// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// Benchmark one dense ONNX input [batch, features] using CPU, AMD MIGraphX or NVIDIA CUDA.
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/ControlService.h>
#include <Framework/ProcessingContext.h>
#include <Framework/runDataProcessing.h>

#include <onnxruntime_c_api.h>
#include <onnxruntime_cxx_api.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <ios>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <dlfcn.h>

namespace
{
// Load GPU telemetry at runtime so CPU builds need no HIP or CUDA headers.
// Memory is device-wide, not process-local; samples are not a transient peak.
class GpuMemory
{
 public:
  explicit GpuMemory(int device, bool cuda)
  {
    const std::vector<const char*> libraries = cuda
                                                 ? std::vector<const char*>{"libcudart.so", "libcudart.so.13", "libcudart.so.12", "libcudart.so.11.0"}
                                                 : std::vector<const char*>{"libamdhip64.so"};
    for (auto name : libraries) {
      library = dlopen(name, RTLD_LAZY | RTLD_LOCAL);
      if (library) {
        break;
      }
    }
    if (!library) {
      throw std::runtime_error(cuda ? "GPU telemetry requires libcudart.so (CUDA 11/12/13) on LD_LIBRARY_PATH"
                                    : "GPU telemetry requires libamdhip64.so on LD_LIBRARY_PATH");
    }
    auto setDevice = reinterpret_cast<int (*)(int)>(dlsym(library, cuda ? "cudaSetDevice" : "hipSetDevice"));
    memoryFunction = cuda ? "cudaMemGetInfo" : "hipMemGetInfo";
    getMemory = reinterpret_cast<int (*)(size_t*, size_t*)>(dlsym(library, memoryFunction));
    if (!setDevice || !getMemory || setDevice(device) != 0) {
      dlclose(library);
      library = nullptr;
      throw std::runtime_error(cuda ? "Cannot select CUDA device or resolve CUDA memory API"
                                    : "Cannot select HIP device or resolve HIP memory API");
    }
  }
  ~GpuMemory()
  {
    if (library) {
      dlclose(library);
    }
  }
  GpuMemory(GpuMemory const&) = delete;
  GpuMemory& operator=(GpuMemory const&) = delete;
  uint64_t used() const
  {
    size_t free = 0, total = 0;
    if (getMemory(&free, &total) != 0 || free > total) {
      throw std::runtime_error(std::string(memoryFunction) + " failed");
    }
    return total - free;
  }

 private:
  const char* memoryFunction = nullptr;
  void* library = nullptr;
  int (*getMemory)(size_t*, size_t*) = nullptr;
};

std::string csvQuote(std::string const& value)
{
  std::string result = "\"";
  for (char c : value) {
    result += c;
    if (c == '"') {
      result += '"';
    }
  }
  return result + '"';
}
} // namespace

using namespace o2::framework;

struct BenchmarkOnnxModel {
  Configurable<std::string> modelPath{"model-path", "model.onnx", "Local ONNX model"};
  Configurable<std::string> provider{"provider", "cpu", "cpu, migraphx or cuda"};
  Configurable<int> gpuId{"gpu-id", 0, "HIP/CUDA logical device ID (after visibility masking)"};
  Configurable<bool> fp16{"fp16", false, "Require FP16 input; does not convert model weights"};
  Configurable<int> batchSize{"batch-size", 1024, "Samples per inference"};
  Configurable<int> features{"features", 5, "Features per sample, in model order"};
  Configurable<int> warmup{"warmup", 10, "Untimed warm-up calls"};
  Configurable<int> iterations{"iterations", 100, "Timed inference calls"};
  Configurable<int> cpuThreads{"cpu-threads", 1, "ORT intra-op CPU threads"};
  Configurable<std::string> outputCsv{"output-csv", "onnx-benchmark.csv", "CSV output (overwritten)"};
  Configurable<std::string> inputFile{"input-file", "", "Optional whitespace-separated feature matrix [batch-size, features]; empty uses synthetic inputs"};
  Configurable<std::string> predictionsFile{"predictions-file", "", "Optional first-output matrix as text (overwritten); extra inference outside timing"};

  void run(ProcessingContext& pc)
  {
    benchmark();
    pc.services().get<ControlService>().endOfStream();
    pc.services().get<ControlService>().readyToQuit(QuitRequest::Me);
  }

  void benchmark()
  {
    if (batchSize.value <= 0 || features.value <= 0 || iterations.value <= 0 || warmup.value < 1 || cpuThreads.value <= 0 || gpuId.value < 0) {
      throw std::runtime_error("Require positive batch, features, iterations, warmup, threads and nonnegative GPU ID");
    }
    bool cuda = provider.value == "cuda";
    bool gpu = cuda || provider.value == "migraphx";
    if (!gpu && provider.value != "cpu") {
      throw std::runtime_error("provider must be cpu, migraphx or cuda");
    }
    if (!gpu && fp16.value) {
      throw std::runtime_error("Use FP32 for the CPU baseline");
    }
    std::unique_ptr<GpuMemory> memory;
    uint64_t before = 0, afterWarmup = 0, afterRun = 0;
    Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "tutorial-benchmark");
    Ort::SessionOptions options;
    options.SetIntraOpNumThreads(cpuThreads.value);
    options.SetInterOpNumThreads(1);
    options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
    if (gpu) {
      const std::string providerName = cuda ? "CUDAExecutionProvider" : "MIGraphXExecutionProvider";
      auto providers = Ort::GetAvailableProviders();
      if (std::find(providers.begin(), providers.end(), providerName) == providers.end()) {
        throw std::runtime_error("This ONNX Runtime has no " + providerName);
      }
      memory = std::make_unique<GpuMemory>(gpuId.value, cuda);
      before = memory->used();
      if (cuda) {
        // V2 options allow disabling TF32 for the FP32 comparison.
        OrtCUDAProviderOptionsV2* raw = nullptr;
        Ort::ThrowOnError(Ort::GetApi().CreateCUDAProviderOptions(&raw));
        std::unique_ptr<OrtCUDAProviderOptionsV2, decltype(Ort::GetApi().ReleaseCUDAProviderOptions)> nvidia(raw, Ort::GetApi().ReleaseCUDAProviderOptions);
        const std::string device = std::to_string(gpuId.value);
        const char* keys[] = {"device_id", "use_tf32"};
        const char* values[] = {device.c_str(), "0"};
        Ort::ThrowOnError(Ort::GetApi().UpdateCUDAProviderOptions(nvidia.get(), keys, values, 2));
        options.AppendExecutionProvider_CUDA_V2(*nvidia);
      } else {
        OrtMIGraphXProviderOptions amd{};
        amd.device_id = gpuId.value;
        amd.migraphx_mem_limit = std::numeric_limits<size_t>::max();
        // Keep the supplied model precision; do not automatically quantize FP32.
        amd.migraphx_fp16_enable = 0;
        options.AppendExecutionProvider_MIGraphX(amd);
      }
      options.AddConfigEntry("session.disable_cpu_ep_fallback", "1");
    }
    Ort::Session session(env, modelPath.value.c_str(), options);
    if (session.GetInputCount() != 1 || session.GetOutputCount() < 1) {
      throw std::runtime_error("Tutorial expects exactly one dense input and at least one output");
    }
    auto inputInfo = session.GetInputTypeInfo(0);
    auto tensorInfo = inputInfo.GetTensorTypeAndShapeInfo();
    auto shape = tensorInfo.GetShape();
    auto type = fp16.value ? ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT16 : ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT;
    if (tensorInfo.GetElementType() != type || shape.size() != 2 ||
        (shape[0] >= 0 && shape[0] != batchSize.value) || (shape[1] >= 0 && shape[1] != features.value)) {
      throw std::runtime_error("Input must match fp16 flag and [batch-size, features]; only dynamic dimensions may be overridden");
    }
    shape = {batchSize.value, features.value};
    const size_t count = static_cast<size_t>(batchSize.value) * features.value;
    if (count > std::numeric_limits<size_t>::max() / sizeof(float)) {
      throw std::runtime_error("Input size overflow");
    }
    // File loading and FP16 conversion are outside the timing loop. Without a
    // file, retain deterministic synthetic inputs for performance-only runs.
    std::vector<float> floats(fp16.value ? 0 : count);
    std::vector<_Float16> halves(fp16.value ? count : 0);
    static_assert(sizeof(_Float16) == 2);
    auto store = [&](size_t i, float value) {
      if (!std::isfinite(value)) {
        throw std::runtime_error("Input features must be finite");
      }
      if (fp16.value) {
        halves[i] = _Float16(value);
        if (!std::isfinite((float)halves[i])) {
          throw std::runtime_error("Input feature overflows FP16");
        }
      } else {
        floats[i] = value;
      }
    };
    if (inputFile.value.empty()) {
      for (size_t i = 0; i < count; ++i) {
        store(i, static_cast<float>(i % 101) / 100.f - 0.5f);
      }
    } else {
      std::ifstream file(inputFile.value);
      if (!file) {
        throw std::runtime_error("Cannot open input-file: " + inputFile.value);
      }
      std::string line;
      for (int row = 0; row < batchSize.value; ++row) {
        if (!std::getline(file, line)) {
          throw std::runtime_error("input-file has fewer rows than batch-size");
        }
        std::istringstream values(line);
        for (int column = 0; column < features.value; ++column) {
          float value = 0;
          if (!(values >> value)) {
            throw std::runtime_error("input-file has a missing or invalid feature at row " + std::to_string(row + 1));
          }
          store(static_cast<size_t>(row) * features.value + column, value);
        }
        values >> std::ws;
        if (!values.eof()) {
          throw std::runtime_error("input-file has extra columns at row " + std::to_string(row + 1));
        }
      }
      file >> std::ws;
      if (!file.eof() || file.bad()) {
        throw std::runtime_error("input-file has extra rows or a read error");
      }
    }
    auto host = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
    void* data = fp16.value ? static_cast<void*>(halves.data()) : static_cast<void*>(floats.data());
    size_t bytes = count * (fp16.value ? sizeof(halves[0]) : sizeof(floats[0]));
    auto input = Ort::Value::CreateTensor(host, data, bytes, shape.data(), shape.size(), type);
    Ort::AllocatorWithDefaultOptions allocator;
    auto inputName = session.GetInputNameAllocated(0, allocator);
    const char* inputNames[] = {inputName.get()};
    std::vector<std::string> names;
    names.reserve(session.GetOutputCount());
    for (size_t i = 0; i < session.GetOutputCount(); ++i) {
      names.emplace_back(session.GetOutputNameAllocated(i, allocator).get());
    }
    std::vector<const char*> outputNames;
    outputNames.reserve(names.size());
    for (auto const& name : names) {
      outputNames.push_back(name.c_str());
    }
    Ort::RunOptions runOptions;
    auto infer = [&]() {
      // Default synchronous Run with CPU input/output includes host/device transfers
      // and output allocation/destruction. No asynchronous Run option is enabled.
      auto outputs = session.Run(runOptions, inputNames, &input, 1, outputNames.data(), outputNames.size());
    };
    for (int i = 0; i < warmup.value; ++i) {
      infer();
    }
    if (memory) {
      afterWarmup = memory->used();
    }
    auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < iterations.value; ++i) {
      infer();
    }
    double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (memory) {
      afterRun = memory->used();
    }
    if (!predictionsFile.value.empty()) {
      // A separate, untimed run leaves benchmark allocation/destruction intact.
      auto outputs = session.Run(runOptions, inputNames, &input, 1, outputNames.data(), outputNames.size());
      auto const& output = outputs.front();
      if (!output.IsTensor()) {
        throw std::runtime_error("Prediction export requires a dense tensor as the first output");
      }
      auto info = output.GetTensorTypeAndShapeInfo();
      auto outputShape = info.GetShape();
      auto outputType = info.GetElementType();
      if (outputShape.size() != 2 || outputShape[0] != batchSize.value || outputShape[1] <= 0 ||
          (outputType != ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT && outputType != ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT16)) {
        throw std::runtime_error("Prediction export requires an FP32 or FP16 first output [batch-size, classes]");
      }
      const float* fp32 = outputType == ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT ? output.GetTensorData<float>() : nullptr;
      const auto* fp16Data = outputType == ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT16 ? output.GetTensorData<_Float16>() : nullptr;
      const size_t columns = static_cast<size_t>(outputShape[1]);
      // Validate before opening the destination, so invalid predictions do not
      // replace an existing file with a partial matrix.
      auto valueAt = [&](size_t i) { return fp32 ? fp32[i] : (float)fp16Data[i]; };
      for (size_t i = 0; i < info.GetElementCount(); ++i) {
        if (!std::isfinite(valueAt(i))) {
          throw std::runtime_error("Prediction export found a nonfinite value");
        }
      }
      std::ofstream predictions(predictionsFile.value, std::ios::trunc);
      if (!predictions) {
        throw std::runtime_error("Cannot open predictions-file: " + predictionsFile.value);
      }
      predictions << std::setprecision(std::numeric_limits<float>::max_digits10);
      for (int row = 0; row < batchSize.value; ++row) {
        for (size_t column = 0; column < columns; ++column) {
          if (column) {
            predictions << ' ';
          }
          predictions << valueAt(static_cast<size_t>(row) * columns + column);
        }
        predictions << '\n';
      }
      predictions.close();
      if (!predictions) {
        throw std::runtime_error("Failed writing predictions-file");
      }
    }
    std::ofstream csv(outputCsv.value, std::ios::trunc);
    if (!csv) {
      throw std::runtime_error("Cannot open output-csv");
    }
    csv << "model,ort_version,provider,precision,gpu_id,batch_size,features,warmup,iterations,cpu_threads,seconds,samples_per_second,gpu_used_before_bytes,gpu_used_after_warmup_bytes,gpu_used_after_run_bytes\n";
    csv << csvQuote(modelPath.value) << ',' << csvQuote(Ort::GetVersionString()) << ',' << provider.value << ','
        << (fp16.value ? "fp16" : "fp32") << ',' << (gpu ? gpuId.value : -1) << ','
        << batchSize.value << ',' << features.value << ',' << warmup.value << ',' << iterations.value << ',' << cpuThreads.value << ','
        << std::setprecision(12) << seconds << ',' << (static_cast<double>(batchSize.value) * iterations.value / seconds) << ',';
    if (memory) {
      csv << before << ',' << afterWarmup << ',' << afterRun;
    } else {
      csv << ",,";
    }
    csv << '\n';
    csv.close();
    if (!csv) {
      throw std::runtime_error("Failed writing output-csv");
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<BenchmarkOnnxModel>(cfgc)};
}
