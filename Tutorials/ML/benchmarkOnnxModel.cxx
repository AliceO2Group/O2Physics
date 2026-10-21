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

/// Benchmark one dense ONNX input [batch, features] using CPU or AMD MIGraphX.
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/ControlService.h>
#include <Framework/ProcessingContext.h>
#include <Framework/runDataProcessing.h>

#include "ML/3rdparty/GPUORTFloat16.h"
#include <onnxruntime_cxx_api.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <dlfcn.h>

namespace
{
// Optional HIP telemetry without introducing a ROCm build dependency for CPU users.
// Memory is device-wide, not process-local; samples are not a transient peak.
class HipMemory
{
 public:
  explicit HipMemory(int device)
  {
    library = dlopen("libamdhip64.so", RTLD_LAZY | RTLD_LOCAL);
    if (!library) {
      throw std::runtime_error("GPU telemetry requires libamdhip64.so on LD_LIBRARY_PATH");
    }
    auto setDevice = reinterpret_cast<int (*)(int)>(dlsym(library, "hipSetDevice"));
    getMemory = reinterpret_cast<int (*)(size_t*, size_t*)>(dlsym(library, "hipMemGetInfo"));
    if (!setDevice || !getMemory || setDevice(device) != 0) {
      dlclose(library);
      library = nullptr;
      throw std::runtime_error("Cannot select HIP device or resolve HIP memory API");
    }
  }
  ~HipMemory()
  {
    if (library) {
      dlclose(library);
    }
  }
  HipMemory(HipMemory const&) = delete;
  HipMemory& operator=(HipMemory const&) = delete;
  uint64_t used() const
  {
    size_t free = 0, total = 0;
    if (getMemory(&free, &total) != 0 || free > total) {
      throw std::runtime_error("hipMemGetInfo failed");
    }
    return total - free;
  }

 private:
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
  Configurable<std::string> provider{"provider", "cpu", "cpu or migraphx"};
  Configurable<int> gpuId{"gpu-id", 0, "HIP logical device ID (after visibility masking)"};
  Configurable<bool> fp16{"fp16", false, "Require FP16 input; does not convert model weights"};
  Configurable<int> batchSize{"batch-size", 1024, "Samples per inference"};
  Configurable<int> features{"features", 5, "Features per sample, in model order"};
  Configurable<int> warmup{"warmup", 10, "Untimed warm-up calls"};
  Configurable<int> iterations{"iterations", 100, "Timed inference calls"};
  Configurable<int> cpuThreads{"cpu-threads", 1, "ORT intra-op CPU threads"};
  Configurable<std::string> outputCsv{"output-csv", "onnx-benchmark.csv", "CSV output (overwritten)"};

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
    bool gpu = provider.value == "migraphx";
    if (!gpu && provider.value != "cpu") {
      throw std::runtime_error("provider must be cpu or migraphx");
    }
    if (!gpu && fp16.value) {
      throw std::runtime_error("Use FP32 for the CPU baseline");
    }
    std::unique_ptr<HipMemory> memory;
    uint64_t before = 0, afterWarmup = 0, afterRun = 0;
    Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "tutorial-benchmark");
    Ort::SessionOptions options;
    options.SetIntraOpNumThreads(cpuThreads.value);
    options.SetInterOpNumThreads(1);
    options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
    if (gpu) {
      auto providers = Ort::GetAvailableProviders();
      if (std::find(providers.begin(), providers.end(), "MIGraphXExecutionProvider") == providers.end()) {
        throw std::runtime_error("This ONNX Runtime has no MIGraphXExecutionProvider");
      }
      memory = std::make_unique<HipMemory>(gpuId.value);
      before = memory->used();
      OrtMIGraphXProviderOptions amd{};
      amd.device_id = gpuId.value;
      amd.migraphx_mem_limit = std::numeric_limits<size_t>::max();
      // Keep the supplied model precision; do not automatically quantize FP32.
      amd.migraphx_fp16_enable = 0;
      options.AppendExecutionProvider_MIGraphX(amd);
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
    // Replace this feature filling with AO2D extraction and training preprocessing later.
    std::vector<float> floats(fp16.value ? 0 : count);
    std::vector<o2::OrtDataType::Float16_t> halves(fp16.value ? count : 0);
    static_assert(sizeof(o2::OrtDataType::Float16_t) == 2);
    for (size_t i = 0; i < count; ++i) {
      float value = static_cast<float>(i % 101) / 100.f - 0.5f;
      if (fp16.value) {
        halves[i] = o2::OrtDataType::Float16_t(value);
      } else {
        floats[i] = value;
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
    for (size_t i = 0; i < session.GetOutputCount(); ++i) {
      names.emplace_back(session.GetOutputNameAllocated(i, allocator).get());
    }
    std::vector<const char*> outputNames;
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
