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

///
/// \file     model.h
///
/// \author   Christian Sonnabend <christian.sonnabend@cern.ch>
///
/// \brief    A general-purpose class for ONNX models
///

#ifndef TOOLS_ML_MODEL_H_
#define TOOLS_ML_MODEL_H_

#include <Framework/Logger.h>

#include <onnx/onnx_pb.h>
#include <onnxruntime_c_api.h>
#include <onnxruntime_cxx_api.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iterator>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace o2
{

namespace ml
{

// ============================================================================
// FusedPreprocLinear custom op
// ----------------------------------------------------------------------------
// Fuses the per-column preprocessing AND the first linear layer into a single
// streaming kernel: for every selected (masked-in) track it computes the K
// features on a small stack and immediately applies W/b, writing one [H] output
// row. No [M, K] feature buffer and no per-op intermediate tensors are ever
// materialised, so the peak memory is just the [M, H] output. The op is generic
// and driven entirely by node attributes (the same spec as PreprocFeature), so
// the original network's layers 1..N are copied verbatim on top of its output.
//
// Op enum (matches OnnxModel::PreprocFeature::Op order):
//   0 Passthrough  1 BroadcastScalar  2 NClsSqrtRecip  3 Mod2  4 GatherNormWhere
// Inputs are a single heterogeneous variadic list: the raw spec inputs (in
// declaration order), then W [K,H] (transB already folded in), then optional b.
struct FusedPreprocLinearKernel {
  int mK = 0, mH = 0, mMaskIdx = 0, mWeightIdx = 0, mBiasIdx = -1;
  std::vector<int64_t> mOps, mA, mB, mScale, mFallback;
  std::vector<float> mConsts;

  explicit FusedPreprocLinearKernel(const OrtKernelInfo* info)
  {
    Ort::ConstKernelInfo ki(info);
    mK = static_cast<int>(ki.GetAttribute<int64_t>("n_features"));
    mH = static_cast<int>(ki.GetAttribute<int64_t>("hidden"));
    mMaskIdx = static_cast<int>(ki.GetAttribute<int64_t>("mask_index"));
    mWeightIdx = static_cast<int>(ki.GetAttribute<int64_t>("weight_index"));
    mBiasIdx = static_cast<int>(ki.GetAttribute<int64_t>("bias_index"));
    mOps = ki.GetAttributes<int64_t>("ops");
    mA = ki.GetAttributes<int64_t>("a_idx");
    mB = ki.GetAttributes<int64_t>("b_idx");
    mScale = ki.GetAttributes<int64_t>("scale_idx");
    mFallback = ki.GetAttributes<int64_t>("fallback_idx");
    mConsts = ki.GetAttributes<float>("consts");
  }

  void Compute(OrtKernelContext* context)
  {
    Ort::KernelContext ctx(context);
    const size_t nIn = ctx.GetInputCount();
    struct In {
      const void* p = nullptr;
      ONNXTensorElementDataType t = ONNX_TENSOR_ELEMENT_DATA_TYPE_UNDEFINED;
    };
    std::vector<In> in(nIn);
    std::vector<Ort::ConstValue> hold;
    hold.reserve(nIn);
    for (size_t i = 0; i < nIn; ++i) {
      hold.emplace_back(ctx.GetInput(i));
      in[i] = {hold[i].GetTensorRawData(), hold[i].GetTensorTypeAndShapeInfo().GetElementType()};
    }

    auto readF = [&](int idx, int64_t n) -> float {
      const In& x = in[idx];
      switch (x.t) {
        case ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT:
          return reinterpret_cast<const float*>(x.p)[n];
        case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT32:
          return static_cast<float>(reinterpret_cast<const int32_t*>(x.p)[n]);
        case ONNX_TENSOR_ELEMENT_DATA_TYPE_UINT8:
          return static_cast<float>(reinterpret_cast<const uint8_t*>(x.p)[n]);
        case ONNX_TENSOR_ELEMENT_DATA_TYPE_INT8:
          return static_cast<float>(reinterpret_cast<const int8_t*>(x.p)[n]);
        default:
          return 0.f;
      }
    };
    auto readI = [&](int idx, int64_t n) -> int64_t {
      const In& x = in[idx];
      if (x.t == ONNX_TENSOR_ELEMENT_DATA_TYPE_INT32) {
        return reinterpret_cast<const int32_t*>(x.p)[n];
      }
      return static_cast<int64_t>(readF(idx, n));
    };

    const bool* mask = reinterpret_cast<const bool*>(in[mMaskIdx].p);
    const int64_t N = hold[mMaskIdx].GetTensorTypeAndShapeInfo().GetShape()[0];
    int64_t M = 0;
    for (int64_t n = 0; n < N; ++n) {
      if (mask[n]) {
        ++M;
      }
    }

    const float* W = reinterpret_cast<const float*>(in[mWeightIdx].p); // [K, H]
    const float* b = (mBiasIdx >= 0) ? reinterpret_cast<const float*>(in[mBiasIdx].p) : nullptr;

    auto out = ctx.GetOutput(0, {M, static_cast<int64_t>(mH)});
    float* o = out.GetTensorMutableData<float>();

    std::vector<float> feat(mK);
    int64_t m = 0;
    for (int64_t n = 0; n < N; ++n) {
      if (!mask[n]) {
        continue;
      }
      for (int k = 0; k < mK; ++k) {
        const float c0 = mConsts[3 * k], c1 = mConsts[3 * k + 1], c2 = mConsts[3 * k + 2];
        switch (mOps[k]) {
          case 0: // Passthrough
            feat[k] = readF(static_cast<int>(mA[k]), n);
            break;
          case 1: // BroadcastScalar (e.g. mass, a [1] scalar)
            feat[k] = readF(static_cast<int>(mA[k]), 0);
            break;
          case 2: { // NClsSqrtRecip: sqrt(numer / (a - b))
            const float ncl = readF(static_cast<int>(mA[k]), n) - readF(static_cast<int>(mB[k]), n);
            const float numer = mScale[k] >= 0 ? readF(static_cast<int>(mScale[k]), 0) : c0;
            feat[k] = std::sqrt(numer / ncl);
            break;
          }
          case 3: { // Mod2: fmod(fmod(a, c0) + c1, c2)
            const float a = readF(static_cast<int>(mA[k]), n);
            feat[k] = std::fmod(std::fmod(a, c0) + c1, c2);
            break;
          }
          case 4: { // GatherNormWhere: idx<0 ? fallback : array[idx]/denom
            const int64_t idx = readI(static_cast<int>(mB[k]), n);
            if (idx < 0) {
              feat[k] = mFallback[k] >= 0 ? readF(static_cast<int>(mFallback[k]), 0) : c1;
            } else {
              const float denom = mScale[k] >= 0 ? readF(static_cast<int>(mScale[k]), 0) : c0;
              feat[k] = readF(static_cast<int>(mA[k]), idx) / denom;
            }
            break;
          }
          default:
            feat[k] = 0.f;
            break;
        }
      }
      float* orow = o + m * mH;
      for (int h = 0; h < mH; ++h) {
        orow[h] = b ? b[h] : 0.f;
      }
      for (int k = 0; k < mK; ++k) {
        const float fk = feat[k];
        const float* wrow = W + static_cast<size_t>(k) * mH;
        for (int h = 0; h < mH; ++h) {
          orow[h] += fk * wrow[h];
        }
      }
      ++m;
    }
  }
};

struct FusedPreprocLinearOp : Ort::CustomOpBase<FusedPreprocLinearOp, FusedPreprocLinearKernel> {
  void* CreateKernel(const OrtApi&, const OrtKernelInfo* info) const { return new FusedPreprocLinearKernel(info); }
  const char* GetName() const { return "FusedPreprocLinear"; }
  size_t GetInputTypeCount() const { return 1; }
  ONNXTensorElementDataType GetInputType(size_t) const { return ONNX_TENSOR_ELEMENT_DATA_TYPE_UNDEFINED; }
  OrtCustomOpInputOutputCharacteristic GetInputCharacteristic(size_t) const { return OrtCustomOpInputOutputCharacteristic::INPUT_OUTPUT_VARIADIC; }
  bool GetVariadicInputHomogeneity() const { return false; } // heterogeneous: mixed input types
  int GetVariadicInputMinArity() const { return 1; }
  size_t GetOutputTypeCount() const { return 1; }
  ONNXTensorElementDataType GetOutputType(size_t) const { return ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT; }
};

class OnnxModel
{

 public:
  OnnxModel() = default;
  ~OnnxModel() = default;

  // Inferencing
  void initModel(const std::string&, const bool = false, const int = 0, const uint64_t = 0, const uint64_t = 0);

  // template methods -- best to define them in header
  template <typename T>
  T* evalModel(std::vector<Ort::Value>& input)
  {
    LOG(debug) << "Input tensor shape: " << printShape(input[0].GetTensorTypeAndShapeInfo().GetShape());
    // assert(input[0].GetTensorTypeAndShapeInfo().GetShape() == getNumInputNodes()); --> Fails build in debug mode, TODO: assertion should be checked somehow

    try {
      const Ort::RunOptions runOptions;
      std::vector<const char*> inputNamesChar(mInputNames.size(), nullptr);
      std::transform(std::begin(mInputNames), std::end(mInputNames), std::begin(inputNamesChar),
                     [&](const std::string& str) { return str.c_str(); });

      std::vector<const char*> outputNamesChar(mOutputNames.size(), nullptr);
      std::transform(std::begin(mOutputNames), std::end(mOutputNames), std::begin(outputNamesChar),
                     [&](const std::string& str) { return str.c_str(); });
      auto outputTensors = mSession->Run(runOptions, inputNamesChar.data(), input.data(), input.size(), outputNamesChar.data(), outputNamesChar.size());
      LOG(debug) << "Number of output tensors: " << outputTensors.size();
      if (outputTensors.size() != mOutputNames.size()) {
        LOG(fatal) << "Number of output tensors: " << outputTensors.size() << " does not agree with the model specified size: " << mOutputNames.size();
      }
      for (std::size_t i = 0; i < outputTensors.size(); i++) {
        LOG(debug) << "Output tensor shape: " << printShape(outputTensors[i].GetTensorTypeAndShapeInfo().GetShape());
        if ((outputTensors[i].GetTensorTypeAndShapeInfo().GetShape() != mOutputShapes[i]) && (mOutputShapes[i][0] != -1)) {
          LOG(fatal) << "Shape of tensor " << i << " does not agree with model specification! Output: " << printShape(outputTensors[i].GetTensorTypeAndShapeInfo().GetShape()) << " model: " << printShape(mOutputShapes[i]);
        }
      }
      T* outputValues = outputTensors.back().GetTensorMutableData<T>();
      return outputValues;
    } catch (const Ort::Exception& exception) {
      LOG(error) << "Error running model inference: " << exception.what();
    }
    return nullptr;
  }

  template <typename T>
  T* evalModel(std::vector<T>& input)
  {
    const int64_t size = input.size();
    assert(size % mInputShapes[0][1] == 0);
    std::vector<int64_t> inputShape{size / mInputShapes[0][1], mInputShapes[0][1]};
    std::vector<Ort::Value> inputTensors;
    Ort::MemoryInfo memInfo =
      Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);
    inputTensors.emplace_back(Ort::Value::CreateTensor<T>(memInfo, input.data(), size, inputShape.data(), inputShape.size()));
    LOG(debug) << "Input shape calculated from vector: " << printShape(inputShape);
    return evalModel<T>(inputTensors);
  }

  // For 2D inputs
  template <typename T>
  T* evalModel(std::vector<std::vector<T>>& input)
  {
    std::vector<Ort::Value> inputTensors;

    Ort::MemoryInfo memInfo = Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);

    for (std::size_t iinput = 0; iinput < input.size(); iinput++) {
      [[maybe_unused]] int totalSize = 1;
      int64_t size = input[iinput].size();
      for (std::size_t idim = 1; idim < mInputShapes[iinput].size(); idim++) {
        totalSize *= mInputShapes[iinput][idim];
      }
      assert(size % totalSize == 0);

      std::vector<int64_t> inputShape{static_cast<int64_t>(size / totalSize)};
      for (std::size_t idim = 1; idim < mInputShapes[iinput].size(); idim++) {
        inputShape.push_back(mInputShapes[iinput][idim]);
      }

      inputTensors.emplace_back(Ort::Value::CreateTensor<T>(memInfo, input[iinput].data(), size, inputShape.data(), inputShape.size()));
    }

    return evalModel<T>(inputTensors);
  }

  // Reset session
  void resetSession()
  {
    mSession.reset(new Ort::Session{*mEnv, modelPath.c_str(), sessionOptions});
  }

  // Declaration of a raw graph input fed directly from an Arrow buffer.
  struct PreprocInput {
    enum class Type { TrackFloat,     // per-track float column           [N]
                      TrackInt32,     // per-track int32 column           [N]
                      TrackUint8,     // per-track uint8 column           [N]
                      TrackInt8,      // per-track int8 column            [N]
                      TrackBool,      // per-track bool mask              [N]
                      CollisionFloat, // per-collision float array        [C]
                      ScalarFloat };  // single scalar (e.g. mass)        [1]
    std::string name;
    Type type;
  };

  // Preprocessing recipe for one network feature (produces a [N] float tensor
  // that feeds column i of the decomposed first layer).
  struct PreprocFeature {
    enum class Op {
      Passthrough,     // feature = a                                   (a: TrackFloat)
      BroadcastScalar, // feature = Expand(a, shape(shapeRef))          (a: ScalarFloat)
      NClsSqrtRecip,   // feature = Sqrt(c0 / (float(a) - float(b)))    (a,b: Track int cols)
      Mod2,            // feature = Mod(Mod(a, c0) + c1, c2)            (a: TrackFloat)
      GatherNormWhere  // feature = Where(b<0, fb, Gather(a, b) / c0)   (a: CollisionFloat, b: TrackInt32)
    };
    Op op;
    std::string a;             // primary input
    std::string b;             // secondary input (NCls: 2nd col; Gather: index)
    std::string shapeRef;      // BroadcastScalar: [N] input to size the Expand
    std::string fallbackInput; // GatherNormWhere: scalar input for the b<0 fallback; "" => use c[1]
    std::string scaleInput;    // NCls: numerator scalar input; Gather: divisor scalar input; "" => use c[0]
    std::array<float, 3> c{};  // op constants
  };

  // Rebuild the model from scratch so the network reads its raw Arrow inputs
  // directly and performs all preprocessing + the first linear layer inside the
  // graph.  Each feature column is produced by a small preprocessing subgraph
  // (`features`) from the raw inputs (`inputDefs`), then the first linear layer is
  // decomposed
  //   layer0 = X @ W + b = sum_i (feat_i[N,1] @ W_row_i[1,H]) + b
  // so no [N, K] interleaving / Concat buffer is ever materialised.  The original
  // layers 1..N are copied verbatim on top of the decomposed layer-0 output.
  // Building a fresh model (rather than augmenting the existing one) is required:
  // the Model Editor can only add nodes, so the original layer-0 Gemm would
  // otherwise remain and collide on the layer-0 output name.
  // If maskInput is non-empty it names a bool [N] input; each feature is then
  // Compress'd by it so the matmul runs only on the selected (valid) rows and the
  // output is the compact set of selected tracks, in order.
  void setupColumnInputs(const std::vector<PreprocInput>& inputDefs,
                         const std::vector<PreprocFeature>& features,
                         const std::string& maskInput = "")
  {
    const int numFeatures = static_cast<int>(features.size());
    if (numFeatures != mInputShapes[0][1]) {
      LOG(fatal) << "setupColumnInputs: expected " << mInputShapes[0][1] << " features, got " << numFeatures;
      return;
    }

    onnx::ModelProto onnxModel;
    {
      std::ifstream ifs(modelPath, std::ios::binary);
      if (!ifs || !onnxModel.ParseFromIstream(&ifs)) {
        LOG(fatal) << "setupColumnInputs: failed to parse ONNX model from " << modelPath;
        return;
      }
    }
    const auto& og = onnxModel.graph();

    auto findInit = [&](const std::string& name) -> const onnx::TensorProto* {
      for (int i = 0; i < og.initializer_size(); ++i) {
        if (og.initializer(i).name() == name) {
          return &og.initializer(i);
        }
      }
      return nullptr;
    };
    auto tensorFloats = [&](const onnx::TensorProto* t, int64_t n) {
      std::vector<float> v(n);
      if (t->raw_data().size() > 0) {
        std::memcpy(v.data(), t->raw_data().data(), n * sizeof(float));
      } else {
        for (int64_t i = 0; i < n; ++i) {
          v[i] = t->float_data(i);
        }
      }
      return v;
    };

    // --- locate the first linear layer and extract its weights. The custom op
    //     reproduces this layer's pre-activation output, so the original layers
    //     1..N are copied on top (the first layer and its input Cast are dropped). ---
    const onnx::NodeProto* first = nullptr;
    for (int i = 0; i < og.node_size(); ++i) {
      const auto& n = og.node(i);
      if (n.op_type() == "Gemm" || n.op_type() == "MatMul") {
        first = &n;
        break;
      }
    }
    if (!first) {
      LOG(fatal) << "setupColumnInputs: no Gemm/MatMul layer found in model";
      return;
    }
    const std::string layer0Out = first->output(0);
    const onnx::TensorProto* wT = findInit(first->input(1));
    if (!wT || wT->dims_size() != 2) {
      LOG(fatal) << "setupColumnInputs: first-layer weight initializer not found or not 2D";
      return;
    }
    bool transB = false;
    if (first->op_type() == "Gemm") {
      for (int i = 0; i < first->attribute_size(); ++i) {
        if (first->attribute(i).name() == "transB") {
          transB = (first->attribute(i).i() != 0);
        }
      }
    }
    const int K = transB ? static_cast<int>(wT->dims(1)) : static_cast<int>(wT->dims(0));
    const int H = transB ? static_cast<int>(wT->dims(0)) : static_cast<int>(wT->dims(1));
    if (K != numFeatures) {
      LOG(fatal) << "setupColumnInputs: first-layer K=" << K << " != numFeatures=" << numFeatures;
      return;
    }
    const std::vector<float> wRaw = tensorFloats(wT, static_cast<int64_t>(K) * H);
    std::vector<float> wKH(static_cast<size_t>(K) * H); // row-major [K, H], transB folded in
    for (int k = 0; k < K; ++k) {
      for (int h = 0; h < H; ++h) {
        wKH[static_cast<size_t>(k) * H + h] =
          transB ? wRaw[static_cast<size_t>(h) * K + k] : wRaw[static_cast<size_t>(k) * H + h];
      }
    }
    std::vector<float> bData;
    if (first->op_type() == "Gemm" && first->input_size() >= 3 && !first->input(2).empty()) {
      if (const onnx::TensorProto* bT = findInit(first->input(2))) {
        bData = tensorFloats(bT, H);
      }
    }

    // --- map each raw input name to its position in the custom node's input list,
    //     and encode the per-feature preprocessing spec for the custom kernel. ---
    std::vector<std::string> rawInputNames;
    rawInputNames.reserve(inputDefs.size());
    std::map<std::string, int> inputIndex;
    for (int i = 0; i < static_cast<int>(inputDefs.size()); ++i) {
      inputIndex[inputDefs[i].name] = i;
      rawInputNames.push_back(inputDefs[i].name);
    }
    auto idxOf = [&](const std::string& name) -> int64_t {
      auto it = inputIndex.find(name);
      return it == inputIndex.end() ? -1 : static_cast<int64_t>(it->second);
    };
    std::vector<int64_t> specOps(K), specA(K), specB(K), specScale(K), specFallback(K);
    std::vector<float> specConsts(static_cast<size_t>(3) * K);
    for (int i = 0; i < K; ++i) {
      const auto& f = features[i];
      specOps[i] = static_cast<int64_t>(f.op);
      specA[i] = idxOf(f.a);
      specB[i] = f.b.empty() ? -1 : idxOf(f.b);
      specScale[i] = f.scaleInput.empty() ? -1 : idxOf(f.scaleInput);
      specFallback[i] = f.fallbackInput.empty() ? -1 : idxOf(f.fallbackInput);
      specConsts[3 * i] = f.c[0];
      specConsts[3 * i + 1] = f.c[1];
      specConsts[3 * i + 2] = f.c[2];
    }

    // --- keep the original layers 1..N: every node reachable backwards from the
    //     graph outputs, except the replaced first layer (and the now-dead input
    //     Cast). They consume the custom node's pre-activation output. ---
    std::set<std::string> live;
    for (int i = 0; i < og.output_size(); ++i) {
      live.insert(og.output(i).name());
    }
    std::vector<bool> keep(og.node_size(), false);
    bool changed = true;
    while (changed) {
      changed = false;
      for (int i = 0; i < og.node_size(); ++i) {
        if (keep[i] || &og.node(i) == first) {
          continue;
        }
        const auto& n = og.node(i);
        bool produces = false;
        for (int o = 0; o < n.output_size(); ++o) {
          if (live.count(n.output(o))) {
            produces = true;
            break;
          }
        }
        if (!produces) {
          continue;
        }
        keep[i] = true;
        changed = true;
        for (int in = 0; in < n.input_size(); ++in) {
          live.insert(n.input(in));
        }
      }
    }

    // --- assemble the rebuilt model as an ONNX protobuf and load it from bytes.
    //     (The Model-Editor session path does not resolve user custom ops; the
    //     standard from-bytes session does.) Layers 1..N are copied verbatim; only
    //     the first Gemm + dead input Cast are replaced by the fused custom node. ---
    onnx::ModelProto nm;
    nm.set_ir_version(onnxModel.ir_version());
    for (const auto& oi : onnxModel.opset_import()) {
      *nm.add_opset_import() = oi;
    }
    {
      auto* oi = nm.add_opset_import();
      oi->set_domain("ai.o2.ml");
      oi->set_version(1);
    }
    onnx::GraphProto* ng = nm.mutable_graph();
    ng->set_name(og.name());

    // raw graph inputs
    for (const auto& pin : inputDefs) {
      auto* vp = ng->add_input();
      vp->set_name(pin.name);
      auto* tt = vp->mutable_type()->mutable_tensor_type();
      int et = onnx::TensorProto::FLOAT;
      const char* dim = "N";
      bool scalar = false;
      switch (pin.type) {
        case PreprocInput::Type::TrackFloat: et = onnx::TensorProto::FLOAT; break;
        case PreprocInput::Type::TrackInt32: et = onnx::TensorProto::INT32; break;
        case PreprocInput::Type::TrackUint8: et = onnx::TensorProto::UINT8; break;
        case PreprocInput::Type::TrackInt8: et = onnx::TensorProto::INT8; break;
        case PreprocInput::Type::TrackBool: et = onnx::TensorProto::BOOL; break;
        case PreprocInput::Type::CollisionFloat: et = onnx::TensorProto::FLOAT; dim = "C"; break;
        case PreprocInput::Type::ScalarFloat: et = onnx::TensorProto::FLOAT; scalar = true; break;
      }
      tt->set_elem_type(et);
      auto* sh = tt->mutable_shape();
      if (scalar) {
        sh->add_dim()->set_dim_value(1);
      } else {
        sh->add_dim()->set_dim_param(dim);
      }
    }

    // W [K,H] (transB folded) and optional b [H] as initializers fed to the node.
    auto addRawInit = [&](const std::string& name, const std::vector<float>& data,
                          const std::vector<int64_t>& shape) {
      auto* t = ng->add_initializer();
      t->set_name(name);
      t->set_data_type(onnx::TensorProto::FLOAT);
      for (const auto d : shape) {
        t->add_dims(d);
      }
      t->set_raw_data(std::string(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(float)));
    };
    addRawInit("_l0_weight", wKH, {static_cast<int64_t>(K), static_cast<int64_t>(H)});
    std::vector<std::string> nodeInputs = rawInputNames;
    const int64_t weightIdx = static_cast<int64_t>(nodeInputs.size());
    nodeInputs.push_back("_l0_weight");
    int64_t biasIdx = -1;
    if (!bData.empty()) {
      addRawInit("_l0_bias", bData, {static_cast<int64_t>(H)});
      biasIdx = static_cast<int64_t>(nodeInputs.size());
      nodeInputs.push_back("_l0_bias");
    }

    // the fused preprocessing + first-layer custom node -> layer0Out [M, H]
    {
      auto* cn = ng->add_node();
      cn->set_op_type("FusedPreprocLinear");
      cn->set_domain("ai.o2.ml");
      cn->set_name("_fused_preproc_linear");
      for (const auto& in : nodeInputs) {
        cn->add_input(in);
      }
      cn->add_output(layer0Out);
      auto addAttrI = [&](const char* nm, int64_t v) {
        auto* a = cn->add_attribute();
        a->set_name(nm);
        a->set_type(onnx::AttributeProto::INT);
        a->set_i(v);
      };
      auto addAttrInts = [&](const char* nm, const std::vector<int64_t>& v) {
        auto* a = cn->add_attribute();
        a->set_name(nm);
        a->set_type(onnx::AttributeProto::INTS);
        for (const auto x : v) {
          a->add_ints(x);
        }
      };
      addAttrI("n_features", K);
      addAttrI("hidden", H);
      addAttrI("mask_index", maskInput.empty() ? -1 : idxOf(maskInput));
      addAttrI("weight_index", weightIdx);
      addAttrI("bias_index", biasIdx);
      addAttrInts("ops", specOps);
      addAttrInts("a_idx", specA);
      addAttrInts("b_idx", specB);
      addAttrInts("scale_idx", specScale);
      addAttrInts("fallback_idx", specFallback);
      auto* ac = cn->add_attribute();
      ac->set_name("consts");
      ac->set_type(onnx::AttributeProto::FLOATS);
      for (const auto x : specConsts) {
        ac->add_floats(x);
      }
    }

    // kept original nodes (layers 1..N) + the initializers they reference
    int keptNodes = 0;
    std::set<std::string> neededInits;
    for (int i = 0; i < og.node_size(); ++i) {
      if (!keep[i]) {
        continue;
      }
      const auto& n = og.node(i);
      *ng->add_node() = n;
      ++keptNodes;
      for (int in = 0; in < n.input_size(); ++in) {
        if (findInit(n.input(in))) {
          neededInits.insert(n.input(in));
        }
      }
    }
    for (const auto& name : neededInits) {
      *ng->add_initializer() = *findInit(name);
    }
    for (int i = 0; i < og.output_size(); ++i) {
      *ng->add_output() = og.output(i);
    }

    std::string modelBytes;
    nm.SerializeToString(&modelBytes);

    // Register the FusedPreprocLinear custom op (once) on the session options.
    if (!mFusedOpRegistered) {
      static FusedPreprocLinearOp fusedOp;
      static Ort::CustomOpDomain fusedDomain("ai.o2.ml");
      static bool opAdded = false;
      if (!opAdded) {
        fusedDomain.Add(&fusedOp);
        opAdded = true;
      }
      sessionOptions.Add(fusedDomain);
      mFusedOpRegistered = true;
    }

    mSession = std::make_shared<Ort::Session>(*mEnv, modelBytes.data(), modelBytes.size(), sessionOptions);

    mInputNames = rawInputNames;
    mInputShapes.assign(rawInputNames.size(), std::vector<int64_t>{-1});
    mNumFeatures = K;

    LOG(info) << "setupColumnInputs: rebuilt model with " << rawInputNames.size()
              << " raw inputs -> fused preprocessing + first layer (K=" << K << ", H=" << H
              << ") in one custom op, " << keptNodes << " downstream nodes kept";
  }

  // Getters & Setters
  Ort::SessionOptions* getSessionOptions() { return &sessionOptions; } // For optimizations in post
  std::shared_ptr<Ort::Session> getSession()
  {
    return mSession;
  }
  int getNumInputNodes() const { return mInputShapes[0][1]; }
  bool hasColumnInputs() const { return mInputShapes.size() > 1 || (mInputShapes.size() == 1 && mInputShapes[0].size() == 1); }
  int getNumColumns() const { return static_cast<int>(mInputNames.size()); }
  int getNumFeatures() const { return mNumFeatures; }
  std::vector<std::vector<int64_t>> getInputShapes() const { return mInputShapes; }
  int getNumOutputNodes() const { return mOutputShapes[0][1]; }
  uint64_t getValidityFrom() const { return validFrom; }
  uint64_t getValidityUntil() const { return validUntil; }
  void setActiveThreads(const int);

 private:
  // Environment variables for the ONNX runtime
  std::shared_ptr<Ort::Env> mEnv = nullptr;
  std::shared_ptr<Ort::Session> mSession = nullptr;
  Ort::SessionOptions sessionOptions;

  // Input & Output specifications of the loaded network
  std::vector<std::string> mInputNames;
  std::vector<std::vector<int64_t>> mInputShapes;
  int mNumFeatures = 0;
  bool mFusedOpRegistered = false;
  std::vector<std::string> mOutputNames;
  std::vector<std::vector<int64_t>> mOutputShapes;

  // Environment settings
  std::string modelPath;
  int activeThreads = 0;
  uint64_t validFrom = 0;
  uint64_t validUntil = 0;

  // Internal function for printing the shape of tensors
  std::string printShape(const std::vector<int64_t>&);
  bool checkHyperloop(const bool = true);
};

} // namespace ml

} // namespace o2

#endif // TOOLS_ML_MODEL_H_
