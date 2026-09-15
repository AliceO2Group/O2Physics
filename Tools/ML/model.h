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

#include <onnxruntime_c_api.h>
#include <onnxruntime_cxx_api.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace o2
{

namespace ml
{

/// \brief Thin wrapper around an ONNX Runtime session (CPU only)
///
/// Inference entry points:
///  - evalModelRaw(): runs the session and returns the owning output tensors (one Ort::Value per model output).
///    This is the zero-copy path: the caller owns the tensors and may read them via GetTensorData<T>()
///    for as long as the returned vector is alive.
///  - evalModel(input): convenience wrapper returning a copy of the *last* model output as std::vector<T>.
///  - evalModel(input, output): same as above but writes into a caller-provided vector. The vector is
///    resized as needed and its capacity is reused across calls, which avoids per-call allocations in
///    hot loops (e.g. batched inference).
///
/// Inputs given as std::vector<T> are wrapped in an Ort::Value without copying. The input vector must therefore
/// stay alive until the call returns (which is always the case for the synchronous calls provided here).
class OnnxModel
{

 public:
  OnnxModel() = default;
  ~OnnxModel() = default;
  // The ONNX session options are not copyable, so the model is move-only
  OnnxModel(const OnnxModel&) = delete;
  OnnxModel& operator=(const OnnxModel&) = delete;
  OnnxModel(OnnxModel&&) = default;
  OnnxModel& operator=(OnnxModel&&) = default;

  // Inferencing
  void initModel(const std::string&, const bool = false, const int = 0, const uint64_t = 0, const uint64_t = 0);

  /// Run the model on already prepared input tensors
  /// \param input one Ort::Value per model input
  /// \return the output tensors of the model (owning). Access to the data via output[i].GetTensorData<T>()
  std::vector<Ort::Value> evalModelRaw(std::vector<Ort::Value>& input);

  /// Run the model and copy the last output tensor into a vector
  /// \param input one Ort::Value per model input
  /// \return flattened content of the last output tensor
  template <typename T>
  std::vector<T> evalModel(std::vector<Ort::Value>& input)
  {
    std::vector<T> output;
    evalModel<T>(input, output);
    return output;
  }

  /// Run the model and copy the last output tensor into the provided vector (allocation is reused between calls)
  /// \param input one Ort::Value per model input
  /// \param output vector which is filled with the flattened content of the last output tensor
  template <typename T>
  void evalModel(std::vector<Ort::Value>& input, std::vector<T>& output)
  {
    const std::vector<Ort::Value> outputTensors = evalModelRaw(input);
    copyLastOutput<T>(outputTensors, output);
  }

  /// Run a single-input model on a flat vector of features (batches are inferred from the model input shape)
  /// \param input flattened input features, size must be a multiple of the number of input nodes
  /// \return flattened content of the last output tensor
  template <typename T>
  std::vector<T> evalModel(std::vector<T>& input)
  {
    std::vector<T> output;
    evalModel<T>(input, output);
    return output;
  }

  /// Run a single-input model on a flat vector of features (batches are inferred from the model input shape)
  /// \param input flattened input features, size must be a multiple of the number of input nodes
  /// \param output vector which is filled with the flattened content of the last output tensor
  template <typename T>
  void evalModel(std::vector<T>& input, std::vector<T>& output)
  {
    if (mInputShapes.size() != 1) {
      LOG(fatal) << "Model has " << mInputShapes.size() << " inputs but a single input vector was provided. Use the std::vector<std::vector<T>> or std::vector<Ort::Value> overload.";
    }
    std::vector<Ort::Value> inputTensors;
    inputTensors.reserve(1);
    addInputTensor<T>(inputTensors, input, 0);
    evalModel<T>(inputTensors, output);
  }

  /// Run a multi-input model: one flat vector of features per model input
  /// \param input one flattened vector per model input
  /// \return flattened content of the last output tensor
  template <typename T>
  std::vector<T> evalModel(std::vector<std::vector<T>>& input)
  {
    std::vector<T> output;
    evalModel<T>(input, output);
    return output;
  }

  /// Run a multi-input model: one flat vector of features per model input
  /// \param input one flattened vector per model input
  /// \param output vector which is filled with the flattened content of the last output tensor
  template <typename T>
  void evalModel(std::vector<std::vector<T>>& input, std::vector<T>& output)
  {
    if (input.size() != mInputShapes.size()) {
      LOG(fatal) << "Model has " << mInputShapes.size() << " inputs but " << input.size() << " input vectors were provided.";
    }
    std::vector<Ort::Value> inputTensors;
    inputTensors.reserve(input.size());
    for (std::size_t iinput = 0; iinput < input.size(); iinput++) {
      addInputTensor<T>(inputTensors, input[iinput], iinput);
    }
    evalModel<T>(inputTensors, output);
  }

  // Reset session
  void resetSession()
  {
    mSession.reset(new Ort::Session{*mEnv, modelPath.c_str(), sessionOptions});
  }

  // Getters & Setters
  Ort::SessionOptions* getSessionOptions() { return &sessionOptions; } // For optimizations in post
  std::shared_ptr<Ort::Session> getSession()
  {
    return mSession;
  }
  int getNumInputNodes() const { return mInputShapes[0][1]; }
  std::vector<std::vector<int64_t>> getInputShapes() const { return mInputShapes; }
  int getNumOutputNodes() const { return mOutputShapes[0][1]; }
  std::vector<std::vector<int64_t>> getOutputShapes() const { return mOutputShapes; }
  uint64_t getValidityFrom() const { return validFrom; }
  uint64_t getValidityUntil() const { return validUntil; }
  void setActiveThreads(const int);

 private:
  // Environment variables for the ONNX runtime
  std::shared_ptr<Ort::Env> mEnv = nullptr;
  std::shared_ptr<Ort::Session> mSession = nullptr;
  Ort::SessionOptions sessionOptions;
  Ort::MemoryInfo mMemInfo{nullptr}; // CPU memory info, created once in initModel

  // Input & Output specifications of the loaded network
  std::vector<std::string> mInputNames;
  std::vector<std::vector<int64_t>> mInputShapes;
  std::vector<std::string> mOutputNames;
  std::vector<std::vector<int64_t>> mOutputShapes;

  // Environment settings
  std::string modelPath;
  int activeThreads = 0;
  uint64_t validFrom = 0;
  uint64_t validUntil = 0;

  /// Derive the tensor shape for input number iinput from the number of provided values
  /// The first dimension is treated as batch dimension; all other dimensions are taken from the model
  std::vector<int64_t> inferInputShape(const std::size_t iinput, const int64_t size) const;

  /// Wrap a flat data vector (no copy) into an Ort::Value and append it to tensors
  template <typename T>
  void addInputTensor(std::vector<Ort::Value>& tensors, std::vector<T>& data, const std::size_t iinput) const
  {
    const std::vector<int64_t> inputShape = inferInputShape(iinput, static_cast<int64_t>(data.size()));
    LOG(debug) << "Input shape calculated from vector: " << printShape(inputShape);
    tensors.emplace_back(Ort::Value::CreateTensor<T>(mMemInfo, data.data(), data.size(), inputShape.data(), inputShape.size()));
  }

  /// Copy the content of the last output tensor into output (reusing its allocation)
  template <typename T>
  void copyLastOutput(const std::vector<Ort::Value>& outputTensors, std::vector<T>& output) const
  {
    if (outputTensors.empty()) {
      LOG(fatal) << "Model returned no output tensors";
    }
    const Ort::Value& tensor = outputTensors.back();
    const auto info = tensor.GetTensorTypeAndShapeInfo();
    if (info.GetElementType() != Ort::TypeToTensorType<T>::type) {
      LOG(fatal) << "Requested output type (ONNX type id " << static_cast<int>(Ort::TypeToTensorType<T>::type) << ") does not match the model output tensor type (ONNX type id " << static_cast<int>(info.GetElementType()) << ")";
    }
    const T* data = tensor.GetTensorData<T>();
    output.assign(data, data + info.GetElementCount());
  }

  // Internal function for printing the shape of tensors
  static std::string printShape(const std::vector<int64_t>&);
  bool checkHyperloop(const bool = true);
};

} // namespace ml

} // namespace o2

#endif // TOOLS_ML_MODEL_H_
