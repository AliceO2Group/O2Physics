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

#include <onnxruntime_cxx_api.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

namespace o2
{

namespace ml
{

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

  // Getters & Setters
  Ort::SessionOptions* getSessionOptions() { return &sessionOptions; } // For optimizations in post
  std::shared_ptr<Ort::Session> getSession()
  {
    return mSession;
  }
  int getNumInputNodes() const { return mInputShapes[0][1]; }
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
