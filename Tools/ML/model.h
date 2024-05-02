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

// C++ and system includes
#if __has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>
#else
#include <onnxruntime_cxx_api.h>
#endif
#include <vector>
#include <string>
#include <memory>
#include <map>

// ROOT includes
#include "TSystem.h"

// O2 includes
#include "Framework/Logger.h"

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
  void initModel(std::string, bool = false, int = 0, uint64_t = 0, uint64_t = 0);

  // template methods -- best to define them in header
  template <typename T>
  T* evalModel(std::vector<Ort::Value>& input)
  {
    LOG(debug) << "Input tensor shape: " << printShape(input[0].GetTensorTypeAndShapeInfo().GetShape());
    // assert(input[0].GetTensorTypeAndShapeInfo().GetShape() == getNumInputNodes()); --> Fails build in debug mode, TODO: assertion should be checked somehow

    try {
#if __has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
      auto outputTensors = mSession->Run(mInputNames, input, mOutputNames);
#else
      Ort::RunOptions runOptions;
      std::vector<const char*> inputNamesChar(mInputNames.size(), nullptr);
      std::transform(std::begin(mInputNames), std::end(mInputNames), std::begin(inputNamesChar),
                     [&](const std::string& str) { return str.c_str(); });

      std::vector<const char*> outputNamesChar(mOutputNames.size(), nullptr);
      std::transform(std::begin(mOutputNames), std::end(mOutputNames), std::begin(outputNamesChar),
                     [&](const std::string& str) { return str.c_str(); });
      auto outputTensors = mSession->Run(runOptions, inputNamesChar.data(), input.data(), input.size(), outputNamesChar.data(), outputNamesChar.size());
#endif
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
    int64_t size = input.size();
    assert(size % mInputShapes[0][1] == 0);
    std::vector<int64_t> inputShape{size / mInputShapes[0][1], mInputShapes[0][1]};
    std::vector<Ort::Value> inputTensors;
#if __has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
    inputTensors.emplace_back(Ort::Experimental::Value::CreateTensor<T>(input.data(), size, inputShape));
#else
    Ort::MemoryInfo mem_info =
      Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);
    inputTensors.emplace_back(Ort::Value::CreateTensor<T>(mem_info, input.data(), size, inputShape.data(), inputShape.size()));
#endif
    LOG(debug) << "Input shape calculated from vector: " << printShape(inputShape);
    return evalModel<T>(inputTensors);
  }

  // Reset session
#if __has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
  void resetSession() { mSession.reset(new Ort::Experimental::Session{*mEnv, modelPath, sessionOptions}); }
#else
  void resetSession() { mSession.reset(new Ort::Session{*mEnv, modelPath.c_str(), sessionOptions}); }
#endif

  // Getters & Setters
  Ort::SessionOptions* getSessionOptions() { return &sessionOptions; } // For optimizations in post
#if __has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
  std::shared_ptr<Ort::Experimental::Session> getSession() { return mSession; }
#else
  std::shared_ptr<Ort::Session> getSession() { return mSession; }
#endif
  int getNumInputNodes() const { return mInputShapes[0][1]; }
  int getNumOutputNodes() const { return mOutputShapes[0][1]; }
  uint64_t getValidityFrom() const { return validFrom; }
  uint64_t getValidityUntil() const { return validUntil; }
  void setActiveThreads(int);

 private:
  // Environment variables for the ONNX runtime
  std::shared_ptr<Ort::Env> mEnv = nullptr;
#if __has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
  std::shared_ptr<Ort::Experimental::Session> mSession = nullptr;
#else
  std::shared_ptr<Ort::Session> mSession = nullptr;
#endif
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
  bool checkHyperloop(bool = true);
};

} // namespace ml

} // namespace o2

#endif // TOOLS_ML_MODEL_H_
