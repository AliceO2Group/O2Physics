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
/// \file     model.cxx
///
/// \author   Christian Sonnabend <christian.sonnabend@cern.ch>
///
/// \brief    A general-purpose class with functions for ONNX model applications
///

#include "Tools/ML/model.h"

#include <Framework/Logger.h>

#include <TSystem.h>

#include <onnxruntime_c_api.h>
#include <onnxruntime_cxx_api.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace o2
{

namespace ml
{

std::string OnnxModel::printShape(const std::vector<int64_t>& v)
{
  if (v.empty()) {
    return "[]";
  }
  std::stringstream ss("");
  for (std::size_t i = 0; i < v.size() - 1; i++)
    ss << v[i] << "x";
  ss << v[v.size() - 1];
  return ss.str();
}

bool OnnxModel::checkHyperloop(const bool verbose)
{
  /// Testing hyperloop core settings
  const char* alienCores = gSystem->Getenv("ALIEN_JDL_CPUCORES");
  const bool alienCoresFound = (alienCores != NULL);
  if (alienCoresFound) {
    if (verbose) {
      LOGP(info, "Hyperloop test/Grid job detected! Number of cores = {}. Setting threads anyway to 1.", alienCores);
    }
    activeThreads = 1;
    sessionOptions.SetIntraOpNumThreads(activeThreads);
  } else {
    if (verbose) {
      LOGP(info, "Not running on Hyperloop.");
    }
  }

  return alienCoresFound;
}

void OnnxModel::initModel(const std::string& localPath, const bool enableOptimizations, const int threads, const uint64_t from, const uint64_t until)
{

  assert(from <= until);

  LOG(info) << "--- ONNX-ML model ---";
  modelPath = localPath;
  activeThreads = threads;

  /// Running on Hyperloop
  if (!checkHyperloop(true)) {
    sessionOptions.SetIntraOpNumThreads(activeThreads);
  }

  /// Enableing optimizations
  if (enableOptimizations) {
    sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED);
  }

  mEnv = std::make_shared<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "onnx-model");
  mSession = std::make_shared<Ort::Session>(*mEnv, modelPath.c_str(), sessionOptions);
  mMemInfo = Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);

  mInputNames.clear();
  mInputShapes.clear();
  mOutputNames.clear();
  mOutputShapes.clear();
  Ort::AllocatorWithDefaultOptions const tmpAllocator;
  for (std::size_t i = 0; i < mSession->GetInputCount(); ++i) {
    mInputNames.push_back(mSession->GetInputNameAllocated(i, tmpAllocator).get());
  }
  for (std::size_t i = 0; i < mSession->GetInputCount(); ++i) {
    mInputShapes.emplace_back(mSession->GetInputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape());
  }
  for (std::size_t i = 0; i < mSession->GetOutputCount(); ++i) {
    mOutputNames.push_back(mSession->GetOutputNameAllocated(i, tmpAllocator).get());
  }
  for (std::size_t i = 0; i < mSession->GetOutputCount(); ++i) {
    mOutputShapes.emplace_back(mSession->GetOutputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape());
  }
  LOG(info) << "Input Nodes:";
  for (std::size_t i = 0; i < mInputNames.size(); i++) {
    LOG(info) << "\t" << mInputNames[i] << " : " << printShape(mInputShapes[i]);
  }

  LOG(info) << "Output Nodes:";
  for (std::size_t i = 0; i < mOutputNames.size(); i++) {
    LOG(info) << "\t" << mOutputNames[i] << " : " << printShape(mOutputShapes[i]);
  }

  validFrom = from;
  validUntil = until;

  LOG(info) << "Model validity - From: " << validFrom << ", Until: " << validUntil;

  LOG(info) << "--- Model initialized! ---";
}

std::vector<int64_t> OnnxModel::inferInputShape(const std::size_t iinput, const int64_t size) const
{
  const std::vector<int64_t>& modelShape = mInputShapes[iinput];

  // Rank-1 input: the whole vector is the tensor
  if (modelShape.size() < 2) {
    return {size};
  }

  // Product of all non-batch dimensions; dynamic dimensions (< 0) cannot be inferred
  int64_t totalSize = 1;
  bool hasDynamicDim = false;
  for (std::size_t idim = 1; idim < modelShape.size(); idim++) {
    if (modelShape[idim] < 0) {
      hasDynamicDim = true;
    } else {
      totalSize *= modelShape[idim];
    }
  }

  if (hasDynamicDim) {
    if (modelShape.size() == 2) {
      // [batch, features] with dynamic feature dimension: interpret the vector as a single sample
      return {1, size};
    }
    LOG(fatal) << "Input " << iinput << " (" << mInputNames[iinput] << ") has dynamic non-batch dimensions (" << printShape(modelShape) << "), the tensor shape cannot be inferred from a flat vector. Please provide std::vector<Ort::Value> inputs instead.";
  }

  if (totalSize <= 0 || size % totalSize != 0) {
    LOG(fatal) << "Size of the input vector (" << size << ") is not a multiple of the model input size (" << totalSize << ") for input " << iinput << " (" << mInputNames[iinput] << ", shape " << printShape(modelShape) << ")";
  }

  std::vector<int64_t> inputShape;
  inputShape.reserve(modelShape.size());
  inputShape.push_back(size / totalSize);
  for (std::size_t idim = 1; idim < modelShape.size(); idim++) {
    inputShape.push_back(modelShape[idim]);
  }
  return inputShape;
}

std::vector<Ort::Value> OnnxModel::evalModelRaw(std::vector<Ort::Value>& input)
{
  if (!mSession) {
    LOG(fatal) << "OnnxModel::evalModel called before initModel()";
  }
  if (input.size() != mInputNames.size()) {
    LOG(fatal) << "Number of input tensors (" << input.size() << ") does not agree with the number of model inputs (" << mInputNames.size() << ")";
  }
  for (std::size_t i = 0; i < input.size(); i++) {
    LOG(debug) << "Input tensor " << i << " shape: " << printShape(input[i].GetTensorTypeAndShapeInfo().GetShape());
  }

  std::vector<const char*> inputNamesChar(mInputNames.size(), nullptr);
  std::transform(std::begin(mInputNames), std::end(mInputNames), std::begin(inputNamesChar),
                 [](const std::string& str) { return str.c_str(); });

  std::vector<const char*> outputNamesChar(mOutputNames.size(), nullptr);
  std::transform(std::begin(mOutputNames), std::end(mOutputNames), std::begin(outputNamesChar),
                 [](const std::string& str) { return str.c_str(); });

  std::vector<Ort::Value> outputTensors;
  try {
    const Ort::RunOptions runOptions;
    outputTensors = mSession->Run(runOptions, inputNamesChar.data(), input.data(), input.size(), outputNamesChar.data(), outputNamesChar.size());
  } catch (const Ort::Exception& exception) {
    LOG(fatal) << "Error running model inference: " << exception.what();
  }

  LOG(debug) << "Number of output tensors: " << outputTensors.size();
  if (outputTensors.size() != mOutputNames.size()) {
    LOG(fatal) << "Number of output tensors: " << outputTensors.size() << " does not agree with the model specified size: " << mOutputNames.size();
  }
  for (std::size_t i = 0; i < outputTensors.size(); i++) {
    const std::vector<int64_t> shape = outputTensors[i].GetTensorTypeAndShapeInfo().GetShape();
    LOG(debug) << "Output tensor " << i << " shape: " << printShape(shape);
    bool shapeOk = (shape.size() == mOutputShapes[i].size());
    for (std::size_t idim = 0; shapeOk && idim < shape.size(); idim++) {
      // dynamic dimensions of the model (< 0) can take any value
      shapeOk = (mOutputShapes[i][idim] < 0) || (shape[idim] == mOutputShapes[i][idim]);
    }
    if (!shapeOk) {
      LOG(fatal) << "Shape of output tensor " << i << " does not agree with model specification! Output: " << printShape(shape) << " model: " << printShape(mOutputShapes[i]);
    }
  }

  return outputTensors;
}

void OnnxModel::setActiveThreads(const int threads)
{
  activeThreads = threads;
  if (!checkHyperloop(false)) {
    sessionOptions.SetIntraOpNumThreads(activeThreads);
  }
}

} // namespace ml

} // namespace o2
