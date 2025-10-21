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

#include <onnxruntime_cxx_api.h>

#include <cassert>
#include <cstddef>
#include <cstdint>
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

void OnnxModel::setActiveThreads(const int threads)
{
  activeThreads = threads;
  if (!checkHyperloop(false)) {
    sessionOptions.SetIntraOpNumThreads(activeThreads);
  }
}

} // namespace ml

} // namespace o2
