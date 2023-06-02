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

// ONNX includes
#include "Tools/ML/model.h"

namespace o2
{

namespace ml
{

std::string OnnxModel::printShape(const std::vector<int64_t>& v)
{
  std::stringstream ss("");
  for (size_t i = 0; i < v.size() - 1; i++)
    ss << v[i] << "x";
  ss << v[v.size() - 1];
  return ss.str();
}

bool OnnxModel::checkHyperloop(bool verbose)
{
  /// Testing hyperloop core settings
  const char* alienCores = gSystem->Getenv("ALIEN_JDL_CPUCORES");
  bool alienCoresFound = (alienCores != NULL);
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

/// Access model from CCDB
/// \param onnxFile is the onnx file name
/// \param ccdbApi is the CCDB API
/// \param pathCCDB is the model path in CCDB
/// \param timestampCCDB is the CCDB timestamp
void OnnxModel::accessModelFromCCDB(std::string onnxFile, o2::ccdb::CcdbApi& ccdbApi, std::string pathCCDB,long timestampCCDB) {
  std::map<std::string, std::string> metadata;
  bool retrieveSuccess = ccdbApi.retrieveBlob(pathCCDB, ".", metadata, timestampCCDB, false, onnxFile);
  if (retrieveSuccess) {
    modelPath = onnxFile;
  } else {
    LOG(fatal) << "Error encountered while accessing the ML model from CCDB! Maybe the ML model doesn't exist yet for this runnumber/timestamp?";
  }
}

void OnnxModel::initModel(std::string localPath, bool enableOptimizations, int threads, uint64_t from, uint64_t until)
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
  mSession = std::make_shared<Ort::Experimental::Session>(*mEnv, modelPath, sessionOptions);

  mInputNames = mSession->GetInputNames();
  mInputShapes = mSession->GetInputShapes();
  mOutputNames = mSession->GetOutputNames();
  mOutputShapes = mSession->GetOutputShapes();

  LOG(info) << "Input Nodes:";
  for (size_t i = 0; i < mInputNames.size(); i++) {
    LOG(info) << "\t" << mInputNames[i] << " : " << printShape(mInputShapes[i]);
  }

  LOG(info) << "Output Nodes:";
  for (size_t i = 0; i < mOutputNames.size(); i++) {
    LOG(info) << "\t" << mOutputNames[i] << " : " << printShape(mOutputShapes[i]);
  }

  validFrom = from;
  validUntil = until;

  LOG(info) << "Model validity - From: " << validFrom << ", Until: " << validUntil;

  LOG(info) << "--- Model initialized! ---";
}

void OnnxModel::initModel(bool enableOptimizations, int threads, uint64_t from, uint64_t until)
{
  initModel(modelPath, enableOptimizations, threads, from, until);
}

void OnnxModel::setActiveThreads(int threads)
{
  activeThreads = threads;
  if (!checkHyperloop(false)) {
    sessionOptions.SetIntraOpNumThreads(activeThreads);
  }
}

} // namespace ml

} // namespace o2
