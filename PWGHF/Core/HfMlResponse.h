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

/// \file HFMLResponse.h
/// \brief Class to compute the ML response for HF-analysis selections
/// \author Fabio Catalano <fabio.catalano@cern.ch>, CERN
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#ifndef PWGHF_CORE_HFMLRESPONSE_H_
#define PWGHF_CORE_HFMLRESPONSE_H_

#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>

#include <map>
#include <string>
#include <vector>

#include "CCDB/CcdbApi.h"
#include "Framework/Array2D.h"

#include "Tools/ML/model.h"

#include "PWGHF/Core/SelectorCuts.h"

namespace o2::analysis
{

enum PIDStatus : uint8_t {
  TpcAndTof = 0,
  TpcOnly,
  TofOnly
};

/// Tag track PID status (TPC and TOF, only TPC, only TOF)
/// \param track is the track
/// \return tag of PID status of the track
template <typename T1>
uint8_t tagPID(const T1& track)
{
  uint8_t tag{0};
  if (track.hasTPC() && track.hasTOF()) {
    SETBIT(tag, PIDStatus::TpcAndTof);
  } else if (track.hasTPC() && !track.hasTOF()) {
    SETBIT(tag, PIDStatus::TpcOnly);
  } else if (!track.hasTPC() && track.hasTOF()) {
    SETBIT(tag, PIDStatus::TofOnly);
  }
  return tag;
}

/// Combine TPC and TOF nSigma
/// \param nSigTpc is the number of sigmas for TPC
/// \param nSigTof is the number of sigmas for TOF
/// \param tagPID is the tag on PID status of the track
/// \return combined nSigma value
template <typename T1, typename T2>
T1 getCombinedNSigma(const T1& nSigTpc, const T1& nSigTof, const T2& tagPID)
{
  if (TESTBIT(tagPID, PIDStatus::TpcAndTof)) {
    return std::sqrt((nSigTpc * nSigTpc + nSigTof * nSigTof) / 2);
  } else if (TESTBIT(tagPID, PIDStatus::TpcOnly)) {
    return std::abs(nSigTpc);
  } else if (TESTBIT(tagPID, PIDStatus::TofOnly)) {
    return std::abs(nSigTof);
  }

  return 0.;
}

template <typename T = float>
class HfMlResponse
{
 public:
  /// Default constructor
  HfMlResponse() = default;
  /// Default destructor
  ~HfMlResponse() = default;

  /// Configure class instance (import configurables)
  /// \param binsLimits is a vector containing bins limits
  /// \param cuts is a LabeledArray containing selections per bin
  /// \param cutDir is a vector telling whether to reject score values greater or smaller than the threshold
  /// \param nClasses is the number of classes for each model
  void configure(const std::vector<double>& binsLimits, const o2::framework::LabeledArray<double>& cuts, const std::vector<int>& cutDir, const uint8_t& nClasses)
  {
    mBinsLimits = binsLimits;
    mCuts = cuts;
    mCutDir = cutDir;
    mNClasses = nClasses;
    mNModels = binsLimits.size() - 1;
    mModels = std::vector<o2::ml::OnnxModel>(mNModels);
    mPaths = std::vector<std::string>(mNModels);
  }

  /// Set models paths to CCDB
  /// \param onnxFiles is a vector of onnx file names, one for each bin
  /// \param ccdbApi is the CCDB API
  /// \param pathCCDB is the model path in CCDB
  /// \param timestampCCDB is the CCDB timestamp
  void setModelPathsCCDB(const std::vector<std::string>& onnxFiles, const o2::ccdb::CcdbApi& ccdbApi, std::string pathCCDB, int64_t timestampCCDB)
  {
    if (onnxFiles.size() != mNModels) {
      LOG(fatal) << "Number of expected models different from the one set! Please check your configurables.";
    }

    uint8_t counterModel{0};
    for (const auto& onnxFile : onnxFiles) {
      std::map<std::string, std::string> metadata;
      bool retrieveSuccess = ccdbApi.retrieveBlob(pathCCDB, ".", metadata, timestampCCDB, false, onnxFile);
      if (retrieveSuccess) {
        mPaths[counterModel] = onnxFile;
      } else {
        LOG(fatal) << "Error encountered while accessing the ML model from CCDB! Maybe the ML model doesn't exist yet for this runnumber/timestamp?";
      }
      ++counterModel;
    }
  }

  /// Set models paths to local or cvmfs
  /// \param onnxFiles is a vector of onnx file names, one for each bin
  void setModelPathsLocal(const std::vector<std::string>& onnxFiles)
  {
    if (onnxFiles.size() != mNModels) {
      LOG(fatal) << "Number of expected models different from the one set! Please check your configurables.";
    }

    mPaths = onnxFiles;
  }

  /// Initialize class instance (initialize OnnxModels)
  /// \param enableOptimizations is a switch no enable optimizations
  /// \param threads is the number of active threads
  void init(bool enableOptimizations = false, int threads = 0)
  {
    uint8_t counterModel{0};
    for (const auto& path : mPaths) {
      mModels[counterModel].initModel(path, enableOptimizations, threads);
      ++counterModel;
    }
  }

  /// Get vector with model predictions
  /// \param input a vector containing the values of features used in the model
  /// \param nModel is the model index
  /// \return model prediction for each class and the selected model
  template <typename T1, typename T2>
  std::vector<T> getModelOutput(T1& input, const T2& nModel)
  {
    T* outputPtr = mModels[nModel].evalModel(input);
    return std::vector<T>{outputPtr, outputPtr + mNClasses};
  }

  /// ML selections
  /// \param input is the input features
  /// \param pt is the candidate transverse momentum
  /// \return boolean telling if model predictions pass the cuts
  template <typename T1, typename T2>
  bool isSelectedMl(T1& input, const T2& pt)
  {
    auto nModel = o2::analysis::findBin(&mBinsLimits, pt);
    auto output = getModelOutput(input, nModel);
    uint8_t iClass{0};
    for (const auto& outputValue : output) {
      uint8_t dir = mCutDir.at(iClass);
      if (dir != hf_cuts_ml::CutDirection::CutNot) {
        if (dir == hf_cuts_ml::CutDirection::CutGreater && outputValue > mCuts.get(nModel, iClass)) {
          return false;
        }
        if (dir == hf_cuts_ml::CutDirection::CutSmaller && outputValue < mCuts.get(nModel, iClass)) {
          return false;
        }
      }
      ++iClass;
    }
    return true;
  }

  /// ML selections
  /// \param input is the input features
  /// \param pt is the candidate transverse momentum
  /// \param output is a container to be filled with model output
  /// \return boolean telling if model predictions pass the cuts
  template <typename T1, typename T2>
  bool isSelectedMl(T1& input, const T2& pt, std::vector<T>& output)
  {
    auto nModel = o2::analysis::findBin(&mBinsLimits, pt);
    output = getModelOutput(input, nModel);
    uint8_t iClass{0};
    for (const auto& outputValue : output) {
      uint8_t dir = mCutDir.at(iClass);
      if (dir != hf_cuts_ml::CutDirection::CutNot) {
        if (dir == hf_cuts_ml::CutDirection::CutGreater && outputValue > mCuts.get(nModel, iClass)) {
          return false;
        }
        if (dir == hf_cuts_ml::CutDirection::CutSmaller && outputValue < mCuts.get(nModel, iClass)) {
          return false;
        }
      }
      ++iClass;
    }
    return true;
  }

 private:
  std::vector<o2::ml::OnnxModel> mModels;         // OnnxModel objects, one for each bin
  uint8_t mNModels = 1;                           // number of bins
  uint8_t mNClasses = 3;                          // number of model classes
  std::vector<double> mBinsLimits = {};           // bin limits of the variable (e.g. pT) used to select which model to use
  std::vector<std::string> mPaths = {""};         // paths to the models, one for each bin
  std::vector<int> mCutDir = {};                  // direction of the cuts on the model scores (no cut is also supported)
  o2::framework::LabeledArray<double> mCuts = {}; // array of cut values to apply on the model scores
};

} // namespace o2::analysis

#endif // PWGHF_CORE_HFMLRESPONSE_H_
