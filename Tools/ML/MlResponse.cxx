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

/// \file MlResponse.h
/// \brief Class to compute the ML response for analysis selections
#include "Tools/ML/MlResponse.h"
#include "Tools/ML/model.h"
#include "CCDB/CcdbApi.h"

namespace o2::analysis
{

template <typename TypeOutputScore>
void MlResponse<TypeOutputScore>::configure(const std::vector<double>& binsLimits, const o2::framework::LabeledArray<double>& cuts, const std::vector<int>& cutDir, const uint8_t& nClasses)
{
  if (cutDir.size() != nClasses) {
    LOG(fatal) << "Number of classes (" << static_cast<int>(nClasses) << ") different from the number of cuts on model scores (" << cutDir.size() << ")! Please check your configurables.";
  }

  this->mBinsLimits = binsLimits;
  mCuts = cuts;
  mCutDir = cutDir;
  mNClasses = nClasses;
  mNModels = binsLimits.size() - 1;
  this->mModels = std::vector<o2::ml::OnnxModel>(mNModels);
  mPaths = std::vector<std::string>(mNModels);
}

/// Set model paths to CCDB
/// \param onnxFiles is a vector of onnx file names, one for each bin
/// \param ccdbApi is the CCDB API
/// \param pathsCCDB is a vector of model paths in CCDB, one for each bin
/// \param timestampCCDB is the CCDB timestamp
/// \note On the CCDB, different models must be stored in different folders
template <typename TypeOutputScore>
void MlResponse<TypeOutputScore>::setModelPathsCCDB(const std::vector<std::string>& onnxFiles, const o2::ccdb::CcdbApi& ccdbApi, const std::vector<std::string>& pathsCCDB, int64_t timestampCCDB)
{
  if (onnxFiles.size() != mNModels) {
    LOG(fatal) << "Number of expected models (" << mNModels << ") different from the one set (" << onnxFiles.size() << ")! Please check your configurables.";
  }
  if (pathsCCDB.size() != mNModels) {
    LOG(fatal) << "Number of expected models (" << mNModels << ") different from the number of CCDB paths (" << pathsCCDB.size() << ")! Please check your configurables.";
  }

  // check that the path is unique for each BDT model (otherwise CCDB download does not work as expected)
  for (auto iThisFile{0}; iThisFile < mNModels; ++iThisFile) {
    for (auto iOtherFile{iThisFile + 1}; iOtherFile < mNModels; ++iOtherFile) {
      if ((pathsCCDB[iThisFile] == pathsCCDB[iOtherFile]) && (onnxFiles[iThisFile] != onnxFiles[iOtherFile])) {
        LOGP(fatal, "More than one model ({} and {}) in the same CCDB directory ({})! Each directory in CCDB can contain only one model. Please check your configurables.", onnxFiles[iThisFile], onnxFiles[iOtherFile], pathsCCDB[iThisFile]);
      }
    }
  }

  for (auto iFile{0}; iFile < mNModels; ++iFile) {
    std::map<std::string, std::string> metadata;
    bool retrieveSuccess = ccdbApi.retrieveBlob(pathsCCDB[iFile], ".", metadata, timestampCCDB, false, onnxFiles[iFile]);
    if (retrieveSuccess) {
      mPaths[iFile] = onnxFiles[iFile];
    } else {
      LOG(fatal) << "Error encountered while accessing the ML model from " << pathsCCDB[iFile] << "! Maybe the ML model doesn't exist yet for this run number or timestamp?";
    }
  }
}
/// Initialize class instance (initialize OnnxModels)
/// \param enableOptimizations is a switch to enable optimizations
/// \param threads is the number of active threads
template <typename TypeOutputScore>
void MlResponse<TypeOutputScore>::init(bool enableOptimizations, int threads)
{
  uint8_t counterModel{0};
  for (const auto& path : mPaths) {
    this->mModels[counterModel].initModel(path, enableOptimizations, threads);
    ++counterModel;
  }
}

template <typename TypeOutputScore>
template <typename T1, typename T2>
std::vector<TypeOutputScore> MlResponse<TypeOutputScore>::getModelOutput(T1& input, const T2& nModel)
{
  if (nModel < 0 || static_cast<std::size_t>(nModel) >= this->mModels.size()) {
    LOG(fatal) << "Model index " << nModel << " is out of range! The number of initialised models is " << this->mModels.size() << ". Please check your configurables.";
  }

  TypeOutputScore* outputPtr = this->mModels[nModel].evalModel(input);
  return std::vector<TypeOutputScore>{outputPtr, outputPtr + mNClasses};
}

template <typename TypeOutputScore>
template <typename T1, typename T2>
bool MlResponse<TypeOutputScore>::isSelectedMl(T1& input, const T2& candVar)
{
  int nModel = findBin(candVar);
  auto output = getModelOutput(input, nModel);
  uint8_t iClass{0};
  for (const auto& outputValue : output) {
    uint8_t dir = mCutDir.at(iClass);
    if (dir != o2::cuts_ml::CutDirection::CutNot) {
      if (dir == o2::cuts_ml::CutDirection::CutGreater && outputValue > mCuts.get(nModel, iClass)) {
        return false;
      }
      if (dir == o2::cuts_ml::CutDirection::CutSmaller && outputValue < mCuts.get(nModel, iClass)) {
        return false;
      }
    }
    ++iClass;
  }
  return true;
}

template <typename TypeOutputScore>
template <typename T1, typename T2>
bool MlResponse<TypeOutputScore>::isSelectedMl(T1& input, const T2& candVar, std::vector<TypeOutputScore>& output)
{
  int nModel = findBin(candVar);
  output = getModelOutput(input, nModel);
  uint8_t iClass{0};
  for (const auto& outputValue : output) {
    uint8_t dir = mCutDir.at(iClass);
    if (dir != o2::cuts_ml::CutDirection::CutNot) {
      if (dir == o2::cuts_ml::CutDirection::CutGreater && outputValue > mCuts.get(nModel, iClass)) {
        return false;
      }
      if (dir == o2::cuts_ml::CutDirection::CutSmaller && outputValue < mCuts.get(nModel, iClass)) {
        return false;
      }
    }
    ++iClass;
  }
  return true;
}

template <typename TypeOutputScore>
void MlResponse<TypeOutputScore>::cacheInputFeaturesIndices(std::vector<std::string> const& cfgInputFeatures)
{
  setAvailableInputFeatures();
  for (const auto& inputFeature : cfgInputFeatures) {
    if (mAvailableInputFeatures.count(inputFeature)) {
      mCachedIndices.emplace_back(mAvailableInputFeatures[inputFeature]);
    } else {
      LOG(fatal) << "Input feature " << inputFeature << " not available! Please check your configurables.";
    }
  }
}

template <typename TypeOutputScore>
void MlResponse<TypeOutputScore>::setModelPathsLocal(const std::vector<std::string>& onnxFiles)
{
  if (onnxFiles.size() != mNModels) {
    LOG(fatal) << "Number of expected models (" << mNModels << ") different from the one set (" << onnxFiles.size() << ")! Please check your configurables.";
  }
  mPaths = onnxFiles;
}

template bool MlResponse<float>::isSelectedMl<std::vector<float>, float>(std::vector<float>&, const float&, std::vector<float>&);
template class MlResponse<float>;
template class MlResponse<double>;
} // namespace o2::analysis
