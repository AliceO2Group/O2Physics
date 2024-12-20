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
/// \author Fabio Catalano <fabio.catalano@cern.ch>, CERN
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#ifndef TOOLS_ML_MLRESPONSE_H_
#define TOOLS_ML_MLRESPONSE_H_

#if __has_include(<onnxruntime/core/session/onnxruntime_cxx_api.h>)
#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>
#else
#include <onnxruntime_cxx_api.h>
#endif

#include <map>
#include <string>
#include <vector>

#include "CCDB/CcdbApi.h"
#include "Framework/Array2D.h"

#include "Tools/ML/model.h"

namespace o2
{
namespace cuts_ml
{
// direction of the cut
enum CutDirection {
  CutGreater = 0, // require score < cut value
  CutSmaller,     // require score > cut value
  CutNot          // do not cut on score
};
} // namespace cuts_ml

namespace analysis
{
// TypeOutputScore is the type of the output score from o2::ml::OnnxModel (float by default)
template <typename TypeOutputScore = float>
class MlResponse
{
 public:
  /// Default constructor
  MlResponse() = default;
  /// Default destructor
  virtual ~MlResponse() = default;

  /// Configure class instance (import configurables)
  /// \param binsLimits is a vector containing bins limits
  /// \param cuts is a LabeledArray containing selections per bin
  /// \param cutDir is a vector telling whether to reject score values greater or smaller than the threshold
  /// \param nClasses is the number of classes for each model
  void configure(const std::vector<double>& binsLimits, const o2::framework::LabeledArray<double>& cuts, const std::vector<int>& cutDir, const uint8_t& nClasses);

  /// Set model paths to CCDB
  /// \param onnxFiles is a vector of onnx file names, one for each bin
  /// \param ccdbApi is the CCDB API
  /// \param pathsCCDB is a vector of model paths in CCDB, one for each bin
  /// \param timestampCCDB is the CCDB timestamp
  /// \note On the CCDB, different models must be stored in different folders
  void setModelPathsCCDB(const std::vector<std::string>& onnxFiles, const o2::ccdb::CcdbApi& ccdbApi, const std::vector<std::string>& pathsCCDB, int64_t timestampCCDB);

  /// Set model paths to local or cvmfs
  /// \param onnxFiles is a vector of onnx file names, one for each bin
  void setModelPathsLocal(const std::vector<std::string>& onnxFiles);

  /// Initialize class instance (initialize OnnxModels)
  /// \param enableOptimizations is a switch to enable optimizations
  /// \param threads is the number of active threads
  void init(bool enableOptimizations = false, int threads = 0);

  /// Method to translate configurable input-feature strings into integers
  /// \param cfgInputFeatures array of input features names
  void cacheInputFeaturesIndices(std::vector<std::string> const& cfgInputFeatures);

  /// Get vector with model predictions
  /// \param input a vector containing the values of features used in the model
  /// \param nModel is the model index
  /// \return model prediction for each class and the selected model
  template <typename T1, typename T2>
  std::vector<TypeOutputScore> getModelOutput(T1& input, const T2& nModel);

  /// ML selections
  /// \param input is the input features
  /// \param candVar is the variable value (e.g. pT) used to select which model to use
  /// \return boolean telling if model predictions pass the cuts
  template <typename T1, typename T2>
  bool isSelectedMl(T1& input, const T2& candVar);

  /// ML selections
  /// \param input is the input features
  /// \param candVar is the variable value (e.g. pT) used to select which model to use
  /// \param output is a container to be filled with model output
  /// \return boolean telling if model predictions pass the cuts
  template <typename T1, typename T2>
  bool isSelectedMl(T1& input, const T2& candVar, std::vector<TypeOutputScore>& output);

 protected:
  std::vector<o2::ml::OnnxModel> mModels;                 // OnnxModel objects, one for each bin
  uint8_t mNModels = 1;                                   // number of bins
  uint8_t mNClasses = 3;                                  // number of model classes
  std::vector<double> mBinsLimits = {};                   // bin limits of the variable (e.g. pT) used to select which model to use
  std::vector<std::string> mPaths = {""};                 // paths to the models, one for each bin
  std::vector<int> mCutDir = {};                          // direction of the cuts on the model scores (no cut is also supported)
  o2::framework::LabeledArray<double> mCuts = {};         // array of cut values to apply on the model scores
  std::map<std::string, uint8_t> mAvailableInputFeatures; // map of available input features
  std::vector<uint8_t> mCachedIndices;                    // vector of index correspondance between configurables and available input features

  virtual void setAvailableInputFeatures() { return; } // method to fill the map of available input features

 private:
  /// Finds matching bin in mBinsLimits
  /// \param value e.g. pT
  /// \return index of the matching bin, used to access mModels
  /// \note Accounts for the offset due to mBinsLimits storing bin limits (same convention as needed to configure a histogram axis)
  template <typename T>
  int findBin(T const& value)
  {
    if (value < mBinsLimits.front()) {
      return -1;
    }
    if (value >= mBinsLimits.back()) {
      return -1;
    }
    return std::distance(mBinsLimits.begin(), std::upper_bound(mBinsLimits.begin(), mBinsLimits.end(), value)) - 1;
  }
};

extern template class MlResponse<float>;
extern template class MlResponse<double>;

} // namespace analysis
} // namespace o2

#endif // TOOLS_ML_MLRESPONSE_H_
