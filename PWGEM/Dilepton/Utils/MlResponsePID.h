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

/// \file MlResponsePID.h
/// \brief Class to compute the ML response for fwdtracks
/// \author Daiki Sekihata <daiki.sekihata@cern.ch>

#ifndef PWGEM_DILEPTON_UTILS_MLRESPONSEPID_H_
#define PWGEM_DILEPTON_UTILS_MLRESPONSEPID_H_

#include "Tools/ML/MlResponse.h"

#include <Framework/Logger.h>

#include <cstdint>
#include <string>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_TRACK(FEATURE) \
  {                             \
    #FEATURE, static_cast<uint8_t>(InputFeaturesPID::FEATURE)}

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER=FEATURE from track
#define CHECK_AND_FILL_TRACK(GETTER)                     \
  case static_cast<uint8_t>(InputFeaturesPID::GETTER): { \
    inputFeature = track.GETTER;                         \
    break;                                               \
  }

namespace o2::analysis
{
// possible input features for ML
enum class InputFeaturesPID : uint8_t {
  tpcInnerParam,
  tpcNClsFound,
  tpcChi2NCl,
  tpcNSigmaEl,
  tofNSigmaEl,
  meanClusterSizeITSobCosTgl,
};

namespace pwgem::dilepton::mlpid
{
struct candidate {
  float tpcInnerParam{0};
  int tpcNClsFound{0};
  float tpcChi2NCl{0};
  float tpcNSigmaEl{0};
  float tofNSigmaEl{0};
  float meanClusterSizeITSobCosTgl{0};
};
} // namespace pwgem::dilepton::mlpid

template <typename TypeOutputScore = float>
class MlResponsePID : public MlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  MlResponsePID() = default;
  /// Default destructor
  virtual ~MlResponsePID() = default;

  template <typename T>
  float return_feature(uint8_t idx, T const& track)
  {
    float inputFeature = 0.;
    switch (idx) {
      CHECK_AND_FILL_TRACK(tpcInnerParam);
      CHECK_AND_FILL_TRACK(tpcNClsFound);
      CHECK_AND_FILL_TRACK(tpcChi2NCl);
      CHECK_AND_FILL_TRACK(tpcNSigmaEl);
      CHECK_AND_FILL_TRACK(tofNSigmaEl);
      CHECK_AND_FILL_TRACK(meanClusterSizeITSobCosTgl);
    }
    return inputFeature;
  }

  /// Method to get the input features vector needed for ML inference
  /// \param track is the single track, \param collision is the collision
  /// \return inputFeatures vector
  template <typename T>
  std::vector<float> getInputFeatures(T const& track)
  {
    std::vector<float> inputFeatures;
    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      float inputFeature = return_feature(idx, track);
      inputFeatures.emplace_back(inputFeature);
    }
    return inputFeatures;
  }

  /// Method to get the value of variable chosen for binning
  /// \param track is the single track, \param collision is the collision
  /// \return binning variable
  template <typename T>
  float getBinningFeature(T const& track)
  {
    return return_feature(mCachedIndexBinning, track);
  }

  void cacheBinningIndex(std::string const& cfgBinningFeature)
  {
    setAvailableInputFeatures();
    if (MlResponse<TypeOutputScore>::mAvailableInputFeatures.count(cfgBinningFeature)) {
      mCachedIndexBinning = MlResponse<TypeOutputScore>::mAvailableInputFeatures[cfgBinningFeature];
    } else {
      LOG(fatal) << "Binning feature " << cfgBinningFeature << " not available! Please check your configurables.";
    }
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_TRACK(tpcInnerParam),
      FILL_MAP_TRACK(tpcNClsFound),
      FILL_MAP_TRACK(tpcChi2NCl),
      FILL_MAP_TRACK(tpcNSigmaEl),
      FILL_MAP_TRACK(tofNSigmaEl),
      FILL_MAP_TRACK(meanClusterSizeITSobCosTgl),
    };
  }

  uint8_t mCachedIndexBinning; // index correspondance between configurable and available input features
};

} // namespace o2::analysis

#undef FILL_MAP_TRACK
#undef CHECK_AND_FILL_TRACK

#endif // PWGEM_DILEPTON_UTILS_MLRESPONSEPID_H_
