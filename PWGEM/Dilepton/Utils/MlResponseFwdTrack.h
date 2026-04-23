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

/// \file MlResponseFwdTrack.h
/// \brief Class to compute the ML response for fwdtracks
/// \author Daiki Sekihata <daiki.sekihata@cern.ch>

#ifndef PWGEM_DILEPTON_UTILS_MLRESPONSEFWDTRACK_H_
#define PWGEM_DILEPTON_UTILS_MLRESPONSEFWDTRACK_H_

#include "Tools/ML/MlResponse.h"

#include <Framework/Logger.h>

#include <cstdint>
#include <string>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_TRACK(FEATURE)                                \
  {                                                               \
    #FEATURE, static_cast<uint8_t>(InputFeaturesFwdTrack::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER=FEATURE from track
#define CHECK_AND_FILL_TRACK(GETTER)                      \
  case static_cast<uint8_t>(InputFeaturesFwdTrack::GETTER): { \
    inputFeature = track.GETTER;                             \
    break;                                                   \
  }

namespace o2::analysis
{
// possible input features for ML
enum class InputFeaturesFwdTrack : uint8_t {
    multFT0C,
    multMFT,
    ptMCHMID,
    rSigned1Pt,
    dEta,
    dPhi,
    dX,
    dY,
    chi2MatchMCHMFT,
    sigmaPhiMFT,
    sigmaTglMFT,
    sigmaPhiMCHMID,
    sigmaTglMCHMID,
};

template <typename TypeOutputScore = float>
class MlResponseFwdTrack : public MlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  MlResponseFwdTrack() = default;
  /// Default destructor
  virtual ~MlResponseFwdTrack() = default;

  template <typename T>
  float return_feature(uint8_t idx, T const& track)
  {
    float inputFeature = 0.;
    switch (idx) {
      CHECK_AND_FILL_TRACK(multFT0C);
      CHECK_AND_FILL_TRACK(multMFT);
      CHECK_AND_FILL_TRACK(ptMCHMID);
      CHECK_AND_FILL_TRACK(rSigned1Pt);
      CHECK_AND_FILL_TRACK(dEta);
      CHECK_AND_FILL_TRACK(dPhi);
      CHECK_AND_FILL_TRACK(dX);
      CHECK_AND_FILL_TRACK(dY);
      CHECK_AND_FILL_TRACK(chi2MatchMCHMFT);
      CHECK_AND_FILL_TRACK(sigmaPhiMFT);
      CHECK_AND_FILL_TRACK(sigmaTglMFT);
      CHECK_AND_FILL_TRACK(sigmaPhiMCHMID);
      CHECK_AND_FILL_TRACK(sigmaTglMCHMID);
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
      FILL_MAP_TRACK(multFT0C),
      FILL_MAP_TRACK(multMFT),
      FILL_MAP_TRACK(ptMCHMID),
      FILL_MAP_TRACK(rSigned1Pt),
      FILL_MAP_TRACK(dEta),
      FILL_MAP_TRACK(dPhi),
      FILL_MAP_TRACK(dX),
      FILL_MAP_TRACK(dY),
      FILL_MAP_TRACK(chi2MatchMCHMFT),
      FILL_MAP_TRACK(sigmaPhiMFT),
      FILL_MAP_TRACK(sigmaTglMFT),
      FILL_MAP_TRACK(sigmaPhiMCHMID),
      FILL_MAP_TRACK(sigmaTglMCHMID),
    };
  }

  uint8_t mCachedIndexBinning; // index correspondance between configurable and available input features
};

} // namespace o2::analysis

#undef FILL_MAP_TRACK
#undef CHECK_AND_FILL_TRACK

#endif // PWGEM_DILEPTON_UTILS_MLRESPONSEFWDTRACK_H_
