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

/// \file MlResponseSCT.h
/// \brief Class to compute the ML response for fwdtracks
/// \author Daiki Sekihata <daiki.sekihata@cern.ch>

#ifndef PWGEM_DILEPTON_UTILS_MLRESPONSESCT_H_
#define PWGEM_DILEPTON_UTILS_MLRESPONSESCT_H_

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
    #FEATURE, static_cast<uint8_t>(InputFeaturesSCT::FEATURE)}

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER=FEATURE from track
#define CHECK_AND_FILL_TRACK(GETTER)                     \
  case static_cast<uint8_t>(InputFeaturesSCT::GETTER): { \
    inputFeature = track.GETTER;                         \
    break;                                               \
  }

namespace o2::analysis
{
// possible input features for ML
enum class InputFeaturesSCT : uint8_t {
  ptL,
  impParXYL,
  impParZL,
  impPar3DL,
  impParXYLinSigma,
  impParZLinSigma,
  impPar3DLinSigma,
  ptH,
  massH,
  tpcNSigmaKa,
  impParXYH,
  impParZH,
  impPar3DH,
  impParXYHinSigma,
  impParZHinSigma,
  impPar3DHinSigma,
  signLH,
  dcaLH,
  logChi2PCA,
  massLH,
  signedMassLH,
  missingPtNuPerpToFD,
  correctedMass,
  cpa,
  cpaXY,
  impParXY,
  impParZ,
  impPar3D,
  impParXYinSigma,
  impParZinSigma,
  impPar3DinSigma,
  decayLengthXY,
  decayLengthZ,
  decayLength3D,
  decayLengthXYinSigma,
  decayLengthZinSigma,
  decayLength3DinSigma,
};

namespace pwgem::dilepton::sct
{
struct candidate {
  // lepton info. (don't use them to avoid bias to single lepton)
  float ptL{0};
  float impParXYL{0};
  float impParZL{0};
  float impPar3DL{0};
  float impParXYLinSigma{0};
  float impParZLinSigma{0};
  float impPar3DLinSigma{0};

  // hadron information
  float ptH{0};
  float massH{0}; // only for V0s and Cascades
  float tpcNSigmaKa{0};
  float impParXYH{0};
  float impParZH{0};
  float impPar3DH{0};
  float impParXYHinSigma{0};
  float impParZHinSigma{0};
  float impPar3DHinSigma{0};

  // LH pair information
  int signLH{0};
  float dcaLH{0};
  float logChi2PCA{0};
  float massLH{0};
  float signedMassLH{0};
  float missingPtNuPerpToFD{0};
  float correctedMass{0};
  float cpa{0};
  float cpaXY{0};
  float impParXY{0};
  float impParZ{0};
  float impPar3D{0};
  float impParXYinSigma{0};
  float impParZinSigma{0};
  float impPar3DinSigma{0};
  float decayLengthXY{0};
  float decayLengthZ{0};
  float decayLength3D{0};
  float decayLengthXYinSigma{0};
  float decayLengthZinSigma{0};
  float decayLength3DinSigma{0};
};
} // namespace pwgem::dilepton::sct

template <typename TypeOutputScore = float>
class MlResponseSCT : public MlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  MlResponseSCT() = default;
  /// Default destructor
  virtual ~MlResponseSCT() = default;

  template <typename T>
  float return_feature(uint8_t idx, T const& track)
  {
    float inputFeature = 0.;
    switch (idx) {
      CHECK_AND_FILL_TRACK(ptL);
      CHECK_AND_FILL_TRACK(impParXYL);
      CHECK_AND_FILL_TRACK(impParZL);
      CHECK_AND_FILL_TRACK(impPar3DL);
      CHECK_AND_FILL_TRACK(impParXYLinSigma);
      CHECK_AND_FILL_TRACK(impParZLinSigma);
      CHECK_AND_FILL_TRACK(impPar3DLinSigma);
      CHECK_AND_FILL_TRACK(ptH);
      CHECK_AND_FILL_TRACK(massH);
      CHECK_AND_FILL_TRACK(tpcNSigmaKa);
      CHECK_AND_FILL_TRACK(impParXYH);
      CHECK_AND_FILL_TRACK(impParZH);
      CHECK_AND_FILL_TRACK(impPar3DH);
      CHECK_AND_FILL_TRACK(impParXYHinSigma);
      CHECK_AND_FILL_TRACK(impParZHinSigma);
      CHECK_AND_FILL_TRACK(impPar3DHinSigma);
      CHECK_AND_FILL_TRACK(signLH);
      CHECK_AND_FILL_TRACK(dcaLH);
      CHECK_AND_FILL_TRACK(logChi2PCA);
      CHECK_AND_FILL_TRACK(massLH);
      CHECK_AND_FILL_TRACK(signedMassLH);
      CHECK_AND_FILL_TRACK(missingPtNuPerpToFD);
      CHECK_AND_FILL_TRACK(correctedMass);
      CHECK_AND_FILL_TRACK(cpa);
      CHECK_AND_FILL_TRACK(cpaXY);
      CHECK_AND_FILL_TRACK(impParXY);
      CHECK_AND_FILL_TRACK(impParZ);
      CHECK_AND_FILL_TRACK(impPar3D);
      CHECK_AND_FILL_TRACK(impParXYinSigma);
      CHECK_AND_FILL_TRACK(impParZinSigma);
      CHECK_AND_FILL_TRACK(impPar3DinSigma);
      CHECK_AND_FILL_TRACK(decayLengthXY);
      CHECK_AND_FILL_TRACK(decayLengthZ);
      CHECK_AND_FILL_TRACK(decayLength3D);
      CHECK_AND_FILL_TRACK(decayLengthXYinSigma);
      CHECK_AND_FILL_TRACK(decayLengthZinSigma);
      CHECK_AND_FILL_TRACK(decayLength3DinSigma);
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
      FILL_MAP_TRACK(ptL),
      FILL_MAP_TRACK(impParXYL),
      FILL_MAP_TRACK(impParZL),
      FILL_MAP_TRACK(impPar3DL),
      FILL_MAP_TRACK(impParXYLinSigma),
      FILL_MAP_TRACK(impParZLinSigma),
      FILL_MAP_TRACK(impPar3DLinSigma),
      FILL_MAP_TRACK(ptH),
      FILL_MAP_TRACK(massH),
      FILL_MAP_TRACK(tpcNSigmaKa),
      FILL_MAP_TRACK(impParXYH),
      FILL_MAP_TRACK(impParZH),
      FILL_MAP_TRACK(impPar3DH),
      FILL_MAP_TRACK(impParXYHinSigma),
      FILL_MAP_TRACK(impParZHinSigma),
      FILL_MAP_TRACK(impPar3DHinSigma),
      FILL_MAP_TRACK(signLH),
      FILL_MAP_TRACK(dcaLH),
      FILL_MAP_TRACK(logChi2PCA),
      FILL_MAP_TRACK(massLH),
      FILL_MAP_TRACK(signedMassLH),
      FILL_MAP_TRACK(missingPtNuPerpToFD),
      FILL_MAP_TRACK(correctedMass),
      FILL_MAP_TRACK(cpa),
      FILL_MAP_TRACK(cpaXY),
      FILL_MAP_TRACK(impParXY),
      FILL_MAP_TRACK(impParZ),
      FILL_MAP_TRACK(impPar3D),
      FILL_MAP_TRACK(impParXYinSigma),
      FILL_MAP_TRACK(impParZinSigma),
      FILL_MAP_TRACK(impPar3DinSigma),
      FILL_MAP_TRACK(decayLengthXY),
      FILL_MAP_TRACK(decayLengthZ),
      FILL_MAP_TRACK(decayLength3D),
      FILL_MAP_TRACK(decayLengthXYinSigma),
      FILL_MAP_TRACK(decayLengthZinSigma),
      FILL_MAP_TRACK(decayLength3DinSigma),
    };
  }

  uint8_t mCachedIndexBinning; // index correspondance between configurable and available input features
};

} // namespace o2::analysis

#undef FILL_MAP_TRACK
#undef CHECK_AND_FILL_TRACK

#endif // PWGEM_DILEPTON_UTILS_MLRESPONSESCT_H_
