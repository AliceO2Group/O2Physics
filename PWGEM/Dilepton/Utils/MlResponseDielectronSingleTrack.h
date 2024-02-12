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

/// \file MlResponseDiletponSingleTrack.h
/// \brief Class to compute the ML response for dielectron analyses at the single track level
/// \author Daniel Samitz <daniel.samitz@cern.ch>, SMI Vienna
///         Elisa Meninno, <elisa.meninno@cern.ch>, SMI Vienna

#ifndef PWGEM_DILEPTON_UTILS_MLRESPONSEDIELECTRONSINGLETRACK_H_
#define PWGEM_DILEPTON_UTILS_MLRESPONSEDIELECTRONSINGLETRACK_H_

#include <map>
#include <string>
#include <vector>

#include "Tools/ML/MlResponse.h"

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_DIELECTRON_SINGLE_TRACK(FEATURE)                                 \
  {                                                                               \
#FEATURE, static_cast < uint8_t>(InputFeaturesDielectronSingleTrack::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER=FEATURE from track
#define CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(GETTER)                 \
  case static_cast<uint8_t>(InputFeaturesDielectronSingleTrack::GETTER): { \
    inputFeatures.emplace_back(track.GETTER());                            \
    break;                                                                 \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER form track and applying a sqrt
#define CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK_SQRT(FEATURE, GETTER)    \
  case static_cast<uint8_t>(InputFeaturesDielectronSingleTrack::FEATURE): { \
    inputFeatures.emplace_back(sqrt(track.GETTER()));                       \
    break;                                                                  \
  }

namespace o2::analysis
{
// possible input features for ML
enum class InputFeaturesDielectronSingleTrack : uint8_t {
  sign,
  pt,
  eta,
  phi,
  dcaXY,
  dcaZ,
  dcaResXY,
  dcaResZ,
  tpcNClsFindable,
  tpcNClsFound,
  tpcNClsCrossedRows,
  tpcChi2NCl,
  tpcInnerParam,
  tpcSignal,
  tpcNSigmaEl,
  tpcNSigmaMu,
  tpcNSigmaPi,
  tpcNSigmaKa,
  tpcNSigmaPr,
  beta,
  tofNSigmaEl,
  tofNSigmaMu,
  tofNSigmaPi,
  tofNSigmaKa,
  tofNSigmaPr,
  itsClusterMap,
  itsChi2NCl
};

template <typename TypeOutputScore = float>
class MlResponseDielectronSingleTrack : public MlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  MlResponseDielectronSingleTrack() = default;
  /// Default destructor
  virtual ~MlResponseDielectronSingleTrack() = default;

  /// Method to get the input features vector needed for ML inference
  /// \param track is the single track
  /// \return inputFeatures vector
  template <typename T>
  std::vector<float> getInputFeatures(T const& track)
  {
    std::vector<float> inputFeatures;

    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      switch (idx) {
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(sign);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(pt);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(eta);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(phi);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(dcaXY);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(dcaZ);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK_SQRT(dcaResXY, cYY);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK_SQRT(dcaResZ, cZZ);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(tpcNClsFindable);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(tpcNClsFound);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(tpcNClsCrossedRows);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(tpcChi2NCl);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(tpcInnerParam);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(tpcSignal);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(tpcNSigmaEl);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(tpcNSigmaMu);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(tpcNSigmaPi);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(tpcNSigmaKa);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(tpcNSigmaPr);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(beta);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(tofNSigmaEl);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(tofNSigmaMu);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(tofNSigmaPi);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(tofNSigmaKa);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(tofNSigmaPr);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(itsClusterMap);
        CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK(itsChi2NCl);
      }
    }

    return inputFeatures;
  }

 protected:
  /// Method to fill the map of available input features
  void setAvailableInputFeatures()
  {
    MlResponse<TypeOutputScore>::mAvailableInputFeatures = {
      FILL_MAP_DIELECTRON_SINGLE_TRACK(sign),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(pt),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(eta),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(phi),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(dcaXY),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(dcaZ),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(dcaResXY),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(dcaResZ),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcNClsFindable),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcNClsFound),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcNClsCrossedRows),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcChi2NCl),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcInnerParam),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcSignal),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcNSigmaEl),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcNSigmaMu),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcNSigmaPi),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcNSigmaKa),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcNSigmaPr),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(beta),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tofNSigmaEl),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tofNSigmaMu),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tofNSigmaPi),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tofNSigmaKa),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tofNSigmaPr),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(itsClusterMap),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(itsChi2NCl)};
  }
};

} // namespace o2::analysis

#undef FILL_MAP_DIELECTRON_SINGLE_TRACK
#undef CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK
#undef CHECK_AND_FILL_VEC_DIELECTRON_SINGLE_TRACK_SQRT

#endif // PWGEM_DILEPTON_UTILS_MLRESPONSEDIELECTRONSINGLETRACK_H_
