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

/// \file MlResponseO2Track.h
/// \brief Class to compute the ML response for dielectron analyses at the single track level
/// \author Daniel Samitz <daniel.samitz@cern.ch>, SMI Vienna
///         Elisa Meninno, <elisa.meninno@cern.ch>, SMI Vienna

#ifndef PWGEM_DILEPTON_UTILS_MLRESPONSEO2TRACK_H_
#define PWGEM_DILEPTON_UTILS_MLRESPONSEO2TRACK_H_

#include "Tools/ML/MlResponse.h"

#include <map>
#include <string>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_O2_TRACK(FEATURE)                                \
  {                                                               \
    #FEATURE, static_cast<uint8_t>(InputFeaturesO2Track::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER=FEATURE from track
#define CHECK_AND_FILL_O2_TRACK(GETTER)                      \
  case static_cast<uint8_t>(InputFeaturesO2Track::GETTER): { \
    inputFeature = track.GETTER();                           \
    break;                                                   \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER=FEATURE from track
#define CHECK_AND_FILL_O2_TRACKPARCOV(FEATURE, GETTER)        \
  case static_cast<uint8_t>(InputFeaturesO2Track::FEATURE): { \
    inputFeature = trackParCov.GETTER();                      \
    break;                                                    \
  }

// mean ITS cluster size
#define CHECK_AND_FILL_O2_TRACK_MEAN_ITSCLUSTER_SIZE(FEATURE, v1, v2)    \
  case static_cast<uint8_t>(InputFeaturesO2Track::FEATURE): {            \
    int nsize = 0;                                                       \
    int ncls = 0;                                                        \
    for (int il = v1; il < v2; il++) {                                   \
      nsize += track.itsClsSizeInLayer(il);                              \
      if (nsize > 0) {                                                   \
        ncls++;                                                          \
      }                                                                  \
    }                                                                    \
    inputFeature = static_cast<float>(nsize) / static_cast<float>(ncls); \
    break;                                                               \
  }

// mean ITS cluster size x cos(lambda)
#define CHECK_AND_FILL_O2_TRACK_MEAN_ITSCLUSTER_SIZE_COS(FEATURE, v1, v2)                                            \
  case static_cast<uint8_t>(InputFeaturesO2Track::FEATURE): {                                                        \
    int nsize = 0;                                                                                                   \
    int ncls = 0;                                                                                                    \
    for (int il = v1; il < v2; il++) {                                                                               \
      nsize += track.itsClsSizeInLayer(il);                                                                          \
      if (nsize > 0) {                                                                                               \
        ncls++;                                                                                                      \
      }                                                                                                              \
    }                                                                                                                \
    inputFeature = static_cast<float>(nsize) / static_cast<float>(ncls) * std::cos(std::atan(trackParCov.getTgl())); \
    break;                                                                                                           \
  }

// TPC+TOF combined nSigma
#define CHECK_AND_FILL_O2_TRACK_TPCTOF(FEATURE, GETTER1, GETTER2, GETTER3)                    \
  case static_cast<uint8_t>(InputFeaturesO2Track::FEATURE): {                                 \
    if (!track.GETTER3()) {                                                                   \
      inputFeature = track.GETTER1();                                                         \
    } else {                                                                                  \
      if (track.GETTER1() > 0) {                                                              \
        inputFeature = sqrt((pow(track.GETTER1(), 2) + pow(track.GETTER2(), 2)) / 2.);        \
      } else {                                                                                \
        inputFeature = (-1) * sqrt((pow(track.GETTER1(), 2) + pow(track.GETTER2(), 2)) / 2.); \
      }                                                                                       \
    }                                                                                         \
    break;                                                                                    \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER from track and applying a sqrt
#define CHECK_AND_FILL_O2_TRACK_SQRT(FEATURE, GETTER)         \
  case static_cast<uint8_t>(InputFeaturesO2Track::FEATURE): { \
    inputFeature = sqrt(track.GETTER());                      \
    break;                                                    \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER1 from track and multiplying with cos(atan(GETTER2))
#define CHECK_AND_FILL_O2_TRACK_COS(FEATURE, GETTER1, GETTER2)             \
  case static_cast<uint8_t>(InputFeaturesO2Track::FEATURE): {              \
    inputFeature = track.GETTER1() * std::cos(std::atan(track.GETTER2())); \
    break;                                                                 \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER1 and GETTER2 from track.
#define CHECK_AND_FILL_O2_TRACK_RELDIFF(FEATURE, GETTER1, GETTER2)                    \
  case static_cast<uint8_t>(InputFeaturesO2Track::FEATURE): {                         \
    inputFeature = (track.GETTER2() - trackParCov.GETTER1()) / trackParCov.GETTER1(); \
    break;                                                                            \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER=FEATURE from collision
#define CHECK_AND_FILL_DIELECTRON_COLLISION(GETTER)          \
  case static_cast<uint8_t>(InputFeaturesO2Track::GETTER): { \
    inputFeature = collision.GETTER();                       \
    break;                                                   \
  }

namespace o2::analysis
{
// possible input features for ML
enum class InputFeaturesO2Track : uint8_t {
  tpcInnerParam,
  reldiffp,
  tpcSignal,
  tpcNSigmaEl,
  // tpcNSigmaMu,
  tpcNSigmaPi,
  tpcNSigmaKa,
  tpcNSigmaPr,
  beta,
  tofNSigmaEl,
  // tofNSigmaMu,
  tofNSigmaPi,
  tofNSigmaKa,
  tofNSigmaPr,
  tpctofNSigmaEl,
  // tpctofNSigmaMu,
  tpctofNSigmaPi,
  tpctofNSigmaKa,
  tpctofNSigmaPr,
  tpcNClsFound,
  tpcNClsPID,
  tpcNClsCrossedRows,
  tpcChi2NCl,
  hasITS,
  hasTPC,
  hasTRD,
  hasTOF,
  tgl,
  p,
  itsClusterSizes,
  meanClusterSizeITS,
  meanClusterSizeITSib,
  meanClusterSizeITSob,
  meanClusterSizeITSCos,
  meanClusterSizeITSibCos,
  meanClusterSizeITSobCos,
  posZ,
  numContrib,
  trackOccupancyInTimeRange,
  ft0cOccupancyInTimeRange,
};

template <typename TypeOutputScore = float>
class MlResponseO2Track : public MlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  MlResponseO2Track() = default;
  /// Default destructor
  virtual ~MlResponseO2Track() = default;

  template <typename T, typename U, typename V>
  float return_feature(uint8_t idx, T const& track, U const& trackParCov, V const& collision)
  {
    float inputFeature = 0.;
    switch (idx) {
      CHECK_AND_FILL_O2_TRACK(tpcInnerParam);
      CHECK_AND_FILL_O2_TRACK_RELDIFF(reldiffp, getP, tpcInnerParam);
      CHECK_AND_FILL_O2_TRACK(tpcSignal);
      CHECK_AND_FILL_O2_TRACK(tpcNSigmaEl);
      // CHECK_AND_FILL_O2_TRACK(tpcNSigmaMu);
      CHECK_AND_FILL_O2_TRACK(tpcNSigmaPi);
      CHECK_AND_FILL_O2_TRACK(tpcNSigmaKa);
      CHECK_AND_FILL_O2_TRACK(tpcNSigmaPr);
      CHECK_AND_FILL_O2_TRACK(beta);
      CHECK_AND_FILL_O2_TRACK(tofNSigmaEl);
      // CHECK_AND_FILL_O2_TRACK(tofNSigmaMu);
      CHECK_AND_FILL_O2_TRACK(tofNSigmaPi);
      CHECK_AND_FILL_O2_TRACK(tofNSigmaKa);
      CHECK_AND_FILL_O2_TRACK(tofNSigmaPr);
      CHECK_AND_FILL_O2_TRACK_TPCTOF(tpctofNSigmaEl, tpcNSigmaEl, tofNSigmaEl, hasTOF);
      // CHECK_AND_FILL_O2_TRACK_TPCTOF(tpctofNSigmaMu, tpcNSigmaMu, tofNSigmaMu, hasTOF);
      CHECK_AND_FILL_O2_TRACK_TPCTOF(tpctofNSigmaPi, tpcNSigmaPi, tofNSigmaPi, hasTOF);
      CHECK_AND_FILL_O2_TRACK_TPCTOF(tpctofNSigmaKa, tpcNSigmaKa, tofNSigmaKa, hasTOF);
      CHECK_AND_FILL_O2_TRACK_TPCTOF(tpctofNSigmaPr, tpcNSigmaPr, tofNSigmaPr, hasTOF);
      CHECK_AND_FILL_O2_TRACK(tpcNClsFound);
      CHECK_AND_FILL_O2_TRACK(tpcNClsPID);
      CHECK_AND_FILL_O2_TRACK(tpcNClsCrossedRows);
      CHECK_AND_FILL_O2_TRACK(tpcChi2NCl);
      CHECK_AND_FILL_O2_TRACK(hasITS);
      CHECK_AND_FILL_O2_TRACK(hasTPC);
      CHECK_AND_FILL_O2_TRACK(hasTRD);
      CHECK_AND_FILL_O2_TRACK(hasTOF);
      CHECK_AND_FILL_O2_TRACKPARCOV(tgl, getTgl);
      CHECK_AND_FILL_O2_TRACKPARCOV(p, getP);
      CHECK_AND_FILL_O2_TRACK(itsClusterSizes);
      CHECK_AND_FILL_O2_TRACK_MEAN_ITSCLUSTER_SIZE(meanClusterSizeITS, 0, 7);
      CHECK_AND_FILL_O2_TRACK_MEAN_ITSCLUSTER_SIZE(meanClusterSizeITSib, 0, 3);
      CHECK_AND_FILL_O2_TRACK_MEAN_ITSCLUSTER_SIZE(meanClusterSizeITSob, 3, 7);
      CHECK_AND_FILL_O2_TRACK_MEAN_ITSCLUSTER_SIZE_COS(meanClusterSizeITSCos, 0, 7);
      CHECK_AND_FILL_O2_TRACK_MEAN_ITSCLUSTER_SIZE_COS(meanClusterSizeITSibCos, 0, 3);
      CHECK_AND_FILL_O2_TRACK_MEAN_ITSCLUSTER_SIZE_COS(meanClusterSizeITSobCos, 3, 7);
      CHECK_AND_FILL_DIELECTRON_COLLISION(posZ);
      CHECK_AND_FILL_DIELECTRON_COLLISION(numContrib);
      CHECK_AND_FILL_DIELECTRON_COLLISION(trackOccupancyInTimeRange);
      CHECK_AND_FILL_DIELECTRON_COLLISION(ft0cOccupancyInTimeRange);
    }
    return inputFeature;
  }

  /// Method to get the input features vector needed for ML inference
  /// \param track is the single track, \param collision is the collision
  /// \return inputFeatures vector
  template <typename T, typename U, typename V>
  std::vector<float> getInputFeatures(T const& track, U const& trackParCov, V const& collision)
  {
    std::vector<float> inputFeatures;
    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      float inputFeature = return_feature(idx, track, trackParCov, collision);
      inputFeatures.emplace_back(inputFeature);
    }
    return inputFeatures;
  }

  /// Method to get the value of variable chosen for binning
  /// \param track is the single track, \param collision is the collision
  /// \return binning variable
  template <typename T, typename U, typename V>
  float getBinningFeature(T const& track, U const& trackParCov, V const& collision)
  {
    return return_feature(mCachedIndexBinning, track, trackParCov, collision);
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
      FILL_MAP_O2_TRACK(tpcInnerParam),
      FILL_MAP_O2_TRACK(reldiffp),
      FILL_MAP_O2_TRACK(tpcSignal),
      FILL_MAP_O2_TRACK(tpcNSigmaEl),
      // FILL_MAP_O2_TRACK(tpcNSigmaMu),
      FILL_MAP_O2_TRACK(tpcNSigmaPi),
      FILL_MAP_O2_TRACK(tpcNSigmaKa),
      FILL_MAP_O2_TRACK(tpcNSigmaPr),
      FILL_MAP_O2_TRACK(beta),
      FILL_MAP_O2_TRACK(tofNSigmaEl),
      // FILL_MAP_O2_TRACK(tofNSigmaMu),
      FILL_MAP_O2_TRACK(tofNSigmaPi),
      FILL_MAP_O2_TRACK(tofNSigmaKa),
      FILL_MAP_O2_TRACK(tofNSigmaPr),
      FILL_MAP_O2_TRACK(tpctofNSigmaEl),
      // FILL_MAP_O2_TRACK(tpctofNSigmaMu),
      FILL_MAP_O2_TRACK(tpctofNSigmaPi),
      FILL_MAP_O2_TRACK(tpctofNSigmaKa),
      FILL_MAP_O2_TRACK(tpctofNSigmaPr),
      FILL_MAP_O2_TRACK(tpcNClsFound),
      FILL_MAP_O2_TRACK(tpcNClsPID),
      FILL_MAP_O2_TRACK(tpcNClsCrossedRows),
      FILL_MAP_O2_TRACK(tpcChi2NCl),
      FILL_MAP_O2_TRACK(hasITS),
      FILL_MAP_O2_TRACK(hasTPC),
      FILL_MAP_O2_TRACK(hasTRD),
      FILL_MAP_O2_TRACK(hasTOF),
      FILL_MAP_O2_TRACK(tgl),
      FILL_MAP_O2_TRACK(p),
      FILL_MAP_O2_TRACK(itsClusterSizes),
      FILL_MAP_O2_TRACK(meanClusterSizeITS),
      FILL_MAP_O2_TRACK(meanClusterSizeITSib),
      FILL_MAP_O2_TRACK(meanClusterSizeITSob),
      FILL_MAP_O2_TRACK(meanClusterSizeITSCos),
      FILL_MAP_O2_TRACK(meanClusterSizeITSibCos),
      FILL_MAP_O2_TRACK(meanClusterSizeITSobCos),
      FILL_MAP_O2_TRACK(posZ),
      FILL_MAP_O2_TRACK(numContrib),
      FILL_MAP_O2_TRACK(trackOccupancyInTimeRange),
      FILL_MAP_O2_TRACK(ft0cOccupancyInTimeRange)};
  }

  uint8_t mCachedIndexBinning; // index correspondance between configurable and available input features
};

} // namespace o2::analysis

#undef FILL_MAP_O2_TRACK
#undef CHECK_AND_FILL_O2_TRACK
#undef CHECK_AND_FILL_O2_TRACKPARCOV
#undef CHECK_AND_FILL_O2_TRACK_MEAN_ITSCLUSTER_SIZE
#undef CHECK_AND_FILL_O2_TRACK_MEAN_ITSCLUSTER_SIZE_COS
#undef CHECK_AND_FILL_O2_TRACK_SQRT
#undef CHECK_AND_FILL_O2_TRACK_COS
#undef CHECK_AND_FILL_O2_TRACK_TPCTOF
#undef CHECK_AND_FILL_O2_TRACK_RELDIFF
#undef CHECK_AND_FILL_DIELECTRON_COLLISION

#endif // PWGEM_DILEPTON_UTILS_MLRESPONSEO2TRACK_H_
