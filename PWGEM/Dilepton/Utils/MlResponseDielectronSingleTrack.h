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

#include "Tools/ML/MlResponse.h"

#include <map>
#include <string>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_DIELECTRON_SINGLE_TRACK(FEATURE)                               \
  {                                                                             \
    #FEATURE, static_cast<uint8_t>(InputFeaturesDielectronSingleTrack::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER=FEATURE from track
#define CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(GETTER)                     \
  case static_cast<uint8_t>(InputFeaturesDielectronSingleTrack::GETTER): { \
    inputFeature = track.GETTER();                                         \
    break;                                                                 \
  }

// TPC+TOF combined nSigma
#define CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK_TPCTOF(FEATURE, GETTER1, GETTER2, GETTER3)     \
  case static_cast<uint8_t>(InputFeaturesDielectronSingleTrack::FEATURE): {                   \
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
#define CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK_SQRT(FEATURE, GETTER)        \
  case static_cast<uint8_t>(InputFeaturesDielectronSingleTrack::FEATURE): { \
    inputFeature = sqrt(track.GETTER());                                    \
    break;                                                                  \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER1 from track and multiplying with cos(atan(GETTER2))
#define CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK_COS(FEATURE, GETTER1, GETTER2) \
  case static_cast<uint8_t>(InputFeaturesDielectronSingleTrack::FEATURE): {   \
    inputFeature = track.GETTER1() * std::cos(std::atan(track.GETTER2()));    \
    break;                                                                    \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER1 and GETTER2 from track.
#define CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK_RELDIFF(FEATURE, GETTER1, GETTER2) \
  case static_cast<uint8_t>(InputFeaturesDielectronSingleTrack::FEATURE): {       \
    inputFeature = (track.GETTER2() - track.GETTER1()) / track.GETTER1();         \
    break;                                                                        \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER=FEATURE from collision
#define CHECK_AND_FILL_DIELECTRON_COLLISION(GETTER)                        \
  case static_cast<uint8_t>(InputFeaturesDielectronSingleTrack::GETTER): { \
    inputFeature = collision.GETTER();                                     \
    break;                                                                 \
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
  tpcNClsFindableMinusFound,
  tpcNClsFindableMinusCrossedRows,
  tpcNClsShared,
  tpcChi2NCl,
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
  itsClusterSizes,
  itsChi2NCl,
  tofChi2,
  detectorMap,
  x,
  alpha,
  y,
  z,
  snp,
  tgl,
  isAssociatedToMPC,
  tpcNClsFound,
  tpcNClsPID,
  tpcNClsCrossedRows,
  tpcCrossedRowsOverFindableCls,
  tpcFoundOverFindableCls,
  tpcFractionSharedCls,
  itsClusterMap,
  itsNCls,
  itsNClsInnerBarrel,
  hasITS,
  hasTPC,
  hasTRD,
  hasTOF,
  signed1Pt,
  p,
  px,
  py,
  pz,
  theta,
  meanClusterSizeITS,
  meanClusterSizeITSib,
  meanClusterSizeITSob,
  meanClusterSizeITSCos,
  meanClusterSizeITSibCos,
  meanClusterSizeITSobCos,
  cYY,
  cZY,
  cZZ,
  cSnpY,
  cSnpZ,
  cSnpSnp,
  cTglY,
  cTglZ,
  cTglSnp,
  cTglTgl,
  c1PtY,
  c1PtZ,
  c1PtSnp,
  c1PtTgl,
  c1Pt21Pt2,
  posX,
  posY,
  posZ,
  numContrib,
  trackOccupancyInTimeRange,
  ft0cOccupancyInTimeRange,
  // covXX,
  // covXY,
  // covXZ,
  // covYY,
  // covYZ,
  // covZZ,
  // chi2,
  multFT0A,
  multFT0C,
  multNTracksPV,
  multNTracksPVeta1,
  multNTracksPVetaHalf,
  isInelGt0,
  isInelGt1,
  multFT0M,
  centFT0M,
  centFT0A,
  centFT0C
};

template <typename TypeOutputScore = float>
class MlResponseDielectronSingleTrack : public MlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  MlResponseDielectronSingleTrack() = default;
  /// Default destructor
  virtual ~MlResponseDielectronSingleTrack() = default;

  template <typename T, typename U>
  float return_feature(uint8_t idx, T const& track, U const& collision)
  {
    float inputFeature = 0.;
    switch (idx) {
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(sign);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(pt);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(eta);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(phi);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(dcaXY);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(dcaZ);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK_SQRT(dcaResXY, cYY);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK_SQRT(dcaResZ, cZZ);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tpcNClsFindable);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tpcNClsFindableMinusFound);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tpcNClsFindableMinusCrossedRows);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tpcNClsShared);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tpcChi2NCl);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tpcInnerParam);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK_RELDIFF(reldiffp, p, tpcInnerParam);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tpcSignal);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tpcNSigmaEl);
      // CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tpcNSigmaMu);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tpcNSigmaPi);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tpcNSigmaKa);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tpcNSigmaPr);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(beta);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tofNSigmaEl);
      // CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tofNSigmaMu);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tofNSigmaPi);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tofNSigmaKa);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tofNSigmaPr);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK_TPCTOF(tpctofNSigmaEl, tpcNSigmaEl, tofNSigmaEl, hasTOF);
      // CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK_TPCTOF(tpctofNSigmaMu, tpcNSigmaMu, tofNSigmaMu, hasTOF);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK_TPCTOF(tpctofNSigmaPi, tpcNSigmaPi, tofNSigmaPi, hasTOF);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK_TPCTOF(tpctofNSigmaKa, tpcNSigmaKa, tofNSigmaKa, hasTOF);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK_TPCTOF(tpctofNSigmaPr, tpcNSigmaPr, tofNSigmaPr, hasTOF);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(itsClusterSizes);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(itsChi2NCl);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tofChi2);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(detectorMap);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(x);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(alpha);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(y);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(z);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(snp);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tgl);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(isAssociatedToMPC);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tpcNClsFound);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tpcNClsPID);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tpcNClsCrossedRows);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tpcCrossedRowsOverFindableCls);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tpcFoundOverFindableCls);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(tpcFractionSharedCls);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(itsClusterMap);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(itsNCls);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(itsNClsInnerBarrel);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(hasITS);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(hasTPC);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(hasTRD);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(hasTOF);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(signed1Pt);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(p);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(px);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(py);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(pz);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(theta);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(meanClusterSizeITS);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(meanClusterSizeITSib);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(meanClusterSizeITSob);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK_COS(meanClusterSizeITSCos, meanClusterSizeITS, tgl);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK_COS(meanClusterSizeITSibCos, meanClusterSizeITSib, tgl);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK_COS(meanClusterSizeITSobCos, meanClusterSizeITSob, tgl);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(cYY);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(cZY);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(cZZ);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(cSnpY);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(cSnpZ);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(cSnpSnp);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(cTglY);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(cTglZ);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(cTglSnp);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(cTglTgl);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(c1PtY);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(c1PtZ);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(c1PtSnp);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(c1PtTgl);
      CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK(c1Pt21Pt2);
      CHECK_AND_FILL_DIELECTRON_COLLISION(posX);
      CHECK_AND_FILL_DIELECTRON_COLLISION(posY);
      CHECK_AND_FILL_DIELECTRON_COLLISION(posZ);
      CHECK_AND_FILL_DIELECTRON_COLLISION(numContrib);
      CHECK_AND_FILL_DIELECTRON_COLLISION(trackOccupancyInTimeRange);
      CHECK_AND_FILL_DIELECTRON_COLLISION(ft0cOccupancyInTimeRange);
      // CHECK_AND_FILL_DIELECTRON_COLLISION(covXX);
      // CHECK_AND_FILL_DIELECTRON_COLLISION(covXY);
      // CHECK_AND_FILL_DIELECTRON_COLLISION(covXZ);
      // CHECK_AND_FILL_DIELECTRON_COLLISION(covYY);
      // CHECK_AND_FILL_DIELECTRON_COLLISION(covYZ);
      // CHECK_AND_FILL_DIELECTRON_COLLISION(covZZ);
      // CHECK_AND_FILL_DIELECTRON_COLLISION(chi2);
      CHECK_AND_FILL_DIELECTRON_COLLISION(multFT0A);
      CHECK_AND_FILL_DIELECTRON_COLLISION(multFT0C);
      CHECK_AND_FILL_DIELECTRON_COLLISION(multNTracksPV);
      CHECK_AND_FILL_DIELECTRON_COLLISION(multNTracksPVeta1);
      CHECK_AND_FILL_DIELECTRON_COLLISION(multNTracksPVetaHalf);
      CHECK_AND_FILL_DIELECTRON_COLLISION(isInelGt0);
      CHECK_AND_FILL_DIELECTRON_COLLISION(isInelGt1);
      CHECK_AND_FILL_DIELECTRON_COLLISION(multFT0M);
      CHECK_AND_FILL_DIELECTRON_COLLISION(centFT0M);
      CHECK_AND_FILL_DIELECTRON_COLLISION(centFT0A);
      CHECK_AND_FILL_DIELECTRON_COLLISION(centFT0C);
    }
    return inputFeature;
  }

  /// Method to get the input features vector needed for ML inference
  /// \param track is the single track, \param collision is the collision
  /// \return inputFeatures vector
  template <typename T, typename U>
  std::vector<float> getInputFeatures(T const& track, U const& collision)
  {
    std::vector<float> inputFeatures;
    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      float inputFeature = return_feature(idx, track, collision);
      inputFeatures.emplace_back(inputFeature);
    }
    return inputFeatures;
  }

  /// Method to get the value of variable chosen for binning
  /// \param track is the single track, \param collision is the collision
  /// \return binning variable
  template <typename T, typename U>
  float getBinningFeature(T const& track, U const& collision)
  {
    return return_feature(mCachedIndexBinning, track, collision);
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
      FILL_MAP_DIELECTRON_SINGLE_TRACK(sign),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(pt),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(eta),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(phi),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(dcaXY),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(dcaZ),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(dcaResXY),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(dcaResZ),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcNClsFindable),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcNClsFindableMinusFound),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcNClsFindableMinusCrossedRows),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcNClsShared),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcChi2NCl),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcInnerParam),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(reldiffp),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcSignal),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcNSigmaEl),
      // FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcNSigmaMu),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcNSigmaPi),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcNSigmaKa),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcNSigmaPr),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(beta),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tofNSigmaEl),
      // FILL_MAP_DIELECTRON_SINGLE_TRACK(tofNSigmaMu),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tofNSigmaPi),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tofNSigmaKa),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tofNSigmaPr),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpctofNSigmaEl),
      // FILL_MAP_DIELECTRON_SINGLE_TRACK(tpctofNSigmaMu),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpctofNSigmaPi),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpctofNSigmaKa),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpctofNSigmaPr),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(itsClusterSizes),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(itsChi2NCl),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tofChi2),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(detectorMap),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(x),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(alpha),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(y),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(z),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(snp),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tgl),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(isAssociatedToMPC),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcNClsFound),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcNClsPID),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcNClsCrossedRows),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcCrossedRowsOverFindableCls),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcFoundOverFindableCls),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(tpcFractionSharedCls),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(itsClusterMap),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(itsNCls),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(itsNClsInnerBarrel),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(hasITS),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(hasTPC),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(hasTRD),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(hasTOF),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(signed1Pt),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(p),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(px),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(py),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(pz),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(theta),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(meanClusterSizeITS),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(meanClusterSizeITSib),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(meanClusterSizeITSob),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(meanClusterSizeITSCos),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(meanClusterSizeITSibCos),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(meanClusterSizeITSobCos),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(cYY),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(cZY),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(cZZ),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(cSnpY),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(cSnpZ),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(cSnpSnp),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(cTglY),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(cTglZ),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(cTglSnp),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(cTglTgl),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(c1PtY),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(c1PtZ),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(c1PtSnp),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(c1PtTgl),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(c1Pt21Pt2),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(posX),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(posY),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(posZ),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(numContrib),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(trackOccupancyInTimeRange),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(ft0cOccupancyInTimeRange),
      // FILL_MAP_DIELECTRON_SINGLE_TRACK(covXX),
      // FILL_MAP_DIELECTRON_SINGLE_TRACK(covXY),
      // FILL_MAP_DIELECTRON_SINGLE_TRACK(covXZ),
      // FILL_MAP_DIELECTRON_SINGLE_TRACK(covYY),
      // FILL_MAP_DIELECTRON_SINGLE_TRACK(covYZ),
      // FILL_MAP_DIELECTRON_SINGLE_TRACK(covZZ),
      // FILL_MAP_DIELECTRON_SINGLE_TRACK(chi2),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(multFT0A),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(multFT0C),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(multNTracksPV),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(multNTracksPVeta1),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(multNTracksPVetaHalf),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(isInelGt0),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(isInelGt1),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(multFT0M),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(centFT0M),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(centFT0A),
      FILL_MAP_DIELECTRON_SINGLE_TRACK(centFT0C)};
  }

  uint8_t mCachedIndexBinning; // index correspondance between configurable and available input features
};

} // namespace o2::analysis

#undef FILL_MAP_DIELECTRON_SINGLE_TRACK
#undef CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK
#undef CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK_SQRT
#undef CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK_COS
#undef CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK_TPCTOF
#undef CHECK_AND_FILL_DIELECTRON_SINGLE_TRACK_RELDIFF
#undef CHECK_AND_FILL_DIELECTRON_COLLISION

#endif // PWGEM_DILEPTON_UTILS_MLRESPONSEDIELECTRONSINGLETRACK_H_
