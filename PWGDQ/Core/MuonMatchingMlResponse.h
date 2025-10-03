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

/// \file MuonMatchingMlResponse.h
/// \brief Class to compute the ML response for MFT-Muon matching
/// \author Maurice Coquet <maurice.louis.coquet@cern.ch>

#ifndef PWGDQ_CORE_MUONMATCHINGMLRESPONSE_H_
#define PWGDQ_CORE_MUONMATCHINGMLRESPONSE_H_

#include "Tools/ML/MlResponse.h"

#include <map>
#include <string>
#include <vector>

// Fill the map of available input features
// the key is the feature's name (std::string)
// the value is the corresponding value in EnumInputFeatures
#define FILL_MAP_MFTMUON_MATCH(FEATURE)                                \
  {                                                                    \
    #FEATURE, static_cast<uint8_t>(InputFeaturesMFTMuonMatch::FEATURE) \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER=FEATURE from track
#define CHECK_AND_FILL_MUON_TRACK(FEATURE, GETTER)                 \
  case static_cast<uint8_t>(InputFeaturesMFTMuonMatch::FEATURE): { \
    inputFeature = muon.GETTER();                                  \
    break;                                                         \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER=FEATURE from track
#define CHECK_AND_FILL_MUONGLOB_TRACK(FEATURE, GETTER)             \
  case static_cast<uint8_t>(InputFeaturesMFTMuonMatch::FEATURE): { \
    inputFeature = muonglob.GETTER();                              \
    break;                                                         \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER=FEATURE from track
#define CHECK_AND_FILL_MFT_TRACK(FEATURE, GETTER)                  \
  case static_cast<uint8_t>(InputFeaturesMFTMuonMatch::FEATURE): { \
    inputFeature = mft.GETTER();                                   \
    break;                                                         \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER=FEATURE from track
#define CHECK_AND_FILL_MUON_COV(FEATURE, GETTER)                   \
  case static_cast<uint8_t>(InputFeaturesMFTMuonMatch::FEATURE): { \
    inputFeature = muoncov.GETTER();                               \
    break;                                                         \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER=FEATURE from track
#define CHECK_AND_FILL_MFT_COV(FEATURE, GETTER)                    \
  case static_cast<uint8_t>(InputFeaturesMFTMuonMatch::FEATURE): { \
    inputFeature = mftcov.GETTER();                                \
    break;                                                         \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER1 and GETTER2 from track.
#define CHECK_AND_FILL_MFTMUON_DIFF(FEATURE, GETTER1, GETTER2)     \
  case static_cast<uint8_t>(InputFeaturesMFTMuonMatch::FEATURE): { \
    inputFeature = (mft.GETTER2() - muon.GETTER1());               \
    break;                                                         \
  }

// Check if the index of mCachedIndices (index associated to a FEATURE)
// matches the entry in EnumInputFeatures associated to this FEATURE
// if so, the inputFeatures vector is filled with the FEATURE's value
// by calling the corresponding GETTER=FEATURE from collision
#define CHECK_AND_FILL_MFTMUON_COLLISION(GETTER)                  \
  case static_cast<uint8_t>(InputFeaturesMFTMuonMatch::GETTER): { \
    inputFeature = collision.GETTER();                            \
    break;                                                        \
  }

namespace o2::analysis
{
// possible input features for ML
enum class InputFeaturesMFTMuonMatch : uint8_t {
  zMatching,
  xMFT,
  yMFT,
  qOverptMFT,
  tglMFT,
  phiMFT,
  dcaXY,
  dcaZ,
  chi2MFT,
  nClustersMFT,
  xMCH,
  yMCH,
  qOverptMCH,
  tglMCH,
  phiMCH,
  nClustersMCH,
  chi2MCH,
  pdca,
  Rabs,
  cXXMFT,
  cXYMFT,
  cYYMFT,
  cPhiYMFT,
  cPhiXMFT,
  cPhiPhiMFT,
  cTglYMFT,
  cTglXMFT,
  cTglPhiMFT,
  cTglTglMFT,
  c1PtYMFT,
  c1PtXMFT,
  c1PtPhiMFT,
  c1PtTglMFT,
  c1Pt21Pt2MFT,
  cXXMCH,
  cXYMCH,
  cYYMCH,
  cPhiYMCH,
  cPhiXMCH,
  cPhiPhiMCH,
  cTglYMCH,
  cTglXMCH,
  cTglPhiMCH,
  cTglTglMCH,
  c1PtYMCH,
  c1PtXMCH,
  c1PtPhiMCH,
  c1PtTglMCH,
  c1Pt21Pt2MCH,
  deltaX,
  deltaY,
  deltaPhi,
  deltaEta,
  deltaPt,
  posX,
  posY,
  posZ,
  numContrib,
  trackOccupancyInTimeRange,
  ft0cOccupancyInTimeRange,
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
  centFT0C,
  chi2MCHMFT,
  chi2GlobMUON
};

template <typename TypeOutputScore = float>
class MlResponseMFTMuonMatch : public MlResponse<TypeOutputScore>
{
 public:
  /// Default constructor
  MlResponseMFTMuonMatch() = default;
  /// Default destructor
  virtual ~MlResponseMFTMuonMatch() = default;

  template <typename T1, typename T2, typename C1, typename C2, typename U>
  float returnFeature(uint8_t idx, T1 const& muon, T2 const& mft, C1 const& muoncov, C2 const& mftcov, U const& collision)
  {
    float inputFeature = 0.;
    switch (idx) {
      CHECK_AND_FILL_MFT_TRACK(zMatching, z);
      CHECK_AND_FILL_MFT_TRACK(xMFT, x);
      CHECK_AND_FILL_MFT_TRACK(yMFT, y);
      CHECK_AND_FILL_MFT_TRACK(qOverptMFT, signed1Pt);
      CHECK_AND_FILL_MFT_TRACK(tglMFT, tgl);
      CHECK_AND_FILL_MFT_TRACK(phiMFT, phi);
      CHECK_AND_FILL_MFT_TRACK(chi2MFT, chi2);
      CHECK_AND_FILL_MFT_TRACK(nClustersMFT, nClusters);
      CHECK_AND_FILL_MUON_TRACK(dcaXY, fwddcaXY);
      CHECK_AND_FILL_MUON_TRACK(dcaZ, fwddcaz);
      CHECK_AND_FILL_MUON_TRACK(xMCH, x);
      CHECK_AND_FILL_MUON_TRACK(yMCH, y);
      CHECK_AND_FILL_MUON_TRACK(qOverptMCH, signed1Pt);
      CHECK_AND_FILL_MUON_TRACK(tglMCH, tgl);
      CHECK_AND_FILL_MUON_TRACK(phiMCH, phi);
      CHECK_AND_FILL_MUON_TRACK(nClustersMCH, nClusters);
      CHECK_AND_FILL_MUON_TRACK(chi2MCH, chi2);
      CHECK_AND_FILL_MUON_TRACK(pdca, pDca);
      CHECK_AND_FILL_MFT_COV(cXXMFT, cXX);
      CHECK_AND_FILL_MFT_COV(cXYMFT, cXY);
      CHECK_AND_FILL_MFT_COV(cYYMFT, cYY);
      CHECK_AND_FILL_MFT_COV(cPhiYMFT, cPhiY);
      CHECK_AND_FILL_MFT_COV(cPhiXMFT, cPhiX);
      CHECK_AND_FILL_MFT_COV(cPhiPhiMFT, cPhiPhi);
      CHECK_AND_FILL_MFT_COV(cTglYMFT, cTglY);
      CHECK_AND_FILL_MFT_COV(cTglXMFT, cTglX);
      CHECK_AND_FILL_MFT_COV(cTglPhiMFT, cTglPhi);
      CHECK_AND_FILL_MFT_COV(cTglTglMFT, cTglTgl);
      CHECK_AND_FILL_MFT_COV(c1PtYMFT, c1PtY);
      CHECK_AND_FILL_MFT_COV(c1PtXMFT, c1PtX);
      CHECK_AND_FILL_MFT_COV(c1PtPhiMFT, c1PtPhi);
      CHECK_AND_FILL_MFT_COV(c1PtTglMFT, c1PtTgl);
      CHECK_AND_FILL_MFT_COV(c1Pt21Pt2MFT, c1Pt21Pt2);
      CHECK_AND_FILL_MUON_COV(cXXMCH, cXX);
      CHECK_AND_FILL_MUON_COV(cXYMCH, cXY);
      CHECK_AND_FILL_MUON_COV(cYYMCH, cYY);
      CHECK_AND_FILL_MUON_COV(cPhiYMCH, cPhiY);
      CHECK_AND_FILL_MUON_COV(cPhiXMCH, cPhiX);
      CHECK_AND_FILL_MUON_COV(cPhiPhiMCH, cPhiPhi);
      CHECK_AND_FILL_MUON_COV(cTglYMCH, cTglY);
      CHECK_AND_FILL_MUON_COV(cTglXMCH, cTglX);
      CHECK_AND_FILL_MUON_COV(cTglPhiMCH, cTglPhi);
      CHECK_AND_FILL_MUON_COV(cTglTglMCH, cTglTgl);
      CHECK_AND_FILL_MUON_COV(c1PtYMCH, c1PtY);
      CHECK_AND_FILL_MUON_COV(c1PtXMCH, c1PtX);
      CHECK_AND_FILL_MUON_COV(c1PtPhiMCH, c1PtPhi);
      CHECK_AND_FILL_MUON_COV(c1PtTglMCH, c1PtTgl);
      CHECK_AND_FILL_MUON_COV(c1Pt21Pt2MCH, c1Pt21Pt2);
      CHECK_AND_FILL_MFTMUON_COLLISION(posX);
      CHECK_AND_FILL_MFTMUON_COLLISION(posY);
      CHECK_AND_FILL_MFTMUON_COLLISION(posZ);
      CHECK_AND_FILL_MFTMUON_COLLISION(numContrib);
      CHECK_AND_FILL_MFTMUON_COLLISION(trackOccupancyInTimeRange);
      CHECK_AND_FILL_MFTMUON_COLLISION(ft0cOccupancyInTimeRange);
      CHECK_AND_FILL_MFTMUON_COLLISION(multFT0A);
      CHECK_AND_FILL_MFTMUON_COLLISION(multFT0C);
      CHECK_AND_FILL_MFTMUON_COLLISION(multNTracksPV);
      CHECK_AND_FILL_MFTMUON_COLLISION(multNTracksPVeta1);
      CHECK_AND_FILL_MFTMUON_COLLISION(multNTracksPVetaHalf);
      CHECK_AND_FILL_MFTMUON_COLLISION(isInelGt0);
      CHECK_AND_FILL_MFTMUON_COLLISION(isInelGt1);
      CHECK_AND_FILL_MFTMUON_COLLISION(multFT0M);
      CHECK_AND_FILL_MFTMUON_COLLISION(centFT0M);
      CHECK_AND_FILL_MFTMUON_COLLISION(centFT0A);
      CHECK_AND_FILL_MFTMUON_COLLISION(centFT0C);
      CHECK_AND_FILL_MUON_TRACK(chi2MCHMFT, chi2MatchMCHMFT);
    }
    return inputFeature;
  }

  template <typename T1, typename T2, typename T3, typename U>
  float returnFeatureGlob(uint8_t idx, T1 const& muonglob, T2 const& muon, T3 const& mft, U const& collision)
  {
    float inputFeature = 0.;
    switch (idx) {
      CHECK_AND_FILL_MFT_TRACK(xMFT, getX);
      CHECK_AND_FILL_MFT_TRACK(yMFT, getY);
      CHECK_AND_FILL_MFT_TRACK(qOverptMFT, getInvQPt);
      CHECK_AND_FILL_MFT_TRACK(tglMFT, getTanl);
      CHECK_AND_FILL_MFT_TRACK(phiMFT, getPhi);
      CHECK_AND_FILL_MFT_TRACK(chi2MFT, getTrackChi2);
      CHECK_AND_FILL_MUON_TRACK(xMCH, getX);
      CHECK_AND_FILL_MUON_TRACK(yMCH, getY);
      CHECK_AND_FILL_MUON_TRACK(qOverptMCH, getInvQPt);
      CHECK_AND_FILL_MUON_TRACK(tglMCH, getTanl);
      CHECK_AND_FILL_MUON_TRACK(phiMCH, getPhi);
      CHECK_AND_FILL_MUON_TRACK(chi2MCH, getTrackChi2);
      CHECK_AND_FILL_MUONGLOB_TRACK(chi2MCHMFT, chi2MatchMCHMFT);
      CHECK_AND_FILL_MUONGLOB_TRACK(chi2GlobMUON, chi2);
      CHECK_AND_FILL_MUONGLOB_TRACK(Rabs, rAtAbsorberEnd);
      // Below are dummy files to remove warning of unused parameters
      CHECK_AND_FILL_MFTMUON_COLLISION(posZ);
    }
    return inputFeature;
  }

  template <typename T1>
  float returnFeatureTest(uint8_t idx, T1 const& muon)
  {
    float inputFeature = 0.;
    switch (idx) {
      CHECK_AND_FILL_MUON_TRACK(chi2MCHMFT, chi2MatchMCHMFT);
    }
    return inputFeature;
  }

  /// Method to get the input features vector needed for ML inference
  /// \param track is the single track, \param collision is the collision
  /// \return inputFeatures vector
  template <typename T1, typename T2, typename C1, typename C2, typename U>
  std::vector<float> getInputFeatures(T1 const& muon, T2 const& mft, C1 const& muoncov, C2 const& mftcov, U const& collision)
  {
    std::vector<float> inputFeatures;
    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      float inputFeature = returnFeature(idx, muon, mft, muoncov, mftcov, collision);
      inputFeatures.emplace_back(inputFeature);
    }
    return inputFeatures;
  }

  template <typename T1>
  std::vector<float> getInputFeaturesTest(T1 const& muon)
  {
    std::vector<float> inputFeatures;
    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      float inputFeature = returnFeatureTest(idx, muon);
      inputFeatures.emplace_back(inputFeature);
    }
    return inputFeatures;
  }

  template <typename T1, typename T2, typename T3, typename U>
  std::vector<float> getInputFeaturesGlob(T1 const& muonglob, T2 const& muon, T3 const& mft, U const& collision)
  {
    std::vector<float> inputFeatures;
    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      float inputFeature = returnFeatureGlob(idx, muonglob, muon, mft, collision);
      inputFeatures.emplace_back(inputFeature);
    }
    return inputFeatures;
  }

  /// Method to get the value of variable chosen for binning
  /// \param track is the single track, \param collision is the collision
  /// \return binning variable
  template <typename T1, typename T2, typename C1, typename C2, typename U>
  float getBinningFeature(T1 const& muon, T2 const& mft, C1 const& muoncov, C2 const& mftcov, U const& collision)
  {
    return returnFeature(mCachedIndexBinning, muon, mft, muoncov, mftcov, collision);
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
      FILL_MAP_MFTMUON_MATCH(zMatching),
      FILL_MAP_MFTMUON_MATCH(xMFT),
      FILL_MAP_MFTMUON_MATCH(yMFT),
      FILL_MAP_MFTMUON_MATCH(qOverptMFT),
      FILL_MAP_MFTMUON_MATCH(tglMFT),
      FILL_MAP_MFTMUON_MATCH(phiMFT),
      FILL_MAP_MFTMUON_MATCH(dcaXY),
      FILL_MAP_MFTMUON_MATCH(dcaZ),
      FILL_MAP_MFTMUON_MATCH(chi2MFT),
      FILL_MAP_MFTMUON_MATCH(nClustersMFT),
      FILL_MAP_MFTMUON_MATCH(xMCH),
      FILL_MAP_MFTMUON_MATCH(yMCH),
      FILL_MAP_MFTMUON_MATCH(qOverptMCH),
      FILL_MAP_MFTMUON_MATCH(tglMCH),
      FILL_MAP_MFTMUON_MATCH(phiMCH),
      FILL_MAP_MFTMUON_MATCH(nClustersMCH),
      FILL_MAP_MFTMUON_MATCH(chi2MCH),
      FILL_MAP_MFTMUON_MATCH(pdca),
      FILL_MAP_MFTMUON_MATCH(Rabs),
      FILL_MAP_MFTMUON_MATCH(cXXMFT),
      FILL_MAP_MFTMUON_MATCH(cXYMFT),
      FILL_MAP_MFTMUON_MATCH(cYYMFT),
      FILL_MAP_MFTMUON_MATCH(cPhiYMFT),
      FILL_MAP_MFTMUON_MATCH(cPhiXMFT),
      FILL_MAP_MFTMUON_MATCH(cPhiPhiMFT),
      FILL_MAP_MFTMUON_MATCH(cTglYMFT),
      FILL_MAP_MFTMUON_MATCH(cTglXMFT),
      FILL_MAP_MFTMUON_MATCH(cTglPhiMFT),
      FILL_MAP_MFTMUON_MATCH(cTglTglMFT),
      FILL_MAP_MFTMUON_MATCH(c1PtYMFT),
      FILL_MAP_MFTMUON_MATCH(c1PtXMFT),
      FILL_MAP_MFTMUON_MATCH(c1PtPhiMFT),
      FILL_MAP_MFTMUON_MATCH(c1PtTglMFT),
      FILL_MAP_MFTMUON_MATCH(c1Pt21Pt2MFT),
      FILL_MAP_MFTMUON_MATCH(cXXMCH),
      FILL_MAP_MFTMUON_MATCH(cXYMCH),
      FILL_MAP_MFTMUON_MATCH(cYYMCH),
      FILL_MAP_MFTMUON_MATCH(cPhiYMCH),
      FILL_MAP_MFTMUON_MATCH(cPhiXMCH),
      FILL_MAP_MFTMUON_MATCH(cPhiPhiMCH),
      FILL_MAP_MFTMUON_MATCH(cTglYMCH),
      FILL_MAP_MFTMUON_MATCH(cTglXMCH),
      FILL_MAP_MFTMUON_MATCH(cTglPhiMCH),
      FILL_MAP_MFTMUON_MATCH(cTglTglMCH),
      FILL_MAP_MFTMUON_MATCH(c1PtYMCH),
      FILL_MAP_MFTMUON_MATCH(c1PtXMCH),
      FILL_MAP_MFTMUON_MATCH(c1PtPhiMCH),
      FILL_MAP_MFTMUON_MATCH(c1PtTglMCH),
      FILL_MAP_MFTMUON_MATCH(c1Pt21Pt2MCH),
      FILL_MAP_MFTMUON_MATCH(chi2MCHMFT),
      FILL_MAP_MFTMUON_MATCH(chi2GlobMUON)};
  }

  uint8_t mCachedIndexBinning; // index correspondance between configurable and available input features
};

} // namespace o2::analysis

#undef FILL_MAP_MFTMUON_MAP
#undef CHECK_AND_FILL_MUON_TRACK
#undef CHECK_AND_FILL_MFT_TRACK
#undef CHECK_AND_FILL_MUON_COV
#undef CHECK_AND_FILL_MFT_COV
#undef CHECK_AND_FILL_MFTMUON_DIFF
#undef CHECK_AND_FILL_MFTMUON_COLLISION

#endif // PWGDQ_CORE_MUONMATCHINGMLRESPONSE_H_
