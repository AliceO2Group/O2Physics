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

#include <Framework/Logger.h>

#include <cstdint>
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
#define CHECK_AND_FILL_FEATURE(FEATURE, GETTER)                    \
  case static_cast<uint8_t>(InputFeaturesMFTMuonMatch::FEATURE): { \
    inputFeature = (GETTER);                                       \
    break;                                                         \
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

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename U>
  float returnFeature(uint8_t idx, T1 const& muon, T2 const& mft, T3 const& mch, T4 const& mftprop, T5 const& mchprop, U const& collision)
  {
    float inputFeature = 0.;
    switch (idx) {
      CHECK_AND_FILL_FEATURE(zMatching, mftprop.getZ());
      CHECK_AND_FILL_FEATURE(xMFT, mftprop.getX());
      CHECK_AND_FILL_FEATURE(yMFT, mftprop.getY());
      CHECK_AND_FILL_FEATURE(qOverptMFT, mftprop.getInvQPt());
      CHECK_AND_FILL_FEATURE(tglMFT, mftprop.getTanl());
      CHECK_AND_FILL_FEATURE(phiMFT, mftprop.getPhi());
      CHECK_AND_FILL_FEATURE(chi2MFT, mft.chi2());
      CHECK_AND_FILL_FEATURE(nClustersMFT, mft.nClusters());
      /*dummy value*/ CHECK_AND_FILL_FEATURE(dcaXY, 0);
      /*dummy value*/ CHECK_AND_FILL_FEATURE(dcaZ, 0);
      CHECK_AND_FILL_FEATURE(xMCH, mchprop.getX());
      CHECK_AND_FILL_FEATURE(yMCH, mchprop.getY());
      CHECK_AND_FILL_FEATURE(qOverptMCH, mchprop.getInvQPt());
      CHECK_AND_FILL_FEATURE(tglMCH, mchprop.getTanl());
      CHECK_AND_FILL_FEATURE(phiMCH, mchprop.getPhi());
      CHECK_AND_FILL_FEATURE(nClustersMCH, mch.nClusters());
      CHECK_AND_FILL_FEATURE(chi2MCH, mch.chi2());
      CHECK_AND_FILL_FEATURE(pdca, muon.pDca());
      CHECK_AND_FILL_FEATURE(cXXMFT, mftprop.getCovariances()(0, 0));
      CHECK_AND_FILL_FEATURE(cXYMFT, mftprop.getCovariances()(0, 1));
      CHECK_AND_FILL_FEATURE(cYYMFT, mftprop.getCovariances()(1, 1));
      CHECK_AND_FILL_FEATURE(cPhiXMFT, mftprop.getCovariances()(0, 2));
      CHECK_AND_FILL_FEATURE(cPhiYMFT, mftprop.getCovariances()(1, 2));
      CHECK_AND_FILL_FEATURE(cPhiPhiMFT, mftprop.getCovariances()(2, 2));
      CHECK_AND_FILL_FEATURE(cTglXMFT, mftprop.getCovariances()(0, 3));
      CHECK_AND_FILL_FEATURE(cTglYMFT, mftprop.getCovariances()(1, 3));
      CHECK_AND_FILL_FEATURE(cTglPhiMFT, mftprop.getCovariances()(2, 3));
      CHECK_AND_FILL_FEATURE(cTglTglMFT, mftprop.getCovariances()(3, 3));
      CHECK_AND_FILL_FEATURE(c1PtXMFT, mftprop.getCovariances()(0, 4));
      CHECK_AND_FILL_FEATURE(c1PtYMFT, mftprop.getCovariances()(1, 4));
      CHECK_AND_FILL_FEATURE(c1PtPhiMFT, mftprop.getCovariances()(2, 4));
      CHECK_AND_FILL_FEATURE(c1PtTglMFT, mftprop.getCovariances()(3, 4));
      CHECK_AND_FILL_FEATURE(c1Pt21Pt2MFT, mftprop.getCovariances()(4, 4));
      CHECK_AND_FILL_FEATURE(cXXMCH, mchprop.getCovariances()(0, 0));
      CHECK_AND_FILL_FEATURE(cXYMCH, mchprop.getCovariances()(0, 1));
      CHECK_AND_FILL_FEATURE(cYYMCH, mchprop.getCovariances()(1, 1));
      CHECK_AND_FILL_FEATURE(cPhiXMCH, mchprop.getCovariances()(0, 2));
      CHECK_AND_FILL_FEATURE(cPhiYMCH, mchprop.getCovariances()(1, 2));
      CHECK_AND_FILL_FEATURE(cPhiPhiMCH, mchprop.getCovariances()(2, 2));
      CHECK_AND_FILL_FEATURE(cTglXMCH, mchprop.getCovariances()(0, 3));
      CHECK_AND_FILL_FEATURE(cTglYMCH, mchprop.getCovariances()(1, 3));
      CHECK_AND_FILL_FEATURE(cTglPhiMCH, mchprop.getCovariances()(2, 3));
      CHECK_AND_FILL_FEATURE(cTglTglMCH, mchprop.getCovariances()(3, 3));
      CHECK_AND_FILL_FEATURE(c1PtXMCH, mchprop.getCovariances()(0, 4));
      CHECK_AND_FILL_FEATURE(c1PtYMCH, mchprop.getCovariances()(1, 4));
      CHECK_AND_FILL_FEATURE(c1PtPhiMCH, mchprop.getCovariances()(2, 4));
      CHECK_AND_FILL_FEATURE(c1PtTglMCH, mchprop.getCovariances()(3, 4));
      CHECK_AND_FILL_FEATURE(c1Pt21Pt2MCH, mchprop.getCovariances()(4, 4));
      CHECK_AND_FILL_FEATURE(posX, collision.posX());
      CHECK_AND_FILL_FEATURE(posY, collision.posY());
      CHECK_AND_FILL_FEATURE(posZ, collision.posZ());
      CHECK_AND_FILL_FEATURE(numContrib, collision.numContrib());
      CHECK_AND_FILL_FEATURE(trackOccupancyInTimeRange, collision.trackOccupancyInTimeRange());
      CHECK_AND_FILL_FEATURE(ft0cOccupancyInTimeRange, collision.ft0cOccupancyInTimeRange());
      CHECK_AND_FILL_FEATURE(multFT0A, collision.multFT0A());
      CHECK_AND_FILL_FEATURE(multFT0C, collision.multFT0C());
      CHECK_AND_FILL_FEATURE(multNTracksPV, collision.multNTracksPV());
      CHECK_AND_FILL_FEATURE(multNTracksPVeta1, collision.multNTracksPVeta1());
      CHECK_AND_FILL_FEATURE(multNTracksPVetaHalf, collision.multNTracksPVetaHalf());
      CHECK_AND_FILL_FEATURE(isInelGt0, collision.isInelGt0());
      CHECK_AND_FILL_FEATURE(isInelGt1, collision.isInelGt1());
      CHECK_AND_FILL_FEATURE(multFT0M, collision.multFT0M());
      CHECK_AND_FILL_FEATURE(centFT0M, collision.centFT0M());
      CHECK_AND_FILL_FEATURE(centFT0A, collision.centFT0A());
      CHECK_AND_FILL_FEATURE(centFT0C, collision.centFT0C());
      CHECK_AND_FILL_FEATURE(chi2MCHMFT, muon.chi2MatchMCHMFT());
    }
    return inputFeature;
  }

  /// Method to get the input features vector needed for ML inference
  /// \param track is the single track, \param collision is the collision
  /// \return inputFeatures vector
  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename U>
  std::vector<float> getInputFeatures(T1 const& muon, T2 const& mft, T3 const& mch, T4 const& mftprop, T5 const& mchprop, U const& collision)
  {
    std::vector<float> inputFeatures;
    for (const auto& idx : MlResponse<TypeOutputScore>::mCachedIndices) {
      float inputFeature = returnFeature(idx, muon, mft, mch, mftprop, mchprop, collision);
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
