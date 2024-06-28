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

//
// Class for dimuon selection
//

#ifndef PWGEM_PHOTONMESON_CORE_DIMUONCUT_H_
#define PWGEM_PHOTONMESON_CORE_DIMUONCUT_H_

#include <algorithm>
#include <set>
#include <vector>
#include <utility>
#include <string>
#include "TNamed.h"
#include "Math/Vector4D.h"

#include "Framework/Logger.h"
#include "Framework/DataTypes.h"
#include "CommonConstants/PhysicsConstants.h"

class DimuonCut : public TNamed
{
 public:
  DimuonCut() = default;
  DimuonCut(const char* name, const char* title) : TNamed(name, title) {}

  enum class DimuonCuts : int {
    // pair cut
    kMass = 0,
    kPairPtRange,
    kPairYRange,
    kPairDCARange,
    // track cut
    kTrackType,
    kTrackPtRange,
    kTrackEtaRange,
    kDCAxy,
    kMFTNCls,
    kMCHMIDNCls,
    kChi2,
    kMatchingChi2MCHMFT,
    kMatchingChi2MCHMID,
    kRabs,
    kPDCA,
    kNCuts
  };

  template <typename TPair>
  bool IsSelected(TPair const& pair) const
  {
    auto t1 = std::get<0>(pair);
    auto t2 = std::get<1>(pair);

    if (!IsSelectedTrack(t1) || !IsSelectedTrack(t2)) {
      return false;
    }

    if (!IsSelectedPair(t1, t2)) {
      return false;
    }

    return true;
  }

  template <typename TTrack1, typename TTrack2>
  bool IsSelectedPair(TTrack1 const& t1, TTrack2 const& t2) const
  {
    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassMuon);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    float dcaX1 = t1.fwdDcaX();
    float dcaY1 = t1.fwdDcaY();
    float cXX1 = t1.cXX();
    float cYY1 = t1.cYY();
    float cXY1 = t1.cXY();
    float dcaX2 = t2.fwdDcaX();
    float dcaY2 = t2.fwdDcaY();
    float cXX2 = t2.cXX();
    float cYY2 = t2.cYY();
    float cXY2 = t2.cXY();

    float dca_xy1 = 999.f;
    float det1 = cXX1 * cYY1 - cXY1 * cXY1;
    if (det1 < 0) {
      dca_xy1 = 999.f;
    } else {
      float chi2 = (dcaX1 * dcaX1 * cYY1 + dcaY1 * dcaY1 * cXX1 - 2. * dcaX1 * dcaY1 * cXY1) / det1;
      dca_xy1 = std::sqrt(std::abs(chi2) / 2.); // in sigma
    }

    float dca_xy2 = 999.f;
    float det2 = cXX2 * cYY2 - cXY2 * cXY2;
    if (det2 < 0) {
      dca_xy2 = 999.f;
    } else {
      float chi2 = (dcaX2 * dcaX2 * cYY2 + dcaY2 * dcaY2 * cXX2 - 2. * dcaX2 * dcaY2 * cXY2) / det2;
      dca_xy2 = std::sqrt(std::abs(chi2) / 2.); // in sigma
    }

    float pair_dca_xy = std::sqrt((std::pow(dca_xy1, 2) + std::pow(dca_xy2, 2)) / 2.);

    if (v12.M() < mMinMass || mMaxMass < v12.M()) {
      return false;
    }

    if (v12.Pt() < mMinPairPt || mMaxPairPt < v12.Pt()) {
      return false;
    }

    if (v12.Rapidity() < mMinPairY || mMaxPairY < v12.Rapidity()) {
      return false;
    }

    if (pair_dca_xy < mMinPairDCAxy || mMaxPairDCAxy < pair_dca_xy) { // in sigma for pair
      return false;
    }
    return true;
  }

  template <typename TTrack>
  bool IsSelectedTrack(TTrack const& track) const
  {
    if (!IsSelectedTrack(track, DimuonCuts::kTrackType)) {
      return false;
    }
    if (!IsSelectedTrack(track, DimuonCuts::kTrackPtRange)) {
      return false;
    }
    if (!IsSelectedTrack(track, DimuonCuts::kTrackEtaRange)) {
      return false;
    }
    if (!IsSelectedTrack(track, DimuonCuts::kDCAxy)) {
      return false;
    }
    if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) && !IsSelectedTrack(track, DimuonCuts::kMFTNCls)) {
      return false;
    }
    if (!IsSelectedTrack(track, DimuonCuts::kMCHMIDNCls)) {
      return false;
    }
    if (!IsSelectedTrack(track, DimuonCuts::kChi2)) {
      return false;
    }
    if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) && !IsSelectedTrack(track, DimuonCuts::kMatchingChi2MCHMFT)) {
      return false;
    }
    if (!IsSelectedTrack(track, DimuonCuts::kMatchingChi2MCHMID)) {
      return false;
    }
    if (!IsSelectedTrack(track, DimuonCuts::kPDCA)) {
      return false;
    }
    if (!IsSelectedTrack(track, DimuonCuts::kRabs)) {
      return false;
    }

    return true;
  }

  template <typename T>
  bool IsSelectedTrack(T const& track, const DimuonCuts& cut) const
  {
    switch (cut) {
      case DimuonCuts::kTrackType:
        return track.trackType() == mTrackType;

      case DimuonCuts::kTrackPtRange:
        return track.pt() > mMinTrackPt && track.pt() < mMaxTrackPt;

      case DimuonCuts::kTrackEtaRange:
        return track.eta() > mMinTrackEta && track.eta() < mMaxTrackEta;

      case DimuonCuts::kDCAxy:
        return mMinDcaXY < std::sqrt(std::pow(track.fwdDcaX(), 2) + std::pow(track.fwdDcaY(), 2)) && std::sqrt(std::pow(track.fwdDcaX(), 2) + std::pow(track.fwdDcaY(), 2)) < mMaxDcaXY;

      case DimuonCuts::kMFTNCls:
        return track.nClustersMFT() >= mMinNClustersMFT;

      case DimuonCuts::kMCHMIDNCls:
        return track.nClusters() >= mMinNClustersMCHMID;

      case DimuonCuts::kChi2:
        return track.chi2() < mMaxChi2;

      case DimuonCuts::kMatchingChi2MCHMFT:
        return track.chi2MatchMCHMFT() < mMaxMatchingChi2MCHMFT;

      case DimuonCuts::kMatchingChi2MCHMID:
        return track.chi2MatchMCHMID() < mMaxMatchingChi2MCHMID;

      case DimuonCuts::kPDCA:
        return track.pDca() < mMaxPDCARabsDep(track.rAtAbsorberEnd());

      case DimuonCuts::kRabs:
        return mMinRabs < track.rAtAbsorberEnd() && track.rAtAbsorberEnd() < mMaxRabs;

      default:
        return false;
    }
  }

  // Setters
  void SetMassRange(float min = 0.f, float max = 1e+10);
  void SetPairPtRange(float minPt = 0.f, float maxPt = 1e10f);
  void SetPairYRange(float minY = -1e10f, float maxY = 1e10f);
  void SetPairDCAxyRange(float min = 0.f, float max = 1e10f); // DCAxy in cm

  void SetTrackType(int track_type); // 0: MFT-MCH-MID (global muon), 3: MCH-MID (standalone muon)
  void SetTrackPtRange(float minPt = 0.f, float maxPt = 1e10f);
  void SetTrackEtaRange(float minEta = -1e10f, float maxEta = 1e10f);
  void SetNClustersMFT(int min, int max);
  void SetNClustersMCHMID(int min, int max);
  void SetChi2(float min, float max);
  void SetMatchingChi2MCHMFT(float min, float max);
  void SetMatchingChi2MCHMID(float min, float max);
  void SetDCAxy(float min, float max); // in cm
  void SetRabs(float min, float max);  // in cm
  void SetMaxPDCARabsDep(std::function<float(float)> RabsDepCut);

 private:
  // pair cuts
  float mMinMass{0.f}, mMaxMass{1e10f};
  float mMinPairPt{0.f}, mMaxPairPt{1e10f};       // range in pT
  float mMinPairY{-1e10f}, mMaxPairY{1e10f};      // range in rapidity
  float mMinPairDCAxy{0.f}, mMaxPairDCAxy{1e10f}; // range in 3D DCA in sigma

  // kinematic cuts
  float mMinTrackPt{0.f}, mMaxTrackPt{1e10f};      // range in pT
  float mMinTrackEta{-1e10f}, mMaxTrackEta{1e10f}; // range in eta

  // track quality cuts
  int mTrackType{3};
  int mMinNClustersMFT{0}, mMaxNClustersMFT{10};                    // min number of TPC clusters
  int mMinNClustersMCHMID{0}, mMaxNClustersMCHMID{16};              // min number of TPC clusters
  float mMinChi2{0.f}, mMaxChi2{1e10f};                             // max tpc fit chi2 per TPC cluster
  float mMinMatchingChi2MCHMFT{0.f}, mMaxMatchingChi2MCHMFT{1e10f}; // max tpc fit chi2 per TPC cluster
  float mMinMatchingChi2MCHMID{0.f}, mMaxMatchingChi2MCHMID{1e10f}; // max tpc fit chi2 per TPC cluster
  std::function<float(float)> mMaxPDCARabsDep{};                    // max pdca in xy plane as function of Rabs

  float mMinRabs{17.6}, mMaxRabs{89.5};
  float mMinDcaXY{0.0f}, mMaxDcaXY{1e10f};

  ClassDef(DimuonCut, 1);
};

#endif // PWGEM_PHOTONMESON_CORE_DIMUONCUT_H_
