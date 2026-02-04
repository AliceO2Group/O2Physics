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

#ifndef PWGEM_DILEPTON_CORE_DIMUONCUT_H_
#define PWGEM_DILEPTON_CORE_DIMUONCUT_H_

#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/DataTypes.h"
#include "Framework/Logger.h"
#include "MathUtils/Utils.h"

#include "Math/Vector4D.h"
#include "TNamed.h"

#include <algorithm>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;

class DimuonCut : public TNamed
{
 public:
  DimuonCut() = default;
  DimuonCut(const char* name, const char* title) : TNamed(name, title) {}

  ~DimuonCut() {}

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
    kTrackPhiRange,
    kDCAxy,
    kMFTNCls,
    kMCHMIDNCls,
    kChi2,
    kMatchingChi2MCHMFT,
    kMatchingChi2MCHMID,
    kChi2MFT,
    kRabs,
    kPDCA,
    kMFTHitMap,
    kDPtDEtaDPhiwrtMCHMID,
    kDr_MatchingChi2MCHMFT_2D,
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

  template <bool dont_require_rapidity = false, typename TTrack1, typename TTrack2>
  bool IsSelectedPair(TTrack1 const& t1, TTrack2 const& t2) const
  {
    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassMuon);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    float dca_xy_t1 = fwdDcaXYinSigma(t1);
    float dca_xy_t2 = fwdDcaXYinSigma(t2);
    float pair_dca_xy = std::sqrt((dca_xy_t1 * dca_xy_t1 + dca_xy_t2 * dca_xy_t2) / 2.);

    if (v12.M() < mMinMass || mMaxMass < v12.M()) {
      return false;
    }

    if (v12.Pt() < mMinPairPt || mMaxPairPt < v12.Pt()) {
      return false;
    }

    if (!dont_require_rapidity && (v12.Rapidity() < mMinPairY || mMaxPairY < v12.Rapidity())) {
      return false;
    }

    if (pair_dca_xy < mMinPairDCAxy || mMaxPairDCAxy < pair_dca_xy) { // in sigma for pair
      return false;
    }

    float deta = v1.Eta() - v2.Eta();
    float dphi = v1.Phi() - v2.Phi();
    o2::math_utils::bringToPMPi(dphi);
    if (mApplydEtadPhi && std::pow(deta / mMinDeltaEta, 2) + std::pow(dphi / mMinDeltaPhi, 2) < 1.f) {
      return false;
    }

    return true;
  }

  template <bool dont_require_pteta = false, typename TTrack>
  bool IsSelectedTrack(TTrack const& track) const
  {
    if (!IsSelectedTrack(track, DimuonCuts::kTrackType)) {
      return false;
    }

    if (!dont_require_pteta) {
      if (!IsSelectedTrack(track, DimuonCuts::kTrackPtRange)) {
        return false;
      }
      if (!IsSelectedTrack(track, DimuonCuts::kTrackEtaRange)) {
        return false;
      }
    }
    if (!IsSelectedTrack(track, DimuonCuts::kTrackPhiRange)) {
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
    if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) && !IsSelectedTrack(track, DimuonCuts::kChi2MFT)) {
      return false;
    }
    if (!IsSelectedTrack(track, DimuonCuts::kPDCA)) {
      return false;
    }
    if (!IsSelectedTrack(track, DimuonCuts::kRabs)) {
      return false;
    }
    if (mApplyMFTHitMap && track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) && !IsSelectedTrack(track, DimuonCuts::kMFTHitMap)) {
      return false;
    }
    if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) && !IsSelectedTrack(track, DimuonCuts::kDPtDEtaDPhiwrtMCHMID)) {
      return false;
    }
    if (track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) && !IsSelectedTrack(track, DimuonCuts::kDr_MatchingChi2MCHMFT_2D)) {
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

      case DimuonCuts::kTrackPhiRange:
        return track.phi() > mMinTrackPhi && track.phi() < mMaxTrackPhi;

      case DimuonCuts::kDCAxy:
        return mMinDcaXY < std::sqrt(std::pow(track.fwdDcaX(), 2) + std::pow(track.fwdDcaY(), 2)) && std::sqrt(std::pow(track.fwdDcaX(), 2) + std::pow(track.fwdDcaY(), 2)) < mMaxDcaXY;

      case DimuonCuts::kMFTNCls:
        return track.nClustersMFT() >= mMinNClustersMFT;

      case DimuonCuts::kMCHMIDNCls:
        return track.nClusters() >= mMinNClustersMCHMID;

      case DimuonCuts::kChi2:
        return track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) ? track.chi2() / (2.f * (track.nClusters() + track.nClustersMFT()) - 5.f) < mMaxChi2 : track.chi2() < mMaxChi2;

      case DimuonCuts::kChi2MFT:
        return track.trackType() == static_cast<uint8_t>(o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) ? track.chi2MFT() / (2.f * track.nClustersMFT() - 5.f) < mMaxChi2MFT : true;

      case DimuonCuts::kMatchingChi2MCHMFT:
        return track.chi2MatchMCHMFT() < mMaxMatchingChi2MCHMFT;

      case DimuonCuts::kMatchingChi2MCHMID:
        return track.chi2MatchMCHMID() < mMaxMatchingChi2MCHMID;

      case DimuonCuts::kPDCA:
        return track.pDca() < mMaxPDCARabsDep(track.rAtAbsorberEnd());

      case DimuonCuts::kRabs:
        return mMinRabs < track.rAtAbsorberEnd() && track.rAtAbsorberEnd() < mMaxRabs;

      case DimuonCuts::kMFTHitMap: {
        std::vector<bool> mftHitMap{checkMFTHitMap<0, 1>(track), checkMFTHitMap<2, 3>(track), checkMFTHitMap<4, 5>(track), checkMFTHitMap<6, 7>(track), checkMFTHitMap<8, 9>(track)};
        for (const auto& iDisk : mRequiredMFTDisks) {
          if (!mftHitMap[iDisk]) {
            return false;
          }
        }
        return true;
      }

      case DimuonCuts::kDPtDEtaDPhiwrtMCHMID:
        return std::fabs(track.ptMatchedMCHMID() - track.pt()) / track.pt() < mMaxReldPtwrtMCHMID && std::sqrt(std::pow((track.etaMatchedMCHMID() - track.eta()) / mMaxdEtawrtMCHMID, 2) + std::pow((track.phiMatchedMCHMID() - track.phi()) / mMaxdPhiwrtMCHMID, 2)) < 1.f;

      case DimuonCuts::kDr_MatchingChi2MCHMFT_2D:
        return mSlope_dr_chi2MatchMFTMCH * track.chi2MatchMCHMFT() + mIntercept_dr_chi2MatchMFTMCH > std::sqrt(std::pow(track.etaMatchedMCHMID() - track.eta(), 2) + std::pow(track.phiMatchedMCHMID() - track.phi(), 2));

      default:
        return false;
    }
  }

  // Setters
  void SetMassRange(float min = 0.f, float max = 1e+10);
  void SetPairPtRange(float minPt = 0.f, float maxPt = 1e10f);
  void SetPairYRange(float minY = -1e10f, float maxY = 1e10f);
  void SetPairDCAxyRange(float min = 0.f, float max = 1e10f); // DCAxy in cm
  void SetMindEtadPhi(bool flag, float min_deta, float min_dphi);

  void SetTrackType(int track_type); // 0: MFT-MCH-MID (global muon), 3: MCH-MID (standalone muon)
  void SetTrackPtRange(float minPt = 0.f, float maxPt = 1e10f);
  void SetTrackEtaRange(float minEta = -1e10f, float maxEta = 1e10f);
  void SetTrackPhiRange(float minPhi = 0.f, float maxPhi = 2.f * M_PI);
  void SetNClustersMFT(int min, int max);
  void SetNClustersMCHMID(int min, int max);
  void SetChi2(float min, float max);
  void SetChi2MFT(float min, float max);
  void SetMatchingChi2MCHMFT(float min, float max);
  void SetMatchingChi2MCHMID(float min, float max);
  void SetDCAxy(float min, float max); // in cm
  void SetRabs(float min, float max);  // in cm
  void SetMaxPDCARabsDep(std::function<float(float)> RabsDepCut);
  void SetMFTHitMap(bool flag, std::vector<int> hitMap);
  void SetMaxdPtdEtadPhiwrtMCHMID(float reldPtMax, float dEtaMax, float dPhiMax); // this is relevant for global muons
  void SetSlopeAndInterceptDRvsChi2MCHMFT(float slope, float intercept);          // this is relevant for global muons

 private:
  // pair cuts
  float mMinMass{0.f}, mMaxMass{1e10f};
  float mMinPairPt{0.f}, mMaxPairPt{1e10f};       // range in pT
  float mMinPairY{-1e10f}, mMaxPairY{1e10f};      // range in rapidity
  float mMinPairDCAxy{0.f}, mMaxPairDCAxy{1e10f}; // range in 3D DCA in sigma
  bool mApplydEtadPhi{false};                     // flag to apply deta, dphi cut between 2 tracks
  float mMinDeltaEta{0.f};
  float mMinDeltaPhi{0.f};

  // kinematic cuts
  float mMinTrackPt{0.f}, mMaxTrackPt{1e10f};        // range in pT
  float mMinTrackEta{-1e10f}, mMaxTrackEta{1e10f};   // range in eta
  float mMinTrackPhi{0.f}, mMaxTrackPhi{2.f * M_PI}; // range in phi

  // track quality cuts
  int mTrackType{3};
  int mMinNClustersMFT{0}, mMaxNClustersMFT{10};                    // min number of MFT clusters
  int mMinNClustersMCHMID{0}, mMaxNClustersMCHMID{20};              // min number of MCH-MID clusters
  float mMinChi2{0.f}, mMaxChi2{1e10f};                             // max chi2 per MFT + MCH cluster
  float mMinChi2MFT{0.f}, mMaxChi2MFT{1e10f};                       // max chi2 per MFT cluster
  float mMinMatchingChi2MCHMFT{0.f}, mMaxMatchingChi2MCHMFT{1e10f}; // max matching chi2 between MCH-MFT
  float mMinMatchingChi2MCHMID{0.f}, mMaxMatchingChi2MCHMID{1e10f}; // max matching chi2 between MCH-MID
  std::function<float(float)> mMaxPDCARabsDep{};                    // max pdca in xy plane as function of Rabs

  float mMinRabs{17.6}, mMaxRabs{89.5};
  float mMinDcaXY{0.0f}, mMaxDcaXY{1e10f};
  float mMaxReldPtwrtMCHMID{1e10f}, mMaxdEtawrtMCHMID{1e10f}, mMaxdPhiwrtMCHMID{1e10f};
  bool mApplyMFTHitMap{false};
  std::vector<int> mRequiredMFTDisks{};
  float mSlope_dr_chi2MatchMFTMCH{-0.15 / 30};
  float mIntercept_dr_chi2MatchMFTMCH{1e+10};

  ClassDef(DimuonCut, 1);
};

#endif // PWGEM_DILEPTON_CORE_DIMUONCUT_H_
