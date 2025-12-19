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

/// \file V0PhotonCut.h
/// \brief Header of class for V0 photon selection.
/// \author D. Sekihata, daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_CORE_V0PHOTONCUT_H_
#define PWGEM_PHOTONMESON_CORE_V0PHOTONCUT_H_

#include "PWGEM/PhotonMeson/Utils/TrackSelection.h"

#include <Framework/ASoA.h>

#include <TMath.h>
#include <TNamed.h>

#include <Rtypes.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <set>
#include <utility>

using namespace o2::pwgem::photonmeson;

class V0PhotonCut : public TNamed
{
 public:
  V0PhotonCut() = default;
  V0PhotonCut(const char* name, const char* title) : TNamed(name, title) {}

  enum class V0PhotonCuts : int {
    // v0 cut
    kMee = 0,
    kV0PtRange,
    kV0EtaRange,
    kAP,
    kPsiPair,
    kPhiV,
    kRxy,
    kCosPA,
    kPCA,
    kChi2KF,
    kRZLine,
    kOnWwireIB,
    kOnWwireOB,
    // leg cut
    kTrackPtRange,
    kTrackEtaRange,
    kTPCNCls,
    kTPCCrossedRows,
    kTPCCrossedRowsOverNCls,
    kTPCFracSharedClusters,
    kTPCChi2NDF,
    kTPCNsigmaEl,
    kTPCNsigmaPi,
    kDCAxy,
    kDCAz,
    kITSNCls,
    kITSChi2NDF,
    kITSClusterSize,
    kRequireITSTPC,
    kRequireITSonly,
    kRequireTPConly,
    kRequireTPCTRD,
    kRequireTPCTOF,
    kNCuts
  };

  template <o2::soa::is_iterator TV0, typename TLeg>
  bool IsSelected(TV0 const& v0) const
  {
    if (!IsSelectedV0(v0, V0PhotonCuts::kV0PtRange)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kV0EtaRange)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kMee)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kAP)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kPsiPair)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kPhiV)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kRxy)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kCosPA)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kPCA)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kChi2KF)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kRZLine)) {
      return false;
    }

    if (mIsOnWwireIB && mIsOnWwireOB) {
      if (!IsSelectedV0(v0, V0PhotonCuts::kOnWwireIB) && !IsSelectedV0(v0, V0PhotonCuts::kOnWwireOB)) {
        return false;
      }
    } else if (mIsOnWwireIB) {
      if (!IsSelectedV0(v0, V0PhotonCuts::kOnWwireIB)) {
        return false;
      }
    } else if (mIsOnWwireOB) {
      if (!IsSelectedV0(v0, V0PhotonCuts::kOnWwireOB)) {
        return false;
      }
    }

    auto pos = v0.template posTrack_as<TLeg>();
    auto ele = v0.template negTrack_as<TLeg>();

    for (auto& track : {pos, ele}) {
      if (!IsSelectedTrack(track, V0PhotonCuts::kTrackPtRange)) {
        return false;
      }
      if (!IsSelectedTrack(track, V0PhotonCuts::kTrackEtaRange)) {
        return false;
      }
      if (!IsSelectedTrack(track, V0PhotonCuts::kDCAxy)) {
        return false;
      }
      if (!IsSelectedTrack(track, V0PhotonCuts::kDCAz)) {
        return false;
      }
      if (!track.hasITS() && !track.hasTPC()) { // track has to be ITSonly or TPConly or ITS-TPC
        return false;
      }
      if (mDisableITSonly && isITSonlyTrack(track)) {
        return false;
      }
      if (mDisableTPConly && isTPConlyTrack(track)) {
        return false;
      }

      if (mRejectITSib) {
        auto hits_ib = std::count_if(its_ib_Requirement.second.begin(), its_ib_Requirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
        auto hits_ob = std::count_if(its_ob_Requirement.second.begin(), its_ob_Requirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
        bool its_ob_only = (hits_ib <= its_ib_Requirement.first) && (hits_ob >= its_ob_Requirement.first);
        if (isITSonlyTrack(track) && !its_ob_only) { // ITSonly tracks should not have any ITSib hits.
          return false;
        }

        auto hits_ob_itstpc = std::count_if(its_ob_Requirement_ITSTPC.second.begin(), its_ob_Requirement_ITSTPC.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
        bool its_ob_only_itstpc = (hits_ib <= its_ib_Requirement.first) && (hits_ob_itstpc >= its_ob_Requirement_ITSTPC.first);
        if (isITSTPCTrack(track) && !its_ob_only_itstpc) { // ITSTPC tracks should not have any ITSib hits.
          return false;
        }
      }

      if (track.hasITS() && !CheckITSCuts(track)) {
        return false;
      }
      if (track.hasTPC() && !CheckTPCCuts(track)) {
        return false;
      }
      if (track.hasITS() && !track.hasTPC() && (track.hasTRD() || track.hasTOF())) { // remove ITS-TRD, ITS-TOF, ITS-TRD-TOF that are unrealistic tracks.
        return false;
      }

      if (mIsOnWwireIB && !CheckITSCuts(track)) { // photon conversion on ibw requires ITS hits.
        return false;
      }

      if (mRequireITSonly && !IsSelectedTrack(track, V0PhotonCuts::kRequireITSonly)) {
        return false;
      }
      if (mRequireITSTPC && !IsSelectedTrack(track, V0PhotonCuts::kRequireITSTPC)) {
        return false;
      }
      if (mRequireTPConly && !IsSelectedTrack(track, V0PhotonCuts::kRequireTPConly)) {
        return false;
      }
      if (mRequireTPCTRD && !IsSelectedTrack(track, V0PhotonCuts::kRequireTPCTRD)) {
        return false;
      }
      if (mRequireTPCTOF && !IsSelectedTrack(track, V0PhotonCuts::kRequireTPCTOF)) {
        return false;
      }
    }
    return true;
  }

  template <typename T>
  bool CheckITSCuts(T const& track) const
  {
    if (!IsSelectedTrack(track, V0PhotonCuts::kITSNCls)) {
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kITSChi2NDF)) {
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kITSClusterSize)) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool CheckTPCCuts(T const& track) const
  {
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCNCls)) {
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCCrossedRows)) {
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCCrossedRowsOverNCls)) {
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCFracSharedClusters)) {
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCChi2NDF)) {
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCNsigmaEl)) {
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCNsigmaPi)) {
      return false;
    }
    return true;
  }

  template <o2::soa::is_iterator T>
  bool IsSelectedV0(T const& v0, const V0PhotonCuts& cut) const
  {
    switch (cut) {
      case V0PhotonCuts::kV0PtRange:
        return v0.pt() >= mMinV0Pt && v0.pt() <= mMaxV0Pt;

      case V0PhotonCuts::kV0EtaRange:
        return v0.eta() >= mMinV0Eta && v0.eta() <= mMaxV0Eta;

      case V0PhotonCuts::kMee:
        return v0.mGamma() < mMaxQt * 2.f;

      case V0PhotonCuts::kAP:
        return pow(v0.alpha() / mMaxAlpha, 2) + pow(v0.qtarm() / mMaxQt, 2) < 1.0;

      case V0PhotonCuts::kPsiPair:
        return true;

      case V0PhotonCuts::kPhiV:
        return true;

      case V0PhotonCuts::kRxy: {
        if (v0.v0radius() < mMinRxy || mMaxRxy < v0.v0radius()) {
          return false;
        }
        return true;
      }

      case V0PhotonCuts::kCosPA:
        return v0.cospa() >= mMinCosPA;

      case V0PhotonCuts::kPCA:
        return v0.pca() <= mMaxPCA;

      case V0PhotonCuts::kChi2KF:
        return v0.chiSquareNDF() <= mMaxChi2KF;

      case V0PhotonCuts::kRZLine:
        return v0.v0radius() > std::fabs(v0.vz()) * std::tan(2 * std::atan(std::exp(-mMaxV0Eta))) - mMaxMarginZ;

      case V0PhotonCuts::kOnWwireIB: {
        const float margin_xy = 1.0; // cm
        // const float margin_z = 20.0; // cm
        // const float rxy_min = 5.506; // cm
        // const float rxy_max = 14.846;         // cm
        // const float z_min = -17.56; // cm
        // const float z_max = +31.15;           // cm
        float x = std::fabs(v0.vx()); // cm, measured secondary vertex of gamma->ee
        float y = v0.vy();            // cm, measured secondary vertex of gamma->ee
        float z = v0.vz();            // cm, measured secondary vertex of gamma->ee

        float rxy = sqrt(x * x + y * y);
        if (rxy < 7.0 || 14.0 < rxy) {
          return false;
        }

        // r = 0.192 * z + 8.88 (cm) expected wire position in RZ plane.TMath::Tan(10.86 * TMath::DegToRad()) = 0.192
        if (rxy > 0.192 * z + 14.0) { // upper limit
          return false;
        }

        float dxy = std::fabs(1.0 * y - x * std::tan(-8.52 * TMath::DegToRad())) / sqrt(pow(1.0, 2) + pow(std::tan(-8.52 * TMath::DegToRad()), 2));
        return !(dxy > margin_xy);
      }
      case V0PhotonCuts::kOnWwireOB: {
        const float margin_xy = 1.0;                                      // cm
        const float rxy_exp = 30.8;                                       // cm
        const float x_exp = rxy_exp * std::cos(-1.3 * TMath::DegToRad()); // cm, expected position x of W wire
        const float y_exp = rxy_exp * std::sin(-1.3 * TMath::DegToRad()); // cm, expected position y of W wire
        // const float z_min = -47.0;                                          // cm
        // const float z_max = +47.0;                                          // cm
        float x = v0.vx(); // cm, measured secondary vertex of gamma->ee
        float y = v0.vy(); // cm, measured secondary vertex of gamma->ee
        // float z = v0.vz();                                                  // cm, measured secondary vertex of gamma->ee

        // float rxy = sqrt(x * x + y * y);
        // if (rxy < 28.0 || 33.0 < rxy) {
        //   return false;
        // }

        float dxy = std::sqrt(std::pow(x - x_exp, 2) + std::pow(y - y_exp, 2));
        return !(dxy > margin_xy);
      }
      default:
        return false;
    }
  }

  // Temporary function to check if track passes a given selection criteria. To be replaced by framework filters.
  template <typename T>
  bool IsSelectedTrack(T const& track, const V0PhotonCuts& cut) const
  {
    switch (cut) {
      case V0PhotonCuts::kTrackPtRange:
        return track.pt() > mMinTrackPt && track.pt() < mMaxTrackPt;

      case V0PhotonCuts::kTrackEtaRange:
        return track.eta() > mMinTrackEta && track.eta() < mMaxTrackEta;

      case V0PhotonCuts::kTPCNCls:
        return track.tpcNClsFound() >= mMinNClustersTPC;

      case V0PhotonCuts::kTPCCrossedRows:
        return track.tpcNClsCrossedRows() >= mMinNCrossedRowsTPC;

      case V0PhotonCuts::kTPCCrossedRowsOverNCls:
        return track.tpcCrossedRowsOverFindableCls() >= mMinNCrossedRowsOverFindableClustersTPC;

      case V0PhotonCuts::kTPCFracSharedClusters:
        return track.tpcFractionSharedCls() < mMaxFracSharedClustersTPC;

      case V0PhotonCuts::kTPCChi2NDF:
        return mMinChi2PerClusterTPC < track.tpcChi2NCl() && track.tpcChi2NCl() < mMaxChi2PerClusterTPC;

      case V0PhotonCuts::kTPCNsigmaEl:
        return track.tpcNSigmaEl() > mMinTPCNsigmaEl && track.tpcNSigmaEl() < mMaxTPCNsigmaEl;

      case V0PhotonCuts::kTPCNsigmaPi:
        return track.tpcNSigmaPi() > mMinTPCNsigmaPi && track.tpcNSigmaPi() < mMaxTPCNsigmaPi;

      case V0PhotonCuts::kDCAxy:
        return std::fabs(track.dcaXY()) < ((mMaxDcaXYPtDep) ? mMaxDcaXYPtDep(track.pt()) : mMaxDcaXY);

      case V0PhotonCuts::kDCAz:
        return std::fabs(track.dcaZ()) < mMaxDcaZ;

      case V0PhotonCuts::kITSNCls:
        return mMinNClustersITS <= track.itsNCls() && track.itsNCls() <= mMaxNClustersITS;

      case V0PhotonCuts::kITSChi2NDF:
        return mMinChi2PerClusterITS < track.itsChi2NCl() && track.itsChi2NCl() < mMaxChi2PerClusterITS;

      case V0PhotonCuts::kITSClusterSize: {
        if (!isITSonlyTrack(track)) {
          return true;
        }
        return mMinMeanClusterSizeITS < track.meanClusterSizeITSob() * std::cos(std::atan(track.tgl())) && track.meanClusterSizeITSob() * std::cos(std::atan(track.tgl())) < mMaxMeanClusterSizeITS;
      }

      case V0PhotonCuts::kRequireITSTPC:
        return isITSTPCTrack(track);

      case V0PhotonCuts::kRequireITSonly:
        return isITSonlyTrack(track);

      case V0PhotonCuts::kRequireTPConly:
        return isTPConlyTrack(track);

      case V0PhotonCuts::kRequireTPCTRD:
        return isTPCTRDTrack(track);

      case V0PhotonCuts::kRequireTPCTOF:
        return isTPCTOFTrack(track);

      default:
        return false;
    }
  }

  // Setters
  void SetV0PtRange(float minPt = 0.f, float maxPt = 1e10f);
  void SetV0EtaRange(float minEta = -1e10f, float maxEta = 1e10f);
  void SetMeeRange(float min = 0.f, float max = 0.1);
  void SetPsiPairRange(float min = -3.15, float max = +3.15);
  void SetPhivPairRange(float min = 0.f, float max = +3.15);
  void SetAPRange(float max_alpha = 0.95, float max_qt = 0.05); // Armenteros Podolanski
  void SetRxyRange(float min = 0.f, float max = 180.f);
  void SetMinCosPA(float min = 0.95);
  void SetMaxPCA(float max = 2.f);
  void SetMaxChi2KF(float max = 1e+10);
  void SetMaxMarginZ(float max = 7.f);
  void SetMaxMeePsiPairDep(std::function<float(float)> psiDepCut);
  void SetOnWwireIB(bool flag = false);
  void SetOnWwireOB(bool flag = false);
  void RejectITSib(bool flag = false);

  void SetTrackPtRange(float minPt = 0.f, float maxPt = 1e10f);
  void SetTrackEtaRange(float minEta = -1e10f, float maxEta = 1e10f);
  void SetMinNClustersTPC(int minNClustersTPC);
  void SetMinNCrossedRowsTPC(int minNCrossedRowsTPC);
  void SetMinNCrossedRowsOverFindableClustersTPC(float minNCrossedRowsOverFindableClustersTPC);
  void SetMaxFracSharedClustersTPC(float max);
  void SetChi2PerClusterTPC(float min, float max);
  void SetNClustersITS(int min, int max);
  void SetChi2PerClusterITS(float min, float max);
  void SetMeanClusterSizeITSob(float min, float max);

  void SetTPCNsigmaElRange(float min = -3, float max = +3);
  void SetTPCNsigmaPiRange(float min = -1e+10, float max = 1e+10);

  void SetMaxDcaXY(float maxDcaXY);
  void SetMaxDcaZ(float maxDcaZ);
  void SetMaxDcaXYPtDep(std::function<float(float)> ptDepCut);
  void SetIsWithinBeamPipe(bool flag);
  void SetRequireITSTPC(bool flag);
  void SetRequireITSonly(bool flag);
  void SetRequireTPConly(bool flag);
  void SetRequireTPCTRD(bool flag);
  void SetRequireTPCTOF(bool flag);
  void SetDisableITSonly(bool flag);
  void SetDisableTPConly(bool flag);

 private:
  static const std::pair<int8_t, std::set<uint8_t>> its_ib_Requirement;
  static const std::pair<int8_t, std::set<uint8_t>> its_ob_Requirement;
  static const std::pair<int8_t, std::set<uint8_t>> its_ob_Requirement_ITSTPC;
  // v0 cuts
  float mMinMee{0.f}, mMaxMee{0.1f};
  float mMinV0Pt{0.f}, mMaxV0Pt{1e10f};      // range in pT
  float mMinV0Eta{-1e10f}, mMaxV0Eta{1e10f}; // range in eta
  float mMaxAlpha{0.95}, mMaxQt{0.05};
  float mMinPsiPair{-3.15}, mMaxPsiPair{+3.15};
  float mMinPhivPair{0.f}, mMaxPhivPair{+3.15};
  float mMinRxy{0.f}, mMaxRxy{180.f};
  float mMinCosPA{0.95};
  float mMaxPCA{2.f};
  float mMaxChi2KF{1e+10};
  float mMaxMarginZ{7.f};
  std::function<float(float)> mMaxMeePsiPairDep{}; // max mee as a function of psipair
  bool mIsOnWwireIB{false};
  bool mIsOnWwireOB{false};
  bool mRejectITSib{false};

  // pid cuts
  float mMinTPCNsigmaEl{-5}, mMaxTPCNsigmaEl{+5};
  float mMinTPCNsigmaPi{-1e+10}, mMaxTPCNsigmaPi{+1e+10};

  // kinematic cuts
  float mMinTrackPt{0.f}, mMaxTrackPt{1e10f};      // range in pT
  float mMinTrackEta{-1e10f}, mMaxTrackEta{1e10f}; // range in eta

  // track quality cuts
  int mMinNClustersTPC{0};                                             // min number of TPC clusters
  int mMinNCrossedRowsTPC{0};                                          // min number of crossed rows in TPC
  float mMinChi2PerClusterTPC{-1e10f}, mMaxChi2PerClusterTPC{1e10f};   // max tpc fit chi2 per TPC cluster
  float mMinNCrossedRowsOverFindableClustersTPC{0.f};                  // min ratio crossed rows / findable clusters
  float mMaxFracSharedClustersTPC{999.f};                              // max ratio shared clusters / clusters in TPC
  int mMinNClustersITS{0}, mMaxNClustersITS{7};                        // range in number of ITS clusters
  float mMinChi2PerClusterITS{-1e10f}, mMaxChi2PerClusterITS{1e10f};   // max its fit chi2 per ITS cluster
  float mMinMeanClusterSizeITS{-1e10f}, mMaxMeanClusterSizeITS{1e10f}; // max <its cluster size> x cos(Lmabda)

  float mMaxDcaXY{1e10f};                       // max dca in xy plane
  float mMaxDcaZ{1e10f};                        // max dca in z direction
  std::function<float(float)> mMaxDcaXYPtDep{}; // max dca in xy plane as function of pT
  bool mRequireITSTPC{false};
  bool mRequireITSonly{false};
  bool mRequireTPConly{false};
  bool mRequireTPCTRD{false};
  bool mRequireTPCTOF{false};
  bool mDisableITSonly{false};
  bool mDisableTPConly{false};

  ClassDef(V0PhotonCut, 2);
};

#endif // PWGEM_PHOTONMESON_CORE_V0PHOTONCUT_H_
