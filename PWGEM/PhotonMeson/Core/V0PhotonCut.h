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
// Class for v0 photon selection
//

#ifndef PWGEM_PHOTONMESON_CORE_V0PHOTONCUT_H_
#define PWGEM_PHOTONMESON_CORE_V0PHOTONCUT_H_

#include <algorithm>
#include <set>
#include <vector>
#include <utility>
#include <string>
#include "Rtypes.h"
#include "TNamed.h"
#include "TMath.h"

#include "PWGEM/PhotonMeson/Utils/TrackSelection.h"
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
    kRZLine,
    kOnWwireIB,
    kOnWwireOB,
    // leg cut
    kTrackPtRange,
    kTrackEtaRange,
    kTPCNCls,
    kTPCCrossedRows,
    kTPCCrossedRowsOverNCls,
    kTPCChi2NDF,
    kTPCNsigmaEl,
    kTPCNsigmaPi,
    kDCAxy,
    kDCAz,
    kITSNCls,
    kITSChi2NDF,
    kIsWithinBeamPipe,
    kRequireITSTPC,
    kRequireITSonly,
    kRequireTPConly,
    kRequireTPCTRD,
    kRequireTPCTOF,
    kRequireTPCTRDTOF,
    kNCuts
  };

  static const char* mCutNames[static_cast<int>(V0PhotonCuts::kNCuts)];

  template <class TLeg, typename TV0>
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

    if (pos.hasITS() && ele.hasITS()) {
      if (v0.mGammaKFSV() > 0.06 && v0.recalculatedVtxR() < 12.f) {
        return false;
      }
    }

    bool isTPConly_pos = isTPConlyTrack(pos);
    bool isTPConly_ele = isTPConlyTrack(ele);
    if (isTPConly_pos && isTPConly_ele) {
      if (v0.mGammaKFSV() > v0.mGammaKFPV()) {
        return false;
      }
    }

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
      if (mIsWithinBP && !IsSelectedTrack(track, V0PhotonCuts::kIsWithinBeamPipe)) {
        return false;
      }
      if (!track.hasITS() && !track.hasTPC()) { // track has to be ITSonly or TPConly or ITS-TPC
        return false;
      }

      bool isITSonly = isITSonlyTrack(track);
      auto hits_ib = std::count_if(its_ib_Requirement.second.begin(), its_ib_Requirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
      auto hits_ob = std::count_if(its_ob_Requirement.second.begin(), its_ob_Requirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
      bool its_ob_only = (hits_ib <= its_ib_Requirement.first) && (hits_ob >= its_ob_Requirement.first);
      if (isITSonly && !its_ob_only) { // ITSonly tracks should not have any ITSib hits.
        return false;
      }

      bool isITSTPC = isITSTPCTrack(track);
      auto hits_ob_itstpc = std::count_if(its_ob_Requirement_ITSTPC.second.begin(), its_ob_Requirement_ITSTPC.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
      bool its_ob_only_itstpc = (hits_ib <= its_ib_Requirement.first) && (hits_ob_itstpc >= its_ob_Requirement_ITSTPC.first);
      if (isITSTPC && !its_ob_only_itstpc) { // ITSonly tracks should not have any ITSib hits.
        return false;
      }

      if (isITSonly) {
        if (!CheckITSCuts(track)) {
          return false;
        }
      } else if (track.hasTPC()) {
        if (!CheckTPCCuts(track)) {
          return false;
        }
      } else { // remove ITS-TRD, ITS-TOF, ITS-TRD-TOF that are unrealistic tracks.
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
      if (mRequireTPCTRDTOF && !IsSelectedTrack(track, V0PhotonCuts::kRequireTPCTRDTOF)) {
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

  template <typename T>
  uint32_t IsSelectedMask(T const& track) const
  {
    uint32_t flag = 0;

    auto setFlag = [&](const V0PhotonCuts& cut) {
      if (IsSelectedTrack(track, cut)) {
        flag |= 1UL << static_cast<int>(cut);
      }
    };

    setFlag(V0PhotonCuts::kV0PtRange);
    setFlag(V0PhotonCuts::kV0EtaRange);
    setFlag(V0PhotonCuts::kTrackPtRange);
    setFlag(V0PhotonCuts::kTrackEtaRange);
    setFlag(V0PhotonCuts::kTPCNCls);
    setFlag(V0PhotonCuts::kTPCCrossedRows);
    setFlag(V0PhotonCuts::kTPCCrossedRowsOverNCls);
    setFlag(V0PhotonCuts::kTPCChi2NDF);
    setFlag(V0PhotonCuts::kDCAxy);
    setFlag(V0PhotonCuts::kDCAz);

    return flag;
  }

  template <typename T>
  bool IsSelectedV0(T const& v0, const V0PhotonCuts& cut) const
  {
    switch (cut) {
      case V0PhotonCuts::kV0PtRange:
        return v0.pt() >= mMinV0Pt && v0.pt() <= mMaxV0Pt;

      case V0PhotonCuts::kV0EtaRange:
        return v0.eta() >= mMinV0Eta && v0.eta() <= mMaxV0Eta;

      case V0PhotonCuts::kMee:
        return v0.mGamma() <= ((mMaxMeePsiPairDep) ? mMaxMeePsiPairDep(abs(v0.psipair())) : mMaxMee);

      case V0PhotonCuts::kAP:
        return pow(v0.alpha() / mMaxAlpha, 2) + pow(v0.qtarm() / mMaxQt, 2) < 1.0;

      case V0PhotonCuts::kPsiPair:
        return v0.psipair() >= mMinPsiPair && v0.psipair() <= mMaxPsiPair;

      case V0PhotonCuts::kPhiV:
        return v0.phiv() >= mMinPhivPair && v0.phiv() <= mMaxPhivPair;

      case V0PhotonCuts::kRxy:
        return v0.recalculatedVtxR() >= mMinRxy && v0.recalculatedVtxR() <= mMaxRxy;

      case V0PhotonCuts::kCosPA:
        return v0.cospa() >= mMinCosPA;

      case V0PhotonCuts::kPCA:
        return v0.pca() <= mMaxPCA;

      case V0PhotonCuts::kRZLine:
        return v0.recalculatedVtxR() > abs(v0.recalculatedVtxZ()) * TMath::Tan(2 * TMath::ATan(TMath::Exp(-mMaxV0Eta))) - mMaxMarginZ; // as long as z recalculation is not fixed use this

      case V0PhotonCuts::kOnWwireIB: {
        const float margin_xy = 1.0; // cm
        // const float margin_z = 20.0; // cm
        // const float rxy_min = 5.506; // cm
        // const float rxy_max = 14.846;         // cm
        // const float z_min = -17.56; // cm
        // const float z_max = +31.15;           // cm
        float x = abs(v0.recalculatedVtxX()); // cm, measured secondary vertex of gamma->ee
        float y = v0.recalculatedVtxY();      // cm, measured secondary vertex of gamma->ee
        // float z = v0.recalculatedVtxZ();      // cm, measured secondary vertex of gamma->ee

        float rxy = sqrt(x * x + y * y);
        if (rxy < 7.5 || 15.0 < rxy) {
          return false;
        }

        // float z_exp = z_min + (rxy - rxy_min) / TMath::Tan(10.86 * TMath::DegToRad()); // cm, expected position rxy of W wire as a function of z
        // if (abs(z - z_exp) > margin_z) {
        //   return false;
        // }

        float dxy = abs(1.0 * y - x * TMath::Tan(-8.52 * TMath::DegToRad())) / sqrt(pow(1.0, 2) + pow(TMath::Tan(-8.52 * TMath::DegToRad()), 2));
        return !(dxy > margin_xy);
      }
      case V0PhotonCuts::kOnWwireOB: {
        const float margin_x = 2.0;                                         // cm
        const float margin_y = 0.5;                                         // cm
        const float margin_z = 5.0;                                         // cm
        const float rxy_exp = 30.8;                                         // cm
        const float x_exp = rxy_exp * TMath::Cos(-1.3 * TMath::DegToRad()); // cm, expected position x of W wire
        const float y_exp = rxy_exp * TMath::Sin(-1.3 * TMath::DegToRad()); // cm, expected position y of W wire
        const float z_min = -47.0;                                          // cm
        const float z_max = +47.0;                                          // cm
        float x = v0.recalculatedVtxX();                                    // cm, measured secondary vertex of gamma->ee
        float y = v0.recalculatedVtxY();                                    // cm, measured secondary vertex of gamma->ee
        float z = v0.recalculatedVtxZ();                                    // cm, measured secondary vertex of gamma->ee
        if (z + margin_z < z_min || z_max < z - margin_z) {
          return false;
        }
        if (abs(x - x_exp) > margin_x) {
          return false;
        }
        if (abs(y - y_exp) > margin_y) {
          return false;
        }
        return true;
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
        return track.pt() >= mMinTrackPt && track.pt() <= mMaxTrackPt;

      case V0PhotonCuts::kTrackEtaRange:
        return track.eta() >= mMinTrackEta && track.eta() <= mMaxTrackEta;

      case V0PhotonCuts::kTPCNCls:
        return track.tpcNClsFound() >= mMinNClustersTPC;

      case V0PhotonCuts::kTPCCrossedRows:
        return track.tpcNClsCrossedRows() >= mMinNCrossedRowsTPC;

      case V0PhotonCuts::kTPCCrossedRowsOverNCls:
        return track.tpcCrossedRowsOverFindableCls() >= mMinNCrossedRowsOverFindableClustersTPC;

      case V0PhotonCuts::kTPCChi2NDF:
        return mMinChi2PerClusterTPC < track.tpcChi2NCl() && track.tpcChi2NCl() < mMaxChi2PerClusterTPC;

      case V0PhotonCuts::kTPCNsigmaEl:
        return track.tpcNSigmaEl() >= mMinTPCNsigmaEl && track.tpcNSigmaEl() <= mMaxTPCNsigmaEl;

      case V0PhotonCuts::kTPCNsigmaPi:
        return track.tpcNSigmaPi() >= mMinTPCNsigmaPi && track.tpcNSigmaPi() <= mMaxTPCNsigmaPi;

      case V0PhotonCuts::kDCAxy:
        return abs(track.dcaXY()) <= ((mMaxDcaXYPtDep) ? mMaxDcaXYPtDep(track.pt()) : mMaxDcaXY);

      case V0PhotonCuts::kDCAz:
        return abs(track.dcaZ()) <= mMaxDcaZ;

      case V0PhotonCuts::kITSNCls:
        return mMinNClustersITS <= track.itsNCls() && track.itsNCls() <= mMaxNClustersITS;

      case V0PhotonCuts::kITSChi2NDF:
        return mMinChi2PerClusterITS < track.itsChi2NCl() && track.itsChi2NCl() < mMaxChi2PerClusterITS;

      case V0PhotonCuts::kIsWithinBeamPipe: {
        // return track.isWithinBeamPipe();
        if (abs(track.y()) > abs(track.x() * TMath::Tan(10.f * TMath::DegToRad())) + 15.f) {
          return false;
        }

        if (track.x() > 82.9 && abs(track.y()) > abs(track.x() * TMath::Tan(10.f * TMath::DegToRad())) + 5.f) {
          return false;
        }

        const float slope = TMath::Tan(2 * TMath::ATan(TMath::Exp(-0.5)));
        return !(track.x() > 82.9 && abs(track.y()) < 40.f && abs(abs(track.z()) - track.x() / slope) < 3.5f && 15.f < abs(track.dcaXY()));
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

      case V0PhotonCuts::kRequireTPCTRDTOF:
        return isTPCTRDTOFTrack(track);

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
  void SetMaxMarginZ(float max = 7.f);
  void SetMaxMeePsiPairDep(std::function<float(float)> psiDepCut);
  void SetOnWwireIB(bool flag = false);
  void SetOnWwireOB(bool flag = false);

  void SetTrackPtRange(float minPt = 0.f, float maxPt = 1e10f);
  void SetTrackEtaRange(float minEta = -1e10f, float maxEta = 1e10f);
  void SetMinNClustersTPC(int minNClustersTPC);
  void SetMinNCrossedRowsTPC(int minNCrossedRowsTPC);
  void SetMinNCrossedRowsOverFindableClustersTPC(float minNCrossedRowsOverFindableClustersTPC);
  void SetChi2PerClusterTPC(float min, float max);
  void SetNClustersITS(int min, int max);
  void SetChi2PerClusterITS(float min, float max);

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
  void SetRequireTPCTRDTOF(bool flag);

  /// @brief Print the track selection
  void print() const;

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
  float mMaxMarginZ{7.f};
  std::function<float(float)> mMaxMeePsiPairDep{}; // max mee as a function of psipair
  bool mIsOnWwireIB{false};
  bool mIsOnWwireOB{false};

  // pid cuts
  float mMinTPCNsigmaEl{-5}, mMaxTPCNsigmaEl{+5};
  float mMinTPCNsigmaPi{-1e+10}, mMaxTPCNsigmaPi{+1e+10};

  // kinematic cuts
  float mMinTrackPt{0.f}, mMaxTrackPt{1e10f};      // range in pT
  float mMinTrackEta{-1e10f}, mMaxTrackEta{1e10f}; // range in eta

  // track quality cuts
  int mMinNClustersTPC{0};                                           // min number of TPC clusters
  int mMinNCrossedRowsTPC{0};                                        // min number of crossed rows in TPC
  float mMinChi2PerClusterTPC{-1e10f}, mMaxChi2PerClusterTPC{1e10f}; // max tpc fit chi2 per TPC cluster
  float mMinNCrossedRowsOverFindableClustersTPC{0.f};                // min ratio crossed rows / findable clusters
  int mMinNClustersITS{0}, mMaxNClustersITS{7};                      // range in number of ITS clusters
  float mMinChi2PerClusterITS{-1e10f}, mMaxChi2PerClusterITS{1e10f}; // max its fit chi2 per ITS cluster

  float mMaxDcaXY{1e10f};                       // max dca in xy plane
  float mMaxDcaZ{1e10f};                        // max dca in z direction
  std::function<float(float)> mMaxDcaXYPtDep{}; // max dca in xy plane as function of pT
  bool mIsWithinBP{false};
  bool mRequireITSTPC{false};
  bool mRequireITSonly{false};
  bool mRequireTPConly{false};
  bool mRequireTPCTRD{false};
  bool mRequireTPCTOF{false};
  bool mRequireTPCTRDTOF{false};

  ClassDef(V0PhotonCut, 1);
};

#endif // PWGEM_PHOTONMESON_CORE_V0PHOTONCUT_H_
