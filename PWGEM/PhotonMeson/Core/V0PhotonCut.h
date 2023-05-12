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

#include <set>
#include <vector>
#include <utility>
#include <string>
#include "Framework/Logger.h"
#include "Framework/DataTypes.h"
#include "Rtypes.h"
#include "TNamed.h"
#include "TMath.h"

class V0PhotonCut : public TNamed
{
 public:
  V0PhotonCut() = default;
  V0PhotonCut(const char* name, const char* title) : TNamed(name, title) {}

  enum class V0PhotonCuts : int {
    // v0 cut
    kMee = 0,
    kPsiPair,
    kRxy,
    kCosPA,
    kPCA,
    kRZLine,
    kOnWwireIB,
    kOnWwireOB,
    // leg cut
    kPtRange,
    kEtaRange,
    kTPCNCls,
    kTPCCrossedRows,
    kTPCCrossedRowsOverNCls,
    kTPCChi2NDF,
    kTPCNsigmaEl,
    kTPCNsigmaPi,
    kDCAxy,
    kDCAz,
    kNCuts
  };

  static const char* mCutNames[static_cast<int>(V0PhotonCuts::kNCuts)];

  // Temporary function to check if track passes selection criteria. To be replaced by framework filters.
  template <class TLeg, typename TV0>
  bool IsSelected(TV0 const& v0) const
  {
    if (!IsSelectedV0(v0, V0PhotonCuts::kMee)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kPsiPair)) {
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
    if (mIsOnWwireIB && !IsSelectedV0(v0, V0PhotonCuts::kOnWwireIB)) {
      return false;
    }
    if (mIsOnWwireOB && !IsSelectedV0(v0, V0PhotonCuts::kOnWwireOB)) {
      return false;
    }

    auto pos = v0.template posTrack_as<TLeg>();
    auto ele = v0.template negTrack_as<TLeg>();
    for (auto& track : {pos, ele}) {
      if (!IsSelectedTrack(track, V0PhotonCuts::kPtRange)) {
        return false;
      }
      if (!IsSelectedTrack(track, V0PhotonCuts::kEtaRange)) {
        return false;
      }
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
      if (!IsSelectedTrack(track, V0PhotonCuts::kDCAxy)) {
        return false;
      }
      if (!IsSelectedTrack(track, V0PhotonCuts::kDCAz)) {
        return false;
      }
    }
    return true;
  }

  // Temporary function to check if track passes and return a flag. To be replaced by framework filters.
  template <typename T>
  uint32_t IsSelectedMask(T const& track) const
  {
    uint32_t flag = 0;

    auto setFlag = [&](const V0PhotonCuts& cut) {
      if (IsSelectedTrack(track, cut)) {
        flag |= 1UL << static_cast<int>(cut);
      }
    };

    setFlag(V0PhotonCuts::kPtRange);
    setFlag(V0PhotonCuts::kEtaRange);
    setFlag(V0PhotonCuts::kTPCNCls);
    setFlag(V0PhotonCuts::kTPCCrossedRows);
    setFlag(V0PhotonCuts::kTPCCrossedRowsOverNCls);
    setFlag(V0PhotonCuts::kTPCChi2NDF);
    setFlag(V0PhotonCuts::kDCAxy);
    setFlag(V0PhotonCuts::kDCAz);

    return flag;
  }

  // Temporary function to check if track passes a given selection criteria. To be replaced by framework filters.
  template <typename T>
  bool IsSelectedV0(T const& v0, const V0PhotonCuts& cut) const
  {
    const float margin = 1.0; // cm
    switch (cut) {
      case V0PhotonCuts::kMee:
        return v0.mGamma() <= ((mMaxMeePsiPairDep) ? mMaxMeePsiPairDep(abs(v0.psipair())) : mMaxMee);

      case V0PhotonCuts::kPsiPair:
        return v0.psipair() >= mMinPsiPair && v0.psipair() <= mMaxPsiPair;

      case V0PhotonCuts::kRxy:
        return v0.recalculatedVtxR() >= mMinRxy && v0.recalculatedVtxR() <= mMaxRxy;

      case V0PhotonCuts::kCosPA:
        return v0.cospa() >= mMinCosPA;

      case V0PhotonCuts::kPCA:
        return v0.pca() <= mMaxPCA;

      case V0PhotonCuts::kRZLine:
        return TMath::Abs(v0.vz()) < 12.0 + v0.recalculatedVtxR() * TMath::Tan(2 * TMath::ATan(TMath::Exp(-mMaxEta))); // as long as z recalculation is not fixed use this

      case V0PhotonCuts::kOnWwireIB: {
        const float rxy_min = 5.506;          // cm
        const float rxy_max = 14.846;         // cm
        const float z_min = -17.56;           // cm
        const float z_max = +31.15;           // cm
        float x = abs(v0.recalculatedVtxX()); // cm, measured secondary vertex of gamma->ee
        float y = v0.recalculatedVtxY();      // cm, measured secondary vertex of gamma->ee
        float z = v0.recalculatedVtxZ();      // cm, measured secondary vertex of gamma->ee
        float rxy = sqrt(x * x + y * y);
        if ((rxy < rxy_min || rxy_max < rxy) || (z < z_min || z_max < z)) {
          return false;
        }
        float x_exp = abs(rxy * TMath::Cos(-8.52 * TMath::DegToRad()));                // cm, expected position x of W wire
        float y_exp = rxy * TMath::Sin(-8.52 * TMath::DegToRad());                     // cm, expected position y of W wire
        float z_exp = z_min + (rxy - rxy_min) / TMath::Tan(10.86 * TMath::DegToRad()); // cm, expected position rxy of W wire as a function of z
        if (abs(x - x_exp) > margin || abs(y - y_exp) > margin || abs(z - z_exp) > margin) {
          return false;
        }
        return true;
      }
      case V0PhotonCuts::kOnWwireOB: {
        const float rxy_min = 30.8 - margin;  // cm
        const float rxy_max = 30.8 + margin;  // cm
        const float z_min = -47.0;            // cm
        const float z_max = +47.0;            // cm
        float x = abs(v0.recalculatedVtxX()); // cm, measured secondary vertex of gamma->ee
        float y = v0.recalculatedVtxY();      // cm, measured secondary vertex of gamma->ee
        float z = v0.recalculatedVtxZ();      // cm, measured secondary vertex of gamma->ee
        float rxy = sqrt(x * x + y * y);
        if ((rxy < rxy_min || rxy_max < rxy) || (z < z_min || z_max < z)) {
          return false;
        }
        float x_exp = abs(rxy * TMath::Cos(-1.3 * TMath::DegToRad())); // cm, expected position x of W wire
        float y_exp = rxy * TMath::Sin(-1.3 * TMath::DegToRad());      // cm, expected position y of W wire
        if (abs(x - x_exp) > margin || abs(y - y_exp) > margin) {
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
      case V0PhotonCuts::kPtRange:
        return track.pt() >= mMinPt && track.pt() <= mMaxPt;

      case V0PhotonCuts::kEtaRange:
        return track.eta() >= mMinEta && track.eta() <= mMaxEta;

      case V0PhotonCuts::kTPCNCls:
        return track.tpcNClsFound() >= mMinNClustersTPC;

      case V0PhotonCuts::kTPCCrossedRows:
        return track.tpcNClsCrossedRows() >= mMinNCrossedRowsTPC;

      case V0PhotonCuts::kTPCCrossedRowsOverNCls:
        return track.tpcCrossedRowsOverFindableCls() >= mMinNCrossedRowsOverFindableClustersTPC;

      case V0PhotonCuts::kTPCChi2NDF:
        return track.tpcChi2NCl() <= mMaxChi2PerClusterTPC;

      case V0PhotonCuts::kTPCNsigmaEl:
        return track.tpcNSigmaEl() >= mMinTPCNsigmaEl && track.tpcNSigmaEl() <= mMaxTPCNsigmaEl;

      case V0PhotonCuts::kTPCNsigmaPi:
        return track.tpcNSigmaPi() >= mMinTPCNsigmaPi && track.tpcNSigmaPi() <= mMaxTPCNsigmaPi;

      case V0PhotonCuts::kDCAxy:
        return abs(track.dcaXY()) <= ((mMaxDcaXYPtDep) ? mMaxDcaXYPtDep(track.pt()) : mMaxDcaXY);

      case V0PhotonCuts::kDCAz:
        return abs(track.dcaZ()) <= mMaxDcaZ;

      default:
        return false;
    }
  }

  // Setters
  void SetMeeRange(float min = 0.f, float max = 0.1);
  void SetPsiPairRange(float min = -3.15, float max = +3.15);
  void SetRxyRange(float min = 0.f, float max = 180.f);
  void SetMinCosPA(float min = 0.95);
  void SetMaxPCA(float max = 2.f);
  void SetMaxMeePsiPairDep(std::function<float(float)> psiDepCut);
  void SetOnWwireIB(bool flag = false);
  void SetOnWwireOB(bool flag = false);

  void SetPtRange(float minPt = 0.f, float maxPt = 1e10f);
  void SetEtaRange(float minEta = -1e10f, float maxEta = 1e10f);
  void SetMinNClustersTPC(int minNClustersTPC);
  void SetMinNCrossedRowsTPC(int minNCrossedRowsTPC);
  void SetMinNCrossedRowsOverFindableClustersTPC(float minNCrossedRowsOverFindableClustersTPC);
  void SetMaxChi2PerClusterTPC(float maxChi2PerClusterTPC);

  void SetTPCNsigmaElRange(float min = -3, float max = +3);
  void SetTPCNsigmaPiRange(float min = -1e+10, float max = 1e+10);

  void SetMaxDcaXY(float maxDcaXY);
  void SetMaxDcaZ(float maxDcaZ);
  void SetMaxDcaXYPtDep(std::function<float(float)> ptDepCut);

  /// @brief Print the track selection
  void print() const;

 private:
  // v0 cuts
  float mMinMee{0.f}, mMaxMee{0.1f};
  float mMinPsiPair{-3.15}, mMaxPsiPair{+3.15};
  float mMinRxy{0.f}, mMaxRxy{180.f};
  float mMinCosPA{0.95};
  float mMaxPCA{2.f};
  std::function<float(float)> mMaxMeePsiPairDep{}; // max mee as a function of psipair
  bool mIsOnWwireIB{false};
  bool mIsOnWwireOB{false};

  // pid cuts
  float mMinTPCNsigmaEl{-5}, mMaxTPCNsigmaEl{+5};
  float mMinTPCNsigmaPi{-1e+10}, mMaxTPCNsigmaPi{+1e+10};

  // kinematic cuts
  float mMinPt{0.f}, mMaxPt{1e10f};      // range in pT
  float mMinEta{-1e10f}, mMaxEta{1e10f}; // range in eta

  // track quality cuts
  int mMinNClustersTPC{0};                            // min number of TPC clusters
  int mMinNCrossedRowsTPC{0};                         // min number of crossed rows in TPC
  float mMaxChi2PerClusterTPC{1e10f};                 // max tpc fit chi2 per TPC cluster
  float mMinNCrossedRowsOverFindableClustersTPC{0.f}; // min ratio crossed rows / findable clusters

  float mMaxDcaXY{1e10f};                       // max dca in xy plane
  float mMaxDcaZ{1e10f};                        // max dca in z direction
  std::function<float(float)> mMaxDcaXYPtDep{}; // max dca in xy plane as function of pT

  ClassDef(V0PhotonCut, 1);
};

#endif // PWGEM_PHOTONMESON_CORE_V0PHOTONCUT_H_
