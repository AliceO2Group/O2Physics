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
// Class for dalitz ee selection
//

#ifndef PWGEM_PHOTONMESON_CORE_DALITZEECUT_H_
#define PWGEM_PHOTONMESON_CORE_DALITZEECUT_H_

#include <algorithm>
#include <set>
#include <vector>
#include <utility>
#include <string>
#include "Framework/Logger.h"
#include "Framework/DataTypes.h"
#include "Rtypes.h"
#include "TNamed.h"
#include "TMath.h"

class DalitzEECut : public TNamed
{
 public:
  DalitzEECut() = default;
  DalitzEECut(const char* name, const char* title) : TNamed(name, title) {}

  enum class DalitzEECuts : int {
    // pair cut
    kMee = 0,
    kPairPtRange,
    kPairEtaRange,
    kPhiV,
    // track cut
    kTrackPtRange,
    kTrackEtaRange,
    kTPCNCls,
    kTPCCrossedRows,
    kTPCCrossedRowsOverNCls,
    kTPCChi2NDF,
    kTPCNsigmaEl,
    kTPCNsigmaMu,
    kTPCNsigmaPi,
    kTPCNsigmaKa,
    kTPCNsigmaPr,
    kTOFNsigmaEl,
    kTOFNsigmaMu,
    kTOFNsigmaPi,
    kTOFNsigmaKa,
    kTOFNsigmaPr,
    kDCAxy,
    kDCAz,
    kITSNCls,
    kITSChi2NDF,
    kNCuts
  };
  static const char* mCutNames[static_cast<int>(DalitzEECuts::kNCuts)];

  enum class PIDSchemes : int {
    kUnDef = -1,
    // for nominal B analysis
    kTOFreq = 0,
    kTOFif = 1,
    kTPChadrej = 2,
    kTPChadrejORTOFreq = 3,
    kTPConly = 4,
  };

  template <class TLeg, typename TPair>
  bool IsSelected(TPair const& pair) const
  {
    if (!IsSelectedPair(pair, DalitzEECuts::kPairPtRange)) {
      return false;
    }
    if (!IsSelectedPair(pair, DalitzEECuts::kPairEtaRange)) {
      return false;
    }
    if (!IsSelectedPair(pair, DalitzEECuts::kMee)) {
      return false;
    }
    if (!IsSelectedPair(pair, DalitzEECuts::kPhiV)) {
      return false;
    }

    auto pos = pair.template posTrack_as<TLeg>();
    auto ele = pair.template negTrack_as<TLeg>();

    for (auto& track : {pos, ele}) {
      if (!track.hasITS() || !track.hasTPC()) { // track has to be ITS-TPC matched track
        return false;
      }

      if (!IsSelectedTrack(track, DalitzEECuts::kTrackPtRange)) {
        return false;
      }
      if (!IsSelectedTrack(track, DalitzEECuts::kTrackEtaRange)) {
        return false;
      }
      if (!IsSelectedTrack(track, DalitzEECuts::kDCAxy)) {
        return false;
      }
      if (!IsSelectedTrack(track, DalitzEECuts::kDCAz)) {
        return false;
      }

      // ITS cuts
      if (!IsSelectedTrack(track, DalitzEECuts::kITSNCls)) {
        return false;
      }
      if (!IsSelectedTrack(track, DalitzEECuts::kITSChi2NDF)) {
        return false;
      }

      // TPC cuts
      if (!IsSelectedTrack(track, DalitzEECuts::kTPCNCls)) {
        return false;
      }
      if (!IsSelectedTrack(track, DalitzEECuts::kTPCCrossedRows)) {
        return false;
      }
      if (!IsSelectedTrack(track, DalitzEECuts::kTPCCrossedRowsOverNCls)) {
        return false;
      }
      if (!IsSelectedTrack(track, DalitzEECuts::kTPCChi2NDF)) {
        return false;
      }

      // PID cuts here.
      if (!PassPID(track)) {
        return false;
      }

      if (mApplyTOFbeta && (mMinTOFbeta < track.beta() && track.beta() < mMaxTOFbeta)) {
        return false;
      }
    }
    return true;
  }

  template <typename T>
  bool PassPID(T const& track) const
  {
    switch (mPIDScheme) {
      case PIDSchemes::kTOFreq:
        return PassTOFreq(track);

      case PIDSchemes::kTOFif:
        return PassTOFif(track);

      case PIDSchemes::kTPChadrej:
        return PassTPChadrej(track);

      case PIDSchemes::kTPChadrejORTOFreq:
        return PassTPChadrej(track) || PassTOFreq(track);

      case PIDSchemes::kTPConly:
        return PassTPConly(track);

      case PIDSchemes::kUnDef:
        return true;

      default:
        return true;
    }
  }

  template <typename T>
  bool PassTOFreq(T const& track) const
  {
    bool is_el_included_TPC = mMinTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < mMaxTPCNsigmaEl;
    bool is_pi_excluded_TPC = track.tpcNSigmaPi() < mMinTPCNsigmaPi || mMaxTPCNsigmaPi < track.tpcNSigmaPi();
    bool is_el_included_TOF = mMinTOFNsigmaEl < track.tofNSigmaEl() && track.tofNSigmaEl() < mMaxTOFNsigmaEl;
    return is_el_included_TPC && is_pi_excluded_TPC && is_el_included_TOF;
  }

  template <typename T>
  bool PassTOFif(T const& track) const
  {
    bool is_el_included_TPC = mMinTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < mMaxTPCNsigmaEl;
    bool is_pi_excluded_TPC = track.tpcNSigmaPi() < mMinTPCNsigmaPi || mMaxTPCNsigmaPi < track.tpcNSigmaPi();
    bool is_el_included_TOF = (track.tpcInnerParam() < mMinPinTOF || track.beta() < 0.0) ? true : mMinTOFNsigmaEl < track.tofNSigmaEl() && track.tofNSigmaEl() < mMaxTOFNsigmaEl;
    return is_el_included_TPC && is_pi_excluded_TPC && is_el_included_TOF;
  }

  template <typename T>
  bool PassTPChadrej(T const& track) const
  {
    bool is_el_included_TPC = mMinTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < mMaxTPCNsigmaEl;
    bool is_mu_excluded_TPC = mMuonExclusionTPC ? track.tpcNSigmaMu() < mMinTPCNsigmaMu || mMaxTPCNsigmaMu < track.tpcNSigmaMu() : true;
    bool is_pi_excluded_TPC = track.tpcNSigmaPi() < mMinTPCNsigmaPi || mMaxTPCNsigmaPi < track.tpcNSigmaPi();
    bool is_ka_excluded_TPC = track.tpcNSigmaKa() < mMinTPCNsigmaKa || mMaxTPCNsigmaKa < track.tpcNSigmaKa();
    bool is_pr_excluded_TPC = track.tpcNSigmaPr() < mMinTPCNsigmaPr || mMaxTPCNsigmaPr < track.tpcNSigmaPr();
    return is_el_included_TPC && is_mu_excluded_TPC && is_pi_excluded_TPC && is_ka_excluded_TPC && is_pr_excluded_TPC;
  }

  template <typename T>
  bool PassTPConly(T const& track) const
  {
    bool is_el_included_TPC = mMinTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < mMaxTPCNsigmaEl;
    bool is_pi_excluded_TPC = track.tpcNSigmaPi() < mMinTPCNsigmaPi || mMaxTPCNsigmaPi < track.tpcNSigmaPi();
    return is_el_included_TPC && is_pi_excluded_TPC;
  }

  template <typename T>
  bool IsSelectedPair(T const& pair, const DalitzEECuts& cut) const
  {
    switch (cut) {
      case DalitzEECuts::kPairPtRange:
        return pair.pt() >= mMinPairPt && pair.pt() <= mMaxPairPt;

      case DalitzEECuts::kPairEtaRange:
        return pair.eta() >= mMinPairEta && pair.eta() <= mMaxPairEta;

      case DalitzEECuts::kMee:
        return mMinMee <= pair.mee() && pair.mee() <= mMaxMee;

      case DalitzEECuts::kPhiV:
        return mMinPhivPair <= pair.phiv() && pair.phiv() <= (mMaxPhivPairMeeDep ? mMaxPhivPairMeeDep(pair.mee()) : mMaxPhivPair);

      default:
        return false;
    }
  }

  template <typename T>
  bool IsSelectedTrack(T const& track, const DalitzEECuts& cut) const
  {
    switch (cut) {
      case DalitzEECuts::kTrackPtRange:
        return track.pt() >= mMinTrackPt && track.pt() <= mMaxTrackPt;

      case DalitzEECuts::kTrackEtaRange:
        return track.eta() >= mMinTrackEta && track.eta() <= mMaxTrackEta;

      case DalitzEECuts::kTPCNCls:
        return track.tpcNClsFound() >= mMinNClustersTPC;

      case DalitzEECuts::kTPCCrossedRows:
        return track.tpcNClsCrossedRows() >= mMinNCrossedRowsTPC;

      case DalitzEECuts::kTPCCrossedRowsOverNCls:
        return track.tpcCrossedRowsOverFindableCls() >= mMinNCrossedRowsOverFindableClustersTPC;

      case DalitzEECuts::kTPCChi2NDF:
        return mMinChi2PerClusterTPC < track.tpcChi2NCl() && track.tpcChi2NCl() < mMaxChi2PerClusterTPC;

      case DalitzEECuts::kDCAxy:
        return abs(track.dcaXY()) <= ((mMaxDcaXYPtDep) ? mMaxDcaXYPtDep(track.pt()) : mMaxDcaXY);

      case DalitzEECuts::kDCAz:
        return abs(track.dcaZ()) <= mMaxDcaZ;

      case DalitzEECuts::kITSNCls:
        return mMinNClustersITS <= track.itsNCls() && track.itsNCls() <= mMaxNClustersITS;

      case DalitzEECuts::kITSChi2NDF:
        return mMinChi2PerClusterITS < track.itsChi2NCl() && track.itsChi2NCl() < mMaxChi2PerClusterITS;

      default:
        return false;
    }
  }

  // Setters
  void SetPairPtRange(float minPt = 0.f, float maxPt = 1e10f);
  void SetPairEtaRange(float minEta = -1e10f, float maxEta = 1e10f);
  void SetMeeRange(float min = 0.f, float max = 0.5);
  void SetMaxPhivPairMeeDep(std::function<float(float)> meeDepCut);

  void SetTrackPtRange(float minPt = 0.f, float maxPt = 1e10f);
  void SetTrackEtaRange(float minEta = -1e10f, float maxEta = 1e10f);
  void SetMinNClustersTPC(int minNClustersTPC);
  void SetMinNCrossedRowsTPC(int minNCrossedRowsTPC);
  void SetMinNCrossedRowsOverFindableClustersTPC(float minNCrossedRowsOverFindableClustersTPC);
  void SetChi2PerClusterTPC(float min, float max);
  void SetNClustersITS(int min, int max);
  void SetChi2PerClusterITS(float min, float max);

  void SetPIDScheme(PIDSchemes scheme);
  void SetMinPinTOF(float min);
  void SetMuonExclusionTPC(bool flag);
  void SetTOFbetaRange(bool flag, float min, float max);
  void SetTPCNsigmaElRange(float min = -1e+10, float max = 1e+10);
  void SetTPCNsigmaMuRange(float min = -1e+10, float max = 1e+10);
  void SetTPCNsigmaPiRange(float min = -1e+10, float max = 1e+10);
  void SetTPCNsigmaKaRange(float min = -1e+10, float max = 1e+10);
  void SetTPCNsigmaPrRange(float min = -1e+10, float max = 1e+10);
  void SetTOFNsigmaElRange(float min = -1e+10, float max = 1e+10);
  void SetTOFNsigmaMuRange(float min = -1e+10, float max = 1e+10);
  void SetTOFNsigmaPiRange(float min = -1e+10, float max = 1e+10);
  void SetTOFNsigmaKaRange(float min = -1e+10, float max = 1e+10);
  void SetTOFNsigmaPrRange(float min = -1e+10, float max = 1e+10);

  void SetMaxDcaXY(float maxDcaXY);
  void SetMaxDcaZ(float maxDcaZ);
  void SetMaxDcaXYPtDep(std::function<float(float)> ptDepCut);

  /// @brief Print the track selection
  void print() const;

 private:
  // pair cuts
  float mMinMee{0.f}, mMaxMee{1e10f};
  float mMinPairPt{0.f}, mMaxPairPt{1e10f};      // range in pT
  float mMinPairEta{-1e10f}, mMaxPairEta{1e10f}; // range in eta
  float mMinPhivPair{0.f}, mMaxPhivPair{+3.2};
  std::function<float(float)> mMaxPhivPairMeeDep{}; // max phiv as a function of mee

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

  float mMaxDcaXY{1.0f};                        // max dca in xy plane
  float mMaxDcaZ{1.0f};                         // max dca in z direction
  std::function<float(float)> mMaxDcaXYPtDep{}; // max dca in xy plane as function of pT

  // pid cuts
  PIDSchemes mPIDScheme{PIDSchemes::kUnDef};
  float mMinPinTOF{0.0f};        // min pin cut for TOF.
  bool mMuonExclusionTPC{false}; // flag to reject muon in TPC for low B
  bool mApplyTOFbeta{false};     // flag to reject hadron contamination with TOF
  float mMinTOFbeta{0.0}, mMaxTOFbeta{0.95};
  float mMinTPCNsigmaEl{-1e+10}, mMaxTPCNsigmaEl{+1e+10};
  float mMinTPCNsigmaMu{-1e+10}, mMaxTPCNsigmaMu{+1e+10};
  float mMinTPCNsigmaPi{-1e+10}, mMaxTPCNsigmaPi{+1e+10};
  float mMinTPCNsigmaKa{-1e+10}, mMaxTPCNsigmaKa{+1e+10};
  float mMinTPCNsigmaPr{-1e+10}, mMaxTPCNsigmaPr{+1e+10};

  float mMinTOFNsigmaEl{-1e+10}, mMaxTOFNsigmaEl{+1e+10};
  float mMinTOFNsigmaMu{-1e+10}, mMaxTOFNsigmaMu{+1e+10};
  float mMinTOFNsigmaPi{-1e+10}, mMaxTOFNsigmaPi{+1e+10};
  float mMinTOFNsigmaKa{-1e+10}, mMaxTOFNsigmaKa{+1e+10};
  float mMinTOFNsigmaPr{-1e+10}, mMaxTOFNsigmaPr{+1e+10};

  ClassDef(DalitzEECut, 1);
};

#endif // PWGEM_PHOTONMESON_CORE_DALITZEECUT_H_
