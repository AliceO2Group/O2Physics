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
    kDCA3Dsigma,
    kDCAxy,
    kDCAz,
    kITSNCls,
    kITSChi2NDF,
    kPrefilter,
    kNCuts
  };
  static const char* mCutNames[static_cast<int>(DalitzEECuts::kNCuts)];

  enum class PIDSchemes : int {
    kUnDef = -1,
    // for nominal B analysis
    kTOFreq = 0,
    kTPChadrej = 1,
    kTPChadrejORTOFreq = 2,
    kTPConly = 3,

    // for low B analysis
    kTOFreq_lowB = 4,
    kTPChadrej_lowB = 5,
    kTPChadrejORTOFreq_lowB = 6,
    kTPConly_lowB = 7,
    kMuon_lowB = 8,
  };

  template <class TLeg, typename TPair>
  bool IsSelected(TPair const& pair) const
  {
    if (!IsSelectedPair(pair)) {
      return false;
    }

    auto pos = pair.template posTrack_as<TLeg>();
    auto ele = pair.template negTrack_as<TLeg>();

    for (auto& track : {pos, ele}) {
      if (!IsSelectedTrack(track)) {
        return false;
      }
    }
    return true;
  }

  bool IsSelectedPair(const float mass, const float phiv) const
  {
    if (mass < mMinMee || mMaxMee < mass) {
      return false;
    }
    if ((phiv < mMinPhivPair || (mMaxPhivPairMeeDep ? mMaxPhivPairMeeDep(mass) : mMaxPhivPair) < phiv) ^ mSelectPC) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool IsSelectedPair(T const& pair) const
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
    return true;
  }

  template <typename T>
  bool IsSelectedTrack(T const& track) const
  {
    if (!track.hasITS() || !track.hasTPC()) { // track has to be ITS-TPC matched track
      return false;
    }

    if (!IsSelectedTrack(track, DalitzEECuts::kTrackPtRange)) {
      return false;
    }
    if (!IsSelectedTrack(track, DalitzEECuts::kTrackEtaRange)) {
      return false;
    }
    if (!IsSelectedTrack(track, DalitzEECuts::kDCA3Dsigma)) {
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

    if (mApplyPF && !IsSelectedTrack(track, DalitzEECuts::kPrefilter)) {
      return false;
    }

    // PID cuts here.
    if (!PassPID(track)) {
      return false;
    }

    if (mApplyTOFbeta && (mMinTOFbeta < track.beta() && track.beta() < mMaxTOFbeta)) {
      return false;
    }

    return true;
  }

  template <typename T>
  bool PassPID(T const& track) const
  {
    switch (mPIDScheme) {
      case PIDSchemes::kTOFreq:
        return PassTOFreq(track);

      case PIDSchemes::kTPChadrej:
        return PassTPChadrej(track);

      case PIDSchemes::kTPChadrejORTOFreq:
        return PassTPChadrej(track) || PassTOFreq(track);

      case PIDSchemes::kTPConly:
        return PassTPConly(track);

      case PIDSchemes::kTOFreq_lowB:
        return PassTOFreq(track);

      case PIDSchemes::kTPChadrej_lowB:
        return PassTPChadrej_lowB(track);

      case PIDSchemes::kTPChadrejORTOFreq_lowB:
        return PassTPChadrej_lowB(track) || PassTOFreq_lowB(track);

      case PIDSchemes::kTPConly_lowB:
        return PassTPConly_lowB(track);

      case PIDSchemes::kMuon_lowB:
        return PassMuon_lowB(track);

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
  bool PassTOFreq_lowB(T const& track) const
  {
    bool is_el_included_TPC = mMinTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < mMaxTPCNsigmaEl;
    bool is_pi_excluded_TPC = track.tpcInnerParam() < 0.4 ? true : (track.tpcNSigmaPi() < mMinTPCNsigmaPi || mMaxTPCNsigmaPi < track.tpcNSigmaPi());
    bool is_el_included_TOF = mMinTOFNsigmaEl < track.tofNSigmaEl() && track.tofNSigmaEl() < mMaxTOFNsigmaEl;
    return is_el_included_TPC && is_pi_excluded_TPC && is_el_included_TOF;
  }

  template <typename T>
  bool PassTPChadrej_lowB(T const& track) const
  {
    bool is_el_included_TPC = mMinTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < mMaxTPCNsigmaEl;
    bool is_mu_excluded_TPC = mMuonExclusionTPC ? track.tpcNSigmaMu() < mMinTPCNsigmaMu || mMaxTPCNsigmaMu < track.tpcNSigmaMu() : true;
    bool is_pi_excluded_TPC = track.tpcInnerParam() < 0.4 ? true : (track.tpcNSigmaPi() < mMinTPCNsigmaPi || mMaxTPCNsigmaPi < track.tpcNSigmaPi());
    bool is_ka_excluded_TPC = track.tpcNSigmaKa() < mMinTPCNsigmaKa || mMaxTPCNsigmaKa < track.tpcNSigmaKa();
    bool is_pr_excluded_TPC = track.tpcNSigmaPr() < mMinTPCNsigmaPr || mMaxTPCNsigmaPr < track.tpcNSigmaPr();
    return is_el_included_TPC && is_mu_excluded_TPC && is_pi_excluded_TPC && is_ka_excluded_TPC && is_pr_excluded_TPC;
  }

  template <typename T>
  bool PassTPConly_lowB(T const& track) const
  {
    bool is_el_included_TPC = mMinTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < mMaxTPCNsigmaEl;
    bool is_pi_excluded_TPC = track.tpcInnerParam() < 0.4 ? true : (track.tpcNSigmaPi() < mMinTPCNsigmaPi || mMaxTPCNsigmaPi < track.tpcNSigmaPi());
    return is_el_included_TPC && is_pi_excluded_TPC;
  }

  template <typename T>
  bool PassMuon_lowB(T const& track) const
  {
    bool is_el_excluded_TPC = track.tpcNSigmaEl() < mMinTPCNsigmaEl || mMaxTPCNsigmaEl < track.tpcNSigmaEl();
    if (!is_el_excluded_TPC) {
      return false;
    }
    if (track.hasTOF()) {
      bool is_mu_included_TPC = mMinTPCNsigmaMu < track.tpcNSigmaMu() && track.tpcNSigmaMu() < mMaxTPCNsigmaMu;
      bool is_mu_included_TOF = mMinTOFNsigmaMu < track.tofNSigmaMu() && track.tofNSigmaMu() < mMaxTOFNsigmaMu;
      bool is_pi_excluded_TOF = track.tofNSigmaPi() < mMinTOFNsigmaPi;
      return is_mu_included_TPC && is_mu_included_TOF && is_pi_excluded_TOF;
    } else if (track.tpcInnerParam() < mMaxPinMuonTPConly) {
      bool is_mu_included_TPC = mMinTPCNsigmaMu < track.tpcNSigmaMu() && track.tpcNSigmaMu() < mMaxTPCNsigmaMu;
      bool is_pi_excluded_TPC = track.tpcNSigmaPi() < mMinTPCNsigmaPi;
      return is_mu_included_TPC && is_pi_excluded_TPC;
    } else {
      return false;
    }
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
        return mMinMee <= pair.mass() && pair.mass() <= mMaxMee;

      case DalitzEECuts::kPhiV:
        return (mMinPhivPair <= pair.phiv() && pair.phiv() <= (mMaxPhivPairMeeDep ? mMaxPhivPairMeeDep(pair.mass()) : mMaxPhivPair)) ^ mSelectPC;

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

      case DalitzEECuts::kDCA3Dsigma: {
        float dca_3d = 999.f;
        float det = track.cYY() * track.cZZ() - track.cZY() * track.cZY();
        if (det < 0) {
          dca_3d = 999.f;
        } else {
          float chi2 = (track.dcaXY() * track.dcaXY() * track.cZZ() + track.dcaZ() * track.dcaZ() * track.cYY() - 2. * track.dcaXY() * track.dcaZ() * track.cZY()) / det;
          dca_3d = std::sqrt(std::abs(chi2) / 2.);
        }
        return mMinDca3D <= dca_3d && dca_3d <= mMaxDca3D; // in sigma for single leg
      }
      case DalitzEECuts::kDCAxy:
        return abs(track.dcaXY()) <= ((mMaxDcaXYPtDep) ? mMaxDcaXYPtDep(track.pt()) : mMaxDcaXY);

      case DalitzEECuts::kDCAz:
        return abs(track.dcaZ()) <= mMaxDcaZ;

      case DalitzEECuts::kITSNCls:
        return mMinNClustersITS <= track.itsNCls() && track.itsNCls() <= mMaxNClustersITS;

      case DalitzEECuts::kITSChi2NDF:
        return mMinChi2PerClusterITS < track.itsChi2NCl() && track.itsChi2NCl() < mMaxChi2PerClusterITS;

      case DalitzEECuts::kPrefilter:
        return track.pfb() <= 0;

      default:
        return false;
    }
  }

  // Setters
  void SetPairPtRange(float minPt = 0.f, float maxPt = 1e10f);
  void SetPairEtaRange(float minEta = -1e10f, float maxEta = 1e10f);
  void SetMeeRange(float min = 0.f, float max = 0.5);
  void SetMaxPhivPairMeeDep(std::function<float(float)> meeDepCut);
  void SelectPhotonConversion(bool flag);

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
  void SetMaxPinMuonTPConly(float max);

  void SetDca3DRange(float min, float max); // in sigma
  void SetMaxDcaXY(float maxDcaXY);         // in cm
  void SetMaxDcaZ(float maxDcaZ);           // in cm
  void SetMaxDcaXYPtDep(std::function<float(float)> ptDepCut);
  void ApplyPrefilter(bool flag);

  // Getters
  bool IsPhotonConversionSelected() const { return mSelectPC; }

  /// @brief Print the track selection
  void print() const;

 private:
  // pair cuts
  float mMinMee{0.f}, mMaxMee{1e10f};
  float mMinPairPt{0.f}, mMaxPairPt{1e10f};      // range in pT
  float mMinPairEta{-1e10f}, mMaxPairEta{1e10f}; // range in eta
  float mMinPhivPair{0.f}, mMaxPhivPair{+3.2};
  std::function<float(float)> mMaxPhivPairMeeDep{}; // max phiv as a function of mee
  bool mSelectPC{false};                            // flag to select photon conversion used in mMaxPhivPairMeeDep

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
  float mMaxPinMuonTPConly{0.2f};                                    // max pin cut for muon ID with TPConly

  float mMinDca3D{0.0f};                        // min dca in 3D in units of sigma
  float mMaxDca3D{1e+10};                       // max dca in 3D in units of sigma
  float mMaxDcaXY{1.0f};                        // max dca in xy plane
  float mMaxDcaZ{1.0f};                         // max dca in z direction
  std::function<float(float)> mMaxDcaXYPtDep{}; // max dca in xy plane as function of pT
  bool mApplyPF{false};

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
