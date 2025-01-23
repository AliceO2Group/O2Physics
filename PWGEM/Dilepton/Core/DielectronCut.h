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
// Class for dielectron selection
//

#ifndef PWGEM_DILEPTON_CORE_DIELECTRONCUT_H_
#define PWGEM_DILEPTON_CORE_DIELECTRONCUT_H_

#include <algorithm>
#include <set>
#include <vector>
#include <utility>
#include <string>
#include "TNamed.h"
#include "Math/Vector4D.h"

#include "PWGEM/Dilepton/Utils/MlResponseDielectronSingleTrack.h"

#include "Framework/Logger.h"
#include "Framework/DataTypes.h"
#include "CommonConstants/PhysicsConstants.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"
#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"

using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;
using namespace o2::aod::pwgem::dilepton::utils::pairutil;

class DielectronCut : public TNamed
{
 public:
  DielectronCut() = default;
  DielectronCut(const char* name, const char* title) : TNamed(name, title) {}
  ~DielectronCut() {}

  enum class DielectronCuts : int {
    // pair cut
    kMee = 0,
    kPairPtRange,
    kPairYRange,
    kPairDCARange,
    kPhiV,
    // track cut
    kTrackPtRange,
    kTrackEtaRange,
    kTrackPhiRange,
    kTPCNCls,
    kTPCCrossedRows,
    kTPCCrossedRowsOverNCls,
    kTPCFracSharedClusters,
    kRelDiffPin,
    kTPCChi2NDF,
    kDCA3Dsigma,
    kDCAxy,
    kDCAz,
    kITSNCls,
    kITSChi2NDF,
    kITSClusterSize,
    kPrefilter,
    kNCuts
  };

  enum class PIDSchemes : int {
    kUnDef = -1,
    kTOFreq = 0,
    kTPChadrej = 1,
    kTPChadrejORTOFreq = 2,
    kTPConly = 3,
    kTOFif = 4,
    kPIDML = 5,
    kTPChadrejORTOFreq_woTOFif = 6
  };

  template <typename T = int, typename TPair>
  bool IsSelected(TPair const& pair) const
  {
    auto t1 = std::get<0>(pair);
    auto t2 = std::get<1>(pair);
    float bz = std::get<2>(pair);

    if (!IsSelectedTrack(t1) || !IsSelectedTrack(t2)) {
      return false;
    }

    if (!IsSelectedPair(t1, t2, bz)) {
      return false;
    }

    return true;
  }

  template <bool dont_require_rapidity = false, typename TTrack1, typename TTrack2>
  bool IsSelectedPair(TTrack1 const& t1, TTrack2 const& t2, const float bz) const
  {
    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    float dca_ee_3d = pairDCAQuadSum(dca3DinSigma(t1), dca3DinSigma(t2));
    float phiv = getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), bz);
    float opAng = getOpeningAngle(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz());

    if (v12.M() < mMinMee || mMaxMee < v12.M()) {
      return false;
    }

    if (!dont_require_rapidity && (v12.Rapidity() < mMinPairY || mMaxPairY < v12.Rapidity())) {
      return false;
    }

    if (mApplyPhiV) {
      if (((mMinPhivPair < phiv && phiv < mMaxPhivPair) && v12.M() < mMaxMeePhiVDep(phiv)) ^ mSelectPC) {
        return false;
      }
    }

    if (dca_ee_3d < mMinPairDCA3D || mMaxPairDCA3D < dca_ee_3d) { // in sigma for pair
      return false;
    }

    if (opAng < mMinOpAng || mMaxOpAng < opAng) {
      return false;
    }

    if (mRequireDiffSides && t1.eta() * t2.eta() > 0.0) {
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

  template <bool dont_require_pteta = false, bool isML = false, typename TTrack, typename TCollision = int>
  bool IsSelectedTrack(TTrack const& track, TCollision const& collision = 0) const
  {
    if (!track.hasITS() || !track.hasTPC()) { // track has to be ITS-TPC matched track
      return false;
    }

    if (!dont_require_pteta) {
      if (!IsSelectedTrack(track, DielectronCuts::kTrackPtRange)) {
        return false;
      }
      if (!IsSelectedTrack(track, DielectronCuts::kTrackEtaRange)) {
        return false;
      }
    }

    if (!IsSelectedTrack(track, DielectronCuts::kTrackPhiRange)) {
      return false;
    }
    if (!IsSelectedTrack(track, DielectronCuts::kDCA3Dsigma)) {
      return false;
    }
    if (!IsSelectedTrack(track, DielectronCuts::kDCAxy)) {
      return false;
    }
    if (!IsSelectedTrack(track, DielectronCuts::kDCAz)) {
      return false;
    }

    // ITS cuts
    if (!IsSelectedTrack(track, DielectronCuts::kITSNCls)) {
      return false;
    }
    if (!IsSelectedTrack(track, DielectronCuts::kITSChi2NDF)) {
      return false;
    }
    if (!IsSelectedTrack(track, DielectronCuts::kITSClusterSize)) {
      return false;
    }

    if (mRequireITSibAny) {
      auto hits_ib = std::count_if(its_ib_any_Requirement.second.begin(), its_ib_any_Requirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
      if (hits_ib < its_ib_any_Requirement.first) {
        return false;
      }
    }

    if (mRequireITSib1st) {
      auto hits_ib = std::count_if(its_ib_1st_Requirement.second.begin(), its_ib_1st_Requirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
      if (hits_ib < its_ib_1st_Requirement.first) {
        return false;
      }
    }

    // TPC cuts
    if (!IsSelectedTrack(track, DielectronCuts::kTPCNCls)) {
      return false;
    }
    if (!IsSelectedTrack(track, DielectronCuts::kTPCCrossedRows)) {
      return false;
    }
    if (!IsSelectedTrack(track, DielectronCuts::kTPCCrossedRowsOverNCls)) {
      return false;
    }
    if (!IsSelectedTrack(track, DielectronCuts::kTPCFracSharedClusters)) {
      return false;
    }
    if (!IsSelectedTrack(track, DielectronCuts::kRelDiffPin)) {
      return false;
    }
    if (!IsSelectedTrack(track, DielectronCuts::kTPCChi2NDF)) {
      return false;
    }

    if (mApplyPF && !IsSelectedTrack(track, DielectronCuts::kPrefilter)) {
      return false;
    }

    // PID cuts
    if constexpr (isML) {
      if (!PassPIDML(track, collision)) {
        return false;
      }
    } else {
      if (!PassPID(track)) {
        return false;
      }
    }

    return true;
  }

  template <typename TTrack, typename TCollision>
  bool PassPIDML(TTrack const& track, TCollision const& collision) const
  {
    /*if (!PassTOFif(track)) { // Allows for pre-selection. But potentially dangerous if analyzers are not aware of it
      return false;
    }*/
    std::vector<float> inputFeatures = mPIDMlResponse->getInputFeatures(track, collision);
    float binningFeature = mPIDMlResponse->getBinningFeature(track, collision);
    return mPIDMlResponse->isSelectedMl(inputFeatures, binningFeature);
  }

  template <typename T>
  bool PassPID(T const& track) const
  {
    switch (mPIDScheme) {
      case static_cast<int>(PIDSchemes::kTOFreq):
        return PassTOFreq(track);

      case static_cast<int>(PIDSchemes::kTPChadrej):
        return PassTPChadrej(track);

      case static_cast<int>(PIDSchemes::kTPChadrejORTOFreq):
        return PassTPChadrej(track) || PassTOFreq(track);

      case static_cast<int>(PIDSchemes::kTPConly):
        return PassTPConly(track);

      case static_cast<int>(PIDSchemes::kTOFif):
        return PassTOFif(track);

      case static_cast<int>(PIDSchemes::kPIDML):
        return true; // don't use kPIDML here.

      case static_cast<int>(PIDSchemes::kTPChadrejORTOFreq_woTOFif):
        return PassTPConlyhadrej(track) || PassTOFreq(track);

      case static_cast<int>(PIDSchemes::kUnDef):
        return true;

      default:
        return true;
    }
  }

  template <typename T>
  bool PassTOFreq(T const& track) const
  {
    bool is_el_included_TPC = mMinTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < mMaxTPCNsigmaEl;
    bool is_pi_excluded_TPC = track.tpcInnerParam() < mMaxPinForPionRejectionTPC ? (track.tpcNSigmaPi() < mMinTPCNsigmaPi || mMaxTPCNsigmaPi < track.tpcNSigmaPi()) : true;
    bool is_el_included_TOF = (mMinTOFNsigmaEl < track.tofNSigmaEl() && track.tofNSigmaEl() < mMaxTOFNsigmaEl) && (track.hasTOF() && track.tofChi2() < mMaxChi2TOF);
    bool is_ka_excluded_ITS = (mMinP_ITSNsigmaKa < track.p() && track.p() < mMaxP_ITSNsigmaKa) ? (track.itsNSigmaKa() < mMinITSNsigmaKa || mMaxITSNsigmaKa < track.itsNSigmaKa()) : true;
    bool is_pr_excluded_ITS = (mMinP_ITSNsigmaPr < track.p() && track.p() < mMaxP_ITSNsigmaPr) ? (track.itsNSigmaPr() < mMinITSNsigmaPr || mMaxITSNsigmaPr < track.itsNSigmaPr()) : true;
    return is_el_included_TPC && is_pi_excluded_TPC && is_el_included_TOF && is_ka_excluded_ITS && is_pr_excluded_ITS;
  }

  template <typename T>
  bool PassTPChadrej(T const& track) const
  {
    bool is_el_included_TPC = mMinTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < mMaxTPCNsigmaEl;
    bool is_mu_excluded_TPC = mMuonExclusionTPC ? track.tpcNSigmaMu() < mMinTPCNsigmaMu || mMaxTPCNsigmaMu < track.tpcNSigmaMu() : true;
    bool is_pi_excluded_TPC = track.tpcInnerParam() < mMaxPinForPionRejectionTPC ? (track.tpcNSigmaPi() < mMinTPCNsigmaPi || mMaxTPCNsigmaPi < track.tpcNSigmaPi()) : true;
    bool is_ka_excluded_TPC = track.tpcNSigmaKa() < mMinTPCNsigmaKa || mMaxTPCNsigmaKa < track.tpcNSigmaKa();
    bool is_pr_excluded_TPC = track.tpcNSigmaPr() < mMinTPCNsigmaPr || mMaxTPCNsigmaPr < track.tpcNSigmaPr();
    bool is_el_included_TOF = track.hasTOF() ? (mMinTOFNsigmaEl < track.tofNSigmaEl() && track.tofNSigmaEl() < mMaxTOFNsigmaEl && track.tofChi2() < mMaxChi2TOF) : true;
    bool is_ka_excluded_ITS = (mMinP_ITSNsigmaKa < track.p() && track.p() < mMaxP_ITSNsigmaKa) ? (track.itsNSigmaKa() < mMinITSNsigmaKa || mMaxITSNsigmaKa < track.itsNSigmaKa()) : true;
    bool is_pr_excluded_ITS = (mMinP_ITSNsigmaPr < track.p() && track.p() < mMaxP_ITSNsigmaPr) ? (track.itsNSigmaPr() < mMinITSNsigmaPr || mMaxITSNsigmaPr < track.itsNSigmaPr()) : true;
    return is_el_included_TPC && is_mu_excluded_TPC && is_pi_excluded_TPC && is_ka_excluded_TPC && is_pr_excluded_TPC && is_el_included_TOF && is_ka_excluded_ITS && is_pr_excluded_ITS;
  }

  template <typename T>
  bool PassTPConly(T const& track) const
  {
    bool is_el_included_TPC = mMinTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < mMaxTPCNsigmaEl;
    bool is_ka_excluded_ITS = (mMinP_ITSNsigmaKa < track.p() && track.p() < mMaxP_ITSNsigmaKa) ? (track.itsNSigmaKa() < mMinITSNsigmaKa || mMaxITSNsigmaKa < track.itsNSigmaKa()) : true;
    bool is_pr_excluded_ITS = (mMinP_ITSNsigmaPr < track.p() && track.p() < mMaxP_ITSNsigmaPr) ? (track.itsNSigmaPr() < mMinITSNsigmaPr || mMaxITSNsigmaPr < track.itsNSigmaPr()) : true;
    return is_el_included_TPC && is_ka_excluded_ITS && is_pr_excluded_ITS;
  }

  template <typename T>
  bool PassTPConlyhadrej(T const& track) const
  {
    bool is_el_included_TPC = mMinTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < mMaxTPCNsigmaEl;
    bool is_mu_excluded_TPC = mMuonExclusionTPC ? track.tpcNSigmaMu() < mMinTPCNsigmaMu || mMaxTPCNsigmaMu < track.tpcNSigmaMu() : true;
    bool is_pi_excluded_TPC = track.tpcInnerParam() < mMaxPinForPionRejectionTPC ? (track.tpcNSigmaPi() < mMinTPCNsigmaPi || mMaxTPCNsigmaPi < track.tpcNSigmaPi()) : true;
    bool is_ka_excluded_TPC = track.tpcNSigmaKa() < mMinTPCNsigmaKa || mMaxTPCNsigmaKa < track.tpcNSigmaKa();
    bool is_pr_excluded_TPC = track.tpcNSigmaPr() < mMinTPCNsigmaPr || mMaxTPCNsigmaPr < track.tpcNSigmaPr();
    bool is_ka_excluded_ITS = (mMinP_ITSNsigmaKa < track.p() && track.p() < mMaxP_ITSNsigmaKa) ? (track.itsNSigmaKa() < mMinITSNsigmaKa || mMaxITSNsigmaKa < track.itsNSigmaKa()) : true;
    bool is_pr_excluded_ITS = (mMinP_ITSNsigmaPr < track.p() && track.p() < mMaxP_ITSNsigmaPr) ? (track.itsNSigmaPr() < mMinITSNsigmaPr || mMaxITSNsigmaPr < track.itsNSigmaPr()) : true;
    return is_el_included_TPC && is_mu_excluded_TPC && is_pi_excluded_TPC && is_ka_excluded_TPC && is_pr_excluded_TPC && is_ka_excluded_ITS && is_pr_excluded_ITS;
  }

  template <typename T>
  bool PassTOFif(T const& track) const
  {
    bool is_el_included_TPC = mMinTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < mMaxTPCNsigmaEl;
    bool is_pi_excluded_TPC = track.tpcInnerParam() < mMaxPinForPionRejectionTPC ? (track.tpcNSigmaPi() < mMinTPCNsigmaPi || mMaxTPCNsigmaPi < track.tpcNSigmaPi()) : true;
    bool is_el_included_TOF = track.hasTOF() ? (mMinTOFNsigmaEl < track.tofNSigmaEl() && track.tofNSigmaEl() < mMaxTOFNsigmaEl && track.tofChi2() < mMaxChi2TOF) : true;
    bool is_ka_excluded_ITS = (mMinP_ITSNsigmaKa < track.p() && track.p() < mMaxP_ITSNsigmaKa) ? (track.itsNSigmaKa() < mMinITSNsigmaKa || mMaxITSNsigmaKa < track.itsNSigmaKa()) : true;
    bool is_pr_excluded_ITS = (mMinP_ITSNsigmaPr < track.p() && track.p() < mMaxP_ITSNsigmaPr) ? (track.itsNSigmaPr() < mMinITSNsigmaPr || mMaxITSNsigmaPr < track.itsNSigmaPr()) : true;
    return is_el_included_TPC && is_pi_excluded_TPC && is_el_included_TOF && is_ka_excluded_ITS && is_pr_excluded_ITS;
  }

  template <typename T>
  bool IsSelectedTrack(T const& track, const DielectronCuts& cut) const
  {
    switch (cut) {
      case DielectronCuts::kTrackPtRange:
        return track.pt() >= mMinTrackPt && track.pt() <= mMaxTrackPt;

      case DielectronCuts::kTrackEtaRange:
        return track.eta() >= mMinTrackEta && track.eta() <= mMaxTrackEta;

      case DielectronCuts::kTrackPhiRange:
        return track.phi() >= mMinTrackPhi && track.phi() <= mMaxTrackPhi;

      case DielectronCuts::kTPCNCls:
        return track.tpcNClsFound() >= mMinNClustersTPC;

      case DielectronCuts::kTPCCrossedRows:
        return track.tpcNClsCrossedRows() >= mMinNCrossedRowsTPC;

      case DielectronCuts::kTPCCrossedRowsOverNCls:
        return track.tpcCrossedRowsOverFindableCls() >= mMinNCrossedRowsOverFindableClustersTPC;

      case DielectronCuts::kTPCFracSharedClusters:
        return track.tpcFractionSharedCls() <= mMaxFracSharedClustersTPC;

      case DielectronCuts::kRelDiffPin:
        return mMinRelDiffPin < (track.tpcInnerParam() - track.p()) / track.p() && (track.tpcInnerParam() - track.p()) / track.p() < mMaxRelDiffPin;

      case DielectronCuts::kTPCChi2NDF:
        return mMinChi2PerClusterTPC < track.tpcChi2NCl() && track.tpcChi2NCl() < mMaxChi2PerClusterTPC;

      case DielectronCuts::kDCA3Dsigma:
        return mMinDca3D <= dca3DinSigma(track) && dca3DinSigma(track) <= mMaxDca3D; // in sigma for single leg

      case DielectronCuts::kDCAxy:
        return abs(track.dcaXY()) <= ((mMaxDcaXYPtDep) ? mMaxDcaXYPtDep(track.pt()) : mMaxDcaXY);

      case DielectronCuts::kDCAz:
        return abs(track.dcaZ()) <= mMaxDcaZ;

      case DielectronCuts::kITSNCls:
        return mMinNClustersITS <= track.itsNCls() && track.itsNCls() <= mMaxNClustersITS;

      case DielectronCuts::kITSChi2NDF:
        return mMinChi2PerClusterITS < track.itsChi2NCl() && track.itsChi2NCl() < mMaxChi2PerClusterITS;

      case DielectronCuts::kITSClusterSize:
        return ((mMinP_ITSClusterSize < track.p() && track.p() < mMaxP_ITSClusterSize) ? (mMinMeanClusterSizeITS < track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())) && track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())) < mMaxMeanClusterSizeITS) : true);

      case DielectronCuts::kPrefilter:
        return track.pfb() <= 0;

      default:
        return false;
    }
  }

  // Setters
  void SetPairPtRange(float minPt = 0.f, float maxPt = 1e10f);
  void SetPairYRange(float minY = -1e10f, float maxY = 1e10f);
  void SetPairDCARange(float min = 0.f, float max = 1e10f); // 3D DCA in sigma
  void SetMeeRange(float min = 0.f, float max = 0.5);
  void SetPairOpAng(float minOpAng = 0.f, float maxOpAng = 1e10f);
  void SetMaxMeePhiVDep(std::function<float(float)> phivDepCut, float min_phiv, float max_phiv);
  void SelectPhotonConversion(bool flag);
  void SetMindEtadPhi(bool flag, float min_deta, float min_dphi);
  void SetRequireDifferentSides(bool flag);

  void SetTrackPtRange(float minPt = 0.f, float maxPt = 1e10f);
  void SetTrackEtaRange(float minEta = -1e10f, float maxEta = 1e10f);
  void SetTrackPhiRange(float minPhi = 0.f, float maxPhi = 2.f * M_PI);
  void SetMinNClustersTPC(int minNClustersTPC);
  void SetMinNCrossedRowsTPC(int minNCrossedRowsTPC);
  void SetMinNCrossedRowsOverFindableClustersTPC(float minNCrossedRowsOverFindableClustersTPC);
  void SetMaxFracSharedClustersTPC(float max);
  void SetRelDiffPin(float min, float max);
  void SetChi2PerClusterTPC(float min, float max);
  void SetNClustersITS(int min, int max);
  void SetChi2PerClusterITS(float min, float max);
  void SetMeanClusterSizeITS(float min, float max, float minP = 0.f, float maxP = 0.f);
  void SetChi2TOF(float min, float max);

  void SetPIDScheme(int scheme);
  void SetMinPinTOF(float min);
  void SetMuonExclusionTPC(bool flag);
  void SetTOFbetaRange(float min, float max);
  void SetTPCNsigmaElRange(float min, float max);
  void SetTPCNsigmaMuRange(float min, float max);
  void SetTPCNsigmaPiRange(float min, float max);
  void SetTPCNsigmaKaRange(float min, float max);
  void SetTPCNsigmaPrRange(float min, float max);
  void SetTOFNsigmaElRange(float min, float max);
  void SetTOFNsigmaMuRange(float min, float max);
  void SetTOFNsigmaPiRange(float min, float max);
  void SetTOFNsigmaKaRange(float min, float max);
  void SetTOFNsigmaPrRange(float min, float max);
  void SetITSNsigmaElRange(float min, float max);
  void SetITSNsigmaMuRange(float min, float max);
  void SetITSNsigmaPiRange(float min, float max);
  void SetITSNsigmaKaRange(float min, float max);
  void SetITSNsigmaPrRange(float min, float max);

  void SetPRangeForITSNsigmaKa(float min, float max);
  void SetPRangeForITSNsigmaPr(float min, float max);

  void SetMaxPinMuonTPConly(float max);
  void SetMaxPinForPionRejectionTPC(float max);
  void RequireITSibAny(bool flag);
  void RequireITSib1st(bool flag);

  void SetTrackDca3DRange(float min, float max); // in sigma
  void SetTrackMaxDcaXY(float maxDcaXY);         // in cm
  void SetTrackMaxDcaZ(float maxDcaZ);           // in cm
  void SetTrackMaxDcaXYPtDep(std::function<float(float)> ptDepCut);
  void ApplyPrefilter(bool flag);
  void ApplyPhiV(bool flag);

  void SetPIDMlResponse(o2::analysis::MlResponseDielectronSingleTrack<float>* mlResponse)
  {
    mPIDMlResponse = mlResponse;
  }

  // Getters
  bool IsPhotonConversionSelected() const { return mSelectPC; }

 private:
  static const std::pair<int8_t, std::set<uint8_t>> its_ib_any_Requirement;
  static const std::pair<int8_t, std::set<uint8_t>> its_ib_1st_Requirement;
  // pair cuts
  float mMinMee{0.f}, mMaxMee{1e10f};
  float mMinPairPt{0.f}, mMaxPairPt{1e10f};       // range in pT
  float mMinPairY{-1e10f}, mMaxPairY{1e10f};      // range in rapidity
  float mMinPairDCA3D{0.f}, mMaxPairDCA3D{1e10f}; // range in 3D DCA in sigma
  float mMinPhivPair{0.f}, mMaxPhivPair{+3.2};
  std::function<float(float)> mMaxMeePhiVDep{}; // max mee as a function of phiv
  bool mSelectPC{false};                        // flag to select photon conversion used in mMaxPhivPairMeeDep
  bool mApplydEtadPhi{false};                   // flag to apply deta, dphi cut between 2 tracks
  float mMinDeltaEta{0.f};
  float mMinDeltaPhi{0.f};
  float mMinOpAng{0.f}, mMaxOpAng{1e10f};
  bool mRequireDiffSides{false}; // flag to require 2 tracks to be from different sides. (A-C combination). If one wants 2 tracks to be in the same side (A-A or C-C), one can simply use track eta cut.

  // kinematic cuts
  float mMinTrackPt{0.f}, mMaxTrackPt{1e10f};        // range in pT
  float mMinTrackEta{-1e10f}, mMaxTrackEta{1e10f};   // range in eta
  float mMinTrackPhi{0.f}, mMaxTrackPhi{2.f * M_PI}; // range in phi

  // track quality cuts
  int mMinNClustersTPC{0};                                           // min number of TPC clusters
  int mMinNCrossedRowsTPC{0};                                        // min number of crossed rows in TPC
  float mMinChi2PerClusterTPC{-1e10f}, mMaxChi2PerClusterTPC{1e10f}; // max tpc fit chi2 per TPC cluster
  float mMinNCrossedRowsOverFindableClustersTPC{0.f};                // min ratio crossed rows / findable clusters
  float mMaxFracSharedClustersTPC{999.f};                            // max ratio shared clusters / clusters in TPC
  float mMinRelDiffPin{-1e10f}, mMaxRelDiffPin{1e10f};               // max relative difference between p at TPC inner wall and p at PV
  int mMinNClustersITS{0}, mMaxNClustersITS{7};                      // range in number of ITS clusters
  float mMinChi2PerClusterITS{-1e10f}, mMaxChi2PerClusterITS{1e10f}; // max its fit chi2 per ITS cluster
  float mMaxPinMuonTPConly{0.2f};                                    // max pin cut for muon ID with TPConly
  float mMaxPinForPionRejectionTPC{1e10f};                           // max pin for pion rejection in TPC
  bool mRequireITSibAny{true};
  bool mRequireITSib1st{false};
  float mMinChi2TOF{-1e10f}, mMaxChi2TOF{1e10f}; // max tof chi2 per

  float mMinDca3D{0.0f};                        // min dca in 3D in units of sigma
  float mMaxDca3D{1e+10};                       // max dca in 3D in units of sigma
  float mMaxDcaXY{1.0f};                        // max dca in xy plane
  float mMaxDcaZ{1.0f};                         // max dca in z direction
  std::function<float(float)> mMaxDcaXYPtDep{}; // max dca in xy plane as function of pT
  bool mApplyPhiV{true};
  bool mApplyPF{false};
  float mMinMeanClusterSizeITS{-1e10f}, mMaxMeanClusterSizeITS{1e10f}; // max <its cluster size> x cos(Lmabda)
  float mMinP_ITSClusterSize{0.0}, mMaxP_ITSClusterSize{0.0};

  // pid cuts
  int mPIDScheme{-1};
  float mMinPinTOF{0.0f};        // min pin cut for TOF.
  bool mMuonExclusionTPC{false}; // flag to reject muon in TPC for low B
  float mMinTOFbeta{-999}, mMaxTOFbeta{999};
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

  float mMinITSNsigmaEl{-1e+10}, mMaxITSNsigmaEl{+1e+10};
  float mMinITSNsigmaMu{-1e+10}, mMaxITSNsigmaMu{+1e+10};
  float mMinITSNsigmaPi{-1e+10}, mMaxITSNsigmaPi{+1e+10};
  float mMinITSNsigmaKa{-1e+10}, mMaxITSNsigmaKa{+1e+10};
  float mMinITSNsigmaPr{-1e+10}, mMaxITSNsigmaPr{+1e+10};
  float mMinP_ITSNsigmaKa{0.0}, mMaxP_ITSNsigmaKa{0.0};
  float mMinP_ITSNsigmaPr{0.0}, mMaxP_ITSNsigmaPr{0.0};

  o2::analysis::MlResponseDielectronSingleTrack<float>* mPIDMlResponse{nullptr};

  ClassDef(DielectronCut, 1);
};

#endif // PWGEM_DILEPTON_CORE_DIELECTRONCUT_H_
