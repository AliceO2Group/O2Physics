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
#include "TNamed.h"
#include "Math/Vector4D.h"

#include "Tools/ML/MlResponse.h"
#include "Tools/ML/model.h"

#include "Framework/Logger.h"
#include "Framework/DataTypes.h"
#include "CommonConstants/PhysicsConstants.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"
#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"

using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;

class DalitzEECut : public TNamed
{
 public:
  DalitzEECut() = default;
  DalitzEECut(const char* name, const char* title) : TNamed(name, title) {}

  enum class DalitzEECuts : int {
    // pair cut
    kMee = 0,
    kPairPtRange,
    kPairYRange,
    kPhiV,
    // track cut
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
    kITSCluserSize,
    kNCuts
  };
  static const char* mCutNames[static_cast<int>(DalitzEECuts::kNCuts)];

  enum class PIDSchemes : int {
    kUnDef = -1,
    kTPConly = 0,
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

  template <typename TTrack1, typename TTrack2>
  bool IsSelectedPair(TTrack1 const& t1, TTrack2 const& t2, const float bz) const
  {
    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
    float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), bz);

    if (v12.M() < mMinMee || mMaxMee < v12.M()) {
      return false;
    }

    if (v12.Rapidity() < mMinPairY || mMaxPairY < v12.Rapidity()) {
      return false;
    }

    if (mApplyPhiV && ((phiv < mMinPhivPair || (mMaxPhivPairMeeDep ? mMaxPhivPairMeeDep(v12.M()) : mMaxPhivPair) < phiv) ^ mSelectPC)) {
      return false;
    }
    return true;
  }

  template <bool isML = false, typename TTrack, typename TCollision = int>
  bool IsSelectedTrack(TTrack const& track, TCollision const& = 0) const
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
    if (!IsSelectedTrack(track, DalitzEECuts::kITSCluserSize)) {
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

    // PID cuts
    if (!PassPID(track)) {
      return false;
    }

    return true;
  }

  template <typename T>
  bool PassPID(T const& track) const
  {
    switch (mPIDScheme) {
      case static_cast<int>(PIDSchemes::kTPConly):
        return PassTPConly(track);

      case static_cast<int>(PIDSchemes::kUnDef):
        return true;

      default:
        return true;
    }
  }

  template <typename T>
  bool PassTPConly(T const& track) const
  {
    bool is_el_included_TPC = mMinTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < mMaxTPCNsigmaEl;
    bool is_pi_excluded_TPC = track.tpcNSigmaPi() < mMinTPCNsigmaPi || mMaxTPCNsigmaPi < track.tpcNSigmaPi();
    return is_el_included_TPC && is_pi_excluded_TPC;
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

      case DalitzEECuts::kITSCluserSize:
        return mMinMeanClusterSizeITS < track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())) && track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())) < mMaxMeanClusterSizeITS;

      default:
        return false;
    }
  }

  // Setters
  void SetPairPtRange(float minPt = 0.f, float maxPt = 1e10f);
  void SetPairYRange(float minY = -1e10f, float maxY = 1e10f);
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
  void SetMeanClusterSizeITS(float min, float max);

  void SetPIDScheme(int scheme);
  void SetTPCNsigmaElRange(float min = -1e+10, float max = 1e+10);
  void SetTPCNsigmaPiRange(float min = -1e+10, float max = 1e+10);
  void RequireITSibAny(bool flag);
  void RequireITSib1st(bool flag);

  void SetMaxDcaXY(float maxDcaXY); // in cm
  void SetMaxDcaZ(float maxDcaZ);   // in cm
  void SetMaxDcaXYPtDep(std::function<float(float)> ptDepCut);
  void ApplyPrefilter(bool flag);
  void ApplyPhiV(bool flag);

  // Getters
  bool IsPhotonConversionSelected() const { return mSelectPC; }

  /// @brief Print the track selection
  // void print() const;

 private:
  static const std::pair<int8_t, std::set<uint8_t>> its_ib_any_Requirement;
  static const std::pair<int8_t, std::set<uint8_t>> its_ib_1st_Requirement;
  // pair cuts
  float mMinMee{0.f}, mMaxMee{1e10f};
  float mMinPairPt{0.f}, mMaxPairPt{1e10f};  // range in pT
  float mMinPairY{-1e10f}, mMaxPairY{1e10f}; // range in rapidity
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
  bool mRequireITSibAny{true};
  bool mRequireITSib1st{false};

  float mMaxDcaXY{1.0f};                        // max dca in xy plane
  float mMaxDcaZ{1.0f};                         // max dca in z direction
  std::function<float(float)> mMaxDcaXYPtDep{}; // max dca in xy plane as function of pT
  bool mApplyPhiV{true};
  float mMinMeanClusterSizeITS{-1e10f}, mMaxMeanClusterSizeITS{1e10f}; // max <its cluster size> x cos(Lmabda)

  // pid cuts
  int mPIDScheme{-1};
  float mMinTPCNsigmaEl{-1e+10}, mMaxTPCNsigmaEl{+1e+10};
  float mMinTPCNsigmaPi{-1e+10}, mMaxTPCNsigmaPi{+1e+10};

  o2::ml::OnnxModel* mPIDModel{nullptr};

  ClassDef(DalitzEECut, 1);
};

#endif // PWGEM_PHOTONMESON_CORE_DALITZEECUT_H_
