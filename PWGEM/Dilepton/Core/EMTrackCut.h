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
// Class for track selection
//

#ifndef PWGEM_DILEPTON_CORE_EMTRACKCUT_H_
#define PWGEM_DILEPTON_CORE_EMTRACKCUT_H_

#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/DataTypes.h"
#include "Framework/Logger.h"

#include "Math/Vector4D.h"
#include "TNamed.h"

#include <algorithm>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;

class EMTrackCut : public TNamed
{
 public:
  EMTrackCut() = default;
  EMTrackCut(const char* name, const char* title) : TNamed(name, title) {}
  ~EMTrackCut() {}

  enum class EMTrackCuts : int {
    // track cut
    kTrackPtRange,
    kTrackEtaRange,
    kTrackPhiRange,
    // kDCAxy,
    // kDCAz,
    // kTPCNCls,
    // kTPCCrossedRows,
    // kTPCCrossedRowsOverNCls,
    // kTPCFracSharedClusters,
    // kTPCChi2NDF,
    // kITSNCls,
    // kITSChi2NDF,
    kTrackBit,
    kNCuts
  };

  template <typename TTrack>
  bool IsSelected(TTrack const& track) const
  {
    // if (!track.hasITS() || !track.hasTPC()) {
    //   return false;
    // }

    if (!IsSelectedTrack(track, EMTrackCuts::kTrackPtRange)) {
      return false;
    }
    if (!IsSelectedTrack(track, EMTrackCuts::kTrackEtaRange)) {
      return false;
    }
    if (!IsSelectedTrack(track, EMTrackCuts::kTrackPhiRange)) {
      return false;
    }

    // if (!IsSelectedTrack(track, EMTrackCuts::kDCAxy)) {
    //   return false;
    // }
    // if (!IsSelectedTrack(track, EMTrackCuts::kDCAz)) {
    //   return false;
    // }

    if (!IsSelectedTrack(track, EMTrackCuts::kTrackBit)) {
      return false;
    }

    //    // ITS cuts
    //    if (!IsSelectedTrack(track, EMTrackCuts::kITSNCls)) {
    //      return false;
    //    }
    //    if (!IsSelectedTrack(track, EMTrackCuts::kITSChi2NDF)) {
    //      return false;
    //    }
    //
    //    if (mRequireITSibAny) {
    //      auto hits_ib = std::count_if(its_ib_any_Requirement.second.begin(), its_ib_any_Requirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
    //      if (hits_ib < its_ib_any_Requirement.first) {
    //        return false;
    //      }
    //    }
    //
    //    if (mRequireITSib1st) {
    //      auto hits_ib = std::count_if(its_ib_1st_Requirement.second.begin(), its_ib_1st_Requirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
    //      if (hits_ib < its_ib_1st_Requirement.first) {
    //        return false;
    //      }
    //    }
    //
    //    // TPC cuts
    //    if (!IsSelectedTrack(track, EMTrackCuts::kTPCNCls)) {
    //      return false;
    //    }
    //    if (!IsSelectedTrack(track, EMTrackCuts::kTPCCrossedRows)) {
    //      return false;
    //    }
    //    if (!IsSelectedTrack(track, EMTrackCuts::kTPCCrossedRowsOverNCls)) {
    //      return false;
    //    }
    //    if (!IsSelectedTrack(track, EMTrackCuts::kTPCFracSharedClusters)) {
    //      return false;
    //    }
    //    if (!IsSelectedTrack(track, EMTrackCuts::kTPCChi2NDF)) {
    //      return false;
    //    }

    return true;
  }

  template <typename T>
  bool IsSelectedTrack(T const& track, const EMTrackCuts& cut) const
  {
    switch (cut) {
      case EMTrackCuts::kTrackPtRange:
        return track.pt() > mMinTrackPt && track.pt() < mMaxTrackPt;

      case EMTrackCuts::kTrackEtaRange:
        return track.eta() > mMinTrackEta && track.eta() < mMaxTrackEta;

      case EMTrackCuts::kTrackPhiRange:
        return track.phi() > mMinTrackPhi && track.phi() < mMaxTrackPhi;

        // case EMTrackCuts::kDCAxy:
        //   return std::fabs(track.dcaXY()) < ((mMaxDcaXYPtDep) ? mMaxDcaXYPtDep(track.pt()) : mMaxDcaXY);

        // case EMTrackCuts::kDCAz:
        //   return std::fabs(track.dcaZ()) < mMaxDcaZ;

      case EMTrackCuts::kTrackBit: {
        // for (int i = 0; i < 10; i++) {
        //   if ((mTrackBit & (1 << i)) > 0 && !((track.trackBit() & (1 << i)) > 0)) {
        //     return false;
        //   }
        // }
        // return true;
        return (track.trackBit() & mTrackBit) >= mTrackBit;
      }

        // case EMTrackCuts::kTPCNCls:
        //   return track.tpcNClsFound() >= mMinNClustersTPC;

        // case EMTrackCuts::kTPCCrossedRows:
        //   return track.tpcNClsCrossedRows() >= mMinNCrossedRowsTPC;

        // case EMTrackCuts::kTPCCrossedRowsOverNCls:
        //   return track.tpcCrossedRowsOverFindableCls() > mMinNCrossedRowsOverFindableClustersTPC;

        // case EMTrackCuts::kTPCFracSharedClusters:
        //   return track.tpcFractionSharedCls() < mMaxFracSharedClustersTPC;

        // case EMTrackCuts::kTPCChi2NDF:
        //   return mMinChi2PerClusterTPC < track.tpcChi2NCl() && track.tpcChi2NCl() < mMaxChi2PerClusterTPC;

        // case EMTrackCuts::kITSNCls:
        //   return mMinNClustersITS <= track.itsNCls() && track.itsNCls() <= mMaxNClustersITS;

        // case EMTrackCuts::kITSChi2NDF:
        //   return mMinChi2PerClusterITS < track.itsChi2NCl() && track.itsChi2NCl() < mMaxChi2PerClusterITS;

      default:
        return false;
    }
  }

  // Setters
  void SetTrackPtRange(float minPt = 0.f, float maxPt = 1e10f);
  void SetTrackEtaRange(float minEta = -1e10f, float maxEta = 1e10f);
  void SetTrackPhiRange(float minPhi = 0.f, float maxPhi = 6.3f);
  void SetMinNClustersTPC(int minNClustersTPC);
  void SetMinNCrossedRowsTPC(int minNCrossedRowsTPC);
  void SetMinNCrossedRowsOverFindableClustersTPC(float minNCrossedRowsOverFindableClustersTPC);
  void SetMaxFracSharedClustersTPC(float max);
  void SetChi2PerClusterTPC(float min, float max);
  void SetNClustersITS(int min, int max);
  void SetChi2PerClusterITS(float min, float max);

  void SetTrackDca3DRange(float min, float max); // in sigma
  void SetTrackMaxDcaXY(float maxDcaXY);         // in cm
  void SetTrackMaxDcaZ(float maxDcaZ);           // in cm
  void SetTrackMaxDcaXYPtDep(std::function<float(float)> ptDepCut);
  void RequireITSibAny(bool flag);
  void RequireITSib1st(bool flag);
  void SetTrackBit(uint16_t bits);

 private:
  static const std::pair<int8_t, std::set<uint8_t>> its_ib_any_Requirement;
  static const std::pair<int8_t, std::set<uint8_t>> its_ib_1st_Requirement;

  // kinematic cuts
  float mMinTrackPt{0.f}, mMaxTrackPt{1e10f};      // range in pT
  float mMinTrackEta{-1e10f}, mMaxTrackEta{1e10f}; // range in eta
  float mMinTrackPhi{0.f}, mMaxTrackPhi{6.3};      // range in phi

  // track quality cuts
  int mMinNClustersTPC{0};                                        // min number of TPC clusters
  int mMinNCrossedRowsTPC{0};                                     // min number of crossed rows in TPC
  float mMinChi2PerClusterTPC{0.f}, mMaxChi2PerClusterTPC{1e10f}; // max tpc fit chi2 per TPC cluster
  float mMinNCrossedRowsOverFindableClustersTPC{0.f};             // min ratio crossed rows / findable clusters
  float mMaxFracSharedClustersTPC{999.f};                         // max ratio shared clusters / clusters in TPC
  int mMinNClustersITS{0}, mMaxNClustersITS{7};                   // range in number of ITS clusters
  float mMinChi2PerClusterITS{0.f}, mMaxChi2PerClusterITS{1e10f}; // max its fit chi2 per ITS cluster
  bool mRequireITSibAny{true};
  bool mRequireITSib1st{false};
  uint16_t mTrackBit{0};

  float mMaxDcaXY{1.0f};                        // max dca in xy plane
  float mMaxDcaZ{1.0f};                         // max dca in z direction
  std::function<float(float)> mMaxDcaXYPtDep{}; // max dca in xy plane as function of pT

  ClassDef(EMTrackCut, 1);
};

#endif // PWGEM_DILEPTON_CORE_EMTRACKCUT_H_
