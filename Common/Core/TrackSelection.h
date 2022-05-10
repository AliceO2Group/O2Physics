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

#ifndef TrackSelection_H
#define TrackSelection_H

#include "Framework/Logger.h"
#include "Framework/DataTypes.h"
#include <set>
#include <vector>
#include "Rtypes.h"

class TrackSelection
{
 public:
  TrackSelection() = default;

  enum class TrackCuts : int {
    kTrackType = 0,
    kPtRange,
    kEtaRange,
    kTPCNCls,
    kTPCCrossedRows,
    kTPCCrossedRowsOverNCls,
    kTPCChi2NDF,
    kTPCRefit,
    kITSNCls,
    kITSChi2NDF,
    kITSRefit,
    kITSHits,
    kGoldenChi2,
    kDCAxy,
    kDCAz,
    kNCuts
  };

  static const std::string mCutNames[static_cast<int>(TrackCuts::kNCuts)];

  // Temporary function to check if track passes selection criteria. To be replaced by framework filters.
  template <typename T>
  bool IsSelected(T const& track)
  {
    if (!IsSelected(track, TrackCuts::kTrackType)) {
      return false;
    }
    if (!IsSelected(track, TrackCuts::kPtRange)) {
      return false;
    }
    if (!IsSelected(track, TrackCuts::kEtaRange)) {
      return false;
    }
    if (!IsSelected(track, TrackCuts::kTPCNCls)) {
      return false;
    }
    if (!IsSelected(track, TrackCuts::kTPCCrossedRows)) {
      return false;
    }
    if (!IsSelected(track, TrackCuts::kTPCCrossedRowsOverNCls)) {
      return false;
    }
    if (!IsSelected(track, TrackCuts::kTPCChi2NDF)) {
      return false;
    }
    if (!IsSelected(track, TrackCuts::kTPCRefit)) {
      return false;
    }
    if (!IsSelected(track, TrackCuts::kITSNCls)) {
      return false;
    }
    if (!IsSelected(track, TrackCuts::kITSChi2NDF)) {
      return false;
    }
    if (!IsSelected(track, TrackCuts::kITSRefit)) {
      return false;
    }
    if (!IsSelected(track, TrackCuts::kITSHits)) {
      return false;
    }
    if (!IsSelected(track, TrackCuts::kGoldenChi2)) {
      return false;
    }
    if (!IsSelected(track, TrackCuts::kDCAxy)) {
      return false;
    }
    if (!IsSelected(track, TrackCuts::kDCAz)) {
      return false;
    }
    return true;
  }

  // Temporary function to check if track passes and return a flag. To be replaced by framework filters.
  template <typename T>
  uint16_t IsSelectedMask(T const& track)
  {
    uint16_t flag = 0;

    auto setFlag = [&](const TrackCuts& cut) {
      if (IsSelected(track, cut)) {
        flag |= 1UL << static_cast<int>(cut);
      }
    };

    setFlag(TrackCuts::kTrackType);
    setFlag(TrackCuts::kPtRange);
    setFlag(TrackCuts::kEtaRange);
    setFlag(TrackCuts::kTPCNCls);
    setFlag(TrackCuts::kTPCCrossedRows);
    setFlag(TrackCuts::kTPCCrossedRowsOverNCls);
    setFlag(TrackCuts::kTPCChi2NDF);
    setFlag(TrackCuts::kTPCRefit);
    setFlag(TrackCuts::kITSNCls);
    setFlag(TrackCuts::kITSChi2NDF);
    setFlag(TrackCuts::kITSRefit);
    setFlag(TrackCuts::kITSHits);
    setFlag(TrackCuts::kGoldenChi2);
    setFlag(TrackCuts::kDCAxy);
    setFlag(TrackCuts::kDCAz);

    return flag;
  }

  // Temporary function to check if track passes a given selection criteria. To be replaced by framework filters.
  template <typename T>
  bool IsSelected(T const& track, const TrackCuts& cut)
  {
    const bool isRun2 = track.trackType() == o2::aod::track::Run2Track || track.trackType() == o2::aod::track::Run2Tracklet;

    switch (cut) {
      case TrackCuts::kTrackType:
        return track.trackType() == mTrackType;

      case TrackCuts::kPtRange:
        return track.pt() >= mMinPt && track.pt() <= mMaxPt;

      case TrackCuts::kEtaRange:
        return track.eta() >= mMinEta && track.eta() <= mMaxEta;

      case TrackCuts::kTPCNCls:
        return track.tpcNClsFound() >= mMinNClustersTPC;

      case TrackCuts::kTPCCrossedRows:
        return track.tpcNClsCrossedRows() >= mMinNCrossedRowsTPC;

      case TrackCuts::kTPCCrossedRowsOverNCls:
        return track.tpcCrossedRowsOverFindableCls() >= mMinNCrossedRowsOverFindableClustersTPC;

      case TrackCuts::kTPCChi2NDF:
        return track.tpcChi2NCl() <= mMaxChi2PerClusterTPC;

      case TrackCuts::kTPCRefit:
        return (mRequireTPCRefit ? (isRun2 ? (track.flags() & o2::aod::track::TPCrefit) : track.hasTPC()) : true);

      case TrackCuts::kITSNCls:
        return track.itsNCls() >= mMinNClustersITS;

      case TrackCuts::kITSChi2NDF:
        return track.itsChi2NCl() <= mMaxChi2PerClusterITS;

      case TrackCuts::kITSRefit:
        return (mRequireITSRefit ? (isRun2 ? (track.flags() & o2::aod::track::ITSrefit) : track.hasITS()) : true);

      case TrackCuts::kITSHits:
        return FulfillsITSHitRequirements(track.itsClusterMap());

      case TrackCuts::kGoldenChi2:
        return (isRun2 && mRequireGoldenChi2) ? (track.flags() & o2::aod::track::GoldenChi2) : true;

      case TrackCuts::kDCAxy:
        return abs(track.dcaXY()) <= ((mMaxDcaXYPtDep) ? mMaxDcaXYPtDep(track.pt()) : mMaxDcaXY);

      case TrackCuts::kDCAz:
        return abs(track.dcaZ()) <= mMaxDcaZ;

      default:
        return false;
    }
  }

  void SetTrackType(o2::aod::track::TrackTypeEnum trackType) { mTrackType = trackType; }
  void SetPtRange(float minPt = 0.f, float maxPt = 1e10f)
  {
    mMinPt = minPt;
    mMaxPt = maxPt;
  }
  void SetEtaRange(float minEta = -1e10f, float maxEta = 1e10f)
  {
    mMinEta = minEta;
    mMaxEta = maxEta;
  }
  void SetRequireITSRefit(bool requireITSRefit = true)
  {
    mRequireITSRefit = requireITSRefit;
  }
  void SetRequireTPCRefit(bool requireTPCRefit = true)
  {
    mRequireTPCRefit = requireTPCRefit;
  }
  void SetRequireGoldenChi2(bool requireGoldenChi2 = true)
  {
    mRequireGoldenChi2 = requireGoldenChi2;
  }
  void SetMinNClustersTPC(int minNClustersTPC)
  {
    mMinNClustersTPC = minNClustersTPC;
  }
  void SetMinNCrossedRowsTPC(int minNCrossedRowsTPC)
  {
    mMinNCrossedRowsTPC = minNCrossedRowsTPC;
  }
  void SetMinNCrossedRowsOverFindableClustersTPC(float minNCrossedRowsOverFindableClustersTPC)
  {
    mMinNCrossedRowsOverFindableClustersTPC = minNCrossedRowsOverFindableClustersTPC;
  }
  void SetMinNClustersITS(int minNClustersITS)
  {
    mMinNClustersITS = minNClustersITS;
  }
  void SetMaxChi2PerClusterTPC(float maxChi2PerClusterTPC)
  {
    mMaxChi2PerClusterTPC = maxChi2PerClusterTPC;
  }
  void SetMaxChi2PerClusterITS(float maxChi2PerClusterITS)
  {
    mMaxChi2PerClusterITS = maxChi2PerClusterITS;
  }
  void SetMaxDcaXY(float maxDcaXY) { mMaxDcaXY = maxDcaXY; }
  void SetMaxDcaZ(float maxDcaZ) { mMaxDcaZ = maxDcaZ; }

  void SetMaxDcaXYPtDep(std::function<float(float)> ptDepCut)
  {
    mMaxDcaXYPtDep = ptDepCut;
  }
  void SetRequireHitsInITSLayers(int8_t minNRequiredHits, std::set<uint8_t> requiredLayers)
  {
    // layer 0 corresponds to the the innermost ITS layer
    mRequiredITSHits.push_back(std::make_pair(minNRequiredHits, requiredLayers));
  }
  void SetRequireNoHitsInITSLayers(std::set<uint8_t> excludedLayers)
  {
    mRequiredITSHits.push_back(std::make_pair(-1, excludedLayers));
  }
  void ResetITSRequirements() { mRequiredITSHits.clear(); }

 private:
  bool FulfillsITSHitRequirements(uint8_t itsClusterMap);

  o2::aod::track::TrackTypeEnum mTrackType{o2::aod::track::TrackTypeEnum::Track};

  // kinematic cuts
  float mMinPt{0.f}, mMaxPt{1e10f};      // range in pT
  float mMinEta{-1e10f}, mMaxEta{1e10f}; // range in eta

  // track quality cuts
  int mMinNClustersTPC{0};                            // min number of TPC clusters
  int mMinNCrossedRowsTPC{0};                         // min number of crossed rows in TPC
  int mMinNClustersITS{0};                            // min number of ITS clusters
  float mMaxChi2PerClusterTPC{1e10f};                 // max tpc fit chi2 per TPC cluster
  float mMaxChi2PerClusterITS{1e10f};                 // max its fit chi2 per ITS cluster
  float mMinNCrossedRowsOverFindableClustersTPC{0.f}; // min ratio crossed rows / findable clusters

  float mMaxDcaXY{1e10f};                       // max dca in xy plane
  float mMaxDcaZ{1e10f};                        // max dca in z direction
  std::function<float(float)> mMaxDcaXYPtDep{}; // max dca in xy plane as function of pT

  bool mRequireITSRefit{false};   // require refit in ITS
  bool mRequireTPCRefit{false};   // require refit in TPC
  bool mRequireGoldenChi2{false}; // require golden chi2 cut (Run 2 only)

  // vector of ITS requirements (minNRequiredHits in specific requiredLayers)
  std::vector<std::pair<int8_t, std::set<uint8_t>>> mRequiredITSHits{};

  ClassDefNV(TrackSelection, 1);
};

#endif
