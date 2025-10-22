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

#include "Common/Core/TrackSelection.h"

#include <Framework/DataTypes.h>
#include <Framework/Logger.h>

#include <algorithm>
#include <cstdint>
#include <functional>
#include <set>
#include <string>
#include <utility>

bool TrackSelection::FulfillsITSHitRequirements(uint8_t itsClusterMap) const
{
  constexpr uint8_t bit = 1;
  for (auto& itsRequirement : mRequiredITSHits) {
    auto hits = std::count_if(itsRequirement.second.begin(), itsRequirement.second.end(), [&](auto&& requiredLayer) { return itsClusterMap & (bit << requiredLayer); });
    if ((itsRequirement.first == -1) && (hits > 0)) {
      return false; // no hits were required in specified layers
    } else if (hits < itsRequirement.first) {
      return false; // not enough hits found in specified layers
    }
  }
  return true;
}

const std::string TrackSelection::mCutNames[static_cast<int>(TrackSelection::TrackCuts::kNCuts)] = {"TrackType", "PtRange", "EtaRange", "TPCNCls", "TPCCrossedRows", "TPCCrossedRowsOverNCls", "TPCChi2NDF", "TPCRefit", "ITSNCls", "ITSChi2NDF", "ITSRefit", "ITSHits", "GoldenChi2", "DCAxy", "DCAz", "TPCFracSharedCls"};

void TrackSelection::SetTrackType(o2::aod::track::TrackTypeEnum trackType)
{
  mTrackType = trackType;
  LOG(info) << "Track selection, set track type: " << static_cast<int>(mTrackType);
}
void TrackSelection::SetPtRange(float minPt, float maxPt)
{
  mMinPt = minPt;
  mMaxPt = maxPt;
  LOG(info) << "Track selection, set pt range: " << mMinPt << " - " << mMaxPt;
}
void TrackSelection::SetEtaRange(float minEta, float maxEta)
{
  mMinEta = minEta;
  mMaxEta = maxEta;
  LOG(info) << "Track selection, set eta range: " << mMinEta << " - " << mMaxEta;
}
void TrackSelection::SetRequireITSRefit(bool requireITSRefit)
{
  mRequireITSRefit = requireITSRefit;
  LOG(info) << "Track selection, set require ITS refit: " << mRequireITSRefit;
}
void TrackSelection::SetRequireTPCRefit(bool requireTPCRefit)
{
  mRequireTPCRefit = requireTPCRefit;
  LOG(info) << "Track selection, set require TPC refit: " << mRequireTPCRefit;
}
void TrackSelection::SetRequireGoldenChi2(bool requireGoldenChi2)
{
  mRequireGoldenChi2 = requireGoldenChi2;
  LOG(info) << "Track selection, set require golden chi2: " << mRequireGoldenChi2;
}
void TrackSelection::SetMinNClustersTPC(int minNClustersTPC)
{
  mMinNClustersTPC = minNClustersTPC;
  LOG(info) << "Track selection, set min N clusters TPC: " << mMinNClustersTPC;
}
void TrackSelection::SetMinNCrossedRowsTPC(int minNCrossedRowsTPC)
{
  mMinNCrossedRowsTPC = minNCrossedRowsTPC;
  LOG(info) << "Track selection, set min N crossed rows TPC: " << mMinNCrossedRowsTPC;
}
void TrackSelection::SetMinNCrossedRowsOverFindableClustersTPC(float minNCrossedRowsOverFindableClustersTPC)
{
  mMinNCrossedRowsOverFindableClustersTPC = minNCrossedRowsOverFindableClustersTPC;
  LOG(info) << "Track selection, set min N crossed rows over findable clusters TPC: " << mMinNCrossedRowsOverFindableClustersTPC;
}
void TrackSelection::SetMaxTPCFractionSharedCls(float maxTPCFractionSharedCls)
{
  mMaxTPCFractionSharedCls = maxTPCFractionSharedCls;
  LOG(info) << "Track selection, set max fraction of shared clusters TPC: " << mMaxTPCFractionSharedCls;
}
void TrackSelection::SetMinNClustersITS(int minNClustersITS)
{
  mMinNClustersITS = minNClustersITS;
  LOG(info) << "Track selection, set min N clusters ITS: " << mMinNClustersITS;
}
void TrackSelection::SetMaxChi2PerClusterTPC(float maxChi2PerClusterTPC)
{
  mMaxChi2PerClusterTPC = maxChi2PerClusterTPC;
  LOG(info) << "Track selection, set max chi2 per cluster TPC: " << mMaxChi2PerClusterTPC;
}
void TrackSelection::SetMaxChi2PerClusterITS(float maxChi2PerClusterITS)
{
  mMaxChi2PerClusterITS = maxChi2PerClusterITS;
  LOG(info) << "Track selection, set max chi2 per cluster ITS: " << mMaxChi2PerClusterITS;
}
void TrackSelection::SetMaxDcaXY(float maxDcaXY)
{
  mMaxDcaXY = maxDcaXY;
  LOG(info) << "Track selection, set max DCA xy: " << mMaxDcaXY;
}
void TrackSelection::SetMaxDcaZ(float maxDcaZ)
{
  mMaxDcaZ = maxDcaZ;
  LOG(info) << "Track selection, set max DCA z: " << mMaxDcaZ;
}

void TrackSelection::SetMaxDcaXYPtDep(std::function<float(float)> ptDepCut)
{
  mMaxDcaXYPtDep = ptDepCut;
  LOG(info) << "Track selection, set max DCA xy pt dep: " << mMaxDcaXYPtDep(1.0);
}

void TrackSelection::SetRequireHitsInITSLayers(int8_t minNRequiredHits, std::set<uint8_t> requiredLayers)
{
  // layer 0 corresponds to the the innermost ITS layer
  mRequiredITSHits.push_back(std::make_pair(minNRequiredHits, requiredLayers));
  LOG(info) << "Track selection, set require hits in ITS layers: " << static_cast<int>(minNRequiredHits);
}
void TrackSelection::SetRequireNoHitsInITSLayers(std::set<uint8_t> excludedLayers)
{
  mRequiredITSHits.push_back(std::make_pair(-1, excludedLayers));
  LOG(info) << "Track selection, set require no hits in ITS layers";
}

void TrackSelection::print() const
{
  LOG(info) << "Track selection:";
  for (int i = 0; i < static_cast<int>(TrackCuts::kNCuts); i++) {
    switch (static_cast<TrackCuts>(i)) {
      case TrackCuts::kTrackType:
        LOG(info) << mCutNames[i].data() << " == " << static_cast<int>(mTrackType);
        break;
      case TrackCuts::kPtRange:
        LOG(info) << mCutNames[i] << " in [" << mMinPt << ", " << mMaxPt << "]";
        break;
      case TrackCuts::kEtaRange:
        LOG(info) << mCutNames[i] << " in [" << mMinEta << ", " << mMaxEta << "]";
        break;
      case TrackCuts::kTPCNCls:
        LOG(info) << mCutNames[i] << " > " << mMinNClustersTPC;
        break;
      case TrackCuts::kTPCCrossedRows:
        LOG(info) << mCutNames[i] << " > " << mMinNCrossedRowsTPC;
        break;
      case TrackCuts::kTPCCrossedRowsOverNCls:
        LOG(info) << mCutNames[i] << " > " << mMinNCrossedRowsOverFindableClustersTPC;
        break;
      case TrackCuts::kTPCChi2NDF:
        LOG(info) << mCutNames[i] << " < " << mMaxChi2PerClusterTPC;
        break;
      case TrackCuts::kTPCRefit:
        LOG(info) << mCutNames[i] << " == " << mRequireTPCRefit;
        break;
      case TrackCuts::kITSNCls:
        LOG(info) << mCutNames[i] << " > " << mMinNClustersITS;
        break;
      case TrackCuts::kITSChi2NDF:
        LOG(info) << mCutNames[i] << " < " << mMaxChi2PerClusterITS;
        break;
      case TrackCuts::kITSRefit:
        LOG(info) << mCutNames[i] << " == " << mRequireITSRefit;
        break;
      case TrackCuts::kITSHits:
        for (const auto& itsRequirement : mRequiredITSHits) {
          LOG(info) << mCutNames[i] << " == " << itsRequirement.first;
        }
        break;
      case TrackCuts::kGoldenChi2:
        LOG(info) << mCutNames[i] << " == " << mRequireGoldenChi2;
        break;
      case TrackCuts::kDCAxy:
        LOG(info) << mCutNames[i] << " < " << mMaxDcaXY;
        break;
      case TrackCuts::kDCAz:
        LOG(info) << mCutNames[i] << " < " << mMaxDcaZ;
        break;
      case TrackCuts::kTPCFracSharedCls:
        LOG(info) << mCutNames[i] << " < " << mMaxTPCFractionSharedCls;
        break;
      default:
        LOG(fatal) << "Cut unknown!";
    }
  }
}
