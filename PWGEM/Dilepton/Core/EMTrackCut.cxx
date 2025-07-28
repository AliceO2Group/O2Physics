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
// Class for track Cut
//

#include "PWGEM/Dilepton/Core/EMTrackCut.h"

#include "Framework/Logger.h"

#include <set>
#include <utility>

ClassImp(EMTrackCut);

const std::pair<int8_t, std::set<uint8_t>> EMTrackCut::its_ib_any_Requirement = {1, {0, 1, 2}}; // hits on any ITS ib layers.
const std::pair<int8_t, std::set<uint8_t>> EMTrackCut::its_ib_1st_Requirement = {1, {0}};       // hit on 1st ITS ib layers.

void EMTrackCut::SetTrackPtRange(float minPt, float maxPt)
{
  mMinTrackPt = minPt;
  mMaxTrackPt = maxPt;
  LOG(info) << "EMTrack Cut, set track pt range: " << mMinTrackPt << " - " << mMaxTrackPt;
}
void EMTrackCut::SetTrackEtaRange(float minEta, float maxEta)
{
  mMinTrackEta = minEta;
  mMaxTrackEta = maxEta;
  LOG(info) << "EMTrack Cut, set track eta range: " << mMinTrackEta << " - " << mMaxTrackEta;
}
void EMTrackCut::SetTrackPhiRange(float minPhi, float maxPhi)
{
  mMinTrackPhi = minPhi;
  mMaxTrackPhi = maxPhi;
  LOG(info) << "EMTrack Cut, set track phi range (rad.): " << mMinTrackPhi << " - " << mMaxTrackPhi;
}
void EMTrackCut::SetMinNClustersTPC(int minNClustersTPC)
{
  mMinNClustersTPC = minNClustersTPC;
  LOG(info) << "EMTrack Cut, set min N clusters TPC: " << mMinNClustersTPC;
}
void EMTrackCut::SetMinNCrossedRowsTPC(int minNCrossedRowsTPC)
{
  mMinNCrossedRowsTPC = minNCrossedRowsTPC;
  LOG(info) << "EMTrack Cut, set min N crossed rows TPC: " << mMinNCrossedRowsTPC;
}
void EMTrackCut::SetMinNCrossedRowsOverFindableClustersTPC(float minNCrossedRowsOverFindableClustersTPC)
{
  mMinNCrossedRowsOverFindableClustersTPC = minNCrossedRowsOverFindableClustersTPC;
  LOG(info) << "EMTrack Cut, set min N crossed rows over findable clusters TPC: " << mMinNCrossedRowsOverFindableClustersTPC;
}
void EMTrackCut::SetMaxFracSharedClustersTPC(float max)
{
  mMaxFracSharedClustersTPC = max;
  LOG(info) << "EMTrack Cut, set max fraction of shared clusters in  TPC: " << mMaxFracSharedClustersTPC;
}
void EMTrackCut::SetChi2PerClusterTPC(float min, float max)
{
  mMinChi2PerClusterTPC = min;
  mMaxChi2PerClusterTPC = max;
  LOG(info) << "EMTrack Cut, set chi2 per cluster TPC range: " << mMinChi2PerClusterTPC << " - " << mMaxChi2PerClusterTPC;
}

void EMTrackCut::SetNClustersITS(int min, int max)
{
  mMinNClustersITS = min;
  mMaxNClustersITS = max;
  LOG(info) << "EMTrack Cut, set N clusters ITS range: " << mMinNClustersITS << " - " << mMaxNClustersITS;
}
void EMTrackCut::SetChi2PerClusterITS(float min, float max)
{
  mMinChi2PerClusterITS = min;
  mMaxChi2PerClusterITS = max;
  LOG(info) << "EMTrack Cut, set chi2 per cluster ITS range: " << mMinChi2PerClusterITS << " - " << mMaxChi2PerClusterITS;
}

void EMTrackCut::SetTrackMaxDcaXY(float maxDcaXY)
{
  mMaxDcaXY = maxDcaXY;
  LOG(info) << "EMTrack Cut, set max DCA xy: " << mMaxDcaXY;
}
void EMTrackCut::SetTrackMaxDcaZ(float maxDcaZ)
{
  mMaxDcaZ = maxDcaZ;
  LOG(info) << "EMTrack Cut, set max DCA z: " << mMaxDcaZ;
}

void EMTrackCut::SetTrackMaxDcaXYPtDep(std::function<float(float)> ptDepCut)
{
  mMaxDcaXYPtDep = ptDepCut;
  LOG(info) << "EMTrack Cut, set max DCA xy pt dep: " << mMaxDcaXYPtDep(1.0);
}

void EMTrackCut::RequireITSibAny(bool flag)
{
  mRequireITSibAny = flag;
  LOG(info) << "EMTrack Cut, require ITS ib any: " << mRequireITSibAny;
}

void EMTrackCut::RequireITSib1st(bool flag)
{
  mRequireITSib1st = flag;
  LOG(info) << "EMTrack Cut, require ITS ib 1st: " << mRequireITSib1st;
}

void EMTrackCut::SetTrackBit(uint16_t bit)
{
  mTrackBit = bit;
  LOG(info) << "EMTrack Cut, require track bits: " << mTrackBit;
}
