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

#include <Framework/Logger.h>

#include <Rtypes.h>

#include <cstdint>

ClassImp(EMTrackCut);

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

// void EMTrackCut::SetTrackMaxDcaXYPtDep(std::function<float(float)> ptDepCut)
// {
//   mMaxDcaXYPtDep = ptDepCut;
//   LOG(info) << "EMTrack Cut, set max DCA xy pt dep: " << mMaxDcaXYPtDep(1.0);
// }

void EMTrackCut::SetTrackBit(uint16_t bit)
{
  mTrackBit = bit;
  LOG(info) << "EMTrack Cut, require track bits: " << mTrackBit;
}
