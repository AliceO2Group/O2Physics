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
// Class for em event selection
//

#include "Framework/Logger.h"
#include "PWGEM/Dilepton/Core/EMEventCut.h"

ClassImp(EMEventCut);

void EMEventCut::SetRequireSel8(bool flag)
{
  mRequireSel8 = flag;
  LOG(info) << "EM Event Cut, require sel8: " << mRequireSel8;
}

void EMEventCut::SetRequireFT0AND(bool flag)
{
  mRequireFT0AND = flag;
  LOG(info) << "EM Event Cut, require FT0AND: " << mRequireFT0AND;
}

void EMEventCut::SetZvtxRange(float min, float max)
{
  mMinZvtx = min;
  mMaxZvtx = max;
  LOG(info) << "EM Event Cut, set z vtx range: " << mMinZvtx << " - " << mMaxZvtx;
}

void EMEventCut::SetOccupancyRange(int min, int max)
{
  mMinOccupancy = min;
  mMaxOccupancy = max;
  LOG(info) << "EM Event Cut, set occupancy range: " << mMinOccupancy << " - " << mMaxOccupancy;
}

void EMEventCut::SetRequireNoTFB(bool flag)
{
  mRequireNoTFB = flag;
  LOG(info) << "EM Event Cut, require No TF border: " << mRequireNoTFB;
}

void EMEventCut::SetRequireNoITSROFB(bool flag)
{
  mRequireNoITSROFB = flag;
  LOG(info) << "EM Event Cut, require No ITS ROF border: " << mRequireNoITSROFB;
}

void EMEventCut::SetRequireNoSameBunchPileup(bool flag)
{
  mRequireNoSameBunchPileup = flag;
  LOG(info) << "EM Event Cut, require No same bunch pileup: " << mRequireNoSameBunchPileup;
}

void EMEventCut::SetRequireVertexITSTPC(bool flag)
{
  mRequireVertexITSTPC = flag;
  LOG(info) << "EM Event Cut, require vertex reconstructed by ITS-TPC matched track: " << mRequireVertexITSTPC;
}

void EMEventCut::SetRequireGoodZvtxFT0vsPV(bool flag)
{
  mRequireGoodZvtxFT0vsPV = flag;
  LOG(info) << "EM Event Cut, require good Zvtx between FT0 vs. PV: " << mRequireGoodZvtxFT0vsPV;
}

void EMEventCut::SetRequireNoCollInTimeRangeStandard(bool flag)
{
  mRequireNoCollInTimeRangeStandard = flag;
  LOG(info) << "EM Event Cut, require No collision in time range standard: " << mRequireNoCollInTimeRangeStandard;
}
void EMEventCut::SetRequireNoCollInTimeRangeNarrow(bool flag)
{
  mRequireNoCollInTimeRangeNarrow = flag;
  LOG(info) << "EM Event Cut, require No collision in time range narrow: " << mRequireNoCollInTimeRangeNarrow;
}
