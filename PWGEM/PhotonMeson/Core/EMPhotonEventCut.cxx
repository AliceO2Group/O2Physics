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
// Class for em photon event selection
//

#include "Framework/Logger.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"

ClassImp(EMPhotonEventCut);

void EMPhotonEventCut::SetRequireSel8(bool flag)
{
  mRequireSel8 = flag;
  LOG(info) << "EM Event Cut, require sel8: " << mRequireSel8;
}

void EMPhotonEventCut::SetRequireFT0AND(bool flag)
{
  mRequireFT0AND = flag;
  LOG(info) << "EM Event Cut, require FT0AND: " << mRequireFT0AND;
}

void EMPhotonEventCut::SetZvtxRange(float min, float max)
{
  mMinZvtx = min;
  mMaxZvtx = max;
  LOG(info) << "EM Event Cut, set z vtx range: " << mMinZvtx << " - " << mMaxZvtx;
}

void EMPhotonEventCut::SetRequireNoTFB(bool flag)
{
  mRequireNoTFB = flag;
  LOG(info) << "EM Event Cut, require No TF border: " << mRequireNoTFB;
}

void EMPhotonEventCut::SetRequireNoITSROFB(bool flag)
{
  mRequireNoITSROFB = flag;
  LOG(info) << "EM Event Cut, require No ITS ROF border: " << mRequireNoITSROFB;
}

void EMPhotonEventCut::SetRequireNoSameBunchPileup(bool flag)
{
  mRequireNoSameBunchPileup = flag;
  LOG(info) << "EM Event Cut, require No same bunch pileup: " << mRequireNoSameBunchPileup;
}

void EMPhotonEventCut::SetRequireVertexITSTPC(bool flag)
{
  mRequireVertexITSTPC = flag;
  LOG(info) << "EM Event Cut, require vertex reconstructed by ITS-TPC matched track: " << mRequireVertexITSTPC;
}

void EMPhotonEventCut::SetRequireVertexTOFmatched(bool flag)
{
  mRequireVertexTOFmatched = flag;
  LOG(info) << "EM Event Cut, require vertex reconstructed by ITS-TPC-TOF matched track: " << mRequireVertexTOFmatched;
}

void EMPhotonEventCut::SetRequireGoodZvtxFT0vsPV(bool flag)
{
  mRequireGoodZvtxFT0vsPV = flag;
  LOG(info) << "EM Event Cut, require good Zvtx between FT0 vs. PV: " << mRequireGoodZvtxFT0vsPV;
}

void EMPhotonEventCut::SetRequireNoCollInTimeRangeStandard(bool flag)
{
  mRequireNoCollInTimeRangeStandard = flag;
  LOG(info) << "EM Event Cut, require No collision in time range standard: " << mRequireNoCollInTimeRangeStandard;
}

void EMPhotonEventCut::SetRequireNoCollInTimeRangeStrict(bool flag)
{
  mRequireNoCollInTimeRangeStrict = flag;
  LOG(info) << "EM Event Cut, require No collision in time range strict: " << mRequireNoCollInTimeRangeStrict;
}
void EMPhotonEventCut::SetRequireNoCollInITSROFStandard(bool flag)
{
  mRequireNoCollInITSROFStandard = flag;
  LOG(info) << "EM Event Cut, require No collision in ITS TOF standard: " << mRequireNoCollInITSROFStandard;
}

void EMPhotonEventCut::SetRequireNoCollInITSROFStrict(bool flag)
{
  mRequireNoCollInITSROFStrict = flag;
  LOG(info) << "EM Event Cut, require No collision in ITS ROF strict: " << mRequireNoCollInITSROFStrict;
}

void EMPhotonEventCut::SetRequireNoHighMultCollInPrevRof(bool flag)
{
  mRequireNoHighMultCollInPrevRof = flag;
  LOG(info) << "EM Event Cut, require No HM collision in previous ITS ROF: " << mRequireNoHighMultCollInPrevRof;
}

void EMPhotonEventCut::SetRequireGoodITSLayer3(bool flag)
{
  mRequireGoodITSLayer3 = flag;
  LOG(info) << "EM Event Cut, require GoodITSLayer3: " << mRequireGoodITSLayer3;
}

void EMPhotonEventCut::SetRequireGoodITSLayer0123(bool flag)
{
  mRequireGoodITSLayer0123 = flag;
  LOG(info) << "EM Event Cut, require GoodITSLayer0123: " << mRequireGoodITSLayer0123;
}

void EMPhotonEventCut::SetRequireGoodITSLayersAll(bool flag)
{
  mRequireGoodITSLayersAll = flag;
  LOG(info) << "EM Event Cut, require GoodITSLayersAll: " << mRequireGoodITSLayersAll;
}

void EMPhotonEventCut::SetRequireEMCReadoutInMB(bool flag)
{
  mRequireEMCReadoutInMB = flag;
  LOG(info) << "EM Photon Event Cut, require the EMC to be read out in an MB collision by checking kTVXinEMC: " << mRequireEMCReadoutInMB;
}

void EMPhotonEventCut::SetRequireEMCHardwareTriggered(bool flag)
{
  mRequireEMCHardwareTriggered = flag;
  LOG(info) << "EM Photon Event Cut, require the EMC to be triggered by requiring kEMC7 or kDMC7: " << mRequireEMCHardwareTriggered;
}
