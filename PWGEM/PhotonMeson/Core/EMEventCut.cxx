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
#include "PWGEM/PhotonMeson/Core/EMEventCut.h"

ClassImp(EMEventCut);

const char* EMEventCut::mCutNames[static_cast<int>(EMEventCut::EMEventCuts::kNCuts)] = {"Sel8", "FT0AND", "Zvtx", "eNoTFB", "RequireNoITSROFB", "NoSameBunchPileup", "GoodVertexITSTPC", "GoodZvtxFT0vsPV"};

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

void EMEventCut::SetRequireIsGoodZvtxFT0vsPV(bool flag)
{
  mRequireGoodZvtxFT0vsPV = flag;
  LOG(info) << "EM Event Cut, require good Zvtx between FT0 vs. PV: " << mRequireGoodZvtxFT0vsPV;
}

void EMEventCut::print() const
{
  LOG(info) << "EM Event Cut:";
  for (int i = 0; i < static_cast<int>(EMEventCuts::kNCuts); i++) {
    switch (static_cast<EMEventCuts>(i)) {
      case EMEventCuts::kFT0AND:
        LOG(info) << mCutNames[i] << " = " << mRequireFT0AND;
        break;
      case EMEventCuts::kZvtx:
        LOG(info) << mCutNames[i] << " in [" << mMinZvtx << ", " << mMaxZvtx << "]";
        break;
      case EMEventCuts::kNoTFB:
        LOG(info) << mCutNames[i] << " = " << mRequireNoTFB;
        break;
      case EMEventCuts::kNoITSROFB:
        LOG(info) << mCutNames[i] << " = " << mRequireNoITSROFB;
        break;

      default:
        LOG(fatal) << "Cut unknown!";
    }
  }
}
