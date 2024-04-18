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

const char* EMEventCut::mCutNames[static_cast<int>(EMEventCut::EMEventCuts::kNCuts)] = {"RequireFT0AND", "Zvtx", "RequireNoTFB", "RequireNoITSROFB"};

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
