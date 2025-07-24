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

#include "Framework/Logger.h"
#include "PWGEM/PhotonMeson/Core/PairCut.h"

ClassImp(PairCut);

const char* PairCut::mCutNames[static_cast<int>(PairCut::PairCuts::kNCuts)] = {"Asym"};

void PairCut::SetAsymRange(float min, float max)
{
  mMinAsym = min;
  mMaxAsym = max;
  LOG(info) << "Pair Cut, set energy asymmetry range: " << mMinAsym << " - " << mMaxAsym;
}

void PairCut::print() const
{
  LOG(info) << "Pair Cut:";
  for (int i = 0; i < static_cast<int>(PairCuts::kNCuts); i++) {
    switch (static_cast<PairCuts>(i)) {
      case PairCuts::kAsym:
        LOG(info) << mCutNames[i] << " in [" << mMinAsym << ", " << mMaxAsym << "]";
        break;
      default:
        LOG(fatal) << "Cut unknown!";
    }
  }
}
