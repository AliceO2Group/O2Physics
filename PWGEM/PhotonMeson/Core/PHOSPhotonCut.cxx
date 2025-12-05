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

/// \file PHOSPhotonCut.cxx
/// \brief Source of class for phos photon selection.
/// \author D. Sekihata, daiki.sekihata@cern.ch

#include "PWGEM/PhotonMeson/Core/PHOSPhotonCut.h"

#include <Framework/Logger.h>

#include <Rtypes.h>

ClassImp(PHOSPhotonCut);

const char* PHOSPhotonCut::mCutNames[static_cast<int>(PHOSPhotonCut::PHOSPhotonCuts::kNCuts)] = {"Energy", "Dispersion", "CPV"};

void PHOSPhotonCut::SetEnergyRange(float min, float max)
{
  mMinEnergy = min;
  mMaxEnergy = max;
  LOG(info) << "PHOS Photon selection, set mee range: " << mMinEnergy << " - " << mMaxEnergy;
}

void PHOSPhotonCut::print() const
{
  LOG(info) << "PHOS Photon Cut:";
  for (int i = 0; i < static_cast<int>(PHOSPhotonCuts::kNCuts); i++) {
    switch (static_cast<PHOSPhotonCuts>(i)) {
      case PHOSPhotonCuts::kEnergy:
        LOG(info) << mCutNames[i] << " in [" << mMinEnergy << ", " << mMaxEnergy << "]";
        break;
      default:
        LOG(fatal) << "Cut unknown!";
    }
  }
}
