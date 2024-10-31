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
