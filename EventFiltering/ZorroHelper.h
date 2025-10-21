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

#ifndef EVENTFILTERING_ZORROHELPER_H_
#define EVENTFILTERING_ZORROHELPER_H_

#include <Rtypes.h>
#include <RtypesCore.h>

struct ZorroHelper {
  ULong64_t bcAOD, bcEvSel, trigMask[2], selMask[2];
  ClassDefNV(ZorroHelper, 1);
};

#endif // EVENTFILTERING_ZORROHELPER_H_
