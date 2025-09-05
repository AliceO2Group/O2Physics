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

#include "Common/CCDB/TriggerAliases.h"

#include <Framework/Logger.h>

#include <cstdint>
#include <string>

std::string aliasLabels[kNaliases] = {
  "kINT7",
  "kEMC7",
  "kINT7inMUON",
  "kMuonSingleLowPt7",
  "kMuonSingleHighPt7",
  "kMuonUnlikeLowPt7",
  "kMuonLikeLowPt7",
  "kCUP8",
  "kCUP9",
  "kMUP10",
  "kMUP11",
  "kINT1",
  "kUnbiased",
  "kDMC7",
  "kEG1",
  "kEJ1",
  "kEG2",
  "kEJ2",
  "kDG1",
  "kDJ1",
  "kDG2",
  "kDJ2",
  "kTVXinTRD",
  "kTVXinEMC",
  "kTVXinPHOS",
  "kTVXinHMP",
  "kPHOS",
  "kALL"};

void TriggerAliases::AddClassIdToAlias(uint32_t aliasId, int classId)
{
  if (classId < 0 || classId > 99) {
    LOGF(fatal, "Invalid classId = %d for aliasId = %d\n", classId, aliasId);
  } else if (classId < 50) {
    mAliasToTriggerMask[aliasId] |= 1ull << classId;
  } else {
    mAliasToTriggerMaskNext50[aliasId] |= 1ull << (classId - 50);
  }
}

void TriggerAliases::Print()
{
  for (const auto& alias : mAliasToTriggerMask) {
    LOGP(info, "alias={} classMask ={}", aliasLabels[alias.first], alias.second);
  }
}
