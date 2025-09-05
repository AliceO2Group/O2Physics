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

#ifndef COMMON_CCDB_TRIGGERALIASES_H_
#define COMMON_CCDB_TRIGGERALIASES_H_

#include <Rtypes.h>
#include <RtypesCore.h>

#include <cstdint>
#include <map>
#include <string>

enum triggerAliases {
  kINT7 = 0,
  kEMC7,
  kINT7inMUON,
  kMuonSingleLowPt7,
  kMuonSingleHighPt7,
  kMuonUnlikeLowPt7,
  kMuonLikeLowPt7,
  kCUP8,
  kCUP9,
  kMUP10,
  kMUP11,
  kINT1,
  kUnbiased,
  kDMC7,
  kEG1,
  kEJ1,
  kEG2,
  kEJ2,
  kDG1,
  kDJ1,
  kDG2,
  kDJ2,
  kTVXinTRD,
  kTVXinEMC,
  kTVXinPHOS,
  kTVXinHMP,
  kPHOS,
  kALL,
  kNaliases
};

extern std::string aliasLabels[kNaliases];

class TriggerAliases
{
 public:
  TriggerAliases() = default;
  ~TriggerAliases() = default;

  void AddAlias(uint32_t aliasId, std::string const& classNames) { mAliasToClassNames[aliasId] = classNames; }
  void AddClassIdToAlias(uint32_t aliasId, int classId);
  const std::map<uint32_t, std::string>& GetAliasToClassNamesMap() const { return mAliasToClassNames; }
  const std::map<uint32_t, ULong64_t>& GetAliasToTriggerMaskMap() const { return mAliasToTriggerMask; }
  const std::map<uint32_t, ULong64_t>& GetAliasToTriggerMaskNext50Map() const { return mAliasToTriggerMaskNext50; }
  void Print();

 private:
  std::map<uint32_t, std::string> mAliasToClassNames;
  std::map<uint32_t, ULong64_t> mAliasToTriggerMask;
  std::map<uint32_t, ULong64_t> mAliasToTriggerMaskNext50;
  ClassDefNV(TriggerAliases, 8)
};

#endif // COMMON_CCDB_TRIGGERALIASES_H_
