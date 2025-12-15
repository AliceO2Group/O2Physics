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

/// \file JetDerivedDataUtilities.h
/// \brief Jet derived data related utilities
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_CORE_JETDERIVEDDATAUTILITIES_H_
#define PWGJE_CORE_JETDERIVEDDATAUTILITIES_H_

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/CCDB/TriggerAliases.h"

#include <Rtypes.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

namespace jetderiveddatautilities
{

static constexpr float mPion = 0.139; // TDatabasePDG::Instance()->GetParticle(211)->Mass(); //can be removed when pion mass becomes default for unidentified tracks

enum JCollisionSel {
  sel8 = 0,
  sel7 = 1,
  selKINT7 = 2,
  selTVX = 3,
  selNoTimeFrameBorder = 4,
  selNoITSROFrameBorder = 5,
  selNoSameBunchPileup = 6,
  selIsGoodZvtxFT0vsPV = 7,
  selNoCollInTimeRangeStandard = 8,
  selNoCollInRofStandard = 9
};

enum JCollisionSubGeneratorId {
  none = -1,
  mbGap = 0
};

template <typename T>
bool commonCollisionSelection(T const& collision, bool skipMBGapEvents = true, bool rctSelection = true, std::string rctLabel = "CBT_hadronPID", bool rejectLimitedAcceptanceRct = false, bool requireZDCRct = false)
{
  if (skipMBGapEvents && collision.getSubGeneratorId() == JCollisionSubGeneratorId::mbGap) {
    return false;
  }
  o2::aod::rctsel::RCTFlagsChecker rctChecker;
  rctChecker.init(rctLabel, requireZDCRct, rejectLimitedAcceptanceRct);
  if (rctSelection && !rctChecker.checkTable(collision)) { // CBT_hadronPID given as default so that TOF is included in RCT selection to benefit from better timing for tracks. Impact of this for inclusive jets should be studied
    return false;
  }
  return true;
}

template <typename T>
bool selectMcCollision(T const& mcCollision, bool skipMBGapEvents = true, bool rctSelection = true, std::string rctLabel = "CBT_hadronPID", bool rejectLimitedAcceptanceRct = false, bool requireZDCRct = false)
{
  return commonCollisionSelection(mcCollision, skipMBGapEvents, rctSelection, rctLabel, rejectLimitedAcceptanceRct, requireZDCRct);
}

template <typename T>
bool selectCollision(T const& collision, const std::vector<int>& eventSelectionMaskBits, bool skipMBGapEvents = true, bool rctSelection = true, std::string rctLabel = "CBT_hadronPID", bool rejectLimitedAcceptanceRct = false, bool requireZDCRct = false)
{

  if (!commonCollisionSelection(collision, skipMBGapEvents, rctSelection, rctLabel, rejectLimitedAcceptanceRct, requireZDCRct)) {
    return false;
  }
  if (eventSelectionMaskBits.size() == 0) {
    return true;
  }
  for (auto eventSelectionMaskBit : eventSelectionMaskBits) {
    if (!(collision.eventSel() & (1 << eventSelectionMaskBit))) {
      return false;
    }
  }
  return true;
}

bool eventSelectionMasksContainSelection(const std::string& eventSelectionMasks, std::string selection)
{
  size_t position = 0;
  while ((position = eventSelectionMasks.find(selection, position)) != std::string::npos) {
    bool validStart = (position == 0 || eventSelectionMasks[position - 1] == '+');
    bool validEnd = (position + selection.length() == eventSelectionMasks.length() || eventSelectionMasks[position + selection.length()] == '+');
    if (validStart && validEnd) {
      return true;
    }
    position += selection.length();
  }
  return false;
}

std::vector<int> initialiseEventSelectionBits(const std::string& eventSelectionMasks)
{
  std::vector<int> eventSelectionMaskBits;
  if (eventSelectionMasksContainSelection(eventSelectionMasks, "sel8")) {
    eventSelectionMaskBits.push_back(JCollisionSel::sel8);
  }
  if (eventSelectionMasksContainSelection(eventSelectionMasks, "sel7")) {
    eventSelectionMaskBits.push_back(JCollisionSel::sel7);
  }
  if (eventSelectionMasksContainSelection(eventSelectionMasks, "selKINT7")) {
    eventSelectionMaskBits.push_back(JCollisionSel::selKINT7);
  }
  if (eventSelectionMasksContainSelection(eventSelectionMasks, "TVX")) {
    eventSelectionMaskBits.push_back(JCollisionSel::selTVX);
  }
  if (eventSelectionMasksContainSelection(eventSelectionMasks, "NoTimeFrameBorder")) {
    eventSelectionMaskBits.push_back(JCollisionSel::selNoTimeFrameBorder);
  }
  if (eventSelectionMasksContainSelection(eventSelectionMasks, "NoITSROFrameBorder")) {
    eventSelectionMaskBits.push_back(JCollisionSel::selNoITSROFrameBorder);
  }
  if (eventSelectionMasksContainSelection(eventSelectionMasks, "NoSameBunchPileup")) {
    eventSelectionMaskBits.push_back(JCollisionSel::selNoSameBunchPileup);
  }
  if (eventSelectionMasksContainSelection(eventSelectionMasks, "IsGoodZvtxFT0vsPV")) {
    eventSelectionMaskBits.push_back(JCollisionSel::selIsGoodZvtxFT0vsPV);
  }
  if (eventSelectionMasksContainSelection(eventSelectionMasks, "NoCollInTimeRangeStandard")) {
    eventSelectionMaskBits.push_back(JCollisionSel::selNoCollInTimeRangeStandard);
  }
  if (eventSelectionMasksContainSelection(eventSelectionMasks, "NoCollInRofStandard")) {
    eventSelectionMaskBits.push_back(JCollisionSel::selNoCollInRofStandard);
  }

  if (eventSelectionMasksContainSelection(eventSelectionMasks, "sel8Full")) {
    eventSelectionMaskBits.push_back(JCollisionSel::sel8);
    eventSelectionMaskBits.push_back(JCollisionSel::selNoSameBunchPileup);
  }
  if (eventSelectionMasksContainSelection(eventSelectionMasks, "sel8FullPbPb")) {
    eventSelectionMaskBits.push_back(JCollisionSel::sel8);
    eventSelectionMaskBits.push_back(JCollisionSel::selNoCollInTimeRangeStandard);
    eventSelectionMaskBits.push_back(JCollisionSel::selNoCollInRofStandard);
  }
  if (eventSelectionMasksContainSelection(eventSelectionMasks, "selUnanchoredMC")) {
    eventSelectionMaskBits.push_back(JCollisionSel::selTVX);
  }
  if (eventSelectionMasksContainSelection(eventSelectionMasks, "selMC")) {
    eventSelectionMaskBits.push_back(JCollisionSel::selTVX);
    eventSelectionMaskBits.push_back(JCollisionSel::selNoTimeFrameBorder);
  }
  if (eventSelectionMasksContainSelection(eventSelectionMasks, "selMCFull")) {
    eventSelectionMaskBits.push_back(JCollisionSel::selTVX);
    eventSelectionMaskBits.push_back(JCollisionSel::selNoTimeFrameBorder);
    eventSelectionMaskBits.push_back(JCollisionSel::selNoSameBunchPileup);
  }
  if (eventSelectionMasksContainSelection(eventSelectionMasks, "selMCFullPbPb")) {
    eventSelectionMaskBits.push_back(JCollisionSel::selTVX);
    eventSelectionMaskBits.push_back(JCollisionSel::selNoCollInTimeRangeStandard);
    eventSelectionMaskBits.push_back(JCollisionSel::selNoCollInRofStandard);
  }
  if (eventSelectionMasksContainSelection(eventSelectionMasks, "sel7KINT7")) {
    eventSelectionMaskBits.push_back(JCollisionSel::sel7);
    eventSelectionMaskBits.push_back(JCollisionSel::selKINT7);
  }
  return eventSelectionMaskBits;
}

template <typename T>
uint16_t setEventSelectionBit(T const& collision)
{
  uint16_t bit = 0;
  if (collision.sel8()) {
    SETBIT(bit, JCollisionSel::sel8);
  }
  if (collision.sel7()) {
    SETBIT(bit, JCollisionSel::sel7);
  }
  if (collision.alias_bit(kINT7)) {
    SETBIT(bit, JCollisionSel::selKINT7);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
    SETBIT(bit, JCollisionSel::selTVX);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
    SETBIT(bit, JCollisionSel::selNoTimeFrameBorder);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
    SETBIT(bit, JCollisionSel::selNoITSROFrameBorder);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
    SETBIT(bit, JCollisionSel::selNoSameBunchPileup);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
    SETBIT(bit, JCollisionSel::selIsGoodZvtxFT0vsPV);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
    SETBIT(bit, JCollisionSel::selNoCollInTimeRangeStandard);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
    SETBIT(bit, JCollisionSel::selNoCollInRofStandard);
  }
  return bit;
}

template <typename T>
bool eventEMCAL(T const& collision)
{
  // Check for EMCAL in readout requires any of the EMCAL trigger classes (including EMC in MB trigger) to fire
  std::array<triggerAliases, 11> selectAliases = {{triggerAliases::kTVXinEMC, triggerAliases::kEMC7, triggerAliases::kDMC7, triggerAliases::kEG1, triggerAliases::kEG2, triggerAliases::kDG1, triggerAliases::kDG2, triggerAliases::kEJ1, triggerAliases::kEJ2, triggerAliases::kDJ1, triggerAliases::kDJ2}};
  bool found = false;
  for (auto alias : selectAliases) {
    if (collision.alias_bit(alias)) {
      found = true;
      break;
    }
  }
  return found;
}

inline const std::string JTriggerMasks = "fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL";

enum JTrigSel {
  noTrigSel = 0,
  JetChLowPt = 1,
  JetChHighPt = 2,
  TrackLowPt = 3,
  TrackHighPt = 4,
  JetD0ChLowPt = 5,
  JetD0ChHighPt = 6,
  JetLcChLowPt = 7,
  JetLcChHighPt = 8,
  EMCALReadout = 9,
  JetFullHighPt = 10,
  JetFullLowPt = 11,
  JetNeutralHighPt = 12,
  JetNeutralLowPt = 13,
  GammaVeryHighPtEMCAL = 14,
  GammaVeryHighPtDCAL = 15,
  GammaHighPtEMCAL = 16,
  GammaHighPtDCAL = 17,
  GammaLowPtEMCAL = 18,
  GammaLowPtDCAL = 19,
  GammaVeryLowPtEMCAL = 20,
  GammaVeryLowPtDCAL = 21
};

template <typename T>
bool selectTrigger(T const& collision, const std::vector<int>& triggerMaskBits)
{
  if (triggerMaskBits.size() == 0) {
    return true;
  }
  for (auto triggerMaskBit : triggerMaskBits) {
    if (collision.triggerSel() & (1 << triggerMaskBit)) {
      return true;
    }
  }
  return false;
}

template <typename T>
bool selectTrigger(T const& collision, int triggerMaskBit)
{
  if (triggerMaskBit == -1) {
    return false;
  }
  return collision.triggerSel() & (1 << triggerMaskBit);
}

bool triggerMasksContainTrigger(const std::string& triggerMasks, std::string trigger)
{
  size_t position = 0;
  while ((position = triggerMasks.find(trigger, position)) != std::string::npos) {
    bool validStart = (position == 0 || triggerMasks[position - 1] == ',');
    bool validEnd = (position + trigger.length() == triggerMasks.length() || triggerMasks[position + trigger.length()] == ',');
    if (validStart && validEnd) {
      return true;
    }
    position += trigger.length();
  }
  return false;
}

std::vector<int> initialiseTriggerMaskBits(const std::string& triggerMasks)
{
  std::vector<int> triggerMaskBits;
  if (triggerMasksContainTrigger(triggerMasks, "fJetChLowPt")) {
    triggerMaskBits.push_back(JTrigSel::JetChLowPt);
  }
  if (triggerMasksContainTrigger(triggerMasks, "fJetChHighPt")) {
    triggerMaskBits.push_back(JTrigSel::JetChHighPt);
  }
  if (triggerMasksContainTrigger(triggerMasks, "fTrackLowPt")) {
    triggerMaskBits.push_back(JTrigSel::TrackLowPt);
  }
  if (triggerMasksContainTrigger(triggerMasks, "fTrackHighPt")) {
    triggerMaskBits.push_back(JTrigSel::TrackHighPt);
  }
  if (triggerMasksContainTrigger(triggerMasks, "fJetD0ChLowPt")) {
    triggerMaskBits.push_back(JTrigSel::JetD0ChLowPt);
  }
  if (triggerMasksContainTrigger(triggerMasks, "fJetD0ChHighPt")) {
    triggerMaskBits.push_back(JTrigSel::JetD0ChHighPt);
  }
  if (triggerMasksContainTrigger(triggerMasks, "fJetLcChLowPt")) {
    triggerMaskBits.push_back(JTrigSel::JetLcChLowPt);
  }
  if (triggerMasksContainTrigger(triggerMasks, "fJetLcChHighPt")) {
    triggerMaskBits.push_back(JTrigSel::JetLcChHighPt);
  }
  if (triggerMasksContainTrigger(triggerMasks, "fEMCALReadout")) {
    triggerMaskBits.push_back(JTrigSel::EMCALReadout);
  }
  if (triggerMasksContainTrigger(triggerMasks, "fJetFullHighPt")) {
    triggerMaskBits.push_back(JTrigSel::JetFullHighPt);
  }
  if (triggerMasksContainTrigger(triggerMasks, "fJetFullLowPt")) {
    triggerMaskBits.push_back(JTrigSel::JetFullLowPt);
  }
  if (triggerMasksContainTrigger(triggerMasks, "fJetNeutralHighPt")) {
    triggerMaskBits.push_back(JTrigSel::JetNeutralHighPt);
  }
  if (triggerMasksContainTrigger(triggerMasks, "fJetNeutralLowPt")) {
    triggerMaskBits.push_back(JTrigSel::JetNeutralLowPt);
  }
  if (triggerMasksContainTrigger(triggerMasks, "fGammaVeryHighPtEMCAL")) {
    triggerMaskBits.push_back(JTrigSel::GammaVeryHighPtEMCAL);
  }
  if (triggerMasksContainTrigger(triggerMasks, "fGammaVeryHighPtDCAL")) {
    triggerMaskBits.push_back(JTrigSel::GammaVeryHighPtDCAL);
  }
  if (triggerMasksContainTrigger(triggerMasks, "fGammaHighPtEMCAL")) {
    triggerMaskBits.push_back(JTrigSel::GammaHighPtEMCAL);
  }
  if (triggerMasksContainTrigger(triggerMasks, "fGammaHighPtDCAL")) {
    triggerMaskBits.push_back(JTrigSel::GammaHighPtDCAL);
  }
  if (triggerMasksContainTrigger(triggerMasks, "fGammaLowPtEMCAL")) {
    triggerMaskBits.push_back(JTrigSel::GammaLowPtEMCAL);
  }
  if (triggerMasksContainTrigger(triggerMasks, "fGammaLowPtDCAL")) {
    triggerMaskBits.push_back(JTrigSel::GammaLowPtDCAL);
  }
  if (triggerMasksContainTrigger(triggerMasks, "fGammaVeryLowPtEMCAL")) {
    triggerMaskBits.push_back(JTrigSel::GammaVeryLowPtEMCAL);
  }
  if (triggerMasksContainTrigger(triggerMasks, "fGammaVeryLowPtDCAL")) {
    triggerMaskBits.push_back(JTrigSel::GammaVeryLowPtDCAL);
  }
  return triggerMaskBits;
}

uint64_t setTriggerSelectionBit(const std::vector<bool>& triggerDecisions)
{
  uint64_t bit = 0;
  for (std::vector<bool>::size_type i = 0; i < triggerDecisions.size(); i++) {
    if (triggerDecisions[i]) {
      SETBIT(bit, i + 1);
    }
  }
  return bit;
}

enum JTrigSelCh {
  noChargedTigger = 0,
  jetChLowPt = 1,
  jetChHighPt = 2,
  trackLowPt = 3,
  trackHighPt = 4
};

template <typename T>
bool selectChargedTrigger(T const& collision, int triggerSelection)
{
  if (triggerSelection == -1) {
    return true;
  }
  return (collision.chargedTriggerSel() & (1 << triggerSelection));
}

int initialiseChargedTriggerSelection(const std::string& triggerSelection)
{
  if (triggerSelection == "jetChLowPt") {
    return JTrigSelCh::jetChLowPt;
  }
  if (triggerSelection == "jetChHighPt") {
    return JTrigSelCh::jetChHighPt;
  }
  if (triggerSelection == "trackLowPt") {
    return JTrigSelCh::trackLowPt;
  }
  if (triggerSelection == "trackHighPt") {
    return JTrigSelCh::trackHighPt;
  }

  return -1;
}

template <typename T>
uint8_t setChargedTriggerSelectionBit(T const& collision)
{

  uint8_t bit = 0;
  if (collision.hasJetChLowPt()) {
    SETBIT(bit, JTrigSelCh::jetChLowPt);
  }
  if (collision.hasJetChHighPt()) {
    SETBIT(bit, JTrigSelCh::jetChHighPt);
  }
  if (collision.hasTrackLowPt()) {
    SETBIT(bit, JTrigSelCh::trackLowPt);
  }
  if (collision.hasTrackHighPt()) {
    SETBIT(bit, JTrigSelCh::trackHighPt);
  }

  return bit;
}

enum JTrigSelFull {
  noFullTrigger = 0,
  emcalReadout = 1,
  jetFullHighPt = 2,
  jetFullLowPt = 3,
  jetNeutralHighPt = 4,
  jetNeutralLowPt = 5,
  gammaVeryHighPtEMCAL = 6,
  gammaVeryHighPtDCAL = 7,
  gammaHighPtEMCAL = 8,
  gammaHighPtDCAL = 9,
  gammaLowPtEMCAL = 10,
  gammaLowPtDCAL = 11,
  gammaVeryLowPtEMCAL = 12,
  gammaVeryLowPtDCAL = 13
};

template <typename T>
bool selectFullTrigger(T const& collision, int triggerSelection)
{
  if (triggerSelection == -1) {
    return true;
  }
  return (collision.fullTriggerSel() & (1 << triggerSelection));
}

int initialiseFullTriggerSelection(const std::string& triggerSelection)
{
  if (triggerSelection == "emcalReadout") {
    return JTrigSelFull::emcalReadout;
  } else if (triggerSelection == "jetFullHighPt") {
    return JTrigSelFull::jetFullHighPt;
  } else if (triggerSelection == "jetFullLowPt") {
    return JTrigSelFull::jetFullLowPt;
  } else if (triggerSelection == "jetNeutralHighPt") {
    return JTrigSelFull::jetNeutralHighPt;
  } else if (triggerSelection == "jetNeutralLowPt") {
    return JTrigSelFull::jetNeutralLowPt;
  } else if (triggerSelection == "gammaVeryHighPtEMCAL") {
    return JTrigSelFull::gammaVeryHighPtEMCAL;
  } else if (triggerSelection == "gammaVeryHighPtDCAL") {
    return JTrigSelFull::gammaVeryHighPtDCAL;
  } else if (triggerSelection == "gammaHighPtEMCAL") {
    return JTrigSelFull::gammaHighPtEMCAL;
  } else if (triggerSelection == "gammaHighPtDCAL") {
    return JTrigSelFull::gammaHighPtDCAL;
  } else if (triggerSelection == "gammaLowPtEMCAL") {
    return JTrigSelFull::gammaLowPtEMCAL;
  } else if (triggerSelection == "gammaLowPtDCAL") {
    return JTrigSelFull::gammaLowPtDCAL;
  } else if (triggerSelection == "gammaVeryLowPtEMCAL") {
    return JTrigSelFull::gammaVeryLowPtEMCAL;
  } else if (triggerSelection == "gammaVeryLowPtDCAL") {
    return JTrigSelFull::gammaVeryLowPtDCAL;
  }
  return -1;
}

template <typename T>
uint32_t setFullTriggerSelectionBit(T const& collision)
{
  uint32_t bit = 0;
  if (collision.hasEMCALinReadout()) {
    SETBIT(bit, JTrigSelFull::emcalReadout);
  }
  if (collision.hasJetFullHighPt()) {
    SETBIT(bit, JTrigSelFull::jetFullHighPt);
  }
  if (collision.hasJetFullLowPt()) {
    SETBIT(bit, JTrigSelFull::jetFullLowPt);
  }
  if (collision.hasJetNeutralHighPt()) {
    SETBIT(bit, JTrigSelFull::jetNeutralHighPt);
  }
  if (collision.hasJetNeutralLowPt()) {
    SETBIT(bit, JTrigSelFull::jetNeutralLowPt);
  }
  if (collision.hasGammaVeryHighPtEMCAL()) {
    SETBIT(bit, JTrigSelFull::gammaVeryHighPtEMCAL);
  }
  if (collision.hasGammaVeryHighPtDCAL()) {
    SETBIT(bit, JTrigSelFull::gammaVeryHighPtDCAL);
  }
  if (collision.hasGammaHighPtEMCAL()) {
    SETBIT(bit, JTrigSelFull::gammaHighPtEMCAL);
  }
  if (collision.hasGammaHighPtDCAL()) {
    SETBIT(bit, JTrigSelFull::gammaHighPtDCAL);
  }
  if (collision.hasGammaLowPtEMCAL()) {
    SETBIT(bit, JTrigSelFull::gammaLowPtEMCAL);
  }
  if (collision.hasGammaLowPtDCAL()) {
    SETBIT(bit, JTrigSelFull::gammaLowPtDCAL);
  }
  if (collision.hasGammaVeryLowPtEMCAL()) {
    SETBIT(bit, JTrigSelFull::gammaVeryLowPtEMCAL);
  }
  if (collision.hasGammaVeryLowPtDCAL()) {
    SETBIT(bit, JTrigSelFull::gammaVeryLowPtDCAL);
  }
  return bit;
}

enum JTrigSelChHF {
  noChargedHFTigger = 0,
  jetD0ChLowPt = 1,
  jetD0ChHighPt = 2,
  jetLcChLowPt = 3,
  jetLcChHighPt = 4
};

template <typename T>
bool selectChargedHFTrigger(T const& collision, int triggerSelection)
{
  if (triggerSelection == -1) {
    return true;
  }
  return (collision.chargedHFTriggerSel() & (1 << triggerSelection));
}

int initialiseChargedHFTriggerSelection(const std::string& triggerSelection)
{
  if (triggerSelection == "jetD0ChLowPt") {
    return JTrigSelChHF::jetD0ChLowPt;
  }
  if (triggerSelection == "jetD0ChHighPt") {
    return JTrigSelChHF::jetD0ChHighPt;
  }
  if (triggerSelection == "jetLcChLowPt") {
    return JTrigSelChHF::jetLcChLowPt;
  }
  if (triggerSelection == "jetLcChHighPt") {
    return JTrigSelChHF::jetLcChHighPt;
  }
  return -1;
}

template <typename T>
uint8_t setChargedHFTriggerSelectionBit(T const& collision)
{

  uint8_t bit = 0;
  if (collision.hasJetD0ChLowPt()) {
    SETBIT(bit, JTrigSelChHF::jetD0ChLowPt);
  }
  if (collision.hasJetD0ChHighPt()) {
    SETBIT(bit, JTrigSelChHF::jetD0ChHighPt);
  }
  if (collision.hasJetLcChLowPt()) {
    SETBIT(bit, JTrigSelChHF::jetLcChLowPt);
  }
  if (collision.hasJetLcChHighPt()) {
    SETBIT(bit, JTrigSelChHF::jetLcChHighPt);
  }
  return bit;
}

enum JTrackSel {
  trackSign = 0, // warning : this number is hardcoded in the sign coloumn in the JTracks table so should not be changed without changing it there too
  globalTrack = 1,
  qualityTrack = 2,
  qualityTrackWDCA = 3,
  hybridTrack = 4
};

template <typename T>
bool applyTrackKinematics(T const& track, float pTMin = 0.15, float pTMax = 100., float EtaMin = -0.9, float EtaMax = 0.9, float PhiMin = -99., float PhiMax = 99.)
{
  if (track.pt() < pTMin || track.pt() > pTMax || track.eta() < EtaMin || track.eta() > EtaMax || track.phi() < PhiMin || track.phi() > PhiMax) {
    return false;
  }
  return true;
}

template <typename T>
bool selectTrack(T const& track, int trackSelection)
{
  if (trackSelection == -1) {
    return true;
  }
  return (track.trackSel() & (1 << trackSelection));
}

int initialiseTrackSelection(const std::string& trackSelection)
{
  if (trackSelection == "globalTracks") {
    return JTrackSel::globalTrack;
  } else if (trackSelection == "QualityTracks") {
    return JTrackSel::qualityTrack;
  } else if (trackSelection == "QualityTracksWDCA") {
    return JTrackSel::qualityTrackWDCA;
  } else if (trackSelection == "hybridTracks") {
    return JTrackSel::hybridTrack;
  }
  return -1;
}

template <typename T>
uint8_t setTrackSelectionBit(T const& track, float trackDCAZ, float maxDCAZ)
{

  uint8_t bit = 0;

  if (track.sign() == 1) {
    SETBIT(bit, JTrackSel::trackSign);
  }
  if (track.isGlobalTrackWoPtEta()) {
    SETBIT(bit, JTrackSel::globalTrack);
  }
  if (track.isQualityTrack()) {
    SETBIT(bit, JTrackSel::qualityTrack);
    if (std::abs(trackDCAZ) < maxDCAZ) {
      SETBIT(bit, JTrackSel::qualityTrackWDCA);
    }
  }
  if (track.trackCutFlagFb5()) {
    SETBIT(bit, JTrackSel::hybridTrack);
  }
  return bit;
}

uint8_t setSingleTrackSelectionBit(int trackSelection)
{
  uint8_t bit = 0;
  if (trackSelection != -1) {
    SETBIT(bit, trackSelection);
  }
  return bit;
}

template <typename T>
float trackEnergy(T const& track, float mass = mPion)
{
  return std::sqrt((track.p() * track.p()) + (mass * mass));
}

template <typename T>
bool selectTrackDcaZ(T const& track, double dcaZmax = 99.)
{
  return std::abs(track.dcaZ()) < dcaZmax;
}

} // namespace jetderiveddatautilities

#endif // PWGJE_CORE_JETDERIVEDDATAUTILITIES_H_
