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

/// \file JetDetivedDataUtilities.h
/// \brief Jet derived data related utilities
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_CORE_JETDERIVEDDATAUTILITIES_H_
#define PWGJE_CORE_JETDERIVEDDATAUTILITIES_H_

#include <string>
#include "Common/CCDB/TriggerAliases.h"
#include "Common/CCDB/EventSelectionParams.h"

namespace jetderiveddatautilities
{

static constexpr float mPion = 0.139; // TDatabasePDG::Instance()->GetParticle(211)->Mass(); //can be removed when pion mass becomes default for unidentified tracks

enum JCollisionSel {
  sel8 = 0,
  sel8WithoutTimeFrameBorderCut = 1,
  sel7 = 2,
  sel7WithoutTimeFrameBorderCut = 3,
  WithoutTimeFrameBorderCut = 4
};

template <typename T>
bool selectCollision(T const& collision, int eventSelection = -1)
{
  if (eventSelection == -1) {
    return true;
  }
  return (collision.eventSel() & (1 << eventSelection));
}

int initialiseEventSelection(std::string eventSelection)
{
  if (eventSelection == "sel8") {
    return JCollisionSel::sel8;
  } else if (eventSelection == "sel8WithoutTimeFrameBorderCut") {
    return JCollisionSel::sel8WithoutTimeFrameBorderCut;
  } else if (eventSelection == "sel7") {
    return JCollisionSel::sel7;
  } else if (eventSelection == "sel7WithoutTimeFrameBorderCut") {
    return JCollisionSel::sel7WithoutTimeFrameBorderCut;
  } else if (eventSelection == "WithoutTimeFrameBorderCut") {
    return JCollisionSel::WithoutTimeFrameBorderCut;
  }
  return -1;
}

template <typename T>
uint16_t setEventSelectionBit(T const& collision)
{

  uint16_t bit = 0;

  if (!collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
    SETBIT(bit, JCollisionSel::WithoutTimeFrameBorderCut);
  }
  if (collision.sel8()) {
    if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      SETBIT(bit, JCollisionSel::sel8);
    } else {
      SETBIT(bit, JCollisionSel::sel8WithoutTimeFrameBorderCut);
    }
  }
  if (collision.sel7()) {
    if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      SETBIT(bit, JCollisionSel::sel7);
    } else {
      SETBIT(bit, JCollisionSel::sel7WithoutTimeFrameBorderCut);
    }
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

enum JTrigSelCh {
  noChargedTigger = 0,
  chargedLow = 1,
  chargedHigh = 2
};

template <typename T>
bool selectChargedTrigger(T const& collision, int triggerSelection)
{
  if (triggerSelection == -1) {
    return true;
  }
  return (collision.chargedTriggerSel() & (1 << triggerSelection));
}

int initialiseChargedTriggerSelection(std::string triggerSelection)
{
  if (triggerSelection == "chargedLow") {
    return JTrigSelCh::chargedLow;
  }
  if (triggerSelection == "chargedHigh") {
    return JTrigSelCh::chargedHigh;
  }
  return -1;
}

template <typename T>
uint8_t setChargedTriggerSelectionBit(T const& collision)
{

  uint8_t bit = 0;
  if (collision.hasJetChLowPt()) {
    SETBIT(bit, JTrigSelCh::chargedLow);
  }
  if (collision.hasJetChHighPt()) {
    SETBIT(bit, JTrigSelCh::chargedHigh);
  }
  return bit;
}

enum JTrigSelFull {
  noFullTrigger = 0,
  fullHigh = 1,
  fullLow = 2,
  neutralHigh = 3,
  neutralLow = 4,
  gammaVeryHighEMCAL = 5,
  gammaHighEMCAL = 6,
  gammaLowEMCAL = 7,
  gammaVeryLowEMCAL = 8,
  gammaVeryHighDCAL = 9,
  gammaHighDCAL = 10,
  gammaLowDCAL = 11,
  gammaVeryLowDCAL = 12
};

template <typename T>
bool selectFullTrigger(T const& collision, int triggerSelection)
{
  if (triggerSelection == -1) {
    return true;
  }
  return (collision.fullTriggerSel() & (1 << triggerSelection));
}

int initialiseFullTriggerSelection(std::string triggerSelection)
{
  if (triggerSelection == "fullHigh") {
    return JTrigSelFull::fullHigh;
  } else if (triggerSelection == "fullLow") {
    return JTrigSelFull::fullLow;
  } else if (triggerSelection == "neutralHigh") {
    return JTrigSelFull::neutralHigh;
  } else if (triggerSelection == "neutralLow") {
    return JTrigSelFull::neutralLow;
  } else if (triggerSelection == "gammaVeryHighEMCAL") {
    return JTrigSelFull::gammaVeryHighEMCAL;
  } else if (triggerSelection == "gammaHighEMCAL") {
    return JTrigSelFull::gammaHighEMCAL;
  } else if (triggerSelection == "gammaLowEMCAL") {
    return JTrigSelFull::gammaLowEMCAL;
  } else if (triggerSelection == "gammaVeryLowEMCAL") {
    return JTrigSelFull::gammaVeryLowEMCAL;
  } else if (triggerSelection == "gammaVeryHighDCAL") {
    return JTrigSelFull::gammaVeryHighDCAL;
  } else if (triggerSelection == "gammaHighDCAL") {
    return JTrigSelFull::gammaHighDCAL;
  } else if (triggerSelection == "gammaLowDCAL") {
    return JTrigSelFull::gammaLowDCAL;
  } else if (triggerSelection == "gammaVeryLowDCAL") {
    return JTrigSelFull::gammaVeryLowDCAL;
  }
  return -1;
}

template <typename T>
uint32_t setFullTriggerSelectionBit(T const& collision)
{
  uint32_t bit = 0;
  if (collision.hasJetFullHighPt()) {
    SETBIT(bit, JTrigSelFull::fullHigh);
  }
  if (collision.hasJetFullLowPt()) {
    SETBIT(bit, JTrigSelFull::fullLow);
  }
  if (collision.hasJetNeutralHighPt()) {
    SETBIT(bit, JTrigSelFull::neutralHigh);
  }
  if (collision.hasJetNeutralLowPt()) {
    SETBIT(bit, JTrigSelFull::neutralLow);
  }
  if (collision.hasGammaVeryHighPtEMCAL()) {
    SETBIT(bit, JTrigSelFull::gammaVeryHighEMCAL);
  }
  if (collision.hasGammaHighPtEMCAL()) {
    SETBIT(bit, JTrigSelFull::gammaHighEMCAL);
  }
  if (collision.hasGammaLowPtEMCAL()) {
    SETBIT(bit, JTrigSelFull::gammaLowEMCAL);
  }
  if (collision.hasGammaVeryLowPtEMCAL()) {
    SETBIT(bit, JTrigSelFull::gammaVeryLowEMCAL);
  }
  if (collision.hasGammaVeryHighPtDCAL()) {
    SETBIT(bit, JTrigSelFull::gammaVeryHighDCAL);
  }
  if (collision.hasGammaHighPtDCAL()) {
    SETBIT(bit, JTrigSelFull::gammaHighDCAL);
  }
  if (collision.hasGammaLowPtDCAL()) {
    SETBIT(bit, JTrigSelFull::gammaLowDCAL);
  }
  if (collision.hasGammaVeryLowPtDCAL()) {
    SETBIT(bit, JTrigSelFull::gammaVeryLowDCAL);
  }
  return bit;
}

enum JTrackSel {
  trackSign = 0, // warning : this number is hardcoded in the sign coloumn in the JTracks table so should not be changed without changing it there too
  globalTrack = 1,
  qualityTrack = 2,
  hybridTrack = 3
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

int initialiseTrackSelection(std::string trackSelection)
{
  if (trackSelection == "globalTracks") {
    return JTrackSel::globalTrack;
  } else if (trackSelection == "QualityTracks") {
    return JTrackSel::qualityTrack;
  } else if (trackSelection == "hybridTracksJE") {
    return JTrackSel::hybridTrack;
  }
  return -1;
}

template <typename T>
uint8_t setTrackSelectionBit(T const& track)
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

} // namespace jetderiveddatautilities

#endif // PWGJE_CORE_JETDERIVEDDATAUTILITIES_H_
