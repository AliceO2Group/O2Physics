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

// see header for a more detailed description.

#include "v0SelectionGroup.h"

#include "PWGLF/Utils/v0SelectionBits.h"

#include <Framework/Logger.h>

#include <cstdint>

// simple checkers, but ensure 64 bit integers
#define bitset(var, nbit) ((var) |= (static_cast<uint64_t>(1) << static_cast<uint64_t>(nbit)))
#define bitcheck(var, nbit) ((var) & (static_cast<uint64_t>(1) << static_cast<uint64_t>(nbit)))

// trivial 64-bit mask checker
bool v0SelectionGroup::verifyMask(uint64_t bitmap, uint64_t mask) const
{
  return (bitmap & mask) == mask;
}

// provides standard masks given current event selection criteria.
void v0SelectionGroup::provideMasks(uint64_t& maskTopological, uint64_t& maskTrackProperties, uint64_t& maskK0ShortSpecific, uint64_t& maskLambdaSpecific, uint64_t& maskAntiLambdaSpecific) const
{
  maskTopological = (static_cast<uint64_t>(1) << v0data::selCosPA) | (static_cast<uint64_t>(1) << v0data::selRadius) | (static_cast<uint64_t>(1) << v0data::selDCANegToPV) | (static_cast<uint64_t>(1) << v0data::selDCAPosToPV) | (static_cast<uint64_t>(1) << v0data::selDCAV0Dau) | (static_cast<uint64_t>(1) << v0data::selRadiusMax);
  maskK0ShortSpecific = (static_cast<uint64_t>(1) << v0data::selK0ShortRapidity) | (static_cast<uint64_t>(1) << v0data::selK0ShortCTau) | (static_cast<uint64_t>(1) << v0data::selK0ShortArmenteros) | (static_cast<uint64_t>(1) << v0data::selConsiderK0Short);
  maskLambdaSpecific = (static_cast<uint64_t>(1) << v0data::selLambdaRapidity) | (static_cast<uint64_t>(1) << v0data::selLambdaCTau) | (static_cast<uint64_t>(1) << v0data::selConsiderLambda);
  maskAntiLambdaSpecific = (static_cast<uint64_t>(1) << v0data::selLambdaRapidity) | (static_cast<uint64_t>(1) << v0data::selLambdaCTau) | (static_cast<uint64_t>(1) << v0data::selConsiderAntiLambda);

  // ask for specific TPC/TOF PID v0data::selections
  maskTrackProperties = 0;
  if (requirePosITSonly) {
    maskTrackProperties = maskTrackProperties | (static_cast<uint64_t>(1) << v0data::selPosItsOnly) | (static_cast<uint64_t>(1) << v0data::selPosGoodITSTrack);
  } else {
    maskTrackProperties = maskTrackProperties | (static_cast<uint64_t>(1) << v0data::selPosGoodTPCTrack) | (static_cast<uint64_t>(1) << v0data::selPosGoodITSTrack);
    // TPC signal is available: ask for positive track PID
    if (TpcPidNsigmaCut < 1e+5) { // safeguard for no cut
      maskK0ShortSpecific = maskK0ShortSpecific | (static_cast<uint64_t>(1) << v0data::selTPCPIDPositivePion);
      maskLambdaSpecific = maskLambdaSpecific | (static_cast<uint64_t>(1) << v0data::selTPCPIDPositiveProton);
      maskAntiLambdaSpecific = maskAntiLambdaSpecific | (static_cast<uint64_t>(1) << v0data::selTPCPIDPositivePion);
    }
    // TOF PID
    if (TofPidNsigmaCutK0Pi < 1e+5) // safeguard for no cut
      maskK0ShortSpecific = maskK0ShortSpecific | (static_cast<uint64_t>(1) << v0data::selTOFNSigmaPositivePionK0Short) | (static_cast<uint64_t>(1) << v0data::selTOFDeltaTPositivePionK0Short);
    if (TofPidNsigmaCutLaPr < 1e+5) // safeguard for no cut
      maskLambdaSpecific = maskLambdaSpecific | (static_cast<uint64_t>(1) << v0data::selTOFNSigmaPositiveProtonLambda) | (static_cast<uint64_t>(1) << v0data::selTOFDeltaTPositiveProtonLambda);
    if (TofPidNsigmaCutLaPi < 1e+5) // safeguard for no cut
      maskAntiLambdaSpecific = maskAntiLambdaSpecific | (static_cast<uint64_t>(1) << v0data::selTOFNSigmaPositivePionLambda) | (static_cast<uint64_t>(1) << v0data::selTOFDeltaTPositivePionLambda);
  }
  if (requireNegITSonly) {
    maskTrackProperties = maskTrackProperties | (static_cast<uint64_t>(1) << v0data::selNegItsOnly) | (static_cast<uint64_t>(1) << v0data::selNegGoodITSTrack);
  } else {
    maskTrackProperties = maskTrackProperties | (static_cast<uint64_t>(1) << v0data::selNegGoodTPCTrack) | (static_cast<uint64_t>(1) << v0data::selNegGoodITSTrack);
    // TPC signal is available: ask for negative track PID
    if (TpcPidNsigmaCut < 1e+5) { // safeguard for no cut
      maskK0ShortSpecific = maskK0ShortSpecific | (static_cast<uint64_t>(1) << v0data::selTPCPIDNegativePion);
      maskLambdaSpecific = maskLambdaSpecific | (static_cast<uint64_t>(1) << v0data::selTPCPIDNegativePion);
      maskAntiLambdaSpecific = maskAntiLambdaSpecific | (static_cast<uint64_t>(1) << v0data::selTPCPIDNegativeProton);
    }
    // TOF PID
    if (TofPidNsigmaCutK0Pi < 1e+5) // safeguard for no cut
      maskK0ShortSpecific = maskK0ShortSpecific | (static_cast<uint64_t>(1) << v0data::selTOFNSigmaNegativePionK0Short) | (static_cast<uint64_t>(1) << v0data::selTOFDeltaTNegativePionK0Short);
    if (TofPidNsigmaCutLaPr < 1e+5) // safeguard for no cut
      maskLambdaSpecific = maskLambdaSpecific | (static_cast<uint64_t>(1) << v0data::selTOFNSigmaNegativePionLambda) | (static_cast<uint64_t>(1) << v0data::selTOFDeltaTNegativePionLambda);
    if (TofPidNsigmaCutLaPi < 1e+5) // safeguard for no cut
      maskAntiLambdaSpecific = maskAntiLambdaSpecific | (static_cast<uint64_t>(1) << v0data::selTOFNSigmaNegativeProtonLambda) | (static_cast<uint64_t>(1) << v0data::selTOFDeltaTNegativeProtonLambda);
  }

  if (skipTPConly) {
    maskK0ShortSpecific = maskK0ShortSpecific | (static_cast<uint64_t>(1) << v0data::selPosNotTPCOnly) | (static_cast<uint64_t>(1) << v0data::selNegNotTPCOnly);
    maskLambdaSpecific = maskLambdaSpecific | (static_cast<uint64_t>(1) << v0data::selPosNotTPCOnly) | (static_cast<uint64_t>(1) << v0data::selNegNotTPCOnly);
    maskAntiLambdaSpecific = maskAntiLambdaSpecific | (static_cast<uint64_t>(1) << v0data::selPosNotTPCOnly) | (static_cast<uint64_t>(1) << v0data::selNegNotTPCOnly);
  }
}

void v0SelectionGroup::PrintSelections() const
{
  LOGF(info, "+++ Phase space selections ++++++++++++++++++++++++");
  LOGF(info, "Rapidity cut ...........................: %.2f", rapidityCut);
  LOGF(info, "Eta (daughters) cut ....................: %.2f", daughterEtaCut);
  LOGF(info, "+++ Topological selections ++++++++++++++++++++++++");
  LOGF(info, "V0 cosine of PA ........................: %.2f", v0cospa);
  LOGF(info, "DCA V0 daughters .......................: %.2f", dcav0dau);
  LOGF(info, "Neg track DCA to PV ....................: %.2f", dcanegtopv);
  LOGF(info, "Pos track DCA to PV ....................: %.2f", dcapostopv);
  LOGF(info, "Minimum radius .........................: %.2f", v0radius);
  LOGF(info, "Maximum radius .........................: %.2f", v0radiusMax);
  LOGF(info, "+++ Track quality selections ++++++++++++++++++++++");
  LOGF(info, "Minimum TPC rows .......................: %i", minTPCrows);
  LOGF(info, "Minimum ITS clusters ...................: %i", minITSclusters);
  LOGF(info, "Skip TPC only ..........................: %s", skipTPConly ? "true" : "false");
  LOGF(info, "Require positive ITS only ..............: %s", requirePosITSonly ? "true" : "false");
  LOGF(info, "Require negative ITS only ..............: %s", requireNegITSonly ? "true" : "false");
  LOGF(info, "+++ Particle identification selections ++++++++++++");
  LOGF(info, "TPC PID Nsigma .........................: %.2f", TpcPidNsigmaCut);
  LOGF(info, "TOF PID Nsigma LaPr ....................: %.2f", TofPidNsigmaCutLaPr);
  LOGF(info, "TOF PID Nsigma LaPi ....................: %.2f", TofPidNsigmaCutLaPi);
  LOGF(info, "TOF PID Nsigma K0Pi ....................: %.2f", TofPidNsigmaCutK0Pi);
  LOGF(info, "+++ Misc selections ++++++++++++++++++++++++++++++");
  LOGF(info, "K0Short lifetime cut ...................: %.2f", lifetimeCutK0Short);
  LOGF(info, "Lambda lifetime cut ....................: %.2f", lifetimeCutLambda);
  LOGF(info, "Armenteros podolanski parameter ........: %.2f", armPodCut);
}
