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

/// \file EventSelectionParams.cxx
/// \brief Event selection parameters
///
/// \author Evgeny Kryshen <evgeny.kryshen@cern.ch> and Igor Altsybeev <Igor.Altsybeev@cern.ch>

// o2-linter: disable=name/workflow-file

#include "EventSelectionParams.h"

#include <cstddef>
#include <cstring>

namespace o2::aod::evsel
{
const char* selectionLabels[kNsel] = {
  "kIsBBV0A",
  "kIsBBV0C",
  "kIsBBFDA",
  "kIsBBFDC",
  "kIsBBT0A",
  "kIsBBT0C",
  "kNoBGV0A",
  "kNoBGV0C",
  "kNoBGFDA",
  "kNoBGFDC",
  "kNoBGT0A",
  "kNoBGT0C",
  "kIsBBZNA",
  "kIsBBZNC",
  "kIsBBZAC",
  "kNoBGZNA",
  "kNoBGZNC",
  "kNoV0MOnVsOfPileup",
  "kNoSPDOnVsOfPileup",
  "kNoV0Casymmetry",
  "kIsGoodTimeRange",
  "kNoIncompleteDAQ",
  "kNoTPCLaserWarmUp",
  "kNoTPCHVdip",
  "kNoPileupFromSPD",
  "kNoV0PFPileup",
  "kNoSPDClsVsTklBG",
  "kNoV0C012vsTklBG",
  "kNoInconsistentVtx",
  "kNoPileupInMultBins",
  "kNoPileupMV",
  "kNoPileupTPC",
  "kIsTriggerTVX",
  "kIsINT1",
  "kNoITSROFrameBorder",
  "kNoTimeFrameBorder",
  "kNoSameBunchPileup",
  "kIsGoodZvtxFT0vsPV",
  "kIsVertexITSTPC",
  "kIsVertexTOFmatched",
  "kIsVertexTRDmatched",
  "kNoCollInTimeRangeNarrow",
  "kNoCollInTimeRangeStrict",
  "kNoCollInTimeRangeStandard",
  "kNoCollInRofStrict",
  "kNoCollInRofStandard",
  "kNoHighMultCollInPrevRof",
  "kIsGoodITSLayer3",
  "kIsGoodITSLayer0123",
  "kIsGoodITSLayersAll"};
} // namespace o2::aod::evsel

using namespace o2::aod::evsel;

EventSelectionParams::EventSelectionParams(int system, int run)
{
  memset(selectionBarrel, 0, sizeof(selectionBarrel));
  memset(selectionMuonWithPileupCuts, 0, sizeof(selectionMuonWithPileupCuts));
  memset(selectionMuonWithoutPileupCuts, 0, sizeof(selectionMuonWithoutPileupCuts));

  if (run == 1) {
    selectionBarrel[kIsBBV0A] = 1;
    selectionBarrel[kIsBBV0C] = 1;
    selectionBarrel[kNoTPCHVdip] = 1;
    selectionBarrel[kNoIncompleteDAQ] = 1;
    selectionBarrel[kNoTPCLaserWarmUp] = 1;

    selectionMuonWithPileupCuts[kIsBBV0A] = 1;
    selectionMuonWithPileupCuts[kIsBBV0C] = 1;

    selectionMuonWithoutPileupCuts[kIsBBV0A] = 1;
    selectionMuonWithoutPileupCuts[kIsBBV0C] = 1;

    if (system == 1) { // additional ZDC checks in pA
      selectionBarrel[kNoBGZNA] = 1;
      selectionMuonWithPileupCuts[kNoBGZNA] = 1;
      selectionMuonWithoutPileupCuts[kNoBGZNA] = 1;
    }
    if (system == 2) { // additional ZDC checks in Ap
      selectionBarrel[kNoBGZNC] = 1;
      selectionMuonWithPileupCuts[kNoBGZNC] = 1;
      selectionMuonWithoutPileupCuts[kNoBGZNC] = 1;
    }
    if (system == 3) { // AA
      selectionBarrel[kIsBBZAC] = 1;
      selectionMuonWithPileupCuts[kIsBBZAC] = 1;
      selectionMuonWithoutPileupCuts[kIsBBZAC] = 1;
    }
  } else if (run == 2) {
    if (system == 0 || system == 1 || system == 2) { // pp, pA or Ap
      // default barrel selection
      selectionBarrel[kIsBBV0A] = 1;
      selectionBarrel[kIsBBV0C] = 1;
      selectionBarrel[kNoV0MOnVsOfPileup] = 1;
      selectionBarrel[kNoSPDOnVsOfPileup] = 1;
      selectionBarrel[kNoV0Casymmetry] = 1;
      selectionBarrel[kIsGoodTimeRange] = 1;
      selectionBarrel[kNoIncompleteDAQ] = 1;
      selectionBarrel[kNoTPCHVdip] = 1;
      selectionBarrel[kNoPileupFromSPD] = 1;
      selectionBarrel[kNoV0PFPileup] = 1;
      selectionBarrel[kNoSPDClsVsTklBG] = 1;
      selectionBarrel[kNoV0C012vsTklBG] = 1;

      // similar to Barrel but without kNoTPCHVdip and kIsGoodTimeRange checks
      selectionMuonWithPileupCuts[kIsBBV0A] = 1;
      selectionMuonWithPileupCuts[kIsBBV0C] = 1;
      selectionMuonWithPileupCuts[kNoV0MOnVsOfPileup] = 1;
      selectionMuonWithPileupCuts[kNoSPDOnVsOfPileup] = 1;
      selectionMuonWithPileupCuts[kNoV0Casymmetry] = 1;
      selectionMuonWithPileupCuts[kNoIncompleteDAQ] = 1;
      selectionMuonWithPileupCuts[kNoPileupFromSPD] = 1;
      selectionMuonWithPileupCuts[kNoV0PFPileup] = 1;
      selectionMuonWithPileupCuts[kNoSPDClsVsTklBG] = 1;
      selectionMuonWithPileupCuts[kNoV0C012vsTklBG] = 1;

      // basic checks for muon analyses without in/out-of-bunch pileup rejection
      selectionMuonWithoutPileupCuts[kIsBBV0A] = 1;
      selectionMuonWithoutPileupCuts[kIsBBV0C] = 1;
      selectionMuonWithoutPileupCuts[kNoV0Casymmetry] = 1;
      selectionMuonWithoutPileupCuts[kNoIncompleteDAQ] = 1;
      selectionMuonWithoutPileupCuts[kNoV0C012vsTklBG] = 1;
    }

    if (system == 1) { // additional ZDC checks in pA
      selectionBarrel[kNoBGZNA] = 1;
      selectionMuonWithPileupCuts[kNoBGZNA] = 1;
      selectionMuonWithoutPileupCuts[kNoBGZNA] = 1;
    }

    if (system == 2) { // additional ZDC checks in Ap
      selectionBarrel[kNoBGZNC] = 1;
      selectionMuonWithPileupCuts[kNoBGZNC] = 1;
      selectionMuonWithoutPileupCuts[kNoBGZNC] = 1;
    }

    if (system == 3) { // AA
      selectionBarrel[kIsBBV0A] = 1;
      selectionBarrel[kIsBBV0C] = 1;
      selectionBarrel[kIsBBZAC] = 1;
      selectionBarrel[kIsGoodTimeRange] = 1;
      selectionBarrel[kNoTPCHVdip] = 1;

      selectionMuonWithPileupCuts[kIsBBV0A] = 1;
      selectionMuonWithPileupCuts[kIsBBV0C] = 1;
      selectionMuonWithPileupCuts[kIsBBZAC] = 1;

      selectionMuonWithoutPileupCuts[kIsBBV0A] = 1;
      selectionMuonWithoutPileupCuts[kIsBBV0C] = 1;
      selectionMuonWithoutPileupCuts[kIsBBZAC] = 1;
    }
  }
}

void EventSelectionParams::disableOutOfBunchPileupCuts()
{
  selectionBarrel[kNoV0MOnVsOfPileup] = 0;
  selectionBarrel[kNoSPDOnVsOfPileup] = 0;
  selectionBarrel[kNoV0Casymmetry] = 0;
  selectionBarrel[kNoV0PFPileup] = 0;

  selectionMuonWithPileupCuts[kNoV0MOnVsOfPileup] = 0;
  selectionMuonWithPileupCuts[kNoSPDOnVsOfPileup] = 0;
  selectionMuonWithPileupCuts[kNoV0Casymmetry] = 0;
  selectionMuonWithPileupCuts[kNoV0PFPileup] = 0;

  selectionMuonWithoutPileupCuts[kNoV0MOnVsOfPileup] = 0;
  selectionMuonWithoutPileupCuts[kNoSPDOnVsOfPileup] = 0;
  selectionMuonWithoutPileupCuts[kNoV0Casymmetry] = 0;
  selectionMuonWithoutPileupCuts[kNoV0PFPileup] = 0;
}

void EventSelectionParams::setOnVsOfParams(float newV0MOnVsOfA, float newV0MOnVsOfB, float newSPDOnVsOfA, float newSPDOnVsOfB)
{
  fV0MOnVsOfA = newV0MOnVsOfA;
  fV0MOnVsOfB = newV0MOnVsOfB;
  fSPDOnVsOfA = newSPDOnVsOfA;
  fSPDOnVsOfB = newSPDOnVsOfB;
}

bool* EventSelectionParams::getSelection(int iSelection)
{
  if (iSelection == 0) {
    return selectionBarrel;
  } else if (iSelection == 1) {
    return selectionMuonWithPileupCuts;
  } else if (iSelection == 2) {
    return selectionMuonWithoutPileupCuts;
  }
  return NULL;
}
