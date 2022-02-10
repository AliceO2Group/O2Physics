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

#ifndef EventSelectionParams_H
#define EventSelectionParams_H

#include <Rtypes.h>
#include <TMath.h>

namespace evsel
{
// Event selection criteria
enum EventSelectionFlags {
  kIsBBV0A = 0,
  kIsBBV0C,
  kIsBBFDA,
  kIsBBFDC,
  kNoBGV0A,
  kNoBGV0C,
  kNoBGFDA,
  kNoBGFDC,
  kIsBBT0A,
  kIsBBT0C,
  kIsBBZNA,
  kIsBBZNC,
  kIsBBZAC,
  kNoBGZNA,
  kNoBGZNC,
  kNoV0MOnVsOfPileup,
  kNoSPDOnVsOfPileup,
  kNoV0Casymmetry,
  kIsGoodTimeRange,
  kNoIncompleteDAQ,
  kNoTPCLaserWarmUp,
  kNoTPCHVdip,
  kNoPileupFromSPD,
  kNoV0PFPileup,
  kNoSPDClsVsTklBG,
  kNoV0C012vsTklBG,
  kNsel
};

extern const char* selectionLabels[kNsel];

} // namespace evsel

using namespace evsel;

class EventSelectionParams
{
 public:
  EventSelectionParams(int system = 0);
  void DisableOutOfBunchPileupCuts();
  void SetOnVsOfParams(float newV0MOnVsOfA, float newV0MOnVsOfB, float newSPDOnVsOfA, float newSPDOnVsOfB);
  bool* GetSelection(int iSelection);

  bool selectionBarrel[kNsel];
  bool selectionMuonWithPileupCuts[kNsel];
  bool selectionMuonWithoutPileupCuts[kNsel];

  // time-of-flight offset
  float fV0ADist = 329.00 / TMath::Ccgs() * 1e9; // ns
  float fV0CDist = 87.15 / TMath::Ccgs() * 1e9;  // ns

  float fFDADist = (1695.30 + 1698.04) / 2. / TMath::Ccgs() * 1e9; // ns
  float fFDCDist = (1952.90 + 1955.90) / 2. / TMath::Ccgs() * 1e9; // ns

  // beam-beam and beam-gas windows
  float fV0ABBlower = +fV0ADist - 9.5;  // ns
  float fV0ABBupper = +fV0ADist + 22.5; // ns
  float fV0ABGlower = -fV0ADist - 2.5;  // ns
  float fV0ABGupper = -fV0ADist + 5.0;  // ns
  float fV0CBBlower = +fV0CDist - 2.5;  // ns
  float fV0CBBupper = +fV0CDist + 22.5; // ns
  float fV0CBGlower = -fV0CDist - 2.5;  // ns
  float fV0CBGupper = -fV0CDist + 2.5;  // ns

  float fFDABBlower = +fFDADist - 2.5; // ns
  float fFDABBupper = +fFDADist + 2.5; // ns
  float fFDABGlower = -fFDADist - 4.0; // ns
  float fFDABGupper = -fFDADist + 4.0; // ns
  float fFDCBBlower = +fFDCDist - 1.5; // ns
  float fFDCBBupper = +fFDCDist + 1.5; // ns
  float fFDCBGlower = -fFDCDist - 2.0; // ns
  float fFDCBGupper = -fFDCDist + 2.0; // ns

  float fZNDifMean = 0;  // ns
  float fZNSumMean = 0;  // ns
  float fZNDifSigma = 2; // ns
  float fZNSumSigma = 2; // ns

  float fZNABBlower = -2.0;  // ns
  float fZNABBupper = 2.0;   // ns
  float fZNCBBlower = -2.0;  // ns
  float fZNCBBupper = 2.0;   // ns
  float fZNABGlower = 5.0;   // ns
  float fZNABGupper = 100.0; // ns
  float fZNCBGlower = 5.0;   // ns
  float fZNCBGupper = 100.0; // ns

  // TODO rough cuts to be adjusted
  float fT0ABBlower = -2.0; // ns
  float fT0ABBupper = 2.0;  // ns
  float fT0CBBlower = -2.0; // ns
  float fT0CBBupper = 2.0;  // ns

  // Default values from AliOADBTriggerAnalysis constructor
  float fSPDClsVsTklA = 65.f;
  float fSPDClsVsTklB = 4.f;
  float fV0C012vsTklA = 150.f;
  float fV0C012vsTklB = 20.f;
  float fV0MOnVsOfA = -59.56f;
  float fV0MOnVsOfB = 5.22f;
  float fSPDOnVsOfA = -5.62f;
  float fSPDOnVsOfB = 0.85f;
  float fV0CasymA = -25.f;
  float fV0CasymB = 0.15f;

  ClassDefNV(EventSelectionParams, 1)
};

#endif
