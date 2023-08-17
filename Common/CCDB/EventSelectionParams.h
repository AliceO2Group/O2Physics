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

#ifndef COMMON_CCDB_EVENTSELECTIONPARAMS_H_
#define COMMON_CCDB_EVENTSELECTIONPARAMS_H_

#include <Rtypes.h>
#include <TMath.h>

namespace o2::aod::evsel
{
// Event selection criteria
enum EventSelectionFlags {
  kIsBBV0A = 0,        // cell-averaged time in V0A in beam-beam window
  kIsBBV0C,            // cell-averaged time in V0C in beam-beam window (for Run 2 only)
  kIsBBFDA,            // cell-averaged time in FDA (or AD in Run2) in beam-beam window
  kIsBBFDC,            // cell-averaged time in FDC (or AD in Run2) in beam-beam window
  kIsBBT0A,            // cell-averaged time in T0A in beam-beam window
  kIsBBT0C,            // cell-averaged time in T0C in beam-beam window
  kNoBGV0A,            // cell-averaged time in V0A in beam-gas window
  kNoBGV0C,            // cell-averaged time in V0C in beam-gas window (for Run 2 only)
  kNoBGFDA,            // cell-averaged time in FDA (AD in Run2) in beam-gas window
  kNoBGFDC,            // cell-averaged time in FDC (AD in Run2) in beam-gas window
  kNoBGT0A,            // cell-averaged time in T0A in beam-gas window
  kNoBGT0C,            // cell-averaged time in T0C in beam-gas window
  kIsBBZNA,            // time in common ZNA channel in beam-beam window
  kIsBBZNC,            // time in common ZNC channel in beam-beam window
  kIsBBZAC,            // time in ZNA and ZNC in beam-beam window - circular cut in ZNA-ZNC plane
  kNoBGZNA,            // time in common ZNA channel is outside of beam-gas window
  kNoBGZNC,            // time in common ZNC channel is outside of beam-gas window
  kNoV0MOnVsOfPileup,  // no out-of-bunch pileup according to online-vs-offline VOM correlation
  kNoSPDOnVsOfPileup,  // no out-of-bunch pileup according to online-vs-offline SPD correlation
  kNoV0Casymmetry,     // no beam-gas according to correlation of V0C multiplicities in V0C3 and V0C012
  kIsGoodTimeRange,    // good time range
  kNoIncompleteDAQ,    // complete event according to DAQ flags
  kNoTPCLaserWarmUp,   // no TPC laser warm-up event (used in Run 1)
  kNoTPCHVdip,         // no TPC HV dip
  kNoPileupFromSPD,    // no pileup according to SPD vertexer
  kNoV0PFPileup,       // no out-of-bunch pileup according to V0 past-future info
  kNoSPDClsVsTklBG,    // no beam-gas according to cluster-vs-tracklet correlation
  kNoV0C012vsTklBG,    // no beam-gas according to V0C012-vs-tracklet correlation
  kNoInconsistentVtx,  // no inconsistency in SPD and Track vertices
  kNoPileupInMultBins, // no pileup according to multiplicity-differential pileup checks
  kNoPileupMV,         // no pileup according to multi-vertexer
  kNoPileupTPC,        // no pileup in TPC
  kIsTriggerTVX,       // FT0 vertex (acceptable FT0C-FT0A time difference) at trigger level
  kIsINT1,             // SPDGFO >= 1 || V0A || V0C
  kNsel                // counter
};

extern const char* selectionLabels[kNsel];

} // namespace o2::aod::evsel

class EventSelectionParams
{
 public:
  explicit EventSelectionParams(int system = 0, int run = 2);
  void DisableOutOfBunchPileupCuts();
  void SetOnVsOfParams(float newV0MOnVsOfA, float newV0MOnVsOfB, float newSPDOnVsOfA, float newSPDOnVsOfB);
  bool* GetSelection(int iSelection);

  bool selectionBarrel[o2::aod::evsel::kNsel];
  bool selectionMuonWithPileupCuts[o2::aod::evsel::kNsel];
  bool selectionMuonWithoutPileupCuts[o2::aod::evsel::kNsel];

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

  float fT0ABBlower = -2.0; // ns
  float fT0ABBupper = 2.0;  // ns
  float fT0CBBlower = -2.0; // ns
  float fT0CBBupper = 2.0;  // ns
  float fT0ABGlower = 32.7; // ns
  float fT0ABGupper = 32.8; // ns
  float fT0CBGlower = 32.7; // ns
  float fT0CBGupper = 32.8; // ns

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

  ClassDefNV(EventSelectionParams, 3)
};

#endif // COMMON_CCDB_EVENTSELECTIONPARAMS_H_
