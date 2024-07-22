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

#include "JEPFlowAnalysis.h"

using namespace o2;
using namespace o2::framework;
using namespace std;

void JEPFlowAnalysis::FillHistograms(const Int_t cBin, Float_t det, Float_t v2, Float_t v3, Float_t v4) {
  switch (cBin) {
    case 0: {
      mHistRegistry->fill(HIST(mCentClasses[0]) + HIST("fV2EP"), v2, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[0]) + HIST("fV3EP"), v3, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[0]) + HIST("fV4EP"), v4, det, 1.);
    } break;
    case 1: {
      mHistRegistry->fill(HIST(mCentClasses[1]) + HIST("fV2EP"), v2, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[1]) + HIST("fV3EP"), v3, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[1]) + HIST("fV4EP"), v4, det, 1.);
    } break;
    case 2: {
      mHistRegistry->fill(HIST(mCentClasses[2]) + HIST("fV2EP"), v2, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[2]) + HIST("fV3EP"), v3, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[2]) + HIST("fV4EP"), v4, det, 1.);
    } break;
    case 3: {
      mHistRegistry->fill(HIST(mCentClasses[3]) + HIST("fV2EP"), v2, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[3]) + HIST("fV3EP"), v3, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[3]) + HIST("fV4EP"), v4, det, 1.);
    } break;
    case 4: {
      mHistRegistry->fill(HIST(mCentClasses[4]) + HIST("fV2EP"), v2, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[4]) + HIST("fV3EP"), v3, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[4]) + HIST("fV4EP"), v4, det, 1.);
    } break;
    case 5: {
      mHistRegistry->fill(HIST(mCentClasses[5]) + HIST("fV2EP"), v2, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[5]) + HIST("fV3EP"), v3, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[5]) + HIST("fV4EP"), v4, det, 1.);
    } break;
    case 6: {
      mHistRegistry->fill(HIST(mCentClasses[6]) + HIST("fV2EP"), v2, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[6]) + HIST("fV3EP"), v3, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[6]) + HIST("fV4EP"), v4, det, 1.);
    } break;
    case 7: {
      mHistRegistry->fill(HIST(mCentClasses[7]) + HIST("fV2EP"), v2, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[7]) + HIST("fV3EP"), v3, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[7]) + HIST("fV4EP"), v4, det, 1.);
    } break;
    default:
      return;
  }
}


void JEPFlowAnalysis::FillVnHistograms(Float_t cent, Float_t det, Float_t pT, Float_t v2, Float_t v3, Float_t v4) {
  mHistRegistry->fill(HIST("fV2EP"), v2, det, pT, cent, 1.);
  mHistRegistry->fill(HIST("fV3EP"), v3, det, pT, cent, 1.);
  mHistRegistry->fill(HIST("fV4EP"), v4, det, pT, cent, 1.);
}

void JEPFlowAnalysis::FillResolutionHistograms(Float_t cent, Float_t det, Float_t harmN, Float_t ResNumA, Float_t ResNumB, Float_t ResDenom) {
  mHistRegistry->fill(HIST("fResNumA"), ResNumA, det, harmN, cent, 1.);
  mHistRegistry->fill(HIST("fResNumB"), ResNumB, det, harmN, cent, 1.);
  mHistRegistry->fill(HIST("fResDenom"), ResDenom, det, harmN, cent, 1.);
}