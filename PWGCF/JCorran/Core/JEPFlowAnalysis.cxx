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

void JEPFlowAnalysis::fillVnHistograms(const Int_t harmN, Float_t cent, Float_t det, Float_t pT, Float_t vn, Float_t vn_sin)
{
  switch (harmN) {
    case 2: {
      mHistRegistry->fill(HIST("fV2EP"), vn, det, pT, cent, 1.);
      mHistRegistry->fill(HIST("fV2EP_sin"), vn_sin, det, pT, cent, 1.);
    } break;
    case 3: {
      mHistRegistry->fill(HIST("fV3EP"), vn, det, pT, cent, 1.);
      mHistRegistry->fill(HIST("fV3EP_sin"), vn_sin, det, pT, cent, 1.);
    } break;
    case 4: {
      mHistRegistry->fill(HIST("fV4EP"), vn, det, pT, cent, 1.);
      mHistRegistry->fill(HIST("fV4EP_sin"), vn_sin, det, pT, cent, 1.);
    } break;
    default:
      break;
  }
}

void JEPFlowAnalysis::fillResolutionHistograms(Float_t cent, Float_t harmN, Float_t ResNumA, Float_t ResNumB, Float_t ResDenom)
{
  mHistRegistry->fill(HIST("fResNumA"), ResNumA, harmN, cent, 1.);
  mHistRegistry->fill(HIST("fResNumB"), ResNumB, harmN, cent, 1.);
  mHistRegistry->fill(HIST("fResDenom"), ResDenom, harmN, cent, 1.);
}
