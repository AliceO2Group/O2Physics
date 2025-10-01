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

// \author  Maxim Virta (maxim.virta@cern.ch)

#ifndef PWGCF_JCORRAN_CORE_JEPFLOWANALYSIS_H_
#define PWGCF_JCORRAN_CORE_JEPFLOWANALYSIS_H_

#include <TComplex.h>

// O2 headers. //
#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;
using namespace std;

class JEPFlowAnalysis
{
 public:
  JEPFlowAnalysis() = default;
  void setHistRegistry(HistogramRegistry* histReg) { mHistRegistry = histReg; }

  void fillHistograms(const Int_t fCentBin, Float_t det, Float_t v2, Float_t v3, Float_t v4);
  void fillVnHistograms(const Int_t harmN, Float_t fCent, Float_t det, Float_t pT, Float_t vn, Float_t vn_sin);
  void fillResolutionHistograms(Float_t fCent, Float_t harmN, Float_t ResNumA, Float_t ResNumB, Float_t ResDenom);
  TComplex Q(const Int_t harmN, const Int_t p);

  void createHistograms()
  {
    if (!mHistRegistry) {
      LOGF(error, "Histogram registry missing. Quitting...");
      return;
    }

    mHistRegistry->add("FullCentrality", "FullCentrality", HistType::kTH1D, {{100, 0., 100.}}, true);
    mHistRegistry->add("fV2EP", "", {HistType::kTHnD, {{200, -1.05, 1.05}, {3, 0.5, 3.5}, {100, 0.2, 12.}, {20, 0., 100.}}}, true);     // x: v2_cos, y: detector, z: pT, t: centrality
    mHistRegistry->add("fV3EP", "", {HistType::kTHnD, {{200, -1.05, 1.05}, {3, 0.5, 3.5}, {100, 0.2, 12.}, {20, 0., 100.}}}, true);     // x: v2_cos, y: detector, z: pT, t: centrality
    mHistRegistry->add("fV4EP", "", {HistType::kTHnD, {{200, -1.05, 1.05}, {3, 0.5, 3.5}, {100, 0.2, 12.}, {20, 0., 100.}}}, true);     // x: v2_cos, y: detector, z: pT, t: centrality
    mHistRegistry->add("fV2EP_sin", "", {HistType::kTHnD, {{200, -1.05, 1.05}, {3, 0.5, 3.5}, {100, 0.2, 12.}, {20, 0., 100.}}}, true); // x: v2_sin, y: detector, z: pT, t: centrality
    mHistRegistry->add("fV3EP_sin", "", {HistType::kTHnD, {{200, -1.05, 1.05}, {3, 0.5, 3.5}, {100, 0.2, 12.}, {20, 0., 100.}}}, true); // x: v2_sin, y: detector, z: pT, t: centrality
    mHistRegistry->add("fV4EP_sin", "", {HistType::kTHnD, {{200, -1.05, 1.05}, {3, 0.5, 3.5}, {100, 0.2, 12.}, {20, 0., 100.}}}, true); // x: v2_sin, y: detector, z: pT, t: centrality
    mHistRegistry->add("fResNumA", "", {HistType::kTH3D, {{100, -1.05, 1.05}, {3, 1.5, 4.5}, {20, 0., 100.}}}, true);                   // x: resolution, y: harmonic, t: centrality
    mHistRegistry->add("fResNumB", "", {HistType::kTH3D, {{100, -1.05, 1.05}, {3, 1.5, 4.5}, {20, 0., 100.}}}, true);                   // x: resolution, y: harmonic, t: centrality
    mHistRegistry->add("fResDenom", "", {HistType::kTH3D, {{100, -1.05, 1.05}, {3, 1.5, 4.5}, {20, 0., 100.}}}, true);                  // x: resolution, y: harmonic, t: centrality
    mHistRegistry->add("phi", "Phi", {HistType::kTH1D, {{100, 0., TMath::TwoPi()}}}, true);
  }

 private:
  HistogramRegistry* mHistRegistry;

  ClassDefNV(JEPFlowAnalysis, 1);
};

#endif // PWGCF_JCORRAN_CORE_JEPFLOWANALYSIS_H_
