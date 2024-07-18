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
#include "Common/Core/EventPlaneHelper.h"

using namespace o2;
using namespace o2::framework;
using namespace std;

class JEPFlowAnalysis
{
public:
  JEPFlowAnalysis() = default;
  void SetHistRegistry(HistogramRegistry* histReg) { mHistRegistry = histReg; }

  void FillHistograms(const Int_t fCentBin, Float_t det, Float_t v2, Float_t v3, Float_t v4);
  TComplex Q(const Int_t harmN, const Int_t p);

  void CreateHistograms() {
    if (!mHistRegistry) {
      LOGF(error, "Histogram registry missing. Quitting...");
      return;
    }
    
    mHistRegistry->add("FullCentrality", "FullCentrality", HistType::kTH1D, {{100, 0., 100.}}, true);
    mHistRegistry->add("Centrality_0/fV2EP", "", {HistType::kTH2D, {{100, -0.15, 0.15},{3,0.5,3.5}}}, true);
    mHistRegistry->add("Centrality_0/fV3EP", "", {HistType::kTH2D, {{100, -0.15, 0.15},{3,0.5,3.5}}}, true);
    mHistRegistry->add("Centrality_0/fV4EP", "", {HistType::kTH2D, {{100, -0.15, 0.15},{3,0.5,3.5}}}, true);
    mHistRegistry->add("Centrality_0/phi", "Phi", {HistType::kTH1D, {{100, 0., TMath::TwoPi()}}}, true);

    for (UInt_t i = 1; i < 8; i++) {
      mHistRegistry->addClone("Centrality_0/", Form("Centrality_%u/", i));
    }
  }

  EventPlaneHelper *epHelp;

private:
  HistogramRegistry* mHistRegistry;

  static constexpr std::string_view mCentClasses[] = {
    "Centrality_0/",
    "Centrality_1/",
    "Centrality_2/",
    "Centrality_3/",
    "Centrality_4/",
    "Centrality_5/",
    "Centrality_6/",
    "Centrality_7/",
    "Centrality_8/",
    "Centrality_9/"};

  ClassDefNV(JEPFlowAnalysis, 1);
};


#endif // PWGCF_JCORRAN_CORE_JEPFLOWANALYSIS_H_