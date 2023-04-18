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
//
//
// Analysis task to produce smeared pt,eta,phi for electrons/muons in dilepton analysis
//    Please write to: daiki.sekihata@cern.ch

#include <iostream>
#include <vector>
#include <TMath.h>
#include <TH1F.h>
#include <THashList.h>
#include <TString.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include "Framework/DataTypes.h"
#include "Framework/HistogramRegistry.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/CCDB/TriggerAliases.h"

using std::cout;
using std::endl;
using std::string;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

// Some definitions
namespace o2::aod
{

namespace emsmearedtrack
{
DECLARE_SOA_COLUMN(PtSmeared, ptSmeared, float);
DECLARE_SOA_COLUMN(EtaSmeared, etaSmeared, float);
DECLARE_SOA_COLUMN(PhiSmeared, phiSmeared, float);
} // namespace emsmearedtrack

DECLARE_SOA_TABLE(SmearedTracks, "AOD", "SMEAREDTRACK",
                  emsmearedtrack::PtSmeared, emsmearedtrack::EtaSmeared, emsmearedtrack::PhiSmeared);
} // namespace o2::aod

struct ApplySmearing {
  Produces<aod::SmearedTracks> smearedtrack;

  template <typename TTracksMC>
  void applySmearing(TTracksMC const& tracksMC)
  {
    for (auto& mctrack : tracksMC) {
      float ptgen = mctrack.pt();
      float etagen = mctrack.eta();
      float phigen = mctrack.phi();

      if (abs(mctrack.pdgCode()) == 11 || abs(mctrack.pdgCode()) == 13) {
        // apply smearing for electrons or muons.
        float ptsmeared = ptgen;
        float etasmeared = etagen;
        float phismeared = phigen;

        smearedtrack(ptsmeared, etasmeared, phismeared);
      } else {
        // don't touch
        smearedtrack(ptgen, etagen, phigen);
      }
    }
  }

  void processMCanalysis(ReducedMCTracks const& tracksMC)
  {
    applySmearing(tracksMC);
  }

  void processCocktail(aod::McParticles_001 const& tracksMC)
  {
    applySmearing(tracksMC);
  }

  void processDummy(aod::McParticles_001 const& tracksMC) {}

  PROCESS_SWITCH(ApplySmearing, processMCanalysis, "Run for MC analysis", false);
  PROCESS_SWITCH(ApplySmearing, processCocktail, "Run for cocktail analysis", false);
  PROCESS_SWITCH(ApplySmearing, processDummy, "Dummy process function", true);
};

struct CheckSmearing {
  using MyTracks = soa::Join<ReducedMCTracks, SmearedTracks>;
  HistogramRegistry registry{
    "registry",
    {
      {"hCorrelation_Pt", "pT correlation", {HistType::kTH2F, {{1000, 0.0f, 10.0f}, {1000, 0.0f, 10.0f}}}},
      {"hCorrelation_Eta", "eta correlation", {HistType::kTH2F, {{200, -1.0f, +1.0f}, {200, -1.0f, +1.0f}}}},
      {"hCorrelation_Phi", "phi correlation", {HistType::kTH2F, {{100, 0.0f, TMath::TwoPi()}, {100, 0.0f, TMath::TwoPi()}}}},
    },
  };

  void processCheck(MyTracks const& tracksMC)
  {
    for (auto& mctrack : tracksMC) {
      if (abs(mctrack.pdgCode()) != 11 && abs(mctrack.pdgCode()) != 13) {
        continue;
      }
      registry.fill(HIST("hCorrelation_Pt"), mctrack.pt(), mctrack.ptSmeared());
      registry.fill(HIST("hCorrelation_Eta"), mctrack.eta(), mctrack.etaSmeared());
      registry.fill(HIST("hCorrelation_Phi"), mctrack.phi(), mctrack.phiSmeared());
    } // end of mctrack loop
  }

  void processDummy(MyTracks const& tracksMC) {}

  PROCESS_SWITCH(CheckSmearing, processCheck, "Run for MC analysis", false);
  PROCESS_SWITCH(CheckSmearing, processDummy, "Dummy process function", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ApplySmearing>(cfgc, TaskName{"apply-smearing"}),
    adaptAnalysisTask<CheckSmearing>(cfgc, TaskName{"check-smearing"})};
}
