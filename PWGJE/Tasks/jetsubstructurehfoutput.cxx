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

// heavy-flavour jet substructure tree filling task (subscribing to jet finder hf and jet substructure hf tasks)
//
// Author: Nima Zardoshti
//

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructureHF.h"
#include "PWGJE/Core/JetFinder.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::analysis::hf_cuts_d0_to_pi_k;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec hfjetsubstructureoutputMode = {
    "hfjetsubstructureoutputMode",
    VariantType::String,
    "",
    {"HF jet substrcture output mode."},
  };
  workflowOptions.push_back(hfjetsubstructureoutputMode);
}

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

template <typename OutputTable>
struct JetSubstructureHFOutputTask {
  Produces<OutputTable> jetSubstructurehfoutputTable;

  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};

  void processData(soa::Join<aod::HFJets, aod::HFJetConstituents, aod::JetSubtructureHFData>::iterator const& jet, // add template back
                   soa::Join<aod::HfCand2Prong, aod::HfSelD0> const& candidates,
                   aod::Tracks const& tracks)
  {
    auto cands = jet.hfcandidates_as<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>();
    auto cand = cands[0];

    auto invMassCand = -1;
    auto invMassCandBar = -1;
    if (cand.isSelD0() >= selectionFlagD0) {
      invMassCand = invMassD0ToPiK(cand);
    }
    if (cand.isSelD0bar() >= selectionFlagD0bar) {
      invMassCandBar = invMassD0barToKPi(cand);
    }

    jetSubstructurehfoutputTable(jet.pt(), jet.phi(), jet.eta(), jet.tracks().size(), cand.pt(), cand.phi(), cand.eta(), yD0(cand), invMassCand, invMassCandBar, jet.zg(), jet.rg(), jet.nsd());
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processData, "HF jet substructure output on data", true);

  void processMCD(soa::Join<aod::MCDetectorLevelHFJets, aod::MCDetectorLevelHFJetConstituents, aod::JetSubtructureHFMCDet>::iterator const& jet,
                  soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec> const& candidates,
                  aod::Tracks const& tracks)
  {

    auto cands = jet.hfcandidates_as<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>>();
    auto cand = cands[0];
    auto invMassCand = -1;
    auto invMassCandBar = -1;
    if (cand.isSelD0() >= selectionFlagD0) {
      invMassCand = invMassD0ToPiK(cand);
    }
    if (cand.isSelD0bar() >= selectionFlagD0bar) {
      invMassCandBar = invMassD0barToKPi(cand);
    }
    jetSubstructurehfoutputTable(jet.pt(), jet.phi(), jet.eta(), jet.tracks().size(), cand.pt(), cand.phi(), cand.eta(), yD0(cand), invMassCand, invMassCandBar, jet.zg(), jet.rg(), jet.nsd());
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processMCD, "HF jet substructure output on MC detector level", true);

  void processMCP(soa::Join<aod::MCParticleLevelHFJets, aod::MCParticleLevelHFJetConstituents, aod::JetSubtructureHFMCGen>::iterator const& jet,
                  // aod::JetConstituentsSub const& constituentsSub,
                  aod::McParticles const& particles)
  {

    auto hfparticles = jet.hfcandidates_as<aod::McParticles>();
    auto hfparticle = hfparticles[0];
    auto Y = RecoDecay::y(array{hfparticle.px(), hfparticle.py(), hfparticle.pz()}, RecoDecay::getMassPDG(hfparticle.pdgCode()));
    auto M = RecoDecay::getMassPDG(hfparticle.pdgCode());
    jetSubstructurehfoutputTable(jet.pt(), jet.phi(), jet.eta(), jet.tracks().size(), hfparticle.pt(), hfparticle.phi(), hfparticle.eta(), Y, M, M, jet.zg(), jet.rg(), jet.nsd());
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processMCP, "HF jet substructure output on MC particle level", true);
};
using JetSubstructureHFOutputDataLevel = JetSubstructureHFOutputTask<aod::JetSubstructureHFOutputData>;
using JetSubstructureHFOutputMCParticleLevel = JetSubstructureHFOutputTask<aod::JetSubstructureHFOutputMCGen>;
using JetSubstructureHFOutputMCDetectorLevel = JetSubstructureHFOutputTask<aod::JetSubstructureHFOutputMCDet>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  auto hfjetsubstructureoutputMode = cfgc.options().get<std::string>("hfjetsubstructureoutputMode");

  if (hfjetsubstructureoutputMode.find("data") != std::string::npos || hfjetsubstructureoutputMode.empty())
    tasks.emplace_back(adaptAnalysisTask<JetSubstructureHFOutputDataLevel>(cfgc,
                                                                           SetDefaultProcesses{{{"processData", true}, {"processMCP", false}, {"processMCD", false}}},
                                                                           TaskName{"jet-substructure-hf-output-data"}));

  if (hfjetsubstructureoutputMode.find("mcp") != std::string::npos || hfjetsubstructureoutputMode.empty())
    tasks.emplace_back(adaptAnalysisTask<JetSubstructureHFOutputMCParticleLevel>(cfgc,
                                                                                 SetDefaultProcesses{{{"processData", false}, {"processMCP", true}, {"processMCD", false}}},
                                                                                 TaskName{"jet-substructure-hf-output-mcp"}));

  if (hfjetsubstructureoutputMode.find("mcd") != std::string::npos || hfjetsubstructureoutputMode.empty())
    tasks.emplace_back(adaptAnalysisTask<JetSubstructureHFOutputMCDetectorLevel>(cfgc,
                                                                                 SetDefaultProcesses{{{"processData", false}, {"processMCP", false}, {"processMCD", true}}},
                                                                                 TaskName{"jet-substructure-hf-output-mcd"}));

  return WorkflowSpec{tasks};
}
