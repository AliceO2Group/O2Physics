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
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/Core/JetFinder.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

template <typename CollisionTable, typename JetTable, typename CandidateTable, typename OutputTable, typename SubstructureOutputTable>
struct JetSubstructureHFOutputTask {
  Produces<OutputTable> jetOutputTable;
  Produces<SubstructureOutputTable> jetSubstructureOutputTable;

  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT cut"};
  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};

  std::vector<double> jetRadiiValues;

  void init(InitContext const&)
  {
    jetRadiiValues = (std::vector<double>)jetRadii;
  }

  Filter jetSelection = aod::jet::pt >= jetPtMin;

  template <typename T, typename U, typename V>
  void fillTables(T const& collision, U const& jet, V const& cands)
  {
    auto cand = cands[0];
    jetOutputTable(collision.globalIndex(), jet.globalIndex(), cand.globalIndex(), jet.pt(), jet.phi(), jet.eta(), jet.r(), jet.tracks().size());
    jetSubstructureOutputTable(jet.globalIndex(), jet.zg(), jet.rg(), jet.nsd());
  }

  void processDummy(aod::Collision const& collision)
  {
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processDummy, "Dummy process function turned on by default", true);

  void processOutput(typename CollisionTable::iterator const& collision, JetTable const& jets, // add template back
                     CandidateTable const& candidates, aod::Tracks const& tracks)
  {
    for (const auto& jet : jets) {
      for (const auto& jetRadiiValue : jetRadiiValues) {
        if (jet.r() == round(jetRadiiValue * 100.0f)) {
          auto cands = jet.template hfcandidates_as<CandidateTable>();
          fillTables(collision, jet, cands);
        }
      }
    }
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutput, "HF jet substructure output on data", false);
};
using JetSubstructureOutputDataD0 = JetSubstructureHFOutputTask<aod::Collisions, soa::Filtered<soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents, aod::D0ChargedJetSubstructures>>, soa::Join<aod::HfCand2Prong, aod::HfSelD0>, aod::D0ChargedJetOutput, aod::D0ChargedJetSubstructureOutput>;
using JetSubstructureOutputMCDetectorLevelD0 = JetSubstructureHFOutputTask<aod::Collisions, soa::Filtered<soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetSubstructures>>, soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>, aod::D0ChargedMCDetectorLevelJetOutput, aod::D0ChargedMCDetectorLevelJetSubstructureOutput>;
using JetSubstructureOutputMCParticleLevelD0 = JetSubstructureHFOutputTask<aod::McCollisions, soa::Filtered<soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelJetSubstructures>>, aod::McParticles, aod::D0ChargedMCParticleLevelJetOutput, aod::D0ChargedMCParticleLevelJetSubstructureOutput>;

// using JetSubstructureOutputDataLc = JetSubstructureHFOutputTask<aod::Collisions, soa::Filtered<soa::Join<aod::LcChargedJets, aod::LcChargedJetConstituents, aod::LcChargedJetSubstructures>>, soa::Join<aod::HfCand3Prong, aod::HfSelLc>, aod::LcChargedJetOutput, aod::LcChargedJetSubstructureOutput>;
//  using JetSubstructureOutputMCDetectorLevelLc = JetSubstructureHFOutputTask<aod::Collisions, soa::Filtered<soa::Join<aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents, aod::LcChargedMCDetectorLevelJetSubstructures>>, soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>,  aod::LcChargedMCDetectorLevelJetOutput, aod::LcChargedMCDetectorLevelJetSubstructureOutput>;
//  using JetSubstructureOutputMCParticleLevelLc = JetSubstructureHFOutputTask<aod::McCollisions, soa::Filtered<soa::Join<aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents, aod::LcChargedMCParticleLevelJetSubstructures>>, aod::McParticles, aod::LcChargedMCParticleLevelJetOutput, aod::LcChargedMCParticleLevelJetSubstructureOutput>;

// using JetSubstructureOutputDataBplus = JetSubstructureHFOutputTask<aod::Collisions, soa::Filtered<soa::Join<aod::BplusChargedJets, aod::BplusChargedJetConstituents, aod::BplusChargedJetSubstructures>> , soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi>, aod::BplusChargedJetOutput, aod::BplusChargedJetSubstructureOutput>;
// using JetSubstructureOutputMCDetectorLevelBplus = JetSubstructureHFOutputTask<aod::Collisions, soa::Filtered<soa::Join<aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents, aod::BplusChargedMCDetectorLevelJetSubstructures>>, soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi, aod::HfCandBplusMcRec>,  aod::BplusChargedMCDetectorLevelJetOutput, aod::BplusChargedMCDetectorLevelJetSubstructureOutput>;
// using JetSubstructureOutputMCParticleLevelBplus = JetSubstructureHFOutputTask<aod::McCollisions, soa::Filtered<soa::Join<aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents, aod::BplusChargedMCParticleLevelJetSubstructures>>, aod::McParticles, aod::BplusChargedMCParticleLevelJetOutput, aod::BplusChargedMCParticleLevelJetSubstructureOutput>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputDataD0>(cfgc,
                                                                    SetDefaultProcesses{},
                                                                    TaskName{"jet-substructure-output-D0-data"}));

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputMCDetectorLevelD0>(cfgc,
                                                                               SetDefaultProcesses{},
                                                                               TaskName{"jet-substructure-output-D0-mcd"}));

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputMCParticleLevelD0>(cfgc,
                                                                               SetDefaultProcesses{},
                                                                               TaskName{"jet-substructure-output-D0-mcp"}));
  /*
   tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputDataLc>(cfgc,
                                                                     SetDefaultProcesses{},
                                                                     TaskName{"jet-substructure-output-Lc-data"}));

     tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputMCDetectorLevelLc>(cfgc,
                                                                                  SetDefaultProcesses{},
                                                                                  TaskName{"jet-substructure-output-Lc-mcd"}));

     tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputMCParticleLevelLc>(cfgc,
                                                                                  SetDefaultProcesses{},
                                                                                  TaskName{"jet-substructure-output-Lc-mcp"}));

     tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputDataBplus>(cfgc,
                                                                          SetDefaultProcesses{},
                                                                          TaskName{"jet-substructure-output-Bplus-data"}));

     tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputMCDetectorLevelBplus>(cfgc,
                                                                                     SetDefaultProcesses{},
                                                                                     TaskName{"jet-substructure-output-Bplus-mcd"}));

     tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputMCParticleLevelBplus>(cfgc,
                                                                                     SetDefaultProcesses{},
                                                                                     TaskName{"jet-substructure-output-Bplus-mcp"}));
   */
  return WorkflowSpec{tasks};
}
