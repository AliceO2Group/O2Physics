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
using namespace o2::aod::hf_cand_2prong;
using namespace o2::analysis::hf_cuts_d0_to_pi_k;

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

template <typename OutputTable, typename SubstructureOutputTable>
struct JetSubstructureHFOutputTask {
  Produces<OutputTable> jetOutputTable;
  Produces<SubstructureOutputTable> jetSubstructureOutputTable;

  template <typename T, typename U>
  void fillTables(T const& jet, U const& cands)
  {
    auto cand = cands[0];
    jetOutputTable(jet.globalIndex(), cand.globalIndex(), jet.pt(), jet.phi(), jet.eta(), jet.tracks().size());
    jetSubstructureOutputTable(jet.globalIndex(), jet.zg(), jet.rg(), jet.nsd());
  }

  void processDummy(aod::Tracks const& track)
  {
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processDummy, "Dummy process function turned on by default", true);

  void processD0Data(soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents, aod::D0ChargedJetSubstructures>::iterator const& jet, // add template back
                     soa::Join<aod::HfCand2Prong, aod::HfSelD0> const& candidates, aod::Tracks const& tracks)
  {
    auto cands = jet.hfcandidates_as<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>();
    fillTables(jet, cands);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processD0Data, "D0 jet substructure output on data", false);

  void processD0MCD(soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetSubstructures>::iterator const& jet, soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec> const& candidates, aod::Tracks const& tracks)
  {
    auto cands = jet.hfcandidates_as<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>>();
    fillTables(jet, cands);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processD0MCD, "D0 jet substructure output on MC detector level", false);

  void processD0MCP(soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelJetSubstructures>::iterator const& jet, aod::McParticles const& particles)
  {
    auto hfparticles = jet.hfcandidates_as<aod::McParticles>();
    fillTables(jet, hfparticles);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processD0MCP, "D0 jet substructure output on MC particle level", false);

  void processLcData(soa::Join<aod::LcChargedJets, aod::LcChargedJetConstituents, aod::LcChargedJetSubstructures>::iterator const& jet, // add template back
                     soa::Join<aod::HfCand3Prong, aod::HfSelLc> const& candidates, aod::Tracks const& tracks)
  {
    auto cands = jet.hfcandidates_as<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>();
    fillTables(jet, cands);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processLcData, "Lc jet substructure output on data", false);

  void processLcMCD(soa::Join<aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents, aod::LcChargedMCDetectorLevelJetSubstructures>::iterator const& jet, soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec> const& candidates, aod::Tracks const& tracks)
  {
    auto cands = jet.hfcandidates_as<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>>();
    fillTables(jet, cands);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processLcMCD, "Lc jet substructure output on MC detector level", false);

  void processLcMCP(soa::Join<aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents, aod::LcChargedMCParticleLevelJetSubstructures>::iterator const& jet, aod::McParticles const& particles)
  {
    auto hfparticles = jet.hfcandidates_as<aod::McParticles>();
    fillTables(jet, hfparticles);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processLcMCP, "Lc jet substructure output on MC particle level", false);

  void processBplusData(soa::Join<aod::BplusChargedJets, aod::BplusChargedJetConstituents, aod::BplusChargedJetSubstructures>::iterator const& jet, // add template back
                        soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi> const& candidates, aod::Tracks const& tracks)
  {
    auto cands = jet.hfcandidates_as<soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi>>();
    fillTables(jet, cands);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processBplusData, "B+ jet substructure output on data", false);

  void processBplusMCD(soa::Join<aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents, aod::BplusChargedMCDetectorLevelJetSubstructures>::iterator const& jet, soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi, aod::HfCandBplusMcRec> const& candidates, aod::Tracks const& tracks)
  {
    auto cands = jet.hfcandidates_as<soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi, aod::HfCandBplusMcRec>>();
    fillTables(jet, cands);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processBplusMCD, "B+ jet substructure output on MC detector level", false);

  void processBplusMCP(soa::Join<aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents, aod::BplusChargedMCParticleLevelJetSubstructures>::iterator const& jet, aod::McParticles const& particles)
  {
    auto hfparticles = jet.hfcandidates_as<aod::McParticles>();
    fillTables(jet, hfparticles);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processBplusMCP, "B+ jet substructure output on MC particle level", false);
};
using JetSubstructureOutputDataD0 = JetSubstructureHFOutputTask<aod::D0ChargedJetOutput, aod::D0ChargedJetSubstructureOutput>;
using JetSubstructureOutputMCDetectorLevelD0 = JetSubstructureHFOutputTask<aod::D0ChargedMCDetectorLevelJetOutput, aod::D0ChargedMCDetectorLevelJetSubstructureOutput>;
using JetSubstructureOutputMCParticleLevelD0 = JetSubstructureHFOutputTask<aod::D0ChargedMCParticleLevelJetOutput, aod::D0ChargedMCParticleLevelJetSubstructureOutput>;
using JetSubstructureOutputDataLc = JetSubstructureHFOutputTask<aod::LcChargedJetOutput, aod::LcChargedJetSubstructureOutput>;
using JetSubstructureOutputMCDetectorLevelLc = JetSubstructureHFOutputTask<aod::LcChargedMCDetectorLevelJetOutput, aod::LcChargedMCDetectorLevelJetSubstructureOutput>;
using JetSubstructureOutputMCParticleLevelLc = JetSubstructureHFOutputTask<aod::LcChargedMCParticleLevelJetOutput, aod::LcChargedMCParticleLevelJetSubstructureOutput>;
using JetSubstructureOutputDataBplus = JetSubstructureHFOutputTask<aod::BplusChargedJetOutput, aod::BplusChargedJetSubstructureOutput>;
using JetSubstructureOutputMCDetectorLevelBplus = JetSubstructureHFOutputTask<aod::BplusChargedMCDetectorLevelJetOutput, aod::BplusChargedMCDetectorLevelJetSubstructureOutput>;
using JetSubstructureOutputMCParticleLevelBplus = JetSubstructureHFOutputTask<aod::BplusChargedMCParticleLevelJetOutput, aod::BplusChargedMCParticleLevelJetSubstructureOutput>;

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

  return WorkflowSpec{tasks};
}
