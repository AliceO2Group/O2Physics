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

// Task to produce a table joinable to the jet tables for MC Particle level event weights
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "PWGJE/DataModel/Jet.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

template <typename MCParticleLevelJetTable, typename MCParticleLevelWeightsTable>

struct JetEventWeightMCPTask {
  Produces<MCParticleLevelWeightsTable> mcParticleLevelWeightsTable;

  void processDummy(aod::JetMcCollisions const&)
  {
  }
  PROCESS_SWITCH(JetEventWeightMCPTask, processDummy, "Dummy process", true);

  void processMCParticleLevelEventWeight(MCParticleLevelJetTable const& jet, aod::JetMcCollisions const&)
  {
    mcParticleLevelWeightsTable(jet.globalIndex(), jet.mcCollision().weight());
  }
  PROCESS_SWITCH(JetEventWeightMCPTask, processMCParticleLevelEventWeight, "Fill event weight tables for particle level MC jets", false);
};

using ChargedMCJetsEventWeight = JetEventWeightMCPTask<aod::ChargedMCParticleLevelJet, aod::ChargedMCParticleLevelJetEventWeights>;
using NeutralMCJetsEventWeight = JetEventWeightMCPTask<aod::NeutralMCParticleLevelJet, aod::NeutralMCParticleLevelJetEventWeights>;
using FullMCJetsEventWeight = JetEventWeightMCPTask<aod::FullMCParticleLevelJet, aod::FullMCParticleLevelJetEventWeights>;
using D0ChargedMCJetsEventWeight = JetEventWeightMCPTask<aod::D0ChargedMCParticleLevelJet, aod::D0ChargedMCParticleLevelJetEventWeights>;
using DplusChargedMCJetsEventWeight = JetEventWeightMCPTask<aod::DplusChargedMCParticleLevelJet, aod::DplusChargedMCParticleLevelJetEventWeights>;
using LcChargedMCJetsEventWeight = JetEventWeightMCPTask<aod::LcChargedMCParticleLevelJet, aod::LcChargedMCParticleLevelJetEventWeights>;
using BplusChargedMCJetsEventWeight = JetEventWeightMCPTask<aod::BplusChargedMCParticleLevelJet, aod::BplusChargedMCParticleLevelJetEventWeights>;
using V0ChargedMCJetsEventWeight = JetEventWeightMCPTask<aod::V0ChargedMCParticleLevelJet, aod::V0ChargedMCParticleLevelJetEventWeights>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(
    adaptAnalysisTask<ChargedMCJetsEventWeight>(cfgc,
                                                SetDefaultProcesses{}, TaskName{"jet-eventweight-mcp-charged"}));

  tasks.emplace_back(
    adaptAnalysisTask<NeutralMCJetsEventWeight>(cfgc,
                                                SetDefaultProcesses{}, TaskName{"jet-eventweight-mcp-neutral"}));

  tasks.emplace_back(
    adaptAnalysisTask<FullMCJetsEventWeight>(cfgc,
                                             SetDefaultProcesses{}, TaskName{"jet-eventweight-mcp-full"}));

  tasks.emplace_back(
    adaptAnalysisTask<D0ChargedMCJetsEventWeight>(cfgc,
                                                  SetDefaultProcesses{}, TaskName{"jet-d0-eventweight-mcp-charged"}));

  tasks.emplace_back(
    adaptAnalysisTask<DplusChargedMCJetsEventWeight>(cfgc,
                                                     SetDefaultProcesses{}, TaskName{"jet-dplus-eventweight-mcp-charged"}));

  tasks.emplace_back(
    adaptAnalysisTask<LcChargedMCJetsEventWeight>(cfgc,
                                                  SetDefaultProcesses{}, TaskName{"jet-lc-eventweight-mcp-charged"}));

  tasks.emplace_back(
    adaptAnalysisTask<BplusChargedMCJetsEventWeight>(cfgc,
                                                     SetDefaultProcesses{}, TaskName{"jet-bplus-eventweight-mcp-charged"}));

  tasks.emplace_back(
    adaptAnalysisTask<V0ChargedMCJetsEventWeight>(cfgc,
                                                  SetDefaultProcesses{}, TaskName{"jet-v0-eventweight-mcp-charged"}));

  return WorkflowSpec{tasks};
}
