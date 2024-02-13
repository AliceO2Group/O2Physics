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

// Task to produce a table joinable to the jet tables for MC Detector level event weights
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

template <typename MCDetectorLevelJetTable, typename MCDetectorLevelWeightsTable, typename MCDetectorLevelEventWiseSubtractedJetTable, typename MCDetectorLevelEventWiseSubtractedWeightsTable>

struct JetEventWeightMCDTask {
  Produces<MCDetectorLevelWeightsTable> mcDetectorLevelWeightsTable;
  Produces<MCDetectorLevelEventWiseSubtractedWeightsTable> mcDetectorLevelEventWiseSubtractedWeightsTable;

  void processDummy(JetCollisions const& collisions)
  {
  }
  PROCESS_SWITCH(JetEventWeightMCDTask, processDummy, "Dummy process", true);

  void processMCDetectorLevelEventWeight(MCDetectorLevelJetTable const& jet, soa::Join<JetCollisions, aod::JMcCollisionLbs> const& collisions, JetMcCollisions const& mcCollisions)
  {
    auto collision = jet.template collision_as<soa::Join<JetCollisions, aod::JMcCollisionLbs>>();
    mcDetectorLevelWeightsTable(jet.globalIndex(), collision.mcCollision().weight());
  }
  PROCESS_SWITCH(JetEventWeightMCDTask, processMCDetectorLevelEventWeight, "Fill event weight tables for detector level MC jets", false);

  void processMCDetectorLevelEventWiseSubtractedEventWeight(MCDetectorLevelEventWiseSubtractedJetTable const& jet, soa::Join<JetCollisions, aod::JMcCollisionLbs> const& collisions, JetMcCollisions const& mcCollisions)
  {
    auto collision = jet.template collision_as<soa::Join<JetCollisions, aod::JMcCollisionLbs>>();
    mcDetectorLevelEventWiseSubtractedWeightsTable(jet.globalIndex(), collision.mcCollision().weight());
  }
  PROCESS_SWITCH(JetEventWeightMCDTask, processMCDetectorLevelEventWiseSubtractedEventWeight, "Fill event weight tables for detector level MC jets", false);
};

using ChargedMCJetsEventWeight = JetEventWeightMCDTask<aod::ChargedMCDetectorLevelJet, aod::ChargedMCDetectorLevelJetEventWeights, aod::ChargedMCDetectorLevelEventWiseSubtractedJet, aod::ChargedMCDetectorLevelEventWiseSubtractedJetEventWeights>;
using NeutralMCJetsEventWeight = JetEventWeightMCDTask<aod::NeutralMCDetectorLevelJet, aod::NeutralMCDetectorLevelJetEventWeights, aod::NeutralMCDetectorLevelEventWiseSubtractedJet, aod::NeutralMCDetectorLevelEventWiseSubtractedJetEventWeights>;
using FullMCJetsEventWeight = JetEventWeightMCDTask<aod::FullMCDetectorLevelJet, aod::FullMCDetectorLevelJetEventWeights, aod::FullMCDetectorLevelEventWiseSubtractedJet, aod::FullMCDetectorLevelEventWiseSubtractedJetEventWeights>;
using D0ChargedMCJetsEventWeight = JetEventWeightMCDTask<aod::D0ChargedMCDetectorLevelJet, aod::D0ChargedMCDetectorLevelJetEventWeights, aod::D0ChargedMCDetectorLevelEventWiseSubtractedJet, aod::D0ChargedMCDetectorLevelEventWiseSubtractedJetEventWeights>;
using LcChargedMCJetsEventWeight = JetEventWeightMCDTask<aod::LcChargedMCDetectorLevelJet, aod::LcChargedMCDetectorLevelJetEventWeights, aod::LcChargedMCDetectorLevelEventWiseSubtractedJet, aod::LcChargedMCDetectorLevelEventWiseSubtractedJetEventWeights>;
using BplusChargedMCJetsEventWeight = JetEventWeightMCDTask<aod::BplusChargedMCDetectorLevelJet, aod::BplusChargedMCDetectorLevelJetEventWeights, aod::BplusChargedMCDetectorLevelEventWiseSubtractedJet, aod::BplusChargedMCDetectorLevelEventWiseSubtractedJetEventWeights>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(
    adaptAnalysisTask<ChargedMCJetsEventWeight>(cfgc,
                                                SetDefaultProcesses{}, TaskName{"jet-eventweight-mcd-charged"}));

  tasks.emplace_back(
    adaptAnalysisTask<NeutralMCJetsEventWeight>(cfgc,
                                                SetDefaultProcesses{}, TaskName{"jet-eventweight-mcd-neutral"}));

  tasks.emplace_back(
    adaptAnalysisTask<FullMCJetsEventWeight>(cfgc,
                                             SetDefaultProcesses{}, TaskName{"jet-eventweight-mcd-full"}));

  tasks.emplace_back(
    adaptAnalysisTask<D0ChargedMCJetsEventWeight>(cfgc,
                                                  SetDefaultProcesses{}, TaskName{"jet-d0-eventweight-mcd-charged"}));

  tasks.emplace_back(
    adaptAnalysisTask<LcChargedMCJetsEventWeight>(cfgc,
                                                  SetDefaultProcesses{}, TaskName{"jet-lc-eventweight-mcd-charged"}));

  tasks.emplace_back(
    adaptAnalysisTask<BplusChargedMCJetsEventWeight>(cfgc,
                                                     SetDefaultProcesses{}, TaskName{"jet-bplus-eventweight-mcd-charged"}));

  return WorkflowSpec{tasks};
}
