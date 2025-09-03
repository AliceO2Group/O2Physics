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

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/runDataProcessing.h>

#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename MCDetectorLevelJetTable, typename MCDetectorLevelWeightsTable, typename MCDetectorLevelEventWiseSubtractedJetTable, typename MCDetectorLevelEventWiseSubtractedWeightsTable>

struct JetEventWeightMCDTask {
  Produces<MCDetectorLevelWeightsTable> mcDetectorLevelWeightsTable;
  Produces<MCDetectorLevelEventWiseSubtractedWeightsTable> mcDetectorLevelEventWiseSubtractedWeightsTable;

  void processDummy(aod::JetCollisions const&)
  {
  }
  PROCESS_SWITCH(JetEventWeightMCDTask, processDummy, "Dummy process", true);

  void processMCDetectorLevelEventWeight(MCDetectorLevelJetTable const& jet, soa::Join<aod::JetCollisions, aod::JMcCollisionLbs> const&, aod::JetMcCollisions const&)
  {
    auto collision = jet.template collision_as<soa::Join<aod::JetCollisions, aod::JMcCollisionLbs>>();
    mcDetectorLevelWeightsTable(jet.globalIndex(), collision.weight());
  }
  PROCESS_SWITCH(JetEventWeightMCDTask, processMCDetectorLevelEventWeight, "Fill event weight tables for detector level MC jets", false);

  void processMCDetectorLevelEventWiseSubtractedEventWeight(MCDetectorLevelEventWiseSubtractedJetTable const& jet, soa::Join<aod::JetCollisions, aod::JMcCollisionLbs> const&, aod::JetMcCollisions const&)
  {
    auto collision = jet.template collision_as<soa::Join<aod::JetCollisions, aod::JMcCollisionLbs>>();
    mcDetectorLevelEventWiseSubtractedWeightsTable(jet.globalIndex(), collision.weight());
  }
  PROCESS_SWITCH(JetEventWeightMCDTask, processMCDetectorLevelEventWiseSubtractedEventWeight, "Fill event weight tables for detector level MC jets", false);
};

using ChargedMCJetsEventWeight = JetEventWeightMCDTask<aod::ChargedMCDetectorLevelJet, aod::ChargedMCDetectorLevelJetEventWeights, aod::ChargedMCDetectorLevelEventWiseSubtractedJet, aod::ChargedMCDetectorLevelEventWiseSubtractedJetEventWeights>;
using NeutralMCJetsEventWeight = JetEventWeightMCDTask<aod::NeutralMCDetectorLevelJet, aod::NeutralMCDetectorLevelJetEventWeights, aod::NeutralMCDetectorLevelEventWiseSubtractedJet, aod::NeutralMCDetectorLevelEventWiseSubtractedJetEventWeights>;
using FullMCJetsEventWeight = JetEventWeightMCDTask<aod::FullMCDetectorLevelJet, aod::FullMCDetectorLevelJetEventWeights, aod::FullMCDetectorLevelEventWiseSubtractedJet, aod::FullMCDetectorLevelEventWiseSubtractedJetEventWeights>;
using D0ChargedMCJetsEventWeight = JetEventWeightMCDTask<aod::D0ChargedMCDetectorLevelJet, aod::D0ChargedMCDetectorLevelJetEventWeights, aod::D0ChargedMCDetectorLevelEventWiseSubtractedJet, aod::D0ChargedMCDetectorLevelEventWiseSubtractedJetEventWeights>;
using DplusChargedMCJetsEventWeight = JetEventWeightMCDTask<aod::DplusChargedMCDetectorLevelJet, aod::DplusChargedMCDetectorLevelJetEventWeights, aod::DplusChargedMCDetectorLevelEventWiseSubtractedJet, aod::DplusChargedMCDetectorLevelEventWiseSubtractedJetEventWeights>;
using DsChargedMCJetsEventWeight = JetEventWeightMCDTask<aod::DsChargedMCDetectorLevelJet, aod::DsChargedMCDetectorLevelJetEventWeights, aod::DsChargedMCDetectorLevelEventWiseSubtractedJet, aod::DsChargedMCDetectorLevelEventWiseSubtractedJetEventWeights>;
using DstarChargedMCJetsEventWeight = JetEventWeightMCDTask<aod::DstarChargedMCDetectorLevelJet, aod::DstarChargedMCDetectorLevelJetEventWeights, aod::DstarChargedMCDetectorLevelEventWiseSubtractedJet, aod::DstarChargedMCDetectorLevelEventWiseSubtractedJetEventWeights>;
using LcChargedMCJetsEventWeight = JetEventWeightMCDTask<aod::LcChargedMCDetectorLevelJet, aod::LcChargedMCDetectorLevelJetEventWeights, aod::LcChargedMCDetectorLevelEventWiseSubtractedJet, aod::LcChargedMCDetectorLevelEventWiseSubtractedJetEventWeights>;
using B0ChargedMCJetsEventWeight = JetEventWeightMCDTask<aod::B0ChargedMCDetectorLevelJet, aod::B0ChargedMCDetectorLevelJetEventWeights, aod::B0ChargedMCDetectorLevelEventWiseSubtractedJet, aod::B0ChargedMCDetectorLevelEventWiseSubtractedJetEventWeights>;
using BplusChargedMCJetsEventWeight = JetEventWeightMCDTask<aod::BplusChargedMCDetectorLevelJet, aod::BplusChargedMCDetectorLevelJetEventWeights, aod::BplusChargedMCDetectorLevelEventWiseSubtractedJet, aod::BplusChargedMCDetectorLevelEventWiseSubtractedJetEventWeights>;
using XicToXiPiPiChargedMCJetsEventWeight = JetEventWeightMCDTask<aod::XicToXiPiPiChargedMCDetectorLevelJet, aod::XicToXiPiPiChargedMCDetectorLevelJetEventWeights, aod::XicToXiPiPiChargedMCDetectorLevelEventWiseSubtractedJet, aod::XicToXiPiPiChargedMCDetectorLevelEventWiseSubtractedJetEventWeights>;
using V0ChargedMCJetsEventWeight = JetEventWeightMCDTask<aod::V0ChargedMCDetectorLevelJet, aod::V0ChargedMCDetectorLevelJetEventWeights, aod::V0ChargedMCDetectorLevelEventWiseSubtractedJet, aod::V0ChargedMCDetectorLevelEventWiseSubtractedJetEventWeights>;

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
    adaptAnalysisTask<DplusChargedMCJetsEventWeight>(cfgc,
                                                     SetDefaultProcesses{}, TaskName{"jet-dplus-eventweight-mcd-charged"}));

  tasks.emplace_back(
    adaptAnalysisTask<DsChargedMCJetsEventWeight>(cfgc,
                                                  SetDefaultProcesses{}, TaskName{"jet-ds-eventweight-mcd-charged"}));

  tasks.emplace_back(
    adaptAnalysisTask<DstarChargedMCJetsEventWeight>(cfgc,
                                                     SetDefaultProcesses{}, TaskName{"jet-dstar-eventweight-mcd-charged"}));

  tasks.emplace_back(
    adaptAnalysisTask<LcChargedMCJetsEventWeight>(cfgc,
                                                  SetDefaultProcesses{}, TaskName{"jet-lc-eventweight-mcd-charged"}));

  tasks.emplace_back(
    adaptAnalysisTask<B0ChargedMCJetsEventWeight>(cfgc,
                                                  SetDefaultProcesses{}, TaskName{"jet-b0-eventweight-mcd-charged"}));

  tasks.emplace_back(
    adaptAnalysisTask<BplusChargedMCJetsEventWeight>(cfgc,
                                                     SetDefaultProcesses{}, TaskName{"jet-bplus-eventweight-mcd-charged"}));

  tasks.emplace_back(
    adaptAnalysisTask<XicToXiPiPiChargedMCJetsEventWeight>(cfgc,
                                                           SetDefaultProcesses{}, TaskName{"jet-xictoxipipi-eventweight-mcd-charged"}));

  tasks.emplace_back(
    adaptAnalysisTask<V0ChargedMCJetsEventWeight>(cfgc,
                                                  SetDefaultProcesses{}, TaskName{"jet-v0-eventweight-mcd-charged"}));

  return WorkflowSpec{tasks};
}
