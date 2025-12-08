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

// jet substructure output B+ charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGHF/DataModel/DerivedTables.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/DataModel/JetSubtraction.h"
#include "PWGJE/Tasks/jetSubstructureHFOutput.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/runDataProcessing.h>

#include <vector>

using JetSubstructureOutputBplus = JetSubstructureHFOutputTask<aod::CollisionsBplus, soa::Join<aod::McCollisionsBplus, aod::HfBplusMcRCollIds>, aod::McCollisionsBplus, aod::CandidatesBplusData, aod::CandidatesBplusMCD, aod::CandidatesBplusMCP, aod::BkgBplusRhos, aod::BkgBplusMcRhos, aod::JTrackBplusSubs, soa::Join<aod::BplusChargedJets, aod::BplusChargedJetConstituents, aod::BplusCJetSSs>, soa::Join<aod::BplusChargedJets, aod::BplusChargedJetConstituents, aod::BplusCJetSSs, aod::BplusChargedJetsMatchedToBplusChargedEventWiseSubtractedJets>, soa::Join<aod::BplusChargedSPs, aod::BplusChargedSPsMatchedToBplusChargedEventWiseSubtractedSPs>, soa::Join<aod::BplusChargedPRs, aod::BplusChargedPRsMatchedToBplusChargedEventWiseSubtractedPRs>, aod::BplusCJetCOs, aod::BplusCJetOs, aod::BplusCJetSSOs, aod::BplusCJetMOs, soa::Join<aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents, aod::BplusCMCDJetSSs, aod::BplusChargedMCDetectorLevelJetsMatchedToBplusChargedMCParticleLevelJets>, soa::Join<aod::BplusChargedMCDetectorLevelSPs, aod::BplusChargedMCDetectorLevelSPsMatchedToBplusChargedMCParticleLevelSPs>, soa::Join<aod::BplusChargedMCDetectorLevelPRs, aod::BplusChargedMCDetectorLevelPRsMatchedToBplusChargedMCParticleLevelPRs>, aod::BplusCMCDJetCOs, aod::BplusCMCDJetOs, aod::BplusCMCDJetSSOs, aod::BplusCMCDJetMOs, soa::Join<aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents, aod::BplusCMCPJetSSs>, soa::Join<aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents, aod::BplusCMCPJetSSs, aod::BplusChargedMCParticleLevelJetsMatchedToBplusChargedMCDetectorLevelJets>, soa::Join<aod::BplusChargedMCParticleLevelSPs, aod::BplusChargedMCParticleLevelSPsMatchedToBplusChargedMCDetectorLevelSPs>, soa::Join<aod::BplusChargedMCParticleLevelPRs, aod::BplusChargedMCParticleLevelPRsMatchedToBplusChargedMCDetectorLevelPRs>, aod::BplusCMCPJetCOs, aod::BplusCMCPJetMCCOs, aod::BplusCMCPJetOs, aod::BplusCMCPJetSSOs, aod::BplusCMCPJetMOs, soa::Join<aod::BplusChargedEventWiseSubtractedJets, aod::BplusChargedEventWiseSubtractedJetConstituents, aod::BplusCEWSJetSSs, aod::BplusChargedEventWiseSubtractedJetsMatchedToBplusChargedJets>, soa::Join<aod::BplusChargedEventWiseSubtractedSPs, aod::BplusChargedEventWiseSubtractedSPsMatchedToBplusChargedSPs>, soa::Join<aod::BplusChargedEventWiseSubtractedPRs, aod::BplusChargedEventWiseSubtractedPRsMatchedToBplusChargedPRs>, aod::BplusCEWSJetCOs, aod::BplusCEWSJetOs, aod::BplusCEWSJetSSOs, aod::BplusCEWSJetMOs, aod::StoredHfBplusCollBase, aod::StoredHfBplusBases, aod::StoredHfBplusPars, aod::StoredHfBplusParEs, aod::StoredHfBplusParD0s, aod::StoredHfBplusSels, aod::StoredHfBplusMls, aod::StoredHfBplusMlD0s, aod::StoredHfBplusMcs, aod::StoredHfBplusMcCollBases, aod::StoredHfBplusMcRCollIds, aod::StoredHfBplusPBases>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputBplus>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructure-bplus-output"}));

  return WorkflowSpec{tasks};
}
