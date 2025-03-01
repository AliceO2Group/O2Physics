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

// jet substructure output D+ charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <vector>
#include "PWGJE/Tasks/jetSubstructureHFOutput.cxx"

using JetSubstructureOutputDplus = JetSubstructureHFOutputTask<aod::CollisionsDplus, soa::Join<aod::McCollisionsDplus, aod::HfDplusMcRCollIds>, aod::CandidatesDplusData, aod::CandidatesDplusMCD, aod::CandidatesDplusMCP, aod::BkgDplusRhos, aod::BkgDplusMcRhos, aod::JTrackDplusSubs, soa::Join<aod::DplusChargedJets, aod::DplusChargedJetConstituents, aod::DplusCJetSSs>, soa::Join<aod::DplusChargedJets, aod::DplusChargedJetConstituents, aod::DplusCJetSSs, aod::DplusChargedJetsMatchedToDplusChargedEventWiseSubtractedJets>, soa::Join<aod::DplusChargedSPs, aod::DplusChargedSPsMatchedToDplusChargedEventWiseSubtractedSPs>, soa::Join<aod::DplusChargedPRs, aod::DplusChargedPRsMatchedToDplusChargedEventWiseSubtractedPRs>, aod::DplusCJetCOs, aod::DplusCJetOs, aod::DplusCJetSSOs, aod::DplusCJetMOs, soa::Join<aod::DplusChargedMCDetectorLevelJets, aod::DplusChargedMCDetectorLevelJetConstituents, aod::DplusCMCDJetSSs, aod::DplusChargedMCDetectorLevelJetsMatchedToDplusChargedMCParticleLevelJets>, soa::Join<aod::DplusChargedMCDetectorLevelSPs, aod::DplusChargedMCDetectorLevelSPsMatchedToDplusChargedMCParticleLevelSPs>, soa::Join<aod::DplusChargedMCDetectorLevelPRs, aod::DplusChargedMCDetectorLevelPRsMatchedToDplusChargedMCParticleLevelPRs>, aod::DplusCMCDJetCOs, aod::DplusCMCDJetOs, aod::DplusCMCDJetSSOs, aod::DplusCMCDJetMOs, soa::Join<aod::DplusChargedMCParticleLevelJets, aod::DplusChargedMCParticleLevelJetConstituents, aod::DplusCMCPJetSSs>, soa::Join<aod::DplusChargedMCParticleLevelJets, aod::DplusChargedMCParticleLevelJetConstituents, aod::DplusCMCPJetSSs, aod::DplusChargedMCParticleLevelJetsMatchedToDplusChargedMCDetectorLevelJets>, soa::Join<aod::DplusChargedMCParticleLevelSPs, aod::DplusChargedMCParticleLevelSPsMatchedToDplusChargedMCDetectorLevelSPs>, soa::Join<aod::DplusChargedMCParticleLevelPRs, aod::DplusChargedMCParticleLevelPRsMatchedToDplusChargedMCDetectorLevelPRs>, aod::DplusCMCPJetCOs, aod::DplusCMCPJetOs, aod::DplusCMCPJetSSOs, aod::DplusCMCPJetMOs, soa::Join<aod::DplusChargedEventWiseSubtractedJets, aod::DplusChargedEventWiseSubtractedJetConstituents, aod::DplusCEWSJetSSs, aod::DplusChargedEventWiseSubtractedJetsMatchedToDplusChargedJets>, soa::Join<aod::DplusChargedEventWiseSubtractedSPs, aod::DplusChargedEventWiseSubtractedSPsMatchedToDplusChargedSPs>, soa::Join<aod::DplusChargedEventWiseSubtractedPRs, aod::DplusChargedEventWiseSubtractedPRsMatchedToDplusChargedPRs>, aod::DplusCEWSJetCOs, aod::DplusCEWSJetOs, aod::DplusCEWSJetSSOs, aod::DplusCEWSJetMOs, aod::StoredHfDplusCollBase, aod::StoredHfDplusBases, aod::StoredHfDplusPars, aod::StoredHfDplusParEs, aod::JDumDplusParDaus, aod::StoredHfDplusSels, aod::StoredHfDplusMls, aod::JDumDplusMlDaus, aod::StoredHfDplusMcs, aod::StoredHfDplusMcCollBases, aod::StoredHfDplusMcRCollIds, aod::StoredHfDplusPBases>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputDplus>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructure-dplus-output"}));

  return WorkflowSpec{tasks};
}
