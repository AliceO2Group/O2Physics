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

#include "PWGJE/Tasks/jetSubstructureHFOutput.cxx"

using JetSubstructureOutputBplus = JetSubstructureHFOutputTask<aod::CollisionsBplus, soa::Join<aod::McCollisionsBplus, aod::HfBplusMcRCollIds>, aod::CandidatesBplusData, aod::CandidatesBplusMCD, aod::CandidatesBplusMCP, aod::JTrackBplusSubs, soa::Join<aod::BplusChargedJets, aod::BplusChargedJetConstituents, aod::BplusCJetSSs>, soa::Join<aod::BplusChargedJets, aod::BplusChargedJetConstituents, aod::BplusCJetSSs, aod::BplusChargedJetsMatchedToBplusChargedEventWiseSubtractedJets>, aod::BplusCJetCOs, aod::BplusCJetOs, aod::BplusCJetSSOs, aod::BplusCJetMOs, soa::Join<aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents, aod::BplusCMCDJetSSs, aod::BplusChargedMCDetectorLevelJetsMatchedToBplusChargedMCParticleLevelJets>, aod::BplusCMCDJetCOs, aod::BplusCMCDJetOs, aod::BplusCMCDJetSSOs, aod::BplusCMCDJetMOs, soa::Join<aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents, aod::BplusCMCPJetSSs, aod::BplusChargedMCParticleLevelJetsMatchedToBplusChargedMCDetectorLevelJets>, aod::BplusCMCPJetCOs, aod::BplusCMCPJetOs, aod::BplusCMCPJetSSOs, aod::BplusCMCPJetMOs, soa::Join<aod::BplusChargedEventWiseSubtractedJets, aod::BplusChargedEventWiseSubtractedJetConstituents, aod::BplusCEWSJetSSs, aod::BplusChargedEventWiseSubtractedJetsMatchedToBplusChargedJets>, aod::BplusCEWSJetCOs, aod::BplusCEWSJetOs, aod::BplusCEWSJetSSOs, aod::BplusCEWSJetMOs, aod::StoredHfBplusCollBase, aod::StoredHfBplusBases, aod::StoredHfBplusPars, aod::StoredHfBplusParEs, aod::StoredHfBplusParD0s, aod::StoredHfBplusSels, aod::StoredHfBplusMls, aod::StoredHfBplusMlD0s, aod::StoredHfBplusMcs, aod::StoredHfBplusMcCollBases, aod::StoredHfBplusMcRCollIds, aod::StoredHfBplusPBases>; // all the 3P tables have been made into Bplus but they might be made common

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputBplus>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructure-bplus-output"}));

  return WorkflowSpec{tasks};
}
