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

// jet substructure output Dielectron charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/Tasks/jetSubstructureHFOutput.cxx"

using JetSubstructureOutputDielectron = JetSubstructureHFOutputTask<aod::CollisionsDielectron, aod::McCollisionsDielectron, aod::CandidatesDielectronData, aod::CandidatesDielectronMCD, aod::CandidatesDielectronMCP, aod::JTrackDielectronSubs, soa::Join<aod::DielectronChargedJets, aod::DielectronChargedJetConstituents, aod::DielectronCJetSSs>, soa::Join<aod::DielectronChargedJets, aod::DielectronChargedJetConstituents, aod::DielectronCJetSSs, aod::DielectronChargedJetsMatchedToDielectronChargedEventWiseSubtractedJets>, aod::DielectronCJetCOs, aod::DielectronCJetOs, aod::DielectronCJetSSOs, aod::DielectronCJetMOs, soa::Join<aod::DielectronChargedMCDetectorLevelJets, aod::DielectronChargedMCDetectorLevelJetConstituents, aod::DielectronCMCDJetSSs, aod::DielectronChargedMCDetectorLevelJetsMatchedToDielectronChargedMCParticleLevelJets>, aod::DielectronCMCDJetCOs, aod::DielectronCMCDJetOs, aod::DielectronCMCDJetSSOs, aod::DielectronCMCDJetMOs, soa::Join<aod::DielectronChargedMCParticleLevelJets, aod::DielectronChargedMCParticleLevelJetConstituents, aod::DielectronCMCPJetSSs, aod::DielectronChargedMCParticleLevelJetsMatchedToDielectronChargedMCDetectorLevelJets>, aod::DielectronCMCPJetCOs, aod::DielectronCMCPJetOs, aod::DielectronCMCPJetSSOs, aod::DielectronCMCPJetMOs, soa::Join<aod::DielectronChargedEventWiseSubtractedJets, aod::DielectronChargedEventWiseSubtractedJetConstituents, aod::DielectronCEWSJetSSs, aod::DielectronChargedEventWiseSubtractedJetsMatchedToDielectronChargedJets>, aod::DielectronCEWSJetCOs, aod::DielectronCEWSJetOs, aod::DielectronCEWSJetSSOs, aod::DielectronCEWSJetMOs, aod::StoredReducedEvents, aod::StoredDielectrons, aod::JDielectron1Dummys, aod::JDielectron2Dummys, aod::JDielectron3Dummys, aod::JDielectron4Dummys, aod::JDielectron5Dummys, aod::JDielectron6Dummys, aod::JDielectron7Dummys, aod::StoredJDielectronMcCollisions, aod::JDielectron8Dummys, aod::StoredJDielectronMcs>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputDielectron>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructure-dielectron-output"}));

  return WorkflowSpec{tasks};
}
