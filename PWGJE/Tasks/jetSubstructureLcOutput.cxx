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

// jet substructure output Lc charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/Tasks/jetSubstructureHFOutput.cxx"

using JetSubstructureOutputLc = JetSubstructureHFOutputTask<aod::CollisionsLc, soa::Join<aod::McCollisionsLc, aod::HfLcMcRCollIds>, aod::CandidatesLcData, aod::CandidatesLcMCD, aod::CandidatesLcMCP, aod::JTrackLcSubs, soa::Join<aod::LcChargedJets, aod::LcChargedJetConstituents, aod::LcCJetSSs>, soa::Join<aod::LcChargedJets, aod::LcChargedJetConstituents, aod::LcCJetSSs, aod::LcChargedJetsMatchedToLcChargedEventWiseSubtractedJets>, aod::LcCJetCOs, aod::LcCJetOs, aod::LcCJetSSOs, aod::LcCJetMOs, soa::Join<aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents, aod::LcCMCDJetSSs, aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCParticleLevelJets>, aod::LcCMCDJetCOs, aod::LcCMCDJetOs, aod::LcCMCDJetSSOs, aod::LcCMCDJetMOs, soa::Join<aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents, aod::LcCMCPJetSSs, aod::LcChargedMCParticleLevelJetsMatchedToLcChargedMCDetectorLevelJets>, aod::LcCMCPJetCOs, aod::LcCMCPJetOs, aod::LcCMCPJetSSOs, aod::LcCMCPJetMOs, soa::Join<aod::LcChargedEventWiseSubtractedJets, aod::LcChargedEventWiseSubtractedJetConstituents, aod::LcCEWSJetSSs, aod::LcChargedEventWiseSubtractedJetsMatchedToLcChargedJets>, aod::LcCEWSJetCOs, aod::LcCEWSJetOs, aod::LcCEWSJetSSOs, aod::LcCEWSJetMOs, aod::StoredHfLcCollBase, aod::StoredHfLcBases, aod::StoredHfLcPars, aod::StoredHfLcParEs, aod::JDumLcParDaus, aod::StoredHfLcSels, aod::StoredHfLcMls, aod::JDumLcMlDaus, aod::StoredHfLcMcs, aod::StoredHfLcMcCollBases, aod::StoredHfLcMcRCollIds, aod::StoredHfLcPBases>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputLc>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructure-lc-output"}));

  return WorkflowSpec{tasks};
}
