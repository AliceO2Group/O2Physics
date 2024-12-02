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

// jet finder B+ charged QA task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/Tasks/jetFinderHFQA.cxx"

using JetFinderBplusQATask = JetFinderHFQATask<aod::BplusChargedJets, aod::BplusChargedJetConstituents, aod::BplusChargedJetsMatchedToBplusChargedEventWiseSubtractedJets, aod::CandidatesBplusData, aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents, aod::BplusChargedMCDetectorLevelJetsMatchedToBplusChargedMCParticleLevelJets, aod::BplusChargedMCDetectorLevelJetEventWeights, aod::CandidatesBplusMCD, aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents, aod::BplusChargedMCParticleLevelJetsMatchedToBplusChargedMCDetectorLevelJets, aod::BplusChargedMCParticleLevelJetEventWeights, aod::BplusChargedEventWiseSubtractedJets, aod::BplusChargedEventWiseSubtractedJetConstituents, aod::BplusChargedEventWiseSubtractedJetsMatchedToBplusChargedJets, aod::CandidatesBplusMCP, aod::JTrackBplusSubs, aod::BkgBplusRhos>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetFinderBplusQATask>(cfgc,
                                                             SetDefaultProcesses{},
                                                             TaskName{"jet-finder-charged-bplus-qa"}));
  return WorkflowSpec{tasks};
}
