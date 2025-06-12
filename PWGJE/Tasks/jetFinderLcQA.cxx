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

// jet finder Lc charged QA task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/Tasks/jetFinderHFQA.cxx"

using JetFinderLcQATask = JetFinderHFQATask<aod::LcChargedJets, aod::LcChargedJetConstituents, aod::LcChargedJetsMatchedToLcChargedEventWiseSubtractedJets, aod::CandidatesLcData, aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents, aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCParticleLevelJets, aod::LcChargedMCDetectorLevelJetEventWeights, aod::CandidatesLcMCD, aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents, aod::LcChargedMCParticleLevelJetsMatchedToLcChargedMCDetectorLevelJets, aod::LcChargedMCParticleLevelJetEventWeights, aod::LcChargedEventWiseSubtractedJets, aod::LcChargedEventWiseSubtractedJetConstituents, aod::LcChargedEventWiseSubtractedJetsMatchedToLcChargedJets, aod::CandidatesLcMCP, aod::JTrackLcSubs, aod::BkgLcRhos>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetFinderLcQATask>(cfgc,
                                                          SetDefaultProcesses{},
                                                          TaskName{"jet-finder-charged-lc-qa"}));
  return WorkflowSpec{tasks};
}
