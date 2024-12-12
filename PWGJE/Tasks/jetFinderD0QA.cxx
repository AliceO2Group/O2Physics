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

// jet finder D0 charged QA task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/Tasks/jetFinderHFQA.cxx"

using JetFinderD0QATask = JetFinderHFQATask<aod::D0ChargedJets, aod::D0ChargedJetConstituents, aod::D0ChargedJetsMatchedToD0ChargedEventWiseSubtractedJets, aod::CandidatesD0Data, aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets, aod::D0ChargedMCDetectorLevelJetEventWeights, aod::CandidatesD0MCD, aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets, aod::D0ChargedMCParticleLevelJetEventWeights, aod::D0ChargedEventWiseSubtractedJets, aod::D0ChargedEventWiseSubtractedJetConstituents, aod::D0ChargedEventWiseSubtractedJetsMatchedToD0ChargedJets, aod::CandidatesD0MCP, aod::JTrackD0Subs, aod::BkgD0Rhos>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetFinderD0QATask>(cfgc,
                                                          SetDefaultProcesses{},
                                                          TaskName{"jet-finder-charged-d0-qa"}));
  return WorkflowSpec{tasks};
}
