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

// jet finder B0 charged QA task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/Tasks/jetFinderHFQA.cxx"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/runDataProcessing.h>

#include <vector>

using JetFinderB0QATask = JetFinderHFQATask<aod::B0ChargedJets, aod::B0ChargedJetConstituents, aod::B0ChargedJetsMatchedToB0ChargedEventWiseSubtractedJets, aod::CandidatesB0Data, aod::B0ChargedMCDetectorLevelJets, aod::B0ChargedMCDetectorLevelJetConstituents, aod::B0ChargedMCDetectorLevelJetsMatchedToB0ChargedMCParticleLevelJets, aod::B0ChargedMCDetectorLevelJetEventWeights, aod::CandidatesB0MCD, aod::B0ChargedMCParticleLevelJets, aod::B0ChargedMCParticleLevelJetConstituents, aod::B0ChargedMCParticleLevelJetsMatchedToB0ChargedMCDetectorLevelJets, aod::B0ChargedMCParticleLevelJetEventWeights, aod::B0ChargedEventWiseSubtractedJets, aod::B0ChargedEventWiseSubtractedJetConstituents, aod::B0ChargedEventWiseSubtractedJetsMatchedToB0ChargedJets, aod::CandidatesB0MCP, aod::JTrackB0Subs, aod::BkgB0Rhos>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetFinderB0QATask>(cfgc,
                                                          SetDefaultProcesses{},
                                                          TaskName{"jet-finder-charged-b0-qa"}));
  return WorkflowSpec{tasks};
}
