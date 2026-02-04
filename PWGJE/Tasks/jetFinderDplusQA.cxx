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

// jet finder D+ charged QA task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubtraction.h"
#include "PWGJE/Tasks/jetFinderHFQA.h"

#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/runDataProcessing.h>

#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using JetFinderDplusQATask = JetFinderHFQATask<aod::DplusChargedJets, aod::DplusChargedJetConstituents, aod::DplusChargedJetsMatchedToDplusChargedEventWiseSubtractedJets, aod::CandidatesDplusData, aod::DplusChargedMCDetectorLevelJets, aod::DplusChargedMCDetectorLevelJetConstituents, aod::DplusChargedMCDetectorLevelJetsMatchedToDplusChargedMCParticleLevelJets, aod::DplusChargedMCDetectorLevelJetEventWeights, aod::CandidatesDplusMCD, aod::DplusChargedMCParticleLevelJets, aod::DplusChargedMCParticleLevelJetConstituents, aod::DplusChargedMCParticleLevelJetsMatchedToDplusChargedMCDetectorLevelJets, aod::DplusChargedMCParticleLevelJetEventWeights, aod::DplusChargedEventWiseSubtractedJets, aod::DplusChargedEventWiseSubtractedJetConstituents, aod::DplusChargedEventWiseSubtractedJetsMatchedToDplusChargedJets, aod::CandidatesDplusMCP, aod::JTrackDplusSubs, aod::BkgDplusRhos>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetFinderDplusQATask>(cfgc,
                                                             SetDefaultProcesses{},
                                                             TaskName{"jet-finder-charged-dplus-qa"}));
  return WorkflowSpec{tasks};
}
