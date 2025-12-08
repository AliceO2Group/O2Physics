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

// jet substructure D* charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/DataModel/JetSubtraction.h"
#include "PWGJE/Tasks/jetSubstructureHF.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/runDataProcessing.h>

#include <vector>

using JetSubstructureDstar = JetSubstructureHFTask<soa::Join<aod::DstarChargedJets, aod::DstarChargedJetConstituents>, soa::Join<aod::DstarChargedMCDetectorLevelJets, aod::DstarChargedMCDetectorLevelJetConstituents>, soa::Join<aod::DstarChargedMCParticleLevelJets, aod::DstarChargedMCParticleLevelJetConstituents>, soa::Join<aod::DstarChargedEventWiseSubtractedJets, aod::DstarChargedEventWiseSubtractedJetConstituents>, aod::CandidatesDstarData, aod::CandidatesDstarMCP, aod::DstarCJetSSs, aod::DstarChargedSPs, aod::DstarChargedPRs, aod::DstarCMCDJetSSs, aod::DstarChargedMCDetectorLevelSPs, aod::DstarChargedMCDetectorLevelPRs, aod::DstarCMCPJetSSs, aod::DstarChargedMCParticleLevelSPs, aod::DstarChargedMCParticleLevelPRs, aod::DstarCEWSJetSSs, aod::DstarChargedEventWiseSubtractedSPs, aod::DstarChargedEventWiseSubtractedPRs, aod::JTrackDstarSubs>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureDstar>(cfgc,
                                                             SetDefaultProcesses{},
                                                             TaskName{"jet-substructure-dstar"}));
  return WorkflowSpec{tasks};
}
