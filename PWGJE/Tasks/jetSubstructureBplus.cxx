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

// jet substructure B+ charged task
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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using JetSubstructureBplus = JetSubstructureHFTask<soa::Join<aod::BplusChargedJets, aod::BplusChargedJetConstituents>, soa::Join<aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents>, soa::Join<aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents>, soa::Join<aod::BplusChargedEventWiseSubtractedJets, aod::BplusChargedEventWiseSubtractedJetConstituents>, aod::CandidatesBplusData, aod::CandidatesBplusMCP, aod::BplusCJetSSs, aod::BplusCJetRs, aod::BplusChargedSPs, aod::BplusChargedPRs, aod::BplusCMCDJetSSs, aod::BplusCMCDJetRs, aod::BplusChargedMCDetectorLevelSPs, aod::BplusChargedMCDetectorLevelPRs, aod::BplusCMCPJetSSs, aod::BplusCMCPJetRs, aod::BplusChargedMCParticleLevelSPs, aod::BplusChargedMCParticleLevelPRs, aod::BplusCEWSJetSSs, aod::BplusChargedEventWiseSubtractedSPs, aod::BplusChargedEventWiseSubtractedPRs, aod::JTrackBplusSubs>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureBplus>(cfgc,
                                                             SetDefaultProcesses{},
                                                             TaskName{"jet-substructure-bplus"}));
  return WorkflowSpec{tasks};
}
