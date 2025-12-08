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

// jet substructure Dielectron charged task
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

using JetSubstructureDielectron = JetSubstructureHFTask<soa::Join<aod::DielectronChargedJets, aod::DielectronChargedJetConstituents>, soa::Join<aod::DielectronChargedMCDetectorLevelJets, aod::DielectronChargedMCDetectorLevelJetConstituents>, soa::Join<aod::DielectronChargedMCParticleLevelJets, aod::DielectronChargedMCParticleLevelJetConstituents>, soa::Join<aod::DielectronChargedEventWiseSubtractedJets, aod::DielectronChargedEventWiseSubtractedJetConstituents>, aod::CandidatesDielectronData, aod::CandidatesDielectronMCP, aod::DielectronCJetSSs, aod::DielectronChargedSPs, aod::DielectronChargedPRs, aod::DielectronCMCDJetSSs, aod::DielectronChargedMCDetectorLevelSPs, aod::DielectronChargedMCDetectorLevelPRs, aod::DielectronCMCPJetSSs, aod::DielectronChargedMCParticleLevelSPs, aod::DielectronChargedMCParticleLevelPRs, aod::DielectronCEWSJetSSs, aod::DielectronChargedEventWiseSubtractedSPs, aod::DielectronChargedEventWiseSubtractedPRs, aod::JTrackDielectronSubs>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureDielectron>(cfgc,
                                                                  SetDefaultProcesses{},
                                                                  TaskName{"jet-substructure-dielectron"}));
  return WorkflowSpec{tasks};
}
