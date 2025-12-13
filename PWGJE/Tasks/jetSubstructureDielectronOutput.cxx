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

#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedDataDQ.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/DataModel/JetSubtraction.h"
#include "PWGJE/Tasks/jetSubstructureHFOutput.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/runDataProcessing.h>

#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using JetSubstructureOutputDielectron = JetSubstructureHFOutputTask<aod::CollisionsDielectron, soa::Join<aod::McCollisionsDielectron, aod::JDielectronMcRCollDummys>, aod::McCollisionsDielectron, aod::CandidatesDielectronData, aod::CandidatesDielectronMCD, aod::CandidatesDielectronMCP, aod::BkgDielectronRhos, aod::BkgDielectronMcRhos, aod::JTrackDielectronSubs, soa::Join<aod::DielectronChargedJets, aod::DielectronChargedJetConstituents, aod::DielectronCJetSSs>, soa::Join<aod::DielectronChargedJets, aod::DielectronChargedJetConstituents, aod::DielectronCJetSSs, aod::DielectronChargedJetsMatchedToDielectronChargedEventWiseSubtractedJets>, soa::Join<aod::DielectronChargedSPs, aod::DielectronChargedSPsMatchedToDielectronChargedEventWiseSubtractedSPs>, soa::Join<aod::DielectronChargedPRs, aod::DielectronChargedPRsMatchedToDielectronChargedEventWiseSubtractedPRs>, aod::DielectronCJetCOs, aod::DielectronCJetOs, aod::DielectronCJetSSOs, aod::DielectronCJetMOs, soa::Join<aod::DielectronChargedMCDetectorLevelJets, aod::DielectronChargedMCDetectorLevelJetConstituents, aod::DielectronCMCDJetSSs, aod::DielectronChargedMCDetectorLevelJetsMatchedToDielectronChargedMCParticleLevelJets>, soa::Join<aod::DielectronChargedMCDetectorLevelSPs, aod::DielectronChargedMCDetectorLevelSPsMatchedToDielectronChargedMCParticleLevelSPs>, soa::Join<aod::DielectronChargedMCDetectorLevelPRs, aod::DielectronChargedMCDetectorLevelPRsMatchedToDielectronChargedMCParticleLevelPRs>, aod::DielectronCMCDJetCOs, aod::DielectronCMCDJetOs, aod::DielectronCMCDJetSSOs, aod::DielectronCMCDJetMOs, soa::Join<aod::DielectronChargedMCParticleLevelJets, aod::DielectronChargedMCParticleLevelJetConstituents, aod::DielectronCMCPJetSSs>, soa::Join<aod::DielectronChargedMCParticleLevelJets, aod::DielectronChargedMCParticleLevelJetConstituents, aod::DielectronCMCPJetSSs, aod::DielectronChargedMCParticleLevelJetsMatchedToDielectronChargedMCDetectorLevelJets>, soa::Join<aod::DielectronChargedMCParticleLevelSPs, aod::DielectronChargedMCParticleLevelSPsMatchedToDielectronChargedMCDetectorLevelSPs>, soa::Join<aod::DielectronChargedMCParticleLevelPRs, aod::DielectronChargedMCParticleLevelPRsMatchedToDielectronChargedMCDetectorLevelPRs>, aod::DielectronCMCPJetCOs, aod::DielectronCMCPJetMCCOs, aod::DielectronCMCPJetOs, aod::DielectronCMCPJetSSOs, aod::DielectronCMCPJetMOs, soa::Join<aod::DielectronChargedEventWiseSubtractedJets, aod::DielectronChargedEventWiseSubtractedJetConstituents, aod::DielectronCEWSJetSSs, aod::DielectronChargedEventWiseSubtractedJetsMatchedToDielectronChargedJets>, soa::Join<aod::DielectronChargedEventWiseSubtractedSPs, aod::DielectronChargedEventWiseSubtractedSPsMatchedToDielectronChargedSPs>, soa::Join<aod::DielectronChargedEventWiseSubtractedPRs, aod::DielectronChargedEventWiseSubtractedPRsMatchedToDielectronChargedPRs>, aod::DielectronCEWSJetCOs, aod::DielectronCEWSJetOs, aod::DielectronCEWSJetSSOs, aod::DielectronCEWSJetMOs, aod::StoredReducedEvents, aod::StoredDielectrons, aod::StoredDielectronsAll, aod::JDielectron1Dummys, aod::JDielectron2Dummys, aod::JDielectron3Dummys, aod::JDielectron4Dummys, aod::JDielectron5Dummys, aod::JDielectron6Dummys, aod::StoredJDielectronMcCollisions, aod::JDielectronMcRCollDummys, aod::StoredJDielectronMcs>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputDielectron>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructure-dielectron-output"}));

  return WorkflowSpec{tasks};
}
