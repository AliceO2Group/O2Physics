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

// jet substructure output D* charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/Tasks/jetSubstructureHFOutput.cxx"

#include "PWGHF/DataModel/DerivedTables.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedDataHF.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/runDataProcessing.h>

#include <vector>

using JetSubstructureOutputDstar = JetSubstructureHFOutputTask<aod::CollisionsDstar, soa::Join<aod::McCollisionsDstar, aod::HfDstarMcRCollIds>, aod::McCollisionsDstar, aod::CandidatesDstarData, aod::CandidatesDstarMCD, aod::CandidatesDstarMCP, aod::BkgDstarRhos, aod::BkgDstarMcRhos, aod::JTrackDstarSubs, soa::Join<aod::DstarChargedJets, aod::DstarChargedJetConstituents, aod::DstarCJetSSs>, soa::Join<aod::DstarChargedJets, aod::DstarChargedJetConstituents, aod::DstarCJetSSs, aod::DstarChargedJetsMatchedToDstarChargedEventWiseSubtractedJets>, soa::Join<aod::DstarChargedSPs, aod::DstarChargedSPsMatchedToDstarChargedEventWiseSubtractedSPs>, soa::Join<aod::DstarChargedPRs, aod::DstarChargedPRsMatchedToDstarChargedEventWiseSubtractedPRs>, aod::DstarCJetCOs, aod::DstarCJetOs, aod::DstarCJetSSOs, aod::DstarCJetMOs, soa::Join<aod::DstarChargedMCDetectorLevelJets, aod::DstarChargedMCDetectorLevelJetConstituents, aod::DstarCMCDJetSSs, aod::DstarChargedMCDetectorLevelJetsMatchedToDstarChargedMCParticleLevelJets>, soa::Join<aod::DstarChargedMCDetectorLevelSPs, aod::DstarChargedMCDetectorLevelSPsMatchedToDstarChargedMCParticleLevelSPs>, soa::Join<aod::DstarChargedMCDetectorLevelPRs, aod::DstarChargedMCDetectorLevelPRsMatchedToDstarChargedMCParticleLevelPRs>, aod::DstarCMCDJetCOs, aod::DstarCMCDJetOs, aod::DstarCMCDJetSSOs, aod::DstarCMCDJetMOs, soa::Join<aod::DstarChargedMCParticleLevelJets, aod::DstarChargedMCParticleLevelJetConstituents, aod::DstarCMCPJetSSs>, soa::Join<aod::DstarChargedMCParticleLevelJets, aod::DstarChargedMCParticleLevelJetConstituents, aod::DstarCMCPJetSSs, aod::DstarChargedMCParticleLevelJetsMatchedToDstarChargedMCDetectorLevelJets>, soa::Join<aod::DstarChargedMCParticleLevelSPs, aod::DstarChargedMCParticleLevelSPsMatchedToDstarChargedMCDetectorLevelSPs>, soa::Join<aod::DstarChargedMCParticleLevelPRs, aod::DstarChargedMCParticleLevelPRsMatchedToDstarChargedMCDetectorLevelPRs>, aod::DstarCMCPJetCOs, aod::DstarCMCPJetMCCOs, aod::DstarCMCPJetOs, aod::DstarCMCPJetSSOs, aod::DstarCMCPJetMOs, soa::Join<aod::DstarChargedEventWiseSubtractedJets, aod::DstarChargedEventWiseSubtractedJetConstituents, aod::DstarCEWSJetSSs, aod::DstarChargedEventWiseSubtractedJetsMatchedToDstarChargedJets>, soa::Join<aod::DstarChargedEventWiseSubtractedSPs, aod::DstarChargedEventWiseSubtractedSPsMatchedToDstarChargedSPs>, soa::Join<aod::DstarChargedEventWiseSubtractedPRs, aod::DstarChargedEventWiseSubtractedPRsMatchedToDstarChargedPRs>, aod::DstarCEWSJetCOs, aod::DstarCEWSJetOs, aod::DstarCEWSJetSSOs, aod::DstarCEWSJetMOs, aod::StoredHfDstarCollBase, aod::StoredHfDstarBases, aod::StoredHfDstarPars, aod::JDumDstarParEs, aod::HfDstarParD0s, aod::StoredHfDstarSels, aod::StoredHfDstarMls, aod::JDumDstarMlDaus, aod::StoredHfDstarMcs, aod::StoredHfDstarMcCollBases, aod::StoredHfDstarMcRCollIds, aod::StoredHfDstarPBases>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputDstar>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructure-dstar-output"}));

  return WorkflowSpec{tasks};
}
