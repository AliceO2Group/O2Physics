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

// jet substructure output D+ charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGHF/DataModel/DerivedTables.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedDataHF.h"
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

using JetSubstructureOutputDplus = JetSubstructureHFOutputTask<aod::CollisionsDplus, soa::Join<aod::McCollisionsDplus, aod::HfDplusMcRCollIds>, aod::McCollisionsDplus, aod::CandidatesDplusData, aod::CandidatesDplusMCD, aod::CandidatesDplusMCP, aod::BkgDplusRhos, aod::BkgDplusMcRhos, aod::JTrackDplusSubs, soa::Join<aod::DplusChargedJets, aod::DplusChargedJetConstituents, aod::DplusCJetSSs>, soa::Join<aod::DplusChargedJets, aod::DplusChargedJetConstituents, aod::DplusCJetSSs, aod::DplusChargedJetsMatchedToDplusChargedEventWiseSubtractedJets>, aod::DplusCJetRs, soa::Join<aod::DplusChargedSPs, aod::DplusChargedSPsMatchedToDplusChargedEventWiseSubtractedSPs>, soa::Join<aod::DplusChargedPRs, aod::DplusChargedPRsMatchedToDplusChargedEventWiseSubtractedPRs>, aod::DplusCJetCOs, aod::DplusCJetOs, aod::DplusCJetSSOs, aod::DplusCJetMOs, aod::DplusCJetROs, soa::Join<aod::DplusChargedMCDetectorLevelJets, aod::DplusChargedMCDetectorLevelJetConstituents, aod::DplusCMCDJetSSs, aod::DplusChargedMCDetectorLevelJetsMatchedToDplusChargedMCParticleLevelJets>, aod::DplusCMCDJetRs, soa::Join<aod::DplusChargedMCDetectorLevelSPs, aod::DplusChargedMCDetectorLevelSPsMatchedToDplusChargedMCParticleLevelSPs>, soa::Join<aod::DplusChargedMCDetectorLevelPRs, aod::DplusChargedMCDetectorLevelPRsMatchedToDplusChargedMCParticleLevelPRs>, aod::DplusCMCDJetCOs, aod::DplusCMCDJetOs, aod::DplusCMCDJetSSOs, aod::DplusCMCDJetMOs, aod::DplusCMCDJetROs, aod::DplusCMCDJetRMOs, soa::Join<aod::DplusChargedMCParticleLevelJets, aod::DplusChargedMCParticleLevelJetConstituents, aod::DplusCMCPJetSSs>, soa::Join<aod::DplusChargedMCParticleLevelJets, aod::DplusChargedMCParticleLevelJetConstituents, aod::DplusCMCPJetSSs, aod::DplusChargedMCParticleLevelJetsMatchedToDplusChargedMCDetectorLevelJets>, aod::DplusCMCPJetRs, soa::Join<aod::DplusChargedMCParticleLevelSPs, aod::DplusChargedMCParticleLevelSPsMatchedToDplusChargedMCDetectorLevelSPs>, soa::Join<aod::DplusChargedMCParticleLevelPRs, aod::DplusChargedMCParticleLevelPRsMatchedToDplusChargedMCDetectorLevelPRs>, aod::DplusCMCPJetCOs, aod::DplusCMCPJetMCCOs, aod::DplusCMCPJetOs, aod::DplusCMCPJetSSOs, aod::DplusCMCPJetMOs, aod::DplusCMCPJetROs, aod::DplusCMCPJetRMOs, soa::Join<aod::DplusChargedEventWiseSubtractedJets, aod::DplusChargedEventWiseSubtractedJetConstituents, aod::DplusCEWSJetSSs, aod::DplusChargedEventWiseSubtractedJetsMatchedToDplusChargedJets>, soa::Join<aod::DplusChargedEventWiseSubtractedSPs, aod::DplusChargedEventWiseSubtractedSPsMatchedToDplusChargedSPs>, soa::Join<aod::DplusChargedEventWiseSubtractedPRs, aod::DplusChargedEventWiseSubtractedPRsMatchedToDplusChargedPRs>, aod::DplusCEWSJetCOs, aod::DplusCEWSJetOs, aod::DplusCEWSJetSSOs, aod::DplusCEWSJetMOs, aod::StoredHfDplusCollBase, aod::StoredHfDplusBases, aod::StoredHfDplusPars, aod::StoredHfDplusParEs, aod::JDumDplusParDaus, aod::StoredHfDplusSels, aod::StoredHfDplusMls, aod::JDumDplusMlDaus, aod::StoredHfDplusMcs, aod::StoredHfDplusMcCollBases, aod::StoredHfDplusMcRCollIds, aod::StoredHfDplusPBases>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputDplus>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructure-dplus-output"}));

  return WorkflowSpec{tasks};
}
