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

// jet substructure output Ds charged task
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

using JetSubstructureOutputDs = JetSubstructureHFOutputTask<aod::CollisionsDs, soa::Join<aod::McCollisionsDs, aod::HfDsMcRCollIds>, aod::McCollisionsDs, aod::CandidatesDsData, aod::CandidatesDsMCD, aod::CandidatesDsMCP, aod::BkgDsRhos, aod::BkgDsMcRhos, aod::JTrackDsSubs, soa::Join<aod::DsChargedJets, aod::DsChargedJetConstituents, aod::DsCJetSSs>, soa::Join<aod::DsChargedJets, aod::DsChargedJetConstituents, aod::DsCJetSSs, aod::DsChargedJetsMatchedToDsChargedEventWiseSubtractedJets>, aod::DsCJetRs, soa::Join<aod::DsChargedSPs, aod::DsChargedSPsMatchedToDsChargedEventWiseSubtractedSPs>, soa::Join<aod::DsChargedPRs, aod::DsChargedPRsMatchedToDsChargedEventWiseSubtractedPRs>, aod::DsCJetCOs, aod::DsCJetOs, aod::DsCJetSSOs, aod::DsCJetMOs, aod::DsCJetROs, soa::Join<aod::DsChargedMCDetectorLevelJets, aod::DsChargedMCDetectorLevelJetConstituents, aod::DsCMCDJetSSs, aod::DsChargedMCDetectorLevelJetsMatchedToDsChargedMCParticleLevelJets>, aod::DsCMCDJetRs, soa::Join<aod::DsChargedMCDetectorLevelSPs, aod::DsChargedMCDetectorLevelSPsMatchedToDsChargedMCParticleLevelSPs>, soa::Join<aod::DsChargedMCDetectorLevelPRs, aod::DsChargedMCDetectorLevelPRsMatchedToDsChargedMCParticleLevelPRs>, aod::DsCMCDJetCOs, aod::DsCMCDJetOs, aod::DsCMCDJetSSOs, aod::DsCMCDJetMOs, aod::DsCMCDJetROs, aod::DsCMCDJetRMOs, soa::Join<aod::DsChargedMCParticleLevelJets, aod::DsChargedMCParticleLevelJetConstituents, aod::DsCMCPJetSSs>, soa::Join<aod::DsChargedMCParticleLevelJets, aod::DsChargedMCParticleLevelJetConstituents, aod::DsCMCPJetSSs, aod::DsChargedMCParticleLevelJetsMatchedToDsChargedMCDetectorLevelJets>, aod::DsCMCPJetRs, soa::Join<aod::DsChargedMCParticleLevelSPs, aod::DsChargedMCParticleLevelSPsMatchedToDsChargedMCDetectorLevelSPs>, soa::Join<aod::DsChargedMCParticleLevelPRs, aod::DsChargedMCParticleLevelPRsMatchedToDsChargedMCDetectorLevelPRs>, aod::DsCMCPJetCOs, aod::DsCMCPJetMCCOs, aod::DsCMCPJetOs, aod::DsCMCPJetSSOs, aod::DsCMCPJetMOs, aod::DsCMCPJetROs, aod::DsCMCPJetRMOs, soa::Join<aod::DsChargedEventWiseSubtractedJets, aod::DsChargedEventWiseSubtractedJetConstituents, aod::DsCEWSJetSSs, aod::DsChargedEventWiseSubtractedJetsMatchedToDsChargedJets>, soa::Join<aod::DsChargedEventWiseSubtractedSPs, aod::DsChargedEventWiseSubtractedSPsMatchedToDsChargedSPs>, soa::Join<aod::DsChargedEventWiseSubtractedPRs, aod::DsChargedEventWiseSubtractedPRsMatchedToDsChargedPRs>, aod::DsCEWSJetCOs, aod::DsCEWSJetOs, aod::DsCEWSJetSSOs, aod::DsCEWSJetMOs, aod::StoredHfDsCollBase, aod::StoredHfDsBases, aod::StoredHfDsPars, aod::StoredHfDsParEs, aod::JDumDsParDaus, aod::StoredHfDsSels, aod::StoredHfDsMls, aod::JDumDsMlDaus, aod::StoredHfDsMcs, aod::StoredHfDsMcCollBases, aod::StoredHfDsMcRCollIds, aod::StoredHfDsPBases>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputDs>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructure-ds-output"}));

  return WorkflowSpec{tasks};
}
