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

// jet substructure output XicToXiPiPi charged task
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

using JetSubstructureOutputXicToXiPiPi = JetSubstructureHFOutputTask<aod::CollisionsXicToXiPiPi, soa::Join<aod::McCollisionsXicToXiPiPi, aod::HfXicToXiPiPiMcRCollIds>, aod::McCollisionsXicToXiPiPi, aod::CandidatesXicToXiPiPiData, aod::CandidatesXicToXiPiPiMCD, aod::CandidatesXicToXiPiPiMCP, aod::BkgXicToXiPiPiRhos, aod::BkgXicToXiPiPiMcRhos, aod::JTrackXicToXiPiPiSubs, soa::Join<aod::XicToXiPiPiChargedJets, aod::XicToXiPiPiChargedJetConstituents, aod::XicToXiPiPiCJetSSs>, soa::Join<aod::XicToXiPiPiChargedJets, aod::XicToXiPiPiChargedJetConstituents, aod::XicToXiPiPiCJetSSs, aod::XicToXiPiPiChargedJetsMatchedToXicToXiPiPiChargedEventWiseSubtractedJets>, soa::Join<aod::XicToXiPiPiChargedSPs, aod::XicToXiPiPiChargedSPsMatchedToXicToXiPiPiChargedEventWiseSubtractedSPs>, soa::Join<aod::XicToXiPiPiChargedPRs, aod::XicToXiPiPiChargedPRsMatchedToXicToXiPiPiChargedEventWiseSubtractedPRs>, aod::XicToXiPiPiCJetCOs, aod::XicToXiPiPiCJetOs, aod::XicToXiPiPiCJetSSOs, aod::XicToXiPiPiCJetMOs, soa::Join<aod::XicToXiPiPiChargedMCDetectorLevelJets, aod::XicToXiPiPiChargedMCDetectorLevelJetConstituents, aod::XicToXiPiPiCMCDJetSSs, aod::XicToXiPiPiChargedMCDetectorLevelJetsMatchedToXicToXiPiPiChargedMCParticleLevelJets>, soa::Join<aod::XicToXiPiPiChargedMCDetectorLevelSPs, aod::XicToXiPiPiChargedMCDetectorLevelSPsMatchedToXicToXiPiPiChargedMCParticleLevelSPs>, soa::Join<aod::XicToXiPiPiChargedMCDetectorLevelPRs, aod::XicToXiPiPiChargedMCDetectorLevelPRsMatchedToXicToXiPiPiChargedMCParticleLevelPRs>, aod::XicToXiPiPiCMCDJetCOs, aod::XicToXiPiPiCMCDJetOs, aod::XicToXiPiPiCMCDJetSSOs, aod::XicToXiPiPiCMCDJetMOs, soa::Join<aod::XicToXiPiPiChargedMCParticleLevelJets, aod::XicToXiPiPiChargedMCParticleLevelJetConstituents, aod::XicToXiPiPiCMCPJetSSs>, soa::Join<aod::XicToXiPiPiChargedMCParticleLevelJets, aod::XicToXiPiPiChargedMCParticleLevelJetConstituents, aod::XicToXiPiPiCMCPJetSSs, aod::XicToXiPiPiChargedMCParticleLevelJetsMatchedToXicToXiPiPiChargedMCDetectorLevelJets>, soa::Join<aod::XicToXiPiPiChargedMCParticleLevelSPs, aod::XicToXiPiPiChargedMCParticleLevelSPsMatchedToXicToXiPiPiChargedMCDetectorLevelSPs>, soa::Join<aod::XicToXiPiPiChargedMCParticleLevelPRs, aod::XicToXiPiPiChargedMCParticleLevelPRsMatchedToXicToXiPiPiChargedMCDetectorLevelPRs>, aod::XicToXiPiPiCMCPJetCOs, aod::XicToXiPiPiCMCPJetMCCOs, aod::XicToXiPiPiCMCPJetOs, aod::XicToXiPiPiCMCPJetSSOs, aod::XicToXiPiPiCMCPJetMOs, soa::Join<aod::XicToXiPiPiChargedEventWiseSubtractedJets, aod::XicToXiPiPiChargedEventWiseSubtractedJetConstituents, aod::XicToXiPiPiCEWSJetSSs, aod::XicToXiPiPiChargedEventWiseSubtractedJetsMatchedToXicToXiPiPiChargedJets>, soa::Join<aod::XicToXiPiPiChargedEventWiseSubtractedSPs, aod::XicToXiPiPiChargedEventWiseSubtractedSPsMatchedToXicToXiPiPiChargedSPs>, soa::Join<aod::XicToXiPiPiChargedEventWiseSubtractedPRs, aod::XicToXiPiPiChargedEventWiseSubtractedPRsMatchedToXicToXiPiPiChargedPRs>, aod::XicToXiPiPiCEWSJetCOs, aod::XicToXiPiPiCEWSJetOs, aod::XicToXiPiPiCEWSJetSSOs, aod::XicToXiPiPiCEWSJetMOs, aod::StoredHfXicToXiPiPiCollBase, aod::StoredHfXicToXiPiPiBases, aod::StoredHfXicToXiPiPiPars, aod::StoredHfXicToXiPiPiParEs, aod::JDumXicToXiPiPiParDaus, aod::StoredHfXicToXiPiPiSels, aod::StoredHfXicToXiPiPiMls, aod::JDumXicToXiPiPiMlDaus, aod::StoredHfXicToXiPiPiMcs, aod::StoredHfXicToXiPiPiMcCollBases, aod::StoredHfXicToXiPiPiMcRCollIds, aod::StoredHfXicToXiPiPiPBases>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputXicToXiPiPi>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructure-xictoxipipi-output"}));

  return WorkflowSpec{tasks};
}
