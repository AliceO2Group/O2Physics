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

// jet substructure output B0 charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGHF/DataModel/DerivedTables.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/DataModel/JetSubtraction.h"
#include "PWGJE/Tasks/jetSubstructureHFOutput.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/runDataProcessing.h>

#include <vector>

using JetSubstructureOutputB0 = JetSubstructureHFOutputTask<aod::CollisionsB0, soa::Join<aod::McCollisionsB0, aod::HfB0McRCollIds>, aod::McCollisionsB0, aod::CandidatesB0Data, aod::CandidatesB0MCD, aod::CandidatesB0MCP, aod::BkgB0Rhos, aod::BkgB0McRhos, aod::JTrackB0Subs, soa::Join<aod::B0ChargedJets, aod::B0ChargedJetConstituents, aod::B0CJetSSs>, soa::Join<aod::B0ChargedJets, aod::B0ChargedJetConstituents, aod::B0CJetSSs, aod::B0ChargedJetsMatchedToB0ChargedEventWiseSubtractedJets>, soa::Join<aod::B0ChargedSPs, aod::B0ChargedSPsMatchedToB0ChargedEventWiseSubtractedSPs>, soa::Join<aod::B0ChargedPRs, aod::B0ChargedPRsMatchedToB0ChargedEventWiseSubtractedPRs>, aod::B0CJetCOs, aod::B0CJetOs, aod::B0CJetSSOs, aod::B0CJetMOs, soa::Join<aod::B0ChargedMCDetectorLevelJets, aod::B0ChargedMCDetectorLevelJetConstituents, aod::B0CMCDJetSSs, aod::B0ChargedMCDetectorLevelJetsMatchedToB0ChargedMCParticleLevelJets>, soa::Join<aod::B0ChargedMCDetectorLevelSPs, aod::B0ChargedMCDetectorLevelSPsMatchedToB0ChargedMCParticleLevelSPs>, soa::Join<aod::B0ChargedMCDetectorLevelPRs, aod::B0ChargedMCDetectorLevelPRsMatchedToB0ChargedMCParticleLevelPRs>, aod::B0CMCDJetCOs, aod::B0CMCDJetOs, aod::B0CMCDJetSSOs, aod::B0CMCDJetMOs, soa::Join<aod::B0ChargedMCParticleLevelJets, aod::B0ChargedMCParticleLevelJetConstituents, aod::B0CMCPJetSSs>, soa::Join<aod::B0ChargedMCParticleLevelJets, aod::B0ChargedMCParticleLevelJetConstituents, aod::B0CMCPJetSSs, aod::B0ChargedMCParticleLevelJetsMatchedToB0ChargedMCDetectorLevelJets>, soa::Join<aod::B0ChargedMCParticleLevelSPs, aod::B0ChargedMCParticleLevelSPsMatchedToB0ChargedMCDetectorLevelSPs>, soa::Join<aod::B0ChargedMCParticleLevelPRs, aod::B0ChargedMCParticleLevelPRsMatchedToB0ChargedMCDetectorLevelPRs>, aod::B0CMCPJetCOs, aod::B0CMCPJetMCCOs, aod::B0CMCPJetOs, aod::B0CMCPJetSSOs, aod::B0CMCPJetMOs, soa::Join<aod::B0ChargedEventWiseSubtractedJets, aod::B0ChargedEventWiseSubtractedJetConstituents, aod::B0CEWSJetSSs, aod::B0ChargedEventWiseSubtractedJetsMatchedToB0ChargedJets>, soa::Join<aod::B0ChargedEventWiseSubtractedSPs, aod::B0ChargedEventWiseSubtractedSPsMatchedToB0ChargedSPs>, soa::Join<aod::B0ChargedEventWiseSubtractedPRs, aod::B0ChargedEventWiseSubtractedPRsMatchedToB0ChargedPRs>, aod::B0CEWSJetCOs, aod::B0CEWSJetOs, aod::B0CEWSJetSSOs, aod::B0CEWSJetMOs, aod::StoredHfB0CollBase, aod::StoredHfB0Bases, aod::StoredHfB0Pars, aod::StoredHfB0ParEs, aod::StoredHfB0ParDpluss, aod::StoredHfB0Sels, aod::StoredHfB0Mls, aod::StoredHfB0MlDpluss, aod::StoredHfB0Mcs, aod::StoredHfB0McCollBases, aod::StoredHfB0McRCollIds, aod::StoredHfB0PBases>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputB0>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructure-b0-output"}));

  return WorkflowSpec{tasks};
}
