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

// jet substructure B0 charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/Tasks/jetSubstructureHF.cxx"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/runDataProcessing.h>

#include <vector>

using JetSubstructureB0 = JetSubstructureHFTask<soa::Join<aod::B0ChargedJets, aod::B0ChargedJetConstituents>, soa::Join<aod::B0ChargedMCDetectorLevelJets, aod::B0ChargedMCDetectorLevelJetConstituents>, soa::Join<aod::B0ChargedMCParticleLevelJets, aod::B0ChargedMCParticleLevelJetConstituents>, soa::Join<aod::B0ChargedEventWiseSubtractedJets, aod::B0ChargedEventWiseSubtractedJetConstituents>, aod::CandidatesB0Data, aod::CandidatesB0MCP, aod::B0CJetSSs, aod::B0ChargedSPs, aod::B0ChargedPRs, aod::B0CMCDJetSSs, aod::B0ChargedMCDetectorLevelSPs, aod::B0ChargedMCDetectorLevelPRs, aod::B0CMCPJetSSs, aod::B0ChargedMCParticleLevelSPs, aod::B0ChargedMCParticleLevelPRs, aod::B0CEWSJetSSs, aod::B0ChargedEventWiseSubtractedSPs, aod::B0ChargedEventWiseSubtractedPRs, aod::JTrackB0Subs>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureB0>(cfgc,
                                                          SetDefaultProcesses{},
                                                          TaskName{"jet-substructure-b0"}));
  return WorkflowSpec{tasks};
}
