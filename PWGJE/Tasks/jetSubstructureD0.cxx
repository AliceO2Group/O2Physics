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

// jet substructure D0 charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <vector>
#include "PWGJE/Tasks/jetSubstructureHF.cxx"

using JetSubstructureD0 = JetSubstructureHFTask<soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents>, soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents>, soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents>, soa::Join<aod::D0ChargedEventWiseSubtractedJets, aod::D0ChargedEventWiseSubtractedJetConstituents>, aod::CandidatesD0Data, aod::CandidatesD0MCP, aod::D0CJetSSs, aod::D0ChargedSPs, aod::D0ChargedPRs, aod::D0CMCDJetSSs, aod::D0ChargedMCDetectorLevelSPs, aod::D0ChargedMCDetectorLevelPRs, aod::D0CMCPJetSSs, aod::D0ChargedMCParticleLevelSPs, aod::D0ChargedMCParticleLevelPRs, aod::D0CEWSJetSSs, aod::D0ChargedEventWiseSubtractedSPs, aod::D0ChargedEventWiseSubtractedPRs, aod::JTrackD0Subs>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureD0>(cfgc,
                                                          SetDefaultProcesses{},
                                                          TaskName{"jet-substructure-d0"}));
  return WorkflowSpec{tasks};
}
