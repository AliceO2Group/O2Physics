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

// jet finder v0 mcp charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/JetFinders/jetfinderv0.cxx"

using JetFinderV0MCParticleLevelCharged = JetFinderV0Task<CandidatesV0Data, CandidatesV0MCD, CandidatesV0MCP, aod::V0ChargedMCParticleLevelJets, aod::V0ChargedMCParticleLevelJetConstituents>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetFinderV0MCParticleLevelCharged>(cfgc,
                                                                          SetDefaultProcesses{{{"processChargedJetsMCP", true}}},
                                                                          TaskName{"jet-finder-v0-mcp-charged"}));

  return WorkflowSpec{tasks};
}
