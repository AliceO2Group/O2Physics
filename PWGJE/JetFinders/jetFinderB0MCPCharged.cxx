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

// jet finder B0 mcp charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/JetFinders/jetFinderHF.cxx"

#include "PWGJE/DataModel/Jet.h"

#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/runDataProcessing.h>

#include <vector>

using JetFinderB0MCParticleLevelCharged = JetFinderHFTask<aod::CandidatesB0Data, aod::CandidatesB0MCD, aod::CandidatesB0MCP, aod::JetTracksSubB0, aod::JetParticlesSubB0, aod::B0ChargedMCParticleLevelJets, aod::B0ChargedMCParticleLevelJetConstituents, aod::B0ChargedMCParticleLevelEventWiseSubtractedJets, aod::B0ChargedMCParticleLevelEventWiseSubtractedJetConstituents>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetFinderB0MCParticleLevelCharged>(cfgc,
                                                                          SetDefaultProcesses{{{"processChargedJetsMCP", true}}},
                                                                          TaskName{"jet-finder-b0-mcp-charged"}));

  return WorkflowSpec{tasks};
}
