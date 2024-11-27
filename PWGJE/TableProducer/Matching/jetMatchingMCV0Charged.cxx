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

// jet matching V0 charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/TableProducer/Matching/jetMatchingMC.cxx"

using V0ChargedJetMatchingMC = JetMatchingMc<soa::Join<aod::V0ChargedMCDetectorLevelJets, aod::V0ChargedMCDetectorLevelJetConstituents>,
                                             soa::Join<aod::V0ChargedMCParticleLevelJets, aod::V0ChargedMCParticleLevelJetConstituents>,
                                             aod::V0ChargedMCDetectorLevelJetsMatchedToV0ChargedMCParticleLevelJets,
                                             aod::V0ChargedMCParticleLevelJetsMatchedToV0ChargedMCDetectorLevelJets,
                                             aod::CandidatesV0MCD,
                                             aod::CandidatesV0MCP,
                                             aod::JDummys>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<V0ChargedJetMatchingMC>(cfgc, TaskName{"jet-matching-mc-v0-ch"}));
  return WorkflowSpec{tasks};
}
