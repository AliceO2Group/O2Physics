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

// D0 substructure matching mc charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <vector>
#include "PWGJE/TableProducer/Matching/Substructure/jetSubstructureMatching.cxx"

using D0ChargedJetSubstructureMatchingMC = JetSubstructureMatching<soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets>,
                                                                   soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets>,
                                                                   aod::D0ChargedMCDetectorLevelSPsMatchedToD0ChargedMCParticleLevelSPs,
                                                                   aod::D0ChargedMCParticleLevelSPsMatchedToD0ChargedMCDetectorLevelSPs,
                                                                   aod::D0ChargedMCDetectorLevelPRsMatchedToD0ChargedMCParticleLevelPRs,
                                                                   aod::D0ChargedMCParticleLevelPRsMatchedToD0ChargedMCDetectorLevelPRs,
                                                                   aod::D0ChargedMCDetectorLevelSPs,
                                                                   aod::D0ChargedMCParticleLevelSPs,
                                                                   aod::D0ChargedMCDetectorLevelPRs,
                                                                   aod::D0ChargedMCParticleLevelPRs,
                                                                   aod::CandidatesD0MCD,
                                                                   aod::CandidatesD0MCP,
                                                                   aod::JetTracksMCD,
                                                                   aod::JetParticles,
                                                                   aod::JDummys>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<D0ChargedJetSubstructureMatchingMC>(cfgc, TaskName{"jet-substructure-matching-mc-d0-ch"}));
  return WorkflowSpec{tasks};
}
