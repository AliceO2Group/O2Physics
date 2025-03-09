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

// substructure matching mc charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <vector>
#include "PWGJE/TableProducer/Matching/Substructure/jetSubstructureMatching.cxx"

using ChargedJetSubstructureMatchingMC = JetSubstructureMatching<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>,
                                                                 soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>,
                                                                 aod::ChargedMCDetectorLevelSPsMatchedToChargedMCParticleLevelSPs,
                                                                 aod::ChargedMCParticleLevelSPsMatchedToChargedMCDetectorLevelSPs,
                                                                 aod::ChargedMCDetectorLevelPRsMatchedToChargedMCParticleLevelPRs,
                                                                 aod::ChargedMCParticleLevelPRsMatchedToChargedMCDetectorLevelPRs,
                                                                 aod::ChargedMCDetectorLevelSPs,
                                                                 aod::ChargedMCParticleLevelSPs,
                                                                 aod::ChargedMCDetectorLevelPRs,
                                                                 aod::ChargedMCParticleLevelPRs,
                                                                 aod::JCollisions,
                                                                 aod::JMcCollisions,
                                                                 aod::JetTracksMCD,
                                                                 aod::JetParticles,
                                                                 aod::JDummys>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<ChargedJetSubstructureMatchingMC>(cfgc, TaskName{"jet-substructure-matching-mc-ch"}));
  return WorkflowSpec{tasks};
}
