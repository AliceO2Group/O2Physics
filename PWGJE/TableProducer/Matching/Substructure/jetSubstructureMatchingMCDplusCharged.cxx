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

// Dplus substructure matching mc charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <vector>
#include "PWGJE/TableProducer/Matching/Substructure/jetSubstructureMatching.cxx"

using DplusChargedJetSubstructureMatchingMC = JetSubstructureMatching<soa::Join<aod::DplusChargedMCDetectorLevelJets, aod::DplusChargedMCDetectorLevelJetConstituents, aod::DplusChargedMCDetectorLevelJetsMatchedToDplusChargedMCParticleLevelJets>,
                                                                      soa::Join<aod::DplusChargedMCParticleLevelJets, aod::DplusChargedMCParticleLevelJetConstituents, aod::DplusChargedMCParticleLevelJetsMatchedToDplusChargedMCDetectorLevelJets>,
                                                                      aod::DplusChargedMCDetectorLevelSPsMatchedToDplusChargedMCParticleLevelSPs,
                                                                      aod::DplusChargedMCParticleLevelSPsMatchedToDplusChargedMCDetectorLevelSPs,
                                                                      aod::DplusChargedMCDetectorLevelPRsMatchedToDplusChargedMCParticleLevelPRs,
                                                                      aod::DplusChargedMCParticleLevelPRsMatchedToDplusChargedMCDetectorLevelPRs,
                                                                      aod::DplusChargedMCDetectorLevelSPs,
                                                                      aod::DplusChargedMCParticleLevelSPs,
                                                                      aod::DplusChargedMCDetectorLevelPRs,
                                                                      aod::DplusChargedMCParticleLevelPRs,
                                                                      aod::CandidatesDplusMCD,
                                                                      aod::CandidatesDplusMCP,
                                                                      aod::JetTracksMCD,
                                                                      aod::JetParticles,
                                                                      aod::JDummys>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<DplusChargedJetSubstructureMatchingMC>(cfgc, TaskName{"jet-substructure-matching-mc-dplus-ch"}));
  return WorkflowSpec{tasks};
}
