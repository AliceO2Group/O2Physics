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

// Dielectron substructure matching mc charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <vector>
#include "PWGJE/TableProducer/Matching/Substructure/jetSubstructureMatching.cxx"

using DielectronChargedJetSubstructureMatchingMC = JetSubstructureMatching<soa::Join<aod::DielectronChargedMCDetectorLevelJets, aod::DielectronChargedMCDetectorLevelJetConstituents, aod::DielectronChargedMCDetectorLevelJetsMatchedToDielectronChargedMCParticleLevelJets>,
                                                                           soa::Join<aod::DielectronChargedMCParticleLevelJets, aod::DielectronChargedMCParticleLevelJetConstituents, aod::DielectronChargedMCParticleLevelJetsMatchedToDielectronChargedMCDetectorLevelJets>,
                                                                           aod::DielectronChargedMCDetectorLevelSPsMatchedToDielectronChargedMCParticleLevelSPs,
                                                                           aod::DielectronChargedMCParticleLevelSPsMatchedToDielectronChargedMCDetectorLevelSPs,
                                                                           aod::DielectronChargedMCDetectorLevelPRsMatchedToDielectronChargedMCParticleLevelPRs,
                                                                           aod::DielectronChargedMCParticleLevelPRsMatchedToDielectronChargedMCDetectorLevelPRs,
                                                                           aod::DielectronChargedMCDetectorLevelSPs,
                                                                           aod::DielectronChargedMCParticleLevelSPs,
                                                                           aod::DielectronChargedMCDetectorLevelPRs,
                                                                           aod::DielectronChargedMCParticleLevelPRs,
                                                                           aod::CandidatesDielectronMCD,
                                                                           aod::CandidatesDielectronMCP,
                                                                           aod::JetTracksMCD,
                                                                           aod::JetParticles,
                                                                           aod::JDummys>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<DielectronChargedJetSubstructureMatchingMC>(cfgc, TaskName{"jet-substructure-matching-mc-dielectron-ch"}));
  return WorkflowSpec{tasks};
}
