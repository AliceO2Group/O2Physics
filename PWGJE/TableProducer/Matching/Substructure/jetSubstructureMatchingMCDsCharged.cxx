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

// Ds substructure matching mc charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/TableProducer/Matching/Substructure/jetSubstructureMatching.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/runDataProcessing.h>

#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using DsChargedJetSubstructureMatchingMC = JetSubstructureMatching<soa::Join<aod::DsChargedMCDetectorLevelJets, aod::DsChargedMCDetectorLevelJetConstituents, aod::DsChargedMCDetectorLevelJetsMatchedToDsChargedMCParticleLevelJets>,
                                                                   soa::Join<aod::DsChargedMCParticleLevelJets, aod::DsChargedMCParticleLevelJetConstituents, aod::DsChargedMCParticleLevelJetsMatchedToDsChargedMCDetectorLevelJets>,
                                                                   aod::DsChargedMCDetectorLevelSPsMatchedToDsChargedMCParticleLevelSPs,
                                                                   aod::DsChargedMCParticleLevelSPsMatchedToDsChargedMCDetectorLevelSPs,
                                                                   aod::DsChargedMCDetectorLevelPRsMatchedToDsChargedMCParticleLevelPRs,
                                                                   aod::DsChargedMCParticleLevelPRsMatchedToDsChargedMCDetectorLevelPRs,
                                                                   aod::DsChargedMCDetectorLevelSPs,
                                                                   aod::DsChargedMCParticleLevelSPs,
                                                                   aod::DsChargedMCDetectorLevelPRs,
                                                                   aod::DsChargedMCParticleLevelPRs,
                                                                   aod::CandidatesDsMCD,
                                                                   aod::CandidatesDsMCP,
                                                                   aod::JetTracksMCD,
                                                                   aod::JetParticles,
                                                                   aod::JDummys>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<DsChargedJetSubstructureMatchingMC>(cfgc, TaskName{"jet-substructure-matching-mc-ds-ch"}));
  return WorkflowSpec{tasks};
}
