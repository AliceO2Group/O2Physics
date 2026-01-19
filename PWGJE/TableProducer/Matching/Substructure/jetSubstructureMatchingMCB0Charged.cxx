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

// B0 substructure matching mc charged task
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

using B0ChargedJetSubstructureMatchingMC = JetSubstructureMatching<soa::Join<aod::B0ChargedMCDetectorLevelJets, aod::B0ChargedMCDetectorLevelJetConstituents, aod::B0ChargedMCDetectorLevelJetsMatchedToB0ChargedMCParticleLevelJets>,
                                                                   soa::Join<aod::B0ChargedMCParticleLevelJets, aod::B0ChargedMCParticleLevelJetConstituents, aod::B0ChargedMCParticleLevelJetsMatchedToB0ChargedMCDetectorLevelJets>,
                                                                   aod::B0ChargedMCDetectorLevelSPsMatchedToB0ChargedMCParticleLevelSPs,
                                                                   aod::B0ChargedMCParticleLevelSPsMatchedToB0ChargedMCDetectorLevelSPs,
                                                                   aod::B0ChargedMCDetectorLevelPRsMatchedToB0ChargedMCParticleLevelPRs,
                                                                   aod::B0ChargedMCParticleLevelPRsMatchedToB0ChargedMCDetectorLevelPRs,
                                                                   aod::B0ChargedMCDetectorLevelSPs,
                                                                   aod::B0ChargedMCParticleLevelSPs,
                                                                   aod::B0ChargedMCDetectorLevelPRs,
                                                                   aod::B0ChargedMCParticleLevelPRs,
                                                                   aod::CandidatesB0MCD,
                                                                   aod::CandidatesB0MCP,
                                                                   aod::JetTracksMCD,
                                                                   aod::JetParticles,
                                                                   aod::JDummys>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<B0ChargedJetSubstructureMatchingMC>(cfgc, TaskName{"jet-substructure-matching-mc-b0-ch"}));
  return WorkflowSpec{tasks};
}
