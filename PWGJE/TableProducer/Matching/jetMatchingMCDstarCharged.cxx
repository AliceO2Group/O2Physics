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

// jet matching mc D* charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/TableProducer/Matching/jetMatchingMC.cxx"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/runDataProcessing.h>

#include <vector>

using DstarChargedJetMatchingMC = JetMatchingMc<soa::Join<aod::DstarChargedMCDetectorLevelJets, aod::DstarChargedMCDetectorLevelJetConstituents>,
                                                soa::Join<aod::DstarChargedMCParticleLevelJets, aod::DstarChargedMCParticleLevelJetConstituents>,
                                                aod::DstarChargedMCDetectorLevelJetsMatchedToDstarChargedMCParticleLevelJets,
                                                aod::DstarChargedMCParticleLevelJetsMatchedToDstarChargedMCDetectorLevelJets,
                                                aod::CandidatesDstarMCD,
                                                aod::CandidatesDstarMCP,
                                                aod::JDummys>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<DstarChargedJetMatchingMC>(cfgc, TaskName{"jet-matching-mc-dstar-ch"}));
  return WorkflowSpec{tasks};
}
