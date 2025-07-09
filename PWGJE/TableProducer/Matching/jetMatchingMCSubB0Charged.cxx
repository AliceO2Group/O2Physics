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

// jet matching mc subtracted B0 charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/TableProducer/Matching/jetMatchingMCSub.cxx"

#include "PWGJE/DataModel/Jet.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/runDataProcessing.h>

#include <vector>

using B0ChargedJetMatchingMCSub = JetMatchingMcSub<soa::Join<aod::B0ChargedMCDetectorLevelJets, aod::B0ChargedMCDetectorLevelJetConstituents>,
                                                   soa::Join<aod::B0ChargedMCDetectorLevelEventWiseSubtractedJets, aod::B0ChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                                   aod::B0ChargedMCDetectorLevelJetsMatchedToB0ChargedMCDetectorLevelEventWiseSubtractedJets,
                                                   aod::B0ChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToB0ChargedMCDetectorLevelJets,
                                                   aod::CandidatesB0MCD>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<B0ChargedJetMatchingMCSub>(cfgc, TaskName{"jet-matching-mc-sub-b0-ch"}));
  return WorkflowSpec{tasks};
}
