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

// jet matching mc subtracted B+ charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/TableProducer/Matching/jetMatchingMCSub.cxx"

using D0ChargedJetMatchingMCSub = JetMatchingMcSub<soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents>,
                                                   soa::Join<aod::D0ChargedMCDetectorLevelEventWiseSubtractedJets, aod::D0ChargedMCDetectorLevelEventWiseSubtractedJetConstituents>,
                                                   aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCDetectorLevelEventWiseSubtractedJets,
                                                   aod::D0ChargedMCDetectorLevelEventWiseSubtractedJetsMatchedToD0ChargedMCDetectorLevelJets,
                                                   aod::CandidatesD0MCD>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<D0ChargedJetMatchingMCSub>(cfgc, TaskName{"jet-matching-mc-sub-d0-ch"}));
  return WorkflowSpec{tasks};
}
