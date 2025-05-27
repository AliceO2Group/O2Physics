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

// substructure matching event-wise subtracted Lc charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <vector>
#include "PWGJE/TableProducer/Matching/Substructure/jetSubstructureMatchingSub.cxx"

using LcChargedJetSubstructureMatchingSub = JetSubstructureMatchingSub<soa::Join<aod::LcChargedJets, aod::LcChargedJetConstituents, aod::LcChargedJetsMatchedToLcChargedEventWiseSubtractedJets>,
                                                                       soa::Join<aod::LcChargedEventWiseSubtractedJets, aod::LcChargedEventWiseSubtractedJetConstituents, aod::LcChargedEventWiseSubtractedJetsMatchedToLcChargedJets>,
                                                                       aod::LcChargedSPsMatchedToLcChargedEventWiseSubtractedSPs,
                                                                       aod::LcChargedEventWiseSubtractedSPsMatchedToLcChargedSPs,
                                                                       aod::LcChargedPRsMatchedToLcChargedEventWiseSubtractedPRs,
                                                                       aod::LcChargedEventWiseSubtractedPRsMatchedToLcChargedPRs,
                                                                       aod::LcChargedSPs,
                                                                       aod::LcChargedEventWiseSubtractedSPs,
                                                                       aod::LcChargedPRs,
                                                                       aod::LcChargedEventWiseSubtractedPRs,
                                                                       aod::CandidatesLcData,
                                                                       aod::JetTracks,
                                                                       aod::JetTracksSubLc,
                                                                       aod::JDummys>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<LcChargedJetSubstructureMatchingSub>(cfgc, TaskName{"jet-substructure-matching-sub-lc-ch"}));
  return WorkflowSpec{tasks};
}
