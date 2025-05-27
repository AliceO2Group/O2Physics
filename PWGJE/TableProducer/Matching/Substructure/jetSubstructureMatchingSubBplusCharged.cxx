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

// substructure matching event-wise subtracted Bplus charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <vector>
#include "PWGJE/TableProducer/Matching/Substructure/jetSubstructureMatchingSub.cxx"

using BplusChargedJetSubstructureMatchingSub = JetSubstructureMatchingSub<soa::Join<aod::BplusChargedJets, aod::BplusChargedJetConstituents, aod::BplusChargedJetsMatchedToBplusChargedEventWiseSubtractedJets>,
                                                                          soa::Join<aod::BplusChargedEventWiseSubtractedJets, aod::BplusChargedEventWiseSubtractedJetConstituents, aod::BplusChargedEventWiseSubtractedJetsMatchedToBplusChargedJets>,
                                                                          aod::BplusChargedSPsMatchedToBplusChargedEventWiseSubtractedSPs,
                                                                          aod::BplusChargedEventWiseSubtractedSPsMatchedToBplusChargedSPs,
                                                                          aod::BplusChargedPRsMatchedToBplusChargedEventWiseSubtractedPRs,
                                                                          aod::BplusChargedEventWiseSubtractedPRsMatchedToBplusChargedPRs,
                                                                          aod::BplusChargedSPs,
                                                                          aod::BplusChargedEventWiseSubtractedSPs,
                                                                          aod::BplusChargedPRs,
                                                                          aod::BplusChargedEventWiseSubtractedPRs,
                                                                          aod::CandidatesBplusData,
                                                                          aod::JetTracks,
                                                                          aod::JetTracksSubBplus,
                                                                          aod::JDummys>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<BplusChargedJetSubstructureMatchingSub>(cfgc, TaskName{"jet-substructure-matching-sub-bplus-ch"}));
  return WorkflowSpec{tasks};
}
