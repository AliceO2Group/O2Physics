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

// jet matching subtracted Dielectron charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/TableProducer/Matching/jetMatchingSub.cxx"

using DielectronChargedJetMatchingSub = JetMatchingSub<soa::Join<aod::DielectronChargedJets, aod::DielectronChargedJetConstituents>,
                                                       soa::Join<aod::DielectronChargedEventWiseSubtractedJets, aod::DielectronChargedEventWiseSubtractedJetConstituents>,
                                                       aod::DielectronChargedJetsMatchedToDielectronChargedEventWiseSubtractedJets,
                                                       aod::DielectronChargedEventWiseSubtractedJetsMatchedToDielectronChargedJets,
                                                       aod::JTrackDielectronSubs,
                                                       aod::CandidatesDielectronData>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<DielectronChargedJetMatchingSub>(cfgc, TaskName{"jet-matching-sub-dielectron-ch"}));
  return WorkflowSpec{tasks};
}
