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

// substructure matching event-wise subtracted B0 charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/TableProducer/Matching/Substructure/jetSubstructureMatchingSub.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/runDataProcessing.h>

#include <vector>

using B0ChargedJetSubstructureMatchingSub = JetSubstructureMatchingSub<soa::Join<aod::B0ChargedJets, aod::B0ChargedJetConstituents, aod::B0ChargedJetsMatchedToB0ChargedEventWiseSubtractedJets>,
                                                                       soa::Join<aod::B0ChargedEventWiseSubtractedJets, aod::B0ChargedEventWiseSubtractedJetConstituents, aod::B0ChargedEventWiseSubtractedJetsMatchedToB0ChargedJets>,
                                                                       aod::B0ChargedSPsMatchedToB0ChargedEventWiseSubtractedSPs,
                                                                       aod::B0ChargedEventWiseSubtractedSPsMatchedToB0ChargedSPs,
                                                                       aod::B0ChargedPRsMatchedToB0ChargedEventWiseSubtractedPRs,
                                                                       aod::B0ChargedEventWiseSubtractedPRsMatchedToB0ChargedPRs,
                                                                       aod::B0ChargedSPs,
                                                                       aod::B0ChargedEventWiseSubtractedSPs,
                                                                       aod::B0ChargedPRs,
                                                                       aod::B0ChargedEventWiseSubtractedPRs,
                                                                       aod::CandidatesB0Data,
                                                                       aod::JetTracks,
                                                                       aod::JetTracksSubB0,
                                                                       aod::JDummys>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<B0ChargedJetSubstructureMatchingSub>(cfgc, TaskName{"jet-substructure-matching-sub-b0-ch"}));
  return WorkflowSpec{tasks};
}
