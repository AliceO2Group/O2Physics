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

// jet matching subtracted XicToXiPiPi charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubtraction.h"
#include "PWGJE/TableProducer/Matching/jetMatchingSub.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/runDataProcessing.h>

#include <vector>

using XicToXiPiPiChargedJetMatchingSub = JetMatchingSub<soa::Join<aod::XicToXiPiPiChargedJets, aod::XicToXiPiPiChargedJetConstituents>,
                                                        soa::Join<aod::XicToXiPiPiChargedEventWiseSubtractedJets, aod::XicToXiPiPiChargedEventWiseSubtractedJetConstituents>,
                                                        aod::XicToXiPiPiChargedJetsMatchedToXicToXiPiPiChargedEventWiseSubtractedJets,
                                                        aod::XicToXiPiPiChargedEventWiseSubtractedJetsMatchedToXicToXiPiPiChargedJets,
                                                        aod::JTrackXicToXiPiPiSubs,
                                                        aod::CandidatesXicToXiPiPiData>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<XicToXiPiPiChargedJetMatchingSub>(cfgc, TaskName{"jet-matching-sub-xictoxipipi-ch"}));
  return WorkflowSpec{tasks};
}
