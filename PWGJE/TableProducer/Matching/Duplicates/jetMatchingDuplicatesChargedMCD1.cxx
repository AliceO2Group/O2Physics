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

// jet matching duplicates charged mcd task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/TableProducer/Matching/Duplicates/jetMatchingDuplicates.cxx"

using Charged1JetMCDMatchingDupliacates = JetMatchingDuplicates<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>,
                                                                soa::Join<aod::Charged1MCDetectorLevelJets, aod::Charged1MCDetectorLevelJetConstituents>,
                                                                aod::ChargedMCDetectorLevelJetsMatchedToCharged1MCDetectorLevelJets,
                                                                aod::Charged1MCDetectorLevelJetsMatchedToChargedMCDetectorLevelJets,
                                                                aod::JTracks,
                                                                aod::JDummys>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<Charged1JetMCDMatchingDupliacates>(cfgc, SetDefaultProcesses{}, TaskName{"jet-matching-mcd-ch-1"}));
  return WorkflowSpec{tasks};
}
