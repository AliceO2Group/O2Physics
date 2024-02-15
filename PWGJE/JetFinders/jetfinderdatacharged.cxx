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

// jet finder data charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/JetFinders/jetfinder.cxx"

using JetFinderDataCharged = JetFinderTask<aod::ChargedJets, aod::ChargedJetConstituents, aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(
    adaptAnalysisTask<JetFinderDataCharged>(cfgc,
                                            SetDefaultProcesses{{{"processChargedJets", true}}}, TaskName{"jet-finder-data-charged"}));

  return WorkflowSpec{tasks};
}
