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

// jet finder D+ data charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/JetFinders/jetFinderHF.cxx"

using JetFinderDplusDataCharged = JetFinderHFTask<aod::CandidatesDplusData, aod::CandidatesDplusMCD, aod::CandidatesDplusMCP, aod::JetTracksSubDplus, aod::JetParticlesSubDplus, aod::DplusChargedJets, aod::DplusChargedJetConstituents, aod::DplusChargedEventWiseSubtractedJets, aod::DplusChargedEventWiseSubtractedJetConstituents>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetFinderDplusDataCharged>(cfgc,
                                                                  SetDefaultProcesses{{{"processChargedJetsData", true}}},
                                                                  TaskName{"jet-finder-dplus-data-charged"}));

  return WorkflowSpec{tasks};
}
