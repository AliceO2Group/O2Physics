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

// jet finder dielectron mcd charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/JetFinders/jetFinderHF.cxx"

using JetFinderDielectronMCDetectorLevelCharged = JetFinderHFTask<aod::CandidatesDielectronData, aod::CandidatesDielectronMCD, aod::CandidatesDielectronMCP, aod::JetTracksSubDielectron, aod::JetParticlesSubDielectron, aod::DielectronChargedMCDetectorLevelJets, aod::DielectronChargedMCDetectorLevelJetConstituents, aod::DielectronChargedMCDetectorLevelEventWiseSubtractedJets, aod::DielectronChargedMCDetectorLevelEventWiseSubtractedJetConstituents>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetFinderDielectronMCDetectorLevelCharged>(cfgc,
                                                                                  SetDefaultProcesses{{{"processChargedJetsMCD", true}}},
                                                                                  TaskName{"jet-finder-dielectron-mcd-charged"}));

  return WorkflowSpec{tasks};
}
