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

// jet substructure Lc charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/Tasks/jetSubstructureHF.cxx"

using JetSubstructureLc = JetSubstructureHFTask<soa::Join<aod::LcChargedJets, aod::LcChargedJetConstituents>, soa::Join<aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents>, soa::Join<aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents>, soa::Join<aod::LcChargedEventWiseSubtractedJets, aod::LcChargedEventWiseSubtractedJetConstituents>, aod::CandidatesLcData, aod::CandidatesLcMCP, aod::LcCJetSSs, aod::LcCMCDJetSSs, aod::LcCMCPJetSSs, aod::LcCEWSJetSSs, aod::JTrackLcSubs>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureLc>(cfgc,
                                                          SetDefaultProcesses{},
                                                          TaskName{"jet-substructure-lc"}));
  return WorkflowSpec{tasks};
}
