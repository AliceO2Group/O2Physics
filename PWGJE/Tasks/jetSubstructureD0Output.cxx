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

// jet substructure output D0 charged task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/Tasks/jetSubstructureHFOutput.cxx"

using JetSubstructureOutputD0 = JetSubstructureHFOutputTask<aod::CollisionsD0, soa::Join<aod::McCollisionsD0, aod::HfD0McRCollIds>, aod::CandidatesD0Data, aod::CandidatesD0MCD, aod::CandidatesD0MCP, aod::JTrackD0Subs, soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents, aod::D0CJetSSs>, soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents, aod::D0CJetSSs, aod::D0ChargedJetsMatchedToD0ChargedEventWiseSubtractedJets>, aod::D0CJetCOs, aod::D0CJetOs, aod::D0CJetSSOs, aod::D0CJetMOs, soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0CMCDJetSSs, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets>, aod::D0CMCDJetCOs, aod::D0CMCDJetOs, aod::D0CMCDJetSSOs, aod::D0CMCDJetMOs, soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0CMCPJetSSs, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets>, aod::D0CMCPJetCOs, aod::D0CMCPJetOs, aod::D0CMCPJetSSOs, aod::D0CMCPJetMOs, soa::Join<aod::D0ChargedEventWiseSubtractedJets, aod::D0ChargedEventWiseSubtractedJetConstituents, aod::D0CEWSJetSSs, aod::D0ChargedEventWiseSubtractedJetsMatchedToD0ChargedJets>, aod::D0CEWSJetCOs, aod::D0CEWSJetOs, aod::D0CEWSJetSSOs, aod::D0CEWSJetMOs, aod::StoredHfD0CollBase, aod::StoredHfD0Bases, aod::StoredHfD0Pars, aod::StoredHfD0ParEs, aod::JDumD0ParDaus, aod::StoredHfD0Sels, aod::StoredHfD0Mls, aod::JDumD0MlDaus, aod::StoredHfD0Mcs, aod::StoredHfD0McCollBases, aod::StoredHfD0McRCollIds, aod::StoredHfD0PBases>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputD0>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructure-d0-output"}));

  return WorkflowSpec{tasks};
}
