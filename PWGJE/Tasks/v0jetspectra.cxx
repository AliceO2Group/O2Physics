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

// jet spectra for v0 fragmentation study
//
/// \author Gijs van Weelden <g.van.weelden@cern.ch>
//

#include "TH1F.h"
#include "TTree.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/RunningWorkflowInfo.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

#include "CommonConstants/PhysicsConstants.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MCDJets = aod::ChargedMCDetectorLevelJets;
using MCDJetsWithConstituents = soa::Join<MCDJets, aod::ChargedMCDetectorLevelJetConstituents>;
using MatchedMCDJets = soa::Join<MCDJets, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>;
using MatchedMCDJetsWithConstituents = soa::Join<MCDJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>;

using MCPJets = aod::ChargedMCParticleLevelJets;
using MatchedMCPJets = soa::Join<MCPJets, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>;
using MCPJetsWithConstituents = soa::Join<MCPJets, aod::ChargedMCParticleLevelJetConstituents>;
using MatchedMCPJetsWithConstituents = soa::Join<MCPJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>;

// V0 jets
using MCDV0Jets = aod::V0ChargedMCDetectorLevelJets;
using MCDV0JetsWithConstituents = soa::Join<MCDV0Jets, aod::V0ChargedMCDetectorLevelJetConstituents>;
using MatchedMCDV0Jets = soa::Join<MCDV0Jets, aod::V0ChargedMCDetectorLevelJetsMatchedToV0ChargedMCParticleLevelJets>;
using MatchedMCDV0JetsWithConstituents = soa::Join<MCDV0Jets, aod::V0ChargedMCDetectorLevelJetConstituents, aod::V0ChargedMCDetectorLevelJetsMatchedToV0ChargedMCParticleLevelJets>;

using MCPV0Jets = aod::V0ChargedMCParticleLevelJets;
using MCPV0JetsWithConstituents = soa::Join<MCPV0Jets, aod::V0ChargedMCParticleLevelJetConstituents>;
using MatchedMCPV0Jets = soa::Join<MCPV0Jets, aod::V0ChargedMCParticleLevelJetsMatchedToV0ChargedMCDetectorLevelJets>;
using MatchedMCPV0JetsWithConstituents = soa::Join<MCPV0Jets, aod::V0ChargedMCParticleLevelJetConstituents, aod::V0ChargedMCParticleLevelJetsMatchedToV0ChargedMCDetectorLevelJets>;

struct V0JetSpectra {
  HistogramRegistry registry{"registry"};

  Configurable<std::string> evSel{"evSel", "sel8WithoutTimeFrameBorderCut", "choose event selection"};
  Configurable<float> vertexZCut{"vertexZCut", 10.f, "vertex z cut"};
  int eventSelection = -1;

  Filter jetCollisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;

  void init(InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(evSel));
    registry.add("jetPtEtaPhi", "Jets; #it{p}_{T}; #eta; #phi", HistType::kTH3D, {{200, 0., 200.}, {20, -1.f, 1.f}, {18 * 8, 0.f, 2. * TMath::Pi()}});
    registry.add("mcpJetPtEtaPhi", "Jets; #it{p}_{T}; #eta; #phi", HistType::kTH3D, {{200, 0., 200.}, {20, -1.f, 1.f}, {18 * 8, 0.f, 2. * TMath::Pi()}});
  }

  template <typename T>
  void fillHistograms(T const& jets, double weight = 1.)
  {
    for (const auto& jet : jets) {
      registry.fill(HIST("jetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi(), weight);
    }
  }
  template <typename T>
  void fillMCPHistograms(T const& jets, double weight = 1.)
  {
    for (const auto& jet : jets) {
      registry.fill(HIST("mcpJetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi(), weight);
    }
  }

  void processData(soa::Filtered<JetCollisions>::iterator const& jcoll, aod::ChargedJets const& chjets, aod::V0ChargedJets const& v0jets)
  {
    if (!jetderiveddatautilities::selectCollision(jcoll, eventSelection)) {
      return;
    }
    if (v0jets.size() > 0) {
      return;
    } else {
      fillHistograms(chjets);
    }
  }
  PROCESS_SWITCH(V0JetSpectra, processData, "Jet spectra for V0 jets or Ch jets if no V0s in data", false);

  void processMCD(soa::Filtered<JetCollisionsMCD>::iterator const& jcoll, JetMcCollisions const&, MCDJets const& chjets, MCDV0Jets const& v0jets)
  {
    if (!jcoll.has_mcCollision()) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(jcoll, eventSelection)) {
      return;
    }
    double weight = jcoll.mcCollision().weight();
    if (v0jets.size() > 0) {
      return;
    } else {
      fillHistograms(chjets, weight);
    }
  }
  PROCESS_SWITCH(V0JetSpectra, processMCD, "Jet spectra for V0 jets or Ch jets if no V0s", false);

  void processMCP(JetMcCollision const& jcoll, MCPJets const& chjets, MCPV0Jets const& v0jets)
  {
    double weight = jcoll.weight();
    if (v0jets.size() > 0) {
      return;
    } else {
      fillMCPHistograms(chjets, weight);
    }
  }
  PROCESS_SWITCH(V0JetSpectra, processMCP, "Jet spectra for V0 jets or Ch jets if no V0s", false);

  void processMcMatched(soa::Filtered<JetCollisionsMCD>::iterator const& jcoll, JetMcCollisions const&, MatchedMCDJets const& chjetsMCD, MatchedMCPJets const& chjetsMCP, MatchedMCDV0Jets const& v0jetsMCD, MatchedMCPV0Jets const& v0jetsMCP)
  {
    if (!jetderiveddatautilities::selectCollision(jcoll, eventSelection)) {
      return;
    }
    // TODO: Need to add checker to only count matched jets (?)
    double weight = jcoll.mcCollision().weight();
    if (v0jetsMCP.size() > 0) {
      return;
    } else {
      fillMCPHistograms(chjetsMCP, weight);
    } // Particle level loop

    if (v0jetsMCD.size() > 0) {
      return;
    } else {
      fillHistograms(chjetsMCD, weight);
    } // Detector level loop
  }
  PROCESS_SWITCH(V0JetSpectra, processMcMatched, "Jet spectra for matched ch jets if no V0s", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<V0JetSpectra>(cfgc, TaskName{"jet-v0-spectra"})};
}
