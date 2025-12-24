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

// jet tutorial task for hands on tutorial session (09/11/2023)
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
//

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct SubstructureDebugTask {
  HistogramRegistry registry{"registry",
                             {{"h_collisions_posZ", "event posZ;posZ;entries", {HistType::kTH1F, {{20, -10.0, 10.0}}}},
                              {"h_mccollisions_posZ", "mc event posZ;posZ;entries", {HistType::kTH1F, {{20, -10.0, 10.0}}}},
                              {"h_jet_pt", "track pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_jet_phi", "jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_jet_nTracks", "jet tracks;n_{track};entries", {HistType::kTH1F, {{80, 0.0, 80.0}}}},
                              {"h_jet_angularity", "jet angularity;angularity;entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
                              {"h_mcpjet_pt", "mcp jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_mcpjet_eta", "mcp jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_mcpjet_phi", "mcp jet #varphi;#varphi_{jeet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_mcpjet_nTracks", "mcp jet tracks,n_{track};entries", {HistType::kTH1F, {{80, 0.0, 80.0}}}},
                              {"h_mcpjet_angularity", "mcp jet angularity;angularity;entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
                              {"h_jet_pt_matched", "jet pT matched;#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}}},
                              {"h_D0injet_pt", "D0 in jet pT;#it{p}_{T,D0} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_D0jet_pt", "D0 jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_D0jet_eta", "D0 jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_D0jet_phi", "D0 jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_D0jet_nTracks", "D0 jet n tracks;n_{track};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_D0jet_angularity", "D0 jet angularity;angularity;entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_D0jet_pt_matched", "D0 jet pT matched;#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}}},
                              {"h_mcpD0injet_pt", "mcp D0 in jet pT;#it{p}_{T,D0} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_mcpD0jet_pt", "mcp D0 jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_mcpD0jet_eta", "mcp D0 jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_mcpD0jet_phi", "mcp D0 jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_mcpD0jet_nTracks", "mcp D0 jet ntracks;n_{track};entries", {HistType::kTH1F, {{80, 0.0, 80.}}}},
                              {"h_mcpD0jet_angularity", "mcp D0 jet angularity;angularity;entries", {HistType::kTH1F, {{50, 0.0, 1.}}}},
                              {"h_Lcinjet_pt", "Lc in jet pT;#it{p}_{T,Lc} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_Lcjet_pt", "Lc jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_Lcjet_eta", "Lc jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_Lcjet_phi", "Lc jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_Lcjet_nTracks", "Lc jet n tracks;n_{track};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_Lcjet_angularity", "Lc jet angularity;angularity;entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_Lcjet_pt_matched", "Lc jet pT matched;#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}}},
                              {"h_mcpLcinjet_pt", "mcp Lc in jet pT;#it{p}_{T,Lc} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_mcpLcjet_pt", "mcp Lc jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_mcpLcjet_eta", "mcp Lc jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_mcpLcjet_phi", "mcp Lc jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_mcpLcjet_nTracks", "mcp Lc jet ntracks;n_{track};entries", {HistType::kTH1F, {{80, 0.0, 80.}}}},
                              {"h_mcpLcjet_angularity", "mcp Lc jet angularity;angularity;entries", {HistType::kTH1F, {{50, 0.0, 1.}}}}}};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  std::vector<int> eventSelection;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
  }

  void processDataCharged(aod::CJetCOs::iterator const& collision, soa::Join<aod::CJetOs, aod::CJetSSOs> const& jets)
  {
    registry.fill(HIST("h_collisions_posZ"), collision.posZ());
    for (auto const& jet : jets) {
      registry.fill(HIST("h_jet_pt"), jet.jetPt());
      registry.fill(HIST("h_jet_phi"), jet.jetPhi());
      registry.fill(HIST("h_jet_eta"), jet.jetEta());
      registry.fill(HIST("h_jet_nTracks"), jet.jetNConstituents());
      registry.fill(HIST("h_jet_angularity"), jet.angularity());
    }
  }
  PROCESS_SWITCH(SubstructureDebugTask, processDataCharged, "jets data", true);

  void processDataD0(aod::D0CJetCOs::iterator const&, soa::Join<aod::D0CJetOs, aod::D0CJetSSOs> const& jets, aod::CandidatesD0Data const&)
  {
    for (auto const& jet : jets) {
      registry.fill(HIST("h_D0jet_pt"), jet.jetPt());
      registry.fill(HIST("h_D0jet_phi"), jet.jetPhi());
      registry.fill(HIST("h_D0jet_eta"), jet.jetEta());
      registry.fill(HIST("h_D0jet_nTracks"), jet.jetNConstituents());
      registry.fill(HIST("h_D0jet_angularity"), jet.angularity());
      for (auto const& candidate : jet.candidates_as<aod::CandidatesD0Data>()) {
        registry.fill(HIST("h_D0injet_pt"), candidate.pt());
      }
    }
  }
  PROCESS_SWITCH(SubstructureDebugTask, processDataD0, "jets D0 data", true);

  void processDataLc(aod::LcCJetCOs::iterator const&, soa::Join<aod::LcCJetOs, aod::LcCJetSSOs> const& jets, aod::CandidatesLcData const&)
  {
    for (auto const& jet : jets) {
      registry.fill(HIST("h_Lcjet_pt"), jet.jetPt());
      registry.fill(HIST("h_Lcjet_phi"), jet.jetPhi());
      registry.fill(HIST("h_Lcjet_eta"), jet.jetEta());
      registry.fill(HIST("h_Lcjet_nTracks"), jet.jetNConstituents());
      registry.fill(HIST("h_Lcjet_angularity"), jet.angularity());
      for (auto const& candidate : jet.candidates_as<aod::CandidatesLcData>()) {
        registry.fill(HIST("h_Lcinjet_pt"), candidate.pt());
      }
    }
  }
  PROCESS_SWITCH(SubstructureDebugTask, processDataLc, "jets Lc data", true);

  void processMCDCharged(aod::CMCDJetCOs::iterator const& collision, soa::Join<aod::CMCDJetOs, aod::CMCDJetSSOs, aod::CMCDJetMOs> const& jets, aod::CMCPJetOs const&)
  {
    registry.fill(HIST("h_collisions_posZ"), collision.posZ());
    for (auto const& jet : jets) {
      registry.fill(HIST("h_jet_pt"), jet.jetPt());
      registry.fill(HIST("h_jet_phi"), jet.jetPhi());
      registry.fill(HIST("h_jet_eta"), jet.jetEta());
      registry.fill(HIST("h_jet_nTracks"), jet.jetNConstituents());
      registry.fill(HIST("h_jet_angularity"), jet.angularity());
      auto const& matchedJets = jet.matchedJetGeo_as<aod::CMCPJetOs>();
      for (auto const& mcpjet : matchedJets) {
        registry.fill(HIST("h_jet_pt_matched"), mcpjet.jetPt(), jet.jetPt());
      }
    }
  }
  PROCESS_SWITCH(SubstructureDebugTask, processMCDCharged, "jets mcd", false);

  void processMCDD0(aod::D0CMCDJetCOs::iterator const&, soa::Join<aod::D0CMCDJetOs, aod::D0CMCDJetSSOs, aod::D0CMCDJetMOs> const& jets, aod::D0CMCPJetOs const&, aod::CandidatesD0MCD const&)
  {
    for (auto const& jet : jets) {
      registry.fill(HIST("h_D0jet_pt"), jet.jetPt());
      registry.fill(HIST("h_D0jet_phi"), jet.jetPhi());
      registry.fill(HIST("h_D0jet_eta"), jet.jetEta());
      registry.fill(HIST("h_D0jet_nTracks"), jet.jetNConstituents());
      registry.fill(HIST("h_D0jet_angularity"), jet.angularity());
      auto const& matchedJets = jet.matchedJetGeo_as<aod::D0CMCPJetOs>();
      for (auto const& mcpjet : matchedJets) {
        registry.fill(HIST("h_D0jet_pt_matched"), mcpjet.jetPt(), jet.jetPt());
      }
      for (auto const& candidate : jet.candidates_as<aod::CandidatesD0MCD>()) {
        registry.fill(HIST("h_D0injet_pt"), candidate.pt());
      }
    }
  }
  PROCESS_SWITCH(SubstructureDebugTask, processMCDD0, "jets D0 mcd", false);

  void processMCDLc(aod::LcCMCDJetCOs::iterator const&, soa::Join<aod::LcCMCDJetOs, aod::LcCMCDJetSSOs, aod::LcCMCDJetMOs> const& jets, aod::LcCMCPJetOs const&, aod::CandidatesLcMCD const&)
  {
    for (auto const& jet : jets) {
      registry.fill(HIST("h_Lcjet_pt"), jet.jetPt());
      registry.fill(HIST("h_Lcjet_phi"), jet.jetPhi());
      registry.fill(HIST("h_Lcjet_eta"), jet.jetEta());
      registry.fill(HIST("h_Lcjet_nTracks"), jet.jetNConstituents());
      registry.fill(HIST("h_Lcjet_angularity"), jet.angularity());
      auto const& matchedJets = jet.matchedJetGeo_as<aod::LcCMCPJetOs>();
      for (auto const& mcpjet : matchedJets) {
        registry.fill(HIST("h_Lcjet_pt_matched"), mcpjet.jetPt(), jet.jetPt());
      }
      for (auto const& candidate : jet.candidates_as<aod::CandidatesLcMCD>()) {
        registry.fill(HIST("h_Lcinjet_pt"), candidate.pt());
      }
    }
  }
  PROCESS_SWITCH(SubstructureDebugTask, processMCDLc, "jets Lc mcd", false);

  void processMCPCharged(aod::CMCPJetCOs::iterator const& mcCollision, soa::Join<aod::CMCPJetOs, aod::CMCPJetSSOs, aod::CMCPJetMOs> const& jets)
  {
    registry.fill(HIST("h_mccollisions_posZ"), mcCollision.posZ());
    for (auto const& jet : jets) {
      registry.fill(HIST("h_mcpjet_pt"), jet.jetPt());
      registry.fill(HIST("h_mcpjet_phi"), jet.jetPhi());
      registry.fill(HIST("h_mcpjet_eta"), jet.jetEta());
      registry.fill(HIST("h_mcpjet_nTracks"), jet.jetNConstituents());
      registry.fill(HIST("h_mcpjet_angularity"), jet.angularity());
    }
  }
  PROCESS_SWITCH(SubstructureDebugTask, processMCPCharged, "jets mcp", false);

  void processMCPD0(aod::D0CMCPJetCOs::iterator const&, soa::Join<aod::D0CMCPJetOs, aod::D0CMCPJetSSOs, aod::D0CMCPJetMOs> const& jets, aod::CandidatesD0MCP const&)
  {
    for (auto const& jet : jets) {
      registry.fill(HIST("h_mcpD0jet_pt"), jet.jetPt());
      registry.fill(HIST("h_mcpD0jet_phi"), jet.jetPhi());
      registry.fill(HIST("h_mcpD0jet_eta"), jet.jetEta());
      registry.fill(HIST("h_mcpD0jet_nTracks"), jet.jetNConstituents());
      registry.fill(HIST("h_mcpD0jet_angularity"), jet.angularity());
      for (auto const& candidate : jet.candidates_as<aod::CandidatesD0MCP>()) {
        registry.fill(HIST("h_mcpD0injet_pt"), candidate.pt());
      }
    }
  }
  PROCESS_SWITCH(SubstructureDebugTask, processMCPD0, "jets D0 mcp", false);

  void processMCPLc(aod::LcCMCPJetCOs::iterator const&, soa::Join<aod::LcCMCPJetOs, aod::LcCMCPJetSSOs, aod::LcCMCPJetMOs> const& jets, aod::CandidatesLcMCP const&)
  {
    for (auto const& jet : jets) {
      registry.fill(HIST("h_mcpLcjet_pt"), jet.jetPt());
      registry.fill(HIST("h_mcpLcjet_phi"), jet.jetPhi());
      registry.fill(HIST("h_mcpLcjet_eta"), jet.jetEta());
      registry.fill(HIST("h_mcpLcjet_nTracks"), jet.jetNConstituents());
      registry.fill(HIST("h_mcpLcjet_angularity"), jet.angularity());
      for (auto const& candidate : jet.candidates_as<aod::CandidatesLcMCP>()) {
        registry.fill(HIST("h_mcpLcinjet_pt"), candidate.pt());
      }
    }
  }
  PROCESS_SWITCH(SubstructureDebugTask, processMCPLc, "jets Lc mcp", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<SubstructureDebugTask>(cfgc, TaskName{"substructure-debug"})}; }
