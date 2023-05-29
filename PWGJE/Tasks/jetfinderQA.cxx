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

// jet finder hf QA task
//
// Authors: Nima Zardoshti

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Common/Core/RecoDecay.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_2prong;

#include "Framework/runDataProcessing.h"

struct JetFinderQATask {

  AxisSpec weightAxis = {{VARIABLE_WIDTH, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0}, "weights"};
  HistogramRegistry registry{"registry",
                             {{"h_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}, true},
                              {"h_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}, true},
                              {"h_jet_phi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}, true},
                              {"h_jet_ntracks", "jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {{100, -0.5, 99.5}}}, true},
                              {"h2_jet_pt_jet_eta", "#it{p}_{T,jet} (GeV/#it{c}); #eta_{jet}", {HistType::kTH2F, {{200, 0.0, 200}, {100, -1.0, 1.0}}}, true},
                              {"h2_jet_pt_jet_phi", "#it{p}_{T,jet} (GeV/#it{c}); #phi_{jet}", {HistType::kTH2F, {{200, 0.0, 200}, {80, -1.0, 7.}}}, true},
                              {"h2_jet_pt_jet_ntracks", "#it{p}_{T,jet} (GeV/#it{c}); N_{jet tracks}", {HistType::kTH2F, {{200, 0.0, 200}, {100, -0.5, 99.5}}}, true},
                              {"h2_jet_r_jet_pt", "#it{R}_{jet}; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{10, 0.05, 1.05}, {200, 0.0, 200}}}, true},
                              {"h2_jet_r_jet_eta", "#it{R}_{jet}; #eta_{jet}", {HistType::kTH2F, {{10, 0.05, 1.05}, {100, -1.0, 1.0}}}, true},
                              {"h2_jet_r_jet_phi", "#it{R}_{jet}; #phi_{jet}", {HistType::kTH2F, {{10, 0.05, 1.05}, {80, -1.0, 7.}}}, true},
                              {"h2_jet_r_jet_ntracks", "#it{R}_{jet}; N_{jet tracks}", {HistType::kTH2F, {{10, 0.05, 1.05}, {100, -0.5, 99.5}}}, true},
                              {"h2_jet_pt_track_pt", "#it{p}_{T,jet} (GeV/#it{c}); #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{200, 0.0, 200}, {200, 0.0, 200.0}}}, true},
                              {"h2_jet_pt_track_eta", "#it{p}_{T,jet} (GeV/#it{c}); #eta_{track}", {HistType::kTH2F, {{200, 0.0, 200}, {100, -1.0, 1.0}}}, true},
                              {"h2_jet_pt_track_phi", "#it{p}_{T,jet} (GeV/#it{c}); #phi_{track}", {HistType::kTH2F, {{200, 0.0, 200}, {80, -1.0, 7.}}}, true},
                              {"h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0.0, 200.0}}}, true},
                              {"h_track_eta", "track #eta;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}, true},
                              {"h_track_phi", "track #phi;#phi_{track};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}, true},
                              {"h_jet_pt_part", "jet pT;#it{p}_{T,jet}^{part}(GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}, true},
                              {"h_jet_eta_part", "jet #eta;#eta_{jet}^{part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}, true},
                              {"h_jet_phi_part", "jet #phi;#phi_{jet}^{part};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}, true},
                              {"h_jet_ntracks_part", "jet N tracks;N_{jet tracks}^{part};entries", {HistType::kTH1F, {{100, -0.5, 99.5}}}, true},
                              {"h2_jet_pt_part_jet_eta_part", "#it{p}_{T,jet}^{part} (GeV/#it{c}); #eta_{jet}", {HistType::kTH2F, {{200, 0.0, 200}, {100, -1.0, 1.0}}}, true},
                              {"h2_jet_pt_part_jet_phi_part", "#it{p}_{T,jet}^{part} (GeV/#it{c}); #phi_{jet}", {HistType::kTH2F, {{200, 0.0, 200}, {80, -1.0, 7.}}}, true},
                              {"h2_jet_pt_part_jet_ntracks_part", "#it{p}_{T,jet}^{part} (GeV/#it{c}); N_{jet tracks}", {HistType::kTH2F, {{200, 0.0, 200}, {100, -0.5, 99.5}}}, true},
                              {"h2_jet_r_part_jet_pt_part", "#it{R}_{jet}^{part}; #it{p}_{T,jet}^{part} (GeV/#it{c})", {HistType::kTH2F, {{10, 0.05, 1.05}, {200, 0.0, 200}}}, true},
                              {"h2_jet_r_part_jet_eta_part", "#it{R}_{jet}^{part}; #eta_{jet}^{part}", {HistType::kTH2F, {{10, 0.05, 1.05}, {100, -1.0, 1.0}}}, true},
                              {"h2_jet_r_part_jet_phi_part", "#it{R}_{jet}^{part}; #phi_{jet}^{part}", {HistType::kTH2F, {{10, 0.05, 1.05}, {80, -1.0, 7.}}}, true},
                              {"h2_jet_r_part_jet_ntracks_part", "#it{R}_{jet}^{part}; N_{jet tracks}^{part}", {HistType::kTH2F, {{10, 0.05, 1.05}, {100, -0.5, 99.5}}}, true},
                              {"h2_jet_pt_part_track_pt_part", "#it{p}_{T,jet}^{part} (GeV/#it{c}); #it{p}_{T,track}^{part} (GeV/#it{c})", {HistType::kTH2F, {{200, 0.0, 200}, {200, 0.0, 200.0}}}, true},
                              {"h2_jet_pt_part_track_eta_part", "#it{p}_{T,jet}^{part} (GeV/#it{c}); #eta_{track}^{part}", {HistType::kTH2F, {{200, 0.0, 200}, {100, -1.0, 1.0}}}, true},
                              {"h2_jet_pt_part_track_phi_part", "#it{p}_{T,jet}^{part} (GeV/#it{c}); #phi_{track}^{part}", {HistType::kTH2F, {{200, 0.0, 200}, {80, -1.0, 7.}}}, true},
                              {"h_track_pt_part", "track pT;#it{p}_{T,track}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0.0, 200.0}}}, true},
                              {"h_track_eta_part", "track #eta;#eta_{track}^{part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}, true},
                              {"h_track_phi_part", "track #phi;#phi_{track}^{part};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}, true},
                              {"h2_jet_pt_part_jet_pt", "#it{p}_{T,jet}^{part} (GeV/#it{c}); #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{200, 0.0, 200}, {200, 0.0, 200}}}, true},
                              {"h2_jet_eta_part_jet_eta", "#eta_{jet}^{part}; #eta_{jet}", {HistType::kTH2F, {{100, -1.0, 1.0}, {100, -1.0, 1.0}}}, true},
                              {"h2_jet_phi_part_jet_phi", "#phi_{jet}^{part}; #phi_{jet}", {HistType::kTH2F, {{80, -1.0, 7.}, {80, -1.0, 7.}}}, true},
                              {"h2_jet_ntracks_part_jet_ntracks", "N_{jet tracks}^{part}; N_{jet tracks}", {HistType::kTH2F, {{100, -0.5, 99.5}, {100, -0.5, 99.5}}}, true},
                              {"h3_jet_pt_part_jet_eta_part_jet_eta", "#it{p}_{T,jet}^{part} (GeV/#it{c}); #eta_{jet}^{part}; #eta_{jet}", {HistType::kTH3F, {{200, 0.0, 200}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}}, true},
                              {"h3_jet_pt_part_jet_phi_part_jet_phi", "#it{p}_{T,jet}^{part} (GeV/#it{c}); #phi_{jet}^{part}; #phi_{jet}", {HistType::kTH3F, {{200, 0.0, 200}, {80, -1.0, 7.}, {80, -1.0, 7.}}}, true},
                              {"h3_jet_pt_part_jet_ntracks_part_jet_ntracks", "#it{p}_{T,jet}^{part} (GeV/#it{c}); N_{jet tracks}^{part}; N_{jet tracks}", {HistType::kTH3F, {{200, 0.0, 200}, {100, -0.5, 99.5}, {100, -0.5, 99.5}}}, true},
                              {"h2_jet_pt_part_jet_pt_diff", "#it{p}_{T,jet}^{part} (GeV/#it{c}); (#it{p}_{T,jet}^{part} (GeV/#it{c}) - #it{p}_{T,jet} (GeV/#it{c})) / #it{p}_{T,jet}^{part} (GeV/#it{c})", {HistType::kTH2F, {{200, 0.0, 200}, {1000, -5.0, 5.0}}}, true},
                              {"h_collision_eventweight_part", "event weight;event weight;entries", {HistType::kTH1F, {weightAxis}}}}};

  void init(o2::framework::InitContext&)
  {
  }

  using JetTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>;

  template <typename T>
  void fillHistograms(T const& jet, float weight = 1.0)
  {

    registry.fill(HIST("h_jet_pt"), jet.pt(), weight);
    registry.fill(HIST("h_jet_eta"), jet.eta(), weight);
    registry.fill(HIST("h_jet_phi"), jet.phi(), weight);
    registry.fill(HIST("h_jet_ntracks"), jet.tracks().size(), weight);

    registry.fill(HIST("h2_jet_pt_jet_eta"), jet.pt(), jet.eta(), weight);
    registry.fill(HIST("h2_jet_pt_jet_phi"), jet.pt(), jet.phi(), weight);
    registry.fill(HIST("h2_jet_pt_jet_ntracks"), jet.pt(), jet.tracks().size(), weight);

    registry.fill(HIST("h2_jet_r_jet_pt"), jet.r() / 100.0, jet.pt(), weight);
    registry.fill(HIST("h2_jet_r_jet_eta"), jet.r() / 100.0, jet.eta(), weight);
    registry.fill(HIST("h2_jet_r_jet_phi"), jet.r() / 100.0, jet.phi(), weight);
    registry.fill(HIST("h2_jet_r_jet_ntracks"), jet.r() / 100.0, jet.tracks().size(), weight);

    for (auto& constituent : jet.template tracks_as<JetTracks>()) {

      registry.fill(HIST("h2_jet_pt_track_pt"), jet.pt(), constituent.pt(), weight);
      registry.fill(HIST("h2_jet_pt_track_eta"), jet.pt(), constituent.eta(), weight);
      registry.fill(HIST("h2_jet_pt_track_phi"), jet.pt(), constituent.phi(), weight);
      registry.fill(HIST("h_track_pt"), constituent.pt(), weight);
      registry.fill(HIST("h_track_eta"), constituent.eta(), weight);
      registry.fill(HIST("h_track_phi"), constituent.phi(), weight);
    }
  }

  template <typename T>
  void fillMCPHistograms(T const& jet, float weight = 1.0)
  {

    registry.fill(HIST("h_jet_pt_part"), jet.pt(), weight);
    registry.fill(HIST("h_jet_eta_part"), jet.eta(), weight);
    registry.fill(HIST("h_jet_phi_part"), jet.phi(), weight);
    registry.fill(HIST("h_jet_ntracks_part"), jet.tracks().size(), weight);

    registry.fill(HIST("h2_jet_pt_part_jet_eta_part"), jet.pt(), jet.eta(), weight);
    registry.fill(HIST("h2_jet_pt_part_jet_phi_part"), jet.pt(), jet.phi(), weight);
    registry.fill(HIST("h2_jet_pt_part_jet_ntracks_part"), jet.pt(), jet.tracks().size(), weight);

    registry.fill(HIST("h2_jet_r_part_jet_pt_part"), jet.r() / 100.0, jet.pt(), weight);
    registry.fill(HIST("h2_jet_r_part_jet_eta_part"), jet.r() / 100.0, jet.eta(), weight);
    registry.fill(HIST("h2_jet_r_part_jet_phi_part"), jet.r() / 100.0, jet.phi(), weight);
    registry.fill(HIST("h2_jet_r_part_jet_ntracks_part"), jet.r() / 100.0, jet.tracks().size(), weight);

    for (auto& constituent : jet.template tracks_as<aod::McParticles>()) {

      registry.fill(HIST("h2_jet_pt_part_track_pt_part"), jet.pt(), constituent.pt(), weight);
      registry.fill(HIST("h2_jet_pt_part_track_eta_part"), jet.pt(), constituent.eta(), weight);
      registry.fill(HIST("h2_jet_pt_part_track_phi_part"), jet.pt(), constituent.phi(), weight);
      registry.fill(HIST("h_track_pt_part"), constituent.pt(), weight);
      registry.fill(HIST("h_track_eta_part"), constituent.eta(), weight);
      registry.fill(HIST("h_track_phi_part"), constituent.phi(), weight);
    }
  }

  template <typename T>
  void fillMCMatchedHistograms(T const& mcdjet, float weight = 1.0)
  {

    if (mcdjet.has_matchedJetGeo()) {
      const auto& mcpjet = mcdjet.template matchedJetGeo_as<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>>();

      registry.fill(HIST("h2_jet_pt_part_jet_pt"), mcpjet.pt(), mcdjet.pt(), weight);
      registry.fill(HIST("h2_jet_eta_part_jet_eta"), mcpjet.eta(), mcdjet.eta(), weight);
      registry.fill(HIST("h2_jet_phi_part_jet_phi"), mcpjet.phi(), mcdjet.phi(), weight);
      registry.fill(HIST("h2_jet_ntracks_part_jet_ntracks"), mcpjet.tracks().size(), mcdjet.tracks().size(), weight);
      registry.fill(HIST("h3_jet_pt_part_jet_eta_part_jet_eta"), mcpjet.pt(), mcpjet.eta(), mcdjet.eta(), weight);
      registry.fill(HIST("h3_jet_pt_part_jet_phi_part_jet_phi"), mcpjet.pt(), mcpjet.phi(), mcdjet.phi(), weight);
      registry.fill(HIST("h3_jet_pt_part_jet_ntracks_part_jet_ntracks"), mcpjet.pt(), mcpjet.tracks().size(), mcdjet.tracks().size(), weight);
      registry.fill(HIST("h2_jet_pt_part_jet_pt_diff"), mcpjet.pt(), (mcpjet.pt() - mcdjet.pt()) / mcpjet.pt(), weight);
    }
  }

  void processDummy(aod::Tracks const& track) {}
  PROCESS_SWITCH(JetFinderQATask, processDummy, "Dummy process function turned on by default", true);

  void processJetsData(soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>::iterator const& jet, JetTracks const& tracks)
  {
    fillHistograms(jet);
  }
  PROCESS_SWITCH(JetFinderQATask, processJetsData, "jet finder QA data", false);

  void processJetsMCD(soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>::iterator const& jet, JetTracks const& tracks)
  {
    fillHistograms(jet);
  }
  PROCESS_SWITCH(JetFinderQATask, processJetsMCD, "jet finder QA mcd", false);

  void processJetsMCDWeighted(soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetEventWeights>::iterator const& jet, JetTracks const& tracks)
  {
    fillHistograms(jet, jet.eventWeight());
  }
  PROCESS_SWITCH(JetFinderQATask, processJetsMCDWeighted, "jet finder QA mcd with weighted events", false);

  void processJetsMCP(soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>::iterator const& jet, aod::McParticles const& particles)
  {

    fillMCPHistograms(jet);
  }
  PROCESS_SWITCH(JetFinderQATask, processJetsMCP, "jet finder QA mcp", false);

  void processJetsMCPWeighted(soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetEventWeights>::iterator const& jet, aod::McParticles const& particles)
  {

    fillMCPHistograms(jet, jet.eventWeight());
  }
  PROCESS_SWITCH(JetFinderQATask, processJetsMCPWeighted, "jet finder QA mcp with weighted events", false);

  void processJetsMCPMCDMatched(aod::Collisions::iterator const& collision,
                                soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets> const& mcdjets,
                                soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets> const& mcpjets,
                                JetTracks const& tracks, aod::McParticles const& particles)
  {
    for (const auto& mcdjet : mcdjets) {
      fillMCMatchedHistograms(mcdjet);
    }
  }
  PROCESS_SWITCH(JetFinderQATask, processJetsMCPMCDMatched, "jet finder QA matched mcp and mcd", false);

  void processJetsMCPMCDMatchedWeighted(aod::Collisions::iterator const& collision,
                                        soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetEventWeights> const& mcdjets,
                                        soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCParticleLevelJetEventWeights> const& mcpjets,
                                        JetTracks const& tracks, aod::McParticles const& particles)
  {
    for (const auto& mcdjet : mcdjets) {
      fillMCMatchedHistograms(mcdjet, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetFinderQATask, processJetsMCPMCDMatchedWeighted, "jet finder QA matched mcp and mcd with weighted events", false);

  void processMCCollisionsWeighted(aod::McCollision const& collision)
  {

    registry.fill(HIST("h_collision_eventweight_part"), collision.weight());
    registry.fill(HIST("h_collision_eventweight_part"), 0.000000003);
    registry.fill(HIST("h_collision_eventweight_part"), 0.0000003);
    registry.fill(HIST("h_collision_eventweight_part"), 0.0003);
    registry.fill(HIST("h_collision_eventweight_part"), 0.03);
  }
  PROCESS_SWITCH(JetFinderQATask, processMCCollisionsWeighted, "collision QA for weighted events", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetFinderQATask>(cfgc, TaskName{"jet-finder-qa"})}; }
