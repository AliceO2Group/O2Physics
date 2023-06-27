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

#include "EventFiltering/filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_2prong;

#include "Framework/runDataProcessing.h"

struct JetFinderHFQATask {

  AxisSpec weightAxis = {{VARIABLE_WIDTH, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0}, "weights"};
  HistogramRegistry registry{"registry",
                             {{"h_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}, true},
                              {"h_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}, true},
                              {"h_jet_phi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}, true},
                              {"h_jet_ntracks", "jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {{100, -0.5, 99.5}}}, true},
                              {"h2_jet_pt_jet_eta", ";#it{p}_{T,jet} (GeV/#it{c}); #eta_{jet}", {HistType::kTH2F, {{200, 0.0, 200}, {100, -1.0, 1.0}}}, true},
                              {"h2_jet_pt_jet_phi", ";#it{p}_{T,jet} (GeV/#it{c}); #phi_{jet}", {HistType::kTH2F, {{200, 0.0, 200}, {80, -1.0, 7.}}}, true},
                              {"h2_jet_pt_jet_ntracks", ";#it{p}_{T,jet} (GeV/#it{c}); N_{jet tracks}", {HistType::kTH2F, {{200, 0.0, 200}, {100, -0.5, 99.5}}}, true},
                              {"h2_jet_r_jet_pt", ";#it{R}_{jet}; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{10, 0.05, 1.05}, {200, 0.0, 200}}}, true},
                              {"h2_jet_r_jet_eta", ";#it{R}_{jet}; #eta_{jet}", {HistType::kTH2F, {{10, 0.05, 1.05}, {100, -1.0, 1.0}}}, true},
                              {"h2_jet_r_jet_phi", ";#it{R}_{jet}; #phi_{jet}", {HistType::kTH2F, {{10, 0.05, 1.05}, {80, -1.0, 7.}}}, true},
                              {"h2_jet_r_jet_ntracks", ";#it{R}_{jet}; N_{jet tracks}", {HistType::kTH2F, {{10, 0.05, 1.05}, {100, -0.5, 99.5}}}, true},
                              {"h2_jet_pt_track_pt", ";#it{p}_{T,jet} (GeV/#it{c}); #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{200, 0.0, 200}, {200, 0.0, 200.0}}}, true},
                              {"h2_jet_pt_track_eta", ";#it{p}_{T,jet} (GeV/#it{c}); #eta_{track}", {HistType::kTH2F, {{200, 0.0, 200}, {100, -1.0, 1.0}}}, true},
                              {"h2_jet_pt_track_phi", ";#it{p}_{T,jet} (GeV/#it{c}); #phi_{track}", {HistType::kTH2F, {{200, 0.0, 200}, {80, -1.0, 7.}}}, true},
                              {"h2_jet_pt_leadingtrack_pt", ";#it{p}_{T,jet} (GeV/#it{c}); #it{p}_{T,leading track} (GeV/#it{c})", {HistType::kTH2F, {{200, 0.0, 200}, {200, 0.0, 200.0}}}, true},
                              {"h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0.0, 200.0}}}, true},
                              {"h_track_eta", "track #eta;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}, true},
                              {"h_track_phi", "track #phi;#phi_{track};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}, true},
                              {"h_candidate_pt", "candidate pT;#it{p}_{T,candidate} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0.0, 200.0}}}},
                              {"h_candidate_eta", "candidate #eta;#eta_{candidate};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_candidate_phi", "candidate #phi;#phi_{candidate};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_candidate_y", "candidate y;y_{candidate};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h2_jet_pt_candidate_pt", ";#it{p}_{T,jet} (GeV/#it{c}); #it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH2F, {{200, 0.0, 200}, {200, 0.0, 200}}}},
                              {"h_jet_pt_part", "jet pT;#it{p}_{T,jet}^{part}(GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}, true},
                              {"h_jet_eta_part", "jet #eta;#eta_{jet}^{part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}, true},
                              {"h_jet_phi_part", "jet #phi;#phi_{jet}^{part};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}, true},
                              {"h_jet_ntracks_part", "jet N tracks;N_{jet tracks}^{part};entries", {HistType::kTH1F, {{100, -0.5, 99.5}}}, true},
                              {"h2_jet_pt_part_jet_eta_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c}); #eta_{jet}", {HistType::kTH2F, {{200, 0.0, 200}, {100, -1.0, 1.0}}}, true},
                              {"h2_jet_pt_part_jet_phi_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c}); #phi_{jet}", {HistType::kTH2F, {{200, 0.0, 200}, {80, -1.0, 7.}}}, true},
                              {"h2_jet_pt_part_jet_ntracks_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c}); N_{jet tracks}", {HistType::kTH2F, {{200, 0.0, 200}, {100, -0.5, 99.5}}}, true},
                              {"h2_jet_r_part_jet_pt_part", ";#it{R}_{jet}^{part}; #it{p}_{T,jet}^{part} (GeV/#it{c})", {HistType::kTH2F, {{10, 0.05, 1.05}, {200, 0.0, 200}}}, true},
                              {"h2_jet_r_part_jet_eta_part", ";#it{R}_{jet}^{part}; #eta_{jet}^{part}", {HistType::kTH2F, {{10, 0.05, 1.05}, {100, -1.0, 1.0}}}, true},
                              {"h2_jet_r_part_jet_phi_part", ";#it{R}_{jet}^{part}; #phi_{jet}^{part}", {HistType::kTH2F, {{10, 0.05, 1.05}, {80, -1.0, 7.}}}, true},
                              {"h2_jet_r_part_jet_ntracks_part", ";#it{R}_{jet}^{part}; N_{jet tracks}^{part}", {HistType::kTH2F, {{10, 0.05, 1.05}, {100, -0.5, 99.5}}}, true},
                              {"h2_jet_pt_part_track_pt_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c}); #it{p}_{T,track}^{part} (GeV/#it{c})", {HistType::kTH2F, {{200, 0.0, 200}, {200, 0.0, 200.0}}}, true},
                              {"h2_jet_pt_part_track_eta_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c}); #eta_{track}^{part}", {HistType::kTH2F, {{200, 0.0, 200}, {100, -1.0, 1.0}}}, true},
                              {"h2_jet_pt_part_track_phi_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c}); #phi_{track}^{part}", {HistType::kTH2F, {{200, 0.0, 200}, {80, -1.0, 7.}}}, true},
                              {"h_track_pt_part", "track pT;#it{p}_{T,track}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0.0, 200.0}}}, true},
                              {"h_track_eta_part", "track #eta;#eta_{track}^{part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}, true},
                              {"h_track_phi_part", "track #phi;#phi_{track}^{part};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}, true},
                              {"h_candidate_pt_part", "candidate pT;#it{p}_{T,candidate}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0.0, 200.0}}}},
                              {"h_candidate_eta_part", "candidate #eta;#eta_{candidate}^{part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_candidate_phi_part", "candidate #phi;#phi_{candidate}^{part};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_candidate_y_part", "candidate y;y_{candidate}^{part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h2_jet_pt_part_candidate_pt_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c}); #it{p}_{T,candidate}^{part} (GeV/#it{c})", {HistType::kTH2F, {{200, 0.0, 200}, {200, 0.0, 200}}}},
                              {"h2_jet_pt_part_jet_pt", ";#it{p}_{T,jet}^{part} (GeV/#it{c}); #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{200, 0.0, 200}, {200, 0.0, 200}}}, true},
                              {"h2_jet_eta_part_jet_eta", ";#eta_{jet}^{part}; #eta_{jet}", {HistType::kTH2F, {{100, -1.0, 1.0}, {100, -1.0, 1.0}}}, true},
                              {"h2_jet_phi_part_jet_phi", ";#phi_{jet}^{part}; #phi_{jet}", {HistType::kTH2F, {{80, -1.0, 7.}, {80, -1.0, 7.}}}, true},
                              {"h2_jet_ntracks_part_jet_ntracks", ";N_{jet tracks}^{part}; N_{jet tracks}", {HistType::kTH2F, {{100, -0.5, 99.5}, {100, -0.5, 99.5}}}, true},
                              {"h2_candidate_pt_part_candidate_pt", ";#it{p}_{T,candidate}^{part} (GeV/#it{c}); #it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH2F, {{200, 0.0, 200}, {200, 0.0, 200}}}},
                              {"h2_candidate_eta_part_candidate_eta", ";#eta_{candidate}^{part}; #eta_{candidate}", {HistType::kTH2F, {{100, -1.0, 1.0}, {100, -1.0, 1.0}}}},
                              {"h2_candidate_phi_part_candidate_phi", ";#phi_{candidate}^{part}; #phi_{candidate}", {HistType::kTH2F, {{80, -1.0, 7.}, {80, -1.0, 7.}}}},
                              {"h2_candidate_y_part_candidate_y", ";#y_{candidate}^{part}; #y_{candidate}", {HistType::kTH2F, {{100, -1.0, 1.0}, {100, -1.0, 1.0}}}},
                              {"h3_jet_pt_part_jet_eta_part_jet_eta", ";#it{p}_{T,jet}^{part} (GeV/#it{c}); #eta_{jet}^{part}; #eta_{jet}", {HistType::kTH3F, {{200, 0.0, 200}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}}, true},
                              {"h3_jet_pt_part_jet_phi_part_jet_phi", ";#it{p}_{T,jet}^{part} (GeV/#it{c}); #phi_{jet}^{part}; #phi_{jet}", {HistType::kTH3F, {{200, 0.0, 200}, {80, -1.0, 7.}, {80, -1.0, 7.}}}, true},
                              {"h3_jet_pt_part_jet_ntracks_part_jet_ntracks", ";#it{p}_{T,jet}^{part} (GeV/#it{c}); N_{jet tracks}^{part}; N_{jet tracks}", {HistType::kTH3F, {{200, 0.0, 200}, {100, -0.5, 99.5}, {100, -0.5, 99.5}}}, true},
                              {"h2_jet_pt_part_jet_pt_diff", ";#it{p}_{T,jet}^{part} (GeV/#it{c}); (#it{p}_{T,jet}^{part} (GeV/#it{c}) - #it{p}_{T,jet} (GeV/#it{c})) / #it{p}_{T,jet}^{part} (GeV/#it{c})", {HistType::kTH2F, {{200, 0.0, 200}, {1000, -5.0, 5.0}}}, true},
                              {"h_collision_trigger_events", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}}, true},
                              {"h_track_pt_MB", "track pT for MB events;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0.0, 200.0}}}, true},
                              {"h_track_eta_MB", "track #eta for MB events;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}, true},
                              {"h_track_phi_MB", "track #phi for MB events;#phi_{track};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}, true},
                              {"h_track_pt_Triggered", "track pT for Triggered events;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0.0, 200.0}}}, true},
                              {"h_track_eta_Triggered", "track #eta for Triggered events;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}, true},
                              {"h_track_phi_Triggered", "track #phi for Triggered events;#phi_{track};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}, true},
                              {"h_collision_eventweight_part", "event weight;event weight;entries", {HistType::kTH1F, {weightAxis}}}}};

  Configurable<float> triggeredJetsRadius{"triggeredJetsRadius", 0.6, "resolution parameter for triggered jets"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.2, 0.3, 0.4, 0.5, 0.6}, "jet resolution parameters"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 100.0, "maximum pT acceptance for tracks"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  std::string trackSelection;
  std::vector<double> minJetPt;
  std::vector<double> jetRadiiValues;

  void init(o2::framework::InitContext&)
  {
    trackSelection = static_cast<std::string>(trackSelections);
    jetRadiiValues = (std::vector<double>)jetRadii;

    for (auto iJetRadius = 0; iJetRadius < jetRadiiValues.size(); iJetRadius++) {
      minJetPt.push_back(0.0);
    }
    auto jetRadiiBins = (std::vector<double>)jetRadii;
    if (jetRadiiBins.size() > 1) {
      jetRadiiBins.push_back(jetRadiiBins[jetRadiiBins.size() - 1] + (TMath::Abs(jetRadiiBins[jetRadiiBins.size() - 1] - jetRadiiBins[jetRadiiBins.size() - 2])));
    } else {
      jetRadiiBins.push_back(jetRadiiBins[jetRadiiBins.size() - 1] + 0.1);
    }
    registry.add("h3_jet_radius_jet_pt_jet_eta", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#phi_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
    registry.add("h3_jet_radius_jet_pt_jet_phi", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#phi_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
    registry.add("h3_jet_radius_jet_eta_jet_phi", "#it{R}_{jet};#eta_{jet};#phi_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {80, -1.0, 7.}}});
    registry.add("h3_jet_radius_part_jet_pt_part_jet_eta_part", ";#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#phi_{jet}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
    registry.add("h3_jet_radius_part_jet_pt_part_jet_phi_part", ";#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#phi_{jet}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
    registry.add("h3_jet_radius_part_jet_eta_part_jet_phi_part", ";#it{R}_{jet}^{part};#eta_{jet}^{part};#phi_{jet}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {80, -1.0, 7.}}});
    registry.add("h3_jet_radius_jet_pt_collision", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});collision trigger status", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {2, -0.5, 1.5}}});
    registry.add("h3_jet_radius_jet_eta_collision", "#it{R}_{jet};#eta_{jet};collision trigger status", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {2, -0.5, 1.5}}});
    registry.add("h3_jet_radius_jet_phi_collision", "#it{R}_{jet};#phi_{jet};collision trigger status", {HistType::kTH3F, {{jetRadiiBins, ""}, {80, -1.0, 7.}, {2, -0.5, 1.5}}});
    registry.add("h2_jet_radius_jet_pT_triggered", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{jetRadiiBins, ""}, {200, 0., 200.}}});
    registry.add("h3_jet_radius_jet_pt_track_pt_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
    registry.add("h3_jet_radius_jet_pt_track_eta_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
    registry.add("h3_jet_radius_jet_pt_track_phi_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
    registry.add("h3_jet_radius_jet_pt_track_pt_Triggered", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
    registry.add("h3_jet_radius_jet_pt_track_eta_Triggered", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
    registry.add("h3_jet_radius_jet_pt_track_phi_Triggered", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
    registry.add("h3_jet_radius_jet_pt_candidate_pt_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
    registry.add("h3_jet_radius_jet_pt_candidate_eta_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
    registry.add("h3_jet_radius_jet_pt_candidate_phi_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
    registry.add("h3_jet_radius_jet_pt_candidate_y_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});y_{candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
    registry.add("h3_jet_radius_jet_pt_candidate_pt_Triggered", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
    registry.add("h3_jet_radius_jet_pt_candidate_eta_Triggered", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
    registry.add("h3_jet_radius_jet_pt_candidate_phi_Triggered", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
    registry.add("h3_jet_radius_jet_pt_candidate_y_Triggered", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});y_{candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
  }

  using JetTracks = soa::Join<aod::Tracks, aod::TrackSelection>;
  using CandidateD0Data = soa::Join<aod::HfCand2Prong, aod::HfSelD0>;
  using CandidateD0MC = soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>;
  using JetParticles2Prong = soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>;
  using JetParticles3Prong = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;
  using JetParticlesBplus = soa::Join<aod::McParticles, aod::HfCandBplusMcGen>;

  template <typename T>
  bool selectTrack(T const& track)
  {
    if (trackSelection == "globalTracks" && !track.isGlobalTrackWoPtEta()) {
      return false;
    }
    if (trackSelection == "globalTracks" && (track.pt() < trackPtMin || track.pt() >= trackPtMax || track.eta() < trackEtaMin || track.eta() >= trackEtaMax)) {
      return false;
    }
    if (trackSelection == "QualityTracks" && !track.isQualityTrack()) {
      return false;
    }
    if (trackSelection == "hybridTracksJE" && !track.trackCutFlagFb5()) { // isQualityTrack
      return false;
    }
    return true;
  }

  template <typename T>
  void fillDataHistograms(T const& jet, float weight = 1.0)
  {

    registry.fill(HIST("h_jet_pt"), jet.pt(), weight);
    registry.fill(HIST("h_jet_eta"), jet.eta(), weight);
    registry.fill(HIST("h_jet_phi"), jet.phi(), weight);
    registry.fill(HIST("h_jet_ntracks"), jet.tracks().size() + jet.hfcandidates().size(), weight);

    registry.fill(HIST("h2_jet_pt_jet_eta"), jet.pt(), jet.eta(), weight);
    registry.fill(HIST("h2_jet_pt_jet_phi"), jet.pt(), jet.phi(), weight);
    registry.fill(HIST("h2_jet_pt_jet_ntracks"), jet.pt(), jet.tracks().size() + jet.hfcandidates().size(), weight);

    registry.fill(HIST("h2_jet_r_jet_pt"), jet.r() / 100.0, jet.pt(), weight);
    registry.fill(HIST("h2_jet_r_jet_eta"), jet.r() / 100.0, jet.eta(), weight);
    registry.fill(HIST("h2_jet_r_jet_phi"), jet.r() / 100.0, jet.phi(), weight);
    registry.fill(HIST("h2_jet_r_jet_ntracks"), jet.r() / 100.0, jet.tracks().size() + jet.hfcandidates().size(), weight);

    registry.fill(HIST("h3_jet_radius_jet_pt_jet_eta"), jet.r() / 100.0, jet.pt(), jet.eta(), weight);
    registry.fill(HIST("h3_jet_radius_jet_pt_jet_phi"), jet.r() / 100.0, jet.pt(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_radius_jet_eta_jet_phi"), jet.r() / 100.0, jet.eta(), jet.phi(), weight);

    for (auto& constituent : jet.template tracks_as<JetTracks>()) {

      registry.fill(HIST("h2_jet_pt_track_pt"), jet.pt(), constituent.pt(), weight);
      registry.fill(HIST("h2_jet_pt_track_eta"), jet.pt(), constituent.eta(), weight);
      registry.fill(HIST("h2_jet_pt_track_phi"), jet.pt(), constituent.phi(), weight);
      registry.fill(HIST("h_track_pt"), constituent.pt(), weight);
      registry.fill(HIST("h_track_eta"), constituent.eta(), weight);
      registry.fill(HIST("h_track_phi"), constituent.phi(), weight);
    }

    for (auto& hfcandidate : jet.template hfcandidates_as<CandidateD0Data>()) {

      registry.fill(HIST("h2_jet_pt_track_pt"), jet.pt(), hfcandidate.pt(), weight);
      registry.fill(HIST("h2_jet_pt_track_eta"), jet.pt(), hfcandidate.eta(), weight);
      registry.fill(HIST("h2_jet_pt_track_phi"), jet.pt(), hfcandidate.phi(), weight);
      registry.fill(HIST("h_track_pt"), hfcandidate.pt(), weight);
      registry.fill(HIST("h_track_eta"), hfcandidate.eta(), weight);
      registry.fill(HIST("h_track_phi"), hfcandidate.phi(), weight);

      registry.fill(HIST("h_candidate_pt"), hfcandidate.pt(), weight);
      registry.fill(HIST("h_candidate_eta"), hfcandidate.eta(), weight);
      registry.fill(HIST("h_candidate_phi"), hfcandidate.phi(), weight);
      registry.fill(HIST("h_candidate_y"), hfcandidate.y(RecoDecay::getMassPDG(421)), weight);
      registry.fill(HIST("h2_jet_pt_candidate_pt"), jet.pt(), hfcandidate.pt(), weight);
    }
  }

  template <typename T>
  void fillMCDHistograms(T const& jet, float weight = 1.0)
  {

    registry.fill(HIST("h_jet_pt"), jet.pt(), weight);
    registry.fill(HIST("h_jet_eta"), jet.eta(), weight);
    registry.fill(HIST("h_jet_phi"), jet.phi(), weight);
    registry.fill(HIST("h_jet_ntracks"), jet.tracks().size() + jet.hfcandidates().size(), weight);

    registry.fill(HIST("h2_jet_pt_jet_eta"), jet.pt(), jet.eta(), weight);
    registry.fill(HIST("h2_jet_pt_jet_phi"), jet.pt(), jet.phi(), weight);
    registry.fill(HIST("h2_jet_pt_jet_ntracks"), jet.pt(), jet.tracks().size() + jet.hfcandidates().size(), weight);

    registry.fill(HIST("h2_jet_r_jet_pt"), jet.r() / 100.0, jet.pt(), weight);
    registry.fill(HIST("h2_jet_r_jet_eta"), jet.r() / 100.0, jet.eta(), weight);
    registry.fill(HIST("h2_jet_r_jet_phi"), jet.r() / 100.0, jet.phi(), weight);
    registry.fill(HIST("h2_jet_r_jet_ntracks"), jet.r() / 100.0, jet.tracks().size() + jet.hfcandidates().size(), weight);

    registry.fill(HIST("h3_jet_radius_jet_pt_jet_eta"), jet.r() / 100.0, jet.pt(), jet.eta(), weight);
    registry.fill(HIST("h3_jet_radius_jet_pt_jet_phi"), jet.r() / 100.0, jet.pt(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_radius_jet_eta_jet_phi"), jet.r() / 100.0, jet.eta(), jet.phi(), weight);

    for (auto& constituent : jet.template tracks_as<JetTracks>()) {

      registry.fill(HIST("h2_jet_pt_track_pt"), jet.pt(), constituent.pt(), weight);
      registry.fill(HIST("h2_jet_pt_track_eta"), jet.pt(), constituent.eta(), weight);
      registry.fill(HIST("h2_jet_pt_track_phi"), jet.pt(), constituent.phi(), weight);
      registry.fill(HIST("h_track_pt"), constituent.pt(), weight);
      registry.fill(HIST("h_track_eta"), constituent.eta(), weight);
      registry.fill(HIST("h_track_phi"), constituent.phi(), weight);
    }

    for (auto& hfcandidate : jet.template hfcandidates_as<CandidateD0MC>()) {

      registry.fill(HIST("h2_jet_pt_track_pt"), jet.pt(), hfcandidate.pt(), weight);
      registry.fill(HIST("h2_jet_pt_track_eta"), jet.pt(), hfcandidate.eta(), weight);
      registry.fill(HIST("h2_jet_pt_track_phi"), jet.pt(), hfcandidate.phi(), weight);
      registry.fill(HIST("h_track_pt"), hfcandidate.pt(), weight);
      registry.fill(HIST("h_track_eta"), hfcandidate.eta(), weight);
      registry.fill(HIST("h_track_phi"), hfcandidate.phi(), weight);

      registry.fill(HIST("h_candidate_pt"), hfcandidate.pt(), weight);
      registry.fill(HIST("h_candidate_eta"), hfcandidate.eta(), weight);
      registry.fill(HIST("h_candidate_phi"), hfcandidate.phi(), weight);
      registry.fill(HIST("h_candidate_y"), hfcandidate.y(RecoDecay::getMassPDG(421)), weight);
      registry.fill(HIST("h2_jet_pt_candidate_pt"), jet.pt(), hfcandidate.pt(), weight);
    }
  }

  template <typename T>
  void fillMCPHistograms(T const& jet, float weight = 1.0)
  {

    registry.fill(HIST("h_jet_pt_part"), jet.pt(), weight);
    registry.fill(HIST("h_jet_eta_part"), jet.eta(), weight);
    registry.fill(HIST("h_jet_phi_part"), jet.phi(), weight);
    registry.fill(HIST("h_jet_ntracks_part"), jet.tracks().size() + jet.hfcandidates().size(), weight);

    registry.fill(HIST("h2_jet_pt_part_jet_eta_part"), jet.pt(), jet.eta(), weight);
    registry.fill(HIST("h2_jet_pt_part_jet_phi_part"), jet.pt(), jet.phi(), weight);
    registry.fill(HIST("h2_jet_pt_part_jet_ntracks_part"), jet.pt(), jet.tracks().size() + jet.hfcandidates().size(), weight);

    registry.fill(HIST("h2_jet_r_part_jet_pt_part"), jet.r() / 100.0, jet.pt(), weight);
    registry.fill(HIST("h2_jet_r_part_jet_eta_part"), jet.r() / 100.0, jet.eta(), weight);
    registry.fill(HIST("h2_jet_r_part_jet_phi_part"), jet.r() / 100.0, jet.phi(), weight);
    registry.fill(HIST("h2_jet_r_part_jet_ntracks_part"), jet.r() / 100.0, jet.tracks().size() + jet.hfcandidates().size(), weight);

    registry.fill(HIST("h3_jet_radius_part_jet_pt_part_jet_eta_part"), jet.r() / 100.0, jet.pt(), jet.eta(), weight);
    registry.fill(HIST("h3_jet_radius_part_jet_pt_part_jet_phi_part"), jet.r() / 100.0, jet.pt(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_radius_part_jet_eta_part_jet_phi_part"), jet.r() / 100.0, jet.eta(), jet.phi(), weight);

    for (auto& constituent : jet.template tracks_as<JetParticles2Prong>()) {

      registry.fill(HIST("h2_jet_pt_part_track_pt_part"), jet.pt(), constituent.pt(), weight);
      registry.fill(HIST("h2_jet_pt_part_track_eta_part"), jet.pt(), constituent.eta(), weight);
      registry.fill(HIST("h2_jet_pt_part_track_phi_part"), jet.pt(), constituent.phi(), weight);
      registry.fill(HIST("h_track_pt_part"), constituent.pt(), weight);
      registry.fill(HIST("h_track_eta_part"), constituent.eta(), weight);
      registry.fill(HIST("h_track_phi_part"), constituent.phi(), weight);
    }

    for (auto& hfcandidate : jet.template hfcandidates_as<JetParticles2Prong>()) {

      registry.fill(HIST("h2_jet_pt_part_track_pt_part"), jet.pt(), hfcandidate.pt(), weight);
      registry.fill(HIST("h2_jet_pt_part_track_eta_part"), jet.pt(), hfcandidate.eta(), weight);
      registry.fill(HIST("h2_jet_pt_part_track_phi_part"), jet.pt(), hfcandidate.phi(), weight);
      registry.fill(HIST("h_track_pt_part"), hfcandidate.pt(), weight);
      registry.fill(HIST("h_track_eta_part"), hfcandidate.eta(), weight);
      registry.fill(HIST("h_track_phi_part"), hfcandidate.phi(), weight);

      registry.fill(HIST("h_candidate_pt_part"), hfcandidate.pt(), weight);
      registry.fill(HIST("h_candidate_eta_part"), hfcandidate.eta(), weight);
      registry.fill(HIST("h_candidate_phi_part"), hfcandidate.phi(), weight);
      registry.fill(HIST("h_candidate_y_part"), hfcandidate.y(), weight);
      registry.fill(HIST("h2_jet_pt_part_candidate_pt_part"), jet.pt(), hfcandidate.pt(), weight);
    }
  }

  template <typename T>
  void fillMCMatchedHistograms(T const& mcdjet, float weight = 1.0)
  {
    if (mcdjet.has_matchedJetCand() && mcdjet.matchedJetCandId() >= 0) {
      const auto& mcpjet = mcdjet.template matchedJetCand_as<soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents>>();

      auto mcdCandPt = 0.0;
      auto mcdCandPhi = 0.0;
      auto mcdCandEta = 0.0;
      auto mcdCandY = 0.0;
      auto mcpCandPt = 0.0;
      auto mcpCandPhi = 0.0;
      auto mcpCandEta = 0.0;
      auto mcpCandY = 0.0;
      for (auto& hfcandidate_mcd : mcdjet.template hfcandidates_as<CandidateD0MC>()) {

        mcdCandPt = hfcandidate_mcd.pt();
        mcdCandPhi = hfcandidate_mcd.phi();
        mcdCandEta = hfcandidate_mcd.eta();
        mcdCandY = hfcandidate_mcd.y(RecoDecay::getMassPDG(421));
      }

      for (auto& hfcandidate_mcp : mcpjet.template hfcandidates_as<JetParticles2Prong>()) {

        mcpCandPt = hfcandidate_mcp.pt();
        mcpCandPhi = hfcandidate_mcp.phi();
        mcpCandEta = hfcandidate_mcp.eta();
        mcpCandY = hfcandidate_mcp.y();
      }

      registry.fill(HIST("h2_jet_pt_part_jet_pt"), mcpjet.pt(), mcdjet.pt(), weight);
      registry.fill(HIST("h2_jet_eta_part_jet_eta"), mcpjet.eta(), mcdjet.eta(), weight);
      registry.fill(HIST("h2_jet_phi_part_jet_phi"), mcpjet.phi(), mcdjet.phi(), weight);
      registry.fill(HIST("h2_jet_ntracks_part_jet_ntracks"), mcpjet.tracks().size() + mcpjet.hfcandidates().size(), mcdjet.tracks().size() + mcdjet.hfcandidates().size(), weight);
      registry.fill(HIST("h2_candidate_pt_part_candidate_pt"), mcpCandPt, mcdCandPt, weight);
      registry.fill(HIST("h2_candidate_eta_part_candidate_eta"), mcpCandEta, mcdCandEta, weight);
      registry.fill(HIST("h2_candidate_phi_part_candidate_phi"), mcpCandPhi, mcdCandPhi, weight);
      registry.fill(HIST("h2_candidate_y_part_candidate_y"), mcpCandY, mcdCandY, weight);
      registry.fill(HIST("h3_jet_pt_part_jet_eta_part_jet_eta"), mcpjet.pt(), mcpjet.eta(), mcdjet.eta(), weight);
      registry.fill(HIST("h3_jet_pt_part_jet_phi_part_jet_phi"), mcpjet.pt(), mcpjet.phi(), mcdjet.phi(), weight);
      registry.fill(HIST("h3_jet_pt_part_jet_ntracks_part_jet_ntracks"), mcpjet.pt(), mcpjet.tracks().size() + mcpjet.hfcandidates().size(), mcdjet.tracks().size() + mcdjet.hfcandidates().size(), weight);
      registry.fill(HIST("h2_jet_pt_part_jet_pt_diff"), mcpjet.pt(), (mcpjet.pt() - mcdjet.pt()) / mcpjet.pt(), weight);
    }
  }

  void processDummy(aod::Tracks const& track) {}
  PROCESS_SWITCH(JetFinderHFQATask, processDummy, "Dummy process function turned on by default", true);

  void processJetsData(soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents>::iterator const& jet, CandidateD0Data const& candidates, JetTracks const& tracks)
  {

    fillDataHistograms(jet);
    auto leadingTrackpT = -1.0;
    for (auto const& track : tracks) {
      if (!selectTrack(track)) {
        continue;
      }
      if (track.pt() > leadingTrackpT)
        leadingTrackpT = track.pt();
    }
    registry.fill(HIST("h2_jet_pt_leadingtrack_pt"), jet.pt(), leadingTrackpT);
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsData, "jet finder HF QA data", false);

  void processJetsMCD(soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents>::iterator const& jet, CandidateD0MC const& candidates, JetTracks const& tracks)
  {
    fillMCDHistograms(jet);
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsMCD, "jet finder HF QA mcd", false);

  void processJetsMCDWeighted(soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetEventWeights>::iterator const& jet, CandidateD0MC const& candidates, JetTracks const& tracks)
  {
    fillMCDHistograms(jet, jet.eventWeight());
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsMCDWeighted, "jet finder HF QA mcd on weighted events", false);

  void processJetsMCP(soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents>::iterator const& jet, JetParticles2Prong const& particles)
  {
    fillMCPHistograms(jet);
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsMCP, "jet finder HF QA mcp", false);

  void processJetsMCPWeighted(soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelJetEventWeights>::iterator const& jet, JetParticles2Prong const& particles)
  {
    fillMCPHistograms(jet, jet.eventWeight());
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsMCPWeighted, "jet finder HF QA mcp on weighted events", false);

  void processJetsMCPMCDMatched(aod::Collisions::iterator const& collision,
                                soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets> const& mcdjets,
                                soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets> const& mcpjets,
                                CandidateD0MC const& candidates,
                                JetTracks const& tracks, JetParticles2Prong const& particles)
  {

    for (const auto& mcdjet : mcdjets) {

      fillMCMatchedHistograms(mcdjet);
    }
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsMCPMCDMatched, "jet finder HF QA matched mcp and mcd", false);

  void processJetsMCPMCDMatchedWeighted(aod::Collisions::iterator const& collision,
                                        soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets, aod::D0ChargedMCDetectorLevelJetEventWeights> const& mcdjets,
                                        soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets, aod::D0ChargedMCParticleLevelJetEventWeights> const& mcpjets,
                                        CandidateD0MC const& candidates,
                                        JetTracks const& tracks, JetParticles2Prong const& particles)
  {

    for (const auto& mcdjet : mcdjets) {

      fillMCMatchedHistograms(mcdjet);
    }
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsMCPMCDMatchedWeighted, "jet finder HF QA matched mcp and mcd on weighted events", false);

  void processMCCollisionsWeighted(aod::McCollision const& collision)
  {

    registry.fill(HIST("h_collision_eventweight_part"), collision.weight());
  }
  PROCESS_SWITCH(JetFinderHFQATask, processMCCollisionsWeighted, "collision QA for weighted events", false);

  void processTriggeredData(soa::Join<aod::Collisions, aod::EvSels, aod::JetFilters>::iterator const& collision,
                            soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents> const& jets,
                            CandidateD0Data const& candidates,
                            JetTracks const& tracks)
  {
    registry.fill(HIST("h_collision_trigger_events"), 0.5); // all events
    if (collision.posZ() > vertexZCut)
      return;
    registry.fill(HIST("h_collision_trigger_events"), 1.5); // all events with z vertex cut
    if (!collision.sel8())
      return;
    registry.fill(HIST("h_collision_trigger_events"), 2.5); // events with sel8()
    if (collision.hasJetChHighPt() >= 1)
      registry.fill(HIST("h_collision_trigger_events"), 3.5); // events with triggered jets

    for (auto iJetRadius = 0; iJetRadius < jetRadiiValues.size(); iJetRadius++) {
      minJetPt[iJetRadius] = 0.0;
    }
    for (auto& jet : jets) {
      for (auto iJetRadius = 0; iJetRadius < jetRadiiValues.size(); iJetRadius++) {
        if (jet.r() == round(jetRadiiValues[iJetRadius] * 100.0f)) {
          for (double pt = jet.pt(); pt > minJetPt[iJetRadius]; pt -= 1.0) {
            registry.fill(HIST("h2_jet_radius_jet_pT_triggered"), jet.r() / 100.0, pt);
          }
          if (jet.pt() > minJetPt[iJetRadius]) {
            minJetPt[iJetRadius] = jet.pt();
          }
          break;
        }
      }

      if ((jet.eta() < (trackEtaMin + jet.r() / 100.0)) || (jet.eta() > (trackEtaMax - jet.r() / 100.0))) {
        continue;
      }
      registry.fill(HIST("h3_jet_radius_jet_pt_collision"), jet.r() / 100.0, jet.pt(), collision.hasJetChHighPt());
      registry.fill(HIST("h3_jet_radius_jet_eta_collision"), jet.r() / 100.0, jet.eta(), collision.hasJetChHighPt());
      registry.fill(HIST("h3_jet_radius_jet_phi_collision"), jet.r() / 100.0, jet.phi(), collision.hasJetChHighPt());

      for (auto& constituent : jet.template tracks_as<JetTracks>()) {
        registry.fill(HIST("h3_jet_radius_jet_pt_track_pt_MB"), jet.r() / 100.0, jet.pt(), constituent.pt());
        registry.fill(HIST("h3_jet_radius_jet_pt_track_eta_MB"), jet.r() / 100.0, jet.pt(), constituent.eta());
        registry.fill(HIST("h3_jet_radius_jet_pt_track_phi_MB"), jet.r() / 100.0, jet.pt(), constituent.phi());
        if (collision.hasJetChHighPt() >= 1) {
          registry.fill(HIST("h3_jet_radius_jet_pt_track_pt_Triggered"), jet.r() / 100.0, jet.pt(), constituent.pt());
          registry.fill(HIST("h3_jet_radius_jet_pt_track_eta_Triggered"), jet.r() / 100.0, jet.pt(), constituent.eta());
          registry.fill(HIST("h3_jet_radius_jet_pt_track_phi_Triggered"), jet.r() / 100.0, jet.pt(), constituent.phi());
        }
      }
      for (auto& hfcandidate : jet.template hfcandidates_as<CandidateD0Data>()) {

        registry.fill(HIST("h3_jet_radius_jet_pt_candidate_pt_MB"), jet.r() / 100.0, jet.pt(), hfcandidate.pt());
        registry.fill(HIST("h3_jet_radius_jet_pt_candidate_eta_MB"), jet.r() / 100.0, jet.pt(), hfcandidate.eta());
        registry.fill(HIST("h3_jet_radius_jet_pt_candidate_phi_MB"), jet.r() / 100.0, jet.pt(), hfcandidate.phi());
        registry.fill(HIST("h3_jet_radius_jet_pt_candidate_y_MB"), jet.r() / 100.0, jet.pt(), hfcandidate.y(RecoDecay::getMassPDG(421)));
        if (collision.hasJetChHighPt() >= 1) {
          registry.fill(HIST("h3_jet_radius_jet_pt_candidate_pt_Triggered"), jet.r() / 100.0, jet.pt(), hfcandidate.pt());
          registry.fill(HIST("h3_jet_radius_jet_pt_candidate_eta_Triggered"), jet.r() / 100.0, jet.pt(), hfcandidate.eta());
          registry.fill(HIST("h3_jet_radius_jet_pt_candidate_phi_Triggered"), jet.r() / 100.0, jet.pt(), hfcandidate.phi());
          registry.fill(HIST("h3_jet_radius_jet_pt_candidate_y_Triggered"), jet.r() / 100.0, jet.pt(), hfcandidate.y(RecoDecay::getMassPDG(421)));
        }
      }
    }

    for (auto& track : tracks) {
      if (!selectTrack(track)) {
        continue;
      }
      registry.fill(HIST("h_track_pt_MB"), track.pt());
      registry.fill(HIST("h_track_eta_MB"), track.eta());
      registry.fill(HIST("h_track_phi_MB"), track.phi());
      if (collision.hasJetChHighPt() >= 1) {
        registry.fill(HIST("h_track_pt_Triggered"), track.pt());
        registry.fill(HIST("h_track_eta_Triggered"), track.eta());
        registry.fill(HIST("h_track_phi_Triggered"), track.phi());
      }
    }
  }

  PROCESS_SWITCH(JetFinderHFQATask, processTriggeredData, "QA for charged jet trigger", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetFinderHFQATask>(cfgc, TaskName{"jet-finder-hf-qa"})}; }
