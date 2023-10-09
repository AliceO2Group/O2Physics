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

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"

#include "EventFiltering/filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

template <typename JetTableData, typename JetConstituentTableData, typename CandidateTableData, typename JetTableMCD, typename JetConstituentTableMCD, typename JetMatchingTableMCDMCP, typename JetTableMCDWeighted, typename CandidateTableMCD, typename JetTableMCP, typename JetConstituentTableMCP, typename JetMatchingTableMCPMCD, typename JetTableMCPWeighted, typename ParticleTableMCP>
struct JetFinderHFQATask {
  HistogramRegistry registry;

  Configurable<float> selectedJetsRadius{"selectedJetsRadius", 0.4, "resolution parameter for histograms without radius"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 100.0, "maximum pT acceptance for tracks"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<int> selectionFlagLcToPKPi{"selectionFlagLcToPKPi", 1, "Selection Flag for Lc->PKPi"};
  Configurable<int> selectionFlagLcToPiPK{"selectionFlagLcToPiPK", 1, "Selection Flag for Lc->PiPK"};
  Configurable<int> selectionFlagBplus{"selectionFlagBplus", 1, "Selection Flag for B+"};

  HfHelper hfHelper;
  std::string trackSelection;
  std::vector<bool> filledJetR;
  std::vector<double> jetRadiiValues;

  void init(o2::framework::InitContext&)
  {
    trackSelection = static_cast<std::string>(trackSelections);
    jetRadiiValues = (std::vector<double>)jetRadii;

    for (std::size_t iJetRadius = 0; iJetRadius < jetRadiiValues.size(); iJetRadius++) {
      filledJetR.push_back(0.0);
    }
    auto jetRadiiBins = (std::vector<double>)jetRadii;
    if (jetRadiiBins.size() > 1) {
      jetRadiiBins.push_back(jetRadiiBins[jetRadiiBins.size() - 1] + (TMath::Abs(jetRadiiBins[jetRadiiBins.size() - 1] - jetRadiiBins[jetRadiiBins.size() - 2])));
    } else {
      jetRadiiBins.push_back(jetRadiiBins[jetRadiiBins.size() - 1] + 0.1);
    }

    if (doprocessJetsData || doprocessJetsMCD || doprocessJetsMCDWeighted) {
      registry.add("h_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_jet_phi", "jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});
      registry.add("h_jet_ntracks", "jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {{200, -0.5, 199.5}}});
      registry.add("h3_jet_r_jet_pt_jet_eta", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_jet_phi", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_eta_jet_phi", "#it{R}_{jet};#eta_{jet};#varphi_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {80, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_jet_ntracks", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});N_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -0.5, 99.5}}});
      registry.add("h3_jet_r_jet_pt_track_pt", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_track_eta", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_track_phi", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_candidate_pt", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_candidate_eta", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_candidate_phi", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_candidate_y", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});y_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_leadingtrack_pt", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c}); #it{p}_{T,leading track} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200.0}, {200, 0.0, 200.0}}});
      registry.add("h3_candidate_invmass_jet_pt_candidate_pt", ";#it{m}_{inv, candidate} (GeV/#it{c}^{2}); #it{p}_{T,jet} (GeV/#it{c}) ;#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{500, 0.0, 5.0}, {200, 0.0, 200.0}, {200, 0.0, 200.0}}});
      registry.add("h3_candidatebar_invmass_jet_pt_candidate_pt", ";#it{m}_{inv, candidate bar} (GeV/#it{c}^{2}); #it{p}_{T,jet} (GeV/#it{c}) ;#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{500, 0.0, 5.0}, {200, 0.0, 200.0}, {200, 0.0, 200.0}}});
    }

    if (doprocessJetsMCP || doprocessJetsMCPWeighted) {
      registry.add("h_jet_pt_part", "jet pT;#it{p}_{T,jet}^{part}(GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_jet_eta_part", "jet #eta;#eta_{jet}^{part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_jet_phi_part", "jet #varphi;#varphi_{jet}^{part};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});
      registry.add("h_jet_ntracks_part", "jet N tracks;N_{jet tracks}^{part};entries", {HistType::kTH1F, {{100, -0.5, 99.5}}});
      registry.add("h3_jet_r_part_jet_pt_part_jet_eta_part", ";#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#eta_{jet}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_part_jet_pt_part_jet_phi_part", ";#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#varphi_{jet}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
      registry.add("h3_jet_r_part_jet_eta_part_jet_phi_part", ";#it{R}_{jet}^{part};#eta_{jet}^{part};#varphi_{jet}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {80, -1.0, 7.}}});
      registry.add("h3_jet_r_part_jet_pt_part_jet_ntracks_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});N_{jet tracks}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -0.5, 99.5}}});
      registry.add("h3_jet_r_part_jet_pt_part_track_pt_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T,jet tracks}^{part} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_part_jet_pt_part_track_eta_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#eta_{jet tracks}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_part_jet_pt_part_track_phi_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#varphi_{jet tracks}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
      registry.add("h3_jet_r_part_jet_pt_part_candidate_pt_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T,candidate}^{part} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_part_jet_pt_part_candidate_eta_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#eta_{candidate}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_part_jet_pt_part_candidate_phi_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#varphi{candidate}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
      registry.add("h3_jet_r_part_jet_pt_part_candidate_y_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});y_{candidate}^{part} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
    }

    if (doprocessJetsMCPMCDMatched || doprocessJetsMCPMCDMatchedWeighted) {
      registry.add("h3_jet_r_jet_pt_part_jet_pt", "#it{R}_{jet};#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_eta_part_jet_eta", "#it{R}_{jet};#eta_{jet}^{part};#eta_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_phi_part_jet_phi", "#it{R}_{jet};#varphi_{jet}^{part};#varphi_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_ntracks_part_jet_ntracks", "#it{R}_{jet};N_{jet tracks}^{part};N_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -0.5, 99.5}, {100, -0.5, 99.5}}});
      registry.add("h3_jet_r_candidate_pt_part_candidate_pt", "#it{R}_{jet};#it{p}_{T,candidate}^{part} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_candidate_eta_part_candidate_eta", "#it{R}_{jet};#eta_{candidate}^{part};#eta_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_candidate_phi_part_candidate_phi", "#it{R}_{jet};#varphi_{candidate}^{part};#varphi_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
      registry.add("h3_jet_r_candidate_y_part_candidate_y", "#it{R}_{jet};#y_{candidate}^{part};#y_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_part_jet_pt_diff", "#it{R}_{jet};#it{p}_{T,jet}^{part} (GeV/#it{c}); (#it{p}_{T,jet}^{part} (GeV/#it{c}) - #it{p}_{T,jet} (GeV/#it{c})) / #it{p}_{T,jet}^{part} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_part_jet_eta_diff", "#it{R}_{jet};#it{p}_{T,jet}^{part} (GeV/#it{c}); (#eta_{jet}^{part} - #eta_{jet}) / #eta_{jet}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_part_jet_phi_diff", "#it{R}_{jet};#it{p}_{T,jet}^{part} (GeV/#it{c}); (#varphi_{jet}^{part} - #varphi_{jet}) / #varphi_{jet}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_pt_part_jet_eta_part_jet_eta", ";#it{p}_{T,jet}^{part} (GeV/#it{c}); #eta_{jet}^{part}; #eta_{jet}", {HistType::kTH3F, {{200, 0.0, 200}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_pt_part_jet_phi_part_jet_phi", ";#it{p}_{T,jet}^{part} (GeV/#it{c}); #varphi_{jet}^{part}; #varphi_{jet}", {HistType::kTH3F, {{200, 0.0, 200}, {80, -1.0, 7.}, {80, -1.0, 7.}}});
      registry.add("h3_jet_pt_part_jet_ntracks_part_jet_ntracks", ";#it{p}_{T,jet}^{part} (GeV/#it{c}); N_{jet tracks}^{part}; N_{jet tracks}", {HistType::kTH3F, {{200, 0.0, 200}, {100, -0.5, 99.5}, {100, -0.5, 99.5}}});
    }

    if (doprocessTriggeredData) {
      registry.add("h_collision_trigger_events", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.add("h_track_pt_MB", "track pT for MB events;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0.0, 200.0}}});
      registry.add("h_track_eta_MB", "track #eta for MB events;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_track_phi_MB", "track #varphi for MB events;#varphi_{track};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});
      registry.add("h_track_pt_Triggered", "track pT for Triggered events;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0.0, 200.0}}});
      registry.add("h_track_eta_Triggered", "track #eta for Triggered events;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_track_phi_Triggered", "track #varphi for Triggered events;#varphi_{track};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_collision", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});collision trigger status", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {2, -0.5, 1.5}}});
      registry.add("h3_jet_r_jet_eta_collision", "#it{R}_{jet};#eta_{jet};collision trigger status", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {2, -0.5, 1.5}}});
      registry.add("h3_jet_r_jet_phi_collision", "#it{R}_{jet};#varphi_{jet};collision trigger status", {HistType::kTH3F, {{jetRadiiBins, ""}, {80, -1.0, 7.}, {2, -0.5, 1.5}}});
      registry.add("h2_jet_r_jet_pT_triggered", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{jetRadiiBins, ""}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_track_pt_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_track_eta_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_track_phi_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_track_pt_Triggered", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_track_eta_Triggered", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_track_phi_Triggered", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_candidate_pt_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_candidate_eta_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_candidate_phi_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_candidate_y_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});y_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_candidate_pt_Triggered", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_candidate_eta_Triggered", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_candidate_phi_Triggered", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_candidate_y_Triggered", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});y_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
    }

    if (doprocessTracks) {
      registry.add("h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.add("h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_track_eta", "track #eta;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_track_phi", "track #varphi;#varphi_{track};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});
    }

    if (doprocessMCCollisionsWeighted) {
      AxisSpec weightAxis = {{VARIABLE_WIDTH, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0}, "weights"};
      registry.add("h_collision_eventweight_part", "event weight;event weight;entries", {HistType::kTH1F, {weightAxis}});
    }
  }

  using JetTracks = soa::Join<aod::Tracks, aod::TrackSelection>;
  using JetTableDataJoined = soa::Join<JetTableData, JetConstituentTableData>;
  using JetTableMCDJoined = soa::Join<JetTableMCD, JetConstituentTableMCD>;
  using JetTableMCDWeightedJoined = soa::Join<JetTableMCD, JetConstituentTableMCD, JetTableMCDWeighted>;
  using JetTableMCDMatchedJoined = soa::Join<JetTableMCD, JetConstituentTableMCD, JetMatchingTableMCDMCP>;
  using JetTableMCDMatchedWeightedJoined = soa::Join<JetTableMCD, JetConstituentTableMCD, JetMatchingTableMCDMCP, JetTableMCDWeighted>;
  using JetTableMCPJoined = soa::Join<JetTableMCP, JetConstituentTableMCP>;
  using JetTableMCPWeightedJoined = soa::Join<JetTableMCP, JetConstituentTableMCP, JetTableMCPWeighted>;
  using JetTableMCPMatchedJoined = soa::Join<JetTableMCP, JetConstituentTableMCP, JetMatchingTableMCPMCD>;
  using JetTableMCPMatchedWeightedJoined = soa::Join<JetTableMCD, JetConstituentTableMCP, JetMatchingTableMCPMCD, JetTableMCPWeighted>;

  template <typename T>
  bool selectTrack(T const& track)
  {
    if (trackSelection == "globalTracks") {
      if (track.pt() >= trackPtMin && track.pt() < trackPtMax && track.eta() >= trackEtaMin && track.eta() <= trackEtaMax) {
        return track.isGlobalTrackWoPtEta();
      } else {
        return false;
      }
    }
    if (trackSelection == "QualityTracks") {
      return track.isQualityTrack();
    }
    if (trackSelection == "hybridTracksJE") { // isQualityTrack
      return track.trackCutFlagFb5();
    }
    return true;
  }

  template <typename T, typename U>
  void fillHistograms(T const& jet, float weight = 1.0)
  {
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("h_jet_ntracks"), jet.tracks().size() + jet.hfcandidates().size(), weight);
    }

    registry.fill(HIST("h3_jet_r_jet_pt_jet_eta"), jet.r() / 100.0, jet.pt(), jet.eta(), weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_phi"), jet.r() / 100.0, jet.pt(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_jet_eta_jet_phi"), jet.r() / 100.0, jet.eta(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_ntracks"), jet.r() / 100.0, jet.pt(), jet.tracks().size() + jet.hfcandidates().size(), weight);

    for (auto& constituent : jet.template tracks_as<JetTracks>()) {

      registry.fill(HIST("h3_jet_r_jet_pt_track_pt"), jet.r() / 100.0, jet.pt(), constituent.pt(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_track_eta"), jet.r() / 100.0, jet.pt(), constituent.eta(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_track_phi"), jet.r() / 100.0, jet.pt(), constituent.phi(), weight);
    }

    for (auto& hfcandidate : jet.template hfcandidates_as<std::decay_t<U>>()) {

      registry.fill(HIST("h3_jet_r_jet_pt_track_pt"), jet.r() / 100.0, jet.pt(), hfcandidate.pt(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_track_eta"), jet.r() / 100.0, jet.pt(), hfcandidate.eta(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_track_phi"), jet.r() / 100.0, jet.pt(), hfcandidate.phi(), weight);

      registry.fill(HIST("h3_jet_r_jet_pt_candidate_pt"), jet.r() / 100.0, jet.pt(), hfcandidate.pt(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_candidate_eta"), jet.r() / 100.0, jet.pt(), hfcandidate.eta(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_candidate_phi"), jet.r() / 100.0, jet.pt(), hfcandidate.phi(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_candidate_y"), jet.r() / 100.0, jet.pt(), hfcandidate.y(RecoDecay::getMassPDG(421)), weight);

      if (jet.r() == round(selectedJetsRadius * 100.0f)) {
        if constexpr (std::is_same_v<std::decay_t<U>, soa::Join<aod::HfCand2Prong, aod::HfSelD0>> || std::is_same_v<std::decay_t<U>, soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>>) {
          if (hfcandidate.isSelD0() >= selectionFlagD0) {
            registry.fill(HIST("h3_candidate_invmass_jet_pt_candidate_pt"), hfHelper.invMassD0ToPiK(hfcandidate), jet.pt(), hfcandidate.pt(), weight);
          }
          if (hfcandidate.isSelD0bar() >= selectionFlagD0bar) {
            registry.fill(HIST("h3_candidatebar_invmass_jet_pt_candidate_pt"), hfHelper.invMassD0barToKPi(hfcandidate), jet.pt(), hfcandidate.pt(), weight);
          }
        }

        if constexpr (std::is_same_v<std::decay_t<U>, soa::Join<aod::HfCand3Prong, aod::HfSelLc>> || std::is_same_v<std::decay_t<U>, soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>>) {
          if (hfcandidate.isSelLcToPKPi() >= selectionFlagLcToPKPi) {
            registry.fill(HIST("h3_candidate_invmass_jet_pt_candidate_pt"), hfHelper.invMassLcToPKPi(hfcandidate), jet.pt(), hfcandidate.pt(), weight);
          }
          if (hfcandidate.isSelLcToPiKP() >= selectionFlagLcToPiPK) {
            registry.fill(HIST("h3_candidatebar_invmass_jet_pt_candidate_pt"), hfHelper.invMassLcToPiKP(hfcandidate), jet.pt(), hfcandidate.pt(), weight);
          }
        }

        if constexpr (std::is_same_v<std::decay_t<U>, soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi>> || std::is_same_v<std::decay_t<U>, soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi, aod::HfCandBplusMcRec>>) {
          if (hfcandidate.isSelBplusToD0Pi() >= selectionFlagBplus) {
            registry.fill(HIST("h3_candidate_invmass_jet_pt_candidate_pt"), hfHelper.invMassBplusToD0Pi(hfcandidate), jet.pt(), hfcandidate.pt(), weight);
            registry.fill(HIST("h3_candidatebar_invmass_jet_pt_candidate_pt"), hfHelper.invMassBplusToD0Pi(hfcandidate), jet.pt(), hfcandidate.pt(), weight);
          }
        }
      }
    }
  }

  template <typename T, typename U>
  void fillMCPHistograms(T const& jet, float weight = 1.0)
  {
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_jet_pt_part"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta_part"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi_part"), jet.phi(), weight);
      registry.fill(HIST("h_jet_ntracks_part"), jet.tracks().size() + jet.hfcandidates().size(), weight);
    }

    registry.fill(HIST("h3_jet_r_part_jet_pt_part_jet_eta_part"), jet.r() / 100.0, jet.pt(), jet.eta(), weight);
    registry.fill(HIST("h3_jet_r_part_jet_pt_part_jet_phi_part"), jet.r() / 100.0, jet.pt(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_part_jet_eta_part_jet_phi_part"), jet.r() / 100.0, jet.eta(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_part_jet_pt_part_jet_ntracks_part"), jet.r() / 100.0, jet.pt(), jet.tracks().size() + jet.hfcandidates().size(), weight);

    for (auto& constituent : jet.template tracks_as<std::decay_t<U>>()) {

      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_pt_part"), jet.r() / 100.0, jet.pt(), constituent.pt(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_eta_part"), jet.r() / 100.0, jet.pt(), constituent.eta(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_phi_part"), jet.r() / 100.0, jet.pt(), constituent.phi(), weight);
    }

    for (auto& hfcandidate : jet.template hfcandidates_as<std::decay_t<U>>()) {

      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_pt_part"), jet.r() / 100.0, jet.pt(), hfcandidate.pt(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_eta_part"), jet.r() / 100.0, jet.pt(), hfcandidate.eta(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_phi_part"), jet.r() / 100.0, jet.pt(), hfcandidate.phi(), weight);

      registry.fill(HIST("h3_jet_r_part_jet_pt_part_candidate_pt_part"), jet.r() / 100.0, jet.pt(), hfcandidate.pt(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_candidate_eta_part"), jet.r() / 100.0, jet.pt(), hfcandidate.eta(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_candidate_phi_part"), jet.r() / 100.0, jet.pt(), hfcandidate.phi(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_candidate_y_part"), jet.r() / 100.0, jet.pt(), hfcandidate.y(), weight);
    }
  }

  template <typename T, typename U, typename M, typename O>
  void fillMCMatchedHistograms(T const& mcdjet, float weight = 1.0)
  {
    // if (mcdjet.has_matchedJetCand() && mcdjet.matchedJetCandId() >= 0) {old version
    //   const auto& mcpjet = mcdjet.template matchedJetCand_as<std::decay_t<U>>();old version
    for (auto& mcpjet : mcdjet.template matchedJetCand_as<std::decay_t<U>>()) {
      auto mcdCandPt = 0.0;
      auto mcdCandPhi = 0.0;
      auto mcdCandEta = 0.0;
      auto mcdCandY = 0.0;
      auto mcpCandPt = 0.0;
      auto mcpCandPhi = 0.0;
      auto mcpCandEta = 0.0;
      auto mcpCandY = 0.0;
      for (auto& hfcandidate_mcd : mcdjet.template hfcandidates_as<std::decay_t<M>>()) {

        mcdCandPt = hfcandidate_mcd.pt();
        mcdCandPhi = hfcandidate_mcd.phi();
        mcdCandEta = hfcandidate_mcd.eta();
        mcdCandY = hfcandidate_mcd.y(RecoDecay::getMassPDG(421));
      }

      for (auto& hfcandidate_mcp : mcpjet.template hfcandidates_as<std::decay_t<O>>()) {

        mcpCandPt = hfcandidate_mcp.pt();
        mcpCandPhi = hfcandidate_mcp.phi();
        mcpCandEta = hfcandidate_mcp.eta();
        mcpCandY = hfcandidate_mcp.y();
      }

      registry.fill(HIST("h3_jet_r_jet_pt_part_jet_pt"), mcdjet.r() / 100.0, mcpjet.pt(), mcdjet.pt(), weight);
      registry.fill(HIST("h3_jet_r_jet_eta_part_jet_eta"), mcdjet.r() / 100.0, mcpjet.eta(), mcdjet.eta(), weight);
      registry.fill(HIST("h3_jet_r_jet_phi_part_jet_phi"), mcdjet.r() / 100.0, mcpjet.phi(), mcdjet.phi(), weight);
      registry.fill(HIST("h3_jet_r_jet_ntracks_part_jet_ntracks"), mcdjet.r() / 100.0, mcpjet.tracks().size() + mcpjet.hfcandidates().size(), mcdjet.tracks().size() + mcdjet.hfcandidates().size(), weight);
      registry.fill(HIST("h3_jet_r_candidate_pt_part_candidate_pt"), mcdjet.r() / 100.0, mcpCandPt, mcdCandPt, weight);
      registry.fill(HIST("h3_jet_r_candidate_eta_part_candidate_eta"), mcdjet.r() / 100.0, mcpCandEta, mcdCandEta, weight);
      registry.fill(HIST("h3_jet_r_candidate_phi_part_candidate_phi"), mcdjet.r() / 100.0, mcpCandPhi, mcdCandPhi, weight);
      registry.fill(HIST("h3_jet_r_candidate_y_part_candidate_y"), mcdjet.r() / 100.0, mcpCandY, mcdCandY, weight);
      registry.fill(HIST("h3_jet_r_jet_pt_part_jet_pt_diff"), mcdjet.r() / 100.0, mcpjet.pt(), (mcpjet.pt() - mcdjet.pt()) / mcpjet.pt(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_part_jet_eta_diff"), mcdjet.r() / 100.0, mcpjet.pt(), (mcpjet.eta() - mcdjet.eta()) / mcpjet.eta(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_part_jet_phi_diff"), mcdjet.r() / 100.0, mcpjet.pt(), (mcpjet.phi() - mcdjet.phi()) / mcpjet.phi(), weight);

      if (mcdjet.r() == round(selectedJetsRadius * 100.0f)) {
        registry.fill(HIST("h3_jet_pt_part_jet_eta_part_jet_eta"), mcpjet.pt(), mcpjet.eta(), mcdjet.eta(), weight);
        registry.fill(HIST("h3_jet_pt_part_jet_phi_part_jet_phi"), mcpjet.pt(), mcpjet.phi(), mcdjet.phi(), weight);
        registry.fill(HIST("h3_jet_pt_part_jet_ntracks_part_jet_ntracks"), mcpjet.pt(), mcpjet.tracks().size() + mcpjet.hfcandidates().size(), mcdjet.tracks().size() + mcdjet.hfcandidates().size(), weight);
      }
    }
  }

  void processDummy(aod::Collisions const& collision)
  {
  }
  PROCESS_SWITCH(JetFinderHFQATask, processDummy, "dummy task", true);

  void processJetsData(typename JetTableDataJoined::iterator const& jet, CandidateTableData const& candidates, JetTracks const& tracks)
  {
    fillHistograms<typename JetTableDataJoined::iterator, CandidateTableData>(jet);
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsData, "jet finder HF QA data", false);

  void processJetsMCD(typename JetTableMCDJoined::iterator const& jet, CandidateTableMCD const& candidates, JetTracks const& tracks)
  {
    fillHistograms<typename JetTableMCDJoined::iterator, CandidateTableMCD>(jet);
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsMCD, "jet finder HF QA mcd", false);

  void processJetsMCDWeighted(typename JetTableMCDWeightedJoined::iterator const& jet, CandidateTableMCD const& candidates, JetTracks const& tracks)
  {
    fillHistograms<typename JetTableMCDWeightedJoined::iterator, CandidateTableMCD>(jet, jet.eventWeight());
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsMCDWeighted, "jet finder HF QA mcd on weighted events", false);

  void processJetsMCP(typename JetTableMCPJoined::iterator const& jet, ParticleTableMCP const& particles)
  {
    fillMCPHistograms<typename JetTableMCPJoined::iterator, ParticleTableMCP>(jet);
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsMCP, "jet finder HF QA mcp", false);

  void processJetsMCPWeighted(typename JetTableMCDWeightedJoined::iterator const& jet, ParticleTableMCP const& particles)
  {
    fillMCPHistograms<typename JetTableMCDWeightedJoined::iterator, ParticleTableMCP>(jet, jet.eventWeight());
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsMCPWeighted, "jet finder HF QA mcp on weighted events", false);

  void processJetsMCPMCDMatched(aod::Collisions::iterator const& collision,
                                JetTableMCDMatchedJoined const& mcdjets,
                                JetTableMCPMatchedJoined const& mcpjets,
                                CandidateTableMCD const& candidates,
                                JetTracks const& tracks, ParticleTableMCP const& particles)
  {

    for (const auto& mcdjet : mcdjets) {

      fillMCMatchedHistograms<typename JetTableMCDMatchedJoined::iterator, JetTableMCPMatchedJoined, CandidateTableMCD, ParticleTableMCP>(mcdjet);
    }
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsMCPMCDMatched, "jet finder HF QA matched mcp and mcd", false);

  void processJetsMCPMCDMatchedWeighted(aod::Collisions::iterator const& collision,
                                        JetTableMCDMatchedWeightedJoined const& mcdjets,
                                        JetTableMCPMatchedWeightedJoined const& mcpjets,
                                        CandidateTableMCD const& candidates,
                                        JetTracks const& tracks, ParticleTableMCP const& particles)
  {

    for (const auto& mcdjet : mcdjets) {

      fillMCMatchedHistograms<typename JetTableMCDMatchedWeightedJoined::iterator, JetTableMCPMatchedWeightedJoined, CandidateTableMCD, ParticleTableMCP>(mcdjet, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsMCPMCDMatchedWeighted, "jet finder HF QA matched mcp and mcd on weighted events", false);

  void processMCCollisionsWeighted(aod::McCollision const& collision)
  {

    registry.fill(HIST("h_collision_eventweight_part"), collision.weight());
  }
  PROCESS_SWITCH(JetFinderHFQATask, processMCCollisionsWeighted, "collision QA for weighted events", false);

  void processTriggeredData(soa::Join<aod::Collisions, aod::EvSels, aod::JetFilters>::iterator const& collision,
                            JetTableDataJoined const& jets,
                            CandidateTableData const& candidates,
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

    for (std::size_t iJetRadius = 0; iJetRadius < jetRadiiValues.size(); iJetRadius++) {
      filledJetR[iJetRadius] = false;
    }
    for (auto& jet : jets) {
      for (std::size_t iJetRadius = 0; iJetRadius < jetRadiiValues.size(); iJetRadius++) {
        if (jet.r() == round(jetRadiiValues[iJetRadius] * 100.0f) && !filledJetR[iJetRadius]) {
          filledJetR[iJetRadius] = true;
          for (double pt = 0.0; pt <= jet.pt(); pt += 1.0) {
            registry.fill(HIST("h2_jet_r_jet_pT_triggered"), jet.r() / 100.0, pt);
          }
          break;
        }
      }

      if ((jet.eta() < (trackEtaMin + jet.r() / 100.0)) || (jet.eta() > (trackEtaMax - jet.r() / 100.0))) {
        continue;
      }
      registry.fill(HIST("h3_jet_r_jet_pt_collision"), jet.r() / 100.0, jet.pt(), collision.hasJetChHighPt());
      registry.fill(HIST("h3_jet_r_jet_eta_collision"), jet.r() / 100.0, jet.eta(), collision.hasJetChHighPt());
      registry.fill(HIST("h3_jet_r_jet_phi_collision"), jet.r() / 100.0, jet.phi(), collision.hasJetChHighPt());

      for (auto& constituent : jet.template tracks_as<JetTracks>()) {
        registry.fill(HIST("h3_jet_r_jet_pt_track_pt_MB"), jet.r() / 100.0, jet.pt(), constituent.pt());
        registry.fill(HIST("h3_jet_r_jet_pt_track_eta_MB"), jet.r() / 100.0, jet.pt(), constituent.eta());
        registry.fill(HIST("h3_jet_r_jet_pt_track_phi_MB"), jet.r() / 100.0, jet.pt(), constituent.phi());
        if (collision.hasJetChHighPt() >= 1) {
          registry.fill(HIST("h3_jet_r_jet_pt_track_pt_Triggered"), jet.r() / 100.0, jet.pt(), constituent.pt());
          registry.fill(HIST("h3_jet_r_jet_pt_track_eta_Triggered"), jet.r() / 100.0, jet.pt(), constituent.eta());
          registry.fill(HIST("h3_jet_r_jet_pt_track_phi_Triggered"), jet.r() / 100.0, jet.pt(), constituent.phi());
        }
      }
      for (auto& hfcandidate : jet.template hfcandidates_as<CandidateTableData>()) {

        registry.fill(HIST("h3_jet_r_jet_pt_candidate_pt_MB"), jet.r() / 100.0, jet.pt(), hfcandidate.pt());
        registry.fill(HIST("h3_jet_r_jet_pt_candidate_eta_MB"), jet.r() / 100.0, jet.pt(), hfcandidate.eta());
        registry.fill(HIST("h3_jet_r_jet_pt_candidate_phi_MB"), jet.r() / 100.0, jet.pt(), hfcandidate.phi());
        registry.fill(HIST("h3_jet_r_jet_pt_candidate_y_MB"), jet.r() / 100.0, jet.pt(), hfcandidate.y(RecoDecay::getMassPDG(421)));
        if (collision.hasJetChHighPt() >= 1) {
          registry.fill(HIST("h3_jet_r_jet_pt_candidate_pt_Triggered"), jet.r() / 100.0, jet.pt(), hfcandidate.pt());
          registry.fill(HIST("h3_jet_r_jet_pt_candidate_eta_Triggered"), jet.r() / 100.0, jet.pt(), hfcandidate.eta());
          registry.fill(HIST("h3_jet_r_jet_pt_candidate_phi_Triggered"), jet.r() / 100.0, jet.pt(), hfcandidate.phi());
          registry.fill(HIST("h3_jet_r_jet_pt_candidate_y_Triggered"), jet.r() / 100.0, jet.pt(), hfcandidate.y(RecoDecay::getMassPDG(421)));
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

  void processTracks(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                     JetTracks const& tracks)
  {

    registry.fill(HIST("h_collisions"), 0.5);
    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    for (auto const& track : tracks) {
      if (!selectTrack(track)) {
        continue;
      }
      registry.fill(HIST("h_track_pt"), track.pt());
      registry.fill(HIST("h_track_eta"), track.eta());
      registry.fill(HIST("h_track_phi"), track.phi());
    }
  }
  PROCESS_SWITCH(JetFinderHFQATask, processTracks, "QA for charged tracks", false);
};

using JetFinderD0QATask = JetFinderHFQATask<aod::D0ChargedJets, aod::D0ChargedJetConstituents, soa::Join<aod::HfCand2Prong, aod::HfSelD0>, aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets, aod::D0ChargedMCDetectorLevelJetEventWeights, soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>, aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets, aod::D0ChargedMCParticleLevelJetEventWeights, soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>>;
using JetFinderLcQATask = JetFinderHFQATask<aod::LcChargedJets, aod::LcChargedJetConstituents, soa::Join<aod::HfCand3Prong, aod::HfSelLc>, aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents, aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCParticleLevelJets, aod::LcChargedMCDetectorLevelJetEventWeights, soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>, aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents, aod::LcChargedMCParticleLevelJetsMatchedToLcChargedMCDetectorLevelJets, aod::LcChargedMCParticleLevelJetEventWeights, soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>;
using JetFinderBplusQATask = JetFinderHFQATask<aod::BplusChargedJets, aod::BplusChargedJetConstituents, soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi>, aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents, aod::BplusChargedMCDetectorLevelJetsMatchedToBplusChargedMCParticleLevelJets, aod::BplusChargedMCDetectorLevelJetEventWeights, soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi, aod::HfCandBplusMcRec>, aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents, aod::BplusChargedMCParticleLevelJetsMatchedToBplusChargedMCDetectorLevelJets, aod::BplusChargedMCParticleLevelJetEventWeights, soa::Join<aod::McParticles, aod::HfCandBplusMcGen>>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetFinderD0QATask>(cfgc,
                                                          SetDefaultProcesses{},
                                                          TaskName{"jet-finder-charged-d0-qa"}));

  tasks.emplace_back(adaptAnalysisTask<JetFinderBplusQATask>(cfgc,
                                                             SetDefaultProcesses{},
                                                             TaskName{"jet-finder-charged-bplus-qa"}));

  tasks.emplace_back(adaptAnalysisTask<JetFinderLcQATask>(cfgc,
                                                          SetDefaultProcesses{},
                                                          TaskName{"jet-finder-charged-lc-qa"}));

  return WorkflowSpec{tasks};
}
