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

// jet finder full+neutral QA task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/Core/HfHelper.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

#include "EventFiltering/filterTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename JetTableData, typename JetConstituentTableData, typename JetTableMCD, typename JetConstituentTableMCD, typename JetMatchingTableMCDMCP, typename JetTableMCDWeighted, typename JetTableMCP, typename JetConstituentTableMCP, typename JetMatchingTableMCPMCD, typename JetTableMCPWeighted>
struct JetFinderFullQATask {
  HistogramRegistry registry;

  Configurable<float> selectedJetsRadius{"selectedJetsRadius", 0.4, "resolution parameter for histograms without radius"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 100.0, "maximum pT acceptance for tracks"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  // cluster level configurables
  Configurable<std::string> clusterDefinitionS{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};
  Configurable<float> clusterEtaMin{"clusterEtaMin", -0.7, "minimum cluster eta"}; // For ECMAL: |eta| < 0.7, phi = 1.40 - 3.26
  Configurable<float> clusterEtaMax{"clusterEtaMax", 0.7, "maximum cluster eta"};  // For ECMAL: |eta| < 0.7, phi = 1.40 - 3.26
  Configurable<float> clusterPhiMin{"clusterPhiMin", -999, "minimum cluster phi"};
  Configurable<float> clusterPhiMax{"clusterPhiMax", 999, "maximum cluster phi"};
  Configurable<float> clusterEnergyMin{"clusterEnergyMin", 0.5, "minimum cluster energy in EMCAL (GeV)"};
  Configurable<float> clusterTimeMin{"clusterTimeMin", -999., "minimum Cluster time (ns)"};
  Configurable<float> clusterTimeMax{"clusterTimeMax", 999., "maximum Cluster time (ns)"};
  Configurable<bool> clusterRejectExotics{"clusterRejectExotics", true, "Reject exotic clusters"};

  std::vector<bool> filledJetR;
  std::vector<double> jetRadiiValues;

  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {

    trackSelection = JetDerivedDataUtilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
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
      registry.add("h_jet_nclusters", "jet N clusters;N_{jet clusters};entries", {HistType::kTH1F, {{200, -0.5, 199.5}}});
      registry.add("h3_jet_r_jet_pt_jet_eta", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_jet_phi", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_eta_jet_phi", "#it{R}_{jet};#eta_{jet};#varphi_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {80, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_jet_ntracks", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});N_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -0.5, 99.5}}});
      registry.add("h3_jet_r_jet_pt_jet_nclusters", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});N_{jet clusters}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -0.5, 99.5}}});
      registry.add("h3_jet_r_jet_pt_jet_area", "#it{R}_{jet}; #it{p}_{T,jet} (GeV/#it{c}); #it{Area}_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, 0., 1.}}});
      registry.add("h3_jet_r_jet_pt_track_pt", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_track_eta", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_track_phi", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_cluster_pt", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,cluster} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_cluster_eta", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{cluster}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_cluster_phi", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{cluster}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_cluster_energy", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});E_{cluster} (GeV)", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.0}}});
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
    }

    if (doprocessJetsMCPMCDMatched || doprocessJetsMCPMCDMatchedWeighted) {
      registry.add("h3_jet_r_jet_pt_part_jet_pt", "#it{R}_{jet};#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_eta_part_jet_eta", "#it{R}_{jet};#eta_{jet}^{part};#eta_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_phi_part_jet_phi", "#it{R}_{jet};#varphi_{jet}^{part};#varphi_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {80, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_ntracks_part_jet_ntracks", "#it{R}_{jet};N_{jet tracks}^{part};N_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -0.5, 99.5}, {100, -0.5, 99.5}}});
      registry.add("h3_jet_r_cluster_pt_part_cluster_pt", "#it{R}_{jet};#it{p}_{T,cluster}^{part} (GeV/#it{c});#it{p}_{T,cluster} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_cluster_eta_part_cluster_eta", "#it{R}_{jet};#eta_{cluster}^{part};#eta_{cluster}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_cluster_phi_part_cluster_phi", "#it{R}_{jet};#varphi_{cluster}^{part};#varphi_{cluster}", {HistType::kTH3F, {{jetRadiiBins, ""}, {80, -1.0, 7.}, {80, -1.0, 7.}}});
      registry.add("h3_jet_r_cluster_energy_part_cluster_energy", "#it{R}_{jet};#E_{cluster}^{part};#E_{cluster}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.0}}});
      registry.add("h3_jet_r_jet_pt_part_jet_pt_diff", "#it{R}_{jet};#it{p}_{T,jet}^{part} (GeV/#it{c}); (#it{p}_{T,jet}^{part} (GeV/#it{c}) - #it{p}_{T,jet} (GeV/#it{c})) / #it{p}_{T,jet}^{part} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_part_jet_eta_diff", "#it{R}_{jet};#it{p}_{T,jet}^{part} (GeV/#it{c}); (#eta_{jet}^{part} - #eta_{jet}) / #eta_{jet}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_part_jet_phi_diff", "#it{R}_{jet};#it{p}_{T,jet}^{part} (GeV/#it{c}); (#varphi_{jet}^{part} - #varphi_{jet}) / #varphi_{jet}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_pt_part_jet_eta_part_jet_eta", ";#it{p}_{T,jet}^{part} (GeV/#it{c}); #eta_{jet}^{part}; #eta_{jet}", {HistType::kTH3F, {{200, 0.0, 200}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_pt_part_jet_phi_part_jet_phi", ";#it{p}_{T,jet}^{part} (GeV/#it{c}); #varphi_{jet}^{part}; #varphi_{jet}", {HistType::kTH3F, {{200, 0.0, 200}, {80, -1.0, 7.}, {80, -1.0, 7.}}});
      registry.add("h3_jet_pt_part_jet_ntracks_part_jet_ntracks", ";#it{p}_{T,jet}^{part} (GeV/#it{c}); N_{jet tracks}^{part}; N_{jet tracks}", {HistType::kTH3F, {{200, 0.0, 200}, {100, -0.5, 99.5}, {100, -0.5, 99.5}}});
    }

    if (doprocessTracks) {
      registry.add("h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.add("h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_track_eta", "track #eta;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_track_phi", "track #varphi;#varphi_{track};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});
      registry.add("h_cluster_pt", "cluster pT;#it{p}_{T,cluster} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_cluster_eta", "cluster #eta;#eta_{cluster};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_cluster_phi", "cluster #varphi;#varphi_{cluster};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});
      registry.add("h_cluster_energy", "cluster E;E_{cluster} (GeV);entries", {HistType::kTH1F, {{200, 0., 200.}}});
    }

    if (doprocessMCCollisionsWeighted) {
      AxisSpec weightAxis = {{VARIABLE_WIDTH, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0}, "weights"};
      registry.add("h_collision_eventweight_part", "event weight;event weight;entries", {HistType::kTH1F, {weightAxis}});
    }
  }

  using JetTracks = aod::JTracks;
  using JetClusters = aod::JClusters;
  using JetParticles = aod::JMcParticles;
  using JetTableDataJoined = soa::Join<JetTableData, JetConstituentTableData>;
  using JetTableMCDJoined = soa::Join<JetTableMCD, JetConstituentTableMCD>;
  using JetTableMCDWeightedJoined = soa::Join<JetTableMCD, JetConstituentTableMCD, JetTableMCDWeighted>;
  using JetTableMCDMatchedJoined = soa::Join<JetTableMCD, JetConstituentTableMCD, JetMatchingTableMCDMCP>;
  using JetTableMCDMatchedWeightedJoined = soa::Join<JetTableMCD, JetConstituentTableMCD, JetMatchingTableMCDMCP, JetTableMCDWeighted>;
  using JetTableMCPJoined = soa::Join<JetTableMCP, JetConstituentTableMCP>;
  using JetTableMCPWeightedJoined = soa::Join<JetTableMCP, JetConstituentTableMCP, JetTableMCPWeighted>;
  using JetTableMCPMatchedJoined = soa::Join<JetTableMCP, JetConstituentTableMCP, JetMatchingTableMCPMCD>;
  using JetTableMCPMatchedWeightedJoined = soa::Join<JetTableMCD, JetConstituentTableMCP, JetMatchingTableMCPMCD, JetTableMCPWeighted>;

  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  aod::EMCALClusterDefinition clusterDefinition = aod::emcalcluster::getClusterDefinitionFromString(clusterDefinitionS.value);
  Filter clusterFilter = (aod::jcluster::definition == static_cast<int>(clusterDefinition) && aod::jcluster::eta > clusterEtaMin && aod::jcluster::eta < clusterEtaMax && aod::jcluster::phi >= clusterPhiMin && aod::jcluster::phi <= clusterPhiMax && aod::jcluster::energy >= clusterEnergyMin && aod::jcluster::time > clusterTimeMin && aod::jcluster::time < clusterTimeMax && (clusterRejectExotics && aod::jcluster::isExotic != true));

  template <typename T>
  void fillHistograms(T const& jet, float weight = 1.0)
  {
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("h_jet_ntracks"), jet.tracks().size(), weight);
      registry.fill(HIST("h_jet_nclusters"), jet.clusters().size(), weight);
    }

    registry.fill(HIST("h3_jet_r_jet_pt_jet_eta"), jet.r() / 100.0, jet.pt(), jet.eta(), weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_phi"), jet.r() / 100.0, jet.pt(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_jet_eta_jet_phi"), jet.r() / 100.0, jet.eta(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_ntracks"), jet.r() / 100.0, jet.pt(), jet.tracks().size(), weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_nclusters"), jet.r() / 100.0, jet.pt(), jet.clusters().size(), weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_area"), jet.r() / 100.0, jet.pt(), jet.area(), weight);

    for (auto& constituent : jet.template tracks_as<JetTracks>()) {

      registry.fill(HIST("h3_jet_r_jet_pt_track_pt"), jet.r() / 100.0, jet.pt(), constituent.pt(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_track_eta"), jet.r() / 100.0, jet.pt(), constituent.eta(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_track_phi"), jet.r() / 100.0, jet.pt(), constituent.phi(), weight);
    }

    for (auto& cluster : jet.template clusters_as<JetClusters>()) {
      double clusterpt = cluster.energy() / std::cosh(cluster.eta());
      registry.fill(HIST("h3_jet_r_jet_pt_cluster_pt"), jet.r() / 100.0, jet.pt(), clusterpt, weight);
      registry.fill(HIST("h3_jet_r_jet_pt_cluster_eta"), jet.r() / 100.0, jet.pt(), cluster.eta(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_cluster_phi"), jet.r() / 100.0, jet.pt(), cluster.phi(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_cluster_energy"), jet.r() / 100.0, jet.pt(), cluster.energy(), weight);
    }
  }

  template <typename T>
  void fillMCPHistograms(T const& jet, float weight = 1.0)
  {
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_jet_pt_part"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta_part"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi_part"), jet.phi(), weight);
      registry.fill(HIST("h_jet_ntracks_part"), jet.tracks().size(), weight);
    }

    registry.fill(HIST("h3_jet_r_part_jet_pt_part_jet_eta_part"), jet.r() / 100.0, jet.pt(), jet.eta(), weight);
    registry.fill(HIST("h3_jet_r_part_jet_pt_part_jet_phi_part"), jet.r() / 100.0, jet.pt(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_part_jet_eta_part_jet_phi_part"), jet.r() / 100.0, jet.eta(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_part_jet_pt_part_jet_ntracks_part"), jet.r() / 100.0, jet.pt(), jet.tracks().size(), weight);

    for (auto& constituent : jet.template tracks_as<JetParticles>()) {

      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_pt_part"), jet.r() / 100.0, jet.pt(), constituent.pt(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_eta_part"), jet.r() / 100.0, jet.pt(), constituent.eta(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_phi_part"), jet.r() / 100.0, jet.pt(), constituent.phi(), weight);
    }
  }

  template <typename T, typename U>
  void fillMCMatchedHistograms(T const& mcdjet, float weight = 1.0)
  {
    for (auto& mcpjet : mcdjet.template matchedJetGeo_as<std::decay_t<U>>()) {

      registry.fill(HIST("h3_jet_r_jet_pt_part_jet_pt"), mcdjet.r() / 100.0, mcpjet.pt(), mcdjet.pt(), weight);
      registry.fill(HIST("h3_jet_r_jet_eta_part_jet_eta"), mcdjet.r() / 100.0, mcpjet.eta(), mcdjet.eta(), weight);
      registry.fill(HIST("h3_jet_r_jet_phi_part_jet_phi"), mcdjet.r() / 100.0, mcpjet.phi(), mcdjet.phi(), weight);
      registry.fill(HIST("h3_jet_r_jet_ntracks_part_jet_ntracks"), mcdjet.r() / 100.0, mcpjet.tracks().size(), mcdjet.tracks().size(), weight); // Note: At truth level tracks include both neutral and charged particles. Is it worth differentiating these in the jetFinder?
      registry.fill(HIST("h3_jet_r_jet_pt_part_jet_pt_diff"), mcdjet.r() / 100.0, mcpjet.pt(), (mcpjet.pt() - mcdjet.pt()) / mcpjet.pt(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_part_jet_eta_diff"), mcdjet.r() / 100.0, mcpjet.pt(), (mcpjet.eta() - mcdjet.eta()) / mcpjet.eta(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_part_jet_phi_diff"), mcdjet.r() / 100.0, mcpjet.pt(), (mcpjet.phi() - mcdjet.phi()) / mcpjet.phi(), weight);

      if (mcdjet.r() == round(selectedJetsRadius * 100.0f)) {
        registry.fill(HIST("h3_jet_pt_part_jet_eta_part_jet_eta"), mcpjet.pt(), mcpjet.eta(), mcdjet.eta(), weight);
        registry.fill(HIST("h3_jet_pt_part_jet_phi_part_jet_phi"), mcpjet.pt(), mcpjet.phi(), mcdjet.phi(), weight);
        registry.fill(HIST("h3_jet_pt_part_jet_ntracks_part_jet_ntracks"), mcpjet.pt(), mcpjet.tracks().size(), mcdjet.tracks().size(), weight);
      }
    }
  }

  void processDummy(aod::JCollisions const& collision)
  {
  }
  PROCESS_SWITCH(JetFinderFullQATask, processDummy, "dummy task", true);

  void processJetsData(typename JetTableDataJoined::iterator const& jet, JetTracks const& tracks, JetClusters const& clusters)
  {
    fillHistograms(jet);
  }
  PROCESS_SWITCH(JetFinderFullQATask, processJetsData, "jet finder HF QA data", false);

  void processJetsMCD(typename JetTableMCDJoined::iterator const& jet, JetTracks const& tracks, JetClusters const& clusters)
  {
    fillHistograms(jet);
  }
  PROCESS_SWITCH(JetFinderFullQATask, processJetsMCD, "jet finder HF QA mcd", false);

  void processJetsMCDWeighted(typename JetTableMCDWeightedJoined::iterator const& jet, JetTracks const& tracks, JetClusters const& clusters)
  {
    fillHistograms(jet, jet.eventWeight());
  }
  PROCESS_SWITCH(JetFinderFullQATask, processJetsMCDWeighted, "jet finder HF QA mcd on weighted events", false);

  void processJetsMCP(typename JetTableMCPJoined::iterator const& jet, JetParticles const& particles)
  {
    fillMCPHistograms(jet);
  }
  PROCESS_SWITCH(JetFinderFullQATask, processJetsMCP, "jet finder HF QA mcp", false);

  void processJetsMCPWeighted(typename JetTableMCDWeightedJoined::iterator const& jet, JetParticles const& particles)
  {
    fillMCPHistograms(jet, jet.eventWeight());
  }
  PROCESS_SWITCH(JetFinderFullQATask, processJetsMCPWeighted, "jet finder HF QA mcp on weighted events", false);

  void processJetsMCPMCDMatched(aod::JCollisions::iterator const& collision,
                                JetTableMCDMatchedJoined const& mcdjets,
                                JetTableMCPMatchedJoined const& mcpjets,
                                JetTracks const& tracks,
                                JetClusters const& clusters,
                                JetParticles const& particles)
  {

    for (const auto& mcdjet : mcdjets) {

      fillMCMatchedHistograms<typename JetTableMCDMatchedJoined::iterator, JetTableMCPMatchedJoined>(mcdjet);
    }
  }
  PROCESS_SWITCH(JetFinderFullQATask, processJetsMCPMCDMatched, "jet finder HF QA matched mcp and mcd", false);

  void processJetsMCPMCDMatchedWeighted(aod::JCollisions::iterator const& collision,
                                        JetTableMCDMatchedWeightedJoined const& mcdjets,
                                        JetTableMCPMatchedWeightedJoined const& mcpjets,
                                        JetTracks const& tracks,
                                        JetClusters const& clusters,
                                        JetParticles const& particles)
  {

    for (const auto& mcdjet : mcdjets) {

      fillMCMatchedHistograms<typename JetTableMCDMatchedWeightedJoined::iterator, JetTableMCPMatchedWeightedJoined>(mcdjet, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetFinderFullQATask, processJetsMCPMCDMatchedWeighted, "jet finder HF QA matched mcp and mcd on weighted events", false);

  void processMCCollisionsWeighted(aod::JMcCollision const& collision)
  {

    registry.fill(HIST("h_collision_eventweight_part"), collision.weight());
  }
  PROCESS_SWITCH(JetFinderFullQATask, processMCCollisionsWeighted, "collision QA for weighted events", false);

  void processTracks(aod::JCollision const& collision,
                     soa::Filtered<JetTracks> const& tracks,
                     soa::Filtered<JetClusters> const& clusters)
  {

    registry.fill(HIST("h_collisions"), 0.5);
    if (!JetDerivedDataUtilities::eventEMCAL(collision)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    for (auto const& track : tracks) {
      if (!JetDerivedDataUtilities::selectTrack(track, trackSelection) || !JetDerivedDataUtilities::applyTrackKinematics(track, trackPtMin, trackPtMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      registry.fill(HIST("h_track_pt"), track.pt());
      registry.fill(HIST("h_track_eta"), track.eta());
      registry.fill(HIST("h_track_phi"), track.phi());
    }
    for (auto const& cluster : clusters) {
      double clusterpt = cluster.energy() / std::cosh(cluster.eta());
      registry.fill(HIST("h_cluster_pt"), clusterpt);
      registry.fill(HIST("h_cluster_eta"), cluster.eta());
      registry.fill(HIST("h_cluster_phi"), cluster.phi());
      registry.fill(HIST("h_cluster_energy"), cluster.energy());
    }
  }
  PROCESS_SWITCH(JetFinderFullQATask, processTracks, "QA for charged tracks", false);
};

using JetFinderFullJetsQATask = JetFinderFullQATask<aod::FullJets, aod::FullJetConstituents, aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents, aod::FullMCDetectorLevelJetsMatchedToFullMCParticleLevelJets, aod::FullMCDetectorLevelJetEventWeights, aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents, aod::FullMCParticleLevelJetsMatchedToFullMCDetectorLevelJets, aod::FullMCParticleLevelJetEventWeights>;
using JetFinderNeutralJetsQATask = JetFinderFullQATask<aod::NeutralJets, aod::NeutralJetConstituents, aod::NeutralMCDetectorLevelJets, aod::NeutralMCDetectorLevelJetConstituents, aod::NeutralMCDetectorLevelJetsMatchedToNeutralMCParticleLevelJets, aod::NeutralMCDetectorLevelJetEventWeights, aod::NeutralMCParticleLevelJets, aod::NeutralMCParticleLevelJetConstituents, aod::NeutralMCParticleLevelJetsMatchedToNeutralMCDetectorLevelJets, aod::NeutralMCParticleLevelJetEventWeights>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetFinderFullJetsQATask>(cfgc,
                                                                SetDefaultProcesses{},
                                                                TaskName{"jet-finder-full-qa"}));

  tasks.emplace_back(adaptAnalysisTask<JetFinderNeutralJetsQATask>(cfgc,
                                                                   SetDefaultProcesses{},
                                                                   TaskName{"jet-finder-neutral-qa"}));

  return WorkflowSpec{tasks};
}
