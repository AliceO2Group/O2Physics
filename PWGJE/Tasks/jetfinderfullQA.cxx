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
#include "PWGJE/Core/JetFindingUtilities.h"

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
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};
  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 100.0, "maximum pT acceptance for tracks"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};

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

  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};

  std::vector<bool> filledJetR;
  std::vector<double> jetRadiiValues;

  int trackSelection = -1;

  Service<o2::framework::O2DatabasePDG> pdgDatabase;

  void init(o2::framework::InitContext&)
  {

    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
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
      registry.add("h_jet_phi", "jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, -1.0, 7.}}});
      registry.add("h_jet_ntracks", "jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {{200, -0.5, 199.5}}});
      registry.add("h_jet_nclusters", "jet N clusters;N_{jet clusters};entries", {HistType::kTH1F, {{200, -0.5, 199.5}}});
      registry.add("h2_centrality_jet_pt", "centrality vs #it{p}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{1200, -10.0, 110.0}, {200, 0., 200.}}});
      registry.add("h2_centrality_jet_eta", "centrality vs #eta_{jet}; centrality; #eta_{jet}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {100, -1.0, 1.0}}});
      registry.add("h2_centrality_jet_phi", "centrality vs #varphi_{jet}; centrality; #varphi_{jet}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {160, -1.0, 7.}}});
      registry.add("h2_centrality_jet_ntracks", "centrality vs N_{jet tracks}; centrality; N_{jet tracks}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {200, -0.5, 199.5}}});
      registry.add("h2_centrality_jet_nclusters", "centrality vs N_{jet clusters}; centrality; N_{jet clusters}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_jet_pt_centrality", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});centrality", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {1200, -10.0, 110.0}}});
      registry.add("h3_jet_r_jet_pt_jet_eta", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_jet_phi", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_eta_jet_phi", "#it{R}_{jet};#eta_{jet};#varphi_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_jet_ntracks", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});N_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_jet_pt_jet_nclusters", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});N_{jet clusters}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_jet_pt_jet_neutralenergyfraction", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});E_{neutral}/E_{total}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {120, 0.0, 1.2}}});
      registry.add("h3_jet_r_jet_pt_jet_area", "#it{R}_{jet}; #it{p}_{T,jet} (GeV/#it{c}); #it{area}_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {300, 0., 3.}}});
      registry.add("h3_jet_r_jet_pt_track_pt", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_track_eta", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_track_phi", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_cluster_pt", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,cluster} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_cluster_eta", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{cluster}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_cluster_phi", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{cluster}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_cluster_energy", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});E_{cluster} (GeV)", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.0}}});
    }

    if (doprocessJetsMCP || doprocessJetsMCPWeighted) {
      registry.add("h_jet_pt_part", "jet pT;#it{p}_{T,jet}^{part}(GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_jet_eta_part", "jet #eta;#eta_{jet}^{part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_jet_phi_part", "jet #varphi;#varphi_{jet}^{part};entries", {HistType::kTH1F, {{160, -1.0, 7.}}});
      registry.add("h_jet_ntracks_part", "jet N tracks;N_{jet tracks}^{part};entries", {HistType::kTH1F, {{200, -0.5, 199.5}}});
      registry.add("h3_jet_r_part_jet_pt_part_jet_eta_part", ";#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#eta_{jet}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_part_jet_pt_part_jet_phi_part", ";#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#varphi_{jet}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_part_jet_eta_part_jet_phi_part", ";#it{R}_{jet}^{part};#eta_{jet}^{part};#varphi_{jet}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_part_jet_pt_part_jet_ntracks_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});N_{jet tracks}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_part_jet_pt_part_jet_neutralenergyfraction_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});E_{neutral}^{part}/E_{total}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {120, 0.0, 1.2}}});
      registry.add("h3_jet_r_part_jet_pt_part_track_pt_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T,jet tracks}^{part} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_part_jet_pt_part_track_eta_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#eta_{jet tracks}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_part_jet_pt_part_track_phi_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#varphi_{jet tracks}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
    }

    if (doprocessJetsMCPMCDMatched || doprocessJetsMCPMCDMatchedWeighted) {
      registry.add("h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c});#it{p}_{T,jet}^{base} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_eta_tag_jet_eta_base_matchedgeo", "#it{R}_{jet};#eta_{jet}^{tag};#eta_{jet}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_phi_tag_jet_phi_base_matchedgeo", "#it{R}_{jet};#varphi_{jet}^{tag};#varphi_{jet}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedgeo", "#it{R}_{jet};N_{jet tracks}^{tag};N_{jet tracks}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, -0.5, 199.5}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_cluster_pt_tag_cluster_pt_base_matchedgeo", "#it{R}_{jet};#it{p}_{T,cluster}^{tag} (GeV/#it{c});#it{p}_{T,cluster}^{base} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_cluster_eta_tag_cluster_eta_base_matchedgeo", "#it{R}_{jet};#eta_{cluster}^{tag};#eta_{cluster}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_cluster_phi_tag_cluster_phi_base_matchedgeo", "#it{R}_{jet};#varphi_{cluster}^{tag};#varphi_{cluster}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {160, -1.0, 7.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_cluster_energy_tag_cluster_energy_base_matchedgeo", "#it{R}_{jet};#E_{cluster}^{tag};#E_{cluster}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedgeo", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_eta_base_diff_matchedgeo", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#eta_{jet}^{tag} - #eta_{jet}^{base}) / #eta_{jet}^{tag}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_phi_base_diff_matchedgeo", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#varphi_{jet}^{tag} - #varphi_{jet}^{base}) / #varphi_{jet}^{tag}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_pt_tag_jet_eta_tag_jet_eta_base_matchedgeo", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); #eta_{jet}^{tag}; #eta_{jet}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_pt_tag_jet_phi_tag_jet_phi_base_matchedgeo", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); #varphi_{jet}^{tag}; #varphi_{jet}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {160, -1.0, 7.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_pt_tag_jet_ntracks_tag_jet_ntracks_base_matchedgeo", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); N_{jet tracks}^{tag}; N_{jet tracks}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {200, -0.5, 199.5}, {200, -0.5, 199.5}}});

      registry.add("h3_jet_r_jet_pt_tag_jet_pt_base_matchedpt", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c});#it{p}_{T,jet}^{base} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_eta_tag_jet_eta_base_matchedpt", "#it{R}_{jet};#eta_{jet}^{tag};#eta_{jet}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_phi_tag_jet_phi_base_matchedpt", "#it{R}_{jet};#varphi_{jet}^{tag};#varphi_{jet}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedpt", "#it{R}_{jet};N_{jet tracks}^{tag};N_{jet tracks}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, -0.5, 199.5}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_cluster_pt_tag_cluster_pt_base_matchedpt", "#it{R}_{jet};#it{p}_{T,cluster}^{tag} (GeV/#it{c});#it{p}_{T,cluster}^{base} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_cluster_eta_tag_cluster_eta_base_matchedpt", "#it{R}_{jet};#eta_{cluster}^{tag};#eta_{cluster}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_cluster_phi_tag_cluster_phi_base_matchedpt", "#it{R}_{jet};#varphi_{cluster}^{tag};#varphi_{cluster}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {160, -1.0, 7.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_cluster_energy_tag_cluster_energy_base_matchedpt", "#it{R}_{jet};#E_{cluster}^{tag};#E_{cluster}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedpt", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_eta_base_diff_matchedpt", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#eta_{jet}^{tag} - #eta_{jet}^{base}) / #eta_{jet}^{tag}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_phi_base_diff_matchedpt", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#varphi_{jet}^{tag} - #varphi_{jet}^{base}) / #varphi_{jet}^{tag}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_pt_tag_jet_eta_tag_jet_eta_base_matchedpt", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); #eta_{jet}^{tag}; #eta_{jet}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_pt_tag_jet_phi_tag_jet_phi_base_matchedpt", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); #varphi_{jet}^{tag}; #varphi_{jet}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {160, -1.0, 7.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_pt_tag_jet_ntracks_tag_jet_ntracks_base_matchedpt", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); N_{jet tracks}^{tag}; N_{jet tracks}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {200, -0.5, 199.5}, {200, -0.5, 199.5}}});

      registry.add("h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeopt", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c});#it{p}_{T,jet}^{base} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_eta_tag_jet_eta_base_matchedgeopt", "#it{R}_{jet};#eta_{jet}^{tag};#eta_{jet}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_phi_tag_jet_phi_base_matchedgeopt", "#it{R}_{jet};#varphi_{jet}^{tag};#varphi_{jet}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedgeopt", "#it{R}_{jet};N_{jet tracks}^{tag};N_{jet tracks}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, -0.5, 199.5}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_cluster_pt_tag_cluster_pt_base_matchedgeopt", "#it{R}_{jet};#it{p}_{T,cluster}^{tag} (GeV/#it{c});#it{p}_{T,cluster}^{base} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_cluster_eta_tag_cluster_eta_base_matchedgeopt", "#it{R}_{jet};#eta_{cluster}^{tag};#eta_{cluster}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_cluster_phi_tag_cluster_phi_base_matchedgeopt", "#it{R}_{jet};#varphi_{cluster}^{tag};#varphi_{cluster}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {160, -1.0, 7.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_cluster_energy_tag_cluster_energy_base_matchedgeopt", "#it{R}_{jet};#E_{cluster}^{tag};#E_{cluster}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedgeopt", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_eta_base_diff_matchedgeopt", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#eta_{jet}^{tag} - #eta_{jet}^{base}) / #eta_{jet}^{tag}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_phi_base_diff_matchedgeopt", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#varphi_{jet}^{tag} - #varphi_{jet}^{base}) / #varphi_{jet}^{tag}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_pt_tag_jet_eta_tag_jet_eta_base_matchedgeopt", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); #eta_{jet}^{tag}; #eta_{jet}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_pt_tag_jet_phi_tag_jet_phi_base_matchedgeopt", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); #varphi_{jet}^{tag}; #varphi_{jet}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {160, -1.0, 7.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_pt_tag_jet_ntracks_tag_jet_ntracks_base_matchedgeopt", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); N_{jet tracks}^{tag}; N_{jet tracks}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {200, -0.5, 199.5}, {200, -0.5, 199.5}}});
    }

    if (doprocessTracks || doprocessTracksWeighted) {
      registry.add("h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.add("h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_track_eta", "track #eta;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_track_phi", "track #varphi;#varphi_{track};entries", {HistType::kTH1F, {{160, -1.0, 7.}}});
      registry.add("h_cluster_pt", "cluster pT;#it{p}_{T,cluster} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_cluster_eta", "cluster #eta;#eta_{cluster};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_cluster_phi", "cluster #varphi;#varphi_{cluster};entries", {HistType::kTH1F, {{160, -1.0, 7.}}});
      registry.add("h_cluster_energy", "cluster E;E_{cluster} (GeV);entries", {HistType::kTH1F, {{200, 0., 200.}}});
      if (doprocessTracksWeighted) {
        registry.add("h_collisions_weighted", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      }
    }

    if (doprocessMCCollisionsWeighted) {
      AxisSpec weightAxis = {{VARIABLE_WIDTH, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0}, "weights"};
      registry.add("h_collision_eventweight_part", "event weight;event weight;entries", {HistType::kTH1F, {weightAxis}});
    }
  }

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
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);
  aod::EMCALClusterDefinition clusterDefinition = aod::emcalcluster::getClusterDefinitionFromString(clusterDefinitionS.value);
  Filter clusterFilter = (aod::jcluster::definition == static_cast<int>(clusterDefinition) && aod::jcluster::eta > clusterEtaMin && aod::jcluster::eta < clusterEtaMax && aod::jcluster::phi >= clusterPhiMin && aod::jcluster::phi <= clusterPhiMax && aod::jcluster::energy >= clusterEnergyMin && aod::jcluster::time > clusterTimeMin && aod::jcluster::time < clusterTimeMax && (clusterRejectExotics && aod::jcluster::isExotic != true));

  template <typename T>
  void fillHistograms(T const& jet, float centrality, float weight = 1.0)
  {

    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }

    float neutralEnergy = 0.0;
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("h_jet_ntracks"), jet.tracks().size(), weight);
      registry.fill(HIST("h_jet_nclusters"), jet.clusters().size(), weight);
      registry.fill(HIST("h2_centrality_jet_pt"), centrality, jet.pt(), weight);
      registry.fill(HIST("h2_centrality_jet_eta"), centrality, jet.eta(), weight);
      registry.fill(HIST("h2_centrality_jet_phi"), centrality, jet.phi(), weight);
      registry.fill(HIST("h2_centrality_jet_ntracks"), centrality, jet.tracks().size(), weight);
      registry.fill(HIST("h2_centrality_jet_nclusters"), centrality, jet.clusters().size(), weight);
    }

    registry.fill(HIST("h3_jet_r_jet_pt_centrality"), jet.r() / 100.0, jet.pt(), centrality, weight);
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
      neutralEnergy += cluster.energy();
      registry.fill(HIST("h3_jet_r_jet_pt_cluster_pt"), jet.r() / 100.0, jet.pt(), clusterpt, weight);
      registry.fill(HIST("h3_jet_r_jet_pt_cluster_eta"), jet.r() / 100.0, jet.pt(), cluster.eta(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_cluster_phi"), jet.r() / 100.0, jet.pt(), cluster.phi(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_cluster_energy"), jet.r() / 100.0, jet.pt(), cluster.energy(), weight);
    }
    registry.fill(HIST("h3_jet_r_jet_pt_jet_neutralenergyfraction"), jet.r() / 100.0, jet.pt(), neutralEnergy / jet.energy(), weight);
  }

  template <typename T>
  void fillMCPHistograms(T const& jet, float weight = 1.0)
  {

    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCP * pTHat) {
      return;
    }

    float neutralEnergy = 0.0;
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
      auto pdgParticle = pdgDatabase->GetParticle(constituent.pdgCode());
      if (pdgParticle->Charge() == 0) {
        neutralEnergy += constituent.e();
      }
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_pt_part"), jet.r() / 100.0, jet.pt(), constituent.pt(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_eta_part"), jet.r() / 100.0, jet.pt(), constituent.eta(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_phi_part"), jet.r() / 100.0, jet.pt(), constituent.phi(), weight);
    }
    registry.fill(HIST("h3_jet_r_part_jet_pt_part_jet_neutralenergyfraction_part"), jet.r() / 100.0, jet.pt(), neutralEnergy / jet.energy(), weight);
  }

  template <typename T, typename U>
  void fillMatchedHistograms(T const& jetBase, float weight = 1.0)
  {

    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jetBase.pt() > pTHatMaxMCD * pTHat) {
      return;
    }

    if (jetBase.has_matchedJetGeo()) {
      for (auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<U>>()) {
        if (jetTag.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }
        registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"), jetBase.r() / 100.0, jetTag.pt(), jetBase.pt(), weight);
        registry.fill(HIST("h3_jet_r_jet_eta_tag_jet_eta_base_matchedgeo"), jetBase.r() / 100.0, jetTag.eta(), jetBase.eta(), weight);
        registry.fill(HIST("h3_jet_r_jet_phi_tag_jet_phi_base_matchedgeo"), jetBase.r() / 100.0, jetTag.phi(), jetBase.phi(), weight);
        registry.fill(HIST("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedgeo"), jetBase.r() / 100.0, jetTag.tracks().size(), jetBase.tracks().size(), weight);
        registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedgeo"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.pt() - jetBase.pt()) / jetTag.pt(), weight);
        registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_eta_base_diff_matchedgeo"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.eta() - jetBase.eta()) / jetTag.eta(), weight);
        registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_phi_base_diff_matchedgeo"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.phi() - jetBase.phi()) / jetTag.phi(), weight);

        if (jetBase.r() == round(selectedJetsRadius * 100.0f)) {
          registry.fill(HIST("h3_jet_pt_tag_jet_eta_tag_jet_eta_matchedgeo"), jetTag.pt(), jetTag.eta(), jetBase.eta(), weight);
          registry.fill(HIST("h3_jet_pt_tag_jet_phi_tag_jet_phi_matchedgeo"), jetTag.pt(), jetTag.phi(), jetBase.phi(), weight);
          registry.fill(HIST("h3_jet_pt_tag_jet_ntracks_tag_jet_ntracks_matchedgeo"), jetTag.pt(), jetTag.tracks().size(), jetBase.tracks().size(), weight);
        }
      }
    }

    if (jetBase.has_matchedJetPt()) {
      for (auto& jetTag : jetBase.template matchedJetPt_as<std::decay_t<U>>()) {
        if (jetTag.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }
        registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_matchedpt"), jetBase.r() / 100.0, jetTag.pt(), jetBase.pt(), weight);
        registry.fill(HIST("h3_jet_r_jet_eta_tag_jet_eta_base_matchedpt"), jetBase.r() / 100.0, jetTag.eta(), jetBase.eta(), weight);
        registry.fill(HIST("h3_jet_r_jet_phi_tag_jet_phi_base_matchedpt"), jetBase.r() / 100.0, jetTag.phi(), jetBase.phi(), weight);
        registry.fill(HIST("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedpt"), jetBase.r() / 100.0, jetTag.tracks().size(), jetBase.tracks().size(), weight);
        registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedpt"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.pt() - jetBase.pt()) / jetTag.pt(), weight);
        registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_eta_base_diff_matchedpt"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.eta() - jetBase.eta()) / jetTag.eta(), weight);
        registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_phi_base_diff_matchedpt"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.phi() - jetBase.phi()) / jetTag.phi(), weight);

        if (jetBase.r() == round(selectedJetsRadius * 100.0f)) {
          registry.fill(HIST("h3_jet_pt_tag_jet_eta_tag_jet_eta_base_matchedpt"), jetTag.pt(), jetTag.eta(), jetBase.eta(), weight);
          registry.fill(HIST("h3_jet_pt_tag_jet_phi_tag_jet_phi_base_matchedpt"), jetTag.pt(), jetTag.phi(), jetBase.phi(), weight);
          registry.fill(HIST("h3_jet_pt_tag_jet_ntracks_tag_jet_ntracks_base_matchedpt"), jetTag.pt(), jetTag.tracks().size(), jetBase.tracks().size(), weight);
        }
      }
    }

    if (jetBase.has_matchedJetGeo() && jetBase.has_matchedJetPt()) {
      for (auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<U>>()) {
        if (jetTag.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }
        if (jetBase.template matchedJetGeo_first_as<std::decay_t<U>>().globalIndex() == jetBase.template matchedJetPt_first_as<std::decay_t<U>>().globalIndex()) { // not a good way to do this
          registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeopt"), jetBase.r() / 100.0, jetTag.pt(), jetBase.pt(), weight);
          registry.fill(HIST("h3_jet_r_jet_eta_tag_jet_eta_base_matchedgeopt"), jetBase.r() / 100.0, jetTag.eta(), jetBase.eta(), weight);
          registry.fill(HIST("h3_jet_r_jet_phi_tag_jet_phi_base_matchedgeopt"), jetBase.r() / 100.0, jetTag.phi(), jetBase.phi(), weight);
          registry.fill(HIST("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedgeopt"), jetBase.r() / 100.0, jetTag.tracks().size(), jetBase.tracks().size(), weight);
          registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedgeopt"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.pt() - jetBase.pt()) / jetTag.pt(), weight);
          registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_eta_base_diff_matchedgeopt"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.eta() - jetBase.eta()) / jetTag.eta(), weight);
          registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_phi_base_diff_matchedgeopt"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.phi() - jetBase.phi()) / jetTag.phi(), weight);

          if (jetBase.r() == round(selectedJetsRadius * 100.0f)) {
            registry.fill(HIST("h3_jet_pt_tag_jet_eta_tag_jet_eta_base_matchedpt"), jetTag.pt(), jetTag.eta(), jetBase.eta(), weight);
            registry.fill(HIST("h3_jet_pt_tag_jet_phi_tag_jet_phi_base_matchedpt"), jetTag.pt(), jetTag.phi(), jetBase.phi(), weight);
            registry.fill(HIST("h3_jet_pt_tag_jet_ntracks_tag_jet_ntracks_base_matchedpt"), jetTag.pt(), jetTag.tracks().size(), jetBase.tracks().size(), weight);
          }
        }
      }
    }
  }

  template <typename T, typename U>
  void fillTrackHistograms(T const& tracks, U const& clusters, float weight = 1.0)
  {
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      registry.fill(HIST("h_track_pt"), track.pt(), weight);
      registry.fill(HIST("h_track_eta"), track.eta(), weight);
      registry.fill(HIST("h_track_phi"), track.phi(), weight);
    }
    for (auto const& cluster : clusters) {
      double clusterpt = cluster.energy() / std::cosh(cluster.eta());
      registry.fill(HIST("h_cluster_pt"), clusterpt, weight);
      registry.fill(HIST("h_cluster_eta"), cluster.eta(), weight);
      registry.fill(HIST("h_cluster_phi"), cluster.phi(), weight);
      registry.fill(HIST("h_cluster_energy"), cluster.energy(), weight);
    }
  }

  void processDummy(JetCollisions const& collisions)
  {
  }
  PROCESS_SWITCH(JetFinderFullQATask, processDummy, "dummy task", true);

  void processJetsData(soa::Filtered<JetCollisions>::iterator const& collision, JetTableDataJoined const& jets, JetTracks const& tracks, JetClusters const& clusters)
  {
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      fillHistograms(jet, collision.centrality());
    }
  }
  PROCESS_SWITCH(JetFinderFullQATask, processJetsData, "jet finder HF QA data", false);

  void processJetsMCD(soa::Filtered<JetCollisions>::iterator const& collision, JetTableMCDJoined const& jets, JetTracks const& tracks, JetClusters const& clusters)
  {
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      fillHistograms(jet, collision.centrality());
    }
  }
  PROCESS_SWITCH(JetFinderFullQATask, processJetsMCD, "jet finder HF QA mcd", false);

  void processJetsMCDWeighted(soa::Filtered<JetCollisions>::iterator const& collision, JetTableMCDWeightedJoined const& jets, JetTracks const& tracks, JetClusters const& clusters)
  {
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      fillHistograms(jet, collision.centrality(), jet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetFinderFullQATask, processJetsMCDWeighted, "jet finder HF QA mcd on weighted events", false);

  void processJetsMCP(typename JetTableMCPJoined::iterator const& jet, JetParticles const& particles)
  {
    if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
      return;
    }
    fillMCPHistograms(jet);
  }
  PROCESS_SWITCH(JetFinderFullQATask, processJetsMCP, "jet finder HF QA mcp", false);

  void processJetsMCPWeighted(typename JetTableMCPWeightedJoined::iterator const& jet, JetParticles const& particles)
  {
    if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
      return;
    }
    fillMCPHistograms(jet, jet.eventWeight());
  }
  PROCESS_SWITCH(JetFinderFullQATask, processJetsMCPWeighted, "jet finder HF QA mcp on weighted events", false);

  void processJetsMCPMCDMatched(JetCollision const& collision,
                                JetTableMCDMatchedJoined const& mcdjets,
                                JetTableMCPMatchedJoined const& mcpjets,
                                JetTracks const& tracks,
                                JetClusters const& clusters,
                                JetParticles const& particles)
  {

    for (const auto& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      fillMatchedHistograms<typename JetTableMCDMatchedJoined::iterator, JetTableMCPMatchedJoined>(mcdjet);
    }
  }
  PROCESS_SWITCH(JetFinderFullQATask, processJetsMCPMCDMatched, "jet finder HF QA matched mcp and mcd", false);

  void processJetsMCPMCDMatchedWeighted(JetCollision const& collision,
                                        JetTableMCDMatchedWeightedJoined const& mcdjets,
                                        JetTableMCPMatchedWeightedJoined const& mcpjets,
                                        JetTracks const& tracks,
                                        JetClusters const& clusters,
                                        JetParticles const& particles)
  {

    for (const auto& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      fillMatchedHistograms<typename JetTableMCDMatchedWeightedJoined::iterator, JetTableMCPMatchedWeightedJoined>(mcdjet, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetFinderFullQATask, processJetsMCPMCDMatchedWeighted, "jet finder HF QA matched mcp and mcd on weighted events", false);

  void processMCCollisionsWeighted(JetMcCollision const& collision)
  {

    registry.fill(HIST("h_collision_eventweight_part"), collision.weight());
  }
  PROCESS_SWITCH(JetFinderFullQATask, processMCCollisionsWeighted, "collision QA for weighted events", false);

  void processTracks(JetCollision const& collision,
                     soa::Filtered<JetTracks> const& tracks,
                     soa::Filtered<JetClusters> const& clusters)
  {
    registry.fill(HIST("h_collisions"), 0.5);
    if (!jetderiveddatautilities::eventEMCAL(collision)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    fillTrackHistograms(tracks, clusters);
  }
  PROCESS_SWITCH(JetFinderFullQATask, processTracks, "QA for charged tracks", false);

  void processTracksWeighted(soa::Join<aod::JCollisions, aod::JMcCollisionLbs>::iterator const& collision,
                             aod::JMcCollisions const& mcCollisions,
                             soa::Filtered<JetTracks> const& tracks,
                             soa::Filtered<JetClusters> const& clusters)
  {
    float eventWeight = collision.mcCollision().weight();
    registry.fill(HIST("h_collisions"), 0.5);
    registry.fill(HIST("h_collisions_weighted"), 0.5, eventWeight);
    if (!jetderiveddatautilities::eventEMCAL(collision)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    registry.fill(HIST("h_collisions_weighted"), 1.5, eventWeight);
    fillTrackHistograms(tracks, clusters, eventWeight);
  }
  PROCESS_SWITCH(JetFinderFullQATask, processTracksWeighted, "QA for charged tracks weighted", false);
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
