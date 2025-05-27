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

// FullJet Spectra in pp
//
//
/// \author Archita Rani Dash <archita.rani.dash@cern.ch>
#include <vector>
#include <iostream>
#include <utility>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/Core/HfHelper.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/EMCALMatchedCollisions.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"

#include "EventFiltering/filterTables.h"

using namespace std;
using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
// using namespace jetderiveddatautilities;

struct FullJetSpectrapp {

  HistogramRegistry registry;

  // Event configurables
  Configurable<float> VertexZCut{"VertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};
  Configurable<bool> doEMCALEventWorkaround{"doEMCALEventWorkaround", false, "apply the workaround to read the EMC trigger bit by requiring a cell content in the EMCAL"};
  Configurable<bool> doMBGapTrigger{"doMBGapTrigger", true, "set to true only when using MB-Gap Trigger JJ MC to reject MB events at the collision and track level"};
  // Configurable<bool> doMBMC{"doMBMC", false, "set to true only when using MB MC"};
  Configurable<bool> checkMcCollisionIsMatched{"checkMcCollisionIsMatched", false, "0: count whole MCcollisions, 1: select MCcollisions which only have their correspond collisions"};

  // Jet configurables
  Configurable<float> selectedJetsRadius{"selectedJetsRadius", 0.4, "resolution parameter for histograms without radius"};
  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};
  Configurable<float> jetpTMin{"jetpTMin", 20.0, "minimum jet pT"};
  Configurable<float> jetpTMax{"jetpTMax", 350., "maximum jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -0.3, "minimum jet eta"}; // each of these jet configurables are for the fiducial emcal cuts
  Configurable<float> jetEtaMax{"jetEtaMax", 0.3, "maximum jet eta"};  // for R = 0.4 (EMCAL eta acceptance: eta_jet = 0.7 - R)
  Configurable<float> jetPhiMin{"jetPhiMin", 1.80, "minimum jet phi"}; // phi_jet_min for R = 0.4 is 1.80
  Configurable<float> jetPhiMax{"jetPhiMax", 2.86, "maximum jet phi"}; // phi_jet_min for R = 0.4 is 2.86
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};

  // Track configurables
  Configurable<float> trackpTMin{"trackpTMin", 0.15, "minimum track pT"};
  Configurable<float> trackpTMax{"trackpTMax", 350., "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.7, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.7, "maximum track eta"};
  Configurable<float> trackPhiMin{"trackPhiMin", 1.396, "minimum track phi"};
  Configurable<float> trackPhiMax{"trackPhiMax", 3.283, "maximum track phi"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<std::string> eventSelections{"eventSelections", "selMCFull", "choose event selection"};
  Configurable<std::string> particleSelections{"particleSelections", "PhysicalPrimary", "set particle selections"};

  // Cluster configurables

  Configurable<std::string> clusterDefinitionS{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};
  Configurable<float> clusterEtaMin{"clusterEtaMin", -0.7, "minimum cluster eta"};
  Configurable<float> clusterEtaMax{"clusterEtaMax", 0.7, "maximum cluster eta"};
  Configurable<float> clusterPhiMin{"clusterPhiMin", 1.396, "minimum cluster phi"};
  Configurable<float> clusterPhiMax{"clusterPhiMax", 3.283, "maximum cluster phi"};
  Configurable<float> clusterEnergyMin{"clusterEnergyMin", 0.3, "minimum cluster energy in EMCAL (GeV)"};
  Configurable<float> clusterTimeMin{"clusterTimeMin", -15., "minimum cluster time (ns)"};
  Configurable<float> clusterTimeMax{"clusterTimeMax", 15., "maximum cluster time (ns)"};
  Configurable<bool> clusterRejectExotics{"clusterRejectExotics", true, "Reject exotic clusters"};

  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 4.0, "exponent of the event weight for the calculation of pTHeventSelectionBitsat"}; // 6 for MB MC and 4 for JJ MC
  Configurable<float> pTHatAbsoluteMin{"pTHatAbsoluteMin", -99.0, "minimum value of pTHat"};

  int trackSelection = -1;
  std::vector<int> eventSelectionBits;
  std::vector<bool> filledJetR;
  std::vector<double> jetRadiiValues;

  std::string particleSelection;

  Service<o2::framework::O2DatabasePDG> pdgDatabase;

  // Add Collision Histograms' Bin Labels for clarity
  void labelCollisionHistograms(HistogramRegistry& registry)
  {
    if (doprocessTracks) {
      auto h_collisions_unweighted = registry.get<TH1>(HIST("h_collisions_unweighted"));
      h_collisions_unweighted->GetXaxis()->SetBinLabel(1, "AllUnweightedDetColl");
      h_collisions_unweighted->GetXaxis()->SetBinLabel(2, "UnweightedCollWithVertexZ");
      h_collisions_unweighted->GetXaxis()->SetBinLabel(3, "EMCAcceptedUnweightedColl");
      h_collisions_unweighted->GetXaxis()->SetBinLabel(4, "UnweightedCollAfterTrackSel");
    }

    if (doprocessTracksWeighted) {
      auto h_collisions_weighted = registry.get<TH1>(HIST("h_collisions_weighted"));
      h_collisions_weighted->GetXaxis()->SetBinLabel(1, "AllWeightedDetColl");
      h_collisions_weighted->GetXaxis()->SetBinLabel(2, "WeightedCollWithVertexZ");
      h_collisions_weighted->GetXaxis()->SetBinLabel(3, "EMCAcceptedWeightedColl");
      h_collisions_weighted->GetXaxis()->SetBinLabel(4, "WeightedCollAfterTrackSel");
    }

    if (doprocessJetsData || doprocessJetsMCD || doprocessJetsMCDWeighted) {
      auto h_Detcollision_counter = registry.get<TH1>(HIST("h_Detcollision_counter"));
      h_Detcollision_counter->GetXaxis()->SetBinLabel(1, "allDetColl");
      h_Detcollision_counter->GetXaxis()->SetBinLabel(2, "DetCollWithVertexZ");
      h_Detcollision_counter->GetXaxis()->SetBinLabel(3, "EMCAcceptedDetColl");
    }

    if (doprocessJetsMCP || doprocessJetsMCPWeighted) {
      auto h_Partcollision_counter = registry.get<TH1>(HIST("h_Partcollision_counter"));
      h_Partcollision_counter->GetXaxis()->SetBinLabel(1, "allMcColl");
      h_Partcollision_counter->GetXaxis()->SetBinLabel(2, "McCollWithVertexZ");
      h_Partcollision_counter->GetXaxis()->SetBinLabel(3, "DetCollWithSize>1");
      h_Partcollision_counter->GetXaxis()->SetBinLabel(4, "EMCAcceptedDetColl");
    }

    if (doprocessJetsMCPMCDMatched || doprocessJetsMCPMCDMatchedWeighted) {
      auto h_Matchedcollision_counter = registry.get<TH1>(HIST("h_Matchedcollision_counter"));
      h_Matchedcollision_counter->GetXaxis()->SetBinLabel(1, "allDetColl");
      h_Matchedcollision_counter->GetXaxis()->SetBinLabel(2, "DetCollWithVertexZ");
      h_Matchedcollision_counter->GetXaxis()->SetBinLabel(3, "EMCAcceptedDetColl");
    }
  }

  void init(o2::framework::InitContext&)
  {
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    particleSelection = static_cast<std::string>(particleSelections);
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

    // Track QA histograms
    if (doprocessTracks || doprocessTracksWeighted) {
      registry.add("h_collisions_unweighted", "event status; event status;entries", {HistType::kTH1F, {{12, 0., 12.0}}});

      registry.add("h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_track_eta", "track #eta;#eta_{track};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_track_phi", "track #varphi;#varphi_{track};entries", {HistType::kTH1F, {{160, 0., 7.}}});
      registry.add("h_track_energy", "track energy;Energy of tracks;entries", {HistType::kTH1F, {{400, 0., 400.}}});
      registry.add("h_track_energysum", "track energy sum;Sum of track energy per event;entries", {HistType::kTH1F, {{400, 0., 400.}}});

      // Cluster QA histograms
      registry.add("h_clusterTime", "Time of cluster", HistType::kTH1F, {{500, -250, 250, "#it{t}_{cls} (ns)"}});
      registry.add("h_cluster_pt", "cluster pT;#it{p}_{T_cluster} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_cluster_eta", "cluster #eta;#eta_{cluster};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_cluster_phi", "cluster #varphi;#varphi_{cluster};entries", {HistType::kTH1F, {{160, 0., 7.}}});
      registry.add("h_cluster_energy", "cluster energy;Energy of cluster;entries", {HistType::kTH1F, {{400, 0., 400.}}});
      registry.add("h_cluster_energysum", "cluster energy sum;Sum of cluster energy per event;entries", {HistType::kTH1F, {{400, 0., 400.}}});

      if (doprocessTracksWeighted) {
        registry.add("h_collisions_weighted", "event status;event status;entries", {HistType::kTH1F, {{12, 0.0, 12.0}}});
      }
    }

    // Jet QA histograms
    if (doprocessJetsData || doprocessJetsMCD || doprocessJetsMCDWeighted) {
      registry.add("h_Detcollision_counter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.}}});

      registry.add("h_full_jet_pt", "#it{p}_{T,jet};#it{p}_{T_jet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_full_jet_pt_pTHatcut", "#it{p}_{T,jet};#it{p}_{T_jet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_full_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_full_jet_phi", "jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., 7.}}});
      registry.add("h_full_jet_clusterTime", "Time of cluster", HistType::kTH1F, {{500, -250, 250, "#it{t}_{cls} (ns)"}});
      registry.add("h2_full_jet_NEF", "#it{p}_{T,jet} vs NEF at Det Level; #it{p}_{T,jet} (GeV/#it{c});NEF", {HistType::kTH2F, {{350, 0., 350.}, {105, 0., 1.05}}});
      registry.add("h2_full_jet_NEF_rejected", "#it{p}_{T,jet} vs NEF at Det Level for rejected events; #it{p}_{T,jet} (GeV/#it{c});NEF", {HistType::kTH2F, {{350, 0., 350.}, {105, 0., 1.05}}});

      registry.add("h_Detjet_ntracks", "#it{p}_{T,track};#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h2_full_jet_chargedconstituents", "Number of charged constituents at Det Level;#it{p}_{T,jet} (GeV/#it{c});N_{ch}", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 100.}}});
      registry.add("h2_full_jet_neutralconstituents", "Number of neutral constituents at Det Level;#it{p}_{T,jet} (GeV/#it{c});N_{ne}", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 100.}}});
      registry.add("h_full_jet_chargedconstituents_pt", "track pT;#it{p}^{T,jet}_{track} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_full_jet_chargedconstituents_eta", "track #eta;#eta^{jet}_{track};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_full_jet_chargedconstituents_phi", "track #varphi;#varphi^{jet}_{track};entries", {HistType::kTH1F, {{160, 0., 7.}}});
      registry.add("h_full_jet_chargedconstituents_energy", "track energy;Energy of tracks;entries", {HistType::kTH1F, {{400, 0., 400.}}});
      registry.add("h_full_jet_chargedconstituents_energysum", "track energy sum;Sum of track energy per event;entries", {HistType::kTH1F, {{400, 0., 400.}}});
      registry.add("h_full_jet_neutralconstituents_pt", "cluster pT;#it{p}^{T,jet}_{cluster} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_full_jet_neutralconstituents_eta", "cluster #eta;#eta^{jet}_{cluster};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_full_jet_neutralconstituents_phi", "cluster #varphi;#varphi^{jet}_{cluster};entries", {HistType::kTH1F, {{160, 0., 7.}}});
      registry.add("h_full_jet_neutralconstituents_energy", "cluster energy;Energy of cluster;entries", {HistType::kTH1F, {{400, 0., 400.}}});
      registry.add("h_full_jet_neutralconstituents_energysum", "cluster energy sum;Sum of cluster energy per event;entries", {HistType::kTH1F, {{400, 0., 400.}}});
      registry.add("h2_full_jettrack_pt", "#it{p}_{T,jet} vs #it{p}_{T,track}; #it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{350, 0., 350.}, {200, 0., 200.}}});
      registry.add("h2_full_jettrack_eta", "jet #eta vs jet_track #eta; #eta_{jet};#eta_{track}", {HistType::kTH2F, {{100, -1., 1.}, {500, -5., 5.}}});
      registry.add("h2_full_jettrack_phi", "jet #varphi vs jet_track #varphi; #varphi_{jet}; #varphi_{track}", {HistType::kTH2F, {{160, 0., 7.}, {160, -1., 7.}}});

      registry.add("h2_track_etaphi", "jet_track #eta vs jet_track #varphi; #eta_{track};#varphi_{track}", {HistType::kTH2F, {{500, -5., 5.}, {160, -1., 7.}}});
      registry.add("h2_jet_etaphi", "jet #eta vs jet #varphi; #eta_{jet};#varphi_{jet}", {HistType::kTH2F, {{100, -1., 1.}, {160, -1., 7.}}});
    }
    if (doprocessJetsMCP || doprocessJetsMCPWeighted) {
      registry.add("h_Partcollision_counter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});

      registry.add("h_full_mcpjet_tablesize", "", {HistType::kTH1F, {{4, 0., 5.}}});
      registry.add("h_full_mcpjet_ntracks", "", {HistType::kTH1F, {{200, -0.5, 200.}}});
      registry.add("h_full_jet_pt_part", "jet pT;#it{p}_{T_jet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_full_jet_eta_part", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_full_jet_phi_part", "jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., 7.}}});
      registry.add("h2_full_jet_NEF_part", "#it{p}_{T,jet} vs NEF at Part Level;#it{p}_{T,jet} (GeV/#it{c});NEF", {HistType::kTH2F, {{350, 0., 350.}, {105, 0., 1.05}}});

      registry.add("h_Partjet_ntracks", "#it{p}_{T,constituent};#it{p}_{T_constituent} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h2_full_jet_chargedconstituents_part", "Number of charged constituents at Part Level;#it{p}_{T,jet} (GeV/#it{c});N_{ch}", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 100.}}});
      registry.add("h2_full_jet_neutralconstituents_part", "Number of neutral constituents at Part Level;#it{p}_{T,jet} (GeV/#it{c});N_{ne}", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 100.}}});
      registry.add("h_full_jet_neutralconstituents_pt_part", "#it{p}_{T} of neutral constituents at Part Level;#it{p}_{T,ne} (GeV/#it{c}); entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_full_jet_neutralconstituents_eta_part", "#eta of neutral constituents at Part Level;#eta_{ne};entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_full_jet_neutralconstituents_phi_part", "#varphi of neutral constituents at Part Level;#varphi_{ne};entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_full_jet_neutralconstituents_energy_part", "neutral constituents' energy;Energy of neutral constituents;entries", {HistType::kTH1F, {{400, 0., 400.}}});
      registry.add("h_full_jet_neutralconstituents_energysum_part", "neutral constituents' energy sum;Sum of neutral constituents' energy per event;entries", {HistType::kTH1F, {{400, 0., 400.}}});

      registry.add("h2_jettrack_pt_part", "#it{p}_{T,jet} vs #it{p}_{T_track}; #it{p}_{T_jet} (GeV/#it{c});#it{p}_{T_track} (GeV/#it{c})", {HistType::kTH2F, {{350, 0., 350.}, {200, 0., 200.}}});
      registry.add("h2_jettrack_eta_part", "jet #eta vs jet_track #eta; #eta_{jet};#eta_{track}", {HistType::kTH2F, {{100, -1., 1.}, {500, -5., 5.}}});
      registry.add("h2_jettrack_phi_part", "jet #varphi vs jet_track #varphi; #varphi_{jet}; #varphi_{track}", {HistType::kTH2F, {{160, 0., 7.}, {160, -1., 7.}}});

      registry.add("h2_track_etaphi_part", "jet_track #eta vs jet_track #varphi; #eta_{track};#varphi_{track}", {HistType::kTH2F, {{500, -5., 5.}, {160, -1., 7.}}});
      registry.add("h2_jet_etaphi_part", "jet #eta vs jet #varphi; #eta_{jet};#varphi_{jet}", {HistType::kTH2F, {{100, -1., 1.}, {160, -1., 7.}}});

      registry.add("h_NOmcpemcalcollisions", "event status;entries", {HistType::kTH1F, {{100, 0., 100.}}});
      registry.add("h_mcpemcalcollisions", "event status;entries", {HistType::kTH1F, {{100, 0., 100.}}});
    }

    if (doprocessJetsMCPMCDMatched || doprocessJetsMCPMCDMatchedWeighted) {
      registry.add("h_Matchedcollision_counter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});

      registry.add("h_full_matchedmcdjet_tablesize", "", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_full_matchedmcpjet_tablesize", "", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_full_matchedmcdjet_ntracks", "", {HistType::kTH1F, {{200, -0.5, 200.}}});
      registry.add("h_full_matchedmcpjet_ntracks", "", {HistType::kTH1F, {{200, -0.5, 200.}}});
      registry.add("h_full_matchedmcdjet_eta", "Matched MCD jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_full_matchedmcdjet_phi", "Matched MCD jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., 7.}}});
      registry.add("h_full_matchedmcpjet_eta", "Matched MCP jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_full_matchedmcpjet_phi", "Matched MCP jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., 7.}}});
      registry.add("h_full_jet_deltaR", "Distance between matched Det Jet and Part Jet; #Delta R; entries", {HistType::kTH1F, {{100, 0., 1.}}});

      registry.add("h2_full_jet_energyscaleDet", "Jet Energy Scale (det); p_{T,det} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH2F, {{400, 0., 400.}, {200, -1., 1.}}});

      registry.add("h2_matchedjet_etaphiDet", "Det jet #eta vs jet #varphi; #eta_{jet};#varphi_{jet}", {HistType::kTH2F, {{100, -1., 1.}, {160, -1., 7.}}});
      registry.add("h2_matchedjet_etaphiPart", "Part jet #eta vs jet #varphi; #eta_{jet};#varphi_{jet}", {HistType::kTH2F, {{100, -1., 1.}, {160, -1., 7.}}});
      registry.add("h2_matchedjet_deltaEtaCorr", "Correlation between Det Eta and Part Eta; #eta_{jet,det}; #eta_{jet,part}", {HistType::kTH2F, {{100, -1., 1.}, {100, -1., 1.}}});
      registry.add("h2_matchedjet_deltaPhiCorr", "Correlation between Det Phi and Part Phi; #varphi_{jet,det}; #varphi_{jet,part}", {HistType::kTH2F, {{160, 0., 7.}, {160, 0., 7.}}});

      registry.add("h2_full_jet_energyscalePart", "Jet Energy Scale (part); p_{T,part} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH2F, {{400, 0., 400.}, {200, -1., 1.}}});
      registry.add("h3_full_jet_energyscalePart", "R dependence of Jet Energy Scale (Part); #it{R}_{jet};p_{T,det} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, 0., 400.}, {200, -1., 1.}}});
      registry.add("h2_full_jet_etaresolutionPart", ";p_{T,part} (GeV/c); (#eta_{jet,det} - #eta_{jet,part})/#eta_{jet,part}", {HistType::kTH2F, {{400, 0., 400.}, {100, -1., 1.}}});
      registry.add("h2_full_jet_phiresolutionPart", ";p_{T,part} (GeV/c); (#varphi_{jet,det} - #varphi_{jet,part})/#varphi_{jet,part}", {HistType::kTH2F, {{400, 0., 400.}, {160, -1., 7.}}});
      registry.add("h2_full_jet_energyscaleChargedPart", "Jet Energy Scale (charged part); p_{T,part} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH2F, {{400, 0., 400.}, {200, -1., 1.}}});
      registry.add("h2_full_jet_energyscaleNeutralPart", "Jet Energy Scale (neutral part); p_{T,part} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH2F, {{400, 0., 400.}, {200, -1., 1.}}});
      registry.add("h2_full_jet_energyscaleChargedVsFullPart", "Jet Energy Scale (charged part, vs. full jet pt); p_{T,part} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH2F, {{400, 0., 400.}, {200, -1., 1.}}});
      registry.add("h2_full_jet_energyscaleNeutralVsFullPart", "Jet Energy Scale (neutral part, vs. full jet pt); p_{T,part} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH2F, {{400, 0., 400.}, {200, -1., 1.}}});
      registry.add("h2_full_fakemcdjets", "Fake MCD Jets; p_{T,det} (GeV/c); NCounts", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 100.}}});
      registry.add("h2_full_fakemcpjets", "Fake MCP Jets; p_{T,part} (GeV/c); NCounts", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 100.}}});
      // Response Matrix
      registry.add("h_full_jet_ResponseMatrix", "Full Jets Response Matrix; p_{T,det} (GeV/c); p_{T,part} (GeV/c)", {HistType::kTH2F, {{350, 0., 350.}, {350, 0., 350.}}});
    }

    if (doprocessCollisionsWeightedWithMultiplicity || doprocessMBCollisionsWithMultiplicity) {
      registry.add("h_FT0Mults_occupancy", "", {HistType::kTH1F, {{3500, 0., 3500.}}});
      registry.add("h2_full_jet_FT0Amplitude", "; FT0C Amplitude; Counts", {HistType::kTH1F, {{3500, 0., 3500.}}});
      registry.add("h2_full_jet_jetpTDetVsFT0Mults", "; p_{T,det} (GeV/c); FT0C Multiplicity", {HistType::kTH2F, {{350, 0., 350.}, {3500, 0., 3500.}}});
      registry.add("h3_full_jet_jetpTDet_FT0Mults_NEF", "; p_{T,det} (GeV/c); FT0C Multiplicity, NEF", {HistType::kTH3F, {{350, 0., 350.}, {3500, 0., 3500.}, {105, 0.0, 1.05}}});
    }

    // Label the histograms
    labelCollisionHistograms(registry);

  } // init

  using EMCCollisionsData = o2::soa::Join<aod::JetCollisions, aod::JEMCCollisionLbs>;   // JetCollisions with EMCAL Collision Labels
  using EMCCollisionsMCD = o2::soa::Join<aod::JetCollisionsMCD, aod::JEMCCollisionLbs>; // where, JetCollisionsMCD = JetCollisions+JMcCollisionLbs

  using FullJetTableDataJoined = soa::Join<aod::FullJets, aod::FullJetConstituents>;
  using JetTableMCDJoined = soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents>;
  using JetTableMCDWeightedJoined = soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents, aod::FullMCDetectorLevelJetEventWeights>;
  using JetTableMCPJoined = soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents>;
  using JetTableMCPWeightedJoined = soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents, aod::FullMCParticleLevelJetEventWeights>;

  using JetTableMCDMatchedJoined = soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents, aod::FullMCDetectorLevelJetsMatchedToFullMCParticleLevelJets>;
  using JetTableMCPMatchedJoined = soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents, aod::FullMCParticleLevelJetsMatchedToFullMCDetectorLevelJets>;

  using JetTableMCDMatchedWeightedJoined = soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents, aod::FullMCDetectorLevelJetsMatchedToFullMCParticleLevelJets, aod::FullMCDetectorLevelJetsMatchedToFullMCParticleLevelJets, aod::FullMCDetectorLevelJetEventWeights>;
  using JetTableMCPMatchedWeightedJoined = soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents, aod::FullMCParticleLevelJetsMatchedToFullMCDetectorLevelJets, aod::FullMCParticleLevelJetsMatchedToFullMCDetectorLevelJets, aod::FullMCParticleLevelJetEventWeights>;

  // Applying some cuts(filters) on collisions, tracks, clusters

  Filter eventCuts = (nabs(aod::jcollision::posZ) < VertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);
  // Filter EMCeventCuts = (nabs(aod::collision::posZ) < VertexZCut && aod::collision::centrality >= centralityMin && aod::collision::centrality < centralityMax);
  Filter trackCuts = (aod::jtrack::pt >= trackpTMin && aod::jtrack::pt < trackpTMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax && aod::jtrack::phi >= trackPhiMin && aod::jtrack::phi <= trackPhiMax);
  aod::EMCALClusterDefinition clusterDefinition = aod::emcalcluster::getClusterDefinitionFromString(clusterDefinitionS.value);
  Filter clusterFilter = (aod::jcluster::definition == static_cast<int>(clusterDefinition) && aod::jcluster::eta > clusterEtaMin && aod::jcluster::eta < clusterEtaMax && aod::jcluster::phi >= clusterPhiMin && aod::jcluster::phi <= clusterPhiMax && aod::jcluster::energy >= clusterEnergyMin && aod::jcluster::time > clusterTimeMin && aod::jcluster::time < clusterTimeMax && (clusterRejectExotics && aod::jcluster::isExotic != true));
  Preslice<JetTableMCPMatchedJoined> JetMCPPerMcCollision = aod::jet::mcCollisionId;
  PresliceUnsorted<soa::Filtered<aod::JetCollisionsMCD>> CollisionsPerMCPCollision = aod::jmccollisionlb::mcCollisionId;

  template <typename T, typename U>
  bool isAcceptedJet(U const& jet)
  {

    if (jetAreaFractionMin > -98.0) {
      if (jet.area() < jetAreaFractionMin * M_PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    if (leadingConstituentPtMin > -98.0) {
      bool isMinleadingConstituent = false;
      for (auto& constituent : jet.template tracks_as<T>()) {
        if (constituent.pt() >= leadingConstituentPtMin) {
          isMinleadingConstituent = true;
          break;
        }
      }

      if (!isMinleadingConstituent) {
        return false;
      }
    }
    return true;
  }
  template <typename T>
  void fillJetHistograms(T const& jet, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin) { // for MCD jets only to remove outliers; setting pTHatMaxMCD = 1 improves purity
      return;
    }

    float neutralEnergy = 0.0;
    double sumtrackE = 0.0;
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_full_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_full_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_full_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("h2_jet_etaphi"), jet.eta(), jet.phi(), weight);

      for (auto& cluster : jet.template clusters_as<aod::JetClusters>()) {
        registry.fill(HIST("h2_full_jet_neutralconstituents"), jet.pt(), jet.clustersIds().size(), weight);

        neutralEnergy += cluster.energy();
        double clusterpt = cluster.energy() / std::cosh(cluster.eta());
        registry.fill(HIST("h_full_jet_clusterTime"), cluster.time(), weight);
        registry.fill(HIST("h_full_jet_neutralconstituents_pt"), clusterpt, weight);
        registry.fill(HIST("h_full_jet_neutralconstituents_eta"), cluster.eta(), weight);
        registry.fill(HIST("h_full_jet_neutralconstituents_phi"), cluster.phi(), weight);
        registry.fill(HIST("h_full_jet_neutralconstituents_energy"), cluster.energy(), weight);
        registry.fill(HIST("h_full_jet_neutralconstituents_energysum"), neutralEnergy, weight);
      }
      auto NEF = neutralEnergy / jet.energy();
      registry.fill(HIST("h2_full_jet_NEF"), jet.pt(), NEF, weight);

      for (auto& jettrack : jet.template tracks_as<aod::JetTracks>()) {
        sumtrackE += jettrack.energy();

        registry.fill(HIST("h_Detjet_ntracks"), jettrack.pt(), weight);
        registry.fill(HIST("h2_full_jet_chargedconstituents"), jet.pt(), jet.tracksIds().size(), weight);
        registry.fill(HIST("h2_full_jettrack_pt"), jet.pt(), jettrack.pt(), weight);
        registry.fill(HIST("h2_full_jettrack_eta"), jet.eta(), jettrack.eta(), weight);
        registry.fill(HIST("h2_full_jettrack_phi"), jet.phi(), jettrack.phi(), weight);

        registry.fill(HIST("h2_track_etaphi"), jettrack.eta(), jettrack.phi(), weight);
        registry.fill(HIST("h_full_jet_chargedconstituents_pt"), jettrack.pt(), weight);
        registry.fill(HIST("h_full_jet_chargedconstituents_eta"), jettrack.eta(), weight);
        registry.fill(HIST("h_full_jet_chargedconstituents_phi"), jettrack.phi(), weight);
        registry.fill(HIST("h_full_jet_chargedconstituents_energy"), jettrack.energy(), weight);
        registry.fill(HIST("h_full_jet_chargedconstituents_energysum"), sumtrackE, weight);
      }
    } // jet.r()
  }

  // check for NEF distribution for rejected events
  template <typename T>
  void fillRejectedJetHistograms(T const& jet, float weight = 1.0)
  {
    float neutralEnergy = 0.0;
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      for (auto& cluster : jet.template clusters_as<aod::JetClusters>()) {
        neutralEnergy += cluster.energy();
      }
      auto NEF = neutralEnergy / jet.energy();
      registry.fill(HIST("h2_full_jet_NEF_rejected"), jet.pt(), NEF, weight);
    } // jet.r()
  }

  template <typename T>
  void fillMCPHistograms(T const& jet, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin) { // MCP outlier rejection
      return;
    }
    float neutralEnergy = 0.0;
    int neutralconsts = 0;
    int chargedconsts = 0;
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_full_mcpjet_tablesize"), jet.size(), weight);
      registry.fill(HIST("h_full_mcpjet_ntracks"), jet.tracksIds().size(), weight);
      registry.fill(HIST("h_full_jet_pt_part"), jet.pt(), weight);
      registry.fill(HIST("h_full_jet_eta_part"), jet.eta(), weight);
      registry.fill(HIST("h_full_jet_phi_part"), jet.phi(), weight);
      registry.fill(HIST("h2_jet_etaphi_part"), jet.eta(), jet.phi(), weight);

      for (auto& constituent : jet.template tracks_as<aod::JetParticles>()) {
        auto pdgParticle = pdgDatabase->GetParticle(constituent.pdgCode());
        if (pdgParticle->Charge() == 0) {
          neutralconsts++;
          neutralEnergy += constituent.e(); // neutral jet constituents at particle level
          double clusterpt = constituent.e() / std::cosh(constituent.eta());
          registry.fill(HIST("h2_full_jet_neutralconstituents_part"), jet.pt(), neutralconsts, weight);
          registry.fill(HIST("h_full_jet_neutralconstituents_pt_part"), clusterpt, weight);
          registry.fill(HIST("h_full_jet_neutralconstituents_eta_part"), constituent.eta(), weight);
          registry.fill(HIST("h_full_jet_neutralconstituents_phi_part"), constituent.phi(), weight);
          registry.fill(HIST("h_full_jet_neutralconstituents_energy_part"), constituent.e(), weight);
          registry.fill(HIST("h_full_jet_neutralconstituents_energysum_part"), neutralEnergy, weight);

        } else {
          chargedconsts++;
          registry.fill(HIST("h2_full_jet_chargedconstituents_part"), jet.pt(), chargedconsts, weight); // charged jet constituents at particle level
          registry.fill(HIST("h2_jettrack_pt_part"), jet.pt(), constituent.pt(), weight);
          registry.fill(HIST("h2_jettrack_eta_part"), jet.eta(), constituent.eta(), weight);
          registry.fill(HIST("h2_jettrack_phi_part"), jet.phi(), constituent.phi(), weight);
          registry.fill(HIST("h2_track_etaphi_part"), constituent.eta(), constituent.phi(), weight);
        }
      } // constituent loop
      auto NEF = neutralEnergy / jet.energy();
      registry.fill(HIST("h2_full_jet_NEF_part"), jet.pt(), NEF, weight);
    }
  }

  template <typename T, typename U>
  void fillTrackHistograms(T const& tracks, U const& clusters, float weight = 1.0)
  {
    double sumtrackE = 0.0;
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (pTHat < pTHatAbsoluteMin) { // Track outlier rejection
      return;
    }
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      sumtrackE += track.energy();
      registry.fill(HIST("h_track_pt"), track.pt(), weight);
      registry.fill(HIST("h_track_eta"), track.eta(), weight);
      registry.fill(HIST("h_track_phi"), track.phi(), weight);
      registry.fill(HIST("h_track_energysum"), sumtrackE, weight);
    }
    double sumclusterE = 0.0;
    for (auto const& cluster : clusters) {
      double clusterpt = cluster.energy() / std::cosh(cluster.eta());
      sumclusterE += cluster.energy();

      registry.fill(HIST("h_clusterTime"), cluster.time());
      registry.fill(HIST("h_cluster_pt"), clusterpt, weight);
      registry.fill(HIST("h_cluster_eta"), cluster.eta(), weight);
      registry.fill(HIST("h_cluster_phi"), cluster.phi(), weight);
      registry.fill(HIST("h_cluster_energy"), cluster.energy(), weight);
      registry.fill(HIST("h_cluster_energysum"), sumclusterE, weight);
    }
  }

  template <typename T, typename U>
  void fillMatchedHistograms(T const& jetBase, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jetBase.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin) {
      return;
    }

    if (jetBase.has_matchedJetGeo()) { // geometrical jet matching only needed for pp - here,matching Base(Det.level) with Tag (Part. level) jets
      registry.fill(HIST("h_full_matchedmcdjet_tablesize"), jetBase.size(), weight);
      registry.fill(HIST("h_full_matchedmcdjet_ntracks"), jetBase.tracksIds().size(), weight);
      registry.fill(HIST("h2_matchedjet_etaphiDet"), jetBase.eta(), jetBase.phi(), weight);

      for (auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<U>>()) {
        if (jetTag.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin) { // MCP outlier rejection
          continue;
        }
        auto deltaEta = jetBase.eta() - jetTag.eta();
        auto deltaPhi = jetBase.phi() - jetTag.phi();
        auto deltaR = jetutilities::deltaR(jetBase, jetTag);

        registry.fill(HIST("h_full_jet_deltaR"), deltaR, weight);
        registry.fill(HIST("h_full_matchedmcpjet_tablesize"), jetTag.size(), weight);
        registry.fill(HIST("h_full_matchedmcpjet_ntracks"), jetTag.tracksIds().size(), weight);
        registry.fill(HIST("h2_matchedjet_etaphiPart"), jetTag.eta(), jetTag.phi(), weight);
        registry.fill(HIST("h2_matchedjet_deltaEtaCorr"), jetBase.eta(), jetTag.eta(), weight);
        registry.fill(HIST("h2_matchedjet_deltaPhiCorr"), jetBase.phi(), jetTag.phi(), weight);

        // JES for fulljets
        registry.fill(HIST("h2_full_jet_energyscaleDet"), jetBase.pt(), (jetBase.pt() - jetTag.pt()) / jetTag.pt(), weight);
        registry.fill(HIST("h2_full_jet_energyscalePart"), jetTag.pt(), (jetBase.pt() - jetTag.pt()) / jetTag.pt(), weight);
        registry.fill(HIST("h3_full_jet_energyscalePart"), jetBase.r() / 100.0, jetTag.pt(), (jetBase.pt() - jetTag.pt()) / jetTag.pt(), weight);
        registry.fill(HIST("h2_full_jet_etaresolutionPart"), jetTag.pt(), deltaEta / jetTag.eta(), weight);
        registry.fill(HIST("h2_full_jet_phiresolutionPart"), jetTag.pt(), deltaPhi / jetTag.phi(), weight);

        // Response Matrix
        registry.fill(HIST("h_full_jet_ResponseMatrix"), jetBase.pt(), jetTag.pt(), weight); // MCD vs MCP jet pT
      } // jetTag
    } // jetBase
  }

  void processDummy(aod::JetCollisions const&)
  {
  }
  PROCESS_SWITCH(FullJetSpectrapp, processDummy, "dummy task", true);

  void processJetsData(soa::Filtered<EMCCollisionsData>::iterator const& collision, FullJetTableDataJoined const& jets, aod::JetTracks const&, aod::JetClusters const&)
  {
    bool eventAccepted = false;
    registry.fill(HIST("h_Detcollision_counter"), 0.5);

    if (fabs(collision.posZ()) > VertexZCut) {
      return;
    }
    registry.fill(HIST("h_Detcollision_counter"), 1.5);
    if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      return;
    }
    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        if (collision.alias_bit(kTVXinEMC)) {
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
      }
    }

    if (!eventAccepted) {
      for (auto const& jet : jets) {
        if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax) || !isAcceptedJet<aod::JetTracks>(jet)) {
          fillRejectedJetHistograms(jet, 1.0);
        }
      }
      return;
    }
    registry.fill(HIST("h_Detcollision_counter"), 2.5);

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (jet.phi() < jetPhiMin || jet.phi() > jetPhiMax) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      fillJetHistograms(jet);
    }
  }
  PROCESS_SWITCH(FullJetSpectrapp, processJetsData, "Full Jets Data", false);

  void processJetsMCD(soa::Filtered<EMCCollisionsMCD>::iterator const& collision, JetTableMCDJoined const& jets, aod::JetTracks const&, aod::JetClusters const&)
  {
    bool eventAccepted = false;

    registry.fill(HIST("h_Detcollision_counter"), 0.5);
    if (fabs(collision.posZ()) > VertexZCut) {
      return;
    }
    registry.fill(HIST("h_Detcollision_counter"), 1.5);
    if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      return;
    }
    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        if (collision.alias_bit(kTVXinEMC)) {
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
      }
    }

    if (!eventAccepted) {
      for (auto const& jet : jets) {
        if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax) || !isAcceptedJet<aod::JetTracks>(jet)) {
          fillRejectedJetHistograms(jet, 1.0);
        }
      }
      return;
    }
    registry.fill(HIST("h_Detcollision_counter"), 2.5);

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (jet.phi() < jetPhiMin || jet.phi() > jetPhiMax) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      fillJetHistograms(jet);
    }
  }
  PROCESS_SWITCH(FullJetSpectrapp, processJetsMCD, "Full Jets at Detector Level", false);

  void processJetsMCDWeighted(soa::Filtered<EMCCollisionsMCD>::iterator const& collision, JetTableMCDWeightedJoined const& jets, aod::JMcCollisions const&, aod::JetTracks const&, aod::JetClusters const&)
  {
    bool eventAccepted = false;

    registry.fill(HIST("h_Detcollision_counter"), 0.5);
    if (fabs(collision.posZ()) > VertexZCut) {
      return;
    }
    registry.fill(HIST("h_Detcollision_counter"), 1.5);
    if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      return;
    }

    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        if (collision.alias_bit(kTVXinEMC)) {
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
      }
    }

    if (!eventAccepted) {
      for (auto const& jet : jets) {
        if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax) || !isAcceptedJet<aod::JetTracks>(jet)) {
          fillRejectedJetHistograms(jet, collision.mcCollision().weight());
        }
      }
      return;
    }
    registry.fill(HIST("h_Detcollision_counter"), 2.5);

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (jet.phi() < jetPhiMin || jet.phi() > jetPhiMax) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }

      // this cut only to be used for calculating Jet Purity and not for Response Matrix
      // this is mainly applied to remove all high weight jets causing big fluctuations
      double pTHat = 10. / (std::pow(collision.mcCollision().weight(), 1.0 / pTHatExponent));
      if (jet.pt() > 1 * pTHat) {
        registry.fill(HIST("h_full_jet_pt_pTHatcut"), jet.pt(), collision.mcCollision().weight());
      }

      fillJetHistograms(jet, collision.mcCollision().weight());
    }
  }
  PROCESS_SWITCH(FullJetSpectrapp, processJetsMCDWeighted, "Full Jets at Detector Level on weighted events", false);

  void processJetsMCP(aod::JetMcCollision const& mccollision, JetTableMCPJoined const& jets, aod::JetParticles const&, soa::SmallGroups<EMCCollisionsMCD> const& collisions)
  {
    bool eventAccepted = false;

    registry.fill(HIST("h_Partcollision_counter"), 0.5);
    if (fabs(mccollision.posZ()) > VertexZCut) {
      return;
    }
    registry.fill(HIST("h_Partcollision_counter"), 1.5);
    if (collisions.size() < 1) {
      return;
    }
    registry.fill(HIST("h_Partcollision_counter"), 2.5);
    for (auto const& collision : collisions) {
      if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
        return;
      }
      if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
        eventAccepted = true;
        // return;
      }
      if (doEMCALEventWorkaround) {
        if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
          eventAccepted = true;
          if (collision.alias_bit(kTVXinEMC)) {
          }
        }
      } else {
        if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
          eventAccepted = true;
        }
      }
    }
    if (!eventAccepted) {
      return;
    }
    registry.fill(HIST("h_Partcollision_counter"), 3.5);

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetParticles>(jet)) {
        continue;
      }
      if (checkMcCollisionIsMatched) { // basically checks if the same collisions are generated at the Part level as those at the Det level
        auto collisionspermcpjet = collisions.sliceBy(CollisionsPerMCPCollision, jet.mcCollisionId());
        if (collisionspermcpjet.size() >= 1 && jetderiveddatautilities::selectCollision(collisionspermcpjet.begin(), eventSelectionBits)) {
          // Now here for every matched collision, I fill the corresponding jet histograms.
          fillMCPHistograms(jet);
        }
      } else {
        fillMCPHistograms(jet);
      }
    }
  }
  PROCESS_SWITCH(FullJetSpectrapp, processJetsMCP, "Full Jets at Particle Level", false);

  void processJetsMCPWeighted(aod::JetMcCollision const& mccollision, JetTableMCPWeightedJoined const& jets, aod::JetParticles const&, soa::SmallGroups<EMCCollisionsMCD> const& collisions)
  {
    bool eventAccepted = false;

    registry.fill(HIST("h_Partcollision_counter"), 0.5);
    if (fabs(mccollision.posZ()) > VertexZCut) {
      return;
    }
    registry.fill(HIST("h_Partcollision_counter"), 1.5);
    if (collisions.size() < 1) {
      return;
    }
    registry.fill(HIST("h_Partcollision_counter"), 2.5);
    for (auto const& collision : collisions) {
      if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
        return;
      }
      if (jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
        eventAccepted = true;
      }
      if (doEMCALEventWorkaround) {
        if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
          eventAccepted = true;
          if (collision.alias_bit(kTVXinEMC)) {
          }
        }
      } else {
        if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
          eventAccepted = true;
        }
      }
    }
    if (!eventAccepted) {
      return;
    }
    registry.fill(HIST("h_Partcollision_counter"), 3.5);

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        return;
      }
      if (!isAcceptedJet<aod::JetParticles>(jet)) {
        return;
      }
      if (doMBGapTrigger && jet.eventWeight() == 1) {
        return;
      }

      if (checkMcCollisionIsMatched) {
        auto collisionspermcpjet = collisions.sliceBy(CollisionsPerMCPCollision, jet.mcCollisionId());

        if (collisionspermcpjet.size() >= 1 && jetderiveddatautilities::selectCollision(collisionspermcpjet.begin(), eventSelectionBits)) {
          fillMCPHistograms(jet, jet.eventWeight());
        }
      } else {
        fillMCPHistograms(jet, jet.eventWeight());
      }
    }
  }
  PROCESS_SWITCH(FullJetSpectrapp, processJetsMCPWeighted, "Full Jets at Particle Level on weighted events", false);

  void processJetsMCPMCDMatched(soa::Filtered<EMCCollisionsMCD>::iterator const& collision, JetTableMCDMatchedJoined const& mcdjets, JetTableMCPMatchedJoined const& mcpjets, aod::JMcCollisions const&, aod::JetTracks const&, aod::JetClusters const&, aod::JetParticles const&)
  {
    bool eventAccepted = false;
    int fakemcdjet = 0;
    int fakemcpjet = 0;
    const auto mcpJetsPerMcCollision = mcpjets.sliceBy(JetMCPPerMcCollision, collision.mcCollisionId());

    registry.fill(HIST("h_Matchedcollision_counter"), 0.5);

    if (fabs(collision.posZ()) > VertexZCut) { // making double sure this condition is satisfied
      return;
    }
    registry.fill(HIST("h_Matchedcollision_counter"), 1.5);

    if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      return;
    }
    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        if (collision.alias_bit(kTVXinEMC)) {
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
      }
    }
    if (!eventAccepted) {
      return;
    }
    registry.fill(HIST("h_Matchedcollision_counter"), 2.5);

    for (const auto& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      // Check if MCD jet is within the EMCAL fiducial region; if not then flag it as a fake jet
      if (mcdjet.phi() < jetPhiMin || mcdjet.phi() > jetPhiMax || mcdjet.eta() < jetEtaMin || mcdjet.eta() > jetEtaMax) {
        fakemcdjet++;
        registry.fill(HIST("h2_full_fakemcdjets"), mcdjet.pt(), fakemcdjet, 1.0);
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      for (auto& mcpjet : mcdjet.template matchedJetGeo_as<JetTableMCPMatchedJoined>()) {
        // apply emcal fiducial cuts to the matched particle level jets
        if (mcpjet.eta() > jetEtaMax || mcpjet.eta() < jetEtaMin || mcpjet.phi() > jetPhiMax || mcpjet.phi() < jetPhiMin) {
          fakemcpjet++;
          registry.fill(HIST("h2_full_fakemcpjets"), mcpjet.pt(), fakemcpjet, 1.0);
          continue;
        }
      } // mcpjet loop
      fillMatchedHistograms<JetTableMCDMatchedJoined::iterator, JetTableMCPMatchedJoined>(mcdjet);
    } // mcdjet loop
  }
  PROCESS_SWITCH(FullJetSpectrapp, processJetsMCPMCDMatched, "Full Jet finder MCP matched to MCD", false);

  void processJetsMCPMCDMatchedWeighted(soa::Filtered<EMCCollisionsMCD>::iterator const& collision, JetTableMCDMatchedWeightedJoined const& mcdjets, JetTableMCPMatchedWeightedJoined const& mcpjets, aod::JMcCollisions const&, aod::JetTracks const&, aod::JetClusters const&, aod::JetParticles const&)
  {
    float eventWeight = collision.mcCollision().weight();

    registry.fill(HIST("h_Matchedcollision_counter"), 0.5);
    if (fabs(collision.posZ()) > VertexZCut) { // making double sure this condition is satisfied
      return;
    }
    registry.fill(HIST("h_Matchedcollision_counter"), 1.5);

    if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      return;
    }
    bool eventAccepted = false;
    int fakemcdjet = 0;
    int fakemcpjet = 0;
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    const auto mcpJetsPerMcCollision = mcpjets.sliceBy(JetMCPPerMcCollision, collision.mcCollisionId());

    for (auto mcpjet : mcpJetsPerMcCollision) {
      if (mcpjet.pt() > pTHatMaxMCP * pTHat) { // outlier rejection for MCP
        return;
      }
    }
    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        if (collision.alias_bit(kTVXinEMC)) {
          registry.fill(HIST("h_Matchedcollision_counter"), eventWeight);
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
        registry.fill(HIST("h_Matchedcollision_counter"), eventWeight);
      }
    }
    if (!eventAccepted) {
      return;
    }
    registry.fill(HIST("h_Matchedcollision_counter"), 2.5);

    for (const auto& mcdjet : mcdjets) {
      // Check if MCD jet is within the EMCAL fiducial region; if not then flag it as a fake jet
      if (mcdjet.phi() < jetPhiMin || mcdjet.phi() > jetPhiMax || mcdjet.eta() < jetEtaMin || mcdjet.eta() > jetEtaMax) {
        fakemcdjet++;
        registry.fill(HIST("h2_full_fakemcdjets"), mcdjet.pt(), fakemcdjet, eventWeight);
        continue;
      }
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }

      for (auto& mcpjet : mcdjet.template matchedJetGeo_as<JetTableMCPMatchedWeightedJoined>()) {
        // apply emcal fiducial cuts to the matched particle level jets
        if (mcpjet.eta() > jetEtaMax || mcpjet.eta() < jetEtaMin || mcpjet.phi() > jetPhiMax || mcpjet.phi() < jetPhiMin) {
          fakemcpjet++;
          registry.fill(HIST("h2_full_fakemcpjets"), mcpjet.pt(), fakemcpjet, eventWeight);
          continue;
        }
        // // If both MCD-MCP matched jet pairs are within the EMCAL fiducial region, fill these histos
        registry.fill(HIST("h_full_matchedmcpjet_eta"), mcpjet.eta(), eventWeight);
        registry.fill(HIST("h_full_matchedmcpjet_phi"), mcpjet.phi(), eventWeight);
      } // mcpjet
      fillMatchedHistograms<JetTableMCDMatchedWeightedJoined::iterator, JetTableMCPMatchedWeightedJoined>(mcdjet, eventWeight);
    } // mcdjet
  }
  PROCESS_SWITCH(FullJetSpectrapp, processJetsMCPMCDMatchedWeighted, "Full Jet finder MCP matched to MCD on weighted events", false);

  void processTracks(soa::Filtered<EMCCollisionsMCD>::iterator const& collision, soa::Filtered<aod::JetTracks> const& tracks, soa::Filtered<aod::JetClusters> const& clusters)
  {
    bool eventAccepted = false;

    registry.fill(HIST("h_collisions_unweighted"), 0.5);
    if (fabs(collision.posZ()) > VertexZCut) {
      return;
    }
    registry.fill(HIST("h_collisions_unweighted"), 1.5);
    if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      return;
    }
    // needed for the workaround to access EMCAL trigger bits. - This is needed for the MC productions in which the EMC trigger bits are missing. (MB MC LHC24f3, for ex.)
    // It first requires for atleast a cell in EMCAL to have energy content.
    // Once it finds a cell content,
    // it then checks if the collision is not an ambiguous collision (i.e. it has to be a unique collision = no bunch pile up)
    // If all of these conditions are satisfied then it checks for the required trigger bit in EMCAL.
    // For LHC22o, since the EMCAL didn't have hardware triggers, one would only require MB trigger (kTVXinEMC) in the EMCAL.

    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        if (collision.alias_bit(kTVXinEMC)) {
          registry.fill(HIST("h_collisions_unweighted"), 4.0); // Tracks with kTVXinEMC
          registry.fill(HIST("h_Detcollision_counter"), 1.0);
        }
      }
    } else {
      // Check if EMCAL was readout with the MB trigger(kTVXinEMC) fired. If not then reject the event and exit the function.
      // This is the default check for the simulations with proper trigger flags not requiring the above workaround.
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
        registry.fill(HIST("h_collisions_unweighted"), 4.0); // Tracks with kTVXinEMC
        registry.fill(HIST("h_Detcollision_counter"), 1.0);
      }
    }

    if (!eventAccepted) {
      registry.fill(HIST("h_collisions_unweighted"), 8.0); // Tracks w/o kTVXinEMC
      return;
    }
    registry.fill(HIST("h_collisions_unweighted"), 2.5);

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      // Fill Accepted events histos
      fillTrackHistograms(tracks, clusters, 1.0);
    }
    registry.fill(HIST("h_collisions_unweighted"), 3.5);
  }
  PROCESS_SWITCH(FullJetSpectrapp, processTracks, "Full Jet tracks", false);

  void processTracksWeighted(soa::Filtered<EMCCollisionsMCD>::iterator const& collision,
                             aod::JMcCollisions const&,
                             soa::Filtered<aod::JetTracks> const& tracks,
                             soa::Filtered<aod::JetClusters> const& clusters)
  {
    bool eventAccepted = false;
    float eventWeight = collision.mcCollision().weight();

    registry.fill(HIST("h_collisions_weighted"), 0.5);
    if (fabs(collision.posZ()) > VertexZCut) {
      return;
    }
    registry.fill(HIST("h_collisions_weighted"), 1.5);
    if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      return;
    }
    if (doMBGapTrigger && eventWeight == 1) {
      return;
    }
    registry.fill(HIST("h_collisions_weighted"), 1.0, eventWeight); // total events

    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        fillTrackHistograms(tracks, clusters, eventWeight);
        if (collision.alias_bit(kTVXinEMC)) {
          registry.fill(HIST("h_collisions_weighted"), 4.0, eventWeight); // TracksWeighted with kTVXinEMC
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
        registry.fill(HIST("h_collisions_weighted"), 4.0, eventWeight); // TracksWeighted with kTVXinEMC
      }
    }

    if (!eventAccepted) {
      registry.fill(HIST("h_collisions_weighted"), 7.0, eventWeight); // TracksWeighted w/o kTVXinEMC
      return;
    }
    registry.fill(HIST("h_collisions_weighted"), 2.5);

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      // Fill Accepted events histos
      fillTrackHistograms(tracks, clusters, 1.0);
    }
    registry.fill(HIST("h_collisions_weighted"), 3.5);
  }
  PROCESS_SWITCH(FullJetSpectrapp, processTracksWeighted, "Full Jet tracks weighted", false);

  void processCollisionsWeightedWithMultiplicity(soa::Filtered<soa::Join<EMCCollisionsMCD, aod::FT0Mults>>::iterator const& collision, JetTableMCDWeightedJoined const& mcdjets, aod::JMcCollisions const&, soa::Filtered<aod::JetTracks> const& tracks, soa::Filtered<aod::JetClusters> const& clusters)
  {
    bool eventAccepted = false;
    float eventWeight = collision.mcCollision().weight();
    float neutralEnergy = 0.0;

    if (fabs(collision.posZ()) > VertexZCut) {
      return;
    }
    if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      return;
    }
    if (doMBGapTrigger && eventWeight == 1) {
      return;
    }

    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        fillTrackHistograms(tracks, clusters, eventWeight);
        if (collision.alias_bit(kTVXinEMC)) {
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
      }
    }

    if (!eventAccepted) {
      return;
    }
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
    }
    registry.fill(HIST("h_FT0Mults_occupancy"), collision.multiplicity());

    for (auto const& mcdjet : mcdjets) {
      float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
      if (mcdjet.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin) { // MCD jets outlier rejection
        return;
      }
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (mcdjet.phi() < jetPhiMin || mcdjet.phi() > jetPhiMax) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      registry.fill(HIST("h2_full_jet_jetpTDetVsFT0Mults"), mcdjet.pt(), collision.multiplicity(), eventWeight);

      for (auto const& cluster : clusters) {
        neutralEnergy += cluster.energy();
      }
      auto NEF = neutralEnergy / mcdjet.energy();
      registry.fill(HIST("h3_full_jet_jetpTDet_FT0Mults_NEF"), mcdjet.pt(), collision.multiplicity(), NEF, eventWeight);
    }
  }
  PROCESS_SWITCH(FullJetSpectrapp, processCollisionsWeightedWithMultiplicity, "Weighted Collisions for Full Jets Multiplicity Studies", false);

  void processMBCollisionsWithMultiplicity(soa::Filtered<soa::Join<EMCCollisionsMCD, aod::FT0Mults>>::iterator const& collision, JetTableMCDJoined const& mcdjets, aod::JMcCollisions const&, soa::Filtered<aod::JetTracks> const& tracks, soa::Filtered<aod::JetClusters> const& clusters)
  {
    bool eventAccepted = false;
    float pTHat = 10. / (std::pow(1.0, 1.0 / pTHatExponent));
    float neutralEnergy = 0.0;

    if (fabs(collision.posZ()) > VertexZCut) {
      return;
    }
    if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      return;
    }

    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        fillTrackHistograms(tracks, clusters, 1.0);
        if (collision.alias_bit(kTVXinEMC)) {
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
      }
    }

    if (!eventAccepted) {
      return;
    }
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
    }
    registry.fill(HIST("h_FT0Mults_occupancy"), collision.multiplicity());

    for (auto const& mcdjet : mcdjets) {
      if (mcdjet.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin) { // MCD (Detector Level) Outlier Rejection
        return;
      }
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (mcdjet.phi() < jetPhiMin || mcdjet.phi() > jetPhiMax) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      registry.fill(HIST("h2_full_jet_jetpTDetVsFT0Mults"), mcdjet.pt(), collision.multiplicity(), 1.0);

      for (auto const& cluster : clusters) {
        neutralEnergy += cluster.energy();
      }
      auto NEF = neutralEnergy / mcdjet.energy();
      registry.fill(HIST("h3_full_jet_jetpTDet_FT0Mults_NEF"), mcdjet.pt(), collision.multiplicity(), NEF, 1.0);
    }
  }
  PROCESS_SWITCH(FullJetSpectrapp, processMBCollisionsWithMultiplicity, "MB Collisions for Full Jets Multiplicity Studies", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FullJetSpectrapp>(cfgc, TaskName{"full-jet-spectra-pp"})};
}
