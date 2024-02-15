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

// jet finder QA task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <cmath>
#include <TRandom3.h>

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "PWGJE/Core/JetDerivedDataUtilities.h"

#include "EventFiltering/filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetFinderQATask {

  HistogramRegistry registry;

  Configurable<float> selectedJetsRadius{"selectedJetsRadius", 0.4, "resolution parameter for histograms without radius"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};
  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 100.0, "maximum pT acceptance for tracks"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};
  Configurable<float> randomConeR{"randomConeR", 0.4, "size of random Cone for estimating background fluctuations"};

  std::vector<bool> filledJetR_Both;
  std::vector<bool> filledJetR_Low;
  std::vector<bool> filledJetR_High;
  std::vector<double> jetRadiiValues;

  int eventSelection = -1;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    jetRadiiValues = (std::vector<double>)jetRadii;

    for (std::size_t iJetRadius = 0; iJetRadius < jetRadiiValues.size(); iJetRadius++) {
      filledJetR_Both.push_back(0.0);
      filledJetR_Low.push_back(0.0);
      filledJetR_High.push_back(0.0);
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
      registry.add("h2_centrality_jet_pt", "centrality vs #it{p}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{1200, -10.0, 110.0}, {200, 0., 200.}}});
      registry.add("h2_centrality_jet_eta", "centrality vs #eta_{jet}; centrality; #eta_{jet}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {100, -1.0, 1.0}}});
      registry.add("h2_centrality_jet_phi", "centrality vs #varphi_{jet}; centrality; #varphi_{jet}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {160, -1.0, 7.}}});
      registry.add("h2_centrality_jet_ntracks", "centrality vs N_{jet tracks}; centrality; N_{jet tracks}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_jet_pt_centrality", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});centrality", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {1200, -10.0, 110.0}}});
      registry.add("h3_jet_r_jet_pt_jet_eta", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_jet_phi", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_eta_jet_phi", "#it{R}_{jet};#eta_{jet};#varphi_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_jet_ntracks", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});N_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_jet_pt_jet_area", "#it{R}_{jet}; #it{p}_{T,jet} (GeV/#it{c}); #it{area}_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {300, 0., 3.}}});
      registry.add("h3_jet_r_jet_pt_track_pt", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_track_eta", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_track_phi", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_leadingtrack_pt", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c}); #it{p}_{T,leading track} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200.0}, {200, 0.0, 200.0}}});
      registry.add("h_jet_phat", "jet #hat{p};#hat{p} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0, 350}}});
      registry.add("h_jet_ptcut", "p_{T} cut;p_{T,jet} (GeV/#it{c});N;entries", {HistType::kTH2F, {{200, 0, 200}, {20, 0, 5}}});
    }

    if (doprocessJetsRhoAreaSubData) {

      registry.add("h_jet_pt_rhoareasubtracted", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{400, -200., 200.}}});
      registry.add("h_jet_eta_rhoareasubtracted", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_jet_phi_rhoareasubtracted", "jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, -1.0, 7.}}});
      registry.add("h_jet_ntracks_rhoareasubtracted", "jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {{200, -0.5, 199.5}}});
      registry.add("h2_centrality_jet_pt_rhoareasubtracted", "centrality vs #it{p}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{1200, -10.0, 110.0}, {400, -200., 200.}}});
      registry.add("h2_centrality_jet_eta_rhoareasubtracted", "centrality vs #eta_{jet}; centrality; #eta_{jet}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {100, -1.0, 1.0}}});
      registry.add("h2_centrality_jet_phi_rhoareasubtracted", "centrality vs #varphi_{jet}; centrality; #varphi_{jet}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {160, -1.0, 7.}}});
      registry.add("h2_centrality_jet_ntracks_rhoareasubtracted", "centrality vs N_{jet tracks}; centrality; N_{jet tracks}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_jet_pt_centrality_rhoareasubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});centrality", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, -200., 200.}, {1200, -10.0, 110.0}}});
      registry.add("h3_jet_r_jet_pt_jet_eta_rhoareasubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, -200., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_jet_phi_rhoareasubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, -200., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_eta_jet_phi_rhoareasubtracted", "#it{R}_{jet};#eta_{jet};#varphi_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_jet_ntracks_rhoareasubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});N_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, -200., 200.}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_jet_pt_jet_area_rhoareasubtracted", "#it{R}_{jet}; #it{p}_{T,jet} (GeV/#it{c}); #it{area}_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, -200., 200.}, {300, 0., 3.}}});
      registry.add("h3_jet_r_jet_pt_track_pt_rhoareasubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, -200., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_track_eta_rhoareasubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, -200., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_track_phi_rhoareasubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, -200., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_jet_pt_rhoareasubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {400, -200., 200.}}});
    }

    if (doprocessEvtWiseConstSubJetsData) {
      registry.add("h_jet_pt_eventwiseconstituentsubtracted", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_jet_eta_eventwiseconstituentsubtracted", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_jet_phi_eventwiseconstituentsubtracted", "jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, -1.0, 7.}}});
      registry.add("h_jet_ntracks_eventwiseconstituentsubtracted", "jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {{200, -0.5, 199.5}}});
      registry.add("h2_centrality_jet_pt_eventwiseconstituentsubtracted", "centrality vs #it{p}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{1200, -10.0, 110.0}, {200, 0., 200.}}});
      registry.add("h2_centrality_jet_eta_eventwiseconstituentsubtracted", "centrality vs #eta_{jet}; centrality; #eta_{jet}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {100, -1.0, 1.0}}});
      registry.add("h2_centrality_jet_phi_eventwiseconstituentsubtracted", "centrality vs #varphi_{jet}; centrality; #varphi_{jet}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {160, -1.0, 7.}}});
      registry.add("h2_centrality_jet_ntracks_eventwiseconstituentsubtracted", "centrality vs N_{jet tracks}; centrality; N_{jet tracks}", {HistType::kTH2F, {{1200, -10.0, 110.0}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_jet_pt_centrality_eventwiseconstituentsubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});centrality", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {1200, -10.0, 110.0}}});
      registry.add("h3_jet_r_jet_pt_jet_eta_eventwiseconstituentsubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_jet_phi_eventwiseconstituentsubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_eta_jet_phi_eventwiseconstituentsubtracted", "#it{R}_{jet};#eta_{jet};#varphi_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_jet_ntracks_eventwiseconstituentsubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});N_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_jet_pt_jet_area_eventwiseconstituentsubtracted", "#it{R}_{jet}; #it{p}_{T,jet} (GeV/#it{c}); #it{area}_{jet}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {300, 0., 3.}}});
      registry.add("h3_jet_r_jet_pt_track_pt_eventwiseconstituentsubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_track_eta_eventwiseconstituentsubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_track_phi_eventwiseconstituentsubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
    }

    if (doprocessRho) {
      registry.add("h2_centrality_ntracks", "; centrality; N_{tracks};", {HistType::kTH2F, {{1100, 0., 110.0}, {10000, 0.0, 10000.0}}});
      registry.add("h2_ntracks_rho", "; N_{tracks}; #it{rho} (GeV/area);", {HistType::kTH2F, {{10000, 0.0, 10000.0}, {400, 0.0, 400.0}}});
      registry.add("h2_ntracks_rhom", "; N_{tracks}; #it{rho}_{m} (GeV/area);", {HistType::kTH2F, {{10000, 0.0, 10000.0}, {100, 0.0, 100.0}}});
      registry.add("h2_centrality_rho", "; centrality; #it{rho} (GeV/area);", {HistType::kTH2F, {{1100, 0., 110.}, {400, 0., 400.0}}});
      registry.add("h2_centrality_rhom", ";centrality; #it{rho}_{m} (GeV/area)", {HistType::kTH2F, {{1100, 0., 110.}, {100, 0., 100.0}}});
    }

    if (doprocessRandomCone) {
      registry.add("h2_centrality_rhorandomcone", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH2F, {{1100, 0., 110.}, {800, -400.0, 400.0}}});
      registry.add("h2_centrality_rhorandomconewithoutleadingjet", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH2F, {{1100, 0., 110.}, {800, -400.0, 400.0}}});
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
      registry.add("h3_jet_r_part_jet_pt_part_track_pt_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T,jet tracks}^{part} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_part_jet_pt_part_track_eta_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#eta_{jet tracks}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_part_jet_pt_part_track_phi_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#varphi_{jet tracks}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h_jet_phat_part", "jet #hat{p};#hat{p} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0, 350}}});
      registry.add("h_jet_ptcut_part", "p_{T} cut;p_{T,jet}^{part} (GeV/#it{c});N;entries", {HistType::kTH2F, {{200, 0, 200}, {20, 0, 5}}});
    }

    if (doprocessJetsMCPMCDMatched || doprocessJetsMCPMCDMatchedWeighted || doprocessJetsSubMatched) {
      registry.add("h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c});#it{p}_{T,jet}^{base} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_eta_tag_jet_eta_base_matchedgeo", "#it{R}_{jet};#eta_{jet}^{tag};#eta_{jet}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_phi_tag_jet_phi_base_matchedgeo", "#it{R}_{jet};#varphi_{jet}^{tag};#varphi_{jet}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedgeo", "#it{R}_{jet};N_{jet tracks}^{tag};N_{jet tracks}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, -0.5, 199.5}, {200, -0.5, 199.5}}});
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
      registry.add("h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedgeopt", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_eta_base_diff_matchedgeopt", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#eta_{jet}^{tag} - #eta_{jet}^{base}) / #eta_{jet}^{tag}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_phi_base_diff_matchedgeopt", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#varphi_{jet}^{tag} - #varphi_{jet}^{base}) / #varphi_{jet}^{tag}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_pt_tag_jet_eta_tag_jet_eta_base_matchedgeopt", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); #eta_{jet}^{tag}; #eta_{jet}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_pt_tag_jet_phi_tag_jet_phi_base_matchedgeopt", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); #varphi_{jet}^{tag}; #varphi_{jet}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {160, -1.0, 7.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_pt_tag_jet_ntracks_tag_jet_ntracks_base_matchedgeopt", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); N_{jet tracks}^{tag}; N_{jet tracks}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {200, -0.5, 199.5}, {200, -0.5, 199.5}}});
    }

    if (doprocessTriggeredData) {
      registry.add("h_collision_trigger_events", "event status;event status;entries", {HistType::kTH1F, {{6, 0.0, 6.0}}});
      registry.add("h_track_pt_MB", "track pT for MB events;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0.0, 200.0}}});
      registry.add("h_track_eta_MB", "track #eta for MB events;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_track_phi_MB", "track #varphi for MB events;#varphi_{track};entries", {HistType::kTH1F, {{160, -1.0, 7.}}});
      registry.add("h_track_pt_Triggered_Low", "track pT for low #it{p}_{T} Triggered events;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0.0, 200.0}}});
      registry.add("h_track_eta_Triggered_Low", "track #eta for low #it{p}_{T} Triggered events;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_track_phi_Triggered_Low", "track #varphi for low #it{p}_{T} Triggered events;#varphi_{track};entries", {HistType::kTH1F, {{160, -1.0, 7.}}});
      registry.add("h_track_pt_Triggered_High", "track pT for high #it{p}_{T} Triggered events;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0.0, 200.0}}});
      registry.add("h_track_eta_Triggered_High", "track #eta for high #it{p}_{T} Triggered events;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_track_phi_Triggered_High", "track #varphi for high #it{p}_{T} Triggered events;#varphi_{track};entries", {HistType::kTH1F, {{160, -1.0, 7.}}});
      registry.add("h_track_pt_Triggered_Both", "track pT for both #it{p}_{T} Triggered events;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0.0, 200.0}}});
      registry.add("h_track_eta_Triggered_Both", "track #eta for both #it{p}_{T} Triggered events;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_track_phi_Triggered_Both", "track #varphi for both #it{p}_{T} Triggered events;#varphi_{track};entries", {HistType::kTH1F, {{160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_collision", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});collision trigger status", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {4, -0.5, 3.5}}});
      registry.add("h3_jet_r_jet_eta_collision", "#it{R}_{jet};#eta_{jet};collision trigger status", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {4, -0.5, 3.5}}});
      registry.add("h3_jet_r_jet_phi_collision", "#it{R}_{jet};#varphi_{jet};collision trigger status", {HistType::kTH3F, {{jetRadiiBins, ""}, {160, -1.0, 7.}, {4, -0.5, 3.5}}});
      registry.add("h2_jet_r_jet_pT_triggered_Low", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{jetRadiiBins, ""}, {200, 0., 200.}}});
      registry.add("h2_jet_r_jet_pT_triggered_High", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{jetRadiiBins, ""}, {200, 0., 200.}}});
      registry.add("h2_jet_r_jet_pT_triggered_Both", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {{jetRadiiBins, ""}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_track_pt_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_track_eta_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_track_phi_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_track_pt_Triggered_Low", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_track_eta_Triggered_Low", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_track_phi_Triggered_Low", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_track_pt_Triggered_High", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_track_eta_Triggered_High", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_track_phi_Triggered_High", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_track_pt_Triggered_Both", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet tracks} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_track_eta_Triggered_Both", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_track_phi_Triggered_Both", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{jet tracks}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
    }

    if (doprocessTracks || doprocessTracksWeighted) {
      registry.add("h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.add("h_centrality_collisions", "centrality vs collisions; centrality, collisions", {HistType::kTH2F, {{1200, -10.0, 110.0}, {4, 0.0, 4.0}}});
      registry.add("h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_track_eta", "track #eta;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_track_phi", "track #varphi;#varphi_{track};entries", {HistType::kTH1F, {{160, -1.0, 7.}}});
      if (doprocessTracksWeighted) {
        registry.add("h_collisions_weighted", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      }
    }
    if (doprocessTracksSub) {

      registry.add("h_track_pt_eventwiseconstituentsubtracted", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_track_eta_eventwiseconstituentsubtracted", "track #eta;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_track_phi_eventwiseconstituentsubtracted", "track #varphi;#varphi_{track};entries", {HistType::kTH1F, {{160, -1.0, 7.}}});
    }

    if (doprocessMCCollisionsWeighted) {
      AxisSpec weightAxis = {{VARIABLE_WIDTH, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0}, "weights"};
      registry.add("h_collision_eventweight_part", "event weight;event weight;entries", {HistType::kTH1F, {weightAxis}});
    }
  }

  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter trackSubCuts = (aod::jtracksub::pt >= trackPtMin && aod::jtracksub::pt < trackPtMax && aod::jtracksub::eta > trackEtaMin && aod::jtracksub::eta < trackEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);

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
  void fillHistograms(T const& jet, float centrality, float weight = 1.0)
  {

    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }

    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("h_jet_ntracks"), jet.tracksIds().size(), weight);
      registry.fill(HIST("h2_centrality_jet_pt"), centrality, jet.pt(), weight);
      registry.fill(HIST("h2_centrality_jet_eta"), centrality, jet.eta(), weight);
      registry.fill(HIST("h2_centrality_jet_phi"), centrality, jet.phi(), weight);
      registry.fill(HIST("h2_centrality_jet_ntracks"), centrality, jet.tracksIds().size(), weight);
    }

    registry.fill(HIST("h3_jet_r_jet_pt_centrality"), jet.r() / 100.0, jet.pt(), centrality, weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_eta"), jet.r() / 100.0, jet.pt(), jet.eta(), weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_phi"), jet.r() / 100.0, jet.pt(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_jet_eta_jet_phi"), jet.r() / 100.0, jet.eta(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_ntracks"), jet.r() / 100.0, jet.pt(), jet.tracksIds().size(), weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_area"), jet.r() / 100.0, jet.pt(), jet.area(), weight);

    for (auto& constituent : jet.template tracks_as<JetTracks>()) {

      registry.fill(HIST("h3_jet_r_jet_pt_track_pt"), jet.r() / 100.0, jet.pt(), constituent.pt(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_track_eta"), jet.r() / 100.0, jet.pt(), constituent.eta(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_track_phi"), jet.r() / 100.0, jet.pt(), constituent.phi(), weight);
    }
  }

  template <typename T>
  void fillRhoAreaSubtractedHistograms(T const& jet, float centrality, float rho, float weight = 1.0)
  {
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_jet_pt_rhoareasubtracted"), jet.pt() - (rho * jet.area()), weight);
      registry.fill(HIST("h_jet_eta_rhoareasubtracted"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi_rhoareasubtracted"), jet.phi(), weight);
      registry.fill(HIST("h_jet_ntracks_rhoareasubtracted"), jet.tracksIds().size(), weight);
      registry.fill(HIST("h2_centrality_jet_pt_rhoareasubtracted"), centrality, jet.pt() - (rho * jet.area()), weight);
      if (jet.pt() - (rho * jet.area()) > 0) {
        registry.fill(HIST("h2_centrality_jet_eta_rhoareasubtracted"), centrality, jet.eta(), weight);
        registry.fill(HIST("h2_centrality_jet_phi_rhoareasubtracted"), centrality, jet.phi(), weight);
        registry.fill(HIST("h2_centrality_jet_ntracks_rhoareasubtracted"), centrality, jet.tracksIds().size(), weight);
      }
    }

    registry.fill(HIST("h3_jet_r_jet_pt_centrality_rhoareasubtracted"), jet.r() / 100.0, jet.pt() - (rho * jet.area()), centrality, weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_eta_rhoareasubtracted"), jet.r() / 100.0, jet.pt() - (rho * jet.area()), jet.eta(), weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_phi_rhoareasubtracted"), jet.r() / 100.0, jet.pt() - (rho * jet.area()), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_jet_eta_jet_phi_rhoareasubtracted"), jet.r() / 100.0, jet.eta(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_ntracks_rhoareasubtracted"), jet.r() / 100.0, jet.pt() - (rho * jet.area()), jet.tracksIds().size(), weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_area_rhoareasubtracted"), jet.r() / 100.0, jet.pt() - (rho * jet.area()), jet.area(), weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_pt_rhoareasubtracted"), jet.r() / 100.0, jet.pt(), jet.pt() - (rho * jet.area()), weight);

    for (auto& constituent : jet.template tracks_as<JetTracks>()) {

      registry.fill(HIST("h3_jet_r_jet_pt_track_pt_rhoareasubtracted"), jet.r() / 100.0, jet.pt() - (rho * jet.area()), constituent.pt(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_track_eta_rhoareasubtracted"), jet.r() / 100.0, jet.pt() - (rho * jet.area()), constituent.eta(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_track_phi_rhoareasubtracted"), jet.r() / 100.0, jet.pt() - (rho * jet.area()), constituent.phi(), weight);
    }
  }

  template <typename T>
  void fillEventWiseConstituentSubtractedHistograms(T const& jet, float centrality, float weight = 1.0)
  {
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_jet_pt_eventwiseconstituentsubtracted"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta_eventwiseconstituentsubtracted"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi_eventwiseconstituentsubtracted"), jet.phi(), weight);
      registry.fill(HIST("h_jet_ntracks_eventwiseconstituentsubtracted"), jet.tracksIds().size(), weight);
      registry.fill(HIST("h2_centrality_jet_pt_eventwiseconstituentsubtracted"), centrality, jet.pt(), weight);
      registry.fill(HIST("h2_centrality_jet_eta_eventwiseconstituentsubtracted"), centrality, jet.eta(), weight);
      registry.fill(HIST("h2_centrality_jet_phi_eventwiseconstituentsubtracted"), centrality, jet.phi(), weight);
      registry.fill(HIST("h2_centrality_jet_ntracks_eventwiseconstituentsubtracted"), centrality, jet.tracksIds().size(), weight);
    }

    registry.fill(HIST("h3_jet_r_jet_pt_centrality_eventwiseconstituentsubtracted"), jet.r() / 100.0, jet.pt(), centrality, weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_eta_eventwiseconstituentsubtracted"), jet.r() / 100.0, jet.pt(), jet.eta(), weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_phi_eventwiseconstituentsubtracted"), jet.r() / 100.0, jet.pt(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_jet_eta_jet_phi_eventwiseconstituentsubtracted"), jet.r() / 100.0, jet.eta(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_ntracks_eventwiseconstituentsubtracted"), jet.r() / 100.0, jet.pt(), jet.tracksIds().size(), weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_area_eventwiseconstituentsubtracted"), jet.r() / 100.0, jet.pt(), jet.area(), weight);

    for (auto& constituent : jet.template tracks_as<JetTracksSub>()) {

      registry.fill(HIST("h3_jet_r_jet_pt_track_pt_eventwiseconstituentsubtracted"), jet.r() / 100.0, jet.pt(), constituent.pt(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_track_eta_eventwiseconstituentsubtracted"), jet.r() / 100.0, jet.pt(), constituent.eta(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_track_phi_eventwiseconstituentsubtracted"), jet.r() / 100.0, jet.pt(), constituent.phi(), weight);
    }
  }

  template <typename T>
  void fillMCPHistograms(T const& jet, float weight = 1.0)
  {

    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCP * pTHat) {
      return;
    }

    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_jet_pt_part"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta_part"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi_part"), jet.phi(), weight);
      registry.fill(HIST("h_jet_ntracks_part"), jet.tracksIds().size(), weight);
    }

    registry.fill(HIST("h3_jet_r_part_jet_pt_part_jet_eta_part"), jet.r() / 100.0, jet.pt(), jet.eta(), weight);
    registry.fill(HIST("h3_jet_r_part_jet_pt_part_jet_phi_part"), jet.r() / 100.0, jet.pt(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_part_jet_eta_part_jet_phi_part"), jet.r() / 100.0, jet.eta(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_part_jet_pt_part_jet_ntracks_part"), jet.r() / 100.0, jet.pt(), jet.tracksIds().size(), weight);

    for (auto& constituent : jet.template tracks_as<JetParticles>()) {

      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_pt_part"), jet.r() / 100.0, jet.pt(), constituent.pt(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_eta_part"), jet.r() / 100.0, jet.pt(), constituent.eta(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_phi_part"), jet.r() / 100.0, jet.pt(), constituent.phi(), weight);
    }
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
        registry.fill(HIST("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedgeo"), jetBase.r() / 100.0, jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
        registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedgeo"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.pt() - jetBase.pt()) / jetTag.pt(), weight);
        registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_eta_base_diff_matchedgeo"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.eta() - jetBase.eta()) / jetTag.eta(), weight);
        registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_phi_base_diff_matchedgeo"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.phi() - jetBase.phi()) / jetTag.phi(), weight);

        if (jetBase.r() == round(selectedJetsRadius * 100.0f)) {
          registry.fill(HIST("h3_jet_pt_tag_jet_eta_tag_jet_eta_base_matchedgeo"), jetTag.pt(), jetTag.eta(), jetBase.eta(), weight);
          registry.fill(HIST("h3_jet_pt_tag_jet_phi_tag_jet_phi_base_matchedgeo"), jetTag.pt(), jetTag.phi(), jetBase.phi(), weight);
          registry.fill(HIST("h3_jet_pt_tag_jet_ntracks_tag_jet_ntracks_base_matchedgeo"), jetTag.pt(), jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
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
        registry.fill(HIST("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedpt"), jetBase.r() / 100.0, jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
        registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedpt"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.pt() - jetBase.pt()) / jetTag.pt(), weight);
        registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_eta_base_diff_matchedpt"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.eta() - jetBase.eta()) / jetTag.eta(), weight);
        registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_phi_base_diff_matchedpt"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.phi() - jetBase.phi()) / jetTag.phi(), weight);

        if (jetBase.r() == round(selectedJetsRadius * 100.0f)) {
          registry.fill(HIST("h3_jet_pt_tag_jet_eta_tag_jet_eta_base_matchedpt"), jetTag.pt(), jetTag.eta(), jetBase.eta(), weight);
          registry.fill(HIST("h3_jet_pt_tag_jet_phi_tag_jet_phi_base_matchedpt"), jetTag.pt(), jetTag.phi(), jetBase.phi(), weight);
          registry.fill(HIST("h3_jet_pt_tag_jet_ntracks_tag_jet_ntracks_base_matchedpt"), jetTag.pt(), jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
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
          registry.fill(HIST("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedgeopt"), jetBase.r() / 100.0, jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
          registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedgeopt"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.pt() - jetBase.pt()) / jetTag.pt(), weight);
          registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_eta_base_diff_matchedgeopt"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.eta() - jetBase.eta()) / jetTag.eta(), weight);
          registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_phi_base_diff_matchedgeopt"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.phi() - jetBase.phi()) / jetTag.phi(), weight);

          if (jetBase.r() == round(selectedJetsRadius * 100.0f)) {
            registry.fill(HIST("h3_jet_pt_tag_jet_eta_tag_jet_eta_base_matchedgeopt"), jetTag.pt(), jetTag.eta(), jetBase.eta(), weight);
            registry.fill(HIST("h3_jet_pt_tag_jet_phi_tag_jet_phi_base_matchedgeopt"), jetTag.pt(), jetTag.phi(), jetBase.phi(), weight);
            registry.fill(HIST("h3_jet_pt_tag_jet_ntracks_tag_jet_ntracks_base_matchedgeopt"), jetTag.pt(), jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
          }
        }
      }
    }
  }

  template <typename T>
  void fillTrackHistograms(T const& tracks, float weight = 1.0)
  {
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      registry.fill(HIST("h_track_pt"), track.pt(), weight);
      registry.fill(HIST("h_track_eta"), track.eta(), weight);
      registry.fill(HIST("h_track_phi"), track.phi(), weight);
    }
  }

  void processJetsData(soa::Filtered<JetCollisions>::iterator const& collision, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets, JetTracks const& tracks)
  {
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(jet)) {
        continue;
      }
      fillHistograms(jet, collision.centrality());
    }
  }
  PROCESS_SWITCH(JetFinderQATask, processJetsData, "jet finder QA data", false);

  void processJetsRhoAreaSubData(soa::Filtered<soa::Join<JetCollisions, aod::BkgChargedRhos>>::iterator const& collision,
                                 soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                                 JetTracks const& tracks)
  {
    for (auto jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(jet)) {
        continue;
      }
      fillRhoAreaSubtractedHistograms(jet, collision.centrality(), collision.rho());
    }
  }
  PROCESS_SWITCH(JetFinderQATask, processJetsRhoAreaSubData, "jet finder QA for rho-area subtracted jets", false);

  void processEvtWiseConstSubJetsData(soa::Filtered<JetCollisions>::iterator const& collision, soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents> const& jets, JetTracksSub const& tracks)
  {
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracksSub>(jet)) {
        continue;
      }
      fillEventWiseConstituentSubtractedHistograms(jet, collision.centrality());
    }
  }
  PROCESS_SWITCH(JetFinderQATask, processEvtWiseConstSubJetsData, "jet finder QA for eventwise constituent-subtracted jets data", false);

  void processJetsSubMatched(soa::Filtered<JetCollisions>::iterator const& collision,
                             soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::ChargedJetsMatchedToChargedEventWiseSubtractedJets> const& jets,
                             soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents, aod::ChargedEventWiseSubtractedJetsMatchedToChargedJets> const& jetSubs,
                             JetTracks const& tracks, JetTracksSub const& tracksSub)
  {
    for (const auto& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(jet)) {
        continue;
      }
      fillMatchedHistograms<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::ChargedJetsMatchedToChargedEventWiseSubtractedJets>::iterator, soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents, aod::ChargedEventWiseSubtractedJetsMatchedToChargedJets>>(jet);
    }
  }
  PROCESS_SWITCH(JetFinderQATask, processJetsSubMatched, "jet finder QA matched unsubtracted and constituent subtracted jets", false);

  void processJetsMCD(soa::Filtered<JetCollisions>::iterator const& collision, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets, JetTracks const& tracks)
  {
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(jet)) {
        continue;
      }
      fillHistograms(jet, collision.centrality());
    }
  }
  PROCESS_SWITCH(JetFinderQATask, processJetsMCD, "jet finder QA mcd", false);

  void processJetsMCDWeighted(soa::Filtered<JetCollisions>::iterator const& collision, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetEventWeights> const& jets, JetTracks const& tracks)
  {
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(jet)) {
        continue;
      }
      double pTHat = 10. / (std::pow(jet.eventWeight(), 1.0 / 6.0));
      for (int N = 1; N < 21; N++) {
        if (jet.pt() < N * 0.25 * pTHat && jet.r() == round(selectedJetsRadius * 100.0f)) {
          registry.fill(HIST("h_jet_ptcut"), jet.pt(), N * 0.25, jet.eventWeight());
        }
      }
      registry.fill(HIST("h_jet_phat"), pTHat);
      fillHistograms(jet, collision.centrality(), jet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetFinderQATask, processJetsMCDWeighted, "jet finder QA mcd with weighted events", false);

  void processJetsMCP(soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>::iterator const& jet, JetParticles const& particles)
  {
    if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
      return;
    }
    if (!isAcceptedJet<JetParticles>(jet)) {
      return;
    }
    fillMCPHistograms(jet);
  }
  PROCESS_SWITCH(JetFinderQATask, processJetsMCP, "jet finder QA mcp", false);

  void processJetsMCPWeighted(soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetEventWeights>::iterator const& jet, JetParticles const& particles)
  {
    if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
      return;
    }
    if (!isAcceptedJet<JetParticles>(jet)) {
      return;
    }
    double pTHat = 10. / (std::pow(jet.eventWeight(), 1.0 / 6.0));
    for (int N = 1; N < 21; N++) {
      if (jet.pt() < N * 0.25 * pTHat && jet.r() == round(selectedJetsRadius * 100.0f)) {
        registry.fill(HIST("h_jet_ptcut_part"), jet.pt(), N * 0.25, jet.eventWeight());
      }
    }
    registry.fill(HIST("h_jet_phat_part"), pTHat);
    fillMCPHistograms(jet, jet.eventWeight());
  }
  PROCESS_SWITCH(JetFinderQATask, processJetsMCPWeighted, "jet finder QA mcp with weighted events", false);

  void processJetsMCPMCDMatched(soa::Filtered<JetCollisions>::iterator const& collision,
                                soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets> const& mcdjets,
                                soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets> const& mcpjets,
                                JetTracks const& tracks, JetParticles const& particles)
  {
    for (const auto& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(mcdjet)) {
        continue;
      }
      fillMatchedHistograms<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>::iterator, soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>(mcdjet);
    }
  }
  PROCESS_SWITCH(JetFinderQATask, processJetsMCPMCDMatched, "jet finder QA matched mcp and mcd", false);

  void processJetsMCPMCDMatchedWeighted(soa::Filtered<JetCollisions>::iterator const& collision,
                                        soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetEventWeights> const& mcdjets,
                                        soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCParticleLevelJetEventWeights> const& mcpjets,
                                        JetTracks const& tracks, JetParticles const& particles)
  {
    for (const auto& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks>(mcdjet)) {
        continue;
      }
      fillMatchedHistograms<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetEventWeights>::iterator, soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCParticleLevelJetEventWeights>>(mcdjet, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetFinderQATask, processJetsMCPMCDMatchedWeighted, "jet finder QA matched mcp and mcd with weighted events", false);

  void processMCCollisionsWeighted(JetMcCollision const& collision)
  {
    registry.fill(HIST("h_collision_eventweight_part"), collision.weight());
  }
  PROCESS_SWITCH(JetFinderQATask, processMCCollisionsWeighted, "collision QA for weighted events", false);

  void processTriggeredData(soa::Join<JetCollisions, aod::JChTrigSels>::iterator const& collision,
                            soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                            soa::Filtered<JetTracks> const& tracks)
  {
    registry.fill(HIST("h_collision_trigger_events"), 0.5); // all events
    if (collision.posZ() > vertexZCut) {
      return;
    }
    registry.fill(HIST("h_collision_trigger_events"), 1.5); // all events with z vertex cut
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("h_collision_trigger_events"), 2.5); // events with sel8()
    if (jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedLow)) {
      registry.fill(HIST("h_collision_trigger_events"), 3.5); // events with high pT triggered jets
    }
    if (jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedHigh)) {
      registry.fill(HIST("h_collision_trigger_events"), 4.5); // events with high pT triggered jets
    }
    if (jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedLow) && jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedHigh)) {
      registry.fill(HIST("h_collision_trigger_events"), 5.5); // events with high pT triggered jets
    }

    for (std::size_t iJetRadius = 0; iJetRadius < jetRadiiValues.size(); iJetRadius++) {
      filledJetR_Both[iJetRadius] = false;
      filledJetR_Low[iJetRadius] = false;
      filledJetR_High[iJetRadius] = false;
    }
    for (auto& jet : jets) {
      for (std::size_t iJetRadius = 0; iJetRadius < jetRadiiValues.size(); iJetRadius++) {
        if (jet.r() == round(jetRadiiValues[iJetRadius] * 100.0f)) {
          if (jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedLow) && !filledJetR_Low[iJetRadius]) {
            filledJetR_Low[iJetRadius] = true;
            for (double pt = 0.0; pt <= jet.pt(); pt += 1.0) {
              registry.fill(HIST("h2_jet_r_jet_pT_triggered_Low"), jet.r() / 100.0, pt);
            }
          }
          if (jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedHigh) && !filledJetR_High[iJetRadius]) {
            filledJetR_High[iJetRadius] = true;
            for (double pt = 0.0; pt <= jet.pt(); pt += 1.0) {
              registry.fill(HIST("h2_jet_r_jet_pT_triggered_High"), jet.r() / 100.0, pt);
            }
          }
          if (jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedLow) && jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedHigh) && !filledJetR_Both[iJetRadius]) {
            filledJetR_Both[iJetRadius] = true;
            for (double pt = 0.0; pt <= jet.pt(); pt += 1.0) {
              registry.fill(HIST("h2_jet_r_jet_pT_triggered_Both"), jet.r() / 100.0, pt);
            }
          }
        }
      }

      if ((jet.eta() < (trackEtaMin + jet.r() / 100.0)) || (jet.eta() > (trackEtaMax - jet.r() / 100.0))) {
        continue;
      }

      registry.fill(HIST("h3_jet_r_jet_pt_collision"), jet.r() / 100.0, jet.pt(), 0.0);
      registry.fill(HIST("h3_jet_r_jet_eta_collision"), jet.r() / 100.0, jet.eta(), 0.0);
      registry.fill(HIST("h3_jet_r_jet_phi_collision"), jet.r() / 100.0, jet.phi(), 0.0);

      if (jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedLow)) {
        registry.fill(HIST("h3_jet_r_jet_pt_collision"), jet.r() / 100.0, jet.pt(), 1.0);
        registry.fill(HIST("h3_jet_r_jet_eta_collision"), jet.r() / 100.0, jet.eta(), 1.0);
        registry.fill(HIST("h3_jet_r_jet_phi_collision"), jet.r() / 100.0, jet.phi(), 1.0);
      }
      if (jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedHigh)) {
        registry.fill(HIST("h3_jet_r_jet_pt_collision"), jet.r() / 100.0, jet.pt(), 2.0);
        registry.fill(HIST("h3_jet_r_jet_eta_collision"), jet.r() / 100.0, jet.eta(), 2.0);
        registry.fill(HIST("h3_jet_r_jet_phi_collision"), jet.r() / 100.0, jet.phi(), 2.0);
      }
      if (jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedLow) && jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedHigh)) {
        registry.fill(HIST("h3_jet_r_jet_pt_collision"), jet.r() / 100.0, jet.pt(), 3.0);
        registry.fill(HIST("h3_jet_r_jet_eta_collision"), jet.r() / 100.0, jet.eta(), 3.0);
        registry.fill(HIST("h3_jet_r_jet_phi_collision"), jet.r() / 100.0, jet.phi(), 3.0);
      }

      for (auto& constituent : jet.template tracks_as<JetTracks>()) {
        registry.fill(HIST("h3_jet_r_jet_pt_track_pt_MB"), jet.r() / 100.0, jet.pt(), constituent.pt());
        registry.fill(HIST("h3_jet_r_jet_pt_track_eta_MB"), jet.r() / 100.0, jet.pt(), constituent.eta());
        registry.fill(HIST("h3_jet_r_jet_pt_track_phi_MB"), jet.r() / 100.0, jet.pt(), constituent.phi());

        if (jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedLow)) {
          registry.fill(HIST("h3_jet_r_jet_pt_track_pt_Triggered_Low"), jet.r() / 100.0, jet.pt(), constituent.pt());
          registry.fill(HIST("h3_jet_r_jet_pt_track_eta_Triggered_Low"), jet.r() / 100.0, jet.pt(), constituent.eta());
          registry.fill(HIST("h3_jet_r_jet_pt_track_phi_Triggered_Low"), jet.r() / 100.0, jet.pt(), constituent.phi());
        }
        if (jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedHigh)) {
          registry.fill(HIST("h3_jet_r_jet_pt_track_pt_Triggered_High"), jet.r() / 100.0, jet.pt(), constituent.pt());
          registry.fill(HIST("h3_jet_r_jet_pt_track_eta_Triggered_High"), jet.r() / 100.0, jet.pt(), constituent.eta());
          registry.fill(HIST("h3_jet_r_jet_pt_track_phi_Triggered_High"), jet.r() / 100.0, jet.pt(), constituent.phi());
        }
        if (jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedLow) && jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedHigh)) {
          registry.fill(HIST("h3_jet_r_jet_pt_track_pt_Triggered_Both"), jet.r() / 100.0, jet.pt(), constituent.pt());
          registry.fill(HIST("h3_jet_r_jet_pt_track_eta_Triggered_Both"), jet.r() / 100.0, jet.pt(), constituent.eta());
          registry.fill(HIST("h3_jet_r_jet_pt_track_phi_Triggered_Both"), jet.r() / 100.0, jet.pt(), constituent.phi());
        }
      }
    }

    for (auto& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      registry.fill(HIST("h_track_pt_MB"), track.pt());
      registry.fill(HIST("h_track_eta_MB"), track.eta());
      registry.fill(HIST("h_track_phi_MB"), track.phi());

      if (jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedLow)) {
        registry.fill(HIST("h_track_pt_Triggered_Low"), track.pt());
        registry.fill(HIST("h_track_eta_Triggered_Low"), track.eta());
        registry.fill(HIST("h_track_phi_Triggered_Low"), track.phi());
      }
      if (jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedHigh)) {
        registry.fill(HIST("h_track_pt_Triggered_High"), track.pt());
        registry.fill(HIST("h_track_eta_Triggered_High"), track.eta());
        registry.fill(HIST("h_track_phi_Triggered_High"), track.phi());
      }
      if (jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedLow) && jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedHigh)) {
        registry.fill(HIST("h_track_pt_Triggered_Both"), track.pt());
        registry.fill(HIST("h_track_eta_Triggered_Both"), track.eta());
        registry.fill(HIST("h_track_phi_Triggered_Both"), track.phi());
      }
    }
  }
  PROCESS_SWITCH(JetFinderQATask, processTriggeredData, "QA for charged jet trigger", false);

  void processTracks(soa::Filtered<JetCollisions>::iterator const& collision,
                     soa::Filtered<JetTracks> const& tracks)
  {
    registry.fill(HIST("h_collisions"), 0.5);
    registry.fill(HIST("h_centrality_collisions"), collision.centrality(), 0.5);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    registry.fill(HIST("h_centrality_collisions"), collision.centrality(), 1.5);
    fillTrackHistograms(tracks);
  }
  PROCESS_SWITCH(JetFinderQATask, processTracks, "QA for charged tracks", false);

  void processTracksWeighted(soa::Join<JetCollisions, aod::JMcCollisionLbs>::iterator const& collision,
                             JetMcCollisions const& mcCollisions,
                             soa::Filtered<JetTracks> const& tracks)
  {
    float eventWeight = collision.mcCollision().weight();
    registry.fill(HIST("h_collisions"), 0.5);
    registry.fill(HIST("h_collisions_weighted"), 0.5, eventWeight);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    registry.fill(HIST("h_collisions_weighted"), 1.5, eventWeight);
    fillTrackHistograms(tracks, eventWeight);
  }
  PROCESS_SWITCH(JetFinderQATask, processTracksWeighted, "QA for charged tracks weighted", false);

  void processTracksSub(soa::Filtered<JetCollisions>::iterator const& collision,
                        soa::Filtered<JetTracksSub> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    for (auto const& track : tracks) {
      registry.fill(HIST("h_track_pt_eventwiseconstituentsubtracted"), track.pt());
      registry.fill(HIST("h_track_eta_eventwiseconstituentsubtracted"), track.eta());
      registry.fill(HIST("h_track_phi_eventwiseconstituentsubtracted"), track.phi());
    }
  }
  PROCESS_SWITCH(JetFinderQATask, processTracksSub, "QA for charged event-wise embedded subtracted tracks", false);

  void processRho(soa::Filtered<soa::Join<JetCollisions, aod::BkgChargedRhos>>::iterator const& collision, soa::Filtered<JetTracks> const& tracks)
  {
    int nTracks = 0;
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        nTracks++;
      }
    }
    registry.fill(HIST("h2_centrality_ntracks"), collision.centrality(), nTracks);
    registry.fill(HIST("h2_ntracks_rho"), nTracks, collision.rho());
    registry.fill(HIST("h2_ntracks_rhom"), nTracks, collision.rhoM());
    registry.fill(HIST("h2_centrality_rho"), collision.centrality(), collision.rho());
    registry.fill(HIST("h2_centrality_rhom"), collision.centrality(), collision.rhoM());
  }
  PROCESS_SWITCH(JetFinderQATask, processRho, "QA for rho-area subtracted jets", false);

  void processRandomCone(soa::Filtered<soa::Join<JetCollisions, aod::BkgChargedRhos>>::iterator const& collision, soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents> const& jets, JetTracksSub const& tracks)
  {

    TRandom3 rand(0);
    float randomConeEta = rand.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
    float randomConePhi = rand.Uniform(0.0, 2 * M_PI);
    float randomConePt = 0;
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, -M_PI);
        float dEta = track.eta() - randomConeEta;
        if (TMath::Sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
          randomConePt += track.pt();
        }
      }
    }
    registry.fill(HIST("h2_centrality_rhorandomcone"), collision.centrality(), randomConePt - M_PI * randomConeR * randomConeR * collision.rho());

    float dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, -M_PI);
    float dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;

    bool jetWasInCone = false;
    while (TMath::Sqrt(dEtaLeadingJet * dEtaLeadingJet + dPhiLeadingJet * dPhiLeadingJet) < randomConeR) {
      jetWasInCone = true;
      randomConeEta = rand.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
      randomConePhi = rand.Uniform(0.0, 2 * M_PI);
      dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, -M_PI);
      dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;
    }
    if (jetWasInCone) {
      randomConePt = 0.0;
      for (auto const& track : tracks) {
        if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
          float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, -M_PI);
          float dEta = track.eta() - randomConeEta;
          if (TMath::Sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
            randomConePt += track.pt();
          }
        }
      }
    }
    registry.fill(HIST("h2_centrality_rhorandomconewithoutleadingjet"), collision.centrality(), randomConePt - M_PI * randomConeR * randomConeR * collision.rho());
  }
  PROCESS_SWITCH(JetFinderQATask, processRandomCone, "QA for random cone estimation of background fluctuations", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetFinderQATask>(cfgc, TaskName{"jet-finder-charged-qa"})}; }
