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
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <TRandom3.h>

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
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetHFUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"

#include "EventFiltering/filterTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename JetTableData, typename JetConstituentTableData, typename JetMatchingTableData, typename CandidateTableData, typename JetTableMCD, typename JetConstituentTableMCD, typename JetMatchingTableMCDMCP, typename JetTableMCDWeighted, typename CandidateTableMCD, typename JetTableMCP, typename JetConstituentTableMCP, typename JetMatchingTableMCPMCD, typename JetTableMCPWeighted, typename JetTableDataSub, typename JetConstituentTableDataSub, typename JetMatchingTableDataSub, typename CandidateTableMCP, typename JetTracksDataSub, typename BkgRhoTable>
struct JetFinderHFQATask {

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
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};
  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};
  Configurable<float> randomConeR{"randomConeR", 0.4, "size of random Cone for estimating background fluctuations"};

  HfHelper hfHelper;
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
      registry.add("h3_jet_r_jet_pt_candidate_pt", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_candidate_eta", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_candidate_phi", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_candidate_y", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});y_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_candidate_invmass_jet_pt_candidate_pt", ";#it{m}_{inv, candidate} (GeV/#it{c}^{2}); #it{p}_{T,jet} (GeV/#it{c}) ;#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{500, 0.0, 5.0}, {200, 0.0, 200.0}, {200, 0.0, 200.0}}});
      registry.add("h3_candidatebar_invmass_jet_pt_candidate_pt", ";#it{m}_{inv, candidate bar} (GeV/#it{c}^{2}); #it{p}_{T,jet} (GeV/#it{c}) ;#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{500, 0.0, 5.0}, {200, 0.0, 200.0}, {200, 0.0, 200.0}}});
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
      registry.add("h3_jet_r_jet_pt_candidate_pt_rhoareasubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, -200., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_candidate_eta_rhoareasubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, -200., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_candidate_phi_rhoareasubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, -200., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_candidate_y_rhoareasubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});y_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, -200., 200.}, {100, -1.0, 1.0}}});
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
      registry.add("h3_jet_r_jet_pt_candidate_pt_eventwiseconstituentsubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_candidate_eta_eventwiseconstituentsubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_candidate_phi_eventwiseconstituentsubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_candidate_y_eventwiseconstituentsubtracted", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});y_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
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
      registry.add("h3_jet_r_part_jet_pt_part_candidate_pt_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#it{p}_{T,candidate}^{part} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_part_jet_pt_part_candidate_eta_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#eta_{candidate}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_part_jet_pt_part_candidate_phi_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});#varphi{candidate}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_part_jet_pt_part_candidate_y_part", "#it{R}_{jet}^{part};#it{p}_{T,jet}^{part} (GeV/#it{c});y_{candidate}^{part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
    }

    if (doprocessJetsMCPMCDMatched || doprocessJetsMCPMCDMatchedWeighted || doprocessJetsSubMatched) {
      registry.add("h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c});#it{p}_{T,jet}^{base} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_eta_tag_jet_eta_base_matchedgeo", "#it{R}_{jet};#eta_{jet}^{tag};#eta_{jet}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_phi_tag_jet_phi_base_matchedgeo", "#it{R}_{jet};#varphi_{jet}^{tag};#varphi_{jet}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedgeo", "#it{R}_{jet};N_{jet tracks}^{tag};N_{jet tracks}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, -0.5, 199.5}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_candidate_pt_tag_candidate_pt_base_matchedgeo", "#it{R}_{jet};#it{p}_{T,candidate}^{tag} (GeV/#it{c});#it{p}_{T,candidate}^{base} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_candidate_eta_tag_candidate_eta_base_matchedgeo", "#it{R}_{jet};#eta_{candidate}^{tag};#eta_{candidate}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_candidate_phi_tag_candidate_phi_base_matchedgeo", "#it{R}_{jet};#varphi_{candidate}^{tag};#varphi_{candidate}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {160, -1.0, 7.}, {160, -1.0, 7.}}});
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
      registry.add("h3_jet_r_candidate_pt_tag_candidate_pt_base_matchedpt", "#it{R}_{jet};#it{p}_{T,candidate}^{tag} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_candidate_eta_tag_candidate_eta_base_matchedpt", "#it{R}_{jet};#eta_{candidate}^{tag};#eta_{candidate}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_candidate_phi_tag_candidate_phi_base_matchedpt", "#it{R}_{jet};#varphi_{candidate}^{tag};#varphi_{candidate}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {160, -1.0, 7.}, {160, -1.0, 7.}}});
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
      registry.add("h3_jet_r_candidate_pt_tag_candidate_pt_base_matchedgeopt", "#it{R}_{jet};#it{p}_{T,candidate}^{tag} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_candidate_eta_tag_candidate_eta_base_matchedgeopt", "#it{R}_{jet};#eta_{candidate}^{tag};#eta_{candidate}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_candidate_phi_tag_candidate_phi_base_matchedgeopt", "#it{R}_{jet};#varphi_{candidate}^{tag};#varphi_{candidate}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {160, -1.0, 7.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedgeopt", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_eta_base_diff_matchedgeopt", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#eta_{jet}^{tag} - #eta_{jet}^{base}) / #eta_{jet}^{tag}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_phi_base_diff_matchedgeopt", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#varphi_{jet}^{tag} - #varphi_{jet}^{base}) / #varphi_{jet}^{tag}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_pt_tag_jet_eta_tag_jet_eta_base_matchedgeopt", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); #eta_{jet}^{tag}; #eta_{jet}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_pt_tag_jet_phi_tag_jet_phi_base_matchedgeopt", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); #varphi_{jet}^{tag}; #varphi_{jet}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {160, -1.0, 7.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_pt_tag_jet_ntracks_tag_jet_ntracks_base_matchedgeopt", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); N_{jet tracks}^{tag}; N_{jet tracks}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {200, -0.5, 199.5}, {200, -0.5, 199.5}}});

      registry.add("h3_jet_r_jet_pt_tag_jet_pt_base_matchedhf", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c});#it{p}_{T,jet}^{base} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_eta_tag_jet_eta_base_matchedhf", "#it{R}_{jet};#eta_{jet}^{tag};#eta_{jet}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_phi_tag_jet_phi_base_matchedhf", "#it{R}_{jet};#varphi_{jet}^{tag};#varphi_{jet}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedhf", "#it{R}_{jet};N_{jet tracks}^{tag};N_{jet tracks}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, -0.5, 199.5}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_candidate_pt_tag_candidate_pt_base_matchedhf", "#it{R}_{jet};#it{p}_{T,candidate}^{tag} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_candidate_eta_tag_candidate_eta_base_matchedhf", "#it{R}_{jet};#eta_{candidate}^{tag};#eta_{candidate}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_candidate_phi_tag_candidate_phi_base_matchedhf", "#it{R}_{jet};#varphi_{candidate}^{tag};#varphi_{candidate}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {160, -1.0, 7.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedhf", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_eta_base_diff_matchedhf", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#eta_{jet}^{tag} - #eta_{jet}^{base}) / #eta_{jet}^{tag}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_phi_base_diff_matchedhf", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#varphi_{jet}^{tag} - #varphi_{jet}^{base}) / #varphi_{jet}^{tag}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_pt_tag_jet_eta_tag_jet_eta_base_matchedhf", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); #eta_{jet}^{tag}; #eta_{jet}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_pt_tag_jet_phi_tag_jet_phi_base_matchedhf", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); #varphi_{jet}^{tag}; #varphi_{jet}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {160, -1.0, 7.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_pt_tag_jet_ntracks_tag_jet_ntracks_base_matchedhf", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); N_{jet tracks}^{tag}; N_{jet tracks}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {200, -0.5, 199.5}, {200, -0.5, 199.5}}});

      registry.add("h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeohf", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c});#it{p}_{T,jet}^{base} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_eta_tag_jet_eta_base_matchedgeohf", "#it{R}_{jet};#eta_{jet}^{tag};#eta_{jet}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_phi_tag_jet_phi_base_matchedgeohf", "#it{R}_{jet};#varphi_{jet}^{tag};#varphi_{jet}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedgeohf", "#it{R}_{jet};N_{jet tracks}^{tag};N_{jet tracks}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, -0.5, 199.5}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_candidate_pt_tag_candidate_pt_base_matchedgeohf", "#it{R}_{jet};#it{p}_{T,candidate}^{tag} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_candidate_eta_tag_candidate_eta_base_matchedgeohf", "#it{R}_{jet};#eta_{candidate}^{tag};#eta_{candidate}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_candidate_phi_tag_candidate_phi_base_matchedgeohf", "#it{R}_{jet};#varphi_{candidate}^{tag};#varphi_{candidate}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {160, -1.0, 7.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedgeohf", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_eta_base_diff_matchedgeohf", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#eta_{jet}^{tag} - #eta_{jet}^{base}) / #eta_{jet}^{tag}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_phi_base_diff_matchedgeohf", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#varphi_{jet}^{tag} - #varphi_{jet}^{base}) / #varphi_{jet}^{tag}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_pt_tag_jet_eta_tag_jet_eta_base_matchedgeohf", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); #eta_{jet}^{tag}; #eta_{jet}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_pt_tag_jet_phi_tag_jet_phi_base_matchedgeohf", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); #varphi_{jet}^{tag}; #varphi_{jet}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {160, -1.0, 7.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_pt_tag_jet_ntracks_tag_jet_ntracks_base_matchedgeohf", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); N_{jet tracks}^{tag}; N_{jet tracks}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {200, -0.5, 199.5}, {200, -0.5, 199.5}}});

      registry.add("h3_jet_r_jet_pt_tag_jet_pt_base_matchedpthf", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c});#it{p}_{T,jet}^{base} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_eta_tag_jet_eta_base_matchedpthf", "#it{R}_{jet};#eta_{jet}^{tag};#eta_{jet}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_phi_tag_jet_phi_base_matchedpthf", "#it{R}_{jet};#varphi_{jet}^{tag};#varphi_{jet}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedpthf", "#it{R}_{jet};N_{jet tracks}^{tag};N_{jet tracks}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, -0.5, 199.5}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_candidate_pt_tag_candidate_pt_base_matchedpthf", "#it{R}_{jet};#it{p}_{T,candidate}^{tag} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_candidate_eta_tag_candidate_eta_base_matchedpthf", "#it{R}_{jet};#eta_{candidate}^{tag};#eta_{candidate}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_candidate_phi_tag_candidate_phi_base_matchedpthf", "#it{R}_{jet};#varphi_{candidate}^{tag};#varphi_{candidate}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {160, -1.0, 7.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedpthf", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_eta_base_diff_matchedpthf", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#eta_{jet}^{tag} - #eta_{jet}^{base}) / #eta_{jet}^{tag}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_phi_base_diff_matchedpthf", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#varphi_{jet}^{tag} - #varphi_{jet}^{base}) / #varphi_{jet}^{tag}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_pt_tag_jet_eta_tag_jet_eta_base_matchedpthf", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); #eta_{jet}^{tag}; #eta_{jet}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_pt_tag_jet_phi_tag_jet_phi_base_matchedpthf", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); #varphi_{jet}^{tag}; #varphi_{jet}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {160, -1.0, 7.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_pt_tag_jet_ntracks_tag_jet_ntracks_base_matchedpthf", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); N_{jet tracks}^{tag}; N_{jet tracks}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {200, -0.5, 199.5}, {200, -0.5, 199.5}}});

      registry.add("h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeopthf", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c});#it{p}_{T,jet}^{base} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_eta_tag_jet_eta_base_matchedgeopthf", "#it{R}_{jet};#eta_{jet}^{tag};#eta_{jet}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_phi_tag_jet_phi_base_matchedgeopthf", "#it{R}_{jet};#varphi_{jet}^{tag};#varphi_{jet}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedgeopthf", "#it{R}_{jet};N_{jet tracks}^{tag};N_{jet tracks}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, -0.5, 199.5}, {200, -0.5, 199.5}}});
      registry.add("h3_jet_r_candidate_pt_tag_candidate_pt_base_matchedgeopthf", "#it{R}_{jet};#it{p}_{T,candidate}^{tag} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_candidate_eta_tag_candidate_eta_base_matchedgeopthf", "#it{R}_{jet};#eta_{candidate}^{tag};#eta_{candidate}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_candidate_phi_tag_candidate_phi_base_matchedgeopthf", "#it{R}_{jet};#varphi_{candidate}^{tag};#varphi_{candidate}^{base}", {HistType::kTH3F, {{jetRadiiBins, ""}, {160, -1.0, 7.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedgeopthf", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_eta_base_diff_matchedgeopthf", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#eta_{jet}^{tag} - #eta_{jet}^{base}) / #eta_{jet}^{tag}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_r_jet_pt_tag_jet_phi_base_diff_matchedgeopthf", "#it{R}_{jet};#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#varphi_{jet}^{tag} - #varphi_{jet}^{base}) / #varphi_{jet}^{tag}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0.0, 200}, {1000, -5.0, 5.0}}});
      registry.add("h3_jet_pt_tag_jet_eta_tag_jet_eta_base_matchedgeopthf", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); #eta_{jet}^{tag}; #eta_{jet}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {100, -1.0, 1.0}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_pt_tag_jet_phi_tag_jet_phi_base_matchedgeopthf", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); #varphi_{jet}^{tag}; #varphi_{jet}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {160, -1.0, 7.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_pt_tag_jet_ntracks_tag_jet_ntracks_base_matchedgeopthf", ";#it{p}_{T,jet}^{tag} (GeV/#it{c}); N_{jet tracks}^{tag}; N_{jet tracks}^{base}", {HistType::kTH3F, {{200, 0.0, 200}, {200, -0.5, 199.5}, {200, -0.5, 199.5}}});
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
      registry.add("h3_jet_r_jet_pt_candidate_pt_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_candidate_eta_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_candidate_phi_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_candidate_y_MB", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});y_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_candidate_pt_Triggered_Low", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_candidate_eta_Triggered_Low", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_candidate_phi_Triggered_Low", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_candidate_y_Triggered_Low", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});y_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_candidate_pt_Triggered_High", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_candidate_eta_Triggered_High", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_candidate_phi_Triggered_High", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_candidate_y_Triggered_High", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});y_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_candidate_pt_Triggered_Both", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,candidate} (GeV/#it{c})", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h3_jet_r_jet_pt_candidate_eta_Triggered_Both", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#eta_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
      registry.add("h3_jet_r_jet_pt_candidate_phi_Triggered_Both", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});#varphi_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {160, -1.0, 7.}}});
      registry.add("h3_jet_r_jet_pt_candidate_y_Triggered_Both", "#it{R}_{jet};#it{p}_{T,jet} (GeV/#it{c});y_{candidate}", {HistType::kTH3F, {{jetRadiiBins, ""}, {200, 0., 200.}, {100, -1.0, 1.0}}});
    }

    if (doprocessTracks || doprocessTracksWeighted) {
      registry.add("h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.add("h_centrality_collisions", "centrality vs collisions; centrality, collisions", {HistType::kTH2F, {{1200, -10.0, 110.0}, {4, 0.0, 4.0}}});
      registry.add("h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_track_eta", "track #eta;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_track_phi", "track #varphi;#varphi_{track};entries", {HistType::kTH1F, {{160, -1.0, 7.}}});
      registry.add("h_track_energy", "track energy;Energy GeV ;entries", {HistType::kTH1F, {{100, 0.0, 100.0}}});
      if (doprocessTracksWeighted) {
        registry.add("h_collisions_weighted", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      }
    }

    if (doprocessTracksSub) {

      registry.add("h_track_pt_eventwiseconstituentsubtracted", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_track_eta_eventwiseconstituentsubtracted", "track #eta;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h_track_phi_eventwiseconstituentsubtracted", "track #varphi;#varphi_{track};entries", {HistType::kTH1F, {{160, -1.0, 7.}}});
      registry.add("h_track_energy_eventwiseconstituentsubtracted", "track energy;Energy GeV ;entries", {HistType::kTH1F, {{100, 0.0, 100.0}}});
    }

    if (doprocessMCCollisionsWeighted) {
      AxisSpec weightAxis = {{VARIABLE_WIDTH, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0}, "weights"};
      registry.add("h_collision_eventweight_part", "event weight;event weight;entries", {HistType::kTH1F, {weightAxis}});
    }

    if (doprocessCandidates) {
      registry.add("h_candidate_invmass", ";#it{m}_{inv, candidate} (GeV/#it{c}^{2});counts", {HistType::kTH1F, {{500, 0.0, 5.0}}});
      registry.add("h_candidate_pt", ";#it{p}_{T,candidate} (GeV/#it{c});counts", {HistType::kTH1F, {{200, 0.0, 200.0}}});
      registry.add("h_candidate_y", ";y_{candidate};counts", {HistType::kTH1F, {{100, -1.0, 1.0}}});
      registry.add("h2_centrality_ncandidates", "centrality vs N_{candidates};centrality;N_{candidates};", {HistType::kTH2F, {{1200, -10.0, 110.0}, {100, 0.0, 100.0}}});
    }
  }

  using JetTableDataJoined = soa::Join<JetTableData, JetConstituentTableData>;
  using JetTableDataMatchedJoined = soa::Join<JetTableData, JetConstituentTableData, JetMatchingTableData>;
  using JetTableDataSubJoined = soa::Join<JetTableDataSub, JetConstituentTableDataSub>;
  using JetTableDataSubMatchedJoined = soa::Join<JetTableDataSub, JetConstituentTableDataSub, JetMatchingTableDataSub>;

  using JetTableMCDJoined = soa::Join<JetTableMCD, JetConstituentTableMCD>;
  using JetTableMCDWeightedJoined = soa::Join<JetTableMCD, JetConstituentTableMCD, JetTableMCDWeighted>;
  using JetTableMCDMatchedJoined = soa::Join<JetTableMCD, JetConstituentTableMCD, JetMatchingTableMCDMCP>;
  using JetTableMCDMatchedWeightedJoined = soa::Join<JetTableMCD, JetConstituentTableMCD, JetMatchingTableMCDMCP, JetTableMCDWeighted>;

  using JetTableMCPJoined = soa::Join<JetTableMCP, JetConstituentTableMCP>;
  using JetTableMCPWeightedJoined = soa::Join<JetTableMCP, JetConstituentTableMCP, JetTableMCPWeighted>;
  using JetTableMCPMatchedJoined = soa::Join<JetTableMCP, JetConstituentTableMCP, JetMatchingTableMCPMCD>;
  using JetTableMCPMatchedWeightedJoined = soa::Join<JetTableMCD, JetConstituentTableMCP, JetMatchingTableMCPMCD, JetTableMCPWeighted>;

  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter trackSubCuts = (aod::jtracksub::pt >= trackPtMin && aod::jtracksub::pt < trackPtMax && aod::jtracksub::eta > trackEtaMin && aod::jtracksub::eta < trackEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);

  // Filter candidateCutsD0 = (aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar);
  // Filter candidateCutsLc = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLcToPKPi || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLcToPiPK);
  // Filter candidateCutsBplus = (aod::hf_sel_candidate_bplus::isSelBplusToD0Pi >= selectionFlagBplus);

  PresliceOptional<soa::Filtered<JetTracksDataSub>> perD0CandidateTracks = aod::bkgd0::candidateId;
  PresliceOptional<soa::Filtered<JetTracksDataSub>> perLcCandidateTracks = aod::bkglc::candidateId;
  PresliceOptional<soa::Filtered<JetTracksDataSub>> perBplusCandidateTracks = aod::bkgbplus::candidateId;

  PresliceOptional<BkgRhoTable> perD0CandidateRhos = aod::bkgd0::candidateId;
  PresliceOptional<BkgRhoTable> perLcCandidateRhos = aod::bkglc::candidateId;
  PresliceOptional<BkgRhoTable> perBplusCandidateRhos = aod::bkgbplus::candidateId;

  template <typename T, typename U, typename V>
  bool isAcceptedJet(V const& jet)
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
      for (auto& hfconstituent : jet.template hfcandidates_as<U>()) {
        if (hfconstituent.pt() >= leadingConstituentPtMin) {
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

  template <typename T, typename U>
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
    registry.fill(HIST("h3_jet_r_jet_pt_jet_ntracks"), jet.r() / 100.0, jet.pt(), jet.tracksIds().size() + jet.hfcandidatesIds().size(), weight);
    registry.fill(HIST("h3_jet_r_jet_pt_jet_area"), jet.r() / 100.0, jet.pt(), jet.area(), weight);

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
      registry.fill(HIST("h3_jet_r_jet_pt_candidate_y"), jet.r() / 100.0, jet.pt(), hfcandidate.y(), weight);

      if (jet.r() == round(selectedJetsRadius * 100.0f)) {
        registry.fill(HIST("h3_candidate_invmass_jet_pt_candidate_pt"), hfcandidate.m(), jet.pt(), hfcandidate.pt(), weight);
        registry.fill(HIST("h3_candidatebar_invmass_jet_pt_candidate_pt"), hfcandidate.m(), jet.pt(), hfcandidate.pt(), weight);
      }
    }
  }

  template <typename T, typename U>
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

    for (auto& hfcandidate : jet.template hfcandidates_as<std::decay_t<U>>()) {
      registry.fill(HIST("h3_jet_r_jet_pt_candidate_pt_rhoareasubtracted"), jet.r() / 100.0, jet.pt() - (rho * jet.area()), hfcandidate.pt(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_candidate_eta_rhoareasubtracted"), jet.r() / 100.0, jet.pt() - (rho * jet.area()), hfcandidate.eta(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_candidate_phi_rhoareasubtracted"), jet.r() / 100.0, jet.pt() - (rho * jet.area()), hfcandidate.phi(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_candidate_y_rhoareasubtracted"), jet.r() / 100.0, jet.pt() - (rho * jet.area()), hfcandidate.y(), weight);
    }
  }

  template <typename T, typename U, typename V>
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

    for (auto& constituent : jet.template tracks_as<V>()) {

      registry.fill(HIST("h3_jet_r_jet_pt_track_pt_eventwiseconstituentsubtracted"), jet.r() / 100.0, jet.pt(), constituent.pt(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_track_eta_eventwiseconstituentsubtracted"), jet.r() / 100.0, jet.pt(), constituent.eta(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_track_phi_eventwiseconstituentsubtracted"), jet.r() / 100.0, jet.pt(), constituent.phi(), weight);
    }

    for (auto& hfcandidate : jet.template hfcandidates_as<std::decay_t<U>>()) {
      registry.fill(HIST("h3_jet_r_jet_pt_candidate_pt_eventwiseconstituentsubtracted"), jet.r() / 100.0, jet.pt(), hfcandidate.pt(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_candidate_eta_eventwiseconstituentsubtracted"), jet.r() / 100.0, jet.pt(), hfcandidate.eta(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_candidate_phi_eventwiseconstituentsubtracted"), jet.r() / 100.0, jet.pt(), hfcandidate.phi(), weight);
      registry.fill(HIST("h3_jet_r_jet_pt_candidate_y_eventwiseconstituentsubtracted"), jet.r() / 100.0, jet.pt(), hfcandidate.y(), weight);
    }
  }

  template <typename T, typename U, typename V>
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
      registry.fill(HIST("h_jet_ntracks_part"), jet.tracksIds().size() + jet.hfcandidatesIds().size(), weight);
    }

    registry.fill(HIST("h3_jet_r_part_jet_pt_part_jet_eta_part"), jet.r() / 100.0, jet.pt(), jet.eta(), weight);
    registry.fill(HIST("h3_jet_r_part_jet_pt_part_jet_phi_part"), jet.r() / 100.0, jet.pt(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_part_jet_eta_part_jet_phi_part"), jet.r() / 100.0, jet.eta(), jet.phi(), weight);
    registry.fill(HIST("h3_jet_r_part_jet_pt_part_jet_ntracks_part"), jet.r() / 100.0, jet.pt(), jet.tracksIds().size() + jet.hfcandidatesIds().size(), weight);

    for (auto& constituent : jet.template tracks_as<std::decay_t<U>>()) {

      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_pt_part"), jet.r() / 100.0, jet.pt(), constituent.pt(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_eta_part"), jet.r() / 100.0, jet.pt(), constituent.eta(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_phi_part"), jet.r() / 100.0, jet.pt(), constituent.phi(), weight);
    }

    for (auto& hfcandidate : jet.template hfcandidates_as<std::decay_t<V>>()) {

      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_pt_part"), jet.r() / 100.0, jet.pt(), hfcandidate.pt(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_eta_part"), jet.r() / 100.0, jet.pt(), hfcandidate.eta(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_track_phi_part"), jet.r() / 100.0, jet.pt(), hfcandidate.phi(), weight);

      registry.fill(HIST("h3_jet_r_part_jet_pt_part_candidate_pt_part"), jet.r() / 100.0, jet.pt(), hfcandidate.pt(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_candidate_eta_part"), jet.r() / 100.0, jet.pt(), hfcandidate.eta(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_candidate_phi_part"), jet.r() / 100.0, jet.pt(), hfcandidate.phi(), weight);
      registry.fill(HIST("h3_jet_r_part_jet_pt_part_candidate_y_part"), jet.r() / 100.0, jet.pt(), hfcandidate.y(), weight);
    }
  }

  template <typename T, typename U, typename M, typename N>
  void fillMatchedHistograms(T const& jetBase, float weight = 1.0)
  {

    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jetBase.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    auto candidateBasePt = 0.0;
    auto candidateBasePhi = 0.0;
    auto candidateBaseEta = 0.0;
    auto candidateTagPt = 0.0;
    auto candidateTagPhi = 0.0;
    auto candidateTagEta = 0.0;

    auto candidateBase = jetBase.template hfcandidates_first_as<std::decay_t<M>>();
    candidateBasePt = candidateBase.pt();
    candidateBasePhi = candidateBase.phi();
    candidateBaseEta = candidateBase.eta();

    if (jetBase.has_matchedJetGeo()) {
      for (auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<U>>()) {
        if (jetTag.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }
        auto candidateTag = jetTag.template hfcandidates_first_as<std::decay_t<N>>();
        candidateTagPt = candidateTag.pt();
        candidateTagPhi = candidateTag.phi();
        candidateTagEta = candidateTag.eta();

        registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeo"), jetBase.r() / 100.0, jetTag.pt(), jetBase.pt(), weight);
        registry.fill(HIST("h3_jet_r_jet_eta_tag_jet_eta_base_matchedgeo"), jetBase.r() / 100.0, jetTag.eta(), jetBase.eta(), weight);
        registry.fill(HIST("h3_jet_r_jet_phi_tag_jet_phi_base_matchedgeo"), jetBase.r() / 100.0, jetTag.phi(), jetBase.phi(), weight);
        registry.fill(HIST("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedgeo"), jetBase.r() / 100.0, jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
        registry.fill(HIST("h3_jet_r_candidate_pt_tag_candidate_pt_base_matchedgeo"), jetBase.r() / 100.0, candidateTagPt, candidateBasePt, weight);
        registry.fill(HIST("h3_jet_r_candidate_eta_tag_candidate_eta_base_matchedgeo"), jetBase.r() / 100.0, candidateTagEta, candidateBaseEta, weight);
        registry.fill(HIST("h3_jet_r_candidate_phi_tag_candidate_phi_base_matchedgeo"), jetBase.r() / 100.0, candidateTagPhi, candidateBasePhi, weight);
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
        auto candidateTag = jetTag.template hfcandidates_first_as<std::decay_t<N>>();
        candidateTagPt = candidateTag.pt();
        candidateTagPhi = candidateTag.phi();
        candidateTagEta = candidateTag.eta();

        registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_matchedpt"), jetBase.r() / 100.0, jetTag.pt(), jetBase.pt(), weight);
        registry.fill(HIST("h3_jet_r_jet_eta_tag_jet_eta_base_matchedpt"), jetBase.r() / 100.0, jetTag.eta(), jetBase.eta(), weight);
        registry.fill(HIST("h3_jet_r_jet_phi_tag_jet_phi_base_matchedpt"), jetBase.r() / 100.0, jetTag.phi(), jetBase.phi(), weight);
        registry.fill(HIST("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedpt"), jetBase.r() / 100.0, jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
        registry.fill(HIST("h3_jet_r_candidate_pt_tag_candidate_pt_base_matchedpt"), jetBase.r() / 100.0, candidateTagPt, candidateBasePt, weight);
        registry.fill(HIST("h3_jet_r_candidate_eta_tag_candidate_eta_base_matchedpt"), jetBase.r() / 100.0, candidateTagEta, candidateBaseEta, weight);
        registry.fill(HIST("h3_jet_r_candidate_phi_tag_candidate_phi_base_matchedpt"), jetBase.r() / 100.0, candidateTagPhi, candidateBasePhi, weight);
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
    if (jetBase.has_matchedJetCand()) {
      for (auto& jetTag : jetBase.template matchedJetCand_as<std::decay_t<U>>()) {
        if (jetTag.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }
        auto candidateTag = jetTag.template hfcandidates_first_as<std::decay_t<N>>();
        candidateTagPt = candidateTag.pt();
        candidateTagPhi = candidateTag.phi();
        candidateTagEta = candidateTag.eta();

        registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_matchedhf"), jetBase.r() / 100.0, jetTag.pt(), jetBase.pt(), weight);
        registry.fill(HIST("h3_jet_r_jet_eta_tag_jet_eta_base_matchedhf"), jetBase.r() / 100.0, jetTag.eta(), jetBase.eta(), weight);
        registry.fill(HIST("h3_jet_r_jet_phi_tag_jet_phi_base_matchedhf"), jetBase.r() / 100.0, jetTag.phi(), jetBase.phi(), weight);
        registry.fill(HIST("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedhf"), jetBase.r() / 100.0, jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
        registry.fill(HIST("h3_jet_r_candidate_pt_tag_candidate_pt_base_matchedhf"), jetBase.r() / 100.0, candidateTagPt, candidateBasePt, weight);
        registry.fill(HIST("h3_jet_r_candidate_eta_tag_candidate_eta_base_matchedhf"), jetBase.r() / 100.0, candidateTagEta, candidateBaseEta, weight);
        registry.fill(HIST("h3_jet_r_candidate_phi_tag_candidate_phi_base_matchedhf"), jetBase.r() / 100.0, candidateTagPhi, candidateBasePhi, weight);
        registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedhf"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.pt() - jetBase.pt()) / jetTag.pt(), weight);
        registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_eta_base_diff_matchedhf"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.eta() - jetBase.eta()) / jetTag.eta(), weight);
        registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_phi_base_diff_matchedhf"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.phi() - jetBase.phi()) / jetTag.phi(), weight);

        if (jetBase.r() == round(selectedJetsRadius * 100.0f)) {
          registry.fill(HIST("h3_jet_pt_tag_jet_eta_tag_jet_eta_base_matchedhf"), jetTag.pt(), jetTag.eta(), jetBase.eta(), weight);
          registry.fill(HIST("h3_jet_pt_tag_jet_phi_tag_jet_phi_base_matchedhf"), jetTag.pt(), jetTag.phi(), jetBase.phi(), weight);
          registry.fill(HIST("h3_jet_pt_tag_jet_ntracks_tag_jet_ntracks_base_matchedhf"), jetTag.pt(), jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
        }
      }
    }
    if (jetBase.has_matchedJetGeo() && jetBase.has_matchedJetPt()) {
      for (auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<U>>()) {
        if (jetTag.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }
        if (jetBase.template matchedJetGeo_first_as<std::decay_t<U>>().globalIndex() == jetBase.template matchedJetPt_first_as<std::decay_t<U>>().globalIndex()) { // not a good way to do this

          auto candidateTag = jetTag.template hfcandidates_first_as<std::decay_t<N>>();
          candidateTagPt = candidateTag.pt();
          candidateTagPhi = candidateTag.phi();
          candidateTagEta = candidateTag.eta();

          registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeopt"), jetBase.r() / 100.0, jetTag.pt(), jetBase.pt(), weight);
          registry.fill(HIST("h3_jet_r_jet_eta_tag_jet_eta_base_matchedgeopt"), jetBase.r() / 100.0, jetTag.eta(), jetBase.eta(), weight);
          registry.fill(HIST("h3_jet_r_jet_phi_tag_jet_phi_base_matchedgeopt"), jetBase.r() / 100.0, jetTag.phi(), jetBase.phi(), weight);
          registry.fill(HIST("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedgeopt"), jetBase.r() / 100.0, jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
          registry.fill(HIST("h3_jet_r_candidate_pt_tag_candidate_pt_base_matchedgeopt"), jetBase.r() / 100.0, candidateTagPt, candidateBasePt, weight);
          registry.fill(HIST("h3_jet_r_candidate_eta_tag_candidate_eta_base_matchedgeopt"), jetBase.r() / 100.0, candidateTagEta, candidateBaseEta, weight);
          registry.fill(HIST("h3_jet_r_candidate_phi_tag_candidate_phi_base_matchedgeopt"), jetBase.r() / 100.0, candidateTagPhi, candidateBasePhi, weight);
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
    if (jetBase.has_matchedJetGeo() && jetBase.has_matchedJetCand()) {
      for (auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<U>>()) {
        if (jetTag.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }
        if (jetBase.template matchedJetGeo_first_as<std::decay_t<U>>().globalIndex() == jetBase.template matchedJetCand_first_as<std::decay_t<U>>().globalIndex()) { // not a good way to do this

          auto candidateTag = jetTag.template hfcandidates_first_as<std::decay_t<N>>();
          candidateTagPt = candidateTag.pt();
          candidateTagPhi = candidateTag.phi();
          candidateTagEta = candidateTag.eta();

          registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeohf"), jetBase.r() / 100.0, jetTag.pt(), jetBase.pt(), weight);
          registry.fill(HIST("h3_jet_r_jet_eta_tag_jet_eta_base_matchedgeohf"), jetBase.r() / 100.0, jetTag.eta(), jetBase.eta(), weight);
          registry.fill(HIST("h3_jet_r_jet_phi_tag_jet_phi_base_matchedgeohf"), jetBase.r() / 100.0, jetTag.phi(), jetBase.phi(), weight);
          registry.fill(HIST("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedgeohf"), jetBase.r() / 100.0, jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
          registry.fill(HIST("h3_jet_r_candidate_pt_tag_candidate_pt_base_matchedgeohf"), jetBase.r() / 100.0, candidateTagPt, candidateBasePt, weight);
          registry.fill(HIST("h3_jet_r_candidate_eta_tag_candidate_eta_base_matchedgeohf"), jetBase.r() / 100.0, candidateTagEta, candidateBaseEta, weight);
          registry.fill(HIST("h3_jet_r_candidate_phi_tag_candidate_phi_base_matchedgeohf"), jetBase.r() / 100.0, candidateTagPhi, candidateBasePhi, weight);
          registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedgeohf"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.pt() - jetBase.pt()) / jetTag.pt(), weight);
          registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_eta_base_diff_matchedgeohf"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.eta() - jetBase.eta()) / jetTag.eta(), weight);
          registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_phi_base_diff_matchedgeohf"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.phi() - jetBase.phi()) / jetTag.phi(), weight);

          if (jetBase.r() == round(selectedJetsRadius * 100.0f)) {
            registry.fill(HIST("h3_jet_pt_tag_jet_eta_tag_jet_eta_base_matchedgeohf"), jetTag.pt(), jetTag.eta(), jetBase.eta(), weight);
            registry.fill(HIST("h3_jet_pt_tag_jet_phi_tag_jet_phi_base_matchedgeohf"), jetTag.pt(), jetTag.phi(), jetBase.phi(), weight);
            registry.fill(HIST("h3_jet_pt_tag_jet_ntracks_tag_jet_ntracks_base_matchedgeohf"), jetTag.pt(), jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
          }
        }
      }
    }
    if (jetBase.has_matchedJetPt() && jetBase.has_matchedJetCand()) {
      for (auto& jetTag : jetBase.template matchedJetPt_as<std::decay_t<U>>()) {
        if (jetTag.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }
        if (jetBase.template matchedJetPt_first_as<std::decay_t<U>>().globalIndex() == jetBase.template matchedJetCand_first_as<std::decay_t<U>>().globalIndex()) { // not a good way to do this

          auto candidateTag = jetTag.template hfcandidates_first_as<std::decay_t<N>>();
          candidateTagPt = candidateTag.pt();
          candidateTagPhi = candidateTag.phi();
          candidateTagEta = candidateTag.eta();

          registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_matchedpthf"), jetBase.r() / 100.0, jetTag.pt(), jetBase.pt(), weight);
          registry.fill(HIST("h3_jet_r_jet_eta_tag_jet_eta_base_matchedpthf"), jetBase.r() / 100.0, jetTag.eta(), jetBase.eta(), weight);
          registry.fill(HIST("h3_jet_r_jet_phi_tag_jet_phi_base_matchedpthf"), jetBase.r() / 100.0, jetTag.phi(), jetBase.phi(), weight);
          registry.fill(HIST("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedpthf"), jetBase.r() / 100.0, jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
          registry.fill(HIST("h3_jet_r_candidate_pt_tag_candidate_pt_base_matchedpthf"), jetBase.r() / 100.0, candidateTagPt, candidateBasePt, weight);
          registry.fill(HIST("h3_jet_r_candidate_eta_tag_candidate_eta_base_matchedpthf"), jetBase.r() / 100.0, candidateTagEta, candidateBaseEta, weight);
          registry.fill(HIST("h3_jet_r_candidate_phi_tag_candidate_phi_base_matchedpthf"), jetBase.r() / 100.0, candidateTagPhi, candidateBasePhi, weight);
          registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedpthf"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.pt() - jetBase.pt()) / jetTag.pt(), weight);
          registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_eta_base_diff_matchedpthf"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.eta() - jetBase.eta()) / jetTag.eta(), weight);
          registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_phi_base_diff_matchedpthf"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.phi() - jetBase.phi()) / jetTag.phi(), weight);

          if (jetBase.r() == round(selectedJetsRadius * 100.0f)) {
            registry.fill(HIST("h3_jet_pt_tag_jet_eta_tag_jet_eta_base_matchedpthf"), jetTag.pt(), jetTag.eta(), jetBase.eta(), weight);
            registry.fill(HIST("h3_jet_pt_tag_jet_phi_tag_jet_phi_base_matchedpthf"), jetTag.pt(), jetTag.phi(), jetBase.phi(), weight);
            registry.fill(HIST("h3_jet_pt_tag_jet_ntracks_tag_jet_ntracks_base_matchedpthf"), jetTag.pt(), jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
          }
        }
      }
    }

    if (jetBase.has_matchedJetGeo() && jetBase.has_matchedJetPt() && jetBase.has_matchedJetCand()) {
      for (auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<U>>()) {
        if (jetTag.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }
        if (jetBase.template matchedJetGeo_first_as<std::decay_t<U>>().globalIndex() == jetBase.template matchedJetPt_first_as<std::decay_t<U>>().globalIndex()) {
          if (jetBase.template matchedJetGeo_first_as<std::decay_t<U>>().globalIndex() == jetBase.template matchedJetCand_first_as<std::decay_t<U>>().globalIndex()) { // not a good way to do this

            auto candidateTag = jetTag.template hfcandidates_first_as<std::decay_t<N>>();
            candidateTagPt = candidateTag.pt();
            candidateTagPhi = candidateTag.phi();
            candidateTagEta = candidateTag.eta();

            registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_matchedgeopthf"), jetBase.r() / 100.0, jetTag.pt(), jetBase.pt(), weight);
            registry.fill(HIST("h3_jet_r_jet_eta_tag_jet_eta_base_matchedgeopthf"), jetBase.r() / 100.0, jetTag.eta(), jetBase.eta(), weight);
            registry.fill(HIST("h3_jet_r_jet_phi_tag_jet_phi_base_matchedgeopthf"), jetBase.r() / 100.0, jetTag.phi(), jetBase.phi(), weight);
            registry.fill(HIST("h3_jet_r_jet_ntracks_tag_jet_ntracks_base_matchedgeopthf"), jetBase.r() / 100.0, jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
            registry.fill(HIST("h3_jet_r_candidate_pt_tag_candidate_pt_base_matchedgeopthf"), jetBase.r() / 100.0, candidateTagPt, candidateBasePt, weight);
            registry.fill(HIST("h3_jet_r_candidate_eta_tag_candidate_eta_base_matchedgeopthf"), jetBase.r() / 100.0, candidateTagEta, candidateBaseEta, weight);
            registry.fill(HIST("h3_jet_r_candidate_phi_tag_candidate_phi_base_matchedgeopthf"), jetBase.r() / 100.0, candidateTagPhi, candidateBasePhi, weight);
            registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_pt_base_diff_matchedgeopthf"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.pt() - jetBase.pt()) / jetTag.pt(), weight);
            registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_eta_base_diff_matchedgeopthf"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.eta() - jetBase.eta()) / jetTag.eta(), weight);
            registry.fill(HIST("h3_jet_r_jet_pt_tag_jet_phi_base_diff_matchedgeopthf"), jetBase.r() / 100.0, jetTag.pt(), (jetTag.phi() - jetBase.phi()) / jetTag.phi(), weight);

            if (jetBase.r() == round(selectedJetsRadius * 100.0f)) {
              registry.fill(HIST("h3_jet_pt_tag_jet_eta_tag_jet_eta_base_matchedgeopthf"), jetTag.pt(), jetTag.eta(), jetBase.eta(), weight);
              registry.fill(HIST("h3_jet_pt_tag_jet_phi_tag_jet_phi_base_matchedgeopthf"), jetTag.pt(), jetTag.phi(), jetBase.phi(), weight);
              registry.fill(HIST("h3_jet_pt_tag_jet_ntracks_tag_jet_ntracks_base_matchedgeopthf"), jetTag.pt(), jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
            }
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
      registry.fill(HIST("h_track_energy"), track.energy(), weight);
    }
  }

  void processDummy(JetCollisions const& collision)
  {
  }
  PROCESS_SWITCH(JetFinderHFQATask, processDummy, "dummy task", true);

  void processJetsData(soa::Filtered<JetCollisions>::iterator const& collision, JetTableDataJoined const& jets, CandidateTableData const& candidates, JetTracks const& tracks)
  {
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks, CandidateTableData>(jet)) {
        continue;
      }
      fillHistograms<typename JetTableDataJoined::iterator, CandidateTableData>(jet, collision.centrality());
    }
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsData, "jet finder HF QA data", false);

  void processJetsRhoAreaSubData(soa::Filtered<JetCollisions>::iterator const& collision,
                                 BkgRhoTable const& bkgRhos,
                                 JetTableDataJoined const& jets,
                                 CandidateTableData const& candidates,
                                 JetTracks const& tracks)
  {
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks, CandidateTableData>(jet)) {
        continue;
      }
      auto const jetCandidate = jet.template hfcandidates_first_as<CandidateTableData>();
      auto bkgRho = jethfutilities::slicedPerCandidate(bkgRhos, jetCandidate, perD0CandidateRhos, perLcCandidateRhos, perBplusCandidateRhos).iteratorAt(0);
      fillRhoAreaSubtractedHistograms<typename JetTableDataJoined::iterator, CandidateTableData>(jet, collision.centrality(), bkgRho.rho());
    }
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsRhoAreaSubData, "jet finder HF QA for rho-area subtracted jets", false);

  void processEvtWiseConstSubJetsData(soa::Filtered<JetCollisions>::iterator const& collision, JetTableDataSubJoined const& jets, CandidateTableData const& candidates, JetTracksDataSub const& tracks)
  {
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracksDataSub, CandidateTableData>(jet)) {
        continue;
      }
      fillEventWiseConstituentSubtractedHistograms<typename JetTableDataSubJoined::iterator, CandidateTableData, JetTracksDataSub>(jet, collision.centrality());
    }
  }
  PROCESS_SWITCH(JetFinderHFQATask, processEvtWiseConstSubJetsData, "jet finder HF QA for eventwise constituent-subtracted jets data", false);

  void processJetsSubMatched(soa::Filtered<JetCollisions>::iterator const& collision,
                             JetTableDataMatchedJoined const& jets,
                             JetTableDataSubMatchedJoined const& jetSubs,
                             JetTracks const& tracks, JetTracksDataSub const& tracksSub, CandidateTableData const& candidates)
  {
    for (const auto& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks, CandidateTableData>(jet)) {
        continue;
      }
      fillMatchedHistograms<typename JetTableDataMatchedJoined::iterator, JetTableDataSubMatchedJoined, CandidateTableData, CandidateTableData>(jet);
    }
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsSubMatched, "jet finder HF QA matched unsubtracted and constituent subtracted jets", false);

  void processJetsMCD(soa::Filtered<JetCollisions>::iterator const& collision, JetTableMCDJoined const& jets, CandidateTableMCD const& candidates, JetTracks const& tracks)
  {
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks, CandidateTableMCD>(jet)) {
        continue;
      }
      fillHistograms<typename JetTableMCDJoined::iterator, CandidateTableMCD>(jet, collision.centrality());
    }
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsMCD, "jet finder HF QA mcd", false);

  void processJetsMCDWeighted(soa::Filtered<JetCollisions>::iterator const& collision, JetTableMCDWeightedJoined const& jets, CandidateTableMCD const& candidates, JetTracks const& tracks)
  {
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks, CandidateTableMCD>(jet)) {
        continue;
      }
      fillHistograms<typename JetTableMCDWeightedJoined::iterator, CandidateTableMCD>(jet, collision.centrality(), jet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsMCDWeighted, "jet finder HF QA mcd on weighted events", false);

  void processJetsMCP(typename JetTableMCPJoined::iterator const& jet, JetParticles const& particles, CandidateTableMCP const& candidates)
  {
    if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
      return;
    }
    if (!isAcceptedJet<JetParticles, CandidateTableMCP>(jet)) {
      return;
    }
    fillMCPHistograms<typename JetTableMCPJoined::iterator, JetParticles, CandidateTableMCP>(jet);
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsMCP, "jet finder HF QA mcp", false);

  void processJetsMCPWeighted(typename JetTableMCPWeightedJoined::iterator const& jet, JetParticles const& particles, CandidateTableMCP const& candidates)
  {
    if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
      return;
    }
    if (!isAcceptedJet<JetParticles, CandidateTableMCP>(jet)) {
      return;
    }
    fillMCPHistograms<typename JetTableMCPWeightedJoined::iterator, JetParticles, CandidateTableMCP>(jet, jet.eventWeight());
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsMCPWeighted, "jet finder HF QA mcp on weighted events", false);

  void processJetsMCPMCDMatched(soa::Filtered<JetCollisions>::iterator const& collision,
                                JetTableMCDMatchedJoined const& mcdjets,
                                JetTableMCPMatchedJoined const& mcpjets,
                                CandidateTableMCD const& candidatesMCD,
                                JetTracks const& tracks, JetParticles const& particles,
                                CandidateTableMCP const& candidatesMCP)
  {

    for (const auto& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks, CandidateTableMCD>(mcdjet)) {
        continue;
      }
      fillMatchedHistograms<typename JetTableMCDMatchedJoined::iterator, JetTableMCPMatchedJoined, CandidateTableMCD, CandidateTableMCP>(mcdjet);
    }
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsMCPMCDMatched, "jet finder HF QA matched mcp and mcd", false);

  void processJetsMCPMCDMatchedWeighted(soa::Filtered<JetCollisions>::iterator const& collision,
                                        JetTableMCDMatchedWeightedJoined const& mcdjets,
                                        JetTableMCPMatchedWeightedJoined const& mcpjets,
                                        CandidateTableMCD const& candidatesMCD,
                                        JetTracks const& tracks, JetParticles const& particles,
                                        CandidateTableMCP const& candidatesMCP)
  {

    for (const auto& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<JetTracks, CandidateTableMCD>(mcdjet)) {
        continue;
      }
      fillMatchedHistograms<typename JetTableMCDMatchedWeightedJoined::iterator, JetTableMCPMatchedWeightedJoined, CandidateTableMCD, CandidateTableMCP>(mcdjet, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetFinderHFQATask, processJetsMCPMCDMatchedWeighted, "jet finder HF QA matched mcp and mcd on weighted events", false);

  void processMCCollisionsWeighted(JetMcCollision const& collision)
  {
    registry.fill(HIST("h_collision_eventweight_part"), collision.weight());
  }
  PROCESS_SWITCH(JetFinderHFQATask, processMCCollisionsWeighted, "collision QA for weighted events", false);

  void processTriggeredData(soa::Join<JetCollisions, aod::JChTrigSels>::iterator const& collision,
                            JetTableDataJoined const& jets,
                            CandidateTableData const& candidates,
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
      for (auto& hfcandidate : jet.template hfcandidates_as<CandidateTableData>()) {

        registry.fill(HIST("h3_jet_r_jet_pt_candidate_pt_MB"), jet.r() / 100.0, jet.pt(), hfcandidate.pt());
        registry.fill(HIST("h3_jet_r_jet_pt_candidate_eta_MB"), jet.r() / 100.0, jet.pt(), hfcandidate.eta());
        registry.fill(HIST("h3_jet_r_jet_pt_candidate_phi_MB"), jet.r() / 100.0, jet.pt(), hfcandidate.phi());
        registry.fill(HIST("h3_jet_r_jet_pt_candidate_y_MB"), jet.r() / 100.0, jet.pt(), hfcandidate.y());

        if (jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedLow)) {
          registry.fill(HIST("h3_jet_r_jet_pt_candidate_pt_Triggered_Low"), jet.r() / 100.0, jet.pt(), hfcandidate.pt());
          registry.fill(HIST("h3_jet_r_jet_pt_candidate_eta_Triggered_Low"), jet.r() / 100.0, jet.pt(), hfcandidate.eta());
          registry.fill(HIST("h3_jet_r_jet_pt_candidate_phi_Triggered_Low"), jet.r() / 100.0, jet.pt(), hfcandidate.phi());
          registry.fill(HIST("h3_jet_r_jet_pt_candidate_y_Triggered_Low"), jet.r() / 100.0, jet.pt(), hfcandidate.y());
        }
        if (jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedHigh)) {
          registry.fill(HIST("h3_jet_r_jet_pt_candidate_pt_Triggered_High"), jet.r() / 100.0, jet.pt(), hfcandidate.pt());
          registry.fill(HIST("h3_jet_r_jet_pt_candidate_eta_Triggered_High"), jet.r() / 100.0, jet.pt(), hfcandidate.eta());
          registry.fill(HIST("h3_jet_r_jet_pt_candidate_phi_Triggered_High"), jet.r() / 100.0, jet.pt(), hfcandidate.phi());
          registry.fill(HIST("h3_jet_r_jet_pt_candidate_y_Triggered_High"), jet.r() / 100.0, jet.pt(), hfcandidate.y());
        }
        if (jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedLow) && jetderiveddatautilities::selectChargedTrigger(collision, jetderiveddatautilities::JTrigSelCh::chargedHigh)) {
          registry.fill(HIST("h3_jet_r_jet_pt_candidate_pt_Triggered_Both"), jet.r() / 100.0, jet.pt(), hfcandidate.pt());
          registry.fill(HIST("h3_jet_r_jet_pt_candidate_eta_Triggered_Both"), jet.r() / 100.0, jet.pt(), hfcandidate.eta());
          registry.fill(HIST("h3_jet_r_jet_pt_candidate_phi_Triggered_Both"), jet.r() / 100.0, jet.pt(), hfcandidate.phi());
          registry.fill(HIST("h3_jet_r_jet_pt_candidate_y_Triggered_Both"), jet.r() / 100.0, jet.pt(), hfcandidate.y());
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

  PROCESS_SWITCH(JetFinderHFQATask, processTriggeredData, "QA for charged jet trigger", false);

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
  PROCESS_SWITCH(JetFinderHFQATask, processTracks, "QA for charged tracks", false);

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
  PROCESS_SWITCH(JetFinderHFQATask, processTracksWeighted, "QA for charged tracks weighted", false);

  void processTracksSub(soa::Filtered<JetCollisions>::iterator const& collision,
                        CandidateTableData const& candidates,
                        soa::Filtered<JetTracksDataSub> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    for (auto const& candidate : candidates) {

      for (auto const& track : jethfutilities::slicedPerCandidate(tracks, candidate, perD0CandidateTracks, perLcCandidateTracks, perBplusCandidateTracks)) {
        registry.fill(HIST("h_track_pt_eventwiseconstituentsubtracted"), track.pt());
        registry.fill(HIST("h_track_eta_eventwiseconstituentsubtracted"), track.eta());
        registry.fill(HIST("h_track_phi_eventwiseconstituentsubtracted"), track.phi());
        registry.fill(HIST("h_track_energy_eventwiseconstituentsubtracted"), track.energy());
      }
      break; // currently only fills it for the first candidate in the event (not pT ordered)
    }
  }
  PROCESS_SWITCH(JetFinderHFQATask, processTracksSub, "QA for charged event-wise embedded subtracted tracks", false);

  void processRho(JetCollision const& collision, CandidateTableData const& candidates, BkgRhoTable const& bkgRhos, soa::Filtered<JetTracks> const& tracks)
  {
    for (auto const& candidate : candidates) {
      auto bkgRho = jethfutilities::slicedPerCandidate(bkgRhos, candidate, perD0CandidateRhos, perLcCandidateRhos, perBplusCandidateRhos).iteratorAt(0);
      int nTracks = 0;
      for (auto const& track : tracks) {
        if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
          nTracks++;
        }
      }
      registry.fill(HIST("h2_centrality_ntracks"), collision.centrality(), nTracks);
      registry.fill(HIST("h2_ntracks_rho"), nTracks, bkgRho.rho());
      registry.fill(HIST("h2_ntracks_rhom"), nTracks, bkgRho.rhoM());
      registry.fill(HIST("h2_centrality_rho"), collision.centrality(), bkgRho.rho());
      registry.fill(HIST("h2_centrality_rhom"), collision.centrality(), bkgRho.rhoM());
      break; // currently only fills it for the first candidate in the event (not pT ordered)
    }
  }
  PROCESS_SWITCH(JetFinderHFQATask, processRho, "QA for rho-area subtracted jets", false);

  void processRandomCone(soa::Filtered<JetCollisions>::iterator const& collision, JetTableDataSubJoined const& jets, CandidateTableData const& candidates, BkgRhoTable const& bkgRhos, soa::Filtered<JetTracks> const& tracks)
  {

    for (auto const& candidate : candidates) {
      auto bkgRho = jethfutilities::slicedPerCandidate(bkgRhos, candidate, perD0CandidateRhos, perLcCandidateRhos, perBplusCandidateRhos).iteratorAt(0);
      TRandom3 rand(0);
      float randomConeEta = rand.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
      float randomConePhi = rand.Uniform(0.0, 2 * M_PI);
      float randomConePt = 0;
      for (auto const& track : tracks) {
        if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
          float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, static_cast<float>(-M_PI));
          float dEta = track.eta() - randomConeEta;
          if (TMath::Sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
            randomConePt += track.pt();
          }
        }
      }
      registry.fill(HIST("h2_centrality_rhorandomcone"), collision.centrality(), randomConePt - M_PI * randomConeR * randomConeR * bkgRho.rho());

      float dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, static_cast<float>(-M_PI));
      float dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;

      bool jetWasInCone = false;
      while (TMath::Sqrt(dEtaLeadingJet * dEtaLeadingJet + dPhiLeadingJet * dPhiLeadingJet) < jets.iteratorAt(0).r() / 100.0 + randomConeR) {
        jetWasInCone = true;
        randomConeEta = rand.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
        randomConePhi = rand.Uniform(0.0, 2 * M_PI);
        dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, static_cast<float>(-M_PI));
        dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;
      }
      if (jetWasInCone) {
        randomConePt = 0.0;
        for (auto const& track : tracks) {
          if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
            float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, static_cast<float>(-M_PI));
            float dEta = track.eta() - randomConeEta;
            if (TMath::Sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
              randomConePt += track.pt();
            }
          }
        }
      }
      registry.fill(HIST("h2_centrality_rhorandomconewithoutleadingjet"), collision.centrality(), randomConePt - M_PI * randomConeR * randomConeR * bkgRho.rho());
      break; // currently only fills it for the first candidate in the event (not pT ordered). Jet is pT ordered so results for excluding leading jet might not be as expected
    }
  }
  PROCESS_SWITCH(JetFinderHFQATask, processRandomCone, "QA for random cone estimation of background fluctuations", false);

  void processCandidates(soa::Filtered<JetCollisions>::iterator const& collision, CandidateTableData const& candidates)
  {

    for (auto const& candidate : candidates) {
      registry.fill(HIST("h_candidate_invmass"), candidate.m());
      registry.fill(HIST("h_candidate_pt"), candidate.pt());
      registry.fill(HIST("h_candidate_y"), candidate.y());
    }
    registry.fill(HIST("h2_centrality_ncandidates"), collision.centrality(), candidates.size());
  }
  PROCESS_SWITCH(JetFinderHFQATask, processCandidates, "HF candidate QA", false);
};

using JetFinderD0QATask = JetFinderHFQATask<aod::D0ChargedJets, aod::D0ChargedJetConstituents, aod::D0ChargedJetsMatchedToD0ChargedEventWiseSubtractedJets, CandidatesD0Data, aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets, aod::D0ChargedMCDetectorLevelJetEventWeights, CandidatesD0MCD, aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets, aod::D0ChargedMCParticleLevelJetEventWeights, aod::D0ChargedEventWiseSubtractedJets, aod::D0ChargedEventWiseSubtractedJetConstituents, aod::D0ChargedEventWiseSubtractedJetsMatchedToD0ChargedJets, CandidatesD0MCP, aod::JTrackD0Subs, aod::BkgD0Rhos>;
// using JetFinderLcQATask = JetFinderHFQATask<aod::LcChargedJets, aod::LcChargedJetConstituents, aod::LcChargedJetsMatchedToLcChargedEventWiseSubtractedJets, CandidatesLcData, aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents, aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCParticleLevelJets, aod::LcChargedMCDetectorLevelJetEventWeights, CandidatesLcMCD, aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents, aod::LcChargedMCParticleLevelJetsMatchedToLcChargedMCDetectorLevelJets, aod::LcChargedMCParticleLevelJetEventWeights, aod::LcChargedEventWiseSubtractedJets, aod::LcChargedEventWiseSubtractedJetConstituents, aod::LcChargedEventWiseSubtractedJetsMatchedToLcChargedJets, CandidatesLcMCP, aod::JTrackLcSubs, aod::BkgLcRhos>;
// using JetFinderBplusQATask = JetFinderHFQATask<aod::BplusChargedJets, aod::BplusChargedJetConstituents, aod::BplusChargedJetsMatchedToBplusChargedEventWiseSubtractedJets, CandidatesBplusData, aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents, aod::BplusChargedMCDetectorLevelJetsMatchedToBplusChargedMCParticleLevelJets, aod::BplusChargedMCDetectorLevelJetEventWeights, CandidatesBplusMCD, aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents, aod::BplusChargedMCParticleLevelJetsMatchedToBplusChargedMCDetectorLevelJets, aod::BplusChargedMCParticleLevelJetEventWeights, aod::BplusChargedEventWiseSubtractedJets, aod::BplusChargedEventWiseSubtractedJetConstituents, aod::BplusChargedEventWiseSubtractedJetsMatchedToBplusChargedJets, CandidatesBplusMCP, aod::JTrackBplusSubs, aod::BkgBplusRhos>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetFinderD0QATask>(cfgc,
                                                          SetDefaultProcesses{},
                                                          TaskName{"jet-finder-charged-d0-qa"}));

  // tasks.emplace_back(adaptAnalysisTask<JetFinderBplusQATask>(cfgc,
  //                                                       SetDefaultProcesses{},
  //                                                       TaskName{"jet-finder-charged-bplus-qa"}));

  //  tasks.emplace_back(adaptAnalysisTask<JetFinderLcQATask>(cfgc,
  //                                                        SetDefaultProcesses{},
  //                                                      TaskName{"jet-finder-charged-lc-qa"}));

  return WorkflowSpec{tasks};
}
