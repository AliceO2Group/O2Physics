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

/// \file jetSpectraCharged.cxx
/// \brief Charged-particle jet spectra task
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, Aimeric Landou <aimeric.landou@cern.ch>, Wenhui Feng <wenhui.feng@cern.ch>

#include <cmath>
#include <TRandom3.h>
#include <THn.h>
#include <THnSparse.h>
#include <string>
#include <vector>

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

struct JetSpectraCharged {

  using McParticleCollision = soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos>;
  using ChargedMCDMatchedJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>;
  using ChargedMCPMatchedJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>;
  using ChargedMCDMatchedJetsWeighted = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetEventWeights>;
  using ChargedMCPMatchedJetsWeighted = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCParticleLevelJetEventWeights>;

  HistogramRegistry registry;

  Configurable<float> selectedJetsRadius{"selectedJetsRadius", 0.2, "resolution parameter for histograms without radius"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 100.0, "maximum pT acceptance for tracks"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};
  Configurable<double> jetPtMax{"jetPtMax", 200., "set jet pT bin max"};
  Configurable<float> jetEtaMin{"jetEtaMin", -0.7, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 0.7, "maximum jet pseudorapidity"};
  Configurable<int> nBinsEta{"nBinsEta", 200, "number of bins for eta axes"};
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};
  Configurable<float> leadingConstituentPtMax{"leadingConstituentPtMax", 9999.0, "maximum pT selection on jet constituent"};
  Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum track occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  Configurable<int> trackOccupancyInTimeRangeMin{"trackOccupancyInTimeRangeMin", -999999, "minimum track occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  Configurable<bool> checkGeoMatched{"checkGeoMatched", true, "0: turn off geometry matching, 1: do geometry matching "};
  Configurable<bool> checkPtMatched{"checkPtMatched", false, "0: turn off pT matching, 1: do pT matching"};
  Configurable<bool> checkGeoPtMatched{"checkGeoPtMatched", false, "0: turn off geometry and pT matching, 1: do geometry and pT matching"};
  Configurable<bool> doMCPJetSubtraction{"doMCPJetSubtraction", false, "0: turn off subtraction for MCP jets, 1: do subtraction for MCP jets"};

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    AxisSpec centralityAxis = {1200, -10., 110., "Centrality"};
    AxisSpec trackPtAxis = {200, -0.5, 199.5, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec trackEtaAxis = {nBinsEta, -1.0, 1.0, "#eta"};
    AxisSpec phiAxis = {160, -1.0, 7.0, "#varphi"};
    AxisSpec jetPtAxis = {200, 0., 200., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec jetPtAxisRhoAreaSub = {400, -200., 200., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec jetEtaAxis = {nBinsEta, -1.0, 1.0, "#eta"};

    if (doprocessQC || doprocessQCWeighted) {
      registry.add("h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.add("h2_centrality_occupancy", "centrality vs occupancy; centrality; occupancy", {HistType::kTH2F, {centralityAxis, {60, 0, 30000}}});
      registry.add("h_collisions_Zvertex", "position of collision ;#it{Z} (cm)", {HistType::kTH1F, {{300, -15.0, 15.0}}});
      registry.add("h_track_pt", "track pT ; #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH1F, {trackPtAxis}});
      registry.add("h2_track_eta_track_phi", "track eta vs. track phi; #eta; #phi; counts", {HistType::kTH2F, {trackEtaAxis, phiAxis}});
      if (doprocessQCWeighted) {
        registry.add("h_collisions_weighted", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      }
    }

    if (doprocessSpectraData || doprocessSpectraMCD || doprocessSpectraMCDWeighted) {
      registry.add("h_jet_pt_rhoareasubtracted", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_jet_eta_rhoareasubtracted", "jet eta;#eta; counts", {HistType::kTH1F, {jetEtaAxis}});
      registry.add("h_jet_phi_rhoareasubtracted", "jet phi;#phi; counts", {HistType::kTH1F, {phiAxis}});
      registry.add("h2_centrality_jet_pt_rhoareasubtracted", "centrality vs. jet pT;centrality; #it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH2F, {centralityAxis, jetPtAxisRhoAreaSub}});
      registry.add("h2_centrality_jet_eta_rhoareasubtracted", "centrality vs. jet eta;centrality; #eta; counts", {HistType::kTH2F, {centralityAxis, jetEtaAxis}});
      registry.add("h2_centrality_jet_phi_rhoareasubtracted", "centrality vs. jet phi;centrality; #varphi; counts", {HistType::kTH2F, {centralityAxis, phiAxis}});
      registry.add("h2_jet_pt_jet_area_rhoareasubtracted", "jet #it{p}_{T,jet} vs. Area_{jet}; #it{p}_{T,jet} (GeV/#it{c}); Area_{jet}", {HistType::kTH2F, {jetPtAxis, {150, 0., 1.5}}});
      registry.add("h2_jet_pt_jet_ntracks_rhoareasubtracted", "jet #it{p}_{T,jet} vs. N_{jet tracks}; #it{p}_{T,jet} (GeV/#it{c}); N_{jet, tracks}", {HistType::kTH2F, {jetPtAxis, {200, -0.5, 199.5}}});
      registry.add("h2_jet_pt_jet_corr_pt_rhoareasubtracted", "jet #it{p}_{T,jet} vs. #it{p}_{T,corr}; #it{p}_{T,jet} (GeV/#it{c});  #it{p}_{T,corr} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxisRhoAreaSub}});
      registry.add("h2_jet_pt_track_pt_rhoareasubtracted", "jet #it{p}_{T,jet} vs. #it{p}_{T,track}; #it{p}_{T,jet} (GeV/#it{c});  #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, trackPtAxis}});
      registry.add("h3_jet_pt_eta_phi_rhoareasubtracted", "jet_pt_eta_phi_rhoareasubtracted", {HistType::kTH3F, {jetPtAxisRhoAreaSub, jetEtaAxis, phiAxis}});
      if (doprocessSpectraMCDWeighted) {
        registry.add("h_jet_phat", "jet #hat{p};#hat{p} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0, 1000}}});
        registry.add("h_jet_phat_weighted", "jet #hat{p};#hat{p} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0, 1000}}});
      }
    }

    if (doprocessSpectraMCP || doprocessSpectraMCPWeighted) {
      registry.add("h_mc_zvertex", "position of collision ;#it{Z} (cm)", {HistType::kTH1F, {{300, -15.0, 15.0}}});
      registry.add("h_mcp_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxis}});
      registry.add("h_mcp_jet_eta", "jet eta;#eta; counts", {HistType::kTH1F, {jetEtaAxis}});
      registry.add("h_mcp_jet_phi", "jet phi;#phi; counts", {HistType::kTH1F, {phiAxis}});
      registry.add("h2_mcp_jet_pt_jet_area", "jet #it{p}_{T,jet} vs. Area_{jet}; #it{p}_{T,jet} (GeV/#it{c}); Area_{jet}", {HistType::kTH2F, {jetPtAxis, {150, 0., 1.5}}});
      registry.add("h2_mcp_jet_pt_jet_ntracks", "jet #it{p}_{T,jet} vs. N_{jet tracks}; #it{p}_{T,jet} (GeV/#it{c}); N_{jet, tracks}", {HistType::kTH2F, {jetPtAxis, {200, -0.5, 199.5}}});
      registry.add("h2_mcp_jet_pt_track_pt", "jet #it{p}_{T,jet} vs. #it{p}_{T,track}; #it{p}_{T,jet} (GeV/#it{c});  #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, trackPtAxis}});
      registry.add("hsparse_mcp_jet_pt_eta_phi", "mcp jet_pt_eta_phi", HistType::kTHnSparseF, {jetPtAxis, jetEtaAxis, phiAxis});
      if (doMCPJetSubtraction) {
        registry.add("h_mcp_rho", "mc collision rho;#rho (GeV/#it{c}); counts", {HistType::kTH1F, {{500, 0.0, 500.0}}});
        registry.add("h_mcp_jet_corr_pt", "mcp jet corr pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
        registry.add("h2_mcp_jet_corr_pt_jet_area", "mcp jet #it{p}_{T,jet} vs. Area_{jet}; #it{p}_{T,jet} (GeV/#it{c}); Area_{jet}", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {150, 0., 1.5}}});
        registry.add("h2_mcp_jet_corr_pt_jet_ntracks", "mcp jet #it{p}_{T,jet} vs. N_{jet tracks}; #it{p}_{T,jet} (GeV/#it{c}); N_{jet, tracks}", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {200, -0.5, 199.5}}});
        registry.add("hsparse_mcp_jet_corr_pt_eta_phi", "mcp jet_corr_pt_eta_phi", HistType::kTHnSparseF, {jetPtAxisRhoAreaSub, jetEtaAxis, phiAxis});
      }
      if (doprocessSpectraMCPWeighted) {
        registry.add("h2_mcp_jet_ptcut_part", "p_{T} cut;p_{T,jet}^{part} (GeV/#it{c});N;entries", {HistType::kTH2F, {{300, 0, 300}, {20, 0, 5}}});
        registry.add("h_mcp_jet_phat_weighted", "jet #hat{p};#hat{p} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0, 1000}}});
      }
    }

    if (doprocessEvtWiseConstSubJetsData || doprocessEvtWiseConstSubJetsMCD) {
      registry.add("h2_centrality_jet_pt_eventwiseconstituentsubtracted", "centrality vs. jet pT;centrality;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH2F, {centralityAxis, jetPtAxis}});
      registry.add("jet_observables_eventwiseconstituentsubtracted", "jet_observables_eventwiseconstituentsubtracted", HistType::kTHnSparseF, {jetPtAxis, jetEtaAxis, phiAxis});
    }

    if (doprocessJetsMatched || doprocessJetsMatchedWeighted) {
      if (checkGeoMatched) {
        registry.add("h2_jet_pt_reco_jet_pt_gen_matchedgeo", "pT reco vs. pT gen;#it{p}_{T,jet}^{reco} (GeV/#it{c});#it{p}_{T,jet}^{gen} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
        registry.add("h2_jet_eta_reco_jet_eta_gen_matchedgeo", "Eta reco vs. Eta gen;#eta_{jet}^{reco};#eta_{jet}^{gen}", {HistType::kTH2F, {jetEtaAxis, jetEtaAxis}});
        registry.add("h2_jet_phi_reco_jet_phi_gen_matchedgeo", "Phi reco vs. Phi gen;#varphi_{jet}^{reco};#varphi_{jet}^{gen}", {HistType::kTH2F, {phiAxis, phiAxis}});
        registry.add("h2_jet_ntracks_reco_jet_ntracks_gen_matchedgeo", "Ntracks reco vs. Ntracks gen;N_{jet tracks}^{reco};N_{jet tracks}^{gen}", {HistType::kTH2F, {{200, -0.5, 199.5}, {200, -0.5, 199.5}}});
        registry.add("h2_jet_pt_gen_jet_pt_diff_matchedgeo", "jet gen pT vs. delta pT / jet gen pt;#it{p}_{T,jet}^{gen} (GeV/#it{c}); (#it{p}_{T,jet}^{gen} (GeV/#it{c}) - #it{p}_{T,jet}^{reco} (GeV/#it{c})) / #it{p}_{T,jet}^{gen} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_reco_jet_pt_diff_matchedgeo", "jet reco pT vs. delta pT / jet reco pt;#it{p}_{T,jet}^{reco} (GeV/#it{c}); (#it{p}_{T,jet}^{reco} (GeV/#it{c}) - #it{p}_{T,jet}^{gen} (GeV/#it{c})) / #it{p}_{T,jet}^{reco} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_gen_jet_pt_ratio_matchedgeo", "jet gen pT vs. jet reco pT / jet gen pt;#it{p}_{T,jet}^{gen} (GeV/#it{c}); #it{p}_{T,jet}^{reco} (GeV/#it{c}) / #it{p}_{T,jet}^{gen} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 5.0}}});
      }
      if (checkPtMatched) {
        registry.add("h2_jet_pt_reco_jet_pt_gen_matchedpt", "pT reco vs. pT gen;#it{p}_{T,jet}^{reco} (GeV/#it{c});#it{p}_{T,jet}^{gen} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
        registry.add("h2_jet_eta_reco_jet_eta_gen_matchedpt", "Eta reco vs. Eta gen;#eta_{jet}^{reco};#eta_{jet}^{gen}", {HistType::kTH2F, {jetEtaAxis, jetEtaAxis}});
        registry.add("h2_jet_phi_reco_jet_phi_gen_matchedgpt", "Phi reco vs. Phi gen;#varphi_{jet}^{reco};#varphi_{jet}^{gen}", {HistType::kTH2F, {phiAxis, phiAxis}});
        registry.add("h2_jet_ntracks_reco_jet_ntracks_gen_matchedpt", "Ntracks reco vs. Ntracks gen;N_{jet tracks}^{reco};N_{jet tracks}^{gen}", {HistType::kTH2F, {{200, -0.5, 199.5}, {200, -0.5, 199.5}}});
        registry.add("h2_jet_pt_gen_jet_pt_diff_matchedpt", "jet gen pT vs. delta pT / jet gen pt;#it{p}_{T,jet}^{gen} (GeV/#it{c}); (#it{p}_{T,jet}^{gen} (GeV/#it{c}) - #it{p}_{T,jet}^{reco} (GeV/#it{c})) / #it{p}_{T,jet}^{gen} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_reco_jet_pt_diff_matchedpt", "jet reco pT vs. delta pT / jet reco pt;#it{p}_{T,jet}^{reco} (GeV/#it{c}); (#it{p}_{T,jet}^{reco} (GeV/#it{c}) - #it{p}_{T,jet}^{gen} (GeV/#it{c})) / #it{p}_{T,jet}^{reco} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_gen_jet_pt_ratio_matchedpt", "jet gen pT vs. jet reco pT / jet gen pt;#it{p}_{T,jet}^{gen} (GeV/#it{c}); #it{p}_{T,jet}^{reco} (GeV/#it{c}) / #it{p}_{T,jet}^{gen} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 5.0}}});
      }
      if (checkGeoPtMatched) {
        registry.add("h2_jet_pt_reco_jet_pt_gen_matchedgeopt", "pT reco vs. pT gen;#it{p}_{T,jet}^{reco} (GeV/#it{c});#it{p}_{T,jet}^{gen} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
        registry.add("h2_jet_eta_reco_jet_eta_gen_matchedgeopt", "Eta reco vs. Eta gen;#eta_{jet}^{reco};#eta_{jet}^{gen}", {HistType::kTH2F, {jetEtaAxis, jetEtaAxis}});
        registry.add("h2_jet_phi_reco_jet_phi_gen_matchedgeopt", "Phi reco vs. Phi gen;#varphi_{jet}^{reco};#varphi_{jet}^{gen}", {HistType::kTH2F, {phiAxis, phiAxis}});
        registry.add("h2_jet_ntracks_reco_jet_ntracks_gen_matchedgeopt", "Ntracks reco vs. Ntracks gen;N_{jet tracks}^{reco};N_{jet tracks}^{gen}", {HistType::kTH2F, {{200, -0.5, 199.5}, {200, -0.5, 199.5}}});
        registry.add("h2_jet_pt_gen_jet_pt_diff_matchedgeopt", "jet gen pT vs. delta pT / jet gen pt;#it{p}_{T,jet}^{gen} (GeV/#it{c}); (#it{p}_{T,jet}^{gen} (GeV/#it{c}) - #it{p}_{T,jet}^{reco} (GeV/#it{c})) / #it{p}_{T,jet}^{gen} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_reco_jet_pt_diff_matchedgeopt", "jet reco pT vs. delta pT / jet reco pt;#it{p}_{T,jet}^{reco} (GeV/#it{c}); (#it{p}_{T,jet}^{reco} (GeV/#it{c}) - #it{p}_{T,jet}^{gen} (GeV/#it{c})) / #it{p}_{T,jet}^{reco} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_gen_jet_pt_ratio_matchedgeopt", "jet gen pT vs. jet reco pT / jet gen pt;#it{p}_{T,jet}^{gen} (GeV/#it{c}); #it{p}_{T,jet}^{reco} (GeV/#it{c}) / #it{p}_{T,jet}^{gen} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 5.0}}});
      }
    }

    if (doprocessJetsMatchedSubtracted) {
      registry.add("h_mc_collisions_matched", "mc collisions status;event status;entries", {HistType::kTH1F, {{5, 0.0, 5.0}}});
      registry.add("h_mcd_events_matched", "mcd event status;event status;entries", {HistType::kTH1F, {{5, 0.0, 5.0}}});
      registry.add("h_mc_rho_matched", "mc collision rho;#rho (GeV/#it{c}); counts", {HistType::kTH1F, {{500, -100.0, 500.0}}});
      registry.add("h2_reco_jet_corr_pt_gen_jet_corr_pt_matchedgeo", "corr pT reco vs. corr cpT gen;#it{p}_{T,jet}^{reco} (GeV/#it{c});#it{p}_{T,jet}^{gen} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, jetPtAxisRhoAreaSub}});
      registry.add("h2_gen_jet_corr_pt_jet_pt_diff_corr_matchedgeo", "jet gen corr pT vs. corr delta pT / jet gen corr pt;#it{p}_{T,jet}^{gen} (GeV/#it{c}); (#it{p}_{T,jet}^{gen} (GeV/#it{c}) - #it{p}_{T,jet}^{reco} (GeV/#it{c})) / #it{p}_{T,jet}^{gen} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {1000, -5.0, 5.0}}});
      registry.add("h2_reco_jet_corr_pt_jet_pt_diff_corr_matchedgeo", "jet reco corr pT vs. corr delta pT / jet reco corr pt;#it{p}_{T,jet}^{reco} (GeV/#it{c}); (#it{p}_{T,jet}^{reco} (GeV/#it{c}) - #it{p}_{T,jet}^{gen} (GeV/#it{c})) / #it{p}_{T,jet}^{reco} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {1000, -5.0, 5.0}}});
      registry.add("h2_gen_jet_corr_pt_jet_corr_pt_ratio_matchedgeo", "jet gen corr pT vs. jet reco corr pT / jet gen corr pt;#it{p}_{T,jet}^{gen} (GeV/#it{c}); #it{p}_{T,jet}^{reco} (GeV/#it{c}) / #it{p}_{T,jet}^{gen} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, {1000, -5.0, 5.0}}});
    }
  }

  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);
  Preslice<ChargedMCDMatchedJets> mcdjetsPerJCollision = o2::aod::jet::collisionId;

  template <typename TTracks, typename TJets>
  bool isAcceptedJet(TJets const& jet)
  {

    if (jetAreaFractionMin > -98.0) {
      if (jet.area() < jetAreaFractionMin * o2::constants::math::PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    bool checkConstituentPt = true;
    bool checkConstituentMinPt = (leadingConstituentPtMin > -98.0);
    bool checkConstituentMaxPt = (leadingConstituentPtMax < 9998.0);
    if (!checkConstituentMinPt && !checkConstituentMaxPt) {
      checkConstituentPt = false;
    }

    if (checkConstituentPt) {
      bool isMinLeadingConstituent = !checkConstituentMinPt;
      bool isMaxLeadingConstituent = true;

      for (const auto& constituent : jet.template tracks_as<TTracks>()) {
        double pt = constituent.pt();

        if (checkConstituentMinPt && pt >= leadingConstituentPtMin) {
          isMinLeadingConstituent = true;
        }
        if (checkConstituentMaxPt && pt > leadingConstituentPtMax) {
          isMaxLeadingConstituent = false;
        }
      }
      return isMinLeadingConstituent && isMaxLeadingConstituent;
    }
    return true;
  }

  template <typename TJets>
  void fillJetHistograms(TJets const& jet, float centrality, float rho, float weight = 1.0)
  {
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      double jetcorrpt = jet.pt() - (rho * jet.area());
      // fill jet histograms after area-based subtraction
      registry.fill(HIST("h_jet_pt_rhoareasubtracted"), jetcorrpt, weight);
      registry.fill(HIST("h_jet_eta_rhoareasubtracted"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi_rhoareasubtracted"), jet.phi(), weight);
      registry.fill(HIST("h2_centrality_jet_pt_rhoareasubtracted"), centrality, jetcorrpt, weight);
      registry.fill(HIST("h2_centrality_jet_eta_rhoareasubtracted"), centrality, jet.eta(), weight);
      registry.fill(HIST("h2_centrality_jet_phi_rhoareasubtracted"), centrality, jet.phi(), weight);
      registry.fill(HIST("h2_jet_pt_jet_corr_pt_rhoareasubtracted"), jet.pt(), jetcorrpt, weight);
      registry.fill(HIST("h3_jet_pt_eta_phi_rhoareasubtracted"), jetcorrpt, jet.eta(), jet.phi(), weight);
      if (jetcorrpt > 0) {
        registry.fill(HIST("h2_jet_pt_jet_area_rhoareasubtracted"), jetcorrpt, jet.area(), weight);
        registry.fill(HIST("h2_jet_pt_jet_ntracks_rhoareasubtracted"), jetcorrpt, jet.tracksIds().size(), weight);
      }
    }

    for (const auto& constituent : jet.template tracks_as<aod::JetTracks>()) {
      registry.fill(HIST("h2_jet_pt_track_pt_rhoareasubtracted"), jet.pt(), constituent.pt(), weight);
    }
  }

  template <typename TJets>
  void fillMCPHistograms(TJets const& jet, float rho = 0.0, float weight = 1.0)
  {
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      // fill mcp jet histograms
      registry.fill(HIST("h_mcp_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_mcp_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_mcp_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("hsparse_mcp_jet_pt_eta_phi"), jet.pt(), jet.eta(), jet.phi(), weight);
      registry.fill(HIST("h2_mcp_jet_pt_jet_area"), jet.pt(), jet.area(), weight);
      registry.fill(HIST("h2_mcp_jet_pt_jet_ntracks"), jet.pt(), jet.tracksIds().size(), weight);
      if (doMCPJetSubtraction) {
        double jetcorrpt = jet.pt() - (rho * jet.area());
        registry.fill(HIST("h_mcp_jet_corr_pt"), jetcorrpt, weight);
        registry.fill(HIST("hsparse_mcp_jet_corr_pt_eta_phi"), jet.pt(), jet.eta(), jet.phi(), weight);
        if (jetcorrpt > 0) {
          registry.fill(HIST("h2_mcp_jet_corr_pt_jet_area"), jetcorrpt, jet.area(), weight);
          registry.fill(HIST("h2_mcp_jet_corr_pt_jet_ntracks"), jetcorrpt, jet.tracksIds().size(), weight);
        }
      }
    }

    for (const auto& constituent : jet.template tracks_as<aod::JetParticles>()) {
      registry.fill(HIST("h2_mcp_jet_pt_track_pt"), jet.pt(), constituent.pt(), weight);
    }
  }

  template <typename TJets>
  void fillEventWiseConstituentSubtractedHistograms(TJets const& jet, float centrality, float weight = 1.0)
  {
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h2_centrality_jet_pt_eventwiseconstituentsubtracted"), centrality, jet.pt(), weight);
      registry.fill(HIST("jet_observables_eventwiseconstituentsubtracted"), jet.pt(), jet.eta(), jet.phi(), weight);
    }
  }

  template <typename TTracks>
  void fillTrackHistograms(TTracks const& track, float weight = 1.0)
  {
    registry.fill(HIST("h_track_pt"), track.pt(), weight);
    registry.fill(HIST("h2_track_eta_track_phi"), track.eta(), track.phi(), weight);
  }

  template <typename TBase, typename TTag>
  void fillMatchedHistograms(TBase const& jetBase, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jetBase.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    // fill geometry matched histograms
    if (checkGeoMatched) {
      if (jetBase.has_matchedJetGeo()) {
        for (const auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<TTag>>()) {
          if (jetTag.pt() > pTHatMaxMCD * pTHat) {
            continue;
          }
          if (jetBase.r() == round(selectedJetsRadius * 100.0f)) {
            double dpt = jetTag.pt() - jetBase.pt();
            registry.fill(HIST("h2_jet_pt_reco_jet_pt_gen_matchedgeo"), jetBase.pt(), jetTag.pt(), weight);
            registry.fill(HIST("h2_jet_eta_reco_jet_eta_gen_matchedgeo"), jetBase.eta(), jetTag.eta(), weight);
            registry.fill(HIST("h2_jet_phi_reco_jet_phi_gen_matchedgeo"), jetBase.phi(), jetTag.phi(), weight);
            registry.fill(HIST("h2_jet_ntracks_reco_jet_ntracks_gen_matchedgeo"), jetBase.tracksIds().size(), jetTag.tracksIds().size(), weight);
            registry.fill(HIST("h2_jet_pt_gen_jet_pt_diff_matchedgeo"), jetTag.pt(), dpt / jetTag.pt(), weight);
            registry.fill(HIST("h2_jet_pt_reco_jet_pt_diff_matchedgeo"), jetBase.pt(), dpt / jetBase.pt(), weight);
            registry.fill(HIST("h2_jet_pt_gen_jet_pt_ratio_matchedgeo"), jetTag.pt(), jetBase.pt() / jetTag.pt(), weight);
          }
        }
      }
    }
    // fill pt matched histograms
    if (checkPtMatched) {
      if (jetBase.has_matchedJetPt()) {
        for (const auto& jetTag : jetBase.template matchedJetPt_as<std::decay_t<TTag>>()) {
          if (jetTag.pt() > pTHatMaxMCD * pTHat) {
            continue;
          }
          if (jetBase.r() == round(selectedJetsRadius * 100.0f)) {
            double dpt = jetTag.pt() - jetBase.pt();
            registry.fill(HIST("h2_jet_pt_reco_jet_pt_gen_matchedpt"), jetBase.pt(), jetTag.pt(), weight);
            registry.fill(HIST("h2_jet_eta_reco_jet_eta_gen_matchedpt"), jetBase.eta(), jetTag.eta(), weight);
            registry.fill(HIST("h2_jet_phi_reco_jet_phi_gen_matchedpt"), jetBase.phi(), jetTag.phi(), weight);
            registry.fill(HIST("h2_jet_ntracks_reco_jet_ntracks_gen_matchedpt"), jetBase.tracksIds().size(), jetTag.tracksIds().size(), weight);
            registry.fill(HIST("h2_jet_pt_gen_jet_pt_diff_matchedpt"), jetTag.pt(), dpt / jetTag.pt(), weight);
            registry.fill(HIST("h2_jet_pt_reco_jet_pt_diff_matchedpt"), jetBase.pt(), dpt / jetBase.pt(), weight);
            registry.fill(HIST("h2_jet_pt_gen_jet_pt_ratio_matchedpt"), jetTag.pt(), jetBase.pt() / jetTag.pt(), weight);
          }
        }
      }
    }
    // fill geometry and pt histograms
    if (checkGeoPtMatched) {
      if (jetBase.has_matchedJetGeo() && jetBase.has_matchedJetPt()) {
        for (const auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<TTag>>()) {
          if (jetTag.pt() > pTHatMaxMCD * pTHat) {
            continue;
          }
          if (jetBase.template matchedJetGeo_first_as<std::decay_t<TTag>>().globalIndex() == jetBase.template matchedJetPt_first_as<std::decay_t<TTag>>().globalIndex()) { // not a good way to do this
            double dpt = jetTag.pt() - jetBase.pt();
            registry.fill(HIST("h2_jet_pt_reco_jet_pt_gen_matchedgeopt"), jetBase.pt(), jetTag.pt(), weight);
            registry.fill(HIST("h2_jet_eta_reco_jet_eta_gen_matchedgeopt"), jetBase.eta(), jetTag.eta(), weight);
            registry.fill(HIST("h2_jet_phi_reco_jet_phi_gen_matchedgeopt"), jetBase.phi(), jetTag.phi(), weight);
            registry.fill(HIST("h2_jet_ntracks_reco_jet_ntracks_gen_matchedgeopt"), jetBase.tracksIds().size(), jetTag.tracksIds().size(), weight);
            registry.fill(HIST("h2_jet_pt_gen_jet_pt_diff_matchedgeopt"), jetTag.pt(), dpt / jetTag.pt(), weight);
            registry.fill(HIST("h2_jet_pt_reco_jet_pt_diff_matchedgeopt"), jetBase.pt(), dpt / jetBase.pt(), weight);
            registry.fill(HIST("h2_jet_pt_gen_jet_pt_ratio_matchedgeopt"), jetTag.pt(), jetBase.pt() / jetTag.pt(), weight);
          }
        }
      }
    }
  }

  template <typename TBase, typename TTag>
  void fillGeoMatchedCorrHistograms(TBase const& jetBase, float rho, float mcrho = 0.0, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jetBase.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    if (jetBase.has_matchedJetGeo()) {
      for (const auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<TTag>>()) {
        if (jetTag.pt() > pTHatMaxMCD * pTHat) {
          continue;
        }
        if (jetBase.r() == round(selectedJetsRadius * 100.0f)) {
          double corrTagjetpt = jetTag.pt() - (mcrho * jetTag.area());
          double corrBasejetpt = jetBase.pt() - (rho * jetBase.area());
          double dcorrpt = corrTagjetpt - corrBasejetpt;
          registry.fill(HIST("h2_reco_jet_corr_pt_gen_jet_corr_pt_matchedgeo"), corrBasejetpt, corrTagjetpt, weight);
          registry.fill(HIST("h2_gen_jet_corr_pt_jet_pt_diff_corr_matchedgeo"), corrTagjetpt, dcorrpt / corrTagjetpt, weight);
          registry.fill(HIST("h2_reco_jet_corr_pt_jet_pt_diff_corr_matchedgeo"), corrBasejetpt, dcorrpt / corrBasejetpt, weight);
          registry.fill(HIST("h2_gen_jet_corr_pt_jet_corr_pt_ratio_matchedgeo"), corrTagjetpt, corrBasejetpt / corrTagjetpt, weight);
        }
      }
    }
  }

  void processQC(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                 soa::Filtered<soa::Join<aod::JetTracks, aod::JTrackExtras>> const& tracks)
  {
    registry.fill(HIST("h_collisions"), 0.5);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collisions"), 2.5);
    registry.fill(HIST("h2_centrality_occupancy"), collision.centrality(), collision.trackOccupancyInTimeRange());
    registry.fill(HIST("h_collisions_Zvertex"), collision.posZ());

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      fillTrackHistograms(track, collision.centrality());
    }
  }
  PROCESS_SWITCH(JetSpectraCharged, processQC, "collisions and track QC for Data and MCD", true);

  void processQCWeighted(soa::Join<aod::JetCollisions, aod::JMcCollisionLbs>::iterator const& collision,
                         aod::JetMcCollisions const&,
                         soa::Filtered<soa::Join<aod::JetTracks, aod::JTrackExtras>> const& tracks)
  {
    float eventWeight = collision.mcCollision().weight();
    registry.fill(HIST("h_collisions"), 0.5);
    registry.fill(HIST("h_collisions_weighted"), 0.5, eventWeight);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    registry.fill(HIST("h_collisions_weighted"), 1.5, eventWeight);
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collisions"), 2.5);
    registry.fill(HIST("h_collisions_weighted"), 2.5, eventWeight);
    registry.fill(HIST("h_collisions_Zvertex"), collision.posZ(), eventWeight);

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      fillTrackHistograms(track, eventWeight);
    }
  }
  PROCESS_SWITCH(JetSpectraCharged, processQCWeighted, "weighted collsions and tracks QC for MC", false);

  void processSpectraData(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision,
                          soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                          aod::JetTracks const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      fillJetHistograms(jet, collision.centrality(), collision.rho());
    }
  }
  PROCESS_SWITCH(JetSpectraCharged, processSpectraData, "jet spectra for Data", false);

  void processSpectraMCD(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision,
                         soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                         aod::JetTracks const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      fillJetHistograms(jet, collision.centrality(), collision.rho());
    }
  }
  PROCESS_SWITCH(JetSpectraCharged, processSpectraMCD, "jet spectra for Data", false);

  void processSpectraMCDWeighted(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision,
                                 soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetEventWeights> const& jets,
                                 aod::JetTracks const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      float jetweight = jet.eventWeight();
      float pTHat = 10. / (std::pow(jetweight, 1.0 / pTHatExponent));
      if (jet.pt() > pTHatMaxMCD * pTHat) {
        return;
      }
      registry.fill(HIST("h_jet_phat"), pTHat);
      registry.fill(HIST("h_jet_phat_weighted"), pTHat, jetweight);
      fillJetHistograms(jet, collision.centrality(), collision.rho(), jetweight);
    }
  }
  PROCESS_SWITCH(JetSpectraCharged, processSpectraMCDWeighted, "jet finder QA mcd with weighted events", false);

  void processSpectraMCP(McParticleCollision::iterator const& mccollision,
                         soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
                         aod::JetParticles const&)
  {
    if (std::abs(mccollision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("h_mc_zvertex"), mccollision.posZ());
    if (doMCPJetSubtraction) {
      registry.fill(HIST("h_mcp_rho"), mccollision.rho());
    }
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetParticles>(jet)) {
        continue;
      }
      if (doMCPJetSubtraction) {
        fillMCPHistograms(jet, mccollision.rho());
      } else {
        fillMCPHistograms(jet);
      }
    }
  }
  PROCESS_SWITCH(JetSpectraCharged, processSpectraMCP, "jet spectra for MC particle level", false);

  void processSpectraMCPWeighted(aod::JetMcCollisions::iterator const& mccollision,
                                 soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetEventWeights> const& jets,
                                 aod::JetParticles const&)
  {
    if (!(std::abs(mccollision.posZ()) < vertexZCut)) {
      return;
    }
    registry.fill(HIST("h_mc_zvertex"), mccollision.posZ());
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetParticles>(jet)) {
        continue;
      }
      float jetweight = jet.eventWeight();
      double pTHat = 10. / (std::pow(jetweight, 1.0 / pTHatExponent));
      for (int N = 1; N < 21; N++) {
        if (jet.pt() < N * 0.25 * pTHat && jet.r() == round(selectedJetsRadius * 100.0f)) {
          registry.fill(HIST("h2_mcp_jet_ptcut_part"), jet.pt(), N * 0.25, jetweight);
        }
      }
      registry.fill(HIST("h_mcp_jet_phat_weighted"), pTHat, jetweight);
      fillMCPHistograms(jet, 0.0, jetweight);
    }
  }
  PROCESS_SWITCH(JetSpectraCharged, processSpectraMCPWeighted, "jet spectra for MC particle level weighted", false);

  void processEvtWiseConstSubJetsData(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                                      soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents> const& jets,
                                      aod::JetTracksSub const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracksSub>(jet)) {
        continue;
      }
      fillEventWiseConstituentSubtractedHistograms(jet, collision.centrality());
    }
  }
  PROCESS_SWITCH(JetSpectraCharged, processEvtWiseConstSubJetsData, "jet spectrum for eventwise constituent-subtracted jets data", false);

  void processEvtWiseConstSubJetsMCD(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                                     soa::Join<aod::ChargedMCDetectorLevelEventWiseSubtractedJets, aod::ChargedMCDetectorLevelEventWiseSubtractedJetConstituents> const& jets,
                                     aod::JetTracksSub const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracksSub>(jet)) {
        continue;
      }
      fillEventWiseConstituentSubtractedHistograms(jet, collision.centrality());
    }
  }
  PROCESS_SWITCH(JetSpectraCharged, processEvtWiseConstSubJetsMCD, "jet spectrum for eventwise constituent-subtracted mcd jets", false);

  void processJetsMatched(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                          ChargedMCDMatchedJets const& mcdjets,
                          ChargedMCPMatchedJets const&,
                          aod::JetTracks const&, aod::JetParticles const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }

    for (const auto& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      fillMatchedHistograms<ChargedMCDMatchedJets::iterator, ChargedMCPMatchedJets>(mcdjet);
    }
  }
  PROCESS_SWITCH(JetSpectraCharged, processJetsMatched, "matched mcp and mcd jets", false);

  void processJetsMatchedWeighted(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                                  ChargedMCDMatchedJetsWeighted const& mcdjets,
                                  ChargedMCPMatchedJetsWeighted const&,
                                  aod::JetTracks const&, aod::JetParticles const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    for (const auto& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      fillMatchedHistograms<ChargedMCDMatchedJetsWeighted::iterator, ChargedMCPMatchedJetsWeighted>(mcdjet, mcdjet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetSpectraCharged, processJetsMatchedWeighted, "matched mcp and mcd jets with weighted events", false);

  void processJetsMatchedSubtracted(McParticleCollision::iterator const& mccollision,
                                    soa::SmallGroups<soa::Join<aod::JetCollisionsMCD, aod::BkgChargedRhos>> const& collisions,
                                    ChargedMCDMatchedJets const& mcdjets,
                                    ChargedMCPMatchedJets const&,
                                    aod::JetTracks const&, aod::JetParticles const&)
  {
    registry.fill(HIST("h_mc_collisions_matched"), 0.5);
    if (mccollision.size() < 1) {
      return;
    }
    registry.fill(HIST("h_mc_collisions_matched"), 1.5);
    if (!(std::abs(mccollision.posZ()) < vertexZCut)) {
      return;
    }
    registry.fill(HIST("h_mc_collisions_matched"), 2.5);
    double mcrho = mccollision.rho();
    registry.fill(HIST("h_mc_rho_matched"), mcrho);
    for (const auto& collision : collisions) {
      registry.fill(HIST("h_mcd_events_matched"), 0.5);
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || !(std::abs(collision.posZ()) < vertexZCut)) {
        continue;
      }
      registry.fill(HIST("h_mcd_events_matched"), 1.5);
      if (collision.centrality() < centralityMin || collision.centrality() > centralityMax) {
        continue;
      }
      registry.fill(HIST("h_mcd_events_matched"), 2.5);
      if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
        return;
      }
      registry.fill(HIST("h_mcd_events_matched"), 3.5);

      auto collmcdjets = mcdjets.sliceBy(mcdjetsPerJCollision, collision.globalIndex());
      for (const auto& mcdjet : collmcdjets) {
        if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
          continue;
        }
        if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
          continue;
        }
        if (doMCPJetSubtraction) {
          fillGeoMatchedCorrHistograms<ChargedMCDMatchedJets::iterator, ChargedMCPMatchedJets>(mcdjet, collision.rho(), mcrho);
        } else {
          fillGeoMatchedCorrHistograms<ChargedMCDMatchedJets::iterator, ChargedMCPMatchedJets>(mcdjet, collision.rho());
        }
      }
    }
  }
  PROCESS_SWITCH(JetSpectraCharged, processJetsMatchedSubtracted, "matched mcp and mcd jets after subtraction", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetSpectraCharged>(cfgc)};
}
