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

// Charged-particle jet spectra task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
/// \author Aimeric Landou <aimeric.landou@cern.ch>
/// \author Wenhui Feng <wenhui.feng@cern.ch>

#include <cmath>
#include <TRandom3.h>
#include <THn.h>
#include <THnSparse.h>
#include <string>

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

struct JetSpectraChargedTask {

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
  Configurable<bool> checkPtMatched{"checkPtMatched", true, "0: turn off pT matching, 1: do pT matching"};
  Configurable<bool> checkGeoPtMatched{"checkGeoPtMatched", true, "0: turn off geometry and pT matching, 1: do geometry and pT matching"};

  int eventSelection = -1;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    AxisSpec centralityAxis = {1200, -10., 110., "Centrality"};
    AxisSpec trackPtAxis = {200, -0.5, 199.5, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec trackEtaAxis = {nBinsEta, trackEtaMin, trackEtaMax, "#eta"};
    AxisSpec PhiAxis = {160, -1.0, 7.0, "#varphi"};
    AxisSpec jetPtAxis = {200, 0., 200., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec jetPtAxisRhoAreaSub = {400, -200., 200., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec jetEtaAxis = {nBinsEta, jetEtaMin, jetEtaMax, "#eta"};

    if (doprocessQC || doprocessQCWeighted) {
      registry.add("h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.add("h2_centrality_occupancy", "centrality vs occupancy; centrality; occupancy", {HistType::kTH2F, {centralityAxis, {60, 0, 30000}}});
      registry.add("h_collisions_vertexZ", "position of collision ;#it{Z} (cm)", {HistType::kTH1F, {{300, -15.0, 15.0}}});
      registry.add("h_track_pt", "track pT ; #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH1F, {trackPtAxis}});
      registry.add("h2_track_eta_track_phi", "track eta vs. track phi; #eta; #phi; counts", {HistType::kTH2F, {trackEtaAxis, PhiAxis}});
      registry.add("h_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h_jet_eta", "jet eta;#eta; counts", {HistType::kTH1F, {jetEtaAxis}});
      registry.add("h_jet_phi", "jet phi;#phi; counts", {HistType::kTH1F, {PhiAxis}});
      if (doprocessQCWeighted) {
        registry.add("h_collisions_weighted", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      }
    }

    if (doprocessSpectraData || doprocessSpectraMCD || doprocessSpectraMCDWeighted) {
      registry.add("h_jet_pt_rhoareasubtracted", "jet pT;#it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h2_centrality_jet_pt_rhoareasubtracted", "centrality vs. jet pT;centrality; #it{p}_{T,jet} (GeV/#it{c}); counts", {HistType::kTH2F, {centralityAxis, jetPtAxisRhoAreaSub}});
      registry.add("h2_centrality_jet_eta_rhoareasubtracted", "centrality vs. jet eta;centrality; #eta; counts", {HistType::kTH2F, {centralityAxis, jetEtaAxis}});
      registry.add("h2_centrality_jet_phi_rhoareasubtracted", "centrality vs. jet phi;centrality; #varphi; counts", {HistType::kTH2F, {centralityAxis, PhiAxis}});
      registry.add("h2_jet_pt_jet_area_rhoareasubtracted", "jet #it{p}_{T,jet} vs. Area_{jet}; #it{p}_{T,jet} (GeV/#it{c}); Area_{jet}", {HistType::kTH2F, {jetPtAxis, {150, 0., 1.5}}});
      registry.add("h2_jet_pt_jet_ntracks_rhoareasubtracted", "jet #it{p}_{T,jet} vs. N_{jet tracks}; #it{p}_{T,jet} (GeV/#it{c}); N_{jet, tracks}", {HistType::kTH2F, {jetPtAxis, {200, -0.5, 199.5}}});
      registry.add("h2_jet_pt_jet_corr_pt_rhoareasubtracted", "jet #it{p}_{T,jet} vs. #it{p}_{T,corr}; #it{p}_{T,jet} (GeV/#it{c});  #it{p}_{T,corr} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxisRhoAreaSub}});
      registry.add("h2_jet_pt_track_pt_rhoareasubtracted", "jet #it{p}_{T,jet} vs. #it{p}_{T,track}; #it{p}_{T,jet} (GeV/#it{c});  #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxisRhoAreaSub, trackPtAxis}});
      registry.add("h3_jet_pt_eta_phi_rhoareasubtracted", "jet_pt_eta_phi_rhoareasubtracted", {HistType::kTH3F, {jetPtAxisRhoAreaSub, jetEtaAxis, PhiAxis}});
      if (doprocessSpectraMCDWeighted) {
        registry.add("h_jet_phat", "jet #hat{p};#hat{p} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0, 1000}}});
        registry.add("h_jet_phat_weighted", "jet #hat{p};#hat{p} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0, 1000}}});
      }
    }

    if (doprocessEvtWiseConstSubJetsData || doprocessEvtWiseConstSubJetsMCD) {
      registry.add("h2_centrality_jet_pt_eventwiseconstituentsubtracted", "centrality vs. jet pT;centrality;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH2F, {centralityAxis, jetPtAxis}});
      registry.add("jet_observables_eventwiseconstituentsubtracted", "jet_observables_eventwiseconstituentsubtracted", HistType::kTHnSparseF, {jetPtAxis, jetEtaAxis, PhiAxis});
    }

    if (doprocessJetsMatched || doprocessJetsMatchedWeighted) {
      if (checkGeoMatched) {
        registry.add("h2_jet_pt_tag_jet_pt_base_matchedgeo", "pT tag vs. pT base;#it{p}_{T,jet}^{tag} (GeV/#it{c});#it{p}_{T,jet}^{base} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
        registry.add("h2_jet_eta_tag_jet_eta_base_matchedgeo", "Eta tag vs. Eta base;#eta_{jet}^{tag};#eta_{jet}^{base}", {HistType::kTH2F, {jetEtaAxis, jetEtaAxis}});
        registry.add("h2_jet_phi_tag_jet_phi_base_matchedgeo", "Phi tag vs. Phi base;#varphi_{jet}^{tag};#varphi_{jet}^{base}", {HistType::kTH2F, {PhiAxis, PhiAxis}});
        registry.add("h2_jet_ntracks_tag_jet_ntracks_base_matchedgeo", "Ntracks tag vs. Ntracks base;N_{jet tracks}^{tag};N_{jet tracks}^{base}", {HistType::kTH2F, {{200, -0.5, 199.5}, {200, -0.5, 199.5}}});
        registry.add("h2_jet_pt_tag_jet_pt_diff_matchedgeo", "jet tag pT vs. delta pT / jet tag pt;#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_base_jet_pt_diff_matchedgeo", "jet base pT vs. delta pT / jet base pt;#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
      }
      if (checkPtMatched) {
        registry.add("h2_jet_pt_tag_jet_pt_base_matchedpt", "pT tag vs. pT base;#it{p}_{T,jet}^{tag} (GeV/#it{c});#it{p}_{T,jet}^{base} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
        registry.add("h2_jet_eta_tag_jet_eta_base_matchedpt", "Eta tag vs. Eta base;#eta_{jet}^{tag};#eta_{jet}^{base}", {HistType::kTH2F, {jetEtaAxis, jetEtaAxis}});
        registry.add("h2_jet_phi_tag_jet_phi_base_matchedpt", "Phi tag vs. Phi base;#varphi_{jet}^{tag};#varphi_{jet}^{base}", {HistType::kTH2F, {PhiAxis, PhiAxis}});
        registry.add("h2_jet_ntracks_tag_jet_ntracks_base_matchedpt", "Ntracks tag vs. Ntracks base;N_{jet tracks}^{tag};N_{jet tracks}^{base}", {HistType::kTH2F, {{200, -0.5, 199.5}, {200, -0.5, 199.5}}});
        registry.add("h2_jet_pt_tag_jet_pt_diff_matchedpt", "jet tag pT vs. delta pT / jet tag pt;#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_base_jet_pt_diff_matchedpt", "jet base pT vs. delta pT / jet base pt;#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
      }
      if (checkGeoPtMatched) {
        registry.add("h2_jet_pt_tag_jet_pt_base_matchedgeopt", "pT tag vs. pT base;#it{p}_{T,jet}^{tag} (GeV/#it{c});#it{p}_{T,jet}^{base} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
        registry.add("h2_jet_eta_tag_jet_eta_base_matchedgeopt", "Eta tag vs. Eta base;#eta_{jet}^{tag};#eta_{jet}^{base}", {HistType::kTH2F, {jetEtaAxis, jetEtaAxis}});
        registry.add("h2_jet_phi_tag_jet_phi_base_matchedgeopt", "Phi tag vs. Phi base;#varphi_{jet}^{tag};#varphi_{jet}^{base}", {HistType::kTH2F, {PhiAxis, PhiAxis}});
        registry.add("h2_jet_ntracks_tag_jet_ntracks_base_matchedgeopt", "Ntracks tag vs. Ntracks base;N_{jet tracks}^{tag};N_{jet tracks}^{base}", {HistType::kTH2F, {{200, -0.5, 199.5}, {200, -0.5, 199.5}}});
        registry.add("h2_jet_pt_tag_jet_pt_diff_matchedgeopt", "jet tag pT vs. delta pT / jet tag pt;#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_base_jet_pt_diff_matchedgeopt", "jet base pT vs. delta pT / jet base pt;#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
      }
    }

    if (doprocessJetsMatchedSubtracted) {
      registry.add("h2_jet_pt_tag_jet_pt_base_corr_matchedgeo", "pT tag vs. corr pT base;#it{p}_{T,jet}^{tag} (GeV/#it{c});#it{p}_{T,jet}^{base} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
      registry.add("h2_jet_pt_tag_jet_pt_diff_corr_matchedgeo", "jet tag pT vs. corr delta pT / jet tag pt;#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
      registry.add("h2_jet_pt_base_jet_pt_diff_corr_matchedgeo", "jet base pT vs. corr delta pT / jet base pt;#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
    }
  }

  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);

  template <typename TTracks, typename TJets>
  bool isAcceptedJet(TJets const& jet)
  {

    if (jetAreaFractionMin > -98.0) {
      if (jet.area() < jetAreaFractionMin * M_PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
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
      Double_t jetcorrpt = jet.pt() - (rho * jet.area());
      // fill jet histograms after area-based subtraction
      registry.fill(HIST("h_jet_pt_rhoareasubtracted"), jetcorrpt, weight);
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

    for (auto& constituent : jet.template tracks_as<aod::JetTracks>()) {
      registry.fill(HIST("h2_jet_pt_track_pt_rhoareasubtracted"), jet.pt(), constituent.pt(), weight);
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
        for (auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<TTag>>()) {
          if (jetTag.pt() > pTHatMaxMCD * pTHat) {
            continue;
          }
          if (jetBase.r() == round(selectedJetsRadius * 100.0f)) {
            Double_t dpt = jetTag.pt() - jetBase.pt();
            registry.fill(HIST("h2_jet_pt_tag_jet_pt_base_matchedgeo"), jetTag.pt(), jetBase.pt(), weight);
            registry.fill(HIST("h2_jet_eta_tag_jet_eta_base_matchedgeo"), jetTag.eta(), jetBase.eta(), weight);
            registry.fill(HIST("h2_jet_phi_tag_jet_phi_base_matchedgeo"), jetTag.phi(), jetBase.phi(), weight);
            registry.fill(HIST("h2_jet_ntracks_tag_jet_ntracks_base_matchedgeo"), jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
            registry.fill(HIST("h2_jet_pt_tag_jet_pt_diff_matchedgeo"), jetTag.pt(), dpt / jetTag.pt(), weight);
            registry.fill(HIST("h2_jet_pt_base_jet_pt_diff_matchedgeo"), jetBase.pt(), dpt / jetBase.pt(), weight);
          }
        }
      }
    }
    // fill pt matched histograms
    if (checkPtMatched) {
      if (jetBase.has_matchedJetPt()) {
        for (auto& jetTag : jetBase.template matchedJetPt_as<std::decay_t<TTag>>()) {
          if (jetTag.pt() > pTHatMaxMCD * pTHat) {
            continue;
          }
          if (jetBase.r() == round(selectedJetsRadius * 100.0f)) {
            Double_t dpt = jetTag.pt() - jetBase.pt();
            registry.fill(HIST("h2_jet_pt_tag_jet_pt_base_matchedpt"), jetTag.pt(), jetBase.pt(), weight);
            registry.fill(HIST("h2_jet_eta_tag_jet_eta_base_matchedpt"), jetTag.eta(), jetBase.eta(), weight);
            registry.fill(HIST("h2_jet_phi_tag_jet_phi_base_matchedpt"), jetTag.phi(), jetBase.phi(), weight);
            registry.fill(HIST("h2_jet_ntracks_tag_jet_ntracks_base_matchedpt"), jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
            registry.fill(HIST("h2_jet_pt_tag_jet_pt_diff_matchedpt"), jetTag.pt(), dpt / jetTag.pt(), weight);
            registry.fill(HIST("h2_jet_pt_base_jet_pt_diff_matchedpt"), jetBase.pt(), dpt / jetBase.pt(), weight);
          }
        }
      }
    }
    // fill geometry and pt histograms:
    if (checkGeoPtMatched) {
      if (jetBase.has_matchedJetGeo() && jetBase.has_matchedJetPt()) {
        for (auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<TTag>>()) {
          if (jetTag.pt() > pTHatMaxMCD * pTHat) {
            continue;
          }
          if (jetBase.template matchedJetGeo_first_as<std::decay_t<TTag>>().globalIndex() == jetBase.template matchedJetPt_first_as<std::decay_t<TTag>>().globalIndex()) { // not a good way to do this
            Double_t dpt = jetTag.pt() - jetBase.pt();
            registry.fill(HIST("h2_jet_pt_tag_jet_pt_base_matchedgeopt"), jetTag.pt(), jetBase.pt(), weight);
            registry.fill(HIST("h2_jet_eta_tag_jet_eta_base_matchedgeopt"), jetTag.eta(), jetBase.eta(), weight);
            registry.fill(HIST("h2_jet_phi_tag_jet_phi_base_matchedgeopt"), jetTag.phi(), jetBase.phi(), weight);
            registry.fill(HIST("h2_jet_ntracks_tag_jet_ntracks_base_matchedgeopt"), jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
            registry.fill(HIST("h2_jet_pt_tag_jet_pt_diff_matchedgeopt"), jetTag.pt(), dpt / jetTag.pt(), weight);
            registry.fill(HIST("h2_jet_pt_base_jet_pt_diff_matchedgeopt"), jetBase.pt(), dpt / jetBase.pt(), weight);
          }
        }
      }
    }
  }

  template <typename TBase, typename TTag>
  void fillGeoMatchedCorrHistograms(TBase const& jetBase, float rho, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jetBase.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    if (jetBase.has_matchedJetGeo()) {
      for (auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<TTag>>()) {
        if (jetTag.pt() > pTHatMaxMCD * pTHat) {
          continue;
        }
        if (jetBase.r() == round(selectedJetsRadius * 100.0f)) {
          Double_t corrBasejetpt = jetBase.pt() - (rho * jetBase.area());
          Double_t dcorrpt = jetTag.pt() - corrBasejetpt;
          registry.fill(HIST("h2_jet_pt_tag_jet_pt_base_corr_matchedgeo"), jetTag.pt(), corrBasejetpt, weight);
          registry.fill(HIST("h2_jet_pt_tag_jet_pt_diff_corr_matchedgeo"), jetTag.pt(), dcorrpt / jetTag.pt(), weight);
          registry.fill(HIST("h2_jet_pt_base_jet_pt_diff_corr_matchedgeo"), corrBasejetpt, dcorrpt / corrBasejetpt, weight);
        }
      }
    }
  }

  void processQC(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                 soa::Filtered<soa::Join<aod::JetTracks, aod::JTrackExtras>> const& tracks,
                 soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets)
  {
    registry.fill(HIST("h_collisions"), 0.5);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collisions"), 2.5);
    registry.fill(HIST("h2_centrality_occupancy"), collision.centrality(), collision.trackOccupancyInTimeRange());
    registry.fill(HIST("h_collisions_vertexZ"), collision.posZ());

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      fillTrackHistograms(track, collision.centrality());
    }

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      registry.fill(HIST("h_jet_pt"), jet.pt());
      registry.fill(HIST("h_jet_eta"), jet.eta());
      registry.fill(HIST("h_jet_phi"), jet.phi());
    }
  }
  PROCESS_SWITCH(JetSpectraChargedTask, processQC, "collisions and track QC for Data and MCD", true);

  void processQCWeighted(soa::Join<aod::JetCollisions, aod::JMcCollisionLbs>::iterator const& collision,
                         aod::JetMcCollisions const&,
                         soa::Filtered<soa::Join<aod::JetTracks, aod::JTrackExtras>> const& tracks,
                         soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetEventWeights> const& jets)
  {
    float eventWeight = collision.mcCollision().weight();
    registry.fill(HIST("h_collisions"), 0.5);
    registry.fill(HIST("h_collisions_weighted"), 0.5, eventWeight);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    registry.fill(HIST("h_collisions_weighted"), 1.5, eventWeight);
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collisions"), 2.5);
    registry.fill(HIST("h_collisions_weighted"), 2.5, eventWeight);
    registry.fill(HIST("h_collisions_vertexZ"), collision.posZ(), eventWeight);

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      fillTrackHistograms(track, eventWeight);
    }

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      registry.fill(HIST("h_jet_pt"), jet.pt(), eventWeight);
      registry.fill(HIST("h_jet_eta"), jet.eta(), eventWeight);
      registry.fill(HIST("h_jet_phi"), jet.phi(), eventWeight);
    }
  }
  PROCESS_SWITCH(JetSpectraChargedTask, processQCWeighted, "weighted collsions and tracks QC for MC", false);

  void processSpectraData(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision,
                          soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                          aod::JetTracks const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
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
  PROCESS_SWITCH(JetSpectraChargedTask, processSpectraData, "jet spectra for Data", false);

  void processSpectraMCD(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision,
                         soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                         aod::JetTracks const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
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
  PROCESS_SWITCH(JetSpectraChargedTask, processSpectraMCD, "jet spectra for Data", false);

  void processSpectraMCDWeighted(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision,
                                 soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetEventWeights> const& jets,
                                 aod::JetTracks const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
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
  PROCESS_SWITCH(JetSpectraChargedTask, processSpectraMCDWeighted, "jet finder QA mcd with weighted events", false);

  void processEvtWiseConstSubJetsData(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                                      soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents> const& jets,
                                      aod::JetTracksSub const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
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
  PROCESS_SWITCH(JetSpectraChargedTask, processEvtWiseConstSubJetsData, "jet spectrum for eventwise constituent-subtracted jets data", false);

  void processEvtWiseConstSubJetsMCD(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                                     soa::Join<aod::ChargedMCDetectorLevelEventWiseSubtractedJets, aod::ChargedMCDetectorLevelEventWiseSubtractedJetConstituents> const& jets,
                                     aod::JetTracksSub const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
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
  PROCESS_SWITCH(JetSpectraChargedTask, processEvtWiseConstSubJetsMCD, "jet spectrum for eventwise constituent-subtracted mcd jets", false);

  void processJetsMatched(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                          ChargedMCDMatchedJets const& mcdjets,
                          ChargedMCPMatchedJets const&,
                          aod::JetTracks const&, aod::JetParticles const&)
  {
    registry.fill(HIST("h_collisions"), 0.5);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collisions"), 2.5);

    for (const auto& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      registry.fill(HIST("h_jet_pt_mcd"), mcdjet.pt());
      registry.fill(HIST("h_jet_eta_mcd"), mcdjet.eta());
      registry.fill(HIST("h_jet_phi_mcd"), mcdjet.phi());
      fillMatchedHistograms<ChargedMCDMatchedJets::iterator, ChargedMCPMatchedJets>(mcdjet);
    }
  }
  PROCESS_SWITCH(JetSpectraChargedTask, processJetsMatched, "matched mcp and mcd jets", false);

  void processJetsMatchedWeighted(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                                  ChargedMCDMatchedJetsWeighted const& mcdjets,
                                  ChargedMCPMatchedJetsWeighted const&,
                                  aod::JetTracks const&, aod::JetParticles const&)
  {
    registry.fill(HIST("h_collisions"), 0.5);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collisions"), 2.5);
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
  PROCESS_SWITCH(JetSpectraChargedTask, processJetsMatchedWeighted, "matched mcp and mcd jets with weighted events", false);

  void processJetsMatchedSubtracted(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision,
                                    ChargedMCDMatchedJets const& mcdjets,
                                    ChargedMCPMatchedJets const&,
                                    aod::JetTracks const&, aod::JetParticles const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
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
      // now only do subtraction for MCD jets, need to add subtraction for MCP jets
      fillGeoMatchedCorrHistograms<ChargedMCDMatchedJets::iterator, ChargedMCPMatchedJets>(mcdjet, collision.rho());
    }
  }
  PROCESS_SWITCH(JetSpectraChargedTask, processJetsMatchedSubtracted, "matched mcp and mcd jets after subtraction", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetSpectraChargedTask>(cfgc, TaskName{"jet-charged-spectra"})};
}
