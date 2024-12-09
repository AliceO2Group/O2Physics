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


struct JetMatchedResponseTask {

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
  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
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

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));

    AxisSpec jetPtAxis = {200, 0., 200., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec jetPtAxisRhoAreaSub = {400, -200., 200., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec jetEtaAxis = {nBinsEta, jetEtaMin, jetEtaMax, "#eta"};
    AxisSpec jetPhiAxis = {160, -1.0, 7.0, "#varphi"};


    if (doprocessJetsMatched || doprocessJetsMatchedWeighted) {
      registry.add("h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.add("h_jet_pt_mcd", "pT base ;#it{p}_{T,jet}^{base} (GeV/#it{c})", {HistType::kTH1F, {jetPtAxis}});
      registry.add("h_jet_eta_mcd", "Eta base ;#eta^{base}", {HistType::kTH1F, {jetEtaAxis}});
      registry.add("h_jet_phi_mcd", "Phi base ;#varphi^{base}", {HistType::kTH1F, {jetPhiAxis}});
      if(checkGeoMatched){
        registry.add("h2_jet_pt_tag_jet_pt_base_matchedgeo", "pT tag vs. pT base;#it{p}_{T,jet}^{tag} (GeV/#it{c});#it{p}_{T,jet}^{base} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
        registry.add("h2_jet_eta_tag_jet_eta_base_matchedgeo", "Eta tag vs. Eta base;#eta_{jet}^{tag};#eta_{jet}^{base}", {HistType::kTH2F, {jetEtaAxis, jetEtaAxis}});
        registry.add("h2_jet_phi_tag_jet_phi_base_matchedgeo", "Phi tag vs. Phi base;#varphi_{jet}^{tag};#varphi_{jet}^{base}", {HistType::kTH2F, {jetPhiAxis, jetPhiAxis}});
        registry.add("h2_jet_ntracks_tag_jet_ntracks_base_matchedgeo", "Ntracks tag vs. Ntracks base;N_{jet tracks}^{tag};N_{jet tracks}^{base}", {HistType::kTH2F, {{200, -0.5, 199.5}, {200, -0.5, 199.5}}});
        registry.add("h2_jet_pt_tag_jet_pt_diff_matchedgeo", "jet tag pT vs. delta pT / jet tag pt;#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_base_jet_pt_diff_matchedgeo", "jet base pT vs. delta pT / jet base pt;#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h3_ptcut_jet_pt_tag_jet_pt_base_matchedgeo", "N;#it{p}_{T,jet}^{tag} (GeV/#it{c});#it{p}_{T,jet}^{base} (GeV/#it{c})", {HistType::kTH3F, {{20, 0., 5.}, jetPtAxis, jetPtAxis}});

      }
      if(checkPtMatched){
        registry.add("h2_jet_pt_tag_jet_pt_base_matchedpt", "pT tag vs. pT base;#it{p}_{T,jet}^{tag} (GeV/#it{c});#it{p}_{T,jet}^{base} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
        registry.add("h2_jet_eta_tag_jet_eta_base_matchedpt", "Eta tag vs. Eta base;#eta_{jet}^{tag};#eta_{jet}^{base}", {HistType::kTH2F, {jetEtaAxis, jetEtaAxis}});
        registry.add("h2_jet_phi_tag_jet_phi_base_matchedpt", "Phi tag vs. Phi base;#varphi_{jet}^{tag};#varphi_{jet}^{base}", {HistType::kTH2F, {jetPhiAxis, jetPhiAxis}});
        registry.add("h2_jet_ntracks_tag_jet_ntracks_base_matchedpt", "Ntracks tag vs. Ntracks base;N_{jet tracks}^{tag};N_{jet tracks}^{base}", {HistType::kTH2F, {{200, -0.5, 199.5}, {200, -0.5, 199.5}}});
        registry.add("h2_jet_pt_tag_jet_pt_diff_matchedpt", "jet tag pT vs. delta pT / jet tag pt;#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_base_jet_pt_diff_matchedpt", "jet base pT vs. delta pT / jet base pt;#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
      }
      if(checkGeoPtMatched){
        registry.add("h2_jet_pt_tag_jet_pt_base_matchedgeopt", "pT tag vs. pT base;#it{p}_{T,jet}^{tag} (GeV/#it{c});#it{p}_{T,jet}^{base} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
        registry.add("h2_jet_eta_tag_jet_eta_base_matchedgeopt", "Eta tag vs. Eta base;#eta_{jet}^{tag};#eta_{jet}^{base}", {HistType::kTH2F, {jetEtaAxis, jetEtaAxis}});
        registry.add("h2_jet_phi_tag_jet_phi_base_matchedgeopt", "Phi tag vs. Phi base;#varphi_{jet}^{tag};#varphi_{jet}^{base}", {HistType::kTH2F, {jetPhiAxis, jetPhiAxis}});
        registry.add("h2_jet_ntracks_tag_jet_ntracks_base_matchedgeopt", "Ntracks tag vs. Ntracks base;N_{jet tracks}^{tag};N_{jet tracks}^{base}", {HistType::kTH2F, {{200, -0.5, 199.5}, {200, -0.5, 199.5}}});
        registry.add("h2_jet_pt_tag_jet_pt_diff_matchedgeopt", "jet tag pT vs. delta pT / jet tag pt;#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
        registry.add("h2_jet_pt_base_jet_pt_diff_matchedgeopt", "jet base pT vs. delta pT / jet base pt;#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
      }
    }

    if(doprocessJetsMatchedSubtracted){
      registry.add("h2_jet_pt_tag_jet_pt_base_corr_matchedgeo", "pT tag vs. corr pT base;#it{p}_{T,jet}^{tag} (GeV/#it{c});#it{p}_{T,jet}^{base} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
      registry.add("h2_jet_pt_tag_jet_pt_diff_corr_matchedgeo", "jet tag pT vs. corr delta pT / jet tag pt;#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});
      registry.add("h2_jet_pt_base_jet_pt_diff_corr_matchedgeo", "jet base pT vs. corr delta pT / jet base pt;#it{p}_{T,jet}^{tag} (GeV/#it{c}); (#it{p}_{T,jet}^{tag} (GeV/#it{c}) - #it{p}_{T,jet}^{base} (GeV/#it{c})) / #it{p}_{T,jet}^{tag} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, {1000, -5.0, 2.0}}});

    }

  }

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

  template <typename T, typename U>
  void fillGeoMatchedHistograms(T const& jetBase, float weight = 1.0)
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
        if (jetBase.r() == round(selectedJetsRadius * 100.0f)) {
          Double_t dpt = jetTag.pt() - jetBase.pt();
          registry.fill(HIST("h2_jet_pt_tag_jet_pt_base_matchedgeo"), jetTag.pt(), jetBase.pt(), weight);
          registry.fill(HIST("h2_jet_eta_tag_jet_eta_base_matchedgeo"), jetTag.eta(), jetBase.eta(), weight);
          registry.fill(HIST("h2_jet_phi_tag_jet_phi_base_matchedgeo"), jetTag.phi(), jetBase.phi(), weight);
          registry.fill(HIST("h2_jet_ntracks_tag_jet_ntracks_base_matchedgeo"), jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
          registry.fill(HIST("h2_jet_pt_tag_jet_pt_diff_matchedgeo"),  jetTag.pt(), dpt / jetTag.pt(), weight);
          registry.fill(HIST("h2_jet_pt_base_jet_pt_diff_matchedgeo"),  jetBase.pt(), dpt / jetBase.pt(), weight);

          for (int N = 1; N < 21; N++) {
            if (jetBase.pt() < N * 0.25 * pTHat && jetTag.pt() < N * 0.25 * pTHat) {
              registry.fill(HIST("h3_ptcut_jet_pt_tag_jet_pt_base_matchedgeo"), N * 0.25, jetTag.pt(), jetBase.pt(), weight);
            }
          }
        }
      }
    }
  }

  template <typename T, typename U>
  void fillPtMatchedHistograms(T const& jetBase, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jetBase.pt() > pTHatMaxMCD * pTHat) {
      return;
    }

    if (jetBase.has_matchedJetPt()) {
      for (auto& jetTag : jetBase.template matchedJetPt_as<std::decay_t<U>>()) {
        if (jetTag.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }
        if (jetBase.r() == round(selectedJetsRadius * 100.0f)) {
          Double_t dpt = jetTag.pt() - jetBase.pt();
          registry.fill(HIST("h2_jet_pt_tag_jet_pt_base_matchedpt"), jetTag.pt(), jetBase.pt(), weight);
          registry.fill(HIST("h2_jet_eta_tag_jet_eta_base_matchedpt"), jetTag.eta(), jetBase.eta(), weight);
          registry.fill(HIST("h2_jet_phi_tag_jet_phi_base_matchedpt"), jetTag.phi(), jetBase.phi(), weight);
          registry.fill(HIST("h2_jet_ntracks_tag_jet_ntracks_base_matchedpt"), jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
          registry.fill(HIST("h2_jet_pt_tag_jet_pt_diff_matchedpt"),  jetTag.pt(), dpt / jetTag.pt(), weight);
          registry.fill(HIST("h2_jet_pt_base_jet_pt_diff_matchedpt"),  jetBase.pt(), dpt / jetBase.pt(), weight);
        }
      }
    }
  }

  template <typename T, typename U>
  void fillGeoPtMatchedHistograms(T const& jetBase, float weight = 1.0)
  {
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jetBase.pt() > pTHatMaxMCD * pTHat) {
      return;
    }

    if (jetBase.has_matchedJetGeo() && jetBase.has_matchedJetPt()) {

      for (auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<U>>()) {
        if (jetTag.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }

        if (jetBase.template matchedJetGeo_first_as<std::decay_t<U>>().globalIndex() == jetBase.template matchedJetPt_first_as<std::decay_t<U>>().globalIndex()) { // not a good way to do this
          Double_t dpt = jetTag.pt() - jetBase.pt();
          registry.fill(HIST("h2_jet_pt_tag_jet_pt_base_matchedgeopt"), jetTag.pt(), jetBase.pt(), weight);
          registry.fill(HIST("h2_jet_eta_tag_jet_eta_base_matchedgeopt"), jetTag.eta(), jetBase.eta(), weight);
          registry.fill(HIST("h2_jet_phi_tag_jet_phi_base_matchedgeopt"), jetTag.phi(), jetBase.phi(), weight);
          registry.fill(HIST("h2_jet_ntracks_tag_jet_ntracks_base_matchedgeopt"), jetTag.tracksIds().size(), jetBase.tracksIds().size(), weight);
          registry.fill(HIST("h2_jet_pt_tag_jet_pt_diff_matchedgeopt"),  jetTag.pt(), dpt / jetTag.pt(), weight);
          registry.fill(HIST("h2_jet_pt_base_jet_pt_diff_matchedgeopt"),  jetBase.pt(), dpt / jetBase.pt(), weight);

        }
      }
    }
  }

  template <typename T, typename U>
  void fillGeoMatchedCorrHistograms(T const& jetBase, float rho, float weight = 1.0)
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
        if (jetBase.r() == round(selectedJetsRadius * 100.0f)) {
          Double_t corrBasejetpt = jetBase.pt() - (rho* jetBase.area());
          Double_t dcorrpt = jetTag.pt() - corrBasejetpt;
          registry.fill(HIST("h2_jet_pt_tag_jet_pt_base_corr_matchedgeo"), jetTag.pt(), corrBasejetpt, weight);
          registry.fill(HIST("h2_jet_pt_tag_jet_pt_diff_corr_matchedgeo"),  jetTag.pt(), dcorrpt / jetTag.pt(), weight);
          registry.fill(HIST("h2_jet_pt_base_jet_pt_diff_corr_matchedgeo"),  corrBasejetpt, dcorrpt / corrBasejetpt, weight);
        }
      }
    }
  }

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
      if(checkGeoMatched){
        fillGeoMatchedHistograms<ChargedMCDMatchedJets::iterator, ChargedMCPMatchedJets>(mcdjet);
      }
      if(checkPtMatched){
        fillPtMatchedHistograms<ChargedMCDMatchedJets::iterator, ChargedMCPMatchedJets>(mcdjet);
      }
      if(checkGeoPtMatched){//select mcd jet has both geometry and pt matched mcp jets, but not a good way to do this
        fillGeoPtMatchedHistograms<ChargedMCDMatchedJets::iterator, ChargedMCPMatchedJets>(mcdjet);
      }

    }
  }
  PROCESS_SWITCH(JetMatchedResponseTask, processJetsMatched, "matched mcp and mcd jets", true);

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
      registry.fill(HIST("h_jet_pt_mcd"), mcdjet.pt(), mcdjet.eventWeight());
      registry.fill(HIST("h_jet_eta_mcd"), mcdjet.eta(), mcdjet.eventWeight());
      registry.fill(HIST("h_jet_phi_mcd"), mcdjet.phi(), mcdjet.eventWeight());
      if(checkGeoMatched){
        fillGeoMatchedHistograms<ChargedMCDMatchedJetsWeighted::iterator, ChargedMCPMatchedJetsWeighted>(mcdjet, mcdjet.eventWeight());
      }
      if(checkPtMatched){
        fillPtMatchedHistograms<ChargedMCDMatchedJetsWeighted::iterator, ChargedMCPMatchedJetsWeighted>(mcdjet, mcdjet.eventWeight());
      }
      if(checkGeoPtMatched){//select mcd jet has both geometry and pt matched mcp jets, but not a good way to do this
        fillGeoPtMatchedHistograms<ChargedMCDMatchedJetsWeighted::iterator, ChargedMCPMatchedJetsWeighted>(mcdjet, mcdjet.eventWeight());
      }

    }
  }
  PROCESS_SWITCH(JetMatchedResponseTask, processJetsMatchedWeighted, "matched mcp and mcd jets with weighted events", false);

  void processJetsMatchedSubtracted(soa::Filtered<soa::Join<aod::JetCollisions,aod::BkgChargedRhos>>::iterator const& collision,
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
  PROCESS_SWITCH(JetMatchedResponseTask, processJetsMatchedSubtracted, "matched mcp and mcd jets after subtraction", false);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
  return WorkflowSpec{adaptAnalysisTask<JetMatchedResponseTask>(cfgc, TaskName{"jet-matched-response"})};
  }
