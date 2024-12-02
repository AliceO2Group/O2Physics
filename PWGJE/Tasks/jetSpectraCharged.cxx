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

struct JetSpectraChargedTask {
  Service<o2::framework::O2DatabasePDG> pdg;

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
  Configurable<double> jetPtMax{"jetPtMax", 200., "set jet pT bin max"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
  Configurable<int> nBinsEta{"nBinsEta", 200, "number of bins for eta axes"};
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};
  Configurable<float> leadingConstituentPtMax{"leadingConstituentPtMax", 9999.0, "maximum pT selection on jet constituent"};
  Configurable<bool> checkMcCollisionIsMatched{"checkMcCollisionIsMatched", false, "0: count whole MCcollisions, 1: select MCcollisions which only have their correspond collisions"};
  Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  Configurable<int> trackOccupancyInTimeRangeMin{"trackOccupancyInTimeRangeMin", -999999, "minimum occupancy of tracks in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};

  int eventSelection = -1;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    
    AxisSpec centralityAxis = {1200, -10., 110., "Centrality"};
    AxisSpec trackPtAxis = {200, -0.5, 199.5, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec trackEtaAxis = {nBinsEta, trackEtaMin, trackEtaMax, "#eta"};

    AxisSpec jetPtAxis = {200, 0., 200., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec jetPtAxisRhoAreaSub = {400, -200., 200., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec jetEtaAxis = {nBinsEta, jetEtaMin, jetEtaMax, "#eta"};
    AxisSpec PhiAxis = {160, -1.0, 7.0, "#varphi"};

    if (doprocessCollisions) {
      registry.add("h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.add("h2_centrality_collisions", "centrality vs collisions; centrality; collisions", {HistType::kTH2F, {centralityAxis, {4, 0.0, 4.0}}});
      registry.add("h2_centrality_occupancy", "centrality vs occupancy; centrality; occupancy", {HistType::kTH2F, {centralityAxis, {60, 0, 3000}}});
      registry.add("h_collisions_vertexZ", "position of collision ;#it{Z} (cm)", {HistType::kTH1F, {{300, -15.0, 15.0}}});
       if (doprocessCollisionsWeighted) {
        registry.add("h_collisions_weighted", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      }
    }

    if (doprocessTracks || doprocessTracksWeighted) {
      registry.add("h_track_pt", "track pT ; #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH1F, {trackPtAxis}});
      registry.add("h2_track_pt_track_dcaxy", "track pT vs track DCA_{xy}; #it{p}_{T,track} (GeV/#it{c}); track DCA_{xy}", {HistType::kTH2F, {trackPtAxis, {200, -0.15, 0.15}}});
      registry.add("h2_occupancy_Ntracks", "occupancy vs Number of selected track; occupancy; N_{tracks}", {HistType::kTH2F, {{60, 0, 30000}, {10000, 0., 10000.}}});
      registry.add("Centrality_track_observables", "Centrality_jet_observables", HistType::kTHnSparseF, {centralityAxis, trackPtAxis, trackEtaAxis, PhiAxis});
    }

    if (doprocessJetsData || doprocessJetsMCD || doprocessJetsMCDWeighted) {
      registry.add("h2_centrality_jet_pt", "centrality vs. #it{p}_{T,jet}; centrality; #it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH2F, {centralityAxis,jetPtAxis}});
      registry.add("h2_jet_pt_jet_ntracks", "jet #it{p}_{T,jet} vs. N_{jet tracks}; #it{p}_{T,jet} (GeV/#it{c}); N_{jet, tracks}", {HistType::kTH2F, {jetPtAxis, {200, -0.5, 199.5}}});
      registry.add("h3_jet_pt_jet_eta_jet_phi", "jet #it{p}_{T,jet} vs. #eta_{jet} vs. #varphi_{jet}; #it{p}_{T,jet} (GeV/#it{c}); #eta_{jet}; #varphi_{jet}", {HistType::kTH3F, {jetPtAxis, jetEtaAxis, PhiAxis}});
      registry.add("h2_jet_pt_track_pt", "jet #it{p}_{T,jet} vs. N_{jet tracks}; #it{p}_{T,jet} (GeV/#it{c}); N_{jet, tracks}", {HistType::kTH2F, {jetPtAxis, trackPtAxis}});
      registry.add("h_jet_phat", "jet #hat{p};#hat{p} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0, 1000}}});
      registry.add("h_jet_phat_weighted", "jet #hat{p};#hat{p} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0, 1000}}});
    }

    if (doprocessJetsRhoAreaSubData || doprocessJetsRhoAreaSubMCD) {
      registry.add("h_jet_pt_rhoareasubtracted", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxisRhoAreaSub}});
      registry.add("h2_centrality_jet_ntracks_rhoareasubtracted", "centrality vs N_{jet tracks}; centrality; N_{jet tracks}", {HistType::kTH2F, {centralityAxis, {200, -0.5, 199.5}}});
      registry.add("h2_jet_pt_jet_area_rhoareasubtracted", "jet #it{p}_{T,jet} vs. Area_{jet}; #it{p}_{T,jet} (GeV/#it{c}); Area_{jet}", {HistType::kTH2F, {jetPtAxis, {150, 0., 1.5}}});
      registry.add("h2_jet_pt_jet_ntracks_rhoareasubtracted", "jet #it{p}_{T,jet} vs. N_{jet tracks}; #it{p}_{T,jet} (GeV/#it{c}); N_{jet, tracks}", {HistType::kTH2F, {jetPtAxis, {200, -0.5, 199.5}}});
      registry.add("h2_jet_pt_jet_corr_pt_rhoareasubtracted", "jet #it{p}_{T,jet} vs. #it{p}_{T,corr}; #it{p}_{T,jet} (GeV/#it{c});  #it{p}_{T,corr} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, jetPtAxisRhoAreaSub}});
      registry.add("h2_occupancy_njets_rhoareasubtracted", "occupancy vs. number of selected jets ; occupancy; N_{jets}", {HistType::kTH2F, {{60, 0, 30000}, {5000, 0, 5000}}});
      registry.add("Centrality_jet_observables_rhoareasubtracted", "Centrality_jet_observables_rhoareasubtracted", HistType::kTHnSparseF, {centralityAxis, jetPtAxisRhoAreaSub, jetEtaAxis, PhiAxis});
    }

    if (doprocessEvtWiseConstSubJetsData || doprocessEvtWiseConstSubJetsMCD) {
      registry.add("h_jet_pt_eventwiseconstituentsubtracted", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxis}});
      registry.add("Centrality_jet_observables_eventwiseconstituentsubtracted", "Centrality_jet_observables_eventwiseconstituentsubtracted", HistType::kTHnSparseF, {centralityAxis, jetPtAxis, jetEtaAxis, PhiAxis});
    }

    if (doprocessRho) {
      registry.add("h2_centrality_ntracks", "; centrality; N_{tracks};", {HistType::kTH2F, {{1100, 0., 110.0}, {10000, 0.0, 10000.0}}});
      registry.add("h2_ntracks_rho", "; N_{tracks}; #it{rho} (GeV/area);", {HistType::kTH2F, {{10000, 0.0, 10000.0}, {400, 0.0, 400.0}}});
      registry.add("h2_centrality_rho", "; centrality; #it{rho} (GeV/area);", {HistType::kTH2F, {{1100, 0., 110.}, {400, 0., 400.0}}});
    }

    if (doprocessJetsMCP || doprocessJetsMCPWeighted) {
      registry.add("h_jet_pt_part", "jet pT;#it{p}_{T,jet}^{part}(GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxis}});
      registry.add("Jet_observables_part", "Jet_observables_part", HistType::kTHnSparseF, {jetPtAxis, jetEtaAxis, PhiAxis});
      registry.add("h2_jet_pt_part_tracks_pt_part", "#it{p}_{T,jet}^{part} vs. constituent #it{p}_{T, tracks}; #it{p}_{T,jet}^{part} (GeV/#it{c}); #it{p}_{T, tracks} (GeV/#it{c})", {HistType::kTH2F, {jetPtAxis, trackPtAxis}});
      registry.add("h_jet_phat_part", "jet #hat{p};#hat{p} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0, 1000}}});
      registry.add("h_jet_phat_part_weighted", "jet #hat{p};#hat{p} (GeV/#it{c});entries", {HistType::kTH1F, {{1000, 0, 1000}}});
    }
    if (doprocessMcparticles || doprocessMcparticlesWeighted){
      registry.add("h_mccollisions_part","event status; entries", {HistType::kTH1F, {{5, 0, 5}}});
      registry.add("h2_centrality_mccollisions_part","event status vs. centrality; centrality; entries", {HistType::kTH2F, {centralityAxis, {5, 0, 5}}});
      registry.add("h2_centrality_particle_pt","centrality vs. particle #it{p}_{T}; centrality; #it{p}_{T, track} (GeV/#it{c})", {HistType::kTH2F, {centralityAxis, {200, 0., 200.}}});
      registry.add("h2_centrality_particle_eta","centrality vs. particle #eta; centrality; #eta_{track}", {HistType::kTH2F, {centralityAxis, trackEtaAxis}});
      registry.add("h2_centrality_particle_phi","centrality vs. particle #varphi; centrality; #varphi_{track}", {HistType::kTH2F, {centralityAxis, PhiAxis}});
      registry.add("h2_centrality_particle_energy","centrality vs. particle energy; centrality; energy GeV", {HistType::kTH2F, {centralityAxis, {100, 0., 100.}}});
      if(doprocessMcparticlesWeighted){
      registry.add("h_mccollisions_part_weighted","event status; event status; entries", {HistType::kTH1F, {{5, 0, 5}}});
      }
    }
  }

  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);
  Filter particlecuts = (aod::jmcparticle::pt >= trackPtMin && aod::jmcparticle::pt < trackPtMax && aod::jmcparticle::eta > trackEtaMin && aod::jmcparticle::eta < trackEtaMax);
  
  PresliceUnsorted<soa::Filtered<aod::JetCollisionsMCD>> CollisionsPerMCPCollision = aod::jmccollisionlb::mcCollisionId;
  // Preslice<aod::JetTracksMCD> tracksPerJCollision = o2::aod::jtracks::collisionId
  
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

  bool isChargedParticle(int code){
    auto p = pdg->GetParticle(code);
    auto charge = 0.;
    if(p != nullptr){
      charge = p->Charge();
    }
    return std::abs(charge) >= 3.;
  }

  template <typename TTracks, typename TJets>
  bool trackIsInJet(TTracks const& track, TJets const& jet)
  {
    for (auto const& constituentId : jet.tracksIds()) {
      if (constituentId == track.globalIndex()) {
        return true;
      }
    }
    return false;
  }

  template <typename TJets>
  void fillHistograms(TJets const& jet, float centrality, float weight = 1.0)
  {

    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCD * pTHat) {
      return;
    }
    registry.fill(HIST("h_jet_phat"), pTHat);
    registry.fill(HIST("h_jet_phat_weighted"), pTHat, weight);

    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h2_centrality_jet_pt"), centrality, jet.pt(), weight);
      registry.fill(HIST("h2_jet_pt_jet_ntracks"), jet.pt(), jet.tracksIds().size(), weight);
      registry.fill(HIST("h3_jet_pt_jet_eta_jet_phi"), jet.pt(), jet.eta(), jet.phi(), weight);
    }


    for (auto& constituent : jet.template tracks_as<aod::JetTracks>()) {
      registry.fill(HIST("h2_jet_pt_track_pt"), jet.pt(), constituent.pt(), weight);
    }
  }

  template <typename TJets>
  void fillRhoAreaSubtractedHistograms(TJets const& jet, float centrality, float occupancy, float rho, float weight = 1.0)
  {
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_jet_pt_rhoareasubtracted"), jet.pt() - (rho * jet.area()), weight);
      registry.fill(HIST("h2_jet_pt_jet_corr_pt_rhoareasubtracted"), jet.pt(), jet.pt() - (rho * jet.area()), weight);
      registry.fill(HIST("Centrality_jet_observables_rhoareasubtracted"), centrality, jet.pt() - (rho * jet.area()), jet.eta(), jet.phi(), weight);
      if (jet.pt() - (rho * jet.area()) > 0) {
        registry.fill(HIST("h2_centrality_jet_ntracks_rhoareasubtracted"), centrality, jet.tracksIds().size(), weight);
        registry.fill(HIST("h2_jet_pt_jet_area_rhoareasubtracted"), jet.pt() - (rho * jet.area()), jet.area(), weight);
        registry.fill(HIST("h2_jet_pt_jet_ntracks_rhoareasubtracted"), jet.pt() - (rho * jet.area()), jet.tracksIds().size(), weight);
      }
    }
  }

  template <typename TJets>
  void fillEventWiseConstituentSubtractedHistograms(TJets const& jet, float centrality, float weight = 1.0)
  {
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_jet_pt_eventwiseconstituentsubtracted"), jet.pt(), weight);
      registry.fill(HIST("Centrality_jet_observables_eventwiseconstituentsubtracted"), centrality, jet.pt(), jet.eta(), jet.phi(), weight);
    }
  }

  template <typename TJets>
  void fillMCPHistograms(TJets const& jet, float weight = 1.0)
  {

    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    if (jet.pt() > pTHatMaxMCP * pTHat) {
      return;
    }
    registry.fill(HIST("h_jet_phat_part"), pTHat);
    registry.fill(HIST("h_jet_phat_part_weighted"), pTHat, weight);

    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_jet_pt_part"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta_part"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi_part"), jet.phi(), weight);
      registry.fill(HIST("h_jet_ntracks_part"), jet.tracksIds().size(), weight);
      
      for (auto& constituent : jet.template tracks_as<aod::JetParticles>()) {
        registry.fill(HIST("h2_jet_pt_part_tracks_pt_part"), jet.pt(), constituent.pt(), weight);
        registry.fill(HIST("h2_jet_pt_part_tracks_eta_part"), jet.pt(), constituent.eta(), weight);
        registry.fill(HIST("h2_jet_pt_part_tracks_phi_part"), jet.pt(), constituent.phi(), weight);
      }
    }
  }

  template <typename TTracks>
  void fillTrackHistograms(TTracks const& tracks, float centrality, float weight = 1.0)
  {
      registry.fill(HIST("h_track_pt"), track.pt(), weight);
      registry.fill(HIST("h2_track_pt_track_dcaxy"), track.pt(), track.dcaXY(), weight);
      registry.fill(HIST("Centrality_track_observables"), centrality, track.pt(), track.eta(), track.phi(), weight);
    
  }

  template<typename TCollision, typename Tmcparticle>
  void fillParticleHistograms(TCollision const& collision, Tmcparticle const& mcparticles, float weight = 1.0){
    for (auto const& mcparticle : mcparticles){
      if (!isChargedParticle(mcparticle.pdgCode())){
        continue;
      }
      registry.fill(HIST("h2_centrality_particle_pt"), collision.centrality(), mcparticle.pt(), weight);
      registry.fill(HIST("h2_centrality_particle_eta"), collision.centrality(), mcparticle.eta(), weight);
      registry.fill(HIST("h2_centrality_particle_phi"), collision.centrality(), mcparticle.phi(), weight);
      registry.fill(HIST("h2_centrality_particle_energy"), collision.centrality(), mcparticle.energy(), weight);
    }
  }
  
  void processCollisions(soa::Filtered<aod::JetCollisions>::iterator const& collision)
  {
    registry.fill(HIST("h_collisions"), 0.5);
    registry.fill(HIST("h2_centrality_collisions"), collision.centrality(), 0.5);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    registry.fill(HIST("h2_centrality_collisions"), collision.centrality(), 1.5);
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    registry.fill(HIST("h_collisions"), 2.5);
    registry.fill(HIST("h2_centrality_collisions"), collision.centrality(), 2.5);
    registry.fill(HIST("h2_centrality_occupancy"), collision.centrality(), collision.trackOccupancyInTimeRange());
    registry.fill(HIST("h_collisions_vertexZ"), collision.posZ());
      
  }
  PROCESS_SWITCH(JetSpectraChargedTask, processCollisions, "collision QA for events", true);
  
  void processCollisionsWeighted(soa::Join<aod::JetCollisions, aod::JMcCollisionLbs>::iterator const& collision,
                                 aod::JetMcCollisions const&)
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
      
  }
  PROCESS_SWITCH(JetSpectraChargedTask, processCollisionsWeighted, "collision for weighted events", true);

  void processTracks(soa::Filtered<aod::JetCollisions>::iterator const& collision,
                     soa::Filtered<soa::Join<aod::JetTracks, aod::JTrackExtras>> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    
    int nTracks = 0;
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      nTracks++;
      fillTrackHistograms(collision.centrality(), tracks);
      
    }
    registry.fill(HIST("h2_occupancy_Ntracks"), collision.trackOccupancyInTimeRange(), nTracks);
  }
  PROCESS_SWITCH(JetSpectraChargedTask, processTracks, "QA for charged tracks", false);

  void processTracksWeighted(soa::Join<aod::JetCollisions, aod::JMcCollisionLbs>::iterator const& collision,
                             aod::JetMcCollisions const&,
                             soa::Filtered<soa::Join<aod::JetTracks, aod::JTrackExtras>> const& tracks)
  {
    float eventWeight = collision.mcCollision().weight();
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
   
    int nTracks = 0;
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      nTracks++;
      fillTrackHistograms(collision.centrality(), tracks);
    }
    registry.fill(HIST("h2_occupancy_Ntracks"), collision.trackOccupancyInTimeRange(), nTracks);
  }
  PROCESS_SWITCH(JetSpectraChargedTask, processTracksWeighted, "QA for charged tracks weighted", false);

  void processJetsData(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets, aod::JetTracks const&)
  {
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
      fillHistograms(jet, collision.centrality());
    }
  }
  PROCESS_SWITCH(JetSpectraChargedTask, processJetsData, "jet finder QA data", false);

  void processJetsRhoAreaSubData(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision,
                                 soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                                 aod::JetTracks const&)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    int nJets = 0;
    for (auto jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      nJets++;
      fillRhoAreaSubtractedHistograms(jet, collision.centrality(), collision.trackOccupancyInTimeRange(), collision.rho());
    }
    registry.fill(HIST("h2_occupancy_njets_rhoareasubtracted"), collision.trackOccupancyInTimeRange(), nJets);
    
  }
  PROCESS_SWITCH(JetSpectraChargedTask, processJetsRhoAreaSubData, "jet finder QA for rho-area subtracted jets", false);

  void processJetsRhoAreaSubMCD(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision,
                                soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
                                aod::JetTracks const&)
  {
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    int nJets = 0;
    for (auto jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      nJets++;
      fillRhoAreaSubtractedHistograms(jet, collision.centrality(), collision.trackOccupancyInTimeRange(), collision.rho());
    }
    registry.fill(HIST("h2_occupancy_njets_rhoareasubtracted"), collision.trackOccupancyInTimeRange(), nJets);
  }
  PROCESS_SWITCH(JetSpectraChargedTask, processJetsRhoAreaSubMCD, "jet finder QA for rho-area subtracted mcd jets", false);

void processEvtWiseConstSubJetsData(soa::Filtered<aod::JetCollisions>::iterator const& collision, 
                                    soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents> const& jets, 
                                    aod::JetTracksSub const&)
  {
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

  void processJetsMCD(soa::Filtered<aod::JetCollisions>::iterator const& collision, 
                      soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets, 
                      aod::JetTracks const&)
  {
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
      fillHistograms(jet, collision.centrality());
    }
  }
  PROCESS_SWITCH(JetSpectraChargedTask, processJetsMCD, "jet finder QA mcd", false);

  void processJetsMCDWeighted(soa::Filtered<aod::JetCollisions>::iterator const& collision, 
                              soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetEventWeights> const& jets, 
                              aod::JetTracks const&)
  {
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
      fillHistograms(jet, collision.centrality(), jet.eventWeight());
    }
  }
  PROCESS_SWITCH(JetSpectraChargedTask, processJetsMCDWeighted, "jet finder QA mcd with weighted events", false);

  void processJetsMCP(soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>::iterator const& jet, 
                      aod::JetParticles const&, 
                      aod::JetMcCollisions const&, 
                      soa::Filtered<aod::JetCollisionsMCD> const& collisions)
  {
    if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
      return;
    }
    if (!isAcceptedJet<aod::JetParticles>(jet)) {
      return;
    }
    if (checkMcCollisionIsMatched) {
      auto collisionspermcpjet = collisions.sliceBy(CollisionsPerMCPCollision, jet.mcCollisionId());
      if (collisionspermcpjet.size() >= 1 && jetderiveddatautilities::selectCollision(collisionspermcpjet.begin(), eventSelection)) {
        fillMCPHistograms(jet);
      }
    } else {
      fillMCPHistograms(jet);
    }
  }
  PROCESS_SWITCH(JetSpectraChargedTask, processJetsMCP, "jet finder QA mcp", false);

  void processJetsMCPWeighted(soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetEventWeights>::iterator const& jet, 
                              aod::JetParticles const&, 
                              aod::JetMcCollisions const&, 
                              soa::Filtered<aod::JetCollisionsMCD> const& collisions)
  {
    if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
      return;
    }
    if (!isAcceptedJet<aod::JetParticles>(jet)) {
      return;
    }
    if (checkMcCollisionIsMatched) {
      auto collisionspermcpjet = collisions.sliceBy(CollisionsPerMCPCollision, jet.mcCollisionId());
      if (collisionspermcpjet.size() >= 1) {
        fillMCPHistograms(jet, jet.eventWeight());
      }
    } else {
      fillMCPHistograms(jet, jet.eventWeight());
    }                                        
  }
  PROCESS_SWITCH(JetSpectraChargedTask, processJetsMCPWeighted, "jet finder QA mcp with weighted events", false);

  void processMcparticles(aod::JetMcCollision const& mcCollision,
                        soa::SmallGroups<aod::JetCollisionsMCD> const& Collisions,
                        soa::Filtered<aod::JetParticles> const& mcparticles)
  {
      registry.fill(HIST("h_mccollisions_part"), 0.5);
      registry.fill(HIST("h2_centrality_mccollisions_part"), Collisions.begin().centrality(), 0.5);

      if (abs(mcCollision.posZ()) > vertexZCut) {
          return;
      }

      if (!abs(Collisions.size() < 1)) {
          return;
      }

      registry.fill(HIST("h_mccollisions_part"), 1.5);
      registry.fill(HIST("h2_centrality_mccollisions_part"), Collisions.begin().centrality(), 1.5);

      fillParticleHistograms(Collisions.begin(), mcparticles);
  }
  PROCESS_SWITCH(JetSpectraChargedTask, processMcparticles, "QA for charged mc particles", false);
    
  void processMcparticlesWeighted(aod::JetMcCollision const& mcCollision,
                                  soa::SmallGroups<aod::JetCollisionsMCD> const& Collisions,
                                  soa::Filtered<aod::JetParticles> const& mcparticles)
  {
      float eventweight = mcCollision.weight();
      registry.fill(HIST("h_mccollisions_part"), 0.5);
      registry.fill(HIST("h_mccollisions_part_weighted"), 0.5, eventweight);
      if(abs(mcCollision.posZ() > vertexZCut)){
        return;
      }
      if(!abs(Collisions.size() < 1)){
        return;
      }
      registry.fill(HIST("h_mccollisions_part"), 1.5);
      registry.fill(HIST("h_mccollisions_part_weighted"), 1.5, eventweight);
      fillParticleHistograms(Collisions.begin(), mcparticles, eventweight);
  }
  PROCESS_SWITCH(JetSpectraChargedTask, processMcparticlesWeighted, "QA for charged mc particles", false);

  void processRho(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision, soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    int nTracks = 0;
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        nTracks++;
      }
    }
    registry.fill(HIST("h2_centrality_ntracks"), collision.centrality(), nTracks);
    registry.fill(HIST("h2_ntracks_rho"), nTracks, collision.rho());
    registry.fill(HIST("h2_centrality_rho"), collision.centrality(), collision.rho());
  }
  PROCESS_SWITCH(JetSpectraChargedTask, processRho, "QA for rho-area subtracted jets", false);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { 
  return WorkflowSpec{adaptAnalysisTask<JetSpectraChargedTask>(cfgc, TaskName{"jet-charged-spectra"})}; 
  }
