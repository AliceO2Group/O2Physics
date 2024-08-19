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

// track efficiency task (global tracks)
//
/// \author Aimeric Landou <aimeric.landou@cern.ch>

#include <cmath>
#include <TRandom3.h>
#include <TMath.h>

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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct TrackEfficiencyJets {
  Service<o2::framework::O2DatabasePDG> pdg;

  using JetParticlesWithOriginal = soa::Join<JetParticles, aod::JMcParticlePIs>;

  HistogramRegistry registry;

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections; other option: uniformTracks"};

  // Tracking efficiency process function configurables:
  Configurable<bool> checkPrimaryPart{"checkPrimaryPart", true, "0: doesn't check mcparticle.isPhysicalPrimary() - 1: checks particle.isPhysicalPrimary()"};
  Configurable<bool> checkCentrality{"checkCentrality", false, ""};
  Configurable<int> acceptSplitCollisions{"acceptSplitCollisions", 0, "0: only look at mcCollisions that are not split; 1: accept split mcCollisions, 2: accept split mcCollisions but only look at the first reco collision associated with it"};
  Configurable<float> trackEtaAcceptanceCountQA{"trackEtaAcceptanceCountQA", 0.9, "eta acceptance"}; // removed from actual cuts for now because all the histograms have an eta axis
  Configurable<float> centralityMin{"centralityMin", -999, ""};
  Configurable<float> centralityMax{"centralityMax", 999, ""};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> trackDcaZmax{"trackDcaZmax", 99, "additional cut on dcaZ to PV for tracks; uniformTracks in particular don't cut on this at all"};

  // Track QA process function configurables:
  Configurable<float> trackQAEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks in the processTracks QA"};
  Configurable<float> trackQAEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks in the processTracks QA"};
  Configurable<float> trackQAPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks in the processTracks QA"};
  Configurable<float> trackQAPtMax{"trackPtMax", 100.0, "maximum pT acceptance for tracks in the processTracks QA"};

  int eventSelection = -1;
  int trackSelection = -1;

  bool isChargedParticle(int code)
  {
    auto p = pdg->GetParticle(code);
    auto charge = 0.;
    if (p != nullptr) {
      charge = p->Charge();
    }
    return std::abs(charge) >= 3.;
  }

  template <typename T, typename U>
  void fillTrackHistograms(T const& collision, U const& tracks, float weight = 1.0)
  {
    for (auto const& track : tracks) {
      if (!(jetderiveddatautilities::selectTrack(track, trackSelection) && jetderiveddatautilities::selectTrackDcaZ(track, trackDcaZmax))) {
        continue;
      }
      registry.fill(HIST("h2_centrality_track_pt"), collision.centrality(), track.pt(), weight);
      registry.fill(HIST("h2_centrality_track_eta"), collision.centrality(), track.eta(), weight);
      registry.fill(HIST("h2_centrality_track_phi"), collision.centrality(), track.phi(), weight);
      registry.fill(HIST("h2_centrality_track_energy"), collision.centrality(), track.energy(), weight);
      registry.fill(HIST("h2_track_pt_track_sigma1overpt"), track.pt(), track.sigma1Pt(), weight);
      registry.fill(HIST("h2_track_pt_track_sigmapt"), track.pt(), track.sigma1Pt() * track.pt(), weight);
      registry.fill(HIST("h2_track_pt_high_track_sigma1overpt"), track.pt(), track.sigma1Pt(), weight);
      registry.fill(HIST("h2_track_pt_high_track_sigmapt"), track.pt(), track.sigma1Pt() * track.pt(), weight);
    }
  }

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    if (doprocessEFficiencyPurity) {

      registry.add("hMcCollCutsCounts", "McColl cuts count checks", {HistType::kTH1F, {{10, 0., 10.}}});
      registry.get<TH1>(HIST("hMcCollCutsCounts"))->GetXaxis()->SetBinLabel(1, "allMcColl");
      registry.get<TH1>(HIST("hMcCollCutsCounts"))->GetXaxis()->SetBinLabel(2, "vertexZ");
      registry.get<TH1>(HIST("hMcCollCutsCounts"))->GetXaxis()->SetBinLabel(3, "noRecoColl");
      registry.get<TH1>(HIST("hMcCollCutsCounts"))->GetXaxis()->SetBinLabel(4, "splitColl");
      registry.get<TH1>(HIST("hMcCollCutsCounts"))->GetXaxis()->SetBinLabel(5, "recoCollEvtSel");
      registry.get<TH1>(HIST("hMcCollCutsCounts"))->GetXaxis()->SetBinLabel(6, "centralityCut");

      registry.add("hMcPartCutsCounts", "McPart cuts count checks", {HistType::kTH1F, {{10, 0., 10.}}});
      registry.get<TH1>(HIST("hMcPartCutsCounts"))->GetXaxis()->SetBinLabel(1, "allPartsInSelMcColl");
      registry.get<TH1>(HIST("hMcPartCutsCounts"))->GetXaxis()->SetBinLabel(2, "isCharged");
      registry.get<TH1>(HIST("hMcPartCutsCounts"))->GetXaxis()->SetBinLabel(3, "isPrimary");
      registry.get<TH1>(HIST("hMcPartCutsCounts"))->GetXaxis()->SetBinLabel(4, "etaAccept"); // not actually applied here but it will give an idea of what will be done in the post processing

      registry.add("hTrackCutsCounts", "Track cuts count checks", {HistType::kTH1F, {{10, 0., 10.}}});
      registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(1, "allTracksInSelColl");
      registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(2, "trackSel");
      registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(3, "hasMcParticle");
      registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(4, "mcPartIsPrimary");
      registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(5, "etaAcc"); // not actually applied here but it will give an idea of what will be done in the post processing

      AxisSpec ptAxis = {500, 0., 10., "#it{p}_{T} (GeV/#it{c})"};
      AxisSpec ptAxisHigh = {18, 10., 100., "#it{p}_{T} (GeV/#it{c})"};
      AxisSpec etaAxis = {100, -1.0, 1.0, "#eta"};
      AxisSpec phiAxis = {200, -1.0, 7., "#phi"};

      // ptAxisLow
      registry.add("h3_particle_pt_particle_eta_particle_phi_mcpartofinterest", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxis, etaAxis, phiAxis}});
      registry.add("h3_particle_pt_particle_eta_particle_phi_mcpart_nonprimary", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxis, etaAxis, phiAxis}});

      registry.add("h3_track_pt_track_eta_track_phi_nonassociatedtrack", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {ptAxis, etaAxis, phiAxis}});
      registry.add("h3_track_pt_track_eta_track_phi_associatedtrack_primary", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {ptAxis, etaAxis, phiAxis}});
      registry.add("h3_track_pt_track_eta_track_phi_associatedtrack_nonprimary", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {ptAxis, etaAxis, phiAxis}});
      registry.add("h3_track_pt_track_eta_track_phi_associatedtrack_split_primary", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {ptAxis, etaAxis, phiAxis}});
      registry.add("h3_track_pt_track_eta_track_phi_associatedtrack_split_nonprimary", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {ptAxis, etaAxis, phiAxis}});

      registry.add("h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxis, etaAxis, phiAxis}});
      registry.add("h3_particle_pt_particle_eta_particle_phi_associatedtrack_nonprimary", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxis, etaAxis, phiAxis}});
      registry.add("h3_particle_pt_particle_eta_particle_phi_associatedtrack_split_primary", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxis, etaAxis, phiAxis}});
      registry.add("h3_particle_pt_particle_eta_particle_phi_associatedtrack_split_nonprimary", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxis, etaAxis, phiAxis}});

      registry.add("h2_particle_pt_track_pt_residual_associatedtrack_primary", "(#it{p}_{T, mcpart} - #it{p}_{T, track}) / #it{p}_{T, mcpart}; #it{p}_{T, mcpart} (GeV/#it{c})", {HistType::kTH2F, {ptAxis, {200, -1., 1.}}});

      // ptAxisHigh
      registry.add("h3_particle_pt_high_particle_eta_particle_phi_mcpartofinterest", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxisHigh, etaAxis, phiAxis}});

      registry.add("h3_track_pt_high_track_eta_track_phi_nonassociatedtrack", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {ptAxisHigh, etaAxis, phiAxis}});
      registry.add("h3_track_pt_high_track_eta_track_phi_associatedtrack_primary", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {ptAxisHigh, etaAxis, phiAxis}});
      registry.add("h3_track_pt_high_track_eta_track_phi_associatedtrack_nonprimary", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {ptAxisHigh, etaAxis, phiAxis}});
      registry.add("h3_track_pt_high_track_eta_track_phi_associatedtrack_split_primary", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {ptAxisHigh, etaAxis, phiAxis}});
      registry.add("h3_track_pt_high_track_eta_track_phi_associatedtrack_split_nonprimary", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {ptAxisHigh, etaAxis, phiAxis}});

      registry.add("h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_primary", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxisHigh, etaAxis, phiAxis}});
      registry.add("h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_nonprimary", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxisHigh, etaAxis, phiAxis}});
      registry.add("h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_split_primary", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxisHigh, etaAxis, phiAxis}});
      registry.add("h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_split_nonprimary", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {ptAxisHigh, etaAxis, phiAxis}});

      registry.add("h2_particle_pt_high_track_pt_high_residual_associatedtrack_primary", "(#it{p}_{T, mcpart} - #it{p}_{T, track}) / #it{p}_{T, mcpart}; #it{p}_{T, mcpart} (GeV/#it{c})", {HistType::kTH2F, {ptAxisHigh, {200, -1., 1.}}});
    }

    if (doprocessTracks || doprocessTracksWeighted) {
      AxisSpec centAxis = {101, 0., 101., "centrality (%)"};
      registry.add("h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      registry.add("h2_centrality_collisions", "centrality vs collisions; centrality; collisions", {HistType::kTH2F, {centAxis, {4, 0.0, 4.0}}});
      registry.add("h2_centrality_track_pt", "centrality vs track pT; centrality; #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {centAxis, {200, 0., 200.}}});
      registry.add("h2_centrality_track_eta", "centrality vs track #eta; centrality; #eta_{track}", {HistType::kTH2F, {centAxis, {100, -1.0, 1.0}}});
      registry.add("h2_centrality_track_phi", "centrality vs track #varphi; centrality; #varphi_{track}", {HistType::kTH2F, {centAxis, {160, -1.0, 7.}}});
      registry.add("h2_centrality_track_energy", "centrality vs track energy; centrality; Energy GeV", {HistType::kTH2F, {centAxis, {100, 0.0, 100.0}}});
      registry.add("h2_track_pt_track_sigmapt", "#sigma(#it{p}_{T})/#it{p}_{T}; #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 10.}, {100000, 0.0, 100.0}}});
      registry.add("h2_track_pt_high_track_sigmapt", "#sigma(#it{p}_{T})/#it{p}_{T}; #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{90, 10., 100.}, {100000, 0.0, 100.0}}});
      registry.add("h2_track_pt_track_sigma1overpt", "#sigma(1/#it{p}_{T}); #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 10.}, {1000, 0.0, 10.0}}});
      registry.add("h2_track_pt_high_track_sigma1overpt", "#sigma(1/#it{p}_{T}); #it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{90, 10., 100.}, {1000, 0.0, 10.0}}});
      if (doprocessTracksWeighted) {
        registry.add("h_collisions_weighted", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}});
      }
    }
  }

  Preslice<JetTracksMCD> tracksPerJCollision = o2::aod::jtrack::collisionId;

  // filters for processTracks QA functions only:
  Filter trackCuts = (aod::jtrack::pt >= trackQAPtMin && aod::jtrack::pt < trackQAPtMax && aod::jtrack::eta > trackQAEtaMin && aod::jtrack::eta < trackQAEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);

  void processEFficiencyPurity(JetMcCollision const& mcCollision,
                               soa::SmallGroups<JetCollisionsMCD> const& collisions, // smallgroups gives only the collisions associated to the current mccollision, thanks to the mccollisionlabel pre-integrated in jetcollisionsmcd
                               soa::Join<JetTracksMCD, aod::JTrackExtras> const& jetTracks,
                               JetParticlesWithOriginal const& jMcParticles)
  {
    // missing:
    //   * constexpr auto hasCentrality = CollisionMCRecTableCentFT0C::template contains<aod::CentFT0Cs>();
    //           if constexpr (hasCentrality) {
    //   * dividing in centrality bins
    // I should maybe introduce the sel8 cuts on the collisoins (reco, but what about mccoll? maybe not htat way included in efficiency)

    registry.fill(HIST("hMcCollCutsCounts"), 0.5); // all mcCollisions

    if (!(abs(mcCollision.posZ()) < vertexZCut)) {
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 1.5); // mcCollision.posZ() condition

    if (collisions.size() < 1) {
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 2.5); // mcCollisions with at least one reconstructed collision

    if (acceptSplitCollisions == 0 && collisions.size() > 1) {
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 3.5); // split mcCollisions condition

    bool hasSel8Coll = false;
    bool centralityCheck = false;
    if (acceptSplitCollisions == 2) {                                                     // check only that the first reconstructed collision passes the check
      if (jetderiveddatautilities::selectCollision(collisions.begin(), eventSelection)) { // Skipping MC events that have not a single selected reconstructed collision ; effect unclear if mcColl is split
        hasSel8Coll = true;
      }
      if (!checkCentrality || ((centralityMin < collisions.begin().centrality()) && (collisions.begin().centrality() < centralityMax))) { // effect unclear if mcColl is split
        centralityCheck = true;
      }
    } else { // check that at least one of the reconstructed collisions passes the checks
      for (auto& collision : collisions) {
        if (jetderiveddatautilities::selectCollision(collision, eventSelection)) { // Skipping MC events that have not a single selected reconstructed collision ; effect unclear if mcColl is split
          hasSel8Coll = true;
        }
        if (!checkCentrality || ((centralityMin < collision.centrality()) && (collision.centrality() < centralityMax))) { // effect unclear if mcColl is split
          centralityCheck = true;
        }
      }
    }
    if (!hasSel8Coll) {
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 4.5); // at least one of the reconstructed collisions associated with this mcCollision is selected

    if (!centralityCheck) {
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 5.5); // at least one of the reconstructed collisions associated with this mcCollision is selected with regard to centrality

    for (auto& jMcParticle : jMcParticles) {
      registry.fill(HIST("hMcPartCutsCounts"), 0.5); // allPartsInSelMcColl

      if (!isChargedParticle(jMcParticle.pdgCode())) {
        continue;
      }
      registry.fill(HIST("hMcPartCutsCounts"), 1.5); // isCharged

      registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_mcpart_nonprimary"), jMcParticle.pt(), jMcParticle.eta(), jMcParticle.phi());

      if (checkPrimaryPart && !jMcParticle.isPhysicalPrimary()) { // global tracks should be mostly primaries
        continue;
      }
      registry.fill(HIST("hMcPartCutsCounts"), 2.5); // isPrimary

      registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"), jMcParticle.pt(), jMcParticle.eta(), jMcParticle.phi());

      registry.fill(HIST("h3_particle_pt_high_particle_eta_particle_phi_mcpartofinterest"), jMcParticle.pt(), jMcParticle.eta(), jMcParticle.phi());

      if ((abs(jMcParticle.eta()) < trackEtaAcceptanceCountQA)) { // removed from actual cuts for now because all the histograms have an eta axis
        registry.fill(HIST("hMcPartCutsCounts"), 3.5);            // etaAccept // not actually applied here but it will give an idea of what will be done in the post processing
      }
    }

    std::vector<int> seenMcParticlesVector; // is reset every mc collision

    int splitCollCounter = 0;
    for (auto& collision : collisions) {
      splitCollCounter++;
      if (acceptSplitCollisions == 2 && splitCollCounter > 1) {
        return;
      }

      if (!jetderiveddatautilities::selectCollision(collision, eventSelection) || !(abs(collision.posZ()) < vertexZCut)) {
        continue;
      }

      auto collTracks = jetTracks.sliceBy(tracksPerJCollision, collision.globalIndex());
      for (auto& track : collTracks) {
        registry.fill(HIST("hTrackCutsCounts"), 0.5);

        if (!(jetderiveddatautilities::selectTrack(track, trackSelection) && jetderiveddatautilities::selectTrackDcaZ(track, trackDcaZmax))) { // if track selection is uniformTrack, dcaZ cuts need to be added as they aren't in the selection so that they can be studied here
          continue;
        }
        registry.fill(HIST("hTrackCutsCounts"), 1.5);

        if (!track.has_mcParticle()) {
          registry.fill(HIST("h3_track_pt_track_eta_track_phi_nonassociatedtrack"), track.pt(), track.eta(), track.phi());

          registry.fill(HIST("h3_track_pt_high_track_eta_track_phi_nonassociatedtrack"), track.pt(), track.eta(), track.phi());
          continue;
        }
        registry.fill(HIST("hTrackCutsCounts"), 2.5);

        auto jMcParticleFromTrack = track.mcParticle_as<JetParticlesWithOriginal>();
        if (!jMcParticleFromTrack.isPhysicalPrimary()) {
          registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrack_nonprimary"), track.pt(), track.eta(), track.phi());
          registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_associatedtrack_nonprimary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi());

          registry.fill(HIST("h3_track_pt_high_track_eta_track_phi_associatedtrack_nonprimary"), track.pt(), track.eta(), track.phi());
          registry.fill(HIST("h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_nonprimary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi());

          if (std::find(seenMcParticlesVector.begin(), seenMcParticlesVector.end(), jMcParticleFromTrack.globalIndex()) != seenMcParticlesVector.end()) {
            registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrack_split_nonprimary"), track.pt(), track.eta(), track.phi());
            registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_associatedtrack_split_nonprimary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi());

            registry.fill(HIST("h3_track_pt_high_track_eta_track_phi_associatedtrack_split_nonprimary"), track.pt(), track.eta(), track.phi());
            registry.fill(HIST("h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_split_nonprimary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi());
          } else {
            seenMcParticlesVector.push_back(jMcParticleFromTrack.globalIndex());
          }

          continue;
        }

        registry.fill(HIST("hTrackCutsCounts"), 3.5);

        registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrack_primary"), track.pt(), track.eta(), track.phi());
        registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi());
        registry.fill(HIST("h2_particle_pt_track_pt_residual_associatedtrack_primary"), jMcParticleFromTrack.pt(), (jMcParticleFromTrack.pt() - track.pt()) / jMcParticleFromTrack.pt());

        registry.fill(HIST("h3_track_pt_high_track_eta_track_phi_associatedtrack_primary"), track.pt(), track.eta(), track.phi());
        registry.fill(HIST("h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_primary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi());
        registry.fill(HIST("h2_particle_pt_high_track_pt_high_residual_associatedtrack_primary"), jMcParticleFromTrack.pt(), (jMcParticleFromTrack.pt() - track.pt()) / jMcParticleFromTrack.pt());

        if (std::find(seenMcParticlesVector.begin(), seenMcParticlesVector.end(), jMcParticleFromTrack.globalIndex()) != seenMcParticlesVector.end()) {
          registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrack_split_primary"), track.pt(), track.eta(), track.phi());
          registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_associatedtrack_split_primary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi());

          registry.fill(HIST("h3_track_pt_high_track_eta_track_phi_associatedtrack_split_primary"), track.pt(), track.eta(), track.phi());
          registry.fill(HIST("h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_split_primary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi());
        } else {
          seenMcParticlesVector.push_back(jMcParticleFromTrack.globalIndex());
        }

        if (abs(jMcParticleFromTrack.eta()) < trackEtaAcceptanceCountQA) { // not actually applied here but it will give an idea of what will be done in the post processing
          registry.fill(HIST("hTrackCutsCounts"), 4.5);
        }
      }
    }
  }
  PROCESS_SWITCH(TrackEfficiencyJets, processEFficiencyPurity, "Histograms for efficiency and purity quantities", true);

  void processTracks(soa::Filtered<JetCollisions>::iterator const& collision,
                     soa::Filtered<soa::Join<JetTracks, aod::JTrackExtras>> const& tracks)
  {
    registry.fill(HIST("h_collisions"), 0.5);
    registry.fill(HIST("h2_centrality_collisions"), collision.centrality(), 0.5);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    registry.fill(HIST("h2_centrality_collisions"), collision.centrality(), 1.5);
    fillTrackHistograms(collision, tracks);
  }
  PROCESS_SWITCH(TrackEfficiencyJets, processTracks, "QA for charged tracks", false);

  void processTracksWeighted(soa::Join<JetCollisions, aod::JMcCollisionLbs>::iterator const& collision,
                             JetMcCollisions const&,
                             soa::Filtered<soa::Join<JetTracks, aod::JTrackExtras>> const& tracks)
  {
    float eventWeight = collision.mcCollision().weight();
    registry.fill(HIST("h_collisions"), 0.5);
    registry.fill(HIST("h_collisions_weighted"), 0.5, eventWeight);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    registry.fill(HIST("h_collisions_weighted"), 1.5, eventWeight);
    fillTrackHistograms(collision, tracks, eventWeight);
  }
  PROCESS_SWITCH(TrackEfficiencyJets, processTracksWeighted, "QA for charged tracks weighted", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<TrackEfficiencyJets>(cfgc, TaskName{"track-efficiency"})}; }
