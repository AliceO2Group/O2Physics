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

  Configurable<bool> checkPrimaryPart{"checkPrimaryPart", true, "0: doesn't check mcparticle.isPhysicalPrimary() - 1: checks particle.isPhysicalPrimary()"};
  Configurable<bool> checkCentrality{"checkCentrality", false, ""};
  Configurable<bool> splitMcCollOK{"splitMcCollOK", false, "false: only look at mcCollisions that are not split; true: accept split mcCollisions"};

  Configurable<float> trackEtaAcceptance{"trackEtaAcceptance", 0.9, "eta acceptance"};
  Configurable<float> centralityMin{"centralityMin", -999, ""};
  Configurable<float> centralityMax{"centralityMax", 999, ""};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> trackDcaZmax{"trackDcaZmax", -99, "additional cut on dcaZ to PV for tracks; uniformTracks in particular don't cut on this at all"};

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

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    registry.add("hMcCollCutsCounts", "McColl cuts count checks", {HistType::kTH1F, {{10, 0., 10.}}});
    registry.get<TH1>(HIST("hMcCollCutsCounts"))->GetXaxis()->SetBinLabel(1, "allMcColl");
    registry.get<TH1>(HIST("hMcCollCutsCounts"))->GetXaxis()->SetBinLabel(2, "vertexZ");
    registry.get<TH1>(HIST("hMcCollCutsCounts"))->GetXaxis()->SetBinLabel(3, "noRecoColl");
    registry.get<TH1>(HIST("hMcCollCutsCounts"))->GetXaxis()->SetBinLabel(4, "splitColl");
    registry.get<TH1>(HIST("hMcCollCutsCounts"))->GetXaxis()->SetBinLabel(5, "recoCollEvtSel");
    registry.get<TH1>(HIST("hMcCollCutsCounts"))->GetXaxis()->SetBinLabel(6, "centralityCut");

    registry.add("hMcPartCutsCounts", "McPart cuts count checks", {HistType::kTH1F, {{10, 0., 10.}}});
    registry.get<TH1>(HIST("hMcPartCutsCounts"))->GetXaxis()->SetBinLabel(1, "allPartsInSelMcColl");
    registry.get<TH1>(HIST("hMcPartCutsCounts"))->GetXaxis()->SetBinLabel(2, "isPrimary");
    registry.get<TH1>(HIST("hMcPartCutsCounts"))->GetXaxis()->SetBinLabel(3, "etaAccept");
    registry.get<TH1>(HIST("hMcPartCutsCounts"))->GetXaxis()->SetBinLabel(4, "isCharged");

    registry.add("hTrackCutsCounts", "Track cuts count checks", {HistType::kTH1F, {{10, 0., 10.}}});
    registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(1, "allTracksInSelColl");
    registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(2, "trackSel");
    registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(3, "hasMcParticle");
    registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(4, "mcPartIsPrimary");
    registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(5, "etaAcc");

    registry.add("h3_particle_pt_particle_eta_particle_phi_mcpartofinterest", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {{500, 0., 10.}, {100, -1.0, 1.0}, {400, -1.0, 7.}}});

    registry.add("h3_track_pt_track_eta_track_phi_nonassociatedtrack", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {{500, 0., 10.}, {100, -1.0, 1.0}, {400, -1.0, 7.}}});
    registry.add("h3_track_pt_track_eta_track_phi_associatedtrack_primary", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {{500, 0., 10.}, {100, -1.0, 1.0}, {400, -1.0, 7.}}});
    registry.add("h3_track_pt_track_eta_track_phi_associatedtrack_nonprimary", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {{500, 0., 10.}, {100, -1.0, 1.0}, {400, -1.0, 7.}}});
    registry.add("h3_track_pt_track_eta_track_phi_associatedtrack_split", "#it{p}_{T, track} (GeV/#it{c}); #eta_{track}; #phi_{track}", {HistType::kTH3F, {{500, 0., 10.}, {100, -1.0, 1.0}, {400, -1.0, 7.}}});

    registry.add("h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {{500, 0., 10.}, {100, -1.0, 1.0}, {400, -1.0, 7.}}});
    registry.add("h3_particle_pt_particle_eta_particle_phi_associatedtrack_nonprimary", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {{500, 0., 10.}, {100, -1.0, 1.0}, {400, -1.0, 7.}}});
    registry.add("h3_particle_pt_particle_eta_particle_phi_associatedtrack_split", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {{500, 0., 10.}, {100, -1.0, 1.0}, {400, -1.0, 7.}}});
  }

  Preslice<JetTracksMCD> tracksPerJCollision = o2::aod::jtrack::collisionId;

  void process(JetMcCollision const& mcCollision,
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

    if (!splitMcCollOK && collisions.size() > 1) {
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 3.5); // split mcCollisions condition

    bool hasSel8Coll = false;
    bool centralityCheck = false;
    for (auto& collision : collisions) {
      if (jetderiveddatautilities::selectCollision(collision, eventSelection)) { // Skipping MC events that have not a single selected reconstructed collision ; effect unclear if mcColl is split
        hasSel8Coll = true;
      }
      if (!checkCentrality || ((centralityMin < collision.centrality()) && (collision.centrality() < centralityMax))) { // effect unclear if mcColl is split
        centralityCheck = true;
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

      if (checkPrimaryPart && !jMcParticle.isPhysicalPrimary()) { // global tracks should be mostly primaries
        continue;
      }
      registry.fill(HIST("hMcPartCutsCounts"), 1.5); // isPrimary

      if (!(abs(jMcParticle.eta()) < trackEtaAcceptance)) {
        continue;
      }
      registry.fill(HIST("hMcPartCutsCounts"), 2.5); // etaAccept

      if (!isChargedParticle(jMcParticle.pdgCode())) {
        continue;
      }
      registry.fill(HIST("hMcPartCutsCounts"), 3.5); // isCharged

      registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"), jMcParticle.pt(), jMcParticle.eta(), jMcParticle.phi());
    }

    std::vector<int> seenMcParticlesVector; // is reset every mc collision

    for (auto& collision : collisions) {

      auto collTracks = jetTracks.sliceBy(tracksPerJCollision, collision.globalIndex());
      for (auto& track : collTracks) {
        if (!jetderiveddatautilities::selectCollision(collision, eventSelection) || !(abs(collision.posZ()) < vertexZCut)) { // selectCollision is mostly here for readability, as the code only looks at the collision associated to a mc collision for which this has already been checked ; collision.posZ() hasn't been checked yet though
          continue;
        }
        registry.fill(HIST("hTrackCutsCounts"), 0.5);

        if (!(jetderiveddatautilities::selectTrack(track, trackSelection) && jetderiveddatautilities::selectTrackDcaZ(track, trackDcaZmax))) { // if track selection is uniformTrack, dcaXY and dcaZ cuts need to be added as they aren't in the selection so that they can be studied here
          continue;
        }
        registry.fill(HIST("hTrackCutsCounts"), 1.5);

        if (!track.has_mcParticle()) {
          registry.fill(HIST("h3_track_pt_track_eta_track_phi_nonassociatedtrack"), track.pt(), track.eta(), track.phi());
          continue;
        }
        registry.fill(HIST("hTrackCutsCounts"), 2.5);

        auto jMcParticleFromTrack = track.mcParticle_as<JetParticlesWithOriginal>();
        if (!jMcParticleFromTrack.isPhysicalPrimary()) {
          registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrack_nonprimary"), track.pt(), track.eta(), track.phi());
          registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_associatedtrack_nonprimary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi());
          continue;
        }

        registry.fill(HIST("hTrackCutsCounts"), 3.5);

        registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrack_primary"), track.pt(), track.eta(), track.phi());
        registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi());

        if (std::find(seenMcParticlesVector.begin(), seenMcParticlesVector.end(), jMcParticleFromTrack.globalIndex()) != seenMcParticlesVector.end()) {
          registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrack_split"), track.pt(), track.eta(), track.phi());
          registry.fill(HIST("h3_particle_pt_particle_eta_particle_phi_associatedtrack_split"), jMcParticleFromTrack.pt(), jMcParticleFromTrack.eta(), jMcParticleFromTrack.phi());
        } else {
          seenMcParticlesVector.push_back(jMcParticleFromTrack.globalIndex());
        }

        if (abs(jMcParticleFromTrack.eta()) < trackEtaAcceptance) { // not actually applied here but it will give an idea of what will be done in the post processing
          registry.fill(HIST("hTrackCutsCounts"), 4.5);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<TrackEfficiencyJets>(cfgc, TaskName{"track-efficiency"})}; }