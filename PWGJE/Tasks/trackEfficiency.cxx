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
  using myJetTracksMCD = o2::soa::Join<JetTracks, o2::aod::JMcTrackLbs>;

  HistogramRegistry registry;

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections; other option: uniformTracks"};
  Configurable<double> etaAcceptance{"etaAcceptance", 0.9, "eta acceptance; would be good to draw it from the tracksel automatically"};
  Configurable<bool> checkProducedByGen{"checkProducedByGen", true, "0: doesn't check mcparticle.producedByGenerator() - 1: checks particle.producedByGenerator()"};
  Configurable<bool> checkPrimaryPart{"checkPrimaryPart", true, "0: doesn't check mcparticle.isPhysicalPrimary() - 1: checks particle.isPhysicalPrimary()"};
  Configurable<bool> checkHepMCStatusCodeFinalState{"checkHepMCStatusCodeFinalState", false, "0: doesn't make sure that mcparticle.getHepMCStatusCode() is equal to 1 - 1: checks particle.getHepMCStatusCode() == 1 ; see https://pythia.org/latest-manual/ParticleProperties.html"};
  Configurable<bool> checkCentrality{"checkCentrality", false, ""};
  Configurable<double> centralityMin{"centralityMin", -999, ""};
  Configurable<double> centralityMax{"centralityMax", 999, ""};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> trackDcaZmax{"trackDcaZmax", 2, "maximum dcaZ to PV acceptance for tracks"};

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

  double maxDcaXYPtDep(double pt) // global track default: see https://github.com/AliceO2Group/O2Physics/blob/master/Common/Core/TrackSelectionDefaults.cxx
  {
    return 0.0105f + 0.0350f / pow(pt, 1.1f); // could also just re-add this passedDCAxy cut to the uniformTrack selection; studying this cut myself sounds tricky as there are 3 variables
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
    registry.get<TH1>(HIST("hMcPartCutsCounts"))->GetXaxis()->SetBinLabel(5, "producedByGen");
    registry.get<TH1>(HIST("hMcPartCutsCounts"))->GetXaxis()->SetBinLabel(6, "isFinalState");

    registry.add("hTrackCutsCounts", "Track cuts count checks", {HistType::kTH1F, {{10, 0., 10.}}});
    registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(1, "allTracksInAllColl");
    registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(2, "hasMcParticle");
    registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(3, "trackSel");
    registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(4, "mcPartIsPrimary");
    registry.get<TH1>(HIST("hTrackCutsCounts"))->GetXaxis()->SetBinLabel(5, "etaAcc");

    registry.add("h3_track_pt_track_eta_track_phi_associatedtrackSelCollSplit", "#it{p}_{T, associatedSplitTrack} (GeV/#it{c}); #eta_{associatedSplitTrack}; #phi_{associatedSplitTrack}", {HistType::kTH3F, {{500, 0., 10.}, {100, -1.0, 1.0}, {400, -1.0, 7.}}});

    registry.add("h3_track_pt_track_eta_track_phi_mcparticles", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {{500, 0., 10.}, {100, -1.0, 1.0}, {400, -1.0, 7.}}});
    registry.add("h3_track_pt_track_eta_track_phi_associatedtrackSelColl", "#it{p}_{T, associatedTrack} (GeV/#it{c}); #eta_{associatedTrack}; #phi_{associatedTrack}", {HistType::kTH3F, {{500, 0., 10.}, {100, -1.0, 1.0}, {400, -1.0, 7.}}});
    registry.add("h3_track_pt_track_eta_track_phi_associatedtrackNonSelColl", "#it{p}_{T, associatedTrack} (GeV/#it{c}); #eta_{associatedTrack}; #phi_{associatedTrack}", {HistType::kTH3F, {{500, 0., 10.}, {100, -1.0, 1.0}, {400, -1.0, 7.}}});

    registry.add("h3_track_pt_track_eta_track_phi_associatedtrackSelCollNonPrimary", "#it{p}_{T, associatedTrack} (GeV/#it{c}); #eta_{associatedTrack}; #phi_{associatedTrack}", {HistType::kTH3F, {{500, 0., 10.}, {100, -1.0, 1.0}, {400, -1.0, 7.}}});
  }

  Preslice<myJetTracksMCD> tracksPerJCollision = o2::aod::jtrack::collisionId;

  void process(JetMcCollision const& mcCollision, 
              soa::SmallGroups<JetCollisionsMCD> const& collisions, //smallgroups gives only the collisions associated to the current mccollision, thanks to the mccollisionlabel pre-integrated in jetcollisionsmcd
              soa::Join<myJetTracksMCD, aod::JTrackExtras> const& jetTracks,
              JetParticlesWithOriginal const& jMcParticles,
              aod::McParticles const&)
  {
    // missing:
    //   * constexpr auto hasCentrality = CollisionMCRecTableCentFT0C::template contains<aod::CentFT0Cs>();
    //           if constexpr (hasCentrality) {
    //   * dividing in centrality bins
    // I should maybe introduce the sel8 cuts on the collisoins (reco, but what about mccoll? maybe not htat way included in efficiency)



    registry.fill(HIST("hMcCollCutsCounts"), 0.5);
    if (!(mcCollision.posZ() < vertexZCut)) {
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 1.5);
    if (collisions.size() < 1) { // Skipping MC events that have no reconstructed collisions
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 2.5);
    if (collisions.size() > 1) { // Skipping MC events that are split
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 3.5);

    bool hasSel8Coll = false;
    bool centralityCheck = false;
    for (auto& collision : collisions){
      if (jetderiveddatautilities::selectCollision(collision, eventSelection)) { // Skipping MC events whose only reco collision isn't sel8
        hasSel8Coll = true; // should I actually put this cut after filling the denominator? for now before to try and get same results as Abhi
      }
      if (!checkCentrality || ((centralityMin < collision.centrality()) && (collision.centrality() < centralityMax))) {
        centralityCheck = true;
      }
    }
    if (!hasSel8Coll) { // checks that at least one of the collisions associated with this mcColl is selected
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 4.5);
    if (!centralityCheck) {
      return;
    }
    registry.fill(HIST("hMcCollCutsCounts"), 5.5);
    
    for (auto& jMcParticle : jMcParticles) {
      auto mcParticle = jMcParticle.mcParticle_as<aod::McParticles>();

      registry.fill(HIST("hMcPartCutsCounts"), 0.5); // allPartsInSelMcColl

      if (checkPrimaryPart && !jMcParticle.isPhysicalPrimary()) { // global tracks should be mostly primaries
        continue;
      }
      registry.fill(HIST("hMcPartCutsCounts"), 1.5); // isPrimary

      if (!(abs(jMcParticle.eta()) < etaAcceptance)) {
        continue;
      }
      registry.fill(HIST("hMcPartCutsCounts"), 2.5); // etaAccept

      if (!isChargedParticle(jMcParticle.pdgCode())) {
        continue;
      }
      registry.fill(HIST("hMcPartCutsCounts"), 3.5); // isCharged

      if (checkProducedByGen && !mcParticle.producedByGenerator()) {
        continue;
      }
      registry.fill(HIST("hMcPartCutsCounts"), 4.5); // producedByGen
      if (checkHepMCStatusCodeFinalState && jMcParticle.getHepMCStatusCode() != 1) { // a priori shouldn't be checked as we don't care about secondaries given they pollute the jet; to be discussed in a jet meeting to be sure
        continue;
      }
      registry.fill(HIST("hMcPartCutsCounts"), 5.5); // isFinalState
      
      registry.fill(HIST("h3_track_pt_track_eta_track_phi_mcparticles"), jMcParticle.pt(), jMcParticle.eta(), jMcParticle.phi());
    }

    std::vector<int> seenMcParticlesVector; // is reset every mc collision

    for (auto& collision : collisions){ // note: only looks at the only collision of the mcCollision as it is checked that collisions.size() == 1; for loop only present because one might want to change the splitting condition later on 

      auto collTracks = jetTracks.sliceBy(tracksPerJCollision, collision.globalIndex());
      for (auto& track : collTracks) {
        registry.fill(HIST("hTrackCutsCounts"), 0.5);

        if (!track.has_mcParticle()) {
          continue;
        }
        registry.fill(HIST("hTrackCutsCounts"), 1.5);

        if (!jetderiveddatautilities::selectTrack(track, trackSelection) || (trackSelections->compare("uniformTracks") == 1 && (track.dcaXY() > maxDcaXYPtDep(track.pt()) || track.dcaZ() > trackDcaZmax))) { /// if track selection is uniformTrack, I need to add dcaXY and dcaZ cuts as they aren't in the selection so that they can be studied here
          continue;
        }
        registry.fill(HIST("hTrackCutsCounts"), 2.5);

        if (!jetderiveddatautilities::selectCollision(collision, eventSelection) || !(collision.posZ() < vertexZCut)) {
          registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrackNonSelColl"), track.mcParticle_as<JetParticlesWithOriginal>().pt(), track.mcParticle_as<JetParticlesWithOriginal>().eta(), track.mcParticle_as<JetParticlesWithOriginal>().phi());
        }

        auto jMcParticleFromTrack = track.mcParticle_as<JetParticlesWithOriginal>();
        // auto mcParticleFromTrack = jMcParticleFromTrack.mcParticle_as<aod::McParticles>();
        if (!jMcParticleFromTrack.isPhysicalPrimary()) {
          registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrackSelCollNonPrimary"), track.mcParticle_as<JetParticlesWithOriginal>().pt(), track.mcParticle_as<JetParticlesWithOriginal>().eta(), track.mcParticle_as<JetParticlesWithOriginal>().phi());
          continue;
        }


        registry.fill(HIST("hTrackCutsCounts"), 3.5);


        registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrackSelColl"), track.mcParticle_as<JetParticlesWithOriginal>().pt(), track.mcParticle_as<JetParticlesWithOriginal>().eta(), track.mcParticle_as<JetParticlesWithOriginal>().phi());

        if (std::find(seenMcParticlesVector.begin(), seenMcParticlesVector.end(), track.mcParticle_as<JetParticlesWithOriginal>().globalIndex()) != seenMcParticlesVector.end()) { 
          registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrackSelCollSplit"), track.mcParticle_as<JetParticlesWithOriginal>().pt(), track.mcParticle_as<JetParticlesWithOriginal>().eta(), track.mcParticle_as<JetParticlesWithOriginal>().phi());
        } else {
          seenMcParticlesVector.push_back(track.mcParticle_as<JetParticlesWithOriginal>().globalIndex());
        }

        if (abs(track.mcParticle_as<JetParticlesWithOriginal>().eta()) < etaAcceptance) { // not actually applied here but it will give an idea of what will be done in the post processing
          registry.fill(HIST("hTrackCutsCounts"), 4.5);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<TrackEfficiencyJets>(cfgc, TaskName{"track-efficiency"})}; }
