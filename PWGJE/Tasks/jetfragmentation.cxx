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

// jet trigger QA task
//
// Author: Gijs van Weelden
//

#include "TH1F.h"
#include "TTree.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/RunningWorkflowInfo.h"

#include "Common/DataModel/EventSelection.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetFinder.h"

#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using McTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>>;
using McTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>;
using McDJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>;
using McPJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>;

struct JetFragmentation {
  HistogramRegistry registry{"registry"};

  std::vector<int> pdgVector = {211, 321, 2212, 111, 130, 310, 311, 3122};
  std::vector<std::string> hadronVector = {"#pi^{#pm}", "#it{K}^{#pm}", "#it{p}^{#pm}", "#pi^{0}", "#it{K}^{0}_{L}", "#it{K}^{0}_{S}", "#it{K}^{0}", "#Lambda^{0}"};

  // Binning
  ConfigurableAxis binJetPt{"binJetPt", {40, 0.f, 200.f}, ""};
  ConfigurableAxis binEta{"binEta", {20, -1.f, 1.f}, ""};
  ConfigurableAxis binPhi{"binPhi", {18 * 8, 0.f, 2. * TMath::Pi()}, ""};
  ConfigurableAxis binZ{"binZ", {20, -5e-3f, 1.f + 5e-3f}, ""};
  ConfigurableAxis binJetR{"binJetR", {6, 0.05f, 0.65f}, ""};
  ConfigurableAxis binTrackPt{"binTrackPt", {200, 0.f, 100.f}, ""};
  ConfigurableAxis binPDG{"binPDG", {static_cast<double>(pdgVector.size()), -0.5f, static_cast<double>(pdgVector.size()) - 0.5f}, ""};
  ConfigurableAxis binVtxZ{"binVtxZ", {200, -20, 20}, ""};

  ConfigurableAxis binDiff{"binDiff", {50, -5.f, 5.f}, ""};
  ConfigurableAxis binRatio{"binRatio", {100, 0.f, 10.f}, ""};        // Ratio of pt, eta, phi
  ConfigurableAxis binMatchDist{"binMatchDist", {10, 0.f, 0.5f}, ""}; // Distance between matched jets

  ConfigurableAxis binCount{"binCount", {1, .5f, 1.5f}, ""};
  ConfigurableAxis jetCount{"jetCount", {50, -.5f, 49.5f}, ""};
  ConfigurableAxis trackCount{"trackCount", {100, -.5f, 99.5f}, ""};

  Preslice<McPJets> perMcPJet = aod::jet::mcCollisionId;

  void init(InitContext& initContext)
  {
    // Axes
    AxisSpec jetPtAxis = {binJetPt, "#it{p}_{T}^{ jet}"}; // Data
    AxisSpec etaAxis = {binEta, "#eta"};
    AxisSpec phiAxis = {binPhi, "#phi"};
    AxisSpec zAxis = {binZ, "#it{z}"};
    AxisSpec detJetPtAxis = {binJetPt, "#it{p}_{T}^{ jet, det}"}; // MC detector level
    AxisSpec detEtaAxis = {binEta, "#eta^{ jet, det}"};
    AxisSpec detPhiAxis = {binPhi, "#phi^{ jet, det}"};
    AxisSpec detZAxis = {binZ, "#it{z}^{ det}"};
    AxisSpec partJetPtAxis = {binJetPt, "#it{p}_{T}^{ jet, part}"}; // MC particle level
    AxisSpec partEtaAxis = {binEta, "#eta^{ jet, part}"};
    AxisSpec partPhiAxis = {binPhi, "#phi^{ jet, part}"};
    AxisSpec partZAxis = {binZ, "#it{z}^{ part}"};

    AxisSpec pdgAxis = {binPDG, ""};
    AxisSpec trackPtAxis = {binTrackPt, "#it{p}_{T}^{tr}"};
    AxisSpec diffAxis = {binDiff, ""};
    AxisSpec ratioAxis = {binRatio, ""};
    AxisSpec vtxZAxis = {binVtxZ, "Collision vertex z (cm)"};
    AxisSpec matchDistAxis = {binMatchDist, "#Delta"};

    // Data
    registry.add("collision/collisionVtxZ", "Collision vertex z (cm)", HistType::kTH1F, {binVtxZ});

    registry.add("tracks/trackPtEtaPhi", "trackPtEtaPhi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});

    registry.add("jets/jetPtEtaPhi", "Jet #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {jetPtAxis, etaAxis, phiAxis});
    registry.add("jets/jetPtTrackPt", "Jet #it{p}_{T}, track #it{p}_{T};#it{p}_{T}^{jet};#it{p}_{T}^{track}", HistType::kTH2F, {jetPtAxis, trackPtAxis});
    registry.add("jets/jetTrackPtEtaPhi", "Tracks in jets #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});
    registry.add("jets/jetPtTrackProj", "Jet #it{p}_{T}, #it{z}; #it{p}_{T}; #it{z}", HistType::kTH2F, {jetPtAxis, zAxis});
    registry.add("jets/jetPtFrag", "Jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}; #it{p}_{T}; #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2F, {jetPtAxis, zAxis});

    // MC particle level
    registry.add("collision/partCollisionVtxZ", "Collision vertex z (cm)", HistType::kTH1F, {binVtxZ});

    registry.add("tracks/partTrackPtEtaPhi", "partTrackPtEtaPhi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});

    registry.add("jets/partJetPtEtaPhi", "Particle level jet #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {partJetPtAxis, partEtaAxis, partPhiAxis});
    registry.add("jets/partJetPtTrackPt", "Particle level jet #it{p}_{T}, track #it{p}_{T};#it{p}_{T}^{jet};#it{p}_{T}^{track}", HistType::kTH2F, {partJetPtAxis, trackPtAxis});
    registry.add("jets/partJetTrackPtEtaPhi", "Particle level tracks in jets #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {trackPtAxis, partEtaAxis, partPhiAxis});
    registry.add("jets/partJetPtTrackProj", "Particle level jet #it{p}_{T}, #it{z}; #it{p}_{T}; #it{z}", HistType::kTH2F, {partJetPtAxis, partZAxis});
    registry.add("jets/partJetPtFrag", "Particle level jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}; #it{p}_{T}; #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2F, {partJetPtAxis, partZAxis});

    // MC detector level
    registry.add("collision/detCollisionVtxZ", "Collision vertex z (cm)", HistType::kTH1F, {binVtxZ});

    registry.add("tracks/detTrackPtEtaPhi", "detTrackPtEtaPhi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});

    registry.add("jets/detJetPtEtaPhi", "Detector level jet #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {detJetPtAxis, detEtaAxis, detPhiAxis});
    registry.add("jets/detJetPtTrackPt", "Detector level jet #it{p}_{T}, track #it{p}_{T};#it{p}_{T}^{jet};#it{p}_{T}^{track}", HistType::kTH2F, {detJetPtAxis, trackPtAxis});
    registry.add("jets/detJetTrackPtEtaPhi", "Detector level tracks in jets #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {trackPtAxis, detEtaAxis, detPhiAxis});
    registry.add("jets/detJetPtTrackProj", "Detector level jet #it{p}_{T}, #it{z}; #it{p}_{T}; #it{z}", HistType::kTH2F, {detJetPtAxis, detZAxis});
    registry.add("jets/detJetPtFrag", "Detector level jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}; #it{p}_{T}; #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2F, {detJetPtAxis, detZAxis});

    // MC particle-detector level matching
    registry.add("collision/matchCollisionVtxZ", "Collision vertex z (cm)", HistType::kTH1F, {binVtxZ});

    registry.add("tracks/matchDetTrackPtEtaPhi", "matchDetTrackPtEtaPhi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});
    registry.add("tracks/matchPartTrackPtEtaPhi", "matchPartTrackPtEtaPhi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});
    registry.add("tracks/matchDetTrackPtPartTrackPt", "matchDetTrackPtPartTrackPt", HistType::kTH2F, {trackPtAxis, trackPtAxis});
    registry.add("tracks/matchDetTrackEtaPartTrackEta", "matchDetTrackEtaPartTrackEta", HistType::kTH2F, {etaAxis, etaAxis});
    registry.add("tracks/matchDetTrackPhiPartTrackPhi", "matchDetTrackPhiPartTrackPhi", HistType::kTH2F, {phiAxis, phiAxis});
    registry.add("tracks/trackResolutionPt", "trackResolutionPt; #Delta #it{p}_{T}^{tr}", HistType::kTH2F, {trackPtAxis, diffAxis});
    registry.add("tracks/trackResolutionEta", "trackResolutionEta; #Delta #eta", HistType::kTH2F, {etaAxis, diffAxis});
    registry.add("tracks/trackResolutionPhi", "trackResolutionPhi; #Delta #phi", HistType::kTH2F, {phiAxis, diffAxis});

    // Detector level jets with a match
    registry.add("jets/matchDetJetPtEtaPhi", "Matched detector level jet #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {detJetPtAxis, detEtaAxis, detPhiAxis});
    registry.add("jets/matchDetJetPtTrackPt", "Matched detector level jet #it{p}_{T}, track #it{p}_{T};#it{p}_{T}^{jet};#it{p}_{T}^{track}", HistType::kTH2F, {detJetPtAxis, trackPtAxis});
    registry.add("jets/matchDetJetTrackPtEtaPhi", "Matched detector level tracks in jets #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {trackPtAxis, detEtaAxis, detPhiAxis});
    registry.add("jets/matchDetJetPtTrackProj", "Matched detector level jet #it{p}_{T}, #it{z}; #it{p}_{T}; #it{z}", HistType::kTH2F, {detJetPtAxis, detZAxis});
    registry.add("jets/matchDetJetPtFrag", "Matched detector level jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}; #it{p}_{T}; #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2F, {detJetPtAxis, detZAxis});
    // Particle level jets with a match
    registry.add("jets/matchPartJetPtEtaPhi", "Matched particle level jet #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {partJetPtAxis, partEtaAxis, partPhiAxis});
    registry.add("jets/matchPartJetPtTrackPt", "Matched particle level jet #it{p}_{T}, track #it{p}_{T};#it{p}_{T}^{jet};#it{p}_{T}^{track}", HistType::kTH2F, {partJetPtAxis, trackPtAxis});
    registry.add("jets/matchPartJetTrackPtEtaPhi", "Matched particle level tracks in jets #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {trackPtAxis, partEtaAxis, partPhiAxis});
    registry.add("jets/matchPartJetPtTrackProj", "Matched particle level jet #it{p}_{T}, #it{z}; #it{p}_{T}; #it{z}", HistType::kTH2F, {partJetPtAxis, partZAxis});
    registry.add("jets/matchPartJetPtFrag", "Matched particle level jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}; #it{p}_{T}; #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2F, {partJetPtAxis, partZAxis});
    // Combined information of matched jets
    registry.add("jets/matchDetJetPtPartJetPt", "matchDetJetPtPartJetPt; #it{p}_{T}^{jet, det}; #it{p}_{T}^{jet, part}", HistType::kTH2F, {detJetPtAxis, partJetPtAxis});
    registry.add("jets/matchDetJetEtaPartJetEta", "matchDetJetEtaPartJetEta; #eta^{jet, det}; #eta^{jet, part}", HistType::kTH2F, {detEtaAxis, partEtaAxis});
    registry.add("jets/matchDetJetPhiPartJetPhi", "matchDetJetPhiPartJetPhi; #phi^{jet, det}; #phi^{jet, part}", HistType::kTH2F, {detPhiAxis, partPhiAxis});
    registry.add("jets/matchPartJetPtResolutionPt", "#Delta = (#it{p}_{T}^{jet, part} - #it{p}_{T}^{jet, det}) / #it{p}_{T}^{jet, part}; #it{p}_{T}^{jet, part}; #Delta", HistType::kTH2F, {partJetPtAxis, diffAxis});
    registry.add("jets/matchPartJetPtResolutionEta", "(#eta^{jet, part} - #eta^{jet, det}) / #eta^{jet, part}; #eta^{jet, part}; #Delta", HistType::kTH3F, {partJetPtAxis, partEtaAxis, diffAxis});
    registry.add("jets/matchPartJetPtResolutionPhi", "(#phi^{jet, part} - #phi^{jet, det}) / #phi^{jet, part}; #phi^{jet, part}; #Delta", HistType::kTH3F, {partJetPtAxis, partPhiAxis, diffAxis});
    registry.add("jets/matchPartJetPtResolutionTrackProj", "Resolution #it{p}^{proj, part} / #it{p}^{jet, part}", HistType::kTH3F, {partJetPtAxis, partZAxis, diffAxis});
    registry.add("jets/matchPartJetPtResolutionChargeFrag", "Resolution #it{p}_{T}^{tr, part} / #it{p}_{T}^{jet, part}", HistType::kTH3F, {partJetPtAxis, partZAxis, diffAxis});
    registry.add("jets/matchPartJetPtMatchDist", "matchJetMatchDist; #it{p}_{T}^{part}; #Delta", HistType::kTH2F, {partJetPtAxis, matchDistAxis});
    registry.add("jets/matchPartJetPtEnergyScale", "jetEnergyScale; #it{p}_{T}^{part}; #it{p}_{T}^{part}/#it{p}_{T}^{det}", HistType::kTH2F, {partJetPtAxis, ratioAxis});

    // Response matrix, fakes, misses
    registry.add("jets/matchDetJetPtTrackProjPartJetPtTrackProj", "Matched; #it{p}_{T}^{jet, det}; #it{p}^{proj, det} / #it{p}^{jet, det}; #it{p}_{T}^{jet, part}; #it{p}^{proj, part} / #it{p}^{jet, part}", HistType::kTHnSparseF, {detJetPtAxis, detZAxis, partJetPtAxis, partZAxis});
    registry.add("jets/fakeDetJetPtTrackProj", "Fakes; #it{p}_{T}^{jet, det}; #it{p}^{proj, det} / #it{p}^{jet, det}", HistType::kTH2F, {detJetPtAxis, detZAxis});
    registry.add("jets/missPartJetPtTrackProj", "Misses; #it{p}_{T}^{jet, part}; #it{p}^{proj, part} / #it{p}^{jet, part}", HistType::kTH2F, {partJetPtAxis, partZAxis});

    registry.add("jets/matchDetJetPtFragPartJetPtFrag", "Matched; #it{p}_{T}^{jet, det}; #it{p}_{T}^{det} / #it{p}_{T}^{jet, det}; #it{p}_{T}^{jet, part}; #it{p}_{T}^{part} / #it{p}_{T}^{jet, part}", HistType::kTHnSparseF, {detJetPtAxis, detZAxis, partJetPtAxis, partZAxis});
    registry.add("jets/fakeDetJetPtFrag", "Fakes; #it{p}_{T}^{jet, det}; #it{p}_{T}^{det} / #it{p}_{T}^{jet, det}", HistType::kTH2F, {detJetPtAxis, detZAxis});
    registry.add("jets/missPartJetPtFrag", "Misses; #it{p}_{T}^{jet, part}; #it{p}_{T}^{part} / #it{p}_{T}^{jet, part}", HistType::kTH2F, {partJetPtAxis, partZAxis});

    // Bookkeeping
    registry.add("partCount", "partCount", HistType::kTH1I, {binCount});
    registry.add("detCount", "detCount", HistType::kTH1I, {binCount});
    registry.add("datCount", "datCount", HistType::kTH1I, {binCount});
    registry.add("partnJetnTrack", "partnJetnTrack; nJets; nTracks", HistType::kTH2I, {jetCount, trackCount});
    registry.add("detnJetnTrack", "detnJetnTrack; nJets; nTracks", HistType::kTH2I, {jetCount, trackCount});
    registry.add("nJetnTrack", "nJetnTrack; nJets; nTracks", HistType::kTH2I, {jetCount, trackCount});
  } // init

  double CheckDphi(double dphi)
  {
    if (dphi > TMath::Pi()) {
      return (dphi - 2 * TMath::Pi());
    } else if (dphi < -1 * TMath::Pi()) {
      return (dphi + 2 * TMath::Pi());
    }
    return dphi;
  }

  void processDummy(aod::Tracks const& track) {}
  PROCESS_SWITCH(JetFragmentation, processDummy, "Dummy process function turned on by default", true);

  void processMcD(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                  McDJets const& jets,
                  aod::Tracks const& tracks)
  {
    int nJets = -1, nTracks = -1;
    nJets = jets.size(), nTracks = tracks.size();
    registry.fill(HIST("detnJetnTrack"), nJets, nTracks);
    for (const auto& track : tracks) {
      if (track.pt() > 0.1) {
        registry.fill(HIST("tracks/detTrackPtEtaPhi"), track.pt(), track.eta(), track.phi());
      }
    }
    for (const auto& jet : jets) {
      registry.fill(HIST("jets/detJetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi());
      for (const auto& track : jet.tracks_as<aod::Tracks>()) {
        double chargeFrag = -1., trackProj = -1.;
        trackProj = track.px() * jet.px() + track.py() * jet.py() + track.pz() * jet.pz();
        trackProj /= (jet.p() * jet.p());
        chargeFrag = track.pt() / jet.pt();

        registry.fill(HIST("jets/detJetPtTrackPt"), jet.pt(), track.pt());
        registry.fill(HIST("jets/detJetTrackPtEtaPhi"), track.pt(), track.eta(), track.phi());
        registry.fill(HIST("jets/detJetPtFrag"), jet.pt(), chargeFrag);
        registry.fill(HIST("jets/detJetPtTrackProj"), jet.pt(), trackProj);
      }
    }
  }
  PROCESS_SWITCH(JetFragmentation, processMcD, "Monte Carlo detector level", false);

  void processMcP(aod::McCollision const& mcCollision, // Add some form of event selection?
                  McPJets const& jets,
                  aod::McParticles const& particles)
  {
    int nJets = -1, nTracks = -1;
    nJets = jets.size(), nTracks = particles.size();
    registry.fill(HIST("partnJetnTrack"), nJets, nTracks);
    for (const auto& particle : particles) {
      if (particle.pt() > 0.1) {
        registry.fill(HIST("tracks/partTrackPtEtaPhi"), particle.pt(), particle.eta(), particle.phi());
      }
    }
    for (const auto& jet : jets) {
      registry.fill(HIST("jets/partJetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi());
      for (const auto& track : jet.tracks_as<aod::McParticles>()) {
        double chargeFrag = -1., trackProj = -1.;
        trackProj = track.px() * jet.px() + track.py() * jet.py() + track.pz() * jet.pz();
        trackProj /= (jet.p() * jet.p());
        chargeFrag = track.pt() / jet.pt();

        registry.fill(HIST("jets/partJetPtTrackPt"), jet.pt(), track.pt());
        registry.fill(HIST("jets/partJetTrackPtEtaPhi"), track.pt(), track.eta(), track.phi());
        registry.fill(HIST("jets/partJetPtFrag"), jet.pt(), chargeFrag);
        registry.fill(HIST("jets/partJetPtTrackProj"), jet.pt(), trackProj);
      }
    }
  }
  PROCESS_SWITCH(JetFragmentation, processMcP, "Monte Carlo particle level", false);

  void processDataRun3(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                       soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                       aod::Tracks const& tracks)
  {
    registry.fill(HIST("datCount"), 1);
    int nJets = -1, nTracks = -1;
    nJets = jets.size(), nTracks = tracks.size();
    registry.fill(HIST("nJetnTrack"), nJets, nTracks);
    for (const auto& track : tracks) {
      if (track.pt() > 0.1) {
        registry.fill(HIST("tracks/trackPtEtaPhi"), track.pt(), track.eta(), track.phi());
      }
    }
    for (const auto& jet : jets) {
      registry.fill(HIST("jets/jetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi());
      for (const auto& track : jet.tracks_as<aod::Tracks>()) {
        double chargeFrag = -1., trackProj = -1.;
        trackProj = track.px() * jet.px() + track.py() * jet.py() + track.pz() * jet.pz();
        trackProj /= (jet.p() * jet.p());
        chargeFrag = track.pt() / jet.pt();

        registry.fill(HIST("jets/jetPtTrackPt"), jet.pt(), track.pt());
        registry.fill(HIST("jets/jetTrackPtEtaPhi"), track.pt(), track.eta(), track.phi());
        registry.fill(HIST("jets/jetPtFrag"), jet.pt(), chargeFrag);
        registry.fill(HIST("jets/jetPtTrackProj"), jet.pt(), trackProj);
      }
    }
  }
  PROCESS_SWITCH(JetFragmentation, processDataRun3, "Run 3 Data", false);

  // Taken from jet-validation
  void processMcDP(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision,
                   //  McDJets const& mcDetJets,
                   soa::Join<McDJets, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets> const& mcDetJets,
                   McTracks const& tracks,
                   aod::McCollisions const& mcCollisions,
                   McPJets const& mcPartJets,
                   aod::McParticles const& mcParticles)
  {
    for (const auto& detJet : mcDetJets) {
      if (detJet.has_matchedJetGeo()) {
        const auto& partJet = detJet.template matchedJetGeo_as<McPJets>();
        double deltaEta = partJet.eta() - detJet.eta();
        double deltaPhi = partJet.phi() - detJet.phi();
        deltaPhi = CheckDphi(deltaPhi);

        registry.fill(HIST("jets/matchDetJetPtEtaPhi"), detJet.pt(), detJet.eta(), detJet.phi());
        registry.fill(HIST("jets/matchPartJetPtEtaPhi"), partJet.pt(), partJet.eta(), partJet.phi());

        registry.fill(HIST("jets/matchPartJetPtMatchDist"), partJet.pt(), TMath::Sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi));
        registry.fill(HIST("jets/matchPartJetPtEnergyScale"), partJet.pt(), detJet.pt() / partJet.pt());
        registry.fill(HIST("jets/matchDetJetPtPartJetPt"), detJet.pt(), partJet.pt());
        registry.fill(HIST("jets/matchDetJetEtaPartJetEta"), detJet.eta(), partJet.eta());
        registry.fill(HIST("jets/matchDetJetPhiPartJetPhi"), detJet.phi(), partJet.phi());
        registry.fill(HIST("jets/matchPartJetPtResolutionPt"), partJet.pt(), (partJet.pt() - detJet.pt()) / partJet.pt());
        registry.fill(HIST("jets/matchPartJetPtResolutionEta"), partJet.pt(), partJet.eta(), (partJet.eta() - detJet.eta()) / partJet.eta());
        registry.fill(HIST("jets/matchPartJetPtResolutionPhi"), partJet.pt(), partJet.phi(), CheckDphi(partJet.phi() - detJet.phi()) / partJet.phi());

        for (const auto& track : detJet.tracks_as<McTracks>()) {
          bool isTrackMatched = false;
          double detChargeFrag = -1., detTrackProj = -1.;
          detTrackProj = track.px() * detJet.px() + track.py() * detJet.py() + track.pz() * detJet.pz();
          detTrackProj /= (detJet.p() * detJet.p());
          detChargeFrag = track.pt() / detJet.pt();

          registry.fill(HIST("jets/matchDetJetTrackPtEtaPhi"), track.pt(), track.eta(), track.phi());
          registry.fill(HIST("jets/matchDetJetPtTrackPt"), detJet.pt(), track.pt());
          registry.fill(HIST("jets/matchDetJetPtTrackProj"), detJet.pt(), detTrackProj);
          registry.fill(HIST("jets/matchDetJetPtFrag"), detJet.pt(), detChargeFrag);

          for (const auto& particle : partJet.tracks_as<aod::McParticles>()) {
            if (track.has_mcParticle() && particle.globalIndex() == track.template mcParticle_as<aod::McParticles>().globalIndex()) {
              isTrackMatched = true;
              double partChargeFrag = -1., partTrackProj = -1.;
              partTrackProj = particle.px() * partJet.px() + particle.py() * partJet.py() + particle.pz() * partJet.pz();
              partTrackProj /= (partJet.p() * partJet.p());
              partChargeFrag = particle.pt() / partJet.pt();

              registry.fill(HIST("jets/matchPartJetTrackPtEtaPhi"), particle.pt(), particle.eta(), particle.phi());
              registry.fill(HIST("jets/matchPartJetPtTrackPt"), partJet.pt(), particle.pt());
              registry.fill(HIST("jets/matchPartJetPtTrackProj"), partJet.pt(), partTrackProj);
              registry.fill(HIST("jets/matchPartJetPtFrag"), partJet.pt(), partChargeFrag);

              registry.fill(HIST("jets/matchPartJetPtResolutionTrackProj"), partJet.pt(), partTrackProj, (partTrackProj - detTrackProj) / partTrackProj);
              registry.fill(HIST("jets/matchPartJetPtResolutionChargeFrag"), partJet.pt(), partChargeFrag, (partChargeFrag - detChargeFrag) / partChargeFrag);
              // Response
              registry.fill(HIST("jets/matchDetJetPtTrackProjPartJetPtTrackProj"), detJet.pt(), detTrackProj, partJet.pt(), partTrackProj);
              registry.fill(HIST("jets/matchDetJetPtFragPartJetPtFrag"), detJet.pt(), detChargeFrag, partJet.pt(), partChargeFrag);
              break; // No need to inspect other particles
            }        // if track has mcParticle and particle is in matched jet
          }          // for particle in matched partJet
          if (!isTrackMatched) {
            registry.fill(HIST("jets/fakeDetJetPtTrackProj"), detJet.pt(), detTrackProj);
            registry.fill(HIST("jets/fakeDetJetPtFrag"), detJet.pt(), detChargeFrag);
          } // if track is not matched
        }   // for detJet tracks
      } else if (!detJet.has_matchedJetGeo()) {
        for (const auto& track : detJet.tracks_as<McTracks>()) {
          double detChargeFrag = -1., detTrackProj = -1.;
          detTrackProj = track.px() * detJet.px() + track.py() * detJet.py() + track.pz() * detJet.pz();
          detTrackProj /= (detJet.p() * detJet.p());
          detChargeFrag = track.pt() / detJet.pt();
          registry.fill(HIST("jets/fakeDetJetPtTrackProj"), detJet.pt(), detTrackProj);
          registry.fill(HIST("jets/fakeDetJetPtFrag"), detJet.pt(), detChargeFrag);
        }
      } // if detJet does not have a match
    }   // for det jet
    // TODO: how to deal with misses?
    //       -> ParticleToDetector table is currently bugged: size does not correspond to ParticleLevelJets (04.06.2023)
  }
  PROCESS_SWITCH(JetFragmentation, processMcDP, "Monte Carlo particle and detector level", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetFragmentation>(cfgc, TaskName{"jet-fragmentation"})};
}
