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

using McTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>;
using McDJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>;
// using MatchedMcDJets = soa::Filtered<soa::Join<McDJets, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>>;
using McPJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>;

struct JetFragmentation {
  HistogramRegistry registry{"registry"};

  std::vector<int> pdgVector = {211, 321, 2212, 111, 130, 310, 311, 3122};
  std::vector<std::string> hadronVector = {"#pi^{#pm}", "#it{K}^{#pm}", "#it{p}^{#pm}", "#pi^{0}", "#it{K}^{0}_{L}", "#it{K}^{0}_{S}", "#it{K}^{0}", "#Lambda^{0}"};

  Configurable<float> matchedDetJetEtaMin{"matchedDetJetEtaMin", -0.5, "minimum matchedDetJet eta"};
  Configurable<float> matchedDetJetEtaMax{"matchedDetJetEtaMax", 0.5, "maximum matchedDetJet eta"};

  // Binning
  ConfigurableAxis binJetPt{"binJetPt", {40, 0.f, 200.f}, ""};
  ConfigurableAxis binEta{"binEta", {20, -1.f, 1.f}, ""};
  ConfigurableAxis binPhi{"binPhi", {18 * 8, 0.f, 2. * TMath::Pi()}, ""};
  ConfigurableAxis binZ{"binZ", {20, -5e-3f, 1.f + 5e-3f}, ""};
  ConfigurableAxis binJetR{"binJetR", {6, 0.05f, 0.65f}, ""};
  ConfigurableAxis binTrackPt{"binTrackPt", {200, 0.f, 100.f}, ""};
  ConfigurableAxis binPDG{"binPDG", {static_cast<double>(pdgVector.size()), -0.5f, static_cast<double>(pdgVector.size()) - 0.5f}, ""};
  ConfigurableAxis binVtxZ{"binVtxZ", {200, -20, 20}, ""};

  ConfigurableAxis binDiff{"binDiff", {51, -5.5f, 5.5f}, ""};
  ConfigurableAxis binRatio{"binRatio", {100, -0.5f, 9.5f}, ""};      // Ratio of pt, eta, phi
  ConfigurableAxis binMatchDist{"binMatchDist", {10, 0.f, 0.5f}, ""}; // Distance between matched jets

  ConfigurableAxis binCount{"binCount", {1, .5f, 1.5f}, ""};
  ConfigurableAxis jetCount{"jetCount", {50, -.5f, 49.5f}, ""};
  ConfigurableAxis trackCount{"trackCount", {100, -.5f, 99.5f}, ""};

  Preslice<McPJets> perMcPJet = aod::jet::mcCollisionId;
  // Filter matchedDetJetFilter = (aod::chargedmcdetectorleveljets::eta >= matchedDetJetEtaMin && aod::chargedmcdetectorleveljets::eta <= matchedDetJetEtaMax);

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
    registry.add("data/collision/collisionVtxZ", "Collision vertex z (cm)", HistType::kTH1F, {binVtxZ});

    registry.add("data/tracks/trackPtEtaPhi", "trackPtEtaPhi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});

    registry.add("data/jets/jetPtEtaPhi", "Jet #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {jetPtAxis, etaAxis, phiAxis});
    registry.add("data/jets/jetPtTrackPt", "Jet #it{p}_{T}, track #it{p}_{T};#it{p}_{T}^{jet};#it{p}_{T}^{track}", HistType::kTH2F, {jetPtAxis, trackPtAxis});
    registry.add("data/jets/jetTrackPtEtaPhi", "Tracks in jets #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});
    registry.add("data/jets/jetPtTrackProj", "Jet #it{p}_{T}, #it{z}; #it{p}_{T}; #it{z}", HistType::kTH2F, {jetPtAxis, zAxis});
    registry.add("data/jets/jetPtFrag", "Jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}; #it{p}_{T}; #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2F, {jetPtAxis, zAxis});

    // MC particle level
    registry.add("particle-level/collision/partCollisionVtxZ", "Collision vertex z (cm)", HistType::kTH1F, {binVtxZ});

    registry.add("particle-level/tracks/partTrackPtEtaPhi", "partTrackPtEtaPhi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});

    registry.add("particle-level/jets/partJetPtEtaPhi", "Particle level jet #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {partJetPtAxis, partEtaAxis, partPhiAxis});
    registry.add("particle-level/jets/partJetPtTrackPt", "Particle level jet #it{p}_{T}, track #it{p}_{T};#it{p}_{T}^{jet};#it{p}_{T}^{track}", HistType::kTH2F, {partJetPtAxis, trackPtAxis});
    registry.add("particle-level/jets/partJetTrackPtEtaPhi", "Particle level tracks in jets #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {trackPtAxis, partEtaAxis, partPhiAxis});
    registry.add("particle-level/jets/partJetPtTrackProj", "Particle level jet #it{p}_{T}, #it{z}; #it{p}_{T}; #it{z}", HistType::kTH2F, {partJetPtAxis, partZAxis});
    registry.add("particle-level/jets/partJetPtFrag", "Particle level jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}; #it{p}_{T}; #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2F, {partJetPtAxis, partZAxis});

    // MC detector level
    registry.add("detector-level/collision/detCollisionVtxZ", "Collision vertex z (cm)", HistType::kTH1F, {binVtxZ});

    registry.add("detector-level/tracks/detTrackPtEtaPhi", "detTrackPtEtaPhi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});

    registry.add("detector-level/jets/detJetPtEtaPhi", "Detector level jet #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {detJetPtAxis, detEtaAxis, detPhiAxis});
    registry.add("detector-level/jets/detJetPtTrackPt", "Detector level jet #it{p}_{T}, track #it{p}_{T};#it{p}_{T}^{jet};#it{p}_{T}^{track}", HistType::kTH2F, {detJetPtAxis, trackPtAxis});
    registry.add("detector-level/jets/detJetTrackPtEtaPhi", "Detector level tracks in jets #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {trackPtAxis, detEtaAxis, detPhiAxis});
    registry.add("detector-level/jets/detJetPtTrackProj", "Detector level jet #it{p}_{T}, #it{z}; #it{p}_{T}; #it{z}", HistType::kTH2F, {detJetPtAxis, detZAxis});
    registry.add("detector-level/jets/detJetPtFrag", "Detector level jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}; #it{p}_{T}; #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2F, {detJetPtAxis, detZAxis});

    // MC particle-detector level matching
    registry.add("matching/collision/matchCollisionVtxZ", "Collision vertex z (cm)", HistType::kTH1F, {binVtxZ});

    registry.add("matching/tracks/matchDetTrackPtEtaPhi", "matchDetTrackPtEtaPhi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});
    registry.add("matching/tracks/matchPartTrackPtEtaPhi", "matchPartTrackPtEtaPhi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});
    registry.add("matching/tracks/matchDetTrackPtPartTrackPt", "matchDetTrackPtPartTrackPt", HistType::kTH2F, {trackPtAxis, trackPtAxis});
    registry.add("matching/tracks/matchDetTrackEtaPartTrackEta", "matchDetTrackEtaPartTrackEta", HistType::kTH2F, {etaAxis, etaAxis});
    registry.add("matching/tracks/matchDetTrackPhiPartTrackPhi", "matchDetTrackPhiPartTrackPhi", HistType::kTH2F, {phiAxis, phiAxis});
    registry.add("matching/tracks/trackResolutionPt", "trackResolutionPt; #Delta #it{p}_{T}^{tr}", HistType::kTH2F, {trackPtAxis, diffAxis});
    registry.add("matching/tracks/trackResolutionEta", "trackResolutionEta; #Delta #eta", HistType::kTH2F, {etaAxis, diffAxis});
    registry.add("matching/tracks/trackResolutionPhi", "trackResolutionPhi; #Delta #phi", HistType::kTH2F, {phiAxis, diffAxis});

    // Detector level jets with a match
    registry.add("matching/jets/matchDetJetPtEtaPhi", "Matched detector level jet #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {detJetPtAxis, detEtaAxis, detPhiAxis});
    registry.add("matching/jets/matchDetJetPtTrackPt", "Matched detector level jet #it{p}_{T}, track #it{p}_{T};#it{p}_{T}^{jet};#it{p}_{T}^{track}", HistType::kTH2F, {detJetPtAxis, trackPtAxis});
    registry.add("matching/jets/matchDetJetTrackPtEtaPhi", "Matched detector level tracks in jets #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {trackPtAxis, detEtaAxis, detPhiAxis});
    registry.add("matching/jets/matchDetJetPtTrackProj", "Matched detector level jet #it{p}_{T}, #it{z}; #it{p}_{T}; #it{z}", HistType::kTH2F, {detJetPtAxis, detZAxis});
    registry.add("matching/jets/matchDetJetPtFrag", "Matched detector level jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}; #it{p}_{T}; #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2F, {detJetPtAxis, detZAxis});
    // Particle level jets with a match
    registry.add("matching/jets/matchPartJetPtEtaPhi", "Matched particle level jet #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {partJetPtAxis, partEtaAxis, partPhiAxis});
    registry.add("matching/jets/matchPartJetPtTrackPt", "Matched particle level jet #it{p}_{T}, track #it{p}_{T};#it{p}_{T}^{jet};#it{p}_{T}^{track}", HistType::kTH2F, {partJetPtAxis, trackPtAxis});
    registry.add("matching/jets/matchPartJetTrackPtEtaPhi", "Matched particle level tracks in jets #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {trackPtAxis, partEtaAxis, partPhiAxis});
    registry.add("matching/jets/matchPartJetPtTrackProj", "Matched particle level jet #it{p}_{T}, #it{z}; #it{p}_{T}; #it{z}", HistType::kTH2F, {partJetPtAxis, partZAxis});
    registry.add("matching/jets/matchPartJetPtFrag", "Matched particle level jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}; #it{p}_{T}; #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2F, {partJetPtAxis, partZAxis});
    // Combined information of matched jets
    registry.add("matching/jets/matchDetJetPtPartJetPt", "matchDetJetPtPartJetPt; #it{p}_{T}^{jet, det}; #it{p}_{T}^{jet, part}", HistType::kTH2F, {detJetPtAxis, partJetPtAxis});
    registry.add("matching/jets/matchPartJetPtDetJetEtaPartJetEta", "matchPartJetPtDetJetEtaPartJetEta; #it{p}_{T}^{jet, part}; #eta^{jet, det}; #eta^{jet, part}", HistType::kTH3F, {partJetPtAxis, detEtaAxis, partEtaAxis});
    registry.add("matching/jets/matchPartJetPtDetJetPhiPartJetPhi", "matchPartJetPtDetJetPhiPartJetPhi; #it{p}_{T}^{jet, part}; #phi^{jet, det}; #phi^{jet, part}", HistType::kTH3F, {partJetPtAxis, detPhiAxis, partPhiAxis});
    registry.add("matching/jets/matchPartJetPtResolutionPt", "#Delta = (#it{p}_{T}^{jet, part} - #it{p}_{T}^{jet, det}) / #it{p}_{T}^{jet, part}; #it{p}_{T}^{jet, part}; #Delta", HistType::kTH2F, {partJetPtAxis, diffAxis});
    registry.add("matching/jets/matchPartJetPtResolutionEta", "(#eta^{jet, part} - #eta^{jet, det}); #eta^{jet, part}; #Delta", HistType::kTH3F, {partJetPtAxis, partEtaAxis, diffAxis});
    registry.add("matching/jets/matchPartJetPtResolutionPhi", "(#phi^{jet, part} - #phi^{jet, det}); #phi^{jet, part}; #Delta", HistType::kTH3F, {partJetPtAxis, partPhiAxis, diffAxis});
    registry.add("matching/jets/matchPartJetPtResolutionTrackProj", "Resolution #it{p}^{proj, part} / #it{p}^{jet, part}", HistType::kTH3F, {partJetPtAxis, partZAxis, diffAxis});
    registry.add("matching/jets/matchPartJetPtResolutionChargeFrag", "Resolution #it{p}_{T}^{tr, part} / #it{p}_{T}^{jet, part}", HistType::kTH3F, {partJetPtAxis, partZAxis, diffAxis});
    registry.add("matching/jets/matchPartJetPtMatchDist", "matchJetMatchDist; #it{p}_{T}^{part}; #Delta", HistType::kTH2F, {partJetPtAxis, matchDistAxis});
    registry.add("matching/jets/matchPartJetPtEnergyScale", "jetEnergyScale; #it{p}_{T}^{part}; #it{p}_{T}^{part}/#it{p}_{T}^{det}", HistType::kTH2F, {partJetPtAxis, ratioAxis});

    // Response matrix, fakes, misses
    registry.add("matching/jets/matchDetJetPtTrackProjPartJetPtTrackProj", "Matched; #it{p}_{T}^{jet, det}; #it{p}^{proj, det} / #it{p}^{jet, det}; #it{p}_{T}^{jet, part}; #it{p}^{proj, part} / #it{p}^{jet, part}", HistType::kTHnSparseF, {detJetPtAxis, detZAxis, partJetPtAxis, partZAxis});
    registry.add("matching/jets/fakeDetJetPtTrackProj", "Fakes; #it{p}_{T}^{jet, det}; #it{p}^{proj, det} / #it{p}^{jet, det}", HistType::kTH2F, {detJetPtAxis, detZAxis});
    registry.add("matching/jets/missPartJetPtTrackProj", "Misses; #it{p}_{T}^{jet, part}; #it{p}^{proj, part} / #it{p}^{jet, part}", HistType::kTH2F, {partJetPtAxis, partZAxis});

    registry.add("matching/jets/matchDetJetPtFragPartJetPtFrag", "Matched; #it{p}_{T}^{jet, det}; #it{p}_{T}^{det} / #it{p}_{T}^{jet, det}; #it{p}_{T}^{jet, part}; #it{p}_{T}^{part} / #it{p}_{T}^{jet, part}", HistType::kTHnSparseF, {detJetPtAxis, detZAxis, partJetPtAxis, partZAxis});
    registry.add("matching/jets/fakeDetJetPtFrag", "Fakes; #it{p}_{T}^{jet, det}; #it{p}_{T}^{det} / #it{p}_{T}^{jet, det}", HistType::kTH2F, {detJetPtAxis, detZAxis});
    registry.add("matching/jets/missPartJetPtFrag", "Misses; #it{p}_{T}^{jet, part}; #it{p}_{T}^{part} / #it{p}_{T}^{jet, part}", HistType::kTH2F, {partJetPtAxis, partZAxis});

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
        registry.fill(HIST("detector-level/tracks/detTrackPtEtaPhi"), track.pt(), track.eta(), track.phi());
      }
    }
    for (const auto& jet : jets) {
      registry.fill(HIST("detector-level/jets/detJetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi());
      for (const auto& track : jet.tracks_as<aod::Tracks>()) {
        double chargeFrag = -1., trackProj = -1.;
        trackProj = track.px() * jet.px() + track.py() * jet.py() + track.pz() * jet.pz();
        trackProj /= (jet.p() * jet.p());
        chargeFrag = track.pt() / jet.pt();

        registry.fill(HIST("detector-level/jets/detJetPtTrackPt"), jet.pt(), track.pt());
        registry.fill(HIST("detector-level/jets/detJetTrackPtEtaPhi"), track.pt(), track.eta(), track.phi());
        registry.fill(HIST("detector-level/jets/detJetPtFrag"), jet.pt(), chargeFrag);
        registry.fill(HIST("detector-level/jets/detJetPtTrackProj"), jet.pt(), trackProj);
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
        registry.fill(HIST("particle-level/tracks/partTrackPtEtaPhi"), particle.pt(), particle.eta(), particle.phi());
      }
    }
    for (const auto& jet : jets) {
      registry.fill(HIST("particle-level/jets/partJetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi());
      for (const auto& track : jet.tracks_as<aod::McParticles>()) {
        double chargeFrag = -1., trackProj = -1.;
        trackProj = track.px() * jet.px() + track.py() * jet.py() + track.pz() * jet.pz();
        trackProj /= (jet.p() * jet.p());
        chargeFrag = track.pt() / jet.pt();

        registry.fill(HIST("particle-level/jets/partJetPtTrackPt"), jet.pt(), track.pt());
        registry.fill(HIST("particle-level/jets/partJetTrackPtEtaPhi"), track.pt(), track.eta(), track.phi());
        registry.fill(HIST("particle-level/jets/partJetPtFrag"), jet.pt(), chargeFrag);
        registry.fill(HIST("particle-level/jets/partJetPtTrackProj"), jet.pt(), trackProj);
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
        registry.fill(HIST("data/tracks/trackPtEtaPhi"), track.pt(), track.eta(), track.phi());
      }
    }
    for (const auto& jet : jets) {
      registry.fill(HIST("data/jets/jetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi());
      for (const auto& track : jet.tracks_as<aod::Tracks>()) {
        double chargeFrag = -1., trackProj = -1.;
        trackProj = track.px() * jet.px() + track.py() * jet.py() + track.pz() * jet.pz();
        trackProj /= (jet.p() * jet.p());
        chargeFrag = track.pt() / jet.pt();

        registry.fill(HIST("data/jets/jetPtTrackPt"), jet.pt(), track.pt());
        registry.fill(HIST("data/jets/jetTrackPtEtaPhi"), track.pt(), track.eta(), track.phi());
        registry.fill(HIST("data/jets/jetPtFrag"), jet.pt(), chargeFrag);
        registry.fill(HIST("data/jets/jetPtTrackProj"), jet.pt(), trackProj);
      }
    }
  }
  PROCESS_SWITCH(JetFragmentation, processDataRun3, "Run 3 Data", false);

  // Taken from jet-validation
  void processMcDP(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision,
                   //  MatchedMcDJets const& mcDetJets,
                   soa::Join<McDJets, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets> const& mcDetJets,
                   McTracks const& tracks,
                   aod::McCollisions const& mcCollisions,
                   McPJets const& mcPartJets,
                   aod::McParticles const& mcParticles)
  {
    for (const auto& detJet : mcDetJets) {
      if (detJet.eta() < matchedDetJetEtaMin || detJet.eta() > matchedDetJetEtaMax) {
        continue; // TODO: should be done in filter
      }
      for (auto& partJet : detJet.template matchedJetGeo_as<McPJets>()) {
        double deltaEta = partJet.eta() - detJet.eta();
        double deltaPhi = partJet.phi() - detJet.phi();
        deltaPhi = CheckDphi(deltaPhi);

        registry.fill(HIST("matching/jets/matchDetJetPtEtaPhi"), detJet.pt(), detJet.eta(), detJet.phi());
        registry.fill(HIST("matching/jets/matchPartJetPtEtaPhi"), partJet.pt(), partJet.eta(), partJet.phi());

        registry.fill(HIST("matching/jets/matchPartJetPtMatchDist"), partJet.pt(), TMath::Sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi));
        registry.fill(HIST("matching/jets/matchPartJetPtEnergyScale"), partJet.pt(), detJet.pt() / partJet.pt());
        registry.fill(HIST("matching/jets/matchDetJetPtPartJetPt"), detJet.pt(), partJet.pt());
        registry.fill(HIST("matching/jets/matchPartJetPtDetJetEtaPartJetEta"), partJet.pt(), detJet.eta(), partJet.eta());
        registry.fill(HIST("matching/jets/matchPartJetPtDetJetPhiPartJetPhi"), partJet.pt(), detJet.phi(), partJet.phi());
        registry.fill(HIST("matching/jets/matchPartJetPtResolutionPt"), partJet.pt(), (partJet.pt() - detJet.pt()));
        registry.fill(HIST("matching/jets/matchPartJetPtResolutionEta"), partJet.pt(), partJet.eta(), (partJet.eta() - detJet.eta()));
        registry.fill(HIST("matching/jets/matchPartJetPtResolutionPhi"), partJet.pt(), partJet.phi(), CheckDphi(partJet.phi() - detJet.phi()));

        for (const auto& track : detJet.tracks_as<McTracks>()) {
          bool isTrackMatched = false;
          double detChargeFrag = -1., detTrackProj = -1.;
          detTrackProj = track.px() * detJet.px() + track.py() * detJet.py() + track.pz() * detJet.pz();
          detTrackProj /= (detJet.p() * detJet.p());
          detChargeFrag = track.pt() / detJet.pt();

          registry.fill(HIST("matching/jets/matchDetJetTrackPtEtaPhi"), track.pt(), track.eta(), track.phi());
          registry.fill(HIST("matching/jets/matchDetJetPtTrackPt"), detJet.pt(), track.pt());
          registry.fill(HIST("matching/jets/matchDetJetPtTrackProj"), detJet.pt(), detTrackProj);
          registry.fill(HIST("matching/jets/matchDetJetPtFrag"), detJet.pt(), detChargeFrag);

          for (const auto& particle : partJet.tracks_as<aod::McParticles>()) {
            if (track.has_mcParticle() && particle.globalIndex() == track.template mcParticle_as<aod::McParticles>().globalIndex()) {
              isTrackMatched = true;
              double partChargeFrag = -1., partTrackProj = -1.;
              partTrackProj = particle.px() * partJet.px() + particle.py() * partJet.py() + particle.pz() * partJet.pz();
              partTrackProj /= (partJet.p() * partJet.p());
              partChargeFrag = particle.pt() / partJet.pt();

              registry.fill(HIST("matching/jets/matchPartJetTrackPtEtaPhi"), particle.pt(), particle.eta(), particle.phi());
              registry.fill(HIST("matching/jets/matchPartJetPtTrackPt"), partJet.pt(), particle.pt());
              registry.fill(HIST("matching/jets/matchPartJetPtTrackProj"), partJet.pt(), partTrackProj);
              registry.fill(HIST("matching/jets/matchPartJetPtFrag"), partJet.pt(), partChargeFrag);

              registry.fill(HIST("matching/jets/matchPartJetPtResolutionTrackProj"), partJet.pt(), partTrackProj, (partTrackProj - detTrackProj));
              registry.fill(HIST("matching/jets/matchPartJetPtResolutionChargeFrag"), partJet.pt(), partChargeFrag, (partChargeFrag - detChargeFrag));
              // Response
              registry.fill(HIST("matching/jets/matchDetJetPtTrackProjPartJetPtTrackProj"), detJet.pt(), detTrackProj, partJet.pt(), partTrackProj);
              registry.fill(HIST("matching/jets/matchDetJetPtFragPartJetPtFrag"), detJet.pt(), detChargeFrag, partJet.pt(), partChargeFrag);
              break; // No need to inspect other particles
            }        // if track has mcParticle and particle is in matched jet
          }          // for particle in matched partJet
          if (!isTrackMatched) {
            registry.fill(HIST("matching/jets/fakeDetJetPtTrackProj"), detJet.pt(), detTrackProj);
            registry.fill(HIST("matching/jets/fakeDetJetPtFrag"), detJet.pt(), detChargeFrag);
          } // if track is not matched
        }   // for detJet tracks
      }
      if (!detJet.has_matchedJetGeo()) {
        for (const auto& track : detJet.tracks_as<McTracks>()) {
          double detChargeFrag = -1., detTrackProj = -1.;
          detTrackProj = track.px() * detJet.px() + track.py() * detJet.py() + track.pz() * detJet.pz();
          detTrackProj /= (detJet.p() * detJet.p());
          detChargeFrag = track.pt() / detJet.pt();
          registry.fill(HIST("matching/jets/fakeDetJetPtTrackProj"), detJet.pt(), detTrackProj);
          registry.fill(HIST("matching/jets/fakeDetJetPtFrag"), detJet.pt(), detChargeFrag);
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
