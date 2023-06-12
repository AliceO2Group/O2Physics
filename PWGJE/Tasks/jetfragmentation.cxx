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

struct JetFragmentation {
  HistogramRegistry registry{"registry"};

  std::vector<int> pdgVector = {211, 321, 2212, 111, 130, 310, 311, 3122};
  std::vector<std::string> hadronVector = {"#pi^{#pm}", "#it{K}^{#pm}", "#it{p}^{#pm}", "#pi^{0}", "#it{K}^{0}_{L}", "#it{K}^{0}_{S}", "#it{K}^{0}", "#Lambda^{0}"};

  // Binning
  ConfigurableAxis binJetPt{"binJetPt", {200, 0.f, 200.f}, ""};
  ConfigurableAxis binEta{"binEta", {100, -1.f, 1.f}, ""};
  ConfigurableAxis binPhi{"binPhi", {18 * 8, 0.f, 2. * TMath::Pi()}, ""};
  ConfigurableAxis binZ{"binZ", {100, -5e-3f, 1.f + 5e-3f}, ""};
  ConfigurableAxis binPDG{"binPDG", {static_cast<double>(pdgVector.size()), -0.5f, static_cast<double>(pdgVector.size()) - 0.5f}, ""};
  ConfigurableAxis binJetR{"binJetR", {6, 0.05f, 0.65f}, ""};
  ConfigurableAxis binTrackPt{"binTrackPt", {200, 0.f, 100.f}, ""};

  ConfigurableAxis binCount{"binCount", {1, .5f, 1.5f}, ""};
  ConfigurableAxis jetCount{"jetCount", {50, -.5f, 49.5f}, ""};
  ConfigurableAxis trackCount{"trackCount", {100, -.5f, 99.5f}, ""};

  void init(InitContext& initContext)
  {
    // Axis
    AxisSpec jetPtAxis = {binJetPt, "#it{p}_{T}^{jet}"};
    AxisSpec etaAxis = {binEta, "#eta"};
    AxisSpec phiAxis = {binPhi, "#phi"};
    AxisSpec zAxis = {binZ, "#it{z}"};
    AxisSpec pdgAxis = {binPDG, ""};
    AxisSpec rAxis = {binJetR, "#it{R}"};
    AxisSpec trackPtAxis = {binTrackPt, "#it{p}_{T}^{tr}"};

    // Data
    registry.add("tracks/trackPtEtaPhi", "trackPtEtaPhi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});

    registry.add("jets/jetPtEtaPhi", "Jet #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {jetPtAxis, etaAxis, phiAxis});
    registry.add("jets/jetPtTrackPt", "Jet #it{p}_{T}, track #it{p}_{T};#it{p}_{T}^{jet};#it{p}_{T}^{track}", HistType::kTH2F, {jetPtAxis, trackPtAxis});
    registry.add("jets/jetTrackPtEtaPhi", "Tracks in jets #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});
    registry.add("jets/jetPtTrackProj", "Jet #it{p}_{T}, #it{z}; #it{p}_{T}; #it{z}", HistType::kTH2F, {jetPtAxis, zAxis});
    registry.add("jets/jetPtFrag", "Jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}; #it{p}_{T}; #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2F, {jetPtAxis, zAxis});

    // MC particle level
    registry.add("tracks/partTrackPtEtaPhi", "partTrackPtEtaPhi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});

    registry.add("jets/partJetPtEtaPhi", "Particle level jet #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {jetPtAxis, etaAxis, phiAxis});
    registry.add("jets/partJetPtTrackPt", "Particle level jet #it{p}_{T}, track #it{p}_{T};#it{p}_{T}^{jet};#it{p}_{T}^{track}", HistType::kTH2F, {jetPtAxis, trackPtAxis});
    registry.add("jets/partJetTrackPtEtaPhi", "Particle level tracks in jets #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});
    registry.add("jets/partJetPtTrackProj", "Particle level jet #it{p}_{T}, #it{z}; #it{p}_{T}; #it{z}", HistType::kTH2F, {jetPtAxis, zAxis});
    registry.add("jets/partJetPtFrag", "Particle level jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}; #it{p}_{T}; #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2F, {jetPtAxis, zAxis});

    // MC detector level
    registry.add("tracks/detTrackPtEtaPhi", "detTrackPtEtaPhi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});

    registry.add("jets/detJetPtEtaPhi", "Detector level jet #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {jetPtAxis, etaAxis, phiAxis});
    registry.add("jets/detJetPtTrackPt", "Detector level jet #it{p}_{T}, track #it{p}_{T};#it{p}_{T}^{jet};#it{p}_{T}^{track}", HistType::kTH2F, {jetPtAxis, trackPtAxis});
    registry.add("jets/detJetTrackPtEtaPhi", "Detector level tracks in jets #it{p}_{T}, #eta, #phi;#it{p}_{T};#eta;#phi", HistType::kTH3F, {trackPtAxis, etaAxis, phiAxis});
    registry.add("jets/detJetPtTrackProj", "Detector level jet #it{p}_{T}, #it{z}; #it{p}_{T}; #it{z}", HistType::kTH2F, {jetPtAxis, zAxis});
    registry.add("jets/detJetPtFrag", "Detector level jet #it{p}_{T}, #it{p}_{T,jet}/#it{p}_{T,tr}; #it{p}_{T}; #it{p}_{T,jet}/#it{p}_{T,tr}", HistType::kTH2F, {jetPtAxis, zAxis});

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

  void processMCD(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                  soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets,
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
  PROCESS_SWITCH(JetFragmentation, processMCD, "Monte Carlo detector level", false);

  void processMCP(aod::McCollision const& mcCollision, // Add some form of event selection?
                  soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
                  aod::McParticles const& particles)
  {
    int nJets = -1, nTracks = -1;
    nJets = jets.size(), nTracks = particles.size();
    registry.fill(HIST("partnJetnTrack"), nJets, nTracks);
    for (const auto& particle : particles) {
      registry.fill(HIST("tracks/partTrackPtEtaPhi"), particle.pt(), particle.eta(), particle.phi());
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
  PROCESS_SWITCH(JetFragmentation, processMCP, "Monte Carlo particle level", false);

  void processDataRun3(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                       soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets,
                       aod::Tracks const& tracks)
  {
    registry.fill(HIST("datCount"), 1);
    int nJets = -1, nTracks = -1;
    nJets = jets.size(), nTracks = tracks.size();
    registry.fill(HIST("partnJetnTrack"), nJets, nTracks);
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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetFragmentation>(cfgc, TaskName{"jet-fragmentation"})};
}
