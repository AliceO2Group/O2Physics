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

// jet tutorial task for hands on tutorial session (09/11/2023)
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
//

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct JetTutorialTask {
  HistogramRegistry registry{"registry",
                             {{"h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, 0.0, 4.0}}}},
                              {"h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_track_eta", "track #eta;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_track_phi", "track #varphi;#varphi_{track};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_track_chi2PerCluster", "track #chi^{2} per cluste ;#chi^{2};entries", {HistType::kTH1F, {{100, 0, 40}}}},
                              {"h_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_jet_phi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_jet_pt_bkgsub", "jet pT bkg sub;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_part_jet_pt", "particle level jet pT;#it{p}_{T,jet part} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_part_jet_eta", "particle level jet #eta;#eta_{jet part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_part_jet_phi", "particle level jet #phi;#phi_{jet part};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_jet_ntracks", "jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {{40, -0.5, 39.5}}}},
                              {"h_jet_angularity", "jet angularity ;#lambda_{1};entries", {HistType::kTH1F, {{5, 0.0, 0.5}}}},
                              {"h_full_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_full_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_full_jet_phi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_full_jet_ntracks", "jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {{40, -0.5, 39.5}}}},
                              {"h_full_jet_nclusters", "jet N clusters;N_{jet clusters};entries", {HistType::kTH1F, {{40, -0.5, 39.5}}}},
                              {"h_full_jet_angularity", "jet angularity ;#lambda_{1};entries", {HistType::kTH1F, {{5, 0.0, 0.5}}}},
                              {"h_part_jet_angularity", "jet angularity ;#lambda_{1};entries", {HistType::kTH1F, {{5, 0.0, 0.5}}}},
                              {"h_recoil_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_recoil_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_recoil_jet_phi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_recoil_jet_dphi", "hadron-jet #Delta#phi;#Delta#phi_{jet,trigger hadron};entries", {HistType::kTH1F, {{40, -2.0, 2.0}}}},
                              {"h_matched_jets_pt", "#it{p}_{T,jet part}; #it{p}_{T,jet det}", {HistType::kTH2F, {{100, 0., 20.}, {100, 0., 20.0}}}},
                              {"h_matched_jets_eta", "#eta_{jet part}; #eta_{jet det}", {HistType::kTH2F, {{100, -1.0, 1.0}, {100, -1.0, 1.0}}}},
                              {"h_matched_jets_phi", "#phi_{jet part}; #phi_{jet det}", {HistType::kTH2F, {{80, -1.0, 7.}, {80, -1.0, 7.}}}}}};

  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  int eventSelection = -1;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    eventSelection = JetDerivedDataUtilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = JetDerivedDataUtilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
  }

  Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f);

  void processCollisions(aod::JCollision const& collision, aod::JTracks const& tracks)
  {

    registry.fill(HIST("h_collisions"), 0.5);
    if (!JetDerivedDataUtilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    for (auto const& track : tracks) {
      if (!JetDerivedDataUtilities::selectTrack(track, trackSelection)) {
        continue;
      }
      registry.fill(HIST("h_track_pt"), track.pt());
      registry.fill(HIST("h_track_eta"), track.eta());
      registry.fill(HIST("h_track_phi"), track.phi());
    }
  }
  PROCESS_SWITCH(JetTutorialTask, processCollisions, "process self contained collisions", true);

  void processCollisionsWithExternalTracks(aod::JCollision const& collision, soa::Join<aod::JTracks, aod::JTrackPIs> const& tracks, soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection> const& originalTracks)
  {

    registry.fill(HIST("h_collisions"), 0.5);
    if (!JetDerivedDataUtilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 1.5);
    for (auto const& track : tracks) {
      if (!JetDerivedDataUtilities::selectTrack(track, trackSelection)) {
        continue;
      }
      registry.fill(HIST("h_track_pt"), track.pt());
      registry.fill(HIST("h_track_eta"), track.eta());
      registry.fill(HIST("h_track_phi"), track.phi());
      auto originalTrack = track.track_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>();
      registry.fill(HIST("h_track_chi2PerCluster"), originalTrack.tpcChi2NCl());
    }
  }
  PROCESS_SWITCH(JetTutorialTask, processCollisionsWithExternalTracks, "process non self contained collisions", true);

  void processDataCharged(soa::Filtered<aod::ChargedJets>::iterator const& jet)
  {
    registry.fill(HIST("h_jet_pt"), jet.pt());
    registry.fill(HIST("h_jet_eta"), jet.eta());
    registry.fill(HIST("h_jet_phi"), jet.phi());
  }
  PROCESS_SWITCH(JetTutorialTask, processDataCharged, "jets data", true);

  void processMCDetectorLevelCharged(soa::Filtered<aod::ChargedMCDetectorLevelJets>::iterator const& jet)
  {
    registry.fill(HIST("h_jet_pt"), jet.pt());
    registry.fill(HIST("h_jet_eta"), jet.eta());
    registry.fill(HIST("h_jet_phi"), jet.phi());
  }
  PROCESS_SWITCH(JetTutorialTask, processMCDetectorLevelCharged, "jets on detector level MC", false);

  void processMCParticleLevel(soa::Filtered<aod::ChargedMCParticleLevelJets>::iterator const& jet)
  {
    registry.fill(HIST("h_part_jet_pt"), jet.pt());
    registry.fill(HIST("h_part_jet_eta"), jet.eta());
    registry.fill(HIST("h_part_jet_phi"), jet.phi());
  }
  PROCESS_SWITCH(JetTutorialTask, processMCParticleLevel, "jets on particle level MC", false);

  void processMCCharged(aod::JCollisions const& collision, soa::Filtered<aod::ChargedMCDetectorLevelJets> const& mcdjets, soa::Filtered<aod::ChargedMCParticleLevelJets> const& mcpjets)
  {
    for (auto& mcdjet : mcdjets) {
      registry.fill(HIST("h_jet_pt"), mcdjet.pt());
      registry.fill(HIST("h_jet_eta"), mcdjet.eta());
      registry.fill(HIST("h_jet_phi"), mcdjet.phi());
    }
    for (auto& mcpjet : mcpjets) {
      registry.fill(HIST("h_part_jet_pt"), mcpjet.pt());
      registry.fill(HIST("h_part_jet_eta"), mcpjet.eta());
      registry.fill(HIST("h_part_jet_phi"), mcpjet.phi());
    }
  }
  PROCESS_SWITCH(JetTutorialTask, processMCCharged, "jets on detector and particle level MC", false);

  using JetMCPTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;
  void processMCChargedMatched(aod::JCollision const& collision,
                               soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>> const& mcdjets,
                               JetMCPTable const& mcpjets,
                               aod::JTracks const& tracks,
                               aod::JMcParticles const& particles)
  {
    for (const auto& mcdjet : mcdjets) {

      for (auto& mcpjet : mcdjet.template matchedJetGeo_as<JetMCPTable>()) {
        // for (auto& mcpjet : mcdjet.template matchedJetPt_as<JetMCPTable>()) {

        registry.fill(HIST("h_matched_jets_pt"), mcpjet.pt(), mcdjet.pt());
        registry.fill(HIST("h_matched_jets_pt"), mcpjet.phi(), mcdjet.phi());
        registry.fill(HIST("h_matched_jets_pt"), mcpjet.eta(), mcdjet.eta());
      }
    }
  }
  PROCESS_SWITCH(JetTutorialTask, processMCChargedMatched, "jet finder QA matched mcp and mcd", false);

  void processDataChargedSubstructure(soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>>::iterator const& jet, aod::JTracks const& tracks)
  {
    // add aditional selection on jet eta
    registry.fill(HIST("h_jet_pt"), jet.pt());
    registry.fill(HIST("h_jet_eta"), jet.eta());
    registry.fill(HIST("h_jet_phi"), jet.phi());
    registry.fill(HIST("h_jet_ntracks"), jet.tracks().size());
    double angularity = 0.0;
    for (auto& jetConstituent : jet.tracks_as<aod::JTracks>()) {
      angularity += jetConstituent.pt() * TMath::Sqrt(TMath::Power(jet.phi() - jetConstituent.phi(), 2.0) + TMath::Power(jet.eta() - jetConstituent.eta(), 2.0));
    }
    registry.fill(HIST("h_jet_angularity"), angularity / (jet.pt() * round(jet.r() / 100.0f)));
  }
  PROCESS_SWITCH(JetTutorialTask, processDataChargedSubstructure, "jet substructure charged jets", false);

  void processMCParticleSubstructure(soa::Filtered<soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents>>::iterator const& jet, aod::JMcParticles const& particles)
  {
    double angularity = 0.0;
    for (auto& jetConstituents : jet.tracks_as<aod::JMcParticles>()) {
      angularity += jetConstituents.pt() * TMath::Sqrt(TMath::Power(jet.phi() - jetConstituents.phi(), 2.0) + TMath::Power(jet.eta() - jetConstituents.eta(), 2.0));
    }
    registry.fill(HIST("h_part_jet_angularity"), angularity / (jet.pt() * round(jet.r() / 100.0f)));
  }
  PROCESS_SWITCH(JetTutorialTask, processMCParticleSubstructure, "jet substructure particle level full jets", false);

  void processDataFull(soa::Filtered<aod::FullJets>::iterator const& jet)
  {
    registry.fill(HIST("h_jet_pt"), jet.pt());
    registry.fill(HIST("h_jet_eta"), jet.eta());
    registry.fill(HIST("h_jet_phi"), jet.phi());
  }
  PROCESS_SWITCH(JetTutorialTask, processDataFull, "jets data", true);

  void processDataFullSubstructure(soa::Filtered<soa::Join<aod::FullJets, aod::FullJetConstituents>>::iterator const& jet, aod::JTracks const& tracks, aod::JClusters const& clusters)
  {
    // add aditional selection on jet eta
    registry.fill(HIST("h_full_jet_pt"), jet.pt());
    registry.fill(HIST("h_full_jet_eta"), jet.eta());
    registry.fill(HIST("h_full_jet_phi"), jet.phi());
    registry.fill(HIST("h_full_jet_ntracks"), jet.tracks().size());
    registry.fill(HIST("h_full_jet_nclusters"), jet.clusters().size());
    double angularity = 0.0;
    for (auto& jetTrack : jet.tracks_as<aod::JTracks>()) {
      angularity += jetTrack.pt() * TMath::Sqrt(TMath::Power(jet.phi() - jetTrack.phi(), 2.0) + TMath::Power(jet.eta() - jetTrack.eta(), 2.0));
    }

    for (auto& jetCluster : jet.clusters_as<aod::JClusters>()) {
      angularity += jetCluster.energy() * TMath::Sqrt(TMath::Power(jet.phi() - jetCluster.phi(), 2.0) + TMath::Power(jet.eta() - jetCluster.eta(), 2.0));
    }

    registry.fill(HIST("h_full_jet_angularity"), angularity / (jet.pt() * round(jet.r() * 100.0f)));
  }
  PROCESS_SWITCH(JetTutorialTask, processDataFullSubstructure, "jet substructure full jets", false);

  void processDataRecoil(aod::JCollision const& collision, soa::Filtered<aod::ChargedJets> const& jets, aod::JTracks const& tracks)
  {
    if (!JetDerivedDataUtilities::selectCollision(collision, eventSelection)) {
      return;
    }
    double leadingTrackpT = 0.0;
    double leadingTrackPhi = 0.0;
    for (auto& track : tracks) {
      if (track.pt() > 6.0 && track.pt() < 10.0) {
        if (track.pt() > leadingTrackpT) {
          leadingTrackpT = track.pt();
          leadingTrackPhi = track.phi();
        }
      }
    }
    if (leadingTrackpT == 0.0)
      return;
    for (auto& jet : jets) {
      if (TMath::Abs(RecoDecay::constrainAngle(RecoDecay::constrainAngle(jet.phi(), -o2::constants::math::PIHalf) - RecoDecay::constrainAngle(leadingTrackPhi, -o2::constants::math::PIHalf), -o2::constants::math::PIHalf) > 0.6)) {
        registry.fill(HIST("h_recoil_jet_pt"), jet.pt());
        registry.fill(HIST("h_recoil_jet_eta"), jet.eta());
        registry.fill(HIST("h_recoil_jet_phi"), jet.phi());
        registry.fill(HIST("h_recoil_jet_dphi"), jet.phi() - leadingTrackPhi);
      }
    }
  }
  PROCESS_SWITCH(JetTutorialTask, processDataRecoil, "hadron-recoil jets", false);

  void processDataBackgroundSubtracted(soa::Join<aod::JCollisions, aod::JCollisionRhos>::iterator const& collision, soa::Filtered<aod::ChargedJets> const& jets)
  {
    for (auto jet : jets) {
      registry.fill(HIST("h_jet_pt"), jet.pt());
      registry.fill(HIST("h_jet_pt_bkgsub"), jet.pt() - (collision.rho() * jet.area()));
      registry.fill(HIST("h_jet_eta"), jet.eta());
      registry.fill(HIST("h_jet_phi"), jet.phi());
    }
  }
  PROCESS_SWITCH(JetTutorialTask, processDataBackgroundSubtracted, "baackground subtracted jets", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetTutorialTask>(cfgc, TaskName{"jet-tutorial"})}; }
