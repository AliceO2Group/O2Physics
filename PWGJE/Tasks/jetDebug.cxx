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

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct JetDebugTask {
  HistogramRegistry registry{"registry",
                             {{"h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
                              {"h_collisions_selected_posz", "event z;event z;entries", {HistType::kTH1F, {{80, -20.0, 20.0}}}},
                              {"h_collisions_selected_posz_original", "event z;event z;entries", {HistType::kTH1F, {{80, -20.0, 20.0}}}},
                              {"h_mccollisions", "mc event status;event status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
                              {"h_tracks_selected_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_tracks_selected_pt_original", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_track_eta", "track #eta;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_track_phi", "track #varphi;#varphi_{track};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_particle_pt", "particle pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_particle_eta", "particle #eta;#eta_{track};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_particle_phi", "particle #varphi;#varphi_{track};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_jet_pt", "track pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_jet_phi", "jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_jet_nTracks", "jet tracks;n_{track};entries", {HistType::kTH1F, {{80, 0.0, 80.0}}}},
                              {"h_jet_angularity", "jet angularity;angularity;entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
                              {"h_mcpjet_pt", "mcp jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_mcpjet_eta", "mcp jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_mcpjet_phi", "mcp jet #varphi;#varphi_{jeet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_mcpjet_nTracks", "mcp jet tracks,n_{track};entries", {HistType::kTH1F, {{80, 0.0, 80.0}}}},
                              {"h_mcpjet_angularity", "mcp jet angularity;angularity;entries", {HistType::kTH1F, {{50, 0.0, 1.0}}}},
                              {"h_jet_pt_matched", "jet pT matched;#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}}},
                              {"h_D0_pt", "D0 pT;#it{p}_{T,D0} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_D0injet_pt", "D0 in jet pT;#it{p}_{T,D0} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_D0jet_pt", "D0 jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_D0jet_eta", "D0 jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_D0jet_phi", "D0 jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_D0jet_nTracks", "D0 jet n tracks;n_{track};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_D0jet_angularity", "D0 jet angularity;angularity;entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_D0jet_pt_matched", "D0 jet pT matched;#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}}},
                              {"h_mcpD0injet_pt", "mcp D0 in jet pT;#it{p}_{T,D0} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_mcpD0_pt", "mcp D0 pT;#it{p}_{T,D0} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_mcpD0jet_pt", "mcp D0 jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_mcpD0jet_eta", "mcp D0 jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_mcpD0jet_phi", "mcp D0 jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_mcpD0jet_nTracks", "mcp D0 jet ntracks;n_{track};entries", {HistType::kTH1F, {{80, 0.0, 80.}}}},
                              {"h_mcpD0jet_angularity", "mcp D0 jet angularity;angularity;entries", {HistType::kTH1F, {{50, 0.0, 1.}}}},
                              {"h_Lc_pt", "Lc pT;#it{p}_{T,Lc} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_Lcinjet_pt", "Lc in jet pT;#it{p}_{T,Lc} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_Lcjet_pt", "Lc jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_Lcjet_eta", "Lc jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_Lcjet_phi", "Lc jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_Lcjet_nTracks", "Lc jet n tracks;n_{track};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_Lcjet_angularity", "Lc jet angularity;angularity;entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_Lcjet_pt_matched", "Lc jet pT matched;#it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}}},
                              {"h_mcpLc_pt", "mcp Lc pT;#it{p}_{T,Lc} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_mcpLcinjet_pt", "mcp Lc in jet pT;#it{p}_{T,Lc} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_mcpLcjet_pt", "mcp Lc jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_mcpLcjet_eta", "mcp Lc jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_mcpLcjet_phi", "mcp Lc jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_mcpLcjet_nTracks", "mcp Lc jet ntracks;n_{track};entries", {HistType::kTH1F, {{80, 0.0, 80.}}}},
                              {"h_mcpLcjet_angularity", "mcp Lc jet angularity;angularity;entries", {HistType::kTH1F, {{50, 0.0, 1.}}}}}};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  std::vector<int> eventSelection;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
  }

  void processCollisionSelection(aod::JetCollisions const& collisions, soa::Join<aod::Collisions, aod::EvSels> const& originalCollisions)
  {
    for (auto const& collision : collisions) {
      if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>("sel8")), false, false)) {
        registry.fill(HIST("h_collisions_selected_posz"), collision.posZ());
      }
    }
    for (auto const& originalCollision : originalCollisions) {
      if (originalCollision.sel8()) {
        registry.fill(HIST("h_collisions_selected_posz_original"), originalCollision.posZ());
      }
    }
  }
  PROCESS_SWITCH(JetDebugTask, processCollisionSelection, "collision selection", true);

  void processTrackSelection(aod::JetTracks const& tracks, soa::Join<aod::Tracks, aod::TrackSelection> const& originalTracks)
  {
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>("globalTracks")))) {
        registry.fill(HIST("h_tracks_selected_pt"), track.pt());
      }
    }
    for (auto const& originalTrack : originalTracks) {
      if (originalTrack.isGlobalTrackWoPtEta()) {
        registry.fill(HIST("h_tracks_selected_pt_original"), originalTrack.pt());
      }
    }
  }
  PROCESS_SWITCH(JetDebugTask, processTrackSelection, "track selection", true);

  void processDataCharged(aod::JetCollision const& collision, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets, aod::JTracks const& tracks)
  {
    registry.fill(HIST("h_collisions"), 1.0);
    if (jetderiveddatautilities::selectCollision(collision, eventSelection, false, false)) {
      registry.fill(HIST("h_collisions"), 2.0);
    }
    for (auto track : tracks) {
      registry.fill(HIST("h_track_pt"), track.pt());
      registry.fill(HIST("h_track_phi"), track.phi());
      registry.fill(HIST("h_track_eta"), track.eta());
    }
    for (auto const& jet : jets) {
      registry.fill(HIST("h_jet_pt"), jet.pt());
      registry.fill(HIST("h_jet_phi"), jet.phi());
      registry.fill(HIST("h_jet_eta"), jet.eta());
      registry.fill(HIST("h_jet_nTracks"), jet.tracksIds().size());
      double angularity = 0.0;
      for (auto const& jetConstituent : jet.tracks_as<aod::JetTracks>()) {
        angularity += (jetConstituent.pt() * TMath::Sqrt(TMath::Power(RecoDecay::constrainAngle(jet.phi() - jetConstituent.phi(), -M_PI), 2.0) + TMath::Power(jet.eta() - jetConstituent.eta(), 2.0))) / (jet.pt() * (jet.r() / 100.f));
      }
      registry.fill(HIST("h_jet_angularity"), angularity);
    }
  }
  PROCESS_SWITCH(JetDebugTask, processDataCharged, "jets data", true);

  void processDataD0(aod::JetCollision const&, soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents> const& jets, aod::JTracks const&, aod::CandidatesD0Data const& candidates)
  {
    for (auto const& candidate : candidates) {
      registry.fill(HIST("h_D0_pt"), candidate.pt());
    }
    for (auto const& jet : jets) {
      registry.fill(HIST("h_D0jet_pt"), jet.pt());
      registry.fill(HIST("h_D0jet_phi"), jet.phi());
      registry.fill(HIST("h_D0jet_eta"), jet.eta());
      registry.fill(HIST("h_D0jet_nTracks"), jet.tracksIds().size() + jet.candidatesIds().size());
      double angularity = 0.0;
      for (auto const& jetConstituent : jet.tracks_as<aod::JetTracks>()) {
        angularity += (jetConstituent.pt() * TMath::Sqrt(TMath::Power(RecoDecay::constrainAngle(jet.phi() - jetConstituent.phi(), -M_PI), 2.0) + TMath::Power(jet.eta() - jetConstituent.eta(), 2.0))) / (jet.pt() * (jet.r() / 100.f));
      }
      for (auto const& jetCandidate : jet.candidates_as<aod::CandidatesD0Data>()) {
        registry.fill(HIST("h_D0injet_pt"), jetCandidate.pt());
        angularity += (jetCandidate.pt() * TMath::Sqrt(TMath::Power(RecoDecay::constrainAngle(jet.phi() - jetCandidate.phi(), -M_PI), 2.0) + TMath::Power(jet.eta() - jetCandidate.eta(), 2.0))) / (jet.pt() * (jet.r() / 100.f));
      }
      registry.fill(HIST("h_D0jet_angularity"), angularity);
    }
  }
  PROCESS_SWITCH(JetDebugTask, processDataD0, "jets D0 data", true);

  void processDataLc(aod::JetCollision const&, soa::Join<aod::LcChargedJets, aod::LcChargedJetConstituents> const& jets, aod::JTracks const&, aod::CandidatesLcData const& candidates)
  {
    for (auto const& candidate : candidates) {
      registry.fill(HIST("h_Lc_pt"), candidate.pt());
    }
    for (auto const& jet : jets) {
      registry.fill(HIST("h_Lcjet_pt"), jet.pt());
      registry.fill(HIST("h_Lcjet_phi"), jet.phi());
      registry.fill(HIST("h_Lcjet_eta"), jet.eta());
      registry.fill(HIST("h_Lcjet_nTracks"), jet.tracksIds().size() + jet.candidatesIds().size());
      double angularity = 0.0;
      for (auto const& jetConstituent : jet.tracks_as<aod::JetTracks>()) {
        angularity += (jetConstituent.pt() * TMath::Sqrt(TMath::Power(RecoDecay::constrainAngle(jet.phi() - jetConstituent.phi(), -M_PI), 2.0) + TMath::Power(jet.eta() - jetConstituent.eta(), 2.0))) / (jet.pt() * (jet.r() / 100.f));
      }
      for (auto const& jetCandidate : jet.candidates_as<aod::CandidatesLcData>()) {
        registry.fill(HIST("h_Lcinjet_pt"), jetCandidate.pt());
        angularity += (jetCandidate.pt() * TMath::Sqrt(TMath::Power(RecoDecay::constrainAngle(jet.phi() - jetCandidate.phi(), -M_PI), 2.0) + TMath::Power(jet.eta() - jetCandidate.eta(), 2.0))) / (jet.pt() * (jet.r() / 100.f));
      }
      registry.fill(HIST("h_Lcjet_angularity"), angularity);
    }
  }
  PROCESS_SWITCH(JetDebugTask, processDataLc, "jets Lc data", true);

  void processMCDCharged(aod::JetCollision const& collision, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets> const& jets, aod::JTracks const& tracks, aod::ChargedMCParticleLevelJets const&)
  {
    registry.fill(HIST("h_collisions"), 1.0);
    if (jetderiveddatautilities::selectCollision(collision, eventSelection, false, false)) {
      registry.fill(HIST("h_collisions"), 2.0);
    }
    for (auto track : tracks) {
      registry.fill(HIST("h_track_pt"), track.pt());
      registry.fill(HIST("h_track_phi"), track.phi());
      registry.fill(HIST("h_track_eta"), track.eta());
    }
    for (auto const& jet : jets) {
      registry.fill(HIST("h_jet_pt"), jet.pt());
      registry.fill(HIST("h_jet_phi"), jet.phi());
      registry.fill(HIST("h_jet_eta"), jet.eta());
      registry.fill(HIST("h_jet_nTracks"), jet.tracksIds().size());
      auto const& matchedJets = jet.matchedJetPt_as<aod::ChargedMCParticleLevelJets>();
      for (auto const& mcpjet : matchedJets) {
        registry.fill(HIST("h_jet_pt_matched"), mcpjet.pt(), jet.pt());
      }
      double angularity = 0.0;
      for (auto const& jetConstituent : jet.tracks_as<aod::JetTracks>()) {
        angularity += (jetConstituent.pt() * TMath::Sqrt(TMath::Power(RecoDecay::constrainAngle(jet.phi() - jetConstituent.phi(), -M_PI), 2.0) + TMath::Power(jet.eta() - jetConstituent.eta(), 2.0))) / (jet.pt() * (jet.r() / 100.f));
      }
      if (std::abs(jet.eta()) < 0.5) {
        registry.fill(HIST("h_jet_angularity"), angularity);
      }
    }
  }
  PROCESS_SWITCH(JetDebugTask, processMCDCharged, "jets mcd", false);

  void processMCDD0(aod::JetCollision const&, soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets> const& jets, aod::JTracks const&, aod::CandidatesD0MCD const& candidates, aod::D0ChargedMCParticleLevelJets const&)
  {
    for (auto const& candidate : candidates) {
      registry.fill(HIST("h_D0_pt"), candidate.pt());
    }
    for (auto const& jet : jets) {
      registry.fill(HIST("h_D0jet_pt"), jet.pt());
      registry.fill(HIST("h_D0jet_phi"), jet.phi());
      registry.fill(HIST("h_D0jet_eta"), jet.eta());
      registry.fill(HIST("h_D0jet_nTracks"), jet.tracksIds().size() + jet.candidatesIds().size());
      auto const& matchedJets = jet.matchedJetCand_as<aod::D0ChargedMCParticleLevelJets>();
      for (auto const& mcpjet : matchedJets) {
        registry.fill(HIST("h_D0jet_pt_matched"), mcpjet.pt(), jet.pt());
      }
      double angularity = 0.0;
      for (auto const& jetConstituent : jet.tracks_as<aod::JetTracks>()) {
        angularity += (jetConstituent.pt() * TMath::Sqrt(TMath::Power(RecoDecay::constrainAngle(jet.phi() - jetConstituent.phi(), -M_PI), 2.0) + TMath::Power(jet.eta() - jetConstituent.eta(), 2.0))) / (jet.pt() * (jet.r() / 100.f));
      }
      for (auto const& jetCandidate : jet.candidates_as<aod::CandidatesD0MCD>()) {
        registry.fill(HIST("h_D0injet_pt"), jetCandidate.pt());
        angularity += (jetCandidate.pt() * TMath::Sqrt(TMath::Power(RecoDecay::constrainAngle(jet.phi() - jetCandidate.phi(), -M_PI), 2.0) + TMath::Power(jet.eta() - jetCandidate.eta(), 2.0))) / (jet.pt() * (jet.r() / 100.f));
      }
      registry.fill(HIST("h_D0jet_angularity"), angularity);
    }
  }
  PROCESS_SWITCH(JetDebugTask, processMCDD0, "jets D0 mcd", false);

  void processMCDLc(aod::JetCollision const&, soa::Join<aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents, aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCParticleLevelJets> const& jets, aod::JTracks const&, aod::CandidatesLcMCD const& candidates, aod::LcChargedMCParticleLevelJets const&)
  {
    for (auto const& candidate : candidates) {
      registry.fill(HIST("h_Lc_pt"), candidate.pt());
    }
    for (auto const& jet : jets) {
      registry.fill(HIST("h_Lcjet_pt"), jet.pt());
      registry.fill(HIST("h_Lcjet_phi"), jet.phi());
      registry.fill(HIST("h_Lcjet_eta"), jet.eta());
      registry.fill(HIST("h_Lcjet_nTracks"), jet.tracksIds().size() + jet.candidatesIds().size());
      auto const& matchedJets = jet.matchedJetCand_as<aod::LcChargedMCParticleLevelJets>();
      for (auto const& mcpjet : matchedJets) {
        registry.fill(HIST("h_Lcjet_pt_matched"), mcpjet.pt(), jet.pt());
      }
      double angularity = 0.0;
      for (auto const& jetConstituent : jet.tracks_as<aod::JetTracks>()) {
        angularity += (jetConstituent.pt() * TMath::Sqrt(TMath::Power(RecoDecay::constrainAngle(jet.phi() - jetConstituent.phi(), -M_PI), 2.0) + TMath::Power(jet.eta() - jetConstituent.eta(), 2.0))) / (jet.pt() * (jet.r() / 100.f));
      }
      for (auto const& jetCandidate : jet.candidates_as<aod::CandidatesLcMCD>()) {
        registry.fill(HIST("h_Lcinjet_pt"), jetCandidate.pt());
        angularity += (jetCandidate.pt() * TMath::Sqrt(TMath::Power(RecoDecay::constrainAngle(jet.phi() - jetCandidate.phi(), -M_PI), 2.0) + TMath::Power(jet.eta() - jetCandidate.eta(), 2.0))) / (jet.pt() * (jet.r() / 100.f));
      }
      registry.fill(HIST("h_Lcjet_angularity"), angularity);
    }
  }
  PROCESS_SWITCH(JetDebugTask, processMCDLc, "jets Lc mcd", false);

  void processMCPCharged(aod::JMcCollision const&, soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets, aod::JMcParticles const& tracks)
  {
    for (auto track : tracks) {
      registry.fill(HIST("h_particle_pt"), track.pt());
      registry.fill(HIST("h_particle_phi"), track.phi());
      registry.fill(HIST("h_particle_eta"), track.eta());
    }
    for (auto const& jet : jets) {
      registry.fill(HIST("h_mcpjet_pt"), jet.pt());
      registry.fill(HIST("h_mcpjet_phi"), jet.phi());
      registry.fill(HIST("h_mcpjet_eta"), jet.eta());
      registry.fill(HIST("h_mcpjet_nTracks"), jet.tracksIds().size());
      double angularity = 0.0;
      for (auto const& jetConstituent : jet.tracks_as<aod::JMcParticles>()) {
        angularity += (jetConstituent.pt() * TMath::Sqrt(TMath::Power(RecoDecay::constrainAngle(jet.phi() - jetConstituent.phi(), -M_PI), 2.0) + TMath::Power(jet.eta() - jetConstituent.eta(), 2.0))) / (jet.pt() * (jet.r() / 100.f));
      }
      registry.fill(HIST("h_mcpjet_angularity"), angularity);
    }
  }
  PROCESS_SWITCH(JetDebugTask, processMCPCharged, "jets mcp", false);

  void processMCPD0(aod::JMcCollision const&, soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents> const& jets, aod::JMcParticles const&, aod::CandidatesD0MCP const& candidates)
  {
    for (auto const& candidate : candidates) {
      registry.fill(HIST("h_mcpD0_pt"), candidate.pt());
    }
    for (auto const& jet : jets) {
      registry.fill(HIST("h_mcpD0jet_pt"), jet.pt());
      registry.fill(HIST("h_mcpD0jet_phi"), jet.phi());
      registry.fill(HIST("h_mcpD0jet_eta"), jet.eta());
      registry.fill(HIST("h_mcpD0jet_nTracks"), jet.tracksIds().size() + jet.candidatesIds().size());
      double angularity = 0.0;
      for (auto const& jetConstituent : jet.tracks_as<aod::JMcParticles>()) {
        angularity += (jetConstituent.pt() * TMath::Sqrt(TMath::Power(RecoDecay::constrainAngle(jet.phi() - jetConstituent.phi(), -M_PI), 2.0) + TMath::Power(jet.eta() - jetConstituent.eta(), 2.0))) / (jet.pt() * (jet.r() / 100.f));
      }
      for (auto const& jetCandidate : jet.candidates_as<aod::CandidatesD0MCP>()) {
        registry.fill(HIST("h_mcpD0injet_pt"), jetCandidate.pt());
        angularity += (jetCandidate.pt() * TMath::Sqrt(TMath::Power(RecoDecay::constrainAngle(jet.phi() - jetCandidate.phi(), -M_PI), 2.0) + TMath::Power(jet.eta() - jetCandidate.eta(), 2.0))) / (jet.pt() * (jet.r() / 100.f));
      }
      registry.fill(HIST("h_mcpD0jet_angularity"), angularity);
    }
  }
  PROCESS_SWITCH(JetDebugTask, processMCPD0, "jets D0 mcp", false);

  void processMCPLc(aod::JMcCollision const&, soa::Join<aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents> const& jets, aod::JMcParticles const&, aod::CandidatesLcMCP const& candidates)
  {
    for (auto const& candidate : candidates) {
      registry.fill(HIST("h_mcpLc_pt"), candidate.pt());
    }
    for (auto const& jet : jets) {
      registry.fill(HIST("h_mcpLcjet_pt"), jet.pt());
      registry.fill(HIST("h_mcpLcjet_phi"), jet.phi());
      registry.fill(HIST("h_mcpLcjet_eta"), jet.eta());
      registry.fill(HIST("h_mcpLcjet_nTracks"), jet.tracksIds().size() + jet.candidatesIds().size());
      double angularity = 0.0;
      for (auto const& jetConstituent : jet.tracks_as<aod::JMcParticles>()) {
        angularity += (jetConstituent.pt() * TMath::Sqrt(TMath::Power(RecoDecay::constrainAngle(jet.phi() - jetConstituent.phi(), -M_PI), 2.0) + TMath::Power(jet.eta() - jetConstituent.eta(), 2.0))) / (jet.pt() * (jet.r() / 100.f));
        for (auto const& jetCandidate : jet.candidates_as<aod::CandidatesLcMCP>()) {
          registry.fill(HIST("h_mcpLcinjet_pt"), jetCandidate.pt());
          angularity += (jetCandidate.pt() * TMath::Sqrt(TMath::Power(RecoDecay::constrainAngle(jet.phi() - jetCandidate.phi(), -M_PI), 2.0) + TMath::Power(jet.eta() - jetCandidate.eta(), 2.0))) / (jet.pt() * (jet.r() / 100.f));
        }
        registry.fill(HIST("h_mcpLcjet_angularity"), angularity);
      }
    }
  }
  PROCESS_SWITCH(JetDebugTask, processMCPLc, "jets Lc mcp", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetDebugTask>(cfgc, TaskName{"jet-debug"})}; }
