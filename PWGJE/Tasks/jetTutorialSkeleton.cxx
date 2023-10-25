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

// jet tutorial task for hands on tutorial session (27/04/2023)
//
// Author: Nima Zardoshti
//

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "PWGJE/DataModel/EMCALClusters.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct JetTutorialSkeletonTask {
  HistogramRegistry registry{"registry",
                             {{"h_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 50.}}}},
                              {"h_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{30, -1.5, 1.5}}}},
                              {"h_jet_phi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{140, -7.0, 7.}}}},
                              {"h_part_jet_pt", "particle level jet pT;#it{p}_{T,jet part} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 50.}}}},
                              {"h_part_jet_eta", "particle level jet #eta;#eta_{jet part};entries", {HistType::kTH1F, {{30, -1.5, 1.5}}}},
                              {"h_part_jet_phi", "particle level jet #phi;#phi_{jet part};entries", {HistType::kTH1F, {{140, -7.0, 7.}}}},
                              {"h_jet_ntracks", "jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {{40, -0.5, 39.5}}}},
                              {"h_jet_angularity", "jet angularity ;#lambda_{1};entries", {HistType::kTH1F, {{5, 0.0, 0.5}}}},
                              {"h_full_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 50.}}}},
                              {"h_full_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{30, -1.5, 1.5}}}},
                              {"h_full_jet_phi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{140, -7.0, 7.}}}},
                              {"h_full_jet_ntracks", "jet N tracks;N_{jet tracks};entries", {HistType::kTH1F, {{40, -0.5, 39.5}}}},
                              {"h_full_jet_nclusters", "jet N clusters;N_{jet clusters};entries", {HistType::kTH1F, {{40, -0.5, 39.5}}}},
                              {"h_full_jet_angularity", "jet angularity ;#lambda_{1};entries", {HistType::kTH1F, {{5, 0.0, 0.5}}}},
                              {"h_part_jet_angularity", "jet angularity ;#lambda_{1};entries", {HistType::kTH1F, {{5, 0.0, 0.5}}}},
                              {"h_recoil_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 50.}}}},
                              {"h_recoil_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{30, -1.5, 1.5}}}},
                              {"h_recoil_jet_phi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{140, -7.0, 7.}}}},
                              {"h_recoil_jet_dphi", "hadron-jet #Delta#phi;#Delta#phi_{jet,trigger hadron};entries", {HistType::kTH1F, {{40, -2.0, 2.0}}}},
                              {"h_matched_jets_pt", "#it{p}_{T,jet part}; #it{p}_{T,jet det}", {HistType::kTH2F, {{100, 0., 20.}, {100, 0., 20.0}}}},
                              {"h_matched_jets_eta", "#eta_{jet part}; #eta_{jet det}", {HistType::kTH2F, {{30, -1.5, 1.5}, {30, -1.5, 1.5}}}},
                              {"h_matched_jets_phi", "#phi_{jet part}; #phi_{jet det}", {HistType::kTH2F, {{140, -7.0, 7.}, {140, -7.0, 7.}}}}}};

  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};

  void init(InitContext const&) {}

  Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f);

  void processDataCharged(soa::Filtered<aod::ChargedJets>::iterator const& jet)
  {
    registry.fill(HIST("h_jet_pt"), jet.pt());
    registry.fill(HIST("h_jet_eta"), jet.eta());
    registry.fill(HIST("h_jet_phi"), jet.phi());
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processDataCharged, "jets data", true);

  void processMCDetectorLevelCharged(soa::Filtered<aod::ChargedMCDetectorLevelJets>::iterator const& jet)
  {
    registry.fill(HIST("h_jet_pt"), jet.pt());
    registry.fill(HIST("h_jet_eta"), jet.eta());
    registry.fill(HIST("h_jet_phi"), jet.phi());
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processMCDetectorLevelCharged, "jets on detector level MC", false);

  void processMCParticleLevel(soa::Filtered<aod::ChargedMCParticleLevelJets>::iterator const& jet)
  {
    registry.fill(HIST("h_part_jet_pt"), jet.pt());
    registry.fill(HIST("h_part_jet_eta"), jet.eta());
    registry.fill(HIST("h_part_jet_phi"), jet.phi());
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processMCParticleLevel, "jets on particle level MC", false);

  void processMCParticleLevelFull(soa::Filtered<aod::FullMCParticleLevelJets>::iterator const& jet)
  {
    registry.fill(HIST("h_part_jet_pt"), jet.pt());
    registry.fill(HIST("h_part_jet_eta"), jet.eta());
    registry.fill(HIST("h_part_jet_phi"), jet.phi());
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processMCParticleLevelFull, "full jets on particle level MC", false);

  void processMCCharged(aod::Collision const& collision, soa::Filtered<aod::ChargedMCDetectorLevelJets> const& MCDjets, soa::Filtered<aod::ChargedMCParticleLevelJets> const& MCPjets)
  {
    for (auto& MCDjet : MCDjets) {
      registry.fill(HIST("h_jet_pt"), MCDjet.pt());
      registry.fill(HIST("h_jet_eta"), MCDjet.eta());
      registry.fill(HIST("h_jet_phi"), MCDjet.phi());
    }
    for (auto& MCPjet : MCPjets) {
      registry.fill(HIST("h_part_jet_pt"), MCPjet.pt());
      registry.fill(HIST("h_part_jet_eta"), MCPjet.eta());
      registry.fill(HIST("h_part_jet_phi"), MCPjet.phi());
    }
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processMCCharged, "jets on detector and particle level MC", false);

  void processDataChargedSubstructure(soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>>::iterator const& jet, aod::Tracks const& tracks)
  {
    registry.fill(HIST("h_jet_pt"), jet.pt());
    registry.fill(HIST("h_jet_eta"), jet.eta());
    registry.fill(HIST("h_jet_phi"), jet.phi());
    registry.fill(HIST("h_jet_ntracks"), jet.tracks().size());
    double angularity = 0.0;
    for (auto& jetConstituent : jet.tracks_as<aod::Tracks>()) {
      angularity += jetConstituent.pt() * TMath::Sqrt(TMath::Power(jetConstituent.phi() - jet.phi(), 2.0) + TMath::Power(jetConstituent.eta() - jet.eta(), 2.0));
    }
    registry.fill(HIST("h_jet_angularity"), angularity / (jet.pt() * jet.r() / 100.0));
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processDataChargedSubstructure, "jet substructure charged jets", false);

  void processDataFullSubstructure(soa::Filtered<soa::Join<aod::FullJets, aod::FullJetConstituents>>::iterator const& jet, aod::Tracks const& tracks, aod::EMCALClusters const& clusters)
  {
    registry.fill(HIST("h_full_jet_pt"), jet.pt());
    registry.fill(HIST("h_full_jet_eta"), jet.eta());
    registry.fill(HIST("h_full_jet_phi"), jet.phi());
    registry.fill(HIST("h_full_jet_ntracks"), jet.tracks().size());
    registry.fill(HIST("h_full_jet_nclusters"), jet.clusters().size());
    double angularity = 0.0;
    for (auto& jetConstituent : jet.tracks_as<aod::Tracks>()) {
      angularity += jetConstituent.pt() * TMath::Sqrt(TMath::Power(jetConstituent.phi() - jet.phi(), 2.0) + TMath::Power(jetConstituent.eta() - jet.eta(), 2.0));
    }
    for (auto& jetConstituent : jet.clusters_as<aod::EMCALClusters>()) {
      angularity += jetConstituent.energy() * TMath::Sqrt(TMath::Power(jetConstituent.phi() - jet.phi(), 2.0) + TMath::Power(jetConstituent.eta() - jet.eta(), 2.0));
    }
    registry.fill(HIST("h_full_jet_angularity"), angularity / (jet.pt() * jet.r() / 100.0));
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processDataFullSubstructure, "jet substructure full jets", false);

  void processMCParticleSubstructure(soa::Filtered<soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents>>::iterator const& jet, aod::McParticles const& particles)
  {
    registry.fill(HIST("h_part_jet_pt"), jet.pt());
    registry.fill(HIST("h_part_jet_eta"), jet.eta());
    registry.fill(HIST("h_part_jet_phi"), jet.phi());
    registry.fill(HIST("h_part_jet_ntracks"), jet.tracks().size());
    double angularity = 0.0;
    for (auto& jetConstituent : jet.tracks_as<aod::McParticles>()) {
      angularity += jetConstituent.pt() * TMath::Sqrt(TMath::Power(jetConstituent.phi() - jet.phi(), 2.0) + TMath::Power(jetConstituent.eta() - jet.eta(), 2.0));
    }
    registry.fill(HIST("h_part_jet_angularity"), angularity / (jet.pt() * jet.r() / 100.0));
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processMCParticleSubstructure, "jet substructure particle level full jets", false);

  void processDataRecoil(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<aod::ChargedJets> const& jets, soa::Join<aod::Tracks, aod::TracksExtra> const& tracks)
  {
    if (!collision.sel8()) {
      return;
    }
    double leadingTrackpT = 0.0;
    double leadingTrackphi = 0.0;
    for (auto& track : tracks) {
      if (track.pt() > 4.0 && track.pt() < 10.0) {
        if (track.pt() > leadingTrackpT) {
          leadingTrackpT = track.pt();
          leadingTrackphi = track.phi();
        }
      }
    }
    if (leadingTrackpT == 0.0) {
      return;
    }
    for (auto& jet : jets) {

      if (jet.phi() - leadingTrackphi > 0.6) {

        registry.fill(HIST("h_recoil_jet_pt"), jet.pt());
        registry.fill(HIST("h_recoil_jet_eta"), jet.eta());
        registry.fill(HIST("h_recoil_jet_phi"), jet.phi());
        registry.fill(HIST("h_recoil_jet_dphi"), jet.phi() - leadingTrackphi);
      }
    }
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processDataRecoil, "hadron-recoil jets", false);

  using MCDJetTable = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>>;
  using MCPJetTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;
  void processMCMatched(aod::Collision const& collision, MCDJetTable const& MCDjets, MCPJetTable const& MCPjets)
  {
    /* for (const auto& MCDjet : MCDjets) {
      if (MCDjet.has_matchedJetGeo()) {
        const auto& MCPjet = MCDjet.matchedJetGeo_as<MCPJetTable>();
        registry.fill(HIST("h_matched_jets_pt"), MCPjet.pt(), MCDjet.pt());
        registry.fill(HIST("h_matched_jets_pt"), MCPjet.phi(), MCDjet.phi());
        registry.fill(HIST("h_matched_jets_pt"), MCPjet.eta(), MCDjet.eta());
      }
    }*/
  }
  PROCESS_SWITCH(JetTutorialSkeletonTask, processMCMatched, "jets matched on detector and particle level MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetTutorialSkeletonTask>(cfgc, TaskName{"jet-tutorial-skeleton"})}; }
