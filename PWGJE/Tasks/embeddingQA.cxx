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

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/StringHelpers.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include <Framework/runDataProcessing.h>

struct EmbeddingQATask {
  HistogramRegistry registry{"registry",
                             {{"h_collision_zdiff", "collision z diff;z_{Data}-z_{MC} (cm);entries", {HistType::kTH1F, {{200, -10., 10.}}}},
                              {"h_tracks_all_pt", "all tracks;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_tracks_embedded_pt", "embedded tracks;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_tracks_sub_pt", "subtracted tracks;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_particles_pt", "particles;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jets_pt", "jets;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jets_sub_pt", "subtracted jets;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jets_mcd_pt", "detector-level jets;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jets_mcp_pt", "particle-level jets;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jet_detector_response", "detector response matrix;#it{p}_{T,jet} particle (GeV/#it{c});#it{p}_{T,jet} detector (GeV/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}}},
                              {"h_jet_response", "response matrix;#it{p}_{T,jet} particle (GeV/#it{c});#it{p}_{T,jet} embedded subtracted (GeV/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}}},
                              {"h_jet_background_response", "background response matrix;#it{p}_{T,jet} detector (GeV/#it{c});#it{p}_{T,jet} embedded (GeV/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}}},
                              {"h_jet_backgroundfluctuations_response", "background fluctuations response matrix;#it{p}_{T,jet} embedded subtracted (GeV/#it{c});#it{p}_{T,jet} embedded (GeV/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}}},
                              {"h_tracks_data_pt", "tracks derived data;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_tracks_mcd_pt", "tracks derived detector level;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_tracks_mcp_pt", "particles derived particle level;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_tracks_data_original_pt", "tracks original data;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_tracks_mcd_original_pt", "tracks original detector level;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_tracks_mcp_original_pt", "particles original particle level;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}}}};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  std::vector<int> eventSelection;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
  }

  Preslice<aod::Collisions> perBCCollisions = aod::collision::bcId;

  void processTracksOriginal(aod::Tracks const& tracks, aod::TracksFrom<o2::aod::Hash<"EMB"_h>> const& tracksMCD, aod::McParticlesFrom<o2::aod::Hash<"EMB"_h>> const& tracksMCP, aod::JetTracks const& tracksEmbedded, aod::JetParticles const& particles, aod::Collisions const& collisions, soa::Join<aod::CollisionsFrom<o2::aod::Hash<"EMB"_h>>, aod::McCollisionLabelsFrom<o2::aod::Hash<"EMB"_h>>> const& mcdCollisions, aod::McCollisionsFrom<o2::aod::Hash<"EMB"_h>> const&, aod::BCs const& bcs, aod::BCsFrom<o2::aod::Hash<"EMB"_h>> const&)
  {
    std::vector<bool> mcdCollisionGoodForEmbedding(mcdCollisions.size(), false);
    for (auto const& mcdCollision : mcdCollisions) {
      auto const& mcdCollisionBC = mcdCollision.bc_as<aod::BCsFrom<o2::aod::Hash<"EMB"_h>>>();
      auto const& mcdGlobalBC = mcdCollisionBC.globalBC();
      for (auto const& bc : bcs) {
        if (bc.globalBC() == mcdGlobalBC) {
          auto const& collisionsPerBC = collisions.sliceBy(perBCCollisions, bc.globalIndex());
          auto const& mcCollision = mcdCollision.mcCollision_as<aod::McCollisionsFrom<o2::aod::Hash<"EMB"_h>>>();
          if (collisionsPerBC.size() > 0) {
            mcdCollisionGoodForEmbedding[mcdCollision.globalIndex()] = true;
            registry.fill(HIST("h_collision_zdiff"), collisionsPerBC.iteratorAt(0).posZ() - mcCollision.posZ());
            break;
          }
        }
      }
    }

    for (auto const& track : tracks) {
      if (track.collisionId() >= 0) {
        registry.fill(HIST("h_tracks_data_original_pt"), track.pt());
      }
    }
    for (auto const& track : tracksMCD) {
      if (track.collisionId() >= 0) {
        if (mcdCollisionGoodForEmbedding[track.collisionId()]) {
          registry.fill(HIST("h_tracks_mcd_original_pt"), track.pt());
        }
      }
    }
    for (auto const& track : tracksMCP) {
      registry.fill(HIST("h_tracks_mcp_original_pt"), track.pt());
    }
    for (auto const& track : tracksEmbedded) {
      if (!(track.trackSel() & (1ULL << jetderiveddatautilities::JTrackSel::embeddedTrack))) {
        registry.fill(HIST("h_tracks_data_pt"), track.pt());
      }
    }
    for (auto const& track : tracksEmbedded) {
      if ((track.trackSel() & (1ULL << jetderiveddatautilities::JTrackSel::embeddedTrack))) {
        registry.fill(HIST("h_tracks_mcd_pt"), track.pt());
      }
    }
    for (auto const& particle : particles) {
      registry.fill(HIST("h_tracks_mcp_pt"), particle.pt());
    }
  }
  PROCESS_SWITCH(EmbeddingQATask, processTracksOriginal, "track QA", true);

  void processTracks(aod::JetCollision const&, aod::JetTracks const& tracks, aod::JetTracksSub const& subTracks)
  {
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        registry.fill(HIST("h_tracks_all_pt"), track.pt());
      }
      if (jetderiveddatautilities::selectTrack(track, trackSelection, true)) {
        registry.fill(HIST("h_tracks_embedded_pt"), track.pt());
      }
    }
    for (auto const& track : subTracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        registry.fill(HIST("h_tracks_sub_pt"), track.pt());
      }
    }
  }
  PROCESS_SWITCH(EmbeddingQATask, processTracks, "track QA", true);

  void processParticles(aod::JetMcCollision const&, aod::JetParticles const& particles)
  {
    for (auto const& particle : particles) {
      registry.fill(HIST("h_particles_pt"), particle.pt());
    }
  }
  PROCESS_SWITCH(EmbeddingQATask, processParticles, "track QA", true);

  void processJets(aod::ChargedJets const& jets, aod::ChargedEventWiseSubtractedJets const& subJets, aod::ChargedMCDetectorLevelJets const& mcdJets, aod::ChargedMCParticleLevelJets const& mcpJets)
  {
    for (auto const& jet : jets) {
      registry.fill(HIST("h_jets_pt"), jet.pt());
    }
    for (auto const& jet : subJets) {
      registry.fill(HIST("h_jets_sub_pt"), jet.pt());
    }
    for (auto const& jet : mcdJets) {
      registry.fill(HIST("h_jets_mcd_pt"), jet.pt());
    }
    for (auto const& jet : mcpJets) {
      registry.fill(HIST("h_jets_mcp_pt"), jet.pt());
    }
  }
  PROCESS_SWITCH(EmbeddingQATask, processJets, "jet QA", true);

  void processResponse(aod::ChargedJets const&, soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetsMatchedToChargedJets> const&, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetsMatchedToChargedEventWiseSubtractedJets> const&, soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets> const& mcpJets)
  {

    for (auto const& mcpJet : mcpJets) {
      for (auto const& mcdJet : mcpJet.matchedJetGeo_as<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetsMatchedToChargedEventWiseSubtractedJets>>()) {
        registry.fill(HIST("h_jet_detector_response"), mcpJet.pt(), mcdJet.pt());
        for (auto const& subJet : mcdJet.matchedJetGeo_as<soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetsMatchedToChargedJets>>()) {
          registry.fill(HIST("h_jet_response"), mcpJet.pt(), subJet.pt());
          for (auto const& jet : subJet.matchedJetGeo_as<aod::ChargedJets>()) {
            registry.fill(HIST("h_jet_background_response"), mcdJet.pt(), jet.pt());
            registry.fill(HIST("h_jet_backgroundfluctuations_response"), subJet.pt(), jet.pt());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(EmbeddingQATask, processResponse, "jet response", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<EmbeddingQATask>(cfgc, TaskName{"embedding-qa"})}; }
