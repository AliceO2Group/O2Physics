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
///
/// \brief Simple QC task to run
/// \author Fabrizio Grosa, fabrizio.grosa@cern.ch (CERN)

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct mcSimpleQc {

  HistogramRegistry registry{
    "registry",
    {
      {"histRecoOverGenCollisions", ";#it{N}_{reco collisions} / #it{N}_{gen collisions};entries", {HistType::kTH1F, {{200, 0., 2.}}}}, //
      {"histDeltaCollisionPosX", ";#it{X}_{reco} - #it{X}_{gen};entries", {HistType::kTH1F, {{1000, -0.1, 0.1}}}},                      //
      {"histDeltaCollisionPosY", ";#it{Y}_{reco} - #it{Y}_{gen};entries", {HistType::kTH1F, {{1000, -0.1, 0.1}}}},                      //
      {"histDeltaCollisionPosZ", ";#it{Z}_{reco} - #it{Z}_{gen};entries", {HistType::kTH1F, {{1000, -1., 1.}}}},                        //
      {"histNumGlobalTracks", ";#it{N}_{global tracks};entries", {HistType::kTH1F, {{200, -0.5, 200.}}}},                               //
      {"histNumPvContributors", ";#it{N}_{PV contributors};entries", {HistType::kTH1F, {{200, -0.5, 200.}}}},                           //
      {"histFracAmbiguousTracks", ";#it{N}_{ambiguous tracks} / #it{N}_{tracks};entries", {HistType::kTH1F, {{100, 0., 1.}}}}           //
    }                                                                                                                                   //
  };

  using CollisionsWithMCLabels = soa::Join<aod::Collisions, aod::McCollisionLabels>;
  using TracksWithSel = soa::Join<Tracks, TracksExtra, TrackSelection>;
  Filter trackFilter = requireGlobalTrackWoDCAInFilter();

  Preslice<soa::Filtered<TracksWithSel>> perCollision = o2::aod::track::collisionId;

  // group tracks according to collision
  void process(CollisionsWithMCLabels const& collisions,
               soa::Filtered<TracksWithSel> const& tracks,
               AmbiguousTracks const& ambTracks,
               McCollisions const& mcCollisions,
               McParticles const&)
  {
    // checks on collisions
    registry.fill(HIST("histRecoOverGenCollisions"), float(collisions.size()) / mcCollisions.size());

    for (const auto& collision : collisions) {
      if (collision.has_mcCollision()) {
        const auto mcCollision = collision.mcCollision();
        registry.fill(HIST("histDeltaCollisionPosX"), collision.posX() - mcCollision.posX());
        registry.fill(HIST("histDeltaCollisionPosY"), collision.posY() - mcCollision.posY());
        registry.fill(HIST("histDeltaCollisionPosZ"), collision.posZ() - mcCollision.posZ());
      }
      registry.fill(HIST("histNumPvContributors"), collision.numContrib());

      auto tracksPerCollision = tracks.sliceBy(perCollision, collision.globalIndex());
      uint nTracks{0u}, nAmbTracks{0u};
      for (const auto& track : tracksPerCollision) {
        nTracks++;
        for (const auto& ambTrack : ambTracks) {
          if (ambTrack.trackId() == track.globalIndex()) {
            nAmbTracks++;
            break;
          }
        }
      }
      registry.fill(HIST("histNumGlobalTracks"), nTracks);
      registry.fill(HIST("histFracAmbiguousTracks"), float(nAmbTracks) / nTracks);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<mcSimpleQc>(cfgc)};
}