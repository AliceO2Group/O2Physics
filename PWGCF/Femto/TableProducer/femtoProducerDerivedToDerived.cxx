// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoProducerDerivedToDerived.cxx
/// \brief Tasks that produces the femto tables from derived data
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/v0Builder.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/Expressions.h"
#include "Framework/InitContext.h"
#include "Framework/runDataProcessing.h"

#include <cstdint>
#include <string>
#include <unordered_map>

using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femto;

struct FemtoProducerDerivedToDerived {

  // setup tables
  using Collisions = Join<FCols, FColMasks>;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  using Tracks = o2::soa::Join<FTracks, FTrackMasks>;
  using Lambdas = o2::soa::Join<FLambdas, FLambdaMasks>;
  using K0shorts = o2::soa::Join<FK0shorts, FK0shortMasks>;

  SliceCache cache;

  // collision builder
  collisionbuilder::ConfCollisionSelection collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  collisionbuilder::CollisionBuilderDerivedToDerivedProducts collisionBuilderProducts;
  collisionbuilder::CollisionBuilderDerivedToDerived collisionBuilder;

  // track builder
  trackbuilder::TrackBuilderDerivedToDerived trackBuilder;
  trackbuilder::TrackBuilderDerivedToDerivedProducts trackBuilderProducts;
  trackbuilder::ConfTrackTablesDerivedToDerived confTrackBuilder;
  trackbuilder::ConfTrackSelection1 trackSelections1;
  trackbuilder::ConfTrackSelection2 trackSelections2;

  Partition<Tracks> trackPartition1 = MAKE_TRACK_PARTITION(trackSelections1);
  Partition<Tracks> trackPartition2 = MAKE_TRACK_PARTITION(trackSelections2);
  Preslice<Tracks> perColTracks = femtobase::stored::fColId;

  // v0 builder
  v0builder::V0BuilderDerivedToDerived v0Builder;
  v0builder::V0BuilderDerivedToDerivedProducts v0BuilderProducts;
  v0builder::ConfV0TablesDerivedToDerived confV0Builder;

  v0builder::ConfLambdaSelection1 lambdaSelection1;
  Partition<Lambdas> lambdaPartition = MAKE_LAMBDA_PARTITION(lambdaSelection1);
  Preslice<Lambdas> perColLambdas = femtobase::stored::fColId;

  v0builder::ConfK0shortSelection1 k0shortSelection1;
  Partition<K0shorts> k0shortPartition = MAKE_K0SHORT_PARTITION(k0shortSelection1);
  Preslice<K0shorts> perColK0shorts = femtobase::stored::fColId;

  std::unordered_map<int64_t, int64_t>
    indexMapTracks; // for mapping tracks to lambdas, cascades and resonances

  void init(InitContext& /*context*/)
  {
    trackBuilder.init(confTrackBuilder);
    v0Builder.init(confV0Builder);

    if ((doprocessTracks + doprocessLambdas + doprocessK0shorts) > 1) {
      LOG(fatal) << "Only one proccess function can be activated";
    }
  }

  // proccess functions
  void processTracks(FilteredCollision const& col, Tracks const& tracks)
  {
    if (trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache)) {
      return;
    }
    indexMapTracks.clear();
    collisionBuilder.processCollision(col, collisionBuilderProducts);
    trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, indexMapTracks, cache, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracks, "Process tracks", true);

  void processLambdas(FilteredCollision const& col, Tracks const& tracks, Lambdas const& lambdas)
  {
    if (trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) && v0Builder.collisionHasTooFewLambdas(col, lambdas, lambdaPartition, cache)) {
      return;
    }
    indexMapTracks.clear();
    if (trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache)) {
      collisionBuilder.processCollision(col, collisionBuilderProducts);
      trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, indexMapTracks, cache, trackBuilderProducts, collisionBuilderProducts);
      v0Builder.processLambdas(col, lambdas, tracks, lambdaPartition, trackBuilder, indexMapTracks, cache, v0BuilderProducts, trackBuilderProducts, collisionBuilderProducts);
    }
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processLambdas, "Process lambdas and tracks", false);

  void processK0shorts(FilteredCollision const& col, Tracks const& tracks, K0shorts const& k0shorts)
  {
    if (trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) && v0Builder.collisionHasTooFewK0shorts(col, k0shorts, k0shortPartition, cache)) {
      return;
    }
    indexMapTracks.clear();
    if (trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache)) {
      collisionBuilder.processCollision(col, collisionBuilderProducts);
      trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, indexMapTracks, cache, trackBuilderProducts, collisionBuilderProducts);
      v0Builder.processK0shorts(col, k0shorts, tracks, k0shortPartition, trackBuilder, indexMapTracks, cache, v0BuilderProducts, trackBuilderProducts, collisionBuilderProducts);
    }
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processK0shorts, "Process k0short and tracks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FemtoProducerDerivedToDerived>(cfgc)};
  return workflow;
}
