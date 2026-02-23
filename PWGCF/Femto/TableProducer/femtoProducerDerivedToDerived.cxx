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
#include "PWGCF/Femto/Core/kinkBuilder.h"
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

using namespace o2::analysis::femto;

struct FemtoProducerDerivedToDerived {

  // setup tables
  using Collisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  using Tracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks>;
  using Lambdas = o2::soa::Join<o2::aod::FLambdas, o2::aod::FLambdaMasks>;
  using K0shorts = o2::soa::Join<o2::aod::FK0shorts, o2::aod::FK0shortMasks>;
  using Sigma = o2::soa::Join<o2::aod::FSigmas, o2::aod::FSigmaMasks>;
  using SigmaPlus = o2::soa::Join<o2::aod::FSigmaPlus, o2::aod::FSigmaPlusMasks>;

  o2::framework::SliceCache cache;

  // collision builder
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  collisionbuilder::CollisionBuilderDerivedToDerivedProducts collisionBuilderProducts;
  collisionbuilder::CollisionBuilderDerivedToDerived collisionBuilder;

  // track builder
  trackbuilder::TrackBuilderDerivedToDerived trackBuilder;
  trackbuilder::TrackBuilderDerivedToDerivedProducts trackBuilderProducts;
  trackbuilder::ConfTrackTablesDerivedToDerived confTrackBuilder;
  trackbuilder::ConfTrackSelection1 trackSelections1;
  trackbuilder::ConfTrackSelection2 trackSelections2;

  o2::framework::Partition<Tracks> trackPartition1 = MAKE_TRACK_PARTITION(trackSelections1);
  o2::framework::Partition<Tracks> trackPartition2 = MAKE_TRACK_PARTITION(trackSelections2);
  o2::framework::Preslice<Tracks> perColTracks = o2::aod::femtobase::stored::fColId;

  // v0 builder
  v0builder::V0BuilderDerivedToDerived v0Builder;
  v0builder::V0BuilderDerivedToDerivedProducts v0BuilderProducts;
  v0builder::ConfV0TablesDerivedToDerived confV0Builder;

  v0builder::ConfLambdaSelection1 lambdaSelection1;
  o2::framework::Partition<Lambdas> lambdaPartition = MAKE_LAMBDA_PARTITION(lambdaSelection1);
  o2::framework::Preslice<Lambdas> perColLambdas = o2::aod::femtobase::stored::fColId;

  v0builder::ConfK0shortSelection1 k0shortSelection1;
  o2::framework::Partition<K0shorts> k0shortPartition = MAKE_K0SHORT_PARTITION(k0shortSelection1);
  o2::framework::Preslice<K0shorts> perColK0shorts = o2::aod::femtobase::stored::fColId;

  // kink builder
  kinkbuilder::KinkBuilderDerivedToDerived kinkBuilder;
  kinkbuilder::KinkBuilderDerivedToDerivedProducts kinkBuilderProducts;
  kinkbuilder::ConfKinkTablesDerivedToDerived confKinkBuilder;

  kinkbuilder::ConfSigmaSelection1 sigmaSelection1;
  o2::framework::Partition<Sigma> sigmaPartition = MAKE_SIGMA_PARTITION(sigmaSelection1);
  o2::framework::Preslice<Sigma> perColSigma = o2::aod::femtobase::stored::fColId;

  kinkbuilder::ConfSigmaPlusSelection1 sigmaPlusSelection1;
  o2::framework::Partition<SigmaPlus> sigmaPlusPartition = MAKE_SIGMAPLUS_PARTITION(sigmaPlusSelection1);
  o2::framework::Preslice<SigmaPlus> perColSigmaPlus = o2::aod::femtobase::stored::fColId;

  void init(o2::framework::InitContext& /*context*/)
  {
    trackBuilder.init(confTrackBuilder);
    v0Builder.init(confV0Builder);
    kinkBuilder.init(confKinkBuilder);

    if ((doprocessTracks + doprocessTracksLambdas + doprocessTracksK0shorts + doprocessTracksSigma + doprocessTracksSigmaPlus) > 1) {
      LOG(fatal) << "Only one proccess function can be activated";
    }
  }

  // proccess functions
  void processTracks(FilteredCollision const& col, Tracks const& tracks)
  {
    if (trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache)) {
      return;
    }
    collisionBuilder.processCollision(col, collisionBuilderProducts);
    trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracks, "Process tracks", true);

  void processTracksLambdas(FilteredCollision const& col, Tracks const& tracks, Lambdas const& lambdas)
  {
    if (trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) || v0Builder.collisionHasTooFewLambdas(col, lambdas, lambdaPartition, cache)) {
      return;
    }
    collisionBuilder.processCollision(col, collisionBuilderProducts);
    trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    v0Builder.processLambdas(col, lambdas, tracks, lambdaPartition, trackBuilder, cache, v0BuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksLambdas, "Process lambdas and tracks", false);

  void processTracksK0shorts(FilteredCollision const& col, Tracks const& tracks, K0shorts const& k0shorts)
  {
    if (trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) || v0Builder.collisionHasTooFewK0shorts(col, k0shorts, k0shortPartition, cache)) {
      return;
    }
    collisionBuilder.processCollision(col, collisionBuilderProducts);
    trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    v0Builder.processK0shorts(col, k0shorts, tracks, k0shortPartition, trackBuilder, cache, v0BuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksK0shorts, "Process k0short and tracks", false);

  void processTracksSigma(FilteredCollision const& col, Tracks const& tracks, Sigma const& sigma)
  {
    if (trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) || kinkBuilder.collisionHasTooFewSigma(col, sigma, sigmaPartition, cache)) {
      return;
    }
    collisionBuilder.processCollision(col, collisionBuilderProducts);
    trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    kinkBuilder.processSigma(col, sigma, tracks, sigmaPartition, trackBuilder, cache, kinkBuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksSigma, "Process sigma and tracks", false);

  void processTracksSigmaPlus(FilteredCollision const& col, Tracks const& tracks, SigmaPlus const& sigmaplus)
  {
    if (trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) || kinkBuilder.collisionHasTooFewSigmaPlus(col, sigmaplus, sigmaPlusPartition, cache)) {
      return;
    }
    collisionBuilder.processCollision(col, collisionBuilderProducts);
    trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    kinkBuilder.processSigmaPlus(col, sigmaplus, tracks, sigmaPlusPartition, trackBuilder, cache, kinkBuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksSigmaPlus, "Process sigmaPlus and tracks", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::framework::WorkflowSpec workflow{adaptAnalysisTask<FemtoProducerDerivedToDerived>(cfgc)};
  return workflow;
}
