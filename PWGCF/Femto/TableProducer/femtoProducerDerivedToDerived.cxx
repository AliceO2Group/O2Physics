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

#include "PWGCF/Femto/Core/cascadeBuilder.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/kinkBuilder.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/v0Builder.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

using namespace o2::analysis::femto;

struct FemtoProducerDerivedToDerived {

  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using FemtoCollision = FemtoCollisions::iterator;

  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks>;
  using FemtoTracksWithMass = o2::soa::Join<FemtoTracks, o2::aod::FTrackMass>;
  using FemtoLambdas = o2::soa::Join<o2::aod::FLambdas, o2::aod::FLambdaMasks>;
  using FemtoK0shorts = o2::soa::Join<o2::aod::FK0shorts, o2::aod::FK0shortMasks>;
  using FemtoXis = o2::soa::Join<o2::aod::FXis, o2::aod::FXiMasks>;
  using FemtoOmegas = o2::soa::Join<o2::aod::FOmegas, o2::aod::FOmegaMasks>;
  using FemtoSigma = o2::soa::Join<o2::aod::FSigmas, o2::aod::FSigmaMasks>;
  using FemtoSigmaPlus = o2::soa::Join<o2::aod::FSigmaPlus, o2::aod::FSigmaPlusMasks>;

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

  o2::framework::Partition<FemtoTracks> trackPartition1 = MAKE_TRACK_PARTITION(trackSelections1);
  o2::framework::Partition<FemtoTracks> trackPartition2 = MAKE_TRACK_PARTITION(trackSelections2);
  o2::framework::Preslice<FemtoTracks> perColTracks = o2::aod::femtobase::stored::fColId;

  // v0 builder
  v0builder::V0BuilderDerivedToDerived v0Builder;
  v0builder::V0BuilderDerivedToDerivedProducts v0BuilderProducts;
  v0builder::ConfV0TablesDerivedToDerived confV0Builder;

  v0builder::ConfLambdaSelection1 lambdaSelection1;
  o2::framework::Partition<FemtoLambdas> lambdaPartition = MAKE_LAMBDA_PARTITION(lambdaSelection1);
  o2::framework::Preslice<FemtoLambdas> perColLambdas = o2::aod::femtobase::stored::fColId;

  v0builder::ConfK0shortSelection1 k0shortSelection1;
  o2::framework::Partition<FemtoK0shorts> k0shortPartition = MAKE_K0SHORT_PARTITION(k0shortSelection1);
  o2::framework::Preslice<FemtoK0shorts> perColK0shorts = o2::aod::femtobase::stored::fColId;

  // kink builder
  kinkbuilder::KinkBuilderDerivedToDerived kinkBuilder;
  kinkbuilder::KinkBuilderDerivedToDerivedProducts kinkBuilderProducts;
  kinkbuilder::ConfKinkTablesDerivedToDerived confKinkBuilder;

  kinkbuilder::ConfSigmaSelection1 sigmaSelection1;
  o2::framework::Partition<FemtoSigma> sigmaPartition = MAKE_SIGMA_PARTITION(sigmaSelection1);
  o2::framework::Preslice<FemtoSigma> perColSigma = o2::aod::femtobase::stored::fColId;

  kinkbuilder::ConfSigmaPlusSelection1 sigmaPlusSelection1;
  o2::framework::Partition<FemtoSigmaPlus> sigmaPlusPartition = MAKE_SIGMAPLUS_PARTITION(sigmaPlusSelection1);
  o2::framework::Preslice<FemtoSigmaPlus> perColSigmaPlus = o2::aod::femtobase::stored::fColId;

  // cascade builder
  cascadebuilder::CascadeBuilderDerivedToDerived cascadeBuilder;
  cascadebuilder::CascadeBuilderDerivedToDerivedProducts cascadeBuilderProducts;
  cascadebuilder::ConfCascadeTablesDerivedToDerived confCascadeBuilder;

  cascadebuilder::ConfXiSelection xiSelection;
  o2::framework::Partition<FemtoXis> xiPartition = MAKE_CASCADE_PARTITION(xiSelection);
  o2::framework::Preslice<FemtoXis> perColXis = o2::aod::femtobase::stored::fColId;

  cascadebuilder::ConfOmegaSelection omegaSelection;
  o2::framework::Partition<FemtoOmegas> omegaPartition = MAKE_CASCADE_PARTITION(omegaSelection);
  o2::framework::Preslice<FemtoOmegas> perColOmegas = o2::aod::femtobase::stored::fColId;

  void init(o2::framework::InitContext& /*context*/)
  {
    trackBuilder.init(confTrackBuilder);
    v0Builder.init(confV0Builder);
    cascadeBuilder.init(confCascadeBuilder);
    kinkBuilder.init(confKinkBuilder);

    if ((static_cast<int>(doprocessTracks) + static_cast<int>(doprocessTracksWithMass) + static_cast<int>(doprocessTracksLambdas) + static_cast<int>(doprocessLambdas) + static_cast<int>(doprocessTracksXis) + static_cast<int>(doprocessTracksOmegas) + static_cast<int>(doprocessTracksK0shorts) + static_cast<int>(doprocessTracksSigma) + static_cast<int>(doprocessTracksSigmaPlus)) > 1) {
      LOG(fatal) << "Only one proccess function can be activated";
    }
  }

  // proccess functions
  void processTracks(FilteredFemtoCollision const& col, FemtoTracks const& tracks)
  {
    if (trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache)) {
      return;
    }
    collisionBuilder.processCollision(col, collisionBuilderProducts);
    trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracks, "Process tracks", true);

  void processTracksWithMass(FilteredFemtoCollision const& col, FemtoTracksWithMass const& tracks)
  {
    if (trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache)) {
      return;
    }
    collisionBuilder.processCollision(col, collisionBuilderProducts);
    trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksWithMass, "Process tracks with mass", false);

  void processTracksLambdas(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoLambdas const& lambdas)
  {
    if (trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) || v0Builder.collisionHasTooFewLambdas(col, lambdas, lambdaPartition, cache)) {
      return;
    }
    collisionBuilder.processCollision(col, collisionBuilderProducts);
    trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    v0Builder.processLambdas(col, lambdas, tracks, lambdaPartition, trackBuilder, cache, v0BuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksLambdas, "Process lambdas and tracks", false);

  void processLambdas(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoLambdas const& lambdas)
  {
    if (v0Builder.collisionHasTooFewLambdas(col, lambdas, lambdaPartition, cache)) {
      return;
    }
    collisionBuilder.processCollision(col, collisionBuilderProducts);
    v0Builder.processLambdas(col, lambdas, tracks, lambdaPartition, trackBuilder, cache, v0BuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processLambdas, "Process lambdas", false);

  void processTracksXis(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoXis const& xis)
  {
    if (trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) || cascadeBuilder.collisionHasTooFewXis(col, xis, xiPartition, cache)) {
      return;
    }
    collisionBuilder.processCollision(col, collisionBuilderProducts);
    trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    cascadeBuilder.processXis(col, xis, tracks, xiPartition, trackBuilder, cache, cascadeBuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksXis, "Process tracks and xis", false);

  void processTracksOmegas(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoOmegas const& omegas)
  {
    if (trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) || cascadeBuilder.collisionHasTooFewOmegas(col, omegas, omegaPartition, cache)) {
      return;
    }
    collisionBuilder.processCollision(col, collisionBuilderProducts);
    trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    cascadeBuilder.processOmegas(col, omegas, tracks, omegaPartition, trackBuilder, cache, cascadeBuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksOmegas, "Process tracks and omegas", false);

  void processTracksK0shorts(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoK0shorts const& k0shorts)
  {
    if (trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) || v0Builder.collisionHasTooFewK0shorts(col, k0shorts, k0shortPartition, cache)) {
      return;
    }
    collisionBuilder.processCollision(col, collisionBuilderProducts);
    trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    v0Builder.processK0shorts(col, k0shorts, tracks, k0shortPartition, trackBuilder, cache, v0BuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksK0shorts, "Process k0short and tracks", false);

  void processTracksSigma(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoSigma const& sigma)
  {
    if (trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) || kinkBuilder.collisionHasTooFewSigma(col, sigma, sigmaPartition, cache)) {
      return;
    }
    collisionBuilder.processCollision(col, collisionBuilderProducts);
    trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    kinkBuilder.processSigma(col, sigma, tracks, sigmaPartition, trackBuilder, cache, kinkBuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksSigma, "Process sigma and tracks", false);

  void processTracksSigmaPlus(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoSigmaPlus const& sigmaplus)
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

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{adaptAnalysisTask<FemtoProducerDerivedToDerived>(context)};
  return workflow;
}
