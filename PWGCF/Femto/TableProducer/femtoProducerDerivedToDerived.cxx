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
#include "PWGCF/Femto/Core/mcBuilder.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/twoTrackResonanceBuilder.h"
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

// kept outside the task struct on purpose: the struct is reflected via a brace-decomposition
// mechanism capped at 100 direct data members (see o2::framework::PresliceGroup for the same
// concern applied to Preslice members below); type aliases don't need to live inside the struct
// at all, so moving them here keeps that budget free for actual data members.
namespace femtotables
{
using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;

using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

using FemtoMcCollisions = o2::soa::Join<FemtoCollisions, o2::aod::FColLabels>;
using FilteredFemtoMcCollisions = o2::soa::Filtered<FemtoMcCollisions>;

using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks>;
using FemtoTracksWithMass = o2::soa::Join<FemtoTracks, o2::aod::FTrackMass>;
// opt-in: carries the packed purity tables and FTrackDcas through to the re-written output --
// needs the upstream producer's TrackTables.produceItsPurity/produceTpcPurity/produceTofPurity/produceTrackDcas enabled
using FemtoTracksWithPurity = o2::soa::Join<FemtoTracks, o2::aod::FItsPurity, o2::aod::FTpcPurity, o2::aod::FTofPurity, o2::aod::FTrackDcas>;
using FemtoLambdas = o2::soa::Join<o2::aod::FLambdas, o2::aod::FLambdaMasks>;
using FemtoK0shorts = o2::soa::Join<o2::aod::FK0shorts, o2::aod::FK0shortMasks>;
using FemtoRho0s = o2::soa::Join<o2::aod::FRho0s, o2::aod::FRho0Masks>;
using FemtoPhis = o2::soa::Join<o2::aod::FPhis, o2::aod::FPhiMasks>;
using FemtoXis = o2::soa::Join<o2::aod::FXis, o2::aod::FXiMasks>;
using FemtoOmegas = o2::soa::Join<o2::aod::FOmegas, o2::aod::FOmegaMasks>;
using FemtoSigma = o2::soa::Join<o2::aod::FSigmas, o2::aod::FSigmaMasks>;
using FemtoSigmaPlus = o2::soa::Join<o2::aod::FSigmaPlus, o2::aod::FSigmaPlusMasks>;
using FemtoSigmaPlusExtended = o2::soa::Join<FemtoSigmaPlus, o2::aod::FSigmaPlusExtras>;
using FemtoSigmaPlusPcm = o2::soa::Join<o2::aod::FSigmaPlusPcm, o2::aod::FSigmaPlusPcmMasks, o2::aod::FSigmaPlusPcmMlScores>;
using FemtoSigmaPlusPcmExtended = o2::soa::Join<FemtoSigmaPlusPcm, o2::aod::FSigmaPlusPcmExtras>;

using FemtoMcTracks = o2::soa::Join<FemtoTracks, o2::aod::FTrackLabels>;
using FemtoMcSigma = o2::soa::Join<FemtoSigma, o2::aod::FSigmaLabels>;
using FemtoMcSigmaPlus = o2::soa::Join<FemtoSigmaPlus, o2::aod::FSigmaPlusLabels>;
using FemtoMcSigmaPlusExtended = o2::soa::Join<FemtoSigmaPlusExtended, o2::aod::FSigmaPlusLabels>;
using FemtoMcSigmaPlusPcm = o2::soa::Join<FemtoSigmaPlusPcm, o2::aod::FSigmaPlusPcmLabels>;
using FemtoMcSigmaPlusPcmExtended = o2::soa::Join<FemtoSigmaPlusPcmExtended, o2::aod::FSigmaPlusPcmLabels>;
using FemtoMcParticles = o2::soa::Join<o2::aod::FMcParticles, o2::aod::FMcMotherLabels>;
} // namespace femtotables

struct FemtoProducerDerivedToDerived {

  o2::framework::SliceCache cache;

  // all fColId Preslice declarations, grouped into a single data member
  struct : public o2::framework::PresliceGroup {
    o2::framework::Preslice<femtotables::FemtoTracks> perColTracks = o2::aod::femtobase::stored::fColId;
    o2::framework::Preslice<femtotables::FemtoTracksWithMass> perColTracksWithMass = o2::aod::femtobase::stored::fColId;
    o2::framework::Preslice<femtotables::FemtoTracksWithPurity> perColTracksWithPurity = o2::aod::femtobase::stored::fColId;
    o2::framework::Preslice<femtotables::FemtoLambdas> perColLambdas = o2::aod::femtobase::stored::fColId;
    o2::framework::Preslice<femtotables::FemtoK0shorts> perColK0shorts = o2::aod::femtobase::stored::fColId;
    o2::framework::Preslice<femtotables::FemtoRho0s> perColRho0s = o2::aod::femtobase::stored::fColId;
    o2::framework::Preslice<femtotables::FemtoPhis> perColPhis = o2::aod::femtobase::stored::fColId;
    o2::framework::Preslice<femtotables::FemtoSigma> perColSigma = o2::aod::femtobase::stored::fColId;
    o2::framework::Preslice<femtotables::FemtoSigmaPlus> perColSigmaPlus = o2::aod::femtobase::stored::fColId;
    o2::framework::Preslice<femtotables::FemtoSigmaPlusExtended> perColSigmaPlusExtended = o2::aod::femtobase::stored::fColId;
    o2::framework::Preslice<femtotables::FemtoSigmaPlusPcm> perColSigmaPlusPcm = o2::aod::femtobase::stored::fColId;
    o2::framework::Preslice<femtotables::FemtoSigmaPlusPcmExtended> perColSigmaPlusPcmExtended = o2::aod::femtobase::stored::fColId;
    o2::framework::Preslice<femtotables::FemtoXis> perColXis = o2::aod::femtobase::stored::fColId;
    o2::framework::Preslice<femtotables::FemtoOmegas> perColOmegas = o2::aod::femtobase::stored::fColId;
    o2::framework::Preslice<femtotables::FemtoMcTracks> perColMcTracks = o2::aod::femtobase::stored::fColId;
    o2::framework::Preslice<femtotables::FemtoMcSigma> perColMcSigma = o2::aod::femtobase::stored::fColId;
    o2::framework::Preslice<femtotables::FemtoMcSigmaPlus> perColMcSigmaPlus = o2::aod::femtobase::stored::fColId;
    o2::framework::Preslice<femtotables::FemtoMcSigmaPlusExtended> perColMcSigmaPlusExtended = o2::aod::femtobase::stored::fColId;
    o2::framework::Preslice<femtotables::FemtoMcSigmaPlusPcm> perColMcSigmaPlusPcm = o2::aod::femtobase::stored::fColId;
    o2::framework::Preslice<femtotables::FemtoMcSigmaPlusPcmExtended> perColMcSigmaPlusPcmExtended = o2::aod::femtobase::stored::fColId;
  } preslices;

  // plain builder helper objects, grouped into a single data member (not Configurable/Partition/
  // Preslice/Produces, so the framework never needs to discover them individually -- free to nest)
  struct {
    collisionbuilder::CollisionBuilderDerivedToDerived collisionBuilder;
    trackbuilder::TrackBuilderDerivedToDerived trackBuilder;
    v0builder::V0BuilderDerivedToDerived v0Builder;
    twotrackresonancebuilder::TwoTrackResonanceBuilderDerivedToDerived resonanceBuilder;
    kinkbuilder::KinkBuilderDerivedToDerived kinkBuilder;
    cascadebuilder::CascadeBuilderDerivedToDerived cascadeBuilder;
    mcbuilder::McBuilderDerivedToDerived mcBuilder;
  } builders;

  // collision builder
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  collisionbuilder::CollisionBuilderDerivedToDerivedProducts collisionBuilderProducts;

  // track builder
  trackbuilder::TrackBuilderDerivedToDerivedProducts trackBuilderProducts;
  trackbuilder::ConfTrackTablesDerivedToDerived confTrackBuilder;
  trackbuilder::ConfTrackSelection1 trackSelections1;
  trackbuilder::ConfTrackSelection2 trackSelections2;

  o2::framework::Partition<femtotables::FemtoTracks> trackPartition1 = MAKE_TRACK_PARTITION(trackSelections1);
  o2::framework::Partition<femtotables::FemtoTracks> trackPartition2 = MAKE_TRACK_PARTITION(trackSelections2);

  o2::framework::Partition<femtotables::FemtoTracksWithMass> trackWithMassPartition1 = MAKE_TRACK_PARTITION_WITH_MASS(trackSelections1);
  o2::framework::Partition<femtotables::FemtoTracksWithMass> trackWithMassPartition2 = MAKE_TRACK_PARTITION_WITH_MASS(trackSelections2);

  o2::framework::Partition<femtotables::FemtoTracksWithPurity> trackWithPurityPartition1 = MAKE_TRACK_PARTITION(trackSelections1);
  o2::framework::Partition<femtotables::FemtoTracksWithPurity> trackWithPurityPartition2 = MAKE_TRACK_PARTITION(trackSelections2);

  // v0 builder
  v0builder::V0BuilderDerivedToDerivedProducts v0BuilderProducts;
  v0builder::ConfV0TablesDerivedToDerived confV0Builder;

  v0builder::ConfLambdaSelection1 lambdaSelection1;
  o2::framework::Partition<femtotables::FemtoLambdas> lambdaPartition = MAKE_LAMBDA_PARTITION(lambdaSelection1);

  v0builder::ConfK0shortSelection1 k0shortSelection1;
  o2::framework::Partition<femtotables::FemtoK0shorts> k0shortPartition = MAKE_K0SHORT_PARTITION(k0shortSelection1);

  // two-track resonance builder
  twotrackresonancebuilder::TwoTrackResonanceBuilderDerivedToDerivedProducts resonanceBuilderProducts;
  twotrackresonancebuilder::ConfTwoTrackResonanceTablesDerivedToDerived confResonanceBuilder;

  twotrackresonancebuilder::ConfRho0Selection rho0Selection;
  o2::framework::Partition<femtotables::FemtoRho0s> rho0Partition = MAKE_RESONANCE_0_PARTITON(rho0Selection);

  twotrackresonancebuilder::ConfPhiSelection phiSelection;
  o2::framework::Partition<femtotables::FemtoPhis> phiPartition = MAKE_RESONANCE_0_PARTITON(phiSelection);

  // kink builder
  kinkbuilder::KinkBuilderDerivedToDerivedProducts kinkBuilderProducts;
  kinkbuilder::ConfKinkTablesDerivedToDerived confKinkBuilder;

  kinkbuilder::ConfSigmaSelection1 sigmaSelection1;
  o2::framework::Partition<femtotables::FemtoSigma> sigmaPartition = MAKE_SIGMA_PARTITION(sigmaSelection1);

  kinkbuilder::ConfSigmaPlusSelection1 sigmaPlusSelection1;
  o2::framework::Partition<femtotables::FemtoSigmaPlus> sigmaPlusPartition = MAKE_SIGMAPLUS_PARTITION(sigmaPlusSelection1);
  o2::framework::Partition<femtotables::FemtoSigmaPlusExtended> sigmaPlusExtendedPartition = MAKE_SIGMAPLUS_PARTITION(sigmaPlusSelection1);

  kinkbuilder::ConfSigmaPlusPcmSelection1 sigmaPlusPcmSelection1;
  o2::framework::Partition<femtotables::FemtoSigmaPlusPcm> sigmaPlusPcmPartition = MAKE_SIGMAPLUSPCM_PARTITION(sigmaPlusPcmSelection1);
  o2::framework::Partition<femtotables::FemtoSigmaPlusPcmExtended> sigmaPlusPcmExtendedPartition = MAKE_SIGMAPLUSPCM_PARTITION(sigmaPlusPcmSelection1);

  // cascade builder
  cascadebuilder::CascadeBuilderDerivedToDerivedProducts cascadeBuilderProducts;
  cascadebuilder::ConfCascadeTablesDerivedToDerived confCascadeBuilder;

  cascadebuilder::ConfXiSelection xiSelection;
  o2::framework::Partition<femtotables::FemtoXis> xiPartition = MAKE_CASCADE_PARTITION(xiSelection);

  cascadebuilder::ConfOmegaSelection omegaSelection;
  o2::framework::Partition<femtotables::FemtoOmegas> omegaPartition = MAKE_CASCADE_PARTITION(omegaSelection);

  // mc builder
  mcbuilder::McBuilderDerivedToDerivedProducts mcBuilderProducts;
  mcbuilder::ConfMcTablesDerivedToDerived confMcBuilder;

  o2::framework::Partition<femtotables::FemtoMcTracks> mcTrackPartition1 = MAKE_TRACK_PARTITION(trackSelections1);
  o2::framework::Partition<femtotables::FemtoMcTracks> mcTrackPartition2 = MAKE_TRACK_PARTITION(trackSelections2);

  o2::framework::Partition<femtotables::FemtoMcSigma> mcSigmaPartition = MAKE_SIGMA_PARTITION(sigmaSelection1);

  o2::framework::Partition<femtotables::FemtoMcSigmaPlus> mcSigmaPlusPartition = MAKE_SIGMAPLUS_PARTITION(sigmaPlusSelection1);
  o2::framework::Partition<femtotables::FemtoMcSigmaPlusExtended> mcSigmaPlusExtendedPartition = MAKE_SIGMAPLUS_PARTITION(sigmaPlusSelection1);

  o2::framework::Partition<femtotables::FemtoMcSigmaPlusPcm> mcSigmaPlusPcmPartition = MAKE_SIGMAPLUSPCM_PARTITION(sigmaPlusPcmSelection1);
  o2::framework::Partition<femtotables::FemtoMcSigmaPlusPcmExtended> mcSigmaPlusPcmExtendedPartition = MAKE_SIGMAPLUSPCM_PARTITION(sigmaPlusPcmSelection1);

  void init(o2::framework::InitContext& /*context*/)
  {
    const int activeProcesses =
      static_cast<int>(doprocessTracks) + static_cast<int>(doprocessTracksWithMass) + static_cast<int>(doprocessTracksWithPurity) +
      static_cast<int>(doprocessTracksLambdas) + static_cast<int>(doprocessTracksWithMassLambdas) +
      static_cast<int>(doprocessLambdas) +
      static_cast<int>(doprocessTracksXis) + static_cast<int>(doprocessTracksOmegas) +
      static_cast<int>(doprocessTracksK0shorts) + static_cast<int>(doprocessTracksSigma) +
      static_cast<int>(doprocessTracksSigmaPlus) + static_cast<int>(doprocessTracksSigmaPlusExtended) +
      static_cast<int>(doprocessTracksSigmaPlusPcm) + static_cast<int>(doprocessTracksSigmaPlusPcmExtended) +
      static_cast<int>(doprocessTracksRho0) + static_cast<int>(doprocessTracksPhi) +
      static_cast<int>(doprocessTracksMc) + static_cast<int>(doprocessTracksSigmaMc) +
      static_cast<int>(doprocessTracksSigmaPlusMc) + static_cast<int>(doprocessTracksSigmaPlusExtendedMc) +
      static_cast<int>(doprocessTracksSigmaPlusPcmMc) + static_cast<int>(doprocessTracksSigmaPlusPcmExtendedMc);
    if (activeProcesses != 1) {
      LOG(fatal) << "Exactly one process function must be activated (got " << activeProcesses << ").";
    }

    builders.trackBuilder.init(confTrackBuilder);
    builders.v0Builder.init(confV0Builder);
    builders.cascadeBuilder.init(confCascadeBuilder);
    builders.kinkBuilder.init(confKinkBuilder);
    builders.resonanceBuilder.init(confResonanceBuilder);
    builders.mcBuilder.init(confMcBuilder);
  }

  // process functions
  void processTracks(femtotables::FilteredFemtoCollision const& col, femtotables::FemtoTracks const& tracks)
  {
    if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache)) {
      return;
    }
    builders.trackBuilder.reset(tracks);
    builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
    builders.trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracks, "Process tracks", true);

  void processTracksWithMass(femtotables::FilteredFemtoCollision const& col, femtotables::FemtoTracksWithMass const& tracks)
  {
    if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, trackWithMassPartition1, trackWithMassPartition2, cache)) {
      return;
    }
    builders.trackBuilder.reset(tracks);
    builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
    builders.trackBuilder.processTracks(col, tracks, trackWithMassPartition1, trackWithMassPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksWithMass, "Process tracks with mass", false);

  void processTracksWithPurity(femtotables::FilteredFemtoCollision const& col, femtotables::FemtoTracksWithPurity const& tracks)
  {
    if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, trackWithPurityPartition1, trackWithPurityPartition2, cache)) {
      return;
    }
    builders.trackBuilder.reset(tracks);
    builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
    builders.trackBuilder.processTracks(col, tracks, trackWithPurityPartition1, trackWithPurityPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksWithPurity, "Process tracks, carrying through purity/fraction tables", false);

  void processTracksLambdas(femtotables::FilteredFemtoCollision const& col, femtotables::FemtoTracks const& tracks, femtotables::FemtoLambdas const& lambdas)
  {
    if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) || builders.v0Builder.collisionHasTooFewLambdas(col, lambdas, lambdaPartition, cache)) {
      return;
    }
    builders.trackBuilder.reset(tracks);
    builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
    builders.trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    builders.v0Builder.processLambdas(col, lambdas, tracks, lambdaPartition, builders.trackBuilder, cache, v0BuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksLambdas, "Process lambdas and tracks", false);

  void processTracksWithMassLambdas(femtotables::FilteredFemtoCollision const& col, femtotables::FemtoTracksWithMass const& tracks, femtotables::FemtoLambdas const& lambdas)
  {
    if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, trackWithMassPartition1, trackWithMassPartition2, cache) || builders.v0Builder.collisionHasTooFewLambdas(col, lambdas, lambdaPartition, cache)) {
      return;
    }
    builders.trackBuilder.reset(tracks);
    builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
    builders.trackBuilder.processTracks(col, tracks, trackWithMassPartition1, trackWithMassPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    builders.v0Builder.processLambdas(col, lambdas, tracks, lambdaPartition, builders.trackBuilder, cache, v0BuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksWithMassLambdas, "Process lambdas and tracks with mass", false);

  void processLambdas(femtotables::FilteredFemtoCollision const& col, femtotables::FemtoTracks const& tracks, femtotables::FemtoLambdas const& lambdas)
  {
    if (builders.v0Builder.collisionHasTooFewLambdas(col, lambdas, lambdaPartition, cache)) {
      return;
    }
    builders.trackBuilder.reset(tracks);
    builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
    builders.v0Builder.processLambdas(col, lambdas, tracks, lambdaPartition, builders.trackBuilder, cache, v0BuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processLambdas, "Process lambdas only (produced track table contains daughers only)", false);

  void processTracksXis(femtotables::FilteredFemtoCollision const& col, femtotables::FemtoTracks const& tracks, femtotables::FemtoXis const& xis)
  {
    if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) || builders.cascadeBuilder.collisionHasTooFewXis(col, xis, xiPartition, cache)) {
      return;
    }
    builders.trackBuilder.reset(tracks);
    builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
    builders.trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    builders.cascadeBuilder.processXis(col, xis, tracks, xiPartition, builders.trackBuilder, cache, cascadeBuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksXis, "Process tracks and xis", false);

  void processTracksOmegas(femtotables::FilteredFemtoCollision const& col, femtotables::FemtoTracks const& tracks, femtotables::FemtoOmegas const& omegas)
  {
    if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) || builders.cascadeBuilder.collisionHasTooFewOmegas(col, omegas, omegaPartition, cache)) {
      return;
    }
    builders.trackBuilder.reset(tracks);
    builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
    builders.trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    builders.cascadeBuilder.processOmegas(col, omegas, tracks, omegaPartition, builders.trackBuilder, cache, cascadeBuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksOmegas, "Process tracks and omegas", false);

  void processTracksK0shorts(femtotables::FilteredFemtoCollision const& col, femtotables::FemtoTracks const& tracks, femtotables::FemtoK0shorts const& k0shorts)
  {
    if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) || builders.v0Builder.collisionHasTooFewK0shorts(col, k0shorts, k0shortPartition, cache)) {
      return;
    }
    builders.trackBuilder.reset(tracks);
    builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
    builders.trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    builders.v0Builder.processK0shorts(col, k0shorts, tracks, k0shortPartition, builders.trackBuilder, cache, v0BuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksK0shorts, "Process k0short and tracks", false);

  void processTracksSigma(femtotables::FilteredFemtoCollision const& col, femtotables::FemtoTracks const& tracks, femtotables::FemtoSigma const& sigma)
  {
    if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) || builders.kinkBuilder.collisionHasTooFewSigma(col, sigma, sigmaPartition, cache)) {
      return;
    }
    builders.trackBuilder.reset(tracks);
    builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
    builders.trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    builders.kinkBuilder.processSigma(col, sigma, tracks, sigmaPartition, builders.trackBuilder, cache, kinkBuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksSigma, "Process sigma and tracks", false);

  void processTracksSigmaPlus(femtotables::FilteredFemtoCollision const& col, femtotables::FemtoTracks const& tracks, femtotables::FemtoSigmaPlus const& sigmaplus)
  {
    if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) || builders.kinkBuilder.collisionHasTooFewSigmaPlus(col, sigmaplus, sigmaPlusPartition, cache)) {
      return;
    }
    builders.trackBuilder.reset(tracks);
    builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
    builders.trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    builders.kinkBuilder.processSigmaPlus(col, sigmaplus, tracks, sigmaPlusPartition, builders.trackBuilder, cache, kinkBuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksSigmaPlus, "Process sigmaPlus and tracks", false);

  void processTracksSigmaPlusExtended(femtotables::FilteredFemtoCollision const& col, femtotables::FemtoTracks const& tracks, femtotables::FemtoSigmaPlusExtended const& sigmaplus)
  {
    if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) || builders.kinkBuilder.collisionHasTooFewSigmaPlus(col, sigmaplus, sigmaPlusExtendedPartition, cache)) {
      return;
    }
    builders.trackBuilder.reset(tracks);
    builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
    builders.trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    builders.kinkBuilder.processSigmaPlus(col, sigmaplus, tracks, sigmaPlusExtendedPartition, builders.trackBuilder, cache, kinkBuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksSigmaPlusExtended, "Process sigmaPlus (with SigmaPlusExtras) and tracks", false);

  void processTracksSigmaPlusPcm(femtotables::FilteredFemtoCollision const& col, femtotables::FemtoTracks const& tracks, femtotables::FemtoSigmaPlusPcm const& sigmapluspcm)
  {
    if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) || builders.kinkBuilder.collisionHasTooFewSigmaPlusPcm(col, sigmapluspcm, sigmaPlusPcmPartition, cache)) {
      return;
    }
    builders.trackBuilder.reset(tracks);
    builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
    builders.trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    builders.kinkBuilder.processSigmaPlusPcm(col, sigmapluspcm, tracks, sigmaPlusPcmPartition, builders.trackBuilder, cache, kinkBuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksSigmaPlusPcm, "Process sigmaPlusPcm (with SigmaPlusPcmMlScores, selected on the BDT score) and tracks", false);

  void processTracksSigmaPlusPcmExtended(femtotables::FilteredFemtoCollision const& col, femtotables::FemtoTracks const& tracks, femtotables::FemtoSigmaPlusPcmExtended const& sigmapluspcm)
  {
    if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) || builders.kinkBuilder.collisionHasTooFewSigmaPlusPcm(col, sigmapluspcm, sigmaPlusPcmExtendedPartition, cache)) {
      return;
    }
    builders.trackBuilder.reset(tracks);
    builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
    builders.trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    builders.kinkBuilder.processSigmaPlusPcm(col, sigmapluspcm, tracks, sigmaPlusPcmExtendedPartition, builders.trackBuilder, cache, kinkBuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksSigmaPlusPcmExtended, "Process sigmaPlusPcm (with SigmaPlusPcmMlScores and SigmaPlusPcmExtras, selected on the BDT score) and tracks", false);

  void processTracksRho0(femtotables::FilteredFemtoCollision const& col, femtotables::FemtoTracks const& tracks, femtotables::FemtoRho0s const& rho0s)
  {
    if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) || builders.resonanceBuilder.collisionHasTooFewRho0s(col, rho0s, rho0Partition, cache)) {
      return;
    }
    builders.trackBuilder.reset(tracks);
    builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
    builders.trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    builders.resonanceBuilder.processRho0s(col, rho0s, tracks, rho0Partition, builders.trackBuilder, cache, resonanceBuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksRho0, "Process rho0 and tracks", false);

  void processTracksPhi(femtotables::FilteredFemtoCollision const& col, femtotables::FemtoTracks const& tracks, femtotables::FemtoPhis const& phis)
  {
    if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, trackPartition1, trackPartition2, cache) || builders.resonanceBuilder.collisionHasTooFewPhis(col, phis, phiPartition, cache)) {
      return;
    }
    builders.trackBuilder.reset(tracks);
    builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
    builders.trackBuilder.processTracks(col, tracks, trackPartition1, trackPartition2, cache, trackBuilderProducts, collisionBuilderProducts);
    builders.resonanceBuilder.processPhis(col, phis, tracks, phiPartition, builders.trackBuilder, cache, resonanceBuilderProducts, trackBuilderProducts, collisionBuilderProducts);
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksPhi, "Process phi and tracks", false);

  // mc processing
  //   // NOTE: unlike the other process functions, this one subscribes to the full
  // collision table and iterates manually. The MC index maps are keyed by source
  // global index and must live for the whole dataframe
  void processTracksMc(femtotables::FilteredFemtoMcCollisions const& cols,
                       femtotables::FemtoMcTracks const& tracks,
                       o2::aod::FMcCols const& mcCols,
                       femtotables::FemtoMcParticles const& mcParticles,
                       o2::aod::FMcMothers const& mcMothers,
                       o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    // once per dataframe
    builders.mcBuilder.reset(mcCols, mcParticles);

    for (auto const& col : cols) {
      if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, mcTrackPartition1, mcTrackPartition2, cache)) {
        continue;
      }
      // track indexMap is per collision by construction, so it resets here
      builders.trackBuilder.reset(tracks);
      builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
      builders.mcBuilder.fillCollisionWithLabel(col, mcCols, mcBuilderProducts);
      builders.trackBuilder.processTracksWithMc(col, tracks, mcTrackPartition1, mcTrackPartition2, cache,
                                                trackBuilderProducts, collisionBuilderProducts,
                                                builders.mcBuilder, mcCols, mcParticles, mcMothers, mcPartonicMothers, mcBuilderProducts);
    }
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksMc, "Process tracks with MC information", false);

  void processTracksSigmaMc(femtotables::FilteredFemtoMcCollisions const& cols,
                            femtotables::FemtoMcTracks const& tracks,
                            femtotables::FemtoMcSigma const& sigma,
                            o2::aod::FMcCols const& mcCols,
                            femtotables::FemtoMcParticles const& mcParticles,
                            o2::aod::FMcMothers const& mcMothers,
                            o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    builders.mcBuilder.reset(mcCols, mcParticles);

    for (auto const& col : cols) {
      if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, mcTrackPartition1, mcTrackPartition2, cache) || builders.kinkBuilder.collisionHasTooFewSigma(col, sigma, mcSigmaPartition, cache)) {
        continue;
      }
      builders.trackBuilder.reset(tracks);
      builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
      builders.mcBuilder.fillCollisionWithLabel(col, mcCols, mcBuilderProducts);
      builders.trackBuilder.processTracksWithMc(col, tracks, mcTrackPartition1, mcTrackPartition2, cache,
                                                trackBuilderProducts, collisionBuilderProducts,
                                                builders.mcBuilder, mcCols, mcParticles, mcMothers, mcPartonicMothers, mcBuilderProducts);
      builders.kinkBuilder.processSigmaWithMc(col, sigma, tracks, mcSigmaPartition, builders.trackBuilder, cache, kinkBuilderProducts, trackBuilderProducts, collisionBuilderProducts,
                                              builders.mcBuilder, mcCols, mcParticles, mcMothers, mcPartonicMothers, mcBuilderProducts);
    }
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksSigmaMc, "Process sigma and tracks with MC information", false);

  void processTracksSigmaPlusMc(femtotables::FilteredFemtoMcCollisions const& cols,
                                femtotables::FemtoMcTracks const& tracks,
                                femtotables::FemtoMcSigmaPlus const& sigmaplus,
                                o2::aod::FMcCols const& mcCols,
                                femtotables::FemtoMcParticles const& mcParticles,
                                o2::aod::FMcMothers const& mcMothers,
                                o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    builders.mcBuilder.reset(mcCols, mcParticles);

    for (auto const& col : cols) {
      if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, mcTrackPartition1, mcTrackPartition2, cache) || builders.kinkBuilder.collisionHasTooFewSigmaPlus(col, sigmaplus, mcSigmaPlusPartition, cache)) {
        continue;
      }
      builders.trackBuilder.reset(tracks);
      builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
      builders.mcBuilder.fillCollisionWithLabel(col, mcCols, mcBuilderProducts);
      builders.trackBuilder.processTracksWithMc(col, tracks, mcTrackPartition1, mcTrackPartition2, cache,
                                                trackBuilderProducts, collisionBuilderProducts,
                                                builders.mcBuilder, mcCols, mcParticles, mcMothers, mcPartonicMothers, mcBuilderProducts);
      builders.kinkBuilder.processSigmaPlusWithMc(col, sigmaplus, tracks, mcSigmaPlusPartition, builders.trackBuilder, cache, kinkBuilderProducts, trackBuilderProducts, collisionBuilderProducts,
                                                  builders.mcBuilder, mcCols, mcParticles, mcMothers, mcPartonicMothers, mcBuilderProducts);
    }
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksSigmaPlusMc, "Process sigmaPlus and tracks with MC information", false);

  void processTracksSigmaPlusExtendedMc(femtotables::FilteredFemtoMcCollisions const& cols,
                                        femtotables::FemtoMcTracks const& tracks,
                                        femtotables::FemtoMcSigmaPlusExtended const& sigmaplus,
                                        o2::aod::FMcCols const& mcCols,
                                        femtotables::FemtoMcParticles const& mcParticles,
                                        o2::aod::FMcMothers const& mcMothers,
                                        o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    builders.mcBuilder.reset(mcCols, mcParticles);

    for (auto const& col : cols) {
      if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, mcTrackPartition1, mcTrackPartition2, cache) || builders.kinkBuilder.collisionHasTooFewSigmaPlus(col, sigmaplus, mcSigmaPlusExtendedPartition, cache)) {
        continue;
      }
      builders.trackBuilder.reset(tracks);
      builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
      builders.mcBuilder.fillCollisionWithLabel(col, mcCols, mcBuilderProducts);
      builders.trackBuilder.processTracksWithMc(col, tracks, mcTrackPartition1, mcTrackPartition2, cache,
                                                trackBuilderProducts, collisionBuilderProducts,
                                                builders.mcBuilder, mcCols, mcParticles, mcMothers, mcPartonicMothers, mcBuilderProducts);
      builders.kinkBuilder.processSigmaPlusWithMc(col, sigmaplus, tracks, mcSigmaPlusExtendedPartition, builders.trackBuilder, cache, kinkBuilderProducts, trackBuilderProducts, collisionBuilderProducts,
                                                  builders.mcBuilder, mcCols, mcParticles, mcMothers, mcPartonicMothers, mcBuilderProducts);
    }
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksSigmaPlusExtendedMc, "Process sigmaPlus (with SigmaPlusExtras) and tracks with MC information", false);

  void processTracksSigmaPlusPcmMc(femtotables::FilteredFemtoMcCollisions const& cols,
                                   femtotables::FemtoMcTracks const& tracks,
                                   femtotables::FemtoMcSigmaPlusPcm const& sigmapluspcm,
                                   o2::aod::FMcCols const& mcCols,
                                   femtotables::FemtoMcParticles const& mcParticles,
                                   o2::aod::FMcMothers const& mcMothers,
                                   o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    builders.mcBuilder.reset(mcCols, mcParticles);

    for (auto const& col : cols) {
      if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, mcTrackPartition1, mcTrackPartition2, cache) || builders.kinkBuilder.collisionHasTooFewSigmaPlusPcm(col, sigmapluspcm, mcSigmaPlusPcmPartition, cache)) {
        continue;
      }
      builders.trackBuilder.reset(tracks);
      builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
      builders.mcBuilder.fillCollisionWithLabel(col, mcCols, mcBuilderProducts);
      builders.trackBuilder.processTracksWithMc(col, tracks, mcTrackPartition1, mcTrackPartition2, cache,
                                                trackBuilderProducts, collisionBuilderProducts,
                                                builders.mcBuilder, mcCols, mcParticles, mcMothers, mcPartonicMothers, mcBuilderProducts);
      builders.kinkBuilder.processSigmaPlusPcmWithMc(col, sigmapluspcm, tracks, mcSigmaPlusPcmPartition, builders.trackBuilder, cache, kinkBuilderProducts, trackBuilderProducts, collisionBuilderProducts,
                                                     builders.mcBuilder, mcCols, mcParticles, mcMothers, mcPartonicMothers, mcBuilderProducts);
    }
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksSigmaPlusPcmMc, "Process sigmaPlusPcm (with SigmaPlusPcmMlScores, selected on the BDT score) and tracks with MC information", false);

  void processTracksSigmaPlusPcmExtendedMc(femtotables::FilteredFemtoMcCollisions const& cols,
                                           femtotables::FemtoMcTracks const& tracks,
                                           femtotables::FemtoMcSigmaPlusPcmExtended const& sigmapluspcm,
                                           o2::aod::FMcCols const& mcCols,
                                           femtotables::FemtoMcParticles const& mcParticles,
                                           o2::aod::FMcMothers const& mcMothers,
                                           o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    builders.mcBuilder.reset(mcCols, mcParticles);

    for (auto const& col : cols) {
      if (builders.trackBuilder.collisionHasTooFewTracks(col, tracks, mcTrackPartition1, mcTrackPartition2, cache) || builders.kinkBuilder.collisionHasTooFewSigmaPlusPcm(col, sigmapluspcm, mcSigmaPlusPcmExtendedPartition, cache)) {
        continue;
      }
      builders.trackBuilder.reset(tracks);
      builders.collisionBuilder.processCollision(col, collisionBuilderProducts);
      builders.mcBuilder.fillCollisionWithLabel(col, mcCols, mcBuilderProducts);
      builders.trackBuilder.processTracksWithMc(col, tracks, mcTrackPartition1, mcTrackPartition2, cache,
                                                trackBuilderProducts, collisionBuilderProducts,
                                                builders.mcBuilder, mcCols, mcParticles, mcMothers, mcPartonicMothers, mcBuilderProducts);
      builders.kinkBuilder.processSigmaPlusPcmWithMc(col, sigmapluspcm, tracks, mcSigmaPlusPcmExtendedPartition, builders.trackBuilder, cache, kinkBuilderProducts, trackBuilderProducts, collisionBuilderProducts,
                                                     builders.mcBuilder, mcCols, mcParticles, mcMothers, mcPartonicMothers, mcBuilderProducts);
    }
  }
  PROCESS_SWITCH(FemtoProducerDerivedToDerived, processTracksSigmaPlusPcmExtendedMc, "Process sigmaPlusPcm (with SigmaPlusPcmMlScores and SigmaPlusPcmExtras, selected on the BDT score) and tracks with MC information", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{adaptAnalysisTask<FemtoProducerDerivedToDerived>(context)};
  return workflow;
}
