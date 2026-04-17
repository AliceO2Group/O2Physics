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

/// \file femtoProducerResonances.cxx
/// \brief Tasks that produces femto tables for all resonances
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/twoTrackResonanceBuilder.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

using namespace o2::analysis::femto;

struct FemtoProducerResonances {

  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks>;

  o2::framework::SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // setup for resonance daughter tracks
  trackbuilder::ConfPionPlusSelection confPionPlusSelection;
  o2::framework::Partition<FemtoTracks> pionPlusPartition = MAKE_TRACK_PARTITION(confPionPlusSelection);

  trackbuilder::ConfPionMinusSelection confPionMinusSelection;
  o2::framework::Partition<FemtoTracks> pionMinusPartition = MAKE_TRACK_PARTITION(confPionMinusSelection);

  trackbuilder::ConfKaonPlusSelection confKaonPlusSelection;
  o2::framework::Partition<FemtoTracks> kaonPlusPartition = MAKE_TRACK_PARTITION(confKaonPlusSelection);

  trackbuilder::ConfKaonMinusSelection confKaonMinusSelection;
  o2::framework::Partition<FemtoTracks> kaonMinusPartition = MAKE_TRACK_PARTITION(confKaonMinusSelection);

  o2::framework::Preslice<FemtoTracks> perColtracks = o2::aod::femtobase::stored::fColId;

  // resonance filters
  twotrackresonancebuilder::ConfRhoFilters confRhoFilter;
  twotrackresonancebuilder::ConfPhiFilters confPhiFilter;
  twotrackresonancebuilder::ConfKstarFilters confKstarFilter;

  // resonance builders
  twotrackresonancebuilder::ConfTwoTrackResonanceTables confTwoTrackResonanceTables;
  twotrackresonancebuilder::TwoTrackResonanceBuilderProducts twoTrackResonanceBuilderProducts;
  twotrackresonancebuilder::TwoTrackResonanceBuilder<modes::TwoTrackResonance::kRho0> rho0Builder;
  twotrackresonancebuilder::TwoTrackResonanceBuilder<modes::TwoTrackResonance::kPhi> phiBuilder;
  twotrackresonancebuilder::TwoTrackResonanceBuilder<modes::TwoTrackResonance::kKstar0> kstar0Builder;
  twotrackresonancebuilder::TwoTrackResonanceBuilder<modes::TwoTrackResonance::kKstar0Bar> kstar0BarBuilder;

  void init(o2::framework::InitContext& context)
  {
    // init builders
    rho0Builder.init(confRhoFilter, confPionPlusSelection, confPionMinusSelection, confTwoTrackResonanceTables, context);
    phiBuilder.init(confPhiFilter, confKaonPlusSelection, confKaonMinusSelection, confTwoTrackResonanceTables, context);
    kstar0Builder.init(confKstarFilter, confKaonPlusSelection, confPionMinusSelection, confTwoTrackResonanceTables, context);
    kstar0BarBuilder.init(confKstarFilter, confPionPlusSelection, confKaonMinusSelection, confTwoTrackResonanceTables, context);
  }

  // proccess functions
  void processRho0(FilteredFemtoCollision const& col, FemtoTracks const& tracks)
  {
    rho0Builder.fillResonances(col, twoTrackResonanceBuilderProducts, pionPlusPartition, pionMinusPartition, tracks, cache);
  }
  PROCESS_SWITCH(FemtoProducerResonances, processRho0, "Build Rho0 candidates", true);

  void processPhi(FilteredFemtoCollision const& col, FemtoTracks const& tracks)
  {
    phiBuilder.fillResonances(col, twoTrackResonanceBuilderProducts, kaonPlusPartition, kaonMinusPartition, tracks, cache);
  }
  PROCESS_SWITCH(FemtoProducerResonances, processPhi, "Build Phi candidates", true);

  void processKstar0(FilteredFemtoCollision const& col, FemtoTracks const& tracks)
  {
    kstar0Builder.fillResonances(col, twoTrackResonanceBuilderProducts, kaonPlusPartition, pionMinusPartition, tracks, cache);
    kstar0BarBuilder.fillResonances(col, twoTrackResonanceBuilderProducts, pionPlusPartition, kaonMinusPartition, tracks, cache);
  }
  PROCESS_SWITCH(FemtoProducerResonances, processKstar0, "Build Kstar0/Kstar0bar candidates", true);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::framework::WorkflowSpec workflow{adaptAnalysisTask<FemtoProducerResonances>(cfgc)};
  return workflow;
}
