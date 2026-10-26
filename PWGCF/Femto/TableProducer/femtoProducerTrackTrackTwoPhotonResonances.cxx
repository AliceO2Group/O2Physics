// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoProducerTrackTrackTwoPhotonResonances.cxx
/// \brief Task that produces femto tables for resonances built from two tracks and a two-photon
///        resonance (e.g. omega(782) -> pi+ pi- pi0)
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/trackTrackTwoPhotonResonanceBuilder.h"
#include "PWGCF/Femto/Core/twoPhotonResonanceBuilder.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

using namespace o2::analysis::femto;

// reads the already-built FTracks (for the charged pion daughters) and FPi0s (for the neutral
// daughter) and combines them into omega(782) -> pi+ pi- pi0 candidates. A 3rd-stage producer,
// chained after femtoProducer.cxx (FTracks) and femtoProducerTwoPhotonResonances.cxx (FPi0s), same
// pattern as femtoProducerTwoPhotonResonances.cxx itself reading femtoProducer.cxx's FPhotons.
struct FemtoProducerTrackTrackTwoPhotonResonances {

  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks>;
  using FemtoPi0s = o2::soa::Join<o2::aod::FPi0s, o2::aod::FPi0Masks>;

  o2::framework::SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // setup for the two charged (pion) daughters. Dedicated pion-species selections, same ones
  // femtoProducerResonances.cxx uses for other two track resonances
  trackbuilder::ConfPionPlusSelection posDauSelection;
  o2::framework::Partition<FemtoTracks> posDauPartition = MAKE_RESONANCE_DAUGHTER_PARTITION(posDauSelection);
  o2::framework::Preslice<FemtoTracks> perColPosDau = o2::aod::femtobase::stored::fColId;

  trackbuilder::ConfPionMinusSelection negDauSelection;
  o2::framework::Partition<FemtoTracks> negDauPartition = MAKE_RESONANCE_DAUGHTER_PARTITION(negDauSelection);
  o2::framework::Preslice<FemtoTracks> perColNegDau = o2::aod::femtobase::stored::fColId;

  // setup for the neutral (pi0) daughter
  twophotonresonancebuilder::ConfPi0Selection neutralDauSelection;
  o2::framework::Partition<FemtoPi0s> neutralDauPartition = MAKE_TWOPHOTONRESONANCE_PARTITION(neutralDauSelection);
  o2::framework::Preslice<FemtoPi0s> perColNeutralDau = o2::aod::femtobase::stored::fColId;

  // resonance filter and builder
  trackstwophotonresonancebuilder::ConfOmegaMesonFilters confOmegaMesonFilter;
  trackstwophotonresonancebuilder::ConfTrackTrackTwoPhotonResonanceTables confOmegaMesonTables;
  trackstwophotonresonancebuilder::TrackTrackTwoPhotonResonanceBuilderProducts omegaMesonBuilderProducts;
  trackstwophotonresonancebuilder::TrackTrackTwoPhotonResonanceBuilder<modes::TrackTrackTwoPhotonResonance::kOmegaMeson> omegaMesonBuilder;

  void init(o2::framework::InitContext& context)
  {
    omegaMesonBuilder.init(confOmegaMesonFilter, posDauSelection, negDauSelection, confOmegaMesonTables, context);
  }

  void processOmegaMesons(FilteredFemtoCollision const& col, FemtoTracks const& /*tracks*/, FemtoPi0s const& /*pi0s*/)
  {
    omegaMesonBuilder.fillResonances(col, omegaMesonBuilderProducts, posDauPartition, negDauPartition, neutralDauPartition, cache);
  }
  PROCESS_SWITCH(FemtoProducerTrackTrackTwoPhotonResonances, processOmegaMesons, "Build OmegaMeson candidates", true);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{adaptAnalysisTask<FemtoProducerTrackTrackTwoPhotonResonances>(context)};
  return workflow;
}
