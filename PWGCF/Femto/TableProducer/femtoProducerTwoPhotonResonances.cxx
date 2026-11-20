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

/// \file femtoProducerTwoPhotonResonances.cxx
/// \brief Tasks that produces femto tables for resonances built from two PCM photons (pi0, eta, ...)
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/photonBuilder.h"
#include "PWGCF/Femto/Core/twoPhotonResonanceBuilder.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

using namespace o2::analysis::femto;

struct FemtoProducerTwoPhotonResonances {

  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoPhotons = o2::soa::Join<o2::aod::FPhotons, o2::aod::FPhotonMasks>;

  o2::framework::SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // setup for resonance daughter photons
  photonbuilder::ConfPhotonSelection1 confPhotonSelection;
  o2::framework::Partition<FemtoPhotons> photonPartition = MAKE_PHOTON_PARTITION(confPhotonSelection);
  o2::framework::Preslice<FemtoPhotons> perColPhotons = o2::aod::femtobase::stored::fColId;

  // resonance filters
  twophotonresonancebuilder::ConfPi0Filters confPi0Filter;
  twophotonresonancebuilder::ConfEtaFilters confEtaFilter;

  // resonance builders
  twophotonresonancebuilder::ConfTwoPhotonResonanceTables confTwoPhotonResonanceTables;
  twophotonresonancebuilder::TwoPhotonResonanceBuilderProducts twoPhotonResonanceBuilderProducts;
  twophotonresonancebuilder::TwoPhotonResonanceBuilder<modes::TwoPhotonResonance::kPi0> pi0Builder;
  twophotonresonancebuilder::TwoPhotonResonanceBuilder<modes::TwoPhotonResonance::kEta> etaBuilder;

  void init(o2::framework::InitContext& context)
  {
    // init builders
    pi0Builder.init(confPi0Filter, confTwoPhotonResonanceTables, context);
    etaBuilder.init(confEtaFilter, confTwoPhotonResonanceTables, context);
  }

  // process functions
  void processPi0(FilteredFemtoCollision const& col, FemtoPhotons const& /*photons*/)
  {
    pi0Builder.fillResonances(col, twoPhotonResonanceBuilderProducts, photonPartition, cache);
  }
  PROCESS_SWITCH(FemtoProducerTwoPhotonResonances, processPi0, "Build Pi0 candidates", true);

  void processEta(FilteredFemtoCollision const& col, FemtoPhotons const& /*photons*/)
  {
    etaBuilder.fillResonances(col, twoPhotonResonanceBuilderProducts, photonPartition, cache);
  }
  PROCESS_SWITCH(FemtoProducerTwoPhotonResonances, processEta, "Build Eta candidates", true);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{adaptAnalysisTask<FemtoProducerTwoPhotonResonances>(context)};
  return workflow;
}
