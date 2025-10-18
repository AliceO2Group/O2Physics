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

/// \file femtoKinkQa.cxx
/// \brief QA task for kinks
/// \author Anton Riedel, TU München, anton.riedel@cern.ch
/// \author Henrik Fribert, TU München, henrik.fribert@cern.ch

#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/kinkBuilder.h"
#include "PWGCF/Femto/Core/kinkHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"
#include "PWGLF/DataModel/LFKinkDecayTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/InitContext.h"
#include "Framework/OutputObjHeader.h"
#include "Framework/runDataProcessing.h"

#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femto;

struct FemtoKinkQa {

  // setup for collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);

  colhistmanager::CollisionHistManager<modes::Mode::kAnalysis_Qa> colHistManager;
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // using Collisions = o2::soa::Join<FUCols, FUColPos, FUColMults, FUColCents>;
  using Collisions = Join<FCols, FColMasks, FColPos, FColSphericities, FColMults>;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  // Define kink/sigma tables (joining tables for comprehensive information)
  using Sigmas = o2::soa::Join<FSigmas, FSigmaMasks, FSigmaExtras>;
  using Tracks = o2::soa::Join<FTracks, FTrackDcas, FTrackExtras, FTrackPids>;

  SliceCache cache;

  // setup for sigmas
  kinkbuilder::ConfSigmaSelection1 confSigmaSelection;

  Partition<Sigmas> sigmaPartition = MAKE_SIGMA_PARTITION(confSigmaSelection);
  Preslice<Sigmas> perColSigmas = aod::femtobase::stored::collisionId;

  kinkhistmanager::ConfSigmaBinning1 confSigmaBinning;
  kinkhistmanager::ConfSigmaQaBinning1 confSigmaQaBinning;
  kinkhistmanager::KinkHistManager<
    kinkhistmanager::PrefixSigmaQa,
    trackhistmanager::PrefixKinkChaDaughterQa,
    modes::Mode::kAnalysis_Qa,
    modes::Kink::kSigma>
    sigmaHistManager;

  // setup for daughters
  trackhistmanager::ConfKinkChaDauBinning confKinkChaDaughterBinning;
  trackhistmanager::ConfKinkChaDauQaBinning confKinkChaDaughterQaBinning;

  HistogramRegistry hRegistry{"FemtoKinkQa", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    auto sigmaHistSpec = kinkhistmanager::makeKinkQaHistSpecMap(confSigmaBinning, confSigmaQaBinning);
    auto chaDauHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confKinkChaDaughterBinning, confKinkChaDaughterQaBinning);

    sigmaHistManager.init(&hRegistry, sigmaHistSpec, chaDauHistSpec);

    auto collisionHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
    colHistManager.init(&hRegistry, collisionHistSpec);
  };

  // Process function for sigma particles from femto tables
  void processSigma(FilteredCollision const& col, Sigmas const& /*sigmas*/, Tracks const& tracks)
  {
    colHistManager.fill(col);
    auto sigmaSlice = sigmaPartition->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
    for (auto const& sigma : sigmaSlice) {
      sigmaHistManager.fill(sigma, tracks);
    }
  }
  PROCESS_SWITCH(FemtoKinkQa, processSigma, "Process sigmas", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoKinkQa>(cfgc),
  };
  return workflow;
}
