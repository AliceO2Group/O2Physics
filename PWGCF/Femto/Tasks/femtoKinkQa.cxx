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

using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femto;

struct FemtoKinkQa {

  // setup for collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);

  colhistmanager::CollisionHistManager<modes::Mode::kAnalysis_Qa> colHistManager;
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::ConfCollisionQaBinning confCollisionQaBinning;

  using FemtoCollisions = o2::soa::Join<FCols, FColMasks, FColPos, FColSphericities, FColMults>;
  using FemtoCollision = FemtoCollisions::iterator;

  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  // Define kink/sigma tables (joining tables for comprehensive information)
  using FemtoSigmas = o2::soa::Join<FSigmas, FSigmaMasks, FSigmaExtras>;
  using FemtoSigmaPlus = o2::soa::Join<FSigmaPlus, FSigmaPlusMasks, FSigmaPlusExtras>;
  using FemtoTracks = o2::soa::Join<FTracks, FTrackDcas, FTrackExtras, FTrackPids>;

  SliceCache cache;

  // setup for sigmas
  kinkbuilder::ConfSigmaSelection1 confSigmaSelection;

  Partition<FemtoSigmas> sigmaPartition = MAKE_SIGMA_PARTITION(confSigmaSelection);
  Preslice<FemtoSigmas> perColSigmas = femtobase::stored::fColId;

  kinkhistmanager::ConfSigmaBinning1 confSigmaBinning;
  kinkhistmanager::ConfSigmaQaBinning1 confSigmaQaBinning;
  kinkhistmanager::KinkHistManager<
    kinkhistmanager::PrefixSigmaQa,
    trackhistmanager::PrefixKinkChaDaughterQa,
    modes::Mode::kAnalysis_Qa,
    modes::Kink::kSigma>
    sigmaHistManager;

  // setup for sigma plus
  kinkbuilder::ConfSigmaPlusSelection1 confSigmaPlusSelection;

  Partition<FemtoSigmaPlus> sigmaPlusPartition = MAKE_SIGMAPLUS_PARTITION(confSigmaPlusSelection);
  Preslice<FemtoSigmaPlus> perColSigmaPlus = femtobase::stored::fColId;

  kinkhistmanager::ConfSigmaPlusBinning1 confSigmaPlusBinning;
  kinkhistmanager::ConfSigmaPlusQaBinning1 confSigmaPlusQaBinning;
  kinkhistmanager::KinkHistManager<
    kinkhistmanager::PrefixSigmaPlusQa,
    trackhistmanager::PrefixKinkChaDaughterQa,
    modes::Mode::kAnalysis_Qa,
    modes::Kink::kSigmaPlus>
    sigmaPlusHistManager;

  // setup for daughters
  trackhistmanager::ConfKinkChaDauBinning confKinkChaDaughterBinning;
  trackhistmanager::ConfKinkChaDauQaBinning confKinkChaDaughterQaBinning;

  HistogramRegistry hRegistry{"FemtoKinkQa", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    // create a map for histogram specs
    auto colHistSpec = colhistmanager::makeColQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
    colHistManager.init(&hRegistry, colHistSpec, confCollisionQaBinning);

    auto chaDauHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confKinkChaDaughterBinning, confKinkChaDaughterQaBinning);

    if ((doprocessSigma + doprocessSigmaPlus > 1)) {
      LOG(fatal) << "Only one process can be activated";
    }

    if (doprocessSigma) {
      auto sigmaHistSpec = kinkhistmanager::makeKinkQaHistSpecMap(confSigmaBinning, confSigmaQaBinning);
      sigmaHistManager.init(&hRegistry, sigmaHistSpec, confSigmaQaBinning, chaDauHistSpec, confKinkChaDaughterQaBinning);
    }

    if (doprocessSigmaPlus) {
      auto sigmaPlusHistSpec = kinkhistmanager::makeKinkQaHistSpecMap(confSigmaPlusBinning, confSigmaPlusQaBinning);
      sigmaPlusHistManager.init(&hRegistry, sigmaPlusHistSpec, confSigmaPlusQaBinning, chaDauHistSpec, confKinkChaDaughterQaBinning);
    }
  };

  void processSigma(FilteredFemtoCollision const& col, FemtoSigmas const& /*sigmas*/, FemtoTracks const& tracks)
  {
    colHistManager.fill(col);
    auto sigmaSlice = sigmaPartition->sliceByCached(femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& sigma : sigmaSlice) {
      sigmaHistManager.fill(sigma, tracks);
    }
  }
  PROCESS_SWITCH(FemtoKinkQa, processSigma, "Process sigmas", true);

  void processSigmaPlus(FilteredFemtoCollision const& col, FemtoSigmaPlus const& /*sigmaplus*/, FemtoTracks const& tracks)
  {
    colHistManager.fill(col);

    auto sigmaplusSlice = sigmaPlusPartition->sliceByCached(femtobase::stored::fColId, col.globalIndex(), cache);

    for (auto const& sp : sigmaplusSlice) {
      sigmaPlusHistManager.fill(sp, tracks);
    }
  }
  PROCESS_SWITCH(FemtoKinkQa, processSigmaPlus, "Process sigma plus", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoKinkQa>(cfgc),
  };
  return workflow;
}
