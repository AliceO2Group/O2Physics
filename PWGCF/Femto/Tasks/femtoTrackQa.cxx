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

/// \file femtoTrackQa.cxx
/// \brief QA task for tracks
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/InitContext.h"
#include "Framework/OutputObjHeader.h"
#include "Framework/runDataProcessing.h"

#include <map>
#include <vector>

using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femto;

struct FemtoTrackQa {

  // setup tables
  using FemtoCollisions = o2::soa::Join<FCols, FColMasks, FColPos, FColSphericities, FColMults>;
  using FemtoCollision = FemtoCollisions::iterator;

  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoTracks = o2::soa::Join<FTracks, FTrackMasks, FTrackDcas, FTrackExtras, FTrackPids>;

  SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::ConfCollisionQaBinning confCollisionQaBinning;
  colhistmanager::CollisionHistManager<modes::Mode::kAnalysis_Qa> colHistManager;

  // setup tracks
  trackbuilder::ConfTrackSelection1 trackSelections;
  trackhistmanager::ConfTrackBinning1 confTrackBinning;
  trackhistmanager::ConfTrackQaBinning1 confTrackQaBinning;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixTrackQa, modes::Mode::kAnalysis_Qa> trackHistManager;

  Partition<FemtoTracks> trackPartition = MAKE_TRACK_PARTITION(trackSelections);
  Preslice<FemtoTracks> perColReco = femtobase::stored::fColId;

  HistogramRegistry hRegistry{"FemtoTrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    // create a map for histogram specs
    auto colHistSpec = colhistmanager::makeColQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
    colHistManager.init(&hRegistry, colHistSpec, confCollisionQaBinning);
    auto trackHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confTrackBinning, confTrackQaBinning);
    trackHistManager.init(&hRegistry, trackHistSpec, confTrackQaBinning, trackSelections.chargeAbs.value);
  };

  void process(FilteredFemtoCollision const& col, FemtoTracks const& tracks)
  {
    colHistManager.fill(col);
    auto trackSlice = trackPartition->sliceByCached(femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& track : trackSlice) {
      trackHistManager.fill(track, tracks);
    }
  }
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoTrackQa>(cfgc),
  };
  return workflow;
}
