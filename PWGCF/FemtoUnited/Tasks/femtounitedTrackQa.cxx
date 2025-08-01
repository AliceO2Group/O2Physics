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

/// \file femtounitedTrackQa.cxx
/// \brief QA task for tracks
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/FemtoUnited/Core/collisionHistManager.h"
#include "PWGCF/FemtoUnited/Core/collisionSelection.h"
#include "PWGCF/FemtoUnited/Core/modes.h"
#include "PWGCF/FemtoUnited/Core/partitions.h"
#include "PWGCF/FemtoUnited/Core/trackHistManager.h"
#include "PWGCF/FemtoUnited/Core/trackSelection.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoCollisionsDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtounited;

struct FemtounitedTrackQa {

  // setup tables
  using Collisions = FUCols;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  using Tracks = o2::soa::Join<FUTracks, FUTrackMasks, FUTrackDCAs, FUTrackExtras, FUTrackPids>;

  SliceCache cache;

  // setup collisions
  collisionselection::ConfCollisionSelection collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::CollisionHistManager<modes::Mode::kANALYSIS_QA> colHistManager;

  // setup tracks
  trackselection::ConfTrackSelection1 trackSelections;
  trackhistmanager::ConfTrackBinning1 confTrackBinning;
  trackhistmanager::ConfTrackQaBinning1 confTrackQaBinning;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixTrackQa, modes::Mode::kANALYSIS_QA> trackHistManager;

  Partition<Tracks> trackPartition = MAKE_TRACK_PARTITION(trackSelections);
  Preslice<Tracks> perColReco = aod::femtobase::stored::collisionId;

  HistogramRegistry hRegistry{"FemtoTrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    // create a map for histogram specs
    auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
    colHistManager.init(&hRegistry, colHistSpec);
    auto trackHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confTrackBinning, confTrackQaBinning);
    trackHistManager.init(&hRegistry, trackHistSpec);
  };

  void process(FilteredCollision const& col, Tracks const& /*tracks*/)
  {
    colHistManager.fill(col);
    auto trackSlice = trackPartition->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
    for (auto const& track : trackSlice) {
      trackHistManager.fill(track);
    }
  }
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtounitedTrackQa>(cfgc),
  };
  return workflow;
}
