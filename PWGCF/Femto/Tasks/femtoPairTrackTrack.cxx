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

/// \file femtoPairTrackTrack.cxx
/// \brief Tasks that computes correlation between two tracks
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/closePairRejection.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/pairBuilder.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femto;

struct FemtoPairTrackTrack {

  // setup tables
  using Collisions = Join<FCols, FColMasks>;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  using Tracks = o2::soa::Join<FTracks, FTrackMasks>;

  SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // setup tracks
  trackbuilder::ConfTrackSelection1 trackSelections1;
  trackhistmanager::ConfTrackBinning1 confTrackBinning1;
  trackbuilder::ConfTrackSelection2 trackSelections2;
  trackhistmanager::ConfTrackBinning2 confTrackBinning2;

  Partition<Tracks> trackPartition1 = MAKE_TRACK_PARTITION(trackSelections1);
  Partition<Tracks> trackPartition2 = MAKE_TRACK_PARTITION(trackSelections2);

  Preslice<Tracks> perColReco = aod::femtobase::stored::fColId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  closepairrejection::ConfCpr confCpr;

  pairbuilder::PairTrackTrackBuilder<
    trackhistmanager::PrefixTrack1,
    trackhistmanager::PrefixTrack2,
    pairhistmanager::PrefixTrackTrackSe,
    pairhistmanager::PrefixTrackTrackMe,
    closepairrejection::PrefixTrackTrackSe,
    closepairrejection::PrefixTrackTrackMe,
    modes::Mode::kAnalysis>
    pairTrackTrackBuilder;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  ColumnBinningPolicy<femtocollisions::PosZ, femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}, true};
  ColumnBinningPolicy<aod::femtocollisions::PosZ, aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}, true};
  ColumnBinningPolicy<aod::femtocollisions::PosZ, aod::femtocollisions::Mult, aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}, true};
  pairhistmanager::ConfMixing confMixing;

  HistogramRegistry hRegistry{"FemtoTrackTrack", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins, confMixing.multBins.value}, true};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    // setup histogram specs
    auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
    auto trackHistSpec1 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning1);
    auto trackHistSpec2 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning2);
    auto pairHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confTrackBinning1, confTrackBinning2);
    auto cprHistSpec = closepairrejection::makeCprHistSpecMap(confCpr);

    pairTrackTrackBuilder.init(&hRegistry, trackSelections1, trackSelections2, confCpr, confMixing, colHistSpec, trackHistSpec1, trackHistSpec2, pairHistSpec, cprHistSpec);
  };

  void processSameEvent(FilteredCollision const& col, Tracks const& tracks)
  {
    pairTrackTrackBuilder.processSameEvent(col, tracks, trackPartition1, trackPartition2, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processSameEvent, "Enable processing same event processing", true);

  void processMixedEvent(FilteredCollisions const& cols, Tracks const& tracks)
  {
    pairTrackTrackBuilder.processMixedEvent(cols, tracks, trackPartition1, trackPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processMixedEvent, "Enable processing mixed event processing", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairTrackTrack>(cfgc),
  };
  return workflow;
}
