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
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/pairBuilder.h"
#include "PWGCF/Femto/Core/pairHistManager.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/BinningPolicy.h"
#include "Framework/Configurable.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/InitContext.h"
#include "Framework/OutputObjHeader.h"
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
  using FemtoCollisions = Join<FCols, FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoTracks = o2::soa::Join<FTracks, FTrackMasks>;

  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, FTrackLabels>;

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

  Partition<FemtoTracks> trackPartition1 = MAKE_TRACK_PARTITION(trackSelections1);
  Partition<FemtoTracks> trackPartition2 = MAKE_TRACK_PARTITION(trackSelections2);

  Preslice<FemtoTracks> perColReco = aod::femtobase::stored::fColId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::ConfPairCuts confPairCuts;

  closepairrejection::ConfCprTrackTrack confCpr;

  pairbuilder::PairTrackTrackBuilder<
    trackhistmanager::PrefixTrack1,
    trackhistmanager::PrefixTrack2,
    pairhistmanager::PrefixTrackTrackSe,
    pairhistmanager::PrefixTrackTrackMe,
    closepairrejection::PrefixTrackTrackSe,
    closepairrejection::PrefixTrackTrackMe>
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
    auto pairHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning);
    auto cprHistSpec = closepairrejection::makeCprHistSpecMap(confCpr);

    if ((doprocessSameEvent + doprocessSameEventMc) > 1 || (doprocessMixedEvent + doprocessMixedEventMc) > 1) {
      LOG(fatal) << "More than 1 same or mixed event process function is activated. Breaking...";
    }

    bool processData = doprocessSameEvent || doprocessMixedEvent;
    bool processMc = doprocessSameEventMc || doprocessMixedEventMc;

    if (processData && processMc) {
      LOG(fatal) << "Both data and mc processing is activated. Breaking...";
    }

    if (processData) {
      pairTrackTrackBuilder.init<modes::Mode::kAnalysis>(&hRegistry, trackSelections1, trackSelections2, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec1, trackHistSpec2, pairHistSpec, cprHistSpec);
    } else {
      pairTrackTrackBuilder.init<modes::Mode::kAnalysis_Mc>(&hRegistry, trackSelections1, trackSelections2, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec1, trackHistSpec2, pairHistSpec, cprHistSpec);
    }
  };

  void processSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks)
  {
    pairTrackTrackBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, trackPartition1, trackPartition2, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processSameEvent, "Enable processing same event processing", true);

  void processSameEventMc(FilteredFemtoCollision const& col, FemtoTracks const& tracks)
  {
    pairTrackTrackBuilder.processSameEvent<modes::Mode::kAnalysis_Mc>(col, tracks, trackPartition1, trackPartition2, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processSameEventMc, "Enable processing same event processing", false);

  void processMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks)
  {
    pairTrackTrackBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, trackPartition1, trackPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processMixedEvent, "Enable processing mixed event processing", true);

  void processMixedEventMc(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks)
  {
    pairTrackTrackBuilder.processMixedEvent<modes::Mode::kAnalysis_Mc>(cols, tracks, trackPartition1, trackPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processMixedEventMc, "Enable processing mixed event processing", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairTrackTrack>(cfgc),
  };
  return workflow;
}
