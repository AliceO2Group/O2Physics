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

/// \file femtoTripletTrackTrackTrack.cxx
/// \brief Tasks that computes correlation between two tracks
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/closeTripletRejection.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/Core/tripletBuilder.h"
#include "PWGCF/Femto/Core/tripletHistManager.h"
#include "PWGCF/Femto/Core/v0Builder.h"
#include "PWGCF/Femto/Core/v0HistManager.h"
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

struct FemtoTripletTrackTrackTrack {

  // setup tables
  using FemtoCollisions = Join<FCols, FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoTracks = o2::soa::Join<FTracks, FTrackMasks>;
  using FemtoLambdas = o2::soa::Join<FLambdas, FLambdaMasks>;
  using FemtoK0shorts = o2::soa::Join<FK0shorts, FK0shortMasks>;

  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, FTrackLabels>;
  using FemtoLambdasWithLabel = o2::soa::Join<FemtoLambdas, FLambdaLabels>;
  using FemtoK0shortsWithLabel = o2::soa::Join<FemtoK0shorts, FK0shortLabels>;

  SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // setup tracks
  trackbuilder::ConfTrackSelection1 confTrackSelections1;
  trackhistmanager::ConfTrackBinning1 confTrackBinning1;
  trackbuilder::ConfTrackSelection2 confTrackSelections2;
  trackhistmanager::ConfTrackBinning2 confTrackBinning2;
  trackbuilder::ConfTrackSelection3 confTrackSelections3;
  trackhistmanager::ConfTrackBinning3 confTrackBinning3;

  Partition<FemtoTracks> trackPartition1 = MAKE_TRACK_PARTITION(confTrackSelections1);
  Partition<FemtoTracks> trackPartition2 = MAKE_TRACK_PARTITION(confTrackSelections2);
  Partition<FemtoTracks> trackPartition3 = MAKE_TRACK_PARTITION(confTrackSelections3);
  Preslice<FemtoTracks> perColtracks = aod::femtobase::stored::fColId;

  Partition<FemtoTracksWithLabel> trackWithLabelPartition1 = MAKE_TRACK_PARTITION(confTrackSelections1);
  Partition<FemtoTracksWithLabel> trackWithLabelPartition2 = MAKE_TRACK_PARTITION(confTrackSelections2);
  Partition<FemtoTracksWithLabel> trackWithLabelPartition3 = MAKE_TRACK_PARTITION(confTrackSelections3);
  Preslice<FemtoTracksWithLabel> perColtracksWithLabel = aod::femtobase::stored::fColId;

  // setup triplets
  triplethistmanager::ConfTripletBinning confTripletBinning;
  triplethistmanager::ConfTripletCuts confTripletCuts;

  closetripletrejection::ConfCtrTrackTrackTrack confCtr;

  tripletbuilder::TripletTrackTrackTrackBuilder<
    trackhistmanager::PrefixTrack1,
    trackhistmanager::PrefixTrack2,
    trackhistmanager::PrefixTrack3,
    triplethistmanager::PrefixTrackTrackTrackSe,
    triplethistmanager::PrefixTrackTrackTrackMe,
    closetripletrejection::PrefixTrack1Track2Se,
    closetripletrejection::PrefixTrack2Track3Se,
    closetripletrejection::PrefixTrack1Track3Se,
    closetripletrejection::PrefixTrack1Track2Me,
    closetripletrejection::PrefixTrack2Track3Me,
    closetripletrejection::PrefixTrack1Track3Me>
    tripletTrackTrackTrackBuilder;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  ColumnBinningPolicy<femtocollisions::PosZ, femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}, true};
  ColumnBinningPolicy<aod::femtocollisions::PosZ, aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}, true};
  ColumnBinningPolicy<aod::femtocollisions::PosZ, aod::femtocollisions::Mult, aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}, true};
  triplethistmanager::ConfMixing confMixing;

  HistogramRegistry hRegistry{"FemtoTrackTrackTrack", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    if ((doprocessSameEvent + doprocessSameEventMc) > 1 || (doprocessMixedEvent + doprocessMixedEventMc) > 1) {
      LOG(fatal) << "More than 1 same or mixed event process function is activated. Breaking...";
    }
    bool processData = doprocessSameEvent || doprocessMixedEvent;
    bool processMc = doprocessSameEventMc || doprocessMixedEventMc;
    if (processData && processMc) {
      LOG(fatal) << "Both data and mc processing is activated. Breaking...";
    }

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins, confMixing.multBins.value}, true};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    // setup histogram specs
    auto ctrHistSpec = closepairrejection::makeCprHistSpecMap(confCtr);

    if (processData) {
      auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
      auto trackHistSpec1 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning1);
      auto trackHistSpec2 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning2);
      auto trackHistSpec3 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning3);
      auto tripletHistSpec = triplethistmanager::makeTripletHistSpecMap(confTripletBinning);
      tripletTrackTrackTrackBuilder.init<modes::Mode::kAnalysis>(&hRegistry, confTrackSelections1, confTrackSelections2, confTrackSelections3, confCtr, confMixing, confTripletBinning, confTripletCuts, colHistSpec, trackHistSpec1, trackHistSpec2, trackHistSpec3, tripletHistSpec, ctrHistSpec);
    } else {
      auto colHistSpec = colhistmanager::makeColMcHistSpecMap(confCollisionBinning);
      auto trackHistSpec1 = trackhistmanager::makeTrackMcHistSpecMap(confTrackBinning1);
      auto trackHistSpec2 = trackhistmanager::makeTrackMcHistSpecMap(confTrackBinning2);
      auto trackHistSpec3 = trackhistmanager::makeTrackMcHistSpecMap(confTrackBinning3);
      auto tripletHistSpec = triplethistmanager::makeTripletMcHistSpecMap(confTripletBinning);
      tripletTrackTrackTrackBuilder.init<modes::Mode::kAnalysis_Mc>(&hRegistry, confTrackSelections1, confTrackSelections2, confTrackSelections3, confCtr, confMixing, confTripletBinning, confTripletCuts, colHistSpec, trackHistSpec1, trackHistSpec2, trackHistSpec3, tripletHistSpec, ctrHistSpec);
    }
    hRegistry.print();
  };

  void processSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks)
  {
    tripletTrackTrackTrackBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, trackPartition1, trackPartition2, trackPartition3, cache);
  }
  PROCESS_SWITCH(FemtoTripletTrackTrackTrack, processSameEvent, "Enable processing same event processing", true);

  void processSameEventMc(FilteredFemtoCollisionWithLabel const& col, FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FMcParticles const& mcParticles, FMcMothers const& mcMothers, FMcPartMoths const& mcPartonicMothers)
  {
    tripletTrackTrackTrackBuilder.processSameEvent<modes::Mode::kAnalysis_Mc>(col, mcCols, tracks, trackWithLabelPartition1, trackWithLabelPartition2, trackWithLabelPartition3, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoTripletTrackTrackTrack, processSameEventMc, "Enable processing same event processing", false);

  void processMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks)
  {
    tripletTrackTrackTrackBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, trackPartition1, trackPartition2, trackPartition3, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoTripletTrackTrackTrack, processMixedEvent, "Enable processing mixed event processing", true);

  void processMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FMcParticles const& mcParticles)
  {
    tripletTrackTrackTrackBuilder.processMixedEvent<modes::Mode::kAnalysis_Mc>(cols, mcCols, tracks, trackWithLabelPartition1, trackWithLabelPartition2, trackWithLabelPartition3, mcParticles, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoTripletTrackTrackTrack, processMixedEventMc, "Enable processing mixed event processing", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoTripletTrackTrackTrack>(cfgc),
  };
  return workflow;
}
