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
#include "PWGCF/Femto/Core/particleCleaner.h"
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

#include <map>
#include <vector>

using namespace o2::analysis::femto;

struct FemtoPairTrackTrack {

  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, o2::aod::FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks>;

  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, o2::aod::FTrackLabels>;

  o2::framework::SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // setup tracks
  trackbuilder::ConfTrackSelection1 confTrackSelections1;
  particlecleaner::ConfTrackCleaner1 confTrackCleaner1;
  trackhistmanager::ConfTrackBinning1 confTrackBinning1;

  trackbuilder::ConfTrackSelection2 confTrackSelections2;
  particlecleaner::ConfTrackCleaner2 confTrackCleaner2;
  trackhistmanager::ConfTrackBinning2 confTrackBinning2;

  o2::framework::Partition<FemtoTracks> trackPartition1 = MAKE_TRACK_PARTITION(confTrackSelections1);
  o2::framework::Partition<FemtoTracks> trackPartition2 = MAKE_TRACK_PARTITION(confTrackSelections2);
  o2::framework::Preslice<FemtoTracks> perColtracks = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoTracksWithLabel> trackWithLabelPartition1 = MAKE_TRACK_PARTITION(confTrackSelections1);
  o2::framework::Partition<FemtoTracksWithLabel> trackWithLabelPartition2 = MAKE_TRACK_PARTITION(confTrackSelections2);
  o2::framework::Preslice<FemtoTracksWithLabel> perColtracksWithLabel = o2::aod::femtobase::stored::fColId;

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
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}, true};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}, true};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult, o2::aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}, true};
  pairhistmanager::ConfMixing confMixing;

  o2::framework::HistogramRegistry hRegistry{"FemtoTrackTrack", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
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
    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> trackHistSpec1;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> trackHistSpec2;
    std::map<pairhistmanager::PairHist, std::vector<o2::framework::AxisSpec>> pairHistSpec;
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpec = closepairrejection::makeCprHistSpecMap(confCpr);

    if (processData) {
      colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
      trackHistSpec1 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning1);
      trackHistSpec2 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning2);
      pairHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning);
      pairTrackTrackBuilder.init<modes::Mode::kAnalysis>(&hRegistry, confTrackSelections1, confTrackSelections2, confTrackCleaner1, confTrackCleaner2, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec1, trackHistSpec2, pairHistSpec, cprHistSpec);
    } else {
      colHistSpec = colhistmanager::makeColMcHistSpecMap(confCollisionBinning);
      trackHistSpec1 = trackhistmanager::makeTrackMcHistSpecMap(confTrackBinning1);
      trackHistSpec2 = trackhistmanager::makeTrackMcHistSpecMap(confTrackBinning2);
      pairHistSpec = pairhistmanager::makePairMcHistSpecMap(confPairBinning);
      pairTrackTrackBuilder.init<modes::Mode::kAnalysis_Mc>(&hRegistry, confTrackSelections1, confTrackSelections2, confTrackCleaner1, confTrackCleaner2, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec1, trackHistSpec2, pairHistSpec, cprHistSpec);
    }
    hRegistry.print();
  };

  void processSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks)
  {
    pairTrackTrackBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, trackPartition1, trackPartition2, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processSameEvent, "Enable processing same event processing", true);

  void processSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackTrackBuilder.processSameEvent<modes::Mode::kAnalysis_Mc>(col, mcCols, tracks, trackWithLabelPartition1, trackWithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processSameEventMc, "Enable processing same event processing", false);

  void processMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks)
  {
    pairTrackTrackBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, trackPartition1, trackPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processMixedEvent, "Enable processing mixed event processing", true);

  void processMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackTrackBuilder.processMixedEvent<modes::Mode::kAnalysis_Mc>(cols, mcCols, tracks, trackWithLabelPartition1, trackWithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processMixedEventMc, "Enable processing mixed event processing", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairTrackTrack>(cfgc),
  };
  return workflow;
}
