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

/// \file femtoPairTrackV0.cxx
/// \brief Tasks that computes correlation between tracks and v0s
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

struct FemtoPairTrackV0 {

  // setup tables
  using Collisions = Join<FCols, FColMasks>;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  using Tracks = o2::soa::Join<FTracks, FTrackMasks>;
  using Lambdas = o2::soa::Join<FLambdas, FLambdaMasks>;
  using K0shorts = o2::soa::Join<FK0shorts, FK0shortMasks>;

  SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // setup tracks
  trackbuilder::ConfTrackSelection1 trackSelection;
  trackhistmanager::ConfTrackBinning1 confTrackBinning;
  Partition<Tracks> trackPartition = MAKE_TRACK_PARTITION(trackSelection);
  Preslice<Tracks> perColTracks = aod::femtobase::stored::fColId;

  // setup for daughters
  trackhistmanager::ConfV0PosDauBinning confPosDauBinning;
  trackhistmanager::ConfV0NegDauBinning confNegDauBinning;

  // setup lambdas
  v0builder::ConfLambdaSelection1 lambdaSelection;
  v0histmanager::ConfLambdaBinning1 confLambdaBinning;
  Partition<Lambdas> lambdaPartition = MAKE_LAMBDA_PARTITION(lambdaSelection);
  Preslice<Lambdas> perColLambdas = aod::femtobase::stored::fColId;

  // setup k0shorts
  v0builder::ConfK0shortSelection1 k0shortSelection;
  v0histmanager::ConfK0shortBinning1 confK0shortBinning;
  Partition<K0shorts> k0shortPartition = MAKE_K0SHORT_PARTITION(k0shortSelection);
  Preslice<K0shorts> perColk0shorts = aod::femtobase::stored::fColId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::ConfPairCuts confPairCuts;

  pairbuilder::PairTrackV0Builder<
    trackhistmanager::PrefixTrack1,
    v0histmanager::PrefixLambda1,
    trackhistmanager::PrefixV01PosDaughter,
    trackhistmanager::PrefixV01NegDaughter,
    pairhistmanager::PrefixTrackV0Se,
    pairhistmanager::PrefixTrackV0Me,
    closepairrejection::PrefixTrackV0DaughterSe,
    closepairrejection::PrefixTrackV0DaughterMe,
    modes::V0::kLambda>
    pairTrackLambdaBuilder;

  pairbuilder::PairTrackV0Builder<
    trackhistmanager::PrefixTrack1,
    v0histmanager::PrefixK0short1,
    trackhistmanager::PrefixV01PosDaughter,
    trackhistmanager::PrefixV01NegDaughter,
    pairhistmanager::PrefixTrackV0Se,
    pairhistmanager::PrefixTrackV0Me,
    closepairrejection::PrefixTrackV0DaughterSe,
    closepairrejection::PrefixTrackV0DaughterMe,
    modes::V0::kK0short>
    pairTrackK0shortBuilder;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  ColumnBinningPolicy<femtocollisions::PosZ, femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}, true};
  ColumnBinningPolicy<aod::femtocollisions::PosZ, aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}, true};
  ColumnBinningPolicy<aod::femtocollisions::PosZ, aod::femtocollisions::Mult, aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}, true};
  pairhistmanager::ConfMixing confMixing;

  HistogramRegistry hRegistry{"FemtoTrackV0", {}, OutputObjHandlingPolicy::AnalysisObject};

  // setup cpr
  closepairrejection::ConfCprTrackV0Daughter confCpr;

  void init(InitContext&)
  {

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins, confMixing.multBins.value}, true};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    // setup histograms
    auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
    auto trackHistSpec = trackhistmanager::makeTrackHistSpecMap(confTrackBinning);
    auto posDauSpec = trackhistmanager::makeTrackHistSpecMap(confPosDauBinning);
    auto negDauSpec = trackhistmanager::makeTrackHistSpecMap(confNegDauBinning);
    auto cprHistSpec = closepairrejection::makeCprHistSpecMap(confCpr);

    // setup for lambda
    if (doprocessLambdaSameEvent || doprocessLambdaMixedEvent) {
      auto lambdaHistSpec = v0histmanager::makeV0HistSpecMap(confLambdaBinning);
      auto pairTrackLambdaHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning);
      pairTrackLambdaBuilder.init<modes::Mode::kAnalysis>(&hRegistry, trackSelection, lambdaSelection, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, lambdaHistSpec, posDauSpec, negDauSpec, pairTrackLambdaHistSpec, cprHistSpec);
    }

    // setup for k0short
    if (doprocessK0shortSameEvent || doprocessK0shortMixedEvent) {
      auto k0shortHistSpec = v0histmanager::makeV0HistSpecMap(confK0shortBinning);
      auto pairTrackK0shortHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning);
      pairTrackK0shortBuilder.init<modes::Mode::kAnalysis>(&hRegistry, trackSelection, lambdaSelection, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, k0shortHistSpec, posDauSpec, negDauSpec, pairTrackK0shortHistSpec, cprHistSpec);
    }

    if (((doprocessLambdaSameEvent || doprocessLambdaMixedEvent) + (doprocessK0shortSameEvent || doprocessK0shortMixedEvent)) > 1) {
      LOG(fatal) << "Can only process lambda-tracks Or k0short-tracks";
    }
  };

  void processLambdaSameEvent(FilteredCollision const& col, Tracks const& tracks, Lambdas const& lambdas)
  {
    pairTrackLambdaBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, trackPartition, lambdas, lambdaPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackV0, processLambdaSameEvent, "Enable processing same event processing for tracks and lambdas", true);

  void processLambdaMixedEvent(FilteredCollisions const& cols, Tracks const& tracks, Lambdas const& /*lambas*/)
  {
    pairTrackLambdaBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, trackPartition, lambdaPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackV0, processLambdaMixedEvent, "Enable processing mixed event processing for tracks and lambdas", true);
  //
  void processK0shortSameEvent(FilteredCollision const& col, Tracks const& tracks, K0shorts const& k0shorts)
  {
    pairTrackK0shortBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, trackPartition, k0shorts, k0shortPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackV0, processK0shortSameEvent, "Enable processing same event processing for tracks and k0shorts", false);

  void processK0shortMixedEvent(FilteredCollisions const& cols, Tracks const& tracks, K0shorts const& /*k0shorts*/)
  {
    pairTrackK0shortBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, trackPartition, k0shortPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackV0, processK0shortMixedEvent, "Enable processing mixed event processing for tracks and k0shorts", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairTrackV0>(cfgc),
  };
  return workflow;
}
