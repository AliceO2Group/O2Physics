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
  trackbuilder::ConfTrackSelection1 confTrackSelection;
  trackhistmanager::ConfTrackBinning1 confTrackBinning;

  Partition<FemtoTracks> trackPartition = MAKE_TRACK_PARTITION(confTrackSelection);
  Preslice<FemtoTracks> perColTracks = aod::femtobase::stored::fColId;

  Partition<FemtoTracksWithLabel> trackWithLabelPartition = MAKE_TRACK_PARTITION(confTrackSelection);
  Preslice<FemtoTracksWithLabel> perColtracksWithLabel = aod::femtobase::stored::fColId;

  // setup for daughters
  trackhistmanager::ConfV0PosDauBinning confPosDauBinning;
  trackhistmanager::ConfV0NegDauBinning confNegDauBinning;

  // setup lambdas
  v0builder::ConfLambdaSelection1 lambdaSelection;
  v0histmanager::ConfLambdaBinning1 confLambdaBinning;

  Partition<FemtoLambdas> lambdaPartition = MAKE_LAMBDA_PARTITION(lambdaSelection);
  Preslice<FemtoLambdas> perColLambdas = aod::femtobase::stored::fColId;

  Partition<FemtoLambdasWithLabel> lambdaWithLabelPartition = MAKE_LAMBDA_PARTITION(lambdaSelection);
  Preslice<FemtoLambdasWithLabel> perColLambdasWithLabel = aod::femtobase::stored::fColId;

  // setup k0shorts
  v0builder::ConfK0shortSelection1 k0shortSelection;
  v0histmanager::ConfK0shortBinning1 confK0shortBinning;

  Partition<FemtoK0shorts> k0shortPartition = MAKE_K0SHORT_PARTITION(k0shortSelection);
  Preslice<FemtoK0shorts> perColk0shorts = aod::femtobase::stored::fColId;

  Partition<FemtoK0shortsWithLabel> k0shortWithLabelPartition = MAKE_K0SHORT_PARTITION(k0shortSelection);
  Preslice<FemtoK0shortsWithLabel> perColk0shortsWithLabel = aod::femtobase::stored::fColId;

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
    bool processData = doprocessLambdaSameEvent || doprocessLambdaMixedEvent || doprocessK0shortSameEvent || doprocessK0shortMixedEvent;
    bool processMc = doprocessLambdaSameEventMc || doprocessLambdaMixedEventMc || doprocessK0shortSameEventMc || doprocessK0shortMixedEventMc;

    if (processData && processMc) {
      LOG(fatal) << "Both data and mc processing is enabled. Breaking...";
    }

    bool processLambda = doprocessLambdaSameEvent || doprocessLambdaSameEventMc || doprocessLambdaMixedEvent || doprocessLambdaMixedEventMc;
    bool processK0short = doprocessK0shortSameEvent || doprocessK0shortSameEventMc || doprocessK0shortMixedEvent || doprocessK0shortMixedEventMc;

    if (processLambda && processK0short) {
      LOG(fatal) << "Both lambda-track and k0short-track processing is enabled. Breaking...";
    }

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins, confMixing.multBins.value}, true};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    auto cprHistSpec = closepairrejection::makeCprHistSpecMap(confCpr);

    if (processData) {
      auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
      auto trackHistSpec = trackhistmanager::makeTrackHistSpecMap(confTrackBinning);
      auto posDauSpec = trackhistmanager::makeTrackHistSpecMap(confPosDauBinning);
      auto negDauSpec = trackhistmanager::makeTrackHistSpecMap(confNegDauBinning);
      auto pairTrackV0HistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning);
      if (processLambda) {
        auto lambdaHistSpec = v0histmanager::makeV0HistSpecMap(confLambdaBinning);
        pairTrackLambdaBuilder.init<modes::Mode::kAnalysis>(&hRegistry, confTrackSelection, lambdaSelection, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, lambdaHistSpec, posDauSpec, negDauSpec, pairTrackV0HistSpec, cprHistSpec);
      } else {
        auto k0shortHistSpec = v0histmanager::makeV0HistSpecMap(confK0shortBinning);
        pairTrackK0shortBuilder.init<modes::Mode::kAnalysis>(&hRegistry, confTrackSelection, lambdaSelection, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, k0shortHistSpec, posDauSpec, negDauSpec, pairTrackV0HistSpec, cprHistSpec);
      }
    } else {
      auto colHistSpec = colhistmanager::makeColMcHistSpecMap(confCollisionBinning);
      auto trackHistSpec = trackhistmanager::makeTrackMcHistSpecMap(confTrackBinning);
      auto posDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confPosDauBinning);
      auto negDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confNegDauBinning);
      auto pairTrackV0HistSpec = pairhistmanager::makePairMcHistSpecMap(confPairBinning);
      if (processLambda) {
        auto lambdaHistSpec = v0histmanager::makeV0McHistSpecMap(confLambdaBinning);
        pairTrackLambdaBuilder.init<modes::Mode::kAnalysis_Mc>(&hRegistry, confTrackSelection, lambdaSelection, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, lambdaHistSpec, posDauSpec, negDauSpec, pairTrackV0HistSpec, cprHistSpec);
      } else {
        auto k0shortHistSpec = v0histmanager::makeV0McHistSpecMap(confK0shortBinning);
        pairTrackK0shortBuilder.init<modes::Mode::kAnalysis_Mc>(&hRegistry, confTrackSelection, lambdaSelection, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, k0shortHistSpec, posDauSpec, negDauSpec, pairTrackV0HistSpec, cprHistSpec);
      }
    }
    hRegistry.print();
  };

  void processLambdaSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoLambdas const& lambdas)
  {
    pairTrackLambdaBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, trackPartition, lambdas, lambdaPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackV0, processLambdaSameEvent, "Enable processing same event processing for tracks and lambdas", true);

  void processLambdaSameEventMc(FilteredFemtoCollisionWithLabel const& col, FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& lambdas, FMcParticles const& mcParticles, FMcMothers const& mcMothers, FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackLambdaBuilder.processSameEvent<modes::Mode::kAnalysis_Mc>(col, mcCols, tracks, trackWithLabelPartition, lambdas, lambdaWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackV0, processLambdaSameEventMc, "Enable processing same event processing for tracks and lambdas with MC information", false);

  void processLambdaMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoLambdas const& /*lambas*/)
  {
    pairTrackLambdaBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, trackPartition, lambdaPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackV0, processLambdaMixedEvent, "Enable processing mixed event processing for tracks and lambdas", true);

  void processLambdaMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& /*lambas*/, FMcParticles const& mcParticles)
  {
    pairTrackLambdaBuilder.processMixedEvent<modes::Mode::kAnalysis_Mc>(cols, mcCols, tracks, trackWithLabelPartition, lambdaWithLabelPartition, mcParticles, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackV0, processLambdaMixedEventMc, "Enable processing mixed event processing for tracks and lambdas with MC information", false);

  void processK0shortSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoK0shorts const& k0shorts)
  {
    pairTrackK0shortBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, trackPartition, k0shorts, k0shortPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackV0, processK0shortSameEvent, "Enable processing same event processing for tracks and k0shorts", false);

  void processK0shortSameEventMc(FilteredFemtoCollisionWithLabel const& col, FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoK0shortsWithLabel const& k0shorts, FMcParticles const& mcParticles, FMcMothers const& mcMothers, FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackK0shortBuilder.processSameEvent<modes::Mode::kAnalysis_Mc>(col, mcCols, tracks, trackWithLabelPartition, k0shorts, k0shortWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackV0, processK0shortSameEventMc, "Enable processing same event processing for tracks and k0shorts with MC information", false);

  void processK0shortMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoK0shorts const& /*k0shorts*/)
  {
    pairTrackK0shortBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, trackPartition, k0shortPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackV0, processK0shortMixedEvent, "Enable processing mixed event processing for tracks and k0shorts", false);

  void processK0shortMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoK0shortsWithLabel const& /*k0shorts*/, FMcParticles const& mcParticles)
  {
    pairTrackK0shortBuilder.processMixedEvent<modes::Mode::kAnalysis_Mc>(cols, mcCols, tracks, trackWithLabelPartition, k0shortWithLabelPartition, mcParticles, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackV0, processK0shortMixedEventMc, "Enable processing mixed event processing for tracks and k0shorts with mc information", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairTrackV0>(cfgc),
  };
  return workflow;
}
