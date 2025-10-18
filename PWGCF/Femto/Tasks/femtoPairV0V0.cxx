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

/// \file femtoPairV0V0.cxx
/// \brief Tasks that computes correlation between two v0s
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

struct FemtoPairV0V0 {

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

  // setup for daughters
  trackhistmanager::ConfV0PosDauBinning confPosDauBinning;
  trackhistmanager::ConfV0NegDauBinning confNegDauBinning;

  // setup lambdas
  v0builder::ConfLambdaSelection1 lambdaSelection;
  v0histmanager::ConfLambdaBinning1 confLambdaBinning;
  Partition<Lambdas> lambdaPartition = MAKE_LAMBDA_PARTITION(lambdaSelection);
  Preslice<Lambdas> perColLambdas = aod::femtobase::stored::collisionId;

  // setup k0shorts
  v0builder::ConfK0shortSelection1 k0shortSelection;
  v0histmanager::ConfK0shortBinning1 confK0shortBinning;
  Partition<K0shorts> k0shortPartition = MAKE_K0SHORT_PARTITION(k0shortSelection);
  Preslice<K0shorts> perColk0shorts = aod::femtobase::stored::collisionId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;

  pairbuilder::PairV0V0Builder<
    v0histmanager::PrefixLambda1,
    trackhistmanager::PrefixV0PosDaughter1,
    trackhistmanager::PrefixV0NegDaughter2,
    v0histmanager::PrefixLambda2,
    trackhistmanager::PrefixV0PosDaughter1,
    trackhistmanager::PrefixV0NegDaughter2,
    pairhistmanager::PrefixV0V0Se,
    pairhistmanager::PrefixV0V0Me,
    closepairrejection::PrefixV0V0PosSe,
    closepairrejection::PrefixV0V0NegSe,
    closepairrejection::PrefixV0V0PosMe,
    closepairrejection::PrefixV0V0NegMe,
    modes::V0::kLambda,
    modes::V0::kLambda,
    modes::Mode::kAnalysis>
    pairLambdaLambdaBuilder;

  pairbuilder::PairV0V0Builder<
    v0histmanager::PrefixK0short1,
    trackhistmanager::PrefixV0PosDaughter1,
    trackhistmanager::PrefixV0NegDaughter2,
    v0histmanager::PrefixK0short2,
    trackhistmanager::PrefixV0PosDaughter1,
    trackhistmanager::PrefixV0NegDaughter2,
    pairhistmanager::PrefixV0V0Se,
    pairhistmanager::PrefixV0V0Me,
    closepairrejection::PrefixV0V0PosSe,
    closepairrejection::PrefixV0V0NegSe,
    closepairrejection::PrefixV0V0PosMe,
    closepairrejection::PrefixV0V0NegMe,
    modes::V0::kK0short,
    modes::V0::kK0short,
    modes::Mode::kAnalysis>
    pairK0shortK0shortBuilder;

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
  closepairrejection::ConfCpr confCpr;

  void init(InitContext&)
  {

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins, confMixing.multBins.value}, true};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    // setup histograms
    auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
    auto posDauSpec = trackhistmanager::makeTrackHistSpecMap(confPosDauBinning);
    auto negDauSpec = trackhistmanager::makeTrackHistSpecMap(confNegDauBinning);
    auto cprHistSpec = closepairrejection::makeCprHistSpecMap(confCpr);

    // setup for lambda
    if (doprocessLambdaLambdaSameEvent || doprocessLambdaLambdaMixedEvent) {
      auto lambdaHistSpec = v0histmanager::makeV0HistSpecMap(confLambdaBinning);
      auto pairLambdaLambdaHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confLambdaBinning, confLambdaBinning);
      pairLambdaLambdaBuilder.init(&hRegistry, lambdaSelection, lambdaSelection, confCpr, confMixing, colHistSpec, lambdaHistSpec, lambdaHistSpec, posDauSpec, negDauSpec, pairLambdaLambdaHistSpec, cprHistSpec);
    }

    // setup for k0short
    if (doprocessK0shortK0shortSameEvent || doprocessK0shortK0shortMixedEvent) {
      auto k0shortHistSpec = v0histmanager::makeV0HistSpecMap(confK0shortBinning);
      auto pairLambdaLambdaHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confLambdaBinning, confLambdaBinning);
      pairLambdaLambdaBuilder.init(&hRegistry, k0shortSelection, k0shortSelection, confCpr, confMixing, colHistSpec, k0shortHistSpec, k0shortHistSpec, posDauSpec, negDauSpec, pairLambdaLambdaHistSpec, cprHistSpec);
    }

    if (((doprocessLambdaLambdaSameEvent || doprocessLambdaLambdaMixedEvent) + (doprocessK0shortK0shortSameEvent || doprocessK0shortK0shortMixedEvent)) > 1) {
      LOG(fatal) << "Can only process lambda-tracks Or k0short-tracks";
    }
  };

  void processLambdaLambdaSameEvent(FilteredCollision const& col, Tracks const& tracks, Lambdas const& lambdas)
  {
    pairLambdaLambdaBuilder.processSameEvent(col, tracks, lambdas, lambdaPartition, lambdaPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaLambdaSameEvent, "Enable processing same event processing for lambdas", true);

  void processLambdaLambdaMixedEvent(FilteredCollisions const& cols, Tracks const& tracks, Lambdas const& /*lambas*/)
  {
    pairLambdaLambdaBuilder.processMixedEvent(cols, tracks, lambdaPartition, lambdaPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaLambdaMixedEvent, "Enable processing mixed event processing for lambdas", true);

  void processK0shortK0shortSameEvent(FilteredCollision const& col, Tracks const& tracks, K0shorts const& k0shorts)
  {
    pairK0shortK0shortBuilder.processSameEvent(col, tracks, k0shorts, k0shortPartition, lambdaPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processK0shortK0shortSameEvent, "Enable processing same event processing for lambdas", true);

  void processK0shortK0shortMixedEvent(FilteredCollisions const& cols, Tracks const& tracks, K0shorts const& /*k0shorts*/)
  {
    pairK0shortK0shortBuilder.processMixedEvent(cols, tracks, k0shortPartition, k0shortPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processK0shortK0shortMixedEvent, "Enable processing mixed event processing for lambdas", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairV0V0>(cfgc),
  };
  return workflow;
}
