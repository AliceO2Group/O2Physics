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

/// \file femtounitedPairTrackV0.cxx
/// \brief Tasks that computes correlation between tracks and lambdas
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/closePairRejection.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/pairCleaner.h"
#include "PWGCF/Femto/Core/pairHistManager.h"
#include "PWGCF/Femto/Core/pairProcessHelpers.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/Core/v0Builder.h"
#include "PWGCF/Femto/Core/v0HistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include <map>
#include <random>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femto;

struct FemtounitedPairTrackV0 {

  // setup tables
  using Collisions = FCols;
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
  colhistmanager::CollisionHistManager<modes::Mode::kAnalysis> colHistManager;

  // setup tracks
  trackbuilder::ConfTrackSelection1 trackSelections;
  trackhistmanager::ConfTrackBinning1 confTrackBinning;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixTrack1, modes::Mode::kAnalysis> trackHistManager;
  Partition<Tracks> trackPartition = MAKE_TRACK_PARTITION(trackSelections);
  Preslice<Tracks> perColTracks = aod::femtobase::stored::collisionId;

  // setup for daughters
  trackhistmanager::ConfV0PosDauBinning confPosDauBinning;
  trackhistmanager::ConfV0NegDauBinning confNegDauBinning;

  // setup lambdas
  v0builder::ConfLambdaSelection1 lambdaSelection;
  v0histmanager::ConfLambdaBinning1 confLambdaBinning;
  v0histmanager::V0HistManager<
    v0histmanager::PrefixLambda,
    trackhistmanager::PrefixV0PosDaughter,
    trackhistmanager::PrefixV0NegDaughter,
    modes::Mode::kAnalysis,
    modes::V0::kLambda>
    lambdaHistManager;
  Partition<Lambdas> lambdaPartition = MAKE_LAMBDA_PARTITION(lambdaSelection);
  Preslice<Lambdas> perColLambdas = aod::femtobase::stored::collisionId;

  // setup k0shorts
  v0builder::ConfK0shortSelection1 k0shortSelection;
  v0histmanager::ConfK0shortBinning1 confK0shortBinning;
  v0histmanager::V0HistManager<
    v0histmanager::PrefixK0short,
    trackhistmanager::PrefixV0PosDaughter,
    trackhistmanager::PrefixV0NegDaughter,
    modes::Mode::kAnalysis,
    modes::V0::kK0short>
    k0shortHistManager;
  Partition<K0shorts> k0shortPartition = MAKE_K0SHORT_PARTITION(k0shortSelection);
  Preslice<K0shorts> perColk0shorts = aod::femtobase::stored::collisionId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::PairHistManager<pairhistmanager::PrefixTrackV0Se, modes::Mode::kAnalysis> pairHistManagerSe;
  pairhistmanager::PairHistManager<pairhistmanager::PrefixTrackV0Me, modes::Mode::kAnalysis> pairHistManagerMe;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  ColumnBinningPolicy<femtocollisions::PosZ, femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}, true};
  ColumnBinningPolicy<aod::femtocollisions::PosZ, aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}, true};
  ColumnBinningPolicy<aod::femtocollisions::PosZ, aod::femtocollisions::Mult, aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}, true};
  pairhistmanager::ConfMixing confMixing;

  HistogramRegistry hRegistry{"FemtoTrackTrack", {}, OutputObjHandlingPolicy::AnalysisObject};

  // setup cpr
  closepairrejection::ConfCpr confCpr;
  closepairrejection::ClosePairRejectionTrackV0<closepairrejection::PrefixTrackPosDauSe, closepairrejection::PrefixTrackNegDauSe> cprSe;
  closepairrejection::ClosePairRejectionTrackV0<closepairrejection::PrefixTrackPosDauMe, closepairrejection::PrefixTrackNegDauMe> cprMe;
  paircleaner::TrackV0PairCleaner pc;

  void init(InitContext&)
  {

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins, confMixing.multBins.value}, true};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    // setup histograms for tracks
    auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
    colHistManager.init(&hRegistry, colHistSpec);

    auto trackHistSpec1 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning);
    trackHistManager.init(&hRegistry, trackHistSpec1);

    // setup for daughters
    auto posDauSpec = trackhistmanager::makeTrackHistSpecMap(confPosDauBinning);
    auto negDauSpec = trackhistmanager::makeTrackHistSpecMap(confNegDauBinning);

    // setup for lambda
    if (doprocessLambdaSameEvent || doprocessLambdaMixedEvent) {
      auto lambdaHistSpec = v0histmanager::makeV0HistSpecMap(confLambdaBinning);
      lambdaHistManager.init(&hRegistry, lambdaHistSpec, posDauSpec, negDauSpec);
    }

    if (((doprocessLambdaSameEvent || doprocessLambdaMixedEvent) + (doprocessK0shortSameEvent || doprocessK0shortMixedEvent)) > 1) {
      LOG(fatal) << "Can only process lambda-tracks Or k0short-tracks";
    }

    // setup for k0short
    if (doprocessK0shortSameEvent || doprocessK0shortMixedEvent) {
      auto k0shortHistSpec = v0histmanager::makeV0HistSpecMap(confK0shortBinning);
      k0shortHistManager.init(&hRegistry, k0shortHistSpec, posDauSpec, negDauSpec);
    }

    // setup histograms for pair
    auto pairHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confTrackBinning, confLambdaBinning);
    pairHistManagerSe.init(&hRegistry, pairHistSpec);
    pairHistManagerSe.setMass(trackSelections.pdgCode.value, lambdaSelection.pdgCode.value);
    pairHistManagerMe.init(&hRegistry, pairHistSpec);
    pairHistManagerMe.setMass(trackSelections.pdgCode.value, lambdaSelection.pdgCode.value);

    // setup histograms for cpr
    auto cprHistSpec = closepairrejection::makeCprHistSpecMap(confCpr);
    cprSe.init(&hRegistry, cprHistSpec, confCpr.detaMax.value, confCpr.dphistarMax.value, confCpr.on.value);
    cprMe.init(&hRegistry, cprHistSpec, confCpr.detaMax.value, confCpr.dphistarMax.value, confCpr.on.value);
  };

  void processLambdaSameEvent(FilteredCollision const& col, Tracks const& tracks, Lambdas const& /*lambdas*/)
  {
    auto trackSlice = trackPartition->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
    auto lambdaSlice = lambdaPartition->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
    if (trackSlice.size() == 0 || lambdaSlice.size() == 0) {
      return;
    }
    colHistManager.fill(col);
    cprSe.setMagField(col.magField());
    pairprocesshelpers::processSameEvent(trackSlice, lambdaSlice, tracks, trackHistManager, lambdaHistManager, pairHistManagerSe, cprSe, pc);
  }
  PROCESS_SWITCH(FemtounitedPairTrackV0, processLambdaSameEvent, "Enable processing same event processing for tracks and lambdas", true);

  void processLambdaMixedEvent(FilteredCollisions const& cols, Tracks const& tracks)
  {
    switch (confMixing.policy.value) {
      case static_cast<int>(pairhistmanager::kVtxMult):
        pairprocesshelpers::processMixedEvent(cols, trackPartition, lambdaPartition, tracks, cache, mixBinsVtxMult, confMixing.depth.value, trackHistManager, lambdaHistManager, pairHistManagerMe, cprMe, pc);
        break;
      case static_cast<int>(pairhistmanager::kVtxCent):
        pairprocesshelpers::processMixedEvent(cols, trackPartition, lambdaPartition, tracks, cache, mixBinsVtxCent, confMixing.depth.value, trackHistManager, lambdaHistManager, pairHistManagerMe, cprMe, pc);
        break;
      case static_cast<int>(pairhistmanager::kVtxMultCent):
        pairprocesshelpers::processMixedEvent(cols, trackPartition, lambdaPartition, tracks, cache, mixBinsVtxMultCent, confMixing.depth.value, trackHistManager, lambdaHistManager, pairHistManagerMe, cprMe, pc);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(FemtounitedPairTrackV0, processLambdaMixedEvent, "Enable processing mixed event processing for tracks and lambdas", true);

  void processK0shortSameEvent(FilteredCollision const& col, Tracks const& tracks, K0shorts const& /*k0shorts*/)
  {
    auto trackSlice = trackPartition->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
    auto k0shortSlice = k0shortPartition->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
    if (trackSlice.size() == 0 || k0shortSlice.size() == 0) {
      return;
    }
    colHistManager.fill(col);
    cprSe.setMagField(col.magField());
    pairprocesshelpers::processSameEvent(trackSlice, k0shortSlice, tracks, trackHistManager, k0shortHistManager, pairHistManagerSe, cprSe, pc);
  }
  PROCESS_SWITCH(FemtounitedPairTrackV0, processK0shortSameEvent, "Enable processing same event processing for tracks and k0shorts", false);

  void processK0shortMixedEvent(FilteredCollisions const& cols, Tracks const& tracks)
  {
    switch (confMixing.policy.value) {
      case static_cast<int>(pairhistmanager::kVtxMult):
        pairprocesshelpers::processMixedEvent(cols, trackPartition, k0shortPartition, tracks, cache, mixBinsVtxMult, confMixing.depth.value, trackHistManager, k0shortHistManager, pairHistManagerMe, cprMe, pc);
        break;
      case static_cast<int>(pairhistmanager::kVtxCent):
        pairprocesshelpers::processMixedEvent(cols, trackPartition, k0shortPartition, tracks, cache, mixBinsVtxCent, confMixing.depth.value, trackHistManager, k0shortHistManager, pairHistManagerMe, cprMe, pc);
        break;
      case static_cast<int>(pairhistmanager::kVtxMultCent):
        pairprocesshelpers::processMixedEvent(cols, trackPartition, k0shortPartition, tracks, cache, mixBinsVtxMultCent, confMixing.depth.value, trackHistManager, k0shortHistManager, pairHistManagerMe, cprMe, pc);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(FemtounitedPairTrackV0, processK0shortMixedEvent, "Enable processing mixed event processing for tracks and k0shorts", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtounitedPairTrackV0>(cfgc),
  };
  return workflow;
}
