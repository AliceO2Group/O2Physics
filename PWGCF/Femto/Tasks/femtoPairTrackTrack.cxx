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
#include "PWGCF/Femto/Core/pairCleaner.h"
#include "PWGCF/Femto/Core/pairHistManager.h"
#include "PWGCF/Femto/Core/pairProcessHelpers.h"
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

struct FemtoPairTrackTrack {

  // setup tables
  using Collisions = FCols;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  using Tracks = o2::soa::Join<FTracks, FTrackMasks>;

  SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::CollisionHistManager<modes::Mode::kAnalysis> colHistManager;

  // setup tracks
  trackbuilder::ConfTrackSelection1 trackSelections1;
  trackhistmanager::ConfTrackBinning1 confTrackBinning1;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixTrack1, modes::Mode::kAnalysis> trackHistManager1;
  Partition<Tracks> trackPartition1 = MAKE_TRACK_PARTITION(trackSelections1);

  trackbuilder::ConfTrackSelection2 trackSelections2;
  trackhistmanager::ConfTrackBinning2 confTrackBinning2;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixTrack2, modes::Mode::kAnalysis> trackHistManager2;
  Partition<Tracks> trackPartition2 = MAKE_TRACK_PARTITION(trackSelections1);

  Preslice<Tracks> perColReco = aod::femtobase::stored::collisionId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::PairHistManager<pairhistmanager::PrefixTrackTrackSe, modes::Mode::kAnalysis> pairHistManagerSe;
  pairhistmanager::PairHistManager<pairhistmanager::PrefixTrackTrackMe, modes::Mode::kAnalysis> pairHistManagerMe;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  ColumnBinningPolicy<femtocollisions::PosZ, femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}, true};
  ColumnBinningPolicy<aod::femtocollisions::PosZ, aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}, true};
  ColumnBinningPolicy<aod::femtocollisions::PosZ, aod::femtocollisions::Mult, aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}, true};
  pairhistmanager::ConfMixing confMixing;

  HistogramRegistry hRegistry{"FemtoTrackTrack", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::mt19937 rng;
  bool mixParticles = false;

  // setup cpr
  closepairrejection::ConfCpr confCpr;
  closepairrejection::ClosePairRejectionTrackTrack<closepairrejection::PrefixTrackTrackSe> cprSe;
  closepairrejection::ClosePairRejectionTrackTrack<closepairrejection::PrefixTrackTrackMe> cprMe;
  paircleaner::TrackTrackPairCleaner pc;

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

    auto trackHistSpec1 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning1);
    trackHistManager1.init(&hRegistry, trackHistSpec1);

    auto trackHistSpec2 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning2);
    trackHistManager2.init(&hRegistry, trackHistSpec2);

    // setup histograms for pair
    auto pairHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confTrackBinning1, confTrackBinning2);
    pairHistManagerSe.init(&hRegistry, pairHistSpec);
    pairHistManagerSe.setMass(trackSelections1.pdgCode.value, trackSelections2.pdgCode.value);
    pairHistManagerMe.init(&hRegistry, pairHistSpec);
    pairHistManagerMe.setMass(trackSelections1.pdgCode.value, trackSelections2.pdgCode.value);

    // setup histograms for cpr
    auto cprHistSpec = closepairrejection::makeCprHistSpecMap(confCpr);
    cprSe.init(&hRegistry, cprHistSpec, confCpr.detaMax.value, confCpr.dphistarMax.value, confCpr.on.value);
    cprMe.init(&hRegistry, cprHistSpec, confCpr.detaMax.value, confCpr.dphistarMax.value, confCpr.on.value);

    // setup rng if necessary
    if (confMixing.seed.value >= 0) {
      uint64_t randomSeed = 0;
      mixParticles = true;
      if (confMixing.seed.value == 0) {
        randomSeed = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
      } else {
        randomSeed = static_cast<uint64_t>(confMixing.seed.value);
      }
      rng = std::mt19937(randomSeed);
    }
  };

  void processSameEvent(FilteredCollision const& col, Tracks const& /*tracks*/)
  {
    if (confMixing.sameSpecies) {
      auto trackSlice1 = trackPartition1->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
      if (trackSlice1.size() == 0) {
        return;
      }
      colHistManager.fill(col);
      cprSe.setMagField(col.magField());
      pairprocesshelpers::processSameEvent(trackSlice1, trackHistManager1, pairHistManagerSe, cprSe, rng, mixParticles);
    } else {
      auto trackSlice1 = trackPartition1->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
      auto trackSlice2 = trackPartition2->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
      if (trackSlice1.size() == 0 || trackSlice2.size() == 0) {
        return;
      }
      colHistManager.fill(col);
      cprSe.setMagField(col.magField());
      pairprocesshelpers::processSameEvent(trackSlice1, trackSlice2, trackHistManager1, trackHistManager2, pairHistManagerSe, cprSe, pc);
    }
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processSameEvent, "Enable processing same event processing", true);

  void processMixedEvent(FilteredCollisions const& cols, Tracks const& /*tracks*/)
  {
    if (confMixing.sameSpecies) {
      switch (confMixing.policy.value) {
        case static_cast<int>(pairhistmanager::kVtxMult):
          pairprocesshelpers::processMixedEvent(cols, trackPartition1, cache, mixBinsVtxMult, confMixing.depth.value, trackHistManager1, trackHistManager2, pairHistManagerMe, cprMe, pc);
          break;
        case static_cast<int>(pairhistmanager::kVtxCent):
          pairprocesshelpers::processMixedEvent(cols, trackPartition1, cache, mixBinsVtxCent, confMixing.depth.value, trackHistManager1, trackHistManager2, pairHistManagerMe, cprMe, pc);
          break;
        case static_cast<int>(pairhistmanager::kVtxMultCent):
          pairprocesshelpers::processMixedEvent(cols, trackPartition1, cache, mixBinsVtxMultCent, confMixing.depth.value, trackHistManager1, trackHistManager2, pairHistManagerMe, cprMe, pc);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    } else {
      switch (confMixing.policy.value) {
        case static_cast<int>(pairhistmanager::kVtxMult):
          pairprocesshelpers::processMixedEvent(cols, trackPartition1, trackPartition2, cache, mixBinsVtxMult, confMixing.depth.value, trackHistManager1, trackHistManager2, pairHistManagerMe, cprMe, pc);
          break;
        case static_cast<int>(pairhistmanager::kVtxCent):
          pairprocesshelpers::processMixedEvent(cols, trackPartition1, trackPartition2, cache, mixBinsVtxCent, confMixing.depth.value, trackHistManager1, trackHistManager2, pairHistManagerMe, cprMe, pc);
          break;
        case static_cast<int>(pairhistmanager::kVtxMultCent):
          pairprocesshelpers::processMixedEvent(cols, trackPartition1, trackPartition2, cache, mixBinsVtxMultCent, confMixing.depth.value, trackHistManager1, trackHistManager2, pairHistManagerMe, cprMe, pc);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    }
  }
  PROCESS_SWITCH(FemtoPairTrackTrack, processMixedEvent, "Enable processing mixed event processing", true);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairTrackTrack>(cfgc),
  };
  return workflow;
}
