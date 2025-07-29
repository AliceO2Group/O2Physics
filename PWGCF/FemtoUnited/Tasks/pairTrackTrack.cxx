// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file pairTrackTrack.cxx
/// \brief Tasks that computes correlation between two tracks
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/FemtoUnited/Core/closePairRejection.h"
#include "PWGCF/FemtoUnited/Core/collisionHistManager.h"
#include "PWGCF/FemtoUnited/Core/collisionSelection.h"
#include "PWGCF/FemtoUnited/Core/modes.h"
#include "PWGCF/FemtoUnited/Core/pairCleaner.h"
#include "PWGCF/FemtoUnited/Core/pairHistManager.h"
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
#include <random>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtounited;

struct PairTrackTrack {

  struct : ConfigurableGroup {
    std::string prefix = std::string("Options");
    Configurable<bool> correlatedPlots{"correlatedPlots", false, "Enable multidimensional histogramms. High memory consumption."};
  } Options;

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
  colhistmanager::CollisionHistManager colHistManager;

  // setup tracks
  trackselection::ConfTrackSelection1 trackSelections1;
  trackhistmanager::ConfTrackBinning1 confTrackBinning1;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixTrack1> trackHistManager1;
  Partition<Tracks> trackPartition1 = MAKE_TRACK_PARTITION(trackSelections1);

  trackselection::ConfTrackSelection2 trackSelections2;
  trackhistmanager::ConfTrackBinning2 confTrackBinning2;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixTrack2> trackHistManager2;
  Partition<Tracks> trackPartition2 = MAKE_TRACK_PARTITION(trackSelections1);

  Preslice<Tracks> perColReco = aod::femtobase::stored::collisionId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::PairHistManager<pairhistmanager::PrefixTrackTrackSe> pairHistManagerSe;
  pairhistmanager::PairHistManager<pairhistmanager::PrefixTrackTrackMe> pairHistManagerMe;

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
  closepairrejection::ClosePairRejection<closepairrejection::PrefixTrackTrackSe> cprSe;
  closepairrejection::ClosePairRejection<closepairrejection::PrefixTrackTrackMe> cprMe;
  paircleaner::PairCleaner pc;

  void init(InitContext&)
  {

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins, confMixing.multBins.value}, true};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    // setup histograms for tracks
    auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
    colHistManager.init<modes::Mode::kANALYSIS>(&hRegistry, colHistSpec);

    auto trackHistSpec1 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning1);
    trackHistManager1.init<modes::Mode::kANALYSIS>(&hRegistry, trackHistSpec1);

    auto trackHistSpec2 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning2);
    trackHistManager2.init<modes::Mode::kANALYSIS>(&hRegistry, trackHistSpec2);

    // setup histograms for pair
    auto pairHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confTrackBinning1, confTrackBinning2);
    pairHistManagerSe.init<modes::Mode::kANALYSIS>(&hRegistry, pairHistSpec);
    pairHistManagerSe.setMass(trackSelections1.pdgCode.value, trackSelections2.pdgCode.value);
    pairHistManagerMe.init<modes::Mode::kANALYSIS>(&hRegistry, pairHistSpec);
    pairHistManagerMe.setMass(trackSelections1.pdgCode.value, trackSelections2.pdgCode.value);

    // setup histograms for cpr
    auto cprHistSpec = closepairrejection::makeCprHistSpecMap(confCpr);
    cprSe.init<modes::Mode::kANALYSIS>(&hRegistry, cprHistSpec);
    cprSe.activate(confCpr.on.value);
    cprSe.setLimits(confCpr.detaMax.value, confCpr.dphistarMax.value);
    cprMe.init<modes::Mode::kANALYSIS>(&hRegistry, cprHistSpec);
    cprMe.activate(confCpr.on.value);
    cprMe.setLimits(confCpr.detaMax, confCpr.dphistarMax);

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
      colHistManager.fill<modes::Mode::kANALYSIS>(col);
      cprSe.setMagField(col.magField());
      pairhistmanager::processSameEvent<modes::Mode::kANALYSIS>(trackSlice1, trackHistManager1, pairHistManagerSe, cprSe, rng, mixParticles);
    } else {
      auto trackSlice1 = trackPartition1->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
      auto trackSlice2 = trackPartition2->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
      if (trackSlice1.size() == 0 || trackSlice2.size() == 0) {
        return;
      }
      colHistManager.fill<modes::Mode::kANALYSIS>(col);
      cprSe.setMagField(col.magField());
      pairhistmanager::processSameEvent<modes::Mode::kANALYSIS>(trackSlice1, trackSlice2, trackHistManager1, trackHistManager2, pairHistManagerSe, cprSe, pc);
    }
  }
  PROCESS_SWITCH(PairTrackTrack, processSameEvent, "Enable processing same event processing", true);

  void processMixedEvent(FilteredCollisions const& cols, Tracks const& /*tracks*/)
  {
    switch (confMixing.policy.value) {
      case static_cast<int>(pairhistmanager::kVtxMult):
        pairhistmanager::processMixedEvent<modes::Mode::kANALYSIS>(cols, trackPartition1, trackPartition2, cache, mixBinsVtxMult, confMixing.depth.value, trackHistManager1, trackHistManager2, pairHistManagerMe, cprMe, pc);
        break;
      case static_cast<int>(pairhistmanager::kVtxCent):
        pairhistmanager::processMixedEvent<modes::Mode::kANALYSIS>(cols, trackPartition1, trackPartition2, cache, mixBinsVtxCent, confMixing.depth.value, trackHistManager1, trackHistManager2, pairHistManagerMe, cprMe, pc);
        break;
      case static_cast<int>(pairhistmanager::kVtxMultCent):
        pairhistmanager::processMixedEvent<modes::Mode::kANALYSIS>(cols, trackPartition1, trackPartition2, cache, mixBinsVtxMultCent, confMixing.depth.value, trackHistManager1, trackHistManager2, pairHistManagerMe, cprMe, pc);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(PairTrackTrack, processMixedEvent, "Enable processing mixed event processing", true);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<PairTrackTrack>(cfgc),
  };
  return workflow;
}
