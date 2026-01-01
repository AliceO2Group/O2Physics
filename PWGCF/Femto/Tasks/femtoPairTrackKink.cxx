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

/// \file femtoPairTrackKink.cxx
/// \brief Tasks that computes correlation between tracks and kinks
/// \author Anton Riedel, TU München, anton.riedel@cern.ch
/// \author Henrik Fribert, TU München, henrik.fribert@cern.ch

#include "PWGCF/Femto/Core/closePairRejection.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/kinkBuilder.h"
#include "PWGCF/Femto/Core/kinkHistManager.h"
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

struct FemtoPairTrackKink {

  // setup tables
  using Collisions = Join<FCols, FColMasks>;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  using Tracks = o2::soa::Join<FTracks, FTrackMasks>;
  using Sigmas = o2::soa::Join<FSigmas, FSigmaMasks>;
  using SigmaPlus = o2::soa::Join<FSigmaPlus, FSigmaPlusMasks>;

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
  trackhistmanager::ConfKinkChaDauBinning confChaDauBinning;

  // setup sigmas
  kinkbuilder::ConfSigmaSelection1 sigmaSelection;
  kinkhistmanager::ConfSigmaBinning1 confSigmaBinning;
  Partition<Sigmas> sigmaPartition = MAKE_SIGMA_PARTITION(sigmaSelection);
  Preslice<Sigmas> perColSigmas = aod::femtobase::stored::fColId;

  // setup for sigma plus
  kinkbuilder::ConfSigmaPlusSelection1 sigmaPlusSelection;
  kinkhistmanager::ConfSigmaPlusBinning1 confSigmaPlusBinning;
  Partition<SigmaPlus> sigmaPlusPartition = MAKE_SIGMAPLUS_PARTITION(sigmaPlusSelection);
  Preslice<SigmaPlus> perColSigmaPlus = aod::femtobase::stored::fColId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::ConfPairCuts confPairCuts;

  pairbuilder::PairTrackKinkBuilder<
    trackhistmanager::PrefixTrack1,
    kinkhistmanager::PrefixSigma1,
    trackhistmanager::PrefixKinkChaDaughter,
    pairhistmanager::PrefixTrackKinkSe,
    pairhistmanager::PrefixTrackKinkMe,
    closepairrejection::PrefixTrackKinkSe,
    closepairrejection::PrefixTrackKinkMe,
    modes::Kink::kSigma>
    pairTrackSigmaBuilder;

  pairbuilder::PairTrackKinkBuilder<
    trackhistmanager::PrefixTrack1,
    kinkhistmanager::PrefixSigmaPlus1,
    trackhistmanager::PrefixKinkChaDaughter,
    pairhistmanager::PrefixTrackKinkSe,
    pairhistmanager::PrefixTrackKinkMe,
    closepairrejection::PrefixTrackKinkSe,
    closepairrejection::PrefixTrackKinkMe,
    modes::Kink::kSigmaPlus>
    pairTrackSigmaPlusBuilder;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  ColumnBinningPolicy<femtocollisions::PosZ, femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}, true};
  ColumnBinningPolicy<aod::femtocollisions::PosZ, aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}, true};
  ColumnBinningPolicy<aod::femtocollisions::PosZ, aod::femtocollisions::Mult, aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}, true};
  pairhistmanager::ConfMixing confMixing;

  HistogramRegistry hRegistry{"FemtoTrackKink", {}, OutputObjHandlingPolicy::AnalysisObject};

  // setup cpr
  closepairrejection::ConfCprTrackKinkDaughter confCpr;

  void init(InitContext&)
  {
    if (((doprocessSigmaSameEvent || doprocessSigmaMixedEvent) + (doprocessSigmaPlusSameEvent || doprocessSigmaPlusMixedEvent)) > 1) {
      LOG(fatal) << "Can only process sigma-tracks Or sigmaplus-tracks";
    }
    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins, confMixing.multBins.value}, true};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    // setup histograms
    auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
    auto trackHistSpec = trackhistmanager::makeTrackHistSpecMap(confTrackBinning);
    auto chaDauSpec = trackhistmanager::makeTrackHistSpecMap(confChaDauBinning);
    auto pairHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning);
    auto cprHistSpec = closepairrejection::makeCprHistSpecMap(confCpr);

    // setup for sigma
    if (doprocessSigmaSameEvent || doprocessSigmaMixedEvent) {
      auto sigmaHistSpec = kinkhistmanager::makeKinkHistSpecMap(confSigmaBinning);
      auto pairTrackSigmaHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning);
      pairTrackSigmaBuilder.init<modes::Mode::kAnalysis>(&hRegistry, trackSelection, sigmaSelection, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, sigmaHistSpec, chaDauSpec, pairTrackSigmaHistSpec, cprHistSpec);
    }

    // setup for sigma plus
    if (doprocessSigmaPlusSameEvent || doprocessSigmaPlusMixedEvent) {
      auto sigmaplusHistSpec = kinkhistmanager::makeKinkHistSpecMap(confSigmaPlusBinning);
      auto pairTrackSigmaPlusHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning);
      pairTrackSigmaPlusBuilder.init<modes::Mode::kAnalysis>(&hRegistry, trackSelection, sigmaPlusSelection, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, sigmaplusHistSpec, chaDauSpec, pairTrackSigmaPlusHistSpec, cprHistSpec);
    }
    hRegistry.print();
  };

  void processSigmaSameEvent(FilteredCollision const& col, Tracks const& tracks, Sigmas const& sigmas)
  {
    pairTrackSigmaBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, trackPartition, sigmas, sigmaPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackKink, processSigmaSameEvent, "Enable processing same event processing for tracks and sigmas", true);

  void processSigmaMixedEvent(FilteredCollisions const& cols, Tracks const& tracks, Sigmas const& /*sigmas*/)
  {
    pairTrackSigmaBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, trackPartition, sigmaPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackKink, processSigmaMixedEvent, "Enable processing mixed event processing for tracks and sigmas", true);
  //
  void processSigmaPlusSameEvent(FilteredCollision const& col, Tracks const& tracks, SigmaPlus const& sigmaplus)
  {
    pairTrackSigmaPlusBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, trackPartition, sigmaplus, sigmaPlusPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackKink, processSigmaPlusSameEvent, "Enable processing same event processing for tracks and sigma plus", false);

  void processSigmaPlusMixedEvent(FilteredCollisions const& cols, Tracks const& tracks, SigmaPlus const& /*sigmaplus*/)
  {
    pairTrackSigmaPlusBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, trackPartition, sigmaPlusPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackKink, processSigmaPlusMixedEvent, "Enable processing mixed event processing for tracks and sigma plus", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairTrackKink>(cfgc),
  };
  return workflow;
}
