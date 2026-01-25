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
  using FemtoCollisions = Join<FCols, FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoTracks = o2::soa::Join<FTracks, FTrackMasks>;
  using FemtoSigmas = o2::soa::Join<FSigmas, FSigmaMasks>;
  using FemtoSigmaPlus = o2::soa::Join<FSigmaPlus, FSigmaPlusMasks>;

  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, FTrackLabels>;
  using FemtoSigmasWithLabel = o2::soa::Join<FemtoSigmas, FSigmaLabels>;
  using FemtoSigmaPlusWithLabel = o2::soa::Join<FemtoSigmaPlus, FSigmaPlusLabels>;

  SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // setup tracks
  trackbuilder::ConfTrackSelection1 confTrackSelection;
  trackhistmanager::ConfTrackBinning1 confTrackBinning;
  particlecleaner::ConfTrackCleaner1 confTrackCleaner;

  Partition<FemtoTracks> trackPartition = MAKE_TRACK_PARTITION(confTrackSelection);
  Preslice<FemtoTracks> perColTracks = aod::femtobase::stored::fColId;

  Partition<FemtoTracksWithLabel> trackWithLabelPartition = MAKE_TRACK_PARTITION(confTrackSelection);
  Preslice<FemtoTracksWithLabel> perColtracksWithLabel = aod::femtobase::stored::fColId;

  // setup for daughters
  trackhistmanager::ConfKinkChaDauBinning confChaDauBinning;

  // setup sigmas
  kinkbuilder::ConfSigmaSelection1 confSigmaSelection;
  kinkhistmanager::ConfSigmaBinning1 confSigmaBinning;
  particlecleaner::ConfSigmaCleaner1 confSigmaCleaner;

  Partition<FemtoSigmas> sigmaPartition = MAKE_SIGMA_PARTITION(confSigmaSelection);
  Preslice<FemtoSigmas> perColSigmas = aod::femtobase::stored::fColId;

  Partition<FemtoSigmasWithLabel> sigmaWithLabelPartition = MAKE_SIGMA_PARTITION(confSigmaSelection);
  Preslice<FemtoSigmasWithLabel> perColSigmasWithLabel = aod::femtobase::stored::fColId;

  // setup for sigma plus
  kinkbuilder::ConfSigmaPlusSelection1 confSigmaPlusSelection;
  kinkhistmanager::ConfSigmaPlusBinning1 confSigmaPlusBinning;
  particlecleaner::ConfSigmaPlusCleaner1 confSigmaPlusCleaner;

  Partition<FemtoSigmaPlus> sigmaPlusPartition = MAKE_SIGMAPLUS_PARTITION(confSigmaPlusSelection);
  Preslice<FemtoSigmaPlus> perColSigmaPlus = aod::femtobase::stored::fColId;

  Partition<FemtoSigmaPlusWithLabel> sigmaPlusWithLabelPartition = MAKE_SIGMAPLUS_PARTITION(confSigmaPlusSelection);
  Preslice<FemtoSigmaPlusWithLabel> perColSigmaPlusWithLabel = aod::femtobase::stored::fColId;

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
    bool processData = doprocessSigmaSameEvent || doprocessSigmaMixedEvent || doprocessSigmaPlusSameEvent || doprocessSigmaPlusMixedEvent;
    bool processMc = doprocessSigmaSameEventMc || doprocessSigmaMixedEventMc || doprocessSigmaPlusSameEventMc || doprocessSigmaPlusMixedEventMc;

    if (processData && processMc) {
      LOG(fatal) << "Both data and mc processing is enabled. Breaking...";
    }

    bool processSigma = doprocessSigmaSameEvent || doprocessSigmaSameEventMc || doprocessSigmaMixedEvent || doprocessSigmaMixedEventMc;
    bool processSigmaPlus = doprocessSigmaPlusSameEvent || doprocessSigmaPlusSameEventMc || doprocessSigmaPlusMixedEvent || doprocessSigmaPlusMixedEventMc;

    if (processSigma && processSigmaPlus) {
      LOG(fatal) << "Both Sigma-track and SigmaPlus-track processing is enabled. Breaking...";
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
      auto chaDauSpec = trackhistmanager::makeTrackHistSpecMap(confChaDauBinning);
      auto pairTrackKinkHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning);
      if (processSigma) {
        auto sigmaHistSpec = kinkhistmanager::makeKinkHistSpecMap(confSigmaBinning);
        pairTrackSigmaBuilder.init<modes::Mode::kAnalysis>(&hRegistry, confTrackSelection, confTrackCleaner, confSigmaSelection, confSigmaCleaner, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, sigmaHistSpec, chaDauSpec, pairTrackKinkHistSpec, cprHistSpec);
      } else {
        auto sigmaplusHistSpec = kinkhistmanager::makeKinkHistSpecMap(confSigmaPlusBinning);
        pairTrackSigmaPlusBuilder.init<modes::Mode::kAnalysis>(&hRegistry, confTrackSelection, confTrackCleaner, confSigmaPlusSelection, confSigmaPlusCleaner, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, sigmaplusHistSpec, chaDauSpec, pairTrackKinkHistSpec, cprHistSpec);
      }
    } else {
      auto colHistSpec = colhistmanager::makeColMcHistSpecMap(confCollisionBinning);
      auto trackHistSpec = trackhistmanager::makeTrackMcHistSpecMap(confTrackBinning);
      auto chaDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confChaDauBinning);
      auto pairTrackKinkHistSpec = pairhistmanager::makePairMcHistSpecMap(confPairBinning);
      if (processSigma) {
        auto sigmaHistSpec = kinkhistmanager::makeKinkMcHistSpecMap(confSigmaBinning);
        pairTrackSigmaBuilder.init<modes::Mode::kAnalysis_Mc>(&hRegistry, confTrackSelection, confTrackCleaner, confSigmaSelection, confSigmaCleaner, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, sigmaHistSpec, chaDauSpec, pairTrackKinkHistSpec, cprHistSpec);
      } else {
        auto sigmaplusHistSpec = kinkhistmanager::makeKinkMcHistSpecMap(confSigmaPlusBinning);
        pairTrackSigmaPlusBuilder.init<modes::Mode::kAnalysis_Mc>(&hRegistry, confTrackSelection, confTrackCleaner, confSigmaPlusSelection, confSigmaPlusCleaner, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, sigmaplusHistSpec, chaDauSpec, pairTrackKinkHistSpec, cprHistSpec);
      }
    }
    hRegistry.print();
  };

  void processSigmaSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoSigmas const& sigmas)
  {
    pairTrackSigmaBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, trackPartition, sigmas, sigmaPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackKink, processSigmaSameEvent, "Enable processing same event processing for tracks and sigmas", true);

  void processSigmaSameEventMc(FilteredFemtoCollisionWithLabel const& col, FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoSigmasWithLabel const& sigmas, FMcParticles const& mcParticles, FMcMothers const& mcMothers, FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackSigmaBuilder.processSameEvent<modes::Mode::kAnalysis_Mc>(col, mcCols, tracks, trackWithLabelPartition, sigmas, sigmaWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackKink, processSigmaSameEventMc, "Enable processing same event processing for tracks and sigmas with MC information", false);

  void processSigmaMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoSigmas const& /*sigmas*/)
  {
    pairTrackSigmaBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, trackPartition, sigmaPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackKink, processSigmaMixedEvent, "Enable processing mixed event processing for tracks and sigmas", true);

  void processSigmaMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoSigmasWithLabel const& /*sigmas*/, FMcParticles const& mcParticles, FMcMothers const& mcMothers, FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackSigmaBuilder.processMixedEvent<modes::Mode::kAnalysis_Mc>(cols, mcCols, tracks, trackWithLabelPartition, sigmaWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackKink, processSigmaMixedEventMc, "Enable processing mixed event processing for tracks and sigmas with MC information", false);

  void processSigmaPlusSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoSigmaPlus const& sigmaplus)
  {
    pairTrackSigmaPlusBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, trackPartition, sigmaplus, sigmaPlusPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackKink, processSigmaPlusSameEvent, "Enable processing same event processing for tracks and sigma plus", false);

  void processSigmaPlusSameEventMc(FilteredFemtoCollisionWithLabel const& col, FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoSigmaPlusWithLabel const& sigmaplus, FMcParticles const& mcParticles, FMcMothers const& mcMothers, FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackSigmaPlusBuilder.processSameEvent<modes::Mode::kAnalysis_Mc>(col, mcCols, tracks, trackWithLabelPartition, sigmaplus, sigmaPlusWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackKink, processSigmaPlusSameEventMc, "Enable processing same event processing for tracks and sigma plus with MC information", false);

  void processSigmaPlusMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoSigmaPlus const& /*sigmaplus*/)
  {
    pairTrackSigmaPlusBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, trackPartition, sigmaPlusPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackKink, processSigmaPlusMixedEvent, "Enable processing mixed event processing for tracks and sigma plus", false);

  void processSigmaPlusMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoSigmaPlusWithLabel const& /*sigmaplus*/, FMcParticles const& mcParticles, FMcMothers const& mcMothers, FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackSigmaPlusBuilder.processMixedEvent<modes::Mode::kAnalysis_Mc>(cols, mcCols, tracks, trackWithLabelPartition, sigmaPlusWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackKink, processSigmaPlusMixedEventMc, "Enable processing mixed event processing for tracks and sigma plus with MC information", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairTrackKink>(cfgc),
  };
  return workflow;
}
