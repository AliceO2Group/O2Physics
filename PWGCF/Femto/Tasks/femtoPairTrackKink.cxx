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

#include <map>
#include <vector>

using namespace o2::analysis::femto;

struct FemtoPairTrackKink {

  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, o2::aod::FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks>;
  using FemtoSigmas = o2::soa::Join<o2::aod::FSigmas, o2::aod::FSigmaMasks>;
  using FemtoSigmaPlus = o2::soa::Join<o2::aod::FSigmaPlus, o2::aod::FSigmaPlusMasks>;

  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, o2::aod::FTrackLabels>;
  using FemtoSigmasWithLabel = o2::soa::Join<FemtoSigmas, o2::aod::FSigmaLabels>;
  using FemtoSigmaPlusWithLabel = o2::soa::Join<FemtoSigmaPlus, o2::aod::FSigmaPlusLabels>;

  o2::framework::SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // setup tracks
  trackbuilder::ConfTrackSelection1 confTrackSelection;
  trackhistmanager::ConfTrackBinning1 confTrackBinning;
  particlecleaner::ConfTrackCleaner1 confTrackCleaner;

  o2::framework::Partition<FemtoTracks> trackPartition = MAKE_TRACK_PARTITION(confTrackSelection);
  o2::framework::Preslice<FemtoTracks> perColTracks = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoTracksWithLabel> trackWithLabelPartition = MAKE_TRACK_PARTITION(confTrackSelection);
  o2::framework::Preslice<FemtoTracksWithLabel> perColtracksWithLabel = o2::aod::femtobase::stored::fColId;

  // setup for daughters
  trackhistmanager::ConfKinkChaDauBinning confChaDauBinning;

  // setup sigmas
  kinkbuilder::ConfSigmaSelection1 confSigmaSelection;
  kinkhistmanager::ConfSigmaBinning1 confSigmaBinning;
  particlecleaner::ConfSigmaCleaner1 confSigmaCleaner;

  o2::framework::Partition<FemtoSigmas> sigmaPartition = MAKE_SIGMA_PARTITION(confSigmaSelection);
  o2::framework::Preslice<FemtoSigmas> perColSigmas = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoSigmasWithLabel> sigmaWithLabelPartition = MAKE_SIGMA_PARTITION(confSigmaSelection);
  o2::framework::Preslice<FemtoSigmasWithLabel> perColSigmasWithLabel = o2::aod::femtobase::stored::fColId;

  // setup for sigma plus
  kinkbuilder::ConfSigmaPlusSelection1 confSigmaPlusSelection;
  kinkhistmanager::ConfSigmaPlusBinning1 confSigmaPlusBinning;
  particlecleaner::ConfSigmaPlusCleaner1 confSigmaPlusCleaner;

  o2::framework::Partition<FemtoSigmaPlus> sigmaPlusPartition = MAKE_SIGMAPLUS_PARTITION(confSigmaPlusSelection);
  o2::framework::Preslice<FemtoSigmaPlus> perColSigmaPlus = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoSigmaPlusWithLabel> sigmaPlusWithLabelPartition = MAKE_SIGMAPLUS_PARTITION(confSigmaPlusSelection);
  o2::framework::Preslice<FemtoSigmaPlusWithLabel> perColSigmaPlusWithLabel = o2::aod::femtobase::stored::fColId;

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
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}, true};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}, true};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult, o2::aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}, true};
  pairhistmanager::ConfMixing confMixing;

  o2::framework::HistogramRegistry hRegistry{"FemtoTrackKink", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  // setup cpr
  closepairrejection::ConfCprTrackKinkDaughter confCpr;

  void init(o2::framework::InitContext&)
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

    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> trackHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> chaDauSpec;
    std::map<kinkhistmanager::KinkHist, std::vector<o2::framework::AxisSpec>> sigmaHistSpec;
    std::map<kinkhistmanager::KinkHist, std::vector<o2::framework::AxisSpec>> sigmaPlusHistSpec;
    std::map<pairhistmanager::PairHist, std::vector<o2::framework::AxisSpec>> pairTrackKinkHistSpec;
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpec = closepairrejection::makeCprHistSpecMap(confCpr);

    if (processData) {
      colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
      trackHistSpec = trackhistmanager::makeTrackHistSpecMap(confTrackBinning);
      chaDauSpec = trackhistmanager::makeTrackHistSpecMap(confChaDauBinning);
      pairTrackKinkHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning);
      if (processSigma) {
        sigmaHistSpec = kinkhistmanager::makeKinkHistSpecMap(confSigmaBinning);
        pairTrackSigmaBuilder.init<modes::Mode::kAnalysis>(&hRegistry, confCollisionBinning, confTrackSelection, confTrackCleaner, confSigmaSelection, confSigmaCleaner, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, sigmaHistSpec, chaDauSpec, pairTrackKinkHistSpec, cprHistSpec);
      } else {
        sigmaPlusHistSpec = kinkhistmanager::makeKinkHistSpecMap(confSigmaPlusBinning);
        pairTrackSigmaPlusBuilder.init<modes::Mode::kAnalysis>(&hRegistry, confCollisionBinning, confTrackSelection, confTrackCleaner, confSigmaPlusSelection, confSigmaPlusCleaner, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, sigmaPlusHistSpec, chaDauSpec, pairTrackKinkHistSpec, cprHistSpec);
      }
    } else {
      colHistSpec = colhistmanager::makeColMcHistSpecMap(confCollisionBinning);
      trackHistSpec = trackhistmanager::makeTrackMcHistSpecMap(confTrackBinning);
      chaDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confChaDauBinning);
      pairTrackKinkHistSpec = pairhistmanager::makePairMcHistSpecMap(confPairBinning);
      if (processSigma) {
        sigmaHistSpec = kinkhistmanager::makeKinkMcHistSpecMap(confSigmaBinning);
        pairTrackSigmaBuilder.init<modes::Mode::kAnalysis_Mc>(&hRegistry, confCollisionBinning, confTrackSelection, confTrackCleaner, confSigmaSelection, confSigmaCleaner, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, sigmaHistSpec, chaDauSpec, pairTrackKinkHistSpec, cprHistSpec);
      } else {
        sigmaPlusHistSpec = kinkhistmanager::makeKinkMcHistSpecMap(confSigmaPlusBinning);
        pairTrackSigmaPlusBuilder.init<modes::Mode::kAnalysis_Mc>(&hRegistry, confCollisionBinning, confTrackSelection, confTrackCleaner, confSigmaPlusSelection, confSigmaPlusCleaner, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, sigmaPlusHistSpec, chaDauSpec, pairTrackKinkHistSpec, cprHistSpec);
      }
    }
    hRegistry.print();
  };

  void processSigmaSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoSigmas const& sigmas)
  {
    pairTrackSigmaBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, trackPartition, sigmas, sigmaPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackKink, processSigmaSameEvent, "Enable processing same event processing for tracks and sigmas", true);

  void processSigmaSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoSigmasWithLabel const& sigmas, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackSigmaBuilder.processSameEvent<modes::Mode::kAnalysis_Mc>(col, mcCols, tracks, trackWithLabelPartition, sigmas, sigmaWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackKink, processSigmaSameEventMc, "Enable processing same event processing for tracks and sigmas with MC information", false);

  void processSigmaMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoSigmas const& /*sigmas*/)
  {
    pairTrackSigmaBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, trackPartition, sigmaPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackKink, processSigmaMixedEvent, "Enable processing mixed event processing for tracks and sigmas", true);

  void processSigmaMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoSigmasWithLabel const& /*sigmas*/, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackSigmaBuilder.processMixedEvent<modes::Mode::kAnalysis_Mc>(cols, mcCols, tracks, trackWithLabelPartition, sigmaWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackKink, processSigmaMixedEventMc, "Enable processing mixed event processing for tracks and sigmas with MC information", false);

  void processSigmaPlusSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoSigmaPlus const& sigmaplus)
  {
    pairTrackSigmaPlusBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, trackPartition, sigmaplus, sigmaPlusPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackKink, processSigmaPlusSameEvent, "Enable processing same event processing for tracks and sigma plus", false);

  void processSigmaPlusSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoSigmaPlusWithLabel const& sigmaplus, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackSigmaPlusBuilder.processSameEvent<modes::Mode::kAnalysis_Mc>(col, mcCols, tracks, trackWithLabelPartition, sigmaplus, sigmaPlusWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackKink, processSigmaPlusSameEventMc, "Enable processing same event processing for tracks and sigma plus with MC information", false);

  void processSigmaPlusMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoSigmaPlus const& /*sigmaplus*/)
  {
    pairTrackSigmaPlusBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, trackPartition, sigmaPlusPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackKink, processSigmaPlusMixedEvent, "Enable processing mixed event processing for tracks and sigma plus", false);

  void processSigmaPlusMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoSigmaPlusWithLabel const& /*sigmaplus*/, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackSigmaPlusBuilder.processMixedEvent<modes::Mode::kAnalysis_Mc>(cols, mcCols, tracks, trackWithLabelPartition, sigmaPlusWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackKink, processSigmaPlusMixedEventMc, "Enable processing mixed event processing for tracks and sigma plus with MC information", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairTrackKink>(cfgc),
  };
  return workflow;
}
