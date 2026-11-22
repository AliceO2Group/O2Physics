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

/// \file femtoPairTrackSigma0.cxx
/// \brief Tasks that computes correlation between a track and a sigma0
/// \author Yakiv Paroviak, CTU in Prague, yakiv.paroviak@cern.ch

#include "PWGCF/Femto/Core/closePairRejection.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/pairBuilder.h"
#include "PWGCF/Femto/Core/pairCleaner.h"
#include "PWGCF/Femto/Core/pairHistManager.h"
#include "PWGCF/Femto/Core/particleCleaner.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/sigma0Builder.h"
#include "PWGCF/Femto/Core/sigma0HistManager.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <map>
#include <vector>

using namespace o2::analysis::femto;

struct FemtoPairTrackSigma0 {

  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, o2::aod::FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks>;
  using FemtoSigma0s = o2::soa::Join<o2::aod::FSigma0s, o2::aod::FSigma0Masks>;

  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, o2::aod::FTrackLabels>;
  using FemtoSigma0sWithLabel = o2::soa::Join<FemtoSigma0s, o2::aod::FSigma0Labels>;

  using FemtoMcParticlesWithLabel = o2::soa::Join<o2::aod::FMcParticles, o2::aod::FMcMotherLabels>;

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
  o2::framework::Preslice<FemtoTracksWithLabel> perColTracksWithLabel = o2::aod::femtobase::stored::fColId;

  // setup for the Lambda daughter tracks of the Sigma0
  trackhistmanager::ConfSigma0PosDauBinning confPosDauBinning;
  trackhistmanager::ConfSigma0NegDauBinning confNegDauBinning;

  // setup sigma0s
  sigma0builder::ConfSigma0SelectionDefault confSigma0Selection;
  sigma0histmanager::ConfSigma0BinningDefault confSigma0Binning;
  particlecleaner::ConfSigma0Cleaner confSigma0Cleaner;

  o2::framework::Partition<FemtoSigma0s> sigma0Partition = MAKE_SIGMA0_PARTITION(confSigma0Selection);
  o2::framework::Preslice<FemtoSigma0s> perColSigma0s = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoSigma0sWithLabel> sigma0WithLabelPartition = MAKE_SIGMA0_PARTITION(confSigma0Selection);
  o2::framework::Preslice<FemtoSigma0sWithLabel> perColSigma0sWithLabel = o2::aod::femtobase::stored::fColId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::ConfPairCuts confPairCuts;
  paircleaner::ConfPairCleanerBinning confPairCleaner;

  pairbuilder::PairTrackSigma0Builder<
    trackhistmanager::PrefixTrack1Se,
    sigma0histmanager::PrefixSigma0Se,
    trackhistmanager::PrefixSigma0PosDaughterSe,
    trackhistmanager::PrefixSigma0NegDaughterSe,
    trackhistmanager::PrefixTrack1Me,
    sigma0histmanager::PrefixSigma0Me,
    trackhistmanager::PrefixSigma0PosDaughterMe,
    trackhistmanager::PrefixSigma0NegDaughterMe,
    pairhistmanager::PrefixTrackSigma0Se,
    pairhistmanager::PrefixTrackSigma0Me,
    closepairrejection::PrefixTrackSigma0DaughterSe,
    closepairrejection::PrefixTrackSigma0DaughterMe>
    pairTrackSigma0Builder;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult, o2::aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}};
  pairhistmanager::ConfMixing confMixing;

  o2::framework::HistogramRegistry hRegistry{"FemtoTrackSigma0", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  // setup cpr (track vs the same-charge Lambda daughter of the Sigma0; photon legs are not in the track table)
  closepairrejection::ConfCprTrackSigma0Daughter confCpr;

  void init(o2::framework::InitContext&)
  {
    bool processData = doprocessSameEvent || doprocessMixedEvent;
    bool processMc = doprocessSameEventMc || doprocessMixedEventMc;
    if (processData && processMc) {
      LOG(fatal) << "Both data and mc processing is activated. Breaking...";
    }

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins.value, confMixing.multBins.value}};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}};

    // setup histogram specs
    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> trackHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> posDauSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> negDauSpec;
    std::map<sigma0histmanager::Sigma0Hist, std::vector<o2::framework::AxisSpec>> sigma0HistSpec;
    std::map<pairhistmanager::PairHist, std::vector<o2::framework::AxisSpec>> pairTrackSigma0HistSpec;
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpec = closepairrejection::makeCprHistSpecMap(confCpr);
    std::map<paircleaner::PairCleanerHist, std::vector<o2::framework::AxisSpec>> pairCleanerHistSpec = paircleaner::makePairCleanerHistSpecMap(confPairCleaner);

    if (processData) {
      colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
      trackHistSpec = trackhistmanager::makeTrackHistSpecMap(confTrackBinning);
      posDauSpec = trackhistmanager::makeTrackHistSpecMap(confPosDauBinning);
      negDauSpec = trackhistmanager::makeTrackHistSpecMap(confNegDauBinning);
      sigma0HistSpec = sigma0histmanager::makeSigma0HistSpecMap(confSigma0Binning);
      pairTrackSigma0HistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);
      pairTrackSigma0Builder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, confTrackSelection, confTrackCleaner, confSigma0Selection, confSigma0Cleaner, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, sigma0HistSpec, posDauSpec, negDauSpec, pairTrackSigma0HistSpec, cprHistSpec, pairCleanerHistSpec);
    } else if (processMc) {
      colHistSpec = colhistmanager::makeColMcHistSpecMap(confCollisionBinning);
      trackHistSpec = trackhistmanager::makeTrackMcHistSpecMap(confTrackBinning);
      posDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confPosDauBinning);
      negDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confNegDauBinning);
      sigma0HistSpec = sigma0histmanager::makeSigma0McHistSpecMap(confSigma0Binning);
      pairTrackSigma0HistSpec = pairhistmanager::makePairMcHistSpecMap(confPairBinning, confMixing);
      pairTrackSigma0Builder.init<modes::Mode::kSe_Reco_Mc, modes::Mode::kMe_Reco_Mc>(&hRegistry, confCollisionBinning, confTrackSelection, confTrackCleaner, confSigma0Selection, confSigma0Cleaner, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, sigma0HistSpec, posDauSpec, negDauSpec, pairTrackSigma0HistSpec, cprHistSpec, pairCleanerHistSpec);
    }
    hRegistry.print();
  };

  void processSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoSigma0s const& sigma0s)
  {
    pairTrackSigma0Builder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, trackPartition, sigma0s, sigma0Partition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackSigma0, processSameEvent, "Enable processing same event processing for tracks and sigma0s", true);

  void processSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoSigma0sWithLabel const& sigma0s, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackSigma0Builder.processSameEvent<modes::Mode::kSe_Reco_Mc>(col, mcCols, tracks, trackWithLabelPartition, sigma0s, sigma0WithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackSigma0, processSameEventMc, "Enable processing same event processing for tracks and sigma0s with mc information", false);

  void processMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoSigma0s const& /*sigma0s*/)
  {
    pairTrackSigma0Builder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, trackPartition, sigma0Partition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackSigma0, processMixedEvent, "Enable processing mixed event processing for tracks and sigma0s", true);

  void processMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoSigma0sWithLabel const& /*sigma0s*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackSigma0Builder.processMixedEvent<modes::Mode::kMe_Reco_Mc>(cols, mcCols, tracks, trackWithLabelPartition, sigma0WithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackSigma0, processMixedEventMc, "Enable processing mixed event processing for tracks and sigma0s with mc information", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairTrackSigma0>(context),
  };
  return workflow;
}
