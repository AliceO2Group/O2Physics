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

/// \file femtoPairTrackTwoTrackResonance.cxx
/// \brief Tasks that computes correlation between tracks and resonances decaying into two tracks
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
#include "PWGCF/Femto/Core/twoTrackResonanceBuilder.h"
#include "PWGCF/Femto/Core/twoTrackResonanceHistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <string>
#include <vector>

using namespace o2::analysis::femto;

struct FemtoPairTrackTwoTrackResonance {

  // setup tables
  using Collisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks>;
  using FemtoPhis = o2::soa::Join<o2::aod::FPhis, o2::aod::FPhiMasks>;
  using FemtoKstar0s = o2::soa::Join<o2::aod::FKstar0s, o2::aod::FKstar0Masks>;
  using FemtoRho0s = o2::soa::Join<o2::aod::FRho0s, o2::aod::FRho0Masks>;

  o2::framework::SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // setup tracks
  trackbuilder::ConfTrackSelection1 trackSelection;
  trackhistmanager::ConfTrackBinning1 confTrackBinning;
  o2::framework::Partition<FemtoTracks> trackPartition = MAKE_TRACK_PARTITION(trackSelection);
  o2::framework::Preslice<FemtoTracks> perColTracks = o2::aod::femtobase::stored::fColId;

  // setup for daughters
  trackhistmanager::ConfResonancePosDauBinning confPosDauBinning;
  trackhistmanager::ConfResonanceNegDauBinning confNegDauBinning;

  // setup phis
  twotrackresonancebuilder::ConfPhiSelection phiSelection;
  twotrackresonancehistmanager::ConfPhiBinning confPhiBinning;
  o2::framework::Partition<FemtoPhis> phiPartition = MAKE_RESONANCE_0_PARTITON(phiSelection);
  o2::framework::Preslice<FemtoPhis> perColPhis = o2::aod::femtobase::stored::fColId;

  // setup kstar0
  twotrackresonancebuilder::ConfKstar0Selection kstar0Selection;
  twotrackresonancehistmanager::ConfKstar0Binning confKstar0Binning;
  o2::framework::Partition<FemtoKstar0s> kstar0Partition = MAKE_RESONANCE_1_PARTITON(kstar0Selection);
  o2::framework::Preslice<FemtoKstar0s> perColKstar0s = o2::aod::femtobase::stored::fColId;

  // rho0s
  twotrackresonancebuilder::ConfRho0Selection rho0Selection;
  twotrackresonancehistmanager::ConfRho0Binning confRho0Binning;
  o2::framework::Partition<FemtoRho0s> rho0Partition = MAKE_RESONANCE_0_PARTITON(rho0Selection);
  o2::framework::Preslice<FemtoRho0s> perColRho0s = o2::aod::femtobase::stored::fColId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::ConfPairCuts confPairCuts;

  // setup for track-phi pairs
  pairbuilder::PairTrackTwoTrackResonanceBuilder<
    trackhistmanager::PrefixTrack1,
    twotrackresonancehistmanager::PrefixPhi,
    trackhistmanager::PrefixResonancePosDaughter,
    trackhistmanager::PrefixResonanceNegDaughter,
    pairhistmanager::PrefixTrackResonanceSe,
    pairhistmanager::PrefixTrackResonanceMe,
    closepairrejection::PrefixTrackTwoTrackResonanceSe,
    closepairrejection::PrefixTrackTwoTrackResonanceMe,
    modes::TwoTrackResonance::kPhi>
    pairTrackPhiBuilder;

  // setup for track-kstar0 pairs
  pairbuilder::PairTrackTwoTrackResonanceBuilder<
    trackhistmanager::PrefixTrack1,
    twotrackresonancehistmanager::PrefixKstar,
    trackhistmanager::PrefixResonancePosDaughter,
    trackhistmanager::PrefixResonanceNegDaughter,
    pairhistmanager::PrefixTrackResonanceSe,
    pairhistmanager::PrefixTrackResonanceMe,
    closepairrejection::PrefixTrackTwoTrackResonanceSe,
    closepairrejection::PrefixTrackTwoTrackResonanceMe,
    modes::TwoTrackResonance::kKstar0>
    pairTrackKstar0Builder;

  // setup for track-rho0 pairs
  pairbuilder::PairTrackTwoTrackResonanceBuilder<
    trackhistmanager::PrefixTrack1,
    twotrackresonancehistmanager::PrefixRho,
    trackhistmanager::PrefixResonancePosDaughter,
    trackhistmanager::PrefixResonanceNegDaughter,
    pairhistmanager::PrefixTrackResonanceSe,
    pairhistmanager::PrefixTrackResonanceMe,
    closepairrejection::PrefixTrackTwoTrackResonanceSe,
    closepairrejection::PrefixTrackTwoTrackResonanceMe,
    modes::TwoTrackResonance::kRho0>
    pairTrackRho0Builder;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}, true};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}, true};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult, o2::aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}, true};
  pairhistmanager::ConfMixing confMixing;

  o2::framework::HistogramRegistry hRegistry{"FemtoTrackTwoTrackResonance", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  // setup cpr
  closepairrejection::ConfCprTrackResonanceDaughter confCpr;

  void init(o2::framework::InitContext&)
  {

    if (((doprocessPhiSameEvent || doprocessPhiMixedEvent) + (doprocessKstar0SameEvent || doprocessKstar0MixedEvent)) + (doprocessRho0SameEvent || doprocessRho0MixedEvent) > 1) {
      LOG(fatal) << "Can only process phi-tracks, rho-tracks or k0*-tracks";
    }

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

    // setup for phi
    if (doprocessPhiSameEvent || doprocessPhiMixedEvent) {
      auto phiHistSpec = twotrackresonancehistmanager::makeTwoTrackResonanceHistSpecMap(confPhiBinning);
      auto pairTrackPhiHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);
      pairTrackPhiBuilder.init<modes::Mode::kSe_Analysis, modes::Mode::kMe_Analysis>(&hRegistry, confCollisionBinning, trackSelection, phiSelection, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, phiHistSpec, posDauSpec, negDauSpec, pairTrackPhiHistSpec, cprHistSpec);
    }

    // setup for kstar0
    if (doprocessKstar0SameEvent || doprocessKstar0MixedEvent) {
      auto kstar0HistSpec = twotrackresonancehistmanager::makeTwoTrackResonanceHistSpecMap(confKstar0Binning);
      auto pairTrackKstar0HistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);
      pairTrackKstar0Builder.init<modes::Mode::kSe_Analysis, modes::Mode::kMe_Analysis>(&hRegistry, confCollisionBinning, trackSelection, kstar0Selection, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, kstar0HistSpec, posDauSpec, negDauSpec, pairTrackKstar0HistSpec, cprHistSpec);
    }

    // setup for kstar0
    if (doprocessRho0SameEvent || doprocessRho0MixedEvent) {
      auto rho0HistSpec = twotrackresonancehistmanager::makeTwoTrackResonanceHistSpecMap(confRho0Binning);
      auto pairTrackRho0HistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);
      pairTrackRho0Builder.init<modes::Mode::kSe_Analysis, modes::Mode::kMe_Analysis>(&hRegistry, confCollisionBinning, trackSelection, rho0Selection, confCpr, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, rho0HistSpec, posDauSpec, negDauSpec, pairTrackRho0HistSpec, cprHistSpec);
    }
  };

  void processPhiSameEvent(FilteredCollision const& col, FemtoTracks const& tracks, FemtoPhis const& phis)
  {
    pairTrackPhiBuilder.processSameEvent<modes::Mode::kSe_Analysis>(col, tracks, trackPartition, phis, phiPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackTwoTrackResonance, processPhiSameEvent, "Enable processing same event processing for tracks and phis", true);

  void processPhiMixedEvent(FilteredCollisions const& cols, FemtoTracks const& tracks, FemtoPhis const& /*phis*/)
  {
    pairTrackPhiBuilder.processMixedEvent<modes::Mode::kMe_Analysis>(cols, tracks, trackPartition, phiPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackTwoTrackResonance, processPhiMixedEvent, "Enable processing mixed event processing for tracks and phis", true);

  void processKstar0SameEvent(FilteredCollision const& col, FemtoTracks const& tracks, FemtoKstar0s const& kstar0s)
  {
    pairTrackKstar0Builder.processSameEvent<modes::Mode::kSe_Analysis>(col, tracks, trackPartition, kstar0s, kstar0Partition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackTwoTrackResonance, processKstar0SameEvent, "Enable processing same event processing for tracks and kstar0s", false);

  void processKstar0MixedEvent(FilteredCollisions const& cols, FemtoTracks const& tracks, FemtoKstar0s const& /*kstar0s*/)
  {
    pairTrackKstar0Builder.processMixedEvent<modes::Mode::kMe_Analysis>(cols, tracks, trackPartition, kstar0Partition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackTwoTrackResonance, processKstar0MixedEvent, "Enable processing mixed event processing for tracks and kstar0s", false);

  void processRho0SameEvent(FilteredCollision const& col, FemtoTracks const& tracks, FemtoRho0s const& rho0s)
  {
    pairTrackRho0Builder.processSameEvent<modes::Mode::kSe_Analysis>(col, tracks, trackPartition, rho0s, rho0Partition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackTwoTrackResonance, processRho0SameEvent, "Enable processing same event processing for tracks and rho0s", false);

  void processRho0MixedEvent(FilteredCollisions const& cols, FemtoTracks const& tracks, FemtoRho0s const& /*rho0s*/)
  {
    pairTrackRho0Builder.processMixedEvent<modes::Mode::kMe_Analysis>(cols, tracks, trackPartition, rho0Partition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackTwoTrackResonance, processRho0MixedEvent, "Enable processing mixed event processing for tracks and rho0s", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairTrackTwoTrackResonance>(cfgc),
  };
  return workflow;
}
