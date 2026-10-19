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

/// \file femtoPairTrackTwoPhotonResonance.cxx
/// \brief Task that computes correlation between tracks and resonances decaying into two photons (pi0, eta)
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/closePairRejection.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/pairBuilder.h"
#include "PWGCF/Femto/Core/pairCleaner.h"
#include "PWGCF/Femto/Core/pairHistManager.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/photonHistManager.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/Core/twoPhotonResonanceBuilder.h"
#include "PWGCF/Femto/Core/twoPhotonResonanceHistManager.h"
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

struct FemtoPairTrackTwoPhotonResonance {

  // setup tables
  using Collisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks>;
  using FemtoPhotons = o2::soa::Join<o2::aod::FPhotons, o2::aod::FPhotonMasks>;
  using FemtoPi0s = o2::soa::Join<o2::aod::FPi0s, o2::aod::FPi0Masks>;
  using FemtoEtas = o2::soa::Join<o2::aod::FEtas, o2::aod::FEtaMasks>;

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

  // setup for photon daughters (shared between pi0 and eta -- unordered, no pos/neg or species asymmetry)
  photonhistmanager::ConfPhotonBinning confPhotonDauBinning;
  photonhistmanager::ConfPhotonQaBinning confPhotonDauQaBinning;

  // setup pi0s
  twophotonresonancebuilder::ConfPi0Selection pi0Selection;
  twophotonresonancehistmanager::ConfPi0Binning confPi0Binning;
  o2::framework::Partition<FemtoPi0s> pi0Partition = MAKE_TWOPHOTONRESONANCE_PARTITION(pi0Selection);
  o2::framework::Preslice<FemtoPi0s> perColPi0s = o2::aod::femtobase::stored::fColId;

  // setup etas
  twophotonresonancebuilder::ConfEtaSelection etaSelection;
  twophotonresonancehistmanager::ConfEtaBinning confEtaBinning;
  o2::framework::Partition<FemtoEtas> etaPartition = MAKE_TWOPHOTONRESONANCE_PARTITION(etaSelection);
  o2::framework::Preslice<FemtoEtas> perColEtas = o2::aod::femtobase::stored::fColId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::ConfPairCuts confPairCuts;
  paircleaner::ConfPairCleanerBinning confPairCleaner;

  // setup for track-pi0 pairs
  pairbuilder::PairTrackTwoPhotonResonanceBuilder<
    trackhistmanager::PrefixTrack1,
    twophotonresonancehistmanager::PrefixPi0,
    photonhistmanager::PrefixTwoPhotonResonanceDau1Qa,
    photonhistmanager::PrefixTwoPhotonResonanceDau2Qa,
    pairhistmanager::PrefixTrackPi0Se,
    pairhistmanager::PrefixTrackPi0Me,
    closepairrejection::PrefixTrackTwoPhotonResonanceSe,
    closepairrejection::PrefixTrackTwoPhotonResonanceMe,
    modes::TwoPhotonResonance::kPi0>
    pairTrackPi0Builder;

  // setup for track-eta pairs
  pairbuilder::PairTrackTwoPhotonResonanceBuilder<
    trackhistmanager::PrefixTrack1,
    twophotonresonancehistmanager::PrefixEta,
    photonhistmanager::PrefixTwoPhotonResonanceDau1Qa,
    photonhistmanager::PrefixTwoPhotonResonanceDau2Qa,
    pairhistmanager::PrefixTrackEtaSe,
    pairhistmanager::PrefixTrackEtaMe,
    closepairrejection::PrefixTrackTwoPhotonResonanceSe,
    closepairrejection::PrefixTrackTwoPhotonResonanceMe,
    modes::TwoPhotonResonance::kEta>
    pairTrackEtaBuilder;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult, o2::aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}};
  pairhistmanager::ConfMixing confMixing;

  o2::framework::HistogramRegistry hRegistry{"FemtoTrackTwoPhotonResonance", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  // setup cpr (no-op for two-photon resonances: PCM photon daughters are not indexed into the track
  // table, so there is no trajectory to compare against a track's -- see closePairRejection.h)
  closepairrejection::ConfCprTrackTwoPhotonResonanceDaughter confCpr;

  void init(o2::framework::InitContext&)
  {

    if (static_cast<int>(doprocessPi0SameEvent || doprocessPi0MixedEvent) + static_cast<int>(doprocessEtaSameEvent || doprocessEtaMixedEvent) > 1) {
      LOG(fatal) << "Can only process SE/ME for pi0-tracks or eta-tracks";
    }

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins.value, confMixing.multBins.value}};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}};

    // setup histograms
    auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
    auto trackHistSpec = trackhistmanager::makeTrackHistSpecMap(confTrackBinning);
    auto dauHistSpec = photonhistmanager::makePhotonHistSpecMap(confPhotonDauBinning);
    auto cprHistSpec = closepairrejection::makeCprHistSpecMap(confCpr);
    auto pairCleanerHistSpec = paircleaner::makePairCleanerHistSpecMap(confPairCleaner);

    // setup for pi0
    if (doprocessPi0SameEvent || doprocessPi0MixedEvent) {
      auto pi0HistSpec = twophotonresonancehistmanager::makeTwoPhotonResonanceHistSpecMap(confPi0Binning);
      auto pairTrackPi0HistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);
      pairTrackPi0Builder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, trackSelection, pi0Selection, confCpr, confMixing, confPairBinning, confPairCuts, confPhotonDauQaBinning, colHistSpec, trackHistSpec, pi0HistSpec, dauHistSpec, pairTrackPi0HistSpec, cprHistSpec, pairCleanerHistSpec);
    }

    // setup for eta
    if (doprocessEtaSameEvent || doprocessEtaMixedEvent) {
      auto etaHistSpec = twophotonresonancehistmanager::makeTwoPhotonResonanceHistSpecMap(confEtaBinning);
      auto pairTrackEtaHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);
      pairTrackEtaBuilder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, trackSelection, etaSelection, confCpr, confMixing, confPairBinning, confPairCuts, confPhotonDauQaBinning, colHistSpec, trackHistSpec, etaHistSpec, dauHistSpec, pairTrackEtaHistSpec, cprHistSpec, pairCleanerHistSpec);
    }

    hRegistry.print();
  };

  void processPi0SameEvent(FilteredCollision const& col, FemtoTracks const& /*tracks*/, FemtoPi0s const& /*pi0s*/, FemtoPhotons const& photons)
  {
    pairTrackPi0Builder.processSameEvent<modes::Mode::kSe_Reco>(col, trackPartition, pi0Partition, photons, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackTwoPhotonResonance, processPi0SameEvent, "Enable same event processing for tracks and pi0s", true);

  void processPi0MixedEvent(FilteredCollisions const& cols, FemtoTracks const& /*tracks*/, FemtoPi0s const& /*pi0s*/, FemtoPhotons const& photons)
  {
    pairTrackPi0Builder.processMixedEvent<modes::Mode::kMe_Reco>(cols, trackPartition, pi0Partition, photons, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackTwoPhotonResonance, processPi0MixedEvent, "Enable mixed event processing for tracks and pi0s", true);

  void processEtaSameEvent(FilteredCollision const& col, FemtoTracks const& /*tracks*/, FemtoEtas const& /*etas*/, FemtoPhotons const& photons)
  {
    pairTrackEtaBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, trackPartition, etaPartition, photons, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackTwoPhotonResonance, processEtaSameEvent, "Enable same event processing for tracks and etas", false);

  void processEtaMixedEvent(FilteredCollisions const& cols, FemtoTracks const& /*tracks*/, FemtoEtas const& /*etas*/, FemtoPhotons const& photons)
  {
    pairTrackEtaBuilder.processMixedEvent<modes::Mode::kMe_Reco>(cols, trackPartition, etaPartition, photons, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackTwoPhotonResonance, processEtaMixedEvent, "Enable mixed event processing for tracks and etas", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairTrackTwoPhotonResonance>(context),
  };
  return workflow;
}
