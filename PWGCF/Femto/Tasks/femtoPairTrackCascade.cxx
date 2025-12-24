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

/// \file femtoPairTrackCascade.cxx
/// \brief Tasks that computes correlation between tracks and cascades
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/cascadeBuilder.h"
#include "PWGCF/Femto/Core/cascadeHistManager.h"
#include "PWGCF/Femto/Core/closePairRejection.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
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

struct FemtoPairTrackCascade {

  // setup tables
  using Collisions = Join<FCols, FColMasks>;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  using Tracks = o2::soa::Join<FTracks, FTrackMasks>;
  using Xis = o2::soa::Join<FXis, FXiMasks>;
  using Omegas = o2::soa::Join<FOmegas, FOmegaMasks>;

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

  // setup for daughters/bachelor
  trackhistmanager::ConfCascadePosDauBinning confPosDauBinning;
  trackhistmanager::ConfCascadeNegDauBinning confNegDauBinning;
  trackhistmanager::ConfCascadeBachelorBinning confBachelorBinning;

  // setup xis
  cascadebuilder::ConfXiSelection xiSelection;
  cascadehistmanager::ConfXiBinning confXiBinning;
  Partition<Xis> xiPartition = MAKE_CASCADE_PARTITION(xiSelection);
  Preslice<Xis> perColXis = aod::femtobase::stored::fColId;

  // setup omegas
  cascadebuilder::ConfOmegaSelection omegaSelection;
  cascadehistmanager::ConfOmegaBinning confOmegaBinning;
  Partition<Omegas> omegaPartition = MAKE_CASCADE_PARTITION(omegaSelection);
  Preslice<Omegas> perColOmegas = aod::femtobase::stored::fColId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::ConfPairCuts confPairCuts;

  pairbuilder::PairTrackCascadeBuilder<
    trackhistmanager::PrefixTrack1,
    cascadehistmanager::PrefixXi,
    trackhistmanager::PrefixCascadeBachelor,
    trackhistmanager::PrefixCascadePosDaughter,
    trackhistmanager::PrefixCascadeNegDaughter,
    pairhistmanager::PrefixTrackCascadeSe,
    pairhistmanager::PrefixTrackCascadeMe,
    closepairrejection::PrefixTrackCascadeBachelorSe,
    closepairrejection::PrefixTrackV0DaughterSe,
    closepairrejection::PrefixTrackCascadeBachelorMe,
    closepairrejection::PrefixTrackV0DaughterMe,
    modes::Cascade::kXi>
    pairTrackXiBuilder;

  pairbuilder::PairTrackCascadeBuilder<
    trackhistmanager::PrefixTrack1,
    cascadehistmanager::PrefixOmega,
    trackhistmanager::PrefixCascadeBachelor,
    trackhistmanager::PrefixCascadePosDaughter,
    trackhistmanager::PrefixCascadeNegDaughter,
    pairhistmanager::PrefixTrackCascadeSe,
    pairhistmanager::PrefixTrackCascadeMe,
    closepairrejection::PrefixTrackCascadeBachelorSe,
    closepairrejection::PrefixTrackV0DaughterSe,
    closepairrejection::PrefixTrackCascadeBachelorMe,
    closepairrejection::PrefixTrackV0DaughterMe,
    modes::Cascade::kOmega>
    pairTrackOmegaBuilder;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  ColumnBinningPolicy<femtocollisions::PosZ, femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}, true};
  ColumnBinningPolicy<aod::femtocollisions::PosZ, aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}, true};
  ColumnBinningPolicy<aod::femtocollisions::PosZ, aod::femtocollisions::Mult, aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}, true};
  pairhistmanager::ConfMixing confMixing;

  HistogramRegistry hRegistry{"FemtoTrackCascade", {}, OutputObjHandlingPolicy::AnalysisObject};

  // setup cpr
  closepairrejection::ConfCprTrackCascadeBachelor confCprBachelor;
  closepairrejection::ConfCprTrackV0Daughter confCprV0Daughter;

  void init(InitContext&)
  {

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins, confMixing.multBins.value}, true};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    // setup histograms
    auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
    auto trackHistSpec = trackhistmanager::makeTrackHistSpecMap(confTrackBinning);
    auto bachelorHistSpec = trackhistmanager::makeTrackHistSpecMap(confBachelorBinning);
    auto posDauSpec = trackhistmanager::makeTrackHistSpecMap(confPosDauBinning);
    auto negDauSpec = trackhistmanager::makeTrackHistSpecMap(confNegDauBinning);
    auto pairHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning);
    auto cprHistSpecBachelor = closepairrejection::makeCprHistSpecMap(confCprBachelor);
    auto cprHistSpecV0Daughter = closepairrejection::makeCprHistSpecMap(confCprV0Daughter);

    // setup for xis
    if (doprocessXiSameEvent || doprocessXiMixedEvent) {
      auto xiHistSpec = cascadehistmanager::makeCascadeHistSpecMap(confXiBinning);
      auto pairTrackXiHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning);
      pairTrackXiBuilder.init<modes::Mode::kAnalysis>(&hRegistry, trackSelection, xiSelection, confCprBachelor, confCprV0Daughter, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, xiHistSpec, bachelorHistSpec, posDauSpec, negDauSpec, pairTrackXiHistSpec, cprHistSpecBachelor, cprHistSpecV0Daughter);
    }

    // setup for omegas
    if (doprocessOmegaSameEvent || doprocessOmegaMixedEvent) {
      auto omegaHistSpec = cascadehistmanager::makeCascadeHistSpecMap(confOmegaBinning);
      auto pairTrackOmegaHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning);
      pairTrackOmegaBuilder.init<modes::Mode::kAnalysis>(&hRegistry, trackSelection, xiSelection, confCprBachelor, confCprV0Daughter, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, omegaHistSpec, bachelorHistSpec, posDauSpec, negDauSpec, pairTrackOmegaHistSpec, cprHistSpecBachelor, cprHistSpecV0Daughter);
    }

    if (((doprocessXiSameEvent || doprocessXiMixedEvent) + (doprocessOmegaSameEvent || doprocessOmegaMixedEvent)) > 1) {
      LOG(fatal) << "Can only process xi-tracks Or omega-tracks";
    }
  };

  void processXiSameEvent(FilteredCollision const& col, Tracks const& tracks, Xis const& xis)
  {
    pairTrackXiBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, trackPartition, xis, xiPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processXiSameEvent, "Enable processing same event processing for tracks and xis", true);

  void processXiMixedEvent(FilteredCollisions const& cols, Tracks const& tracks, Xis const& /*xis*/)
  {
    pairTrackXiBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, trackPartition, xiPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processXiMixedEvent, "Enable processing mixed event processing for tracks and xis", true);

  void processOmegaSameEvent(FilteredCollision const& col, Tracks const& tracks, Omegas const& omegas)
  {
    pairTrackOmegaBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, trackPartition, omegas, omegaPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processOmegaSameEvent, "Enable processing same event processing for tracks and omegas", false);

  void processOmegaMixedEvent(FilteredCollisions const& cols, Tracks const& tracks, Omegas const& /*omegas*/)
  {
    pairTrackOmegaBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, trackPartition, omegaPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processOmegaMixedEvent, "Enable processing mixed event processing for tracks and omegas", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairTrackCascade>(cfgc),
  };
  return workflow;
}
