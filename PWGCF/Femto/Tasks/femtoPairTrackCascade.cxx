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

#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femto;

struct FemtoPairTrackCascade {

  // setup tables
  using FemtoCollisions = Join<FCols, FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoTracks = o2::soa::Join<FTracks, FTrackMasks>;
  using FemtoXis = o2::soa::Join<FXis, FXiMasks>;
  using FemtoOmegas = o2::soa::Join<FOmegas, FOmegaMasks>;

  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, FTrackLabels>;
  using FemtoXisWithLabel = o2::soa::Join<FemtoXis, FLambdaLabels>;
  using FemtoOmegasWithLabel = o2::soa::Join<FemtoOmegas, FK0shortLabels>;

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

  // setup for daughters/bachelor
  trackhistmanager::ConfCascadePosDauBinning confPosDauBinning;
  trackhistmanager::ConfCascadeNegDauBinning confNegDauBinning;
  trackhistmanager::ConfCascadeBachelorBinning confBachelorBinning;

  // setup xis
  cascadebuilder::ConfXiSelection confXiSelection;
  cascadehistmanager::ConfXiBinning confXiBinning;
  particlecleaner::ConfXiCleaner1 confXiCleaner;

  Partition<FemtoXis> xiPartition = MAKE_CASCADE_PARTITION(confXiSelection);
  Preslice<FemtoXis> perColXis = aod::femtobase::stored::fColId;

  Partition<FemtoXisWithLabel> xiWithLabelPartition = MAKE_CASCADE_PARTITION(confXiSelection);
  Preslice<FemtoXisWithLabel> perColXisWithLabel = aod::femtobase::stored::fColId;

  // setup omegas
  cascadebuilder::ConfOmegaSelection confOmegaSelection;
  cascadehistmanager::ConfOmegaBinning confOmegaBinning;
  particlecleaner::ConfOmegaCleaner1 confOmegaCleaner;

  Partition<FemtoOmegas> omegaPartition = MAKE_CASCADE_PARTITION(confOmegaSelection);
  Preslice<FemtoOmegas> perColOmegas = aod::femtobase::stored::fColId;

  Partition<FemtoOmegasWithLabel> omegaWithLabelPartition = MAKE_CASCADE_PARTITION(confOmegaSelection);
  Preslice<FemtoOmegasWithLabel> perColOmegasWithLabel = aod::femtobase::stored::fColId;

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
    bool processData = doprocessXiSameEvent || doprocessXiMixedEvent || doprocessOmegaSameEvent || doprocessOmegaMixedEvent;
    bool processMc = doprocessXiSameEventMc || doprocessXiMixedEventMc || doprocessOmegaSameEventMc || doprocessOmegaMixedEventMc;

    if (processData && processMc) {
      LOG(fatal) << "Both data and mc processing is enabled. Breaking...";
    }

    bool processXi = doprocessXiSameEvent || doprocessXiSameEventMc || doprocessXiMixedEvent || doprocessXiMixedEventMc;
    bool processOmega = doprocessOmegaSameEvent || doprocessOmegaSameEventMc || doprocessOmegaMixedEvent || doprocessOmegaMixedEventMc;

    if (processXi && processOmega) {
      LOG(fatal) << "Both xi-track and omega-track processing is enabled. Breaking...";
    }

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins, confMixing.multBins.value}, true};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    auto cprHistSpecBachelor = closepairrejection::makeCprHistSpecMap(confCprBachelor);
    auto cprHistSpecV0Daughter = closepairrejection::makeCprHistSpecMap(confCprV0Daughter);

    if (processData) {
      auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
      auto trackHistSpec = trackhistmanager::makeTrackHistSpecMap(confTrackBinning);
      auto bachelorHistSpec = trackhistmanager::makeTrackHistSpecMap(confBachelorBinning);
      auto posDauSpec = trackhistmanager::makeTrackHistSpecMap(confPosDauBinning);
      auto negDauSpec = trackhistmanager::makeTrackHistSpecMap(confNegDauBinning);
      auto pairHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning);
      if (processXi) {
        auto xiHistSpec = cascadehistmanager::makeCascadeHistSpecMap(confXiBinning);
        auto pairTrackXiHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning);
        pairTrackXiBuilder.init<modes::Mode::kAnalysis>(&hRegistry, confTrackSelection, confTrackCleaner, confXiSelection, confXiCleaner, confCprBachelor, confCprV0Daughter, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, xiHistSpec, bachelorHistSpec, posDauSpec, negDauSpec, pairTrackXiHistSpec, cprHistSpecBachelor, cprHistSpecV0Daughter);
      } else {
        auto omegaHistSpec = cascadehistmanager::makeCascadeHistSpecMap(confOmegaBinning);
        auto pairTrackOmegaHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning);
        pairTrackOmegaBuilder.init<modes::Mode::kAnalysis>(&hRegistry, confTrackSelection, confTrackCleaner, confOmegaSelection, confOmegaCleaner, confCprBachelor, confCprV0Daughter, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, omegaHistSpec, bachelorHistSpec, posDauSpec, negDauSpec, pairTrackOmegaHistSpec, cprHistSpecBachelor, cprHistSpecV0Daughter);
      }
    } else {
      auto colHistSpec = colhistmanager::makeColMcHistSpecMap(confCollisionBinning);
      auto trackHistSpec = trackhistmanager::makeTrackMcHistSpecMap(confTrackBinning);
      auto bachelorHistSpec = trackhistmanager::makeTrackMcHistSpecMap(confBachelorBinning);
      auto posDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confPosDauBinning);
      auto negDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confNegDauBinning);
      auto pairHistSpec = pairhistmanager::makePairMcHistSpecMap(confPairBinning);
      if (processXi) {
        auto xiHistSpec = cascadehistmanager::makeCascadeMcHistSpecMap(confXiBinning);
        auto pairTrackXiHistSpec = pairhistmanager::makePairMcHistSpecMap(confPairBinning);
        pairTrackXiBuilder.init<modes::Mode::kAnalysis_Mc>(&hRegistry, confTrackSelection, confTrackCleaner, confXiSelection, confXiCleaner, confCprBachelor, confCprV0Daughter, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, xiHistSpec, bachelorHistSpec, posDauSpec, negDauSpec, pairTrackXiHistSpec, cprHistSpecBachelor, cprHistSpecV0Daughter);
      } else {
        auto omegaHistSpec = cascadehistmanager::makeCascadeMcHistSpecMap(confOmegaBinning);
        auto pairTrackOmegaHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning);
        pairTrackOmegaBuilder.init<modes::Mode::kAnalysis_Mc>(&hRegistry, confTrackSelection, confTrackCleaner, confOmegaSelection, confOmegaCleaner, confCprBachelor, confCprV0Daughter, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, omegaHistSpec, bachelorHistSpec, posDauSpec, negDauSpec, pairTrackOmegaHistSpec, cprHistSpecBachelor, cprHistSpecV0Daughter);
      }
    }
  };

  void processXiSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoXis const& xis)
  {
    pairTrackXiBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, trackPartition, xis, xiPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processXiSameEvent, "Enable processing same event processing for tracks and xis", true);

  void processXiSameEventMc(FilteredFemtoCollisionWithLabel const& col, FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoXisWithLabel const& xis, FMcParticles const& mcParticles, FMcMothers const& mcMothers, FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackXiBuilder.processSameEvent<modes::Mode::kAnalysis_Mc>(col, mcCols, tracks, trackWithLabelPartition, xis, xiWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processXiSameEventMc, "Enable processing same event processing for tracks and xis with mc information", false);

  void processXiMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoXis const& /*xis*/)
  {
    pairTrackXiBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, trackPartition, xiPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processXiMixedEvent, "Enable processing mixed event processing for tracks and xis", true);

  void processXiMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoXisWithLabel const& /*xis*/, FMcParticles const& mcParticles, FMcMothers const& mcMothers, FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackXiBuilder.processMixedEvent<modes::Mode::kAnalysis_Mc>(cols, mcCols, tracks, trackWithLabelPartition, xiWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processXiMixedEventMc, "Enable processing mixed event processing for tracks and xis with mc information", false);

  void processOmegaSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoOmegas const& omegas)
  {
    pairTrackOmegaBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, trackPartition, omegas, omegaPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processOmegaSameEvent, "Enable processing same event processing for tracks and omegas", false);

  void processOmegaSameEventMc(FilteredFemtoCollisionWithLabel const& col, FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoOmegasWithLabel const& omegas, FMcParticles const& mcParticles, FMcMothers const& mcMothers, FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackOmegaBuilder.processSameEvent<modes::Mode::kAnalysis_Mc>(col, mcCols, tracks, trackWithLabelPartition, omegas, omegaWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processOmegaSameEventMc, "Enable processing same event processing for tracks and omegas with mc information", false);

  void processOmegaMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoOmegas const& /*omegas*/)
  {
    pairTrackOmegaBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, trackPartition, omegaPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processOmegaMixedEvent, "Enable processing mixed event processing for tracks and omegas", false);

  void processOmegaMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoXisWithLabel const& /*xis*/, FMcParticles const& mcParticles, FMcMothers const& mcMothers, FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackOmegaBuilder.processMixedEvent<modes::Mode::kAnalysis_Mc>(cols, mcCols, tracks, trackWithLabelPartition, omegaWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processOmegaMixedEventMc, "Enable processing mixed event processing for tracks and omegas with mc information", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairTrackCascade>(cfgc),
  };
  return workflow;
}
