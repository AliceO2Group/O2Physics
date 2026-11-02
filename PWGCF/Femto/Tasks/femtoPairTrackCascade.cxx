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
#include "PWGCF/Femto/Core/pairCleaner.h"
#include "PWGCF/Femto/Core/pairHistManager.h"
#include "PWGCF/Femto/Core/particleCleaner.h"
#include "PWGCF/Femto/Core/partitions.h"
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

struct FemtoPairTrackCascade {

  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, o2::aod::FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks>;
  using FemtoXis = o2::soa::Join<o2::aod::FXis, o2::aod::FXiMasks>;
  using FemtoOmegas = o2::soa::Join<o2::aod::FOmegas, o2::aod::FOmegaMasks>;

  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, o2::aod::FTrackLabels, o2::aod::FItsPurity, o2::aod::FTpcPurity, o2::aod::FTofPurity, o2::aod::FTrackDcas>;
  using FemtoXisWithLabel = o2::soa::Join<FemtoXis, o2::aod::FXiLabels, o2::aod::FXiFrac>;
  using FemtoOmegasWithLabel = o2::soa::Join<FemtoOmegas, o2::aod::FOmegaLabels, o2::aod::FOmegaFrac>;

  // reco-only + FitVars (no MC labels, but carries the opt-in purity/fraction tables)
  using FemtoTracksWithFitVars = o2::soa::Join<FemtoTracks, o2::aod::FItsPurity, o2::aod::FTpcPurity, o2::aod::FTofPurity, o2::aod::FTrackDcas>;
  using FemtoXisWithFitVars = o2::soa::Join<FemtoXis, o2::aod::FXiFrac>;
  using FemtoOmegasWithFitVars = o2::soa::Join<FemtoOmegas, o2::aod::FOmegaFrac>;

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
  o2::framework::Preslice<FemtoTracksWithLabel> perColtracksWithLabel = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoTracksWithFitVars> trackWithFitVarsPartition = MAKE_TRACK_PARTITION(confTrackSelection);
  o2::framework::Preslice<FemtoTracksWithFitVars> perColTracksWithFitVars = o2::aod::femtobase::stored::fColId;

  // setup for daughters/bachelor
  trackhistmanager::ConfCascadePosDauBinning confPosDauBinning;
  trackhistmanager::ConfCascadeNegDauBinning confNegDauBinning;
  trackhistmanager::ConfCascadeBachelorBinning confBachelorBinning;

  // setup xis
  cascadebuilder::ConfXiSelection confXiSelection;
  cascadehistmanager::ConfXiBinning confXiBinning;
  particlecleaner::ConfXiCleaner1 confXiCleaner;

  o2::framework::Partition<FemtoXis> xiPartition = MAKE_CASCADE_PARTITION(confXiSelection);
  o2::framework::Preslice<FemtoXis> perColXis = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoXisWithLabel> xiWithLabelPartition = MAKE_CASCADE_PARTITION(confXiSelection);
  o2::framework::Preslice<FemtoXisWithLabel> perColXisWithLabel = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoXisWithFitVars> xiWithFitVarsPartition = MAKE_CASCADE_PARTITION(confXiSelection);
  o2::framework::Preslice<FemtoXisWithFitVars> perColXisWithFitVars = o2::aod::femtobase::stored::fColId;

  // setup omegas
  cascadebuilder::ConfOmegaSelection confOmegaSelection;
  cascadehistmanager::ConfOmegaBinning confOmegaBinning;
  particlecleaner::ConfOmegaCleaner1 confOmegaCleaner;

  o2::framework::Partition<FemtoOmegas> omegaPartition = MAKE_CASCADE_PARTITION(confOmegaSelection);
  o2::framework::Preslice<FemtoOmegas> perColOmegas = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoOmegasWithLabel> omegaWithLabelPartition = MAKE_CASCADE_PARTITION(confOmegaSelection);
  o2::framework::Preslice<FemtoOmegasWithLabel> perColOmegasWithLabel = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoOmegasWithFitVars> omegaWithFitVarsPartition = MAKE_CASCADE_PARTITION(confOmegaSelection);
  o2::framework::Preslice<FemtoOmegasWithFitVars> perColOmegasWithFitVars = o2::aod::femtobase::stored::fColId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::ConfPairCuts confPairCuts;
  paircleaner::ConfPairCleanerBinning confPairCleaner;

  pairbuilder::PairTrackCascadeBuilder<
    trackhistmanager::PrefixTrack1Se,
    cascadehistmanager::PrefixXiSe,
    trackhistmanager::PrefixCascadeBachelorSe,
    trackhistmanager::PrefixCascadePosDaughterSe,
    trackhistmanager::PrefixCascadeNegDaughterSe,
    trackhistmanager::PrefixTrack1Me,
    cascadehistmanager::PrefixXiMe,
    trackhistmanager::PrefixCascadeBachelorMe,
    trackhistmanager::PrefixCascadePosDaughterMe,
    trackhistmanager::PrefixCascadeNegDaughterMe,
    pairhistmanager::PrefixTrackCascadeSe,
    pairhistmanager::PrefixTrackCascadeMe,
    closepairrejection::PrefixTrackCascadeBachelorSe,
    closepairrejection::PrefixTrackV0DaughterSe,
    closepairrejection::PrefixTrackCascadeBachelorMe,
    closepairrejection::PrefixTrackV0DaughterMe,
    modes::Cascade::kXi>
    pairTrackXiBuilder;

  pairbuilder::PairTrackCascadeBuilder<
    trackhistmanager::PrefixTrack1Se,
    cascadehistmanager::PrefixOmegaSe,
    trackhistmanager::PrefixCascadeBachelorSe,
    trackhistmanager::PrefixCascadePosDaughterSe,
    trackhistmanager::PrefixCascadeNegDaughterSe,
    trackhistmanager::PrefixTrack1Me,
    cascadehistmanager::PrefixOmegaMe,
    trackhistmanager::PrefixCascadeBachelorMe,
    trackhistmanager::PrefixCascadePosDaughterMe,
    trackhistmanager::PrefixCascadeNegDaughterMe,
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
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult, o2::aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}};
  pairhistmanager::ConfMixing confMixing;

  o2::framework::HistogramRegistry hRegistry{"FemtoTrackCascade", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  // setup cpr
  closepairrejection::ConfCprTrackCascadeBachelor confCprBachelor;
  closepairrejection::ConfCprTrackV0Daughter confCprV0Daughter;

  void init(o2::framework::InitContext&)
  {
    bool processData = doprocessXiSameEvent || doprocessXiMixedEvent || doprocessOmegaSameEvent || doprocessOmegaMixedEvent;
    bool processFitVars = doprocessXiSameEventWithFitVars || doprocessXiMixedEventWithFitVars || doprocessOmegaSameEventWithFitVars || doprocessOmegaMixedEventWithFitVars;
    bool processMc = doprocessXiSameEventMc || doprocessXiMixedEventMc || doprocessOmegaSameEventMc || doprocessOmegaMixedEventMc;

    if (static_cast<int>(processData) + static_cast<int>(processFitVars) + static_cast<int>(processMc) > 1) {
      LOG(fatal) << "More than one of data, FitVars or mc processing is activated. Breaking...";
    }

    bool processXi = doprocessXiSameEvent || doprocessXiSameEventWithFitVars || doprocessXiSameEventMc || doprocessXiMixedEvent || doprocessXiMixedEventWithFitVars || doprocessXiMixedEventMc;
    bool processOmega = doprocessOmegaSameEvent || doprocessOmegaSameEventWithFitVars || doprocessOmegaSameEventMc || doprocessOmegaMixedEvent || doprocessOmegaMixedEventWithFitVars || doprocessOmegaMixedEventMc;

    if (processXi && processOmega) {
      LOG(fatal) << "Both xi-track and omega-track processing is enabled. Breaking...";
    }

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins.value, confMixing.multBins.value}};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}};

    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> trackHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> bachelorHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> posDauSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> negDauSpec;
    std::map<cascadehistmanager::CascadeHist, std::vector<o2::framework::AxisSpec>> xiHistSpec;
    std::map<cascadehistmanager::CascadeHist, std::vector<o2::framework::AxisSpec>> omegaHistSpec;
    std::map<pairhistmanager::PairHist, std::vector<o2::framework::AxisSpec>> pairTrackCascadeHistSpec;

    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpecBachelor = closepairrejection::makeCprHistSpecMap(confCprBachelor);
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpecV0Daughter = closepairrejection::makeCprHistSpecMap(confCprV0Daughter);
    std::map<paircleaner::PairCleanerHist, std::vector<o2::framework::AxisSpec>> pairCleanerHistSpec = paircleaner::makePairCleanerHistSpecMap(confPairCleaner);

    if (processFitVars) {
      colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
      trackHistSpec = trackhistmanager::makeTrackHistSpecMap(confTrackBinning);
      bachelorHistSpec = trackhistmanager::makeTrackHistSpecMap(confBachelorBinning);
      posDauSpec = trackhistmanager::makeTrackHistSpecMap(confPosDauBinning);
      negDauSpec = trackhistmanager::makeTrackHistSpecMap(confNegDauBinning);
      pairTrackCascadeHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);
      if (processXi) {
        xiHistSpec = cascadehistmanager::makeCascadeHistSpecMap(confXiBinning);
        pairTrackXiBuilder.init<modes::Mode::kSe_Reco_FitVars, modes::Mode::kMe_Reco_FitVars>(&hRegistry, confCollisionBinning, confTrackSelection, confTrackCleaner, confXiSelection, confXiCleaner, confCprBachelor, confCprV0Daughter, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, xiHistSpec, bachelorHistSpec, posDauSpec, negDauSpec, pairTrackCascadeHistSpec, cprHistSpecBachelor, cprHistSpecV0Daughter, pairCleanerHistSpec);
      } else {
        omegaHistSpec = cascadehistmanager::makeCascadeHistSpecMap(confOmegaBinning);
        pairTrackOmegaBuilder.init<modes::Mode::kSe_Reco_FitVars, modes::Mode::kMe_Reco_FitVars>(&hRegistry, confCollisionBinning, confTrackSelection, confTrackCleaner, confOmegaSelection, confOmegaCleaner, confCprBachelor, confCprV0Daughter, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, omegaHistSpec, bachelorHistSpec, posDauSpec, negDauSpec, pairTrackCascadeHistSpec, cprHistSpecBachelor, cprHistSpecV0Daughter, pairCleanerHistSpec);
      }
    } else if (processData) {
      colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
      trackHistSpec = trackhistmanager::makeTrackHistSpecMap(confTrackBinning);
      bachelorHistSpec = trackhistmanager::makeTrackHistSpecMap(confBachelorBinning);
      posDauSpec = trackhistmanager::makeTrackHistSpecMap(confPosDauBinning);
      negDauSpec = trackhistmanager::makeTrackHistSpecMap(confNegDauBinning);
      pairTrackCascadeHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);
      if (processXi) {
        xiHistSpec = cascadehistmanager::makeCascadeHistSpecMap(confXiBinning);
        pairTrackXiBuilder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, confTrackSelection, confTrackCleaner, confXiSelection, confXiCleaner, confCprBachelor, confCprV0Daughter, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, xiHistSpec, bachelorHistSpec, posDauSpec, negDauSpec, pairTrackCascadeHistSpec, cprHistSpecBachelor, cprHistSpecV0Daughter, pairCleanerHistSpec);
      } else {
        omegaHistSpec = cascadehistmanager::makeCascadeHistSpecMap(confOmegaBinning);
        pairTrackCascadeHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);
        pairTrackOmegaBuilder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, confTrackSelection, confTrackCleaner, confOmegaSelection, confOmegaCleaner, confCprBachelor, confCprV0Daughter, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, omegaHistSpec, bachelorHistSpec, posDauSpec, negDauSpec, pairTrackCascadeHistSpec, cprHistSpecBachelor, cprHistSpecV0Daughter, pairCleanerHistSpec);
      }
    } else {
      colHistSpec = colhistmanager::makeColMcHistSpecMap(confCollisionBinning);
      trackHistSpec = trackhistmanager::makeTrackMcHistSpecMap(confTrackBinning);
      bachelorHistSpec = trackhistmanager::makeTrackMcHistSpecMap(confBachelorBinning);
      posDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confPosDauBinning);
      negDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confNegDauBinning);
      pairTrackCascadeHistSpec = pairhistmanager::makePairMcHistSpecMap(confPairBinning, confMixing);
      if (processXi) {
        xiHistSpec = cascadehistmanager::makeCascadeMcHistSpecMap(confXiBinning);
        pairTrackXiBuilder.init<modes::Mode::kSe_Reco_Mc_FitVars, modes::Mode::kMe_Reco_Mc_FitVars>(&hRegistry, confCollisionBinning, confTrackSelection, confTrackCleaner, confXiSelection, confXiCleaner, confCprBachelor, confCprV0Daughter, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, xiHistSpec, bachelorHistSpec, posDauSpec, negDauSpec, pairTrackCascadeHistSpec, cprHistSpecBachelor, cprHistSpecV0Daughter, pairCleanerHistSpec);
      } else {
        omegaHistSpec = cascadehistmanager::makeCascadeMcHistSpecMap(confOmegaBinning);
        pairTrackOmegaBuilder.init<modes::Mode::kSe_Reco_Mc_FitVars, modes::Mode::kMe_Reco_Mc_FitVars>(&hRegistry, confCollisionBinning, confTrackSelection, confTrackCleaner, confOmegaSelection, confOmegaCleaner, confCprBachelor, confCprV0Daughter, confMixing, confPairBinning, confPairCuts, colHistSpec, trackHistSpec, omegaHistSpec, bachelorHistSpec, posDauSpec, negDauSpec, pairTrackCascadeHistSpec, cprHistSpecBachelor, cprHistSpecV0Daughter, pairCleanerHistSpec);
      }
    }
    hRegistry.print();
  };

  void processXiSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoXis const& xis)
  {
    pairTrackXiBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, trackPartition, xis, xiPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processXiSameEvent, "Enable processing same event processing for tracks and xis", true);

  void processXiSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoXisWithLabel const& xis, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackXiBuilder.processSameEvent<modes::Mode::kSe_Reco_Mc_FitVars>(col, mcCols, tracks, trackWithLabelPartition, xis, xiWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processXiSameEventMc, "Enable processing same event processing for tracks and xis with mc information", false);

  void processXiMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoXis const& /*xis*/)
  {
    pairTrackXiBuilder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, trackPartition, xiPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processXiMixedEvent, "Enable processing mixed event processing for tracks and xis", true);

  void processXiMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoXisWithLabel const& /*xis*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackXiBuilder.processMixedEvent<modes::Mode::kMe_Reco_Mc_FitVars>(cols, mcCols, tracks, trackWithLabelPartition, xiWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processXiMixedEventMc, "Enable processing mixed event processing for tracks and xis with mc information", false);

  void processXiSameEventWithFitVars(FilteredFemtoCollision const& col, FemtoTracksWithFitVars const& tracks, FemtoXisWithFitVars const& xis)
  {
    pairTrackXiBuilder.processSameEvent<modes::Mode::kSe_Reco_FitVars>(col, tracks, trackWithFitVarsPartition, xis, xiWithFitVarsPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processXiSameEventWithFitVars, "Enable processing same event processing for tracks and xis, carrying through FitVars tables", false);

  void processXiMixedEventWithFitVars(FilteredFemtoCollisions const& cols, FemtoTracksWithFitVars const& tracks, FemtoXisWithFitVars const& /*xis*/)
  {
    pairTrackXiBuilder.processMixedEvent<modes::Mode::kMe_Reco_FitVars>(cols, tracks, trackWithFitVarsPartition, xiWithFitVarsPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processXiMixedEventWithFitVars, "Enable processing mixed event processing for tracks and xis, carrying through FitVars tables", false);

  void processOmegaSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoOmegas const& omegas)
  {
    pairTrackOmegaBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, trackPartition, omegas, omegaPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processOmegaSameEvent, "Enable processing same event processing for tracks and omegas", false);

  void processOmegaSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoOmegasWithLabel const& omegas, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackOmegaBuilder.processSameEvent<modes::Mode::kSe_Reco_Mc_FitVars>(col, mcCols, tracks, trackWithLabelPartition, omegas, omegaWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processOmegaSameEventMc, "Enable processing same event processing for tracks and omegas with mc information", false);

  void processOmegaMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoOmegas const& /*omegas*/)
  {
    pairTrackOmegaBuilder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, trackPartition, omegaPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processOmegaMixedEvent, "Enable processing mixed event processing for tracks and omegas", false);

  void processOmegaMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoOmegasWithLabel const& /*omegas*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairTrackOmegaBuilder.processMixedEvent<modes::Mode::kMe_Reco_Mc_FitVars>(cols, mcCols, tracks, trackWithLabelPartition, omegaWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processOmegaMixedEventMc, "Enable processing mixed event processing for tracks and omegas with mc information", false);

  void processOmegaSameEventWithFitVars(FilteredFemtoCollision const& col, FemtoTracksWithFitVars const& tracks, FemtoOmegasWithFitVars const& omegas)
  {
    pairTrackOmegaBuilder.processSameEvent<modes::Mode::kSe_Reco_FitVars>(col, tracks, trackWithFitVarsPartition, omegas, omegaWithFitVarsPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processOmegaSameEventWithFitVars, "Enable processing same event processing for tracks and omegas, carrying through FitVars tables", false);

  void processOmegaMixedEventWithFitVars(FilteredFemtoCollisions const& cols, FemtoTracksWithFitVars const& tracks, FemtoOmegasWithFitVars const& /*omegas*/)
  {
    pairTrackOmegaBuilder.processMixedEvent<modes::Mode::kMe_Reco_FitVars>(cols, tracks, trackWithFitVarsPartition, omegaWithFitVarsPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairTrackCascade, processOmegaMixedEventWithFitVars, "Enable processing mixed event processing for tracks and omegas, carrying through FitVars tables", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairTrackCascade>(context),
  };
  return workflow;
}
