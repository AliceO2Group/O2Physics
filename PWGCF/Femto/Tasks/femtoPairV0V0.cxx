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

/// \file femtoPairV0V0.cxx
/// \brief Tasks that computes correlation between two v0s
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/closePairRejection.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/pairBuilder.h"
#include "PWGCF/Femto/Core/pairHistManager.h"
#include "PWGCF/Femto/Core/particleCleaner.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/Core/v0Builder.h"
#include "PWGCF/Femto/Core/v0HistManager.h"
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

struct FemtoPairV0V0 {

  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, o2::aod::FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks>;
  using FemtoLambdas = o2::soa::Join<o2::aod::FLambdas, o2::aod::FLambdaMasks>;
  using FemtoK0shorts = o2::soa::Join<o2::aod::FK0shorts, o2::aod::FK0shortMasks>;

  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, o2::aod::FTrackLabels>;
  using FemtoLambdasWithLabel = o2::soa::Join<FemtoLambdas, o2::aod::FLambdaLabels>;
  using FemtoK0shortsWithLabel = o2::soa::Join<FemtoK0shorts, o2::aod::FK0shortLabels>;
  //
  o2::framework::SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // setup for daughters
  trackhistmanager::ConfV0PosDauBinning confPosDauBinning;
  trackhistmanager::ConfV0NegDauBinning confNegDauBinning;

  // setup lambdas
  v0builder::ConfLambdaSelection1 confLambdaSelection;
  particlecleaner::ConfLambdaCleaner1 confLambdaCleaner;
  v0histmanager::ConfLambdaBinning1 confLambdaBinning;

  o2::framework::Partition<FemtoLambdas> lambdaPartition = MAKE_LAMBDA_PARTITION(confLambdaSelection);
  o2::framework::Preslice<FemtoLambdas> perColLambdas = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoLambdasWithLabel> lambdaWithLabelPartition = MAKE_LAMBDA_PARTITION(confLambdaSelection);
  o2::framework::Preslice<FemtoLambdasWithLabel> perCollambdasWithLabel = o2::aod::femtobase::stored::fColId;

  // setup k0shorts
  v0builder::ConfK0shortSelection1 confK0shortSelection;
  particlecleaner::ConfK0shortCleaner1 confK0shortCleaner;
  v0histmanager::ConfK0shortBinning1 confK0shortBinning;

  o2::framework::Partition<FemtoK0shorts> k0shortPartition = MAKE_K0SHORT_PARTITION(confK0shortSelection);
  o2::framework::Preslice<FemtoK0shorts> perColk0shorts = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoK0shortsWithLabel> k0shortWithLabelPartition = MAKE_K0SHORT_PARTITION(confK0shortSelection);
  o2::framework::Preslice<FemtoK0shortsWithLabel> perColk0shortsWithLabel = o2::aod::femtobase::stored::fColId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::ConfPairCuts confPairCuts;

  pairbuilder::PairV0V0Builder<
    v0histmanager::PrefixLambda1,
    trackhistmanager::PrefixV01PosDaughter,
    trackhistmanager::PrefixV01NegDaughter,
    v0histmanager::PrefixLambda2,
    trackhistmanager::PrefixV02PosDaughter,
    trackhistmanager::PrefixV02NegDaughter,
    pairhistmanager::PrefixV0V0Se,
    pairhistmanager::PrefixV0V0Me,
    closepairrejection::PrefixV0V0PosSe,
    closepairrejection::PrefixV0V0NegSe,
    closepairrejection::PrefixV0V0PosMe,
    closepairrejection::PrefixV0V0NegMe,
    modes::V0::kLambda,
    modes::V0::kLambda>
    pairLambdaLambdaBuilder;

  pairbuilder::PairV0V0Builder<
    v0histmanager::PrefixK0short1,
    trackhistmanager::PrefixV01PosDaughter,
    trackhistmanager::PrefixV01NegDaughter,
    v0histmanager::PrefixK0short2,
    trackhistmanager::PrefixV02PosDaughter,
    trackhistmanager::PrefixV02NegDaughter,
    pairhistmanager::PrefixV0V0Se,
    pairhistmanager::PrefixV0V0Me,
    closepairrejection::PrefixV0V0PosSe,
    closepairrejection::PrefixV0V0NegSe,
    closepairrejection::PrefixV0V0PosMe,
    closepairrejection::PrefixV0V0NegMe,
    modes::V0::kK0short,
    modes::V0::kK0short>
    pairK0shortK0shortBuilder;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}, true};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}, true};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult, o2::aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}, true};
  pairhistmanager::ConfMixing confMixing;

  o2::framework::HistogramRegistry hRegistry{"FemtoV0V0", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  // setup cpr
  closepairrejection::ConfCprV0DaugherV0DaughterPos confCprPos;
  closepairrejection::ConfCprV0DaugherV0DaughterNeg confCprNeg;

  void init(o2::framework::InitContext&)
  {
    // TODO: implement lambda-k0short
    bool processData = doprocessLambdaLambdaSameEvent || doprocessLambdaLambdaSameEvent || doprocessK0shortK0shortSameEvent || doprocessK0shortK0shortSameEvent;
    bool processMc = doprocessLambdaLambdaSameEventMc || doprocessLambdaLambdaSameEventMc || doprocessK0shortK0shortSameEventMc || doprocessK0shortK0shortSameEventMc;

    if (processData && processMc) {
      LOG(fatal) << "Both data and mc processing is enabled. Breaking...";
    }

    bool processLambdaLambda = doprocessLambdaLambdaSameEvent || doprocessLambdaLambdaMixedEvent || doprocessLambdaLambdaSameEventMc || doprocessLambdaLambdaMixedEventMc;
    bool processK0shortK0short = doprocessK0shortK0shortSameEvent || doprocessK0shortK0shortMixedEvent || doprocessK0shortK0shortSameEventMc || doprocessK0shortK0shortMixedEventMc;

    if (processLambdaLambda && processK0shortK0short) {
      LOG(fatal) << "Both lambda-lambda and k0short-k0short processing is enabled. Breaking...";
    }

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins, confMixing.multBins.value}, true};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    // setup histograms
    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> trackHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> posDauSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> negDauSpec;
    std::map<v0histmanager::V0Hist, std::vector<o2::framework::AxisSpec>> lambdaHistSpec;
    std::map<v0histmanager::V0Hist, std::vector<o2::framework::AxisSpec>> k0shortHistSpec;
    std::map<pairhistmanager::PairHist, std::vector<o2::framework::AxisSpec>> pairV0V0HistSpec;
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpecPos = closepairrejection::makeCprHistSpecMap(confCprPos);
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpecNeg = closepairrejection::makeCprHistSpecMap(confCprNeg);

    if (processData) {
      colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
      posDauSpec = trackhistmanager::makeTrackHistSpecMap(confPosDauBinning);
      negDauSpec = trackhistmanager::makeTrackHistSpecMap(confNegDauBinning);
      if (processLambdaLambda) {
        lambdaHistSpec = v0histmanager::makeV0HistSpecMap(confLambdaBinning);
        pairV0V0HistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning);
        pairLambdaLambdaBuilder.init<modes::Mode::kAnalysis>(&hRegistry, confCollisionBinning, confLambdaSelection, confLambdaSelection, confLambdaCleaner, confLambdaCleaner, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, lambdaHistSpec, lambdaHistSpec, posDauSpec, negDauSpec, pairV0V0HistSpec, cprHistSpecPos, cprHistSpecNeg);
      }

      // setup for k0short
      if (doprocessK0shortK0shortSameEvent || doprocessK0shortK0shortMixedEvent) {
        k0shortHistSpec = v0histmanager::makeV0HistSpecMap(confK0shortBinning);
        pairV0V0HistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning);
        pairK0shortK0shortBuilder.init<modes::Mode::kAnalysis>(&hRegistry, confCollisionBinning, confK0shortSelection, confK0shortSelection, confK0shortCleaner, confK0shortCleaner, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, k0shortHistSpec, k0shortHistSpec, posDauSpec, negDauSpec, pairV0V0HistSpec, cprHistSpecPos, cprHistSpecNeg);
      }
    } else {
      colHistSpec = colhistmanager::makeColMcHistSpecMap(confCollisionBinning);
      posDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confPosDauBinning);
      negDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confNegDauBinning);
      if (processLambdaLambda) {
        lambdaHistSpec = v0histmanager::makeV0McHistSpecMap(confLambdaBinning);
        pairV0V0HistSpec = pairhistmanager::makePairMcHistSpecMap(confPairBinning);
        pairLambdaLambdaBuilder.init<modes::Mode::kAnalysis_Qa>(&hRegistry, confCollisionBinning, confLambdaSelection, confLambdaSelection, confLambdaCleaner, confLambdaCleaner, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, lambdaHistSpec, lambdaHistSpec, posDauSpec, negDauSpec, pairV0V0HistSpec, cprHistSpecPos, cprHistSpecNeg);
      }

      // setup for k0short
      if (doprocessK0shortK0shortSameEvent || doprocessK0shortK0shortMixedEvent) {
        k0shortHistSpec = v0histmanager::makeV0McHistSpecMap(confK0shortBinning);
        pairV0V0HistSpec = pairhistmanager::makePairMcHistSpecMap(confPairBinning);
        pairK0shortK0shortBuilder.init<modes::Mode::kAnalysis_Qa>(&hRegistry, confCollisionBinning, confK0shortSelection, confK0shortSelection, confK0shortCleaner, confK0shortCleaner, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, k0shortHistSpec, k0shortHistSpec, posDauSpec, negDauSpec, pairV0V0HistSpec, cprHistSpecPos, cprHistSpecNeg);
      }
    }
  };

  void processLambdaLambdaSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoLambdas const& lambdas)
  {
    pairLambdaLambdaBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, lambdas, lambdaPartition, lambdaPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaLambdaSameEvent, "Enable processing same event processing for lambda-lambda", true);

  void processLambdaLambdaSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& lambdas, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairLambdaLambdaBuilder.processSameEvent<modes::Mode::kAnalysis_Mc>(col, mcCols, tracks, lambdas, lambdaWithLabelPartition, lambdaWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaLambdaSameEventMc, "Enable processing same event processing for lambda-lambda with mc information", false);

  void processLambdaLambdaMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoLambdas const& lambdas)
  {
    pairLambdaLambdaBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, lambdas, lambdaPartition, lambdaPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaLambdaMixedEvent, "Enable processing mixed event processing for lambda-lambda", true);

  void processLambdaLambdaMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& /*lambdas*/, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairLambdaLambdaBuilder.processMixedEvent<modes::Mode::kAnalysis_Mc>(cols, mcCols, tracks, lambdaWithLabelPartition, lambdaWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaLambdaMixedEventMc, "Enable processing mixed event processing for lambda-lambda with mc information", false);

  void processK0shortK0shortSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoK0shorts const& k0shorts)
  {
    pairK0shortK0shortBuilder.processSameEvent<modes::Mode::kAnalysis>(col, tracks, k0shorts, k0shortPartition, k0shortPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processK0shortK0shortSameEvent, "Enable processing same event processing for k0short-k0short", false);

  void processK0shortK0shortSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoK0shortsWithLabel const& k0shorts, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairK0shortK0shortBuilder.processSameEvent<modes::Mode::kAnalysis_Mc>(col, mcCols, tracks, k0shorts, k0shortWithLabelPartition, k0shortWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processK0shortK0shortSameEventMc, "Enable processing same event processing for k0short-k0short with mc information", false);

  void processK0shortK0shortMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoK0shorts const& k0shorts)
  {
    pairK0shortK0shortBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, tracks, k0shorts, k0shortPartition, k0shortPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processK0shortK0shortMixedEvent, "Enable processing mixed event processing for k0short-k0short", false);

  void processK0shortK0shortMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoK0shortsWithLabel const& /*k0shorts*/, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairK0shortK0shortBuilder.processMixedEvent<modes::Mode::kAnalysis>(cols, mcCols, tracks, k0shortWithLabelPartition, k0shortWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processK0shortK0shortMixedEventMc, "Enable processing mixed event processing for k0short-k0short with mc information", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairV0V0>(cfgc),
  };
  return workflow;
}
