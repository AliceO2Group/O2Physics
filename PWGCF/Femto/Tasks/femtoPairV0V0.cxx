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

// Each pair slot (1 and 2) has its own selection, cleaner, binning and partition, so that
// any combination of v0s can be correlated:
//   lambda-lambda     : LambdaSelection1 == LambdaSelection2,   Mixing.sameSpecies = true
//   lambda-antilambda : LambdaSelection1/2 differ only in sign, Mixing.sameSpecies = false
//   k0short-k0short   : K0shortSelection1 == K0shortSelection2, Mixing.sameSpecies = true
//   lambda-k0short    : LambdaSelection1 and K0shortSelection2, Mixing.sameSpecies = false
// Slot N always uses selection N, cleaner N, binning N and partition N.
struct FemtoPairV0V0 {

  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, o2::aod::FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoTracks = o2::aod::FTracks;
  using FemtoLambdas = o2::soa::Join<o2::aod::FLambdas, o2::aod::FLambdaMasks>;
  using FemtoK0shorts = o2::soa::Join<o2::aod::FK0shorts, o2::aod::FK0shortMasks>;

  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, o2::aod::FTrackLabels>;
  using FemtoLambdasWithLabel = o2::soa::Join<FemtoLambdas, o2::aod::FLambdaLabels>;
  using FemtoK0shortsWithLabel = o2::soa::Join<FemtoK0shorts, o2::aod::FK0shortLabels>;

  using FemtoMcParticlesWithLabel = o2::soa::Join<o2::aod::FMcParticles, o2::aod::FMcMotherLabels>;

  o2::framework::SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // setup for daughters
  // separate binning per v0 slot, so that e.g. lambda (proton/pion daughters) and
  // k0short (pion/pion daughters) can be binned independently
  trackhistmanager::ConfV01PosDauBinning confV01PosDauBinning;
  trackhistmanager::ConfV01NegDauBinning confV01NegDauBinning;
  trackhistmanager::ConfV02PosDauBinning confV02PosDauBinning;
  trackhistmanager::ConfV02NegDauBinning confV02NegDauBinning;

  // setup lambdas
  v0builder::ConfLambdaSelection1 confLambdaSelection1;
  v0builder::ConfLambdaSelection2 confLambdaSelection2;
  particlecleaner::ConfLambdaCleaner1 confLambdaCleaner1;
  particlecleaner::ConfLambdaCleaner2 confLambdaCleaner2;
  v0histmanager::ConfLambdaBinning1 confLambdaBinning1;
  v0histmanager::ConfLambdaBinning2 confLambdaBinning2;

  o2::framework::Partition<FemtoLambdas> lambdaPartition1 = MAKE_LAMBDA_PARTITION(confLambdaSelection1);
  o2::framework::Partition<FemtoLambdas> lambdaPartition2 = MAKE_LAMBDA_PARTITION(confLambdaSelection2);
  o2::framework::Preslice<FemtoLambdas> perColLambdas = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoLambdasWithLabel> lambdaWithLabelPartition1 = MAKE_LAMBDA_PARTITION(confLambdaSelection1);
  o2::framework::Partition<FemtoLambdasWithLabel> lambdaWithLabelPartition2 = MAKE_LAMBDA_PARTITION(confLambdaSelection2);
  o2::framework::Preslice<FemtoLambdasWithLabel> perColLambdasWithLabel = o2::aod::femtobase::stored::fColId;

  // setup k0shorts
  v0builder::ConfK0shortSelection1 confK0shortSelection1;
  v0builder::ConfK0shortSelection2 confK0shortSelection2;
  particlecleaner::ConfK0shortCleaner1 confK0shortCleaner1;
  particlecleaner::ConfK0shortCleaner2 confK0shortCleaner2;
  v0histmanager::ConfK0shortBinning1 confK0shortBinning1;
  v0histmanager::ConfK0shortBinning2 confK0shortBinning2;

  o2::framework::Partition<FemtoK0shorts> k0shortPartition1 = MAKE_K0SHORT_PARTITION(confK0shortSelection1);
  o2::framework::Partition<FemtoK0shorts> k0shortPartition2 = MAKE_K0SHORT_PARTITION(confK0shortSelection2);
  o2::framework::Preslice<FemtoK0shorts> perColK0shorts = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoK0shortsWithLabel> k0shortWithLabelPartition1 = MAKE_K0SHORT_PARTITION(confK0shortSelection1);
  o2::framework::Partition<FemtoK0shortsWithLabel> k0shortWithLabelPartition2 = MAKE_K0SHORT_PARTITION(confK0shortSelection2);
  o2::framework::Preslice<FemtoK0shortsWithLabel> perColK0shortsWithLabel = o2::aod::femtobase::stored::fColId;

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

  pairbuilder::PairV0V0Builder<
    v0histmanager::PrefixLambda1,
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
    modes::V0::kLambda,
    modes::V0::kK0short>
    pairLambdaK0shortBuilder;

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
    bool processData = doprocessLambdaLambdaSameEvent || doprocessLambdaLambdaMixedEvent ||
                       doprocessK0shortK0shortSameEvent || doprocessK0shortK0shortMixedEvent ||
                       doprocessLambdaK0shortSameEvent || doprocessLambdaK0shortMixedEvent;
    bool processMc = doprocessLambdaLambdaSameEventMc || doprocessLambdaLambdaMixedEventMc ||
                     doprocessK0shortK0shortSameEventMc || doprocessK0shortK0shortMixedEventMc ||
                     doprocessLambdaK0shortSameEventMc || doprocessLambdaK0shortMixedEventMc;

    if (processData && processMc) {
      LOG(fatal) << "Both data and mc processing is enabled. Breaking...";
    }
    if (!processData && !processMc) {
      LOG(fatal) << "Neither data nor mc processing is enabled. Breaking...";
    }

    bool processLambdaLambda = doprocessLambdaLambdaSameEvent || doprocessLambdaLambdaMixedEvent ||
                               doprocessLambdaLambdaSameEventMc || doprocessLambdaLambdaMixedEventMc;
    bool processK0shortK0short = doprocessK0shortK0shortSameEvent || doprocessK0shortK0shortMixedEvent ||
                                 doprocessK0shortK0shortSameEventMc || doprocessK0shortK0shortMixedEventMc;
    bool processLambdaK0short = doprocessLambdaK0shortSameEvent || doprocessLambdaK0shortMixedEvent ||
                                doprocessLambdaK0shortSameEventMc || doprocessLambdaK0shortMixedEventMc;

    if (static_cast<int>(processLambdaLambda) + static_cast<int>(processK0shortK0short) + static_cast<int>(processLambdaK0short) > 1) {
      LOG(fatal) << "Only one of lambda-lambda, k0short-k0short or lambda-k0short processing can be enabled. Breaking...";
    }

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins.value, confMixing.multBins.value}, true};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    // setup histograms
    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> posDauSpec1;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> negDauSpec1;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> posDauSpec2;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> negDauSpec2;
    std::map<v0histmanager::V0Hist, std::vector<o2::framework::AxisSpec>> lambdaHistSpec1;
    std::map<v0histmanager::V0Hist, std::vector<o2::framework::AxisSpec>> lambdaHistSpec2;
    std::map<v0histmanager::V0Hist, std::vector<o2::framework::AxisSpec>> k0shortHistSpec1;
    std::map<v0histmanager::V0Hist, std::vector<o2::framework::AxisSpec>> k0shortHistSpec2;
    std::map<pairhistmanager::PairHist, std::vector<o2::framework::AxisSpec>> pairV0V0HistSpec;
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpecPos = closepairrejection::makeCprHistSpecMap(confCprPos);
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpecNeg = closepairrejection::makeCprHistSpecMap(confCprNeg);

    if (processData) {
      colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
      posDauSpec1 = trackhistmanager::makeTrackHistSpecMap(confV01PosDauBinning);
      negDauSpec1 = trackhistmanager::makeTrackHistSpecMap(confV01NegDauBinning);
      posDauSpec2 = trackhistmanager::makeTrackHistSpecMap(confV02PosDauBinning);
      negDauSpec2 = trackhistmanager::makeTrackHistSpecMap(confV02NegDauBinning);
      pairV0V0HistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);

      if (processLambdaLambda) {
        lambdaHistSpec1 = v0histmanager::makeV0HistSpecMap(confLambdaBinning1);
        lambdaHistSpec2 = v0histmanager::makeV0HistSpecMap(confLambdaBinning2);
        pairLambdaLambdaBuilder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, confLambdaSelection1, confLambdaSelection2, confLambdaCleaner1, confLambdaCleaner2, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, lambdaHistSpec1, lambdaHistSpec2, posDauSpec1, negDauSpec1, posDauSpec2, negDauSpec2, pairV0V0HistSpec, cprHistSpecPos, cprHistSpecNeg);
      }

      if (processK0shortK0short) {
        k0shortHistSpec1 = v0histmanager::makeV0HistSpecMap(confK0shortBinning1);
        k0shortHistSpec2 = v0histmanager::makeV0HistSpecMap(confK0shortBinning2);
        pairK0shortK0shortBuilder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, confK0shortSelection1, confK0shortSelection2, confK0shortCleaner1, confK0shortCleaner2, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, k0shortHistSpec1, k0shortHistSpec2, posDauSpec1, negDauSpec1, posDauSpec2, negDauSpec2, pairV0V0HistSpec, cprHistSpecPos, cprHistSpecNeg);
      }

      if (processLambdaK0short) {
        lambdaHistSpec1 = v0histmanager::makeV0HistSpecMap(confLambdaBinning1);
        k0shortHistSpec2 = v0histmanager::makeV0HistSpecMap(confK0shortBinning2);
        pairLambdaK0shortBuilder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, confLambdaSelection1, confK0shortSelection2, confLambdaCleaner1, confK0shortCleaner2, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, lambdaHistSpec1, k0shortHistSpec2, posDauSpec1, negDauSpec1, posDauSpec2, negDauSpec2, pairV0V0HistSpec, cprHistSpecPos, cprHistSpecNeg);
      }
    } else {
      colHistSpec = colhistmanager::makeColMcHistSpecMap(confCollisionBinning);
      posDauSpec1 = trackhistmanager::makeTrackMcHistSpecMap(confV01PosDauBinning);
      negDauSpec1 = trackhistmanager::makeTrackMcHistSpecMap(confV01NegDauBinning);
      posDauSpec2 = trackhistmanager::makeTrackMcHistSpecMap(confV02PosDauBinning);
      negDauSpec2 = trackhistmanager::makeTrackMcHistSpecMap(confV02NegDauBinning);
      pairV0V0HistSpec = pairhistmanager::makePairMcHistSpecMap(confPairBinning, confMixing);

      if (processLambdaLambda) {
        lambdaHistSpec1 = v0histmanager::makeV0McHistSpecMap(confLambdaBinning1);
        lambdaHistSpec2 = v0histmanager::makeV0McHistSpecMap(confLambdaBinning2);
        pairLambdaLambdaBuilder.init<modes::Mode::kSe_Reco_Mc, modes::Mode::kMe_Reco_Mc>(&hRegistry, confCollisionBinning, confLambdaSelection1, confLambdaSelection2, confLambdaCleaner1, confLambdaCleaner2, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, lambdaHistSpec1, lambdaHistSpec2, posDauSpec1, negDauSpec1, posDauSpec2, negDauSpec2, pairV0V0HistSpec, cprHistSpecPos, cprHistSpecNeg);
      }

      if (processK0shortK0short) {
        k0shortHistSpec1 = v0histmanager::makeV0McHistSpecMap(confK0shortBinning1);
        k0shortHistSpec2 = v0histmanager::makeV0McHistSpecMap(confK0shortBinning2);
        pairK0shortK0shortBuilder.init<modes::Mode::kSe_Reco_Mc, modes::Mode::kMe_Reco_Mc>(&hRegistry, confCollisionBinning, confK0shortSelection1, confK0shortSelection2, confK0shortCleaner1, confK0shortCleaner2, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, k0shortHistSpec1, k0shortHistSpec2, posDauSpec1, negDauSpec1, posDauSpec2, negDauSpec2, pairV0V0HistSpec, cprHistSpecPos, cprHistSpecNeg);
      }

      if (processLambdaK0short) {
        lambdaHistSpec1 = v0histmanager::makeV0McHistSpecMap(confLambdaBinning1);
        k0shortHistSpec2 = v0histmanager::makeV0McHistSpecMap(confK0shortBinning2);
        pairLambdaK0shortBuilder.init<modes::Mode::kSe_Reco_Mc, modes::Mode::kMe_Reco_Mc>(&hRegistry, confCollisionBinning, confLambdaSelection1, confK0shortSelection2, confLambdaCleaner1, confK0shortCleaner2, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, lambdaHistSpec1, k0shortHistSpec2, posDauSpec1, negDauSpec1, posDauSpec2, negDauSpec2, pairV0V0HistSpec, cprHistSpecPos, cprHistSpecNeg);
      }
    }
  }

  void processLambdaLambdaSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoLambdas const& lambdas)
  {
    pairLambdaLambdaBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, lambdas, lambdaPartition1, lambdaPartition2, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaLambdaSameEvent, "Enable processing same event processing for lambda-lambda", true);

  void processLambdaLambdaSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& lambdas, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairLambdaLambdaBuilder.processSameEvent<modes::Mode::kSe_Reco_Mc>(col, mcCols, tracks, lambdas, lambdaWithLabelPartition1, lambdaWithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaLambdaSameEventMc, "Enable processing same event processing for lambda-lambda with mc information", false);

  void processLambdaLambdaMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/)
  {
    pairLambdaLambdaBuilder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, lambdaPartition1, lambdaPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaLambdaMixedEvent, "Enable processing mixed event processing for lambda-lambda", true);

  void processLambdaLambdaMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& /*lambdas*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairLambdaLambdaBuilder.processMixedEvent<modes::Mode::kMe_Reco_Mc>(cols, mcCols, tracks, lambdaWithLabelPartition1, lambdaWithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaLambdaMixedEventMc, "Enable processing mixed event processing for lambda-lambda with mc information", false);

  void processK0shortK0shortSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoK0shorts const& k0shorts)
  {
    pairK0shortK0shortBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, k0shorts, k0shortPartition1, k0shortPartition2, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processK0shortK0shortSameEvent, "Enable processing same event processing for k0short-k0short", false);

  void processK0shortK0shortSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoK0shortsWithLabel const& k0shorts, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairK0shortK0shortBuilder.processSameEvent<modes::Mode::kSe_Reco_Mc>(col, mcCols, tracks, k0shorts, k0shortWithLabelPartition1, k0shortWithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processK0shortK0shortSameEventMc, "Enable processing same event processing for k0short-k0short with mc information", false);

  void processK0shortK0shortMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoK0shorts const& /*k0shorts*/)
  {
    pairK0shortK0shortBuilder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, k0shortPartition1, k0shortPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processK0shortK0shortMixedEvent, "Enable processing mixed event processing for k0short-k0short", false);

  void processK0shortK0shortMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoK0shortsWithLabel const& /*k0shorts*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairK0shortK0shortBuilder.processMixedEvent<modes::Mode::kMe_Reco_Mc>(cols, mcCols, tracks, k0shortWithLabelPartition1, k0shortWithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processK0shortK0shortMixedEventMc, "Enable processing mixed event processing for k0short-k0short with mc information", false);

  void processLambdaK0shortSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoLambdas const& lambdas, FemtoK0shorts const& /*k0shorts*/)
  {
    pairLambdaK0shortBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, lambdas, lambdaPartition1, k0shortPartition2, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaK0shortSameEvent, "Enable processing same event processing for lambda-k0short", false);

  void processLambdaK0shortSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& lambdas, FemtoK0shortsWithLabel const& /*k0shorts*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairLambdaK0shortBuilder.processSameEvent<modes::Mode::kSe_Reco_Mc>(col, mcCols, tracks, lambdas, lambdaWithLabelPartition1, k0shortWithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaK0shortSameEventMc, "Enable processing same event processing for lambda-k0short with mc information", false);

  void processLambdaK0shortMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/, FemtoK0shorts const& /*k0shorts*/)
  {
    pairLambdaK0shortBuilder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, lambdaPartition1, k0shortPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaK0shortMixedEvent, "Enable processing mixed event processing for lambda-k0short", false);

  void processLambdaK0shortMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& /*lambdas*/, FemtoK0shortsWithLabel const& /*k0shorts*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairLambdaK0shortBuilder.processMixedEvent<modes::Mode::kMe_Reco_Mc>(cols, mcCols, tracks, lambdaWithLabelPartition1, k0shortWithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaK0shortMixedEventMc, "Enable processing mixed event processing for lambda-k0short with mc information", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairV0V0>(context),
  };
  return workflow;
}
