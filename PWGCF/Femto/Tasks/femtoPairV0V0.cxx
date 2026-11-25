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
#include "PWGCF/Femto/Core/collisionTagger.h"
#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/pairBuilder.h"
#include "PWGCF/Femto/Core/pairCleaner.h"
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
#include <Framework/DeviceSpec.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/RunningWorkflowInfo.h>
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

  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, o2::aod::FTrackLabels, o2::aod::FItsPurity, o2::aod::FTpcPurity, o2::aod::FTofPurity, o2::aod::FTrackDcas>;
  using FemtoLambdasWithLabel = o2::soa::Join<FemtoLambdas, o2::aod::FLambdaLabels, o2::aod::FLambdaFrac>;
  using FemtoK0shortsWithLabel = o2::soa::Join<FemtoK0shorts, o2::aod::FK0shortLabels, o2::aod::FK0shortFrac>;

  using FemtoMcParticlesWithLabel = o2::soa::Join<o2::aod::FMcParticles, o2::aod::FMcMotherLabels>;

  // reco-only + FitVars (no MC labels, but carries the opt-in purity/fraction tables)
  using FemtoTracksWithFitVars = o2::soa::Join<FemtoTracks, o2::aod::FItsPurity, o2::aod::FTpcPurity, o2::aod::FTofPurity, o2::aod::FTrackDcas>;
  using FemtoLambdasWithFitVars = o2::soa::Join<FemtoLambdas, o2::aod::FLambdaFrac>;
  using FemtoK0shortsWithFitVars = o2::soa::Join<FemtoK0shorts, o2::aod::FK0shortFrac>;

  // opt-in: carries FColTags, only needed by the processXxxMixedEventWithTagging functions below.
  using FemtoCollisionsWithTags = o2::soa::Join<FemtoCollisions, o2::aod::FColTags>;
  using FilteredFemtoCollisionsWithTags = o2::soa::Filtered<FemtoCollisionsWithTags>;

  // opt-in: carries both FColLabels and FColTags, for the processXxxMixedEventMcWithTagging
  // functions below -- same collisionFilter auto-apply reasoning as FemtoCollisionsWithTags above.
  using FemtoCollisionsWithLabelAndTags = o2::soa::Join<FemtoCollisionsWithLabel, o2::aod::FColTags>;
  using FilteredFemtoCollisionsWithLabelAndTags = o2::soa::Filtered<FemtoCollisionsWithLabelAndTags>;

  o2::framework::SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // opt-in, only needed by the processXxxMixedEventWithTagging functions below -- three separate
  // instances since lambda-lambda/k0short-k0short/lambda-k0short are mutually exclusive on one
  // instance, but tagged independently (see init()).
  collisiontagger::ConfTagging<collisiontagger::role::LambdaLambdaPair> confTaggingLambdaLambda;
  collisiontagger::ConfTagging<collisiontagger::role::K0shortK0shortPair> confTaggingK0shortK0short;
  collisiontagger::ConfTagging<collisiontagger::role::LambdaK0shortPair> confTaggingLambdaK0short;

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

  o2::framework::Partition<FemtoLambdasWithFitVars> lambdaWithFitVarsPartition1 = MAKE_LAMBDA_PARTITION(confLambdaSelection1);
  o2::framework::Partition<FemtoLambdasWithFitVars> lambdaWithFitVarsPartition2 = MAKE_LAMBDA_PARTITION(confLambdaSelection2);
  o2::framework::Preslice<FemtoLambdasWithFitVars> perColLambdasWithFitVars = o2::aod::femtobase::stored::fColId;

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

  o2::framework::Partition<FemtoK0shortsWithFitVars> k0shortWithFitVarsPartition1 = MAKE_K0SHORT_PARTITION(confK0shortSelection1);
  o2::framework::Partition<FemtoK0shortsWithFitVars> k0shortWithFitVarsPartition2 = MAKE_K0SHORT_PARTITION(confK0shortSelection2);
  o2::framework::Preslice<FemtoK0shortsWithFitVars> perColK0shortsWithFitVars = o2::aod::femtobase::stored::fColId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::ConfPairCuts confPairCuts;
  paircleaner::ConfPairCleanerBinning confPairCleaner;

  pairbuilder::PairV0V0Builder<
    v0histmanager::PrefixLambda1Se,
    trackhistmanager::PrefixV01PosDaughterSe,
    trackhistmanager::PrefixV01NegDaughterSe,
    v0histmanager::PrefixLambda2Se,
    trackhistmanager::PrefixV02PosDaughterSe,
    trackhistmanager::PrefixV02NegDaughterSe,
    v0histmanager::PrefixLambda1Me,
    trackhistmanager::PrefixV01PosDaughterMe,
    trackhistmanager::PrefixV01NegDaughterMe,
    v0histmanager::PrefixLambda2Me,
    trackhistmanager::PrefixV02PosDaughterMe,
    trackhistmanager::PrefixV02NegDaughterMe,
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
    v0histmanager::PrefixK0short1Se,
    trackhistmanager::PrefixV01PosDaughterSe,
    trackhistmanager::PrefixV01NegDaughterSe,
    v0histmanager::PrefixK0short2Se,
    trackhistmanager::PrefixV02PosDaughterSe,
    trackhistmanager::PrefixV02NegDaughterSe,
    v0histmanager::PrefixK0short1Me,
    trackhistmanager::PrefixV01PosDaughterMe,
    trackhistmanager::PrefixV01NegDaughterMe,
    v0histmanager::PrefixK0short2Me,
    trackhistmanager::PrefixV02PosDaughterMe,
    trackhistmanager::PrefixV02NegDaughterMe,
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
    v0histmanager::PrefixLambda1Se,
    trackhistmanager::PrefixV01PosDaughterSe,
    trackhistmanager::PrefixV01NegDaughterSe,
    v0histmanager::PrefixK0short2Se,
    trackhistmanager::PrefixV02PosDaughterSe,
    trackhistmanager::PrefixV02NegDaughterSe,
    v0histmanager::PrefixLambda1Me,
    trackhistmanager::PrefixV01PosDaughterMe,
    trackhistmanager::PrefixV01NegDaughterMe,
    v0histmanager::PrefixK0short2Me,
    trackhistmanager::PrefixV02PosDaughterMe,
    trackhistmanager::PrefixV02NegDaughterMe,
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
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult, o2::aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}};
  pairhistmanager::ConfMixing confMixing;

  o2::framework::HistogramRegistry hRegistry{"FemtoV0V0", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  // setup cpr
  closepairrejection::ConfCprV0DaugherV0DaughterPos confCprPos;
  closepairrejection::ConfCprV0DaugherV0DaughterNeg confCprNeg;

  void init(o2::framework::InitContext& context)
  {
    bool processData = doprocessLambdaLambdaSameEvent || doprocessLambdaLambdaMixedEvent || doprocessLambdaLambdaMixedEventWithTagging ||
                       doprocessK0shortK0shortSameEvent || doprocessK0shortK0shortMixedEvent || doprocessK0shortK0shortMixedEventWithTagging ||
                       doprocessLambdaK0shortSameEvent || doprocessLambdaK0shortMixedEvent || doprocessLambdaK0shortMixedEventWithTagging;
    bool processFitVars = doprocessLambdaLambdaSameEventWithFitVars || doprocessLambdaLambdaMixedEventWithFitVars ||
                          doprocessK0shortK0shortSameEventWithFitVars || doprocessK0shortK0shortMixedEventWithFitVars ||
                          doprocessLambdaK0shortSameEventWithFitVars || doprocessLambdaK0shortMixedEventWithFitVars;
    bool processMc = doprocessLambdaLambdaSameEventMc || doprocessLambdaLambdaMixedEventMc ||
                     doprocessK0shortK0shortSameEventMc || doprocessK0shortK0shortMixedEventMc ||
                     doprocessLambdaK0shortSameEventMc || doprocessLambdaK0shortMixedEventMc ||
                     doprocessLambdaLambdaMixedEventMcWithTagging || doprocessK0shortK0shortMixedEventMcWithTagging || doprocessLambdaK0shortMixedEventMcWithTagging;

    if (static_cast<int>(processData) + static_cast<int>(processFitVars) + static_cast<int>(processMc) > 1) {
      LOG(fatal) << "More than one of data, FitVars or mc processing is activated. Breaking...";
    }
    if (!processData && !processFitVars && !processMc) {
      LOG(fatal) << "None of data, FitVars or mc processing is activated. Breaking...";
    }

    bool processLambdaLambda = doprocessLambdaLambdaSameEvent || doprocessLambdaLambdaMixedEvent || doprocessLambdaLambdaMixedEventWithTagging ||
                               doprocessLambdaLambdaSameEventWithFitVars || doprocessLambdaLambdaMixedEventWithFitVars ||
                               doprocessLambdaLambdaSameEventMc || doprocessLambdaLambdaMixedEventMc || doprocessLambdaLambdaMixedEventMcWithTagging;
    bool processK0shortK0short = doprocessK0shortK0shortSameEvent || doprocessK0shortK0shortMixedEvent || doprocessK0shortK0shortMixedEventWithTagging ||
                                 doprocessK0shortK0shortSameEventWithFitVars || doprocessK0shortK0shortMixedEventWithFitVars ||
                                 doprocessK0shortK0shortSameEventMc || doprocessK0shortK0shortMixedEventMc || doprocessK0shortK0shortMixedEventMcWithTagging;
    bool processLambdaK0short = doprocessLambdaK0shortSameEvent || doprocessLambdaK0shortMixedEvent || doprocessLambdaK0shortMixedEventWithTagging ||
                                doprocessLambdaK0shortSameEventWithFitVars || doprocessLambdaK0shortMixedEventWithFitVars ||
                                doprocessLambdaK0shortSameEventMc || doprocessLambdaK0shortMixedEventMc || doprocessLambdaK0shortMixedEventMcWithTagging;

    if (static_cast<int>(processLambdaLambda) + static_cast<int>(processK0shortK0short) + static_cast<int>(processLambdaK0short) > 1) {
      LOG(fatal) << "Only one of lambda-lambda, k0short-k0short or lambda-k0short processing can be enabled. Breaking...";
    }

    // Each pair type discovers its own tag bit independently (see Core/collisionTagger.h).
    if (doprocessLambdaLambdaMixedEventWithTagging || doprocessK0shortK0shortMixedEventWithTagging || doprocessLambdaK0shortMixedEventWithTagging ||
        doprocessLambdaLambdaMixedEventMcWithTagging || doprocessK0shortK0shortMixedEventMcWithTagging || doprocessLambdaK0shortMixedEventMcWithTagging) {
      auto const& workflows = context.services().get<o2::framework::RunningWorkflowInfo const>();
      auto const& ownDevice = context.services().get<o2::framework::DeviceSpec const>();
      auto roles = collisiontagger::discoverTaggedRoles(workflows);
      if (doprocessLambdaLambdaMixedEventWithTagging || doprocessLambdaLambdaMixedEventMcWithTagging) {
        int bit = collisiontagger::findOwnBit(roles, ownDevice.name, confTaggingLambdaLambda.Role);
        if (bit < 0) {
          LOG(fatal) << "A processLambdaLambdaMixedEventWithTagging* function is enabled, but no femto-collision-tagger device "
                     << "(with processV0V0 enabled) tagging this device (" << ownDevice.name << ") was found in this workflow.";
        }
        confTaggingLambdaLambda.internalTagMask.value = static_cast<datatypes::CollisionTagType>(1) << bit;
        LOG(info) << "Collision tagger bit for " << ownDevice.name << " (LambdaLambda): " << bit;
      }
      if (doprocessK0shortK0shortMixedEventWithTagging || doprocessK0shortK0shortMixedEventMcWithTagging) {
        int bit = collisiontagger::findOwnBit(roles, ownDevice.name, confTaggingK0shortK0short.Role);
        if (bit < 0) {
          LOG(fatal) << "A processK0shortK0shortMixedEventWithTagging* function is enabled, but no femto-collision-tagger device "
                     << "(with processV0V0 enabled) tagging this device (" << ownDevice.name << ") was found in this workflow.";
        }
        confTaggingK0shortK0short.internalTagMask.value = static_cast<datatypes::CollisionTagType>(1) << bit;
        LOG(info) << "Collision tagger bit for " << ownDevice.name << " (K0shortK0short): " << bit;
      }
      if (doprocessLambdaK0shortMixedEventWithTagging || doprocessLambdaK0shortMixedEventMcWithTagging) {
        int bit = collisiontagger::findOwnBit(roles, ownDevice.name, confTaggingLambdaK0short.Role);
        if (bit < 0) {
          LOG(fatal) << "A processLambdaK0shortMixedEventWithTagging* function is enabled, but no femto-collision-tagger device "
                     << "(with processV0V0 enabled) tagging this device (" << ownDevice.name << ") was found in this workflow.";
        }
        confTaggingLambdaK0short.internalTagMask.value = static_cast<datatypes::CollisionTagType>(1) << bit;
        LOG(info) << "Collision tagger bit for " << ownDevice.name << " (LambdaK0short): " << bit;
      }
    }

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins.value, confMixing.multBins.value}};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}};

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
    std::map<paircleaner::PairCleanerHist, std::vector<o2::framework::AxisSpec>> pairCleanerHistSpec = paircleaner::makePairCleanerHistSpecMap(confPairCleaner);

    if (processFitVars) {
      colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
      posDauSpec1 = trackhistmanager::makeTrackHistSpecMap(confV01PosDauBinning);
      negDauSpec1 = trackhistmanager::makeTrackHistSpecMap(confV01NegDauBinning);
      posDauSpec2 = trackhistmanager::makeTrackHistSpecMap(confV02PosDauBinning);
      negDauSpec2 = trackhistmanager::makeTrackHistSpecMap(confV02NegDauBinning);
      pairV0V0HistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);

      if (processLambdaLambda) {
        lambdaHistSpec1 = v0histmanager::makeV0HistSpecMap(confLambdaBinning1);
        lambdaHistSpec2 = v0histmanager::makeV0HistSpecMap(confLambdaBinning2);
        pairLambdaLambdaBuilder.init<modes::Mode::kSe_Reco_FitVars, modes::Mode::kMe_Reco_FitVars>(&hRegistry, confCollisionBinning, confLambdaSelection1, confLambdaSelection2, confLambdaCleaner1, confLambdaCleaner2, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, lambdaHistSpec1, lambdaHistSpec2, posDauSpec1, negDauSpec1, posDauSpec2, negDauSpec2, pairV0V0HistSpec, cprHistSpecPos, cprHistSpecNeg, pairCleanerHistSpec);
      }

      if (processK0shortK0short) {
        k0shortHistSpec1 = v0histmanager::makeV0HistSpecMap(confK0shortBinning1);
        k0shortHistSpec2 = v0histmanager::makeV0HistSpecMap(confK0shortBinning2);
        pairK0shortK0shortBuilder.init<modes::Mode::kSe_Reco_FitVars, modes::Mode::kMe_Reco_FitVars>(&hRegistry, confCollisionBinning, confK0shortSelection1, confK0shortSelection2, confK0shortCleaner1, confK0shortCleaner2, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, k0shortHistSpec1, k0shortHistSpec2, posDauSpec1, negDauSpec1, posDauSpec2, negDauSpec2, pairV0V0HistSpec, cprHistSpecPos, cprHistSpecNeg, pairCleanerHistSpec);
      }

      if (processLambdaK0short) {
        lambdaHistSpec1 = v0histmanager::makeV0HistSpecMap(confLambdaBinning1);
        k0shortHistSpec2 = v0histmanager::makeV0HistSpecMap(confK0shortBinning2);
        pairLambdaK0shortBuilder.init<modes::Mode::kSe_Reco_FitVars, modes::Mode::kMe_Reco_FitVars>(&hRegistry, confCollisionBinning, confLambdaSelection1, confK0shortSelection2, confLambdaCleaner1, confK0shortCleaner2, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, lambdaHistSpec1, k0shortHistSpec2, posDauSpec1, negDauSpec1, posDauSpec2, negDauSpec2, pairV0V0HistSpec, cprHistSpecPos, cprHistSpecNeg, pairCleanerHistSpec);
      }
    } else if (processData) {
      colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
      posDauSpec1 = trackhistmanager::makeTrackHistSpecMap(confV01PosDauBinning);
      negDauSpec1 = trackhistmanager::makeTrackHistSpecMap(confV01NegDauBinning);
      posDauSpec2 = trackhistmanager::makeTrackHistSpecMap(confV02PosDauBinning);
      negDauSpec2 = trackhistmanager::makeTrackHistSpecMap(confV02NegDauBinning);
      pairV0V0HistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);

      if (processLambdaLambda) {
        lambdaHistSpec1 = v0histmanager::makeV0HistSpecMap(confLambdaBinning1);
        lambdaHistSpec2 = v0histmanager::makeV0HistSpecMap(confLambdaBinning2);
        pairLambdaLambdaBuilder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, confLambdaSelection1, confLambdaSelection2, confLambdaCleaner1, confLambdaCleaner2, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, lambdaHistSpec1, lambdaHistSpec2, posDauSpec1, negDauSpec1, posDauSpec2, negDauSpec2, pairV0V0HistSpec, cprHistSpecPos, cprHistSpecNeg, pairCleanerHistSpec);
      }

      if (processK0shortK0short) {
        k0shortHistSpec1 = v0histmanager::makeV0HistSpecMap(confK0shortBinning1);
        k0shortHistSpec2 = v0histmanager::makeV0HistSpecMap(confK0shortBinning2);
        pairK0shortK0shortBuilder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, confK0shortSelection1, confK0shortSelection2, confK0shortCleaner1, confK0shortCleaner2, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, k0shortHistSpec1, k0shortHistSpec2, posDauSpec1, negDauSpec1, posDauSpec2, negDauSpec2, pairV0V0HistSpec, cprHistSpecPos, cprHistSpecNeg, pairCleanerHistSpec);
      }

      if (processLambdaK0short) {
        lambdaHistSpec1 = v0histmanager::makeV0HistSpecMap(confLambdaBinning1);
        k0shortHistSpec2 = v0histmanager::makeV0HistSpecMap(confK0shortBinning2);
        pairLambdaK0shortBuilder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, confLambdaSelection1, confK0shortSelection2, confLambdaCleaner1, confK0shortCleaner2, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, lambdaHistSpec1, k0shortHistSpec2, posDauSpec1, negDauSpec1, posDauSpec2, negDauSpec2, pairV0V0HistSpec, cprHistSpecPos, cprHistSpecNeg, pairCleanerHistSpec);
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
        pairLambdaLambdaBuilder.init<modes::Mode::kSe_Reco_Mc_FitVars, modes::Mode::kMe_Reco_Mc_FitVars>(&hRegistry, confCollisionBinning, confLambdaSelection1, confLambdaSelection2, confLambdaCleaner1, confLambdaCleaner2, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, lambdaHistSpec1, lambdaHistSpec2, posDauSpec1, negDauSpec1, posDauSpec2, negDauSpec2, pairV0V0HistSpec, cprHistSpecPos, cprHistSpecNeg, pairCleanerHistSpec);
      }

      if (processK0shortK0short) {
        k0shortHistSpec1 = v0histmanager::makeV0McHistSpecMap(confK0shortBinning1);
        k0shortHistSpec2 = v0histmanager::makeV0McHistSpecMap(confK0shortBinning2);
        pairK0shortK0shortBuilder.init<modes::Mode::kSe_Reco_Mc_FitVars, modes::Mode::kMe_Reco_Mc_FitVars>(&hRegistry, confCollisionBinning, confK0shortSelection1, confK0shortSelection2, confK0shortCleaner1, confK0shortCleaner2, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, k0shortHistSpec1, k0shortHistSpec2, posDauSpec1, negDauSpec1, posDauSpec2, negDauSpec2, pairV0V0HistSpec, cprHistSpecPos, cprHistSpecNeg, pairCleanerHistSpec);
      }

      if (processLambdaK0short) {
        lambdaHistSpec1 = v0histmanager::makeV0McHistSpecMap(confLambdaBinning1);
        k0shortHistSpec2 = v0histmanager::makeV0McHistSpecMap(confK0shortBinning2);
        pairLambdaK0shortBuilder.init<modes::Mode::kSe_Reco_Mc_FitVars, modes::Mode::kMe_Reco_Mc_FitVars>(&hRegistry, confCollisionBinning, confLambdaSelection1, confK0shortSelection2, confLambdaCleaner1, confK0shortCleaner2, confCprPos, confCprNeg, confMixing, confPairBinning, confPairCuts, colHistSpec, lambdaHistSpec1, k0shortHistSpec2, posDauSpec1, negDauSpec1, posDauSpec2, negDauSpec2, pairV0V0HistSpec, cprHistSpecPos, cprHistSpecNeg, pairCleanerHistSpec);
      }
    }
    hRegistry.print();
  }

  void processLambdaLambdaSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoLambdas const& lambdas)
  {
    pairLambdaLambdaBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, lambdas, lambdaPartition1, lambdaPartition2, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaLambdaSameEvent, "Enable processing same event processing for lambda-lambda", true);

  void processLambdaLambdaSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& lambdas, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairLambdaLambdaBuilder.processSameEvent<modes::Mode::kSe_Reco_Mc_FitVars>(col, mcCols, tracks, lambdas, lambdaWithLabelPartition1, lambdaWithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaLambdaSameEventMc, "Enable processing same event processing for lambda-lambda with mc information", false);

  // Opt-in counterpart of processLambdaLambdaMixedEvent() below -- same "Partition declared and
  // bound inside process()" idiom as femtoPairTrackTrack.cxx's processMixedEventWithTagging.
  void processLambdaLambdaMixedEventWithTagging(FilteredFemtoCollisionsWithTags const& cols, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/)
  {
    o2::framework::Partition<FilteredFemtoCollisionsWithTags> taggedCollisions = ncheckbit(o2::aod::femtocollisions::collisionTag, confTaggingLambdaLambda.internalTagMask);
    taggedCollisions.bindTable(cols);
    pairLambdaLambdaBuilder.processMixedEvent<modes::Mode::kMe_Reco>(taggedCollisions, tracks, lambdaPartition1, lambdaPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaLambdaMixedEventWithTagging, "Enable processing mixed event processing for lambda-lambda, pre-filtered using the collision tagger", false);

  void processLambdaLambdaMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/)
  {
    pairLambdaLambdaBuilder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, lambdaPartition1, lambdaPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaLambdaMixedEvent, "Enable processing mixed event processing for lambda-lambda", true);

  void processLambdaLambdaMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& /*lambdas*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairLambdaLambdaBuilder.processMixedEvent<modes::Mode::kMe_Reco_Mc_FitVars>(cols, mcCols, tracks, lambdaWithLabelPartition1, lambdaWithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaLambdaMixedEventMc, "Enable processing mixed event processing for lambda-lambda with mc information", false);

  void processLambdaLambdaSameEventWithFitVars(FilteredFemtoCollision const& col, FemtoTracksWithFitVars const& tracks, FemtoLambdasWithFitVars const& lambdas)
  {
    pairLambdaLambdaBuilder.processSameEvent<modes::Mode::kSe_Reco_FitVars>(col, tracks, lambdas, lambdaWithFitVarsPartition1, lambdaWithFitVarsPartition2, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaLambdaSameEventWithFitVars, "Enable processing same event processing for lambda-lambda, carrying through FitVars tables", false);

  void processLambdaLambdaMixedEventWithFitVars(FilteredFemtoCollisions const& cols, FemtoTracksWithFitVars const& tracks, FemtoLambdasWithFitVars const& /*lambdas*/)
  {
    pairLambdaLambdaBuilder.processMixedEvent<modes::Mode::kMe_Reco_FitVars>(cols, tracks, lambdaWithFitVarsPartition1, lambdaWithFitVarsPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaLambdaMixedEventWithFitVars, "Enable processing mixed event processing for lambda-lambda, carrying through FitVars tables", false);

  // Same tagging filter as processLambdaLambdaMixedEventWithTagging above, on a collision table
  // that also carries FColLabels, and the MC-aware processMixedEvent overload (same as processLambdaLambdaMixedEventMc).
  void processLambdaLambdaMixedEventMcWithTagging(FilteredFemtoCollisionsWithLabelAndTags const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& /*lambdas*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    o2::framework::Partition<FilteredFemtoCollisionsWithLabelAndTags> taggedCollisions = ncheckbit(o2::aod::femtocollisions::collisionTag, confTaggingLambdaLambda.internalTagMask);
    taggedCollisions.bindTable(cols);
    pairLambdaLambdaBuilder.processMixedEvent<modes::Mode::kMe_Reco_Mc_FitVars>(taggedCollisions, mcCols, tracks, lambdaWithLabelPartition1, lambdaWithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaLambdaMixedEventMcWithTagging, "Enable processing mixed event processing for lambda-lambda, pre-filtered using the collision tagger, with mc information", false);

  void processK0shortK0shortSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoK0shorts const& k0shorts)
  {
    pairK0shortK0shortBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, k0shorts, k0shortPartition1, k0shortPartition2, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processK0shortK0shortSameEvent, "Enable processing same event processing for k0short-k0short", false);

  void processK0shortK0shortSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoK0shortsWithLabel const& k0shorts, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairK0shortK0shortBuilder.processSameEvent<modes::Mode::kSe_Reco_Mc_FitVars>(col, mcCols, tracks, k0shorts, k0shortWithLabelPartition1, k0shortWithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processK0shortK0shortSameEventMc, "Enable processing same event processing for k0short-k0short with mc information", false);

  // Opt-in counterpart of processK0shortK0shortMixedEvent() below -- same idiom as
  // processLambdaLambdaMixedEventWithTagging above.
  void processK0shortK0shortMixedEventWithTagging(FilteredFemtoCollisionsWithTags const& cols, FemtoTracks const& tracks, FemtoK0shorts const& /*k0shorts*/)
  {
    o2::framework::Partition<FilteredFemtoCollisionsWithTags> taggedCollisions = ncheckbit(o2::aod::femtocollisions::collisionTag, confTaggingK0shortK0short.internalTagMask);
    taggedCollisions.bindTable(cols);
    pairK0shortK0shortBuilder.processMixedEvent<modes::Mode::kMe_Reco>(taggedCollisions, tracks, k0shortPartition1, k0shortPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processK0shortK0shortMixedEventWithTagging, "Enable processing mixed event processing for k0short-k0short, pre-filtered using the collision tagger", false);

  void processK0shortK0shortMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoK0shorts const& /*k0shorts*/)
  {
    pairK0shortK0shortBuilder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, k0shortPartition1, k0shortPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processK0shortK0shortMixedEvent, "Enable processing mixed event processing for k0short-k0short", false);

  void processK0shortK0shortMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoK0shortsWithLabel const& /*k0shorts*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairK0shortK0shortBuilder.processMixedEvent<modes::Mode::kMe_Reco_Mc_FitVars>(cols, mcCols, tracks, k0shortWithLabelPartition1, k0shortWithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processK0shortK0shortMixedEventMc, "Enable processing mixed event processing for k0short-k0short with mc information", false);

  void processK0shortK0shortSameEventWithFitVars(FilteredFemtoCollision const& col, FemtoTracksWithFitVars const& tracks, FemtoK0shortsWithFitVars const& k0shorts)
  {
    pairK0shortK0shortBuilder.processSameEvent<modes::Mode::kSe_Reco_FitVars>(col, tracks, k0shorts, k0shortWithFitVarsPartition1, k0shortWithFitVarsPartition2, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processK0shortK0shortSameEventWithFitVars, "Enable processing same event processing for k0short-k0short, carrying through FitVars tables", false);

  void processK0shortK0shortMixedEventWithFitVars(FilteredFemtoCollisions const& cols, FemtoTracksWithFitVars const& tracks, FemtoK0shortsWithFitVars const& /*k0shorts*/)
  {
    pairK0shortK0shortBuilder.processMixedEvent<modes::Mode::kMe_Reco_FitVars>(cols, tracks, k0shortWithFitVarsPartition1, k0shortWithFitVarsPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processK0shortK0shortMixedEventWithFitVars, "Enable processing mixed event processing for k0short-k0short, carrying through FitVars tables", false);

  // Same tagging filter as processK0shortK0shortMixedEventWithTagging above, on a collision table
  // that also carries FColLabels, and the MC-aware processMixedEvent overload (same as processK0shortK0shortMixedEventMc).
  void processK0shortK0shortMixedEventMcWithTagging(FilteredFemtoCollisionsWithLabelAndTags const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoK0shortsWithLabel const& /*k0shorts*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    o2::framework::Partition<FilteredFemtoCollisionsWithLabelAndTags> taggedCollisions = ncheckbit(o2::aod::femtocollisions::collisionTag, confTaggingK0shortK0short.internalTagMask);
    taggedCollisions.bindTable(cols);
    pairK0shortK0shortBuilder.processMixedEvent<modes::Mode::kMe_Reco_Mc_FitVars>(taggedCollisions, mcCols, tracks, k0shortWithLabelPartition1, k0shortWithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processK0shortK0shortMixedEventMcWithTagging, "Enable processing mixed event processing for k0short-k0short, pre-filtered using the collision tagger, with mc information", false);

  void processLambdaK0shortSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoLambdas const& lambdas, FemtoK0shorts const& /*k0shorts*/)
  {
    pairLambdaK0shortBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, lambdas, lambdaPartition1, k0shortPartition2, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaK0shortSameEvent, "Enable processing same event processing for lambda-k0short", false);

  void processLambdaK0shortSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& lambdas, FemtoK0shortsWithLabel const& /*k0shorts*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairLambdaK0shortBuilder.processSameEvent<modes::Mode::kSe_Reco_Mc_FitVars>(col, mcCols, tracks, lambdas, lambdaWithLabelPartition1, k0shortWithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaK0shortSameEventMc, "Enable processing same event processing for lambda-k0short with mc information", false);

  // Opt-in counterpart of processLambdaK0shortMixedEvent() below -- same idiom as
  // processLambdaLambdaMixedEventWithTagging above.
  void processLambdaK0shortMixedEventWithTagging(FilteredFemtoCollisionsWithTags const& cols, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/, FemtoK0shorts const& /*k0shorts*/)
  {
    o2::framework::Partition<FilteredFemtoCollisionsWithTags> taggedCollisions = ncheckbit(o2::aod::femtocollisions::collisionTag, confTaggingLambdaK0short.internalTagMask);
    taggedCollisions.bindTable(cols);
    pairLambdaK0shortBuilder.processMixedEvent<modes::Mode::kMe_Reco>(taggedCollisions, tracks, lambdaPartition1, k0shortPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaK0shortMixedEventWithTagging, "Enable processing mixed event processing for lambda-k0short, pre-filtered using the collision tagger", false);

  void processLambdaK0shortMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/, FemtoK0shorts const& /*k0shorts*/)
  {
    pairLambdaK0shortBuilder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, lambdaPartition1, k0shortPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaK0shortMixedEvent, "Enable processing mixed event processing for lambda-k0short", false);

  void processLambdaK0shortMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& /*lambdas*/, FemtoK0shortsWithLabel const& /*k0shorts*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairLambdaK0shortBuilder.processMixedEvent<modes::Mode::kMe_Reco_Mc_FitVars>(cols, mcCols, tracks, lambdaWithLabelPartition1, k0shortWithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaK0shortMixedEventMc, "Enable processing mixed event processing for lambda-k0short with mc information", false);

  void processLambdaK0shortSameEventWithFitVars(FilteredFemtoCollision const& col, FemtoTracksWithFitVars const& tracks, FemtoLambdasWithFitVars const& lambdas, FemtoK0shortsWithFitVars const& /*k0shorts*/)
  {
    pairLambdaK0shortBuilder.processSameEvent<modes::Mode::kSe_Reco_FitVars>(col, tracks, lambdas, lambdaWithFitVarsPartition1, k0shortWithFitVarsPartition2, cache);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaK0shortSameEventWithFitVars, "Enable processing same event processing for lambda-k0short, carrying through FitVars tables", false);

  void processLambdaK0shortMixedEventWithFitVars(FilteredFemtoCollisions const& cols, FemtoTracksWithFitVars const& tracks, FemtoLambdasWithFitVars const& /*lambdas*/, FemtoK0shortsWithFitVars const& /*k0shorts*/)
  {
    pairLambdaK0shortBuilder.processMixedEvent<modes::Mode::kMe_Reco_FitVars>(cols, tracks, lambdaWithFitVarsPartition1, k0shortWithFitVarsPartition2, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaK0shortMixedEventWithFitVars, "Enable processing mixed event processing for lambda-k0short, carrying through FitVars tables", false);

  // Same tagging filter as processLambdaK0shortMixedEventWithTagging above, on a collision table
  // that also carries FColLabels, and the MC-aware processMixedEvent overload (same as processLambdaK0shortMixedEventMc).
  void processLambdaK0shortMixedEventMcWithTagging(FilteredFemtoCollisionsWithLabelAndTags const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& /*lambdas*/, FemtoK0shortsWithLabel const& /*k0shorts*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    o2::framework::Partition<FilteredFemtoCollisionsWithLabelAndTags> taggedCollisions = ncheckbit(o2::aod::femtocollisions::collisionTag, confTaggingLambdaK0short.internalTagMask);
    taggedCollisions.bindTable(cols);
    pairLambdaK0shortBuilder.processMixedEvent<modes::Mode::kMe_Reco_Mc_FitVars>(taggedCollisions, mcCols, tracks, lambdaWithLabelPartition1, k0shortWithLabelPartition2, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0V0, processLambdaK0shortMixedEventMcWithTagging, "Enable processing mixed event processing for lambda-k0short, pre-filtered using the collision tagger, with mc information", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairV0V0>(context),
  };
  return workflow;
}
