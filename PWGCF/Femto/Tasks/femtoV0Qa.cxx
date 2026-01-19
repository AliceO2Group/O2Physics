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

/// \file femtoV0Qa.cxx
/// \brief QA task for v0s
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/particleCleaner.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/Core/v0Builder.h"
#include "PWGCF/Femto/Core/v0HistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/InitContext.h"
#include "Framework/OutputObjHeader.h"
#include "Framework/runDataProcessing.h"

#include <string>
#include <vector>

using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femto;

struct FemtoV0Qa {

  // setup tables
  using FemtoCollisions = o2::soa::Join<FCols, FColMasks, FColPos, FColSphericities, FColMults>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoLambdas = o2::soa::Join<FLambdas, FLambdaMasks, FLambdaExtras>;
  using FemtoK0shorts = o2::soa::Join<FK0shorts, FK0shortMasks, FK0shortExtras>;
  using FemtoTracks = o2::soa::Join<FTracks, FTrackDcas, FTrackExtras, FTrackPids>;

  using FemtoLambdasWithLabel = o2::soa::Join<FemtoLambdas, FLambdaLabels>;
  using FemtoK0shortsWithLabel = o2::soa::Join<FemtoK0shorts, FK0shortLabels>;
  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, FTrackLabels>;

  SliceCache cache;

  // setup for collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::CollisionHistManager colHistManager;
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::ConfCollisionQaBinning confCollisionQaBinning;

  // setup for lambdas
  particlecleaner::ConfLambdaCleaner1 confLambdaCleaner;
  v0builder::ConfLambdaSelection1 confLambdaSelection;

  Partition<FemtoLambdas> lambdaPartition = MAKE_LAMBDA_PARTITION(confLambdaSelection);
  Preslice<FemtoLambdas> perColLambdas = femtobase::stored::fColId;

  Partition<FemtoLambdasWithLabel> lambdaWithLabelPartition = MAKE_LAMBDA_PARTITION(confLambdaSelection);
  Preslice<FemtoLambdasWithLabel> perColLambdasWithLabel = femtobase::stored::fColId;

  particlecleaner::ParticleCleaner lambdaCleaner;

  v0histmanager::ConfLambdaBinning1 confLambdaBinning;
  v0histmanager::ConfLambdaQaBinning1 confLambdaQaBinning;
  v0histmanager::V0HistManager<
    v0histmanager::PrefixLambdaQa,
    trackhistmanager::PrefixV0PosDaughterQa,
    trackhistmanager::PrefixV0NegDaughterQa,
    modes::V0::kLambda>
    lambdaHistManager;

  // setup for k0shorts
  particlecleaner::ConfK0shortCleaner1 confK0shortCleaner;
  v0builder::ConfK0shortSelection1 confK0shortSelection;

  Partition<FemtoK0shorts> k0shortPartition = MAKE_K0SHORT_PARTITION(confK0shortSelection);
  Preslice<FemtoK0shorts> perColK0shorts = femtobase::stored::fColId;

  Partition<FemtoK0shortsWithLabel> k0shortWithLabelPartition = MAKE_K0SHORT_PARTITION(confK0shortSelection);
  Preslice<FemtoK0shortsWithLabel> perColK0shortsWithLabel = femtobase::stored::fColId;

  particlecleaner::ParticleCleaner k0shortCleaner;

  v0histmanager::ConfK0shortBinning1 confK0shortBinning;
  v0histmanager::ConfK0shortQaBinning1 confK0shortQaBinning;
  v0histmanager::V0HistManager<
    v0histmanager::PrefixK0shortQa,
    trackhistmanager::PrefixV0PosDaughterQa,
    trackhistmanager::PrefixV0NegDaughterQa,
    modes::V0::kK0short>
    k0shortHistManager;

  // setup for daughters
  trackhistmanager::ConfV0PosDauBinning confV0PosDaughterBinning;
  trackhistmanager::ConfV0PosDauQaBinning confV0PosDaughterQaBinning;

  trackhistmanager::ConfV0NegDauBinning confV0NegDaughterBinning;
  trackhistmanager::ConfV0NegDauQaBinning confV0NegDaughterQaBinning;

  HistogramRegistry hRegistry{"FemtoV0Qa", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    if ((doprocessLambda + doprocessLambdaMc + doprocessK0short + doprocessK0shortMc) > 1) {
      LOG(fatal) << "Only one process can be activated";
    }
    bool processData = doprocessLambda || doprocessK0short;

    lambdaCleaner.init(confLambdaCleaner);
    k0shortCleaner.init(confK0shortCleaner);
    if (processData) {
      auto colHistSpec = colhistmanager::makeColQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
      colHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, colHistSpec, confCollisionQaBinning);
      auto posDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confV0PosDaughterBinning, confV0PosDaughterQaBinning);
      auto negDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confV0NegDaughterBinning, confV0NegDaughterQaBinning);
      if (doprocessLambda) {
        auto lambdaHistSpec = v0histmanager::makeV0QaHistSpecMap(confLambdaBinning, confLambdaQaBinning);
        lambdaHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, lambdaHistSpec, confLambdaSelection, confLambdaQaBinning, posDaughterHistSpec, confV0PosDaughterQaBinning, negDaughterHistSpec, confV0NegDaughterQaBinning);
      }
      if (doprocessK0short) {
        auto k0shortHistSpec = v0histmanager::makeV0QaHistSpecMap(confK0shortBinning, confK0shortQaBinning);
        k0shortHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, k0shortHistSpec, confK0shortSelection, confK0shortQaBinning, posDaughterHistSpec, confV0PosDaughterQaBinning, negDaughterHistSpec, confV0NegDaughterQaBinning);
      }
    } else {
      auto colHistSpec = colhistmanager::makeColMcQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
      colHistManager.init<modes::Mode::kAnalysis_Qa_Mc>(&hRegistry, colHistSpec, confCollisionQaBinning);
      auto posDaughterHistSpec = trackhistmanager::makeTrackMcQaHistSpecMap(confV0PosDaughterBinning, confV0PosDaughterQaBinning);
      auto negDaughterHistSpec = trackhistmanager::makeTrackMcQaHistSpecMap(confV0NegDaughterBinning, confV0NegDaughterQaBinning);
      if (doprocessLambdaMc) {
        auto lambdaHistSpec = v0histmanager::makeV0McQaHistSpecMap(confLambdaBinning, confLambdaQaBinning);
        lambdaHistManager.init<modes::Mode::kAnalysis_Qa_Mc>(&hRegistry, lambdaHistSpec, confLambdaSelection, confLambdaQaBinning, posDaughterHistSpec, confV0PosDaughterQaBinning, negDaughterHistSpec, confV0NegDaughterQaBinning);
      }
      if (doprocessK0shortMc) {
        auto k0shortHistSpec = v0histmanager::makeV0McQaHistSpecMap(confK0shortBinning, confK0shortQaBinning);
        k0shortHistManager.init<modes::Mode::kAnalysis_Qa_Mc>(&hRegistry, k0shortHistSpec, confK0shortSelection, confK0shortQaBinning, posDaughterHistSpec, confV0PosDaughterQaBinning, negDaughterHistSpec, confV0NegDaughterQaBinning);
      }
    }
    hRegistry.print();
  };

  void processK0short(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoK0shorts const& /*k0shorts*/)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa>(col);
    auto k0shortSlice = k0shortPartition->sliceByCached(femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& k0short : k0shortSlice) {
      k0shortHistManager.fill<modes::Mode::kAnalysis_Qa>(k0short, tracks);
    }
  }
  PROCESS_SWITCH(FemtoV0Qa, processK0short, "Process k0shorts", false);

  void processK0shortMc(FilteredFemtoCollisionWithLabel const& col, FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoK0shortsWithLabel const& /*k0shorts*/, FMcParticles const& mcParticles, FMcMothers const& mcMothers, FMcPartMoths const& mcPartonicMothers)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(col, mcCols);
    auto k0shortSlice = k0shortWithLabelPartition->sliceByCached(femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& k0short : k0shortSlice) {
      if (!k0shortCleaner.isClean(k0short, mcParticles, mcMothers, mcPartonicMothers)) {
        continue;
      }
      k0shortHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(k0short, tracks, mcParticles, mcMothers, mcPartonicMothers);
    }
  }
  PROCESS_SWITCH(FemtoV0Qa, processK0shortMc, "Process k0shorts with MC information", false);

  void processLambda(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa>(col);
    auto lambdaSlice = lambdaPartition->sliceByCached(femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& lambda : lambdaSlice) {
      lambdaHistManager.fill<modes::Mode::kAnalysis_Qa>(lambda, tracks);
    }
  }
  PROCESS_SWITCH(FemtoV0Qa, processLambda, "Process lambdas", true);

  void processLambdaMc(FilteredFemtoCollisionWithLabel const& col, FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& /*lambdas*/, FMcParticles const& mcParticles, FMcMothers const& mcMothers, FMcPartMoths const& mcPartonicMothers)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(col, mcCols);
    auto lambdaSlice = lambdaWithLabelPartition->sliceByCached(femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& lambda : lambdaSlice) {
      if (!lambdaCleaner.isClean(lambda, mcParticles, mcMothers, mcPartonicMothers)) {
        continue;
      }
      lambdaHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(lambda, tracks, mcParticles, mcMothers, mcPartonicMothers);
    }
  }
  PROCESS_SWITCH(FemtoV0Qa, processLambdaMc, "Process lambdas", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoV0Qa>(cfgc),
  };
  return workflow;
}
