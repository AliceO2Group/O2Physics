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

#include <map>
#include <vector>

using namespace o2::analysis::femto;

struct FemtoV0Qa {

  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks, o2::aod::FColPos, o2::aod::FColSphericities, o2::aod::FColMults>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, o2::aod::FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoLambdas = o2::soa::Join<o2::aod::FLambdas, o2::aod::FLambdaMasks, o2::aod::FLambdaExtras>;
  using FemtoK0shorts = o2::soa::Join<o2::aod::FK0shorts, o2::aod::FK0shortMasks, o2::aod::FK0shortExtras>;
  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackDcas, o2::aod::FTrackExtras, o2::aod::FTrackPids>;

  using FemtoLambdasWithLabel = o2::soa::Join<FemtoLambdas, o2::aod::FLambdaLabels>;
  using FemtoK0shortsWithLabel = o2::soa::Join<FemtoK0shorts, o2::aod::FK0shortLabels>;
  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, o2::aod::FTrackLabels>;

  o2::framework::SliceCache cache;

  // setup for collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::CollisionHistManager colHistManager;
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::ConfCollisionQaBinning confCollisionQaBinning;

  // setup for lambdas
  particlecleaner::ConfLambdaCleaner1 confLambdaCleaner;
  v0builder::ConfLambdaSelection1 confLambdaSelection;

  o2::framework::Partition<FemtoLambdas> lambdaPartition = MAKE_LAMBDA_PARTITION(confLambdaSelection);
  o2::framework::Preslice<FemtoLambdas> perColLambdas = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoLambdasWithLabel> lambdaWithLabelPartition = MAKE_LAMBDA_PARTITION(confLambdaSelection);
  o2::framework::Preslice<FemtoLambdasWithLabel> perColLambdasWithLabel = o2::aod::femtobase::stored::fColId;

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

  o2::framework::Partition<FemtoK0shorts> k0shortPartition = MAKE_K0SHORT_PARTITION(confK0shortSelection);
  o2::framework::Preslice<FemtoK0shorts> perColK0shorts = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoK0shortsWithLabel> k0shortWithLabelPartition = MAKE_K0SHORT_PARTITION(confK0shortSelection);
  o2::framework::Preslice<FemtoK0shortsWithLabel> perColK0shortsWithLabel = o2::aod::femtobase::stored::fColId;

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

  o2::framework::HistogramRegistry hRegistry{"FemtoV0Qa", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    if ((doprocessLambda + doprocessLambdaMc + doprocessK0short + doprocessK0shortMc) > 1) {
      LOG(fatal) << "Only one process can be activated";
    }
    bool processData = doprocessLambda || doprocessK0short;

    lambdaCleaner.init(confLambdaCleaner);
    k0shortCleaner.init(confK0shortCleaner);

    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> posDaughterHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> negDaughterHistSpec;
    std::map<v0histmanager::V0Hist, std::vector<o2::framework::AxisSpec>> lambdaHistSpec;
    std::map<v0histmanager::V0Hist, std::vector<o2::framework::AxisSpec>> k0shortHistSpec;

    if (processData) {
      colHistSpec = colhistmanager::makeColQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
      colHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, colHistSpec, confCollisionQaBinning);
      posDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confV0PosDaughterBinning, confV0PosDaughterQaBinning);
      negDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confV0NegDaughterBinning, confV0NegDaughterQaBinning);
      if (doprocessLambda) {
        lambdaHistSpec = v0histmanager::makeV0QaHistSpecMap(confLambdaBinning, confLambdaQaBinning);
        lambdaHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, lambdaHistSpec, confLambdaSelection, confLambdaQaBinning, posDaughterHistSpec, confV0PosDaughterQaBinning, negDaughterHistSpec, confV0NegDaughterQaBinning);
      }
      if (doprocessK0short) {
        k0shortHistSpec = v0histmanager::makeV0QaHistSpecMap(confK0shortBinning, confK0shortQaBinning);
        k0shortHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, k0shortHistSpec, confK0shortSelection, confK0shortQaBinning, posDaughterHistSpec, confV0PosDaughterQaBinning, negDaughterHistSpec, confV0NegDaughterQaBinning);
      }
    } else {
      colHistSpec = colhistmanager::makeColMcQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
      colHistManager.init<modes::Mode::kAnalysis_Qa_Mc>(&hRegistry, colHistSpec, confCollisionQaBinning);
      posDaughterHistSpec = trackhistmanager::makeTrackMcQaHistSpecMap(confV0PosDaughterBinning, confV0PosDaughterQaBinning);
      negDaughterHistSpec = trackhistmanager::makeTrackMcQaHistSpecMap(confV0NegDaughterBinning, confV0NegDaughterQaBinning);
      if (doprocessLambdaMc) {
        lambdaHistSpec = v0histmanager::makeV0McQaHistSpecMap(confLambdaBinning, confLambdaQaBinning);
        lambdaHistManager.init<modes::Mode::kAnalysis_Qa_Mc>(&hRegistry, lambdaHistSpec, confLambdaSelection, confLambdaQaBinning, posDaughterHistSpec, confV0PosDaughterQaBinning, negDaughterHistSpec, confV0NegDaughterQaBinning);
      }
      if (doprocessK0shortMc) {
        k0shortHistSpec = v0histmanager::makeV0McQaHistSpecMap(confK0shortBinning, confK0shortQaBinning);
        k0shortHistManager.init<modes::Mode::kAnalysis_Qa_Mc>(&hRegistry, k0shortHistSpec, confK0shortSelection, confK0shortQaBinning, posDaughterHistSpec, confV0PosDaughterQaBinning, negDaughterHistSpec, confV0NegDaughterQaBinning);
      }
    }
    hRegistry.print();
  };

  void processK0short(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoK0shorts const& /*k0shorts*/)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa>(col);
    auto k0shortSlice = k0shortPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& k0short : k0shortSlice) {
      k0shortHistManager.fill<modes::Mode::kAnalysis_Qa>(k0short, tracks);
    }
  }
  PROCESS_SWITCH(FemtoV0Qa, processK0short, "Process k0shorts", false);

  void processK0shortMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoK0shortsWithLabel const& /*k0shorts*/, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(col, mcCols);
    auto k0shortSlice = k0shortWithLabelPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
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
    auto lambdaSlice = lambdaPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& lambda : lambdaSlice) {
      lambdaHistManager.fill<modes::Mode::kAnalysis_Qa>(lambda, tracks);
    }
  }
  PROCESS_SWITCH(FemtoV0Qa, processLambda, "Process lambdas", true);

  void processLambdaMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& /*lambdas*/, o2::aod::FMcParticles const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(col, mcCols);
    auto lambdaSlice = lambdaWithLabelPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& lambda : lambdaSlice) {
      if (!lambdaCleaner.isClean(lambda, mcParticles, mcMothers, mcPartonicMothers)) {
        continue;
      }
      lambdaHistManager.fill<modes::Mode::kAnalysis_Qa_Mc>(lambda, tracks, mcParticles, mcMothers, mcPartonicMothers);
    }
  }
  PROCESS_SWITCH(FemtoV0Qa, processLambdaMc, "Process lambdas", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoV0Qa>(cfgc),
  };
  return workflow;
}
