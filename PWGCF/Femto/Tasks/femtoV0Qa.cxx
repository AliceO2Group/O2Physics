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
#include <string>
#include <vector>

using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femto;

struct FemtoV0Qa {

  using FemtoCollisions = o2::soa::Join<FCols, FColMasks, FColPos, FColSphericities, FColMults>;

  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoLambdas = o2::soa::Join<FLambdas, FLambdaMasks, FLambdaExtras>;
  using FemtoK0shorts = o2::soa::Join<FK0shorts, FK0shortMasks, FK0shortExtras>;
  using FemtoTracks = o2::soa::Join<FTracks, FTrackDcas, FTrackExtras, FTrackPids>;

  SliceCache cache;

  // setup for collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::CollisionHistManager<modes::Mode::kAnalysis_Qa> colHistManager;
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::ConfCollisionQaBinning confCollisionQaBinning;

  // setup for lambdas
  v0builder::ConfLambdaSelection1 confLambdaSelection;

  Partition<FemtoLambdas> lambdaPartition = MAKE_LAMBDA_PARTITION(confLambdaSelection);
  Preslice<FemtoLambdas> perColLambdas = femtobase::stored::fColId;

  v0histmanager::ConfLambdaBinning1 confLambdaBinning;
  v0histmanager::ConfLambdaQaBinning1 confLambdaQaBinning;
  v0histmanager::V0HistManager<
    v0histmanager::PrefixLambdaQa,
    trackhistmanager::PrefixV0PosDaughterQa,
    trackhistmanager::PrefixV0NegDaughterQa,
    modes::Mode::kAnalysis_Qa,
    modes::V0::kLambda>
    lambdaHistManager;

  // setup for k0shorts
  v0builder::ConfK0shortSelection1 confK0shortSelection;

  Partition<FemtoK0shorts> k0shortPartition = MAKE_K0SHORT_PARTITION(confK0shortSelection);
  Preslice<FemtoK0shorts> perColK0shorts = femtobase::stored::fColId;

  v0histmanager::ConfK0shortBinning1 confK0shortBinning;
  v0histmanager::ConfK0shortQaBinning1 confK0shortQaBinning;
  v0histmanager::V0HistManager<
    v0histmanager::PrefixK0shortQa,
    trackhistmanager::PrefixV0PosDaughterQa,
    trackhistmanager::PrefixV0NegDaughterQa,
    modes::Mode::kAnalysis_Qa,
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
    // create a map for histogram specs
    auto colHistSpec = colhistmanager::makeColQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
    colHistManager.init(&hRegistry, colHistSpec, confCollisionQaBinning);

    auto posDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confV0PosDaughterBinning, confV0PosDaughterQaBinning);
    auto negDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confV0NegDaughterBinning, confV0NegDaughterQaBinning);

    if ((doprocessK0short + doprocessLambda) > 1) {
      LOG(fatal) << "Only one process can be activated";
    }

    if (doprocessLambda) {
      auto lambdaHistSpec = v0histmanager::makeV0QaHistSpecMap(confLambdaBinning, confLambdaQaBinning);
      lambdaHistManager.init(&hRegistry, lambdaHistSpec, confLambdaQaBinning, posDaughterHistSpec, confV0PosDaughterQaBinning, negDaughterHistSpec, confV0NegDaughterQaBinning);
    }

    if (doprocessK0short) {
      auto k0shortHistSpec = v0histmanager::makeV0QaHistSpecMap(confK0shortBinning, confK0shortQaBinning);
      k0shortHistManager.init(&hRegistry, k0shortHistSpec, confK0shortQaBinning, posDaughterHistSpec, confV0PosDaughterQaBinning, negDaughterHistSpec, confV0NegDaughterQaBinning);
    }
  };

  void processK0short(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoK0shorts const& /*k0shorts*/)
  {
    colHistManager.fill(col);
    auto k0shortSlice = k0shortPartition->sliceByCached(femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& k0short : k0shortSlice) {
      k0shortHistManager.fill(k0short, tracks);
    }
  }
  PROCESS_SWITCH(FemtoV0Qa, processK0short, "Process k0shorts", false);

  void processLambda(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/)
  {
    colHistManager.fill(col);
    auto lambdaSlice = lambdaPartition->sliceByCached(femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& lambda : lambdaSlice) {
      lambdaHistManager.fill(lambda, tracks);
    }
  }
  PROCESS_SWITCH(FemtoV0Qa, processLambda, "Process lambdas", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoV0Qa>(cfgc),
  };
  return workflow;
}
