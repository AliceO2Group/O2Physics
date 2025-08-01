// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file v0Qa.cxx
/// \brief QA task for v0s
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/FemtoUnited/Core/collisionHistManager.h"
#include "PWGCF/FemtoUnited/Core/collisionSelection.h"
#include "PWGCF/FemtoUnited/Core/modes.h"
#include "PWGCF/FemtoUnited/Core/partitions.h"
#include "PWGCF/FemtoUnited/Core/trackHistManager.h"
#include "PWGCF/FemtoUnited/Core/v0HistManager.h"
#include "PWGCF/FemtoUnited/Core/v0Selection.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoCollisionsDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoV0sDerived.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtounited;

struct V0Qa {

  // setup for collisions
  colhistmanager::CollisionHistManager<modes::Mode::kANALYSIS_QA> colHistManager;
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  collisionselection::ConfCollisionSelection collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);

  // using Collisions = o2::soa::Join<FUCols, FUColPos, FUColMults, FUColCents>;
  using Collisions = FUCols;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  using Lambdas = o2::soa::Join<FULambdas, FULambdaMasks, FULambdaExtras>;
  using K0shorts = o2::soa::Join<FUK0shorts, FUK0shortMasks, FUK0shortExtras>;
  using Tracks = o2::soa::Join<FUTracks, FUTrackDCAs, FUTrackExtras, FUTrackPids>;

  SliceCache cache;

  // setup for lambdas
  v0selection::ConfLambdaSelection1 confLambdaSelection;

  Partition<Lambdas> lambdaPartition = MAKE_LAMBDA_PARTITION(confLambdaSelection);
  Preslice<Lambdas> perColLambdas = aod::femtobase::stored::collisionId;

  v0histmanager::ConfLambdaBinning1 confLambdaBinning;
  v0histmanager::ConfLambdaQaBinning1 confLambdaQaBinning;
  v0histmanager::V0HistManager<
    v0histmanager::PrefixLambdaQa,
    trackhistmanager::PrefixV0PosDaughterQa,
    trackhistmanager::PrefixV0NegDaughterQa,
    modes::Mode::kANALYSIS_QA,
    modes::V0::kLambda>
    lambdaHistManager;

  // setup for k0shorts
  v0selection::ConfK0shortSelection1 confK0shortSelection;

  Partition<K0shorts> k0shortPartition = MAKE_K0SHORT_PARTITION(confK0shortSelection);
  Preslice<K0shorts> perColK0shorts = aod::femtobase::stored::collisionId;

  v0histmanager::ConfK0shortBinning1 confK0shortBinning;
  v0histmanager::ConfK0shortQaBinning1 confK0shortQaBinning;
  v0histmanager::V0HistManager<
    v0histmanager::PrefixK0shortQa,
    trackhistmanager::PrefixV0PosDaughterQa,
    trackhistmanager::PrefixV0NegDaughterQa,
    modes::Mode::kANALYSIS_QA,
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
    auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
    colHistManager.init(&hRegistry, colHistSpec);

    auto posDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confV0PosDaughterBinning, confV0PosDaughterQaBinning);
    auto negDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confV0NegDaughterBinning, confV0NegDaughterQaBinning);

    if ((doprocessK0short + doprocessLambda) > 1) {
      LOG(fatal) << "Only one process can be activated";
    }

    if (doprocessLambda) {
      auto lambdaHistSpec = v0histmanager::makeV0QaHistSpecMap(confLambdaBinning, confLambdaQaBinning);
      lambdaHistManager.init(&hRegistry, lambdaHistSpec, posDaughterHistSpec, negDaughterHistSpec);
    }

    if (doprocessK0short) {
      auto k0shortHistSpec = v0histmanager::makeV0QaHistSpecMap(confK0shortBinning, confK0shortQaBinning);
      k0shortHistManager.init(&hRegistry, k0shortHistSpec, posDaughterHistSpec, negDaughterHistSpec);
    }
  };

  void processK0short(FilteredCollision const& col, K0shorts const& /*k0shorts*/, Tracks const& tracks)
  {
    colHistManager.fill(col);
    auto k0shortSlice = k0shortPartition->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
    for (auto const& k0short : k0shortSlice) {
      k0shortHistManager.fill(k0short, tracks);
    }
  }
  PROCESS_SWITCH(V0Qa, processK0short, "Process k0shorts", false);

  void processLambda(FilteredCollision const& col, Lambdas const& /*lambdas*/, Tracks const& tracks)
  {
    colHistManager.fill(col);
    auto lambdaSlice = lambdaPartition->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
    for (auto const& lambda : lambdaSlice) {
      lambdaHistManager.fill(lambda, tracks);
    }
  }
  PROCESS_SWITCH(V0Qa, processLambda, "Process lambdas", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<V0Qa>(cfgc),
  };
  return workflow;
}
