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

/// \file femtoTwotrackresonanceQa.cxx
/// \brief Qa task for two track resonances
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/Core/twoTrackResonanceBuilder.h"
#include "PWGCF/Femto/Core/twoTrackResonanceHistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

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
using namespace o2::analysis::femto;

struct FemtoTwotrackresonanceQa {

  // setup tables
  using Collisions = FCols;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  using Phis = o2::soa::Join<FPhis, FPhiMasks>;
  using Rho0s = o2::soa::Join<FRho0s, FRho0Masks>;
  using Kstar0s = o2::soa::Join<FKstar0s, FKstar0Masks>;
  using Tracks = o2::soa::Join<FTracks, FTrackDcas, FTrackExtras, FTrackPids>;

  SliceCache cache;

  // setup for collisions
  colhistmanager::CollisionHistManager<modes::Mode::kAnalysis_Qa> colHistManager;
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  collisionbuilder::ConfCollisionFilter collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);

  // setup for phis
  twotrackresonancebuilder::ConfPhiSelection confPhiSelection;
  Partition<Phis> phiPartition = MAKE_RESONANCE_0_PARTITON(confPhiSelection);
  Preslice<Phis> perColPhis = aod::femtobase::stored::collisionId;

  twotrackresonancehistmanager::ConfPhiBinning confPhiBinning;
  twotrackresonancehistmanager::TwoTrackResonanceHistManager<
    twotrackresonancehistmanager::PrefixPhi,
    trackhistmanager::PrefixResonancePosDaughterQa,
    trackhistmanager::PrefixResonanceNegDaughterQa,
    modes::Mode::kAnalysis_Qa,
    modes::TwoTrackResonance::kPhi>
    phiHistManager;

  // setup for rho0s
  twotrackresonancebuilder::ConfRho0Selection confRho0Selection;
  Partition<Rho0s> rho0Partition = MAKE_RESONANCE_0_PARTITON(confRho0Selection);
  Preslice<Rho0s> perColRhos = aod::femtobase::stored::collisionId;

  twotrackresonancehistmanager::ConfRho0Binning confRho0Binning;
  twotrackresonancehistmanager::TwoTrackResonanceHistManager<
    twotrackresonancehistmanager::PrefixRho,
    trackhistmanager::PrefixResonancePosDaughterQa,
    trackhistmanager::PrefixResonanceNegDaughterQa,
    modes::Mode::kAnalysis_Qa,
    modes::TwoTrackResonance::kRho0>
    rho0HistManager;

  //  setup for kstar0s
  twotrackresonancebuilder::ConfKstar0Selection confKstar0Selection;
  Partition<Kstar0s> kstar0Partition = MAKE_RESONANCE_1_PARTITON(confKstar0Selection);
  Preslice<Kstar0s> perColKstars = aod::femtobase::stored::collisionId;

  twotrackresonancehistmanager::ConfKstar0Binning confKstar0Binning;
  twotrackresonancehistmanager::TwoTrackResonanceHistManager<
    twotrackresonancehistmanager::PrefixKstar,
    trackhistmanager::PrefixResonancePosDaughterQa,
    trackhistmanager::PrefixResonanceNegDaughterQa,
    modes::Mode::kAnalysis_Qa,
    modes::TwoTrackResonance::kKstar0>
    kstar0HistManager;

  // setup for daughters
  trackhistmanager::ConfResonancePosDauBinning confPosDaughterBinning;
  trackhistmanager::ConfResonancePosDauQaBinning confPosDaughterQaBinning;
  trackhistmanager::ConfResonanceNegDauBinning confNegDaughterBinning;
  trackhistmanager::ConfResonanceNegDauQaBinning confNegDaughterQaBinning;

  HistogramRegistry hRegistry{"ResonanceQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    // create a map for histogram specs
    auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
    colHistManager.init(&hRegistry, colHistSpec);

    auto posDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confPosDaughterBinning, confPosDaughterQaBinning);
    auto negDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confNegDaughterBinning, confNegDaughterQaBinning);

    if ((doprocessPhis + doprocessRho0s + doprocessKstar0s) > 1) {
      LOG(fatal) << "Only one process can be activated";
    }

    if (doprocessPhis) {
      auto phiHistSpec = twotrackresonancehistmanager::makeTwoTrackResonanceQaHistSpecMap(confPhiBinning);
      phiHistManager.init(&hRegistry, phiHistSpec, posDaughterHistSpec, negDaughterHistSpec);
    }
    if (doprocessRho0s) {
      auto rho0HistSpec = twotrackresonancehistmanager::makeTwoTrackResonanceQaHistSpecMap(confRho0Binning);
      rho0HistManager.init(&hRegistry, rho0HistSpec, posDaughterHistSpec, negDaughterHistSpec);
    }

    if (doprocessKstar0s) {
      auto kstar0HistSpec = twotrackresonancehistmanager::makeTwoTrackResonanceQaHistSpecMap(confKstar0Binning);
      kstar0HistManager.init(&hRegistry, kstar0HistSpec, posDaughterHistSpec, negDaughterHistSpec);
    }
  };

  void processPhis(FilteredCollision const& col, Phis const& /*phis*/, Tracks const& tracks)
  {
    colHistManager.fill(col);
    auto phiSlice = phiPartition->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
    for (auto const& phi : phiSlice) {
      phiHistManager.fill(phi, tracks);
    }
  };
  PROCESS_SWITCH(FemtoTwotrackresonanceQa, processPhis, "Process Phis", true);

  void processRho0s(FilteredCollision const& col, Rho0s const& /*rho0s*/, Tracks const& tracks)
  {
    colHistManager.fill(col);
    auto rho0Slice = rho0Partition->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
    for (auto const& rho0 : rho0Slice) {
      rho0HistManager.fill(rho0, tracks);
    }
  };
  PROCESS_SWITCH(FemtoTwotrackresonanceQa, processRho0s, "Process Rho0s", false);

  void processKstar0s(FilteredCollision const& col, Kstar0s const& /*kstar0s*/, Tracks const& tracks)
  {
    colHistManager.fill(col);
    auto kstar0Slice = kstar0Partition->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
    for (auto const& kstar0 : kstar0Slice) {
      kstar0HistManager.fill(kstar0, tracks);
    }
  };
  PROCESS_SWITCH(FemtoTwotrackresonanceQa, processKstar0s, "Process Kstar0s", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoTwotrackresonanceQa>(cfgc),
  };
  return workflow;
}
