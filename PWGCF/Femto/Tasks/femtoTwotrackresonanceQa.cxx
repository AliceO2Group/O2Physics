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
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/Core/twoTrackResonanceBuilder.h"
#include "PWGCF/Femto/Core/twoTrackResonanceHistManager.h"
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

struct FemtoTwotrackresonanceQa {

  // setup tables
  using FemtoCollisions = o2::soa::Join<FCols, FColMasks, FColPos, FColSphericities, FColMults>;
  using FemtoCollision = FemtoCollisions::iterator;

  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoPhis = o2::soa::Join<FPhis, FPhiMasks>;
  using FemtoRho0s = o2::soa::Join<FRho0s, FRho0Masks>;
  using FemtoKstar0s = o2::soa::Join<FKstar0s, FKstar0Masks>;
  using FemtoTracks = o2::soa::Join<FTracks, FTrackDcas, FTrackExtras, FTrackPids>;

  SliceCache cache;

  // setup for collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::CollisionHistManager colHistManager;
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::ConfCollisionQaBinning confCollisionQaBinning;

  // setup for phis
  twotrackresonancebuilder::ConfPhiSelection confPhiSelection;
  Partition<FemtoPhis> phiPartition = MAKE_RESONANCE_0_PARTITON(confPhiSelection);
  Preslice<FemtoPhis> perColPhis = femtobase::stored::fColId;

  twotrackresonancehistmanager::ConfPhiBinning confPhiBinning;
  twotrackresonancehistmanager::TwoTrackResonanceHistManager<
    twotrackresonancehistmanager::PrefixPhi,
    trackhistmanager::PrefixResonancePosDaughterQa,
    trackhistmanager::PrefixResonanceNegDaughterQa,
    modes::TwoTrackResonance::kPhi>
    phiHistManager;

  // setup for rho0s
  twotrackresonancebuilder::ConfRho0Selection confRho0Selection;
  Partition<FemtoRho0s> rho0Partition = MAKE_RESONANCE_0_PARTITON(confRho0Selection);
  Preslice<FemtoRho0s> perColRhos = femtobase::stored::fColId;

  twotrackresonancehistmanager::ConfRho0Binning confRho0Binning;
  twotrackresonancehistmanager::TwoTrackResonanceHistManager<
    twotrackresonancehistmanager::PrefixRho,
    trackhistmanager::PrefixResonancePosDaughterQa,
    trackhistmanager::PrefixResonanceNegDaughterQa,
    modes::TwoTrackResonance::kRho0>
    rho0HistManager;

  //  setup for kstar0s
  twotrackresonancebuilder::ConfKstar0Selection confKstar0Selection;
  Partition<FemtoKstar0s> kstar0Partition = MAKE_RESONANCE_1_PARTITON(confKstar0Selection);
  Preslice<FemtoKstar0s> perColKstars = femtobase::stored::fColId;

  twotrackresonancehistmanager::ConfKstar0Binning confKstar0Binning;
  twotrackresonancehistmanager::TwoTrackResonanceHistManager<
    twotrackresonancehistmanager::PrefixKstar,
    trackhistmanager::PrefixResonancePosDaughterQa,
    trackhistmanager::PrefixResonanceNegDaughterQa,
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
    auto colHistSpec = colhistmanager::makeColQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
    colHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, colHistSpec, confCollisionBinning, confCollisionQaBinning);

    auto posDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confPosDaughterBinning, confPosDaughterQaBinning);
    auto negDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confNegDaughterBinning, confNegDaughterQaBinning);

    if ((doprocessPhis + doprocessRho0s + doprocessKstar0s) > 1) {
      LOG(fatal) << "Only one process can be activated";
    }

    if (doprocessPhis) {
      auto phiHistSpec = twotrackresonancehistmanager::makeTwoTrackResonanceQaHistSpecMap(confPhiBinning);
      phiHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, phiHistSpec, confPhiSelection, posDaughterHistSpec, confPosDaughterQaBinning, negDaughterHistSpec, confNegDaughterQaBinning);
    }
    if (doprocessRho0s) {
      auto rho0HistSpec = twotrackresonancehistmanager::makeTwoTrackResonanceQaHistSpecMap(confRho0Binning);
      rho0HistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, rho0HistSpec, confRho0Selection, posDaughterHistSpec, confPosDaughterQaBinning, negDaughterHistSpec, confNegDaughterQaBinning);
    }

    if (doprocessKstar0s) {
      auto kstar0HistSpec = twotrackresonancehistmanager::makeTwoTrackResonanceQaHistSpecMap(confKstar0Binning);
      kstar0HistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, kstar0HistSpec, confKstar0Selection, posDaughterHistSpec, confPosDaughterQaBinning, negDaughterHistSpec, confNegDaughterQaBinning);
    }
  };

  void processPhis(FilteredFemtoCollision const& col, FemtoPhis const& /*phis*/, FemtoTracks const& tracks)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa>(col, 0, 0, 0);
    auto phiSlice = phiPartition->sliceByCached(femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& phi : phiSlice) {
      phiHistManager.fill<modes::Mode::kAnalysis_Qa>(phi, tracks);
    }
  };
  PROCESS_SWITCH(FemtoTwotrackresonanceQa, processPhis, "Process Phis", true);

  void processRho0s(FilteredFemtoCollision const& col, FemtoRho0s const& /*rho0s*/, FemtoTracks const& tracks)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa>(col, 0, 0, 0);
    auto rho0Slice = rho0Partition->sliceByCached(femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& rho0 : rho0Slice) {
      rho0HistManager.fill<modes::Mode::kAnalysis_Qa>(rho0, tracks);
    }
  };
  PROCESS_SWITCH(FemtoTwotrackresonanceQa, processRho0s, "Process Rho0s", false);

  void processKstar0s(FilteredFemtoCollision const& col, FemtoKstar0s const& /*kstar0s*/, FemtoTracks const& tracks)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa>(col, 0, 0, 0);
    auto kstar0Slice = kstar0Partition->sliceByCached(femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& kstar0 : kstar0Slice) {
      kstar0HistManager.fill<modes::Mode::kAnalysis_Qa>(kstar0, tracks);
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
