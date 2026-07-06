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

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <string>
#include <vector>

using namespace o2::analysis::femto;

struct FemtoTwotrackresonanceQa {

  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks, o2::aod::FColPos, o2::aod::FColSphericities, o2::aod::FColMults>;
  using FemtoCollision = FemtoCollisions::iterator;

  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoPhis = o2::soa::Join<o2::aod::FPhis, o2::aod::FPhiMasks>;
  using FemtoRho0s = o2::soa::Join<o2::aod::FRho0s, o2::aod::FRho0Masks>;
  using FemtoKstar0s = o2::soa::Join<o2::aod::FKstar0s, o2::aod::FKstar0Masks>;
  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMass, o2::aod::FTrackDcas, o2::aod::FTrackExtras, o2::aod::FTrackPids>;

  o2::framework::SliceCache cache;

  // setup for collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::CollisionHistManager colHistManager;
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::ConfCollisionQaBinning confCollisionQaBinning;

  // setup for phis
  twotrackresonancebuilder::ConfPhiSelection confPhiSelection;
  o2::framework::Partition<FemtoPhis> phiPartition = MAKE_RESONANCE_0_PARTITON(confPhiSelection);
  o2::framework::Preslice<FemtoPhis> perColPhis = o2::aod::femtobase::stored::fColId;

  twotrackresonancehistmanager::ConfPhiBinning confPhiBinning;
  twotrackresonancehistmanager::TwoTrackResonanceHistManager<
    twotrackresonancehistmanager::PrefixPhi,
    trackhistmanager::PrefixResonancePosDaughterQa,
    trackhistmanager::PrefixResonanceNegDaughterQa,
    modes::TwoTrackResonance::kPhi>
    phiHistManager;

  // setup for rho0s
  twotrackresonancebuilder::ConfRho0Selection confRho0Selection;
  o2::framework::Partition<FemtoRho0s> rho0Partition = MAKE_RESONANCE_0_PARTITON(confRho0Selection);
  o2::framework::Preslice<FemtoRho0s> perColRhos = o2::aod::femtobase::stored::fColId;

  twotrackresonancehistmanager::ConfRho0Binning confRho0Binning;
  twotrackresonancehistmanager::TwoTrackResonanceHistManager<
    twotrackresonancehistmanager::PrefixRho,
    trackhistmanager::PrefixResonancePosDaughterQa,
    trackhistmanager::PrefixResonanceNegDaughterQa,
    modes::TwoTrackResonance::kRho0>
    rho0HistManager;

  //  setup for kstar0s
  twotrackresonancebuilder::ConfKstar0Selection confKstar0Selection;
  o2::framework::Partition<FemtoKstar0s> kstar0Partition = MAKE_RESONANCE_1_PARTITON(confKstar0Selection);
  o2::framework::Preslice<FemtoKstar0s> perColKstars = o2::aod::femtobase::stored::fColId;

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

  o2::framework::HistogramRegistry hRegistry{"ResonanceQA", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    // create a map for histogram specs
    auto colHistSpec = colhistmanager::makeColQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
    colHistManager.init<modes::Mode::kReco_Qa>(&hRegistry, colHistSpec, confCollisionBinning, confCollisionQaBinning);

    auto posDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confPosDaughterBinning, confPosDaughterQaBinning);
    auto negDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confNegDaughterBinning, confNegDaughterQaBinning);

    if ((doprocessPhis + doprocessRho0s + doprocessKstar0s) > 1) {
      LOG(fatal) << "Only one process can be activated";
    }

    if (doprocessPhis) {
      auto phiHistSpec = twotrackresonancehistmanager::makeTwoTrackResonanceQaHistSpecMap(confPhiBinning);
      phiHistManager.init<modes::Mode::kReco_Qa>(&hRegistry, phiHistSpec, confPhiSelection, posDaughterHistSpec, confPosDaughterQaBinning, negDaughterHistSpec, confNegDaughterQaBinning);
    }
    if (doprocessRho0s) {
      auto rho0HistSpec = twotrackresonancehistmanager::makeTwoTrackResonanceQaHistSpecMap(confRho0Binning);
      rho0HistManager.init<modes::Mode::kReco_Qa>(&hRegistry, rho0HistSpec, confRho0Selection, posDaughterHistSpec, confPosDaughterQaBinning, negDaughterHistSpec, confNegDaughterQaBinning);
    }

    if (doprocessKstar0s) {
      auto kstar0HistSpec = twotrackresonancehistmanager::makeTwoTrackResonanceQaHistSpecMap(confKstar0Binning);
      kstar0HistManager.init<modes::Mode::kReco_Qa>(&hRegistry, kstar0HistSpec, confKstar0Selection, posDaughterHistSpec, confPosDaughterQaBinning, negDaughterHistSpec, confNegDaughterQaBinning);
    }
  };

  void processPhis(FilteredFemtoCollision const& col, FemtoPhis const& /*phis*/, FemtoTracks const& tracks)
  {
    auto phiSlice = phiPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (phiSlice.size() == 0) {
      return;
    }
    colHistManager.fill<modes::Mode::kReco_Qa>(col);
    for (auto const& phi : phiSlice) {
      phiHistManager.fill<modes::Mode::kReco_Qa>(phi, tracks);
    }
  };
  PROCESS_SWITCH(FemtoTwotrackresonanceQa, processPhis, "Process Phis", true);

  void processRho0s(FilteredFemtoCollision const& col, FemtoRho0s const& /*rho0s*/, FemtoTracks const& tracks)
  {
    auto rho0Slice = rho0Partition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (rho0Slice.size() == 0) {
      return;
    }
    colHistManager.fill<modes::Mode::kReco_Qa>(col);
    for (auto const& rho0 : rho0Slice) {
      rho0HistManager.fill<modes::Mode::kReco_Qa>(rho0, tracks);
    }
  };
  PROCESS_SWITCH(FemtoTwotrackresonanceQa, processRho0s, "Process Rho0s", false);

  void processKstar0s(FilteredFemtoCollision const& col, FemtoKstar0s const& /*kstar0s*/, FemtoTracks const& tracks)
  {
    auto kstar0Slice = kstar0Partition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (kstar0Slice.size() == 0) {
      return;
    }
    colHistManager.fill<modes::Mode::kReco_Qa>(col);
    for (auto const& kstar0 : kstar0Slice) {
      kstar0HistManager.fill<modes::Mode::kReco_Qa>(kstar0, tracks);
    }
  };
  PROCESS_SWITCH(FemtoTwotrackresonanceQa, processKstar0s, "Process Kstar0s", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoTwotrackresonanceQa>(cfgc),
  };
  return workflow;
}
