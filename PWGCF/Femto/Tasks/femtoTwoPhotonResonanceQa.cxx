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

/// \file femtoTwoPhotonResonanceQa.cxx
/// \brief QA task for resonances built from two photons (pi0, eta, ...)
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/photonHistManager.h"
#include "PWGCF/Femto/Core/twoPhotonResonanceBuilder.h"
#include "PWGCF/Femto/Core/twoPhotonResonanceHistManager.h"
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

using namespace o2::analysis::femto;

struct FemtoTwoPhotonResonanceQa {

  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks, o2::aod::FColPos, o2::aod::FColSphericities, o2::aod::FColMults, o2::aod::FColCents>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoPi0s = o2::soa::Join<o2::aod::FPi0s, o2::aod::FPi0Masks>;
  using FemtoEtas = o2::soa::Join<o2::aod::FEtas, o2::aod::FEtaMasks>;
  using FemtoPhotons = o2::soa::Join<o2::aod::FPhotons, o2::aod::FPhotonMasks, o2::aod::FPhotonExtras>;

  o2::framework::SliceCache cache;

  // setup for collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::CollisionHistManager colHistManager;
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::ConfCollisionQaBinning confCollisionQaBinning;

  // setup for pi0s
  twophotonresonancebuilder::ConfPi0Selection confPi0Selection;
  o2::framework::Partition<FemtoPi0s> pi0Partition = MAKE_TWOPHOTONRESONANCE_PARTITION(confPi0Selection);
  o2::framework::Preslice<FemtoPi0s> perColPi0s = o2::aod::femtobase::stored::fColId;

  twophotonresonancehistmanager::ConfPi0Binning confPi0Binning;
  twophotonresonancehistmanager::TwoPhotonResonanceHistManager<
    twophotonresonancehistmanager::PrefixPi0,
    photonhistmanager::PrefixTwoPhotonResonanceDau1Qa,
    photonhistmanager::PrefixTwoPhotonResonanceDau2Qa,
    modes::TwoPhotonResonance::kPi0>
    pi0HistManager;

  // setup for etas
  twophotonresonancebuilder::ConfEtaSelection confEtaSelection;
  o2::framework::Partition<FemtoEtas> etaPartition = MAKE_TWOPHOTONRESONANCE_PARTITION(confEtaSelection);
  o2::framework::Preslice<FemtoEtas> perColEtas = o2::aod::femtobase::stored::fColId;

  twophotonresonancehistmanager::ConfEtaBinning confEtaBinning;
  twophotonresonancehistmanager::TwoPhotonResonanceHistManager<
    twophotonresonancehistmanager::PrefixEta,
    photonhistmanager::PrefixTwoPhotonResonanceDau1Qa,
    photonhistmanager::PrefixTwoPhotonResonanceDau2Qa,
    modes::TwoPhotonResonance::kEta>
    etaHistManager;

  // setup for daughters (shared between pi0 and eta -- unordered, no pos/neg or species asymmetry)
  photonhistmanager::ConfPhotonBinning confPhotonDauBinning;
  photonhistmanager::ConfPhotonQaBinning confPhotonDauQaBinning;

  o2::framework::HistogramRegistry hRegistry{"FemtoTwoPhotonResonanceQa", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    auto colHistSpec = colhistmanager::makeColQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
    colHistManager.init<modes::Mode::kReco_Qa>(&hRegistry, colHistSpec, confCollisionBinning, confCollisionQaBinning);

    auto dauHistSpec = photonhistmanager::makePhotonQaHistSpecMap(confPhotonDauBinning, confPhotonDauQaBinning);

    if ((static_cast<int>(doprocessPi0s) + static_cast<int>(doprocessEtas)) > 1) {
      LOG(fatal) << "Only one process can be activated";
    }

    if (doprocessPi0s) {
      auto pi0HistSpec = twophotonresonancehistmanager::makeTwoPhotonResonanceQaHistSpecMap(confPi0Binning);
      pi0HistManager.init<modes::Mode::kReco_Qa>(&hRegistry, pi0HistSpec, dauHistSpec, confPhotonDauQaBinning);
    }
    if (doprocessEtas) {
      auto etaHistSpec = twophotonresonancehistmanager::makeTwoPhotonResonanceQaHistSpecMap(confEtaBinning);
      etaHistManager.init<modes::Mode::kReco_Qa>(&hRegistry, etaHistSpec, dauHistSpec, confPhotonDauQaBinning);
    }
  }

  void processPi0s(FilteredFemtoCollision const& col, FemtoPi0s const& /*pi0s*/, FemtoPhotons const& photons)
  {
    auto pi0Slice = pi0Partition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (pi0Slice.size() == 0) {
      return;
    }
    colHistManager.fill<modes::Mode::kReco_Qa>(col);
    for (auto const& pi0 : pi0Slice) {
      pi0HistManager.fill<modes::Mode::kReco_Qa>(pi0, photons);
    }
  }
  PROCESS_SWITCH(FemtoTwoPhotonResonanceQa, processPi0s, "Process Pi0s", true);

  void processEtas(FilteredFemtoCollision const& col, FemtoEtas const& /*etas*/, FemtoPhotons const& photons)
  {
    auto etaSlice = etaPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (etaSlice.size() == 0) {
      return;
    }
    colHistManager.fill<modes::Mode::kReco_Qa>(col);
    for (auto const& eta : etaSlice) {
      etaHistManager.fill<modes::Mode::kReco_Qa>(eta, photons);
    }
  }
  PROCESS_SWITCH(FemtoTwoPhotonResonanceQa, processEtas, "Process Etas", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoTwoPhotonResonanceQa>(context),
  };
  return workflow;
}
