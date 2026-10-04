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

/// \file femtoPhotonQa.cxx
/// \brief QA task for photons (PCM)
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/photonBuilder.h"
#include "PWGCF/Femto/Core/photonHistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
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

struct FemtoPhotonQa {

  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks, o2::aod::FColPos, o2::aod::FColSphericities, o2::aod::FColMults, o2::aod::FColCents>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoPhotons = o2::soa::Join<o2::aod::FPhotons, o2::aod::FPhotonMasks, o2::aod::FPhotonExtras>;

  o2::framework::SliceCache cache;

  // setup for collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::CollisionHistManager colHistManager;
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::ConfCollisionQaBinning confCollisionQaBinning;

  // setup for photons
  photonbuilder::ConfPhotonSelection1 confPhotonSelection;

  o2::framework::Partition<FemtoPhotons> photonPartition = MAKE_PHOTON_PARTITION(confPhotonSelection);
  o2::framework::Preslice<FemtoPhotons> perColPhotons = o2::aod::femtobase::stored::fColId;

  photonhistmanager::ConfPhotonBinning confPhotonBinning;
  photonhistmanager::ConfPhotonQaBinning confPhotonQaBinning;
  photonhistmanager::PhotonHistManager<photonhistmanager::PrefixPhotonQa> photonHistManager;

  o2::framework::HistogramRegistry hRegistry{"FemtoPhotonQa", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    auto colHistSpec = colhistmanager::makeColQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
    colHistManager.init<modes::Mode::kReco_Qa>(&hRegistry, colHistSpec, confCollisionBinning, confCollisionQaBinning);

    auto photonHistSpec = photonhistmanager::makePhotonQaHistSpecMap(confPhotonBinning, confPhotonQaBinning);
    photonHistManager.init<modes::Mode::kReco_Qa>(&hRegistry, photonHistSpec, confPhotonQaBinning);

    hRegistry.print();
  }

  void processPhotons(FilteredFemtoCollision const& col, FemtoPhotons const& /*photons*/)
  {
    auto photonSlice = photonPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (photonSlice.size() == 0) {
      return;
    }
    colHistManager.fill<modes::Mode::kReco_Qa>(col);
    for (auto const& photon : photonSlice) {
      photonHistManager.fill<modes::Mode::kReco_Qa>(photon);
    }
  }
  PROCESS_SWITCH(FemtoPhotonQa, processPhotons, "Process photons", true);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPhotonQa>(context),
  };
  return workflow;
}
