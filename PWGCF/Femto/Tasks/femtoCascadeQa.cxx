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

/// \file femtoCascadeQa.cxx
/// \brief Tasks for Qa of cascades
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/cascadeBuilder.h"
#include "PWGCF/Femto/Core/cascadeHistManager.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
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

struct FemtoCascadeQa {

  // setup tables
  using FemtoCollisions = o2::soa::Join<FCols, FColMasks, FColPos, FColSphericities, FColMults>;
  using FemtoCollision = FemtoCollisions::iterator;

  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoXis = o2::soa::Join<FXis, FXiMasks, FXiExtras>;
  using FemtoOmegas = o2::soa::Join<FOmegas, FOmegaMasks, FOmegaExtras>;
  using FemtoTracks = o2::soa::Join<FTracks, FTrackDcas, FTrackExtras, FTrackPids>;

  SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::CollisionHistManager colHistManager;
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::ConfCollisionQaBinning confCollisionQaBinning;

  // setup for xis
  cascadebuilder::ConfXiSelection confXiSelection;
  Partition<FemtoXis> xiPartition = MAKE_CASCADE_PARTITION(confXiSelection);
  Preslice<FemtoXis> preColXis = femtobase::stored::fColId;

  cascadehistmanager::ConfXiBinning confXiBinning;
  cascadehistmanager::ConfXiQaBinning confXiQaBinning;
  cascadehistmanager::CascadeHistManager<
    cascadehistmanager::PrefixXiQa,
    trackhistmanager::PrefixCascadeBachelorQa,
    trackhistmanager::PrefixV0PosDaughterQa,
    trackhistmanager::PrefixV0NegDaughterQa,
    modes::Cascade::kXi>
    xiHistManager;

  // setup for omegas
  cascadebuilder::ConfOmegaSelection confOmegaSelection;
  Partition<FemtoOmegas> omegaPartition = MAKE_CASCADE_PARTITION(confOmegaSelection);
  Preslice<FemtoOmegas> preColOmegas = femtobase::stored::fColId;

  cascadehistmanager::ConfOmegaBinning confOmegaBinning;
  cascadehistmanager::ConfOmegaQaBinning confOmegaQaBinning;
  cascadehistmanager::CascadeHistManager<
    cascadehistmanager::PrefixOmegaQa,
    trackhistmanager::PrefixCascadeBachelorQa,
    trackhistmanager::PrefixV0PosDaughterQa,
    trackhistmanager::PrefixV0NegDaughterQa,
    modes::Cascade::kOmega>
    omegaHistManager;

  // setup for daughters/bachelor
  trackhistmanager::ConfCascadePosDauBinning confPosDaughterBinning;
  trackhistmanager::ConfCascadePosDauQaBinning confPosDaughterQaBinning;
  trackhistmanager::ConfCascadeNegDauBinning confNegDaughterBinning;
  trackhistmanager::ConfCascadeNegDauQaBinning confNegDaughterQaBinning;
  trackhistmanager::ConfCascadeBachelorBinning confBachelorBinning;
  trackhistmanager::ConfCascadeBachelorQaBinning confBachelorQaBinning;

  HistogramRegistry hRegistry{"FemtoCascadeQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    if ((doprocessXis + doprocessOmegas) > 1) {
      LOG(fatal) << "Only one process can be activated";
    }

    auto colHistSpec = colhistmanager::makeColQaHistSpecMap(confCollisionBinning, confCollisionQaBinning);
    colHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, colHistSpec, confCollisionQaBinning);
    auto bachelorHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confBachelorBinning, confBachelorQaBinning);
    auto posDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confPosDaughterBinning, confPosDaughterQaBinning);
    auto negDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confNegDaughterBinning, confNegDaughterQaBinning);

    if (doprocessXis) {
      auto xiHistSpec = cascadehistmanager::makeCascadeQaHistSpecMap(confXiBinning, confXiQaBinning);
      xiHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, xiHistSpec, confXiSelection, confXiQaBinning, bachelorHistSpec, confBachelorQaBinning, posDaughterHistSpec, confPosDaughterQaBinning, negDaughterHistSpec, confNegDaughterQaBinning);
    }

    if (doprocessOmegas) {
      auto omegaHistSpec = cascadehistmanager::makeCascadeQaHistSpecMap(confOmegaBinning, confOmegaQaBinning);
      omegaHistManager.init<modes::Mode::kAnalysis_Qa>(&hRegistry, omegaHistSpec, confOmegaSelection, confOmegaQaBinning, bachelorHistSpec, confBachelorQaBinning, posDaughterHistSpec, confPosDaughterQaBinning, negDaughterHistSpec, confNegDaughterQaBinning);
    }
    hRegistry.print();
  };

  void processXis(FilteredFemtoCollision const& col, FemtoXis const& /*xis*/, FemtoTracks const& tracks)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa>(col);
    auto xiSlice = xiPartition->sliceByCached(femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& xi : xiSlice) {
      xiHistManager.fill<modes::Mode::kAnalysis_Qa>(xi, tracks);
    }
  }
  PROCESS_SWITCH(FemtoCascadeQa, processXis, "Process Xis", true);

  void processOmegas(FilteredFemtoCollision const& col, FemtoOmegas const& /*omegas*/, FemtoTracks const& tracks)
  {
    colHistManager.fill<modes::Mode::kAnalysis_Qa>(col);
    auto omegaSlice = omegaPartition->sliceByCached(femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& omega : omegaSlice) {
      omegaHistManager.fill<modes::Mode::kAnalysis_Qa>(omega, tracks);
    }
  }
  PROCESS_SWITCH(FemtoCascadeQa, processOmegas, "Process Omegas", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoCascadeQa>(cfgc),
  };
  return workflow;
}
