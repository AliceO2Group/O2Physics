// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
// ========================
//
// This code produces event selection table for PWG-EM.
//    Please write to: daiki.sekihata@cern.ch

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using MyCollisions_Cent = soa::Join<MyCollisions, aod::Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;

using MyCollisionsMC = soa::Join<MyCollisions, aod::McCollisionLabels>;
using MyCollisionsMC_Cent = soa::Join<MyCollisionsMC, aod::Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;

struct EMEventSelection {
  Produces<o2::aod::EMEvSels> emevsel;

  // Configurables
  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0.f, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};

  Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
  Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
  Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
  Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border in event cut"};
  Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", true, "require no ITS readout frame border in event cut"};
  Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
  Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
  Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. track occupancy"};
  Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. track occupancy"};
  Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. occupancy"};
  Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. occupancy"};
  Configurable<bool> cfgRequireNoCollInTimeRangeStandard{"cfgRequireNoCollInTimeRangeStandard", false, "require no collision in time range standard"};

  void init(InitContext&) {}

  template <typename TCollision>
  bool isSelectedEvent(TCollision const& collision)
  {
    if constexpr (std::is_same_v<std::decay_t<TCollision>, MyCollisionsMC::iterator> || std::is_same_v<std::decay_t<TCollision>, MyCollisionsMC_Cent::iterator>) {
      if (!collision.has_mcCollision()) {
        return false;
      }
    }

    if (fabs(collision.posZ()) > cfgZvtxMax) {
      return false;
    }

    if (cfgRequireFT0AND && !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      return false;
    }

    if (cfgRequireNoTFB && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }

    if (cfgRequireNoITSROFB && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }

    if (cfgRequireNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }

    if (cfgRequireGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }

    if (cfgRequireNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }

    if (!(cfgTrackOccupancyMin <= collision.trackOccupancyInTimeRange() && collision.trackOccupancyInTimeRange() < cfgTrackOccupancyMax)) {
      return false;
    }

    if (!(cfgFT0COccupancyMin <= collision.ft0cOccupancyInTimeRange() && collision.ft0cOccupancyInTimeRange() < cfgFT0COccupancyMax)) {
      return false;
    }

    if constexpr (std::is_same_v<std::decay_t<TCollision>, MyCollisions_Cent::iterator>) {
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        return false;
      }
    }

    return true;
  }

  template <typename TCollisions>
  void processEventSelection(TCollisions const& collisions)
  {
    for (auto& collision : collisions) {
      emevsel(isSelectedEvent(collision));
    } // end of collision loop
  } // end of process

  PROCESS_SWITCH_FULL(EMEventSelection, processEventSelection<MyCollisions>, processEventSelection, "event selection", true);
  PROCESS_SWITCH_FULL(EMEventSelection, processEventSelection<MyCollisions_Cent>, processEventSelection_Cent, "event selection with cent", false);
  PROCESS_SWITCH_FULL(EMEventSelection, processEventSelection<MyCollisionsMC>, processEventSelectionMC, "event selection MC", false);
  PROCESS_SWITCH_FULL(EMEventSelection, processEventSelection<MyCollisionsMC_Cent>, processEventSelectionMC_Cent, "event selection MC with cent", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<EMEventSelection>(cfgc, TaskName{"em-event-selection"})};
}
