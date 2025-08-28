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

#include <string>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/CCDB/RCTSelectionFlags.h"
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
  Configurable<float> cfgCentMin{"cfgCentMin", -1.f, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};

  // for RCT
  Configurable<bool> cfgRequireGoodRCT{"cfgRequireGoodRCT", false, "require good detector flag in run condtion table"};
  Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "CBT_hadronPID", "select 1 [CBT, CBT_hadronPID, CBT_muon_glo] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
  Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "set ZDC flag for PbPb"};
  Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};

  Configurable<float> cfgZvtxMin{"cfgZvtxMin", -1e+10, "min. Zvtx"};
  Configurable<float> cfgZvtxMax{"cfgZvtxMax", 1e+10, "max. Zvtx"};
  Configurable<bool> cfgRequireSel8{"cfgRequireSel8", false, "require sel8 in event cut"};
  Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", false, "require FT0AND in event cut"};
  Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
  Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
  Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
  Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
  Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. track occupancy"};
  Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. track occupancy"};
  Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. occupancy"};
  Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. occupancy"};
  Configurable<bool> cfgRequireNoCollInTimeRangeStandard{"cfgRequireNoCollInTimeRangeStandard", false, "require no collision in time range standard"};

  Configurable<bool> cfgRequireTVXinEMC{"cfgRequireTVXinEMC", false, "require kTVXinEMC (only for EMC analyses)"};

  o2::aod::rctsel::RCTFlagsChecker rctChecker;

  void init(InitContext&)
  {
    rctChecker.init(cfgRCTLabel.value, cfgCheckZDC.value, cfgTreatLimitedAcceptanceAsBad.value);
  }

  template <typename TCollision>
  bool isSelectedEvent(TCollision const& collision)
  {
    if constexpr (std::is_same_v<std::decay_t<TCollision>, MyCollisionsMC::iterator> || std::is_same_v<std::decay_t<TCollision>, MyCollisionsMC_Cent::iterator>) {
      if (!collision.has_mcCollision()) {
        return false;
      }
    }

    if (collision.posZ() < cfgZvtxMin || cfgZvtxMax < collision.posZ()) {
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

    if (cfgRequireTVXinEMC && !collision.alias_bit(triggerAliases::kTVXinEMC)) {
      return false;
    }

    if constexpr (std::is_same_v<std::decay_t<TCollision>, MyCollisions_Cent::iterator>) {
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        return false;
      }
    }

    if (cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
      // LOGF(info, "rejected by RCT flag");
      return false;
    }

    return true;
  }

  template <typename TCollisions>
  void processEventSelection(TCollisions const& collisions)
  {
    for (const auto& collision : collisions) {
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
