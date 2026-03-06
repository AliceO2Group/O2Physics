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

/// \file   flowZdcEnergy.cxx
/// \author Kegang Xiong (kegang.xiong@cern.ch)
/// \since  03/2026
/// \brief  A try to use the znc energy.

#include "PWGCF/DataModel/SPTableZDC.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"

#include "TF1.h"
#include "TPDGCode.h"

#include <algorithm>
#include <map>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;
// using namespace o2::analysis;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct flowZdcEnergy {

  RCTFlagsChecker rctChecker;

  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgEvtUseRCTFlagChecker, bool, false, "Evt sel: use RCT flag checker");
    O2_DEFINE_CONFIGURABLE(cfgEvtRCTFlagCheckerLabel, std::string, "CBT_hadronPID", "Evt sel: RCT flag checker label (CBT, CBT_hadronPID)"); // all Labels can be found in Common/CCDB/RCTSelectionFlags.h
    O2_DEFINE_CONFIGURABLE(cfgEvtRCTFlagCheckerZDCCheck, bool, false, "Evt sel: RCT flag checker ZDC check");
    O2_DEFINE_CONFIGURABLE(cfgEvtRCTFlagCheckerLimitAcceptAsBad, bool, false, "Evt sel: RCT flag checker treat Limited Acceptance As Bad");
  } rctFlags;

  struct : ConfigurableGroup {
    // Additional event selections
    O2_DEFINE_CONFIGURABLE(cfgMaxOccupancy, int, 10000, "Maximum occupancy of selected events");
    O2_DEFINE_CONFIGURABLE(cfgCentMin, float, 0, "Minimum cenrality for selected events");
    O2_DEFINE_CONFIGURABLE(cfgCentMax, float, 90, "Maximum cenrality for selected events");
  } EvSel;

  // Configurables containing vector
  O2_DEFINE_CONFIGURABLE(cfgVtxZ, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgEtaMax, float, 0.8f, "Maximum track #eta")
  O2_DEFINE_CONFIGURABLE(cfgPtMin, float, 0.2f, "Minimum track #P_{t}")
  O2_DEFINE_CONFIGURABLE(cfgPtMax, float, 10.0f, "Maximum track #P_{t}")
  O2_DEFINE_CONFIGURABLE(cfgDcaXYMax, float, 0.2f, "Maximum DCAxy")
  O2_DEFINE_CONFIGURABLE(cfgDcaZMax, float, 2.0f, "Maximum DCAz")

  // Filter trackFilter = nabs(aod::track::eta) < cfgEtaMax && aod::track::pt > cfgPtMin&& aod::track::pt < cfgPtMax && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && nabs(aod::track::dcaXY) < cfgDcaXYMax&& nabs(aod::track::dcaZ) < cfgDcaZMax;

  // using UsedTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;
  using ZDCCollisions = soa::Join<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>, aod::SPTableZDC>;
  using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

  //  Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    rctChecker.init(rctFlags.cfgEvtRCTFlagCheckerLabel, rctFlags.cfgEvtRCTFlagCheckerZDCCheck, rctFlags.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    AxisSpec axisCent = {100, 0, 100, "Centrality(%)"};
    AxisSpec axisEnergy = {300, 0, 300, "Energy"};
    AxisSpec axisRescaledDiff = {40, -2, 2, "(#E_{A}-#E_{C}) / (#E_{A}+#E_{C})"};

    registry.add("hEnergyWithCent_ZNA_Common", "", {HistType::kTH2D, {axisEnergy, axisCent}});
    registry.add("hEnergyWithCent_ZNC_Common", "", {HistType::kTH2D, {axisEnergy, axisCent}});
    registry.add("hEnergyWithCent_RescaledDiff", "", {HistType::kTH2D, {axisRescaledDiff, axisCent}});
    registry.add("hEnergyWithCent_ZNA_1", "", {HistType::kTH2D, {axisEnergy, axisCent}});
    registry.add("hEnergyWithCent_ZNA_2", "", {HistType::kTH2D, {axisEnergy, axisCent}});
    registry.add("hEnergyWithCent_ZNA_3", "", {HistType::kTH2D, {axisEnergy, axisCent}});
    registry.add("hEnergyWithCent_ZNA_4", "", {HistType::kTH2D, {axisEnergy, axisCent}});
    registry.add("hEnergyWithCent_ZNC_1", "", {HistType::kTH2D, {axisEnergy, axisCent}});
    registry.add("hEnergyWithCent_ZNC_2", "", {HistType::kTH2D, {axisEnergy, axisCent}});
    registry.add("hEnergyWithCent_ZNC_3", "", {HistType::kTH2D, {axisEnergy, axisCent}});
    registry.add("hEnergyWithCent_ZNC_4", "", {HistType::kTH2D, {axisEnergy, axisCent}});
    registry.add("hEnergyWithCent_ZNA_SumSectors", "", {HistType::kTH2D, {axisEnergy, axisCent}});
    registry.add("hEnergyWithCent_ZNC_SumSectors", "", {HistType::kTH2D, {axisEnergy, axisCent}});
    registry.add("hEnergyWithCent_RescaledSumDiff", "", {HistType::kTH2D, {axisRescaledDiff, axisCent}});

    
  }

  void process(ZDCCollisions::iterator const& collision, BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcs*/)
  {

    double centrality = collision.centFT0C();
    // event selection
    if (centrality < EvSel.cfgCentMin || centrality > EvSel.cfgCentMax || !collision.sel8() || std::abs(collision.posZ()) > cfgVtxZ) {
      return;
    }

    const auto& foundBC = collision.foundBC_as<BCsRun3>();
    if (!foundBC.has_zdc()) {
      return;
    }

    const auto& zdcCol = foundBC.zdc();

    double SumEnergyZNA = zdcCol.energySectorZNA()[0] + zdcCol.energySectorZNA()[1] + zdcCol.energySectorZNA()[2] + zdcCol.energySectorZNA()[3];
    double SumEnergyZNC = zdcCol.energySectorZNC()[0] + zdcCol.energySectorZNC()[1] + zdcCol.energySectorZNC()[2] + zdcCol.energySectorZNC()[3];
    double commonDen = zdcCol.energyCommonZNA() + zdcCol.energyCommonZNC();
    double sumDen = SumEnergyZNA + SumEnergyZNC;
    
    if (commonDen > 1e-3) {
      registry.fill(HIST("hEnergyWithCent_RescaledDiff"),
                    (zdcCol.energyCommonZNA() - zdcCol.energyCommonZNC()) / commonDen,
                    centrality);
    }
    if (sumDen > 1e-3) {
      registry.fill(HIST("hEnergyWithCent_RescaledSumDiff"),
                    (SumEnergyZNA - SumEnergyZNC) / sumDen,
                    centrality);
    }
    registry.fill(HIST("hEnergyWithCent_ZNA_Common"), zdcCol.energyCommonZNA(), centrality);
    registry.fill(HIST("hEnergyWithCent_ZNC_Common"), zdcCol.energyCommonZNC(), centrality);
    registry.fill(HIST("hEnergyWithCent_ZNA_1"), zdcCol.energySectorZNA()[0], centrality);
    registry.fill(HIST("hEnergyWithCent_ZNA_2"), zdcCol.energySectorZNA()[1], centrality);
    registry.fill(HIST("hEnergyWithCent_ZNA_3"), zdcCol.energySectorZNA()[2], centrality);
    registry.fill(HIST("hEnergyWithCent_ZNA_4"), zdcCol.energySectorZNA()[3], centrality);
    registry.fill(HIST("hEnergyWithCent_ZNC_1"), zdcCol.energySectorZNC()[0], centrality);
    registry.fill(HIST("hEnergyWithCent_ZNC_2"), zdcCol.energySectorZNC()[1], centrality);
    registry.fill(HIST("hEnergyWithCent_ZNC_3"), zdcCol.energySectorZNC()[2], centrality);
    registry.fill(HIST("hEnergyWithCent_ZNC_4"), zdcCol.energySectorZNC()[3], centrality);
    registry.fill(HIST("hEnergyWithCent_ZNA_SumSectors"), SumEnergyZNA, centrality);
    registry.fill(HIST("hEnergyWithCent_ZNC_SumSectors"), SumEnergyZNC, centrality);
    
  
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<flowZdcEnergy>(cfgc)};
}
