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
/// \author Kegang Xiong
/// \since  03/2026
/// \brief  Study ZDC energy observables versus centrality for Run 2 / Run 3.

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"

#include <chrono>
#include <cmath>
#include <cstdint>

using namespace o2;
using namespace o2::framework;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct flowZdcEnergy {

  struct : ConfigurableGroup{
             O2_DEFINE_CONFIGURABLE(cfgCentMin, float, 0.f, "Minimum centrality for selected events")
               O2_DEFINE_CONFIGURABLE(cfgCentMax, float, 90.f, "Maximum centrality for selected events")
                 O2_DEFINE_CONFIGURABLE(cfgVtxZ, float, 10.f, "Accepted z-vertex range")} evsel;

  ConfigurableAxis axisCent{"axisCent", {90, 0, 90}, "Centrality (%)"};
  ConfigurableAxis axisMult{"axisMult", {100, 0, 100000}, "Multiplicity"};
  ConfigurableAxis axisEnergy{"axisEnergy", {300, 0, 300}, "Energy"};
  ConfigurableAxis axisRescaledDiff{"axisRescaledDiff", {400, -1, 1}, "(EA-EC)/(EA+EC)"};

  // Event counter bins
  enum SelectionCriteria : uint8_t {
    kAllEvents = 0,
    kSeln,
    kZvtx,
    kCentrality,
    kBCHasZDC,
    kSelectedZDC,
    kNSelections
  };

  Service<ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry registry{"registry"};

  // Run 3
  using CollisionsRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>;
  using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
  // Run 2
  using CollisionsRun2 = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentRun2V0Ms>;
  using BCsRun2 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run2MatchedToBCSparse>;

  void init(InitContext const&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(
                 std::chrono::system_clock::now().time_since_epoch())
                 .count();
    ccdb->setCreatedNotAfter(now);

    registry.add("hEventCount", "Event counter;Selection;Events", {HistType::kTH1D, {{kNSelections, 0, kNSelections}}});
    auto hCount = registry.get<TH1>(HIST("hEventCount"));
    hCount->GetXaxis()->SetBinLabel(kAllEvents + 1, "All events");
    hCount->GetXaxis()->SetBinLabel(kSeln + 1, "Sel7/8");
    hCount->GetXaxis()->SetBinLabel(kZvtx + 1, "Zvtx");
    hCount->GetXaxis()->SetBinLabel(kCentrality + 1, "Centrality");
    hCount->GetXaxis()->SetBinLabel(kBCHasZDC + 1, "BC has ZDC");
    hCount->GetXaxis()->SetBinLabel(kSelectedZDC + 1, "Selected ZDC");

    registry.add("hCentrality", "", {HistType::kTH1D, {axisCent}});
    registry.add("hMultiplicity", "", {HistType::kTH1D, {axisMult}});

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

  // Helper: event selection
  template <typename TCollision>
  bool acceptEvent(TCollision const& collision, float centrality, const int runmode)
  {
    registry.fill(HIST("hEventCount"), kAllEvents);

    if (runmode == 2 && !collision.sel7()) {
      return false;
    }
    if (runmode == 3 && !collision.sel8()) {
      return false;
    }
    registry.fill(HIST("hEventCount"), kSeln);

    if (std::abs(collision.posZ()) > evsel.cfgVtxZ) {
      return false;
    }
    registry.fill(HIST("hEventCount"), kZvtx);

    if (centrality < evsel.cfgCentMin || centrality > evsel.cfgCentMax) {
      return false;
    }
    registry.fill(HIST("hEventCount"), kCentrality);

    return true;
  }

  // Helper: fill ZDC observables
  template <typename TCollision, typename TBCs>
  void fillZDCObservables(TCollision const& collision, float centrality)
  {
    const auto& foundBC = collision.template foundBC_as<TBCs>();
    if (!foundBC.has_zdc()) {
      return;
    }
    registry.fill(HIST("hEventCount"), kBCHasZDC);

    const auto& zdc = foundBC.zdc();
    if (zdc.energyCommonZNA() <= 1.f || zdc.energyCommonZNC() <= 1.f) {
      return;
    }
    registry.fill(HIST("hEventCount"), kSelectedZDC);

    const float energyCommonZNA = zdc.energyCommonZNA();
    const float energyCommonZNC = zdc.energyCommonZNC();
    const float energySectorZNA1 = zdc.energySectorZNA()[0];
    const float energySectorZNA2 = zdc.energySectorZNA()[1];
    const float energySectorZNA3 = zdc.energySectorZNA()[2];
    const float energySectorZNA4 = zdc.energySectorZNA()[3];
    const float energySectorZNC1 = zdc.energySectorZNC()[0];
    const float energySectorZNC2 = zdc.energySectorZNC()[1];
    const float energySectorZNC3 = zdc.energySectorZNC()[2];
    const float energySectorZNC4 = zdc.energySectorZNC()[3];

    const float sumEnergyZNA = energySectorZNA1 + energySectorZNA2 + energySectorZNA3 + energySectorZNA4;
    const float sumEnergyZNC = energySectorZNC1 + energySectorZNC2 + energySectorZNC3 + energySectorZNC4;

    const float commonDen = energyCommonZNA + energyCommonZNC;
    const float sumDen = sumEnergyZNA + sumEnergyZNC;

    registry.fill(HIST("hEnergyWithCent_ZNA_Common"), energyCommonZNA, centrality);
    registry.fill(HIST("hEnergyWithCent_ZNC_Common"), energyCommonZNC, centrality);
    registry.fill(HIST("hEnergyWithCent_ZNA_1"), energySectorZNA1, centrality);
    registry.fill(HIST("hEnergyWithCent_ZNA_2"), energySectorZNA2, centrality);
    registry.fill(HIST("hEnergyWithCent_ZNA_3"), energySectorZNA3, centrality);
    registry.fill(HIST("hEnergyWithCent_ZNA_4"), energySectorZNA4, centrality);
    registry.fill(HIST("hEnergyWithCent_ZNC_1"), energySectorZNC1, centrality);
    registry.fill(HIST("hEnergyWithCent_ZNC_2"), energySectorZNC2, centrality);
    registry.fill(HIST("hEnergyWithCent_ZNC_3"), energySectorZNC3, centrality);
    registry.fill(HIST("hEnergyWithCent_ZNC_4"), energySectorZNC4, centrality);
    registry.fill(HIST("hEnergyWithCent_ZNA_SumSectors"), sumEnergyZNA, centrality);
    registry.fill(HIST("hEnergyWithCent_ZNC_SumSectors"), sumEnergyZNC, centrality);

    if (commonDen > 1.e-6f) {
      registry.fill(HIST("hEnergyWithCent_RescaledDiff"), (energyCommonZNA - energyCommonZNC) / commonDen, centrality);
    }
    if (sumDen > 1.e-6f) {
      registry.fill(HIST("hEnergyWithCent_RescaledSumDiff"), (sumEnergyZNA - sumEnergyZNC) / sumDen, centrality);
    }
  }

  // Run 3 process
  void processRun3(CollisionsRun3::iterator const& collision,
                   BCsRun3 const&,
                   aod::Zdcs const&)
  {
    const float centrality = collision.centFT0C();
    const float multi = collision.multFT0C();

    if (!acceptEvent(collision, centrality, 3)) {
      return;
    }
    registry.fill(HIST("hCentrality"), centrality);
    registry.fill(HIST("hMultiplicity"), multi);

    fillZDCObservables<CollisionsRun3::iterator, BCsRun3>(collision, centrality);
  }

  // Run 2 process
  void processRun2(CollisionsRun2::iterator const& collision,
                   BCsRun2 const&,
                   aod::Zdcs const&)
  {
    const float centrality = collision.centRun2V0M();
    const float multi = collision.multFV0M();

    if (!acceptEvent(collision, centrality, 2)) {
      return;
    }
    registry.fill(HIST("hCentrality"), centrality);
    registry.fill(HIST("hMultiplicity"), multi);

    fillZDCObservables<CollisionsRun2::iterator, BCsRun2>(collision, centrality);
  }

  // Process switches
  PROCESS_SWITCH(flowZdcEnergy, processRun3, "Process Run 3 data", true);
  PROCESS_SWITCH(flowZdcEnergy, processRun2, "Process Run 2 data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<flowZdcEnergy>(cfgc)};
}
 