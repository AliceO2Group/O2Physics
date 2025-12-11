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

/// \file zdcTaskLightIons.cxx
/// \brief Task for ZDC in light ions
/// \author chiara.oppedisano@cern.ch

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/ZDCLightIons.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/runDataProcessing.h>

#include <TH2.h>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::evsel;

using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms>;

struct ZdcTaskLightIons {

  Produces<aod::ZDCLightIons> zdcTableLI;

  // Configurable parameters
  Configurable<int> nBinsTiming{"nBinsTiming", 200, "n bins for debunching histo"};
  Configurable<bool> tdcCut{"tdcCut", true, "Flag for TDC cut"};
  Configurable<float> tdcZNmincut{"tdcZNmincut", -2.5, "Min. ZN TDC cut value"};
  Configurable<float> tdcZNmaxcut{"tdcZNmaxcut", 2.5, "Max. ZN TDC cut value"};
  //
  // Event selections
  Configurable<float> cfgEvSelVtxZ{"cfgEvSelVtxZ", 10, "Event selection: zVtx"};
  Configurable<bool> cfgEvSelSel8{"cfgEvSelSel8", true, "Event selection: sel8"};
  Configurable<bool> cfgEvSelsDoOccupancySel{"cfgEvSelsDoOccupancySel", true, "Event selection: do occupancy selection"};
  Configurable<float> cfgEvSelsMaxOccupancy{"cfgEvSelsMaxOccupancy", 10000, "Event selection: set max occupancy"};
  Configurable<bool> cfgEvSelsNoSameBunchPileupCut{"cfgEvSelsNoSameBunchPileupCut", true, "Event selection: no same bunch pileup cut"};
  Configurable<bool> cfgEvSelsIsGoodZvtxFT0vsPV{"cfgEvSelsIsGoodZvtxFT0vsPV", true, "Event selection: is good ZVTX FT0 vs PV"};
  Configurable<bool> cfgEvSelsNoCollInTimeRangeStandard{"cfgEvSelsNoCollInTimeRangeStandard", true, "Event selection: no collision in time range standard"};
  Configurable<bool> cfgEvSelsIsVertexITSTPC{"cfgEvSelsIsVertexITSTPC", true, "Event selection: is vertex ITSTPC"};
  Configurable<bool> cfgEvSelsIsGoodITSLayersAll{"cfgEvSelsIsGoodITSLayersAll", true, "Event selection: is good ITS layers all"};
  //
  HistogramRegistry registry{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  enum SelectionCriteria {
    evSel_zvtx,
    evSel_kIsTriggerTVX,
    evSel_sel8,
    evSel_occupancy,
    evSel_kNoSameBunchPileup,
    evSel_kIsGoodZvtxFT0vsPV,
    evSel_kNoCollInTimeRangeStandard,
    evSel_kNoTimeFrameBorder,
    evSel_kNoITSROFrameBorder,
    evSel_kIsGoodITSLayersAll,
    evSel_allEvents,
    nEventSelections
  };

  void init(InitContext const&)
  {
    registry.add("zdcDebunchHist", "ZN sum vs. diff; ZNA-ZNC (ns); ZNA+ZNC (ns)", {HistType::kTH2D, {{nBinsTiming, -20., 20.}, {nBinsTiming, -20., 20.}}});

    // if (doprocessALICEcoll) {
    registry.add("hEventCount", "Number of events; Cut; # of selected events ", {HistType::kTH1D, {{nEventSelections, 0, nEventSelections}}});
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_allEvents + 1, "All events");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_zvtx + 1, "vtxZ");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_sel8 + 1, "sel8");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kIsTriggerTVX + 1, "kIsTriggerTVX");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_occupancy + 1, "kOccupancy");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kNoSameBunchPileup + 1, "kNoSameBunchPileup");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kIsGoodZvtxFT0vsPV + 1, "kIsGoodZvtxFT0vsPV");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kNoCollInTimeRangeStandard + 1, "kNoCollInTimeRangeStandard");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kNoTimeFrameBorder + 1, "kNoTimeFrameBorder");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kNoITSROFrameBorder + 1, "kNoITSROFrameBorder");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kIsGoodITSLayersAll + 1, "kkIsGoodITSLayersAll");
    //}
  }

  template <typename TCollision>
  uint8_t eventSelected(TCollision collision)
  {
    uint8_t selectionBits = 0;
    bool selected = false;

    registry.fill(HIST("hEventCount"), evSel_allEvents);

    selected = std::fabs(collision.posZ()) < cfgEvSelVtxZ;
    if (selected) {
      selectionBits |= (uint8_t)(0x1u << evSel_zvtx);
      registry.fill(HIST("hEventCount"), evSel_zvtx);
    }

    selected = collision.sel8();
    if (selected) {
      selectionBits |= (uint8_t)(0x1u << evSel_sel8);
      registry.fill(HIST("hEventCount"), evSel_sel8);
    }

    selected = collision.selection_bit(kIsTriggerTVX);
    if (selected) {
      selectionBits |= (uint8_t)(0x1u << evSel_kIsTriggerTVX);
      registry.fill(HIST("hEventCount"), evSel_kIsTriggerTVX);
    }

    auto occupancy = collision.trackOccupancyInTimeRange();
    selected = occupancy <= cfgEvSelsMaxOccupancy;
    if (selected) {
      selectionBits |= (uint8_t)(0x1u << evSel_occupancy);
      registry.fill(HIST("hEventCount"), evSel_occupancy);
    }

    selected = collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup);
    if (selected) {
      selectionBits |= (uint8_t)(0x1u << evSel_kNoSameBunchPileup);
      registry.fill(HIST("hEventCount"), evSel_kNoSameBunchPileup);
    }

    selected = collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV);
    if (selected) {
      selectionBits |= (uint8_t)(0x1u << evSel_kIsGoodZvtxFT0vsPV);
      registry.fill(HIST("hEventCount"), evSel_kIsGoodZvtxFT0vsPV);
    }

    selected = collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard);
    if (selected) {
      selectionBits |= (uint8_t)(0x1u << evSel_kNoCollInTimeRangeStandard);
      registry.fill(HIST("hEventCount"), evSel_kNoCollInTimeRangeStandard);
    }

    selected = collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder);
    if (selected) {
      selectionBits |= (uint8_t)(0x1u << evSel_kNoTimeFrameBorder);
      registry.fill(HIST("hEventCount"), evSel_kNoTimeFrameBorder);
    }

    selected = collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder);
    if (selected) {
      selectionBits |= (uint8_t)(0x1u << evSel_kNoITSROFrameBorder);
      registry.fill(HIST("hEventCount"), evSel_kNoITSROFrameBorder);
    }

    selected = collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll);
    if (selected) {
      selectionBits |= (uint8_t)(0x1u << evSel_kIsGoodITSLayersAll);
      registry.fill(HIST("hEventCount"), evSel_kIsGoodITSLayersAll);
    }

    return selectionBits;
  }

  void processZDCBC(BCsRun3 const& bcs, aod::Zdcs const& /*zdcs*/)
  {
    for (const auto& bc : bcs) {
      if (bc.has_zdc()) {

        auto tdcZNA = bc.zdc().timeZNA();
        auto tdcZNC = bc.zdc().timeZNC();
        auto tdcZPA = bc.zdc().timeZPA();
        auto tdcZPC = bc.zdc().timeZPC();
        auto tdcZEM1 = bc.zdc().timeZEM1();
        auto tdcZEM2 = bc.zdc().timeZEM2();
        //
        double zna = bc.zdc().amplitudeZNA();
        double znc = bc.zdc().amplitudeZNC();
        double zpa = bc.zdc().amplitudeZPA();
        double zpc = bc.zdc().amplitudeZPC();
        double zem1 = bc.zdc().amplitudeZEM1();
        double zem2 = bc.zdc().amplitudeZEM2();
        //
        double pmcZNA = bc.zdc().energyCommonZNA();
        double pmcZNC = bc.zdc().energyCommonZNC();
        double pmcZPA = bc.zdc().energyCommonZPA();
        double pmcZPC = bc.zdc().energyCommonZPC();
        double pmqZNC[4] = {
          0,
          0,
          0,
          0,
        };
        double pmqZNA[4] = {
          0,
          0,
          0,
          0,
        };
        const int noofZNsectors = 4;
        for (int itow = 0; itow < noofZNsectors; itow++) {
          pmqZNA[itow] = (bc.zdc().energySectorZNA())[itow];
          pmqZNC[itow] = (bc.zdc().energySectorZNC())[itow];
        }

        bool isZNChit = false, isZNAhit = false;
        if (tdcCut) { // a narrow TDC window is set
          if ((tdcZNC >= tdcZNmincut) && (tdcZNC <= tdcZNmaxcut)) {
            isZNChit = true;
          }
          if ((tdcZNA >= tdcZNmincut) && (tdcZNA <= tdcZNmaxcut)) {
            isZNAhit = true;
          }
        } else { // if no window on TDC is set
          if (pmcZNC > 0.) {
            isZNChit = true;
          }
          if (pmcZNA > 0.) {
            isZNAhit = true;
          }
        }
        if (isZNChit && isZNAhit) {
          registry.get<TH2>(HIST("zdcDebunchHist"))->Fill(zna - znc, zna + znc);
        }

        zdcTableLI(tdcZNA, zna, pmcZNA, pmqZNA[0], pmqZNA[1], pmqZNA[2], pmqZNA[3],
                   tdcZNC, znc, pmcZNC, pmqZNC[0], pmqZNC[1], pmqZNC[2], pmqZNC[3],
                   tdcZPA, zpa, pmcZPA,
                   tdcZPC, zpc, pmcZPC,
                   tdcZEM1, zem1, tdcZEM2, zem2,
                   -1, -1, -1,
                   -1.,
                   -1, -1, -1,
                   bc.timestamp(),
                   -1);
      }
    }
  }
  /// name, description, function pointer, default value
  /// note that it has to be declared after the function, so that the pointer is known
  PROCESS_SWITCH(ZdcTaskLightIons, processZDCBC, "Processing ZDC 4 auto-triggered events", true);

  void processALICEcoll(ColEvSels const& cols, BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcs*/, aod::FT0s const& /*ft0s*/, aod::FV0As const& /*fv0*/)
  {
    // collision-based event selection
    for (auto const& collision : cols) {

      const auto& foundBC = collision.foundBC_as<BCsRun3>();
      uint8_t evSelection = eventSelected(collision);
      //
      // sel8
      if (!collision.sel8())
        continue;
      // vertex cut
      if (std::fabs(collision.posZ()) > cfgEvSelVtxZ)
        continue;

      auto zv = collision.posZ();
      auto centralityFT0C = collision.centFT0C();
      auto centralityFT0A = collision.centFT0A();
      auto centralityFT0M = collision.centFT0M();

      // FT0
      float multFT0A = 0.;
      float multFT0C = 0.;
      if (foundBC.has_ft0()) {
        for (auto const& amplitude : foundBC.ft0().amplitudeA()) {
          multFT0A += amplitude;
        }
        for (auto const& amplitude : foundBC.ft0().amplitudeC()) {
          multFT0C += amplitude;
        }
      }
      // FV0
      float multV0A = 0;
      if (foundBC.has_fv0a()) {
        for (auto const& amplitude : foundBC.fv0a().amplitude()) {
          multV0A += amplitude;
        }
      }

      if (foundBC.has_zdc()) {
        const auto& zdc = foundBC.zdc();

        auto tdcZNA = zdc.timeZNA();
        auto tdcZNC = zdc.timeZNC();
        auto tdcZPA = zdc.timeZPA();
        auto tdcZPC = zdc.timeZPC();
        auto tdcZEM1 = zdc.timeZEM1();
        auto tdcZEM2 = zdc.timeZEM2();
        //
        double zna = zdc.amplitudeZNA();
        double znc = zdc.amplitudeZNC();
        double zpa = zdc.amplitudeZPA();
        double zpc = zdc.amplitudeZPC();
        double zem1 = zdc.amplitudeZEM1();
        double zem2 = zdc.amplitudeZEM2();
        //
        double pmcZNA = zdc.energyCommonZNA();
        double pmcZNC = zdc.energyCommonZNC();
        double pmcZPA = zdc.energyCommonZPA();
        double pmcZPC = zdc.energyCommonZPC();
        double pmqZNC[4] = {
          0,
          0,
          0,
          0,
        };
        double pmqZNA[4] = {
          0,
          0,
          0,
          0,
        };
        const int noofZNsectors = 4;
        for (int itow = 0; itow < noofZNsectors; itow++) {
          pmqZNA[itow] = (zdc.energySectorZNA())[itow];
          pmqZNC[itow] = (zdc.energySectorZNC())[itow];
        }

        bool isZNChit = false, isZNAhit = false;
        if (tdcCut) { // a narrow TDC window is set
          if ((tdcZNC >= tdcZNmincut) && (tdcZNC <= tdcZNmaxcut)) {
            isZNChit = true;
          }
          if ((tdcZNA >= tdcZNmincut) && (tdcZNA <= tdcZNmaxcut)) {
            isZNAhit = true;
          }
        } else { // if no window on TDC is set
          if (pmcZNC > 0.) {
            isZNChit = true;
          }
          if (pmcZNA > 0.) {
            isZNAhit = true;
          }
        }
        if (isZNChit && isZNAhit) {
          registry.get<TH2>(HIST("zdcDebunchHist"))->Fill(zna - znc, zna + znc);
        }

        zdcTableLI(tdcZNA, zna, pmcZNA, pmqZNA[0], pmqZNA[1], pmqZNA[2], pmqZNA[3],
                   tdcZNC, znc, pmcZNC, pmqZNC[0], pmqZNC[1], pmqZNC[2], pmqZNC[3],
                   tdcZPA, zpa, pmcZPA,
                   tdcZPC, zpc, pmcZPC,
                   tdcZEM1, zem1, tdcZEM2, zem2,
                   multFT0A, multFT0C, multV0A,
                   zv,
                   centralityFT0C, centralityFT0A, centralityFT0M,
                   foundBC.timestamp(),
                   evSelection);
      }
    }
  }
  /// name, description, function pointer, default value
  /// note that it has to be declared after the function, so that the pointer is known
  PROCESS_SWITCH(ZdcTaskLightIons, processALICEcoll, "Processing ZDC for ALICE collisions", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) // o2-linter: disable=name/file-cpp
{
  return WorkflowSpec{
    adaptAnalysisTask<ZdcTaskLightIons>(cfgc)};
}
