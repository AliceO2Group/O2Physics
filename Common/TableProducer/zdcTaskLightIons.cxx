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

/// \file ZDCTaskLightIons.cxx
/// \brief Task for ZDC
/// \author chiara.oppedisano@cern.ch

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/ZDCLightIons.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "TH1F.h"
#include "TH2F.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::evsel;
using namespace o2::aod::track;

using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms>;

struct ZdcTaskLightIons {

  Produces<aod::ZDCLightIons> zdcTableLI;

  // Configurable parameters
  Configurable<int> nBinsTiming{"nBinsTiming", 200, "n bins for debunching histo"};
  Configurable<bool> tdcCut{"tdcCut", true, "Flag for TDC cut"};
  Configurable<float> tdcZNmincut{"tdcZNmincut", -2.5, "Min. ZN TDC cut value"};
  Configurable<float> tdcZNmaxcut{"tdcZNmaxcut", 2.5, "Max. ZN TDC cut value"};
  //
  // Event selections
  Configurable<float> etaRange{"etaRange", 1.0f, "Eta range to consider"};
  Configurable<float> cfgEvSelVtxZ{"cfgEvSelVtxZ", 10, "Event selection: zVtx"};
  Configurable<bool> cfgEvSelSel8{"cfgEvSelSel8", true, "Event selection: sel8"};
  Configurable<bool> isApplySameBunchPileup{"isApplySameBunchPileup", true, "Enable SameBunchPileup cut"};
  Configurable<bool> isApplyGoodZvtxFT0vsPV{"isApplyGoodZvtxFT0vsPV", true, "Enable GoodZvtxFT0vsPV cut"};
  Configurable<bool> isApplyVertexITSTPC{"isApplyVertexITSTPC", true, "Enable VertexITSTPC cut"};
  Configurable<float> nPVtracksCut{"nPVtracksCut", 1.0f, "Apply extra nPVtracks cut"};
  Configurable<float> cfgEvSelsMaxOccupancy{"cfgEvSelsMaxOccupancy", 10000, "Event selection: set max occupancy"};
  Configurable<bool> isApplyNoCollInTimeRangeStandard{"isApplyNoCollInTimeRangeStandard", true, "Enable NoCollInTimeRangeStandard cut"};
  Configurable<bool> isApplyNoCollInRofStandard{"isApplyNoCollInRofStandard", true, "Enable NoCollInRofStandard cut"};
  Configurable<bool> isApplyFT0CbasedOccupancy{"isApplyFT0CbasedOccupancy", true, "Enable FT0CbasedOccupancy cut"};
  //
  HistogramRegistry registry{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    if (doprocessZDCautoTrig) {
      registry.add("debunchTrue", "ZN sum vs. diff; ZNA-ZNC (ns); ZNA+ZNC (ns)", {HistType::kTH2D, {{nBinsTiming, -40., 40.}, {nBinsTiming, -40., 40.}}});
    }
    if (doprocessALICEcoll) {
      registry.add("EventHist", "EventSel; sel.; N.events", {HistType::kTH1D, {{12, 0.5, 12.5}}});
      registry.add("debunchcheck", "ZN sum vs. diff; ZNA-ZNC (ns); ZNA+ZNC (ns)", {HistType::kTH2D, {{nBinsTiming, -20., 20.}, {nBinsTiming, -20., 20.}}});
    }

    auto hstat = registry.get<TH1>(HIST("EventHist"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All events");
    x->SetBinLabel(2, "sel8");
    x->SetBinLabel(3, "vtxZ");
    x->SetBinLabel(4, "kNoSameBunchPileup");
    x->SetBinLabel(5, "kIsGoodZvtxFT0vsPV");
    x->SetBinLabel(6, "kIsVertexITSTPC");
    x->SetBinLabel(7, "kOccupancy");
    x->SetBinLabel(8, "kIsGoodITSLayersAll");
    x->SetBinLabel(9, "NoCollInTimeRangeStandard");
    x->SetBinLabel(10, "NoCollInRofStandard");
  }

  template <typename CheckTrack>
  bool isTrackSelected(CheckTrack const& track)
  {
    if (std::abs(track.eta()) > etaRange) {
      return false;
    }
    /*if (isApplyExtraPhiCut && ((track.phi() > 3.07666 && track.phi() < 3.12661) || track.phi() <= 0.03 || track.phi() >= 6.253)) {
      return false;
    }*/
    return true;
  }

  template <typename TCollision>
  bool isEventSelected(TCollision collision)
  {
    registry.fill(HIST("EventHist"), 1);

    if (!collision.sel8()) {
      return false;
    }
    registry.fill(HIST("EventHist"), 2);

    if (std::fabs(collision.posZ()) > cfgEvSelVtxZ) {
      return false;
    }
    registry.fill(HIST("EventHist"), 3);

    if (isApplySameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    registry.fill(HIST("EventHist"), 4);

    if (isApplyGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    registry.fill(HIST("EventHist"), 5);

    if (isApplyVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    registry.fill(HIST("EventHist"), 6);

    auto occupancy = collision.trackOccupancyInTimeRange();
    if (occupancy > cfgEvSelsMaxOccupancy) {
      return false;
    }
    registry.fill(HIST("EventHist"), 7);

    if (!collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }
    registry.fill(HIST("EventHist"), 8);

    if (isApplyNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    registry.fill(HIST("EventHist"), 9);

    if (isApplyNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    registry.fill(HIST("EventHist"), 10);

    return true;
  }

  void processZDCautoTrig(aod::Zdc const& zdc)
  {
    // auto-triggered events for ZDC

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
      registry.get<TH1>(HIST("debunchcheck"))->Fill(zna - znc, zna + znc);
    }

    zdcTableLI(tdcZNA, zna, pmcZNA, pmqZNA[0], pmqZNA[1], pmqZNA[2], pmqZNA[3],
               tdcZNC, znc, pmcZNC, pmqZNC[0], pmqZNC[1], pmqZNC[2], pmqZNC[3],
               tdcZPA, zpa, tdcZPC, zpc, tdcZEM1, zem1, tdcZEM2, zem2,
               -1, -1, -1,
               -1, -1.,
               -1, -1, -1,
               -1);
  }
  /// name, description, function pointer, default value
  /// note that it has to be declared after the function, so that the pointer is known
  PROCESS_SWITCH(ZdcTaskLightIons, processZDCautoTrig, "Processing ZDC 4 auto-triggered events", true);

  void processALICEcoll(ColEvSels const& cols, BCsRun3 const& /*bcs*/, aod::Tracks const& tracks, aod::Zdcs const& /*zdcs*/)
  {
    // collision-based event selection
    for (auto const& collision : cols) {

      const auto& foundBC = collision.foundBC_as<BCsRun3>();
      uint8_t evSelection = isEventSelected(collision);
      auto zv = collision.posZ();

      auto nTracks = 0;
      for (const auto& track : tracks) {
        if (!isTrackSelected(track)) {
          continue;
        }
        nTracks++;
      }

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
          registry.get<TH1>(HIST("debunchcheck"))->Fill(zna - znc, zna + znc);
        }

        zdcTableLI(tdcZNA, zna, pmcZNA, pmqZNA[0], pmqZNA[1], pmqZNA[2], pmqZNA[3],
                   tdcZNC, znc, pmcZNC, pmqZNC[0], pmqZNC[1], pmqZNC[2], pmqZNC[3],
                   tdcZPA, zpa, tdcZPC, zpc, tdcZEM1, zem1, tdcZEM2, zem2,
                   multFT0A, multFT0C, multV0A,
                   nTracks, zv,
                   centralityFT0C, centralityFT0A, centralityFT0M,
                   evSelection);
      }
    }
  }
  /// name, description, function pointer, default value
  /// note that it has to be declared after the function, so that the pointer is known
  PROCESS_SWITCH(ZdcTaskLightIons, processALICEcoll, "Processing ZDC for ALICE collisions", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ZdcTaskLightIons>(cfgc)};
}
