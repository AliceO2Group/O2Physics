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

/// \file zdcExtraTableProducer.cxx
/// \brief Task creating table with ZDC PMTs energies and calculated centroid (Q-vector) to be used for spectator plane measurement
/// \author Chiara Oppedisano <chiara.oppedisano@cern.ch>, INFN Torino
/// \author Uliana Dmitrieva <uliana.dmitrieva@cern.ch>, INFN Torino

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/ZDCExtra.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TRandom3.h>

#include <cstdint>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::evsel;

using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;

struct ZdcExtraTableProducer {

  Produces<aod::ZdcExtras> zdcextras;

  // Configurable parameters
  //
  Configurable<int> nBins{"nBins", 400, "n bins"};
  Configurable<float> maxZN{"maxZN", 399.5, "Max ZN signal"};
  Configurable<bool> tdcCut{"tdcCut", false, "Flag for TDC cut"};
  Configurable<float> tdcZNmincut{"tdcZNmincut", -2.5, "Min ZN TDC cut"};
  Configurable<float> tdcZNmaxcut{"tdcZNmaxcut", 2.5, "Max ZN TDC cut"};
  Configurable<bool> cfgUsePMC{"cfgUsePMC", true, "Use common PM (true) or sum of PMs (false) "};
  // Event selections
  Configurable<bool> cfgEvSelSel8{"cfgEvSelSel8", true, "Event selection: sel8"};
  Configurable<float> cfgEvSelVtxZ{"cfgEvSelVtxZ", 10, "Event selection: zVtx"};
  Configurable<bool> cfgEvSelsDoOccupancySel{"cfgEvSelsDoOccupancySel", true, "Event selection: do occupancy selection"};
  Configurable<float> cfgEvSelsMaxOccupancy{"cfgEvSelsMaxOccupancy", 10000, "Event selection: set max occupancy"};
  Configurable<bool> cfgEvSelsNoSameBunchPileupCut{"cfgEvSelsNoSameBunchPileupCut", true, "Event selection: no same bunch pileup cut"};
  Configurable<bool> cfgEvSelsIsGoodZvtxFT0vsPV{"cfgEvSelsIsGoodZvtxFT0vsPV", true, "Event selection: is good ZVTX FT0 vs PV"};
  Configurable<bool> cfgEvSelsNoCollInTimeRangeStandard{"cfgEvSelsNoCollInTimeRangeStandard", true, "Event selection: no collision in time range standard"};
  Configurable<bool> cfgEvSelsIsVertexITSTPC{"cfgEvSelsIsVertexITSTPC", true, "Event selection: is vertex ITSTPC"};
  Configurable<bool> cfgEvSelsIsGoodITSLayersAll{"cfgEvSelsIsGoodITSLayersAll", true, "Event selection: is good ITS layers all"};
  // Calibration settings
  Configurable<float> cfgCalibrationDownscaling{"cfgCalibrationDownscaling", 1.f, "Percentage of events to be saved to derived table"};

  HistogramRegistry registry{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  enum SelectionCriteria {
    evSel_zvtx,
    evSel_sel8,
    evSel_occupancy,
    evSel_kNoSameBunchPileup,
    evSel_kIsGoodZvtxFT0vsPV,
    evSel_kNoCollInTimeRangeStandard,
    evSel_kIsVertexITSTPC,
    evSel_kIsGoodITSLayersAll,
    evSel_allEvents,
    nEventSelections
  };

  void init(InitContext const&)
  {
    registry.add("ZNApmc", "ZNApmc; ZNA PMC; Entries", {HistType::kTH1F, {{nBins, -0.5, maxZN}}});
    registry.add("ZNCpmc", "ZNCpmc; ZNC PMC; Entries", {HistType::kTH1F, {{nBins, -0.5, maxZN}}});
    registry.add("ZNApm1", "ZNApm1; ZNA PM1; Entries", {HistType::kTH1F, {{nBins, -0.5, maxZN}}});
    registry.add("ZNApm2", "ZNApm2; ZNA PM2; Entries", {HistType::kTH1F, {{nBins, -0.5, maxZN}}});
    registry.add("ZNApm3", "ZNApm3; ZNA PM3; Entries", {HistType::kTH1F, {{nBins, -0.5, maxZN}}});
    registry.add("ZNApm4", "ZNApm4; ZNA PM4; Entries", {HistType::kTH1F, {{nBins, -0.5, maxZN}}});
    registry.add("ZNCpm1", "ZNCpm1; ZNC PM1; Entries", {HistType::kTH1F, {{nBins, -0.5, maxZN}}});
    registry.add("ZNCpm2", "ZNCpm2; ZNC PM2; Entries", {HistType::kTH1F, {{nBins, -0.5, maxZN}}});
    registry.add("ZNCpm3", "ZNCpm3; ZNC PM3; Entries", {HistType::kTH1F, {{nBins, -0.5, maxZN}}});
    registry.add("ZNCpm4", "ZNCpm4; ZNC PM4; Entries", {HistType::kTH1F, {{nBins, -0.5, maxZN}}});
    registry.add("ZNAsumq", "ZNAsumq; ZNA uncalib. sum PMQ; Entries", {HistType::kTH1F, {{nBins, -0.5, maxZN}}});
    registry.add("ZNCsumq", "ZNCsumq; ZNC uncalib. sum PMQ; Entries", {HistType::kTH1F, {{nBins, -0.5, maxZN}}});

    registry.add("ZNACentroid", "ZNA Centroid; X; Y", {HistType::kTH2F, {{50, -1.5, 1.5}, {50, -1.5, 1.5}}});
    registry.add("ZNCCentroid", "ZNC Centroid; X; Y", {HistType::kTH2F, {{50, -1.5, 1.5}, {50, -1.5, 1.5}}});

    registry.add("hEventCount", "Number of Event; Cut; #Events Passed Cut", {HistType::kTH1D, {{nEventSelections, 0, nEventSelections}}});
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_allEvents + 1, "All events");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_zvtx + 1, "vtxZ");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_sel8 + 1, "Sel8");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_occupancy + 1, "kOccupancy");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kNoSameBunchPileup + 1, "kNoSameBunchPileup");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kIsGoodZvtxFT0vsPV + 1, "kIsGoodZvtxFT0vsPV");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kNoCollInTimeRangeStandard + 1, "kNoCollInTimeRangeStandard");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kIsVertexITSTPC + 1, "kIsVertexITSTPC");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kIsGoodITSLayersAll + 1, "kIsGoodITSLayersAll");
  }

  template <typename TCollision>
  uint8_t eventSelected(TCollision collision)
  {
    uint8_t selectionBits = 0;
    bool selected;

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

    selected = collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC);
    if (selected) {
      selectionBits |= (uint8_t)(0x1u << evSel_kIsVertexITSTPC);
      registry.fill(HIST("hEventCount"), evSel_kIsVertexITSTPC);
    }

    selected = collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll);
    if (selected) {
      selectionBits |= (uint8_t)(0x1u << evSel_kIsGoodITSLayersAll);
      registry.fill(HIST("hEventCount"), evSel_kIsGoodITSLayersAll);
    }

    return selectionBits;
  }

  void process(ColEvSels const& cols, BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcs*/)
  {
    // collision-based event selection
    int nTowers = 4; // number of ZDC towers

    for (auto const& collision : cols) {
      const auto& foundBC = collision.foundBC_as<BCsRun3>();
      if (foundBC.has_zdc()) {
        const auto& zdc = foundBC.zdc();

        uint8_t evSelection = eventSelected(collision);

        float centrality = collision.centFT0C();

        // To assure that ZN have a genuine signal (tagged by the relative TDC)
        // we can check that the amplitude is >0 or that ADC is NOT very negative (-inf)

        double pmcZNC = zdc.energyCommonZNC();
        double pmcZNA = zdc.energyCommonZNA();
        bool isZNChit = false, isZNAhit = false;
        //
        double tdcZNC = zdc.timeZNC();
        double tdcZNA = zdc.timeZNA();
        // OR we can select a narrow window in both ZN TDCs using the configurable parameters
        if (tdcCut) { // a narrow TDC window is set
          if ((tdcZNC >= tdcZNmincut) && (tdcZNC <= tdcZNmaxcut)) {
            isZNChit = true;
          }
          if ((tdcZNA >= tdcZNmincut) && (tdcZNA <= tdcZNmaxcut)) {
            isZNAhit = true;
          }
        } else { // if no window on TDC is set
          if (pmcZNC > -1.) {
            isZNChit = true;
          }
          if (pmcZNA > -1.) {
            isZNAhit = true;
          }
        }
        //
        double sumZNC = 0;
        double sumZNA = 0;
        double pmqZNC[4] = {};
        double pmqZNA[4] = {};
        //
        if (isZNChit) {
          for (int it = 0; it < nTowers; it++) {
            pmqZNC[it] = (zdc.energySectorZNC())[it];
            sumZNC += pmqZNC[it];
          }
          registry.get<TH1>(HIST("ZNCpmc"))->Fill(pmcZNC);
          registry.get<TH1>(HIST("ZNCpm1"))->Fill(pmqZNC[0]);
          registry.get<TH1>(HIST("ZNCpm2"))->Fill(pmqZNC[1]);
          registry.get<TH1>(HIST("ZNCpm3"))->Fill(pmqZNC[2]);
          registry.get<TH1>(HIST("ZNCpm4"))->Fill(pmqZNC[3]);
          registry.get<TH1>(HIST("ZNCsumq"))->Fill(sumZNC);
        }
        if (isZNAhit) {
          for (int it = 0; it < nTowers; it++) {
            pmqZNA[it] = (zdc.energySectorZNA())[it];
            sumZNA += pmqZNA[it];
          }
          //
          registry.get<TH1>(HIST("ZNApmc"))->Fill(pmcZNA);
          registry.get<TH1>(HIST("ZNApm1"))->Fill(pmqZNA[0]);
          registry.get<TH1>(HIST("ZNApm2"))->Fill(pmqZNA[1]);
          registry.get<TH1>(HIST("ZNApm3"))->Fill(pmqZNA[2]);
          registry.get<TH1>(HIST("ZNApm4"))->Fill(pmqZNA[3]);
          registry.get<TH1>(HIST("ZNAsumq"))->Fill(sumZNA);
        }

        // Q-vectors (centroid) calculation
        // kBeamEne --  LHC Run 3 Pb-Pb collision energy (5.36 TeV per nucleon pair)
        constexpr float kBeamEne = 5.36 * 0.5;

        // Provide coordinates of centroid over ZN (side C) front face
        constexpr float X[4] = {-1.75, 1.75, -1.75, 1.75};
        constexpr float Y[4] = {-1.75, -1.75, 1.75, 1.75};
        constexpr float kAlpha = 0.395; // saturation correction

        float numXZNC = 0., numYZNC = 0., denZNC = 0.;
        float numXZNA = 0., numYZNA = 0., denZNA = 0.;

        // Calculate weighted sums of the x and y coordinates
        constexpr int kNTowers = 4; // number of ZDC towers
        for (int i = 0; i < kNTowers; i++) {
          if (pmqZNC[i] > 0.) {
            float wZNC = std::pow(pmqZNC[i], kAlpha);
            numXZNC -= X[i] * wZNC; // numerator x (minus sign due to opposite orientation of ZNC)
            numYZNC += Y[i] * wZNC; // numerator y
            denZNC += wZNC;         // denominator
          }
          if (pmqZNA[i] > 0.) {
            float wZNA = std::pow(pmqZNA[i], kAlpha);
            numXZNA += X[i] * wZNA; // numerator x
            numYZNA += Y[i] * wZNA; // numerator y
            denZNA += wZNA;         // denominator
          }
        }
        // Calculate centroid coordinates (in cm) with correction factor c depending on the number of spectator nucleons (nSpec)

        float zncCommon = 0;
        float znaCommon = 0;

        // Use sum of PMTs (cfgUsePMC == false) when common PMT is saturated
        if (cfgUsePMC) {
          zncCommon = pmcZNC;
          znaCommon = pmcZNA;
        } else {
          zncCommon = sumZNC;
          znaCommon = sumZNA;
        }

        float centroidZNC[2], centroidZNA[2];

        if (denZNC != 0.) {
          float nSpecnC = zncCommon / kBeamEne;
          float cZNC = 1.89358 - 0.71262 / (nSpecnC + 0.71789);
          centroidZNC[0] = cZNC * numXZNC / denZNC;
          centroidZNC[1] = cZNC * numYZNC / denZNC;
        } else {
          centroidZNC[0] = 999.;
          centroidZNC[1] = 999.;
        }
        //
        if (denZNA != 0.) {
          float nSpecnA = znaCommon / kBeamEne;
          float cZNA = 1.89358 - 0.71262 / (nSpecnA + 0.71789);
          centroidZNA[0] = cZNA * numXZNA / denZNA;
          centroidZNA[1] = cZNA * numYZNA / denZNA;
        } else {
          centroidZNA[0] = 999.;
          centroidZNA[1] = 999.;
        }
        registry.get<TH2>(HIST("ZNCCentroid"))->Fill(centroidZNC[0], centroidZNC[1]);
        registry.get<TH2>(HIST("ZNACentroid"))->Fill(centroidZNA[0], centroidZNA[1]);

        auto vz = collision.posZ();
        auto vx = collision.posX();
        auto vy = collision.posY();

        if ((isZNAhit || isZNChit) && (gRandom->Uniform() < cfgCalibrationDownscaling)) {
          zdcextras(pmcZNA, pmqZNA[0], pmqZNA[1], pmqZNA[2], pmqZNA[3], tdcZNA, centroidZNA[0], centroidZNA[1], pmcZNC, pmqZNC[0], pmqZNC[1], pmqZNC[2], pmqZNC[3], tdcZNC, centroidZNC[0], centroidZNC[1], centrality, vx, vy, vz, foundBC.timestamp(), foundBC.runNumber(), evSelection);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ZdcExtraTableProducer>(cfgc)};
}
