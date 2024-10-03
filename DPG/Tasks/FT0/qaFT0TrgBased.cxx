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

/// \file   qaFT0TrgBased.cxx
/// \author Uliana Dmitrieva uliana.dmitrieva@cern.ch
/// \brief  FT0 QA

#include <bitset>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "DataFormatsFT0/Digit.h"

#include "CommonConstants/LHCConstants.h"
#include "DataFormatsFIT/Triggers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"
#include "Framework/InitContext.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;

using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::BcSels, aod::Timestamps,
                                       aod::Run3MatchedToBCSparse>;

struct qaFT0TrgBased {

  // Configurable<bool> isMC{"isMC", 0, "0 - data, 1 - MC"};
  Configurable<int> selection{"selection", 0, "trigger: 0 - no sel, 8 - sel8"};
  Configurable<bool> isLowFlux{"isLowFlux", 1,
                               "1 - low flux (pp, pPb), 0 - high flux (PbPb)"};

  // Event selection conditions
  // Caution required if order is changed.
  enum eConditions {
    kAll,
    kHasFT0,
    kOrAFT0,
    kOrCFT0,
    kOrACFT0,
    k0TVX,
    kCent,
    kSemiCent,
    kMinBias,
    kNConditions
  };

  // Event selection condition names
  // Caution required if names are changed. See FillConditionhistos macro below
  // (+ post processing scripts)
  inline static std::map<eConditions, const char*> mapConditionNames = {
    {kAll, "All"},
    {kHasFT0, "HasFT0"},
    {kOrAFT0, "FT0OrA"},
    {kOrCFT0, "FT0OrC"},
    {kOrACFT0, "FT0OrAC"},
    {k0TVX, "0TVX"},
    {kCent, "Cent"},
    {kSemiCent, "SemiCent"},
    {kMinBias, "MinBias"}};

  // Observables
  enum eObservables {
    kT0A,
    kT0C,
    kT0AC,
    kT0A_V0A,
    kT0res,
    kAmpT0A,
    kAmpT0C,
    kMultT0A,
    kMultT0C,
    kMultT0AC,
    kMultT0A_V0A,
    kMultT0A_ZNA,
    kMultT0C_ZNC,
    kMultZNA_ZNC,
    kT0vertex,
    kPV,
    kT0V_PV_diff,
    kT0V_PV,
    kContrib,
    kNObservables
  };

  // Observable names
  inline static std::map<eObservables, const char*> mapObservableNames = {
    {kT0A, "hT0A"},
    {kT0C, "hT0C"},
    {kT0AC, "hT0AC"},
    {kT0A_V0A, "hT0A_V0A"},
    {kT0res, "hT0res"},
    {kAmpT0A, "hAmpT0A"},
    {kAmpT0C, "hAmpT0C"},
    {kMultT0A, "hMultT0A"},
    {kMultT0C, "hMultT0C"},
    {kMultT0AC, "hMultT0AC"},
    {kMultT0A_V0A, "hMultT0A_V0A"},
    {kMultT0A_ZNA, "hMultT0A_ZNA"},
    {kMultT0C_ZNC, "hMultT0C_ZNC"},
    {kMultZNA_ZNC, "hMultZNA_ZNC"},
    {kT0vertex, "hT0vertex"},
    {kPV, "hPV"},
    {kT0V_PV_diff, "hT0V_PV_diff"},
    {kT0V_PV, "hT0V_PV"},
    {kContrib, "hContrib"}};

  HistogramRegistry histos{
    "Histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {

    if (kNConditions != mapConditionNames.size()) {
      LOG(fatal)
        << "Number of conditions does not match number of condition names";
    }
    if (kNObservables != mapObservableNames.size()) {
      LOG(fatal)
        << "Number of observables does not match number of observable names";
    }

    const AxisSpec axisTime{500, -5., 5., "collision time (ns)"};
    const AxisSpec axisColTimeRes{isLowFlux ? 500 : 2000, -0.5, 0.5,
                                  "(T0A - T0C)/2 (ns)"};
    const AxisSpec axisVertex{300, -30., 30.};

    const AxisSpec axisAmpT0A{isLowFlux ? 500 : 2000, 0.,
                              isLowFlux ? 500. : 2000.,
                              "T0A amplitude (# ADC channels)"};

    const AxisSpec axisAmpT0C{isLowFlux ? 500 : 2000, 0.,
                              isLowFlux ? 500. : 2000.,
                              "T0C amplitude (# ADC channels)"};

    const AxisSpec axisMultT0A{5000, 0., isLowFlux ? 10000. : 200000.,
                               "T0A multiplicity (# ADC channels)"};
    const AxisSpec axisMultT0C{1000, 0., isLowFlux ? 2000. : 80000.,
                               "T0C multiplicity (# ADC channels)"};
    const AxisSpec axisMultT0AC{12000, 0., isLowFlux ? 12000. : 240000.,
                                "T0AC multiplicity (# ADC channels)"};

    const AxisSpec axisMultV0A{2000, 0., isLowFlux ? 40000. : 200000.,
                               "V0A multiplicity"};

    const AxisSpec axisMultZNA{isLowFlux ? 500 : 1500, 0.,
                               isLowFlux ? 50. : 300., "E_{ZNA} (TeV)"};

    const AxisSpec axisMultZNC{isLowFlux ? 500 : 1500, 0.,
                               isLowFlux ? 50. : 300., "E_{ZNC} (TeV)"};

    const AxisSpec axisNcontrib{isLowFlux ? 150 : 4500, 0.,
                                isLowFlux ? 150. : 4500., "# contributors"};

    const AxisSpec axisCounter{
      2, 0.5, 1.5, ""}; // should be deleted as soon as BC task rewritten

    const AxisSpec axisEvSelStats{kNConditions, 0, kNConditions};

    // Lambda that creates one histogram for 'observable' per condition,
    // excluding kAll condition
    bool isBCs = 0;

    auto const makeConditionHistos = [&](bool isBCs, eObservables observable,
                                         HistType histType,
                                         const std::vector<AxisSpec>& axis) {
      for (int iCondition = 1; iCondition < kNConditions; iCondition++) {

        const std::string prefix = isBCs ? "BCs/" : "Colls/";

        const std::string condition =
          mapConditionNames[static_cast<eConditions>(iCondition)];
        const std::string histoName = prefix + mapObservableNames[observable] +
                                      std::string("/") + condition;
        histos.add(histoName, histoName.c_str(), histType, axis);
      }
    };

    // Collisions task
    if (doprocessCollisions) {

      isBCs = 0;

      makeConditionHistos(isBCs, kT0A, kTH1F, {axisTime});
      makeConditionHistos(isBCs, kT0C, kTH1F, {axisTime});
      makeConditionHistos(isBCs, kT0AC, kTH1F, {axisTime});
      makeConditionHistos(isBCs, kT0A_V0A, kTH2F, {axisTime, axisTime});

      makeConditionHistos(isBCs, kT0res, kTH1F, {axisColTimeRes});

      makeConditionHistos(isBCs, kAmpT0A, kTH1F, {axisAmpT0A});
      makeConditionHistos(isBCs, kAmpT0C, kTH1F, {axisAmpT0C});

      makeConditionHistos(isBCs, kMultT0A, kTH1F, {axisMultT0A});
      makeConditionHistos(isBCs, kMultT0C, kTH1F, {axisMultT0C});
      makeConditionHistos(isBCs, kMultT0AC, kTH1F, {axisMultT0AC});

      makeConditionHistos(isBCs, kMultT0A_V0A, kTH2F,
                          {axisMultT0A, axisMultV0A});

      makeConditionHistos(isBCs, kMultT0A_ZNA, kTH2F,
                          {axisMultT0A, axisMultZNA});
      makeConditionHistos(isBCs, kMultT0C_ZNC, kTH2F,
                          {axisMultT0C, axisMultZNC});

      makeConditionHistos(isBCs, kT0vertex, kTH1F, {axisVertex});

      histos.add("Colls/hPV/All", "PV;primary vertex (cm);counts", kTH1F,
                 {axisVertex});
      makeConditionHistos(isBCs, kPV, kTH1F, {axisVertex});

      makeConditionHistos(isBCs, kT0V_PV_diff, kTH1F, {axisVertex});
      makeConditionHistos(isBCs, kT0V_PV, kTH2F, {axisVertex, axisVertex});

      histos.add("Colls/hContrib/All", ";(# contrubutors);counts", kTH1F,
                 {axisNcontrib});
      makeConditionHistos(isBCs, kContrib, kTH1F, {axisNcontrib});

      // Histogram for storing event selection statistics
      auto hc = histos.add<TH1>("Colls/EventSelectionStats",
                                "Event selection statistics", kTH1F,
                                {axisEvSelStats});

      for (int iCondition = 0; iCondition < kNConditions; iCondition++) {
        hc->GetXaxis()->SetBinLabel(
          iCondition + 1,
          mapConditionNames[static_cast<eConditions>(iCondition)]);
      }
    }

    //_______________________________________________________________________________________________________________________
    // BC task

    if (doprocessBCs) {
      isBCs = 1;

      makeConditionHistos(isBCs, kT0A, kTH1F, {axisTime});
      makeConditionHistos(isBCs, kT0C, kTH1F, {axisTime});
      makeConditionHistos(isBCs, kT0AC, kTH1F, {axisTime});
      makeConditionHistos(isBCs, kT0A_V0A, kTH2F, {axisTime, axisTime});

      makeConditionHistos(isBCs, kAmpT0A, kTH1F, {axisAmpT0A});
      makeConditionHistos(isBCs, kAmpT0C, kTH1F, {axisAmpT0C});

      makeConditionHistos(isBCs, kMultT0A, kTH1F, {axisMultT0A});
      makeConditionHistos(isBCs, kMultT0C, kTH1F, {axisMultT0C});
      makeConditionHistos(isBCs, kMultT0AC, kTH1F, {axisMultT0AC});

      makeConditionHistos(isBCs, kMultT0A_ZNA, kTH2F,
                          {axisMultT0A, axisMultZNA});
      makeConditionHistos(isBCs, kMultT0C_ZNC, kTH2F,
                          {axisMultT0C, axisMultZNC});

      makeConditionHistos(isBCs, kT0vertex, kTH1F, {axisVertex});

      histos.add("BCs/hMultZNA_ZNC/All", ";(# contrubutors);counts", kTH2F,
                 {axisMultZNA, axisMultZNC}); // just to check ZDC
      histos.add("BCs/hMultZNA_ZNC_EM/All", ";(# contrubutors);counts", kTH2F,
                 {axisMultZNA, axisMultZNC}); // just to check ZDC

      // Histogram for storing event selection statistics
      auto hbc = histos.add<TH1>("BCs/EventSelectionStats",
                                 "Event selection statistics", kTH1F,
                                 {axisEvSelStats});

      for (int iCondition = 0; iCondition < kNConditions; iCondition++) {
        hbc->GetXaxis()->SetBinLabel(
          iCondition + 1,
          mapConditionNames[static_cast<eConditions>(iCondition)]);
      }
    }
  }

  void
    processCollisions(soa::Join<aod::Collisions, aod::EvSels, aod::Mults,
                                aod::FT0sCorrected>::iterator const& collision,
                      aod::FT0s const&, aod::FV0As const&,
                      aod::Zdcs const&)
  {

    if (selection == 8 && !collision.sel8()) {
      return;
    }

    histos.fill(HIST("Colls/EventSelectionStats"), kAll);

    const bool hasFT0 = collision.has_foundFT0();
    const bool hasZDC = (collision.foundZDCId() >= 0);
    const bool hasFV0 = collision.has_foundFV0();

    bool OrAFT0 = false;
    bool OrCFT0 = false;
    bool OrACFT0 = false;
    bool TVX = false;
    bool Cent = false;
    bool SemiCent = false;
    bool MinBias = false;

    float multFT0A = 0.f;
    float multFT0C = 0.f;
    float multFT0M = 0.f;
    float multFV0 = 0.f;

    // is FT0
    if (hasFT0) {

      histos.fill(HIST("Colls/EventSelectionStats"), kHasFT0);

      auto ft0 = collision.foundFT0();

      std::bitset<8> triggers = ft0.triggerMask();
      OrAFT0 = triggers[o2::ft0::Triggers::bitA];
      OrCFT0 = triggers[o2::ft0::Triggers::bitC];
      OrACFT0 = (OrAFT0 && OrCFT0);
      TVX = triggers[o2::ft0::Triggers::bitVertex];
      Cent = triggers[o2::ft0::Triggers::bitCen];
      SemiCent = triggers[o2::ft0::Triggers::bitSCen];
      MinBias = (TVX && (Cent || SemiCent));

      if (OrAFT0) {
        histos.fill(HIST("Colls/EventSelectionStats"), kOrAFT0);
      }
      if (OrCFT0) {
        histos.fill(HIST("Colls/EventSelectionStats"), kOrCFT0);
      }
      if (OrACFT0) {
        histos.fill(HIST("Colls/EventSelectionStats"), kOrACFT0);
      }
      if (TVX) {
        histos.fill(HIST("Colls/EventSelectionStats"), k0TVX);
      }
      if (Cent) {
        histos.fill(HIST("Colls/EventSelectionStats"), kCent);
      }
      if (SemiCent) {
        histos.fill(HIST("Colls/EventSelectionStats"), kSemiCent);
      }
      if (MinBias) {
        histos.fill(HIST("Colls/EventSelectionStats"), kMinBias);
      }

      // Macro for filling histograms for 'observable' based on conditions for
      // CollisionsTask
#define FillConditionHistos(observable, ...)                           \
  if (hasFT0) {                                                        \
    histos.fill(HIST("Colls/" observable "/HasFT0"), __VA_ARGS__);     \
    if (OrAFT0) {                                                      \
      histos.fill(HIST("Colls/" observable "/FT0OrA"), __VA_ARGS__);   \
    }                                                                  \
    if (OrCFT0) {                                                      \
      histos.fill(HIST("Colls/" observable "/FT0OrC"), __VA_ARGS__);   \
    }                                                                  \
    if (OrACFT0) {                                                     \
      histos.fill(HIST("Colls/" observable "/FT0OrAC"), __VA_ARGS__);  \
    }                                                                  \
    if (TVX) {                                                         \
      histos.fill(HIST("Colls/" observable "/0TVX"), __VA_ARGS__);     \
    }                                                                  \
    if (Cent) {                                                        \
      histos.fill(HIST("Colls/" observable "/Cent"), __VA_ARGS__);     \
    }                                                                  \
    if (SemiCent) {                                                    \
      histos.fill(HIST("Colls/" observable "/SemiCent"), __VA_ARGS__); \
    }                                                                  \
    if (MinBias) {                                                     \
      histos.fill(HIST("Colls/" observable "/MinBias"), __VA_ARGS__);  \
    }                                                                  \
  }

      if (ft0.isValidTimeA()) {
        FillConditionHistos("hT0A", ft0.timeA());
      }
      if (ft0.isValidTimeC()) {
        FillConditionHistos("hT0C", ft0.timeC());
      }

      if (ft0.isValidTimeA() && ft0.isValidTimeC()) {
        FillConditionHistos("hT0AC", ft0.collTime());
        FillConditionHistos("hT0vertex", ft0.posZ());
        FillConditionHistos("hT0V_PV", ft0.posZ(), collision.posZ());
        FillConditionHistos("hT0V_PV_diff", ft0.posZ() - collision.posZ());
        FillConditionHistos("hPV", ft0.posZ());
      }

      if (collision.t0CCorrectedValid() && collision.t0ACorrectedValid()) {
        FillConditionHistos("hT0res", collision.t0resolution());
      }

      for (std::size_t i_a = 0; i_a < ft0.amplitudeA().size(); i_a++) {

        float amplitudeA = ft0.amplitudeA()[i_a];
        //   uint8_t channel = ft0.channelA()[i_a];
        FillConditionHistos("hAmpT0A", amplitudeA);
      }

      for (std::size_t i_c = 0; i_c < ft0.amplitudeC().size(); i_c++) {
        float amplitudeC = ft0.amplitudeC()[i_c];
        //  uint8_t channel = ft0.channelC()[i_c];
        FillConditionHistos("hAmpT0C", amplitudeC);
      }

      multFT0A = collision.multFT0A();
      multFT0C = collision.multFT0C();
      multFT0M = multFT0A + multFT0C;

      FillConditionHistos("hMultT0A", multFT0A);
      FillConditionHistos("hMultT0C", multFT0C);
      FillConditionHistos("hMultT0AC", multFT0M);

      // Multiplicity correlations with ZDC

      float multZNA = 0;
      float multZNC = 0;

      // is ZDC
      if (hasZDC) {

        multZNA = collision.multZNA();
        multZNC = collision.multZNC();

        FillConditionHistos("hMultT0A_ZNA", multFT0A, multZNA);
        FillConditionHistos("hMultT0C_ZNC", multFT0C, multZNC);
      } // end of if (collision.foundZDCId() >= 0)

      // is FV0
      if (hasFV0) {
        auto fv0 = collision.foundFV0();

        FillConditionHistos("hT0A_V0A", fv0.time(), ft0.timeA());

        multFV0 = collision.multFV0A();
        FillConditionHistos("hMultT0A_V0A", multFT0A, multFV0);

      } // end of if (collision.has_foundFV0())

    } // end of if (collision.has_foundFT0())

    // number of contributers used for vertex calculation
    int nContrib = collision.numContrib();
    histos.fill(HIST("Colls/hContrib/All"), nContrib);
    FillConditionHistos("hContrib", nContrib);

    histos.fill(HIST("Colls/hPV/All"), collision.posZ());

#undef FillConditionHistos

  } // end of processCollsions()

  PROCESS_SWITCH(qaFT0TrgBased, processCollisions, "per-collision analysis",
                 true);

  void processBCs(BCsWithRun3Matchings::iterator const& bc, aod::FV0As const&,
                  aod::FT0s const&, aod::Zdcs const&)
  {

    histos.fill(HIST("BCs/EventSelectionStats"), kAll);

    bool OrAFT0 = false;
    bool OrCFT0 = false;
    bool OrACFT0 = false;
    bool TVX = false;
    bool Cent = false;
    bool SemiCent = false;
    bool MinBias = false;

    float multFT0A = 0.f;
    float multFT0C = 0.f;
    float multFT0M = 0.f;
    float multZNA = 0;
    float multZNC = 0;

    const bool hasFT0 = bc.has_ft0();
    const bool hasZDC = bc.has_zdc();
    const bool hasFV0 = bc.has_fv0a();

    // isFT0
    if (hasFT0) {

      auto ft0 = bc.ft0();

      histos.fill(HIST("BCs/EventSelectionStats"), kHasFT0);

      std::bitset<8> triggers = ft0.triggerMask();
      OrAFT0 = triggers[o2::ft0::Triggers::bitA];
      OrCFT0 = triggers[o2::ft0::Triggers::bitC];
      OrACFT0 = (OrAFT0 && OrCFT0);
      TVX = triggers[o2::ft0::Triggers::bitVertex];
      Cent = triggers[o2::ft0::Triggers::bitCen];
      SemiCent = triggers[o2::ft0::Triggers::bitSCen];
      MinBias = (TVX && (Cent || SemiCent));

      if (OrAFT0) {
        histos.fill(HIST("BCs/EventSelectionStats"), kOrAFT0);
      }
      if (OrCFT0) {
        histos.fill(HIST("BCs/EventSelectionStats"), kOrCFT0);
      }
      if (OrACFT0) {
        histos.fill(HIST("BCs/EventSelectionStats"), kOrACFT0);
      }
      if (TVX) {
        histos.fill(HIST("BCs/EventSelectionStats"), k0TVX);
      }
      if (Cent) {
        histos.fill(HIST("BCs/EventSelectionStats"), kCent);
      }
      if (SemiCent) {
        histos.fill(HIST("BCs/EventSelectionStats"), kSemiCent);
      }
      if (MinBias) {
        histos.fill(HIST("BCs/EventSelectionStats"), kMinBias);
      }

      // Macro for filling histograms for 'observable' based on conditions for
      // BCsTask
#define FillConditionHistosBCs(observable, ...)                      \
  if (hasFT0) {                                                      \
    histos.fill(HIST("BCs/" observable "/HasFT0"), __VA_ARGS__);     \
    if (OrAFT0) {                                                    \
      histos.fill(HIST("BCs/" observable "/FT0OrA"), __VA_ARGS__);   \
    }                                                                \
    if (OrCFT0) {                                                    \
      histos.fill(HIST("BCs/" observable "/FT0OrC"), __VA_ARGS__);   \
    }                                                                \
    if (OrACFT0) {                                                   \
      histos.fill(HIST("BCs/" observable "/FT0OrAC"), __VA_ARGS__);  \
    }                                                                \
    if (TVX) {                                                       \
      histos.fill(HIST("BCs/" observable "/0TVX"), __VA_ARGS__);     \
    }                                                                \
    if (Cent) {                                                      \
      histos.fill(HIST("BCs/" observable "/Cent"), __VA_ARGS__);     \
    }                                                                \
    if (SemiCent) {                                                  \
      histos.fill(HIST("BCs/" observable "/SemiCent"), __VA_ARGS__); \
    }                                                                \
    if (MinBias) {                                                   \
      histos.fill(HIST("BCs/" observable "/MinBias"), __VA_ARGS__);  \
    }                                                                \
  }

      if (ft0.isValidTimeA()) {
        FillConditionHistosBCs("hT0A", ft0.timeA());
      }
      if (ft0.isValidTimeC()) {
        FillConditionHistosBCs("hT0C", ft0.timeC());
      }

      if (ft0.isValidTimeA() && ft0.isValidTimeC()) {
        FillConditionHistosBCs("hT0AC", ft0.collTime());
        FillConditionHistosBCs("hT0vertex", ft0.posZ());
      }

      for (std::size_t i_a = 0; i_a < ft0.amplitudeA().size(); i_a++) {

        float amplitudeA = ft0.amplitudeA()[i_a];
        //   uint8_t channel = ft0.channelA()[i_a];
        FillConditionHistosBCs("hAmpT0A", amplitudeA);
      }

      for (std::size_t i_c = 0; i_c < ft0.amplitudeC().size(); i_c++) {
        float amplitudeC = ft0.amplitudeC()[i_c];
        //  uint8_t channel = ft0.channelC()[i_c];
        FillConditionHistosBCs("hAmpT0C", amplitudeC);
      }

      multFT0A = ft0.sumAmpA();
      multFT0C = ft0.sumAmpC();
      multFT0M = multFT0A + multFT0C;

      if (multFT0A > 0) {
        FillConditionHistosBCs("hMultT0A", multFT0A);
      }
      if (multFT0C > 0) {
        FillConditionHistosBCs("hMultT0C", multFT0C);
      }
      if (multFT0M > 0) {
        FillConditionHistosBCs("hMultT0AC", multFT0M);
      }

      // Multiplicity correlations with ZDC
      if (hasZDC) {
        multZNA = bc.zdc().energyCommonZNA();
        multZNC = bc.zdc().energyCommonZNC();

        FillConditionHistosBCs("hMultT0A_ZNA", multFT0A, multZNA);
        FillConditionHistosBCs("hMultT0C_ZNC", multFT0C, multZNC);

      } // end if(hasZDC)

      // is FV0
      if (hasFV0) {
        auto fv0 = bc.fv0a();

        FillConditionHistosBCs("hT0A_V0A", fv0.time(), ft0.timeA());

      } // end of if (hasFV0)

    } // bc.hasFT0()
      //

    // Check only ZDC (without FIT)
    bool isZEM1 = false;
    bool isZEM2 = false;

    if (hasZDC) {

      isZEM1 = bc.zdc().timeZEM1() > -std::numeric_limits<float>::infinity();
      isZEM2 = bc.zdc().timeZEM2() > -std::numeric_limits<float>::infinity();

      multZNA = bc.zdc().energyCommonZNA();
      multZNC = bc.zdc().energyCommonZNC();

      histos.fill(HIST("BCs/hMultZNA_ZNC/All"), multZNA, multZNC);

      if (!isZEM1 && !isZEM2) {

        histos.fill(HIST("BCs/hMultZNA_ZNC_EM/All"), multZNA, multZNC);
      }
    }

  } //  end of processBCs()

  PROCESS_SWITCH(qaFT0TrgBased, processBCs, "per-BC analysis", true);
}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qaFT0TrgBased>(cfgc, TaskName{"ft0-qa-trg-based"})};
}
