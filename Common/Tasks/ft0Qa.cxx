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

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"

#include <DataFormatsFT0/Digit.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <bitset>
#include <cstddef>

using namespace o2;
using namespace o2::framework;

using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::BcSels, aod::Timestamps,
                                       aod::Run3MatchedToBCSparse>;

struct ft0QaTask {

  // Configurable<bool> isMC{"isMC", 0, "0 - data, 1 - MC"};
  Configurable<int> selection{"selection", 0, "trigger: 0 - no sel, 8 - sel8"};
  Configurable<bool> isLowFlux{"isLowFlux", 1,
                               "1 - low flux (pp, pPb), 0 - high flux (PbPb)"};

  HistogramRegistry histos{
    "Histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {

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

    const AxisSpec axisMultT0A{isLowFlux ? 10000 : 200000, 0.,
                               isLowFlux ? 10000. : 200000.,
                               "T0A multiplicity (# ADC channels)"};
    const AxisSpec axisMultT0C{isLowFlux ? 2000 : 70000, 0.,
                               isLowFlux ? 2000. : 70000.,
                               "T0C multiplicity (# ADC channels)"};
    const AxisSpec axisMultT0AC{isLowFlux ? 12000 : 27000, 0.,
                                isLowFlux ? 12000. : 270000.,
                                "T0AC multiplicity (# ADC channels)"};
    const AxisSpec axisMultV0A{isLowFlux ? 40000 : 200000, 0.,
                               isLowFlux ? 40000. : 200000.,
                               "V0A multiplicity (# ADC channels)"};

    const AxisSpec axisNcontrib{isLowFlux ? 150 : 4500, 0.,
                                isLowFlux ? 150. : 4500., "# contributors"};

    const AxisSpec axisCounter{2, 0.5, 1.5, ""};

    // Collisions task
    //  FT0 time
    histos.add("hT0A", "T0A;T0A time (ns);counts", kTH1F, {axisTime});
    histos.add("hT0C", "T0C;T0C time (ns);counts", kTH1F, {axisTime});
    histos.add("hT0AC", "T0AC;T0AC time (ns);counts", kTH1F, {axisTime});
    histos.add("hT0res", "FT0 resolution", kTH1F, {axisColTimeRes});
    histos.add("hColTime", "", kTH1F, {axisTime});

    // FT0 vertex
    histos.add("hT0vertex", "FT0 vertex;FT0 vertex (cm);counts", kTH1F,
               {axisVertex});
    histos.add("hPV", "PV;primary vertex (cm);counts", kTH1F, {axisVertex});
    histos.add("hT0vertexDiff", "FT0V - PV;FT0 vertex -  PV (cm);counts", kTH1F,
               {axisVertex});
    histos.add("hVertex_T0_PV",
               "PV vs. FT0V;FT0 vertex (cm);primary vertex (cm)", kTH2F,
               {axisVertex, axisVertex});
    histos.add("hVertex_T0_PV_nC20",
               "PV vs. FT0V;FT0 vertex (cm);primary vertex (cm)", kTH2F,
               {axisVertex, axisVertex});
    histos.add("hT0vertexDiff_nC20", "FT0V - PV;FT0 vertex -  PV (cm);counts",
               kTH1F, {axisVertex});
    histos.add("hVertex_T0_PV_TVX",
               "PV vs. FT0V;FT0 vertex (cm);primary vertex (cm)", kTH2F,
               {axisVertex, axisVertex});
    histos.add("hT0AC_T0Vertex_TVX",
               "time vs. FT0V;FT0 vertex (cm);T0AC time (ns)", kTH2F,
               {axisVertex, axisTime});

    histos.add("hPV_nContrib",
               "PV vs. Ncontributers;primary vertex (cm);(# contrubutors)",
               kTH2F, {axisVertex, axisNcontrib});

    // FT0 amplitude and multiplicity
    histos.add("hAmpT0A", "amplitude T0A;#ADC channels;counts", kTH1F,
               {axisAmpT0A});
    histos.add("hAmpT0C", "amplitude T0C;#ADC channels;counts", kTH1F,
               {axisAmpT0C});

    histos.add("hMultT0A", ";(# ADC channels);counts", kTH1F, {axisMultT0A});
    histos.add("hMultT0AOrA", "OrA trig;(# ADC channels);counts", kTH1F,
               {axisMultT0A});
    histos.add("hMultT0ATVX", "TVX trig;(# ADC channels);counts", kTH1F,
               {axisMultT0A});
    histos.add("hMultT0AOrAC", "OrA && OrC trig;(# ADC channels);counts", kTH1F,
               {axisMultT0A});
    histos.add("hMultT0ACent", "Central trig;(# ADC channels);counts", kTH1F,
               {axisMultT0A});
    histos.add("hMultT0ASCent", "SemiCentral trig;(# ADC channels);counts",
               kTH1F, {axisMultT0A});

    histos.add("hMultT0C", ";(# ADC channels);counts", kTH1F, {axisMultT0C});
    histos.add("hMultT0COrC", "OrC trig;(# ADC channels);counts", kTH1F,
               {axisMultT0C});
    histos.add("hMultT0CTVX", "TVX trig;(# ADC channels);counts", kTH1F,
               {axisMultT0C});
    histos.add("hMultT0COrAC", "OrA && OrC trig;(# ADC channels);counts", kTH1F,
               {axisMultT0C});
    histos.add("hMultT0CCent", "Central trig;(# ADC channels);counts", kTH1F,
               {axisMultT0C});
    histos.add("hMultT0CSCent", "SemiCentral trig;(# ADC channels);counts",
               kTH1F, {axisMultT0C});

    histos.add("hMultT0AC", ";(# ADC channels);counts", kTH1F, {axisMultT0AC});
    histos.add("hMultT0ACTVX", "TVX trig;(# ADC channels);counts", kTH1F,
               {axisMultT0AC});
    histos.add("hMultT0ACOrAC", "OrA && OrC trig;(# ADC channels);counts",
               kTH1F, {axisMultT0AC});
    histos.add("hMultT0ACCent", "Central trig;(# ADC channels);counts", kTH1F,
               {axisMultT0AC});
    histos.add("hMultT0ACSCent", "SemiCentral trig;(# ADC channels);counts",
               kTH1F, {axisMultT0AC});

    // Number of contributers to PV for FT0 triggers
    histos.add("hContrib", ";(# contrubutors);counts", kTH1F, {axisNcontrib});
    histos.add("hContribOrAC",
               "Ncontributers with T0AC;(# contrubutors);counts", kTH1F,
               {axisNcontrib});
    histos.add("hContribOrA", "Ncontributers with T0A;(# contrubutors);counts",
               kTH1F, {axisNcontrib});
    histos.add("hContribOrC", "Ncontributers with T0C;(# contrubutors);counts",
               kTH1F, {axisNcontrib});
    histos.add("hContribTVX", "Ncontributers with TVX;(# contrubutors);counts",
               kTH1F, {axisNcontrib});
    histos.add("hContribCent",
               "Ncontributers with bitCen;(# contrubutors);counts", kTH1F,
               {axisNcontrib});
    histos.add("hContribSemiCent",
               "Ncontributers with bitSCen;(# contrubutors);counts", kTH1F,
               {axisNcontrib});

    // FV0
    histos.add("hContribV0", "Ncontributers with V0A", kTH1F, {axisNcontrib});
    histos.add("hT0V0time", "T0A vs V0 time;V0 time (ns);T0A time (ns)", kTH2F,
               {axisTime, axisTime});

    //_______________________________________________________________________________________________________________________
    // BC task

    histos.add("hBcCounterAll", "", kTH1F, {axisCounter});
    histos.add("hBcCounterFT0", "", kTH1F, {axisCounter});
    histos.add("hBcCounterTVX", "", kTH1F, {axisCounter});

    histos.add("hBcT0A", "T0A;T0A time (ns);counts", kTH1F, {axisTime});
    histos.add("hBcT0C", "T0C;T0C time (ns);counts", kTH1F, {axisTime});
    histos.add("hBcT0AC", "T0AC;T0AC time (ns);counts", kTH1F, {axisTime});
    histos.add("hBcT0ACdiff", "(T0A -T0C)/2;(T0A -T0C)/2 (ns);counts", kTH1F,
               {axisColTimeRes});

    histos.add("hBcT0ATVX", "T0A TVX;T0A time (ns);counts", kTH1F, {axisTime});
    histos.add("hBcT0CTVX", "T0C TVX;T0C time (ns);counts", kTH1F, {axisTime});
    histos.add("hBcT0ACTVX", "T0AC TVX;T0AC time (ns);counts", kTH1F,
               {axisTime});
    histos.add("hBcT0ACdiffTVX", "(T0A -T0C)/2;(T0A -T0C)/2 (ns);counts", kTH1F,
               {axisColTimeRes});

    // FT0 amplitude and multiplicity
    histos.add("hBcAmpT0A", "amplitude T0A;#ADC channels;counts", kTH1F,
               {axisAmpT0A});
    histos.add("hBcAmpT0C", "amplitude T0C;#ADC channels;counts", kTH1F,
               {axisAmpT0C});

    histos.add("hBcMultT0A", ";(# ADC channels);counts", kTH1F, {axisMultT0A});
    histos.add("hBcMultT0AOrA", "OrA trig;(# ADC channels);counts", kTH1F,
               {axisMultT0A});
    histos.add("hBcMultT0ATVX", "TVX trig;(# ADC channels);counts", kTH1F,
               {axisMultT0A});
    histos.add("hBcMultT0AOrAC", "OrA && OrC trig;(# ADC channels);counts",
               kTH1F, {axisMultT0A});
    histos.add("hBcMultT0ACent", "Central trig;(# ADC channels);counts", kTH1F,
               {axisMultT0A});
    histos.add("hBcMultT0ASCent", "SemiCentral trig;(# ADC channels);counts",
               kTH1F, {axisMultT0A});

    histos.add("hBcMultT0C", ";(# ADC channels);counts", kTH1F, {axisMultT0C});
    histos.add("hBcMultT0COrC", "OrC trig;(# ADC channels);counts", kTH1F,
               {axisMultT0C});
    histos.add("hBcMultT0CTVX", "TVX trig;(# ADC channels);counts", kTH1F,
               {axisMultT0C});
    histos.add("hBcMultT0COrAC", "OrA && OrC trig;(# ADC channels);counts",
               kTH1F, {axisMultT0C});
    histos.add("hBcMultT0CCent", "Central trig;(# ADC channels);counts", kTH1F,
               {axisMultT0C});
    histos.add("hBcMultT0CSCent", "SemiCentral trig;(# ADC channels);counts",
               kTH1F, {axisMultT0C});

    histos.add("hBcMultT0AC", ";(# ADC channels);counts", kTH1F,
               {axisMultT0AC});
    histos.add("hBcMultT0ACTVX", "TVX trig;(# ADC channels);counts", kTH1F,
               {axisMultT0AC});
    histos.add("hBcMultT0ACOrAC", "OrA && OrC trig;(# ADC channels);counts",
               kTH1F, {axisMultT0AC});
    histos.add("hBcMultT0ACCent", "Central trig;(# ADC channels);counts", kTH1F,
               {axisMultT0AC});
    histos.add("hBcMultT0ACSCent", "SemiCentral trig;(# ADC channels);counts",
               kTH1F, {axisMultT0AC});
  }

  void
    processCollisions(soa::Join<aod::Collisions, aod::EvSels, aod::Mults,
                                aod::FT0sCorrected>::iterator const& collision,
                      aod::FT0s const&, aod::FV0As const&)
  {

    if (selection == 8 && !collision.sel8()) {
      return;
    }

    float multFT0A = 0.f;
    float multFT0C = 0.f;
    float multFT0M = 0.f;

    bool ora = false;
    bool orc = false;
    bool tvx = false;
    bool cent = false;
    bool semicent = false;

    // number of contributers used for vertex calculation
    int nContrib = collision.numContrib();
    histos.fill(HIST("hContrib"), nContrib);

    histos.fill(HIST("hPV_nContrib"), collision.posZ(), nContrib);

    //  collision time
    float colTime = collision.collisionTime();
    histos.fill(HIST("hColTime"), colTime);

    // is FT0
    if (collision.has_foundFT0()) {

      auto ft0 = collision.foundFT0();

      std::bitset<8> triggers = ft0.triggerMask();
      ora = triggers[o2::ft0::Triggers::bitA];
      orc = triggers[o2::ft0::Triggers::bitC];
      tvx = triggers[o2::ft0::Triggers::bitVertex];
      cent = triggers[o2::ft0::Triggers::bitCen];
      semicent = triggers[o2::ft0::Triggers::bitSCen];

      // FT0 multiplicity calculation
      if (tvx) {
        histos.fill(HIST("hContribTVX"), nContrib);
      }
      if (ora) {
        histos.fill(HIST("hContribOrA"), nContrib);
      }
      if (orc) {
        histos.fill(HIST("hContribOrC"), nContrib);
      }
      if (ora && orc) {
        histos.fill(HIST("hContribOrAC"), nContrib);
      }
      if (cent) {
        histos.fill(HIST("hContribCent"), nContrib);
      }
      if (semicent) {
        histos.fill(HIST("hContribSemiCent"), nContrib);
      }

      for (std::size_t i_a = 0; i_a < ft0.amplitudeA().size(); i_a++) {

        float amplitudeA = ft0.amplitudeA()[i_a];
        //   uint8_t channel = ft0.channelA()[i_a];

        histos.fill(HIST("hAmpT0A"), amplitudeA);
      }

      for (std::size_t i_c = 0; i_c < ft0.amplitudeC().size(); i_c++) {
        float amplitudeC = ft0.amplitudeC()[i_c];
        //  uint8_t channel = ft0.channelC()[i_c];
        histos.fill(HIST("hAmpT0C"), amplitudeC);
      }

      multFT0A = collision.multFT0A();
      multFT0C = collision.multFT0C();
      multFT0M = multFT0A + multFT0C;

      // Multiplicities (no triggers)

      histos.fill(HIST("hMultT0A"), multFT0A);
      histos.fill(HIST("hMultT0C"), multFT0C);

      histos.fill(HIST("hMultT0AC"), multFT0M);

      // Multiplicities (incl. triggers)
      if (tvx) {
        histos.fill(HIST("hMultT0ATVX"), multFT0A);
        histos.fill(HIST("hMultT0CTVX"), multFT0C);
        histos.fill(HIST("hMultT0ACTVX"), multFT0M);
      }

      if (ora) {
        histos.fill(HIST("hMultT0AOrA"), multFT0A);
      }

      if (orc) {
        histos.fill(HIST("hMultT0COrC"), multFT0C);
      }

      if (ora && orc) {
        histos.fill(HIST("hMultT0AOrAC"), multFT0A);
        histos.fill(HIST("hMultT0COrAC"), multFT0C);
        histos.fill(HIST("hMultT0ACOrAC"), multFT0M);
      }

      if (cent) {
        histos.fill(HIST("hMultT0ACent"), multFT0A);
        histos.fill(HIST("hMultT0CCent"), multFT0C);
        histos.fill(HIST("hMultT0ACCent"), multFT0M);
      }

      if (semicent) {
        histos.fill(HIST("hMultT0ASCent"), multFT0A);
        histos.fill(HIST("hMultT0CSCent"), multFT0C);
        histos.fill(HIST("hMultT0ACSCent"), multFT0M);
      }

      // vertex corrected FT0A and FT0C times
      if (collision.t0ACorrectedValid()) {
        histos.fill(HIST("hT0A"), collision.t0ACorrected());
      }
      if (collision.t0CCorrectedValid()) {
        histos.fill(HIST("hT0C"), collision.t0CCorrected());
      }

      if (collision.t0CCorrectedValid() && collision.t0ACorrectedValid()) {
        histos.fill(HIST("hT0AC"), collision.t0AC());
        histos.fill(HIST("hT0vertex"), ft0.posZ());
        histos.fill(HIST("hVertex_T0_PV"), ft0.posZ(), collision.posZ());
        histos.fill(HIST("hPV"), collision.posZ());
        histos.fill(HIST("hT0res"), collision.t0resolution());
        histos.fill(HIST("hT0vertexDiff"), ft0.posZ() - collision.posZ());

        if (nContrib > 20) {

          histos.fill(HIST("hVertex_T0_PV_nC20"), ft0.posZ(), collision.posZ());
          histos.fill(HIST("hT0vertexDiff_nC20"),
                      ft0.posZ() - collision.posZ());
        }

        if (tvx) {
          histos.fill(HIST("hVertex_T0_PV_TVX"), ft0.posZ(), collision.posZ());
          histos.fill(HIST("hT0AC_T0Vertex_TVX"), ft0.posZ(), collision.t0AC());
        }
      }

    } // end of if (collision.has_foundFT0())

    // is FV0
    if (collision.has_foundFV0()) {
      auto fv0 = collision.foundFV0();

      histos.fill(HIST("hContribV0"), nContrib);

      if (collision.has_foundFT0()) {
        auto ft0 = collision.foundFT0();
        histos.fill(HIST("hT0V0time"), fv0.time(), ft0.timeA());
      }
    } // end of if (collision.has_foundFV0())

  } // end of processCollsions()

  PROCESS_SWITCH(ft0QaTask, processCollisions, "per-collision analysis", true);

  // soa::Join<aod::BCs, aod::BcSels, aod::Mults >::iterator const &collision,
  // aod::FT0s const &ft0s, aod::FV0As const &fv0s

  void processBCs(BCsWithRun3Matchings::iterator const& bc, aod::FV0As const&,
                  aod::FT0s const&)
  {

    histos.fill(HIST("hBcCounterAll"), 1);

    // float multFV0A = 0.f;
    float multFT0A = 0.f;
    float multFT0C = 0.f;
    float multFT0M = 0.f;

    float timeA = -999.0;
    float timeC = -999.0;
    float timeAC = -999.0;
    float timeACdiff = -999.0;

    bool ora = false;
    bool orc = false;
    bool tvx = false;
    bool cent = false;
    bool semicent = false;

    if (bc.has_ft0()) {

      histos.fill(HIST("hBcCounterFT0"), 1);
      auto ft0 = bc.ft0();

      std::bitset<8> triggers = ft0.triggerMask();
      ora = triggers[o2::ft0::Triggers::bitA];
      orc = triggers[o2::ft0::Triggers::bitC];
      tvx = triggers[o2::ft0::Triggers::bitVertex];
      cent = triggers[o2::ft0::Triggers::bitCen];
      semicent = triggers[o2::ft0::Triggers::bitSCen];

      for (auto amplitude : ft0.amplitudeA()) {
        multFT0A += amplitude;
        histos.fill(HIST("hBcAmpT0A"), amplitude);
      }
      for (auto amplitude : ft0.amplitudeC()) {
        multFT0C += amplitude;
        histos.fill(HIST("hBcAmpT0C"), amplitude);
      }

      multFT0M = multFT0A + multFT0C;

      // Multiplicities (no triggers)

      histos.fill(HIST("hBcMultT0A"), multFT0A);
      histos.fill(HIST("hBcMultT0C"), multFT0C);

      histos.fill(HIST("hBcMultT0AC"), multFT0M);

      timeA = ft0.timeA();
      timeC = ft0.timeC();
      timeAC = (ft0.timeA() + ft0.timeC()) / 2;
      timeACdiff = (ft0.timeA() - ft0.timeC()) / 2;

      if (ora) {
        histos.fill(HIST("hBcMultT0AOrA"), multFT0A);

        histos.fill(HIST("hBcT0A"), timeA);
      }

      if (orc) {
        histos.fill(HIST("hBcMultT0COrC"), multFT0C);

        histos.fill(HIST("hBcT0C"), timeC);
      }

      if (ora && orc) {
        histos.fill(HIST("hBcMultT0AOrAC"), multFT0A);
        histos.fill(HIST("hBcMultT0COrAC"), multFT0C);
        histos.fill(HIST("hBcMultT0ACOrAC"), multFT0M);

        histos.fill(HIST("hBcT0AC"), timeAC);
        histos.fill(HIST("hBcT0ACdiff"), timeACdiff);
      }

      if (tvx) {
        //    continue;
        histos.fill(HIST("hBcCounterTVX"), 1);
        histos.fill(HIST("hBcMultT0ACTVX"), multFT0M);
        histos.fill(HIST("hBcMultT0ATVX"), multFT0A);
        histos.fill(HIST("hBcMultT0CTVX"), multFT0C);

        histos.fill(HIST("hBcT0ATVX"), timeA);
        histos.fill(HIST("hBcT0CTVX"), timeC);
        histos.fill(HIST("hBcT0ACTVX"), timeAC);
        histos.fill(HIST("hBcT0ACdiffTVX"), timeACdiff);
      }

      if (cent) {
        histos.fill(HIST("hBcMultT0ACent"), multFT0A);
        histos.fill(HIST("hBcMultT0CCent"), multFT0C);
        histos.fill(HIST("hBcMultT0ACCent"), multFT0M);
      }

      if (semicent) {
        histos.fill(HIST("hBcMultT0ASCent"), multFT0A);
        histos.fill(HIST("hBcMultT0CSCent"), multFT0C);
        histos.fill(HIST("hBcMultT0ACSCent"), multFT0M);
      }
    }

    // if (bc.has_fv0a()) {
    //   auto fv0a = bc.fv0a();
    //   for (auto amplitude : fv0a.amplitude()) {
    //     multFV0A += amplitude;
    //   }
    // }
  }
  PROCESS_SWITCH(ft0QaTask, processBCs, "per-BC analysis", true);

}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<ft0QaTask>(cfgc, TaskName{"ft0-qa"})};
}
