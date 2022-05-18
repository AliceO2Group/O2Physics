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
// This code calculates output histograms for centrality calibration
// as well as vertex-Z dependencies of raw variables (either for calibration
// of vtx-Z dependencies or for the calibration of those).
//
// This task is not strictly necessary in a typical analysis workflow,
// except for centrality calibration! The necessary task is the multiplicity
// tables.

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "TH1F.h"
#include "TH2F.h"
using namespace o2;
using namespace o2::framework;

struct FT0Qa {
  OutputObj<TH1F> hT0A{TH1F("hT0A", "T0A; ns", 200, -2, 2)};
  OutputObj<TH1F> hT0C{TH1F("hT0C", "T0C; ns", 200, -2, 2)};
  OutputObj<TH1F> hT0AC{TH1F("hT0AC", "(T0C+T0A)/2; ns", 200, -2, 2)};
  OutputObj<TH1F> hT0res{TH1F("hT0res", "T0 resolution; ns", 200, -0.5, 0.5)};
  OutputObj<TH1F> hT0Vertex{TH1F("hT0vertex", "T0 vertex; cm", 200, -30, 30.)};
  OutputObj<TH1F> hPV{TH1F("hPV", "Primary vertex; cm", 200, -30, 30.)};
  OutputObj<TH1F> hT0VertexDiff{TH1F("hT0vertexDiff", "T0 vertex -  PV; cm", 200, -30, 30.)};
  OutputObj<TH2F> hVertex_T0_PV{TH2F(" hVertex_T0_PV", "T0 vertex vs Primary Vertex; T0 vertex, cm; PV, cm", 200, -30, 30., 200, -30, 30.)};
  OutputObj<TH2F> hResMult{TH2F("hResMult", "T0 resolution vs event multiplicity", 100, 0., 100, 100, -0.5, 0.5)};
  OutputObj<TH2F> hResColMult{TH2F("hResCol", "Col. time resolution vs event multiplicity", 100, 0., 100, 100, -0.5, 0.5)};
  OutputObj<TH1F> hColTime{TH1F("hColTime", "ColTime, ns; Coiilisions, ns", 500, -5, 5.)};
  OutputObj<TH2F> hT0V0mult{TH2F("hT0V0mult", "T0A vs V0 multiplicity;V0Mult, #ADC channels;T0Mmult, #ADC channels", 200, 0., 500., 300, 0., 1500.)};
  OutputObj<TH2F> hT0V0time{TH2F("hT0V0time", "T0A vs V0 time ;V0time;T0A", 200, -2, 2., 200, -2, 2)};
  OutputObj<TH1F> hNcontrib{TH1F("hContrib", "Ncontributers;#contributers", 100, -0.5, 99.5)};
  OutputObj<TH1F> hNcontribAC{TH1F("hContribAC", "Ncontributers with T0AC;#contributers", 100, -0.5, 99.5)};
  OutputObj<TH1F> hNcontribA{TH1F("hContribA", "Ncontributers with T0A;#contributers", 100, -0.5, 99.5)};
  OutputObj<TH1F> hNcontribC{TH1F("hContribC", "Ncontributers with T0C;#contributers", 100, -0.5, 99.5)};
  OutputObj<TH1F> hNcontribV0{TH1F("hContribV0", "Ncontributers with V0A;#contributers", 100, -0.5, 99.5)};
  OutputObj<TH1F> hAmpV0{TH1F("hAmpV0", "amplitude V0A;#ADC channels", 500, 0, 3000)};
  OutputObj<TH1F> hAmpT0{TH1F("hAmpT0", "amplitude T0A;#ADC channels", 500, 0, 500)};

  // Configurable<bool> isMC{"isMC", 0, "0 - data, 1 - MC"};

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::FT0sCorrected>::iterator const& col, aod::FT0s const& ft0s, aod::FV0As const& fv0s)
  {
    float sumAmpFT0 = 0;
    float sumAmpFV0 = 0;
    int innerFV0 = 24;
    float fv0Time = -5000;
    uint16_t nContrib = col.numContrib();
    hNcontrib->Fill(nContrib);
    float colTime = col.collisionTime();
    hColTime->Fill(float(colTime));
    if (col.has_foundFV0()) {
      auto fv0 = col.foundFV0();
      fv0Time = fv0.time();
      hNcontribV0->Fill(nContrib);
      for (std::size_t ich = 0; ich < fv0.amplitude().size(); ich++) {
        hAmpV0->Fill(fv0.amplitude()[ich]);
        if (int(fv0.channel()[ich]) > innerFV0)
          continue;
        sumAmpFV0 += fv0.amplitude()[ich];
      }
    }
    if (col.has_foundFT0()) {
      auto ft0 = col.foundFT0();
      if (col.t0ACorrectedValid()) {
        hT0A->Fill(col.t0ACorrected());
        hNcontribA->Fill(nContrib);
        if (col.has_foundFV0()) {
          hT0V0time->Fill(fv0Time, ft0.timeA());
        }
        for (auto amplitude : ft0.amplitudeA()) {
          sumAmpFT0 += amplitude;
          hAmpT0->Fill(amplitude);
        }
      }
      if (col.has_foundFV0()) {
        if (sumAmpFV0 > 0 && sumAmpFT0 > 0) {
          hT0V0mult->Fill(sumAmpFV0, sumAmpFT0);
          LOG(debug) << "V0 amp  " << sumAmpFV0 << " T0 amp " << sumAmpFT0;
        }
      }
      if (col.t0CCorrectedValid()) {
        hT0C->Fill(col.t0CCorrected());
        hNcontribC->Fill(nContrib);
      }
      if (col.t0CCorrectedValid() && col.t0ACorrectedValid()) {
        LOGF(debug, "multFV0A=%f;  multFT0A=%f; PV=%f; T0vertex=%f; T0A=%f; T0C=%f; T0AC=%f; ColTime=%f", col.multFV0A(), col.multFT0A(), col.posZ(), ft0.posZ(), col.t0ACorrected(), col.t0CCorrected(), col.t0AC(), colTime);
        hT0AC->Fill(col.t0AC());
        hT0Vertex->Fill(ft0.posZ());
        hVertex_T0_PV->Fill(ft0.posZ(), col.posZ());
        hPV->Fill(col.posZ());
        hT0res->Fill(col.t0resolution());
        hT0VertexDiff->Fill(ft0.posZ() - col.posZ());
        hResMult->Fill(nContrib, (col.t0ACorrected() - col.t0CCorrected()) / 2.);
        hResColMult->Fill(nContrib, colTime);
        hNcontribAC->Fill(nContrib);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FT0Qa>(cfgc, TaskName{"ft0-qa"})};
}
