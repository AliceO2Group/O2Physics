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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "TH1F.h"
#include "TH2F.h"
#include "DataFormatsFT0/Digit.h"
#include <bitset>

using namespace o2;
using namespace o2::framework;

using CollissionsWithEvSelAndMultsAndFT0VertCorr = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::FT0sCorrected>;

struct FT0qanew {

  OutputObj<TH1F> hAmpT0A{TH1F("hAmpT0A", "amplitude T0A;#ADC channels", 10000, 0, 10000)};
  OutputObj<TH1F> hAmpT0C{TH1F("hAmpT0C", "amplitude T0C;#ADC channels", 10000, 0, 10000)};

  OutputObj<TH1F> hMultT0A{TH1F("hMultT0A", "multiplicity T0A; #ADC channels", 100100, -100, 100000)};
  OutputObj<TH1F> hMultT0AOrA{TH1F("hMultT0AOrA", "multiplicity T0A OrA; #ADC channels", 100100, -100, 100000)};
  OutputObj<TH1F> hMultT0ATVX{TH1F("hMultT0ATVX", "multiplicity T0A TVX; #ADC channels", 100100, -100, 100000)};
  OutputObj<TH1F> hMultT0AOrAC{TH1F("hMultT0AOrAC", "multiplicity T0A OrAC; #ADC channels", 100100, -100, 100000)};
  OutputObj<TH1F> hMultT0ACent{TH1F("hMultT0ACent", "multiplicity T0A Cent trig; #ADC channels", 100100, -100, 100000)};
  OutputObj<TH1F> hMultT0ASCent{TH1F("hMultT0ASCent", "multiplicity T0A SemiCent trig; #ADC channels", 100100, -100, 100000)};

  OutputObj<TH1F> hMultT0C{TH1F("hMultT0C", "multiplicity T0C; #ADC channels", 100100, -100, 100000)};
  OutputObj<TH1F> hMultT0COrC{TH1F("hMultT0COrC", "multiplicity T0C OrC; #ADC channels", 100100, -100, 100000)};
  OutputObj<TH1F> hMultT0CTVX{TH1F("hMultT0CTVX", "multiplicity T0C TVX; #ADC channels", 100100, -100, 100000)};
  OutputObj<TH1F> hMultT0COrAC{TH1F("hMultT0COrAC", "multiplicity T0C OrAC; #ADC channels", 100100, -100, 100000)};
  OutputObj<TH1F> hMultT0CCent{TH1F("hMultT0CCent", "multiplicity T0C Cent trig; #ADC channels", 100100, -100, 100000)};
  OutputObj<TH1F> hMultT0CSCent{TH1F("hMultT0CSCent", "multiplicity T0A SemiCent trig; #ADC channels", 100100, -100, 100000)};

  OutputObj<TH1F> hAmpT0A_3{TH1F("hAmpT0A_3", "amplitudeT0A; #ADC channels", 5000, 0, 10000)};
  OutputObj<TH1F> hAmpT0A_3_TVX{TH1F("hAmpT0A_3_TVX", "amplitudeT0A TVX; #ADC channels", 5000, 0, 10000)};
  OutputObj<TH1F> hAmpT0A_3_Cent{TH1F("hAmpT0A_3_Cent", "amplitudeT0A Cent trig; #ADC channels", 5000, 0, 10000)};
  OutputObj<TH1F> hAmpT0A_3_SCent{TH1F("hAmpT0A_3_SCent", "amplitudeT0A SemiCent trig; #ADC channels", 5000, 0, 10000)};

  OutputObj<TH1F> hAmpT0C_97{TH1F("hAmpT0C_97", "amplitudeT0C; #ADC channels", 5000, 0, 10000)};
  OutputObj<TH1F> hAmpT0C_97_TVX{TH1F("hAmpT0C_97_TVX", "amplitudeT0C TVX; #ADC channels", 5000, 0, 10000)};
  OutputObj<TH1F> hAmpT0C_97_Cent{TH1F("hAmpT0C_97_Cent", "amplitudeT0C Cent trig; #ADC channels", 5000, 0, 10000)};
  OutputObj<TH1F> hAmpT0C_97_SCent{TH1F("hAmpT0C_97_SCent", "amplitudeT0A SemiCent trig; #ADC channels", 5000, 0, 10000)};

  OutputObj<TH1F> hAmpT0A_17{TH1F("hAmpT0A_17", "amplitudeT0A; #ADC channels", 5000, 0, 10000)};
  OutputObj<TH1F> hAmpT0A_17_TVX{TH1F("hAmpT0A_17_TVX", "amplitudeT0A TVX; #ADC channels", 5000, 0, 10000)};
  OutputObj<TH1F> hAmpT0A_17_Cent{TH1F("hAmpT0A_17_Cent", "amplitudeT0A Cent trig; #ADC channels", 5000, 0, 10000)};
  OutputObj<TH1F> hAmpT0A_17_SCent{TH1F("hAmpT0A_17_SCent", "amplitudeT0A SemiCent trig; #ADC channels", 5000, 0, 10000)};

  OutputObj<TH1F> hAmpT0C_117{TH1F("hAmpT0C_117", "amplitudeT0C; #ADC channels", 5000, 0, 10000)};
  OutputObj<TH1F> hAmpT0C_117_TVX{TH1F("hAmpT0C_117_TVX", "amplitudeT0C TVX; #ADC channels", 5000, 0, 10000)};
  OutputObj<TH1F> hAmpT0C_117_Cent{TH1F("hAmpT0C_117_Cent", "amplitudeT0C Cent trig; #ADC channels", 5000, 0, 10000)};
  OutputObj<TH1F> hAmpT0C_117_SCent{TH1F("hAmpT0C_117_SCent", "amplitudeT0A SemiCent trig; #ADC channels", 5000, 0, 10000)};

  OutputObj<TH1F> hMultT0AC{TH1F("hMultT0AC", "multiplicity T0A+C; #ADC channels", 500100, -100, 500000)};
  OutputObj<TH1F> hMultT0ACOrAC{TH1F("hMultT0ACOrAC", "multiplicity T0A+C OrAC; #ADC channels", 500100, -100, 500000)};
  OutputObj<TH1F> hMultT0ACTVX{TH1F("hMultT0ACTVX", "multiplicity T0A+C TVX; #ADC channels", 500100, -100, 500000)};
  OutputObj<TH1F> hMultT0ACCent{TH1F("hMultT0ACCent", "multiplicity T0A+C Cent trig; #ADC channels", 500100, -100, 500000)};
  OutputObj<TH1F> hMultT0ACSCent{TH1F("hMultT0ACSCent", "multiplicity T0A+C SemiCent trig; #ADC channels", 500100, -100, 500000)};

  OutputObj<TH1F> hNcontribOrAC{TH1F("hContribOrAC", "Ncontributers with T0AC;#contributers", 120, -0.5, 119.5)};
  OutputObj<TH1F> hNcontribOrA{TH1F("hContribOrA", "Ncontributers with T0A;#contributers", 120, -0.5, 119.5)};
  OutputObj<TH1F> hNcontribOrC{TH1F("hContribOrC", "Ncontributers with T0C;#contributers", 120, -0.5, 119.5)};

  OutputObj<TH1F> hNcontribTVX{TH1F("hContribTVX", "Ncontributers with TVX;#contributers", 120, -0.5, 119.5)};
  OutputObj<TH1F> hNcontribCent{TH1F("hContribCent", "Ncontributers with bitCen;#contributers", 120, -0.5, 119.5)};
  OutputObj<TH1F> hNcontribSemiCent{TH1F("hContribSemiCent", "Ncontributers with bitSCen; #contributers", 120, -0.5, 119.5)};
  OutputObj<TH1F> hNcontrib{TH1F("hContrib", "Ncontributers;#contributers", 120, -0.5, 119.5)};

  OutputObj<TH1F> hT0A{TH1F("hT0A", "T0A; ns", 200, -2, 2)};
  OutputObj<TH1F> hT0C{TH1F("hT0C", "T0C; ns", 200, -2, 2)};
  OutputObj<TH1F> hT0AC{TH1F("hT0AC", "(T0C+T0A)/2; ns", 200, -2, 2)};

  OutputObj<TH1F> hT0res{TH1F("hT0res", "T0 resolution; (T0A - T0C)/2 (ns)", 400, -0.5, 0.5)};
  OutputObj<TH1F> hColTime{TH1F("hColTime", "Colission time (ns); Coiilisions; (ns)", 500, -5, 5.)};

  OutputObj<TH1F> hT0Vertex{TH1F("hT0vertex", "T0 vertex; (cm)", 200, -30, 30.)};
  OutputObj<TH1F> hPV{TH1F("hPV", "Primary vertex; (cm)", 200, -30, 30.)};
  OutputObj<TH1F> hT0VertexDiff{TH1F("hT0vertexDiff", "FT0 vertex -  PV; (cm)", 200, -30, 30.)};
  OutputObj<TH2F> hVertex_T0_PV{TH2F(" hVertex_T0_PV", "T0 vertex vs Primary Vertex; FT0 vertex (cm); PV (cm)", 200, -30, 30., 200, -30, 30.)};

  void process(CollissionsWithEvSelAndMultsAndFT0VertCorr const& collisions, aod::FV0As const&, aod::FT0s const&)
  {

    for (auto& collision : collisions) {

      //  if (!collision.sel8()) continue;

      float multFT0A = 0.f;
      float multFT0C = 0.f;
      float multFT0AC = 0.f;

      float ampA_3 = 0.f;
      float ampA_17 = 0.f;
      float ampC_97 = 0.f;
      float ampC_117 = 0.f;

      bool ora = false;
      bool orc = false;
      bool tvx = false;
      bool cent = false;
      bool semicent = false;

      // number of contributers used for vertex calculation
      uint16_t nContrib = collision.numContrib();
      hNcontrib->Fill(nContrib);

      //  collision time ?!
      float colTime = collision.collisionTime();
      hColTime->Fill(float(colTime));

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
          hNcontribTVX->Fill(nContrib);
        }
        if (ora) {
          hNcontribOrA->Fill(nContrib);
        }
        if (orc) {
          hNcontribOrC->Fill(nContrib);
        }
        if (ora && orc) {
          hNcontribOrAC->Fill(nContrib);
        }
        if (cent) {
          hNcontribCent->Fill(nContrib);
        }
        if (semicent) {
          hNcontribSemiCent->Fill(nContrib);
        }

        for (std::size_t i_a = 0; i_a < ft0.amplitudeA().size(); i_a++) {

          float amplitudeA = ft0.amplitudeA()[i_a];
          uint8_t channel = ft0.channelA()[i_a];

          hAmpT0A->Fill(amplitudeA);

          if (channel == 3) {
            ampA_3 += amplitudeA;
          }
          if (channel == 17) {
            ampA_17 += amplitudeA;
          }
        }

        for (std::size_t i_c = 0; i_c < ft0.amplitudeC().size(); i_c++) {
          float amplitudeC = ft0.amplitudeC()[i_c];
          uint8_t channel = ft0.channelC()[i_c];
          hAmpT0C->Fill(amplitudeC);

          if (channel == 97 - 96) {
            ampC_97 += amplitudeC;
          }
          if (channel == 117 - 96) {
            ampC_117 += amplitudeC;
          }
        }

        multFT0A = collision.multFT0A();
        multFT0C = collision.multFT0C();
        multFT0AC = multFT0A + multFT0C;

        // Multiplicities (no triggers)

        hMultT0A->Fill(multFT0A);
        hMultT0C->Fill(multFT0C);

        hMultT0AC->Fill(multFT0AC);

        if (ampA_3) {
          hAmpT0A_3->Fill(ampA_3);
        }
        if (ampA_17) {
          hAmpT0A_17->Fill(ampA_17);
        }
        if (ampC_97) {
          hAmpT0C_97->Fill(ampC_97);
        }
        if (ampC_117) {
          hAmpT0C_117->Fill(ampC_117);
        }

        // Multiplicities (incl. triggers)
        if (tvx) {
          hMultT0ATVX->Fill(multFT0A);
          hMultT0CTVX->Fill(multFT0C);
          hMultT0ACTVX->Fill(multFT0AC);

          if (ampA_3) {
            hAmpT0A_3_TVX->Fill(ampA_3);
          }
          if (ampA_17) {
            hAmpT0A_17_TVX->Fill(ampA_17);
          }
          if (ampC_97) {
            hAmpT0C_97_TVX->Fill(ampC_97);
          }
          if (ampC_117) {
            hAmpT0C_117_TVX->Fill(ampC_117);
          }
        }
        if (ora) {
          hMultT0AOrA->Fill(multFT0A);
        }

        if (orc) {
          hMultT0COrC->Fill(multFT0C);
        }
        if (ora && orc) {
          hMultT0AOrAC->Fill(multFT0A);
          hMultT0COrAC->Fill(multFT0C);
          hMultT0ACOrAC->Fill(multFT0AC);
        }

        if (cent) {
          hMultT0ACent->Fill(multFT0A);
          hMultT0CCent->Fill(multFT0C);
          hMultT0ACCent->Fill(multFT0AC);

          if (ampA_3) {
            hAmpT0A_3_Cent->Fill(ampA_3);
          }
          if (ampA_17) {
            hAmpT0A_17_Cent->Fill(ampA_17);
          }
          if (ampC_97) {
            hAmpT0C_97_Cent->Fill(ampC_97);
          }
          if (ampC_117) {
            hAmpT0C_117_Cent->Fill(ampC_117);
          }
        }

        if (semicent) {
          hMultT0ASCent->Fill(multFT0A);
          hMultT0CSCent->Fill(multFT0C);
          hMultT0ACSCent->Fill(multFT0AC);

          if (ampA_3) {
            hAmpT0A_3_SCent->Fill(ampA_3);
          }
          if (ampA_17) {
            hAmpT0A_17_SCent->Fill(ampA_17);
          }
          if (ampC_97) {
            hAmpT0C_97_SCent->Fill(ampC_97);
          }
          if (ampC_117) {
            hAmpT0C_117_SCent->Fill(ampC_117);
          }
        }

        // vertex corrected FT0A and FT0C times
        if (collision.t0ACorrectedValid()) {
          hT0A->Fill(collision.t0ACorrected());
        }
        if (collision.t0CCorrectedValid()) {
          hT0C->Fill(collision.t0CCorrected());
        }

        if (collision.t0CCorrectedValid() && collision.t0ACorrectedValid()) {
          hT0AC->Fill(collision.t0AC());
          hT0Vertex->Fill(ft0.posZ());
          hVertex_T0_PV->Fill(ft0.posZ(), collision.posZ());
          hPV->Fill(collision.posZ());
          hT0res->Fill(collision.t0resolution());
          hT0VertexDiff->Fill(ft0.posZ() - collision.posZ());
        }

      } // end of if (collision.has_foundFT0())

    } // end of collisions loop
  }   // end of procces()
};    // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FT0qanew>(cfgc, TaskName{"ft0-qa-ext"})};
}
