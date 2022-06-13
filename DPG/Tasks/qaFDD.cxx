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
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "DataFormatsFDD/Digit.h"
#include "DataFormatsFIT/Triggers.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace o2;
using namespace o2::framework;
int nBCsPerOrbit = 3564;
struct fddQA {
  Configurable<double> minOrbit{"minOrbit", 0, "minimum orbit"};
  Configurable<int> nOrbits{"nOrbits", 10000, "number of orbits"};
  Configurable<int> refBC{"refBC", 1238, "reference bc"};
  Configurable<int> nOrbitsPerTF{"nOrbitsPerTF", 256, "reference bc"};

  OutputObj<TH1F> hBcCol{TH1F("hBcCol", ";;", nBCsPerOrbit, 0., double(nBCsPerOrbit))};
  OutputObj<TH1F> hBcFDD{TH1F("hBcFDD", ";;", nBCsPerOrbit, 0., double(nBCsPerOrbit))};
  OutputObj<TH2F> hBcOrbitColl{TH2F("hBcOrbitColl", "Orbit vs BC [Collision];Orbit;BC", double(nOrbitsPerTF),
                                    0, double(nOrbitsPerTF), nBCsPerOrbit, 0., double(nBCsPerOrbit))};
  OutputObj<TH2F> hBcOrbitFDD{TH2F("hBcOrbitFDD", "Orbit vs BC [FDD];Orbit;BC", double(nOrbitsPerTF),
                                   0, double(nOrbitsPerTF), nBCsPerOrbit, 0., double(nBCsPerOrbit))};
  OutputObj<TH1F> hTFDDA{TH1F("hTFDDA", "Time (FDDA); ns", 2000, -20, 20)};
  OutputObj<TH1F> hTFDDC{TH1F("hTFDDC", " Time (FDDC); ns", 2000, -20, 20)};
  OutputObj<TH1F> hTFDDAC{TH1F("hTFDDAC", " Time (FDDA+FDDC)/2; ns", 2000, -20, 20)};
  OutputObj<TH1F> hChFDDA{TH1F("hChFDDA", "FDDA; Charge in ADC;", 5010, -10, 5000)};
  OutputObj<TH1F> hChFDDC{TH1F("hChFDDC", "FDDC; Charge in ADC;", 5010, -10, 5000)};
  OutputObj<TH1F> hTotalChargeFDDAC{TH1F("hTotalChargeFDDAC", "FDDC; Charge in ADC;", 8010, -10, 8000)};
  OutputObj<TH1F> hNcontribColl{TH1F("hNcontribColl", "Ncontributers in Coll TABLE;#contributers", 100, -0.5, 99.5)};
  OutputObj<TH1F> hNcontribFDD{TH1F("hNcontribFDD", "Ncontributers in FDD;#contributers", 100, -0.5, 99.5)};
  OutputObj<TH1F> hNcontribFDDAC{TH1F("hNcontribFDDAC", "Ncontributers in FDD A and C;#contributers", 100, -0.5, 99.5)};
  OutputObj<TH1F> hNcontribFDDorAC{TH1F("hNcontribFDDorAC", "Ncontributers in FDD A or C;#contributers", 100, -0.5, 99.5)};
  OutputObj<TH1F> hNcontribFDDA{TH1F("hNcontribFDDA", "Ncontributers in FDDA;#contributers", 100, -0.5, 99.5)};
  OutputObj<TH1F> hNcontribFDDC{TH1F("hNcontribFDDC", "Ncontributers with FDDC;#contributers", 100, -0.5, 99.5)};

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator const& collision, aod::FDDs const& fdds, aod::BCs const&)
  {
    float multFDDA = 0.f;
    float multFDDC = 0.f;
    float totalCharge = 0.f;
    auto bc = collision.bc_as<aod::BCs>();
    uint64_t globalBC = bc.globalBC();
    uint64_t orbit = globalBC % nOrbitsPerTF;
    int localBC = globalBC % nBCsPerOrbit;
    hBcCol->Fill(localBC);
    hBcOrbitColl->Fill(orbit, localBC);
    int nContributors = collision.numContrib();
    hNcontribColl->Fill(nContributors);
    if (collision.has_foundFDD()) {
      auto fdd = collision.foundFDD();
      hBcFDD->Fill(localBC);
      hBcOrbitFDD->Fill(orbit, localBC);
      std::bitset<8> fddTriggers = fdd.triggerMask();
      bool orA = fddTriggers[o2::fdd::Triggers::bitA];
      bool orC = fddTriggers[o2::fdd::Triggers::bitC];
      // bool vertex = fddTriggers[o2::fdd::Triggers::bitVertex];
      // bool central = fddTriggers[o2::fdd::Triggers::bitCen];
      // bool semiCentral = fddTriggers[o2::fdd::Triggers::bitSCen];
      // bool outputsAreBlocked = fddTriggers[o2::fdd::Triggers::bitOutputsAreBlocked];
      // bool dataIsValid = fddTriggers[o2::fdd::Triggers::bitDataIsValid];
      /*std::cout<<" orA "<< orA <<" orC "<< orC << " vertex "<<vertex<<" Central "<<central<<" SemiCentral "
                <<semiCentral<<" OutputsAreBlocked "<<outputsAreBlocked<<" DataIsValid "<<dataIsValid<<std::endl;*/
      hNcontribFDD->Fill(nContributors);
      if (orA) {
        hNcontribFDDA->Fill(nContributors);
        hTFDDA->Fill(fdd.timeA());
      }
      if (orC) {
        hNcontribFDDC->Fill(nContributors);
        hTFDDC->Fill(fdd.timeC());
      }
      if (orA && orC) {
        hNcontribFDDAC->Fill(nContributors);
        hTFDDAC->Fill((fdd.timeA() + fdd.timeC()) / 2.0);
      }
      if (orA || orC) {
        hNcontribFDDorAC->Fill(nContributors);
      }
      for (auto amplitude : fdd.chargeA()) {
        multFDDA += amplitude;
      }
      for (auto amplitude : fdd.chargeC()) {
        multFDDC += amplitude;
      }
      totalCharge = multFDDA + multFDDC;
      hChFDDA->Fill(multFDDA);
      hChFDDC->Fill(multFDDC);
      hTotalChargeFDDAC->Fill(totalCharge);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<fddQA>(cfgc, TaskName{"fdd-qa"})};
}