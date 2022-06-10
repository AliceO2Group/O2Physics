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
#include "Common/CCDB/EventSelectionParams.h"
#include <CCDB/BasicCCDBManager.h>
#include "Framework/HistogramRegistry.h"
#include "CommonDataFormat/BunchFilling.h"
#include "TH1F.h"
#include "TH2F.h"
using namespace o2::framework;
using namespace o2;
using namespace evsel;

using BCsRun2 = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps, aod::BcSels, aod::Run2MatchedToBCSparse>;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

struct EventSelectionQaTask {
  Configurable<bool> isMC{"isMC", 0, "0 - data, 1 - MC"};
  Configurable<double> minGlobalBC{"minGlobalBC", 0, "minimum global bc"};
  Configurable<int> nGlobalBCs{"nGlobalBCs", 100000, "number of global bcs"};
  Configurable<double> minOrbit{"minOrbit", 0, "minimum orbit"};
  Configurable<int> nOrbits{"nOrbits", 10000, "number of orbits"};
  Configurable<int> refBC{"refBC", 1238, "reference bc"};
  Configurable<bool> isLowFlux{"isLowFlux", 1, "1 - low flux (pp, pPb), 0 - high flux (PbPb)"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  bool* applySelection = NULL;
  int nBCsPerOrbit = 3564;
  int lastRunNumber = -1;

  std::bitset<o2::constants::lhc::LHCMaxBunches> beamPatternA;
  std::bitset<o2::constants::lhc::LHCMaxBunches> beamPatternC;
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternA;
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternC;
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternB;

  void init(InitContext&)
  {
    // ccdb->setURL("http://ccdb-test.cern.ch:8080");
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    float maxMultV0M = isLowFlux ? 40000 : 40000;
    float maxMultV0A = isLowFlux ? 30000 : 30000;
    float maxMultV0C = isLowFlux ? 30000 : 30000;
    float maxMultT0A = isLowFlux ? 2000 : 100000;
    float maxMultT0C = isLowFlux ? 2000 : 100000;
    float maxMultFDA = isLowFlux ? 5000 : 10000;
    float maxMultFDC = isLowFlux ? 5000 : 10000;
    float maxMultZNA = isLowFlux ? 1000 : 200000;
    float maxMultZNC = isLowFlux ? 1000 : 200000;
    const AxisSpec axisTime{700, -35, 35};
    histos.add("hTimeV0Aall", "All bcs;V0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeV0Call", "All bcs;V0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNAall", "All bcs;ZNA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNCall", "All bcs;ZNC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Aall", "All bcs;T0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Call", "All bcs;T0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDAall", "All bcs;FDA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDCall", "All bcs;FDC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZACall", "All bcs; ZNC-ZNA time (ns); ZNC+ZNA time (ns)", kTH2F, {{100, -5, 5}, {100, -5, 5}});
    histos.add("hTimeV0Abga", "BeamA-only bcs;V0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeV0Cbga", "BeamA-only bcs;V0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNAbga", "BeamA-only bcs;ZNA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNCbga", "BeamA-only bcs;ZNC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Abga", "BeamA-only bcs;T0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Cbga", "BeamA-only bcs;T0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDAbga", "BeamA-only bcs;FDA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDCbga", "BeamA-only bcs;FDC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeV0Abgc", "BeamC-only bcs;V0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeV0Cbgc", "BeamC-only bcs;V0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNAbgc", "BeamC-only bcs;ZNA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNCbgc", "BeamC-only bcs;ZNC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Abgc", "BeamC-only bcs;T0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Cbgc", "BeamC-only bcs;T0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDAbgc", "BeamC-only bcs;FDA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDCbgc", "BeamC-only bcs;FDC time (ns);Entries", kTH1F, {axisTime});

    histos.add("hTimeV0Aref", "Reference bcs;V0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeV0Cref", "Reference bcs;V0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNAref", "Reference bcs;ZNA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNCref", "Reference bcs;ZNC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Aref", "Reference bcs;T0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Cref", "Reference bcs;T0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDAref", "Reference bcs;FDA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDCref", "Reference bcs;FDC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZACref", "Reference bcs; ZNC-ZNA time (ns); ZNC+ZNA time (ns)", kTH2F, {{100, -5, 5}, {100, -5, 5}});
    histos.add("hTimeV0Acol", "All events;V0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeV0Ccol", "All events;V0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNAcol", "All events;ZNA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNCcol", "All events;ZNC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Acol", "All events;T0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Ccol", "All events;T0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDAcol", "All events;FDA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDCcol", "All events;FDC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZACcol", "All events; ZNC-ZNA time (ns); ZNC+ZNA time (ns)", kTH2F, {{100, -5, 5}, {100, -5, 5}});
    histos.add("hTimeV0Aacc", "Accepted events;V0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeV0Cacc", "Accepted events;V0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNAacc", "Accepted events;ZNA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNCacc", "Accepted events;ZNC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Aacc", "Accepted events;T0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Cacc", "Accepted events;T0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDAacc", "Accepted events;FDA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDCacc", "Accepted events;FDC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZACacc", "Accepted events; ZNC-ZNA time (ns); ZNC+ZNA time (ns)", kTH2F, {{100, -5, 5}, {100, -5, 5}});
    histos.add("hSPDClsVsTklCol", "All events;n tracklets;n clusters", kTH2F, {{200, 0., isLowFlux ? 200. : 6000.}, {100, 0., isLowFlux ? 100. : 20000.}});
    histos.add("hV0C012vsTklCol", "All events;n tracklets;V0C012 multiplicity", kTH2F, {{150, 0., 150.}, {150, 0., 600.}});
    histos.add("hV0MOnVsOfCol", "All events;Offline V0M;Online V0M", kTH2F, {{200, 0., isLowFlux ? 1000. : 50000.}, {400, 0., isLowFlux ? 8000. : 40000.}});
    histos.add("hSPDOnVsOfCol", "All events;Offline FOR;Online FOR", kTH2F, {{300, 0., isLowFlux ? 300. : 1200.}, {300, 0., isLowFlux ? 300. : 1200.}});
    histos.add("hV0C3vs012Col", "All events;V0C012 multiplicity;V0C3 multiplicity", kTH2F, {{200, 0., 800.}, {300, 0., 300.}});
    histos.add("hSPDClsVsTklAcc", "Accepted events;n tracklets;n clusters", kTH2F, {{200, 0., isLowFlux ? 200. : 6000.}, {100, 0., isLowFlux ? 100. : 20000.}});
    histos.add("hV0C012vsTklAcc", "Accepted events;n tracklets;V0C012 multiplicity", kTH2F, {{150, 0., 150.}, {150, 0., 600.}});
    histos.add("hV0MOnVsOfAcc", "Accepted events;Offline V0M;Online V0M", kTH2F, {{200, 0., isLowFlux ? 1000. : 50000.}, {400, 0., isLowFlux ? 8000. : 40000.}});
    histos.add("hSPDOnVsOfAcc", "Accepted events;Offline FOR;Online FOR", kTH2F, {{300, 0., isLowFlux ? 300. : 1200.}, {300, 0., isLowFlux ? 300. : 1200.}});
    histos.add("hV0C3vs012Acc", "Accepted events;V0C012 multiplicity;V0C3 multiplicity", kTH2F, {{200, 0., 800.}, {300, 0., 300.}});

    histos.add("hColCounterAll", "", kTH1F, {{kNaliases, 0, kNaliases}});
    histos.add("hColCounterAcc", "", kTH1F, {{kNaliases, 0, kNaliases}});
    histos.add("hBcCounterAll", "", kTH1F, {{kNaliases, 0, kNaliases}});
    histos.add("hSelCounter", "", kTH1F, {{kNsel, 0, kNsel}});
    histos.add("hSelMask", "", kTH1F, {{kNsel, 0, kNsel}});

    histos.add("hGlobalBcAll", ";;", kTH1F, {{nGlobalBCs, minGlobalBC, minGlobalBC + nGlobalBCs}});
    histos.add("hGlobalBcCol", ";;", kTH1F, {{nGlobalBCs, minGlobalBC, minGlobalBC + nGlobalBCs}});
    histos.add("hGlobalBcFT0", ";;", kTH1F, {{nGlobalBCs, minGlobalBC, minGlobalBC + nGlobalBCs}});
    histos.add("hGlobalBcFV0", ";;", kTH1F, {{nGlobalBCs, minGlobalBC, minGlobalBC + nGlobalBCs}});
    histos.add("hGlobalBcFDD", ";;", kTH1F, {{nGlobalBCs, minGlobalBC, minGlobalBC + nGlobalBCs}});
    histos.add("hOrbitAll", ";;", kTH1F, {{nOrbits, minOrbit, minOrbit + nOrbits}});
    histos.add("hOrbitCol", ";;", kTH1F, {{nOrbits, minOrbit, minOrbit + nOrbits}});
    histos.add("hOrbitFT0", ";;", kTH1F, {{nOrbits, minOrbit, minOrbit + nOrbits}});
    histos.add("hOrbitFV0", ";;", kTH1F, {{nOrbits, minOrbit, minOrbit + nOrbits}});
    histos.add("hOrbitFDD", ";;", kTH1F, {{nOrbits, minOrbit, minOrbit + nOrbits}});
    histos.add("hBcAll", ";;", kTH1F, {{nBCsPerOrbit, 0., double(nBCsPerOrbit)}});
    histos.add("hBcCol", ";;", kTH1F, {{nBCsPerOrbit, 0., double(nBCsPerOrbit)}});
    histos.add("hBcFT0", ";;", kTH1F, {{nBCsPerOrbit, 0., double(nBCsPerOrbit)}});
    histos.add("hBcFV0", ";;", kTH1F, {{nBCsPerOrbit, 0., double(nBCsPerOrbit)}});
    histos.add("hBcFDD", ";;", kTH1F, {{nBCsPerOrbit, 0., double(nBCsPerOrbit)}});

    histos.add("hMultV0Aall", "All bcs;V0A multiplicity;Entries", kTH1F, {{10000, 0., maxMultV0M}});
    histos.add("hMultV0Call", "All bcs;V0C multiplicity;Entries", kTH1F, {{10000, 0., maxMultV0C}});
    histos.add("hMultZNAall", "All bcs;ZNA multiplicity;Entries", kTH1F, {{1000, 0., maxMultZNA}});
    histos.add("hMultZNCall", "All bcs;ZNC multiplicity;Entries", kTH1F, {{1000, 0., maxMultZNC}});
    histos.add("hMultT0Aall", "All bcs;T0A multiplicity;Entries", kTH1F, {{1000, 0., maxMultT0A}});
    histos.add("hMultT0Call", "All bcs;T0C multiplicity;Entries", kTH1F, {{1000, 0., maxMultT0C}});
    histos.add("hMultFDAall", "All bcs;FDA multiplicity;Entries", kTH1F, {{1000, 0., maxMultFDA}});
    histos.add("hMultFDCall", "All bcs;FDC multiplicity;Entries", kTH1F, {{1000, 0., maxMultFDC}});
    histos.add("hMultV0Aref", "Reference bcs;V0A multiplicity;Entries", kTH1F, {{10000, 0., maxMultV0M}});
    histos.add("hMultV0Cref", "Reference bcs;V0C multiplicity;Entries", kTH1F, {{10000, 0., maxMultV0C}});
    histos.add("hMultZNAref", "Reference bcs;ZNA multiplicity;Entries", kTH1F, {{1000, 0., maxMultZNA}});
    histos.add("hMultZNCref", "Reference bcs;ZNC multiplicity;Entries", kTH1F, {{1000, 0., maxMultZNC}});
    histos.add("hMultT0Aref", "Reference bcs;T0A multiplicity;Entries", kTH1F, {{1000, 0., maxMultT0A}});
    histos.add("hMultT0Cref", "Reference bcs;T0C multiplicity;Entries", kTH1F, {{1000, 0., maxMultT0C}});
    histos.add("hMultFDAref", "Reference bcs;FDA multiplicity;Entries", kTH1F, {{1000, 0., maxMultFDA}});
    histos.add("hMultFDCref", "Reference bcs;FDC multiplicity;Entries", kTH1F, {{1000, 0., maxMultFDC}});
    histos.add("hMultV0Mcol", "All events;V0M multiplicity;Entries", kTH1F, {{10000, 0., maxMultV0M}});
    histos.add("hMultV0Acol", "All events;V0A multiplicity;Entries", kTH1F, {{10000, 0., maxMultV0A}});
    histos.add("hMultV0Ccol", "All events;V0C multiplicity;Entries", kTH1F, {{10000, 0., maxMultV0C}});
    histos.add("hMultZNAcol", "All events;ZNA multiplicity;Entries", kTH1F, {{1000, 0., maxMultZNA}});
    histos.add("hMultZNCcol", "All events;ZNC multiplicity;Entries", kTH1F, {{1000, 0., maxMultZNC}});
    histos.add("hMultT0Acol", "All events;T0A multiplicity;Entries", kTH1F, {{1000, 0., maxMultT0A}});
    histos.add("hMultT0Ccol", "All events;T0C multiplicity;Entries", kTH1F, {{1000, 0., maxMultT0C}});
    histos.add("hMultFDAcol", "All events;FDA multiplicity;Entries", kTH1F, {{1000, 0., maxMultFDA}});
    histos.add("hMultFDCcol", "All events;FDC multiplicity;Entries", kTH1F, {{1000, 0., maxMultFDC}});
    histos.add("hMultV0Macc", "Accepted events;V0M multiplicity;Entries", kTH1F, {{10000, 0., maxMultV0M}});
    histos.add("hMultV0Aacc", "Accepted events;V0A multiplicity;Entries", kTH1F, {{10000, 0., maxMultV0A}});
    histos.add("hMultV0Cacc", "Accepted events;V0C multiplicity;Entries", kTH1F, {{10000, 0., maxMultV0C}});
    histos.add("hMultZNAacc", "Accepted events;ZNA multiplicity;Entries", kTH1F, {{1000, 0., maxMultZNA}});
    histos.add("hMultZNCacc", "Accepted events;ZNC multiplicity;Entries", kTH1F, {{1000, 0., maxMultZNC}});
    histos.add("hMultT0Aacc", "Accepted events;T0A multiplicity;Entries", kTH1F, {{1000, 0., maxMultT0A}});
    histos.add("hMultT0Cacc", "Accepted events;T0C multiplicity;Entries", kTH1F, {{1000, 0., maxMultT0C}});
    histos.add("hMultFDAacc", "Accepted events;FDA multiplicity;Entries", kTH1F, {{1000, 0., maxMultFDA}});
    histos.add("hMultFDCacc", "Accepted events;FDC multiplicity;Entries", kTH1F, {{1000, 0., maxMultFDC}});

    histos.add("hColTimeRes", ";collision time resolution (ns)", kTH1F, {{7000, 0, 7000}});
    histos.add("hColTimeResVsNcontrib", ";n contributors; collision time resolution (ns)", kTH2F, {{100, 0, 100}, {7000, 0, 7000}});
    histos.add("hColTimeResVsNcontribITSonly", ";n contributors; collision time resolution (ns)", kTH2F, {{100, 0, 100}, {7000, 0, 7000}});
    histos.add("hColTimeResVsNcontribWithTOF", ";n contributors; collision time resolution (ns)", kTH2F, {{100, 0, 100}, {7000, 0, 7000}});
    histos.add("hColBcDiffVsNcontrib", ";n contributors; collision bc difference", kTH2F, {{100, 0, 100}, {600, -300, 300}});
    histos.add("hColBcDiffVsNcontribITSonly", ";n contributors; collision bc difference", kTH2F, {{100, 0, 100}, {600, -300, 300}});
    histos.add("hColBcDiffVsNcontribWithTOF", ";n contributors; collision bc difference", kTH2F, {{100, 0, 100}, {600, -300, 300}});
    histos.add("hDFstartOrbit", "", kTH1F, {{100000, 0., 1e+9}});
    histos.add("hFT0sPerDF", "", kTH1F, {{100000, 0., 1e+9}});

    histos.add("hNcontribCol", ";n contributors;", kTH1F, {{100, 0, 100.}});
    histos.add("hNcontribAcc", ";n contributors;", kTH1F, {{100, 0, 100.}});

    // MC histograms
    histos.add("hGlobalBcColMC", ";;", kTH1F, {{nGlobalBCs, minGlobalBC, minGlobalBC + nGlobalBCs}});
    histos.add("hOrbitColMC", ";;", kTH1F, {{nOrbits, minOrbit, minOrbit + nOrbits}});
    histos.add("hBcColMC", "", kTH1F, {{nBCsPerOrbit, 0., double(nBCsPerOrbit)}});
    histos.add("hVertexXMC", "", kTH1F, {{1000, -1, 1}});
    histos.add("hVertexYMC", "", kTH1F, {{1000, -1, 1}});
    histos.add("hVertexZMC", "", kTH1F, {{1000, -10, 10}});

    for (int i = 0; i < kNsel; i++) {
      histos.get<TH1>(HIST("hSelCounter"))->GetXaxis()->SetBinLabel(i + 1, selectionLabels[i]);
      histos.get<TH1>(HIST("hSelMask"))->GetXaxis()->SetBinLabel(i + 1, selectionLabels[i]);
    }
    for (int i = 0; i < kNaliases; i++) {
      histos.get<TH1>(HIST("hColCounterAll"))->GetXaxis()->SetBinLabel(i + 1, aliasLabels[i]);
      histos.get<TH1>(HIST("hColCounterAcc"))->GetXaxis()->SetBinLabel(i + 1, aliasLabels[i]);
      histos.get<TH1>(HIST("hBcCounterAll"))->GetXaxis()->SetBinLabel(i + 1, aliasLabels[i]);
    }
  }

  void processRun2(
    soa::Join<aod::Collisions, aod::EvSels> cols,
    BCsRun2 const& bcs,
    aod::Zdcs const& zdcs,
    aod::FV0As const& fv0as,
    aod::FV0Cs const& fv0cs,
    aod::FT0s const& ft0s,
    aod::FDDs const& fdds)
  {
    if (!applySelection) {
      auto first_bc = bcs.iteratorAt(0);
      EventSelectionParams* par = ccdb->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", first_bc.timestamp());
      applySelection = par->GetSelection(0);
      for (int i = 0; i < kNsel; i++) {
        histos.get<TH1>(HIST("hSelMask"))->SetBinContent(i + 1, applySelection[i]);
      }
    }

    // bc-based event selection qa
    for (auto& bc : bcs) {
      for (int iAlias = 0; iAlias < kNaliases; iAlias++) {
        histos.fill(HIST("hBcCounterAll"), iAlias, bc.alias()[iAlias]);
      }
    }

    // collision-based event selection qa
    for (auto& col : cols) {
      for (int iAlias = 0; iAlias < kNaliases; iAlias++) {
        if (!col.alias()[iAlias]) {
          continue;
        }
        histos.fill(HIST("hColCounterAll"), iAlias, 1);
        if (!col.sel7()) {
          continue;
        }
        histos.fill(HIST("hColCounterAcc"), iAlias, 1);
      }

      // further checks just on minimum bias triggers
      if (!isMC && !col.alias()[kINT7]) {
        continue;
      }
      for (int i = 0; i < kNsel; i++) {
        histos.fill(HIST("hSelCounter"), i, col.selection()[i]);
      }

      auto bc = col.bc_as<BCsRun2>();
      uint64_t globalBC = bc.globalBC();
      uint64_t orbit = globalBC / nBCsPerOrbit;
      int localBC = globalBC % nBCsPerOrbit;
      histos.fill(HIST("hGlobalBcAll"), globalBC);
      histos.fill(HIST("hOrbitAll"), orbit);
      histos.fill(HIST("hBcAll"), localBC);
      if (col.selection()[kIsBBV0A] || col.selection()[kIsBBV0C]) {
        histos.fill(HIST("hGlobalBcFV0"), globalBC);
        histos.fill(HIST("hOrbitFV0"), orbit);
        histos.fill(HIST("hBcFV0"), localBC);
      }
      if (col.selection()[kIsBBT0A] || col.selection()[kIsBBT0C]) {
        histos.fill(HIST("hGlobalBcFT0"), globalBC);
        histos.fill(HIST("hOrbitFT0"), orbit);
        histos.fill(HIST("hBcFT0"), localBC);
      }
      if (col.selection()[kIsBBFDA] || col.selection()[kIsBBFDC]) {
        histos.fill(HIST("hGlobalBcFDD"), globalBC);
        histos.fill(HIST("hOrbitFDD"), orbit);
        histos.fill(HIST("hBcFDD"), localBC);
      }

      float timeZNA = bc.has_zdc() ? bc.zdc().timeZNA() : -999.f;
      float timeZNC = bc.has_zdc() ? bc.zdc().timeZNC() : -999.f;
      float timeV0A = bc.has_fv0a() ? bc.fv0a().time() : -999.f;
      float timeV0C = bc.has_fv0c() ? bc.fv0c().time() : -999.f;
      float timeT0A = bc.has_ft0() ? bc.ft0().timeA() : -999.f;
      float timeT0C = bc.has_ft0() ? bc.ft0().timeC() : -999.f;
      float timeFDA = bc.has_fdd() ? bc.fdd().timeA() : -999.f;
      float timeFDC = bc.has_fdd() ? bc.fdd().timeC() : -999.f;
      float znSum = timeZNA + timeZNC;
      float znDif = timeZNA - timeZNC;
      float ofSPD = bc.spdFiredChipsL0() + bc.spdFiredChipsL1();
      float onSPD = bc.spdFiredFastOrL0() + bc.spdFiredFastOrL1();
      float multV0A = col.multRingV0A()[0] + col.multRingV0A()[1] + col.multRingV0A()[2] + col.multRingV0A()[3];
      float multV0C = col.multRingV0C()[0] + col.multRingV0C()[1] + col.multRingV0C()[2] + col.multRingV0C()[3];
      float multV0M = multV0A + multV0C;
      float multRingV0C3 = col.multRingV0C()[3];
      float multRingV0C012 = multV0C - multRingV0C3;
      float onV0M = bc.v0TriggerChargeA() + bc.v0TriggerChargeC();
      float ofV0M = multV0A + multV0C - col.multRingV0A()[0];
      int nTracklets = col.nTracklets();
      int spdClusters = col.spdClusters();

      float multT0A = 0;
      for (auto amplitude : bc.ft0().amplitudeA()) {
        multT0A += amplitude;
      }
      float multT0C = 0;
      for (auto amplitude : bc.ft0().amplitudeC()) {
        multT0C += amplitude;
      }
      float multFDA = 0;
      for (auto amplitude : bc.fdd().chargeA()) {
        multFDA += amplitude;
      }
      float multFDC = 0;
      for (auto amplitude : bc.fdd().chargeC()) {
        multFDC += amplitude;
      }
      float multZNA = bc.zdc().energyCommonZNA();
      float multZNC = bc.zdc().energyCommonZNC();

      histos.fill(HIST("hMultV0Mcol"), multV0M);
      histos.fill(HIST("hMultV0Acol"), multV0A);
      histos.fill(HIST("hMultV0Ccol"), multV0C);
      histos.fill(HIST("hMultZNAcol"), multZNA);
      histos.fill(HIST("hMultZNCcol"), multZNC);
      histos.fill(HIST("hMultT0Acol"), multT0A);
      histos.fill(HIST("hMultT0Ccol"), multT0C);
      histos.fill(HIST("hMultFDAcol"), multFDA);
      histos.fill(HIST("hMultFDCcol"), multFDC);

      histos.fill(HIST("hTimeV0Acol"), timeV0A);
      histos.fill(HIST("hTimeV0Ccol"), timeV0C);
      histos.fill(HIST("hTimeZNAcol"), timeZNA);
      histos.fill(HIST("hTimeZNCcol"), timeZNC);
      histos.fill(HIST("hTimeT0Acol"), timeT0A);
      histos.fill(HIST("hTimeT0Ccol"), timeT0C);
      histos.fill(HIST("hTimeFDAcol"), timeFDA);
      histos.fill(HIST("hTimeFDCcol"), timeFDC);
      histos.fill(HIST("hTimeZACcol"), znDif, znSum);
      histos.fill(HIST("hSPDClsVsTklCol"), nTracklets, spdClusters);
      histos.fill(HIST("hSPDOnVsOfCol"), ofSPD, onSPD);
      histos.fill(HIST("hV0MOnVsOfCol"), ofV0M, onV0M);
      histos.fill(HIST("hV0C3vs012Col"), multRingV0C012, multRingV0C3);
      histos.fill(HIST("hV0C012vsTklCol"), nTracklets, multRingV0C012);

      // filling plots for accepted events
      if (!col.sel7()) {
        continue;
      }

      histos.fill(HIST("hMultV0Macc"), multV0M);
      histos.fill(HIST("hMultV0Aacc"), multV0A);
      histos.fill(HIST("hMultV0Cacc"), multV0C);
      histos.fill(HIST("hMultZNAacc"), multZNA);
      histos.fill(HIST("hMultZNCacc"), multZNC);
      histos.fill(HIST("hMultT0Aacc"), multT0A);
      histos.fill(HIST("hMultT0Cacc"), multT0C);
      histos.fill(HIST("hMultFDAacc"), multFDA);
      histos.fill(HIST("hMultFDCacc"), multFDC);

      histos.fill(HIST("hTimeV0Aacc"), timeV0A);
      histos.fill(HIST("hTimeV0Cacc"), timeV0C);
      histos.fill(HIST("hTimeZNAacc"), timeZNA);
      histos.fill(HIST("hTimeZNCacc"), timeZNC);
      histos.fill(HIST("hTimeT0Aacc"), timeT0A);
      histos.fill(HIST("hTimeT0Cacc"), timeT0C);
      histos.fill(HIST("hTimeFDAacc"), timeFDA);
      histos.fill(HIST("hTimeFDCacc"), timeFDC);
      histos.fill(HIST("hTimeZACacc"), znDif, znSum);
      histos.fill(HIST("hSPDClsVsTklAcc"), nTracklets, spdClusters);
      histos.fill(HIST("hSPDOnVsOfAcc"), ofSPD, onSPD);
      histos.fill(HIST("hV0MOnVsOfAcc"), ofV0M, onV0M);
      histos.fill(HIST("hV0C3vs012Acc"), multRingV0C012, multRingV0C3);
      histos.fill(HIST("hV0C012vsTklAcc"), nTracklets, multRingV0C012);
    }
  }
  PROCESS_SWITCH(EventSelectionQaTask, processRun2, "Process Run2 event selection QA", true);

  void processRun3(
    soa::Join<aod::Collisions, aod::EvSels> cols,
    aod::FullTracks const& tracks,
    BCsRun3 const& bcs,
    aod::Zdcs const& zdcs,
    aod::FV0As const& fv0as,
    aod::FT0s const& ft0s,
    aod::FDDs const& fdds)
  {
    int runNumber = bcs.iteratorAt(0).runNumber();
    if (runNumber != lastRunNumber && runNumber >= 500000) { // using BC filling scheme for data only
      lastRunNumber = runNumber;
      auto bf = ccdb->getForTimeStamp<o2::BunchFilling>("GLO/GRP/BunchFilling", bcs.iteratorAt(0).timestamp());
      beamPatternA = bf->getBeamPattern(0);
      beamPatternC = bf->getBeamPattern(1);
      bcPatternA = beamPatternA & ~beamPatternC;
      bcPatternC = ~beamPatternA & beamPatternC;
      bcPatternB = beamPatternA & beamPatternC;
    }

    // background studies
    for (auto& bc : bcs) {
      int localBC = bc.globalBC() % nBCsPerOrbit;
      float timeV0A = bc.has_fv0a() ? bc.fv0a().time() : -999.f;
      float timeT0A = bc.has_ft0() ? bc.ft0().timeA() : -999.f;
      float timeT0C = bc.has_ft0() ? bc.ft0().timeC() : -999.f;
      float timeFDA = bc.has_fdd() ? bc.fdd().timeA() : -999.f;
      float timeFDC = bc.has_fdd() ? bc.fdd().timeC() : -999.f;
      if (bcPatternA[localBC] == 1 || bcPatternA[(localBC + 1) % nBCsPerOrbit]) {
        histos.fill(HIST("hTimeV0Abga"), timeV0A);
        histos.fill(HIST("hTimeT0Abga"), timeT0A);
        histos.fill(HIST("hTimeT0Cbga"), timeT0C);
        histos.fill(HIST("hTimeFDAbga"), timeFDA);
        histos.fill(HIST("hTimeFDCbga"), timeFDC);
      }
      if (bcPatternA[(localBC + 5) % nBCsPerOrbit]) {
        histos.fill(HIST("hTimeFDAbga"), timeFDA);
        histos.fill(HIST("hTimeFDCbga"), timeFDC);
      }
      if (bcPatternC[localBC] == 1 || bcPatternC[(localBC + 1) % nBCsPerOrbit]) {
        histos.fill(HIST("hTimeV0Abgc"), timeV0A);
        histos.fill(HIST("hTimeT0Abgc"), timeT0A);
        histos.fill(HIST("hTimeT0Cbgc"), timeT0C);
        histos.fill(HIST("hTimeFDAbgc"), timeFDA);
        histos.fill(HIST("hTimeFDCbgc"), timeFDC);
      }
      if (bcPatternC[(localBC + 5) % nBCsPerOrbit]) {
        histos.fill(HIST("hTimeFDAbgc"), timeFDA);
        histos.fill(HIST("hTimeFDCbgc"), timeFDC);
      }
    }

    // per-DF info to deduce FT0 rate
    if (bcs.size() > 0) {
      uint64_t orbit = bcs.iteratorAt(0).globalBC() / nBCsPerOrbit;
      histos.fill(HIST("hDFstartOrbit"), orbit);
      histos.fill(HIST("hFT0sPerDF"), orbit, ft0s.size());
    }

    // bc-based event selection qa
    for (auto& bc : bcs) {
      for (int iAlias = 0; iAlias < kNaliases; iAlias++) {
        histos.fill(HIST("hBcCounterAll"), iAlias, bc.alias()[iAlias]);
      }
      uint64_t globalBC = bc.globalBC();
      uint64_t orbit = globalBC / nBCsPerOrbit;
      int localBC = globalBC % nBCsPerOrbit;
      float timeZNA = bc.has_zdc() ? bc.zdc().timeZNA() : -999.f;
      float timeZNC = bc.has_zdc() ? bc.zdc().timeZNC() : -999.f;
      float timeV0A = bc.has_fv0a() ? bc.fv0a().time() : -999.f;
      float timeT0A = bc.has_ft0() ? bc.ft0().timeA() : -999.f;
      float timeT0C = bc.has_ft0() ? bc.ft0().timeC() : -999.f;
      float timeFDA = bc.has_fdd() ? bc.fdd().timeA() : -999.f;
      float timeFDC = bc.has_fdd() ? bc.fdd().timeC() : -999.f;
      histos.fill(HIST("hTimeV0Aall"), timeV0A);
      histos.fill(HIST("hTimeZNAall"), timeZNA);
      histos.fill(HIST("hTimeZNCall"), timeZNC);
      histos.fill(HIST("hTimeT0Aall"), timeT0A);
      histos.fill(HIST("hTimeT0Call"), timeT0C);
      histos.fill(HIST("hTimeFDAall"), timeFDA);
      histos.fill(HIST("hTimeFDCall"), timeFDC);
      if (bcPatternB[localBC]) {
        histos.fill(HIST("hTimeV0Aref"), timeV0A);
        histos.fill(HIST("hTimeZNAref"), timeZNA);
        histos.fill(HIST("hTimeZNCref"), timeZNC);
        histos.fill(HIST("hTimeT0Aref"), timeT0A);
        histos.fill(HIST("hTimeT0Cref"), timeT0C);
        histos.fill(HIST("hTimeFDAref"), timeFDA);
        histos.fill(HIST("hTimeFDCref"), timeFDC);
      }

      histos.fill(HIST("hGlobalBcAll"), globalBC);
      histos.fill(HIST("hOrbitAll"), orbit);
      histos.fill(HIST("hBcAll"), localBC);

      // FV0
      if (bc.has_fv0a()) {
        histos.fill(HIST("hGlobalBcFV0"), globalBC);
        histos.fill(HIST("hOrbitFV0"), orbit);
        histos.fill(HIST("hBcFV0"), localBC);
        float multV0A = 0;
        for (auto amplitude : bc.fv0a().amplitude()) {
          multV0A += amplitude;
        }
        histos.fill(HIST("hMultV0Aall"), multV0A);
        if (localBC == refBC) {
          histos.fill(HIST("hMultV0Aref"), multV0A);
        }
      }

      // FT0
      if (bc.has_ft0()) {
        histos.fill(HIST("hGlobalBcFT0"), globalBC);
        histos.fill(HIST("hOrbitFT0"), orbit);
        histos.fill(HIST("hBcFT0"), localBC);
        float multT0A = 0;
        for (auto amplitude : bc.ft0().amplitudeA()) {
          multT0A += amplitude;
        }
        float multT0C = 0;
        for (auto amplitude : bc.ft0().amplitudeC()) {
          multT0C += amplitude;
        }
        histos.fill(HIST("hMultT0Aall"), multT0A);
        histos.fill(HIST("hMultT0Call"), multT0C);
        if (localBC == refBC) {
          histos.fill(HIST("hMultT0Aref"), multT0A);
          histos.fill(HIST("hMultT0Cref"), multT0C);
        }
      }

      // FDD
      if (bc.has_fdd()) {
        histos.fill(HIST("hGlobalBcFDD"), globalBC);
        histos.fill(HIST("hOrbitFDD"), orbit);
        histos.fill(HIST("hBcFDD"), localBC);
        float multFDA = 0;
        for (auto amplitude : bc.fdd().chargeA()) {
          multFDA += amplitude;
        }
        float multFDC = 0;
        for (auto amplitude : bc.fdd().chargeC()) {
          multFDC += amplitude;
        }
        histos.fill(HIST("hMultFDAall"), multFDA);
        histos.fill(HIST("hMultFDCall"), multFDC);
        if (localBC == refBC) {
          histos.fill(HIST("hMultFDAref"), multFDA);
          histos.fill(HIST("hMultFDCref"), multFDC);
        }
      }

      // ZDC
      if (bc.has_zdc()) {
        float multZNA = bc.zdc().energyCommonZNA();
        float multZNC = bc.zdc().energyCommonZNC();
        histos.fill(HIST("hMultZNAall"), multZNA);
        histos.fill(HIST("hMultZNCall"), multZNC);
        if (localBC == refBC) {
          histos.fill(HIST("hMultZNAref"), multZNA);
          histos.fill(HIST("hMultZNCref"), multZNC);
        }
      }
    }

    // collision-based event selection qa
    for (auto& col : cols) {
      for (int iAlias = 0; iAlias < kNaliases; iAlias++) {
        if (!col.alias()[iAlias]) {
          continue;
        }
        histos.fill(HIST("hColCounterAll"), iAlias, 1);
        if (!col.sel8()) {
          continue;
        }
        histos.fill(HIST("hColCounterAcc"), iAlias, 1);
      }

      for (int i = 0; i < kNsel; i++) {
        histos.fill(HIST("hSelCounter"), i, col.selection()[i]);
      }

      auto bc = col.bc_as<BCsRun3>();
      uint64_t globalBC = bc.globalBC();
      uint64_t orbit = globalBC / nBCsPerOrbit;
      int localBC = globalBC % nBCsPerOrbit;
      histos.fill(HIST("hGlobalBcCol"), globalBC);
      histos.fill(HIST("hOrbitCol"), orbit);
      histos.fill(HIST("hBcCol"), localBC);

      auto tracksGrouped = tracks.sliceBy(aod::track::collisionId, col.globalIndex());
      int nTPCtracks = 0;
      int nTOFtracks = 0;
      for (auto& track : tracksGrouped) {
        if (!track.isPVContributor()) {
          continue;
        }
        nTPCtracks += track.hasTPC();
        nTOFtracks += track.hasTOF();
      }

      auto foundBC = col.foundBC_as<BCsRun3>();
      uint64_t foundGlobalBC = bc.globalBC();
      int nContributors = col.numContrib();
      float timeRes = col.collisionTimeRes();
      int bcDiff = int(globalBC - foundGlobalBC);
      histos.fill(HIST("hColTimeRes"), timeRes);
      histos.fill(HIST("hColBcDiffVsNcontrib"), nContributors, bcDiff);
      histos.fill(HIST("hColTimeResVsNcontrib"), nContributors, timeRes);
      if (nTPCtracks == 0) {
        histos.fill(HIST("hColBcDiffVsNcontribITSonly"), nContributors, bcDiff);
        histos.fill(HIST("hColTimeResVsNcontribITSonly"), nContributors, timeRes);
      }
      if (nTOFtracks > 0) {
        histos.fill(HIST("hColBcDiffVsNcontribWithTOF"), nContributors, bcDiff);
        histos.fill(HIST("hColTimeResVsNcontribWithTOF"), nContributors, timeRes);
      }
      histos.fill(HIST("hNcontribCol"), nContributors);

      float timeZNA = foundBC.has_zdc() ? foundBC.zdc().timeZNA() : -999.f;
      float timeZNC = foundBC.has_zdc() ? foundBC.zdc().timeZNC() : -999.f;
      float timeV0A = foundBC.has_fv0a() ? foundBC.fv0a().time() : -999.f;
      float timeT0A = foundBC.has_ft0() ? foundBC.ft0().timeA() : -999.f;
      float timeT0C = foundBC.has_ft0() ? foundBC.ft0().timeC() : -999.f;
      float timeFDA = foundBC.has_fdd() ? foundBC.fdd().timeA() : -999.f;
      float timeFDC = foundBC.has_fdd() ? foundBC.fdd().timeC() : -999.f;
      float znSum = timeZNA + timeZNC;
      float znDif = timeZNA - timeZNC;

      histos.fill(HIST("hTimeV0Acol"), timeV0A);
      histos.fill(HIST("hTimeZNAcol"), timeZNA);
      histos.fill(HIST("hTimeZNCcol"), timeZNC);
      histos.fill(HIST("hTimeT0Acol"), timeT0A);
      histos.fill(HIST("hTimeT0Ccol"), timeT0C);
      histos.fill(HIST("hTimeFDAcol"), timeFDA);
      histos.fill(HIST("hTimeFDCcol"), timeFDC);
      histos.fill(HIST("hTimeZACcol"), znDif, znSum);

      // FT0
      float multT0A = 0;
      float multT0C = 0;
      if (foundBC.has_ft0()) {
        for (auto amplitude : foundBC.ft0().amplitudeA()) {
          multT0A += amplitude;
        }
        for (auto amplitude : foundBC.ft0().amplitudeC()) {
          multT0C += amplitude;
        }
      }
      // FV0
      float multV0A = 0;
      if (foundBC.has_fv0a()) {
        for (auto amplitude : foundBC.fv0a().amplitude()) {
          multV0A += amplitude;
        }
      }
      // FDD
      float multFDA = 0;
      float multFDC = 0;
      if (foundBC.has_fdd()) {
        for (auto amplitude : foundBC.fdd().chargeA()) {
          multFDA += amplitude;
        }
        for (auto amplitude : foundBC.fdd().chargeC()) {
          multFDC += amplitude;
        }
      }
      histos.fill(HIST("hMultT0Acol"), multT0A);
      histos.fill(HIST("hMultT0Ccol"), multT0C);
      histos.fill(HIST("hMultV0Acol"), multV0A);
      histos.fill(HIST("hMultFDAcol"), multFDA);
      histos.fill(HIST("hMultFDCcol"), multFDC);

      // filling plots for accepted events
      if (!col.sel8()) {
        continue;
      }

      histos.fill(HIST("hTimeV0Aacc"), timeV0A);
      histos.fill(HIST("hTimeZNAacc"), timeZNA);
      histos.fill(HIST("hTimeZNCacc"), timeZNC);
      histos.fill(HIST("hTimeT0Aacc"), timeT0A);
      histos.fill(HIST("hTimeT0Cacc"), timeT0C);
      histos.fill(HIST("hTimeFDAacc"), timeFDA);
      histos.fill(HIST("hTimeFDCacc"), timeFDC);
      histos.fill(HIST("hTimeZACacc"), znDif, znSum);
      histos.fill(HIST("hMultT0Aacc"), multT0A);
      histos.fill(HIST("hMultT0Cacc"), multT0C);
      histos.fill(HIST("hMultV0Aacc"), multV0A);
      histos.fill(HIST("hMultFDAacc"), multFDA);
      histos.fill(HIST("hMultFDCacc"), multFDC);
      histos.fill(HIST("hNcontribAcc"), nContributors);
    }
  }
  PROCESS_SWITCH(EventSelectionQaTask, processRun3, "Process Run3 event selection QA", false);

  void processMCRun3(aod::McCollisions const& mcCols, BCsRun3 const& bcs)
  {
    for (auto& mcCol : mcCols) {
      auto bc = mcCol.bc_as<BCsRun3>();
      uint64_t globalBC = bc.globalBC();
      uint64_t orbit = globalBC / nBCsPerOrbit;
      int localBC = globalBC % nBCsPerOrbit;
      histos.fill(HIST("hGlobalBcColMC"), globalBC);
      histos.fill(HIST("hOrbitColMC"), orbit);
      histos.fill(HIST("hBcColMC"), localBC);
      histos.fill(HIST("hVertexXMC"), mcCol.posX());
      histos.fill(HIST("hVertexYMC"), mcCol.posY());
      histos.fill(HIST("hVertexZMC"), mcCol.posZ());
    }
  }
  PROCESS_SWITCH(EventSelectionQaTask, processMCRun3, "Process Run3 MC event selection QA", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<EventSelectionQaTask>(cfgc)};
}
