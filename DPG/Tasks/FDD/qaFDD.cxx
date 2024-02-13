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
#include "DataFormatsFDD/Digit.h"
#include "DataFormatsFIT/Triggers.h"
#include "Common/DataModel/FT0Corrected.h"

#include "CCDB/CcdbApi.h"
#include "CommonDataFormat/BunchFilling.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
const int nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;
using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;
using BCPattern = std::bitset<o2::constants::lhc::LHCMaxBunches>;
using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::BcSels, aod::Timestamps, aod::Run3MatchedToBCSparse>;

// const int gdeltaBC =5;

struct fddQA {
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  BCPattern CollidingBunch;
  double timeFrameInMs;
  int newRunNumber = -999;
  int oldRunNumber = -999;
  int nTF = 0;
  int nOrAFDD, nOrCFDD, nVertexFDD;
  int nOrAFTO, nOrCFTO, nVertexFTO;
  std::array<int, nBCsPerOrbit> RateVertexPerBCFDD, RateVertexPerBCFT0;

  Configurable<double> minOrbit{"minOrbit", 0, "minimum orbit"};
  Configurable<int> nOrbits{"nOrbits", 10000, "number of orbits"};
  Configurable<int> refBC{"refBC", 1238, "reference bc"};
  Configurable<int> nOrbitsPerTF{"nOrbitsPerTF", 256, "reference bc"};
  Configurable<bool> doZdcCorrela{"doZdcCorrela", true, "switch on the correlation plots for FIT and FDD"};

  OutputObj<TH2F> h2ChargeFT0CvsFV0A{
    TH2F("h2ChargeFT0CvsFV0A", "FT0 C  Vs FV0 A; FT0C-total charge (ADC); FV0A-total charge (ADC)", 800, 0.,
         8000., 10000, 0., 100000.)};
  OutputObj<TH1F> hBcCol{TH1F("hBcCol", ";;", nBCsPerOrbit, 0., static_cast<double>(nBCsPerOrbit))};
  OutputObj<TH1F> hBcFDD{TH1F("hBcFDD", ";;", nBCsPerOrbit, 0., static_cast<double>(nBCsPerOrbit))};
  OutputObj<TH2F> hBcOrbitColl{TH2F("hBcOrbitColl", "Orbit vs BC [Collision];Orbit;BC", static_cast<double>(nOrbitsPerTF),
                                    0, static_cast<double>(nOrbitsPerTF), nBCsPerOrbit, 0., static_cast<double>(nBCsPerOrbit))};
  OutputObj<TH2F> hBcOrbitFDD{TH2F("hBcOrbitFDD", "Orbit vs BC [FDD];Orbit;BC", static_cast<double>(nOrbitsPerTF),
                                   0, static_cast<double>(nOrbitsPerTF), nBCsPerOrbit, 0., static_cast<double>(nBCsPerOrbit))};
  OutputObj<TH1F> hTFDDA{TH1F("hTFDDA", "Time (FDDA); ns", 2000, -20, 20)};
  OutputObj<TH1F> hTFDDC{TH1F("hTFDDC", " Time (FDDC); ns", 2000, -20, 20)};
  OutputObj<TH1F> hTFDDAC{TH1F("hTFDDAC", " Time (FDDA+FDDC)/2; ns", 2000, -20, 20)};
  OutputObj<TH1F> hChFDDA{TH1F("hChFDDA", "FDDA; Charge in ADC;", 5010, -10, 5000)};
  OutputObj<TH1F> hChFDDC{TH1F("hChFDDC", "FDDC; Charge in ADC;", 5010, -10, 5000)};
  OutputObj<TH1F> hTotalChargeFDDAC{TH1F("hTotalChargeFDDAC", "FDDC; Charge in ADC;", 8010, -10, 8000)};
  OutputObj<TH1F> hNcontribColl{TH1F("hNcontribColl", "Ncontributers in Coll TABLE;#contributors", 100, -0.5, 99.5)};
  OutputObj<TH1F> hNcontribFDD{TH1F("hNcontribFDD", "Ncontributers in FDD;#contributors", 100, -0.5, 99.5)};
  OutputObj<TH1F> hNcontribFDDAC{
    TH1F("hNcontribFDDAC", "Ncontributers in FDD A and C;#contributors", 100, -0.5, 99.5)};
  OutputObj<TH1F> hNcontribFDDorAC{
    TH1F("hNcontribFDDorAC", "Ncontributers in FDD A or C;#contributors", 100, -0.5, 99.5)};
  OutputObj<TH1F> hNcontribFDDA{TH1F("hNcontribFDDA", "Ncontributers in FDDA;#contributors", 100, -0.5, 99.5)};
  OutputObj<TH1F> hNcontribFDDC{TH1F("hNcontribFDDC", "Ncontributers with FDDC;#contributors", 100, -0.5, 99.5)};

  HistogramRegistry registry;

  void init(InitContext&)
  {
    const AxisSpec axisMultT0M{1000, 0., 250000., "FT0M multiplicity"};
    const AxisSpec axisMultFDDM{1000, 0., 50000., "FDDM multiplicity"};
    const AxisSpec axisMultFDDA{1000, 0., 50000., "FDDA multiplicity"};
    const AxisSpec axisMultFDDC{1000, 0., 50000., "FDDC multiplicity"};
    const AxisSpec axisMultT0A{1000, 0., 250000., "FT0A multiplicity"};
    const AxisSpec axisMultT0C{1000, 0., 250000., "FT0C multiplicity"};
    const AxisSpec axisMultV0A{1000, 0., 250000., "FV0A multiplicity"};
    const AxisSpec axisBC{nBCsPerOrbit, 0., nBCsPerOrbit, "Bunch cross axis"};
    const AxisSpec axisSignalZN{400, -10., 4000, "ZN signal"};
    const AxisSpec axisSignalZNA{400, -10., 4000, "ZNA signal"};
    const AxisSpec axisSignalZNC{400, -10., 4000, "ZNC signal"};
    const AxisSpec axisColTime{3300, -16.5, 16.5};
    const AxisSpec axisVertex{1000, -50., 50.};
    const AxisSpec axisNfiredFT0{220, 0., 220.};
    const AxisSpec axisNfiredFV0{50, 0., 50.};

    const AxisSpec avgTimeNS{300, -15., 15, "Time in  ns"};
    const AxisSpec avgTimeSumAndDiff{300, -15., 15, "Time in  ns"};

    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    registry.add("FT0FV0/h2ChargeFT0CvsFV0A",
                 "charge correlation between FT0 and FV0;FT0C-total charge (ADC); FV0A-total charge (ADC)",
                 {HistType::kTH2F, {axisMultT0M, axisMultV0A}});
    registry.add("FT0FV0/h2ChargeFT0CvsFV0ACollidingBC",
                 "charge correlation between FT0 and FV0 (colliding BC);FT0C-total charge (ADC); FV0A-total charge (ADC)",
                 {HistType::kTH2F, {axisMultT0M, axisMultV0A}});
    registry.add("FT0FV0/bcFV0CH", "BC distribution (FV0CH);BC ID; Counts", {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FT0FV0/bcFV0CHandFT0VX", "BC distribution (FV0CH & FT0VX);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FT0FV0/bcFT0CE", "BC distribution (FT0CE);BC ID; Counts", {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FT0FV0/bcFT0VX", "BC distribution (FT0VX);BC ID; Counts", {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FT0FV0/bcFT0VXandFT0SC", "BC distribution (FT0VXandFT0SC);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FT0FV0/bcFT0VXandFT0SCorFT0CE", "BC distribution (FT0VXandFT0SCorFT0CE);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FT0FV0/h1Vertex", "FT0 vertex counter;; Counts", {HistType::kTH1F, {{1, 0, 1}}});
    registry.add("FT0FV0/bcFT0CEandFT0VX", "BC distribution (FT0CE & FT0VX);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FT0FV0Table/bcFT0VTX", "BC distribution (FT0VX);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FT0FV0Table/bcFT0CE", "BC distribution (FT0CE);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FT0FV0Table/bcFT0SC", "BC distribution (FT0SC);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});

    // hist for FIT correlation from BC table
    registry.add("FT0FV0FromBCTable/BC/bcFT0VTX", "BC distribution (FT0Vtx);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FT0FV0FromBCTable/BC/bcFT0Cen", "BC distribution (FT0Cen);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FT0FV0FromBCTable/BC/bcFT0SCen", "BC distribution (FT0SCen);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FT0FV0FromBCTable/BC/bcFT0OrA", "BC distribution (FT0OrA);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FT0FV0FromBCTable/BC/bcFT0OrC", "BC distribution (FT0OrC);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FT0FV0FromBCTable/BC/bcFT0VTXandFV0CH", "BC distribution (FT0VTX and FV0CH);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FT0FV0FromBCTable/BC/bcFT0CenandFT0VTX", "BC distribution (FT0VTX and FT0Cen);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});

    registry.add("FITFromBCTable/BC/bcFT0VtxandCenOrSCenCollBC", "BC distribution (FT0VTX and (FT0Cen or FT0SCen));BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FITFromBCTable/BC/bcFT0VtxandCenOrSCen", "BC distribution (FT0VTX and (FT0Cen or FT0SCen));BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FITFromBCTable/BC/bcFT0VtxandCenCollBC", "BC distribution (FT0VTX and FT0Cen);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FITFromBCTable/BC/bcFT0VtxandFV0CHCollBC", "BC distribution (FT0VTX and FV0SCH);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FITFromBCTable/BC/bcFT0VtxandCenOrSCenAndZDCCollBC", "BC distribution (FT0VTX and (FT0Cen or FT0SCen) and ZDC);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FITFromBCTable/BC/bcFT0VtxandCenOrSCenAndZDC", "BC distribution (FT0VTX and (FT0Cen or FT0SCen) and ZDC);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FITFromBCTable/BC/bcFT0VtxandFT0CenandZDCCollBC", "BC distribution (FT0VTX and FV0SCH);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FITFromBCTable/BC/bcFT0VtxandFV0CHandZDCCollBC", "BC distribution (FT0VTX and FV0SCH);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});

    // Total Charge FT0 for different triggers
    registry.add("FITFromBCTable/Charge1D/chargeFT0Vtx", "Amplitude distribution (FT0VTX);FT0M multiplicity; Counts",
                 {HistType::kTH1F, {axisMultT0M}});
    registry.add("FITFromBCTable/Charge1D/chargeFT0SCen", "Amplitude distribution (SCen);FT0M multiplicity; Counts",
                 {HistType::kTH1F, {axisMultT0M}});
    registry.add("FITFromBCTable/Charge1D/chargeFT0Cen", "Amplitude distribution (Cen);FT0M multiplicity; Counts",
                 {HistType::kTH1F, {axisMultT0M}});
    registry.add("FITFromBCTable/Charge1D/chargeFT0OrA", "Amplitude distribution (OrA);FT0M multiplicity; Counts",
                 {HistType::kTH1F, {axisMultT0M}});
    registry.add("FITFromBCTable/Charge1D/chargeFT0OrC", "Amplitude distribution (OrC);FT0M multiplicity; Counts",
                 {HistType::kTH1F, {axisMultT0M}});

    // Total Charge FV0 for different triggers
    registry.add("FITFromBCTable/Charge1D/chargeFV0OrA", "Amplitude distribution (OrA);FV0A multiplicity; Counts",
                 {HistType::kTH1F, {axisMultV0A}});
    registry.add("FITFromBCTable/Charge1D/chargeFV0Charge", "Amplitude distribution (NCharge);FV0A multiplicity; Counts",
                 {HistType::kTH1F, {axisMultV0A}});
    registry.add("FITFromBCTable/Charge1D/chargeFV0NCh", "Amplitude distribution (Nchan);FV0A multiplicity; Counts",
                 {HistType::kTH1F, {axisMultV0A}});
    registry.add("FITFromBCTable/Charge1D/chargeFV0Inner", "Amplitude distribution (Inner);FV0A multiplicity; Counts",
                 {HistType::kTH1F, {axisMultV0A}});
    registry.add("FITFromBCTable/Charge1D/chargeFV0Outer", "Amplitude distribution (Outer);FV0A multiplicity; Counts",
                 {HistType::kTH1F, {axisMultV0A}});

    registry.add("FITFromBCTable/Charge1D/chargeFV0OrAWithFT0vTX", "Amplitude distribution (OrA WithFT0vTX);FV0A multiplicity; Counts",
                 {HistType::kTH1F, {axisMultV0A}});

    // Number of fired channels for FT0 for different triggers
    registry.add("FITFromBCTable/FiredChan/FT0Vtx", "NFired Channels (FT0VTX);FT0-fired chanels; Counts",
                 {HistType::kTH1F, {axisNfiredFT0}});
    registry.add("FITFromBCTable/FiredChan/FT0SCen", "NFired Channels (SCen);FT0-fired chanels; Counts",
                 {HistType::kTH1F, {axisNfiredFT0}});
    registry.add("FITFromBCTable/FiredChan/FT0Cen", "NFired Channels (Cen);FT0-fired chanels; Counts",
                 {HistType::kTH1F, {axisNfiredFT0}});
    registry.add("FITFromBCTable/FiredChan/FT0OrA", "NFired Channels (OrA);FT0-fired chanels; Counts",
                 {HistType::kTH1F, {axisNfiredFT0}});
    registry.add("FITFromBCTable/FiredChan/FT0OrC", "NFired Channels (OrC);FT0-fired chanels; Counts",
                 {HistType::kTH1F, {axisNfiredFT0}});

    // Number of fired channels for FV0 for different triggers
    registry.add("FITFromBCTable/FiredChan/FV0Inner", "NFired Channels (FV0 Inner);FV0-fired chanels; Counts",
                 {HistType::kTH1F, {axisNfiredFV0}});
    registry.add("FITFromBCTable/FiredChan/FV0Outer", "NFired Channels (Outer);FV0-fired chanels; Counts",
                 {HistType::kTH1F, {axisNfiredFV0}});
    registry.add("FITFromBCTable/FiredChan/FV0NChan", "NFired Channels (NChan);FV0-fired chanels; Counts",
                 {HistType::kTH1F, {axisNfiredFV0}});
    registry.add("FITFromBCTable/FiredChan/FV0OrA", "NFired Channels (OrA);FV0-fired chanels; Counts",
                 {HistType::kTH1F, {axisNfiredFV0}});
    registry.add("FITFromBCTable/FiredChan/FV0Charge", "NFired Channels (Charge);FV0-fired chanels; Counts",
                 {HistType::kTH1F, {axisNfiredFV0}});

    // Fired Chanels Vs Total Charge
    registry.add("FITFromBCTable/chargeFT0MVsFiredChannels", "(chargeFT0M and nFired (FT0 Vtx));Charge; nFired",
                 {HistType::kTH2F, {axisMultT0M, axisNfiredFT0}});

    // FDD, FT0 time plot
    registry.add("FITFromBCTable/Time/hVertexVsCollTimeFT0Vtx", "Vertex vs. Coll Time; FT0 vertex (cm);Collision time (ns)", kTH2F, {axisVertex, axisColTime});
    registry.add("FITFromBCTable/Time/hVertexVsCollTimeFT0VtxCollBC", "Vertex vs. Coll Time (Colliding BC); FT0 vertex (cm);Collision time (ns)", kTH2F, {axisVertex, axisColTime});
    registry.add("FITFromBCTable/Time/hVertexVsCollTimeFDDVtx", "Vertex vs. Coll Time; FDD vertex (cm);Collision time (ns)", kTH2F, {axisVertex, axisColTime});
    registry.add("FITFromBCTable/Time/hVertexVsCollTimeFDDVtxCollBC", "Vertex vs. Coll Time (Colliding BC); FDD vertex (cm);Collision time (ns)", kTH2F, {axisVertex, axisColTime});

    // FT0 and FV0 charge
    registry.add("FT0FV0FromBCTable/chargeFT0MandFV0A", "(chargeFT0M and FV0A without trigger);BC ID; Counts",
                 {HistType::kTH2F, {axisMultV0A, axisMultT0M}});
    registry.add("FT0FV0FromBCTable/chargeFT0MVtxandFV0A", "(FT0M (Vtx) and FV0A (orA));BC ID; Counts",
                 {HistType::kTH2F, {axisMultV0A, axisMultT0M}});
    registry.add("FT0FV0FromBCTable/chargeFT0MVtxandFV0ACollBC",
                 "(FT0M (Vtx) and FV0A (orA) (Colliding BC));BC ID; Counts",
                 {HistType::kTH2F, {axisMultV0A, axisMultT0M}});

    registry.add("FT0FV0FromBCTable/chargeFT0AVtxandFV0A", "(FT0A (Vtx) and FV0A (orA));BC ID; Counts",
                 {HistType::kTH2F, {axisMultV0A, axisMultT0A}});
    registry.add("FT0FV0FromBCTable/chargeFT0CVtxandFV0A", "(FT0C (Vtx) and FV0A (orA));BC ID; Counts",
                 {HistType::kTH2F, {axisMultV0A, axisMultT0C}});

    registry.add("FT0FV0FromBCTable/chargeFT0AVtxandFV0ACollBC", "(FT0A (Vtx) and FV0A (orA): CollBC);BC ID; Counts",
                 {HistType::kTH2F, {axisMultV0A, axisMultT0A}});
    registry.add("FT0FV0FromBCTable/chargeFT0CVtxandFV0ACollBC", "(FT0C (Vtx) and FV0A (orA): CollBC);BC ID; Counts",
                 {HistType::kTH2F, {axisMultV0A, axisMultT0C}});

    // FT0 and FDD charge
    registry.add("FDDFT0FromBCTable/chargeFT0MVtxandFDDMVtx",
                 "(chargeFT0M (Vtx) and FDDM (Vtx) without trigger);BC ID; Counts",
                 {HistType::kTH2F, {axisMultFDDM, axisMultT0M}});
    registry.add("FDDFT0FromBCTable/chargeFT0AVtxandFDDAVtx", "(FT0A (Vtx) and FDDA (Vtx));BC ID; Counts",
                 {HistType::kTH2F, {axisMultFDDA, axisMultT0A}});
    registry.add("FDDFT0FromBCTable/chargeFT0CVtxandFDDCVtx",
                 "(FT0C (Vtx) and FDDC (Vtx) (Colliding BC));BC ID; Counts",
                 {HistType::kTH2F, {axisMultFDDC, axisMultT0C}});

    registry.add("FT0FV0FromBCTable/BC/bcFV0OrA", "BC distribution (FV0OrA);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});
    registry.add("FT0FV0FromBCTable/BC/bcFV0CH", "BC distribution (FV0Ch);BC ID; Counts",
                 {HistType::kTH1F, {{3564, 0, 3564}}});

    if (doZdcCorrela) {
      // Check if the process function for ZDCCollCorrela is enabled
      registry.add("FT0FV0FromBCTable/ZDCFIT/ZNvsFT0Mcorrel", "ZNvsFT0Mcorrel",
                   {HistType::kTH2F, {axisMultT0M, axisSignalZN}});
      registry.add("FT0FV0FromBCTable/ZDCFIT/ZNvsFT0MVTXcorrel", "ZNvsFT0MVTXcorrel",
                   {HistType::kTH2F, {axisMultT0M, axisSignalZN}});
      registry.add("FT0FV0FromBCTable/ZDCFIT/ZNvsFT0Mminbiascorrel", "ZNvsFT0MVTXcorrel (FT0 (cen||scen)&& Vtx)",
                   {HistType::kTH2F, {axisMultT0M, axisSignalZN}});
      registry.add("FT0FV0FromBCTable/ZDCFIT/ZNAvsFT0AVTXcorrel", "ZNAvsFT0AVTXcorrel",
                   {HistType::kTH2F, {axisMultT0A, axisSignalZNA}});
      registry.add("FT0FV0FromBCTable/ZDCFIT/ZNCvsFT0AVTXcorrel", "ZNAvsFT0AVTXcorrel",
                   {HistType::kTH2F, {axisMultT0A, axisSignalZNC}});
      registry.add("FT0FV0FromBCTable/ZDCFIT/ZNCvsFT0CVTXcorrel", "ZNCvsFT0CVTXcorrel",
                   {HistType::kTH2F, {axisMultT0C, axisSignalZNC}});
      registry.add("FT0FV0FromBCTable/ZDCFIT/ZNAvsFV0AWithFT0VTXcorrel", "ZNCvsFV0AWithFT0VTXcorrel",
                   {HistType::kTH2F, {axisMultV0A, axisSignalZNA}});
      registry.add("FT0FV0FromBCTable/ZDCFIT/ZNCvsFV0AWithFT0VTXcorrel", "ZNCvsFV0AWithFT0VTXcorrel",
                   {HistType::kTH2F, {axisMultV0A, axisSignalZNC}});
      registry.add("FT0FV0FromBCTable/ZDCFIT/ZNvsFV0AWithFT0VTXcorrel", "ZNCvsFV0AWithFT0VTXcorrel",
                   {HistType::kTH2F, {axisMultV0A, axisSignalZN}});

      registry.add("FDDFT0FromBCTable/BC/FDDVTX", "FDDVTX (BC distribution)", {HistType::kTH1F, {axisBC}});
      registry.add("FDDFT0FromBCTable/BC/FDDVTXCollBC", "FDDVTX (BC distribution-CollBC)", {HistType::kTH1F, {axisBC}});
      registry.add("FDDFT0FromBCTable/BC/FDDOrA", "FDDOrA (BC distribution)", {HistType::kTH1F, {axisBC}});
      registry.add("FDDFT0FromBCTable/BC/FDDOrC", "FDDOrC (BC distribution)", {HistType::kTH1F, {axisBC}});
      registry.add("FDDFT0FromBCTable/BC/FDDFT0VTX", "FDDFT0VTX (BC distribution with FT0 vtx)",
                   {HistType::kTH1F, {axisBC}});
      registry.add("FDDFT0FromBCTable/BC/FDDVTXandFT0VTX", "FDDVTXandFT0VTX (BC distribution)",
                   {HistType::kTH1F, {axisBC}});
    }
  }

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::FDDs const& fdds,
               aod::BCs const&)
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

  PROCESS_SWITCH(fddQA, process,
                 "Process FDD and FT0 info", true);

  void
    processCorr(soa::Join<aod::Collisions, aod::EvSels, aod::FT0sCorrected>::iterator const& col, aod::FT0s const& ft0s,
                aod::FV0As const& fv0s, aod::Zdcs const& zdcs, aod::BCs const&)
  {
    float sumAmpFT0C = 0;
    float sumAmpFV0 = 0;

    std::bitset<8> fT0Triggers;
    std::bitset<8> fV0Triggers;
    bool isVetexFT0 = kFALSE, isCentralFT0 = kFALSE, isSemiCentralFT0 = kFALSE;
    bool isFV0OrA = kFALSE;
    bool isFV0TrgNCh = kFALSE;

    // Trigger rates for FDD
    newRunNumber = col.bc().runNumber();
    // LOG(info)<< " newRunNumber  "<<newRunNumber;
    uint64_t ts{};

    uint64_t globalBC = col.bc().globalBC();
    // uint64_t orbit = globalBC % nOrbitsPerTF;
    int localBC = globalBC % nBCsPerOrbit;

    if (newRunNumber != oldRunNumber) {
      std::map<string, string> metadataRCT, headers;
      headers = ccdbApi.retrieveHeaders(Form("RCT/Info/RunInformation/%i", newRunNumber), metadataRCT, -1);
      ts = atol(headers["SOR"].c_str());

      LOG(info) << " newRunNumber  " << newRunNumber << " time stamp " << ts;
      oldRunNumber = newRunNumber;
      std::map<std::string, std::string> mapMetadata;
      std::map<std::string, std::string> mapHeader;
      auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", ts);
      CollidingBunch = grplhcif->getBunchFilling().getBCPattern();
      for (int i = 0; i < static_cast<int>(CollidingBunch.size()); i++) {
        if (CollidingBunch.test(i))
          LOG(info) << i << "  ";
      }
    } // new run number

    if (col.has_foundFV0()) {
      auto fv0 = col.foundFV0();
      fV0Triggers = fv0.triggerMask();
      isFV0OrA = fV0Triggers[o2::fit::Triggers::bitA];
      isFV0TrgNCh = fV0Triggers[o2::fit::Triggers::bitTrgCharge];

      for (std::size_t ich = 0; ich < fv0.amplitude().size(); ich++) {
        sumAmpFV0 += fv0.amplitude()[ich];
      }
      if (isFV0TrgNCh) {
        int localBCFV0 = fv0.bc().globalBC() % nBCsPerOrbit;
        registry.get<TH1>(HIST("FT0FV0/bcFV0CH"))->Fill(localBCFV0);
      }
    } // fv0

    if (col.has_foundFT0()) {
      auto ft0 = col.foundFT0();
      int localBCFT0 = ft0.bc().globalBC() % nBCsPerOrbit;
      fT0Triggers = ft0.triggerMask();
      isVetexFT0 = fT0Triggers[o2::fit::Triggers::bitVertex];
      isCentralFT0 = fT0Triggers[o2::fit::Triggers::bitCen];
      isSemiCentralFT0 = fT0Triggers[o2::fit::Triggers::bitSCen];

      if (col.t0CCorrectedValid()) {
        for (auto amplitude : ft0.amplitudeC()) {
          sumAmpFT0C += amplitude;
        }
      }

      bool isColl = CollidingBunch.test(localBC);
      if (col.has_foundFV0() && isFV0TrgNCh) {
        if (isVetexFT0) {
          registry.get<TH2>(HIST("FT0FV0/h2ChargeFT0CvsFV0A"))->Fill(sumAmpFT0C, sumAmpFV0);
          registry.get<TH1>(HIST("FT0FV0/bcFV0CHandFT0VX"))->Fill(localBCFT0);
        }
      } // fv0 and ft0
      if (isVetexFT0) {
        registry.get<TH1>(HIST("FT0FV0/bcFT0VX"))->Fill(localBCFT0);
        registry.get<TH1>(HIST("FT0FV0/h1Vertex"))->Fill(1);
        if (isCentralFT0)
          registry.get<TH1>(HIST("FT0FV0/bcFT0CEandFT0VX"))->Fill(localBCFT0);
        if (isSemiCentralFT0)
          registry.get<TH1>(HIST("FT0FV0/bcFT0VXandFT0SC"))->Fill(localBCFT0);
        if (isCentralFT0 || isSemiCentralFT0)
          registry.get<TH1>(HIST("FT0FV0/bcFT0VXandFT0SCorFT0CE"))->Fill(localBCFT0);
      }

      if (isCentralFT0) {
        registry.get<TH1>(HIST("FT0FV0/bcFT0CE"))->Fill(localBCFT0);
      }

      if (col.has_foundFV0() && isFV0OrA && isVetexFT0) {
        h2ChargeFT0CvsFV0A->Fill(sumAmpFT0C, sumAmpFV0);
        if (isColl)
          registry.get<TH2>(HIST("FT0FV0/h2ChargeFT0CvsFV0ACollidingBC"))->Fill(sumAmpFT0C, sumAmpFV0);
      }

      // correlation of FT0 with ZDC
      /*   if(col.has_foundZDC()) {
             h2ChargeFT0CvsFV0A->Fill(sumAmpFT0C, sumAmpFV0);
             if(isColl) registry.get<TH2>(HIST("ZDCFIT/ZNvsFT0correl"))->Fill(sumAmpFT0C, col.foundZDC().amplitudeZNA() +  col.foundZDC().amplitudeZNC());
         }
         */
    } // ft0
  }

  PROCESS_SWITCH(fddQA, processCorr,
                 "Process FDD time correlation", false);

  void processFT0(aod::FT0s const& ft0s, aod::FV0As const& fv0s, aod::BCs const&)
  {
    for (auto& ft0 : ft0s) {
      int localBCFT0 = ft0.bc().globalBC() % nBCsPerOrbit;
      std::bitset<8> triggerMASK = ft0.triggerMask();
      auto isVetexFT0 = triggerMASK[o2::fit::Triggers::bitVertex];
      auto isSCenFT0 = triggerMASK[o2::fit::Triggers::bitSCen];
      auto isCenFT0 = triggerMASK[o2::fit::Triggers::bitCen];
      if (isVetexFT0) {
        registry.get<TH1>(HIST("FT0FV0Table/bcFT0VTX"))->Fill(localBCFT0);
      }
      if (isCenFT0) {
        registry.get<TH1>(HIST("FT0FV0Table/bcFT0CE"))->Fill(localBCFT0);
      }
      if (isSCenFT0) {
        registry.get<TH1>(HIST("FT0FV0Table/bcFT0SC"))->Fill(localBCFT0);
      }
    }
  }

  PROCESS_SWITCH(fddQA, processFT0,
                 "Process FT0", false);

  void
    processFITFromBC(BCsWithRun3Matchings::iterator const& bc, aod::FV0As const&, aod::FT0s const&, aod::FDDs const&,
                     aod::Zdcs const& zdcs)
  {
    bool Ora = false;
    bool Orc = false;
    bool Tvx = false;
    bool Cent = false;
    bool SemiCentral = false;

    float multFT0C = 0.f;
    float multFT0A = 0.f;
    float multFV0A = 0.f;
    float multFT0M = 0.f;
    newRunNumber = bc.runNumber();

    if (newRunNumber != oldRunNumber) {
      uint64_t ts{};
      std::map<string, string> metadataRCT, headers;
      headers = ccdbApi.retrieveHeaders(Form("RCT/Info/RunInformation/%i", newRunNumber), metadataRCT, -1);
      ts = atol(headers["SOR"].c_str());

      LOG(info) << " newRunNumber  " << newRunNumber << " time stamp " << ts;
      oldRunNumber = newRunNumber;
      std::map<std::string, std::string> mapMetadata;
      std::map<std::string, std::string> mapHeader;
      auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", ts);
      CollidingBunch = grplhcif->getBunchFilling().getBCPattern();
      for (int i = 0; i < static_cast<int>(CollidingBunch.size()); i++) {
        if (CollidingBunch.test(i)) {
          LOG(info) << i << "  ";
        }
      }
    } // new run number

    if (bc.has_ft0()) {
      auto ft0 = bc.ft0();
      std::bitset<8> triggers = ft0.triggerMask();
      Ora = triggers[o2::fit::Triggers::bitA];
      Orc = triggers[o2::fit::Triggers::bitC];
      Tvx = triggers[o2::fit::Triggers::bitVertex];
      Cent = triggers[o2::fit::Triggers::bitCen];
      SemiCentral = triggers[o2::fit::Triggers::bitSCen];
      int localBCFT0 = bc.globalBC() % nBCsPerOrbit;
      // calculate charge
      for (auto amplitude : ft0.amplitudeA()) {
        multFT0A += amplitude;
      }
      for (auto amplitude : ft0.amplitudeC()) {
        multFT0C += amplitude;
      }
      multFT0M = multFT0A + multFT0C;

      int nfiredA_FT0 = 0;
      int nfiredC_FT0 = 0;
      int nFired_FT0 = 0;
      nfiredA_FT0 = ft0.channelA().size();
      nfiredC_FT0 = ft0.channelC().size();
      nFired_FT0 = nfiredA_FT0 + nfiredC_FT0;

      if (SemiCentral) {
        registry.get<TH1>(HIST("FITFromBCTable/Charge1D/chargeFT0SCen"))->Fill(multFT0M);
        registry.get<TH1>(HIST("FITFromBCTable/FiredChan/FT0SCen"))->Fill(nFired_FT0);
      }

      if (Cent) {
        registry.get<TH1>(HIST("FITFromBCTable/Charge1D/chargeFT0Cen"))->Fill(multFT0M);
        registry.get<TH1>(HIST("FITFromBCTable/FiredChan/FT0Cen"))->Fill(nFired_FT0);
      }

      if (Ora) {
        registry.get<TH1>(HIST("FITFromBCTable/Charge1D/chargeFT0OrA"))->Fill(multFT0M);
        registry.get<TH1>(HIST("FT0FV0FromBCTable/BC/bcFT0OrA"))->Fill(localBCFT0);
        registry.get<TH1>(HIST("FITFromBCTable/FiredChan/FT0OrA"))->Fill(nFired_FT0);
      }
      if (Orc) {
        registry.get<TH1>(HIST("FT0FV0FromBCTable/BC/bcFT0OrC"))->Fill(localBCFT0);
        registry.get<TH1>(HIST("FITFromBCTable/Charge1D/chargeFT0OrC"))->Fill(multFT0M);
        registry.get<TH1>(HIST("FITFromBCTable/FiredChan/FT0OrC"))->Fill(nFired_FT0);
      }
      if (Tvx) {
        registry.get<TH1>(HIST("FITFromBCTable/FiredChan/FT0Vtx"))->Fill(nFired_FT0);
        registry.get<TH1>(HIST("FITFromBCTable/Charge1D/chargeFT0Vtx"))->Fill(multFT0M);
        registry.get<TH1>(HIST("FT0FV0FromBCTable/BC/bcFT0VTX"))->Fill(localBCFT0);
        registry.get<TH2>(HIST("FITFromBCTable/Time/hVertexVsCollTimeFT0Vtx"))->Fill(((ft0.timeC() - ft0.timeA()) / 2) * o2::constants::physics::LightSpeedCm2NS, (ft0.timeA() + ft0.timeC()) / 2.0);
        if (CollidingBunch.test(localBCFT0)) {
          registry.get<TH2>(HIST("FITFromBCTable/Time/hVertexVsCollTimeFT0VtxCollBC"))->Fill(((ft0.timeC() - ft0.timeA()) / 2) * o2::constants::physics::LightSpeedCm2NS, (ft0.timeA() + ft0.timeC()) / 2.0);
        }
      }
      if (Cent) {
        registry.get<TH1>(HIST("FT0FV0FromBCTable/BC/bcFT0Cen"))->Fill(localBCFT0);
      }
      if (SemiCentral) {
        registry.get<TH1>(HIST("FT0FV0FromBCTable/BC/bcFT0SCen"))->Fill(localBCFT0);
      }
      if (Cent && Tvx) {
        registry.get<TH1>(HIST("FT0FV0FromBCTable/BC/bcFT0CenandFT0VTX"))->Fill(localBCFT0);
        if (CollidingBunch.test(localBCFT0)) {
          registry.get<TH1>(HIST("FITFromBCTable/BC/bcFT0VtxandCenCollBC"))->Fill(localBCFT0);
        }
      }
      if (Tvx && (SemiCentral || Cent)) {
        registry.get<TH1>(HIST("FITFromBCTable/BC/bcFT0VtxandCenOrSCen"))->Fill(localBCFT0);
        if (CollidingBunch.test(localBCFT0)) {
          registry.get<TH1>(HIST("FITFromBCTable/BC/bcFT0VtxandCenOrSCenCollBC"))->Fill(localBCFT0);
        }
      }
      if (bc.has_fv0a()) {
        auto fv0 = bc.fv0a();
        std::bitset<8> fV0Triggers = fv0.triggerMask();

        for (auto amplitude : fv0.amplitude()) {
          multFV0A += amplitude;
        }
        bool isFV0TrgNCh = fV0Triggers[o2::fit::Triggers::bitTrgCharge];
        bool isFV0OrA = fV0Triggers[o2::fit::Triggers::bitA];
        // here combine information from FT0
        registry.get<TH2>(HIST("FT0FV0FromBCTable/chargeFT0MandFV0A"))->Fill(multFV0A, multFT0M);
        if (isFV0TrgNCh && Tvx) {
          registry.get<TH1>(HIST("FT0FV0FromBCTable/BC/bcFT0VTXandFV0CH"))->Fill(localBCFT0);
          if (CollidingBunch.test(localBCFT0)) {
            if (bc.has_zdc()) {
              registry.get<TH1>(HIST("FITFromBCTable/BC/bcFT0VtxandFV0CHandZDCCollBC"))->Fill(localBCFT0);
            }
            registry.get<TH1>(HIST("FITFromBCTable/BC/bcFT0VtxandFV0CHCollBC"))->Fill(localBCFT0);
          }
        }
        if (Tvx && isFV0OrA) {
          registry.get<TH2>(HIST("FT0FV0FromBCTable/chargeFT0AVtxandFV0A"))->Fill(multFV0A, multFT0A);
          registry.get<TH2>(HIST("FT0FV0FromBCTable/chargeFT0CVtxandFV0A"))->Fill(multFV0A, multFT0C);
          registry.get<TH1>(HIST("FITFromBCTable/Charge1D/chargeFV0OrAWithFT0vTX"))->Fill(multFV0A);
          registry.get<TH2>(HIST("FT0FV0FromBCTable/chargeFT0MVtxandFV0A"))->Fill(multFV0A, multFT0M);
          if (CollidingBunch.test(localBCFT0)) {
            registry.get<TH2>(HIST("FT0FV0FromBCTable/chargeFT0MVtxandFV0ACollBC"))->Fill(multFV0A, multFT0M);
            registry.get<TH2>(HIST("FT0FV0FromBCTable/chargeFT0AVtxandFV0ACollBC"))->Fill(multFV0A, multFT0A);
            registry.get<TH2>(HIST("FT0FV0FromBCTable/chargeFT0CVtxandFV0ACollBC"))->Fill(multFV0A, multFT0C);
          }
        }
      } // fv0

      if (bc.has_fdd()) {
        float multFDDA = 0;
        float multFDDC = 0;
        float multFDDM = 0;

        auto fdd = bc.fdd();
        for (auto amplitude : fdd.chargeA()) {
          multFDDA += amplitude;
        }

        for (auto amplitude : fdd.chargeC()) {
          multFDDC += amplitude;
        }
        multFDDM = multFDDA + multFDDC;
        auto localBC = bc.globalBC() % nBCsPerOrbit;
        std::bitset<8> fddTriggers = fdd.triggerMask();
        bool FDDVtx = fddTriggers[o2::fit::Triggers::bitVertex];

        if (Tvx) {
          registry.get<TH1>(HIST("FDDFT0FromBCTable/BC/FDDFT0VTX"))->Fill(localBC);
          if (FDDVtx) {
            registry.get<TH2>(HIST("FDDFT0FromBCTable/chargeFT0MVtxandFDDMVtx"))->Fill(multFDDM, multFT0M);
            registry.get<TH2>(HIST("FDDFT0FromBCTable/chargeFT0AVtxandFDDAVtx"))->Fill(multFDDA, multFT0A);
            registry.get<TH2>(HIST("FDDFT0FromBCTable/chargeFT0CVtxandFDDCVtx"))->Fill(multFDDC, multFT0C);
            registry.get<TH1>(HIST("FDDFT0FromBCTable/BC/FDDVTXandFT0VTX"))->Fill(localBC);
          }
        }
      } // FDD

      if (bc.has_zdc()) {
        registry.get<TH2>(HIST("FT0FV0FromBCTable/ZDCFIT/ZNvsFT0Mcorrel"))->Fill(multFT0M, bc.zdc().amplitudeZNA() + bc.zdc().amplitudeZNC());
        if (Tvx && (SemiCentral || Cent)) {
          registry.get<TH2>(HIST("FT0FV0FromBCTable/ZDCFIT/ZNvsFT0Mminbiascorrel"))->Fill(multFT0M, bc.zdc().amplitudeZNA() + bc.zdc().amplitudeZNC());
          registry.get<TH1>(HIST("FITFromBCTable/BC/bcFT0VtxandCenOrSCenAndZDC"))->Fill(localBCFT0);
          if (CollidingBunch.test(localBCFT0)) {
            registry.get<TH1>(HIST("FITFromBCTable/BC/bcFT0VtxandCenOrSCenAndZDCCollBC"))->Fill(localBCFT0);
          }
        }
        if (Tvx && Cent && CollidingBunch.test(localBCFT0)) {
          registry.get<TH1>(HIST("FITFromBCTable/BC/bcFT0VtxandFT0CenandZDCCollBC"))->Fill(localBCFT0);
        }

        if (Tvx) {
          if (CollidingBunch.test(localBCFT0)) {
            registry.get<TH2>(HIST("FT0FV0FromBCTable/ZDCFIT/ZNvsFT0MVTXcorrel"))->Fill(multFT0M, bc.zdc().amplitudeZNA() + bc.zdc().amplitudeZNC());
            registry.get<TH2>(HIST("FT0FV0FromBCTable/ZDCFIT/ZNAvsFT0AVTXcorrel"))->Fill(multFT0A, bc.zdc().amplitudeZNA());
            registry.get<TH2>(HIST("FT0FV0FromBCTable/ZDCFIT/ZNCvsFT0AVTXcorrel"))->Fill(multFT0A, bc.zdc().amplitudeZNC());
            registry.get<TH2>(HIST("FT0FV0FromBCTable/ZDCFIT/ZNCvsFT0CVTXcorrel"))->Fill(multFT0C, bc.zdc().amplitudeZNC());
          }
        }
        if (bc.has_fv0a()) {
          if (Tvx) {
            registry.get<TH2>(HIST("FT0FV0FromBCTable/ZDCFIT/ZNvsFV0AWithFT0VTXcorrel"))->Fill(multFV0A, bc.zdc().amplitudeZNA() + bc.zdc().amplitudeZNC());
          }
          if (Tvx) {
            registry.get<TH2>(HIST("FT0FV0FromBCTable/ZDCFIT/ZNAvsFV0AWithFT0VTXcorrel"))->Fill(multFV0A, bc.zdc().amplitudeZNA());
          }
          if (Tvx) {
            registry.get<TH2>(HIST("FT0FV0FromBCTable/ZDCFIT/ZNCvsFV0AWithFT0VTXcorrel"))->Fill(multFV0A, bc.zdc().amplitudeZNC());
          }
        }
      }
    }
    if (bc.has_fdd()) {
      auto fdd = bc.fdd();
      auto localBC = bc.globalBC() % nBCsPerOrbit;
      std::bitset<8> fddTriggers = fdd.triggerMask();
      bool FDDVtx = fddTriggers[o2::fit::Triggers::bitVertex];
      bool FDDOrA = fddTriggers[o2::fit::Triggers::bitA];
      bool FDDOrC = fddTriggers[o2::fit::Triggers::bitC];

      if (FDDOrA) {
        registry.get<TH1>(HIST("FDDFT0FromBCTable/BC/FDDOrA"))->Fill(localBC);
      }
      if (FDDOrC) {
        registry.get<TH1>(HIST("FDDFT0FromBCTable/BC/FDDOrC"))->Fill(localBC);
      }

      if (FDDVtx) {
        registry.get<TH1>(HIST("FDDFT0FromBCTable/BC/FDDVTX"))->Fill(localBC);
        if (CollidingBunch.test(localBC)) {
          registry.get<TH1>(HIST("FDDFT0FromBCTable/BC/FDDVTXCollBC"))->Fill(localBC);
        }
        registry.get<TH2>(HIST("FITFromBCTable/Time/hVertexVsCollTimeFDDVtx"))->Fill(((fdd.timeC() - fdd.timeA()) / 2) * o2::constants::physics::LightSpeedCm2NS, (fdd.timeA() + fdd.timeC()) / 2.0);
        if (CollidingBunch.test(localBC)) {

          registry.get<TH2>(HIST("FITFromBCTable/Time/hVertexVsCollTimeFDDVtxCollBC"))->Fill(((fdd.timeC() - fdd.timeA()) / 2) * o2::constants::physics::LightSpeedCm2NS, (fdd.timeA() + fdd.timeC()) / 2.0);
        }
      }
    }

    if (bc.has_fv0a()) {
      auto fv0 = bc.fv0a();
      float multFV0ATrg = 0;
      for (auto amplitude : fv0.amplitude()) {
        multFV0ATrg += amplitude;
      }

      int NchannelFired = 0;
      NchannelFired = fv0.channel().size();

      std::bitset<8> fV0Triggers = fv0.triggerMask();
      bool isFV0OrA = fV0Triggers[o2::fit::Triggers::bitA];
      bool isFV0TrgNCh = fV0Triggers[o2::fit::Triggers::bitTrgCharge];
      bool isNChan = fV0Triggers[o2::fit::Triggers::bitTrgNchan];
      bool isFV0Inner = fV0Triggers[o2::fit::Triggers::bitAIn];
      bool isFV0Outer = fV0Triggers[o2::fit::Triggers::bitAOut];

      int localBCFV0 = bc.globalBC() % nBCsPerOrbit;
      if (isFV0OrA) {
        registry.get<TH1>(HIST("FT0FV0FromBCTable/BC/bcFV0OrA"))->Fill(localBCFV0);
        registry.get<TH1>(HIST("FITFromBCTable/Charge1D/chargeFV0OrA"))->Fill(multFV0ATrg);
        registry.get<TH1>(HIST("FITFromBCTable/FiredChan/FV0OrA"))->Fill(NchannelFired);
      }
      if (isFV0TrgNCh) {
        registry.get<TH1>(HIST("FT0FV0FromBCTable/BC/bcFV0CH"))->Fill(localBCFV0);
        registry.get<TH1>(HIST("FITFromBCTable/Charge1D/chargeFV0Charge"))->Fill(multFV0ATrg);
        registry.get<TH1>(HIST("FITFromBCTable/FiredChan/FV0Charge"))->Fill(NchannelFired);
      }
      if (isNChan) {
        registry.get<TH1>(HIST("FITFromBCTable/FiredChan/FV0NChan"))->Fill(NchannelFired);
        registry.get<TH1>(HIST("FITFromBCTable/Charge1D/chargeFV0NCh"))->Fill(multFV0ATrg);
      }
      if (isFV0Inner) {
        registry.get<TH1>(HIST("FITFromBCTable/Charge1D/chargeFV0Inner"))->Fill(multFV0ATrg);
        registry.get<TH1>(HIST("FITFromBCTable/FiredChan/FV0Inner"))->Fill(NchannelFired);
      }
      if (isFV0Outer) {
        registry.get<TH1>(HIST("FITFromBCTable/Charge1D/chargeFV0Outer"))->Fill(multFV0ATrg);
        registry.get<TH1>(HIST("FITFromBCTable/FiredChan/FV0Outer"))->Fill(NchannelFired);
      }
    }
  }
  PROCESS_SWITCH(fddQA, processFITFromBC, "Process FT0", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<fddQA>(cfgc, TaskName{"fdd-qa"})};
}
