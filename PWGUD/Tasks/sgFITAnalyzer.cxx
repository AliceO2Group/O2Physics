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
// \Single Gap Event Analyzer
// \author Sasha Bylinkin, alexander.bylinkin@gmail.com
// \since  April 2023

#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/SGTrackSelector.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Common/DataModel/PIDResponseTOF.h"

#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace bcNtr
{
DECLARE_SOA_COLUMN(BC, bc, uint64_t);
DECLARE_SOA_COLUMN(PV, pv, int);
DECLARE_SOA_COLUMN(GS, gs, int);
} // namespace bcNtr
namespace o2::aod
{
DECLARE_SOA_TABLE(BcPvGs, "AOD", "BCPVGS",
                  bcNtr::BC, bcNtr::PV, bcNtr::GS);
}
struct SGFITAnalyzer {
  Produces<o2::aod::BcPvGs> bcNtr;
  SGSelector sgSelector;
  Service<o2::framework::O2DatabasePDG> pdg;

  // configurables
  Configurable<bool> verbose{"Verbose", {}, "Additional print outs"};
  ConfigurableAxis ptAxis{"ptAxis", {250, 0.0, 2.5}, "p_T axis"};
  // ConfigurableAxis BCAxis{"BCAxis", {1000000000000, 0.5, 1000000000000.5}, "BC axis"};
  ConfigurableAxis BCAxis{"BCAxis", {100000000000, 500000000000.5, 600000000000.5}, "BC axis"};
  ConfigurableAxis etaAxis{"etaAxis", {300, -1.5, 1.5}, ""};
  ConfigurableAxis sigTPCAxis{"sigTPCAxis", {100, -100.0, 100.0}, ""};
  ConfigurableAxis sigTOFAxis{"sigTOFAxis", {100, -100.0, 100.0}, ""};
  ConfigurableAxis multAxis{"multAxis", {51, -.5, 50.5}, ""};
  ConfigurableAxis FitAxis{"FitAxis", {2000, -0.5, 3999.5}, ""};
  ConfigurableAxis ZDCAxis{"ZDCAxis", {1000, -2.5, 199.5}, ""};
  Configurable<float> FV0_cut{"FV0", 100., "FV0A threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};
  Configurable<float> FT0A_cut{"FT0A", 100., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 50., "FT0C threshold"};
  Configurable<float> FDDA_cut{"FDDA", 10000., "FDDA threshold"};
  Configurable<float> FDDC_cut{"FDDC", 10000., "FDDC threshold"};
  // Track Selections
  Configurable<float> PV_cut{"PV_cut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcaZ_cut{"dcaZ_cut", 2.0, "dcaZ cut"};
  Configurable<float> dcaXY_cut{"dcaXY_cut", 0.0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2_cut{"tpcChi2_cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindable_cut{"tpcNClsFindable_cut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2_cut{"itsChi2_cut", 36, "Max itsChi2NCl"};
  Configurable<float> eta_cut{"eta_cut", 0.9, "Track Pseudorapidity"};
  Configurable<float> pt_cut{"pt_cut", 0.1, "Track Pt"};
  Configurable<std::string> outputFileName{"outputFileName", "AnalysisResults.root", "Output file name"};
  // initialize histogram registry
  HistogramRegistry registry{
    "registry",
    {}};

  void init(InitContext&)
  {
    const AxisSpec axispt{ptAxis, "p_{T}"};
    const AxisSpec axismeanpt{ptAxis, "<p_{T}>"};
    const AxisSpec axisBC{BCAxis, "BC"};
    const AxisSpec axiseta{etaAxis, "#eta"};
    const AxisSpec axismult{multAxis, "N_{tracks}"};
    const AxisSpec axisfit{FitAxis, "FIT Amplitude"};
    const AxisSpec axiszdc{ZDCAxis, "ZDC Amplitude"};
    // Collision histograms
    registry.add("collisions/BC", "Relative BC number; Relative BC; Collisions", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
    registry.add("collisions/BC_PVCA", "Global BC; Global BC; Multiplicity", {HistType::kTH2F, {axisBC, axismult}});
    registry.add("collisions/BC_PVCC", "Global BC; Global BC; Multiplicity", {HistType::kTH2F, {axisBC, axismult}});
    registry.add("collisions/BC_PVCAC", "Global BC; Global BC; Multiplicity", {HistType::kTH2F, {axisBC, axismult}});
    registry.add("collisions/multiplicityAll", "Multiplicity of all tracks; Tracks; Tracks", {HistType::kTH1F, {{axismult}}});
    registry.add("collisions/multiplicityPVC", "Multiplicity of PV contributors; PV contributors; Tracks", {HistType::kTH1F, {{axismult}}});
    registry.add("collisions/multiplicityPVCA", "Multiplicity of PV contributors A-side; PV contributors; Tracks", {HistType::kTH1F, {{axismult}}});
    registry.add("collisions/multiplicityPVCC", "Multiplicity of PV contributors C-side; PV contributors; Tracks", {HistType::kTH1F, {{axismult}}});
    registry.add("collisions/multiplicityPVCAC", "Multiplicity of PV contributors AC-side; PV contributors; Tracks", {HistType::kTH1F, {{axismult}}});
    registry.add("collisions/multiplicityTPVC", "Multiplicity of PV contributors; PV contributors; Tracks", {HistType::kTH1F, {{axismult}}});
    registry.add("collisions/multiplicityTPVCA", "Multiplicity of PV contributors A-side; PV contributors; Tracks", {HistType::kTH1F, {{axismult}}});
    registry.add("collisions/multiplicityTPVCC", "Multiplicity of PV contributors C-side; PV contributors; Tracks", {HistType::kTH1F, {{axismult}}});
    registry.add("collisions/multiplicityTPVCAC", "Multiplicity of PV contributors AC-side; PV contributors; Tracks", {HistType::kTH1F, {{axismult}}});
    registry.add("collisions/multiplicityZ0PVCA", "Multiplicity of PV contributors A-side; PV contributors; Tracks", {HistType::kTH1F, {{axismult}}});
    registry.add("collisions/multiplicityZ0PVCC", "Multiplicity of PV contributors C-side; PV contributors; Tracks", {HistType::kTH1F, {{axismult}}});
    registry.add("collisions/multiplicityZ0PVCAC", "Multiplicity of PV contributors AC-side; PV contributors; Tracks", {HistType::kTH1F, {{axismult}}});
    registry.add("collisions/multiplicityZ1PVCA", "Multiplicity of PV contributors A-side; PV contributors; Tracks", {HistType::kTH1F, {{axismult}}});
    registry.add("collisions/multiplicityZ1PVCC", "Multiplicity of PV contributors C-side; PV contributors; Tracks", {HistType::kTH1F, {{axismult}}});
    registry.add("collisions/multiplicityZ1PVCAC", "Multiplicity of PV contributors AC-side; PV contributors; Tracks", {HistType::kTH1F, {{axismult}}});
    registry.add("collisions/GapSide", "Gap Side: A, C, A+C", {HistType::kTH1F, {{3, -0.5, 2.5}}});
    registry.add("collisions/TrueGapSide", "Gap Side: A, C, A+C", {HistType::kTH1F, {{4, -1.5, 2.5}}});
    registry.add("collisions/2D/multiplicityVsMeanPtPVCA", "Multiplicity of PV contributors A-side; PV contributors; mean #{p}_{T}", {HistType::kTH2F, {{axismult}, {axismeanpt}}});
    registry.add("collisions/2D/multiplicityVsEtaPVCA", "Multiplicity of PV contributors A-side; PV contributors; #eta", {HistType::kTH2F, {{axismult}, {axiseta}}});
    registry.add("collisions/2D/multiplicityVsMeanPtPVCC", "Multiplicity of PV contributors C-side; PV contributors; mean #{p}_{T}", {HistType::kTH2F, {{axismult}, {axismeanpt}}});
    registry.add("collisions/2D/multiplicityVsEtaPVCC", "Multiplicity of PV contributors C-side; PV contributors; #eta", {HistType::kTH2F, {{axismult}, {axiseta}}});
    registry.add("collisions/2D/multiplicityVsMeanPtPVCAC", "Multiplicity of PV contributors AC-side; PV contributors; mean #{p}_{T}", {HistType::kTH2F, {{axismult}, {axismeanpt}}});
    registry.add("collisions/2D/multiplicityVsEtaPVCAC", "Multiplicity of PV contributors AC-side; PV contributors; #eta", {HistType::kTH2F, {{axismult}, {axiseta}}});

    // track histograms
    registry.add("tracks/QCAll", "Track QC of all tracks; Hit in detector; Tracks", {HistType::kTH1F, {{5, -0.5, 4.5}}});
    registry.add("tracks/QCPVC", "Track QC of PV contributors; Hit in detector; Tracks", {HistType::kTH1F, {{5, -0.5, 4.5}}});
    registry.add("tracks/ptAll", "track pt of all tracks; p_{T} [GeV/c]; Tracks", {HistType::kTH1F, {axispt}});
    registry.add("tracks/ptPVC", "track pt of PV contributors; p_{T} [GeV/c]; Tracks", {HistType::kTH1F, {axispt}});
    registry.add("tracks/etaA", "track pt of PV contributors; p_{T} [GeV/c]; Tracks", {HistType::kTH1F, {axiseta}});
    registry.add("tracks/etaC", "track pt of PV contributors; p_{T} [GeV/c]; Tracks", {HistType::kTH1F, {axiseta}});
    registry.add("tracks/etaAC", "track pt of PV contributors; p_{T} [GeV/c]; Tracks", {HistType::kTH1F, {axiseta}});
    registry.add("tracks/etaApv", "track pt of PV contributors; p_{T} [GeV/c]; Tracks", {HistType::kTH1F, {axiseta}});
    registry.add("tracks/etaCpv", "track pt of PV contributors; p_{T} [GeV/c]; Tracks", {HistType::kTH1F, {axiseta}});
    registry.add("tracks/etaACpv", "track pt of PV contributors; p_{T} [GeV/c]; Tracks", {HistType::kTH1F, {axiseta}});
    registry.add("tracks/eta2Apv", "track pt of PV contributors; p_{T} [GeV/c]; Tracks", {HistType::kTH1F, {axiseta}});
    registry.add("tracks/eta2Cpv", "track pt of PV contributors; p_{T} [GeV/c]; Tracks", {HistType::kTH1F, {axiseta}});
    registry.add("tracks/eta2ACpv", "track pt of PV contributors; p_{T} [GeV/c]; Tracks", {HistType::kTH1F, {axiseta}});
    registry.add("tracks/etavsptA", "track eta versus pt of all tracks; eta; p_{T} [GeV/c]; Tracks", {HistType::kTH2F, {axiseta, axispt}});
    registry.add("tracks/etavsptC", "track eta versus pt of PV contributors; eta; p_{T} [GeV/c]; Tracks", {HistType::kTH2F, {axiseta, axispt}});
    registry.add("tracks/etavsptAC", "track eta versus pt of PV contributors; eta; p_{T} [GeV/c]; Tracks", {HistType::kTH2F, {axiseta, axispt}});

    const AxisSpec axisp{ptAxis, "momentum axis"};
    const AxisSpec axisTPCsig{sigTPCAxis, "TPC signal"};
    const AxisSpec axisTOFsig{sigTOFAxis, "TOF signal"};
    registry.add("tracks/TPCSignalvspAll", "TPC signal versus track momentum of all tracks; Track momentum [GeV/c]; TPC signal [arb. units]; Tracks", {HistType::kTH2F, {axisp, axisTPCsig}});
    registry.add("tracks/TPCSignalvspPVC", "TPC signal versus track momentum of PV contributors; Track momentum [GeV/c]; TPC signal [arb. units]; Tracks", {HistType::kTH2F, {axisp, axisTPCsig}});
    registry.add("tracks/TOFSignalvspAll", "TOF signal versus track momentum of all tracks; Track momentum [GeV/c]; TOF signal [arb. units]; Tracks", {HistType::kTH2F, {axisp, axisTOFsig}});
    registry.add("tracks/TOFSignalvspPVC", "TOF signal versus track momentum of PV contributors; Track momentum [GeV/c]; TOF signal [arb. units]; Tracks", {HistType::kTH2F, {axisp, axisTOFsig}});

    // FIT histograms
    registry.add("FIT/BBFV0A", "Beam-beam in V0A; BC relative to associated BC; Collisions", {HistType::kTH1F, {{32, -16.5, 15.5}}});
    registry.add("FIT/BBFT0A", "Beam-beam in T0A; BC relative to associated BC; Collisions", {HistType::kTH1F, {{32, -16.5, 15.5}}});
    registry.add("FIT/BBFT0C", "Beam-beam in T0C; BC relative to associated BC; Collisions", {HistType::kTH1F, {{32, -16.5, 15.5}}});
    registry.add("FIT/BBFDDA", "Beam-beam in FDDA; BC relative to associated BC; Collisions", {HistType::kTH1F, {{32, -16.5, 15.5}}});
    registry.add("FIT/BBFDDC", "Beam-beam in FDDA; BC relative to associated BC; Collisions", {HistType::kTH1F, {{32, -16.5, 15.5}}});
    registry.add("ZDC/AZNA", "Amplitude ZNA, A Gap", {HistType::kTH1F, {{axiszdc}}});
    registry.add("ZDC/CZNA", "Amplitude ZNA, C Gap", {HistType::kTH1F, {{axiszdc}}});
    registry.add("ZDC/ACZNA", "Amplitude ZNA, AC Gap", {HistType::kTH1F, {{axiszdc}}});
    registry.add("ZDC/ACZNA_CR", "Amplitude ZNA, AC Gap", {HistType::kTH1F, {{axiszdc}}});
    registry.add("ZDC/AZNA_CR", "Amplitude ZNA, AC Gap", {HistType::kTH1F, {{axiszdc}}});
    registry.add("ZDC/CZNA_CR", "Amplitude ZNA, AC Gap", {HistType::kTH1F, {{axiszdc}}});
    registry.add("ZDC/ACZNA_CM", "Amplitude ZNA, AC Gap", {HistType::kTH1F, {{axiszdc}}});
    registry.add("ZDC/AZNA_CM", "Amplitude ZNA, AC Gap", {HistType::kTH1F, {{axiszdc}}});
    registry.add("ZDC/CZNA_CM", "Amplitude ZNA, AC Gap", {HistType::kTH1F, {{axiszdc}}});
    registry.add("ZDC/AZNC", "Amplitude ZNC, A Gap", {HistType::kTH1F, {{axiszdc}}});
    registry.add("ZDC/CZNC", "Amplitude ZNC, C Gap", {HistType::kTH1F, {{axiszdc}}});
    registry.add("ZDC/ACZNC", "Amplitude ZNC, AC Gap", {HistType::kTH1F, {{axiszdc}}});
    registry.add("ZDC/ACZNC_CR", "Amplitude ZNC, AC Gap", {HistType::kTH1F, {{axiszdc}}});
    registry.add("ZDC/AZNC_CR", "Amplitude ZNC, AC Gap", {HistType::kTH1F, {{axiszdc}}});
    registry.add("ZDC/CZNC_CR", "Amplitude ZNC, AC Gap", {HistType::kTH1F, {{axiszdc}}});
    registry.add("ZDC/ACZNC_CM", "Amplitude ZNC, AC Gap", {HistType::kTH1F, {{axiszdc}}});
    registry.add("ZDC/AZNC_CM", "Amplitude ZNC, AC Gap", {HistType::kTH1F, {{axiszdc}}});
    registry.add("ZDC/CZNC_CM", "Amplitude ZNC, AC Gap", {HistType::kTH1F, {{axiszdc}}});
    registry.add("ZDC/tAZNA", "Time ZNA", {HistType::kTH1F, {{100, -19.5, 19.5}}});
    registry.add("ZDC/tAZNC", "Time ZNC", {HistType::kTH1F, {{100, -19.5, 19.5}}});
    registry.add("ZDC/tCZNA", "Time ZNA", {HistType::kTH1F, {{100, -19.5, 19.5}}});
    registry.add("ZDC/tCZNC", "Time ZNC", {HistType::kTH1F, {{100, -19.5, 19.5}}});
    registry.add("ZDC/tACZNA", "Time ZNA", {HistType::kTH1F, {{100, -19.5, 19.5}}});
    registry.add("ZDC/tACZNC", "Time ZNC", {HistType::kTH1F, {{100, -19.5, 19.5}}});
    registry.add("FIT/AFV0A", "Amplitude FV0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFV0A0", "Amplitude FV0A 0n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFV0A1", "Amplitude FV0A 1n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFV0A2", "Amplitude FV0A 2n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFV0A3", "Amplitude FV0A 3n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFV0A4", "Amplitude FV0A 4n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFT0A", "Amplitude FT0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFT0A0", "Amplitude FT0A 0n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFT0A1", "Amplitude FT0A 1n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFT0A2", "Amplitude FT0A 2n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFT0A3", "Amplitude FT0A 3n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFT0A4", "Amplitude FT0A 4n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFT0C", "Amplitude FT0C", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFT0C0", "Amplitude FT0C 0n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFT0C00", "Amplitude FT0C 0n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFT0C1", "Amplitude FT0C 1n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFT0C2", "Amplitude FT0C 2n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFT0C3", "Amplitude FT0C 3n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFT0C4", "Amplitude FT0C 4n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/tAFT0C", "Amplitude FT0C", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/zAFT0C", "Amplitude FT0C", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFDDA", "Amplitude FDDA", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFDDC", "Amplitude FDDC", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/zAFDDC", "Amplitude FDDC", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFDDC0", "Amplitude FDDC 0n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFDDC1", "Amplitude FDDC 1n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFDDC2", "Amplitude FDDC 2n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFDDC3", "Amplitude FDDC 3n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFDDC4", "Amplitude FDDC 4n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFV0A", "Amplitude FV0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/zCFV0A", "Amplitude FV0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFV0A0", "Amplitude FV0A 0n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFV0A1", "Amplitude FV0A 1n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFV0A2", "Amplitude FV0A 2n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFV0A3", "Amplitude FV0A 3n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFV0A4", "Amplitude FV0A 4n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFT0A", "Amplitude FT0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFT0A0", "Amplitude FT0A 0n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFT0A00", "Amplitude FT0A 0n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFT0A1", "Amplitude FT0A 1n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFT0A2", "Amplitude FT0A 2n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFT0A3", "Amplitude FT0A 3n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFT0A4", "Amplitude FT0A 4n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/tCFT0A", "Amplitude FT0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/zCFT0A", "Amplitude FT0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFT0C", "Amplitude FT0C", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFDDA", "Amplitude FDDA", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/zCFDDA", "Amplitude FDDA", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFDDA0", "Amplitude FDDA 0n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFDDA1", "Amplitude FDDA 1n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFDDA2", "Amplitude FDDA 2n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFDDA3", "Amplitude FDDA 3n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFDDA4", "Amplitude FDDA 4n", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFDDC", "Amplitude FDDC", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ACFV0A_CR", "Amplitude FV0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ACFT0A_CR", "Amplitude FT0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ACFT0C_CR", "Amplitude FT0C", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ACFDDA_CR", "Amplitude FDDA", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ACFDDC_CR", "Amplitude FDDC", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFV0A_CR", "Amplitude FV0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFT0A_CR", "Amplitude FT0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFT0C_CR", "Amplitude FT0C", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFDDA_CR", "Amplitude FDDA", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFDDC_CR", "Amplitude FDDC", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFV0A_CR", "Amplitude FV0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFT0A_CR", "Amplitude FT0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFT0C_CR", "Amplitude FT0C", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFDDA_CR", "Amplitude FDDA", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFDDC_CR", "Amplitude FDDC", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ACFV0A_CM", "Amplitude FV0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ACFT0A_CM", "Amplitude FT0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ACFT0C_CM", "Amplitude FT0C", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ACFDDA_CM", "Amplitude FDDA", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ACFDDC_CM", "Amplitude FDDC", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFV0A_CM", "Amplitude FV0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFT0A_CM", "Amplitude FT0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFT0C_CM", "Amplitude FT0C", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFDDA_CM", "Amplitude FDDA", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFDDC_CM", "Amplitude FDDC", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFV0A_CM", "Amplitude FV0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFT0A_CM", "Amplitude FT0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFT0C_CM", "Amplitude FT0C", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFDDA_CM", "Amplitude FDDA", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFDDC_CM", "Amplitude FDDC", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ACFV0A", "Amplitude FV0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ACFT0A", "Amplitude FT0A", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ACFT0C", "Amplitude FT0C", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ACFDDA", "Amplitude FDDA", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ACFDDC", "Amplitude FDDC", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ACFITC", "Amplitude FITC", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFITC", "Amplitude FITC", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFITC", "Amplitude FITC", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/ACFITA", "Amplitude FITA", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/AFITA", "Amplitude FITA", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/CFITA", "Amplitude FITA", {HistType::kTH1F, {{axisfit}}});
    registry.add("FIT/TFV0A", "Time FV0A", {HistType::kTH1F, {{100, -49.5, 50.5}}});
    registry.add("FIT/TFT0A", "Time FT0A", {HistType::kTH1F, {{100, -49.5, 50.5}}});
    registry.add("FIT/TFT0C", "Time FT0C", {HistType::kTH1F, {{100, -49.5, 50.5}}});
    registry.add("FIT/TFDDA", "Time FDDA", {HistType::kTH1F, {{100, -49.5, 50.5}}});
    registry.add("FIT/TFDDC", "Time FDDC", {HistType::kTH1F, {{100, -49.5, 50.5}}});
    registry.add("ZDC/AZNAC", "ZNA vs ZNC, A Gap", {HistType::kTH2F, {{axiszdc}, {axiszdc}}});
    registry.add("ZDC/CZNAC", "ZNA vs ZNC, C Gap", {HistType::kTH2F, {{axiszdc}, {axiszdc}}});
    registry.add("ZDC/ACZNAC", "ZNA vs ZNC, AC Gap", {HistType::kTH2F, {{axiszdc}, {axiszdc}}});
    registry.add("ZDC/AZNAFT0A", "ZNA vs FT0A, A Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/AZNAFV0A", "ZNA vs FV0A, A Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/AZNCFV0A", "ZNC vs FV0A, A Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/CZNAFT0A", "ZNA vs FT0A, C Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/CZNAFV0A", "ZNA vs FV0A, C Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/CZNCFV0A", "ZNC vs FV0A, C Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/ACZNAFT0A", "ZNA vs FT0A, AC Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/ACZNAFV0A", "ZNA vs FV0A, AC Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/ACZNCFV0A", "ZNC vs FV0A, AC Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/AZNAFT0C", "ZNA vs FT0C, A Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/CZNAFT0C", "ZNA vs FT0C, C Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/ACZNAFT0C", "ZNA vs FT0C, AC Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/AZNCFT0A", "ZNC vs FT0A, A Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/CZNCFT0A", "ZNC vs FT0A, C Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/ACZNCFT0A", "ZNC vs FT0A, AC Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/AZNCFT0C", "ZNC vs FT0C, A Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/CZNCFT0C", "ZNC vs FT0C, C Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/ACZNCFT0C", "ZNC vs FT0C, AC Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/AZNAFDDA", "ZNA vs FDDA, A Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/CZNAFDDA", "ZNA vs FDDA, C Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/ACZNAFDDA", "ZNA vs FDDA, AC Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/AZNAFDDC", "ZNA vs FDDC, A Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/CZNAFDDC", "ZNA vs FDDC, C Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/ACZNAFDDC", "ZNA vs FDDC, AC Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/AZNCFDDA", "ZNC vs FDDA, A Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/CZNCFDDA", "ZNC vs FDDA, C Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/ACZNCFDDA", "ZNC vs FDDA, AC Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/AZNCFDDC", "ZNC vs FDDC, A Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/CZNCFDDC", "ZNC vs FDDC, C Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("ZDC/ACZNCFDDC", "ZNC vs FDDC, AC Gap", {HistType::kTH2F, {{axiszdc}, {axisfit}}});
    registry.add("FIT/TAFV0A", "Time vs Amp FV0A", {HistType::kTH2F, {{40, -5.5, 34.5}, {axisfit}}});
    registry.add("FIT/TAFT0A", "Time vs Amp FT0A", {HistType::kTH2F, {{40, -5.5, 34.5}, {axisfit}}});
    registry.add("FIT/TTFT0", "Time FT0A vs Time FT0C", {HistType::kTH2F, {{40, -5.5, 34.5}, {40, -5.5, 34.5}}});
    registry.add("FIT/TTFT01", "Time FT0A vs Time FT0C", {HistType::kTH2F, {{40, -5.5, 34.5}, {40, -5.5, 34.5}}});
    registry.add("FIT/TTFT02", "Time FT0A vs Time FT0C", {HistType::kTH2F, {{40, -5.5, 34.5}, {40, -5.5, 34.5}}});
    registry.add("FIT/ACFT0", "Amp FITA vs Amp FITC", {HistType::kTH2F, {{100, -.5, 99.5}, {100, -.5, 99.5}}});
    registry.add("FIT/AFT0", "Amp FITA vs Amp FITC", {HistType::kTH2F, {{100, -.5, 99.5}, {100, -.5, 99.5}}});
    registry.add("FIT/CFT0", "Amp FITA vs Amp FITC", {HistType::kTH2F, {{100, -.5, 99.5}, {100, -.5, 99.5}}});
    registry.add("FIT/AFT0AC", "Amp FT0A vs Amp FT0C, A Gap", {HistType::kTH2F, {{axisfit}, {axisfit}}});
    registry.add("FIT/CFT0AC", "Amp FT0A vs Amp FT0C, C Gap", {HistType::kTH2F, {{axisfit}, {axisfit}}});
    registry.add("FIT/ACFT0AC", "Amp FT0A vs Amp FT0C, AC Gap", {HistType::kTH2F, {{axisfit}, {axisfit}}});
    registry.add("FIT/AFTV0A", "Amp FT0A vs Amp FV0A, A Gap", {HistType::kTH2F, {{axisfit}, {axisfit}}});
    registry.add("FIT/CFTV0A", "Amp FT0A vs Amp FV0A, C Gap", {HistType::kTH2F, {{axisfit}, {axisfit}}});
    registry.add("FIT/ACFTV0A", "Amp FT0A vs Amp FV0A, AC Gap", {HistType::kTH2F, {{axisfit}, {axisfit}}});
    registry.add("FIT/AFTV0C", "Amp FV0A vs Amp FT0C, A Gap", {HistType::kTH2F, {{axisfit}, {axisfit}}});
    registry.add("FIT/CFTV0C", "Amp FV0A vs Amp FT0C, C Gap", {HistType::kTH2F, {{axisfit}, {axisfit}}});
    registry.add("FIT/ACFTV0C", "Amp FV0A vs Amp FT0C, AC Gap", {HistType::kTH2F, {{axisfit}, {axisfit}}});
    registry.add("FIT/AFDV0A", "Amp FDDA vs Amp FV0A, A Gap", {HistType::kTH2F, {{axisfit}, {axisfit}}});
    registry.add("FIT/CFDV0A", "Amp FDDA vs Amp FV0A, C Gap", {HistType::kTH2F, {{axisfit}, {axisfit}}});
    registry.add("FIT/ACFDV0A", "Amp FDDA vs Amp FV0A, AC Gap", {HistType::kTH2F, {{axisfit}, {axisfit}}});
    registry.add("FIT/AFDT0A", "Amp FDDA vs Amp FT0A, A Gap", {HistType::kTH2F, {{axisfit}, {axisfit}}});
    registry.add("FIT/CFDT0A", "Amp FDDA vs Amp FT0A, C Gap", {HistType::kTH2F, {{axisfit}, {axisfit}}});
    registry.add("FIT/ACFDT0A", "Amp FDDA vs Amp FT0A, AC Gap", {HistType::kTH2F, {{axisfit}, {axisfit}}});
    registry.add("FIT/TAFT0C", "Time vs Amp FT0C", {HistType::kTH2F, {{40, -5.5, 34.5}, {axisfit}}});
    registry.add("FIT/TAFDDA", "Time vs Amp FDDA", {HistType::kTH2F, {{40, -5.5, 34.5}, {axisfit}}});
    registry.add("FIT/TAFDDC", "Time vs Amp FDDA", {HistType::kTH2F, {{40, -5.5, 34.5}, {axisfit}}});
    registry.add("FIT/TCFV0A", "Time vs Amp FV0A", {HistType::kTH2F, {{40, -5.5, 34.5}, {axisfit}}});
    registry.add("FIT/TCFT0A", "Time vs Amp FT0A", {HistType::kTH2F, {{40, -5.5, 34.5}, {axisfit}}});
    registry.add("FIT/TCFT0C", "Time vs Amp FT0C", {HistType::kTH2F, {{40, -5.5, 34.5}, {axisfit}}});
    registry.add("FIT/TCFDDA", "Time vs Amp FDDA", {HistType::kTH2F, {{40, -5.5, 34.5}, {axisfit}}});
    registry.add("FIT/TCFDDC", "Time vs Amp FDDC", {HistType::kTH2F, {{40, -5.5, 34.5}, {axisfit}}});
    registry.add("FIT/MFV0A", "Track number vs Amp FV0A", {HistType::kTH2F, {{axismult}, {axisfit}}});
    registry.add("FIT/MFT0A", "Track number vs Amp FT0A", {HistType::kTH2F, {{axismult}, {axisfit}}});
    registry.add("FIT/MFT0C", "Track number vs Amp FT0C", {HistType::kTH2F, {{axismult}, {axisfit}}});
    registry.add("FIT/MFDDA", "Track number vs Amp FDDA", {HistType::kTH2F, {{axismult}, {axisfit}}});
    registry.add("FIT/MFDDC", "Track number vs Amp FDDC", {HistType::kTH2F, {{axismult}, {axisfit}}});
    registry.add("ZDC/MAZNA", "Track number vs Amp FV0A", {HistType::kTH2F, {{axismult}, {axiszdc}}});
    registry.add("ZDC/MAZNC", "Track number vs Amp FT0A", {HistType::kTH2F, {{axismult}, {axiszdc}}});
    registry.add("ZDC/MCZNC", "Track number vs Amp FT0C", {HistType::kTH2F, {{axismult}, {axiszdc}}});
    registry.add("ZDC/MCZNA", "Track number vs Amp FDDA", {HistType::kTH2F, {{axismult}, {axiszdc}}});
    registry.add("ZDC/MACZNA", "Track number vs Amp FDDC", {HistType::kTH2F, {{axismult}, {axiszdc}}});
    registry.add("ZDC/MACZNC", "Track number vs Amp FDDC", {HistType::kTH2F, {{axismult}, {axiszdc}}});
  }

  // define data types
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>; // UDCollisions
  // using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions>; // UDCollisions
  using UDCollisionFull = UDCollisionsFull::iterator;
  // using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags>;
  using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  //  using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksFlags>;

  void process(UDCollisionFull const& dgcand, UDTracksFull const& dgtracks)
  {
    if (verbose) {
      LOGF(info, "<UDTutorial01sg> DG candidate %d", dgcand.globalIndex());
    }

    const float mpion = pdg->Mass(211);
    const float mmuon = pdg->Mass(13);
    // fill collision histograms
    registry.get<TH1>(HIST("collisions/GapSide"))->Fill(dgcand.gapSide(), 1.);
    // int truegapSide = sgSelector.trueGap(dgcand, FV0_cut, ZDC_cut);
    float FIT_cut[5] = {FV0_cut, FT0A_cut, FT0C_cut, FDDA_cut, FDDC_cut};
    // int truegapSide = sgSelector.trueGap(collision, *FIT_cut, ZDC_cut);
    int truegapSide = sgSelector.trueGap(dgcand, FIT_cut[0], FIT_cut[1], FIT_cut[2], ZDC_cut);
    int gs = truegapSide;
    registry.get<TH1>(HIST("collisions/TrueGapSide"))->Fill(truegapSide, 1.);
    // select PV contributors
    Partition<UDTracksFull> PVContributors = aod::udtrack::isPVContributor == true;
    PVContributors.bindTable(dgtracks);
    //    if (PVContributors.size() > 50)
    //     return;
    registry.get<TH1>(HIST("collisions/multiplicityPVC"))->Fill(PVContributors.size(), 1.);
    int pv = PVContributors.size();
    bool tof = false;
    // relative BC number
    auto bcnum = dgcand.globalBC() % o2::constants::lhc::LHCMaxBunches;
    uint64_t bc = dgcand.globalBC();
    registry.get<TH1>(HIST("collisions/BC"))->Fill(bcnum, 1.);
    // fill track histograms
    if (verbose) {
      LOGF(info, "<UDTutorial01sg>   Number of tracks %d", dgtracks.size());
      LOGF(info, "<UDTutorial01sg>   Number of PV contributors %d", PVContributors.size());
    }
    std::vector<float> parameters = {PV_cut, dcaZ_cut, dcaXY_cut, tpcChi2_cut, tpcNClsFindable_cut, itsChi2_cut, eta_cut, pt_cut};
    // check rho0 signals
    bool coh_rho0 = false;
    bool coh_jpsi = false;
    TLorentzVector rho, jpsi;
    std::vector<TLorentzVector> goodTracks;
    std::vector<TLorentzVector> muonTracks;
    float sign = 0;
    for (auto t : dgtracks) {
      TLorentzVector a;
      TLorentzVector b;
      a.SetXYZM(t.px(), t.py(), t.pz(), mpion);
      b.SetXYZM(t.px(), t.py(), t.pz(), mmuon);
      if (trackselector(t, parameters)) {
        sign += t.sign();
        goodTracks.push_back(a);
        if (std::abs(t.tpcNSigmaMu()) < 3)
          muonTracks.push_back(b);
      }
    }
    if (goodTracks.size() == 2) {
      for (auto pion : goodTracks) {
        rho += pion;
      }
      if (sign == 0 && TMath::Abs(rho.Rapidity()) < .9 && rho.M() > .5 && rho.M() < 1.2 && rho.Pt() < 0.1)
        coh_rho0 = true;
      if (muonTracks.size() == 2) {
        for (auto muon : muonTracks) {
          jpsi += muon;
        }
        if (sign == 0 && TMath::Abs(jpsi.Rapidity()) < .9 && jpsi.M() > 2.8 && jpsi.M() < 3.35 && jpsi.Pt() < 0.1)
          coh_jpsi = true;
      }
    }
    int pva = 0;
    double avPtPVa = 0;
    std::vector<double> vecEtaPVa;
    int pvc = 0;
    double avPtPVc = 0;
    std::vector<double> vecEtaPVc;
    int pvac = 0;
    double avPtPVac = 0;
    std::vector<double> vecEtaPVac;
    int z0pva = 0;
    int z0pvc = 0;
    int z0pvac = 0;
    int z1pva = 0;
    int z1pvc = 0;
    int z1pvac = 0;
    registry.get<TH1>(HIST("FIT/TFT0A"))->Fill(dgcand.timeFT0A(), 1.);
    registry.get<TH1>(HIST("FIT/TFT0C"))->Fill(dgcand.timeFT0C(), 1.);
    registry.get<TH1>(HIST("FIT/TFV0A"))->Fill(dgcand.timeFV0A(), 1.);
    registry.get<TH1>(HIST("FIT/TFDDA"))->Fill(dgcand.timeFDDA(), 1.);
    registry.get<TH1>(HIST("FIT/TFDDC"))->Fill(dgcand.timeFDDC(), 1.);
    //    if (truegapSide == 0) {
    registry.get<TH2>(HIST("FIT/TAFT0A"))->Fill(dgcand.timeFT0A(), dgcand.totalFT0AmplitudeA());
    if (truegapSide == 0)
      registry.get<TH2>(HIST("FIT/TTFT0"))->Fill(dgcand.timeFT0A(), dgcand.timeFT0C());
    if (truegapSide == 1)
      registry.get<TH2>(HIST("FIT/TTFT01"))->Fill(dgcand.timeFT0A(), dgcand.timeFT0C());
    if (truegapSide == 2)
      registry.get<TH2>(HIST("FIT/TTFT02"))->Fill(dgcand.timeFT0A(), dgcand.timeFT0C());
    if (truegapSide == 0)
      registry.get<TH2>(HIST("FIT/AFT0AC"))->Fill(dgcand.totalFT0AmplitudeA(), dgcand.totalFT0AmplitudeC());
    if (truegapSide == 1)
      registry.get<TH2>(HIST("FIT/CFT0AC"))->Fill(dgcand.totalFT0AmplitudeA(), dgcand.totalFT0AmplitudeC());
    if (truegapSide == 2)
      registry.get<TH2>(HIST("FIT/ACFT0AC"))->Fill(dgcand.totalFT0AmplitudeA(), dgcand.totalFT0AmplitudeC());
    if (truegapSide == 0)
      registry.get<TH2>(HIST("FIT/AFTV0A"))->Fill(dgcand.totalFT0AmplitudeA(), dgcand.totalFV0AmplitudeA());
    if (truegapSide == 1)
      registry.get<TH2>(HIST("FIT/CFTV0A"))->Fill(dgcand.totalFT0AmplitudeA(), dgcand.totalFV0AmplitudeA());
    if (truegapSide == 2)
      registry.get<TH2>(HIST("FIT/ACFTV0A"))->Fill(dgcand.totalFT0AmplitudeA(), dgcand.totalFV0AmplitudeA());
    if (truegapSide == 0)
      registry.get<TH2>(HIST("FIT/AFTV0C"))->Fill(dgcand.totalFV0AmplitudeA(), dgcand.totalFT0AmplitudeC());
    if (truegapSide == 1)
      registry.get<TH2>(HIST("FIT/CFTV0C"))->Fill(dgcand.totalFV0AmplitudeA(), dgcand.totalFT0AmplitudeC());
    if (truegapSide == 2)
      registry.get<TH2>(HIST("FIT/ACFTV0C"))->Fill(dgcand.totalFV0AmplitudeA(), dgcand.totalFT0AmplitudeC());
    if (truegapSide == 0)
      registry.get<TH2>(HIST("FIT/AFDV0A"))->Fill(dgcand.totalFDDAmplitudeA(), dgcand.totalFV0AmplitudeA());
    if (truegapSide == 1)
      registry.get<TH2>(HIST("FIT/CFDV0A"))->Fill(dgcand.totalFDDAmplitudeA(), dgcand.totalFV0AmplitudeA());
    if (truegapSide == 2)
      registry.get<TH2>(HIST("FIT/ACFDV0A"))->Fill(dgcand.totalFDDAmplitudeA(), dgcand.totalFV0AmplitudeA());
    if (truegapSide == 0)
      registry.get<TH2>(HIST("FIT/AFDT0A"))->Fill(dgcand.totalFDDAmplitudeA(), dgcand.totalFT0AmplitudeA());
    if (truegapSide == 1)
      registry.get<TH2>(HIST("FIT/CFDT0A"))->Fill(dgcand.totalFDDAmplitudeA(), dgcand.totalFT0AmplitudeA());
    if (truegapSide == 2)
      registry.get<TH2>(HIST("FIT/ACFDT0A"))->Fill(dgcand.totalFDDAmplitudeA(), dgcand.totalFT0AmplitudeA());
    float totalA = dgcand.totalFDDAmplitudeA() + dgcand.totalFT0AmplitudeA() + dgcand.totalFV0AmplitudeA();
    float totalC = dgcand.totalFDDAmplitudeC() + dgcand.totalFT0AmplitudeC();
    float zna = -1;
    float znc = -1;
    int an = 0;
    int cn = 0;
    if (dgcand.energyCommonZNC() > 0)
      znc = dgcand.energyCommonZNC();
    if (dgcand.energyCommonZNA() > 0)
      zna = dgcand.energyCommonZNA();
    if (zna > 0 && zna < 4)
      an = 1;
    else if (zna > 4 && zna < 10)
      an = 2;
    else if (zna > 10 && zna < 30)
      an = 3;
    else if (zna > 30)
      an = 4;
    if (znc > 0 && znc < 4)
      cn = 1;
    else if (znc > 4 && znc < 10)
      cn = 2;
    else if (znc > 10 && znc < 30)
      cn = 3;
    else if (znc > 30)
      cn = 4;
    // if (zna >0 && znc> 0)      LOGF(info, "ZNA  %f, ZNC %f", zna, znc);
    registry.get<TH2>(HIST("FIT/TAFT0C"))->Fill(dgcand.timeFT0C(), dgcand.totalFT0AmplitudeC());
    registry.get<TH2>(HIST("FIT/TAFV0A"))->Fill(dgcand.timeFV0A(), dgcand.totalFV0AmplitudeA());
    registry.get<TH2>(HIST("FIT/TAFDDA"))->Fill(dgcand.timeFDDA(), dgcand.totalFDDAmplitudeA());
    registry.get<TH2>(HIST("FIT/TAFDDC"))->Fill(dgcand.timeFDDC(), dgcand.totalFDDAmplitudeC());
    //    }
    if (truegapSide == 1) {
      registry.get<TH2>(HIST("FIT/TCFT0A"))->Fill(dgcand.timeFT0A(), dgcand.totalFT0AmplitudeA());
      registry.get<TH2>(HIST("FIT/TCFT0C"))->Fill(dgcand.timeFT0C(), dgcand.totalFT0AmplitudeC());
      registry.get<TH2>(HIST("FIT/TCFV0A"))->Fill(dgcand.timeFV0A(), dgcand.totalFV0AmplitudeA());
      registry.get<TH2>(HIST("FIT/TCFDDA"))->Fill(dgcand.timeFDDA(), dgcand.totalFDDAmplitudeA());
      registry.get<TH2>(HIST("FIT/TCFDDC"))->Fill(dgcand.timeFDDC(), dgcand.totalFDDAmplitudeC());
      registry.get<TH1>(HIST("ZDC/tCZNA"))->Fill(dgcand.timeZNA(), 1.);
      registry.get<TH1>(HIST("ZDC/tCZNC"))->Fill(dgcand.timeZNC(), 1.);
    }
    if (truegapSide == 0) {
      if (coh_rho0) {
        registry.get<TH1>(HIST("ZDC/AZNA_CR"))->Fill(zna, 1.);
        registry.get<TH1>(HIST("ZDC/AZNC_CR"))->Fill(znc, 1.);
        registry.get<TH1>(HIST("FIT/AFT0A_CR"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/AFT0C_CR"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
        registry.get<TH1>(HIST("FIT/AFV0A_CR"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/AFDDA_CR"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/AFDDC_CR"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
      }
      if (coh_jpsi) {
        registry.get<TH1>(HIST("ZDC/AZNA_CM"))->Fill(zna, 1.);
        registry.get<TH1>(HIST("ZDC/AZNC_CM"))->Fill(znc, 1.);
        registry.get<TH1>(HIST("FIT/AFT0A_CM"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/AFT0C_CM"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
        registry.get<TH1>(HIST("FIT/AFV0A_CM"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/AFDDA_CM"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/AFDDC_CM"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
      }
      registry.get<TH1>(HIST("ZDC/tAZNA"))->Fill(dgcand.timeZNA(), 1.);
      registry.get<TH1>(HIST("ZDC/tAZNC"))->Fill(dgcand.timeZNC(), 1.);
      registry.get<TH1>(HIST("ZDC/AZNA"))->Fill(zna, 1.);
      registry.get<TH1>(HIST("ZDC/AZNC"))->Fill(znc, 1.);
      registry.get<TH2>(HIST("ZDC/AZNAC"))->Fill(zna, znc);
      registry.get<TH2>(HIST("ZDC/AZNAFT0A"))->Fill(zna, dgcand.totalFT0AmplitudeA());
      registry.get<TH2>(HIST("ZDC/AZNAFT0C"))->Fill(zna, dgcand.totalFT0AmplitudeC());
      registry.get<TH2>(HIST("ZDC/AZNCFT0A"))->Fill(znc, dgcand.totalFT0AmplitudeA());
      registry.get<TH2>(HIST("ZDC/AZNCFT0C"))->Fill(znc, dgcand.totalFT0AmplitudeC());
      registry.get<TH2>(HIST("ZDC/AZNAFDDA"))->Fill(zna, dgcand.totalFDDAmplitudeA());
      registry.get<TH2>(HIST("ZDC/AZNAFDDC"))->Fill(zna, dgcand.totalFDDAmplitudeC());
      registry.get<TH2>(HIST("ZDC/AZNCFDDA"))->Fill(znc, dgcand.totalFDDAmplitudeA());
      registry.get<TH2>(HIST("ZDC/AZNCFDDC"))->Fill(znc, dgcand.totalFDDAmplitudeC());
      registry.get<TH2>(HIST("ZDC/AZNAFV0A"))->Fill(zna, dgcand.totalFV0AmplitudeA());
      registry.get<TH2>(HIST("ZDC/AZNCFV0A"))->Fill(znc, dgcand.totalFV0AmplitudeA());
      registry.get<TH2>(HIST("FIT/AFT0"))->Fill(totalA, totalC);
      registry.get<TH1>(HIST("FIT/AFITA"))->Fill(totalA, 1.);
      registry.get<TH1>(HIST("FIT/AFITC"))->Fill(totalC, 1.);
      registry.get<TH1>(HIST("FIT/AFT0A"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
      registry.get<TH1>(HIST("FIT/AFT0C"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
      if (an && cn)
        registry.get<TH1>(HIST("FIT/zAFT0C"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
      if (cn == 0) {
        registry.get<TH1>(HIST("FIT/AFT0A0"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/AFV0A0"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/AFT0C0"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
        if (an == 0)
          registry.get<TH1>(HIST("FIT/AFT0C00"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
        registry.get<TH1>(HIST("FIT/AFDDC0"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
      } else if (cn == 1) {
        registry.get<TH1>(HIST("FIT/AFT0A1"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/AFV0A1"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/AFT0C1"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
        registry.get<TH1>(HIST("FIT/AFDDC1"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
      } else if (cn == 2) {
        registry.get<TH1>(HIST("FIT/AFT0A2"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/AFV0A2"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/AFT0C2"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
        registry.get<TH1>(HIST("FIT/AFDDC2"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
      } else if (cn == 3) {
        registry.get<TH1>(HIST("FIT/AFT0A3"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/AFV0A3"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/AFT0C3"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
        registry.get<TH1>(HIST("FIT/AFDDC3"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
      } else if (cn == 4) {
        registry.get<TH1>(HIST("FIT/AFT0A4"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/AFV0A4"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/AFT0C4"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
        registry.get<TH1>(HIST("FIT/AFDDC4"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
      }
      if (an && cn)
        registry.get<TH1>(HIST("FIT/zAFDDC"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);

      registry.get<TH1>(HIST("FIT/AFV0A"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
      registry.get<TH1>(HIST("FIT/AFDDA"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
      registry.get<TH1>(HIST("FIT/AFDDC"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
      registry.get<TH2>(HIST("FIT/MFT0C"))->Fill(PVContributors.size(), dgcand.totalFT0AmplitudeC());
      registry.get<TH2>(HIST("FIT/MFDDC"))->Fill(PVContributors.size(), dgcand.totalFDDAmplitudeC());
      registry.get<TH2>(HIST("ZDC/MAZNA"))->Fill(PVContributors.size(), zna);
      registry.get<TH2>(HIST("ZDC/MAZNC"))->Fill(PVContributors.size(), znc);
    }
    if (truegapSide == 1) {
      if (coh_rho0) {
        registry.get<TH1>(HIST("ZDC/CZNA_CR"))->Fill(zna, 1.);
        registry.get<TH1>(HIST("ZDC/CZNC_CR"))->Fill(znc, 1.);
        registry.get<TH1>(HIST("FIT/CFT0A_CR"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/CFT0C_CR"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
        registry.get<TH1>(HIST("FIT/CFV0A_CR"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/CFDDA_CR"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/CFDDC_CR"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
      }
      if (coh_jpsi) {
        registry.get<TH1>(HIST("ZDC/CZNA_CM"))->Fill(zna, 1.);
        registry.get<TH1>(HIST("ZDC/CZNC_CM"))->Fill(znc, 1.);
        registry.get<TH1>(HIST("FIT/CFT0A_CM"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/CFT0C_CM"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
        registry.get<TH1>(HIST("FIT/CFV0A_CM"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/CFDDA_CM"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/CFDDC_CM"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
      }
      registry.get<TH1>(HIST("ZDC/CZNA"))->Fill(zna, 1.);
      registry.get<TH1>(HIST("ZDC/CZNC"))->Fill(znc, 1.);
      registry.get<TH2>(HIST("ZDC/CZNAC"))->Fill(zna, znc);
      registry.get<TH2>(HIST("ZDC/CZNAFT0A"))->Fill(zna, dgcand.totalFT0AmplitudeA());
      registry.get<TH2>(HIST("ZDC/CZNAFT0C"))->Fill(zna, dgcand.totalFT0AmplitudeC());
      registry.get<TH2>(HIST("ZDC/CZNCFT0A"))->Fill(znc, dgcand.totalFT0AmplitudeA());
      registry.get<TH2>(HIST("ZDC/CZNCFT0C"))->Fill(znc, dgcand.totalFT0AmplitudeC());
      registry.get<TH2>(HIST("ZDC/CZNAFDDA"))->Fill(zna, dgcand.totalFDDAmplitudeA());
      registry.get<TH2>(HIST("ZDC/CZNAFDDC"))->Fill(zna, dgcand.totalFDDAmplitudeC());
      registry.get<TH2>(HIST("ZDC/CZNCFDDA"))->Fill(znc, dgcand.totalFDDAmplitudeA());
      registry.get<TH2>(HIST("ZDC/CZNCFDDC"))->Fill(znc, dgcand.totalFDDAmplitudeC());
      registry.get<TH2>(HIST("ZDC/CZNAFV0A"))->Fill(zna, dgcand.totalFV0AmplitudeA());
      registry.get<TH2>(HIST("ZDC/CZNCFV0A"))->Fill(znc, dgcand.totalFV0AmplitudeA());
      registry.get<TH2>(HIST("FIT/CFT0"))->Fill(totalA, totalC);
      registry.get<TH1>(HIST("FIT/CFITA"))->Fill(totalA, 1.);
      registry.get<TH1>(HIST("FIT/CFITC"))->Fill(totalC, 1.);
      registry.get<TH2>(HIST("FIT/MFT0A"))->Fill(PVContributors.size(), dgcand.totalFT0AmplitudeA());
      registry.get<TH2>(HIST("FIT/MFV0A"))->Fill(PVContributors.size(), dgcand.totalFV0AmplitudeA());
      registry.get<TH2>(HIST("FIT/MFDDA"))->Fill(PVContributors.size(), dgcand.totalFDDAmplitudeA());
      registry.get<TH2>(HIST("ZDC/MCZNA"))->Fill(PVContributors.size(), zna);
      registry.get<TH2>(HIST("ZDC/MCZNC"))->Fill(PVContributors.size(), znc);
      registry.get<TH1>(HIST("FIT/CFT0C"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
      registry.get<TH1>(HIST("FIT/CFV0A"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
      registry.get<TH1>(HIST("FIT/CFT0A"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
      registry.get<TH1>(HIST("FIT/CFDDA"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
      if (an && cn) {
        registry.get<TH1>(HIST("FIT/zCFV0A"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/zCFT0A"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/zCFDDA"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
      }
      if (an == 0) {
        registry.get<TH1>(HIST("FIT/CFT0A0"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
        if (cn == 0)
          registry.get<TH1>(HIST("FIT/CFT0A00"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/CFV0A0"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/CFDDA0"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
      } else if (an == 1) {
        registry.get<TH1>(HIST("FIT/CFT0A1"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/CFV0A1"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/CFDDA1"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
      } else if (an == 2) {
        registry.get<TH1>(HIST("FIT/CFT0A2"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/CFV0A2"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/CFDDA2"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
      } else if (an == 3) {
        registry.get<TH1>(HIST("FIT/CFT0A3"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/CFV0A3"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/CFDDA3"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
      } else if (an == 4) {
        registry.get<TH1>(HIST("FIT/CFT0A4"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/CFV0A4"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/CFDDA4"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
      }

      registry.get<TH1>(HIST("FIT/CFDDA"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
      registry.get<TH1>(HIST("FIT/CFDDC"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
    }
    if (truegapSide == 2) {
      registry.get<TH1>(HIST("ZDC/tACZNA"))->Fill(dgcand.timeZNA(), 1.);
      registry.get<TH1>(HIST("ZDC/tACZNC"))->Fill(dgcand.timeZNC(), 1.);
      registry.get<TH1>(HIST("ZDC/ACZNA"))->Fill(zna, 1.);
      registry.get<TH1>(HIST("ZDC/ACZNC"))->Fill(znc, 1.);
      registry.get<TH2>(HIST("ZDC/ACZNAC"))->Fill(zna, znc);
      registry.get<TH2>(HIST("ZDC/ACZNAFT0A"))->Fill(zna, dgcand.totalFT0AmplitudeA());
      registry.get<TH2>(HIST("ZDC/ACZNAFT0C"))->Fill(zna, dgcand.totalFT0AmplitudeC());
      registry.get<TH2>(HIST("ZDC/ACZNCFT0A"))->Fill(znc, dgcand.totalFT0AmplitudeA());
      registry.get<TH2>(HIST("ZDC/ACZNCFT0C"))->Fill(znc, dgcand.totalFT0AmplitudeC());
      registry.get<TH2>(HIST("ZDC/ACZNAFDDA"))->Fill(zna, dgcand.totalFDDAmplitudeA());
      registry.get<TH2>(HIST("ZDC/ACZNAFDDC"))->Fill(zna, dgcand.totalFDDAmplitudeC());
      registry.get<TH2>(HIST("ZDC/ACZNCFDDA"))->Fill(znc, dgcand.totalFDDAmplitudeA());
      registry.get<TH2>(HIST("ZDC/ACZNCFDDC"))->Fill(znc, dgcand.totalFDDAmplitudeC());
      registry.get<TH2>(HIST("ZDC/ACZNAFV0A"))->Fill(zna, dgcand.totalFV0AmplitudeA());
      registry.get<TH2>(HIST("ZDC/ACZNCFV0A"))->Fill(znc, dgcand.totalFV0AmplitudeA());
      registry.get<TH2>(HIST("FIT/ACFT0"))->Fill(totalA, totalC);
      registry.get<TH1>(HIST("FIT/ACFITA"))->Fill(totalA, 1.);
      registry.get<TH1>(HIST("FIT/ACFITC"))->Fill(totalC, 1.);
      if (coh_rho0) {
        registry.get<TH1>(HIST("ZDC/ACZNA_CR"))->Fill(zna, 1.);
        registry.get<TH1>(HIST("ZDC/ACZNC_CR"))->Fill(znc, 1.);
        registry.get<TH1>(HIST("FIT/ACFT0A_CR"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/ACFT0C_CR"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
        registry.get<TH1>(HIST("FIT/ACFV0A_CR"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/ACFDDA_CR"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/ACFDDC_CR"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
      }
      if (coh_jpsi) {
        registry.get<TH1>(HIST("ZDC/ACZNA_CM"))->Fill(zna, 1.);
        registry.get<TH1>(HIST("ZDC/ACZNC_CM"))->Fill(znc, 1.);
        registry.get<TH1>(HIST("FIT/ACFT0A_CM"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/ACFT0C_CM"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
        registry.get<TH1>(HIST("FIT/ACFV0A_CM"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/ACFDDA_CM"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
        registry.get<TH1>(HIST("FIT/ACFDDC_CM"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
      }
      registry.get<TH1>(HIST("FIT/ACFT0A"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
      registry.get<TH1>(HIST("FIT/ACFT0C"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
      registry.get<TH1>(HIST("FIT/ACFV0A"))->Fill(dgcand.totalFV0AmplitudeA(), 1.);
      registry.get<TH1>(HIST("FIT/ACFDDA"))->Fill(dgcand.totalFDDAmplitudeA(), 1.);
      registry.get<TH1>(HIST("FIT/ACFDDC"))->Fill(dgcand.totalFDDAmplitudeC(), 1.);
      registry.get<TH2>(HIST("ZDC/MACZNA"))->Fill(PVContributors.size(), zna);
      registry.get<TH2>(HIST("ZDC/MACZNC"))->Fill(PVContributors.size(), znc);
    }
    for (auto track : dgtracks) {
      registry.get<TH1>(HIST("tracks/QCAll"))->Fill(0., 1.);
      registry.get<TH1>(HIST("tracks/QCAll"))->Fill(1., track.hasITS() * 1.);
      registry.get<TH1>(HIST("tracks/QCAll"))->Fill(2., track.hasTPC() * 1.);
      registry.get<TH1>(HIST("tracks/QCAll"))->Fill(3., track.hasTRD() * 1.);
      registry.get<TH1>(HIST("tracks/QCAll"))->Fill(4., track.hasTOF() * 1.);
      if (track.tofChi2() > -10)
        tof = true;
      auto vtrk = TVector3(track.px(), track.py(), track.pz());
      registry.get<TH1>(HIST("tracks/ptAll"))->Fill(track.pt(), 1.);
      if (truegapSide == 0)
        registry.get<TH1>(HIST("tracks/etaA"))->Fill(vtrk.Eta(), 1.);
      if (truegapSide == 1)
        registry.get<TH1>(HIST("tracks/etaC"))->Fill(vtrk.Eta(), 1.);
      if (truegapSide == 2)
        registry.get<TH1>(HIST("tracks/etaAC"))->Fill(vtrk.Eta(), 1.);

      auto signalTPC = track.tpcSignal() * track.sign();
      registry.get<TH2>(HIST("tracks/TPCSignalvspAll"))->Fill(vtrk.Mag(), signalTPC, 1.);
      auto signalTOF = track.tofSignal() * track.sign() / 1.E3;
      registry.get<TH2>(HIST("tracks/TOFSignalvspAll"))->Fill(vtrk.Mag(), signalTOF, 1.);
      if (track.isPVContributor()) {
        if (truegapSide == 0)
          registry.get<TH2>(HIST("tracks/etavsptA"))->Fill(vtrk.Eta(), track.pt(), 1.);
        if (truegapSide == 1)
          registry.get<TH2>(HIST("tracks/etavsptC"))->Fill(vtrk.Eta(), track.pt(), 1.);
        if (truegapSide == 2)
          registry.get<TH2>(HIST("tracks/etavsptAC"))->Fill(vtrk.Eta(), track.pt(), 1.);
        registry.get<TH1>(HIST("tracks/QCPVC"))->Fill(0., 1.);
        registry.get<TH1>(HIST("tracks/QCPVC"))->Fill(1., track.hasITS() * 1.);
        registry.get<TH1>(HIST("tracks/QCPVC"))->Fill(2., track.hasTPC() * 1.);
        registry.get<TH1>(HIST("tracks/QCPVC"))->Fill(3., track.hasTRD() * 1.);
        registry.get<TH1>(HIST("tracks/QCPVC"))->Fill(4., track.hasTOF() * 1.);
        registry.get<TH1>(HIST("tracks/ptPVC"))->Fill(track.pt(), 1.);
        //        registry.get<TH2>(HIST("tracks/etavsptPVC"))->Fill(vtrk.Eta(), track.pt(), 1.);
        if (truegapSide == 0) {
          pva++;
          avPtPVa += track.pt();
          vecEtaPVa.push_back(vtrk.Eta());
          registry.get<TH1>(HIST("tracks/etaApv"))->Fill(vtrk.Eta(), 1.);
          if (!an)
            registry.get<TH1>(HIST("tracks/eta2Apv"))->Fill(vtrk.Eta(), 1.);
          if (!an)
            z0pva++;
          else if (an == 4)
            z1pva++;
        }
        if (truegapSide == 1) {
          pvc++;
          avPtPVc += track.pt();
          vecEtaPVc.push_back(vtrk.Eta());
          registry.get<TH1>(HIST("tracks/etaCpv"))->Fill(vtrk.Eta(), 1.);
          if (!cn)
            registry.get<TH1>(HIST("tracks/eta2Cpv"))->Fill(vtrk.Eta(), 1.);
          if (!cn)
            z0pvc++;
          else if (cn == 4)
            z1pvc++;
        }
        if (truegapSide == 2) {
          pvac++;
          avPtPVac += track.pt();
          vecEtaPVac.push_back(vtrk.Eta());
          registry.get<TH1>(HIST("tracks/etaACpv"))->Fill(vtrk.Eta(), 1.);
          if (!an && !cn)
            registry.get<TH1>(HIST("tracks/eta2ACpv"))->Fill(vtrk.Eta(), 1.);
          if (!an && !cn)
            z0pvac++;
          else if (an >= 3 && cn >= 3)
            z1pvac++;
        }
        registry.get<TH2>(HIST("tracks/TPCSignalvspPVC"))->Fill(vtrk.Mag(), signalTPC, 1.);
        registry.get<TH2>(HIST("tracks/TOFSignalvspPVC"))->Fill(vtrk.Mag(), signalTOF, 1.);
      }
    }
    if (pva) {
      registry.get<TH1>(HIST("collisions/multiplicityPVCA"))->Fill(pva, 1.);
      registry.get<TH2>(HIST("collisions/2D/multiplicityVsMeanPtPVCA"))->Fill(pva, avPtPVa / pva, 1.);
      for (auto& element : vecEtaPVa)
        registry.get<TH2>(HIST("collisions/2D/multiplicityVsEtaPVCA"))->Fill(pva, element, 1.);
    }
    if (pvc) {
      registry.get<TH1>(HIST("collisions/multiplicityPVCC"))->Fill(pvc, 1.);
      registry.get<TH2>(HIST("collisions/2D/multiplicityVsMeanPtPVCC"))->Fill(pvc, avPtPVc / pvc, 1.);
      for (auto& element : vecEtaPVc)
        registry.get<TH2>(HIST("collisions/2D/multiplicityVsEtaPVCC"))->Fill(pvc, element, 1.);
    }
    if (pvac) {
      registry.get<TH1>(HIST("collisions/multiplicityPVCAC"))->Fill(pvac, 1.);
      registry.get<TH2>(HIST("collisions/2D/multiplicityVsMeanPtPVCAC"))->Fill(pvac, avPtPVac / pvac, 1.);
      for (auto& element : vecEtaPVac)
        registry.get<TH2>(HIST("collisions/2D/multiplicityVsEtaPVCAC"))->Fill(pvac, element, 1.);
    }
    if (pva)
      registry.get<TH1>(HIST("collisions/multiplicityZ0PVCA"))->Fill(z0pva, 1.);
    if (pvc)
      registry.get<TH1>(HIST("collisions/multiplicityZ0PVCC"))->Fill(z0pvc, 1.);
    if (pvac)
      registry.get<TH1>(HIST("collisions/multiplicityZ0PVCAC"))->Fill(z0pvac, 1.);
    if (pva)
      registry.get<TH1>(HIST("collisions/multiplicityZ1PVCA"))->Fill(z1pva, 1.);
    if (pvc)
      registry.get<TH1>(HIST("collisions/multiplicityZ1PVCC"))->Fill(z1pvc, 1.);
    if (pvac)
      registry.get<TH1>(HIST("collisions/multiplicityZ1PVCAC"))->Fill(z1pvac, 1.);
    if (tof) {
      if (truegapSide == 0)
        registry.get<TH1>(HIST("FIT/tAFT0C"))->Fill(dgcand.totalFT0AmplitudeC(), 1.);
      if (truegapSide == 1)
        registry.get<TH1>(HIST("FIT/tCFT0A"))->Fill(dgcand.totalFT0AmplitudeA(), 1.);
      registry.get<TH1>(HIST("collisions/multiplicityTPVC"))->Fill(PVContributors.size(), 1.);
      if (pva)
        registry.get<TH1>(HIST("collisions/multiplicityTPVCA"))->Fill(pva, 1.);
      if (pvc)
        registry.get<TH1>(HIST("collisions/multiplicityTPVCC"))->Fill(pvc, 1.);
      if (pvac)
        registry.get<TH1>(HIST("collisions/multiplicityTPVCAC"))->Fill(pvac, 1.);
    }
    // fill FIT histograms
    for (auto bit = 0; bit < 33; bit++) {
      registry.get<TH1>(HIST("FIT/BBFV0A"))->Fill(bit - 16, TESTBIT(dgcand.bbFV0Apf(), bit));
      registry.get<TH1>(HIST("FIT/BBFT0A"))->Fill(bit - 16, TESTBIT(dgcand.bbFT0Apf(), bit));
      registry.get<TH1>(HIST("FIT/BBFT0C"))->Fill(bit - 16, TESTBIT(dgcand.bbFT0Cpf(), bit));
      registry.get<TH1>(HIST("FIT/BBFDDA"))->Fill(bit - 16, TESTBIT(dgcand.bbFDDApf(), bit));
      registry.get<TH1>(HIST("FIT/BBFDDC"))->Fill(bit - 16, TESTBIT(dgcand.bbFDDCpf(), bit));
    }
    // Fill Table here
    bcNtr(bc, pv, gs);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SGFITAnalyzer>(cfgc, TaskName{"sgfitanalyzer"}),
  };
}
