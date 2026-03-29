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

/// \file nchCumulantsId.cxx
/// \brief Event by Event conserved charges fluctuations
/// \author Pravata Panigrahi <pravata.panigrahi@cern.ch> :: Sadhana Dash(sadhana@phy.iitb.ac.in)

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include <TPDGCode.h>

#include <algorithm>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics; // for constants
using namespace std;

#define ID_BIT_PI 0 // Identificationi bits for PID checks
#define ID_BIT_KA 1
#define ID_BIT_PR 2
#define ID_BIT_EL 3
#define ID_BIT_DE 4

#define BITSET(mask, ithBit) ((mask) |= (1 << (ithBit)))  // avoid name bitset as std::bitset is already there
#define BITCHECK(mask, ithBit) ((mask) & (1 << (ithBit))) // bit check will return int value, not bool, use BITCHECK != 0 in Analysi

enum PidEnum {
  kCh = 0,
  kPi, // dont use kPion, kKaon, as these enumeration
  kKa, // are already defined in $ROOTSYS/root/include/TPDGCode.h
  kPr,
  kEl,
  kDe
};

enum ChargeEnum {
  kPos = 0,
  kNeg = 1
};

static constexpr std::string_view PidDire[] = {
  "Ch/",
  "Pi/",
  "Ka/",
  "Pr/",
  "El/",
  "De/"};

static constexpr std::string_view ChargeDire[] = {
  "Pos/",
  "Neg/"};

std::string getModifiedStr(const std::string& myString)
{
  size_t pos = myString.rfind('/');
  if (pos != std::string::npos) {
    std::string subString = myString.substr(0, pos); // remove "/" from end of the string
    return subString;
  } else {
    return myString;
  }
}

struct NchCumulantsId {

  HistogramRegistry hist{"hist", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry recoTracks{"recoTracks", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry genAnalysis{"genAnalysis", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry recoAnalysis{"recoAnalysis", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry purityAnalysis{"purityAnalysis", {}, OutputObjHandlingPolicy::AnalysisObject};

  // PDG data base
  Service<o2::framework::O2DatabasePDG> pdgDB;

  Configurable<float> cfgCutPosZ{"cfgCutPosZ", 10.0, "cut for vertex Z"};
  Configurable<float> cfgCutDcaXY{"cfgCutDcaXY", 0.12, "cut for dcaXY"};
  Configurable<float> cfgCutDcaZ{"cfgCutDcaZ", 0.3, "cut for dcaZ"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "cut for eta"};
  Configurable<float> cfgCutPtMax{"cfgCutPtMax", 3.0, "max cut for pT"};
  Configurable<float> cfgCutPtMin{"cfgCutPtMin", 0.15, "min cut for pT"};

  Configurable<bool> checkCollPosZMc{"checkCollPosZMc", false, "checkCollPosZMc"};
  Configurable<bool> flagUnusedVariableError{"flagUnusedVariableError", false, "flagUnusedVariableError"};

  // Configurables for particle Identification
  Configurable<bool> cfgId01CheckVetoCut{"cfgId01CheckVetoCut", false, "cfgId01CheckVetoCut"};
  Configurable<bool> cfgId02DoElRejection{"cfgId02DoElRejection", true, "cfgId02DoElRejection"};
  Configurable<bool> cfgId03DoDeRejection{"cfgId03DoDeRejection", false, "cfgId03DoDeRejection"};
  Configurable<bool> cfgId04DoPdependentId{"cfgId04DoPdependentId", true, "cfgId04DoPdependentId"};
  Configurable<bool> cfgId05DoTpcInnerParamId{"cfgId05DoTpcInnerParamId", false, "cfgId05DoTpcInnerParamId"};

  Configurable<float> cfgIdPi01ThrPforTOF{"cfgIdPi01ThrPforTOF", 0.7, "cfgIdPi01ThrPforTOF"};
  Configurable<int> cfgIdPi02IdCutTypeLowP{"cfgIdPi02IdCutTypeLowP", 0, "cfgIdPi02IdCutTypeLowP"};
  Configurable<float> cfgIdPi03NSigmaTPCLowP{"cfgIdPi03NSigmaTPCLowP", 2.0, "cfgIdPi03NSigmaTPCLowP"};
  Configurable<float> cfgIdPi04NSigmaTOFLowP{"cfgIdPi04NSigmaTOFLowP", 2.0, "cfgIdPi04NSigmaTOFLowP"};
  Configurable<float> cfgIdPi05NSigmaRadLowP{"cfgIdPi05NSigmaRadLowP", 4.0, "cfgIdPi05NSigmaRadLowP"};
  Configurable<int> cfgIdPi06IdCutTypeHighP{"cfgIdPi06IdCutTypeHighP", 0, "cfgIdPi06IdCutTypeHighP"};
  Configurable<float> cfgIdPi07NSigmaTPCHighP{"cfgIdPi07NSigmaTPCHighP", 2.0, "cfgIdPi07NSigmaTPCHighP"};
  Configurable<float> cfgIdPi08NSigmaTOFHighP{"cfgIdPi08NSigmaTOFHighP", 2.0, "cfgIdPi08NSigmaTOFHighP"};
  Configurable<float> cfgIdPi09NSigmaRadHighP{"cfgIdPi09NSigmaRadHighP", 4.0, "cfgIdPi09NSigmaRadHighP"};

  Configurable<float> cfgIdKa01ThrPforTOF{"cfgIdKa01ThrPforTOF", 0.8, "cfgIdKa01ThrPforTOF"};
  Configurable<int> cfgIdKa02IdCutTypeLowP{"cfgIdKa02IdCutTypeLowP", 0, "cfgIdKa02IdCutTypeLowP"};
  Configurable<float> cfgIdKa03NSigmaTPCLowP{"cfgIdKa03NSigmaTPCLowP", 2.0, "cfgIdKa03NSigmaTPCLowP"};
  Configurable<float> cfgIdKa04NSigmaTOFLowP{"cfgIdKa04NSigmaTOFLowP", 2.0, "cfgIdKa04NSigmaTOFLowP"};
  Configurable<float> cfgIdKa05NSigmaRadLowP{"cfgIdKa05NSigmaRadLowP", 4.0, "cfgIdKa05NSigmaRadLowP"};
  Configurable<int> cfgIdKa06IdCutTypeHighP{"cfgIdKa06IdCutTypeHighP", 0, "cfgIdKa06IdCutTypeHighP"};
  Configurable<float> cfgIdKa07NSigmaTPCHighP{"cfgIdKa07NSigmaTPCHighP", 2.0, "cfgIdKa07NSigmaTPCHighP"};
  Configurable<float> cfgIdKa08NSigmaTOFHighP{"cfgIdKa08NSigmaTOFHighP", 2.0, "cfgIdKa08NSigmaTOFHighP"};
  Configurable<float> cfgIdKa09NSigmaRadHighP{"cfgIdKa09NSigmaRadHighP", 4.0, "cfgIdKa09NSigmaRadHighP"};

  Configurable<float> cfgIdPr01ThrPforTOF{"cfgIdPr01ThrPforTOF", 0.8, "cfgIdPr01ThrPforTOF"};
  Configurable<int> cfgIdPr02IdCutTypeLowP{"cfgIdPr02IdCutTypeLowP", 0, "cfgIdPr02IdCutTypeLowP"};
  Configurable<float> cfgIdPr03NSigmaTPCLowP{"cfgIdPr03NSigmaTPCLowP", 2.0, "cfgIdPr03NSigmaTPCLowP"};
  Configurable<float> cfgIdPr04NSigmaTOFLowP{"cfgIdPr04NSigmaTOFLowP", 2.0, "cfgIdPr04NSigmaTOFLowP"};
  Configurable<float> cfgIdPr05NSigmaRadLowP{"cfgIdPr05NSigmaRadLowP", 4.0, "cfgIdPr05NSigmaRadLowP"};
  Configurable<int> cfgIdPr06IdCutTypeHighP{"cfgIdPr06IdCutTypeHighP", 0, "cfgIdPr06IdCutTypeHighP"};
  Configurable<float> cfgIdPr07NSigmaTPCHighP{"cfgIdPr07NSigmaTPCHighP", 2.0, "cfgIdPr07NSigmaTPCHighP"};
  Configurable<float> cfgIdPr08NSigmaTOFHighP{"cfgIdPr08NSigmaTOFHighP", 2.0, "cfgIdPr08NSigmaTOFHighP"};
  Configurable<float> cfgIdPr09NSigmaRadHighP{"cfgIdPr09NSigmaRadHighP", 4.0, "cfgIdPr09NSigmaRadHighP"};

  struct : ConfigurableGroup {
    Configurable<float> cfgVetoId01PiTPC{"cfgVetoId01PiTPC", 3.0, "cfgVetoId01PiTPC"};
    Configurable<float> cfgVetoId02PiTOF{"cfgVetoId02PiTOF", 3.0, "cfgVetoId02PiTOF"};
    Configurable<float> cfgVetoId03KaTPC{"cfgVetoId03KaTPC", 3.0, "cfgVetoId03KaTPC"};
    Configurable<float> cfgVetoId04KaTOF{"cfgVetoId04KaTOF", 3.0, "cfgVetoId04KaTOF"};
    Configurable<float> cfgVetoId05PrTPC{"cfgVetoId05PrTPC", 3.0, "cfgVetoId05PrTPC"};
    Configurable<float> cfgVetoId06PrTOF{"cfgVetoId06PrTOF", 3.0, "cfgVetoId06PrTOF"};
    Configurable<float> cfgVetoId07ElTPC{"cfgVetoId07ElTPC", 3.0, "cfgVetoId07ElTPC"};
    Configurable<float> cfgVetoId08ElTOF{"cfgVetoId08ElTOF", 3.0, "cfgVetoId08ElTOF"};
    Configurable<float> cfgVetoId09DeTPC{"cfgVetoId09DeTPC", 3.0, "cfgVetoId09DeTPC"};
    Configurable<float> cfgVetoId10DeTOF{"cfgVetoId10DeTOF", 3.0, "cfgVetoId10DeTOF"};
  } cfgVetoIdCut;

  struct : ConfigurableGroup {
    Configurable<std::string> cfgCCDB01URL{"cfgCCDB01URL", "http://ccdb-test.cern.ch:8080", "cfgCCDB01URL"};
    Configurable<std::string> cfgCCDB02Path{"cfgCCDB02Path", "Users/p/ppanigra/NchCumulants/correctionWeights/", "cfgCCDB02Path"};
    Configurable<int> cfgCCDB03SOR{"cfgCCDB03SOR", 1, "cfgCCDB03SOR"};
    Configurable<std::string> cfgCCDB04StrHistP{"cfgCCDB04StrHistP", "h01_p", "cfgCCDB04StrHistP"};
    Configurable<std::string> cfgCCDB05StrHistPt{"cfgCCDB05StrHistPt", "h02_pt", "cfgCCDB05StrHistPt"};
    Configurable<std::string> cfgCCDB06StrHistPtEta{"cfgCCDB06StrHistPtEta", "h20_pt_eta", "cfgCCDB06StrHistPtEta"};
  } cfgCCDB;

  struct : ConfigurableGroup {
    Configurable<bool> printDebugMessages{"printDebugMessages", false, "printDebugMessages"};
    Configurable<bool> resetHistograms{"resetHistograms", false, "resetHistograms"};
  } cfgDebug;

  // Efficieny containing Histograms.
  TH2F* hPtEtaForBinSearch = nullptr;
  std::vector<std::array<TH2F*, kNeg + 1>> hPtEtaForEffCorrection{kDe + 1, std::array<TH2F*, kNeg + 1>{}};

  void init(InitContext const&)
  {
    auto& mgr = o2::ccdb::BasicCCDBManager::instance();
    mgr.setURL(cfgCCDB.cfgCCDB01URL);
    mgr.setCaching(true);
    auto ccdbObj = mgr.getForTimeStamp<TList>(cfgCCDB.cfgCCDB02Path, cfgCCDB.cfgCCDB03SOR);
    if (!ccdbObj) {
      if (cfgDebug.printDebugMessages)
        LOG(info) << "DEBUG :: CCDB OBJECT NOT FOUND";
    } else {
      if (cfgDebug.printDebugMessages)
        LOG(info) << "DEBUG :: CCDB OBJECT FOUND";
    }

    ccdbObj->Print();

    hPtEtaForBinSearch = reinterpret_cast<TH2F*>(ccdbObj->FindObject("hPtEta"));
    if (cfgDebug.printDebugMessages)
      LOG(info) << "DEBUG :: Obj Name = " << hPtEtaForBinSearch->GetName() << " :: entries = " << hPtEtaForBinSearch->GetEntries();
    std::string name = "";
    for (int i = 0; i <= kDe; i++) {
      for (int j = 0; j < (kNeg + 1); j++) {
        name = "hPtEta" + getModifiedStr(static_cast<std::string>(PidDire[i])) + getModifiedStr(static_cast<std::string>(ChargeDire[j]));
        hPtEtaForEffCorrection[i][j] = reinterpret_cast<TH2F*>(ccdbObj->FindObject(name.c_str()));
        if (cfgDebug.printDebugMessages)
          LOG(info) << "DEBUG :: Obj Name = " << hPtEtaForEffCorrection[i][j]->GetName() << " :: entries = " << hPtEtaForBinSearch->GetEntries();
      }
    }

    mgr.setURL("http://alice-ccdb.cern.ch"); // RESET the URL otherwise the other process functions which contains ccdb lookups will fail

    // QA check axes
    const AxisSpec axisEvents{1, 0, 1, "Counts"};
    const AxisSpec axisEta{100, -1., +1., "#eta"};
    const AxisSpec axisRapidity{200, -5, 5, "Rapidity (y)"};
    const AxisSpec axisPt{100, 0., 5., "p_{T} (GeV/c)"};
    const AxisSpec axisP{100, 0., 5., "p (GeV/c)"};
    const AxisSpec axisTPCInnerParam{100, 0, 3, "P_innerParam_Gev"};
    const AxisSpec axisdEdx(100, 20, 500, {"#frac{dE}{dx}"});
    const AxisSpec axisVtxZ{80, -20., 20., "V_{Z} (cm)"};
    const AxisSpec axisDCAz{200, -3., 3., "DCA_{Z} (cm)"};
    const AxisSpec axisDCAxy{200, -3., 3., "DCA_{XY} (cm)"};
    const AxisSpec axisMultFT0(150, 0, 1500, "MultFT0");
    const AxisSpec axisCent(103, -1., 102., "FT0C(%)");
    const AxisSpec axisPhi(80, -1, 7, "phi");

    const AxisSpec axisTOFBeta = {40, -2.0, 2.0, "tofBeta"};
    const AxisSpec axisTPCSignal = {100, -1, 1000, "tpcSignal"};
    const AxisSpec axisTPCNSigma = {200, -10.0, 10.0, "n#sigma_{TPC}"};
    const AxisSpec axisTOFNSigma = {200, -10.0, 10.0, "n#sigma_{TOF}"};
    const AxisSpec axisTOFExpMom = {200, 0.0f, 10.0f, "#it{p}_{tofExpMom} (GeV/#it{c})"};

    const AxisSpec axisNch(100, -50, 50, "Net_charge_dN");
    const AxisSpec axisPosCh(101, -1, 100, "Pos_charge");
    const AxisSpec axisNegCh(101, -1, 100, "Neg_charge");
    const AxisSpec axisNt(201, -1, 200, "Mult_midRap_Nch");
    const AxisSpec axisPrCh(101, -1, 100, "Pr_charge");
    const AxisSpec axisAPrCh(101, -1, 100, "APr_charge");
    const AxisSpec axisKaCh(101, -1, 100, "Ka_charge");
    const AxisSpec axisAKaCh(101, -1, 100, "AKa_charge");
    const AxisSpec axisPiCh(101, -1, 100, "Pion_Positive");
    const AxisSpec axisAPiCh(101, -1, 100, "Pion_Negative");

    HistogramConfigSpec qnHist1({HistType::kTHnSparseD, {axisNch, axisPosCh, axisNegCh, axisPrCh, axisAPrCh, axisKaCh, axisAKaCh, axisNt, axisCent}});
    HistogramConfigSpec qnHist2({HistType::kTHnSparseD, {axisNch, axisPosCh, axisNegCh, axisPiCh, axisAPiCh, axisKaCh, axisAKaCh, axisNt, axisCent}});

    HistogramConfigSpec histPPt({HistType::kTH2F, {axisP, axisPt}});
    HistogramConfigSpec histPTpcInnerParam({HistType::kTH2F, {axisP, axisTPCInnerParam}});
    HistogramConfigSpec histPTpcSignal({HistType::kTH2F, {axisP, axisTPCSignal}});
    HistogramConfigSpec histTpcInnerParamTpcSignal({HistType::kTH2F, {axisTPCInnerParam, axisTPCSignal}});
    HistogramConfigSpec histPBeta({HistType::kTH2F, {axisP, axisTOFBeta}});
    HistogramConfigSpec histTpcInnerParamBeta({HistType::kTH2F, {axisTPCInnerParam, axisTOFBeta}});
    HistogramConfigSpec histPTpcNSigma({HistType::kTH2F, {axisP, axisTPCNSigma}});
    HistogramConfigSpec histPtTpcNSigma({HistType::kTH2F, {axisPt, axisTPCNSigma}});
    HistogramConfigSpec histTpcInnerParamTpcNSigma({HistType::kTH2F, {axisTPCInnerParam, axisTPCNSigma}});
    HistogramConfigSpec histTofExpMomTpcNSigma({HistType::kTH2F, {axisTOFExpMom, axisTPCNSigma}});
    HistogramConfigSpec histPTofNSigma({HistType::kTH2F, {axisP, axisTOFNSigma}});
    HistogramConfigSpec histPtTofNSigma({HistType::kTH2F, {axisPt, axisTOFNSigma}});
    HistogramConfigSpec histTpcInnerParamTofNSigma({HistType::kTH2F, {axisTPCInnerParam, axisTOFNSigma}});
    HistogramConfigSpec histTofExpMomTofNSigma({HistType::kTH2F, {axisTOFExpMom, axisTOFNSigma}});
    HistogramConfigSpec histTpcNSigmaTofNSigma({HistType::kTH2F, {axisTPCNSigma, axisTOFNSigma}});

    HistogramConfigSpec histPtMc({HistType::kTH1F, {axisPt}});

    // QA check histos

    hist.add("QA/events/preSel/h_VtxZ", "V_{Z}", kTH1D, {axisVtxZ});
    hist.add("QA/events/preSel/h_Counts", "Counts", kTH1D, {axisEvents});
    hist.add("QA/events/preSel/multFT0", "multFT0", kTH1F, {axisMultFT0});
    hist.add("QA/events/preSel/centFT0", "centFT0", kTH1F, {axisCent});
    hist.addClone("QA/events/preSel/", "QA/events/postSel/");
    hist.add("QA/events/postSel/net_charge", "net_charge", kTH1F, {axisNch});
    hist.add("QA/events/postSel/Nt_centFT", "Mid_rap_Mult_VS_Cent", kTH2D, {{axisCent}, {axisNt}});

    hist.add("QA/tracks/preSel/h_P", "p (Gev/c)", kTH1D, {axisP});
    hist.add("QA/tracks/preSel/h_P_InnerParameter", "p_InnerParameter (Gev/c)", kTH1D, {axisTPCInnerParam});
    hist.add("QA/tracks/preSel/h_Pt", "p_{T} (TPC & TPC+TOF)", kTH1D, {axisPt});
    hist.add("QA/tracks/preSel/h_Eta", "#eta ", kTH1D, {axisEta});
    hist.add("QA/tracks/preSel/h_phi", "#phi ", kTH1D, {axisPhi});
    hist.add("QA/tracks/preSel/h2_Pt_DcaZ", "DCA_{z}", kTH2D, {{axisPt}, {axisDCAz}});
    hist.add("QA/tracks/preSel/h2_Pt_DcaXY", "DCA_{xy}", kTH2D, {{axisPt}, {axisDCAxy}});
    hist.add("QA/tracks/preSel/h2_p_DcaZ", "DCA_{z}", kTH2D, {{axisP}, {axisDCAz}});
    hist.add("QA/tracks/preSel/h2_p_DcaXY", "DCA_{xy}", kTH2D, {{axisP}, {axisDCAxy}});
    hist.add("QA/tracks/preSel/dE_dx1", "dE/dx vs p", kTH2F, {axisP, axisdEdx});
    hist.add("QA/tracks/preSel/dE_dx2", "dE/dx vs innerparam", kTH2F, {axisTPCInnerParam, axisdEdx});
    hist.add("QA/tracks/preSel/p_pt", "p_vs_pT", kTH2F, {axisP, axisPt});
    hist.add("QA/tracks/preSel/p_pInnerParameter", "p_vs_innerparameter", kTH2F, {axisP, axisTPCInnerParam});

    // tofBeta
    hist.add("QA/tracks/preSel/p_beta", "p_beta", histPBeta);
    hist.add("QA/tracks/preSel/tpcInnerParam_beta", "tpcInnerParam_beta", histTpcInnerParamBeta);

    // Look at Pion
    hist.add("QA/tracks/preSel/Pi/NoId/p_tpcNSigma", "p_tpcNSigma", histPTpcNSigma);
    hist.add("QA/tracks/preSel/Pi/NoId/pt_tpcNSigma", "pt_tpcNSigma", histPtTpcNSigma);
    hist.add("QA/tracks/preSel/Pi/NoId/tpcInnerParam_tpcNSigma", "tpcInnerParam_tpcNSigma", histTpcInnerParamTpcNSigma);
    hist.add("QA/tracks/preSel/Pi/NoId/tofExpMom_tpcNSigma", "tofExpMom_tpcNSigma", histTofExpMomTpcNSigma);
    hist.add("QA/tracks/preSel/Pi/NoId/p_tofNSigma", "p_tofNSigma", histPTofNSigma);
    hist.add("QA/tracks/preSel/Pi/NoId/pt_tofNSigma", "pt_tofNSigma", histPtTofNSigma);
    hist.add("QA/tracks/preSel/Pi/NoId/tpcInnerParam_tofNSigma", "tpcInnerParam_tofNSigma", histTpcInnerParamTofNSigma);
    hist.add("QA/tracks/preSel/Pi/NoId/tofExpMom_tofNSigma", "tofExpMom_tofNSigma", histTofExpMomTofNSigma);
    hist.add("QA/tracks/preSel/Pi/NoId/tpcNSigma_tofNSigma", "tpcNSigma_tofNSigma", histTpcNSigmaTofNSigma);

    hist.addClone("QA/tracks/preSel/Pi/", "QA/tracks/preSel/Ka/");
    hist.addClone("QA/tracks/preSel/Pi/", "QA/tracks/preSel/Pr/");

    hist.addClone("QA/tracks/preSel/", "QA/tracks/postSel/");

    hist.add("QA/tracks/Idfd/Pi/tpcId/p_pt", "p_pt", histPPt);
    hist.add("QA/tracks/Idfd/Pi/tpcId/p_tpcInnerParam", "p_tpcInnerParam", histPTpcInnerParam);
    // tpcSignal
    hist.add("QA/tracks/Idfd/Pi/tpcId/p_tpcSignal", "p_tpcSignal", histPTpcSignal);
    hist.add("QA/tracks/Idfd/Pi/tpcId/tpcInnerParam_tpcSignal", "tpcInnerParam_tpcSignal", histTpcInnerParamTpcSignal);
    // tofBeta
    hist.add("QA/tracks/Idfd/Pi/tpcId/p_beta", "p_beta", histPBeta);
    hist.add("QA/tracks/Idfd/Pi/tpcId/tpcInnerParam_beta", "tpcInnerParam_beta", histTpcInnerParamBeta);
    // Look at Pion
    hist.add("QA/tracks/Idfd/Pi/tpcId/p_tpcNSigma", "p_tpcNSigma", histPTpcNSigma);
    hist.add("QA/tracks/Idfd/Pi/tpcId/pt_tpcNSigma", "pt_tpcNSigma", histPtTpcNSigma);
    hist.add("QA/tracks/Idfd/Pi/tpcId/tpcInnerParam_tpcNSigma", "tpcInnerParam_tpcNSigma", histTpcInnerParamTpcNSigma);
    hist.add("QA/tracks/Idfd/Pi/tpcId/tofExpMom_tpcNSigma", "tofExpMom_tpcNSigma", histTofExpMomTpcNSigma);
    hist.add("QA/tracks/Idfd/Pi/tpcId/p_tofNSigma", "p_tofNSigma", histPTofNSigma);
    hist.add("QA/tracks/Idfd/Pi/tpcId/pt_tofNSigma", "pt_tofNSigma", histPtTofNSigma);
    hist.add("QA/tracks/Idfd/Pi/tpcId/tpcInnerParam_tofNSigma", "tpcInnerParam_tofNSigma", histTpcInnerParamTofNSigma);
    hist.add("QA/tracks/Idfd/Pi/tpcId/tofExpMom_tofNSigma", "tofExpMom_tofNSigma", histTofExpMomTofNSigma);
    hist.add("QA/tracks/Idfd/Pi/tpcId/tpcNSigma_tofNSigma", "tpcNSigma_tofNSigma", histTpcNSigmaTofNSigma);

    hist.addClone("QA/tracks/Idfd/Pi/tpcId/", "QA/tracks/Idfd/Pi/tpctofId/");
    hist.addClone("QA/tracks/Idfd/Pi/", "QA/tracks/Idfd/Ka/");
    hist.addClone("QA/tracks/Idfd/Pi/", "QA/tracks/Idfd/Pr/");

    hist.add("sparse1", "sparse1", qnHist1);
    hist.add("sparse2", "sparse2", qnHist2);

    hist.add("sim/gen/sparse1", "sparse1", qnHist1);
    hist.add("sim/gen/sparse2", "sparse2", qnHist2);

    hist.add("sim/reco/sparse1", "sparse1", qnHist1);
    hist.add("sim/reco/sparse2", "sparse2", qnHist2);

    hist.add("sim/purity/sparse1", "sparse1", qnHist1);
    hist.add("sim/purity/sparse2", "sparse2", qnHist2);

    recoAnalysis.add("recoAnalysis/Pi/Pos/h12_p", "p", kTH1F, {axisP});
    recoAnalysis.add("recoAnalysis/Pi/Pos/h13_pt", "pt", kTH1F, {axisPt});
    recoAnalysis.add("recoAnalysis/Pi/Pos/h14_eta", "eta", kTH1F, {axisEta});
    recoAnalysis.add("recoAnalysis/Pi/Pos/h15_phi", "phi", kTH1F, {axisPhi});
    recoAnalysis.add("recoAnalysis/Pi/Pos/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    recoAnalysis.add("recoAnalysis/Pi/Pos/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});
    recoAnalysis.add("recoAnalysis/Pi/Neg/h12_p", "p", kTH1F, {axisP});
    recoAnalysis.add("recoAnalysis/Pi/Neg/h13_pt", "pt", kTH1F, {axisPt});
    recoAnalysis.add("recoAnalysis/Pi/Neg/h14_eta", "eta", kTH1F, {axisEta});
    recoAnalysis.add("recoAnalysis/Pi/Neg/h15_phi", "phi", kTH1F, {axisPhi});
    recoAnalysis.add("recoAnalysis/Pi/Neg/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    recoAnalysis.add("recoAnalysis/Pi/Neg/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});

    recoAnalysis.addClone("recoAnalysis/Pi/", "recoAnalysis/Ka/");
    recoAnalysis.addClone("recoAnalysis/Pi/", "recoAnalysis/Pr/");
    // recoAnalysis.addClone("recoAnalysis/Pi/", "recoAnalysis/El/");
    // recoAnalysis.addClone("recoAnalysis/Pi/", "recoAnalysis/De/");
    recoAnalysis.add("recoAnalysis/Charge/Pos/h12_p", "p", kTH1F, {axisP});
    recoAnalysis.add("recoAnalysis/Charge/Pos/h13_pt", "pt", kTH1F, {axisPt});
    recoAnalysis.add("recoAnalysis/Charge/Pos/h14_eta", "eta", kTH1F, {axisEta});
    recoAnalysis.add("recoAnalysis/Charge/Pos/h15_phi", "phi", kTH1F, {axisPhi});
    recoAnalysis.add("recoAnalysis/Charge/Pos/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    recoAnalysis.add("recoAnalysis/Charge/Pos/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});
    recoAnalysis.add("recoAnalysis/Charge/Neg/h12_p", "p", kTH1F, {axisP});
    recoAnalysis.add("recoAnalysis/Charge/Neg/h13_pt", "pt", kTH1F, {axisPt});
    recoAnalysis.add("recoAnalysis/Charge/Neg/h14_eta", "eta", kTH1F, {axisEta});
    recoAnalysis.add("recoAnalysis/Charge/Neg/h15_phi", "phi", kTH1F, {axisPhi});
    recoAnalysis.add("recoAnalysis/Charge/Neg/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    recoAnalysis.add("recoAnalysis/Charge/Neg/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});

    genAnalysis.add("genAnalysis/Pi/Pos/h12_p", "p", kTH1F, {axisP});
    genAnalysis.add("genAnalysis/Pi/Pos/h13_pt", "pt", kTH1F, {axisPt});
    genAnalysis.add("genAnalysis/Pi/Pos/h14_eta", "eta", kTH1F, {axisEta});
    genAnalysis.add("genAnalysis/Pi/Pos/h15_phi", "phi", kTH1F, {axisPhi});
    genAnalysis.add("genAnalysis/Pi/Pos/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    genAnalysis.add("genAnalysis/Pi/Pos/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});
    genAnalysis.add("genAnalysis/Pi/Neg/h12_p", "p", kTH1F, {axisP});
    genAnalysis.add("genAnalysis/Pi/Neg/h13_pt", "pt", kTH1F, {axisPt});
    genAnalysis.add("genAnalysis/Pi/Neg/h14_eta", "eta", kTH1F, {axisEta});
    genAnalysis.add("genAnalysis/Pi/Neg/h15_phi", "phi", kTH1F, {axisPhi});
    genAnalysis.add("genAnalysis/Pi/Neg/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    genAnalysis.add("genAnalysis/Pi/Neg/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});

    genAnalysis.addClone("genAnalysis/Pi/", "genAnalysis/Ka/");
    genAnalysis.addClone("genAnalysis/Pi/", "genAnalysis/Pr/");
    // genAnalysis.addClone("genAnalysis/Pi/", "genAnalysis/El/");
    // genAnalysis.addClone("genAnalysis/Pi/", "genAnalysis/De/");

    genAnalysis.add("genAnalysis/Charge/Pos/h12_p", "p", kTH1F, {axisP});
    genAnalysis.add("genAnalysis/Charge/Pos/h13_pt", "pt", kTH1F, {axisPt});
    genAnalysis.add("genAnalysis/Charge/Pos/h14_eta", "eta", kTH1F, {axisEta});
    genAnalysis.add("genAnalysis/Charge/Pos/h15_phi", "phi", kTH1F, {axisPhi});
    genAnalysis.add("genAnalysis/Charge/Pos/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    genAnalysis.add("genAnalysis/Charge/Pos/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});
    genAnalysis.add("genAnalysis/Charge/Neg/h12_p", "p", kTH1F, {axisP});
    genAnalysis.add("genAnalysis/Charge/Neg/h13_pt", "pt", kTH1F, {axisPt});
    genAnalysis.add("genAnalysis/Charge/Neg/h14_eta", "eta", kTH1F, {axisEta});
    genAnalysis.add("genAnalysis/Charge/Neg/h15_phi", "phi", kTH1F, {axisPhi});
    genAnalysis.add("genAnalysis/Charge/Neg/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    genAnalysis.add("genAnalysis/Charge/Neg/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});

    purityAnalysis.add("purityAnalysis/Charge/Pos/h12_p", "p", kTH1F, {axisP});
    purityAnalysis.add("purityAnalysis/Charge/Pos/h13_pt", "pt", kTH1F, {axisPt});
    purityAnalysis.add("purityAnalysis/Charge/Pos/h14_eta", "eta", kTH1F, {axisEta});
    purityAnalysis.add("purityAnalysis/Charge/Pos/h15_phi", "phi", kTH1F, {axisPhi});
    purityAnalysis.add("purityAnalysis/Charge/Pos/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    purityAnalysis.add("purityAnalysis/Charge/Pos/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});
    purityAnalysis.add("purityAnalysis/Charge/Neg/h12_p", "p", kTH1F, {axisP});
    purityAnalysis.add("purityAnalysis/Charge/Neg/h13_pt", "pt", kTH1F, {axisPt});
    purityAnalysis.add("purityAnalysis/Charge/Neg/h14_eta", "eta", kTH1F, {axisEta});
    purityAnalysis.add("purityAnalysis/Charge/Neg/h15_phi", "phi", kTH1F, {axisPhi});
    purityAnalysis.add("purityAnalysis/Charge/Neg/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    purityAnalysis.add("purityAnalysis/Charge/Neg/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});

    purityAnalysis.add("purityAnalysis/Pi/Pos/h12_p", "p", kTH1F, {axisP});
    purityAnalysis.add("purityAnalysis/Pi/Pos/h13_pt", "pt", kTH1F, {axisPt});
    purityAnalysis.add("purityAnalysis/Pi/Pos/h14_eta", "eta", kTH1F, {axisEta});
    purityAnalysis.add("purityAnalysis/Pi/Pos/h15_phi", "phi", kTH1F, {axisPhi});
    purityAnalysis.add("purityAnalysis/Pi/Pos/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    purityAnalysis.add("purityAnalysis/Pi/Pos/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});
    purityAnalysis.add("purityAnalysis/Pi/Neg/h12_p", "p", kTH1F, {axisP});
    purityAnalysis.add("purityAnalysis/Pi/Neg/h13_pt", "pt", kTH1F, {axisPt});
    purityAnalysis.add("purityAnalysis/Pi/Neg/h14_eta", "eta", kTH1F, {axisEta});
    purityAnalysis.add("purityAnalysis/Pi/Neg/h15_phi", "phi", kTH1F, {axisPhi});
    purityAnalysis.add("purityAnalysis/Pi/Neg/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    purityAnalysis.add("purityAnalysis/Pi/Neg/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});

    purityAnalysis.addClone("purityAnalysis/Pi/", "purityAnalysis/Ka/");
    purityAnalysis.addClone("purityAnalysis/Pi/", "purityAnalysis/Pr/");
    // purityAnalysis.addClone("purityAnalysis/Pi/", "purityAnalysis/El/");
    // purityAnalysis.addClone("purityAnalysis/Pi/", "purityAnalysis/De/");

  } // init ends

  static constexpr std::string_view HistRegDire2[] = {
    "v0Table/Full/",
    "v0Table/postK0sCheck/",
    "v0Table/postMassCut/",
    "v0Table/postSelectionCut/",
    "recoK0s/PreSel/",
    "recoK0s/PostSel/",
    "recoTracks/PreSel/",
    "recoTracks/PostSel/",
    "recoAnalysis/",
    "genAnalysis/",
    "purityAnalysis/"};

  enum HistRegEnum2 {
    v0Full = 0,
    v0PostK0sCheck,
    v0PostMassCut,
    v0PostSelectionCut,
    recoK0sPreSel,
    recoK0sPostSel,
    recoTracksPreSel,
    recoTracksPostSel,
    recoAnalysisDir,
    genAnalysisDir,
    purityAnalysisDir
  };

  enum IdentificationType {
    kTPCidentified = 0,
    kTOFidentified,
    kTPCTOFidentified,
    kUnidentified
  };

  enum TpcTofCutType {
    kRectangularCut = 0,
    kCircularCut,
    kEllipsoidalCut
  };

  enum DetEnum {
    tpcId = 0,
    tofId,
    tpctofId,
    NoId
  };

  static constexpr std::string_view DetDire[] = {
    "tpcId/",
    "tofId/",
    "tpctofId/",
    "NoId/"};

  enum HistRegEnum {
    qaEventPreSel = 0,
    qaEventPostSel,
    qaTracksPreSel,
    qaTracksPostSel,
    qaTracksIdfd,
  };

  static constexpr std::string_view HistRegDire[] = {
    "QA/events/preSel/",
    "QA/events/postSel/",
    "QA/tracks/preSel/",
    "QA/tracks/postSel/",
    "QA/tracks/Idfd/"};

  enum ProcessTypeEnum {
    doDataProcessing = 0,
    doRecoProcessing,
    doPurityProcessing,
    doGenProcessing,
    doSimProcessing
  };

  // particle identifications
  //  tpc Selections

  template <int pidMode, typename T>
  bool vetoIdOthersTPC(const T& track)
  {
    if (pidMode != kPi) {
      if (std::fabs(track.tpcNSigmaPi()) < cfgVetoIdCut.cfgVetoId01PiTPC)
        return false;
    }
    if (pidMode != kKa) {
      if (std::fabs(track.tpcNSigmaKa()) < cfgVetoIdCut.cfgVetoId03KaTPC)
        return false;
    }
    if (pidMode != kPr) {
      if (std::fabs(track.tpcNSigmaPr()) < cfgVetoIdCut.cfgVetoId05PrTPC)
        return false;
    }
    if (cfgId02DoElRejection) {
      if (pidMode != kEl) {
        if (std::fabs(track.tpcNSigmaEl()) < cfgVetoIdCut.cfgVetoId07ElTPC)
          return false;
      }
    }
    if (cfgId03DoDeRejection) {
      if (pidMode != kDe) {
        if (std::fabs(track.tpcNSigmaDe()) < cfgVetoIdCut.cfgVetoId09DeTPC)
          return false;
      }
    }
    return true;
  }

  template <int pidMode, typename T>
  bool vetoIdOthersTOF(const T& track)
  {
    if (pidMode != kPi) {
      if (std::fabs(track.tofNSigmaPi()) < cfgVetoIdCut.cfgVetoId02PiTOF)
        return false;
    }
    if (pidMode != kKa) {
      if (std::fabs(track.tofNSigmaKa()) < cfgVetoIdCut.cfgVetoId04KaTOF)
        return false;
    }
    if (pidMode != kPr) {
      if (std::fabs(track.tofNSigmaPr()) < cfgVetoIdCut.cfgVetoId06PrTOF)
        return false;
    }
    if (cfgId02DoElRejection) {
      if (pidMode != kEl) {
        if (std::fabs(track.tofNSigmaEl()) < cfgVetoIdCut.cfgVetoId08ElTOF)
          return false;
      }
    }
    if (cfgId03DoDeRejection) {
      if (pidMode != kDe) {
        if (std::fabs(track.tofNSigmaDe()) < cfgVetoIdCut.cfgVetoId10DeTOF)
          return false;
      }
    }
    return true;
  }

  template <int pidMode, typename T>
  bool vetoIdOthersTPCTOF(const T& track)
  {
    if (pidMode != kPi) {
      if (std::fabs(track.tpcNSigmaPi()) < cfgVetoIdCut.cfgVetoId01PiTPC && std::fabs(track.tofNSigmaPi()) < cfgVetoIdCut.cfgVetoId02PiTOF)
        return false;
    }
    if (pidMode != kKa) {
      if (std::fabs(track.tpcNSigmaKa()) < cfgVetoIdCut.cfgVetoId03KaTPC && std::fabs(track.tofNSigmaKa()) < cfgVetoIdCut.cfgVetoId04KaTOF)
        return false;
    }
    if (pidMode != kPr) {
      if (std::fabs(track.tpcNSigmaPr()) < cfgVetoIdCut.cfgVetoId05PrTPC && std::fabs(track.tofNSigmaPr()) < cfgVetoIdCut.cfgVetoId06PrTOF)
        return false;
    }
    if (cfgId02DoElRejection) {
      if (pidMode != kEl) {
        if (std::fabs(track.tpcNSigmaEl()) < cfgVetoIdCut.cfgVetoId07ElTPC && std::fabs(track.tofNSigmaEl()) < cfgVetoIdCut.cfgVetoId08ElTOF)
          return false;
      }
    }
    if (cfgId03DoDeRejection) {
      if (pidMode != kDe) {
        if (std::fabs(track.tpcNSigmaDe()) < cfgVetoIdCut.cfgVetoId09DeTPC && std::fabs(track.tofNSigmaDe()) < cfgVetoIdCut.cfgVetoId10DeTOF)
          return false;
      }
    }
    return true;
  }

  template <int pidMode, typename T>
  bool selIdRectangularCut(const T& track, const float& nSigmaTPC, const float& nSigmaTOF)
  {
    switch (pidMode) {
      case kPi:
        if (std::fabs(track.tpcNSigmaPi()) < nSigmaTPC &&
            std::fabs(track.tofNSigmaPi()) < nSigmaTOF) {
          return true;
        }
        break;
      case kKa:
        if (std::fabs(track.tpcNSigmaKa()) < nSigmaTPC &&
            std::fabs(track.tofNSigmaKa()) < nSigmaTOF) {
          return true;
        }
        break;
      case kPr:
        if (std::fabs(track.tpcNSigmaPr()) < nSigmaTPC &&
            std::fabs(track.tofNSigmaPr()) < nSigmaTOF) {
          return true;
        }
        break;
      default:
        return false;
        break;
    }
    return false;
  }

  template <int pidMode, typename T>
  bool selIdEllipsoidalCut(const T& track, const float& nSigmaTPC, const float& nSigmaTOF)
  {
    switch (pidMode) {
      case kPi:
        if (std::pow(track.tpcNSigmaPi() / nSigmaTPC, 2) + std::pow(track.tofNSigmaPi() / nSigmaTOF, 2) < 1.0)
          return true;
        break;
      case kKa:
        if (std::pow(track.tpcNSigmaKa() / nSigmaTPC, 2) + std::pow(track.tofNSigmaKa() / nSigmaTOF, 2) < 1.0)
          return true;
        break;
      case kPr:
        if (std::pow(track.tpcNSigmaPr() / nSigmaTPC, 2) + std::pow(track.tofNSigmaPr() / nSigmaTOF, 2) < 1.0)
          return true;
        break;
      default:
        return false;
        break;
    }
    return false;
  }

  template <int pidMode, typename T>
  bool selIdCircularCut(const T& track, const float& nSigmaSquaredRad)
  {
    switch (pidMode) {
      case kPi:
        if (std::pow(track.tpcNSigmaPi(), 2) + std::pow(track.tofNSigmaPi(), 2) < nSigmaSquaredRad)
          return true;
        break;
      case kKa:
        if (std::pow(track.tpcNSigmaKa(), 2) + std::pow(track.tofNSigmaKa(), 2) < nSigmaSquaredRad)
          return true;
        break;
      case kPr:
        if (std::pow(track.tpcNSigmaPr(), 2) + std::pow(track.tofNSigmaPr(), 2) < nSigmaSquaredRad)
          return true;
        break;
      default:
        return false;
        break;
    }
    return false;
  }

  template <typename T>
  bool checkReliableTOF(const T& track)
  {
    if (track.hasTOF())
      return true; // which check makes the information of TOF relaiable? should track.beta() be checked?
    else
      return false;
  }

  template <int pidMode, typename T>
  bool idTPC(const T& track, const float& nSigmaTPC)
  {
    if (cfgId01CheckVetoCut && !vetoIdOthersTPC<pidMode>(track))
      return false;
    switch (pidMode) {
      case kPi:
        if (std::fabs(track.tpcNSigmaPi()) < nSigmaTPC)
          return true;
        break;
      case kKa:
        if (std::fabs(track.tpcNSigmaKa()) < nSigmaTPC)
          return true;
        break;
      case kPr:
        if (std::fabs(track.tpcNSigmaPr()) < nSigmaTPC)
          return true;
        break;
      default:
        return false;
        break;
    }
    return false;
  }

  template <int pidMode, typename T>
  bool idTPCTOF(const T& track, const int& pidCutType, const float& nSigmaTPC, const float& nSigmaTOF, const float& nSigmaSquaredRad)
  {
    if (cfgId01CheckVetoCut && !vetoIdOthersTPCTOF<pidMode>(track))
      return false;
    if (pidCutType == kRectangularCut) {
      return selIdRectangularCut<pidMode>(track, nSigmaTPC, nSigmaTOF);
    } else if (pidCutType == kCircularCut) {
      return selIdCircularCut<pidMode>(track, nSigmaSquaredRad);
    } else if (pidCutType == kEllipsoidalCut) {
      return selIdEllipsoidalCut<pidMode>(track, nSigmaTPC, nSigmaTOF);
    }
    return false;
  }

  template <int pidMode, typename T>
  bool selPdependent(const T& track, int& IdMethod, const float& cfgIdThrPforTOF,
                     const int& idCutTypeLowP, const float& nSigmaTPCLowP, const float& nSigmaTOFLowP, const float& nSigmaRadLowP,
                     const int& idCutTypeHighP, const float& nSigmaTPCHighP, const float& nSigmaTOFHighP, const float& nSigmaRadHighP)
  {
    if (track.p() < cfgIdThrPforTOF) {
      if (checkReliableTOF(track)) {
        if (idTPCTOF<pidMode>(track, idCutTypeLowP, nSigmaTPCLowP, nSigmaTOFLowP, nSigmaRadLowP)) {
          IdMethod = kTPCTOFidentified;
          return true;
        }
        return false;
      } else {
        if (idTPC<pidMode>(track, nSigmaTPCLowP)) {
          IdMethod = kTPCidentified;
          return true;
        }
        return false;
      }
    } else {
      if (checkReliableTOF(track)) {
        if (idTPCTOF<pidMode>(track, idCutTypeHighP, nSigmaTPCHighP, nSigmaTOFHighP, nSigmaRadHighP)) {
          IdMethod = kTPCTOFidentified;
          return true;
        }
        return false;
      }
      return false;
    }
  }

  //
  //______________________________Identification Functions________________________________________________________________
  // Pion
  template <typename T>
  bool selPion(const T& track, int& IdMethod)
  {
    if (cfgId04DoPdependentId) {
      return selPdependent<kPi>(track, IdMethod,
                                cfgIdPi01ThrPforTOF, cfgIdPi02IdCutTypeLowP, cfgIdPi03NSigmaTPCLowP, cfgIdPi04NSigmaTOFLowP, cfgIdPi05NSigmaRadLowP,
                                cfgIdPi06IdCutTypeHighP, cfgIdPi07NSigmaTPCHighP, cfgIdPi08NSigmaTOFHighP, cfgIdPi09NSigmaRadHighP);
    }
    return false;
  }

  // Kaon
  template <typename T>
  bool selKaon(const T& track, int& IdMethod)
  {
    if (cfgId04DoPdependentId) {
      return selPdependent<kKa>(track, IdMethod,
                                cfgIdKa01ThrPforTOF, cfgIdKa02IdCutTypeLowP, cfgIdKa03NSigmaTPCLowP, cfgIdKa04NSigmaTOFLowP, cfgIdKa05NSigmaRadLowP,
                                cfgIdKa06IdCutTypeHighP, cfgIdKa07NSigmaTPCHighP, cfgIdKa08NSigmaTOFHighP, cfgIdKa09NSigmaRadHighP);
    }
    return false;
  }

  // Proton
  template <typename T>
  bool selProton(const T& track, int& IdMethod)
  {
    if (cfgId04DoPdependentId) {
      return selPdependent<kPr>(track, IdMethod,
                                cfgIdPr01ThrPforTOF, cfgIdPr02IdCutTypeLowP, cfgIdPr03NSigmaTPCLowP, cfgIdPr04NSigmaTOFLowP, cfgIdPr05NSigmaRadLowP,
                                cfgIdPr06IdCutTypeHighP, cfgIdPr07NSigmaTPCHighP, cfgIdPr08NSigmaTOFHighP, cfgIdPr09NSigmaRadHighP);
    }
    return false;
  }

  template <int Mode, int pidMode, int detMode, typename H, typename T>
  void fillIdentificationQA(H histReg, const T& track)
  {
    float tpcNSigmaVal = -999, tofNSigmaVal = -999;
    switch (pidMode) {
      case kPi:
        tpcNSigmaVal = track.tpcNSigmaPi();
        tofNSigmaVal = track.tofNSigmaPi();
        break;
      case kKa:
        tpcNSigmaVal = track.tpcNSigmaKa();
        tofNSigmaVal = track.tofNSigmaKa();
        break;
      case kPr:
        tpcNSigmaVal = track.tpcNSigmaPr();
        tofNSigmaVal = track.tofNSigmaPr();
        break;
      case kEl:
        tpcNSigmaVal = track.tpcNSigmaEl();
        tofNSigmaVal = track.tofNSigmaEl();
        break;
      case kDe:
        tpcNSigmaVal = track.tpcNSigmaDe();
        tofNSigmaVal = track.tofNSigmaDe();
        break;
      default:
        tpcNSigmaVal = -999, tofNSigmaVal = -999;
        break;
    }

    // NSigma
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("p_tpcNSigma"), track.p(), tpcNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("pt_tpcNSigma"), track.pt(), tpcNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("tpcInnerParam_tpcNSigma"), track.tpcInnerParam(), tpcNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("tofExpMom_tpcNSigma"), track.tofExpMom(), tpcNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("p_tofNSigma"), track.p(), tofNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("pt_tofNSigma"), track.pt(), tofNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("tpcInnerParam_tofNSigma"), track.tpcInnerParam(), tofNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("tofExpMom_tofNSigma"), track.tofExpMom(), tofNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST(DetDire[detMode]) + HIST("tpcNSigma_tofNSigma"), tpcNSigmaVal, tofNSigmaVal);
  }

  template <int mode, typename T>
  void fillCollQA(const T& col, const int& nCh, const int& nT)
  {
    hist.fill(HIST(HistRegDire[mode]) + HIST("h_VtxZ"), col.posZ());
    hist.fill(HIST(HistRegDire[mode]) + HIST("h_Counts"), 0.5);
    hist.fill(HIST(HistRegDire[mode]) + HIST("multFT0"), col.multFT0C());
    hist.fill(HIST(HistRegDire[mode]) + HIST("centFT0"), col.centFT0M());
    if (mode == qaEventPostSel) {
      hist.fill(HIST(HistRegDire[mode]) + HIST("net_charge"), nCh);
      hist.fill(HIST(HistRegDire[mode]) + HIST("Nt_centFT"), col.centFT0M(), nT);
    }
  }

  template <int mode, typename T>
  void fillTrackQA(const T& track)
  {
    hist.fill(HIST(HistRegDire[mode]) + HIST("h_P"), track.p());
    hist.fill(HIST(HistRegDire[mode]) + HIST("h_P_InnerParameter"), track.tpcInnerParam());
    hist.fill(HIST(HistRegDire[mode]) + HIST("h_Pt"), track.pt());
    hist.fill(HIST(HistRegDire[mode]) + HIST("h_Eta"), track.eta());
    hist.fill(HIST(HistRegDire[mode]) + HIST("h_phi"), track.phi());
    hist.fill(HIST(HistRegDire[mode]) + HIST("h2_Pt_DcaZ"), track.pt(), track.dcaZ());
    hist.fill(HIST(HistRegDire[mode]) + HIST("h2_Pt_DcaXY"), track.pt(), track.dcaXY());
    hist.fill(HIST(HistRegDire[mode]) + HIST("h2_p_DcaZ"), track.p(), track.dcaZ());
    hist.fill(HIST(HistRegDire[mode]) + HIST("h2_p_DcaXY"), track.p(), track.dcaXY());
    hist.fill(HIST(HistRegDire[mode]) + HIST("dE_dx1"), track.p(), track.tpcSignal());
    hist.fill(HIST(HistRegDire[mode]) + HIST("dE_dx2"), track.tpcInnerParam(), track.tpcSignal());
    hist.fill(HIST(HistRegDire[mode]) + HIST("p_pt"), track.p(), track.pt());
    hist.fill(HIST(HistRegDire[mode]) + HIST("p_pInnerParameter"), track.p(), track.tpcInnerParam());

    // tofBeta
    hist.fill(HIST(HistRegDire[mode]) + HIST("p_beta"), track.p(), track.beta());
    hist.fill(HIST(HistRegDire[mode]) + HIST("tpcInnerParam_beta"), track.tpcInnerParam(), track.beta());

    // Look at Pion
    fillIdentificationQA<mode, kPi, NoId>(hist, track); // Look at Pion
    fillIdentificationQA<mode, kKa, NoId>(hist, track); // Look at Kaon
    fillIdentificationQA<mode, kPr, NoId>(hist, track); // Look at Proton
  }

  template <HistRegEnum2 DIR, PidEnum P, ChargeEnum C, typename TrackType>
  void fillGenTrackQA(HistogramRegistry& r, const TrackType& t)
  {
    auto base =
      HIST(HistRegDire2[DIR]) +
      HIST(PidDire[P]) +
      HIST(ChargeDire[C]);
    r.fill(base + HIST("h12_p"), t.p());
    r.fill(base + HIST("h13_pt"), t.pt());
    r.fill(base + HIST("h14_eta"), t.eta());
    r.fill(base + HIST("h15_phi"), t.phi());
    r.fill(base + HIST("h16_rapidity"), t.y());
    r.fill(base + HIST("h20_pt_eta"), t.pt(), t.eta());
  }

  template <HistRegEnum2 DIR, PidEnum P, ChargeEnum C, typename TrackType>
  void fillRecoTrackQA(HistogramRegistry& r, const TrackType& t)
  {
    auto base =
      HIST(HistRegDire2[DIR]) +
      HIST(PidDire[P]) +
      HIST(ChargeDire[C]);
    r.fill(base + HIST("h12_p"), t.p());
    r.fill(base + HIST("h13_pt"), t.pt());
    r.fill(base + HIST("h14_eta"), t.eta());
    r.fill(base + HIST("h15_phi"), t.phi());
    r.fill(base + HIST("h16_rapidity"), t.y());
    r.fill(base + HIST("h20_pt_eta"), t.pt(), t.eta());
  }

  template <HistRegEnum2 DIR, PidEnum P, ChargeEnum C, typename TrackType>
  void fillPurityTrackQA(HistogramRegistry& r, const TrackType& t)
  {
    auto base =
      HIST(HistRegDire2[DIR]) +
      HIST(PidDire[P]) +
      HIST(ChargeDire[C]);
    r.fill(base + HIST("h12_p"), t.p());
    r.fill(base + HIST("h13_pt"), t.pt());
    r.fill(base + HIST("h14_eta"), t.eta());
    r.fill(base + HIST("h15_phi"), t.phi());
    r.fill(base + HIST("h16_rapidity"), t.y());
    r.fill(base + HIST("h20_pt_eta"), t.pt(), t.eta());
  }

  template <typename T, typename H>
  void executeTrackAnalysisPart(const T& track, const int& trackIdTag, float& nP, float& nM,
                                const int& idMethodPi, const bool& trackIsPion, float& nAPi, float& nPi,
                                const int& idMethodKa, const bool& trackIsKaon, float& nAKa, float& nKa,
                                const int& idMethodPr, const bool& trackIsProton, float& nPr, float& nAPr, H& recoAnalysis)
  {
    if (flagUnusedVariableError)
      LOG(info) << trackIdTag << idMethodPi << ":" << idMethodKa << ":" << idMethodPr;
    if (track.sign() > 0) {
      // fillRecoTrackQA<kPos>(recoAnalysis, track);
      nP++;
      recoAnalysis.fill(HIST("recoAnalysis/Charge/Pos/h12_p"), track.p());
      recoAnalysis.fill(HIST("recoAnalysis/Charge/Pos/h13_pt"), track.pt());
      recoAnalysis.fill(HIST("recoAnalysis/Charge/Pos/h14_eta"), track.eta());
      recoAnalysis.fill(HIST("recoAnalysis/Charge/Pos/h15_phi"), track.phi());
      recoAnalysis.fill(HIST("recoAnalysis/Charge/Pos/h16_rapidity"), track.y());
      recoAnalysis.fill(HIST("recoAnalysis/Charge/Pos/h20_pt_eta"), track.pt(), track.eta());
    }
    if (track.sign() < 0) {
      // fillRecoTrackQA<kNeg>(recoAnalysis, track);
      nM++;
      recoAnalysis.fill(HIST("recoAnalysis/Charge/Neg/h12_p"), track.p());
      recoAnalysis.fill(HIST("recoAnalysis/Charge/Neg/h13_pt"), track.pt());
      recoAnalysis.fill(HIST("recoAnalysis/Charge/Neg/h14_eta"), track.eta());
      recoAnalysis.fill(HIST("recoAnalysis/Charge/Neg/h15_phi"), track.phi());
      recoAnalysis.fill(HIST("recoAnalysis/Charge/Neg/h16_rapidity"), track.y());
      recoAnalysis.fill(HIST("recoAnalysis/Charge/Neg/h20_pt_eta"), track.pt(), track.eta());
    }

    if (trackIsPion) {
      // if (idMethodPi == kTPCidentified) {
      //   fillIdentificationQA<recoAnalysisDir, kPi, tpcId, true>(hist, track); // set hist as recoAnalysis after tpcId etc add true
      // } else if (idMethodPi == kTPCTOFidentified) {
      //   fillIdentificationQA<recoAnalysisDir, kPi, tpctofId, true>(hist, track);
      // } else if (idMethodPi == kUnidentified) {
      //   fillIdentificationQA<recoAnalysisDir, kPi, NoId, true>(hist, track);
      // }
      if (track.sign() > 0) {
        nPi++;
        fillRecoTrackQA<recoAnalysisDir, kPi, kPos>(recoAnalysis, track);
      } else if (track.sign() < 0) {
        nAPi++;
        fillRecoTrackQA<recoAnalysisDir, kPi, kNeg>(recoAnalysis, track);
      }
      // fillRecoTrackQA<recoAnalysisDir, kPi>(recoAnalysis, track);
    }
    if (trackIsKaon) {
      // if (idMethodKa == kTPCidentified) {
      //   fillIdentificationQA<recoAnalysisDir, kKa, tpcId, true>(hist, track);
      // } else if (idMethodKa == kTPCTOFidentified) {
      //   fillIdentificationQA<recoAnalysisDir, kKa, tpctofId, true>(hist, track);
      // } else if (idMethodKa == kUnidentified) {
      //   fillIdentificationQA<recoAnalysisDir, kKa, NoId, true>(hist, track);
      // }
      if (track.sign() > 0) {
        nKa++;
        fillRecoTrackQA<recoAnalysisDir, kKa, kPos>(recoAnalysis, track);
      } else if (track.sign() < 0) {
        nAKa++;
        fillRecoTrackQA<recoAnalysisDir, kKa, kNeg>(recoAnalysis, track);
      }
      // fillRecoTrackQA<recoAnalysisDir, kKa>(recoAnalysis, track);
    }
    if (trackIsProton) {
      // if (idMethodPr == kTPCidentified) {
      //   fillIdentificationQA<recoAnalysisDir, kPr, tpcId, true>(hist, track);
      // } else if (idMethodPr == kTPCTOFidentified) {
      //   fillIdentificationQA<recoAnalysisDir, kPr, tpctofId, true>(hist, track);
      // } else if (idMethodPr == kUnidentified) {
      //   fillIdentificationQA<recoAnalysisDir, kPr, NoId, true>(hist, track);
      // }
      if (track.sign() > 0) {
        nPr++;
        fillRecoTrackQA<recoAnalysisDir, kPr, kPos>(recoAnalysis, track);
      } else if (track.sign() < 0) {
        nAPr++;
        fillRecoTrackQA<recoAnalysisDir, kPr, kNeg>(recoAnalysis, track);
      }
      // fillRecoTrackQA<recoAnalysisDir, kPr>(recoAnalysis, track);
    }

    // recoAnalysis.fill(HIST("recoAnalysis/SelectedTrack_IdentificationTag"), trackIdTag);
  }

  using MyAllTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra,
                                aod::TracksDCA, aod::pidTOFbeta, aod::pidTOFmass, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl,
                                aod::pidTPCFullDe, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullDe>;
  using MyCollisions = soa::Join<aod::Collisions, aod::Mults, aod::CentFT0Ms>;

  using MyTracksWithMclabels = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTOFbeta, aod::pidTOFmass,
                                         aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullDe,
                                         aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullDe, aod::McTrackLabels>;
  using MyCollisionsWithMcLabels = soa::Join<aod::Collisions, aod::Mults, aod::CentFT0Ms, aod::McCollisionLabels>;

  // tracks and collision filters
  Filter col = aod::evsel::sel8 == true;
  Filter colFilter = nabs(aod::collision::posZ) < cfgCutPosZ;
  Filter trackFilter = requireGlobalTrackInFilter();
  Filter trackPt = (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax);
  Filter trackDCAxy = nabs(aod::track::dcaXY) < cfgCutDcaXY;
  Filter trackDCAz = nabs(aod::track::dcaZ) < cfgCutDcaZ;
  Filter tracketa = nabs(aod::track::eta) < cfgCutEta;

  using MyFilteredCol = soa::Filtered<MyCollisions>;
  using MyFilteredTracks = soa::Filtered<MyAllTracks>;

  using MyFilteredColsWithMcLabels = soa::Filtered<MyCollisionsWithMcLabels>;
  using MyFilteredTracksWithMcLabels = soa::Filtered<MyTracksWithMclabels>;

  // manual sliceby
  SliceCache cache;

  Preslice<MyFilteredTracks> tracksPerCollisionPreslice = o2::aod::track::collisionId;
  Preslice<MyFilteredTracksWithMcLabels> mctracksPerCollisionPreslice = o2::aod::track::collisionId;

  template <int analysisType, typename C, typename T>
  void executeAnalysis(const C& collisions, const T& tracks)
  {
    float nP = 0;
    float nM = 0;
    float nCh = 0;
    float nT = 0;
    float nPr = 0;
    float nAPr = 0;
    float nKa = 0;
    float nAKa = 0;
    float nPi = 0;
    float nAPi = 0;

    bool trackIsPion = false;
    bool trackIsKaon = false;
    bool trackIsProton = false;

    int trackIdTag = 0;
    int idMethodPi = kUnidentified;
    int idMethodKa = kUnidentified;
    int idMethodPr = kUnidentified;

    int ptEtaBin = -1;

    if constexpr (analysisType == doDataProcessing) {
      for (const auto& col : collisions) {
        nP = 0;
        nM = 0;
        nCh = 0;
        nT = 0;
        nPr = 0;
        nAPr = 0;
        nKa = 0;
        nAKa = 0;
        nPi = 0;
        nAPi = 0;

        // group tracks manually with corresponding collision using col id;
        const uint64_t collIdx = col.globalIndex();
        const auto tracksTablePerColl = tracks.sliceBy(tracksPerCollisionPreslice, collIdx);

        for (const auto& track : tracksTablePerColl) {

          ptEtaBin = hPtEtaForBinSearch->FindBin(track.pt(), track.eta()); // Find Track Bin for efficiency correction

          fillTrackQA<qaTracksPostSel>(track);

          if (track.sign() == 1) {
            // nP++;
            nP += hPtEtaForEffCorrection[kCh][kPos]->GetBinContent(ptEtaBin);
          }
          if (track.sign() == -1) {
            nM += hPtEtaForEffCorrection[kCh][kNeg]->GetBinContent(ptEtaBin);
          }

          int idMethod;
          // pion
          if (selPion(track, idMethod)) {
            if (track.sign() == 1) {
              nPi += hPtEtaForEffCorrection[kPi][kPos]->GetBinContent(ptEtaBin);
            }
            if (track.sign() == -1) {
              nAPi += hPtEtaForEffCorrection[kPi][kNeg]->GetBinContent(ptEtaBin);
            }

            if (idMethod == kTPCidentified)
              fillIdentificationQA<qaTracksIdfd, kPi, tpcId>(hist, track);
            if (idMethod == kTPCTOFidentified)
              fillIdentificationQA<qaTracksIdfd, kPi, tpctofId>(hist, track);
          }
          // kaon
          if (selKaon(track, idMethod)) {
            if (track.sign() == 1) {
              nKa += hPtEtaForEffCorrection[kKa][kPos]->GetBinContent(ptEtaBin);
            }
            if (track.sign() == -1) {
              nAKa += hPtEtaForEffCorrection[kKa][kNeg]->GetBinContent(ptEtaBin);
            }

            if (idMethod == kTPCidentified)
              fillIdentificationQA<qaTracksIdfd, kKa, tpcId>(hist, track);
            if (idMethod == kTPCTOFidentified)
              fillIdentificationQA<qaTracksIdfd, kKa, tpctofId>(hist, track);
          }
          // proton
          if (selProton(track, idMethod)) {
            if (track.sign() == 1) {
              nPr += hPtEtaForEffCorrection[kPr][kPos]->GetBinContent(ptEtaBin);
            }
            if (track.sign() == -1) {
              nAPr += hPtEtaForEffCorrection[kPr][kNeg]->GetBinContent(ptEtaBin);
            }

            if (idMethod == kTPCidentified)
              fillIdentificationQA<qaTracksIdfd, kPr, tpcId>(hist, track);
            if (idMethod == kTPCTOFidentified)
              fillIdentificationQA<qaTracksIdfd, kPr, tpctofId>(hist, track);
          }
        } // track loop ends
        nCh = nP - nM;
        nT = nP + nM;

        fillCollQA<qaEventPostSel>(col, nCh, nT);

        hist.fill(HIST("sparse1"), nCh, nP, nM, nPr, nAPr, nKa, nAKa, nT, col.centFT0M());
        hist.fill(HIST("sparse2"), nCh, nP, nM, nPi, nAPi, nKa, nAKa, nT, col.centFT0M());

      } // collision loop ends
    } else if constexpr (analysisType == doRecoProcessing || analysisType == doPurityProcessing) {

      for (const auto& col : collisions) {

        if (!col.has_mcCollision()) {
          LOG(warning) << "No MC collision for this collision, skip...";
          continue;
        }
        float nP = 0;
        float nM = 0;
        float nCh = 0;
        float nT = 0;
        float nPr = 0;
        float nAPr = 0;
        float nKa = 0;
        float nAKa = 0;
        float nPi = 0;
        float nAPi = 0;
        // group tracks manually with corresponding collision using col id;
        const uint64_t collIdx = col.globalIndex();
        const auto tracksTablePerColl = tracks.sliceBy(mctracksPerCollisionPreslice, collIdx);

        for (const auto& track : tracksTablePerColl) {

          if (!track.has_mcParticle()) {
            LOG(warning) << "No MC Particle for this track, skip...";
            continue;
          }

          auto mcPart = track.mcParticle();

          fillTrackQA<qaTracksPostSel>(track);

          // Do Proper Track Identification
          trackIsPion = false;
          trackIsKaon = false;
          trackIsProton = false;

          trackIdTag = 0;
          idMethodPi = kUnidentified;
          idMethodKa = kUnidentified;
          idMethodPr = kUnidentified;

          if (selPion(track, idMethodPi)) {
            trackIsPion = true;
            BITSET(trackIdTag, ID_BIT_PI);
          }
          if (selKaon(track, idMethodKa)) {
            trackIsKaon = true;
            BITSET(trackIdTag, ID_BIT_KA);
          }
          if (selProton(track, idMethodPr)) {
            trackIsProton = true;
            BITSET(trackIdTag, ID_BIT_PR);
          }

          if constexpr (analysisType == doPurityProcessing) {
            if (trackIsPion) {
              if (track.sign() > 0 && mcPart.pdgCode() != kPiPlus) {
                trackIsPion = false;
              }
              if (track.sign() < 0 && mcPart.pdgCode() != kPiMinus) {
                trackIsPion = false;
              }
            }
            if (trackIsKaon) {
              if (track.sign() > 0 && mcPart.pdgCode() != kKPlus) {
                trackIsKaon = false;
              }
              if (track.sign() < 0 && mcPart.pdgCode() != kKMinus) {
                trackIsKaon = false;
              }
            }
            if (trackIsProton) {
              if (track.sign() > 0 && mcPart.pdgCode() != kProton) {
                trackIsProton = false;
              }
              if (track.sign() < 0 && mcPart.pdgCode() != kProtonBar) {
                trackIsProton = false;
              }
            }
          }

          executeTrackAnalysisPart(track, trackIdTag, nP, nM,
                                   idMethodPi, trackIsPion, nAPi, nPi,
                                   idMethodKa, trackIsKaon, nAKa, nKa,
                                   idMethodPr, trackIsProton, nPr, nAPr, recoAnalysis);

        } // track itteration ends

        nCh = nP - nM;
        nT = nP + nM;

        fillCollQA<qaEventPostSel>(col, nCh, nT);

        hist.fill(HIST("sparse1"), nCh, nP, nM, nPr, nAPr, nKa, nAKa, nT, col.centFT0M());
        hist.fill(HIST("sparse2"), nCh, nP, nM, nPi, nAPi, nKa, nAKa, nT, col.centFT0M());

      } // collision ends
    }
  }

  void processRun3(MyFilteredCol const& collisions, MyFilteredTracks const& tracks)
  {
    executeAnalysis<doDataProcessing>(collisions, tracks);
  }
  PROCESS_SWITCH(NchCumulantsId, processRun3, "Process for Run-3", false);

  void processMcRecco(MyFilteredColsWithMcLabels const& collisions, MyFilteredTracksWithMcLabels const& tracks, aod::McParticles const&)
  {
    executeAnalysis<doRecoProcessing>(collisions, tracks);
  }
  PROCESS_SWITCH(NchCumulantsId, processMcRecco, "Process for MC Recco", false);

  void processPurity(MyFilteredColsWithMcLabels const& collisions, MyFilteredTracksWithMcLabels const& tracks, aod::McParticles const&)
  {
    executeAnalysis<doPurityProcessing>(collisions, tracks);
  }
  PROCESS_SWITCH(NchCumulantsId, processPurity, "Process for MC generated purity", false);
  // process fun is over

  // process for truth and response matrix starts now

  Preslice<aod::McParticles> mcTracksPerMcCollisionPreslice = o2::aod::mcparticle::mcCollisionId;

  using MyMcCollisions = aod::McCollisions;

  void processGen(MyMcCollisions const&, MyFilteredColsWithMcLabels const& collisions, aod::McParticles const& mcParticles)
  {
    for (const auto& col : collisions) {
      if (!col.has_mcCollision()) {
        LOG(warning) << "No MC collision for this event, skip...";
        continue;
      }
      const auto& mcColl = col.mcCollision();

      // ---- apply same Vz cut as data/reco ----
      if (std::abs(mcColl.posZ()) > cfgCutPosZ) {
        continue; // reject GEN events outside acceptance
      }
      // slice MCParticles for this MC collision
      const auto mcTracksTablePerMcColl = mcParticles.sliceBy(mcTracksPerMcCollisionPreslice, mcColl.globalIndex());

      float nP = 0;
      float nM = 0;
      float nCh = 0;
      float nT = 0;
      float nPr = 0;
      float nAPr = 0;
      float nKa = 0;
      float nAKa = 0;
      float nPi = 0;
      float nAPi = 0;
      for (const auto& mcTrack : mcTracksTablePerMcColl) {
        if (!mcTrack.isPhysicalPrimary()) {
          continue;
        }
        // pt / eta acceptance (same as RECO acceptance)
        if (mcTrack.pt() <= cfgCutPtMin ||
            mcTrack.pt() >= cfgCutPtMax ||
            std::abs(mcTrack.eta()) >= cfgCutEta) {
          continue;
        }

        int pdg = mcTrack.pdgCode();

        if (pdg == kPiPlus || pdg == kKPlus || pdg == kProton || pdg == kPositron || pdg == kMuonPlus || pdg == kDeuteron) {
          nP++;

          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h12_p"), mcTrack.p());
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h13_pt"), mcTrack.pt());
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h14_eta"), mcTrack.eta());
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h15_phi"), mcTrack.phi());
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h16_rapidity"), mcTrack.y());
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h20_pt_eta"), mcTrack.pt(), mcTrack.eta());
        } else if (pdg == kPiMinus || pdg == kKMinus || pdg == kProtonBar || pdg == kElectron || pdg == kMuonMinus || pdg == -kDeuteron) {
          nM++;

          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h12_p"), mcTrack.p());
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h13_pt"), mcTrack.pt());
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h14_eta"), mcTrack.eta());
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h15_phi"), mcTrack.phi());
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h16_rapidity"), mcTrack.y());
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h20_pt_eta"), mcTrack.pt(), mcTrack.eta());
        }
        // ----- Pions -----
        if (pdg == kPiPlus) {
          nPi++;
          fillGenTrackQA<genAnalysisDir, kPi, kPos>(genAnalysis, mcTrack);
        } else if (pdg == kPiMinus) {
          nAPi++;
          fillGenTrackQA<genAnalysisDir, kPi, kNeg>(genAnalysis, mcTrack);
        } else if (pdg == kKPlus) {
          nKa++;
          fillGenTrackQA<genAnalysisDir, kKa, kPos>(genAnalysis, mcTrack);
        } else if (pdg == kKMinus) {
          nAKa++;
          fillGenTrackQA<genAnalysisDir, kKa, kNeg>(genAnalysis, mcTrack);
        } else if (pdg == kProton) {
          nPr++;
          fillGenTrackQA<genAnalysisDir, kPr, kPos>(genAnalysis, mcTrack);
        } else if (pdg == kProtonBar) {
          nAPr++;
          fillGenTrackQA<genAnalysisDir, kPr, kNeg>(genAnalysis, mcTrack);
        }
      }
      nCh = nP - nM;
      nT = nP + nM;

      hist.fill(HIST("sparse1"), nCh, nP, nM, nPr, nAPr, nKa, nAKa, nT, col.centFT0M());
      hist.fill(HIST("sparse2"), nCh, nP, nM, nPi, nAPi, nKa, nAKa, nT, col.centFT0M());
    } // collision ends
  }
  PROCESS_SWITCH(NchCumulantsId, processGen, "Process for MC generated or truth", false);

  void processSim(MyFilteredColsWithMcLabels const& collisions, MyFilteredTracksWithMcLabels const& tracks, aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    if (flagUnusedVariableError)
      LOG(info) << mcCollisions.size();
    bool trackIsPion = false;
    bool trackIsKaon = false;
    bool trackIsProton = false;

    int trackIdTag = 0;
    int idMethodPi = kUnidentified;
    int idMethodKa = kUnidentified;
    int idMethodPr = kUnidentified;

    int ptEtaBin = -1;

    for (auto const& col : collisions) {

      if (!col.has_mcCollision()) {
        LOG(warning) << "No MC collision for this collision, skip...";
        continue;
      }
      auto mcCollision = col.mcCollision();

      if (checkCollPosZMc && std::abs(mcCollision.posZ()) > cfgCutPosZ)
        continue;

      // slice reco tracks to this collision
      const uint64_t collIdx = col.globalIndex();
      const auto tracksTablePerColl = tracks.sliceBy(mctracksPerCollisionPreslice, collIdx);

      // slice mc particles to mc collisions
      const auto mcTracksTablePerMcColl = mcParticles.sliceBy(mcTracksPerMcCollisionPreslice, mcCollision.globalIndex());

      // Denominator -- Generator level(truth)

      float nPGen = 0, nMGen = 0, nChGen = 0, nTGen = 0;
      float nPrGen = 0, nAPrGen = 0;
      float nKaGen = 0, nAKaGen = 0;
      float nPiGen = 0, nAPiGen = 0;

      for (const auto& mcTrack : mcTracksTablePerMcColl) {
        if (!mcTrack.isPhysicalPrimary())
          continue;
        if (mcTrack.pt() <= cfgCutPtMin ||
            mcTrack.pt() >= cfgCutPtMax ||
            std::abs(mcTrack.eta()) >= cfgCutEta)
          continue;

        int pdg = mcTrack.pdgCode();

        if (pdg == kPiPlus || pdg == kKPlus || pdg == kProton || pdg == kPositron || pdg == kMuonPlus || pdg == kDeuteron) {
          nPGen++;

          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h12_p"), mcTrack.p());
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h13_pt"), mcTrack.pt());
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h14_eta"), mcTrack.eta());
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h15_phi"), mcTrack.phi());
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h16_rapidity"), mcTrack.y());
          genAnalysis.fill(HIST("genAnalysis/Charge/Pos/h20_pt_eta"), mcTrack.pt(), mcTrack.eta());
        } else if (pdg == kPiMinus || pdg == kKMinus || pdg == kProtonBar || pdg == kElectron || pdg == kMuonMinus || pdg == -kDeuteron) {
          nMGen++;

          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h12_p"), mcTrack.p());
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h13_pt"), mcTrack.pt());
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h14_eta"), mcTrack.eta());
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h15_phi"), mcTrack.phi());
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h16_rapidity"), mcTrack.y());
          genAnalysis.fill(HIST("genAnalysis/Charge/Neg/h20_pt_eta"), mcTrack.pt(), mcTrack.eta());
        }

        if (pdg == kPiPlus) {
          nPiGen++;
          fillGenTrackQA<genAnalysisDir, kPi, kPos>(genAnalysis, mcTrack);
        } else if (pdg == kPiMinus) {
          nAPiGen++;
          fillGenTrackQA<genAnalysisDir, kPi, kNeg>(genAnalysis, mcTrack);
        } else if (pdg == kKPlus) {
          nKaGen++;
          fillGenTrackQA<genAnalysisDir, kKa, kPos>(genAnalysis, mcTrack);
        } else if (pdg == kKMinus) {
          nAKaGen++;
          fillGenTrackQA<genAnalysisDir, kKa, kNeg>(genAnalysis, mcTrack);
        } else if (pdg == kProton) {
          nPrGen++;
          fillGenTrackQA<genAnalysisDir, kPr, kPos>(genAnalysis, mcTrack);
        } else if (pdg == kProtonBar) {
          nAPrGen++;
          fillGenTrackQA<genAnalysisDir, kPr, kNeg>(genAnalysis, mcTrack);
        }

      } // particle/track loop for gen ends
      nChGen = nPGen - nMGen;
      nTGen = nPGen + nMGen;

      // ── Fill GEN sparse (denominator) ────────────────────────
      hist.fill(HIST("sim/gen/sparse1"), nChGen, nPGen, nMGen,
                nPrGen, nAPrGen, nKaGen, nAKaGen, nTGen,
                col.centFT0M());
      hist.fill(HIST("sim/gen/sparse2"), nChGen, nPGen, nMGen,
                nPiGen, nAPiGen, nKaGen, nAKaGen, nTGen,
                col.centFT0M());
      //
      // Numerator - Reconstructed + truth matched
      // reco->selFunc passed, no pdg
      // purity -> selFunc matched with pdg matched
      //
      float nPRec = 0, nMRec = 0, nChRec = 0, nTRec = 0;
      float nPrRec = 0, nAPrRec = 0;
      float nKaRec = 0, nAKaRec = 0;
      float nPiRec = 0, nAPiRec = 0;

      // purity counters — separate from reco
      float nPPur = 0, nMPur = 0, nChPur = 0, nTPur = 0;
      float nPrPur = 0, nAPrPur = 0;
      float nKaPur = 0, nAKaPur = 0;
      float nPiPur = 0, nAPiPur = 0;

      for (const auto& track : tracksTablePerColl) {

        if (!track.has_mcParticle()) {
          LOG(warning) << "No MC Particle for this track, skip...";
          continue;
        }
        auto mcPart = track.mcParticle();
        if (!mcPart.isPhysicalPrimary())
          continue;
        int pdg = mcPart.pdgCode();

        fillTrackQA<qaTracksPostSel>(track);

        trackIsPion = false;
        trackIsKaon = false;
        trackIsProton = false;
        trackIdTag = 0;
        idMethodPi = kUnidentified;
        idMethodKa = kUnidentified;
        idMethodPr = kUnidentified;

        if (selPion(track, idMethodPi)) {
          trackIsPion = true;
          BITSET(trackIdTag, ID_BIT_PI);
        }
        if (selKaon(track, idMethodKa)) {
          trackIsKaon = true;
          BITSET(trackIdTag, ID_BIT_KA);
        }
        if (selProton(track, idMethodPr)) {
          trackIsProton = true;
          BITSET(trackIdTag, ID_BIT_PR);
        }

        // bool isKnownCharged = (std::abs(pdg) == kPiPlus  ||
        //                          std::abs(pdg) == kKPlus   ||
        //                          std::abs(pdg) == kProton  ||
        //                          std::abs(pdg) == kElectron ||
        //                          std::abs(pdg) == kMuonMinus ||
        //                          std::abs(pdg) == kDeuteron);

        ptEtaBin = hPtEtaForBinSearch->FindBin(track.pt(), track.eta()); // Find Track Bin for efficiency correction

        if (track.sign() > 0) {
          nPRec += hPtEtaForEffCorrection[kCh][kPos]->GetBinContent(ptEtaBin);
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Pos/h12_p"), track.p());
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Pos/h13_pt"), track.pt());
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Pos/h14_eta"), track.eta());
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Pos/h15_phi"), track.phi());
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Pos/h16_rapidity"), track.y());
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Pos/h20_pt_eta"), track.pt(), track.eta());
        } else if (track.sign() < 0) {
          nMRec += hPtEtaForEffCorrection[kCh][kNeg]->GetBinContent(ptEtaBin);
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Neg/h12_p"), track.p());
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Neg/h13_pt"), track.pt());
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Neg/h14_eta"), track.eta());
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Neg/h15_phi"), track.phi());
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Neg/h16_rapidity"), track.y());
          recoAnalysis.fill(HIST("recoAnalysis/Charge/Neg/h20_pt_eta"), track.pt(), track.eta());
        }

        // species reco — sel passes, PDG not checked (raw reco)
        if (trackIsPion) {
          if (track.sign() > 0) {
            nPiRec += hPtEtaForEffCorrection[kPi][kPos]->GetBinContent(ptEtaBin);
            fillRecoTrackQA<recoAnalysisDir, kPi, kPos>(recoAnalysis, track);
          } else if (track.sign() < 0) {
            nAPiRec += hPtEtaForEffCorrection[kPi][kNeg]->GetBinContent(ptEtaBin);
            fillRecoTrackQA<recoAnalysisDir, kPi, kNeg>(recoAnalysis, track);
          }
        } else if (trackIsKaon) {
          if (track.sign() > 0) {
            nKaRec += hPtEtaForEffCorrection[kKa][kPos]->GetBinContent(ptEtaBin);
            fillRecoTrackQA<recoAnalysisDir, kKa, kPos>(recoAnalysis, track);
          } else if (track.sign() < 0) {
            nAKaRec += hPtEtaForEffCorrection[kKa][kNeg]->GetBinContent(ptEtaBin);
            fillRecoTrackQA<recoAnalysisDir, kKa, kNeg>(recoAnalysis, track);
          }
        } else if (trackIsProton) {
          if (track.sign() > 0) {
            nPrRec += hPtEtaForEffCorrection[kPr][kPos]->GetBinContent(ptEtaBin);
            fillRecoTrackQA<recoAnalysisDir, kPr, kPos>(recoAnalysis, track);
          } else if (track.sign() < 0) {
            nAPrRec += hPtEtaForEffCorrection[kPr][kNeg]->GetBinContent(ptEtaBin);
            fillRecoTrackQA<recoAnalysisDir, kPr, kNeg>(recoAnalysis, track);
          }
        }
        // purity check - check pdg aginst sign
        bool purityPion = false;
        if (trackIsPion) {
          if (track.sign() > 0 && pdg == kPiPlus)
            purityPion = true;
          if (track.sign() < 0 && pdg == kPiMinus)
            purityPion = true;
        }

        bool purityKaon = false;
        if (trackIsKaon) {
          if (track.sign() > 0 && pdg == kKPlus)
            purityKaon = true;
          if (track.sign() < 0 && pdg == kKMinus)
            purityKaon = true;
        }

        bool purityProton = false;
        if (trackIsProton) {
          if (track.sign() > 0 && pdg == kProton)
            purityProton = true;
          if (track.sign() < 0 && pdg == kProtonBar)
            purityProton = true;
        }

        // charge purity — track.sign() + isKnownCharged + PDG sign consistency
        bool pdgPositive = (pdg == kPiPlus || pdg == kKPlus || pdg == kProton ||
                            pdg == kPositron || pdg == kMuonPlus || pdg == kDeuteron);
        bool pdgNegative = (pdg == kPiMinus || pdg == kKMinus || pdg == kProtonBar ||
                            pdg == kElectron || pdg == kMuonMinus || pdg == -kDeuteron);

        if (track.sign() > 0 && pdgPositive) {
          nPPur++;
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Pos/h12_p"), track.p());
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Pos/h13_pt"), track.pt());
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Pos/h14_eta"), track.eta());
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Pos/h15_phi"), track.phi());
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Pos/h16_rapidity"), track.y());
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Pos/h20_pt_eta"), track.pt(), track.eta());
        } else if (track.sign() < 0 && pdgNegative) {
          nMPur++;
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Neg/h12_p"), track.p());
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Neg/h13_pt"), track.pt());
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Neg/h14_eta"), track.eta());
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Neg/h15_phi"), track.phi());
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Neg/h16_rapidity"), track.y());
          purityAnalysis.fill(HIST("purityAnalysis/Charge/Neg/h20_pt_eta"), track.pt(), track.eta());
        }

        // species purity fill
        if (purityPion) {
          if (track.sign() > 0) {
            nPiPur++;
            fillPurityTrackQA<purityAnalysisDir, kPi, kPos>(purityAnalysis, track);
          } else {
            nAPiPur++;
            fillPurityTrackQA<purityAnalysisDir, kPi, kNeg>(purityAnalysis, track);
          }
        } else if (purityKaon) {
          if (track.sign() > 0) {
            nKaPur++;
            fillPurityTrackQA<purityAnalysisDir, kKa, kPos>(purityAnalysis, track);
          } else {
            nAKaPur++;
            fillPurityTrackQA<purityAnalysisDir, kKa, kNeg>(purityAnalysis, track);
          }
        } else if (purityProton) {
          if (track.sign() > 0) {
            nPrPur++;
            fillPurityTrackQA<purityAnalysisDir, kPr, kPos>(purityAnalysis, track);
          } else {
            nAPrPur++;
            fillPurityTrackQA<purityAnalysisDir, kPr, kNeg>(purityAnalysis, track);
          }
        }
      }
      nChRec = nPRec - nMRec;
      nTRec = nPRec + nMRec;
      nChPur = nPPur - nMPur;
      nTPur = nPPur + nMPur;

      // ── fill reco histos ─────────────────────────────────────
      hist.fill(HIST("sim/reco/sparse1"), nChRec, nPRec, nMRec,
                nPrRec, nAPrRec, nKaRec, nAKaRec, nTRec, col.centFT0M());
      hist.fill(HIST("sim/reco/sparse2"), nChRec, nPRec, nMRec,
                nPiRec, nAPiRec, nKaRec, nAKaRec, nTRec, col.centFT0M());

      // ── fill purity histos ───────────────────────────────────
      hist.fill(HIST("sim/purity/sparse1"), nChPur, nPPur, nMPur,
                nPrPur, nAPrPur, nKaPur, nAKaPur, nTPur, col.centFT0M());
      hist.fill(HIST("sim/purity/sparse2"), nChPur, nPPur, nMPur,
                nPiPur, nAPiPur, nKaPur, nAKaPur, nTPur, col.centFT0M());

    } // common collision loop ends
  } // process sim ends
  PROCESS_SWITCH(NchCumulantsId, processSim, "Process Sim: Gen + Reco + Purity", true);
}; // structure ends

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NchCumulantsId>(cfgc)};
}
