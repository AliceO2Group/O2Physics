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

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include <algorithm>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics; // for constants
using namespace std;

struct NchCumulantsId {

  HistogramRegistry hist{"hist", {}, OutputObjHandlingPolicy::AnalysisObject};

  // PDG data base
  Service<o2::framework::O2DatabasePDG> pdgDB;

  Configurable<float> cfgCutPosZ{"cfgCutPosZ", 10.0, "cut for vertex Z"};
  Configurable<float> cfgCutDcaXY{"cfgCutDcaXY", 0.12, "cut for dcaXY"};
  Configurable<float> cfgCutDcaZ{"cfgCutDcaZ", 0.3, "cut for dcaZ"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "cut for eta"};
  Configurable<float> cfgCutPtMax{"cfgCutPtMax", 3.0, "max cut for pT"};
  Configurable<float> cfgCutPtMin{"cfgCutPtMin", 0.15, "min cut for pT"};

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

  void init(InitContext const&)
  {
    // QA check axes
    const AxisSpec axisEvents{1, 0, 1, "Counts"};
    const AxisSpec axisEta{100, -1., +1., "#eta"};
    const AxisSpec axisPt{100, 0., 3., "p_{T} (GeV/c)"};
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
    const AxisSpec axisPosCh(1001, -1, 1000, "Pos_charge");
    const AxisSpec axisNegCh(1001, -1, 1000, "Neg_charge");
    const AxisSpec axisNt(5001, -1, 5000, "Mult_midRap_Nch");
    const AxisSpec axisPrCh(101, -1, 100, "Pr_charge");
    const AxisSpec axisAPrCh(101, -1, 100, "APr_charge");
    const AxisSpec axisKaCh(101, -1, 100, "Ka_charge");
    const AxisSpec axisAKaCh(101, -1, 100, "AKa_charge");
    const AxisSpec axisPiCh(1001, -1, 1000, "Pion_Positive");
    const AxisSpec axisAPiCh(1001, -1, 1000, "Pion_Negative");

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
  } // init ends

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
    qaTracksIdfd
  };

  static constexpr std::string_view HistRegDire[] = {
    "QA/events/preSel/",
    "QA/events/postSel/",
    "QA/tracks/preSel/",
    "QA/tracks/postSel/",
    "QA/tracks/Idfd/"};

  enum PidEnum {
    kPi = 0, // dont use kPion, kKaon, as these enumeration
    kKa,     // are already defined in $ROOTSYS/root/include/TPDGCode.h
    kPr,
    kEl,
    kDe
  };

  static constexpr std::string_view PidDire[] = {
    "Pi/",
    "Ka/",
    "Pr/",
    "El/",
    "De/"};

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

  using MyAllTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi,
                                aod::pidTOFPi, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCPr, aod::pidTOFPr, aod::pidTOFbeta, aod::pidTOFmass,
                                aod::pidTPCEl, aod::pidTOFEl, aod::pidTPCDe, aod::pidTOFDe>;
  using MyCollisions = soa::Join<aod::Collisions, aod::Mults, aod::CentFT0Ms>;
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

  // manual sliceby
  SliceCache cache;
  Preslice<MyFilteredTracks> tracksPerCollisionPreslice = o2::aod::track::collisionId;

  void process(MyFilteredCol const& collisions, MyFilteredTracks const& tracks)
  {
    for (const auto& col : collisions) {

      int nP = 0;
      int nM = 0;
      int nCh = 0;
      int nT = 0;
      int nPr = 0;
      int nAPr = 0;
      int nKa = 0;
      int nAKa = 0;
      int nPi = 0;
      int nAPi = 0;
      // group tracks manually with corresponding collision using col id;
      const uint64_t collIdx = col.globalIndex();
      const auto tracksTablePerColl = tracks.sliceBy(tracksPerCollisionPreslice, collIdx);

      for (const auto& track : tracksTablePerColl) {

        fillTrackQA<qaTracksPostSel>(track);

        if (track.sign() == 1) {
          nP++;
        }
        if (track.sign() == -1) {
          nM++;
        }

        int idMethod;
        // pion
        if (selPion(track, idMethod)) {
          if (track.sign() == 1) {
            nPi++;
          }
          if (track.sign() == -1) {
            nAPi++;
          }

          if (idMethod == kTPCidentified)
            fillIdentificationQA<qaTracksIdfd, kPi, tpcId>(hist, track);
          if (idMethod == kTPCTOFidentified)
            fillIdentificationQA<qaTracksIdfd, kPi, tpctofId>(hist, track);
        }
        // kaon
        if (selKaon(track, idMethod)) {
          if (track.sign() == 1) {
            nKa++;
          }
          if (track.sign() == -1) {
            nAKa++;
          }

          if (idMethod == kTPCidentified)
            fillIdentificationQA<qaTracksIdfd, kKa, tpcId>(hist, track);
          if (idMethod == kTPCTOFidentified)
            fillIdentificationQA<qaTracksIdfd, kKa, tpctofId>(hist, track);
        }
        // proton
        if (selProton(track, idMethod)) {
          if (track.sign() == 1) {
            nPr++;
          }
          if (track.sign() == -1) {
            nAPr++;
          }

          if (idMethod == kTPCidentified)
            fillIdentificationQA<qaTracksIdfd, kPr, tpcId>(hist, track);
          if (idMethod == kTPCTOFidentified)
            fillIdentificationQA<qaTracksIdfd, kPr, tpctofId>(hist, track);
        }

      } // track itteration ends
      nCh = nP - nM;
      nT = nP + nM;

      fillCollQA<qaEventPostSel>(col, nCh, nT);

      hist.fill(HIST("sparse1"), nCh, nP, nM, nPr, nAPr, nKa, nAKa, nT, col.centFT0M());
      hist.fill(HIST("sparse2"), nCh, nP, nM, nPi, nAPi, nKa, nAKa, nT, col.centFT0M());

    } // collision ends
  } // process ends
}; // structure ends

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NchCumulantsId>(cfgc)};
}
