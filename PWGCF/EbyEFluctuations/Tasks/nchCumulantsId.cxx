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

/// \file NchCumulantsId.cxx
/// \brief Event by Event conserved charges fluctuations
///        it is meant to be a blank page for further developments.
/// \author Pravata Panigrahi <pravata.panigrahi@cern.ch>:: Sadhana Dash (sadhana@phy.iitb.ac.in) and Rahul Verma (rahul.verma@iitb.ac.in)
#include <algorithm>
#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "TLorentzVector.h"

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

    hist.addClone("QA/tracks/Idfd/Pi/tpcId/", "QA/tracks/Idfd/Pi/tofId/");
    hist.addClone("QA/tracks/Idfd/Pi/", "QA/tracks/Idfd/Ka/");
    hist.addClone("QA/tracks/Idfd/Pi/", "QA/tracks/Idfd/Pr/");

    hist.add("sparse1", "sparse1", qnHist1);
    hist.add("sparse2", "sparse2", qnHist2);
  } // init ends

  // particle identifications
  //  tpc Selections
  template <typename T>
  bool selPionTPCInnerParam(T track)
  {
    if (std::abs(track.tpcNSigmaEl()) > 3.0 && std::abs(track.tpcNSigmaKa()) > 3.0 && std::abs(track.tpcNSigmaPr()) > 3.0 && std::abs(track.tpcNSigmaDe()) > 3.0) {
      if (0.05 <= track.tpcInnerParam() && track.tpcInnerParam() < 0.70 && std::abs(track.tpcNSigmaPi()) < 3.0) {
        return true;
      }
      if (0.70 <= track.tpcInnerParam() && std::abs(track.tpcNSigmaPi()) < 2.0) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selKaonTPCInnerParam(T track)
  {
    // p dependent cuts
    if (std::abs(track.tpcNSigmaEl()) > 3.0 && std::abs(track.tpcNSigmaPi()) > 3.0 && std::abs(track.tpcNSigmaPr()) > 3.0 && std::abs(track.tpcNSigmaDe()) > 3.0) {
      if (0.05 <= track.tpcInnerParam() && track.tpcInnerParam() < 0.70 && std::abs(track.tpcNSigmaKa()) < 3.0) {
        return true;
      }
      if (0.70 <= track.tpcInnerParam() && std::abs(track.tpcNSigmaKa()) < 2.0) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selProtonTPCInnerParam(T track)
  {
    if (std::abs(track.tpcNSigmaEl()) > 3.0 && std::abs(track.tpcNSigmaPi()) > 3.0 && std::abs(track.tpcNSigmaKa()) > 3.0 && std::abs(track.tpcNSigmaDe()) > 3.0) {
      if (0.05 <= track.tpcInnerParam() && track.tpcInnerParam() < 1.60 && std::abs(track.tpcNSigmaPr()) < 3.0) {
        return true;
      }
      if (1.60 <= track.tpcInnerParam() && std::abs(track.tpcNSigmaPr()) < 2.0) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selDeuteronTPCInnerParam(T track)
  {
    if (std::abs(track.tpcNSigmaEl()) > 3.0 && std::abs(track.tpcNSigmaPi()) > 3.0 && std::abs(track.tpcNSigmaKa()) > 3.0 && std::abs(track.tpcNSigmaPr()) > 3.0) {
      if (0.05 <= track.tpcInnerParam() && track.tpcInnerParam() < 1.80 && std::abs(track.tpcNSigmaDe()) < 3.0) {
        return true;
      }
      if (1.80 <= track.tpcInnerParam() && std::abs(track.tpcNSigmaDe()) < 2.0) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selElectronTPCInnerParam(T track)
  {
    if (track.tpcNSigmaEl() < 3.0 && track.tpcNSigmaPi() > 3.0 && track.tpcNSigmaKa() > 3.0 && track.tpcNSigmaPr() > 3.0 && track.tpcNSigmaDe() > 3.0) {
      return true;
    }
    return false;
  }
  //

  // TOF Selections
  // Pion
  template <typename T>
  bool selPionTOF(T track)
  {
    if (track.p() <= 0.75 && std::abs(track.tpcNSigmaPi()) < 3.0 && std::abs(track.tofNSigmaPi()) < 3.0 && std::abs(track.tofNSigmaEl()) > 3.0) {
      return true;
    } else if (0.75 < track.p() // after p = 0.75, Pi and Ka lines of nSigma 3.0 will start intersecting
               && std::abs(track.tpcNSigmaPi()) < 2.0 && std::abs(track.tofNSigmaPi()) < 2.0 && std::abs(track.tofNSigmaEl()) > 3.0 && std::abs(track.tofNSigmaKa()) > 3.0 && std::abs(track.tofNSigmaPr()) > 3.0 && std::abs(track.tofNSigmaDe()) > 3.0) {
      return true;
    }
    return false;
  }

  // Kaon
  template <typename T>
  bool selKaonTOF(T track)
  {
    if (track.p() <= 0.75 && std::abs(track.tpcNSigmaKa()) < 3.0 && std::abs(track.tofNSigmaKa()) < 3.0) {
      return true;
    }
    if (0.75 < track.p() && track.p() <= 1.30 // after 0.75 Pi and Ka lines of nSigma 3.0 will start intersecting
        && std::abs(track.tpcNSigmaKa()) < 3.0 && std::abs(track.tofNSigmaKa()) < 3.0 && std::abs(track.tofNSigmaPi()) > 3.0 && std::abs(track.tofNSigmaEl()) > 3.0) {
      return true;
    }
    if (1.30 < track.p() // after 1.30 Pr and Ka lines of nSigma 3.0 will start intersecting
        && std::abs(track.tpcNSigmaKa()) < 2.0 && std::abs(track.tofNSigmaKa()) < 2.0 && std::abs(track.tofNSigmaEl()) > 3.0 && std::abs(track.tofNSigmaPi()) > 3.0 && std::abs(track.tofNSigmaPr()) > 3.0 && std::abs(track.tofNSigmaDe()) > 3.0) {
      return true;
    }
    return false;
  }

  // Proton
  template <typename T>
  bool selProtonTOF(T track)
  {
    if (track.p() <= 1.30 && std::abs(track.tpcNSigmaPr()) < 3.0 && std::abs(track.tofNSigmaPr()) < 3.0) {
      return true;
    }
    if (1.30 < track.p() && track.p() <= 3.10                                                                                                                                                                                                     // after 1.30 Pr and Ka lines of nSigma 3.0 will start intersecting
        && std::abs(track.tpcNSigmaPr()) < 3.0 && std::abs(track.tofNSigmaPr()) < 3.0 && std::abs(track.tofNSigmaEl()) > 3.0 && std::abs(track.tofNSigmaPi()) > 3.0 && std::abs(track.tofNSigmaKa()) > 3.0 && std::abs(track.tofNSigmaDe()) > 3.0 // Some Deuteron contamination is still coming in p dependent cuts
    ) {
      return true;
    }
    if (3.10 < track.p() // after 3.10 Pr and De lines of nSigma 3.0 will start intersecting
        && std::abs(track.tpcNSigmaPr()) < 2.0 && std::abs(track.tofNSigmaPr()) < 2.0 && std::abs(track.tofNSigmaEl()) > 3.0 && std::abs(track.tofNSigmaPi()) > 3.0 && std::abs(track.tofNSigmaKa()) > 3.0 && std::abs(track.tofNSigmaDe()) > 3.0) {
      return true;
    }
    return false;
  }

  // Deuteron
  template <typename T>
  bool selDeuteronTOF(T track)
  {
    if (track.p() <= 3.10 && std::abs(track.tpcNSigmaDe()) < 3.0 && std::abs(track.tofNSigmaDe()) < 3.0 && std::abs(track.tofNSigmaEl()) > 3.0 && std::abs(track.tofNSigmaPi()) > 3.0 && std::abs(track.tofNSigmaKa()) > 3.0 && std::abs(track.tofNSigmaPr()) > 3.0) {
      return true;
    }
    if (3.10 < track.p() // after 3.10 De and Pr lines of nSigma 3.0 will start intersecting
        && std::abs(track.tpcNSigmaDe()) < 2.0 && std::abs(track.tofNSigmaDe()) < 2.0 && std::abs(track.tofNSigmaEl()) > 3.0 && std::abs(track.tofNSigmaPi()) > 3.0 && std::abs(track.tofNSigmaKa()) > 3.0 && std::abs(track.tofNSigmaPr()) > 3.0) {
      return true;
    }
    return false;
  }

  // Electron
  template <typename T>
  bool selElectronTOF(T track)
  {
    if (
      (std::pow(track.tpcNSigmaEl(), 2) + std::pow(track.tofNSigmaEl(), 2)) < 9.00 && std::abs(track.tofNSigmaPi()) > 3.0 && std::abs(track.tofNSigmaKa()) > 3.0 && std::abs(track.tofNSigmaPr()) > 3.0 && std::abs(track.tofNSigmaDe()) > 3.0) {
      return true;
    }
    return false;
  }
  //

  // SelectionFunctions
  // Pion
  template <typename T>
  bool selPion(T track, int& IdMethod)
  {
    if (!track.hasTOF() && selPionTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() < 0.0 && selPionTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && selPionTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && !selPionTPCInnerParam(track)) {
      IdMethod = 1;
      return selPionTOF(track);
    }
    return false;
  }

  // Kaon
  template <typename T>
  bool selKaon(T track, int& IdMethod)
  {
    if (!track.hasTOF() && selKaonTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() < 0.0 && selKaonTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && selKaonTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && !selKaonTPCInnerParam(track)) {
      IdMethod = 1;
      return selKaonTOF(track);
    }
    return false;
  }

  // Proton
  template <typename T>
  bool selProton(T track, int& IdMethod)
  {
    if (!track.hasTOF() && selProtonTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() < 0.0 && selProtonTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && selProtonTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && !selProtonTPCInnerParam(track)) {
      IdMethod = 1;
      return selProtonTOF(track);
    }
    return false;
  }

  // Deuteron
  template <typename T>
  bool selDeuteron(T track, int& IdMethod)
  {
    if (!track.hasTOF() && selDeuteronTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() < 0.0 && selDeuteronTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && selDeuteronTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && !selDeuteronTPCInnerParam(track)) {
      IdMethod = 1;
      return selDeuteronTOF(track);
    }
    return false;
  }

  // Electron
  template <typename T>
  bool selElectron(T track, int& IdMethod)
  {
    if (!track.hasTOF() && selElectronTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() < 0.0 && selElectronTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && selElectronTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && !selElectronTPCInnerParam(track)) {
      IdMethod = 1;
      return selElectronTOF(track);
    }
    return false;
  }
  //

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

  enum DetEnum {
    tpcId = 0,
    tofId,
    NoId
  };

  static constexpr std::string_view DetDire[] = {
    "tpcId/",
    "tofId/",
    "NoId/"};

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
  Filter trackPt = (aod::track::pt > 0.15f) && (aod::track::pt < 2.0f);
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

          if (idMethod == 0)
            fillIdentificationQA<qaTracksIdfd, kPi, tpcId>(hist, track);
          if (idMethod == 1)
            fillIdentificationQA<qaTracksIdfd, kPi, tofId>(hist, track);
        }
        // kaon
        if (selKaon(track, idMethod)) {
          if (track.sign() == 1) {
            nKa++;
          }
          if (track.sign() == -1) {
            nAKa++;
          }

          if (idMethod == 0)
            fillIdentificationQA<qaTracksIdfd, kKa, tpcId>(hist, track);
          if (idMethod == 1)
            fillIdentificationQA<qaTracksIdfd, kKa, tofId>(hist, track);
        }
        // proton
        if (selProton(track, idMethod)) {
          if (track.sign() == 1) {
            nPr++;
          }
          if (track.sign() == -1) {
            nAPr++;
          }

          if (idMethod == 0)
            fillIdentificationQA<qaTracksIdfd, kPr, tpcId>(hist, track);
          if (idMethod == 1)
            fillIdentificationQA<qaTracksIdfd, kPr, tofId>(hist, track);
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
