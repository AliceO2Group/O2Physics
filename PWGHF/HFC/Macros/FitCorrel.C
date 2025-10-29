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

/// \file DhCorrelationExtraction.cxx
/// \brief class for D-h correlation extraction
/// \author Samuele Cattaruzzi <samuele.cattaruzzi@cern.ch>
/// \author Swapnesh Santosh Khade <swapnesh.santosh.khade@cern.ch>

#include "DhCorrelationFitter.h"

#include <TAttMarker.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TText.h>

#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>

#include <Rtypes.h>
#include <RtypesCore.h>

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace rapidjson;

bool removeNSPeakLowPt = false;

template <typename ValueType>
void readArray(const Value& jsonArray, vector<ValueType>& output)
{
  for (const auto* it = jsonArray.Begin(); it != jsonArray.End(); it++) {
    auto value = it->template Get<ValueType>();
    output.emplace_back(value);
  }
}

void setTH1HistoStyle(TH1D*& histo, TString hTitle, TString hXaxisTitle, TString hYaxisTitle,
                      Style_t markerStyle, Color_t markerColor, Double_t markerSize,
                      Color_t lineColor, Int_t lineWidth, Float_t hTitleXaxisOffset = 1.3, Float_t hTitleYaxisOffset = 1.3,
                      Float_t hTitleXaxisSize = 0.045, Float_t hTitleYaxisSize = 0.045, Float_t hLabelXaxisSize = 0.045, Float_t hLabelYaxisSize = 0.045,
                      Bool_t centerXaxisTitle = false, Bool_t centerYaxisTitle = false);
void setTH1HistoStyle(TH1F*& histo, TString hTitle, TString hXaxisTitle, TString hYaxisTitle,
                      Style_t markerStyle, Color_t markerColor, Double_t markerSize,
                      Color_t lineColor, Int_t lineWidth, Float_t hTitleXaxisOffset = 1.3, Float_t hTitleYaxisOffset = 1.3,
                      Float_t hTitleXaxisSize = 0.045, Float_t hTitleYaxisSize = 0.045, Float_t hLabelXaxisSize = 0.045, Float_t hLabelYaxisSize = 0.045,
                      Bool_t centerXaxisTitle = false, Bool_t centerYaxisTitle = false);

void fitCorrelDs(const TString cfgFileName = "config_CorrAnalysis.json")
{
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.005);
  gStyle->SetPadBottomMargin(0.2);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLineWidth(2);
  gStyle->SetCanvasDefH(1126);
  gStyle->SetCanvasDefW(1840);

  // Load config
  FILE* configFile = fopen(cfgFileName.Data(), "r");
  Document config;
  char readBuffer[65536];
  FileReadStream is(configFile, readBuffer, sizeof(readBuffer));
  config.ParseStream(is);
  fclose(configFile);

  string codeNameAnalysis = config["CodeName"].GetString();
  gSystem->Exec(Form("rm -rf Output_CorrelationFitting_%s_Root/ Output_CorrelationFitting_%s_png/", codeNameAnalysis.data(), codeNameAnalysis.data()));
  gSystem->Exec(Form("mkdir Output_CorrelationFitting_%s_Root/ Output_CorrelationFitting_%s_png/", codeNameAnalysis.data(), codeNameAnalysis.data()));

  string inputFileNameFit = config["InputFileNameFitCorr"].GetString();
  const TString inFileName = Form("Output_CorrelationExtraction_%s_Root/%s", codeNameAnalysis.data(), inputFileNameFit.data());

  bool const isReflected = config["IsRiflected"].GetBool();
  bool const drawSystematicErrors = config["DrawSystematics"].GetBool();
  bool const sameSystematics = config["SameSystematics"].GetBool();
  bool const shiftBaseUp = config["ShiftBaseUp"].GetBool();
  bool const shiftBaseDown = config["ShiftBaseDown"].GetBool();

  std::vector<double> binsPtCandIntervalsVec;
  std::vector<double> binsPtHadIntervals;
  std::vector<int> fitFunc;

  const Value& ptCandValue = config["binsPtCandIntervals"];
  readArray(ptCandValue, binsPtCandIntervalsVec);

  const Value& ptHadValue = config["binsPtHadIntervals"];
  readArray(ptHadValue, binsPtHadIntervals);

  const int nBinsPtCand = binsPtCandIntervalsVec.size() - 1;
  const int nBinsPtHad = binsPtHadIntervals.size() - 1;

  double binsPtCandIntervals[nBinsPtCand + 1];
  for (int i = 0; i < nBinsPtCand + 1; i++) {
    binsPtCandIntervals[i] = binsPtCandIntervalsVec[i];
  }

  const Value& fitFuncValue = config["FitFunction"];
  readArray(fitFuncValue, fitFunc);

  int const fixBase = config["FixBaseline"].GetInt();
  int const fixMean = config["FixMean"].GetInt();

  int const nBaselinePoints = config["nBaselinePoints"].GetInt();
  vector<int> pointsForBaselineVec;
  const Value& pointsForBaselineValue = config["binsForBaseline"];
  readArray(pointsForBaselineValue, pointsForBaselineVec);
  if (pointsForBaselineVec.size() != nBaselinePoints) {
    cout << "ERROR: size of the vector pointsForBaseline is different from the number of nBaselinePoints" << endl;
    return;
  }
  int pointsForBaseline[nBaselinePoints];
  for (int i = 0; i < nBaselinePoints; i++) {
    pointsForBaseline[i] = pointsForBaselineVec[i];
  }

  std::cout << "=========================== " << std::endl;
  std::cout << "Input variables from config" << std::endl;
  for (int iBinPtCand = 0; iBinPtCand < nBinsPtCand; iBinPtCand++) {
    std::cout << "iPt = " << iBinPtCand + 1 << "  FitFunction    = " << fitFunc[iBinPtCand] << std::endl;
  }
  std::cout << "FixBaseline    = " << fixBase << std::endl;
  std::cout << "FixMean    = " << fixMean << std::endl;
  std::cout << "=========================== " << std::endl;
  std::cout << " " << std::endl;

  // TODO: reflections
  bool const refl = false;

  // Input file
  auto* inFile = new TFile(inFileName.Data());
  auto* inFileSystematicErrors = new TFile("OutputSystematicUncertainties/SystematicUncertaintesAngCorrMerged.root");
  auto* inFileFitSystematicErrors = new TFile("OutputSystematicUncertainties/SystematicUncertaintesFitPhysObsMerged.root");

  // Canvas
  TCanvas* canvasCorrPhi[nBinsPtHad];

  // Histograms
  TH1D* hCorrPhi[nBinsPtCand][nBinsPtHad];
  TH1F* hSystematicErrors[nBinsPtCand][nBinsPtHad];
  TH1D* hSystematicErrorsPlot[nBinsPtCand][nBinsPtHad];

  const int nBinsPtD = 5;
  if (nBinsPtD != nBinsPtCand) {
    std::cout << "[ERROR]: nBinsPtD != nBinsPtCand" << std::endl;
    return;
  }
  double const systUncCorrelatedDs[nBinsPtD] = {20, 20, 20, 10}; // % (just the MC Closure uncertainty to put in the plot)

  // DhCorrelationFitter
  const double fMin{-0.5 * TMath::Pi()}, fMax{1.5 * TMath::Pi()}; // limits for the fitting function
  DhCorrelationFitter* corrFitter[nBinsPtHad][nBinsPtCand];

  // Input parameters for fitting
  const int npars{10}; // PED     NSY  NSM    NSW    ASY   ASM     ASW   BETA    v2D     v2h
  const Double_t vals[npars] = {3., 2., 0., 0.5, 2., 3.14, 0.3, 2., 0.1, 0.1};
  const Double_t lowBounds[npars] = {0., 0., -1., 0., 0., 2., 0., 0.5, 0., 0.};
  const Double_t uppBounds[npars] = {9999., 999., 1., 3.14 / 3., 999., 4., 3.14 / 2., 3.5, 0.5, 0.5};
  const Double_t v2AssocPart[nBinsPtD] = {0.15, 0.15, 0.15, 0.15};
  const Double_t v2Dmeson[nBinsPtD] = {0.175, 0.09, 0.04, 0.04};

  // Output histograms
  TH1D* hBaselin[nBinsPtHad];
  TH1D* hNSYield[nBinsPtHad];
  TH1D* hNSSigma[nBinsPtHad];
  TH1D* hASYield[nBinsPtHad];
  TH1D* hASSigma[nBinsPtHad];
  TH1D* hBeta[nBinsPtHad];
  TH1D* hNSYieldBinCount[nBinsPtHad];
  TH1D* hASYieldBinCount[nBinsPtHad];

  // extract TH1D and prepare fit
  for (int iBinPtCand = 0; iBinPtCand < nBinsPtCand; iBinPtCand++) {
    for (int iBinPtHad = 0; iBinPtHad < nBinsPtHad; iBinPtHad++) {
      if (isReflected) {
        hCorrPhi[iBinPtCand][iBinPtHad] = reinterpret_cast<TH1D*>(inFile->Get(Form("hCorrectedCorrReflected_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", binsPtCandIntervals[iBinPtCand], binsPtCandIntervals[iBinPtCand + 1], binsPtHadIntervals[iBinPtHad], binsPtHadIntervals[iBinPtHad + 1])));
      } else {
        hCorrPhi[iBinPtCand][iBinPtHad] = reinterpret_cast<TH1D*>(inFile->Get(Form("hCorrectedCorr_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", binsPtCandIntervals[iBinPtCand], binsPtCandIntervals[iBinPtCand + 1], binsPtHadIntervals[iBinPtHad], binsPtHadIntervals[iBinPtHad + 1])));
      }

      corrFitter[iBinPtHad][iBinPtCand] = new DhCorrelationFitter(reinterpret_cast<TH1F*>(hCorrPhi[iBinPtCand][iBinPtHad]), fMin, fMax);
      corrFitter[iBinPtHad][iBinPtCand]->setHistoIsReflected(refl);
      corrFitter[iBinPtHad][iBinPtCand]->setFixBaseline(fixBase);
      corrFitter[iBinPtHad][iBinPtCand]->setBaselineUpOrDown(shiftBaseUp, shiftBaseDown);
      corrFitter[iBinPtHad][iBinPtCand]->setPointsForBaseline(nBaselinePoints, pointsForBaseline);
      corrFitter[iBinPtHad][iBinPtCand]->setv2(v2AssocPart[iBinPtCand], v2Dmeson[iBinPtCand]);
      corrFitter[iBinPtHad][iBinPtCand]->setReflectedCorrHisto(isReflected);

      corrFitter[iBinPtHad][iBinPtCand]->setFixMean(fixMean);
      corrFitter[iBinPtHad][iBinPtCand]->setPtRanges(binsPtCandIntervals[iBinPtCand], binsPtCandIntervals[iBinPtCand + 1], binsPtHadIntervals[iBinPtHad], binsPtHadIntervals[iBinPtHad + 1]);
      corrFitter[iBinPtHad][iBinPtCand]->setExternalValsAndBounds(npars, vals, lowBounds, uppBounds); // these are starting points and limits...
    }
  }

  // Plots and fit
  for (int iBinPtHad = 0; iBinPtHad < nBinsPtHad; iBinPtHad++) {
    canvasCorrPhi[iBinPtHad] = new TCanvas(Form("CanvasCorrPhi_PtBinAssoc%d", iBinPtHad + 1), Form("CorrPhiDs_PtBinAssoc%d", iBinPtHad + 1));
    if (nBinsPtCand <= 4) {
      canvasCorrPhi[iBinPtHad]->Divide(2, 2);
    }
    if (nBinsPtCand > 4 && nBinsPtCand <= 6) {
      canvasCorrPhi[iBinPtHad]->Divide(3, 2);
    }
    // histograms with fir parameters
    hBaselin[iBinPtHad] = new TH1D(Form("hBaselin_PtBinAssoc%d", iBinPtHad + 1), "", nBinsPtCand, binsPtCandIntervals);
    hNSYield[iBinPtHad] = new TH1D(Form("hNSYield_PtBinAssoc%d", iBinPtHad + 1), "", nBinsPtCand, binsPtCandIntervals);
    hNSSigma[iBinPtHad] = new TH1D(Form("hNSSigma_PtBinAssoc%d", iBinPtHad + 1), "", nBinsPtCand, binsPtCandIntervals);
    hASYield[iBinPtHad] = new TH1D(Form("hASYield_PtBinAssoc%d", iBinPtHad + 1), "", nBinsPtCand, binsPtCandIntervals);
    hASSigma[iBinPtHad] = new TH1D(Form("hASSigma_PtBinAssoc%d", iBinPtHad + 1), "", nBinsPtCand, binsPtCandIntervals);
    hBeta[iBinPtHad] = new TH1D(Form("hBeta_PtBinAssoc%d", iBinPtHad + 1), "", nBinsPtCand, binsPtCandIntervals);
    hNSYieldBinCount[iBinPtHad] = new TH1D(Form("hNSYieldBinCount_PtBinAssoc%d", iBinPtHad + 1), "", nBinsPtCand, binsPtCandIntervals);
    hASYieldBinCount[iBinPtHad] = new TH1D(Form("hASYieldBinCount_PtBinAssoc%d", iBinPtHad + 1), "", nBinsPtCand, binsPtCandIntervals);

    for (int iBinPtCand = 0; iBinPtCand < nBinsPtCand; iBinPtCand++) {
      setTH1HistoStyle(hCorrPhi[iBinPtCand][iBinPtHad], "", "#Delta#phi [rad]", "#frac{1}{N_{D_{s}}}#frac{dN^{assoc}}{d#Delta#phi} [rad^{-1}]", kFullCircle, kRed + 1, 1.4, kRed + 1, 3);

      canvasCorrPhi[iBinPtHad]->cd(iBinPtCand + 1);
      canvasCorrPhi[iBinPtHad]->SetTickx();
      canvasCorrPhi[iBinPtHad]->SetTicky();
      hCorrPhi[iBinPtCand][iBinPtHad]->SetStats(false);
      hCorrPhi[iBinPtCand][iBinPtHad]->SetMinimum(0);
      // hCorrPhi[iBinPtCand][iBinPtHad] -> Draw();

      // draw systematic errors
      int const nBinsPhi = hCorrPhi[iBinPtCand][iBinPtHad]->GetNbinsX();
      if (drawSystematicErrors) {
        hSystematicErrors[iBinPtCand][iBinPtHad] = reinterpret_cast<TH1F*>(inFileSystematicErrors->Get(Form("hSystematicErrorsMerged_PtBin%d_PtBinAssoc%d", iBinPtCand + 1, iBinPtHad + 1)));
        hSystematicErrorsPlot[iBinPtCand][iBinPtHad] = reinterpret_cast<TH1D*>(hCorrPhi[iBinPtCand][iBinPtHad]->Clone(Form("hSystematicErrorsPlot_PtBin%d_PtBinAssoc%d", iBinPtCand + 1, iBinPtHad + 1)));
        for (int iPhi = 0; iPhi < nBinsPhi; iPhi++) {
          if (sameSystematics) {
            hSystematicErrorsPlot[iBinPtCand][iBinPtHad]->SetBinError(iPhi + 1, hSystematicErrors[iBinPtCand][iBinPtHad]->GetBinContent(1) * hCorrPhi[iBinPtCand][iBinPtHad]->GetBinContent(iPhi + 1));
          } else {
            hSystematicErrorsPlot[iBinPtCand][iBinPtHad]->SetBinError(iPhi + 1, hSystematicErrors[iBinPtCand][iBinPtHad]->GetBinContent(iPhi + 1) * hCorrPhi[iBinPtCand][iBinPtHad]->GetBinContent(iPhi + 1));
          }
          hSystematicErrorsPlot[iBinPtCand][iBinPtHad]->SetBinContent(iPhi + 1, hCorrPhi[iBinPtCand][iBinPtHad]->GetBinContent(iPhi + 1));
        }
        setTH1HistoStyle(hSystematicErrorsPlot[iBinPtCand][iBinPtHad], "", "#Delta#phi [rad]", "#frac{1}{N_{D_{s}}}#frac{dN^{assoc}}{d#Delta#phi} [rad^{-1}]", kFullCircle, kRed + 1, 1.4, kRed + 1, 2);
        hSystematicErrorsPlot[iBinPtCand][iBinPtHad]->SetLineColor(kRed - 4);
        hSystematicErrorsPlot[iBinPtCand][iBinPtHad]->SetFillStyle(0);
        // hSystematicErrorsPlot[iBinPtCand][iBinPtHad] -> Draw("E2same");
      }
      // hCorrPhi[iBinPtCand][iBinPtHad] -> Draw("same");

      // Fit
      corrFitter[iBinPtHad][iBinPtCand]->setFuncType(static_cast<DhCorrelationFitter::FunctionType>(fitFunc[iBinPtCand]));
      corrFitter[iBinPtHad][iBinPtCand]->fitting(kTRUE, kTRUE); // the first term is for drawing the fit functions, the second argument is useExternalParams

      TF1* fFit = corrFitter[iBinPtHad][iBinPtCand]->getFitFunction();

      // Title of the histogram
      auto* pttext = new TPaveText(0.15, 0.9, 0.85, 0.95, "NDC");
      pttext->SetFillStyle(0);
      pttext->SetBorderSize(0);
      TText* tpT = pttext->AddText(0., 0.8, Form("%.0f < p_{T}^{D_{s}} < %.0f GeV/c, p_{T}^{assoc} > %.1f GeV/c", binsPtCandIntervals[iBinPtCand], binsPtCandIntervals[iBinPtCand + 1], binsPtHadIntervals[iBinPtHad]));
      // pttext -> Draw("same");

      // Fill the histograms with the fit parameters
      hBaselin[iBinPtHad]->SetBinContent(iBinPtCand + 1, corrFitter[iBinPtHad][iBinPtCand]->getPedestal());
      hBaselin[iBinPtHad]->SetBinError(iBinPtCand + 1, corrFitter[iBinPtHad][iBinPtCand]->getPedestalError());
      if (iBinPtCand == 0 && removeNSPeakLowPt) {
        hNSYield[iBinPtHad]->SetBinContent(iBinPtCand + 1, -1);
        hNSYield[iBinPtHad]->SetBinError(iBinPtCand + 1, 0);

        hNSSigma[iBinPtHad]->SetBinContent(iBinPtCand + 1, -1);
        hNSSigma[iBinPtHad]->SetBinError(iBinPtCand + 1, 0);

        hBeta[iBinPtHad]->SetBinContent(iBinPtCand + 1, -1);
        hBeta[iBinPtHad]->SetBinError(iBinPtCand + 1, 0);
      } else {
        hNSYield[iBinPtHad]->SetBinContent(iBinPtCand + 1, corrFitter[iBinPtHad][iBinPtCand]->getNsYield());
        hNSYield[iBinPtHad]->SetBinError(iBinPtCand + 1, corrFitter[iBinPtHad][iBinPtCand]->getNsYieldError());

        if (fitFunc[iBinPtCand] != 5 && fitFunc[iBinPtCand] != 6) {
          hNSSigma[iBinPtHad]->SetBinContent(iBinPtCand + 1, corrFitter[iBinPtHad][iBinPtCand]->getNsSigma());
          hNSSigma[iBinPtHad]->SetBinError(iBinPtCand + 1, corrFitter[iBinPtHad][iBinPtCand]->getNsSigmaError());
        } else {
          hNSSigma[iBinPtHad]->SetBinContent(iBinPtCand + 1, TMath::Sqrt(1. / corrFitter[iBinPtHad][iBinPtCand]->getNsSigma()));
          Double_t const errrel = corrFitter[iBinPtHad][iBinPtCand]->getNsSigmaError() / corrFitter[iBinPtHad][iBinPtCand]->getNsSigma() / 2.;
          hNSSigma[iBinPtHad]->SetBinError(iBinPtCand + 1, errrel * TMath::Sqrt(1. / corrFitter[iBinPtHad][iBinPtCand]->getNsSigma()));
        }
      }
      hNSYieldBinCount[iBinPtHad]->SetBinContent(iBinPtCand + 1, corrFitter[iBinPtHad][iBinPtCand]->getBinCountingNsYield());
      hNSYieldBinCount[iBinPtHad]->SetBinError(iBinPtCand + 1, corrFitter[iBinPtHad][iBinPtCand]->getBinCountingNsYieldErr());

      hASYield[iBinPtHad]->SetBinContent(iBinPtCand + 1, corrFitter[iBinPtHad][iBinPtCand]->getAsYield());
      hASYield[iBinPtHad]->SetBinError(iBinPtCand + 1, corrFitter[iBinPtHad][iBinPtCand]->getAsYieldError());

      hASYieldBinCount[iBinPtHad]->SetBinContent(iBinPtCand + 1, corrFitter[iBinPtHad][iBinPtCand]->getBinCountingAsYield());
      hASYieldBinCount[iBinPtHad]->SetBinError(iBinPtCand + 1, corrFitter[iBinPtHad][iBinPtCand]->getBinCountingAsYieldErr());
      if (fitFunc[iBinPtCand] != 5 && fitFunc[iBinPtCand] != 6) {
        hASSigma[iBinPtHad]->SetBinContent(iBinPtCand + 1, corrFitter[iBinPtHad][iBinPtCand]->getAsSigma());
        hASSigma[iBinPtHad]->SetBinError(iBinPtCand + 1, corrFitter[iBinPtHad][iBinPtCand]->getAsSigmaError());
      } else {
        hASSigma[iBinPtHad]->SetBinContent(iBinPtCand + 1, TMath::Sqrt(1. / corrFitter[iBinPtHad][iBinPtCand]->getAsSigma()));
        Double_t const errrel = corrFitter[iBinPtHad][iBinPtCand]->getAsSigmaError() / corrFitter[iBinPtHad][iBinPtCand]->getAsSigma() / 2.;
        hASSigma[iBinPtHad]->SetBinError(iBinPtCand + 1, errrel * TMath::Sqrt(1. / corrFitter[iBinPtHad][iBinPtCand]->getAsSigma()));
      }
      if (fitFunc[iBinPtCand] == 4) { // param beta for gen. gauss
        hBeta[iBinPtHad]->SetBinContent(iBinPtCand + 1, corrFitter[iBinPtHad][iBinPtCand]->getBeta());
        hBeta[iBinPtHad]->SetBinError(iBinPtCand + 1, corrFitter[iBinPtHad][iBinPtCand]->getBetaError());
      }

      // Draw
      auto* tCorrUncDs = new TPaveText(0.413, 0.311, 0.877, 0.392, "NDC");
      tCorrUncDs->SetFillStyle(0);
      tCorrUncDs->SetBorderSize(0);
      tCorrUncDs->SetTextSize(0.05);
      tCorrUncDs->SetTextFont(42);
      tCorrUncDs->SetTextAlign(13);
      tCorrUncDs->SetTextColor(kRed + 1);
      tCorrUncDs->AddText(0., 0., Form("#splitline{+%.0f%%}{#minus%.0f%%}", systUncCorrelatedDs[iBinPtCand], systUncCorrelatedDs[iBinPtCand]));

      auto* tScaleUnc = new TPaveText(0.501, 0.292, 0.968, 0.372, "NDC");
      tScaleUnc->SetFillStyle(0);
      tScaleUnc->SetBorderSize(0);
      tScaleUnc->SetTextSize(0.05);
      tScaleUnc->SetTextFont(42);
      tScaleUnc->SetTextAlign(13);
      tScaleUnc->SetTextColor(kBlack);
      tScaleUnc->AddText(0., 0., "corr. unc.");

      if (drawSystematicErrors) {
        hSystematicErrorsPlot[iBinPtCand][iBinPtHad]->Draw("E2same");
      }
      hCorrPhi[iBinPtCand][iBinPtHad]->Draw("same");
      pttext->Draw("same");
      if (drawSystematicErrors) {
        tCorrUncDs->Draw("same");
        tScaleUnc->Draw("same");
      }
    }
    canvasCorrPhi[iBinPtHad]->SaveAs(Form("Output_CorrelationFitting_%s_png/CorrPhiDs_PtBinAssoc%d.png", codeNameAnalysis.data(), iBinPtHad + 1));
    canvasCorrPhi[iBinPtHad]->SaveAs(Form("Output_CorrelationFitting_%s_Root/CorrPhiDs_PtBinAssoc%d.root", codeNameAnalysis.data(), iBinPtHad + 1));
  }

  // histogram with fit parameter and errors
  auto* outFile = new TFile(Form("Output_CorrelationFitting_%s_Root/CorrPhiDs_FinalPlots.root", codeNameAnalysis.data()), "RECREATE");
  outFile->cd();
  for (int iBinPtHad = 0; iBinPtHad < nBinsPtHad; iBinPtHad++) {
    hBaselin[iBinPtHad]->Write();
    hNSYield[iBinPtHad]->Write();
    hNSSigma[iBinPtHad]->Write();
    hASYield[iBinPtHad]->Write();
    hASSigma[iBinPtHad]->Write();
    hBeta[iBinPtHad]->Write();
    hNSYieldBinCount[iBinPtHad]->Write();
    hASYieldBinCount[iBinPtHad]->Write();
  }
  outFile->Close();

  // Draw plots
  for (int iBinPtHad = 0; iBinPtHad < nBinsPtHad; iBinPtHad++) {
    auto* c1 = new TCanvas(Form("NS_yield_PtAssoc%d", iBinPtHad + 1), Form("NS_yield_PtAssoc%d", iBinPtHad + 1));
    auto* c2 = new TCanvas(Form("AS_yield_PtAssoc%d", iBinPtHad + 1), Form("AS_yield_PtAssoc%d", iBinPtHad + 1));
    auto* c3 = new TCanvas(Form("NS_sigma_PtAssoc%d", iBinPtHad + 1), Form("AS_sigma_PtAssoc%d", iBinPtHad + 1));
    auto* c4 = new TCanvas(Form("AS_sigma_PtAssoc%d", iBinPtHad + 1), Form("AS_sigma_PtAssoc%d", iBinPtHad + 1));
    auto* c5 = new TCanvas(Form("Baseline_PtAssoc%d", iBinPtHad + 1), Form("Baseline_PtAssoc%d", iBinPtHad + 1));
    setTH1HistoStyle(hBaselin[iBinPtHad], Form("p_{T}^{assoc} > %.1f GeV/c", binsPtHadIntervals[iBinPtHad]), "p_{T} (GeV/c)", "Baseline", kFullSquare, kBlue, 1.8, kBlue, 2);
    setTH1HistoStyle(hNSYield[iBinPtHad], Form("p_{T}^{assoc} > %.1f GeV/c", binsPtHadIntervals[iBinPtHad]), "p_{T} (GeV/c)", "Y^{NS}", kFullSquare, kRed, 1.8, kRed, 2);
    setTH1HistoStyle(hASYield[iBinPtHad], Form("p_{T}^{assoc} > %.1f GeV/c", binsPtHadIntervals[iBinPtHad]), "p_{T} (GeV/c)", "Y^{AS}", kFullSquare, kMagenta, 1.8, kMagenta, 2);
    setTH1HistoStyle(hNSSigma[iBinPtHad], Form("p_{T}^{assoc} > %.1f GeV/c", binsPtHadIntervals[iBinPtHad]), "p_{T} (GeV/c)", "#sigma_{NS}", kFullSquare, kOrange + 8, 1.8, kOrange + 8, 2);
    setTH1HistoStyle(hASSigma[iBinPtHad], Form("p_{T}^{assoc} > %.1f GeV/c", binsPtHadIntervals[iBinPtHad]), "p_{T} (GeV/c)", "#sigma_{AS}", kFullSquare, kViolet - 5, 1.8, kViolet - 5, 2);
    c1->cd();
    hNSYield[iBinPtHad]->SetMinimum(0);
    hNSYield[iBinPtHad]->Draw();
    c2->cd();
    hASYield[iBinPtHad]->SetMinimum(0);
    hASYield[iBinPtHad]->Draw();
    c3->cd();
    hNSSigma[iBinPtHad]->SetMinimum(0);
    hNSSigma[iBinPtHad]->Draw();
    c4->cd();
    hASSigma[iBinPtHad]->SetMinimum(0);
    hASSigma[iBinPtHad]->Draw();
    c5->cd();
    hBaselin[iBinPtHad]->SetMinimum(0);
    hBaselin[iBinPtHad]->Draw();

    if (drawSystematicErrors) {
      TH1F* hBaselinSyst = reinterpret_cast<TH1F*>(inFileFitSystematicErrors->Get(Form("hSystematicErrorsBaselinMerged_PtBinAssoc%d", iBinPtHad + 1)));
      TH1F* hNSYieldSyst = reinterpret_cast<TH1F*>(inFileFitSystematicErrors->Get(Form("hSystematicErrorsNSYieldMerged_PtBinAssoc%d", iBinPtHad + 1)));
      TH1F* hNSSigmaSyst = reinterpret_cast<TH1F*>(inFileFitSystematicErrors->Get(Form("hSystematicErrorsNSSigmaMerged_PtBinAssoc%d", iBinPtHad + 1)));
      TH1F* hASYieldSyst = reinterpret_cast<TH1F*>(inFileFitSystematicErrors->Get(Form("hSystematicErrorsASYieldMerged_PtBinAssoc%d", iBinPtHad + 1)));
      TH1F* hASSigmaSyst = reinterpret_cast<TH1F*>(inFileFitSystematicErrors->Get(Form("hSystematicErrorsASSigmaMerged_PtBinAssoc%d", iBinPtHad + 1)));

      for (int iBinPtCand = 0; iBinPtCand < nBinsPtCand; iBinPtCand++) {
        hBaselinSyst->SetBinError(iBinPtCand + 1, hBaselinSyst->GetBinContent(iBinPtCand + 1) * hBaselin[iBinPtHad]->GetBinContent(iBinPtCand + 1));
        hNSYieldSyst->SetBinError(iBinPtCand + 1, hNSYieldSyst->GetBinContent(iBinPtCand + 1) * hNSYield[iBinPtHad]->GetBinContent(iBinPtCand + 1));
        hNSSigmaSyst->SetBinError(iBinPtCand + 1, hNSSigmaSyst->GetBinContent(iBinPtCand + 1) * hNSSigma[iBinPtHad]->GetBinContent(iBinPtCand + 1));
        hASYieldSyst->SetBinError(iBinPtCand + 1, hASYieldSyst->GetBinContent(iBinPtCand + 1) * hASYield[iBinPtHad]->GetBinContent(iBinPtCand + 1));
        hASSigmaSyst->SetBinError(iBinPtCand + 1, hASSigmaSyst->GetBinContent(iBinPtCand + 1) * hASSigma[iBinPtHad]->GetBinContent(iBinPtCand + 1));

        hBaselinSyst->SetBinContent(iBinPtCand + 1, hBaselin[iBinPtHad]->GetBinContent(iBinPtCand + 1));
        hNSYieldSyst->SetBinContent(iBinPtCand + 1, hNSYield[iBinPtHad]->GetBinContent(iBinPtCand + 1));
        hNSSigmaSyst->SetBinContent(iBinPtCand + 1, hNSSigma[iBinPtHad]->GetBinContent(iBinPtCand + 1));
        hASYieldSyst->SetBinContent(iBinPtCand + 1, hASYield[iBinPtHad]->GetBinContent(iBinPtCand + 1));
        hASSigmaSyst->SetBinContent(iBinPtCand + 1, hASSigma[iBinPtHad]->GetBinContent(iBinPtCand + 1));
      }

      c1->cd();
      setTH1HistoStyle(hNSYieldSyst, Form("p_{T}^{assoc} > %.1f GeV/c", binsPtHadIntervals[iBinPtHad]), "p_{T} (GeV/c)", "Y^{NS}", kFullSquare, kRed, 1.8, kRed, 2);
      hNSYieldSyst->SetFillStyle(0);
      hNSYieldSyst->Draw("E2same");

      c2->cd();
      setTH1HistoStyle(hASYieldSyst, Form("p_{T}^{assoc} > %.1f GeV/c", binsPtHadIntervals[iBinPtHad]), "p_{T} (GeV/c)", "Y^{AS}", kFullSquare, kMagenta, 1.8, kMagenta, 2);
      hASYieldSyst->SetFillStyle(0);
      hASYieldSyst->Draw("E2same");

      c3->cd();
      setTH1HistoStyle(hNSSigmaSyst, Form("p_{T}^{assoc} > %.1f GeV/c", binsPtHadIntervals[iBinPtHad]), "p_{T} (GeV/c)", "#sigma_{NS}", kFullSquare, kOrange + 8, 1.8, kOrange + 8, 2);
      hNSSigmaSyst->SetFillStyle(0);
      hNSSigmaSyst->Draw("E2same");

      c4->cd();
      setTH1HistoStyle(hASSigmaSyst, Form("p_{T}^{assoc} > %.1f GeV/c", binsPtHadIntervals[iBinPtHad]), "p_{T} (GeV/c)", "#sigma_{AS}", kFullSquare, kViolet - 5, 1.8, kViolet - 5, 2);
      hASSigmaSyst->SetFillStyle(0);
      hASSigmaSyst->Draw("E2same");

      c5->cd();
      setTH1HistoStyle(hBaselinSyst, Form("p_{T}^{assoc} > %.1f GeV/c", binsPtHadIntervals[iBinPtHad]), "p_{T} (GeV/c)", "Baseline", kFullSquare, kBlue, 1.8, kBlue, 2);
      hBaselinSyst->SetFillStyle(0);
      hBaselinSyst->Draw("E2same");
    }

    c1->SaveAs(Form("Output_CorrelationFitting_%s_png/NearSideYield_PtAssoc%d.png", codeNameAnalysis.data(), iBinPtHad + 1));
    c2->SaveAs(Form("Output_CorrelationFitting_%s_png/AwaySideYield_PtAssoc%d.png", codeNameAnalysis.data(), iBinPtHad + 1));
    c3->SaveAs(Form("Output_CorrelationFitting_%s_png/NearSideSigma_PtAssoc%d.png", codeNameAnalysis.data(), iBinPtHad + 1));
    c4->SaveAs(Form("Output_CorrelationFitting_%s_png/AwaySideSigma_PtAssoc%d.png", codeNameAnalysis.data(), iBinPtHad + 1));
    c5->SaveAs(Form("Output_CorrelationFitting_%s_png/Baseline_PtAssoc%d.png", codeNameAnalysis.data(), iBinPtHad + 1));
    c1->SaveAs(Form("Output_CorrelationFitting_%s_Root/NearSideYield_PtBinAssoc%d.root", codeNameAnalysis.data(), iBinPtHad + 1));
    c2->SaveAs(Form("Output_CorrelationFitting_%s_Root/AwaySideYield_PtBinAssoc%d.root", codeNameAnalysis.data(), iBinPtHad + 1));
    c3->SaveAs(Form("Output_CorrelationFitting_%s_Root/NearSideSigma_PtBinAssoc%d.root", codeNameAnalysis.data(), iBinPtHad + 1));
    c4->SaveAs(Form("Output_CorrelationFitting_%s_Root/AwaySideSigma_PtBinAssoc%d.root", codeNameAnalysis.data(), iBinPtHad + 1));
    c5->SaveAs(Form("Output_CorrelationFitting_%s_Root/Baseline_PtBinAssoc%d.root", codeNameAnalysis.data(), iBinPtHad + 1));
  }
}

void setTH1HistoStyle(TH1D*& histo, TString hTitle, TString hXaxisTitle, TString hYaxisTitle,
                      Style_t markerStyle, Color_t markerColor, Double_t markerSize,
                      Color_t lineColor, Int_t lineWidth, Float_t hTitleXaxisOffset, Float_t hTitleYaxisOffset,
                      Float_t hTitleXaxisSize, Float_t hTitleYaxisSize, Float_t hLabelXaxisSize, Float_t hLabelYaxisSize,
                      Bool_t centerXaxisTitle, Bool_t centerYaxisTitle)
{

  histo->SetTitle(hTitle.Data());
  histo->GetXaxis()->SetTitle(hXaxisTitle.Data());
  histo->GetYaxis()->SetTitle(hYaxisTitle.Data());
  histo->SetMarkerStyle(markerStyle);
  histo->SetMarkerColor(markerColor);
  histo->SetMarkerSize(markerSize);
  histo->SetLineColor(lineColor);
  histo->SetLineWidth(lineWidth);
  histo->GetXaxis()->SetTitleOffset(hTitleXaxisOffset);
  histo->GetYaxis()->SetTitleOffset(hTitleYaxisOffset);
  histo->GetXaxis()->SetTitleSize(hTitleXaxisSize);
  histo->GetYaxis()->SetTitleSize(hTitleYaxisSize);
  histo->GetXaxis()->SetLabelSize(hLabelXaxisSize);
  histo->GetYaxis()->SetLabelSize(hLabelYaxisSize);
  histo->GetXaxis()->CenterTitle(centerXaxisTitle);
  histo->GetYaxis()->CenterTitle(centerYaxisTitle);
}

void setTH1HistoStyle(TH1F*& histo, TString hTitle, TString hXaxisTitle, TString hYaxisTitle,
                      Style_t markerStyle, Color_t markerColor, Double_t markerSize,
                      Color_t lineColor, Int_t lineWidth, Float_t hTitleXaxisOffset, Float_t hTitleYaxisOffset,
                      Float_t hTitleXaxisSize, Float_t hTitleYaxisSize, Float_t hLabelXaxisSize, Float_t hLabelYaxisSize,
                      Bool_t centerXaxisTitle, Bool_t centerYaxisTitle)
{

  histo->SetTitle(hTitle.Data());
  histo->GetXaxis()->SetTitle(hXaxisTitle.Data());
  histo->GetYaxis()->SetTitle(hYaxisTitle.Data());
  histo->SetMarkerStyle(markerStyle);
  histo->SetMarkerColor(markerColor);
  histo->SetMarkerSize(markerSize);
  histo->SetLineColor(lineColor);
  histo->SetLineWidth(lineWidth);
  histo->GetXaxis()->SetTitleOffset(hTitleXaxisOffset);
  histo->GetYaxis()->SetTitleOffset(hTitleYaxisOffset);
  histo->GetXaxis()->SetTitleSize(hTitleXaxisSize);
  histo->GetYaxis()->SetTitleSize(hTitleYaxisSize);
  histo->GetXaxis()->SetLabelSize(hLabelXaxisSize);
  histo->GetYaxis()->SetLabelSize(hLabelYaxisSize);
  histo->GetXaxis()->CenterTitle(centerXaxisTitle);
  histo->GetYaxis()->CenterTitle(centerYaxisTitle);
}
