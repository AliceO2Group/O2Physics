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

/// \file FitCorrel.C
/// \brief Macro to perform the azimuthal correlation fit
/// \usage .L DhCorrelationFitter.cxx+
/// \usage .x FitCorrel.C("config-file-name")
/// \author Samuele Cattaruzzi <samuele.cattaruzzi@cern.ch>
/// \author Swapnesh Santosh Khade <swapnesh.santosh.khade@cern.ch>

#include "Riostream.h"
#include <TROOT.h>
#include <TF1.h>
#include <TH1D.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TSystem.h>
#include <TStyle.h>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include "DhCorrelationFitter.h"

using namespace rapidjson;

template <typename ValueType>
void readArray(const Value& jsonArray, std::vector<ValueType>& output)
{
  for (auto it = jsonArray.Begin(); it != jsonArray.End(); it++) {
    auto value = it->template Get<ValueType>();
    output.emplace_back(value);
  }
}

void SetTH1HistoStyle(TH1D*& histo, TString hTitle, TString hXaxisTitle, TString hYaxisTitle,
                      Style_t markerStyle, Color_t markerColor, Double_t markerSize,
                      Color_t lineColor, Int_t lineWidth, Float_t hTitleXaxisOffset = 1., Float_t hTitleYaxisOffset = 1.,
                      Float_t hTitleXaxisSize = 0.060, Float_t hTitleYaxisSize = 0.060, Float_t hLabelXaxisSize = 0.060, Float_t hLabelYaxisSize = 0.060,
                      Bool_t centerXaxisTitle = false, Bool_t centerYaxisTitle = false);

void FitCorrel(TString cfgFileName = "config_CorrAnalysis.json")
{
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
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

  string CodeNameAnalysis = config["CodeName"].GetString();
  gSystem->Exec(Form("rm -rf Output_CorrelationFitting_%s_Root/ Output_CorrelationFitting_%s_png/", CodeNameAnalysis.data(), CodeNameAnalysis.data()));
  gSystem->Exec(Form("mkdir Output_CorrelationFitting_%s_Root/ Output_CorrelationFitting_%s_png/", CodeNameAnalysis.data(), CodeNameAnalysis.data()));

  const TString inFileName = Form("Output_CorrelationExtraction_%s_Root/ExtractCorrelationsResults.root", CodeNameAnalysis.data());

  vector<double> binsPtCandIntervals;
  vector<double> binsPtHadIntervals;

  const Value& PtCandValue = config["binsPtCandIntervals"];
  readArray(PtCandValue, binsPtCandIntervals);

  const Value& PtHadValue = config["binsPtHadIntervals"];
  readArray(PtHadValue, binsPtHadIntervals);

  int fitFunc = config["FitFunction"].GetInt();
  int fixBase = config["FixBaseline"].GetInt();
  int fixMean = config["FixMean"].GetInt();

  std::cout << "=========================== " << std::endl;
  std::cout << "Input variables from config" << std::endl;
  std::cout << "FitFunction    = " << fitFunc << std::endl;
  std::cout << "FixBaseline    = " << fixBase << std::endl;
  std::cout << "FixMean    = " << fixMean << std::endl;
  std::cout << "=========================== " << std::endl;
  std::cout << " " << std::endl;

  // TODO: reflections
  bool refl = false;

  const int nBinsPtCand = binsPtCandIntervals.size() - 1;
  const int nBinsPtHad = binsPtHadIntervals.size() - 1;

  // Input file
  TFile* inFile = new TFile(inFileName.Data());

  // Canvas
  TCanvas* CanvasCorrPhi[nBinsPtHad];

  // Histograms
  TH1D* hCorrPhi[nBinsPtCand][nBinsPtHad];

  int nBinsPhi;
  double baselineFromThreePoints[nBinsPtCand][nBinsPtHad], baselineFromThreePointsError[nBinsPtCand][nBinsPtHad];

  // DhCorrelationFitter
  const double fMin{-0.5 * TMath::Pi()}, fMax{1.5 * TMath::Pi()}; // limits for the fitting function
  DhCorrelationFitter* corrFitter[nBinsPtHad][nBinsPtCand];

  // Input parameters for fitting
  const int npars{8}; // PED     NSY  NSM    NSW    ASY   ASM     ASW   BETA
  Double_t vals[npars] = {3., 2., 0., 0.5, 2., 3.14, 0.3, 2.};
  Double_t lowBounds[npars] = {0., 0., -1., 0., 0., 2., 0., 0.5};
  Double_t uppBounds[npars] = {9999., 999., 1., 3.14 / 3., 999., 4., 3.14 / 2., 3.5};

  const int nBaselinePoints{8};
  Int_t pointsForBaseline[nBaselinePoints] = {1, 2, 13, 14, 15, 16, 31, 32};

  // extract TH1D and prepare fit
  for (int iBinPtCand = 0; iBinPtCand < nBinsPtCand; iBinPtCand++) {
    for (int iBinPtHad = 0; iBinPtHad < nBinsPtHad; iBinPtHad++) {
      hCorrPhi[iBinPtCand][iBinPtHad] = reinterpret_cast<TH1D*>(inFile->Get(Form("hCorrectedCorr_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", binsPtCandIntervals[iBinPtCand], binsPtCandIntervals[iBinPtCand + 1], binsPtHadIntervals[iBinPtHad], binsPtHadIntervals[iBinPtHad + 1])));

      corrFitter[iBinPtHad][iBinPtCand] = new DhCorrelationFitter(reinterpret_cast<TH1F*>(hCorrPhi[iBinPtCand][iBinPtHad], fMin, fMax));
      corrFitter[iBinPtHad][iBinPtCand]->SetHistoIsReflected(refl);
      corrFitter[iBinPtHad][iBinPtCand]->SetFuncType(static_cast<DhCorrelationFitter::FunctionType>(fitFunc));
      corrFitter[iBinPtHad][iBinPtCand]->SetFixBaseline(fixBase);
      corrFitter[iBinPtHad][iBinPtCand]->SetPointsForBaseline(nBaselinePoints, pointsForBaseline);

      corrFitter[iBinPtHad][iBinPtCand]->SetFixMean(fixMean);
      corrFitter[iBinPtHad][iBinPtCand]->SetPtRanges(binsPtCandIntervals[iBinPtCand], binsPtCandIntervals[iBinPtCand + 1], binsPtHadIntervals[iBinPtHad], binsPtHadIntervals[iBinPtHad + 1]);
      corrFitter[iBinPtHad][iBinPtCand]->SetExternalValsAndBounds(npars, vals, lowBounds, uppBounds); // these are starting points and limits...
    }
  }

  // Plots and fit
  for (int iBinPtHad = 0; iBinPtHad < nBinsPtHad; iBinPtHad++) {
    CanvasCorrPhi[iBinPtHad] = new TCanvas(Form("CanvasCorrPhi_PtBinAssoc%d", iBinPtHad + 1), Form("CorrPhiDs_PtBinAssoc%d", iBinPtHad + 1));
    if (nBinsPtCand <= 4) {
      CanvasCorrPhi[iBinPtHad]->Divide(2, 2);
    }
    if (nBinsPtCand > 4 && nBinsPtCand <= 6) {
      CanvasCorrPhi[iBinPtHad]->Divide(2, 3);
    }
    for (int iBinPtCand = 0; iBinPtCand < nBinsPtCand; iBinPtCand++) {
      SetTH1HistoStyle(hCorrPhi[iBinPtCand][iBinPtHad], "", "#Delta#phi [rad]", "#frac{1}{N_{D_{s}}}#frac{dN^{assoc}}{d#Delta#phi} [rad^{-1}]", kFullCircle, kRed + 1, 1.4, kRed + 1, 3);

      CanvasCorrPhi[iBinPtHad]->cd(iBinPtCand + 1);
      hCorrPhi[iBinPtCand][iBinPtHad]->Draw();

      // Fit
      corrFitter[iBinPtHad][iBinPtCand]->Fitting(kTRUE, kTRUE); // the first term is for drawing the fit functions, the second argument is useExternalParams

      TF1* fFit = corrFitter[iBinPtHad][iBinPtCand]->GetFitFunction();

      // Title of the histogram
      TPaveText* pttext = new TPaveText(0.15, 0.9, 0.85, 0.95, "NDC");
      pttext->SetFillStyle(0);
      pttext->SetBorderSize(0);
      TText* tpT = pttext->AddText(0., 0.8, Form("%.0f < p_{T}^{D_{s}} < %.0f GeV/c, p_{T}^{assoc} > 0.3 GeV/c", binsPtCandIntervals[iBinPtCand], binsPtCandIntervals[iBinPtCand + 1]));
      pttext->Draw("same");
    }
    CanvasCorrPhi[iBinPtHad]->SaveAs(Form("Output_CorrelationFitting_%s_png/CorrPhi_PtBinAssoc%d.png", CodeNameAnalysis.data(), iBinPtHad + 1));
    CanvasCorrPhi[iBinPtHad]->SaveAs(Form("Output_CorrelationFitting_%s_Root/CorrPhi_PtBinAssoc%d.root", CodeNameAnalysis.data(), iBinPtHad + 1));
  }

  return;
}

void SetTH1HistoStyle(TH1D*& histo, TString hTitle, TString hXaxisTitle, TString hYaxisTitle,
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

  return;
}
