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
/// \file saveCorrelation.C
/// \brief
/// \author ALICE

#include "multCalibrator.h"
#include "multGlauberNBDFitter.h"

#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <iostream>

//________________________________________________________________
Double_t FastIntegrate(TF1* f1, Double_t a, Double_t b, Int_t n = 5)
{
  // Do fast integration with N sampling points
  const Int_t nc = n;
  Double_t x[nc], y[nc];
  Double_t lWidth = (b - a) / ((double)(n - 1));
  for (Int_t ii = 0; ii < n; ii++) {
    x[ii] = a + ((double)(ii)) * lWidth;
    y[ii] = f1->Eval(x[ii]);
  }
  // Now go via trapezoids, please (this probably has a name)
  Double_t lIntegral = 0;
  for (Int_t ii = 0; ii < n - 1; ii++) {
    lIntegral += 0.5 * lWidth * (y[ii] + y[ii + 1]);
  }
  return lIntegral / (b - a);
}

Double_t GetBoundaryForPercentile(TH1* histo, Double_t lPercentileRequested)
{
  // This function returns the boundary for a specific percentile.
  Double_t lReturnValue = 0.0;
  Double_t lPercentile = 100.0 - lPercentileRequested;

  const Long_t lNBins = histo->GetNbinsX();
  Double_t lCountDesired = lPercentile * histo->GetEntries() / 100;
  Long_t lCount = 0;
  for (Long_t ibin = 1; ibin < lNBins; ibin++) {
    lCount += histo->GetBinContent(ibin);
    if (lCount >= lCountDesired) {
      // Found bin I am looking for!
      Double_t lWidth = histo->GetBinWidth(ibin);
      Double_t lLeftPercentile = 100. * (lCount - histo->GetBinContent(ibin)) / histo->GetEntries();
      Double_t lRightPercentile = 100. * lCount / histo->GetEntries();

      Double_t lProportion = (lPercentile - lLeftPercentile) / (lRightPercentile - lLeftPercentile);

      lReturnValue = histo->GetBinLowEdge(ibin) + lProportion * lWidth;
      break;
    }
  }
  return lReturnValue;
}

/// @brief master glauber fit function
/// @param lInputFileName input file name (from the grid, typically centrality-studies task)
/// @param histogramName histogram name: histogram to use within the input file
/// @param ancestorMode ancestor mode: 0: truncation, 1: rounding, 2: effective / non-integer (default: 2)
/// @param lFreek free k: keep k value free (default Pb-Pb: fixed at 1.5)
/// @param use_dMu_dNanc use dMu/dNanc: allow for a varying production of Nch vs ancestor if Nancestor is large. Models detector saturation in effective manner.
/// @param lFreef free f: keep f value free (default Pb-Pb: fixed at 0.8)
/// @param lfvalue f value: the value to use for fixed f
/// @param outputFile name of output file
int runGlauberFit(TString lInputFileName = "AnalysisResultsLHC24ar.root", TString histogramName = "hFT0C_BCs", int ancestorMode = 2, Bool_t lFreek = kFALSE, Bool_t use_dMu_dNanc = kFALSE, Bool_t lFreef = kFALSE, Float_t lfvalue = 0.800, TString outputFile = "output.root")
{
  gStyle->SetLineScalePS(1);
  gStyle->SetOptStat(0);

  cout << "Starting!" << endl;
  TFile* file = new TFile(lInputFileName.Data(), "READ");
  if (!file)
    cout << "Problem with file!" << endl;
  TH1F* hV0Mfine = 0x0;

  hV0Mfine = (TH1F*)file->Get(Form("centrality-study/%s", histogramName.Data()));

  // disregard bin zero
  cout << "Received bin zero content: " << hV0Mfine->GetBinContent(0) << ", will set to zero..." << endl;
  hV0Mfine->SetBinContent(0, 0);

  if (!hV0Mfine)
    cout << "Problem with histogram!" << endl;

  cout << "Input histogram has been received successfully! Information: " << endl;

  cout << "Counts: " << hV0Mfine->GetEntries() << endl;
  cout << "NbinsX: " << hV0Mfine->GetNbinsX() << endl;
  cout << "MaxX: " << hV0Mfine->GetBinLowEdge(hV0Mfine->GetNbinsX() + 1) << endl;

  cout << "Creating output file..." << endl;
  TString lProcessedFileName = lInputFileName.Data();
  TString lkMode = "fixedK";
  if (lFreek)
    lkMode = "freeK";
  TString lSaturationMode = "fixedMu";
  if (use_dMu_dNanc)
    lSaturationMode = "freeMu";

  TFile* fOutput = new TFile(outputFile.Data(), "RECREATE");

  TH1F* hV0M = (TH1F*)hV0Mfine->Clone("hV0M");
  TH1F* hV0MUltraFine = (TH1F*)hV0Mfine->Clone("hV0MUltraFine");
  hV0M->SetName("hV0M");
  hV0M->SetTitle("");

  //____________________________________________
  // maximum fit range estimate (avoid tails)
  // may need adjusting
  Double_t lFitRangeMax = GetBoundaryForPercentile(hV0Mfine, 0.008);
  cout << "Fit range max estimated from histogram: " << lFitRangeMax << endl;

  // minimum fit range estimate (guess region that may be unfittable)
  // may need adjusting
  Double_t maxPercent = 0.01;
  Double_t fractionOfMax = 0.012;
  Double_t lFitRange = fractionOfMax * GetBoundaryForPercentile(hV0Mfine, maxPercent);

  // adjust if low mult (Ntracks, typically)
  Double_t maxRangeForTracks = 10000;
  Double_t fractionOfMaxBroader = 0.02;
  if (lFitRangeMax < maxRangeForTracks)
    lFitRange = fractionOfMaxBroader * GetBoundaryForPercentile(hV0Mfine, maxPercent);

  cout << "Fit range min estimated from histogram: " << lFitRange << endl;

  //____________________________________________
  // rebinning matters
  int rebinFactor = 20;
  if (lFitRangeMax < maxRangeForTracks)
    rebinFactor = 1;

  cout << "Creating rebinned histogram with rebin factor: " << rebinFactor << endl;
  hV0M->Rebin(rebinFactor);

  //____________________________________________
  // simple plots for inspection
  TCanvas* c1 = new TCanvas("c1", "", 1300, 900);
  c1->SetFrameFillStyle(0);
  c1->SetFillStyle(0);
  c1->Divide(1, 2);
  c1->cd(1)->SetFrameFillStyle(0);
  c1->cd(1)->SetFillStyle(0);
  c1->cd(2)->SetFrameFillStyle(0);
  c1->cd(2)->SetFillStyle(0);

  c1->cd(1);
  c1->cd(1)->SetLogy();
  c1->cd(1)->SetTicks(1, 1);
  c1->cd(1)->SetPad(0, 0.5, 1, 1);
  c1->cd(2)->SetPad(0, 0.0, 1, .5);

  c1->cd(1)->SetBottomMargin(0.001);
  c1->cd(1)->SetRightMargin(0.25);
  c1->cd(1)->SetTopMargin(0.02);
  c1->cd(1)->SetLeftMargin(0.07);

  c1->cd(2)->SetBottomMargin(0.14);
  c1->cd(2)->SetRightMargin(0.25);
  c1->cd(2)->SetTopMargin(0.001);
  c1->cd(2)->SetLeftMargin(0.07);
  c1->cd(2)->SetTicks(1, 1);
  c1->cd(1);

  hV0M->GetXaxis()->SetRangeUser(0, lFitRangeMax);
  hV0M->GetYaxis()->SetRangeUser(0.25, hV0M->GetMaximum() * 3);
  hV0M->SetLineColor(kBlack);
  hV0M->SetMarkerStyle(20);
  hV0M->SetMarkerColor(kBlack);
  hV0M->SetMarkerSize(0.5);
  hV0M->GetYaxis()->SetTitleSize(0.07);
  hV0M->GetYaxis()->SetLabelSize(0.05);
  hV0M->GetYaxis()->SetTitle("Count");
  hV0M->GetYaxis()->SetTitleOffset(0.5);
  hV0M->GetXaxis()->SetLabelSize(0.05);
  hV0M->GetXaxis()->SetTitleSize(0.06);
  hV0M->GetXaxis()->SetTitle("FT0A+C Amplitude");
  hV0M->GetYaxis()->SetTickLength(0.015);
  hV0M->SetStats(0);
  hV0M->Draw("E");

  // Stand back! Imma gonna do GLAUBER FITTIN'
  multGlauberNBDFitter* g = new multGlauberNBDFitter("lglau");
  g->SetAncestorMode(ancestorMode);

  // Step 1: open the (Npart, Ncoll) pair information, provide
  TFile* fbasefile = new TFile("basehistos.root", "READ");
  TH2D* hNpNc = (TH2D*)fbasefile->Get("hNpNc");
  // return to proper scope
  fOutput->cd();
  g->SetNpartNcollCorrelation(hNpNc);
  g->SetInputV0M(hV0M);
  g->SetFitRange(lFitRange, lFitRangeMax);
  // Step 3: go for it ...
  g->SetNorm(1.53527e+08);
  TString lString = "REM0";
  g->SetFitOptions(lString.Data());
  g->SetFitNpx(100000);

  TF1* fitfunc = g->GetGlauberNBD();

  //____________________________________________
  //
  // set initial fit parameters here
  // may require manual tuning depending on data!

  Double_t guessedMu = lFitRangeMax / 53968.4 * 0.175 * 3.53971e+02;
  cout << "Guessed GlauberNBD mu value: " << guessedMu << endl;

  fitfunc->SetParameter(0, guessedMu); // mu value
  fitfunc->SetParLimits(0, 0.25 * guessedMu, guessedMu * 2);

  if (!lFreek) {
    fitfunc->FixParameter(1, 1.5); // k value
  } else {
    fitfunc->SetParLimits(1, 0.01, 35);
    fitfunc->SetParameter(1, 1.5);
  }
  if (!lFreef) {
    fitfunc->FixParameter(2, lfvalue); // f value
  } else {
    fitfunc->SetParLimits(2, 0.55, 0.97);
    fitfunc->SetParameter(2, 0.800);
  }
  fitfunc->SetParLimits(3, 0.1e+5, 800e+10); // normalization
  fitfunc->SetParameter(3, 0.80832e+08);

  // dMu/dNanc
  fitfunc->SetParameter(4, 0);
  if (!use_dMu_dNanc) {
    fitfunc->FixParameter(4, 0);
  }

  // fitfunc->SetParameter(4,-1.15443e-01);
  // fitfunc->SetParameter(4,-4.55957e-02);

  // dk/dNanc
  fitfunc->FixParameter(5, 0);
  // fitfunc->SetParameter(5,1.63590e-03);

  // d2Mu/dNanc2
  fitfunc->FixParameter(6, 0.0);
  // fitfunc->SetParameter(6,4.02271e-05);

  fitfunc->FixParameter(7, 0.0);
  // fitfunc->SetParameter(7,-1.24349e-06);

  //____________________________________________
  // handle internals for fitting: needs to
  // be done before the fit is attempted!
  g->InitializeNpNc();
  g->InitAncestor();

  cout << "WILL NOW ATTEMPT GLAUBER FIT" << endl;
  cout << "This will take a while. Please wait..." << endl;
  Int_t lFitStatus = 0;
  lFitStatus = g->DoFit();
  Int_t lAttempts = 1;
  Int_t lMaxAttempts = 10;
  while (lAttempts < lMaxAttempts && lFitStatus == 0) {
    // insist on fitting until it works
    cout << "Attempting fit again (" << lAttempts << " attempt)..." << endl;
    lFitStatus = g->DoFit();
  }
  cout << "Final fit status: " << lFitStatus << endl;

  gStyle->SetOptStat(0);
  // Do a ratio plot
  TH1D* hGlauber = (TH1D*)hV0MUltraFine->Clone("hGlauber");
  TH1D* hRatio = (TH1D*)hV0MUltraFine->Clone("hRatio");
  hGlauber->Reset();

  c1->cd(1);
  fitfunc->SetLineColor(kRed);
  fitfunc->SetLineWidth(2);
  fitfunc->Draw("same");

  cout << "Calculating glauber function histogram with the same binning as data input... please wait..." << endl;
  for (Int_t ii = 1; ii < hGlauber->GetNbinsX() + 1; ii++) {
    Double_t lFuncVal = FastIntegrate(fitfunc, hGlauber->GetBinLowEdge(ii), hGlauber->GetBinLowEdge(ii + 1), 4);
    hGlauber->SetBinContent(ii, lFuncVal);
    Int_t printEveryThisManyBins = 500;
    if (ii % printEveryThisManyBins == 0) {
      cout << "At integration #" << ii << "/" << hGlauber->GetNbinsX() + 1 << "..." << endl;
    }
  }
  cout << "Glauber function evaluated. Should go quickly now." << endl;

  c1->cd(2);
  Float_t lLoRangeRatio = 0.35;
  Float_t lHiRangeRatio = 1.65;
  hRatio->Divide(hGlauber);
  hRatio->GetYaxis()->SetTitle("Data/Fit");
  hRatio->GetXaxis()->SetTitle("FT0C amplitude");
  hRatio->GetYaxis()->SetTitleSize(0.055);
  hRatio->GetYaxis()->SetTitleOffset(0.7);
  hRatio->GetXaxis()->SetTitleSize(0.055);
  hRatio->GetYaxis()->SetLabelSize(0.045);
  hRatio->GetXaxis()->SetLabelSize(0.045);
  hRatio->GetYaxis()->SetRangeUser(lLoRangeRatio, lHiRangeRatio);
  hRatio->GetXaxis()->SetRangeUser(0, lFitRangeMax);
  hRatio->SetMarkerStyle(20);
  hRatio->SetMarkerColor(kGray + 2);
  hRatio->SetLineColor(kGray + 2);
  // hRatio->SetMarkerSize(1.0);
  hRatio->SetMarkerSize(.7);
  hRatio->SetStats(0);
  // hRatioWide->SetStats(0);

  hRatio->Draw("hist");

  TLine* line = new TLine(0, 1, lFitRangeMax, 1);
  line->SetLineStyle(7);
  line->SetLineColor(kGray + 1);
  line->Draw();

  TLine* lFitRangeLine = new TLine(lFitRange, lLoRangeRatio, lFitRange, 0.9);
  lFitRangeLine->SetLineColor(kBlue);
  lFitRangeLine->SetLineWidth(1);
  lFitRangeLine->SetLineStyle(2);
  lFitRangeLine->Draw();

  TH1D* hRatioGrayed = (TH1D*)hRatio->Clone("hRatioGrayed");
  hRatioGrayed->SetMarkerColor(kGray + 2);
  hRatioGrayed->SetLineColor(kGray + 2);
  hRatioGrayed->Draw("same");

  hRatio->SetLineWidth(1);
  hRatio->Draw("same hist");
  // hRatioWide->Draw("same");

  c1->cd(1);
  TLatex* lat = new TLatex();
  lat->SetNDC();
  Float_t lPosText = 0.76;
  Float_t lYShift = 0.25;
  lat->SetTextSize(0.042);

  // save the glauber parameters explicitly
  TH1D* hGlauberParameters = new TH1D("hGlauberParameters", "", 10, 0, 10);
  TH1D* hGlauberFitRange = new TH1D("hGlauberFitRange", "", 10, 0, 10);

  // fitfunc
  hGlauberParameters->SetBinContent(1, fitfunc->GetParameter(0));
  hGlauberParameters->SetBinContent(2, fitfunc->GetParameter(1));
  hGlauberParameters->SetBinContent(3, fitfunc->GetParameter(2));
  hGlauberParameters->SetBinContent(4, fitfunc->GetParameter(3));
  hGlauberParameters->SetBinContent(5, fitfunc->GetParameter(4));
  hGlauberParameters->Write();

  Double_t lLoRangeGlauber, lHiRangeGlauber;
  fitfunc->GetRange(lLoRangeGlauber, lHiRangeGlauber);
  hGlauberFitRange->SetBinContent(1, lLoRangeGlauber);
  hGlauberFitRange->SetBinContent(2, lHiRangeGlauber);
  hGlauberFitRange->Write();

  hRatio->Write();
  fOutput->Write();

  return 0;
}
