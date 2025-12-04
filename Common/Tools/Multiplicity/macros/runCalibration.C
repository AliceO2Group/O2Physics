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

#include <iostream>

/// @brief function to calibrate centrality
/// @param lInputFileName name of input file.
/// @param anchorPointPercentage anchor point percentage to use
/// @param matchRange width of region in which data/glauber matching is to be done in rolling anchoring test
/// @param doNpartNcoll wether or not to attempt calculating Npart, Ncoll in centrality bins
void runCalibration(TString lInputFileName = "results/AR_544122_glauberNBD_ancestorMode2_hFT0C_BCs.root", double anchorPointPercentage = 90.0, double matchRange = 200.0, bool doNpartNcoll = false)
{
  TFile* file = new TFile(lInputFileName.Data(), "READ");
  file->ls();

  TH1F* hData = (TH1F*)file->Get("hV0MUltraFine");
  TH1F* hGlauberParameters = (TH1F*)file->Get("hGlauberParameters");
  TH1F* hGlauberFitRange = (TH1F*)file->Get("hGlauberFitRange");
  hData->SetName("hData");
  TH1F* hStitched = (TH1F*)hData->Clone("hStitched");
  TH1F* hFit = (TH1F*)file->Get("hGlauber");

  TCanvas* c1 = new TCanvas("c1", "", 800, 600);
  c1->SetLeftMargin(0.17);
  c1->SetBottomMargin(0.17);
  c1->SetRightMargin(0.15);
  c1->SetTopMargin(0.05);
  c1->SetTicks(1, 1);
  c1->SetLogz();
  c1->SetFrameFillStyle(0);
  c1->SetFillStyle(0);

  cout << "Data bin width: " << hData->GetBinWidth(1) << endl;
  cout << "Fit bin width: " << hFit->GetBinWidth(1) << endl;
  cout << "Match range to use: " << matchRange << endl;

  //____________________________________________
  double anchorPointFraction = anchorPointPercentage / 100.f;
  double anchorPoint = -1; // the anchor point value in raw

  //____________________________________________
  // doing partial integration up to certain point for finding anchor point bin
  for (int ii = 1; ii < hData->GetNbinsX() + 1; ii++) {
    // renormalize data curve
    int bin1 = ii + 1;
    int bin2 = hData->FindBin(hData->GetBinLowEdge(ii + 1) + matchRange + 1e-3);
    double matchRangeData = hData->Integral(bin1, bin2);
    double matchRangeFit = hFit->Integral(bin1, bin2);

    // rescale fit to match in the vicinity of the region we're at
    hFit->Scale(matchRangeData / matchRangeFit);

    double integralFit = hFit->Integral(1, ii);
    double integralData = hData->Integral(ii + 1, hData->GetNbinsX() + 1);
    double integralAll = integralFit + integralData;

    cout << "at bin #" << ii << ", integrated up to " << hData->GetBinLowEdge(ii + 1) << " fraction above this value is: " << integralData / integralAll << endl;
    anchorPoint = hData->GetBinLowEdge(ii + 1);

    if (integralData / integralAll < anchorPointFraction)
      break;
  }

  //____________________________________________
  for (int ii = 1; ii < hData->GetNbinsX() + 1; ii++) {
    // renormalize data curve
    if (hData->GetBinCenter(ii) < anchorPoint)
      hStitched->SetBinContent(ii, hFit->GetBinContent(ii));
  }

  cout << "Anchor point determined to be: " << anchorPoint << endl;
  cout << "Preparing stitched histogram ... " << endl;

  hFit->SetLineColor(kRed);
  hStitched->SetLineColor(kBlue);

  hData->GetYaxis()->SetTitleSize(0.055);
  hData->GetXaxis()->SetTitleSize(0.055);
  hData->GetYaxis()->SetLabelSize(0.04);
  hData->GetXaxis()->SetLabelSize(0.04);
  hData->SetTitle("");
  hData->Draw("hist");
  hFit->Draw("hist same");
  hStitched->Draw("hist same");

  // All fine, let's try the calibrator
  multCalibrator* lCalib = new multCalibrator("lCalib");
  lCalib->SetAnchorPointPercentage(100.0f);
  lCalib->SetAnchorPointRaw(-1e-6);

  // Set standard Pb-Pb boundaries
  lCalib->SetStandardOnePercentBoundaries();

  TString calibFileName = lInputFileName.Data();
  calibFileName.ReplaceAll("glauberNBD", "calibration");
  TFile* fileCalib = new TFile(calibFileName.Data(), "RECREATE");

  TH1F* hCalib = lCalib->GetCalibrationHistogram(hStitched, "hCalib");

  TCanvas* c2 = new TCanvas("c2", "", 800, 600);
  c2->SetLeftMargin(0.17);
  c2->SetBottomMargin(0.17);
  c2->SetRightMargin(0.15);
  c2->SetTopMargin(0.05);
  c2->SetTicks(1, 1);
  //  c2->SetLogz();
  c2->SetFrameFillStyle(0);
  c2->SetFillStyle(0);

  hCalib->GetYaxis()->SetTitleSize(0.055);
  hCalib->GetXaxis()->SetTitleSize(0.055);
  hCalib->GetYaxis()->SetLabelSize(0.04);
  hCalib->GetXaxis()->SetLabelSize(0.04);
  hCalib->SetTitle("");
  hCalib->Draw();

  fileCalib->cd();

  hData->Write();
  hCalib->Write();
  hStitched->Write();
  hFit->Write();

  if (doNpartNcoll) {
    cout << "Will now attempt to calculate % -> Np, Nc map..." << endl;

    TProfile* hProfileNpart = new TProfile("hProfileNpart", "", 100, 0, 100);
    TProfile* hProfileNcoll = new TProfile("hProfileNcoll", "", 100, 0, 100);
    TH2F* h2dNpart = new TH2F("h2dNpart", "", 100, 0, 100, 500, -0.5f, 499.5f);
    TH2F* h2dNcoll = new TH2F("h2dNcoll", "", 100, 0, 100, 3000, -0.5f, 2999.5);

    // Replay
    multGlauberNBDFitter* g = new multGlauberNBDFitter("lglau");
    TF1* fitfunc = g->GetGlauberNBD();

    // Step 1: open the (Npart, Ncoll) pair information, provide
    TFile* fbasefile = new TFile("basehistos.root", "READ");
    TH2D* hNpNc = (TH2D*)fbasefile->Get("hNpNc");
    g->SetNpartNcollCorrelation(hNpNc);
    g->InitializeNpNc();

    fitfunc->SetParameter(0, hGlauberParameters->GetBinContent(1));
    fitfunc->SetParameter(1, hGlauberParameters->GetBinContent(2));
    fitfunc->SetParameter(2, hGlauberParameters->GetBinContent(3));
    fitfunc->SetParameter(3, hGlauberParameters->GetBinContent(4));
    fitfunc->SetParameter(4, hGlauberParameters->GetBinContent(5));

    Double_t lMax = hData->GetBinLowEdge(hData->GetNbinsX() + 1);

    // uncomment if Np Nc needed -> Warning, slow!
    g->CalculateAvNpNc(hProfileNpart, hProfileNcoll, h2dNpart, h2dNcoll, hCalib, 0, lMax);

    hProfileNpart->Write();
    hProfileNcoll->Write();
    h2dNpart->Write();
    h2dNcoll->Write();
  }

  fileCalib->Write();
}
