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

/// \file evalMchTrackingEfficiency.cxx
/// \brief Macro to evaluate the MCH tracking efficiency using the results from taskMuonTrkEfficiency.cxx
///
/// \author Zaida Conesa del Valle <zaida.conesa.del.valle@cern.ch>
/// \author Andrea Tavira García <tavira-garcia@ijclab.in2p3.fr>
#include <vector>

#include "TFile.h"
#include "TDirectoryFile.h"
#include "THn.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

double computeEfficiencyPerChamber(THnF* hnf, int iAxis, int iCh, double binLimits[2]);
double computeEfficiencyPerChamber(THnF* hnf, const int iAxis[3], int iCh, double binLimits[3][2]);
double computeStationEfficiency(double effPerChamber[4], int iStation);
double computeTotalEfficiency(double effPerStation[5]);

void getBinLimits(THnF* hnf, int nBins, int iAxis, std::vector<double>& limits);
void setHistoMinPt(THnF* hnf, int iAxis, float minPt);
void setStyleCanvas(TCanvas* c);
void setHistoStylePerChamber(TH1F* h, int iCh);

//
// Main function to evaluate the muon tracking efficiency
//
void evalMchTrackingEfficiency(const char* inFile = "AnalysisResults.root", const char* outFileName = "evalMchTrkEff.root", const char* suffix = "", const float minPt = 0.5)
{
  // read ouput from the taskComputeMchEfficiency
  TFile* file = TFile::Open(inFile, "read");
  TDirectoryFile* dir = reinterpret_cast<TDirectoryFile*>(file->Get("task-muon-mch-trk-efficiency"));
  THnF* hHitsEtaPtPhi = reinterpret_cast<THnF*>(dir->Get("hHitsEtaPtPhi"));

  // retrieve the eta/pt/phi axis binning configuration
  const int etaAx = 1, ptAx = 2, phiAx = 3;
  const int listAx[3] = {etaAx, ptAx, phiAx};

  int nBinsEta = hHitsEtaPtPhi->GetAxis(etaAx)->GetNbins();
  int nBinsPt = hHitsEtaPtPhi->GetAxis(ptAx)->GetNbins();
  int nBinsPhi = hHitsEtaPtPhi->GetAxis(phiAx)->GetNbins();

  std::vector<double> limitsEta, limitsPt, limitsPhi;

  getBinLimits(hHitsEtaPtPhi, nBinsEta, etaAx, limitsEta); // these are values!
  getBinLimits(hHitsEtaPtPhi, nBinsPt, ptAx, limitsPt);    // these are values!
  getBinLimits(hHitsEtaPtPhi, nBinsPhi, phiAx, limitsPhi); // these are values!

  int nBinMinPt = hHitsEtaPtPhi->GetAxis(ptAx)->FindBin(minPt);

  // Define the output file
  TFile* outFile = new TFile(outFileName, "recreate");

  // define histograms efficiency per chamber
  //  vs eta, pt, phi
  int const kChamber = 10, kStation = 5;
  TH1F *hEffPerChamberEta[kChamber], *hEffPerChamberPt[kChamber], *hEffPerChamberPhi[kChamber];

  for (int iCh = 0; iCh < kChamber; iCh++) {
    hEffPerChamberEta[iCh] = new TH1F(Form("hEffPerChamberEta_%d", iCh), Form("hEff Chamber %d vs. #eta; #eta ; Efficiency", iCh), nBinsEta, limitsEta.data());
    hEffPerChamberPt[iCh] = new TH1F(Form("hEffPerChamberPt_%d", iCh), Form("hEff Chamber %d vs. #it{p}_{T}; #it{p}_{T} (GeV/#it{c}) ; Efficiency", iCh), nBinsPt, limitsPt.data());
    hEffPerChamberPhi[iCh] = new TH1F(Form("hEffPerChamberPhi_%d", iCh), Form("hEff Chamber %d vs. #varphi; #varphi (rad.) ; Efficiency", iCh), nBinsPhi, limitsPhi.data());
  }

  // define histogram for integrated efficiency per chamber
  // as well as the one per station
  // the total integrated value is stored in both histos as last quantity for reference
  TH1F* hEffIntegratedChamber = new TH1F("hEffIntegratedChamber", "integrated efficiency per chamber; ;Efficiency", kChamber + 1, -0.5, kChamber + 0.5);
  const char* hChNames[kChamber + 1];
  for (int i = 0; i < kChamber; i++) {
    hEffIntegratedChamber->GetXaxis()->SetBinLabel(i + 1, Form("Chamber %d", i));
  }
  hEffIntegratedChamber->GetXaxis()->SetBinLabel(kChamber + 1, "Total");

  // Define histogram for integrated efficiency per station
  TH1F* hEffIntegratedStation = new TH1F("hEffIntegratedStation", "integrated efficiency per station; ;Efficiency", kStation, -0.5, kStation + 0.5);
  const char* hStNames[kStation];
  for (int i = 0; i < 3; i++) {
    hEffIntegratedStation->GetXaxis()->SetBinLabel(i + 1, Form("Station %d", i));
  }
  hEffIntegratedStation->GetXaxis()->SetBinLabel(4, "Station 3 & 4");
  hEffIntegratedStation->GetXaxis()->SetBinLabel(kStation, "Total");

  // define also the eta-phi efficiency per chamber
  TH2F* hEffPerChamberEtaPhi[kChamber];
  for (int iCh = 0; iCh < kChamber; iCh++) {
    hEffPerChamberEtaPhi[iCh] = new TH2F(Form("hEffPerChamberEtaPhi_%d", iCh), Form("hEff Chamber %d vs. #eta vs. #varphi; #eta ; #varphi ; Efficiency", iCh), nBinsEta, limitsEta.data(), nBinsPhi, limitsPhi.data());
  }

  double effIntChamber[kChamber];
  double effIntStation[kStation];
  double effIntegrated = 0.;

  // Calculation of the efficiency per chamber for each variable
  for (int iCh = 0; iCh < 10; iCh++) {
    double binLimits[2] = {0., 0.};

    // Calculation vs. eta
    //   check min pt interval
    setHistoMinPt(hHitsEtaPtPhi, ptAx, minPt);
    for (int ikBin = 1; ikBin <= nBinsEta; ikBin++) {
      binLimits[0] = limitsEta[ikBin - 1]; // these are NOT bins
      binLimits[1] = limitsEta[ikBin];     // these are NOT bins
      double eff = computeEfficiencyPerChamber(hHitsEtaPtPhi, etaAx, iCh, binLimits);
      int hBin = hEffPerChamberEta[iCh]->FindBin(binLimits[0] + (binLimits[1] - binLimits[0]) / 2.);
      hEffPerChamberEta[iCh]->SetBinContent(hBin, eff);
    }
    // Calculation vs pt
    for (int ikBin = 1; ikBin <= nBinsPt; ikBin++) {
      binLimits[0] = limitsPt[ikBin - 1];
      binLimits[1] = limitsPt[ikBin];
      double eff = computeEfficiencyPerChamber(hHitsEtaPtPhi, ptAx, iCh, binLimits);
      int hBin = hEffPerChamberPt[iCh]->FindBin(binLimits[0] + (binLimits[1] - binLimits[0]) / 2.);
      hEffPerChamberPt[iCh]->SetBinContent(hBin, eff);
    }
    // Calculation vs phi
    //   check min pt interval
    setHistoMinPt(hHitsEtaPtPhi, ptAx, minPt);
    for (int ikBin = 1; ikBin <= nBinsPhi; ikBin++) {
      binLimits[0] = limitsPhi[ikBin - 1];
      binLimits[1] = limitsPhi[ikBin];
      double eff = computeEfficiencyPerChamber(hHitsEtaPtPhi, phiAx, iCh, binLimits);
      int hBin = hEffPerChamberPhi[iCh]->FindBin(binLimits[0] + (binLimits[1] - binLimits[0]) / 2.);
      hEffPerChamberPhi[iCh]->SetBinContent(hBin, eff);
    }
    // Calculation of the integrated quantity per chamber
    binLimits[0] = limitsEta[0];
    binLimits[1] = limitsEta[nBinsEta];
    effIntChamber[iCh] = computeEfficiencyPerChamber(hHitsEtaPtPhi, etaAx, iCh, binLimits);
    hEffIntegratedChamber->SetBinContent(iCh + 1, effIntChamber[iCh]);

  } // end calculation efficiency per chamber

  // Do integrated calculation per station
  double tmp[4] = {effIntChamber[0], effIntChamber[1], 0., 0.};
  effIntStation[0] = computeStationEfficiency(tmp, 0); // Station 1
  tmp[0] = effIntChamber[2];
  tmp[1] = effIntChamber[3];
  effIntStation[1] = computeStationEfficiency(tmp, 1); // Station 2
  tmp[0] = effIntChamber[4];
  tmp[1] = effIntChamber[5];
  effIntStation[2] = computeStationEfficiency(tmp, 2); // Station 3
  tmp[0] = effIntChamber[6];
  tmp[1] = effIntChamber[7];
  tmp[2] = effIntChamber[8];
  tmp[3] = effIntChamber[9];
  effIntStation[3] = computeStationEfficiency(tmp, 3); // Station 4 & 5

  for (int i = 1; i <= 4; i++) {
    hEffIntegratedStation->SetBinContent(i, effIntStation[i - 1]);
  }

  // Do integrated total calculation
  effIntegrated = computeTotalEfficiency(effIntStation);
  hEffIntegratedStation->SetBinContent(kStation, effIntegrated);
  hEffIntegratedChamber->SetBinContent(kChamber + 1, effIntegrated);

  // Calculation of the 2D eta-phi efficiency per chamber
  for (int iCh = 0; iCh < 10; iCh++) {
    double binLimits[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
    binLimits[ptAx - 1][0] = minPt;                                               // these are values!!
    binLimits[ptAx - 1][1] = hHitsEtaPtPhi->GetAxis(ptAx)->GetBinUpEdge(nBinsPt); // these are values!!
    // Loop over eta
    for (int ikBin = 1; ikBin <= nBinsEta; ikBin++) {
      binLimits[etaAx - 1][0] = limitsEta[ikBin - 1];
      binLimits[etaAx - 1][1] = limitsEta[ikBin];

      // Loop over phi
      for (int isbin = 1; isbin <= nBinsPhi; isbin++) {
        binLimits[phiAx - 1][0] = limitsPhi[isbin - 1];
        binLimits[phiAx - 1][1] = limitsPhi[isbin];
        double eff = computeEfficiencyPerChamber(hHitsEtaPtPhi, listAx, iCh, binLimits);
        int xBin = hEffPerChamberEtaPhi[iCh]->GetXaxis()->FindBin(binLimits[etaAx - 1][0] + (binLimits[etaAx - 1][1] - binLimits[etaAx - 1][0]) / 2.);
        int yBin = hEffPerChamberEtaPhi[iCh]->GetYaxis()->FindBin(binLimits[phiAx - 1][0] + (binLimits[phiAx - 1][1] - binLimits[phiAx - 1][0]) / 2.);
        hEffPerChamberEtaPhi[iCh]->SetBinContent(xBin, yBin, eff);
      }
    }
  }

  // Drawing output values
  double yPlotLimits[2] = {0, 1.1};
  gStyle->SetOptStat(0);
  TLegend* legEffChamber = new TLegend(0.6, 0.2, 0.8, 0.7);

  TCanvas* cEffChEta = new TCanvas(Form("cEffChEta%s", suffix), Form("MCH tracking efficiency per chamber vs. eta %s", suffix));
  setStyleCanvas(cEffChEta);
  hEffPerChamberEta[0]->GetYaxis()->SetRangeUser(yPlotLimits[0], yPlotLimits[1]);
  hEffPerChamberEta[0]->Draw();

  for (int iCh = 0; iCh < 10; iCh++) {
    setHistoStylePerChamber(hEffPerChamberEta[iCh], iCh);
    hEffPerChamberEta[iCh]->Draw("hsame");
    legEffChamber->AddEntry(hEffPerChamberEta[iCh], Form("Chamber %d", iCh), "l");
  }
  legEffChamber->Draw();

  TCanvas* cEffChPt = new TCanvas(Form("cEffChPt%s", suffix), Form("MCH tracking efficiency per chamber vs. pT %s", suffix));
  setStyleCanvas(cEffChPt);
  hEffPerChamberPt[0]->GetYaxis()->SetRangeUser(yPlotLimits[0], yPlotLimits[1]);
  hEffPerChamberPt[0]->Draw();
  for (int iCh = 0; iCh < 10; iCh++) {
    setHistoStylePerChamber(hEffPerChamberPt[iCh], iCh);
    hEffPerChamberPt[iCh]->Draw("hsame");
  }
  legEffChamber->Draw();

  TCanvas* cEffChPhi = new TCanvas(Form("cEffChPhi%s", suffix), Form("MCH tracking efficiency per chamber vs. phi %s", suffix));
  setStyleCanvas(cEffChPhi);
  hEffPerChamberPhi[0]->GetYaxis()->SetRangeUser(yPlotLimits[0], yPlotLimits[1]);
  hEffPerChamberPhi[0]->Draw();
  for (int iCh = 0; iCh < 10; iCh++) {
    setHistoStylePerChamber(hEffPerChamberPhi[iCh], iCh);
    hEffPerChamberPhi[iCh]->Draw("hsame");
  }
  legEffChamber->Draw();

  TCanvas* cEffIntegrated = new TCanvas(Form("cEffIntegrated%s", suffix), Form("MCH tracking integrated efficiency %s", suffix));
  setStyleCanvas(cEffIntegrated);
  setHistoStylePerChamber(hEffIntegratedChamber, 10);
  setHistoStylePerChamber(hEffIntegratedStation, 10);
  cEffIntegrated->Divide(1, 2);
  cEffIntegrated->cd(1);
  hEffIntegratedChamber->GetYaxis()->SetRangeUser(yPlotLimits[0], yPlotLimits[1]);
  hEffIntegratedChamber->SetMarkerSize(1.8);
  hEffIntegratedChamber->Draw("h,text");
  cEffIntegrated->cd(2);
  hEffIntegratedStation->SetMarkerSize(1.8);
  hEffIntegratedStation->Draw("h,text");

  TCanvas* cEffEtaPhi = new TCanvas(Form("cEffEtaPhi%s", suffix), Form("MCH tracking eta-phi efficiency %s", suffix));
  cEffEtaPhi->Divide(2, 5);
  for (int iCh = 0; iCh < 10; iCh++) {
    cEffEtaPhi->cd(iCh + 1);
    hEffPerChamberEtaPhi[iCh]->Draw("colz");
  }

  outFile->Write();
  file->Close();
}

// Function to retrieve the axis limits from the THn
void getBinLimits(THnF* hnf, int nBins, int iAxis, std::vector<double>& limits)
{
  TAxis* ax = hnf->GetAxis(iAxis);
  for (int i = 1; i <= nBins; i++) {
    limits.push_back(ax->GetBinLowEdge(i));
  }
  limits.push_back(ax->GetBinUpEdge(nBins));
  return;
}

// Evaluate the efficiency per Chamber
// Eff(i) = N(i+j) / ( N(i+j) + N(0+j) )
double computeEfficiencyPerChamber(THnF* hnf, int iAxis, int iCh, double binLimits[2])
{
  // Set range of study
  hnf->GetAxis(iAxis)->SetRangeUser(binLimits[0], binLimits[1]);

  // Project onto the Nhits axis
  TH1F* htmp = reinterpret_cast<TH1F*>(hnf->Projection(0));
  double NhitPairing[16];
  for (int i = 0; i < 16; i++) {
    NhitPairing[i] = 0;
    NhitPairing[i] = htmp->GetBinContent(i + 1);
  }
  double eff = 0.;
  if (iCh == 0 && (NhitPairing[0] + NhitPairing[2]) > 0.) { // Chamber 0 : St 1
    eff = NhitPairing[0] / (NhitPairing[0] + NhitPairing[2]);
  } else if (iCh == 1 && NhitPairing[0] > 0.) { // Chamber 1 : St 1
    eff = NhitPairing[0] / (NhitPairing[0] + NhitPairing[1]);
  } else if (iCh == 2 && NhitPairing[3] > 0.) { // Chamber 2 : St 2
    eff = NhitPairing[3] / (NhitPairing[3] + NhitPairing[5]);
  } else if (iCh == 3 && NhitPairing[3] > 0.) { // Chamber 3 : St 2
    eff = NhitPairing[3] / (NhitPairing[3] + NhitPairing[4]);
  } else if (iCh == 4 && NhitPairing[6] > 0.) { // Chamber 4 : St 3
    eff = NhitPairing[6] / (NhitPairing[6] + NhitPairing[8]);
  } else if (iCh == 5 && NhitPairing[6] > 0.) { // Chamber 5 : St 3
    eff = NhitPairing[6] / (NhitPairing[6] + NhitPairing[7]);
  } else if (iCh == 6 && NhitPairing[9] > 0.) { // Chamber 6 : St 4
    eff = NhitPairing[9] / (NhitPairing[9] + NhitPairing[11]);
  } else if (iCh == 7 && NhitPairing[9] > 0.) { // Chamber 7 : St 4
    eff = NhitPairing[9] / (NhitPairing[9] + NhitPairing[10]);
  } else if (iCh == 8 && NhitPairing[12] > 0.) { // Chamber 8 : St 5
    eff = NhitPairing[12] / (NhitPairing[12] + NhitPairing[14]);
  } else if (iCh == 9 && NhitPairing[12] > 0.) { // Chamber 9 : St 5
    eff = NhitPairing[12] / (NhitPairing[12] + NhitPairing[13]);
  }

  delete htmp;

  // reset to full range
  int nBins = hnf->GetAxis(iAxis)->GetNbins();
  hnf->GetAxis(iAxis)->SetRange(1, nBins);

  return eff;
}

// Evaluate the efficiency per Chamber
// Eff(i) = N(i+j) / ( N(i+j) + N(0+j) )
double computeEfficiencyPerChamber(THnF* hnf, const int iAxis[3], int iCh, double binLimits[3][2])
{
  // Set range of study
  for (int i = 0; i < 3; i++) {
    hnf->GetAxis(iAxis[i])->SetRangeUser(binLimits[i][0], binLimits[i][1]);
  }

  // Project onto the Nhits axis
  TH1F* htmp = reinterpret_cast<TH1F*>(hnf->Projection(0));
  double NhitPairing[16];
  for (int i = 0; i < 16; i++) {
    NhitPairing[i] = 0;
    NhitPairing[i] = htmp->GetBinContent(i + 1);
  }
  double eff = 0.;
  if (iCh == 0 && (NhitPairing[0] + NhitPairing[2]) > 0.) { // Chamber 0 : St 1
    eff = NhitPairing[0] / (NhitPairing[0] + NhitPairing[2]);
  } else if (iCh == 1 && NhitPairing[0] > 0.) { // Chamber 1 : St 1
    eff = NhitPairing[0] / (NhitPairing[0] + NhitPairing[1]);
  } else if (iCh == 2 && NhitPairing[3] > 0.) { // Chamber 2 : St 2
    eff = NhitPairing[3] / (NhitPairing[3] + NhitPairing[5]);
  } else if (iCh == 3 && NhitPairing[3] > 0.) { // Chamber 3 : St 2
    eff = NhitPairing[3] / (NhitPairing[3] + NhitPairing[4]);
  } else if (iCh == 4 && NhitPairing[6] > 0.) { // Chamber 4 : St 3
    eff = NhitPairing[6] / (NhitPairing[6] + NhitPairing[8]);
  } else if (iCh == 5 && NhitPairing[6] > 0.) { // Chamber 5 : St 3
    eff = NhitPairing[6] / (NhitPairing[6] + NhitPairing[7]);
  } else if (iCh == 6 && NhitPairing[9] > 0.) { // Chamber 6 : St 4
    eff = NhitPairing[9] / (NhitPairing[9] + NhitPairing[11]);
  } else if (iCh == 7 && NhitPairing[9] > 0.) { // Chamber 7 : St 4
    eff = NhitPairing[9] / (NhitPairing[9] + NhitPairing[10]);
  } else if (iCh == 8 && NhitPairing[12] > 0.) { // Chamber 8 : St 5
    eff = NhitPairing[12] / (NhitPairing[12] + NhitPairing[14]);
  } else if (iCh == 9 && NhitPairing[12] > 0.) { // Chamber 9 : St 5
    eff = NhitPairing[12] / (NhitPairing[12] + NhitPairing[13]);
  }

  delete htmp;

  // reset to full range
  for (int i = 0; i < 3; i++) {
    int nBins = hnf->GetAxis(iAxis[i])->GetNbins();
    hnf->GetAxis(iAxis[i])->SetRange(1, nBins);
  }

  return eff;
}

// Compute the efficiency per station
// Stations 1, 2, 3:               Eff = 1 - ( 1 - E(i) ) * ( 1 - (E(j) )
// Stations 4 & 5 together:   Eff = product(i=6...9) E(i) +  sum(i=6...9) [ ( 1 - E(i) )* product(j=6...9, i!=i) E(j) ]
double computeStationEfficiency(double effPerChamber[4], int iStation)
{
  double effSt = 0.;
  if (iStation < 3) { // Calculation for Station 1, 2, 3
    effSt = 1. - (1. - effPerChamber[0]) * (1. - effPerChamber[1]);
  } else { // Calculation for Station 4 & 5
    double product = 1., sum = 0.;
    for (int i = 0; i < 4; i++) {
      product *= effPerChamber[i];
    }
    for (int i = 0; i < 4; i++) {
      double tmpProduct = 1.;
      for (int j = 0; j < 4; j++) {
        if (i != j) {
          tmpProduct *= effPerChamber[j];
        }
      }
      sum += (1. - effPerChamber[i]) * tmpProduct;
    }
    effSt = product + sum;
  }
  return effSt;
}

// Total efficiency corresponds to the product of the efficiency per station
double computeTotalEfficiency(double effPerStation[5])
{
  double effTot = 1.;
  for (int i = 0; i < 4; i++) {
    effTot *= effPerStation[i];
  }
  return effTot;
}

// Set style for histos
void setHistoStylePerChamber(TH1F* h, int iCh)
{
  int iColor[11] = {kRed, kMagenta, kBlue, kAzure + 10, kGreen, kGreen + 3, kOrange + 7, kOrange + 4, kGray + 1, kGray + 3, kBlack};
  h->SetLineColor(iColor[iCh]);
  h->SetMarkerColor(iColor[iCh]);
}
// Set style for canvas
void setStyleCanvas(TCanvas* c)
{
  c->SetTickx();
  c->SetTicky();
}

void setHistoMinPt(THnF* hnf, int iAxis, float minPt)
{
  // Get the axis and its number of bins
  TAxis* axis = hnf->GetAxis(iAxis);
  int nBinsPt = axis->GetNbins();

  // Find the bin corresponding to minPt
  int nBinMinPt = axis->FindBin(minPt); // FindBin already gives the correct bin index

  // Set range using bin indices
  axis->SetRange(nBinMinPt, nBinsPt);
}
