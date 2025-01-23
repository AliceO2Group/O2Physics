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
/// \brief Class for D-h correlation extraction
/// \author Samuele Cattaruzzi <samuele.cattaruzzi@cern.ch>
/// \author Swapnesh Santosh Khade <swapnesh.santosh.khade@cern.ch>

#include "DhCorrelationExtraction.h"

#include <cstdio>
#include <iostream>

DhCorrelationExtraction::DhCorrelationExtraction() : // default constructor
                                                     fFileMass(0x0),
                                                     fFileSE(0x0),
                                                     fFileME(0x0),
                                                     fDirMass(0x0),
                                                     fDirSE(0x0),
                                                     fDirME(0x0),
                                                     fCorrectedCorrHisto(0x0),
                                                     fDmesonSpecies(kDsToKKPi),
                                                     fDmesonLabel("Ds"),
                                                     fNpools(9),
                                                     fDeltaEtaMin(-1.),
                                                     fDeltaEtaMax(1.),
                                                     fCorrectPoolsSeparately(kTRUE),
                                                     fSubtractSoftPiME(kFALSE),
                                                     fFileNameSE(""),
                                                     fFileNameME(""),
                                                     fDirNameSE(""),
                                                     fDirNameME(""),
                                                     fMassHistoNameSgn(""),
                                                     fMassHistoNameBkg(""),
                                                     fMassHistoNameSBs(""),
                                                     fSECorrelSignalRegionName(""),
                                                     fSECorrelSidebandsName(""),
                                                     fMECorrelSignalRegionName(""),
                                                     fMECorrelSidebandsName(""),
                                                     fBkgScaleFactor(1.),
                                                     fSgnYieldNorm(1.),
                                                     fRebin2Dhisto(kFALSE),
                                                     fRebinAxisDeltaEta(1),
                                                     fRebinAxisDeltaPhi(1),
                                                     fDebug(0)
{
}

DhCorrelationExtraction::DhCorrelationExtraction(const DhCorrelationExtraction& source) : // copy constructor
                                                                                          fFileMass(source.fFileMass),
                                                                                          fFileSE(source.fFileSE),
                                                                                          fFileME(source.fFileME),
                                                                                          fDirMass(source.fDirMass),
                                                                                          fDirSE(source.fDirSE),
                                                                                          fDirME(source.fDirME),
                                                                                          fCorrectedCorrHisto(source.fCorrectedCorrHisto),
                                                                                          fDmesonSpecies(source.fDmesonSpecies),
                                                                                          fDmesonLabel(source.fDmesonLabel),
                                                                                          fNpools(source.fNpools),
                                                                                          fDeltaEtaMin(source.fDeltaEtaMin),
                                                                                          fDeltaEtaMax(source.fDeltaEtaMax),
                                                                                          fCorrectPoolsSeparately(source.fCorrectPoolsSeparately),
                                                                                          fSubtractSoftPiME(source.fSubtractSoftPiME),
                                                                                          fFileNameSE(source.fFileNameSE),
                                                                                          fFileNameME(source.fFileNameME),
                                                                                          fDirNameSE(source.fDirNameSE),
                                                                                          fDirNameME(source.fDirNameME),
                                                                                          fMassHistoNameSgn(source.fMassHistoNameSgn),
                                                                                          fMassHistoNameBkg(source.fMassHistoNameBkg),
                                                                                          fMassHistoNameSBs(source.fMassHistoNameSBs),
                                                                                          fSECorrelSignalRegionName(source.fSECorrelSignalRegionName),
                                                                                          fSECorrelSidebandsName(source.fSECorrelSidebandsName),
                                                                                          fMECorrelSignalRegionName(source.fMECorrelSignalRegionName),
                                                                                          fMECorrelSidebandsName(source.fMECorrelSidebandsName),
                                                                                          fBkgScaleFactor(source.fBkgScaleFactor),
                                                                                          fSgnYieldNorm(source.fSgnYieldNorm),
                                                                                          fRebin2Dhisto(source.fRebin2Dhisto),
                                                                                          fRebinAxisDeltaEta(source.fRebinAxisDeltaEta),
                                                                                          fRebinAxisDeltaPhi(source.fRebinAxisDeltaPhi),
                                                                                          fDebug(source.fDebug)
{
}

DhCorrelationExtraction::~DhCorrelationExtraction()
// destructor
{
}

Bool_t DhCorrelationExtraction::SetDmesonSpecie(DmesonSpecie k)
{

  if (k < 0 || k > 3) {
    printf("[ERROR] D meson specie not correctly set!\n");
    return kFALSE;
  } else if (k == 0) {
    fDmesonLabel = "Dzero";
  } else if (k == 1) {
    fDmesonLabel = "Dplus";
  } else if (k == 2) {
    fDmesonLabel = "Ds";
  } else {
    fDmesonLabel = "Dstar";
  }

  fDmesonSpecies = k;
  return kTRUE;
}

Bool_t DhCorrelationExtraction::ExtractCorrelations(Double_t PtCandMin, Double_t PtCandMax, Double_t PtHadMin, Double_t PtHadMax, TString codeName)
{

  if (fSubtractSoftPiME) {
    printf("[INFO] Fake softPi subtraction in ME via extraction code is enabled!\n");
  }

  if (!fCorrectPoolsSeparately)
    fNpools = 1; // single histogram with integrated pools

  // Histograms definition
  TH2D* hSE_Sign[fNpools];
  TH2D* hME_Sign[fNpools];
  TH2D* hME_Sign_SoftPi[fNpools];
  TH2D* hSE_Sideb[fNpools];
  TH2D* hME_Sideb[fNpools];
  TH2D* hME_Sideb_SoftPi[fNpools];

  TH2D* hCorr_Sign[fNpools];
  TH2D* hCorr_Sideb[fNpools];

  TH2D* h2D_Sign;
  TH2D* h2D_Sideb;
  TH2D* h2D_Subtr;

  TH1D* h1D_Sign;
  TH1D* h1D_Sideb;
  TH1D* h1D_Subtr;
  TH1D* h1D_SubtrNorm;

  // if (fIntegratePtBins && iBinPtHad>0) continue;

  for (int iPool = 0; iPool < fNpools; iPool++) {

    // Retrieve 2D plots for SE and ME, signal and bkg regions, for each pTbin and pool
    hSE_Sign[iPool] = GetCorrelHisto(kSE, kSign, iPool, PtCandMin, PtCandMax, PtHadMin, PtHadMax);
    hME_Sign[iPool] = GetCorrelHisto(kME, kSign, iPool, PtCandMin, PtCandMax, PtHadMin, PtHadMax);
    hSE_Sideb[iPool] = GetCorrelHisto(kSE, kSideb, iPool, PtCandMin, PtCandMax, PtHadMin, PtHadMax);
    hME_Sideb[iPool] = GetCorrelHisto(kME, kSideb, iPool, PtCandMin, PtCandMax, PtHadMin, PtHadMax);

    hSE_Sign[iPool]->Sumw2();
    hME_Sign[iPool]->Sumw2();
    hSE_Sideb[iPool]->Sumw2();
    hME_Sideb[iPool]->Sumw2();

    // rebin axes deltaEta and deltaPhi
    if (fRebin2Dhisto) {
      hSE_Sign[iPool]->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi); // Xaxis: deltaEta, Yaxis: deltaPhi
      hSE_Sideb[iPool]->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi);
      hME_Sign[iPool]->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi);
      hME_Sideb[iPool]->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi);
      if (fSubtractSoftPiME) {
        hME_Sideb_SoftPi[iPool]->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi);
      }
    }

    if (fDebug >= 1) {
      TCanvas* c = new TCanvas(Form("cSEME_Original_%d_%1.1fto%1.1f", iPool, PtHadMin, PtHadMax), Form("cSEME_Original_%s_pool%d_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), iPool, PtCandMin, PtCandMax, PtHadMin, PtHadMax), 100, 100, 1600, 900);
      c->Divide(2, 2);
      c->cd(1);
      hSE_Sign[iPool]->SetMinimum(0);
      hSE_Sign[iPool]->Draw("lego2");
      c->cd(2);
      hME_Sign[iPool]->SetMinimum(0);
      hME_Sign[iPool]->Draw("lego2");
      c->cd(3);
      hSE_Sideb[iPool]->SetMinimum(0);
      hSE_Sideb[iPool]->Draw("lego2");
      c->cd(4);
      hME_Sideb[iPool]->SetMinimum(0);
      hME_Sideb[iPool]->Draw("lego2");
      c->SaveAs(Form("Output_CorrelationExtraction_%s_png/CorrSEandME_Original_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
      c->SaveAs(Form("Output_CorrelationExtraction_%s_Root/CorrSEandME_Original_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
    }

    // Scale bkg plots by ratio of signal region/sidebands
    hSE_Sideb[iPool]->Scale(fBkgScaleFactor);
    hME_Sideb[iPool]->Scale(fBkgScaleFactor);
    hSE_Sideb[iPool]->SetEntries(hSE_Sideb[iPool]->GetEntries() * fBkgScaleFactor);
    hME_Sideb[iPool]->SetEntries(hME_Sideb[iPool]->GetEntries() * fBkgScaleFactor);

    if (fSubtractSoftPiME) {
      hME_Sideb_SoftPi[iPool]->Scale(fBkgScaleFactor);
      hME_Sideb_SoftPi[iPool]->SetEntries(hME_Sideb_SoftPi[iPool]->GetEntries() * fBkgScaleFactor);
    }

    // Normalize ME plots for the entries in (deltaEta, deltaPhi) = (0, 0)
    NormalizeMEplot(hME_Sign[iPool], hME_Sign_SoftPi[iPool]);
    NormalizeMEplot(hME_Sideb[iPool], hME_Sideb_SoftPi[iPool]);

    // Apply Event Mixing Correction
    hCorr_Sign[iPool] = reinterpret_cast<TH2D*>(hSE_Sign[iPool]->Clone(Form("hCorr_Sign_Pool%d", iPool)));
    hCorr_Sign[iPool]->Sumw2();
    hCorr_Sign[iPool]->Divide(hME_Sign[iPool]);

    hCorr_Sideb[iPool] = reinterpret_cast<TH2D*>(hSE_Sideb[iPool]->Clone(Form("hCorr_Sideb_Pool%d", iPool)));
    hCorr_Sideb[iPool]->Sumw2();
    hCorr_Sideb[iPool]->Divide(hME_Sideb[iPool]);

    Double_t N_SEsign = 0, N_SEsideb = 0, N_sign = 0, N_sideb = 0;
    for (int i = 1; i <= hCorr_Sign[iPool]->GetXaxis()->GetNbins(); i++) {
      for (int j = 1; j <= hCorr_Sign[iPool]->GetYaxis()->GetNbins(); j++) {
        N_SEsign += hSE_Sign[iPool]->GetBinContent(i, j);
        N_SEsideb += hSE_Sideb[iPool]->GetBinContent(i, j);
        N_sign += hCorr_Sign[iPool]->GetBinContent(i, j);
        N_sideb += hCorr_Sideb[iPool]->GetBinContent(i, j);
      }
    }
    hSE_Sign[iPool]->SetEntries(N_SEsign);
    hSE_Sideb[iPool]->SetEntries(N_SEsideb);
    hCorr_Sign[iPool]->SetEntries(N_sign);
    hCorr_Sideb[iPool]->SetEntries(N_sideb);

    if (fDebug >= 1) {
      TCanvas* c = new TCanvas(Form("cSEME_%d_%1.1fto%1.1f", iPool, PtHadMin, PtHadMax), Form("cSEME_%s_pool%d_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), iPool, PtCandMin, PtCandMax, PtHadMin, PtHadMax), 100, 100, 1600, 900);
      c->Divide(3, 2);
      c->cd(1);
      hSE_Sign[iPool]->SetMinimum(0);
      hSE_Sign[iPool]->Draw("lego2");
      c->cd(2);
      hME_Sign[iPool]->SetMinimum(0);
      hME_Sign[iPool]->Draw("lego2");
      c->cd(3);
      hCorr_Sign[iPool]->SetMinimum(0);
      hCorr_Sign[iPool]->Draw("lego2");
      c->cd(4);
      hSE_Sideb[iPool]->SetMinimum(0);
      hSE_Sideb[iPool]->Draw("lego2");
      c->cd(5);
      hME_Sideb[iPool]->SetMinimum(0);
      hME_Sideb[iPool]->Draw("lego2");
      c->cd(6);
      hCorr_Sideb[iPool]->SetMinimum(0);
      hCorr_Sideb[iPool]->Draw("lego2");
      c->SaveAs(Form("Output_CorrelationExtraction_%s_png/CorrSEandME_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
      c->SaveAs(Form("Output_CorrelationExtraction_%s_Root/CorrSEandME_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
    }

    // Pools integration
    if (iPool == 0) {
      h2D_Sign = reinterpret_cast<TH2D*>(hCorr_Sign[0]->Clone("h2D_Sign"));
      h2D_Sideb = reinterpret_cast<TH2D*>(hCorr_Sideb[0]->Clone("h2D_Sideb"));
      h2D_Sign->Sumw2();
      h2D_Sideb->Sumw2();
    } else {
      h2D_Sign->Add(hCorr_Sign[iPool]);
      h2D_Sideb->Add(hCorr_Sideb[iPool]);
    }
  } // end pool loop

  // Draw 2D plots (Signal region and Sidebands)
  TCanvas* c2D = new TCanvas(Form("c2D_IntPools_PtHad%.0fto%.0f", PtHadMin, PtHadMax), Form("c2D_%s_IntPools_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), PtHadMin, PtHadMax), 100, 100, 1500, 800);
  SetTH2HistoStyle(h2D_Sign, Form("Signal region, %.0f < p^{%s}_{T} < %.0f GeV/c, %.0f < p^{assoc}_{T} < %.0f GeV/c", PtCandMin, fDmesonLabel.Data(), PtCandMax, PtHadMin, PtHadMax), "#Delta#eta", "#Delta#phi [rad]", "entries", 1.6, 1.6, 1.6, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04);
  SetTH2HistoStyle(h2D_Sideb, Form("Sideband region, %.0f < p^{%s}_{T} < %.0f GeV/c, %.0f < p^{assoc}_{T} < %.0f GeV/c", PtCandMin, fDmesonLabel.Data(), PtCandMax, PtHadMin, PtHadMax), "#Delta#eta", "#Delta#phi [rad]", "#frac{Y_{Bkg}}{Y_{SB}} entries", 1.6, 1.6, 1.6, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04);
  c2D->Divide(2, 1);
  c2D->cd(1);
  h2D_Sign->SetMinimum(0);
  h2D_Sign->Draw("lego2");
  c2D->cd(2);
  h2D_Sideb->SetMinimum(0);
  h2D_Sideb->Draw("lego2");
  c2D->SaveAs(Form("Output_CorrelationExtraction_%s_png/h2D_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
  c2D->SaveAs(Form("Output_CorrelationExtraction_%s_Root/h2D_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));

  // Bkg subtraction (2D plot)
  TCanvas* c2D_Sub = new TCanvas(Form("c2D_Subtr_IntPools_PtHAd%.0fto%.0f", PtHadMin, PtHadMax), Form("c2D_%s_Subtr_IntPools_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), PtHadMin, PtHadMax), 100, 100, 1500, 800);
  h2D_Subtr = reinterpret_cast<TH2D*>(h2D_Sign->Clone("h2D_Subtr"));
  h2D_Subtr->Sumw2();
  h2D_Subtr->Add(h2D_Sideb, -1);
  h2D_Subtr->SetEntries(h2D_Sign->GetEntries() - h2D_Sideb->GetEntries());
  h2D_Subtr->SetTitle("Signal region after sideb. subt. corr. - 2D");
  h2D_Subtr->Draw("lego2");
  c2D_Sub->SaveAs(Form("Output_CorrelationExtraction_%s_png/h2D_%s_Subtr_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
  c2D_Sub->SaveAs(Form("Output_CorrelationExtraction_%s_Root/h2D_%s_Subtr_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));

  // 1D projection
  h1D_Sign = reinterpret_cast<TH1D*>(h2D_Sign->ProjectionY("h1D_Sign")); // projection on deltaPhi axis
  h1D_Sideb = reinterpret_cast<TH1D*>(h2D_Sideb->ProjectionY("h1D_Sideb"));
  h1D_Sign->SetTitle("Signal region correlations");
  h1D_Sideb->SetTitle("Sidebands correlations");
  h1D_Sign->Scale(1. / h1D_Sign->GetXaxis()->GetBinWidth(1));
  h1D_Sideb->Scale(1. / h1D_Sideb->GetXaxis()->GetBinWidth(1));

  // Bkg subtraction (1D plot)
  h1D_Subtr = reinterpret_cast<TH1D*>(h1D_Sign->Clone("h1D_Subtr"));
  h1D_Subtr->Sumw2();
  h1D_Subtr->Add(h1D_Sideb, -1);
  h1D_Subtr->SetEntries(h1D_Sign->GetEntries() - h1D_Sideb->GetEntries());
  h1D_Subtr->SetTitle("Signal region after sideb. subt. corr.");

  // Draw 1D plots (Signal region, Sidebands, S-SB (subtr.))
  TCanvas* c1D = new TCanvas(Form("c1D_IntPools_%.0fto%.0f", PtHadMin, PtHadMax), Form("c1D_%s_IntPools_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), PtHadMin, PtHadMax), 100, 100, 1600, 500);
  c1D->Divide(3, 1);
  c1D->cd(1);
  h1D_Sign->Draw();
  c1D->cd(2);
  h1D_Sideb->Draw();
  c1D->cd(3);
  h1D_Subtr->Draw();
  c1D->SaveAs(Form("Output_CorrelationExtraction_%s_png/h1D_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
  c1D->SaveAs(Form("Output_CorrelationExtraction_%s_Root/h1D_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));

  // Apply normalization to number of triggers
  h1D_SubtrNorm = reinterpret_cast<TH1D*>(h1D_Subtr->Clone("h1D_SubtrNorm"));
  h1D_SubtrNorm->Sumw2();
  h1D_SubtrNorm->Scale(1. / fSgnYieldNorm);
  h1D_SubtrNorm->SetTitle("Signal region after sideb. subt. corr. - Normalized to # of triggers");

  fCorrectedCorrHisto = reinterpret_cast<TH1D*>(h1D_SubtrNorm->Clone(Form("hCorrectedCorr_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax)));

  // Draw 1D plots (Signal region, normalized)
  TCanvas* cFinal = new TCanvas(Form("cFinal_%.0fto%.0f", PtHadMin, PtHadMax), Form("cFinal_%s_IntPools_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), PtHadMin, PtHadMax), 100, 100, 1200, 700);
  h1D_SubtrNorm->Draw();
  cFinal->SaveAs(Form("Output_CorrelationExtraction_%s_png/AzimCorrDistr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
  cFinal->SaveAs(Form("Output_CorrelationExtraction_%s_Root/AzimCorrDistr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));

  return kTRUE;
}

Bool_t DhCorrelationExtraction::ReadInputSEandME()
{

  fFileSE = TFile::Open(fFileNameSE.Data());
  if (!fFileSE) {
    std::cout << "[ERROR] File " << fFileNameSE << " cannot be opened! check your file path!";
    return kFALSE;
  }

  fFileME = TFile::Open(fFileNameME.Data());
  if (!fFileME) {
    std::cout << "[ERROR] File " << fFileNameME << " cannot be opened! check your file path!";
    return kFALSE;
  }

  fDirSE = reinterpret_cast<TDirectoryFile*>(fFileSE->Get(fDirNameSE.Data()));
  fDirME = reinterpret_cast<TDirectoryFile*>(fFileME->Get(fDirNameME.Data()));

  std::cout << "===================== " << std::endl;
  std::cout << "Read inputs SE and ME" << std::endl;
  std::cout << "TFile SE    = " << fFileNameSE << std::endl;
  std::cout << "TFile ME    = " << fFileNameME << std::endl;
  std::cout << "TDir SE    = " << fDirNameSE << std::endl;
  std::cout << "TDir ME    = " << fDirNameME << std::endl;
  std::cout << "===================== " << std::endl;
  std::cout << " " << std::endl;

  return kTRUE;
}

Bool_t DhCorrelationExtraction::ReadInputInvMass()
{

  fFileMass = TFile::Open(fFileNameMass.Data());
  if (!fFileMass) {
    std::cout << "[ERROR] File " << fFileNameMass << " cannot be opened! check your file path!";
    return kFALSE;
  }

  std::cout << "===================== " << std::endl;
  std::cout << "Read inputs inv. mass" << std::endl;
  std::cout << "TFile Mass    = " << fFileNameMass << std::endl;
  std::cout << "Histo Sgn Yield    = " << fMassHistoNameSgn << std::endl;
  std::cout << "Histo Bkg Yield    = " << fMassHistoNameBkg << std::endl;
  std::cout << "Histo SBs Yield    = " << fMassHistoNameSBs << std::endl;
  std::cout << "===================== " << std::endl;
  std::cout << " " << std::endl;

  return kTRUE;
}

TH2D* DhCorrelationExtraction::GetCorrelHisto(Int_t SEorME, Int_t SorSB, Int_t pool, Double_t PtCandMin, Double_t PtCandMax, Double_t PtHadMin, Double_t PtHadMax)
{
  // TODO: Subtraction of softpion
  TH2D* h2D = new TH2D(); // pointer to be returned

  THnSparseD* hSparse = 0x0;

  if (SEorME == kSE) { // Same Event
    if (SorSB == kSign) {
      hSparse = reinterpret_cast<THnSparseD*>(fDirSE->Get(fSECorrelSignalRegionName.Data()));
    } else {
      hSparse = reinterpret_cast<THnSparseD*>(fDirSE->Get(fSECorrelSidebandsName.Data()));
    }
  } else { // Mixed Event
    if (SorSB == kSign) {
      hSparse = reinterpret_cast<THnSparseD*>(fDirME->Get(fMECorrelSignalRegionName.Data()));
    } else {
      hSparse = reinterpret_cast<THnSparseD*>(fDirME->Get(fMECorrelSidebandsName.Data()));
    }
  }
  Int_t binExtPtCandMin = (Int_t)hSparse->GetAxis(2)->FindBin(PtCandMin + 0.01); // axis2: ptCand, the 0.01 to avoid bin edges!
  Int_t binExtPtCandMax = (Int_t)hSparse->GetAxis(2)->FindBin(PtCandMax - 0.01);
  Int_t binExtPtHadMin = (Int_t)hSparse->GetAxis(3)->FindBin(PtHadMin + 0.01); // axis3: ptHad
  Int_t binExtPtHadMax = (Int_t)hSparse->GetAxis(3)->FindBin(PtHadMax - 0.01);
  Int_t binExtPoolMin;
  Int_t binExtPoolMax;
  if (fCorrectPoolsSeparately) {
    binExtPoolMin = (Int_t)hSparse->GetAxis(4)->FindBin(pool + 1.01); // axis4: pool bin
    binExtPoolMax = (Int_t)hSparse->GetAxis(4)->FindBin(pool + 1.99);
  } else { // merge all pools in one
    binExtPoolMin = 1;
    binExtPoolMax = (Int_t)hSparse->GetAxis(4)->GetNbins();
    // cout << "binExtPoolMax:" << binExtPoolMax <<endl;
  }
  // possibility to select a certain eta region
  Int_t binExtEtaMin = (Int_t)hSparse->GetAxis(1)->FindBin(fDeltaEtaMin + 0.0001);
  Int_t binExtEtaMax = (Int_t)hSparse->GetAxis(1)->FindBin(fDeltaEtaMax - 0.0001);
  if (binExtEtaMax > hSparse->GetAxis(1)->GetNbins())
    binExtEtaMax = hSparse->GetAxis(1)->GetNbins();
  if (binExtEtaMin < 1)
    binExtEtaMin = 1;

  hSparse->GetAxis(1)->SetRange(binExtEtaMin, binExtEtaMax);       // axis1: deltaEta
  hSparse->GetAxis(2)->SetRange(binExtPtCandMin, binExtPtCandMax); // axis2: ptCand
  hSparse->GetAxis(3)->SetRange(binExtPtHadMin, binExtPtHadMax);   // axis3: ptHad
  // hSparse -> GetAxis(4) -> SetRange(binExtPoolMin, binExtPoolMax); // axis4: pool bin

  h2D = reinterpret_cast<TH2D*>(hSparse->Projection(0, 1)); // axis0: deltaPhi, axis1: deltaEta
  if (SEorME == kSE) {                                      // Same Event
    if (SorSB == kSign) {
      h2D->SetName(Form("hCorr_SE_Sig_2D_PtCandBin%d_PtHadBin%d_iPool%d", binExtPtCandMin, binExtPtHadMin, pool));
    } else {
      h2D->SetName(Form("hCorr_SE_Sideb_2D_PtCandBin%d_PtHadBin%d_iPool%d", binExtPtCandMin, binExtPtHadMin, pool));
    }
  } else { // Mixed Event
    if (SorSB == kSign) {
      h2D->SetName(Form("hCorr_ME_Sig_2D_PtCandBin%d_PtHadBin%d_iPool%d", binExtPtCandMin, binExtPtHadMin, pool));
    } else {
      h2D->SetName(Form("hCorr_ME_Sideb_2D_PtCandBin%d_PtHadBin%d_iPool%d", binExtPtCandMin, binExtPtHadMin, pool));
    }
  }

  return h2D;
}

void DhCorrelationExtraction::GetSignalAndBackgroundForNorm(Double_t PtCandMin, Double_t PtCandMax)
{

  // using results obtained from HFInvariantMassFitter.cxx class
  TH1F* hMassFitSgnYield = reinterpret_cast<TH1F*>(fFileMass->Get(fMassHistoNameSgn.Data()));
  TH1F* hMassFitBkgYield = reinterpret_cast<TH1F*>(fFileMass->Get(fMassHistoNameBkg.Data()));
  TH1F* hMassFitSBsYield = reinterpret_cast<TH1F*>(fFileMass->Get(fMassHistoNameSBs.Data()));

  Int_t PtCandBin = hMassFitSgnYield->FindBin(PtCandMin + 0.01);
  if (PtCandBin != hMassFitSgnYield->FindBin(PtCandMax - 0.01))
    std::cout << "[ERROR] Pt bin in invariant mass histogram not univocally defined " << std::endl;

  Float_t SgnYield = hMassFitSgnYield->GetBinContent(PtCandBin);
  Float_t BkgYield = hMassFitBkgYield->GetBinContent(PtCandBin);
  Float_t SBsYield = hMassFitSBsYield->GetBinContent(PtCandBin);

  std::cout << "================================= " << std::endl;
  std::cout << "Getting invariant mass parameters " << std::endl;
  std::cout << "Pt cand " << PtCandMin << " - " << PtCandMax << std::endl;
  std::cout << "Signal yield    = " << SgnYield << std::endl;
  std::cout << "Bkg yield    = " << BkgYield << std::endl;
  std::cout << "Sideband yield    = " << SBsYield << std::endl;
  std::cout << "================================= " << std::endl;
  std::cout << " " << std::endl;

  SetSignalYieldforNorm(SgnYield);
  SetBkgScaleFactor(BkgYield / SBsYield);

  return;
}

void DhCorrelationExtraction::NormalizeMEplot(TH2D*& histoME, TH2D*& histoMEsoftPi)
{

  Double_t bin0phi = histoME->GetYaxis()->FindBin(0.);
  Double_t bin0eta = histoME->GetXaxis()->FindBin(0.);

  // evaluate the normalization (from ALL tracks, including possible fake softpions) -> **histoME indeed includes bin1+bin2 of THnSparse, i.e. all the tracks**
  Double_t factorNorm = 0;
  for (int in = -1; in <= 0; in++) {
    factorNorm += histoME->GetBinContent(bin0eta, bin0phi + in);
  }
  for (int in = -1; in <= 0; in++) {
    factorNorm += histoME->GetBinContent(bin0eta - 1, bin0phi + in);
  }
  factorNorm /= 4.;

  if (fSubtractSoftPiME) {
    histoME->Add(histoMEsoftPi, -1); // remove the tracks compatible with soft pion (if requested)
  }

  // apply the normalization
  histoME->Scale(1. / factorNorm);

  return;
}

void DhCorrelationExtraction::SetTH1HistoStyle(TH1D*& histo, TString hTitle, TString hXaxisTitle, TString hYaxisTitle,
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

void DhCorrelationExtraction::SetTH2HistoStyle(TH2D*& histo, TString hTitle, TString hXaxisTitle, TString hYaxisTitle, TString hZaxisTitle,
                                               Float_t hTitleXaxisOffset, Float_t hTitleYaxisOffset, Float_t hTitleZaxisOffset,
                                               Float_t hTitleXaxisSize, Float_t hTitleYaxisSize, Float_t hTitleZaxisSize,
                                               Float_t hLabelXaxisSize, Float_t hLabelYaxisSize, Float_t hLabelZaxisSize,
                                               Bool_t centerXaxisTitle, Bool_t centerYaxisTitle)
{

  histo->SetTitle(hTitle.Data());
  histo->GetXaxis()->SetTitle(hXaxisTitle.Data());
  histo->GetYaxis()->SetTitle(hYaxisTitle.Data());
  histo->GetZaxis()->SetTitle(hZaxisTitle.Data());
  histo->GetXaxis()->SetTitleOffset(hTitleXaxisOffset);
  histo->GetYaxis()->SetTitleOffset(hTitleYaxisOffset);
  histo->GetZaxis()->SetTitleOffset(hTitleZaxisOffset);
  histo->GetXaxis()->SetTitleSize(hTitleXaxisSize);
  histo->GetYaxis()->SetTitleSize(hTitleYaxisSize);
  histo->GetZaxis()->SetTitleSize(hTitleZaxisSize);
  histo->GetXaxis()->SetLabelSize(hLabelXaxisSize);
  histo->GetYaxis()->SetLabelSize(hLabelYaxisSize);
  histo->GetZaxis()->SetLabelSize(hLabelZaxisSize);
  histo->GetXaxis()->CenterTitle(centerXaxisTitle);
  histo->GetYaxis()->CenterTitle(centerYaxisTitle);

  return;
}
