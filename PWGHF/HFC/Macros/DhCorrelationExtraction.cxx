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

#include "DhCorrelationExtraction.h"

#include <TCanvas.h>
#include <TDirectoryFile.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TString.h>

#include <RtypesCore.h>

#include <cstdio>
#include <iostream>

DhCorrelationExtraction::DhCorrelationExtraction() : // default constructor
                                                     fFileMass(0x0),
                                                     fFileSE(0x0),
                                                     fFileME(0x0),
                                                     fFileFDTemplate(0x0),
                                                     fFileFDPromptFrac(0x0),
                                                     fFileSecPart(0x0),
                                                     fDirMass(0x0),
                                                     fDirSE(0x0),
                                                     fDirME(0x0),
                                                     fDirSecPart(0x0),
                                                     fFilePromptMc(0x0),
                                                     fFileNonPromptMc(0x0),
                                                     fCorrectedCorrHisto(0x0),
                                                     fCorrectedCorrHisto_BaselineSubtr(0x0),
                                                     fCorrectedCorrHisto_Reflected(0x0),
                                                     fCorrectedCorrHisto_Reflected_BaselineSubtr(0x0),
                                                     fDmesonSpecies(kDsToKKPi),
                                                     fDmesonLabel("Ds"),
                                                     fNpools(9),
                                                     fDeltaEtaMin(-1.),
                                                     fDeltaEtaMax(1.),
                                                     fCorrectPoolsSeparately(kTRUE),
                                                     fSubtractSoftPiME(kFALSE),
                                                     fFileNameSE(""),
                                                     fFileNameME(""),
                                                     fFileSecPartName(""),
                                                     fFilePromptMcRecName(""),
                                                     fFileNonPromptMcRecName(""),
                                                     fDirNameSE(""),
                                                     fDirNameME(""),
                                                     fDirSecPartName(""),
                                                     fMassHistoNameSgn(""),
                                                     fMassHistoNameBkg(""),
                                                     fMassHistoNameSBs(""),
                                                     fSECorrelSignalRegionName(""),
                                                     fSECorrelSidebandsName(""),
                                                     fSECorrelSidebandLeftName(""),
                                                     fSECorrelSidebandRightName(""),
                                                     fMECorrelSignalRegionName(""),
                                                     fMECorrelSidebandsName(""),
                                                     fMECorrelSidebandLeftName(""),
                                                     fMECorrelSidebandRightName(""),
                                                     fFileFDTemplateName(""),
                                                     fFileFDPromptFracName(""),
                                                     fHistoFDTemplatePromptName(""),
                                                     fHistoFDTemplateNonPromptName(""),
                                                     fHistoFDPromptFracName(""),
                                                     fHistoPrimaryPartName(""),
                                                     fHistoAllPartName(""),
                                                     fBkgScaleFactor(1.),
                                                     fSgnYieldNorm(1.),
                                                     fBkgYield(1.),
                                                     fRebinAngCorr(kFALSE),
                                                     fRebinFDCorr(kFALSE),
                                                     fRebinSecPart(kFALSE),
                                                     fSidebandDivided(kFALSE),
                                                     fUseSidebLeft(kFALSE),
                                                     fUseSidebRight(kFALSE),
                                                     fRebinAxisDeltaEta(1),
                                                     fRebinAxisDeltaPhi(1),
                                                     fBinPtCand(0),
                                                     fBinPtHad(0),
                                                     fDebug(0),
                                                     fFDsubtraction(0),
                                                     fSecPartContamination(0),
                                                     fCorrBiasBtoD(0)
{
}

DhCorrelationExtraction::DhCorrelationExtraction(const DhCorrelationExtraction& source) : // copy constructor
                                                                                          fFileMass(source.fFileMass),
                                                                                          fFileSE(source.fFileSE),
                                                                                          fFileME(source.fFileME),
                                                                                          fFileFDTemplate(source.fFileFDTemplate),
                                                                                          fFileFDPromptFrac(source.fFileFDPromptFrac),
                                                                                          fFileSecPart(source.fFileSecPart),
                                                                                          fFilePromptMc(source.fFilePromptMc),
                                                                                          fFileNonPromptMc(source.fFileNonPromptMc),
                                                                                          fDirMass(source.fDirMass),
                                                                                          fDirSecPart(source.fDirSecPart),
                                                                                          fCorrectedCorrHisto(source.fCorrectedCorrHisto),
                                                                                          fCorrectedCorrHisto_BaselineSubtr(source.fCorrectedCorrHisto_BaselineSubtr),
                                                                                          fCorrectedCorrHisto_Reflected(source.fCorrectedCorrHisto_Reflected),
                                                                                          fCorrectedCorrHisto_Reflected_BaselineSubtr(source.fCorrectedCorrHisto_Reflected_BaselineSubtr),
                                                                                          fDmesonSpecies(source.fDmesonSpecies),
                                                                                          fDmesonLabel(source.fDmesonLabel),
                                                                                          fNpools(source.fNpools),
                                                                                          fDeltaEtaMin(source.fDeltaEtaMin),
                                                                                          fDeltaEtaMax(source.fDeltaEtaMax),
                                                                                          fCorrectPoolsSeparately(source.fCorrectPoolsSeparately),
                                                                                          fSubtractSoftPiME(source.fSubtractSoftPiME),
                                                                                          fFileNameSE(source.fFileNameSE),
                                                                                          fFileNameME(source.fFileNameME),
                                                                                          fFileSecPartName(source.fFileSecPartName),
                                                                                          fFilePromptMcRecName(source.fFilePromptMcRecName),
                                                                                          fFileNonPromptMcRecName(source.fFileNonPromptMcRecName),
                                                                                          fDirNameSE(source.fDirNameSE),
                                                                                          fDirNameME(source.fDirNameME),
                                                                                          fDirSecPartName(source.fDirSecPartName),
                                                                                          fMassHistoNameSgn(source.fMassHistoNameSgn),
                                                                                          fMassHistoNameBkg(source.fMassHistoNameBkg),
                                                                                          fMassHistoNameSBs(source.fMassHistoNameSBs),
                                                                                          fSECorrelSignalRegionName(source.fSECorrelSignalRegionName),
                                                                                          fSECorrelSidebandsName(source.fSECorrelSidebandsName),
                                                                                          fSECorrelSidebandLeftName(source.fSECorrelSidebandLeftName),
                                                                                          fSECorrelSidebandRightName(source.fSECorrelSidebandRightName),
                                                                                          fMECorrelSignalRegionName(source.fMECorrelSignalRegionName),
                                                                                          fMECorrelSidebandsName(source.fMECorrelSidebandsName),
                                                                                          fMECorrelSidebandLeftName(source.fMECorrelSidebandLeftName),
                                                                                          fMECorrelSidebandRightName(source.fMECorrelSidebandRightName),
                                                                                          fFileFDTemplateName(source.fFileFDTemplateName),
                                                                                          fFileFDPromptFracName(source.fFileFDPromptFracName),
                                                                                          fHistoFDTemplatePromptName(source.fHistoFDTemplatePromptName),
                                                                                          fHistoFDTemplateNonPromptName(source.fHistoFDTemplateNonPromptName),
                                                                                          fHistoFDPromptFracName(source.fHistoFDPromptFracName),
                                                                                          fHistoPrimaryPartName(source.fHistoPrimaryPartName),
                                                                                          fHistoAllPartName(source.fHistoAllPartName),
                                                                                          fBkgScaleFactor(source.fBkgScaleFactor),
                                                                                          fSgnYieldNorm(source.fSgnYieldNorm),
                                                                                          fBkgYield(source.fBkgYield),
                                                                                          fRebinAngCorr(source.fRebinAngCorr),
                                                                                          fRebinFDCorr(source.fRebinFDCorr),
                                                                                          fRebinSecPart(source.fRebinSecPart),
                                                                                          fSidebandDivided(source.fSidebandDivided),
                                                                                          fUseSidebLeft(source.fUseSidebLeft),
                                                                                          fUseSidebRight(source.fUseSidebRight),
                                                                                          fRebinAxisDeltaEta(source.fRebinAxisDeltaEta),
                                                                                          fRebinAxisDeltaPhi(source.fRebinAxisDeltaPhi),
                                                                                          fBinPtCand(source.fBinPtCand),
                                                                                          fBinPtHad(source.fBinPtHad),
                                                                                          fDebug(source.fDebug),
                                                                                          fFDsubtraction(source.fFDsubtraction),
                                                                                          fSecPartContamination(source.fSecPartContamination),
                                                                                          fCorrBiasBtoD(source.fCorrBiasBtoD)
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

  TH2D* h2D_FDTemplatePrompt;
  TH2D* h2D_FDTemplateNonPrompt;

  TH1D* h1D_Sign;
  TH1D* h1D_Sideb;
  TH1D* h1D_Subtr;
  TH1D* h1D_SignNorm;
  TH1D* h1D_SidebNorm;
  TH1D* h1D_SubtrNorm;
  TH1D* h1D_FDTemplatePrompt;
  TH1D* h1D_FDTemplateNonPrompt;
  TH1D* h1D_TemplateTotal;
  TH1D* h1D_SubtrFDNorm;
  TH1D* h1D_PrimaryPartCorr;
  TH1D* h1D_AllPartCorr;
  TH1D* h1D_SecPartFrac;
  TH1D* h1D_SubtrNorm_SecPart;
  TH1D* h1D_BaselineSubtr;
  TH1D* h1D_ReflCorr;
  TH1D* h1D_ReflCorr_BaselineSubtr;
  TH1D* hModul;
  TH1D* hBeforeModulCorr;

  Double_t FDPromptFrac;

  // if (fIntegratePtBins && iBinPtHad>0) continue;

  for (int iPool = 0; iPool < fNpools; iPool++) {
    // Retrieve 2D plots for SE and ME, signal and bkg regions, for each pTbin and pool
    hSE_Sign[iPool] = GetCorrelHisto(kSE, kSign, iPool, PtCandMin, PtCandMax, PtHadMin, PtHadMax);
    std::cout << "Got SE histogram signal region" << std::endl;
    hME_Sign[iPool] = GetCorrelHisto(kME, kSign, iPool, PtCandMin, PtCandMax, PtHadMin, PtHadMax);
    std::cout << "Got ME histogram signal region" << std::endl;
    hSE_Sideb[iPool] = GetCorrelHisto(kSE, kSideb, iPool, PtCandMin, PtCandMax, PtHadMin, PtHadMax);
    std::cout << "Got SE histogram sdeband region" << std::endl;
    hME_Sideb[iPool] = GetCorrelHisto(kME, kSideb, iPool, PtCandMin, PtCandMax, PtHadMin, PtHadMax);
    std::cout << "Got ME histogram sdeband region" << std::endl;

    hSE_Sign[iPool]->Sumw2();
    hME_Sign[iPool]->Sumw2();
    hSE_Sideb[iPool]->Sumw2();
    hME_Sideb[iPool]->Sumw2();

    // rebin axes deltaEta and deltaPhi
    if (fRebinAngCorr) {
      hSE_Sign[iPool]->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi); // Xaxis: deltaEta, Yaxis: deltaPhi
      hSE_Sideb[iPool]->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi);
      hME_Sign[iPool]->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi);
      hME_Sideb[iPool]->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi);
      if (fSubtractSoftPiME) {
        hME_Sideb_SoftPi[iPool]->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi);
      }
      std::cout << "SE and ME histograms rebinned" << std::endl;
    }

    if (fDebug >= 1) {
      TCanvas* c = new TCanvas(Form("cSE_Original_%d_%1.1fto%1.1f", iPool, PtHadMin, PtHadMax), Form("cSE_Original_%s_pool%d_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), iPool, PtCandMin, PtCandMax, PtHadMin, PtHadMax), 100, 100, 1600, 900);
      c->Divide(2, 1);
      c->cd(1);
      hSE_Sign[iPool]->SetMinimum(0);
      hSE_Sign[iPool]->Draw("lego2");
      c->cd(2);
      hSE_Sideb[iPool]->SetMinimum(0);
      hSE_Sideb[iPool]->Draw("lego2");
      c->SaveAs(Form("Output_CorrelationExtraction_%s_png/CorrSE_Original_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
      c->SaveAs(Form("Output_CorrelationExtraction_%s_Root/CorrSE_Original_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
    }

    if (fDebug >= 1) {
      TCanvas* c = new TCanvas(Form("cME_Original_%d_%1.1fto%1.1f", iPool, PtHadMin, PtHadMax), Form("cME_Original_%s_pool%d_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), iPool, PtCandMin, PtCandMax, PtHadMin, PtHadMax), 100, 100, 1600, 900);
      c->Divide(2, 1);
      c->cd(1);
      hME_Sign[iPool]->SetMinimum(0);
      hME_Sign[iPool]->Draw("lego2");
      c->cd(2);
      hME_Sideb[iPool]->SetMinimum(0);
      hME_Sideb[iPool]->Draw("lego2");
      c->SaveAs(Form("Output_CorrelationExtraction_%s_png/CorrME_Original_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
      c->SaveAs(Form("Output_CorrelationExtraction_%s_Root/CorrME_Original_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
    }

    // Scale bkg plots by ratio of signal region/sidebands
    hSE_Sideb[iPool]->Scale(fBkgScaleFactor);
    hME_Sideb[iPool]->Scale(fBkgScaleFactor); // when normalised this factor should cancel out
    std::cout << "[INFO] fBkgScaleFactor    = " << fBkgScaleFactor << std::endl;
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
      c->SaveAs(Form("Output_CorrelationExtraction_%s_png/CorrSEandME_%s_Canvas_PtCand%.0fto%.0f_Pool%d_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, iPool, PtHadMin, PtHadMax));
      c->SaveAs(Form("Output_CorrelationExtraction_%s_Root/CorrSEandME_%s_Canvas_PtCand%.0fto%.0f_Pool%d_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, iPool, PtHadMin, PtHadMax));
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

  // Get FD correlations for FD subtraction
  if (fFDsubtraction) {
    h2D_FDTemplatePrompt = GetFDTemplateHisto(kPrompt, PtCandMin, PtCandMax, PtHadMin, PtHadMax);
    h2D_FDTemplateNonPrompt = GetFDTemplateHisto(kFD, PtCandMin, PtCandMax, PtHadMin, PtHadMax);
    // h1D_BaselineSubtr
    FDPromptFrac = GetFDPromptFrac(PtCandMin, PtCandMax, PtHadMin, PtHadMax);

    if (fRebinFDCorr) {
      h2D_FDTemplatePrompt->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi);
      h2D_FDTemplateNonPrompt->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi);
    }

    if (fDebug >= 1) {
      TCanvas* c = new TCanvas(Form("cFDTemplate_PtCand%.0fto%.0f_PtHad%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax), Form("cFDTemplate_%s_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
      c->Divide(2, 1);
      c->cd(1);
      h2D_FDTemplatePrompt->SetMinimum(0);
      h2D_FDTemplatePrompt->Draw("lego2");
      c->cd(2);
      h2D_FDTemplateNonPrompt->SetMinimum(0);
      h2D_FDTemplateNonPrompt->Draw("lego2");
      c->SaveAs(Form("Output_CorrelationExtraction_%s_png/CorrFDTemplate_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
    }
  }

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
  TCanvas* c1D = new TCanvas(Form("c1D_IntPools_PtCand%.0fto%.0f_PtHad%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax), Form("c1D_%s_IntPools_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), PtHadMin, PtHadMax), 100, 100, 1600, 500);
  c1D->Divide(3, 1);
  c1D->cd(1);
  h1D_Sign->Draw();
  c1D->cd(2);
  h1D_Sideb->Draw();
  c1D->cd(3);
  h1D_Subtr->Draw();
  c1D->SaveAs(Form("Output_CorrelationExtraction_%s_png/h1D_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
  c1D->SaveAs(Form("Output_CorrelationExtraction_%s_Root/h1D_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));

  if (fDebug >= 1) {
    h1D_SignNorm = reinterpret_cast<TH1D*>(h1D_Sign->Clone("h1D_Sign_Norm"));
    h1D_SidebNorm = reinterpret_cast<TH1D*>(h1D_Sideb->Clone("h1D_Sideb_Norm"));
    h1D_SignNorm->Scale(1. / (fSgnYieldNorm + fBkgYield));
    // h1D_SidebNorm -> Scale(1./fBkgYield);
    h1D_SidebNorm->Scale(1. / fBkgScaleFactor);
    h1D_SidebNorm->Scale(1. / fSBYield);
    h1D_SignNorm->SetMarkerStyle(kFullCircle);
    h1D_SignNorm->SetMarkerSize(1.2);
    h1D_SignNorm->SetLineColor(kRed);
    h1D_SignNorm->SetMarkerColor(kRed);
    h1D_SignNorm->SetLineWidth(2);
    h1D_SidebNorm->SetMinimum(0);
    h1D_SidebNorm->SetMarkerStyle(kFullSquare);
    h1D_SidebNorm->SetMarkerSize(1.2);
    h1D_SidebNorm->SetLineColor(kBlue);
    h1D_SidebNorm->SetMarkerColor(kBlue);
    h1D_SidebNorm->SetLineWidth(2);
    h1D_SidebNorm->SetTitle(Form("%.0f < p_{T} < %.0f", PtCandMin, PtCandMax));
    TCanvas* c = new TCanvas(Form("c_IntPools_PtCand%.0fto%.0f_PtHad%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax), "");
    c->cd();
    h1D_SidebNorm->Draw();
    h1D_SignNorm->Draw("same");
    c->SaveAs(Form("Output_CorrelationExtraction_%s_png/ComparisonSignalSidebCorr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
    c->SaveAs(Form("Output_CorrelationExtraction_%s_Root/ComparisonSignalSidebCorr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
  }
  // Apply normalization to number of triggers
  h1D_SubtrNorm = reinterpret_cast<TH1D*>(h1D_Subtr->Clone("h1D_SubtrNorm"));
  h1D_SubtrNorm->Sumw2();
  h1D_SubtrNorm->Scale(1. / fSgnYieldNorm);
  h1D_SubtrNorm->SetTitle("Signal region after sideb. subt. corr. - Normalized to # of triggers");

  // Correction for bias B to D topologies
  if (fCorrBiasBtoD) {
    hModul = EvaluateMCClosModulations(PtCandMin, PtCandMax, PtHadMin, PtHadMax);
    TCanvas* c1D_corrBbias = new TCanvas(Form("c1D_corrBbias_IntPools_PtCand%.0fto%.0f_PtHad%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax), Form("c1D_corrBbias_%s_IntPools_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), PtHadMin, PtHadMax), 100, 100, 1600, 500);
    c1D_corrBbias->cd();
    hBeforeModulCorr = reinterpret_cast<TH1D*>(h1D_SubtrNorm->Clone("hBeforeModulCorr"));
    hBeforeModulCorr->SetLineColor(kViolet - 3);
    hBeforeModulCorr->GetYaxis()->SetRangeUser(0., 5.);
    hBeforeModulCorr->Draw();
    h1D_SubtrNorm->Multiply(hModul);
    h1D_SubtrNorm->Draw("same");
    c1D_corrBbias->SaveAs(Form("Output_CorrelationExtraction_%s_png/ComparisonCorrBiasBtoD_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
    c1D_corrBbias->SaveAs(Form("Output_CorrelationExtraction_%s_Root/ComparisonCorrBiasBtoD_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));

    TFile* file = new TFile(Form("Output_CorrelationExtraction_%s_Root/SystematicCorrBiasBtoD_%s_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax), "RECREATE"); // Open file in write mode
    TH1D* h1D_SubtrNorm_Clone = reinterpret_cast<TH1D*>(h1D_SubtrNorm->Clone("h1D_SubtrNorm_Clone"));
    h1D_SubtrNorm_Clone = ReflectCorrHistogram(h1D_SubtrNorm_Clone);
    hBeforeModulCorr = ReflectCorrHistogram(hBeforeModulCorr);
    TH1D* hSystematicCorrBiasBtoD = reinterpret_cast<TH1D*>(h1D_SubtrNorm_Clone->Clone("hSystematicCorrBiasBtoD"));
    hSystematicCorrBiasBtoD->Add(h1D_SubtrNorm_Clone, hBeforeModulCorr, 1, -1);
    // Set bin contents to absolute values
    for (int i = 1; i <= hSystematicCorrBiasBtoD->GetNbinsX(); ++i) {
      hSystematicCorrBiasBtoD->SetBinContent(i, std::abs(hSystematicCorrBiasBtoD->GetBinContent(i)) / TMath::Sqrt(12));
      hSystematicCorrBiasBtoD->SetBinError(i, 0.);
    }
    hSystematicCorrBiasBtoD->Write();
    file->Close();
  }

  // Secondary particle contamination
  if (fSecPartContamination) {
    h1D_PrimaryPartCorr = GetCorrelHistoSecondaryPart(kPrimaryPart, PtCandMin, PtCandMax, PtHadMin, PtHadMax);
    h1D_AllPartCorr = GetCorrelHistoSecondaryPart(kAllPart, PtCandMin, PtCandMax, PtHadMin, PtHadMax);
    h1D_PrimaryPartCorr->Sumw2();
    h1D_AllPartCorr->Sumw2();
    if (fRebinSecPart) {
      h1D_PrimaryPartCorr->RebinX(fRebinAxisDeltaPhi); // Xaxis: deltaPhi
      h1D_AllPartCorr->RebinX(fRebinAxisDeltaPhi);
      std::cout << "Secondary particle histogram rebinned" << std::endl;
    }
    h1D_SecPartFrac = reinterpret_cast<TH1D*>(h1D_PrimaryPartCorr->Clone(Form("hCorrRatio_PtD%.0fto%.0f_PtHad%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax)));
    h1D_SecPartFrac->Sumw2();
    h1D_SecPartFrac->Divide(h1D_PrimaryPartCorr, h1D_AllPartCorr, 1., 1., "B");

    TCanvas* c1D = new TCanvas(Form("c1D_CorrPrimaryPart_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax), Form("c1D_%s_CorrPrimaryPart_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), PtHadMin, PtHadMax));
    c1D->cd();
    SetTH1HistoStyle(h1D_SecPartFrac, Form("%.0f < p_{T} < %.0f GeV/c", PtCandMin, PtCandMax), "#Delta#phi [rad]", "#frac{primary part.}{part. selected}");
    h1D_SecPartFrac->Draw();
    c1D->SaveAs(Form("Output_CorrelationExtraction_%s_png/CorrPrimaryPartRatio_%s_Canvas_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
    c1D->SaveAs(Form("Output_CorrelationExtraction_%s_Root/CorrPrimaryPartRatio_%s_Canvas_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));

    h1D_SubtrNorm_SecPart = reinterpret_cast<TH1D*>(h1D_SubtrNorm->Clone("h1D_SubtrNorm_SecPart"));
    h1D_SubtrNorm_SecPart->Sumw2();
    Int_t nBinsPhi = h1D_SubtrNorm_SecPart->GetNbinsX();
    if (nBinsPhi != h1D_SecPartFrac->GetNbinsX()) {
      std::cout << "[ERROR]: nBinsPhi different between h1D_SubtrNorm and h1D_SecPartFrac" << std::endl;
      return kFALSE;
    }
    h1D_SubtrNorm_SecPart->Multiply(h1D_SecPartFrac);
  }

  // FD Subtraction
  if (fFDsubtraction) {
    h1D_FDTemplatePrompt = reinterpret_cast<TH1D*>(h2D_FDTemplatePrompt->ProjectionY("h1D_FDTemplatePrompt"));
    h1D_FDTemplateNonPrompt = reinterpret_cast<TH1D*>(h2D_FDTemplateNonPrompt->ProjectionY("h1D_FDTemplateNonPrompt"));

    h1D_FDTemplatePrompt->Scale(1. / h1D_FDTemplatePrompt->GetXaxis()->GetBinWidth(1));
    h1D_FDTemplateNonPrompt->Scale(1. / h1D_FDTemplateNonPrompt->GetXaxis()->GetBinWidth(1));

    h1D_TemplateTotal = reinterpret_cast<TH1D*>(h1D_FDTemplatePrompt->Clone("h1D_TemplateTotal"));
    h1D_TemplateTotal->Sumw2();
    h1D_TemplateTotal->Scale(FDPromptFrac);
    h1D_TemplateTotal->Add(h1D_FDTemplateNonPrompt, 1 - FDPromptFrac);

    if (fDebug >= 1) {
      TCanvas* c = new TCanvas(Form("cFDTemplate_1D_%1.1fto%1.1f", PtHadMin, PtHadMax), Form("cFDTemplate_%s_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
      c->cd();
      h1D_TemplateTotal->SetMinimum(0);
      h1D_FDTemplateNonPrompt->SetMinimum(0);
      h1D_TemplateTotal->SetMarkerColor(kGreen);
      h1D_TemplateTotal->SetLineColor(kGreen);
      h1D_TemplateTotal->SetLineWidth(2);
      h1D_TemplateTotal->SetMarkerStyle(kFullCircle);
      h1D_FDTemplatePrompt->SetMarkerColor(kRed);
      h1D_FDTemplatePrompt->SetLineColor(kRed);
      h1D_FDTemplatePrompt->SetLineWidth(2);
      h1D_FDTemplatePrompt->SetMarkerStyle(kFullCircle);
      h1D_FDTemplateNonPrompt->SetMarkerColor(kBlue);
      h1D_FDTemplateNonPrompt->SetLineColor(kBlue);
      h1D_FDTemplateNonPrompt->SetLineWidth(2);
      h1D_FDTemplateNonPrompt->SetMarkerStyle(kFullCircle);
      h1D_FDTemplateNonPrompt->Draw();
      h1D_FDTemplatePrompt->Draw("same");
      h1D_TemplateTotal->Draw("same");
      TLegend* lFD = new TLegend();
      lFD->AddEntry(h1D_TemplateTotal, "Total template");
      lFD->AddEntry(h1D_FDTemplatePrompt, "Prompt Template");
      lFD->AddEntry(h1D_FDTemplateNonPrompt, "Non prompt template");
      lFD->Draw("same");
      c->SaveAs(Form("Output_CorrelationExtraction_%s_png/CorrFDTemplate_1D_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
      c->SaveAs(Form("Output_CorrelationExtraction_%s_Root/CorrFDTemplate_1D_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
    }

    Double_t BaselineFD = CalculateBaseline(h1D_TemplateTotal, kTRUE);
    Double_t BaselineData;
    if (fSecPartContamination) {
      BaselineData = CalculateBaseline(h1D_SubtrNorm_SecPart, kTRUE);
    } else {
      BaselineData = CalculateBaseline(h1D_SubtrNorm, kTRUE);
    }

    std::cout << "===================== " << std::endl;
    std::cout << "Baseline FD: " << BaselineFD << std::endl;
    std::cout << "Baseline Data: " << BaselineData << std::endl;
    std::cout << "===================== " << std::endl;
    std::cout << " " << std::endl;

    Double_t Baselinediff = BaselineData - BaselineFD;
    TH1D* hBaselineDiff = reinterpret_cast<TH1D*>(h1D_FDTemplateNonPrompt->Clone("hBaselineDiff"));
    for (int iBin = 0; iBin < hBaselineDiff->GetNbinsX(); iBin++) {
      hBaselineDiff->SetBinContent(iBin + 1, Baselinediff);
    }
    h1D_FDTemplateNonPrompt->Add(hBaselineDiff);
    h1D_TemplateTotal->Add(hBaselineDiff);
    if (fSecPartContamination) {
      h1D_SubtrFDNorm = reinterpret_cast<TH1D*>(h1D_SubtrNorm_SecPart->Clone("h1D_SubtrFDNorm"));
    } else {
      h1D_SubtrFDNorm = reinterpret_cast<TH1D*>(h1D_SubtrNorm->Clone("h1D_SubtrFDNorm"));
    }
    h1D_FDTemplateNonPrompt->Scale(1 - FDPromptFrac);
    h1D_SubtrFDNorm->Add(h1D_FDTemplateNonPrompt, -1);
    h1D_SubtrFDNorm->Scale(1. / FDPromptFrac);

    if (fDebug >= 1) {
      TCanvas* c1 = new TCanvas(Form("cFDTemplateSubtr_%1.1fto%1.1f", PtHadMin, PtHadMax), Form("cFDTemplateSubtr_%s_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax), 100, 100, 1600, 900);
      c1->cd();
      h1D_SubtrNorm->SetLineColor(kRed);
      h1D_SubtrNorm_SecPart->SetLineColor(kOrange);
      h1D_FDTemplateNonPrompt->SetLineColor(kBlue);
      h1D_SubtrFDNorm->SetLineColor(kGreen);
      h1D_TemplateTotal->SetLineColor(kMagenta);
      h1D_SubtrNorm->SetMinimum(0);
      h1D_SubtrNorm_SecPart->SetMinimum(0);
      h1D_FDTemplateNonPrompt->SetMinimum(0);
      h1D_SubtrFDNorm->SetMinimum(0);
      // h1D_SubtrNorm -> GetYaxis() -> SetRangeUser(0., 8.);
      h1D_SubtrNorm->SetMarkerStyle(kFullCircle);
      h1D_SubtrNorm->SetMarkerSize(1.2);
      h1D_SubtrNorm->SetMarkerColor(kRed);
      h1D_SubtrNorm->SetLineWidth(2);
      h1D_SubtrNorm_SecPart->SetMarkerStyle(kFullCircle);
      h1D_SubtrNorm_SecPart->SetMarkerSize(1.2);
      h1D_SubtrNorm_SecPart->SetMarkerColor(kOrange);
      h1D_SubtrNorm_SecPart->SetLineWidth(2);
      h1D_SubtrFDNorm->SetMarkerStyle(kFullCircle);
      h1D_SubtrFDNorm->SetMarkerSize(1.2);
      h1D_SubtrFDNorm->SetMarkerColor(kGreen);
      h1D_SubtrFDNorm->SetLineWidth(2);
      h1D_SubtrNorm->GetYaxis()->SetTitle("#frac{1}{N_{D}} #frac{dN^{assoc. part}}{d#Delta#phi}");
      h1D_SubtrNorm_SecPart->GetYaxis()->SetTitle("#frac{1}{N_{D}} #frac{dN^{assoc. part}}{d#Delta#phi}");
      if (fSecPartContamination) {
        h1D_SubtrNorm_SecPart->Draw();
      } else {
        h1D_SubtrNorm->Draw();
      }
      // h1D_FDTemplateNonPrompt -> Draw("same");
      h1D_SubtrFDNorm->Draw("same");
      h1D_TemplateTotal->Draw("same");
      c1->SaveAs(Form("Output_CorrelationExtraction_%s_png/CorrFDTemplateSubtr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
      c1->SaveAs(Form("Output_CorrelationExtraction_%s_Root/CorrFDTemplateSubtr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
    }
  }

  if (fFDsubtraction) {
    fCorrectedCorrHisto = reinterpret_cast<TH1D*>(h1D_SubtrFDNorm->Clone(Form("hCorrectedCorr_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax)));
  } else if (fSecPartContamination) {
    fCorrectedCorrHisto = reinterpret_cast<TH1D*>(h1D_SubtrNorm_SecPart->Clone(Form("hCorrectedCorr_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax)));
  } else {
    fCorrectedCorrHisto = reinterpret_cast<TH1D*>(h1D_SubtrNorm->Clone(Form("hCorrectedCorr_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax)));
  }

  std::cout << "Analysis steps completed - baseline subtraction missing" << std::endl;

  // Draw 1D plots (Signal region, normalized)
  TCanvas* cFinal = new TCanvas(Form("cFinal_%.0fto%.0f", PtHadMin, PtHadMax), Form("cFinal_%s_IntPools_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), PtHadMin, PtHadMax), 100, 100, 1200, 700);
  h1D_SubtrNorm->SetLineColor(kBlue + 1);
  h1D_SubtrNorm->SetMarkerColor(kBlue + 1);
  h1D_SubtrNorm->SetMarkerStyle(kFullCircle);
  h1D_SubtrNorm->SetMinimum(0);
  h1D_SubtrNorm->Draw();
  if (fSecPartContamination) {
    h1D_SubtrNorm_SecPart->SetLineColor(kRed + 1);
    h1D_SubtrNorm_SecPart->SetMarkerColor(kRed + 1);
    h1D_SubtrNorm_SecPart->SetMarkerStyle(kFullCircle);
    h1D_SubtrNorm_SecPart->Draw("same");
  }
  if (fFDsubtraction) {
    h1D_SubtrFDNorm->SetLineColor(kGreen + 2);
    h1D_SubtrFDNorm->SetMarkerColor(kGreen + 2);
    h1D_SubtrFDNorm->SetMarkerStyle(kFullCircle);
    h1D_SubtrFDNorm->Draw("same");
  }
  if (fFDsubtraction)
    h1D_TemplateTotal->Draw("same");
  TLegend* lFinal = new TLegend();
  lFinal->AddEntry(h1D_SubtrNorm, "Corr. after bkg subtr.");
  if (fFDsubtraction)
    lFinal->AddEntry(h1D_TemplateTotal, "CR Mode 2 total template");
  if (fSecPartContamination) {
    lFinal->AddEntry(h1D_SubtrNorm_SecPart, "Corr. after sec. part. correction");
  }
  if (fFDsubtraction) {
    lFinal->AddEntry(h1D_SubtrFDNorm, "Corr. FD subtr.");
  }
  lFinal->Draw("same");
  cFinal->SaveAs(Form("Output_CorrelationExtraction_%s_png/AzimCorrDistr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
  cFinal->SaveAs(Form("Output_CorrelationExtraction_%s_Root/AzimCorrDistr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));

  // Baseline subtraction
  Double_t BaselineData, BaselineDataErr;
  TH1D* hBaseline = reinterpret_cast<TH1D*>(h1D_SubtrNorm->Clone("hBaseline"));
  hBaseline->Sumw2();
  if (fFDsubtraction) {
    BaselineData = CalculateBaseline(h1D_SubtrFDNorm, kTRUE, kFALSE); // introduced kFALSE
    BaselineDataErr = CalculateBaselineError(h1D_SubtrFDNorm, kTRUE, kFALSE);
    for (int iBin = 0; iBin < hBaseline->GetNbinsX(); iBin++) {
      hBaseline->SetBinContent(iBin + 1, BaselineData);
      hBaseline->SetBinError(iBin + 1, BaselineDataErr);
    }
    h1D_BaselineSubtr = reinterpret_cast<TH1D*>(h1D_SubtrFDNorm->Clone("h1D_BaselineSubtr"));
    h1D_BaselineSubtr->Add(hBaseline, -1.);
  } else if (fSecPartContamination) {
    BaselineData = CalculateBaseline(h1D_SubtrNorm_SecPart, kTRUE, kFALSE);
    BaselineDataErr = CalculateBaselineError(h1D_SubtrNorm_SecPart, kTRUE, kFALSE);
    for (int iBin = 0; iBin < hBaseline->GetNbinsX(); iBin++) {
      hBaseline->SetBinContent(iBin + 1, BaselineData);
      hBaseline->SetBinError(iBin + 1, BaselineDataErr);
    }
    h1D_BaselineSubtr = reinterpret_cast<TH1D*>(h1D_SubtrNorm_SecPart->Clone("h1D_BaselineSubtr"));
    h1D_BaselineSubtr->Add(hBaseline, -1.);
  } else {
    BaselineData = CalculateBaseline(h1D_SubtrNorm, kTRUE, kFALSE);
    BaselineDataErr = CalculateBaselineError(h1D_SubtrNorm, kTRUE, kFALSE);
    for (int iBin = 0; iBin < hBaseline->GetNbinsX(); iBin++) {
      hBaseline->SetBinContent(iBin + 1, BaselineData);
      hBaseline->SetBinError(iBin + 1, BaselineDataErr);
    }
    h1D_BaselineSubtr = reinterpret_cast<TH1D*>(h1D_SubtrNorm->Clone("h1D_BaselineSubtr"));
    h1D_BaselineSubtr->Add(hBaseline, -1.);
  }

  fCorrectedCorrHisto_BaselineSubtr = reinterpret_cast<TH1D*>(h1D_BaselineSubtr->Clone(Form("hCorrectedCorrBaselineSubtr_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax)));

  TCanvas* cFinal_BaselineSubtr = new TCanvas(Form("cFinal_BaselineSubtr_%.0fto%.0f", PtHadMin, PtHadMax), Form("cFinal_BaselineSubtr_%s_IntPools_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), PtHadMin, PtHadMax), 100, 100, 1200, 700);
  h1D_BaselineSubtr->SetMarkerColor(kOrange + 8);
  h1D_BaselineSubtr->SetLineColor(kOrange + 8);
  h1D_BaselineSubtr->GetYaxis()->SetRangeUser(-0.2, 8.);
  h1D_BaselineSubtr->Draw();
  if (fFDsubtraction) {
    h1D_SubtrFDNorm->Draw("same");
  } else if (fSecPartContamination) {
    h1D_SubtrNorm_SecPart->Draw("same");
  } else {
    h1D_SubtrNorm->Draw("same");
  }
  hBaseline->SetMarkerColor(kPink - 6);
  hBaseline->SetMarkerStyle(kFullSquare);
  hBaseline->SetLineColor(kPink - 6);
  hBaseline->Draw("same");

  cFinal_BaselineSubtr->SaveAs(Form("Output_CorrelationExtraction_%s_png/AzimCorrDistr_BaselineSubtr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
  cFinal_BaselineSubtr->SaveAs(Form("Output_CorrelationExtraction_%s_Root/AzimCorrDistr_BaselineSubtr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));

  // Reflected histograms
  if (fFDsubtraction) {
    h1D_ReflCorr = ReflectCorrHistogram(h1D_SubtrFDNorm);
  } else if (fSecPartContamination) {
    h1D_ReflCorr = ReflectCorrHistogram(h1D_SubtrNorm_SecPart);
  } else {
    h1D_ReflCorr = ReflectCorrHistogram(h1D_SubtrNorm);
  }

  /* used as control using Run2 reflection function
  if (fFDsubtraction) {
    h1D_ReflCorr = ReflectHistoRun2(h1D_SubtrFDNorm, 0.5);
  } else if (fSecPartContamination) {
    h1D_ReflCorr = ReflectHistoRun2(h1D_SubtrNorm_SecPart, 0.5);
  } else {
    h1D_ReflCorr = ReflectHistoRun2(h1D_SubtrNorm, 0.5);
  }*/

  TCanvas* cFinal_Reflected = new TCanvas(Form("cFinal_Reflected_%.0fto%.0f", PtHadMin, PtHadMax), Form("cFinal_Reflected_%s_IntPools_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), PtHadMin, PtHadMax), 100, 100, 1200, 700);
  cFinal_Reflected->cd();
  SetTH1HistoStyle(h1D_ReflCorr, Form("%.0f < p_{T} < %.0f GeV/c", PtCandMin, PtCandMax), "#Delta#phi [rad]", "#frac{1}{N_{D}}#frac{dN^{assoc}}{d#Delta#phi} [rad^{-1}]", kFullCircle, kOrange + 8, 1.6, kOrange + 8, 3);
  h1D_ReflCorr->SetMinimum(0);
  h1D_ReflCorr->Draw();
  cFinal_Reflected->SaveAs(Form("Output_CorrelationExtraction_%s_png/AzimCorrDistr_Reflected_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
  cFinal_Reflected->SaveAs(Form("Output_CorrelationExtraction_%s_Root/AzimCorrDistr_Reflected_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));

  // Reflected histograms baseline subtracted
  TH1D* hBaseline_Refl = reinterpret_cast<TH1D*>(h1D_ReflCorr->Clone("hBaseline_Refl"));
  hBaseline_Refl->Sumw2();
  BaselineData = CalculateBaseline(h1D_ReflCorr, kFALSE, kTRUE);
  BaselineDataErr = CalculateBaselineError(h1D_ReflCorr, kFALSE, kTRUE);

  for (int iBin = 0; iBin < hBaseline_Refl->GetNbinsX(); iBin++) {
    hBaseline_Refl->SetBinContent(iBin + 1, BaselineData);
    hBaseline_Refl->SetBinError(iBin + 1, BaselineDataErr);
  }
  h1D_ReflCorr_BaselineSubtr = reinterpret_cast<TH1D*>(h1D_ReflCorr->Clone("h1D_ReflCorr_BaselineSubtr"));
  h1D_ReflCorr_BaselineSubtr->Sumw2();
  h1D_ReflCorr_BaselineSubtr->Add(hBaseline_Refl, -1.);

  TF1* fConstZero = new TF1("fConstZero", "[0]", 0., TMath::Pi());
  fConstZero->SetParameter(0, 0.);
  fConstZero->SetLineColor(kMagenta);
  fConstZero->SetLineStyle(9);
  fConstZero->SetLineWidth(4);
  fConstZero->SetTitle("");

  TCanvas* cFinal_Reflected_BaselineSubtr = new TCanvas(Form("cFinal_Reflected_BaselineSubtr_%.0fto%.0f", PtHadMin, PtHadMax), Form("cFinal_Reflected_BaselineSubtr_%s_IntPools_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), PtHadMin, PtHadMax), 100, 100, 1200, 700);
  SetTH1HistoStyle(h1D_ReflCorr_BaselineSubtr, Form("%.0f < p_{T} < %.0f GeV/c", PtCandMin, PtCandMax), "#Delta#phi [rad]", "#frac{1}{N_{D}}#frac{dN^{assoc}}{d#Delta#phi} [rad^{-1}]", kFullCircle, kRed + 1, 1.6, kRed + 1, 3);
  hBaseline_Refl->SetMarkerColor(kOrange);
  hBaseline_Refl->SetMarkerStyle(kFullSquare);
  hBaseline_Refl->SetLineColor(kOrange);
  cFinal_Reflected_BaselineSubtr->cd();
  h1D_ReflCorr->SetMinimum(-0.8);
  h1D_ReflCorr->SetStats(0);
  hBaseline_Refl->SetStats(0);
  h1D_ReflCorr->Draw();
  hBaseline_Refl->Draw("same");
  h1D_ReflCorr_BaselineSubtr->SetStats(0);
  h1D_ReflCorr_BaselineSubtr->Draw("same"); // then keep just this
  fConstZero->Draw("same");
  cFinal_Reflected_BaselineSubtr->SaveAs(Form("Output_CorrelationExtraction_%s_png/AzimCorrDistr_Reflected_BaselineSubtr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));
  cFinal_Reflected_BaselineSubtr->SaveAs(Form("Output_CorrelationExtraction_%s_Root/AzimCorrDistr_Reflected_BaselineSubtr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax));

  fCorrectedCorrHisto_Reflected = reinterpret_cast<TH1D*>(h1D_ReflCorr->Clone(Form("hCorrectedCorrReflected_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax)));
  fCorrectedCorrHisto_Reflected_BaselineSubtr = reinterpret_cast<TH1D*>(h1D_ReflCorr_BaselineSubtr->Clone(Form("hCorrectedCorrReflected_BaselineSubtr_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax)));

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

Bool_t DhCorrelationExtraction::ReadInputFDSubtr()
{

  fFileFDTemplate = TFile::Open(fFileFDTemplateName.Data());
  fFileFDPromptFrac = TFile::Open(fFileFDPromptFracName.Data());
  if (!fFileFDTemplate) {
    std::cout << "[ERROR] File " << fFileFDTemplateName << " cannot be opened! check your file path!" << std::endl;
    return kFALSE;
  }
  if (!fFileFDPromptFrac) {
    std::cout << "[ERROR] File " << fFileFDPromptFracName << " cannot be opened! check your file path!" << std::endl;
    return kFALSE;
  }

  std::cout << "===================== " << std::endl;
  std::cout << "Read inputs FD template" << std::endl;
  std::cout << "TFile FD template    = " << fFileFDTemplateName << std::endl;
  std::cout << "TFile FD Prompt Frac    = " << fFileFDPromptFracName << std::endl;
  std::cout << "Histo FD template Prompt    = " << fHistoFDTemplatePromptName << std::endl;
  std::cout << "Histo FD template Non Prompt     = " << fHistoFDTemplateNonPromptName << std::endl;
  std::cout << "Histo FD Prompt Frac     = " << fHistoFDPromptFracName << std::endl;
  std::cout << "===================== " << std::endl;
  std::cout << " " << std::endl;

  return kTRUE;
}

Bool_t DhCorrelationExtraction::ReadInputSecondaryPartContamination()
{

  fFileSecPart = TFile::Open(fFileSecPartName.Data());
  if (!fFileSecPart) {
    std::cout << "[ERROR] File " << fFileSecPartName << " cannot be opened! check your file path!" << std::endl;
    return kFALSE;
  }

  fDirSecPart = reinterpret_cast<TDirectoryFile*>(fFileSecPart->Get(fDirSecPartName.Data()));

  if (!fDirSecPart) {
    std::cout << "[ERROR] Directory " << fDirSecPart << " cannot be opened! check your file path!" << std::endl;
    return kFALSE;
  }

  std::cout << "===================== " << std::endl;
  std::cout << "Read inputs SE and ME" << std::endl;
  std::cout << "TFile Sec. part.    = " << fFileSecPartName << std::endl;
  std::cout << "TDir Sec. part.    = " << fDirSecPartName << std::endl;
  std::cout << "===================== " << std::endl;
  std::cout << " " << std::endl;

  return kTRUE;
}

TH1D* DhCorrelationExtraction::EvaluateMCClosModulations(Double_t PtCandMin, Double_t PtCandMax, Double_t PtHadMin, Double_t PtHadMax)
{

  TH1D* hModul = new TH1D();

  fFilePromptMc = TFile::Open(fFilePromptMcRecName.Data());
  fFileNonPromptMc = TFile::Open(fFileNonPromptMcRecName.Data());

  if (!fFilePromptMc) {
    std::cout << "[ERROR] File prompt MC rec cannot be opened! check your file path!" << std::endl;
  }
  if (!fFileNonPromptMc) {
    std::cout << "[ERROR] File non-prompt MC rec cannot be opened! check your file path!" << std::endl;
  }

  // TODO: generalise this part
  TH1D* hRecPrompt = reinterpret_cast<TH1D*>(fFilePromptMc->Get(Form("h1D_Rec_iPtD%d_iPtAssoc%d", fBinPtCand, fBinPtHad)));
  TH1D* hRecNonPrompt = reinterpret_cast<TH1D*>(fFileNonPromptMc->Get(Form("h1D_Rec_iPtD%d_iPtAssoc%d", fBinPtCand, fBinPtHad)));
  TH1D* hGenPrompt = reinterpret_cast<TH1D*>(fFilePromptMc->Get(Form("h1D_Gen_iPtD%d_iPtAssoc%d", fBinPtCand, fBinPtHad)));
  TH1D* hGenNonPrompt = reinterpret_cast<TH1D*>(fFileNonPromptMc->Get(Form("h1D_Gen_iPtD%d_iPtAssoc%d", fBinPtCand, fBinPtHad)));

  printf("[INFO] Bin cand %d - Bin had %d \n", fBinPtCand, fBinPtHad);

  // hRecPrompt = ReflectCorrHistogram(hRecPrompt);
  // hRecNonPrompt = ReflectCorrHistogram(hRecNonPrompt);
  // hGenPrompt = ReflectCorrHistogram(hGenPrompt);
  // hGenNonPrompt = ReflectCorrHistogram(hGenNonPrompt);

  hRecNonPrompt->Sumw2();
  hRecNonPrompt->Sumw2();
  hGenPrompt->Sumw2();
  hGenNonPrompt->Sumw2();

  TH1D* hRatioNonPrompt = reinterpret_cast<TH1D*>(hRecNonPrompt->Clone("hRatioNonPrompt"));
  hRatioNonPrompt->Sumw2();
  hRatioNonPrompt->Divide(hRecNonPrompt, hGenNonPrompt, 1., 1., "B");
  hModul = reinterpret_cast<TH1D*>(hRatioNonPrompt->Clone("hModul"));

  TF1* funFit = new TF1("funFit", "[0]", TMath::Pi() * 3. / 8., TMath::Pi() * 3 / 2);
  hRatioNonPrompt->Fit(funFit, "R");
  Double_t fitVal = funFit->GetParameter(0);

  TCanvas* cRatio_MCClosure = new TCanvas(Form("cRatio_MCClosure_PtCand%.0fto%.0f_Pthad%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax), Form("cRatio_MCClosure_PtCand%.0fto%.0f_Pthad%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax), 100, 100, 1200, 700);
  cRatio_MCClosure->cd();
  hRatioNonPrompt->GetYaxis()->SetRangeUser(0.2, 1.8);
  hRatioNonPrompt->Draw();

  Double_t FPrompt = GetFDPromptFrac(PtCandMin, PtCandMax, PtHadMin, PtHadMax);
  Double_t relAmplC[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t relAmplB[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t recoKineVal[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t modul[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  for (int iBin = 0; iBin < hRatioNonPrompt->GetNbinsX(); iBin++) {
    if (iBin > 1 && iBin < 13) {
      recoKineVal[iBin - 2] = hRatioNonPrompt->GetBinContent(iBin + 1) - (fitVal - 1);
      relAmplC[iBin - 2] = hRecPrompt->GetBinContent(iBin + 1) / (hRecPrompt->GetBinContent(iBin + 1) * FPrompt + hRecNonPrompt->GetBinContent(iBin + 1) * (1 - FPrompt));
      relAmplB[iBin - 2] = hRecNonPrompt->GetBinContent(iBin + 1) / (hRecPrompt->GetBinContent(iBin + 1) * FPrompt + hRecNonPrompt->GetBinContent(iBin + 1) * (1 - FPrompt));
      modul[iBin - 2] = relAmplC[iBin - 2] * FPrompt + relAmplB[iBin - 2] * (1 - FPrompt) / recoKineVal[iBin - 2];
      hModul->SetBinContent(iBin + 1, modul[iBin - 2]);
      hModul->SetBinError(iBin + 1, 0.);

      printf("[INFO] Bin%d MODUL = %1.5f\t (Reco/Kine-fitVal = %1.4f, FPrompt = %1.3f, Ampl_ratio C,B = %1.4f, %1.4f)\n", iBin + 1, modul[iBin - 2], recoKineVal[iBin - 2], FPrompt, relAmplC[iBin - 2], relAmplB[iBin - 2]);
    } else {
      hModul->SetBinContent(iBin + 1, 1.);
      hModul->SetBinError(iBin + 1, 0.);
    }
  }

  hModul->SetLineColor(kMagenta);
  hModul->Draw("same");

  cRatio_MCClosure->SaveAs(Form("Output_CorrelationExtraction_Thin2023_FullAnalysis_CentralPoints_png/Ratio_MCClosure_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", PtCandMin, PtCandMax, PtHadMin, PtHadMax));

  return hModul;
}

TH2D* DhCorrelationExtraction::GetCorrelHisto(Int_t SEorME, Int_t SorSB, Int_t pool, Double_t PtCandMin, Double_t PtCandMax, Double_t PtHadMin, Double_t PtHadMax)
{
  // TODO: Subtraction of softpion
  TH2D* h2D = new TH2D(); // pointer to be returned

  THnSparseD* hSparse = 0x0;
  if (SEorME == kSE) { // Same Event
    if (SorSB == kSign) {
      hSparse = reinterpret_cast<THnSparseD*>(fDirSE->Get(fSECorrelSignalRegionName.Data()));
    } else if (!fSidebandDivided) {
      hSparse = reinterpret_cast<THnSparseD*>(fDirSE->Get(fSECorrelSidebandsName.Data()));
    } else if (fSidebandDivided) {
      if (fUseSidebLeft && !fUseSidebRight) {
        hSparse = reinterpret_cast<THnSparseD*>(fDirSE->Get(fSECorrelSidebandLeftName.Data()));
      } else if (!fUseSidebLeft && fUseSidebRight) {
        hSparse = reinterpret_cast<THnSparseD*>(fDirSE->Get(fSECorrelSidebandRightName.Data()));
      } else if (fUseSidebLeft && fUseSidebRight) {
        hSparse = reinterpret_cast<THnSparseD*>(fDirSE->Get(fSECorrelSidebandLeftName.Data()));
        hSparse->SetName("hSparse");
        THnSparseD* hSparseRightSideb = reinterpret_cast<THnSparseD*>(fDirSE->Get(fSECorrelSidebandRightName.Data()));
        hSparse->Add(hSparseRightSideb, 1.);
      }
    }
  } else { // Mixed Event
    if (SorSB == kSign) {
      hSparse = reinterpret_cast<THnSparseD*>(fDirME->Get(fMECorrelSignalRegionName.Data()));
    } else if (!fSidebandDivided) {
      hSparse = reinterpret_cast<THnSparseD*>(fDirME->Get(fMECorrelSidebandsName.Data()));
    } else if (fSidebandDivided) {
      if (fUseSidebLeft && !fUseSidebRight) {
        hSparse = reinterpret_cast<THnSparseD*>(fDirME->Get(fMECorrelSidebandLeftName.Data()));
      } else if (!fUseSidebLeft && fUseSidebRight) {
        hSparse = reinterpret_cast<THnSparseD*>(fDirME->Get(fMECorrelSidebandRightName.Data()));
      } else if (fUseSidebLeft && fUseSidebRight) {
        hSparse = reinterpret_cast<THnSparseD*>(fDirME->Get(fMECorrelSidebandLeftName.Data()));
        hSparse->SetName("hSparse");
        THnSparseD* hSparseRightSideb = reinterpret_cast<THnSparseD*>(fDirME->Get(fMECorrelSidebandRightName.Data()));
        hSparse->Add(hSparseRightSideb, 1.);
      }
    }
  }
  /*else if (fSidebandDivided) { // Mixed Event
    if (SorSB == kSign) { hSparse = reinterpret_cast<THnSparseD*> fDirME -> Get(fMECorrelSignalRegionName.Data());
    } else if (!fSidebandDivided) { hSparse = reinterpret_cast<THnSparseD*> fDirME -> Get(fMECorrelSidebandsName.Data());
    } else {
      if (fUseSidebLeft && !fUseSidebRight) {
        hSparse = reinterpret_cast<THnSparseD*> fDirME -> Get(fMECorrelSidebandLeftName.Data());
      } else if (!fUseSidebLeft && fUseSidebRight) {
        hSparse = reinterpret_cast<THnSparseD*> fDirME -> Get(fMECorrelSidebandRightName.Data());
      } else if (fUseSidebLeft && fUseSidebRight) {
        hSparse = reinterpret_cast<THnSparseD*> fDirME -> Get(fMECorrelSidebandLeftName.Data());
        hSparse -> SetName("hSparse");
        THnSparseD *hSparseRightSideb = reinterpret_cast<THnSparseD*> fDirME -> Get(fMECorrelSidebandRightName.Data());
        hSparse -> Add(hSparseRightSideb, 1.);
      }
    }
  }*/

  Int_t binExtPtCandMin = (Int_t)hSparse->GetAxis(2)->FindBin(PtCandMin + 0.01); // axis2: ptCand, the 0.01 to avoid bin edges!
  Int_t binExtPtCandMax = (Int_t)hSparse->GetAxis(2)->FindBin(PtCandMax - 0.01);
  Int_t binExtPtHadMin = (Int_t)hSparse->GetAxis(3)->FindBin(PtHadMin + 0.01); // axis3: ptHad
  Int_t binExtPtHadMax = (Int_t)hSparse->GetAxis(3)->FindBin(PtHadMax - 0.01);
  Int_t binExtPoolMin;
  Int_t binExtPoolMax;
  if (fCorrectPoolsSeparately) {
    binExtPoolMin = (Int_t)hSparse->GetAxis(4)->FindBin(pool + 0.01); // axis4: pool bin
    binExtPoolMax = (Int_t)hSparse->GetAxis(4)->FindBin(pool + 0.99);
  } else { // merge all pools in one
    binExtPoolMin = 1;
    binExtPoolMax = (Int_t)hSparse->GetAxis(4)->GetNbins();
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
  hSparse->GetAxis(4)->SetRange(binExtPoolMin, binExtPoolMax);     // axis4: pool bin
  h2D = reinterpret_cast<TH2D*>(hSparse->Projection(0, 1));        // axis0: deltaPhi, axis1: deltaEta
  if (SEorME == kSE) {                                             // Same Event
    if (SorSB == kSign) {
      h2D->SetName(Form("hCorr_SE_Sig_2D_PtCandBin%d_PtHadBin%d_iPool%d", binExtPtCandMin, binExtPtHadMin, pool));
    } else if (!fSidebandDivided) {
      h2D->SetName(Form("hCorr_SE_Sideb_2D_PtCandBin%d_PtHadBin%d_iPool%d", binExtPtCandMin, binExtPtHadMin, pool));
    } else {
      if (fUseSidebLeft && !fUseSidebRight) {
        h2D->SetName(Form("hCorr_SE_Sideb_Left_2D_PtCandBin%d_PtHadBin%d_iPool%d", binExtPtCandMin, binExtPtHadMin, pool));
      } else if (!fUseSidebLeft && fUseSidebRight) {
        h2D->SetName(Form("hCorr_SE_Sideb_Right_2D_PtCandBin%d_PtHadBin%d_iPool%d", binExtPtCandMin, binExtPtHadMin, pool));
      } else if (fUseSidebLeft && fUseSidebRight) {
        h2D->SetName(Form("hCorr_SE_Sideb_LeftAndRight_2D_PtCandBin%d_PtHadBin%d_iPool%d", binExtPtCandMin, binExtPtHadMin, pool));
      }
    }
  } else { // Mixed Event
    if (SorSB == kSign) {
      h2D->SetName(Form("hCorr_ME_Sig_2D_PtCandBin%d_PtHadBin%d_iPool%d", binExtPtCandMin, binExtPtHadMin, pool));
    } else if (!fSidebandDivided) {
      h2D->SetName(Form("hCorr_SE_Sideb_2D_PtCandBin%d_PtHadBin%d_iPool%d", binExtPtCandMin, binExtPtHadMin, pool));
    } else {
      if (fUseSidebLeft && !fUseSidebRight) {
        h2D->SetName(Form("hCorr_ME_Sideb_Left_2D_PtCandBin%d_PtHadBin%d_iPool%d", binExtPtCandMin, binExtPtHadMin, pool));
      } else if (!fUseSidebLeft && fUseSidebRight) {
        h2D->SetName(Form("hCorr_ME_Sideb_Right_2D_PtCandBin%d_PtHadBin%d_iPool%d", binExtPtCandMin, binExtPtHadMin, pool));
      } else if (fUseSidebLeft && fUseSidebRight) {
        h2D->SetName(Form("hCorr_ME_Sideb_LeftAndRight_2D_PtCandBin%d_PtHadBin%d_iPool%d", binExtPtCandMin, binExtPtHadMin, pool));
      }
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
  TH1F* hMassFitSBLYield = reinterpret_cast<TH1F*>(fFileMass->Get("hBackgroundSidebandLeft"));
  TH1F* hMassFitSBRYield = reinterpret_cast<TH1F*>(fFileMass->Get("hBackgroundSidebandRight"));

  Int_t PtCandBin = hMassFitSgnYield->FindBin(PtCandMin + 0.01);
  if (PtCandBin != hMassFitSgnYield->FindBin(PtCandMax - 0.01))
    std::cout << "[ERROR] Pt bin in invariant mass histogram not univocally defined " << std::endl;

  Float_t SgnYield = hMassFitSgnYield->GetBinContent(PtCandBin);
  Float_t BkgYield = hMassFitBkgYield->GetBinContent(PtCandBin);
  Float_t SBsYield = hMassFitSBsYield->GetBinContent(PtCandBin);
  Float_t SBLYield = hMassFitSBLYield->GetBinContent(PtCandBin);
  Float_t SBRYield = hMassFitSBRYield->GetBinContent(PtCandBin);

  std::cout << "================================= " << std::endl;
  std::cout << "Getting invariant mass parameters " << std::endl;
  std::cout << "Pt cand " << PtCandMin << " - " << PtCandMax << std::endl;
  std::cout << "Signal yield    = " << SgnYield << std::endl;
  std::cout << "Bkg yield    = " << BkgYield << std::endl;
  std::cout << "Sideband yield    = " << SBsYield << std::endl;
  std::cout << "Sideband left yield    = " << SBLYield << std::endl;
  std::cout << "Sideband right yield    = " << SBRYield << std::endl;
  std::cout << "================================= " << std::endl;
  std::cout << " " << std::endl;

  SetSignalYieldforNorm(SgnYield);
  SetBkgYield(BkgYield);
  if (fUseSidebLeft && fUseSidebRight) {
    SetBkgScaleFactor(BkgYield / SBsYield);
    SetSBYield(SBsYield);
  } else if (fUseSidebLeft && !fUseSidebRight) {
    SetBkgScaleFactor(BkgYield / SBLYield);
    SetSBYield(SBLYield);
  } else if (!fUseSidebLeft && fUseSidebRight) {
    SetBkgScaleFactor(BkgYield / SBRYield);
    SetSBYield(SBRYield);
  }

  return;
}

TH2D* DhCorrelationExtraction::GetFDTemplateHisto(Int_t PromptOrFD, Double_t PtCandMin, Double_t PtCandMax, Double_t PtHadMin, Double_t PtHadMax)
{

  TH2D* h2D = new TH2D(); // pointer to be returned

  if (PromptOrFD == kPrompt) {
    h2D = reinterpret_cast<TH2D*>(fFileFDTemplate->Get(Form("%s%.0f_%.0f_ptassoc%.0f_%.0f", fHistoFDTemplatePromptName.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax)));
  } else {
    h2D = reinterpret_cast<TH2D*>(fFileFDTemplate->Get(Form("%s%.0f_%.0f_ptassoc%.0f_%.0f", fHistoFDTemplateNonPromptName.Data(), PtCandMin, PtCandMax, PtHadMin, PtHadMax)));
  }

  Int_t binExtEtaMin = (Int_t)h2D->GetXaxis()->FindBin(fDeltaEtaMin + 0.000001);
  Int_t binExtEtaMax = (Int_t)h2D->GetXaxis()->FindBin(fDeltaEtaMax - 0.000001);
  if (binExtEtaMax > h2D->GetXaxis()->GetNbins())
    binExtEtaMax = h2D->GetXaxis()->GetNbins();
  if (binExtEtaMin < 1)
    binExtEtaMin = 1;

  h2D->GetXaxis()->SetRange(binExtEtaMin, binExtEtaMax);
  if (PromptOrFD == kPrompt) {
    h2D->SetName(Form("hFDTemplatePrompt_2D_PtCand%.0fto%.0f_PtHad%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax));
  } else {
    h2D->SetName(Form("hFDTemplateNonPrompt_2D_PtCand%.0fto%.0f_PtHad%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax));
  }
  h2D->GetYaxis()->SetTitle("#Delta#phi (rad)");
  h2D->GetXaxis()->SetTitle("#Delta#eta");

  return h2D;
}

TH1D* DhCorrelationExtraction::GetCorrelHistoSecondaryPart(Int_t PartType, Double_t PtCandMin, Double_t PtCandMax, Double_t PtHadMin, Double_t PtHadMax)
{

  TH1D* h1D = new TH1D(); // pointer to be returned

  THnSparseD* hSparse = 0x0;

  if (PartType == kPrimaryPart) { // primary particles
    hSparse = reinterpret_cast<THnSparseD*>(fDirSecPart->Get(fHistoPrimaryPartName.Data()));
  } else { // all selected particles
    hSparse = reinterpret_cast<THnSparseD*>(fDirSecPart->Get(fHistoAllPartName.Data()));
  }
  Int_t binExtPtCandMin = (Int_t)hSparse->GetAxis(2)->FindBin(PtCandMin + 0.01); // axis2: ptCand, the 0.01 to avoid bin edges!
  Int_t binExtPtCandMax = (Int_t)hSparse->GetAxis(2)->FindBin(PtCandMax - 0.01);
  Int_t binExtPtHadMin = (Int_t)hSparse->GetAxis(3)->FindBin(PtHadMin + 0.01); // axis3: ptHad
  Int_t binExtPtHadMax = (Int_t)hSparse->GetAxis(3)->FindBin(PtHadMax - 0.01);
  Int_t binExtPoolMin;
  Int_t binExtPoolMax;
  if (PartType == kAllPart) {
    binExtPoolMin = 1;
    binExtPoolMax = (Int_t)hSparse->GetAxis(4)->GetNbins();
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
  if (PartType == kAllPart) {
    hSparse->GetAxis(4)->SetRange(binExtPoolMin, binExtPoolMax); // axis4: pool bin
  }

  h1D = reinterpret_cast<TH1D*>(hSparse->Projection(0)); // axis0: deltaPhi
  if (PartType == kPrimaryPart) {                        // primary particles
    h1D->SetName(Form("hPrimaryPartCorr_PtD%.0fto%.0f_PtHad%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax));
  } else { // all selected particles
    h1D->SetName(Form("hAllPartCorr_PtD%.0fto%.0f_PtHad%.0fto%.0f", PtCandMin, PtCandMax, PtHadMin, PtHadMax));
  }

  return h1D;
}

TH1D* DhCorrelationExtraction::ReflectCorrHistogram(TH1D*& histo)
{

  // nBinsPhi must be a multple of 4 in order to reflect correcty the histogram
  Int_t nBinsPhi = histo->GetNbinsX();
  Int_t nBinsPhiRefl = nBinsPhi / 2;
  Int_t bin0Phi = nBinsPhi / 4 + 1;
  Int_t binPiPhi = 3 * nBinsPhi / 4;

  TH1D* h1D = new TH1D("h1D_Reflected", "", nBinsPhiRefl, 0., TMath::Pi()); // pointer to be returned
  h1D->Sumw2();
  // TH1D* h1D = reinterpret_cast<TH1D*> histo -> Clone("h1D_Reflected");
  // h1D -> GetXaxis() -> SetRange(bin0Phi, binPiPhi);

  // reflection
  Double_t reflectedContent, reflectedContentError;
  for (int iBin = 0; iBin < nBinsPhiRefl / 2; iBin++) {
    reflectedContent = (histo->GetBinContent(bin0Phi - iBin - 1) + histo->GetBinContent(bin0Phi + iBin)) / 2;
    std::cout << "[INFO] Paired bins for reflection: " << bin0Phi - iBin - 1 << " - " << bin0Phi + iBin << std::endl;
    std::cout << "[INFO] Bin filled: " << iBin + 1 << std::endl;
    reflectedContentError = 0.5 * TMath::Sqrt(TMath::Power(histo->GetBinError(iBin + 1), 2) + TMath::Power(histo->GetBinError(bin0Phi + iBin), 2));
    h1D->SetBinContent(iBin + 1, reflectedContent);
    h1D->SetBinError(iBin + 1, reflectedContentError);
  }
  for (int iBin = nBinsPhiRefl / 2; iBin < nBinsPhiRefl; iBin++) {
    reflectedContent = (histo->GetBinContent(bin0Phi + iBin) + histo->GetBinContent(binPiPhi + 2 * bin0Phi - iBin - 2)) / 2;
    reflectedContentError = 0.5 * TMath::Sqrt(TMath::Power(histo->GetBinError(bin0Phi + iBin), 2) + TMath::Power(histo->GetBinError(binPiPhi + 2 * bin0Phi - iBin - 2), 2));
    std::cout << "[INFO] Paired bins for reflection: " << bin0Phi + iBin << " - " << binPiPhi + 2 * bin0Phi - iBin - 2 << std::endl;
    std::cout << "[INFO] Bin filled: " << iBin + 1 << std::endl;
    h1D->SetBinContent(iBin + 1, reflectedContent);
    h1D->SetBinError(iBin + 1, reflectedContentError);
  }

  return h1D;
}

TH1D* DhCorrelationExtraction::ReflectHistoRun2(TH1D* h, Double_t scale)
{

  TH1D* h2 = new TH1D(Form("%sReflected", h->GetName()), Form("%sReflected", h->GetName()), h->GetNbinsX() / 2., 0., TMath::Pi());
  for (Int_t j = 1; j <= h->GetNbinsX(); j++) {
    Double_t x = h->GetBinCenter(j);
    Double_t y0 = h->GetBinContent(j);
    Double_t ey0 = h->GetBinError(j);
    Int_t j2;
    if (x > 0 && x < TMath::Pi()) {
      j2 = h2->FindBin(x);
    } else if (x < 0) {
      j2 = h2->FindBin(-1. * x);
    } else if (x > TMath::Pi()) {
      j2 = h2->FindBin(2. * TMath::Pi() - x);
    } else {
      printf("Point %d excluded \n", j);
      continue;
    }
    Double_t y = h2->GetBinContent(j2);
    Double_t ey = h2->GetBinError(j2);
    h2->SetBinContent(j2, (y + y0));
    h2->SetBinError(j2, TMath::Sqrt(ey0 * ey0 + ey * ey));
  }
  h2->Scale(scale);

  return h2;
}

Double_t DhCorrelationExtraction::GetFDPromptFrac(Double_t PtCandMin, Double_t PtCandMax, Double_t PtHadMin, Double_t PtHadMax)
{

  TH1D* h1D = new TH1D();
  h1D = reinterpret_cast<TH1D*>(fFileFDPromptFrac->Get(fHistoFDPromptFracName.Data()));

  Int_t binPtCandMin = h1D->GetXaxis()->FindBin(PtCandMin + 0.01);
  Int_t binPtCandMax = h1D->GetXaxis()->FindBin(PtCandMax - 0.01);
  Double_t PromptFraction;
  if (binPtCandMin == binPtCandMax) {
    PromptFraction = h1D->GetBinContent(binPtCandMin);
  } else {
    std::cout << "[ERROR] Different bin obtained from PtCandMin and PtCandMax";
    return 0.;
  }

  return PromptFraction;
}

void DhCorrelationExtraction::NormalizeMEplot(TH2D*& histoME, TH2D*& histoMEsoftPi)
{

  Int_t bin0phi = histoME->GetYaxis()->FindBin(0.);
  Int_t bin0eta = histoME->GetXaxis()->FindBin(0.);

  // evaluate the normalization (from ALL tracks, including possible fake softpions) -> **histoME indeed includes bin1+bin2 of THnSparse, i.e. all the tracks**
  Double_t factorNorm = 0;
  for (int in = -1; in <= 0; in++) {
    factorNorm += histoME->GetBinContent(bin0eta, bin0phi + in);
  }
  for (int in = -1; in <= 0; in++) {
    factorNorm += histoME->GetBinContent(bin0eta - 1, bin0phi + in);
  }
  factorNorm /= 4.;

  std::cout << "bin 0 phi: " << bin0phi << std::endl;
  std::cout << "bin 0 eta: " << bin0eta << std::endl;
  std::cout << "Factor norm. ME: " << factorNorm << std::endl;
  std::cout << "Bin content (0,0) ME: " << histoME->GetBinContent(bin0eta, bin0phi) << std::endl;

  if (fSubtractSoftPiME)
    histoME->Add(histoMEsoftPi, -1); // remove the tracks compatible with soft pion (if requested)

  // apply the normalization
  histoME->Scale(1. / factorNorm);

  return;
}

Double_t DhCorrelationExtraction::CalculateBaseline(TH1D*& histo, Bool_t totalRange, Bool_t reflected)
{

  // total range = 2*Pi
  // half range = Pi , for histogram reflected under symmetric assumption

  Double_t baseline, errBaseline;
  Int_t nBinsPhi = histo->GetNbinsX();
  Int_t binPhiHalf = nBinsPhi / 2;
  Int_t binPhiHalfMinus1 = nBinsPhi / 2 - 1;
  Int_t binPhiHalfPlus1 = nBinsPhi / 2 + 1;
  Int_t binPhiHalfPlus2 = nBinsPhi / 2 + 1;

  if (totalRange) {
    // baseline evaluated considering: the two first points, the last two points and four points in the middle (corresponding to the outer points)
    if (nBinsPhi >= 32) {
      baseline =
        ((histo->GetBinContent(1)) * (1. / TMath::Power(histo->GetBinError(1), 2)) +
         (histo->GetBinContent(2)) * (1. / TMath::Power(histo->GetBinError(2), 2)) +
         (histo->GetBinContent(binPhiHalfMinus1)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfMinus1), 2)) +
         (histo->GetBinContent(binPhiHalf)) * (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
         (histo->GetBinContent(binPhiHalfPlus1)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
         (histo->GetBinContent(binPhiHalfPlus2)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus2), 2)) +
         (histo->GetBinContent(nBinsPhi - 1)) * (1. / TMath::Power(histo->GetBinError(nBinsPhi - 1), 2)) +
         (histo->GetBinContent(nBinsPhi)) * (1. / TMath::Power(histo->GetBinError(nBinsPhi), 2))) /
        ((1. / TMath::Power(histo->GetBinError(1), 2)) +
         (1. / TMath::Power(histo->GetBinError(2), 2)) +
         (1. / TMath::Power(histo->GetBinError(binPhiHalfMinus1), 2)) +
         (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
         (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
         (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus2), 2)) +
         (1. / TMath::Power(histo->GetBinError(nBinsPhi - 1), 2)) +
         (1. / TMath::Power(histo->GetBinError(nBinsPhi), 2)));
    } else {
      baseline =
        ((histo->GetBinContent(1)) * (1. / TMath::Power(histo->GetBinError(1), 2)) +
         (histo->GetBinContent(binPhiHalf)) * (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
         (histo->GetBinContent(binPhiHalfPlus1)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
         (histo->GetBinContent(nBinsPhi)) * (1. / TMath::Power(histo->GetBinError(nBinsPhi), 2))) /
        ((1. / TMath::Power(histo->GetBinError(1), 2)) +
         (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
         (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
         (1. / TMath::Power(histo->GetBinError(nBinsPhi), 2)));
    }
  } else {
    if (reflected) {
      baseline =
        ((histo->GetBinContent(binPhiHalfMinus1)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfMinus1), 2)) +
         (histo->GetBinContent(binPhiHalf)) * (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
         (histo->GetBinContent(binPhiHalfPlus1)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
         (histo->GetBinContent(binPhiHalfPlus2)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus2), 2))) /
        ((1. / TMath::Power(histo->GetBinError(binPhiHalfMinus1), 2)) +
         (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
         (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
         (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus2), 2)));
    } else {
      // baseline evaluated using the 4 middle points in the transverese region
      if (nBinsPhi >= 32) {
        baseline =
          ((histo->GetBinContent(binPhiHalfMinus1)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfMinus1), 2)) +
           (histo->GetBinContent(binPhiHalf)) * (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
           (histo->GetBinContent(binPhiHalfPlus1)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
           (histo->GetBinContent(binPhiHalfPlus2)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus2), 2))) /
          ((1. / TMath::Power(histo->GetBinError(binPhiHalfMinus1), 2)) +
           (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
           (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
           (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus2), 2)));
      } else {
        baseline =
          ((histo->GetBinContent(binPhiHalf)) * (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
           (histo->GetBinContent(binPhiHalfPlus1)) * (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2))) /
          ((1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
           (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)));
      }
    }
  }

  return baseline;
}

Double_t DhCorrelationExtraction::CalculateBaselineError(TH1D*& histo, Bool_t totalRange, Bool_t reflected)
{

  // total range = 2*Pi
  // half range = Pi , for histogram reflected under symmetric assumption

  Double_t errBaseline;
  Int_t nBinsPhi = histo->GetNbinsX();
  Int_t binPhiHalf = nBinsPhi / 2;
  Int_t binPhiHalfMinus1 = nBinsPhi / 2 - 1;
  Int_t binPhiHalfPlus1 = nBinsPhi / 2 + 1;
  Int_t binPhiHalfPlus2 = nBinsPhi / 2 + 1;

  if (totalRange) {
    // baseline evaluated considering: the two first points, the last two points and four points in the middle (corresponding to the outer points)
    if (nBinsPhi >= 32) {
      errBaseline = 1. /
                    TMath::Sqrt((1. / TMath::Power(histo->GetBinError(1), 2)) +
                                (1. / TMath::Power(histo->GetBinError(2), 2)) +
                                (1. / TMath::Power(histo->GetBinError(binPhiHalfMinus1), 2)) +
                                (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
                                (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
                                (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus2), 2)) +
                                (1. / TMath::Power(histo->GetBinError(nBinsPhi - 1), 2)) +
                                (1. / TMath::Power(histo->GetBinError(nBinsPhi), 2)));
    } else { // fon nBinsPhi = 16 (rebin 4)
      errBaseline = 1. /
                    TMath::Sqrt((1. / TMath::Power(histo->GetBinError(1), 2)) +
                                (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
                                (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
                                (1. / TMath::Power(histo->GetBinError(nBinsPhi), 2)));
    }
  } else {
    // baseline evaluated using the 4 middle points in the transverese region
    if (reflected) {
      errBaseline = 1. /
                    TMath::Sqrt((1. / TMath::Power(histo->GetBinError(binPhiHalfMinus1), 2)) +
                                (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
                                (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
                                (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus2), 2)));
    } else {
      if (nBinsPhi >= 32) {
        errBaseline = 1. /
                      TMath::Sqrt((1. / TMath::Power(histo->GetBinError(binPhiHalfMinus1), 2)) +
                                  (1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
                                  (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)) +
                                  (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus2), 2)));
      } else {
        errBaseline = 1. /
                      TMath::Sqrt((1. / TMath::Power(histo->GetBinError(binPhiHalf), 2)) +
                                  (1. / TMath::Power(histo->GetBinError(binPhiHalfPlus1), 2)));
      }
    }
  }

  return errBaseline;
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
