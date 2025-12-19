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

#include <TAttMarker.h>
#include <TCanvas.h>
#include <TDirectoryFile.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TLegend.h>
#include <TMath.h>
#include <TString.h>

#include <Rtypes.h>
#include <RtypesCore.h>

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

DhCorrelationExtraction::DhCorrelationExtraction() : // default constructor
                                                     fFileMass(nullptr),
                                                     fFileSE(nullptr),
                                                     fFileME(nullptr),
                                                     fFileFDTemplate(nullptr),
                                                     fFileFDPromptFrac(nullptr),
                                                     fFileSecPart(nullptr),
                                                     fFilePromptMc(nullptr),
                                                     fFileNonPromptMc(nullptr),
                                                     fDirMass(nullptr),
                                                     fDirSE(nullptr),
                                                     fDirME(nullptr),
                                                     fDirSecPart(nullptr),
                                                     fCorrectedCorrHisto(nullptr),
                                                     fCorrectedCorrHistoBaselineSubtr(nullptr),
                                                     fCorrectedCorrHistoReflected(nullptr),
                                                     fCorrectedCorrHistoReflectedBaselineSubtr(nullptr),
                                                     fDmesonSpecies(kDsToKKPi),
                                                     fDmesonLabel("Ds"),
                                                     fFileNameSE(""),
                                                     fFileNameME(""),
                                                     fFileSecPartName(""),
                                                     fFileFDTemplateName(""),
                                                     fFileFDPromptFracName(""),
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
                                                     fHistoFDTemplatePromptName(""),
                                                     fHistoFDTemplateNonPromptName(""),
                                                     fHistoFDPromptFracName(""),
                                                     fHistoPrimaryPartName(""),
                                                     fHistoAllPartName(""),
                                                     fNpools(9),
                                                     fRebinAxisDeltaEta(1),
                                                     fRebinAxisDeltaPhi(1),
                                                     fDebug(0),
                                                     fBinPtCand(0),
                                                     fBinPtHad(0),
                                                     fDeltaEtaMin(-1.),
                                                     fDeltaEtaMax(1.),
                                                     fBkgScaleFactor(1.),
                                                     fSgnYieldNorm(1.),
                                                     fBkgYield(1.),
                                                     fCorrectPoolsSeparately(kTRUE),
                                                     fSubtractSoftPiME(kFALSE),
                                                     fRebinAngCorr(kFALSE),
                                                     fRebinFDCorr(kFALSE),
                                                     fRebinSecPart(kFALSE),
                                                     fSidebandDivided(kFALSE),
                                                     fUseSidebLeft(kFALSE),
                                                     fUseSidebRight(kFALSE),
                                                     fFDsubtraction(false),
                                                     fSecPartContamination(false),
                                                     fCorrBiasBtoD(false)
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
                                                                                          fDirSE(source.fDirSE),
                                                                                          fDirME(source.fDirME),
                                                                                          fDirSecPart(source.fDirSecPart),
                                                                                          fCorrectedCorrHisto(source.fCorrectedCorrHisto),
                                                                                          fCorrectedCorrHistoBaselineSubtr(source.fCorrectedCorrHistoBaselineSubtr),
                                                                                          fCorrectedCorrHistoReflected(source.fCorrectedCorrHistoReflected),
                                                                                          fCorrectedCorrHistoReflectedBaselineSubtr(source.fCorrectedCorrHistoReflectedBaselineSubtr),
                                                                                          fDmesonSpecies(source.fDmesonSpecies),
                                                                                          fDmesonLabel(source.fDmesonLabel),
                                                                                          fFileNameSE(source.fFileNameSE),
                                                                                          fFileNameME(source.fFileNameME),
                                                                                          fFileSecPartName(source.fFileSecPartName),
                                                                                          fFileFDTemplateName(source.fFileFDTemplateName),
                                                                                          fFileFDPromptFracName(source.fFileFDPromptFracName),
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
                                                                                          fHistoFDTemplatePromptName(source.fHistoFDTemplatePromptName),
                                                                                          fHistoFDTemplateNonPromptName(source.fHistoFDTemplateNonPromptName),
                                                                                          fHistoFDPromptFracName(source.fHistoFDPromptFracName),
                                                                                          fHistoPrimaryPartName(source.fHistoPrimaryPartName),
                                                                                          fHistoAllPartName(source.fHistoAllPartName),
                                                                                          fNpools(source.fNpools),
                                                                                          fRebinAxisDeltaEta(source.fRebinAxisDeltaEta),
                                                                                          fRebinAxisDeltaPhi(source.fRebinAxisDeltaPhi),
                                                                                          fDebug(source.fDebug),
                                                                                          fBinPtCand(source.fBinPtCand),
                                                                                          fBinPtHad(source.fBinPtHad),
                                                                                          fDeltaEtaMin(source.fDeltaEtaMin),
                                                                                          fDeltaEtaMax(source.fDeltaEtaMax),
                                                                                          fBkgScaleFactor(source.fBkgScaleFactor),
                                                                                          fSgnYieldNorm(source.fSgnYieldNorm),
                                                                                          fBkgYield(source.fBkgYield),
                                                                                          fCorrectPoolsSeparately(source.fCorrectPoolsSeparately),
                                                                                          fSubtractSoftPiME(source.fSubtractSoftPiME),
                                                                                          fRebinAngCorr(source.fRebinAngCorr),
                                                                                          fRebinFDCorr(source.fRebinFDCorr),
                                                                                          fRebinSecPart(source.fRebinSecPart),
                                                                                          fSidebandDivided(source.fSidebandDivided),
                                                                                          fUseSidebLeft(source.fUseSidebLeft),
                                                                                          fUseSidebRight(source.fUseSidebRight),
                                                                                          fFDsubtraction(source.fFDsubtraction),
                                                                                          fSecPartContamination(source.fSecPartContamination),
                                                                                          fCorrBiasBtoD(source.fCorrBiasBtoD)
{
}

DhCorrelationExtraction::~DhCorrelationExtraction()
  // destructor
  = default;

Bool_t DhCorrelationExtraction::setDmesonSpecie(DmesonSpecie k)
{

  if (k < 0 || k > 3) {
    printf("[ERROR] D meson specie not correctly set!\n");
    return kFALSE;
  }
  if (k == 0) {
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

Bool_t DhCorrelationExtraction::extractCorrelations(Double_t ptCandMin, Double_t ptCandMax, Double_t ptHadMin, Double_t ptHadMax, TString codeName)
{

  if (fSubtractSoftPiME) {
    printf("[INFO] Fake softPi subtraction in ME via extraction code is enabled!\n");
  }

  if (!fCorrectPoolsSeparately) {
    fNpools = 1; // single histogram with integrated pools
  }

  // Histograms definition
  TH2D* hSeSign[fNpools];
  TH2D* hMeSign[fNpools];
  TH2D* hMeSignSoftPi[fNpools];
  TH2D* hSeSideb[fNpools];
  TH2D* hMeSideb[fNpools];
  TH2D* hMeSidebSoftPi[fNpools];

  TH2D* hCorrSign[fNpools];
  TH2D* hCorrSideb[fNpools];

  TH2D* h2DSign;
  TH2D* h2DSideb;
  TH2D* h2DSubtr;

  TH2D* h2DFdTemplatePrompt;
  TH2D* h2DFdTemplateNonPrompt;

  TH1D* h1DSign;
  TH1D* h1DSideb;
  TH1D* h1DSubtr;
  TH1D* h1DSignNorm;
  TH1D* h1DSidebNorm;
  TH1D* h1DSubtrNorm;
  TH1D* h1DFdTemplatePrompt;
  TH1D* h1DFdTemplateNonPrompt;
  TH1D* h1DTemplateTotal;
  TH1D* h1DSubtrFdNorm;
  TH1D* h1DPrimaryPartCorr;
  TH1D* h1DAllPartCorr;
  TH1D* h1DSecPartFrac;
  TH1D* h1DSubtrNormSecPart;
  TH1D* h1DBaselineSubtr;
  TH1D* h1DReflCorr;
  TH1D* h1DReflCorrBaselineSubtr;
  TH1D* hModul;
  TH1D* hBeforeModulCorr;

  Double_t fdPromptFrac;

  // if (fIntegratePtBins && iBinPtHad>0) continue;

  for (int iPool = 0; iPool < fNpools; iPool++) {
    // Retrieve 2D plots for SE and ME, signal and bkg regions, for each pTbin and pool
    hSeSign[iPool] = getCorrelHisto(kSE, kSign, iPool, ptCandMin, ptCandMax, ptHadMin, ptHadMax);
    std::cout << "Got SE histogram signal region" << std::endl;
    hMeSign[iPool] = getCorrelHisto(kME, kSign, iPool, ptCandMin, ptCandMax, ptHadMin, ptHadMax);
    std::cout << "Got ME histogram signal region" << std::endl;
    hSeSideb[iPool] = getCorrelHisto(kSE, kSideb, iPool, ptCandMin, ptCandMax, ptHadMin, ptHadMax);
    std::cout << "Got SE histogram sdeband region" << std::endl;
    hMeSideb[iPool] = getCorrelHisto(kME, kSideb, iPool, ptCandMin, ptCandMax, ptHadMin, ptHadMax);
    std::cout << "Got ME histogram sdeband region" << std::endl;

    hSeSign[iPool]->Sumw2();
    hMeSign[iPool]->Sumw2();
    hSeSideb[iPool]->Sumw2();
    hMeSideb[iPool]->Sumw2();

    // rebin axes deltaEta and deltaPhi
    if (fRebinAngCorr) {
      hSeSign[iPool]->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi); // Xaxis: deltaEta, Yaxis: deltaPhi
      hSeSideb[iPool]->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi);
      hMeSign[iPool]->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi);
      hMeSideb[iPool]->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi);
      if (fSubtractSoftPiME) {
        hMeSidebSoftPi[iPool]->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi);
      }
      std::cout << "SE and ME histograms rebinned" << std::endl;
    }

    if (fDebug >= 1) {
      auto* c = new TCanvas(Form("cSE_Original_%d_%1.1fto%1.1f", iPool, ptHadMin, ptHadMax), Form("cSE_Original_%s_pool%d_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), iPool, ptCandMin, ptCandMax, ptHadMin, ptHadMax), 100, 100, 1600, 900);
      c->Divide(2, 1);
      c->cd(1);
      hSeSign[iPool]->SetMinimum(0);
      hSeSign[iPool]->Draw("lego2");
      c->cd(2);
      hSeSideb[iPool]->SetMinimum(0);
      hSeSideb[iPool]->Draw("lego2");
      c->SaveAs(Form("Output_CorrelationExtraction_%s_png/CorrSE_Original_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
      c->SaveAs(Form("Output_CorrelationExtraction_%s_Root/CorrSE_Original_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
    }

    if (fDebug >= 1) {
      auto* c = new TCanvas(Form("cME_Original_%d_%1.1fto%1.1f", iPool, ptHadMin, ptHadMax), Form("cME_Original_%s_pool%d_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), iPool, ptCandMin, ptCandMax, ptHadMin, ptHadMax), 100, 100, 1600, 900);
      c->Divide(2, 1);
      c->cd(1);
      hMeSign[iPool]->SetMinimum(0);
      hMeSign[iPool]->Draw("lego2");
      c->cd(2);
      hMeSideb[iPool]->SetMinimum(0);
      hMeSideb[iPool]->Draw("lego2");
      c->SaveAs(Form("Output_CorrelationExtraction_%s_png/CorrME_Original_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
      c->SaveAs(Form("Output_CorrelationExtraction_%s_Root/CorrME_Original_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
    }

    // Scale bkg plots by ratio of signal region/sidebands
    hSeSideb[iPool]->Scale(fBkgScaleFactor);
    hMeSideb[iPool]->Scale(fBkgScaleFactor); // when normalised this factor should cancel out
    std::cout << "[INFO] fBkgScaleFactor    = " << fBkgScaleFactor << std::endl;
    hSeSideb[iPool]->SetEntries(hSeSideb[iPool]->GetEntries() * fBkgScaleFactor);
    hMeSideb[iPool]->SetEntries(hMeSideb[iPool]->GetEntries() * fBkgScaleFactor);

    if (fSubtractSoftPiME) {
      hMeSidebSoftPi[iPool]->Scale(fBkgScaleFactor);
      hMeSidebSoftPi[iPool]->SetEntries(hMeSidebSoftPi[iPool]->GetEntries() * fBkgScaleFactor);
    }

    // Normalize ME plots for the entries in (deltaEta, deltaPhi) = (0, 0)
    normalizeMePlot(hMeSign[iPool], hMeSignSoftPi[iPool]);
    normalizeMePlot(hMeSideb[iPool], hMeSidebSoftPi[iPool]);

    // Apply Event Mixing Correction
    hCorrSign[iPool] = reinterpret_cast<TH2D*>(hSeSign[iPool]->Clone(Form("hCorr_Sign_Pool%d", iPool)));
    hCorrSign[iPool]->Sumw2();
    hCorrSign[iPool]->Divide(hMeSign[iPool]);

    hCorrSideb[iPool] = reinterpret_cast<TH2D*>(hSeSideb[iPool]->Clone(Form("hCorr_Sideb_Pool%d", iPool)));
    hCorrSideb[iPool]->Sumw2();
    hCorrSideb[iPool]->Divide(hMeSideb[iPool]);

    Double_t nSEsign = 0, nSEsideb = 0, nSign = 0, nSideb = 0;
    for (int i = 1; i <= hCorrSign[iPool]->GetXaxis()->GetNbins(); i++) {
      for (int j = 1; j <= hCorrSign[iPool]->GetYaxis()->GetNbins(); j++) {
        nSEsign += hSeSign[iPool]->GetBinContent(i, j);
        nSEsideb += hSeSideb[iPool]->GetBinContent(i, j);
        nSign += hCorrSign[iPool]->GetBinContent(i, j);
        nSideb += hCorrSideb[iPool]->GetBinContent(i, j);
      }
    }
    hSeSign[iPool]->SetEntries(nSEsign);
    hSeSideb[iPool]->SetEntries(nSEsideb);
    hCorrSign[iPool]->SetEntries(nSign);
    hCorrSideb[iPool]->SetEntries(nSideb);

    if (fDebug >= 1) {
      auto* c = new TCanvas(Form("cSEME_%d_%1.1fto%1.1f", iPool, ptHadMin, ptHadMax), Form("cSEME_%s_pool%d_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), iPool, ptCandMin, ptCandMax, ptHadMin, ptHadMax), 100, 100, 1600, 900);
      c->Divide(3, 2);
      c->cd(1);
      hSeSign[iPool]->SetMinimum(0);
      hSeSign[iPool]->Draw("lego2");
      c->cd(2);
      hMeSign[iPool]->SetMinimum(0);
      hMeSign[iPool]->Draw("lego2");
      c->cd(3);
      hCorrSign[iPool]->SetMinimum(0);
      hCorrSign[iPool]->Draw("lego2");
      c->cd(4);
      hSeSideb[iPool]->SetMinimum(0);
      hSeSideb[iPool]->Draw("lego2");
      c->cd(5);
      hMeSideb[iPool]->SetMinimum(0);
      hMeSideb[iPool]->Draw("lego2");
      c->cd(6);
      hCorrSideb[iPool]->SetMinimum(0);
      hCorrSideb[iPool]->Draw("lego2");
      c->SaveAs(Form("Output_CorrelationExtraction_%s_png/CorrSEandME_%s_Canvas_PtCand%.0fto%.0f_Pool%d_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, iPool, ptHadMin, ptHadMax));
      c->SaveAs(Form("Output_CorrelationExtraction_%s_Root/CorrSEandME_%s_Canvas_PtCand%.0fto%.0f_Pool%d_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, iPool, ptHadMin, ptHadMax));
    }

    // Pools integration
    if (iPool == 0) {
      h2DSign = reinterpret_cast<TH2D*>(hCorrSign[0]->Clone("h2D_Sign"));
      h2DSideb = reinterpret_cast<TH2D*>(hCorrSideb[0]->Clone("h2D_Sideb"));
      h2DSign->Sumw2();
      h2DSideb->Sumw2();
    } else {
      h2DSign->Add(hCorrSign[iPool]);
      h2DSideb->Add(hCorrSideb[iPool]);
    }
  } // end pool loop

  // Draw 2D plots (Signal region and Sidebands)
  auto* c2D = new TCanvas(Form("c2D_IntPools_PtHad%.0fto%.0f", ptHadMin, ptHadMax), Form("c2D_%s_IntPools_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), ptHadMin, ptHadMax), 100, 100, 1500, 800);
  setTH2HistoStyle(h2DSign, Form("Signal region, %.0f < p^{%s}_{T} < %.0f GeV/c, %.0f < p^{assoc}_{T} < %.0f GeV/c", ptCandMin, fDmesonLabel.Data(), ptCandMax, ptHadMin, ptHadMax), "#Delta#eta", "#Delta#phi [rad]", "entries", 1.6, 1.6, 1.6, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04);
  setTH2HistoStyle(h2DSideb, Form("Sideband region, %.0f < p^{%s}_{T} < %.0f GeV/c, %.0f < p^{assoc}_{T} < %.0f GeV/c", ptCandMin, fDmesonLabel.Data(), ptCandMax, ptHadMin, ptHadMax), "#Delta#eta", "#Delta#phi [rad]", "#frac{Y_{Bkg}}{Y_{SB}} entries", 1.6, 1.6, 1.6, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04);
  c2D->Divide(2, 1);
  c2D->cd(1);
  h2DSign->SetMinimum(0);
  h2DSign->Draw("lego2");
  c2D->cd(2);
  h2DSideb->SetMinimum(0);
  h2DSideb->Draw("lego2");
  c2D->SaveAs(Form("Output_CorrelationExtraction_%s_png/h2D_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
  c2D->SaveAs(Form("Output_CorrelationExtraction_%s_Root/h2D_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));

  // Get FD correlations for FD subtraction
  if (fFDsubtraction) {
    h2DFdTemplatePrompt = getFdTemplateHisto(kPrompt, ptCandMin, ptCandMax, ptHadMin, ptHadMax);
    h2DFdTemplateNonPrompt = getFdTemplateHisto(kFD, ptCandMin, ptCandMax, ptHadMin, ptHadMax);
    // h1D_BaselineSubtr
    fdPromptFrac = getFdPromptFrac(ptCandMin, ptCandMax, ptHadMin, ptHadMax);

    if (fRebinFDCorr) {
      h2DFdTemplatePrompt->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi);
      h2DFdTemplateNonPrompt->Rebin2D(fRebinAxisDeltaEta, fRebinAxisDeltaPhi);
    }

    if (fDebug >= 1) {
      auto* c = new TCanvas(Form("cFDTemplate_PtCand%.0fto%.0f_PtHad%.0fto%.0f", ptCandMin, ptCandMax, ptHadMin, ptHadMax), Form("cFDTemplate_%s_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
      c->Divide(2, 1);
      c->cd(1);
      h2DFdTemplatePrompt->SetMinimum(0);
      h2DFdTemplatePrompt->Draw("lego2");
      c->cd(2);
      h2DFdTemplateNonPrompt->SetMinimum(0);
      h2DFdTemplateNonPrompt->Draw("lego2");
      c->SaveAs(Form("Output_CorrelationExtraction_%s_png/CorrFDTemplate_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
    }
  }

  // Bkg subtraction (2D plot)
  auto* c2DSub = new TCanvas(Form("c2D_Subtr_IntPools_PtHAd%.0fto%.0f", ptHadMin, ptHadMax), Form("c2D_%s_Subtr_IntPools_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), ptHadMin, ptHadMax), 100, 100, 1500, 800);
  h2DSubtr = reinterpret_cast<TH2D*>(h2DSign->Clone("h2D_Subtr"));
  h2DSubtr->Sumw2();
  h2DSubtr->Add(h2DSideb, -1);
  h2DSubtr->SetEntries(h2DSign->GetEntries() - h2DSideb->GetEntries());
  h2DSubtr->SetTitle("Signal region after sideb. subt. corr. - 2D");
  h2DSubtr->Draw("lego2");
  c2DSub->SaveAs(Form("Output_CorrelationExtraction_%s_png/h2D_%s_Subtr_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
  c2DSub->SaveAs(Form("Output_CorrelationExtraction_%s_Root/h2D_%s_Subtr_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));

  // 1D projection
  h1DSign = h2DSign->ProjectionY("h1D_Sign"); // projection on deltaPhi axis
  h1DSideb = h2DSideb->ProjectionY("h1D_Sideb");
  h1DSign->SetTitle("Signal region correlations");
  h1DSideb->SetTitle("Sidebands correlations");
  h1DSign->Scale(1. / h1DSign->GetXaxis()->GetBinWidth(1));
  h1DSideb->Scale(1. / h1DSideb->GetXaxis()->GetBinWidth(1));

  // Bkg subtraction (1D plot)
  h1DSubtr = reinterpret_cast<TH1D*>(h1DSign->Clone("h1D_Subtr"));
  h1DSubtr->Sumw2();
  h1DSubtr->Add(h1DSideb, -1);
  h1DSubtr->SetEntries(h1DSign->GetEntries() - h1DSideb->GetEntries());
  h1DSubtr->SetTitle("Signal region after sideb. subt. corr.");

  // Draw 1D plots (Signal region, Sidebands, S-SB (subtr.))
  auto* c1D = new TCanvas(Form("c1D_IntPools_PtCand%.0fto%.0f_PtHad%.0fto%.0f", ptCandMin, ptCandMax, ptHadMin, ptHadMax), Form("c1D_%s_IntPools_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), ptHadMin, ptHadMax), 100, 100, 1600, 500);
  c1D->Divide(3, 1);
  c1D->cd(1);
  h1DSign->Draw();
  c1D->cd(2);
  h1DSideb->Draw();
  c1D->cd(3);
  h1DSubtr->Draw();
  c1D->SaveAs(Form("Output_CorrelationExtraction_%s_png/h1D_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
  c1D->SaveAs(Form("Output_CorrelationExtraction_%s_Root/h1D_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));

  if (fDebug >= 1) {
    h1DSignNorm = reinterpret_cast<TH1D*>(h1DSign->Clone("h1D_Sign_Norm"));
    h1DSidebNorm = reinterpret_cast<TH1D*>(h1DSideb->Clone("h1D_Sideb_Norm"));
    h1DSignNorm->Scale(1. / (fSgnYieldNorm + fBkgYield));
    // h1D_SidebNorm -> Scale(1./fBkgYield);
    h1DSidebNorm->Scale(1. / fBkgScaleFactor);
    h1DSidebNorm->Scale(1. / fSBYield);
    h1DSignNorm->SetMarkerStyle(kFullCircle);
    h1DSignNorm->SetMarkerSize(1.2);
    h1DSignNorm->SetLineColor(kRed);
    h1DSignNorm->SetMarkerColor(kRed);
    h1DSignNorm->SetLineWidth(2);
    h1DSidebNorm->SetMinimum(0);
    h1DSidebNorm->SetMarkerStyle(kFullSquare);
    h1DSidebNorm->SetMarkerSize(1.2);
    h1DSidebNorm->SetLineColor(kBlue);
    h1DSidebNorm->SetMarkerColor(kBlue);
    h1DSidebNorm->SetLineWidth(2);
    h1DSidebNorm->SetTitle(Form("%.0f < p_{T} < %.0f", ptCandMin, ptCandMax));
    auto* c = new TCanvas(Form("c_IntPools_PtCand%.0fto%.0f_PtHad%.0fto%.0f", ptCandMin, ptCandMax, ptHadMin, ptHadMax), "");
    c->cd();
    h1DSidebNorm->Draw();
    h1DSignNorm->Draw("same");
    c->SaveAs(Form("Output_CorrelationExtraction_%s_png/ComparisonSignalSidebCorr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
    c->SaveAs(Form("Output_CorrelationExtraction_%s_Root/ComparisonSignalSidebCorr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
  }
  // Apply normalization to number of triggers
  h1DSubtrNorm = reinterpret_cast<TH1D*>(h1DSubtr->Clone("h1D_SubtrNorm"));
  h1DSubtrNorm->Sumw2();
  h1DSubtrNorm->Scale(1. / fSgnYieldNorm);
  h1DSubtrNorm->SetTitle("Signal region after sideb. subt. corr. - Normalized to # of triggers");

  // Correction for bias B to D topologies
  if (fCorrBiasBtoD) {
    hModul = evaluateMcClosModulations(ptCandMin, ptCandMax, ptHadMin, ptHadMax);
    auto* c1DCorrBbias = new TCanvas(Form("c1D_corrBbias_IntPools_PtCand%.0fto%.0f_PtHad%.0fto%.0f", ptCandMin, ptCandMax, ptHadMin, ptHadMax), Form("c1D_corrBbias_%s_IntPools_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), ptHadMin, ptHadMax), 100, 100, 1600, 500);
    c1DCorrBbias->cd();
    hBeforeModulCorr = reinterpret_cast<TH1D*>(h1DSubtrNorm->Clone("hBeforeModulCorr"));
    hBeforeModulCorr->SetLineColor(kViolet - 3);
    hBeforeModulCorr->GetYaxis()->SetRangeUser(0., 5.);
    hBeforeModulCorr->Draw();
    h1DSubtrNorm->Multiply(hModul);
    h1DSubtrNorm->Draw("same");
    c1DCorrBbias->SaveAs(Form("Output_CorrelationExtraction_%s_png/ComparisonCorrBiasBtoD_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
    c1DCorrBbias->SaveAs(Form("Output_CorrelationExtraction_%s_Root/ComparisonCorrBiasBtoD_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));

    auto* file = new TFile(Form("Output_CorrelationExtraction_%s_Root/SystematicCorrBiasBtoD_%s_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax), "RECREATE"); // Open file in write mode
    TH1D* h1DSubtrNormClone = reinterpret_cast<TH1D*>(h1DSubtrNorm->Clone("h1D_SubtrNorm_Clone"));
    h1DSubtrNormClone = reflectCorrHistogram(h1DSubtrNormClone);
    hBeforeModulCorr = reflectCorrHistogram(hBeforeModulCorr);
    TH1D* hSystematicCorrBiasBtoD = reinterpret_cast<TH1D*>(h1DSubtrNormClone->Clone("hSystematicCorrBiasBtoD"));
    hSystematicCorrBiasBtoD->Add(h1DSubtrNormClone, hBeforeModulCorr, 1, -1);
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
    h1DPrimaryPartCorr = getCorrelHistoSecondaryPart(kPrimaryPart, ptCandMin, ptCandMax, ptHadMin, ptHadMax);
    h1DAllPartCorr = getCorrelHistoSecondaryPart(kAllPart, ptCandMin, ptCandMax, ptHadMin, ptHadMax);
    h1DPrimaryPartCorr->Sumw2();
    h1DAllPartCorr->Sumw2();
    if (fRebinSecPart) {
      h1DPrimaryPartCorr->RebinX(fRebinAxisDeltaPhi); // Xaxis: deltaPhi
      h1DAllPartCorr->RebinX(fRebinAxisDeltaPhi);
      std::cout << "Secondary particle histogram rebinned" << std::endl;
    }
    h1DSecPartFrac = reinterpret_cast<TH1D*>(h1DPrimaryPartCorr->Clone(Form("hCorrRatio_PtD%.0fto%.0f_PtHad%.0fto%.0f", ptCandMin, ptCandMax, ptHadMin, ptHadMax)));
    h1DSecPartFrac->Sumw2();
    h1DSecPartFrac->Divide(h1DPrimaryPartCorr, h1DAllPartCorr, 1., 1., "B");

    auto* c1D = new TCanvas(Form("c1D_CorrPrimaryPart_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", ptCandMin, ptCandMax, ptHadMin, ptHadMax), Form("c1D_%s_CorrPrimaryPart_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), ptHadMin, ptHadMax));
    c1D->cd();
    setTH1HistoStyle(h1DSecPartFrac, Form("%.0f < p_{T} < %.0f GeV/c", ptCandMin, ptCandMax), "#Delta#phi [rad]", "#frac{primary part.}{part. selected}");
    h1DSecPartFrac->Draw();
    c1D->SaveAs(Form("Output_CorrelationExtraction_%s_png/CorrPrimaryPartRatio_%s_Canvas_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
    c1D->SaveAs(Form("Output_CorrelationExtraction_%s_Root/CorrPrimaryPartRatio_%s_Canvas_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));

    h1DSubtrNormSecPart = reinterpret_cast<TH1D*>(h1DSubtrNorm->Clone("h1D_SubtrNorm_SecPart"));
    h1DSubtrNormSecPart->Sumw2();
    Int_t const nBinsPhi = h1DSubtrNormSecPart->GetNbinsX();
    if (nBinsPhi != h1DSecPartFrac->GetNbinsX()) {
      std::cout << "[ERROR]: nBinsPhi different between h1D_SubtrNorm and h1D_SecPartFrac" << std::endl;
      return kFALSE;
    }
    h1DSubtrNormSecPart->Multiply(h1DSecPartFrac);
  }

  // FD Subtraction
  if (fFDsubtraction) {
    h1DFdTemplatePrompt = h2DFdTemplatePrompt->ProjectionY("h1D_FDTemplatePrompt");
    h1DFdTemplateNonPrompt = h2DFdTemplateNonPrompt->ProjectionY("h1D_FDTemplateNonPrompt");

    h1DFdTemplatePrompt->Scale(1. / h1DFdTemplatePrompt->GetXaxis()->GetBinWidth(1));
    h1DFdTemplateNonPrompt->Scale(1. / h1DFdTemplateNonPrompt->GetXaxis()->GetBinWidth(1));

    h1DTemplateTotal = reinterpret_cast<TH1D*>(h1DFdTemplatePrompt->Clone("h1D_TemplateTotal"));
    h1DTemplateTotal->Sumw2();
    h1DTemplateTotal->Scale(fdPromptFrac);
    h1DTemplateTotal->Add(h1DFdTemplateNonPrompt, 1 - fdPromptFrac);

    if (fDebug >= 1) {
      auto* c = new TCanvas(Form("cFDTemplate_1D_%1.1fto%1.1f", ptHadMin, ptHadMax), Form("cFDTemplate_%s_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
      c->cd();
      h1DTemplateTotal->SetMinimum(0);
      h1DFdTemplateNonPrompt->SetMinimum(0);
      h1DTemplateTotal->SetMarkerColor(kGreen);
      h1DTemplateTotal->SetLineColor(kGreen);
      h1DTemplateTotal->SetLineWidth(2);
      h1DTemplateTotal->SetMarkerStyle(kFullCircle);
      h1DFdTemplatePrompt->SetMarkerColor(kRed);
      h1DFdTemplatePrompt->SetLineColor(kRed);
      h1DFdTemplatePrompt->SetLineWidth(2);
      h1DFdTemplatePrompt->SetMarkerStyle(kFullCircle);
      h1DFdTemplateNonPrompt->SetMarkerColor(kBlue);
      h1DFdTemplateNonPrompt->SetLineColor(kBlue);
      h1DFdTemplateNonPrompt->SetLineWidth(2);
      h1DFdTemplateNonPrompt->SetMarkerStyle(kFullCircle);
      h1DFdTemplateNonPrompt->Draw();
      h1DFdTemplatePrompt->Draw("same");
      h1DTemplateTotal->Draw("same");
      auto* lFD = new TLegend();
      lFD->AddEntry(h1DTemplateTotal, "Total template");
      lFD->AddEntry(h1DFdTemplatePrompt, "Prompt Template");
      lFD->AddEntry(h1DFdTemplateNonPrompt, "Non prompt template");
      lFD->Draw("same");
      c->SaveAs(Form("Output_CorrelationExtraction_%s_png/CorrFDTemplate_1D_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
      c->SaveAs(Form("Output_CorrelationExtraction_%s_Root/CorrFDTemplate_1D_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
    }

    Double_t const baselineFd = calculateBaseline(h1DTemplateTotal, kTRUE);
    Double_t baselineData;
    if (fSecPartContamination) {
      baselineData = calculateBaseline(h1DSubtrNormSecPart, kTRUE);
    } else {
      baselineData = calculateBaseline(h1DSubtrNorm, kTRUE);
    }

    std::cout << "===================== " << std::endl;
    std::cout << "Baseline FD: " << baselineFd << std::endl;
    std::cout << "Baseline Data: " << baselineData << std::endl;
    std::cout << "===================== " << std::endl;
    std::cout << " " << std::endl;

    Double_t const baselinediff = baselineData - baselineFd;
    TH1D* hBaselineDiff = reinterpret_cast<TH1D*>(h1DFdTemplateNonPrompt->Clone("hBaselineDiff"));
    for (int iBin = 0; iBin < hBaselineDiff->GetNbinsX(); iBin++) {
      hBaselineDiff->SetBinContent(iBin + 1, baselinediff);
    }
    h1DFdTemplateNonPrompt->Add(hBaselineDiff);
    h1DTemplateTotal->Add(hBaselineDiff);
    if (fSecPartContamination) {
      h1DSubtrFdNorm = reinterpret_cast<TH1D*>(h1DSubtrNormSecPart->Clone("h1D_SubtrFDNorm"));
    } else {
      h1DSubtrFdNorm = reinterpret_cast<TH1D*>(h1DSubtrNorm->Clone("h1D_SubtrFDNorm"));
    }
    h1DFdTemplateNonPrompt->Scale(1 - fdPromptFrac);
    h1DSubtrFdNorm->Add(h1DFdTemplateNonPrompt, -1);
    h1DSubtrFdNorm->Scale(1. / fdPromptFrac);

    if (fDebug >= 1) {
      auto* c1 = new TCanvas(Form("cFDTemplateSubtr_%1.1fto%1.1f", ptHadMin, ptHadMax), Form("cFDTemplateSubtr_%s_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax), 100, 100, 1600, 900);
      c1->cd();
      h1DSubtrNorm->SetLineColor(kRed);
      h1DSubtrNormSecPart->SetLineColor(kOrange);
      h1DFdTemplateNonPrompt->SetLineColor(kBlue);
      h1DSubtrFdNorm->SetLineColor(kGreen);
      h1DTemplateTotal->SetLineColor(kMagenta);
      h1DSubtrNorm->SetMinimum(0);
      h1DSubtrNormSecPart->SetMinimum(0);
      h1DFdTemplateNonPrompt->SetMinimum(0);
      h1DSubtrFdNorm->SetMinimum(0);
      // h1D_SubtrNorm -> GetYaxis() -> SetRangeUser(0., 8.);
      h1DSubtrNorm->SetMarkerStyle(kFullCircle);
      h1DSubtrNorm->SetMarkerSize(1.2);
      h1DSubtrNorm->SetMarkerColor(kRed);
      h1DSubtrNorm->SetLineWidth(2);
      h1DSubtrNormSecPart->SetMarkerStyle(kFullCircle);
      h1DSubtrNormSecPart->SetMarkerSize(1.2);
      h1DSubtrNormSecPart->SetMarkerColor(kOrange);
      h1DSubtrNormSecPart->SetLineWidth(2);
      h1DSubtrFdNorm->SetMarkerStyle(kFullCircle);
      h1DSubtrFdNorm->SetMarkerSize(1.2);
      h1DSubtrFdNorm->SetMarkerColor(kGreen);
      h1DSubtrFdNorm->SetLineWidth(2);
      h1DSubtrNorm->GetYaxis()->SetTitle("#frac{1}{N_{D}} #frac{dN^{assoc. part}}{d#Delta#phi}");
      h1DSubtrNormSecPart->GetYaxis()->SetTitle("#frac{1}{N_{D}} #frac{dN^{assoc. part}}{d#Delta#phi}");
      if (fSecPartContamination) {
        h1DSubtrNormSecPart->Draw();
      } else {
        h1DSubtrNorm->Draw();
      }
      // h1D_FDTemplateNonPrompt -> Draw("same");
      h1DSubtrFdNorm->Draw("same");
      h1DTemplateTotal->Draw("same");
      c1->SaveAs(Form("Output_CorrelationExtraction_%s_png/CorrFDTemplateSubtr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
      c1->SaveAs(Form("Output_CorrelationExtraction_%s_Root/CorrFDTemplateSubtr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
    }
  }

  if (fFDsubtraction) {
    fCorrectedCorrHisto = reinterpret_cast<TH1D*>(h1DSubtrFdNorm->Clone(Form("hCorrectedCorr_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", ptCandMin, ptCandMax, ptHadMin, ptHadMax)));
  } else if (fSecPartContamination) {
    fCorrectedCorrHisto = reinterpret_cast<TH1D*>(h1DSubtrNormSecPart->Clone(Form("hCorrectedCorr_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", ptCandMin, ptCandMax, ptHadMin, ptHadMax)));
  } else {
    fCorrectedCorrHisto = reinterpret_cast<TH1D*>(h1DSubtrNorm->Clone(Form("hCorrectedCorr_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", ptCandMin, ptCandMax, ptHadMin, ptHadMax)));
  }

  std::cout << "Analysis steps completed - baseline subtraction missing" << std::endl;

  // Draw 1D plots (Signal region, normalized)
  auto* cFinal = new TCanvas(Form("cFinal_%.0fto%.0f", ptHadMin, ptHadMax), Form("cFinal_%s_IntPools_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), ptHadMin, ptHadMax), 100, 100, 1200, 700);
  h1DSubtrNorm->SetLineColor(kBlue + 1);
  h1DSubtrNorm->SetMarkerColor(kBlue + 1);
  h1DSubtrNorm->SetMarkerStyle(kFullCircle);
  h1DSubtrNorm->SetMinimum(0);
  h1DSubtrNorm->Draw();
  if (fSecPartContamination) {
    h1DSubtrNormSecPart->SetLineColor(kRed + 1);
    h1DSubtrNormSecPart->SetMarkerColor(kRed + 1);
    h1DSubtrNormSecPart->SetMarkerStyle(kFullCircle);
    h1DSubtrNormSecPart->Draw("same");
  }
  if (fFDsubtraction) {
    h1DSubtrFdNorm->SetLineColor(kGreen + 2);
    h1DSubtrFdNorm->SetMarkerColor(kGreen + 2);
    h1DSubtrFdNorm->SetMarkerStyle(kFullCircle);
    h1DSubtrFdNorm->Draw("same");
  }
  if (fFDsubtraction) {
    h1DTemplateTotal->Draw("same");
  }
  auto* lFinal = new TLegend();
  lFinal->AddEntry(h1DSubtrNorm, "Corr. after bkg subtr.");
  if (fFDsubtraction) {
    lFinal->AddEntry(h1DTemplateTotal, "CR Mode 2 total template");
  }
  if (fSecPartContamination) {
    lFinal->AddEntry(h1DSubtrNormSecPart, "Corr. after sec. part. correction");
  }
  if (fFDsubtraction) {
    lFinal->AddEntry(h1DSubtrFdNorm, "Corr. FD subtr.");
  }
  lFinal->Draw("same");
  cFinal->SaveAs(Form("Output_CorrelationExtraction_%s_png/AzimCorrDistr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
  cFinal->SaveAs(Form("Output_CorrelationExtraction_%s_Root/AzimCorrDistr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));

  // Baseline subtraction
  Double_t baselineData, baselineDataErr;
  TH1D* hBaseline = reinterpret_cast<TH1D*>(h1DSubtrNorm->Clone("hBaseline"));
  hBaseline->Sumw2();
  if (fFDsubtraction) {
    baselineData = calculateBaseline(h1DSubtrFdNorm, kTRUE, kFALSE); // introduced kFALSE
    baselineDataErr = calculateBaselineError(h1DSubtrFdNorm, kTRUE, kFALSE);
    for (int iBin = 0; iBin < hBaseline->GetNbinsX(); iBin++) {
      hBaseline->SetBinContent(iBin + 1, baselineData);
      hBaseline->SetBinError(iBin + 1, baselineDataErr);
    }
    h1DBaselineSubtr = reinterpret_cast<TH1D*>(h1DSubtrFdNorm->Clone("h1D_BaselineSubtr"));
    h1DBaselineSubtr->Add(hBaseline, -1.);
  } else if (fSecPartContamination) {
    baselineData = calculateBaseline(h1DSubtrNormSecPart, kTRUE, kFALSE);
    baselineDataErr = calculateBaselineError(h1DSubtrNormSecPart, kTRUE, kFALSE);
    for (int iBin = 0; iBin < hBaseline->GetNbinsX(); iBin++) {
      hBaseline->SetBinContent(iBin + 1, baselineData);
      hBaseline->SetBinError(iBin + 1, baselineDataErr);
    }
    h1DBaselineSubtr = reinterpret_cast<TH1D*>(h1DSubtrNormSecPart->Clone("h1D_BaselineSubtr"));
    h1DBaselineSubtr->Add(hBaseline, -1.);
  } else {
    baselineData = calculateBaseline(h1DSubtrNorm, kTRUE, kFALSE);
    baselineDataErr = calculateBaselineError(h1DSubtrNorm, kTRUE, kFALSE);
    for (int iBin = 0; iBin < hBaseline->GetNbinsX(); iBin++) {
      hBaseline->SetBinContent(iBin + 1, baselineData);
      hBaseline->SetBinError(iBin + 1, baselineDataErr);
    }
    h1DBaselineSubtr = reinterpret_cast<TH1D*>(h1DSubtrNorm->Clone("h1D_BaselineSubtr"));
    h1DBaselineSubtr->Add(hBaseline, -1.);
  }

  fCorrectedCorrHistoBaselineSubtr = reinterpret_cast<TH1D*>(h1DBaselineSubtr->Clone(Form("hCorrectedCorrBaselineSubtr_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", ptCandMin, ptCandMax, ptHadMin, ptHadMax)));

  auto* cFinalBaselineSubtr = new TCanvas(Form("cFinal_BaselineSubtr_%.0fto%.0f", ptHadMin, ptHadMax), Form("cFinal_BaselineSubtr_%s_IntPools_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), ptHadMin, ptHadMax), 100, 100, 1200, 700);
  h1DBaselineSubtr->SetMarkerColor(kOrange + 8);
  h1DBaselineSubtr->SetLineColor(kOrange + 8);
  h1DBaselineSubtr->GetYaxis()->SetRangeUser(-0.2, 8.);
  h1DBaselineSubtr->Draw();
  if (fFDsubtraction) {
    h1DSubtrFdNorm->Draw("same");
  } else if (fSecPartContamination) {
    h1DSubtrNormSecPart->Draw("same");
  } else {
    h1DSubtrNorm->Draw("same");
  }
  hBaseline->SetMarkerColor(kPink - 6);
  hBaseline->SetMarkerStyle(kFullSquare);
  hBaseline->SetLineColor(kPink - 6);
  hBaseline->Draw("same");

  cFinalBaselineSubtr->SaveAs(Form("Output_CorrelationExtraction_%s_png/AzimCorrDistr_BaselineSubtr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
  cFinalBaselineSubtr->SaveAs(Form("Output_CorrelationExtraction_%s_Root/AzimCorrDistr_BaselineSubtr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));

  // Reflected histograms
  if (fFDsubtraction) {
    h1DReflCorr = reflectCorrHistogram(h1DSubtrFdNorm);
  } else if (fSecPartContamination) {
    h1DReflCorr = reflectCorrHistogram(h1DSubtrNormSecPart);
  } else {
    h1DReflCorr = reflectCorrHistogram(h1DSubtrNorm);
  }

  /* used as control using Run2 reflection function
  if (fFDsubtraction) {
    h1D_ReflCorr = ReflectHistoRun2(h1D_SubtrFDNorm, 0.5);
  } else if (fSecPartContamination) {
    h1D_ReflCorr = ReflectHistoRun2(h1D_SubtrNorm_SecPart, 0.5);
  } else {
    h1D_ReflCorr = ReflectHistoRun2(h1D_SubtrNorm, 0.5);
  }*/

  auto* cFinalReflected = new TCanvas(Form("cFinal_Reflected_%.0fto%.0f", ptHadMin, ptHadMax), Form("cFinal_Reflected_%s_IntPools_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), ptHadMin, ptHadMax), 100, 100, 1200, 700);
  cFinalReflected->cd();
  setTH1HistoStyle(h1DReflCorr, Form("%.0f < p_{T} < %.0f GeV/c", ptCandMin, ptCandMax), "#Delta#phi [rad]", "#frac{1}{N_{D}}#frac{dN^{assoc}}{d#Delta#phi} [rad^{-1}]", kFullCircle, kOrange + 8, 1.6, kOrange + 8, 3);
  h1DReflCorr->SetMinimum(0);
  h1DReflCorr->Draw();
  cFinalReflected->SaveAs(Form("Output_CorrelationExtraction_%s_png/AzimCorrDistr_Reflected_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
  cFinalReflected->SaveAs(Form("Output_CorrelationExtraction_%s_Root/AzimCorrDistr_Reflected_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));

  // Reflected histograms baseline subtracted
  TH1D* hBaselineRefl = reinterpret_cast<TH1D*>(h1DReflCorr->Clone("hBaseline_Refl"));
  hBaselineRefl->Sumw2();
  baselineData = calculateBaseline(h1DReflCorr, kFALSE, kTRUE);
  baselineDataErr = calculateBaselineError(h1DReflCorr, kFALSE, kTRUE);

  for (int iBin = 0; iBin < hBaselineRefl->GetNbinsX(); iBin++) {
    hBaselineRefl->SetBinContent(iBin + 1, baselineData);
    hBaselineRefl->SetBinError(iBin + 1, baselineDataErr);
  }
  h1DReflCorrBaselineSubtr = reinterpret_cast<TH1D*>(h1DReflCorr->Clone("h1D_ReflCorr_BaselineSubtr"));
  h1DReflCorrBaselineSubtr->Sumw2();
  h1DReflCorrBaselineSubtr->Add(hBaselineRefl, -1.);

  TF1* fConstZero = new TF1("fConstZero", "[0]", 0., TMath::Pi());
  fConstZero->SetParameter(0, 0.);
  fConstZero->SetLineColor(kMagenta);
  fConstZero->SetLineStyle(9);
  fConstZero->SetLineWidth(4);
  fConstZero->SetTitle("");

  auto* cFinalReflectedBaselineSubtr = new TCanvas(Form("cFinal_Reflected_BaselineSubtr_%.0fto%.0f", ptHadMin, ptHadMax), Form("cFinal_Reflected_BaselineSubtr_%s_IntPools_PtAssoc%.0fto%.0f", fDmesonLabel.Data(), ptHadMin, ptHadMax), 100, 100, 1200, 700);
  setTH1HistoStyle(h1DReflCorrBaselineSubtr, Form("%.0f < p_{T} < %.0f GeV/c", ptCandMin, ptCandMax), "#Delta#phi [rad]", "#frac{1}{N_{D}}#frac{dN^{assoc}}{d#Delta#phi} [rad^{-1}]", kFullCircle, kRed + 1, 1.6, kRed + 1, 3);
  hBaselineRefl->SetMarkerColor(kOrange);
  hBaselineRefl->SetMarkerStyle(kFullSquare);
  hBaselineRefl->SetLineColor(kOrange);
  cFinalReflectedBaselineSubtr->cd();
  h1DReflCorr->SetMinimum(-0.8);
  h1DReflCorr->SetStats(false);
  hBaselineRefl->SetStats(false);
  h1DReflCorr->Draw();
  hBaselineRefl->Draw("same");
  h1DReflCorrBaselineSubtr->SetStats(false);
  h1DReflCorrBaselineSubtr->Draw("same"); // then keep just this
  fConstZero->Draw("same");
  cFinalReflectedBaselineSubtr->SaveAs(Form("Output_CorrelationExtraction_%s_png/AzimCorrDistr_Reflected_BaselineSubtr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));
  cFinalReflectedBaselineSubtr->SaveAs(Form("Output_CorrelationExtraction_%s_Root/AzimCorrDistr_Reflected_BaselineSubtr_%s_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.root", codeName.Data(), fDmesonLabel.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax));

  fCorrectedCorrHistoReflected = reinterpret_cast<TH1D*>(h1DReflCorr->Clone(Form("hCorrectedCorrReflected_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", ptCandMin, ptCandMax, ptHadMin, ptHadMax)));
  fCorrectedCorrHistoReflectedBaselineSubtr = reinterpret_cast<TH1D*>(h1DReflCorrBaselineSubtr->Clone(Form("hCorrectedCorrReflected_BaselineSubtr_PtCand%.0fto%.0f_PtAssoc%.0fto%.0f", ptCandMin, ptCandMax, ptHadMin, ptHadMax)));

  return kTRUE;
}

Bool_t DhCorrelationExtraction::readInputSeAndMe()
{

  fFileSE = TFile::Open(fFileNameSE.Data());
  if (fFileSE == nullptr) {
    std::cout << "[ERROR] File " << fFileNameSE << " cannot be opened! check your file path!";
    return kFALSE;
  }

  fFileME = TFile::Open(fFileNameME.Data());
  if (fFileME == nullptr) {
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

Bool_t DhCorrelationExtraction::readInputInvMass()
{

  fFileMass = TFile::Open(fFileNameMass.Data());
  if (fFileMass == nullptr) {
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

Bool_t DhCorrelationExtraction::readInputFdSubtr()
{

  fFileFDTemplate = TFile::Open(fFileFDTemplateName.Data());
  fFileFDPromptFrac = TFile::Open(fFileFDPromptFracName.Data());
  if (fFileFDTemplate == nullptr) {
    std::cout << "[ERROR] File " << fFileFDTemplateName << " cannot be opened! check your file path!" << std::endl;
    return kFALSE;
  }
  if (fFileFDPromptFrac == nullptr) {
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

Bool_t DhCorrelationExtraction::readInputSecondaryPartContamination()
{

  fFileSecPart = TFile::Open(fFileSecPartName.Data());
  if (fFileSecPart == nullptr) {
    std::cout << "[ERROR] File " << fFileSecPartName << " cannot be opened! check your file path!" << std::endl;
    return kFALSE;
  }

  fDirSecPart = reinterpret_cast<TDirectoryFile*>(fFileSecPart->Get(fDirSecPartName.Data()));

  if (fDirSecPart == nullptr) {
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

TH1D* DhCorrelationExtraction::evaluateMcClosModulations(Double_t ptCandMin, Double_t ptCandMax, Double_t ptHadMin, Double_t ptHadMax)
{

  TH1D* hModul = nullptr;

  fFilePromptMc = TFile::Open(fFilePromptMcRecName.Data());
  if (fFilePromptMc == nullptr) {
    throw std::runtime_error("[ERROR] File prompt MC rec cannot be opened! check your file path!");
  }
  fFileNonPromptMc = TFile::Open(fFileNonPromptMcRecName.Data());
  if (fFileNonPromptMc == nullptr) {
    throw std::runtime_error("[ERROR] File non-prompt MC rec cannot be opened! check your file path!");
  }

  // TODO: generalise this part
  TH1D* hRecPrompt = reinterpret_cast<TH1D*>(fFilePromptMc->Get(Form("h1D_Rec_iPtD%d_iPtAssoc%d", fBinPtCand, fBinPtHad)));
  TH1D* hRecNonPrompt = reinterpret_cast<TH1D*>(fFileNonPromptMc->Get(Form("h1D_Rec_iPtD%d_iPtAssoc%d", fBinPtCand, fBinPtHad)));
  TH1D* hGenPrompt = reinterpret_cast<TH1D*>(fFilePromptMc->Get(Form("h1D_Gen_iPtD%d_iPtAssoc%d", fBinPtCand, fBinPtHad)));
  TH1D* hGenNonPrompt = reinterpret_cast<TH1D*>(fFileNonPromptMc->Get(Form("h1D_Gen_iPtD%d_iPtAssoc%d", fBinPtCand, fBinPtHad)));

  printf("[INFO] Bin cand %d - Bin had %d \n", fBinPtCand, fBinPtHad);

  // hRecPrompt = reflectCorrHistogram(hRecPrompt);
  // hRecNonPrompt = reflectCorrHistogram(hRecNonPrompt);
  // hGenPrompt = reflectCorrHistogram(hGenPrompt);
  // hGenNonPrompt = reflectCorrHistogram(hGenNonPrompt);

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
  Double_t const fitVal = funFit->GetParameter(0);

  auto* cRatioMcClosure = new TCanvas(Form("cRatio_MCClosure_PtCand%.0fto%.0f_Pthad%.0fto%.0f", ptCandMin, ptCandMax, ptHadMin, ptHadMax), Form("cRatio_MCClosure_PtCand%.0fto%.0f_Pthad%.0fto%.0f", ptCandMin, ptCandMax, ptHadMin, ptHadMax), 100, 100, 1200, 700);
  cRatioMcClosure->cd();
  hRatioNonPrompt->GetYaxis()->SetRangeUser(0.2, 1.8);
  hRatioNonPrompt->Draw();

  Double_t const fPrompt = getFdPromptFrac(ptCandMin, ptCandMax, ptHadMin, ptHadMax);
  Double_t relAmplC[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t relAmplB[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t recoKineVal[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t modul[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  for (int iBin = 0; iBin < hRatioNonPrompt->GetNbinsX(); iBin++) {
    if (iBin > 1 && iBin < 13) {
      recoKineVal[iBin - 2] = hRatioNonPrompt->GetBinContent(iBin + 1) - (fitVal - 1);
      relAmplC[iBin - 2] = hRecPrompt->GetBinContent(iBin + 1) / (hRecPrompt->GetBinContent(iBin + 1) * fPrompt + hRecNonPrompt->GetBinContent(iBin + 1) * (1 - fPrompt));
      relAmplB[iBin - 2] = hRecNonPrompt->GetBinContent(iBin + 1) / (hRecPrompt->GetBinContent(iBin + 1) * fPrompt + hRecNonPrompt->GetBinContent(iBin + 1) * (1 - fPrompt));
      modul[iBin - 2] = relAmplC[iBin - 2] * fPrompt + relAmplB[iBin - 2] * (1 - fPrompt) / recoKineVal[iBin - 2];
      hModul->SetBinContent(iBin + 1, modul[iBin - 2]);
      hModul->SetBinError(iBin + 1, 0.);

      printf("[INFO] Bin%d MODUL = %1.5f\t (Reco/Kine-fitVal = %1.4f, FPrompt = %1.3f, Ampl_ratio C,B = %1.4f, %1.4f)\n", iBin + 1, modul[iBin - 2], recoKineVal[iBin - 2], fPrompt, relAmplC[iBin - 2], relAmplB[iBin - 2]);
    } else {
      hModul->SetBinContent(iBin + 1, 1.);
      hModul->SetBinError(iBin + 1, 0.);
    }
  }

  hModul->SetLineColor(kMagenta);
  hModul->Draw("same");

  cRatioMcClosure->SaveAs(Form("Output_CorrelationExtraction_Thin2023_FullAnalysis_CentralPoints_png/Ratio_MCClosure_Canvas_PtCand%.0fto%.0f_PoolInt_PtAssoc%.0fto%.0f.png", ptCandMin, ptCandMax, ptHadMin, ptHadMax));

  return hModul;
}

TH2D* DhCorrelationExtraction::getCorrelHisto(Int_t sEorMe, Int_t sorSb, Int_t pool, Double_t ptCandMin, Double_t ptCandMax, Double_t ptHadMin, Double_t ptHadMax)
{
  // TODO: Subtraction of softpion
  TH2D* h2D = nullptr; // pointer to be returned

  THnSparseD* hSparse = nullptr;
  if (sEorMe == kSE) { // Same Event
    if (sorSb == kSign) {
      hSparse = reinterpret_cast<THnSparseD*>(fDirSE->Get(fSECorrelSignalRegionName.Data()));
    } else if (!fSidebandDivided) {
      hSparse = reinterpret_cast<THnSparseD*>(fDirSE->Get(fSECorrelSidebandsName.Data()));
    } else {
      if (fUseSidebLeft && !fUseSidebRight) {
        hSparse = reinterpret_cast<THnSparseD*>(fDirSE->Get(fSECorrelSidebandLeftName.Data()));
      } else if (!fUseSidebLeft && fUseSidebRight) {
        hSparse = reinterpret_cast<THnSparseD*>(fDirSE->Get(fSECorrelSidebandRightName.Data()));
      } else if (fUseSidebLeft && fUseSidebRight) {
        hSparse = reinterpret_cast<THnSparseD*>(fDirSE->Get(fSECorrelSidebandLeftName.Data()));
        hSparse->SetName("hSparse");
        auto* hSparseRightSideb = reinterpret_cast<THnSparseD*>(fDirSE->Get(fSECorrelSidebandRightName.Data()));
        hSparse->Add(hSparseRightSideb, 1.);
      }
    }
  } else { // Mixed Event
    if (sorSb == kSign) {
      hSparse = reinterpret_cast<THnSparseD*>(fDirME->Get(fMECorrelSignalRegionName.Data()));
    } else if (!fSidebandDivided) {
      hSparse = reinterpret_cast<THnSparseD*>(fDirME->Get(fMECorrelSidebandsName.Data()));
    } else {
      if (fUseSidebLeft && !fUseSidebRight) {
        hSparse = reinterpret_cast<THnSparseD*>(fDirME->Get(fMECorrelSidebandLeftName.Data()));
      } else if (!fUseSidebLeft && fUseSidebRight) {
        hSparse = reinterpret_cast<THnSparseD*>(fDirME->Get(fMECorrelSidebandRightName.Data()));
      } else if (fUseSidebLeft && fUseSidebRight) {
        hSparse = reinterpret_cast<THnSparseD*>(fDirME->Get(fMECorrelSidebandLeftName.Data()));
        hSparse->SetName("hSparse");
        auto* hSparseRightSideb = reinterpret_cast<THnSparseD*>(fDirME->Get(fMECorrelSidebandRightName.Data()));
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

  Int_t const binExtPtCandMin = hSparse->GetAxis(2)->FindBin(ptCandMin + 0.01); // axis2: ptCand, the 0.01 to avoid bin edges!
  Int_t const binExtPtCandMax = hSparse->GetAxis(2)->FindBin(ptCandMax - 0.01);
  Int_t const binExtPtHadMin = hSparse->GetAxis(3)->FindBin(ptHadMin + 0.01); // axis3: ptHad
  Int_t const binExtPtHadMax = hSparse->GetAxis(3)->FindBin(ptHadMax - 0.01);
  Int_t binExtPoolMin;
  Int_t binExtPoolMax;
  if (fCorrectPoolsSeparately) {
    binExtPoolMin = hSparse->GetAxis(4)->FindBin(pool + 0.01); // axis4: pool bin
    binExtPoolMax = hSparse->GetAxis(4)->FindBin(pool + 0.99);
  } else { // merge all pools in one
    binExtPoolMin = 1;
    binExtPoolMax = hSparse->GetAxis(4)->GetNbins();
  }
  // possibility to select a certain eta region
  Int_t binExtEtaMin = hSparse->GetAxis(1)->FindBin(fDeltaEtaMin + 0.0001);
  Int_t binExtEtaMax = hSparse->GetAxis(1)->FindBin(fDeltaEtaMax - 0.0001);
  if (binExtEtaMax > hSparse->GetAxis(1)->GetNbins()) {
    binExtEtaMax = hSparse->GetAxis(1)->GetNbins();
  }
  if (binExtEtaMin < 1) {
    binExtEtaMin = 1;
  }
  hSparse->GetAxis(1)->SetRange(binExtEtaMin, binExtEtaMax);       // axis1: deltaEta
  hSparse->GetAxis(2)->SetRange(binExtPtCandMin, binExtPtCandMax); // axis2: ptCand
  hSparse->GetAxis(3)->SetRange(binExtPtHadMin, binExtPtHadMax);   // axis3: ptHad
  hSparse->GetAxis(4)->SetRange(binExtPoolMin, binExtPoolMax);     // axis4: pool bin
  h2D = hSparse->Projection(0, 1);                                 // axis0: deltaPhi, axis1: deltaEta
  if (sEorMe == kSE) {                                             // Same Event
    if (sorSb == kSign) {
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
    if (sorSb == kSign) {
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

void DhCorrelationExtraction::getSignalAndBackgroundForNorm(Double_t ptCandMin, Double_t ptCandMax)
{

  // using results obtained from HFInvariantMassFitter.cxx class
  TH1F* hMassFitSgnYield = reinterpret_cast<TH1F*>(fFileMass->Get(fMassHistoNameSgn.Data()));
  TH1F* hMassFitBkgYield = reinterpret_cast<TH1F*>(fFileMass->Get(fMassHistoNameBkg.Data()));
  TH1F* hMassFitSBsYield = reinterpret_cast<TH1F*>(fFileMass->Get(fMassHistoNameSBs.Data()));
  TH1F* hMassFitSBLYield = reinterpret_cast<TH1F*>(fFileMass->Get("hBackgroundSidebandLeft"));
  TH1F* hMassFitSBRYield = reinterpret_cast<TH1F*>(fFileMass->Get("hBackgroundSidebandRight"));

  Int_t const ptCandBin = hMassFitSgnYield->FindBin(ptCandMin + 0.01);
  if (ptCandBin != hMassFitSgnYield->FindBin(ptCandMax - 0.01)) {
    std::cout << "[ERROR] Pt bin in invariant mass histogram not univocally defined " << std::endl;
  }

  Float_t const sgnYield = hMassFitSgnYield->GetBinContent(ptCandBin);
  Float_t const bkgYield = hMassFitBkgYield->GetBinContent(ptCandBin);
  Float_t const sBsYield = hMassFitSBsYield->GetBinContent(ptCandBin);
  Float_t const sblYield = hMassFitSBLYield->GetBinContent(ptCandBin);
  Float_t const sbrYield = hMassFitSBRYield->GetBinContent(ptCandBin);

  std::cout << "================================= " << std::endl;
  std::cout << "Getting invariant mass parameters " << std::endl;
  std::cout << "Pt cand " << ptCandMin << " - " << ptCandMax << std::endl;
  std::cout << "Signal yield    = " << sgnYield << std::endl;
  std::cout << "Bkg yield    = " << bkgYield << std::endl;
  std::cout << "Sideband yield    = " << sBsYield << std::endl;
  std::cout << "Sideband left yield    = " << sblYield << std::endl;
  std::cout << "Sideband right yield    = " << sbrYield << std::endl;
  std::cout << "================================= " << std::endl;
  std::cout << " " << std::endl;

  setSignalYieldforNorm(sgnYield);
  setBkgYield(bkgYield);
  if (fUseSidebLeft && fUseSidebRight) {
    setBkgScaleFactor(bkgYield / sBsYield);
    setSbYield(sBsYield);
  } else if (fUseSidebLeft && !fUseSidebRight) {
    setBkgScaleFactor(bkgYield / sblYield);
    setSbYield(sblYield);
  } else if (!fUseSidebLeft && fUseSidebRight) {
    setBkgScaleFactor(bkgYield / sbrYield);
    setSbYield(sbrYield);
  }
}

TH2D* DhCorrelationExtraction::getFdTemplateHisto(Int_t promptOrFd, Double_t ptCandMin, Double_t ptCandMax, Double_t ptHadMin, Double_t ptHadMax)
{

  TH2D* h2D = nullptr; // pointer to be returned

  if (promptOrFd == kPrompt) {
    h2D = reinterpret_cast<TH2D*>(fFileFDTemplate->Get(Form("%s%.0f_%.0f_ptassoc%.0f_%.0f", fHistoFDTemplatePromptName.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax)));
  } else {
    h2D = reinterpret_cast<TH2D*>(fFileFDTemplate->Get(Form("%s%.0f_%.0f_ptassoc%.0f_%.0f", fHistoFDTemplateNonPromptName.Data(), ptCandMin, ptCandMax, ptHadMin, ptHadMax)));
  }

  Int_t binExtEtaMin = h2D->GetXaxis()->FindBin(fDeltaEtaMin + 0.000001);
  Int_t binExtEtaMax = h2D->GetXaxis()->FindBin(fDeltaEtaMax - 0.000001);
  if (binExtEtaMax > h2D->GetXaxis()->GetNbins()) {
    binExtEtaMax = h2D->GetXaxis()->GetNbins();
  }
  if (binExtEtaMin < 1) {
    binExtEtaMin = 1;
  }

  h2D->GetXaxis()->SetRange(binExtEtaMin, binExtEtaMax);
  if (promptOrFd == kPrompt) {
    h2D->SetName(Form("hFDTemplatePrompt_2D_PtCand%.0fto%.0f_PtHad%.0fto%.0f", ptCandMin, ptCandMax, ptHadMin, ptHadMax));
  } else {
    h2D->SetName(Form("hFDTemplateNonPrompt_2D_PtCand%.0fto%.0f_PtHad%.0fto%.0f", ptCandMin, ptCandMax, ptHadMin, ptHadMax));
  }
  h2D->GetYaxis()->SetTitle("#Delta#phi (rad)");
  h2D->GetXaxis()->SetTitle("#Delta#eta");

  return h2D;
}

TH1D* DhCorrelationExtraction::getCorrelHistoSecondaryPart(Int_t partType, Double_t ptCandMin, Double_t ptCandMax, Double_t ptHadMin, Double_t ptHadMax)
{

  TH1D* h1D = nullptr; // pointer to be returned

  THnSparseD* hSparse = nullptr;

  if (partType == kPrimaryPart) { // primary particles
    hSparse = reinterpret_cast<THnSparseD*>(fDirSecPart->Get(fHistoPrimaryPartName.Data()));
  } else { // all selected particles
    hSparse = reinterpret_cast<THnSparseD*>(fDirSecPart->Get(fHistoAllPartName.Data()));
  }
  Int_t const binExtPtCandMin = hSparse->GetAxis(2)->FindBin(ptCandMin + 0.01); // axis2: ptCand, the 0.01 to avoid bin edges!
  Int_t const binExtPtCandMax = hSparse->GetAxis(2)->FindBin(ptCandMax - 0.01);
  Int_t const binExtPtHadMin = hSparse->GetAxis(3)->FindBin(ptHadMin + 0.01); // axis3: ptHad
  Int_t const binExtPtHadMax = hSparse->GetAxis(3)->FindBin(ptHadMax - 0.01);
  Int_t binExtPoolMin;
  Int_t binExtPoolMax;
  if (partType == kAllPart) {
    binExtPoolMin = 1;
    binExtPoolMax = hSparse->GetAxis(4)->GetNbins();
  }
  // possibility to select a certain eta region
  Int_t binExtEtaMin = hSparse->GetAxis(1)->FindBin(fDeltaEtaMin + 0.0001);
  Int_t binExtEtaMax = hSparse->GetAxis(1)->FindBin(fDeltaEtaMax - 0.0001);
  if (binExtEtaMax > hSparse->GetAxis(1)->GetNbins()) {
    binExtEtaMax = hSparse->GetAxis(1)->GetNbins();
  }
  if (binExtEtaMin < 1) {
    binExtEtaMin = 1;
  }

  hSparse->GetAxis(1)->SetRange(binExtEtaMin, binExtEtaMax);       // axis1: deltaEta
  hSparse->GetAxis(2)->SetRange(binExtPtCandMin, binExtPtCandMax); // axis2: ptCand
  hSparse->GetAxis(3)->SetRange(binExtPtHadMin, binExtPtHadMax);   // axis3: ptHad
  if (partType == kAllPart) {
    hSparse->GetAxis(4)->SetRange(binExtPoolMin, binExtPoolMax); // axis4: pool bin
  }

  h1D = hSparse->Projection(0);   // axis0: deltaPhi
  if (partType == kPrimaryPart) { // primary particles
    h1D->SetName(Form("hPrimaryPartCorr_PtD%.0fto%.0f_PtHad%.0fto%.0f", ptCandMin, ptCandMax, ptHadMin, ptHadMax));
  } else { // all selected particles
    h1D->SetName(Form("hAllPartCorr_PtD%.0fto%.0f_PtHad%.0fto%.0f", ptCandMin, ptCandMax, ptHadMin, ptHadMax));
  }

  return h1D;
}

TH1D* DhCorrelationExtraction::reflectCorrHistogram(TH1D*& histo)
{

  // nBinsPhi must be a multple of 4 in order to reflect correcty the histogram
  Int_t const nBinsPhi = histo->GetNbinsX();
  Int_t const nBinsPhiRefl = nBinsPhi / 2;
  Int_t const bin0Phi = nBinsPhi / 4 + 1;
  Int_t const binPiPhi = 3 * nBinsPhi / 4;

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

TH1D* DhCorrelationExtraction::reflectHistoRun2(TH1D* h, Double_t scale)
{

  TH1D* h2 = new TH1D(Form("%sReflected", h->GetName()), Form("%sReflected", h->GetName()), h->GetNbinsX() / 2., 0., TMath::Pi());
  for (Int_t j = 1; j <= h->GetNbinsX(); j++) {
    Double_t const x = h->GetBinCenter(j);
    Double_t const y0 = h->GetBinContent(j);
    Double_t const ey0 = h->GetBinError(j);
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
    Double_t const y = h2->GetBinContent(j2);
    Double_t const ey = h2->GetBinError(j2);
    h2->SetBinContent(j2, (y + y0));
    h2->SetBinError(j2, TMath::Sqrt(ey0 * ey0 + ey * ey));
  }
  h2->Scale(scale);

  return h2;
}

Double_t DhCorrelationExtraction::getFdPromptFrac(Double_t ptCandMin, Double_t ptCandMax, Double_t /*ptHadMin*/, Double_t /*ptHadMax*/)
{

  TH1D* h1D = reinterpret_cast<TH1D*>(fFileFDPromptFrac->Get(fHistoFDPromptFracName.Data()));

  Int_t const binPtCandMin = h1D->GetXaxis()->FindBin(ptCandMin + 0.01);
  Int_t const binPtCandMax = h1D->GetXaxis()->FindBin(ptCandMax - 0.01);
  Double_t promptFraction;
  if (binPtCandMin == binPtCandMax) {
    promptFraction = h1D->GetBinContent(binPtCandMin);
  } else {
    std::cout << "[ERROR] Different bin obtained from PtCandMin and PtCandMax";
    return 0.;
  }

  return promptFraction;
}

void DhCorrelationExtraction::normalizeMePlot(TH2D*& histoME, TH2D*& histoMEsoftPi) const
{

  Int_t const bin0phi = histoME->GetYaxis()->FindBin(0.);
  Int_t const bin0eta = histoME->GetXaxis()->FindBin(0.);

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

  if (fSubtractSoftPiME) {
    histoME->Add(histoMEsoftPi, -1); // remove the tracks compatible with soft pion (if requested)
  }

  // apply the normalization
  histoME->Scale(1. / factorNorm);
}

Double_t DhCorrelationExtraction::calculateBaseline(TH1D*& histo, Bool_t totalRange, Bool_t reflected)
{

  // total range = 2*Pi
  // half range = Pi , for histogram reflected under symmetric assumption

  Double_t baseline;
  Int_t const nBinsPhi = histo->GetNbinsX();
  Int_t const binPhiHalf = nBinsPhi / 2;
  Int_t const binPhiHalfMinus1 = nBinsPhi / 2 - 1;
  Int_t const binPhiHalfPlus1 = nBinsPhi / 2 + 1;
  Int_t const binPhiHalfPlus2 = nBinsPhi / 2 + 1;

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

Double_t DhCorrelationExtraction::calculateBaselineError(TH1D*& histo, Bool_t totalRange, Bool_t reflected)
{

  // total range = 2*Pi
  // half range = Pi , for histogram reflected under symmetric assumption

  Double_t errBaseline;
  Int_t const nBinsPhi = histo->GetNbinsX();
  Int_t const binPhiHalf = nBinsPhi / 2;
  Int_t const binPhiHalfMinus1 = nBinsPhi / 2 - 1;
  Int_t const binPhiHalfPlus1 = nBinsPhi / 2 + 1;
  Int_t const binPhiHalfPlus2 = nBinsPhi / 2 + 1;

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

void DhCorrelationExtraction::setTH1HistoStyle(TH1D*& histo, TString hTitle, TString hXaxisTitle, TString hYaxisTitle,
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

void DhCorrelationExtraction::setTH2HistoStyle(TH2D*& histo, TString hTitle, TString hXaxisTitle, TString hYaxisTitle, TString hZaxisTitle,
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
}
