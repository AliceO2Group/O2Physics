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

/// \file DhCorrelationFitter.cxx
/// \brief Class to perform the azimuthal correlation fit
/// \author Samuele Cattaruzzi <samuele.cattaruzzi@cern.ch>
/// \author Swapnesh Santosh Khade <swapnesh.santosh.khade@cern.ch>

#include "DhCorrelationFitter.h"

#include <cstdio>
#include <iostream>
#include <sstream>

#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TF2.h>
#include <TMatrixD.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <Riostream.h>
#include <TBufferFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TPaveLabel.h>
#include <TVirtualPad.h>
#include <TMath.h>
#include <TLatex.h>
#include <TColor.h>
#include <TClass.h>
#include <TVirtualFitter.h>
#include <TMinuit.h>

using namespace std;

DhCorrelationFitter::DhCorrelationFitter() : // default constructor
                                             fIsReflected(kFALSE),
                                             fTypeOfFitFunc(kConstwoGaus),
                                             fFixBase(0),
                                             fFixMean(0),
                                             fMinCandPt(0.),
                                             fMaxCandPt(99.),
                                             fMinAssoPt(0.),
                                             fMaxAssoPt(99.),
                                             fNpars(0),
                                             fExtParsVals(0x0),
                                             fExtParsLowBounds(0x0),
                                             fExtParsUppBounds(0x0),
                                             fUseExternalPars(kFALSE),
                                             fNbasleinePoints(0),
                                             fBinsBaseline(0x0),
                                             fHist(0x0),
                                             fMinCorr(0),
                                             fMaxCorr(0),
                                             fBaseline(0.),
                                             fErrBaseline(0.),
                                             fFit(0x0),
                                             fGausNS(0x0),
                                             fGausAS(0x0),
                                             fPed(0x0),
                                             fNSyieldBinCount(0.),
                                             fErrNSyieldBinCount(0.),
                                             fASyieldBinCount(0.),
                                             fErrASyieldBinCount(0.)
{
}

DhCorrelationFitter::DhCorrelationFitter(TH1F* histoToFit, Double_t min, Double_t max) : // standard constructor
                                                                                         fIsReflected(kFALSE),
                                                                                         fTypeOfFitFunc(kConstwoGaus),
                                                                                         fFixBase(0),
                                                                                         fFixMean(0),
                                                                                         fMinCandPt(0.),
                                                                                         fMaxCandPt(99.),
                                                                                         fMinAssoPt(0.),
                                                                                         fMaxAssoPt(99.),
                                                                                         fNpars(0),
                                                                                         fExtParsVals(0x0),
                                                                                         fExtParsLowBounds(0x0),
                                                                                         fExtParsUppBounds(0x0),
                                                                                         fUseExternalPars(kFALSE),
                                                                                         fNbasleinePoints(0),
                                                                                         fBinsBaseline(0x0),
                                                                                         fHist(0x0),
                                                                                         fMinCorr(0.),
                                                                                         fMaxCorr(0.),
                                                                                         fBaseline(0.),
                                                                                         fErrBaseline(0.),
                                                                                         fFit(0x0),
                                                                                         fGausNS(0x0),
                                                                                         fGausAS(0x0),
                                                                                         fPed(0x0),
                                                                                         fNSyieldBinCount(0.),
                                                                                         fErrNSyieldBinCount(0.),
                                                                                         fASyieldBinCount(0.),
                                                                                         fErrASyieldBinCount(0.)
{
  fHist = histoToFit;
  fMinCorr = min;
  fMaxCorr = max;
}

DhCorrelationFitter::DhCorrelationFitter(const DhCorrelationFitter& source) : // copy constructor
                                                                              fIsReflected(source.fIsReflected),
                                                                              fTypeOfFitFunc(source.fTypeOfFitFunc),
                                                                              fFixBase(source.fFixBase),
                                                                              fFixMean(source.fFixMean),
                                                                              fMinCandPt(source.fMinCandPt),
                                                                              fMaxCandPt(source.fMaxCandPt),
                                                                              fMinAssoPt(source.fMinAssoPt),
                                                                              fMaxAssoPt(source.fMaxAssoPt),
                                                                              fNpars(source.fNpars),
                                                                              fExtParsVals(source.fExtParsVals),
                                                                              fExtParsLowBounds(source.fExtParsLowBounds),
                                                                              fExtParsUppBounds(source.fExtParsUppBounds),
                                                                              fUseExternalPars(source.fUseExternalPars),
                                                                              fNbasleinePoints(source.fNbasleinePoints),
                                                                              fBinsBaseline(source.fBinsBaseline),
                                                                              fHist(source.fHist),
                                                                              fMinCorr(source.fMinCorr),
                                                                              fMaxCorr(source.fMaxCorr),
                                                                              fBaseline(source.fBaseline),
                                                                              fErrBaseline(source.fErrBaseline),
                                                                              fFit(source.fFit),
                                                                              fGausNS(source.fGausNS),
                                                                              fGausAS(source.fGausAS),
                                                                              fPed(source.fPed),
                                                                              fNSyieldBinCount(source.fNSyieldBinCount),
                                                                              fErrNSyieldBinCount(source.fErrNSyieldBinCount),
                                                                              fASyieldBinCount(source.fASyieldBinCount),
                                                                              fErrASyieldBinCount(source.fErrASyieldBinCount)
{
}

DhCorrelationFitter::~DhCorrelationFitter()
// destructor
{
  Info("DhCorrelationFitter.cxx", "Destructor is calling");

  if (fHist) {
    delete fHist;
    fHist = 0;
  }
  if (fFit) {
    delete fFit;
    fFit = 0;
  }
  if (fGausNS) {
    delete fGausNS;
    fGausNS = 0;
  }
  // if (fGausNS2) {delete fGausNS2; fGausNS2 = 0;}
  if (fPed) {
    delete fPed;
    fPed = 0;
  }
}

DhCorrelationFitter& DhCorrelationFitter::operator=(const DhCorrelationFitter& cfit)
// assignment operator
{
  if (&cfit == this)
    return *this;

  fIsReflected = cfit.fIsReflected;
  fTypeOfFitFunc = cfit.fTypeOfFitFunc;
  fFixBase = cfit.fFixBase;
  fFixMean = cfit.fFixMean;
  fMinCandPt = cfit.fMinCandPt;
  fMaxCandPt = cfit.fMaxCandPt;
  fMinAssoPt = cfit.fMinAssoPt;
  fMaxAssoPt = cfit.fMaxAssoPt;
  fNpars = cfit.fNpars;
  fExtParsVals = cfit.fExtParsVals;
  fExtParsLowBounds = cfit.fExtParsLowBounds;
  fExtParsUppBounds = cfit.fExtParsUppBounds;
  fUseExternalPars = cfit.fUseExternalPars;
  fNbasleinePoints = cfit.fNbasleinePoints;
  fBinsBaseline = cfit.fBinsBaseline;
  fHist = cfit.fHist;
  fMinCorr = cfit.fMinCorr;
  fMaxCorr = cfit.fMaxCorr;
  fBaseline = cfit.fBaseline;
  fErrBaseline = cfit.fErrBaseline;
  fFit = cfit.fFit;
  fGausNS = cfit.fGausNS;
  fGausAS = cfit.fGausAS;
  fPed = cfit.fPed;
  fNSyieldBinCount = cfit.fNSyieldBinCount;
  fErrNSyieldBinCount = cfit.fErrNSyieldBinCount;
  fASyieldBinCount = cfit.fASyieldBinCount;
  fErrASyieldBinCount = cfit.fErrASyieldBinCount;

  return *this;
}

void DhCorrelationFitter::SetExternalValsAndBounds(Int_t nPars, Double_t* vals, Double_t* lowBounds, Double_t* uppBounds)
{

  fNpars = nPars;

  fExtParsVals = new Double_t[fNpars];
  fExtParsLowBounds = new Double_t[fNpars];
  fExtParsUppBounds = new Double_t[fNpars];

  for (int i = 0; i < fNpars; i++) {
    fExtParsVals[i] = vals[i];
    fExtParsLowBounds[i] = lowBounds[i];
    fExtParsUppBounds[i] = uppBounds[i];
  }

  return;
}

void DhCorrelationFitter::Fitting(Bool_t drawSplitTerm, Bool_t useExternalPars)
{
  // -> fFixBase = 0 : baseline free
  //             = 1 : fix the baseline to the minimum of the histogram
  //             < 0 : fix the baseline to the weighted average of the abs(fFixBaseline) lower points
  //             = 2 : zyam at pi/2. Fix the baseline averaging the 2 points around +-pi/2 value
  //             = 3 : fix the baseline to the weighted average of the points passed through the function SetPointsForBaseline()
  //
  // -> fFixMean = 0 : NS & AS mean free
  //             = 1 : NS mean fixed to 0, AS mean free
  //             = 2 : AS mean fixed to pi, NS mean free
  //             = 3 : NS mean fixed to 0, AS mean to pi

  if (useExternalPars)
    fUseExternalPars = kTRUE;

  if (fFixBase != 0 && fFixBase != 6) {
    Printf("[INFO] DhCorrelationFitter::Fitting, Finding baseline");
    FindBaseline();
  }
  Printf("[INFO] DhCorrelationFitter::Fitting, Setting Function");
  SetFitFunction();

  if (fFixBase != 0) {
    fFit->FixParameter(0, fBaseline);
  }
  if (fFixMean == 1 || fFixMean == 3) {
    fFit->FixParameter(2, 0.);
  }
  if (fFixMean == 2 || fFixMean == 3) {
    if (fTypeOfFitFunc != 0)
      fFit->FixParameter(5, TMath::Pi());
  }
  Printf("[INFO] DhCorrelationFitter::Fitting, Fitting");
  TVirtualFitter::SetMaxIterations(20000);
  TFitResultPtr fitptr = fHist->Fit(fFit, "RIMES", "", fMinCorr, fMaxCorr);
  TMatrixD cor = fitptr->GetCorrelationMatrix();
  TMatrixD cov = fitptr->GetCovarianceMatrix();
  printf("[INFO] Correlation Matrix - The final one! \n");
  cor.Print();
  gMinuit->mnmatu(1);
  printf("[INFO] Covariance Matrix - The final one! \n");
  cov.Print();
  if (fFixBase == 0) {
    fBaseline = fFit->GetParameter(0);
    fErrBaseline = fFit->GetParError(0);
  }
  Printf("[INFO] DhCorrelationFitter::Fitting, Calculating yields with BC");
  CalculateYieldsAboveBaseline();
  fHist->SetTitle(";#Delta#varphi (rad); #frac{1}{N_{D}}#frac{dN^{assoc}}{d#Delta#varphi} (rad^{-1})");
  Printf("[INFO] DhCorrelationFitter::Fitting, Now drawing, if requested");
  SetSingleTermsForDrawing(drawSplitTerm);
}

void DhCorrelationFitter::SetFitFunction()
{
  // -> fitFunc = 1: const+ G NS + G AS (w/o periodicity)
  //            = 2: const+ G NS + G AS  (w/ periodicity)
  //            = 3: const+ yieldNS*[fact*(G NS)+(1- fact)*(G2 NS)] + yieldAS*(G AS)  (w/ periodicity)
  //            = 4: const +yieldNS*(G NS) + yieldAS*[fact*(G AS)+(1- fact)*(G2 AS)]   (w/ periodicity)
  //            = 5: v2 modulation (no gaussian terms)
  //            = 6: v2 modulation + G NS + G AS  (w/ periodicity)
  //            = 7: const+ GenG NS + G AS  (w/ periodicity)
  //            = 8: const+ GenG fixBeta NS + G AS  (w/ periodicity)
  //            = 9: const+ GenG constrBeta NS + G AS  (w/ periodicity)

  if (fFit) {
    delete fFit;
    delete fGausNS;
    // delete fGausNS2;
    delete fGausAS;
    // delete fGausAS2;
    delete fPed;
  }
  switch (fTypeOfFitFunc) {
    case 1:
      fFit = new TF1("ConstwoGaus", "[0]+[1]/TMath::Sqrt(2.*TMath::Pi())/[3]*TMath::Exp(-(x-[2])*(x-[2])/2./([3]*[3]))+[4]/TMath::Sqrt(2.*TMath::Pi())/[6]*TMath::Exp(-(x-[5])*(x-[5])/2./([6]*[6]))", fMinCorr, fMaxCorr);
      fGausNS = new TF1("fGausNS", "[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-[1])*(x-[1])/2./([2]*[2]))", fMinCorr, fMaxCorr);
      fGausAS = new TF1("fGausAS", "[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-[1])*(x-[1])/2./([2]*[2]))", fMinCorr, fMaxCorr);
      fPed = new TF1("fPed", "[0]", fMinCorr, fMaxCorr);

      if (!fUseExternalPars) {
        fFit->SetParLimits(0, 0, 9999.);
        fFit->SetParLimits(1, 0, 999.);
        fFit->SetParLimits(2, -1, 1);
        fFit->SetParLimits(3, 0, 3.14 / 3.);
        fFit->SetParLimits(4, 0, 999.);
        fFit->SetParLimits(5, 2., 4);
        fFit->SetParLimits(6, 0, 3.14 / 2.);

        fFit->SetParameter(0, 3);
        fFit->SetParameter(1, 2);
        fFit->SetParameter(2, 0.);
        fFit->SetParameter(3, 0.3);
        fFit->SetParameter(4, 2);
        fFit->SetParameter(5, 3.14);
        fFit->SetParameter(6, 0.3);

      } else {
        for (int i = 0; i < fNpars; i++) {
          fFit->SetParameter(i, fExtParsVals[i]);
          fFit->SetParLimits(i, fExtParsLowBounds[i], fExtParsUppBounds[i]);
        }
      }

      fFit->SetParName(0, "ped");
      fFit->SetParName(1, "NS Y");
      fFit->SetParName(2, "NS mean");
      fFit->SetParName(3, "NS #sigma");
      fFit->SetParName(4, "AS Y");
      fFit->SetParName(5, "AS mean");
      fFit->SetParName(6, "AS #sigma");

      break;

    case 2:
      fFit = new TF1("TwoGausPeriodicity", "[0]+[1]/TMath::Sqrt(2.*TMath::Pi())/[3]*TMath::Exp(-(x-[2])*(x-[2])/2./([3]*[3]))+[4]/TMath::Sqrt(2.*TMath::Pi())/[6]*TMath::Exp(-(x-[5])*(x-[5])/2./([6]*[6]))+[1]/TMath::Sqrt(2.*TMath::Pi())/[3]*TMath::Exp(-(x-2.*TMath::Pi()-[2])*(x-2.*TMath::Pi()-[2])/2./([3]*[3]))+[1]/TMath::Sqrt(2.*TMath::Pi())/[3]*TMath::Exp(-(x+2.*TMath::Pi()-[2])*(x+2.*TMath::Pi()-[2])/2./([3]*[3]))+[4]/TMath::Sqrt(2.*TMath::Pi())/[6]*TMath::Exp(-(x+2.*TMath::Pi()-[5])*(x+2.*TMath::Pi()-[5])/2./([6]*[6]))+[4]/TMath::Sqrt(2.*TMath::Pi())/[6]*TMath::Exp(-(x-2.*TMath::Pi()-[5])*(x-2.*TMath::Pi()-[5])/2./([6]*[6]))", fMinCorr, fMaxCorr);
      fGausNS = new TF1("fGausNSper", "[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-[1])*(x-[1])/2./([2]*[2]))+[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-2.*TMath::Pi()-[1])*(x-2.*TMath::Pi()-[1])/2./([2]*[2]))+[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x+2.*TMath::Pi()-[1])*(x+2.*TMath::Pi()-[1])/2./([2]*[2]))", fMinCorr, fMaxCorr);
      fGausAS = new TF1("fGausASper", "[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-[1])*(x-[1])/2./([2]*[2]))+[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-2.*TMath::Pi()-[1])*(x-2.*TMath::Pi()-[1])/2./([2]*[2]))+[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x+2.*TMath::Pi()-[1])*(x+2.*TMath::Pi()-[1])/2./([2]*[2]))", fMinCorr, fMaxCorr);
      fPed = new TF1("fPed", "[0]", fMinCorr, fMaxCorr);

      if (!fUseExternalPars) {
        fFit->SetParLimits(0, 0., 999.);
        fFit->SetParLimits(1, 0, 999.);
        fFit->SetParLimits(2, -0.55, 0.55);
        fFit->SetParLimits(3, 0, 3.14 / 3.);
        fFit->SetParLimits(4, 0, 999.);
        fFit->SetParLimits(5, 2.85, 3.55);
        fFit->SetParLimits(6, 0, 3.14 / 2.);

        fFit->SetParameter(0, 0.6);
        fFit->SetParameter(1, 3);
        fFit->SetParameter(2, 0.);
        fFit->SetParameter(3, 0.3);
        fFit->SetParameter(4, 2);
        fFit->SetParameter(5, TMath::Pi());
        fFit->SetParameter(6, 0.3);

      } else {
        for (int i = 0; i < fNpars; i++) {
          fFit->SetParameter(i, fExtParsVals[i]);
          fFit->SetParLimits(i, fExtParsLowBounds[i], fExtParsUppBounds[i]);
        }
      }

      fFit->SetParName(0, "ped");
      fFit->SetParName(1, "NS Y");
      fFit->SetParName(2, "NS mean");
      fFit->SetParName(3, "NS #sigma");
      fFit->SetParName(4, "AS Y");
      fFit->SetParName(5, "AS mean");
      fFit->SetParName(6, "AS #sigma");

      break;
  }
}

void DhCorrelationFitter::SetPointsForBaseline(Int_t nBaselinePoints, Int_t* binsBaseline)
{

  fNbasleinePoints = nBaselinePoints;

  fBinsBaseline = new Int_t[fNbasleinePoints];

  for (int i = 0; i < fNbasleinePoints; i++) {
    fBinsBaseline[i] = binsBaseline[i];
  }

  return;
}

Double_t DhCorrelationFitter::FindBaseline()
{

  // baseline free
  if (fFixBase == 0) {
    Printf("[INFO] DhCorrelationFitter::FindBasline(). The baseline option is set to free baseline: now the full fit will be done. Beware!");
    Fitting(); // TODO: not sure
    return fBaseline;
  }

  // fix the baseline to the minimum of the histogram
  if (fFixBase == 1) {
    Double_t min = 1.e10;
    Int_t iBin = -1;
    for (Int_t j = 1; j <= fHist->GetNbinsX(); j++) {
      if (fHist->GetBinContent(j) < min) {
        min = fHist->GetBinContent(j);
        iBin = j;
      }
    }
    fBaseline = min;
    fErrBaseline = fHist->GetBinError(iBin);

    return fBaseline;
  }

  // fix the baseline to the weighted average of the abs(fFixBaseline) lower points
  if (fFixBase < 0) {
    Int_t npointsAv = TMath::Abs(fFixBase);
    Int_t* ind = new Int_t[fHist->GetNbinsX()];
    Float_t* hval = new Float_t[fHist->GetNbinsX()];
    for (Int_t k = 1; k <= fHist->GetNbinsX(); k++) {
      hval[k - 1] = fHist->GetBinContent(k);
    }
    Double_t errAv = 0., Av = 0.;
    TMath::Sort(fHist->GetNbinsX(), hval, ind, kFALSE); //  KFALSE -> increasing order
    // Average of abs(fFixBase) lower points
    for (Int_t k = 0; k < npointsAv; k++) {
      if (fHist->GetBinError(ind[k] + 1) == 0.) // in case of null entries which induce a crash. Could bias the basline in upward direction!
      {
        printf("[WARNING] Null entries found in histogram to be fit! These points are been excluded from baseline evaluation...\n");
        npointsAv++;
        continue;
      }
      Av += fHist->GetBinContent(ind[k] + 1) / (fHist->GetBinError(ind[k] + 1) * fHist->GetBinError(ind[k] + 1));
      errAv += 1. / (fHist->GetBinError(ind[k] + 1) * fHist->GetBinError(ind[k] + 1));
    }
    Av /= errAv;
    errAv = TMath::Sqrt(1. / errAv);
    printf("[RESULT] Average fBaseline: %.3f +- %.3f", Av, errAv);
    fBaseline = Av;
    fErrBaseline = errAv;
    return fBaseline;
  }

  // zyam at pi/2. Fix the baseline averaging the 2 points around +-pi/2 value
  if (fFixBase == 2) {
    Double_t errAv = 0., Av = 0.;
    Int_t binPhi = fHist->FindBin(TMath::Pi() / 2.);
    Av += fHist->GetBinContent(binPhi) / (fHist->GetBinError(binPhi) * fHist->GetBinError(binPhi));
    errAv += 1. / (fHist->GetBinError(binPhi) * fHist->GetBinError(binPhi));
    if (!fIsReflected) {
      binPhi = fHist->FindBin(-TMath::Pi() / 2.);
      if (binPhi < 1)
        binPhi = 1;
      Av += fHist->GetBinContent(binPhi) / (fHist->GetBinError(binPhi) * fHist->GetBinError(binPhi));
      errAv += 1. / (fHist->GetBinError(binPhi) * fHist->GetBinError(binPhi));
    } else {
      printf("[INFO] Reflected histo: only the point at +pi/2 used to evaluate baseline");
    }
    Av /= errAv;
    errAv = TMath::Sqrt(1. / errAv);
    printf("[RESULT] Average fBaseline: %.3f +- %.3f \n", Av, errAv);
    fBaseline = Av;
    fErrBaseline = errAv;

    return fBaseline;
  }

  // fix the baseline to the weighted average of the points passed through the function SetPointsForBaseline()
  if (fFixBase == 3) {
    if (fNbasleinePoints == 0) {
      printf("[ERROR] No baseline points set for the baseline evaluation, SetPointsForBaseline(Int_t nBaselinePoints, Double_t* valsBaseline). Returning -1");
      return -1;
    }
    Double_t errAv = 0., Av = 0.;
    for (int i = 0; i < fNbasleinePoints; i++) {
      Av += fHist->GetBinContent(fBinsBaseline[i]) / (fHist->GetBinError(fBinsBaseline[i]) * fHist->GetBinError(fBinsBaseline[i]));
      errAv += 1. / (fHist->GetBinError(fBinsBaseline[i]) * fHist->GetBinError(fBinsBaseline[i]));
    }
    Av /= errAv;
    errAv = TMath::Sqrt(1. / errAv);
    printf("[RESULT] Average fBaseline: %.3f +- %.3f \n", Av, errAv);
    fBaseline = Av;
    fErrBaseline = errAv;

    return fBaseline;
  }

  Printf("[ERROR] DhCorrelationFitter::FindBaseline - WRONG BASELINE OPTION SET. Returning -1");
  return -1.;
}

void DhCorrelationFitter::CalculateYieldsAboveBaseline()
{

  fNSyieldBinCount = 0;
  fErrNSyieldBinCount = 0;
  fASyieldBinCount = 0;
  fErrASyieldBinCount = 0;
  cout << "[RESULT] Baseline: " << fBaseline << " +- " << fErrBaseline << endl;
  Int_t binMinNS = fHist->FindBin(-1.5); // slightly more than -pi/2
  if (binMinNS < 1)
    binMinNS = 1;                              // with this, it is ok even in the case of a reflected fHist (range 0 - pi)
  Int_t binMaxNS = fHist->FindBin(1.5);        // slightly less than +pi/2
  Int_t binMinAS = fHist->FindBin(1.6);        // slightly more than +pi/2
  Int_t binMaxAS = fHist->FindBin(3.14 + 1.5); // slightly less than +3pi/2
  if (binMaxAS > fHist->GetNbinsX())
    binMaxNS = fHist->GetNbinsX(); // with this, it is ok even in the case of a reflected fHist (range 0 - pi) TODO: maybe binMaxAS
  // Near Side Yield from bin counting
  for (Int_t bmNS = binMinNS; bmNS <= binMaxNS; bmNS++) {
    fNSyieldBinCount += (fHist->GetBinContent(bmNS) - fBaseline) * fHist->GetBinWidth(bmNS);
    fErrNSyieldBinCount += (fHist->GetBinError(bmNS) * fHist->GetBinError(bmNS)) * fHist->GetBinWidth(bmNS) * fHist->GetBinWidth(bmNS);
  }
  fErrNSyieldBinCount = TMath::Sqrt(fErrNSyieldBinCount);

  // Away Side Yield from bin counting
  for (Int_t bmAS = binMinAS; bmAS <= binMaxAS; bmAS++) {
    fASyieldBinCount += (fHist->GetBinContent(bmAS) - fBaseline) * fHist->GetBinWidth(bmAS);
    fErrASyieldBinCount += (fHist->GetBinError(bmAS) * fHist->GetBinError(bmAS)) * fHist->GetBinWidth(bmAS) * fHist->GetBinWidth(bmAS);
  }
  fErrASyieldBinCount = TMath::Sqrt(fErrASyieldBinCount);

  printf("[RESULT] Bin counting results: NS Yield = %.3f +- %.3f \n[RESULT] Bin counting results: AS Yield: %.3f +- %.3f \n", fNSyieldBinCount, fErrNSyieldBinCount, fASyieldBinCount, fErrASyieldBinCount);

  return;
}

void DhCorrelationFitter::SetSingleTermsForDrawing(Bool_t draw)
{
  Double_t* par = 0;
  if (fTypeOfFitFunc == 1 || fTypeOfFitFunc == 2) {
    par = new Double_t[7];
  } else {
    Printf("[ERROR] DhCorrelationFitter::SetSingleTermsForDrawing, wrong type of function");
    return;
  }

  fFit->GetParameters(par);
  fFit->SetLineWidth(4);

  fPed->SetParameter(0, par[0]);
  fPed->SetLineColor(6); // pink
  fPed->SetLineStyle(9);
  fPed->SetLineWidth(4);

  fGausNS->SetParameter(0, par[1]);
  fGausNS->SetParameter(1, par[2]);
  fGausNS->SetParameter(2, par[3]);
  fGausAS->SetParameter(0, par[4]);
  fGausAS->SetParameter(1, par[5]);
  fGausAS->SetParameter(2, par[6]);

  fGausNS->SetLineStyle(9);
  fGausNS->SetLineColor(kBlue);
  fGausAS->SetLineStyle(9);
  fGausAS->SetLineColor(kGreen);
  fGausNS->SetLineWidth(4);
  fGausAS->SetLineWidth(4);

  TPaveText* pvStatTests1 = new TPaveText(0.51, 0.58, 0.85, 0.90, "NDC");
  pvStatTests1->SetFillStyle(0);
  pvStatTests1->SetTextSize(0.035);
  pvStatTests1->SetBorderSize(0);
  TText *t0, *t1, *t2, *t3, *t4, *t5, *t6;
  t0 = pvStatTests1->AddText(0., 1.00, Form("#chi^{2}/ndf = %.1f/%d ", fFit->GetChisquare(), fFit->GetNDF()));
  t1 = pvStatTests1->AddText(0., 0.80, Form("Ped = %.3f#pm%.3f ", fBaseline, fErrBaseline));
  t2 = pvStatTests1->AddText(0., 0.65, Form("NS Y = %.3f#pm%.3f ", fFit->GetParameter("NS Y"), fFit->GetParError(fFit->GetParNumber("NS Y"))));
  t3 = pvStatTests1->AddText(0., 0.50, Form("NS #sigma = %.3f#pm%.3f ", fFit->GetParameter("NS #sigma"), fFit->GetParError(fFit->GetParNumber("NS #sigma"))));
  t4 = pvStatTests1->AddText(0., 0.35, Form("AS Y = %.3f#pm%.3f ", fFit->GetParameter("AS Y"), fFit->GetParError(fFit->GetParNumber("AS Y"))));
  t5 = pvStatTests1->AddText(0., 0.20, Form("AS #sigma = %.3f#pm%.3f ", fFit->GetParameter("AS #sigma"), fFit->GetParError(fFit->GetParNumber("AS #sigma"))));

  if (draw) {
    fFit->Draw("same");
    fPed->Draw("same");
    fGausAS->Draw("same");
    fGausNS->Draw("same");
    pvStatTests1->Draw("same");
  }
}
