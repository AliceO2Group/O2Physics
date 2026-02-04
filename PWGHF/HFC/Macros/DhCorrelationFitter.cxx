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
/// \brief class for for performing the fit of azimuthal correlations
/// \author Samuele Cattaruzzi <samuele.cattaruzzi@cern.ch>
/// \author Swapnesh Santosh Khade <swapnesh.santosh.khade@cern.ch>

#include "DhCorrelationFitter.h"

#include <TError.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMath.h>
#include <TMathBase.h>
#include <TMatrixD.h> // IWYU pragma: keep (do not replace with TMatrixDfwd.h)
#include <TMatrixDfwd.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TString.h>
#include <TVirtualFitter.h>
#include <TVirtualPad.h>

#include <Rtypes.h>
#include <RtypesCore.h>

#include <cstdio>
#include <iostream>

DhCorrelationFitter::DhCorrelationFitter() : // default constructor
                                             fHist(nullptr),
                                             fFit(nullptr),
                                             fGausNS(nullptr),
                                             fGausAS(nullptr),
                                             fPed(nullptr),
                                             fBaseTransvReg(nullptr),
                                             fIsReflected(kFALSE),
                                             fUseExternalPars(kFALSE),
                                             fShiftBaselineUp(kFALSE),
                                             fShiftBaselineDown(kFALSE),
                                             fIsTotal(kTRUE),
                                             fTypeOfFitFunc(kConstwoGaus),
                                             fFixBase(0),
                                             fFixMean(0),
                                             fNpars(0),
                                             fNbasleinePoints(0),
                                             fBinsBaseline(nullptr),
                                             fMinCorr(0),
                                             fMaxCorr(0),
                                             fMinCandPt(0.),
                                             fMaxCandPt(99.),
                                             fMinAssoPt(0.),
                                             fMaxAssoPt(99.),
                                             fBaseline(0.),
                                             fErrBaseline(0.),
                                             fNSyieldBinCount(0.),
                                             fErrNSyieldBinCount(0.),
                                             fASyieldBinCount(0.),
                                             fErrASyieldBinCount(0.),
                                             fv2AssocPart(0.),
                                             fv2Dmeson(0.),
                                             fExtParsVals(nullptr),
                                             fExtParsLowBounds(nullptr),
                                             fExtParsUppBounds(nullptr)
{
}

DhCorrelationFitter::DhCorrelationFitter(TH1F* histoToFit, Double_t min, Double_t max) : // standard constructor
                                                                                         fHist(histoToFit),
                                                                                         fFit(nullptr),
                                                                                         fGausNS(nullptr),
                                                                                         fGausAS(nullptr),
                                                                                         fPed(nullptr),
                                                                                         fBaseTransvReg(nullptr),
                                                                                         fIsReflected(kFALSE),
                                                                                         fUseExternalPars(kFALSE),
                                                                                         fShiftBaselineUp(kFALSE),
                                                                                         fShiftBaselineDown(kFALSE),
                                                                                         fIsTotal(kTRUE),
                                                                                         fTypeOfFitFunc(kConstwoGaus),
                                                                                         fFixBase(0),
                                                                                         fFixMean(0),
                                                                                         fNpars(0),
                                                                                         fNbasleinePoints(0),
                                                                                         fBinsBaseline(nullptr),
                                                                                         fMinCorr(min),
                                                                                         fMaxCorr(max),
                                                                                         fMinCandPt(0.),
                                                                                         fMaxCandPt(99.),
                                                                                         fMinAssoPt(0.),
                                                                                         fMaxAssoPt(99.),
                                                                                         fBaseline(0.),
                                                                                         fErrBaseline(0.),
                                                                                         fNSyieldBinCount(0.),
                                                                                         fErrNSyieldBinCount(0.),
                                                                                         fASyieldBinCount(0.),
                                                                                         fErrASyieldBinCount(0.),
                                                                                         fv2AssocPart(0.),
                                                                                         fv2Dmeson(0.),
                                                                                         fExtParsVals(nullptr),
                                                                                         fExtParsLowBounds(nullptr),
                                                                                         fExtParsUppBounds(nullptr)
{
}

DhCorrelationFitter::DhCorrelationFitter(const DhCorrelationFitter& source)

  = default;

DhCorrelationFitter::~DhCorrelationFitter()
// destructor
{
  Info("DhCorrelationFitter.cxx", "Destructor is calling");

  if (fHist != nullptr) {
    delete fHist;
    fHist = nullptr;
  }
  if (fFit != nullptr) {
    delete fFit;
    fFit = nullptr;
  }
  if (fGausNS != nullptr) {
    delete fGausNS;
    fGausNS = nullptr;
  }
  // if (fGausNS2) {delete fGausNS2; fGausNS2 = 0;}
  if (fPed != nullptr) {
    delete fPed;
    fPed = nullptr;
  }
}

DhCorrelationFitter& DhCorrelationFitter::operator=(const DhCorrelationFitter& cfit)
// assignment operator
{
  if (&cfit == this) {
    return *this;
  }

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
  fShiftBaselineUp = cfit.fShiftBaselineUp;
  fShiftBaselineDown = cfit.fShiftBaselineDown;
  fIsTotal = cfit.fIsTotal;
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
  fBaseTransvReg = cfit.fBaseTransvReg;
  fv2AssocPart = cfit.fv2AssocPart;
  fv2Dmeson = cfit.fv2Dmeson;
  fNSyieldBinCount = cfit.fNSyieldBinCount;
  fErrNSyieldBinCount = cfit.fErrNSyieldBinCount;
  fASyieldBinCount = cfit.fASyieldBinCount;
  fErrASyieldBinCount = cfit.fErrASyieldBinCount;

  return *this;
}

void DhCorrelationFitter::setExternalValsAndBounds(Int_t nPars, const Double_t* vals, const Double_t* lowBounds, const Double_t* uppBounds)
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
}

void DhCorrelationFitter::fitting(Bool_t drawSplitTerm, Bool_t useExternalPars)
{
  // -> fFixBase = 0 : baseline free
  //             = 1 : fix the baseline to the minimum of the histogram
  //             < 0 : fix the baseline to the weighted average of the abs(fFixBaseline) lower points
  //             = 2 : zyam at pi/2. Fix the baseline averaging the 2 points around +-pi/2 value
  //             = 3 : fix the baseline to the weighted average of the points passed through the function SetPointsForBaseline()
  //             = 4 : fix the baseline to the weighted average of the points in the transverse region in default configuration (i.e. for 32 bins the first 2, the last 2 and the 4 middle ones)
  //
  // -> fFixMean = 0 : NS & AS mean free
  //             = 1 : NS mean fixed to 0, AS mean free
  //             = 2 : AS mean fixed to pi, NS mean free
  //             = 3 : NS mean fixed to 0, AS mean to pi

  if (useExternalPars) {
    fUseExternalPars = kTRUE;
  }

  if (fFixBase != 0 && fFixBase != 6) {
    Printf("[INFO] DhCorrelationFitter::Fitting, Finding baseline");
    findBaseline();
  }
  if (fFixBase == 0) {
    // set initial value of the fBaseline
    fBaseline = calculateBaseline(fHist, fIsTotal);
  }
  Printf("[INFO] DhCorrelationFitter::Fitting, Setting Function");
  if (fTypeOfFitFunc == 7) { // case for v2 modulation
    fitBaselineWv2();        // to contrain the B parameter in the fit function for the pedestal
    Printf("[INFO] B parameter for v2 fit: %.3f", fBaseline);
  }
  setFitFunction();

  if (fFixBase != 0) {
    fFit->FixParameter(0, fBaseline);
  }
  if (fFixMean == 1 || fFixMean == 3) {
    fFit->FixParameter(2, 0.);
  }
  if (fFixMean == 2 || fFixMean == 3) {
    if (fTypeOfFitFunc != 0 && fTypeOfFitFunc != 3) {
      fFit->FixParameter(5, TMath::Pi());
    }
    if (fTypeOfFitFunc == 3 || fTypeOfFitFunc == 6) {
      fFit->FixParameter(2, TMath::Pi());
    }
  }

  Printf("[INFO] DhCorrelationFitter::Fitting, Fitting");
  TVirtualFitter::SetMaxIterations(20000);
  TFitResultPtr const fitptr = fHist->Fit(fFit, "RIMES", "", fMinCorr, fMaxCorr);
  TMatrixD const cor = fitptr->GetCorrelationMatrix();
  TMatrixD const cov = fitptr->GetCovarianceMatrix();
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
  calculateYieldsAboveBaseline();
  fHist->SetTitle(";#Delta#varphi (rad); #frac{1}{N_{D}}#frac{dN^{assoc}}{d#Delta#varphi} (rad^{-1})");
  Printf("[INFO] DhCorrelationFitter::Fitting, Now drawing, if requested");
  setSingleTermsForDrawing(drawSplitTerm);

  // NS yield from bin counting
  double fNSyield = 0.;
  double fNSyieldErr = 0.;
  double baselinBinCount = 0;
  for (int iBin = 1; iBin <= 6; iBin++) {                                                // first six bins
    fNSyield += 2 * (fHist->GetBinContent(iBin) - fBaseline) * fHist->GetBinWidth(iBin); // x2 due to the fatct the histogram is reflected
    fNSyieldErr += 4 * (TMath::Power(fHist->GetBinError(iBin), 2) + TMath::Power(fErrBaseline, 2)) * (fHist->GetBinWidth(iBin) * fHist->GetBinWidth(iBin));
    baselinBinCount += fBaseline * fHist->GetBinWidth(iBin);
  }
  fNSyieldErr = TMath::Sqrt(fNSyieldErr);

  // AS yield from bin counting
  double fASyield = 0.;
  double fASyieldErr = 0.;
  for (int iBin = 11; iBin <= 16; iBin++) { // last six bins
    fASyield += 2 * (fHist->GetBinContent(iBin) - fBaseline) * fHist->GetBinWidth(iBin);
    fASyieldErr += 4 * (TMath::Power(fHist->GetBinError(iBin), 2) + TMath::Power(fErrBaseline, 2)) * (fHist->GetBinWidth(iBin) * fHist->GetBinWidth(iBin));
  }
  fASyieldErr = TMath::Sqrt(fASyieldErr);

  printf("[RESULT MINE] Bin counting results: NS Yield = %.3f +- %.3f \n[RESULT MINE] Bin counting results: AS Yield: %.3f +- %.3f \n[RESULT MINE] baseline = %.3f \n", fNSyield, fNSyieldErr, fASyield, fASyieldErr, baselinBinCount);
}

void DhCorrelationFitter::setFitFunction()
{
  // -> fitFunc = 1: const + G NS + G AS (w/o periodicity)
  //            = 2: const + G NS + G AS  (w/ periodicity)
  //            = 3: const + G AS
  //            = 4: const + GenG NS + G AS  (w/ periodicity)
  //            = 5: const + VonMises NS + VonMises AS (w/periodicity)
  //            = 6: const + VonMises AS
  //            = 7: baseline w v2 modulation + G NS + G AS  (w/ periodicity)

  if (fFit != nullptr) {
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

    case 3:
      fFit = new TF1("OneGausPeriodicity", "[0]+[1]/TMath::Sqrt(2.*TMath::Pi())/[3]*TMath::Exp(-(x-[2])*(x-[2])/2./([3]*[3]))+[1]/TMath::Sqrt(2.*TMath::Pi())/[3]*TMath::Exp(-(x+2.*TMath::Pi()-[2])*(x+2.*TMath::Pi()-[2])/2./([3]*[3]))+[1]/TMath::Sqrt(2.*TMath::Pi())/[3]*TMath::Exp(-(x-2.*TMath::Pi()-[2])*(x-2.*TMath::Pi()-[2])/2./([3]*[3]))", fMinCorr, fMaxCorr);
      fGausAS = new TF1("fGausASper", "[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-[1])*(x-[1])/2./([2]*[2]))+[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-2.*TMath::Pi()-[1])*(x-2.*TMath::Pi()-[1])/2./([2]*[2]))+[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x+2.*TMath::Pi()-[1])*(x+2.*TMath::Pi()-[1])/2./([2]*[2]))", fMinCorr, fMaxCorr);
      fPed = new TF1("fPed", "[0]", fMinCorr, fMaxCorr);

      // TODO: add possibility to use external parameters
      // if(!fUseExternalPars) {
      fFit->SetParLimits(0, 0, 9999.);
      fFit->SetParLimits(1, 0, 999.);
      fFit->SetParLimits(2, 2., 4.);
      fFit->SetParLimits(3, 0, 3.14 / 2.);

      fFit->SetParameter(0, 3);
      fFit->SetParameter(1, 2);
      fFit->SetParameter(2, 3.14);
      fFit->SetParameter(3, 0.3);

      /*} else {
        for(int i=0; i<fNpars; i++) {
          fFit -> SetParameter(i, fExtParsVals[i]);
          fFit -> SetParLimits(i, fExtParsLowBounds[i], fExtParsUppBounds[i]);
        }
      }  */

      fFit->SetParName(0, "ped");
      fFit->SetParName(1, "AS Y");
      fFit->SetParName(2, "AS mean");
      fFit->SetParName(3, "AS #sigma");

      break;

    case 4:
      fFit = new TF1("kModifNSGausPeriodicity", "[0]+[1]*([7]*TMath::Sqrt(TMath::Gamma(3./[7]))/(2.*[3]*TMath::Power(TMath::Gamma(1./[7]),3./2.))*TMath::Exp(-TMath::Power(TMath::Abs(x-[2])*TMath::Sqrt(TMath::Gamma(3./[7]))/([3]*TMath::Sqrt(TMath::Gamma(1./[7]))),[7])))+[4]/TMath::Sqrt(2.*TMath::Pi())/[6]*TMath::Exp(-(x-[5])*(x-[5])/2./([6]*[6]))+[1]*([7]*TMath::Sqrt(TMath::Gamma(3./[7]))/(2.*[3]*TMath::Power(TMath::Gamma(1./[7]),3./2.))*TMath::Exp(-TMath::Power(TMath::Abs(x-2*TMath::Pi()-[2])*TMath::Sqrt(TMath::Gamma(3./[7]))/([3]*TMath::Sqrt(TMath::Gamma(1./[7]))),[7])))+[1]*([7]*TMath::Sqrt(TMath::Gamma(3./[7]))/(2.*[3]*TMath::Power(TMath::Gamma(1./[7]),3./2.))*TMath::Exp(-TMath::Power(TMath::Abs(x+2*TMath::Pi()-[2])*TMath::Sqrt(TMath::Gamma(3./[7]))/([3]*TMath::Sqrt(TMath::Gamma(1./[7]))),[7])))+[4]/TMath::Sqrt(2.*TMath::Pi())/[6]*TMath::Exp(-(x-2.*TMath::Pi()-[5])*(x-2.*TMath::Pi()-[5])/2./([6]*[6]))+[4]/TMath::Sqrt(2.*TMath::Pi())/[6]*TMath::Exp(-(x+2.*TMath::Pi()-[5])*(x+2.*TMath::Pi()-[5])/2./([6]*[6]))", fMinCorr, fMaxCorr);
      fGausNS = new TF1("fModGausNSper", "[0]*([3]*TMath::Sqrt(TMath::Gamma(3./[3]))/(2.*[2]*TMath::Power(TMath::Gamma(1./[3]),3./2.))*TMath::Exp(-TMath::Power(TMath::Abs(x-[1])*TMath::Sqrt(TMath::Gamma(3./[3]))/([2]*TMath::Sqrt(TMath::Gamma(1./[3]))),[3])))+[0]*([3]*TMath::Sqrt(TMath::Gamma(3./[3]))/(2.*[2]*TMath::Power(TMath::Gamma(1./[3]),3./2.))*TMath::Exp(-TMath::Power(TMath::Abs(x-2*TMath::Pi()-[1])*TMath::Sqrt(TMath::Gamma(3./[3]))/([2]*TMath::Sqrt(TMath::Gamma(1./[3]))),[3])))+[0]*([3]*TMath::Sqrt(TMath::Gamma(3./[3]))/(2.*[2]*TMath::Power(TMath::Gamma(1./[3]),3./2.))*TMath::Exp(-TMath::Power(TMath::Abs(x+2*TMath::Pi()-[1])*TMath::Sqrt(TMath::Gamma(3./[3]))/([2]*TMath::Sqrt(TMath::Gamma(1./[3]))),[3])))", fMinCorr, fMaxCorr);
      fGausAS = new TF1("fGausASper", "[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-[1])*(x-[1])/2./([2]*[2]))+[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-2.*TMath::Pi()-[1])*(x-2.*TMath::Pi()-[1])/2./([2]*[2]))+[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x+2.*TMath::Pi()-[1])*(x+2.*TMath::Pi()-[1])/2./([2]*[2]))", fMinCorr, fMaxCorr);
      fPed = new TF1("fPed", "[0]", fMinCorr, fMaxCorr);

      fFit->SetParLimits(0, 0., 999.);
      fFit->SetParLimits(1, 0.005, 25.);
      fFit->SetParLimits(2, -0.55, 0.55);
      fFit->SetParLimits(3, 0, 0.8);
      fFit->SetParLimits(4, 0.005, 25.);
      fFit->SetParLimits(5, 2.85, 3.55);
      fFit->SetParLimits(6, 0.05, 3.14 / 2.);
      fFit->SetParLimits(7, 0.5, 3.5);

      // default starting pars
      fFit->SetParameter(0, 1.);
      fFit->SetParameter(1, 1.);
      fFit->SetParameter(2, 0.);
      fFit->SetParameter(3, 0.3);
      fFit->SetParameter(4, 0.25);
      fFit->SetParameter(5, TMath::Pi());
      fFit->SetParameter(6, 0.3);
      fFit->SetParameter(7, 2);

      fFit->SetParName(0, "ped");
      fFit->SetParName(1, "NS Y");
      fFit->SetParName(2, "NS mean");
      fFit->SetParName(3, "NS #sigma");
      fFit->SetParName(4, "AS Y");
      fFit->SetParName(5, "AS mean");
      fFit->SetParName(6, "AS #sigma");
      fFit->SetParName(7, "NS shape par"); // beta of the gen. gaussian
      break;

    case 5:
      fFit = new TF1("kVonMises", "[0] +[1]/(2*TMath::Pi()*TMath::BesselI0([3]))*TMath::Exp([3]*TMath::Cos(x- 2*TMath::Pi() - [2])) + [4]/(2*TMath::Pi()*TMath::BesselI0([6]))*TMath::Exp([6]*TMath::Cos(x- 2*TMath::Pi()-[5]))", fMinCorr, fMaxCorr);
      fGausNS = new TF1("fVonMisesNS", "[0]/(2*TMath::Pi()*TMath::BesselI0([2]))*TMath::Exp([2]*TMath::Cos(x- 2*TMath::Pi() - [1]))", fMinCorr, fMaxCorr);
      fGausAS = new TF1("fVonMisesAS", "[0]/(2*TMath::Pi()*TMath::BesselI0([2]))*TMath::Exp([2]*TMath::Cos(x- 2*TMath::Pi()-[1]))", fMinCorr, fMaxCorr);
      fPed = new TF1("fPed", "[0]", fMinCorr, fMaxCorr);

      fFit->SetParLimits(0, 0., 999.);
      fFit->SetParLimits(1, 0.005, 25.);
      fFit->SetParLimits(2, -0.55, 0.55);
      fFit->SetParLimits(3, 0, 15.);
      fFit->SetParLimits(4, 0.005, 25.);
      fFit->SetParLimits(5, 2.85, 3.55);
      fFit->SetParLimits(6, 0., 15.);

      // default starting pars
      fFit->SetParameter(0, 3.5);
      fFit->SetParameter(1, 0.8);
      fFit->SetParameter(2, 0.);
      fFit->SetParameter(3, 1.);
      fFit->SetParameter(4, 1.5);
      fFit->SetParameter(5, TMath::Pi());
      fFit->SetParameter(6, 1.);

      fFit->SetParName(0, "ped");
      fFit->SetParName(1, "NS Y");
      fFit->SetParName(2, "NS mean");
      fFit->SetParName(3, "NS #sigma");
      fFit->SetParName(4, "AS Y");
      fFit->SetParName(5, "AS mean");
      fFit->SetParName(6, "AS #sigma");

      break;

    case 6:
      fFit = new TF1("kSingleVonMises", "[0] +[1]/(2*TMath::Pi()*TMath::BesselI0([3]))*TMath::Exp([3]*TMath::Cos(x- 2*TMath::Pi() - [2]))", fMinCorr, fMaxCorr);
      fGausAS = new TF1("fVonMisesAS", "[0]/(2*TMath::Pi()*TMath::BesselI0([2]))*TMath::Exp([2]*TMath::Cos(x- 2*TMath::Pi()-[1]))", fMinCorr, fMaxCorr);
      fPed = new TF1("fPed", "[0]", fMinCorr, fMaxCorr);

      fFit->SetParLimits(0, 0., 999.);
      fFit->SetParLimits(1, 0.005, 25.);
      fFit->SetParLimits(2, 2.85, 3.55);
      fFit->SetParLimits(3, 0., 15.);

      // default starting pars
      fFit->SetParameter(0, 3.5);
      fFit->SetParameter(1, 1.5);
      fFit->SetParameter(2, TMath::Pi());
      fFit->SetParameter(3, 1.);

      fFit->SetParName(0, "ped");
      fFit->SetParName(1, "AS Y");
      fFit->SetParName(2, "AS mean");
      fFit->SetParName(3, "AS #sigma");

      break;

    case 7: // case 2 Gaus w periodicity + v2 modulation

      fFit = new TF1("kTwoGausPeriodicityPlusV2modulation", "[1]/TMath::Sqrt(2.*TMath::Pi())/[3]*TMath::Exp(-(x-[2])*(x-[2])/2./([3]*[3]))+[4]/TMath::Sqrt(2.*TMath::Pi())/[6]*TMath::Exp(-(x-[5])*(x-[5])/2./([6]*[6]))+[1]/TMath::Sqrt(2.*TMath::Pi())/[3]*TMath::Exp(-(x-2.*TMath::Pi()-[2])*(x-2.*TMath::Pi()-[2])/2./([3]*[3]))+[1]/TMath::Sqrt(2.*TMath::Pi())/[3]*TMath::Exp(-(x+2.*TMath::Pi()-[2])*(x+2.*TMath::Pi()-[2])/2./([3]*[3]))+[4]/TMath::Sqrt(2.*TMath::Pi())/[6]*TMath::Exp(-(x+2.*TMath::Pi()-[5])*(x+2.*TMath::Pi()-[5])/2./([6]*[6]))+[4]/TMath::Sqrt(2.*TMath::Pi())/[6]*TMath::Exp(-(x-2.*TMath::Pi()-[5])*(x-2.*TMath::Pi()-[5])/2./([6]*[6]))+[0]*(1+2*[7]*[8]*TMath::Cos(2*x))", fMinCorr, fMaxCorr);
      fGausNS = new TF1("fGausNSper", "[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-[1])*(x-[1])/2./([2]*[2]))+[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-2.*TMath::Pi()-[1])*(x-2.*TMath::Pi()-[1])/2./([2]*[2]))+[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x+2.*TMath::Pi()-[1])*(x+2.*TMath::Pi()-[1])/2./([2]*[2]))", fMinCorr, fMaxCorr);
      fGausAS = new TF1("fGausASper", "[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-[1])*(x-[1])/2./([2]*[2]))+[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x-2.*TMath::Pi()-[1])*(x-2.*TMath::Pi()-[1])/2./([2]*[2]))+[0]/TMath::Sqrt(2.*TMath::Pi())/[2]*TMath::Exp(-(x+2.*TMath::Pi()-[1])*(x+2.*TMath::Pi()-[1])/2./([2]*[2]))", fMinCorr, fMaxCorr);
      fPed = new TF1("fPedv2Mod", "[0]*(1+2*[1]*[2]*TMath::Cos(2*x))", fMinCorr, fMaxCorr);

      fFit->SetParLimits(0, 0., 999.);
      fFit->SetParLimits(1, 0, 999.);
      fFit->SetParLimits(2, -0.55, 0.55);
      fFit->SetParLimits(3, 0, 3.14 / 3.);
      fFit->SetParLimits(4, 0, 999.);
      fFit->SetParLimits(5, 2.85, 3.55);
      fFit->SetParLimits(6, 0, 3.14 / 2.);
      fFit->SetParLimits(7, -1, 1);
      fFit->SetParLimits(8, -1, 1);

      fFit->FixParameter(0, fBaseline);
      fFit->SetParameter(1, 3);
      fFit->SetParameter(2, 0.);
      fFit->SetParameter(3, 0.3);
      fFit->SetParameter(4, 2);
      fFit->SetParameter(5, TMath::Pi());
      fFit->SetParameter(6, 0.3);
      fFit->SetParameter(7, 0);
      fFit->SetParameter(8, 0);

      if (fUseExternalPars) { // overwrites previous configuration
        for (int i = 0; i < fNpars; i++) {
          fFit->SetParameter(i, fExtParsVals[i]);
          fFit->SetParLimits(i, fExtParsLowBounds[i], fExtParsUppBounds[i]);
        }
      }

      fFit->FixParameter(7, fv2AssocPart);
      fFit->FixParameter(8, fv2Dmeson);

      fPed->FixParameter(0, fBaseline);
      fPed->FixParameter(1, fv2AssocPart);
      fPed->FixParameter(2, fv2Dmeson);

      fFit->SetParName(0, "ped");
      fFit->SetParName(1, "NS Y");
      fFit->SetParName(2, "NS mean");
      fFit->SetParName(3, "NS #sigma");
      fFit->SetParName(4, "AS Y");
      fFit->SetParName(5, "AS mean");
      fFit->SetParName(6, "AS #sigma");
      fFit->SetParName(7, "v_{2} hadron");
      fFit->SetParName(8, "v_{2} D meson");
      break;
  }
}

void DhCorrelationFitter::setPointsForBaseline(Int_t nBaselinePoints, const Int_t* binsBaseline)
{

  fNbasleinePoints = nBaselinePoints;

  fBinsBaseline = new Int_t[fNbasleinePoints];

  for (int i = 0; i < fNbasleinePoints; i++) {
    fBinsBaseline[i] = binsBaseline[i];
  }
}

Double_t DhCorrelationFitter::findBaseline()
{

  // baseline free
  if (fFixBase == 0) {
    Printf("[INFO] DhCorrelationFitter::FindBasline(). The baseline option is set to free baseline: now the full fit will be done. Beware!");
    fitting(); // TODO: not sure
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

    if (fShiftBaselineUp) {
      fBaseline += fErrBaseline;
      printf("[INFO] Shift baseline up of its statistical uncertainty");
    }

    if (fShiftBaselineDown) {
      fBaseline -= fErrBaseline;
      printf("[INFO] Shift baseline down of its statistical uncertainty");
    }

    return fBaseline;
  }

  // fix the baseline to the weighted average of the abs(fFixBaseline) lower points
  if (fFixBase < 0) {
    Int_t npointsAv = TMath::Abs(fFixBase);
    auto* ind = new Int_t[fHist->GetNbinsX()];
    auto* hval = new Float_t[fHist->GetNbinsX()];
    for (Int_t k = 1; k <= fHist->GetNbinsX(); k++) {
      hval[k - 1] = fHist->GetBinContent(k);
    }
    Double_t errAv = 0., av = 0.;
    TMath::Sort(fHist->GetNbinsX(), hval, ind, kFALSE); //  KFALSE -> increasing order
    delete[] hval;
    // Average of abs(fFixBase) lower points
    for (Int_t k = 0; k < npointsAv; k++) {
      if (fHist->GetBinError(ind[k] + 1) == 0.) // in case of null entries which induce a crash. Could bias the basline in upward direction!
      {
        printf("[WARNING] Null entries found in histogram to be fit! These points are been excluded from baseline evaluation...\n");
        npointsAv++;
        continue;
      }
      av += fHist->GetBinContent(ind[k] + 1) / (fHist->GetBinError(ind[k] + 1) * fHist->GetBinError(ind[k] + 1));
      errAv += 1. / (fHist->GetBinError(ind[k] + 1) * fHist->GetBinError(ind[k] + 1));
    }
    delete[] ind;
    av /= errAv;
    errAv = TMath::Sqrt(1. / errAv);
    printf("[RESULT] Average fBaseline: %.3f +- %.3f", av, errAv);
    fBaseline = av;
    fErrBaseline = errAv;

    if (fShiftBaselineUp) {
      fBaseline += fErrBaseline;
      printf("[INFO] Shift baseline up of its statistical uncertainty");
    }

    if (fShiftBaselineDown) {
      fBaseline -= fErrBaseline;
      printf("[INFO] Shift baseline down of its statistical uncertainty");
    }

    return fBaseline;
  }

  // zyam at pi/2. Fix the baseline averaging the 2 points around +-pi/2 value
  if (fFixBase == 2) {
    Double_t errAv = 0., av = 0.;
    Int_t binPhi = fHist->FindBin(TMath::Pi() / 2.);
    av += fHist->GetBinContent(binPhi) / (fHist->GetBinError(binPhi) * fHist->GetBinError(binPhi));
    errAv += 1. / (fHist->GetBinError(binPhi) * fHist->GetBinError(binPhi));
    if (!fIsReflected) {
      binPhi = fHist->FindBin(-TMath::Pi() / 2.);
      if (binPhi < 1) {
        binPhi = 1;
      }
      av += fHist->GetBinContent(binPhi) / (fHist->GetBinError(binPhi) * fHist->GetBinError(binPhi));
      errAv += 1. / (fHist->GetBinError(binPhi) * fHist->GetBinError(binPhi));
    } else {
      printf("[INFO] Reflected histo: only the point at +pi/2 used to evaluate baseline");
    }
    av /= errAv;
    errAv = TMath::Sqrt(1. / errAv);
    printf("[RESULT] Average fBaseline: %.3f +- %.3f \n", av, errAv);
    fBaseline = av;
    fErrBaseline = errAv;

    if (fShiftBaselineUp) {
      fBaseline += fErrBaseline;
      printf("[INFO] Shift baseline up of its statistical uncertainty");
    }

    if (fShiftBaselineDown) {
      fBaseline -= fErrBaseline;
      printf("[INFO] Shift baseline down of its statistical uncertainty");
    }

    return fBaseline;
  }

  // fix the baseline to the weighted average of the points passed through the function SetPointsForBaseline()
  if (fFixBase == 3) {
    if (fNbasleinePoints == 0) {
      printf("[ERROR] No baseline points set for the baseline evaluation, SetPointsForBaseline(Int_t nBaselinePoints, Double_t* valsBaseline). Returning -1");
      return -1;
    }
    Double_t errAv = 0., av = 0.;
    for (int i = 0; i < fNbasleinePoints; i++) {
      av += fHist->GetBinContent(fBinsBaseline[i]) / (fHist->GetBinError(fBinsBaseline[i]) * fHist->GetBinError(fBinsBaseline[i]));
      errAv += 1. / (fHist->GetBinError(fBinsBaseline[i]) * fHist->GetBinError(fBinsBaseline[i]));
    }
    av /= errAv;
    errAv = TMath::Sqrt(1. / errAv);
    printf("[RESULT] Average fBaseline: %.3f +- %.3f \n", av, errAv);
    fBaseline = av;
    fErrBaseline = errAv;

    if (fShiftBaselineUp) {
      fBaseline += fErrBaseline;
      printf("[INFO] Shift baseline up of its statistical uncertainty \n");
    }

    if (fShiftBaselineDown) {
      fBaseline -= fErrBaseline;
      printf("[INFO] Shift baseline down of its statistical uncertainty \n");
    }

    return fBaseline;
  }

  if (fFixBase == 4) {
    fBaseline = calculateBaseline(fHist, fIsTotal); // TODO: add the option for total range/ reflected range to pass in input
    fErrBaseline = calculateBaselineError(fHist, fIsTotal);

    if (fShiftBaselineUp) {
      fBaseline += fErrBaseline;
      printf("[INFO] Shift baseline up of its statistical uncertainty \n");
    }

    if (fShiftBaselineDown) {
      fBaseline -= fErrBaseline;
      printf("[INFO] Shift baseline down of its statistical uncertainty \n");
    }

    return fBaseline;
  }

  Printf("[ERROR] DhCorrelationFitter::FindBaseline - WRONG BASELINE OPTION SET. Returning -1");
  return -1.;
}

void DhCorrelationFitter::fitBaselineWv2()
{

  fBaseTransvReg = new TF1("fBaseTransvReg", [](const double* x, const double* p) {
    double const xx = x[0]; // x value
    if ((xx >= -TMath::Pi()/2 && xx <= -3*TMath::Pi()/8) || (xx >= 3*TMath::Pi()/8 && xx <= 5*TMath::Pi()/8) || (xx >= 11*TMath::Pi()/8 && xx <= 3*TMath::Pi()/2)) {
        // Gaussian example: p[0] = amplitude, p[1] = mean, p[2] = sigma
        return p[0]*(1+2*p[1]*p[2]*TMath::Cos(2*xx));
    }
    return 0.; }, -TMath::Pi() / 2, 3 * TMath::Pi() / 2, 3); // Function valid for [0,10], with 3 parameters

  fBaseTransvReg->FixParameter(1, fv2AssocPart);
  fBaseTransvReg->FixParameter(2, fv2Dmeson);

  TFitResultPtr const rFit = fHist->Fit(fBaseTransvReg, "RIMES", "", fMinCorr, fMaxCorr);
  fBaseline = fBaseTransvReg->GetParameter(0);
}

void DhCorrelationFitter::calculateYieldsAboveBaseline()
{

  fNSyieldBinCount = 0;
  fErrNSyieldBinCount = 0;
  fASyieldBinCount = 0;
  fErrASyieldBinCount = 0;
  std::cout << "[RESULT] Baseline: " << fBaseline << " +- " << fErrBaseline << std::endl;
  Int_t binMinNS = fHist->FindBin(-1.5); // slightly more than -pi/2
  if (binMinNS < 1) {
    binMinNS = 1; // with this, it is ok even in the case of a reflected fHist (range 0 - pi)
  }
  Int_t const binMaxNS = 6;  // fHist -> FindBin(1.5); // slightly less than +pi/2
  Int_t const binMinAS = 11; // fHist -> FindBin(1.6); // slightly more than +pi/2
  Int_t binMaxAS = 16;       // fHist -> FindBin(3.14+1.5); // slightly less than +3pi/2
  if (binMaxAS > fHist->GetNbinsX()) {
    binMaxAS = fHist->GetNbinsX(); // with this, it is ok even in the case of a reflected fHist (range 0 - pi)
  }
  std::cout << "N bins : " << fHist->GetNbinsX() << std::endl;
  std::cout << "binMinNS : " << binMinNS << std::endl;
  std::cout << "binMaxNS : " << binMaxNS << std::endl;
  std::cout << "binMinAS : " << binMinAS << std::endl;
  std::cout << "binMaxAS : " << binMaxAS << std::endl;
  // Near Side Yield from bin counting
  for (Int_t bmNS = binMinNS; bmNS <= binMaxNS; bmNS++) {
    fNSyieldBinCount += 2 * (fHist->GetBinContent(bmNS) - fBaseline) * fHist->GetBinWidth(bmNS);
    fErrNSyieldBinCount += 4 * (fHist->GetBinError(bmNS) * fHist->GetBinError(bmNS)) * fHist->GetBinWidth(bmNS) * fHist->GetBinWidth(bmNS);
  }
  fErrNSyieldBinCount = TMath::Sqrt(fErrNSyieldBinCount);

  // Away Side Yield from bin counting
  for (Int_t bmAS = binMinAS; bmAS <= binMaxAS; bmAS++) {
    fASyieldBinCount += 2 * (fHist->GetBinContent(bmAS) - fBaseline) * fHist->GetBinWidth(bmAS);
    fErrASyieldBinCount += 4 * (fHist->GetBinError(bmAS) * fHist->GetBinError(bmAS)) * fHist->GetBinWidth(bmAS) * fHist->GetBinWidth(bmAS);
  }
  fErrASyieldBinCount = TMath::Sqrt(fErrASyieldBinCount);

  printf("[RESULT] Bin counting results: NS Yield = %.3f +- %.3f \n[RESULT] Bin counting results: AS Yield: %.3f +- %.3f \n", fNSyieldBinCount, fErrNSyieldBinCount, fASyieldBinCount, fErrASyieldBinCount);
}

Double_t DhCorrelationFitter::calculateBaseline(TH1F*& histo, Bool_t totalRange)
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
    printf("[INFO] Using total deltaPhi range \n");
    // baseline evaluated considering: the two first points, the last two points and four points in the middle (corresponding to the outer points)
    if (nBinsPhi >= 32) {
      printf("[INFO] nBinsPhi >= 32 \n");
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
      printf("[INFO] nBinsPhi < 32 \n");
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
    printf("[INFO] Using reflected deltaPhi range \n");
    // baseline evaluated using the 4 middle points in the transverese region
    if (nBinsPhi >= 16) {
      printf("[INFO] 4 central points in the transverse region for baseline \n");
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

  return baseline;
}

Double_t DhCorrelationFitter::calculateBaselineError(TH1F*& histo, Bool_t totalRange)
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

  return errBaseline;
}

void DhCorrelationFitter::setSingleTermsForDrawing(Bool_t draw)
{
  Double_t* par = nullptr;
  if (fTypeOfFitFunc == 1 || fTypeOfFitFunc == 2 || fTypeOfFitFunc == 5) {
    par = new Double_t[7];
  } else if (fTypeOfFitFunc == 3 || fTypeOfFitFunc == 6) {
    par = new Double_t[4];
  } else if (fTypeOfFitFunc == 4) {
    par = new Double_t[8];
  } else if (fTypeOfFitFunc == 7) {
    par = new Double_t[9];
  } else {
    Printf("[ERROR] DhCorrelationFitter::SetSingleTermsForDrawing, wrong type of function");
    return;
  }

  if (fTypeOfFitFunc == 3 || fTypeOfFitFunc == 6) {
    fFit->GetParameters(par);
    fFit->SetLineWidth(4);

    fPed->SetParameter(0, par[0]);
    fPed->SetLineColor(6); // pink
    fPed->SetLineStyle(9);
    fPed->SetLineWidth(4);

    fGausAS->SetParameter(0, par[1]);
    fGausAS->SetParameter(1, par[2]);
    fGausAS->SetParameter(2, par[3]);

    fGausAS->SetLineStyle(9);
    fGausAS->SetLineColor(kGreen);
    fGausAS->SetLineWidth(4);

    auto* pvStatTests1 = new TPaveText(0.51, 0.58, 0.85, 0.90, "NDC");
    pvStatTests1->SetFillStyle(0);
    pvStatTests1->SetTextSize(0.045);
    pvStatTests1->SetBorderSize(0);
    // TText *t0, *t1, *t2, *t3;
    // t0 = pvStatTests1->AddText(0., 1.00, Form("#chi^{2}/ndf = %.1f/%d ", fFit->GetChisquare(), fFit->GetNDF()));
    // t1 = pvStatTests1->AddText(0., 0.80, Form("Ped = %.3f#pm%.3f ", fBaseline, fErrBaseline));
    // t2 = pvStatTests1->AddText(0., 0.65, Form("AS Y = %.3f#pm%.3f ", fFit->GetParameter("AS Y"), fFit->GetParError(fFit->GetParNumber("AS Y"))));
    // t3 = pvStatTests1->AddText(0., 0.50, Form("AS #sigma = %.3f#pm%.3f ", fFit->GetParameter("AS #sigma"), fFit->GetParError(fFit->GetParNumber("AS #sigma"))));

    if (draw) {
      fFit->Draw("same");
      fPed->Draw("same");
      fGausAS->Draw("same");
      pvStatTests1->Draw("same");
    }
  } else if (fTypeOfFitFunc == 4) {
    fFit->GetParameters(par);
    fFit->SetLineWidth(4);

    fPed->SetParameter(0, par[0]);
    fPed->SetLineColor(6); // pink
    fPed->SetLineStyle(9);
    fPed->SetLineWidth(4);

    fGausNS->SetParameter(0, par[1]);
    fGausNS->SetParameter(1, par[2]);
    fGausNS->SetParameter(2, par[3]);
    fGausNS->SetParameter(3, par[7]);
    fGausAS->SetParameter(0, par[4]);
    fGausAS->SetParameter(1, par[5]);
    fGausAS->SetParameter(2, par[6]);

    fGausNS->SetLineStyle(9);
    fGausNS->SetLineColor(kBlue);
    fGausAS->SetLineStyle(9);
    fGausAS->SetLineColor(kGreen);
    fGausNS->SetLineWidth(4);
    fGausAS->SetLineWidth(4);

    auto* pvStatTests1 = new TPaveText(0.51, 0.58, 0.85, 0.90, "NDC");
    pvStatTests1->SetFillStyle(0);
    pvStatTests1->SetTextSize(0.045);
    pvStatTests1->SetBorderSize(0);
    // TText *t0, *t1, *t2, *t3, *t4, *t5, *t6;
    // t0 = pvStatTests1->AddText(0., 1.00, Form("#chi^{2}/ndf = %.1f/%d ", fFit->GetChisquare(), fFit->GetNDF()));
    // t1 = pvStatTests1->AddText(0., 0.80, Form("Ped = %.3f#pm%.3f ", fBaseline, fErrBaseline));
    // t2 = pvStatTests1->AddText(0., 0.65, Form("NS Y = %.3f#pm%.3f ", fFit->GetParameter("NS Y"), fFit->GetParError(fFit->GetParNumber("NS Y"))));
    // t3 = pvStatTests1->AddText(0., 0.50, Form("NS #sigma = %.3f#pm%.3f ", fFit->GetParameter("NS #sigma"), fFit->GetParError(fFit->GetParNumber("NS #sigma"))));
    // t4 = pvStatTests1->AddText(0., 0.35, Form("AS Y = %.3f#pm%.3f ", fFit->GetParameter("AS Y"), fFit->GetParError(fFit->GetParNumber("AS Y"))));
    // t5 = pvStatTests1->AddText(0., 0.20, Form("AS #sigma = %.3f#pm%.3f ", fFit->GetParameter("AS #sigma"), fFit->GetParError(fFit->GetParNumber("AS #sigma"))));
    // t6 = pvStatTests1->AddText(0., 0.05, Form("#beta = %.3f#pm%.3f ", fFit->GetParameter("NS shape par"), fFit->GetParError(fFit->GetParNumber("NS shape par"))));

    if (draw) {
      fFit->Draw("same");
      fPed->Draw("same");
      fGausAS->Draw("same");
      fGausNS->Draw("same");
      pvStatTests1->Draw("same");
    }
  } else if (fTypeOfFitFunc == 7) {
    fFit->GetParameters(par);
    fFit->SetLineWidth(4);

    fBaseTransvReg->SetLineColor(15);
    fBaseTransvReg->SetLineStyle(9);
    fBaseTransvReg->SetLineWidth(4);

    fPed->SetParameter(0, par[0]);
    fPed->SetParameter(1, par[7]);
    fPed->SetParameter(2, par[8]);
    fPed->SetLineColor(6); // pink
    fPed->SetLineStyle(9);
    fPed->SetLineWidth(4);

    fGausNS->SetParameter(0, par[1]);
    fGausNS->SetParameter(1, par[2]);
    fGausNS->SetParameter(2, par[3]);
    // fGausNS -> SetParameter(3, par[7]);
    fGausAS->SetParameter(0, par[4]);
    fGausAS->SetParameter(1, par[5]);
    fGausAS->SetParameter(2, par[6]);

    fGausNS->SetLineStyle(9);
    fGausNS->SetLineColor(kBlue);
    fGausAS->SetLineStyle(9);
    fGausAS->SetLineColor(kGreen);
    fGausNS->SetLineWidth(4);
    fGausAS->SetLineWidth(4);

    auto* pvStatTests1 = new TPaveText(0.51, 0.58, 0.85, 0.90, "NDC");
    pvStatTests1->SetFillStyle(0);
    pvStatTests1->SetTextSize(0.045);
    pvStatTests1->SetBorderSize(0);
    // TText *t0, *t1, *t2, *t3, *t4, *t5, *t6, *t7, *t8;
    // t0 = pvStatTests1->AddText(0., 1.00, Form("#chi^{2}/ndf = %.1f/%d ", fFit->GetChisquare(), fFit->GetNDF()));
    // t2 = pvStatTests1->AddText(0., 0.80, Form("NS Y = %.3f#pm%.3f ", fFit->GetParameter("NS Y"), fFit->GetParError(fFit->GetParNumber("NS Y"))));
    // t3 = pvStatTests1->AddText(0., 0.65, Form("NS #sigma = %.3f#pm%.3f ", fFit->GetParameter("NS #sigma"), fFit->GetParError(fFit->GetParNumber("NS #sigma"))));
    // t4 = pvStatTests1->AddText(0., 0.50, Form("AS Y = %.3f#pm%.3f ", fFit->GetParameter("AS Y"), fFit->GetParError(fFit->GetParNumber("AS Y"))));
    // t5 = pvStatTests1->AddText(0., 0.35, Form("AS #sigma = %.3f#pm%.3f ", fFit->GetParameter("AS #sigma"), fFit->GetParError(fFit->GetParNumber("AS #sigma"))));
    // t6 = pvStatTests1 -> AddText(0., 0.20, Form("#beta = %.3f#pm%.3f ", fFit -> GetParameter("NS shape par"), fFit -> GetParError(fFit->GetParNumber("NS shape par"))));

    auto* pvStatTests2 = new TPaveText(0.51, 0.28, 0.85, 0.60, "NDC");
    pvStatTests2->SetFillStyle(0);
    pvStatTests2->SetTextSize(0.045);
    pvStatTests2->SetBorderSize(0);
    // t1 = pvStatTests2->AddText(0., 1.00, Form("Ped = %.3f#pm%.3f ", fFit->GetParameter("ped"), fErrBaseline /*fFit -> GetParError(fFit->GetParNumber("ped")*/));
    // t7 = pvStatTests2->AddText(0., 0.65, Form("v_{2}^{hadron} = %.3f#pm%.3f ", fFit->GetParameter("v_{2} hadron"), fFit->GetParError(fFit->GetParNumber("v_{2} hadron"))));
    // t8 = pvStatTests2->AddText(0., 0.35, Form("v_{2}^{D} = %.3f#pm%.3f ", fFit->GetParameter("v_{2} D meson"), fFit->GetParError(fFit->GetParNumber("v_{2} D meson"))));

    if (draw) {
      fFit->Draw("same");
      fPed->Draw("same");
      fBaseTransvReg->Draw("same");
      fGausAS->Draw("same");
      fGausNS->Draw("same");
      pvStatTests1->Draw("same");
      pvStatTests2->Draw("same");
    }
  } else {
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

    auto* pvStatTests1 = new TPaveText(0.51, 0.58, 0.85, 0.90, "NDC");
    pvStatTests1->SetFillStyle(0);
    pvStatTests1->SetTextSize(0.045);
    pvStatTests1->SetBorderSize(0);
    // TText *t0, *t1, *t2, *t3, *t4, *t5, *t6;
    // t0 = pvStatTests1->AddText(0., 1.00, Form("#chi^{2}/ndf = %.1f/%d ", fFit->GetChisquare(), fFit->GetNDF()));
    // t1 = pvStatTests1->AddText(0., 0.80, Form("Ped = %.3f#pm%.3f ", fBaseline, fErrBaseline));
    // t2 = pvStatTests1->AddText(0., 0.65, Form("NS Y = %.3f#pm%.3f ", fFit->GetParameter("NS Y"), fFit->GetParError(fFit->GetParNumber("NS Y"))));
    // t3 = pvStatTests1->AddText(0., 0.50, Form("NS #sigma = %.3f#pm%.3f ", fFit->GetParameter("NS #sigma"), fFit->GetParError(fFit->GetParNumber("NS #sigma"))));
    // t4 = pvStatTests1->AddText(0., 0.35, Form("AS Y = %.3f#pm%.3f ", fFit->GetParameter("AS Y"), fFit->GetParError(fFit->GetParNumber("AS Y"))));
    // t5 = pvStatTests1->AddText(0., 0.20, Form("AS #sigma = %.3f#pm%.3f ", fFit->GetParameter("AS #sigma"), fFit->GetParError(fFit->GetParNumber("AS #sigma"))));

    if (draw) {
      fFit->Draw("same");
      fPed->Draw("same");
      fGausAS->Draw("same");
      fGausNS->Draw("same");
      pvStatTests1->Draw("same");
    }
  }
  delete[] par;
}
