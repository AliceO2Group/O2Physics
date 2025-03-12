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

/// \file DhCorrelationFitter.h
/// \brief Class to perform the azimuthal correlation fit
/// \author Samuele Cattaruzzi <samuele.cattaruzzi@cern.ch>
/// \author Swapnesh Santosh Khade <swapnesh.santosh.khade@cern.ch>

#ifndef PWGHF_HFC_MACROS_DHCORRELATIONFITTER_H_
#define PWGHF_HFC_MACROS_DHCORRELATIONFITTER_H_

#include <cstdio>

#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLatex.h>

class DhCorrelationFitter
{

 public:
  enum FunctionType { kConstwoGaus = 1,
                      kTwoGausPeriodicity = 2 };

  /// Constructors
  DhCorrelationFitter();
  DhCorrelationFitter(TH1F* histoToFit, Double_t min, Double_t max);
  virtual ~DhCorrelationFitter();
  DhCorrelationFitter(const DhCorrelationFitter& source);
  DhCorrelationFitter& operator=(const DhCorrelationFitter& cfit);

  /// Setters
  void SetHistoIsReflected(Bool_t isrefl) { fIsReflected = isrefl; }
  void SetFuncType(FunctionType fitType) { fTypeOfFitFunc = fitType; }
  void SetFixBaseline(Int_t fixBase) { fFixBase = fixBase; }
  void SetFixMean(Int_t fixMean) { fFixMean = fixMean; }
  void SetPtRanges(Double_t PtCandMin, Double_t PtCandMax, Double_t PtAssocMin, Double_t PtAssocMax)
  {
    fMinCandPt = PtCandMin;
    fMaxCandPt = PtCandMax;
    fMinAssoPt = PtAssocMin;
    fMaxAssoPt = PtAssocMax;
  }
  void SetExternalValsAndBounds(Int_t nPars, Double_t* vals, Double_t* lowBounds, Double_t* uppBounds);
  void SetPointsForBaseline(Int_t nBaselinePoints, Int_t* binsBaseline);

  /// Functions for fitting
  void Fitting(Bool_t drawSplitTerm = kTRUE, Bool_t useExternalPars = kFALSE);
  void SetFitFunction();
  void CalculateYieldsAboveBaseline();
  void SetSingleTermsForDrawing(Bool_t draw);
  Double_t FindBaseline();

  /// Getters
  Double_t GetNSSigma() { return fFit->GetParameter("NS #sigma"); } // TODO: case kConstThreeGausPeriodicity
  Double_t GetASSigma() { return fFit->GetParameter("AS #sigma"); } // TODO: case kConstThreeGausPeriodicity
  Double_t GetNSYield() { return fFit->GetParameter("NS Y"); }
  Double_t GetASYield() { return fFit->GetParameter("AS Y"); }
  Double_t GetBeta() { return fFit->GetParameter(7); }
  Double_t GetPedestal() { return fBaseline; }
  Double_t Getv2hadron() { return fFit->GetParameter("v_{2} hadron"); }
  Double_t Getv2Dmeson() { return fFit->GetParameter("v_{2} D meson"); }
  Double_t GetNSSigmaError() { return fFit->GetParError(fFit->GetParNumber("NS #sigma")); } // TODO: case kConstThreeGausPeriodicity
  Double_t GetASSigmaError() { return fFit->GetParError(fFit->GetParNumber("AS #sigma")); } // TODO: case kConstThreeGausPeriodicityAS
  Double_t GetNSYieldError() { return fFit->GetParError(fFit->GetParNumber("NS Y")); }
  Double_t GetASYieldError() { return fFit->GetParError(fFit->GetParNumber("AS Y")); }
  Double_t GetBetaError() { return fFit->GetParError(7); }
  Double_t GetPedestalError() { return fErrBaseline; }
  Double_t Getv2hadronError() { return fFit->GetParError(fFit->GetParNumber("v_{2} hadron")); }
  Double_t Getv2DmesonError() { return fFit->GetParError(fFit->GetParNumber("v_{2} D meson")); }
  Double_t GetBinCountingNSYield() { return fNSyieldBinCount; }
  Double_t GetBinCountingASYield() { return fASyieldBinCount; }
  Double_t GetBinCountingNSYieldErr() { return fErrNSyieldBinCount; }
  Double_t GetBinCountingASYieldErr() { return fErrASyieldBinCount; }
  TF1* GetFitFunction()
  {
    if (!fFit) {
      printf("[ERROR] DhCorrelationFitter::GetFitFunction, No fit function");
      return NULL;
    }
    return fFit;
  }

 private:
  TH1F* fHist; // 1D azimuthal correlation histogram

  TF1* fFit;    // Total fit function
  TF1* fGausNS; // Near-Side (NS) Gaussian
  TF1* fGausAS; // Away-Side (AS) Gaussian
  TF1* fPed;    // Baseline function

  Bool_t fIsReflected;
  Bool_t fUseExternalPars; // To use external fit parameters initial values and bounds

  FunctionType fTypeOfFitFunc; // Type of fit function

  Int_t fFixBase;         // Fix baseline or keep it as free parameter (see Fitting())
  Int_t fFixMean;         // Fix means or keep them as free parameters (see Fitting())
  Int_t fNpars;           // Number of parameters of the fit function
  Int_t fNbasleinePoints; // Number of points passed as inputs for the baseline estimation

  Int_t* fBinsBaseline; // Bin values for calculating the baseline

  Double_t fMinCorr;            // Min value of correlation histogram
  Double_t fMaxCorr;            // Max value of correlation histogram
  Double_t fMinCandPt;          // Min pt value of candidate
  Double_t fMaxCandPt;          // Max pt value of candidate
  Double_t fMinAssoPt;          // Min pt value of assoc. particles
  Double_t fMaxAssoPt;          // Max pt value of assoc. particles
  Double_t fBaseline;           // Baseline value
  Double_t fErrBaseline;        // Baseline error
  Double_t fNSyieldBinCount;    // NS Yield from bin counting
  Double_t fErrNSyieldBinCount; // NS Yield error from bin counting
  Double_t fASyieldBinCount;    // AS Yield from bin counting
  Double_t fErrASyieldBinCount; // AS Yield error from bin counting

  Double_t* fExtParsVals;      // Fit parameters initial values
  Double_t* fExtParsLowBounds; // Fit parameters lower bounds
  Double_t* fExtParsUppBounds; // Fit parameters upper bounds
};

#endif // PWGHF_HFC_MACROS_DHCORRELATIONFITTER_H_
