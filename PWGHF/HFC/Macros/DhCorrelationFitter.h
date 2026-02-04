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

#include <TF1.h>
#include <TH1.h>

#include <RtypesCore.h>

#include <cstdio>

class DhCorrelationFitter
{

 public:
  enum FunctionType { kConstwoGaus = 1,
                      kTwoGausPeriodicity = 2,
                      kSingleGaus = 3,
                      kGenGaus = 4,
                      kVonMises = 5,
                      kSingleVonMises = 6,
                      kTwoGausPeriodicityPlusV2modulation = 7 };

  /// Constructors
  DhCorrelationFitter();
  DhCorrelationFitter(TH1F* histoToFit, Double_t min, Double_t max);
  virtual ~DhCorrelationFitter();
  DhCorrelationFitter(const DhCorrelationFitter& source);
  DhCorrelationFitter& operator=(const DhCorrelationFitter& cfit);

  /// Setters
  void setHistoIsReflected(Bool_t isrefl) { fIsReflected = isrefl; }
  void setFuncType(FunctionType fitType) { fTypeOfFitFunc = fitType; }
  void setFixBaseline(Int_t fixBase) { fFixBase = fixBase; }
  void setFixMean(Int_t fixMean) { fFixMean = fixMean; }
  void setPtRanges(Double_t ptCandMin, Double_t ptCandMax, Double_t ptAssocMin, Double_t ptAssocMax)
  {
    fMinCandPt = ptCandMin;
    fMaxCandPt = ptCandMax;
    fMinAssoPt = ptAssocMin;
    fMaxAssoPt = ptAssocMax;
  }
  void setExternalValsAndBounds(Int_t nPars, const Double_t* vals, const Double_t* lowBounds, const Double_t* uppBounds);
  void setPointsForBaseline(Int_t nBaselinePoints, const Int_t* binsBaseline);
  void setReflectedCorrHisto(Bool_t isReflected) { fIsTotal = !isReflected; }
  void setBaselineUpOrDown(Bool_t baseUp, Bool_t baseDown)
  {
    fShiftBaselineUp = baseUp;
    fShiftBaselineDown = baseDown;
  }
  void setv2(Double_t v2AssocPart, Double_t v2Dmeson)
  {
    fv2AssocPart = v2AssocPart;
    fv2Dmeson = v2Dmeson;
  }

  /// Functions for fitting
  void fitting(Bool_t drawSplitTerm = kTRUE, Bool_t useExternalPars = kFALSE);
  void setFitFunction();
  void calculateYieldsAboveBaseline();
  void fitBaselineWv2();
  Double_t calculateBaseline(TH1F*& histo, Bool_t totalRange = kTRUE);
  Double_t calculateBaselineError(TH1F*& histo, Bool_t totalRange = kTRUE);
  void setSingleTermsForDrawing(Bool_t draw);
  Double_t findBaseline();

  /// Getters
  Double_t getNsSigma() { return fFit->GetParameter("NS #sigma"); } // TODO: case kConstThreeGausPeriodicity
  Double_t getAsSigma() { return fFit->GetParameter("AS #sigma"); } // TODO: case kConstThreeGausPeriodicity
  Double_t getNsYield() { return fFit->GetParameter("NS Y"); }
  Double_t getAsYield() { return fFit->GetParameter("AS Y"); }
  Double_t getBeta() { return fFit->GetParameter(7); }
  [[nodiscard]] Double_t getPedestal() const { return fBaseline; }
  Double_t getv2hadron() { return fFit->GetParameter("v_{2} hadron"); }
  Double_t getv2Dmeson() { return fFit->GetParameter("v_{2} D meson"); }
  Double_t getNsSigmaError() { return fFit->GetParError(fFit->GetParNumber("NS #sigma")); } // TODO: case kConstThreeGausPeriodicity
  Double_t getAsSigmaError() { return fFit->GetParError(fFit->GetParNumber("AS #sigma")); } // TODO: case kConstThreeGausPeriodicityAS
  Double_t getNsYieldError() { return fFit->GetParError(fFit->GetParNumber("NS Y")); }
  Double_t getAsYieldError() { return fFit->GetParError(fFit->GetParNumber("AS Y")); }
  Double_t getBetaError() { return fFit->GetParError(7); }
  [[nodiscard]] Double_t getPedestalError() const { return fErrBaseline; }
  Double_t getv2hadronError() { return fFit->GetParError(fFit->GetParNumber("v_{2} hadron")); }
  Double_t getv2DmesonError() { return fFit->GetParError(fFit->GetParNumber("v_{2} D meson")); }
  [[nodiscard]] Double_t getBinCountingNsYield() const { return fNSyieldBinCount; }
  [[nodiscard]] Double_t getBinCountingAsYield() const { return fASyieldBinCount; }
  [[nodiscard]] Double_t getBinCountingNsYieldErr() const { return fErrNSyieldBinCount; }
  [[nodiscard]] Double_t getBinCountingAsYieldErr() const { return fErrASyieldBinCount; }
  TF1* getFitFunction()
  {
    if (fFit == nullptr) {
      printf("[ERROR] DhCorrelationFitter::GetFitFunction, No fit function");
      return nullptr;
    }
    return fFit;
  }

 private:
  TH1F* fHist; // 1D azimuthal correlation histogram

  TF1* fFit;           // Total fit function
  TF1* fGausNS;        // Near-Side (NS) Gaussian
  TF1* fGausAS;        // Away-Side (AS) Gaussian
  TF1* fPed;           // Baseline function
  TF1* fBaseTransvReg; // Baseline function with v2

  Bool_t fIsReflected;       // To use if reflected azimuthal correlation are given as input
  Bool_t fUseExternalPars;   // To use external fit parameters initial values and bounds
  Bool_t fShiftBaselineUp;   // To shift the baseline up of its statistical uncertainty
  Bool_t fShiftBaselineDown; // To shift baseline down of its statistical uncertainty
  Bool_t fIsTotal;           // Total range of 2*pi in the azimuthal correlation distribution

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
  Double_t fv2AssocPart;        // v2 associated particles
  Double_t fv2Dmeson;           // v2 of D mesons

  Double_t* fExtParsVals;      // Fit parameters initial values
  Double_t* fExtParsLowBounds; // Fit parameters lower bounds
  Double_t* fExtParsUppBounds; // Fit parameters upper bounds
};

#endif // PWGHF_HFC_MACROS_DHCORRELATIONFITTER_H_
