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
#ifndef COMMON_TOOLS_MULTIPLICITY_MULTGLAUBERNBDFITTER_H_
#define COMMON_TOOLS_MULTIPLICITY_MULTGLAUBERNBDFITTER_H_

#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TNamed.h>
#include <TProfile.h>
#include <TString.h>

#include <Rtypes.h>
#include <RtypesCore.h>

class multGlauberNBDFitter : public TNamed
{

 public:
  // basic functionality
  multGlauberNBDFitter();
  explicit multGlauberNBDFitter(const char* name, const char* title = "Glauber+NBD fitter");
  ~multGlauberNBDFitter();

  // Master fitter function
  Double_t ProbDistrib(Double_t* x, Double_t* par);

  void InitAncestor();

  // Do Fit: where everything happens
  Bool_t DoFit();

  // Set input characteristics: the 2D plot with Npart, Nanc
  Bool_t SetNpartNcollCorrelation(TH2* hNpNc);

  // Set main input to be fitted (the V0M distribution)
  Bool_t SetInputV0M(TH1* hV0M);

  // Interface to get funtions if asked to
  TF1* GetNBD();
  TF1* GetGlauberNBD();

  // Helper
  Bool_t InitializeNpNc();

  // Interface for debug
  void SetAncestorMode(Int_t lAncMode = 0) { fAncestorMode = lAncMode; }
  Int_t GetAncestorMode() { return fAncestorMode; }
  TH1D* GetAncestorHistogram() { return fhNanc; }

  // Interface to set vals
  void SetMu(Double_t lVal) { fMu = lVal; }
  void Setk(Double_t lVal) { fk = lVal; }
  void Setf(Double_t lVal) { ff = lVal; }
  void SetNorm(Double_t lVal) { fnorm = lVal; }

  // Interface to get vals
  Double_t GetMu() { return fMu; }
  Double_t Getk() { return fk; }
  Double_t Getf() { return ff; }
  Double_t GetNorm() { return fnorm; }

  void SetFitRange(Double_t lMin, Double_t lMax);
  void SetFitOptions(TString lOpt);
  void SetFitNpx(Long_t lNpx);

  // For ancestor mode 2
  Double_t ContinuousNBD(Double_t n, Double_t mu, Double_t k);

  // For estimating Npart, Ncoll in multiplicity bins
  // also viable: eccentricity, impact parameter, ancestor cross-check plot
  void CalculateAvNpNc(TProfile* lNPartProf, TProfile* lNCollProf, TH2F* lNPart2DPlot, TH2F* lNColl2DPlot, TH1F* hPercentileMap, Double_t lLoRange = -1, Double_t lHiRange = -1, TH3D* lNpNcEcc = 0x0, TH2F* lEcc2DPlot = 0x0, TH3D* lNpNcB = 0x0, TH2F* lB2DPlot = 0x0, TH2F* lNancestor2DPlot = 0x0, Double_t fProbabilityCutoff = -1);

  // void    Print(Option_t *option="") const;

 private:
  // This function serves as the (analytical) NBD
  TF1* fNBD;

  // This function is the key fitting function
  TF1* fGlauberNBD;

  // Reference histo
  TH1D* fhNanc; // basic ancestor distribution
  TH2* fhNpNc;  // correlation between Npart and Ncoll
  TH1* fhV0M;   // basic ancestor distribution

  // Fitting utilities
  Bool_t ffChanged;
  Double_t fCurrentf;

  // 0: truncation, 1: rounding, 2: analytical continuation
  Int_t fAncestorMode;

  // Buffer for (Npart, Ncoll) pairs in memory
  Double_t* fNpart;
  Double_t* fNcoll;
  Long_t* fContent;
  Long_t fNNpNcPairs; // number of pairs to use
  Long_t fMaxNpNcPairs;

  // The actual output: mu, k, f, norm
  Double_t fMu;
  Double_t fdMu; // variable mu option
  Double_t fk;
  Double_t ff;
  Double_t fnorm;

  TString fFitOptions;
  Long_t fFitNpx;

  ClassDef(multGlauberNBDFitter, 1);
};
#endif // COMMON_TOOLS_MULTIPLICITY_MULTGLAUBERNBDFITTER_H_
