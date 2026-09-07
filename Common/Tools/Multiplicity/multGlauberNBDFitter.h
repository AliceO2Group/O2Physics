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

/// \file multGlauberNBDFitter.h
/// \brief Helper class to perform Glauber fits
/// \author David Dobrigkeit Chinellato

#ifndef COMMON_TOOLS_MULTIPLICITY_MULTGLAUBERNBDFITTER_H_
#define COMMON_TOOLS_MULTIPLICITY_MULTGLAUBERNBDFITTER_H_

#include <Framework/Logger.h>

#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TNamed.h>
#include <TProfile.h>
#include <TString.h>

#include <Rtypes.h>

class multGlauberNBDFitter : public TNamed
{
 public:
  // basic functionality
  multGlauberNBDFitter();
  explicit multGlauberNBDFitter(const char* name, const char* title = "Glauber+NBD fitter");
  ~multGlauberNBDFitter() override;

  enum class AncestorMode {
    Truncated = 0, // Nancestors = (int) truncation
    Rounded,       // Nancestors = rounded to nearest integer
    Continuous     // Nancestors = kept as float, uses ContinuousNBD
  };

  enum class NBDFitterMode {
    Glauber = 0, // Uses Npart & Ncoll from TGlauberMC
    Trento       // Uses trento entropy
  };

  enum class FitPar {
    mu = 0,
    k,
    f,
    norm,
    dMu
  };

  enum class NBDPar {
    p = 0,
    k
  };

  template <typename TEnum>
  static constexpr int Index(TEnum e)
  {
    return static_cast<int>(e);
  }

  // Master fitter function
  double GlauberProbDistrib(const double* x, const double* par);
  double TrentoProbDistrib(const double* x, const double* par);

  void InitAncestor();

  // Do Fit: where everything happens
  bool DoFit();
  bool DoGlauberFit();
  bool DoTrentoFit();

  // Set input characteristics: the 2D plot with Npart, Nanc
  bool SetNpartNcollCorrelation(TH2* hNpNc);
  // Set input characteristics: the 1D plot with trento entropy
  bool SetNSources(TH1* hNSources);

  // Set main input to be fitted (the V0M distribution)
  bool SetInputV0M(TH1* hV0M);

  // Interface to get funtions if asked to
  TF1* GetNBD();
  TF1* GetGlauberNBD();

  // Helper
  bool InitializeNpNc();
  void InitGlauberNBD(const float mu, const float k, const float f, const float norm);
  void InitTrentoNBD(const float mu = 45, const float k = 1.5, const float norm = 100);

  // Interface for debug
  void SetAncestorMode(AncestorMode lAncMode = AncestorMode::Truncated) { fAncestorMode = lAncMode; }
  void SetAncestorMode(int lAncMode = Index(AncestorMode::Truncated))
  {
    if (lAncMode < Index(AncestorMode::Truncated) || lAncMode > Index(AncestorMode::Continuous)) {
      LOG(fatal) << "Invalid ancestor mode: " << lAncMode;
      return;
    }
    fAncestorMode = static_cast<AncestorMode>(lAncMode);
  }
  int GetAncestorMode() { return static_cast<int>(fAncestorMode); }
  TH1D* GetAncestorHistogram() { return fhNanc; }

  // Interface to set vals
  void SetMu(const double lVal) { fMu = lVal; }
  void Setk(const double lVal) { fk = lVal; }
  void Setf(const double lVal) { ff = lVal; }
  void SetNorm(const double lVal) { fnorm = lVal; }

  // Interface to get vals
  double GetMu() { return fMu; }
  double Getk() { return fk; }
  double Getf() { return ff; }
  double GetNorm() { return fnorm; }

  void SetFitRange(const double lMin, const double lMax);
  void SetFitOptions(const TString& lOpt);
  void SetFitNpx(const long lNpx);

  // For ancestor mode 2
  double ContinuousNBD(const double n, const double mu, const double k);

  // For estimating Npart, Ncoll in multiplicity bins
  // also viable: eccentricity, impact parameter, ancestor cross-check plot
  void CalculateAvNpNc(TProfile* lNPartProf, TProfile* lNCollProf, TH2F* lNPart2DPlot, TH2F* lNColl2DPlot, TH1F* hPercentileMap, double lLoRange = -1, double lHiRange = -1, TH3D* lNpNcEcc = nullptr, TH2F* lEcc2DPlot = nullptr, TH3D* lNpNcB = nullptr, TH2F* lB2DPlot = nullptr, TH2F* lNancestor2DPlot = nullptr, double fProbabilityCutoff = -1);

 private:
  // This function serves as the (analytical) NBD
  TF1* fNBD;

  // This function is the key fitting function
  TF1* fGlauberNBD;
  TF1* fTrentoNBD;

  // This is the fitting mode that we use
  NBDFitterMode fNBDFitterMode;

  // Reference histo
  TH1D* fhNanc;    // basic ancestor distribution
  TH2* fhNpNc;     // correlation between Npart and Ncoll
  TH1* fhNSources; // Trento entropy
  TH1* fhV0M;      // basic ancestor distribution

  // Fitting utilities
  bool ffChanged;
  double fCurrentf;

  // 0: truncation, 1: rounding, 2: analytical continuation
  AncestorMode fAncestorMode;

  // Buffer for (Npart, Ncoll) pairs in memory
  double* fNpart;
  double* fNcoll;
  long* fContent;
  long fNNpNcPairs; // number of pairs to use
  long fMaxNpNcPairs;

  // The actual output: mu, k, f, norm
  double fMu;
  double fdMu; // variable mu option
  double fk;
  double ff;
  double fnorm;

  TString fFitOptions;
  long fFitNpx;

  ClassDefOverride(multGlauberNBDFitter, 1);
};
#endif // COMMON_TOOLS_MULTIPLICITY_MULTGLAUBERNBDFITTER_H_
