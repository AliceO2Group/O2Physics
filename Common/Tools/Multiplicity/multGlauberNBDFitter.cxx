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

/// \file multGlauberNBDFitter.cxx
/// \brief Helper class to perform Glauber fits
/// \author David Dobrigkeit Chinellato

/**********************************************
 *
 *   Class meant to do Glauber+NBD fits
 *
 *   This class makes full use of analytical
 *   properties of the Negative Binomial Function
 *
 *   Only the Glauber component is MC, while
 *   the NBD is evaluated probabilistically
 *
 *  - bugs, comments, suggestions, complaints?
 *  - Feel free to write to:
 *     david.dobrigkeit.chinellato@cern.ch
 *
 **********************************************/

#include "multGlauberNBDFitter.h"

#include <Framework/Logger.h>

#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TNamed.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TVirtualFitter.h>

#include <Rtypes.h>

#include <cmath>

ClassImp(multGlauberNBDFitter);

multGlauberNBDFitter::multGlauberNBDFitter() : multGlauberNBDFitter("multGlauberNBDFitter") {}
multGlauberNBDFitter::multGlauberNBDFitter(const char* name, const char* title) : TNamed(name, title),
                                                                                  fNBD(new TF1("fNBD", "ROOT::Math::negative_binomial_pdf(x,[0],[1])", 0, 45000)),
                                                                                  fGlauberNBD(nullptr),
                                                                                  fTrentoNBD(nullptr),
                                                                                  fNBDFitterMode(NBDFitterMode::Glauber),
                                                                                  fhNanc(nullptr),
                                                                                  fhNpNc(nullptr),
                                                                                  fhNSources(nullptr),
                                                                                  fhV0M(nullptr),
                                                                                  ffChanged(true),
                                                                                  fCurrentf(-1),
                                                                                  fAncestorMode(AncestorMode::Continuous),
                                                                                  fNpart(nullptr),
                                                                                  fNcoll(nullptr),
                                                                                  fContent(nullptr),
                                                                                  fNNpNcPairs(-1),
                                                                                  fMaxNpNcPairs(1000000),
                                                                                  fMu(45),
                                                                                  fdMu(0.0),
                                                                                  fk(1.5),
                                                                                  ff(0.8),
                                                                                  fnorm(100),
                                                                                  fFitOptions("R0"),
                                                                                  fFitNpx(5000)
{
  // Named constructor
  fNpart = new double[fMaxNpNcPairs];
  fNcoll = new double[fMaxNpNcPairs];
  fContent = new int64_t[fMaxNpNcPairs];

  // NBD
  fNBD->SetNpx(45000);

  // master function
  InitGlauberNBD(fMu, fk, ff, fnorm);
}

//________________________________________________________________
multGlauberNBDFitter::~multGlauberNBDFitter()
{
  // Destructor
  delete fNBD;
  delete fhNanc;
  delete fhNpNc;
  delete[] fNpart;
  delete[] fNcoll;
  delete[] fContent;
}

void multGlauberNBDFitter::InitGlauberNBD(const float mu, const float k, const float f, const float norm)
{
  fNBDFitterMode = NBDFitterMode::Glauber;
  fGlauberNBD = new TF1("fGlauberNBD", this, &multGlauberNBDFitter::GlauberProbDistrib,
                        0, 50000, 5, "multGlauberNBDFitter", "GlauberProbDistrib");
  fGlauberNBD->SetParameter(Index(FitPar::mu), mu);
  fGlauberNBD->SetParameter(Index(FitPar::k), k);
  fGlauberNBD->SetParameter(Index(FitPar::f), f);
  fGlauberNBD->SetParameter(Index(FitPar::norm), norm);

  fGlauberNBD->SetParName(Index(FitPar::mu), "mu");
  fGlauberNBD->SetParName(Index(FitPar::k), "k");
  fGlauberNBD->SetParName(Index(FitPar::f), "f");
  fGlauberNBD->SetParName(Index(FitPar::norm), "norm");
  fGlauberNBD->SetParName(Index(FitPar::dMu), "dMu/dNanc");
}

void multGlauberNBDFitter::InitTrentoNBD(const float mu, const float k, const float norm)
{
  fNBDFitterMode = NBDFitterMode::Trento;
  fGlauberNBD = nullptr;
  fTrentoNBD = new TF1("fTrentoNBD", this, &multGlauberNBDFitter::TrentoProbDistrib,
                       0, 50000, 5, "multGlauberNBDFitter", "TrentoProbDistrib");
  fTrentoNBD->SetParameter(Index(FitPar::mu), mu);
  fTrentoNBD->SetParameter(Index(FitPar::k), k);
  fTrentoNBD->SetParameter(Index(FitPar::norm), norm);
  fTrentoNBD->FixParameter(Index(FitPar::dMu), 0);

  fTrentoNBD->SetParName(Index(FitPar::mu), "mu");
  fTrentoNBD->SetParName(Index(FitPar::k), "k");
  fTrentoNBD->SetParName(Index(FitPar::f), "f");
  fTrentoNBD->SetParName(Index(FitPar::norm), "norm");
  fTrentoNBD->SetParName(Index(FitPar::dMu), "dMu/dNanc");
}

//______________________________________________________
double multGlauberNBDFitter::GlauberProbDistrib(const double* x, const double* par)
// Master fitter function
{
  double lMultValue = x[0];
  double lProbability = 0.0;
  ffChanged = true;
  static constexpr double Almost0 = 1.e-13;
  // Comment this line in order to make the code evaluate Nancestor all the time
  if (std::abs(fCurrentf - par[Index(FitPar::f)]) < Almost0) {
    ffChanged = false;
  }

  //______________________________________________________
  // Recalculate the ancestor distribution in case f changed
  if (ffChanged) {
    fCurrentf = par[Index(FitPar::f)];
    fhNanc->Reset();

    for (int ibin = 0; ibin < fNNpNcPairs; ibin++) {
      double lOption0 = static_cast<int>(fNpart[ibin] * par[Index(FitPar::f)] + fNcoll[ibin] * (1.0 - par[Index(FitPar::f)]));
      double lOption1 = std::floor(fNpart[ibin] * par[Index(FitPar::f)] + fNcoll[ibin] * (1.0 - par[Index(FitPar::f)]) + 0.5);
      double lOption2 = (fNpart[ibin] * par[Index(FitPar::f)] + fNcoll[ibin] * (1.0 - par[Index(FitPar::f)]));
      if (fAncestorMode == AncestorMode::Truncated) {
        fhNanc->Fill(lOption0, fContent[ibin]);
      }
      if (fAncestorMode == AncestorMode::Rounded) {
        fhNanc->Fill(lOption1, fContent[ibin]);
      }
      if (fAncestorMode == AncestorMode::Continuous) {
        fhNanc->Fill(lOption2, fContent[ibin]);
      }
    }
    if (fhNanc->Integral() < 1) {
      LOG(info) << "ERROR: ANCESTOR HISTOGRAM EMPTY";
      LOG(info) << "Will not do anything. Call InitNpNc if you want to plot without fitting";
      return 0;
    }
    fhNanc->Scale(1. / fhNanc->Integral());
  }
  //______________________________________________________
  // Actually evaluate function
  int lStartBin = fhNanc->FindBin(0.0) + 1;
  for (int64_t iNanc = lStartBin; iNanc < fhNanc->GetNbinsX() + 1; ++iNanc) {
    double lNancestors = fhNanc->GetBinCenter(iNanc);
    double lNancestorCount = fhNanc->GetBinContent(iNanc);

    // allow for variable mu in case requested
    double lThisMu = lNancestors * (par[Index(FitPar::mu)] + par[Index(FitPar::dMu)] * lNancestors);
    double lThisk = lNancestors * par[Index(FitPar::k)];
    double lpval = std::pow(1.0 + lThisMu / lThisk, -1);
    fNBD->SetParameter(Index(NBDPar::k), lThisk);
    fNBD->SetParameter(Index(NBDPar::p), lpval);

    double lMult = 0.0;
    static constexpr double MultTolerance = 1e-6;
    if (lMultValue > MultTolerance) {
      lMult = fAncestorMode != AncestorMode::Continuous ? fNBD->Eval(lMultValue) : ContinuousNBD(lMultValue, lThisMu, lThisk);
    }
    lProbability += lNancestorCount * lMult;
  }
  //______________________________________________________
  return par[Index(FitPar::norm)] * lProbability;
}

double multGlauberNBDFitter::TrentoProbDistrib(const double* x, const double* par)
// Master fitter function
{
  double lMultValue = x[0];
  double lProbability = 0.0;
  //______________________________________________________
  // Actually ealuate function
  for (int64_t iNSrc = 1; iNSrc < fhNSources->GetNbinsX() + 1; ++iNSrc) {
    double lNsources = fhNSources->GetBinCenter(iNSrc);
    double lThisMu = lNsources * par[Index(FitPar::mu)];
    double lThisk = lNsources * par[Index(FitPar::k)];
    double lpval = std::pow(1 + lThisMu / lThisk, -1);
    fNBD->SetParameter(Index(NBDPar::k), lThisk);
    fNBD->SetParameter(Index(NBDPar::p), lpval);
    double lMult = fNBD->Eval(lMultValue);
    lProbability += fhNSources->GetBinContent(fhNSources->FindBin(iNSrc)) * lMult;
  }
  //______________________________________________________
  return par[Index(FitPar::norm)] * lProbability;
}

//________________________________________________________________
bool multGlauberNBDFitter::SetNpartNcollCorrelation(TH2* hNpNc)
{
  bool lReturnValue = true;
  if (hNpNc) {
    fhNpNc = hNpNc;
  } else {
    lReturnValue = false;
  }
  return lReturnValue;
}

//________________________________________________________________
bool multGlauberNBDFitter::SetNSources(TH1* hNSources)
{
  bool lReturnValue = true;
  if (hNSources) {
    fhNSources = hNSources;
  } else {
    lReturnValue = false;
  }
  return lReturnValue;
}

//________________________________________________________________
bool multGlauberNBDFitter::SetInputV0M(TH1* hV0M)
{
  bool lReturnValue = true;
  if (hV0M) {
    fhV0M = hV0M;
  } else {
    lReturnValue = false;
  }
  return lReturnValue;
}

//________________________________________________________________
TF1* multGlauberNBDFitter::GetNBD()
{
  return fNBD;
}

//________________________________________________________________
TF1* multGlauberNBDFitter::GetGlauberNBD()
{
  return fGlauberNBD;
}

TF1* multGlauberNBDFitter::GetTrentoNBD()
{
  return fTrentoNBD;
}

//________________________________________________________________
void multGlauberNBDFitter::SetFitRange(const double lMin, const double lMax)
{
  if (fGlauberNBD) {
    fGlauberNBD->SetRange(lMin, lMax);
  }
  if (fTrentoNBD) {
    fTrentoNBD->SetRange(lMin, lMax);
  }
}

//________________________________________________________________
void multGlauberNBDFitter::SetFitOptions(const TString& lOpt)
{
  fFitOptions = std::move(lOpt);
}

//________________________________________________________________
void multGlauberNBDFitter::SetFitNpx(const int64_t lNpx)
{
  fFitNpx = lNpx;
}

//________________________________________________________________
void multGlauberNBDFitter::InitAncestor()
{
  if (!fhNanc) {
    fhNanc = new TH1D("fhNanc", "", fAncestorMode == AncestorMode::Continuous ? 10000 : 1000, -0.5, 999.5);
  }
}

//________________________________________________________________
bool multGlauberNBDFitter::DoFit()
{
  // Try very hard, please
  TVirtualFitter::SetMaxIterations(5000000);
  auto timer = new TStopwatch();
  timer->Start(true);

  bool lReturnValue = false;
  switch (fNBDFitterMode) {
    case NBDFitterMode::Glauber:
      LOG(info) << "Doing Glauber fit!!";
      lReturnValue = DoGlauberFit();
      break;
    case NBDFitterMode::Trento:
      LOG(info) << "Doing Trento fit!!";
      lReturnValue = DoTrentoFit();
      break;

    default:
      break;
  }

  timer->Stop();
  double lTotalTime = timer->RealTime();
  if (lReturnValue) {
    LOG(info) << "---> Fitting succeeded after " << lTotalTime << " seconds";
  } else {
    LOG(info) << "---> Fitting failed after " << lTotalTime << " seconds";
  }
  return lReturnValue;
}

bool multGlauberNBDFitter::DoGlauberFit()
{
  InitAncestor();
  if (!InitializeNpNc()) {
    LOG(info) << "---> Initialization of Npart x Ncoll correlation info failed!";
    return false;
  }

  if (fAncestorMode == AncestorMode::Truncated) {
    LOG(info) << "---> Config: Nancestors will be truncated";
  }
  if (fAncestorMode == AncestorMode::Rounded) {
    LOG(info) << "---> Config: Nancestors will be rounded";
  }
  if (fAncestorMode == AncestorMode::Continuous) {
    LOG(info) << "---> Config: Nancestors will be taken as float";
  }
  LOG(info) << "---> Now fitting, please wait...";

  fGlauberNBD->SetNpx(fFitNpx);
  TFitResultPtr fitptr;
  fFitOptions.Append("S");
  fitptr = fhV0M->Fit("fGlauberNBD", fFitOptions.Data());

  fMu = fGlauberNBD->GetParameter(Index(FitPar::mu));
  fk = fGlauberNBD->GetParameter(Index(FitPar::k));
  ff = fGlauberNBD->GetParameter(Index(FitPar::f));
  fnorm = fGlauberNBD->GetParameter(Index(FitPar::norm));
  fdMu = fGlauberNBD->GetParameter(Index(FitPar::dMu));
  return fitptr.Get()->IsValid();
}

bool multGlauberNBDFitter::DoTrentoFit()
{
  fTrentoNBD->SetNpx(fFitNpx);
  TFitResultPtr fitptr;
  fFitOptions.Append("S");
  fitptr = fhV0M->Fit("fTrentoNBD", fFitOptions.Data());
  fMu = fTrentoNBD->GetParameter(Index(FitPar::mu));
  fk = fTrentoNBD->GetParameter(Index(FitPar::k));
  fnorm = fTrentoNBD->GetParameter(Index(FitPar::norm));
  return fitptr.Get()->IsValid();
}

//________________________________________________________________
bool multGlauberNBDFitter::InitializeNpNc()
{
  // This function initializes fhNpNc
  // Warning: X == Npart, Y == Ncoll
  bool lReturnValue = false;
  if (fhNpNc) {
    fNNpNcPairs = 0;
    // Sweep all allowed values of Npart, Ncoll; find counters
    static constexpr int NBinsX = 500;
    static constexpr int NBinsY = 3000;
    for (int xbin = 1; xbin < NBinsX; xbin++) {
      for (int ybin = 1; ybin < NBinsY; ybin++) {
        if (fhNpNc->GetBinContent(fhNpNc->FindBin(xbin, ybin)) != 0) {
          fNpart[fNNpNcPairs] = xbin;
          fNcoll[fNNpNcPairs] = ybin;
          fContent[fNNpNcPairs] = static_cast<int64_t>(fhNpNc->GetBinContent(fhNpNc->FindBin(xbin, ybin)));
          fNNpNcPairs++;
        }
      }
    }
    LOG(info) << "Initialized with number of (Npart, Ncoll) pairs: " << fNNpNcPairs;
    lReturnValue = true;
  } else {
    LOG(info) << "Failed to initialize! Please provide input histogram with (Npart, Ncoll) info!";
    LOG(info) << "Please remember to call SetNpartNcollCorrelation before doing fit!";
  }
  return lReturnValue;
}

//________________________________________________________________
double multGlauberNBDFitter::ContinuousNBD(const double n, const double mu, const double k)
{
  // Adaptation of the negative binomial distribution
  // for non-integer arguments: analytical continuation
  //
  // This function would actually also be fine with integers;
  // in fact it is equivalent to that if 'n' is typecast as
  // an integer prior to use

  double F{};
  double f{};

  static constexpr double NumStabilityThreshold = 100.0;
  if (n + k > NumStabilityThreshold) {
    // log method for handling large numbers
    F = std::lgamma(n + k) - std::lgamma(n + 1.) - std::lgamma(k);
    f = n * std::log(mu / k) - (n + k) * std::log1p(mu / k);
    F = F + f;
    F = std::exp(F);
  } else {
    F = std::tgamma(n + k) / (std::tgamma(n + 1.) * std::tgamma(k));
    f = n * std::log(mu / k) - (n + k) * std::log1p(mu / k);
    f = std::exp(f);
    F *= f;
  }
  return F;
}

void multGlauberNBDFitter::CalculateAvNpNc(TProfile* lNPartProf, TProfile* lNCollProf, TH2F* lNPart2DPlot, TH2F* lNColl2DPlot, TH1F* hPercentileMap, double lLoRange, double lHiRange, TH3D* lNpNcEcc, TH2F* lEcc2DPlot, TH3D* lNpNcB, TH2F* lB2DPlot, TH2F* lNancestor2DPlot, double fProbabilityCutoff)
{
  LOG(info) << "Calculating <Npart>, <Ncoll> in centrality bins...";
  LOG(info) << "Range to calculate: " << lLoRange << " to " << lHiRange;

  LOG(info) << "Acquiring values from the fit function...";

  fMu = fGlauberNBD->GetParameter(Index(FitPar::mu));
  fk = fGlauberNBD->GetParameter(Index(FitPar::k));
  ff = fGlauberNBD->GetParameter(Index(FitPar::f));
  fnorm = fGlauberNBD->GetParameter(Index(FitPar::norm));
  fdMu = fGlauberNBD->GetParameter(Index(FitPar::dMu));

  LOG(info) << "Please inspect now: ";
  LOG(info) << "Glauber NBD mu ............: " << fMu;
  LOG(info) << "Glauber NBD k .............: " << fk;
  LOG(info) << "Glauber NBD f .............: " << ff;
  LOG(info) << "Glauber NBD norm ..........: " << fnorm;
  LOG(info) << "Glauber NBD dmu/dNanc .....: " << fdMu;

  // 2-fold nested loop:
  //  + looping over all Nancestor combinations
  //  + looping over all possible final multiplicities
  //  ^---> final product already multiplicity-binned

  //______________________________________________________
  if (lLoRange < -1 && lHiRange < -1) {
    fGlauberNBD->GetRange(lLoRange, lHiRange);
  }
  // bypass to zero
  for (int ibin = 0; ibin < fNNpNcPairs; ibin++) {
    static constexpr int LogInterval = 200;
    if (ibin % LogInterval == 0) {
      LOG(info) << "At NpNc pair #" << ibin << " of " << fNNpNcPairs << "...";
    }
    double lNAncestors0 = static_cast<int>(fNpart[ibin] * ff + fNcoll[ibin] * (1.0 - ff));
    double lNAncestors1 = std::floor(fNpart[ibin] * ff + fNcoll[ibin] * (1.0 - ff) + 0.5);
    double lNAncestors2 = (fNpart[ibin] * ff + fNcoll[ibin] * (1.0 - ff));

    // define ancestors officially
    double lNancestors = lNAncestors0;
    if (fAncestorMode == AncestorMode::Rounded) {
      lNancestors = lNAncestors1;
    }
    if (fAncestorMode == AncestorMode::Continuous) {
      lNancestors = lNAncestors2;
    }

    // eccentricity handling
    TH1D* hEccentricity = nullptr;
    if (lNpNcEcc) {
      // locate the histogram that corresponds to the eccentricity distribution in this NpNc pair
      lNpNcEcc->GetXaxis()->SetRange(lNpNcEcc->GetXaxis()->FindBin(fNpart[ibin]), lNpNcEcc->GetXaxis()->FindBin(fNpart[ibin]));
      lNpNcEcc->GetYaxis()->SetRange(lNpNcEcc->GetYaxis()->FindBin(fNcoll[ibin]), lNpNcEcc->GetYaxis()->FindBin(fNcoll[ibin]));
      hEccentricity = dynamic_cast<TH1D*>(lNpNcEcc->Project3D("z"));
      hEccentricity->SetName(Form("hEccentricity_%i", ibin));

      // normalize into unitary fractions
      double eccIntegral = hEccentricity->Integral(1, hEccentricity->GetNbinsX() + 1);
      static constexpr double EccTolerance = 1e-6;
      if (eccIntegral > EccTolerance) { // no counts
        hEccentricity->Scale(1. / eccIntegral);
      } else {
        hEccentricity->Scale(0.0);
      }
    }

    // impact parameter handling
    TH1D* hImpactParameter = nullptr;
    if (lNpNcB) {
      // locate the histogram that corresponds to the eccentricity distribution in this NpNc pair
      lNpNcB->GetXaxis()->SetRange(lNpNcB->GetXaxis()->FindBin(fNpart[ibin]), lNpNcB->GetXaxis()->FindBin(fNpart[ibin]));
      lNpNcB->GetYaxis()->SetRange(lNpNcB->GetYaxis()->FindBin(fNcoll[ibin]), lNpNcB->GetYaxis()->FindBin(fNcoll[ibin]));
      hImpactParameter = dynamic_cast<TH1D*>(lNpNcB->Project3D("z"));
      hImpactParameter->SetName(Form("hImpactParameter_%i", ibin));

      // normalize into unitary fractions
      double bIntegral = hImpactParameter->Integral(1, hImpactParameter->GetNbinsX() + 1);
      static constexpr double ImpactParTolerance = 1e-6;
      if (bIntegral > ImpactParTolerance) { // no counts
        hImpactParameter->Scale(1. / bIntegral);
      } else {
        hImpactParameter->Scale(0.0);
      }
    }

    for (int64_t lMultValue = 1; lMultValue < lHiRange; lMultValue++) {
      double lNancestorCount = fContent[ibin];
      double lThisMu = lNancestors * fMu;
      double lThisk = lNancestors * fk;
      double lpval = std::pow(1 + lThisMu / lThisk, -1);
      fNBD->SetParameter(1, lThisk);
      fNBD->SetParameter(0, lpval);
      double lMult = 0.0;
      static constexpr double MultTolerance = 1e-6;
      if (lMultValue > MultTolerance) {
        lMult = fAncestorMode != AncestorMode::Continuous ? fNBD->Eval(lMultValue) : ContinuousNBD(lMultValue, lThisMu, lThisk);
      }
      double lProbability = lNancestorCount * lMult;

      if (lProbability < fProbabilityCutoff) {
        continue; // skip if probability of contributing too small
      }

      double lMultValueToFill = lMultValue;
      if (hPercentileMap) {
        lMultValueToFill = hPercentileMap->GetBinContent(hPercentileMap->FindBin(lMultValue));
      }
      lNPartProf->Fill(lMultValueToFill, fNpart[ibin], lProbability);
      lNCollProf->Fill(lMultValueToFill, fNcoll[ibin], lProbability);
      if (lNancestor2DPlot) {
        // fill cross-check histogram with lNancestorCount at lNancestors value
        lNancestor2DPlot->Fill(lMultValueToFill, lNancestors, lProbability);
      }
      if (lNPart2DPlot) {
        lNPart2DPlot->Fill(lMultValueToFill, fNpart[ibin], lProbability);
      }
      if (lNColl2DPlot) {
        lNColl2DPlot->Fill(lMultValueToFill, fNcoll[ibin], lProbability);
      }
      if (lNpNcEcc) {
        // collapse the entire eccentricity distribution for this combo
        for (int ib = 1; ib < hEccentricity->GetNbinsX() + 1; ib++) {
          lEcc2DPlot->Fill(lMultValueToFill, hEccentricity->GetBinCenter(ib), lProbability * hEccentricity->GetBinContent(ib));
        }
      }
      if (lNpNcB) {
        // collapse the entire impact parameter distribution for this combo
        for (int ib = 1; ib < hImpactParameter->GetNbinsX() + 1; ib++) {
          lB2DPlot->Fill(lMultValueToFill, hImpactParameter->GetBinCenter(ib), lProbability * hImpactParameter->GetBinContent(ib));
        }
      }
    }
  }
}
