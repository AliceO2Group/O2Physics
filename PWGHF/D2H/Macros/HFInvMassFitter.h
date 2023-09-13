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

/// \file HFInvMassFitter.h
/// \brief HFInvMassFitter class
///
/// \author Zhen Zhang <zhenz@cern.ch>
/// \author Mingyu Zhang <mingyu.zang@cern.ch>
/// \author Xinye Peng  <xinye.peng@cern.ch>
/// \author Biao Zhang <biao.zhang@cern.ch>

#ifndef PWGHF_D2H_MACROS_HFINVMASSFITTER_H_
#define PWGHF_D2H_MACROS_HFINVMASSFITTER_H_

#include <string> // std::string

#include <RooWorkspace.h>
#include <TCanvas.h>
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TH1.h>
#include <TH1F.h>
#include <TNamed.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TVirtualFitter.h>

using namespace RooFit;

class TF1;
class TH1F;

class HFInvMassFitter : public TNamed
{
 public:
  enum TypeOfBkgPdf {
    Expo = 0,
    Poly1 = 1,
    Poly2 = 2,
    Pow = 3,
    PowExpo = 4,
    Poly3 = 5,
    NoBkg = 6
  };
  enum TypeOfSgnPdf {
    SingleGaus = 0,
    DoubleGaus = 1,
    DoubleGausSigmaRatioPar = 2,
    GausSec = 3
  };
  enum TypeOfReflPdf {
    SingleGausRefl = 0,
    DoubleGausRefl = 1,
    Poly3Refl = 2,
    Poly6Refl = 3
  };
  HFInvMassFitter();
  HFInvMassFitter(const TH1F* histoToFit, Double_t minValue, Double_t maxValue, Int_t fitTypeBkg = Expo, Int_t fitTypeSgn = SingleGaus);
  ~HFInvMassFitter();
  void setHistogramForFit(const TH1F* histoToFit)
  {
    if (mHistoInvMass) {
      delete mHistoInvMass;
    }
    mHistoInvMass = reinterpret_cast<TH1F*>(histoToFit->Clone("mHistoInvMass"));
    mHistoInvMass->SetDirectory(0);
  }
  void setUseLikelihoodFit() { mFitOption = "L,E"; }
  void setUseChi2Fit() { mFitOption = "Chi2"; }
  void setFitOption(TString opt) { mFitOption = opt.Data(); }
  RooAbsPdf* createBackgroundFitFunction(RooWorkspace* w1);
  RooAbsPdf* createSignalFitFunction(RooWorkspace* w1);
  RooAbsPdf* createReflectionFitFunction(RooWorkspace* w1);

  void setFitRange(Double_t minValue, Double_t maxValue)
  {
    mMinMass = minValue;
    mMaxMass = maxValue;
  }
  void setFitFunctions(Int_t fitTypeBkg, Int_t fitTypeSgn)
  {
    mTypeOfBkgPdf = fitTypeBkg;
    mTypeOfSgnPdf = fitTypeSgn;
  }
  void setSigmaLimit(Double_t sigmaValue, Double_t sigmaLimit)
  {
    mSigmaValue = sigmaValue;
    mParamSgn = sigmaLimit;
  }
  void setParticlePdgMass(Double_t mass) { mMassParticle = mass; }
  Double_t getParticlePdgMass() { return mMassParticle; }
  void setInitialGaussianMean(Double_t mean)
  {
    mMass = mean;
    mSecMass = mean;
  }
  void setInitialGaussianSigma(Double_t sigma)
  {
    mSigmaSgn = sigma;
    mSecSigma = sigma;
  }
  void setInitialSecondGaussianSigma(Double_t sigma) { mSigmaSgnDoubleGaus = sigma; }
  void setInitialFracDoubleGaus(Double_t frac) { mFracDoubleGaus = frac; }
  void setInitialRatioDoubleGausSigma(Double_t fracSigma) { mRatioDoubleGausSigma = fracSigma; }
  void setFixGaussianMean(Double_t mean)
  {
    setInitialGaussianMean(mean);
    mFixedMean = kTRUE;
  }
  void setBoundGaussianMean(Double_t mean, Double_t meanLowLimit, Double_t meanUpLimit)
  {
    if (mean < meanLowLimit ||
        mean > meanUpLimit) {
      cout << "Invalid Gaussian mean limmit!" << endl;
    }
    setInitialGaussianMean(mean);
    mMassLowLimit = meanLowLimit;
    mMassUpLimit = meanUpLimit;
    mBoundMean = kTRUE;
  }
  void setBoundReflGausMean(Double_t mean, Double_t meanLowLimit, Double_t meanUpLimit)
  {
    if (mean < meanLowLimit ||
        mean > meanUpLimit) {
      cout << "Invalid Gaussian mean limmit for reflection!" << endl;
    }
    setInitialGaussianMean(mean);
    mMassReflLowLimit = meanLowLimit;
    mMassReflUpLimit = meanUpLimit;
    mBoundReflMean = kTRUE;
  }
  void setFixGaussianSigma(Double_t sigma)
  {
    setInitialGaussianSigma(sigma);
    mFixedSigma = kTRUE;
  }
  void setBoundGausSigma(Double_t sigma, Double_t sigmaLimit)
  {
    setInitialGaussianSigma(sigma);
    setSigmaLimit(sigma, sigmaLimit);
    mBoundSigma = kTRUE;
  }
  void setFixSecondGaussianSigma(Double_t sigma)
  {
    if (mTypeOfSgnPdf != DoubleGaus) {
      cout << "Fit type should be 2Gaus!" << endl;
    }
    setInitialSecondGaussianSigma(sigma);
    mFixedSigmaDoubleGaus = kTRUE;
  }
  void setFixFrac2Gaus(Double_t frac)
  {
    if (mTypeOfSgnPdf != DoubleGaus &&
        mTypeOfSgnPdf != DoubleGausSigmaRatioPar) {
      cout << "Fit type should be 2Gaus or 2GausSigmaRatio!" << endl;
    }
    setInitialFracDoubleGaus(frac);
    mFixedFracDoubleGaus = kTRUE;
  }
  void setFixRatioToGausSigma(Double_t sigmaFrac)
  {
    if (mTypeOfSgnPdf != DoubleGausSigmaRatioPar) {
      cout << "Fit type should be set to k2GausSigmaRatioPar!" << endl;
    }
    setInitialRatioDoubleGausSigma(sigmaFrac);
    mFixedRatioDoubleGausSigma = kTRUE;
  }
  void setFixSignalYield(Double_t yield) { mFixedRawYield = yield; }
  void setNumberOfSigmaForSidebands(Double_t numberOfSigma) { mNSigmaForSidebands = numberOfSigma; }
  void plotBkg(RooAbsPdf* mFunc);
  void plotRefl(RooAbsPdf* mFunc);
  void setReflFuncFixed();
  void doFit(Bool_t draw = kTRUE);
  void setInitialReflOverSgn(Double_t reflOverSgn) { mReflOverSgn = reflOverSgn; }
  void setFixReflOverSgn(Double_t reflOverSgn)
  {
    setInitialReflOverSgn(reflOverSgn);
    mFixReflOverSgn = kTRUE;
  }
  void setTemplateReflections(const TH1* histoRefl, Int_t fitTypeRefl = DoubleGaus)
  {
    if (!histoRefl) {
      mEnableReflections = kFALSE;
    }
    mHistoTemplateRefl = reinterpret_cast<TH1F*>(histoRefl->Clone("mHistoTemplateRefl"));
  }
  Double_t getChiSquareOverNDF() const { return mChiSquareOverNdf; }
  Double_t getRawYield() const { return mRawYield; }
  Double_t getRawYieldError() const { return mRawYieldErr; }
  Double_t getBkgYield() const { return mBkgYield; }
  Double_t getBkgYieldError() const { return mBkgYieldErr; }
  Double_t getSignificance() const { return mSignificance; }
  Double_t getSignificanceError() const { return mSignificanceErr; }
  Double_t getMean() const { return mRooMeanSgn->getVal(); }
  Double_t getMeanUncertainty() const { return mRooMeanSgn->getError(); }
  Double_t getSigma() const { return mRooSigmaSgn->getVal(); }
  Double_t getSigmaUncertainty() const { return mRooSigmaSgn->getError(); }
  Double_t getReflOverSig() const
  {
    if (mReflPdf) {
      return mReflOverSgn;
    } else {
      return 0;
    }
  }
  void calculateSignal(Double_t& signal, Double_t& signalErr) const;
  void calculateBackground(Double_t& bkg, Double_t& bkgErr) const;
  void calculateSignificance(Double_t& significance, Double_t& significanceErr) const;
  void checkForSignal(Double_t& estimatedSignal);
  void drawFit(TVirtualPad* c, Int_t writeFitInfo = 2);
  void drawResidual(TVirtualPad* c);
  void drawReflection(TVirtualPad* c);

 private:
  HFInvMassFitter(const HFInvMassFitter& source);
  HFInvMassFitter& operator=(const HFInvMassFitter& source);
  void fillWorkspace(RooWorkspace& w);

  TH1F* mHistoInvMass; // histogram to fit
  TString mFitOption;
  Double_t mMinMass;                 // lower mass limit
  Double_t mMaxMass;                 // upper mass limit
  Int_t mTypeOfBkgPdf;               // background fit function
  Int_t mTypeOfSgnPdf;               // signal fit function
  Int_t mTypeOfReflPdf;              // reflection fit function
  Double_t mMassParticle;            // pdg value of particle mass
  Double_t mMass;                    /// signal gaussian mean value
  Double_t mMassLowLimit;            /// lower limit of the allowed mass range
  Double_t mMassUpLimit;             /// upper limit of the allowed mass range
  Double_t mMassReflLowLimit;        /// lower limit of the allowed mass range for reflection
  Double_t mMassReflUpLimit;         /// upper limit of the allowed mass range for reflection
  Double_t mSecMass;                 /// Second peak mean value
  Double_t mMassErr;                 /// uncertainty on signal gaussian mean value
  Double_t mSigmaSgn;                /// signal gaussian sigma
  Double_t mSecSigma;                /// Second peak gaussian sigma
  Int_t mNSigmaForSidebands;         /// number of sigmas to veto the signal peak
  Int_t mNSigmaForSgn;               /// number of sigmas to veto the signal peak
  Double_t mSigmaSgnErr;             /// uncertainty on signal gaussian sigma
  Double_t mSigmaSgnDoubleGaus;      /// signal 2gaussian sigma
  Double_t mFixedMean;               /// switch for fix mean of gaussian
  Bool_t mBoundMean;                 /// switch for bound mean of guassian
  Bool_t mBoundReflMean;             /// switch for bound mean of guassian for reflection
  Bool_t mFixedSigma;                /// fix sigma or not
  Bool_t mFixedSigmaDoubleGaus;      /// fix sigma of 2gaussian or not
  Bool_t mBoundSigma;                /// set bound sigma or not
  Double_t mSigmaValue;              /// value of sigma
  Double_t mParamSgn;                /// +/- range variation of bound Sigma of gaussian in %
  Double_t mFracDoubleGaus;          /// initialization for fraction of 2nd gaussian in case of k2Gaus or k2GausSigmaRatioPar
  Double_t mFixedRawYield;           /// initialization for raw yield
  Bool_t mFixedFracDoubleGaus;       /// switch for fixed fraction of 2nd gaussian in case of k2Gaus or k2GausSigmaRatioPar
  Double_t mRatioDoubleGausSigma;    /// initialization for ratio between two gaussian sigmas in case of k2GausSigmaRatioPar
  Bool_t mFixedRatioDoubleGausSigma; /// switch for fixed ratio between two gaussian sigmas in case of k2GausSigmaRatioPar
  Double_t mReflOverSgn;             /// reflection/signal
  Bool_t mEnableReflections;         /// flag use/not use reflections
  Double_t mRawYield;                /// signal gaussian integral
  Double_t mRawYieldErr;             /// err on signal gaussian integral
  Double_t mBkgYield;                /// background
  Double_t mBkgYieldErr;             /// err on background
  Double_t mSignificance;            /// significance
  Double_t mSignificanceErr;         /// err on significance
  Double_t mChiSquareOverNdf;        /// chi2/ndf
  Bool_t mFixReflOverSgn;            /// switch for fix refl/signal
  RooRealVar* mRooMeanSgn;           /// mean for gaussian of signal
  RooRealVar* mRooSigmaSgn;          /// sigma for gaussian of signal
  RooAbsPdf* mSgnPdf;                /// signal fit function
  RooAbsPdf* mBkgPdf;                /// background fit function
  RooAbsPdf* mReflPdf;               /// reflection fit function
  RooRealVar* mRooNSgn;              /// total Signal fit function integral
  RooRealVar* mRooNBkg;              /// total background fit function integral
  RooRealVar* mRooNRefl;             /// total reflection fit function integral
  RooAbsPdf* mTotalPdf;              /// total fit function
  RooPlot* mInvMassFrame;            /// frame of mass
  RooPlot* mReflFrame;               /// reflection frame
  RooPlot* mReflOnlyFrame;           /// reflection frame plot on reflection only
  RooPlot* mResidualFrame;           /// residual frame
  RooPlot* mResidualFrameForCalulation;
  RooRealVar* mass;         /// mass
  RooWorkspace* mWorkspace; /// workspace
  Double_t mIntegralHisto;  /// integral of histogram to fit
  Double_t mIntegralBkg;    /// integral of background fit function
  Double_t mIntegralSgn;    /// integral of signal fit function
  TH1F* mHistoTemplateRefl; /// reflection histogram

  ClassDef(HFInvMassFitter, 1);
};

#endif // PWGHF_D2H_MACROS_HFINVMASSFITTER_H_
