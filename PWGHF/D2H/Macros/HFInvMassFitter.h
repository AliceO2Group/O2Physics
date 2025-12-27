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
/// \author Oleksii Lubynets <oleksii.lubynets@cern.ch>
/// \author Phil Stahlhut <phil.lennart.stahlhut@cern.ch>

#ifndef PWGHF_D2H_MACROS_HFINVMASSFITTER_H_
#define PWGHF_D2H_MACROS_HFINVMASSFITTER_H_

#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <TF1.h>
#include <TH1.h>
#include <TNamed.h>
#include <TVirtualPad.h>

#include <Rtypes.h>
#include <RtypesCore.h>

#include <array>
#include <string>
#include <vector>

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
    NoBkg = 6,
    NTypesOfBkgPdf
  };
  std::array<std::string, NTypesOfBkgPdf> namesOfBkgPdf{"bkgFuncExpo", "bkgFuncPoly1", "bkgFuncPoly2", "bkgFuncPow", "bkgFuncPowExpo", "bkgFuncPoly3"};
  enum TypeOfSgnPdf {
    SingleGaus = 0,
    DoubleGaus = 1,
    DoubleGausSigmaRatioPar = 2,
    GausSec = 3,
    NTypesOfSgnPdf
  };
  enum TypeOfReflPdf {
    SingleGausRefl = 0,
    DoubleGausRefl = 1,
    Poly3Refl = 2,
    Poly6Refl = 3,
    NTypesOfReflPdf
  };
  std::array<std::string, NTypesOfReflPdf> namesOfReflPdf{"reflFuncGaus", "reflFuncDoubleGaus", "reflFuncPoly3", "reflFuncPoly6"};
  HFInvMassFitter() = delete;
  HFInvMassFitter(TH1* histoToFit, double minValue, double maxValue, int fitTypeBkg = Expo, int fitTypeSgn = SingleGaus);
  ~HFInvMassFitter() override;
  void setHistogramForFit(TH1* histoToFit);
  void setUseLikelihoodFit() { mFitOption = "L,E"; }
  void setUseChi2Fit() { mFitOption = "Chi2"; }
  void setFitOption(const std::string& opt) { mFitOption = opt; }
  RooAbsPdf* createBackgroundFitFunction(RooWorkspace* w1) const;
  RooAbsPdf* createSignalFitFunction(RooWorkspace* w1);
  RooAbsPdf* createReflectionFitFunction(RooWorkspace* w1) const;

  void setFitRange(double minValue, double maxValue);
  void setFitFunctions(int fitTypeBkg, int fitTypeSgn);
  void setSigmaLimit(double sigmaValue, double sigmaLimit);
  void setParticlePdgMass(double mass) { mMassParticle = mass; }
  [[nodiscard]] double getParticlePdgMass() const { return mMassParticle; }
  void setInitialGaussianMean(double mean);
  void setInitialGaussianSigma(double sigma);
  void setInitialSecondGaussianSigma(double sigma) { mSigmaSgnDoubleGaus = sigma; }
  void setInitialFracDoubleGaus(double frac) { mFracDoubleGaus = frac; }
  void setInitialRatioDoubleGausSigma(double fracSigma) { mRatioDoubleGausSigma = fracSigma; }
  void setFixGaussianMean(double mean);
  void setBoundGaussianMean(double mean, double meanLowLimit, double meanUpLimit);
  void setBoundReflGausMean(double mean, double meanLowLimit, double meanUpLimit);
  void setFixGaussianSigma(double sigma);
  void setBoundGausSigma(double sigma, double sigmaLimit);
  void setFixSecondGaussianSigma(double sigma);
  void setFixFrac2Gaus(double frac);
  void setFixRatioToGausSigma(double sigmaFrac);
  void setFixSignalYield(double yield) { mFixedRawYield = yield; }
  void setNumberOfSigmaForSidebands(double numberOfSigma) { mNSigmaForSidebands = numberOfSigma; }
  void plotBkg(RooAbsPdf* mFunc, Color_t color = kRed);
  void plotRefl(RooAbsPdf* mFunc);
  void setReflFuncFixed();
  void doFit();
  void setInitialReflOverSgn(double reflOverSgn) { mReflOverSgn = reflOverSgn; }
  void setFixReflOverSgn(double reflOverSgn);
  void setTemplateReflections(TH1* histoRefl);
  void setDrawBgPrefit(bool value = true) { mDrawBgPrefit = value; }
  void setHighlightPeakRegion(bool value = true) { mHighlightPeakRegion = value; }
  [[nodiscard]] double getChiSquareOverNDFTotal() const { return mChiSquareOverNdfTotal; }
  [[nodiscard]] double getChiSquareOverNDFBkg() const { return mChiSquareOverNdfBkg; }
  [[nodiscard]] double getRawYield() const { return mRawYield; }
  [[nodiscard]] double getRawYieldError() const { return mRawYieldErr; }
  [[nodiscard]] double getRawYieldCounted() const { return mRawYieldCounted; }
  [[nodiscard]] double getRawYieldCountedError() const { return mRawYieldCountedErr; }
  [[nodiscard]] double getBkgYield() const { return mBkgYield; }
  [[nodiscard]] double getBkgYieldError() const { return mBkgYieldErr; }
  [[nodiscard]] double getSignificance() const { return mSignificance; }
  [[nodiscard]] double getSignificanceError() const { return mSignificanceErr; }
  [[nodiscard]] double getMean() const { return mRooMeanSgn->getVal(); }
  [[nodiscard]] double getMeanUncertainty() const { return mRooMeanSgn->getError(); }
  [[nodiscard]] double getSigma() const { return mRooSigmaSgn->getVal(); }
  [[nodiscard]] double getSigmaUncertainty() const { return mRooSigmaSgn->getError(); }
  [[nodiscard]] double getSecSigma() const { return mRooSecSigmaSgn->getVal(); }
  [[nodiscard]] double getSecSigmaUncertainty() const { return mRooSecSigmaSgn->getError(); }
  [[nodiscard]] double getFracDoubleGaus() const { return mRooFracDoubleGaus->getVal(); }
  [[nodiscard]] double getFracDoubleGausUncertainty() const { return mRooFracDoubleGaus->getError(); }
  [[nodiscard]] double getReflOverSig() const { return mReflPdf != nullptr ? mReflOverSgn : 0.; }
  void calculateSignal(double& signal, double& signalErr) const;
  void countSignal(double& signal, double& signalErr) const;
  void calculateBackground(double& bkg, double& bkgErr) const;
  void calculateSignificance(double& significance, double& significanceErr) const;
  void checkForSignal(double& estimatedSignal);
  void calculateFitToDataRatio() const;
  void drawFit(TVirtualPad* c, const std::vector<std::string>& plotLabels, bool writeParInfo = true);
  void drawResidual(TVirtualPad* c);
  void drawRatio(TVirtualPad* c);
  void drawReflection(TVirtualPad* c);

 private:
  HFInvMassFitter(const HFInvMassFitter& source);
  HFInvMassFitter& operator=(const HFInvMassFitter& source);
  void fillWorkspace(RooWorkspace& w) const;
  void highlightPeakRegion(const RooPlot* plot, Color_t color = kGray + 1, Width_t width = 1, Style_t style = 2) const;

  TH1* mHistoInvMass; // histogram to fit
  std::string mFitOption;
  double mMinMass;                 // lower mass limit
  double mMaxMass;                 // upper mass limit
  int mTypeOfBkgPdf;               // background fit function
  int mTypeOfSgnPdf;               // signal fit function
  int mTypeOfReflPdf;              // reflection fit function
  double mMassParticle;            // pdg value of particle mass
  double mMass;                    /// signal gaussian mean value
  double mMassLowLimit;            /// lower limit of the allowed mass range
  double mMassUpLimit;             /// upper limit of the allowed mass range
  double mMassReflLowLimit;        /// lower limit of the allowed mass range for reflection
  double mMassReflUpLimit;         /// upper limit of the allowed mass range for reflection
  double mSecMass;                 /// Second peak mean value
  double mSigmaSgn;                /// signal gaussian sigma
  double mSecSigma;                /// Second peak gaussian sigma
  int mNSigmaForSidebands;         /// number of sigmas to veto the signal peak
  int mNSigmaForSgn;               /// number of sigmas to veto the signal peak
  double mSigmaSgnErr;             /// uncertainty on signal gaussian sigma
  double mSigmaSgnDoubleGaus;      /// signal 2gaussian sigma
  bool mFixedMean;                 /// switch for fix mean of gaussian
  bool mBoundMean;                 /// switch for bound mean of guassian
  bool mBoundReflMean;             /// switch for bound mean of guassian for reflection
  bool mFixedSigma;                /// fix sigma or not
  bool mFixedSigmaDoubleGaus;      /// fix sigma of 2gaussian or not
  bool mBoundSigma;                /// set bound sigma or not
  double mSigmaValue;              /// value of sigma
  double mParamSgn;                /// +/- range variation of bound Sigma of gaussian in %
  double mFracDoubleGaus;          /// initialization for fraction of 2nd gaussian in case of k2Gaus or k2GausSigmaRatioPar
  double mFixedRawYield;           /// initialization for raw yield
  bool mFixedFracDoubleGaus;       /// switch for fixed fraction of 2nd gaussian in case of k2Gaus or k2GausSigmaRatioPar
  double mRatioDoubleGausSigma;    /// initialization for ratio between two gaussian sigmas in case of k2GausSigmaRatioPar
  bool mFixedRatioDoubleGausSigma; /// switch for fixed ratio between two gaussian sigmas in case of k2GausSigmaRatioPar
  double mReflOverSgn;             /// reflection/signal
  bool mEnableReflections;         /// flag use/not use reflections
  double mRawYield;                /// signal gaussian integral
  double mRawYieldErr;             /// err on signal gaussian integral
  double mRawYieldCounted;         /// signal gaussian integral evaluated via bin counting
  double mRawYieldCountedErr;      /// err on signal gaussian integral evaluated via bin counting
  double mBkgYield;                /// background
  double mBkgYieldErr;             /// err on background
  double mSignificance;            /// significance
  double mSignificanceErr;         /// err on significance
  double mChiSquareOverNdfTotal;   /// chi2/ndf of the total fit
  double mChiSquareOverNdfBkg;     /// chi2/ndf of the background (sidebands) pre-fit
  bool mFixReflOverSgn;            /// switch for fix refl/signal
  RooRealVar* mRooMeanSgn;         /// mean for gaussian of signal
  RooRealVar* mRooSigmaSgn;        /// sigma for gaussian of signal
  RooRealVar* mRooSecSigmaSgn;     /// second sigma for composite gaussian of signal
  RooRealVar* mRooFracDoubleGaus;  /// fraction of second gaussian for composite gaussian of signal
  RooAbsPdf* mSgnPdf;              /// signal fit function
  RooAbsPdf* mBkgPdf;              /// background fit function
  RooAbsPdf* mReflPdf;             /// reflection fit function
  RooRealVar* mRooNSgn;            /// total Signal fit function integral
  RooRealVar* mRooNBkg;            /// total background fit function integral
  RooRealVar* mRooNRefl;           /// total reflection fit function integral
  RooAbsPdf* mTotalPdf;            /// total fit function
  RooPlot* mInvMassFrame;          /// frame of mass
  RooPlot* mReflFrame;             /// reflection frame
  RooPlot* mReflOnlyFrame;         /// reflection frame plot on reflection only
  RooPlot* mResidualFrame;         /// residual frame
  RooPlot* mRatioFrame;            /// fit/data ratio frame
  RooPlot* mResidualFrameForCalculation;
  RooWorkspace* mWorkspace;  /// workspace
  double mIntegralHisto;     /// integral of histogram to fit
  double mIntegralBkg;       /// integral of background fit function
  double mIntegralSgn;       /// integral of signal fit function
  TH1* mHistoTemplateRefl;   /// reflection histogram
  bool mDrawBgPrefit;        /// draw background after fitting the sidebands
  bool mHighlightPeakRegion; /// draw vertical lines showing the peak region (usually +- 3 sigma)

  ClassDefOverride(HFInvMassFitter, 1);
};

#endif // PWGHF_D2H_MACROS_HFINVMASSFITTER_H_
