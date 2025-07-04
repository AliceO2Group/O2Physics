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

/// \file HFInvMassFitter.cxx
/// \brief HFInvMassFitter class
///
/// \author Zhen Zhang <zhenz@cern.ch>
/// \author Mingyu Zhang <mingyu.zang@cern.ch>
/// \author Xinye Peng  <xinye.peng@cern.ch>
/// \author Biao Zhang <biao.zhang@cern.ch>
/// \author Oleksii Lubynets <oleksii.lubynets@cern.ch>

#include "HFInvMassFitter.h"

#include <RooAddPdf.h>
#include <RooDataHist.h>
#include <RooExponential.h>
#include <RooFitResult.h>
#include <RooFormulaVar.h>
#include <RooGamma.h>
#include <RooGaussian.h>
#include <RooGenericPdf.h>
#include <RooGlobalFunc.h>
#include <RooHist.h>
#include <RooPlot.h>
#include <RooPolynomial.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <TColor.h>
#include <TDatabasePDG.h>
#include <TLine.h>
#include <TNamed.h>
#include <TPaveText.h>
#include <TString.h>
#include <TStyle.h>
#include <TVirtualPad.h>

#include <Rtypes.h>
#include <RtypesCore.h>

#include <array>
#include <cmath>
#include <cstring>

using namespace RooFit;

ClassImp(HFInvMassFitter);

HFInvMassFitter::HFInvMassFitter() : TNamed(),
                                     mHistoInvMass(0x0),
                                     mFitOption("L,E"),
                                     mMinMass(0),
                                     mMaxMass(5),
                                     mTypeOfBkgPdf(Expo),
                                     mMassParticle(TDatabasePDG::Instance()->GetParticle("D0")->Mass()),
                                     mTypeOfSgnPdf(SingleGaus),
                                     mTypeOfReflPdf(1),
                                     mMass(1.865),
                                     mSecMass(1.969),
                                     mSigmaSgn(0.012),
                                     mSecSigma(0.006),
                                     mNSigmaForSidebands(4.),
                                     mNSigmaForSgn(3.),
                                     mSigmaSgnErr(0.),
                                     mSigmaSgnDoubleGaus(0.012),
                                     mFixedMean(kFALSE),
                                     mBoundMean(kFALSE),
                                     mBoundReflMean(kFALSE),
                                     mRooMeanSgn(0x0),
                                     mRooSigmaSgn(0x0),
                                     mMassLowLimit(0),
                                     mMassUpLimit(0),
                                     mMassReflLowLimit(0),
                                     mMassReflUpLimit(0),
                                     mFixedSigma(kFALSE),
                                     mFixedSigmaDoubleGaus(kFALSE),
                                     mBoundSigma(kFALSE),
                                     mSigmaValue(0.012),
                                     mParamSgn(0.1),
                                     mFracDoubleGaus(0.2),
                                     mFixedRawYield(-1.),
                                     mFixedFracDoubleGaus(kFALSE),
                                     mRatioDoubleGausSigma(0.),
                                     mFixedRatioDoubleGausSigma(kFALSE),
                                     mReflOverSgn(0),
                                     mEnableReflections(kFALSE),
                                     mRawYield(0),
                                     mRawYieldErr(0),
                                     mBkgYield(0),
                                     mBkgYieldErr(0),
                                     mSignificance(0),
                                     mSignificanceErr(0),
                                     mChiSquareOverNdf(0),
                                     mSgnPdf(0x0),
                                     mBkgPdf(0x0),
                                     mReflPdf(0x0),
                                     mIntegralHisto(0),
                                     mIntegralBkg(0),
                                     mIntegralSgn(0),
                                     mRooNSgn(0x0),
                                     mRooNBkg(0x0),
                                     mRooNRefl(0x0),
                                     mTotalPdf(0x0),
                                     mInvMassFrame(0x0),
                                     mReflFrame(0x0),
                                     mReflOnlyFrame(0x0),
                                     mResidualFrame(0x0),
                                     mWorkspace(0x0),
                                     mHistoTemplateRefl(0x0),
                                     mDrawBgPrefit(kFALSE),
                                     mHighlightPeakRegion(kFALSE)
{
  // default constructor
}

HFInvMassFitter::HFInvMassFitter(const TH1* histoToFit, Double_t minValue, Double_t maxValue, Int_t fitTypeBkg, Int_t fitTypeSgn) : TNamed(),
                                                                                                                                    mHistoInvMass(0x0),
                                                                                                                                    mFitOption("L,E"),
                                                                                                                                    mMinMass(minValue),
                                                                                                                                    mMaxMass(maxValue),
                                                                                                                                    mTypeOfBkgPdf(fitTypeBkg),
                                                                                                                                    mMassParticle(TDatabasePDG::Instance()->GetParticle("D0")->Mass()),
                                                                                                                                    mTypeOfSgnPdf(fitTypeSgn),
                                                                                                                                    mTypeOfReflPdf(1),
                                                                                                                                    mMass(1.865),
                                                                                                                                    mSecMass(1.969),
                                                                                                                                    mSigmaSgn(0.012),
                                                                                                                                    mSecSigma(0.006),
                                                                                                                                    mNSigmaForSidebands(3.),
                                                                                                                                    mNSigmaForSgn(3.),
                                                                                                                                    mSigmaSgnErr(0.),
                                                                                                                                    mSigmaSgnDoubleGaus(0.012),
                                                                                                                                    mFixedMean(kFALSE),
                                                                                                                                    mBoundMean(kFALSE),
                                                                                                                                    mBoundReflMean(kFALSE),
                                                                                                                                    mRooMeanSgn(0x0),
                                                                                                                                    mRooSigmaSgn(0x0),
                                                                                                                                    mMassLowLimit(0),
                                                                                                                                    mMassUpLimit(0),
                                                                                                                                    mMassReflLowLimit(0),
                                                                                                                                    mMassReflUpLimit(0),
                                                                                                                                    mFixedSigma(kFALSE),
                                                                                                                                    mFixedSigmaDoubleGaus(kFALSE),
                                                                                                                                    mBoundSigma(kFALSE),
                                                                                                                                    mSigmaValue(0.012),
                                                                                                                                    mParamSgn(0.1),
                                                                                                                                    mFracDoubleGaus(0.2),
                                                                                                                                    mFixedRawYield(-1.),
                                                                                                                                    mFixedFracDoubleGaus(kFALSE),
                                                                                                                                    mRatioDoubleGausSigma(0.),
                                                                                                                                    mFixedRatioDoubleGausSigma(kFALSE),
                                                                                                                                    mReflOverSgn(0),
                                                                                                                                    mEnableReflections(kFALSE),
                                                                                                                                    mRawYield(0),
                                                                                                                                    mRawYieldErr(0),
                                                                                                                                    mBkgYield(0),
                                                                                                                                    mBkgYieldErr(0),
                                                                                                                                    mSignificance(0),
                                                                                                                                    mSignificanceErr(0),
                                                                                                                                    mChiSquareOverNdf(0),
                                                                                                                                    mSgnPdf(0x0),
                                                                                                                                    mBkgPdf(0x0),
                                                                                                                                    mReflPdf(0x0),
                                                                                                                                    mIntegralHisto(0),
                                                                                                                                    mIntegralBkg(0),
                                                                                                                                    mIntegralSgn(0),
                                                                                                                                    mRooNSgn(0x0),
                                                                                                                                    mRooNBkg(0x0),
                                                                                                                                    mRooNRefl(0x0),
                                                                                                                                    mTotalPdf(0x0),
                                                                                                                                    mInvMassFrame(0x0),
                                                                                                                                    mReflFrame(0x0),
                                                                                                                                    mReflOnlyFrame(0x0),
                                                                                                                                    mResidualFrame(0x0),
                                                                                                                                    mWorkspace(0x0),
                                                                                                                                    mHistoTemplateRefl(0x0),
                                                                                                                                    mDrawBgPrefit(kFALSE),
                                                                                                                                    mHighlightPeakRegion(kFALSE)
{
  // standard constructor
  mHistoInvMass = dynamic_cast<TH1*>(histoToFit->Clone(histoToFit->GetTitle()));
  mHistoInvMass->SetDirectory(0);
}

HFInvMassFitter::~HFInvMassFitter()
{

  /// destructor
  delete mHistoInvMass;
  delete mHistoTemplateRefl;
  delete mRooMeanSgn;
  delete mRooSigmaSgn;
  delete mSgnPdf;
  delete mBkgPdf;
  delete mReflPdf;
  delete mTotalPdf;
  delete mRooNSgn;
  delete mRooNBkg;
  delete mRooNRefl;
  delete mInvMassFrame;
  delete mReflFrame;
  delete mReflOnlyFrame;
  delete mResidualFrame;
  delete mWorkspace;
}

void HFInvMassFitter::doFit()
{
  mIntegralHisto = mHistoInvMass->Integral(mHistoInvMass->FindBin(mMinMass), mHistoInvMass->FindBin(mMaxMass));
  mWorkspace = new RooWorkspace("mWorkspace");
  fillWorkspace(*mWorkspace);
  RooRealVar* mass = mWorkspace->var("mass");
  RooDataHist dataHistogram("dataHistogram", "data", *mass, Import(*mHistoInvMass));

  if (mTypeOfBkgPdf == NoBkg) { // MC
    mass->setRange("signal", mMass - 3. * mSigmaSgn, mMass + 3. * mSigmaSgn);
  } else {
    if (mTypeOfSgnPdf == GausSec) { // Second Peak fit range
      mass->setRange("SBL", mMinMass, mMass - mNSigmaForSidebands * mSigmaSgn);
      mass->setRange("SBR", mMass + mNSigmaForSidebands * mSigmaSgn, mSecMass - mNSigmaForSidebands * mSecSigma);
      mass->setRange("SEC", mSecMass + mNSigmaForSidebands * mSecSigma, mMaxMass);
      mass->setRange("signal", mSecMass - mNSigmaForSidebands * mSecSigma, mSecMass + mNSigmaForSidebands * mSecSigma);
    } else { // Single Peak fit range
      mass->setRange("SBL", mMinMass, mMass - mNSigmaForSidebands * mSigmaSgn);
      mass->setRange("SBR", mMass + mNSigmaForSidebands * mSigmaSgn, mMaxMass);
      mass->setRange("signal", mMass - mNSigmaForSgn * mSigmaSgn, mMass + mNSigmaForSgn * mSigmaSgn);
    }
  }
  mass->setRange("bkg", mMass - 4 * mSigmaSgn, mMass + 4 * mSigmaSgn);
  mass->setRange("full", mMinMass, mMaxMass);
  mInvMassFrame = mass->frame(Title(Form("%s", mHistoInvMass->GetTitle()))); // define the frame to plot
  dataHistogram.plotOn(mInvMassFrame, Name("data_c"));                       // plot data histogram on the frame

  // define number of background and background fit function
  mRooNBkg = new RooRealVar("mRooNBkg", "number of background", 0.3 * mIntegralHisto, 0., 1.2 * mIntegralHisto); // background yield
  RooAbsPdf* bkgPdf = createBackgroundFitFunction(mWorkspace);                                                   // Create background pdf
  RooAbsPdf* sgnPdf = createSignalFitFunction(mWorkspace);                                                       // Create signal pdf

  // fit MC or Data
  if (mTypeOfBkgPdf == NoBkg) {                                                                                // MC
    mRooNSgn = new RooRealVar("mRooNSig", "number of signal", 0.3 * mIntegralHisto, 0., 1.2 * mIntegralHisto); // signal yield
    mTotalPdf = new RooAddPdf("mMCFunc", "MC fit function", RooArgList(*sgnPdf), RooArgList(*mRooNSgn));       // create total pdf
    if (!strcmp(mFitOption.Data(), "Chi2")) {
      mTotalPdf->chi2FitTo(dataHistogram, Range("signal"));
    } else {
      mTotalPdf->fitTo(dataHistogram, Range("signal"));
    }
    RooAbsReal* signalIntegralMc = mTotalPdf->createIntegral(*mass, NormSet(*mass), Range("signal")); // sig yield from fit
    mIntegralSgn = signalIntegralMc->getValV();
    calculateSignal(mRawYield, mRawYieldErr);        // calculate signal and signal error
    mTotalPdf->plotOn(mInvMassFrame, Name("Tot_c")); // plot total function
  } else {                                           // data
    mBkgPdf = new RooAddPdf("mBkgPdf", "background fit function", RooArgList(*bkgPdf), RooArgList(*mRooNBkg));
    if (mTypeOfSgnPdf == GausSec) { // two peak fit
      if (!strcmp(mFitOption.Data(), "Chi2")) {
        mBkgPdf->chi2FitTo(dataHistogram, Range("SBL,SBR,SEC"), Save());
      } else {
        mBkgPdf->fitTo(dataHistogram, Range("SBL,SBR,SEC"), Save());
      }
    } else { // single peak fit
      if (!strcmp(mFitOption.Data(), "Chi2")) {
        mBkgPdf->chi2FitTo(dataHistogram, Range("SBL,SBR"), Save());
      } else {
        mBkgPdf->fitTo(dataHistogram, Range("SBL,SBR"), Save());
      }
    }
    RooAbsPdf* mBkgPdfPrefit{nullptr};
    if (mDrawBgPrefit) {
      mBkgPdfPrefit = dynamic_cast<RooAbsPdf*>(mBkgPdf->Clone());
      mBkgPdfPrefit->plotOn(mInvMassFrame, Range("full"), Name("Bkg_c_prefit"), LineColor(kGray));
      delete mBkgPdfPrefit;
    }

    // estimate signal yield
    RooAbsReal* bkgIntegral = mBkgPdf->createIntegral(*mass, NormSet(*mass), Range("bkg")); // bkg integral
    mIntegralBkg = bkgIntegral->getValV();                                                  // fraction of BG's integral in "bkg" range out of that in "full" range (which is 1 by construction). Not an absolute value.
    Double_t estimatedSignal;
    checkForSignal(estimatedSignal);              // SIG's absolute integral in "bkg" range
    calculateBackground(mBkgYield, mBkgYieldErr); // BG's absolute integral in "bkg" range

    mRooNSgn = new RooRealVar("mNSgn", "number of signal", 0.3 * estimatedSignal, 0., 1.2 * estimatedSignal); // estimated signal yield
    if (mFixedRawYield > 0) {
      mRooNSgn->setVal(mFixedRawYield); // fixed signal yield
      mRooNSgn->setConstant(kTRUE);
    }
    mSgnPdf = new RooAddPdf("mSgnPdf", "signal fit function", RooArgList(*sgnPdf), RooArgList(*mRooNSgn));
    // create reflection template and fit to reflection
    if (mHistoTemplateRefl) {
      RooAbsPdf* reflPdf = createReflectionFitFunction(mWorkspace); // create reflection pdf
      RooDataHist reflHistogram("reflHistogram", "refl for fit", *mass, Import(*mHistoTemplateRefl));
      mReflFrame = mass->frame();
      mReflOnlyFrame = mass->frame(Title(Form("%s", mHistoTemplateRefl->GetTitle())));
      reflHistogram.plotOn(mReflOnlyFrame);
      mRooNRefl = new RooRealVar("mNRefl", "number of reflection", 0.5 * mHistoTemplateRefl->Integral(), 0, mHistoTemplateRefl->Integral());
      RooAddPdf reflFuncTemp("reflFuncTemp", "template reflection fit function", RooArgList(*reflPdf), RooArgList(*mRooNRefl));
      if (!strcmp(mFitOption.Data(), "Chi2")) {
        reflFuncTemp.chi2FitTo(reflHistogram);
      } else {
        reflFuncTemp.fitTo(reflHistogram);
      }
      reflFuncTemp.plotOn(mReflOnlyFrame);

      mRooNRefl->setVal(mReflOverSgn * estimatedSignal);
      mRooNRefl->setConstant(kTRUE);
      setReflFuncFixed(); // fix reflection pdf parameter
      mTotalPdf = new RooAddPdf("mTotalPdf", "background + signal + reflection fit function", RooArgList(*bkgPdf, *sgnPdf, *reflPdf), RooArgList(*mRooNBkg, *mRooNSgn, *mRooNRefl));
      if (!strcmp(mFitOption.Data(), "Chi2")) {
        mTotalPdf->chi2FitTo(dataHistogram);
      } else {
        mTotalPdf->fitTo(dataHistogram);
      }
      mTotalPdf->plotOn(mInvMassFrame, Name("Tot_c"));
      mReflPdf = new RooAddPdf("mReflPdf", "reflection fit function", RooArgList(*reflPdf), RooArgList(*mRooNRefl));
      RooAddPdf reflBkgPdf("reflBkgPdf", "reflBkgPdf", RooArgList(*bkgPdf, *reflPdf), RooArgList(*mRooNBkg, *mRooNRefl));
      reflBkgPdf.plotOn(mInvMassFrame, Normalization(1.0, RooAbsReal::RelativeExpected), LineStyle(7), LineColor(kRed + 1), Name("ReflBkg_c"));
      plotBkg(mTotalPdf);                                              // plot bkg pdf in total pdf
      plotRefl(mTotalPdf);                                             // plot reflection in total pdf
      mChiSquareOverNdf = mInvMassFrame->chiSquare("Tot_c", "data_c"); // calculate reduced chi2 / NDF

      // plot residual distribution
      RooHist* residualHistogram = mInvMassFrame->residHist("data_c", "ReflBkg_c");
      mResidualFrame = mass->frame(Title("Residual Distribution"));
      mResidualFrame->addPlotable(residualHistogram, "p");
      mSgnPdf->plotOn(mResidualFrame, Normalization(1.0, RooAbsReal::RelativeExpected), LineColor(kBlue));
    } else {
      mTotalPdf = new RooAddPdf("mTotalPdf", "background + signal pdf", RooArgList(*bkgPdf, *sgnPdf), RooArgList(*mRooNBkg, *mRooNSgn));
      if (!strcmp(mFitOption.Data(), "Chi2")) {
        mTotalPdf->chi2FitTo(dataHistogram);
      } else {
        mTotalPdf->fitTo(dataHistogram);
      }
      plotBkg(mTotalPdf);
      mTotalPdf->plotOn(mInvMassFrame, Name("Tot_c"), LineColor(kBlue));
      mSgnPdf->plotOn(mInvMassFrame, Normalization(1.0, RooAbsReal::RelativeExpected), DrawOption("F"), FillColor(TColor::GetColorTransparent(kBlue, 0.2)), VLines());
      mChiSquareOverNdf = mInvMassFrame->chiSquare("Tot_c", "data_c"); // calculate reduced chi2 / DNF
      // plot residual distribution
      mResidualFrame = mass->frame(Title("Residual Distribution"));
      RooHist* residualHistogram = mInvMassFrame->residHist("data_c", "Bkg_c");
      mResidualFrame->addPlotable(residualHistogram, "P");
      mSgnPdf->plotOn(mResidualFrame, Normalization(1.0, RooAbsReal::RelativeExpected), LineColor(kBlue));
    }
    mass->setRange("bkgForSignificance", mRooMeanSgn->getVal() - mNSigmaForSgn * mRooSigmaSgn->getVal(), mRooMeanSgn->getVal() + mNSigmaForSgn * mRooSigmaSgn->getVal());
    bkgIntegral = mBkgPdf->createIntegral(*mass, NormSet(*mass), Range("bkgForSignificance"));
    mIntegralBkg = bkgIntegral->getValV();
    calculateBackground(mBkgYield, mBkgYieldErr);

    RooAbsReal* sgnIntegral = mSgnPdf->createIntegral(*mass, NormSet(*mass), Range("signal"));
    mIntegralSgn = sgnIntegral->getValV();
    calculateSignal(mRawYield, mRawYieldErr);
    countSignal(mRawYieldCounted, mRawYieldCountedErr);
    calculateSignificance(mSignificance, mSignificanceErr);
  }
}

void HFInvMassFitter::fillWorkspace(RooWorkspace& workspace) const
{
  // Declare observable variable
  RooRealVar mass("mass", "mass", mMinMass, mMaxMass, "GeV/c^{2}");
  // bkg expo
  RooRealVar tau("tau", "tau", -1, -5., 5.);
  RooAbsPdf* bkgFuncExpo = new RooExponential("bkgFuncExpo", "background fit function", mass, tau);
  workspace.import(*bkgFuncExpo);
  delete bkgFuncExpo;
  // bkg poly1
  RooRealVar polyParam0("polyParam0", "Parameter of Poly function", 0.5, -5., 5.);
  RooRealVar polyParam1("polyParam1", "Parameter of Poly function", 0.2, -5., 5.);
  RooAbsPdf* bkgFuncPoly1 = new RooPolynomial("bkgFuncPoly1", "background fit function", mass, RooArgSet(polyParam0, polyParam1));
  workspace.import(*bkgFuncPoly1);
  delete bkgFuncPoly1;
  // bkg poly2
  RooRealVar polyParam2("polyParam2", "Parameter of Poly function", 0.2, -5., 5.);
  RooAbsPdf* bkgFuncPoly2 = new RooPolynomial("bkgFuncPoly2", "background fit function", mass, RooArgSet(polyParam0, polyParam1, polyParam2));
  workspace.import(*bkgFuncPoly2);
  delete bkgFuncPoly2;
  // bkg poly3
  RooRealVar polyParam3("polyParam3", "Parameter of Poly function", 0.2, -1., 1.);
  RooAbsPdf* bkgFuncPoly3 = new RooPolynomial("bkgFuncPoly3", "background pdf", mass, RooArgSet(polyParam0, polyParam1, polyParam2, polyParam3));
  workspace.import(*bkgFuncPoly3);
  delete bkgFuncPoly3;
  // bkg power law
  RooRealVar powParam1("powParam1", "Parameter of Pow function", TDatabasePDG::Instance()->GetParticle("pi+")->Mass());
  RooRealVar powParam2("powParam2", "Parameter of Pow function", 1., -10, 10);
  RooAbsPdf* bkgFuncPow = new RooGenericPdf("bkgFuncPow", "bkgFuncPow", "(mass-powParam1)^powParam2", RooArgSet(mass, powParam1, powParam2));
  workspace.import(*bkgFuncPow);
  delete bkgFuncPow;
  // pow * exp
  RooRealVar powExpoParam1("powExpoParam1", "Parameter of PowExpo function", 1 / 2);
  RooRealVar powExpoParam2("powExpoParam2", "Parameter of PowExpo function", 1, -10, 10);
  RooRealVar massPi("massPi", "mass of pion", TDatabasePDG::Instance()->GetParticle("pi+")->Mass());
  RooFormulaVar powExpoParam3("powExpoParam3", "powExpoParam1 + 1", RooArgList(powExpoParam1));
  RooFormulaVar powExpoParam4("powExpoParam4", "1./powExpoParam2", RooArgList(powExpoParam2));
  RooAbsPdf* bkgFuncPowExpo = new RooGamma("bkgFuncPowExpo", "background pdf", mass, powExpoParam3, powExpoParam4, massPi);
  workspace.import(*bkgFuncPowExpo);
  delete bkgFuncPowExpo;

  // signal pdf
  RooRealVar mean("mean", "mean for signal fit", mMass, 0, 5);
  if (mBoundMean) {
    mean.setMax(mMassUpLimit);
    mean.setMin(mMassLowLimit);
  }
  // signal Gaussian
  if (mFixedMean) {
    mean.setVal(mMass);
    mean.setConstant(kTRUE);
  }
  RooRealVar sigma("sigma", "sigma for signal", mSigmaSgn, mSigmaSgn - 0.01, mSigmaSgn + 0.01);
  if (mFixedSigma) {
    sigma.setVal(mSigmaSgn);
    sigma.setConstant(kTRUE);
  }
  if (mBoundSigma) {
    sigma.setMax(mSigmaSgn * (1 + mParamSgn));
    sigma.setMin(mSigmaSgn * (1 - mParamSgn));
  }
  RooAbsPdf* sgnFuncGaus = new RooGaussian("sgnFuncGaus", "signal pdf", mass, mean, sigma);
  workspace.import(*sgnFuncGaus);
  delete sgnFuncGaus;
  // signal double Gaussian
  RooRealVar sigmaDoubleGaus("sigmaDoubleGaus", "sigma2Gaus", mSigmaSgn, mSigmaSgn - 0.01, mSigmaSgn + 0.01);
  if (mBoundSigma) {
    sigmaDoubleGaus.setMax(mSigmaSgn * (1 + mParamSgn));
    sigmaDoubleGaus.setMin(mSigmaSgn * (1 - mParamSgn));
  }
  if (mFixedSigma) {
    sigma.setVal(mSigmaSgn);
    sigma.setConstant(kTRUE);
  }
  if (mFixedSigmaDoubleGaus) {
    sigmaDoubleGaus.setVal(mSigmaSgnDoubleGaus);
    sigmaDoubleGaus.setConstant(kTRUE);
  }
  RooGaussian gaus1("gaus1", "gaus1", mass, mean, sigma);
  RooGaussian gaus2("gaus2", "gaus2", mass, mean, sigmaDoubleGaus);
  RooRealVar fracDoubleGaus("fracDoubleGaus", "frac of two gauss", mFracDoubleGaus, 0, 1.);
  if (mFixedFracDoubleGaus) {
    fracDoubleGaus.setVal(mFracDoubleGaus);
    fracDoubleGaus.setConstant(kTRUE);
  }
  RooAbsPdf* sgnFuncDoubleGaus = new RooAddPdf("sgnFuncDoubleGaus", "signal pdf", RooArgList(gaus1, gaus2), fracDoubleGaus);
  workspace.import(*sgnFuncDoubleGaus);
  delete sgnFuncDoubleGaus;
  // double Gaussian ratio
  RooRealVar ratio("ratio", "ratio of sigma12", mRatioDoubleGausSigma, 0, 10);
  if (mFixedSigma) {
    sigma.setVal(mSigmaSgn);
    sigma.setConstant(kTRUE);
  }
  if (mFixedRatioDoubleGausSigma) {
    ratio.setVal(mRatioDoubleGausSigma);
    ratio.setConstant(kTRUE);
  }
  if (mBoundSigma) {
    sigma.setMax(mSigmaSgn * (1 + mParamSgn));
    sigma.setMin(mSigmaSgn * (1 - mParamSgn));
  }
  RooRealVar sigmaDoubleGausRatio("sigmaDoubleGausRatio", "sigmaDoubleGausRatio", sigma.getVal() * ratio.getVal());
  RooGaussian gausRatio1("gausRatio1", "gausratio1", mass, mean, sigma);
  RooGaussian gausRatio2("gausRatio2", "gausratio2", mass, mean, sigmaDoubleGausRatio);
  RooRealVar fracDoubleGausRatio("fracDoubleGausRatio", "fraction of two gauss ratio", 0.5, 0, 1.);
  if (mFixedFracDoubleGaus) {
    fracDoubleGausRatio.setVal(mFracDoubleGaus);
    fracDoubleGausRatio.setConstant(kTRUE);
  }
  RooAbsPdf* sgnFuncGausRatio = new RooAddPdf("sgnFuncGausRatio", "signal pdf", RooArgList(gausRatio1, gausRatio2), fracDoubleGausRatio);
  workspace.import(*sgnFuncGausRatio);
  delete sgnFuncGausRatio;
  // double peak for Ds
  RooRealVar meanSec("meanSec", "mean for second peak fit", mSecMass, mMinMass, mMaxMass);
  RooRealVar sigmaSec("sigmaSec", "sigmaSec", mSecSigma, mSecSigma - 0.005, mSecSigma + 0.01);
  if (mFixedMean) {
    meanSec.setVal(mSecMass);
    meanSec.setConstant(kTRUE);
  }
  if (mBoundMean) {
    meanSec.setMax(mMassUpLimit);
    meanSec.setMin(mMassLowLimit);
  }
  if (mFixedSigma) {
    sigmaSec.setVal(mSecSigma);
    sigmaSec.setConstant(kTRUE);
  }
  if (mBoundSigma) {
    sigmaSec.setMax(mSecSigma * (1 + mParamSgn));
    sigmaSec.setMin(mSecSigma * (1 - mParamSgn));
  }
  RooGaussian gausSec1("gausSec1", "gausSec1", mass, mean, sigmaSec);
  RooGaussian gausSec2("gausSec2", "gausSec2", mass, meanSec, sigmaSec);
  RooRealVar fracSec("fracSec", "frac of two peak", 0.5, 0, 1.);
  RooAbsPdf* sgnFuncDoublePeak = new RooAddPdf("sgnFuncDoublePeak", "signal pdf", RooArgList(gausSec1, gausSec2), fracSec);
  workspace.import(*sgnFuncDoublePeak);
  delete sgnFuncDoublePeak;
  // reflection Gaussian
  RooRealVar meanRefl("meanRefl", "mean for reflection", mMass, 0.0, mMass + 0.05);
  if (mBoundReflMean) {
    meanRefl.setMax(mMassReflUpLimit);
    meanRefl.setMin(mMassReflLowLimit);
  }
  RooRealVar sigmaRefl("sigmaRefl", "sigma for reflection", 0.012, 0, 0.25);
  RooAbsPdf* reflFuncGaus = new RooGaussian("reflFuncGaus", "reflection pdf", mass, meanRefl, sigmaRefl);
  workspace.import(*reflFuncGaus);
  delete reflFuncGaus;
  // reflection double Gaussian
  RooRealVar meanReflDoubleGaus("meanReflDoubleGaus", "mean for reflection double gaussian", mMass, 0.0, mMass + 0.05);
  if (mBoundReflMean) {
    meanReflDoubleGaus.setMax(mMassReflUpLimit);
    meanReflDoubleGaus.setMin(mMassReflLowLimit);
  }
  RooRealVar sigmaReflDoubleGaus("sigmaReflDoubleGaus", "sigmaReflDoubleGaus", 0.012, 0.0, 0.25);
  RooGaussian gausRefl1("gausRefl1", "gausRefl1", mass, meanRefl, sigmaRefl);
  RooGaussian gausRefl2("gausRefl2", "gausRefl2", mass, meanReflDoubleGaus, sigmaReflDoubleGaus);
  RooRealVar fracRefl("fracRefl", "frac of two gauss", 0.5, 0, 1.);
  RooAbsPdf* reflFuncDoubleGaus = new RooAddPdf("reflFuncDoubleGaus", "reflection pdf", RooArgList(gausRefl1, gausRefl2), fracRefl);
  workspace.import(*reflFuncDoubleGaus);
  delete reflFuncDoubleGaus;
  // reflection poly3
  RooRealVar polyReflParam0("polyReflParam0", "polyReflParam0", 0.5, -1., 1.);
  RooRealVar polyReflParam1("polyReflParam1", "polyReflParam1", 0.2, -1., 1.);
  RooRealVar polyReflParam2("polyReflParam2", "polyReflParam2", 0.2, -1., 1.);
  RooRealVar polyReflParam3("polyReflParam3", "polyReflParam3", 0.2, -1., 1.);
  RooAbsPdf* reflFuncPoly3 = new RooPolynomial("reflFuncPoly3", "reflection PDF", mass, RooArgSet(polyReflParam0, polyReflParam1, polyReflParam2, polyReflParam3));
  workspace.import(*reflFuncPoly3);
  delete reflFuncPoly3;
  // reflection poly6
  RooRealVar polyReflParam4("polyReflParam4", "polyReflParam4", 0.2, -1., 1.);
  RooRealVar polyReflParam5("polyReflParam5", "polyReflParam5", 0.2, -1., 1.);
  RooRealVar polyReflParam6("polyReflParam6", "polyReflParam6", 0.2, -1., 1.);
  RooAbsPdf* reflFuncPoly6 = new RooPolynomial("reflFuncPoly6", "reflection pdf", mass, RooArgSet(polyReflParam0, polyReflParam1, polyReflParam2, polyReflParam3, polyReflParam4, polyReflParam5, polyReflParam6));
  workspace.import(*reflFuncPoly6);
  delete reflFuncPoly6;
}
// draw fit output
void HFInvMassFitter::drawFit(TVirtualPad* pad, Int_t writeFitInfo)
{
  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  pad->cd();
  if (writeFitInfo > 0) {
    TPaveText* textInfoLeft = new TPaveText(0.12, 0.65, 0.47, 0.89, "NDC");
    TPaveText* textInfoRight = new TPaveText(0.6, 0.7, 1., .87, "NDC");
    textInfoLeft->SetBorderSize(0);
    textInfoLeft->SetFillStyle(0);
    textInfoRight->SetBorderSize(0);
    textInfoRight->SetFillStyle(0);
    textInfoRight->SetTextColor(kBlue);
    textInfoLeft->AddText(Form("S = %.0f #pm %.0f ", mRawYield, mRawYieldErr));
    textInfoLeft->AddText(Form("S_{count} = %.0f #pm %.0f ", mRawYieldCounted, mRawYieldCountedErr));
    if (mTypeOfBkgPdf != NoBkg) {
      textInfoLeft->AddText(Form("B (%d#sigma) = %.0f #pm %.0f", mNSigmaForSidebands, mBkgYield, mBkgYieldErr));
      textInfoLeft->AddText(Form("S/B (%d#sigma) = %.4g ", mNSigmaForSidebands, mRawYield / mBkgYield));
    }
    if (mReflPdf) {
      textInfoLeft->AddText(Form("Refl/Sig =  %.3f #pm %.3f ", mReflOverSgn, 0.0));
    }
    if (mTypeOfBkgPdf != NoBkg) {
      textInfoLeft->AddText(Form("Signif (%d#sigma) = %.1f #pm %.1f ", mNSigmaForSidebands, mSignificance, mSignificanceErr));
      textInfoLeft->AddText(Form("#chi^{2} / ndf  =  %.3f", mChiSquareOverNdf));
    }
    if (mFixedMean) {
      textInfoRight->AddText(Form("mean(fixed) = %.3f #pm %.3f", mRooMeanSgn->getVal(), mRooMeanSgn->getError()));
    } else {
      textInfoRight->AddText(Form("mean(free) = %.3f #pm %.3f", mRooMeanSgn->getVal(), mRooMeanSgn->getError()));
    }
    if (mFixedSigma) {
      textInfoRight->AddText(Form("sigma(fixed) = %.3f #pm %.3f", mRooSigmaSgn->getVal(), mRooSigmaSgn->getError()));
    } else {
      textInfoRight->AddText(Form("sigma(free) = %.3f #pm %.3f", mRooSigmaSgn->getVal(), mRooSigmaSgn->getError()));
    }
    mInvMassFrame->addObject(textInfoLeft);
    mInvMassFrame->addObject(textInfoRight);
    mInvMassFrame->GetYaxis()->SetTitleOffset(1.8);
    gPad->SetLeftMargin(0.15);
    mInvMassFrame->GetYaxis()->SetTitle(Form("%s", mHistoInvMass->GetYaxis()->GetTitle()));
    mInvMassFrame->GetXaxis()->SetTitle(Form("%s", mHistoInvMass->GetXaxis()->GetTitle()));
    mInvMassFrame->Draw();
    highlightPeakRegion(mInvMassFrame);
    if (mHistoTemplateRefl) {
      mReflFrame->Draw("same");
    }
  }
}

// draw residual distribution on canvas
void HFInvMassFitter::drawResidual(TVirtualPad* pad)
{
  pad->cd();
  mResidualFrame->GetYaxis()->SetTitle("");
  TPaveText* textInfo = new TPaveText(0.12, 0.65, 0.47, .89, "NDC");
  textInfo->SetBorderSize(0);
  textInfo->SetFillStyle(0);
  textInfo->SetTextColor(kBlue);
  textInfo->AddText(Form("S = %.0f #pm %.0f ", mRawYield, mRawYieldErr));
  textInfo->AddText(Form("S_{count} = %.0f #pm %.0f ", mRawYieldCounted, mRawYieldCountedErr));
  textInfo->AddText(Form("mean = %.3f #pm %.3f", mRooMeanSgn->getVal(), mRooMeanSgn->getError()));
  textInfo->AddText(Form("sigma = %.3f #pm %.3f", mRooSigmaSgn->getVal(), mRooSigmaSgn->getError()));
  mResidualFrame->addObject(textInfo);
  mResidualFrame->Draw();
  highlightPeakRegion(mResidualFrame);
}

// draw peak region with vertical lines
void HFInvMassFitter::highlightPeakRegion(const RooPlot* plot, Color_t color, Width_t width, Style_t style) const
{
  if (!mHighlightPeakRegion)
    return;
  double yMin = plot->GetMinimum();
  double yMax = plot->GetMaximum();
  const Double_t mean = mRooMeanSgn->getVal();
  const Double_t sigma = mRooSigmaSgn->getVal();
  const Double_t minForSgn = mean - mNSigmaForSidebands * sigma;
  const Double_t maxForSgn = mean + mNSigmaForSidebands * sigma;
  TLine* leftLine = new TLine(minForSgn, yMin, minForSgn, yMax);
  TLine* rightLine = new TLine(maxForSgn, yMin, maxForSgn, yMax);
  for (const auto& line : std::array<TLine*, 2>{leftLine, rightLine}) {
    line->SetLineColor(color);
    line->SetLineWidth(width);
    line->SetLineStyle(style);
    line->Draw();
  }
}

// draw reflection distribution on canvas
void HFInvMassFitter::drawReflection(TVirtualPad* pad)
{
  pad->cd();
  mReflOnlyFrame->GetYaxis()->SetTitle("");
  mReflOnlyFrame->Draw();
}

// calculate signal yield via bin counting
void HFInvMassFitter::countSignal(Double_t& signal, Double_t& signalErr) const
{
  const Double_t mean = mRooMeanSgn->getVal();
  const Double_t sigma = mRooSigmaSgn->getVal();
  const Double_t minForSgn = mean - mNSigmaForSidebands * sigma;
  const Double_t maxForSgn = mean + mNSigmaForSidebands * sigma;
  const Int_t binForMinSgn = mHistoInvMass->FindBin(minForSgn);
  const Int_t binForMaxSgn = mHistoInvMass->FindBin(maxForSgn);
  const Double_t binForMinSgnUpperEdge = mHistoInvMass->GetBinLowEdge(binForMinSgn + 1);
  const Double_t binForMaxSgnLowerEdge = mHistoInvMass->GetBinLowEdge(binForMaxSgn);
  const Double_t binForMinSgnFraction = (binForMinSgnUpperEdge - minForSgn) / mHistoInvMass->GetBinWidth(binForMinSgn);
  const Double_t binForMaxSgnFraction = (maxForSgn - binForMaxSgnLowerEdge) / mHistoInvMass->GetBinWidth(binForMaxSgn);

  Double_t sum = 0;
  sum += mHistoInvMass->GetBinContent(binForMinSgn) * binForMinSgnFraction;
  for (Int_t iBin = binForMinSgn + 1; iBin <= binForMaxSgn - 1; iBin++) {
    sum += mHistoInvMass->GetBinContent(iBin);
  }
  sum += mHistoInvMass->GetBinContent(binForMaxSgn) * binForMaxSgnFraction;

  Double_t bkg, errBkg;
  calculateBackground(bkg, errBkg);

  signal = sum - bkg;
  signalErr = std::sqrt(sum + errBkg * errBkg); // sum error squared is equal to sum
}

// calculate signal yield
void HFInvMassFitter::calculateSignal(Double_t& signal, Double_t& errSignal) const
{
  signal = mRooNSgn->getVal();
  errSignal = mRooNSgn->getError();
}

// calculate background yield
void HFInvMassFitter::calculateBackground(Double_t& bkg, Double_t& errBkg) const
{
  bkg = mRooNBkg->getVal() * mIntegralBkg;
  errBkg = mRooNBkg->getError() * mIntegralBkg;
}

// calculate significance
void HFInvMassFitter::calculateSignificance(Double_t& significance, Double_t& errSignificance) const
{
  Double_t signal, errSignal;
  calculateSignal(signal, errSignal);
  Double_t bkg, errBkg;
  calculateBackground(bkg, errBkg);
  Double_t sgnErrSquare = errSignal * errSignal;
  Double_t bkgErrSquare = errBkg * errBkg;
  Double_t totalSgnBkg = signal + bkg;
  significance = signal / std::sqrt(signal + bkg);
  errSignificance = significance * std::sqrt((sgnErrSquare + bkgErrSquare) / (mNSigmaForSidebands * totalSgnBkg * totalSgnBkg) + (bkg / totalSgnBkg) * (sgnErrSquare / signal / signal));
}

// estimate Signal
void HFInvMassFitter::checkForSignal(Double_t& estimatedSignal)
{
  Double_t minForSgn = mMass - 4 * mSigmaSgn;
  Double_t maxForSgn = mMass + 4 * mSigmaSgn;
  Int_t binForMinSgn = mHistoInvMass->FindBin(minForSgn);
  Int_t binForMaxSgn = mHistoInvMass->FindBin(maxForSgn);

  Double_t sum = 0;
  for (Int_t i = binForMinSgn; i <= binForMaxSgn; i++) {
    sum += mHistoInvMass->GetBinContent(i);
  }
  Double_t bkg, errBkg;
  calculateBackground(bkg, errBkg);
  estimatedSignal = sum - bkg;
}

// Create Background Fit Function
RooAbsPdf* HFInvMassFitter::createBackgroundFitFunction(RooWorkspace* workspace) const
{
  RooAbsPdf* bkgPdf;
  switch (mTypeOfBkgPdf) {
    case 0: // exponential
    {
      bkgPdf = workspace->pdf("bkgFuncExpo");
    } break;
    case 1: // poly1
    {
      bkgPdf = workspace->pdf("bkgFuncPoly1");
    } break;
    case 2: {
      bkgPdf = workspace->pdf("bkgFuncPoly2");
    } break;
    case 3: {
      bkgPdf = workspace->pdf("bkgFuncPow");
    } break;
    case 4: {
      bkgPdf = workspace->pdf("bkgFuncPowExpo");
    } break;
    case 5: {
      bkgPdf = workspace->pdf("bkgFuncPoly3");
    } break;
    case 6: // MC
      break;
    default:
      break;
  }
  return bkgPdf;
}

// Create Signal Fit Function
RooAbsPdf* HFInvMassFitter::createSignalFitFunction(RooWorkspace* workspace)
{
  RooAbsPdf* sgnPdf;
  switch (mTypeOfSgnPdf) {
    case 0: {
      sgnPdf = workspace->pdf("sgnFuncGaus");
      mRooSigmaSgn = workspace->var("sigma");
      mRooMeanSgn = workspace->var("mean");
    } break;
    case 1: {
      sgnPdf = workspace->pdf("sgnFuncDoubleGaus");
      mRooSigmaSgn = workspace->var("sigmaDoubleGaus");
      mRooMeanSgn = workspace->var("mean");
    } break;
    case 2: {
      sgnPdf = workspace->pdf("sgnFuncGausRatio");
      mRooSigmaSgn = workspace->var("sigmaDoubleGausRatio");
      mRooMeanSgn = workspace->var("mean");
    } break;
    case 3: {
      sgnPdf = workspace->pdf("sgnFuncDoublePeak");
      mRooSigmaSgn = workspace->var("sigmaSec");
      mRooMeanSgn = workspace->var("meanSec");
    } break;
    default:
      break;
  }
  return sgnPdf;
}

// Create Reflection Fit Function
RooAbsPdf* HFInvMassFitter::createReflectionFitFunction(RooWorkspace* workspace) const
{
  RooAbsPdf* reflPdf;
  switch (mTypeOfReflPdf) {
    case 0: {
      reflPdf = workspace->pdf("reflFuncGaus");
    } break;
    case 1: {
      reflPdf = workspace->pdf("reflFuncDoubleGaus");
    } break;
    case 2: {
      reflPdf = workspace->pdf("reflFuncPoly3");
    } break;
    case 3: {
      reflPdf = workspace->pdf("reflFuncPoly6");
    } break;
    default:
      break;
  }
  return reflPdf;
}

// Plot Bkg components of fTotFunction
void HFInvMassFitter::plotBkg(RooAbsPdf* pdf, Color_t color)
{
  switch (mTypeOfBkgPdf) {
    case 0:
      pdf->plotOn(mInvMassFrame, Components("bkgFuncExpo"), Name("Bkg_c"), LineColor(color));
      break;
    case 1:
      pdf->plotOn(mInvMassFrame, Components("bkgFuncPoly1"), Name("Bkg_c"), LineColor(color));
      break;
    case 2:
      pdf->plotOn(mInvMassFrame, Components("bkgFuncPoly2"), Name("Bkg_c"), LineColor(color));
      break;
    case 3:
      pdf->plotOn(mInvMassFrame, Components("bkgFuncPow"), Name("Bkg_c"), LineColor(color));
      break;
    case 4:
      pdf->plotOn(mInvMassFrame, Components("bkgFuncPowExp"), Name("Bkg_c"), LineColor(color));
      break;
    case 5:
      pdf->plotOn(mInvMassFrame, Components("bkgFuncPoly3"), Name("Bkg_c"), LineColor(color));
      break;
    case 6:
      break;
    default:
      break;
  }
}

// Plot Refl distribution on canvas
void HFInvMassFitter::plotRefl(RooAbsPdf* pdf)
{
  switch (mTypeOfReflPdf) {
    case 0:
      pdf->plotOn(mInvMassFrame, Components("reflFuncGaus"), Name("Refl_c"), LineColor(kGreen));
      break;
    case 1:
      pdf->plotOn(mInvMassFrame, Components("reflFuncDoubleGaus"), Name("Refl_c"), LineColor(kGreen));
      break;
    case 2:
      pdf->plotOn(mInvMassFrame, Components("reflFuncPoly3"), Name("Refl_c"), LineColor(kGreen));
      break;
    case 3:
      pdf->plotOn(mInvMassFrame, Components("reflFuncPoly6"), Name("Refl_c"), LineColor(kGreen));
      break;
    default:
      break;
  }
}

// Fix reflection pdf
void HFInvMassFitter::setReflFuncFixed()
{
  switch (mTypeOfReflPdf) {
    case 0: // exponential
    {
      RooRealVar* meanRefl = mWorkspace->var("meanRefl");
      RooRealVar* sigmaRefl = mWorkspace->var("sigmaRefl");
      meanRefl->setConstant(kTRUE);
      sigmaRefl->setConstant(kTRUE);
    } break;
    case 1: // poly1
    {
      RooRealVar* meanRefl = mWorkspace->var("meanRefl");
      RooRealVar* sigmaRefl = mWorkspace->var("sigmaRefl");
      RooRealVar* meanReflDoubleGaus = mWorkspace->var("meanReflDoubleGaus");
      RooRealVar* sigmaReflDoubleGaus = mWorkspace->var("sigmaReflDoubleGaus");
      RooRealVar* fracRefl = mWorkspace->var("fracRefl");
      meanRefl->setConstant(kTRUE);
      sigmaRefl->setConstant(kTRUE);
      meanReflDoubleGaus->setConstant(kTRUE);
      sigmaReflDoubleGaus->setConstant(kTRUE);
      fracRefl->setConstant(kTRUE);
    } break;
    case 2: {
      RooRealVar* polyReflParam0 = mWorkspace->var("polyReflParam0");
      RooRealVar* polyReflParam1 = mWorkspace->var("polyReflParam1");
      RooRealVar* polyReflParam2 = mWorkspace->var("polyReflParam2");
      RooRealVar* polyReflParam3 = mWorkspace->var("polyReflParam3");
      polyReflParam0->setConstant(kTRUE);
      polyReflParam1->setConstant(kTRUE);
      polyReflParam2->setConstant(kTRUE);
      polyReflParam3->setConstant(kTRUE);
    } break;
    case 3: {
      RooRealVar* polyReflParam0 = mWorkspace->var("polyReflParam0");
      RooRealVar* polyReflParam1 = mWorkspace->var("polyReflParam1");
      RooRealVar* polyReflParam2 = mWorkspace->var("polyReflParam2");
      RooRealVar* polyReflParam3 = mWorkspace->var("polyReflParam3");
      RooRealVar* polyReflParam4 = mWorkspace->var("polyReflParam4");
      RooRealVar* polyReflParam5 = mWorkspace->var("polyReflParam5");
      RooRealVar* polyReflParam6 = mWorkspace->var("polyReflParam6");
      polyReflParam0->setConstant(kTRUE);
      polyReflParam1->setConstant(kTRUE);
      polyReflParam2->setConstant(kTRUE);
      polyReflParam3->setConstant(kTRUE);
      polyReflParam4->setConstant(kTRUE);
      polyReflParam5->setConstant(kTRUE);
      polyReflParam6->setConstant(kTRUE);
    } break;
    default:
      break;
  }
}
