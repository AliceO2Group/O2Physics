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

/// \file runMassFitter.C
/// \brief HFInvMassFitter class steering macro
///
/// \author Zhen Zhang <zhenz@cern.ch>
/// \author Mingyu Zhang <mingyu.zang@cern.ch>
/// \author Xinye Peng  <xinye.peng@cern.ch>
/// \author Biao Zhang <biao.zhang@cern.ch>
/// \author Oleksii Lubynets <oleksii.lubynets@cern.ch>
/// \author Phil Stahlhut <phil.lennart.stahlhut@cern.ch>

#if !defined(__CINT__) || defined(__CLING__)

#include "HFInvMassFitter.h"

#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TVirtualPad.h>

#include <rapidjson/document.h> // if .h file not found, please include your local rapidjson/document.h and rapidjson/filereadstream.h here
#include <rapidjson/filereadstream.h>

#include <Rtypes.h>
#include <RtypesCore.h>

#include <cmath>
#include <cstdio> // for fclose
#include <functional>
#include <map>
#include <stdexcept>
#include <string> // std::string
#include <tuple>
#include <vector> // std::vector

#endif

using namespace rapidjson;

int runMassFitter(const TString& configFileName = "config_massfitter.json");

template <typename ValueType>
void readArray(const Value& jsonArray, std::vector<ValueType>& output)
{
  for (const auto* it = jsonArray.Begin(); it != jsonArray.End(); it++) {
    auto value = it->template Get<ValueType>();
    output.emplace_back(value);
  }
}

void parseStringArray(const Value& jsonArray, std::vector<std::string>& output)
{
  size_t const arrayLength = jsonArray.Size();
  for (size_t i = 0; i < arrayLength; i++) {
    if (jsonArray[i].IsString()) {
      output.emplace_back(jsonArray[i].GetString());
    }
  }
}

void divideCanvas(TCanvas* c, int nSliceVarBins);
void setHistoStyle(TH1* histo, Color_t color = kBlack, Size_t markerSize = 1);

int runMassFitter(const TString& configFileName)
{
  // load config
  FILE* configFile = fopen(configFileName.Data(), "r");
  if (configFile == nullptr) {
    throw std::runtime_error("ERROR: Missing configuration json file: " + configFileName);
  }

  Document config;
  char readBuffer[65536];
  FileReadStream is(configFile, readBuffer, sizeof(readBuffer));
  config.ParseStream(is);
  fclose(configFile);

  Bool_t const isMc = config["IsMC"].GetBool();
  Bool_t const writeSignalPar = config["WriteSignalPar"].GetBool();
  TString const particleName = config["Particle"].GetString();
  TString const collisionSystem = config["CollisionSystem"].GetString();
  TString const inputFileName = config["InFileName"].GetString();
  TString const reflFileName = config["ReflFileName"].GetString();
  TString outputFileName = config["OutFileName"].GetString();

  std::vector<std::string> inputHistoName;
  std::vector<std::string> promptHistoName;
  std::vector<std::string> fdHistoName;
  std::vector<std::string> reflHistoName;
  std::vector<std::string> promptSecPeakHistoName;
  std::vector<std::string> fdSecPeakHistoName;
  TString sliceVarName;
  TString sliceVarUnit;
  std::vector<double> sliceVarMin;
  std::vector<double> sliceVarMax;
  std::vector<double> massMin;
  std::vector<double> massMax;
  std::vector<double> fixMeanManual;
  std::vector<double> fixSigmaManual;
  std::vector<double> fixSecondSigmaManual;
  std::vector<double> fixFracDoubleGausManual;
  std::vector<int> nRebin;
  std::vector<int> bkgFuncConfig;
  std::vector<int> sgnFuncConfig;

  const Value& inputHistoNameValue = config["InputHistoName"];
  parseStringArray(inputHistoNameValue, inputHistoName);

  const Value& promptHistoNameValue = config["PromptHistoName"];
  parseStringArray(promptHistoNameValue, promptHistoName);

  const Value& fdHistoNameValue = config["FDHistoName"];
  parseStringArray(fdHistoNameValue, fdHistoName);

  const Value& reflHistoNameValue = config["ReflHistoName"];
  parseStringArray(reflHistoNameValue, reflHistoName);

  const Value& promptSecPeakHistoNameValue = config["PromptSecPeakHistoName"];
  parseStringArray(promptSecPeakHistoNameValue, promptSecPeakHistoName);

  const Value& fdSecPeakHistoNameValue = config["FDSecPeakHistoName"];
  parseStringArray(fdSecPeakHistoNameValue, fdSecPeakHistoName);

  const bool fixMean = config["FixMean"].GetBool();
  const std::string meanFile = config["MeanFile"].GetString();

  const Value& fixMeanManualValue = config["FixMeanManual"];
  readArray(fixMeanManualValue, fixMeanManual);

  const bool fixSigma = config["FixSigma"].GetBool();
  const std::string sigmaFile = config["SigmaFile"].GetString();

  const Value& fixSigmaManualValue = config["FixSigmaManual"];
  readArray(fixSigmaManualValue, fixSigmaManual);

  const bool fixSecondSigma = config["FixSecondSigma"].GetBool();
  const std::string secondSigmaFile = config["SecondSigmaFile"].GetString();

  const Value& fixSecondSigmaManualValue = config["FixSecondSigmaManual"];
  readArray(fixSecondSigmaManualValue, fixSecondSigmaManual);

  const bool fixFracDoubleGaus = config["FixFracDoubleGaus"].GetBool();
  const std::string fracDoubleGausFile = config["FracDoubleGausFile"].GetString();

  const Value& fixFracDoubleGausManualValue = config["FixFracDoubleGausManual"];
  readArray(fixFracDoubleGausManualValue, fixFracDoubleGausManual);

  sliceVarName = config["SliceVarName"].GetString();
  sliceVarUnit = config["SliceVarUnit"].GetString();

  const Value& sliceVarMinValue = config["SliceVarMin"];
  readArray(sliceVarMinValue, sliceVarMin);

  const Value& sliceVarMaxValue = config["SliceVarMax"];
  readArray(sliceVarMaxValue, sliceVarMax);

  const Value& massMinValue = config["MassMin"];
  readArray(massMinValue, massMin);

  const Value& massMaxValue = config["MassMax"];
  readArray(massMaxValue, massMax);

  const Value& rebinValue = config["Rebin"];
  readArray(rebinValue, nRebin);

  bool const includeSecPeak = config["InclSecPeak"].GetBool();
  bool const useLikelihood = config["UseLikelihood"].GetBool();

  const Value& bkgFuncValue = config["BkgFunc"];
  readArray(bkgFuncValue, bkgFuncConfig);

  const Value& sgnFuncValue = config["SgnFunc"];
  readArray(sgnFuncValue, sgnFuncConfig);

  const bool enableRefl = config["EnableRefl"].GetBool();

  const bool drawBgPrefit = config["drawBgPrefit"].GetBool();
  const bool highlightPeakRegion = config["highlightPeakRegion"].GetBool();

  const unsigned int nSliceVarBins = sliceVarMin.size();
  std::vector<int> bkgFunc(nSliceVarBins);
  std::vector<int> sgnFunc(nSliceVarBins);
  std::vector<double> sliceVarLimits(nSliceVarBins + 1);

  for (unsigned int iSliceVar = 0; iSliceVar < nSliceVarBins; iSliceVar++) {
    sliceVarLimits[iSliceVar] = sliceVarMin[iSliceVar];
    sliceVarLimits[iSliceVar + 1] = sliceVarMax[iSliceVar];

    if (bkgFuncConfig[iSliceVar] < 0 || bkgFuncConfig[iSliceVar] >= HFInvMassFitter::NTypesOfBkgPdf) {
      throw std::runtime_error("ERROR: only Expo, Poly1, Poly2, Pow and PowEx background functions supported! Exit");
    }
    bkgFunc[iSliceVar] = bkgFuncConfig[iSliceVar];

    if (sgnFuncConfig[iSliceVar] < 0 || sgnFuncConfig[iSliceVar] >= HFInvMassFitter::NTypesOfSgnPdf) {
      throw std::runtime_error("ERROR: only SingleGaus, DoubleGaus and DoubleGausSigmaRatioPar signal functions supported! Exit");
    }
    sgnFunc[iSliceVar] = sgnFuncConfig[iSliceVar];
  }

  std::map<std::string, std::tuple<std::string, std::string, std::string>> particles{
    {"Dplus", {"K#pi#pi", "D+", "D^{+} #rightarrow K^{-}#pi^{+}#pi^{+} + c.c."}},
    {"D0", {"K#pi", "D0", "D^{0} #rightarrow K^{-}#pi^{+} + c.c."}},
    {"Ds", {"KK#pi", "D_s+", "D_{s}^{+} #rightarrow K^{-}K^{+}#pi^{+} + c.c."}},
    {"LcToPKPi", {"pK#pi", "Lambda_c+", "#Lambda_{c}^{+} #rightarrow pK^{-}#pi^{+} + c.c."}},
    {"LcToPK0s", {"pK^{0}_{s}", "Lambda_c+", "#Lambda_{c}^{+} #rightarrow pK^{0}_{s} + c.c."}},
    {"Dstar", {"D^{0}pi^{+}", "D*+", "D^{*+} #rightarrow D^{0}#pi^{+} + c.c."}},
    {"XicToXiPiPi", {"#Xi#pi#pi", "Xi_c+", "#Xi_{c}^{+} #rightarrow #Xi^{-}#pi^{+}#pi^{+} + c.c."}}};
  if (particles.find(particleName.Data()) == particles.end()) {
    throw std::runtime_error("ERROR: only Dplus, D0, Ds, LcToPKPi, LcToPK0s, Dstar and XicToXiPiPi particles supported! Exit");
  }
  const auto& particleTuple = particles[particleName.Data()];
  const TString massAxisTitle = "#it{M}(" + std::get<0>(particleTuple) + ") (GeV/#it{c}^{2})";
  const double massPDG = TDatabasePDG::Instance()->GetParticle(std::get<1>(particleTuple).c_str())->Mass();
  const std::vector<std::string> plotLabels = {std::get<2>(particleTuple), collisionSystem.Data()};

  // load inv-mass histograms
  auto* inputFile = TFile::Open(inputFileName.Data());
  if ((inputFile == nullptr) || !inputFile->IsOpen()) {
    return -1;
  }

  TFile* inputFileRefl = nullptr;
  if (enableRefl) {
    inputFileRefl = TFile::Open(reflFileName.Data());
    if ((inputFileRefl == nullptr) || !inputFileRefl->IsOpen()) {
      return -1;
    }
  }

  std::vector<TH1*> hMassSgn(nSliceVarBins);
  std::vector<TH1*> hMassRefl(nSliceVarBins);
  std::vector<TH1*> hMass(nSliceVarBins);

  for (unsigned int iSliceVar = 0; iSliceVar < nSliceVarBins; iSliceVar++) {
    if (!isMc) {
      hMass[iSliceVar] = inputFile->Get<TH1>(inputHistoName[iSliceVar].data());
      if (enableRefl) {
        hMassRefl[iSliceVar] = inputFileRefl->Get<TH1>(reflHistoName[iSliceVar].data());
        if (hMassRefl[iSliceVar] == nullptr) {
          throw std::runtime_error("ERROR: MC reflection histogram not found! Exit!");
        }
        hMassSgn[iSliceVar] = inputFileRefl->Get<TH1>(fdHistoName[iSliceVar].data());
        if (hMassSgn[iSliceVar] == nullptr) {
          throw std::runtime_error("ERROR: MC prompt or FD histogram not found! Exit!");
        }
        hMassSgn[iSliceVar]->Add(inputFileRefl->Get<TH1>(promptHistoName[iSliceVar].data()));
      }
    } else {
      hMass[iSliceVar] = inputFile->Get<TH1>(promptHistoName[iSliceVar].data());
      hMass[iSliceVar]->Add(inputFile->Get<TH1>(fdHistoName[iSliceVar].data()));
      if (includeSecPeak) {
        hMass[iSliceVar]->Add(inputFile->Get<TH1>(promptSecPeakHistoName[iSliceVar].data()));
        hMass[iSliceVar]->Add(inputFile->Get<TH1>(fdSecPeakHistoName[iSliceVar].data()));
      }
    }
    if (hMass[iSliceVar] == nullptr) {
      throw std::runtime_error("ERROR: input histogram for fit not found! Exit!");
    }
    hMass[iSliceVar]->SetDirectory(nullptr);
  }
  inputFile->Close();

  // define output histos
  auto* hRawYieldsSignal = new TH1D("hRawYieldsSignal", ";" + sliceVarName + "(" + sliceVarUnit + ");raw yield", nSliceVarBins, sliceVarLimits.data());
  auto* hRawYieldsSignalCounted = new TH1D("hRawYieldsSignalCounted", ";" + sliceVarName + "(" + sliceVarUnit + ");raw yield via bin count", nSliceVarBins, sliceVarLimits.data());
  auto* hRawYieldsBkg = new TH1D("hRawYieldsBkg", ";" + sliceVarName + "(" + sliceVarUnit + ");Background (3#sigma)", nSliceVarBins, sliceVarLimits.data());
  auto* hRawYieldsSgnOverBkg = new TH1D("hRawYieldsSgnOverBkg", ";" + sliceVarName + "(" + sliceVarUnit + ");S/B (3#sigma)", nSliceVarBins, sliceVarLimits.data());
  auto* hRawYieldsSignificance = new TH1D("hRawYieldsSignificance", ";" + sliceVarName + "(" + sliceVarUnit + ");significance (3#sigma)", nSliceVarBins, sliceVarLimits.data());
  auto* hRawYieldsChiSquareBkg = new TH1D("hRawYieldsChiSquareBkg", ";" + sliceVarName + "(" + sliceVarUnit + ");#chi^{2}/#it{ndf}", nSliceVarBins, sliceVarLimits.data());
  auto* hRawYieldsChiSquareTotal = new TH1D("hRawYieldsChiSquareTotal", ";" + sliceVarName + "(" + sliceVarUnit + ");#chi^{2}/#it{ndf}", nSliceVarBins, sliceVarLimits.data());
  auto* hReflectionOverSignal = new TH1D("hReflectionOverSignal", ";" + sliceVarName + "(" + sliceVarUnit + ");Refl/Signal", nSliceVarBins, sliceVarLimits.data());
  auto* hRawYieldsMean = new TH1D("hRawYieldsMean", ";" + sliceVarName + "(" + sliceVarUnit + ");mean (GeV/#it{c}^{2})", nSliceVarBins, sliceVarLimits.data());
  auto* hRawYieldsSigma = new TH1D("hRawYieldsSigma", ";" + sliceVarName + "(" + sliceVarUnit + ");width (GeV/#it{c}^{2})", nSliceVarBins, sliceVarLimits.data());
  auto* hRawYieldsSecSigma = new TH1D("hRawYieldsSecSigma", ";" + sliceVarName + "(" + sliceVarUnit + ");width (GeV/#it{c}^{2})", nSliceVarBins, sliceVarLimits.data());
  auto* hRawYieldsFracDoubleGaus = new TH1D("hRawYieldsFracDoubleGaus", ";" + sliceVarName + "(" + sliceVarUnit + ");fraction of double gaussian", nSliceVarBins, sliceVarLimits.data());

  const Int_t nConfigsToSave = 6;
  auto* hFitConfig = new TH2F("hfitConfig", "Fit Configurations", nConfigsToSave, 0, 6, nSliceVarBins, sliceVarLimits.data());
  const char* hFitConfigXLabel[nConfigsToSave] = {"mass min", "mass max", "rebin num", "fix sigma", "bkg func", "sgn func"};
  hFitConfig->SetStats(false);
  hFitConfig->LabelsDeflate("X");
  hFitConfig->LabelsDeflate("Y");
  hFitConfig->LabelsOption("v");
  for (int i = 0; i < nConfigsToSave; i++) {
    hFitConfig->GetXaxis()->SetBinLabel(i + 1, hFitConfigXLabel[i]);
  }

  setHistoStyle(hRawYieldsSignal);
  setHistoStyle(hRawYieldsSignalCounted);
  setHistoStyle(hRawYieldsBkg);
  setHistoStyle(hRawYieldsSgnOverBkg);
  setHistoStyle(hRawYieldsSignificance);
  setHistoStyle(hRawYieldsChiSquareBkg);
  setHistoStyle(hRawYieldsChiSquareTotal);
  setHistoStyle(hReflectionOverSignal, kRed + 1);
  setHistoStyle(hRawYieldsMean);
  setHistoStyle(hRawYieldsSigma);
  setHistoStyle(hRawYieldsSecSigma);
  setHistoStyle(hRawYieldsFracDoubleGaus);

  auto getHistToFix = [&nSliceVarBins](bool const& isFix, std::vector<double> const& fixManual, std::string const& fixFileName, std::string const& var) -> TH1* {
    TH1* histToFix = nullptr;
    if (isFix) {
      if (fixManual.empty()) {
        auto* fixInputFile = TFile::Open(fixFileName.data());
        if (fixInputFile == nullptr) {
          throw std::runtime_error("Cannot open file for fixed " + var);
        }
        const std::string histName = "hRawYields" + var;
        histToFix = fixInputFile->Get<TH1>(histName.data());
        histToFix->SetDirectory(nullptr);
        if (static_cast<unsigned int>(histToFix->GetNbinsX()) != nSliceVarBins) {
          throw std::runtime_error("Different number of bins for this analysis and histo for fixed " + var);
        }
        fixInputFile->Close();
      }
    }
    return histToFix;
  };

  TH1* hSigmaToFix = getHistToFix(fixSigma, fixSigmaManual, sigmaFile, "Sigma");
  TH1* hMeanToFix = getHistToFix(fixMean, fixMeanManual, meanFile, "Mean");
  TH1* hSecondSigmaToFix = getHistToFix(fixSecondSigma, fixSecondSigmaManual, secondSigmaFile, "SecSigma");
  TH1* hFracDoubleGausToFix = getHistToFix(fixFracDoubleGaus, fixFracDoubleGausManual, fracDoubleGausFile, "FracDoubleGaus");

  // fit histograms

  std::vector<TH1*> hMassForFit(nSliceVarBins);
  std::vector<TH1*> hMassForRefl(nSliceVarBins);
  std::vector<TH1*> hMassForSgn(nSliceVarBins);

  Int_t canvasSize[2] = {1920, 1080};
  if (nSliceVarBins == 1) {
    canvasSize[0] = 500;
    canvasSize[1] = 500;
  }

  Int_t const nCanvasesMax = 20; // do not put more than 20 bins per canvas to make them visible
  const Int_t nCanvases = std::ceil(static_cast<float>(nSliceVarBins) / nCanvasesMax);
  std::vector<TCanvas*> canvasMass(nCanvases);
  std::vector<TCanvas*> canvasResiduals(nCanvases);
  std::vector<TCanvas*> canvasRatio(nCanvases);
  std::vector<TCanvas*> canvasRefl(nCanvases);
  for (int iCanvas = 0; iCanvas < nCanvases; iCanvas++) {
    const int nPads = (nCanvases == 1) ? nSliceVarBins : nCanvasesMax;
    canvasMass[iCanvas] = new TCanvas(Form("canvasMass%d", iCanvas), Form("canvasMass%d", iCanvas),
                                      canvasSize[0], canvasSize[1]);
    divideCanvas(canvasMass[iCanvas], nPads);

    canvasResiduals[iCanvas] =
      new TCanvas(Form("canvasResiduals%d", iCanvas), Form("canvasResiduals%d", iCanvas), canvasSize[0], canvasSize[1]);
    divideCanvas(canvasResiduals[iCanvas], nPads);

    canvasRatio[iCanvas] = new TCanvas(Form("canvasRatio%d", iCanvas), Form("canvasRatio%d", iCanvas),
                                       canvasSize[0], canvasSize[1]);
    divideCanvas(canvasRatio[iCanvas], nPads);

    if (enableRefl) {
      canvasRefl[iCanvas] = new TCanvas(Form("canvasRefl%d", iCanvas), Form("canvasRefl%d", iCanvas),
                                        canvasSize[0], canvasSize[1]);
      divideCanvas(canvasRefl[iCanvas], nPads);
    }
  }

  for (unsigned int iSliceVar = 0; iSliceVar < nSliceVarBins; iSliceVar++) {
    const Int_t iCanvas = std::floor(static_cast<float>(iSliceVar) / nCanvasesMax);

    hMassForFit[iSliceVar] = hMass[iSliceVar]->Rebin(nRebin[iSliceVar]);
    TString const ptTitle =
      Form("%0.2f < " + sliceVarName + " < %0.2f " + sliceVarUnit, sliceVarMin[iSliceVar], sliceVarMax[iSliceVar]);
    hMassForFit[iSliceVar]->SetTitle(Form("%s;%s;Counts per %0.1f MeV/#it{c}^{2}",
                                          ptTitle.Data(), massAxisTitle.Data(),
                                          hMassForFit[iSliceVar]->GetBinWidth(1) * 1000));
    hMassForFit[iSliceVar]->SetName(Form("MassForFit%d", iSliceVar));

    if (enableRefl) {
      hMassForRefl[iSliceVar] = hMassRefl[iSliceVar]->Rebin(nRebin[iSliceVar]);
      hMassForSgn[iSliceVar] = hMassSgn[iSliceVar]->Rebin(nRebin[iSliceVar]);
    }

    Double_t reflOverSgn = 0;

    if (isMc) {
      HFInvMassFitter* massFitter;
      massFitter = new HFInvMassFitter(hMassForFit[iSliceVar], massMin[iSliceVar], massMax[iSliceVar], HFInvMassFitter::NoBkg, sgnFunc[iSliceVar]);
      massFitter->setDrawBgPrefit(drawBgPrefit);
      massFitter->setHighlightPeakRegion(highlightPeakRegion);
      massFitter->setInitialGaussianMean(massPDG);
      massFitter->setParticlePdgMass(massPDG);
      massFitter->setBoundGaussianMean(massPDG, 0.8 * massPDG, 1.2 * massPDG);
      massFitter->doFit();

      if (nSliceVarBins > 1) {
        canvasMass[iCanvas]->cd(iSliceVar - nCanvasesMax * iCanvas + 1);
      } else {
        canvasMass[iCanvas]->cd();
      }

      massFitter->drawFit(gPad, plotLabels, writeSignalPar);

      const Double_t rawYield = massFitter->getRawYield();
      const Double_t rawYieldErr = massFitter->getRawYieldError();
      const Double_t rawYieldCounted = massFitter->getRawYieldCounted();
      const Double_t rawYieldCountedErr = massFitter->getRawYieldCountedError();
      const Double_t reducedChiSquareBkg = massFitter->getChiSquareOverNDFBkg();
      const Double_t reducedChiSquareTotal = massFitter->getChiSquareOverNDFTotal();
      const Double_t mean = massFitter->getMean();
      const Double_t meanErr = massFitter->getMeanUncertainty();
      const Double_t sigma = massFitter->getSigma();
      const Double_t sigmaErr = massFitter->getSigmaUncertainty();

      hRawYieldsSignal->SetBinContent(iSliceVar + 1, rawYield);
      hRawYieldsSignal->SetBinError(iSliceVar + 1, rawYieldErr);
      hRawYieldsSignalCounted->SetBinContent(iSliceVar + 1, rawYieldCounted);
      hRawYieldsSignalCounted->SetBinError(iSliceVar + 1, rawYieldCountedErr);
      hRawYieldsChiSquareBkg->SetBinContent(iSliceVar + 1, reducedChiSquareBkg);
      hRawYieldsChiSquareBkg->SetBinError(iSliceVar + 1, 0.);
      hRawYieldsChiSquareTotal->SetBinContent(iSliceVar + 1, reducedChiSquareTotal);
      hRawYieldsChiSquareTotal->SetBinError(iSliceVar + 1, 0.);
      hRawYieldsMean->SetBinContent(iSliceVar + 1, mean);
      hRawYieldsMean->SetBinError(iSliceVar + 1, meanErr);
      hRawYieldsSigma->SetBinContent(iSliceVar + 1, sigma);
      hRawYieldsSigma->SetBinError(iSliceVar + 1, sigmaErr);

      if (sgnFunc[iSliceVar] != HFInvMassFitter::SingleGaus) {
        const Double_t secSigma = massFitter->getSecSigma();
        const Double_t secSigmaErr = massFitter->getSecSigmaUncertainty();
        hRawYieldsSecSigma->SetBinContent(iSliceVar + 1, secSigma);
        hRawYieldsSecSigma->SetBinError(iSliceVar + 1, secSigmaErr);
      }
      if (sgnFunc[iSliceVar] == HFInvMassFitter::DoubleGaus || sgnFunc[iSliceVar] == HFInvMassFitter::DoubleGausSigmaRatioPar) {
        const Double_t fracDoubleGaus = massFitter->getFracDoubleGaus();
        const Double_t fracDoubleGausErr = massFitter->getFracDoubleGausUncertainty();
        hRawYieldsFracDoubleGaus->SetBinContent(iSliceVar + 1, fracDoubleGaus);
        hRawYieldsFracDoubleGaus->SetBinError(iSliceVar + 1, fracDoubleGausErr);
      }
    } else {
      HFInvMassFitter* massFitter;
      massFitter = new HFInvMassFitter(hMassForFit[iSliceVar], massMin[iSliceVar], massMax[iSliceVar],
                                       bkgFunc[iSliceVar], sgnFunc[iSliceVar]);
      massFitter->setDrawBgPrefit(drawBgPrefit);
      massFitter->setHighlightPeakRegion(highlightPeakRegion);
      massFitter->setInitialGaussianMean(massPDG);
      massFitter->setParticlePdgMass(massPDG);
      massFitter->setBoundGaussianMean(massPDG, 0.8 * massPDG, 1.2 * massPDG);
      if (useLikelihood) {
        massFitter->setUseLikelihoodFit();
      }

      auto setFixedValue = [&iSliceVar](bool const& isFix, std::vector<double> const& fixManual, const TH1* histToFix, std::function<void(Double_t)> setFunc, std::string const& var) -> void {
        if (isFix) {
          if (fixManual.empty()) {
            setFunc(histToFix->GetBinContent(iSliceVar + 1));
            printf("*****************************\n");
            printf("FIXED %s: %f\n", var.data(), histToFix->GetBinContent(iSliceVar + 1));
            printf("*****************************\n");
          } else {
            setFunc(fixManual[iSliceVar]);
            printf("*****************************\n");
            printf("FIXED %s: %f\n", var.data(), fixManual[iSliceVar]);
            printf("*****************************\n");
          }
        }
      };

      setFixedValue(fixMean, fixMeanManual, hMeanToFix, std::bind(&HFInvMassFitter::setFixGaussianMean, massFitter, std::placeholders::_1), "MEAN");
      setFixedValue(fixSigma, fixSigmaManual, hSigmaToFix, std::bind(&HFInvMassFitter::setFixGaussianSigma, massFitter, std::placeholders::_1), "SIGMA");
      setFixedValue(fixSecondSigma, fixSecondSigmaManual, hSecondSigmaToFix, std::bind(&HFInvMassFitter::setFixSecondGaussianSigma, massFitter, std::placeholders::_1), "SECOND SIGMA");
      setFixedValue(fixFracDoubleGaus, fixFracDoubleGausManual, hFracDoubleGausToFix, std::bind(&HFInvMassFitter::setFixFrac2Gaus, massFitter, std::placeholders::_1), "FRAC DOUBLE GAUS");

      if (enableRefl) {
        reflOverSgn = hMassForSgn[iSliceVar]->Integral(hMassForSgn[iSliceVar]->FindBin(massMin[iSliceVar] * 1.0001), hMassForSgn[iSliceVar]->FindBin(massMax[iSliceVar] * 0.999));
        reflOverSgn = hMassForRefl[iSliceVar]->Integral(hMassForRefl[iSliceVar]->FindBin(massMin[iSliceVar] * 1.0001), hMassForRefl[iSliceVar]->FindBin(massMax[iSliceVar] * 0.999)) / reflOverSgn;
        massFitter->setFixReflOverSgn(reflOverSgn);
        massFitter->setTemplateReflections(hMassRefl[iSliceVar]);
      }

      massFitter->doFit();

      const double rawYield = massFitter->getRawYield();
      const double rawYieldErr = massFitter->getRawYieldError();
      const double rawYieldCounted = massFitter->getRawYieldCounted();
      const double rawYieldCountedErr = massFitter->getRawYieldCountedError();
      const double bkg = massFitter->getBkgYield();
      const double bkgErr = massFitter->getBkgYieldError();
      const double significance = massFitter->getSignificance();
      const double significanceErr = massFitter->getSignificanceError();
      const double reducedChiSquareBkg = massFitter->getChiSquareOverNDFBkg();
      const double reducedChiSquareTotal = massFitter->getChiSquareOverNDFTotal();
      const double mean = massFitter->getMean();
      const double meanErr = massFitter->getMeanUncertainty();
      const double sigma = massFitter->getSigma();
      const double sigmaErr = massFitter->getSigmaUncertainty();

      hRawYieldsSignal->SetBinContent(iSliceVar + 1, rawYield);
      hRawYieldsSignal->SetBinError(iSliceVar + 1, rawYieldErr);
      hRawYieldsSignalCounted->SetBinContent(iSliceVar + 1, rawYieldCounted);
      hRawYieldsSignalCounted->SetBinError(iSliceVar + 1, rawYieldCountedErr);
      hRawYieldsBkg->SetBinContent(iSliceVar + 1, bkg);
      hRawYieldsBkg->SetBinError(iSliceVar + 1, bkgErr);
      hRawYieldsSgnOverBkg->SetBinContent(iSliceVar + 1, rawYield / bkg);
      hRawYieldsSgnOverBkg->SetBinError(iSliceVar + 1, rawYield / bkg * std::sqrt(rawYieldErr / rawYield * rawYieldErr / rawYield + bkgErr / bkg * bkgErr / bkg));
      hRawYieldsSignificance->SetBinContent(iSliceVar + 1, significance);
      hRawYieldsSignificance->SetBinError(iSliceVar + 1, significanceErr);
      hRawYieldsChiSquareBkg->SetBinContent(iSliceVar + 1, reducedChiSquareBkg);
      hRawYieldsChiSquareBkg->SetBinError(iSliceVar + 1, 1.e-20);
      hRawYieldsChiSquareTotal->SetBinContent(iSliceVar + 1, reducedChiSquareTotal);
      hRawYieldsChiSquareTotal->SetBinError(iSliceVar + 1, 1.e-20);
      hRawYieldsMean->SetBinContent(iSliceVar + 1, mean);
      hRawYieldsMean->SetBinError(iSliceVar + 1, meanErr);
      hRawYieldsSigma->SetBinContent(iSliceVar + 1, sigma);
      hRawYieldsSigma->SetBinError(iSliceVar + 1, sigmaErr);

      if (sgnFunc[iSliceVar] != HFInvMassFitter::SingleGaus) {
        const double secSigma = massFitter->getSecSigma();
        const double secSigmaErr = massFitter->getSecSigmaUncertainty();
        hRawYieldsSecSigma->SetBinContent(iSliceVar + 1, secSigma);
        hRawYieldsSecSigma->SetBinError(iSliceVar + 1, secSigmaErr);
      }
      if (sgnFunc[iSliceVar] == HFInvMassFitter::DoubleGaus || sgnFunc[iSliceVar] == HFInvMassFitter::DoubleGausSigmaRatioPar) {
        const double fracDoubleGaus = massFitter->getFracDoubleGaus();
        const double fracDoubleGausErr = massFitter->getFracDoubleGausUncertainty();
        hRawYieldsFracDoubleGaus->SetBinContent(iSliceVar + 1, fracDoubleGaus);
        hRawYieldsFracDoubleGaus->SetBinError(iSliceVar + 1, fracDoubleGausErr);
      }

      if (enableRefl) {
        hReflectionOverSignal->SetBinContent(iSliceVar + 1, reflOverSgn);
        if (nSliceVarBins > 1) {
          canvasRefl[iCanvas]->cd(iSliceVar - nCanvasesMax * iCanvas + 1);
        } else {
          canvasRefl[iCanvas]->cd();
        }
        massFitter->drawReflection(gPad);
        canvasRefl[iCanvas]->Modified();
        canvasRefl[iCanvas]->Update();
      }

      if (nSliceVarBins > 1) {
        canvasMass[iCanvas]->cd(iSliceVar - nCanvasesMax * iCanvas + 1);
      } else {
        canvasMass[iCanvas]->cd();
      }
      massFitter->drawFit(gPad, plotLabels, writeSignalPar);
      canvasMass[iCanvas]->Modified();
      canvasMass[iCanvas]->Update();

      if (bkgFunc[iSliceVar] != HFInvMassFitter::NoBkg) {
        if (nSliceVarBins > 1) {
          canvasResiduals[iCanvas]->cd(iSliceVar - nCanvasesMax * iCanvas + 1);
        } else {
          canvasResiduals[iCanvas]->cd();
        }
        massFitter->drawResidual(gPad);
        canvasResiduals[iCanvas]->Modified();
        canvasResiduals[iCanvas]->Update();

        if (nSliceVarBins > 1) {
          canvasRatio[iCanvas]->cd(iSliceVar - nCanvasesMax * iCanvas + 1);
        } else {
          canvasRatio[iCanvas]->cd();
        }
        massFitter->drawRatio(gPad);
        canvasRatio[iCanvas]->Modified();
        canvasRatio[iCanvas]->Update();
      }
    }

    hFitConfig->SetBinContent(1, iSliceVar + 1, massMin[iSliceVar]);
    hFitConfig->SetBinContent(2, iSliceVar + 1, massMax[iSliceVar]);
    hFitConfig->SetBinContent(3, iSliceVar + 1, nRebin[iSliceVar]);
    if (fixSigma) {
      if (fixSigmaManual.empty()) {
        hFitConfig->SetBinContent(4, iSliceVar + 1, hSigmaToFix->GetBinContent(iSliceVar + 1));
      } else {
        hFitConfig->SetBinContent(4, iSliceVar + 1, fixSigmaManual[iSliceVar]);
      }
    }
    hFitConfig->SetBinContent(5, iSliceVar + 1, bkgFuncConfig[iSliceVar]);
    hFitConfig->SetBinContent(6, iSliceVar + 1, sgnFuncConfig[iSliceVar]);
  }

  // save output histograms
  TFile outputFile(outputFileName.Data(), "recreate");
  for (int iCanvas = 0; iCanvas < nCanvases; iCanvas++) {
    canvasMass[iCanvas]->Write();
    if (!isMc) {
      canvasResiduals[iCanvas]->Write();
      canvasRatio[iCanvas]->Write();
      if (enableRefl) {
        canvasRefl[iCanvas]->Write();
      }
    }
  }

  for (unsigned int iSliceVar = 0; iSliceVar < nSliceVarBins; iSliceVar++) {
    hMass[iSliceVar]->Write();
  }
  hRawYieldsSignal->Write();
  hRawYieldsSignalCounted->Write();
  hRawYieldsBkg->Write();
  hRawYieldsSgnOverBkg->Write();
  hRawYieldsSignificance->Write();
  hRawYieldsChiSquareBkg->Write();
  hRawYieldsChiSquareTotal->Write();
  hRawYieldsMean->Write();
  hRawYieldsSigma->Write();
  hRawYieldsSecSigma->Write();
  hRawYieldsFracDoubleGaus->Write();
  if (enableRefl) {
    hReflectionOverSignal->Write();
  }
  hFitConfig->Write();

  outputFile.Close();

  outputFileName.ReplaceAll(".root", ".pdf");
  TString outputFileNameResidual = outputFileName;
  outputFileNameResidual.ReplaceAll(".pdf", "_Residuals.pdf");
  TString outputFileRatio = outputFileName;
  outputFileRatio.ReplaceAll(".pdf", "_Ratio.pdf");
  for (int iCanvas = 0; iCanvas < nCanvases; iCanvas++) {
    if (iCanvas == 0 && nCanvases > 1) {
      canvasMass[iCanvas]->SaveAs(Form("%s[", outputFileName.Data()));
    }
    canvasMass[iCanvas]->SaveAs(outputFileName.Data());
    if (iCanvas == nCanvases - 1 && nCanvases > 1) {
      canvasMass[iCanvas]->SaveAs(Form("%s]", outputFileName.Data()));
    }
    if (!isMc) {
      // residuals
      if (iCanvas == 0 && nCanvases > 1) {
        canvasResiduals[iCanvas]->SaveAs(Form("%s[", outputFileNameResidual.Data()));
      }
      canvasResiduals[iCanvas]->SaveAs(outputFileNameResidual.Data());
      if (iCanvas == nCanvases - 1 && nCanvases > 1) {
        canvasResiduals[iCanvas]->SaveAs(Form("%s]", outputFileNameResidual.Data()));
      }
      // ratio
      if (iCanvas == 0 && nCanvases > 1) {
        canvasRatio[iCanvas]->SaveAs(Form("%s[", outputFileRatio.Data()));
      }
      canvasRatio[iCanvas]->SaveAs(outputFileRatio.Data());
      if (iCanvas == nCanvases - 1 && nCanvases > 1) {
        canvasRatio[iCanvas]->SaveAs(Form("%s]", outputFileRatio.Data()));
      }
    }
  }
  return 0;
}

void setHistoStyle(TH1* histo, Color_t color, Size_t markerSize)
{
  histo->SetStats(kFALSE);
  histo->SetMarkerSize(markerSize);
  histo->SetMarkerStyle(20);
  histo->SetLineWidth(2);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
}

void divideCanvas(TCanvas* canvas, int nSliceVarBins)
{
  int nCols = std::ceil(std::sqrt(nSliceVarBins));
  int nRows = std::ceil(static_cast<double>(nSliceVarBins) / nCols);
  canvas->Divide(nCols, nRows);
}

int main(int argc, const char* argv[])
{
  if (argc == 1) {
    throw std::runtime_error("Not enough arguments. Please use\n./runMassFitter configFileName");
  }

  const std::string configFileName = argv[1];

  runMassFitter(configFileName);

  return 0;
}
