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

#if !defined(__CINT__) || defined(__CLING__)

#include "HFInvMassFitter.h"

// if .h file not found, please include your local rapidjson/document.h and rapidjson/filereadstream.h here
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH2F.h>

#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>

#include <cstdio> // for fclose
#include <map>
#include <stdexcept>
#include <string> // std::string
#include <utility>
#include <vector> // std::vector

#endif

using namespace rapidjson;

int runMassFitter(const TString& configFileName = "config_massfitter.json");

template <typename ValueType>
void readArray(const Value& jsonArray, std::vector<ValueType>& output)
{
  for (auto it = jsonArray.Begin(); it != jsonArray.End(); it++) {
    auto value = it->template Get<ValueType>();
    output.emplace_back(value);
  }
}

void parseStringArray(const Value& jsonArray, std::vector<std::string>& output)
{
  size_t arrayLength = jsonArray.Size();
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
  if (!configFile) {
    throw std::runtime_error("ERROR: Missing configuration json file: " + configFileName);
  }

  Document config;
  char readBuffer[65536];
  FileReadStream is(configFile, readBuffer, sizeof(readBuffer));
  config.ParseStream(is);
  fclose(configFile);

  Bool_t isMc = config["IsMC"].GetBool();
  TString inputFileName = config["InFileName"].GetString();
  TString reflFileName = config["ReflFileName"].GetString();
  TString outputFileName = config["OutFileName"].GetString();
  TString particleName = config["Particle"].GetString();

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
  std::vector<double> fixSigmaManual;
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

  bool fixSigma = config["FixSigma"].GetBool();
  std::string sigmaFile = config["SigmaFile"].GetString();

  bool fixMean = config["FixMean"].GetBool();
  std::string meanFile = config["MeanFile"].GetString();

  const Value& fixSigmaManualValue = config["FixSigmaManual"];
  readArray(fixSigmaManualValue, fixSigmaManual);

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

  bool includeSecPeak = config["InclSecPeak"].GetBool();
  bool useLikelihood = config["UseLikelihood"].GetBool();

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

  std::map<std::string, std::pair<std::string, std::string>> particles{
    {"Dplus", {"K#pi#pi", "D+"}},
    {"D0", {"K#pi", "D0"}},
    {"Ds", {"KK#pi", "D_s+"}},
    {"LcToPKPi", {"pK#pi", "Lambda_c+"}},
    {"LcToPK0s", {"pK^{0}_{s}", "Lambda_c+"}},
    {"Dstar", {"D^{0}pi^{+}", "D*+"}}};
  if (particles.find(particleName.Data()) == particles.end()) {
    throw std::runtime_error("ERROR: only Dplus, D0, Ds, LcToPKPi, LcToPK0s and Dstar particles supported! Exit");
  }
  const TString massAxisTitle = "#it{M}(" + particles[particleName.Data()].first + ") (GeV/#it{c}^{2})";
  const double massPDG = TDatabasePDG::Instance()->GetParticle(particles[particleName.Data()].second.c_str())->Mass();

  // load inv-mass histograms
  auto inputFile = TFile::Open(inputFileName.Data());
  if (!inputFile || !inputFile->IsOpen()) {
    return -1;
  }

  TFile* inputFileRefl = nullptr;
  if (enableRefl) {
    inputFileRefl = TFile::Open(reflFileName.Data());
    if (!inputFileRefl || !inputFileRefl->IsOpen()) {
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
        hMassSgn[iSliceVar] = inputFileRefl->Get<TH1>(fdHistoName[iSliceVar].data());
        hMassSgn[iSliceVar]->Add(inputFileRefl->Get<TH1>(promptHistoName[iSliceVar].data()));
        if (!hMassRefl[iSliceVar]) {
          throw std::runtime_error("ERROR: MC reflection histogram not found! Exit!");
        }
        if (!hMassSgn[iSliceVar]) {
          throw std::runtime_error("ERROR: MC prompt or FD histogram not found! Exit!");
        }
      }
    } else {
      hMass[iSliceVar] = inputFile->Get<TH1>(promptHistoName[iSliceVar].data());
      hMass[iSliceVar]->Add(inputFile->Get<TH1>(fdHistoName[iSliceVar].data()));
      if (includeSecPeak) {
        hMass[iSliceVar]->Add(inputFile->Get<TH1>(promptSecPeakHistoName[iSliceVar].data()));
        hMass[iSliceVar]->Add(inputFile->Get<TH1>(fdSecPeakHistoName[iSliceVar].data()));
      }
    }
    if (!hMass[iSliceVar]) {
      throw std::runtime_error("ERROR: input histogram for fit not found! Exit!");
    }
    hMass[iSliceVar]->SetDirectory(nullptr);
  }
  inputFile->Close();

  // define output histos
  auto hRawYieldsSignal = new TH1D("hRawYieldsSignal", ";" + sliceVarName + "(" + sliceVarUnit + ");raw yield",
                                   nSliceVarBins, sliceVarLimits.data());
  auto hRawYieldsSignalCounted = new TH1D("hRawYieldsSignalCounted", ";" + sliceVarName + "(" + sliceVarUnit + ");raw yield via bin count",
                                          nSliceVarBins, sliceVarLimits.data());
  auto hRawYieldsSigma = new TH1D(
    "hRawYieldsSigma", ";" + sliceVarName + "(" + sliceVarUnit + ");width (GeV/#it{c}^{2})",
    nSliceVarBins, sliceVarLimits.data());
  auto hRawYieldsMean = new TH1D(
    "hRawYieldsMean", ";" + sliceVarName + "(" + sliceVarUnit + ");mean (GeV/#it{c}^{2})",
    nSliceVarBins, sliceVarLimits.data());
  auto hRawYieldsSignificance = new TH1D(
    "hRawYieldsSignificance",
    ";" + sliceVarName + "(" + sliceVarUnit + ");significance (3#sigma)", nSliceVarBins, sliceVarLimits.data());
  auto hRawYieldsSgnOverBkg =
    new TH1D("hRawYieldsSgnOverBkg", ";" + sliceVarName + "(" + sliceVarUnit + ");S/B (3#sigma)",
             nSliceVarBins, sliceVarLimits.data());
  auto hRawYieldsBkg =
    new TH1D("hRawYieldsBkg", ";" + sliceVarName + "(" + sliceVarUnit + ");Background (3#sigma)",
             nSliceVarBins, sliceVarLimits.data());
  auto hRawYieldsChiSquareBkg =
    new TH1D("hRawYieldsChiSquareBkg",
             ";" + sliceVarName + "(" + sliceVarUnit + ");#chi^{2}/#it{ndf}", nSliceVarBins, sliceVarLimits.data());
  auto hRawYieldsChiSquareTotal =
    new TH1D("hRawYieldsChiSquareTotal",
             ";" + sliceVarName + "(" + sliceVarUnit + ");#chi^{2}/#it{ndf}", nSliceVarBins, sliceVarLimits.data());
  auto hReflectionOverSignal =
    new TH1D("hReflectionOverSignal", ";" + sliceVarName + "(" + sliceVarUnit + ");Refl/Signal",
             nSliceVarBins, sliceVarLimits.data());

  const Int_t nConfigsToSave = 6;
  auto hFitConfig = new TH2F("hfitConfig", "Fit Configurations", nConfigsToSave, 0, 6, nSliceVarBins, sliceVarLimits.data());
  const char* hFitConfigXLabel[nConfigsToSave] = {"mass min", "mass max", "rebin num", "fix sigma", "bkg func", "sgn func"};
  hFitConfig->SetStats(0);
  hFitConfig->LabelsDeflate("X");
  hFitConfig->LabelsDeflate("Y");
  hFitConfig->LabelsOption("v");
  for (int i = 0; i < nConfigsToSave; i++) {
    hFitConfig->GetXaxis()->SetBinLabel(i + 1, hFitConfigXLabel[i]);
  }

  setHistoStyle(hRawYieldsSignal);
  setHistoStyle(hRawYieldsSignalCounted);
  setHistoStyle(hRawYieldsSigma);
  setHistoStyle(hRawYieldsMean);
  setHistoStyle(hRawYieldsSignificance);
  setHistoStyle(hRawYieldsSgnOverBkg);
  setHistoStyle(hRawYieldsBkg);
  setHistoStyle(hRawYieldsChiSquareBkg);
  setHistoStyle(hRawYieldsChiSquareTotal);
  setHistoStyle(hReflectionOverSignal, kRed + 1);

  TH1* hSigmaToFix = nullptr;
  if (fixSigma) {
    if (fixSigmaManual.empty()) {
      auto inputFileSigma = TFile::Open(sigmaFile.data());
      if (!inputFileSigma) {
        return -2;
      }
      hSigmaToFix = inputFileSigma->Get<TH1>("hRawYieldsSigma");
      hSigmaToFix->SetDirectory(0);
      if (static_cast<unsigned int>(hSigmaToFix->GetNbinsX()) != nSliceVarBins) {
        printf("WARNING: Different number of bins for this analysis and histo for fix sigma!\n");
      }
      inputFileSigma->Close();
    }
  }

  TH1* hMeanToFix = nullptr;
  if (fixMean) {
    auto inputFileMean = TFile::Open(meanFile.data());
    if (!inputFileMean) {
      return -3;
    }
    hMeanToFix = inputFileMean->Get<TH1>("hRawYieldsMean");
    hMeanToFix->SetDirectory(0);
    if (static_cast<unsigned int>(hMeanToFix->GetNbinsX()) != nSliceVarBins) {
      printf("WARNING: Different number of bins for this analysis and histo for fix mean\n");
    }
    inputFileMean->Close();
  }

  // fit histograms

  std::vector<TH1*> hMassForFit(nSliceVarBins);
  std::vector<TH1*> hMassForRefl(nSliceVarBins);
  std::vector<TH1*> hMassForSgn(nSliceVarBins);

  Int_t canvasSize[2] = {1920, 1080};
  if (nSliceVarBins == 1) {
    canvasSize[0] = 500;
    canvasSize[1] = 500;
  }

  Int_t nCanvasesMax = 20; // do not put more than 20 bins per canvas to make them visible
  const Int_t nCanvases = std::ceil(static_cast<float>(nSliceVarBins) / nCanvasesMax);
  std::vector<TCanvas*> canvasMass(nCanvases);
  std::vector<TCanvas*> canvasResiduals(nCanvases);
  std::vector<TCanvas*> canvasRefl(nCanvases);
  for (int iCanvas = 0; iCanvas < nCanvases; iCanvas++) {
    const int nPads = (nCanvases == 1) ? nSliceVarBins : nCanvasesMax;
    canvasMass[iCanvas] = new TCanvas(Form("canvasMass%d", iCanvas), Form("canvasMass%d", iCanvas),
                                      canvasSize[0], canvasSize[1]);
    divideCanvas(canvasMass[iCanvas], nPads);

    canvasResiduals[iCanvas] =
      new TCanvas(Form("canvasResiduals%d", iCanvas), Form("canvasResiduals%d", iCanvas), canvasSize[0], canvasSize[1]);
    divideCanvas(canvasResiduals[iCanvas], nPads);
    canvasRefl[iCanvas] = new TCanvas(Form("canvasRefl%d", iCanvas), Form("canvasRefl%d", iCanvas),
                                      canvasSize[0], canvasSize[1]);
    divideCanvas(canvasRefl[iCanvas], nPads);
  }

  for (unsigned int iSliceVar = 0; iSliceVar < nSliceVarBins; iSliceVar++) {
    const Int_t iCanvas = std::floor(static_cast<float>(iSliceVar) / nCanvasesMax);

    hMassForFit[iSliceVar] = static_cast<TH1*>(hMass[iSliceVar]->Rebin(nRebin[iSliceVar]));
    TString ptTitle =
      Form("%0.2f < " + sliceVarName + " < %0.2f " + sliceVarUnit, sliceVarMin[iSliceVar], sliceVarMax[iSliceVar]);
    hMassForFit[iSliceVar]->SetTitle(Form("%s;%s;Counts per %0.1f MeV/#it{c}^{2}",
                                          ptTitle.Data(), massAxisTitle.Data(),
                                          hMassForFit[iSliceVar]->GetBinWidth(1) * 1000));
    hMassForFit[iSliceVar]->SetName(Form("MassForFit%d", iSliceVar));

    if (enableRefl) {
      hMassForRefl[iSliceVar] = static_cast<TH1*>(hMassRefl[iSliceVar]->Rebin(nRebin[iSliceVar]));
      hMassForSgn[iSliceVar] = static_cast<TH1*>(hMassSgn[iSliceVar]->Rebin(nRebin[iSliceVar]));
    }

    Double_t reflOverSgn = 0;
    double markerSize = 1.;
    constexpr int NSliceVarBinsLarge = 15;
    if (nSliceVarBins > NSliceVarBinsLarge) {
      markerSize = 0.5;
    }

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

      massFitter->drawFit(gPad);

      const Double_t rawYield = massFitter->getRawYield();
      const Double_t rawYieldErr = massFitter->getRawYieldError();
      const Double_t rawYieldCounted = massFitter->getRawYieldCounted();
      const Double_t rawYieldCountedErr = massFitter->getRawYieldCountedError();

      const Double_t sigma = massFitter->getSigma();
      const Double_t sigmaErr = massFitter->getSigmaUncertainty();
      const Double_t mean = massFitter->getMean();
      const Double_t meanErr = massFitter->getMeanUncertainty();
      const Double_t reducedChiSquareBkg = massFitter->getChiSquareOverNDFBkg();
      const Double_t reducedChiSquareTotal = massFitter->getChiSquareOverNDFTotal();

      hRawYieldsSignal->SetBinContent(iSliceVar + 1, rawYield);
      hRawYieldsSignal->SetBinError(iSliceVar + 1, rawYieldErr);
      hRawYieldsSignalCounted->SetBinContent(iSliceVar + 1, rawYieldCounted);
      hRawYieldsSignalCounted->SetBinError(iSliceVar + 1, rawYieldCountedErr);
      hRawYieldsSigma->SetBinContent(iSliceVar + 1, sigma);
      hRawYieldsSigma->SetBinError(iSliceVar + 1, sigmaErr);
      hRawYieldsMean->SetBinContent(iSliceVar + 1, mean);
      hRawYieldsMean->SetBinError(iSliceVar + 1, meanErr);
      hRawYieldsChiSquareBkg->SetBinContent(iSliceVar + 1, reducedChiSquareBkg);
      hRawYieldsChiSquareBkg->SetBinError(iSliceVar + 1, 0.);
      hRawYieldsChiSquareTotal->SetBinContent(iSliceVar + 1, reducedChiSquareTotal);
      hRawYieldsChiSquareTotal->SetBinError(iSliceVar + 1, 0.);
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
      if (fixMean) {
        massFitter->setFixGaussianMean(hMeanToFix->GetBinContent(iSliceVar + 1));
      }
      if (fixSigma) {
        if (fixSigmaManual.empty()) {
          massFitter->setFixGaussianSigma(hSigmaToFix->GetBinContent(iSliceVar + 1));
          printf("*****************************\n");
          printf("FIXED SIGMA: %f\n", hSigmaToFix->GetBinContent(iSliceVar + 1));
          printf("*****************************\n");
        } else if (!fixSigmaManual.empty()) {
          massFitter->setFixGaussianSigma(fixSigmaManual[iSliceVar]);
          printf("*****************************\n");
          printf("FIXED SIGMA: %f\n", fixSigmaManual[iSliceVar]);
          printf("*****************************\n");
        } else {
          printf("WARNING: impossible to fix sigma! Wrong fix sigma file or value!\n");
        }
      }

      if (enableRefl) {
        reflOverSgn = hMassForSgn[iSliceVar]->Integral(hMassForSgn[iSliceVar]->FindBin(massMin[iSliceVar] * 1.0001), hMassForSgn[iSliceVar]->FindBin(massMax[iSliceVar] * 0.999));
        reflOverSgn = hMassForRefl[iSliceVar]->Integral(hMassForRefl[iSliceVar]->FindBin(massMin[iSliceVar] * 1.0001), hMassForRefl[iSliceVar]->FindBin(massMax[iSliceVar] * 0.999)) / reflOverSgn;
        massFitter->setFixReflOverSgn(reflOverSgn);
        massFitter->setTemplateReflections(hMassRefl[iSliceVar], HFInvMassFitter::DoubleGaus);
      }

      massFitter->doFit();

      const double rawYield = massFitter->getRawYield();
      const double rawYieldErr = massFitter->getRawYieldError();
      const double rawYieldCounted = massFitter->getRawYieldCounted();
      const double rawYieldCountedErr = massFitter->getRawYieldCountedError();
      const double sigma = massFitter->getSigma();
      const double sigmaErr = massFitter->getSigmaUncertainty();
      const double mean = massFitter->getMean();
      const double meanErr = massFitter->getMeanUncertainty();
      const double reducedChiSquareBkg = massFitter->getChiSquareOverNDFBkg();
      const double reducedChiSquareTotal = massFitter->getChiSquareOverNDFTotal();
      const double significance = massFitter->getSignificance();
      const double significanceErr = massFitter->getSignificanceError();
      const double bkg = massFitter->getBkgYield();
      const double bkgErr = massFitter->getBkgYieldError();

      hRawYieldsSignal->SetBinContent(iSliceVar + 1, rawYield);
      hRawYieldsSignal->SetBinError(iSliceVar + 1, rawYieldErr);
      hRawYieldsSignalCounted->SetBinContent(iSliceVar + 1, rawYieldCounted);
      hRawYieldsSignalCounted->SetBinError(iSliceVar + 1, rawYieldCountedErr);
      hRawYieldsSigma->SetBinContent(iSliceVar + 1, sigma);
      hRawYieldsSigma->SetBinError(iSliceVar + 1, sigmaErr);
      hRawYieldsMean->SetBinContent(iSliceVar + 1, mean);
      hRawYieldsMean->SetBinError(iSliceVar + 1, meanErr);
      hRawYieldsSignificance->SetBinContent(iSliceVar + 1, significance);
      hRawYieldsSignificance->SetBinError(iSliceVar + 1, significanceErr);
      hRawYieldsSgnOverBkg->SetBinContent(iSliceVar + 1, rawYield / bkg);
      hRawYieldsSgnOverBkg->SetBinError(iSliceVar + 1, rawYield / bkg * std::sqrt(rawYieldErr / rawYield * rawYieldErr / rawYield + bkgErr / bkg * bkgErr / bkg));
      hRawYieldsBkg->SetBinContent(iSliceVar + 1, bkg);
      hRawYieldsBkg->SetBinError(iSliceVar + 1, bkgErr);
      hRawYieldsChiSquareBkg->SetBinContent(iSliceVar + 1, reducedChiSquareBkg);
      hRawYieldsChiSquareBkg->SetBinError(iSliceVar + 1, 1.e-20);
      hRawYieldsChiSquareTotal->SetBinContent(iSliceVar + 1, reducedChiSquareTotal);
      hRawYieldsChiSquareTotal->SetBinError(iSliceVar + 1, 1.e-20);
      if (enableRefl) {
        hReflectionOverSignal->SetBinContent(iSliceVar + 1, reflOverSgn);
      }

      if (enableRefl) {
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
      massFitter->drawFit(gPad);
      canvasMass[iCanvas]->Modified();
      canvasMass[iCanvas]->Update();

      if (nSliceVarBins > 1) {
        canvasResiduals[iCanvas]->cd(iSliceVar - nCanvasesMax * iCanvas + 1);
      } else {
        canvasResiduals[iCanvas]->cd();
      }
      massFitter->drawResidual(gPad);
      canvasResiduals[iCanvas]->Modified();
      canvasResiduals[iCanvas]->Update();
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
      canvasRefl[iCanvas]->Write();
    }
  }

  for (unsigned int iSliceVar = 0; iSliceVar < nSliceVarBins; iSliceVar++) {
    hMass[iSliceVar]->Write();
  }
  hRawYieldsSignal->Write();
  hRawYieldsSignalCounted->Write();
  hRawYieldsSigma->Write();
  hRawYieldsMean->Write();
  hRawYieldsSignificance->Write();
  hRawYieldsSgnOverBkg->Write();
  hRawYieldsBkg->Write();
  hRawYieldsChiSquareBkg->Write();
  hRawYieldsChiSquareTotal->Write();
  if (enableRefl) {
    hReflectionOverSignal->Write();
  }
  hFitConfig->Write();

  outputFile.Close();

  outputFileName.ReplaceAll(".root", ".pdf");
  TString outputFileNameResidual = outputFileName;
  outputFileNameResidual.ReplaceAll(".pdf", "_Residuals.pdf");
  for (int iCanvas = 0; iCanvas < nCanvases; iCanvas++) {
    if (iCanvas == 0 && nCanvases > 1) {
      canvasMass[iCanvas]->SaveAs(Form("%s[", outputFileName.Data()));
    }
    canvasMass[iCanvas]->SaveAs(outputFileName.Data());
    if (iCanvas == nCanvases - 1 && nCanvases > 1) {
      canvasMass[iCanvas]->SaveAs(Form("%s]", outputFileName.Data()));
    }
    if (!isMc) {
      if (iCanvas == 0 && nCanvases > 1) {
        canvasResiduals[iCanvas]->SaveAs(Form("%s[", outputFileNameResidual.Data()));
      }
      canvasResiduals[iCanvas]->SaveAs(outputFileNameResidual.Data());
      if (iCanvas == nCanvases - 1 && nCanvases > 1) {
        canvasResiduals[iCanvas]->SaveAs(Form("%s]", outputFileNameResidual.Data()));
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
  const int rectangularSideMin = std::floor(std::sqrt(nSliceVarBins));
  constexpr int RectangularSidesDiffMax = 2;
  for (int rectangularSidesDiff = 0; rectangularSidesDiff < RectangularSidesDiffMax; ++rectangularSidesDiff) {
    if (rectangularSideMin * (rectangularSideMin + rectangularSidesDiff) >= nSliceVarBins) {
      canvas->Divide(rectangularSideMin + rectangularSidesDiff, rectangularSideMin);
    }
  }
}

int main(int argc, char* argv[])
{
  if (argc == 1) {
    throw std::runtime_error("Not enough arguments. Please use\n./runMassFitter configFileName");
  }

  const std::string configFileName = argv[1];

  runMassFitter(configFileName);

  return 0;
}
