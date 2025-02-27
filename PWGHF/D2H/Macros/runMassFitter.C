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

#if !defined(__CINT__) || defined(__CLING__)

#include "HFInvMassFitter.h"

#include <iostream> // std::cout
#include <string>   // std::string
#include <vector>   // std::vector

#include <Riostream.h>
#include <TROOT.h>

// if .h file not found, please include your local rapidjson/document.h and rapidjson/filereadstream.h here
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>

#endif

using namespace std;
using namespace rapidjson;

int runMassFitter(TString configFileName = "config_massfitter.json");

template <typename ValueType>
void readArray(const Value& jsonArray, std::vector<ValueType>& output)
{
  for (auto it = jsonArray.Begin(); it != jsonArray.End(); it++) {
    auto value = it->template Get<ValueType>();
    output.emplace_back(value);
  }
}

void parseStringArray(const Value& jsonArray, std::vector<string>& output)
{
  size_t arrayLength = jsonArray.Size();
  for (size_t i = 0; i < arrayLength; i++) {
    if (jsonArray[i].IsString()) {
      output.emplace_back(jsonArray[i].GetString());
    }
  }
}

void divideCanvas(TCanvas* c, int nPtBins);
void setHistoStyle(TH1* histo, int color = kBlack, double markerSize = 1.);

int runMassFitter(TString configFileName)
{
  // load config
  FILE* configFile = fopen(configFileName.Data(), "r");
  if (!configFile) {
    cerr << "ERROR: Missing configuration json file: " << configFileName << endl;
    return -1;
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

  vector<string> inputHistoName;
  vector<string> promptHistoName;
  vector<string> fdHistoName;
  vector<string> reflHistoName;
  vector<string> promptSecPeakHistoName;
  vector<string> fdSecPeakHistoName;
  vector<double> ptMin;
  vector<double> ptMax;
  vector<double> massMin;
  vector<double> massMax;
  vector<double> fixSigmaManual;
  vector<int> nRebin;
  vector<int> bkgFuncConfig;
  vector<int> sgnFuncConfig;

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
  parseStringArray(promptSecPeakHistoNameValue, promptSecPeakHistoName);

  bool fixSigma = config["FixSigma"].GetBool();
  string sigmaFile = config["SigmaFile"].GetString();
  double sigmaMultFactor =
    config["SigmaMultFactor"].GetDouble();
  bool fixMean = config["FixMean"].GetBool();
  string meanFile = config["MeanFile"].GetString();

  const Value& fixSigmaManualValue = config["FixSigmaManual"];
  readArray(fixSigmaManualValue, fixSigmaManual);

  const Value& ptMinValue = config["PtMin"];
  readArray(ptMinValue, ptMin);

  const Value& ptMaxValue = config["PtMax"];
  readArray(ptMaxValue, ptMax);

  const Value& massMinValue = config["MassMin"];
  readArray(massMinValue, massMin);

  const Value& massMaxValue = config["MassMax"];
  readArray(massMaxValue, massMax);

  const Value& rebinValue = config["Rebin"];
  readArray(rebinValue, nRebin);

  bool includeSecPeak = config["InclSecPeak"].GetBool();
  string sigmaSecPeak = config["SigmaSecPeak"].GetString();
  string sigmaFileSecPeak =
    config["SigmaFileSecPeak"].GetString();
  double sigmaMultFactorSecPeak =
    config["SigmaMultFactorSecPeak"].GetDouble();
  bool fixSigmaToFirstPeak =
    config["FixSigmaToFirstPeak"].GetBool();
  bool useLikelihood = config["UseLikelihood"].GetBool();

  const Value& bkgFuncValue = config["BkgFunc"];
  readArray(bkgFuncValue, bkgFuncConfig);

  const Value& sgnFuncValue = config["SgnFunc"];
  readArray(sgnFuncValue, sgnFuncConfig);

  bool fixSigmaRatio = config["FixSigmaRatio"].GetBool();
  TString sigmaRatioFile = config["SigmaRatioFile"].GetString();
  bool boundMean = config["BoundMean"].GetBool();
  bool enableRefl = config["EnableRefl"].GetBool();

  const unsigned int nPtBins = ptMin.size();
  int bkgFunc[nPtBins], sgnFunc[nPtBins];
  double ptLimits[nPtBins + 1];

  for (unsigned int iPt = 0; iPt < nPtBins; iPt++) {
    ptLimits[iPt] = ptMin[iPt];
    ptLimits[iPt + 1] = ptMax[iPt];

    if (bkgFuncConfig[iPt] == 0) {
      bkgFunc[iPt] = HFInvMassFitter::Expo;
    } else if (bkgFuncConfig[iPt] == 1) {
      bkgFunc[iPt] = HFInvMassFitter::Poly1;
    } else if (bkgFuncConfig[iPt] == 2) {
      bkgFunc[iPt] = HFInvMassFitter::Poly2;
    } else if (bkgFuncConfig[iPt] == 3) {
      bkgFunc[iPt] = HFInvMassFitter::Pow;
    } else if (bkgFuncConfig[iPt] == 4) {
      bkgFunc[iPt] = HFInvMassFitter::PowExpo;
    } else if (bkgFuncConfig[iPt] == 5) {
      bkgFunc[iPt] = HFInvMassFitter::Poly3;
    } else if (bkgFuncConfig[iPt] == 6) {
      bkgFunc[iPt] = HFInvMassFitter::NoBkg;
    } else {
      cerr << "ERROR: only Expo, Poly1, Poly2, Pow and PowEx background "
              "functions supported! Exit"
           << endl;
      return -1;
    }

    if (sgnFuncConfig[iPt] == 0) {
      sgnFunc[iPt] = HFInvMassFitter::SingleGaus;
    } else if (sgnFuncConfig[iPt] == 1) {
      sgnFunc[iPt] = HFInvMassFitter::DoubleGaus;
    } else if (sgnFuncConfig[iPt] == 2) {
      sgnFunc[iPt] = HFInvMassFitter::DoubleGausSigmaRatioPar;
    } else {
      cerr << "ERROR: only SingleGaus, DoubleGaus and DoubleGausSigmaRatioPar signal "
              "functions supported! Exit"
           << endl;
      return -1;
    }
  }

  TString massAxisTitle = "";
  if (particleName == "Dplus") {
    massAxisTitle = "#it{M}(K#pi#pi) (GeV/#it{c}^{2})";
  } else if (particleName == "D0") {
    massAxisTitle = "#it{M}(K#pi) (GeV/#it{c}^{2})";
  } else if (particleName == "Ds") {
    massAxisTitle = "#it{M}(KK#pi) (GeV/#it{c}^{2})";
  } else if (particleName == "LcToPKPi") {
    massAxisTitle = "#it{M}(pK#pi) (GeV/#it{c}^{2})";
  } else if (particleName == "LcToPK0s") {
    massAxisTitle = "#it{M}(pK^{0}_{s}) (GeV/#it{c}^{2})";
  } else if (particleName == "Dstar") {
    massAxisTitle = "#it{M}(pi^{+}) (GeV/#it{c}^{2})";
  } else {
    cerr << "ERROR: only Dplus, D0, Ds, LcToPKPi, LcToPK0s and Dstar particles supported! Exit" << endl;
    return -1;
  }

  // load inv-mass histograms
  auto inputFile = TFile::Open(inputFileName.Data());
  if (!inputFile || !inputFile->IsOpen()) {
    return -1;
  }

  TFile* inputFileRefl = NULL;
  if (enableRefl) {
    inputFileRefl = TFile::Open(reflFileName.Data());
    if (!inputFileRefl || !inputFileRefl->IsOpen()) {
      return -1;
    }
  }

  TH1F* hMassSgn[nPtBins];
  TH1F* hMassRefl[nPtBins];
  TH1F* hMass[nPtBins];

  for (unsigned int iPt = 0; iPt < nPtBins; iPt++) {
    if (!isMc) {
      hMass[iPt] = static_cast<TH1F*>(inputFile->Get(inputHistoName[iPt].data()));
      if (enableRefl) {
        hMassRefl[iPt] = static_cast<TH1F*>(inputFileRefl->Get(reflHistoName[iPt].data()));
        hMassSgn[iPt] = static_cast<TH1F*>(inputFileRefl->Get(fdHistoName[iPt].data()));
        hMassSgn[iPt]->Add(static_cast<TH1F*>(inputFileRefl->Get(promptHistoName[iPt].data())));
        if (!hMassRefl[iPt]) {
          cerr << "ERROR: MC reflection histogram not found! Exit!" << endl;
          return -1;
        }
        if (!hMassSgn[iPt]) {
          cerr << "ERROR: MC prompt or FD histogram not found! Exit!" << endl;
          return -1;
        }
      }
    } else {
      hMass[iPt] = static_cast<TH1F*>(inputFile->Get(promptHistoName[iPt].data()));
      hMass[iPt]->Add(static_cast<TH1F*>(inputFile->Get(fdHistoName[iPt].data())));
      if (includeSecPeak) {
        hMass[iPt]->Add(static_cast<TH1F*>(inputFile->Get(promptSecPeakHistoName[iPt].data())));
        hMass[iPt]->Add(static_cast<TH1F*>(inputFile->Get(fdSecPeakHistoName[iPt].data())));
      }
    }
    if (!hMass[iPt]) {
      cerr << "ERROR: input histogram for fit not found! Exit!" << endl;
      return -1;
    }
    hMass[iPt]->SetDirectory(0);
  }
  inputFile->Close();

  // define output histos
  auto hRawYields = new TH1D("hRawYields", ";#it{p}_{T} (GeV/#it{c});raw yield",
                             nPtBins, ptLimits);
  auto hRawYieldsSigma = new TH1D(
    "hRawYieldsSigma", ";#it{p}_{T} (GeV/#it{c});width (GeV/#it{c}^{2})",
    nPtBins, ptLimits);
  auto hRawYieldsSigmaRatio = new TH1D(
    "hRawYieldsSigmaRatio",
    ";#it{p}_{T} (GeV/#it{c});ratio #sigma_{1}/#sigma_{2}", nPtBins, ptLimits);
  auto hRawYieldsSigma2 = new TH1D(
    "hRawYieldsSigma2", ";#it{p}_{T} (GeV/#it{c});width (GeV/#it{c}^{2})",
    nPtBins, ptLimits);
  auto hRawYieldsMean = new TH1D(
    "hRawYieldsMean", ";#it{p}_{T} (GeV/#it{c});mean (GeV/#it{c}^{2})",
    nPtBins, ptLimits);
  auto hRawYieldsFracGaus2 = new TH1D(
    "hRawYieldsFracGaus2",
    ";#it{p}_{T} (GeV/#it{c});second-gaussian fraction", nPtBins, ptLimits);
  auto hRawYieldsSignificance = new TH1D(
    "hRawYieldsSignificance",
    ";#it{p}_{T} (GeV/#it{c});significance (3#sigma)", nPtBins, ptLimits);
  auto hRawYieldsSgnOverBkg =
    new TH1D("hRawYieldsSgnOverBkg", ";#it{p}_{T} (GeV/#it{c});S/B (3#sigma)",
             nPtBins, ptLimits);
  auto hRawYieldsSignal =
    new TH1D("hRawYieldsSignal", ";#it{p}_{T} (GeV/#it{c});Signal (3#sigma)",
             nPtBins, ptLimits);
  auto hRawYieldsBkg =
    new TH1D("hRawYieldsBkg", ";#it{p}_{T} (GeV/#it{c});Background (3#sigma)",
             nPtBins, ptLimits);
  auto hRawYieldsChiSquare =
    new TH1D("hRawYieldsChiSquare",
             ";#it{p}_{T} (GeV/#it{c});#chi^{2}/#it{ndf}", nPtBins, ptLimits);
  auto hRawYieldsSecondPeak = new TH1D(
    "hRawYieldsSecondPeak", ";#it{p}_{T} (GeV/#it{c});raw yield second peak",
    nPtBins, ptLimits);
  auto hRawYieldsMeanSecondPeak =
    new TH1D("hRawYieldsMeanSecondPeak",
             ";#it{p}_{T} (GeV/#it{c});mean second peak (GeV/#it{c}^{2})",
             nPtBins, ptLimits);
  auto hRawYieldsSigmaSecondPeak =
    new TH1D("hRawYieldsSigmaSecondPeak",
             ";#it{p}_{T} (GeV/#it{c});width second peak (GeV/#it{c}^{2})",
             nPtBins, ptLimits);
  auto hRawYieldsSignificanceSecondPeak =
    new TH1D("hRawYieldsSignificanceSecondPeak",
             ";#it{p}_{T} (GeV/#it{c});signficance second peak (3#sigma)",
             nPtBins, ptLimits);
  auto hRawYieldsSigmaRatioSecondFirstPeak =
    new TH1D("hRawYieldsSigmaRatioSecondFirstPeak",
             ";#it{p}_{T} (GeV/#it{c});width second peak / width first peak",
             nPtBins, ptLimits);
  auto hRawYieldsSoverBSecondPeak = new TH1D(
    "hRawYieldsSoverBSecondPeak",
    ";#it{p}_{T} (GeV/#it{c});S/B second peak (3#sigma)", nPtBins, ptLimits);
  auto hRawYieldsSignalSecondPeak = new TH1D(
    "hRawYieldsSignalSecondPeak",
    ";#it{p}_{T} (GeV/#it{c});Signal second peak (3#sigma)", nPtBins, ptLimits);
  auto hRawYieldsBkgSecondPeak =
    new TH1D("hRawYieldsBkgSecondPeak",
             ";#it{p}_{T} (GeV/#it{c});Background second peak (3#sigma)",
             nPtBins, ptLimits);
  auto hReflectionOverSignal =
    new TH1D("hReflectionOverSignal", ";#it{p}_{T} (GeV/#it{c});Refl/Signal",
             nPtBins, ptLimits);

  const Int_t nConfigsToSave = 6;
  auto hFitConfig = new TH2F("hfitConfig", "Fit Configurations", nConfigsToSave, 0, 6, nPtBins, ptLimits);
  const char* hFitConfigXLabel[nConfigsToSave] = {"mass min", "mass max", "rebin num", "fix sigma", "bkg func", "sgn func"};
  hFitConfig->SetStats(0);
  hFitConfig->LabelsDeflate("X");
  hFitConfig->LabelsDeflate("Y");
  hFitConfig->LabelsOption("v");
  for (int i = 0; i < nConfigsToSave; i++) {
    hFitConfig->GetXaxis()->SetBinLabel(i + 1, hFitConfigXLabel[i]);
  }

  setHistoStyle(hRawYields);
  setHistoStyle(hRawYieldsSigma);
  setHistoStyle(hRawYieldsSigma2);
  setHistoStyle(hRawYieldsMean);
  setHistoStyle(hRawYieldsFracGaus2);
  setHistoStyle(hRawYieldsSignificance);
  setHistoStyle(hRawYieldsSgnOverBkg);
  setHistoStyle(hRawYieldsSignal);
  setHistoStyle(hRawYieldsBkg);
  setHistoStyle(hRawYieldsChiSquare);
  setHistoStyle(hRawYieldsSecondPeak, kRed + 1);
  setHistoStyle(hRawYieldsMeanSecondPeak, kRed + 1);
  setHistoStyle(hRawYieldsSigmaSecondPeak, kRed + 1);
  setHistoStyle(hRawYieldsSignificanceSecondPeak, kRed + 1);
  setHistoStyle(hRawYieldsSigmaRatioSecondFirstPeak, kRed + 1);
  setHistoStyle(hRawYieldsSoverBSecondPeak, kRed + 1);
  setHistoStyle(hRawYieldsSignalSecondPeak, kRed + 1);
  setHistoStyle(hRawYieldsBkgSecondPeak, kRed + 1);
  setHistoStyle(hReflectionOverSignal, kRed + 1);

  TH1D* hSigmaToFix = NULL;
  if (fixSigma) {
    if (fixSigmaManual.empty()) {
      auto inputFileSigma = TFile::Open(sigmaFile.data());
      if (!inputFileSigma) {
        return -2;
      }
      hSigmaToFix = static_cast<TH1D*>(inputFileSigma->Get("hRawYieldsSigma"));
      hSigmaToFix->SetDirectory(0);
      if (static_cast<unsigned int>(hSigmaToFix->GetNbinsX()) != nPtBins) {
        cout << "WARNING: Different number of bins for this analysis and histo for fix sigma!" << endl;
      }
      inputFileSigma->Close();
    }
  }

  TH1D* hMeanToFix = NULL;
  if (fixMean) {
    auto inputFileMean = TFile::Open(meanFile.data());
    if (!inputFileMean) {
      return -3;
    }
    hMeanToFix = static_cast<TH1D*>(inputFileMean->Get("hRawYieldsMean"));
    hMeanToFix->SetDirectory(0);
    if (static_cast<unsigned int>(hMeanToFix->GetNbinsX()) != nPtBins) {
      cout << "WARNING: Different number of bins for this analysis and histo for fix mean" << endl;
    }
    inputFileMean->Close();
  }

  // fit histograms

  TH1F* hMassForFit[nPtBins];
  TH1F* hMassForRefl[nPtBins];
  TH1F* hMassForSgn[nPtBins];

  Int_t canvasSize[2] = {1920, 1080};
  if (nPtBins == 1) {
    canvasSize[0] = 500;
    canvasSize[1] = 500;
  }

  Int_t nCanvasesMax = 20; // do not put more than 20 bins per canvas to make them visible
  const Int_t nCanvases = ceil((float)nPtBins / nCanvasesMax);
  TCanvas *canvasMass[nCanvases], *canvasResiduals[nCanvases], *canvasRefl[nCanvases];
  for (int iCanvas = 0; iCanvas < nCanvases; iCanvas++) {
    int nPads = (nCanvases == 1) ? nPtBins : nCanvasesMax;
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

  for (unsigned int iPt = 0; iPt < nPtBins; iPt++) {
    Int_t iCanvas = floor((float)iPt / nCanvasesMax);

    hMassForFit[iPt] = reinterpret_cast<TH1F*>(hMass[iPt]->Rebin(nRebin[iPt]));
    TString ptTitle =
      Form("%0.1f < #it{p}_{T} < %0.1f GeV/#it{c}", ptMin[iPt], ptMax[iPt]);
    hMassForFit[iPt]->SetTitle(Form("%s;%s;Counts per %0.f MeV/#it{c}^{2}",
                                    ptTitle.Data(), massAxisTitle.Data(),
                                    hMassForFit[iPt]->GetBinWidth(1) * 1000));
    hMassForFit[iPt]->SetName(Form("MassForFit%d", iPt));

    if (enableRefl) {
      hMassForRefl[iPt] =
        reinterpret_cast<TH1F*>(hMassRefl[iPt]->Rebin(nRebin[iPt]));
      hMassForSgn[iPt] =
        reinterpret_cast<TH1F*>(hMassSgn[iPt]->Rebin(nRebin[iPt]));
    }

    Double_t reflOverSgn = 0;
    double markerSize = 1.;
    if (nPtBins > 15) {
      markerSize = 0.5;
    }

    if (isMc) {
      HFInvMassFitter* massFitter;
      massFitter = new HFInvMassFitter(hMassForFit[iPt], massMin[iPt], massMax[iPt], HFInvMassFitter::NoBkg, sgnFunc[iPt]);
      massFitter->doFit(false);

      if (nPtBins > 1) {
        canvasMass[iCanvas]->cd(iPt - nCanvasesMax * iCanvas + 1);
      } else {
        canvasMass[iCanvas]->cd();
      }

      massFitter->drawFit(gPad);

      Double_t rawYield = massFitter->getRawYield();
      Double_t rawYieldErr = massFitter->getRawYieldError();

      Double_t sigma = massFitter->getSigma();
      Double_t sigmaErr = massFitter->getSigmaUncertainty();
      Double_t mean = massFitter->getMean();
      Double_t meanErr = massFitter->getMeanUncertainty();
      Double_t reducedChiSquare = massFitter->getChiSquareOverNDF();

      hRawYields->SetBinContent(iPt + 1, rawYield);
      hRawYields->SetBinError(iPt + 1, rawYieldErr);
      hRawYieldsSigma->SetBinContent(iPt + 1, sigma);
      hRawYieldsSigma->SetBinError(iPt + 1, sigmaErr);
      hRawYieldsMean->SetBinContent(iPt + 1, mean);
      hRawYieldsMean->SetBinError(iPt + 1, meanErr);
      hRawYieldsChiSquare->SetBinContent(iPt + 1, reducedChiSquare);
      hRawYieldsChiSquare->SetBinError(iPt + 1, 0.);
    } else {
      HFInvMassFitter* massFitter;
      massFitter = new HFInvMassFitter(hMassForFit[iPt], massMin[iPt], massMax[iPt],
                                       bkgFunc[iPt], sgnFunc[iPt]);
      if (useLikelihood) {
        massFitter->setUseLikelihoodFit();
      }
      if (fixMean) {
        massFitter->setFixGaussianMean(hMeanToFix->GetBinContent(iPt + 1));
      }
      if (fixSigma) {
        if (fixSigmaManual.empty()) {
          massFitter->setFixGaussianSigma(hSigmaToFix->GetBinContent(iPt + 1));
          cout << "*****************************"
               << "\n"
               << "FIXED SIGMA: " << hSigmaToFix->GetBinContent(iPt + 1) << "\n"
               << "*****************************" << endl;
        } else if (!fixSigmaManual.empty()) {
          massFitter->setFixGaussianSigma(fixSigmaManual[iPt]);
          cout << "*****************************"
               << "\n"
               << "FIXED SIGMA: " << fixSigmaManual[iPt] << "\n"
               << "*****************************" << endl;
        } else {
          cout << "WARNING: impossible to fix sigma! Wrong fix sigma file or value!" << endl;
        }
      }

      if (enableRefl) {
        reflOverSgn = hMassForSgn[iPt]->Integral(hMassForSgn[iPt]->FindBin(massMin[iPt] * 1.0001), hMassForSgn[iPt]->FindBin(massMax[iPt] * 0.999));
        reflOverSgn = hMassForRefl[iPt]->Integral(hMassForRefl[iPt]->FindBin(massMin[iPt] * 1.0001), hMassForRefl[iPt]->FindBin(massMax[iPt] * 0.999)) / reflOverSgn;
        massFitter->setFixReflOverSgn(reflOverSgn);
        massFitter->setTemplateReflections(hMassRefl[iPt], HFInvMassFitter::DoubleGaus);
      }

      massFitter->doFit(false);

      double rawYield = massFitter->getRawYield();
      double rawYieldErr = massFitter->getRawYieldError();
      double sigma = massFitter->getSigma();
      double sigmaErr = massFitter->getSigmaUncertainty();
      double mean = massFitter->getMean();
      double meanErr = massFitter->getMeanUncertainty();
      double reducedChiSquare = massFitter->getChiSquareOverNDF();
      double significance = massFitter->getSignificance();
      double significanceErr = massFitter->getSignificanceError();
      double bkg = massFitter->getBkgYield();
      double bkgErr = massFitter->getBkgYieldError();

      hRawYields->SetBinContent(iPt + 1, rawYield);
      hRawYields->SetBinError(iPt + 1, rawYieldErr);
      hRawYieldsSigma->SetBinContent(iPt + 1, sigma);
      hRawYieldsSigma->SetBinError(iPt + 1, sigmaErr);
      hRawYieldsMean->SetBinContent(iPt + 1, mean);
      hRawYieldsMean->SetBinError(iPt + 1, meanErr);
      hRawYieldsSignificance->SetBinContent(iPt + 1, significance);
      hRawYieldsSignificance->SetBinError(iPt + 1, significanceErr);
      hRawYieldsSgnOverBkg->SetBinContent(iPt + 1, rawYield / bkg);
      hRawYieldsSgnOverBkg->SetBinError(iPt + 1, rawYield / bkg * std::sqrt(rawYieldErr / rawYield * rawYieldErr / rawYield + bkgErr / bkg * bkgErr / bkg));
      hRawYieldsSignal->SetBinContent(iPt + 1, rawYield);
      hRawYieldsSignal->SetBinError(iPt + 1, rawYieldErr);
      hRawYieldsBkg->SetBinContent(iPt + 1, bkg);
      hRawYieldsBkg->SetBinError(iPt + 1, bkgErr);
      hRawYieldsChiSquare->SetBinContent(iPt + 1, reducedChiSquare);
      hRawYieldsChiSquare->SetBinError(iPt + 1, 1.e-20);
      if (enableRefl) {
        hReflectionOverSignal->SetBinContent(iPt + 1, reflOverSgn);
      }

      if (enableRefl) {
        if (nPtBins > 1) {
          canvasRefl[iCanvas]->cd(iPt - nCanvasesMax * iCanvas + 1);
        } else {
          canvasRefl[iCanvas]->cd();
        }
        massFitter->drawReflection(gPad);
        canvasRefl[iCanvas]->Modified();
        canvasRefl[iCanvas]->Update();
      }

      if (nPtBins > 1) {
        canvasMass[iCanvas]->cd(iPt - nCanvasesMax * iCanvas + 1);
      } else {
        canvasMass[iCanvas]->cd();
      }
      massFitter->drawFit(gPad);
      canvasMass[iCanvas]->Modified();
      canvasMass[iCanvas]->Update();

      if (nPtBins > 1) {
        canvasResiduals[iCanvas]->cd(iPt - nCanvasesMax * iCanvas + 1);
      } else {
        canvasResiduals[iCanvas]->cd();
      }
      massFitter->drawResidual(gPad);
      canvasResiduals[iCanvas]->Modified();
      canvasResiduals[iCanvas]->Update();
    }

    hFitConfig->SetBinContent(1, iPt + 1, massMin[iPt]);
    hFitConfig->SetBinContent(2, iPt + 1, massMax[iPt]);
    hFitConfig->SetBinContent(3, iPt + 1, nRebin[iPt]);
    if (fixSigma) {
      if (fixSigmaManual.empty()) {
        hFitConfig->SetBinContent(4, iPt + 1, hSigmaToFix->GetBinContent(iPt + 1));
      } else {
        hFitConfig->SetBinContent(4, iPt + 1, fixSigmaManual[iPt]);
      }
    }
    hFitConfig->SetBinContent(5, iPt + 1, bkgFuncConfig[iPt]);
    hFitConfig->SetBinContent(6, iPt + 1, sgnFuncConfig[iPt]);
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

  for (unsigned int iPt = 0; iPt < nPtBins; iPt++) {
    hMass[iPt]->Write();
  }
  hRawYields->Write();
  hRawYieldsSigma->Write();
  hRawYieldsMean->Write();
  hRawYieldsSignificance->Write();
  hRawYieldsSgnOverBkg->Write();
  hRawYieldsSignal->Write();
  hRawYieldsBkg->Write();
  hRawYieldsChiSquare->Write();
  hRawYieldsSigma2->Write();
  hRawYieldsFracGaus2->Write();
  hRawYieldsSecondPeak->Write();
  hRawYieldsMeanSecondPeak->Write();
  hRawYieldsSigmaSecondPeak->Write();
  hRawYieldsSignificanceSecondPeak->Write();
  hRawYieldsSigmaRatioSecondFirstPeak->Write();
  hRawYieldsSoverBSecondPeak->Write();
  hRawYieldsSignalSecondPeak->Write();
  hRawYieldsBkgSecondPeak->Write();
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

void setHistoStyle(TH1* histo, int color, double markerSize)
{
  histo->SetStats(kFALSE);
  histo->SetMarkerSize(markerSize);
  histo->SetMarkerStyle(20);
  histo->SetLineWidth(2);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
}

void divideCanvas(TCanvas* canvas, int nPtBins)
{
  if (nPtBins < 2) {
    canvas->cd();
  } else if (nPtBins == 2 || nPtBins == 3) {
    canvas->Divide(nPtBins, 1);
  } else if (nPtBins == 4 || nPtBins == 6 || nPtBins == 8) {
    canvas->Divide(nPtBins / 2, 2);
  } else if (nPtBins == 5 || nPtBins == 7) {
    canvas->Divide((nPtBins + 1) / 2, 2);
  } else if (nPtBins == 9 || nPtBins == 12 || nPtBins == 15) {
    canvas->Divide(nPtBins / 3, 3);
  } else if (nPtBins == 10 || nPtBins == 11) {
    canvas->Divide(4, 3);
  } else if (nPtBins == 13 || nPtBins == 14) {
    canvas->Divide(5, 3);
  } else if (nPtBins > 15 && nPtBins <= 20 && nPtBins % 4 == 0) {
    canvas->Divide(nPtBins / 4, 4);
  } else if (nPtBins > 15 && nPtBins <= 20 && nPtBins % 4 != 0) {
    canvas->Divide(5, 4);
  } else if (nPtBins == 21) {
    canvas->Divide(7, 3);
  } else if (nPtBins > 21 && nPtBins <= 25) {
    canvas->Divide(5, 5);
  } else if (nPtBins > 25 && nPtBins % 2 == 0) {
    canvas->Divide(nPtBins / 2, 2);
  } else {
    canvas->Divide((nPtBins + 1) / 2, 2);
  }
}
