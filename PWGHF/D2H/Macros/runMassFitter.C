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
#include <cstdio>
#include <functional>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#endif

using namespace rapidjson;

void runMassFitter(const std::string& configFileName = "config_massfitter.json");

TFile* openFileWithNullptrCheck(const std::string& fileName, const std::string& option = "read");

template <typename T>
T* getObjectWithNullPtrCheck(TFile* fileIn, const std::string& objectName);

template <typename T>
T readJsonField(const Document& config, const std::string& fieldName, const T& defaultValue);

template <typename T>
T readJsonField(const Document& config, const std::string& fieldName);

template <typename T>
void readJsonVector(std::vector<T>& vec, const Document& config, const std::string& fieldName, bool isRequired = false);

void readJsonVectorHistogram(std::vector<double>& vec, const Document& config, const std::string& fileNameFieldName, const std::string& histoNameFieldName);

void divideCanvas(TCanvas* c, int nSliceVarBins);

void setHistoStyle(TH1* histo, Color_t color = kBlack, Size_t markerSize = 1);

void runMassFitter(const std::string& configFileName)
{
  // load config
  FILE* configFile = fopen(configFileName.c_str(), "r");
  if (configFile == nullptr) {
    throw std::runtime_error("ERROR: Missing configuration json file: " + configFileName);
  }

  Document config;
  char readBuffer[65536];
  FileReadStream is(configFile, readBuffer, sizeof(readBuffer));
  config.ParseStream(is);
  fclose(configFile);

  bool const isMc = readJsonField<bool>(config, "IsMC");
  bool const writeSignalPar = readJsonField<bool>(config, "WriteSignalPar", true);
  std::string const particleName = readJsonField<std::string>(config, "Particle");
  std::string const collisionSystem = readJsonField<std::string>(config, "CollisionSystem", "");
  std::string const inputFileName = readJsonField<std::string>(config, "InFileName");
  std::string const reflFileName = readJsonField<std::string>(config, "ReflFileName", "");
  TString outputFileName = readJsonField<std::string>(config, "OutFileName", "mInvFit.root");

  std::vector<std::string> inputHistoName;
  std::vector<std::string> promptHistoName;
  std::vector<std::string> fdHistoName;
  std::vector<std::string> reflHistoName;
  std::vector<std::string> promptSecPeakHistoName;
  std::vector<std::string> fdSecPeakHistoName;
  std::vector<double> sliceVarMin;
  std::vector<double> sliceVarMax;
  std::vector<double> massMin;
  std::vector<double> massMax;
  std::vector<double> fixMeanManual;
  std::vector<double> fixSigmaManual;
  std::vector<double> fixSecondSigmaManual;
  std::vector<double> fixFracDoubleGausManual;
  std::vector<int> nRebin;
  std::vector<int> bkgFunc;
  std::vector<int> sgnFunc;

  readJsonVector(inputHistoName, config, "InputHistoName", true);
  readJsonVector(promptHistoName, config, "PromptHistoName");
  readJsonVector(fdHistoName, config, "FDHistoName");
  readJsonVector(reflHistoName, config, "ReflHistoName");
  readJsonVector(promptSecPeakHistoName, config, "PromptSecPeakHistoName");
  readJsonVector(fdSecPeakHistoName, config, "FDSecPeakHistoName");

  const bool fixMean = readJsonField<bool>(config, "FixMean", false);
  const std::string meanFile = readJsonField<std::string>(config, "MeanFile", "");
  readJsonVector(fixMeanManual, config, "FixMeanManual");

  const bool fixSigma = readJsonField<bool>(config, "FixSigma", false);
  const std::string sigmaFile = readJsonField<std::string>(config, "SigmaFile", "");
  readJsonVector(fixSigmaManual, config, "FixSigmaManual");

  const bool fixSecondSigma = readJsonField<bool>(config, "FixSecondSigma", false);
  const std::string secondSigmaFile = readJsonField<std::string>(config, "SecondSigmaFile", "");
  readJsonVector(fixSecondSigmaManual, config, "FixSecondSigmaManual");

  const bool fixFracDoubleGaus = readJsonField<bool>(config, "FixFracDoubleGaus", false);
  const std::string fracDoubleGausFile = readJsonField<std::string>(config, "FracDoubleGausFile", "");
  readJsonVector(fixFracDoubleGausManual, config, "FixFracDoubleGausManual");

  const TString sliceVarName = readJsonField<std::string>(config, "SliceVarName");
  const TString sliceVarUnit = readJsonField<std::string>(config, "SliceVarUnit");

  readJsonVector(sliceVarMin, config, "SliceVarMin", true);
  readJsonVector(sliceVarMax, config, "SliceVarMax", true);

  readJsonVector(massMin, config, "MassMin", true);
  readJsonVector(massMax, config, "MassMax", true);

  readJsonVector(nRebin, config, "Rebin", true);

  bool const includeSecPeak = readJsonField<bool>(config, "InclSecPeak", false);
  bool const useLikelihood = readJsonField<bool>(config, "UseLikelihood");

  readJsonVector(bkgFunc, config, "BkgFunc", true);
  readJsonVector(sgnFunc, config, "SgnFunc", true);

  const bool enableRefl = readJsonField<bool>(config, "EnableRefl", false);
  const bool drawBgPrefit = readJsonField<bool>(config, "drawBgPrefit", true);
  const bool highlightPeakRegion = readJsonField<bool>(config, "highlightPeakRegion", true);

  const int nSliceVarBins = static_cast<int>(sliceVarMin.size());
  std::vector<double> sliceVarLimits(nSliceVarBins + 1);

  auto checkVectorSize = [&](const auto& vec, const std::string& name = "", const bool isEmptyOk = false) {
    if (vec.size() != static_cast<size_t>(nSliceVarBins)) {
      if (isEmptyOk && vec.empty()) {
        return;
      }
      throw std::runtime_error("ERROR: inconsistent vector size for " + name + "! Exit");
    }
  };

  checkVectorSize(inputHistoName, "inputHistoName", true);
  checkVectorSize(promptHistoName, "promptHistoName", true);
  checkVectorSize(fdHistoName, "fdHistoName", true);
  checkVectorSize(reflHistoName, "reflHistoName", true);
  checkVectorSize(promptSecPeakHistoName, "promptSecPeakHistoName", true);
  checkVectorSize(fdSecPeakHistoName, "fdSecPeakHistoName", true);
  checkVectorSize(sliceVarMin, "sliceVarMin");
  checkVectorSize(sliceVarMax, "sliceVarMax");
  checkVectorSize(massMin, "massMin");
  checkVectorSize(massMax, "massMax");
  checkVectorSize(fixMeanManual, "fixMeanManual", true);
  checkVectorSize(fixSigmaManual, "fixSigmaManual", true);
  checkVectorSize(fixSecondSigmaManual, "fixSecondSigmaManual", true);
  checkVectorSize(fixFracDoubleGausManual, "fixFracDoubleGausManual", true);
  checkVectorSize(nRebin, "nRebin");
  checkVectorSize(bkgFunc, "bkgFunc");
  checkVectorSize(sgnFunc, "sgnFunc");

  for (int iSliceVar = 0; iSliceVar < nSliceVarBins; iSliceVar++) {
    sliceVarLimits[iSliceVar] = sliceVarMin[iSliceVar];

    if (bkgFunc[iSliceVar] < 0 || bkgFunc[iSliceVar] >= HFInvMassFitter::NTypesOfBkgPdf) {
      throw std::runtime_error("ERROR: only Expo, Poly1, Poly2, Pow and PowEx background functions supported! Exit");
    }
    if (isMc && bkgFunc[iSliceVar] != HFInvMassFitter::NoBkg) {
      throw std::runtime_error("ERROR: in MC mode the background function must be NoBkg! Exit");
    }
    if (sgnFunc[iSliceVar] < 0 || sgnFunc[iSliceVar] >= HFInvMassFitter::NTypesOfSgnPdf) {
      throw std::runtime_error("ERROR: only SingleGaus, DoubleGaus and DoubleGausSigmaRatioPar signal functions supported! Exit");
    }
  }
  sliceVarLimits[nSliceVarBins] = sliceVarMax[nSliceVarBins - 1];

  struct decayInfo {
    std::string decayProducts;
    std::string pdgName;
    std::string decayFormulaLhs;
    std::string decayFormulaRhs;
  };

  std::map<std::string, decayInfo> particles{
    {"Dplus", {"K#pi#pi", "D+", "D^{+}", "K^{-}#pi^{+}#pi^{+}"}},
    {"D0", {"K#pi", "D0", "D^{0}", "K^{-}#pi^{+}"}},
    {"Ds", {"KK#pi", "D_s+", "D_{s}^{+}", "K^{-}K^{+}#pi^{+}"}},
    {"LcToPKPi", {"pK#pi", "Lambda_c+", "#Lambda_{c}^{+}", "pK^{-}#pi^{+}"}},
    {"LcToPK0s", {"pK^{0}_{s}", "Lambda_c+", "#Lambda_{c}^{+}", "pK^{0}_{s}"}},
    {"Dstar", {"D^{0}pi^{+}", "D*+", "D^{*+}", "D^{0}#pi^{+}"}},
    {"XicToXiPiPi", {"#Xi#pi#pi", "Xi_c+", "#Xi_{c}^{+}", "#Xi^{-}#pi^{+}#pi^{+}"}}};
  if (particles.find(particleName.c_str()) == particles.end()) {
    throw std::runtime_error("ERROR: only Dplus, D0, Ds, LcToPKPi, LcToPK0s, Dstar and XicToXiPiPi particles supported! Exit");
  }
  const auto& particle = particles[particleName.c_str()];
  const std::string massAxisTitle = "#it{M}(" + particle.decayProducts + ") (GeV/#it{c}^{2})";
  const double massPDG = TDatabasePDG::Instance()->GetParticle(particle.pdgName.c_str())->Mass();
  const std::vector<std::string> plotLabels = {(particle.decayFormulaLhs + " #rightarrow " + particle.decayFormulaRhs + " + c.c.").c_str(), collisionSystem};

  // load inv-mass histograms
  auto* inputFile = openFileWithNullptrCheck(inputFileName);

  TFile* inputFileRefl = nullptr;
  if (enableRefl) {
    inputFileRefl = openFileWithNullptrCheck(reflFileName);
  }

  std::vector<TH1*> hMassSgn(nSliceVarBins);
  std::vector<TH1*> hMassRefl(nSliceVarBins);
  std::vector<TH1*> hMass(nSliceVarBins);

  for (int iSliceVar = 0; iSliceVar < nSliceVarBins; iSliceVar++) {
    if (!isMc) {
      hMass[iSliceVar] = getObjectWithNullPtrCheck<TH1>(inputFile, inputHistoName[iSliceVar]);
      if (enableRefl) {
        hMassRefl[iSliceVar] = getObjectWithNullPtrCheck<TH1>(inputFileRefl, reflHistoName[iSliceVar]);
        hMassSgn[iSliceVar] = getObjectWithNullPtrCheck<TH1>(inputFileRefl, fdHistoName[iSliceVar]);
        hMassSgn[iSliceVar]->Add(getObjectWithNullPtrCheck<TH1>(inputFileRefl, promptHistoName[iSliceVar]));
      }
    } else {
      hMass[iSliceVar] = getObjectWithNullPtrCheck<TH1>(inputFile, promptHistoName[iSliceVar]);
      hMass[iSliceVar]->Add(getObjectWithNullPtrCheck<TH1>(inputFile, fdHistoName[iSliceVar]));
      if (includeSecPeak) {
        hMass[iSliceVar]->Add(getObjectWithNullPtrCheck<TH1>(inputFile, promptSecPeakHistoName[iSliceVar]));
        hMass[iSliceVar]->Add(getObjectWithNullPtrCheck<TH1>(inputFile, fdSecPeakHistoName[iSliceVar]));
      }
    }
    hMass[iSliceVar]->SetDirectory(nullptr);
  }
  inputFile->Close();
  if (enableRefl) {
    inputFileRefl->Close();
  }

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

  enum {
    ConfigMassMin = 1,
    ConfigMassMax,
    ConfigNRebin,
    ConfigFixSigma,
    ConfigBkgFunc,
    ConfigSgnFunc,
    NConfigsToSave
  };
  auto* hFitConfig = new TH2F("hfitConfig", "Fit Configurations", NConfigsToSave - 1, 0, NConfigsToSave - 1, nSliceVarBins, sliceVarLimits.data());
  const char* hFitConfigXLabel[NConfigsToSave - 1] = {"mass min", "mass max", "rebin num", "fix sigma", "bkg func", "sgn func"};
  hFitConfig->SetStats(false);
  for (int i = 0; i < NConfigsToSave - 1; i++) {
    hFitConfig->GetXaxis()->SetBinLabel(i + 1, hFitConfigXLabel[i]);
  }
  hFitConfig->LabelsDeflate("X");
  hFitConfig->LabelsDeflate("Y");
  hFitConfig->LabelsOption("v");

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
        auto* fixInputFile = openFileWithNullptrCheck(fixFileName);
        const std::string histName = "hRawYields" + var;
        histToFix = getObjectWithNullPtrCheck<TH1>(fixInputFile, histName);
        histToFix->SetDirectory(nullptr);
        if (histToFix->GetNbinsX() != nSliceVarBins) {
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

  int canvasSize[2] = {1920, 1080};
  if (nSliceVarBins == 1) {
    canvasSize[0] = 500;
    canvasSize[1] = 500;
  }

  int constexpr nCanvasesMax = 20; // do not put more than 20 bins per canvas to make them visible
  const int nCanvases = std::ceil(static_cast<float>(nSliceVarBins) / nCanvasesMax);
  std::vector<TCanvas*> canvasMass(nCanvases);
  std::vector<TCanvas*> canvasResiduals(nCanvases);
  std::vector<TCanvas*> canvasRatio(nCanvases);
  std::vector<TCanvas*> canvasRefl(nCanvases);
  for (int iCanvas = 0; iCanvas < nCanvases; iCanvas++) {
    const int nPads = (nCanvases == 1) ? nSliceVarBins : nCanvasesMax;
    canvasMass[iCanvas] = new TCanvas(Form("canvasMass%d", iCanvas), Form("canvasMass%d", iCanvas), canvasSize[0], canvasSize[1]);
    divideCanvas(canvasMass[iCanvas], nPads);

    canvasResiduals[iCanvas] = new TCanvas(Form("canvasResiduals%d", iCanvas), Form("canvasResiduals%d", iCanvas), canvasSize[0], canvasSize[1]);
    divideCanvas(canvasResiduals[iCanvas], nPads);

    canvasRatio[iCanvas] = new TCanvas(Form("canvasRatio%d", iCanvas), Form("canvasRatio%d", iCanvas), canvasSize[0], canvasSize[1]);
    divideCanvas(canvasRatio[iCanvas], nPads);

    if (enableRefl) {
      canvasRefl[iCanvas] = new TCanvas(Form("canvasRefl%d", iCanvas), Form("canvasRefl%d", iCanvas), canvasSize[0], canvasSize[1]);
      divideCanvas(canvasRefl[iCanvas], nPads);
    }
  }

  for (int iSliceVar = 0; iSliceVar < nSliceVarBins; iSliceVar++) {
    const int iCanvas = std::floor(static_cast<float>(iSliceVar) / nCanvasesMax);

    hMassForFit[iSliceVar] = hMass[iSliceVar]->Rebin(nRebin[iSliceVar]);
    TString const ptTitle =
      Form("%0.2f < " + sliceVarName + " < %0.2f " + sliceVarUnit, sliceVarMin[iSliceVar], sliceVarMax[iSliceVar]);
    hMassForFit[iSliceVar]->SetTitle(Form("%s;%s;Counts per %0.1f MeV/#it{c}^{2}",
                                          ptTitle.Data(), massAxisTitle.c_str(),
                                          hMassForFit[iSliceVar]->GetBinWidth(1) * 1000));
    hMassForFit[iSliceVar]->SetName(Form("MassForFit%d", iSliceVar));

    if (enableRefl) {
      hMassForRefl[iSliceVar] = hMassRefl[iSliceVar]->Rebin(nRebin[iSliceVar]);
      hMassForSgn[iSliceVar] = hMassSgn[iSliceVar]->Rebin(nRebin[iSliceVar]);
    }

    double reflOverSgn = 0;

    HFInvMassFitter* massFitter = new HFInvMassFitter(hMassForFit[iSliceVar], massMin[iSliceVar], massMax[iSliceVar], bkgFunc[iSliceVar], sgnFunc[iSliceVar]);
    massFitter->setDrawBgPrefit(drawBgPrefit);
    massFitter->setHighlightPeakRegion(highlightPeakRegion);
    massFitter->setInitialGaussianMean(massPDG);
    massFitter->setParticlePdgMass(massPDG);
    massFitter->setBoundGaussianMean(massPDG, 0.8 * massPDG, 1.2 * massPDG);
    if (useLikelihood) {
      massFitter->setUseLikelihoodFit();
    } else {
      massFitter->setUseChi2Fit();
    }

    auto setFixedValue = [&iSliceVar](bool const& isFix, std::vector<double> const& fixManual, const TH1* histToFix, std::function<void(double)> setFunc, std::string const& var) -> void {
      if (isFix) {
        if (fixManual.empty() && histToFix == nullptr) {
          throw std::runtime_error("Histogram to fix " + var + " is null while isFix==true and fixManual is empty");
        }
        const auto valueToFix = fixManual.empty() ? histToFix->GetBinContent(iSliceVar + 1) : fixManual[iSliceVar];
        setFunc(valueToFix);
        printf("*****************************\n");
        printf("FIXED %s: %f\n", var.data(), valueToFix);
        printf("*****************************\n");
      }
    };

    setFixedValue(fixMean, fixMeanManual, hMeanToFix, std::bind(&HFInvMassFitter::setFixGaussianMean, massFitter, std::placeholders::_1), "MEAN");
    setFixedValue(fixSigma, fixSigmaManual, hSigmaToFix, std::bind(&HFInvMassFitter::setFixGaussianSigma, massFitter, std::placeholders::_1), "SIGMA");
    setFixedValue(fixSecondSigma, fixSecondSigmaManual, hSecondSigmaToFix, std::bind(&HFInvMassFitter::setFixSecondGaussianSigma, massFitter, std::placeholders::_1), "SECOND SIGMA");
    setFixedValue(fixFracDoubleGaus, fixFracDoubleGausManual, hFracDoubleGausToFix, std::bind(&HFInvMassFitter::setFixFrac2Gaus, massFitter, std::placeholders::_1), "FRAC DOUBLE GAUS");

    if (!isMc && enableRefl) {
      reflOverSgn = hMassForSgn[iSliceVar]->Integral(hMassForSgn[iSliceVar]->FindBin(massMin[iSliceVar] * 1.0001), hMassForSgn[iSliceVar]->FindBin(massMax[iSliceVar] * 0.999));
      reflOverSgn = hMassForRefl[iSliceVar]->Integral(hMassForRefl[iSliceVar]->FindBin(massMin[iSliceVar] * 1.0001), hMassForRefl[iSliceVar]->FindBin(massMax[iSliceVar] * 0.999)) / reflOverSgn;
      massFitter->setFixReflOverSgn(reflOverSgn);
      massFitter->setTemplateReflections(hMassRefl[iSliceVar]);
    }

    massFitter->doFit();

    auto drawOnCanvas = [&](std::vector<TCanvas*>& canvas, std::function<void()> drawer) {
      if (nSliceVarBins > 1) {
        canvas[iCanvas]->cd(iSliceVar - nCanvasesMax * iCanvas + 1);
      } else {
        canvas[iCanvas]->cd();
      }
      drawer();
      canvas[iCanvas]->Modified();
      canvas[iCanvas]->Update();
    };

    drawOnCanvas(canvasMass, [&]() { massFitter->drawFit(gPad, plotLabels, writeSignalPar); });
    drawOnCanvas(canvasRatio, [&]() { massFitter->drawRatio(gPad); });
    if (bkgFunc[iSliceVar] != HFInvMassFitter::NoBkg) {
      drawOnCanvas(canvasResiduals, [&]() { massFitter->drawResidual(gPad); });
    }
    if (enableRefl) {
      drawOnCanvas(canvasRefl, [&]() { massFitter->drawReflection(gPad); });
    }

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
    hReflectionOverSignal->SetBinContent(iSliceVar + 1, reflOverSgn);

    if (sgnFunc[iSliceVar] != HFInvMassFitter::SingleGaus) { // TODO foresee DSCB and Voigt cases
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

    hFitConfig->SetBinContent(ConfigMassMin, iSliceVar + 1, massMin[iSliceVar]);
    hFitConfig->SetBinContent(ConfigMassMax, iSliceVar + 1, massMax[iSliceVar]);
    hFitConfig->SetBinContent(ConfigNRebin, iSliceVar + 1, nRebin[iSliceVar]);
    if (fixSigma) {
      const auto valueToFix = fixSigmaManual.empty() ? hSigmaToFix->GetBinContent(iSliceVar + 1) : fixSigmaManual[iSliceVar];
      hFitConfig->SetBinContent(ConfigFixSigma, iSliceVar + 1, valueToFix);
    }
    hFitConfig->SetBinContent(ConfigBkgFunc, iSliceVar + 1, bkgFunc[iSliceVar]);
    hFitConfig->SetBinContent(ConfigSgnFunc, iSliceVar + 1, sgnFunc[iSliceVar]);
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

  for (int iSliceVar = 0; iSliceVar < nSliceVarBins; iSliceVar++) {
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
    const std::string printingBracket = nCanvases == 1 ? "" : iCanvas == 0             ? "("
                                                            : iCanvas == nCanvases - 1 ? ")"
                                                                                       : "";
    canvasMass[iCanvas]->Print(Form("%s%s", outputFileName.Data(), printingBracket.c_str()), "pdf");
    canvasRatio[iCanvas]->Print(Form("%s%s", outputFileRatio.Data(), printingBracket.c_str()), "pdf");
    if (!isMc) {
      canvasResiduals[iCanvas]->Print(Form("%s%s", outputFileNameResidual.Data(), printingBracket.c_str()), "pdf");
    }
  }
}

void setHistoStyle(TH1* histo, Color_t color, Size_t markerSize)
{
  histo->SetStats(false);
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

TFile* openFileWithNullptrCheck(const std::string& fileName, const std::string& option)
{
  TFile* file = TFile::Open(fileName.c_str(), option.c_str());
  if (file == nullptr || file->IsZombie()) {
    throw std::runtime_error("openFileWithNullptrCheck(): Cannot open file " + fileName);
  }
  return file;
}

template <typename T>
T* getObjectWithNullPtrCheck(TFile* fileIn, const std::string& objectName)
{
  T* ptr = fileIn->Get<T>(objectName.c_str());
  if (ptr == nullptr) {
    throw std::runtime_error("getObjectWithNullptrCheck() - object " + objectName + " in file " + fileIn->GetName() + " is missing");
  }
  return ptr;
}

template <typename T>
T getJsonValue(const Value& value)
{
  if constexpr (std::is_same_v<std::decay_t<T>, std::string>) {
    return value.GetString();
  } else if constexpr (std::is_same_v<std::decay_t<T>, bool>) {
    return value.GetBool();
  } else if constexpr (std::is_same_v<std::decay_t<T>, int>) {
    return value.GetInt();
  } else if constexpr (std::is_same_v<std::decay_t<T>, double>) {
    return value.GetDouble();
  }
  throw std::runtime_error("getJsonValue(): unsupported type!");
  return T();
}

template <typename T>
T readJsonField(const Document& config, const std::string& fieldName, const T& defaultValue)
{
  if (!config.HasMember(fieldName.c_str())) {
    return defaultValue;
  }
  const auto& value = config[fieldName.c_str()];
  return getJsonValue<T>(value);
}

template <typename T>
T readJsonField(const Document& config, const std::string& fieldName)
{
  if (!config.HasMember(fieldName.c_str())) {
    throw std::runtime_error("readJsonField(): missing field " + fieldName);
  }
  return readJsonField<T>(config, fieldName, T());
}

template <typename T>
void readJsonVector(std::vector<T>& vec, const Document& config, const std::string& fieldName, const bool isRequired)
{
  if (!vec.empty()) {
    throw std::runtime_error("readJsonVector(): vector is not empty!");
  }
  if (config.HasMember(fieldName.c_str())) {
    const Value& jsonArray = config[fieldName.c_str()];
    for (auto it = jsonArray.Begin(); it != jsonArray.End(); it++) {
      vec.push_back(getJsonValue<T>(*it));
    }
  } else if (isRequired) {
    throw std::runtime_error("readJsonVector(): missing required field " + fieldName);
  }
}

void readJsonVectorHistogram(std::vector<double>& vec, const Document& config, const std::string& fileNameFieldName, const std::string& histoNameFieldName)
{
  if (!vec.empty()) {
    throw std::runtime_error("readJsonVectorHistogram(): vector is not empty!");
  }
  const auto fileName = readJsonField<std::string>(config, fileNameFieldName);
  const auto histoName = readJsonField<std::string>(config, histoNameFieldName);
  if (fileName.empty() || histoName.empty()) {
    return;
  }
  TFile* inputFile = openFileWithNullptrCheck(fileName);
  TH1* histo = getObjectWithNullPtrCheck<TH1>(inputFile, histoName);
  for (int iBin = 1; iBin <= histo->GetNbinsX(); iBin++) {
    vec.push_back(histo->GetBinContent(iBin));
  }
  inputFile->Close();
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
