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

/// \file GetRawYields.C
/// \brief Macro for fitting charm nesons and invariant-mass spectra
///
/// \author Grazia Luparello <gluparel@cern.ch>

#if !defined(__CINT__) || defined(__CLING__)

#include <string>
#include <vector>

// if .h file not found, please include your local rapidjson/document.h and rapidjson/filereadstream.h here
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>

#include "TROOT.h"
#include "Riostream.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TDatabasePDG.h"
#include "TH2F.h"
#include "HFMassFitter.h"

#endif

using namespace std;
using namespace rapidjson;

enum { kDplus,
       kD0,
       kDs,
       kLcToPKPi,
       kLcToPK0s,
       kDstar };

//__________________________________________________________________________________________________________________
int GetRawYields(bool isMC = false, TString infilename = "AnalysisResults.root", TString inDirectory = "hf-correlator-ds-hadrons", TString invMassHistoName2D = "hMassDsVsPt", TString refFileName = "", TString cfgfilename = "config_massfitter.json", TString outFileName = "RawYieldsDsToPKPi.root", bool doSideband = false);
double SingleGaus(double* m, double* pars);
double DoublePeakSingleGaus(double* x, double* pars);
double DoubleGaus(double* m, double* pars);
double DoublePeakDoubleGaus(double* m, double* pars);
double DoubleGausSigmaRatio(double* m, double* pars);
void SetHistoStyle(TH1* histo, int color = kBlack, double markersize = 1.);
void SetStyle();
void DivideCanvas(TCanvas* c, int nPtBins);
TH1D* RebinHisto(TH1* hOrig, Int_t reb, Int_t firstUse = -1);

//__________________________________________________________________________________________________________________
template <typename ValueType>
void readArray(const Value& jsonArray, vector<ValueType>& output)
{
  for (auto it = jsonArray.Begin(); it != jsonArray.End(); it++) {
    auto value = it->template Get<ValueType>();
    output.emplace_back(value);
  }
}
//__________________________________________________________________________________________________________________
int GetRawYields(bool isMC, TString infilename, TString inDirectory, TString invMassHistoName2D, TString refFileName, TString cfgfilename, TString outFileName, bool doSideband)
{
  SetStyle();

  //=======================
  // LOAD CONFIGURATION FROM JSON FILE
  FILE* configFile = fopen(cfgfilename.Data(), "r");

  Document config;
  char readBuffer[65536];
  FileReadStream is(configFile, readBuffer, sizeof(readBuffer));
  config.ParseStream(is);
  fclose(configFile);

  TString ParticleName = config["Particle"].GetString();
  cout << "Particle to be analysed" << ParticleName << endl;

  int particle;
  if (ParticleName == "Dplus") {
    particle = kDplus;
  } else if (ParticleName == "Ds") {
    particle = kDs;
  } else if (ParticleName == "D0") {
    particle = kD0;
  } else if (ParticleName == "LcToPKPi") {
    particle = kLcToPKPi;
  } else if (ParticleName == "LcToPK0s") {
    particle = kLcToPK0s;
  } else if (ParticleName == "Dstar") {
    particle = kDstar;
  } else {
    cerr << "ERROR: only Dplus, Ds, D0, Dstar and Lc are supported! Exit";
    return -1;
  }

  vector<double> PtMin;
  vector<double> PtMax;
  vector<double> MassMin;
  vector<double> MassMax;
  vector<double> Rebin;
  vector<int> bkgFuncConfig;
  vector<int> sgnFuncConfig;
  vector<int> fixSigma;
  vector<double> fixSigmaManually;
  vector<int> fixMean;
  vector<double> fixMeanManually;
  vector<int> InclSecPeak;     // For Ds only
  vector<double> SigmaSecPeak; // For Ds only
  vector<double> MinSBLeft;
  vector<double> MaxSBLeft;
  vector<double> MinSBRight;
  vector<double> MaxSBRight;

  const Value& ptMinValue = config["PtMin"];
  readArray(ptMinValue, PtMin);

  cout << "Number of PtBins for fit = " << PtMin.size() << endl;
  const unsigned int nPtBins = PtMin.size();

  const Value& ptMaxValue = config["PtMax"];
  readArray(ptMaxValue, PtMax);
  if (PtMax.size() != nPtBins) {
    cout << "ERROR: size of the vector PtMax is different from the size of vector PtMin" << endl;
    return -1;
  }

  const Value& massMinValue = config["MassMin"];
  readArray(massMinValue, MassMin);
  if (MassMin.size() != nPtBins) {
    cout << "ERROR: size of the vector MassMin is different from the number of PtBins" << endl;
    return -1;
  }

  const Value& massMaxValue = config["MassMax"];
  readArray(massMaxValue, MassMax);
  if (MassMax.size() != nPtBins) {
    cout << "ERROR: size of the vector MassMax is different from the number of PtBins" << endl;
    return -1;
  }

  const Value& optFixSigma = config["FixSigma"];
  readArray(optFixSigma, fixSigma);
  if (fixSigma.size() == 1)
    fixSigma.assign(PtMin.size(), fixSigma[0]);
  if (fixSigma.size() != nPtBins) {
    cout << "ERROR: size of the vector fixSigma is different from the number of PtBins. You have to correct the length of the vector in the config file or use just one entry" << endl;
    return -1;
  }

  const Value& fixSigmaManual = config["FixSigmaManual"];
  readArray(fixSigmaManual, fixSigmaManually);
  if (fixSigmaManually.size() != 0 && fixSigmaManually.size() != nPtBins) {
    cout << "ERROR: size of the vector fixSigmaManual is different from the number of PtBins" << endl;
    return -1;
  }

  string infilenameSigma = config["SigmaFile"].GetString();
  bool fixSigmaRatio = config["FixSigmaRatio"].GetBool(); // used only if SgnFunc=k2GausSigmaRatioPar
  string infilenameSigmaRatio = config["SigmaRatioFile"].GetString();

  bool isSigmaMultFromUnc = false;
  double sigmaMult = 1.;
  double sigmaMultFromConfig = config["SigmaMultFactor"].GetDouble();
  if (sigmaMultFromConfig != 0)
    sigmaMult = sigmaMultFromConfig;

  string sigmaMultFromUnc = config["SigmaMultFactorUnc"].GetString();
  if (sigmaMultFromUnc != "" && (sigmaMultFromUnc == "MinusUnc" || sigmaMultFromUnc == "PlusUnc")) {
    cout << sigmaMultFromUnc << endl;
    isSigmaMultFromUnc = true;
  }

  const Value& optFixMean = config["FixMean"];
  readArray(optFixMean, fixMean);
  if (fixMean.size() == 1)
    fixMean.assign(PtMin.size(), fixMean[0]);
  if (fixMean.size() != nPtBins) {
    cout << "ERROR: size of the vector fixMean is different from the number of PtBins. You have to correct the length of the vector in the config file or use just one entry" << endl;
    return -1;
  }

  const Value& fixMeanManual = config["FixMeanManual"];
  readArray(fixMeanManual, fixMeanManually);
  if (fixMeanManually.size() != 0 && fixMeanManually.size() != nPtBins) {
    cout << "ERROR: size of the vector fixMean is different from the number of PtBins" << endl;
    return -1;
  }

  string meanFile = config["MeanFile"].GetString();
  bool boundMean = config["BoundMean"].GetBool();

  bool enableRefl = config["EnableRefl"].GetBool();
  bool UseLikelihood = config["UseLikelihood"].GetBool();

  const Value& rebinValue = config["Rebin"];
  readArray(rebinValue, Rebin);
  if (Rebin.size() != nPtBins) {
    cout << "ERROR: size of the vector Rebin is different from the number of PtBins" << endl;
    return -1;
  }

  const Value& optInclSecPeak = config["InclSecPeak"];
  readArray(optInclSecPeak, InclSecPeak);
  if (InclSecPeak.size() == 1)
    InclSecPeak.assign(PtMin.size(), InclSecPeak[0]);
  if (InclSecPeak.size() != nPtBins) {
    cout << "ERROR: size of the vector InclSecPeak is different from the number of PtBins. You have to correct the length of the vector in the config file or use just one entry" << endl;
    return -1;
  }
  string infilenameSigmaSecPeak = config["SigmaFileSecPeak"].GetString();

  const Value& sigmaSecPeak = config["SigmaSecPeak"];
  readArray(sigmaSecPeak, SigmaSecPeak);
  if (SigmaSecPeak.size() != 0 && SigmaSecPeak.size() != nPtBins) {
    cout << "ERROR: size of the vector MassMax is different from the number of PtBins" << endl;
    return -1;
  }

  double sigmaMultSecPeak = config["SigmaMultFactorSecPeak"].GetDouble();
  bool fixSigmaToFirstPeak = config["FixSigmaToFirstPeak"].GetBool();

  if (doSideband) {
    const Value& minSBL = config["MinSidebandLeft"];
    readArray(minSBL, MinSBLeft);
    if (MinSBLeft.size() != nPtBins) {
      cout << "ERROR: size of the vector MinSBLeft is different from the number of PtBins" << endl;
      return -1;
    }

    const Value& maxSBL = config["MaxSidebandLeft"];
    readArray(maxSBL, MaxSBLeft);
    if (MaxSBLeft.size() != nPtBins) {
      cout << "ERROR: size of the vector MaxSBLeft is different from the number of PtBins" << endl;
      return -1;
    }

    const Value& minSBR = config["MinSidebandRight"];
    readArray(minSBR, MinSBRight);
    if (MinSBRight.size() != nPtBins) {
      cout << "ERROR: size of the vector MinSBRight is different from the number of PtBins" << endl;
      return -1;
    }

    const Value& maxSBR = config["MaxSidebandRight"];
    readArray(maxSBR, MaxSBRight);
    if (MaxSBRight.size() != nPtBins) {
      cout << "ERROR: size of the vector MaxSBRight is different from the number of PtBins" << endl;
      return -1;
    }
  }

  const Value& bkgFuncValue = config["BkgFunc"];
  readArray(bkgFuncValue, bkgFuncConfig);
  if (bkgFuncConfig.size() != nPtBins) {
    cout << "ERROR: size of the vector bkgFuncConfig is different from the number of PtBins" << endl;
    return -1;
  }

  const Value& sgnFuncValue = config["SgnFunc"];
  readArray(sgnFuncValue, sgnFuncConfig);
  if (sgnFuncConfig.size() != nPtBins) {
    cout << "ERROR: size of the vector signFuncConfig is different from the number of PtBins" << endl;
    return -1;
  }

  // const unsigned int nPtBins = PtMin.size();
  int BkgFunc[nPtBins], SgnFunc[nPtBins], degPol[nPtBins];
  double PtLims[nPtBins + 1];
  for (unsigned int iPt = 0; iPt < nPtBins; iPt++) {
    PtLims[iPt] = PtMin[iPt];
    PtLims[iPt + 1] = PtMax[iPt];

    degPol[iPt] = -1;
    if (bkgFuncConfig[iPt] == HFMassFitter::kExpo)
      BkgFunc[iPt] = HFMassFitter::kExpo;
    else if (bkgFuncConfig[iPt] == HFMassFitter::kLin)
      BkgFunc[iPt] = HFMassFitter::kLin;
    else if (bkgFuncConfig[iPt] == HFMassFitter::kPol2)
      BkgFunc[iPt] = HFMassFitter::kPol2;
    else if (bkgFuncConfig[iPt] == 6) // Pol3
    {
      BkgFunc[iPt] = 6;
      degPol[iPt] = 3;
      if (PtMin.size() > 1 && InclSecPeak[iPt] == 1) {
        cerr << "Pol3 and Pol4 fits work only with one bin if you have the secondary peak! Exit!" << endl;
        return -1;
      }
    } else if (bkgFuncConfig[iPt] == 7) // Pol3
    {
      BkgFunc[iPt] = 6;
      degPol[iPt] = 4;
      if (PtMin.size() > 1 && InclSecPeak[iPt] == 1) {
        cerr << "Pol3 and Pol4 fits work only with one bin if you have the secondary peak! Exit!" << endl;
        return -1;
      }
    } else if (bkgFuncConfig[iPt] == HFMassFitter::kPow)
      BkgFunc[iPt] = HFMassFitter::kPow;
    else if (bkgFuncConfig[iPt] == HFMassFitter::kPowEx)
      BkgFunc[iPt] = HFMassFitter::kPowEx;
    else {
      cerr << "ERROR: only kExpo, kLin, kPol2, kPol3, and kPol4, kPow and kPowEx background functions supported! Exit" << endl;
      return -1;
    }

    if (sgnFuncConfig[iPt] == HFMassFitter::kGaus)
      SgnFunc[iPt] = HFMassFitter::kGaus;
    else if (sgnFuncConfig[iPt] == HFMassFitter::k2Gaus)
      SgnFunc[iPt] = HFMassFitter::k2Gaus;
    else if (sgnFuncConfig[iPt] == HFMassFitter::k2GausSigmaRatioPar)
      SgnFunc[iPt] = HFMassFitter::k2GausSigmaRatioPar;
    else {
      cerr << "ERROR: only kGaus, k2Gaus and k2GausSigmaRatioPar signal functions supported! Exit" << endl;
      return -1;
    }
  }

  //=======================
  TString massaxistit = "";
  if (particle == kDplus)
    massaxistit = "#it{M}(K#pi#pi) (GeV/#it{c}^{2})";
  else if (particle == kD0)
    massaxistit = "#it{M}(K#pi) (GeV/#it{c}^{2})";
  else if (particle == kDs)
    massaxistit = "#it{M}(KK#pi) (GeV/#it{c}^{2})";
  else if (particle == kLcToPKPi)
    massaxistit = "#it{M}(pK#pi) (GeV/#it{c}^{2})";
  else if (particle == kLcToPK0s)
    massaxistit = "#it{M}(pK^{0}_{s}) (GeV/#it{c}^{2})";
  else if (particle == kDstar)
    massaxistit = "#it{M}(pi^{+}) (GeV/#it{c}^{2})";

  //=======================
  // Load Invariant Mass Histogram from train output
  auto infile = TFile::Open(infilename.Data());
  if (!infile || !infile->IsOpen())
    return -1;
  TFile* infileref = NULL;
  if (enableRefl) {
    infileref = TFile::Open(refFileName.Data());
    if (!infileref || !infileref->IsOpen())
      return -1;
  }
  TDirectoryFile* dir = (TDirectoryFile*)infile->Get(inDirectory.Data());
  TH2F* hMassVsPt = (TH2F*)dir->Get(invMassHistoName2D.Data());
  hMassVsPt->Draw();

  TH1F* hMassSig[nPtBins];
  TH1F* hMassRef[nPtBins];
  TH1F* hMass[nPtBins];
  TH1F* hEv = NULL; // TO BE ADDED

  TH1F* hPtAxis = (TH1F*)hMassVsPt->ProjectionY("hMassVsPt_py"); // prima _px
  for (unsigned int iPt = 0; iPt < nPtBins; iPt++) {
    Int_t binExtPtMin = hPtAxis->FindBin(PtMin[iPt] + 0.001);
    Int_t binExtPtMax = hPtAxis->FindBin(PtMax[iPt] - 0.001);
    hMass[iPt] = (TH1F*)hMassVsPt->ProjectionX(Form("hMass_binPt%d", iPt), binExtPtMin, binExtPtMax);
    if (enableRefl) {
      // TO BE CHECKED WITH REF FILES
      hMassRef[iPt] = static_cast<TH1F*>(infileref->Get(Form("hVarReflMass_%0.f_%0.f", PtMin[iPt] * 10, PtMax[iPt] * 10)));
      hMassSig[iPt] = static_cast<TH1F*>(infileref->Get(Form("hFDMass_%0.f_%0.f", PtMin[iPt] * 10, PtMax[iPt] * 10)));
      hMassSig[iPt]->Add(static_cast<TH1F*>(infileref->Get(Form("hPromptMass_%0.f_%0.f", PtMin[iPt] * 10, PtMax[iPt] * 10))));
    }
    hMass[iPt]->SetDirectory(0);
    SetHistoStyle(hMass[iPt]);
  }
  // infile ->Close();

  //=======================
  // DEFINE OUTPUT HISTOS
  auto hRawYields = new TH1D("hRawYields", ";#it{p}_{T} (GeV/#it{c});raw yield", nPtBins, PtLims);
  auto hRawYieldsSigma = new TH1D("hRawYieldsSigma", ";#it{p}_{T} (GeV/#it{c});width (GeV/#it{c}^{2})", nPtBins, PtLims);
  auto hRawYieldsSigmaRatio = new TH1D("hRawYieldsSigmaRatio", ";#it{p}_{T} (GeV/#it{c});ratio #sigma_{1}/#sigma_{2}", nPtBins, PtLims);
  auto hRawYieldsSigma2 = new TH1D("hRawYieldsSigma2", ";#it{p}_{T} (GeV/#it{c});width (GeV/#it{c}^{2})", nPtBins, PtLims);
  auto hRawYieldsMean = new TH1D("hRawYieldsMean", ";#it{p}_{T} (GeV/#it{c});mean (GeV/#it{c}^{2})", nPtBins, PtLims);
  auto hRawYieldsFracGaus2 = new TH1D("hRawYieldsFracGaus2", ";#it{p}_{T} (GeV/#it{c});second-gaussian fraction", nPtBins, PtLims);
  auto hRawYieldsSignificance = new TH1D("hRawYieldsSignificance", ";#it{p}_{T} (GeV/#it{c});significance (3#sigma)", nPtBins, PtLims);
  auto hRawYieldsSoverB = new TH1D("hRawYieldsSoverB", ";#it{p}_{T} (GeV/#it{c});S/B (3#sigma)", nPtBins, PtLims);
  auto hRawYieldsSignal = new TH1D("hRawYieldsSignal", ";#it{p}_{T} (GeV/#it{c});Signal (3#sigma)", nPtBins, PtLims);
  auto hRawYieldsBkg = new TH1D("hRawYieldsBkg", ";#it{p}_{T} (GeV/#it{c});Background (3#sigma)", nPtBins, PtLims);
  auto hRawYieldsChiSquare = new TH1D("hRawYieldsChiSquare", ";#it{p}_{T} (GeV/#it{c});#chi^{2}/#it{ndf}", nPtBins, PtLims);
  auto hRawYieldsSecondPeak = new TH1D("hRawYieldsSecondPeak", ";#it{p}_{T} (GeV/#it{c});raw yield second peak", nPtBins, PtLims);
  auto hRawYieldsMeanSecondPeak = new TH1D("hRawYieldsMeanSecondPeak", ";#it{p}_{T} (GeV/#it{c});mean second peak (GeV/#it{c}^{2})", nPtBins, PtLims);
  auto hRawYieldsSigmaSecondPeak = new TH1D("hRawYieldsSigmaSecondPeak", ";#it{p}_{T} (GeV/#it{c});width second peak (GeV/#it{c}^{2})", nPtBins, PtLims);
  auto hRawYieldsSignificanceSecondPeak = new TH1D("hRawYieldsSignificanceSecondPeak", ";#it{p}_{T} (GeV/#it{c});signficance second peak (3#sigma)", nPtBins, PtLims);
  auto hRawYieldsSigmaRatioSecondFirstPeak = new TH1D("hRawYieldsSigmaRatioSecondFirstPeak", ";#it{p}_{T} (GeV/#it{c});width second peak / width first peak", nPtBins, PtLims);
  auto hRawYieldsSoverBSecondPeak = new TH1D("hRawYieldsSoverBSecondPeak", ";#it{p}_{T} (GeV/#it{c});S/B second peak (3#sigma)", nPtBins, PtLims);
  auto hRawYieldsSignalSecondPeak = new TH1D("hRawYieldsSignalSecondPeak", ";#it{p}_{T} (GeV/#it{c});Signal second peak (3#sigma)", nPtBins, PtLims);
  auto hRawYieldsBkgSecondPeak = new TH1D("hRawYieldsBkgSecondPeak", ";#it{p}_{T} (GeV/#it{c});Background second peak (3#sigma)", nPtBins, PtLims);
  auto hRawYieldsTrue = new TH1D("hRawYieldsTrue", ";#it{p}_{T} (GeV/#it{c});true signal", nPtBins, PtLims);
  auto hRawYieldsSecondPeakTrue = new TH1D("hRawYieldsSecondPeakTrue", ";#it{p}_{T} (GeV/#it{c});true signal second peak", nPtBins, PtLims);
  auto hRelDiffRawYieldsFitTrue = new TH1D("hRelDiffRawYieldsFitTrue", ";#it{p}_{T} (GeV/#it{c}); (Y_{fit} - Y_{true}) / Y_{true}", nPtBins, PtLims);
  auto hRelDiffRawYieldsSecondPeakFitTrue = new TH1D("hRelDiffRawYieldsSecondPeakFitTrue", ";#it{p}_{T} (GeV/#it{c});(Y_{fit} - Y_{true}) / Y_{true} second peak", nPtBins, PtLims);
  auto hBackgroundSidebands = new TH1D("hBackgroundSideband", ";#it{p}_{T} (GeV/#it{c});Background sideband", nPtBins, PtLims);

  SetHistoStyle(hRawYields);
  SetHistoStyle(hRawYieldsSigma);
  SetHistoStyle(hRawYieldsSigma2);
  SetHistoStyle(hRawYieldsMean);
  SetHistoStyle(hRawYieldsFracGaus2);
  SetHistoStyle(hRawYieldsSignificance);
  SetHistoStyle(hRawYieldsSoverB);
  SetHistoStyle(hRawYieldsSignal);
  SetHistoStyle(hRawYieldsBkg);
  SetHistoStyle(hRawYieldsChiSquare);
  SetHistoStyle(hRawYieldsSecondPeak, kRed + 1);
  SetHistoStyle(hRawYieldsMeanSecondPeak, kRed + 1);
  SetHistoStyle(hRawYieldsSigmaSecondPeak, kRed + 1);
  SetHistoStyle(hRawYieldsSignificanceSecondPeak, kRed + 1);
  SetHistoStyle(hRawYieldsSigmaRatioSecondFirstPeak, kRed + 1);
  SetHistoStyle(hRawYieldsSoverBSecondPeak, kRed + 1);
  SetHistoStyle(hRawYieldsSignalSecondPeak, kRed + 1);
  SetHistoStyle(hRawYieldsBkgSecondPeak, kRed + 1);
  SetHistoStyle(hRawYieldsTrue);
  SetHistoStyle(hRawYieldsSecondPeakTrue, kRed + 1);
  SetHistoStyle(hRelDiffRawYieldsFitTrue);
  SetHistoStyle(hRelDiffRawYieldsSecondPeakFitTrue, kRed + 1);
  SetHistoStyle(hRawYieldsSigmaRatio, kRed + 1);
  SetHistoStyle(hBackgroundSidebands, kRed + 1);

  //=======================
  // READ HISTOS TO FIX VALUES - Note that also the possibility to manually define the fixed values in the config file
  TH1D* hSigmaToFix = NULL;
  if (accumulate(fixSigma.begin(), fixSigma.end(), 0) > 0 && infilenameSigma.compare("")) {
    auto infileSigma = TFile::Open(infilenameSigma.data());
    if (!infileSigma)
      return -2;
    hSigmaToFix = static_cast<TH1D*>(infileSigma->Get("hRawYieldsSigma"));
    hSigmaToFix->SetDirectory(0);
    if (static_cast<unsigned int>(hSigmaToFix->GetNbinsX()) != nPtBins)
      cout << "WARNING: Different number of bins for this analysis and histo for fix sigma" << endl;
    infileSigma->Close();
  }

  TH1D* hSigmaToFix1 = NULL;
  TH1D* hSigmaToFix2 = NULL;
  if (fixSigmaRatio) {
    std::cout << "Load sigma ratio from file " << infilenameSigmaRatio.data() << std::endl;
    auto infileSigmaRatio = TFile::Open(infilenameSigmaRatio.data());
    if (!infileSigmaRatio)
      return -2;
    hSigmaToFix1 = static_cast<TH1D*>(infileSigmaRatio->Get("hRawYieldsSigma"));
    hSigmaToFix1->SetDirectory(0);
    if (static_cast<unsigned int>(hSigmaToFix1->GetNbinsX()) != nPtBins)
      cout << "WARNING: Different number of bins for this analysis and histo for fix sigma ratio" << endl;

    hSigmaToFix2 = static_cast<TH1D*>(infileSigmaRatio->Get("hRawYieldsSigma2"));
    hSigmaToFix2->SetDirectory(0);
    if (static_cast<unsigned int>(hSigmaToFix2->GetNbinsX()) != nPtBins)
      cout << "WARNING: Different number of bins for this analysis and histo for fix sigma ratio" << endl;

    infileSigmaRatio->Close();
  }

  TH1D* hMeanToFix = NULL;
  if (accumulate(fixMean.begin(), fixMean.end(), 0) > 0 && meanFile.compare("")) {
    auto infileMean = TFile::Open(meanFile.data());
    if (!infileMean)
      return -3;
    hMeanToFix = static_cast<TH1D*>(infileMean->Get("hRawYieldsMean"));
    hMeanToFix->SetDirectory(0);
    if (static_cast<unsigned int>(hMeanToFix->GetNbinsX()) != nPtBins)
      cout << "WARNING: Different number of bins for this analysis and histo for fix mean" << endl;
    infileMean->Close();
  }

  TH1D* hSigmaFirstPeakMC = NULL;
  TH1D* hSigmaToFixSecPeak = NULL;

  /*
  //TO BE UNDERSTOOD
  bool InclSecPeakGlobal = false;
  for(int isp:InclSecPeak){
      if(isp != 0){
          InclSecPeakGlobal = true;
          break;
      }
  }

  TFile *infileSigmaSecPeak = NULL;
  if(InclSecPeakGlobal){
      infileSigmaSecPeak = TFile::Open(infilenameSigmaSecPeak.data());
  }
  if(!infileSigmaSecPeak && fixSigmaToFirstPeak)
      return -2;
  if(infileSigmaSecPeak) {
      hSigmaFirstPeakMC = static_cast<TH1D*>(infileSigmaSecPeak->Get("hRawYieldsSigma"));
      hSigmaToFixSecPeak = static_cast<TH1D*>(infileSigmaSecPeak->Get("hRawYieldsSigmaSecondPeak"));
      hSigmaFirstPeakMC->SetDirectory(0);
      hSigmaToFixSecPeak->SetDirectory(0);
      if(static_cast<unsigned int>(hSigmaFirstPeakMC->GetNbinsX()) != nPtBins ||
         static_cast<unsigned int>(hSigmaToFixSecPeak->GetNbinsX()) != nPtBins)
          cout << "WARNING: Different number of bins for this analysis and histos for fix sigma" << endl;
      infileSigmaSecPeak->Close();
}*/

  //=======================
  // FIT HISTOS
  double massDplus = TDatabasePDG::Instance()->GetParticle(411)->Mass();
  double massDs = TDatabasePDG::Instance()->GetParticle(431)->Mass();
  double massD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  double massLc = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
  double dmassDstar = TDatabasePDG::Instance()->GetParticle(413)->Mass() - TDatabasePDG::Instance()->GetParticle(421)->Mass();
  double massForFit = -1.;

  if (particle == kDplus)
    massForFit = massDplus;
  else if (particle == kDs)
    massForFit = massDs;
  else if (particle == kD0)
    massForFit = massD0;
  else if (particle == kLcToPKPi || particle == kLcToPK0s)
    massForFit = massLc;
  else if (particle == kDstar)
    massForFit = dmassDstar;

  TH1F* hMassForFit[nPtBins];
  TH1F* hMassForRel[nPtBins];
  TH1F* hMassForSig[nPtBins];

  //=======================
  // SET UP CANVAS
  int canvSize[2] = {1920, 1080};
  if (nPtBins == 1) {
    canvSize[0] = 500;
    canvSize[1] = 500;
  }

  int nMaxCanvases = 20; // do not put more than 20 bins per canvas to make them visible
  const int nCanvases = ceil((float)nPtBins / nMaxCanvases);
  TCanvas *cMass[nCanvases], *cResiduals[nCanvases];
  for (int iCanv = 0; iCanv < nCanvases; iCanv++) {
    int nPads = (nCanvases == 1) ? nPtBins : nMaxCanvases;
    cMass[iCanv] = new TCanvas(Form("cMass%d", iCanv), Form("cMass%d", iCanv), canvSize[0], canvSize[1]);
    DivideCanvas(cMass[iCanv], nPads);
    cResiduals[iCanv] = new TCanvas(Form("cResiduals%d", iCanv), Form("cResiduals%d", iCanv), canvSize[0], canvSize[1]);
    DivideCanvas(cResiduals[iCanv], nPads);
  }

  //=======================
  // LOOP on Pt BINS
  for (unsigned int iPt = 0; iPt < nPtBins; iPt++) {

    int iCanv = floor((float)iPt / nMaxCanvases);

    hMassForFit[iPt] = reinterpret_cast<TH1F*>(RebinHisto(hMass[iPt], Rebin[iPt]));
    TString pttitle = Form("%0.1f < #it{p}_{T} < %0.1f GeV/#it{c}", PtMin[iPt], PtMax[iPt]);
    hMassForFit[iPt]->SetTitle(Form("%s;%s;Counts per %0.f MeV/#it{c}^{2}", pttitle.Data(), massaxistit.Data(), hMassForFit[iPt]->GetBinWidth(1) * 1000));
    hMassForFit[iPt]->SetName(Form("MassForFit%d", iPt));
    if (enableRefl) {
      hMassForRel[iPt] = reinterpret_cast<TH1F*>(RebinHisto(hMassRef[iPt], Rebin[iPt]));
      hMassForSig[iPt] = reinterpret_cast<TH1F*>(RebinHisto(hMassSig[iPt], Rebin[iPt]));
    }
    double markerSize = 1.;
    if (nPtBins > 15)
      markerSize = 0.5;
    SetHistoStyle(hMassForFit[iPt], kBlack, markerSize);

    if (isMC) {                                                            // MC
      int parRawYield = 0, parMean = 1., parSigma1 = 2, parSigmaRatio = 3; // always the same
      int parSigma2 = -1, parFrac2Gaus = -1, parRawYieldSecPeak = -1, parMeanSecPeak = -1, parSigmaSecPeak = -1;
      TF1* massFunc = NULL;
      if (SgnFunc[iPt] == HFMassFitter::kGaus) {
        if (!(InclSecPeak[iPt] == 1 && particle == kDs)) {
          // if(!(InclSecPeak && particle==kDs)) {
          massFunc = new TF1(Form("massFunc%d", iPt), SingleGaus, MassMin[iPt], MassMax[iPt], 3);
          massFunc->SetParameters(hMassForFit[iPt]->Integral() * hMassForFit[iPt]->GetBinWidth(1), massForFit, 0.010);
        } else {
          massFunc = new TF1(Form("massFunc%d", iPt), DoublePeakSingleGaus, MassMin[iPt], MassMax[iPt], 6);
          massFunc->SetParameters(hMassForFit[iPt]->Integral() * hMassForFit[iPt]->GetBinWidth(1), massForFit, 0.010, hMassForFit[iPt]->Integral() * hMassForFit[iPt]->GetBinWidth(1), massDplus, 0.010);
          parRawYieldSecPeak = 3;
          parMeanSecPeak = 4;
          parSigmaSecPeak = 5;
        }
      } else if (SgnFunc[iPt] == HFMassFitter::k2Gaus) {
        parSigma2 = 3;
        parFrac2Gaus = 4;
        if (!(InclSecPeak[iPt] == 1 && particle == kDs)) {
          massFunc = new TF1(Form("massFunc%d", iPt), DoubleGaus, MassMin[iPt], MassMax[iPt], 5);
          massFunc->SetParameters(hMassForFit[iPt]->Integral() * hMassForFit[iPt]->GetBinWidth(1), massForFit, 0.010, 0.030, 0.9);
        } else {
          massFunc = new TF1(Form("massFunc%d", iPt), DoublePeakDoubleGaus, MassMin[iPt], MassMax[iPt], 8);
          massFunc->SetParameters(hMassForFit[iPt]->Integral() * hMassForFit[iPt]->GetBinWidth(1), massForFit, 0.010, 0.030, 0.9, hMassForFit[iPt]->Integral() * hMassForFit[iPt]->GetBinWidth(1), massDplus, 0.010);
          parRawYieldSecPeak = 5;
          parMeanSecPeak = 6;
          parSigmaSecPeak = 7;
        }
      } else if (SgnFunc[iPt] == HFMassFitter::k2GausSigmaRatioPar) {
        massFunc = new TF1(Form("massFunc%d", iPt), DoubleGausSigmaRatio, MassMin[iPt], MassMax[iPt], 5);
        massFunc->SetParameters(hMassForFit[iPt]->Integral() * hMassForFit[iPt]->GetBinWidth(1),
                                massForFit, 0.010, 2, 0.5);
        massFunc->SetParLimits(3, 1, 1e6);
        massFunc->SetParLimits(4, 0, 1);
      }

      if (nPtBins > 1)
        cMass[iCanv]->cd(iPt - nMaxCanvases * iCanv + 1);
      else
        cMass[iCanv]->cd();
      hMassForFit[iPt]->Fit(massFunc, "E"); // fit with chi2

      double rawyield = massFunc->GetParameter(parRawYield);
      double rawyielderr = massFunc->GetParError(parRawYield);
      double sigma = massFunc->GetParameter(parSigma1);
      double sigmaerr = massFunc->GetParError(parSigma1);
      double mean = massFunc->GetParameter(parMean);
      double meanerr = massFunc->GetParError(parMean);
      double redchi2 = massFunc->GetChisquare() / massFunc->GetNDF();

      hRawYields->SetBinContent(iPt + 1, rawyield);
      hRawYields->SetBinError(iPt + 1, rawyielderr);
      hRawYieldsSigma->SetBinContent(iPt + 1, sigma);
      hRawYieldsSigma->SetBinError(iPt + 1, sigmaerr);
      hRawYieldsMean->SetBinContent(iPt + 1, mean);
      hRawYieldsMean->SetBinError(iPt + 1, meanerr);
      hRawYieldsChiSquare->SetBinContent(iPt + 1, redchi2);
      hRawYieldsChiSquare->SetBinError(iPt + 1, 0.);

      hRawYieldsTrue->SetBinContent(iPt + 1, hMassForFit[iPt]->Integral());
      hRawYieldsTrue->SetBinError(iPt + 1, TMath::Sqrt(hMassForFit[iPt]->Integral()));
      hRelDiffRawYieldsFitTrue->SetBinContent(iPt + 1, rawyield - hMassForFit[iPt]->Integral());
      hRelDiffRawYieldsFitTrue->SetBinError(iPt + 1, TMath::Sqrt(rawyielderr * rawyielderr + hMassForFit[iPt]->Integral()));

      if (InclSecPeak[iPt] == 1 && particle == kDs) {
        double rawyieldSecPeak = massFunc->GetParameter(parRawYieldSecPeak);
        double rawyieldSecPeakerr = massFunc->GetParError(parRawYieldSecPeak);
        double sigmasecondpeak = massFunc->GetParameter(parSigmaSecPeak);
        double sigmasecondpeakerr = massFunc->GetParError(parSigmaSecPeak);
        double meansecondpeak = massFunc->GetParameter(parMeanSecPeak);
        double meansecondpeakerr = massFunc->GetParError(parMeanSecPeak);
        hRawYieldsSecondPeak->SetBinContent(iPt + 1, rawyieldSecPeak);
        hRawYieldsSecondPeak->SetBinError(iPt + 1, rawyieldSecPeakerr);
        hRawYieldsMeanSecondPeak->SetBinContent(iPt + 1, meansecondpeak);
        hRawYieldsMeanSecondPeak->SetBinError(iPt + 1, meansecondpeakerr);
        hRawYieldsSigmaSecondPeak->SetBinContent(iPt + 1, sigmasecondpeak);
        hRawYieldsSigmaSecondPeak->SetBinError(iPt + 1, sigmasecondpeakerr);
        hRawYieldsSigmaRatioSecondFirstPeak->SetBinContent(iPt + 1, sigmasecondpeak / sigma);
        hRawYieldsSigmaRatioSecondFirstPeak->SetBinError(iPt + 1, TMath::Sqrt(sigmaerr * sigmaerr / (sigma * sigma) + sigmasecondpeakerr * sigmasecondpeakerr / (sigmasecondpeak * sigmasecondpeak)) * sigmasecondpeak / sigma); // neglected correlation between parameters

        hRawYieldsSecondPeakTrue->SetBinContent(iPt + 1, rawyield);
        hRelDiffRawYieldsSecondPeakFitTrue->SetBinContent(iPt + 1, rawyield);
      }
      if (SgnFunc[iPt] == HFMassFitter::k2Gaus) {
        double sigma2 = massFunc->GetParameter(parSigma2);
        double sigma2err = massFunc->GetParError(parSigma2);
        double frac2gaus = massFunc->GetParameter(parFrac2Gaus);
        double frac2gauserr = massFunc->GetParError(parFrac2Gaus);
        hRawYieldsSigma2->SetBinContent(iPt + 1, sigma2);
        hRawYieldsSigma2->SetBinError(iPt + 1, sigma2err);
        hRawYieldsFracGaus2->SetBinContent(iPt + 1, frac2gaus);
        hRawYieldsFracGaus2->SetBinError(iPt + 1, frac2gauserr);
      } else if (SgnFunc[iPt] == HFMassFitter::k2GausSigmaRatioPar) {
        double sigmaRatio1_2 = massFunc->GetParameter(parSigmaRatio);
        double sigmaRatio1_2err = massFunc->GetParError(parSigmaRatio);
        double sigma2 = sigma / sigmaRatio1_2;
        double sigma2err = 0; // todo: add error propagation (correlation?)
        hRawYieldsSigmaRatio->SetBinContent(iPt + 1, sigmaRatio1_2);
        hRawYieldsSigmaRatio->SetBinError(iPt + 1, sigmaRatio1_2err);
        hRawYieldsSigma2->SetBinContent(iPt + 1, sigma2);
        hRawYieldsSigma2->SetBinError(iPt + 1, sigma2err);
      }
    } else { // data
      auto massFitter = new HFMassFitter(hMassForFit[iPt], MassMin[iPt], MassMax[iPt], BkgFunc[iPt], SgnFunc[iPt]);
      if (degPol[iPt] > 0)
        massFitter->SetPolDegreeForBackgroundFit(degPol[iPt]);
      if (UseLikelihood)
        massFitter->SetUseLikelihoodFit();
      if (fixMean[iPt]) {
        if (hMeanToFix != NULL)
          massFitter->SetFixGaussianMean(hMeanToFix->GetBinContent(iPt + 1));
        else
          massFitter->SetFixGaussianMean(fixMeanManually[iPt]);
      }
      if (boundMean)
        massFitter->SetBoundGaussianMean(massForFit, MassMin[iPt], MassMax[iPt]);
      else
        massFitter->SetInitialGaussianMean(massForFit);
      if (fixSigma[iPt]) {
        if (!isSigmaMultFromUnc)
          if (hSigmaToFix)
            massFitter->SetFixGaussianSigma(hSigmaToFix->GetBinContent(iPt + 1) * sigmaMult);
          else
            massFitter->SetFixGaussianSigma(fixSigmaManually[iPt]);
        else {
          if (sigmaMultFromUnc == "MinusUnc")
            massFitter->SetFixGaussianSigma(hSigmaToFix->GetBinContent(iPt + 1) - hSigmaToFix->GetBinError(iPt + 1));
          else if (sigmaMultFromUnc == "PlusUnc")
            massFitter->SetFixGaussianSigma(hSigmaToFix->GetBinContent(iPt + 1) + hSigmaToFix->GetBinError(iPt + 1));
          else
            cout << "WARNING: impossible to fix sigma! Wrong mult factor set in config file!" << endl;
        }
      } else {
        if (hSigmaToFix)
          massFitter->SetInitialGaussianSigma(hSigmaToFix->GetBinContent(iPt + 1) * sigmaMult);
        else if (particle == kDstar)
          massFitter->SetInitialGaussianSigma(0.001);
        else
          cout << "Initialise sigma" << endl;
        massFitter->SetInitialGaussianSigma(0.008);
      }
      if (fixSigmaRatio) {
        massFitter->SetFixRatio2GausSigma(
          hSigmaToFix1->GetBinContent(iPt + 1) / hSigmaToFix2->GetBinContent(iPt + 1));
      }
      if (InclSecPeak[iPt] && particle == kDs) {
        if (hSigmaToFixSecPeak) {
          massFitter->IncludeSecondGausPeak(massDplus, false, hSigmaToFixSecPeak->GetBinContent(iPt + 1) * sigmaMultSecPeak, true);
          if (fixSigmaToFirstPeak) {
            // fix D+ peak to sigmaMC(D+)/sigmaMC(Ds+)*sigmaData(Ds+)
            massFitter->MassFitter(false);
            double sigmaFirstPeak = massFitter->GetSigma();
            double sigmaRatioMC = hSigmaToFixSecPeak->GetBinContent(iPt + 1) / hSigmaFirstPeakMC->GetBinContent(iPt + 1);
            massFitter->IncludeSecondGausPeak(massDplus, false, sigmaRatioMC * sigmaFirstPeak, true);
          }
        } else {
          massFitter->IncludeSecondGausPeak(massDplus, false, SigmaSecPeak[iPt], true);
        }
      }

      if (enableRefl) {
        double rOverS = hMassForSig[iPt]->Integral(hMassForSig[iPt]->FindBin(MassMin[iPt] * 1.0001), hMassForSig[iPt]->FindBin(MassMax[iPt] * 0.999));
        rOverS = hMassForRel[iPt]->Integral(hMassForRel[iPt]->FindBin(MassMin[iPt] * 1.0001), hMassForRel[iPt]->FindBin(MassMax[iPt] * 0.999)) / rOverS;
        massFitter->SetFixReflOverS(rOverS);
        massFitter->SetTemplateReflections(hMassRef[iPt], "2gaus", MassMin[iPt], MassMax[iPt]);
      }

      massFitter->MassFitter(false);

      //=======================
      // RETRIEVE FITTER PARAMETERS
      double rawyield = massFitter->GetRawYield();
      double rawyielderr = massFitter->GetRawYieldError();
      double sigma = massFitter->GetSigma();
      double sigmaerr = massFitter->GetSigmaUncertainty();
      double mean = massFitter->GetMean();
      double meanerr = massFitter->GetMeanUncertainty();
      double redchi2 = massFitter->GetReducedChiSquare();
      double signif = 0., signiferr = 0.;
      double sgn = 0., sgnerr = 0.;
      double bkg = 0., bkgerr = 0.;
      massFitter->Significance(3, signif, signiferr);
      massFitter->Signal(3, sgn, sgnerr);
      massFitter->Background(3, bkg, bkgerr);

      hRawYields->SetBinContent(iPt + 1, rawyield);
      hRawYields->SetBinError(iPt + 1, rawyielderr);
      hRawYieldsSigma->SetBinContent(iPt + 1, sigma);
      hRawYieldsSigma->SetBinError(iPt + 1, sigmaerr);
      hRawYieldsMean->SetBinContent(iPt + 1, mean);
      hRawYieldsMean->SetBinError(iPt + 1, meanerr);
      hRawYieldsSignificance->SetBinContent(iPt + 1, signif);
      hRawYieldsSignificance->SetBinError(iPt + 1, signiferr);
      hRawYieldsSoverB->SetBinContent(iPt + 1, sgn / bkg);
      hRawYieldsSoverB->SetBinError(iPt + 1, sgn / bkg * TMath::Sqrt(sgnerr / sgn * sgnerr / sgn + bkgerr / bkg * bkgerr / bkg));
      hRawYieldsSignal->SetBinContent(iPt + 1, sgn);
      hRawYieldsSignal->SetBinError(iPt + 1, sgnerr);
      hRawYieldsBkg->SetBinContent(iPt + 1, bkg);
      hRawYieldsBkg->SetBinError(iPt + 1, bkgerr);
      hRawYieldsChiSquare->SetBinContent(iPt + 1, redchi2);
      hRawYieldsChiSquare->SetBinError(iPt + 1, 1.e-20);

      /*            for(int iS = 0; iS < nMassWindows; iS++)
            {
                massFitter->Significance(nSigma4SandB[iS],signif,signiferr);
                massFitter->Signal(nSigma4SandB[iS],sgn,sgnerr);
                massFitter->Background(nSigma4SandB[iS],bkg,bkgerr);

                hRawYieldsSignalDiffSigma[iS]->SetBinContent(iPt+1,sgn);
                hRawYieldsSignalDiffSigma[iS]->SetBinError(iPt+1,sgnerr);
                hRawYieldsBkgDiffSigma[iS]->SetBinContent(iPt+1,bkg);
                hRawYieldsBkgDiffSigma[iS]->SetBinError(iPt+1,bkgerr);
                hRawYieldsSoverBDiffSigma[iS]->SetBinContent(iPt+1,sgn/bkg);
                hRawYieldsSoverBDiffSigma[iS]->SetBinError(iPt+1,sgn/bkg*TMath::Sqrt(sgnerr/sgn*sgnerr/sgn+bkgerr/bkg*bkgerr/bkg));
                hRawYieldsSignifDiffSigma[iS]->SetBinContent(iPt+1,signif);
                hRawYieldsSignifDiffSigma[iS]->SetBinError(iPt+1,signiferr);
    }*/

      TF1* fTotFunc = massFitter->GetMassFunc();
      TF1* fBkgFunc = massFitter->GetBackgroundRecalcFunc();

      double parFrac2Gaus = -1, parsecondsigma = -1;
      if (SgnFunc[iPt] == HFMassFitter::k2Gaus) {
        if (!(InclSecPeak[iPt] == 1 && particle == kDs)) {
          parFrac2Gaus = fTotFunc->GetNpar() - 2;
          parsecondsigma = fTotFunc->GetNpar() - 1;
        } else {
          parFrac2Gaus = fTotFunc->GetNpar() - 5;
          parsecondsigma = fTotFunc->GetNpar() - 4;
        }

        double sigma2 = fTotFunc->GetParameter(parsecondsigma);
        double sigma2err = fTotFunc->GetParError(parsecondsigma);
        double frac2gaus = fTotFunc->GetParameter(parFrac2Gaus);
        double frac2gauserr = fTotFunc->GetParError(parFrac2Gaus);
        hRawYieldsSigma2->SetBinContent(iPt + 1, sigma2);
        hRawYieldsSigma2->SetBinError(iPt + 1, sigma2err);
        hRawYieldsFracGaus2->SetBinContent(iPt + 1, frac2gaus);
        hRawYieldsFracGaus2->SetBinError(iPt + 1, frac2gauserr);
      } else if (SgnFunc[iPt] == HFMassFitter::k2GausSigmaRatioPar) {
        int parRatioSigma = fTotFunc->GetNpar() - 1;
        double sigmaRatio1_2 = fTotFunc->GetParameter(parRatioSigma);
        double sigmaRatio1_2err = fTotFunc->GetParError(parRatioSigma);
        double sigma2 = sigma / sigmaRatio1_2;
        double sigma2err = 0; // todo: add error propagation (correlation?)
        hRawYieldsSigmaRatio->SetBinContent(iPt + 1, sigmaRatio1_2);
        hRawYieldsSigmaRatio->SetBinError(iPt + 1, sigmaRatio1_2err);
        hRawYieldsSigma2->SetBinContent(iPt + 1, sigma2);
        hRawYieldsSigma2->SetBinError(iPt + 1, sigma2err);
      }

      if (InclSecPeak[iPt] == 1 && particle == kDs) {
        int paryieldSecPeak = fTotFunc->GetNpar() - 3;
        int parMeansecondpeak = fTotFunc->GetNpar() - 2;
        int parSigmasecondpeak = fTotFunc->GetNpar() - 1;

        double rawyieldSecPeak = fTotFunc->GetParameter(paryieldSecPeak) / hMassForFit[iPt]->GetBinWidth(1);
        double rawyieldSecPeakerr = fTotFunc->GetParError(paryieldSecPeak) / hMassForFit[iPt]->GetBinWidth(1);
        double meansecondpeak = fTotFunc->GetParameter(parMeansecondpeak);
        double meansecondpeakerr = fTotFunc->GetParError(parMeansecondpeak);
        double sigmasecondpeak = fTotFunc->GetParameter(parSigmasecondpeak);
        double sigmasecondpeakerr = fTotFunc->GetParError(parSigmasecondpeak);

        double bkgSecPeak = fBkgFunc->Integral(meansecondpeak - 3 * sigmasecondpeak, meansecondpeak + 3 * sigmasecondpeak) / hMassForFit[iPt]->GetBinWidth(1);
        double bkgSecPeakerr = TMath::Sqrt(bkgSecPeak);
        double signalSecPeak = fTotFunc->Integral(meansecondpeak - 3 * sigmasecondpeak, meansecondpeak + 3 * sigmasecondpeak) / hMassForFit[iPt]->GetBinWidth(1) - bkgSecPeak;
        double signalSecPeakerr = TMath::Sqrt(signalSecPeak + bkgSecPeak);
        double signifSecPeak = -1., signifSecPeakerr = -1.;
        // AliVertexingHFUtils::ComputeSignificance(signalSecPeak,signalSecPeakerr,bkgSecPeak,bkgSecPeakerr,signifSecPeak,signifSecPeakerr);
        // TO BE FIXED. At the moment there is no method to calc this sigificance.

        hRawYieldsSecondPeak->SetBinContent(iPt + 1, rawyieldSecPeak);
        hRawYieldsSecondPeak->SetBinError(iPt + 1, rawyieldSecPeakerr);
        hRawYieldsMeanSecondPeak->SetBinContent(iPt + 1, meansecondpeak);
        hRawYieldsMeanSecondPeak->SetBinError(iPt + 1, meansecondpeakerr);
        hRawYieldsSigmaSecondPeak->SetBinContent(iPt + 1, sigmasecondpeak);
        hRawYieldsSigmaSecondPeak->SetBinError(iPt + 1, sigmasecondpeakerr);
        hRawYieldsSignificanceSecondPeak->SetBinContent(iPt + 1, signifSecPeak);
        hRawYieldsSignificanceSecondPeak->SetBinError(iPt + 1, signifSecPeakerr);
        hRawYieldsSigmaRatioSecondFirstPeak->SetBinContent(iPt + 1, sigmasecondpeak / sigma);
        hRawYieldsSigmaRatioSecondFirstPeak->SetBinError(iPt + 1, TMath::Sqrt(sigmaerr * sigmaerr / (sigma * sigma) + sigmasecondpeakerr * sigmasecondpeakerr / (sigmasecondpeak * sigmasecondpeak)) * sigmasecondpeak / sigma); // neglected correlation between parameters
        hRawYieldsSoverBSecondPeak->SetBinContent(iPt + 1, signalSecPeak / bkgSecPeak);
        hRawYieldsSoverBSecondPeak->SetBinError(iPt + 1, signalSecPeak / bkgSecPeak * TMath::Sqrt(signalSecPeakerr / signalSecPeak * signalSecPeakerr / signalSecPeak + bkgSecPeakerr / bkgSecPeak * bkgSecPeakerr / bkgSecPeak));
        hRawYieldsSignalSecondPeak->SetBinContent(iPt + 1, signalSecPeak);
        hRawYieldsSignalSecondPeak->SetBinError(iPt + 1, signalSecPeakerr);
        hRawYieldsBkgSecondPeak->SetBinContent(iPt + 1, bkgSecPeak);
        hRawYieldsBkgSecondPeak->SetBinError(iPt + 1, bkgSecPeakerr);
      }

      if (doSideband) {
        double bkgSB = fBkgFunc->Integral(MinSBLeft[iPt], MaxSBLeft[iPt]) / (Double_t)hMassForFit[iPt]->GetBinWidth(1) + fBkgFunc->Integral(MinSBRight[iPt], MaxSBRight[iPt]) / (Double_t)hMassForFit[iPt]->GetBinWidth(1);
        double bkgSBErr = TMath::Sqrt(bkgSB);
        // relative error evaluation: from histo

        int minLeft = hMassForFit[iPt]->FindBin(MinSBLeft[iPt]);
        int maxLeft = hMassForFit[iPt]->FindBin(MaxSBLeft[iPt]);
        int minRight = hMassForFit[iPt]->FindBin(MinSBRight[iPt]);
        int maxRight = hMassForFit[iPt]->FindBin(MaxSBRight[iPt]);

        double intB = hMassForFit[iPt]->Integral(minLeft, maxLeft) + hMassForFit[iPt]->Integral(minRight, maxRight);
        double sum2 = 0;
        for (int i = minLeft; i <= maxLeft; i++) {
          sum2 += hMassForFit[iPt]->GetBinError(i) * hMassForFit[iPt]->GetBinError(i);
        }
        for (int i = minRight; i <= maxRight; i++) {
          sum2 += hMassForFit[iPt]->GetBinError(i) * hMassForFit[iPt]->GetBinError(i);
        }
        double intBerr = TMath::Sqrt(sum2);
        double errbackground = intBerr / intB * bkgSB;
        cout << "Check integrals bkg sidebands: " << fBkgFunc->Integral(MinSBLeft[iPt], MaxSBLeft[iPt]) / (Double_t)hMassForFit[iPt]->GetBinWidth(1) << " " << fBkgFunc->Integral(MinSBRight[iPt], MaxSBRight[iPt]) / (Double_t)hMassForFit[iPt]->GetBinWidth(1) << " " << intB << endl;
        cout << "Check errors bkg sidebands: " << errbackground << "  " << bkgSBErr << endl;
        hBackgroundSidebands->SetBinContent(iPt + 1, bkgSB);
        hBackgroundSidebands->SetBinError(iPt + 1, bkgSBErr);
      }

      if (nPtBins > 1)
        cMass[iCanv]->cd(iPt - nMaxCanvases * iCanv + 1);
      else
        cMass[iCanv]->cd();

      hMassForFit[iPt]->GetYaxis()->SetRangeUser(hMassForFit[iPt]->GetMinimum() * 0.95, hMassForFit[iPt]->GetMaximum() * 1.2);
      massFitter->DrawHere(gPad);

      if (!isMC) {
        // residuals
        if (nPtBins > 1)
          cResiduals[iCanv]->cd(iPt - nMaxCanvases * iCanv + 1);
        else
          cResiduals[iCanv]->cd();
        massFitter->DrawHistoMinusFit(gPad);
      }
    }
    cMass[iCanv]->Modified();
    cMass[iCanv]->Update();
    cResiduals[iCanv]->Modified();
    cResiduals[iCanv]->Update();
  }

  // save output histos
  TFile outFile(outFileName.Data(), "recreate");
  for (int iCanv = 0; iCanv < nCanvases; iCanv++) {
    cMass[iCanv]->Write();
    if (!isMC)
      cResiduals[iCanv]->Write();
  }
  for (unsigned int iPt = 0; iPt < nPtBins; iPt++)
    hMass[iPt]->Write();
  hRawYields->Write();
  hRawYieldsSigma->Write();
  hRawYieldsMean->Write();
  hRawYieldsSignificance->Write();
  hRawYieldsSoverB->Write();
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
  hRawYieldsTrue->Write();
  hRawYieldsSecondPeakTrue->Write();
  hRelDiffRawYieldsFitTrue->Write();
  hRelDiffRawYieldsSecondPeakFitTrue->Write();
  hRawYieldsSigmaRatio->Write();
  hBackgroundSidebands->Write();
  // hEv->Write();
  if (!isMC) {
    /*TDirectoryFile dir("SandBDiffNsigma", "SandBDiffNsigma");
      dir.Write();
      dir.cd();
 for(int iS = 0; iS < nMassWindows; iS++)
      {
          hRawYieldsSignalDiffSigma[iS]->Write();
          hRawYieldsBkgDiffSigma[iS]->Write();
          hRawYieldsSoverBDiffSigma[iS]->Write();
          hRawYieldsSignifDiffSigma[iS]->Write();
    }
      dir.Close();*/
  }
  outFile.Close();

  outFileName.ReplaceAll(".root", ".pdf");
  TString outFileNameRes = outFileName;
  outFileNameRes.ReplaceAll(".pdf", "_Residuals.pdf");
  for (int iCanv = 0; iCanv < nCanvases; iCanv++) {
    if (iCanv == 0 && nCanvases > 1)
      cMass[iCanv]->SaveAs(Form("%s[", outFileName.Data()));
    cMass[iCanv]->SaveAs(outFileName.Data());
    if (iCanv == nCanvases - 1 && nCanvases > 1)
      cMass[iCanv]->SaveAs(Form("%s]", outFileName.Data()));

    if (!isMC) {
      if (iCanv == 0 && nCanvases > 1)
        cResiduals[iCanv]->SaveAs(Form("%s[", outFileNameRes.Data()));
      cResiduals[iCanv]->SaveAs(outFileNameRes.Data());
      if (iCanv == nCanvases - 1 && nCanvases > 1)
        cResiduals[iCanv]->SaveAs(Form("%s]", outFileNameRes.Data()));
    }
  }

  return 0;
}

//__________________________________________________________________________________________________________________
double SingleGaus(double* m, double* pars)
{
  double norm = pars[0], mean = pars[1], sigma = pars[2];

  return norm * TMath::Gaus(m[0], mean, sigma, true);
}

//__________________________________________________________________________________________________________________
double DoubleGaus(double* m, double* pars)
{
  double norm = pars[0], mean = pars[1], sigma1 = pars[2], sigma1_2 = pars[3], fg = pars[4];

  return norm * ((1 - fg) * TMath::Gaus(m[0], mean, sigma1, true) + fg * TMath::Gaus(m[0], mean, sigma1_2, true));
}

//__________________________________________________________________________________________________________________
double DoubleGausSigmaRatio(double* m, double* pars)
{
  double norm = pars[0], mean = pars[1], sigma1 = pars[2], sigma1OverSigma2 = pars[3], fg = pars[4];

  return norm * ((1 - fg) * TMath::Gaus(m[0], mean, sigma1, true) + fg * TMath::Gaus(m[0], mean, sigma1 / sigma1OverSigma2, true));
}

//__________________________________________________________________________________________________________________
double DoublePeakSingleGaus(double* m, double* pars)
{
  double norm1 = pars[0], mean1 = pars[1], sigma1 = pars[2]; // Ds peak
  double norm2 = pars[3], mean2 = pars[4], sigma2 = pars[5]; // Dplus peak

  return norm1 * TMath::Gaus(m[0], mean1, sigma1, true) + norm2 * TMath::Gaus(m[0], mean2, sigma2, true);
}

//__________________________________________________________________________________________________________________
double DoublePeakDoubleGaus(double* m, double* pars)
{
  double norm1 = pars[0], mean = pars[1], sigma1 = pars[2], sigma1_2 = pars[3], fg = pars[4]; // Ds peak
  double norm2 = pars[5], mean2 = pars[6], sigma2 = pars[7];                                  // Dplus peak

  return norm1 * ((1 - fg) * TMath::Gaus(m[0], mean, sigma1, true) + fg * TMath::Gaus(m[0], mean, sigma1_2, true)) + norm2 * TMath::Gaus(m[0], mean2, sigma2, true);
}

//__________________________________________________________________________________________________________________
void SetHistoStyle(TH1* histo, int color, double markersize)
{
  histo->SetStats(kFALSE);
  histo->SetMarkerSize(markersize);
  histo->SetMarkerStyle(20);
  histo->SetLineWidth(2);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
}

//__________________________________________________________________________________________________________________
void DivideCanvas(TCanvas* c, int nPtBins)
{
  if (nPtBins < 2)
    c->cd();
  else if (nPtBins == 2 || nPtBins == 3)
    c->Divide(nPtBins, 1);
  else if (nPtBins == 4 || nPtBins == 6 || nPtBins == 8)
    c->Divide(nPtBins / 2, 2);
  else if (nPtBins == 5 || nPtBins == 7)
    c->Divide((nPtBins + 1) / 2, 2);
  else if (nPtBins == 9 || nPtBins == 12 || nPtBins == 15)
    c->Divide(nPtBins / 3, 3);
  else if (nPtBins == 10 || nPtBins == 11)
    c->Divide(4, 3);
  else if (nPtBins == 13 || nPtBins == 14)
    c->Divide(5, 3);
  else if (nPtBins > 15 && nPtBins <= 20 && nPtBins % 4 == 0)
    c->Divide(nPtBins / 4, 4);
  else if (nPtBins > 15 && nPtBins <= 20 && nPtBins % 4 != 0)
    c->Divide(5, 4);
  else if (nPtBins == 21)
    c->Divide(7, 3);
  else if (nPtBins > 21 && nPtBins <= 25)
    c->Divide(5, 5);
  else if (nPtBins > 25 && nPtBins % 2 == 0)
    c->Divide(nPtBins / 2, 2);
  else
    c->Divide((nPtBins + 1) / 2, 2);
}

//__________________________________________________________________________________________________________________
void SetStyle()
{
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetTitleOffset(1.4, "y");
  gStyle->SetTitleOffset(1.2, "x");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TGaxis::SetMaxDigits(3);
}

//__________________________________________________________________________________________________________________
TH1D* RebinHisto(TH1* hOrig, Int_t reb, Int_t firstUse)
{
  /// Rebin histogram, from bin firstUse to lastUse
  /// Use all bins if firstUse=-1
  /// If ngroup is not an exact divider of the number of bins,
  ///  the bin width is kept as reb*original width
  ///  and the range of rebinned histogram is adapted

  Int_t nBinOrig = hOrig->GetNbinsX();
  Int_t firstBinOrig = 1;
  Int_t lastBinOrig = nBinOrig;
  Int_t nBinOrigUsed = nBinOrig;
  Int_t nBinFinal = nBinOrig / reb;
  if (firstUse >= 1) {
    firstBinOrig = firstUse;
    nBinFinal = (nBinOrig - firstUse + 1) / reb;
    nBinOrigUsed = nBinFinal * reb;
    lastBinOrig = firstBinOrig + nBinOrigUsed - 1;
  } else {
    Int_t exc = nBinOrigUsed % reb;
    if (exc != 0) {
      nBinOrigUsed -= exc;
      lastBinOrig = firstBinOrig + nBinOrigUsed - 1;
    }
  }

  printf("Rebin from %d bins to %d bins -- Used bins=%d in range %d-%d\n", nBinOrig, nBinFinal, nBinOrigUsed, firstBinOrig, lastBinOrig);
  Float_t lowLim = hOrig->GetXaxis()->GetBinLowEdge(firstBinOrig);
  Float_t hiLim = hOrig->GetXaxis()->GetBinUpEdge(lastBinOrig);
  TH1D* hRebin = new TH1D(Form("%s-rebin", hOrig->GetName()), hOrig->GetTitle(), nBinFinal, lowLim, hiLim);
  Int_t lastSummed = firstBinOrig - 1;
  for (Int_t iBin = 1; iBin <= nBinFinal; iBin++) {
    Float_t sum = 0.;
    Float_t sume2 = 0.;
    for (Int_t iOrigBin = 0; iOrigBin < reb; iOrigBin++) {
      sum += hOrig->GetBinContent(lastSummed + 1);
      sume2 += (hOrig->GetBinError(lastSummed + 1) * hOrig->GetBinError(lastSummed + 1));
      lastSummed++;
    }
    hRebin->SetBinContent(iBin, sum);
    hRebin->SetBinError(iBin, TMath::Sqrt(sume2));
  }
  return hRebin;
}
