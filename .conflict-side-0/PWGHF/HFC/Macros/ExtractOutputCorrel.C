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

/// \file ExtractOutputCorrel.C
/// \brief Macro to perform the correlation extraction
/// \usage .L DhCorrelationExtraction.cxx+
/// \usage .x ExtractOutputCorrel.C("config-file-name")
/// \author Samuele Cattaruzzi <samuele.cattaruzzi@cern.ch>
/// \author Swapnesh Santosh Khade <swapnesh.santosh.khade@cern.ch>

#include "Riostream.h"
#include <TROOT.h>
#include <TStyle.h>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include "DhCorrelationExtraction.h"

using namespace rapidjson;

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

void SetInputCorrelNames(DhCorrelationExtraction* plotter, TString pathFileMass, TString pathFileSE, TString pathFileME, TString dirSE, TString dirME, TString histoNameCorrSignal, TString histoNameCorrSideba);
void SetInputHistoInvMassNames(DhCorrelationExtraction* plotter, std::vector<string> inputMassNames);

void ExtractOutputCorrel(TString cfgFileName = "config_CorrAnalysis.json")
{
  // gStyle -> SetOptStat(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLineWidth(2);
  gStyle->SetCanvasDefH(1126);
  gStyle->SetCanvasDefW(1840);

  // Load config
  FILE* configFile = fopen(cfgFileName.Data(), "r");
  Document config;
  char readBuffer[65536];
  FileReadStream is(configFile, readBuffer, sizeof(readBuffer));
  config.ParseStream(is);
  fclose(configFile);

  string CodeNameAnalysis = config["CodeName"].GetString();
  gSystem->Exec(Form("rm -rf Output_CorrelationExtraction_%s_Root/ Output_CorrelationExtraction_%s_png/", CodeNameAnalysis.data(), CodeNameAnalysis.data()));
  gSystem->Exec(Form("mkdir Output_CorrelationExtraction_%s_Root/ Output_CorrelationExtraction_%s_png/", CodeNameAnalysis.data(), CodeNameAnalysis.data()));

  string pathFileSE = config["pathFileSE"].GetString();
  string pathFileME = config["pathFileME"].GetString();
  string pathFileMass = config["pathFileMass"].GetString();

  string dirSE = config["InputDirSE"].GetString();
  string dirME = config["InputDirME"].GetString();
  string histoNameCorrSignal = config["InputHistoCorrSignalName"].GetString();
  string histoNameCorrSideba = config["InputHistoCorrSidebaName"].GetString();

  vector<string> InputHistoMassName;

  const Value& inputMassNames = config["InputHistoMassName"];
  parseStringArray(inputMassNames, InputHistoMassName);

  cout << InputHistoMassName[0].data() << endl;
  cout << InputHistoMassName[1].data() << endl;
  cout << InputHistoMassName[2].data() << endl;

  vector<double> binsPtCandIntervals;
  vector<double> binsPtHadIntervals;
  vector<double> deltaEtaInterval;

  const Value& PtCandValue = config["binsPtCandIntervals"];
  readArray(PtCandValue, binsPtCandIntervals);

  const Value& PtHadValue = config["binsPtHadIntervals"];
  readArray(PtHadValue, binsPtHadIntervals);

  const Value& deltaEtaValue = config["deltaEtaInterval"];
  readArray(deltaEtaValue, deltaEtaInterval);
  double deltaEtaMin = deltaEtaInterval[0];
  double deltaEtaMax = deltaEtaInterval[1];

  int specie = config["DmesonSpecie"].GetInt();

  int npools = config["NumberOfPools"].GetInt();
  bool poolByPool = config["CorrectPoolsSeparately"].GetBool();

  std::cout << "=========================== " << std::endl;
  std::cout << "Input variables from config" << std::endl;
  std::cout << "deltaEtaMin    = " << deltaEtaMin << std::endl;
  std::cout << "deltaEtaMax    = " << deltaEtaMax << std::endl;
  std::cout << "DmesonSpecie    = " << specie << std::endl;
  std::cout << "nPools    = " << npools << std::endl;
  std::cout << "poolByPool    = " << poolByPool << std::endl;
  std::cout << "=========================== " << std::endl;
  std::cout << " " << std::endl;

  const int nBinsPtCand = binsPtCandIntervals.size() - 1;
  const int nBinsPtHad = binsPtHadIntervals.size() - 1;

  TH1D* hCorrectedCorrel[nBinsPtCand][nBinsPtHad];

  // Create and set the correlation plotter class
  DhCorrelationExtraction* plotter = new DhCorrelationExtraction();

  Bool_t flagSpecie = plotter->SetDmesonSpecie(static_cast<DhCorrelationExtraction::DmesonSpecie>(specie));
  plotter->SetNpools(npools);
  plotter->SetCorrectPoolsSeparately(poolByPool); // kTRUE = pool.by-pool extraction and correction; kFALSE = merged ME pools
  plotter->SetDeltaEtaRange(deltaEtaMin, deltaEtaMax);
  plotter->SetSubtractSoftPiInMEdistr(kFALSE);
  plotter->SetRebin2DcorrelHisto(2, 2); // Xaxis: deltaEta, Yaxis: deltaPhi
  plotter->SetDebugLevel(1);

  if (!flagSpecie)
    cout << "[ERROR] Wrong D meson flag" << endl;

  // Set the input file config
  SetInputCorrelNames(plotter, pathFileMass, pathFileSE, pathFileME, dirSE, dirME, histoNameCorrSignal, histoNameCorrSideba);
  SetInputHistoInvMassNames(plotter, InputHistoMassName);
  Bool_t readSEandME = plotter->ReadInputSEandME();
  Bool_t readInvMass = plotter->ReadInputInvMass();
  if (readSEandME)
    cout << "Files SE and ME read correctly" << endl;
  if (readInvMass)
    cout << "Files inv. mass read correctly" << endl;

  // Loop over candidate pt and assoc. particle pt
  for (int iBinPtCand = 0; iBinPtCand < nBinsPtCand; iBinPtCand++) {
    plotter->GetSignalAndBackgroundForNorm(binsPtCandIntervals[iBinPtCand], binsPtCandIntervals[iBinPtCand + 1]);
    for (int iBinPtHad = 0; iBinPtHad < nBinsPtHad; iBinPtHad++) {
      plotter->ExtractCorrelations(binsPtCandIntervals[iBinPtCand], binsPtCandIntervals[iBinPtCand + 1], binsPtHadIntervals[iBinPtHad], binsPtHadIntervals[iBinPtHad + 1], CodeNameAnalysis);
      hCorrectedCorrel[iBinPtCand][iBinPtHad] = reinterpret_cast<TH1D*>(plotter->GetCorrectedCorrHisto());
    }
  }

  // output file
  TFile* outFile = new TFile(Form("Output_CorrelationExtraction_%s_Root/ExtractCorrelationsResults.root", CodeNameAnalysis.data()), "RECREATE");
  outFile->cd();
  for (int iBinPtCand = 0; iBinPtCand < nBinsPtCand; iBinPtCand++) {
    for (int iBinPtHad = 0; iBinPtHad < nBinsPtHad; iBinPtHad++) {
      hCorrectedCorrel[iBinPtCand][iBinPtHad]->Write();
    }
  }
  outFile->Close();

  return;
}

void SetInputCorrelNames(DhCorrelationExtraction* plotter, TString pathFileMass, TString pathFileSE, TString pathFileME, TString dirSE, TString dirME, TString histoNameCorrSignal, TString histoNameCorrSideba)
{

  // paths
  plotter->SetInputFilenameMass(pathFileMass.Data());
  plotter->SetInputFilenameSE(pathFileSE.Data());
  plotter->SetInputFilenameME(pathFileME.Data());
  plotter->SetDirNameSE(dirSE.Data());
  plotter->SetDirNameME(dirME.Data());
  plotter->SetSECorrelHistoSignalName(histoNameCorrSignal.Data());
  plotter->SetSECorrelHistoSidebandName(histoNameCorrSideba.Data());
  plotter->SetMECorrelHistoSignalName(histoNameCorrSignal.Data());
  plotter->SetMECorrelHistoSidebandName(histoNameCorrSideba.Data());

  return;
}

void SetInputHistoInvMassNames(DhCorrelationExtraction* plotter, std::vector<std::string> inputMassNames)
{ // to use if sgn and bkg extraction is done apart

  plotter->SetMassHistoNameSgn(inputMassNames[0].data());
  plotter->SetMassHistoNameBkg(inputMassNames[1].data());
  plotter->SetMassHistoNameSBs(inputMassNames[2].data());

  return;
}
