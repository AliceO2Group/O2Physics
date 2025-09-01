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

/// \file DhCorrelationExtraction.cxx
/// \brief class for D-h correlation extraction
/// \usage .L DhCorrelationExtraction.cxx+
/// \usage .x ExtractOutputCorrel.C("config-file-name")
/// \author Samuele Cattaruzzi <samuele.cattaruzzi@cern.ch>
/// \author Swapnesh Santosh Khade <swapnesh.santosh.khade@cern.ch>

#include "DhCorrelationExtraction.h"
#include "Riostream.h"

#include <TROOT.h>
#include <TStyle.h>

#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

using namespace rapidjson;

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

void SetInputCorrelNames(DhCorrelationExtraction* plotter, TString pathFileSE, TString pathFileME, TString dirSE, TString dirME, TString histoNameCorrSignal, TString histoNameCorrSideba, TString histoNameCorrSidebaLeft, TString histoNameCorrSidebaRight);
void SetInputHistoInvMassNames(DhCorrelationExtraction* plotter, TString pathFileMass, std::vector<std::string> inputMassNames);
void SetInputHistoFDSubtraction(DhCorrelationExtraction* plotter, TString pathFileFDTemplate, TString pathFileFDPromptFrac, TString histoNameFDTemplatePrompt, TString histoNameFDTemplateNonPrompt, TString histoNameRawFracPrompt);
void SetInputHistoSecPart(DhCorrelationExtraction* plotter, TString pathFileSecPart, TString dirSecPartName, TString histoNamePrimaryPart, TString histoNameAllPart);
void SetInputHistoBiasBtoD(DhCorrelationExtraction* plotter, TString pathfFilePromptMcRec, TString pathfFileNonPromptMcRec);

void ExtractOutputCorrel_Ds(const TString cfgFileName = "config_CorrAnalysis.json")
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
  string pathFileFDTemplate = config["pathFileFDTemplate"].GetString();
  string pathFileFDPromptFrac = config["pathFileFDPromptFrac"].GetString();
  string pathFileSecPart = config["pathFileSecPart"].GetString();
  string pathfFilePromptMcRec = config["pathfFilePromptMcRec"].GetString();
  string pathfFileNonPromptMcRec = config["pathfFileNonPromptMcRec"].GetString();

  string dirSE = config["InputDirSE"].GetString();
  string dirME = config["InputDirME"].GetString();
  string dirSecPart = config["InputDirSecPart"].GetString();
  string histoNameCorrSignal = config["InputHistoCorrSignalName"].GetString();
  string histoNameCorrSideba = config["InputHistoCorrSidebaName"].GetString();
  string histoNameCorrSidebaLeft = config["InputHistoCorrSidebaLeftName"].GetString();
  string histoNameCorrSidebaRight = config["InputHistoCorrSidebaRightName"].GetString();
  string histoNameFDTemplatePrompt = config["InputHistoFDTemplatePrompt"].GetString();
  string histoNameFDTemplateNonPrompt = config["InputHistoFDTemplateNonPrompt"].GetString();
  string histoNameRawFracPrompt = config["InputHistoFDPromptFrac"].GetString();
  string histoNamePrimaryPart = config["InputHistoPrimaryPart"].GetString();
  string histoNameAllPart = config["InputHistoAllPart"].GetString();

  std::vector<std::string> InputHistoMassName;

  const Value& inputMassNames = config["InputHistoMassName"];
  parseStringArray(inputMassNames, InputHistoMassName);

  std::cout << InputHistoMassName[0].data() << std::endl;
  std::cout << InputHistoMassName[1].data() << std::endl;
  std::cout << InputHistoMassName[2].data() << std::endl;

  std::vector<double> binsPtCandIntervals;
  std::vector<double> binsPtHadIntervals;
  std::vector<double> deltaEtaInterval;

  const Value& PtCandValue = config["binsPtCandIntervals"];
  readArray(PtCandValue, binsPtCandIntervals);

  const Value& PtHadValue = config["binsPtHadIntervals"];
  readArray(PtHadValue, binsPtHadIntervals);

  const Value& deltaEtaValue = config["deltaEtaInterval"];
  readArray(deltaEtaValue, deltaEtaInterval);
  double deltaEtaMin = deltaEtaInterval[0];
  double deltaEtaMax = deltaEtaInterval[1];

  int specie = config["DmesonSpecie"].GetInt();
  bool rebinAngCorr = config["RebinAngCorr"].GetBool();
  bool rebinFDCorr = config["RebinFDCorr"].GetBool();
  bool rebinSecPart = config["RebinSecPart"].GetBool();
  int rebinDeltaPhi = config["nRebinDeltaPhi"].GetInt();
  int rebinDeltaEta = config["nRebinDeltaEta"].GetInt();

  int npools = config["NumberOfPools"].GetInt();
  bool poolByPool = config["CorrectPoolsSeparately"].GetBool();
  bool applySecPartCorr = config["ApplySecPartCorr"].GetBool();
  bool applyBiasBtoDCorr = config["ApplyBiasBtoDCorr"].GetBool();
  bool applyFDCorr = config["ApplyFDCorr"].GetBool();
  bool isDividedSideb = config["IsDividedSideb"].GetBool();
  bool useSidebLeft = config["UseSidebLeft"].GetBool();
  bool useSidebRight = config["UseSidebRight"].GetBool();

  if (useSidebLeft && useSidebLeft) {
    std::cout << "Using left and right" << std::endl;
  }

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
  TH1D* hCorrectedCorrel_BaselineSubtr[nBinsPtCand][nBinsPtHad];
  TH1D* hCorrectedCorrel_Reflected[nBinsPtCand][nBinsPtHad];
  TH1D* hCorrectedCorrel_Reflected_BaselineSubtr[nBinsPtCand][nBinsPtHad];

  // Create and set the correlation plotter class
  DhCorrelationExtraction* plotter = new DhCorrelationExtraction();

  Bool_t flagSpecie = plotter->SetDmesonSpecie(static_cast<DhCorrelationExtraction::DmesonSpecie>(specie));
  plotter->SetNpools(npools);
  plotter->SetCorrectPoolsSeparately(poolByPool); // kTRUE = pool.by-pool extraction and correction; kFALSE = merged ME pools
  plotter->SetFDSubtraction(applyFDCorr);
  plotter->SetSecPartContamination(applySecPartCorr);
  plotter->SetDeltaEtaRange(deltaEtaMin, deltaEtaMax);
  plotter->SetSubtractSoftPiInMEdistr(kFALSE);
  plotter->SetRebinOptions(rebinAngCorr, rebinFDCorr, rebinSecPart);
  plotter->SetRebin2DcorrelHisto(rebinDeltaEta, rebinDeltaPhi); // Xaxis: deltaEta, Yaxis: deltaPhi
  plotter->SetCorrBiasBtoD(applyBiasBtoDCorr);
  plotter->SetDebugLevel(1);

  if (!flagSpecie)
    std::cout << "[ERROR] Wrong D meson flag" << std::endl;

  // Set the input file config
  SetInputCorrelNames(plotter, pathFileSE, pathFileME, dirSE, dirME, histoNameCorrSignal, histoNameCorrSideba, histoNameCorrSidebaLeft, histoNameCorrSidebaRight);
  SetInputHistoInvMassNames(plotter, pathFileMass, InputHistoMassName);
  if (applyFDCorr)
    SetInputHistoFDSubtraction(plotter, pathFileFDTemplate, pathFileFDPromptFrac, histoNameFDTemplatePrompt, histoNameFDTemplateNonPrompt, histoNameRawFracPrompt);
  if (applySecPartCorr)
    SetInputHistoSecPart(plotter, pathFileSecPart, dirSecPart, histoNamePrimaryPart, histoNameAllPart);
  if (applyBiasBtoDCorr)
    SetInputHistoBiasBtoD(plotter, pathfFilePromptMcRec, pathfFileNonPromptMcRec);
  Bool_t readSEandME = plotter->ReadInputSEandME();
  if (readSEandME)
    std::cout << "Files SE and ME read correctly" << std::endl;
  Bool_t readInvMass = plotter->ReadInputInvMass();
  if (readInvMass)
    std::cout << "Files inv. mass read correctly" << std::endl;
  if (applyFDCorr) {
    Bool_t readFDSubtr = plotter->ReadInputFDSubtr();
    if (readFDSubtr)
      std::cout << "Files for FD subtr. read correctly" << std::endl;
  }
  if (applySecPartCorr) {
    Bool_t readSecPart = plotter->ReadInputSecondaryPartContamination();
    if (readSecPart)
      std::cout << "Files for secondary part. contamination read correctly" << std::endl;
  }

  // Loop over candidate pt and assoc. particle pt
  for (int iBinPtCand = 0; iBinPtCand < nBinsPtCand; iBinPtCand++) {
    plotter->SetDividedSidebands(isDividedSideb, useSidebLeft, useSidebRight);
    plotter->GetSignalAndBackgroundForNorm(binsPtCandIntervals[iBinPtCand], binsPtCandIntervals[iBinPtCand + 1]);
    for (int iBinPtHad = 0; iBinPtHad < nBinsPtHad; iBinPtHad++) {
      plotter->SetBinCandAndHad(iBinPtCand + 1, iBinPtHad + 1);
      plotter->ExtractCorrelations(binsPtCandIntervals[iBinPtCand], binsPtCandIntervals[iBinPtCand + 1], binsPtHadIntervals[iBinPtHad], binsPtHadIntervals[iBinPtHad + 1], CodeNameAnalysis);
      hCorrectedCorrel[iBinPtCand][iBinPtHad] = (TH1D*)plotter->GetCorrectedCorrHisto();
      hCorrectedCorrel_BaselineSubtr[iBinPtCand][iBinPtHad] = (TH1D*)plotter->GetCorrectedCorrHisto_BaselineSubtr();
      hCorrectedCorrel_Reflected[iBinPtCand][iBinPtHad] = (TH1D*)plotter->GetCorrectedCorrHisto_Reflected();
      hCorrectedCorrel_Reflected_BaselineSubtr[iBinPtCand][iBinPtHad] = (TH1D*)plotter->GetCorrectedCorrHisto_Reflected_BaselineSubtr();
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

  // output file baseline subtr.
  TFile* outFile_BaselineSubtr = new TFile(Form("Output_CorrelationExtraction_%s_Root/ExtractCorrelationsResults_BaselineSubtr.root", CodeNameAnalysis.data()), "RECREATE");
  outFile_BaselineSubtr->cd();
  for (int iBinPtCand = 0; iBinPtCand < nBinsPtCand; iBinPtCand++) {
    for (int iBinPtHad = 0; iBinPtHad < nBinsPtHad; iBinPtHad++) {
      hCorrectedCorrel_BaselineSubtr[iBinPtCand][iBinPtHad]->Write();
    }
  }
  outFile_BaselineSubtr->Close();

  // output file reflected
  TFile* outFile_Reflected = new TFile(Form("Output_CorrelationExtraction_%s_Root/ExtractCorrelationsResults_Reflected.root", CodeNameAnalysis.data()), "RECREATE");
  outFile_Reflected->cd();
  for (int iBinPtCand = 0; iBinPtCand < nBinsPtCand; iBinPtCand++) {
    for (int iBinPtHad = 0; iBinPtHad < nBinsPtHad; iBinPtHad++) {
      hCorrectedCorrel_Reflected[iBinPtCand][iBinPtHad]->Write();
    }
  }
  outFile_Reflected->Close();

  // output file reflected baseline subtr.
  TFile* outFile_Reflected_BaselineSubtr = new TFile(Form("Output_CorrelationExtraction_%s_Root/ExtractCorrelationsResults_Reflected_BaselineSubtr.root", CodeNameAnalysis.data()), "RECREATE");
  outFile_Reflected_BaselineSubtr->cd();
  for (int iBinPtCand = 0; iBinPtCand < nBinsPtCand; iBinPtCand++) {
    for (int iBinPtHad = 0; iBinPtHad < nBinsPtHad; iBinPtHad++) {
      hCorrectedCorrel_Reflected_BaselineSubtr[iBinPtCand][iBinPtHad]->Write();
    }
  }
  outFile_Reflected_BaselineSubtr->Close();

  return;
}

void SetInputCorrelNames(DhCorrelationExtraction* plotter, TString pathFileSE, TString pathFileME, TString dirSE, TString dirME, TString histoNameCorrSignal, TString histoNameCorrSideba, TString histoNameCorrSidebaLeft, TString histoNameCorrSidebaRight)
{

  // Ds paths
  plotter->SetInputFilenameSE(pathFileSE.Data());
  plotter->SetInputFilenameME(pathFileME.Data());
  plotter->SetDirNameSE(dirSE.Data());
  plotter->SetDirNameME(dirME.Data());
  plotter->SetSECorrelHistoSignalName(histoNameCorrSignal.Data());
  plotter->SetSECorrelHistoSidebandName(histoNameCorrSideba.Data());
  plotter->SetMECorrelHistoSignalName(histoNameCorrSignal.Data());
  plotter->SetMECorrelHistoSidebandName(histoNameCorrSideba.Data());
  plotter->SetSECorrelHistoSidebandLeftName(histoNameCorrSidebaLeft.Data());
  plotter->SetMECorrelHistoSidebandLeftName(histoNameCorrSidebaLeft.Data());
  plotter->SetSECorrelHistoSidebandRightName(histoNameCorrSidebaRight.Data());
  plotter->SetMECorrelHistoSidebandRightName(histoNameCorrSidebaRight.Data());

  return;
}

void SetInputHistoInvMassNames(DhCorrelationExtraction* plotter, TString pathFileMass, std::vector<std::string> inputMassNames)
{ // to use if sgn and bkg extraction is done apart

  plotter->SetInputFilenameMass(pathFileMass.Data());
  plotter->SetMassHistoNameSgn(inputMassNames[0].data());
  plotter->SetMassHistoNameBkg(inputMassNames[1].data());
  plotter->SetMassHistoNameSBs(inputMassNames[2].data());

  return;
}

void SetInputHistoFDSubtraction(DhCorrelationExtraction* plotter, TString pathFileFDTemplate, TString pathFileFDPromptFrac, TString histoNameFDTemplatePrompt, TString histoNameFDTemplateNonPrompt, TString histoNameRawFracPrompt)
{

  plotter->SetInputFilenameFDTemplate(pathFileFDTemplate.Data());
  plotter->SetInputFilenameFDPromptFrac(pathFileFDPromptFrac.Data());
  plotter->SetInputHistoNameFDTemplatePrompt(histoNameFDTemplatePrompt.Data());
  plotter->SetInputHistoNameFDTemplateNonPrompt(histoNameFDTemplateNonPrompt.Data());
  plotter->SetInputHistoNameFDPromptFrac(histoNameRawFracPrompt.Data());

  return;
}

void SetInputHistoSecPart(DhCorrelationExtraction* plotter, TString pathFileSecPart, TString dirSecPartName, TString histoNamePrimaryPart, TString histoNameAllPart)
{

  plotter->SetInputFilenameSecPart(pathFileSecPart.Data());
  plotter->SetDirNameSecPart(dirSecPartName.Data());
  plotter->SetHistoSecPartName(histoNamePrimaryPart.Data(), histoNameAllPart.Data());

  return;
}

void SetInputHistoBiasBtoD(DhCorrelationExtraction* plotter, TString pathfFilePromptMcRec, TString pathfFileNonPromptMcRec)
{

  plotter->SetInputFilenameBiasBtoD(pathfFilePromptMcRec.Data(), pathfFileNonPromptMcRec.Data());

  return;
}
