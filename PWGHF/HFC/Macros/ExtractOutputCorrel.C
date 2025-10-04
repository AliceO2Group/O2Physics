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

#include <TH1.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>

#include <RtypesCore.h>

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

using namespace rapidjson;

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

void setInputCorrelNames(DhCorrelationExtraction* plotter, TString pathFileSE, TString pathFileME, TString dirSE, TString dirME, TString histoNameCorrSignal, TString histoNameCorrSideba, TString histoNameCorrSidebaLeft, TString histoNameCorrSidebaRight);
void setInputHistoInvMassNames(DhCorrelationExtraction* plotter, TString pathFileMass, std::vector<std::string> inputMassNames);
void setInputHistoFdSubtraction(DhCorrelationExtraction* plotter, TString pathFileFDTemplate, TString pathFileFDPromptFrac, TString histoNameFDTemplatePrompt, TString histoNameFDTemplateNonPrompt, TString histoNameRawFracPrompt);
void setInputHistoSecPart(DhCorrelationExtraction* plotter, TString pathFileSecPart, TString dirSecPartName, TString histoNamePrimaryPart, TString histoNameAllPart);
void setInputHistoBiasBtoD(DhCorrelationExtraction* plotter, TString pathfFilePromptMcRec, TString pathfFileNonPromptMcRec);

void extractOutputCorrelDs(const TString cfgFileName = "config_CorrAnalysis.json")
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

  std::string codeNameAnalysis = config["CodeName"].GetString();
  gSystem->Exec(Form("rm -rf Output_CorrelationExtraction_%s_Root/ Output_CorrelationExtraction_%s_png/", codeNameAnalysis.data(), codeNameAnalysis.data()));
  gSystem->Exec(Form("mkdir Output_CorrelationExtraction_%s_Root/ Output_CorrelationExtraction_%s_png/", codeNameAnalysis.data(), codeNameAnalysis.data()));

  std::string const pathFileSE = config["pathFileSE"].GetString();
  std::string const pathFileME = config["pathFileME"].GetString();
  std::string const pathFileMass = config["pathFileMass"].GetString();
  std::string const pathFileFDTemplate = config["pathFileFDTemplate"].GetString();
  std::string const pathFileFDPromptFrac = config["pathFileFDPromptFrac"].GetString();
  std::string const pathFileSecPart = config["pathFileSecPart"].GetString();
  std::string const pathfFilePromptMcRec = config["pathfFilePromptMcRec"].GetString();
  std::string const pathfFileNonPromptMcRec = config["pathfFileNonPromptMcRec"].GetString();

  std::string const dirSE = config["InputDirSE"].GetString();
  std::string const dirME = config["InputDirME"].GetString();
  std::string const dirSecPart = config["InputDirSecPart"].GetString();
  std::string const histoNameCorrSignal = config["InputHistoCorrSignalName"].GetString();
  std::string const histoNameCorrSideba = config["InputHistoCorrSidebaName"].GetString();
  std::string const histoNameCorrSidebaLeft = config["InputHistoCorrSidebaLeftName"].GetString();
  std::string const histoNameCorrSidebaRight = config["InputHistoCorrSidebaRightName"].GetString();
  std::string const histoNameFDTemplatePrompt = config["InputHistoFDTemplatePrompt"].GetString();
  std::string const histoNameFDTemplateNonPrompt = config["InputHistoFDTemplateNonPrompt"].GetString();
  std::string const histoNameRawFracPrompt = config["InputHistoFDPromptFrac"].GetString();
  std::string const histoNamePrimaryPart = config["InputHistoPrimaryPart"].GetString();
  std::string const histoNameAllPart = config["InputHistoAllPart"].GetString();

  std::vector<std::string> inputHistoMassName;

  const Value& inputMassNames = config["InputHistoMassName"];
  parseStringArray(inputMassNames, inputHistoMassName);

  std::cout << inputHistoMassName[0].data() << std::endl;
  std::cout << inputHistoMassName[1].data() << std::endl;
  std::cout << inputHistoMassName[2].data() << std::endl;

  std::vector<double> binsPtCandIntervals;
  std::vector<double> binsPtHadIntervals;
  std::vector<double> deltaEtaInterval;

  const Value& ptCandValue = config["binsPtCandIntervals"];
  readArray(ptCandValue, binsPtCandIntervals);

  const Value& ptHadValue = config["binsPtHadIntervals"];
  readArray(ptHadValue, binsPtHadIntervals);

  const Value& deltaEtaValue = config["deltaEtaInterval"];
  readArray(deltaEtaValue, deltaEtaInterval);
  double const deltaEtaMin = deltaEtaInterval[0];
  double const deltaEtaMax = deltaEtaInterval[1];

  int const specie = config["DmesonSpecie"].GetInt();
  bool const rebinAngCorr = config["RebinAngCorr"].GetBool();
  bool const rebinFDCorr = config["RebinFDCorr"].GetBool();
  bool const rebinSecPart = config["RebinSecPart"].GetBool();
  int const rebinDeltaPhi = config["nRebinDeltaPhi"].GetInt();
  int const rebinDeltaEta = config["nRebinDeltaEta"].GetInt();

  int const npools = config["NumberOfPools"].GetInt();
  bool const poolByPool = config["CorrectPoolsSeparately"].GetBool();
  bool const applySecPartCorr = config["ApplySecPartCorr"].GetBool();
  bool const applyBiasBtoDCorr = config["ApplyBiasBtoDCorr"].GetBool();
  bool const applyFDCorr = config["ApplyFDCorr"].GetBool();
  bool const isDividedSideb = config["IsDividedSideb"].GetBool();
  bool const useSidebLeft = config["UseSidebLeft"].GetBool();
  bool const useSidebRight = config["UseSidebRight"].GetBool();

  if (useSidebLeft && useSidebRight) {
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
  TH1D* hCorrectedCorrelBaselineSubtr[nBinsPtCand][nBinsPtHad];
  TH1D* hCorrectedCorrelReflected[nBinsPtCand][nBinsPtHad];
  TH1D* hCorrectedCorrelReflectedBaselineSubtr[nBinsPtCand][nBinsPtHad];

  // Create and set the correlation plotter class
  auto* plotter = new DhCorrelationExtraction();

  Bool_t const flagSpecie = plotter->setDmesonSpecie(static_cast<DhCorrelationExtraction::DmesonSpecie>(specie));
  plotter->setNpools(npools);
  plotter->setCorrectPoolsSeparately(poolByPool); // kTRUE = pool.by-pool extraction and correction; kFALSE = merged ME pools
  plotter->setFdSubtraction(applyFDCorr);
  plotter->setSecPartContamination(applySecPartCorr);
  plotter->setDeltaEtaRange(deltaEtaMin, deltaEtaMax);
  plotter->setSubtractSoftPiInMEdistr(kFALSE);
  plotter->setRebinOptions(rebinAngCorr, rebinFDCorr, rebinSecPart);
  plotter->setRebin2DcorrelHisto(rebinDeltaEta, rebinDeltaPhi); // Xaxis: deltaEta, Yaxis: deltaPhi
  plotter->setCorrBiasBtoD(applyBiasBtoDCorr);
  plotter->setDebugLevel(1);

  if (!flagSpecie) {
    std::cout << "[ERROR] Wrong D meson flag" << std::endl;
  }

  // Set the input file config
  setInputCorrelNames(plotter, pathFileSE, pathFileME, dirSE, dirME, histoNameCorrSignal, histoNameCorrSideba, histoNameCorrSidebaLeft, histoNameCorrSidebaRight);
  setInputHistoInvMassNames(plotter, pathFileMass, inputHistoMassName);
  if (applyFDCorr) {
    setInputHistoFdSubtraction(plotter, pathFileFDTemplate, pathFileFDPromptFrac, histoNameFDTemplatePrompt, histoNameFDTemplateNonPrompt, histoNameRawFracPrompt);
  }
  if (applySecPartCorr) {
    setInputHistoSecPart(plotter, pathFileSecPart, dirSecPart, histoNamePrimaryPart, histoNameAllPart);
  }
  if (applyBiasBtoDCorr) {
    setInputHistoBiasBtoD(plotter, pathfFilePromptMcRec, pathfFileNonPromptMcRec);
  }
  Bool_t const readSEandME = plotter->readInputSeAndMe();
  if (readSEandME) {
    std::cout << "Files SE and ME read correctly" << std::endl;
  }
  Bool_t const readInvMass = plotter->readInputInvMass();
  if (readInvMass) {
    std::cout << "Files inv. mass read correctly" << std::endl;
  }
  if (applyFDCorr) {
    Bool_t const readFDSubtr = plotter->readInputFdSubtr();
    if (readFDSubtr) {
      std::cout << "Files for FD subtr. read correctly" << std::endl;
    }
  }
  if (applySecPartCorr) {
    Bool_t const readSecPart = plotter->readInputSecondaryPartContamination();
    if (readSecPart) {
      std::cout << "Files for secondary part. contamination read correctly" << std::endl;
    }
  }

  // Loop over candidate pt and assoc. particle pt
  for (int iBinPtCand = 0; iBinPtCand < nBinsPtCand; iBinPtCand++) {
    plotter->setDividedSidebands(isDividedSideb, useSidebLeft, useSidebRight);
    plotter->getSignalAndBackgroundForNorm(binsPtCandIntervals[iBinPtCand], binsPtCandIntervals[iBinPtCand + 1]);
    for (int iBinPtHad = 0; iBinPtHad < nBinsPtHad; iBinPtHad++) {
      plotter->setBinCandAndHad(iBinPtCand + 1, iBinPtHad + 1);
      plotter->extractCorrelations(binsPtCandIntervals[iBinPtCand], binsPtCandIntervals[iBinPtCand + 1], binsPtHadIntervals[iBinPtHad], binsPtHadIntervals[iBinPtHad + 1], codeNameAnalysis);
      hCorrectedCorrel[iBinPtCand][iBinPtHad] = plotter->getCorrectedCorrHisto();
      hCorrectedCorrelBaselineSubtr[iBinPtCand][iBinPtHad] = plotter->getCorrectedCorrHistoBaselineSubtr();
      hCorrectedCorrelReflected[iBinPtCand][iBinPtHad] = plotter->getCorrectedCorrHistoReflected();
      hCorrectedCorrelReflectedBaselineSubtr[iBinPtCand][iBinPtHad] = plotter->getCorrectedCorrHistoReflectedBaselineSubtr();
    }
  }

  // output file
  auto* outFile = new TFile(Form("Output_CorrelationExtraction_%s_Root/ExtractCorrelationsResults.root", codeNameAnalysis.data()), "RECREATE");
  outFile->cd();
  for (int iBinPtCand = 0; iBinPtCand < nBinsPtCand; iBinPtCand++) {
    for (int iBinPtHad = 0; iBinPtHad < nBinsPtHad; iBinPtHad++) {
      hCorrectedCorrel[iBinPtCand][iBinPtHad]->Write();
    }
  }
  outFile->Close();

  // output file baseline subtr.
  auto* outFileBaselineSubtr = new TFile(Form("Output_CorrelationExtraction_%s_Root/ExtractCorrelationsResults_BaselineSubtr.root", codeNameAnalysis.data()), "RECREATE");
  outFileBaselineSubtr->cd();
  for (int iBinPtCand = 0; iBinPtCand < nBinsPtCand; iBinPtCand++) {
    for (int iBinPtHad = 0; iBinPtHad < nBinsPtHad; iBinPtHad++) {
      hCorrectedCorrelBaselineSubtr[iBinPtCand][iBinPtHad]->Write();
    }
  }
  outFileBaselineSubtr->Close();

  // output file reflected
  auto* outFileReflected = new TFile(Form("Output_CorrelationExtraction_%s_Root/ExtractCorrelationsResults_Reflected.root", codeNameAnalysis.data()), "RECREATE");
  outFileReflected->cd();
  for (int iBinPtCand = 0; iBinPtCand < nBinsPtCand; iBinPtCand++) {
    for (int iBinPtHad = 0; iBinPtHad < nBinsPtHad; iBinPtHad++) {
      hCorrectedCorrelReflected[iBinPtCand][iBinPtHad]->Write();
    }
  }
  outFileReflected->Close();

  // output file reflected baseline subtr.
  auto* outFileReflectedBaselineSubtr = new TFile(Form("Output_CorrelationExtraction_%s_Root/ExtractCorrelationsResults_Reflected_BaselineSubtr.root", codeNameAnalysis.data()), "RECREATE");
  outFileReflectedBaselineSubtr->cd();
  for (int iBinPtCand = 0; iBinPtCand < nBinsPtCand; iBinPtCand++) {
    for (int iBinPtHad = 0; iBinPtHad < nBinsPtHad; iBinPtHad++) {
      hCorrectedCorrelReflectedBaselineSubtr[iBinPtCand][iBinPtHad]->Write();
    }
  }
  outFileReflectedBaselineSubtr->Close();
}

void setInputCorrelNames(DhCorrelationExtraction* plotter, TString pathFileSE, TString pathFileME, TString dirSE, TString dirME, TString histoNameCorrSignal, TString histoNameCorrSideba, TString histoNameCorrSidebaLeft, TString histoNameCorrSidebaRight)
{

  // Ds paths
  plotter->setInputFilenameSe(pathFileSE.Data());
  plotter->setInputFilenameMe(pathFileME.Data());
  plotter->setDirNameSe(dirSE.Data());
  plotter->setDirNameMe(dirME.Data());
  plotter->setSeCorrelHistoSignalName(histoNameCorrSignal.Data());
  plotter->setSeCorrelHistoSidebandName(histoNameCorrSideba.Data());
  plotter->setMeCorrelHistoSignalName(histoNameCorrSignal.Data());
  plotter->setMeCorrelHistoSidebandName(histoNameCorrSideba.Data());
  plotter->setSeCorrelHistoSidebandLeftName(histoNameCorrSidebaLeft.Data());
  plotter->setMeCorrelHistoSidebandLeftName(histoNameCorrSidebaLeft.Data());
  plotter->setSeCorrelHistoSidebandRightName(histoNameCorrSidebaRight.Data());
  plotter->setMeCorrelHistoSidebandRightName(histoNameCorrSidebaRight.Data());
}

void setInputHistoInvMassNames(DhCorrelationExtraction* plotter, TString pathFileMass, std::vector<std::string> inputMassNames)
{ // to use if sgn and bkg extraction is done apart

  plotter->setInputFilenameMass(pathFileMass.Data());
  plotter->setMassHistoNameSgn(inputMassNames[0].data());
  plotter->setMassHistoNameBkg(inputMassNames[1].data());
  plotter->setMassHistoNameSBs(inputMassNames[2].data());
}

void setInputHistoFdSubtraction(DhCorrelationExtraction* plotter, TString pathFileFDTemplate, TString pathFileFDPromptFrac, TString histoNameFDTemplatePrompt, TString histoNameFDTemplateNonPrompt, TString histoNameRawFracPrompt)
{

  plotter->setInputFilenameFdTemplate(pathFileFDTemplate.Data());
  plotter->setInputFilenameFdPromptFrac(pathFileFDPromptFrac.Data());
  plotter->setInputHistoNameFdTemplatePrompt(histoNameFDTemplatePrompt.Data());
  plotter->setInputHistoNameFdTemplateNonPrompt(histoNameFDTemplateNonPrompt.Data());
  plotter->setInputHistoNameFdPromptFrac(histoNameRawFracPrompt.Data());
}

void setInputHistoSecPart(DhCorrelationExtraction* plotter, TString pathFileSecPart, TString dirSecPartName, TString histoNamePrimaryPart, TString histoNameAllPart)
{

  plotter->setInputFilenameSecPart(pathFileSecPart.Data());
  plotter->setDirNameSecPart(dirSecPartName.Data());
  plotter->setHistoSecPartName(histoNamePrimaryPart.Data(), histoNameAllPart.Data());
}

void setInputHistoBiasBtoD(DhCorrelationExtraction* plotter, TString pathfFilePromptMcRec, TString pathfFileNonPromptMcRec)
{

  plotter->setInputFilenameBiasBtoD(pathfFilePromptMcRec.Data(), pathfFileNonPromptMcRec.Data());
}
