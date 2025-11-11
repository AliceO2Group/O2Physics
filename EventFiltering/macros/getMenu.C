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

#include "CCDB/BasicCCDBManager.h"

#include <TAxis.h>
#include <TFile.h>
#include <TGrid.h>
#include <TH1.h>
#include <TSystem.h>

#include <iostream>
#include <regex>
#include <string>

void getMenu(int runNumber, std::string baseCCDBPath = "Users/m/mpuccio/EventFiltering/OTS/Chunked/")
{
  auto& ccdb = o2::ccdb::BasicCCDBManager::instance();
  TH1* counters = ccdb.getForRun<TH1>(baseCCDBPath + "FilterCounters", runNumber);
  TAxis* axis = counters->GetXaxis();

  std::vector<std::string> binLabels(axis->GetNbins() - 2); // skip first and last bins
  std::cout << "Menu for run " << runNumber << ":\n";
  for (int i = 2; i < axis->GetNbins(); ++i) {
    binLabels[i - 1] = axis->GetBinLabel(i);
    std::cout << "Id " << i - 2 << ": " << axis->GetBinLabel(i) << "\n";
  }
}

void getMenu(std::string period)
{
  std::regex pattern(R"(LHC(\d{2})[A-Za-z]{1,2})");
  std::smatch match;

  int year{2000};
  if (std::regex_match(period, match, pattern)) {
    std::cout << "Year = " << match[1] << std::endl;
    year += std::stoi(match[1]);
  } else {
    std::cout << "Invalid format" << std::endl;
  }
  gSystem->Exec(Form("alien_find /alice/data/%i/%s/ ctf_skim_full/AnalysisResults_fullrun.root > list_tmp.txt", year, period.data()));

  std::ifstream file("list_tmp.txt");
  if (!file) {
    std::cerr << "Error: could not open file\n";
    return;
  }

  std::string firstLine;
  if (!std::getline(file, firstLine)) {
    std::cerr << "Error: file is empty or read failed\n";
    return;
  }

  std::cout << "First line: " << firstLine << '\n';
  TGrid::Connect("alien://");
  TFile *scalersFile = TFile::Open((std::string("alien://") + firstLine).data(), "READ");
  TH1D* counters = (TH1D*)scalersFile->Get("central-event-filter-task/scalers/mFiltered");
  TAxis* axis = counters->GetXaxis();

  std::vector<std::string> binLabels(axis->GetNbins() - 2); // skip first and last bins
  for (int i = 2; i < axis->GetNbins(); ++i) {
    binLabels[i - 1] = axis->GetBinLabel(i);
    std::cout << "Id " << i - 2 << ": " << axis->GetBinLabel(i) << "\n";
  }
  gSystem->Exec("rm list_tmp.txt");
}
