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
#include <set>
#include <sstream>
#include <string>
#include <vector>

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

std::vector<std::string> getMenuForPeriod(std::string period)
{
  std::regex pattern(R"(LHC(\d{2})[A-Za-z]{1,2})");
  std::smatch match;

  int year{2000};
  if (!std::regex_match(period, match, pattern)) {
    std::cout << "Invalid format for period: " << period << std::endl;
    return {};
  }

  year += std::stoi(match[1]);
  gSystem->Exec(Form("alien_find /alice/data/%i/%s/ ctf_skim_full/AnalysisResults_fullrun.root > list_tmp_%s.txt", year, period.data(), period.data()));

  std::ifstream file(Form("list_tmp_%s.txt", period.data()));
  if (!file) {
    std::cerr << "Error: could not open file for period " << period << "\n";
    return {};
  }

  std::string firstLine;
  if (!std::getline(file, firstLine)) {
    std::cerr << "Error: file is empty or read failed for period " << period << "\n";
    return {};
  }

  TGrid::Connect("alien://");
  TFile* scalersFile = TFile::Open((std::string("alien://") + firstLine).data(), "READ");
  TH1D* counters = (TH1D*)scalersFile->Get("central-event-filter-task/scalers/mFiltered");
  TAxis* axis = counters->GetXaxis();

  std::vector<std::string> binLabels(axis->GetNbins() - 2);
  for (int i = 2; i < axis->GetNbins(); ++i) {
    binLabels[i - 2] = axis->GetBinLabel(i);
  }

  scalersFile->Close();
  delete scalersFile;
  gSystem->Exec(Form("rm list_tmp_%s.txt", period.data()));

  return binLabels;
}

void getMenu(std::string periods)
{
  std::stringstream ss(periods);
  std::string period;
  std::vector<std::string> periodList;

  // Parse comma-separated periods
  while (std::getline(ss, period, ',')) {
    // Trim whitespace
    period.erase(0, period.find_first_not_of(" \t"));
    period.erase(period.find_last_not_of(" \t") + 1);
    periodList.push_back(period);
  }

  std::map<std::vector<std::string>, std::vector<std::string>> menuGroups;

  // Get menus for each period
  for (const auto& p : periodList) {
    auto menu = getMenuForPeriod(p);
    if (!menu.empty()) {
      menuGroups[menu].push_back(p);
    }
  }

  // Report different menus
  int menuId = 1;
  for (const auto& [menu, periods] : menuGroups) {
    std::cout << "\n=== Menu " << menuId++ << " (periods: ";
    for (size_t i = 0; i < periods.size(); ++i) {
      std::cout << periods[i];
      if (i < periods.size() - 1)
        std::cout << ", ";
    }
    std::cout << ") ===\n";

    for (size_t i = 0; i < menu.size(); ++i) {
      std::cout << "Id " << i << ": " << menu[i] << "\n";
    }
  }
}
