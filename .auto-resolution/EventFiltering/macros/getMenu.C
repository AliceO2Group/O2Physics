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
#include <TH1.h>

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
