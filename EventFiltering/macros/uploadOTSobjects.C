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
//

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <array>

#include <ROOT/RDataFrame.hxx>
#include "TFile.h"
#include "TKey.h"
#include "TH1.h"

#include "CCDB/BasicCCDBManager.h"

const std::string kBaseCCDBPath = "Users/r/rlietava/EventFiltering/OTS/";

#pragma link C++ class std::vector < std::array < uint64_t, 2>> + ;
struct bcInfo {
  bcInfo() = default;
  ULong64_t bcAOD, bcEvSel, trigMask, selMask;
  void print() const;
};
void uploadOTSobjects(std::string inputList)
{
  o2::ccdb::CcdbApi api;
  api.init("http://alice-ccdb.cern.ch");

  std::ifstream file(inputList.data());
  std::string path;
  std::map<std::string, std::string> metadata;
  while (std::getline(file, path)) {
    auto pos = path.find("/5") + 1; /// in the path at some point there is the run number
    const std::string runString = path.substr(pos, 6);
    std::cout << "Processing run " << runString << std::endl;
    const int runNumber = std::stoi(runString);
    metadata["runNumber"] = runString;
    std::pair<int64_t, int64_t> duration = o2::ccdb::BasicCCDBManager::getRunDuration(api, runNumber);
    duration.first -= 10000;   // subtract 3 minutes from the run start
    duration.second += 180000; // add 3 minutes to the run duration
    TFile scalersFile((path + "/AnalysisResults_fullrun.root").data(), "READ");
    TH1* scalers = (TH1*)scalersFile.Get("central-event-filter-task/scalers/mScalers");
    TH1* filters = (TH1*)scalersFile.Get("central-event-filter-task/scalers/mFiltered");
    api.storeAsTFile(scalers, kBaseCCDBPath + "FilterCounters", metadata, duration.first, duration.second);
    api.storeAsTFile(filters, kBaseCCDBPath + "SelectionCounters", metadata, duration.first, duration.second);
    TH1* hCounterTVX = (TH1*)scalersFile.Get("bc-selection-task/hCounterTVX");
    api.storeAsTFile(hCounterTVX, kBaseCCDBPath + "InspectedTVX", metadata, duration.first, duration.second);

    std::vector<std::array<uint64_t, 2>> bcRanges, filterBitMask, selectionBitMask;
    TFile bcRangesFile((path + "/bcRanges_fullrun.root").data(), "READ");
    int Nmax = 0;
    for (auto key : *(bcRangesFile.GetListOfKeys())) {
      TTree* cefpTree = (TTree*)bcRangesFile.Get(Form("%s/selectedBC", key->GetName()));
      if (!cefpTree)
        continue;
      bcInfo bci;
      cefpTree->SetBranchAddress("bcAO2D", &bci.bcAOD);
      cefpTree->SetBranchAddress("bcEvSel", &bci.bcEvSel);
      cefpTree->SetBranchAddress("selMask", &bci.selMask);
      cefpTree->SetBranchAddress("triMask", &bci.trigMask);
      for (int i = 0; i < cefpTree->GetEntries(); i++) {
        if ((i < Nmax) || (Nmax == 0)) {
          cefpTree->GetEntry(i);
          // Check consistency
          if (~bci.trigMask & bci.selMask) {
            std::cout << "ERROR selMask is not subset of trigMask:";
            // bcAO2D.print();
          }
          bcRanges.push_back({bci.bcAOD, bci.bcEvSel});
          filterBitMask.push_back({bci.trigMask, 0ull});
          selectionBitMask.push_back({bci.selMask, 0ull});
        }
      }
    }
    api.storeAsTFileAny(&bcRanges, kBaseCCDBPath + "SelectedBCs", metadata, duration.first, duration.second);
    api.storeAsTFileAny(&filterBitMask, kBaseCCDBPath + "FilterBitMasks", metadata, duration.first, duration.second);
    api.storeAsTFileAny(&selectionBitMask, kBaseCCDBPath + "SelectionBitMasks", metadata, duration.first, duration.second);
  }
}
