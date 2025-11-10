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

#include "Common/Core/ZorroHelper.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/LHCConstants.h"

#include "TFile.h"
#include "TGrid.h"
#include "TH1.h"
#include "TKey.h"
#include "TSystem.h"
#include "TTree.h"

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

constexpr uint32_t chunkSize = 1000000;

void uploadOTSobjects(std::string inputList, std::string passName, bool useAlien, bool chunkedProcessing = true)
{
  const std::string kBaseCCDBPath = "EventFiltering/Zorro/";
  std::string baseCCDBpath = passName.empty() ? kBaseCCDBPath : kBaseCCDBPath + passName + "/";
  if (useAlien) {
    TGrid::Connect("alien://");
  }
  o2::ccdb::CcdbApi api;
  api.init("http://alice-ccdb.cern.ch");
  auto& ccdb = o2::ccdb::BasicCCDBManager::instance();

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
    duration.first -= 60000;  // subtract 1 minutes from the run start
    duration.second += 60000; // add 1 minutes to the run duration
    std::cout << ">>> Begin - end timestamps for the upload: " << duration.first << " - " << duration.second << std::endl;
    auto ctp = ccdb.getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", duration.first / 2 + duration.second / 2);
    auto orbitResetTimestamp = (*ctp)[0];
    path = useAlien ? "alien://" + path : path;
    std::unique_ptr<TFile> scalersFile{TFile::Open((path + "/AnalysisResults_fullrun.root").data(), "READ")};
    TH1* scalers = static_cast<TH1*>(scalersFile->Get("central-event-filter-task/scalers/mScalers"));
    TH1* filters = static_cast<TH1*>(scalersFile->Get("central-event-filter-task/scalers/mFiltered"));
    api.storeAsTFile(scalers, baseCCDBpath + "FilterCounters", metadata, duration.first, duration.second + 1);
    api.storeAsTFile(filters, baseCCDBpath + "SelectionCounters", metadata, duration.first, duration.second + 1);
    TH1* hCounterTVX = static_cast<TH1*>(scalersFile->Get("bc-selection-task/hCounterTVX"));
    if (!hCounterTVX) {
      hCounterTVX = static_cast<TH1*>(scalersFile->Get("lumi-task/hCounterTVX"));
      if (!hCounterTVX) {
        std::cout << "No hCounterTVX histogram found in the file, skipping upload for run " << runString << std::endl;
        continue;
      }
    }
    api.storeAsTFile(hCounterTVX, baseCCDBpath + "InspectedTVX", metadata, duration.first, duration.second + 1);

    std::vector<ZorroHelper> zorroHelpers;
    std::unique_ptr<TFile> bcRangesFile{TFile::Open((path + "/bcRanges_fullrun.root").data(), "READ")};
    int Nmax = 0;
    for (auto key : *(bcRangesFile->GetListOfKeys())) {
      TTree* cefpTree = static_cast<TTree*>(bcRangesFile->Get(Form("%s/selectedBC", key->GetName())));
      if (!cefpTree)
        continue;
      ZorroHelper bci;
      cefpTree->SetBranchAddress("bcAO2D", &bci.bcAOD);
      cefpTree->SetBranchAddress("bcEvSel", &bci.bcEvSel);
      if (cefpTree->GetBranch("selMask") && cefpTree->GetBranch("triMask")) {
        cefpTree->SetBranchAddress("selMask", &bci.selMask[0]);
        cefpTree->SetBranchAddress("triMask", &bci.trigMask[0]);
      } else {
        cefpTree->SetBranchAddress("selMask0", &bci.selMask[0]);
        cefpTree->SetBranchAddress("triMask0", &bci.trigMask[0]);
        cefpTree->SetBranchAddress("selMask1", &bci.selMask[1]);
        cefpTree->SetBranchAddress("triMask1", &bci.trigMask[1]);
      }
      for (int i = 0; i < cefpTree->GetEntries(); i++) {
        if ((i < Nmax) || (Nmax == 0)) {
          cefpTree->GetEntry(i);
          zorroHelpers.push_back(bci);
        }
      }
    }
    std::sort(zorroHelpers.begin(), zorroHelpers.end(), [](const ZorroHelper& a, const ZorroHelper& b) {
      return a.bcAOD < b.bcAOD;
    });
    if (!chunkedProcessing) {
      api.storeAsTFileAny(&zorroHelpers, baseCCDBpath + "ZorroHelpers", metadata, duration.first, duration.second + 1);
      std::cout << std::endl;
    } else {
      uint32_t helperIndex{0};
      auto startTS = duration.first;
      int64_t endTS = 0;
      while (helperIndex < zorroHelpers.size()) {
        std::vector<ZorroHelper> chunk;
        uint32_t endIndex = helperIndex + chunkSize;
        while (true) {
          if (endIndex >= zorroHelpers.size() - 2) {
            endTS = duration.second;
            chunk.insert(chunk.begin(), zorroHelpers.begin() + helperIndex, zorroHelpers.end());
            break;
          }
          auto bcEnd{zorroHelpers[endIndex].bcAOD > zorroHelpers[endIndex].bcEvSel ? zorroHelpers[endIndex].bcAOD : zorroHelpers[endIndex].bcEvSel};
          const auto& nextHelper = zorroHelpers[endIndex + 1];
          auto bcNext{nextHelper.bcAOD < nextHelper.bcEvSel ? nextHelper.bcAOD : nextHelper.bcEvSel};
          if (bcNext - bcEnd > 2001 / o2::constants::lhc::LHCBunchSpacingMUS) { /// ensure a gap of 2ms between chunks
            chunk.insert(chunk.begin(), zorroHelpers.begin() + helperIndex, zorroHelpers.begin() + endIndex + 1);
            endTS = (orbitResetTimestamp + int64_t(bcEnd * o2::constants::lhc::LHCBunchSpacingNS * 1e-3)) / 1000 + 1;
            break;
          }
          bcEnd = nextHelper.bcAOD > nextHelper.bcEvSel ? nextHelper.bcAOD : nextHelper.bcEvSel;
          endIndex++;
        }
        std::cout << ">>> Chunk " << helperIndex << " - " << helperIndex + chunk.size() << " : " << startTS << " - " << endTS << " \t" << (endTS - startTS) * 1.e-3 << std::endl;
        api.storeAsTFileAny(&chunk, baseCCDBpath + "ZorroHelpers", metadata, startTS, endTS + 1);
        startTS = endTS + 1;
        helperIndex += chunk.size();
      }
    }
  }
}

void uploadOTSobjects(std::string periodName)
{
  int year = 2000 + std::stoi(periodName.substr(3, 2));
  gSystem->Exec(Form("alien_find /alice/data/%i/%s/ ctf_skim_full/AnalysisResults_fullrun.root | sed 's:/AnalysisResults_fullrun\\.root::' > list_%s.txt", year, periodName.data(), periodName.data()));
  uploadOTSobjects(Form("list_%s.txt", periodName.data()), "", true, true);
}
