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

const std::string kBaseCCDBPath = "Users/m/mpuccio/EventFiltering/OTS/";

#pragma link C++ class std::vector<std::array<uint64_t, 2>> + ;

void uploadOTSobjects(std::string inputList) {
  o2::ccdb::CcdbApi api;
  api.init("http://alice-ccdb.cern.ch");

  std::ifstream file(inputList.data());
  std::string path;
  std::map<std::string, std::string> metadata;
  while (std::getline(file, path)) {
    auto pos = path.find("/5") + 1;  /// in the path at some point there is the run number
    const std::string runString = path.substr(pos, 6);
    std::cout << "Processing run " << runString << std::endl;
    const int runNumber = std::stoi(runString);
    metadata["runNumber"] = runString;
    std::pair<int64_t, int64_t> duration = o2::ccdb::BasicCCDBManager::getRunDuration(api, runNumber);
    duration.first -= 10000; // subtract 3 minutes from the run start
    duration.second += 180000; // add 3 minutes to the run duration
    TFile scalersFile((path + "AnalysisResults_fullrun.root").data(), "READ");
    TH1* scalers = (TH1*)scalersFile.Get("central-event-filter-task/scalers/mScalers");
    TH1* filters = (TH1*)scalersFile.Get("central-event-filter-task/scalers/mFiltered");
    api.storeAsTFile(scalers, kBaseCCDBPath + "FilterCounters", metadata, duration.first, duration.second);
    api.storeAsTFile(filters, kBaseCCDBPath + "SelectionCounters", metadata, duration.first, duration.second);
    TH1* hCounterTVX = (TH1*)scalersFile.Get("bc-selection-task/hCounterTVX");
    api.storeAsTFile(hCounterTVX, kBaseCCDBPath + "InspectedTVX", metadata, duration.first, duration.second);

    std::vector<std::array<uint64_t, 2>> bcRanges, filterBitMask, selectionBitMask;
    TFile bcRangesFile((path + "bcRanges_fullrun.root").data(), "READ");
    TKey* key;
    auto klst = bcRangesFile.GetListOfKeys();
    TIter nextkey(klst);
    while ((key = (TKey*)nextkey())) {
      std::string kcl(key->GetClassName());
      if (kcl == "TDirectoryFile") {
        ROOT::RDataFrame bcRangesFrame(std::string(key->GetName()) + "/selectedBC", &bcRangesFile);
        auto bcAO2D = bcRangesFrame.Take<uint64_t>("bcAO2D");
        auto bcEvSel = bcRangesFrame.Take<uint64_t>("bcEvSel");
        auto selMask = bcRangesFrame.Take<uint64_t>("selMask");
        auto triMask = bcRangesFrame.Take<uint64_t>("triMask");
        for (size_t i = 0; i < bcAO2D->size(); i++) {
          bcRanges.push_back({bcAO2D->at(i), bcEvSel->at(i)});
          filterBitMask.push_back({bcAO2D->at(i), 0ull});
          selectionBitMask.push_back({bcAO2D->at(i), 0ull});
        }
      }
    }
    api.storeAsTFileAny(&bcRanges, kBaseCCDBPath + "SelectedBCs", metadata, duration.first, duration.second);
    api.storeAsTFileAny(&filterBitMask, kBaseCCDBPath + "FilterBitMasks", metadata, duration.first, duration.second);
    api.storeAsTFileAny(&selectionBitMask, kBaseCCDBPath + "SelectionBitMasks", metadata, duration.first, duration.second);
  }
}
