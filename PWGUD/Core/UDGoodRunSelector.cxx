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

#include "PWGUD/Core/UDGoodRunSelector.h"

#include "Framework/Logger.h"

#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"

#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

class TFile;

using namespace rapidjson;

UDGoodRunSelector::UDGoodRunSelector(std::string const& goodRunsFile)
{
  init(goodRunsFile);
}

UDGoodRunSelector::~UDGoodRunSelector()
{
  clear();
}

void UDGoodRunSelector::clear()
{
  mgoodRuns.clear();
  mrunMap.clear();
  mrnMin = -1;
  mrnMax = -1;
  misActive = false;
}

void UDGoodRunSelector::Print()
{
  LOGF(info, "RunSelector: isActive %i", misActive);
  for (const auto& [period, runNumbers] : mrunMap) {
    LOGF(info, "  run numbers of period \"%s\"", period);
    for (const auto& runNumber : runNumbers) {
      LOGF(info, "    %i", runNumber);
    }
  }
  LOGF(info, "  all run numbers");
  for (const auto& runNumber : mgoodRuns) {
    LOGF(info, "    %i", runNumber);
  }
}

bool UDGoodRunSelector::isGoodRun(int runNumber)
{
  // search for runNumber in mgoodRuns
  if (!misActive) {
    return true;
  } else {
    return std::find(mgoodRuns.begin(), mgoodRuns.end(), runNumber) != mgoodRuns.end();
  }
}

std::vector<int> UDGoodRunSelector::goodRuns(std::string runPeriod)
{
  auto it = mrunMap.find(runPeriod.c_str());
  if (it != mrunMap.end()) {
    return it->second;
  } else {
    return std::vector<int>{};
  }
}

bool UDGoodRunSelector::init(std::string const& goodRunsFile)
{
  // parse goodRunsFile and fill mgoodRuns
  // goodRunsFile is expected to be a json file with the structure like
  //  {
  //    "goodRuns": [
  //      {
  //        "period": "LHC22e",
  //        "runlist":  [ 519041, 519043, 519045 ]
  //      },
  //      {
  //        "period": "LHC22f",
  //        "runlist":  [ 520143 ]
  //      }
  //    ]
  //  }
  misActive = false;
  if (goodRunsFile.empty()) {
    LOGF(info, "goodRuns was not specified!");
    return true;
  }

  // open the file
  FILE* fjson = fopen(goodRunsFile.c_str(), "r");
  if (!fjson) {
    LOGF(info, "Could not open goodRuns file %s", goodRunsFile);
    return false;
  }

  // create streamer
  char readBuffer[65536];
  FileReadStream jsonStream(fjson, readBuffer, sizeof(readBuffer));

  // parse the json file
  Document jsonDocument;
  jsonDocument.ParseStream(jsonStream);

  // is it a proper json document?
  if (jsonDocument.HasParseError()) {
    LOGF(error, "Check the goodRuns file! There is a problem with the format!");
    return false;
  }

  // goodRuns
  std::string runPeriod{""};
  int runNumber;

  const char* itemName = "goodRuns";
  if (!jsonDocument.HasMember(itemName)) {
    LOGF(error, "Check the goodRuns file! Item %s is missing!", itemName);
    return false;
  }

  const Value& item0 = jsonDocument[itemName];
  if (!item0.IsArray()) {
    LOGF(error, "Check the goodRuns file! Item %s must be an array!", itemName);
    return false;
  }

  // loop over all runPeriods
  for (auto& item1 : item0.GetArray()) {
    if (!item1.IsObject()) {
      LOGF(error, "Check the goodRuns file! %s must be objects!", itemName);
      return false;
    }
    // name of the period
    itemName = "period";
    if (item1.HasMember(itemName)) {
      if (item1[itemName].IsString()) {
        runPeriod = item1[itemName].GetString();
      } else {
        LOGF(error, "Check the goodRuns file!%s must be a string!", itemName);
        return false;
      }
    }
    // runNumbers
    itemName = "runlist";
    if (item1.HasMember(itemName)) {
      if (!item1[itemName].IsArray()) {
        LOGF(error, "Check the goodRuns file! Item %s must be an array!", itemName);
        return false;
      }
      // update goodRuns and mrunMap
      for (auto& item2 : item1[itemName].GetArray()) {
        runNumber = item2.GetInt();
        if (runNumber < mrnMin) {
          mrnMin = runNumber;
        } else if (runNumber > mrnMax) {
          mrnMax = runNumber;
        }
        mgoodRuns.push_back(runNumber);
        mrunMap[runPeriod].push_back(runNumber);
        misActive = true;
      }
    }
  }
  // make mgoodRuns unique
  std::sort(mgoodRuns.begin(), mgoodRuns.end());
  auto last = std::unique(mgoodRuns.begin(), mgoodRuns.end());
  mgoodRuns.erase(last, mgoodRuns.end());

  // clean up
  fclose(fjson);

  return misActive;
}

// =============================================================================
