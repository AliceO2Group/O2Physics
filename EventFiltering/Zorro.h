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
// Zero Obstacles Results Retriever for Offline trigger selections

#ifndef EVENTFILTERING_ZORRO_H_
#define EVENTFILTERING_ZORRO_H_

#include <bitset>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "TH1D.h"
#include "TH2D.h"
#include "CommonDataFormat/IRFrame.h"
#include "Framework/HistogramRegistry.h"
#include "ZorroHelper.h"
#include "ZorroSummary.h"

namespace o2
{
namespace ccdb
{
class BasicCCDBManager;
};
}; // namespace o2

class Zorro
{
 public:
  Zorro() = default;
  std::vector<int> initCCDB(o2::ccdb::BasicCCDBManager* ccdb, int runNumber, uint64_t timestamp, std::string tois, int bcTolerance = 500);
  std::bitset<128> fetch(uint64_t bcGlobalId, uint64_t tolerance = 100);
  bool isSelected(uint64_t bcGlobalId, uint64_t tolerance = 100, TH2* toiHisto = nullptr);
  bool isNotSelectedByAny(uint64_t bcGlobalId, uint64_t tolerance = 100);

  void populateHistRegistry(o2::framework::HistogramRegistry& histRegistry, int runNumber, std::string folderName = "Zorro");
  void populateExternalHists(int runNumber, TH2* zorroHisto = nullptr, TH2* toiHisto = nullptr);

  TH1D* getScalers() const { return mScalers; }
  TH1D* getSelections() const { return mSelections; }
  TH1D* getInspectedTVX() const { return mInspectedTVX; }
  std::bitset<128> getLastResult() const { return mLastResult; }
  std::vector<int> getTOIcounters() const { return mTOIcounts; }
  std::vector<bool> getTriggerOfInterestResults(uint64_t bcGlobalId, uint64_t tolerance = 100);
  std::vector<bool> getTriggerOfInterestResults() const;

  void setCCDBpath(std::string path) { mBaseCCDBPath = path; }
  void setBaseCCDBPath(std::string path) { mBaseCCDBPath = path; }
  void setBCtolerance(int tolerance) { mBCtolerance = tolerance; }

  ZorroSummary* getZorroSummary() { return &mZorroSummary; }

 private:
  void setupHelpers(int64_t timestamp);

  ZorroSummary mZorroSummary{"ZorroSummary", "ZorroSummary"};

  std::string mBaseCCDBPath = "Users/m/mpuccio/EventFiltering/OTS/Chunked/";
  int mRunNumber = 0;
  std::pair<int64_t, int64_t> mRunDuration;
  int64_t mOrbitResetTimestamp = 0;
  TH1* mAnalysedTriggers;           /// Accounting for all triggers in the current run
  TH1* mAnalysedTriggersOfInterest; /// Accounting for triggers of interest in the current run

  std::vector<int> mRunNumberHistos;
  std::vector<TH1*> mAnalysedTriggersList;           /// Per run histograms
  std::vector<TH1*> mAnalysedTriggersOfInterestList; /// Per run histograms

  int mBCtolerance = 100;
  uint64_t mLastBCglobalId = 0;
  uint64_t mLastSelectedIdx = 0;
  TH1D* mScalers = nullptr;
  TH1D* mSelections = nullptr;
  TH1D* mInspectedTVX = nullptr;
  std::bitset<128> mLastResult;
  std::vector<bool> mAccountedBCranges; /// Avoid double accounting of inspected BC ranges
  std::vector<o2::dataformats::IRFrame> mBCranges;
  std::vector<ZorroHelper>* mZorroHelpers = nullptr;
  std::vector<std::string> mTOIs;
  std::vector<int> mTOIidx;
  std::vector<int> mTOIcounts;
  o2::ccdb::BasicCCDBManager* mCCDB = nullptr;
};

#endif // EVENTFILTERING_ZORRO_H_
