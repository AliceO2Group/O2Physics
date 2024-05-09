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

#include "Zorro.h"

#include <map>

#include "TH1D.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonDataFormat/InteractionRecord.h"

using o2::InteractionRecord;

std::vector<int> Zorro::initCCDB(o2::ccdb::BasicCCDBManager* ccdb, int runNumber, uint64_t timestamp, std::string tois, int bcRange)
{
  if (mRunNumber == runNumber) {
    return mTOIidx;
  }
  mCCDB = ccdb;
  mRunNumber = runNumber;
  mBCtolerance = bcRange;
  std::map<std::string, std::string> metadata;
  metadata["runNumber"] = std::to_string(runNumber);
  mScalers = mCCDB->getSpecific<TH1D>(mBaseCCDBPath + "FilterCounters", timestamp, metadata);
  mSelections = mCCDB->getSpecific<TH1D>(mBaseCCDBPath + "SelectionCounters", timestamp, metadata);
  mInspectedTVX = mCCDB->getSpecific<TH1D>(mBaseCCDBPath + "InspectedTVX", timestamp, metadata);
  auto selectedBCs = mCCDB->getSpecific<std::vector<std::array<uint64_t, 2>>>(mBaseCCDBPath + "SelectedBCs", timestamp, metadata);
  mBCranges.clear();
  for (auto bc : *selectedBCs) {
    mBCranges.emplace_back(InteractionRecord::long2IR(std::min(bc[0], bc[1])), InteractionRecord::long2IR(std::max(bc[0], bc[1])));
  }
  std::sort(mBCranges.begin(), mBCranges.end(), [](const auto& a, const auto& b) { return a.getMin() < b.getMin(); });

  mSelectionBitMask = mCCDB->getSpecific<std::vector<std::array<uint64_t, 2>>>(mBaseCCDBPath + "SelectionBitMask", timestamp, metadata);
  mFilterBitMask = mCCDB->getSpecific<std::vector<std::array<uint64_t, 2>>>(mBaseCCDBPath + "FilterBitMask", timestamp, metadata);

  mLastBCglobalId = 0;
  mLastSelectedIdx = 0;
  mTOIs.clear();
  mTOIidx.clear();
  size_t pos = 0;
  while ((pos = tois.find(",")) != std::string::npos) {
    std::string token = tois.substr(0, pos);
    // Trim leading and trailing whitespaces from the token
    token.erase(0, token.find_first_not_of(" "));
    token.erase(token.find_last_not_of(" ") + 1);
    int bin = mScalers->GetXaxis()->FindBin(token.c_str()) - 2;
    mTOIs.push_back(token);
    mTOIidx.push_back(bin);
    tois.erase(0, pos + 1);
  }
  mTOIcounts.resize(mTOIs.size(), 0);
  return mTOIidx;
}

std::bitset<128> Zorro::fetch(uint64_t bcGlobalId, uint64_t tolerance)
{
  std::bitset<128> result;
  o2::dataformats::IRFrame bcFrame{InteractionRecord::long2IR(bcGlobalId) - tolerance, InteractionRecord::long2IR(bcGlobalId) + tolerance};
  if (bcGlobalId < mLastBCglobalId) {
    mLastSelectedIdx = 0;
  }
  for (size_t i = mLastSelectedIdx; i < mBCranges.size(); i++) {
    if (!bcFrame.getOverlap(mBCranges[i]).isZeroLength()) {
      for (int iMask{0}; iMask < 2; ++iMask) {
        for (int iTOI{0}; iTOI < 64; ++iTOI) {
          result.set(iMask * 64 + iTOI, mFilterBitMask->at(i)[iMask] & (1ull << iTOI));
        }
      }
      mLastSelectedIdx = i;
      return result;
    }
  }
  return result;
}

bool Zorro::isSelected(uint64_t bcGlobalId, uint64_t tolerance)
{
  int lastSelectedIdx = mLastSelectedIdx;
  std::bitset<128> result = fetch(bcGlobalId, tolerance);
  for (size_t i{0}; i < mTOIidx.size(); ++i) {
    if (mTOIidx[i] < 0) {
      continue;
    } else if (result.test(mTOIidx[i])) {
      mTOIcounts[i] += (lastSelectedIdx != mLastSelectedIdx); /// Avoid double counting
      return true;
    }
  }
  return false;
}