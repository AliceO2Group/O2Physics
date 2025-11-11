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

#include "Common/Core/ZorroHelper.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <CommonDataFormat/IRFrame.h>
#include <CommonDataFormat/InteractionRecord.h>
#include <CommonUtils/StringUtils.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>

#include <TH1.h>
#include <TH2.h>
#include <TString.h>

#include <RtypesCore.h>

#include <algorithm>
#include <bitset>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

using o2::InteractionRecord;

namespace
{
int findBin(TH1* hist, const std::string& label)
{ // Find bin by label, avoiding the axis extention from the native ROOT implementation
  for (int iBin{1}; iBin <= hist->GetNbinsX(); ++iBin) {
    if (label == hist->GetXaxis()->GetBinLabel(iBin)) {
      return iBin;
    }
  }
  return -1;
}
} // namespace

void Zorro::populateHistRegistry(o2::framework::HistogramRegistry& histRegistry, int runNumber, std::string folderName)
{
  int runId{-1};
  for (size_t i{0}; i < mRunNumberHistos.size(); ++i) {
    if (mRunNumberHistos[i] == runNumber) {
      runId = i;
      break;
    }
  }
  if (runId > -1) {
    /// Support jobs running on non-continuous run numbers
    mAnalysedTriggers = mAnalysedTriggersList[runId];
    mAnalysedTriggersOfInterest = mAnalysedTriggersOfInterestList[runId];
    return;
  }
  if (mSelections) {
    mAnalysedTriggers = histRegistry.add<TH1>((folderName + "/" + std::to_string(runNumber) + "/" + "AnalysedTriggers").data(), "", o2::framework::HistType::kTH1D, {{mSelections->GetNbinsX() - 2, -0.5, mSelections->GetNbinsX() - 2.5}}).get();
    for (int iBin{2}; iBin < mSelections->GetNbinsX(); ++iBin) { // Exclude first and last bins as they are total number of analysed and selected events, respectively
      mAnalysedTriggers->GetXaxis()->SetBinLabel(iBin - 1, mSelections->GetXaxis()->GetBinLabel(iBin));
    }
    std::shared_ptr<TH1> selections = histRegistry.add<TH1>((folderName + "/" + std::to_string(runNumber) + "/" + "Selections").data(), "", o2::framework::HistType::kTH1D, {{mSelections->GetNbinsX(), -0.5, static_cast<double>(mSelections->GetNbinsX() - 0.5)}});
    selections->SetBit(TH1::kIsAverage);
    for (int iBin{1}; iBin <= mSelections->GetNbinsX(); ++iBin) {
      selections->GetXaxis()->SetBinLabel(iBin, mSelections->GetXaxis()->GetBinLabel(iBin));
      selections->SetBinContent(iBin, mSelections->GetBinContent(iBin));
      selections->SetBinError(iBin, mSelections->GetBinError(iBin));
    }
  }
  if (mScalers) {
    std::shared_ptr<TH1> scalers = histRegistry.add<TH1>((folderName + "/" + std::to_string(runNumber) + "/" + "Scalers").data(), "", o2::framework::HistType::kTH1D, {{mScalers->GetNbinsX(), -0.5, static_cast<double>(mScalers->GetNbinsX() - 0.5)}});
    scalers->SetBit(TH1::kIsAverage);
    for (int iBin{1}; iBin <= mScalers->GetNbinsX(); ++iBin) {
      scalers->GetXaxis()->SetBinLabel(iBin, mScalers->GetXaxis()->GetBinLabel(iBin));
      scalers->SetBinContent(iBin, mScalers->GetBinContent(iBin));
      scalers->SetBinError(iBin, mScalers->GetBinError(iBin));
    }
  }
  if (mInspectedTVX) {
    std::shared_ptr<TH1> inspectedTVX = histRegistry.add<TH1>((folderName + "/" + std::to_string(runNumber) + "/" + "InspectedTVX").data(), "", o2::framework::HistType::kTH1D, {{mInspectedTVX->GetNbinsX(), -0.5, static_cast<double>(mInspectedTVX->GetNbinsX() - 0.5)}});
    inspectedTVX->SetBit(TH1::kIsAverage);
    for (int iBin{1}; iBin <= mInspectedTVX->GetNbinsX(); ++iBin) {
      inspectedTVX->GetXaxis()->SetBinLabel(iBin, mInspectedTVX->GetXaxis()->GetBinLabel(iBin));
      inspectedTVX->SetBinContent(iBin, mInspectedTVX->GetBinContent(iBin));
      inspectedTVX->SetBinError(iBin, mInspectedTVX->GetBinError(iBin));
    }
  }
  if (mTOIs.size()) {
    mAnalysedTriggersOfInterest = histRegistry.add<TH1>((folderName + "/" + std::to_string(runNumber) + "/" + "AnalysedTriggersOfInterest").data(), "", o2::framework::HistType::kTH1D, {{static_cast<int>(mTOIs.size()), -0.5, static_cast<double>(mTOIs.size() - 0.5)}}).get();
    for (size_t i{0}; i < mTOIs.size(); ++i) {
      mAnalysedTriggersOfInterest->GetXaxis()->SetBinLabel(i + 1, mTOIs[i].data());
    }
  }
  mAnalysedTriggersList.push_back(mAnalysedTriggers);
  mAnalysedTriggersOfInterestList.push_back(mAnalysedTriggersOfInterest);
  mRunNumberHistos.push_back(runNumber);
}

void Zorro::populateExternalHists(int runNumber, TH2* ZorroHisto, TH2* ToiHisto)
{
  // x-axis is run number, y-axis is same as ZorroSummary
  int runId{-1};
  for (size_t i{0}; i < mRunNumberHistos.size(); ++i) {
    if (mRunNumberHistos[i] == runNumber) {
      runId = i;
      break;
    }
  }
  if (runId > -1) {
    return;
  }
  // if the summary histogram is not set, create a new one
  if (!ZorroHisto) {
    LOGF(info, "Summary histogram not set, creating a new one");
    ZorroHisto = new TH2D("Zorro", "Zorro", 1, -0.5, 0.5, 1 + mTOIs.size() * 2, -0.5, mTOIs.size() * 2 - 0.5);
    ZorroHisto->SetBit(TH1::kIsAverage);
  }
  if (!ToiHisto) {
    LOGF(info, "TOI histogram not set, creating a new one");
    ToiHisto = new TH2D("TOI", "TOI", 1, -0.5, 0.5, mTOIs.size() * 2, -0.5, mTOIs.size() * 2 - 0.5);
  }
  // if it is the first run, initialize the histogram
  if (mRunNumberHistos.size() == 0) {
    ZorroHisto->SetBins(1, -0.5, 0.5, 1 + mTOIs.size() * 2, -0.5, mTOIs.size() * 2 - 0.5);
    ZorroHisto->SetBit(TH1::kIsAverage);
    ZorroHisto->GetXaxis()->SetBinLabel(1, Form("%d", runNumber));
    ZorroHisto->GetYaxis()->SetBinLabel(1, "inspected TVX");
    for (size_t i{0}; i < mTOIs.size(); ++i) {
      ZorroHisto->GetYaxis()->SetBinLabel(i + 2, Form("%s selections", mTOIs[i].data()));
      ZorroHisto->GetYaxis()->SetBinLabel(i + 2 + mTOIs.size(), Form("%s scalers", mTOIs[i].data()));
    }
    // TOI histogram
    ToiHisto->SetBins(1, -0.5, 0.5, mTOIs.size() * 2, -0.5, mTOIs.size() * 2 - 0.5);
    ToiHisto->GetXaxis()->SetBinLabel(1, Form("%d", runNumber));
    for (size_t i{0}; i < mTOIs.size(); ++i) {
      ToiHisto->GetYaxis()->SetBinLabel(i * 2 + 1, mTOIs[i].data());
      ToiHisto->GetYaxis()->SetBinLabel(i * 2 + 2, Form("%s AnalysedTriggers", mTOIs[i].data()));
    }
  }
  if (mInspectedTVX) {
    ZorroHisto->Fill(Form("%d", runNumber), "inspected TVX", mInspectedTVX->GetBinContent(1));
    ZorroHisto->SetBinError(mRunNumberHistos.size() + 1, 1, mInspectedTVX->GetBinError(1));
  }
  if (mSelections) {
    mAnalysedTriggers = new TH1D("AnalysedTriggers", "", mSelections->GetNbinsX() - 2, -0.5, mSelections->GetNbinsX() - 2.5);
    for (int iBin{2}; iBin < mSelections->GetNbinsX(); ++iBin) { // Exclude first and last bins as they are total number of analysed and selected events, respectively
      mAnalysedTriggers->GetXaxis()->SetBinLabel(iBin - 1, mSelections->GetXaxis()->GetBinLabel(iBin));
    }
    for (size_t i{0}; i < mTOIs.size(); ++i) {
      int bin = findBin(mSelections, mTOIs[i]);
      ZorroHisto->Fill(Form("%d", runNumber), Form("%s selections", mTOIs[i].data()), mSelections->GetBinContent(bin));
      ZorroHisto->SetBinError(mRunNumberHistos.size() + 1, i + 2, mSelections->GetBinError(bin));
    }
  }
  if (mScalers) {
    for (size_t i{0}; i < mTOIs.size(); ++i) {
      int bin = findBin(mScalers, mTOIs[i]);
      ZorroHisto->Fill(Form("%d", runNumber), Form("%s scalers", mTOIs[i].data()), mScalers->GetBinContent(bin));
      ZorroHisto->SetBinError(mRunNumberHistos.size() + 1, i + 2 + mTOIs.size(), mScalers->GetBinError(bin));
    }
  }

  mRunNumberHistos.push_back(runNumber);
}

std::vector<int> Zorro::initCCDB(o2::ccdb::BasicCCDBManager* ccdb, int runNumber, uint64_t timestamp, std::string tois, int bcRange)
{
  if (mRunNumber == runNumber) {
    return mTOIidx;
  }
  mCCDB = ccdb;
  mRunNumber = runNumber;
  mBCtolerance = bcRange;
  auto ctp = ccdb->getForRun<std::vector<Long64_t>>("CTP/Calib/OrbitReset", runNumber, false);
  mOrbitResetTimestamp = (*ctp)[0];
  mScalers = mCCDB->getForRun<TH1D>(mBaseCCDBPath + "FilterCounters", runNumber, true);
  mSelections = mCCDB->getForRun<TH1D>(mBaseCCDBPath + "SelectionCounters", runNumber, true);
  mInspectedTVX = mCCDB->getForRun<TH1D>(mBaseCCDBPath + "InspectedTVX", runNumber, true);
  setupHelpers(timestamp);
  mLastBCglobalId = 0;
  mLastSelectedIdx = 0;
  mTOIs.clear();
  mTOIidx.clear();
  std::vector<std::string> tokens = o2::utils::Str::tokenize(tois, ','); // tokens are trimmed
  for (auto const& token : tokens) {
    int bin = findBin(mSelections, token) - 2;
    mTOIs.push_back(token);
    mTOIidx.push_back(bin);
  }
  mTOIcounts.resize(mTOIs.size(), 0);
  mATcounts.resize(mSelections->GetNbinsX() - 2, 0);
  LOGF(info, "Zorro initialized for run %d, triggers of interest:", runNumber);
  for (size_t i{0}; i < mTOIs.size(); ++i) {
    LOGF(info, ">>> %s : %i", mTOIs[i].data(), mTOIidx[i]);
  }
  mZorroSummary.setupTOIs(mTOIs.size(), mTOIs);
  std::vector<double> toiCounters(mTOIs.size(), 0.);
  for (size_t i{0}; i < mTOIs.size(); ++i) {
    toiCounters[i] = mSelections->GetBinContent(mTOIidx[i] + 2);
  }
  mZorroSummary.setupRun(runNumber, mInspectedTVX->GetBinContent(1), toiCounters);

  return mTOIidx;
}

std::bitset<128> Zorro::fetch(uint64_t bcGlobalId, uint64_t tolerance)
{
  mLastResult.reset();
  if (bcGlobalId < mBCranges.front().getMin().toLong() - tolerance || bcGlobalId > mBCranges.back().getMax().toLong() + tolerance) {
    setupHelpers((mOrbitResetTimestamp + static_cast<int64_t>(bcGlobalId * o2::constants::lhc::LHCBunchSpacingNS * 1e-3)) / 1000);
  }

  o2::dataformats::IRFrame bcFrame{InteractionRecord::long2IR(bcGlobalId) - tolerance, InteractionRecord::long2IR(bcGlobalId) + tolerance};
  if (bcGlobalId < mLastBCglobalId) { /// Handle the possible discontinuity in the BC processed by the analyses
    mLastSelectedIdx = 0;
  }
  uint64_t lastSelectedIdx = mLastSelectedIdx;
  mLastBCglobalId = bcGlobalId;
  for (size_t i = mLastSelectedIdx; i < mBCranges.size(); i++) {
    if (!mBCranges[i].isOutside(bcFrame)) {
      for (int iMask{0}; iMask < 2; ++iMask) {
        for (int iTOI{0}; iTOI < 64; ++iTOI) {
          if (mZorroHelpers->at(i).selMask[iMask] & (1ull << iTOI)) {
            mLastResult.set(iMask * 64 + iTOI, 1);
            if (!mAccountedBCranges[i]) {
              mATcounts[iMask * 64 + iTOI]++;
              if (mAnalysedTriggers) {
                mAnalysedTriggers->Fill(iMask * 64 + iTOI);
              }
            }
          }
        }
      }
      mAccountedBCranges[i] = true;
      mLastSelectedIdx = mLastSelectedIdx == lastSelectedIdx-- ? i : mLastSelectedIdx; /// Decrease lastSelectedIdx to make sure this check is valid only in its first instance
    } else if (mBCranges[i].getMax() < bcFrame.getMin()) {
      mLastSelectedIdx = i;
    } else if (mBCranges[i].getMin() > bcFrame.getMax()) {
      break;
    }
  }
  return mLastResult;
}

bool Zorro::isSelected(uint64_t bcGlobalId, uint64_t tolerance, TH2* ToiHisto)
{
  uint64_t lastSelectedIdx = mLastSelectedIdx;
  fetch(bcGlobalId, tolerance);
  bool retVal{false};
  for (size_t i{0}; i < mTOIidx.size(); ++i) {
    if (mTOIidx[i] < 0) {
      continue;
    } else if (mLastResult.test(mTOIidx[i])) {
      if (ToiHisto && mAnalysedTriggers) {
        int binX = ToiHisto->GetXaxis()->FindBin(Form("%d", mRunNumber));
        int binY = ToiHisto->GetYaxis()->FindBin(Form("%s AnalysedTriggers", mTOIs[i].data()));
        ToiHisto->SetBinContent(binX, binY, mAnalysedTriggers->GetBinContent(mAnalysedTriggers->GetXaxis()->FindBin(mTOIs[i].data())));
      }
      mTOIcounts[i] += (lastSelectedIdx != mLastSelectedIdx); /// Avoid double counting
      if (mAnalysedTriggersOfInterest && lastSelectedIdx != mLastSelectedIdx) {
        mAnalysedTriggersOfInterest->Fill(i);
        mZorroSummary.increaseTOIcounter(mRunNumber, i);
      }
      if (ToiHisto && lastSelectedIdx != mLastSelectedIdx) {
        ToiHisto->Fill(Form("%d", mRunNumber), Form("%s", mTOIs[i].data()), 1);
      }
      retVal = true;
    }
  }
  return retVal;
}

std::vector<bool> Zorro::getTriggerOfInterestResults(uint64_t bcGlobalId, uint64_t tolerance)
{
  fetch(bcGlobalId, tolerance);
  return getTriggerOfInterestResults();
}

std::vector<bool> Zorro::getTriggerOfInterestResults() const
{
  std::vector<bool> results(mTOIidx.size(), false);
  for (size_t i{0}; i < mTOIidx.size(); ++i) {
    if (mTOIidx[i] < 0) {
      continue;
    } else if (mLastResult.test(mTOIidx[i])) {
      results[i] = true;
    }
  }
  return results;
}

bool Zorro::isNotSelectedByAny(uint64_t bcGlobalId, uint64_t tolerance)
{
  fetch(bcGlobalId, tolerance);
  return mLastResult.none();
}

void Zorro::setupHelpers(int64_t timestamp)
{
  if (mCCDB->isCachedObjectValid(mBaseCCDBPath + "ZorroHelpers", timestamp)) {
    return;
  }
  mZorroHelpers = mCCDB->getSpecific<std::vector<ZorroHelper>>(mBaseCCDBPath + "ZorroHelpers", timestamp, {{"runNumber", std::to_string(mRunNumber)}});
  std::sort(mZorroHelpers->begin(), mZorroHelpers->end(), [](const auto& a, const auto& b) { return std::min(a.bcAOD, a.bcEvSel) < std::min(b.bcAOD, b.bcEvSel); });
  mBCranges.clear();
  mAccountedBCranges.clear();
  for (const auto& helper : *mZorroHelpers) {
    mBCranges.emplace_back(InteractionRecord::long2IR(std::min(helper.bcAOD, helper.bcEvSel)), InteractionRecord::long2IR(std::max(helper.bcAOD, helper.bcEvSel)));
  }
  mAccountedBCranges.resize(mBCranges.size(), false);
}
