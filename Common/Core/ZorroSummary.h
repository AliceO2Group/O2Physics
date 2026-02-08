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

#ifndef EVENTFILTERING_ZORROSUMMARY_H_
#define EVENTFILTERING_ZORROSUMMARY_H_

#include <TNamed.h>

#include <Rtypes.h>
#include <RtypesCore.h>

#include <string>
#include <unordered_map>
#include <vector>

class ZorroSummary : public TNamed
{
 public:
  ZorroSummary() = default;
  ZorroSummary(const char* name, const char* objTitle) : TNamed(name, objTitle) {}
  virtual ~ZorroSummary() = default;   // NOLINT: Making this override breaks compilation for unknown reason
  virtual void Copy(TObject& c) const; // NOLINT: Making this override breaks compilation for unknown reason
  virtual Long64_t Merge(TCollection* list);

  void setupTOIs(int ntois, const std::vector<std::string>& toinames)
  {
    mNtois = ntois;
    if (toinames.size() == 0) {
      return;
    }
    mTOInames = toinames[0];
    for (size_t i = 1; i < toinames.size(); i++) {
      mTOInames += "," + toinames[i];
    }
  }
  void setupRun(int runNumber, double tvxCountes, const std::vector<double>& toiCounters)
  {
    if (mRunNumber == runNumber) {
      return;
    }
    mRunNumber = runNumber;
    mTVXcounters[runNumber] = tvxCountes;
    mTOIcounters[runNumber] = toiCounters;
    mAnalysedTOIcounters.try_emplace(runNumber, std::vector<ULong64_t>(mNtois, 0ull));
    mCurrentAnalysedTOIcounters = &mAnalysedTOIcounters[runNumber];
  }
  double getNormalisationFactor(int toiId) const;
  std::unordered_map<int, std::vector<double>> getPerRunNormalisationFactors() const;

  void increaseTOIcounter(int runNumber, int toiId)
  {
    if (runNumber != mRunNumber) {
      return;
    }
    mCurrentAnalysedTOIcounters->at(toiId)++;
  }

  const auto& getTOInames() const { return mTOInames; }
  const auto& getTOIcounters() const { return mTOIcounters; }
  const auto& getTVXcounters() const { return mTVXcounters; }
  const auto& getAnalysedTOIcounters() const { return mAnalysedTOIcounters; }

 private:
  int mRunNumber = 0;                                            //! Run currently being analysed
  std::vector<ULong64_t>* mCurrentAnalysedTOIcounters = nullptr; //! Analysed TOI counters for the current run

  int mNtois = 0;
  std::string mTOInames;
  std::unordered_map<int, std::vector<ULong64_t>> mAnalysedTOIcounters;
  std::unordered_map<int, std::vector<double>> mTOIcounters;
  std::unordered_map<int, double> mTVXcounters;

  ClassDef(ZorroSummary, 1);
};

#endif // EVENTFILTERING_ZORROSUMMARY_H_
