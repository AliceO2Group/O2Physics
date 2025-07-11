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

#include "ZorroSummary.h"

#include <TCollection.h>
#include <TObject.h>

#include <RtypesCore.h>

#include <cstddef>

void ZorroSummary::Copy(TObject& c) const
{
  static_cast<ZorroSummary&>(c) = *this;
}

Long64_t ZorroSummary::Merge(TCollection* list)
{
  if (!list) {
    return 0;
  }
  int n = 1;
  if (list->IsEmpty()) {
    return n;
  }

  for (auto* obj : *list) {
    auto* entry = dynamic_cast<ZorroSummary*>(obj);
    if (!entry || entry->getTOInames() != mTOInames) {
      continue;
    }
    n++;
    auto& analysedToiCounters = entry->getAnalysedTOIcounters();
    for (const auto& [runNumber, currentAnalysedToiCounters] : analysedToiCounters) {
      if (mAnalysedTOIcounters.find(runNumber) == mAnalysedTOIcounters.end()) {
        mAnalysedTOIcounters[runNumber] = currentAnalysedToiCounters;
        mTVXcounters[runNumber] = entry->getTVXcounters().at(runNumber);
        mTOIcounters[runNumber] = entry->getTOIcounters().at(runNumber);
      } else {
        auto& thisCounters = mAnalysedTOIcounters[runNumber];
        for (std::size_t i = 0; i < thisCounters.size(); ++i) {
          thisCounters[i] += currentAnalysedToiCounters[i];
        }
      }
    }
  }
  return n;
}

double ZorroSummary::getNormalisationFactor(int toiId) const
{
  double totalTOI{0.}, totalTVX{0.};
  ULong64_t totalAnalysedTOI{0};
  for (const auto& [runNumber, toiCounters] : mTOIcounters) {
    totalTOI += toiCounters.at(toiId);
  }
  for (const auto& [runNumber, tvxCounters] : mTVXcounters) {
    totalTVX += tvxCounters;
  }
  for (const auto& [runNumber, analysedTOIcounters] : mAnalysedTOIcounters) {
    totalAnalysedTOI += analysedTOIcounters.at(toiId);
  }

  return totalTVX * totalAnalysedTOI / totalTOI;
}
