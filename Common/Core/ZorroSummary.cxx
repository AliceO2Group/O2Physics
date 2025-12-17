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

std::unordered_map<int, std::vector<double>> ZorroSummary::getPerRunNormalisationFactors() const
{
  std::unordered_map<int, std::vector<double>> runNormalisations;
  for (const auto& [runNumber, analysedTOIcounters] : mAnalysedTOIcounters) {
    const double tvxCount = mTVXcounters.at(runNumber);
    for (int toiId = 0; toiId < mNtois; ++toiId) {
      double toiCount = mTOIcounters.at(runNumber).at(toiId);
      ULong64_t analysedTOI = analysedTOIcounters.at(toiId);
      runNormalisations[runNumber].push_back(tvxCount * analysedTOI / toiCount);
    }
  }
  return runNormalisations;
}

double ZorroSummary::getNormalisationFactor(int toiId) const
{
  double normalisationFactor{0.};
  auto runNormalisations = getPerRunNormalisationFactors();
  for (const auto& [runNumber, normalisations] : runNormalisations) {
    normalisationFactor += normalisations.at(toiId);
  }
  return normalisationFactor;
}
