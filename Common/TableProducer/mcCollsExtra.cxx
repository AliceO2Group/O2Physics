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
// Quick and dirty task to correlate MC <-> data
//

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <numeric>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using FullCollisions = soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As, aod::FT0Mults>;

// simple checkers
#define bitset(var, nbit) ((var) |= (1 << (nbit)))
#define bitcheck(var, nbit) ((var) & (1 << (nbit)))

struct mcCollisionExtra {
  Produces<aod::McCollsExtra> mcCollsExtra;
  Produces<aod::McCollContexts> mcCollContexts; // collision context with respect to its neighbours

  Configurable<int> pdgCodeOfInterest{"pdgCodeOfInterest", 3312, "PDG Code of interest"};
  Configurable<bool> pdgCodeAbsolute{"pdgCodeAbsolute", true, "if true, accept +/- pdgCodeOfInterest"};
  Configurable<float> poiEtaWindow{"poiEtaWindow", 0.8, "PDG code requirement within this eta window"};

  // For manual sliceBy
  Preslice<aod::McParticle> perMcCollision = aod::mcparticle::mcCollisionId;

  template <typename T>
  std::vector<std::size_t> sort_indices(const std::vector<T>& v)
  {
    std::vector<std::size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(),
                     [&v](std::size_t i1, std::size_t i2) { return v[i1] < v[i2]; });
    return idx;
  }

  void processNoCentrality(aod::McCollision const&, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions)
  {
    int biggestNContribs = -1;
    int bestCollisionIndex = -1;
    for (auto& collision : collisions) {
      if (biggestNContribs < collision.numContrib()) {
        biggestNContribs = collision.numContrib();
        bestCollisionIndex = collision.globalIndex();
      }
    }
    mcCollsExtra(collisions.size(), bestCollisionIndex, 0.0f);
  }
  void processWithCentrality(aod::McCollision const&, soa::SmallGroups<FullCollisions> const& collisions)
  {
    int biggestNContribs = -1;
    int bestCollisionIndex = -1;
    float bestCollisionCentFT0C = 100.5f;
    for (auto& collision : collisions) {
      if (biggestNContribs < collision.numContrib()) {
        biggestNContribs = collision.numContrib();
        bestCollisionIndex = collision.globalIndex();
        bestCollisionCentFT0C = collision.centFT0C();
      }
    }
    mcCollsExtra(collisions.size(), bestCollisionIndex, bestCollisionCentFT0C);
  }
  void processMcContexts(aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticlesUngrouped, FullCollisions const& collisions)
  {
    std::vector<float> mcCollisionTimes;
    std::vector<bool> mcCollisionHasPoI;
    for (auto& mcCollision : mcCollisions) {
      bool PoIpresent = false;
      auto mcParticles = mcParticlesUngrouped.sliceBy(perMcCollision, mcCollision.globalIndex());
      for (auto& mcParticle : mcParticles) {
        if (std::abs(mcParticle.eta()) < poiEtaWindow) {
          if (mcParticle.pdgCode() == pdgCodeOfInterest) {
            PoIpresent = true;
          }
          if (mcParticle.pdgCode() == -pdgCodeOfInterest && pdgCodeAbsolute) {
            PoIpresent = true;
          }
        }
      }
      mcCollisionTimes.emplace_back(mcCollision.t());
      mcCollisionHasPoI.emplace_back(PoIpresent);
    }
    // sort mcCollisions according to time
    auto sortedIndices = sort_indices(mcCollisionTimes);
    for (auto& collision : collisions) {
      uint16_t forwardHistory = 0, backwardHistory = 0;
      if (!collision.has_mcCollision()) {
        mcCollContexts(forwardHistory, backwardHistory);
        continue;
      }
      auto mcCollision = collision.mcCollision();
      auto iter = std::find(sortedIndices.begin(), sortedIndices.end(), mcCollision.index());
      if (iter != sortedIndices.end()) {
        auto index = std::distance(iter, sortedIndices.begin());
        for (auto iMcColl = index + 1; iMcColl < index + 17; iMcColl++) {
          if (iMcColl >= std::ssize(sortedIndices))
            continue;
          if (mcCollisionHasPoI[sortedIndices[iMcColl]])
            bitset(forwardHistory, iMcColl - index - 1);
        }
        for (int iMcColl = index - 1; iMcColl > index - 17; iMcColl--) {
          if (iMcColl <= 0)
            continue;
          if (mcCollisionHasPoI[sortedIndices[iMcColl]])
            bitset(backwardHistory, index + 1 - iMcColl);
        }
      }
      mcCollContexts(forwardHistory, backwardHistory);
    }
  }

  PROCESS_SWITCH(mcCollisionExtra, processNoCentrality, "process as if real data", false);
  PROCESS_SWITCH(mcCollisionExtra, processWithCentrality, "process as if real data", true);
  PROCESS_SWITCH(mcCollisionExtra, processMcContexts, "process mc collision contexts", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mcCollisionExtra>(cfgc)};
}
