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
///
/// \file eventMixingValidation.cxx
/// \brief Validation tasks for event mixing.
/// \author Karwowska Maja
/// \since

#include <vector>
#include "Framework/AnalysisDataModel.h"
#include "Framework/SliceCache.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct MixedEventsEmptyTables {
  SliceCache cache;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  // Dummy filter to enforce empty tables
  Filter trackFilter = (aod::track::x > -0.8f) && (aod::track::x < -0.8f);
  using MyTracks = soa::Filtered<aod::Tracks>;

  std::vector<double> xBins{VARIABLE_WIDTH, -0.064, -0.062, -0.060, 0.066, 0.068, 0.070, 0.072};
  std::vector<double> yBins{VARIABLE_WIDTH, -0.320, -0.301, -0.300, 0.330, 0.340, 0.350, 0.360};
  using BinningType = ColumnBinningPolicy<aod::collision::PosX, aod::collision::PosY>;
  BinningType binningOnPositions{{xBins, yBins}, true};                                 // true is for 'ignore overflows' (true by default)
  SameKindPair<aod::Collisions, MyTracks, BinningType> pair{binningOnPositions, 5, -1, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

  void process(aod::Collisions const& collisions, MyTracks const& tracks)
  {
    LOGF(info, "Input data Collisions %d, Tracks %d ", collisions.size(), tracks.size());

    int count = 0;
    for (const auto& [c1, tracks1, c2, tracks2] : pair) {
      LOGF(info, "Mixed event collisions: (%d, %d)", c1.globalIndex(), c2.globalIndex());
      count++;
      if (count == 10)
        break;

      // Example of using tracks from mixed events -- iterate over all track pairs from the two collisions
      int trackCount = 0;
      for (const auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
        LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d)", t1.index(), t2.index(), c1.index(), c2.index(), t1.collision().index(), t2.collision().index());
        trackCount++;
        if (trackCount == 10)
          break;
      }
    }
  }
};

struct MixedEventsJoinedTracks {
  SliceCache cache;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  std::vector<double> xBins{VARIABLE_WIDTH, -0.064, -0.062, -0.060, 0.066, 0.068, 0.070, 0.072};
  std::vector<double> yBins{VARIABLE_WIDTH, -0.320, -0.301, -0.300, 0.330, 0.340, 0.350, 0.360};
  using BinningType = ColumnBinningPolicy<aod::collision::PosX, aod::collision::PosY>;
  BinningType binningOnPositions{{xBins, yBins}, true};                                        // true is for 'ignore overflows' (true by default)
  SameKindPair<aod::Collisions, aod::FullTracks, BinningType> pair{binningOnPositions, 5, -1, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

  void process(aod::Collisions const& collisions, aod::FullTracks const& tracks)
  {
    LOGF(info, "Input data Collisions %d, Tracks %d ", collisions.size(), tracks.size());

    int count = 0;
    for (const auto& [c1, tracks1, c2, tracks2] : pair) {
      LOGF(info, "Mixed event collisions: (%d, %d)", c1.globalIndex(), c2.globalIndex());
      count++;
      if (count == 100)
        break;

      // Example of using tracks from mixed events -- iterate over all track pairs from the two collisions
      int trackCount = 0;
      for (const auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
        LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d)", t1.index(), t2.index(), c1.index(), c2.index(), t1.collision().index(), t2.collision().index());
        trackCount++;
        if (trackCount == 10)
          break;
      }
    }
  }
};

// It should not compile
// struct MixedEventsBadSubscription {
//  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
//  std::vector<double> xBins{VARIABLE_WIDTH, -0.064, -0.062, -0.060, 0.066, 0.068, 0.070, 0.072};
//  std::vector<double> yBins{VARIABLE_WIDTH, -0.320, -0.301, -0.300, 0.330, 0.340, 0.350, 0.360};
//  using BinningType = ColumnBinningPolicy<aod::collision::PosX, aod::collision::PosY>;
//  BinningType binningOnPositions{{xBins, yBins}, true}; // true is for 'ignore overflows' (true by default)
//  SameKindPair<aod::Collisions, aod::Cascades, BinningType> pair{binningOnPositions, 5, -1}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored
//
//  void process(aod::Collisions const& collisions, aod::Tracks const& tracks)
//  {
//    int count = 0;
//    for (const auto& [c1, tracks1, c2, tracks2] : pair) {
//      LOGF(info, "Mixed event collisions: (%d, %d)", c1.globalIndex(), c2.globalIndex());
//      count++;
//      if (count == 10)
//        break;
//
//      // Example of using tracks from mixed events -- iterate over all track pairs from the two collisions
//      int trackCount = 0;
//      for (const auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
//        LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d)", t1.index(), t2.index(), c1.index(), c2.index(), t1.collision().index(), t2.collision().index());
//        trackCount++;
//        if (trackCount == 10)
//          break;
//      }
//    }
//  }
//};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MixedEventsEmptyTables>(cfgc),
    adaptAnalysisTask<MixedEventsJoinedTracks>(cfgc),
    // adaptAnalysisTask<MixedEventsBadSubscription>(cfgc), // Should not compile
  };
}
