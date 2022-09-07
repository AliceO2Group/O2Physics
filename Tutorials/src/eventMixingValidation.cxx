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
/// \brief Validation tasks for event mixing.
/// \author
/// \since

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
  // Dummy filter to enforce empty tables
  Filter trackFilter = (aod::track::x > -0.8f) && (aod::track::x < -0.8f);
  using myTracks = soa::Filtered<aod::Tracks>;

  std::vector<double> xBins{VARIABLE_WIDTH, -0.064, -0.062, -0.060, 0.066, 0.068, 0.070, 0.072};
  std::vector<double> yBins{VARIABLE_WIDTH, -0.320, -0.301, -0.300, 0.330, 0.340, 0.350, 0.360};
  using BinningType = ColumnBinningPolicy<aod::collision::PosX, aod::collision::PosY>;
  BinningType binningOnPositions{{xBins, yBins}, true};                                 // true is for 'ignore overflows' (true by default)
  SameKindPair<aod::Collisions, myTracks, BinningType> pair{binningOnPositions, 5, -1}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

  void process(aod::Collisions const& collisions, myTracks const& tracks)
  {
    LOGF(info, "Input data Collisions %d, Tracks %d ", collisions.size(), tracks.size());

    int count = 0;
    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      LOGF(info, "Mixed event collisions: (%d, %d)", c1.globalIndex(), c2.globalIndex());
      count++;
      if (count == 10)
        break;

      // Example of using tracks from mixed events -- iterate over all track pairs from the two collisions
      int trackCount = 0;
      for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
        LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d)", t1.index(), t2.index(), c1.index(), c2.index(), t1.collision().index(), t2.collision().index());
        trackCount++;
        if (trackCount == 10)
          break;
      }
    }
  }
};

struct MixedEventsJoinedTracks {
  std::vector<double> xBins{VARIABLE_WIDTH, -0.064, -0.062, -0.060, 0.066, 0.068, 0.070, 0.072};
  std::vector<double> yBins{VARIABLE_WIDTH, -0.320, -0.301, -0.300, 0.330, 0.340, 0.350, 0.360};
  using BinningType = ColumnBinningPolicy<aod::collision::PosX, aod::collision::PosY>;
  BinningType binningOnPositions{{xBins, yBins}, true};                                        // true is for 'ignore overflows' (true by default)
  SameKindPair<aod::Collisions, aod::FullTracks, BinningType> pair{binningOnPositions, 5, -1}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

  void process(aod::Collisions const& collisions, aod::FullTracks const& tracks)
  {
    LOGF(info, "Input data Collisions %d, Tracks %d ", collisions.size(), tracks.size());

    int count = 0;
    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      LOGF(info, "Mixed event collisions: (%d, %d)", c1.globalIndex(), c2.globalIndex());
      count++;
      if (count == 100)
        break;

      // Example of using tracks from mixed events -- iterate over all track pairs from the two collisions
      int trackCount = 0;
      for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
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
//  std::vector<double> xBins{VARIABLE_WIDTH, -0.064, -0.062, -0.060, 0.066, 0.068, 0.070, 0.072};
//  std::vector<double> yBins{VARIABLE_WIDTH, -0.320, -0.301, -0.300, 0.330, 0.340, 0.350, 0.360};
//  using BinningType = ColumnBinningPolicy<aod::collision::PosX, aod::collision::PosY>;
//  BinningType binningOnPositions{{xBins, yBins}, true}; // true is for 'ignore overflows' (true by default)
//  SameKindPair<aod::Collisions, aod::Cascades, BinningType> pair{binningOnPositions, 5, -1}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored
//
//  void process(aod::Collisions const& collisions, aod::Tracks const& tracks)
//  {
//    int count = 0;
//    for (auto& [c1, tracks1, c2, tracks2] : pair) {
//      LOGF(info, "Mixed event collisions: (%d, %d)", c1.globalIndex(), c2.globalIndex());
//      count++;
//      if (count == 10)
//        break;
//
//      // Example of using tracks from mixed events -- iterate over all track pairs from the two collisions
//      int trackCount = 0;
//      for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
//        LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d)", t1.index(), t2.index(), c1.index(), c2.index(), t1.collision().index(), t2.collision().index());
//        trackCount++;
//        if (trackCount == 10)
//          break;
//      }
//    }
//  }
//};

struct MixedEventsCounters {
  // This task shows how to extract variables needed for weighted correlations.
  // NOTE: The same number of currentWindowNeighbours is returned for all kinds of block combinations.
  // Strictly upper: the first collision will is paired with exactly currentWindowNeighbours other collisions.
  // Upper: the first collision is paired with (currentWindowNeighbours + 1) collisions, including itself.
  // Full: (currentWindowNeighbours + 1) pairs with the first collision in the first position (c1)
  //       + there are other combinations with the first collision at other positions.

  std::vector<double> xBins{VARIABLE_WIDTH, -0.064, -0.062, -0.060, 0.066, 0.068, 0.070, 0.072};
  std::vector<double> yBins{VARIABLE_WIDTH, -0.320, -0.301, -0.300, 0.330, 0.340, 0.350, 0.360};
  using BinningType = ColumnBinningPolicy<aod::collision::PosX, aod::collision::PosY>;
  BinningType binningOnPositions{{xBins, yBins}, true};
  // true is for 'ignore overflows' (true by default)
  SameKindPair<aod::Collisions, aod::Tracks, BinningType> pair{binningOnPositions, 5, -1}; // indicates that 5 events should be mixed and under / overflow(-1) to be ignored

  void process(aod::Collisions const& collisions, aod::Tracks const& tracks)
  {
    LOGF(info, "Input data Collisions %d, Tracks %d ", collisions.size(), tracks.size());

    for (auto it = pair.begin(); it != pair.end(); it++) {
      auto& [c1, tracks1, c2, tracks2] = *it;
      LOGF(info, "Mixed event collisions: (%d, %d), is it first pair with %d: %d, number of collisions mixed with the first: %d", c1.globalIndex(), c2.globalIndex(), c1.globalIndex(), it.isNewWindow(), it.currentWindowNeighbours());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MixedEventsEmptyTables>(cfgc),
    adaptAnalysisTask<MixedEventsJoinedTracks>(cfgc),
    // adaptAnalysisTask<MixedEventsBadSubscription>(cfgc), // Should not compile
    adaptAnalysisTask<MixedEventsCounters>(cfgc),
  };
}
