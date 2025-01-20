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
/// \brief Example tasks for event mixing.
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

namespace o2::aod
{
namespace hash
{
DECLARE_SOA_COLUMN(Bin, bin, int);
} // namespace hash
DECLARE_SOA_TABLE(MixingHashes, "AOD", "HASH", hash::Bin);

} // namespace o2::aod

struct MixedEvents {
  SliceCache cache;
  std::vector<double> xBins{VARIABLE_WIDTH, -0.064, -0.062, -0.060, 0.066, 0.068, 0.070, 0.072};
  std::vector<double> yBins{VARIABLE_WIDTH, -0.320, -0.301, -0.300, 0.330, 0.340, 0.350, 0.360};
  using BinningType = ColumnBinningPolicy<aod::collision::PosX, aod::collision::PosY>;
  BinningType binningOnPositions{{xBins, yBins}, true};                                    // true is for 'ignore overflows' (true by default)
  SameKindPair<aod::Collisions, aod::Tracks, BinningType> pair{binningOnPositions, 5, -1, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

  void process(aod::Collisions const& collisions, aod::Tracks const& tracks)
  {
    LOGF(info, "Input data Collisions %d, Tracks %d ", collisions.size(), tracks.size());

    int count = 0;
    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      // Example of using getBin() to check collisions bin
      // NOTE that it is a bit different for FlexibleBinning -- check in respective example
      int bin = binningOnPositions.getBin({c1.posX(), c1.posY()});

      LOGF(info, "Mixed event collisions: (%d, %d) from bin %d", c1.globalIndex(), c2.globalIndex(), bin);
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

struct MixedEventsInsideProcess {
  SliceCache cache;
  std::vector<double> xBins{VARIABLE_WIDTH, -0.064, -0.062, -0.060, 0.066, 0.068, 0.070, 0.072};
  std::vector<double> yBins{VARIABLE_WIDTH, -0.320, -0.301, -0.300, 0.330, 0.340, 0.350, 0.360};
  using BinningType = ColumnBinningPolicy<aod::collision::PosX, aod::collision::PosY>;
  BinningType binningOnPositions{{xBins, yBins}, true}; // true is for 'ignore overflows' (true by default)

  void process(aod::Collisions& collisions, aod::Tracks& tracks)
  {
    LOGF(info, "Input data Collisions %d, Tracks %d ", collisions.size(), tracks.size());

    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<aod::Collisions, aod::Tracks, BinningType> pair{binningOnPositions, 5, -1, collisions, tracksTuple, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

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

struct MixedEventsFilteredTracks {
  SliceCache cache;
  Filter trackFilter = aod::track::eta < 1.0f;
  using myTracks = soa::Filtered<aod::Tracks>;

  std::vector<double> xBins{VARIABLE_WIDTH, -0.064, -0.062, -0.060, 0.066, 0.068, 0.070, 0.072};
  std::vector<double> yBins{VARIABLE_WIDTH, -0.320, -0.301, -0.300, 0.330, 0.340, 0.350, 0.360};
  using BinningType = ColumnBinningPolicy<aod::collision::PosX, aod::collision::PosY>;
  BinningType binningOnPositions{{xBins, yBins}, true};                                 // true is for 'ignore overflows' (true by default)
  SameKindPair<aod::Collisions, myTracks, BinningType> pair{binningOnPositions, 5, -1, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

  void process(aod::Collisions const& collisions, myTracks const& tracks)
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
        LOGF(info, "Mixed event tracks pair: (%d, %d) (%.2f, %.2f) < 1.0f from events (%d, %d), track event: (%d, %d)", t1.index(), t2.index(), t1.eta(), t2.eta(), c1.index(), c2.index(), t1.collision().index(), t2.collision().index());
        trackCount++;
        if (trackCount == 10)
          break;
      }
    }
  }
};

struct MixedEventsJoinedCollisions {
  SliceCache cache;
  using aodCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>;
  std::vector<double> xBins{VARIABLE_WIDTH, -0.064, -0.062, -0.060, 0.066, 0.068, 0.070, 0.072};
  std::vector<double> yBins{VARIABLE_WIDTH, -0.320, -0.301, -0.300, 0.330, 0.340, 0.350, 0.360};
  using BinningType = ColumnBinningPolicy<aod::collision::PosX, aod::collision::PosY>;
  BinningType binningOnPositions{{xBins, yBins}, true};                                  // true is for 'ignore overflows' (true by default)
  SameKindPair<aodCollisions, aod::Tracks, BinningType> pair{binningOnPositions, 5, -1, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

  void process(aodCollisions& collisions, aod::Tracks const& tracks, aod::BCsWithTimestamps const&)
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

struct MixedEventsDynamicColumns {
  SliceCache cache;
  using aodCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
  std::vector<double> zBins{7, -7, 7};
  std::vector<double> multBins{VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1};
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFV0M<aod::mult::MultFV0A, aod::mult::MultFV0C>>;
  BinningType corrBinning{{zBins, multBins}, true};                               // true is for 'ignore overflows' (true by default)
  SameKindPair<aodCollisions, aod::Tracks, BinningType> pair{corrBinning, 5, -1, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

  void process(aodCollisions& collisions, aod::Tracks const& tracks)
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

struct MixedEventsVariousKinds {
  SliceCache cache;
  std::vector<double> xBins{VARIABLE_WIDTH, -0.064, -0.062, -0.060, 0.066, 0.068, 0.070, 0.072};
  std::vector<double> yBins{VARIABLE_WIDTH, -0.320, -0.301, -0.300, 0.330, 0.340, 0.350, 0.360};
  using BinningType = ColumnBinningPolicy<aod::collision::PosX, aod::collision::PosY>;
  BinningType binningOnPositions{{xBins, yBins}, true};                                      // true is for 'ignore overflows' (true by default)
  Pair<aod::Collisions, aod::Tracks, aod::V0s, BinningType> pair{binningOnPositions, 5, -1, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

  void process(aod::Collisions const& collisions, aod::Tracks const& tracks, aod::V0s const& v0s)
  {
    LOGF(info, "Input data Collisions %d, Tracks %d V0s %d", collisions.size(), tracks.size(), v0s.size());

    int count = 0;
    // tracks1 is an aod::Tracks table of tracks belonging to collision c1 (aod::Collision::iterator)
    // tracks2 is an aod::V0s table of V0s belonging to collision c2 (aod::Collision::iterator)
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

struct MixedEventsTriple {
  SliceCache cache;
  std::vector<double> xBins{VARIABLE_WIDTH, -0.064, -0.062, -0.060, 0.066, 0.068, 0.070, 0.072};
  std::vector<double> yBins{VARIABLE_WIDTH, -0.320, -0.301, -0.300, 0.330, 0.340, 0.350, 0.360};
  std::vector<double> zBins{7, -7, 7};
  using BinningType = ColumnBinningPolicy<aod::collision::PosX, aod::collision::PosY, aod::collision::PosZ>;
  BinningType binningOnPositions{{xBins, yBins, zBins}, true};                                 // true is for 'ignore overflows' (true by default)
  SameKindTriple<aod::Collisions, aod::Tracks, BinningType> triple{binningOnPositions, 5, -1, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

  void process(aod::Collisions const& collisions, aod::Tracks const& tracks)
  {
    LOGF(info, "Input data Collisions %d, Tracks %d", collisions.size(), tracks.size());

    int count = 0;
    for (auto& [c1, tracks1, c2, tracks2, c3, tracks3] : triple) {
      LOGF(info, "Mixed event collisions: (%d, %d, %d)", c1.globalIndex(), c2.globalIndex(), c3.globalIndex());
      count++;
      if (count == 100)
        break;

      // Example of using tracks from mixed events -- iterate over all track triplets from the three collisions
      int trackCount = 0;
      for (auto& [t1, t2, t3] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2, tracks3))) {
        LOGF(info, "Mixed event tracks triple: (%d, %d, %d) from events (%d, %d, %d), track event: (%d, %d, %d)", t1.index(), t2.index(), t3.index(), c1.index(), c2.index(), c3.index(), t1.collision().index(), t2.collision().index(), t3.collision().index());
        trackCount++;
        if (trackCount == 10)
          break;
      }
    }
  }
};

struct MixedEventsTripleVariousKinds {
  SliceCache cache;
  std::vector<double> xBins{VARIABLE_WIDTH, -0.064, -0.062, -0.060, 0.066, 0.068, 0.070, 0.072};
  std::vector<double> yBins{VARIABLE_WIDTH, -0.320, -0.301, -0.300, 0.330, 0.340, 0.350, 0.360};
  std::vector<double> zBins{7, -7, 7};
  using BinningType = ColumnBinningPolicy<aod::collision::PosX, aod::collision::PosY, aod::collision::PosZ>;
  BinningType binningOnPositions{{xBins, yBins, zBins}, true};                                                // true is for 'ignore overflows' (true by default)
  Triple<aod::Collisions, aod::Tracks, aod::V0s, aod::Tracks, BinningType> triple{binningOnPositions, 5, -1, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

  void process(aod::Collisions const& collisions, aod::Tracks const& tracks, aod::V0s const& v0s)
  {
    LOGF(info, "Input data Collisions %d, Tracks %d V0s %d", collisions.size(), tracks.size(), v0s.size());

    int count = 0;
    // tracks1 is an aod::Tracks table of tracks belonging to collision c1 (aod::Collision::iterator)
    // tracks2 is an aod::V0s table of V0s belonging to collision c2 (aod::Collision::iterator)
    // tracks3 is an aod::Tracks table of tracks belonging to collision c3 (aod::Collision::iterator)
    for (auto& [c1, tracks1, c2, tracks2, c3, tracks3] : triple) {
      LOGF(info, "Mixed event collisions: (%d, %d, %d)", c1.globalIndex(), c2.globalIndex(), c3.globalIndex());
      count++;
      if (count == 100)
        break;

      // Example of using tracks from mixed events -- iterate over all track triplets from the three collisions
      int trackCount = 0;
      for (auto& [t1, t2, t3] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2, tracks3))) {
        LOGF(info, "Mixed event tracks triple: (%d, %d, %d) from events (%d, %d, %d), track event: (%d, %d, %d)", t1.index(), t2.index(), t3.index(), c1.index(), c2.index(), c3.index(), t1.collision().index(), t2.collision().index(), t3.collision().index());
        trackCount++;
        if (trackCount == 10)
          break;
      }
    }
  }
};

struct HashTask {
  std::vector<float> xBins{-0.064f, -0.062f, -0.060f, 0.066f, 0.068f, 0.070f, 0.072};
  std::vector<float> yBins{-0.320f, -0.301f, -0.300f, 0.330f, 0.340f, 0.350f, 0.360};
  Produces<aod::MixingHashes> hashes;

  // Calculate hash for an element based on 2 properties and their bins.
  int getHash(const std::vector<float>& xBins, const std::vector<float>& yBins, float colX, float colY)
  {
    // underflow
    if (colX < xBins[0] || colY < yBins[0])
      return -1;

    for (unsigned int i = 1; i < xBins.size(); i++) {
      if (colX < xBins[i]) {
        for (unsigned int j = 1; j < yBins.size(); j++) {
          if (colY < yBins[j]) {
            return (i - 1) + (j - 1) * (xBins.size() - 1);
          }
        }
        // overflow for yBins only
        return -1;
      }
    }

    // overflow
    return -1;
  }

  void process(aod::Collisions const& collisions)
  {
    for (auto& collision : collisions) {
      int hash = getHash(xBins, yBins, collision.posX(), collision.posY());
      // LOGF(info, "Collision: %d (%f, %f, %f) hash: %d", collision.index(), collision.posX(), collision.posY(), collision.posZ(), hash);
      hashes(hash);
    }
  }
};

struct MixedEventsWithHashTask {
  SliceCache cache;
  using myCollisions = soa::Join<aod::MixingHashes, aod::Collisions>;
  NoBinningPolicy<aod::hash::Bin> hashBin;
  SameKindPair<myCollisions, aod::Tracks, NoBinningPolicy<aod::hash::Bin>> pair{hashBin, 5, -1, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

  void process(myCollisions& collisions, aod::Tracks& tracks)
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
        LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d)", t1.globalIndex(), t2.globalIndex(), c1.globalIndex(), c2.globalIndex());
        trackCount++;
        if (trackCount == 10)
          break;
      }
    }
  }
};

struct MixedEventsPartitionedTracks {
  SliceCache cache;
  Filter trackFilter = (aod::track::x > -0.8f) && (aod::track::x < 0.8f) && (aod::track::y > 1.0f);
  using myTracks = soa::Filtered<aod::Tracks>;

  Configurable<float> philow{"phiLow", 1.0f, "lowest phi"};
  Configurable<float> phiup{"phiUp", 2.0f, "highest phi"};

  std::vector<double> xBins{VARIABLE_WIDTH, -0.064, -0.062, -0.060, 0.066, 0.068, 0.070, 0.072};
  std::vector<double> yBins{VARIABLE_WIDTH, -0.320, -0.301, -0.300, 0.330, 0.340, 0.350, 0.360};
  using BinningType = ColumnBinningPolicy<aod::collision::PosX, aod::collision::PosY>;
  BinningType binningOnPositions{{xBins, yBins}, true};                                 // true is for 'ignore overflows' (true by default)
  SameKindPair<aod::Collisions, myTracks, BinningType> pair{binningOnPositions, 5, -1, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

  void process(aod::Collisions const& collisions, myTracks const& tracks)
  {
    LOGF(info, "Input data Collisions %d, Tracks %d ", collisions.size(), tracks.size());

    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      Partition<myTracks> leftPhi1 = aod::track::phi < philow;
      Partition<myTracks> leftPhi2 = aod::track::phi < philow;
      Partition<myTracks> rightPhi1 = aod::track::phi >= phiup;
      Partition<myTracks> rightPhi2 = aod::track::phi >= phiup;
      leftPhi1.bindTable(tracks1);
      leftPhi2.bindTable(tracks2);
      rightPhi1.bindTable(tracks1);
      rightPhi2.bindTable(tracks2);

      LOGF(info, "Mixed event collisions: (%d, %d), tracks: (%d, %d), left phis: (%d, %d), right phis: (%d, %d)", c1.globalIndex(), c2.globalIndex(), tracks1.size(), tracks2.size(), leftPhi1.size(), leftPhi2.size(), rightPhi1.size(), rightPhi2.size());

      // Example of using tracks from mixed events -- iterate over all track pairs from the two partitions from two collisions
      for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(leftPhi1, leftPhi2))) {
        if (t1.phi() >= (float)philow || t2.phi() >= (float)philow) {
          LOGF(info, "WRONG Mixed event left tracks pair: (%d, %d) from events (%d, %d), phi: (%.3f. %.3f) < %.3f", t1.index(), t2.index(), c1.index(), c2.index(), t1.phi(), t2.phi(), (float)philow);
        }
      }
      for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(rightPhi1, rightPhi2))) {
        if (t1.phi() < (float)phiup || t2.phi() < (float)phiup) {
          LOGF(info, "WRONG Mixed event right tracks pair: (%d, %d) from events (%d, %d), phi: (%.3f. %.3f) >= %.3f", t1.index(), t2.index(), c1.index(), c2.index(), t1.phi(), t2.phi(), (float)phiup);
        }
      }
    }
  }
};

struct MixedEventsLambdaBinning {
  SliceCache cache;
  Preslice<aod::Tracks> perCol = aod::track::collisionId;

  ConfigurableAxis axisVertex{"axisVertex", {14, -7, 7}, "vertex axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0.0, 2.750, 5.250, 7.750, 12.750, 17.750, 22.750, 27.750, 32.750, 37.750, 42.750, 47.750, 52.750, 57.750, 62.750, 67.750, 72.750, 77.750, 82.750, 87.750, 92.750, 97.750, 250.1}, "multiplicity axis for histograms"};

  void process(aod::Collisions& collisions, aod::Tracks& tracks)
  {
    LOGF(info, "Input data Collisions %d, Tracks %d ", collisions.size(), tracks.size());

    auto getTracksSize =
      [&tracks, this](aod::Collision const& col) {
        auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, col.globalIndex(), this->cache); // it's cached, so slicing/grouping happens only once
        return associatedTracks.size();
      };

    using BinningType = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;
    BinningType binningWithLambda{{getTracksSize}, {axisVertex, axisMultiplicity}, true};

    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<aod::Collisions, aod::Tracks, BinningType> pair{binningWithLambda, 5, -1, collisions, tracksTuple, &cache}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

    int count = 0;
    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      // NOTE that getBin() with FlexibleBinningPolicy needs explicit tuple construction
      int bin = binningWithLambda.getBin(std::tuple(getTracksSize(c1), c1.posZ()));

      LOGF(info, "Mixed event collisions: (%d, %d) from bin %d z: (%.3f, %.3f), tracks size (%d, %d)", c1.globalIndex(), c2.globalIndex(), bin, c1.posZ(), c2.posZ(), tracks1.size(), tracks2.size());
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

struct MixedEventsCounters {
  // Please read carefully the tutorial description in O2 documentation to understand the behaviour of currentWindowNeighbours().
  // https://aliceo2group.github.io/analysis-framework/docs/tutorials/eventMixing.html#mixedeventscounters

  SliceCache cache;

  std::vector<double> xBins{VARIABLE_WIDTH, -0.064, -0.062, -0.060, 0.066, 0.068, 0.070, 0.072};
  std::vector<double> yBins{VARIABLE_WIDTH, -0.320, -0.301, -0.300, 0.330, 0.340, 0.350, 0.360};
  using BinningType = ColumnBinningPolicy<aod::collision::PosX, aod::collision::PosY>;
  BinningType binningOnPositions{{xBins, yBins}, true};
  // true is for 'ignore overflows' (true by default)
  SameKindPair<aod::Collisions, aod::Tracks, BinningType> pair{binningOnPositions, 5, -1, &cache}; // indicates that 5 events should be mixed and under / overflow(-1) to be ignored

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
    adaptAnalysisTask<MixedEvents>(cfgc),
    adaptAnalysisTask<MixedEventsInsideProcess>(cfgc),
    adaptAnalysisTask<MixedEventsFilteredTracks>(cfgc),
    adaptAnalysisTask<MixedEventsJoinedCollisions>(cfgc),
    adaptAnalysisTask<MixedEventsDynamicColumns>(cfgc),
    adaptAnalysisTask<MixedEventsVariousKinds>(cfgc),
    adaptAnalysisTask<MixedEventsTriple>(cfgc),
    adaptAnalysisTask<MixedEventsTripleVariousKinds>(cfgc),
    adaptAnalysisTask<HashTask>(cfgc),
    adaptAnalysisTask<MixedEventsWithHashTask>(cfgc),
    adaptAnalysisTask<MixedEventsPartitionedTracks>(cfgc),
    adaptAnalysisTask<MixedEventsLambdaBinning>(cfgc),
    adaptAnalysisTask<MixedEventsCounters>(cfgc),
  };
}
