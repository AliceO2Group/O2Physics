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
/// \brief Example task for event mixing.
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"

namespace o2::aod
{
namespace hash
{
DECLARE_SOA_COLUMN(Bin, bin, int);
} // namespace hash
DECLARE_SOA_TABLE(Hashes, "AOD", "HASH", hash::Bin);

using Hash = Hashes::iterator;
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct HashTask {
  // NOTE: The bins need to be ordered for the hash calculations to make sense
  std::vector<float> xBins{-0.064f, -0.062f, -0.060f, 0.066f, 0.068f, 0.070f, 0.072f};
  std::vector<float> yBins{-0.320f, -0.301f, -0.300f, 0.330f, 0.340f, 0.350f, 0.360f};
  Produces<aod::Hashes> hashes;

  // Calculate hash for an element based on 2 properties and their bins.
  int getHash(std::vector<float> const& xBins, std::vector<float> const& yBins, float colX, float colY)
  {
    if (colX < xBins[0] || colY < yBins[0]) {
      return -1;
    }
    for (unsigned int i = 1; i < xBins.size(); i++) {
      if (colX < xBins[i]) {
        for (unsigned int j = 1; j < yBins.size(); j++) {
          if (colY < yBins[j]) {
            return i + j * (xBins.size() + 1);
          }
        }
        return -1;
      }
    }
    return -1;
  }

  void process(aod::Collisions const& collisions)
  {
    for (auto& collision : collisions) {
      int hash = getHash(xBins, yBins, collision.posX(), collision.posY());
      LOGF(info, "Collision: %d (%f, %f, %f) hash: %d", collision.index(), collision.posX(), collision.posY(), collision.posZ(), hash);
      hashes(hash);
    }
  }
};

struct MixedEventsPartitionedTracks {
  Filter trackFilter = (aod::track::x > -0.8f) && (aod::track::x < 0.8f) && (aod::track::y > 1.0f);
  using myTracks = soa::Filtered<aod::Tracks>;

  Configurable<float> philow{"phiLow", 1.0f, "lowest phi"};
  Configurable<float> phiup{"phiUp", 2.0f, "highest phi"};

  void process(aod::Hashes const& hashes, aod::Collisions& collisions, myTracks& tracks)
  {
    collisions.bindExternalIndices(&tracks);
    auto tracksTuple = std::make_tuple(tracks);
    GroupSlicer slicer(collisions, tracksTuple);

    // Strictly upper categorised collisions
    for (auto& [c1, c2] : selfCombinations("fBin", 5, -1, join(hashes, collisions), join(hashes, collisions))) {
      // LOGF(info, "Collisions bin: %d pair: %d (%f, %f, %f), %d (%f, %f, %f)", c1.bin(), c1.index(), c1.posX(), c1.posY(), c1.posZ(), c2.index(), c2.posX(), c2.posY(), c2.posZ());

      auto it1 = slicer.begin();
      auto it2 = slicer.begin();
      for (auto& slice : slicer) {
        if (slice.groupingElement().index() == c1.index()) {
          it1 = slice;
          break;
        }
      }
      for (auto& slice : slicer) {
        if (slice.groupingElement().index() == c2.index()) {
          it2 = slice;
          break;
        }
      }
      auto tracks1 = std::get<myTracks>(it1.associatedTables());
      tracks1.bindExternalIndices(&collisions);
      auto tracks2 = std::get<myTracks>(it2.associatedTables());
      tracks2.bindExternalIndices(&collisions);

      Partition<myTracks> leftPhi1 = aod::track::phi < philow;
      leftPhi1.bindTable(tracks1);
      Partition<myTracks> leftPhi2 = aod::track::phi < philow;
      leftPhi2.bindTable(tracks2);
      Partition<myTracks> rightPhi1 = aod::track::phi >= phiup;
      rightPhi1.bindTable(tracks1);
      Partition<myTracks> rightPhi2 = aod::track::phi >= phiup;
      rightPhi2.bindTable(tracks2);

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

struct MixedEvents {
  SameKindPair<aod::Hashes, aod::Collisions, aod::Tracks> pair{"fBin", 5, -1};

  // Collisions must be first, not hashes!
  void process(aod::Collisions const& collisions, aod::Hashes const& hashes, aod::Tracks const& tracks)
  {
    LOGF(info, "Input data Collisions %d, Tracks %d ", collisions.size(), tracks.size());

    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      LOGF(info, "Mixed event bin: %d collisions: (%d, %d), tracks: (%d, %d)", c1.bin(), c1.globalIndex(), c2.globalIndex(), tracks1.size(), tracks2.size());
      for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
        LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d)", t1.index(), t2.index(), c1.index(), c2.index(), t1.collision().index(), t2.collision().index());
      }
    }
  }
};

struct MixedEventsInsideProcess {
  void process(aod::Collisions& collisions, aod::Hashes& hashes, aod::Tracks& tracks)
  {
    LOGF(info, "Input data Collisions %d, Tracks %d ", collisions.size(), tracks.size());

    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<aod::Hashes, aod::Collisions, aod::Tracks> pair{"fBin", 5, -1, hashes, collisions, tracksTuple};

    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      LOGF(info, "Mixed event bin: %d collisions: (%d, %d), tracks: (%d, %d)", c1.bin(), c1.globalIndex(), c2.globalIndex(), tracks1.size(), tracks2.size());
      for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
        LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d)", t1.globalIndex(), t2.globalIndex(), c1.globalIndex(), c2.globalIndex());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HashTask>(cfgc),
    adaptAnalysisTask<MixedEvents>(cfgc),
    adaptAnalysisTask<MixedEventsInsideProcess>(cfgc),
    adaptAnalysisTask<MixedEventsPartitionedTracks>(cfgc),
  };
}
