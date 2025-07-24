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
/// \brief Iterate over pair of tracks with combinations interface, without and with binning.
//         Iterate over pairs of collisions with configurable binning.
//         Iterate over pairs of tracks with binning predefined by a hash task.
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/Multiplicity.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;

namespace o2::aod
{
namespace hash
{
DECLARE_SOA_COLUMN(Bin, bin, int);
} // namespace hash
DECLARE_SOA_TABLE(MixingHashes, "AOD", "HASH", hash::Bin);

} // namespace o2::aod

struct TrackCombinations {
  void process(aod::Tracks const& tracks)
  {
    int count = 0;
    // Strictly upper tracks
    for (auto& [t0, t1] : combinations(tracks, tracks)) {
      LOGF(info, "Tracks pair: %d %d", t0.index(), t1.index());
      count++;
      if (count > 100)
        break;
    }
  }
};

struct BinnedTrackCombinations {
  std::vector<double> xBins{VARIABLE_WIDTH, -0.064, -0.062, -0.060, 0.066, 0.068, 0.070, 0.072};
  std::vector<double> yBins{VARIABLE_WIDTH, -0.320, -0.301, -0.300, 0.330, 0.340, 0.350, 0.360};
  ColumnBinningPolicy<aod::track::X, aod::track::Y> trackBinning{{xBins, yBins}, true};

  void process(aod::Tracks const& tracks)
  {
    int count = 0;
    // Strictly upper tracks binned by x and y position
    for (auto& [t0, t1] : selfCombinations(trackBinning, 5, -1, tracks, tracks)) {
      int bin = trackBinning.getBin({t0.x(), t0.y()});
      LOGF(info, "Tracks bin: %d pair: %d (%f, %f, %f), %d (%f, %f, %f)", bin, t0.index(), t0.x(), t0.y(), t0.z(), t1.index(), t1.x(), t1.y(), t1.z());
      count++;
      if (count > 100)
        break;
    }
  }
};

struct ConfigurableBinnedCollisionCombinations {
  ConfigurableAxis axisVertex{"axisVertex", {7, -7, 7}, "vertex axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity / centrality axis for histograms"};

  void process(soa::Join<aod::Collisions, aod::Mults> const& collisions)
  {
    using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFV0M<aod::mult::MultFV0A, aod::mult::MultFV0C>>;
    BinningType colBinning{{axisVertex, axisMultiplicity}, true};

    int count = 0;
    // Strictly upper tracks binned by x and y position
    for (auto& [c0, c1] : selfCombinations(colBinning, 5, -1, collisions, collisions)) {
      int bin = colBinning.getBin({c0.posZ(), c0.multFV0M()});
      LOGF(info, "Collision bin: %d pair: %d (%f, %f), %d (%f, %f)", bin, c0.index(), c0.posZ(), c0.multFV0M(), c1.index(), c1.posZ(), c1.multFV0A());
      count++;
      if (count > 100)
        break;
    }
  }
};

struct BinnedTrackPartitionsCombinations {
  std::vector<double> xBins{VARIABLE_WIDTH, -0.064, -0.062, -0.060, 0.066, 0.068, 0.070, 0.072};
  std::vector<double> yBins{VARIABLE_WIDTH, -0.320, -0.301, -0.300, 0.330, 0.340, 0.350, 0.360};
  ColumnBinningPolicy<aod::track::X, aod::track::Y> trackBinning{{xBins, yBins}, true};
  Configurable<float> philow{"phiLow", 1.0f, "lowest phi"};

  void process(aod::Tracks const& tracks)
  {
    Partition<aod::Tracks> leftPhi = aod::track::phi < philow;
    Partition<aod::Tracks> rightPhi = aod::track::phi >= philow;
    leftPhi.bindTable(tracks);
    rightPhi.bindTable(tracks);

    int count = 0;
    // Strictly upper tracks binned by x and y position
    for (auto& [t0, t1] : selfCombinations(trackBinning, 5, -1, leftPhi, rightPhi)) {
      int bin = trackBinning.getBin({t0.x(), t0.y()});
      LOGF(info, "Tracks bin: %d pair: %d (%f, %f, %f) phi %f, %d (%f, %f, %f) phi %f", bin, t0.index(), t0.x(), t0.y(), t0.z(), t0.phi(), t1.index(), t1.x(), t1.y(), t1.z(), t1.phi());
      count++;
      if (count > 100)
        break;
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
            return i + j * (xBins.size() + 1);
          }
        }
        // overflow for yBins only
        return -1;
      }
    }

    // overflow
    return -1;
  }

  void process(aod::Tracks const& tracks)
  {
    for (auto& track : tracks) {
      int hash = getHash(xBins, yBins, track.x(), track.y());
      LOGF(info, "Track: %d (%f, %f, %f) hash: %d", track.index(), track.x(), track.y(), track.z(), hash);
      hashes(hash);
    }
  }
};

struct BinnedTrackCombinationsWithHashTable {
  NoBinningPolicy<aod::hash::Bin> hashBin;

  void process(soa::Join<aod::MixingHashes, aod::Tracks> const& hashedTracks)
  {
    int count = 0;
    // Strictly upper categorised tracks
    for (auto& [t0, t1] : selfCombinations(hashBin, 5, -1, hashedTracks, hashedTracks)) {
      int bin = hashBin.getBin({t0.bin()});
      LOGF(info, "Tracks bin: %d from policy: %d  pair: %d (%f, %f, %f), %d (%f, %f, %f)", t0.bin(), bin, t0.index(), t0.x(), t0.y(), t0.z(), t1.index(), t1.x(), t1.y(), t1.z());
      count++;
      if (count > 100)
        break;
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TrackCombinations>(cfgc),
    adaptAnalysisTask<BinnedTrackCombinations>(cfgc),
    adaptAnalysisTask<ConfigurableBinnedCollisionCombinations>(cfgc),
    adaptAnalysisTask<BinnedTrackPartitionsCombinations>(cfgc),
    adaptAnalysisTask<HashTask>(cfgc),
    adaptAnalysisTask<BinnedTrackCombinationsWithHashTable>(cfgc),
  };
}
