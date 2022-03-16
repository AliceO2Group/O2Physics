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
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;

struct TrackCombinations {
  void process(aod::Tracks const& tracks)
  {
    // Strictly upper tracks
    for (auto& [t0, t1] : combinations(tracks, tracks)) {
      LOGF(info, "Tracks pair: %d %d", t0.index(), t1.index());
    }
  }
};

struct BinnedTrackCombinations {
  std::vector<double> xBins{-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5};
  std::vector<double> yBins{-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5};
  BinningPolicy<aod::track::X, aod::track::Y> trackBinning{{xBins, yBins}, true};

  void process(aod::Tracks const& tracks)
  {
    // Strictly upper tracks binned by x and y position
    for (auto& [t0, t1] : selfCombinations(trackBinning, 5, -1, tracks, tracks)) {
      LOGF(info, "Tracks bin: %d pair: %d (%f, %f, %f), %d (%f, %f, %f)", t0.bin(), t0.index(), t0.x(), t0.y(), t0.z(), t1.index(), t1.x(), t1.y(), t1.z());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TrackCombinations>(cfgc),
    adaptAnalysisTask<BinnedTrackCombinations>(cfgc),
  };
}
