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

/// \file mftTracksConverter.cxx
/// \brief Converts MFTTracks table from version 000 to 001

/// This task allows for the conversion of the MFTTracks table from version 000 to 001.
/// The conversion is needed because the table has been extended with the MFTClusterSizesAndTracksFlags column
/// while nClusters and isCA columns are evaluated dynamically from it.
/// In the converter a dummy MFTClusterSizesAndTracksFlags is created and filled according to the number of clusters

/// \author L.Micheletti <luca.micheletti@cern.ch>

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

#include <cstdint>

using namespace o2;
using namespace o2::framework;

struct MftTracksConverter {
  Produces<aod::StoredMFTTracks_001> mftTracks_001;
  void process(aod::MFTTracks_000 const& mftTracks_000)
  {

    for (const auto& track0 : mftTracks_000) {
      uint64_t mftClusterSizesAndTrackFlags = 0;
      int8_t nClusters = track0.nClusters();

      for (int layer = 0; layer < 10; ++layer) {
        mftClusterSizesAndTrackFlags &= ~(0x3fULL << (layer * 6));
        mftClusterSizesAndTrackFlags |= (layer < nClusters) ? (1ULL << (layer * 6)) : 0;
      }

      mftTracks_001(track0.collisionId(),
                    track0.x(),
                    track0.y(),
                    track0.z(),
                    track0.phi(),
                    track0.tgl(),
                    track0.signed1Pt(),
                    mftClusterSizesAndTrackFlags,
                    track0.chi2(),
                    track0.trackTime(),
                    track0.trackTimeRes());
    }
  }
};

/// Spawn the extended table for MFTTracks001 to avoid the call to the internal spawner and a consequent circular dependency
struct MFTTracksSpawner {
  Spawns<aod::MFTTracks_001> mftTracks_001;
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MftTracksConverter>(cfgc),
    adaptAnalysisTask<MFTTracksSpawner>(cfgc),
  };
}
