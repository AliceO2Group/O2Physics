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

/// \file tracksExtraConverter.cxx
/// \brief Converts TracksExtra table from version 000 to 001

/// This task allows for the conversion of the TracksExtra table from version 000 to 001.
/// The conversion is needed because the table has been extended with the ITSClusterSize column
/// and the ITSClusterMap column is evaluated dynamically from it.
/// In the converter a dummy ITSClusterSize column is filled with overflows if a hit in the layer is present

/// \author F.Mazzaschi <fmazzasc@cern.ch>

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

#include <cstdint>

using namespace o2;
using namespace o2::framework;

struct TracksExtraConverter {
  Produces<aod::StoredTracksExtra_001> tracksExtra_001;
  void process(aod::TracksExtra_000 const& tracksExtra_000)
  {

    for (const auto& track0 : tracksExtra_000) {

      uint32_t itsClusterSizes = 0;
      for (int layer = 0; layer < 7; layer++) {
        if (track0.itsClusterMap() & (1 << layer)) {
          itsClusterSizes |= (0xf << (layer * 4));
        }
      }

      tracksExtra_001(track0.tpcInnerParam(),
                      track0.flags(),
                      itsClusterSizes,
                      track0.tpcNClsFindable(),
                      track0.tpcNClsFindableMinusFound(),
                      track0.tpcNClsFindableMinusCrossedRows(),
                      track0.tpcNClsShared(),
                      track0.trdPattern(),
                      track0.itsChi2NCl(),
                      track0.tpcChi2NCl(),
                      track0.trdChi2(),
                      track0.tofChi2(),
                      track0.tpcSignal(),
                      track0.trdSignal(),
                      track0.length(),
                      track0.tofExpMom(),
                      track0.trackEtaEmcal(),
                      track0.trackPhiEmcal(),
                      track0.trackTime(),
                      track0.trackTimeRes());
    }
  }
};

/// Spawn the extended table for TracksExtra001 to avoid the call to the internal spawner and a consequent circular dependency
struct TracksExtraSpawner {
  Spawns<aod::TracksExtra_001> tracksExtra_001;
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TracksExtraConverter>(cfgc),
    adaptAnalysisTask<TracksExtraSpawner>(cfgc),
  };
}
