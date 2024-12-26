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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

using namespace o2;
using namespace o2::framework;

struct TracksExtraV002Converter {
  Produces<aod::StoredTracksExtra_002> tracksExtra_002;

  void init(InitContext const&)
  {
    if (doprocessV000ToV002 == false && doprocessV001ToV002 == false) {
      LOGF(fatal, "Neither processV000ToV002 nor processV001ToV002 is enabled. Please choose one!");
    }
    if (doprocessV000ToV002 == true && doprocessV001ToV002 == true) {
      LOGF(fatal, "Both processV000ToV002 and processV001ToV002 are enabled. Please choose only one!");
    }
  }

  void processV000ToV002(aod::TracksExtra_000 const& tracksExtra_000)
  {

    for (const auto& track0 : tracksExtra_000) {

      uint32_t itsClusterSizes = 0;
      for (int layer = 0; layer < 7; layer++) {
        if (track0.itsClusterMap() & (1 << layer)) {
          itsClusterSizes |= (0xf << (layer * 4));
        }
      }

      int8_t TPCNClsFindableMinusPID = 0;

      tracksExtra_002(track0.tpcInnerParam(),
                      track0.flags(),
                      itsClusterSizes,
                      track0.tpcNClsFindable(),
                      track0.tpcNClsFindableMinusFound(),
                      TPCNClsFindableMinusPID,
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
  PROCESS_SWITCH(TracksExtraV002Converter, processV000ToV002, "process v000-to-v002 conversion", false);

  void processV001ToV002(aod::TracksExtra_001 const& tracksExtra_001)
  {

    for (const auto& track1 : tracksExtra_001) {

      int8_t TPCNClsFindableMinusPID = 0;

      tracksExtra_002(track1.tpcInnerParam(),
                      track1.flags(),
                      track1.itsClusterSizes(),
                      track1.tpcNClsFindable(),
                      track1.tpcNClsFindableMinusFound(),
                      TPCNClsFindableMinusPID,
                      track1.tpcNClsFindableMinusCrossedRows(),
                      track1.tpcNClsShared(),
                      track1.trdPattern(),
                      track1.itsChi2NCl(),
                      track1.tpcChi2NCl(),
                      track1.trdChi2(),
                      track1.tofChi2(),
                      track1.tpcSignal(),
                      track1.trdSignal(),
                      track1.length(),
                      track1.tofExpMom(),
                      track1.trackEtaEmcal(),
                      track1.trackPhiEmcal(),
                      track1.trackTime(),
                      track1.trackTimeRes());
    }
  }
  PROCESS_SWITCH(TracksExtraV002Converter, processV001ToV002, "process v001-to-v002 conversion", true);
};

/// Spawn the extended table for TracksExtra002 to avoid the call to the internal spawner and a consequent circular dependency
struct TracksExtraSpawner {
  Spawns<aod::TracksExtra_002> tracksExtra_002;
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TracksExtraV002Converter>(cfgc),
    adaptAnalysisTask<TracksExtraSpawner>(cfgc),
  };
}
