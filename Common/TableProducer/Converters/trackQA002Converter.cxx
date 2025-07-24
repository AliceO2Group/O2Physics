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
#include <limits>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

using namespace o2;
using namespace o2::framework;

struct TrackQAConverter002 {
  Produces<aod::TracksQA_002> tracksQA_002;

  void process000(aod::TracksQA_000 const& tracksQA_000)
  {
    for (const auto& trackQA : tracksQA_000) {
      tracksQA_002(
        trackQA.trackId(),
        trackQA.tpcTime0(),
        trackQA.tpcdcaR(),
        trackQA.tpcdcaZ(),
        trackQA.tpcClusterByteMask(),
        trackQA.tpcdEdxMax0R(),
        trackQA.tpcdEdxMax1R(),
        trackQA.tpcdEdxMax2R(),
        trackQA.tpcdEdxMax3R(),
        trackQA.tpcdEdxTot0R(),
        trackQA.tpcdEdxTot1R(),
        trackQA.tpcdEdxTot2R(),
        trackQA.tpcdEdxTot3R(),
        // dummy values, not available in _000
        std::numeric_limits<int8_t>::min(),  // deltaRefContParamY
        std::numeric_limits<int8_t>::min(),  // deltaRefContParamZ
        std::numeric_limits<int8_t>::min(),  // deltaRefContParamSnp
        std::numeric_limits<int8_t>::min(),  // deltaRefContParamTgl
        std::numeric_limits<int8_t>::min(),  // deltaRefContParamQ2Pt
        std::numeric_limits<int8_t>::min(),  // deltaRefGloParamY
        std::numeric_limits<int8_t>::min(),  // deltaRefGloParamZ
        std::numeric_limits<int8_t>::min(),  // deltaRefGloParamSnp
        std::numeric_limits<int8_t>::min(),  // deltaRefGloParamTgl
        std::numeric_limits<int8_t>::min(),  // deltaRefGloParamQ2Pt
        std::numeric_limits<int8_t>::min(),  // dTofdX
        std::numeric_limits<int8_t>::min()); // dTofdY
    }
  }
  PROCESS_SWITCH(TrackQAConverter002, process000, "process v000-to-v002 conversion", false);

  void process001(aod::TracksQA_001 const& tracksQA_001)
  {
    for (const auto& trackQA : tracksQA_001) {
      tracksQA_002(
        trackQA.trackId(),
        trackQA.tpcTime0(),
        trackQA.tpcdcaR(),
        trackQA.tpcdcaZ(),
        trackQA.tpcClusterByteMask(),
        trackQA.tpcdEdxMax0R(),
        trackQA.tpcdEdxMax1R(),
        trackQA.tpcdEdxMax2R(),
        trackQA.tpcdEdxMax3R(),
        trackQA.tpcdEdxTot0R(),
        trackQA.tpcdEdxTot1R(),
        trackQA.tpcdEdxTot2R(),
        trackQA.tpcdEdxTot3R(),
        trackQA.deltaRefContParamY(),
        trackQA.deltaRefITSParamZ(),
        trackQA.deltaRefContParamSnp(),
        trackQA.deltaRefContParamTgl(),
        trackQA.deltaRefContParamQ2Pt(),
        trackQA.deltaRefGloParamY(),
        trackQA.deltaRefGloParamZ(),
        trackQA.deltaRefGloParamSnp(),
        trackQA.deltaRefGloParamTgl(),
        trackQA.deltaRefGloParamQ2Pt(),
        // dummy values, not available in _001
        std::numeric_limits<int8_t>::min(),  // dTofdX
        std::numeric_limits<int8_t>::min()); // dTofdY
    }
  }
  PROCESS_SWITCH(TrackQAConverter002, process001, "process v001-to-v002 conversion", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TrackQAConverter002>(cfgc),
  };
}
