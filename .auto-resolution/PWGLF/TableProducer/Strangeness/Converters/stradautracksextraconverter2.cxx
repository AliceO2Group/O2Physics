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
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"

using namespace o2;
using namespace o2::framework;

// Converts daughter TracksExtra from 1 to 2
struct stradautracksextraconverter2 {
  Produces<aod::DauTrackExtras_002> dauTrackExtras_002;

  void process(aod::DauTrackExtras_001 const& dauTrackExtras_001)
  {
    for (auto& values : dauTrackExtras_001) {
      const int maxFindable = 130; // synthetic findable to ensure range is ok
      int findableMinusFound = maxFindable - values.tpcClusters();
      int findableMinusCrossedRows = maxFindable - values.tpcCrossedRows();
      dauTrackExtras_002(values.itsChi2PerNcl(),
                         values.detectorMap(),
                         values.itsClusterSizes(),
                         static_cast<uint8_t>(maxFindable),              // findable (unknown in old format)
                         static_cast<int8_t>(findableMinusFound),        // findable minus found: we know found
                         static_cast<int8_t>(findableMinusCrossedRows)); // findable minus crossed rows: we know crossed rows
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<stradautracksextraconverter2>(cfgc)};
}
