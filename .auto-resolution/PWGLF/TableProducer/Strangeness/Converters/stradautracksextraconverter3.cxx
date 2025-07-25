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

// Converts daughter TracksExtra from 2 to 3
struct stradautracksextraconverter3 {
  Produces<aod::DauTrackExtras_003> dauTrackExtras_003;

  void process(aod::DauTrackExtras_002 const& dauTrackExtras_002)
  {
    for (auto& values : dauTrackExtras_002) {
      dauTrackExtras_003(values.itsChi2PerNcl(),
                         -1 /* dummy tpcChi2PerNcl value */,
                         values.detectorMap(),
                         values.itsClusterSizes(),
                         values.tpcNClsFindable(),
                         values.tpcNClsFindableMinusFound(),
                         values.tpcNClsFindableMinusCrossedRows(),
                         -1 /* dummy tpcNClsShared value */);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<stradautracksextraconverter3>(cfgc)};
}
