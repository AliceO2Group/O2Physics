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

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/AggregatedRunInfo.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::evsel;

// Converts straevselsextrasconverter1 converts StraEvSelExtras_000 into StraEvSelExtras_001
struct straevselextrasconverter {
  Produces<aod::StraEvSelExtras_001> straEvSelExtras_001;

  void process(soa::Join<aod::StraEvSels_005, aod::StraEvSelExtras_000> const& straEvSels_005)
  {
    for (auto& values : straEvSels_005) {
      straEvSelExtras_001(values.multZNA(),
                          values.multZNC(),
                          values.multZEM1(),
                          values.multZEM2(),
                          values.multZPA(),
                          values.multZPC(),
                          values.multNTracksITSTPC(),
                          values.multAllTracksTPCOnly(),
                          values.multAllTracksITSTPC(),
                          values.trackOccupancyInTimeRange(),
                          values.ft0cOccupancyInTimeRange(),
                          values.timeFDDA(),
                          values.timeFDDC(),
                          values.timeFV0A(),
                          values.timeFT0A(),
                          values.timeFT0C(),
                          values.triggerMaskFT0(),
                          values.gapSide(),
                          values.totalFT0AmplitudeA(),
                          values.totalFT0AmplitudeC(),
                          values.totalFV0AmplitudeA(),
                          values.totalFDDAmplitudeA(),
                          values.totalFDDAmplitudeC(),
                          values.timeZNA(),
                          values.timeZNC(),
                          values.energyCommonZNA(),
                          values.energyCommonZNC());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<straevselextrasconverter>(cfgc)};
}
