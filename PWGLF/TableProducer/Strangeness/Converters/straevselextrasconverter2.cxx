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

// Produce dummy StraEvSelExtras for analysis subscribing to StraEvSelExtras but not saved in the strangeness derived data (typically when running over pp strangeness derived data)
struct straevselextrasconverter2 {
  Produces<aod::StraEvSelExtras> straEvSelExtras;

  void process(aod::StraEvSels const& straEvSels)
  {
    for (int ii = 0; ii < straEvSels.size(); ii++) {
      straEvSelExtras(-999.,  // dummy multZNA,
                      -999.,  // dummy multZNC,
                      -999.,  // dummy multZEM1,
                      -999.,  // dummy multZEM2,
                      -999.,  // dummy multZPA,
                      -999.,  // dummy multZPC,
                      -999.,  // dummy multNTracksITSTPC,
                      -999.,  // dummy multAllTracksTPCOnly,
                      -999.,  // dummy multAllTracksITSTPC,
                      -999.,  // dummy trackOccupancyInTimeRange,
                      -999.,  // dummy ft0cOccupancyInTimeRange,
                      -999.,  // dummy timeFDDA,
                      -999.,  // dummy timeFDDC,
                      -999.,  // dummy timeFV0A,
                      -999.,  // dummy timeFT0A,
                      -999.,  // dummy timeFT0C,
                      0,      // dummy triggerMaskFT0,
                      -999,   // dummy gapSide,
                      -999.,  // dummy totalFT0AmplitudeA,
                      -999.,  // dummy totalFT0AmplitudeC,
                      -999.,  // dummy totalFV0AmplitudeA,
                      -999.,  // dummy totalFDDAmplitudeA,
                      -999.,  // dummy totalFDDAmplitudeC,
                      -999.,  // dummy timeZNA,
                      -999.,  // dummy timeZNC,
                      -999.,  // dummy energyCommonZNA,
                      -999.); // dummy energyCommonZNC);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<straevselextrasconverter2>(cfgc)};
}
