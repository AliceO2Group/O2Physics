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
//
/// \file straevselextrasconverter.cxx
/// \brief Converts straevselsextrasconverter1 converts StraEvSelExtras_000 into StraEvSelExtras_001
///
/// \author David Dobrigkeit Chinellato <david.dobrigkeit.chinellato@cern.ch>, Austrian Academy of Sciences & MBI
/// \author Romain Schotter <romain.schotter@cern.ch>, Austrian Academy of Sciences & MBI
//
//__________________________________________________
//

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/AggregatedRunInfo.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::evsel;

struct straevselextrasconverter {
  Produces<aod::StraEvSelExtras_001> straEvSelExtras_001;

  void init(InitContext&)
  {
    LOGF(info, "Initializing now: cross-checking correctness...");
    if (doprocessAll + doprocessStraEvSelsOnly > 1) {
      LOGF(fatal, "You have enabled more than one process function. Please check your configuration! Aborting now.");
    }
  }

  void processAll(soa::Join<aod::StraEvSels_005, aod::StraEvSelExtras_000> const& straEvSels_005)
  {
    straEvSelExtras_001.reserve(straEvSels_005.size());
    for (const auto& values : straEvSels_005) {
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

  void processStraEvSelsOnly(aod::StraEvSels_005 const& straEvSels_005)
  {
    straEvSelExtras_001.reserve(straEvSels_005.size());
    for (const auto& values : straEvSels_005) {
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
                          -999., // dummy timeFDDA,
                          -999., // dummy timeFDDC,
                          -999., // dummy timeFV0A,
                          -999., // dummy timeFT0A,
                          -999., // dummy timeFT0C,
                          0,     // dummy triggerMaskFT0,
                          values.gapSide(),
                          values.totalFT0AmplitudeA(),
                          values.totalFT0AmplitudeC(),
                          values.totalFV0AmplitudeA(),
                          values.totalFDDAmplitudeA(),
                          values.totalFDDAmplitudeC(),
                          -999., // dummy timeZNA,
                          -999., // dummy timeZNC,
                          values.energyCommonZNA(),
                          values.energyCommonZNC());
    }
  }
  PROCESS_SWITCH(straevselextrasconverter, processAll, "Store StraEvSels005 and StraEvSelExtras000 into StraEvSelExtras001", true);
  PROCESS_SWITCH(straevselextrasconverter, processStraEvSelsOnly, "Store only StraEvSels005 into StraEvSelExtras001. Other columns are set to dummy values", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<straevselextrasconverter>(cfgc)};
}
