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
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "Common/DataModel/Multiplicity.h"

using namespace o2;
using namespace o2::framework;

struct MultDqExtraConverter {

  Produces<aod::ReducedEventsMultPV_001> reducedEventsMultPV_001;
  void process(aod::ReducedEventsMultPV_000 const& reducedEventsMultPV_000)
  {
    for (const auto& r : reducedEventsMultPV_000) {
      reducedEventsMultPV_001(r.multNTracksHasITS(),
                              r.multNTracksHasTPC(), r.multNTracksHasTOF(), r.multNTracksHasTRD(), r.multNTracksITSOnly(),
                              r.multNTracksTPCOnly(), r.multNTracksITSTPC(), 0.0f, 0.0f, 0.0f,
                              r.trackOccupancyInTimeRange());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MultDqExtraConverter>(cfgc)};
}
