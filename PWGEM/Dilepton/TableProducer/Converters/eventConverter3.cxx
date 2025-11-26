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
// ========================
//
// This code runs loop over ULS ee pars for virtual photon QC.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct eventConverter3 {
  Produces<aod::EMEvents_003> event_003;

  void process(aod::EMEvents_002 const& collisions)
  {
    for (const auto& collision : collisions) {
      event_003(
        collision.globalIndex(),
        collision.runNumber(),
        collision.globalBC(),
        collision.alias_raw(),
        collision.selection_raw(),
        collision.rct_raw(),
        collision.timestamp(),
        collision.posZ(),
        collision.numContrib(),
        collision.trackOccupancyInTimeRange(),
        collision.ft0cOccupancyInTimeRange());
    } // end of collision loop
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<eventConverter3>(cfgc, TaskName{"event-converter3"})};
}
