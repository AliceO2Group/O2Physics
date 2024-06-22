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
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

using namespace o2;
using namespace o2::framework;

// Converts EMEvents version 000 to 001
struct emeventsconverter {
  Produces<aod::EMEvents_001> events_001;

  void process(aod::EMEvents_000 const& events)
  {
    for (auto& event : events) {
      events_001(
        event.collisionId(),
        event.runNumber(),
        0,
        event.sel8(),
        event.alias_raw(),
        event.selection_raw(),
        0,
        event.ncollsPerBC(),
        event.posX(),
        event.posY(),
        event.posZ(),
        event.numContrib(),
        0 /*dummy occupancy value*/);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<emeventsconverter>(cfgc, TaskName{"emevents-converter"})};
}
