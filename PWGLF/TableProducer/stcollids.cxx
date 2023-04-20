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

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct STCollIdTask {
  Produces<aod::TrackedCascadeColls> trackedCascadeColls;

  void init(InitContext const&) {}

  void processTrackedCascades(aod::TrackedCascades const& trackedCascades, aod::Tracks const& tracks)
  {
    for (const auto& trackedCascade : trackedCascades)
      trackedCascadeColls(trackedCascade.track().collisionId());
  }
  PROCESS_SWITCH(STCollIdTask, processTrackedCascades, "process cascades from strangeness tracking", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<STCollIdTask>(cfgc)};
}