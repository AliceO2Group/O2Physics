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

/// \file stCollIds.cxx
/// \brief Task to create joinable tables with collision IDs of strangeness-tracked objects
///
/// \author Jochen Klein

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct StCollIds {
  Produces<aod::TrackedCascadeColls> trackedCascadeColls;
  Produces<aod::TrackedV0Colls> trackedV0Colls;
  Produces<aod::Tracked3BodyColls> tracked3BodyColls;

  void init(InitContext const&) {}

  void processTrackedCascades(aod::TrackedCascades const& trackedCascades, aod::TracksIU const& /*tracks*/)
  {
    for (const auto& trackedCascade : trackedCascades)
      trackedCascadeColls(trackedCascade.track_as<aod::TracksIU>().collisionId());
  }
  PROCESS_SWITCH(StCollIds, processTrackedCascades, "process cascades from strangeness tracking", true);

  void processTrackedV0s(aod::TrackedV0s const& trackedV0s, aod::TracksIU const& /*tracks*/)
  {
    for (const auto& trackedV0 : trackedV0s)
      trackedV0Colls(trackedV0.track_as<aod::TracksIU>().collisionId());
  }
  PROCESS_SWITCH(StCollIds, processTrackedV0s, "process V0s from strangeness tracking", true);

  void processTracked3Bodys(aod::Tracked3Bodys const& tracked3Bodys, aod::TracksIU const& /*tracks*/)
  {
    for (const auto& tracked3Body : tracked3Bodys)
      tracked3BodyColls(tracked3Body.track_as<aod::TracksIU>().collisionId());
  }
  PROCESS_SWITCH(StCollIds, processTracked3Bodys, "process cascades from strangeness tracking", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<StCollIds>(cfgc)};
}
