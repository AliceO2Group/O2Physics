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
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;

// Converts V0 and cascade version 000 to 001
// Build indices to group V0s and cascades to collisions

struct WeakDecayIndicesV0 {
  Produces<aod::V0s_001> v0s_001;

  void process(aod::V0s_000 const& v0s, aod::Tracks const&)
  {
    for (auto& v0 : v0s) {
      if (v0.posTrack().collisionId() != v0.negTrack().collisionId()) {
        LOGF(fatal, "V0 %d has inconsistent collision information (%d, %d)", v0.globalIndex(), v0.posTrack().collisionId(), v0.negTrack().collisionId());
      }
      v0s_001(v0.posTrack().collisionId(), v0.posTrackId(), v0.negTrackId());
    }
  }
};

// NOTE These tasks have to be split because for the cascades, V0s and not V0s_000 are needed
struct WeakDecayIndicesCascades {
  Produces<aod::Cascades_001> cascades_001;

  void process(aod::V0s const&, aod::Cascades_000 const& cascades, aod::Tracks const&)
  {
    for (auto& cascade : cascades) {
      if (cascade.bachelor().collisionId() != cascade.v0().posTrack().collisionId() || cascade.v0().posTrack().collisionId() != cascade.v0().negTrack().collisionId()) {
        LOGF(fatal, "Cascade %d has inconsistent collision information (%d, %d, %d) track ids %d %d %d", cascade.globalIndex(), cascade.bachelor().collisionId(),
             cascade.v0().posTrack().collisionId(), cascade.v0().negTrack().collisionId(), cascade.bachelorId(), cascade.v0().posTrackId(), cascade.v0().negTrackId());
      }
      cascades_001(cascade.bachelor().collisionId(), cascade.v0Id(), cascade.bachelorId());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<WeakDecayIndicesV0>(cfgc),
    adaptAnalysisTask<WeakDecayIndicesCascades>(cfgc),
  };
}
