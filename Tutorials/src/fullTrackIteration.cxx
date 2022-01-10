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
///
/// \brief Joined tables can be used as argument to the process function.
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

using namespace o2;
using namespace o2::framework;

struct UseJoins {
  void process(soa::Join<aod::Tracks, aod::TracksExtra> const& fullTracks)
  {
    for (auto& track : fullTracks) {
      LOGF(info, "%d, %f %f", track.globalIndex(), track.alpha(), track.tpcSignal());
    }
  }
};

struct LoopAmbiguousTracks {
  void process(aod::AmbiguousTracks const& tracks, aod::BCs const&)
  {
    for (auto& track : tracks) {
      LOGF(info, "We look at track %d which has %d possible BCs", track.globalIndex(), track.bc().size());
      for (auto& bc : track.bc()) {
        LOGF(info, "  BC %d with global BC %lld", bc.globalIndex(), bc.globalBC());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UseJoins>(cfgc),
    adaptAnalysisTask<LoopAmbiguousTracks>(cfgc),
  };
}
