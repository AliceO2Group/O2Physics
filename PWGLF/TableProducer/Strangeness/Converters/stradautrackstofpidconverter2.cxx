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
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;

// converts DauTrackTOFPIDs_000 to _001
struct stradautrackstofpidconverter2 {
  Produces<aod::DauTrackTOFPIDs_001> dautracktofpids;
  Produces<aod::StraEvTimes> straEvTimes;

  void process(aod::StraCollisions const& collisions, soa::Join<aod::DauTrackExtras, aod::DauTrackTOFPIDs_000> const& dauTracks, soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras> const& v0s)
  {
    // create new TOFPIDs
    for (int ii = 0; ii < dauTracks.size(); ii++) {
      auto dauTrack = dauTracks.rawIteratorAt(ii);
      dautracktofpids(-1, -1, dauTrack.tofSignal(), dauTrack.tofEvTime(), dauTrack.length(), 0.0f);
    }

    // fill EvTimes (created simultaneously, so done in the same converter)
    // recover as much as possible from the previous format
    std::vector<double> collisionEventTime(collisions.size(), 0.0);
    std::vector<int> collisionNtracks(collisions.size(), 0);

    for (const auto& v0 : v0s) {
      auto posTrackTOF = dauTracks.rawIteratorAt(v0.posTrackExtraId());
      auto negTrackTOF = dauTracks.rawIteratorAt(v0.negTrackExtraId());
      if (posTrackTOF.hasTOF()) {
        collisionEventTime[v0.straCollisionId()] += posTrackTOF.tofEvTime();
        collisionNtracks[v0.straCollisionId()]++;
      }
      if (negTrackTOF.hasTOF()) {
        collisionEventTime[v0.straCollisionId()] += negTrackTOF.tofEvTime();
        collisionNtracks[v0.straCollisionId()]++;
      }
    }
    for (const auto& collision : collisions) {
      if (collisionNtracks[collision.globalIndex()] > 0) {
        collisionEventTime[collision.globalIndex()] /= static_cast<double>(collisionNtracks[collision.globalIndex()]);
      } else {
        collisionEventTime[collision.globalIndex()] = -1e+6; // undefined
      }
      straEvTimes(collisionEventTime[collision.globalIndex()]);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<stradautrackstofpidconverter2>(cfgc)};
}
