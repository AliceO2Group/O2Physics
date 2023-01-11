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

/// \file trackToCollisionAssociator.cxx
/// \brief Associates tracks to collisions considering ambiguities
///
/// \author

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

struct HfTrackToCollisionAssociation {
  Produces<HfTrackAssoc> association;

  using TracksWithSel = soa::Join<Tracks, TracksExtra, TrackSelection>;
  Filter trackFilter = requireGlobalTrackWoDCAInFilter();

  void process(Collisions const& collisions,
               soa::Filtered<TracksWithSel> const& tracks,
               BCs const& bcs)
  {
    auto trackBegin = tracks.begin();
    const auto bOffsetMax = 241;
    for (const auto& collision : collisions) {
      float collTime = collision.collisionTime();
      float collTimeRes2 = collision.collisionTimeRes() * collision.collisionTimeRes();
      uint64_t collBC = collision.bc().globalBC();
      bool iteratorMoved = false;
      for (auto track = trackBegin; track != tracks.end(); ++track) {
        if (!track.has_collision()) {
          continue;
        }
        int64_t bcOffset = (int64_t)track.collision().bc().globalBC() - (int64_t)collBC;
        float deltaTime = track.trackTime() - collTime + bcOffset * 25.f;
        LOGP(debug, "collision time={}, collision time res={}, track time={}, track time res={}, bc collision={}, bc track={}, delta time={}", collTime, collision.collisionTimeRes(), track.trackTime(), track.trackTimeRes(), collBC, track.collision().bc().globalBC(), deltaTime);
        if (!iteratorMoved && bcOffset > -bOffsetMax) {
          trackBegin.setCursor(track.filteredIndex());
          iteratorMoved = true;
          LOGP(debug, "Moving iterator begin {}", track.globalIndex());
        } else if (bcOffset > bOffsetMax) {
          LOGP(debug, "Stopping iterator {}", track.globalIndex());
          break;
        }
        float sigmaTimeRes2 = collTimeRes2 + track.trackTimeRes() * track.trackTimeRes();
        if (deltaTime * deltaTime < 9 * sigmaTimeRes2) {
          LOGP(debug, "Filling track id {} for coll id {}", track.globalIndex(), collision.globalIndex());
          association(collision.globalIndex(), track.globalIndex());
        }
      }
    }
    LOGP(info, "Table association size {}, table track size {}", association.lastIndex(), tracks.size());
  }
};

//________________________________________________________________________________________________________________________
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfTrackToCollisionAssociation>(cfgc));

  return workflow;
}
