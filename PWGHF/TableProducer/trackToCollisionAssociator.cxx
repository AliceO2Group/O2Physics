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
  Produces<HfTrackCompColls> reverseIndices;

  Configurable<float> nSigmaForTimeCompat{"nSigmaForTimeCompat", 3.f, "number of sigmas for time compatibility"};

  using TracksWithSel = soa::Join<Tracks, TracksExtra, TrackSelection>;
  Filter trackFilter = requireGlobalTrackWoDCAInFilter();

  void processAssocWithTime(Collisions const& collisions,
                            TracksWithSel const& tracksUnfiltered,
                            soa::Filtered<TracksWithSel> const& tracks,
                            BCs const& bcs)
  {
    std::vector<int>** collsPerTrack = new std::vector<int>*[tracksUnfiltered.size()];
    memset(collsPerTrack, 0x0, sizeof(collsPerTrack) * tracksUnfiltered.size());

    // loop over collisions to find time-compatible tracks
    const float nSigmaForTimeCompat2 = nSigmaForTimeCompat * nSigmaForTimeCompat;
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
        float trackTime{0.};
        float trackTimeRes2{0.};
        if (track.isPVContributor()) {
          trackTime = track.collision().collisionTime(); // if PV contributor, we assume the time to be the one of the collision
          trackTimeRes2 = 25.f;                          // 1 BC
        } else {
          trackTime = track.trackTime();
          trackTimeRes2 = track.trackTimeRes() * track.trackTimeRes();
        }
        float deltaTime = trackTime - collTime + bcOffset * 25.f;
        float sigmaTimeRes2 = collTimeRes2 + trackTimeRes2;
        LOGP(debug, "collision time={}, collision time res={}, track time={}, track time res={}, bc collision={}, bc track={}, delta time={}", collTime, collision.collisionTimeRes(), track.trackTime(), track.trackTimeRes(), collBC, track.collision().bc().globalBC(), deltaTime);
        if (!iteratorMoved && bcOffset > -bOffsetMax) {
          trackBegin.setCursor(track.filteredIndex());
          iteratorMoved = true;
          LOGP(debug, "Moving iterator begin {}", track.globalIndex());
        } else if (bcOffset > bOffsetMax) {
          LOGP(debug, "Stopping iterator {}", track.globalIndex());
          break;
        }
        if (deltaTime * deltaTime < nSigmaForTimeCompat2 * sigmaTimeRes2) {
          auto collIdx = collision.globalIndex();
          auto trackIdx = track.globalIndex();
          LOGP(debug, "Filling track id {} for coll id {}", trackIdx, collIdx);
          if (collsPerTrack[trackIdx] == nullptr) {
            collsPerTrack[trackIdx] = new std::vector<int>;
          }
          collsPerTrack[trackIdx]->push_back(collIdx);
          association(collIdx, trackIdx);
        }
      }
    }

    // create reverse index track to collisions
    std::vector<int> empty{};
    for (const auto& track : tracksUnfiltered) {
      LOGP(debug, "Track id {}", track.globalIndex());
      if (collsPerTrack[track.globalIndex()] == nullptr) {
        reverseIndices(empty);
      } else {
        // for (const auto& collId : *collsPerTrack[track.globalIndex()]) {
        //   LOGP(info, "  -> Coll id {}", collId);
        // }
        reverseIndices(*collsPerTrack[track.globalIndex()]);
      }
    }

    for (int iTrack{0u}; iTrack < tracksUnfiltered.size(); ++iTrack) {
      delete collsPerTrack[iTrack];
    }
    delete[] collsPerTrack;
  }

  PROCESS_SWITCH(HfTrackToCollisionAssociation, processAssocWithTime, "Use track-to-collision association based on time", false);

  Preslice<Tracks> perCollision = o2::aod::track::collisionId;

  void processStandardAssoc(Collisions const& collisions,
                            Tracks const& tracks)
  {
    for (const auto& collision : collisions) { // we do it for all tracks, to be compatible with Run2 analyses
      uint64_t collIdx = collision.globalIndex();
      auto tracksPerCollision = tracks.sliceBy(perCollision, collIdx);
      for (const auto& track : tracksPerCollision) {
        association(collIdx, track.globalIndex());
      }
    }

    std::vector<int> empty{};
    for (const auto& track : tracks) {
      if (!track.has_collision()) {
        reverseIndices(empty);
      } else {
        reverseIndices(std::vector<int>{track.collisionId()});
      }
    }
  }

  PROCESS_SWITCH(HfTrackToCollisionAssociation, processStandardAssoc, "Use standard track-to-collision association", true);
};

//________________________________________________________________________________________________________________________
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfTrackToCollisionAssociation>(cfgc));

  return workflow;
}
