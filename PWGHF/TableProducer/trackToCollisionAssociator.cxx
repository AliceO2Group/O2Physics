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
  Produces<HfTrackAssocExtra> associationExtra;
  Produces<HfCompColls> reverseIndices;

  Configurable<float> nSigmaForTimeCompat{"nSigmaForTimeCompat", 3.f, "number of sigmas for time compatibility"};
  Configurable<float> timeMargin{"timeMargin", 0., "time margin in ns added to uncertainty because of uncalibrated TPC"};

  Filter trackFilter = requireGlobalTrackWoDCAInFilter();
  using TracksWithSel = soa::Join<Tracks, TracksExtra, TrackSelection>;
  using TracksWithSelFilter = soa::Filtered<TracksWithSel>;

  void processAssocWithTime(Collisions const& collisions,
                            TracksWithSel const& tracksUnfiltered,
                            TracksWithSelFilter const& tracks,
                            BCs const& bcs)
  {
    std::vector<std::unique_ptr<std::vector<int>>> collsPerTrack(tracksUnfiltered.size());

    // loop over collisions to find time-compatible tracks
    auto trackBegin = tracks.begin();
    const auto bOffsetMax = 241; // 6 mus (ITS)
    for (const auto& collision : collisions) {
      const float collTime = collision.collisionTime();
      const float collTimeRes2 = collision.collisionTimeRes() * collision.collisionTimeRes();
      uint64_t collBC = collision.bc().globalBC();
      bool iteratorMoved = false;
      for (auto track = trackBegin; track != tracks.end(); ++track) {
        if (!track.has_collision()) {
          continue;
        }
        const int64_t bcOffset = (int64_t)track.collision().bc().globalBC() - (int64_t)collBC;
        float trackTime{0.};
        float trackTimeRes{0.};
        if (track.isPVContributor()) {
          trackTime = track.collision().collisionTime(); // if PV contributor, we assume the time to be the one of the collision
          trackTimeRes = 25.f;                           // 1 BC
        } else {
          trackTime = track.trackTime();
          trackTimeRes = track.trackTimeRes();
        }
        const float deltaTime = trackTime - collTime + bcOffset * 25.f;
        float sigmaTimeRes2 = collTimeRes2 + trackTimeRes * trackTimeRes;
        LOGP(debug, "collision time={}, collision time res={}, track time={}, track time res={}, bc collision={}, bc track={}, delta time={}", collTime, collision.collisionTimeRes(), track.trackTime(), track.trackTimeRes(), collBC, track.collision().bc().globalBC(), deltaTime);
        // optimization to avoid looping over the full track list each time. This assumes that tracks are sorted by BCs (which they should be because collisions are sorted by BCs)
        if (!iteratorMoved && bcOffset > -bOffsetMax) {
          trackBegin.setCursor(track.filteredIndex());
          iteratorMoved = true;
          LOGP(debug, "Moving iterator begin {}", track.globalIndex());
        } else if (bcOffset > bOffsetMax) {
          LOGP(debug, "Stopping iterator {}", track.globalIndex());
          break;
        }

        float thresholdTime = 0.;
        if (track.isPVContributor()) {
          thresholdTime = trackTimeRes;
        } else if (TESTBIT(track.flags(), o2::aod::track::TrackTimeResIsRange)) {
          thresholdTime = std::sqrt(sigmaTimeRes2) + timeMargin;
        } else {
          thresholdTime = nSigmaForTimeCompat * std::sqrt(sigmaTimeRes2) + timeMargin;
        }

        if (std::abs(deltaTime) < thresholdTime) {
          const auto collIdx = collision.globalIndex();
          const auto trackIdx = track.globalIndex();
          LOGP(debug, "Filling track id {} for coll id {}", trackIdx, collIdx);
          if (collsPerTrack[trackIdx] == nullptr) {
            collsPerTrack[trackIdx] = std::make_unique<std::vector<int>>();
          }
          collsPerTrack[trackIdx].get()->push_back(collIdx);
          association(collIdx, trackIdx);
        }
      }
    }

    // create reverse index track to collisions
    std::vector<int> empty{};
    for (const auto& track : tracksUnfiltered) {

      const auto trackId = track.globalIndex();
      if (collsPerTrack[trackId] == nullptr) {
        reverseIndices(empty);
      } else {
        // for (const auto& collId : *collsPerTrack[trackId]) {
        //   LOGP(info, "  -> Coll id {}", collId);
        // }
        reverseIndices(*collsPerTrack[trackId].get());
      }
    }
  }

  PROCESS_SWITCH(HfTrackToCollisionAssociation, processAssocWithTime, "Use track-to-collision association based on time", true);

  Partition<TracksWithSel> pvContributors = ((aod::track::flags & (uint32_t)aod::track::PVContributor) == (uint32_t)aod::track::PVContributor);

  void processAssocWithAmb(Collisions const& collisions,
                           TracksWithSel const& tracksUnfiltered,
                           TracksWithSelFilter const& tracks,
                           AmbiguousTracks const& ambTracks,
                           BCs const& bcs)
  {
    std::vector<uint64_t> collIds{}, trackIds{}, trackTypes{};
    std::vector<std::unique_ptr<std::vector<int>>> collsPerTrack(tracksUnfiltered.size());

    // associate collisions with tracks as default in the AO2Ds
    for (const auto& track : tracks) {
      const auto trackId = track.globalIndex();
      const auto collId = track.collisionId();
      collIds.push_back(collId);
      trackIds.push_back(trackId);
      trackTypes.push_back(hf_track_association::eTrackType::Regular);
      if (collsPerTrack[trackId] == nullptr) {
        collsPerTrack[trackId] = std::make_unique<std::vector<int>>();
      }
      collsPerTrack[trackId].get()->push_back(collId);
    }

    // associate collisions with ambiguous tracks as in the AO2Ds
    for (const auto& ambTrack : ambTracks) {
      auto track = ambTrack.track_as<TracksWithSel>();
      if (!track.isGlobalTrackWoDCA()) {
        continue;
      }

      for (const auto& collision : collisions) {
        const auto collId = collision.globalIndex();
        if (collId == track.collisionId()) {
          continue;
        }
        const uint64_t mostProbableBc = collision.bc().globalBC();
        for (const auto& bc : ambTrack.bc()) {
          if (bc.globalBC() == mostProbableBc) {
            const uint64_t trackId = track.globalIndex();
            collIds.push_back(collId);
            trackIds.push_back(trackId);
            trackTypes.push_back(hf_track_association::eTrackType::Ambiguous);
            if (collsPerTrack[trackId] == nullptr) {
              collsPerTrack[trackId] = std::make_unique<std::vector<int>>();
            }
            collsPerTrack[trackId].get()->push_back(collId);
          }
        }
      }
    }

    // associate collisions with PV contributors in in-bunch pileup events
    for (auto& bc : bcs) {
      std::vector<uint64_t> repeatCollIds{};
      for (const auto& collision : collisions) {
        /// count how many coll. have the current bc as most probable
        if (bc.globalBC() == collision.bc().globalBC()) {
          repeatCollIds.push_back(collision.globalIndex());
        }
      } /// end loop on collisiuons

      /// we do not want bcs with only 1 coll
      if (repeatCollIds.size() < 2) {
        continue;
      }

      for (const auto& repCollId : repeatCollIds) {
        const auto& pvContrPerCollision = pvContributors->sliceByCached(aod::track::collisionId, repCollId);
        for (const auto& pvContr : pvContrPerCollision) {
          if (pvContr.isGlobalTrackWoDCA()) {
            for (const auto& repCollId2 : repeatCollIds) {
              if (repCollId2 != repCollId) {
                const uint64_t trackId = pvContr.globalIndex();
                collIds.push_back(repCollId);
                trackIds.push_back(trackId);
                trackTypes.push_back(hf_track_association::eTrackType::PVContributor);
                if (collsPerTrack[trackId] == nullptr) {
                  collsPerTrack[trackId] = std::make_unique<std::vector<int>>();
                }
                collsPerTrack[trackId].get()->push_back(repCollId);
              }
            }
          }
        }
      }
    }

    std::vector<uint64_t> sortedCollIds(collIds.size());
    std::iota(sortedCollIds.begin(), sortedCollIds.end(), 0);
    std::sort(sortedCollIds.begin(), sortedCollIds.end(), [collIds](uint64_t i, uint64_t j) { return collIds[i] < collIds[j]; });

    for (const auto& index : sortedCollIds) {
      association(collIds[index], trackIds[index]);
      associationExtra(trackTypes[index]);
    }

    // create reverse index track to collisions
    std::vector<int> empty{};
    for (const auto& track : tracksUnfiltered) {

      const auto trackId = track.globalIndex();
      if (collsPerTrack[trackId] == nullptr) {
        reverseIndices(empty);
      } else {
        // for (const auto& collId : *collsPerTrack[trackId]) {
        //   LOGP(info, "  -> Coll id {}", collId);
        // }
        reverseIndices(*collsPerTrack[trackId].get());
      }
    }
  }

  PROCESS_SWITCH(HfTrackToCollisionAssociation, processAssocWithAmb, "Use track-to-collision association based on ambiguous tracks", false);

  Preslice<Tracks> perCollision = o2::aod::track::collisionId;

  void processStandardAssoc(Collisions const& collisions,
                            Tracks const& tracks)
  {
    for (const auto& collision : collisions) { // we do it for all tracks, to be compatible with Run2 analyses
      const uint64_t collIdx = collision.globalIndex();
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

  PROCESS_SWITCH(HfTrackToCollisionAssociation, processStandardAssoc, "Use standard track-to-collision association", false);
};

//________________________________________________________________________________________________________________________
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfTrackToCollisionAssociation>(cfgc));

  return workflow;
}
