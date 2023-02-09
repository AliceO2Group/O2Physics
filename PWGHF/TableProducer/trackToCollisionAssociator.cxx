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
/// \author Jan Fiete Grosse-Oetringhaus <jan.fiete.grosse-oetringhaus@cern.ch>, CERN
/// \author Fabrizio Grosa <fgrosa@cern.ch>, CERN
/// \author Mattia Faggin <mfaggin@cern.ch>, University and INFN Padova

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

  Configurable<float> nSigmaForTimeCompat{"nSigmaForTimeCompat", 4.f, "number of sigmas for time compatibility"};
  Configurable<float> timeMargin{"timeMargin", 0.f, "time margin in ns added to uncertainty because of uncalibrated TPC"};
  Configurable<bool> applyTrackSelForRun2{"applyMinimalTrackSelForRun2", false, "flag to apply minimal track selection for Run2 in case of standard association"};
  Configurable<bool> applyIsGlobalTrack{"applyIsGlobalTrack", true, "flag to apply global track w/o DCA selection"};
  Configurable<bool> debug{"debug", false, "fill a table with flag to keep track of kind of track (PV contributor or ambiguous)"};

  Filter trackFilter = (applyIsGlobalTrack == false) || requireGlobalTrackWoDCAInFilter();
  using TracksWithSel = soa::Join<Tracks, TracksExtra, TrackSelection>;
  using TracksWithSelFilter = soa::Filtered<TracksWithSel>;

  void init(InitContext const&)
  {
    std::array<int, 3> doProcess = {doprocessAssocWithTime, doprocessAssocWithAmb, doprocessStandardAssoc};
    int doProcessSum = std::accumulate(doProcess.begin(), doProcess.end(), 0);
    if (doProcessSum != 1) {
      LOGP(fatal, "Exactly one process function should be enabled! Exit");
    }

    if (applyIsGlobalTrack && applyTrackSelForRun2) {
      LOGP(fatal, "You cannot apply simultaneously Run2 track selections and isGlobalWoDCA! Exit");
    }
  }

  void processAssocWithTime(Collisions const& collisions,
                            TracksWithSel const& tracksUnfiltered,
                            TracksWithSelFilter const& tracks,
                            BCs const& bcs)
  {
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
        int trackType = hf_track_association::TrackTypes::Ambiguous;
        if (track.isPVContributor()) {
          trackTime = track.collision().collisionTime(); // if PV contributor, we assume the time to be the one of the collision
          trackTimeRes = 25.f;                           // 1 BC
          trackType = hf_track_association::TrackTypes::PVContributor;
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
          association(collIdx, trackIdx);
          if (debug) {
            associationExtra(trackType);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(HfTrackToCollisionAssociation, processAssocWithTime, "Use track-to-collision association based on time", false);

  Partition<TracksWithSel> pvContributors = ((aod::track::flags & (uint32_t)aod::track::PVContributor) == (uint32_t)aod::track::PVContributor);

  void processAssocWithAmb(Collisions const& collisions,
                           TracksWithSel const& tracksUnfiltered,
                           TracksWithSelFilter const& tracks,
                           AmbiguousTracks const& ambTracks,
                           BCs const& bcs)
  {
    std::vector<uint64_t> collIds{}, trackIds{}, trackTypes{};

    // associate collisions with tracks as default in the AO2Ds
    for (const auto& track : tracks) {
      if (!track.has_collision()) {
        continue;
      }
      const auto trackId = track.globalIndex();
      const auto collId = track.collisionId();
      collIds.push_back(collId);
      trackIds.push_back(trackId);
      trackTypes.push_back(hf_track_association::TrackTypes::Regular);
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
            trackTypes.push_back(hf_track_association::TrackTypes::Ambiguous);
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
                collIds.push_back(repCollId2);
                trackIds.push_back(trackId);
                trackTypes.push_back(hf_track_association::TrackTypes::PVContributor);
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
      if (debug) {
        associationExtra(trackTypes[index]);
      }
    }
  }

  PROCESS_SWITCH(HfTrackToCollisionAssociation, processAssocWithAmb, "Use track-to-collision association based on ambiguous tracks", false);

  Preslice<Tracks> perCollision = o2::aod::track::collisionId;

  void processStandardAssoc(Collisions const& collisions,
                            TracksWithSel const& tracks)
  {
    for (const auto& collision : collisions) { // we do it for all tracks, to be compatible with Run2 analyses
      const uint64_t collIdx = collision.globalIndex();
      auto tracksPerCollision = tracks.sliceBy(perCollision, collIdx);
      for (const auto& track : tracksPerCollision) {
        bool hasGoodQuality = true;
        if (applyIsGlobalTrack && !track.isGlobalTrackWoDCA()) {
          hasGoodQuality = false;
        } else if (applyTrackSelForRun2) {
          unsigned char itsClusterMap = track.itsClusterMap();
          if (!(track.tpcNClsFound() >= 50 && track.flags() & o2::aod::track::ITSrefit && track.flags() & o2::aod::track::TPCrefit && (TESTBIT(itsClusterMap, 0) || TESTBIT(itsClusterMap, 1)))) {
            hasGoodQuality = false;
          }
        }
        if (hasGoodQuality) {
          association(collIdx, track.globalIndex());
        }
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
