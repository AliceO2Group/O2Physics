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

#include "Common/DataModel/CollisionAssociation.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/LHCConstants.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

struct TrackToCollisionAssociation {
  Produces<TrackAssoc> association;

  Configurable<float> nSigmaForTimeCompat{"nSigmaForTimeCompat", 4.f, "number of sigmas for time compatibility"};
  Configurable<float> timeMargin{"timeMargin", 0.f, "time margin in ns added to uncertainty because of uncalibrated TPC"};
  Configurable<bool> applyMinimalTrackSelForRun2{"applyMinimalTrackSelForRun2", false, "flag to apply minimal track selection for Run 2 in case of standard association"};
  Configurable<bool> applyIsGlobalTrackWoDCA{"applyIsGlobalTrackWoDCA", true, "flag to apply global track w/o DCA selection"};
  Configurable<bool> usePVAssociation{"usePVAssociation", true, "if the track is a PV contributor, use the collision time for it"};
  Configurable<bool> includeUnassigned{"includeUnassigned", false, "consider also tracks which are not assigned to any collision"};

  Filter trackFilter = (applyIsGlobalTrackWoDCA == false) || requireGlobalTrackWoDCAInFilter();
  using TracksWithSel = soa::Join<Tracks, TracksExtra, TrackSelection>;
  using TracksWithSelFilter = soa::Filtered<TracksWithSel>;

  void init(InitContext const&)
  {
    if (doprocessAssocWithTime == doprocessStandardAssoc) {
      LOGP(fatal, "Exactly one process function should be enabled!");
    }

    if (applyIsGlobalTrackWoDCA && applyMinimalTrackSelForRun2) {
      LOGP(fatal, "You cannot apply simultaneously Run 2 track selections and applyIsGlobalTrackWoDCA!");
    }
  }

  void processAssocWithTime(Collisions const& collisions,
                            TracksWithSelFilter const& tracks,
                            AmbiguousTracks const& ambiguousTracks,
                            BCs const& bcs)
  {
    // cache globalBC
    std::vector<uint64_t> globalBC;
    for (const auto& track : tracks) {
      if (track.has_collision()) {
        globalBC.push_back(track.collision().bc().globalBC());
      } else {
        for (const auto& ambTrack : ambiguousTracks) {
          if (ambTrack.trackId() == track.globalIndex()) {
            globalBC.push_back(ambTrack.bc().begin().globalBC());
            break;
          }
        }
      }
    }

    // loop over collisions to find time-compatible tracks
    auto trackBegin = tracks.begin();
    constexpr auto bOffsetMax = 241; // 6 mus (ITS)
    for (const auto& collision : collisions) {
      const float collTime = collision.collisionTime();
      const float collTimeRes2 = collision.collisionTimeRes() * collision.collisionTimeRes();
      uint64_t collBC = collision.bc().globalBC();
      // bool iteratorMoved = false;
      for (auto track = trackBegin; track != tracks.end(); ++track) {
        if (!includeUnassigned && !track.has_collision()) {
          continue;
        }
        const int64_t bcOffset = (int64_t)globalBC[track.filteredIndex()] - (int64_t)collBC;
        if (std::abs(bcOffset) > bOffsetMax) {
          continue;
        }

        float trackTime{0.};
        float trackTimeRes{0.};
        if (usePVAssociation && track.isPVContributor()) {
          trackTime = track.collision().collisionTime();    // if PV contributor, we assume the time to be the one of the collision
          trackTimeRes = constants::lhc::LHCBunchSpacingNS; // 1 BC
        } else {
          trackTime = track.trackTime();
          trackTimeRes = track.trackTimeRes();
        }
        const float deltaTime = trackTime - collTime + bcOffset * constants::lhc::LHCBunchSpacingNS;
        float sigmaTimeRes2 = collTimeRes2 + trackTimeRes * trackTimeRes;
        LOGP(debug, "collision time={}, collision time res={}, track time={}, track time res={}, bc collision={}, bc track={}, delta time={}", collTime, collision.collisionTimeRes(), track.trackTime(), track.trackTimeRes(), collBC, globalBC[track.filteredIndex()], deltaTime);

        // optimization to avoid looping over the full track list each time. This assumes that tracks are sorted by BCs (which they should be because collisions are sorted by BCs)
        // NOTE this does not work anymore if includeUnassigned is set as the unassigned blocks can be somewhere (and we can have merged DFs, too)
        // if (!iteratorMoved && bcOffset > -bOffsetMax) {
        //   trackBegin.setCursor(track.filteredIndex());
        //   iteratorMoved = true;
        //   LOGP(debug, "Moving iterator begin {}", track.globalIndex());
        // } else if (bcOffset > bOffsetMax) {
        //   LOGP(debug, "Stopping iterator {}", track.globalIndex());
        //   break;
        // }

        float thresholdTime = 0.;
        if (usePVAssociation && track.isPVContributor()) {
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
        }
      }
    }
  }

  PROCESS_SWITCH(TrackToCollisionAssociation, processAssocWithTime, "Use track-to-collision association based on time", false);

  void processStandardAssoc(Collision const& collision,
                            TracksWithSel const& tracksPerCollision)
  {
    // we do it for all tracks, to be compatible with Run 2 analyses
    for (const auto& track : tracksPerCollision) {
      bool hasGoodQuality = true;
      if (applyIsGlobalTrackWoDCA && !track.isGlobalTrackWoDCA()) {
        hasGoodQuality = false;
      } else if (applyMinimalTrackSelForRun2) {
        unsigned char itsClusterMap = track.itsClusterMap();
        if (!(track.tpcNClsFound() >= 50 && track.flags() & o2::aod::track::ITSrefit && track.flags() & o2::aod::track::TPCrefit && (TESTBIT(itsClusterMap, 0) || TESTBIT(itsClusterMap, 1)))) {
          hasGoodQuality = false;
        }
      }
      if (hasGoodQuality) {
        association(collision.globalIndex(), track.globalIndex());
      }
    }
  }

  PROCESS_SWITCH(TrackToCollisionAssociation, processStandardAssoc, "Use standard track-to-collision association", true);
};

//________________________________________________________________________________________________________________________
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TrackToCollisionAssociation>(cfgc)};
}
