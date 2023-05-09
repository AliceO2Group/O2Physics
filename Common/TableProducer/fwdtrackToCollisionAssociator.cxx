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

/// \file fwdtrackToCollisionAssociator.cxx
/// \brief Associates fwd and MFT tracks to collisions considering ambiguities
/// \author Sarah Herrmann <sarah.herrmann@cern.ch>, IP2I Lyon
/// \author Maurice Coquet <maurice.louis.coquet@cern.ch>, CEA-Saclay/Irfu

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

struct FwdTrackToCollisionAssociation {
  Produces<FwdTrackAssoc> fwdassociation;
  Produces<FwdTrkCompColls> fwdreverseIndices;
  Produces<MFTTrackAssoc> mftassociation;
  Produces<MFTTrkCompColls> mftreverseIndices;

  Configurable<float> nSigmaForTimeCompat{"nSigmaForTimeCompat", 4.f, "number of sigmas for time compatibility"};
  Configurable<float> timeMargin{"timeMargin", 0.f, "time margin in ns added to uncertainty because of uncalibrated TPC"};
  Configurable<bool> includeUnassigned{"includeUnassigned", false, "consider also tracks which are not assigned to any collision"};
  Configurable<bool> fillTableOfCollIdsPerTrack{"fillTableOfCollIdsPerTrack", false, "fill additional table with vector of collision ids per track"};

  Preslice<FwdTracks> muonsPerCollisions = aod::fwdtrack::collisionId;
  Preslice<MFTTracks> mftsPerCollisions = aod::fwdtrack::collisionId;

  void init(InitContext const&)
  {
  }

  template <typename TTracks, typename Slice, typename Assoc, typename RevIndices>
  void runStandardAssoc(Collisions const& collisions,
                        TTracks const& tracks, Slice perCollisions, Assoc association, RevIndices reverseIndices)
  {
    // we do it for all tracks, to be compatible with Run 2 analyses
    for (const auto& collision : collisions) {
      auto tracksThisCollision = tracks.sliceBy(perCollisions, collision.globalIndex());
      for (const auto& track : tracksThisCollision) {
        association(collision.globalIndex(), track.globalIndex());
      }
    }

    // create reverse index track to collisions if enabled
    std::vector<int> empty{};
    if (fillTableOfCollIdsPerTrack) {
      for (const auto& track : tracks) {
        if (track.has_collision()) {
          reverseIndices(std::vector<int>{track.collisionId()});
        } else {
          reverseIndices(empty);
        }
      }
    }
  }

  template <typename TTracks, typename TAmbiTracks, typename Assoc, typename RevIndices>
  void runAssocWithTime(Collisions const& collisions,
                        TTracks const& tracks,
                        TAmbiTracks const& ambiguousTracks,
                        BCs const& bcs,
                        Assoc association, RevIndices reverseIndices)
  {
    // cache globalBC
    std::vector<uint64_t> globalBC;
    for (const auto& track : tracks) {
      if (track.has_collision()) {
        globalBC.push_back(track.collision().bc().globalBC());
      } else {
        for (const auto& ambTrack : ambiguousTracks) {
          if (ambTrack.template getId<TTracks>() == track.globalIndex()) {
            globalBC.push_back(ambTrack.bc().begin().globalBC());
            break;
          }
        }
      }
    }

    // define vector of vectors to store indices of compatible collisions per track
    std::vector<std::unique_ptr<std::vector<int>>> collsPerTrack(tracks.size());

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
        float trackTime = track.trackTime();
        float trackTimeRes = track.trackTimeRes();

        const float deltaTime = trackTime - collTime + bcOffset * constants::lhc::LHCBunchSpacingNS;
        float sigmaTimeRes2 = collTimeRes2 + trackTimeRes * trackTimeRes;
        LOGP(debug, "collision time={}, collision time res={}, track time={}, track time res={}, bc collision={}, bc track={}, delta time={}", collTime, collision.collisionTimeRes(), track.trackTime(), track.trackTimeRes(), collBC, globalBC[track.filteredIndex()], deltaTime);

        float thresholdTime = 0.;

        thresholdTime = nSigmaForTimeCompat * std::sqrt(sigmaTimeRes2) + timeMargin;

        if (std::abs(deltaTime) < thresholdTime) {
          const auto collIdx = collision.globalIndex();
          const auto trackIdx = track.globalIndex();
          LOGP(debug, "Filling track id {} for coll id {}", trackIdx, collIdx);
          association(collIdx, trackIdx);
          if (fillTableOfCollIdsPerTrack) {
            if (collsPerTrack[trackIdx] == nullptr) {
              collsPerTrack[trackIdx] = std::make_unique<std::vector<int>>();
            }
            collsPerTrack[trackIdx].get()->push_back(collIdx);
          }
        }
      }
    }
    // create reverse index track to collisions if enabled
    if (fillTableOfCollIdsPerTrack) {
      std::vector<int> empty{};
      for (const auto& track : tracks) {

        const auto trackId = track.globalIndex();
        if (collsPerTrack[trackId] == nullptr) {
          reverseIndices(empty);
        } else {
          reverseIndices(*collsPerTrack[trackId].get());
        }
      }
    }
  }

  void processFwdAssocWithTime(Collisions const& collisions,
                               FwdTracks const& muons,
                               AmbiguousFwdTracks const& ambiTracksFwd,
                               BCs const& bcs)
  {
    runAssocWithTime(collisions, muons, ambiTracksFwd, bcs, fwdassociation, fwdreverseIndices);
  }
  PROCESS_SWITCH(FwdTrackToCollisionAssociation, processFwdAssocWithTime, "Use fwdtrack-to-collision association based on time", true);

  void processFwdStandardAssoc(Collisions const& collisions,
                               FwdTracks const& muons)
  {
    runStandardAssoc(collisions, muons, muonsPerCollisions, fwdassociation, fwdreverseIndices);
  }
  PROCESS_SWITCH(FwdTrackToCollisionAssociation, processFwdStandardAssoc, "Use standard fwdtrack-to-collision association", false);

  void processMFTAssocWithTime(Collisions const& collisions,
                               MFTTracks const& tracks,
                               AmbiguousMFTTracks const& ambiguousTracks,
                               BCs const& bcs)
  {
    runAssocWithTime(collisions, tracks, ambiguousTracks, bcs, mftassociation, mftreverseIndices);
  }
  PROCESS_SWITCH(FwdTrackToCollisionAssociation, processMFTAssocWithTime, "Use MFTtrack-to-collision association based on time", true);

  void processMFTStandardAssoc(Collisions const& collisions,
                               MFTTracks const& tracks)
  {
    runStandardAssoc(collisions, tracks, mftsPerCollisions, mftassociation, mftreverseIndices);
  }
  PROCESS_SWITCH(FwdTrackToCollisionAssociation, processMFTStandardAssoc, "Use standard mfttrack-to-collision association", false);
};

//________________________________________________________________________________________________________________________
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<FwdTrackToCollisionAssociation>(cfgc)};
}
