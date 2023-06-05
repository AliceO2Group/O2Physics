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

/// \file CollisionAssociation.h
/// \brief Function and table definitions for track to collision associators
/// \author Jan Fiete Grosse-Oetringhaus <jan.fiete.grosse-oetringhaus@cern.ch>, CERN
/// \author Fabrizio Grosa <fgrosa@cern.ch>, CERN
/// \author Mattia Faggin <mfaggin@cern.ch>, University and INFN Padova
/// \author Sarah Herrmann <sarah.herrmann@cern.ch>, IP2I Lyon
/// \author Maurice Coquet <maurice.louis.coquet@cern.ch>, CEA-Saclay/Irfu

#ifndef COMMON_DATAMODEL_COLLISIONASSOCIATION_H_
#define COMMON_DATAMODEL_COLLISIONASSOCIATION_H_

#include "CommonConstants/LHCConstants.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{

namespace track_association
{
enum TrackSelection {
  CentralBarrelRun2 = -1,
  None = 0,
  GlobalTrackWoDCA = 1,
  QualityTracksITS = 2
};

DECLARE_SOA_INDEX_COLUMN(Collision, collision);            //! Collision index
DECLARE_SOA_INDEX_COLUMN(Track, track);                    //! Track index
DECLARE_SOA_ARRAY_INDEX_COLUMN(Collision, compatibleColl); //! Array of collision indices
} // namespace track_association

DECLARE_SOA_TABLE(TrackAssoc, "AOD", "TRACKASSOC", //! Table for track-to-collision association for e.g. HF vertex finding - tracks can appear for several collisions
                  track_association::CollisionId,
                  track_association::TrackId);

DECLARE_SOA_TABLE(TrackCompColls, "AOD", "TRKCOMPCOLLS", //! Table with vectors of collision indices stored per track
                  track_association::CollisionIds);

namespace fwdtrack_association
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);            //! Collision index
DECLARE_SOA_INDEX_COLUMN(FwdTrack, fwdtrack);              //! FwdTrack index
DECLARE_SOA_INDEX_COLUMN(MFTTrack, mfttrack);              //! MFTTrack index
DECLARE_SOA_ARRAY_INDEX_COLUMN(Collision, compatibleColl); //! Array of collision indices
} // namespace fwdtrack_association

DECLARE_SOA_TABLE(FwdTrackAssoc, "AOD", "FWDTRACKASSOC", //! Table for fwdtrack-to-collision association
                  fwdtrack_association::CollisionId,
                  fwdtrack_association::FwdTrackId);

DECLARE_SOA_TABLE(FwdTrkCompColls, "AOD", "FWDTRKCOMPCOLL", //! Table with vectors of collision indices stored per fwdtrack
                  fwdtrack_association::CollisionIds);

DECLARE_SOA_TABLE(MFTTrackAssoc, "AOD", "MFTTRACKASSOC", //! Table for mfttrack-to-collision association
                  fwdtrack_association::CollisionId,
                  fwdtrack_association::MFTTrackId);

DECLARE_SOA_TABLE(MFTTrkCompColls, "AOD", "MFTTRKCOMPCOLL", //! Table with vectors of collision indices stored per mfttrack
                  fwdtrack_association::CollisionIds, o2::soa::Marker<1>);

template <bool isCentralBarrel, typename TTracks, typename Slice, typename Assoc, typename RevIndices>
void runStandardAssoc(Collisions const& collisions,
                      TTracks const& tracks,
                      Slice perCollisions,
                      Assoc association,
                      RevIndices reverseIndices,
                      const int& setTrackSelections,
                      const bool fillTableOfCollIdsPerTrack)
{
  // we do it for all tracks, to be compatible with Run 2 analyses
  for (const auto& collision : collisions) {
    auto tracksThisCollision = tracks.sliceBy(perCollisions, collision.globalIndex());
    for (const auto& track : tracksThisCollision) {
      if constexpr (isCentralBarrel) {
        bool hasGoodQuality = true;
        switch (setTrackSelections) {
          case track_association::CentralBarrelRun2: {
            unsigned char itsClusterMap = track.itsClusterMap();
            if (!(track.tpcNClsFound() >= 50 && track.flags() & o2::aod::track::ITSrefit && track.flags() & o2::aod::track::TPCrefit && (TESTBIT(itsClusterMap, 0) || TESTBIT(itsClusterMap, 1)))) {
              hasGoodQuality = false;
            }
            break;
          }
          case track_association::None: {
            break;
          }
          case track_association::GlobalTrackWoDCA: {
            if (!track.isGlobalTrackWoDCA()) {
              hasGoodQuality = false;
            }
            break;
          }
          case track_association::QualityTracksITS: {
            if (!track.isQualityTrackITS()) {
              hasGoodQuality = false;
            }
            break;
          }
        }
        if (hasGoodQuality) {
          association(collision.globalIndex(), track.globalIndex());
        }
      } else {
        association(collision.globalIndex(), track.globalIndex());
      }
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

template <bool isCentralBarrel, typename TTracksUnfiltered, typename TTracks, typename TAmbiTracks, typename Assoc, typename RevIndices>
void runAssocWithTime(Collisions const& collisions,
                      TTracksUnfiltered const& tracksUnfiltered,
                      TTracks const& tracks,
                      TAmbiTracks const& ambiguousTracks,
                      BCs const& bcs,
                      Assoc association,
                      RevIndices reverseIndices,
                      const float& nSigmaForTimeCompat,
                      const float& timeMargin,
                      const bool usePVAssociation,
                      const bool includeUnassigned,
                      const bool fillTableOfCollIdsPerTrack)
{
  // cache globalBC
  std::vector<uint64_t> globalBC;
  for (const auto& track : tracks) {
    if (track.has_collision()) {
      globalBC.push_back(track.collision().bc().globalBC());
    } else {
      for (const auto& ambTrack : ambiguousTracks) {
        if constexpr (isCentralBarrel) {
          if (ambTrack.trackId() == track.globalIndex()) {
            globalBC.push_back(ambTrack.bc().begin().globalBC());
            break;
          }
        } else {
          if (ambTrack.template getId<TTracks>() == track.globalIndex()) {
            globalBC.push_back(ambTrack.bc().begin().globalBC());
            break;
          }
        }
      }
    }
  }

  // define vector of vectors to store indices of compatible collisions per track
  std::vector<std::unique_ptr<std::vector<int>>> collsPerTrack(tracksUnfiltered.size());

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
      if constexpr (isCentralBarrel) {
        if (usePVAssociation && track.isPVContributor()) {
          thresholdTime = trackTimeRes;
        } else if (TESTBIT(track.flags(), o2::aod::track::TrackTimeResIsRange)) {
          thresholdTime = std::sqrt(sigmaTimeRes2) + timeMargin;
        } else {
          thresholdTime = nSigmaForTimeCompat * std::sqrt(sigmaTimeRes2) + timeMargin;
        }
      } else {
        thresholdTime = nSigmaForTimeCompat * std::sqrt(sigmaTimeRes2) + timeMargin;
      }

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
    for (const auto& track : tracksUnfiltered) {

      const auto trackId = track.globalIndex();
      if (collsPerTrack[trackId] == nullptr) {
        reverseIndices(empty);
      } else {
        reverseIndices(*collsPerTrack[trackId].get());
      }
    }
  }
}

} // namespace o2::aod

#endif // COMMON_DATAMODEL_COLLISIONASSOCIATION_H_
