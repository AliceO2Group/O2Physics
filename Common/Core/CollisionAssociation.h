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
/// \brief Class for track to collision associators
/// \author Jan Fiete Grosse-Oetringhaus <jan.fiete.grosse-oetringhaus@cern.ch>, CERN
/// \author Fabrizio Grosa <fgrosa@cern.ch>, CERN
/// \author Mattia Faggin <mfaggin@cern.ch>, University and INFN Padova
/// \author Sarah Herrmann <sarah.herrmann@cern.ch>, IP2I Lyon
/// \author Maurice Coquet <maurice.louis.coquet@cern.ch>, CEA-Saclay/Irfu

#ifndef COMMON_CORE_COLLISIONASSOCIATION_H_
#define COMMON_CORE_COLLISIONASSOCIATION_H_

#include <vector>
#include <memory>

#include "CommonConstants/LHCConstants.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

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

} // namespace track_association
} // namespace o2::aod

template <bool isCentralBarrel>
class CollisionAssociation
{
 public:
  /// Default constructor
  CollisionAssociation() = default;

  void setNumSigmaForTimeCompat(float nSigma) { mNumSigmaForTimeCompat = nSigma; }
  void setTimeMargin(float timeMargin) { mTimeMargin = timeMargin; }
  void setTrackSelectionOptionForStdAssoc(int option) { mTrackSelection = option; }
  void setUsePvAssociation(bool enable = true) { mUsePvAssociation = enable; }
  void setIncludeUnassigned(bool enable = true) { mIncludeUnassigned = enable; }
  void setFillTableOfCollIdsPerTrack(bool fill = true) { mFillTableOfCollIdsPerTrack = fill; }
  void setBcWindow(int bcWindow = 115) { mBcWindowForOneSigma = bcWindow; }

  template <typename TTracks, typename Slice, typename Assoc, typename RevIndices>
  void runStandardAssoc(o2::aod::Collisions const& collisions,
                        TTracks const& tracks,
                        Slice& perCollisions,
                        Assoc& association,
                        RevIndices& reverseIndices)
  {
    // we do it for all tracks, to be compatible with Run 2 analyses
    for (const auto& collision : collisions) {
      auto tracksThisCollision = tracks.sliceBy(perCollisions, collision.globalIndex());
      for (const auto& track : tracksThisCollision) {
        if constexpr (isCentralBarrel) {
          bool hasGoodQuality = true;
          switch (mTrackSelection) {
            case o2::aod::track_association::TrackSelection::CentralBarrelRun2: {
              unsigned char itsClusterMap = track.itsClusterMap();
              if (!(track.tpcNClsFound() >= 50 && track.flags() & o2::aod::track::ITSrefit && track.flags() & o2::aod::track::TPCrefit && (TESTBIT(itsClusterMap, 0) || TESTBIT(itsClusterMap, 1)))) {
                hasGoodQuality = false;
              }
              break;
            }
            case o2::aod::track_association::TrackSelection::None: {
              break;
            }
            case o2::aod::track_association::TrackSelection::GlobalTrackWoDCA: {
              if (!track.isGlobalTrackWoDCA()) {
                hasGoodQuality = false;
              }
              break;
            }
            case o2::aod::track_association::TrackSelection::QualityTracksITS: {
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
    if (mFillTableOfCollIdsPerTrack) {
      for (const auto& track : tracks) {
        if (track.has_collision()) {
          reverseIndices(std::vector<int>{track.collisionId()});
        } else {
          reverseIndices(empty);
        }
      }
    }
  }

  template <typename TTracksUnfiltered, typename TTracks, typename TAmbiTracks, typename Assoc, typename RevIndices>
  void runAssocWithTime(o2::aod::Collisions const& collisions,
                        TTracksUnfiltered const& tracksUnfiltered,
                        TTracks const& tracks,
                        TAmbiTracks const& ambiguousTracks,
                        o2::aod::BCs const& bcs,
                        Assoc& association,
                        RevIndices& reverseIndices)
  {
    // cache globalBC
    std::vector<uint64_t> globalBC;
    for (const auto& track : tracks) {
      if (track.has_collision()) {
        globalBC.push_back(track.collision().bc().globalBC());
      } else {
        for (const auto& ambTrack : ambiguousTracks) {
          if constexpr (isCentralBarrel) { // FIXME: to be removed as soon as it is possible to use getId<Table>() for joined tables
            if (ambTrack.trackId() == track.globalIndex()) {
              if (!ambTrack.has_bc() || ambTrack.bc().size() == 0) {
                globalBC.push_back(-1);
                break;
              }
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
    auto bOffsetMax = mBcWindowForOneSigma * mNumSigmaForTimeCompat + mTimeMargin / o2::constants::lhc::LHCBunchSpacingNS;
    for (const auto& collision : collisions) {
      const float collTime = collision.collisionTime();
      const float collTimeRes2 = collision.collisionTimeRes() * collision.collisionTimeRes();
      uint64_t collBC = collision.bc().globalBC();
      // bool iteratorMoved = false;
      for (auto track = trackBegin; track != tracks.end(); ++track) {
        if (!mIncludeUnassigned && !track.has_collision()) {
          continue;
        }

        float trackTime = track.trackTime();
        if (globalBC[track.filteredIndex()] < 0) {
          continue;
        }
        const int64_t bcOffsetWindow = (int64_t)globalBC[track.filteredIndex()] + trackTime / o2::constants::lhc::LHCBunchSpacingNS - (int64_t)collBC;
        if (std::abs(bcOffsetWindow) > bOffsetMax) {
          continue;
        }

        float trackTimeRes = track.trackTimeRes();
        if constexpr (isCentralBarrel) {
          if (mUsePvAssociation && track.isPVContributor()) {
            trackTime = track.collision().collisionTime();    // if PV contributor, we assume the time to be the one of the collision
            trackTimeRes = o2::constants::lhc::LHCBunchSpacingNS; // 1 BC
          }
        }

        const int64_t bcOffset = (int64_t)globalBC[track.filteredIndex()] - (int64_t)collBC;
        const float deltaTime = trackTime - collTime + bcOffset * o2::constants::lhc::LHCBunchSpacingNS;
        float sigmaTimeRes2 = collTimeRes2 + trackTimeRes * trackTimeRes;
        LOGP(debug, "collision time={}, collision time res={}, track time={}, track time res={}, bc collision={}, bc track={}, delta time={}", collTime, collision.collisionTimeRes(), track.trackTime(), track.trackTimeRes(), collBC, globalBC[track.filteredIndex()], deltaTime);

        // optimization to avoid looping over the full track list each time. This assumes that tracks are sorted by BCs (which they should be because collisions are sorted by BCs)
        // NOTE this does not work anymore if mIncludeUnassigned is set as the unassigned blocks can be somewhere (and we can have merged DFs, too)
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
          if (mUsePvAssociation && track.isPVContributor()) {
            thresholdTime = trackTimeRes;
          } else if (TESTBIT(track.flags(), o2::aod::track::TrackTimeResIsRange)) {
            thresholdTime = std::sqrt(sigmaTimeRes2) + mTimeMargin;
          } else {
            thresholdTime = mNumSigmaForTimeCompat * std::sqrt(sigmaTimeRes2) + mTimeMargin;
          }
        } else {
          thresholdTime = mNumSigmaForTimeCompat * std::sqrt(sigmaTimeRes2) + mTimeMargin;
        }

        if (std::abs(deltaTime) < thresholdTime) {
          const auto collIdx = collision.globalIndex();
          const auto trackIdx = track.globalIndex();
          LOGP(debug, "Filling track id {} for coll id {}", trackIdx, collIdx);
          association(collIdx, trackIdx);
          if (mFillTableOfCollIdsPerTrack) {
            if (collsPerTrack[trackIdx] == nullptr) {
              collsPerTrack[trackIdx] = std::make_unique<std::vector<int>>();
            }
            collsPerTrack[trackIdx].get()->push_back(collIdx);
          }
        }
      }
    }
    // create reverse index track to collisions if enabled
    if (mFillTableOfCollIdsPerTrack) {
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

 private:
  float mNumSigmaForTimeCompat{4.};                                                  // number of sigma for time compatibility
  float mTimeMargin{500.};                                                           // additional time margin in ns
  int mTrackSelection{o2::aod::track_association::TrackSelection::GlobalTrackWoDCA}; // track selection for central barrel tracks (standard association only)
  bool mUsePvAssociation{true};                                                      // use the information of PV contributors
  bool mIncludeUnassigned{true};                                                     // include tracks that were originally not assigned to any collision
  bool mFillTableOfCollIdsPerTrack{false};                                           // fill additional table with vectors of compatible collisions per track
  int mBcWindowForOneSigma{115};                                                     // BC window to be multiplied by the number of sigmas to define maximum window to be considered
};

#endif // COMMON_CORE_COLLISIONASSOCIATION_H_
