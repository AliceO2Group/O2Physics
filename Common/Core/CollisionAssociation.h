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

#include <CommonConstants/LHCConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/DataTypes.h>
#include <Framework/Logger.h>

#include <Rtypes.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <utility>
#include <vector>

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
              int minTpcNClsFound{50};
              if (!(track.tpcNClsFound() >= minTpcNClsFound && track.flags() & o2::aod::track::ITSrefit && track.flags() & o2::aod::track::TPCrefit && (TESTBIT(itsClusterMap, 0) || TESTBIT(itsClusterMap, 1)))) {
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
    // cache globalBC and track time in BC for optimization
    std::vector<int64_t> globalBC;
    std::vector<int64_t> trackBCCache;
    std::vector<std::pair<typename TTracks::iterator, typename TTracks::iterator>> trackIterationWindows; // continous regions in which we can count on increasing globalBC numbers
    globalBC.reserve(tracks.size());
    trackBCCache.reserve(tracks.size());
    auto trackBegin = tracks.begin();
    int lastCollisionId = 0;
    if (tracks.size() > 0) {
      lastCollisionId = trackBegin.collisionId();
    }
    auto track = trackBegin;
    for (; track != tracks.end(); ++track) {
      int64_t trackBC = -1;
      if (track.has_collision()) {
        trackBC = track.collision().bc().globalBC();
      } else if (mIncludeUnassigned) {
        for (const auto& ambTrack : ambiguousTracks) {
          if constexpr (isCentralBarrel) { // FIXME: to be removed as soon as it is possible to use getId<Table>() for joined tables
            if (ambTrack.trackId() == track.globalIndex()) {
              // special check to avoid crashes (in particular on some MC datasets)
              // related to shifts in ambiguous tracks association to bc slices (off by 1) - see https://mattermost.web.cern.ch/alice/pl/g9yaaf3tn3g4pgn7c1yex9copy
              if (ambTrack.bcIds()[0] >= bcs.size() || ambTrack.bcIds()[1] >= bcs.size()) {
                break;
              }
              if (!ambTrack.has_bc() || ambTrack.bc().size() == 0) {
                break;
              }
              trackBC = ambTrack.bc().begin().globalBC();
              break;
            }
          } else {
            if (ambTrack.template getId<TTracks>() == track.globalIndex()) {
              trackBC = ambTrack.bc().begin().globalBC();
              break;
            }
          }
        }
      }
      globalBC.push_back(trackBC);
      trackBCCache.push_back(trackBC + track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS);
      // find uniform blocks
      if ((track.collisionId() < lastCollisionId) || (lastCollisionId < 0 && track.collisionId() >= 0)) {
        if (lastCollisionId >= 0 || mIncludeUnassigned) {
          LOGP(debug, "Found track block from {} to {}, current id {}, last id {}", trackBegin.filteredIndex(), track.filteredIndex() - 1, track.collisionId(), lastCollisionId);
          trackIterationWindows.push_back(std::make_pair(trackBegin, track));
        }
        trackBegin = track;
      }
      lastCollisionId = track.collisionId();
    }
    // trackIterationWindows.push_back(std::make_pair(trackBegin, tracks.end()));
    if (lastCollisionId >= 0 || mIncludeUnassigned) {
      LOGP(debug, "Found track block from {} to {}", trackBegin.filteredIndex(), tracks.size() - 1);
      trackIterationWindows.push_back(std::make_pair(trackBegin, track));
    }

    // define vector of vectors to store indices of compatible collisions per track
    std::vector<std::unique_ptr<std::vector<int>>> collsPerTrack(tracksUnfiltered.size());

    // loop over collisions to find time-compatible tracks
    int64_t bcOffsetMax = mBcWindowForOneSigma * mNumSigmaForTimeCompat + mTimeMargin / o2::constants::lhc::LHCBunchSpacingNS;
    for (const auto& collision : collisions) {
      const float collTime = collision.collisionTime();
      const float collTimeRes2 = collision.collisionTimeRes() * collision.collisionTimeRes();
      uint64_t collBC = collision.bc().globalBC();

      // This is done per block to allow optimization below. Within each block the globalBC increase continously
      for (auto& iterationWindow : trackIterationWindows) { // o2-linter: disable=const-ref-in-for-loop (iterationWindow is modified)
        bool iteratorMoved = false;
        const bool isAssignedTrackWindow = (iterationWindow.first != iterationWindow.second) ? iterationWindow.first.has_collision() : false;
        for (auto trackInWindow = iterationWindow.first; trackInWindow != iterationWindow.second; ++trackInWindow) {
          int64_t trackBC = globalBC[trackInWindow.filteredIndex()];
          if (trackBC < 0) {
            continue;
          }

          // Optimization to avoid looping over the full track list each time. This builds on that tracks are sorted by BCs (which they should be because collisions are sorted by BCs)
          const int64_t bcOffset = trackBC - static_cast<int64_t>(collBC);
          if constexpr (isCentralBarrel) {
            // only for blocks with collision association
            if (isAssignedTrackWindow) {
              constexpr int margin = 200;
              if (!iteratorMoved && bcOffset > -bcOffsetMax - margin) {
                iterationWindow.first.setCursor(trackInWindow.filteredIndex());
                iteratorMoved = true;
                LOGP(debug, "Moving iterator begin {}", trackInWindow.filteredIndex());
              } else if (bcOffset > bcOffsetMax + margin) {
                LOGP(debug, "Stopping iterator {}", trackInWindow.filteredIndex());
                break;
              }
            }
          }

          int64_t bcOffsetWindow = trackBCCache[trackInWindow.filteredIndex()] - static_cast<int64_t>(collBC);
          if (std::abs(bcOffsetWindow) > bcOffsetMax) {
            continue;
          }

          float trackTime = 0;
          float trackTimeRes = 0;
          if constexpr (isCentralBarrel) {
            if (mUsePvAssociation && trackInWindow.isPVContributor()) {
              trackTime = trackInWindow.collision().collisionTime(); // if PV contributor, we assume the time to be the one of the collision
              trackTimeRes = o2::constants::lhc::LHCBunchSpacingNS;  // 1 BC
            } else {
              trackTime = trackInWindow.trackTime();
              trackTimeRes = trackInWindow.trackTimeRes();
            }
          } else {
            trackTime = trackInWindow.trackTime();
            trackTimeRes = trackInWindow.trackTimeRes();
          }

          const float deltaTime = trackTime - collTime + bcOffset * o2::constants::lhc::LHCBunchSpacingNS;
          float sigmaTimeRes2 = collTimeRes2 + trackTimeRes * trackTimeRes;
          LOGP(debug, "collision time={}, collision time res={}, track time={}, track time res={}, bc collision={}, bc track={}, delta time={}", collTime, collision.collisionTimeRes(), trackInWindow.trackTime(), trackInWindow.trackTimeRes(), collBC, trackBC, deltaTime);

          float thresholdTime = 0.;
          if constexpr (isCentralBarrel) {
            if (mUsePvAssociation && trackInWindow.isPVContributor()) {
              thresholdTime = trackTimeRes;
            } else if (TESTBIT(trackInWindow.flags(), o2::aod::track::TrackTimeResIsRange)) {
              // the track time resolution is a range, not a gaussian resolution
              thresholdTime = trackTimeRes + mNumSigmaForTimeCompat * std::sqrt(collTimeRes2) + mTimeMargin;
            } else {
              thresholdTime = mNumSigmaForTimeCompat * std::sqrt(sigmaTimeRes2) + mTimeMargin;
            }
          } else {
            // the track is not a central track
            if constexpr (TTracks::template contains<o2::aod::MFTTracks>()) {
              // then the track is an MFT track, or an MFT track with additionnal joined info
              // in this case TrackTimeResIsRange
              thresholdTime = trackTimeRes + mNumSigmaForTimeCompat * std::sqrt(collTimeRes2) + mTimeMargin;
            } else if constexpr (TTracks::template contains<o2::aod::FwdTracks>()) {
              // the track is a fwd track, with a gaussian time resolution
              thresholdTime = mNumSigmaForTimeCompat * std::sqrt(sigmaTimeRes2) + mTimeMargin;
            }
          }

          if (std::abs(deltaTime) < thresholdTime) {
            const auto collIdx = collision.globalIndex();
            const auto trackIdx = trackInWindow.globalIndex();
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
    }
    // create reverse index track to collisions if enabled
    if (mFillTableOfCollIdsPerTrack) {
      std::vector<int> empty{};
      for (const auto& trackUnfiltered : tracksUnfiltered) {

        const auto trackId = trackUnfiltered.globalIndex();
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
