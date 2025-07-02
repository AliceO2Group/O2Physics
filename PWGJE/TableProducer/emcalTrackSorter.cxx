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
///
/// EMCAL track sorting and selection task
///
/// \file emcalTrackSorter.cxx
/// \brief Task that selects tracks and sorts them to BCid to be later used in track matching
/// \author Marvin Hemmer (marvin.hemmer@cern.ch) Goethe-University
///

#include "PWGJE/DataModel/EMCALClusters.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <cstdint>
#include <span>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using MyGlobTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TrackSelection, o2::aod::TracksDCA, o2::aod::TracksCov>;
using CollEventSels = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
using FilteredTracks = o2::soa::Filtered<MyGlobTracks>;
using FilteredCollisions = o2::soa::Filtered<CollEventSels>;

struct EmcalTrackSorter {
  Produces<o2::aod::SortedTracks> sortedTracks;

  // Options for the track selection
  Configurable<float> trackMinPt{"trackMinPt", 0.3f, "Minimum pT for tracks to be selected."};
  Configurable<float> trackDCAz{"trackDCAz", 1.f, "Maximum DCAz of a track to be selected."};
  Configurable<float> trackDCAxy{"trackDCAxy", 1.f, "Maximum DCAxy of a track to be selected."};
  Configurable<bool> useOnlyGlobal{"useOnlyGlobal", true, "Select only global tracks."};

  // Filter for the tracks
  const float trackNotOnEMCal = -900;
  Filter trackFilter = (aod::track::pt >= trackMinPt) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t)true)) && (aod::track::trackEtaEmcal > trackNotOnEMCal) && (aod::track::trackPhiEmcal > trackNotOnEMCal) && (nabs(aod::track::dcaXY) < trackDCAxy) && (nabs(aod::track::dcaZ) < trackDCAz) && (aod::track::collisionId >= 0);
  Filter collisionFilter = aod::evsel::foundBCId != -1;

  // QA
  o2::framework::HistogramRegistry mHistManager{"EmcalTrackSorterQAHistograms"};

  // Preslice
  Preslice<o2::aod::Tracks> perCollision = o2::aod::track::collisionId;

  struct MySortedTracks {
    int bcID;
    int trackID;
    float eta;
    float phi;
  };

  void init(InitContext const&)
  {
    LOG(debug) << "Start init!";

    // Setup QA hists.
    // using O2HistType = o2::framework::HistType;
    // o2::framework::AxisSpec etaAxis{160, -0.8, 0.8, "#it{#eta}"};
    // o2::framework::AxisSpec phiAxis{72, 0, 2 * 3.14159, "#it{#varphi} (rad)"};

    // mHistManager.add("hTrackEtaPhiEMCal", "hTrackEtaPhiEMCal", O2HistType::kTH2D, {etaAxis, phiAxis});
    LOG(debug) << "Completed init!";
  }

  //  Appears to need the BC to be accessed to be available in the collision table...
  void processFull(FilteredCollisions const& cols, FilteredTracks const& tracks)
  {
    LOG(debug) << "Starting process full.";

    // outer vector is size as number of ALL calocells in current DF, inside vector will store trackIDs
    std::vector<MySortedTracks> vecSortedTracks;
    vecSortedTracks.reserve(tracks.size());

    for (const auto& col : cols) {
      if (!col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        continue;
      }
      auto groupedTracks = tracks.sliceBy(perCollision, col.globalIndex());
      auto foundBC = col.foundBCId();
      // loop over all selected tracks
      for (const auto& track : groupedTracks) {
        // fill the bcID, trackID, eta and phi values at EMCal surface to the vector
        vecSortedTracks.emplace_back(foundBC, track.globalIndex(), track.trackEtaEmcal(), RecoDecay::constrainAngle(track.trackPhiEmcal()));
      } // end of loop over tracks
    } // end of loop over collisions
    // sort the tracks of this DF according to their bcID
    std::sort(vecSortedTracks.begin(), vecSortedTracks.end(), [](const MySortedTracks& a, const MySortedTracks& b) {
      return a.bcID < b.bcID;
    });
    for (const auto& vecSortedTracksIter : vecSortedTracks) {
      // mHistManager.fill(HIST("hTrackEtaPhiEMCal"), vecSortedTracksIter.eta, vecSortedTracksIter.phi);
      sortedTracks(vecSortedTracksIter.bcID, vecSortedTracksIter.trackID, vecSortedTracksIter.eta, vecSortedTracksIter.phi);
    }
    LOG(debug) << "Ending process full.";
  }
  PROCESS_SWITCH(EmcalTrackSorter, processFull, "run full track sorting analysis", true);

  //  Appears to need the BC to be accessed to be available in the collision table...
  void processNormalSort(CollEventSels const&, FilteredTracks const& tracks)
  {
    LOG(debug) << "Starting process full.";

    std::vector<MySortedTracks> vecSortedTracks;
    vecSortedTracks.reserve(tracks.size());

    int64_t prevColId = -1;
    int32_t cachedBCID = -2;
    bool cachedCollisionOk = true;

    for (const auto& track : tracks) {
      const auto col = track.collision_as<CollEventSels>();
      const int64_t colID = col.globalIndex();

      if (colID != prevColId) {
        cachedBCID = col.foundBCId();
        cachedCollisionOk = col.selection_bit(o2::aod::evsel::kNoSameBunchPileup);
        prevColId = colID;
      }

      if (cachedCollisionOk && cachedBCID >= 0) {
        vecSortedTracks.emplace_back(cachedBCID, track.globalIndex(), track.trackEtaEmcal(), RecoDecay::constrainAngle(track.trackPhiEmcal()));
      }
    }

    std::sort(vecSortedTracks.begin(), vecSortedTracks.end(),
              [](const MySortedTracks& a, const MySortedTracks& b) {
                return a.bcID < b.bcID;
              });

    sortedTracks.reserve(vecSortedTracks.size());
    for (const auto& sTrack : vecSortedTracks) {
      // mHistManager.fill(HIST("hTrackEtaPhiEMCal"), sTrack.eta, sTrack.phi);
      sortedTracks(sTrack.bcID, sTrack.trackID, sTrack.eta, sTrack.phi);
    }

    LOG(debug) << "Ending process full.";
  }
  PROCESS_SWITCH(EmcalTrackSorter, processNormalSort, "run full track sorting analysis by looping over tracks not collisions", false);

  // Appears to need the BC to be accessed to be available in the collision table...
  void processBoostSort(FilteredCollisions const& collisions, FilteredTracks const& tracks)
  {

    LOG(debug) << "Starting process full.";

    if (!tracks.size()) {
      LOG(debug) << "No tracks found.";
      return;
    }

    std::vector<MySortedTracks> vecSortedTracks;
    vecSortedTracks.reserve(tracks.size());

    int32_t cachedBCID = -1;

    auto track = tracks.begin();
    const auto trackEnd = tracks.end();
    int32_t currentCollId = track.collisionId();

    auto selection = collisions.getSelectedRows();
    for (auto& collisionId : selection) {
      // tracks are behind collision, move ahead in track table
      while (track != trackEnd && collisionId > currentCollId) {
        if (++track != trackEnd) {
          currentCollId = track.collisionId();
        }
      }
      // tracks are ahead of collisions, move to next collision
      if (collisionId < currentCollId) {
        continue;
      }

      // we reached the first track for this collision, get the BCId
      cachedBCID = track.collision_as<FilteredCollisions>().foundBCId();
      // loop over all tracks in current collision
      while (track != trackEnd && collisionId == currentCollId) {
        vecSortedTracks.emplace_back(cachedBCID, track.globalIndex(), track.trackEtaEmcal(), RecoDecay::constrainAngle(track.trackPhiEmcal()));
        if (++track != trackEnd) {
          currentCollId = track.collisionId();
        }
      }
    }

    std::sort(vecSortedTracks.begin(), vecSortedTracks.end(),
              [](const MySortedTracks& a, const MySortedTracks& b) {
                return a.bcID < b.bcID;
              });

    sortedTracks.reserve(vecSortedTracks.size());
    for (const auto& sTrack : vecSortedTracks) {
      // mHistManager.fill(HIST("hTrackEtaPhiEMCal"), sTrack.eta, sTrack.phi);
      sortedTracks(sTrack.bcID, sTrack.trackID, sTrack.eta, sTrack.phi);
    }

    LOG(debug) << "Ending process full.";
  }
  PROCESS_SWITCH(EmcalTrackSorter, processBoostSort, "run full track sorting analysis by looping over tracks not collisions", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<EmcalTrackSorter>(cfgc)};
}
