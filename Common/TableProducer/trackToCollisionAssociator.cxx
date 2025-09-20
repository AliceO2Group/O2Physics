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

#include "Common/Core/CollisionAssociation.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

struct TrackToCollisionAssociation {

  Produces<TrackAssoc> association;
  Produces<TrackCompColls> reverseIndices;

  Configurable<float> nSigmaForTimeCompat{"nSigmaForTimeCompat", 4.f, "number of sigmas for time compatibility"};
  Configurable<float> timeMargin{"timeMargin", 0.f, "time margin in ns added to uncertainty because of uncalibrated TPC"};
  Configurable<int> setTrackSelections{"setTrackSelections", 1, "flag to apply track selections: -1=minimal track selection for Run 2 (standard association); 0=none; 1=global track w/o DCA selection; 2=only ITS quality"};
  Configurable<bool> usePVAssociation{"usePVAssociation", true, "if the track is a PV contributor, use the collision time for it"};
  Configurable<bool> includeUnassigned{"includeUnassigned", false, "consider also tracks which are not assigned to any collision"};
  Configurable<bool> fillTableOfCollIdsPerTrack{"fillTableOfCollIdsPerTrack", false, "fill additional table with vector of collision ids per track"};
  Configurable<int> bcWindowForOneSigma{"bcWindowForOneSigma", 60, "BC window to be multiplied by the number of sigmas to define maximum window to be considered"};

  CollisionAssociation<true> collisionAssociator;

  Filter trackFilter = (setTrackSelections.node() == 0) ||                                        // no track selections
                       ((setTrackSelections.node() == 1) && requireGlobalTrackWoDCAInFilter()) || // global track selections w/o dca
                       ((setTrackSelections.node() == 2) && requireQualityTracksITSInFilter());   // ITS-quality selections only
  using TracksWithSel = soa::Join<Tracks, TracksExtra, TrackSelection>;
  using TracksWithSelFilter = soa::Filtered<TracksWithSel>;
  Preslice<TracksWithSel> tracksPerCollisions = aod::track::collisionId;

  void init(InitContext const&)
  {
    if (doprocessAssocWithTime == doprocessStandardAssoc) {
      LOGP(fatal, "Exactly one process function between standard and time-based association should be enabled!");
    }

    // set options in track-to-collision association
    collisionAssociator.setNumSigmaForTimeCompat(nSigmaForTimeCompat);
    collisionAssociator.setTimeMargin(timeMargin);
    collisionAssociator.setTrackSelectionOptionForStdAssoc(setTrackSelections);
    collisionAssociator.setUsePvAssociation(usePVAssociation);
    collisionAssociator.setIncludeUnassigned(includeUnassigned);
    collisionAssociator.setFillTableOfCollIdsPerTrack(fillTableOfCollIdsPerTrack);
    collisionAssociator.setBcWindow(bcWindowForOneSigma);
  }

  void processAssocWithTime(Collisions const& collisions, TracksWithSel const& tracksUnfiltered, TracksWithSelFilter const& tracks, AmbiguousTracks const& ambiguousTracks, BCs const& bcs)
  {
    collisionAssociator.runAssocWithTime(collisions, tracksUnfiltered, tracks, ambiguousTracks, bcs, association, reverseIndices);
  }
  PROCESS_SWITCH(TrackToCollisionAssociation, processAssocWithTime, "Use track-to-collision association based on time", true);

  void processStandardAssoc(Collisions const& collisions, TracksWithSel const& tracksUnfiltered)
  {
    collisionAssociator.runStandardAssoc(collisions, tracksUnfiltered, tracksPerCollisions, association, reverseIndices);
  }
  PROCESS_SWITCH(TrackToCollisionAssociation, processStandardAssoc, "Use standard track-to-collision association", false);
};

//________________________________________________________________________________________________________________________
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TrackToCollisionAssociation>(cfgc)};
}
