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

  void processFwdAssocWithTime(Collisions const& collisions,
                               FwdTracks const& muons,
                               AmbiguousFwdTracks const& ambiTracksFwd,
                               BCs const& bcs)
  {
    runAssocWithTime<false>(collisions, muons, muons, ambiTracksFwd, bcs, fwdassociation, fwdreverseIndices, nSigmaForTimeCompat, timeMargin, 0, includeUnassigned, fillTableOfCollIdsPerTrack);
  }
  PROCESS_SWITCH(FwdTrackToCollisionAssociation, processFwdAssocWithTime, "Use fwdtrack-to-collision association based on time", true);

  void processFwdStandardAssoc(Collisions const& collisions,
                               FwdTracks const& muons)
  {
    runStandardAssoc<false>(collisions, muons, muonsPerCollisions, fwdassociation, fwdreverseIndices, 0, fillTableOfCollIdsPerTrack);
  }
  PROCESS_SWITCH(FwdTrackToCollisionAssociation, processFwdStandardAssoc, "Use standard fwdtrack-to-collision association", false);

  void processMFTAssocWithTime(Collisions const& collisions,
                               MFTTracks const& tracks,
                               AmbiguousMFTTracks const& ambiguousTracks,
                               BCs const& bcs)
  {
    runAssocWithTime<false>(collisions, tracks, tracks, ambiguousTracks, bcs, mftassociation, mftreverseIndices, nSigmaForTimeCompat, timeMargin, 0, includeUnassigned, fillTableOfCollIdsPerTrack);
  }
  PROCESS_SWITCH(FwdTrackToCollisionAssociation, processMFTAssocWithTime, "Use MFTtrack-to-collision association based on time", true);

  void processMFTStandardAssoc(Collisions const& collisions,
                               MFTTracks const& tracks)
  {
    runStandardAssoc<false>(collisions, tracks, mftsPerCollisions, mftassociation, mftreverseIndices, 0, fillTableOfCollIdsPerTrack);
  }
  PROCESS_SWITCH(FwdTrackToCollisionAssociation, processMFTStandardAssoc, "Use standard mfttrack-to-collision association", false);
};

//________________________________________________________________________________________________________________________
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<FwdTrackToCollisionAssociation>(cfgc)};
}
