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

/// \file alice3TrackingTranslator.cxx
///
/// \brief Translator task to convert tracking software to the AO2D format digestible with the O2Physics analysis framework
///
/// \author Nicolò Jacazio, Universita del Piemonte Orientale (IT)
///

#include "ALICE3/DataModel/collisionAlice3.h"
#include "ALICE3/DataModel/tracksAlice3.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

struct Alice3TrackingTranslator {
  o2::framework::Produces<o2::aod::Collisions> tableCollisions;
  o2::framework::Produces<o2::aod::McCollisionLabels> tableMcCollisionLabels;
  o2::framework::Produces<o2::aod::StoredTracks> tableStoredTracks;
  o2::framework::Produces<o2::aod::TracksExtension> tableTracksExtension;
  o2::framework::Produces<o2::aod::StoredTracksCov> tableStoredTracksCov;
  o2::framework::Produces<o2::aod::TracksCovExtension> tableTracksCovExtension;
  o2::framework::Produces<o2::aod::McTrackLabels> tableMcTrackLabels;
  o2::framework::Produces<o2::aod::TracksDCA> tableTracksDCA;
  o2::framework::Produces<o2::aod::TracksDCACov> tableTracksDCACov;
  o2::framework::Produces<o2::aod::CollisionsAlice3> tableCollisionsAlice3;
  o2::framework::Produces<o2::aod::TracksAlice3> tableTracksAlice3;
  o2::framework::Produces<o2::aod::TracksExtraA3> tableTracksExtraA3;

  o2::framework::Produces<o2::aod::StoredTracksExtra_002> tableStoredTracksExtra;
  o2::framework::Produces<o2::aod::TrackSelection> tableTrackSelection;
  o2::framework::Produces<o2::aod::TrackSelectionExtension> tableTrackSelectionExtension;

  void init(o2::framework::InitContext&)
  {
    // Initialization if needed
    LOG(info) << "Alice3TrackingTranslator init called";
  }

  void process(o2::aod::BCs const&)
  {
    LOG(info) << "Alice3TrackingTranslator process called";
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::framework::WorkflowSpec w;
  w.push_back(adaptAnalysisTask<Alice3TrackingTranslator>(cfgc));
  return w;
}
