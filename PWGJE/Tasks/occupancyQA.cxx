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

// QA task for occupancy in the jet framework
//
/// \author Aimeric Landou <aimeric.landou@cern.ch>

#include <cmath>
#include <TRandom3.h>
#include <string>

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "PWGJE/Core/JetDerivedDataUtilities.h"

#include "EventFiltering/filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct OccupancyQATask {

  HistogramRegistry registry;

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 100.0, "maximum pT acceptance for tracks"};
  // Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  // Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  // Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};
  // Configurable<bool> checkMcCollisionIsMatched{"checkMcCollisionIsMatched", false, "0: count whole MCcollisions, 1: select MCcollisions which only have their correspond collisions"};

  int eventSelection = -1;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    AxisSpec occupancyAxis = {142, -1.5, 14000.5, "occupancy"};
    AxisSpec nTracksAxis = {16001, -1., 16000, "n tracks"};

    if (doprocessEventsJetData) {
      registry.add("h2_occupancy_ntracksall_presel8", "occupancy vs N_{tracks}; occupancy; N_{tracks}", {HistType::kTH2I, {occupancyAxis, nTracksAxis}});
      registry.add("h2_occupancy_ntracksall_postsel8", "occupancy vs N_{tracks}; occupancy; N_{tracks}", {HistType::kTH2I, {occupancyAxis, nTracksAxis}});
      registry.add("h2_occupancy_ntrackssel_presel8", "occupancy vs N_{tracks}; occupancy; N_{tracks}", {HistType::kTH2I, {occupancyAxis, nTracksAxis}});
      registry.add("h2_occupancy_ntrackssel_postsel8", "occupancy vs N_{tracks}; occupancy; N_{tracks}", {HistType::kTH2I, {occupancyAxis, nTracksAxis}});
      registry.add("h2_occupancy_ntracksselptetacuts_presel8", "occupancy vs N_{tracks}; occupancy; N_{tracks}", {HistType::kTH2I, {occupancyAxis, nTracksAxis}});
      registry.add("h2_occupancy_ntracksselptetacuts_postsel8", "occupancy vs N_{tracks}; occupancy; N_{tracks}", {HistType::kTH2I, {occupancyAxis, nTracksAxis}});
    }
  }

  // Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);
  Filter eventCutsAOD = (nabs(aod::collision::posZ) < vertexZCut);

  void processEventsJetData(soa::Filtered<aod::JetCollisions>::iterator const& collision, aod::JetTracks const& tracks)
  {

    // Before sel8-like selection
    int occupancy = collision.trackOccupancyInTimeRange();
    int nTracksAll = tracks.size();
    int nTracksAllAcceptanceAndSelected = 0;
    int nTracksInAcceptanceAndSelected = 0;
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        nTracksAllAcceptanceAndSelected += 1;
        if (track.pt() >= trackPtMin && track.pt() < trackPtMax && track.eta() > trackEtaMin && track.eta() < trackEtaMax) {
          nTracksInAcceptanceAndSelected += 1;
        }
      }
    }

    registry.fill(HIST("h2_occupancy_ntracksall_presel8"), occupancy, nTracksAll);
    registry.fill(HIST("h2_occupancy_ntrackssel_presel8"), occupancy, nTracksAllAcceptanceAndSelected);
    registry.fill(HIST("h2_occupancy_ntracksselptetacuts_presel8"), occupancy, nTracksInAcceptanceAndSelected);
    if (jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      registry.fill(HIST("h2_occupancy_ntracksall_postsel8"), occupancy, nTracksAll);
      registry.fill(HIST("h2_occupancy_ntrackssel_postsel8"), occupancy, nTracksAllAcceptanceAndSelected);
      registry.fill(HIST("h2_occupancy_ntracksselptetacuts_postsel8"), occupancy, nTracksInAcceptanceAndSelected);
    }
  }
  PROCESS_SWITCH(OccupancyQATask, processEventsJetData, "occupancy QA on jet derived data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<OccupancyQATask>(cfgc, TaskName{"occupancy-qa"})}; }