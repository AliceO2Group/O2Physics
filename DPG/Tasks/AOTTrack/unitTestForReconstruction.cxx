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
/// \file unitTestForReconstruction.cxx
///
/// \brief Unit test for validating the reconstruction software
/// \author Alberto Caliva (alberto.caliva@cern.ch), Catalin-Lucian Ristea (catalin.ristea@cern.ch)
/// \since September 9, 2025

#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataTypes.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"

#include <cmath>
#include <memory>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::constants::math;

struct UnitTestForReconstruction {

  // Histogram registry
  HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // global IDs of events to be inspected
  Configurable<std::vector<int>> eventNr{"eventNr", {1, 5, 12, 44, 76, 99, 102, 115, 180, 220}, "eventNr"};
  std::unordered_set<int> eventSet;

  void init(InitContext const&)
  {
    // Define histogram to monitor event counts at different selection stages
    registryData.add("eventCounter", "eventCounter", HistType::kTH1F, {{10, 0, 10, ""}});

    // Define histogram for the transverse momentum spectrum of reconstructed charged tracks
    registryData.add("ptChargedTracks", "ptChargedTracks", HistType::kTH2F, {{11, 0, 11, "event"}, {1000, 0, 10, "#it{p}_{T} (GeV/#it{c})"}});

    // Fast lookup set from configurable event list
    eventSet = std::unordered_set<int>(eventNr->begin(), eventNr->end());
  }

  // Process Data
  void processData(o2::aod::Collisions const& collisions, o2::aod::Tracks const& tracks)
  {
    // Event index
    int eventIndex = 0;
    static constexpr int indexAllEvts = 0;

    // Loop over reconstructed events
    for (const auto& collision : collisions) {

      // Event counter: before event selection
      registryData.fill(HIST("eventCounter"), 0.5);

      // Check if event global index is in the list of events to process
      int ev = collision.globalIndex();
      if (eventSet.count(ev)) {

        // Increment event index
        eventIndex++;

        // Fill event counter
        registryData.fill(HIST("eventCounter"), 1.5);

        // Loop over reconstructed tracks
        for (auto const& track : tracks) {
          registryData.fill(HIST("ptChargedTracks"), indexAllEvts, track.pt());
          registryData.fill(HIST("ptChargedTracks"), eventIndex, track.pt());
        }
      }
    }
  }
  PROCESS_SWITCH(UnitTestForReconstruction, processData, "Process Data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<UnitTestForReconstruction>(cfgc)};
}
