// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file MIDEfficiency.cxx
/// \brief task to prepare the tables and convert the AO2D forward track into MID track to feed to the efficiency calculator
///
/// @param midEfficiency MID efficiency converter
/// Struct for writing the table and convert the data
/// to MID tracks needed to compute the efficiency of the MID RPCs
///
/// \author Luca Quaglia <luca.quaglia@cern.ch>
///

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/CallbackService.h"

// O2 physics classes
#include "PWGDQ/DataModel/ReducedInfoTables.h"

// O2
#include "DataFormatsMID/Track.h" //MID track from O2
#include "Framework/Variant.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/CompletionPolicyHelpers.h"

#include "MIDEfficiency/Efficiency.h"
#include "MIDBase/DetectorParameters.h"
#include "MIDBase/Mapping.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;

struct midEfficiency {

  // Histogram registry
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurables for histogram axes
  Configurable<int> nBinsLocal{"nBinsLocal", 936, "N bins in local board counts histo"};
  Configurable<int> nBinsRPC{"nBinsRPC", 72, "N bins in RPC counts histo"};
  Configurable<int> nBinsPlane{"nBinsPlane", 4, "N bins in plane counts histo"};
  Configurable<bool> createRootFile{"createRootFile", false, "if true it creates the mid-reco.root file for debug purposes"};

  // Vector of MID tracks to pass to the efficiency calculator
  std::vector<o2::mid::Track> dummyTrack;
  // MID track placeholder for processing
  o2::mid::Track trk;
  // MID mapping for LB calculation
  o2::mid::Mapping mapping;

  void init(o2::framework::InitContext const&)
  {

    LOGF(debug, "Initialization starting");

    // Axes definition
    const AxisSpec axisLocalBoards{nBinsLocal, 0.5, 936.5, "Local board"};
    const AxisSpec axisRPCs{nBinsRPC, -0.5, 71.5, "RPC"};
    const AxisSpec axisPlanes{nBinsPlane, -0.5, 3.5, "Plane"};

    LOGF(debug, "Creating histograms");

    // Same names as O2 task
    // Local boards
    histos.add("nFiredBPperBoard", "nFiredBPperBoard", kTH1F, {axisLocalBoards});
    histos.add("nFiredNBPperBoard", "nFiredNBPperBoard", kTH1F, {axisLocalBoards});
    histos.add("nFiredBothperBoard", "nFiredBothperBoard", kTH1F, {axisLocalBoards});
    histos.add("nTotperBoard", "nTotperBoard", kTH1F, {axisLocalBoards});
    // RPCs
    histos.add("nFiredBPperRPC", "nFiredBPperRPC", kTH1F, {axisRPCs});
    histos.add("nFiredNBPperRPC", "nFiredNBPperRPC", kTH1F, {axisRPCs});
    histos.add("nFiredBothperRPC", "nFiredBothperRPC", kTH1F, {axisRPCs});
    histos.add("nTotperRPC", "nTotperRPC", kTH1F, {axisRPCs});
    // Planes
    histos.add("nFiredBPperPlane", "nFiredBPperPlane", kTH1F, {axisPlanes});
    histos.add("nFiredNBPperPlane", "nFiredNBPperPlane", kTH1F, {axisPlanes});
    histos.add("nFiredBothperPlane", "nFiredBothperPlane", kTH1F, {axisPlanes});
    histos.add("nTotperPlane", "nTotperPlane", kTH1F, {axisPlanes});

  } // end of init

  template <typename TEvent, typename Muons>
  void runMidEffCounters(TEvent const& /*event*/, Muons const& muons)
  {
    LOGF(debug, "Calling process function");

    // Loop over all forward tracks
    // LOGP(info, "collision index = {} ,  nTracks = {}", event.globalIndex(), muons.size());
    for (auto& track : muons) {

      LOGF(debug, "Processing a track");

      trk.setEfficiencyWord(track.midBoards());

      auto deIdMT11 = trk.getFiredDEId();
      auto isRight = o2::mid::detparams::isRightSide(deIdMT11);
      auto rpcLine = o2::mid::detparams::getRPCLine(deIdMT11);
      auto effFlag = trk.getEfficiencyFlag();

      if (effFlag < 0) {
        continue;
      }

      // Loop on the four planes and fill histograms accordingly
      for (int ich = 0; ich < 4; ++ich) {

        bool isFiredBP = trk.isFiredChamber(ich, 0);
        bool isFiredNBP = trk.isFiredChamber(ich, 1);
        // Plane
        // Fill all counts - plane
        histos.fill(HIST("nTotperPlane"), ich);
        if (isFiredBP)
          histos.fill(HIST("nFiredBPperPlane"), ich); // BP - Plane
        if (isFiredNBP)
          histos.fill(HIST("nFiredNBPperPlane"), ich); // NBP - Plane
        if (isFiredBP && isFiredNBP)
          histos.fill(HIST("nFiredBothperPlane"), ich); // Both planes - plane

        if (effFlag < 2) {
          continue;
        }

        // Get RPC id
        auto deId = o2::mid::detparams::getDEId(isRight, ich, rpcLine);

        // Fill all counts - RPC
        histos.fill(HIST("nTotperRPC"), deId);
        if (isFiredBP)
          histos.fill(HIST("nFiredBPperRPC"), deId); // BP - RPC
        if (isFiredNBP)
          histos.fill(HIST("nFiredNBPperRPC"), deId); // NBP - RPC
        if (isFiredBP && isFiredNBP)
          histos.fill(HIST("nFiredBothperRPC"), deId); // Both planes - RPC

        if (effFlag < 3) {
          continue;
        }

        // Get LB ID
        auto firedColumn = trk.getFiredColumnId(); // Get fired column - needed for LB calculation
        auto firedLine = trk.getFiredLineId();     // Get fired line - needed for LB calculation

        auto LB = ich * o2::mid::detparams::NLocalBoards + mapping.getBoardId(firedLine, firedColumn, deId);

        histos.fill(HIST("nTotperBoard"), LB);

        if (isFiredBP)
          histos.fill(HIST("nFiredBPperBoard"), LB);
        if (isFiredNBP)
          histos.fill(HIST("nFiredNBPperBoard"), LB);
        if (isFiredBP && isFiredNBP)
          histos.fill(HIST("nFiredBothperBoard"), LB);
      }
    }

  } // end of runMidEffCounters

  void processMidEffCounter(aod::ReducedEvents::iterator const& event, MyMuonTracks const& muons)
  {
    runMidEffCounters(event, muons); // call efficiency calculator function
  }

  PROCESS_SWITCH(midEfficiency, processMidEffCounter, "process reconstructed information", true);
}; // End of struct midEfficiency

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<midEfficiency>(cfgc)};
}