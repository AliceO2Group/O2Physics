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
/// \author Luca Quaglia <luca.quaglia@cern.ch

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/CallbackService.h"
#include "Framework/ASoAHelpers.h"
#include <iostream>

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
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;

struct midEfficiency {

  // Histogram registry
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurables for histogram axes
  // Centrality
  Configurable<int> nBinsCentrality{"nBinsCentrality", 9, "N bins for centrality histo in PbPb"};
  Configurable<float> minCentrality{"minCentrality", 0., "minimum of axis in centrality histo in PbPb"};
  Configurable<float> maxCentrality{"maxCentrality", 90., "maximum of axis in centrality histo in PbPb"};
  Configurable<bool> isPbPb{"isPbPb", false, "If true, the task will be used to run on PbPb data and will enable the filling of THnSparse for centrality studies"};

  // pt
  Configurable<int> nBinsPt{"nBinsPt", 2000, "N bins for pt histo"};
  Configurable<float> minPt{"minPt", 0., "minimum of pt axis"}; // GeV/c
  Configurable<float> maxPt{"maxPt", 20., "maximum of pt axis"};

  // eta
  Configurable<int> nBinsEta{"nBinsEta", 500, "N bins for eta histo"};
  Configurable<float> minEta{"minEta", -5., "minimum of eta axis"}; //
  Configurable<float> maxEta{"maxEta", 5., "maximum of eta axis"};

  // phi
  Configurable<int> nBinsPhi{"nBinsPhi", 500, "N bins for phi histo"};
  Configurable<float> minPhi{"minPhi", -2. * TMath::Pi(), "minimum of phi axis"}; //
  Configurable<float> maxPhi{"maxPhi", 2 * TMath::Pi(), "maximum of phi axis"};

  // MID track placeholder for processing
  o2::mid::Track trk;
  // MID mapping for LB calculation
  o2::mid::Mapping mapping;

  // Filter only for MCH-MID matched tracks
  Filter muonTrackType = aod::fwdtrack::trackType == uint8_t(3);

  void init(o2::framework::InitContext const& /*ic*/)
  {

    LOGF(debug, "Initialization starting");

    // Axes definition
    const AxisSpec axisLocalBoards{936, 0.5, 936.5, "Local board"}; // These are not defined as configurable since they are fixed
    const AxisSpec axisRPCs{72, -0.5, 71.5, "RPC"};                 // These are not defined as configurable since they are fixed
    const AxisSpec axisPlanes{4, -0.5, 3.5, "Plane"};               // These are not defined as configurable since they are fixed
    const AxisSpec axisTrackType{5, -0.5, 4.5, "Muon track type"};  // These are not defined as configurable since they are fixed

    const AxisSpec axisCent{nBinsCentrality, minCentrality, maxCentrality, "Centrality"};
    const AxisSpec axisPt{nBinsPt, minPt, maxPt, "track p_{t} [GeV/c]"};
    const AxisSpec axisEta{nBinsEta, minEta, maxEta, "track #eta"};
    const AxisSpec axisPhi{nBinsPhi, minPhi, maxPhi, "track #phi"};

    LOGF(debug, "Creating histograms");

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
    // Centrality test
    histos.add("hCentr", "hCentr", kTH1F, {axisCent});

    // Track type
    histos.add("hTrackType", "hTrackType", kTH1F, {axisTrackType});

    // If this is true -> PbPb data, add THnSparse with centrality
    if (isPbPb) {
      // Local boards
      histos.add("hSparseCentFiredBPperBoard", "THn for centrality studies", HistType::kTHnSparseF, {axisLocalBoards, axisCent, axisPt, axisEta, axisPhi});
      histos.add("hSparseCentFiredNBPperBoard", "THn for centrality studies", HistType::kTHnSparseF, {axisLocalBoards, axisCent, axisPt, axisEta, axisPhi});
      histos.add("hSparseCentFiredBothperBoard", "THn for centrality studies", HistType::kTHnSparseF, {axisLocalBoards, axisCent, axisPt, axisEta, axisPhi});
      histos.add("hSparseCentFiredTotperBoard", "THn for centrality studies", HistType::kTHnSparseF, {axisLocalBoards, axisCent, axisPt, axisEta, axisPhi});

      // RPCs
      histos.add("hSparseCentFiredBPperRPC", "THn for centrality studies", HistType::kTHnSparseF, {axisRPCs, axisCent, axisPt, axisEta, axisPhi});
      histos.add("hSparseCentFiredNBPperRPC", "THn for centrality studies", HistType::kTHnSparseF, {axisRPCs, axisCent, axisPt, axisEta, axisPhi});
      histos.add("hSparseCentFiredBothperRPC", "THn for centrality studies", HistType::kTHnSparseF, {axisRPCs, axisCent, axisPt, axisEta, axisPhi});
      histos.add("hSparseCentFiredTotperRPC", "THn for centrality studies", HistType::kTHnSparseF, {axisRPCs, axisCent, axisPt, axisEta, axisPhi});

      // Planes
      histos.add("hSparseCentFiredBPperPlane", "THn for centrality studies", HistType::kTHnSparseF, {axisPlanes, axisCent, axisPt, axisEta, axisPhi});
      histos.add("hSparseCentFiredNBPperPlane", "THn for centrality studies", HistType::kTHnSparseF, {axisPlanes, axisCent, axisPt, axisEta, axisPhi});
      histos.add("hSparseCentFiredBothperPlane", "THn for centrality studies", HistType::kTHnSparseF, {axisPlanes, axisCent, axisPt, axisEta, axisPhi});
      histos.add("hSparseCentFiredTotperPlane", "THn for centrality studies", HistType::kTHnSparseF, {axisPlanes, axisCent, axisPt, axisEta, axisPhi});
    } else { // THnSparse without centrality in pp
      // Local boards
      histos.add("hSparseCentFiredBPperBoard", "THn for centrality studies", HistType::kTHnSparseF, {axisLocalBoards, axisPt, axisEta, axisPhi});
      histos.add("hSparseCentFiredNBPperBoard", "THn for centrality studies", HistType::kTHnSparseF, {axisLocalBoards, axisPt, axisEta, axisPhi});
      histos.add("hSparseCentFiredBothperBoard", "THn for centrality studies", HistType::kTHnSparseF, {axisLocalBoards, axisPt, axisEta, axisPhi});
      histos.add("hSparseCentFiredTotperBoard", "THn for centrality studies", HistType::kTHnSparseF, {axisLocalBoards, axisPt, axisEta, axisPhi});

      // RPCs
      histos.add("hSparseCentFiredBPperRPC", "THn for centrality studies", HistType::kTHnSparseF, {axisRPCs, axisPt, axisEta, axisPhi});
      histos.add("hSparseCentFiredNBPperRPC", "THn for centrality studies", HistType::kTHnSparseF, {axisRPCs, axisPt, axisEta, axisPhi});
      histos.add("hSparseCentFiredBothperRPC", "THn for centrality studies", HistType::kTHnSparseF, {axisRPCs, axisPt, axisEta, axisPhi});
      histos.add("hSparseCentFiredTotperRPC", "THn for centrality studies", HistType::kTHnSparseF, {axisRPCs, axisPt, axisEta, axisPhi});

      // Planes
      histos.add("hSparseCentFiredBPperPlane", "THn for centrality studies", HistType::kTHnSparseF, {axisPlanes, axisPt, axisEta, axisPhi});
      histos.add("hSparseCentFiredNBPperPlane", "THn for centrality studies", HistType::kTHnSparseF, {axisPlanes, axisPt, axisEta, axisPhi});
      histos.add("hSparseCentFiredBothperPlane", "THn for centrality studies", HistType::kTHnSparseF, {axisPlanes, axisPt, axisEta, axisPhi});
      histos.add("hSparseCentFiredTotperPlane", "THn for centrality studies", HistType::kTHnSparseF, {axisPlanes, axisPt, axisEta, axisPhi});
    }
  } // end of init

  template <typename TEvent, typename Muons>
  void runMidEffCounters(TEvent const& event, Muons const& muons)
  {
    LOGF(debug, "Calling process function");

    float cent = event.centFT0C();

    if (isPbPb)
      histos.fill(HIST("hCentr"), cent); // Fill centrality histo

    // Loop over all forward tracks
    for (auto& track : muons) {

      LOGF(debug, "Processing a track");

      histos.fill(HIST("hTrackType"), track.trackType());

      trk.setEfficiencyWord(track.midBoards());

      auto deIdMT11 = trk.getFiredDEId();
      auto isRight = o2::mid::detparams::isRightSide(deIdMT11);
      auto rpcLine = o2::mid::detparams::getRPCLine(deIdMT11);
      auto effFlag = trk.getEfficiencyFlag();

      float pt = track.pt();
      float eta = track.eta();
      float phi = track.phi();

      if (effFlag < 0) {
        continue;
      }

      // Loop on the four planes and fill histograms accordingly
      for (int ich = 0; ich < 4; ++ich) {

        // Check if BP/NBP has been fired by the track
        bool isFiredBP = trk.isFiredChamber(ich, 0);
        bool isFiredNBP = trk.isFiredChamber(ich, 1);

        // Plane
        histos.fill(HIST("nTotperPlane"), ich); // All counts - plane
        if (isPbPb)
          histos.fill(HIST("hSparseCentFiredTotperPlane"), ich, cent, pt, eta, phi);
        else
          histos.fill(HIST("hSparseCentFiredTotperPlane"), ich, pt, eta, phi);

        if (isFiredBP) {
          histos.fill(HIST("nFiredBPperPlane"), ich); // BP - Plane
          if (isPbPb)
            histos.fill(HIST("hSparseCentFiredBPperPlane"), ich, cent, pt, eta, phi);
          else
            histos.fill(HIST("hSparseCentFiredBPperPlane"), ich, pt, eta, phi);
        }
        if (isFiredNBP) {
          histos.fill(HIST("nFiredNBPperPlane"), ich); // NBP - Plane
          if (isPbPb)
            histos.fill(HIST("hSparseCentFiredNBPperPlane"), ich, cent, pt, eta, phi);
          else
            histos.fill(HIST("hSparseCentFiredNBPperPlane"), ich, pt, eta, phi);
        }
        if (isFiredBP && isFiredNBP) {
          histos.fill(HIST("nFiredBothperPlane"), ich); // Both planes - plane
          if (isPbPb)
            histos.fill(HIST("hSparseCentFiredBothperPlane"), ich, cent, pt, eta, phi);
          else
            histos.fill(HIST("hSparseCentFiredBothperPlane"), ich, pt, eta, phi);
        }

        if (effFlag < 2) {
          continue;
        }

        // Get RPC id
        auto deId = o2::mid::detparams::getDEId(isRight, ich, rpcLine);

        // RPC
        histos.fill(HIST("nTotperRPC"), deId); // All counts - RPC
        if (isPbPb)
          histos.fill(HIST("hSparseCentFiredTotperRPC"), deId, cent, pt, eta, phi);
        else
          histos.fill(HIST("hSparseCentFiredTotperRPC"), deId, pt, eta, phi);

        if (isFiredBP) {
          histos.fill(HIST("nFiredBPperRPC"), deId); // BP - RPC
          if (isPbPb)
            histos.fill(HIST("hSparseCentFiredBPperRPC"), deId, cent, pt, eta, phi);
          else
            histos.fill(HIST("hSparseCentFiredBPperRPC"), deId, pt, eta, phi);
        }
        if (isFiredNBP) {
          histos.fill(HIST("nFiredNBPperRPC"), deId); // NBP - RPC
          if (isPbPb)
            histos.fill(HIST("hSparseCentFiredNBPperRPC"), deId, cent, pt, eta, phi);
          else
            histos.fill(HIST("hSparseCentFiredNBPperRPC"), deId, pt, eta, phi);
        }
        if (isFiredBP && isFiredNBP) {
          histos.fill(HIST("nFiredBothperRPC"), deId); // Both planes - RPC
          if (isPbPb)
            histos.fill(HIST("hSparseCentFiredBothperRPC"), deId, cent, pt, eta, phi);
          else
            histos.fill(HIST("hSparseCentFiredBothperRPC"), deId, pt, eta, phi);
        }

        if (effFlag < 3) {
          continue;
        }

        // Get fired column and line -> needed for LB calculation
        auto firedColumn = trk.getFiredColumnId();
        auto firedLine = trk.getFiredLineId();
        // Get LB ID
        auto LB = ich * o2::mid::detparams::NLocalBoards + mapping.getBoardId(firedLine, firedColumn, deId);

        // LB
        histos.fill(HIST("nTotperBoard"), LB); // All counts - LB
        if (isPbPb)
          histos.fill(HIST("hSparseCentFiredTotperBoard"), LB, cent, pt, eta, phi);
        else
          histos.fill(HIST("hSparseCentFiredTotperBoard"), LB, pt, eta, phi);

        if (isFiredBP) {
          histos.fill(HIST("nFiredBPperBoard"), LB); // BP - LB
          if (isPbPb)
            histos.fill(HIST("hSparseCentFiredBPperBoard"), LB, cent, pt, eta, phi);
          else
            histos.fill(HIST("hSparseCentFiredBPperBoard"), LB, pt, eta, phi);
        }
        if (isFiredNBP) {
          histos.fill(HIST("nFiredNBPperBoard"), LB); // NBP - LB
          if (isPbPb)
            histos.fill(HIST("hSparseCentFiredNBPperBoard"), LB, cent, pt, eta, phi);
          else
            histos.fill(HIST("hSparseCentFiredNBPperBoard"), LB, pt, eta, phi);
        }
        if (isFiredBP && isFiredNBP) {
          histos.fill(HIST("nFiredBothperBoard"), LB); // Both Planes - LB
          if (isPbPb)
            histos.fill(HIST("hSparseCentFiredBothperBoard"), LB, cent, pt, eta, phi);
          else
            histos.fill(HIST("hSparseCentFiredBothperBoard"), LB, pt, eta, phi);
        }
      }
    }

  } // end of runMidEffCounters

  // void processMidEffCounter(aod::ReducedEvents::iterator const& event, soa::Filtered<MyMuonTracks> const& muons)
  void processMidEffCounter(MyEvents::iterator const& event, soa::Filtered<MyMuonTracks> const& muons)
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
