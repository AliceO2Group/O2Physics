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

#include <climits>
#include <cstdlib>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "EMCALBase/Geometry.h"
#include "EMCALCalib/BadChannelMap.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/EMCALMatchedTracks.h"
#include "DataFormatsEMCAL/Cell.h"
#include "DataFormatsEMCAL/Constants.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"

#include "CommonDataFormat/InteractionRecord.h"

// \struct EmcalMatchedTracksTask
/// \brief Simple table producer task for EMCal matched tracks
/// \author Marvin Hemmer <marvin.hemmer@cern.ch>
/// \since 16.11.2023
///
/// This task is meant to be used for monitoring the matching between global tracks and EMCal clusters
/// properties, such as:
/// - cluster energy over track momentum
/// - difference in eta
/// - difference in phi
using namespace o2::framework;
using namespace o2::framework::expressions;

using myCollisions = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
using myCollision = myCollisions::iterator;
using myBCs = o2::soa::Join<o2::aod::BCsWithTimestamps, o2::aod::BcSels>;
using selectedClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;
using myTracksPID = o2::soa::Filtered<o2::soa::Join<o2::aod::pidTPCFullEl, o2::aod::pidTPCFullPi, o2::aod::pidTOFFullEl, o2::aod::pidTOFFullPi, o2::aod::FullTracks, o2::aod::TracksCov, o2::aod::TrackSelection>>;

struct EmcalMatchedTracksTask {
  Produces<o2::aod::EmcalMTs> clustersmatchedtracks;

  HistogramRegistry mHistManager{"EmcalMatchedTracksHistograms", {}, OutputObjHandlingPolicy::AnalysisObject};
  o2::emcal::Geometry* mGeometry = nullptr;

  Preslice<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalclustercell::emcalclusterId;

  // configurable parameters
  Configurable<bool> mDoEventSel{"doEventSel", 0, "demand kINT7"};
  Configurable<double> mVertexCut{"vertexCut", 10., "apply z-vertex cut with value in cm"};
  Configurable<int> mClusterDefinition{"clusterDefinition", 10, "cluster definition to be selected, e.g. 10=kV3Default"};
  Configurable<bool> hasPropagatedTracks{"hasPropagatedTracks", false, "temporary flag, only set to true when running over data which has the tracks propagated to EMCal/PHOS!"};
  Configurable<float> minTime{"minTime", -25., "Minimum cluster time for time cut"};
  Configurable<float> maxTime{"maxTime", +20., "Maximum cluster time for time cut"};
  Configurable<float> minM02{"minM02", 0.1, "Minimum M02 for M02 cut"};
  Configurable<float> maxM02{"maxM02", 0.9, "Maximum M02 for M02 cut"};
  Configurable<float> minTrackPt{"minTrackPt", 0.15, "Minimum pT for tracks"};
  Configurable<float> minDEta{"minDEta", 0.1, "Minimum dEta between track and cluster"};
  Configurable<float> minDPhi{"minDPhi", 0.1, "Minimum dPhi between track and cluster"};

  std::vector<int> mSelectBCIDs;

  /// \brief Create output histograms and initialize geometry
  void init(InitContext const&)
  {
    // create histograms
    using o2HistType = HistType;

    // load geometry just in case we need it
    mGeometry = o2::emcal::Geometry::GetInstanceFromRunNumber(300000);

    // event properties
    mHistManager.add("eventsAll", "Number of events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventsSelected", "Number of events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventVertexZAll", "z-vertex of event (all events)", o2HistType::kTH1F, {{200, -20, 20}});
    mHistManager.add("eventVertexZSelected", "z-vertex of event (selected events)", o2HistType::kTH1F, {{200, -20, 20}});
  }

  // define cluster filter. It selects only those clusters which are of the type
  Filter clusterDefinitionSelection = (o2::aod::emcalcluster::definition == mClusterDefinition) && (o2::aod::emcalcluster::time >= minTime) && (o2::aod::emcalcluster::time <= maxTime) && (o2::aod::emcalcluster::m02 > minM02) && (o2::aod::emcalcluster::m02 < maxM02);
  // cut on tracks
  Filter trackSelection = (o2::aod::track::pt >= minTrackPt);

  /// \brief Process EMCAL clusters that are matched to a collisions
  void processCollisions(myCollision const& theCollision, myBCs const& bcs, selectedClusters const& clusters, o2::aod::EMCALMatchedTracks const& matchedtracks, myTracksPID const& alltracks)
  {
    mHistManager.fill(HIST("eventsAll"), 1);

    o2::InteractionRecord eventIR;
    eventIR.setFromLong(theCollision.bc_as<myBCs>().globalBC());

    // do event selection if mDoEventSel is specified
    // currently the event selection is hard coded to kINT7
    // but other selections are possible that are defined in TriggerAliases.h
    bool isSelected = true;
    if (mDoEventSel) {
      if (theCollision.bc().runNumber() < 300000) {
        if (!theCollision.alias_bit(kINT7)) {
          isSelected = false;
        }
      } else {
        if (!theCollision.alias_bit(kTVXinEMC)) {
          isSelected = false;
        }
      }
    }
    if (!isSelected) {
      LOG(debug) << "Event not selected because it is not kINT7 or does not have EMCAL in readout, skipping";
      return;
    }
    mHistManager.fill(HIST("eventVertexZAll"), theCollision.posZ());
    if (mVertexCut > 0 && TMath::Abs(theCollision.posZ()) > mVertexCut) {
      LOG(debug) << "Event not selected because of z-vertex cut z= " << theCollision.posZ() << " > " << mVertexCut << " cm, skipping";
      return;
    }
    mHistManager.fill(HIST("eventsSelected"), 1);
    mHistManager.fill(HIST("eventVertexZSelected"), theCollision.posZ());

    // skip events with no clusters
    if (clusters.size() == 0) {
      return;
    }
    // loop over all clusters from accepted collision
    for (const auto& cluster : clusters) {
      auto tracksofcluster = matchedtracks.sliceBy(perClusterMatchedTracks, cluster.globalIndex());
      // skip clusters with no matched tracks!
      if (tracksofcluster.size() == 0) {
        continue;
      } else if (tracksofcluster.size() == 1) {
        if (fabs(tracksofcluster.iteratorAt(0).track_as<myTracksPID>().trackEtaEmcal() - cluster.eta()) >= minDEta || fabs(tracksofcluster.iteratorAt(0).track_as<myTracksPID>().trackPhiEmcal() - cluster.phi()) >= minDPhi) {
          continue;
        }
        clustersmatchedtracks(eventIR.orbit, theCollision.bc_as<myBCs>().timestamp(), theCollision.bc_as<myBCs>().runNumber(),
                              tracksofcluster.iteratorAt(0).track_as<myTracksPID>().x(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().alpha(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().p(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().signed1Pt(),
                              tracksofcluster.iteratorAt(0).track_as<myTracksPID>().y(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().z(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().snp(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().tgl(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().pt(),
                              tracksofcluster.iteratorAt(0).track_as<myTracksPID>().sigmaY(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().sigmaZ(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().sigmaSnp(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().sigmaTgl(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().sigma1Pt(),
                              tracksofcluster.iteratorAt(0).track_as<myTracksPID>().eta(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().phi(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().trackEtaEmcal(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().trackPhiEmcal(),
                              tracksofcluster.iteratorAt(0).track_as<myTracksPID>().trackEtaEmcal() - cluster.eta(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().trackPhiEmcal() - cluster.phi(),
                              tracksofcluster.iteratorAt(0).track_as<myTracksPID>().itsNCls(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().tofExpMom(),
                              tracksofcluster.iteratorAt(0).track_as<myTracksPID>().tpcNSigmaEl(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().tpcNSigmaPi(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().tofNSigmaEl(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().tofNSigmaPi(),
                              0.f, 0.f, 0.f, 0.f,
                              0.f, 0.f, 0.f, 0.f, 0.f,
                              0.f, 0.f, 0.f, 0.f, 0.f,
                              0.f, 0.f, 0.f, 0.f,
                              0.f, 0.f,
                              0, 0.f,
                              0.f, 0.f, 0.f, 0.f,
                              cluster.energy(), cluster.eta(), cluster.phi(), cluster.m02(), cluster.nCells(), cluster.time());

      } else if (tracksofcluster.size() >= 2) { // do at least two cluster pre matched
        if (fabs(tracksofcluster.iteratorAt(0).track_as<myTracksPID>().trackEtaEmcal() - cluster.eta()) >= minDEta || fabs(tracksofcluster.iteratorAt(0).track_as<myTracksPID>().trackPhiEmcal() - cluster.phi()) >= minDPhi) {
          // if not even the first track is within tighter matching window, go to next cluster
          continue;
        } else if (fabs(tracksofcluster.iteratorAt(1).track_as<myTracksPID>().trackEtaEmcal() - cluster.eta()) >= minDEta || fabs(tracksofcluster.iteratorAt(1).track_as<myTracksPID>().trackPhiEmcal() - cluster.phi()) >= minDPhi) {
          // if only the first track is within tighter matching window, just write that one to table
          clustersmatchedtracks(eventIR.orbit, theCollision.bc_as<myBCs>().timestamp(), theCollision.bc_as<myBCs>().runNumber(),
                                tracksofcluster.iteratorAt(0).track_as<myTracksPID>().x(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().alpha(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().p(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().signed1Pt(),
                                tracksofcluster.iteratorAt(0).track_as<myTracksPID>().y(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().z(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().snp(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().tgl(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().pt(),
                                tracksofcluster.iteratorAt(0).track_as<myTracksPID>().sigmaY(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().sigmaZ(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().sigmaSnp(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().sigmaTgl(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().sigma1Pt(),
                                tracksofcluster.iteratorAt(0).track_as<myTracksPID>().eta(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().phi(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().trackEtaEmcal(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().trackPhiEmcal(),
                                tracksofcluster.iteratorAt(0).track_as<myTracksPID>().trackEtaEmcal() - cluster.eta(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().trackPhiEmcal() - cluster.phi(),
                                tracksofcluster.iteratorAt(0).track_as<myTracksPID>().itsNCls(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().tofExpMom(),
                                tracksofcluster.iteratorAt(0).track_as<myTracksPID>().tpcNSigmaEl(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().tpcNSigmaPi(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().tofNSigmaEl(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().tofNSigmaPi(),
                                0.f, 0.f, 0.f, 0.f,
                                0.f, 0.f, 0.f, 0.f, 0.f,
                                0.f, 0.f, 0.f, 0.f, 0.f,
                                0.f, 0.f, 0.f, 0.f,
                                0.f, 0.f,
                                0, 0.f,
                                0.f, 0.f, 0.f, 0.f,
                                cluster.energy(), cluster.eta(), cluster.phi(), cluster.m02(), cluster.nCells(), cluster.time());
        } else {
          // if both the first and second track are within tighter matching window, write both down
          clustersmatchedtracks(eventIR.orbit, theCollision.bc_as<myBCs>().timestamp(), theCollision.bc_as<myBCs>().runNumber(),
                                tracksofcluster.iteratorAt(0).track_as<myTracksPID>().x(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().alpha(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().p(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().signed1Pt(),
                                tracksofcluster.iteratorAt(0).track_as<myTracksPID>().y(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().z(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().snp(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().tgl(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().pt(),
                                tracksofcluster.iteratorAt(0).track_as<myTracksPID>().sigmaY(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().sigmaZ(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().sigmaSnp(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().sigmaTgl(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().sigma1Pt(),
                                tracksofcluster.iteratorAt(0).track_as<myTracksPID>().eta(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().phi(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().trackEtaEmcal(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().trackPhiEmcal(),
                                tracksofcluster.iteratorAt(0).track_as<myTracksPID>().trackEtaEmcal() - cluster.eta(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().trackPhiEmcal() - cluster.phi(),
                                tracksofcluster.iteratorAt(0).track_as<myTracksPID>().itsNCls(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().tofExpMom(),
                                tracksofcluster.iteratorAt(0).track_as<myTracksPID>().tpcNSigmaEl(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().tpcNSigmaPi(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().tofNSigmaEl(), tracksofcluster.iteratorAt(0).track_as<myTracksPID>().tofNSigmaPi(),
                                tracksofcluster.iteratorAt(1).track_as<myTracksPID>().x(), tracksofcluster.iteratorAt(1).track_as<myTracksPID>().alpha(), tracksofcluster.iteratorAt(1).track_as<myTracksPID>().p(), tracksofcluster.iteratorAt(1).track_as<myTracksPID>().signed1Pt(),
                                tracksofcluster.iteratorAt(1).track_as<myTracksPID>().y(), tracksofcluster.iteratorAt(1).track_as<myTracksPID>().z(), tracksofcluster.iteratorAt(1).track_as<myTracksPID>().snp(), tracksofcluster.iteratorAt(1).track_as<myTracksPID>().tgl(), tracksofcluster.iteratorAt(1).track_as<myTracksPID>().pt(),
                                tracksofcluster.iteratorAt(1).track_as<myTracksPID>().sigmaY(), tracksofcluster.iteratorAt(1).track_as<myTracksPID>().sigmaZ(), tracksofcluster.iteratorAt(1).track_as<myTracksPID>().sigmaSnp(), tracksofcluster.iteratorAt(1).track_as<myTracksPID>().sigmaTgl(), tracksofcluster.iteratorAt(1).track_as<myTracksPID>().sigma1Pt(),
                                tracksofcluster.iteratorAt(1).track_as<myTracksPID>().eta(), tracksofcluster.iteratorAt(1).track_as<myTracksPID>().phi(), tracksofcluster.iteratorAt(1).track_as<myTracksPID>().trackEtaEmcal(), tracksofcluster.iteratorAt(1).track_as<myTracksPID>().trackPhiEmcal(),
                                tracksofcluster.iteratorAt(1).track_as<myTracksPID>().trackEtaEmcal() - cluster.eta(), tracksofcluster.iteratorAt(1).track_as<myTracksPID>().trackPhiEmcal() - cluster.phi(),
                                tracksofcluster.iteratorAt(1).track_as<myTracksPID>().itsNCls(), tracksofcluster.iteratorAt(1).track_as<myTracksPID>().tofExpMom(),
                                tracksofcluster.iteratorAt(1).track_as<myTracksPID>().tpcNSigmaEl(), tracksofcluster.iteratorAt(1).track_as<myTracksPID>().tpcNSigmaPi(), tracksofcluster.iteratorAt(1).track_as<myTracksPID>().tofNSigmaEl(), tracksofcluster.iteratorAt(1).track_as<myTracksPID>().tofNSigmaPi(), cluster.energy(), cluster.eta(), cluster.phi(), cluster.m02(), cluster.nCells(), cluster.time());
        }
      }
    } // end of loop over clusters
  }
  PROCESS_SWITCH(EmcalMatchedTracksTask, processCollisions, "Process clusters from collision", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<EmcalMatchedTracksTask>(cfgc, TaskName{"emc-matchedtracks-writer"})};
  return workflow;
}
