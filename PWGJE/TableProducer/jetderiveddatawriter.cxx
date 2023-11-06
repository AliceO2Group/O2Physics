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

/// \file jetderiveddatawriter.cxx
/// \brief Task to skim jet framework tables (JetCollisions, JetTracks, JetClusters, ...)
/// while adjusting indices accordingly
///
/// \author Jochen Klein <jochen.klein@cern.ch>
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <algorithm>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetDerivedDataWriter {

  Configurable<float> jetPtMin{"jetPtMin", 0.0, "Minimum jet pt to accept event"};

  Produces<aod::StoredJCollisions> storedJetCollisions;
  Produces<aod::StoredJTracks> storedJetTracks;

  Produces<aod::JBCs> storedJBCsTable;
  Produces<aod::JBCPIs> storedJBCParentIndexTable;
  Produces<aod::JCollisions> storedJCollisionsTable;
  Produces<aod::JCollisionPIs> storedJCollisionsParentIndexTable;
  Produces<aod::JCollisionBCs> storedJCollisionsBunchCrossingIndexTable;
  Produces<aod::JChTrigSels> storedJChargedTriggerSelsTable;
  Produces<aod::JFullTrigSels> storedJFullTriggerSelsTable;
  Produces<aod::JMcCollisionLbs> storedJMcCollisionsLabelTable;
  Produces<aod::JMcCollisions> storedJMcCollisionsTable;
  Produces<aod::JMcCollisionPIs> storedJMcCollisionsParentIndexTable;
  Produces<aod::JTracks> storedJTracksTable;
  Produces<aod::JTrackPIs> storedJTracksParentIndexTable;
  Produces<aod::JMcTrackLbs> storedJMcTracksLabelTable;
  Produces<aod::JMcParticles> storedJMcParticlesTable;
  Produces<aod::JMcParticlePIs> storedJParticlesParentIndexTable;
  Produces<aod::JClusters> storedJClustersTable;
  Produces<aod::JClusterTracks> storedJClustersMatchedTracksTable;

  std::vector<bool> collFlag;
  std::vector<int> bcIndicies;
  std::map<int32_t, int32_t> bcMapping;

  bool acceptCollision(aod::JCollision const& collision)
  {
    return true;
  }

  // explicit process function used to reset acceptance flags
  void process(aod::JCollisions const& collisions)
  {
    LOGF(info, "processing %d collisions", collisions.size());
    collFlag.reserve(collisions.size());
    std::fill(collFlag.begin(), collFlag.end(), false);
  }

  template <typename Jets>
  void processJets(Jets& jets)
  {
    for (const auto& jet : jets) {
      if (jet.pt() > jetPtMin) {
        collFlag[jet.collisionId()] = true;
      }
    }
  }
// TODO: replace with PROCESS_SWITCH_FULL when available
#define PROCESS_SWITCH_JKL(_Class_, _Method_, _Name_, _Help_, _Default_) \
  decltype(ProcessConfigurable{&_Class_ ::_Method_, #_Name_, _Default_, _Help_}) do##_Name_ = ProcessConfigurable{&_Class_ ::_Method_, #_Name_, _Default_, _Help_};
  PROCESS_SWITCH_JKL(JetDerivedDataWriter, processJets<aod::ChargedJets>, processChargedJets, "process charged jets", true);
  PROCESS_SWITCH_JKL(JetDerivedDataWriter, processJets<aod::NeutralJets>, processNeutralJets, "process neutral jets", false);
  PROCESS_SWITCH_JKL(JetDerivedDataWriter, processJets<aod::D0ChargedJets>, processD0ChargedJets, "process D0 charged jets", false);
  PROCESS_SWITCH_JKL(JetDerivedDataWriter, processJets<aod::LcChargedJets>, processLcChargedJets, "process Lc charged jets", false);

  void processData(soa::Join<aod::JCollision, aod::JCollisionPIs, aod::JCollisionBCs, aod::JChTrigSels, aod::JFullTrigSels>::iterator const& collision, aod::JBCs bcs, soa::Join<aod::JTracks, aod::JTrackPIs> const& tracks, soa::Join<aod::JClusters, aod::JClusterTracks> const& clusters)
  {
    std::map<int32_t, int32_t> trackMapping;

    if (collFlag[collision.globalIndex()]) {

      auto bc = collision.bc_as<aod::JBCs>;
      if (std::find(bcIndicies.begin(), bcIndicies.end(), bc.globalIndex()) == bcIndicies.end()) {
        storedJBCsTable(bc.runNumber(), bc.globalBC(), bc.timestamp());
        storedJBCParentIndexTable(bc.bcId());
        bcIndicies.push_back(bc.globalIndex());
        bcMapping.insert(make_pair(bc.globalIndex(), storedJBCsTable.lastIndex()));
      }

      storedJCollisionsTable(collision.posZ(), collision.eventSel(), collision.alias_raw());
      storedJCollisionsParentIndexTable(collision.collisionId());
      int32_t storedBCID = -1;
      auto JBCIndex = bcMapping.find(collision.bcId());
      if (JBCIndex != bcMapping.end()) {
      storedBCID = JBCIndex->second);
      }
      storedJCollisionsBunchCrossingIndexTable(storedBCID);
      storedJChargedTriggerSelsTable(collision.chargedTriggerSel());
      storedJFullTriggerSelsTable(collision.fullTriggerSel());

      for (const auto& track : tracks) {
      storedJTracksTable(storedJCollisionsTable.lastIndex(), track.pt(), track.eta(), track.phi(), track.energy(), track.trackSel());
      storedJTracksParentIndexTable(track.trackId());
      trackMapping.insert(make_pair(track.globalIndex(), storedJTracksTable.lastIndex()));
      }

      for (const auto& cluster : clusters) {
      storedJClustersTable(storedJCollisionsTable.lastIndex(), cluster.id(), cluster.energy(), cluster.coreEnergy(), cluster.rawEnergy(),
                           cluster.eta(), cluster.phi(), cluster.m02(), cluster.m20(), cluster.nCells(), cluster.time(), cluster.isExotic(), cluster.distanceToBadChannel(),
                           cluster.nlm(), cluster.definition(), cluster.leadingCellEnergy(), cluster.subleadingCellEnergy(), cluster.leadingCellNumber(), cluster.subleadingCellNumber());

      std::vector<int> clusterStoredJTrackIDs;
      for (const auto& clusterTrack : cluster.matchedTracks_as<soa::Join<aod::JTracks, aod::JTrackPIs>>) {
        auto JtrackIndex = trackMapping.find(clusterTrack.globalIndex());
        if (JtrackIndex != trackMapping.end()) {
          clusterStoredJTrackIDs.push_back(JtrackIndex->second);
        }
      }
      jClustersMatchedTracksTable(clusterStoredJTrackIDs);
      }
    }
  }
  // process switch for output writing must be last
  // to run after all jet selections
  PROCESS_SWITCH(JetDerivedDataWriter, processData, "write out data output tables", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetDerivedDataWriter>(cfgc, TaskName{"jet-deriveddata-writer"}));

  return WorkflowSpec{tasks};
}
