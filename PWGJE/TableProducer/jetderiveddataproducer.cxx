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

// task to produce a self contained data format for jet analyses from the full AO2D
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <vector>
#include <string>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/RecoDecay.h"
#include "PWGJE/DataModel/EMCALClusters.h"

#include "EventFiltering/filterTables.h"

#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetDerivedDataProducerTask {
  Produces<aod::JDummys> jDummysTable;
  Produces<aod::JBCs> jBCsTable;
  Produces<aod::JBCPIs> jBCParentIndexTable;
  Produces<aod::JCollisions> jCollisionsTable;
  Produces<aod::JCollisionPIs> jCollisionsParentIndexTable;
  Produces<aod::JCollisionBCs> jCollisionsBunchCrossingIndexTable;
  Produces<aod::JMcCollisionLbs> jMcCollisionsLabelTable;
  Produces<aod::JMcCollisions> jMcCollisionsTable;
  Produces<aod::JMcCollisionPIs> jMcCollisionsParentIndexTable;
  Produces<aod::JTracks> jTracksTable;
  Produces<aod::JTrackPIs> jTracksParentIndexTable;
  Produces<aod::JMcTrackLbs> jMcTracksLabelTable;
  Produces<aod::JMcParticles> jMcParticlesTable;
  Produces<aod::JMcParticlePIs> jParticlesParentIndexTable;
  Produces<aod::JClusters> jClustersTable;
  Produces<aod::JClusterPIs> jClustersParentIndexTable;
  Produces<aod::JClusterTracks> jClustersMatchedTracksTable;

  Preslice<aod::EMCALClusterCells> perClusterCells = aod::emcalclustercell::emcalclusterId;
  Preslice<aod::EMCALMatchedTracks> perClusterTracks = aod::emcalclustercell::emcalclusterId;

  void init(InitContext const&)
  {
  }

  void processBunchCossings(soa::Join<aod::BCs, aod::Timestamps>::iterator const& bc)
  {
    jBCsTable(bc.runNumber(), bc.globalBC(), bc.timestamp());
    jBCParentIndexTable(bc.globalIndex());
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processBunchCossings, "produces derived bunch crossing table", false);

  void processCollisions(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision)
  {
    jCollisionsTable(collision.posZ(), JetDerivedDataUtilities::setEventSelectionBit(collision), collision.alias_raw());
    jCollisionsParentIndexTable(collision.globalIndex());
    jCollisionsBunchCrossingIndexTable(collision.bcId());
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processCollisions, "produces derived collision tables", true);

  void processMcCollisionLabels(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision)
  {

    if (collision.has_mcCollision()) {
      jMcCollisionsLabelTable(collision.mcCollisionId());
    } else {
      jMcCollisionsLabelTable(-1);
    }
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processMcCollisionLabels, "produces derived MC collision labels table", false);

  void processMcCollisions(aod::McCollision const& McCollision)
  {
    jMcCollisionsTable(McCollision.posZ(), McCollision.weight());
    jMcCollisionsParentIndexTable(McCollision.globalIndex());
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processMcCollisions, "produces derived MC collision table", false);

  void processTracks(soa::Join<aod::Tracks, aod::TrackSelection>::iterator const& track)
  {
    jTracksTable(track.collisionId(), track.pt(), track.eta(), track.phi(), JetDerivedDataUtilities::trackEnergy(track), track.sign(), JetDerivedDataUtilities::setTrackSelectionBit(track));
    jTracksParentIndexTable(track.globalIndex());
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processTracks, "produces derived track table", true);

  void processMcTrackLabels(soa::Join<aod::Tracks, aod::McTrackLabels>::iterator const& track)
  {
    if (track.has_mcParticle()) {
      jMcTracksLabelTable(track.mcParticleId());
    } else {
      jMcTracksLabelTable(-1);
    }
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processMcTrackLabels, "produces derived track labels table", false);

  void processParticles(aod::McParticle const& particle)
  {
    std::vector<int> mothersId;
    if (particle.has_mothers()) {
      auto mothersIdTemps = particle.mothersIds();
      for (auto mothersIdTemp : mothersIdTemps) {
        mothersId.push_back(mothersIdTemp);
      }
    }
    int daughtersId[2] = {-1, -1};
    auto i = 0;
    if (particle.has_daughters()) {
      for (auto daughterId : particle.daughtersIds()) {
        if (i > 1) {
          break;
        }
        daughtersId[i] = daughterId;
        i++;
      }
    }
    jMcParticlesTable(particle.mcCollisionId(), particle.pt(), particle.eta(), particle.phi(), particle.y(), particle.e(), particle.pdgCode(), particle.getGenStatusCode(), particle.getHepMCStatusCode(), particle.isPhysicalPrimary(), mothersId, daughtersId);
    jParticlesParentIndexTable(particle.globalIndex());
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processParticles, "produces derived parrticle table", false);

  void processClusters(aod::Collision const& collision, aod::EMCALClusters const& clusters, aod::EMCALClusterCells const& cells, aod::Calos const& calos, aod::EMCALMatchedTracks const& matchedTracks, aod::Tracks const& tracks)
  {

    for (auto cluster : clusters) {

      auto const clusterCells = cells.sliceBy(perClusterCells, cluster.globalIndex());

      float leadingCellEnergy = -1.0;
      float subleadingCellEnergy = -1.0;
      float cellAmplitude = -1.0;
      int leadingCellNumber = -1;
      int subleadingCellNumber = -1;
      int cellNumber = -1;
      for (auto const& clutserCell : clusterCells) {
        cellAmplitude = clutserCell.calo().amplitude();
        cellNumber = clutserCell.calo().cellNumber();
        if (cellAmplitude > subleadingCellEnergy) {
          subleadingCellEnergy = cellAmplitude;
          subleadingCellNumber = cellNumber;
        }
        if (subleadingCellEnergy > leadingCellEnergy) {
          std::swap(leadingCellEnergy, subleadingCellEnergy);
          std::swap(leadingCellNumber, subleadingCellNumber);
        }
      }

      jClustersTable(cluster.collisionId(), cluster.id(), cluster.energy(), cluster.coreEnergy(), cluster.rawEnergy(), cluster.eta(), cluster.phi(), cluster.m02(), cluster.m20(), cluster.nCells(), cluster.time(), cluster.isExotic(), cluster.distanceToBadChannel(), cluster.nlm(), cluster.definition(), leadingCellEnergy, subleadingCellEnergy, leadingCellNumber, subleadingCellNumber);
      jClustersParentIndexTable(cluster.globalIndex());

      auto const clusterTracks = matchedTracks.sliceBy(perClusterTracks, cluster.globalIndex());
      std::vector<int> clusterTrackIDs;
      for (const auto& clusterTrack : clusterTracks) {
        clusterTrackIDs.push_back(clusterTrack.trackId());
      }
      jClustersMatchedTracksTable(clusterTrackIDs);
    }
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processClusters, "produces derived cluster tables", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetDerivedDataProducerTask>(cfgc, TaskName{"jet-deriveddata-producer"})};
}
