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

  Configurable<float> chargedJetPtMin{"chargedJetPtMin", 0.0, "Minimum charged jet pt to accept event"};
  Configurable<float> chargedMCPJetPtMin{"chargedMCPJetPtMin", 0.0, "Minimum charged mcp jet pt to accept event"};
  Configurable<float> neutralJetPtMin{"neutralJetPtMin", 0.0, "Minimum charged jet pt to accept event"};
  Configurable<float> fullJetPtMin{"fullJetPtMin", 0.0, "Minimum full jet pt to accept event"};
  Configurable<float> chargedD0JetPtMin{"chargedD0JetPtMin", 0.0, "Minimum charged D0 jet pt to accept event"};
  Configurable<float> chargedLcJetPtMin{"chargedLcJetPtMin", 0.0, "Minimum charged Lc jet pt to accept event"};

  Produces<aod::StoredJDummys> storedJDummysTable;
  Produces<aod::StoredJBCs> storedJBCsTable;
  Produces<aod::StoredJBCPIs> storedJBCParentIndexTable;
  Produces<aod::StoredJCollisions> storedJCollisionsTable;
  Produces<aod::StoredJCollisionPIs> storedJCollisionsParentIndexTable;
  Produces<aod::StoredJCollisionBCs> storedJCollisionsBunchCrossingIndexTable;
  Produces<aod::StoredJChTrigSels> storedJChargedTriggerSelsTable;
  Produces<aod::StoredJFullTrigSels> storedJFullTriggerSelsTable;
  Produces<aod::StoredJMcCollisionLbs> storedJMcCollisionsLabelTable;
  Produces<aod::StoredJMcCollisions> storedJMcCollisionsTable;
  Produces<aod::StoredJMcCollisionPIs> storedJMcCollisionsParentIndexTable;
  Produces<aod::StoredJTracks> storedJTracksTable;
  Produces<aod::StoredJTrackPIs> storedJTracksParentIndexTable;
  Produces<aod::StoredJMcTrackLbs> storedJMcTracksLabelTable;
  Produces<aod::StoredJMcParticles> storedJMcParticlesTable;
  Produces<aod::StoredJMcParticlePIs> storedJParticlesParentIndexTable;
  Produces<aod::StoredJClusters> storedJClustersTable;
  Produces<aod::StoredJClusterPIs> storedJClustersParentIndexTable;
  Produces<aod::StoredJClusterTracks> storedJClustersMatchedTracksTable;

  PresliceUnsorted<soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::JCollisionBCs, aod::JChTrigSels, aod::JFullTrigSels, aod::JMcCollisionLbs>> CollisionsPerMcCollision = aod::jmccollisionlb::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::JMcParticles, aod::JMcParticlePIs>> ParticlesPerMcCollision = aod::jmcparticle::mcCollisionId;
  Preslice<soa::Join<aod::JTracks, aod::JTrackPIs, aod::JMcTrackLbs>> TracksPerCollision = aod::jtrack::collisionId;
  Preslice<soa::Join<aod::JClusters, aod::JClusterPIs, aod::JClusterTracks>> ClustersPerCollision = aod::jcluster::collisionId;

  std::vector<bool> collisionFlag;
  std::vector<bool> McCollisionFlag;
  std::vector<int> bcIndicies;
  std::map<int32_t, int32_t> bcMapping;

  bool acceptCollision(aod::JCollision const& collision)
  {
    return true;
  }

  void processCollisions(aod::JCollisions const& collisions)
  {
    collisionFlag.reserve(collisions.size());
    std::fill(collisionFlag.begin(), collisionFlag.end(), false);
  }

  void processMcCollisions(aod::JMcCollisions const& Mccollisions)
  {
    McCollisionFlag.reserve(Mccollisions.size());
    std::fill(McCollisionFlag.begin(), McCollisionFlag.end(), false);
  }

  template <typename T>
  void processJets(T& jets)
  {
    float jetPtMin = 0.0;
    if constexpr (std::is_same_v<std::decay_t<T>, aod::ChargedJets> || std::is_same_v<std::decay_t<T>, aod::ChargedMCDetectorLevelJets>) {
      jetPtMin = chargedJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::ChargedMCParticleLevelJets>) {
      jetPtMin = chargedMCPJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::NeutralJets>) {
      jetPtMin = neutralJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::FullJets>) {
      jetPtMin = fullJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::D0ChargedJets>) {
      jetPtMin = chargedD0JetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::LcChargedJets>) {
      jetPtMin = chargedLcJetPtMin;
    } else {
      jetPtMin = 0.0;
    }

    for (const auto& jet : jets) {
      if (jet.pt() >= jetPtMin) {
        if constexpr (std::is_same_v<std::decay_t<T>, aod::ChargedMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::NeutralMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::FullMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::D0ChargedMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::LcChargedMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::BplusChargedMCParticleLevelJets>) {
          McCollisionFlag[jet.mcCollisionId()] = true;
        } else {
          collisionFlag[jet.collisionId()] = true;
        }
      }
    }
  }
// TODO: replace with PROCESS_SWITCH_FULL when available
#define PROCESS_SWITCH_JKL(_Class_, _Method_, _Name_, _Help_, _Default_) \
  decltype(ProcessConfigurable{&_Class_ ::_Method_, #_Name_, _Default_, _Help_}) do##_Name_ = ProcessConfigurable{&_Class_ ::_Method_, #_Name_, _Default_, _Help_};
  PROCESS_SWITCH(JetDerivedDataWriter, processCollisions, "setup the writing for data and MCD", true);
  PROCESS_SWITCH(JetDerivedDataWriter, processMcCollisions, "setup the writing for MCP", false);
  PROCESS_SWITCH_JKL(JetDerivedDataWriter, processJets<aod::ChargedJets>, processChargedJets, "process charged jets", true);
  PROCESS_SWITCH_JKL(JetDerivedDataWriter, processJets<aod::ChargedMCDetectorLevelJets>, processChargedMCDJets, "process charged mcd jets", false);
  PROCESS_SWITCH_JKL(JetDerivedDataWriter, processJets<aod::ChargedMCParticleLevelJets>, processChargedMCPJets, "process charged mcp jets", false);
  PROCESS_SWITCH_JKL(JetDerivedDataWriter, processJets<aod::NeutralJets>, processNeutralJets, "process neutral jets", false);
  PROCESS_SWITCH_JKL(JetDerivedDataWriter, processJets<aod::FullJets>, processFullJets, "process full jets", false);
  PROCESS_SWITCH_JKL(JetDerivedDataWriter, processJets<aod::D0ChargedJets>, processD0ChargedJets, "process D0 charged jets", false);
  PROCESS_SWITCH_JKL(JetDerivedDataWriter, processJets<aod::LcChargedJets>, processLcChargedJets, "process Lc charged jets", false);

  void processDummy(aod::JDummys const& Dummys)
  {
    storedJDummysTable(1);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDummy, "write out dummy output table", true);

  void processData(soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::JCollisionBCs, aod::JChTrigSels, aod::JFullTrigSels>::iterator const& collision, soa::Join<aod::JBCs, aod::JBCPIs> const& bcs, soa::Join<aod::JTracks, aod::JTrackPIs> const& tracks, soa::Join<aod::JClusters, aod::JClusterPIs, aod::JClusterTracks> const& clusters)
  {
    std::map<int32_t, int32_t> trackMapping;

    if (collisionFlag[collision.globalIndex()]) {

      auto bc = collision.bc_as<soa::Join<aod::JBCs, aod::JBCPIs>>();
      if (std::find(bcIndicies.begin(), bcIndicies.end(), bc.globalIndex()) == bcIndicies.end()) {
        storedJBCsTable(bc.runNumber(), bc.globalBC(), bc.timestamp());
        storedJBCParentIndexTable(bc.bcId());
        bcIndicies.push_back(bc.globalIndex());
        bcMapping.insert(std::make_pair(bc.globalIndex(), storedJBCsTable.lastIndex()));
      }

      storedJCollisionsTable(collision.posZ(), collision.eventSel(), collision.alias_raw());
      storedJCollisionsParentIndexTable(collision.collisionId());
      int32_t storedBCID = -1;
      auto JBCIndex = bcMapping.find(collision.bcId());
      if (JBCIndex != bcMapping.end()) {
        storedBCID = JBCIndex->second;
      }
      storedJCollisionsBunchCrossingIndexTable(storedBCID);
      storedJChargedTriggerSelsTable(collision.chargedTriggerSel());
      storedJFullTriggerSelsTable(collision.fullTriggerSel());

      for (const auto& track : tracks) {
        storedJTracksTable(storedJCollisionsTable.lastIndex(), track.pt(), track.eta(), track.phi(), track.energy(), track.sign(), track.trackSel());
        storedJTracksParentIndexTable(track.trackId());
        trackMapping.insert(std::make_pair(track.globalIndex(), storedJTracksTable.lastIndex()));
      }

      for (const auto& cluster : clusters) {
        storedJClustersTable(storedJCollisionsTable.lastIndex(), cluster.id(), cluster.energy(), cluster.coreEnergy(), cluster.rawEnergy(),
                             cluster.eta(), cluster.phi(), cluster.m02(), cluster.m20(), cluster.nCells(), cluster.time(), cluster.isExotic(), cluster.distanceToBadChannel(),
                             cluster.nlm(), cluster.definition(), cluster.leadingCellEnergy(), cluster.subleadingCellEnergy(), cluster.leadingCellNumber(), cluster.subleadingCellNumber());
        storedJClustersParentIndexTable(cluster.clusterId());

        std::vector<int> clusterStoredJTrackIDs;
        for (const auto& clusterTrack : cluster.matchedTracks_as<soa::Join<aod::JTracks, aod::JTrackPIs>>()) {
          auto JtrackIndex = trackMapping.find(clusterTrack.globalIndex());
          if (JtrackIndex != trackMapping.end()) {
            clusterStoredJTrackIDs.push_back(JtrackIndex->second);
          }
        }
        storedJClustersMatchedTracksTable(clusterStoredJTrackIDs);
      }
    }
  }
  // process switch for output writing must be last
  // to run after all jet selections
  PROCESS_SWITCH(JetDerivedDataWriter, processData, "write out data output tables", true);

  void processMC(soa::Join<aod::JMcCollisions, aod::JMcCollisionPIs> const& mcCollisions, soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::JCollisionBCs, aod::JChTrigSels, aod::JFullTrigSels, aod::JMcCollisionLbs> const& collisions, soa::Join<aod::JBCs, aod::JBCPIs> const& bcs, soa::Join<aod::JTracks, aod::JTrackPIs, aod::JMcTrackLbs> const& tracks, soa::Join<aod::JClusters, aod::JClusterPIs, aod::JClusterTracks> const& clusters, soa::Join<aod::JMcParticles, aod::JMcParticlePIs> const& particles)
  {

    std::map<int32_t, int32_t> paticleMapping;
    std::map<int32_t, int32_t> mcCollisionMapping;
    int particleTableIndex = 0;
    for (auto mcCollision : mcCollisions) {
      bool collisionSelected = false;
      const auto collisionsPerMcCollision = collisions.sliceBy(CollisionsPerMcCollision, mcCollision.globalIndex());
      for (auto collision : collisionsPerMcCollision) {
        if (collisionFlag[collision.globalIndex()]) {
          collisionSelected = true;
        }
      }

      if (McCollisionFlag[mcCollision.globalIndex()] || collisionSelected) {

        const auto particlesPerMcCollision = particles.sliceBy(ParticlesPerMcCollision, mcCollision.globalIndex());

        storedJMcCollisionsTable(mcCollision.posZ(), mcCollision.weight());
        storedJMcCollisionsParentIndexTable(mcCollision.mcCollisionId());
        mcCollisionMapping.insert(std::make_pair(mcCollision.globalIndex(), storedJMcCollisionsTable.lastIndex()));

        for (auto particle : particlesPerMcCollision) {
          paticleMapping.insert(std::make_pair(particle.globalIndex(), particleTableIndex));
          particleTableIndex++;
        }
        for (auto particle : particlesPerMcCollision) {

          std::vector<int> mothersId;
          if (particle.has_mothers()) {
            auto mothersIdTemps = particle.mothersIds();
            for (auto mothersIdTemp : mothersIdTemps) {

              auto JMotherIndex = paticleMapping.find(mothersIdTemp);
              if (JMotherIndex != paticleMapping.end()) {
                mothersId.push_back(JMotherIndex->second);
              }
            }
          }
          int daughtersId[2] = {-1, -1};
          auto i = 0;
          if (particle.has_daughters()) {
            for (auto daughterId : particle.daughtersIds()) {
              if (i > 1) {
                break;
              }
              auto JDaughterIndex = paticleMapping.find(daughterId);
              if (JDaughterIndex != paticleMapping.end()) {
                daughtersId[i] = JDaughterIndex->second;
              }
              i++;
            }
          }
          storedJMcParticlesTable(storedJMcCollisionsTable.lastIndex(), particle.pt(), particle.eta(), particle.phi(), particle.y(), particle.e(), particle.pdgCode(), particle.getGenStatusCode(), particle.getHepMCStatusCode(), particle.isPhysicalPrimary(), mothersId, daughtersId);
          storedJParticlesParentIndexTable(particle.mcParticleId());
        }
      }
    }

    for (auto mcCollision : mcCollisions) {
      bool collisionSelected = false;
      const auto collisionsPerMcCollision = collisions.sliceBy(CollisionsPerMcCollision, mcCollision.globalIndex());
      for (auto collision : collisionsPerMcCollision) {
        if (collisionFlag[collision.globalIndex()]) {
          collisionSelected = true;
        }
      }

      if (McCollisionFlag[mcCollision.globalIndex()] || collisionSelected) {

        for (auto collision : collisionsPerMcCollision) {
          std::map<int32_t, int32_t> trackMapping;

          auto bc = collision.bc_as<soa::Join<aod::JBCs, aod::JBCPIs>>();
          if (std::find(bcIndicies.begin(), bcIndicies.end(), bc.globalIndex()) == bcIndicies.end()) {
            storedJBCsTable(bc.runNumber(), bc.globalBC(), bc.timestamp());
            storedJBCParentIndexTable(bc.bcId());
            bcIndicies.push_back(bc.globalIndex());
            bcMapping.insert(std::make_pair(bc.globalIndex(), storedJBCsTable.lastIndex()));
          }

          storedJCollisionsTable(collision.posZ(), collision.eventSel(), collision.alias_raw());
          storedJCollisionsParentIndexTable(collision.collisionId());

          auto JMcCollisionIndex = mcCollisionMapping.find(mcCollision.globalIndex());
          if (JMcCollisionIndex != mcCollisionMapping.end()) {
            storedJMcCollisionsLabelTable(JMcCollisionIndex->second);
          }
          int32_t storedBCID = -1;
          auto JBCIndex = bcMapping.find(collision.bcId());
          if (JBCIndex != bcMapping.end()) {
            storedBCID = JBCIndex->second;
          }
          storedJCollisionsBunchCrossingIndexTable(storedBCID);
          storedJChargedTriggerSelsTable(collision.chargedTriggerSel());
          storedJFullTriggerSelsTable(collision.fullTriggerSel());

          const auto tracksPerCollision = tracks.sliceBy(TracksPerCollision, collision.globalIndex());
          for (const auto& track : tracksPerCollision) {
            storedJTracksTable(storedJCollisionsTable.lastIndex(), track.pt(), track.eta(), track.phi(), track.energy(), track.sign(), track.trackSel());
            storedJTracksParentIndexTable(track.trackId());

            if (track.has_mcParticle()) {
              auto JParticleIndex = paticleMapping.find(track.mcParticleId());
              if (JParticleIndex != paticleMapping.end()) {
                storedJMcTracksLabelTable(JParticleIndex->second);
              } else {
                storedJMcTracksLabelTable(-1); // this can happen because there are some tracks that are reconstucted in a wrong collision, but their original McCollision did not pass the required cuts so that McParticle is not saved. These are very few but we should look into them further and see what to do about them
              }
            } else {
              storedJMcTracksLabelTable(-1);
            }
            trackMapping.insert(std::make_pair(track.globalIndex(), storedJTracksTable.lastIndex()));
          }

          const auto clustersPerCollision = clusters.sliceBy(ClustersPerCollision, collision.globalIndex());
          for (const auto& cluster : clustersPerCollision) {
            storedJClustersTable(storedJCollisionsTable.lastIndex(), cluster.id(), cluster.energy(), cluster.coreEnergy(), cluster.rawEnergy(),
                                 cluster.eta(), cluster.phi(), cluster.m02(), cluster.m20(), cluster.nCells(), cluster.time(), cluster.isExotic(), cluster.distanceToBadChannel(),
                                 cluster.nlm(), cluster.definition(), cluster.leadingCellEnergy(), cluster.subleadingCellEnergy(), cluster.leadingCellNumber(), cluster.subleadingCellNumber());
            storedJClustersParentIndexTable(cluster.clusterId());

            std::vector<int> clusterStoredJTrackIDs;
            for (const auto& clusterTrack : cluster.matchedTracks_as<soa::Join<aod::JTracks, aod::JTrackPIs>>()) {
              auto JtrackIndex = trackMapping.find(clusterTrack.globalIndex());
              if (JtrackIndex != trackMapping.end()) {
                clusterStoredJTrackIDs.push_back(JtrackIndex->second);
              }
            }
            storedJClustersMatchedTracksTable(clusterStoredJTrackIDs);
          }
        }
      }
    }
  }
  // process switch for output writing must be last
  // to run after all jet selections
  PROCESS_SWITCH(JetDerivedDataWriter, processMC, "write out data output tables for mc", false);

  void processMCP(soa::Join<aod::JMcCollisions, aod::JMcCollisionPIs> const& mcCollisions, soa::Join<aod::JMcParticles, aod::JMcParticlePIs> const& particles)
  {

    int particleTableIndex = 0;
    for (auto mcCollision : mcCollisions) {
      if (McCollisionFlag[mcCollision.globalIndex()]) { // you can also check if any of its detector level counterparts are correct
        std::map<int32_t, int32_t> paticleMapping;

        storedJMcCollisionsTable(mcCollision.posZ(), mcCollision.weight());
        storedJMcCollisionsParentIndexTable(mcCollision.mcCollisionId());

        const auto particlesPerMcCollision = particles.sliceBy(ParticlesPerMcCollision, mcCollision.globalIndex());

        for (auto particle : particlesPerMcCollision) {
          paticleMapping.insert(std::make_pair(particle.globalIndex(), particleTableIndex));
          particleTableIndex++;
        }
        for (auto particle : particlesPerMcCollision) {

          std::vector<int> mothersId;
          int daughtersId[2];
          if (particle.has_mothers()) {
            for (auto const& mother : particle.template mothers_as<soa::Join<aod::JMcParticles, aod::JMcParticlePIs>>()) {

              auto JMotherIndex = paticleMapping.find(mother.globalIndex());
              if (JMotherIndex != paticleMapping.end()) {
                mothersId.push_back(JMotherIndex->second);
              }
            }
          }
          auto i = 0;
          if (particle.has_daughters()) {
            for (auto const& daughter : particle.template daughters_as<soa::Join<aod::JMcParticles, aod::JMcParticlePIs>>()) {
              if (i > 1) {
                break;
              }
              auto JDaughterIndex = paticleMapping.find(daughter.globalIndex());
              if (JDaughterIndex != paticleMapping.end()) {
                daughtersId[i] = JDaughterIndex->second;
              }
              i++;
            }
          }
          storedJMcParticlesTable(storedJMcCollisionsTable.lastIndex(), particle.pt(), particle.eta(), particle.phi(), particle.y(), particle.e(), particle.pdgCode(), particle.getGenStatusCode(), particle.getHepMCStatusCode(), particle.isPhysicalPrimary(), mothersId, daughtersId);
          storedJParticlesParentIndexTable(particle.mcParticleId());
        }
      }
    }
  }
  // process switch for output writing must be last
  // to run after all jet selections
  PROCESS_SWITCH(JetDerivedDataWriter, processMCP, "write out data output tables for mcp", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetDerivedDataWriter>(cfgc, TaskName{"jet-deriveddata-writer"}));

  return WorkflowSpec{tasks};
}
