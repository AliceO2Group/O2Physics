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

#include <MathUtils/Utils.h>
#include <algorithm>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"

#include "PWGJE/Core/JetHFUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetDerivedDataWriter {

  Configurable<float> chargedJetPtMin{"chargedJetPtMin", 0.0, "Minimum charged jet pt to accept event"};
  Configurable<float> chargedEventWiseSubtractedJetPtMin{"chargedEventWiseSubtractedJetPtMin", 0.0, "Minimum charged event-wise subtracted jet pt to accept event"};
  Configurable<float> chargedMCPJetPtMin{"chargedMCPJetPtMin", 0.0, "Minimum charged mcp jet pt to accept event"};
  Configurable<float> neutralJetPtMin{"neutralJetPtMin", 0.0, "Minimum charged jet pt to accept event"};
  Configurable<float> fullJetPtMin{"fullJetPtMin", 0.0, "Minimum full jet pt to accept event"};
  Configurable<float> chargedEventWiseSubtractedD0JetPtMin{"chargedEventWiseSubtractedD0JetPtMin", 0.0, "Minimum charged event-wise subtracted D0 jet pt to accept event"};
  Configurable<float> chargedD0JetPtMin{"chargedD0JetPtMin", 0.0, "Minimum charged D0 jet pt to accept event"};
  Configurable<float> chargedEventWiseSubtractedLcJetPtMin{"chargedEventWiseSubtractedLcJetPtMin", 0.0, "Minimum charged event-wise subtracted Lc jet pt to accept event"};
  Configurable<float> chargedLcJetPtMin{"chargedLcJetPtMin", 0.0, "Minimum charged Lc jet pt to accept event"};

  Configurable<bool> performTrackSelection{"performTrackSelection", true, "only save tracks that pass one of the track selections"};
  Configurable<bool> saveBCsTable{"saveBCsTable", true, "save the bunch crossing table to the output"};
  Configurable<bool> saveClustersTable{"saveClustersTable", true, "save the clusters table to the output"};
  Configurable<bool> saveD0Table{"saveD0Table", false, "save the D0 table to the output"};
  Configurable<bool> saveLcTable{"saveLcTable", false, "save the Lc table to the output"};

  Produces<aod::StoredCollisionCounts> storedCollisionCountsTable;
  Produces<aod::StoredJDummys> storedJDummysTable;
  Produces<aod::StoredJBCs> storedJBCsTable;
  Produces<aod::StoredJBCPIs> storedJBCParentIndexTable;
  Produces<aod::StoredJCollisions> storedJCollisionsTable;
  Produces<aod::StoredJCollisionPIs> storedJCollisionsParentIndexTable;
  Produces<aod::StoredJCollisionBCs> storedJCollisionsBunchCrossingIndexTable;
  Produces<aod::StoredJChTrigSels> storedJChargedTriggerSelsTable;
  Produces<aod::StoredJFullTrigSels> storedJFullTriggerSelsTable;
  Produces<aod::StoredJChHFTrigSels> storedJChargedHFTriggerSelsTable;
  Produces<aod::StoredJMcCollisionLbs> storedJMcCollisionsLabelTable;
  Produces<aod::StoredJMcCollisions> storedJMcCollisionsTable;
  Produces<aod::StoredJMcCollisionPIs> storedJMcCollisionsParentIndexTable;
  Produces<aod::StoredJTracks> storedJTracksTable;
  Produces<aod::StoredJTrackExtras> storedJTracksExtraTable;
  Produces<aod::StoredJTrackPIs> storedJTracksParentIndexTable;
  Produces<aod::StoredJMcTrackLbs> storedJMcTracksLabelTable;
  Produces<aod::StoredJMcParticles> storedJMcParticlesTable;
  Produces<aod::StoredJMcParticlePIs> storedJParticlesParentIndexTable;
  Produces<aod::StoredJClusters> storedJClustersTable;
  Produces<aod::StoredJClusterPIs> storedJClustersParentIndexTable;
  Produces<aod::StoredJClusterTracks> storedJClustersMatchedTracksTable;

  Produces<aod::StoredHfD0CollBases> storedD0CollisionsTable;
  Produces<aod::StoredHfD0Bases> storedD0sTable;
  Produces<aod::StoredHfD0Pars> storedD0ParsTable;
  Produces<aod::StoredHfD0ParEs> storedD0ParExtrasTable;
  Produces<aod::StoredHfD0Sels> storedD0SelsTable;
  Produces<aod::StoredHfD0Mls> storedD0MlsTable;
  Produces<aod::StoredHfD0Mcs> storedD0McsTable;
  Produces<aod::StoredJD0Ids> storedD0IdsTable;
  Produces<aod::StoredHfD0PBases> storedD0ParticlesTable;
  Produces<aod::StoredJD0PIds> storedD0ParticleIdsTable;

  Produces<aod::StoredHf3PCollBases> storedLcCollisionsTable;
  Produces<aod::StoredHf3PBases> storedLcsTable;
  Produces<aod::StoredHf3PPars> storedLcParsTable;
  Produces<aod::StoredHf3PParEs> storedLcParExtrasTable;
  Produces<aod::StoredHf3PSels> storedLcSelsTable;
  Produces<aod::StoredHf3PMls> storedLcMlsTable;
  Produces<aod::StoredHf3PMcs> storedLcMcsTable;
  Produces<aod::StoredJLcIds> storedLcIdsTable;
  Produces<aod::StoredHf3PPBases> storedLcParticlesTable;
  Produces<aod::StoredJLcPIds> storedLcParticleIdsTable;

  PresliceUnsorted<soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::JCollisionBCs, aod::JChTrigSels, aod::JFullTrigSels, aod::JChHFTrigSels, aod::JMcCollisionLbs>> CollisionsPerMcCollision = aod::jmccollisionlb::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::JMcParticles, aod::JMcParticlePIs>> ParticlesPerMcCollision = aod::jmcparticle::mcCollisionId;
  Preslice<soa::Join<aod::JTracks, aod::JTrackPIs, aod::JMcTrackLbs>> TracksPerCollision = aod::jtrack::collisionId;
  Preslice<soa::Join<aod::JClusters, aod::JClusterPIs, aod::JClusterTracks>> ClustersPerCollision = aod::jcluster::collisionId;
  Preslice<CandidatesD0MCD> D0sPerCollision = aod::jd0indices::collisionId;
  Preslice<CandidatesLcMCD> LcsPerCollision = aod::jlcindices::collisionId;

  std::vector<bool> collisionFlag;
  std::vector<bool> McCollisionFlag;
  std::vector<int> bcIndicies;

  uint32_t precisionPositionMask;
  uint32_t precisionMomentumMask;

  void init(InitContext&)
  {
    precisionPositionMask = 0xFFFFFC00; // 13 bits
    precisionMomentumMask = 0xFFFFFC00; // 13 bits  this is currently keept at 13 bits wihich gives roughly a resolution of 1/8000. This can be increased to 15 bits if really needed
  }

  bool acceptCollision(aod::JCollision const&)
  {
    return true;
  }

  void processCollisions(aod::JCollisions const& collisions)
  {
    collisionFlag.resize(collisions.size());
    std::fill(collisionFlag.begin(), collisionFlag.end(), false);
  }

  void processMcCollisions(aod::JMcCollisions const& Mccollisions)
  {
    McCollisionFlag.resize(Mccollisions.size());
    std::fill(McCollisionFlag.begin(), McCollisionFlag.end(), false);
  }

  template <typename T>
  void processJets(T& jets)
  {
    float jetPtMin = 0.0;
    if constexpr (std::is_same_v<std::decay_t<T>, aod::ChargedJets> || std::is_same_v<std::decay_t<T>, aod::ChargedMCDetectorLevelJets>) {
      jetPtMin = chargedJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::ChargedEventWiseSubtractedJets>) {
      jetPtMin = chargedEventWiseSubtractedJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::ChargedMCParticleLevelJets>) {
      jetPtMin = chargedMCPJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::NeutralJets>) {
      jetPtMin = neutralJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::FullJets>) {
      jetPtMin = fullJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::D0ChargedJets>) {
      jetPtMin = chargedD0JetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::D0ChargedEventWiseSubtractedJets>) {
      jetPtMin = chargedEventWiseSubtractedD0JetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::LcChargedJets>) {
      jetPtMin = chargedLcJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::LcChargedEventWiseSubtractedJets>) {
      jetPtMin = chargedEventWiseSubtractedLcJetPtMin;
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
  PROCESS_SWITCH(JetDerivedDataWriter, processCollisions, "setup the writing for data and MCD", true);
  PROCESS_SWITCH(JetDerivedDataWriter, processMcCollisions, "setup the writing for MCP", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::ChargedJets>, processChargedJets, "process charged jets", true);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::ChargedEventWiseSubtractedJets>, processChargedEventWiseSubtractedJets, "process charged event-wise subtracted jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::ChargedMCDetectorLevelJets>, processChargedMCDJets, "process charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::ChargedMCParticleLevelJets>, processChargedMCPJets, "process charged mcp jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::NeutralJets>, processNeutralJets, "process neutral jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::FullJets>, processFullJets, "process full jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::D0ChargedJets>, processD0ChargedJets, "process D0 charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::D0ChargedEventWiseSubtractedJets>, processD0ChargedEventWiseSubtractedJets, "process D0 event-wise subtracted charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::LcChargedJets>, processLcChargedJets, "process Lc charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::LcChargedEventWiseSubtractedJets>, processLcChargedEventWiseSubtractedJets, "process Lc event-wise subtracted charged jets", false);

  void processDummyTable(aod::JDummys const&)
  {
    storedJDummysTable(1);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDummyTable, "write out dummy output table", true);

  void processCollisionCounting(aod::JCollisions const& collisions, aod::CollisionCounts const& collisionCounts)
  {
    int readCollisionCounter = 0;
    int writtenCollisionCounter = 0;
    for (const auto& collision : collisions) {
      readCollisionCounter++;
      if (collisionFlag[collision.globalIndex()]) {
        writtenCollisionCounter++;
      }
    }
    std::vector<int> previousReadCounts;
    std::vector<int> previousWrittenCounts;
    int iPreviousDataFrame = 0;
    for (const auto& collisionCount : collisionCounts) {
      auto readCollisionCounterSpan = collisionCount.readCounts();
      auto writtenCollisionCounterSpan = collisionCount.writtenCounts();
      if (iPreviousDataFrame == 0) {
        std::copy(readCollisionCounterSpan.begin(), readCollisionCounterSpan.end(), std::back_inserter(previousReadCounts));
        std::copy(writtenCollisionCounterSpan.begin(), writtenCollisionCounterSpan.end(), std::back_inserter(previousWrittenCounts));
      } else {
        for (unsigned int i = 0; i < previousReadCounts.size(); i++) {
          previousReadCounts[i] += readCollisionCounterSpan[i];
          previousWrittenCounts[i] += writtenCollisionCounterSpan[i];
        }
      }
      iPreviousDataFrame++;
    }
    previousReadCounts.push_back(readCollisionCounter);
    previousWrittenCounts.push_back(writtenCollisionCounter);
    storedCollisionCountsTable(previousReadCounts, previousWrittenCounts);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processCollisionCounting, "write out collision counting output table", false);

  void processData(soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::JCollisionBCs, aod::JChTrigSels, aod::JFullTrigSels, aod::JChHFTrigSels>::iterator const& collision, soa::Join<aod::JBCs, aod::JBCPIs> const&, soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs> const& tracks, soa::Join<aod::JClusters, aod::JClusterPIs, aod::JClusterTracks> const& clusters, aod::HfD0CollBases const&, CandidatesD0Data const& D0s, aod::Hf3PCollBases const&, CandidatesLcData const& Lcs)
  {
    std::map<int32_t, int32_t> bcMapping;
    std::map<int32_t, int32_t> trackMapping;

    if (collisionFlag[collision.globalIndex()]) {
      if (saveBCsTable) {
        auto bc = collision.bc_as<soa::Join<aod::JBCs, aod::JBCPIs>>();
        if (std::find(bcIndicies.begin(), bcIndicies.end(), bc.globalIndex()) == bcIndicies.end()) {
          storedJBCsTable(bc.runNumber(), bc.globalBC(), bc.timestamp());
          storedJBCParentIndexTable(bc.bcId());
          bcIndicies.push_back(bc.globalIndex());
          bcMapping.insert(std::make_pair(bc.globalIndex(), storedJBCsTable.lastIndex()));
        }
      }

      storedJCollisionsTable(collision.posX(), collision.posY(), collision.posZ(), collision.multiplicity(), collision.centrality(), collision.eventSel(), collision.alias_raw());
      storedJCollisionsParentIndexTable(collision.collisionId());
      if (saveBCsTable) {
        int32_t storedBCID = -1;
        auto JBCIndex = bcMapping.find(collision.bcId());
        if (JBCIndex != bcMapping.end()) {
          storedBCID = JBCIndex->second;
        }
        storedJCollisionsBunchCrossingIndexTable(storedBCID);
      }
      storedJChargedTriggerSelsTable(collision.chargedTriggerSel());
      storedJFullTriggerSelsTable(collision.fullTriggerSel());
      storedJChargedHFTriggerSelsTable(collision.chargedHFTriggerSel());

      for (const auto& track : tracks) {
        if (performTrackSelection && !(track.trackSel() & ~(1 << jetderiveddatautilities::JTrackSel::trackSign))) { // skips tracks that pass no selections. This might cause a problem with tracks matched with clusters. We should generate a track selection purely for cluster matched tracks so that they are kept
          continue;
        }
        storedJTracksTable(storedJCollisionsTable.lastIndex(), o2::math_utils::detail::truncateFloatFraction(track.pt(), precisionMomentumMask), o2::math_utils::detail::truncateFloatFraction(track.eta(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.phi(), precisionPositionMask), track.trackSel());
        storedJTracksExtraTable(o2::math_utils::detail::truncateFloatFraction(track.dcaXY(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.dcaZ(), precisionPositionMask));
        storedJTracksParentIndexTable(track.trackId());
        trackMapping.insert(std::make_pair(track.globalIndex(), storedJTracksTable.lastIndex()));
      }
      if (saveClustersTable) {
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

      if (saveD0Table) {
        int nD0InCollision = 0;
        int32_t collisionD0Index = -1;
        for (const auto& D0 : D0s) {
          if (nD0InCollision == 0) {
            jethfutilities::fillD0CollisionTable(D0.hfD0CollBase_as<aod::HfD0CollBases>(), storedD0CollisionsTable, collisionD0Index);
          }
          nD0InCollision++;

          int32_t D0Index = -1;
          jethfutilities::fillD0CandidateTable<false>(D0, collisionD0Index, storedD0sTable, storedD0ParsTable, storedD0ParExtrasTable, storedD0SelsTable, storedD0MlsTable, storedD0McsTable, D0Index);

          int32_t prong0Id = -1;
          int32_t prong1Id = -1;
          auto JtrackIndex = trackMapping.find(D0.prong0Id());
          if (JtrackIndex != trackMapping.end()) {
            prong0Id = JtrackIndex->second;
          }
          JtrackIndex = trackMapping.find(D0.prong1Id());
          if (JtrackIndex != trackMapping.end()) {
            prong1Id = JtrackIndex->second;
          }
          storedD0IdsTable(storedJCollisionsTable.lastIndex(), prong0Id, prong1Id);
        }
      }

      if (saveLcTable) {
        int nLcInCollision = 0;
        int32_t collisionLcIndex = -1;
        for (const auto& Lc : Lcs) {
          if (nLcInCollision == 0) {
            jethfutilities::fillLcCollisionTable(Lc.hf3PCollBase_as<aod::Hf3PCollBases>(), storedLcCollisionsTable, collisionLcIndex);
          }
          nLcInCollision++;

          int32_t LcIndex = -1;
          jethfutilities::fillLcCandidateTable<false>(Lc, collisionLcIndex, storedLcsTable, storedLcParsTable, storedLcParExtrasTable, storedLcSelsTable, storedLcMlsTable, storedLcMcsTable, LcIndex);

          int32_t prong0Id = -1;
          int32_t prong1Id = -1;
          int32_t prong2Id = -1;
          auto JtrackIndex = trackMapping.find(Lc.prong0Id());
          if (JtrackIndex != trackMapping.end()) {
            prong0Id = JtrackIndex->second;
          }
          JtrackIndex = trackMapping.find(Lc.prong1Id());
          if (JtrackIndex != trackMapping.end()) {
            prong1Id = JtrackIndex->second;
          }
          JtrackIndex = trackMapping.find(Lc.prong2Id());
          if (JtrackIndex != trackMapping.end()) {
            prong2Id = JtrackIndex->second;
          }
          storedLcIdsTable(storedJCollisionsTable.lastIndex(), prong0Id, prong1Id, prong2Id);
        }
      }
    }
  }
  // process switch for output writing must be last
  // to run after all jet selections
  PROCESS_SWITCH(JetDerivedDataWriter, processData, "write out data output tables", false);

  void processMC(soa::Join<aod::JMcCollisions, aod::JMcCollisionPIs> const& mcCollisions, soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::JCollisionBCs, aod::JChTrigSels, aod::JFullTrigSels, aod::JChHFTrigSels, aod::JMcCollisionLbs> const& collisions, soa::Join<aod::JBCs, aod::JBCPIs> const&, soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs, aod::JMcTrackLbs> const& tracks, soa::Join<aod::JClusters, aod::JClusterPIs, aod::JClusterTracks> const& clusters, soa::Join<aod::JMcParticles, aod::JMcParticlePIs> const& particles, aod::HfD0CollBases const&, CandidatesD0MCD const& D0s, CandidatesD0MCP const& D0Particles, aod::Hf3PCollBases const&, CandidatesLcMCD const& Lcs, CandidatesLcMCP const& LcParticles)
  {
    std::map<int32_t, int32_t> bcMapping;
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

        storedJMcCollisionsTable(mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(), mcCollision.weight());
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
          storedJMcParticlesTable(storedJMcCollisionsTable.lastIndex(), o2::math_utils::detail::truncateFloatFraction(particle.pt(), precisionMomentumMask), o2::math_utils::detail::truncateFloatFraction(particle.eta(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(particle.phi(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(particle.y(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(particle.e(), precisionMomentumMask), particle.pdgCode(), particle.getGenStatusCode(), particle.getHepMCStatusCode(), particle.isPhysicalPrimary(), mothersId, daughtersId);
          storedJParticlesParentIndexTable(particle.mcParticleId());
        }

        if (saveD0Table) {
          for (const auto& D0Particle : D0Particles) {
            int32_t D0ParticleIndex = -1;
            jethfutilities::fillD0CandidateMcTable(D0Particle, storedD0ParticlesTable, D0ParticleIndex);
            int32_t d0ParticleId = -1;
            auto JParticleIndex = paticleMapping.find(D0Particle.mcParticleId());
            if (JParticleIndex != paticleMapping.end()) {
              d0ParticleId = JParticleIndex->second;
            }
            storedD0ParticleIdsTable(storedJMcCollisionsTable.lastIndex(), d0ParticleId);
          }
        }

        if (saveLcTable) {
          for (const auto& LcParticle : LcParticles) {
            int32_t LcParticleIndex = -1;
            jethfutilities::fillLcCandidateMcTable(LcParticle, storedLcParticlesTable, LcParticleIndex);
            int32_t LcParticleId = -1;
            auto JParticleIndex = paticleMapping.find(LcParticle.mcParticleId());
            if (JParticleIndex != paticleMapping.end()) {
              LcParticleId = JParticleIndex->second;
            }
            storedLcParticleIdsTable(storedJMcCollisionsTable.lastIndex(), LcParticleId);
          }
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
          if (saveBCsTable) {
            auto bc = collision.bc_as<soa::Join<aod::JBCs, aod::JBCPIs>>();
            if (std::find(bcIndicies.begin(), bcIndicies.end(), bc.globalIndex()) == bcIndicies.end()) {
              storedJBCsTable(bc.runNumber(), bc.globalBC(), bc.timestamp());
              storedJBCParentIndexTable(bc.bcId());
              bcIndicies.push_back(bc.globalIndex());
              bcMapping.insert(std::make_pair(bc.globalIndex(), storedJBCsTable.lastIndex()));
            }
          }

          storedJCollisionsTable(collision.posX(), collision.posY(), collision.posZ(), collision.multiplicity(), collision.centrality(), collision.eventSel(), collision.alias_raw());
          storedJCollisionsParentIndexTable(collision.collisionId());

          auto JMcCollisionIndex = mcCollisionMapping.find(mcCollision.globalIndex());
          if (JMcCollisionIndex != mcCollisionMapping.end()) {
            storedJMcCollisionsLabelTable(JMcCollisionIndex->second);
          }
          if (saveBCsTable) {
            int32_t storedBCID = -1;
            auto JBCIndex = bcMapping.find(collision.bcId());
            if (JBCIndex != bcMapping.end()) {
              storedBCID = JBCIndex->second;
            }
            storedJCollisionsBunchCrossingIndexTable(storedBCID);
          }
          storedJChargedTriggerSelsTable(collision.chargedTriggerSel());
          storedJFullTriggerSelsTable(collision.fullTriggerSel());
          storedJChargedHFTriggerSelsTable(collision.chargedHFTriggerSel());

          const auto tracksPerCollision = tracks.sliceBy(TracksPerCollision, collision.globalIndex());
          for (const auto& track : tracksPerCollision) {
            if (performTrackSelection && !(track.trackSel() & ~(1 << jetderiveddatautilities::JTrackSel::trackSign))) { // skips tracks that pass no selections. This might cause a problem with tracks matched with clusters. We should generate a track selection purely for cluster matched tracks so that they are kept
              continue;
            }
            storedJTracksTable(storedJCollisionsTable.lastIndex(), o2::math_utils::detail::truncateFloatFraction(track.pt(), precisionMomentumMask), o2::math_utils::detail::truncateFloatFraction(track.eta(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.phi(), precisionPositionMask), track.trackSel());
            storedJTracksExtraTable(o2::math_utils::detail::truncateFloatFraction(track.dcaXY(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.dcaZ(), precisionPositionMask));
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
          if (saveClustersTable) {
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

          if (saveD0Table) {
            const auto d0sPerCollision = D0s.sliceBy(D0sPerCollision, collision.globalIndex());
            int nD0InCollision = 0;
            int32_t collisionD0Index = -1;
            for (const auto& D0 : d0sPerCollision) {
              if (nD0InCollision == 0) {
                jethfutilities::fillD0CollisionTable(D0.hfD0CollBase_as<aod::HfD0CollBases>(), storedD0CollisionsTable, collisionD0Index);
              }
              nD0InCollision++;

              int32_t D0Index = -1;
              jethfutilities::fillD0CandidateTable<true>(D0, collisionD0Index, storedD0sTable, storedD0ParsTable, storedD0ParExtrasTable, storedD0SelsTable, storedD0MlsTable, storedD0McsTable, D0Index);

              int32_t prong0Id = -1;
              int32_t prong1Id = -1;
              auto JtrackIndex = trackMapping.find(D0.prong0Id());
              if (JtrackIndex != trackMapping.end()) {
                prong0Id = JtrackIndex->second;
              }
              JtrackIndex = trackMapping.find(D0.prong1Id());
              if (JtrackIndex != trackMapping.end()) {
                prong1Id = JtrackIndex->second;
              }
              storedD0IdsTable(storedJCollisionsTable.lastIndex(), prong0Id, prong1Id);
            }
          }

          if (saveLcTable) {
            const auto lcsPerCollision = Lcs.sliceBy(LcsPerCollision, collision.globalIndex());
            int nLcInCollision = 0;
            int32_t collisionLcIndex = -1;
            for (const auto& Lc : lcsPerCollision) {
              if (nLcInCollision == 0) {
                jethfutilities::fillLcCollisionTable(Lc.hf3PCollBase_as<aod::Hf3PCollBases>(), storedLcCollisionsTable, collisionLcIndex);
              }
              nLcInCollision++;

              int32_t LcIndex = -1;
              jethfutilities::fillLcCandidateTable<true>(Lc, collisionLcIndex, storedLcsTable, storedLcParsTable, storedLcParExtrasTable, storedLcSelsTable, storedLcMlsTable, storedLcMcsTable, LcIndex);

              int32_t prong0Id = -1;
              int32_t prong1Id = -1;
              int32_t prong2Id = -1;
              auto JtrackIndex = trackMapping.find(Lc.prong0Id());
              if (JtrackIndex != trackMapping.end()) {
                prong0Id = JtrackIndex->second;
              }
              JtrackIndex = trackMapping.find(Lc.prong1Id());
              if (JtrackIndex != trackMapping.end()) {
                prong1Id = JtrackIndex->second;
              }
              JtrackIndex = trackMapping.find(Lc.prong2Id());
              if (JtrackIndex != trackMapping.end()) {
                prong2Id = JtrackIndex->second;
              }
              storedLcIdsTable(storedJCollisionsTable.lastIndex(), prong0Id, prong1Id, prong2Id);
            }
          }
        }
      }
    }
  }
  // process switch for output writing must be last
  // to run after all jet selections
  PROCESS_SWITCH(JetDerivedDataWriter, processMC, "write out data output tables for mc", false);

  void processMCP(soa::Join<aod::JMcCollisions, aod::JMcCollisionPIs> const& mcCollisions, soa::Join<aod::JMcParticles, aod::JMcParticlePIs> const& particles, CandidatesD0MCP const& D0Particles, CandidatesLcMCP const& LcParticles)
  {

    int particleTableIndex = 0;
    for (auto mcCollision : mcCollisions) {
      if (McCollisionFlag[mcCollision.globalIndex()]) { // you can also check if any of its detector level counterparts are correct
        std::map<int32_t, int32_t> paticleMapping;

        storedJMcCollisionsTable(mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(), mcCollision.weight());
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
          storedJMcParticlesTable(storedJMcCollisionsTable.lastIndex(), o2::math_utils::detail::truncateFloatFraction(particle.pt(), precisionMomentumMask), o2::math_utils::detail::truncateFloatFraction(particle.eta(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(particle.phi(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(particle.y(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(particle.e(), precisionMomentumMask), particle.pdgCode(), particle.getGenStatusCode(), particle.getHepMCStatusCode(), particle.isPhysicalPrimary(), mothersId, daughtersId);
          storedJParticlesParentIndexTable(particle.mcParticleId());
        }

        if (saveD0Table) {
          for (const auto& D0Particle : D0Particles) {
            int32_t D0ParticleIndex = -1;
            jethfutilities::fillD0CandidateMcTable(D0Particle, storedD0ParticlesTable, D0ParticleIndex);
            int32_t d0ParticleId = -1;
            auto JParticleIndex = paticleMapping.find(D0Particle.mcParticleId());
            if (JParticleIndex != paticleMapping.end()) {
              d0ParticleId = JParticleIndex->second;
            }
            storedD0ParticleIdsTable(storedJMcCollisionsTable.lastIndex(), d0ParticleId);
          }
        }
        if (saveLcTable) {
          for (const auto& LcParticle : LcParticles) {
            int32_t LcParticleIndex = -1;
            jethfutilities::fillLcCandidateMcTable(LcParticle, storedLcParticlesTable, LcParticleIndex);
            int32_t LcParticleId = -1;
            auto JParticleIndex = paticleMapping.find(LcParticle.mcParticleId());
            if (JParticleIndex != paticleMapping.end()) {
              LcParticleId = JParticleIndex->second;
            }
            storedLcParticleIdsTable(storedJMcCollisionsTable.lastIndex(), LcParticleId);
          }
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
