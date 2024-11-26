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
/// \brief Task to skim jet framework tables (aod::JetCollisions, aod::JetTracks, aod::JetClusters, ...)
/// while adjusting indices accordingly
///
/// \author Jochen Klein <jochen.klein@cern.ch>
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <MathUtils/Utils.h>
#include <algorithm>

#include <TRandom3.h>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"

#include "PWGJE/Core/JetHFUtilities.h"
#include "PWGJE/Core/JetDQUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetDerivedDataWriter {

  struct : ConfigurableGroup {
    Configurable<float> thresholdChargedJetPtMin{"thresholdChargedJetPtMin", 0.0, "Minimum charged jet pt to accept event"};
    Configurable<float> thresholdChargedEventWiseSubtractedJetPtMin{"thresholdChargedEventWiseSubtractedJetPtMin", 0.0, "Minimum charged event-wise subtracted jet pt to accept event"};
    Configurable<float> thresholdChargedMCPJetPtMin{"thresholdChargedMCPJetPtMin", 0.0, "Minimum charged mcp jet pt to accept event"};
    Configurable<float> thresholdNeutralJetPtMin{"thresholdNeutralJetPtMin", 0.0, "Minimum neutral jet pt to accept event"};
    Configurable<float> thresholdNeutralMCPJetPtMin{"thresholdNeutralMCPJetPtMin", 0.0, "Minimum neutal mcp jet pt to accept event"};
    Configurable<float> thresholdFullJetPtMin{"thresholdFullJetPtMin", 0.0, "Minimum full jet pt to accept event"};
    Configurable<float> thresholdFullMCPJetPtMin{"thresholdFullMCPJetPtMin", 0.0, "Minimum full mcp jet pt to accept event"};
    Configurable<float> thresholdChargedD0JetPtMin{"thresholdChargedD0JetPtMin", 0.0, "Minimum charged D0 jet pt to accept event"};
    Configurable<float> thresholdChargedEventWiseSubtractedD0JetPtMin{"thresholdChargedEventWiseSubtractedD0JetPtMin", 0.0, "Minimum charged event-wise subtracted D0 jet pt to accept event"};
    Configurable<float> thresholdChargedD0MCPJetPtMin{"thresholdChargedD0MCPJetPtMin", 0.0, "Minimum charged D0 mcp jet pt to accept event"};
    Configurable<float> thresholdChargedLcJetPtMin{"thresholdChargedLcJetPtMin", 0.0, "Minimum charged Lc jet pt to accept event"};
    Configurable<float> thresholdChargedEventWiseSubtractedLcJetPtMin{"thresholdChargedEventWiseSubtractedLcJetPtMin", 0.0, "Minimum charged event-wise subtracted Lc jet pt to accept event"};
    Configurable<float> thresholdChargedLcMCPJetPtMin{"thresholdChargedLcMCPJetPtMin", 0.0, "Minimum charged Lc mcp jet pt to accept event"};
    Configurable<float> thresholdChargedDielectronJetPtMin{"thresholdChargedDielectronJetPtMin", 0.0, "Minimum charged Dielectron jet pt to accept event"};
    Configurable<float> thresholdChargedEventWiseSubtractedDielectronJetPtMin{"thresholdChargedEventWiseSubtractedDielectronJetPtMin", 0.0, "Minimum charged event-wise subtracted Dielectron jet pt to accept event"};
    Configurable<float> thresholdChargedDielectronMCPJetPtMin{"thresholdChargedDielectronMCPJetPtMin", 0.0, "Minimum charged Dielectron mcp jet pt to accept event"};
    Configurable<float> thresholdTriggerTrackPtMin{"thresholdTriggerTrackPtMin", 0.0, "Minimum trigger track pt to accept event"};
    Configurable<float> thresholdClusterEnergyMin{"thresholdClusterEnergyMin", 0.0, "Minimum cluster energy to accept event"};
    Configurable<int> downscaleFactor{"downscaleFactor", 1, "random downscale of selected events"};

    Configurable<float> vertexZCut{"vertexZCut", 10.0, "z-vertex cut on event"};
    Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
    Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};
    Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
    Configurable<bool> performTrackSelection{"performTrackSelection", true, "only save tracks that pass one of the track selections"};
    Configurable<float> trackPtSelectionMin{"trackPtSelectionMin", 0.15, "only save tracks that have a pT larger than this pT"};
    Configurable<float> trackEtaSelectionMax{"trackEtaSelectionMax", 0.9, "only save tracks that have an eta smaller than this eta"};
    Configurable<bool> saveBCsTable{"saveBCsTable", true, "save the bunch crossing table to the output"};
    Configurable<bool> saveClustersTable{"saveClustersTable", false, "save the clusters table to the output"};
    Configurable<bool> saveD0Table{"saveD0Table", false, "save the D0 table to the output"};
    Configurable<bool> saveLcTable{"saveLcTable", false, "save the Lc table to the output"};
    Configurable<bool> saveDielectronTable{"saveDielectronTable", false, "save the Dielectron table to the output"};

    Configurable<std::string> triggerMasks{"triggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};
  } config;

  struct : ProducesGroup {
    Produces<aod::StoredJDummys> storedJDummysTable;
    Produces<aod::StoredJBCs> storedJBCsTable;
    Produces<aod::StoredJBCPIs> storedJBCParentIndexTable;
    Produces<aod::StoredJCollisions> storedJCollisionsTable;
    Produces<aod::StoredJCollisionPIs> storedJCollisionsParentIndexTable;
    Produces<aod::StoredJCollisionBCs> storedJCollisionsBunchCrossingIndexTable;
    Produces<aod::StoredJEMCCollisionLbs> storedJCollisionsEMCalLabelTable;
    Produces<aod::StoredJMcCollisionLbs> storedJMcCollisionsLabelTable;
    Produces<aod::StoredJMcCollisions> storedJMcCollisionsTable;
    Produces<aod::StoredJMcCollisionPIs> storedJMcCollisionsParentIndexTable;
    Produces<aod::StoredJTracks> storedJTracksTable;
    Produces<aod::StoredJTrackExtras> storedJTracksExtraTable;
    Produces<aod::StoredJEMCTracks> storedJTracksEMCalTable;
    Produces<aod::StoredJTrackPIs> storedJTracksParentIndexTable;
    Produces<aod::StoredJMcTrackLbs> storedJMcTracksLabelTable;
    Produces<aod::StoredJMcParticles> storedJMcParticlesTable;
    Produces<aod::StoredJMcParticlePIs> storedJParticlesParentIndexTable;
    Produces<aod::StoredJClusters> storedJClustersTable;
    Produces<aod::StoredJClusterPIs> storedJClustersParentIndexTable;
    Produces<aod::StoredJClusterTracks> storedJClustersMatchedTracksTable;
    Produces<aod::StoredJMcClusterLbs> storedJMcClustersLabelTable;

    Produces<aod::StoredHfD0CollBases> storedD0CollisionsTable;
    Produces<aod::StoredJD0CollisionIds> storedD0CollisionIdsTable;
    Produces<aod::StoredHfD0Bases> storedD0sTable;
    Produces<aod::StoredHfD0Pars> storedD0ParsTable;
    Produces<aod::StoredHfD0ParEs> storedD0ParExtrasTable;
    Produces<aod::StoredHfD0Sels> storedD0SelsTable;
    Produces<aod::StoredHfD0Mls> storedD0MlsTable;
    Produces<aod::StoredHfD0Mcs> storedD0McsTable;
    Produces<aod::StoredJD0Ids> storedD0IdsTable;
    Produces<aod::StoredHfD0McCollBases> storedD0McCollisionsTable;
    Produces<aod::StoredJD0McCollisionIds> storedD0McCollisionIdsTable;
    Produces<aod::StoredHfD0McRCollIds> storedD0McCollisionsMatchingTable;
    Produces<aod::StoredHfD0PBases> storedD0ParticlesTable;
    Produces<aod::StoredJD0PIds> storedD0ParticleIdsTable;

    Produces<aod::StoredHf3PCollBases> storedLcCollisionsTable;
    Produces<aod::StoredJLcCollisionIds> storedLcCollisionIdsTable;
    Produces<aod::StoredHf3PBases> storedLcsTable;
    Produces<aod::StoredHf3PPars> storedLcParsTable;
    Produces<aod::StoredHf3PParEs> storedLcParExtrasTable;
    Produces<aod::StoredHf3PSels> storedLcSelsTable;
    Produces<aod::StoredHf3PMls> storedLcMlsTable;
    Produces<aod::StoredHf3PMcs> storedLcMcsTable;
    Produces<aod::StoredJLcIds> storedLcIdsTable;
    Produces<aod::StoredHf3PMcCollBases> storedLcMcCollisionsTable;
    Produces<aod::StoredJLcMcCollisionIds> storedLcMcCollisionIdsTable;
    Produces<aod::StoredHf3PMcRCollIds> storedLcMcCollisionsMatchingTable;
    Produces<aod::StoredHf3PPBases> storedLcParticlesTable;
    Produces<aod::StoredJLcPIds> storedLcParticleIdsTable;

    Produces<aod::StoredReducedEvents> storedDielectronCollisionsTable;
    Produces<aod::StoredJDielectronCollisionIds> storedDielectronCollisionIdsTable;
    Produces<aod::StoredDielectrons> storedDielectronsTable;
    Produces<aod::StoredJDielectronIds> storedDielectronIdsTable;
    Produces<aod::StoredJDielectronMcCollisions> storedDielectronMcCollisionsTable;
    Produces<aod::StoredJDielectronMcCollisionIds> storedDielectronMcCollisionIdsTable;
    // Produces<aod::StoredHfD0McRCollIds> storedD0McCollisionsMatchingTable; //this doesnt exist for Dileptons yet
    Produces<aod::StoredJDielectronMcs> storedDielectronParticlesTable;
    Produces<aod::StoredJDielectronMcIds> storedDielectronParticleIdsTable;
  } products;

  PresliceUnsorted<soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::JCollisionBCs, aod::JChTrigSels, aod::JFullTrigSels, aod::JChHFTrigSels, aod::JMcCollisionLbs>> CollisionsPerMcCollision = aod::jmccollisionlb::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::JMcParticles, aod::JMcParticlePIs>> ParticlesPerMcCollision = aod::jmcparticle::mcCollisionId;
  Preslice<soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs, aod::JMcTrackLbs>> TracksPerCollision = aod::jtrack::collisionId;
  Preslice<soa::Join<aod::JClusters, aod::JClusterPIs, aod::JClusterTracks>> ClustersPerCollision = aod::jcluster::collisionId;
  Preslice<soa::Join<aod::McCollisionsD0, aod::HfD0McRCollIds>> D0McCollisionsPerMcCollision = aod::jd0indices::mcCollisionId;
  Preslice<soa::Join<aod::McCollisionsLc, aod::Hf3PMcRCollIds>> LcMcCollisionsPerMcCollision = aod::jlcindices::mcCollisionId;
  Preslice<aod::McCollisionsDielectron> DielectronMcCollisionsPerMcCollision = aod::jdielectronindices::mcCollisionId;
  Preslice<aod::CollisionsD0> D0CollisionsPerCollision = aod::jd0indices::collisionId;
  Preslice<aod::CollisionsLc> LcCollisionsPerCollision = aod::jlcindices::collisionId;
  Preslice<aod::CollisionsDielectron> DielectronCollisionsPerCollision = aod::jdielectronindices::collisionId;
  Preslice<aod::CandidatesD0MCD> D0sPerCollision = aod::jd0indices::collisionId;
  Preslice<aod::CandidatesLcMCD> LcsPerCollision = aod::jlcindices::collisionId;
  Preslice<aod::CandidatesDielectronMCD> DielectronsPerCollision = aod::jdielectronindices::collisionId;
  PresliceUnsorted<aod::JEMCTracks> EMCTrackPerTrack = aod::jemctrack::trackId;

  std::vector<bool> collisionFlag;
  std::vector<bool> McCollisionFlag;
  std::vector<int32_t> bcIndicies;

  uint32_t precisionPositionMask;
  uint32_t precisionMomentumMask;

  TRandom3 randomNumber;

  std::vector<int> triggerMaskBits;
  void init(InitContext&)
  {
    precisionPositionMask = 0xFFFFFC00; // 13 bits
    precisionMomentumMask = 0xFFFFFC00; // 13 bits  this is currently keept at 13 bits wihich gives roughly a resolution of 1/8000. This can be increased to 15 bits if really needed
    randomNumber.SetSeed(0);
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(config.triggerMasks);
  }

  bool acceptCollision(aod::JCollision const&)
  {
    return true;
  }

  void processSetupCollisions(aod::JCollisions const& collisions)
  {
    collisionFlag.clear();
    collisionFlag.resize(collisions.size());
    std::fill(collisionFlag.begin(), collisionFlag.end(), false);
  }

  void processSetupMcCollisions(aod::JMcCollisions const& McCollisions)
  {
    McCollisionFlag.clear();
    McCollisionFlag.resize(McCollisions.size());
    std::fill(McCollisionFlag.begin(), McCollisionFlag.end(), false);
  }

  void processSetupAllCollisionsWithDownscaling(aod::JCollisions const& collisions)
  {
    collisionFlag.clear();
    collisionFlag.resize(collisions.size());
    for (const auto& collision : collisions) {
      if (randomNumber.Integer(config.downscaleFactor) == 0) {
        collisionFlag[collision.globalIndex()] = true;
      } else {
        collisionFlag[collision.globalIndex()] = false;
      }
    }
  }

  void processSetupAllMcCollisionsWithDownscaling(aod::JMcCollisions const& McCollisions)
  {
    McCollisionFlag.clear();
    McCollisionFlag.resize(McCollisions.size());
    for (const auto& mcCollision : McCollisions) {
      if (randomNumber.Integer(config.downscaleFactor) == 0) {
        McCollisionFlag[mcCollision.globalIndex()] = true;
      } else {
        McCollisionFlag[mcCollision.globalIndex()] = false;
      }
    }
  }

  template <typename T>
  void processDoDownscaling(T const& collisions)
  {
    for (const auto& collision : collisions) {
      if constexpr (std::is_same_v<std::decay_t<T>, aod::JCollisions>) {
        if (collisionFlag[collision.globalIndex()] && randomNumber.Integer(config.downscaleFactor) != 0) {
          collisionFlag[collision.globalIndex()] = false;
        }
      }
      if constexpr (std::is_same_v<std::decay_t<T>, aod::JMcCollisions>) {
        if (McCollisionFlag[collision.globalIndex()] && randomNumber.Integer(config.downscaleFactor) != 0) {
          McCollisionFlag[collision.globalIndex()] = false;
        }
      }
    }
  }

  void processSetupEventTriggering(aod::JCollision const& collision)
  {
    if (jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      collisionFlag[collision.globalIndex()] = true;
    }
  }

  void processDoCollisionSelections(aod::JCollision const& collision)
  { // can also add event selection like sel8 but goes a little against the derived data idea
    if (collision.centrality() < config.centralityMin || collision.centrality() >= config.centralityMax || collision.trackOccupancyInTimeRange() > config.trackOccupancyInTimeRangeMax || std::abs(collision.posZ()) > config.vertexZCut) {
      collisionFlag[collision.globalIndex()] = false;
    }
  }

  template <typename T>
  void processSelectionObjects(T& selectionObjects)
  {
    float selectionObjectPtMin = 0.0;
    if constexpr (std::is_same_v<std::decay_t<T>, aod::ChargedJets> || std::is_same_v<std::decay_t<T>, aod::ChargedMCDetectorLevelJets>) {
      selectionObjectPtMin = config.thresholdChargedJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::ChargedEventWiseSubtractedJets> || std::is_same_v<std::decay_t<T>, aod::ChargedMCDetectorLevelEventWiseSubtractedJets>) {
      selectionObjectPtMin = config.thresholdChargedEventWiseSubtractedJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::ChargedMCParticleLevelJets>) {
      selectionObjectPtMin = config.thresholdChargedMCPJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::NeutralJets> || std::is_same_v<std::decay_t<T>, aod::NeutralMCDetectorLevelJets>) {
      selectionObjectPtMin = config.thresholdNeutralJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::NeutralMCParticleLevelJets>) {
      selectionObjectPtMin = config.thresholdNeutralMCPJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::FullJets> || std::is_same_v<std::decay_t<T>, aod::FullMCDetectorLevelJets>) {
      selectionObjectPtMin = config.thresholdFullJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::FullMCParticleLevelJets>) {
      selectionObjectPtMin = config.thresholdFullMCPJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::D0ChargedJets> || std::is_same_v<std::decay_t<T>, aod::D0ChargedMCDetectorLevelJets>) {
      selectionObjectPtMin = config.thresholdChargedD0JetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::D0ChargedEventWiseSubtractedJets> || std::is_same_v<std::decay_t<T>, aod::D0ChargedMCDetectorLevelEventWiseSubtractedJets>) {
      selectionObjectPtMin = config.thresholdChargedEventWiseSubtractedD0JetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::D0ChargedMCParticleLevelJets>) {
      selectionObjectPtMin = config.thresholdChargedD0MCPJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::LcChargedJets> || std::is_same_v<std::decay_t<T>, aod::LcChargedMCDetectorLevelJets>) {
      selectionObjectPtMin = config.thresholdChargedLcJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::LcChargedEventWiseSubtractedJets> || std::is_same_v<std::decay_t<T>, aod::LcChargedMCDetectorLevelEventWiseSubtractedJets>) {
      selectionObjectPtMin = config.thresholdChargedEventWiseSubtractedLcJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::LcChargedMCParticleLevelJets>) {
      selectionObjectPtMin = config.thresholdChargedLcMCPJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::DielectronChargedJets> || std::is_same_v<std::decay_t<T>, aod::DielectronChargedMCDetectorLevelJets>) {
      selectionObjectPtMin = config.thresholdChargedDielectronJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::DielectronChargedEventWiseSubtractedJets> || std::is_same_v<std::decay_t<T>, aod::DielectronChargedMCDetectorLevelEventWiseSubtractedJets>) {
      selectionObjectPtMin = config.thresholdChargedEventWiseSubtractedDielectronJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::DielectronChargedMCParticleLevelJets>) {
      selectionObjectPtMin = config.thresholdChargedDielectronMCPJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::JTracks>) {
      selectionObjectPtMin = config.thresholdTriggerTrackPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::JClusters>) {
      selectionObjectPtMin = config.thresholdClusterEnergyMin;
    } else {
      selectionObjectPtMin = 0.0;
    }
    for (const auto& selectionObject : selectionObjects) {
      bool isTriggerObject = false;
      if constexpr (std::is_same_v<std::decay_t<T>, aod::JClusters>) {
        if (selectionObject.energy() >= selectionObjectPtMin) {
          isTriggerObject = true;
        }
      } else {
        if constexpr (std::is_same_v<std::decay_t<T>, aod::JTracks>) {
          if (config.performTrackSelection && !(selectionObject.trackSel() & ~(1 << jetderiveddatautilities::JTrackSel::trackSign))) {
            continue;
          }
        }
        if (selectionObject.pt() >= selectionObjectPtMin) {
          isTriggerObject = true;
        }
      }
      if (isTriggerObject) {
        if constexpr (std::is_same_v<std::decay_t<T>, aod::ChargedMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::NeutralMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::FullMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::D0ChargedMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::LcChargedMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::BplusChargedMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::DielectronChargedMCParticleLevelJets>) {
          if (selectionObject.mcCollisionId() >= 0) {
            McCollisionFlag[selectionObject.mcCollisionId()] = true;
          }
        } else {
          if (selectionObject.collisionId() >= 0) {
            collisionFlag[selectionObject.collisionId()] = true;
          }
        }
      }
    }
  }
  // Todo : Check memory consumption of having so many Process Switches
  PROCESS_SWITCH(JetDerivedDataWriter, processSetupCollisions, "setup the writing for data and MCD based on collisions", true);
  PROCESS_SWITCH(JetDerivedDataWriter, processSetupMcCollisions, "setup the writing for MCP based on mcCollisions", false);
  PROCESS_SWITCH(JetDerivedDataWriter, processSetupAllCollisionsWithDownscaling, "setup the writing of untriggered collisions with downscaling", false);
  PROCESS_SWITCH(JetDerivedDataWriter, processSetupAllMcCollisionsWithDownscaling, "setup the writing of untriggered mccollisions with downscaling", false);
  PROCESS_SWITCH(JetDerivedDataWriter, processSetupEventTriggering, "process software triggers", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::ChargedJets>, processSelectingChargedJets, "process charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::ChargedEventWiseSubtractedJets>, processSelectingChargedEventWiseSubtractedJets, "process charged event-wise subtracted jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::ChargedMCDetectorLevelJets>, processSelectingChargedMCDJets, "process charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::ChargedMCDetectorLevelEventWiseSubtractedJets>, processSelectingChargedMCDetectorLevelEventWiseSubtractedJets, "process charged event-wise subtracted mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::ChargedMCParticleLevelJets>, processSelectingChargedMCPJets, "process charged mcp jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::NeutralJets>, processSelectingNeutralJets, "process neutral jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::NeutralMCDetectorLevelJets>, processSelectingNeutralMCDJets, "process neutral mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::NeutralMCParticleLevelJets>, processSelectingNeutralMCPJets, "process neutral mcp jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::FullJets>, processSelectingFullJets, "process full jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::FullMCDetectorLevelJets>, processSelectingFullMCDJets, "process full mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::FullMCParticleLevelJets>, processSelectingFullMCPJets, "process full mcp jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::D0ChargedJets>, processSelectingD0ChargedJets, "process D0 charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::D0ChargedEventWiseSubtractedJets>, processSelectingD0ChargedEventWiseSubtractedJets, "process D0 event-wise subtracted charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::D0ChargedMCDetectorLevelJets>, processSelectingD0ChargedMCDJets, "process D0 charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::D0ChargedMCDetectorLevelEventWiseSubtractedJets>, processSelectingD0ChargedMCDetectorLevelEventWiseSubtractedJets, "process D0 event-wise subtracted charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::D0ChargedMCParticleLevelJets>, processSelectingD0ChargedMCPJets, "process D0 charged mcp jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::LcChargedJets>, processSelectingLcChargedJets, "process Lc charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::LcChargedEventWiseSubtractedJets>, processSelectingLcChargedEventWiseSubtractedJets, "process Lc event-wise subtracted charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::LcChargedMCDetectorLevelJets>, processSelectingLcChargedMCDJets, "process Lc charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::LcChargedMCDetectorLevelEventWiseSubtractedJets>, processSelectingLcChargedMCDetectorLevelEventWiseSubtractedJets, "process Lc event-wise subtracted charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::LcChargedMCParticleLevelJets>, processSelectingLcChargedMCPJets, "process Lc charged mcp jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::DielectronChargedJets>, processSelectingDielectronChargedJets, "process Dielectron charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::DielectronChargedEventWiseSubtractedJets>, processSelectingDielectronChargedEventWiseSubtractedJets, "process Dielectron event-wise subtracted charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::DielectronChargedMCDetectorLevelJets>, processSelectingDielectronChargedMCDJets, "process Dielectron charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::DielectronChargedMCDetectorLevelEventWiseSubtractedJets>, processSelectingDielectronChargedMCDetectorLevelEventWiseSubtractedJets, "process Dielectron event-wise subtracted charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::DielectronChargedMCParticleLevelJets>, processSelectingDielectronChargedMCPJets, "process Dielectron charged mcp jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::JClusters>, processSelectingClusters, "process EMCal clusters", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processSelectionObjects<aod::JTracks>, processSelectingTracks, "process high pt tracks", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processDoDownscaling<aod::JCollisions>, processCollisionDownscaling, "process downsaling of triggered collisions", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processDoDownscaling<aod::JMcCollisions>, processMcCollisionDownscaling, "process downsaling of triggered mccollisions", false);
  PROCESS_SWITCH(JetDerivedDataWriter, processDoCollisionSelections, "process event selections for saved events", false);

  void processStoreDummyTable(aod::JDummys const&)
  {
    products.storedJDummysTable(1);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processStoreDummyTable, "write out dummy output table", true);

  void processStoreData(soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::JCollisionBCs, aod::JEMCCollisionLbs>::iterator const& collision, soa::Join<aod::JBCs, aod::JBCPIs> const&, soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs> const& tracks, aod::JEMCTracks const& emcTracks, soa::Join<aod::JClusters, aod::JClusterPIs, aod::JClusterTracks> const& clusters, aod::CollisionsD0 const& D0Collisions, aod::CandidatesD0Data const& D0s, aod::CollisionsLc const& LcCollisions, aod::CandidatesLcData const& Lcs, aod::CollisionsDielectron const& DielectronCollisions, aod::CandidatesDielectronData const& Dielectrons)
  {
    std::map<int32_t, int32_t> bcMapping;
    std::map<int32_t, int32_t> trackMapping;

    if (collisionFlag[collision.globalIndex()]) {
      if (config.saveBCsTable) {
        auto bc = collision.bc_as<soa::Join<aod::JBCs, aod::JBCPIs>>();
        if (std::find(bcIndicies.begin(), bcIndicies.end(), bc.globalIndex()) == bcIndicies.end()) {
          products.storedJBCsTable(bc.runNumber(), bc.globalBC(), bc.timestamp(), bc.alias_raw(), bc.selection_raw());
          products.storedJBCParentIndexTable(bc.bcId());
          bcIndicies.push_back(bc.globalIndex());
          bcMapping.insert(std::make_pair(bc.globalIndex(), products.storedJBCsTable.lastIndex()));
        }
      }

      products.storedJCollisionsTable(collision.posX(), collision.posY(), collision.posZ(), collision.multiplicity(), collision.centrality(), collision.trackOccupancyInTimeRange(), collision.eventSel(), collision.alias_raw(), collision.triggerSel());
      products.storedJCollisionsParentIndexTable(collision.collisionId());
      if (config.saveBCsTable) {
        int32_t storedBCID = -1;
        auto JBCIndex = bcMapping.find(collision.bcId());
        if (JBCIndex != bcMapping.end()) {
          storedBCID = JBCIndex->second;
        }
        products.storedJCollisionsBunchCrossingIndexTable(storedBCID);
      }
      if (config.saveClustersTable) {
        products.storedJCollisionsEMCalLabelTable(collision.isAmbiguous(), collision.isEmcalReadout());
      }

      for (const auto& track : tracks) {
        if (config.performTrackSelection && !(track.trackSel() & ~(1 << jetderiveddatautilities::JTrackSel::trackSign))) { // skips tracks that pass no selections. This might cause a problem with tracks matched with clusters. We should generate a track selection purely for cluster matched tracks so that they are kept. This includes also the track pT selction.
          continue;
        }
        if (track.pt() < config.trackPtSelectionMin || std::abs(track.eta()) > config.trackEtaSelectionMax) {
          continue;
        }
        products.storedJTracksTable(products.storedJCollisionsTable.lastIndex(), o2::math_utils::detail::truncateFloatFraction(track.pt(), precisionMomentumMask), o2::math_utils::detail::truncateFloatFraction(track.eta(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.phi(), precisionPositionMask), track.trackSel());
        products.storedJTracksExtraTable(o2::math_utils::detail::truncateFloatFraction(track.dcaX(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.dcaY(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.dcaZ(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.dcaXY(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.dcaXYZ(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.sigmadcaZ(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.sigmadcaXY(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.sigmadcaXYZ(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.sigma1Pt(), precisionMomentumMask));
        products.storedJTracksParentIndexTable(track.trackId());
        trackMapping.insert(std::make_pair(track.globalIndex(), products.storedJTracksTable.lastIndex()));
      }
      if (config.saveClustersTable) {
        for (const auto& cluster : clusters) {
          products.storedJClustersTable(products.storedJCollisionsTable.lastIndex(), cluster.id(), cluster.energy(), cluster.coreEnergy(), cluster.rawEnergy(),
                                        cluster.eta(), cluster.phi(), cluster.m02(), cluster.m20(), cluster.nCells(), cluster.time(), cluster.isExotic(), cluster.distanceToBadChannel(),
                                        cluster.nlm(), cluster.definition(), cluster.leadingCellEnergy(), cluster.subleadingCellEnergy(), cluster.leadingCellNumber(), cluster.subleadingCellNumber());
          products.storedJClustersParentIndexTable(cluster.clusterId());

          std::vector<int32_t> clusterStoredJTrackIDs;
          for (const auto& clusterTrack : cluster.matchedTracks_as<soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs>>()) {
            auto JtrackIndex = trackMapping.find(clusterTrack.globalIndex());
            if (JtrackIndex != trackMapping.end()) {
              clusterStoredJTrackIDs.push_back(JtrackIndex->second);
              auto emcTracksPerTrack = emcTracks.sliceBy(EMCTrackPerTrack, clusterTrack.globalIndex());
              auto emcTrackPerTrack = emcTracksPerTrack.iteratorAt(0);
              products.storedJTracksEMCalTable(JtrackIndex->second, emcTrackPerTrack.etaEmcal(), emcTrackPerTrack.phiEmcal());
            }
          }
          products.storedJClustersMatchedTracksTable(clusterStoredJTrackIDs);
        }
      }

      if (config.saveD0Table) {
        int32_t collisionD0Index = -1;
        for (const auto& D0Collision : D0Collisions) { // should only ever be one
          jethfutilities::fillD0CollisionTable(D0Collision, products.storedD0CollisionsTable, collisionD0Index);
          products.storedD0CollisionIdsTable(products.storedJCollisionsTable.lastIndex());
        }
        for (const auto& D0 : D0s) {
          int32_t D0Index = -1;
          jethfutilities::fillD0CandidateTable<false>(D0, collisionD0Index, products.storedD0sTable, products.storedD0ParsTable, products.storedD0ParExtrasTable, products.storedD0SelsTable, products.storedD0MlsTable, products.storedD0McsTable, D0Index);

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
          products.storedD0IdsTable(products.storedJCollisionsTable.lastIndex(), prong0Id, prong1Id);
        }
      }

      if (config.saveLcTable) {
        int32_t collisionLcIndex = -1;
        for (const auto& LcCollision : LcCollisions) { // should only ever be one
          jethfutilities::fillLcCollisionTable(LcCollision, products.storedLcCollisionsTable, collisionLcIndex);
          products.storedLcCollisionIdsTable(products.storedJCollisionsTable.lastIndex());
        }
        for (const auto& Lc : Lcs) {
          int32_t LcIndex = -1;
          jethfutilities::fillLcCandidateTable<false>(Lc, collisionLcIndex, products.storedLcsTable, products.storedLcParsTable, products.storedLcParExtrasTable, products.storedLcSelsTable, products.storedLcMlsTable, products.storedLcMcsTable, LcIndex);

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
          products.storedLcIdsTable(products.storedJCollisionsTable.lastIndex(), prong0Id, prong1Id, prong2Id);
        }
      }
      if (config.saveDielectronTable) {
        int32_t collisionDielectronIndex = -1;
        for (const auto& DielectronCollision : DielectronCollisions) { // should only ever be one
          jetdqutilities::fillDielectronCollisionTable(DielectronCollision, products.storedDielectronCollisionsTable, collisionDielectronIndex);
          products.storedDielectronCollisionIdsTable(products.storedJCollisionsTable.lastIndex());
        }
        for (const auto& Dielectron : Dielectrons) {
          int32_t DielectronIndex = -1;
          jetdqutilities::fillDielectronCandidateTable(Dielectron, collisionDielectronIndex, products.storedDielectronsTable, DielectronIndex);

          int32_t prong0Id = -1;
          int32_t prong1Id = -1;
          auto JtrackIndex = trackMapping.find(Dielectron.prong0Id());
          if (JtrackIndex != trackMapping.end()) {
            prong0Id = JtrackIndex->second;
          }
          JtrackIndex = trackMapping.find(Dielectron.prong1Id());
          if (JtrackIndex != trackMapping.end()) {
            prong1Id = JtrackIndex->second;
          }
          products.storedDielectronIdsTable(products.storedJCollisionsTable.lastIndex(), prong0Id, prong1Id);
        }
      }
    }
  }
  // process switch for output writing must be last
  // to run after all jet selections
  PROCESS_SWITCH(JetDerivedDataWriter, processStoreData, "write out data output tables", false);

  void processStoreMC(soa::Join<aod::JMcCollisions, aod::JMcCollisionPIs> const& mcCollisions, soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::JCollisionBCs, aod::JMcCollisionLbs, aod::JEMCCollisionLbs> const& collisions, soa::Join<aod::JBCs, aod::JBCPIs> const&, soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs, aod::JMcTrackLbs> const& tracks, aod::JEMCTracks const& emcTracks, soa::Join<aod::JClusters, aod::JClusterPIs, aod::JClusterTracks, aod::JMcClusterLbs> const& clusters, soa::Join<aod::JMcParticles, aod::JMcParticlePIs> const& particles, aod::CollisionsD0 const& D0Collisions, aod::CandidatesD0MCD const& D0s, soa::Join<aod::McCollisionsD0, aod::HfD0McRCollIds> const& D0McCollisions, aod::CandidatesD0MCP const& D0Particles, aod::CollisionsLc const& LcCollisions, aod::CandidatesLcMCD const& Lcs, soa::Join<aod::McCollisionsLc, aod::Hf3PMcRCollIds> const& LcMcCollisions, aod::CandidatesLcMCP const& LcParticles, aod::CollisionsDielectron const& DielectronCollisions, aod::CandidatesDielectronMCD const& Dielectrons, aod::McCollisionsDielectron const& DielectronMcCollisions, aod::CandidatesDielectronMCP const& DielectronParticles)
  {
    std::map<int32_t, int32_t> bcMapping;
    std::map<int32_t, int32_t> paticleMapping;
    std::map<int32_t, int32_t> mcCollisionMapping;
    std::map<int32_t, int32_t> D0CollisionMapping;
    std::map<int32_t, int32_t> LcCollisionMapping;
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

        products.storedJMcCollisionsTable(mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(), mcCollision.weight());
        products.storedJMcCollisionsParentIndexTable(mcCollision.mcCollisionId());
        mcCollisionMapping.insert(std::make_pair(mcCollision.globalIndex(), products.storedJMcCollisionsTable.lastIndex()));

        for (auto particle : particlesPerMcCollision) {
          paticleMapping.insert(std::make_pair(particle.globalIndex(), particleTableIndex));
          particleTableIndex++;
        }
        for (auto particle : particlesPerMcCollision) {

          std::vector<int32_t> mothersId;
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
          products.storedJMcParticlesTable(products.storedJMcCollisionsTable.lastIndex(), o2::math_utils::detail::truncateFloatFraction(particle.pt(), precisionMomentumMask), o2::math_utils::detail::truncateFloatFraction(particle.eta(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(particle.phi(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(particle.y(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(particle.e(), precisionMomentumMask), particle.pdgCode(), particle.getGenStatusCode(), particle.getHepMCStatusCode(), particle.isPhysicalPrimary(), mothersId, daughtersId);
          products.storedJParticlesParentIndexTable(particle.mcParticleId());
        }

        if (config.saveD0Table) {
          const auto d0McCollisionsPerMcCollision = D0McCollisions.sliceBy(D0McCollisionsPerMcCollision, mcCollision.globalIndex());
          int32_t mcCollisionD0Index = -1;
          for (const auto& d0McCollisionPerMcCollision : d0McCollisionsPerMcCollision) { // should only ever be one
            jethfutilities::fillD0McCollisionTable(d0McCollisionPerMcCollision, products.storedD0McCollisionsTable, mcCollisionD0Index);
            products.storedD0McCollisionIdsTable(products.storedJMcCollisionsTable.lastIndex());
          }
          for (const auto& D0Particle : D0Particles) {
            int32_t D0ParticleIndex = -1;
            jethfutilities::fillD0CandidateMcTable(D0Particle, mcCollisionD0Index, products.storedD0ParticlesTable, D0ParticleIndex);
            int32_t D0ParticleId = -1;
            auto JParticleIndex = paticleMapping.find(D0Particle.mcParticleId());
            if (JParticleIndex != paticleMapping.end()) {
              D0ParticleId = JParticleIndex->second;
            }
            products.storedD0ParticleIdsTable(products.storedJMcCollisionsTable.lastIndex(), D0ParticleId);
          }
        }

        if (config.saveLcTable) {
          const auto lcMcCollisionsPerMcCollision = LcMcCollisions.sliceBy(LcMcCollisionsPerMcCollision, mcCollision.globalIndex());
          int32_t mcCollisionLcIndex = -1;
          for (const auto& lcMcCollisionPerMcCollision : lcMcCollisionsPerMcCollision) { // should only ever be one
            jethfutilities::fillLcMcCollisionTable(lcMcCollisionPerMcCollision, products.storedLcMcCollisionsTable, mcCollisionLcIndex);
            products.storedLcMcCollisionIdsTable(products.storedJMcCollisionsTable.lastIndex());
          }
          for (const auto& LcParticle : LcParticles) {
            int32_t LcParticleIndex = -1;
            jethfutilities::fillLcCandidateMcTable(LcParticle, mcCollisionLcIndex, products.storedLcParticlesTable, LcParticleIndex);
            int32_t LcParticleId = -1;
            auto JParticleIndex = paticleMapping.find(LcParticle.mcParticleId());
            if (JParticleIndex != paticleMapping.end()) {
              LcParticleId = JParticleIndex->second;
            }
            products.storedLcParticleIdsTable(products.storedJMcCollisionsTable.lastIndex(), LcParticleId);
          }
        }
        if (config.saveDielectronTable) {
          const auto dielectronMcCollisionsPerMcCollision = DielectronMcCollisions.sliceBy(DielectronMcCollisionsPerMcCollision, mcCollision.globalIndex());
          int32_t mcCollisionDielectronIndex = -1;
          for (const auto& dielectronMcCollisionPerMcCollision : dielectronMcCollisionsPerMcCollision) { // should only ever be one
            jetdqutilities::fillDielectronMcCollisionTable(dielectronMcCollisionPerMcCollision, products.storedDielectronMcCollisionsTable, mcCollisionDielectronIndex);
            products.storedDielectronMcCollisionIdsTable(products.storedJMcCollisionsTable.lastIndex());
          }
          for (const auto& DielectronParticle : DielectronParticles) {
            int32_t DielectronParticleIndex = -1;
            jetdqutilities::fillDielectronCandidateMcTable(DielectronParticle, mcCollisionDielectronIndex, products.storedDielectronParticlesTable, DielectronParticleIndex);
            int32_t DielectronParticleId = -1;
            auto JParticleIndex = paticleMapping.find(DielectronParticle.mcParticleId());
            if (JParticleIndex != paticleMapping.end()) {
              DielectronParticleId = JParticleIndex->second;
            }
            std::vector<int32_t> DielectronMothersId;
            int DielectronDaughtersId[2];
            if (DielectronParticle.has_mothers()) {
              for (auto const& DielectronMother : DielectronParticle.template mothers_as<soa::Join<aod::JMcParticles, aod::JMcParticlePIs>>()) {
                auto JDielectronMotherIndex = paticleMapping.find(DielectronMother.globalIndex());
                if (JDielectronMotherIndex != paticleMapping.end()) {
                  DielectronMothersId.push_back(JDielectronMotherIndex->second);
                }
              }
            }
            auto i = 0;
            if (DielectronParticle.has_daughters()) {
              for (auto const& DielectronDaughter : DielectronParticle.template daughters_as<soa::Join<aod::JMcParticles, aod::JMcParticlePIs>>()) {
                if (i > 1) {
                  break;
                }
                auto JDielectronDaughterIndex = paticleMapping.find(DielectronDaughter.globalIndex());
                if (JDielectronDaughterIndex != paticleMapping.end()) {
                  DielectronDaughtersId[i] = JDielectronDaughterIndex->second;
                }
                i++;
              }
            }
            products.storedDielectronParticleIdsTable(products.storedJMcCollisionsTable.lastIndex(), DielectronParticleId, DielectronMothersId, DielectronDaughtersId);
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
          if (config.saveBCsTable) {
            auto bc = collision.bc_as<soa::Join<aod::JBCs, aod::JBCPIs>>();
            if (std::find(bcIndicies.begin(), bcIndicies.end(), bc.globalIndex()) == bcIndicies.end()) {
              products.storedJBCsTable(bc.runNumber(), bc.globalBC(), bc.timestamp(), bc.alias_raw(), bc.selection_raw());
              products.storedJBCParentIndexTable(bc.bcId());
              bcIndicies.push_back(bc.globalIndex());
              bcMapping.insert(std::make_pair(bc.globalIndex(), products.storedJBCsTable.lastIndex()));
            }
          }

          products.storedJCollisionsTable(collision.posX(), collision.posY(), collision.posZ(), collision.multiplicity(), collision.centrality(), collision.trackOccupancyInTimeRange(), collision.eventSel(), collision.alias_raw(), collision.triggerSel());
          products.storedJCollisionsParentIndexTable(collision.collisionId());

          auto JMcCollisionIndex = mcCollisionMapping.find(mcCollision.globalIndex());
          if (JMcCollisionIndex != mcCollisionMapping.end()) {
            products.storedJMcCollisionsLabelTable(JMcCollisionIndex->second);
          }
          if (config.saveBCsTable) {
            int32_t storedBCID = -1;
            auto JBCIndex = bcMapping.find(collision.bcId());
            if (JBCIndex != bcMapping.end()) {
              storedBCID = JBCIndex->second;
            }
            products.storedJCollisionsBunchCrossingIndexTable(storedBCID);
          }
          if (config.saveClustersTable) {
            products.storedJCollisionsEMCalLabelTable(collision.isAmbiguous(), collision.isEmcalReadout());
          }

          const auto tracksPerCollision = tracks.sliceBy(TracksPerCollision, collision.globalIndex());
          for (const auto& track : tracksPerCollision) {
            if (config.performTrackSelection && !(track.trackSel() & ~(1 << jetderiveddatautilities::JTrackSel::trackSign))) { // skips tracks that pass no selections. This might cause a problem with tracks matched with clusters. We should generate a track selection purely for cluster matched tracks so that they are kept
              continue;
            }
            if (track.pt() < config.trackPtSelectionMin || std::abs(track.eta()) > config.trackEtaSelectionMax) {
              continue;
            }
            products.storedJTracksTable(products.storedJCollisionsTable.lastIndex(), o2::math_utils::detail::truncateFloatFraction(track.pt(), precisionMomentumMask), o2::math_utils::detail::truncateFloatFraction(track.eta(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.phi(), precisionPositionMask), track.trackSel());
            products.storedJTracksExtraTable(o2::math_utils::detail::truncateFloatFraction(track.dcaX(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.dcaY(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.dcaZ(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.dcaXY(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.dcaXYZ(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.sigmadcaZ(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.sigmadcaXY(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.sigmadcaXYZ(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.sigma1Pt(), precisionMomentumMask));
            products.storedJTracksParentIndexTable(track.trackId());

            if (track.has_mcParticle()) {
              auto JParticleIndex = paticleMapping.find(track.mcParticleId());
              if (JParticleIndex != paticleMapping.end()) {
                products.storedJMcTracksLabelTable(JParticleIndex->second);
              } else {
                products.storedJMcTracksLabelTable(-1); // this can happen because there are some tracks that are reconstucted in a wrong collision, but their original McCollision did not pass the required cuts so that McParticle is not saved. These are very few but we should look into them further and see what to do about them
              }
            } else {
              products.storedJMcTracksLabelTable(-1);
            }
            trackMapping.insert(std::make_pair(track.globalIndex(), products.storedJTracksTable.lastIndex()));
          }
          if (config.saveClustersTable) {
            const auto clustersPerCollision = clusters.sliceBy(ClustersPerCollision, collision.globalIndex());
            for (const auto& cluster : clustersPerCollision) {
              products.storedJClustersTable(products.storedJCollisionsTable.lastIndex(), cluster.id(), cluster.energy(), cluster.coreEnergy(), cluster.rawEnergy(),
                                            cluster.eta(), cluster.phi(), cluster.m02(), cluster.m20(), cluster.nCells(), cluster.time(), cluster.isExotic(), cluster.distanceToBadChannel(),
                                            cluster.nlm(), cluster.definition(), cluster.leadingCellEnergy(), cluster.subleadingCellEnergy(), cluster.leadingCellNumber(), cluster.subleadingCellNumber());
              products.storedJClustersParentIndexTable(cluster.clusterId());

              std::vector<int32_t> clusterStoredJTrackIDs;
              for (const auto& clusterTrack : cluster.matchedTracks_as<soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs>>()) {
                auto JtrackIndex = trackMapping.find(clusterTrack.globalIndex());
                if (JtrackIndex != trackMapping.end()) {
                  clusterStoredJTrackIDs.push_back(JtrackIndex->second);
                  const auto emcTracksPerTrack = emcTracks.sliceBy(EMCTrackPerTrack, clusterTrack.globalIndex());
                  auto emcTrackPerTrack = emcTracksPerTrack.iteratorAt(0);
                  products.storedJTracksEMCalTable(JtrackIndex->second, emcTrackPerTrack.etaEmcal(), emcTrackPerTrack.phiEmcal());
                }
              }
              products.storedJClustersMatchedTracksTable(clusterStoredJTrackIDs);

              std::vector<int32_t> clusterStoredJParticleIDs;
              for (const auto& clusterParticleId : cluster.mcParticleIds()) {
                auto JParticleIndex = paticleMapping.find(clusterParticleId);
                if (JParticleIndex != paticleMapping.end()) {
                  clusterStoredJParticleIDs.push_back(JParticleIndex->second);
                }
              }
              std::vector<float> amplitudeA;
              auto amplitudeASpan = cluster.amplitudeA();
              std::copy(amplitudeASpan.begin(), amplitudeASpan.end(), std::back_inserter(amplitudeA));
              products.storedJMcClustersLabelTable(clusterStoredJParticleIDs, amplitudeA);
            }
          }

          if (config.saveD0Table) {
            const auto d0CollisionsPerCollision = D0Collisions.sliceBy(D0CollisionsPerCollision, collision.globalIndex());
            int32_t collisionD0Index = -1;
            for (const auto& d0CollisionPerCollision : d0CollisionsPerCollision) { // should only ever be one
              jethfutilities::fillD0CollisionTable(d0CollisionPerCollision, products.storedD0CollisionsTable, collisionD0Index);
              products.storedD0CollisionIdsTable(products.storedJCollisionsTable.lastIndex());
              D0CollisionMapping.insert(std::make_pair(d0CollisionPerCollision.globalIndex(), products.storedD0CollisionsTable.lastIndex()));
            }
            const auto d0sPerCollision = D0s.sliceBy(D0sPerCollision, collision.globalIndex());
            for (const auto& D0 : d0sPerCollision) {
              int32_t D0Index = -1;
              jethfutilities::fillD0CandidateTable<true>(D0, collisionD0Index, products.storedD0sTable, products.storedD0ParsTable, products.storedD0ParExtrasTable, products.storedD0SelsTable, products.storedD0MlsTable, products.storedD0McsTable, D0Index);

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
              products.storedD0IdsTable(products.storedJCollisionsTable.lastIndex(), prong0Id, prong1Id);
            }
          }

          if (config.saveLcTable) {
            const auto lcCollisionsPerCollision = LcCollisions.sliceBy(LcCollisionsPerCollision, collision.globalIndex());
            int32_t collisionLcIndex = -1;
            for (const auto& lcCollisionPerCollision : lcCollisionsPerCollision) { // should only ever be one
              jethfutilities::fillLcCollisionTable(lcCollisionPerCollision, products.storedLcCollisionsTable, collisionLcIndex);
              products.storedLcCollisionIdsTable(products.storedJCollisionsTable.lastIndex());
              LcCollisionMapping.insert(std::make_pair(lcCollisionPerCollision.globalIndex(), products.storedLcCollisionsTable.lastIndex()));
            }
            const auto lcsPerCollision = Lcs.sliceBy(LcsPerCollision, collision.globalIndex());
            for (const auto& Lc : lcsPerCollision) {
              int32_t LcIndex = -1;
              jethfutilities::fillLcCandidateTable<true>(Lc, collisionLcIndex, products.storedLcsTable, products.storedLcParsTable, products.storedLcParExtrasTable, products.storedLcSelsTable, products.storedLcMlsTable, products.storedLcMcsTable, LcIndex);

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
              products.storedLcIdsTable(products.storedJCollisionsTable.lastIndex(), prong0Id, prong1Id, prong2Id);
            }
          }

          if (config.saveDielectronTable) {
            const auto dielectronCollisionsPerCollision = DielectronCollisions.sliceBy(DielectronCollisionsPerCollision, collision.globalIndex());
            int32_t collisionDielectronIndex = -1;
            for (const auto& dielectronCollisionPerCollision : dielectronCollisionsPerCollision) { // should only ever be one
              jetdqutilities::fillDielectronCollisionTable(dielectronCollisionPerCollision, products.storedDielectronCollisionsTable, collisionDielectronIndex);
              products.storedDielectronCollisionIdsTable(products.storedJCollisionsTable.lastIndex());
              // D0CollisionMapping.insert(std::make_pair(d0CollisionPerCollision.globalIndex(), products.storedD0CollisionsTable.lastIndex())); // if DielectronMCCollisions are indexed to Dielectron Collisions then this can be added
            }
            const auto dielectronsPerCollision = Dielectrons.sliceBy(DielectronsPerCollision, collision.globalIndex());
            for (const auto& Dielectron : dielectronsPerCollision) {
              int32_t DielectronIndex = -1;
              jetdqutilities::fillDielectronCandidateTable(Dielectron, collisionDielectronIndex, products.storedDielectronsTable, DielectronIndex);

              int32_t prong0Id = -1;
              int32_t prong1Id = -1;
              auto JtrackIndex = trackMapping.find(Dielectron.prong0Id());
              if (JtrackIndex != trackMapping.end()) {
                prong0Id = JtrackIndex->second;
              }
              JtrackIndex = trackMapping.find(Dielectron.prong1Id());
              if (JtrackIndex != trackMapping.end()) {
                prong1Id = JtrackIndex->second;
              }
              products.storedDielectronIdsTable(products.storedJCollisionsTable.lastIndex(), prong0Id, prong1Id);
            }
          }
        }

        if (config.saveD0Table) {
          const auto d0McCollisionsPerMcCollision = D0McCollisions.sliceBy(D0McCollisionsPerMcCollision, mcCollision.globalIndex());
          for (const auto& d0McCollisionPerMcCollision : d0McCollisionsPerMcCollision) { // should just be one
            std::vector<int32_t> d0CollisionIDs;
            for (auto const& d0CollisionPerMcCollision : d0McCollisionPerMcCollision.template hfCollBases_as<aod::CollisionsD0>()) {
              auto d0CollisionIndex = D0CollisionMapping.find(d0CollisionPerMcCollision.globalIndex());
              if (d0CollisionIndex != D0CollisionMapping.end()) {
                d0CollisionIDs.push_back(d0CollisionIndex->second);
              }
            }
            products.storedD0McCollisionsMatchingTable(d0CollisionIDs);
          }
        }

        if (config.saveLcTable) {
          const auto lcMcCollisionsPerMcCollision = LcMcCollisions.sliceBy(LcMcCollisionsPerMcCollision, mcCollision.globalIndex());
          for (const auto& lcMcCollisionPerMcCollision : lcMcCollisionsPerMcCollision) { // should just be one
            std::vector<int32_t> lcCollisionIDs;
            for (auto const& lcCollisionPerMcCollision : lcMcCollisionPerMcCollision.template hfCollBases_as<aod::CollisionsLc>()) {
              auto lcCollisionIndex = LcCollisionMapping.find(lcCollisionPerMcCollision.globalIndex());
              if (lcCollisionIndex != LcCollisionMapping.end()) {
                lcCollisionIDs.push_back(lcCollisionIndex->second);
              }
            }
            products.storedLcMcCollisionsMatchingTable(lcCollisionIDs);
          }
        }
      }
    }
  }
  // process switch for output writing must be last
  // to run after all jet selections
  PROCESS_SWITCH(JetDerivedDataWriter, processStoreMC, "write out data output tables for mc", false);

  void processStoreMCP(soa::Join<aod::JMcCollisions, aod::JMcCollisionPIs> const& mcCollisions, soa::Join<aod::JMcParticles, aod::JMcParticlePIs> const& particles, aod::McCollisionsD0 const& D0McCollisions, aod::CandidatesD0MCP const& D0Particles, aod::McCollisionsLc const& LcMcCollisions, aod::CandidatesLcMCP const& LcParticles, aod::McCollisionsDielectron const& DielectronMcCollisions, aod::CandidatesDielectronMCP const& DielectronParticles)
  {

    int particleTableIndex = 0;
    for (auto mcCollision : mcCollisions) {
      if (McCollisionFlag[mcCollision.globalIndex()]) { // you can also check if any of its detector level counterparts are correct
        std::map<int32_t, int32_t> paticleMapping;

        products.storedJMcCollisionsTable(mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(), mcCollision.weight());
        products.storedJMcCollisionsParentIndexTable(mcCollision.mcCollisionId());

        const auto particlesPerMcCollision = particles.sliceBy(ParticlesPerMcCollision, mcCollision.globalIndex());

        for (auto particle : particlesPerMcCollision) {
          paticleMapping.insert(std::make_pair(particle.globalIndex(), particleTableIndex));
          particleTableIndex++;
        }
        for (auto particle : particlesPerMcCollision) {
          std::vector<int32_t> mothersId;
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
          products.storedJMcParticlesTable(products.storedJMcCollisionsTable.lastIndex(), o2::math_utils::detail::truncateFloatFraction(particle.pt(), precisionMomentumMask), o2::math_utils::detail::truncateFloatFraction(particle.eta(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(particle.phi(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(particle.y(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(particle.e(), precisionMomentumMask), particle.pdgCode(), particle.getGenStatusCode(), particle.getHepMCStatusCode(), particle.isPhysicalPrimary(), mothersId, daughtersId);
          products.storedJParticlesParentIndexTable(particle.mcParticleId());
        }

        if (config.saveD0Table) {
          const auto d0McCollisionsPerMcCollision = D0McCollisions.sliceBy(D0McCollisionsPerMcCollision, mcCollision.globalIndex());
          int32_t mcCollisionD0Index = -1;
          for (const auto& d0McCollisionPerMcCollision : d0McCollisionsPerMcCollision) { // should only ever be one
            jethfutilities::fillD0McCollisionTable(d0McCollisionPerMcCollision, products.storedD0McCollisionsTable, mcCollisionD0Index);
            products.storedD0McCollisionIdsTable(products.storedJMcCollisionsTable.lastIndex());
          }
          for (const auto& D0Particle : D0Particles) {
            int32_t D0ParticleIndex = -1;
            jethfutilities::fillD0CandidateMcTable(D0Particle, mcCollisionD0Index, products.storedD0ParticlesTable, D0ParticleIndex);
            int32_t D0ParticleId = -1;
            auto JParticleIndex = paticleMapping.find(D0Particle.mcParticleId());
            if (JParticleIndex != paticleMapping.end()) {
              D0ParticleId = JParticleIndex->second;
            }
            products.storedD0ParticleIdsTable(products.storedJMcCollisionsTable.lastIndex(), D0ParticleId);
          }
        }
        if (config.saveLcTable) {
          const auto lcMcCollisionsPerMcCollision = LcMcCollisions.sliceBy(LcMcCollisionsPerMcCollision, mcCollision.globalIndex());
          int32_t mcCollisionLcIndex = -1;
          for (const auto& lcMcCollisionPerMcCollision : lcMcCollisionsPerMcCollision) { // should only ever be one
            jethfutilities::fillLcMcCollisionTable(lcMcCollisionPerMcCollision, products.storedLcMcCollisionsTable, mcCollisionLcIndex);
            products.storedLcMcCollisionIdsTable(products.storedJMcCollisionsTable.lastIndex());
          }
          for (const auto& LcParticle : LcParticles) {
            int32_t LcParticleIndex = -1;
            jethfutilities::fillLcCandidateMcTable(LcParticle, mcCollisionLcIndex, products.storedLcParticlesTable, LcParticleIndex);
            int32_t LcParticleId = -1;
            auto JParticleIndex = paticleMapping.find(LcParticle.mcParticleId());
            if (JParticleIndex != paticleMapping.end()) {
              LcParticleId = JParticleIndex->second;
            }
            products.storedLcParticleIdsTable(products.storedJMcCollisionsTable.lastIndex(), LcParticleId);
          }
        }
        if (config.saveDielectronTable) {
          const auto dielectronMcCollisionsPerMcCollision = DielectronMcCollisions.sliceBy(DielectronMcCollisionsPerMcCollision, mcCollision.globalIndex());
          int32_t mcCollisionDielectronIndex = -1;
          for (const auto& dielectronMcCollisionPerMcCollision : dielectronMcCollisionsPerMcCollision) { // should only ever be one
            jetdqutilities::fillDielectronMcCollisionTable(dielectronMcCollisionPerMcCollision, products.storedDielectronMcCollisionsTable, mcCollisionDielectronIndex);
            products.storedDielectronMcCollisionIdsTable(products.storedJMcCollisionsTable.lastIndex());
          }
          for (const auto& DielectronParticle : DielectronParticles) {
            int32_t DielectronParticleIndex = -1;
            jetdqutilities::fillDielectronCandidateMcTable(DielectronParticle, mcCollisionDielectronIndex, products.storedDielectronParticlesTable, DielectronParticleIndex);
            int32_t DielectronParticleId = -1;
            auto JParticleIndex = paticleMapping.find(DielectronParticle.mcParticleId());
            if (JParticleIndex != paticleMapping.end()) {
              DielectronParticleId = JParticleIndex->second;
            }
            std::vector<int32_t> DielectronMothersId;
            int DielectronDaughtersId[2];
            if (DielectronParticle.has_mothers()) {
              for (auto const& DielectronMother : DielectronParticle.template mothers_as<soa::Join<aod::JMcParticles, aod::JMcParticlePIs>>()) {
                auto JDielectronMotherIndex = paticleMapping.find(DielectronMother.globalIndex());
                if (JDielectronMotherIndex != paticleMapping.end()) {
                  DielectronMothersId.push_back(JDielectronMotherIndex->second);
                }
              }
            }
            auto i = 0;
            if (DielectronParticle.has_daughters()) {
              for (auto const& DielectronDaughter : DielectronParticle.template daughters_as<soa::Join<aod::JMcParticles, aod::JMcParticlePIs>>()) {
                if (i > 1) {
                  break;
                }
                auto JDielectronDaughterIndex = paticleMapping.find(DielectronDaughter.globalIndex());
                if (JDielectronDaughterIndex != paticleMapping.end()) {
                  DielectronDaughtersId[i] = JDielectronDaughterIndex->second;
                }
                i++;
              }
            }
            products.storedDielectronParticleIdsTable(products.storedJMcCollisionsTable.lastIndex(), DielectronParticleId, DielectronMothersId, DielectronDaughtersId);
          }
        }
      }
    }
  }
  // process switch for output writing must be last
  // to run after all jet selections
  PROCESS_SWITCH(JetDerivedDataWriter, processStoreMCP, "write out data output tables for mcp", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetDerivedDataWriter>(cfgc, TaskName{"jet-deriveddata-writer"}));

  return WorkflowSpec{tasks};
}
