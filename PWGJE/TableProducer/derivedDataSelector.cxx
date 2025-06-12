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

/// \file jetderiveddataselector.cxx
/// \brief Task to store decision of which events to skim for producing jet framework tables (aod::JetCollisions, aod::JetTracks, aod::JetClusters, ...)
/// while adjusting indices accordingly
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
/// \author Jochen Klein <jochen.klein@cern.ch>

#include <MathUtils/Utils.h>
#include <algorithm>
#include <string>
#include <vector>

#include <TRandom3.h>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetReducedDataSelector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetDerivedDataSelector {

  Produces<aod::JCollisionSelections> collisionSelectionsTable;
  Produces<aod::JMcCollisionSelections> mcCollisionSelectionsTable;

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
    Configurable<float> thresholdChargedDplusJetPtMin{"thresholdChargedDplusJetPtMin", 0.0, "Minimum charged Dplus jet pt to accept event"};
    Configurable<float> thresholdChargedEventWiseSubtractedDplusJetPtMin{"thresholdChargedEventWiseSubtractedDplusJetPtMin", 0.0, "Minimum charged event-wise subtracted Dplus jet pt to accept event"};
    Configurable<float> thresholdChargedDplusMCPJetPtMin{"thresholdChargedDplusMCPJetPtMin", 0.0, "Minimum charged Dplus mcp jet pt to accept event"};
    Configurable<float> thresholdChargedLcJetPtMin{"thresholdChargedLcJetPtMin", 0.0, "Minimum charged Lc jet pt to accept event"};
    Configurable<float> thresholdChargedEventWiseSubtractedLcJetPtMin{"thresholdChargedEventWiseSubtractedLcJetPtMin", 0.0, "Minimum charged event-wise subtracted Lc jet pt to accept event"};
    Configurable<float> thresholdChargedLcMCPJetPtMin{"thresholdChargedLcMCPJetPtMin", 0.0, "Minimum charged Lc mcp jet pt to accept event"};
    Configurable<float> thresholdChargedBplusJetPtMin{"thresholdChargedBplusJetPtMin", 0.0, "Minimum charged Bplus jet pt to accept event"};
    Configurable<float> thresholdChargedEventWiseSubtractedBplusJetPtMin{"thresholdChargedEventWiseSubtractedBplusJetPtMin", 0.0, "Minimum charged event-wise subtracted Bplus jet pt to accept event"};
    Configurable<float> thresholdChargedBplusMCPJetPtMin{"thresholdChargedBplusMCPJetPtMin", 0.0, "Minimum charged Bplus mcp jet pt to accept event"};
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

    Configurable<std::string> triggerMasks{"triggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};
  } config;

  std::vector<bool> collisionFlag;
  std::vector<bool> McCollisionFlag;

  TRandom3 randomNumber;

  std::vector<int> triggerMaskBits;
  void init(InitContext&)
  {
    randomNumber.SetSeed(0);
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(config.triggerMasks);
  }

  PresliceUnsorted<soa::Join<aod::JCollisions, aod::JMcCollisionLbs>> CollisionsPerMcCollision = aod::jmccollisionlb::mcCollisionId;

  void processSetupCollisions(aod::JCollisions const& collisions)
  {
    collisionFlag.clear();
    collisionFlag.resize(collisions.size(), false);
  }

  void processSetupMcCollisions(aod::JMcCollisions const& mcCollisions)
  {
    McCollisionFlag.clear();
    McCollisionFlag.resize(mcCollisions.size(), false);
  }

  void processSelectMcCollisionsPerCollision(aod::JMcCollisions const& mcCollisions, soa::Join<aod::JCollisions, aod::JMcCollisionLbs> const& collisions)
  {
    for (auto mcCollision : mcCollisions) {
      const auto collisionsPerMcCollision = collisions.sliceBy(CollisionsPerMcCollision, mcCollision.globalIndex());
      for (auto collision : collisionsPerMcCollision) {
        if (collisionFlag[collision.globalIndex()]) {
          McCollisionFlag[mcCollision.globalIndex()] = true;
        }
      }
    }
  }

  void processSelectCollisionsPerMcCollision(soa::Join<aod::JCollisions, aod::JMcCollisionLbs>::iterator const& collision)
  {
    if (McCollisionFlag[collision.mcCollisionId()]) {
      collisionFlag[collision.globalIndex()] = true;
    }
  }

  void processSetupAllCollisionsWithDownscaling(aod::JCollisions const& collisions)
  {
    collisionFlag.clear();
    collisionFlag.resize(collisions.size(), false);
    for (const auto& collision : collisions) {
      if (randomNumber.Integer(config.downscaleFactor) == 0) {
        collisionFlag[collision.globalIndex()] = true;
      }
    }
  }

  void processSetupAllMcCollisionsWithDownscaling(aod::JMcCollisions const& mcCollisions)
  {
    McCollisionFlag.clear();
    McCollisionFlag.resize(mcCollisions.size(), false);
    for (const auto& mcCollision : mcCollisions) {
      if (randomNumber.Integer(config.downscaleFactor) == 0) {
        McCollisionFlag[mcCollision.globalIndex()] = true;
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
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::DplusChargedJets> || std::is_same_v<std::decay_t<T>, aod::DplusChargedMCDetectorLevelJets>) {
      selectionObjectPtMin = config.thresholdChargedDplusJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::DplusChargedEventWiseSubtractedJets> || std::is_same_v<std::decay_t<T>, aod::DplusChargedMCDetectorLevelEventWiseSubtractedJets>) {
      selectionObjectPtMin = config.thresholdChargedEventWiseSubtractedDplusJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::DplusChargedMCParticleLevelJets>) {
      selectionObjectPtMin = config.thresholdChargedDplusMCPJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::LcChargedJets> || std::is_same_v<std::decay_t<T>, aod::LcChargedMCDetectorLevelJets>) {
      selectionObjectPtMin = config.thresholdChargedLcJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::LcChargedEventWiseSubtractedJets> || std::is_same_v<std::decay_t<T>, aod::LcChargedMCDetectorLevelEventWiseSubtractedJets>) {
      selectionObjectPtMin = config.thresholdChargedEventWiseSubtractedLcJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::LcChargedMCParticleLevelJets>) {
      selectionObjectPtMin = config.thresholdChargedLcMCPJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::BplusChargedJets> || std::is_same_v<std::decay_t<T>, aod::BplusChargedMCDetectorLevelJets>) {
      selectionObjectPtMin = config.thresholdChargedBplusJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::BplusChargedEventWiseSubtractedJets> || std::is_same_v<std::decay_t<T>, aod::BplusChargedMCDetectorLevelEventWiseSubtractedJets>) {
      selectionObjectPtMin = config.thresholdChargedEventWiseSubtractedBplusJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::BplusChargedMCParticleLevelJets>) {
      selectionObjectPtMin = config.thresholdChargedBplusMCPJetPtMin;
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
          if (selectionObject.pt() < config.trackPtSelectionMin || std::abs(selectionObject.eta()) > config.trackEtaSelectionMax) {
            continue;
          }
        }
        if (selectionObject.pt() >= selectionObjectPtMin) {
          isTriggerObject = true;
        }
      }
      if (isTriggerObject) {
        if constexpr (std::is_same_v<std::decay_t<T>, aod::ChargedMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::NeutralMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::FullMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::D0ChargedMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::DplusChargedMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::LcChargedMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::BplusChargedMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::DielectronChargedMCParticleLevelJets>) {
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
  PROCESS_SWITCH(JetDerivedDataSelector, processSetupCollisions, "setup the writing for data and MCD based on collisions", true);
  PROCESS_SWITCH(JetDerivedDataSelector, processSetupMcCollisions, "setup the writing for MCP based on mcCollisions", false);
  PROCESS_SWITCH(JetDerivedDataSelector, processSetupAllCollisionsWithDownscaling, "setup the writing of untriggered collisions with downscaling", false);
  PROCESS_SWITCH(JetDerivedDataSelector, processSetupAllMcCollisionsWithDownscaling, "setup the writing of untriggered mccollisions with downscaling", false);
  PROCESS_SWITCH(JetDerivedDataSelector, processSetupEventTriggering, "process software triggers", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::ChargedJets>, processSelectingChargedJets, "process charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::ChargedEventWiseSubtractedJets>, processSelectingChargedEventWiseSubtractedJets, "process charged event-wise subtracted jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::ChargedMCDetectorLevelJets>, processSelectingChargedMCDJets, "process charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::ChargedMCDetectorLevelEventWiseSubtractedJets>, processSelectingChargedMCDetectorLevelEventWiseSubtractedJets, "process charged event-wise subtracted mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::ChargedMCParticleLevelJets>, processSelectingChargedMCPJets, "process charged mcp jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::NeutralJets>, processSelectingNeutralJets, "process neutral jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::NeutralMCDetectorLevelJets>, processSelectingNeutralMCDJets, "process neutral mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::NeutralMCParticleLevelJets>, processSelectingNeutralMCPJets, "process neutral mcp jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::FullJets>, processSelectingFullJets, "process full jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::FullMCDetectorLevelJets>, processSelectingFullMCDJets, "process full mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::FullMCParticleLevelJets>, processSelectingFullMCPJets, "process full mcp jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::D0ChargedJets>, processSelectingD0ChargedJets, "process D0 charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::D0ChargedEventWiseSubtractedJets>, processSelectingD0ChargedEventWiseSubtractedJets, "process D0 event-wise subtracted charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::D0ChargedMCDetectorLevelJets>, processSelectingD0ChargedMCDJets, "process D0 charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::D0ChargedMCDetectorLevelEventWiseSubtractedJets>, processSelectingD0ChargedMCDetectorLevelEventWiseSubtractedJets, "process D0 event-wise subtracted charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::D0ChargedMCParticleLevelJets>, processSelectingD0ChargedMCPJets, "process D0 charged mcp jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::DplusChargedJets>, processSelectingDplusChargedJets, "process Dplus charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::DplusChargedEventWiseSubtractedJets>, processSelectingDplusChargedEventWiseSubtractedJets, "process Dplus event-wise subtracted charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::DplusChargedMCDetectorLevelJets>, processSelectingDplusChargedMCDJets, "process Dplus charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::DplusChargedMCDetectorLevelEventWiseSubtractedJets>, processSelectingDplusChargedMCDetectorLevelEventWiseSubtractedJets, "process Dplus event-wise subtracted charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::DplusChargedMCParticleLevelJets>, processSelectingDplusChargedMCPJets, "process Dplus charged mcp jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::LcChargedJets>, processSelectingLcChargedJets, "process Lc charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::LcChargedEventWiseSubtractedJets>, processSelectingLcChargedEventWiseSubtractedJets, "process Lc event-wise subtracted charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::LcChargedMCDetectorLevelJets>, processSelectingLcChargedMCDJets, "process Lc charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::LcChargedMCDetectorLevelEventWiseSubtractedJets>, processSelectingLcChargedMCDetectorLevelEventWiseSubtractedJets, "process Lc event-wise subtracted charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::LcChargedMCParticleLevelJets>, processSelectingLcChargedMCPJets, "process Lc charged mcp jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::BplusChargedJets>, processSelectingBplusChargedJets, "process Bplus charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::BplusChargedEventWiseSubtractedJets>, processSelectingBplusChargedEventWiseSubtractedJets, "process Bplus event-wise subtracted charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::BplusChargedMCDetectorLevelJets>, processSelectingBplusChargedMCDJets, "process Bplus charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::BplusChargedMCDetectorLevelEventWiseSubtractedJets>, processSelectingBplusChargedMCDetectorLevelEventWiseSubtractedJets, "process Bplus event-wise subtracted charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::BplusChargedMCParticleLevelJets>, processSelectingBplusChargedMCPJets, "process Bplus charged mcp jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::DielectronChargedJets>, processSelectingDielectronChargedJets, "process Dielectron charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::DielectronChargedEventWiseSubtractedJets>, processSelectingDielectronChargedEventWiseSubtractedJets, "process Dielectron event-wise subtracted charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::DielectronChargedMCDetectorLevelJets>, processSelectingDielectronChargedMCDJets, "process Dielectron charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::DielectronChargedMCDetectorLevelEventWiseSubtractedJets>, processSelectingDielectronChargedMCDetectorLevelEventWiseSubtractedJets, "process Dielectron event-wise subtracted charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::DielectronChargedMCParticleLevelJets>, processSelectingDielectronChargedMCPJets, "process Dielectron charged mcp jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::JClusters>, processSelectingClusters, "process EMCal clusters", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processSelectionObjects<aod::JTracks>, processSelectingTracks, "process high pt tracks", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processDoDownscaling<aod::JCollisions>, processCollisionDownscaling, "process downsaling of triggered collisions", false);
  PROCESS_SWITCH_FULL(JetDerivedDataSelector, processDoDownscaling<aod::JMcCollisions>, processMcCollisionDownscaling, "process downsaling of triggered mccollisions", false);
  PROCESS_SWITCH(JetDerivedDataSelector, processDoCollisionSelections, "process event selections for saved events", false);
  PROCESS_SWITCH(JetDerivedDataSelector, processSelectMcCollisionsPerCollision, "select McCollisions due to a triggered reconstructed collision", false);
  PROCESS_SWITCH(JetDerivedDataSelector, processSelectCollisionsPerMcCollision, "select collisions due to a triggered McCollision", false);

  void processStoreCollisionDecision(aod::JCollision const& collision)
  {
    if (collisionFlag[collision.globalIndex()]) {
      collisionSelectionsTable(true);
    } else {
      collisionSelectionsTable(false);
    }
  }
  PROCESS_SWITCH(JetDerivedDataSelector, processStoreCollisionDecision, "write out decision of storing collision", true);

  void processStoreMcCollisionDecision(aod::JMcCollision const& mcCollision)
  {
    if (McCollisionFlag[mcCollision.globalIndex()]) {
      mcCollisionSelectionsTable(true);
    } else {
      mcCollisionSelectionsTable(false);
    }
  }
  PROCESS_SWITCH(JetDerivedDataSelector, processStoreMcCollisionDecision, "write out decision of storing mcCollision", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetDerivedDataSelector>(cfgc, TaskName{"jet-deriveddata-selector"}));

  return WorkflowSpec{tasks};
}
