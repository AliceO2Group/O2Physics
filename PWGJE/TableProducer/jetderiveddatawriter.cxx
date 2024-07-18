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
#include "PWGJE/Core/JetDQUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetDerivedDataWriter {

  struct : ConfigurableGroup {
    Configurable<float> chargedJetPtMin{"chargedJetPtMin", 0.0, "Minimum charged jet pt to accept event"};
    Configurable<float> chargedEventWiseSubtractedJetPtMin{"chargedEventWiseSubtractedJetPtMin", 0.0, "Minimum charged event-wise subtracted jet pt to accept event"};
    Configurable<float> chargedMCPJetPtMin{"chargedMCPJetPtMin", 0.0, "Minimum charged mcp jet pt to accept event"};
    Configurable<float> neutralJetPtMin{"neutralJetPtMin", 0.0, "Minimum neutral jet pt to accept event"};
    Configurable<float> neutralMCPJetPtMin{"neutralMCPJetPtMin", 0.0, "Minimum neutal mcp jet pt to accept event"};
    Configurable<float> fullJetPtMin{"fullJetPtMin", 0.0, "Minimum full jet pt to accept event"};
    Configurable<float> fullMCPJetPtMin{"fullMCPJetPtMin", 0.0, "Minimum full mcp jet pt to accept event"};
    Configurable<float> chargedD0JetPtMin{"chargedD0JetPtMin", 0.0, "Minimum charged D0 jet pt to accept event"};
    Configurable<float> chargedEventWiseSubtractedD0JetPtMin{"chargedEventWiseSubtractedD0JetPtMin", 0.0, "Minimum charged event-wise subtracted D0 jet pt to accept event"};
    Configurable<float> chargedD0MCPJetPtMin{"chargedD0MCPJetPtMin", 0.0, "Minimum charged D0 mcp jet pt to accept event"};
    Configurable<float> chargedLcJetPtMin{"chargedLcJetPtMin", 0.0, "Minimum charged Lc jet pt to accept event"};
    Configurable<float> chargedEventWiseSubtractedLcJetPtMin{"chargedEventWiseSubtractedLcJetPtMin", 0.0, "Minimum charged event-wise subtracted Lc jet pt to accept event"};
    Configurable<float> chargedLcMCPJetPtMin{"chargedLcMCPJetPtMin", 0.0, "Minimum charged Lc mcp jet pt to accept event"};
    Configurable<float> chargedDielectronJetPtMin{"chargedDielectronJetPtMin", 0.0, "Minimum charged Dielectron jet pt to accept event"};
    Configurable<float> chargedEventWiseSubtractedDielectronJetPtMin{"chargedEventWiseSubtractedDielectronJetPtMin", 0.0, "Minimum charged event-wise subtracted Dielectron jet pt to accept event"};
    Configurable<float> chargedDielectronMCPJetPtMin{"chargedDielectronMCPJetPtMin", 0.0, "Minimum charged Dielectron mcp jet pt to accept event"};

    Configurable<bool> performTrackSelection{"performTrackSelection", true, "only save tracks that pass one of the track selections"};
    Configurable<bool> saveBCsTable{"saveBCsTable", true, "save the bunch crossing table to the output"};
    Configurable<bool> saveClustersTable{"saveClustersTable", true, "save the clusters table to the output"};
    Configurable<bool> saveD0Table{"saveD0Table", false, "save the D0 table to the output"};
    Configurable<bool> saveLcTable{"saveLcTable", false, "save the Lc table to the output"};
    Configurable<bool> saveDielectronTable{"saveDielectronTable", false, "save the Dielectron table to the output"};

    Configurable<std::string> eventSelectionForCounting{"eventSelectionForCounting", "sel8", "choose event selection for collision counter"};
  } config;

  struct : ProducesGroup {
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
  Preslice<soa::Join<McCollisionsD0, aod::HfD0McRCollIds>> D0McCollisionsPerMcCollision = aod::jd0indices::mcCollisionId;
  Preslice<soa::Join<McCollisionsLc, aod::Hf3PMcRCollIds>> LcMcCollisionsPerMcCollision = aod::jlcindices::mcCollisionId;
  Preslice<McCollisionsDielectron> DielectronMcCollisionsPerMcCollision = aod::jdielectronindices::mcCollisionId;
  Preslice<CollisionsD0> D0CollisionsPerCollision = aod::jd0indices::collisionId;
  Preslice<CollisionsLc> LcCollisionsPerCollision = aod::jlcindices::collisionId;
  Preslice<CollisionsDielectron> DielectronCollisionsPerCollision = aod::jdielectronindices::collisionId;
  Preslice<CandidatesD0MCD> D0sPerCollision = aod::jd0indices::collisionId;
  Preslice<CandidatesLcMCD> LcsPerCollision = aod::jlcindices::collisionId;
  Preslice<CandidatesDielectronMCD> DielectronsPerCollision = aod::jdielectronindices::collisionId;

  std::vector<bool> collisionFlag;
  std::vector<bool> McCollisionFlag;
  std::vector<int> bcIndicies;

  uint32_t precisionPositionMask;
  uint32_t precisionMomentumMask;

  int eventSelection = -1;
  void init(InitContext&)
  {
    precisionPositionMask = 0xFFFFFC00; // 13 bits
    precisionMomentumMask = 0xFFFFFC00; // 13 bits  this is currently keept at 13 bits wihich gives roughly a resolution of 1/8000. This can be increased to 15 bits if really needed
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(config.eventSelectionForCounting));
  }

  bool acceptCollision(aod::JCollision const&)
  {
    return true;
  }

  void processCollisions(aod::JCollisions const& collisions)
  {
    collisionFlag.clear();
    collisionFlag.resize(collisions.size());
    std::fill(collisionFlag.begin(), collisionFlag.end(), false);
  }

  void processMcCollisions(aod::JMcCollisions const& Mccollisions)
  {
    McCollisionFlag.clear();
    McCollisionFlag.resize(Mccollisions.size());
    std::fill(McCollisionFlag.begin(), McCollisionFlag.end(), false);
  }

  template <typename T>
  void processJets(T& jets)
  {
    float jetPtMin = 0.0;
    if constexpr (std::is_same_v<std::decay_t<T>, aod::ChargedJets> || std::is_same_v<std::decay_t<T>, aod::ChargedMCDetectorLevelJets>) {
      jetPtMin = config.chargedJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::ChargedEventWiseSubtractedJets> || std::is_same_v<std::decay_t<T>, aod::ChargedMCDetectorLevelEventWiseSubtractedJets>) {
      jetPtMin = config.chargedEventWiseSubtractedJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::ChargedMCParticleLevelJets>) {
      jetPtMin = config.chargedMCPJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::NeutralJets> || std::is_same_v<std::decay_t<T>, aod::NeutralMCDetectorLevelJets>) {
      jetPtMin = config.neutralJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::NeutralMCParticleLevelJets>) {
      jetPtMin = config.neutralMCPJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::FullJets> || std::is_same_v<std::decay_t<T>, aod::FullMCDetectorLevelJets>) {
      jetPtMin = config.fullJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::FullMCParticleLevelJets>) {
      jetPtMin = config.fullMCPJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::D0ChargedJets> || std::is_same_v<std::decay_t<T>, aod::D0ChargedMCDetectorLevelJets>) {
      jetPtMin = config.chargedD0JetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::D0ChargedEventWiseSubtractedJets> || std::is_same_v<std::decay_t<T>, aod::D0ChargedMCDetectorLevelEventWiseSubtractedJets>) {
      jetPtMin = config.chargedEventWiseSubtractedD0JetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::D0ChargedMCParticleLevelJets>) {
      jetPtMin = config.chargedD0MCPJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::LcChargedJets> || std::is_same_v<std::decay_t<T>, aod::LcChargedMCDetectorLevelJets>) {
      jetPtMin = config.chargedLcJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::LcChargedEventWiseSubtractedJets> || std::is_same_v<std::decay_t<T>, aod::LcChargedMCDetectorLevelEventWiseSubtractedJets>) {
      jetPtMin = config.chargedEventWiseSubtractedLcJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::LcChargedMCParticleLevelJets>) {
      jetPtMin = config.chargedLcMCPJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::DielectronChargedJets> || std::is_same_v<std::decay_t<T>, aod::DielectronChargedMCDetectorLevelJets>) {
      jetPtMin = config.chargedDielectronJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::DielectronChargedEventWiseSubtractedJets> || std::is_same_v<std::decay_t<T>, aod::DielectronChargedMCDetectorLevelEventWiseSubtractedJets>) {
      jetPtMin = config.chargedEventWiseSubtractedDielectronJetPtMin;
    } else if constexpr (std::is_same_v<std::decay_t<T>, aod::DielectronChargedMCParticleLevelJets>) {
      jetPtMin = config.chargedDielectronMCPJetPtMin;
    } else {
      jetPtMin = 0.0;
    }
    for (const auto& jet : jets) {
      if (jet.pt() >= jetPtMin) {
        if constexpr (std::is_same_v<std::decay_t<T>, aod::ChargedMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::NeutralMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::FullMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::D0ChargedMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::LcChargedMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::BplusChargedMCParticleLevelJets> || std::is_same_v<std::decay_t<T>, aod::DielectronChargedMCParticleLevelJets>) {
          McCollisionFlag[jet.mcCollisionId()] = true;
        } else {
          collisionFlag[jet.collisionId()] = true;
        }
      }
    }
  }
  // Todo : Check memory consumption of having so many Process Switches
  PROCESS_SWITCH(JetDerivedDataWriter, processCollisions, "setup the writing for data and MCD", true);
  PROCESS_SWITCH(JetDerivedDataWriter, processMcCollisions, "setup the writing for MCP", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::ChargedJets>, processChargedJets, "process charged jets", true);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::ChargedEventWiseSubtractedJets>, processChargedEventWiseSubtractedJets, "process charged event-wise subtracted jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::ChargedMCDetectorLevelJets>, processChargedMCDJets, "process charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::ChargedMCDetectorLevelEventWiseSubtractedJets>, processChargedMCDetectorLevelEventWiseSubtractedJets, "process charged event-wise subtracted mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::ChargedMCParticleLevelJets>, processChargedMCPJets, "process charged mcp jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::NeutralJets>, processNeutralJets, "process neutral jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::NeutralMCDetectorLevelJets>, processNeutralMCDJets, "process neutral mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::NeutralMCParticleLevelJets>, processNeutralMCPJets, "process neutral mcp jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::FullJets>, processFullJets, "process full jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::FullMCDetectorLevelJets>, processFullMCDJets, "process full mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::FullMCParticleLevelJets>, processFullMCPJets, "process full mcp jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::D0ChargedJets>, processD0ChargedJets, "process D0 charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::D0ChargedEventWiseSubtractedJets>, processD0ChargedEventWiseSubtractedJets, "process D0 event-wise subtracted charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::D0ChargedMCDetectorLevelJets>, processD0ChargedMCDJets, "process D0 charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::D0ChargedMCDetectorLevelEventWiseSubtractedJets>, processD0ChargedMCDetectorLevelEventWiseSubtractedJets, "process D0 event-wise subtracted charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::D0ChargedMCParticleLevelJets>, processD0ChargedMCPJets, "process D0 charged mcp jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::LcChargedJets>, processLcChargedJets, "process Lc charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::LcChargedEventWiseSubtractedJets>, processLcChargedEventWiseSubtractedJets, "process Lc event-wise subtracted charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::LcChargedMCDetectorLevelJets>, processLcChargedMCDJets, "process Lc charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::LcChargedMCDetectorLevelEventWiseSubtractedJets>, processLcChargedMCDetectorLevelEventWiseSubtractedJets, "process Lc event-wise subtracted charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::LcChargedMCParticleLevelJets>, processLcChargedMCPJets, "process Lc charged mcp jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::DielectronChargedJets>, processDielectronChargedJets, "process Dielectron charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::DielectronChargedEventWiseSubtractedJets>, processDielectronChargedEventWiseSubtractedJets, "process Dielectron event-wise subtracted charged jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::DielectronChargedMCDetectorLevelJets>, processDielectronChargedMCDJets, "process Dielectron charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::DielectronChargedMCDetectorLevelEventWiseSubtractedJets>, processDielectronChargedMCDetectorLevelEventWiseSubtractedJets, "process Dielectron event-wise subtracted charged mcd jets", false);
  PROCESS_SWITCH_FULL(JetDerivedDataWriter, processJets<aod::DielectronChargedMCParticleLevelJets>, processDielectronChargedMCPJets, "process Dielectron charged mcp jets", false);

  void processDummyTable(aod::JDummys const&)
  {
    products.storedJDummysTable(1);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processDummyTable, "write out dummy output table", true);

  void processCollisionCounting(aod::JCollisions const& collisions, aod::CollisionCounts const& collisionCounts)
  {
    int readCollisionCounter = 0;
    int readSelectedCollisionCounter = 0;
    int writtenCollisionCounter = 0;
    for (const auto& collision : collisions) {
      readCollisionCounter++;
      if (jetderiveddatautilities::selectCollision(collision, eventSelection)) {
        readSelectedCollisionCounter++;
      }
      if (collisionFlag[collision.globalIndex()]) {
        writtenCollisionCounter++;
      }
    }
    std::vector<int> previousReadCounts;
    std::vector<int> previousReadSelectedCounts;
    std::vector<int> previousWrittenCounts;
    int iPreviousDataFrame = 0;
    for (const auto& collisionCount : collisionCounts) {
      auto readCollisionCounterSpan = collisionCount.readCounts();
      auto readSelectedCollisionCounterSpan = collisionCount.readSelectedCounts();
      auto writtenCollisionCounterSpan = collisionCount.writtenCounts();
      if (iPreviousDataFrame == 0) {
        std::copy(readCollisionCounterSpan.begin(), readCollisionCounterSpan.end(), std::back_inserter(previousReadCounts));
        std::copy(readSelectedCollisionCounterSpan.begin(), readSelectedCollisionCounterSpan.end(), std::back_inserter(previousReadSelectedCounts));
        std::copy(writtenCollisionCounterSpan.begin(), writtenCollisionCounterSpan.end(), std::back_inserter(previousWrittenCounts));
      } else {
        for (unsigned int i = 0; i < previousReadCounts.size(); i++) {
          previousReadCounts[i] += readCollisionCounterSpan[i];
          previousReadSelectedCounts[i] += readSelectedCollisionCounterSpan[i];
          previousWrittenCounts[i] += writtenCollisionCounterSpan[i];
        }
      }
      iPreviousDataFrame++;
    }
    previousReadCounts.push_back(readCollisionCounter);
    previousReadSelectedCounts.push_back(readSelectedCollisionCounter);
    previousWrittenCounts.push_back(writtenCollisionCounter);
    products.storedCollisionCountsTable(previousReadCounts, previousReadSelectedCounts, previousWrittenCounts);
  }
  PROCESS_SWITCH(JetDerivedDataWriter, processCollisionCounting, "write out collision counting output table", false);

  void processData(soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::JCollisionBCs, aod::JChTrigSels, aod::JFullTrigSels, aod::JChHFTrigSels>::iterator const& collision, soa::Join<aod::JBCs, aod::JBCPIs> const&, soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs> const& tracks, soa::Join<aod::JClusters, aod::JClusterPIs, aod::JClusterTracks> const& clusters, CollisionsD0 const& D0Collisions, CandidatesD0Data const& D0s, CollisionsLc const& LcCollisions, CandidatesLcData const& Lcs, CollisionsDielectron const& DielectronCollisions, CandidatesDielectronData const& Dielectrons)
  {
    std::map<int32_t, int32_t> bcMapping;
    std::map<int32_t, int32_t> trackMapping;

    if (collisionFlag[collision.globalIndex()]) {
      if (config.saveBCsTable) {
        auto bc = collision.bc_as<soa::Join<aod::JBCs, aod::JBCPIs>>();
        if (std::find(bcIndicies.begin(), bcIndicies.end(), bc.globalIndex()) == bcIndicies.end()) {
          products.storedJBCsTable(bc.runNumber(), bc.globalBC(), bc.timestamp());
          products.storedJBCParentIndexTable(bc.bcId());
          bcIndicies.push_back(bc.globalIndex());
          bcMapping.insert(std::make_pair(bc.globalIndex(), products.storedJBCsTable.lastIndex()));
        }
      }

      products.storedJCollisionsTable(collision.posX(), collision.posY(), collision.posZ(), collision.multiplicity(), collision.centrality(), collision.eventSel(), collision.alias_raw());
      products.storedJCollisionsParentIndexTable(collision.collisionId());
      if (config.saveBCsTable) {
        int32_t storedBCID = -1;
        auto JBCIndex = bcMapping.find(collision.bcId());
        if (JBCIndex != bcMapping.end()) {
          storedBCID = JBCIndex->second;
        }
        products.storedJCollisionsBunchCrossingIndexTable(storedBCID);
      }
      products.storedJChargedTriggerSelsTable(collision.chargedTriggerSel());
      products.storedJFullTriggerSelsTable(collision.fullTriggerSel());
      products.storedJChargedHFTriggerSelsTable(collision.chargedHFTriggerSel());

      for (const auto& track : tracks) {
        if (config.performTrackSelection && !(track.trackSel() & ~(1 << jetderiveddatautilities::JTrackSel::trackSign))) { // skips tracks that pass no selections. This might cause a problem with tracks matched with clusters. We should generate a track selection purely for cluster matched tracks so that they are kept
          continue;
        }
        products.storedJTracksTable(products.storedJCollisionsTable.lastIndex(), o2::math_utils::detail::truncateFloatFraction(track.pt(), precisionMomentumMask), o2::math_utils::detail::truncateFloatFraction(track.eta(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.phi(), precisionPositionMask), track.trackSel());
        products.storedJTracksExtraTable(o2::math_utils::detail::truncateFloatFraction(track.dcaXY(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.dcaZ(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.sigma1Pt(), precisionMomentumMask));
        products.storedJTracksParentIndexTable(track.trackId());
        trackMapping.insert(std::make_pair(track.globalIndex(), products.storedJTracksTable.lastIndex()));
      }
      if (config.saveClustersTable) {
        for (const auto& cluster : clusters) {
          products.storedJClustersTable(products.storedJCollisionsTable.lastIndex(), cluster.id(), cluster.energy(), cluster.coreEnergy(), cluster.rawEnergy(),
                                        cluster.eta(), cluster.phi(), cluster.m02(), cluster.m20(), cluster.nCells(), cluster.time(), cluster.isExotic(), cluster.distanceToBadChannel(),
                                        cluster.nlm(), cluster.definition(), cluster.leadingCellEnergy(), cluster.subleadingCellEnergy(), cluster.leadingCellNumber(), cluster.subleadingCellNumber());
          products.storedJClustersParentIndexTable(cluster.clusterId());

          std::vector<int> clusterStoredJTrackIDs;
          for (const auto& clusterTrack : cluster.matchedTracks_as<soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs>>()) {
            auto JtrackIndex = trackMapping.find(clusterTrack.globalIndex());
            if (JtrackIndex != trackMapping.end()) {
              clusterStoredJTrackIDs.push_back(JtrackIndex->second);
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
  PROCESS_SWITCH(JetDerivedDataWriter, processData, "write out data output tables", false);

  void processMC(soa::Join<aod::JMcCollisions, aod::JMcCollisionPIs> const& mcCollisions, soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::JCollisionBCs, aod::JChTrigSels, aod::JFullTrigSels, aod::JChHFTrigSels, aod::JMcCollisionLbs> const& collisions, soa::Join<aod::JBCs, aod::JBCPIs> const&, soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs, aod::JMcTrackLbs> const& tracks, soa::Join<aod::JClusters, aod::JClusterPIs, aod::JClusterTracks, aod::JMcClusterLbs> const& clusters, soa::Join<aod::JMcParticles, aod::JMcParticlePIs> const& particles, CollisionsD0 const& D0Collisions, CandidatesD0MCD const& D0s, soa::Join<McCollisionsD0, aod::HfD0McRCollIds> const& D0McCollisions, CandidatesD0MCP const& D0Particles, CollisionsLc const& LcCollisions, CandidatesLcMCD const& Lcs, soa::Join<McCollisionsLc, aod::Hf3PMcRCollIds> const& LcMcCollisions, CandidatesLcMCP const& LcParticles, CollisionsDielectron const& DielectronCollisions, CandidatesDielectronMCD const& Dielectrons, McCollisionsDielectron const& DielectronMcCollisions, CandidatesDielectronMCP const& DielectronParticles)
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
            std::vector<int> DielectronMothersId;
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
              products.storedJBCsTable(bc.runNumber(), bc.globalBC(), bc.timestamp());
              products.storedJBCParentIndexTable(bc.bcId());
              bcIndicies.push_back(bc.globalIndex());
              bcMapping.insert(std::make_pair(bc.globalIndex(), products.storedJBCsTable.lastIndex()));
            }
          }

          products.storedJCollisionsTable(collision.posX(), collision.posY(), collision.posZ(), collision.multiplicity(), collision.centrality(), collision.eventSel(), collision.alias_raw());
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
          products.storedJChargedTriggerSelsTable(collision.chargedTriggerSel());
          products.storedJFullTriggerSelsTable(collision.fullTriggerSel());
          products.storedJChargedHFTriggerSelsTable(collision.chargedHFTriggerSel());

          const auto tracksPerCollision = tracks.sliceBy(TracksPerCollision, collision.globalIndex());
          for (const auto& track : tracksPerCollision) {
            if (config.performTrackSelection && !(track.trackSel() & ~(1 << jetderiveddatautilities::JTrackSel::trackSign))) { // skips tracks that pass no selections. This might cause a problem with tracks matched with clusters. We should generate a track selection purely for cluster matched tracks so that they are kept
              continue;
            }
            products.storedJTracksTable(products.storedJCollisionsTable.lastIndex(), o2::math_utils::detail::truncateFloatFraction(track.pt(), precisionMomentumMask), o2::math_utils::detail::truncateFloatFraction(track.eta(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.phi(), precisionPositionMask), track.trackSel());
            products.storedJTracksExtraTable(o2::math_utils::detail::truncateFloatFraction(track.dcaXY(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.dcaZ(), precisionPositionMask), o2::math_utils::detail::truncateFloatFraction(track.sigma1Pt(), precisionMomentumMask));
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

              std::vector<int> clusterStoredJTrackIDs;
              for (const auto& clusterTrack : cluster.matchedTracks_as<soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs>>()) {
                auto JtrackIndex = trackMapping.find(clusterTrack.globalIndex());
                if (JtrackIndex != trackMapping.end()) {
                  clusterStoredJTrackIDs.push_back(JtrackIndex->second);
                }
              }
              products.storedJClustersMatchedTracksTable(clusterStoredJTrackIDs);

              std::vector<int> clusterStoredJParticleIDs;
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
            for (auto const& d0CollisionPerMcCollision : d0McCollisionPerMcCollision.template hfCollBases_as<CollisionsD0>()) {
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
            for (auto const& lcCollisionPerMcCollision : lcMcCollisionPerMcCollision.template hfCollBases_as<CollisionsLc>()) {
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
  PROCESS_SWITCH(JetDerivedDataWriter, processMC, "write out data output tables for mc", false);

  void processMCP(soa::Join<aod::JMcCollisions, aod::JMcCollisionPIs> const& mcCollisions, soa::Join<aod::JMcParticles, aod::JMcParticlePIs> const& particles, McCollisionsD0 const& D0McCollisions, CandidatesD0MCP const& D0Particles, McCollisionsLc const& LcMcCollisions, CandidatesLcMCP const& LcParticles, McCollisionsDielectron const& DielectronMcCollisions, CandidatesDielectronMCP const& DielectronParticles)
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
            std::vector<int> DielectronMothersId;
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
  PROCESS_SWITCH(JetDerivedDataWriter, processMCP, "write out data output tables for mcp", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetDerivedDataWriter>(cfgc, TaskName{"jet-deriveddata-writer"}));

  return WorkflowSpec{tasks};
}
