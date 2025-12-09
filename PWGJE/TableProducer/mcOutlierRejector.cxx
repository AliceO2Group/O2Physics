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

// Task to produce a table joinable to the jcollision table which contains an outlier rejector
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct McOutlierRejectorTask {
  Produces<aod::JCollisionOutliers> collisionOutliers;
  Produces<aod::JMcCollisionOutliers> mcCollisionOutliers;
  Preslice<soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs>> perColParticle = aod::jmccollision::mcCollisionId;

  Configurable<bool> checkmcCollisionForCollision{"checkmcCollisionForCollision", true, "additionally reject collision based on mcCollision"};
  Configurable<float> ptHatMax{"ptHatMax", 4.0, "maximum factor of pt hat the leading jet in the event is allowed"};
  Configurable<float> ptTrackMaxMinBias{"ptTrackMaxMinBias", 20.0, "maximum pt for track originating from minimum bias event"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  std::vector<bool> collisionFlag;
  std::vector<bool> mcCollisionFlag;

  int trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

  void processSetupCollisionSelection(aod::JCollisions const& collisions)
  {
    collisionFlag.clear();
    collisionFlag.resize(collisions.size(), false);
  }
  PROCESS_SWITCH(McOutlierRejectorTask, processSetupCollisionSelection, "Setup collision processing", true);

  void processSetupMcCollisionSelection(aod::JMcCollisions const& mcCollisions)
  {
    mcCollisionFlag.clear();
    mcCollisionFlag.resize(mcCollisions.size(), false);
  }
  PROCESS_SWITCH(McOutlierRejectorTask, processSetupMcCollisionSelection, "Setup MC Collision processing", true);

  template <typename T>
  void collisionSelection(int32_t collisionIndex, int32_t mcCollisionId, T const& selectionObjects, float ptHard, std::vector<bool>& flagArray, std::optional<std::reference_wrapper<const soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs>>> mcCollisionsOpt = std::nullopt)
  {
    if (selectionObjects.size() != 0) {
      float selectionObjectPt = 0.0;
      if constexpr (std::is_same_v<std::decay_t<T>, aod::JetTracksMCD> || std::is_same_v<std::decay_t<T>, aod::JetParticles>) {
        for (auto selectionObject : selectionObjects) {
          selectionObjectPt = selectionObject.pt();
          // may be slow - could save only MC particle then check difference only for tracks IDd as outliers?
          if constexpr (std::is_same_v<std::decay_t<T>, aod::JetTracksMCD>) { // tracks
            if (!jetderiveddatautilities::selectTrack(selectionObject, trackSelection)) {
              continue;
            }
            auto& mcCollisions = mcCollisionsOpt.value().get();
            auto mcParticle = selectionObject.template mcParticle_as<soa::Join<aod::JetParticles, aod::JMcParticlePIs>>();
            auto mcCollision = mcCollisions.sliceBy(perColParticle, mcParticle.mcCollisionId());
            int subGenID = mcCollision.begin().subGeneratorId();
            int diffCollisionID = mcParticle.mcCollisionId() - mcCollisionId;
            if (diffCollisionID != 0 &&
                selectionObjectPt > ptHatMax * ptHard) {
              if (subGenID != jetderiveddatautilities::JCollisionSubGeneratorId::mbGap && selectionObjectPt > ptHatMax * ptHard) {
                flagArray[collisionIndex] = true;
              }
              if (subGenID == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap && selectionObjectPt > ptTrackMaxMinBias) {
                flagArray[collisionIndex] = true;
              }
            }
          } else { // particles
            if (selectionObjectPt > ptHatMax * ptHard) {
              flagArray[collisionIndex] = true;
            }
          }
        }
      } else { // jets
        selectionObjectPt = selectionObjects.iteratorAt(0).pt();
        if (selectionObjectPt > ptHatMax * ptHard) {
          flagArray[collisionIndex] = true; // Currently if running multiple different jet finders, then a single type of jet can veto an event for others. Decide if this is the best way
        }
      }
    }
  }

  template <typename T>
  void processSelectionObjects(aod::JetCollisionMCD const& collision, T const& selectionObjects, aod::JetMcCollisions const&)
  {
    auto mcCollision = collision.mcCollision_as<aod::JetMcCollisions>();
    collisionSelection(collision.globalIndex(), 1., selectionObjects, mcCollision.ptHard(), collisionFlag);
  }

  template <typename T>
  void processSelectionObjectsTracks(aod::JetCollisionMCD const& collision, T const& selectionObjects, soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs> const& mcCollisions, soa::Join<aod::JetParticles, aod::JMcParticlePIs> const&)
  {
    auto mcCollision = collision.mcCollision_as<soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs>>();
    collisionSelection(collision.globalIndex(), collision.mcCollisionId(), selectionObjects, mcCollision.ptHard(), collisionFlag, mcCollisions);
  }

  template <typename T>
  void processSelectionMcObjects(aod::JetMcCollision const& mcCollision, T const& selectionMcObjects)
  {
    collisionSelection(mcCollision.globalIndex(), 1., selectionMcObjects, mcCollision.ptHard(), mcCollisionFlag);
  }

  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionObjects<aod::ChargedMCDetectorLevelJets>, processSelectingChargedMCDetectorLevelJets, "process mc detector level charged jets", true);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionObjects<aod::NeutralMCDetectorLevelJets>, processSelectingNeutralMCDetectorLevelJets, "process mc detector level neutral jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionObjects<aod::FullMCDetectorLevelJets>, processSelectingFullMCDetectorLevelJets, "process mc detector level full jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionObjects<aod::D0ChargedMCDetectorLevelJets>, processSelectingD0ChargedMCDetectorLevelJets, "process mc detector level D0 charged jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionObjects<aod::DplusChargedMCDetectorLevelJets>, processSelectingDplusChargedMCDetectorLevelJets, "process mc detector level Dplus charged jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionObjects<aod::DsChargedMCDetectorLevelJets>, processSelectingDsChargedMCDetectorLevelJets, "process mc detector level Ds charged jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionObjects<aod::DstarChargedMCDetectorLevelJets>, processSelectingDstarChargedMCDetectorLevelJets, "process mc detector level Dstar charged jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionObjects<aod::LcChargedMCDetectorLevelJets>, processSelectingLcChargedMCDetectorLevelJets, "process mc detector level Lc charged jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionObjects<aod::B0ChargedMCDetectorLevelJets>, processSelectingB0ChargedMCDetectorLevelJets, "process mc detector level B0 charged jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionObjects<aod::BplusChargedMCDetectorLevelJets>, processSelectingBplusChargedMCDetectorLevelJets, "process mc detector level Bplus charged jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionObjects<aod::XicToXiPiPiChargedMCDetectorLevelJets>, processSelectingXicToXiPiPiChargedMCDetectorLevelJets, "process mc detector level XicToXiPiPi charged jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionObjects<aod::DielectronChargedMCDetectorLevelJets>, processSelectingDielectronChargedMCDetectorLevelJets, "process mc detector level Dielectron charged jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionObjectsTracks<aod::JetTracksMCD>, processSelectingTracks, "process tracks", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionMcObjects<aod::ChargedMCParticleLevelJets>, processSelectingChargedMCParticleLevelJets, "process mc particle level charged jets", true);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionMcObjects<aod::NeutralMCParticleLevelJets>, processSelectingNeutralMCParticleLevelJets, "process mc particle level neutral jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionMcObjects<aod::FullMCParticleLevelJets>, processSelectingFullMCParticleLevelJets, "process mc particle level full jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionMcObjects<aod::D0ChargedMCParticleLevelJets>, processSelectingD0ChargedMCParticleLevelJets, "process mc particle level D0 charged jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionMcObjects<aod::DplusChargedMCParticleLevelJets>, processSelectingDplusChargedMCParticleLevelJets, "process mc particle level Dplus charged jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionMcObjects<aod::DsChargedMCParticleLevelJets>, processSelectingDsChargedMCParticleLevelJets, "process mc particle level Ds charged jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionMcObjects<aod::DstarChargedMCParticleLevelJets>, processSelectingDstarChargedMCParticleLevelJets, "process mc particle level Dstar charged jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionMcObjects<aod::LcChargedMCParticleLevelJets>, processSelectingLcChargedMCParticleLevelJets, "process mc particle level Lc charged jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionMcObjects<aod::B0ChargedMCParticleLevelJets>, processSelectingB0ChargedMCParticleLevelJets, "process mc particle level B0 charged jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionMcObjects<aod::BplusChargedMCParticleLevelJets>, processSelectingBplusChargedMCParticleLevelJets, "process mc particle level Bplus charged jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionMcObjects<aod::XicToXiPiPiChargedMCParticleLevelJets>, processSelectingXicToXiPiPiChargedMCParticleLevelJets, "process mc particle level XicToXiPiPi charged jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionMcObjects<aod::DielectronChargedMCParticleLevelJets>, processSelectingDielectronChargedMCParticleLevelJets, "process mc particle level Dielectron charged jets", false);
  PROCESS_SWITCH_FULL(McOutlierRejectorTask, processSelectionMcObjects<aod::JetParticles>, processSelectingParticles, "process mc particles", false);

  void processStoreCollisionDecision(aod::JetCollisionMCD const& collision)
  {
    bool rejectCollision = collisionFlag[collision.globalIndex()];
    if (!rejectCollision && checkmcCollisionForCollision && collision.has_mcCollision()) {
      rejectCollision = mcCollisionFlag[collision.mcCollisionId()];
    }
    collisionOutliers(rejectCollision);
  }
  PROCESS_SWITCH(McOutlierRejectorTask, processStoreCollisionDecision, "write out decision of rejecting collision", true);

  void processStoreMcCollisionDecision(aod::JetMcCollision const& mcCollision)
  {
    mcCollisionOutliers(mcCollisionFlag[mcCollision.globalIndex()]);
  }
  PROCESS_SWITCH(McOutlierRejectorTask, processStoreMcCollisionDecision, "write out decision of rejecting mcCollision", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<McOutlierRejectorTask>(cfgc, TaskName{"mc-outlier-rejector"})}; }
