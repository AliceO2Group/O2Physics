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

/// \file luminosityproducer.cxx
/// \brief Task to produce tables needed for normalisation and luminosity calculation
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "JetDerivedDataUtilities.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/CCDB/EventSelectionParams.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct LuminosityProducer {

  Produces<aod::StoredBCCounts> bcCountsTable;
  Produces<aod::StoredCollisionCounts> collisionCountsTable;

  Configurable<float> vertexZCutForCounting{"vertexZCutForCounting", 10.0, "choose z-vertex cut for collision counter"};
  Configurable<std::string> customEventSelections{"customEventSelections", "sel8", "choose custom event selection to be added"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", true, "decide to run over MB gap events or not"};
  Configurable<bool> applyRCTSelections{"applyRCTSelections", true, "decide to apply RCT selections"};

  void init(InitContext&)
  {
  }

  void processBCCountingNonDerived(aod::JBCs const& bcs)
  {
    int bcCounter = 0;
    int bcWithTVXCounter = 0;
    int bcWithTVXAndNoTFBCounter = 0;
    int bcWithTVXAndNoTFBAndNoITSROFBCounter = 0;
    for (const auto& bc : bcs) {
      bcCounter++;
      if (bc.selection_bit(aod::evsel::EventSelectionFlags::kIsTriggerTVX)) {
        bcWithTVXCounter++;
        if (bc.selection_bit(aod::evsel::EventSelectionFlags::kNoTimeFrameBorder)) {
          bcWithTVXAndNoTFBCounter++;
          if (bc.selection_bit(aod::evsel::EventSelectionFlags::kNoITSROFrameBorder)) {
            bcWithTVXAndNoTFBAndNoITSROFBCounter++;
          }
        }
      }
    }
    bcCountsTable(bcCounter, bcWithTVXCounter, bcWithTVXAndNoTFBCounter, bcWithTVXAndNoTFBAndNoITSROFBCounter);
  }
  PROCESS_SWITCH(LuminosityProducer, processBCCountingNonDerived, "write out bc counting output table for running on full AO2D", true);

  void processBCCountingDerived(aod::BCCounts const& bcCounts)
  {
    int bcCounter = 0;
    int bcWithTVXCounter = 0;
    int bcWithTVXAndNoTFBCounter = 0;
    int bcWithTVXAndNoTFBAndNoITSROFBCounter = 0;

    for (const auto& bcCount : bcCounts) {
      bcCounter += bcCount.counts();
      bcWithTVXCounter += bcCount.countsWithTVX();
      bcWithTVXAndNoTFBCounter += bcCount.countsWithTVXAndNoTFB();
      bcWithTVXAndNoTFBAndNoITSROFBCounter += bcCount.countsWithTVXAndNoTFBAndNoITSROFB();
    }
    bcCountsTable(bcCounter, bcWithTVXCounter, bcWithTVXAndNoTFBCounter, bcWithTVXAndNoTFBAndNoITSROFBCounter);
  }
  PROCESS_SWITCH(LuminosityProducer, processBCCountingDerived, "write out bc counting output table for running on derived data", false);

  void processCollisionCountingNonDerived(aod::JetCollisions const& collisions)
  {

    int collisionCounter = 0;
    int collisionWithTVXCounter = 0;
    int collisionWithTVXAndZVertexAndSel8Counter = 0;
    int collisionWithTVXAndZVertexAndSel8FullCounter = 0;
    int collisionWithTVXAndZVertexAndSel8FullPbPbCounter = 0;
    int collisionWithTVXAndZVertexAndSelMCCounter = 0;
    int collisionWithTVXAndZVertexAndSelMCFullCounter = 0;
    int collisionWithTVXAndZVertexAndSelMCFullPbPbCounter = 0;
    int collisionWithTVXAndZVertexAndSelUnanchoredMCCounter = 0;
    int collisionWithTVXAndZVertexAndSelTVXCounter = 0; // redundant but we keep it
    int collisionWithTVXAndZVertexAndSel7Counter = 0;
    int collisionWithTVXAndZVertexAndSel7KINT7Counter = 0;
    int collisionWithCustomCounter = 0;
    for (const auto& collision : collisions) {
      collisionCounter++;
      if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("TVX"), skipMBGapEvents, false)) { // asuumes all selections include the TVX trigger but for this step does not include the rct flags
        collisionWithTVXCounter++;
        if (std::abs(collision.posZ()) > vertexZCutForCounting) {
          continue;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("sel8"), skipMBGapEvents, applyRCTSelections)) {
          collisionWithTVXAndZVertexAndSel8Counter++;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("TVX"), skipMBGapEvents, applyRCTSelections)) {
          collisionWithTVXAndZVertexAndSelTVXCounter++;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("sel7"), skipMBGapEvents, applyRCTSelections)) {
          collisionWithTVXAndZVertexAndSel7Counter++;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("sel7KINT7"), skipMBGapEvents, applyRCTSelections)) {
          collisionWithTVXAndZVertexAndSel7KINT7Counter++;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("sel8Full"), skipMBGapEvents, applyRCTSelections)) {
          collisionWithTVXAndZVertexAndSel8FullCounter++;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("sel8FullPbPb"), skipMBGapEvents, applyRCTSelections)) {
          collisionWithTVXAndZVertexAndSel8FullPbPbCounter++;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("selMC"), skipMBGapEvents, applyRCTSelections)) {
          collisionWithTVXAndZVertexAndSelMCCounter++;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("selMCFull"), skipMBGapEvents, applyRCTSelections)) {
          collisionWithTVXAndZVertexAndSelMCFullCounter++;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("selMCFullPbPb"), skipMBGapEvents, applyRCTSelections)) {
          collisionWithTVXAndZVertexAndSelMCFullPbPbCounter++;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("selUnanchoredMC"), skipMBGapEvents, applyRCTSelections)) {
          collisionWithTVXAndZVertexAndSelUnanchoredMCCounter++;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(customEventSelections)), skipMBGapEvents, applyRCTSelections)) {
          collisionWithCustomCounter++;
        }
      }
    }

    collisionCountsTable(collisionCounter, collisionWithTVXCounter, collisionWithTVXAndZVertexAndSel8Counter, collisionWithTVXAndZVertexAndSel8FullCounter, collisionWithTVXAndZVertexAndSel8FullPbPbCounter, collisionWithTVXAndZVertexAndSelMCCounter, collisionWithTVXAndZVertexAndSelMCFullCounter, collisionWithTVXAndZVertexAndSelMCFullPbPbCounter, collisionWithTVXAndZVertexAndSelUnanchoredMCCounter, collisionWithTVXAndZVertexAndSelTVXCounter, collisionWithTVXAndZVertexAndSel7Counter, collisionWithTVXAndZVertexAndSel7KINT7Counter, collisionWithCustomCounter);
  }
  PROCESS_SWITCH(LuminosityProducer, processCollisionCountingNonDerived, "write out collision counting output table for running on full AO2D", true);

  void processStoreCollisionCountingDerived(aod::CollisionCounts const& collisionCounts)
  {
    int collisionCounter = 0;
    int collisionWithTVXCounter = 0;
    int collisionWithTVXAndZVertexAndSel8Counter = 0;
    int collisionWithTVXAndZVertexAndSel8FullCounter = 0;
    int collisionWithTVXAndZVertexAndSel8FullPbPbCounter = 0;
    int collisionWithTVXAndZVertexAndSelMCCounter = 0;
    int collisionWithTVXAndZVertexAndSelMCFullCounter = 0;
    int collisionWithTVXAndZVertexAndSelMCFullPbPbCounter = 0;
    int collisionWithTVXAndZVertexAndSelUnanchoredMCCounter = 0;
    int collisionWithTVXAndZVertexAndSelTVXCounter = 0; // redundant but we keep it
    int collisionWithTVXAndZVertexAndSel7Counter = 0;
    int collisionWithTVXAndZVertexAndSel7KINT7Counter = 0;
    int collisionWithCustomCounter = 0;
    for (const auto& collisionCount : collisionCounts) {
      collisionCounter += collisionCount.counts();
      collisionWithTVXCounter += collisionCount.countsWithTVX();
      collisionWithTVXAndZVertexAndSel8Counter += collisionCount.countsWithTVXAndZVertexAndSel8();
      collisionWithTVXAndZVertexAndSel8FullCounter += collisionCount.countsWithTVXAndZVertexAndSel8Full();
      collisionWithTVXAndZVertexAndSel8FullPbPbCounter += collisionCount.countsWithTVXAndZVertexAndSel8FullPbPb();
      collisionWithTVXAndZVertexAndSelMCCounter += collisionCount.countsWithTVXAndZVertexAndSelMC();
      collisionWithTVXAndZVertexAndSelMCFullCounter += collisionCount.countsWithTVXAndZVertexAndSelMCFull();
      collisionWithTVXAndZVertexAndSelMCFullPbPbCounter += collisionCount.countsWithTVXAndZVertexAndSelMCFullPbPb();
      collisionWithTVXAndZVertexAndSelUnanchoredMCCounter += collisionCount.countsWithTVXAndZVertexAndSelUnanchoredMC();
      collisionWithTVXAndZVertexAndSelTVXCounter += collisionCount.countsWithTVXAndZVertexAndSelTVX();
      collisionWithTVXAndZVertexAndSel7Counter += collisionCount.countsWithTVXAndZVertexAndSel7();
      collisionWithTVXAndZVertexAndSel7KINT7Counter += collisionCount.countsWithTVXAndZVertexAndSel7KINT7();
      collisionWithCustomCounter += collisionCount.countsWithCustomSelection();
    }

    collisionCountsTable(collisionCounter, collisionWithTVXCounter, collisionWithTVXAndZVertexAndSel8Counter, collisionWithTVXAndZVertexAndSel8FullCounter, collisionWithTVXAndZVertexAndSel8FullPbPbCounter, collisionWithTVXAndZVertexAndSelMCCounter, collisionWithTVXAndZVertexAndSelMCFullCounter, collisionWithTVXAndZVertexAndSelMCFullPbPbCounter, collisionWithTVXAndZVertexAndSelUnanchoredMCCounter, collisionWithTVXAndZVertexAndSelTVXCounter, collisionWithTVXAndZVertexAndSel7Counter, collisionWithTVXAndZVertexAndSel7KINT7Counter, collisionWithCustomCounter);
  }
  PROCESS_SWITCH(LuminosityProducer, processStoreCollisionCountingDerived, "write out collision counting output table for running on derived data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<LuminosityProducer>(cfgc, TaskName{"jet-luminosity-producer"}));

  return WorkflowSpec{tasks};
}
