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

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct LuminosityProducer {

  Produces<aod::StoredBCCounts> storedBCCountsTable;
  Produces<aod::StoredCollisionCounts> storedCollisionCountsTable;

  Configurable<float> vertexZCutForCounting{"vertexZCutForCounting", 10.0, "choose z-vertex cut for collision counter"};
  Configurable<std::string> customEventSelections{"customEventSelections", "sel8", "choose custom event selection to be added"};

  void init(InitContext&)
  {
  }

  void processStoreBCCounting(aod::JBCs const& bcs, aod::BCCounts const& bcCounts)
  {
    int readBCCounter = 0;
    int readBCWithTVXCounter = 0;
    int readBCWithTVXAndNoTFBCounter = 0;
    int readBCWithTVXAndNoTFBAndNoITSROFBCounter = 0;
    for (const auto& bc : bcs) {
      readBCCounter++;
      if (bc.selection_bit(aod::evsel::EventSelectionFlags::kIsTriggerTVX)) {
        readBCWithTVXCounter++;
        if (bc.selection_bit(aod::evsel::EventSelectionFlags::kNoTimeFrameBorder)) {
          readBCWithTVXAndNoTFBCounter++;
          if (bc.selection_bit(aod::evsel::EventSelectionFlags::kNoITSROFrameBorder)) {
            readBCWithTVXAndNoTFBAndNoITSROFBCounter++;
          }
        }
      }
    }
    std::vector<int> previousReadCounts;
    std::vector<int> previousReadCountsWithTVX;
    std::vector<int> previousReadCountsWithTVXAndNoTFB;
    std::vector<int> previousReadCountsWithTVXAndNoTFBAndNoITSROFB;
    int iPreviousDataFrame = 0;
    for (const auto& bcCount : bcCounts) {
      auto readBCCounterSpan = bcCount.readCounts();
      auto readBCWithTVXCounterSpan = bcCount.readCountsWithTVX();
      auto readBCWithTVXAndNoTFBCounterSpan = bcCount.readCountsWithTVXAndNoTFB();
      auto readBCWithTVXAndNoTFBAndNoITSROFBCounterSpan = bcCount.readCountsWithTVXAndNoTFBAndNoITSROFB();
      if (iPreviousDataFrame == 0) {
        std::copy(readBCCounterSpan.begin(), readBCCounterSpan.end(), std::back_inserter(previousReadCounts));
        std::copy(readBCWithTVXCounterSpan.begin(), readBCWithTVXCounterSpan.end(), std::back_inserter(previousReadCountsWithTVX));
        std::copy(readBCWithTVXAndNoTFBCounterSpan.begin(), readBCWithTVXAndNoTFBCounterSpan.end(), std::back_inserter(previousReadCountsWithTVXAndNoTFB));
        std::copy(readBCWithTVXAndNoTFBAndNoITSROFBCounterSpan.begin(), readBCWithTVXAndNoTFBAndNoITSROFBCounterSpan.end(), std::back_inserter(previousReadCountsWithTVXAndNoTFBAndNoITSROFB));
      } else {
        for (unsigned int i = 0; i < previousReadCounts.size(); i++) { // in principle we only care about the first element, but might be interesting information to keep
          previousReadCounts[i] += readBCCounterSpan[i];
          previousReadCountsWithTVX[i] += readBCWithTVXCounterSpan[i];
          previousReadCountsWithTVXAndNoTFB[i] += readBCWithTVXAndNoTFBCounterSpan[i];
          previousReadCountsWithTVXAndNoTFBAndNoITSROFB[i] += readBCWithTVXAndNoTFBAndNoITSROFBCounterSpan[i];
        }
      }
      iPreviousDataFrame++;
    }
    previousReadCounts.push_back(readBCCounter);
    previousReadCountsWithTVX.push_back(readBCWithTVXCounter);
    previousReadCountsWithTVXAndNoTFB.push_back(readBCWithTVXAndNoTFBCounter);
    previousReadCountsWithTVXAndNoTFBAndNoITSROFB.push_back(readBCWithTVXAndNoTFBAndNoITSROFBCounter);
    storedBCCountsTable(previousReadCounts, previousReadCountsWithTVX, previousReadCountsWithTVXAndNoTFB, previousReadCountsWithTVXAndNoTFBAndNoITSROFB);
  }
  PROCESS_SWITCH(LuminosityProducer, processStoreBCCounting, "write out bc counting output table", true);

  void processStoreCollisionCounting(aod::JetCollisions const& collisions, aod::CollisionCounts const& collisionCounts)
  {
    int readCollisionCounter = 0;
    int readCollisionWithTVXCounter = 0;
    int readCollisionWithTVXAndZVertexAndSel8Counter = 0;
    int readCollisionWithTVXAndZVertexAndSel8FullCounter = 0;
    int readCollisionWithTVXAndZVertexAndSel8FullPbPbCounter = 0;
    int readCollisionWithTVXAndZVertexAndSelMCCounter = 0;
    int readCollisionWithTVXAndZVertexAndSelMCFullCounter = 0;
    int readCollisionWithTVXAndZVertexAndSelMCFullPbPbCounter = 0;
    int readCollisionWithTVXAndZVertexAndSelUnanchoredMCCounter = 0;
    int readCollisionWithTVXAndZVertexAndSelTVXCounter = 0; // redundant but we keep it
    int readCollisionWithTVXAndZVertexAndSel7Counter = 0;
    int readCollisionWithTVXAndZVertexAndSel7KINT7Counter = 0;
    int readCollisionWithCustomCounter = 0;
    for (const auto& collision : collisions) {
      readCollisionCounter++;
      if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("TVX"))) { // asuumes all selections include the TVX trigger
        readCollisionWithTVXCounter++;
        if (std::abs(collision.posZ()) > vertexZCutForCounting) {
          continue;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("sel8"))) {
          readCollisionWithTVXAndZVertexAndSel8Counter++;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("TVX"))) {
          readCollisionWithTVXAndZVertexAndSelTVXCounter++;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("sel7"))) {
          readCollisionWithTVXAndZVertexAndSel7Counter++;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("sel7KINT7"))) {
          readCollisionWithTVXAndZVertexAndSel7KINT7Counter++;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("sel8Full"))) {
          readCollisionWithTVXAndZVertexAndSel8FullCounter++;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("sel8FullPbPb"))) {
          readCollisionWithTVXAndZVertexAndSel8FullPbPbCounter++;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("selMC"))) {
          readCollisionWithTVXAndZVertexAndSelMCCounter++;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("selMCFull"))) {
          readCollisionWithTVXAndZVertexAndSelMCFullCounter++;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("selMCFullPbPb"))) {
          readCollisionWithTVXAndZVertexAndSelMCFullPbPbCounter++;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits("selUnanchoredMC"))) {
          readCollisionWithTVXAndZVertexAndSelUnanchoredMCCounter++;
        }
        if (jetderiveddatautilities::selectCollision(collision, jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(customEventSelections)))) {
          readCollisionWithCustomCounter++;
        }
      }
    }
    std::vector<int> previousReadCounts;
    std::vector<int> previousReadCountsWithTVX;
    std::vector<int> previousReadCountsWithTVXAndZVertexAndSel8;
    std::vector<int> previousReadCountsWithTVXAndZVertexAndSel8Full;
    std::vector<int> previousReadCountsWithTVXAndZVertexAndSel8FullPbPb;
    std::vector<int> previousReadCountsWithTVXAndZVertexAndSelMC;
    std::vector<int> previousReadCountsWithTVXAndZVertexAndSelMCFull;
    std::vector<int> previousReadCountsWithTVXAndZVertexAndSelMCFullPbPb;
    std::vector<int> previousReadCountsWithTVXAndZVertexAndSelUnanchoredMC;
    std::vector<int> previousReadCountsWithTVXAndZVertexAndSelTVX;
    std::vector<int> previousReadCountsWithTVXAndZVertexAndSel7;
    std::vector<int> previousReadCountsWithTVXAndZVertexAndSel7KINT7;
    std::vector<int> previousReadCountsWithCustom;

    int iPreviousDataFrame = 0;
    for (const auto& collisionCount : collisionCounts) {
      auto readCollisionCounterSpan = collisionCount.readCounts();
      auto readCollisionWithTVXCounterSpan = collisionCount.readCountsWithTVX();
      auto readCollisionWithTVXAndZVertexAndSel8CounterSpan = collisionCount.readCountsWithTVXAndZVertexAndSel8();
      auto readCollisionWithTVXAndZVertexAndSel8FullCounterSpan = collisionCount.readCountsWithTVXAndZVertexAndSel8Full();
      auto readCollisionWithTVXAndZVertexAndSel8FullPbPbCounterSpan = collisionCount.readCountsWithTVXAndZVertexAndSel8FullPbPb();
      auto readCollisionWithTVXAndZVertexAndSelMCCounterSpan = collisionCount.readCountsWithTVXAndZVertexAndSelMC();
      auto readCollisionWithTVXAndZVertexAndSelMCFullCounterSpan = collisionCount.readCountsWithTVXAndZVertexAndSelMCFull();
      auto readCollisionWithTVXAndZVertexAndSelMCFullPbPbCounterSpan = collisionCount.readCountsWithTVXAndZVertexAndSelMCFullPbPb();
      auto readCollisionWithTVXAndZVertexAndSelUnanchoredMCCounterSpan = collisionCount.readCountsWithTVXAndZVertexAndSelUnanchoredMC();
      auto readCollisionWithTVXAndZVertexAndSelTVXCounterSpan = collisionCount.readCountsWithTVXAndZVertexAndSelTVX();
      auto readCollisionWithTVXAndZVertexAndSel7CounterSpan = collisionCount.readCountsWithTVXAndZVertexAndSel7();
      auto readCollisionWithTVXAndZVertexAndSel7KINT7CounterSpan = collisionCount.readCountsWithTVXAndZVertexAndSel7KINT7();
      auto readCollisionWithCustomCounterSpan = collisionCount.readCountsWithCustom();

      if (iPreviousDataFrame == 0) {
        std::copy(readCollisionCounterSpan.begin(), readCollisionCounterSpan.end(), std::back_inserter(previousReadCounts));
        std::copy(readCollisionWithTVXCounterSpan.begin(), readCollisionWithTVXCounterSpan.end(), std::back_inserter(previousReadCountsWithTVX));
        std::copy(readCollisionWithTVXAndZVertexAndSel8CounterSpan.begin(), readCollisionWithTVXAndZVertexAndSel8CounterSpan.end(), std::back_inserter(previousReadCountsWithTVXAndZVertexAndSel8));
        std::copy(readCollisionWithTVXAndZVertexAndSel8FullCounterSpan.begin(), readCollisionWithTVXAndZVertexAndSel8FullCounterSpan.end(), std::back_inserter(previousReadCountsWithTVXAndZVertexAndSel8Full));
        std::copy(readCollisionWithTVXAndZVertexAndSel8FullPbPbCounterSpan.begin(), readCollisionWithTVXAndZVertexAndSel8FullPbPbCounterSpan.end(), std::back_inserter(previousReadCountsWithTVXAndZVertexAndSel8FullPbPb));
        std::copy(readCollisionWithTVXAndZVertexAndSelMCCounterSpan.begin(), readCollisionWithTVXAndZVertexAndSelMCCounterSpan.end(), std::back_inserter(previousReadCountsWithTVXAndZVertexAndSelMC));
        std::copy(readCollisionWithTVXAndZVertexAndSelMCFullCounterSpan.begin(), readCollisionWithTVXAndZVertexAndSelMCFullCounterSpan.end(), std::back_inserter(previousReadCountsWithTVXAndZVertexAndSelMCFull));
        std::copy(readCollisionWithTVXAndZVertexAndSelMCFullPbPbCounterSpan.begin(), readCollisionWithTVXAndZVertexAndSelMCFullPbPbCounterSpan.end(), std::back_inserter(previousReadCountsWithTVXAndZVertexAndSelMCFullPbPb));
        std::copy(readCollisionWithTVXAndZVertexAndSelUnanchoredMCCounterSpan.begin(), readCollisionWithTVXAndZVertexAndSelUnanchoredMCCounterSpan.end(), std::back_inserter(previousReadCountsWithTVXAndZVertexAndSelUnanchoredMC));
        std::copy(readCollisionWithTVXAndZVertexAndSelTVXCounterSpan.begin(), readCollisionWithTVXAndZVertexAndSelTVXCounterSpan.end(), std::back_inserter(previousReadCountsWithTVXAndZVertexAndSelTVX));
        std::copy(readCollisionWithTVXAndZVertexAndSel7CounterSpan.begin(), readCollisionWithTVXAndZVertexAndSel7CounterSpan.end(), std::back_inserter(previousReadCountsWithTVXAndZVertexAndSel7));
        std::copy(readCollisionWithTVXAndZVertexAndSel7KINT7CounterSpan.begin(), readCollisionWithTVXAndZVertexAndSel7KINT7CounterSpan.end(), std::back_inserter(previousReadCountsWithTVXAndZVertexAndSel7KINT7));
        std::copy(readCollisionWithCustomCounterSpan.begin(), readCollisionWithCustomCounterSpan.end(), std::back_inserter(previousReadCountsWithCustom));

      } else {
        for (unsigned int i = 0; i < previousReadCounts.size(); i++) { // in principle we only care about the first element, but might be interesting information to keep
          previousReadCounts[i] += readCollisionCounterSpan[i];
          previousReadCountsWithTVX[i] += readCollisionWithTVXCounterSpan[i];
          previousReadCountsWithTVXAndZVertexAndSel8[i] += readCollisionWithTVXAndZVertexAndSel8CounterSpan[i];
          previousReadCountsWithTVXAndZVertexAndSel8Full[i] += readCollisionWithTVXAndZVertexAndSel8FullCounterSpan[i];
          previousReadCountsWithTVXAndZVertexAndSel8FullPbPb[i] += readCollisionWithTVXAndZVertexAndSel8FullPbPbCounterSpan[i];
          previousReadCountsWithTVXAndZVertexAndSelMC[i] += readCollisionWithTVXAndZVertexAndSelMCCounterSpan[i];
          previousReadCountsWithTVXAndZVertexAndSelMCFull[i] += readCollisionWithTVXAndZVertexAndSelMCFullCounterSpan[i];
          previousReadCountsWithTVXAndZVertexAndSelMCFullPbPb[i] += readCollisionWithTVXAndZVertexAndSelMCFullPbPbCounterSpan[i];
          previousReadCountsWithTVXAndZVertexAndSelUnanchoredMC[i] += readCollisionWithTVXAndZVertexAndSelUnanchoredMCCounterSpan[i];
          previousReadCountsWithTVXAndZVertexAndSelTVX[i] += readCollisionWithTVXAndZVertexAndSelTVXCounterSpan[i];
          previousReadCountsWithTVXAndZVertexAndSel7[i] += readCollisionWithTVXAndZVertexAndSel7CounterSpan[i];
          previousReadCountsWithTVXAndZVertexAndSel7KINT7[i] += readCollisionWithTVXAndZVertexAndSel7KINT7CounterSpan[i];
          previousReadCountsWithCustom[i] += readCollisionWithCustomCounterSpan[i];
        }
      }
      iPreviousDataFrame++;
    }
    previousReadCounts.push_back(readCollisionCounter);
    previousReadCountsWithTVX.push_back(readCollisionWithTVXCounter);
    previousReadCountsWithTVXAndZVertexAndSel8.push_back(readCollisionWithTVXAndZVertexAndSel8Counter);
    previousReadCountsWithTVXAndZVertexAndSel8Full.push_back(readCollisionWithTVXAndZVertexAndSel8FullCounter);
    previousReadCountsWithTVXAndZVertexAndSel8FullPbPb.push_back(readCollisionWithTVXAndZVertexAndSel8FullPbPbCounter);
    previousReadCountsWithTVXAndZVertexAndSelMC.push_back(readCollisionWithTVXAndZVertexAndSelMCCounter);
    previousReadCountsWithTVXAndZVertexAndSelMCFull.push_back(readCollisionWithTVXAndZVertexAndSelMCFullCounter);
    previousReadCountsWithTVXAndZVertexAndSelMCFullPbPb.push_back(readCollisionWithTVXAndZVertexAndSelMCFullPbPbCounter);
    previousReadCountsWithTVXAndZVertexAndSelUnanchoredMC.push_back(readCollisionWithTVXAndZVertexAndSelUnanchoredMCCounter);
    previousReadCountsWithTVXAndZVertexAndSelTVX.push_back(readCollisionWithTVXAndZVertexAndSelTVXCounter);
    previousReadCountsWithTVXAndZVertexAndSel7.push_back(readCollisionWithTVXAndZVertexAndSel7Counter);
    previousReadCountsWithTVXAndZVertexAndSel7KINT7.push_back(readCollisionWithTVXAndZVertexAndSel7KINT7Counter);
    previousReadCountsWithCustom.push_back(readCollisionWithCustomCounter);

    storedCollisionCountsTable(previousReadCounts, previousReadCountsWithTVX, previousReadCountsWithTVXAndZVertexAndSel8, previousReadCountsWithTVXAndZVertexAndSel8Full, previousReadCountsWithTVXAndZVertexAndSel8FullPbPb, previousReadCountsWithTVXAndZVertexAndSelMC, previousReadCountsWithTVXAndZVertexAndSelMCFull, previousReadCountsWithTVXAndZVertexAndSelMCFullPbPb, previousReadCountsWithTVXAndZVertexAndSelUnanchoredMC, previousReadCountsWithTVXAndZVertexAndSelTVX, previousReadCountsWithTVXAndZVertexAndSel7, previousReadCountsWithTVXAndZVertexAndSel7KINT7, previousReadCountsWithCustom);
  }
  PROCESS_SWITCH(LuminosityProducer, processStoreCollisionCounting, "write out collision counting output table", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<LuminosityProducer>(cfgc, TaskName{"jet-luminosity-producer"}));

  return WorkflowSpec{tasks};
}
