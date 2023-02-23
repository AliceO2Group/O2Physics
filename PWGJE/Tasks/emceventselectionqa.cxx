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

// Monitoring task for EMCAL event selection
//
// Author: Markus Fasel

#include <unordered_map>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using bcEvSels = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels>;
using collEventSels = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;

struct EmcEventSelectionQA {
  o2::framework::HistogramRegistry mHistManager{"EMCALEventSelectionQAHistograms"};

  // Require EMCAL cells (CALO type 1)
  Filter emccellfilter = aod::calo::caloType == 1;

  void init(o2::framework::InitContext const&)
  {
    using o2HistType = o2::framework::HistType;
    using o2Axis = o2::framework::AxisSpec;

    o2Axis matchingAxis{3, -0.5, 2.5, "matchingStatus", "Matching status"}, // 0, no vertex,1 vertex found , 2 multiple vertices found
      bcAxis{4001, -0.5, 4000.5, "bcid", "BC ID"};

    mHistManager.add("hCollisionMatching", "hCollisionMatching", o2HistType::kTH1F, {matchingAxis});
    mHistManager.add("hCollisionMatchingReadout", "hCollisionMatching", o2HistType::kTH1F, {matchingAxis});
    mHistManager.add("hBCCollisions", "Bunch crossings of found collisions", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCEmcalReadout", "Bunch crossings with EMCAL trigger from CTP", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCTVX", "Bunch crossings with FIT TVX trigger from CTP", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCEmcalCellContent", "Bunch crossings with non-0 EMCAL cell content", o2HistType::kTH1F, {bcAxis});
  }

  Preslice<collEventSels> perFoundBC = aod::evsel::foundBCId;

  void process(bcEvSels const& bcs, collEventSels const& collisions, soa::Filtered<aod::Calos> const& cells)
  {
    std::unordered_map<uint64_t, int> cellGlobalBCs;
    // Build map of number of cells for corrected BCs using global BCs
    // used later in the determination whether a BC has EMC cell content (for speed reason)
    for (const auto& cell : cells) {
      auto globalbcid = cell.bc_as<bcEvSels>().globalBC();
      auto found = cellGlobalBCs.find(globalbcid);
      if (found != cellGlobalBCs.end()) {
        found->second++;
      } else {
        cellGlobalBCs.insert(std::pair<uint64_t, int>(globalbcid, 1));
      }
    }

    for (const auto& bc : bcs) {
      bool isEMCALreadout = false;
      auto bcID = bc.globalBC() % 3564;

      if (bc.runNumber() > 300000) {
        // in case of run3 not all BCs contain EMCAL data, require trigger selection also for min. bias
        // in addition select also L0/L1 triggers as triggers with EMCAL in reaodut
        if (bc.alias()[kTVXinEMC] || bc.alias()[kEMC7] || bc.alias()[kEG1] || bc.alias()[kEG2] || bc.alias()[kEJ1] || bc.alias()[kEJ2]) {
          isEMCALreadout = true;
        }
      } else {
        // run1/2: rely on trigger cluster, runlist must contain only runs with EMCAL in readout
        // Select min. bias trigger and EMCAL L0/L1 triggers
        if (bc.alias()[kINT7] || bc.alias()[kEMC7] || bc.alias()[kEG1] || bc.alias()[kEG2] || bc.alias()[kEJ1] || bc.alias()[kEJ2]) {
          isEMCALreadout = true;
        }
      }

      // Monitoring BCs with EMCAL trigger / readout / FIT trigger
      if (isEMCALreadout) {
        mHistManager.fill(HIST("hBCEmcalReadout"), bcID);
      }

      if (bc.selection()[kIsTriggerTVX]) {
        mHistManager.fill(HIST("hBCTVX"), bcID);
      }

      // lookup number of cells for global BC of this BC
      // avoid iteration over cell table for speed reason
      auto found = cellGlobalBCs.find(bc.globalBC());
      if (found != cellGlobalBCs.end()) {
        // require at least 1 cell for global BC
        if (found->second > 0) {
          mHistManager.fill(HIST("hBCEmcalCellContent"), bcID);
        }
      }

      auto collisionsGrouped = collisions.sliceBy(perFoundBC, bc.globalIndex());
      int collisionStatus = -1;
      if (!collisionsGrouped.size()) {
        collisionStatus = 0;
      } else if (collisionsGrouped.size() == 1) {
        collisionStatus = 1;
      } else {
        collisionStatus = 2;
      }
      if (collisionStatus >= 0) {
        mHistManager.fill(HIST("hCollisionMatching"), collisionStatus);
        if (isEMCALreadout) {
          mHistManager.fill(HIST("hCollisionMatchingReadout"), collisionStatus);
        }
      }
      if (collisionStatus > 0) {
        mHistManager.fill(HIST("hBCCollisions"), bcID);
      }
    }
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{
    o2::framework::adaptAnalysisTask<EmcEventSelectionQA>(cfgc)};
}
