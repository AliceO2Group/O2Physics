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

/// \file emcEventSelectionQA.cxx
/// \brief Monitoring task for EMCAL event selection
/// \author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laoratory

#include "PWGJE/Core/utilsBcSelEMC.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/EventSelection.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <cstdint>
#include <unordered_map>
#include <utility>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::emc_evsel;

using BCEvSels = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels>;
using CollEventSels = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
using FilteredCells = o2::soa::Filtered<aod::Calos>;

struct EmcEventSelectionQA {

  EMCEventSelection emcEvSel; // event selection and monitoring
  o2::framework::HistogramRegistry mHistManager{"EMCALEventSelectionQAHistograms"};

  // Require EMCAL cells (CALO type 1)
  Filter emccellfilter = aod::calo::caloType == 1;

  const int mRun3MinNumber = 300000;

  void init(o2::framework::InitContext const&)
  {
    using O2HistType = o2::framework::HistType;
    using O2Axis = o2::framework::AxisSpec;

    O2Axis matchingAxis{3, -0.5, 2.5, "Matching Status (0, 1, 2+ collisions)", "Matching status"}, // 0, no vertex,1 vertex found , 2 multiple vertices found
      bcAxis{4001, -0.5, 4000.5, "bcid", "BC ID"},
      amplitudeAxisLarge{1000, 0., 100., "amplitudeLarge", "Amplitude (GeV)"},
      timeAxisLarge{1500, -600, 900, "celltime", "#it{t}_{cell} (ns)"};

    mHistManager.add("hCollisionMatching", "Collision Status", O2HistType::kTH1F, {matchingAxis});
    mHistManager.add("hCollisionMatchingReadout", "Collision Status EMCAL Readout", O2HistType::kTH1F, {matchingAxis});
    mHistManager.add("hCollisionMatchingMB", "Collision Status EMCAL MB", O2HistType::kTH1F, {matchingAxis});
    mHistManager.add("hCollisionMatching0EMC", "Collision Status EMCAL L0 trigger", O2HistType::kTH1F, {matchingAxis});
    mHistManager.add("hCollisionMatching0DMC", "Collision Status DCAL L0 tr", O2HistType::kTH1F, {matchingAxis});
    mHistManager.add("hCollisionMatchingEG1", "Collision Status EG1 trigger", O2HistType::kTH1F, {matchingAxis});
    mHistManager.add("hCollisionMatchingDG1", "Collision Status DG1 trigger", O2HistType::kTH1F, {matchingAxis});
    mHistManager.add("hCollisionMatchingEG2", "Collision Status EG2 trigger", O2HistType::kTH1F, {matchingAxis});
    mHistManager.add("hCollisionMatchingDG2", "Collision Status DG2 trigger", O2HistType::kTH1F, {matchingAxis});
    mHistManager.add("hCollisionMatchingEJ1", "Collision Status EJ1 trigger", O2HistType::kTH1F, {matchingAxis});
    mHistManager.add("hCollisionMatchingDJ1", "Collision Status DJ1 trigger", O2HistType::kTH1F, {matchingAxis});
    mHistManager.add("hCollisionMatchingEJ2", "Collision Status EJ2 trigger", O2HistType::kTH1F, {matchingAxis});
    mHistManager.add("hCollisionMatchingDJ2", "Collision Status DJ2 trigger", O2HistType::kTH1F, {matchingAxis});
    mHistManager.add("hBCCollisions", "Bunch crossings of found collisions", O2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCEmcalReadout", "Bunch crossings with EMCAL trigger from CTP", O2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCEmcalMB", "Bunch crossings with EMCAL MB from CTP", O2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCEmcal0EMC", "Bunch crossings with EMCAL L0 from CTP", O2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCEmcal0DMC", "Bunch crossings with DCAL L0 from CTP", O2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCEmcalEG1", "Bunch crossings with EG1 trigger from CTP", O2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCEmcalDG1", "Bunch crossings with DG1 trigger from CTP", O2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCEmcalEG2", "Bunch crossings with EG1 trigger from CTP", O2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCEmcalDG2", "Bunch crossings with DG2 trigger from CTP", O2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCEmcalEJ1", "Bunch crossings with EJ1 trigger from CTP", O2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCEmcalDJ1", "Bunch crossings with DJ1 trigger from CTP", O2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCEmcalEJ2", "Bunch crossings with EJ2 trigger from CTP", O2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCEmcalDJ2", "Bunch crossings with DJ2 trigger from CTP", O2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCTVX", "Bunch crossings with FIT TVX trigger from CTP", O2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCNoTVXEmcalReadout", "Bunch crossings with no FIT TVX trigger from CTP but with EMCal im readout", O2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCEmcalCellContent", "Bunch crossings with non-0 EMCAL cell content", O2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCCollisionCounter_TVX", "Number of BCs with a certain number of rec. colls", O2HistType::kTH2F, {bcAxis, matchingAxis});
    mHistManager.add("hBCEMCalReadoutAndEmcalCellContent", "Bunch crossings with EMCAL trigger from CTP and non-0 EMCAL cell content", O2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCNotEMCalReadoutButEmcalCellContent", "Bunch crossings without EMCAL trigger from CTP but with non-0 EMCAL cell content", O2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCNotAcceptedButEMCalReadout", "Bunch crossings with EMCAL trigger from CTP but not accpeted due to BC selection", O2HistType::kTH1F, {bcAxis});
    mHistManager.add("hBCNotAcceptedButEmcalCellContent", "Bunch crossings with non-0 EMCAL cell content but not accpeted due to BC selection", O2HistType::kTH1F, {bcAxis});
    mHistManager.add("hAmplitudevsCellTimeNoReadout", "Amplitude vs cell time for bunch crossings without EMCAL trigger from CTP but with non-0 EMCAL cell content", O2HistType::kTH2D, {timeAxisLarge, amplitudeAxisLarge});

    initCollisionHistogram(mHistManager.get<TH1>(HIST("hCollisionMatching")).get());
    initCollisionHistogram(mHistManager.get<TH1>(HIST("hCollisionMatchingReadout")).get());
    initCollisionHistogram(mHistManager.get<TH1>(HIST("hCollisionMatchingMB")).get());
    initCollisionHistogram(mHistManager.get<TH1>(HIST("hCollisionMatching0EMC")).get());
    initCollisionHistogram(mHistManager.get<TH1>(HIST("hCollisionMatching0DMC")).get());
    initCollisionHistogram(mHistManager.get<TH1>(HIST("hCollisionMatchingEG1")).get());
    initCollisionHistogram(mHistManager.get<TH1>(HIST("hCollisionMatchingDG1")).get());
    initCollisionHistogram(mHistManager.get<TH1>(HIST("hCollisionMatchingEG2")).get());
    initCollisionHistogram(mHistManager.get<TH1>(HIST("hCollisionMatchingDG2")).get());
    initCollisionHistogram(mHistManager.get<TH1>(HIST("hCollisionMatchingEJ1")).get());
    initCollisionHistogram(mHistManager.get<TH1>(HIST("hCollisionMatchingDJ1")).get());
    initCollisionHistogram(mHistManager.get<TH1>(HIST("hCollisionMatchingEJ2")).get());
    initCollisionHistogram(mHistManager.get<TH1>(HIST("hCollisionMatchingDJ2")).get());

    emcEvSel.addHistograms(mHistManager); // collision monitoring
  }

  PresliceUnsorted<CollEventSels> perFoundBC = aod::evsel::foundBCId;
  Preslice<FilteredCells> cellsPerFoundBC = aod::calo::bcId;

  void process(BCEvSels const& bcs, CollEventSels const& collisions, soa::Filtered<aod::Calos> const& cells)
  {
    std::unordered_map<uint64_t, int> cellGlobalBCs;
    // Build map of number of cells for corrected BCs using global BCs
    // used later in the determination whether a BC has EMC cell content (for speed reason)
    for (const auto& cell : cells) {
      auto globalbcid = cell.bc_as<BCEvSels>().globalBC();
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

      // get bitmask with bc selection info
      const auto rejectionMask = emcEvSel.getEMCCollisionRejectionMask<true, BCEvSels::iterator>(bc);
      // monitor the satisfied event selections
      emcEvSel.fillHistograms(rejectionMask);

      if (bc.runNumber() > mRun3MinNumber) {
        // in case of run3 not all BCs contain EMCAL data, require trigger selection also for min. bias
        // in addition select also L0/L1 triggers as triggers with EMCAL in reaodut
        if (bc.alias_bit(kTVXinEMC) || bc.alias_bit(kEMC7) || bc.alias_bit(kEG1) || bc.alias_bit(kEG2) || bc.alias_bit(kDG1) || bc.alias_bit(kDG2) || bc.alias_bit(kEJ1) || bc.alias_bit(kEJ2) || bc.alias_bit(kDJ1) || bc.alias_bit(kDJ2)) {
          isEMCALreadout = true;
        }
      } else {
        // run1/2: rely on trigger cluster, runlist must contain only runs with EMCAL in readout
        // Select min. bias trigger and EMCAL L0/L1 triggers
        if (bc.alias_bit(kINT7) || bc.alias_bit(kEMC7) || bc.alias_bit(kEG1) || bc.alias_bit(kEG2) || bc.alias_bit(kEJ1) || bc.alias_bit(kEJ2)) {
          isEMCALreadout = true;
        }
      }

      // lookup number of cells for global BC of this BC
      // avoid iteration over cell table for speed reason
      auto found = cellGlobalBCs.find(bc.globalBC());

      if (rejectionMask != 0) {
        // at least one event selection not satisfied --> reject the candidate
        continue;
      } else {
        if (isEMCALreadout) {
          mHistManager.fill(HIST("hBCNotAcceptedButEMCalReadout"), bcID);
        }
        if (found != cellGlobalBCs.end()) {
          // require at least 1 cell for global BC
          if (found->second > 0) {
            mHistManager.fill(HIST("hBCNotAcceptedButEmcalCellContent"), bcID);
          }
        }
      }

      // Monitoring BCs with EMCAL trigger / readout / FIT trigger
      if (isEMCALreadout) {
        mHistManager.fill(HIST("hBCEmcalReadout"), bcID);
        // various triggers
        if (bc.alias_bit(kTVXinEMC)) {
          mHistManager.fill(HIST("hBCEmcalMB"), bcID);
        }
        if (bc.alias_bit(kEMC7)) {
          mHistManager.fill(HIST("hBCEmcal0EMC"), bcID);
        }
        if (bc.alias_bit(kDMC7)) {
          mHistManager.fill(HIST("hBCEmcal0DMC"), bcID);
        }
        if (bc.alias_bit(kEG1)) {
          mHistManager.fill(HIST("hBCEmcalEG1"), bcID);
        }
        if (bc.alias_bit(kDG1)) {
          mHistManager.fill(HIST("hBCEmcalDG1"), bcID);
        }
        if (bc.alias_bit(kEG2)) {
          mHistManager.fill(HIST("hBCEmcalEG2"), bcID);
        }
        if (bc.alias_bit(kDG2)) {
          mHistManager.fill(HIST("hBCEmcalDG2"), bcID);
        }
        if (bc.alias_bit(kEJ1)) {
          mHistManager.fill(HIST("hBCEmcalEJ1"), bcID);
        }
        if (bc.alias_bit(kDJ1)) {
          mHistManager.fill(HIST("hBCEmcalDJ1"), bcID);
        }
        if (bc.alias_bit(kEJ2)) {
          mHistManager.fill(HIST("hBCEmcalEJ2"), bcID);
        }
        if (bc.alias_bit(kDJ2)) {
          mHistManager.fill(HIST("hBCEmcalDJ2"), bcID);
        }
      }

      // lookup number of cells for global BC of this BC
      // avoid iteration over cell table for speed reason
      if (found != cellGlobalBCs.end()) {
        // require at least 1 cell for global BC
        if (found->second > 0) {
          mHistManager.fill(HIST("hBCEmcalCellContent"), bcID);
          if (isEMCALreadout) {
            mHistManager.fill(HIST("hBCEMCalReadoutAndEmcalCellContent"), bcID);
          } else {
            mHistManager.fill(HIST("hBCNotEMCalReadoutButEmcalCellContent"), bcID);
            auto cellsInBC = cells.sliceBy(cellsPerFoundBC, bc.globalIndex());
            for (const auto& cell : cellsInBC) {
              mHistManager.fill(HIST("hAmplitudevsCellTimeNoReadout"), cell.time(), cell.amplitude());
            }
          }
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

      if (bc.selection_bit(aod::evsel::kIsTriggerTVX)) {
        mHistManager.fill(HIST("hBCTVX"), bcID);
        mHistManager.fill(HIST("hBCCollisionCounter_TVX"), bcID, collisionStatus);
      }

      if (collisionStatus >= 0) {
        mHistManager.fill(HIST("hCollisionMatching"), collisionStatus);
        if (isEMCALreadout) {
          mHistManager.fill(HIST("hCollisionMatchingReadout"), collisionStatus);
          if (!bc.selection_bit(aod::evsel::kIsTriggerTVX)) {
            mHistManager.fill(HIST("hBCNoTVXEmcalReadout"), bcID);
          }
          // various triggers
          if (bc.alias_bit(kTVXinEMC)) {
            mHistManager.fill(HIST("hCollisionMatchingMB"), collisionStatus);
          }
          if (bc.alias_bit(kEMC7)) {
            mHistManager.fill(HIST("hCollisionMatching0EMC"), collisionStatus);
          }
          if (bc.alias_bit(kDMC7)) {
            mHistManager.fill(HIST("hCollisionMatching0DMC"), collisionStatus);
          }
          if (bc.alias_bit(kEG1)) {
            mHistManager.fill(HIST("hCollisionMatchingEG1"), collisionStatus);
          }
          if (bc.alias_bit(kDG1)) {
            mHistManager.fill(HIST("hCollisionMatchingDG1"), collisionStatus);
          }
          if (bc.alias_bit(kEG2)) {
            mHistManager.fill(HIST("hCollisionMatchingEG2"), collisionStatus);
          }
          if (bc.alias_bit(kDG2)) {
            mHistManager.fill(HIST("hCollisionMatchingDG2"), collisionStatus);
          }
          if (bc.alias_bit(kEJ1)) {
            mHistManager.fill(HIST("hCollisionMatchingEJ1"), collisionStatus);
          }
          if (bc.alias_bit(kDJ1)) {
            mHistManager.fill(HIST("hCollisionMatchingDJ1"), collisionStatus);
          }
          if (bc.alias_bit(kEJ2)) {
            mHistManager.fill(HIST("hCollisionMatchingEJ2"), collisionStatus);
          }
          if (bc.alias_bit(kDJ2)) {
            mHistManager.fill(HIST("hCollisionMatchingDJ2"), collisionStatus);
          }
        }
      }
      if (collisionStatus > 0) {
        mHistManager.fill(HIST("hBCCollisions"), bcID);
      }
    }
  }

  void initCollisionHistogram(TH1* hist)
  {
    // Beautify collision type lables
    hist->GetXaxis()->SetTitle("Collision status");
    hist->GetXaxis()->SetBinLabel(1, "No vertex");
    hist->GetXaxis()->SetBinLabel(2, "1 vertex");
    hist->GetXaxis()->SetBinLabel(3, "Pileup");
    hist->GetYaxis()->SetTitle("Number of BCs");
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{
    o2::framework::adaptAnalysisTask<EmcEventSelectionQA>(cfgc)};
}
