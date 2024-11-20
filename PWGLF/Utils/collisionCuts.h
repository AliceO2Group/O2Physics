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

/// \file collisionCuts.h
/// \brief Traditional event selection cuts for O2 analysis
///
/// Simply copied from FemtoDreamCollisionSelection.h
/// original author: Laura Serksnyte, TU München
///
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>

#ifndef PWGLF_UTILS_COLLISIONCUTS_H_
#define PWGLF_UTILS_COLLISIONCUTS_H_

#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"
#include "Common/DataModel/EventSelection.h"

using namespace o2::framework;

namespace o2::analysis
{

class CollisonCuts
{
 public:
  virtual ~CollisonCuts() = default;

  enum EvtSel {
    kAllEvent = 0,
    kFlagZvertex,
    kFlagTrigerTVX,
    kFlagTimeFrameBorder,
    kFlagITSROFrameBorder,
    kFlagSel8,
    kFlagVertexITSTPC,
    kFlagBunchPileup,
    kFlagZvtxFT0vsPV,
    kFlagOccupancy,
    kAllpassed
  };

  inline int binLabel(int index)
  {
    return index + 1;
  }

  /// \brief Pass the selection criteria to the class
  /// \param zvtxMax Maximal value of the z-vertex
  /// \param checkTrigger whether or not to check for the trigger alias
  /// \param trig Requested trigger alias
  /// \param checkOffline whether or not to check for offline selection criteria
  void setCuts(float zvtxMax, bool checkTrigger, int trig, bool checkOffline, bool checkRun3, bool triggerTVXsel = false, int trackOccupancyInTimeRangeMax = -1)
  {
    mCutsSet = true;
    mZvtxMax = zvtxMax;
    mCheckTrigger = checkTrigger;
    mTrigger = trig;
    mCheckOffline = checkOffline;
    mTriggerTVXselection = triggerTVXsel;
    mCheckIsRun3 = checkRun3;
    mApplyTFBorderCut = false;
    mApplyITSTPCvertex = false;
    mApplyZvertexTimedifference = false;
    mApplyPileupRejection = false;
    mApplyNoITSROBorderCut = false;
    mtrackOccupancyInTimeRangeMax = trackOccupancyInTimeRangeMax;
  }

  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  void init(HistogramRegistry* registry)
  {
    if (!mCutsSet) {
      LOGF(error, "Event selection not set - quitting!");
    }
    mHistogramRegistry = registry;
    mHistogramRegistry->add("Event/posZ", "; vtx_{z} (cm); Entries", kTH1F, {{250, -12.5, 12.5}});       // z-vertex histogram after event selections
    mHistogramRegistry->add("Event/posZ_noCut", "; vtx_{z} (cm); Entries", kTH1F, {{250, -12.5, 12.5}}); // z-vertex histogram before all selections
    if (mCheckIsRun3) {
      mHistogramRegistry->add("Event/CentFV0A", "; vCentV0A; Entries", kTH1F, {{110, 0, 110}});
      mHistogramRegistry->add("Event/CentFT0M", "; vCentT0M; Entries", kTH1F, {{110, 0, 110}});
      mHistogramRegistry->add("Event/CentFT0C", "; vCentT0C; Entries", kTH1F, {{110, 0, 110}});
      mHistogramRegistry->add("Event/CentFT0A", "; vCentT0A; Entries", kTH1F, {{110, 0, 110}});
      mHistogramRegistry->add("Event/posZ_ITSOnly", "; vtx_{z} (cm); Entries", kTH1F, {{250, -12.5, 12.5}});
      mHistogramRegistry->add("Event/posZ_ITSTPC", "; vtx_{z} (cm); Entries", kTH1F, {{250, -12.5, 12.5}});
      mHistogramRegistry->add("Event/trackOccupancyInTimeRange_noCut", "; Occupancy; Entries", kTH1F, {{500, 0., 20000.}});
    } else {
      mHistogramRegistry->add("Event/CentRun2V0M", "; vCentV0M; Entries", kTH1F, {{110, 0, 110}});
    }
    mHistogramRegistry->add("CollCutCounts", "; ; Entries", kTH1F, {{11, 0., 11.}});
    mHistogramRegistry->get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(binLabel(EvtSel::kAllEvent), "all");
    mHistogramRegistry->get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(binLabel(EvtSel::kFlagZvertex), "Zvtx");
    mHistogramRegistry->get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(binLabel(EvtSel::kFlagTrigerTVX), "IsTriggerTVX");
    mHistogramRegistry->get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(binLabel(EvtSel::kFlagTimeFrameBorder), "NoTimeFrameBorder");
    mHistogramRegistry->get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(binLabel(EvtSel::kFlagITSROFrameBorder), "NoITSROFrameBorder");
    mHistogramRegistry->get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(binLabel(EvtSel::kFlagSel8), "sel8");
    mHistogramRegistry->get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(binLabel(EvtSel::kFlagVertexITSTPC), "IsVertexITSTPC");
    mHistogramRegistry->get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(binLabel(EvtSel::kFlagBunchPileup), "NoSameBunchPileup");
    mHistogramRegistry->get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(binLabel(EvtSel::kFlagZvtxFT0vsPV), "IsGoodZvtxFT0vsPV");
    mHistogramRegistry->get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(binLabel(EvtSel::kFlagOccupancy), "LowOccupancy");
    mHistogramRegistry->get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(binLabel(EvtSel::kAllpassed), "Allpassed");
  }

  /// Print some debug information
  void printCuts()
  {
    LOGF(info, "Debug information for Collison Cuts \n Max. z-vertex: %f \n Check trigger: %d \n Trigger: %d \n Check offline: %d \n Check Run3: %d \n Trigger TVX selection: %d \n Apply time frame border cut: %d \n Apply ITS-TPC vertex: %d \n Apply Z-vertex time difference: %d \n Apply Pileup rejection: %d \n Apply NoITSRO frame border cut: %d \n Track occupancy in time range max: %d",
         mZvtxMax, mCheckTrigger, mTrigger, mCheckOffline, mCheckIsRun3, mTriggerTVXselection, mApplyTFBorderCut, mApplyITSTPCvertex, mApplyZvertexTimedifference, mApplyPileupRejection, mApplyNoITSROBorderCut, mtrackOccupancyInTimeRangeMax);
  }

  /// Set MB selection
  void setTriggerTVX(bool triggerTVXsel) { mTriggerTVXselection = triggerTVXsel; }

  /// Scan the trigger alias of the event
  void setInitialTriggerScan(bool checkTrigger) { mInitialTriggerScan = checkTrigger; }

  /// Set the time frame border cut
  void setApplyTFBorderCut(bool applyTFBorderCut) { mApplyTFBorderCut = applyTFBorderCut; }

  /// Set the ITS-TPC matching cut
  void setApplyITSTPCvertex(bool applyITSTPCvertex) { mApplyITSTPCvertex = applyITSTPCvertex; }

  /// Set the Z-vertex time difference cut
  void setApplyZvertexTimedifference(bool applyZvertexTimedifference) { mApplyZvertexTimedifference = applyZvertexTimedifference; }

  /// Set the Pileup rejection
  void setApplyPileupRejection(bool applyPileupRejection) { mApplyPileupRejection = applyPileupRejection; }

  /// Set the NoITSRO frame border cut
  void setApplyNoITSROBorderCut(bool applyNoITSROBorderCut) { mApplyNoITSROBorderCut = applyNoITSROBorderCut; }

  /// Check whether the collisions fulfills the specified selections
  /// \tparam T type of the collision
  /// \param col Collision
  /// \return whether or not the collisions fulfills the specified selections
  template <typename T>
  bool isSelected(T const& col)
  {
    mHistogramRegistry->fill(HIST("Event/posZ_noCut"), col.posZ());
    if (mCheckIsRun3) {
      mHistogramRegistry->fill(HIST("Event/trackOccupancyInTimeRange_noCut"), col.trackOccupancyInTimeRange());
    }
    mHistogramRegistry->fill(HIST("CollCutCounts"), EvtSel::kAllEvent);
    if (std::abs(col.posZ()) > mZvtxMax) {
      LOGF(debug, "Vertex out of range");
      return false;
    }
    mHistogramRegistry->fill(HIST("CollCutCounts"), EvtSel::kFlagZvertex);
    if (mCheckIsRun3) { // Run3 case
      if (!col.selection_bit(aod::evsel::kIsTriggerTVX) && mTriggerTVXselection) {
        LOGF(debug, "Offline selection TVX failed (Run3)");
        return false;
      }
      mHistogramRegistry->fill(HIST("CollCutCounts"), EvtSel::kFlagTrigerTVX);
      if (!col.selection_bit(aod::evsel::kNoTimeFrameBorder) && mApplyTFBorderCut) {
        LOGF(debug, "Time frame border cut failed");
        return false;
      }
      mHistogramRegistry->fill(HIST("CollCutCounts"), EvtSel::kFlagTimeFrameBorder);
      if (!col.selection_bit(aod::evsel::kNoITSROFrameBorder) && mApplyNoITSROBorderCut) {
        LOGF(debug, "NoITSRO frame border cut failed");
        return false;
      }
      mHistogramRegistry->fill(HIST("CollCutCounts"), EvtSel::kFlagITSROFrameBorder);
      if (!col.sel8() && mCheckOffline) {
        LOGF(debug, "Offline selection failed (Run3)");
        return false;
      }
      mHistogramRegistry->fill(HIST("CollCutCounts"), EvtSel::kFlagSel8);
      if (!col.selection_bit(o2::aod::evsel::kIsVertexITSTPC) && mApplyITSTPCvertex) {
        LOGF(debug, "ITS-TPC matching cut failed");
        return false;
      }
      mHistogramRegistry->fill(HIST("CollCutCounts"), EvtSel::kFlagVertexITSTPC);
      if (!col.selection_bit(o2::aod::evsel::kNoSameBunchPileup) && mApplyPileupRejection) {
        LOGF(debug, "Pileup rejection failed");
        return false;
      }
      mHistogramRegistry->fill(HIST("CollCutCounts"), EvtSel::kFlagBunchPileup);
      if (!col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV) && mApplyZvertexTimedifference) {
        LOGF(debug, "Z-vertex time difference cut failed");
        return false;
      }
      mHistogramRegistry->fill(HIST("CollCutCounts"), EvtSel::kFlagZvtxFT0vsPV);
      if (mtrackOccupancyInTimeRangeMax > 0 && col.trackOccupancyInTimeRange() > mtrackOccupancyInTimeRangeMax) {
        LOGF(debug, "trackOccupancyInTimeRange selection failed");
        return false;
      }
      mHistogramRegistry->fill(HIST("CollCutCounts"), EvtSel::kFlagOccupancy);
    } else { // Run2 case
      if (mCheckOffline && !col.sel7()) {
        LOGF(debug, "Offline selection failed (sel7)");
        return false;
      }
      mHistogramRegistry->fill(HIST("CollCutCounts"), EvtSel::kFlagTrigerTVX);
    }
    if (mCheckTrigger && !col.alias_bit(mTrigger)) {
      LOGF(debug, "Trigger selection failed");
      if (mInitialTriggerScan) { // Print out the trigger bits
        LOGF(debug, "Trigger scan initialized");
        for (int i = 0; i < kNaliases; i++) {
          if (col.alias_bit(i)) {
            LOGF(debug, "Trigger %d fired", i);
          }
        }
        mInitialTriggerScan = false;
      }
      return false;
    }
    mHistogramRegistry->fill(HIST("CollCutCounts"), EvtSel::kAllpassed);
    return true;
  }

  /// Some basic QA of the event
  /// \tparam T type of the collision
  /// \param col Collision
  template <typename T>
  void fillQA(T const& col)
  {
    if (mHistogramRegistry) {
      mHistogramRegistry->fill(HIST("Event/posZ"), col.posZ());
      if (!col.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        mHistogramRegistry->fill(HIST("Event/posZ_ITSOnly"), col.posZ());
      } else {
        mHistogramRegistry->fill(HIST("Event/posZ_ITSTPC"), col.posZ());
      }
      mHistogramRegistry->fill(HIST("Event/CentFV0A"), col.centFV0A());
      mHistogramRegistry->fill(HIST("Event/CentFT0M"), col.centFT0M());
      mHistogramRegistry->fill(HIST("Event/CentFT0C"), col.centFT0C());
      mHistogramRegistry->fill(HIST("Event/CentFT0A"), col.centFT0A());
    }
  }

  /// Some basic QA of the event
  /// \tparam T type of the collision
  /// \param col Collision
  template <typename T>
  void fillQARun2(T const& col)
  {
    if (mHistogramRegistry) {
      mHistogramRegistry->fill(HIST("Event/posZ"), col.posZ());
      mHistogramRegistry->fill(HIST("Event/CentRun2V0M"), col.centRun2V0M());
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry = nullptr; ///< For QA output
  bool mCutsSet = false;                           ///< Protection against running without cuts
  bool mCheckTrigger = false;                      ///< Check for trigger
  bool mTriggerTVXselection = false;               ///< Check for trigger TVX selection
  bool mCheckOffline = false;                      ///< Check for offline criteria (might change)
  bool mCheckIsRun3 = false;                       ///< Check if running on Pilot Beam
  bool mInitialTriggerScan = false;                ///< Check trigger when the event is first selected
  bool mApplyTFBorderCut = false;                  ///< Apply time frame border cut
  bool mApplyITSTPCvertex = false;                 ///< Apply at least one ITS-TPC track for vertexing
  bool mApplyZvertexTimedifference = false;        ///< removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference.
  bool mApplyPileupRejection = false;              ///< Pileup rejection
  bool mApplyNoITSROBorderCut = false;             ///< Apply NoITSRO frame border cut
  int mTrigger = kINT7;                            ///< Trigger to check for
  float mZvtxMax = 999.f;                          ///< Maximal deviation from nominal z-vertex (cm)
  int mtrackOccupancyInTimeRangeMax = -1;          ///< Maximum trackOccupancyInTimeRange cut (-1 no cut)
};
} // namespace o2::analysis

#endif // PWGLF_UTILS_COLLISIONCUTS_H_
