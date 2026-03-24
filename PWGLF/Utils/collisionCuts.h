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
/// \author Hirak Kumar Koley <hirak.koley@cern.ch>

#ifndef PWGLF_UTILS_COLLISIONCUTS_H_
#define PWGLF_UTILS_COLLISIONCUTS_H_

#include "PWGLF/Utils/collisionCutsGroup.h"

#include "Common/DataModel/EventSelection.h"

#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"

#include <map>
#include <string>
#include <vector>

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
    kFlagVertexITSTPC,
    kFlagZvtxFT0vsPV,
    kFlagVertexTOFmatched,
    kFlagVertexTRDmatched,
    kFlagBunchPileup,
    kNoCollInTimeRangeStandard,
    kNoCollInTimeRangeNarrow,
    kNoCollInTimeRangeStrict,
    kNoCollInRofStandard,
    kNoCollInRofStrict,
    kNoHighMultCollInPrevRof,
    kIsGoodITSLayersAll,
    kIsGoodITSLayer3,
    kIsGoodITSLayer0123,
    kFlagBBT0A,
    kFlagBBT0C,
    kFlagSel8,
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
  /// \param checkOffline whether or not to check for offline selection criteria
  void setCuts(float zvtxMax, bool checkTrigger, bool checkOffline, bool checkRun3, bool triggerTVXsel = false, int trackOccupancyInTimeRangeMax = -1, int trackOccupancyInTimeRangeMin = -1)
  {
    mCutsSet = true;
    mZvtxMax = zvtxMax;
    mCheckTrigger = checkTrigger;
    mCheckOffline = checkOffline;
    mTriggerTVXselection = triggerTVXsel;
    mCheckIsRun3 = checkRun3;
    mApplyTFBorderCut = false;
    mApplyITSTPCvertex = false;
    mApplyZvertexTimedifference = false;
    mApplyPileupRejection = false;
    mApplyNoITSROBorderCut = false;
    mApplyCollInTimeRangeNarrow = false;
    mApplyCollInTimeRangeStandard = false;
    mApplyGoodITSLayersAll = false;
    mtrackOccupancyInTimeRangeMax = trackOccupancyInTimeRangeMax;
    mtrackOccupancyInTimeRangeMin = trackOccupancyInTimeRangeMin;
    mApplyRun2AliEventCuts = true;
    mApplyRun2INELgtZERO = false;
  }

  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  void init(o2::framework::HistogramRegistry* registry)
  {
    if (!mCutsSet) {
      LOGF(error, "Event selection not set - quitting!");
    }
    for (int i = 0; i < kNaliases; i++) {
      bitList.push_back(1 << i); // BIT(i)
    }
    mHistogramRegistry = registry;
    mHistogramRegistry->add("Event/posZ", "; vtx_{z} (cm); Entries", o2::framework::kTH1F, {{250, -12.5, 12.5}});       // z-vertex histogram after event selections
    mHistogramRegistry->add("Event/posZ_noCut", "; vtx_{z} (cm); Entries", o2::framework::kTH1F, {{250, -12.5, 12.5}}); // z-vertex histogram before all selections
    if (mCheckIsRun3) {
      mHistogramRegistry->add("Event/CentFT0M", "; FT0M Percentile; Entries", o2::framework::kTH1F, {{110, 0, 110}});
      mHistogramRegistry->add("Event/CentFT0C", "; FT0C Percentile; Entries", o2::framework::kTH1F, {{110, 0, 110}});
      mHistogramRegistry->add("Event/CentFT0A", "; FT0A Percentile; Entries", o2::framework::kTH1F, {{110, 0, 110}});
      mHistogramRegistry->add("Event/posZ_ITSOnly", "; vtx_{z} (cm); Entries", o2::framework::kTH1F, {{250, -12.5, 12.5}});
      mHistogramRegistry->add("Event/posZ_ITSTPC", "; vtx_{z} (cm); Entries", o2::framework::kTH1F, {{250, -12.5, 12.5}});
      mHistogramRegistry->add("Event/trackOccupancyInTimeRange_noCut", "; Occupancy; Entries", o2::framework::kTH1F, {{500, 0., 20000.}});
    } else {
      mHistogramRegistry->add("Event/CentRun2V0M", "; vCentV0M; Entries", o2::framework::kTH1F, {{110, 0, 110}});
    }
    mHistogramRegistry->add("CollCutCounts", "; ; Entries", o2::framework::kTH1F, {{kAllpassed + 1, 0, kAllpassed + 1}});
    mHistogramRegistry->get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(binLabel(EvtSel::kAllEvent), "all");
    mHistogramRegistry->get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(binLabel(EvtSel::kFlagZvertex), "Zvtx");

    // Automatically set labels from metadata registry
    for (const auto& [evtSelFlag, metadata] : sSelectionRegistry) {
      mHistogramRegistry->get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(binLabel(metadata.histogramBin), metadata.label);
    }

    mHistogramRegistry->get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(binLabel(EvtSel::kFlagSel8), "sel8");
    mHistogramRegistry->get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(binLabel(EvtSel::kFlagOccupancy), "LowOccupancy");
    mHistogramRegistry->get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(binLabel(EvtSel::kFlagVertexTOFmatched), "VertexTOFmatched");
    mHistogramRegistry->get<TH1>(HIST("CollCutCounts"))->GetXaxis()->SetBinLabel(binLabel(EvtSel::kAllpassed), "Allpassed");
  }

  /// Print some debug information
  void printCuts()
  {
    LOGF(info, "Debug information for Collison Cuts");
    LOGF(info, "Max. z-vertex: %f", mZvtxMax);
    LOGF(info, "Check trigger: %d", mCheckTrigger);
    LOGF(info, "Check offline: %d", mCheckOffline);
    LOGF(info, "Check Run3: %d", mCheckIsRun3);
    if (mCheckIsRun3) {
      LOGF(info, "Trigger TVX selection: %d", mTriggerTVXselection);
      LOGF(info, "Apply time frame border cut: %d", mApplyTFBorderCut);
      LOGF(info, "Apply ITS-TPC vertex: %d", mApplyITSTPCvertex);
      LOGF(info, "Apply NoCollInTimeRangeNarrow: %d", mApplyCollInTimeRangeNarrow);
      LOGF(info, "Apply Z-vertex time difference: %d", mApplyZvertexTimedifference);
      LOGF(info, "Apply Pileup rejection: %d", mApplyPileupRejection);
      LOGF(info, "Apply NoITSRO frame border cut: %d", mApplyNoITSROBorderCut);
      LOGF(info, "Track occupancy in time range max: %d", mtrackOccupancyInTimeRangeMax);
      LOGF(info, "Track occupancy in time range min: %d", mtrackOccupancyInTimeRangeMin);
      LOGF(info, "Apply NoCollInTimeRangeStandard: %d", mApplyCollInTimeRangeStandard);
    } else {
      LOGF(info, "Apply Run2 AliEventCuts: %d", mApplyRun2AliEventCuts);
      LOGF(info, "Apply Run2 INELgtZERO: %d", mApplyRun2INELgtZERO);
    }
  }

  /// Set MB selection
  void setTriggerTVX(bool triggerTVXsel)
  {
    mTriggerTVXselection = triggerTVXsel;
    setSelection(EvtSel::kFlagTrigerTVX, triggerTVXsel);
  }

  /// Set the time frame border cut
  void setApplyTFBorderCut(bool applyTFBorderCut)
  {
    mApplyTFBorderCut = applyTFBorderCut;
    setSelection(EvtSel::kFlagTimeFrameBorder, applyTFBorderCut);
  }

  /// Set the ITS-TPC matching cut
  void setApplyITSTPCvertex(bool applyITSTPCvertex)
  {
    mApplyITSTPCvertex = applyITSTPCvertex;
    setSelection(EvtSel::kFlagVertexITSTPC, applyITSTPCvertex);
  }

  /// Set the NoCollInTimeRangeNarrow cut
  void setApplyCollInTimeRangeNarrow(bool applyCollInTimeRangeNarrow)
  {
    mApplyCollInTimeRangeNarrow = applyCollInTimeRangeNarrow;
    setSelection(EvtSel::kNoCollInTimeRangeNarrow, applyCollInTimeRangeNarrow);
  }

  /// Set the Z-vertex time difference cut
  void setApplyZvertexTimedifference(bool applyZvertexTimedifference)
  {
    mApplyZvertexTimedifference = applyZvertexTimedifference;
    setSelection(EvtSel::kFlagZvtxFT0vsPV, applyZvertexTimedifference);
  }

  /// Set the Pileup rejection
  void setApplyPileupRejection(bool applyPileupRejection)
  {
    mApplyPileupRejection = applyPileupRejection;
    setSelection(EvtSel::kFlagBunchPileup, applyPileupRejection);
  }

  /// Set the NoITSRO frame border cut
  void setApplyNoITSROBorderCut(bool applyNoITSROBorderCut)
  {
    mApplyNoITSROBorderCut = applyNoITSROBorderCut;
    setSelection(EvtSel::kFlagITSROFrameBorder, applyNoITSROBorderCut);
  }

  /// Set the track occupancy in time range cut
  void setTrackOccupancyInTimeRange(int trackOccupancyInTimeRangeMax, int trackOccupancyInTimeRangeMin)
  {
    mtrackOccupancyInTimeRangeMax = trackOccupancyInTimeRangeMax;
    mtrackOccupancyInTimeRangeMin = trackOccupancyInTimeRangeMin;
  }
  /// Set the NoCollInTimeRangeStandard cut
  void setApplyCollInTimeRangeStandard(bool applyCollInTimeRangeStandard)
  {
    mApplyCollInTimeRangeStandard = applyCollInTimeRangeStandard;
    setSelection(EvtSel::kNoCollInTimeRangeStandard, applyCollInTimeRangeStandard);
  }

  /// Set the GoodITSLayersAll cut
  void setApplyGoodITSLayersAll(bool applyGoodITSLayersAll)
  {
    mApplyGoodITSLayersAll = applyGoodITSLayersAll;
    setSelection(EvtSel::kIsGoodITSLayersAll, applyGoodITSLayersAll);
  }

  /// Set the GoodITSLayer3 cut
  void setApplyGoodITSLayer3(bool applyGoodITSLayer3)
  {
    setSelection(EvtSel::kIsGoodITSLayer3, applyGoodITSLayer3);
  }

  /// Set the GoodITSLayer0123 cut
  void setApplyGoodITSLayer0123(bool applyGoodITSLayer0123)
  {
    setSelection(EvtSel::kIsGoodITSLayer0123, applyGoodITSLayer0123);
  }

  /// Set the NoCollInTimeRangeStrict cut
  void setApplyCollInTimeRangeStrict(bool apply)
  {
    mApplyCollInTimeRangeStrict = apply;
    setSelection(EvtSel::kNoCollInTimeRangeStrict, apply);
  }

  /// Set the NoCollInRofStandard cut
  void setApplyCollInRofStandard(bool apply)
  {
    mApplyCollInRofStandard = apply;
    setSelection(EvtSel::kNoCollInRofStandard, apply);
  }

  /// Set the NoCollInRofStrict cut
  void setApplyCollInRofStrict(bool apply)
  {
    mApplyCollInRofStrict = apply;
    setSelection(EvtSel::kNoCollInRofStrict, apply);
  }

  /// Set the NoHighMultCollInPrevRof cut
  void setApplyHighMultCollInPrevRof(bool apply)
  {
    mApplyHighMultCollInPrevRof = apply;
    setSelection(EvtSel::kNoHighMultCollInPrevRof, apply);
  }

  /// Set the VertexTOFmatched cut
  void setApplyVertexTOFmatched(bool apply)
  {
    mApplyVertexTOFmatched = apply;
    setSelection(EvtSel::kFlagVertexTOFmatched, apply);
  }

  /// Set the VertexTRDmatched cut
  void setApplyVertexTRDmatched(bool apply)
  {
    mApplyVertexTRDmatched = apply;
    setSelection(EvtSel::kFlagVertexTRDmatched, apply);
  }

  /// Set the BBT0A cut
  void setApplyBBT0A(bool apply)
  {
    mApplyBBT0A = apply;
    setSelection(EvtSel::kFlagBBT0A, apply);
  }

  /// Set the BBT0C cut
  void setApplyBBT0C(bool apply)
  {
    mApplyBBT0C = apply;
    setSelection(EvtSel::kFlagBBT0C, apply);
  }

  /// Set the Run2 AliEventCuts cut
  void setApplyRun2AliEventCuts(bool applyRun2AliEventCuts) { mApplyRun2AliEventCuts = applyRun2AliEventCuts; }

  /// Set the Run2 INELgtZERO cut
  void setApplyRun2INELgtZERO(bool applyRun2INELgtZERO) { mApplyRun2INELgtZERO = applyRun2INELgtZERO; }

  /// Overload that accepts a CollisionCutsGroup and isRun3 flag
  /// \param settings The CollisionCutsGroup containing all selection parameters
  /// \param checkRun3 Whether this is Run3 data (true) or Run2 (false)
  void setCuts(const CollisionCutsGroup& settings, bool checkRun3)
  {
    mCutsSet = true;
    mZvtxMax = settings.getZvtxMax();
    mCheckTrigger = settings.getCheckTrigger();
    mCheckOffline = settings.getCheckOffline();
    mCheckIsRun3 = checkRun3;
    mApplyTFBorderCut = false;
    mApplyITSTPCvertex = false;
    mApplyZvertexTimedifference = false;
    mApplyPileupRejection = false;
    mApplyNoITSROBorderCut = false;
    mApplyCollInTimeRangeNarrow = false;
    mApplyCollInTimeRangeStandard = false;
    mApplyGoodITSLayersAll = false;
    mtrackOccupancyInTimeRangeMax = settings.getOccupancyInTimeRangeMax();
    mtrackOccupancyInTimeRangeMin = settings.getOccupancyInTimeRangeMin();
    mApplyRun2AliEventCuts = settings.getApplyRun2AliEventCuts();
    mApplyRun2INELgtZERO = settings.getApplyRun2INELgtZERO();
    mTriggerTVXselection = false;

    // Apply all other settings via delegation
    settings.applyToCollisonCuts(*this);
  }

  /// Apply settings from a CollisionCutsGroup object (alternative to setCuts)
  /// \param settings The CollisionCutsGroup containing all selection parameters
  void applySettings(const CollisionCutsGroup& settings)
  {
    settings.applyToCollisonCuts(*this);
  }

  /// Generic selection control (NEW API)
  /// \param evtSelFlag Event selection flag from EvtSel enum
  /// \param enable Whether to enable (true) or disable (false) this selection
  void setSelection(int evtSelFlag, bool enable)
  {
    mEnabledSelections[evtSelFlag] = enable;
  }

  /// Get current state of a selection
  /// \param evtSelFlag Event selection flag from EvtSel enum
  /// \return Whether the selection is enabled
  bool getSelection(int evtSelFlag) const
  {
    auto it = mEnabledSelections.find(evtSelFlag);
    return it != mEnabledSelections.end() ? it->second : false;
  }

  /// Enable multiple selections at once
  /// \param selections Vector of pairs (evtSelFlag, enabled)
  void setSelections(const std::vector<std::pair<int, bool>>& selections)
  {
    for (const auto& [flag, enabled] : selections) {
      setSelection(flag, enabled);
    }
  }

  /// Check whether the collisions fulfills the specified selections
  /// \tparam T type of the collision
  /// \param col Collision
  /// \return whether or not the collisions fulfills the specified selections
  template <typename T>
  bool isSelected(T const& col, const bool QA = true)
  {
    if (QA) {
      mHistogramRegistry->fill(HIST("Event/posZ_noCut"), col.posZ());
      if (mCheckIsRun3) {
        mHistogramRegistry->fill(HIST("Event/trackOccupancyInTimeRange_noCut"), col.trackOccupancyInTimeRange());
      }
      mHistogramRegistry->fill(HIST("CollCutCounts"), EvtSel::kAllEvent);
    }
    if (std::abs(col.posZ()) > mZvtxMax) {
      LOGF(debug, "Vertex out of range");
      return false;
    }
    if (mInitialColBitScan) {
      for (const auto& bit : bitList) {
        if (col.selection_bit(bit)) {
          LOGF(info, "Trigger %d fired", bit);
        }
      }
      mInitialColBitScan = false;
    }
    if (QA) {
      mHistogramRegistry->fill(HIST("CollCutCounts"), EvtSel::kFlagZvertex);
    }
    if (mCheckIsRun3) { // Run3 case
      // Metadata-driven selection loop for registry-based cuts
      for (const auto& [evtSelFlag, enabled] : mEnabledSelections) {
        if (!enabled) {
          continue; // Skip disabled selections
        }

        auto it = sSelectionRegistry.find(evtSelFlag);
        if (it == sSelectionRegistry.end()) {
          LOGF(warning, "Unknown selection flag %d in mEnabledSelections", evtSelFlag);
          continue;
        }

        const auto& metadata = it->second;

        // Check the selection_bit
        if (!col.selection_bit(metadata.evselFlag)) {
          LOGF(debug, "%s", metadata.debugMessage);
          return false;
        }

        // Fill QA histogram
        if (QA) {
          mHistogramRegistry->fill(HIST("CollCutCounts"), metadata.histogramBin);
        }
      }

      // Special case: sel8 (not using selection_bit)
      if (!col.sel8() && mCheckOffline) {
        LOGF(debug, "Offline selection failed (Run3)");
        return false;
      }
      if (QA) {
        mHistogramRegistry->fill(HIST("CollCutCounts"), EvtSel::kFlagSel8);
      }

      // Special case: Occupancy cuts (custom logic)
      if (mtrackOccupancyInTimeRangeMax > 0 && col.trackOccupancyInTimeRange() > mtrackOccupancyInTimeRangeMax) {
        LOGF(debug, "trackOccupancyInTimeRange selection failed");
        return false;
      }
      if (mtrackOccupancyInTimeRangeMin > 0 && col.trackOccupancyInTimeRange() < mtrackOccupancyInTimeRangeMin) {
        LOGF(debug, "trackOccupancyInTimeRange selection failed");
        return false;
      }
      if (QA) {
        mHistogramRegistry->fill(HIST("CollCutCounts"), EvtSel::kFlagOccupancy);
      }
    } else { // Run2 case
      if (mCheckOffline && !col.sel7()) {
        LOGF(debug, "Offline selection failed (sel7)");
        return false;
      }
      auto bc = col.template bc_as<BCsWithRun2Info>();
      if (!(bc.eventCuts() & BIT(aod::Run2EventCuts::kAliEventCutsAccepted)) && !mApplyRun2AliEventCuts) {
        LOGF(debug, "Offline selection failed (AliEventCuts)");
        return false;
      }
      mHistogramRegistry->fill(HIST("CollCutCounts"), EvtSel::kFlagSel8);
      if (!(bc.eventCuts() & BIT(aod::Run2EventCuts::kINELgtZERO)) && !mApplyRun2INELgtZERO) {
        LOGF(debug, "INELgtZERO selection failed");
        return false;
      }
      if (QA) {
        mHistogramRegistry->fill(HIST("CollCutCounts"), EvtSel::kAllpassed);
      }
    }
    if (QA) {
      mHistogramRegistry->fill(HIST("CollCutCounts"), EvtSel::kAllpassed);
    }
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
  /// Metadata structure for event selections
  struct SelectionMetadata {
    int evselFlag;            ///< o2::aod::evsel::EventSelectionFlags value
    const char* label;        ///< Histogram label
    const char* debugMessage; ///< Debug message for LOGF
    int histogramBin;         ///< Position in EvtSel enum (for histogram bin)
  };

  /// Central registry mapping EvtSel enum to metadata (auto-generated from EventSelectionFlagsMapping.def)
  inline static const std::map<int, SelectionMetadata> sSelectionRegistry = {
  // AUTO-GENERATED from EventSelectionFlagsMapping.def
#define EVSEL_FLAG(enumVal, member, defaultVal, evtSelEnum, setter, getter, label, desc) \
  {EvtSel::evtSelEnum, {aod::evsel::enumVal, label, desc " selection failed", EvtSel::evtSelEnum}},
#include "PWGLF/Utils/EventSelectionFlagsMapping.def"
#undef EVSEL_FLAG
  };

  /// Dynamic selection state: which selections are enabled
  std::map<int, bool> mEnabledSelections;

  using BCsWithRun2Info = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;
  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr; ///< For QA output
  std::vector<int> bitList;
  bool mCutsSet = false;                      ///< Protection against running without cuts
  bool mInitialColBitScan = true;             ///< Scan for collision bit
  bool mCheckTrigger = false;                 ///< Check for trigger
  bool mTriggerTVXselection = false;          ///< Check for trigger TVX selection
  bool mCheckOffline = false;                 ///< Check for offline criteria (might change)
  bool mCheckIsRun3 = false;                  ///< Check if running on Pilot Beam
  bool mApplyTFBorderCut = false;             ///< Apply time frame border cut
  bool mApplyITSTPCvertex = false;            ///< Apply at least one ITS-TPC track for vertexing
  bool mApplyCollInTimeRangeNarrow = false;   ///< Apply NoCollInTimeRangeNarrow selection
  bool mApplyZvertexTimedifference = false;   ///< removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference.
  bool mApplyPileupRejection = false;         ///< Pileup rejection
  bool mApplyNoITSROBorderCut = false;        ///< Apply NoITSRO frame border cut
  bool mApplyCollInTimeRangeStandard = false; ///< Apply NoCollInTimeRangeStandard selection
  bool mApplyGoodITSLayersAll = false;        ///< Apply GoodITSLayersAll selection (cuts time intervals with dead ITS staves)
  bool mApplyCollInTimeRangeStrict = false;   ///< Apply NoCollInTimeRangeStrict selection
  bool mApplyCollInRofStandard = false;       ///< Apply NoCollInRofStandard selection
  bool mApplyCollInRofStrict = false;         ///< Apply NoCollInRofStrict selection
  bool mApplyHighMultCollInPrevRof = false;   ///< Apply NoHighMultCollInPrevRof selection
  bool mApplyVertexTOFmatched = false;        ///< Apply vertex TOF matched selection
  bool mApplyVertexTRDmatched = false;        ///< Apply vertex TRD matched selection
  bool mApplyBBT0A = false;                   ///< Apply T0A beam-beam timing selection
  bool mApplyBBT0C = false;                   ///< Apply T0C beam-beam timing selection
  bool mApplyRun2AliEventCuts = true;         ///< Apply Run2 AliEventCuts
  bool mApplyRun2INELgtZERO = false;          ///< Apply Run2 INELgtZERO selection
  float mZvtxMax = 999.f;                     ///< Maximal deviation from nominal z-vertex (cm)
  int mtrackOccupancyInTimeRangeMax = -1;     ///< Maximum trackOccupancyInTimeRange cut (-1 no cut)
  int mtrackOccupancyInTimeRangeMin = -1;     ///< Minimum trackOccupancyInTimeRange cut (-1 no cut)
};
} // namespace o2::analysis

#endif // PWGLF_UTILS_COLLISIONCUTS_H_
