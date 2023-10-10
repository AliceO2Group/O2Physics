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

using namespace o2::framework;

namespace o2::analysis
{

class CollisonCuts
{
 public:
  virtual ~CollisonCuts() = default;

  /// \brief Pass the selection criteria to the class
  /// \param zvtxMax Maximal value of the z-vertex
  /// \param checkTrigger whether or not to check for the trigger alias
  /// \param trig Requested trigger alias
  /// \param checkOffline whether or not to check for offline selection criteria
  void setCuts(float zvtxMax, bool checkTrigger, int trig, bool checkOffline, bool checkRun3)
  {
    mCutsSet = true;
    mZvtxMax = zvtxMax;
    mCheckTrigger = checkTrigger;
    mTrigger = trig;
    mCheckOffline = checkOffline;
    mCheckIsRun3 = checkRun3;
  }

  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  void init(HistogramRegistry* registry)
  {
    if (!mCutsSet) {
      LOGF(error, "Event selection not set - quitting!");
    }
    mHistogramRegistry = registry;
    mHistogramRegistry->add("Event/posZ", "; vtx_{z} (cm); Entries", kTH1F, {{250, -12.5, 12.5}});
    mHistogramRegistry->add("Event/CentFV0A", "; vCentV0A; Entries", kTH1F, {{110, 0, 110}});
    mHistogramRegistry->add("Event/CentFT0M", "; vCentT0M; Entries", kTH1F, {{110, 0, 110}});
    mHistogramRegistry->add("Event/CentFT0C", "; vCentT0C; Entries", kTH1F, {{110, 0, 110}});
    mHistogramRegistry->add("Event/CentFT0A", "; vCentT0A; Entries", kTH1F, {{110, 0, 110}});
    mHistogramRegistry->add("Event/MultFT0M", "; FT0M signal; Entries", kTH1F, {{100000, 0, 100000}});
    mHistogramRegistry->add("Event/MultFT0C", "; FT0C signal; Entries", kTH1F, {{100000, 0, 100000}});
    mHistogramRegistry->add("Event/MultFT0A", "; FT0A signal; Entries", kTH1F, {{100000, 0, 100000}});
  }

  /// Print some debug information
  void printCuts()
  {
    LOGF(info, "Debug information for Collison Cuts \n Max. z-vertex: %f \n Check trigger: %d \n Trigger: %d \n Check offline: %d \n", mZvtxMax, mCheckTrigger, mTrigger, mCheckOffline);
  }

  /// Scan the trigger alias of the event
  void setInitialTriggerScan(bool checkTrigger) { mInitialTriggerScan = checkTrigger; }

  /// Check whether the collisions fulfills the specified selections
  /// \tparam T type of the collision
  /// \param col Collision
  /// \return whether or not the collisions fulfills the specified selections
  template <typename T>
  bool isSelected(T const& col)
  {
    if (std::abs(col.posZ()) > mZvtxMax) {
      LOGF(debug, "Vertex out of range");
      return false;
    }
    if (mCheckIsRun3) { // Run3 case
      if (mCheckOffline && !col.sel8()) {
        LOGF(debug, "Offline selection failed (Run3)");
        return false;
      }
    } else { // Run2 case
      if (mCheckOffline && !col.sel7()) {
        LOGF(debug, "Offline selection failed (sel7)");
        return false;
      }
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
      mHistogramRegistry->fill(HIST("Event/CentFV0A"), col.centFV0A());
      mHistogramRegistry->fill(HIST("Event/CentFT0M"), col.centFT0M());
      mHistogramRegistry->fill(HIST("Event/CentFT0C"), col.centFT0C());
      mHistogramRegistry->fill(HIST("Event/CentFT0A"), col.centFT0A());
      mHistogramRegistry->fill(HIST("Event/MultFT0M"), col.multFT0M());
      mHistogramRegistry->fill(HIST("Event/MultFT0C"), col.multFT0C());
      mHistogramRegistry->fill(HIST("Event/MultFT0A"), col.multFT0A());
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry = nullptr; ///< For QA output
  bool mCutsSet = false;                           ///< Protection against running without cuts
  bool mCheckTrigger = false;                      ///< Check for trigger
  bool mCheckOffline = false;                      ///< Check for offline criteria (might change)
  bool mCheckIsRun3 = false;                       ///< Check if running on Pilot Beam
  bool mInitialTriggerScan = false;                ///< Check trigger when the event is first selected
  int mTrigger = kINT7;                            ///< Trigger to check for
  float mZvtxMax = 999.f;                          ///< Maximal deviation from nominal z-vertex (cm)
};
} // namespace o2::analysis

#endif // PWGLF_UTILS_COLLISIONCUTS_H_
