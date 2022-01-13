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

/// \file FemtoDreamCollisionSelection.h
/// \brief FemtoDreamCollisionSelection - event selection within the o2femtodream framework
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#ifndef ANALYSIS_TASKS_PWGCF_FEMTODREAM_FEMTODREAMCOLLISIONSELECTION_H_
#define ANALYSIS_TASKS_PWGCF_FEMTODREAM_FEMTODREAMCOLLISIONSELECTION_H_

#include "Common/CCDB/TriggerAliases.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"

#include <string>
#include <iostream>

using namespace o2::framework;

namespace o2::analysis::femtoDream
{

/// \class FemtoDreamCollisionSelection
/// \brief Small selection class to check whether a given collision fulfills the specified selections
class FemtoDreamCollisionSelection
{
 public:
  /// Destructor
  virtual ~FemtoDreamCollisionSelection() = default;

  /// Pass the selection criteria to the class
  /// \param zvtxMax Maximal value of the z-vertex
  /// \param checkTrigger whether or not to check for the trigger alias
  /// \param trig Requested trigger alias
  /// \param checkOffline whether or not to check for offline selection criteria
  void setCuts(float zvtxMax, bool checkTrigger, int trig, bool checkOffline, bool checkRun3)
  {
    mCutsSet = true;
    mZvtxMax = zvtxMax;
    mCheckTrigger = checkTrigger;
    mTrigger = static_cast<triggerAliases>(trig);
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
    mHistogramRegistry->add("Event/zvtxhist", "; vtx_{z} (cm); Entries", kTH1F, {{300, -12.5, 12.5}});
    mHistogramRegistry->add("Event/MultV0M", "; vMultV0M; Entries", kTH1F, {{600, 0, 600}});
    mHistogramRegistry->add("Event/MultT0A", "; vMultT0A; Entries", kTH1F, {{600, 0, 600}});
  }

  /// Print some debug information
  void printCuts()
  {
    std::cout << "Debug information for FemtoDreamCollisionSelection \n Max. z-vertex: " << mZvtxMax << "\n Check trigger: " << mCheckTrigger << "\n Trigger: " << mTrigger << "\n Check offline: " << mCheckOffline << "\n";
  }

  /// Check whether the collisions fulfills the specified selections
  /// \tparam T type of the collision
  /// \param col Collision
  /// \return whether or not the collisions fulfills the specified selections
  template <typename T>
  bool isSelected(T const& col)
  {
    bool pass = false;
    if (!mCheckIsRun3) {
      if (std::abs(col.posZ()) > mZvtxMax) {
        return pass;
      }
      if (mCheckTrigger && col.alias()[mTrigger] != 1) {
        return pass;
      }
      if (mCheckOffline && col.sel7() != 1) {
        return pass;
      }
      pass = true;
    } else if (mCheckIsRun3) {
      if (!col.sel8()) {
        return pass;
      }
      pass = true;
    }
    return pass;
  }

  /// Some basic QA of the event
  /// \tparam T type of the collision
  /// \param col Collision
  template <typename T>
  void fillQA(T const& col)
  {
    if (mHistogramRegistry) {
      mHistogramRegistry->fill(HIST("Event/zvtxhist"), col.posZ());
      mHistogramRegistry->fill(HIST("Event/MultV0M"), col.multV0M());
      mHistogramRegistry->fill(HIST("Event/MultT0A"), col.multT0A());
    }
  }

  /// \todo to be implemented!
  /// Compute the sphericity of an event
  /// Important here is that the filter on tracks does not interfere here!
  /// In Run 2 we used here global tracks within |eta| < 0.8
  /// \tparam T1 type of the collision
  /// \tparam T2 type of the tracks
  /// \param col Collision
  /// \param tracks All tracks
  /// \return value of the sphericity of the event
  template <typename T1, typename T2>
  float computeSphericity(T1 const& col, T2 const& tracks)
  {
    return 2.f;
  }

 private:
  HistogramRegistry* mHistogramRegistry = nullptr; ///< For QA output
  bool mCutsSet = false;                           ///< Protection against running without cuts
  bool mCheckTrigger = false;                      ///< Check for trigger
  bool mCheckOffline = false;                      ///< Check for offline criteria (might change)
  bool mCheckIsRun3 = false;                       ///< Check if running on Pilot Beam
  triggerAliases mTrigger = kINT7;                 ///< Trigger to check for
  float mZvtxMax = 999.f;                          ///< Maximal deviation from nominal z-vertex (cm)
};
} // namespace o2::analysis::femtoDream

#endif /* ANALYSIS_TASKS_PWGCF_FEMTODREAM_FEMTODREAMCOLLISIONSELECTION_H_ */
