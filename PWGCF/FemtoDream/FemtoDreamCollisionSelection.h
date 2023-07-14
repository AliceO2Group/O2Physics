// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
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

#ifndef PWGCF_FEMTODREAM_FEMTODREAMCOLLISIONSELECTION_H_
#define PWGCF_FEMTODREAM_FEMTODREAMCOLLISIONSELECTION_H_

#include <string>
#include <iostream>
#include "Common/CCDB/TriggerAliases.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"

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
    mHistogramRegistry->add("Event/MultV0M", "; vMultV0M; Entries", kTH1F, {{16384, 0, 32768}});
    mHistogramRegistry->add("Event/MultT0M", "; vMultT0M; Entries", kTH1F, {{4096, 0, 8192}});
    mHistogramRegistry->add("Event/MultNTracksPV", "; vMultNTracksPV; Entries", kTH1F, {{120, 0, 120}});
    mHistogramRegistry->add("Event/MultNTracklets", "; vMultNTrackslets; Entries", kTH1F, {{300, 0, 300}});
    mHistogramRegistry->add("Event/MultTPC", "; vMultTPC; Entries", kTH1I, {{600, 0, 600}});
  }

  /// Print some debug information
  void printCuts()
  {
    LOG(info) << "Debug information for FemtoDreamCollisionSelection";
    LOG(info) << "Max. z-vertex: " << mZvtxMax;
    LOG(info) << "Check trigger: " << mCheckTrigger;
    LOG(info) << "Trigger: " << mTrigger;
    LOG(info) << " Check offline: " << mCheckOffline;
  }

  /// Check whether the collisions fulfills the specified selections
  /// \tparam T type of the collision
  /// \param col Collision
  /// \return whether or not the collisions fulfills the specified selections
  template <typename C>
  bool isSelectedCollision(C const& col)
  {
    if (std::abs(col.posZ()) > mZvtxMax) {
      return false;
    }
    if (mCheckIsRun3) {
      if (mCheckOffline && !col.sel8()) {
        return false;
      }
    } else {
      if (mCheckTrigger && !col.alias_bit(mTrigger)) {
        return false;
      }
      if (mCheckOffline && !col.sel7()) {
        return false;
      }
    }
    return true;
  }

  template <typename C, typename T, typename TC>
  bool isEmptyCollision(C const& col, T const& tracks, TC& trackCuts)
  {
    // check if there is no selected track
    for (auto const& track : tracks) {
      if (trackCuts.isSelectedMinimal(track)) {
        return false;
      }
    }
    return true;
  }

  template <typename C, typename V, typename VC, typename T>
  bool isEmptyCollision(C const& col, V const& V0s, VC& V0Cuts, T const& Tracks)
  {
    // check if there is no selected V0
    for (auto const& V0 : V0s) {
      auto postrack = V0.template posTrack_as<T>();
      auto negtrack = V0.template negTrack_as<T>();
      if (V0Cuts.isSelectedMinimal(col, V0, postrack, negtrack)) {
        return false;
      }
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
      mHistogramRegistry->fill(HIST("Event/zvtxhist"), col.posZ());
      mHistogramRegistry->fill(HIST("Event/MultT0M"), col.multFT0M());
      mHistogramRegistry->fill(HIST("Event/MultNTracksPV"), col.multNTracksPV());
      mHistogramRegistry->fill(HIST("Event/MultNTracklets"), col.multTracklets());
      mHistogramRegistry->fill(HIST("Event/MultTPC"), col.multTPC());
      if (mCheckIsRun3) {
        mHistogramRegistry->fill(HIST("Event/MultV0M"), col.multFV0M());
      } else {
        mHistogramRegistry->fill(HIST("Event/MultV0M"), 0.5 * (col.multFV0M())); // in AliPhysics, the VOM was defined by (V0A + V0C)/2.
      }
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

#endif // PWGCF_FEMTODREAM_FEMTODREAMCOLLISIONSELECTION_H_
