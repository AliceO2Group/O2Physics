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

#ifndef PWGCF_FEMTODREAM_CORE_FEMTODREAMCOLLISIONSELECTION_H_
#define PWGCF_FEMTODREAM_CORE_FEMTODREAMCOLLISIONSELECTION_H_

#include <string>
#include <iostream>
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/EventSelection.h"
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
  void setCuts(float zvtxMax, bool checkTrigger, int trig, bool checkOffline, bool addCheckOffline, bool checkRun3)
  {
    mCutsSet = true;
    mZvtxMax = zvtxMax;
    mCheckTrigger = checkTrigger;
    mTrigger = static_cast<triggerAliases>(trig);
    mCheckOffline = checkOffline;
    mAddCheckOffline = addCheckOffline;
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
    mHistogramRegistry->add("Event/Zvtx", "; vtx_{z} (cm); Entries", kTH1F, {{300, -12.5, 12.5}});
    mHistogramRegistry->add("Event/MultPercentile", "; Multiplicity Percentile; Entries", kTH1F, {{100, 0, 100}});
    mHistogramRegistry->add("Event/MultPercentileVSMultNTracksPV", "; Multiplicity Percentile; MultNTracks", kTH2F, {{100, 0, 100}, {200, 0, 200}});
    mHistogramRegistry->add("Event/MultNTracksPV", "; MultNTracksPV; Entries", kTH1F, {{200, 0, 200}});
    mHistogramRegistry->add("Event/MultNTracklets", "; MultNTrackslets; Entries", kTH1F, {{300, 0, 300}});
    mHistogramRegistry->add("Event/MultTPC", "; MultTPC; Entries", kTH1F, {{600, 0, 600}});
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
      // all checks additional to sel8 are pending to be included in sel8, check them explicitly now and remove them once they have been added to sel8
      // kIsGoodZvtxFT0vsPV can be a dangerous cut because the default event selection value is rather tight
      // Remeber to open the cut (~4cm) with custom event selection task on hyperloop
      if (mAddCheckOffline && (!col.selection_bit(aod::evsel::kNoTimeFrameBorder) || !col.selection_bit(o2::aod::evsel::kNoITSROFrameBorder) || !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup) || !col.selection_bit(o2::aod::evsel::kIsVertexITSTPC) || !col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
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
  bool isEmptyCollision(C const& /*col*/, T const& tracks, TC& trackCuts)
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
  bool isEmptyCollision(C const& col, V const& V0s, VC& V0Cuts, T const& /*Tracks*/)
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

  template <typename C, typename Casc, typename CascC, typename T>
  bool isCollisionWithoutTrkCasc(C const& col, Casc const& Cascades, CascC& CascadeCuts, T const& /*Tracks*/)
  {
    // check if there is no selected Cascade
    for (auto const& Cascade : Cascades) {
      auto postrack = Cascade.template posTrack_as<T>();
      auto negtrack = Cascade.template negTrack_as<T>();
      auto bachtrack = Cascade.template bachelor_as<T>();
      if (CascadeCuts.isSelectedMinimal(col, Cascade, postrack, negtrack, bachtrack)) {
        return false;
      }
    }
    return true;
  }

  /// Some basic QA of the event
  /// \tparam T type of the collision
  /// \param col Collision
  template <typename T>
  void fillQA(T const& col, float cent)
  {
    if (mHistogramRegistry) {
      mHistogramRegistry->fill(HIST("Event/Zvtx"), col.posZ());
      mHistogramRegistry->fill(HIST("Event/MultNTracksPV"), col.multNTracksPV());
      mHistogramRegistry->fill(HIST("Event/MultTPC"), col.multTPC());
      if (mCheckIsRun3) {
        mHistogramRegistry->fill(HIST("Event/MultPercentile"), cent);
        mHistogramRegistry->fill(HIST("Event/MultPercentileVSMultNTracksPV"), cent, col.multNTracksPV());
      } else {
        mHistogramRegistry->fill(HIST("Event/MultNTracklets"), col.multTracklets());
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
  float computeSphericity(T1 const& /*col*/, T2 const& /*tracks*/)
  {
    return 2.f;
  }

 private:
  HistogramRegistry* mHistogramRegistry = nullptr; ///< For QA output
  bool mCutsSet = false;                           ///< Protection against running without cuts
  bool mCheckTrigger = false;                      ///< Check for trigger
  bool mCheckOffline = false;                      ///< Check for offline criteria (might change)
  bool mAddCheckOffline = false;                   ///< Additional check for offline criteria (added to sel8 soon)
  bool mCheckIsRun3 = false;                       ///< Check if running on Pilot Beam
  triggerAliases mTrigger = kINT7;                 ///< Trigger to check for
  float mZvtxMax = 999.f;                          ///< Maximal deviation from nominal z-vertex (cm)
};
} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_CORE_FEMTODREAMCOLLISIONSELECTION_H_
