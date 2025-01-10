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

/// \file FemtoUniverseCollisionSelection.h
/// \brief FemtoUniverseCollisionSelection - event selection within the o2femtouniverse framework
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch
/// \author Pritam Chakraborty, WUT Warsaw, pritam.chakraborty@pw.edu.pl

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSECOLLISIONSELECTION_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSECOLLISIONSELECTION_H_

#include <string>
#include "Common/CCDB/TriggerAliases.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"

using namespace o2::framework;

namespace o2::analysis::femto_universe
{

/// \class FemtoUniverseCollisionSelection
/// \brief Small selection class to check whether a given collision fulfills the specified selections
class FemtoUniverseCollisionSelection
{
 public:
  /// Destructor
  virtual ~FemtoUniverseCollisionSelection() = default;

  /// Pass the selection criteria to the class
  /// \param zvtxMax Maximal value of the z-vertex
  /// \param checkTrigger Whether or not to check for the trigger alias
  /// \param trig Requested trigger alias
  /// \param checkOffline Whether or not to check for offline selection criteria
  /// \param checkRun3 To check for the Run3 data
  /// \param centmin Minimum value of centrality selection
  /// \param centmax Maximum value of centrality selection
  void setCuts(float zvtxMax, bool checkTrigger, int trig, bool checkOffline, bool checkRun3, float centmin, float centmax)
  // void setCuts(float zvtxMax, bool checkTrigger, int trig, bool checkOffline, bool checkRun3)
  {
    mCutsSet = true;
    mZvtxMax = zvtxMax;
    mCheckTrigger = checkTrigger;
    mTrigger = static_cast<triggerAliases>(trig);
    mCheckOffline = checkOffline;
    mCheckIsRun3 = checkRun3;
    mCentMin = centmin;
    mCentMax = centmax;
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
    mHistogramRegistry->add("Event/Sphericity", "; Sphericity; Entries", kTH1I, {{200, 0, 3}});
  }

  /// Print some debug information
  void printCuts()
  {
    LOG(info) << "Debug information for FemtoUniverseCollisionSelection";
    LOG(info) << "Max. z-vertex: " << mZvtxMax;
    LOG(info) << "Check trigger: " << mCheckTrigger;
    LOG(info) << "Trigger: " << mTrigger;
    LOG(info) << " Check offline: " << mCheckOffline;
    LOG(info) << " Minimum Centrality: " << mCentMin;
    LOG(info) << " Maximum Centrality: " << mCentMax;
  }

  /// Check whether the collisions fulfills the specified selections
  /// \tparam T type of the collision
  /// \param col Collision
  /// \return whether or not the collisions fulfills the specified selections
  template <typename T>
  bool isSelected(T const& col)
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

  /// Check whether the collisions fulfills the specified selections for Run3
  /// \tparam T type of the collision
  /// \param col Collision
  /// \return whether or not the collisions fulfills the specified selections
  template <typename T>
  bool isSelectedRun3(T const& col)
  {
    if (std::abs(col.posZ()) > mZvtxMax) {
      return false;
    }
    if (mCheckOffline && !col.sel8()) {
      return false;
    }
    if ((col.centFT0C() < mCentMin) || (col.centFT0C() > mCentMax)) {
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

  /// Compute the sphericity of an event
  /// \tparam T1 type of the collision
  /// \tparam T2 type of the tracks
  /// \param col Collision
  /// \param tracks All tracks
  /// \return value of the sphericity of the event
  template <typename T1, typename T2>
  float computeSphericity(T1 const& /*col*/, T2 const& tracks)
  {
    double kS00 = 0;
    double kS11 = 0;
    double kS10 = 0;
    double sumPt = 0;
    int partNumber = 0;
    double spher = 0;

    for (const auto& p : tracks) {
      double phi = p.phi();
      double pT = p.pt();
      double px = pT * std::cos(phi);
      double py = pT * std::sin(phi);

      kS00 = kS00 + px * px / pT;
      kS11 = kS11 + py * py / pT;
      kS10 = kS10 + px * py / pT;
      sumPt = sumPt + pT;
      partNumber++;
    }

    if (sumPt != 0) {
      kS00 = kS00 / sumPt;
      kS11 = kS11 / sumPt;
      kS10 = kS10 / sumPt;

      double lambda1 = (kS00 + kS11 + std::sqrt((kS00 + kS11) * (kS00 + kS11) - 4.0 * (kS00 * kS11 - kS10 * kS10))) / 2.0;
      double lambda2 = (kS00 + kS11 - std::sqrt((kS00 + kS11) * (kS00 + kS11) - 4.0 * (kS00 * kS11 - kS10 * kS10))) / 2.0;

      if ((lambda1 + lambda2) != 0 && partNumber > 2) {
        spher = 2 * lambda2 / (lambda1 + lambda2);
      } else {
        spher = 2;
      }
    } else {
      spher = 2;
    }

    if (mHistogramRegistry) {
      mHistogramRegistry->fill(HIST("Event/Sphericity"), spher);
    }

    return spher;
  }

 private:
  HistogramRegistry* mHistogramRegistry = nullptr; ///< For QA output
  bool mCutsSet = false;                           ///< Protection against running without cuts
  bool mCheckTrigger = false;                      ///< Check for trigger
  bool mCheckOffline = false;                      ///< Check for offline criteria (might change)
  bool mCheckIsRun3 = false;                       ///< Check if running on Pilot Beam
  triggerAliases mTrigger = kINT7;                 ///< Trigger to check for
  float mZvtxMax = 999.f;                          ///< Maximal deviation from nominal z-vertex (cm)
  float mCentMin = 0.0;                            ///< Minimum centrality value
  float mCentMax = 100.0;                          ///< Maximum centrality value
};
} // namespace o2::analysis::femto_universe

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSECOLLISIONSELECTION_H_
