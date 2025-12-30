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

/// \file collisionCutsGroup.h
/// \brief Configurable group for collision event selection cuts
///
/// This class follows the v0SelectionGroup pattern to enable automatic
/// configurable registration in O2 analysis tasks via ROOT serialization.
///
/// Usage in analysis tasks:
///   Configurable<collisionCutsGroup> collCuts{"collCuts", {}, "Collision event selection"};
///   o2::analysis::CollisonCuts colCutsChecker;
///
///   void init(InitContext const&) {
///     colCutsChecker.init(&registry);
///     colCutsChecker.setCuts(collCuts, /*isRun3*/ true);
///   }
///
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>

#ifndef PWGLF_UTILS_COLLISIONCUTSGROUP_H_
#define PWGLF_UTILS_COLLISIONCUTSGROUP_H_

#include "Common/CCDB/EventSelectionParams.h"

#include <Rtypes.h>

#include <iosfwd>

// Forward declaration to avoid circular dependency
namespace o2::analysis
{
class CollisonCuts;
}

/// \class collisionCutsGroup
/// \brief ROOT-serializable container for collision event selection parameters
///
/// This class encapsulates all collision selection parameters as member variables,
/// enabling automatic configurable registration via the DPL framework's ROOT
/// dictionary introspection. This eliminates the need for manual configurable
/// declarations in each analysis task.
///
/// Event selection flags are synchronized with o2::aod::evsel::EventSelectionFlags
/// via X-macros defined in EventSelectionFlagsMapping.def. To add new flags, simply
/// add an entry to that file - all getters, members, setters, and metadata are
/// auto-generated.
class CollisionCutsGroup
{
 public:
  /// Default constructor with standard default values
  CollisionCutsGroup();

  /// Apply all settings to a CollisonCuts object
  /// \param cuts The CollisonCuts object to configure
  void applyToCollisonCuts(o2::analysis::CollisonCuts& cuts) const;

  /// Print all selection settings for debugging
  void printSelections() const;

  // ===== Getters for non-flag members =====
  float getZvtxMax() const { return zvtxMax; }
  int getOccupancyInTimeRangeMax() const { return occupancyInTimeRangeMax; }
  int getOccupancyInTimeRangeMin() const { return occupancyInTimeRangeMin; }
  bool getCheckTrigger() const { return checkTrigger; }
  bool getCheckOffline() const { return checkOffline; }

  // ===== AUTO-GENERATED GETTERS for event selection flags =====
  // Generated from EventSelectionFlagsMapping.def
#define EVSEL_FLAG(enumVal, member, defaultVal, evtSelEnum, setter, getter, label, desc) \
  bool get##getter() const { return member; }
#include "EventSelectionFlagsMapping.def"
#undef EVSEL_FLAG

  // ===== Getters for Run2-specific flags (not in EventSelectionFlags) =====
  bool getApplyRun2AliEventCuts() const { return applyRun2AliEventCuts; }
  bool getApplyRun2INELgtZERO() const { return applyRun2INELgtZERO; }

 private:
  // ===== Basic collision parameters =====
  float zvtxMax;               ///< Maximum z-vertex position (cm)
  int occupancyInTimeRangeMax; ///< Maximum track occupancy (-1 = no cut)
  int occupancyInTimeRangeMin; ///< Minimum track occupancy (-1 = no cut)

  // ===== General event selection flags =====
  bool checkTrigger; ///< Check for trigger
  bool checkOffline; ///< Check for offline selection

  // ===== AUTO-GENERATED event selection flag members =====
  // Generated from EventSelectionFlagsMapping.def
#define EVSEL_FLAG(enumVal, member, defaultVal, evtSelEnum, setter, getter, label, desc) \
  bool member;                            ///< desc
#include "EventSelectionFlagsMapping.def" // NOLINT(build/include)
#undef EVSEL_FLAG

  // ===== Run2-specific selections (not in EventSelectionFlags) =====
  bool applyRun2AliEventCuts; ///< Run2 AliEventCuts
  bool applyRun2INELgtZERO;   ///< Run2 INEL>0 selection

  ClassDefNV(CollisionCutsGroup, 2);
};

#endif // PWGLF_UTILS_COLLISIONCUTSGROUP_H_
