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

/// \file collisionCutsGroup.cxx
/// \brief Implementation of CollisionCutsGroup class
///
/// Event selection flags are auto-generated from EventSelectionFlagsMapping.def
/// using X-macros for member initialization, setter calls, and logging.
///
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>
/// o2-linter: disable=name/workflow-file (Library implementation file, not a workflow)

#include "PWGLF/Utils/collisionCutsGroup.h"

#include "PWGLF/Utils/collisionCuts.h"

#include "Framework/Logger.h"

CollisionCutsGroup::CollisionCutsGroup()
  : zvtxMax{10.0f},
    occupancyInTimeRangeMax{-1},
    occupancyInTimeRangeMin{-1},
    checkTrigger{false},
    checkOffline{true},
// AUTO-GENERATED member initialization from EventSelectionFlagsMapping.def
#define EVSEL_FLAG(enumVal, member, defaultVal, evtSelEnum, setter, getter, label, desc) \
  member{defaultVal},
#include "EventSelectionFlagsMapping.def"
#undef EVSEL_FLAG
    // Manual Run2 flag initialization
    applyRun2AliEventCuts{true},
    applyRun2INELgtZERO{false}
{
  // Constructor - member initialization done in initializer list
}

void CollisionCutsGroup::applyToCollisonCuts(o2::analysis::CollisonCuts& cuts) const
{
  // AUTO-GENERATED setter method calls from EventSelectionFlagsMapping.def
#define EVSEL_FLAG(enumVal, member, defaultVal, evtSelEnum, setter, getter, label, desc) \
  cuts.setter(member);
#include "EventSelectionFlagsMapping.def"  // NOLINT(build/include)
#undef EVSEL_FLAG

  // Manual setter calls for non-flag members
  cuts.setTrackOccupancyInTimeRange(occupancyInTimeRangeMax, occupancyInTimeRangeMin);
  cuts.setApplyRun2AliEventCuts(applyRun2AliEventCuts);
  cuts.setApplyRun2INELgtZERO(applyRun2INELgtZERO);
}

void CollisionCutsGroup::printSelections() const
{
  LOGF(info, "");
  LOGF(info, "+++ Collision Selection Settings ++++++++++++++++++++++++");
  LOGF(info, "");
  LOGF(info, "Basic Parameters:");
  LOGF(info, "  Max z-vertex (cm) ......................: %.2f", zvtxMax);
  LOGF(info, "  Occupancy max ..........................: %d", occupancyInTimeRangeMax);
  LOGF(info, "  Occupancy min ..........................: %d", occupancyInTimeRangeMin);
  LOGF(info, "");
  LOGF(info, "General Flags:");
  LOGF(info, "  Check trigger ..........................: %s", checkTrigger ? "true" : "false");
  LOGF(info, "  Check offline ..........................: %s", checkOffline ? "true" : "false");
  LOGF(info, "");
  LOGF(info, "Run 3 Selections (from EventSelectionFlags):");

  // AUTO-GENERATED logging from EventSelectionFlagsMapping.def
#define EVSEL_FLAG(enumVal, member, defaultVal, evtSelEnum, setter, getter, label, desc) \
  LOGF(info, "  %-40s: %s", desc, member ? "true" : "false");
#include "EventSelectionFlagsMapping.def"  // NOLINT(build/include)
#undef EVSEL_FLAG

  LOGF(info, "");
  LOGF(info, "Run 2 Selections:");
  LOGF(info, "  Apply Run2 AliEventCuts ................: %s", applyRun2AliEventCuts ? "true" : "false");
  LOGF(info, "  Apply Run2 INEL>0 ......................: %s", applyRun2INELgtZERO ? "true" : "false");
  LOGF(info, "");
  LOGF(info, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
  LOGF(info, "");
}
