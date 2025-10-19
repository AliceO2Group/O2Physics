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

/// \file utilsTpcSkimsTableCreator.h
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>
/// \author Christian Sonnabend <christian.sonnabend@cern.ch>
/// \author Jeremy Wilkinson <jeremy.wilkinson@cern.ch>
/// \author Ana Marin <ana.marin@cern.ch>
/// \author Oleksii Lubynets <oleksii.lubynets@cern.ch>
/// \brief  Helper functions used both in V0 and TOF structs of tpcSkimsTableCreator.cxx

#ifndef DPG_TASKS_TPC_UTILSTPCSKIMSTABLECREATOR_H_
#define DPG_TASKS_TPC_UTILSTPCSKIMSTABLECREATOR_H_

namespace o2::dpg_tpcskimstablecreator
{
enum {
  TrackSelectionNoCut = 0,
  TrackSelectionGlobalTrack,
  TrackSelectionTrackWoPtEta,
  TrackSelectionGlobalTrackWoDCA,
  TrackSelectionQualityTracks,
  TrackSelectionInAcceptanceTracks
};

enum {
  EventSelectionNo = 0,
  EventSelectionRun2,
  EventSelectionRun3
};

/// Event selection
template <typename CollisionType>
bool isEventSelected(const CollisionType& collision, const int applyEvSel)
{
  if (applyEvSel == EventSelectionRun2) {
    if (!collision.sel7()) {
      return false;
    }
  } else if (applyEvSel == EventSelectionRun3) {
    if (!collision.sel8()) {
      return false;
    }
  }
  return true;
};

} // namespace o2::dpg_tpcskimstablecreator
#endif // DPG_TASKS_TPC_UTILSTPCSKIMSTABLECREATOR_H_
