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

#ifndef PWGMM_MULT_CORE_INCLUDE_SELECTIONS_H_
#define PWGMM_MULT_CORE_INCLUDE_SELECTIONS_H_

#include "Common/DataModel/TrackSelectionTables.h"

namespace pwgmm::mult
{

// default quality criteria for tracks with ITS contribution
static constexpr o2::aod::track::TrackSelectionFlags::flagtype trackSelectionITS =
  o2::aod::track::TrackSelectionFlags::kITSNCls | o2::aod::track::TrackSelectionFlags::kITSChi2NDF |
  o2::aod::track::TrackSelectionFlags::kITSHits;

// default quality criteria for tracks with TPC contribution
static constexpr o2::aod::track::TrackSelectionFlags::flagtype trackSelectionTPC =
  o2::aod::track::TrackSelectionFlags::kTPCNCls |
  o2::aod::track::TrackSelectionFlags::kTPCCrossedRowsOverNCls |
  o2::aod::track::TrackSelectionFlags::kTPCChi2NDF;

// default standard DCA cuts
static constexpr o2::aod::track::TrackSelectionFlags::flagtype trackSelectionDCA =
  o2::aod::track::TrackSelectionFlags::kDCAz | o2::aod::track::TrackSelectionFlags::kDCAxy;

// default standard transversal-only DCA cuts
static constexpr o2::aod::track::TrackSelectionFlags::flagtype trackSelectionDCAXYonly =
  o2::aod::track::TrackSelectionFlags::kDCAxy;
} // namespace pwgmm::mult

#endif // PWGMM_MULT_CORE_INCLUDE_SELECTIONS_H_
