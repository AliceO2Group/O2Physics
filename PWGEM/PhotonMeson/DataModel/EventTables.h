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

/// \file EventTables.h
/// \brief This header provides the table definitions to store photon meson event tables
/// \author Marvin Hemmer (marvin.hemmer@cern.ch) - Goethe University Frankfurt

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>

#include <Rtypes.h>

#include <cmath>
#include <cstdint>

#ifndef PWGEM_PHOTONMESON_DATAMODEL_EVENTTABLES_H_
#define PWGEM_PHOTONMESON_DATAMODEL_EVENTTABLES_H_

namespace o2::aod
{

namespace pmevsel
{
// Event selection criteria. See O2Physics/Common/CCDB/EventSelectionParams.h
enum EventSelectionFlags {
  kIsTriggerTVX = 0,          // FT0 vertex (acceptable FT0C-FT0A time difference) at trigger level
  kNoITSROFrameBorder,        // bunch crossing is far from ITS RO Frame border
  kNoTimeFrameBorder,         // bunch crossing is far from Time Frame borders
  kNoSameBunchPileup,         // reject collisions in case of pileup with another collision in the same foundBC
  kIsGoodZvtxFT0vsPV,         // small difference between z-vertex from PV and from FT0
  kIsVertexITSTPC,            // at least one ITS-TPC track (reject vertices built from ITS-only tracks)
  kIsVertexTOFmatched,        // at least one of vertex contributors is matched to TOF
  kIsVertexTRDmatched,        // at least one of vertex contributors is matched to TRD
  kNoCollInTimeRangeNarrow,   // no other collisions in specified time range (narrower than Strict)
  kNoCollInTimeRangeStrict,   // no other collisions in specified time range
  kNoCollInTimeRangeStandard, // no other collisions in specified time range with per-collision multiplicity above threshold
  kNoCollInRofStrict,         // no other collisions in this Readout Frame
  kNoCollInRofStandard,       // no other collisions in this Readout Frame with per-collision multiplicity above threshold
  kNoHighMultCollInPrevRof,   // veto an event if FT0C amplitude in previous ITS ROF is above threshold
  kIsGoodITSLayer3,           // number of inactive chips on ITS layer 3 is below maximum allowed value
  kIsGoodITSLayer0123,        // numbers of inactive chips on ITS layers 0-3 are below maximum allowed values
  kIsGoodITSLayersAll,        // numbers of inactive chips on all ITS layers are below maximum allowed values
  kNsel                       // counter
};

DECLARE_SOA_BITMAP_COLUMN(Selection, selection, 32); //! Bitmask of selection flags
DECLARE_SOA_DYNAMIC_COLUMN(Sel8, sel8, [](uint32_t selection_bit) -> bool { return (selection_bit & BIT(o2::aod::pmevsel::kIsTriggerTVX)) && (selection_bit & BIT(o2::aod::pmevsel::kNoTimeFrameBorder)) && (selection_bit & BIT(o2::aod::pmevsel::kNoITSROFrameBorder)); });

} // namespace pmevsel

namespace pmevent
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int);

DECLARE_SOA_DYNAMIC_COLUMN(Sel8, sel8, [](uint64_t selection_bit) -> bool { return (selection_bit & BIT(o2::aod::evsel::kIsTriggerTVX)) && (selection_bit & BIT(o2::aod::evsel::kNoTimeFrameBorder)) && (selection_bit & BIT(o2::aod::evsel::kNoITSROFrameBorder)); });
} // namespace pmevent

DECLARE_SOA_TABLE(PMEvents, "AOD", "PMEVENT", //!   Main event information table
                  o2::soa::Index<>, pmevent::CollisionId, bc::RunNumber, bc::GlobalBC, evsel::Selection, evsel::Rct, timestamp::Timestamp,
                  collision::PosZ,
                  collision::NumContrib, evsel::NumTracksInTimeRange, evsel::SumAmpFT0CInTimeRange, pmevent::Sel8<evsel::Selection>);

using PMEvent = PMEvents::iterator;

} // namespace o2::aod

#endif // PWGEM_PHOTONMESON_DATAMODEL_EVENTTABLES_H_
