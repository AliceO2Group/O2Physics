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
/// \brief This header provides the table definitions to store photon meson event tables and ccdb tables
/// \author Marvin Hemmer (marvin.hemmer@cern.ch) - Goethe University Frankfurt

#ifndef PWGEM_PHOTONMESON_DATAMODEL_EVENTTABLES_H_
#define PWGEM_PHOTONMESON_DATAMODEL_EVENTTABLES_H_

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"

#include <DataFormatsParameters/GRPMagField.h>
#include <EMCALCalib/BadChannelMap.h>
#include <Framework/AnalysisDataModel.h>

#include <TBufferFile.h> // IWYU pragma: keep
#include <TClass.h>      // IWYU pragma: keep

#include <Rtypes.h> // for BIT

#include <cstdint>
#include <vector>

namespace o2::aod
{

namespace pmevsel
{

DECLARE_SOA_DYNAMIC_COLUMN(Sel8, sel8, [](uint64_t selection_bit, int runNumber) -> bool {
  return (selection_bit & BIT(o2::aod::evsel::kIsTriggerTVX)) && (selection_bit & BIT(o2::aod::evsel::kNoTimeFrameBorder)) && (runNumber < 568873 ? (selection_bit & BIT(o2::aod::evsel::kNoITSROFrameBorder)) : true); // o2-linter: disable=magic-number (hard-coded run range to indicate 2026 datataking)
});

// Enum used for filling table for event-norm purposes
enum EventAcceptanceBits {
  kAll = 0, // o2-linter: disable=magic-number (enum)
  kHasMCColl,
  kGoodZVtx,
  kIsFT0AND,
  kNoTFB,
  kITSROFB,
  kNoSameBunchPileUp,
  kGoodZVtxFTOPV,
  kNoCollInTimeRange,
  kGoodTrackOccupancy,
  kGoodFT0Occupancy,
  kTVXInEMC,
  kGoodCent,
  kGoodRCT,
  kGoodSel8,
  kSize
};

} // namespace pmevsel

namespace pmevent
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int);

} // namespace pmevent

DECLARE_SOA_TABLE(PMEvents, "AOD", "PMEVENT", //!   Main event information table
                  o2::soa::Index<>, pmevent::CollisionId, bc::RunNumber, bc::GlobalBC, evsel::Selection, evsel::Rct, timestamp::Timestamp,
                  collision::PosZ,
                  collision::NumContrib, evsel::NumTracksInTimeRange, evsel::SumAmpFT0CInTimeRange, pmevsel::Sel8<evsel::Selection, bc::RunNumber>);

using PMEvent = PMEvents::iterator;

// Tables for event selection and event bookkeeping

DECLARE_SOA_COLUMN(IsSelected, isSelected, bool); //! MB event selection info
DECLARE_SOA_TABLE(PMEvSels, "AOD", "PMEVSEL",     //! joinable to o2::aod::Collisions
                  IsSelected);
using PMEvSel = PMEvSels::iterator;

DECLARE_SOA_COLUMN(EventSelectionBit, eventSelectionBit, std::vector<uint64_t>); //! Event selection info stored in binned data for each DF
DECLARE_SOA_TABLE(PMEvSelBits, "AOD", "PMEVSELBITS",                             //! produces binned data that can be loaded in analysis task for event counting
                  EventSelectionBit);
using PMEvSelBit = PMEvSelBits::iterator;

namespace ccdbPcm
{
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
DECLARE_SOA_CCDB_COLUMN(GRPMagField, grpMagField, o2::parameters::GRPMagField, "GLO/Config/GRPMagField"); //!
} // namespace ccdbPcm

/// Full table — join with aod::BCsWithTimestamps to obtain all four objects.
DECLARE_SOA_TIMESTAMPED_TABLE(PcmObjects, aod::Timestamps, o2::aod::timestamp::Timestamp, 0, "PCMOBJECTS", //!
                              ccdbPcm::GRPMagField);

namespace em::ccdbMagField
{
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
DECLARE_SOA_CCDB_COLUMN(GRPMagField, grpMagField, o2::parameters::GRPMagField, "GLO/Config/GRPMagField"); //!
} // namespace em::ccdbMagField

/// Full table — join with aod::BCsWithTimestamps to obtain all four objects.
DECLARE_SOA_TIMESTAMPED_TABLE(EmMagFields, aod::PMEvents, o2::aod::timestamp::Timestamp, 0, "EMMAGFIELDS", //!
                              em::ccdbMagField::GRPMagField);

namespace em::ccdbEmcal
{
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
DECLARE_SOA_CCDB_COLUMN(BadChannelMap, badChannelMap, o2::emcal::BadChannelMap, "EMC/Calib/BadChannelMap"); //! EMCal BadChannelMap for EmEvents
} // namespace em::ccdbEmcal

/// Full table — join with aod::BCsWithTimestamps to obtain all four objects.
DECLARE_SOA_TIMESTAMPED_TABLE(EmEmcalObjects, aod::PMEvents, o2::aod::timestamp::Timestamp, 0, "EMEMCALOBJECTS", //!
                              em::ccdbEmcal::BadChannelMap);

} // namespace o2::aod

#endif // PWGEM_PHOTONMESON_DATAMODEL_EVENTTABLES_H_
