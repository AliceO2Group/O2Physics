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

/// \file EventSelection.h
/// \brief Definitions of event selection tables
///
/// \author Evgeny Kryshen <evgeny.kryshen@cern.ch> and Igor Altsybeev <Igor.Altsybeev@cern.ch>

#ifndef COMMON_DATAMODEL_EVENTSELECTION_H_
#define COMMON_DATAMODEL_EVENTSELECTION_H_

#include <Framework/AnalysisDataModel.h>
#include <Framework/Logger.h>

#include <algorithm>
#include <cstdint>
#include <ranges>
#include <vector>

namespace o2::aod
{

// Bits in eventCuts bitmask in Run2BCInfos table
// Must be consistent with EventSelectionCut enum in the Run2 converter
enum Run2EventCuts {
  kINELgtZERO = 0,
  kPileupInMultBins,
  kConsistencySPDandTrackVertices,
  kTrackletsVsClusters,
  kNonZeroNContribs,
  kIncompleteDAQ,
  kPileUpMV,
  kTPCPileUp,
  kTimeRangeCut,
  kEMCALEDCut,
  kAliEventCutsAccepted,
  kIsPileupFromSPD,
  kIsV0PFPileup,
  kIsTPCHVdip,
  kIsTPCLaserWarmUp,
  kTRDHCO, // Offline TRD cosmic trigger decision
  kTRDHJT, // Offline TRD jet trigger decision
  kTRDHSE, // Offline TRD single electron trigger decision
  kTRDHQU, // Offline TRD quarkonium trigger decision
  kTRDHEE  // Offline TRD single-electron-in-EMCAL-acceptance trigger decision
};

namespace evsel
{
DECLARE_SOA_BITMAP_COLUMN(Alias, alias, 32);                                //! Bitmask of fired trigger aliases (see TriggerAliases.h for definitions)
DECLARE_SOA_BITMAP_COLUMN(Selection, selection, 64);                        //! Bitmask of selection flags (see EventSelectionParams.h for definitions)
DECLARE_SOA_BITMAP_COLUMN(Rct, rct, 32);                                    //! Bitmask of RCT flags
DECLARE_SOA_COLUMN(Sel7, sel7, bool);                                       //! Event selection decision based on V0A & V0C
DECLARE_SOA_COLUMN(Sel8, sel8, bool);                                       //! Event selection decision based on TVX
DECLARE_SOA_INDEX_COLUMN_FULL(FoundBC, foundBC, int, BCs, "_foundBC");      //! BC entry index in BCs table (-1 if doesn't exist)
DECLARE_SOA_INDEX_COLUMN_FULL(FoundFT0, foundFT0, int, FT0s, "_foundFT0");  //! FT0 entry index in FT0s table (-1 if doesn't exist)
DECLARE_SOA_INDEX_COLUMN_FULL(FoundFV0, foundFV0, int, FV0As, "_foundFV0"); //! FV0 entry index in FV0As table (-1 if doesn't exist)
DECLARE_SOA_INDEX_COLUMN_FULL(FoundFDD, foundFDD, int, FDDs, "_foundFDD");  //! FDD entry index in FDDs table (-1 if doesn't exist)
DECLARE_SOA_INDEX_COLUMN_FULL(FoundZDC, foundZDC, int, Zdcs, "_foundZDC");  //! ZDC entry index in ZDCs table (-1 if doesn't exist)
DECLARE_SOA_COLUMN(NumTracksInTimeRange, trackOccupancyInTimeRange, int);   //! Occupancy in specified time interval by a number of tracks from nearby collisions // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(SumAmpFT0CInTimeRange, ft0cOccupancyInTimeRange, float); //! Occupancy in specified time interval by a sum of FT0C amplitudes from nearby collisions // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(OccupancyMedianTime, occupancyMedianTime, float);        //! Median time for the occupancy in specified time interval (the delta-time wrt to the collision at which occupancy crosses 1/2 of its full value)
} // namespace evsel

// bc-joinable event selection decisions
DECLARE_SOA_TABLE(BcSels, "AOD", "BCSEL", //!
                  evsel::Alias, evsel::Selection, evsel::Rct, evsel::FoundFT0Id, evsel::FoundFV0Id, evsel::FoundFDDId, evsel::FoundZDCId);
using BcSel = BcSels::iterator;

// collision-joinable event selection decisions
DECLARE_SOA_TABLE(EvSels, "AOD", "EVSEL", //!
                  evsel::Alias,
                  evsel::Selection,
                  evsel::Rct,
                  evsel::Sel7,
                  evsel::Sel8,
                  evsel::FoundBCId,
                  evsel::FoundFT0Id,
                  evsel::FoundFV0Id,
                  evsel::FoundFDDId,
                  evsel::FoundZDCId,
                  evsel::NumTracksInTimeRange,
                  evsel::SumAmpFT0CInTimeRange,
                  evsel::OccupancyMedianTime);
using EvSel = EvSels::iterator;

// occupancy wrt BCs: same estimators and flags as in EvSels (occupancy from other collisions nearby in time),
// calculated wrt a BC; collisions found in this BC are excluded (their number is stored in numCollsInBC)
// N.B.: not related to the occupancy tables in OccupancyTables.h
namespace bcocc
{
enum BcOccupancyFlags {         // bits of occupancyFlags
  kNoCollInTimeRangeNarrow = 0, // no other collisions in specified time range (narrower than Strict)
  kNoCollInTimeRangeStrict,     // no other collisions in specified time range
  kNoCollInTimeRangeStandard,   // no other collisions in specified time range with per-collision multiplicity above threshold
  kNoCollInRofStrict,           // no other collisions in this ITS Readout Frame
  kNoCollInRofStandard,         // no other collisions in this ITS Readout Frame with per-collision multiplicity above threshold (no veto on collisions with close vZ)
  kNoHighMultCollInPrevRof,     // FT0C amplitude in previous ITS ROF is below threshold
  kUpcMode                      // ITS ROF of this BC is in UPC mode: no collisions in the ROF, or collisions reconstructed in the UPC mode (set for all BCs, also without calculated occupancy)
};
DECLARE_SOA_COLUMN(NumTracksInTimeRange, trackOccupancyInTimeRange, int);   //! Occupancy in specified time interval by a number of tracks from nearby collisions (-1 if not available) // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(SumAmpFT0CInTimeRange, ft0cOccupancyInTimeRange, float); //! Occupancy in specified time interval by a sum of FT0C amplitudes from nearby collisions (-1 if not available) // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(OccupancyMedianTime, occupancyMedianTime, float);        //! Median time for the occupancy in specified time interval (the delta-time wrt to the BC at which occupancy crosses 1/2 of its full value)
DECLARE_SOA_BITMAP_COLUMN(OccupancyFlags, occupancyFlags, 8);               //! Bitmask of occupancy flags (see BcOccupancyFlags), check with occupancyFlags_bit(bit)
DECLARE_SOA_COLUMN(NumCollsInBC, numCollsInBC, int16_t);                    //! Number of collisions found in this BC (excluded from the occupancy)
DECLARE_SOA_COLUMN(NumCollsInRof, numCollsInRof, int16_t);                  //! Number of collisions found in the ITS ROF of this BC (within the TF), including the collisions in this BC (-1 if the ROF parameters are not available)
} // namespace bcocc

// bc-joinable occupancy
DECLARE_SOA_TABLE(BcOccs, "AOD", "BCOCC", //!
                  bcocc::NumTracksInTimeRange,
                  bcocc::SumAmpFT0CInTimeRange,
                  bcocc::OccupancyMedianTime,
                  bcocc::OccupancyFlags,
                  bcocc::NumCollsInBC,
                  bcocc::NumCollsInRof);
using BcOcc = BcOccs::iterator;

// occupancy wrt every BC of the TFs present in the DF (not joinable): for each TF (in increasing order) all its BCs, in increasing order of global BC
// same estimators and flags as in BcOccs
DECLARE_SOA_TABLE(BcOccsAll, "AOD", "BCOCCALL", //!
                  bc::GlobalBC,
                  bcocc::NumTracksInTimeRange,
                  bcocc::SumAmpFT0CInTimeRange,
                  bcocc::OccupancyMedianTime,
                  bcocc::OccupancyFlags,
                  bcocc::NumCollsInBC,
                  bcocc::NumCollsInRof);
using BcOccAll = BcOccsAll::iterator;

namespace bcocc
{
// fast lookup of the aod::BcOccsAll entry for a given global BC
// usage (within the process function, the index keeps a pointer to the table):
//   aod::bcocc::BcOccAllIndex bcOccAllIndexed(bcOccsAll);
//   auto occ = bcOccAllIndexed[bc.globalBC()]; // iterator of aod::BcOccsAll
//   int occByTracks = occ.trackOccupancyInTimeRange();
template <typename TBcOccsAll>
class BcOccAllIndex
{
 public:
  explicit BcOccAllIndex(TBcOccsAll const& bcOccsAll) : mTable(&bcOccsAll)
  {
    const int64_t nRows = static_cast<int64_t>(bcOccsAll.size());
    // offset globalBC - row: constant along consecutive BCs, changes only where TFs are missing in the DF
    auto offsetAt = [&bcOccsAll](int64_t row) { return static_cast<int64_t>(bcOccsAll.iteratorAt(row).globalBC()) - row; };
    for (int64_t first = 0, end = 0; first < nRows; first = end) {
      // the offset is non-decreasing (global BCs are strictly increasing): the segment ends at the first row with a larger offset
      const int64_t offset = offsetAt(first);
      // binary search over the rows [first, nRows) for the end of the segment: the first row with a different offset
      end = *std::ranges::partition_point(std::views::iota(first, nRows), [&](int64_t row) { return offsetAt(row) == offset; });
      mSegments.push_back({offset + first, offset + end, offset});
    }
  }

  // entry of aod::BcOccsAll for this global BC (fatal if the BC is not in the table)
  auto operator[](uint64_t globalBC) const
  {
    const int64_t row = rowOf(globalBC);
    if (row < 0) {
      LOGP(fatal, "BcOccAllIndex: global BC {} is not in aod::BcOccsAll ({} rows): its TF is not in this DF, or the occupancy is not available for this DF", globalBC, mTable->size());
    }
    return mTable->iteratorAt(row);
  }

  // helper function to check if a global BC is contained in the table
  bool contains(uint64_t globalBC) const { return rowOf(globalBC) >= 0; }

 private:
  // row of aod::BcOccsAll for this global BC, -1 if the BC is not in the table
  int64_t rowOf(uint64_t globalBC) const
  {
    const auto bc = static_cast<int64_t>(globalBC);
    for (const auto& seg : mSegments) { // a few segments: one per range of consecutive TFs in the DF
      if (bc >= seg.firstGlobalBC && bc < seg.endGlobalBC) {
        return bc - seg.offset;
      }
    }
    return -1;
  }

  struct Segment {
    int64_t firstGlobalBC; // first global BC of the segment
    int64_t endGlobalBC;   // last global BC of the segment + 1
    int64_t offset;        // globalBC - row, constant within the segment
  };
  TBcOccsAll const* mTable;       // not owned: the index is valid as long as the table (within the process function)
  std::vector<Segment> mSegments; // in increasing order of global BC
};
} // namespace bcocc
} // namespace o2::aod

#endif // COMMON_DATAMODEL_EVENTSELECTION_H_
