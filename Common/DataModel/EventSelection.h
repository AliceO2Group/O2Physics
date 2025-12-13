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

#include "Common/CCDB/EventSelectionParams.h" // FIXME: remove
#include "Common/CCDB/RCTSelectionFlags.h"    // FIXME: remove
#include "Common/CCDB/TriggerAliases.h"       // FIXME: remove

#include <Framework/AnalysisDataModel.h>

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
} // namespace o2::aod

#endif // COMMON_DATAMODEL_EVENTSELECTION_H_
