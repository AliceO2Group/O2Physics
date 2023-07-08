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

///
/// \file   A3DecayFinderTables.h
/// \since  04/07/2023
/// \brief  Set of tables for ALICE 3 decay finder
///

#ifndef ALICE3_DATAMODEL_A3DECAYFINDERTABLES_H_
#define ALICE3_DATAMODEL_A3DECAYFINDERTABLES_H_

// O2 includes
#include "Framework/AnalysisDataModel.h"

enum a3selectionBit : uint32_t { kDCAxy = 0,
                                 kInnerTOFPion,
                                 kInnerTOFKaon,
                                 kInnerTOFProton,
                                 kOuterTOFPion,
                                 kOuterTOFKaon,
                                 kOuterTOFProton,
                                 kRICHPion,
                                 kRICHKaon,
                                 kRICHProton,
                                 kTruePion,
                                 kTrueKaon,
                                 kTrueProton,
                                 kTruePiPlusFromD,
                                 kTrueKaPlusFromD,
                                 kTruePiMinusFromD,
                                 kTrueKaMinusFromD,
                                 kTruePiPlusFromLc,
                                 kTrueKaPlusFromLc,
                                 kTruePrPlusFromLc,
                                 kTruePiMinusFromLc,
                                 kTrueKaMinusFromLc,
                                 kTruePrMinusFromLc };

namespace o2::aod
{
namespace a3DecayMap
{
DECLARE_SOA_COLUMN(DecayMap, decayMap, uint32_t); //! simple map to process passing / not passing criteria
} // namespace a3DecayMap
DECLARE_SOA_TABLE(Alice3DecayMaps, "AOD", "ALICE3DECAYMAP",
                  a3DecayMap::DecayMap);

using Alice3DecayMap = Alice3DecayMaps::iterator;

} // namespace o2::aod

#endif // ALICE3_DATAMODEL_A3DECAYFINDERTABLES_H_
