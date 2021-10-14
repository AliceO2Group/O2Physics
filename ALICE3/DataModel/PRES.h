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
/// \file   PRES.h
/// \author A. Marin
/// \note   Based on MID.h
/// \brief  Set of tables for the ALICE3 PREShower information
///

#ifndef O2_ANALYSIS_ALICE3_PRES_H_
#define O2_ANALYSIS_ALICE3_PRES_H_

// O2 includes
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace alice3pres
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);    //!
DECLARE_SOA_INDEX_COLUMN(Track, track);            //!
DECLARE_SOA_COLUMN(PRESIsElectron, presIsElectron, uint8_t); //! FIXME: To be changed to bool once bool columns are groupable.
} // namespace alice3mid

DECLARE_SOA_TABLE(PRESs, "AOD", "PRES", //!
                  o2::soa::Index<>,
                  alice3pres::CollisionId,
                  alice3pres::TrackId,
                  alice3pres::PRESIsElectron);

using PRES = PRESs::iterator;

} // namespace o2::aod

#endif // O2_ANALYSIS_ALICE3_PRES_H_
