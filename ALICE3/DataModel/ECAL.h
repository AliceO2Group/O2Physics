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
/// \file   ECAL.h
/// \author Nicolo' Jacazio
/// \since  14/09/2021
/// \brief  Set of tables for the ALICE3 ECAL information
///

#ifndef O2_ANALYSIS_ALICE3_ECAL_H_
#define O2_ANALYSIS_ALICE3_ECAL_H_

// O2 includes
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace alice3ecal
{
DECLARE_SOA_INDEX_COLUMN(Track, track);    //! Index to travel from track to ECAL
DECLARE_SOA_COLUMN(Energy, energy, float); //! Signal in ECAL
DECLARE_SOA_COLUMN(PosX, posX, float);     //! Error on the ECAL signal
DECLARE_SOA_COLUMN(PosY, posY, float);     //! signal - exp. signal for electrons
DECLARE_SOA_COLUMN(PosZ, posZ, float);     //! signal - exp. signal for muons
} // namespace alice3ecal

DECLARE_SOA_TABLE(ECALs, "AOD", "A3ECAL", //! Table for the ALICE3 ECAL detector
                  o2::soa::Index<>,
                  alice3ecal::TrackId,
                  alice3ecal::Energy,
                  alice3ecal::PosX,
                  alice3ecal::PosY,
                  alice3ecal::PosZ);

using ECAL = ECALs::iterator;

} // namespace o2::aod

#endif // O2_ANALYSIS_ALICE3_ECAL_H_
