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
DECLARE_SOA_INDEX_COLUMN(Collision, collision);   //! Index to travel from ECAL to collision
DECLARE_SOA_INDEX_COLUMN(Track, track);           //! Index to travel from ECAL to track
DECLARE_SOA_INDEX_COLUMN(McParticle, mcparticle); //! Index to travel from ECAL to particle
DECLARE_SOA_COLUMN(E, e, double);                 //! Signal in ECAL
DECLARE_SOA_COLUMN(Px, px, double);               //! Px
DECLARE_SOA_COLUMN(Py, py, double);               //! Py
DECLARE_SOA_COLUMN(Pz, pz, double);               //! Pz
DECLARE_SOA_COLUMN(PosZ, posZ, float);            //! Position in Z
DECLARE_SOA_COLUMN(PosPhi, posPhi, float);        //! Position in phi
} // namespace alice3ecal

DECLARE_SOA_TABLE(ECALs, "AOD", "A3ECAL", //! Table for the ALICE3 ECAL detector
                  o2::soa::Index<>,
                  alice3ecal::CollisionId,
                  alice3ecal::TrackId,
                  alice3ecal::McParticleId,
                  alice3ecal::E,
                  alice3ecal::Px,
                  alice3ecal::Py,
                  alice3ecal::Pz,
                  alice3ecal::PosPhi,
                  alice3ecal::PosZ);

using ECAL = ECALs::iterator;

} // namespace o2::aod

#endif // O2_ANALYSIS_ALICE3_ECAL_H_
