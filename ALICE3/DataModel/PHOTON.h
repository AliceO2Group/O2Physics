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
/// \file   PHOTON.h
/// \note   Based on ECAL.h
/// \author Ana Marin
/// \since  11/10/2021
/// \brief  Set of tables for the ALICE3 Photon Conversions information
///

#ifndef O2_ANALYSIS_ALICE3_PHOTON_H_
#define O2_ANALYSIS_ALICE3_PHOTON_H_

// O2 includes
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace alice3photon
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);   //! Index to travel from PHOTON to collision
DECLARE_SOA_INDEX_COLUMN(Track, track);           //! Index to travel from PHOTON to track
DECLARE_SOA_INDEX_COLUMN(McParticle, mcparticle); //! Index to travel from PHOTON to particle
DECLARE_SOA_COLUMN(PX, px, float);                //! Px
DECLARE_SOA_COLUMN(PY, py, float);                //! Py
DECLARE_SOA_COLUMN(PZ, pz, float);                //! Pz
} // namespace alice3photon

DECLARE_SOA_TABLE(PHOTONs, "AOD", "PHOTONCONV", //! Table for the ALICE3 Photon detector
                  o2::soa::Index<>,
                  alice3photon::CollisionId,
                  alice3photon::TrackId,
                  alice3photon::McParticleId,
                  alice3photon::PX,
                  alice3photon::PY,
                  alice3photon::PZ);

using PHOTON = PHOTONs::iterator;

} // namespace o2::aod

#endif // O2_ANALYSIS_ALICE3_PHOTON_H_
