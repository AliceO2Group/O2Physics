// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
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
/// \brief  Table definitions for defining whether a collision is matched to a buch crossing including either one or multiple collisions
/// \file   EMCALMatchedCollisions.h
/// \author Nicolas Strangmann <nicolas.strangmann@cern.ch>
/// \since  2023-12-16

#ifndef PWGJE_DATAMODEL_EMCALMATCHEDCOLLISIONS_H_
#define PWGJE_DATAMODEL_EMCALMATCHEDCOLLISIONS_H_

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h> // IWYU pragma: keep

namespace o2::aod
{
namespace emcalcollisionmatch
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);       //! collisionID used as index for matched collisions
DECLARE_SOA_COLUMN(Ambiguous, ambiguous, bool);       //! boolean stating whether the collision is ambiguous (in a BC with multiple collisions)
DECLARE_SOA_COLUMN(IsEMCReadout, isemcreadout, bool); //! boolean stating whether the EMCal was readout in that collision (based on whether the BC contains at least one cell)
} // namespace emcalcollisionmatch

DECLARE_SOA_TABLE(EMCALMatchedCollisions, "AOD", "EMCALMCS",                                                                              //!
                  o2::soa::Index<>, emcalcollisionmatch::CollisionId, emcalcollisionmatch::Ambiguous, emcalcollisionmatch::IsEMCReadout); //

using EMCALMatchedCollision = EMCALMatchedCollisions::iterator;

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_EMCALMATCHEDCOLLISIONS_H_
