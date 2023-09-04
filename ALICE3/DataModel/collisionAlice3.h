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
/// \file   collisionExtra.h
/// \author David Dobrigkeit Chinellato
/// \since  11/05/2023
/// \brief  Table for ALICE 3 collision-related info
///

#ifndef ALICE3_DATAMODEL_COLLISIONALICE3_H_
#define ALICE3_DATAMODEL_COLLISIONALICE3_H_

// O2 includes
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace collision_alice3
{
DECLARE_SOA_COLUMN(MultDensity, multDensity, float); //! midrapidity Nch density
} // namespace collision_alice3
DECLARE_SOA_TABLE(CollisionsAlice3, "AOD", "COLLALICE3",
                  collision_alice3::MultDensity);

using CollisionAlice3 = CollisionsAlice3::iterator;

} // namespace o2::aod

#endif // ALICE3_DATAMODEL_COLLISIONALICE3_H_
