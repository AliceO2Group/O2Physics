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

// Helper tables to correlate collisions if needed

#ifndef COMMON_DATAMODEL_MCCOLLISIONEXTRA_H_
#define COMMON_DATAMODEL_MCCOLLISIONEXTRA_H_

#include <Framework/AnalysisDataModel.h>

#include <cmath>
#include <cstdint>

namespace o2::aod
{
namespace mccollisionprop
{
DECLARE_SOA_COLUMN(NumRecoCollision, numRecoCollision, int);             //! stores N times this PV was recoed
DECLARE_SOA_COLUMN(BestCollisionIndex, bestCollisionIndex, int);         //! stores N times this PV was recoed
DECLARE_SOA_COLUMN(BestCollisionCentFT0C, bestCollisionCentFT0C, float); //! stores best FT0C centrality

// collision MC context (neighbours contain PoI?)
DECLARE_SOA_COLUMN(ForwardCollisionMap, forwardCollisionMap, uint32_t);   //! stores bitmap telling if PoI found in collisions after this one (bits forward in time)
DECLARE_SOA_COLUMN(BackwardCollisionMap, backwardCollisionMap, uint32_t); //! stores bitmap telling if PoI found in collisions before this one (bits backward in time)
} // namespace mccollisionprop
DECLARE_SOA_TABLE(McCollsExtra, "AOD", "MCCOLLSEXTRA",
                  mccollisionprop::NumRecoCollision, mccollisionprop::BestCollisionIndex, mccollisionprop::BestCollisionCentFT0C);
DECLARE_SOA_TABLE(McCollContexts, "AOD", "MCCOLLCONTEXT",
                  mccollisionprop::ForwardCollisionMap, mccollisionprop::BackwardCollisionMap);
} // namespace o2::aod

#endif // COMMON_DATAMODEL_MCCOLLISIONEXTRA_H_
