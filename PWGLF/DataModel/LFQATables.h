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

#ifndef PWGLF_DATAMODEL_LFQATABLES_H_
#define PWGLF_DATAMODEL_LFQATABLES_H_

#include <cmath>
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace mccollisionprop
{
DECLARE_SOA_COLUMN(HasRecoCollision, hasRecoCollision, int);     //! stores N times this PV was recoed
DECLARE_SOA_COLUMN(BestCollisionIndex, bestCollisionIndex, int); //! stores N times this PV was recoed
} // namespace mccollisionprop
DECLARE_SOA_TABLE(McCollsExtra, "AOD", "MCCOLLSEXTRA",
                  mccollisionprop::HasRecoCollision, mccollisionprop::BestCollisionIndex);
} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFQATABLES_H_