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

#ifndef COMMON_DATAMODEL_STRANGENESSTRACKING_H_
#define COMMON_DATAMODEL_STRANGENESSTRACKING_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
DECLARE_SOA_TABLE(TrackedCascadeColls, "AOD", "TRACASCCOLL", //! Table joinable with TrackedCascades containing collision ids
                  track::CollisionId, o2::soa::Marker<1>);

using TrackedCascadeColl = TrackedCascadeColls::iterator;

DECLARE_SOA_TABLE(TrackedV0Colls, "AOD", "TRAV0COLL", //! Table joinable with TrackedV0s containing collision ids
                  track::CollisionId, o2::soa::Marker<2>);
using TrackedV0Coll = TrackedV0Colls::iterator;

DECLARE_SOA_TABLE(Tracked3BodyColls, "AOD", "TRA3BODYCOLL", //! Table joinable with Tracked3Bodys containing collision ids
                  track::CollisionId, o2::soa::Marker<3>);
using Tracked3BodyColl = Tracked3BodyColls::iterator;
} // namespace o2::aod

#endif // COMMON_DATAMODEL_STRANGENESSTRACKING_H_
