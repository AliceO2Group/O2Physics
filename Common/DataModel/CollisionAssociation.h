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
#ifndef COMMON_DATAMODEL_COLLISIONASSOCIATION_H_
#define COMMON_DATAMODEL_COLLISIONASSOCIATION_H_

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{

namespace track_association
{

enum TrackTypes { Regular = 0,
                  Ambiguous,
                  PVContributor };

DECLARE_SOA_INDEX_COLUMN(Collision, collision);   //! Collision index
DECLARE_SOA_INDEX_COLUMN(Track, track);           //! Track index
DECLARE_SOA_COLUMN(TrackType, trackType, int8_t); //! Track type
} // namespace track_association

DECLARE_SOA_TABLE(TrackAssoc, "AOD", "TRACKASSOC", //! Table for track-to-collision association for e.g. HF vertex finding - tracks can appear for several collisions
                  track_association::CollisionId,
                  track_association::TrackId);

DECLARE_SOA_TABLE(TrackAssocExtra, "AOD", "TRACKASSOCEX", //!
                  track_association::TrackType);
} // namespace o2::aod

#endif // COMMON_DATAMODEL_COLLISIONASSOCIATION_H_
