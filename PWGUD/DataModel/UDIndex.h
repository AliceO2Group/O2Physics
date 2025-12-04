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
/// \file   UDIndex.h
/// \author Roman Lavička
/// \since  2025-04-15
/// \brief  A table to store indices for association of MC truth to MC reco
///

#ifndef PWGUD_DATAMODEL_UDINDEX_H_
#define PWGUD_DATAMODEL_UDINDEX_H_

#include "UDTables.h"

#include "Framework/AnalysisDataModel.h"
namespace o2::aod
{
namespace udidx
{
DECLARE_SOA_ARRAY_INDEX_COLUMN(UDTrack, udtracks);
DECLARE_SOA_ARRAY_INDEX_COLUMN(UDCollision, udcollisions);
} // namespace udidx
DECLARE_SOA_TABLE(UDMcParticlesToUDTracks, "AOD", "UDP2UDT", udidx::UDTrackIds);
DECLARE_SOA_TABLE(UDMcCollisionsToUDCollisions, "AOD", "UDMCC2UDC", udidx::UDCollisionIds);
} // namespace o2::aod
#endif // PWGUD_DATAMODEL_UDINDEX_H_
