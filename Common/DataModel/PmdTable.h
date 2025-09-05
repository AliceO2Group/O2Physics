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
/// \file PmdTable.h
///
/// \brief pmd index table define
/// \author Abhi Modak (abhi.modak@cern.ch)
/// \since May 17, 2025

#ifndef COMMON_DATAMODEL_PMDTABLE_H_
#define COMMON_DATAMODEL_PMDTABLE_H_

#include <Framework/AnalysisDataModel.h>

namespace o2::aod
{
namespace pmdtrack
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_ARRAY_INDEX_COLUMN(Collision, collisions);
DECLARE_SOA_INDEX_COLUMN(BC, bc);
DECLARE_SOA_SLICE_INDEX_COLUMN(Pmd, pmd);
} // namespace pmdtrack

DECLARE_SOA_INDEX_TABLE_USER(PMDTracksIndex, BCs, "PMDTRKIDX", pmdtrack::CollisionId, pmdtrack::BCId, pmdtrack::PmdIdSlice);
} // namespace o2::aod

#endif // COMMON_DATAMODEL_PMDTABLE_H_
