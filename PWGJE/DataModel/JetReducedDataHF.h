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
/// \brief Table definitions for reduced data model for hf jets
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_DATAMODEL_JETREDUCEDDATAHF_H_
#define PWGJE_DATAMODEL_JETREDUCEDDATAHF_H_

#include <cmath>
#include <vector>
#include "Framework/AnalysisDataModel.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/JetReducedData.h"

namespace o2::aod
{

namespace jd0indices
{
DECLARE_SOA_INDEX_COLUMN(JCollision, collision);
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, JTracks, "_0");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, JTracks, "_1");
DECLARE_SOA_INDEX_COLUMN(JMcCollision, mcCollision);
DECLARE_SOA_INDEX_COLUMN(JMcParticle, mcParticle);
} // namespace jd0indices

DECLARE_SOA_TABLE_STAGED(JD0CollisionIds, "JD0COLLID",
                         jd0indices::JCollisionId);

DECLARE_SOA_TABLE_STAGED(JD0McCollisionIds, "JD0MCCOLLID",
                         jd0indices::JMcCollisionId);

DECLARE_SOA_TABLE_STAGED(JD0Ids, "JD0ID",
                         jd0indices::JCollisionId,
                         jd0indices::Prong0Id,
                         jd0indices::Prong1Id);

DECLARE_SOA_TABLE_STAGED(JD0PIds, "JD0PID",
                         jd0indices::JMcCollisionId,
                         jd0indices::JMcParticleId);

namespace jlcindices
{
DECLARE_SOA_INDEX_COLUMN(JCollision, collision);
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, JTracks, "_0");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, JTracks, "_1");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong2, prong2, int, JTracks, "_2");
DECLARE_SOA_INDEX_COLUMN(JMcCollision, mcCollision);
DECLARE_SOA_INDEX_COLUMN(JMcParticle, mcParticle);
} // namespace jlcindices

DECLARE_SOA_TABLE_STAGED(JLcCollisionIds, "JLCCOLLID",
                         jlcindices::JCollisionId);

DECLARE_SOA_TABLE_STAGED(JLcMcCollisionIds, "JLCMCCOLLID",
                         jlcindices::JMcCollisionId);

DECLARE_SOA_TABLE_STAGED(JLcIds, "JLCID",
                         jlcindices::JCollisionId,
                         jlcindices::Prong0Id,
                         jlcindices::Prong1Id,
                         jlcindices::Prong2Id);

DECLARE_SOA_TABLE_STAGED(JLcPIds, "JLCPID",
                         jlcindices::JMcCollisionId,
                         jlcindices::JMcParticleId);

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETREDUCEDDATAHF_H_
