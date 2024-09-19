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

DECLARE_SOA_TABLE(JD0CollisionIds, "AOD", "JD0COLLID",
                  jd0indices::JCollisionId);

DECLARE_SOA_TABLE(StoredJD0CollisionIds, "AOD1", "JD0COLLID",
                  jd0indices::JCollisionId,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JD0McCollisionIds, "AOD", "JD0MCCOLLID",
                  jd0indices::JMcCollisionId);

DECLARE_SOA_TABLE(StoredJD0McCollisionIds, "AOD1", "JD0MCCOLLID",
                  jd0indices::JMcCollisionId,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JD0Ids, "AOD", "JD0ID",
                  jd0indices::JCollisionId,
                  jd0indices::Prong0Id,
                  jd0indices::Prong1Id);

DECLARE_SOA_TABLE(StoredJD0Ids, "AOD1", "JD0ID",
                  jd0indices::JCollisionId,
                  jd0indices::Prong0Id,
                  jd0indices::Prong1Id,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JD0PIds, "AOD", "JD0PID",
                  jd0indices::JMcCollisionId,
                  jd0indices::JMcParticleId);

DECLARE_SOA_TABLE(StoredJD0PIds, "AOD1", "JD0PID",
                  jd0indices::JMcCollisionId,
                  jd0indices::JMcParticleId,
                  o2::soa::Marker<1>);

namespace jlcindices
{
DECLARE_SOA_INDEX_COLUMN(JCollision, collision);
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, JTracks, "_0");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, JTracks, "_1");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong2, prong2, int, JTracks, "_2");
DECLARE_SOA_INDEX_COLUMN(JMcCollision, mcCollision);
DECLARE_SOA_INDEX_COLUMN(JMcParticle, mcParticle);
} // namespace jlcindices

DECLARE_SOA_TABLE(JLcCollisionIds, "AOD", "JLCCOLLID",
                  jlcindices::JCollisionId);

DECLARE_SOA_TABLE(StoredJLcCollisionIds, "AOD1", "JLCCOLLID",
                  jlcindices::JCollisionId,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JLcMcCollisionIds, "AOD", "JLCMCCOLLID",
                  jlcindices::JMcCollisionId);

DECLARE_SOA_TABLE(StoredJLcMcCollisionIds, "AOD1", "JLCMCCOLLID",
                  jlcindices::JMcCollisionId,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JLcIds, "AOD", "JLCID",
                  jlcindices::JCollisionId,
                  jlcindices::Prong0Id,
                  jlcindices::Prong1Id,
                  jlcindices::Prong2Id);

DECLARE_SOA_TABLE(StoredJLcIds, "AOD1", "JLCID",
                  jlcindices::JCollisionId,
                  jlcindices::Prong0Id,
                  jlcindices::Prong1Id,
                  jlcindices::Prong2Id,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JLcPIds, "AOD", "JLCPID",
                  jlcindices::JMcCollisionId,
                  jlcindices::JMcParticleId);

DECLARE_SOA_TABLE(StoredJLcPIds, "AOD1", "JLCPID",
                  jlcindices::JMcCollisionId,
                  jlcindices::JMcParticleId,
                  o2::soa::Marker<1>);

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETREDUCEDDATAHF_H_
