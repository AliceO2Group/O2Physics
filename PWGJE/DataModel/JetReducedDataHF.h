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

#include "PWGJE/DataModel/EMCALClusters.h" // IWYU pragma: keep
#include "PWGJE/DataModel/JetReducedData.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h> // IWYU pragma: keep

#include <sys/types.h>

#include <cmath>

namespace o2::aod
{

constexpr uint JMarkerD0 = 1;
constexpr uint JMarkerDplus = 2;
constexpr uint JMarkerLc = 3;
constexpr uint JMarkerBplus = 4;
constexpr uint JMarkerDielectron = 5;
constexpr uint JMarkerDstar = 6;
constexpr uint JMarkerB0 = 7;
constexpr uint JMarkerDs = 8;
constexpr uint JMarkerXicToXiPiPi = 9;

namespace jcandidateindices
{
DECLARE_SOA_INDEX_COLUMN(JCollision, collision);
DECLARE_SOA_INDEX_COLUMN(JMcCollision, mcCollision);
DECLARE_SOA_INDEX_COLUMN(JMcParticle, mcParticle);
} // namespace jcandidateindices

namespace jd0indices
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, JTracks, "_0");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, JTracks, "_1");
} // namespace jd0indices

DECLARE_SOA_TABLE_STAGED(JD0CollisionIds, "JD0COLLID",
                         jcandidateindices::JCollisionId,
                         o2::soa::Marker<JMarkerD0>);

DECLARE_SOA_TABLE_STAGED(JD0McCollisionIds, "JD0MCCOLLID",
                         jcandidateindices::JMcCollisionId,
                         o2::soa::Marker<JMarkerD0>);

DECLARE_SOA_TABLE_STAGED(JD0Ids, "JD0ID",
                         jcandidateindices::JCollisionId,
                         jd0indices::Prong0Id,
                         jd0indices::Prong1Id);

DECLARE_SOA_TABLE_STAGED(JD0PIds, "JD0PID",
                         jcandidateindices::JMcCollisionId,
                         jcandidateindices::JMcParticleId,
                         o2::soa::Marker<JMarkerD0>);

namespace jdummyd0
{
DECLARE_SOA_COLUMN(DummyD0, dummyD0, bool);
} // namespace jdummyd0

DECLARE_SOA_TABLE(JDumD0ParDaus, "AOD", "JDUMD0PARDAU",
                  jdummyd0::DummyD0,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JDumD0MlDaus, "AOD", "JDumD0MLDAU",
                  jdummyd0::DummyD0,
                  o2::soa::Marker<2>);

namespace jdplusindices
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, JTracks, "_0");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, JTracks, "_1");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong2, prong2, int, JTracks, "_2");
} // namespace jdplusindices

DECLARE_SOA_TABLE_STAGED(JDplusCollisionIds, "JDPCOLLID",
                         jcandidateindices::JCollisionId,
                         o2::soa::Marker<JMarkerDplus>);

DECLARE_SOA_TABLE_STAGED(JDplusMcCollisionIds, "JDPMCCOLLID",
                         jcandidateindices::JMcCollisionId,
                         o2::soa::Marker<JMarkerDplus>);

DECLARE_SOA_TABLE_STAGED(JDplusIds, "JDPID",
                         jcandidateindices::JCollisionId,
                         jdplusindices::Prong0Id,
                         jdplusindices::Prong1Id,
                         jdplusindices::Prong2Id);

DECLARE_SOA_TABLE_STAGED(JDplusPIds, "JDPPID",
                         jcandidateindices::JMcCollisionId,
                         jcandidateindices::JMcParticleId,
                         o2::soa::Marker<JMarkerDplus>);

namespace jdummydplus
{

DECLARE_SOA_COLUMN(DummyDplus, dummyDplus, bool);

} // namespace jdummydplus
DECLARE_SOA_TABLE(JDumDplusParDaus, "AOD", "JDUMDPPARDAU",
                  jdummydplus::DummyDplus,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JDumDplusMlDaus, "AOD", "JDUMDPMLDAU",
                  jdummydplus::DummyDplus,
                  o2::soa::Marker<2>);

// might have to update!!!!!!!!!!!!!!!!!!!!!!

namespace jdsindices
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, JTracks, "_0");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, JTracks, "_1");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong2, prong2, int, JTracks, "_2");
} // namespace jdsindices

DECLARE_SOA_TABLE_STAGED(JDsCollisionIds, "JDSCOLLID",
                         jcandidateindices::JCollisionId,
                         o2::soa::Marker<JMarkerDs>);

DECLARE_SOA_TABLE_STAGED(JDsMcCollisionIds, "JDSMCCOLLID",
                         jcandidateindices::JMcCollisionId,
                         o2::soa::Marker<JMarkerDs>);

DECLARE_SOA_TABLE_STAGED(JDsIds, "JDSID",
                         jcandidateindices::JCollisionId,
                         jdsindices::Prong0Id,
                         jdsindices::Prong1Id,
                         jdsindices::Prong2Id);

DECLARE_SOA_TABLE_STAGED(JDsPIds, "JDSPID",
                         jcandidateindices::JMcCollisionId,
                         jcandidateindices::JMcParticleId,
                         o2::soa::Marker<JMarkerDs>);

namespace jdummyds
{

DECLARE_SOA_COLUMN(DummyDs, dummyDs, bool);

} // namespace jdummyds
DECLARE_SOA_TABLE(JDumDsParDaus, "AOD", "JDUMDSPARDAU",
                  jdummyds::DummyDs,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JDumDsMlDaus, "AOD", "JDUMDSMLDAU",
                  jdummyds::DummyDs,
                  o2::soa::Marker<2>);

namespace jdstarindices
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, JTracks, "_0");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, JTracks, "_1");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong2, prong2, int, JTracks, "_2");
} // namespace jdstarindices

DECLARE_SOA_TABLE_STAGED(JDstarCollisionIds, "JDSTCOLLID",
                         jcandidateindices::JCollisionId,
                         o2::soa::Marker<JMarkerDstar>);

DECLARE_SOA_TABLE_STAGED(JDstarMcCollisionIds, "JDSTMCCOLLID",
                         jcandidateindices::JMcCollisionId,
                         o2::soa::Marker<JMarkerDstar>);

DECLARE_SOA_TABLE_STAGED(JDstarIds, "JDSTID",
                         jcandidateindices::JCollisionId,
                         jdstarindices::Prong0Id,
                         jdstarindices::Prong1Id,
                         jdstarindices::Prong2Id);

DECLARE_SOA_TABLE_STAGED(JDstarPIds, "JDSTPID",
                         jcandidateindices::JMcCollisionId,
                         jcandidateindices::JMcParticleId,
                         o2::soa::Marker<JMarkerDstar>);

namespace jdummydstar
{

DECLARE_SOA_COLUMN(DummyDstar, dummyDstar, bool);

} // namespace jdummydstar
DECLARE_SOA_TABLE(JDumDstarParEs, "AOD", "JDUMDSTPARE",
                  jdummydstar::DummyDstar,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JDumDstarMlDaus, "AOD", "JDUMDSTMLDAU",
                  jdummydstar::DummyDstar,
                  o2::soa::Marker<2>);

namespace jlcindices
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, JTracks, "_0");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, JTracks, "_1");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong2, prong2, int, JTracks, "_2");
} // namespace jlcindices

DECLARE_SOA_TABLE_STAGED(JLcCollisionIds, "JLCCOLLID",
                         jcandidateindices::JCollisionId,
                         o2::soa::Marker<JMarkerLc>);

DECLARE_SOA_TABLE_STAGED(JLcMcCollisionIds, "JLCMCCOLLID",
                         jcandidateindices::JMcCollisionId,
                         o2::soa::Marker<JMarkerLc>);

DECLARE_SOA_TABLE_STAGED(JLcIds, "JLCID",
                         jcandidateindices::JCollisionId,
                         jlcindices::Prong0Id,
                         jlcindices::Prong1Id,
                         jlcindices::Prong2Id);

DECLARE_SOA_TABLE_STAGED(JLcPIds, "JLCPID",
                         jcandidateindices::JMcCollisionId,
                         jcandidateindices::JMcParticleId,
                         o2::soa::Marker<JMarkerLc>);

namespace jdummylc
{

DECLARE_SOA_COLUMN(DummyLc, dummyLc, bool);

} // namespace jdummylc
DECLARE_SOA_TABLE(JDumLcParDaus, "AOD", "JDUMLCPARDAU",
                  jdummylc::DummyLc,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JDumLcMlDaus, "AOD", "JDUMLCMLDAU",
                  jdummylc::DummyLc,
                  o2::soa::Marker<2>);

namespace jb0indices
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, JTracks, "_0");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, JTracks, "_1");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong2, prong2, int, JTracks, "_2");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong3, prong3, int, JTracks, "_3");
} // namespace jb0indices

DECLARE_SOA_TABLE_STAGED(JB0CollisionIds, "JB0COLLID",
                         jcandidateindices::JCollisionId,
                         o2::soa::Marker<JMarkerB0>);

DECLARE_SOA_TABLE_STAGED(JB0McCollisionIds, "JB0MCCOLLID",
                         jcandidateindices::JMcCollisionId,
                         o2::soa::Marker<JMarkerB0>);

DECLARE_SOA_TABLE_STAGED(JB0Ids, "JB0ID",
                         jcandidateindices::JCollisionId,
                         jb0indices::Prong0Id,
                         jb0indices::Prong1Id,
                         jb0indices::Prong2Id,
                         jb0indices::Prong3Id);

DECLARE_SOA_TABLE_STAGED(JB0PIds, "JB0PID",
                         jcandidateindices::JMcCollisionId,
                         jcandidateindices::JMcParticleId,
                         o2::soa::Marker<JMarkerB0>);

namespace jbplusindices
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, JTracks, "_0");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, JTracks, "_1");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong2, prong2, int, JTracks, "_2");
} // namespace jbplusindices

DECLARE_SOA_TABLE_STAGED(JBplusCollisionIds, "JBPCOLLID",
                         jcandidateindices::JCollisionId,
                         o2::soa::Marker<JMarkerBplus>);

DECLARE_SOA_TABLE_STAGED(JBplusMcCollisionIds, "JBPMCCOLLID",
                         jcandidateindices::JMcCollisionId,
                         o2::soa::Marker<JMarkerBplus>);

DECLARE_SOA_TABLE_STAGED(JBplusIds, "JBPID",
                         jcandidateindices::JCollisionId,
                         jbplusindices::Prong0Id,
                         jbplusindices::Prong1Id,
                         jbplusindices::Prong2Id);

DECLARE_SOA_TABLE_STAGED(JBplusPIds, "JBPPID",
                         jcandidateindices::JMcCollisionId,
                         jcandidateindices::JMcParticleId,
                         o2::soa::Marker<JMarkerBplus>);

namespace jxictoxipipiindices
{
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, JTracks, "_0");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, JTracks, "_1");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong2, prong2, int, JTracks, "_2");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong3, prong3, int, JTracks, "_3");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong4, prong4, int, JTracks, "_4");
} // namespace jxictoxipipiindices

DECLARE_SOA_TABLE_STAGED(JXicToXiPiPiCollisionIds, "JXICXPPCOLLID",
                         jcandidateindices::JCollisionId,
                         o2::soa::Marker<JMarkerXicToXiPiPi>);

DECLARE_SOA_TABLE_STAGED(JXicToXiPiPiMcCollisionIds, "JXICXPPMCCOLLID",
                         jcandidateindices::JMcCollisionId,
                         o2::soa::Marker<JMarkerXicToXiPiPi>);

DECLARE_SOA_TABLE_STAGED(JXicToXiPiPiIds, "JXICXPPID",
                         jcandidateindices::JCollisionId,
                         jxictoxipipiindices::Prong0Id,
                         jxictoxipipiindices::Prong1Id,
                         jxictoxipipiindices::Prong2Id,
                         jxictoxipipiindices::Prong3Id,
                         jxictoxipipiindices::Prong4Id);

DECLARE_SOA_TABLE_STAGED(JXicToXiPiPiPIds, "JXICXPPPID",
                         jcandidateindices::JMcCollisionId,
                         jcandidateindices::JMcParticleId,
                         o2::soa::Marker<JMarkerXicToXiPiPi>);

namespace jdummyxictoxipipi
{
DECLARE_SOA_COLUMN(DummyXicToXiPiPi, dummyXicToXiPiPi, bool);

} // namespace jdummyxictoxipipi
DECLARE_SOA_TABLE(JDumXicToXiPiPiParDaus, "AOD", "JDUMXICXPPPARDAU",
                  jdummyxictoxipipi::DummyXicToXiPiPi,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JDumXicToXiPiPiMlDaus, "AOD", "JDUMXICXPPMLDAU",
                  jdummyxictoxipipi::DummyXicToXiPiPi,
                  o2::soa::Marker<2>);

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETREDUCEDDATAHF_H_
