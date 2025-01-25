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
/// \brief Table definitions for reduced data model for dq jets
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_DATAMODEL_JETREDUCEDDATADQ_H_
#define PWGJE_DATAMODEL_JETREDUCEDDATADQ_H_

#include <cmath>
#include <vector>
#include "Framework/AnalysisDataModel.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGJE/DataModel/JetReducedData.h"

namespace o2::aod
{

DECLARE_SOA_TABLE_STAGED(JDielectronMcCollisions, "JDIELMCCOLL",
                         o2::soa::Index<>,
                         jmccollision::PosX,
                         jmccollision::PosY,
                         jmccollision::PosZ);

namespace jdielectronindices
{
DECLARE_SOA_INDEX_COLUMN(JCollision, collision);
DECLARE_SOA_INDEX_COLUMN_CUSTOM(JDielectronMcCollision, dielectronmccollision, "JDIELMCCOLLS");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, JTracks, "_0");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, JTracks, "_1");
DECLARE_SOA_INDEX_COLUMN(JMcCollision, mcCollision);
DECLARE_SOA_INDEX_COLUMN(JMcParticle, mcParticle);
} // namespace jdielectronindices

DECLARE_SOA_TABLE_STAGED(JDielectronCollisionIds, "JDIELCOLLID",
                         jdielectronindices::JCollisionId);

DECLARE_SOA_TABLE_STAGED(JDielectronMcCollisionIds, "JDIELMCCOLLID",
                         jdielectronindices::JMcCollisionId);

DECLARE_SOA_TABLE_STAGED(JDielectronIds, "JDIELID",
                         jdielectronindices::JCollisionId,
                         jdielectronindices::Prong0Id,
                         jdielectronindices::Prong1Id);

namespace jdielectronmc
{
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(E, e, float);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(PdgCode, pdgCode, int);
DECLARE_SOA_COLUMN(GenStatusCode, getGenStatusCode, int); // TODO : We can look at combining this with the two below
DECLARE_SOA_COLUMN(HepMCStatusCode, getHepMCStatusCode, int);
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool);
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(Mothers, mothers);
DECLARE_SOA_SELF_SLICE_INDEX_COLUMN(Daughters, daughters);
DECLARE_SOA_COLUMN(DecayFlag, decayFlag, int8_t);
DECLARE_SOA_COLUMN(Origin, origin, int);
DECLARE_SOA_DYNAMIC_COLUMN(Px, px,
                           [](float pt, float phi) -> float { return pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py,
                           [](float pt, float phi) -> float { return pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz,
                           [](float pt, float eta) -> float { return pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p,
                           [](float pt, float eta) -> float { return pt * std::cosh(eta); });
} // namespace jdielectronmc

DECLARE_SOA_TABLE_STAGED(JDielectronMcs, "JDIELMC",
                         o2::soa::Index<>,
                         jdielectronindices::JDielectronMcCollisionId,
                         jdielectronmc::Pt,
                         jdielectronmc::Eta,
                         jdielectronmc::Phi,
                         jdielectronmc::Y,
                         jdielectronmc::E,
                         jdielectronmc::M,
                         jdielectronmc::PdgCode,
                         jdielectronmc::GenStatusCode,
                         jdielectronmc::HepMCStatusCode,
                         jdielectronmc::IsPhysicalPrimary,
                         jdielectronmc::DecayFlag,
                         jdielectronmc::Origin,
                         jdielectronmc::Px<jdielectronmc::Pt, jdielectronmc::Phi>,
                         jdielectronmc::Py<jdielectronmc::Pt, jdielectronmc::Phi>,
                         jdielectronmc::Pz<jdielectronmc::Pt, jdielectronmc::Eta>,
                         jdielectronmc::P<jdielectronmc::Pt, jdielectronmc::Eta>);

using JDielectronMc = JDielectronMcs::iterator;
using StoredJDielectronMc = StoredJDielectronMcs::iterator;

DECLARE_SOA_TABLE_STAGED(JDielectronMcIds, "JDIELMCID",
                         jdielectronindices::JMcCollisionId,
                         jdielectronindices::JMcParticleId,
                         jdielectronmc::MothersIds,
                         jdielectronmc::DaughtersIdSlice);

namespace jdummydq
{

DECLARE_SOA_COLUMN(DummyDQ, dummyDQ, bool);

} // namespace jdummydq
DECLARE_SOA_TABLE(JDielectron1Dummys, "AOD", "JDIEL1DUMMY",
                  jdummydq::DummyDQ,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JDielectron2Dummys, "AOD", "JDIEL2DUMMY",
                  jdummydq::DummyDQ,
                  o2::soa::Marker<2>);

DECLARE_SOA_TABLE(JDielectron3Dummys, "AOD", "JDIEL3DUMMY",
                  jdummydq::DummyDQ,
                  o2::soa::Marker<3>);

DECLARE_SOA_TABLE(JDielectron4Dummys, "AOD", "JDIEL4DUMMY",
                  jdummydq::DummyDQ,
                  o2::soa::Marker<4>);

DECLARE_SOA_TABLE(JDielectron5Dummys, "AOD", "JDIEL5DUMMY",
                  jdummydq::DummyDQ,
                  o2::soa::Marker<5>);

DECLARE_SOA_TABLE(JDielectron6Dummys, "AOD", "JDIEL6DUMMY",
                  jdummydq::DummyDQ,
                  o2::soa::Marker<6>);

DECLARE_SOA_TABLE(JDielectron7Dummys, "AOD", "JDIEL7DUMMY",
                  jdummydq::DummyDQ,
                  o2::soa::Marker<7>);

DECLARE_SOA_TABLE(JDielectron8Dummys, "AOD", "JDIEL8DUMMY",
                  jdummydq::DummyDQ,
                  o2::soa::Marker<8>);

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETREDUCEDDATADQ_H_
