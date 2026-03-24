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

#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetReducedDataHF.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h> // IWYU pragma: keep

#include <cmath>
#include <cstdint>

namespace o2::aod
{
namespace jdielectronmccollision
{
DECLARE_SOA_COLUMN(DummyDQ, dummyDQ, bool);
} // namespace jdielectronmccollision

DECLARE_SOA_TABLE_STAGED(JDielectronMcCollisions, "JDIELMCCOLL",
                         o2::soa::Index<>,
                         jmccollision::PosX,
                         jmccollision::PosY,
                         jmccollision::PosZ);

DECLARE_SOA_TABLE_STAGED(JDielectronMcRCollDummys, "JDIELMCRCOLLDUM",
                         jdielectronmccollision::DummyDQ);

namespace jdielectronindices
{
DECLARE_SOA_INDEX_COLUMN_CUSTOM(JDielectronMcCollision, dielectronmccollision, "JDIELMCCOLLS");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, JTracks, "_0");
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, JTracks, "_1");
} // namespace jdielectronindices

DECLARE_SOA_TABLE_STAGED(JDielectronCollisionIds, "JDIELCOLLID",
                         jcandidateindices::JCollisionId,
                         o2::soa::Marker<JMarkerDielectron>);

DECLARE_SOA_TABLE_STAGED(JDielectronMcCollisionIds, "JDIELMCCOLLID",
                         jcandidateindices::JMcCollisionId,
                         o2::soa::Marker<JMarkerDielectron>);

DECLARE_SOA_TABLE_STAGED(JDielectronIds, "JDIELID",
                         jcandidateindices::JCollisionId,
                         jdielectronindices::Prong0Id,
                         jdielectronindices::Prong1Id,
                         o2::soa::Marker<JMarkerDielectron>);

namespace jdielectronmc
{
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(E, e, float);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Vx, vx, float);
DECLARE_SOA_COLUMN(Vy, vy, float);
DECLARE_SOA_COLUMN(Vz, vz, float);
DECLARE_SOA_COLUMN(Vt, vt, float);
DECLARE_SOA_COLUMN(PdgCode, pdgCode, int);
DECLARE_SOA_COLUMN(StatusCode, statusCode, int);
DECLARE_SOA_COLUMN(Flags, flags, uint8_t);
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
DECLARE_SOA_DYNAMIC_COLUMN(Energy, energy,
                           [](float e) -> float { return e; });

DECLARE_SOA_DYNAMIC_COLUMN(ProducedByGenerator, producedByGenerator, //! True if particle produced by the generator (==TMCProcess::kPrimary); False if by the transport code
                           [](uint8_t flags) -> bool { return (flags & o2::aod::mcparticle::enums::ProducedByTransport) == 0x0; });
DECLARE_SOA_DYNAMIC_COLUMN(FromBackgroundEvent, fromBackgroundEvent, //! Particle from background event
                           [](uint8_t flags) -> bool { return (flags & o2::aod::mcparticle::enums::FromBackgroundEvent) == o2::aod::mcparticle::enums::FromBackgroundEvent; });
DECLARE_SOA_DYNAMIC_COLUMN(GetProcess, getProcess, //! The VMC physics code (as int) that generated this particle (see header TMCProcess.h in ROOT)
                           [](uint8_t flags, int statusCode) -> int { if ((flags & o2::aod::mcparticle::enums::ProducedByTransport) == 0x0) { return 0 /*TMCProcess::kPrimary*/; } else { return statusCode; } });
DECLARE_SOA_DYNAMIC_COLUMN(GetGenStatusCode, getGenStatusCode, //! The native status code put by the generator, or -1 if a particle produced during transport
                           [](uint8_t flags, int statusCode) -> int { if ((flags & o2::aod::mcparticle::enums::ProducedByTransport) == 0x0) { return o2::mcgenstatus::getGenStatusCode(statusCode); } else { return -1; } });
DECLARE_SOA_DYNAMIC_COLUMN(GetHepMCStatusCode, getHepMCStatusCode, //! The HepMC status code put by the generator, or -1 if a particle produced during transport
                           [](uint8_t flags, int statusCode) -> int { if ((flags & o2::aod::mcparticle::enums::ProducedByTransport) == 0x0) { return o2::mcgenstatus::getHepMCStatusCode(statusCode); } else { return -1; } });
DECLARE_SOA_DYNAMIC_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, //! True if particle is considered a physical primary according to the ALICE definition
                           [](uint8_t flags) -> bool { return (flags & o2::aod::mcparticle::enums::PhysicalPrimary) == o2::aod::mcparticle::enums::PhysicalPrimary; });
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
                         jdielectronmc::Vx,
                         jdielectronmc::Vy,
                         jdielectronmc::Vz,
                         jdielectronmc::Vt,
                         jdielectronmc::PdgCode,
                         jdielectronmc::StatusCode,
                         jdielectronmc::Flags,
                         jdielectronmc::DecayFlag,
                         jdielectronmc::Origin,
                         jdielectronmc::Px<jdielectronmc::Pt, jdielectronmc::Phi>,
                         jdielectronmc::Py<jdielectronmc::Pt, jdielectronmc::Phi>,
                         jdielectronmc::Pz<jdielectronmc::Pt, jdielectronmc::Eta>,
                         jdielectronmc::P<jdielectronmc::Pt, jdielectronmc::Eta>,
                         jdielectronmc::Energy<jdielectronmc::E>,
                         jdielectronmc::ProducedByGenerator<jdielectronmc::Flags>,
                         jdielectronmc::FromBackgroundEvent<jdielectronmc::Flags>,
                         jdielectronmc::GetProcess<jdielectronmc::Flags, jdielectronmc::StatusCode>,
                         jdielectronmc::GetGenStatusCode<jdielectronmc::Flags, jdielectronmc::StatusCode>,
                         jdielectronmc::GetHepMCStatusCode<jdielectronmc::Flags, jdielectronmc::StatusCode>,
                         jdielectronmc::IsPhysicalPrimary<jdielectronmc::Flags>);

using JDielectronMc = JDielectronMcs::iterator;
using StoredJDielectronMc = StoredJDielectronMcs::iterator;

DECLARE_SOA_TABLE_STAGED(JDielectronMcIds, "JDIELMCID",
                         jcandidateindices::JMcCollisionId,
                         jcandidateindices::JMcParticleId,
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

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETREDUCEDDATADQ_H_
