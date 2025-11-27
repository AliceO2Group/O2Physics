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
/// \brief Table definitions for reduced data model for lf jets
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_DATAMODEL_JETREDUCEDDATAV0_H_
#define PWGJE_DATAMODEL_JETREDUCEDDATAV0_H_

#include "PWGJE/DataModel/JetReducedData.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h> // IWYU pragma: keep

#include <cmath>
#include <cstdint>

namespace o2::aod
{

DECLARE_SOA_TABLE_STAGED(JV0McCollisions, "JV0MCCOLL",
                         o2::soa::Index<>,
                         jmccollision::PosX,
                         jmccollision::PosY,
                         jmccollision::PosZ,
                         o2::soa::Marker<3>);

namespace jv0indices
{
DECLARE_SOA_INDEX_COLUMN(JCollision, collision);
DECLARE_SOA_INDEX_COLUMN(JV0McCollision, v0mccollision);
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, JTracks, "_0");
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, JTracks, "_1");
DECLARE_SOA_INDEX_COLUMN(JMcCollision, mcCollision);
DECLARE_SOA_INDEX_COLUMN(JMcParticle, mcParticle);
} // namespace jv0indices

DECLARE_SOA_TABLE_STAGED(JV0CollisionIds, "JV0COLLID",
                         jv0indices::JCollisionId);

DECLARE_SOA_TABLE_STAGED(JV0McCollisionIds, "JV0MCCOLLID",
                         jv0indices::JMcCollisionId);

DECLARE_SOA_TABLE(JV0Ids, "AOD", "JV0ID",
                  jv0indices::JCollisionId,
                  jv0indices::PosTrackId,
                  jv0indices::NegTrackId);

namespace jv0mc
{
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(E, e, float);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(PdgCode, pdgCode, int);
DECLARE_SOA_COLUMN(StatusCode, statusCode, int);
DECLARE_SOA_COLUMN(Flags, flags, uint8_t);
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(Mothers, mothers);
DECLARE_SOA_SELF_SLICE_INDEX_COLUMN(Daughters, daughters);
DECLARE_SOA_COLUMN(DecayFlag, decayFlag, int8_t);
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
} // namespace jv0mc

DECLARE_SOA_TABLE_STAGED(JV0Mcs, "JV0MC",
                         o2::soa::Index<>,
                         jv0indices::JV0McCollisionId,
                         jv0mc::Pt,
                         jv0mc::Eta,
                         jv0mc::Phi,
                         jv0mc::Y,
                         jv0mc::E,
                         jv0mc::M,
                         jv0mc::PdgCode,
                         jv0mc::StatusCode,
                         jv0mc::Flags,
                         jv0mc::DecayFlag,
                         jv0mc::Px<jv0mc::Pt, jv0mc::Phi>,
                         jv0mc::Py<jv0mc::Pt, jv0mc::Phi>,
                         jv0mc::Pz<jv0mc::Pt, jv0mc::Eta>,
                         jv0mc::P<jv0mc::Pt, jv0mc::Eta>,
                         jv0mc::Energy<jv0mc::E>,
                         jv0mc::ProducedByGenerator<jv0mc::Flags>,
                         jv0mc::FromBackgroundEvent<jv0mc::Flags>,
                         jv0mc::GetProcess<jv0mc::Flags, jv0mc::StatusCode>,
                         jv0mc::GetGenStatusCode<jv0mc::Flags, jv0mc::StatusCode>,
                         jv0mc::GetHepMCStatusCode<jv0mc::Flags, jv0mc::StatusCode>,
                         jv0mc::IsPhysicalPrimary<jv0mc::Flags>);

using JV0Mc = JV0Mcs::iterator;
using StoredJV0Mc = StoredJV0Mcs::iterator;

DECLARE_SOA_TABLE_STAGED(JV0McIds, "JV0MCID",
                         jv0indices::JMcCollisionId,
                         jv0indices::JMcParticleId,
                         jv0mc::MothersIds,
                         jv0mc::DaughtersIdSlice);

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETREDUCEDDATAV0_H_
