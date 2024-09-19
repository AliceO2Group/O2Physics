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

#include <cmath>
#include <vector>
#include "Framework/AnalysisDataModel.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/JetReducedData.h"

namespace o2::aod
{

DECLARE_SOA_TABLE(JV0McCollisions, "AOD", "JV0MCCOLL",
                  o2::soa::Index<>,
                  jmccollision::PosX,
                  jmccollision::PosY,
                  jmccollision::PosZ,
                  o2::soa::Marker<1>);

using JMcCollision = JMcCollisions::iterator;

DECLARE_SOA_TABLE(StoredJV0McCollisions, "AOD1", "JV0MCCOLL",
                  o2::soa::Index<>,
                  jmccollision::PosX,
                  jmccollision::PosY,
                  jmccollision::PosZ,
                  o2::soa::Marker<2>);

namespace jv0indices
{
DECLARE_SOA_INDEX_COLUMN(JCollision, collision);
DECLARE_SOA_INDEX_COLUMN(JV0McCollision, v0mccollision);
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, JTracks, "_0");
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, JTracks, "_1");
DECLARE_SOA_INDEX_COLUMN(JMcCollision, mcCollision);
DECLARE_SOA_INDEX_COLUMN(JMcParticle, mcParticle);
} // namespace jv0indices

DECLARE_SOA_TABLE(JV0CollisionIds, "AOD", "JV0COLLID",
                  jv0indices::JCollisionId);

DECLARE_SOA_TABLE(StoredJV0CollisionIds, "AOD1", "JV0COLLID",
                  jv0indices::JCollisionId,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JV0McCollisionIds, "AOD", "JV0MCCOLLID",
                  jv0indices::JMcCollisionId);

DECLARE_SOA_TABLE(StoredJV0McCollisionIds, "AOD1", "JV0MCCOLLID",
                  jv0indices::JMcCollisionId,
                  o2::soa::Marker<1>);

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
DECLARE_SOA_COLUMN(GenStatusCode, getGenStatusCode, int); // TODO : We can look at combining this with the two below
DECLARE_SOA_COLUMN(HepMCStatusCode, getHepMCStatusCode, int);
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool);
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
} // namespace jv0mc

DECLARE_SOA_TABLE(JV0Mcs, "AOD", "JV0MC",
                  o2::soa::Index<>,
                  jv0indices::JV0McCollisionId,
                  jv0mc::Pt,
                  jv0mc::Eta,
                  jv0mc::Phi,
                  jv0mc::Y,
                  jv0mc::E,
                  jv0mc::M,
                  jv0mc::PdgCode,
                  jv0mc::GenStatusCode,
                  jv0mc::HepMCStatusCode,
                  jv0mc::IsPhysicalPrimary,
                  jv0mc::DecayFlag,
                  jv0mc::Px<jv0mc::Pt, jv0mc::Phi>,
                  jv0mc::Py<jv0mc::Pt, jv0mc::Phi>,
                  jv0mc::Pz<jv0mc::Pt, jv0mc::Eta>,
                  jv0mc::P<jv0mc::Pt, jv0mc::Eta>);

using JV0Mc = JV0Mcs::iterator;

DECLARE_SOA_TABLE(StoredJV0Mcs, "AOD1", "JV0MC",
                  o2::soa::Index<>,
                  jv0indices::JV0McCollisionId,
                  jv0mc::Pt,
                  jv0mc::Eta,
                  jv0mc::Phi,
                  jv0mc::Y,
                  jv0mc::E,
                  jv0mc::M,
                  jv0mc::PdgCode,
                  jv0mc::GenStatusCode,
                  jv0mc::HepMCStatusCode,
                  jv0mc::IsPhysicalPrimary,
                  jv0mc::DecayFlag,
                  jv0mc::Px<jv0mc::Pt, jv0mc::Phi>,
                  jv0mc::Py<jv0mc::Pt, jv0mc::Phi>,
                  jv0mc::Pz<jv0mc::Pt, jv0mc::Eta>,
                  jv0mc::P<jv0mc::Pt, jv0mc::Eta>,
                  o2::soa::Marker<1>);

using StoredJV0Mc = StoredJV0Mcs::iterator;

DECLARE_SOA_TABLE(JV0McIds, "AOD", "JV0MCID",
                  jv0indices::JMcCollisionId,
                  jv0indices::JMcParticleId,
                  jv0mc::MothersIds,
                  jv0mc::DaughtersIdSlice);

DECLARE_SOA_TABLE(StoredJV0McIds, "AOD1", "JV0MCID",
                  jv0indices::JMcCollisionId,
                  jv0indices::JMcParticleId,
                  jv0mc::MothersIds,
                  jv0mc::DaughtersIdSlice,
                  o2::soa::Marker<1>);

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETREDUCEDDATAV0_H_
