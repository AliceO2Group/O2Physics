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

#ifndef PWGJE_DATAMODEL_JETREDUCEDDATALF_H_
#define PWGJE_DATAMODEL_JETREDUCEDDATALF_H_

#include <cmath>
#include "Framework/AnalysisDataModel.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/JetReducedData.h"

namespace o2::aod
{

namespace jv0indices
{
DECLARE_SOA_INDEX_COLUMN(JCollision, collision);
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, JTracks, "_0");
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, JTracks, "_1");
DECLARE_SOA_INDEX_COLUMN(JMcCollision, mcCollision);
DECLARE_SOA_INDEX_COLUMN(JMcParticle, mcParticle);
} // namespace jv0indices

DECLARE_SOA_TABLE(JV0Ids, "AOD", "JV0ID",
                  jv0indices::JCollisionId,
                  jv0indices::PosTrackId,
                  jv0indices::NegTrackId);

DECLARE_SOA_TABLE(StoredJV0Ids, "AOD1", "JV0ID",
                  jv0indices::JCollisionId,
                  jv0indices::PosTrackId,
                  jv0indices::NegTrackId,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JV0PIds, "AOD", "JV0PID",
                  jv0indices::JMcCollisionId,
                  jv0indices::JMcParticleId);

DECLARE_SOA_TABLE(StoredJV0PIds, "AOD1", "JV0PID",
                  jv0indices::JMcCollisionId,
                  jv0indices::JMcParticleId,
                  o2::soa::Marker<1>);

namespace jv0mcparticle
{
DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollisionParent);
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticleParent);
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
} // namespace jv0mcparticle

DECLARE_SOA_TABLE(JV0McParticles, "AOD", "JV0MCPARTICLE",
                  o2::soa::Index<>,
                  jv0mcparticle::McCollisionId,
                  jv0mcparticle::McParticleId,
                  jv0mcparticle::Pt,
                  jv0mcparticle::Eta,
                  jv0mcparticle::Phi,
                  jv0mcparticle::Y,
                  jv0mcparticle::E,
                  jv0mcparticle::M,
                  jv0mcparticle::PdgCode,
                  jv0mcparticle::GenStatusCode,
                  jv0mcparticle::HepMCStatusCode,
                  jv0mcparticle::IsPhysicalPrimary,
                  jv0mcparticle::MothersIds,
                  jv0mcparticle::DaughtersIdSlice,
                  jv0mcparticle::DecayFlag,
                  jv0mcparticle::Px<jv0mcparticle::Pt, jv0mcparticle::Phi>,
                  jv0mcparticle::Py<jv0mcparticle::Pt, jv0mcparticle::Phi>,
                  jv0mcparticle::Pz<jv0mcparticle::Pt, jv0mcparticle::Eta>,
                  jv0mcparticle::P<jv0mcparticle::Pt, jv0mcparticle::Eta>);

using JV0McParticle = JV0McParticles::iterator;

DECLARE_SOA_TABLE(StoredJV0McParticles, "AOD1", "JV0MCPARTICLE",
                  o2::soa::Index<>,
                  jv0mcparticle::Pt,
                  jv0mcparticle::Eta,
                  jv0mcparticle::Phi,
                  jv0mcparticle::Y,
                  jv0mcparticle::E,
                  jv0mcparticle::M,
                  jv0mcparticle::PdgCode,
                  jv0mcparticle::GenStatusCode,
                  jv0mcparticle::HepMCStatusCode,
                  jv0mcparticle::IsPhysicalPrimary,
                  jv0mcparticle::MothersIds,
                  jv0mcparticle::DaughtersIdSlice,
                  jv0mcparticle::DecayFlag,
                  jv0mcparticle::Px<jv0mcparticle::Pt, jv0mcparticle::Phi>,
                  jv0mcparticle::Py<jv0mcparticle::Pt, jv0mcparticle::Phi>,
                  jv0mcparticle::Pz<jv0mcparticle::Pt, jv0mcparticle::Eta>,
                  jv0mcparticle::P<jv0mcparticle::Pt, jv0mcparticle::Eta>,
                  o2::soa::Marker<1>);

using StoredJV0McParticle = StoredJV0McParticles::iterator;

DECLARE_SOA_TABLE(JV0McPIs, "AOD", "JV0MCPI",
                  jv0mcparticle::McCollisionId, jv0mcparticle::McParticleId);

DECLARE_SOA_TABLE(StoredJV0McPIs, "AOD1", "JV0MCPI",
                  jv0mcparticle::McCollisionId, jv0mcparticle::McParticleId,
                  o2::soa::Marker<1>);

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETREDUCEDDATALF_H_
