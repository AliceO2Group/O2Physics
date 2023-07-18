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
/// \brief Table definitions for reduced data model for jets
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_DATAMODEL_JETREDUCEDDATA_H_
#define PWGJE_DATAMODEL_JETREDUCEDDATA_H_

#include <cmath>
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{

namespace jcollision
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(EventSel, eventSel, uint8_t);
DECLARE_SOA_BITMAP_COLUMN(Alias, alias, 32);
DECLARE_SOA_COLUMN(ChargedTriggerSel, chargedTriggerSel, uint16_t);
DECLARE_SOA_COLUMN(FullTriggerSel, fullTriggerSel, uint16_t);
} // namespace jcollision

DECLARE_SOA_TABLE(JCollisions, "AOD", "JCollisions",
                  o2::soa::Index<>,
                  jcollision::PosZ,
                  jcollision::EventSel,
                  jcollision::Alias);

using JCollision = JCollisions::iterator;

DECLARE_SOA_TABLE(StoredJCollisions, "AOD", "JCollisions",
                  o2::soa::Index<>,
                  jcollision::PosZ,
                  jcollision::EventSel,
                  jcollision::Alias,
                  o2::soa::Marker<1>);

using StoredJCollision = StoredJCollisions::iterator;

DECLARE_SOA_TABLE(JCollisionPIs, "AOD", "JCollisionPIs",
                  jcollision::CollisionId);

DECLARE_SOA_TABLE(StoredJCollisionPIs, "AOD", "JCollisionPIs",
                  jcollision::CollisionId,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JChTrigSels, "AOD", "JChrgTrigSels",
                  jcollision::ChargedTriggerSel);

DECLARE_SOA_TABLE(StoredJChTrigSels, "AOD", "JChargTrigSels",
                  jcollision::ChargedTriggerSel,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JFullTrigSels, "AOD", "JFullTrigSels",
                  jcollision::FullTriggerSel);

DECLARE_SOA_TABLE(StoredJFullTrigSels, "AOD", "JFullTrigSels",
                  jcollision::FullTriggerSel,
                  o2::soa::Marker<1>);

namespace jmccollision
{
DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(Weight, weight, float);
} // namespace jmccollision
DECLARE_SOA_TABLE(JMcCollisions, "AOD", "JMcCollisions",
                  o2::soa::Index<>,
                  jmccollision::PosZ,
                  jmccollision::Weight);

using JMcCollision = JMcCollisions::iterator;

DECLARE_SOA_TABLE(StoredJMcCollisions, "AOD", "JMcCollisions",
                  o2::soa::Index<>,
                  jmccollision::PosZ,
                  jmccollision::Weight,
                  o2::soa::Marker<1>);

using StoredJMcCollision = StoredJMcCollisions::iterator;

DECLARE_SOA_TABLE(JMcCollisionPIs, "AOD", "JMcCollisionPIs",
                  jmccollision::McCollisionId);

DECLARE_SOA_TABLE(StoredJMcCollisionPIs, "AOD", "JMcCollisionPIs",
                  jmccollision::McCollisionId,
                  o2::soa::Marker<1>);

namespace jmccollisionlb
{
DECLARE_SOA_INDEX_COLUMN(JMcCollision, mcCollision);
}

DECLARE_SOA_TABLE(JMcCollisionLbs, "AOD", "JMcCollisionLbs",
                  jmccollisionlb::JMcCollisionId);

DECLARE_SOA_TABLE(StoredJMcCollisionLbs, "AOD", "JMcCollisionLbs",
                  jmccollisionlb::JMcCollisionId,
                  o2::soa::Marker<1>);

namespace jtrack
{
DECLARE_SOA_INDEX_COLUMN(JCollision, collision);
DECLARE_SOA_INDEX_COLUMN(Track, track);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Energy, energy, float);
DECLARE_SOA_COLUMN(TrackSel, trackSel, uint8_t);
DECLARE_SOA_DYNAMIC_COLUMN(Px, px,
                           [](float pt, float phi) -> float { return pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py,
                           [](float pt, float phi) -> float { return pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz,
                           [](float pt, float eta) -> float { return pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p,
                           [](float pt, float eta) -> float { return pt * std::cosh(eta); });
} // namespace jtrack

DECLARE_SOA_TABLE(JTracks, "AOD", "JTracks",
                  o2::soa::Index<>,
                  jtrack::JCollisionId,
                  jtrack::Pt,
                  jtrack::Eta,
                  jtrack::Phi,
                  jtrack::Energy,
                  jtrack::TrackSel,
                  jtrack::Px<jtrack::Pt, jtrack::Phi>,
                  jtrack::Py<jtrack::Pt, jtrack::Phi>,
                  jtrack::Pz<jtrack::Pt, jtrack::Eta>,
                  jtrack::P<jtrack::Pt, jtrack::Eta>);

using JTrack = JTracks::iterator;

DECLARE_SOA_TABLE(StoredJTracks, "AOD", "JTracks",
                  o2::soa::Index<>,
                  jtrack::JCollisionId,
                  jtrack::Pt,
                  jtrack::Eta,
                  jtrack::Phi,
                  jtrack::Energy,
                  jtrack::TrackSel,
                  jtrack::Px<jtrack::Pt, jtrack::Phi>,
                  jtrack::Py<jtrack::Pt, jtrack::Phi>,
                  jtrack::Pz<jtrack::Pt, jtrack::Eta>,
                  jtrack::P<jtrack::Pt, jtrack::Eta>,
                  o2::soa::Marker<1>);

using StoredJTrack = StoredJTracks::iterator;

DECLARE_SOA_TABLE(JTrackPIs, "AOD", "JTrackPIs",
                  jtrack::TrackId);

DECLARE_SOA_TABLE(StoredJTrackPIs, "AOD", "JTrackPIs",
                  jtrack::TrackId,
                  o2::soa::Marker<1>);

namespace jmcparticle
{
DECLARE_SOA_INDEX_COLUMN(JMcCollision, collision);
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(Energy, energy, float);
DECLARE_SOA_COLUMN(PdgCode, pdgCode, int);
DECLARE_SOA_COLUMN(GenStatusCode, getGenStatusCode, int); //TODO : We can look at combining this with the two below
DECLARE_SOA_COLUMN(HepMCStatusCode, getHepMCStatusCode, int);
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, uint8_t);
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(Mothers, mothers);
DECLARE_SOA_SELF_SLICE_INDEX_COLUMN(Daughters, daughters);
DECLARE_SOA_DYNAMIC_COLUMN(Px, px,
                           [](float pt, float phi) -> float { return pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py,
                           [](float pt, float phi) -> float { return pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz,
                           [](float pt, float eta) -> float { return pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p,
                           [](float pt, float eta) -> float { return pt * std::cosh(eta); });
} // namespace jmcparticle

DECLARE_SOA_TABLE(JMcParticles, "AOD", "JMcParticles",
                  o2::soa::Index<>,
                  jmcparticle::JMcCollisionId,
                  jmcparticle::Pt,
                  jmcparticle::Eta,
                  jmcparticle::Phi,
                  jmcparticle::Y,
                  jmcparticle::Energy,
                  jmcparticle::PdgCode,
                  jmcparticle::GenStatusCode,
                  jmcparticle::HepMCStatusCode,
                  jmcparticle::IsPhysicalPrimary,
                  jmcparticle::MothersIds,
                  jmcparticle::DaughtersIdSlice,
                  jmcparticle::Px<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Py<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Pz<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::P<jmcparticle::Pt, jmcparticle::Eta>);

using JMcParticle = JMcParticles::iterator;

DECLARE_SOA_TABLE(StoredJMcParticles, "AOD", "JMcParticles",
                  o2::soa::Index<>,
                  jmcparticle::JMcCollisionId,
                  jmcparticle::Pt,
                  jmcparticle::Eta,
                  jmcparticle::Phi,
                  jmcparticle::Y,
                  jmcparticle::Energy,
                  jmcparticle::PdgCode,
                  jmcparticle::GenStatusCode,
                  jmcparticle::HepMCStatusCode,
                  jmcparticle::IsPhysicalPrimary,
                  jmcparticle::MothersIds,
                  jmcparticle::DaughtersIdSlice,
                  jmcparticle::Px<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Py<jmcparticle::Pt, jmcparticle::Phi>,
                  jmcparticle::Pz<jmcparticle::Pt, jmcparticle::Eta>,
                  jmcparticle::P<jmcparticle::Pt, jmcparticle::Eta>,
                  o2::soa::Marker<1>);

using StoredJMcParticle = StoredJMcParticles::iterator;

DECLARE_SOA_TABLE(JMcParticlePIs, "AOD", "JMcParticlePIs",
                  jmcparticle::McParticleId);

DECLARE_SOA_TABLE(StoredJMcParticlePIs, "AOD", "JMcParticlePIs",
                  jmcparticle::McParticleId,
                  o2::soa::Marker<1>);

namespace jmctracklb
{
DECLARE_SOA_INDEX_COLUMN(JMcParticle, mcParticle);
}

DECLARE_SOA_TABLE(JMcTrackLbs, "AOD", "JMcTrackLbs", //! Table joined to the track table containing the MC index
                  jmctracklb::JMcParticleId);

DECLARE_SOA_TABLE(StoredJMcTrackLbs, "AOD", "JMcTrackLbs", //! Table joined to the track table containing the MC index
                  jmctracklb::JMcParticleId,
                  o2::soa::Marker<1>);

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETREDUCEDDATA_H_