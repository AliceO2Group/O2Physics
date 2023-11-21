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
#include "PWGJE/DataModel/EMCALClusters.h"

namespace o2::aod
{

namespace jbc
{
DECLARE_SOA_INDEX_COLUMN(BC, bc);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(GlobalBC, globalBC, uint64_t);
DECLARE_SOA_COLUMN(Timestamp, timestamp, uint64_t);
} // namespace jbc

DECLARE_SOA_TABLE(JBCs, "AOD", "JBCs",
                  o2::soa::Index<>,
                  jbc::RunNumber,
                  jbc::GlobalBC,
                  jbc::Timestamp);

using JBC = JBCs::iterator;

DECLARE_SOA_TABLE(StoredJBCs, "DYN", "JBCs",
                  o2::soa::Index<>,
                  jbc::RunNumber,
                  jbc::GlobalBC,
                  jbc::Timestamp,
                  o2::soa::Marker<1>);

using StoredJBC = StoredJBCs::iterator;

DECLARE_SOA_TABLE(JBCPIs, "AOD", "JBCPIs",
                  jbc::BCId);

DECLARE_SOA_TABLE(StoredJBCPIs, "DYN", "JBCPIs",
                  jbc::BCId,
                  o2::soa::Marker<1>);

namespace jcollision
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_INDEX_COLUMN(JBC, bc);
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

DECLARE_SOA_TABLE(StoredJCollisions, "DYN", "JCollisions",
                  o2::soa::Index<>,
                  jcollision::PosZ,
                  jcollision::EventSel,
                  jcollision::Alias,
                  o2::soa::Marker<1>);

using StoredJCollision = StoredJCollisions::iterator;

DECLARE_SOA_TABLE(JCollisionPIs, "AOD", "JCollisionPIs",
                  jcollision::CollisionId);

DECLARE_SOA_TABLE(StoredJCollisionPIs, "DYN", "JCollisionPIs",
                  jcollision::CollisionId,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JChTrigSels, "AOD", "JChrgTrigSels",
                  jcollision::ChargedTriggerSel);

DECLARE_SOA_TABLE(StoredJChTrigSels, "DYN", "JChargTrigSels",
                  jcollision::ChargedTriggerSel,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JFullTrigSels, "AOD", "JFullTrigSels",
                  jcollision::FullTriggerSel);

DECLARE_SOA_TABLE(StoredJFullTrigSels, "DYN", "JFullTrigSels",
                  jcollision::FullTriggerSel,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JCollisionBCs, "AOD", "JCollisionBcs",
                  jcollision::JBCId);

DECLARE_SOA_TABLE(StoredJCollisionBCs, "DYN", "JCollisionBcs",
                  jcollision::JBCId,
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

DECLARE_SOA_TABLE(StoredJMcCollisions, "DYN", "JMcCollisions",
                  o2::soa::Index<>,
                  jmccollision::PosZ,
                  jmccollision::Weight,
                  o2::soa::Marker<1>);

using StoredJMcCollision = StoredJMcCollisions::iterator;

DECLARE_SOA_TABLE(JMcCollisionPIs, "AOD", "JMcCollisionPIs",
                  jmccollision::McCollisionId);

DECLARE_SOA_TABLE(StoredJMcCollisionPIs, "DYN", "JMcCollisionPIs",
                  jmccollision::McCollisionId,
                  o2::soa::Marker<1>);

namespace jmccollisionlb
{
DECLARE_SOA_INDEX_COLUMN(JMcCollision, mcCollision);
}

DECLARE_SOA_TABLE(JMcCollisionLbs, "AOD", "JMcCollisionLbs",
                  jmccollisionlb::JMcCollisionId);

DECLARE_SOA_TABLE(StoredJMcCollisionLbs, "DYN", "JMcCollisionLbs",
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
DECLARE_SOA_COLUMN(Sign, sign, float);
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
                  jtrack::Sign,
                  jtrack::TrackSel,
                  jtrack::Px<jtrack::Pt, jtrack::Phi>,
                  jtrack::Py<jtrack::Pt, jtrack::Phi>,
                  jtrack::Pz<jtrack::Pt, jtrack::Eta>,
                  jtrack::P<jtrack::Pt, jtrack::Eta>);

using JTrack = JTracks::iterator;

DECLARE_SOA_TABLE(StoredJTracks, "DYN", "JTracks",
                  o2::soa::Index<>,
                  jtrack::JCollisionId,
                  jtrack::Pt,
                  jtrack::Eta,
                  jtrack::Phi,
                  jtrack::Energy,
                  jtrack::Sign,
                  jtrack::TrackSel,
                  jtrack::Px<jtrack::Pt, jtrack::Phi>,
                  jtrack::Py<jtrack::Pt, jtrack::Phi>,
                  jtrack::Pz<jtrack::Pt, jtrack::Eta>,
                  jtrack::P<jtrack::Pt, jtrack::Eta>,
                  o2::soa::Marker<1>);

using StoredJTrack = StoredJTracks::iterator;

DECLARE_SOA_TABLE(JTrackPIs, "AOD", "JTrackPIs",
                  jtrack::TrackId);

DECLARE_SOA_TABLE(StoredJTrackPIs, "DYN", "JTrackPIs",
                  jtrack::TrackId,
                  o2::soa::Marker<1>);

namespace jmcparticle
{
DECLARE_SOA_INDEX_COLUMN(JMcCollision, mcCollision);
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(E, e, float);
DECLARE_SOA_COLUMN(PdgCode, pdgCode, int);
DECLARE_SOA_COLUMN(GenStatusCode, getGenStatusCode, int); // TODO : We can look at combining this with the two below
DECLARE_SOA_COLUMN(HepMCStatusCode, getHepMCStatusCode, int);
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isPhysicalPrimary, bool);
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
                  jmcparticle::E,
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

DECLARE_SOA_TABLE(StoredJMcParticles, "DYN", "JMcParticles",
                  o2::soa::Index<>,
                  jmcparticle::JMcCollisionId,
                  jmcparticle::Pt,
                  jmcparticle::Eta,
                  jmcparticle::Phi,
                  jmcparticle::Y,
                  jmcparticle::E,
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

DECLARE_SOA_TABLE(StoredJMcParticlePIs, "DYN", "JMcParticlePIs",
                  jmcparticle::McParticleId,
                  o2::soa::Marker<1>);

namespace jmctracklb
{
DECLARE_SOA_INDEX_COLUMN(JMcParticle, mcParticle);
}

DECLARE_SOA_TABLE(JMcTrackLbs, "AOD", "JMcTrackLbs", //! Table joined to the track table containing the MC index
                  jmctracklb::JMcParticleId);

DECLARE_SOA_TABLE(StoredJMcTrackLbs, "DYN", "JMcTrackLbs", //! Table joined to the track table containing the MC index
                  jmctracklb::JMcParticleId,
                  o2::soa::Marker<1>);

namespace jcluster
{
DECLARE_SOA_INDEX_COLUMN(JCollision, collision);                       //! collisionID used as index for matched clusters
DECLARE_SOA_INDEX_COLUMN(EMCALCluster, cluster);                       //! cluster ID of original cluster
DECLARE_SOA_ARRAY_INDEX_COLUMN(JTrack, matchedTracks);                 // ! array of indices to tracks matched to the cluster
DECLARE_SOA_COLUMN(ID, id, int);                                       //! cluster ID identifying cluster in event
DECLARE_SOA_COLUMN(Energy, energy, float);                             //! cluster energy (GeV)
DECLARE_SOA_COLUMN(CoreEnergy, coreEnergy, float);                     //! cluster core energy (GeV)
DECLARE_SOA_COLUMN(RawEnergy, rawEnergy, float);                       //! raw cluster energy (GeV)
DECLARE_SOA_COLUMN(Eta, eta, float);                                   //! cluster pseudorapidity (calculated using vertex)
DECLARE_SOA_COLUMN(Phi, phi, float);                                   //! cluster azimuthal angle (calculated using vertex)
DECLARE_SOA_COLUMN(M02, m02, float);                                   //! shower shape long axis
DECLARE_SOA_COLUMN(M20, m20, float);                                   //! shower shape short axis
DECLARE_SOA_COLUMN(NCells, nCells, int);                               //! number of cells in cluster
DECLARE_SOA_COLUMN(Time, time, float);                                 //! cluster time (ns)
DECLARE_SOA_COLUMN(IsExotic, isExotic, bool);                          //! flag to mark cluster as exotic
DECLARE_SOA_COLUMN(DistanceToBadChannel, distanceToBadChannel, float); //! distance to bad channel
DECLARE_SOA_COLUMN(NLM, nlm, int);                                     //! number of local maxima
DECLARE_SOA_COLUMN(Definition, definition, int);                       //! cluster definition, see EMCALClusterDefinition.h
DECLARE_SOA_COLUMN(LeadingCellEnergy, leadingCellEnergy, float);       //! energy of leading cell in the cluster
DECLARE_SOA_COLUMN(SubleadingCellEnergy, subleadingCellEnergy, float); //! energy of leading cell in the cluster
DECLARE_SOA_COLUMN(LeadingCellNumber, leadingCellNumber, int);         //! energy of leading cell in the cluster
DECLARE_SOA_COLUMN(SubleadingCellNumber, subleadingCellNumber, int);   //! energy of leading cell in the cluster

} // namespace jcluster

DECLARE_SOA_TABLE(JClusters, "AOD", "JClusters", //!
                  o2::soa::Index<>, jcluster::JCollisionId, jcluster::ID, jcluster::Energy,
                  jcluster::CoreEnergy, jcluster::RawEnergy, jcluster::Eta, jcluster::Phi,
                  jcluster::M02, jcluster::M20, jcluster::NCells, jcluster::Time,
                  jcluster::IsExotic, jcluster::DistanceToBadChannel, jcluster::NLM, jcluster::Definition,
                  jcluster::LeadingCellEnergy, jcluster::SubleadingCellEnergy, jcluster::LeadingCellNumber, jcluster::SubleadingCellNumber);

using JCluster = JClusters::iterator;

DECLARE_SOA_TABLE(StoredJClusters, "DYN", "JClusters",
                  o2::soa::Index<>, jcluster::JCollisionId, jcluster::ID, jcluster::Energy,
                  jcluster::CoreEnergy, jcluster::RawEnergy, jcluster::Eta, jcluster::Phi,
                  jcluster::M02, jcluster::M20, jcluster::NCells, jcluster::Time,
                  jcluster::IsExotic, jcluster::DistanceToBadChannel, jcluster::NLM, jcluster::Definition,
                  jcluster::LeadingCellEnergy, jcluster::SubleadingCellEnergy, jcluster::LeadingCellNumber, jcluster::SubleadingCellNumber,
                  o2::soa::Marker<1>);

using StoredJCluster = StoredJClusters::iterator;

DECLARE_SOA_TABLE(JClusterPIs, "AOD", "JClusterPIs",
                  jcluster::EMCALClusterId);

DECLARE_SOA_TABLE(StoredJClusterPIs, "DYN", "JClusterPIs",
                  jcluster::EMCALClusterId,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JClusterTracks, "AOD", "JClusterTracks", //!
                  jcluster::JTrackIds);

DECLARE_SOA_TABLE(StoredJClusterTracks, "DYN", "JClusterTracks", //!
                  jcluster::JTrackIds,
                  o2::soa::Marker<1>);

namespace jdummy
{

DECLARE_SOA_COLUMN(Dummy, dummy, bool);

} // namespace jdummy
DECLARE_SOA_TABLE(JDummys, "AOD", "JDummys",
                  o2::soa::Index<>,
                  jdummy::Dummy);

DECLARE_SOA_TABLE(StoredJDummys, "DYN", "JDummys",
                  o2::soa::Index<>,
                  jdummy::Dummy,
                  o2::soa::Marker<1>);

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETREDUCEDDATA_H_
