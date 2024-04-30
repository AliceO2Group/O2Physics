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
#include <vector>
#include "Framework/AnalysisDataModel.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

namespace o2::aod
{

namespace jbc
{
DECLARE_SOA_INDEX_COLUMN(BC, bc);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(GlobalBC, globalBC, uint64_t);
DECLARE_SOA_COLUMN(Timestamp, timestamp, uint64_t);
} // namespace jbc

DECLARE_SOA_TABLE(JBCs, "AOD", "JBC",
                  o2::soa::Index<>,
                  jbc::RunNumber,
                  jbc::GlobalBC,
                  jbc::Timestamp);

using JBC = JBCs::iterator;

DECLARE_SOA_TABLE(StoredJBCs, "AOD1", "JBC",
                  o2::soa::Index<>,
                  jbc::RunNumber,
                  jbc::GlobalBC,
                  jbc::Timestamp,
                  o2::soa::Marker<1>);

using StoredJBC = StoredJBCs::iterator;

DECLARE_SOA_TABLE(JBCPIs, "AOD", "JBCPI",
                  jbc::BCId);

DECLARE_SOA_TABLE(StoredJBCPIs, "AOD1", "JBCPI",
                  jbc::BCId,
                  o2::soa::Marker<1>);

namespace jcollision
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_INDEX_COLUMN(JBC, bc);
DECLARE_SOA_COLUMN(PosX, posX, float);
DECLARE_SOA_COLUMN(PosY, posY, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float);
DECLARE_SOA_COLUMN(Centrality, centrality, float);
DECLARE_SOA_COLUMN(EventSel, eventSel, uint8_t);
DECLARE_SOA_BITMAP_COLUMN(Alias, alias, 32);
DECLARE_SOA_COLUMN(ChargedTriggerSel, chargedTriggerSel, uint8_t);
DECLARE_SOA_COLUMN(FullTriggerSel, fullTriggerSel, uint32_t);
DECLARE_SOA_COLUMN(ChargedHFTriggerSel, chargedHFTriggerSel, uint8_t);
DECLARE_SOA_COLUMN(ReadCounts, readCounts, std::vector<int>);
DECLARE_SOA_COLUMN(WrittenCounts, writtenCounts, std::vector<int>);
} // namespace jcollision

DECLARE_SOA_TABLE(JCollisions, "AOD", "JCOLLISION",
                  o2::soa::Index<>,
                  jcollision::PosX,
                  jcollision::PosY,
                  jcollision::PosZ,
                  jcollision::Multiplicity,
                  jcollision::Centrality,
                  jcollision::EventSel,
                  jcollision::Alias);

using JCollision = JCollisions::iterator;

DECLARE_SOA_TABLE(StoredJCollisions, "AOD1", "JCOLLISION",
                  o2::soa::Index<>,
                  jcollision::PosX,
                  jcollision::PosY,
                  jcollision::PosZ,
                  jcollision::Multiplicity,
                  jcollision::Centrality,
                  jcollision::EventSel,
                  jcollision::Alias,
                  o2::soa::Marker<1>);

using StoredJCollision = StoredJCollisions::iterator;

DECLARE_SOA_TABLE(JCollisionPIs, "AOD", "JCOLLISIONPI",
                  jcollision::CollisionId);

DECLARE_SOA_TABLE(StoredJCollisionPIs, "AOD1", "JCOLLISIONPI",
                  jcollision::CollisionId,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JChTrigSels, "AOD", "JCHTRIGSEL",
                  jcollision::ChargedTriggerSel);

DECLARE_SOA_TABLE(StoredJChTrigSels, "AOD1", "JCHTRIGSEL",
                  jcollision::ChargedTriggerSel,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JFullTrigSels, "AOD", "JFULLTRIGSEL",
                  jcollision::FullTriggerSel);

DECLARE_SOA_TABLE(StoredJFullTrigSels, "AOD1", "JFULLTRIGSEL",
                  jcollision::FullTriggerSel,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JChHFTrigSels, "AOD", "JCHHFTRIGSEL",
                  jcollision::ChargedHFTriggerSel);

DECLARE_SOA_TABLE(StoredJChHFTrigSels, "AOD1", "JCHHFTRIGSEL",
                  jcollision::ChargedHFTriggerSel,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JCollisionBCs, "AOD", "JCOLLISIONBC",
                  jcollision::JBCId);

DECLARE_SOA_TABLE(StoredJCollisionBCs, "AOD1", "JCOLLISIONBC",
                  jcollision::JBCId,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(CollisionCounts, "AOD", "COLLCOUNT",
                  jcollision::ReadCounts,
                  jcollision::WrittenCounts);

DECLARE_SOA_TABLE(StoredCollisionCounts, "AOD1", "COLLCOUNT",
                  jcollision::ReadCounts,
                  jcollision::WrittenCounts,
                  o2::soa::Marker<1>);

namespace jmccollision
{
DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision);
DECLARE_SOA_COLUMN(PosX, posX, float);
DECLARE_SOA_COLUMN(PosY, posY, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(Weight, weight, float);
} // namespace jmccollision
DECLARE_SOA_TABLE(JMcCollisions, "AOD", "JMCCOLLISION",
                  o2::soa::Index<>,
                  jmccollision::PosX,
                  jmccollision::PosY,
                  jmccollision::PosZ,
                  jmccollision::Weight);

using JMcCollision = JMcCollisions::iterator;

DECLARE_SOA_TABLE(StoredJMcCollisions, "AOD1", "JMCCOLLISION",
                  o2::soa::Index<>,
                  jmccollision::PosX,
                  jmccollision::PosY,
                  jmccollision::PosZ,
                  jmccollision::Weight,
                  o2::soa::Marker<1>);

using StoredJMcCollision = StoredJMcCollisions::iterator;

DECLARE_SOA_TABLE(JMcCollisionPIs, "AOD", "JMCCOLLISIONPI",
                  jmccollision::McCollisionId);

DECLARE_SOA_TABLE(StoredJMcCollisionPIs, "AOD1", "JMCCOLLISIONPI",
                  jmccollision::McCollisionId,
                  o2::soa::Marker<1>);

namespace jmccollisionlb
{
DECLARE_SOA_INDEX_COLUMN(JMcCollision, mcCollision);
}

DECLARE_SOA_TABLE(JMcCollisionLbs, "AOD", "JMCCOLLISIONLB",
                  jmccollisionlb::JMcCollisionId);

DECLARE_SOA_TABLE(StoredJMcCollisionLbs, "AOD1", "JMCCOLLISIONLB",
                  jmccollisionlb::JMcCollisionId,
                  o2::soa::Marker<1>);

namespace jtrack
{
DECLARE_SOA_INDEX_COLUMN(JCollision, collision);
DECLARE_SOA_INDEX_COLUMN(Track, track);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(DCAXY, dcaXY, float);
DECLARE_SOA_COLUMN(DCAZ, dcaZ, float);
DECLARE_SOA_COLUMN(TrackSel, trackSel, uint8_t);
DECLARE_SOA_DYNAMIC_COLUMN(Px, px,
                           [](float pt, float phi) -> float { return pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py,
                           [](float pt, float phi) -> float { return pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz,
                           [](float pt, float eta) -> float { return pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p,
                           [](float pt, float eta) -> float { return pt * std::cosh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(Energy, energy,
                           [](float pt, float eta) -> float { return std::sqrt((pt * std::cosh(eta) * pt * std::cosh(eta)) + (jetderiveddatautilities::mPion * jetderiveddatautilities::mPion)); });
DECLARE_SOA_DYNAMIC_COLUMN(Sign, sign,
                           [](uint8_t trackSel) -> int { if (trackSel & (1 << jetderiveddatautilities::JTrackSel::trackSign)){ return 1;} else{return -1;} });
} // namespace jtrack

DECLARE_SOA_TABLE(JTracks, "AOD", "JTRACK",
                  o2::soa::Index<>,
                  jtrack::JCollisionId,
                  jtrack::Pt,
                  jtrack::Eta,
                  jtrack::Phi,
                  jtrack::TrackSel,
                  jtrack::Px<jtrack::Pt, jtrack::Phi>,
                  jtrack::Py<jtrack::Pt, jtrack::Phi>,
                  jtrack::Pz<jtrack::Pt, jtrack::Eta>,
                  jtrack::P<jtrack::Pt, jtrack::Eta>,
                  jtrack::Energy<jtrack::Pt, jtrack::Eta>,
                  jtrack::Sign<jtrack::TrackSel>);

using JTrack = JTracks::iterator;

DECLARE_SOA_TABLE(StoredJTracks, "AOD1", "JTRACK",
                  o2::soa::Index<>,
                  jtrack::JCollisionId,
                  jtrack::Pt,
                  jtrack::Eta,
                  jtrack::Phi,
                  jtrack::TrackSel,
                  jtrack::Px<jtrack::Pt, jtrack::Phi>,
                  jtrack::Py<jtrack::Pt, jtrack::Phi>,
                  jtrack::Pz<jtrack::Pt, jtrack::Eta>,
                  jtrack::P<jtrack::Pt, jtrack::Eta>,
                  jtrack::Energy<jtrack::Pt, jtrack::Eta>,
                  jtrack::Sign<jtrack::TrackSel>,
                  o2::soa::Marker<1>);

using StoredJTrack = StoredJTracks::iterator;

DECLARE_SOA_TABLE(JTrackExtras, "AOD", "JTRACKEXTRA",
                  jtrack::DCAXY,
                  jtrack::DCAZ);

DECLARE_SOA_TABLE(StoredJTrackExtras, "AOD1", "JTRACKEXTRA",
                  jtrack::DCAXY,
                  jtrack::DCAZ,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JTrackPIs, "AOD", "JTRACKPI",
                  jtrack::TrackId);

DECLARE_SOA_TABLE(StoredJTrackPIs, "AOD1", "JTRACKPI",
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
DECLARE_SOA_DYNAMIC_COLUMN(Energy, energy,
                           [](float e) -> float { return e; });
} // namespace jmcparticle

DECLARE_SOA_TABLE(JMcParticles, "AOD", "JMCPARTICLE",
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
                  jmcparticle::Energy<jmcparticle::E>);

using JMcParticle = JMcParticles::iterator;

DECLARE_SOA_TABLE(StoredJMcParticles, "AOD1", "JMCPARTICLE",
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
                  jmcparticle::Energy<jmcparticle::E>,
                  o2::soa::Marker<1>);

using StoredJMcParticle = StoredJMcParticles::iterator;

DECLARE_SOA_TABLE(JMcParticlePIs, "AOD", "JMCPARTICLEPI",
                  jmcparticle::McParticleId);

DECLARE_SOA_TABLE(StoredJMcParticlePIs, "AOD1", "JMCPARTICLEPI",
                  jmcparticle::McParticleId,
                  o2::soa::Marker<1>);

namespace jmctracklb
{
DECLARE_SOA_INDEX_COLUMN(JMcParticle, mcParticle);
}

DECLARE_SOA_TABLE(JMcTrackLbs, "AOD", "JMCTRACKLB", //! Table joined to the track table containing the MC index
                  jmctracklb::JMcParticleId);

DECLARE_SOA_TABLE(StoredJMcTrackLbs, "AOD1", "JMCTRACKLB", //! Table joined to the track table containing the MC index
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

DECLARE_SOA_TABLE(JClusters, "AOD", "JCLUSTER", //!
                  o2::soa::Index<>, jcluster::JCollisionId, jcluster::ID, jcluster::Energy,
                  jcluster::CoreEnergy, jcluster::RawEnergy, jcluster::Eta, jcluster::Phi,
                  jcluster::M02, jcluster::M20, jcluster::NCells, jcluster::Time,
                  jcluster::IsExotic, jcluster::DistanceToBadChannel, jcluster::NLM, jcluster::Definition,
                  jcluster::LeadingCellEnergy, jcluster::SubleadingCellEnergy, jcluster::LeadingCellNumber, jcluster::SubleadingCellNumber);

using JCluster = JClusters::iterator;

DECLARE_SOA_TABLE(StoredJClusters, "AOD1", "JCLUSTER",
                  o2::soa::Index<>, jcluster::JCollisionId, jcluster::ID, jcluster::Energy,
                  jcluster::CoreEnergy, jcluster::RawEnergy, jcluster::Eta, jcluster::Phi,
                  jcluster::M02, jcluster::M20, jcluster::NCells, jcluster::Time,
                  jcluster::IsExotic, jcluster::DistanceToBadChannel, jcluster::NLM, jcluster::Definition,
                  jcluster::LeadingCellEnergy, jcluster::SubleadingCellEnergy, jcluster::LeadingCellNumber, jcluster::SubleadingCellNumber,
                  o2::soa::Marker<1>);

using StoredJCluster = StoredJClusters::iterator;

DECLARE_SOA_TABLE(JClusterPIs, "AOD", "JCLUSTERPI",
                  jcluster::EMCALClusterId);

DECLARE_SOA_TABLE(StoredJClusterPIs, "AOD1", "JCLUSTERPI",
                  jcluster::EMCALClusterId,
                  o2::soa::Marker<1>);

DECLARE_SOA_TABLE(JClusterTracks, "AOD", "JCLUSTERTRACK", //!
                  jcluster::JTrackIds);

DECLARE_SOA_TABLE(StoredJClusterTracks, "AOD1", "JCLUSTERTRACK", //!
                  jcluster::JTrackIds,
                  o2::soa::Marker<1>);

namespace jdummy
{

DECLARE_SOA_COLUMN(Dummy, dummy, bool);

} // namespace jdummy
DECLARE_SOA_TABLE(JDummys, "AOD", "JDUMMY",
                  o2::soa::Index<>,
                  jdummy::Dummy);

DECLARE_SOA_TABLE(StoredJDummys, "AOD1", "JDUMMY",
                  o2::soa::Index<>,
                  jdummy::Dummy,
                  o2::soa::Marker<1>);

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETREDUCEDDATA_H_
