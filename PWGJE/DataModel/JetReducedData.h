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
DECLARE_SOA_BITMAP_COLUMN(Alias, alias, 32);
DECLARE_SOA_BITMAP_COLUMN(Selection, selection, 64);
DECLARE_SOA_COLUMN(ReadCounts, readCounts, std::vector<int>);
DECLARE_SOA_COLUMN(ReadCountsWithTVX, readCountsWithTVX, std::vector<int>);
DECLARE_SOA_COLUMN(ReadCountsWithTVXAndNoTFB, readCountsWithTVXAndNoTFB, std::vector<int>);
DECLARE_SOA_COLUMN(ReadCountsWithTVXAndNoTFBAndNoITSROFB, readCountsWithTVXAndNoTFBAndNoITSROFB, std::vector<int>);
} // namespace jbc

DECLARE_SOA_TABLE_STAGED(JBCs, "JBC",
                         o2::soa::Index<>,
                         jbc::RunNumber,
                         jbc::GlobalBC,
                         jbc::Timestamp,
                         jbc::Alias,
                         jbc::Selection);

using JBC = JBCs::iterator;
using StoredJBC = StoredJBCs::iterator;

DECLARE_SOA_TABLE_STAGED(JBCPIs, "JBCPI",
                         jbc::BCId);

DECLARE_SOA_TABLE_STAGED(BCCounts, "BCCOUNT",
                         jbc::ReadCounts,
                         jbc::ReadCountsWithTVX,
                         jbc::ReadCountsWithTVXAndNoTFB,
                         jbc::ReadCountsWithTVXAndNoTFBAndNoITSROFB);

namespace jcollision
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_INDEX_COLUMN(JBC, bc);
DECLARE_SOA_COLUMN(PosX, posX, float);
DECLARE_SOA_COLUMN(PosY, posY, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float);
DECLARE_SOA_COLUMN(Centrality, centrality, float);
DECLARE_SOA_COLUMN(Weight, weight, float);
DECLARE_SOA_COLUMN(SubGeneratorId, subGeneratorId, int);
DECLARE_SOA_COLUMN(EventSel, eventSel, uint16_t);
DECLARE_SOA_BITMAP_COLUMN(Alias, alias, 32);
DECLARE_SOA_COLUMN(TrackOccupancyInTimeRange, trackOccupancyInTimeRange, int);
DECLARE_SOA_COLUMN(TriggerSel, triggerSel, uint64_t);
DECLARE_SOA_COLUMN(ChargedTriggerSel, chargedTriggerSel, uint8_t);
DECLARE_SOA_COLUMN(FullTriggerSel, fullTriggerSel, uint32_t);
DECLARE_SOA_COLUMN(ChargedHFTriggerSel, chargedHFTriggerSel, uint8_t);
DECLARE_SOA_COLUMN(ReadCounts, readCounts, std::vector<int>);
DECLARE_SOA_COLUMN(ReadCountsWithTVX, readCountsWithTVX, std::vector<int>);
DECLARE_SOA_COLUMN(ReadCountsWithTVXAndZVertexAndSel8, readCountsWithTVXAndZVertexAndSel8, std::vector<int>);
DECLARE_SOA_COLUMN(ReadCountsWithTVXAndZVertexAndSel8Full, readCountsWithTVXAndZVertexAndSel8Full, std::vector<int>);
DECLARE_SOA_COLUMN(ReadCountsWithTVXAndZVertexAndSel8FullPbPb, readCountsWithTVXAndZVertexAndSel8FullPbPb, std::vector<int>);
DECLARE_SOA_COLUMN(ReadCountsWithTVXAndZVertexAndSelMC, readCountsWithTVXAndZVertexAndSelMC, std::vector<int>);
DECLARE_SOA_COLUMN(ReadCountsWithTVXAndZVertexAndSelMCFull, readCountsWithTVXAndZVertexAndSelMCFull, std::vector<int>);
DECLARE_SOA_COLUMN(ReadCountsWithTVXAndZVertexAndSelMCFullPbPb, readCountsWithTVXAndZVertexAndSelMCFullPbPb, std::vector<int>);
DECLARE_SOA_COLUMN(ReadCountsWithTVXAndZVertexAndSelUnanchoredMC, readCountsWithTVXAndZVertexAndSelUnanchoredMC, std::vector<int>);
DECLARE_SOA_COLUMN(ReadCountsWithTVXAndZVertexAndSelTVX, readCountsWithTVXAndZVertexAndSelTVX, std::vector<int>);
DECLARE_SOA_COLUMN(ReadCountsWithTVXAndZVertexAndSel7, readCountsWithTVXAndZVertexAndSel7, std::vector<int>);
DECLARE_SOA_COLUMN(ReadCountsWithTVXAndZVertexAndSel7KINT7, readCountsWithTVXAndZVertexAndSel7KINT7, std::vector<int>);
DECLARE_SOA_COLUMN(ReadCountsWithCustom, readCountsWithCustom, std::vector<int>);
DECLARE_SOA_COLUMN(IsAmbiguous, isAmbiguous, bool);
DECLARE_SOA_COLUMN(IsEMCALReadout, isEmcalReadout, bool);
} // namespace jcollision

DECLARE_SOA_TABLE_STAGED(JCollisions, "JCOLLISION",
                         o2::soa::Index<>,
                         jcollision::PosX,
                         jcollision::PosY,
                         jcollision::PosZ,
                         jcollision::Multiplicity,
                         jcollision::Centrality,
                         jcollision::TrackOccupancyInTimeRange,
                         jcollision::EventSel,
                         jcollision::Alias,
                         jcollision::TriggerSel);

using JCollision = JCollisions::iterator;
using StoredJCollision = StoredJCollisions::iterator;

DECLARE_SOA_TABLE_STAGED(JCollisionMcInfos, "JCOLLISIONMCINFO",
                         jcollision::Weight,
                         jcollision::SubGeneratorId);

DECLARE_SOA_TABLE_STAGED(JEMCCollisionLbs, "JEMCCOLLISIONLB",
                         jcollision::IsAmbiguous,
                         jcollision::IsEMCALReadout);
using JEMCCollisionLb = JEMCCollisionLbs::iterator;
using StoredJEMCCollisionLb = StoredJEMCCollisionLbs::iterator;

DECLARE_SOA_TABLE_STAGED(JCollisionPIs, "JCOLLISIONPI",
                         jcollision::CollisionId);

DECLARE_SOA_TABLE_STAGED(JCollisionBCs, "JCOLLISIONBC",
                         jcollision::JBCId);

DECLARE_SOA_TABLE(JChTrigSels, "AOD", "JCHTRIGSEL",
                  jcollision::ChargedTriggerSel);

DECLARE_SOA_TABLE(JFullTrigSels, "AOD", "JFULLTRIGSEL",
                  jcollision::FullTriggerSel);

DECLARE_SOA_TABLE(JChHFTrigSels, "AOD", "JCHHFTRIGSEL",
                  jcollision::ChargedHFTriggerSel);

DECLARE_SOA_TABLE_STAGED(CollisionCounts, "COLLCOUNT",
                         jcollision::ReadCounts,
                         jcollision::ReadCountsWithTVX,
                         jcollision::ReadCountsWithTVXAndZVertexAndSel8,
                         jcollision::ReadCountsWithTVXAndZVertexAndSel8Full,
                         jcollision::ReadCountsWithTVXAndZVertexAndSel8FullPbPb,
                         jcollision::ReadCountsWithTVXAndZVertexAndSelMC,
                         jcollision::ReadCountsWithTVXAndZVertexAndSelMCFull,
                         jcollision::ReadCountsWithTVXAndZVertexAndSelMCFullPbPb,
                         jcollision::ReadCountsWithTVXAndZVertexAndSelUnanchoredMC,
                         jcollision::ReadCountsWithTVXAndZVertexAndSelTVX,
                         jcollision::ReadCountsWithTVXAndZVertexAndSel7,
                         jcollision::ReadCountsWithTVXAndZVertexAndSel7KINT7,
                         jcollision::ReadCountsWithCustom);

namespace jmccollision
{
DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision);
DECLARE_SOA_COLUMN(PosX, posX, float);
DECLARE_SOA_COLUMN(PosY, posY, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(Weight, weight, float);
DECLARE_SOA_COLUMN(SubGeneratorId, subGeneratorId, int);
} // namespace jmccollision
DECLARE_SOA_TABLE_STAGED(JMcCollisions, "JMCCOLLISION",
                         o2::soa::Index<>,
                         jmccollision::PosX,
                         jmccollision::PosY,
                         jmccollision::PosZ,
                         jmccollision::Weight,
                         jmccollision::SubGeneratorId);

using JMcCollision = JMcCollisions::iterator;
using StoredJMcCollision = StoredJMcCollisions::iterator;

DECLARE_SOA_TABLE_STAGED(JMcCollisionPIs, "JMCCOLLISIONPI",
                         jmccollision::McCollisionId);

namespace jmccollisionlb
{
DECLARE_SOA_INDEX_COLUMN(JMcCollision, mcCollision);
}

DECLARE_SOA_TABLE_STAGED(JMcCollisionLbs, "JMCCOLLISIONLB",
                         jmccollisionlb::JMcCollisionId);

namespace jtrack
{
DECLARE_SOA_INDEX_COLUMN(JCollision, collision);
DECLARE_SOA_INDEX_COLUMN(Track, track);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(DCAX, dcaX, float);
DECLARE_SOA_COLUMN(DCAY, dcaY, float);
DECLARE_SOA_COLUMN(DCAZ, dcaZ, float);
DECLARE_SOA_COLUMN(DCAXY, dcaXY, float);
DECLARE_SOA_COLUMN(DCAXYZ, dcaXYZ, float);
DECLARE_SOA_COLUMN(SigmaDCAZ, sigmadcaZ, float);
DECLARE_SOA_COLUMN(SigmaDCAXY, sigmadcaXY, float);
DECLARE_SOA_COLUMN(SigmaDCAXYZ, sigmadcaXYZ, float);
DECLARE_SOA_COLUMN(Sigma1Pt, sigma1Pt, float);
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

DECLARE_SOA_TABLE_STAGED(JTracks, "JTRACK",
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
using StoredJTrack = StoredJTracks::iterator;

DECLARE_SOA_TABLE_STAGED(JTrackExtras, "JTRACKEXTRA",
                         jtrack::DCAX,
                         jtrack::DCAY,
                         jtrack::DCAZ,
                         jtrack::DCAXY,
                         jtrack::DCAXYZ,
                         jtrack::SigmaDCAZ,
                         jtrack::SigmaDCAXY,
                         jtrack::SigmaDCAXYZ,
                         jtrack::Sigma1Pt);

DECLARE_SOA_TABLE_STAGED(JTrackPIs, "JTRACKPI",
                         jtrack::TrackId);

namespace jemctrack
{
DECLARE_SOA_INDEX_COLUMN(JTrack, track);
DECLARE_SOA_COLUMN(EtaEMCAL, etaEmcal, float);
DECLARE_SOA_COLUMN(PhiEMCAL, phiEmcal, float);
} // namespace jemctrack

DECLARE_SOA_TABLE_STAGED(JEMCTracks, "JEMCTrack",
                         jemctrack::JTrackId,
                         jemctrack::EtaEMCAL,
                         jemctrack::PhiEMCAL);

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

DECLARE_SOA_TABLE_STAGED(JMcParticles, "JMCPARTICLE",
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
using StoredJMcParticle = StoredJMcParticles::iterator;

DECLARE_SOA_TABLE_STAGED(JMcParticlePIs, "JMCPARTICLEPI",
                         jmcparticle::McParticleId);

namespace jmctracklb
{
DECLARE_SOA_INDEX_COLUMN(JMcParticle, mcParticle);
}

DECLARE_SOA_TABLE_STAGED(JMcTrackLbs, "JMCTRACKLB", //! Table joined to the track table containing the MC index
                         jmctracklb::JMcParticleId);

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

DECLARE_SOA_TABLE_STAGED(JClusters, "JCLUSTER", //!
                         o2::soa::Index<>, jcluster::JCollisionId, jcluster::ID, jcluster::Energy,
                         jcluster::CoreEnergy, jcluster::RawEnergy, jcluster::Eta, jcluster::Phi,
                         jcluster::M02, jcluster::M20, jcluster::NCells, jcluster::Time,
                         jcluster::IsExotic, jcluster::DistanceToBadChannel, jcluster::NLM, jcluster::Definition,
                         jcluster::LeadingCellEnergy, jcluster::SubleadingCellEnergy, jcluster::LeadingCellNumber, jcluster::SubleadingCellNumber);

using JCluster = JClusters::iterator;
using StoredJCluster = StoredJClusters::iterator;

DECLARE_SOA_TABLE_STAGED(JClusterPIs, "JCLUSTERPI",
                         jcluster::EMCALClusterId);

DECLARE_SOA_TABLE_STAGED(JClusterTracks, "JCLUSTERTRACK", //!
                         jcluster::JTrackIds);

namespace jclusterhadroniccorrection
{
DECLARE_SOA_COLUMN(EnergyCorrectedOneTrack1, energyCorrectedOneTrack1, float);   //! with hadronic correction fraction (100%) for one matched track
DECLARE_SOA_COLUMN(EnergyCorrectedOneTrack2, energyCorrectedOneTrack2, float);   //! with hadronic correction fraction (70%) for one matched track - systematic studies
DECLARE_SOA_COLUMN(EnergyCorrectedAllTracks1, energyCorrectedAllTracks1, float); //! with hadronic correction fraction (100%) for all matched tracks
DECLARE_SOA_COLUMN(EnergyCorrectedAllTracks2, energyCorrectedAllTracks2, float); //! with hadronic correction fraction (70%) for all matched tracks - for systematic studies
} // namespace jclusterhadroniccorrection

DECLARE_SOA_TABLE_STAGED(JClustersCorrectedEnergies, "JCLUSTERCORRE",            //! if this table changes it needs to be reflected in FastJetUtilities.h!!
                         jclusterhadroniccorrection::EnergyCorrectedOneTrack1,   // corrected cluster energy for 1 matched track (f = 100%)
                         jclusterhadroniccorrection::EnergyCorrectedOneTrack2,   // corrected cluster energy for 1 matched track (f = 70%)
                         jclusterhadroniccorrection::EnergyCorrectedAllTracks1,  // corrected cluster energy for all matched tracks (f = 100%)
                         jclusterhadroniccorrection::EnergyCorrectedAllTracks2); // corrected cluster energy for all matched tracks (f = 70%)

namespace jmcclusterlb
{
DECLARE_SOA_ARRAY_INDEX_COLUMN(JMcParticle, mcParticles);
DECLARE_SOA_COLUMN(AmplitudeA, amplitudeA, std::vector<float>);
} // namespace jmcclusterlb

DECLARE_SOA_TABLE_STAGED(JMcClusterLbs, "JMCCLUSTERLB", //!
                         jmcclusterlb::JMcParticleIds, jmcclusterlb::AmplitudeA);

namespace jdummy
{

DECLARE_SOA_COLUMN(Dummy, dummy, bool);

} // namespace jdummy
DECLARE_SOA_TABLE_STAGED(JDummys, "JDUMMY",
                         o2::soa::Index<>,
                         jdummy::Dummy);

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_JETREDUCEDDATA_H_
