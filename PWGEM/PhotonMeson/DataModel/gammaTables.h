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

#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/CaloClusters.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include <TMath.h>

#ifndef PWGEM_PHOTONMESON_DATAMODEL_GAMMATABLES_H_
#define PWGEM_PHOTONMESON_DATAMODEL_GAMMATABLES_H_

// todo: declare more columns in this file dynamic or expression, atm I save a lot of redundant information
namespace o2::aod
{
namespace gammatrackreco
{
DECLARE_SOA_COLUMN(Eta, eta, float);                                                     //! Pseudorapidity
DECLARE_SOA_COLUMN(P, p, float);                                                         //! Total momentum in GeV/c
DECLARE_SOA_COLUMN(Phi, phi, float);                                                     //! Azimuthal angle
DECLARE_SOA_COLUMN(Pt, pt, float);                                                       //! Transversal momentum in GeV/c
DECLARE_SOA_COLUMN(PositivelyCharged, positivelyCharged, bool);                          //! True for positively charged track
DECLARE_SOA_COLUMN(TpcCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, float); //! Ratio  crossed rows over findable clusters
DECLARE_SOA_COLUMN(TpcFoundOverFindableCls, tpcFoundOverFindableCls, float);             //! Ratio of found over findable clusters
DECLARE_SOA_COLUMN(TpcNClsCrossedRows, tpcNClsCrossedRows, float);                       //! Number of crossed TPC Rows

} // namespace gammatrackreco

DECLARE_SOA_TABLE(V0DaughterTracks, "AOD", "V0TRACKS",
                  o2::soa::Index<>,
                  v0data::V0Id,
                  track::DcaXY,
                  gammatrackreco::Eta,
                  gammatrackreco::P,
                  gammatrackreco::Phi,
                  gammatrackreco::Pt,
                  gammatrackreco::PositivelyCharged,
                  gammatrackreco::TpcCrossedRowsOverFindableCls,
                  gammatrackreco::TpcFoundOverFindableCls,
                  gammatrackreco::TpcNClsCrossedRows,
                  pidtpc::TPCNSigmaEl,
                  pidtpc::TPCNSigmaPi,
                  track::TPCSignal);

namespace MCTracksTrue
{
DECLARE_SOA_COLUMN(SameMother, sameMother, bool); // Do the tracks have the same mother particle?
} // namespace MCTracksTrue

DECLARE_SOA_TABLE(V0DaughterMcParticles, "AOD", "MCTRACKTRUE",
                  mcparticle::PdgCode,
                  mcparticle::Px,
                  mcparticle::Py,
                  mcparticle::Pz,
                  MCTracksTrue::SameMother);

// Index table to associate V0DaughterTracks to the corresponding V0MCDaughterParticles if an entrie exists. This table contains an index column which contains the corresponding row number or -1, if there is no corresponding V0MCDaughterParticles
// This table is one column that is joinable with V0DaughterTracks as far as i understand
namespace MCParticleTrueIndex
{
DECLARE_SOA_INDEX_COLUMN(V0DaughterMcParticle, v0DaughterMcParticle);
} // namespace MCParticleTrueIndex

// DECLARE_SOA_INDEX_TABLE_USER(MCTrackIndex, V0MCDaughterParticles, "MCTRACKINDEX", MCParticleTrueIndex::V0DaughterTrackId);
DECLARE_SOA_TABLE(MCParticleIndex, "AOD", "MCPARTICLEINDEX", MCParticleTrueIndex::V0DaughterMcParticleId);

namespace gammarecalculated
{
DECLARE_SOA_COLUMN(RecalculatedVtxX, recalculatedVtxX, float); //! Recalculated conversion point
DECLARE_SOA_COLUMN(RecalculatedVtxY, recalculatedVtxY, float); //! Recalculated conversion point
DECLARE_SOA_COLUMN(RecalculatedVtxZ, recalculatedVtxZ, float); //! Recalculated conversion point
DECLARE_SOA_DYNAMIC_COLUMN(RecalculatedVtxR, recalculatedVtxR, [](float x, float y) { return TMath::Sqrt(x * x + y * y); });
} // namespace gammarecalculated

DECLARE_SOA_TABLE(V0Recalculated, "AOD", "V0RECALCULATED",
                  gammarecalculated::RecalculatedVtxX,
                  gammarecalculated::RecalculatedVtxY,
                  gammarecalculated::RecalculatedVtxZ,
                  gammarecalculated::RecalculatedVtxR<o2::aod::gammarecalculated::RecalculatedVtxX, o2::aod::gammarecalculated::RecalculatedVtxY>);

namespace gammamctrue
{
DECLARE_SOA_COLUMN(P, p, float); //! Absolute momentum in GeV/c
} // namespace gammamctrue

DECLARE_SOA_TABLE(McDaughterTrue, "AOD", "MCDAUTRUE",
                  gammamctrue::P);

namespace gammamctrue
{
DECLARE_SOA_COLUMN(Gamma, gamma, int64_t);       //! Used as reference for the daughters
DECLARE_SOA_COLUMN(NDaughters, nDaughters, int); //! Number of daughters
DECLARE_SOA_COLUMN(Eta, eta, float);             //! Pseudorapidity
DECLARE_SOA_COLUMN(Phi, phi, float);             //! Angle phi in rad
DECLARE_SOA_COLUMN(Pt, pt, float);               //! Transversal momentum in GeV/c
DECLARE_SOA_COLUMN(Y, y, float);                 //! Rapidity

DECLARE_SOA_COLUMN(ConversionX, conversionX, float); //! x of conversion point in cm
DECLARE_SOA_COLUMN(ConversionY, conversionY, float); //! y of conversion point in cm
DECLARE_SOA_COLUMN(ConversionZ, conversionZ, float); //! z of conversion point in cm
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);       //! 2d radius of conversion point
// DECLARE_SOA_INDEX_COLUMN(McDaughterTrue, mcDaughterTrue);
DECLARE_SOA_INDEX_COLUMN_FULL(McDaughterTrueOne, mcDaughterTrueOne, int, McDaughterTrue, "_One"); // this is a reference that points to the entry in the McDaughterTrues table
DECLARE_SOA_INDEX_COLUMN_FULL(McDaughterTrueTwo, mcDaughterTrueTwo, int, McDaughterTrue, "_Two"); // this is a reference that points to the entry in the McDaughterTrues table
} // namespace gammamctrue

DECLARE_SOA_TABLE(McGammasTrue, "AOD", "MCGATRUE",
                  o2::soa::Index<>,
                  mcparticle::McCollisionId,
                  gammamctrue::Gamma,
                  v0data::V0Id, // reference to reconstructed v0 (if its a task with reconstucted info)
                  mcparticle::PdgCode, mcparticle::StatusCode, mcparticle::Flags,
                  mcparticle::Px, mcparticle::Py, mcparticle::Pz,
                  mcparticle::Vx, mcparticle::Vy, mcparticle::Vz, mcparticle::Vt,
                  gammamctrue::NDaughters,
                  // todo: make those expression columns in an extended table
                  gammamctrue::Eta, gammamctrue::Phi, gammamctrue::P, gammamctrue::Pt, gammamctrue::Y,
                  gammamctrue::ConversionX, gammamctrue::ConversionY, gammamctrue::ConversionZ,
                  gammamctrue::V0Radius,

                  // Index columns
                  gammamctrue::McDaughterTrueOneId,
                  gammamctrue::McDaughterTrueTwoId,

                  // Dynamic columns
                  mcparticle::ProducedByGenerator<mcparticle::Flags>,
                  mcparticle::FromBackgroundEvent<mcparticle::Flags>,
                  mcparticle::GetGenStatusCode<mcparticle::Flags, mcparticle::StatusCode>,
                  mcparticle::GetProcess<mcparticle::Flags, mcparticle::StatusCode>,
                  mcparticle::IsPhysicalPrimary<mcparticle::Flags>);

namespace gammacaloreco
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                        //! collisionID used as index for matched clusters
DECLARE_SOA_INDEX_COLUMN(BC, bc);                                      //! bunch crossing ID used as index for ambiguous clusters
DECLARE_SOA_COLUMN(ID, id, int);                                       //! cluster ID identifying cluster in event
DECLARE_SOA_COLUMN(Energy, energy, float);                             //! cluster energy (GeV)
DECLARE_SOA_COLUMN(CoreEnergy, coreEnergy, float);                     //! cluster core energy (GeV)
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
// DECLARE_SOA_INDEX_COLUMN(Calo, calo);                                  //! linked to calo cells
// DECLARE_SOA_INDEX_COLUMN(Track, track); //! linked to Track table only for tracks that were matched
// DECLARE_SOA_EXPRESSION_COLUMN(Pt, pt, float, //! transverse momentum of a photon candidate 2
//                               (gammacaloreco::energy * 2.f) / (nexp(gammacaloreco::eta) + nexp(gammacaloreco::eta * -1.f)));
// DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float pz, float E) { return atanh(pz / E); });  //! pseudorapidity of the cluster
// DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float px, float py) { return atan2(py, px); }); //! phi angle of the cluster
} // namespace gammacaloreco

DECLARE_SOA_TABLE(SkimEMCClusters, "AOD", "SKIMEMCCLUSTERS", //!
                  o2::soa::Index<>, gammacaloreco::CollisionId, gammacaloreco::ID, gammacaloreco::Energy, gammacaloreco::CoreEnergy,
                  gammacaloreco::Eta, gammacaloreco::Phi, gammacaloreco::M02, gammacaloreco::M20, gammacaloreco::NCells, gammacaloreco::Time,
                  gammacaloreco::IsExotic, gammacaloreco::DistanceToBadChannel, gammacaloreco::NLM, gammacaloreco::Definition);

} // namespace o2::aod
#endif // PWGEM_PHOTONMESON_DATAMODEL_GAMMATABLES_H_
