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
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/CaloClusters.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include <TMath.h>

#ifndef PWGEM_PHOTONMESON_DATAMODEL_GAMMATABLES_H_
#define PWGEM_PHOTONMESON_DATAMODEL_GAMMATABLES_H_

// todo: declare more columns in this file dynamic or expression, atm I save a lot of redundant information
namespace o2::aod
{

namespace emreducedevent
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision); //!
DECLARE_SOA_COLUMN(Tag, tag, uint64_t);         //!  Bit-field for storing event information (e.g. high level info, cut decisions)
DECLARE_SOA_COLUMN(NgammaPCM, ngpcm, int);
DECLARE_SOA_COLUMN(NgammaPHOS, ngphos, int);
DECLARE_SOA_COLUMN(NgammaEMC, ngemc, int);

} // namespace emreducedevent
DECLARE_SOA_TABLE(EMReducedEvents, "AOD", "EMREDUCEDEVENT", //!   Main event information table
                  o2::soa::Index<>, emreducedevent::CollisionId, emreducedevent::Tag, bc::RunNumber, evsel::Sel8,
                  collision::PosX, collision::PosY, collision::PosZ,
                  collision::NumContrib, collision::CollisionTime, collision::CollisionTimeRes,
                  emreducedevent::NgammaPCM, emreducedevent::NgammaPHOS, emreducedevent::NgammaEMC);
using EMReducedEvent = EMReducedEvents::iterator;

namespace v0leg
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);   //!
DECLARE_SOA_INDEX_COLUMN(Track, track);           //!
DECLARE_SOA_COLUMN(Sign, sign, int);              //!
DECLARE_SOA_COLUMN(IsAmbTrack, isAmbTrack, bool); //!

} // namespace v0leg
// reconstructed v0 information
DECLARE_SOA_TABLE(V0Legs, "AOD", "V0LEG", //!
                  o2::soa::Index<>, v0leg::CollisionId,
                  v0leg::TrackId, v0leg::Sign, v0leg::IsAmbTrack,
                  track::Pt, track::Eta, track::Phi,
                  track::DcaXY, track::DcaZ,
                  track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows,
                  track::TPCChi2NCl, track::TPCInnerParam,
                  track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr,
                  track::ITSClusterMap, track::ITSChi2NCl, track::DetectorMap,

                  // dynamic column
                  track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::ITSNCls<track::ITSClusterMap>,
                  track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>,
                  track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>);

// iterators
using V0Leg = V0Legs::iterator;

namespace v0photon
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                         //!
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, V0Legs, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, V0Legs, "_Neg"); //!
DECLARE_SOA_COLUMN(Vx, vx, float);                                      //!
DECLARE_SOA_COLUMN(Vy, vy, float);                                      //!
DECLARE_SOA_COLUMN(Vz, vz, float);                                      //!
DECLARE_SOA_COLUMN(PxPosAtSV, pxposatsv, float);                        //!
DECLARE_SOA_COLUMN(PyPosAtSV, pyposatsv, float);                        //!
DECLARE_SOA_COLUMN(PzPosAtSV, pzposatsv, float);                        //!
DECLARE_SOA_COLUMN(PxNegAtSV, pxnegatsv, float);                        //!
DECLARE_SOA_COLUMN(PyNegAtSV, pynegatsv, float);                        //!
DECLARE_SOA_COLUMN(PzNegAtSV, pznegatsv, float);                        //!
DECLARE_SOA_COLUMN(CosPA, cospa, float);                                //!
DECLARE_SOA_COLUMN(PCA, pca, float);                                    //!

DECLARE_SOA_DYNAMIC_COLUMN(Px, px, [](float pxpos, float pxneg) -> float { return pxpos + pxneg; });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, [](float pypos, float pyneg) -> float { return pypos + pyneg; });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, [](float pzpos, float pzneg) -> float { return pzpos + pzneg; });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float pxpos, float pypos, float pxneg, float pyneg) -> float { return RecoDecay::sqrtSumOfSquares(pxpos + pxneg, pypos + pyneg); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::eta(array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float pxpos, float pypos, float pxneg, float pyneg) -> float { return RecoDecay::phi(pxpos + pxneg, pypos + pyneg); });

} // namespace v0photon
// reconstructed v0 information
DECLARE_SOA_TABLE(V0Photons, "AOD", "V0PHOTON", //!
                  o2::soa::Index<>, v0photon::CollisionId, v0photon::PosTrackId, v0photon::NegTrackId,
                  v0photon::Vx, v0photon::Vy, v0photon::Vz,
                  v0photon::PxPosAtSV, v0photon::PyPosAtSV, v0photon::PzPosAtSV,
                  v0photon::PxNegAtSV, v0photon::PyNegAtSV, v0photon::PzNegAtSV,
                  v0photon::CosPA, v0photon::PCA,

                  // dynamic column
                  v0photon::Px<v0photon::PxPosAtSV, v0photon::PxNegAtSV>,
                  v0photon::Py<v0photon::PyPosAtSV, v0photon::PyNegAtSV>,
                  v0photon::Pz<v0photon::PzPosAtSV, v0photon::PzNegAtSV>,
                  v0photon::Pt<v0photon::PxPosAtSV, v0photon::PyPosAtSV, v0photon::PxNegAtSV, v0photon::PyNegAtSV>,
                  v0photon::Eta<v0photon::PxPosAtSV, v0photon::PyPosAtSV, v0photon::PzPosAtSV, v0photon::PxNegAtSV, v0photon::PyNegAtSV, v0photon::PzNegAtSV>,
                  v0photon::Phi<v0photon::PxPosAtSV, v0photon::PyPosAtSV, v0photon::PxNegAtSV, v0photon::PyNegAtSV>,
                  v0data::V0Radius<v0photon::Vx, v0photon::Vy>,
                  v0data::Alpha<v0photon::PxPosAtSV, v0photon::PyPosAtSV, v0photon::PzPosAtSV, v0photon::PxNegAtSV, v0photon::PyNegAtSV, v0photon::PzNegAtSV>,
                  v0data::QtArm<v0photon::PxPosAtSV, v0photon::PyPosAtSV, v0photon::PzPosAtSV, v0photon::PxNegAtSV, v0photon::PyNegAtSV, v0photon::PzNegAtSV>,
                  v0data::PsiPair<v0photon::PxPosAtSV, v0photon::PyPosAtSV, v0photon::PzPosAtSV, v0photon::PxNegAtSV, v0photon::PyNegAtSV, v0photon::PzNegAtSV>,

                  // invariant mass mee
                  v0data::MGamma<v0photon::PxPosAtSV, v0photon::PyPosAtSV, v0photon::PzPosAtSV, v0photon::PxNegAtSV, v0photon::PyNegAtSV, v0photon::PzNegAtSV>);
// iterators
using V0Photon = V0Photons::iterator;

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
DECLARE_SOA_TABLE(SkimEMCClusters, "AOD", "SKIMEMCCLUSTERS", //! table of skimmed EMCal clusters
                  o2::soa::Index<>, gammacaloreco::CollisionId, gammacaloreco::ID, gammacaloreco::Energy, gammacaloreco::CoreEnergy,
                  gammacaloreco::Eta, gammacaloreco::Phi, gammacaloreco::M02, gammacaloreco::M20, gammacaloreco::NCells, gammacaloreco::Time,
                  gammacaloreco::IsExotic, gammacaloreco::DistanceToBadChannel, gammacaloreco::NLM, gammacaloreco::Definition);

namespace phoscluster
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                                     //! collisionID used as index for matched clusters
DECLARE_SOA_INDEX_COLUMN(BC, bc);                                                   //! bunch crossing ID used as index for ambiguous clusters
DECLARE_SOA_INDEX_COLUMN_FULL(MatchedTrack, matchedTrack, int, Tracks, "_Matched"); //! matched track index
DECLARE_SOA_COLUMN(E, e, float);                                                    //! cluster energy (GeV)
DECLARE_SOA_COLUMN(X, x, float);                                                    //! cluster hit position in ALICE global coordinate
DECLARE_SOA_COLUMN(Y, y, float);                                                    //! cluster hit position in ALICE global coordinate
DECLARE_SOA_COLUMN(Z, z, float);                                                    //! cluster hit position in ALICE global coordinate
DECLARE_SOA_COLUMN(M02, m02, float);                                                //! shower shape long axis
DECLARE_SOA_COLUMN(M20, m20, float);                                                //! shower shape short axis
DECLARE_SOA_COLUMN(NCells, nCells, int);                                            //! number of cells in cluster
DECLARE_SOA_COLUMN(Time, time, float);                                              //! cluster time (ns)
DECLARE_SOA_COLUMN(DistanceToBadChannel, distanceToBadChannel, float);              //! distance to bad channel in cm
DECLARE_SOA_COLUMN(NLM, nlm, int);                                                  //! number of local maxima
// DECLARE_SOA_COLUMN(TrackEta, tracketa, float);                                      //! eta of the matched track
// DECLARE_SOA_COLUMN(TrackPhi, trackphi, float);                                      //! phi of the matched track
// DECLARE_SOA_COLUMN(TrackP, trackp, float);                                          //! momentum of the matched track
// DECLARE_SOA_COLUMN(TrackPt, trackpt, float);                                        //! pt of the matched track
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, [](float e, float x, float y, float z, float m) -> float { return x / RecoDecay::sqrtSumOfSquares(x, y, z) * sqrt(e * e - m * m); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, [](float e, float x, float y, float z, float m) -> float { return y / RecoDecay::sqrtSumOfSquares(x, y, z) * sqrt(e * e - m * m); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, [](float e, float x, float y, float z, float m) -> float { return z / RecoDecay::sqrtSumOfSquares(x, y, z) * sqrt(e * e - m * m); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float e, float x, float y, float z, float m) -> float { return RecoDecay::sqrtSumOfSquares(x, y) / RecoDecay::sqrtSumOfSquares(x, y, z) * sqrt(e * e - m * m); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float x, float y, float z) -> float { return RecoDecay::eta(array{x, y, z}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float x, float y) -> float { return RecoDecay::phi(x, y); });
} // namespace phoscluster

DECLARE_SOA_TABLE(PHOSClusters, "AOD", "PHOSCLUSTERS", //!
                  o2::soa::Index<>, phoscluster::CollisionId, phoscluster::BCId, phoscluster::MatchedTrackId,
                  phoscluster::E, phoscluster::X, phoscluster::Y, phoscluster::Z,
                  phoscluster::M02, phoscluster::M20, phoscluster::NCells,
                  phoscluster::Time, phoscluster::DistanceToBadChannel, phoscluster::NLM,
                  // phoscluster::TrackEta, phoscluster::TrackPhi, phoscluster::TrackP, phoscluster::TrackPt,
                  // dynamic column
                  phoscluster::Px<phoscluster::E, phoscluster::X, phoscluster::Y, phoscluster::Z>,
                  phoscluster::Py<phoscluster::E, phoscluster::X, phoscluster::Y, phoscluster::Z>,
                  phoscluster::Pz<phoscluster::E, phoscluster::X, phoscluster::Y, phoscluster::Z>,
                  phoscluster::Pt<phoscluster::E, phoscluster::X, phoscluster::Y, phoscluster::Z>,
                  phoscluster::Eta<phoscluster::X, phoscluster::Y, phoscluster::Z>,
                  phoscluster::Phi<phoscluster::X, phoscluster::Y>);
using PHOSCluster = PHOSClusters::iterator;

using SkimEMCCluster = SkimEMCClusters::iterator;

namespace gammacellsreco
{
DECLARE_SOA_INDEX_COLUMN(SkimEMCCluster, skimEMCCluster); //! SkimEMCClusterID
DECLARE_SOA_INDEX_COLUMN(Calo, calo);                     //! CaloID (CellID)
// DECLARE_SOA_INDEX_COLUMN(Track, track);                   //! CaloID (CellID)
DECLARE_SOA_COLUMN(TrackEta, tracketa, float); //! eta of the matched track
DECLARE_SOA_COLUMN(TrackPhi, trackphi, float); //! phi of the matched track
DECLARE_SOA_COLUMN(TrackP, trackp, float);     //! momentum of the matched track
DECLARE_SOA_COLUMN(TrackPt, trackpt, float);   //! pt of the matched track
} // namespace gammacellsreco

DECLARE_SOA_TABLE(SkimEMCCells, "AOD", "SKIMEMCCELLS",                                         //! table of link between skimmed EMCal clusters and their cells
                  o2::soa::Index<>, gammacellsreco::SkimEMCClusterId, gammacellsreco::CaloId); //!

DECLARE_SOA_TABLE(SkimEMCMTs, "AOD", "SKIMEMCMTS", //! table of link between skimmed EMCal clusters and their matched tracks
                  o2::soa::Index<>, gammacellsreco::SkimEMCClusterId, gammacellsreco::TrackEta,
                  gammacellsreco::TrackPhi, gammacellsreco::TrackP, gammacellsreco::TrackPt);

namespace gammareco
{
DECLARE_SOA_COLUMN(Method, method, int);                                         //! cut bit for PCM photon candidates
DECLARE_SOA_INDEX_COLUMN_FULL(SkimmedPCM, skimmedPCM, int, V0Photons, "");       //! reference to the gamma in the skimmed PCM table
DECLARE_SOA_INDEX_COLUMN_FULL(SkimmedPHOS, skimmedPHOS, int, PHOSClusters, "");  //! reference to the gamma in the skimmed PHOS table
DECLARE_SOA_INDEX_COLUMN_FULL(SkimmedEMC, skimmedEMC, int, SkimEMCClusters, ""); //! reference to the gamma in the skimmed EMCal table
DECLARE_SOA_COLUMN(PCMCutBit, pcmcutbit, uint64_t);                              //! cut bit for PCM photon candidates
DECLARE_SOA_COLUMN(PHOSCutBit, phoscutbit, uint64_t);                            //! cut bit for PHOS photon candidates
DECLARE_SOA_COLUMN(EMCCutBit, emccutbit, uint64_t);                              //! cut bit for EMCal photon candidates
} // namespace gammareco
DECLARE_SOA_TABLE(SkimGammas, "AOD", "SKIMGAMMAS", //! table of all gamma candidates (PCM, EMCal and PHOS) after cuts
                  o2::soa::Index<>, gammacaloreco::CollisionId, gammareco::Method,
                  gammacaloreco::Energy, gammacaloreco::Eta, gammacaloreco::Phi,
                  gammareco::SkimmedEMCId, gammareco::SkimmedPHOSId);
DECLARE_SOA_TABLE(SkimPCMCuts, "AOD", "SKIMPCMCUTS",                                  //! table of link between skimmed PCM photon candidates and their cuts
                  o2::soa::Index<>, gammareco::SkimmedPCMId, gammareco::PCMCutBit);   //!
DECLARE_SOA_TABLE(SkimPHOSCuts, "AOD", "SKIMPHOSCUTS",                                //! table of link between skimmed PHOS photon candidates and their cuts
                  o2::soa::Index<>, gammareco::SkimmedPHOSId, gammareco::PHOSCutBit); //!
DECLARE_SOA_TABLE(SkimEMCCuts, "AOD", "SKIMEMCCUTS",                                  //! table of link between skimmed EMCal photon candidates and their cuts
                  o2::soa::Index<>, gammareco::SkimmedEMCId, gammareco::EMCCutBit);   //!
} // namespace o2::aod

#endif // PWGEM_PHOTONMESON_DATAMODEL_GAMMATABLES_H_
