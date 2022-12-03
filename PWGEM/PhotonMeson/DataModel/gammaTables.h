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

#ifndef PWGEM_PHOTONMESON_DATAMODEL_GAMMATABLES_H_
#define PWGEM_PHOTONMESON_DATAMODEL_GAMMATABLES_H_

#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include <TMath.h>

// todo: declare more columns in this file dynamic or expression, atm I save a lot of redundant information
namespace o2::aod
{

namespace emreducedevent
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision); //!
DECLARE_SOA_COLUMN(Tag, tag, uint64_t);         //!  Bit-field for storing event information (e.g. high level info, cut decisions)
} // namespace emreducedevent
DECLARE_SOA_TABLE(EMReducedEvents, "AOD", "EMREDUCEDEVENT", //!   Main event information table
                  o2::soa::Index<>, emreducedevent::CollisionId, emreducedevent::Tag, bc::RunNumber, evsel::Sel8,
                  collision::PosX, collision::PosY, collision::PosZ,
                  collision::NumContrib, collision::CollisionTime, collision::CollisionTimeRes);
using EMReducedEvent = EMReducedEvents::iterator;

namespace v0leg
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);           //!
DECLARE_SOA_INDEX_COLUMN(EMReducedEvent, emreducedevent); //!
DECLARE_SOA_INDEX_COLUMN(Track, track);                   //!
DECLARE_SOA_COLUMN(Sign, sign, int);                      //!
DECLARE_SOA_COLUMN(IsAmbTrack, isAmbTrack, bool);         //!
// DECLARE_SOA_COLUMN(NumBC, numBC, unsigned int);                                            //!
// DECLARE_SOA_COLUMN(NumColl, numColl, unsigned int);                                            //!

} // namespace v0leg
// reconstructed v0 information
DECLARE_SOA_TABLE(V0Legs, "AOD", "V0LEG", //!
                  o2::soa::Index<>, v0leg::CollisionId, v0leg::EMReducedEventId,
                  v0leg::TrackId, v0leg::Sign, v0leg::IsAmbTrack,
                  // v0leg::NumBC, v0leg::NumColl,
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
DECLARE_SOA_INDEX_COLUMN(EMReducedEvent, emreducedevent);               //!
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
                  o2::soa::Index<>, v0photon::CollisionId, v0photon::EMReducedEventId, v0photon::PosTrackId, v0photon::NegTrackId,
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
} // namespace o2::aod

#endif // PWGEM_PHOTONMESON_DATAMODEL_GAMMATABLES_H_
