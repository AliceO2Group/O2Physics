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

/// \file LFResonanceTables.h
/// \brief Definitions of tables of resonance decay candidates
///
/// Inspired by StrangenessTables.h, FemtoDerived.h
///
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>

#ifndef O2_ANALYSIS_LFRESONANCETABLES_H_
#define O2_ANALYSIS_LFRESONANCETABLES_H_

#include <cmath>

#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
/// Resonance Collisions
namespace resocollision
{
DECLARE_SOA_COLUMN(MultV0M, multV0M, float);       //! V0M multiplicity
DECLARE_SOA_COLUMN(Sphericity, sphericity, float); //! Sphericity of the event
} // namespace resocollision
DECLARE_SOA_TABLE(ResoCollisions, "AOD", "RESOCOL",
                  o2::soa::Index<>,
                  o2::aod::collision::PosX,
                  o2::aod::collision::PosY,
                  o2::aod::collision::PosZ,
                  resocollision::MultV0M,
                  resocollision::Sphericity,
                  timestamp::Timestamp);
using ResoCollision = ResoCollisions::iterator;

// Resonance Daughters
// inspired from PWGCF/DataModel/FemtoDerived.h
namespace resodaughter
{
/// Distinuishes the different daughter types
enum DaughterType {
  kTrack,         //! Track
  kV0,            //! V0
  kCascade,       //! Cascade
  kNDaughterTypes //! Number of Daughter types
};
static constexpr std::string_view DaughterTypeName[kNDaughterTypes] = {"Tracks", "V0", "Cascade"}; //! Naming of the different particle types

enum PDGtype {
    kPion = BIT(0),
    kKaon = BIT(1),
    kProton = BIT(2)
};

DECLARE_SOA_INDEX_COLUMN(ResoCollision, resoCollision);
DECLARE_SOA_COLUMN(Pt, pt, float);                                   //! p_T (GeV/c)
DECLARE_SOA_COLUMN(Eta, eta, float);                                 //! Eta
DECLARE_SOA_COLUMN(Phi, phi, float);                                 //! Phi
DECLARE_SOA_COLUMN(PartType, partType, uint8_t);                     //! Type of the particle, according to resodaughter::ParticleType
DECLARE_SOA_COLUMN(TempFitVar, tempFitVar, float);                   //! Observable for the template fitting (Track: DCA_xy, V0: CPA)
DECLARE_SOA_COLUMN(Indices, indices, int[2]);                        //! Field for the track indices to remove auto-correlations
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                              //! Sign of the track charge
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, uint8_t); //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);                             //! DCA_xy
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);                               //! DCA_z
DECLARE_SOA_COLUMN(X, x, float);                                     //! x position of the track
DECLARE_SOA_COLUMN(Alpha, alpha, float);                             //! alpha position of the track
DECLARE_SOA_COLUMN(TPCPIDselection, tpcPIDselection, uint8_t);       //! TPC PID selection
DECLARE_SOA_COLUMN(TOFPIDselection, tofPIDselection, uint8_t);       //! TOF PID selection
DECLARE_SOA_COLUMN(TPCnSigmaPi, tpcNSigmaPi, float);                 //! Pion TPC nSigma
DECLARE_SOA_COLUMN(TPCnSigmaKa, tpcNSigmaKa, float);                 //! Kaon TPC nSigma
DECLARE_SOA_COLUMN(TPCnSigmaPr, tpcNSigmaPr, float);                 //! Proton TPC nSigma
DECLARE_SOA_COLUMN(TOFnSigmaPi, tofNSigmaPi, float);                 //! Pion TOF nSigma
DECLARE_SOA_COLUMN(TOFnSigmaKa, tofNSigmaKa, float);                 //! Kaon TOF nSigma
DECLARE_SOA_COLUMN(TOFnSigmaPr, tofNSigmaPr, float);                 //! Proton TOF nSigma
DECLARE_SOA_COLUMN(DaughDCA, daughDCA, float);                       //! DCA between daughters
DECLARE_SOA_COLUMN(MLambda, mLambda, float);                         //! The invariant mass of V0 candidate, assuming lambda
DECLARE_SOA_COLUMN(MAntiLambda, mAntiLambda, float);                 //! The invariant mass of V0 candidate, assuming antilambda
DECLARE_SOA_COLUMN(TransRadius, transRadius, float);                 //! Transverse radius of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxX, decayVtxX, float);                     //! X position of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxY, decayVtxY, float);                     //! Y position of the decay vertex
DECLARE_SOA_COLUMN(DecayVtxZ, decayVtxZ, float);                     //! Z position of the decay vertex

DECLARE_SOA_DYNAMIC_COLUMN(Theta, theta, //! Compute the theta of the track
                           [](float eta) -> float {
                             return 2.f * std::atan(std::exp(-eta));
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, //! Compute the momentum in x in GeV/c
                           [](float pt, float phi) -> float {
                             return pt * std::sin(phi);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //! Compute the momentum in y in GeV/c
                           [](float pt, float phi) -> float {
                             return pt * std::cos(phi);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //! Compute the momentum in z in GeV/c
                           [](float pt, float eta) -> float {
                             return pt * std::sinh(eta);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //! Compute the overall momentum in GeV/c
                           [](float pt, float eta) -> float {
                             return pt * std::cosh(eta);
                           });
} // namespace resodaughter
DECLARE_SOA_TABLE(ResoDaughters, "AOD", "RESODAUGHTERS",
                  o2::soa::Index<>,
                  resodaughter::ResoCollisionId,
                  resodaughter::Pt,
                  resodaughter::Eta,
                  resodaughter::Phi,
                  resodaughter::PartType,
                  resodaughter::TempFitVar,
                  resodaughter::Indices,
                  resodaughter::Sign,
                  resodaughter::TPCNClsCrossedRows,
                  resodaughter::DcaXY,
                  resodaughter::DcaZ,
                  resodaughter::X,
                  resodaughter::Alpha,
                  resodaughter::TPCPIDselection,
                  resodaughter::TOFPIDselection,
                  resodaughter::TPCnSigmaPi,
                  resodaughter::TPCnSigmaKa,
                  resodaughter::TPCnSigmaPr,
                  resodaughter::TOFnSigmaPi,
                  resodaughter::TOFnSigmaKa,
                  resodaughter::TOFnSigmaPr,
                  resodaughter::DaughDCA,
                  resodaughter::MLambda,
                  resodaughter::MAntiLambda,
                  resodaughter::TransRadius,
                  resodaughter::DecayVtxX,
                  resodaughter::DecayVtxY,
                  resodaughter::DecayVtxZ,
                  resodaughter::Theta<resodaughter::Eta>,
                  resodaughter::Px<resodaughter::Pt, resodaughter::Phi>,
                  resodaughter::Py<resodaughter::Pt, resodaughter::Phi>,
                  resodaughter::Pz<resodaughter::Pt, resodaughter::Eta>,
                  resodaughter::P<resodaughter::Pt, resodaughter::Eta>);
using ResoDaughter = ResoDaughters::iterator;

/// Reconstruction of track-track decay resonance candidates
///
namespace reso2trktrkdata
{
// Needed to have shorter table that does not rely on existing one (filtering!)
DECLARE_SOA_INDEX_COLUMN_FULL(Track1, track1, int, Tracks, "_1"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(Track2, track2, int, Tracks, "_2"); //!
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                   //!

// General two track properties: position, momentum
DECLARE_SOA_COLUMN(Trk1Sign, trk1Sign, int); //! track1 sign
DECLARE_SOA_COLUMN(Trk2Sign, trk2Sign, int); //! track2 sign
DECLARE_SOA_COLUMN(PxTrk1, pxTrk1, float);   //! track1 px at min
DECLARE_SOA_COLUMN(PyTrk1, pyTrk1, float);   //! track1 py at min
DECLARE_SOA_COLUMN(PzTrk1, pzTrk1, float);   //! track1 pz at min
DECLARE_SOA_COLUMN(PxTrk2, pxTrk2, float);   //! track2 px at min
DECLARE_SOA_COLUMN(PyTrk2, pyTrk2, float);   //! track2 py at min
DECLARE_SOA_COLUMN(PzTrk2, pzTrk2, float);   //! track2 pz at min

// Saved from finding: DCAs
DECLARE_SOA_COLUMN(DCArTrk1ToPV, dcarTrk1toPV, float); //! DCA track1 to PVr
DECLARE_SOA_COLUMN(DCArTrk2ToPV, dcarTrk2toPV, float); //! DCA track2 to PVr
DECLARE_SOA_COLUMN(DCAzTrk1ToPV, dcazTrk1toPV, float); //! DCA track1 to PVz
DECLARE_SOA_COLUMN(DCAzTrk2ToPV, dcazTrk2toPV, float); //! DCA track2 to PVz

// Derived expressions
// Momenta
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //! resonance pT
                           [](float pxTrk1, float pyTrk1, float pxTrk2, float pyTrk2) -> float { return RecoDecay::sqrtSumOfSquares(pxTrk1 + pxTrk2, pyTrk1 + pyTrk2); });

// Psi pair angle: angle between the plane defined by the electron and positron momenta and the xy plane
DECLARE_SOA_DYNAMIC_COLUMN(PsiPair, psiPair, //! psi pair angle
                           [](float pxTrk1, float pyTrk1, float pzTrk1, float pxTrk2, float pyTrk2, float pzTrk2) {
                             auto clipToPM1 = [](float x) { return x < -1.f ? -1.f : (x > 1.f ? 1.f : x); };
                             float ptot2 = RecoDecay::p2(pxTrk1, pyTrk1, pzTrk1) * RecoDecay::p2(pxTrk2, pyTrk2, pzTrk2);
                             float argcos = RecoDecay::dotProd(array{pxTrk1, pyTrk1, pzTrk1}, array{pxTrk2, pyTrk2, pzTrk2}) / std::sqrt(ptot2);
                             float thetaPos = std::atan2(RecoDecay::sqrtSumOfSquares(pxTrk1, pyTrk1), pzTrk1);
                             float thetaNeg = std::atan2(RecoDecay::sqrtSumOfSquares(pxTrk2, pyTrk2), pzTrk2);
                             float argsin = (thetaNeg - thetaPos) / std::acos(clipToPM1(argcos));
                             return std::asin(clipToPM1(argsin));
                           });

DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, //! Resonance eta
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::eta(array{Px, Py, Pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, //! Resonance phi
                           [](float Px, float Py) -> float { return RecoDecay::phi(Px, Py); });

DECLARE_SOA_DYNAMIC_COLUMN(Track2Pt, track2Pt, //! Track2 daughter pT
                           [](float pxTrk2, float pyTrk2) -> float { return RecoDecay::sqrtSumOfSquares(pxTrk2, pyTrk2); });
DECLARE_SOA_DYNAMIC_COLUMN(Track1Pt, track1Pt, //! Track1 daughter pT
                           [](float pxTrk1, float pyTrk1) -> float { return RecoDecay::sqrtSumOfSquares(pxTrk1, pyTrk1); });
DECLARE_SOA_DYNAMIC_COLUMN(Track2Eta, track2Eta, //! Track2 daughter eta
                           [](float pxTrk2, float pyTrk2, float pzTrk2) -> float { return RecoDecay::eta(array{pxTrk2, pyTrk2, pzTrk2}); });
DECLARE_SOA_DYNAMIC_COLUMN(Track2Phi, track2Phi, //! Track2 daughter phi
                           [](float pxTrk2, float pyTrk2) -> float { return RecoDecay::phi(pxTrk2, pyTrk2); });
DECLARE_SOA_DYNAMIC_COLUMN(Track1Eta, track1Eta, //! Track1 daughter eta
                           [](float pxTrk1, float pyTrk1, float pzTrk1) -> float { return RecoDecay::eta(array{pxTrk1, pyTrk1, pzTrk1}); });
DECLARE_SOA_DYNAMIC_COLUMN(Track1Phi, track1Phi, //! Track1 daughter phi
                           [](float pxTrk1, float pyTrk1) -> float { return RecoDecay::phi(pxTrk1, pyTrk1); });

DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //! Resonance px
                              float, 1.f * aod::reso2trktrkdata::pxTrk1 + 1.f * aod::reso2trktrkdata::pxTrk2);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //! Resonance py
                              float, 1.f * aod::reso2trktrkdata::pyTrk1 + 1.f * aod::reso2trktrkdata::pyTrk2);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //! Resonance pz
                              float, 1.f * aod::reso2trktrkdata::pzTrk1 + 1.f * aod::reso2trktrkdata::pzTrk2);
} // namespace reso2trktrkdata

DECLARE_SOA_TABLE_FULL(StoredReso2TrackTrackDatas, "Reso2TrackTrackDatas", "AOD", "Reso2TTDATA", //!
                       o2::soa::Index<>, reso2trktrkdata::Track1Id, reso2trktrkdata::Track2Id, reso2trktrkdata::CollisionId,
                       reso2trktrkdata::Trk1Sign, reso2trktrkdata::Trk2Sign,
                       reso2trktrkdata::PxTrk1, reso2trktrkdata::PyTrk1, reso2trktrkdata::PzTrk1,
                       reso2trktrkdata::PxTrk2, reso2trktrkdata::PyTrk2, reso2trktrkdata::PzTrk2,

                       // Topological selections
                       reso2trktrkdata::DCArTrk1ToPV, reso2trktrkdata::DCArTrk2ToPV,
                       reso2trktrkdata::DCAzTrk1ToPV, reso2trktrkdata::DCAzTrk2ToPV,

                       // Dynamic columns
                       reso2trktrkdata::Pt<reso2trktrkdata::PxTrk1, reso2trktrkdata::PyTrk1, reso2trktrkdata::PxTrk2, reso2trktrkdata::PyTrk2>,
                       reso2trktrkdata::PsiPair<reso2trktrkdata::PxTrk1, reso2trktrkdata::PyTrk1, reso2trktrkdata::PzTrk1, reso2trktrkdata::PxTrk2, reso2trktrkdata::PyTrk2, reso2trktrkdata::PzTrk2>,

                       // Longitudinal
                       reso2trktrkdata::Eta<reso2trktrkdata::Px, reso2trktrkdata::Py, reso2trktrkdata::Pz>,
                       reso2trktrkdata::Phi<reso2trktrkdata::Px, reso2trktrkdata::Py>,
                       reso2trktrkdata::Track2Pt<reso2trktrkdata::PxTrk2, reso2trktrkdata::PyTrk2>,
                       reso2trktrkdata::Track1Pt<reso2trktrkdata::PxTrk1, reso2trktrkdata::PyTrk1>,
                       reso2trktrkdata::Track2Eta<reso2trktrkdata::PxTrk2, reso2trktrkdata::PyTrk2, reso2trktrkdata::PzTrk2>,
                       reso2trktrkdata::Track2Phi<reso2trktrkdata::PxTrk2, reso2trktrkdata::PyTrk2>,
                       reso2trktrkdata::Track1Eta<reso2trktrkdata::PxTrk1, reso2trktrkdata::PyTrk1, reso2trktrkdata::PzTrk1>,
                       reso2trktrkdata::Track1Phi<reso2trktrkdata::PxTrk1, reso2trktrkdata::PyTrk1>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(Reso2TrackTrackDatas, StoredReso2TrackTrackDatas, "Reso2TTDATAEXT", //!
                                reso2trktrkdata::Px, reso2trktrkdata::Py, reso2trktrkdata::Pz);     // the table name has here to be the one with EXT which is not nice and under study

using Reso2TrackTrackData = Reso2TrackTrackDatas::iterator;

using Reso2TracksExt = soa::Join<aod::FullTracks, aod::TracksDCA>;
using Reso2TracksMC = soa::Join<aod::FullTracks, McTrackLabels>;
using Reso2TracksPID = soa::Join<aod::FullTracks, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;
using Reso2TracksPIDExt = soa::Join<Reso2TracksPID, aod::TracksDCA>;

/// Reconstruction of track-v0 decay resonance candidates
///
namespace reso2trkv0data
{
// Needed to have shorter table that does not rely on existing one (filtering!)
DECLARE_SOA_INDEX_COLUMN_FULL(Track, track, int, Tracks, "_track"); //!
DECLARE_SOA_INDEX_COLUMN(V0, v0);                                   //!
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                     //!

// General V0, track properties: position, momentum
DECLARE_SOA_COLUMN(Sign, sign, int);         //!
DECLARE_SOA_COLUMN(PxPos, pxpos, float);     //!
DECLARE_SOA_COLUMN(PyPos, pypos, float);     //!
DECLARE_SOA_COLUMN(PzPos, pzpos, float);     //!
DECLARE_SOA_COLUMN(PxNeg, pxneg, float);     //!
DECLARE_SOA_COLUMN(PyNeg, pyneg, float);     //!
DECLARE_SOA_COLUMN(PzNeg, pzneg, float);     //!
DECLARE_SOA_COLUMN(PxTrack, pxtrack, float); //!
DECLARE_SOA_COLUMN(PyTrack, pytrack, float); //!
DECLARE_SOA_COLUMN(PzTrack, pztrack, float); //!
DECLARE_SOA_COLUMN(XV0, xV0, float);         //!
DECLARE_SOA_COLUMN(YV0, yV0, float);         //!
DECLARE_SOA_COLUMN(ZV0, zV0, float);         //!

// Saved from finding: DCAs
DECLARE_SOA_COLUMN(DCAV0Daughters, dcaV0daughters, float); //!
DECLARE_SOA_COLUMN(DCAPosToPV, dcapostopv, float);         //!
DECLARE_SOA_COLUMN(DCANegToPV, dcanegtopv, float);         //!
DECLARE_SOA_COLUMN(DCATrackToPV, dcatracktopv, float);     //!

// Derived expressions
// Momenta
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //!
                           [](float Px, float Py) -> float { return RecoDecay::sqrtSumOfSquares(Px, Py); });

// Length quantities
DECLARE_SOA_DYNAMIC_COLUMN(V0Radius, v0radius, //!
                           [](float xV0, float yV0) -> float { return RecoDecay::sqrtSumOfSquares(xV0, yV0); });

// CosPAs
DECLARE_SOA_DYNAMIC_COLUMN(V0CosPA, v0cosPA, //!
                           [](float XV0, float YV0, float ZV0, float PxV0, float PyV0, float PzV0, float pvX, float pvY, float pvZ) -> float { return RecoDecay::cpa(array{pvX, pvY, pvZ}, array{XV0, YV0, ZV0}, array{PxV0, PyV0, PzV0}); });
DECLARE_SOA_DYNAMIC_COLUMN(DCAV0ToPV, dcav0topv, //!
                           [](float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ) -> float { return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz)); });

// Calculated on the fly with mass assumption + dynamic tables
DECLARE_SOA_DYNAMIC_COLUMN(MLambda, mLambda, //! mass under lambda hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(array{array{pxpos, pypos, pzpos}, array{pxneg, pyneg, pzneg}}, array{RecoDecay::getMassPDG(kProton), RecoDecay::getMassPDG(kPiPlus)}); });
DECLARE_SOA_DYNAMIC_COLUMN(MAntiLambda, mAntiLambda, //! mass under antilambda hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(array{array{pxpos, pypos, pzpos}, array{pxneg, pyneg, pzneg}}, array{RecoDecay::getMassPDG(kPiPlus), RecoDecay::getMassPDG(kProton)}); });
DECLARE_SOA_DYNAMIC_COLUMN(MK0Short, mK0Short, //! mass under K0short hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(array{array{pxpos, pypos, pzpos}, array{pxneg, pyneg, pzneg}}, array{RecoDecay::getMassPDG(kPiPlus), RecoDecay::getMassPDG(kPiPlus)}); });

DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, //!
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::eta(array{Px, Py, Pz}); });

DECLARE_SOA_EXPRESSION_COLUMN(PxV0, pxV0, //!
                              float, 1.f * aod::reso2trkv0data::pxpos + 1.f * aod::reso2trkv0data::pxneg);
DECLARE_SOA_EXPRESSION_COLUMN(PyV0, pyV0, //!
                              float, 1.f * aod::reso2trkv0data::pypos + 1.f * aod::reso2trkv0data::pyneg);
DECLARE_SOA_EXPRESSION_COLUMN(PzV0, pzV0, //!
                              float, 1.f * aod::reso2trkv0data::pzpos + 1.f * aod::reso2trkv0data::pzneg);
DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //!
                              float, 1.f * aod::reso2trkv0data::pxpos + 1.f * aod::reso2trkv0data::pxneg + 1.f * aod::reso2trkv0data::pxtrack);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //!
                              float, 1.f * aod::reso2trkv0data::pypos + 1.f * aod::reso2trkv0data::pyneg + 1.f * aod::reso2trkv0data::pytrack);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //!
                              float, 1.f * aod::reso2trkv0data::pzpos + 1.f * aod::reso2trkv0data::pzneg + 1.f * aod::reso2trkv0data::pztrack);
} // namespace reso2trkv0data

DECLARE_SOA_TABLE(Reso2TrackV0Data, "AOD", "Reso2TVDATA", //!
                  o2::soa::Index<>, reso2trkv0data::TrackId, reso2trkv0data::V0Id, reso2trkv0data::CollisionId,

                  // General V0, track properties: position, momentum
                  reso2trkv0data::Sign,
                  reso2trkv0data::PxPos, reso2trkv0data::PyPos, reso2trkv0data::PzPos,
                  reso2trkv0data::PxNeg, reso2trkv0data::PyNeg, reso2trkv0data::PzNeg,
                  reso2trkv0data::PxTrack, reso2trkv0data::PyTrack, reso2trkv0data::PzTrack,
                  reso2trkv0data::XV0, reso2trkv0data::YV0, reso2trkv0data::ZV0,

                  // Saved from finding: DCAs
                  reso2trkv0data::DCAV0Daughters,
                  reso2trkv0data::DCAPosToPV, reso2trkv0data::DCANegToPV, reso2trkv0data::DCATrackToPV,

                  // Dynamic columns
                  reso2trkv0data::Pt<reso2trkv0data::Px, reso2trkv0data::Py>,
                  reso2trkv0data::V0Radius<reso2trkv0data::XV0, reso2trkv0data::YV0>,
                  reso2trkv0data::V0CosPA<reso2trkv0data::XV0, reso2trkv0data::YV0, reso2trkv0data::ZV0, reso2trkv0data::PxV0, reso2trkv0data::PyV0, reso2trkv0data::PzV0>,
                  reso2trkv0data::DCAV0ToPV<reso2trkv0data::XV0, reso2trkv0data::YV0, reso2trkv0data::ZV0, reso2trkv0data::PxV0, reso2trkv0data::PyV0, reso2trkv0data::PzV0>,

                  // V0 Invariant masses
                  reso2trkv0data::MLambda<reso2trkv0data::PxPos, reso2trkv0data::PyPos, reso2trkv0data::PzPos, reso2trkv0data::PxNeg, reso2trkv0data::PyNeg, reso2trkv0data::PzNeg>,
                  reso2trkv0data::MAntiLambda<reso2trkv0data::PxPos, reso2trkv0data::PyPos, reso2trkv0data::PzPos, reso2trkv0data::PxNeg, reso2trkv0data::PyNeg, reso2trkv0data::PzNeg>,
                  reso2trkv0data::MK0Short<reso2trkv0data::PxPos, reso2trkv0data::PyPos, reso2trkv0data::PzPos, reso2trkv0data::PxNeg, reso2trkv0data::PyNeg, reso2trkv0data::PzNeg>,

                  reso2trkv0data::Eta<reso2trkv0data::Px, reso2trkv0data::Py, reso2trkv0data::Pz>);

using Reso2TrackV0DataOrigin = Reso2TrackV0Data;

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(Reso2TrackV0DataExt, Reso2TrackV0DataOrigin, "Reso2TVDATAEXT", //!
                                reso2trkv0data::PxV0, reso2trkv0data::PyV0, reso2trkv0data::PzV0,
                                reso2trkv0data::Px, reso2trkv0data::Py, reso2trkv0data::Pz); // the table name has here to be the one with EXT which is not nice and under study

using Reso2TrackV0DataFull = Reso2TrackV0DataExt;

} // namespace o2::aod
#endif // O2_ANALYSIS_LFRESONANCETABLES_H_