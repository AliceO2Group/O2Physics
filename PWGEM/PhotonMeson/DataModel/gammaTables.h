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

#include <vector>

#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/CaloClusters.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
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
DECLARE_SOA_COLUMN(IsPHOSCPVReadout, isPHOSCPVreadout, bool);
DECLARE_SOA_COLUMN(IsEMCReadout, isEMCreadout, bool);

} // namespace emreducedevent
DECLARE_SOA_TABLE(EMReducedEvents, "AOD", "EMREDUCEDEVENT", //!   Main event information table
                  o2::soa::Index<>, emreducedevent::CollisionId, emreducedevent::Tag, bc::RunNumber, bc::TriggerMask, evsel::Sel8,
                  emreducedevent::IsPHOSCPVReadout, emreducedevent::IsEMCReadout,
                  collision::PosX, collision::PosY, collision::PosZ,
                  collision::NumContrib, collision::CollisionTime, collision::CollisionTimeRes,
                  mult::MultTPC, mult::MultFV0A, mult::MultFV0C, mult::MultFT0A, mult::MultFT0C,
                  mult::MultFDDA, mult::MultFDDC, mult::MultZNA, mult::MultZNC, mult::MultTracklets, mult::MultNTracksPV, mult::MultNTracksPVeta1,
                  cent::CentFT0M, cent::CentFT0A, cent::CentFT0C, cent::CentNTPV,
                  emreducedevent::NgammaPCM, emreducedevent::NgammaPHOS, emreducedevent::NgammaEMC);
using EMReducedEvent = EMReducedEvents::iterator;

namespace emreducedmcevent
{
DECLARE_SOA_COLUMN(PosX, posX, float); //!
DECLARE_SOA_COLUMN(PosY, posY, float); //!
DECLARE_SOA_COLUMN(PosZ, posZ, float); //!
} // namespace emreducedmcevent
DECLARE_SOA_TABLE(EMReducedMCEvents, "AOD", "EMMCEVENT", //!   MC event information table
                  o2::soa::Index<>, mccollision::GeneratorsID,
                  emreducedmcevent::PosX, emreducedmcevent::PosY, emreducedmcevent::PosZ,
                  mccollision::T, mccollision::Weight, mccollision::ImpactParameter);
using EMReducedMCEvent = EMReducedMCEvents::iterator;

namespace emmceventlabel
{
DECLARE_SOA_INDEX_COLUMN(EMReducedMCEvent, emreducedmcevent); //! MC collision
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);                 //! Bit mask to indicate collision mismatches (bit ON means mismatch). Bit 15: indicates negative label
} // namespace emmceventlabel

DECLARE_SOA_TABLE(EMReducedMCEventLabels, "AOD", "EMMCEVENTLABEL", //! Table joined to the EMReducedEvents table containing the MC index
                  emmceventlabel::EMReducedMCEventId, emmceventlabel::McMask);
using EMReducedMCEventLabel = EMReducedMCEventLabels::iterator;

namespace emmcparticle
{
DECLARE_SOA_INDEX_COLUMN(EMReducedMCEvent, emreducedmcevent);                             //!
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(Mother0, mother0, int, "EMMCParticles_Mother0");       //! Track index of the first mother
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(Mother1, mother1, int, "EMMCParticles_Mother1");       //! Track index of the last mother
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(Daughter0, daughter0, int, "EMMCParticles_Daughter0"); //! Track index of the first daughter
DECLARE_SOA_SELF_INDEX_COLUMN_FULL(Daughter1, daughter1, int, "EMMCParticles_Daughter1"); //! Track index of the last daughter
DECLARE_SOA_SELF_ARRAY_INDEX_COLUMN(Mothers, mothers);                                    //! Mother tracks (possible empty) array. Iterate over mcParticle.mothers_as<aod::McParticles>())
DECLARE_SOA_SELF_SLICE_INDEX_COLUMN(Daughters, daughters);                                //! Daughter tracks (possibly empty) slice. Check for non-zero with mcParticle.has_daughters(). Iterate over mcParticle.daughters_as<aod::McParticles>())
DECLARE_SOA_COLUMN(Pt, pt, float);                                                        //!
DECLARE_SOA_COLUMN(Eta, eta, float);                                                      //!
DECLARE_SOA_COLUMN(Phi, phi, float);                                                      //!
DECLARE_SOA_COLUMN(E, e, float);                                                          //!
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, [](float pt, float phi) -> float { return pt * std::cos(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, [](float pt, float phi) -> float { return pt * std::sin(phi); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, [](float pt, float eta) -> float { return pt * std::sinh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float pt, float eta) -> float { return pt * std::cosh(eta); });
DECLARE_SOA_DYNAMIC_COLUMN(Y, y, //! Particle rapidity
                           [](float pt, float eta, float e) -> float {
                             float pz = pt * std::sinh(eta);
                             if ((e - pz) > static_cast<float>(1e-7)) {
                               return 0.5f * std::log((e + pz) / (e - pz));
                             } else {
                               return -999.0f;
                             }
                           });
} // namespace emmcparticle
// NOTE: This table is nearly identical to the one from Framework (except that it points to the event ID, not the BC id)
//       This table contains all MC truth tracks (both v0 and calos)
DECLARE_SOA_TABLE_FULL(EMMCParticles, "EMMCParticles", "AOD", "EMMCPARTICLE", //!  MC track information (on disk)
                       o2::soa::Index<>, emmcparticle::EMReducedMCEventId,
                       mcparticle::PdgCode, mcparticle::StatusCode, mcparticle::Flags,
                       emmcparticle::MothersIds, emmcparticle::DaughtersIdSlice,
                       mcparticle::Weight,
                       emmcparticle::Pt, emmcparticle::Eta, emmcparticle::Phi, emmcparticle::E,
                       mcparticle::Vx, mcparticle::Vy, mcparticle::Vz, mcparticle::Vt,

                       // dynamic column
                       emmcparticle::Px<emmcparticle::Pt, emmcparticle::Phi>,
                       emmcparticle::Py<emmcparticle::Pt, emmcparticle::Phi>,
                       emmcparticle::Pz<emmcparticle::Pt, emmcparticle::Eta>,
                       emmcparticle::P<emmcparticle::Pt, emmcparticle::Eta>,
                       emmcparticle::Y<emmcparticle::Pt, emmcparticle::Eta, emmcparticle::E>,
                       mcparticle::ProducedByGenerator<mcparticle::Flags>,
                       mcparticle::FromBackgroundEvent<mcparticle::Flags>,
                       mcparticle::GetGenStatusCode<mcparticle::Flags, mcparticle::StatusCode>,
                       mcparticle::GetProcess<mcparticle::Flags, mcparticle::StatusCode>,
                       mcparticle::IsPhysicalPrimary<mcparticle::Flags>);

using EMMCParticle = EMMCParticles::iterator;

namespace v0legmclabel
{
DECLARE_SOA_INDEX_COLUMN(EMMCParticle, emmcparticle); //!
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);
} // namespace v0legmclabel

// NOTE: MC labels. This table has one entry for each reconstructed track (joinable with v0leg table)
DECLARE_SOA_TABLE(V0LegMCLabels, "AOD", "V0LEGMCLABEL", //!
                  v0legmclabel::EMMCParticleId, v0legmclabel::McMask);
using V0LegMCLabel = V0LegMCLabels::iterator;

namespace emprimarytrackmclabel
{
DECLARE_SOA_INDEX_COLUMN(EMMCParticle, emmcparticle); //!
DECLARE_SOA_COLUMN(McMask, mcMask, uint16_t);
} // namespace emprimarytrackmclabel

// NOTE: MC labels. This table has one entry for each reconstructed track (joinable with emprimarytrack table)
DECLARE_SOA_TABLE(EMPrimaryTrackMCLabels, "AOD", "EMPRMTRKMCLABEL", //!
                  emprimarytrackmclabel::EMMCParticleId, emprimarytrackmclabel::McMask);
using EMPrimaryTrackMCLabel = EMPrimaryTrackMCLabels::iterator;

namespace v0leg
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);   //!
DECLARE_SOA_INDEX_COLUMN(Track, track);           //!
DECLARE_SOA_COLUMN(Sign, sign, int);              //!
DECLARE_SOA_COLUMN(IsAmbTrack, isAmbTrack, bool); //!

} // namespace v0leg
DECLARE_SOA_TABLE(V0Legs, "AOD", "V0LEG", //!
                  o2::soa::Index<>, v0leg::CollisionId,
                  v0leg::TrackId, v0leg::Sign, v0leg::IsAmbTrack,
                  track::Pt, track::Eta, track::Phi, track::P,
                  track::DcaXY, track::DcaZ,
                  track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows,
                  track::TPCChi2NCl, track::TPCInnerParam,
                  track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaPi,
                  track::ITSClusterMap, track::ITSChi2NCl, track::DetectorMap,
                  track::X, track::Y, track::Z, track::Snp, track::Tgl, track::Alpha, track::Signed1Pt,

                  // dynamic column
                  track::IsWithinBeamPipe<track::X>,
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
DECLARE_SOA_INDEX_COLUMN(EMReducedEvent, emreducedevent);               //!
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
DECLARE_SOA_DYNAMIC_COLUMN(E, e, [](float pxpos, float pxneg, float pypos, float pyneg, float pzpos, float pzneg, float m = 0) -> float { return RecoDecay::sqrtSumOfSquares(pxpos + pxneg, pypos + pyneg, pzpos + pzneg, m); }); //! energy of v0 photn, mass to be given as argument when getter is called!
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float pxpos, float pypos, float pxneg, float pyneg) -> float { return RecoDecay::sqrtSumOfSquares(pxpos + pxneg, pypos + pyneg); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::eta(std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float pxpos, float pypos, float pxneg, float pyneg) -> float { return RecoDecay::phi(pxpos + pxneg, pypos + pyneg); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::sqrtSumOfSquares(pxpos + pxneg, pypos + pyneg, pzpos + pzneg); });
DECLARE_SOA_DYNAMIC_COLUMN(PPos, ppos, [](float px, float py, float pz) -> float { return RecoDecay::sqrtSumOfSquares(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(PNeg, pneg, [](float px, float py, float pz) -> float { return RecoDecay::sqrtSumOfSquares(px, py, pz); });

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
                  v0photon::E<v0photon::PxPosAtSV, v0photon::PxNegAtSV, v0photon::PyPosAtSV, v0photon::PyNegAtSV, v0photon::PzPosAtSV, v0photon::PzNegAtSV>,
                  v0photon::Pt<v0photon::PxPosAtSV, v0photon::PyPosAtSV, v0photon::PxNegAtSV, v0photon::PyNegAtSV>,
                  v0photon::Eta<v0photon::PxPosAtSV, v0photon::PyPosAtSV, v0photon::PzPosAtSV, v0photon::PxNegAtSV, v0photon::PyNegAtSV, v0photon::PzNegAtSV>,
                  v0photon::Phi<v0photon::PxPosAtSV, v0photon::PyPosAtSV, v0photon::PxNegAtSV, v0photon::PyNegAtSV>,
                  v0photon::P<v0photon::PxPosAtSV, v0photon::PyPosAtSV, v0photon::PzPosAtSV, v0photon::PxNegAtSV, v0photon::PyNegAtSV, v0photon::PzNegAtSV>,
                  v0photon::PPos<v0photon::PxPosAtSV, v0photon::PyPosAtSV, v0photon::PzPosAtSV>,
                  v0photon::PNeg<v0photon::PxNegAtSV, v0photon::PyNegAtSV, v0photon::PzNegAtSV>,
                  v0data::V0Radius<v0photon::Vx, v0photon::Vy>,
                  v0data::Alpha<v0photon::PxPosAtSV, v0photon::PyPosAtSV, v0photon::PzPosAtSV, v0photon::PxNegAtSV, v0photon::PyNegAtSV, v0photon::PzNegAtSV>,
                  v0data::QtArm<v0photon::PxPosAtSV, v0photon::PyPosAtSV, v0photon::PzPosAtSV, v0photon::PxNegAtSV, v0photon::PyNegAtSV, v0photon::PzNegAtSV>,
                  v0data::PsiPair<v0photon::PxPosAtSV, v0photon::PyPosAtSV, v0photon::PzPosAtSV, v0photon::PxNegAtSV, v0photon::PyNegAtSV, v0photon::PzNegAtSV>,

                  // invariant mass mee
                  v0data::MGamma<v0photon::PxPosAtSV, v0photon::PyPosAtSV, v0photon::PzPosAtSV, v0photon::PxNegAtSV, v0photon::PyNegAtSV, v0photon::PzNegAtSV>);
// iterators
using V0Photon = V0Photons::iterator;

DECLARE_SOA_TABLE(V0EMReducedEventIds, "AOD", "V0EMEVENTID", v0photon::EMReducedEventId); // To be joined with V0Photons table at analysis level.
// iterators
using V0EMReducedEventId = V0EMReducedEventIds::iterator;

namespace v0photonkf
{
DECLARE_SOA_INDEX_COLUMN(EMReducedEvent, emreducedevent);               //!
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                         //!
DECLARE_SOA_INDEX_COLUMN(V0Photon, v0photon);                           //!
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, V0Legs, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, V0Legs, "_Neg"); //!
DECLARE_SOA_COLUMN(Vx, vx, float);                                      //!
DECLARE_SOA_COLUMN(Vy, vy, float);                                      //!
DECLARE_SOA_COLUMN(Vz, vz, float);                                      //!
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(MGamma, mGamma, float);
DECLARE_SOA_COLUMN(MGammaKFPV, mGammaKFPV, float);
DECLARE_SOA_COLUMN(MGammaKFSV, mGammaKFSV, float);
DECLARE_SOA_COLUMN(CosPA, cospa, float); //!
DECLARE_SOA_COLUMN(PCA, pca, float);     //!

DECLARE_SOA_COLUMN(Alpha, alpha, float);
DECLARE_SOA_COLUMN(QtArm, qtarm, float);
DECLARE_SOA_COLUMN(PsiPair, psipair, float);
DECLARE_SOA_COLUMN(PhiV, phiv, float);
DECLARE_SOA_COLUMN(ChiSquareNDF, chiSquareNDF, float); // Chi2 / NDF of the reconstructed V0

DECLARE_SOA_DYNAMIC_COLUMN(E, e, [](float px, float py, float pz, float m = 0) -> float { return RecoDecay::sqrtSumOfSquares(px, py, pz, m); }); //! energy of v0 photn, mass to be given as argument when getter is called!
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float px, float py) -> float { return RecoDecay::sqrtSumOfSquares(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float px, float py, float pz) -> float { return RecoDecay::eta(std::array{px, py, pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float px, float py) -> float { return RecoDecay::phi(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, [](float px, float py, float pz) -> float { return RecoDecay::sqrtSumOfSquares(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(V0Radius, v0radius, [](float vx, float vy) -> float { return RecoDecay::sqrtSumOfSquares(vx, vy); });

} // namespace v0photonkf
DECLARE_SOA_TABLE(V0PhotonsKF, "AOD", "V0PHOTONKF", //!
                  o2::soa::Index<>, v0photonkf::CollisionId, v0photonkf::V0PhotonId, v0photon::PosTrackId, v0photon::NegTrackId,
                  v0photonkf::Vx, v0photonkf::Vy, v0photonkf::Vz,
                  v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz,
                  v0photonkf::MGamma, v0photonkf::MGammaKFPV, v0photonkf::MGammaKFSV,
                  v0photonkf::CosPA, v0photonkf::PCA,
                  v0photonkf::Alpha, v0photonkf::QtArm, v0photonkf::PsiPair, v0photonkf::PhiV,
                  v0photonkf::ChiSquareNDF,
                  // dynamic column
                  v0photonkf::E<v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz>,
                  v0photonkf::Pt<v0photonkf::Px, v0photonkf::Py>,
                  v0photonkf::Eta<v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz>,
                  v0photonkf::Phi<v0photonkf::Px, v0photonkf::Py>,
                  v0photonkf::P<v0photonkf::Px, v0photonkf::Py, v0photonkf::Pz>,
                  v0photonkf::V0Radius<v0photonkf::Vx, v0photonkf::Vy>);
// iterators
using V0PhotonKF = V0PhotonsKF::iterator;

DECLARE_SOA_TABLE(V0KFEMReducedEventIds, "AOD", "V0KFEMEVENTID", v0photonkf::EMReducedEventId); // To be joined with V0PhotonsKF table at analysis level.
// iterators
using V0KFEMReducedEventId = V0KFEMReducedEventIds::iterator;

namespace emprimarytrack
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int); //!
DECLARE_SOA_COLUMN(TrackId, trackId, int);         //!
DECLARE_SOA_COLUMN(Sign, sign, int);               //!
} // namespace emprimarytrack
DECLARE_SOA_TABLE(EMPrimaryTracks, "AOD", "EMPRIMARYTRACK", //!
                  o2::soa::Index<>, emprimarytrack::CollisionId,
                  emprimarytrack::TrackId, emprimarytrack::Sign,
                  track::Pt, track::Eta, track::Phi, track::DcaXY, track::DcaZ,
                  track::TPCNClsFindable, track::TPCNClsFindableMinusFound, track::TPCNClsFindableMinusCrossedRows,
                  track::TPCChi2NCl, track::TPCInnerParam,
                  track::TPCSignal, pidtpc::TPCNSigmaEl, pidtpc::TPCNSigmaMu, pidtpc::TPCNSigmaPi, pidtpc::TPCNSigmaKa, pidtpc::TPCNSigmaPr,
                  pidtofbeta::Beta, pidtof::TOFNSigmaEl, pidtof::TOFNSigmaMu, pidtof::TOFNSigmaPi, pidtof::TOFNSigmaKa, pidtof::TOFNSigmaPr,
                  track::ITSClusterMap, track::ITSChi2NCl, track::DetectorMap, track::Signed1Pt, track::CYY, track::CZZ,

                  // dynamic column
                  track::TPCNClsFound<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::TPCNClsCrossedRows<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCCrossedRowsOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusCrossedRows>,
                  track::TPCFoundOverFindableCls<track::TPCNClsFindable, track::TPCNClsFindableMinusFound>,
                  track::ITSNCls<track::ITSClusterMap>,
                  track::HasITS<track::DetectorMap>, track::HasTPC<track::DetectorMap>,
                  track::HasTRD<track::DetectorMap>, track::HasTOF<track::DetectorMap>);

// iterators
using EMPrimaryTrack = EMPrimaryTracks::iterator;

namespace dalitzee
{
DECLARE_SOA_INDEX_COLUMN(EMReducedEvent, emreducedevent);                        //!
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, EMPrimaryTracks, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, EMPrimaryTracks, "_Neg"); //!
DECLARE_SOA_COLUMN(CollisionId, collisionId, int);                               //!
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Mee, mee, float);
DECLARE_SOA_COLUMN(PhiV, phiv, float);
DECLARE_SOA_COLUMN(DCAeeXY, dcaeeXY, float);
DECLARE_SOA_COLUMN(DCAeeZ, dcaeeZ, float);
DECLARE_SOA_COLUMN(Sign, sign, int);                                                                                                     //!
DECLARE_SOA_DYNAMIC_COLUMN(Energy, e, [](float pt, float eta, float m) { return RecoDecay::sqrtSumOfSquares(pt * std::cosh(eta), m); }); // e = sqrt(p*p + m*m)
} // namespace dalitzee
DECLARE_SOA_TABLE(DalitzEEs, "AOD", "DALITZEE", //!
                  o2::soa::Index<>, dalitzee::CollisionId, dalitzee::PosTrackId, dalitzee::NegTrackId,
                  dalitzee::Pt, dalitzee::Eta, dalitzee::Phi, dalitzee::Mee, dalitzee::PhiV, dalitzee::DCAeeXY, dalitzee::DCAeeZ, dalitzee::Sign,
                  dalitzee::Energy<o2::aod::dalitzee::Pt, o2::aod::dalitzee::Eta, o2::aod::dalitzee::Mee>);
// iterators
using DalitzEE = DalitzEEs::iterator;

DECLARE_SOA_TABLE(DalitzEEEMReducedEventIds, "AOD", "EEEMEVENTID", dalitzee::EMReducedEventId); // To be joined with DalitzEEs table at analysis level.
// iterators
using DalitzEEEMReducedEventId = DalitzEEEMReducedEventIds::iterator;

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

namespace v0Recalculations
{
DECLARE_SOA_COLUMN(RecalculatedVtxX, recalculatedVtxX, float); //! Recalculated conversion point
DECLARE_SOA_COLUMN(RecalculatedVtxY, recalculatedVtxY, float); //! Recalculated conversion point
DECLARE_SOA_COLUMN(RecalculatedVtxZ, recalculatedVtxZ, float); //! Recalculated conversion point
DECLARE_SOA_DYNAMIC_COLUMN(RecalculatedVtxR, recalculatedVtxR, [](float x, float y) { return sqrt(x * x + y * y); });
} // namespace v0Recalculations

DECLARE_SOA_TABLE(V0Recalculation, "AOD", "V0RECALC",
                  v0Recalculations::RecalculatedVtxX,
                  v0Recalculations::RecalculatedVtxY,
                  v0Recalculations::RecalculatedVtxZ,
                  v0Recalculations::RecalculatedVtxR<o2::aod::v0Recalculations::RecalculatedVtxX, o2::aod::v0Recalculations::RecalculatedVtxY>);

namespace gammamctrue
{
DECLARE_SOA_COLUMN(P, p, float); //! Absolute momentum in GeV/c
} // namespace gammamctrue

DECLARE_SOA_TABLE(McDaughterTrue, "AOD", "MCDAUTRUE",
                  gammamctrue::P);

namespace gammamctrue
{
DECLARE_SOA_INDEX_COLUMN_FULL(V0Photon, v0photon, int, V0Photons, "_V0Photon"); //!
DECLARE_SOA_COLUMN(Gamma, gamma, int64_t);                                      //! Used as reference for the daughters
DECLARE_SOA_COLUMN(NDaughters, nDaughters, int);                                //! Number of daughters
DECLARE_SOA_COLUMN(Eta, eta, float);                                            //! Pseudorapidity
DECLARE_SOA_COLUMN(Phi, phi, float);                                            //! Angle phi in rad
DECLARE_SOA_COLUMN(Pt, pt, float);                                              //! Transversal momentum in GeV/c
DECLARE_SOA_COLUMN(Y, y, float);                                                //! Rapidity

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
                  // v0data::V0Id, // reference to reconstructed v0 (if its a task with reconstucted info)
                  gammamctrue::V0PhotonId, // reference to reconstructed v0 (if its a task with reconstucted info)
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

namespace skimmedcluster
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                        //! collisionID used as index for matched clusters
DECLARE_SOA_INDEX_COLUMN(BC, bc);                                      //! bunch crossing ID used as index for ambiguous clusters
DECLARE_SOA_COLUMN(ID, id, int);                                       //! cluster ID identifying cluster in event
DECLARE_SOA_COLUMN(E, e, float);                                       //! cluster energy (GeV)
DECLARE_SOA_COLUMN(Eta, eta, float);                                   //! cluster pseudorapidity (calculated using vertex)
DECLARE_SOA_COLUMN(Phi, phi, float);                                   //! cluster azimuthal angle (calculated using vertex)
DECLARE_SOA_COLUMN(M02, m02, float);                                   //! shower shape long axis
DECLARE_SOA_COLUMN(M20, m20, float);                                   //! shower shape short axis
DECLARE_SOA_COLUMN(NCells, nCells, int);                               //! number of cells in cluster
DECLARE_SOA_COLUMN(Time, time, float);                                 //! cluster time (ns)
DECLARE_SOA_COLUMN(DistanceToBadChannel, distanceToBadChannel, float); //! distance to bad channel
DECLARE_SOA_COLUMN(NLM, nlm, int);                                     //! number of local maxima
} // namespace skimmedcluster

namespace emccluster
{
DECLARE_SOA_INDEX_COLUMN(EMReducedEvent, emreducedevent);                                                                     //!
DECLARE_SOA_COLUMN(CoreEnergy, coreEnergy, float);                                                                            //! cluster core energy (GeV)
DECLARE_SOA_COLUMN(Time, time, float);                                                                                        //! cluster time (ns)
DECLARE_SOA_COLUMN(IsExotic, isExotic, bool);                                                                                 //! flag to mark cluster as exotic
DECLARE_SOA_COLUMN(Definition, definition, int);                                                                              //! cluster definition, see EMCALClusterDefinition.h
DECLARE_SOA_ARRAY_INDEX_COLUMN(Track, track);                                                                                 //! TrackIds
DECLARE_SOA_COLUMN(TrackEta, tracketa, std::vector<float>);                                                                   //! eta values of the matched tracks
DECLARE_SOA_COLUMN(TrackPhi, trackphi, std::vector<float>);                                                                   //! phi values of the matched tracks
DECLARE_SOA_COLUMN(TrackP, trackp, std::vector<float>);                                                                       //! momentum values of the matched tracks
DECLARE_SOA_COLUMN(TrackPt, trackpt, std::vector<float>);                                                                     //! pt values of the matched tracks
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float e, float eta, float m = 0) -> float { return sqrt(e * e - m * m) / cosh(eta); }); //! cluster pt, mass to be given as argument when getter is called!
} // namespace emccluster
DECLARE_SOA_TABLE(SkimEMCClusters, "AOD", "SKIMEMCCLUSTERS", //! table of skimmed EMCal clusters
                  o2::soa::Index<>, skimmedcluster::CollisionId, skimmedcluster::BCId, skimmedcluster::E, emccluster::CoreEnergy,
                  skimmedcluster::Eta, skimmedcluster::Phi, skimmedcluster::M02, skimmedcluster::M20, skimmedcluster::NCells, skimmedcluster::Time,
                  emccluster::IsExotic, skimmedcluster::DistanceToBadChannel, skimmedcluster::NLM, emccluster::Definition,
                  emccluster::TrackIds, emccluster::TrackEta, emccluster::TrackPhi, emccluster::TrackP, emccluster::TrackPt,
                  // dynamic column
                  emccluster::Pt<skimmedcluster::E, skimmedcluster::Eta>);
using SkimEMCCluster = SkimEMCClusters::iterator;

DECLARE_SOA_TABLE(EMCEMReducedEventIds, "AOD", "EMCEMEVENTID", emccluster::EMReducedEventId); // To be joined with SkimEMCClusters table at analysis level.
// iterators
using EMCEMReducedEventId = EMCEMReducedEventIds::iterator;

namespace phoscluster
{
DECLARE_SOA_INDEX_COLUMN(EMReducedEvent, emreducedevent);                           //!
DECLARE_SOA_INDEX_COLUMN_FULL(MatchedTrack, matchedTrack, int, Tracks, "_Matched"); //! matched track index
DECLARE_SOA_COLUMN(X, x, float);                                                    //! cluster hit position in ALICE global coordinate
DECLARE_SOA_COLUMN(Y, y, float);                                                    //! cluster hit position in ALICE global coordinate
DECLARE_SOA_COLUMN(Z, z, float);                                                    //! cluster hit position in ALICE global coordinate
DECLARE_SOA_COLUMN(CellX, cellx, int);                                              //! cell index x of cluster hit position
DECLARE_SOA_COLUMN(CellZ, cellz, int);                                              //! cell index z of cluster hit position
// DECLARE_SOA_COLUMN(TrackEta, tracketa, float);                                      //! eta of the matched track
// DECLARE_SOA_COLUMN(TrackPhi, trackphi, float);                                      //! phi of the matched track
// DECLARE_SOA_COLUMN(TrackP, trackp, float);                                          //! momentum of the matched track
// DECLARE_SOA_COLUMN(TrackPt, trackpt, float);                                        //! pt of the matched track
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, [](float e, float x, float y, float z, float m = 0) -> float { return x / RecoDecay::sqrtSumOfSquares(x, y, z) * sqrt(e * e - m * m); });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, [](float e, float x, float y, float z, float m = 0) -> float { return y / RecoDecay::sqrtSumOfSquares(x, y, z) * sqrt(e * e - m * m); });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, [](float e, float x, float y, float z, float m = 0) -> float { return z / RecoDecay::sqrtSumOfSquares(x, y, z) * sqrt(e * e - m * m); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, [](float e, float x, float y, float z, float m = 0) -> float { return RecoDecay::sqrtSumOfSquares(x, y) / RecoDecay::sqrtSumOfSquares(x, y, z) * sqrt(e * e - m * m); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, [](float x, float y, float z) -> float { return RecoDecay::eta(std::array{x, y, z}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, [](float x, float y) -> float { return RecoDecay::phi(x, y); });
} // namespace phoscluster

DECLARE_SOA_TABLE(PHOSClusters, "AOD", "PHOSCLUSTERS", //!
                  o2::soa::Index<>, skimmedcluster::CollisionId, phoscluster::MatchedTrackId,
                  skimmedcluster::E, phoscluster::X, phoscluster::Y, phoscluster::Z,
                  skimmedcluster::M02, skimmedcluster::M20, skimmedcluster::NCells,
                  skimmedcluster::Time, skimmedcluster::DistanceToBadChannel, skimmedcluster::NLM,
                  calocluster::Module, phoscluster::CellX, phoscluster::CellZ,
                  // phoscluster::TrackEta, phoscluster::TrackPhi, phoscluster::TrackP, phoscluster::TrackPt,
                  // dynamic column
                  phoscluster::Px<skimmedcluster::E, phoscluster::X, phoscluster::Y, phoscluster::Z>,
                  phoscluster::Py<skimmedcluster::E, phoscluster::X, phoscluster::Y, phoscluster::Z>,
                  phoscluster::Pz<skimmedcluster::E, phoscluster::X, phoscluster::Y, phoscluster::Z>,
                  phoscluster::Pt<skimmedcluster::E, phoscluster::X, phoscluster::Y, phoscluster::Z>,
                  phoscluster::Eta<phoscluster::X, phoscluster::Y, phoscluster::Z>,
                  phoscluster::Phi<phoscluster::X, phoscluster::Y>);
using PHOSCluster = PHOSClusters::iterator;

DECLARE_SOA_TABLE(PHOSEMReducedEventIds, "AOD", "PHOSEMEVENTID", phoscluster::EMReducedEventId); // To be joined with PHOSClusters table at analysis level.
// iterators
using PHOSEMReducedEventId = PHOSEMReducedEventIds::iterator;

namespace caloextra
{
DECLARE_SOA_INDEX_COLUMN_FULL(Cluster, cluster, int, SkimEMCClusters, ""); //! reference to the gamma in the skimmed EMCal table
DECLARE_SOA_INDEX_COLUMN_FULL(Cell, cell, int, Calos, "");                 //! reference to the gamma in the skimmed EMCal table
// DECLARE_SOA_INDEX_COLUMN(Track, track);                   //! TrackID
DECLARE_SOA_COLUMN(TrackEta, tracketa, float); //! eta of the matched track
DECLARE_SOA_COLUMN(TrackPhi, trackphi, float); //! phi of the matched track
DECLARE_SOA_COLUMN(TrackP, trackp, float);     //! momentum of the matched track
DECLARE_SOA_COLUMN(TrackPt, trackpt, float);   //! pt of the matched track
} // namespace caloextra

DECLARE_SOA_TABLE(SkimEMCCells, "AOD", "SKIMEMCCELLS",                        //! table of link between skimmed EMCal clusters and their cells
                  o2::soa::Index<>, caloextra::ClusterId, caloextra::CellId); //!
using SkimEMCCell = SkimEMCCells::iterator;

DECLARE_SOA_TABLE(SkimEMCMTs, "AOD", "SKIMEMCMTS", //! table of link between skimmed EMCal clusters and their matched tracks
                  o2::soa::Index<>, caloextra::ClusterId, caloextra::TrackEta,
                  caloextra::TrackPhi, caloextra::TrackP, caloextra::TrackPt);
using SkimEMCMT = SkimEMCMTs::iterator;

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
                  o2::soa::Index<>, skimmedcluster::CollisionId, gammareco::Method,
                  skimmedcluster::E, skimmedcluster::Eta, skimmedcluster::Phi,
                  gammareco::SkimmedEMCId, gammareco::SkimmedPHOSId);
DECLARE_SOA_TABLE(SkimPCMCuts, "AOD", "SKIMPCMCUTS",                                  //! table of link between skimmed PCM photon candidates and their cuts
                  o2::soa::Index<>, gammareco::SkimmedPCMId, gammareco::PCMCutBit);   //!
DECLARE_SOA_TABLE(SkimPHOSCuts, "AOD", "SKIMPHOSCUTS",                                //! table of link between skimmed PHOS photon candidates and their cuts
                  o2::soa::Index<>, gammareco::SkimmedPHOSId, gammareco::PHOSCutBit); //!
DECLARE_SOA_TABLE(SkimEMCCuts, "AOD", "SKIMEMCCUTS",                                  //! table of link between skimmed EMCal photon candidates and their cuts
                  o2::soa::Index<>, gammareco::SkimmedEMCId, gammareco::EMCCutBit);   //!
} // namespace o2::aod

#endif // PWGEM_PHOTONMESON_DATAMODEL_GAMMATABLES_H_
