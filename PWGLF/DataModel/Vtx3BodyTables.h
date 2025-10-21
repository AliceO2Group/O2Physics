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

/// \file Vtx3BodyTables.h
/// \brief Definitions of analysis tables for 3body decayed hypertriton
/// \author Yuanzhe Wang <yuanzhe.wang@cern.ch>
/// \author Carolina Reetz <c.reetz@cern.ch>

#ifndef PWGLF_DATAMODEL_VTX3BODYTABLES_H_
#define PWGLF_DATAMODEL_VTX3BODYTABLES_H_

#include "Common/Core/RecoDecay.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisDataModel.h"

#include <cmath>

namespace o2::aod
{
namespace vtx3body
{
// indices
DECLARE_SOA_INDEX_COLUMN_FULL(TrackPr, trackPr, int, Tracks, "_pr"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(TrackPi, trackPi, int, Tracks, "_pi"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(TrackDe, trackDe, int, Tracks, "_de"); //!
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                      //!
DECLARE_SOA_INDEX_COLUMN(Decay3Body, decay3body);                    //!

// General 3 body Vtx properties
DECLARE_SOA_COLUMN(Mass, mass, float); //! candidate mass (with H3L or Anti-H3L mass hypothesis depending on deuteron charge)
DECLARE_SOA_COLUMN(Sign, sign, float); //! candidate sign
DECLARE_SOA_COLUMN(X, x, float);       //! decay position X
DECLARE_SOA_COLUMN(Y, y, float);       //! decay position Y
DECLARE_SOA_COLUMN(Z, z, float);       //! decay position Z
DECLARE_SOA_COLUMN(Px, px, float);     //! momentum X
DECLARE_SOA_COLUMN(Py, py, float);     //! momentum Y
DECLARE_SOA_COLUMN(Pz, pz, float);     //! momentum Z
DECLARE_SOA_COLUMN(Chi2, chi2, float); //! KFParticle: chi2geo/ndf or chi2topo/ndf of vertex fit, DCA fitter: Chi2AtPCACandidate value

// daughter properties
DECLARE_SOA_COLUMN(MassV0, massV0, float);       //! V0 mass (with H3L or Anti-H3L mass hypothesis depending on deuteron charge)
DECLARE_SOA_COLUMN(PxTrackPr, pxTrackPr, float); //! track0 px at min
DECLARE_SOA_COLUMN(PyTrackPr, pyTrackPr, float); //! track0 py at min
DECLARE_SOA_COLUMN(PzTrackPr, pzTrackPr, float); //! track0 pz at min
DECLARE_SOA_COLUMN(PxTrackPi, pxTrackPi, float); //! track1 px at min
DECLARE_SOA_COLUMN(PyTrackPi, pyTrackPi, float); //! track1 py at min
DECLARE_SOA_COLUMN(PzTrackPi, pzTrackPi, float); //! track1 pz at min
DECLARE_SOA_COLUMN(PxTrackDe, pxTrackDe, float); //! track2 px at min
DECLARE_SOA_COLUMN(PyTrackDe, pyTrackDe, float); //! track2 py at min
DECLARE_SOA_COLUMN(PzTrackDe, pzTrackDe, float); //! track2 pz at min

// DCAs to PV
DECLARE_SOA_COLUMN(DCAXYTrackPrToPV, dcaXYtrackPrToPv, float);         //! DCAXY of proton to PV (computed with KFParticle)
DECLARE_SOA_COLUMN(DCAXYTrackPiToPV, dcaXYtrackPiToPv, float);         //! DCAXY of pion to PV (computed with KFParticle)
DECLARE_SOA_COLUMN(DCAXYTrackDeToPV, dcaXYtrackDeToPv, float);         //! DCAXY of deuteron to PV (computed with KFParticle)
DECLARE_SOA_COLUMN(DCATrackPrToPV, dcaTrackPrToPv, float);             //! DCA of proton to PV (computed with KFParticle)
DECLARE_SOA_COLUMN(DCATrackPiToPV, dcaTrackPiToPv, float);             //! DCA of pion to PV (computed with KFParticle)
DECLARE_SOA_COLUMN(DCATrackDeToPV, dcaTrackDeToPv, float);             //! DCA of deuteron to PV (computed with KFParticle)
DECLARE_SOA_COLUMN(DCAXYTrackPrToPVProp, dcaXYtrackPrToPvProp, float); //! DCAXY of proton to PV (propagated with O2 Propagator)
DECLARE_SOA_COLUMN(DCAXYTrackPiToPVProp, dcaXYtrackPiToPvProp, float); //! DCAXY of pion to PV (propagated with O2 Propagator)
DECLARE_SOA_COLUMN(DCAXYTrackDeToPVProp, dcaXYtrackDeToPvProp, float); //! DCAXY of deuteron to PV (propagated with O2 Propagator)
DECLARE_SOA_COLUMN(DCATrackPrToPVProp, dcaTrackPrToPvProp, float);     //! DCA of proton to PV (propagated with O2 Propagator)
DECLARE_SOA_COLUMN(DCATrackPiToPVProp, dcaTrackPiToPvProp, float);     //! DCA of pion to PV (propagated with O2 Propagator)
DECLARE_SOA_COLUMN(DCATrackDeToPVProp, dcaTrackDeToPvProp, float);     //! DCA of deuteron to PV (propagated with O2 Propagator)

// DCAs to SV
DECLARE_SOA_COLUMN(DCATrackPrToSV, dcaTrackPrToSv, float);           //! DCA of proton to SV
DECLARE_SOA_COLUMN(DCATrackPiToSV, dcaTrackPiToSv, float);           //! DCA of pion to SV
DECLARE_SOA_COLUMN(DCATrackDeToSV, dcaTrackDeToSv, float);           //! DCA of deuteron to SV
DECLARE_SOA_COLUMN(DCAVtxToDaughtersAv, dcaVtxToDaughtersAv, float); //! Quadratic sum of DCA between daughters at SV

// CosPA
DECLARE_SOA_COLUMN(CosPA, cosPA, float); //! Cosine of pointing angle of the 3body candidate

// Ct
DECLARE_SOA_COLUMN(Ct, ct, float); //! Reconstruction Ct of 3body candidate

// Strangeness tracking
DECLARE_SOA_COLUMN(TrackedClSize, trackedClSize, float); //! Average ITS cluster size of strangeness tracked 3body

// PID
DECLARE_SOA_COLUMN(TPCNSigmaPr, tpcNSigmaPr, float);         //! nsigma proton of TPC PID of the proton daughter
DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcNSigmaPi, float);         //! nsigma pion of TPC PID of the pion daughter
DECLARE_SOA_COLUMN(TPCNSigmaDe, tpcNSigmaDe, float);         //! nsigma deuteron of TPC PID of the bachelor daughter
DECLARE_SOA_COLUMN(TPCNSigmaPiBach, tpcNSigmaPiBach, float); //! nsigma pion of TPC PID of the bachelor daughter
DECLARE_SOA_COLUMN(TOFNSigmaDe, tofNSigmaDe, float);         //! nsigma deuteron of TOF PID of the bachelor daughter
DECLARE_SOA_COLUMN(PIDTrackingDe, pidTrackingDe, uint32_t);  //! PID during tracking of bachelor daughter

// Daughter track quality
DECLARE_SOA_COLUMN(TPCNClTrackPr, tpcNClTrackPr, int); //! Number of TPC clusters of proton daughter
DECLARE_SOA_COLUMN(TPCNClTrackPi, tpcNClTrackPi, int); //! Number of TPC clusters of pion daughter
DECLARE_SOA_COLUMN(TPCNClTrackDe, tpcNClTrackDe, int); //! Number of TPC clusters of deuteron daughter
DECLARE_SOA_COLUMN(ITSClSizePr, itsClsizePr, double);  //! average ITS cluster size of proton daughter
DECLARE_SOA_COLUMN(ITSClSizePi, itsClsizePi, double);  //! average ITS cluster size of pion daughter
DECLARE_SOA_COLUMN(ITSClSizeDe, itsClsizeDe, double);  //! average ITS cluster size of deuteron daughter

// Covariance matrices
DECLARE_SOA_COLUMN(CovProton, covProton, float[21]);     //! covariance matrix elements of proton daughter track
DECLARE_SOA_COLUMN(CovPion, covPion, float[21]);         //! covariance matrix elements of pion daughter track
DECLARE_SOA_COLUMN(CovDeuteron, covDeuteron, float[21]); //! covariance matrix elements of deuteron daughter track
DECLARE_SOA_COLUMN(VtxCovMat, vtxCovMat, float[21]);     //! covariance matrix elements of candidate

// Monte Carlo info
DECLARE_SOA_COLUMN(GenPx, genPx, float);                // generated Px of the hypertriton in GeV/c
DECLARE_SOA_COLUMN(GenPy, genPy, float);                // generated Py of the hypertriton in GeV/c
DECLARE_SOA_COLUMN(GenPz, genPz, float);                // generated Pz of the hypertriton in GeV/c
DECLARE_SOA_COLUMN(GenX, genX, float);                  // generated decay vtx position X of the hypertriton
DECLARE_SOA_COLUMN(GenY, genY, float);                  // generated decay vtx position Y of the hypertriton
DECLARE_SOA_COLUMN(GenZ, genZ, float);                  // generated decay vtx position Z of the hypertriton
DECLARE_SOA_COLUMN(GenCt, genCt, float);                // generated Ct of the hypertriton
DECLARE_SOA_COLUMN(GenPhi, genPhi, float);              // generated Phi of the hypertriton
DECLARE_SOA_COLUMN(GenEta, genEta, float);              // Eta of the hypertriton
DECLARE_SOA_COLUMN(GenRap, genRap, float);              // generated rapidity of the hypertriton
DECLARE_SOA_COLUMN(GenPPr, genPPr, float);              //! generated momentum proton daughter particle
DECLARE_SOA_COLUMN(GenPPi, genPPi, float);              //! generated momentum pion daughter particle
DECLARE_SOA_COLUMN(GenPDe, genPDe, float);              //! generated momentum deuteron daughter particle
DECLARE_SOA_COLUMN(GenPtPr, genPtPr, float);            //! generated transverse momentum proton daughter particle
DECLARE_SOA_COLUMN(GenPtPi, genPtPi, float);            //! generated transverse momentum pion daughter particle
DECLARE_SOA_COLUMN(GenPtDe, genPtDe, float);            //! generated transverse momentum deuteron daughter particle
DECLARE_SOA_COLUMN(IsTrueH3L, isTrueH3l, bool);         //! flag for true hypertriton candidate
DECLARE_SOA_COLUMN(IsTrueAntiH3L, isTrueAntiH3l, bool); //! flag for true anti-hypertriton candidate
DECLARE_SOA_COLUMN(MotherPdgCode, motherPdgCode, int);  //! PDG code of the mother particle
DECLARE_SOA_COLUMN(PrPdgCode, prPdgCode, int);          //! MC particle proton PDG code
DECLARE_SOA_COLUMN(PiPdgCode, piPdgCode, int);          //! MC particle pion PDG code
DECLARE_SOA_COLUMN(DePdgCode, dePdgCode, int);          //! MC particle deuteron PDG code
DECLARE_SOA_COLUMN(IsDePrimary, isDePrimary, bool);     //! flag for deuteron daughter primary
DECLARE_SOA_COLUMN(IsSurvEvSel, isSurvEvSel, int);      //! flag if reco collision survived event selection
DECLARE_SOA_COLUMN(IsReco, isreco, int);                //! flag if candidate was reconstructed

// Derived expressions
// Momenta
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //! 3 body pT in GeV/c
                           [](float px, float py) -> float { return RecoDecay::sqrtSumOfSquares(px, py); });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //! 3 body total momentum in GeV/c
                           [](float px, float py, float pz) -> float { return RecoDecay::sqrtSumOfSquares(px, py, pz); });
DECLARE_SOA_DYNAMIC_COLUMN(GenPt, genPt, //! 3 body pT in GeV/c
                           [](float genPx, float genPy) -> float { return RecoDecay::sqrtSumOfSquares(genPx, genPy); });
DECLARE_SOA_DYNAMIC_COLUMN(GenP, genP, //! 3 body total momentum in GeV/c
                           [](float genPx, float genPy, float genPz) -> float { return RecoDecay::sqrtSumOfSquares(genPx, genPy, genPz); });

// Length quantities
DECLARE_SOA_DYNAMIC_COLUMN(VtxRadius, vtxradius, //! 3 body decay radius (2D, centered at zero)
                           [](float x, float y) -> float { return RecoDecay::sqrtSumOfSquares(x, y); });
DECLARE_SOA_DYNAMIC_COLUMN(GenRadius, genRadius, //! 3 body decay radius (2D, centered at zero)
                           [](float genX, float genY) -> float { return RecoDecay::sqrtSumOfSquares(genX, genY); });

// Distance Over To Mom
DECLARE_SOA_DYNAMIC_COLUMN(DistOverTotMom, distovertotmom, //! PV to 3 body decay distance over total momentum
                           [](float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ) {
                             float P = RecoDecay::sqrtSumOfSquares(Px, Py, Pz);
                             return std::sqrt(std::pow(X - pvX, 2) + std::pow(Y - pvY, 2) + std::pow(Z - pvZ, 2)) / (P + 1E-10);
                           });

// Dca to PV
DECLARE_SOA_DYNAMIC_COLUMN(DCAVtxToPV, dcavtxtopv, //! DCA of 3 body vtx to PV
                           [](float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ) -> float { return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz)); });

// Eta and Phi
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, //! 3 body vtx eta
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::eta(std::array{Px, Py, Pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, //! 3 body vtx phi
                           [](float Px, float Py) -> float { return RecoDecay::phi(Px, Py); });

// Rapidity
DECLARE_SOA_DYNAMIC_COLUMN(Rap, rap, //! 3 body vtx y with hypertriton or antihypertriton hypothesis
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::y(std::array{Px, Py, Pz}, o2::constants::physics::MassHyperTriton); });

// Kinematic information of daughter tracks
DECLARE_SOA_DYNAMIC_COLUMN(TrackPrPt, trackPrPt, //! daughter0 pT
                           [](float pxTrackPr, float pyTrackPr) -> float { return RecoDecay::sqrtSumOfSquares(pxTrackPr, pyTrackPr); });
DECLARE_SOA_DYNAMIC_COLUMN(TrackPiPt, trackPiPt, //! daughter1 pT
                           [](float pxTrackPi, float pyTrackPi) -> float { return RecoDecay::sqrtSumOfSquares(pxTrackPi, pyTrackPi); });
DECLARE_SOA_DYNAMIC_COLUMN(TrackDePt, trackDePt, //! daughter2 pT
                           [](float pxTrackDe, float pyTrackDe) -> float { return RecoDecay::sqrtSumOfSquares(pxTrackDe, pyTrackDe); });
DECLARE_SOA_DYNAMIC_COLUMN(TrackPrEta, trackPrEta, //! daughter0 eta
                           [](float pxTrackPr, float pyTrackPr, float pzTrackPr) -> float { return RecoDecay::eta(std::array{pxTrackPr, pyTrackPr, pzTrackPr}); });
DECLARE_SOA_DYNAMIC_COLUMN(TrackPrPhi, trackPrPhi, //! daughter0 phi
                           [](float pxTrackPr, float pyTrackPr) -> float { return RecoDecay::phi(pxTrackPr, pyTrackPr); });
DECLARE_SOA_DYNAMIC_COLUMN(TrackPiEta, trackPiEta, //! daughter1 eta
                           [](float pxTrackPi, float pyTrackPi, float pzTrackPi) -> float { return RecoDecay::eta(std::array{pxTrackPi, pyTrackPi, pzTrackPi}); });
DECLARE_SOA_DYNAMIC_COLUMN(TrackPiPhi, trackPiPhi, //! daughter1 phi
                           [](float pxTrackPi, float pyTrackPi) -> float { return RecoDecay::phi(pxTrackPi, pyTrackPi); });
DECLARE_SOA_DYNAMIC_COLUMN(TrackDeEta, trackDeEta, //! daughter2 eta
                           [](float pxTrackDe, float pyTrackDe, float pzTrackDe) -> float { return RecoDecay::eta(std::array{pxTrackDe, pyTrackDe, pzTrackDe}); });
DECLARE_SOA_DYNAMIC_COLUMN(TrackDePhi, trackDePhi, //! daughter2 phi
                           [](float pxTrackDe, float pyTrackDe) -> float { return RecoDecay::phi(pxTrackDe, pyTrackDe); });
} // namespace vtx3body

// index table
DECLARE_SOA_TABLE(Decay3BodyIndices, "AOD", "3BodyINDEX", //!
                  o2::soa::Index<>,
                  vtx3body::Decay3BodyId,
                  vtx3body::TrackPrId, vtx3body::TrackPiId, vtx3body::TrackDeId,
                  vtx3body::CollisionId);

// reconstructed candidate table for analysis
DECLARE_SOA_TABLE(Vtx3BodyDatas, "AOD", "VTX3BODYDATA", //!
                  o2::soa::Index<>,
                  vtx3body::Sign,
                  vtx3body::Mass, vtx3body::MassV0,
                  vtx3body::X, vtx3body::Y, vtx3body::Z,
                  vtx3body::Px, vtx3body::Py, vtx3body::Pz,
                  vtx3body::Chi2,
                  vtx3body::TrackedClSize,
                  vtx3body::PxTrackPr, vtx3body::PyTrackPr, vtx3body::PzTrackPr,
                  vtx3body::PxTrackPi, vtx3body::PyTrackPi, vtx3body::PzTrackPi,
                  vtx3body::PxTrackDe, vtx3body::PyTrackDe, vtx3body::PzTrackDe,
                  vtx3body::DCAXYTrackPrToPV, vtx3body::DCAXYTrackPiToPV, vtx3body::DCAXYTrackDeToPV,
                  vtx3body::DCATrackPrToPV, vtx3body::DCATrackPiToPV, vtx3body::DCATrackDeToPV,
                  vtx3body::DCAXYTrackPrToPVProp, vtx3body::DCAXYTrackPiToPVProp, vtx3body::DCAXYTrackDeToPVProp,
                  vtx3body::DCATrackPrToPVProp, vtx3body::DCATrackPiToPVProp, vtx3body::DCATrackDeToPVProp,
                  vtx3body::DCATrackPrToSV, vtx3body::DCATrackPiToSV, vtx3body::DCATrackDeToSV,
                  vtx3body::DCAVtxToDaughtersAv,
                  vtx3body::CosPA, vtx3body::Ct,
                  vtx3body::TPCNSigmaPr, vtx3body::TPCNSigmaPi, vtx3body::TPCNSigmaDe, vtx3body::TPCNSigmaPiBach,
                  vtx3body::TOFNSigmaDe,
                  vtx3body::ITSClSizePr, vtx3body::ITSClSizePi, vtx3body::ITSClSizeDe,
                  vtx3body::TPCNClTrackPr, vtx3body::TPCNClTrackPi, vtx3body::TPCNClTrackDe,
                  vtx3body::PIDTrackingDe,

                  // Dynamic columns
                  vtx3body::P<vtx3body::Px, vtx3body::Py, vtx3body::Pz>,
                  vtx3body::Pt<vtx3body::Px, vtx3body::Py>,
                  vtx3body::VtxRadius<vtx3body::X, vtx3body::Y>,
                  vtx3body::DistOverTotMom<vtx3body::X, vtx3body::Y, vtx3body::Z, vtx3body::Px, vtx3body::Py, vtx3body::Pz>,
                  vtx3body::DCAVtxToPV<vtx3body::X, vtx3body::Y, vtx3body::Z, vtx3body::Px, vtx3body::Py, vtx3body::Pz>,

                  // Longitudinal
                  vtx3body::Rap<vtx3body::Px, vtx3body::Py, vtx3body::Pz>,
                  vtx3body::Eta<vtx3body::Px, vtx3body::Py, vtx3body::Pz>,
                  vtx3body::Phi<vtx3body::Px, vtx3body::Py>,
                  vtx3body::TrackPrPt<vtx3body::PxTrackPr, vtx3body::PyTrackPr>,
                  vtx3body::TrackPrEta<vtx3body::PxTrackPr, vtx3body::PyTrackPr, vtx3body::PzTrackPr>,
                  vtx3body::TrackPrPhi<vtx3body::PxTrackPr, vtx3body::PyTrackPr>,
                  vtx3body::TrackPiPt<vtx3body::PxTrackPi, vtx3body::PyTrackPi>,
                  vtx3body::TrackPiEta<vtx3body::PxTrackPi, vtx3body::PyTrackPi, vtx3body::PzTrackPi>,
                  vtx3body::TrackPiPhi<vtx3body::PxTrackPi, vtx3body::PyTrackPi>,
                  vtx3body::TrackDePt<vtx3body::PxTrackDe, vtx3body::PyTrackDe>,
                  vtx3body::TrackDeEta<vtx3body::PxTrackDe, vtx3body::PyTrackDe, vtx3body::PzTrackDe>,
                  vtx3body::TrackDePhi<vtx3body::PxTrackDe, vtx3body::PyTrackDe>);

// covariance matrix table
DECLARE_SOA_TABLE(Vtx3BodyCovs, "AOD", "VTX3BODYCOV", //!
                  vtx3body::CovProton, vtx3body::CovPion, vtx3body::CovDeuteron,
                  vtx3body::VtxCovMat);

// MC candidate table for analysis
DECLARE_SOA_TABLE(McVtx3BodyDatas, "AOD", "MC3BODYDATA", //!
                  o2::soa::Index<>,
                  vtx3body::Sign,
                  vtx3body::Mass, vtx3body::MassV0,
                  vtx3body::X, vtx3body::Y, vtx3body::Z,
                  vtx3body::Px, vtx3body::Py, vtx3body::Pz,
                  vtx3body::Chi2,
                  vtx3body::TrackedClSize,
                  vtx3body::PxTrackPr, vtx3body::PyTrackPr, vtx3body::PzTrackPr,
                  vtx3body::PxTrackPi, vtx3body::PyTrackPi, vtx3body::PzTrackPi,
                  vtx3body::PxTrackDe, vtx3body::PyTrackDe, vtx3body::PzTrackDe,
                  vtx3body::DCAXYTrackPrToPV, vtx3body::DCAXYTrackPiToPV, vtx3body::DCAXYTrackDeToPV,
                  vtx3body::DCATrackPrToPV, vtx3body::DCATrackPiToPV, vtx3body::DCATrackDeToPV,
                  vtx3body::DCAXYTrackPrToPVProp, vtx3body::DCAXYTrackPiToPVProp, vtx3body::DCAXYTrackDeToPVProp,
                  vtx3body::DCATrackPrToPVProp, vtx3body::DCATrackPiToPVProp, vtx3body::DCATrackDeToPVProp,
                  vtx3body::DCATrackPrToSV, vtx3body::DCATrackPiToSV, vtx3body::DCATrackDeToSV,
                  vtx3body::DCAVtxToDaughtersAv,
                  vtx3body::CosPA, vtx3body::Ct,
                  vtx3body::TPCNSigmaPr, vtx3body::TPCNSigmaPi, vtx3body::TPCNSigmaDe, vtx3body::TPCNSigmaPiBach,
                  vtx3body::TOFNSigmaDe,
                  vtx3body::ITSClSizePr, vtx3body::ITSClSizePi, vtx3body::ITSClSizeDe,
                  vtx3body::TPCNClTrackPr, vtx3body::TPCNClTrackPi, vtx3body::TPCNClTrackDe,
                  vtx3body::PIDTrackingDe,

                  // Monte Carlo information
                  vtx3body::GenPx, vtx3body::GenPy, vtx3body::GenPz,
                  vtx3body::GenX, vtx3body::GenY, vtx3body::GenZ,
                  vtx3body::GenCt,
                  vtx3body::GenPhi, vtx3body::GenEta, vtx3body::GenRap,
                  vtx3body::GenPPr, vtx3body::GenPPi, vtx3body::GenPDe,
                  vtx3body::GenPtPr, vtx3body::GenPtPi, vtx3body::GenPtDe,
                  vtx3body::IsTrueH3L, vtx3body::IsTrueAntiH3L,
                  vtx3body::IsReco,
                  vtx3body::MotherPdgCode,
                  vtx3body::PrPdgCode, vtx3body::PiPdgCode, vtx3body::DePdgCode,
                  vtx3body::IsDePrimary,
                  vtx3body::IsSurvEvSel,

                  // Dynamic columns
                  vtx3body::P<vtx3body::Px, vtx3body::Py, vtx3body::Pz>,
                  vtx3body::Pt<vtx3body::Px, vtx3body::Py>,
                  vtx3body::GenP<vtx3body::GenPx, vtx3body::GenPy, vtx3body::GenPz>,
                  vtx3body::GenPt<vtx3body::GenPx, vtx3body::GenPy>,
                  vtx3body::VtxRadius<vtx3body::X, vtx3body::Y>,
                  vtx3body::GenRadius<vtx3body::GenX, vtx3body::GenY>,
                  vtx3body::DistOverTotMom<vtx3body::X, vtx3body::Y, vtx3body::Z, vtx3body::Px, vtx3body::Py, vtx3body::Pz>,
                  vtx3body::DCAVtxToPV<vtx3body::X, vtx3body::Y, vtx3body::Z, vtx3body::Px, vtx3body::Py, vtx3body::Pz>,

                  // Longitudinal
                  vtx3body::Rap<vtx3body::Px, vtx3body::Py, vtx3body::Pz>,
                  vtx3body::Eta<vtx3body::Px, vtx3body::Py, vtx3body::Pz>,
                  vtx3body::Phi<vtx3body::Px, vtx3body::Py>,
                  vtx3body::TrackPrPt<vtx3body::PxTrackPr, vtx3body::PyTrackPr>,
                  vtx3body::TrackPrEta<vtx3body::PxTrackPr, vtx3body::PyTrackPr, vtx3body::PzTrackPr>,
                  vtx3body::TrackPrPhi<vtx3body::PxTrackPr, vtx3body::PyTrackPr>,
                  vtx3body::TrackPiPt<vtx3body::PxTrackPi, vtx3body::PyTrackPi>,
                  vtx3body::TrackPiEta<vtx3body::PxTrackPi, vtx3body::PyTrackPi, vtx3body::PzTrackPi>,
                  vtx3body::TrackPiPhi<vtx3body::PxTrackPi, vtx3body::PyTrackPi>,
                  vtx3body::TrackDePt<vtx3body::PxTrackDe, vtx3body::PyTrackDe>,
                  vtx3body::TrackDeEta<vtx3body::PxTrackDe, vtx3body::PyTrackDe, vtx3body::PzTrackDe>,
                  vtx3body::TrackDePhi<vtx3body::PxTrackDe, vtx3body::PyTrackDe>);

// Define joins
using Vtx3BodyDatasCovs = soa::Join<Vtx3BodyDatas, Vtx3BodyCovs>;
using Vtx3BodyDatasCovsIndexed = soa::Join<Vtx3BodyDatas, Vtx3BodyCovs, Decay3BodyIndices>;

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_VTX3BODYTABLES_H_
