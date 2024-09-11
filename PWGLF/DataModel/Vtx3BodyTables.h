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

#ifndef PWGLF_DATAMODEL_VTX3BODYTABLES_H_
#define PWGLF_DATAMODEL_VTX3BODYTABLES_H_

#include <cmath>
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/RecoDecay.h"
#include "CommonConstants/PhysicsConstants.h"

namespace o2::aod
{
namespace vtx3body
{
DECLARE_SOA_INDEX_COLUMN_FULL(Track0, track0, int, Tracks, "_0"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(Track1, track1, int, Tracks, "_1"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(Track2, track2, int, Tracks, "_2"); //!
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                   //!
DECLARE_SOA_INDEX_COLUMN(Decay3Body, decay3body);                 //!

// General 3 body Vtx properties: position, momentum
DECLARE_SOA_COLUMN(PxTrack0, pxtrack0, float); //! track0 px at min
DECLARE_SOA_COLUMN(PyTrack0, pytrack0, float); //! track0 py at min
DECLARE_SOA_COLUMN(PzTrack0, pztrack0, float); //! track0 pz at min
DECLARE_SOA_COLUMN(PxTrack1, pxtrack1, float); //! track1 px at min
DECLARE_SOA_COLUMN(PyTrack1, pytrack1, float); //! track1 py at min
DECLARE_SOA_COLUMN(PzTrack1, pztrack1, float); //! track1 pz at min
DECLARE_SOA_COLUMN(PxTrack2, pxtrack2, float); //! track2 px at min
DECLARE_SOA_COLUMN(PyTrack2, pytrack2, float); //! track2 py at min
DECLARE_SOA_COLUMN(PzTrack2, pztrack2, float); //! track2 pz at min
DECLARE_SOA_COLUMN(X, x, float);               //! decay position X
DECLARE_SOA_COLUMN(Y, y, float);               //! decay position Y
DECLARE_SOA_COLUMN(Z, z, float);               //! decay position Z

// Saved from finding: DCAs
DECLARE_SOA_COLUMN(DCAVtxDaughters, dcaVtxdaughters, float); //! DCA among daughters
DECLARE_SOA_COLUMN(DCAXYTrack0ToPV, dcaXYtrack0topv, float); //! DCAXY of prong0 to PV
DECLARE_SOA_COLUMN(DCAXYTrack1ToPV, dcaXYtrack1topv, float); //! DCAXY of prong1 to PV
DECLARE_SOA_COLUMN(DCAXYTrack2ToPV, dcaXYtrack2topv, float); //! DCAXY of prong2 to PV
DECLARE_SOA_COLUMN(DCATrack0ToPV, dcatrack0topv, float);     //! DCA of prong0 to PV
DECLARE_SOA_COLUMN(DCATrack1ToPV, dcatrack1topv, float);     //! DCA of prong1 to PV
DECLARE_SOA_COLUMN(DCATrack2ToPV, dcatrack2topv, float);     //! DCA of prong2 to PV

// Recalculated TOF PID information of bachelor
DECLARE_SOA_COLUMN(TOFNSigmaBachDe, tofNSigmaBachDe, float); //! Recalculated Nsigma seperation with TOF for deuteron

// Derived expressions
// Momenta
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //! 3 body p
                           [](float pxtrack0, float pytrack0, float pztrack0, float pxtrack1, float pytrack1, float pztrack1, float pxtrack2, float pytrack2, float pztrack2) -> float { return RecoDecay::sqrtSumOfSquares(pxtrack0 + pxtrack1 + pxtrack2, pytrack0 + pytrack1 + pytrack2, pztrack0 + pztrack1 + pztrack2); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //! 3 body pT
                           [](float pxtrack0, float pytrack0, float pxtrack1, float pytrack1, float pxtrack2, float pytrack2) -> float { return RecoDecay::sqrtSumOfSquares(pxtrack0 + pxtrack1 + pxtrack2, pytrack0 + pytrack1 + pytrack2); });

// Length quantities
DECLARE_SOA_DYNAMIC_COLUMN(VtxRadius, vtxradius, //! 3 body decay radius (2D, centered at zero)
                           [](float x, float y) -> float { return RecoDecay::sqrtSumOfSquares(x, y); });

// Distance Over To Mom
DECLARE_SOA_DYNAMIC_COLUMN(DistOverTotMom, distovertotmom, //! PV to 3 body decay distance over total momentum
                           [](float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ) {
                             float P = RecoDecay::sqrtSumOfSquares(Px, Py, Pz);
                             return std::sqrt(std::pow(X - pvX, 2) + std::pow(Y - pvY, 2) + std::pow(Z - pvZ, 2)) / (P + 1E-10);
                           });

// CosPA
DECLARE_SOA_DYNAMIC_COLUMN(VtxCosPA, vtxcosPA, //! 3 body vtx CosPA
                           [](float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ) -> float { return RecoDecay::cpa(std::array{pvX, pvY, pvZ}, std::array{X, Y, Z}, std::array{Px, Py, Pz}); });
// Dca to PV
DECLARE_SOA_DYNAMIC_COLUMN(DCAVtxToPV, dcavtxtopv, //! DCA of 3 body vtx to PV
                           [](float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ) -> float { return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz)); });

// Eta and Phi
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, //! 3 body vtx eta
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::eta(std::array{Px, Py, Pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, //! 3 body vtx phi
                           [](float Px, float Py) -> float { return RecoDecay::phi(Px, Py); });

// Calculated on the fly with mother particle hypothesis
DECLARE_SOA_DYNAMIC_COLUMN(MHypertriton, mHypertriton, //! mass under Hypertriton hypothesis
                           [](float pxtrack0, float pytrack0, float pztrack0, float pxtrack1, float pytrack1, float pztrack1, float pxtrack2, float pytrack2, float pztrack2) -> float { return RecoDecay::m(std::array{std::array{pxtrack0, pytrack0, pztrack0}, std::array{pxtrack1, pytrack1, pztrack1}, std::array{pxtrack2, pytrack2, pztrack2}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged, o2::constants::physics::MassDeuteron}); });
DECLARE_SOA_DYNAMIC_COLUMN(MAntiHypertriton, mAntiHypertriton, //! mass under antiHypertriton hypothesis
                           [](float pxtrack0, float pytrack0, float pztrack0, float pxtrack1, float pytrack1, float pztrack1, float pxtrack2, float pytrack2, float pztrack2) -> float { return RecoDecay::m(std::array{std::array{pxtrack0, pytrack0, pztrack0}, std::array{pxtrack1, pytrack1, pztrack1}, std::array{pxtrack2, pytrack2, pztrack2}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron}); });

DECLARE_SOA_DYNAMIC_COLUMN(MHyperHelium4, mHyperHelium4, //! mass under HyperHelium4 hypothesis
                           [](float pxtrack0, float pytrack0, float pztrack0, float pxtrack1, float pytrack1, float pztrack1, float pxtrack2, float pytrack2, float pztrack2) -> float { return RecoDecay::m(std::array{std::array{pxtrack0, pytrack0, pztrack0}, std::array{pxtrack1, pytrack1, pztrack1}, std::array{pxtrack2, pytrack2, pztrack2}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged, o2::constants::physics::MassHelium3}); });
DECLARE_SOA_DYNAMIC_COLUMN(MAntiHyperHelium4, mAntiHyperHelium4, //! mass under antiHyperHelium4 hypothesis
                           [](float pxtrack0, float pytrack0, float pztrack0, float pxtrack1, float pytrack1, float pztrack1, float pxtrack2, float pytrack2, float pztrack2) -> float { return RecoDecay::m(std::array{std::array{pxtrack0, pytrack0, pztrack0}, std::array{pxtrack1, pytrack1, pztrack1}, std::array{pxtrack2, pytrack2, pztrack2}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton, o2::constants::physics::MassHelium3}); });

DECLARE_SOA_DYNAMIC_COLUMN(YHypertriton, yHypertriton, //! 3 body vtx y with hypertriton or antihypertriton hypothesis
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::y(std::array{Px, Py, Pz}, o2::constants::physics::MassHyperTriton); });
DECLARE_SOA_DYNAMIC_COLUMN(YHyperHelium4, yHyperHelium4, //! 3 body vtx y with hyperhelium4 or antihyperhelium4 hypothesis
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::y(std::array{Px, Py, Pz}, o2::constants::physics::MassHyperHelium4); });

// kinematic information of daughter tracks
DECLARE_SOA_DYNAMIC_COLUMN(Track0Pt, track0pt, //! daughter0 pT
                           [](float pxtrack0, float pytrack0) -> float { return RecoDecay::sqrtSumOfSquares(pxtrack0, pytrack0); });
DECLARE_SOA_DYNAMIC_COLUMN(Track1Pt, track1pt, //! daughter1 pT
                           [](float pxtrack1, float pytrack1) -> float { return RecoDecay::sqrtSumOfSquares(pxtrack1, pytrack1); });
DECLARE_SOA_DYNAMIC_COLUMN(Track2Pt, track2pt, //! daughter2 pT
                           [](float pxtrack2, float pytrack2) -> float { return RecoDecay::sqrtSumOfSquares(pxtrack2, pytrack2); });
DECLARE_SOA_DYNAMIC_COLUMN(Track0Eta, track0eta, //! daughter0 eta
                           [](float pxtrack0, float pytrack0, float pztrack0) -> float { return RecoDecay::eta(std::array{pxtrack0, pytrack0, pztrack0}); });
DECLARE_SOA_DYNAMIC_COLUMN(Track0Phi, track0phi, //! daughter0 phi
                           [](float pxtrack0, float pytrack0) -> float { return RecoDecay::phi(pxtrack0, pytrack0); });
DECLARE_SOA_DYNAMIC_COLUMN(Track1Eta, track1eta, //! daughter1 eta
                           [](float pxtrack1, float pytrack1, float pztrack1) -> float { return RecoDecay::eta(std::array{pxtrack1, pytrack1, pztrack1}); });
DECLARE_SOA_DYNAMIC_COLUMN(Track1Phi, track1phi, //! daughter1 phi
                           [](float pxtrack1, float pytrack1) -> float { return RecoDecay::phi(pxtrack1, pytrack1); });
DECLARE_SOA_DYNAMIC_COLUMN(Track2Eta, track2eta, //! daughter2 eta
                           [](float pxtrack2, float pytrack2, float pztrack2) -> float { return RecoDecay::eta(std::array{pxtrack2, pytrack2, pztrack2}); });
DECLARE_SOA_DYNAMIC_COLUMN(Track2Phi, track2phi, //! daughter2 phi
                           [](float pxtrack2, float pytrack2) -> float { return RecoDecay::phi(pxtrack2, pytrack2); });

DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //! 3 body vtx px
                              float, 1.f * aod::vtx3body::pxtrack0 + 1.f * aod::vtx3body::pxtrack1 + 1.f * aod::vtx3body::pxtrack2);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //! 3 body vtx py
                              float, 1.f * aod::vtx3body::pytrack0 + 1.f * aod::vtx3body::pytrack1 + 1.f * aod::vtx3body::pytrack2);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //! 3 body vtx pz
                              float, 1.f * aod::vtx3body::pztrack0 + 1.f * aod::vtx3body::pztrack1 + 1.f * aod::vtx3body::pztrack2);
} // namespace vtx3body

DECLARE_SOA_TABLE_FULL(StoredVtx3BodyDatas, "Vtx3BodyDatas", "AOD", "Vtx3BodyDATA", //!
                       o2::soa::Index<>, vtx3body::Track0Id, vtx3body::Track1Id, vtx3body::Track2Id, vtx3body::CollisionId, vtx3body::Decay3BodyId,
                       vtx3body::X, vtx3body::Y, vtx3body::Z,
                       vtx3body::PxTrack0, vtx3body::PyTrack0, vtx3body::PzTrack0,
                       vtx3body::PxTrack1, vtx3body::PyTrack1, vtx3body::PzTrack1,
                       vtx3body::PxTrack2, vtx3body::PyTrack2, vtx3body::PzTrack2,
                       vtx3body::DCAVtxDaughters,
                       vtx3body::DCAXYTrack0ToPV, vtx3body::DCAXYTrack1ToPV, vtx3body::DCAXYTrack2ToPV,
                       vtx3body::DCATrack0ToPV, vtx3body::DCATrack1ToPV, vtx3body::DCATrack2ToPV,
                       vtx3body::TOFNSigmaBachDe,

                       // Dynamic columns
                       vtx3body::P<vtx3body::PxTrack0, vtx3body::PyTrack0, vtx3body::PzTrack0, vtx3body::PxTrack1, vtx3body::PyTrack1, vtx3body::PzTrack1, vtx3body::PxTrack2, vtx3body::PyTrack2, vtx3body::PzTrack2>,
                       vtx3body::Pt<vtx3body::PxTrack0, vtx3body::PyTrack0, vtx3body::PxTrack1, vtx3body::PyTrack1, vtx3body::PxTrack2, vtx3body::PyTrack2>,
                       vtx3body::VtxRadius<vtx3body::X, vtx3body::Y>,
                       vtx3body::DistOverTotMom<vtx3body::X, vtx3body::Y, vtx3body::Z, vtx3body::Px, vtx3body::Py, vtx3body::Pz>,
                       vtx3body::VtxCosPA<vtx3body::X, vtx3body::Y, vtx3body::Z, vtx3body::Px, vtx3body::Py, vtx3body::Pz>,
                       vtx3body::DCAVtxToPV<vtx3body::X, vtx3body::Y, vtx3body::Z, vtx3body::Px, vtx3body::Py, vtx3body::Pz>,

                       // Invariant masses
                       vtx3body::MHypertriton<vtx3body::PxTrack0, vtx3body::PyTrack0, vtx3body::PzTrack0, vtx3body::PxTrack1, vtx3body::PyTrack1, vtx3body::PzTrack1, vtx3body::PxTrack2, vtx3body::PyTrack2, vtx3body::PzTrack2>,
                       vtx3body::MAntiHypertriton<vtx3body::PxTrack0, vtx3body::PyTrack0, vtx3body::PzTrack0, vtx3body::PxTrack1, vtx3body::PyTrack1, vtx3body::PzTrack1, vtx3body::PxTrack2, vtx3body::PyTrack2, vtx3body::PzTrack2>,
                       vtx3body::MHyperHelium4<vtx3body::PxTrack0, vtx3body::PyTrack0, vtx3body::PzTrack0, vtx3body::PxTrack1, vtx3body::PyTrack1, vtx3body::PzTrack1, vtx3body::PxTrack2, vtx3body::PyTrack2, vtx3body::PzTrack2>,
                       vtx3body::MAntiHyperHelium4<vtx3body::PxTrack0, vtx3body::PyTrack0, vtx3body::PzTrack0, vtx3body::PxTrack1, vtx3body::PyTrack1, vtx3body::PzTrack1, vtx3body::PxTrack2, vtx3body::PyTrack2, vtx3body::PzTrack2>,

                       // Longitudinal
                       vtx3body::YHypertriton<vtx3body::Px, vtx3body::Py, vtx3body::Pz>,
                       vtx3body::YHyperHelium4<vtx3body::Px, vtx3body::Py, vtx3body::Pz>,
                       vtx3body::Eta<vtx3body::Px, vtx3body::Py, vtx3body::Pz>,
                       vtx3body::Phi<vtx3body::Px, vtx3body::Py>,
                       vtx3body::Track0Pt<vtx3body::PxTrack0, vtx3body::PyTrack0>,
                       vtx3body::Track0Eta<vtx3body::PxTrack0, vtx3body::PyTrack0, vtx3body::PzTrack0>,
                       vtx3body::Track0Phi<vtx3body::PxTrack0, vtx3body::PyTrack0>,
                       vtx3body::Track1Pt<vtx3body::PxTrack1, vtx3body::PyTrack1>,
                       vtx3body::Track1Eta<vtx3body::PxTrack1, vtx3body::PyTrack1, vtx3body::PzTrack1>,
                       vtx3body::Track1Phi<vtx3body::PxTrack1, vtx3body::PyTrack1>,
                       vtx3body::Track2Pt<vtx3body::PxTrack2, vtx3body::PyTrack2>,
                       vtx3body::Track2Eta<vtx3body::PxTrack2, vtx3body::PyTrack2, vtx3body::PzTrack2>,
                       vtx3body::Track2Phi<vtx3body::PxTrack2, vtx3body::PyTrack2>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(Vtx3BodyDatas, StoredVtx3BodyDatas, "Vtx3BodyDATAEXT", //!
                                vtx3body::Px, vtx3body::Py, vtx3body::Pz);

using Vtx3BodyData = Vtx3BodyDatas::iterator;
namespace vtx3body
{
DECLARE_SOA_INDEX_COLUMN(Vtx3BodyData, vtx3BodyData); //! Index to Vtx3BodyData entry
}

DECLARE_SOA_TABLE(Decay3BodyDataLink, "AOD", "DECAY3BODYLINK", //! Joinable table with Decay3bodys which links to Vtx3BodyData which is not produced for all entries
                  vtx3body::Vtx3BodyDataId);

using Decay3BodysLinked = soa::Join<Decay3Bodys, Decay3BodyDataLink>;
using Decay3BodyLinked = Decay3BodysLinked::iterator;

// Definition of labels for Vtx3BodyDatas
namespace mcvtx3bodylabel
{
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle); //! MC particle for Vtx3BodyDatas
} // namespace mcvtx3bodylabel

DECLARE_SOA_TABLE(McVtx3BodyLabels, "AOD", "MCVTXLABEL", //! Table joinable with Vtx3BodyData containing the MC labels
                  mcvtx3bodylabel::McParticleId);
using McVtx3BodyLabel = McVtx3BodyLabels::iterator;

// Definition of labels for Decay3Bodys // Full table, joinable with Decay3Bodys (CAUTION: NOT WITH Vtx3BodyDATA)
namespace mcfullvtx3bodylabel
{
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle); //! MC particle for Decay3Bodys
} // namespace mcfullvtx3bodylabel

DECLARE_SOA_TABLE(McFullVtx3BodyLabels, "AOD", "MCFULLVTXLABEL", //! Table joinable with Decay3Bodys
                  mcfullvtx3bodylabel::McParticleId);
using McFullVtx3BodyLabel = McFullVtx3BodyLabels::iterator;

// output table for ML studies
namespace hyp3body
{
// collision
DECLARE_SOA_COLUMN(Centrality, centrality, float); //! centrality
DECLARE_SOA_COLUMN(XPV, xpv, float);               //! primary vertex X
DECLARE_SOA_COLUMN(YPV, ypv, float);               //! primary vertex Y
DECLARE_SOA_COLUMN(ZPV, zpv, float);               //! primary vertex Z
// reconstruced candidate
DECLARE_SOA_COLUMN(IsMatter, isMatter, bool); //! bool: true for matter
DECLARE_SOA_COLUMN(M, m, float);              //! invariant mass
DECLARE_SOA_COLUMN(P, p, float);              //! p
DECLARE_SOA_COLUMN(Pt, pt, float);            //! pT
DECLARE_SOA_COLUMN(Ct, ct, float);            //! ct
DECLARE_SOA_COLUMN(X, x, float);              //! decay position X
DECLARE_SOA_COLUMN(Y, y, float);              //! decay position Y
DECLARE_SOA_COLUMN(Z, z, float);              //! decay position Z
DECLARE_SOA_COLUMN(CosPA, cospa, float);
DECLARE_SOA_COLUMN(DCADaughters, dcaDaughters, float); //! DCA among daughters
DECLARE_SOA_COLUMN(DCACandToPV, dcaCandtopv, float);   //! DCA of the reconstructed track to pv
// kinematic infomation of daughter tracks
DECLARE_SOA_COLUMN(PProton, pProton, float);         //! p of the proton daughter
DECLARE_SOA_COLUMN(PtProton, ptProton, float);       //! pT of the proton daughter
DECLARE_SOA_COLUMN(EtaProton, etaProton, float);     //! eta of the proton daughter
DECLARE_SOA_COLUMN(PhiProton, phiProton, float);     //! phi of the proton daughter
DECLARE_SOA_COLUMN(PPion, pPion, float);             //! p of the pion daughter
DECLARE_SOA_COLUMN(PtPion, ptPion, float);           //! pT of the pion daughter
DECLARE_SOA_COLUMN(EtaPion, etaPion, float);         //! eta of the pion daughter
DECLARE_SOA_COLUMN(PhiPion, phiPion, float);         //! phi of the pion daughter
DECLARE_SOA_COLUMN(PBachelor, pBachelor, float);     //! p of the bachelor daughter
DECLARE_SOA_COLUMN(PtBachelor, ptBachelor, float);   //! pT of the bachelor daughter
DECLARE_SOA_COLUMN(EtaBachelor, etaBachelor, float); //! eta of the bachelor daughter
DECLARE_SOA_COLUMN(PhiBachelor, phiBachelor, float); //! phi of the bachelor daughter
// track quality
DECLARE_SOA_COLUMN(TPCNclusProton, tpcNclusProton, uint8_t);             //! number of TPC clusters of the proton daughter
DECLARE_SOA_COLUMN(TPCNclusPion, tpcNclusPion, uint8_t);                 //! number of TPC clusters of the pion daughter
DECLARE_SOA_COLUMN(TPCNclusBachelor, tpcNclusBachelor, uint8_t);         //! number of TPC clusters of the bachelor daughter
DECLARE_SOA_COLUMN(ITSNclusSizeProton, itsNclusSizeProton, uint8_t);     //! average ITS cluster size of the proton daughter
DECLARE_SOA_COLUMN(ITSNclusSizePion, itsNclusSizePion, uint8_t);         //! average ITS cluster size of the pion daughter
DECLARE_SOA_COLUMN(ITSNclusSizeBachelor, itsNclusSizeBachelor, uint8_t); //! average ITS cluster size of the bachelor daughter
// PID
DECLARE_SOA_COLUMN(TPCNSigmaProton, tpcNSigmaProton, float);     //! nsigma of TPC PID of the proton daughter
DECLARE_SOA_COLUMN(TPCNSigmaPion, tpcNSigmaPion, float);         //! nsigma of TPC PID of the pion daughter
DECLARE_SOA_COLUMN(TPCNSigmaBachelor, tpcNSigmaBachelor, float); //! nsigma of TPC PID of the bachelor daughter
DECLARE_SOA_COLUMN(TOFNSigmaBachelor, tofNSigmaBachelor, float); //! nsigma of TOF PID of the bachelor daughter
// DCA to PV
DECLARE_SOA_COLUMN(DCAXYProtonToPV, dcaxyProtontoPV, float);     //! DCAXY of the proton daughter to pv
DECLARE_SOA_COLUMN(DCAXYPionToPV, dcaxyPiontoPV, float);         //! DCAXY of the pion daughter to pv
DECLARE_SOA_COLUMN(DCAXYBachelorToPV, dcaxyBachelortoPV, float); //! DCAXY of the bachelor daughter to pv
DECLARE_SOA_COLUMN(DCAProtonToPV, dcaProtontoPV, float);         //! DCA of the proton daughter to pv
DECLARE_SOA_COLUMN(DCAPionToPV, dcaPiontoPV, float);             //! DCA of the pion daughter to pv
DECLARE_SOA_COLUMN(DCABachelorToPV, dcaBachelortoPV, float);     //! DCA of the bachelor daughter to pv
// for MC
DECLARE_SOA_COLUMN(GenP, genP, float);                                    // P of the hypertriton
DECLARE_SOA_COLUMN(GenPt, genPt, float);                                  // pT of the hypertriton
DECLARE_SOA_COLUMN(GenCt, genCt, float);                                  // ct of the hypertriton
DECLARE_SOA_COLUMN(GenPhi, genPhi, float);                                // Phi of the hypertriton
DECLARE_SOA_COLUMN(GenEta, genEta, float);                                // Eta of the hypertriton
DECLARE_SOA_COLUMN(GenRapidity, genRapidity, float);                      // Rapidity of the hypertriton
DECLARE_SOA_COLUMN(IsReco, isReco, bool);                                 // bool: true for reco
DECLARE_SOA_COLUMN(IsSignal, isSignal, bool);                             // bool: true for signal
DECLARE_SOA_COLUMN(PdgCode, pdgCode, int);                                // pdgCode of the mcparticle, -1 for fake pair
DECLARE_SOA_COLUMN(SurvivedEventSelection, survivedEventSelection, bool); // bool: true for survived event selection
} // namespace hyp3body

// output table for data
DECLARE_SOA_TABLE(Hyp3BodyCands, "AOD", "HYP3BODYCANDS",
                  o2::soa::Index<>,
                  hyp3body::Centrality,
                  hyp3body::XPV, hyp3body::YPV, hyp3body::ZPV,
                  // secondary vertex and reconstruced candidate
                  hyp3body::IsMatter,
                  hyp3body::M,
                  hyp3body::P,
                  hyp3body::Pt,
                  hyp3body::Ct,
                  hyp3body::X, hyp3body::Y, hyp3body::Z,
                  hyp3body::CosPA,
                  hyp3body::DCADaughters,
                  hyp3body::DCACandToPV,
                  // daughter tracks
                  hyp3body::PProton, hyp3body::PtProton, hyp3body::EtaProton, hyp3body::PhiProton,
                  hyp3body::PPion, hyp3body::PtPion, hyp3body::EtaPion, hyp3body::PhiPion,
                  hyp3body::PBachelor, hyp3body::PtBachelor, hyp3body::EtaBachelor, hyp3body::PhiBachelor,
                  hyp3body::TPCNclusProton, hyp3body::TPCNclusPion, hyp3body::TPCNclusBachelor,
                  hyp3body::ITSNclusSizeProton, hyp3body::ITSNclusSizePion, hyp3body::ITSNclusSizeBachelor,
                  hyp3body::TPCNSigmaProton, hyp3body::TPCNSigmaPion, hyp3body::TPCNSigmaBachelor,
                  hyp3body::TOFNSigmaBachelor,
                  hyp3body::DCAXYProtonToPV, hyp3body::DCAXYPionToPV, hyp3body::DCAXYBachelorToPV,
                  hyp3body::DCAProtonToPV, hyp3body::DCAPionToPV, hyp3body::DCABachelorToPV);

// output table for MC
DECLARE_SOA_TABLE(MCHyp3BodyCands, "AOD", "MCHYP3BODYCANDS",
                  o2::soa::Index<>,
                  hyp3body::Centrality,
                  hyp3body::XPV, hyp3body::YPV, hyp3body::ZPV,
                  // secondary vertex and reconstruced candidate
                  hyp3body::IsMatter,
                  hyp3body::M,
                  hyp3body::P,
                  hyp3body::Pt,
                  hyp3body::Ct,
                  hyp3body::X, hyp3body::Y, hyp3body::Z,
                  hyp3body::CosPA,
                  hyp3body::DCADaughters,
                  hyp3body::DCACandToPV,
                  // daughter tracks
                  hyp3body::PProton, hyp3body::PtProton, hyp3body::EtaProton, hyp3body::PhiProton,
                  hyp3body::PPion, hyp3body::PtPion, hyp3body::EtaPion, hyp3body::PhiPion,
                  hyp3body::PBachelor, hyp3body::PtBachelor, hyp3body::EtaBachelor, hyp3body::PhiBachelor,
                  hyp3body::TPCNclusProton, hyp3body::TPCNclusPion, hyp3body::TPCNclusBachelor,
                  hyp3body::ITSNclusSizeProton, hyp3body::ITSNclusSizePion, hyp3body::ITSNclusSizeBachelor,
                  hyp3body::TPCNSigmaProton, hyp3body::TPCNSigmaPion, hyp3body::TPCNSigmaBachelor,
                  hyp3body::TOFNSigmaBachelor,
                  hyp3body::DCAXYProtonToPV, hyp3body::DCAXYPionToPV, hyp3body::DCAXYBachelorToPV,
                  hyp3body::DCAProtonToPV, hyp3body::DCAPionToPV, hyp3body::DCABachelorToPV,
                  // MC information
                  hyp3body::GenP,
                  hyp3body::GenPt,
                  hyp3body::GenCt,
                  hyp3body::GenPhi,
                  hyp3body::GenEta,
                  hyp3body::GenRapidity,
                  hyp3body::IsSignal,
                  hyp3body::IsReco,
                  hyp3body::PdgCode,
                  hyp3body::SurvivedEventSelection);

//______________________________________________________
// DATAMODEL FOR KFPARTICLE DECAY3BODYS

namespace kfvtx3body
{
// General 3 body Vtx properties: mass, momentum, charge
DECLARE_SOA_COLUMN(Mass, mass, float); //! candidate mass (PID hypothesis depending on bachelor charge)
DECLARE_SOA_COLUMN(Px, px, float);     //! candidate px at decay position
DECLARE_SOA_COLUMN(Py, py, float);     //! candidate py at decay position
DECLARE_SOA_COLUMN(Pz, pz, float);     //! candidate pz at decay position
DECLARE_SOA_COLUMN(Pt, pt, float);     //! candidate pt at decay position
DECLARE_SOA_COLUMN(Sign, sign, float); //! candidate sign

// topological properties
DECLARE_SOA_COLUMN(VtxCosPAKF, vtxcosPAkf, float);             //! 3 body vtx CosPA from KFParticle (using kfpPV)
DECLARE_SOA_COLUMN(VtxCosXYPAKF, vtxcosxyPAkf, float);         //! 3 body vtx CosPA from KFParticle (using kfpPV)
DECLARE_SOA_COLUMN(VtxCosPAKFtopo, vtxcosPAkftopo, float);     //! 3 body vtx CosPA from KFParticle after topological constraint (using kfpPV)
DECLARE_SOA_COLUMN(VtxCosXYPAKFtopo, vtxcosxyPAkftopo, float); //! 3 body vtx CosPA from KFParticle after topological constraint (using kfpPV)
DECLARE_SOA_COLUMN(DCAVtxToPVKF, dcavtxtopvkf, float);         //! 3 body vtx DCA to PV from KFParticle (using kfpPV)
DECLARE_SOA_COLUMN(DCAXYVtxToPVKF, dcaxyvtxtopvkf, float);     //! 3 body vtx DCAxy to PV from KFParticle (using kfpPV)
DECLARE_SOA_COLUMN(DecayLKF, decaylkf, float);                 //! 3 body vtx decay length from KFParticle (using kfpPV after topological constraint)
DECLARE_SOA_COLUMN(DecayLXYKF, decaylxykf, float);             //! 3 body vtx decay length XY from KFParticle (using kfpPV after topological constraint)
DECLARE_SOA_COLUMN(DecayLDeltaL, decayldeltal, float);         //! 3 body vtx l/dl from KFParticle (using kfpPV after topological constraint)
DECLARE_SOA_COLUMN(Chi2geoNDF, chi2geondf, float);             //! 3 body vtx chi2geo from geometrical KFParticle fit
DECLARE_SOA_COLUMN(Chi2topoNDF, chi2topondf, float);           //! 3 body vtx chi2topo from KFParticle topological constraint to the PV (using kfpPV)

// daughters
DECLARE_SOA_COLUMN(DCATrack0ToPVKF, dcatrack0topvkf, float);           //! DCA of proton prong to PV from KFParticle
DECLARE_SOA_COLUMN(DCATrack1ToPVKF, dcatrack1topvkf, float);           //! DCA of pion prong to PV from KFParticle
DECLARE_SOA_COLUMN(DCATrack2ToPVKF, dcatrack2topvkf, float);           //! DCA of deuteron prong to PV from KFParticle
DECLARE_SOA_COLUMN(DCAxyTrack0ToPVKF, dcaxytrack0topvkf, float);       //! DCAxy of proton prong to PV from KFParticle
DECLARE_SOA_COLUMN(DCAxyTrack1ToPVKF, dcaxytrack1topvkf, float);       //! DCAxy of pion prong to PV from KFParticle
DECLARE_SOA_COLUMN(DCAxyTrack2ToPVKF, dcaxytrack2topvkf, float);       //! DCAxy of deuteron prong to PV from KFParticle
DECLARE_SOA_COLUMN(DCAxyTrack0ToSVKF, dcaxytrack0tosvkf, float);       //! DCAxy of proton prong to SV from KFParticle
DECLARE_SOA_COLUMN(DCAxyTrack1ToSVKF, dcaxytrack1tosvkf, float);       //! DCAxy of pion prong to SV from KFParticle
DECLARE_SOA_COLUMN(DCAxyTrack2ToSVKF, dcaxytrack2tosvkf, float);       //! DCAxy of deuteron prong to SV from KFParticle
DECLARE_SOA_COLUMN(DCAxyTrack0ToTrack1, dcaxytrack0totrack1kf, float); //! DCAxy of proton prong to pion from KFParticle
DECLARE_SOA_COLUMN(DCAxyTrack0ToTrack2, dcaxytrack0totrack2kf, float); //! DCAxy of proton prong to deuteron from KFParticle
DECLARE_SOA_COLUMN(DCAxyTrack1ToTrack2, dcaxytrack1totrack2kf, float); //! DCAxy of pion prong to deuteron from KFParticle
DECLARE_SOA_COLUMN(DCAVtxDaughtersKF, dcavtxdaughterskf, float);       //! sum of DCAs between daughters in 3D from KFParticle
DECLARE_SOA_COLUMN(Track0Sign, track0sign, float);                     //! sign of proton daughter track
DECLARE_SOA_COLUMN(Track1Sign, track1sign, float);                     //! sign of pion daughter track
DECLARE_SOA_COLUMN(Track2Sign, track2sign, float);                     //! sign of deuteron daughter track

// V0
DECLARE_SOA_COLUMN(MassV0, massv0, float);         //! proton, pion vertex mass
DECLARE_SOA_COLUMN(Chi2MassV0, chi2massv0, float); //! chi2 of proton, pion mass constraint to Lambda mass

} // namespace kfvtx3body

DECLARE_SOA_TABLE(StoredKFVtx3BodyDatas, "AOD", "KFVTX3BODYDATA", //!
                                                                  // indices
                  o2::soa::Index<>, vtx3body::CollisionId, vtx3body::Track0Id, vtx3body::Track1Id, vtx3body::Track2Id, vtx3body::Decay3BodyId,

                  // hypertriton candidate
                  kfvtx3body::Mass,
                  vtx3body::X, vtx3body::Y, vtx3body::Z,
                  kfvtx3body::Px, kfvtx3body::Py, kfvtx3body::Pz, kfvtx3body::Pt,
                  kfvtx3body::Sign,
                  kfvtx3body::DCAVtxToPVKF, kfvtx3body::DCAXYVtxToPVKF,
                  kfvtx3body::VtxCosPAKF, kfvtx3body::VtxCosXYPAKF,
                  kfvtx3body::VtxCosPAKFtopo, kfvtx3body::VtxCosXYPAKFtopo,
                  kfvtx3body::DecayLKF, kfvtx3body::DecayLXYKF, kfvtx3body::DecayLDeltaL,
                  kfvtx3body::Chi2geoNDF, kfvtx3body::Chi2topoNDF,

                  // V0
                  kfvtx3body::MassV0, kfvtx3body::Chi2MassV0,

                  // daughters
                  vtx3body::PxTrack0, vtx3body::PyTrack0, vtx3body::PzTrack0, // proton
                  vtx3body::PxTrack1, vtx3body::PyTrack1, vtx3body::PzTrack1, // pion
                  vtx3body::PxTrack2, vtx3body::PyTrack2, vtx3body::PzTrack2, // deuteron
                  kfvtx3body::DCATrack0ToPVKF, kfvtx3body::DCATrack1ToPVKF, kfvtx3body::DCATrack2ToPVKF, kfvtx3body::DCAxyTrack0ToPVKF, kfvtx3body::DCAxyTrack1ToPVKF, kfvtx3body::DCAxyTrack2ToPVKF,
                  kfvtx3body::DCAxyTrack0ToSVKF, kfvtx3body::DCAxyTrack1ToSVKF, kfvtx3body::DCAxyTrack2ToSVKF,
                  kfvtx3body::DCAxyTrack0ToTrack1, kfvtx3body::DCAxyTrack0ToTrack2, kfvtx3body::DCAxyTrack1ToTrack2,
                  kfvtx3body::DCAVtxDaughtersKF,
                  vtx3body::DCATrack0ToPV, vtx3body::DCATrack1ToPV, vtx3body::DCATrack2ToPV,
                  kfvtx3body::Track0Sign, kfvtx3body::Track1Sign, kfvtx3body::Track2Sign, // track sing: proton, pion, deuteron
                  vtx3body::TOFNSigmaBachDe,

                  // dynamic columns
                  vtx3body::VtxRadius<vtx3body::X, vtx3body::Y>,
                  vtx3body::DistOverTotMom<vtx3body::X, vtx3body::Y, vtx3body::Z, kfvtx3body::Px, kfvtx3body::Py, kfvtx3body::Pz>,
                  vtx3body::VtxCosPA<vtx3body::X, vtx3body::Y, vtx3body::Z, kfvtx3body::Px, kfvtx3body::Py, kfvtx3body::Pz>,
                  vtx3body::DCAVtxToPV<vtx3body::X, vtx3body::Y, vtx3body::Z, kfvtx3body::Px, kfvtx3body::Py, kfvtx3body::Pz>,

                  // Longitudinal
                  vtx3body::YHypertriton<kfvtx3body::Px, kfvtx3body::Py, kfvtx3body::Pz>,
                  vtx3body::Eta<kfvtx3body::Px, kfvtx3body::Py, kfvtx3body::Pz>,
                  vtx3body::Phi<kfvtx3body::Px, kfvtx3body::Py>,
                  vtx3body::Track0Pt<vtx3body::PxTrack0, vtx3body::PyTrack0>,                      // proton pT
                  vtx3body::Track0Eta<vtx3body::PxTrack0, vtx3body::PyTrack0, vtx3body::PzTrack0>, // proton eta
                  vtx3body::Track0Phi<vtx3body::PxTrack0, vtx3body::PyTrack0>,                     // proton phi
                  vtx3body::Track1Pt<vtx3body::PxTrack1, vtx3body::PyTrack1>,                      // pion pT
                  vtx3body::Track1Eta<vtx3body::PxTrack1, vtx3body::PyTrack1, vtx3body::PzTrack1>, // pion eta
                  vtx3body::Track1Phi<vtx3body::PxTrack1, vtx3body::PyTrack1>,                     // pion phi
                  vtx3body::Track2Pt<vtx3body::PxTrack2, vtx3body::PyTrack2>,                      // deuteron pT
                  vtx3body::Track2Eta<vtx3body::PxTrack2, vtx3body::PyTrack2, vtx3body::PzTrack2>, // deuteron eta
                  vtx3body::Track2Phi<vtx3body::PxTrack2, vtx3body::PyTrack2>);                    // deuteron phi

using StoredKFVtx3BodyData = StoredKFVtx3BodyDatas::iterator;
namespace kfvtx3body
{
DECLARE_SOA_INDEX_COLUMN(StoredKFVtx3BodyData, storedkfvtx3BodyData); //! Index to KFVtx3BodyData entry
}

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_VTX3BODYTABLES_H_
