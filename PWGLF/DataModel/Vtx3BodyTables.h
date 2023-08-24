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
DECLARE_SOA_COLUMN(DCAVtxDaughters, dcaVtxdaughters, float); //! DCA between  daughters
DECLARE_SOA_COLUMN(DCATrack0ToPV, dcatrack0topv, float);     //! DCA prong0 to PV
DECLARE_SOA_COLUMN(DCATrack1ToPV, dcatrack1topv, float);     //! DCA prong1 to PV
DECLARE_SOA_COLUMN(DCATrack2ToPV, dcatrack2topv, float);     //! DCA prong2 to PV

// Derived expressions
// Momenta
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
DECLARE_SOA_DYNAMIC_COLUMN(DCAVtxToPV, dcavtxtopv, //! DCA of 3 body vtx to PV
                           [](float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ) -> float { return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz)); });

// Calculated on the fly with mass assumption + dynamic tables
DECLARE_SOA_DYNAMIC_COLUMN(MHypertriton, mHypertriton, //! mass under Hypertriton hypothesis
                           [](float pxtrack0, float pytrack0, float pztrack0, float pxtrack1, float pytrack1, float pztrack1, float pxtrack2, float pytrack2, float pztrack2) -> float { return RecoDecay::m(std::array{std::array{pxtrack0, pytrack0, pztrack0}, std::array{pxtrack1, pytrack1, pztrack1}, std::array{pxtrack2, pytrack2, pztrack2}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged, o2::constants::physics::MassDeuteron}); });
DECLARE_SOA_DYNAMIC_COLUMN(MAntiHypertriton, mAntiHypertriton, //! mass under antiHypertriton hypothesis
                           [](float pxtrack0, float pytrack0, float pztrack0, float pxtrack1, float pytrack1, float pztrack1, float pxtrack2, float pytrack2, float pztrack2) -> float { return RecoDecay::m(std::array{std::array{pxtrack0, pytrack0, pztrack0}, std::array{pxtrack1, pytrack1, pztrack1}, std::array{pxtrack2, pytrack2, pztrack2}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron}); });

DECLARE_SOA_DYNAMIC_COLUMN(YHypertriton, yHypertriton,                                                                                                           //! 3 body vtx y with hypertriton or antihypertriton hypothesis
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::y(std::array{Px, Py, Pz}, o2::constants::physics::MassHyperTriton); }); // here MassHyperTriton = 2.992
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta,                                                                                                                             //! 3 body vtx eta
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::eta(std::array{Px, Py, Pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, //! 3 body vtx phi
                           [](float Px, float Py) -> float { return RecoDecay::phi(Px, Py); });

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
                       vtx3body::DCATrack0ToPV, vtx3body::DCATrack1ToPV, vtx3body::DCATrack2ToPV,

                       // Dynamic columns
                       vtx3body::Pt<vtx3body::PxTrack0, vtx3body::PyTrack0, vtx3body::PxTrack1, vtx3body::PyTrack1, vtx3body::PxTrack2, vtx3body::PyTrack2>,
                       vtx3body::VtxRadius<vtx3body::X, vtx3body::Y>,
                       vtx3body::DistOverTotMom<vtx3body::X, vtx3body::Y, vtx3body::Z, vtx3body::Px, vtx3body::Py, vtx3body::Pz>,
                       vtx3body::VtxCosPA<vtx3body::X, vtx3body::Y, vtx3body::Z, vtx3body::Px, vtx3body::Py, vtx3body::Pz>,
                       vtx3body::DCAVtxToPV<vtx3body::X, vtx3body::Y, vtx3body::Z, vtx3body::Px, vtx3body::Py, vtx3body::Pz>,

                       // Invariant masses
                       vtx3body::MHypertriton<vtx3body::PxTrack0, vtx3body::PyTrack0, vtx3body::PzTrack0, vtx3body::PxTrack1, vtx3body::PyTrack1, vtx3body::PzTrack1, vtx3body::PxTrack2, vtx3body::PyTrack2, vtx3body::PzTrack2>,
                       vtx3body::MAntiHypertriton<vtx3body::PxTrack0, vtx3body::PyTrack0, vtx3body::PzTrack0, vtx3body::PxTrack1, vtx3body::PyTrack1, vtx3body::PzTrack1, vtx3body::PxTrack2, vtx3body::PyTrack2, vtx3body::PzTrack2>,

                       // Longitudinal
                       vtx3body::YHypertriton<vtx3body::Px, vtx3body::Py, vtx3body::Pz>,
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

} // namespace o2::aod
#endif // PWGLF_DATAMODEL_VTX3BODYTABLES_H_
