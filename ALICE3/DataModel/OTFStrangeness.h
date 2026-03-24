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
/// \file   OTFStrangeness.h
/// \author David Dobrigkeit Chinellato
/// \since  05/08/2024
/// \brief  Set of tables for the ALICE3 strangeness information
///

#ifndef ALICE3_DATAMODEL_OTFSTRANGENESS_H_
#define ALICE3_DATAMODEL_OTFSTRANGENESS_H_

// O2 includes
#include "Common/Core/RecoDecay.h"

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace otfcascade
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                                     //!
DECLARE_SOA_INDEX_COLUMN_FULL(CascadeTrack, cascadeTrack, int, Tracks, "_Cascade"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, Tracks, "_Pos");             //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, Tracks, "_Neg");             //!
DECLARE_SOA_INDEX_COLUMN_FULL(BachTrack, bachTrack, int, Tracks, "_Bach");          //!

// topo vars
DECLARE_SOA_COLUMN(DcaV0Daughters, dcaV0Daughters, float);
DECLARE_SOA_COLUMN(DcaCascadeDaughters, dcaCascadeDaughters, float);
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);
DECLARE_SOA_COLUMN(CascRadius, cascRadius, float);
DECLARE_SOA_COLUMN(CascRadiusMC, cascRadiusMC, float);
DECLARE_SOA_COLUMN(MLambda, mLambda, float);
DECLARE_SOA_COLUMN(MXi, mXi, float);

// strangeness tracking
DECLARE_SOA_COLUMN(FindableClusters, findableClusters, int);
DECLARE_SOA_COLUMN(FoundClusters, foundClusters, int);

} // namespace otfcascade
DECLARE_SOA_TABLE(UpgradeCascades, "AOD", "UPGRADECASCADES",
                  o2::soa::Index<>,
                  otfcascade::CollisionId,
                  otfcascade::CascadeTrackId,
                  otfcascade::PosTrackId,
                  otfcascade::NegTrackId,
                  otfcascade::BachTrackId,
                  otfcascade::DcaV0Daughters,
                  otfcascade::DcaCascadeDaughters,
                  otfcascade::V0Radius,
                  otfcascade::CascRadius,
                  otfcascade::CascRadiusMC,
                  otfcascade::MLambda,
                  otfcascade::MXi,
                  otfcascade::FindableClusters,
                  otfcascade::FoundClusters);

using UpgradeCascade = UpgradeCascades::iterator;

namespace otfv0
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                         //!
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, Tracks, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, Tracks, "_Neg"); //!
DECLARE_SOA_INDEX_COLUMN(V0, v0);                                       //! index of the mc particle corresponding to the V0

// topo vars
DECLARE_SOA_COLUMN(DcaV0Daughters, dcaV0Daughters, float);
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);
DECLARE_SOA_COLUMN(MLambda, mLambda, float);
DECLARE_SOA_COLUMN(MAntiLambda, mAntiLambda, float);
DECLARE_SOA_COLUMN(MK0, mK0, float);

// kinematics
DECLARE_SOA_COLUMN(Pt, pt, float);

} // namespace otfv0
DECLARE_SOA_TABLE(UpgradeV0s, "AOD", "UPGRADEV0S",
                  o2::soa::Index<>,
                  otfv0::CollisionId,
                  otfv0::V0Id,
                  otfv0::PosTrackId,
                  otfv0::NegTrackId,
                  otfv0::DcaV0Daughters,
                  otfv0::V0Radius,
                  otfv0::MLambda,
                  otfv0::MAntiLambda,
                  otfv0::MK0,
                  otfv0::Pt);

using UpgradeV0 = UpgradeV0s::iterator;

namespace candidatev0
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                         //!
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, Tracks, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, Tracks, "_Neg"); //!
DECLARE_SOA_INDEX_COLUMN(V0, v0);                                       //!

// Label to MC particle
DECLARE_SOA_INDEX_COLUMN_FULL(McParticle, mcParticle, int, McParticles, ""); //! label to the MC particle corresponding to the V0

// General V0 properties: position, momentum
DECLARE_SOA_COLUMN(PosX, posX, float);   //! positive track X at min
DECLARE_SOA_COLUMN(NegX, negX, float);   //! negative track X at min
DECLARE_SOA_COLUMN(PxPos, pxPos, float); //! positive track px at min
DECLARE_SOA_COLUMN(PyPos, pyPos, float); //! positive track py at min
DECLARE_SOA_COLUMN(PzPos, pzPos, float); //! positive track pz at min
DECLARE_SOA_COLUMN(PxNeg, pxNeg, float); //! negative track px at min
DECLARE_SOA_COLUMN(PyNeg, pyNeg, float); //! negative track py at min
DECLARE_SOA_COLUMN(PzNeg, pzNeg, float); //! negative track pz at min
DECLARE_SOA_COLUMN(X, x, float);         //! decay position X
DECLARE_SOA_COLUMN(Y, y, float);         //! decay position Y
DECLARE_SOA_COLUMN(Z, z, float);         //! decay position Z

// topo vars
DECLARE_SOA_COLUMN(DcaV0Daughters, dcaV0Daughters, float);
DECLARE_SOA_COLUMN(CosPA, cosPA, float);
DECLARE_SOA_COLUMN(DcaPosToPV, dcaPosToPV, float);
DECLARE_SOA_COLUMN(DcaNegToPV, dcaNegToPV, float);
DECLARE_SOA_COLUMN(DcaV0ToPV, dcaV0ToPV, float);

//______________________________________________________
// DYNAMIC COLUMNS

DECLARE_SOA_DYNAMIC_COLUMN(Px, px, //! V0 px
                           [](float pxPos, float pxNeg) -> float { return pxPos + pxNeg; });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //! V0 py
                           [](float pyPos, float pyNeg) -> float { return pyPos + pyNeg; });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //! V0 pz
                           [](float pzPos, float pzNeg) -> float { return pzPos + pzNeg; });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //! Transverse momentum in GeV/c
                           [](float pxPos, float pyPos, float pxNeg, float pyNeg) -> float {
                             return RecoDecay::sqrtSumOfSquares(pxPos + pxNeg, pyPos + pyNeg);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //! Total momentum in GeV/c
                           [](float pxPos, float pyPos, float pzPos, float pxNeg, float pyNeg, float pzNeg) -> float {
                             return RecoDecay::sqrtSumOfSquares(pxPos + pxNeg, pyPos + pyNeg, pzPos + pzNeg);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, //! Phi in the range [0, 2pi)
                           [](float pxPos, float pyPos, float pxNeg, float pyNeg) -> float { return RecoDecay::phi(pxPos + pxNeg, pyPos + pyNeg); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, //! Pseudorapidity, conditionally defined to avoid FPEs
                           [](float pxPos, float pyPos, float pzPos, float pxNeg, float pyNeg, float pzNeg) -> float {
                             return RecoDecay::eta(std::array{pxPos + pxNeg, pyPos + pyNeg, pzPos + pzNeg});
                           });
// Length quantities
DECLARE_SOA_DYNAMIC_COLUMN(V0radius, v0radius, //! V0 decay radius (2D, centered at zero)
                           [](float x, float y) -> float { return RecoDecay::sqrtSumOfSquares(x, y); });

// Distance Over To Mom
DECLARE_SOA_DYNAMIC_COLUMN(DistOverTotMom, distOverTotMom, //! PV to V0decay distance over total momentum
                           [](float X, float Y, float Z, float pxPos, float pyPos, float pzPos, float pxNeg, float pyNeg, float pzNeg, float pvX, float pvY, float pvZ) {
                             float p = RecoDecay::sqrtSumOfSquares(pxPos + pxNeg, pyPos + pyNeg, pzPos + pzNeg);
                             return std::sqrt(std::pow(X - pvX, 2) + std::pow(Y - pvY, 2) + std::pow(Z - pvZ, 2)) / (p + 1E-10);
                           });

// Armenteros-Podolanski variables
DECLARE_SOA_DYNAMIC_COLUMN(Alpha, alpha, //! Armenteros Alpha
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             float momTot = RecoDecay::p(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
                             float lQlNeg = RecoDecay::dotProd(std::array{pxneg, pyneg, pzneg}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}) / momTot;
                             float lQlPos = RecoDecay::dotProd(std::array{pxpos, pypos, pzpos}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}) / momTot;
                             return (lQlPos - lQlNeg) / (lQlPos + lQlNeg); // alphav0
                           });

DECLARE_SOA_DYNAMIC_COLUMN(QtArm, qtArm, //! Armenteros Qt
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             float momTot = RecoDecay::p2(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
                             float dp = RecoDecay::dotProd(std::array{pxneg, pyneg, pzneg}, std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg});
                             return std::sqrt(RecoDecay::p2(pxneg, pyneg, pzneg) - dp * dp / momTot); // qtarm
                           });
// Mass assumption
DECLARE_SOA_DYNAMIC_COLUMN(MLambda, mLambda, //! mass under lambda hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float {
                             return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(MAntiLambda, mAntiLambda, //! mass under antilambda hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float {
                             return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(MK0Short, mK0Short, //! mass under K0short hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float {
                             return RecoDecay::m(std::array{std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged});
                           });
// Rapidity
DECLARE_SOA_DYNAMIC_COLUMN(YK0Short, yK0Short, //! V0 y with K0short hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float {
                             return RecoDecay::y(std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}, o2::constants::physics::MassKaonNeutral);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(YLambda, yLambda, //! V0 y with lambda or antilambda hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float {
                             return RecoDecay::y(std::array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}, o2::constants::physics::MassLambda);
                           });
// Daughter track momenta
DECLARE_SOA_DYNAMIC_COLUMN(NegativePt, negativePt, //! negative daughter pT
                           [](float pxneg, float pyneg) -> float { return RecoDecay::sqrtSumOfSquares(pxneg, pyneg); });
DECLARE_SOA_DYNAMIC_COLUMN(PositivePt, positivePt, //! positive daughter pT
                           [](float pxpos, float pypos) -> float { return RecoDecay::sqrtSumOfSquares(pxpos, pypos); });
DECLARE_SOA_DYNAMIC_COLUMN(NegativeEta, negativeEta, //! negative daughter eta
                           [](float PxNeg, float PyNeg, float PzNeg) -> float { return RecoDecay::eta(std::array{PxNeg, PyNeg, PzNeg}); });
DECLARE_SOA_DYNAMIC_COLUMN(NegativePhi, negativePhi, //! negative daughter phi
                           [](float PxNeg, float PyNeg) -> float { return RecoDecay::phi(PxNeg, PyNeg); });
DECLARE_SOA_DYNAMIC_COLUMN(PositiveEta, positiveEta, //! positive daughter eta
                           [](float PxPos, float PyPos, float PzPos) -> float { return RecoDecay::eta(std::array{PxPos, PyPos, PzPos}); });
DECLARE_SOA_DYNAMIC_COLUMN(PositivePhi, positivePhi, //! positive daughter phi
                           [](float PxPos, float PyPos) -> float { return RecoDecay::phi(PxPos, PyPos); });
} // namespace candidatev0
DECLARE_SOA_TABLE(V0CandidateIndices, "AOD", "V0CANDIDATEINDEX", //! index table
                  o2::soa::Index<>, candidatev0::CollisionId, candidatev0::PosTrackId, candidatev0::NegTrackId, candidatev0::McParticleId);

DECLARE_SOA_TABLE(V0CandidateCores, "AOD", "V0CANDIDATECORE",
                  o2::soa::Index<>,
                  candidatev0::X, candidatev0::Y, candidatev0::Z,
                  candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos,
                  candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg,
                  candidatev0::DcaV0Daughters, candidatev0::DcaPosToPV, candidatev0::DcaNegToPV,
                  candidatev0::CosPA, candidatev0::DcaV0ToPV,
                  candidatev0::Px<candidatev0::PxPos, candidatev0::PxNeg>,
                  candidatev0::Py<candidatev0::PyPos, candidatev0::PyNeg>,
                  candidatev0::Pz<candidatev0::PzPos, candidatev0::PzNeg>,
                  candidatev0::Pt<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PxNeg, candidatev0::PyNeg>,
                  candidatev0::P<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos, candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::Phi<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PxNeg, candidatev0::PyNeg>,
                  candidatev0::Eta<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos, candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::V0radius<candidatev0::X, candidatev0::Y>,
                  candidatev0::DistOverTotMom<candidatev0::X, candidatev0::Y, candidatev0::Z, candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos, candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::Alpha<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos, candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::QtArm<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos, candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::MLambda<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos, candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::MAntiLambda<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos, candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::MK0Short<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos, candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::YK0Short<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos, candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::YLambda<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos, candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::NegativePt<candidatev0::PxNeg, candidatev0::PyNeg>,
                  candidatev0::PositivePt<candidatev0::PxPos, candidatev0::PyPos>,
                  candidatev0::NegativeEta<candidatev0::PxNeg, candidatev0::PyNeg, candidatev0::PzNeg>,
                  candidatev0::NegativePhi<candidatev0::PxNeg, candidatev0::PyNeg>,
                  candidatev0::PositiveEta<candidatev0::PxPos, candidatev0::PyPos, candidatev0::PzPos>,
                  candidatev0::PositivePhi<candidatev0::PxPos, candidatev0::PyPos>);

using V0CandidateCore = V0CandidateCores::iterator;
} // namespace o2::aod
#endif // ALICE3_DATAMODEL_OTFSTRANGENESS_H_
