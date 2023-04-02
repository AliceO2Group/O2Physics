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
#ifndef O2_ANALYSIS_STRANGENESSTABLES_H_
#define O2_ANALYSIS_STRANGENESSTABLES_H_

#include "Framework/AnalysisDataModel.h"
#include "Common/Core/RecoDecay.h"
#include "CommonConstants/PhysicsConstants.h"
#include <cmath>

namespace o2::aod
{
namespace v0data
{
// Needed to have shorter table that does not rely on existing one (filtering!)
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, Tracks, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, Tracks, "_Neg"); //!
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                         //!
DECLARE_SOA_INDEX_COLUMN(V0, v0);                                       //!

// General V0 properties: position, momentum
DECLARE_SOA_COLUMN(PosX, posX, float);   //! positive track X at min
DECLARE_SOA_COLUMN(NegX, negX, float);   //! negative track X at min
DECLARE_SOA_COLUMN(PxPos, pxpos, float); //! positive track px at min
DECLARE_SOA_COLUMN(PyPos, pypos, float); //! positive track py at min
DECLARE_SOA_COLUMN(PzPos, pzpos, float); //! positive track pz at min
DECLARE_SOA_COLUMN(PxNeg, pxneg, float); //! negative track px at min
DECLARE_SOA_COLUMN(PyNeg, pyneg, float); //! negative track py at min
DECLARE_SOA_COLUMN(PzNeg, pzneg, float); //! negative track pz at min
DECLARE_SOA_COLUMN(X, x, float);         //! decay position X
DECLARE_SOA_COLUMN(Y, y, float);         //! decay position Y
DECLARE_SOA_COLUMN(Z, z, float);         //! decay position Z

// Saved from finding: DCAs
DECLARE_SOA_COLUMN(DCAV0Daughters, dcaV0daughters, float); //! DCA between V0 daughters
DECLARE_SOA_COLUMN(DCAPosToPV, dcapostopv, float);         //! DCA positive prong to PV
DECLARE_SOA_COLUMN(DCANegToPV, dcanegtopv, float);         //! DCA negative prong to PV

// Saved from finding: covariance matrix of parent track (on request)
DECLARE_SOA_COLUMN(PositionCovMat, positionCovMat, float[6]); //! covariance matrix elements
DECLARE_SOA_COLUMN(MomentumCovMat, momentumCovMat, float[6]); //! covariance matrix elements

// Derived expressions
// Momenta
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //! V0 pT
                           [](float pxpos, float pypos, float pxneg, float pyneg) -> float { return RecoDecay::sqrtSumOfSquares(pxpos + pxneg, pypos + pyneg); });
// Account for rigidity in case of hypertriton
DECLARE_SOA_DYNAMIC_COLUMN(PtHypertriton, ptHypertriton, //! V0 pT
                           [](float pxpos, float pypos, float pxneg, float pyneg) -> float { return RecoDecay::sqrtSumOfSquares(2.0f * pxpos + pxneg, 2.0f * pypos + pyneg); });
DECLARE_SOA_DYNAMIC_COLUMN(PtAntiHypertriton, ptAntiHypertriton, //! V0 pT
                           [](float pxpos, float pypos, float pxneg, float pyneg) -> float { return RecoDecay::sqrtSumOfSquares(pxpos + 2.0f * pxneg, pypos + 2.0f * pyneg); });

// Length quantities
DECLARE_SOA_DYNAMIC_COLUMN(V0Radius, v0radius, //! V0 decay radius (2D, centered at zero)
                           [](float x, float y) -> float { return RecoDecay::sqrtSumOfSquares(x, y); });

// Distance Over To Mom
DECLARE_SOA_DYNAMIC_COLUMN(DistOverTotMom, distovertotmom, //! PV to V0decay distance over total momentum
                           [](float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ) {
                             float P = RecoDecay::sqrtSumOfSquares(Px, Py, Pz);
                             return std::sqrt(std::pow(X - pvX, 2) + std::pow(Y - pvY, 2) + std::pow(Z - pvZ, 2)) / (P + 1E-10);
                           });

// CosPA
DECLARE_SOA_DYNAMIC_COLUMN(V0CosPA, v0cosPA, //! V0 CosPA
                           [](float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ) -> float { return RecoDecay::cpa(array{pvX, pvY, pvZ}, array{X, Y, Z}, array{Px, Py, Pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(DCAV0ToPV, dcav0topv, //! DCA of V0 to PV
                           [](float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ) -> float { return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz)); });

// Armenteros-Podolanski variables
DECLARE_SOA_DYNAMIC_COLUMN(Alpha, alpha, //! Armenteros Alpha
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             float momTot = RecoDecay::p(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
                             float lQlNeg = RecoDecay::dotProd(array{pxneg, pyneg, pzneg}, array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}) / momTot;
                             float lQlPos = RecoDecay::dotProd(array{pxpos, pypos, pzpos}, array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}) / momTot;
                             return (lQlPos - lQlNeg) / (lQlPos + lQlNeg); // alphav0
                           });

DECLARE_SOA_DYNAMIC_COLUMN(QtArm, qtarm, //! Armenteros Qt
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             float momTot = RecoDecay::p2(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
                             float dp = RecoDecay::dotProd(array{pxneg, pyneg, pzneg}, array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg});
                             return std::sqrt(RecoDecay::p2(pxneg, pyneg, pzneg) - dp * dp / momTot); // qtarm
                           });

// Psi pair angle: angle between the plane defined by the electron and positron momenta and the xy plane
DECLARE_SOA_DYNAMIC_COLUMN(PsiPair, psipair, //! psi pair angle
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             auto clipToPM1 = [](float x) { return x < -1.f ? -1.f : (x > 1.f ? 1.f : x); };
                             float ptot2 = RecoDecay::p2(pxpos, pypos, pzpos) * RecoDecay::p2(pxneg, pyneg, pzneg);
                             float argcos = RecoDecay::dotProd(array{pxpos, pypos, pzpos}, array{pxneg, pyneg, pzneg}) / std::sqrt(ptot2);
                             float thetaPos = std::atan2(RecoDecay::sqrtSumOfSquares(pxpos, pypos), pzpos);
                             float thetaNeg = std::atan2(RecoDecay::sqrtSumOfSquares(pxneg, pyneg), pzneg);
                             float argsin = (thetaNeg - thetaPos) / std::acos(clipToPM1(argcos));
                             return std::asin(clipToPM1(argsin));
                           });

// calculate the fraction of the pos/neg momentum of the V0 momentum
DECLARE_SOA_DYNAMIC_COLUMN(PFracPos, pfracpos,
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             float ppos = RecoDecay::sqrtSumOfSquares(pxpos, pypos, pzpos);
                             float PV0 = RecoDecay::sqrtSumOfSquares(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
                             return (ppos / PV0);
                           });

DECLARE_SOA_DYNAMIC_COLUMN(PFracNeg, pfracneg,
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             float pneg = RecoDecay::sqrtSumOfSquares(pxneg, pyneg, pzneg);
                             float PV0 = RecoDecay::sqrtSumOfSquares(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
                             return (pneg / PV0);
                           });

// Calculated on the fly with mass assumption + dynamic tables
DECLARE_SOA_DYNAMIC_COLUMN(MLambda, mLambda, //! mass under lambda hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(array{array{pxpos, pypos, pzpos}, array{pxneg, pyneg, pzneg}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged}); });
DECLARE_SOA_DYNAMIC_COLUMN(MAntiLambda, mAntiLambda, //! mass under antilambda hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(array{array{pxpos, pypos, pzpos}, array{pxneg, pyneg, pzneg}}, array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton}); });
DECLARE_SOA_DYNAMIC_COLUMN(MK0Short, mK0Short, //! mass under K0short hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(array{array{pxpos, pypos, pzpos}, array{pxneg, pyneg, pzneg}}, array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged}); });
DECLARE_SOA_DYNAMIC_COLUMN(MGamma, mGamma, //! mass under gamma hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(array{array{pxpos, pypos, pzpos}, array{pxneg, pyneg, pzneg}}, array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron}); });
// Account for rigidity in case of hypertriton
DECLARE_SOA_DYNAMIC_COLUMN(MHypertriton, mHypertriton, //! mass under hypertriton  hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(array{array{2.0f * pxpos, 2.0f * pypos, 2.0f * pzpos}, array{pxneg, pyneg, pzneg}}, array{o2::constants::physics::MassHelium3, o2::constants::physics::MassPionCharged}); });
DECLARE_SOA_DYNAMIC_COLUMN(MAntiHypertriton, mAntiHypertriton, //! mass under antihypertriton hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(array{array{pxpos, pypos, pzpos}, array{2.0f * pxneg, 2.0f * pyneg, 2.0f * pzneg}}, array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassHelium3}); });

DECLARE_SOA_DYNAMIC_COLUMN(YK0Short, yK0Short, //! V0 y with K0short hypothesis
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::y(array{Px, Py, Pz}, o2::constants::physics::MassKaonNeutral); });
DECLARE_SOA_DYNAMIC_COLUMN(YLambda, yLambda, //! V0 y with lambda or antilambda hypothesis
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::y(array{Px, Py, Pz}, o2::constants::physics::MassLambda); });
DECLARE_SOA_DYNAMIC_COLUMN(YHypertriton, yHypertriton, //! V0 y with hypertriton hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::y(array{2.0f * pxpos + pxneg, 2.0f * pypos + pyneg, 2.0f * pzpos + pzneg}, o2::constants::physics::MassHyperTriton); });
DECLARE_SOA_DYNAMIC_COLUMN(YAntiHypertriton, yAntiHypertriton, //! V0 y with antihypertriton hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::y(array{pxpos + 2.0f * pxneg, pypos + 2.0f * pyneg, pzpos + 2.0f * pzneg}, o2::constants::physics::MassHyperTriton); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, //! V0 eta
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::eta(array{Px, Py, Pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, //! V0 phi
                           [](float Px, float Py) -> float { return RecoDecay::phi(Px, Py); });

DECLARE_SOA_DYNAMIC_COLUMN(NegativePt, negativept, //! negative daughter pT
                           [](float pxneg, float pyneg) -> float { return RecoDecay::sqrtSumOfSquares(pxneg, pyneg); });
DECLARE_SOA_DYNAMIC_COLUMN(PositivePt, positivept, //! positive daughter pT
                           [](float pxpos, float pypos) -> float { return RecoDecay::sqrtSumOfSquares(pxpos, pypos); });
DECLARE_SOA_DYNAMIC_COLUMN(NegativeEta, negativeeta, //! negative daughter eta
                           [](float PxNeg, float PyNeg, float PzNeg) -> float { return RecoDecay::eta(array{PxNeg, PyNeg, PzNeg}); });
DECLARE_SOA_DYNAMIC_COLUMN(NegativePhi, negativephi, //! negative daughter phi
                           [](float PxNeg, float PyNeg) -> float { return RecoDecay::phi(PxNeg, PyNeg); });
DECLARE_SOA_DYNAMIC_COLUMN(PositiveEta, positiveeta, //! positive daughter eta
                           [](float PxPos, float PyPos, float PzPos) -> float { return RecoDecay::eta(array{PxPos, PyPos, PzPos}); });
DECLARE_SOA_DYNAMIC_COLUMN(PositivePhi, positivephi, //! positive daughter phi
                           [](float PxPos, float PyPos) -> float { return RecoDecay::phi(PxPos, PyPos); });

DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //! V0 px
                              float, 1.f * aod::v0data::pxpos + 1.f * aod::v0data::pxneg);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //! V0 py
                              float, 1.f * aod::v0data::pypos + 1.f * aod::v0data::pyneg);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //! V0 pz
                              float, 1.f * aod::v0data::pzpos + 1.f * aod::v0data::pzneg);
} // namespace v0data

DECLARE_SOA_TABLE_FULL(StoredV0Datas, "V0Datas", "AOD", "V0DATA", //!
                       o2::soa::Index<>, v0data::PosTrackId, v0data::NegTrackId, v0data::CollisionId, v0data::V0Id,
                       v0data::PosX, v0data::NegX,
                       v0data::X, v0data::Y, v0data::Z,
                       v0data::PxPos, v0data::PyPos, v0data::PzPos,
                       v0data::PxNeg, v0data::PyNeg, v0data::PzNeg,
                       v0data::DCAV0Daughters, v0data::DCAPosToPV, v0data::DCANegToPV,

                       // Dynamic columns
                       v0data::Pt<v0data::PxPos, v0data::PyPos, v0data::PxNeg, v0data::PyNeg>,
                       v0data::PtHypertriton<v0data::PxPos, v0data::PyPos, v0data::PxNeg, v0data::PyNeg>,
                       v0data::PtAntiHypertriton<v0data::PxPos, v0data::PyPos, v0data::PxNeg, v0data::PyNeg>,
                       v0data::V0Radius<v0data::X, v0data::Y>,
                       v0data::DistOverTotMom<v0data::X, v0data::Y, v0data::Z, v0data::Px, v0data::Py, v0data::Pz>,
                       v0data::V0CosPA<v0data::X, v0data::Y, v0data::Z, v0data::Px, v0data::Py, v0data::Pz>,
                       v0data::DCAV0ToPV<v0data::X, v0data::Y, v0data::Z, v0data::Px, v0data::Py, v0data::Pz>,
                       v0data::Alpha<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::QtArm<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::PsiPair<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::PFracPos<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::PFracNeg<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,

                       // Invariant masses
                       v0data::MLambda<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::MAntiLambda<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::MK0Short<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::MGamma<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::MHypertriton<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::MAntiHypertriton<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,

                       // Longitudinal
                       v0data::YK0Short<v0data::Px, v0data::Py, v0data::Pz>,
                       v0data::YLambda<v0data::Px, v0data::Py, v0data::Pz>,
                       v0data::YHypertriton<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::YAntiHypertriton<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::Eta<v0data::Px, v0data::Py, v0data::Pz>,
                       v0data::Phi<v0data::Px, v0data::Py>,
                       v0data::NegativePt<v0data::PxNeg, v0data::PyNeg>,
                       v0data::PositivePt<v0data::PxPos, v0data::PyPos>,
                       v0data::NegativeEta<v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                       v0data::NegativePhi<v0data::PxNeg, v0data::PyNeg>,
                       v0data::PositiveEta<v0data::PxPos, v0data::PyPos, v0data::PzPos>,
                       v0data::PositivePhi<v0data::PxPos, v0data::PyPos>);

DECLARE_SOA_TABLE_FULL(V0Covs, "V0Covs", "AOD", "V0COVS", //!
                       v0data::PositionCovMat, v0data::MomentumCovMat);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(V0Datas, StoredV0Datas, "V0DATAEXT", //!
                                v0data::Px, v0data::Py, v0data::Pz); // the table name has here to be the one with EXT which is not nice and under study

using V0Data = V0Datas::iterator;
namespace v0data
{
DECLARE_SOA_INDEX_COLUMN(V0Data, v0Data); //! Index to V0Data entry
}

DECLARE_SOA_TABLE(V0DataLink, "AOD", "V0DATALINK", //! Joinable table with V0s which links to V0Data which is not produced for all entries
                  v0data::V0DataId);

using V0sLinked = soa::Join<V0s, V0DataLink>;
using V0Linked = V0sLinked::iterator;

// helper for building
namespace v0tag
{
// Global bool
DECLARE_SOA_COLUMN(IsInteresting, isInteresting, bool); //! will this be built or not?

// MC association bools
DECLARE_SOA_COLUMN(IsTrueGamma, isTrueGamma, bool);                     //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueK0Short, isTrueK0Short, bool);                 //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueLambda, isTrueLambda, bool);                   //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueAntiLambda, isTrueAntiLambda, bool);           //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueHypertriton, isTrueHypertriton, bool);         //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueAntiHypertriton, isTrueAntiHypertriton, bool); //! PDG checked correctly in MC

// dE/dx compatibility bools
DECLARE_SOA_COLUMN(IsGammaCandidate, isGammaCandidate, bool);                     //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsK0ShortCandidate, isK0ShortCandidate, bool);                 //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsLambdaCandidate, isLambdaCandidate, bool);                   //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsAntiLambdaCandidate, isAntiLambdaCandidate, bool);           //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsHypertritonCandidate, isHypertritonCandidate, bool);         //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsAntiHypertritonCandidate, isAntiHypertritonCandidate, bool); //! compatible with dE/dx hypotheses
} // namespace v0tag
DECLARE_SOA_TABLE(V0Tags, "AOD", "V0TAGS",
                  v0tag::IsInteresting,
                  v0tag::IsTrueGamma,
                  v0tag::IsTrueK0Short,
                  v0tag::IsTrueLambda,
                  v0tag::IsTrueAntiLambda,
                  v0tag::IsTrueHypertriton,
                  v0tag::IsTrueAntiHypertriton,
                  v0tag::IsGammaCandidate,
                  v0tag::IsK0ShortCandidate,
                  v0tag::IsLambdaCandidate,
                  v0tag::IsAntiLambdaCandidate,
                  v0tag::IsHypertritonCandidate,
                  v0tag::IsAntiHypertritonCandidate);

namespace cascdata
{
// Necessary for full filtering functionality
DECLARE_SOA_INDEX_COLUMN(V0, v0);                                   //!
DECLARE_SOA_INDEX_COLUMN(Cascade, cascade);                         //!
DECLARE_SOA_INDEX_COLUMN_FULL(Bachelor, bachelor, int, Tracks, ""); //!
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                     //!
// General cascade properties: position, momentum
DECLARE_SOA_COLUMN(Sign, sign, int);         //!
DECLARE_SOA_COLUMN(PxPos, pxpos, float);     //!
DECLARE_SOA_COLUMN(PyPos, pypos, float);     //!
DECLARE_SOA_COLUMN(PzPos, pzpos, float);     //!
DECLARE_SOA_COLUMN(PxNeg, pxneg, float);     //!
DECLARE_SOA_COLUMN(PyNeg, pyneg, float);     //!
DECLARE_SOA_COLUMN(PzNeg, pzneg, float);     //!
DECLARE_SOA_COLUMN(PxBach, pxbach, float);   //!
DECLARE_SOA_COLUMN(PyBach, pybach, float);   //!
DECLARE_SOA_COLUMN(PzBach, pzbach, float);   //!
DECLARE_SOA_COLUMN(X, x, float);             //!
DECLARE_SOA_COLUMN(Y, y, float);             //!
DECLARE_SOA_COLUMN(Z, z, float);             //!
DECLARE_SOA_COLUMN(Xlambda, xlambda, float); //!
DECLARE_SOA_COLUMN(Ylambda, ylambda, float); //!
DECLARE_SOA_COLUMN(Zlambda, zlambda, float); //!

// Saved from finding: DCAs
DECLARE_SOA_COLUMN(DCAV0Daughters, dcaV0daughters, float);     //!
DECLARE_SOA_COLUMN(DCACascDaughters, dcacascdaughters, float); //!
DECLARE_SOA_COLUMN(DCAPosToPV, dcapostopv, float);             //!
DECLARE_SOA_COLUMN(DCANegToPV, dcanegtopv, float);             //!
DECLARE_SOA_COLUMN(DCABachToPV, dcabachtopv, float);           //!
DECLARE_SOA_COLUMN(DCACascToPV, dcacasctopv, float);           //!

// Saved from finding: covariance matrix of parent track (on request)
DECLARE_SOA_COLUMN(PositionCovMat, positionCovMat, float[6]); //! covariance matrix elements
DECLARE_SOA_COLUMN(MomentumCovMat, momentumCovMat, float[6]); //! covariance matrix elements

// Derived expressions
// Momenta
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //!
                           [](float Px, float Py) -> float { return RecoDecay::sqrtSumOfSquares(Px, Py); });

// Length quantities
DECLARE_SOA_DYNAMIC_COLUMN(V0Radius, v0radius, //!
                           [](float xlambda, float ylambda) -> float { return RecoDecay::sqrtSumOfSquares(xlambda, ylambda); });
DECLARE_SOA_DYNAMIC_COLUMN(CascRadius, cascradius, //!
                           [](float x, float y) -> float { return RecoDecay::sqrtSumOfSquares(x, y); });

// CosPAs
DECLARE_SOA_DYNAMIC_COLUMN(V0CosPA, v0cosPA, //!
                           [](float Xlambda, float Ylambda, float Zlambda, float PxLambda, float PyLambda, float PzLambda, float pvX, float pvY, float pvZ) -> float { return RecoDecay::cpa(array{pvX, pvY, pvZ}, array{Xlambda, Ylambda, Zlambda}, array{PxLambda, PyLambda, PzLambda}); });
DECLARE_SOA_DYNAMIC_COLUMN(CascCosPA, casccosPA, //!
                           [](float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ) -> float { return RecoDecay::cpa(array{pvX, pvY, pvZ}, array{X, Y, Z}, array{Px, Py, Pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(DCAV0ToPV, dcav0topv, //!
                           [](float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ) -> float { return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz)); });

// Calculated on the fly with mass assumption + dynamic tables
DECLARE_SOA_DYNAMIC_COLUMN(MLambda, mLambda, //!
                           [](int charge, float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::m(array{array{pxpos, pypos, pzpos}, array{pxneg, pyneg, pzneg}}, charge < 0 ? array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged} : array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton}); });
// Calculated on the fly with mass assumption + dynamic tables

DECLARE_SOA_DYNAMIC_COLUMN(MXi, mXi, //!
                           [](float pxbach, float pybach, float pzbach, float PxLambda, float PyLambda, float PzLambda) -> float { return RecoDecay::m(array{array{pxbach, pybach, pzbach}, array{PxLambda, PyLambda, PzLambda}}, array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassLambda}); });
DECLARE_SOA_DYNAMIC_COLUMN(MOmega, mOmega, //!
                           [](float pxbach, float pybach, float pzbach, float PxLambda, float PyLambda, float PzLambda) -> float { return RecoDecay::m(array{array{pxbach, pybach, pzbach}, array{PxLambda, PyLambda, PzLambda}}, array{o2::constants::physics::MassKaonCharged, o2::constants::physics::MassLambda}); });

DECLARE_SOA_DYNAMIC_COLUMN(YXi, yXi, //!
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::y(array{Px, Py, Pz}, o2::constants::physics::MassXiMinus); });
DECLARE_SOA_DYNAMIC_COLUMN(YOmega, yOmega, //!
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::y(array{Px, Py, Pz}, o2::constants::physics::MassOmegaMinus); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, //!
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::eta(array{Px, Py, Pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, //! cascade phi
                           [](float Px, float Py) -> float { return RecoDecay::phi(Px, Py); });
} // namespace cascdata

namespace cascdataext
{
DECLARE_SOA_EXPRESSION_COLUMN(PxLambda, pxlambda, //!
                              float, 1.f * aod::cascdata::pxpos + 1.f * aod::cascdata::pxneg);
DECLARE_SOA_EXPRESSION_COLUMN(PyLambda, pylambda, //!
                              float, 1.f * aod::cascdata::pypos + 1.f * aod::cascdata::pyneg);
DECLARE_SOA_EXPRESSION_COLUMN(PzLambda, pzlambda, //!
                              float, 1.f * aod::cascdata::pzpos + 1.f * aod::cascdata::pzneg);
DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //!
                              float, 1.f * aod::cascdata::pxpos + 1.f * aod::cascdata::pxneg + 1.f * aod::cascdata::pxbach);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //!
                              float, 1.f * aod::cascdata::pypos + 1.f * aod::cascdata::pyneg + 1.f * aod::cascdata::pybach);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //!
                              float, 1.f * aod::cascdata::pzpos + 1.f * aod::cascdata::pzneg + 1.f * aod::cascdata::pzbach);
} // namespace cascdataext

DECLARE_SOA_TABLE(StoredCascDatas, "AOD", "CASCDATA", //!
                  o2::soa::Index<>, cascdata::V0Id, cascdata::CascadeId, cascdata::BachelorId, cascdata::CollisionId,

                  cascdata::Sign,
                  cascdata::X, cascdata::Y, cascdata::Z,
                  cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda,
                  cascdata::PxPos, cascdata::PyPos, cascdata::PzPos,
                  cascdata::PxNeg, cascdata::PyNeg, cascdata::PzNeg,
                  cascdata::PxBach, cascdata::PyBach, cascdata::PzBach,
                  cascdata::DCAV0Daughters, cascdata::DCACascDaughters,
                  cascdata::DCAPosToPV, cascdata::DCANegToPV, cascdata::DCABachToPV, cascdata::DCACascToPV,

                  // Dynamic columns
                  cascdata::Pt<cascdataext::Px, cascdataext::Py>,
                  cascdata::V0Radius<cascdata::Xlambda, cascdata::Ylambda>,
                  cascdata::CascRadius<cascdata::X, cascdata::Y>,
                  cascdata::V0CosPA<cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda, cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda>,
                  cascdata::CascCosPA<cascdata::X, cascdata::Y, cascdata::Z, cascdataext::Px, cascdataext::Py, cascdataext::Pz>,
                  cascdata::DCAV0ToPV<cascdata::Xlambda, cascdata::Ylambda, cascdata::Zlambda, cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda>,

                  // Invariant masses
                  cascdata::MLambda<cascdata::Sign, cascdata::PxPos, cascdata::PyPos, cascdata::PzPos, cascdata::PxNeg, cascdata::PyNeg, cascdata::PzNeg>,
                  cascdata::MXi<cascdata::PxBach, cascdata::PyBach, cascdata::PzBach, cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda>,
                  cascdata::MOmega<cascdata::PxBach, cascdata::PyBach, cascdata::PzBach, cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda>,
                  // Longitudinal
                  cascdata::YXi<cascdataext::Px, cascdataext::Py, cascdataext::Pz>,
                  cascdata::YOmega<cascdataext::Px, cascdataext::Py, cascdataext::Pz>,
                  cascdata::Eta<cascdataext::Px, cascdataext::Py, cascdataext::Pz>,
                  cascdata::Phi<cascdataext::Px, cascdataext::Py>);

DECLARE_SOA_TABLE_FULL(CascCovs, "CascCovs", "AOD", "CASCCOVS", //!
                       cascdata::PositionCovMat, cascdata::MomentumCovMat);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(CascDatas, StoredCascDatas, "CascDATAEXT", //!
                                cascdataext::PxLambda, cascdataext::PyLambda, cascdataext::PzLambda,
                                cascdataext::Px, cascdataext::Py, cascdataext::Pz);

using CascData = CascDatas::iterator;

// For compatibility with previous table declarations
using CascDataFull = CascDatas;
using CascDataExt = CascDatas;

namespace cascdata
{
DECLARE_SOA_INDEX_COLUMN(CascData, cascData); //! Index to CascData entry
}

DECLARE_SOA_TABLE(CascDataLink, "AOD", "CASCDATALINK", //! Joinable table with Cascades which links to CascData which is not produced for all entries
                  cascdata::CascDataId);

using CascadesLinked = soa::Join<Cascades, CascDataLink>;
using CascadeLinked = CascadesLinked::iterator;

namespace casctag
{
DECLARE_SOA_COLUMN(IsInteresting, isInteresting, bool); //! will this be built or not?

// MC association bools
DECLARE_SOA_COLUMN(IsTrueXiMinus, isTrueXiMinus, bool);       //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueXiPlus, isTrueXiPlus, bool);         //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueOmegaMinus, isTrueOmegaMinus, bool); //! PDG checked correctly in MC
DECLARE_SOA_COLUMN(IsTrueOmegaPlus, isTrueOmegaPlus, bool);   //! PDG checked correctly in MC

// dE/dx compatibility bools
DECLARE_SOA_COLUMN(IsXiMinusCandidate, isXiMinusCandidate, bool);       //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsXiPlusCandidate, isXiPlusCandidate, bool);         //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsOmegaMinusCandidate, isOmegaMinusCandidate, bool); //! compatible with dE/dx hypotheses
DECLARE_SOA_COLUMN(IsOmegaPlusCandidate, isOmegaPlusCandidate, bool);   //! compatible with dE/dx hypotheses
} // namespace casctag
DECLARE_SOA_TABLE(CascTags, "AOD", "CASCTAGS",
                  casctag::IsInteresting,
                  casctag::IsTrueXiMinus,
                  casctag::IsTrueXiPlus,
                  casctag::IsTrueOmegaMinus,
                  casctag::IsTrueOmegaPlus,
                  casctag::IsXiMinusCandidate,
                  casctag::IsXiPlusCandidate,
                  casctag::IsOmegaMinusCandidate,
                  casctag::IsOmegaPlusCandidate);

// Definition of labels for V0s
namespace mcv0label
{
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle); //! MC particle for V0
} // namespace mcv0label

DECLARE_SOA_TABLE(McV0Labels, "AOD", "MCV0LABEL", //! Table joinable to V0data containing the MC labels
                  mcv0label::McParticleId);
using McV0Label = McV0Labels::iterator;

// Definition of labels for cascades
namespace mccasclabel
{
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle); //! MC particle for V0
} // namespace mccasclabel

DECLARE_SOA_TABLE(McCascLabels, "AOD", "MCCASCLABEL", //! Table joinable to V0data containing the MC labels
                  mccasclabel::McParticleId);
using McCascLabel = McCascLabels::iterator;

} // namespace o2::aod

#endif // O2_ANALYSIS_STRANGENESSTABLES_H_
