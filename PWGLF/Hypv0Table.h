#ifndef O2_ANALYSIS_HYPV0TABLE_H_
#define O2_ANALYSIS_HYPV0TABLE_H_

#include "Common/Core/RecoDecay.h"
#include <cmath>

namespace o2::aod
{
  namespace hypv0data{

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

// Derived expressions
// Momenta
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //! V0 pT
                           [](float pxpos, float pypos, float pxneg, float pyneg) -> float { return RecoDecay::sqrtSumOfSquares(pxpos + pxneg, pypos + pyneg); });

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
                           [](float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ) -> float { return RecoDecay::CPA(array{pvX, pvY, pvZ}, array{X, Y, Z}, array{Px, Py, Pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(DCAV0ToPV, dcav0topv, //! DCA of V0 to PV
                           [](float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ) -> float { return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz)); });

// Armenteros-Podolanski variables
DECLARE_SOA_DYNAMIC_COLUMN(Alpha, alpha, //! Armenteros Alpha
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             float momTot = RecoDecay::P(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
                             float lQlNeg = RecoDecay::dotProd(array{pxneg, pyneg, pzneg}, array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}) / momTot;
                             float lQlPos = RecoDecay::dotProd(array{pxpos, pypos, pzpos}, array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg}) / momTot;
                             return (lQlPos - lQlNeg) / (lQlPos + lQlNeg); // alphav0
                           });

DECLARE_SOA_DYNAMIC_COLUMN(QtArm, qtarm, //! Armenteros Qt
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             float momTot = RecoDecay::P2(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
                             float dp = RecoDecay::dotProd(array{pxneg, pyneg, pzneg}, array{pxpos + pxneg, pypos + pyneg, pzpos + pzneg});
                             return std::sqrt(RecoDecay::P2(pxneg, pyneg, pzneg) - dp * dp / momTot); // qtarm
                           });

// Psi pair angle: angle between the plane defined by the electron and positron momenta and the xy plane
DECLARE_SOA_DYNAMIC_COLUMN(PsiPair, psipair, //! psi pair angle
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) {
                             auto clipToPM1 = [](float x) { return x < -1.f ? -1.f : (x > 1.f ? 1.f : x); };
                             float ptot2 = RecoDecay::P2(pxpos, pypos, pzpos) * RecoDecay::P2(pxneg, pyneg, pzneg);
                             float argcos = RecoDecay::dotProd(array{pxpos, pypos, pzpos}, array{pxneg, pyneg, pzneg}) / std::sqrt(ptot2);
                             float thetaPos = std::atan2(RecoDecay::sqrtSumOfSquares(pxpos, pypos), pzpos);
                             float thetaNeg = std::atan2(RecoDecay::sqrtSumOfSquares(pxneg, pyneg), pzneg);
                             float argsin = (thetaNeg - thetaPos) / std::acos(clipToPM1(argcos));
                             return std::asin(clipToPM1(argsin));
                           });

// Calculated on the fly with mass assumption + dynamic tables
DECLARE_SOA_DYNAMIC_COLUMN(MHypertriton, mHypertriton, //! mass under Hypertriton hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::M(array{array{pxpos, pypos, pzpos}, array{pxneg, pyneg, pzneg}}, array{2.80839, RecoDecay::getMassPDG(kPiPlus)}); });
DECLARE_SOA_DYNAMIC_COLUMN(MAntiHypertriton, mAntiHypertriton, //! mass under antihypertriton hypothesis
                           [](float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg) -> float { return RecoDecay::M(array{array{pxpos, pypos, pzpos}, array{pxneg, pyneg, pzneg}}, array{RecoDecay::getMassPDG(kPiPlus), 2.80839}); });

DECLARE_SOA_DYNAMIC_COLUMN(YHypertriton, yHypertriton, //! V0 y with hypertriton or antihypertriton hypothesis
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::Y(array{Px, Py, Pz}, 2.991); });
DECLARE_SOA_DYNAMIC_COLUMN(Eta, eta, //! V0 eta
                           [](float Px, float Py, float Pz) -> float { return RecoDecay::Eta(array{Px, Py, Pz}); });
DECLARE_SOA_DYNAMIC_COLUMN(Phi, phi, //! V0 phi
                           [](float Px, float Py) -> float { return RecoDecay::Phi(Px, Py); });


DECLARE_SOA_DYNAMIC_COLUMN(NegativePt, negativept, //! negative daughter pT
                           [](float pxneg, float pyneg) -> float { return RecoDecay::sqrtSumOfSquares(pxneg, pyneg); });
DECLARE_SOA_DYNAMIC_COLUMN(PositivePt, positivept, //! positive daughter pT
                           [](float pxpos, float pypos) -> float { return RecoDecay::sqrtSumOfSquares(pxpos, pypos); });
DECLARE_SOA_DYNAMIC_COLUMN(NegativeEta, negativeeta, //! negative daughter eta
                           [](float PxNeg, float PyNeg, float PzNeg) -> float { return RecoDecay::Eta(array{PxNeg, PyNeg, PzNeg}); });
DECLARE_SOA_DYNAMIC_COLUMN(NegativePhi, negativephi, //! negative daughter phi
                           [](float PxNeg, float PyNeg) -> float { return RecoDecay::Phi(PxNeg, PyNeg); });
DECLARE_SOA_DYNAMIC_COLUMN(PositiveEta, positiveeta, //! positive daughter eta
                           [](float PxPos, float PyPos, float PzPos) -> float { return RecoDecay::Eta(array{PxPos, PyPos, PzPos}); });
DECLARE_SOA_DYNAMIC_COLUMN(PositivePhi, positivephi, //! positive daughter phi
                           [](float PxPos, float PyPos) -> float { return RecoDecay::Phi(PxPos, PyPos); });
DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //! V0 px
                              float, 1.f * aod::hypv0data::pxpos + 1.f * aod::hypv0data::pxneg);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //! V0 py
                              float, 1.f * aod::hypv0data::pypos + 1.f * aod::hypv0data::pyneg);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //! V0 pz
                              float, 1.f * aod::hypv0data::pzpos + 1.f * aod::hypv0data::pzneg);
}

DECLARE_SOA_TABLE(HypStoredV0Datas,"AOD", "HypV0DATA", //!
//DECLARE_SOA_TABLE_FULL(HypStoredV0Datas, "HypV0Datas", "AOD", "HypV0DATA", //!
                       o2::soa::Index<>, hypv0data::PosTrackId, hypv0data::NegTrackId, hypv0data::CollisionId, hypv0data::V0Id,
                       hypv0data::PosX, hypv0data::NegX,
                       hypv0data::X, hypv0data::Y, hypv0data::Z,
                       hypv0data::PxPos, hypv0data::PyPos, hypv0data::PzPos,
                       hypv0data::PxNeg, hypv0data::PyNeg, hypv0data::PzNeg,
                       hypv0data::DCAV0Daughters, hypv0data::DCAPosToPV, hypv0data::DCANegToPV,

                       // Dynamic columns
                       hypv0data::Pt<hypv0data::PxPos, hypv0data::PyPos, hypv0data::PxNeg, hypv0data::PyNeg>,
                       hypv0data::V0Radius<hypv0data::X, hypv0data::Y>,
                       hypv0data::DistOverTotMom<hypv0data::X, hypv0data::Y, hypv0data::Z, hypv0data::Px, hypv0data::Py, hypv0data::Pz>,
                       hypv0data::V0CosPA<hypv0data::X, hypv0data::Y, hypv0data::Z, hypv0data::Px, hypv0data::Py, hypv0data::Pz>,
                       hypv0data::DCAV0ToPV<hypv0data::X, hypv0data::Y, hypv0data::Z, hypv0data::Px, hypv0data::Py, hypv0data::Pz>,
                       hypv0data::Alpha<hypv0data::PxPos, hypv0data::PyPos, hypv0data::PzPos, hypv0data::PxNeg, hypv0data::PyNeg, hypv0data::PzNeg>,
                       hypv0data::QtArm<hypv0data::PxPos, hypv0data::PyPos, hypv0data::PzPos, hypv0data::PxNeg, hypv0data::PyNeg, hypv0data::PzNeg>,
                       hypv0data::PsiPair<hypv0data::PxPos, hypv0data::PyPos, hypv0data::PzPos, hypv0data::PxNeg, hypv0data::PyNeg, hypv0data::PzNeg>,

                       // Invariant masses
                       hypv0data::MHypertriton<hypv0data::PxPos, hypv0data::PyPos, hypv0data::PzPos, hypv0data::PxNeg, hypv0data::PyNeg, hypv0data::PzNeg>,
                       hypv0data::MAntiHypertriton<hypv0data::PxPos, hypv0data::PyPos, hypv0data::PzPos, hypv0data::PxNeg, hypv0data::PyNeg, hypv0data::PzNeg>,

                       // Longitudinal
                       hypv0data::YHypertriton<hypv0data::Px, hypv0data::Py, hypv0data::Pz>,
                       hypv0data::Eta<hypv0data::Px, hypv0data::Py, hypv0data::Pz>,
                       hypv0data::Phi<hypv0data::Px, hypv0data::Py>,
                       hypv0data::NegativePt<hypv0data::PxNeg, hypv0data::PyNeg>,
                       hypv0data::PositivePt<hypv0data::PxPos, hypv0data::PyPos>,
                       hypv0data::NegativeEta<hypv0data::PxNeg, hypv0data::PyNeg, hypv0data::PzNeg>,
                       hypv0data::NegativePhi<hypv0data::PxNeg, hypv0data::PyNeg>,
                       hypv0data::PositiveEta<hypv0data::PxPos, hypv0data::PyPos, hypv0data::PzPos>,
                       hypv0data::PositivePhi<hypv0data::PxPos, hypv0data::PyPos>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HypV0Datas, HypStoredV0Datas, "HYPV0DATAEXT", //!
                                hypv0data::Px, hypv0data::Py, hypv0data::Pz);
using HypV0Data = HypV0Datas::iterator;
namespace hypv0data
{
DECLARE_SOA_INDEX_COLUMN(HypV0Data, Hypv0Data); // Inde to HypV0Data entry
}

}

#endif
