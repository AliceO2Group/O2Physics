// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoMath.h
/// \brief shared kinematic math (kstar, Q3, mt, kt, spherical-harmonics kinematics, ...)
///        used by the pair and triplet histogram managers
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_FEMTOMATH_H_
#define PWGCF_FEMTO_CORE_FEMTOMATH_H_

#include "PWGCF/Femto/Core/modes.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/Logger.h>

#include <Math/Vector3D.h> // IWYU pragma: keep (do not replace with Math/Vector3Dfwd.h)
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

namespace o2::analysis::femto::femtomath
{

// SH reference frame selectors for computeShKinematics(), mirrors ConfPairBinning.shFrame
constexpr int ShFrameLcmsNonIdentical = 0;
constexpr int ShFrameLcmsIdentical = 1;
constexpr int ShFramePrf = 2;

inline float getKt(ROOT::Math::PtEtaPhiMVector const& part1, ROOT::Math::PtEtaPhiMVector const& part2)
{
  auto sum = (part1 + part2);
  double kt = 0.5 * sum.Pt();
  return static_cast<float>(kt);
}

inline float getMt(ROOT::Math::PtEtaPhiMVector const& part1, ROOT::Math::PtEtaPhiMVector const& part2, modes::TransverseMassType mtType = modes::TransverseMassType::kAveragePdgMass)
{
  auto sum = part1 + part2;
  double mt = 0;
  double mass = 0;
  switch (mtType) {
    case modes::TransverseMassType::kAveragePdgMass:
      mass = 0.5 * (part1.M() + part2.M());
      mt = std::hypot(0.5 * sum.Pt(), mass);
      break;
    case modes::TransverseMassType::kReducedPdgMass:
      mass = 2. * (part1.M() * part2.M()) / (part1.M() + part2.M());
      mt = std::hypot(0.5 * sum.Pt(), mass);
      break;
    case modes::TransverseMassType::kMt4Vector:
      mt = 0.5 * sum.Mt();
      break;
    default:
      LOG(fatal) << "Invalid transverse mass type, breaking...";
  }
  return static_cast<float>(mt);
}

inline float getMinv(ROOT::Math::PtEtaPhiMVector const& part1, ROOT::Math::PtEtaPhiMVector const& part2)
{
  return static_cast<float>((part1 + part2).M());
}

inline float getKstar(ROOT::Math::PtEtaPhiMVector const& part1, ROOT::Math::PtEtaPhiMVector const& part2)
{
  // Use Cartesian 4-vectors: addition/M2() become pure arithmetic
  const ROOT::Math::PxPyPzEVector p1(part1);
  const ROOT::Math::PxPyPzEVector p2(part2);

  // Mandelstam s = (p1 + p2)^2
  const double s = (p1 + p2).M2();
  const double m1sq = p1.M2();
  const double m2sq = p2.M2();

  // Källen function λ(s, m1^2, m2^2) = (s - m1^2 - m2^2)² - 4*m1^2*m2^2
  const double kallen = (s - m1sq - m2sq) * (s - m1sq - m2sq) - 4.0 * m1sq * m2sq;

  // k* = 0.5 * sqrt(λ/s)
  return static_cast<float>(0.5 * std::sqrt(std::max(0.0, kallen) / s));
}

// Alpha: azimuthal opening angle between the two particles' transverse momenta,
// i.e. phi2 - phi1 wrapped to [-pi, pi). Purely a lab-frame angular variable (no
// boosting involved) -- small |alpha| means the pair sits close together in phi,
// which is where detector effects like track splitting/merging and two-track
// resolution losses tend to show up.
inline float getAlpha(ROOT::Math::PtEtaPhiMVector const& part1, ROOT::Math::PtEtaPhiMVector const& part2)
{
  return static_cast<float>(RecoDecay::constrainAngle(part2.phi() - part1.phi(), -o2::constants::math::PI));
}

// BetaLab: lab-frame angle between the pair's relative 3-momentum (p1 - p2) and its
// total 3-momentum (p1 + p2, i.e. the kT direction). Note the (p1 - p2) here is the
// raw 3-vector difference, not the Lorentz-invariant k* -- this is a lab-frame proxy
// for how the pair's separation axis is oriented relative to its direction of motion,
// without needing to boost to the pair rest frame.
inline float getBetaLab(ROOT::Math::PtEtaPhiMVector const& part1, ROOT::Math::PtEtaPhiMVector const& part2)
{
  ROOT::Math::XYZVector p1(part1.Px(), part1.Py(), part1.Pz());
  ROOT::Math::XYZVector p2(part2.Px(), part2.Py(), part2.Pz());

  ROOT::Math::XYZVector kstar = p1 - p2;
  ROOT::Math::XYZVector kt = p1 + p2;

  return static_cast<float>(std::acos(kstar.Dot(kt) / std::sqrt(kstar.Mag2() * kt.Mag2())));
}

inline std::tuple<float, float, float> computeBertschPrattLCMS(ROOT::Math::PtEtaPhiMVector const& part1, ROOT::Math::PtEtaPhiMVector const& part2)
{
  const ROOT::Math::PxPyPzEVector p1(part1);
  const ROOT::Math::PxPyPzEVector p2(part2);
  const ROOT::Math::PxPyPzEVector pSum = p1 + p2;
  const ROOT::Math::PxPyPzEVector pDiff = p1 - p2;

  const double tPx = pSum.Px();
  const double tPy = pSum.Py();
  const double tPz = pSum.Pz();
  const double tE = pSum.E();

  const double tPt = std::sqrt(tPx * tPx + tPy * tPy);
  const double tMt = std::sqrt(tE * tE - tPz * tPz);

  static constexpr double MinTransverseMomentum = 1e-9;
  if (tPt < MinTransverseMomentum || tMt < MinTransverseMomentum) {
    return {0.0f, 0.0f, 0.0f};
  }

  const double betaL = tPz / tE;
  const double gammaL = tE / tMt;

  // out/side/long are linear projections of the momentum, so projecting the
  // difference p1-p2 directly is equivalent to (and half the divisions of)
  // computing each particle's out/side/long component and subtracting.
  const double qOut = (pDiff.Px() * tPx + pDiff.Py() * tPy) / tPt;
  const double qSide = (pDiff.Py() * tPx - pDiff.Px() * tPy) / tPt;
  const double qLong = gammaL * (pDiff.Pz() - betaL * pDiff.E());

  return {static_cast<float>(qOut), static_cast<float>(qSide), static_cast<float>(qLong)};
}

// Return the bin index for value given ascending bin edges, or -1 if out of range.
// edges = {e0, e1, ..., eN} defines N bins [e0,e1), [e1,e2), ..., [e_{N-1},eN).
inline int findShBin(double value, std::vector<double> const& edges)
{
  static constexpr std::size_t MinEdgesForOneBin = 2;
  if (edges.size() < MinEdgesForOneBin || value < edges.front() || value >= edges.back()) {
    return -1;
  }
  for (std::size_t i = 0; i < edges.size() - 1; ++i) {
    if (value >= edges[i] && value < edges[i + 1]) {
      return static_cast<int>(i);
    }
  }
  return -1;
}

// Kinematics feeding the spherical-harmonics decomposition, ported 1:1 from
// FemtoUniverseMath::newpairfunc so results are numerically comparable.
// Returns {kv, out, side, long}, selected by shFrame:
//   ShFrameLcmsNonIdentical: {kstar,  fDKOut,  fDKSide,  fDKLong}   (LCMS components of particle 1)
//   ShFrameLcmsIdentical:    {qinv,   outLCMS, sideLCMS, longLCMS}  (LCMS q-differences)
//   ShFramePrf:              {|qPRF|, outPRF,  sidePRF,  longPRF}   (matches FemtoUniverse isIdenPRF=true)
inline std::tuple<float, float, float, float> computeShKinematics(ROOT::Math::PtEtaPhiMVector const& part1,
                                                                  ROOT::Math::PtEtaPhiMVector const& part2,
                                                                  int shFrame)
{
  const ROOT::Math::PxPyPzEVector p1(part1);
  const ROOT::Math::PxPyPzEVector p2(part2);
  const ROOT::Math::PxPyPzEVector sum = p1 + p2;

  const double tPx = sum.Px();
  const double tPy = sum.Py();
  const double tPz = sum.Pz();
  const double tE = sum.E();
  const double tPtSq = tPx * tPx + tPy * tPy;
  const double tMtSq = tE * tE - tPz * tPz;

  static constexpr double MinTransverseMomentum = 1e-9;
  if (tPtSq < MinTransverseMomentum || tMtSq < MinTransverseMomentum) {
    return {0.f, 0.f, 0.f, 0.f};
  }
  const double tPt = std::sqrt(tPtSq);
  const double tMt = std::sqrt(tMtSq);
  const double tM = std::sqrt(std::max(0.0, tMtSq - tPtSq));

  const double beta = tPz / tE;
  const double gamma = tE / tMt;

  // LCMS components of particle 1 (fDKOut/Side/Long in newpairfunc)
  const double e1 = p1.E();
  const double fDKOut = (p1.Px() * tPx + p1.Py() * tPy) / tPt;
  const double fDKSide = (-p1.Px() * tPy + p1.Py() * tPx) / tPt;
  const double fDKLong = gamma * (p1.Pz() - beta * e1);
  const double fDE = gamma * (e1 - beta * p1.Pz());

  // LCMS components of particle 2
  const double e2 = p2.E();
  const double px2 = (p2.Px() * tPx + p2.Py() * tPy) / tPt;
  const double py2 = (p2.Py() * tPx - p2.Px() * tPy) / tPt;
  const double pz2 = gamma * (p2.Pz() - beta * e2);

  const double outLCMS = fDKOut - px2;
  const double sideLCMS = fDKSide - py2;
  const double longLCMS = fDKLong - pz2;

  if (shFrame == ShFrameLcmsIdentical) {
    // LCMS identical: Ylm <- LCMS q-differences, axis <- qinv
    const double qinv = std::sqrt(outLCMS * outLCMS + sideLCMS * sideLCMS + longLCMS * longLCMS);
    return {static_cast<float>(qinv),
            static_cast<float>(outLCMS), static_cast<float>(sideLCMS), static_cast<float>(longLCMS)};
  }
  if (shFrame == ShFramePrf) {
    // PRF: Ylm <- PRF q-differences, axis <- |q_PRF| (matches FemtoUniverse isIdenPRF=true)
    const double pE2LCMS = gamma * (e2 - beta * p2.Pz());
    const double betaOut = tPt / tMt;
    const double gammaOut = tMt / tM;
    const double outPRF = gammaOut * (outLCMS - betaOut * (fDE - pE2LCMS));
    const double sidePRF = sideLCMS;
    const double longPRF = longLCMS;
    const double qPRF = std::sqrt(outPRF * outPRF + sidePRF * sidePRF + longPRF * longPRF);
    return {static_cast<float>(qPRF),
            static_cast<float>(outPRF), static_cast<float>(sidePRF), static_cast<float>(longPRF)};
  }
  // LCMS non-identical (default): Ylm <- LCMS components of particle 1, axis <- k*
  const double betaOut = tPt / tMt;
  const double gammaOut = tMt / tM;
  const double fKOut = gammaOut * (fDKOut - betaOut * fDE);
  const double kstar = std::sqrt(fKOut * fKOut + fDKSide * fDKSide + fDKLong * fDKLong);
  return {static_cast<float>(kstar),
          static_cast<float>(fDKOut), static_cast<float>(fDKSide), static_cast<float>(fDKLong)};
}

template <typename T1>
inline float getPairPhiFromEventPlaneAngle(const ROOT::Math::PtEtaPhiMVector& part1, const ROOT::Math::PtEtaPhiMVector& part2, const T1& col)
{
  const ROOT::Math::PtEtaPhiMVector trackSum = part1 + part2;
  float phi = RecoDecay::constrainAngle(trackSum.Phi() - col.eventPlaneAngle());
  return phi;
}

// triplet-only helpers

inline ROOT::Math::PxPyPzEVector getqij(ROOT::Math::PxPyPzEVector const& vi, ROOT::Math::PxPyPzEVector const& vj)
{
  auto trackSum = vi + vj;
  auto trackDifference = vi - vj;
  const double s = trackSum.M2();
  const double scaling = (s != 0.0) ? (vi.M2() - vj.M2()) / s : 0.0;
  return trackDifference - scaling * trackSum;
}

inline float getQ3(ROOT::Math::PtEtaPhiMVector const& part1, ROOT::Math::PtEtaPhiMVector const& part2, ROOT::Math::PtEtaPhiMVector const& part3)
{
  // upfront conversion to PxPyPzEVector
  const ROOT::Math::PxPyPzEVector p1(part1);
  const ROOT::Math::PxPyPzEVector p2(part2);
  const ROOT::Math::PxPyPzEVector p3(part3);
  auto q12 = getqij(p1, p2);
  auto q23 = getqij(p2, p3);
  auto q31 = getqij(p3, p1);
  double q = q12.M2() + q23.M2() + q31.M2();
  return static_cast<float>(std::sqrt(std::max(0.0, -q)));
}

} // namespace o2::analysis::femto::femtomath

#endif // PWGCF_FEMTO_CORE_FEMTOMATH_H_
