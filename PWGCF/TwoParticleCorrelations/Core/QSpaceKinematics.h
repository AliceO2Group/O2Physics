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

/// \file QSpaceKinematics.h
/// \brief Relative momentum of a particle pair in the pair rest frame (Q_inv, Q_long, Q_out, Q_side)
/// \author Lakshitha Fernando <lakshithafernando@wayne.edu>
/// \since Nov/08/2026
///
/// Definitions follow Pratt and Cheng, PRC 68 (2003) 014907 and Pruneau et al., PRC 107 (2023) 014902.
/// The projection steps follow newpairfunc() in PWGCF/FemtoUniverse/Core/FemtoUniverseMath.h.

#ifndef PWGCF_TWOPARTICLECORRELATIONS_CORE_QSPACEKINEMATICS_H_
#define PWGCF_TWOPARTICLECORRELATIONS_CORE_QSPACEKINEMATICS_H_

#include <algorithm>
#include <cmath>

namespace o2::analysis::qspace
{

/// Pair variables in GeV/c (mass in GeV/c^2). The components are signed: exchanging the two particles flips them.
struct PairQ {
  double qInv = 0.;  ///< from Lorentz invariants, never from the components
  double qLong = 0.; ///< along the beam axis
  double qOut = 0.;  ///< along the pair transverse momentum
  double qSide = 0.; ///< transverse, perpendicular to out; (out, side, long) is right-handed
  double mass = 0.;  ///< pair invariant mass
};

/// Computes the pair variables of particle 1 relative to particle 2.
/// Sequence: Bertsch-Pratt axes, boost along the beam to the frame where the pair has P_z = 0 (LCMS),
/// then boost along out to the pair rest frame, so that qLong^2 + qOut^2 + qSide^2 = qInv^2.
/// For equal masses qInv = 2 k*, with k* the momentum of either particle in the pair rest frame.
inline PairQ computePairQ(double pt1, double eta1, double phi1, double mass1,
                          double pt2, double eta2, double phi2, double mass2)
{
  const double px1 = pt1 * std::cos(phi1);
  const double py1 = pt1 * std::sin(phi1);
  const double pz1 = pt1 * std::sinh(eta1);
  const double e1 = std::sqrt(px1 * px1 + py1 * py1 + pz1 * pz1 + mass1 * mass1);
  const double px2 = pt2 * std::cos(phi2);
  const double py2 = pt2 * std::sin(phi2);
  const double pz2 = pt2 * std::sinh(eta2);
  const double e2 = std::sqrt(px2 * px2 + py2 * py2 + pz2 * pz2 + mass2 * mass2);

  const double sumPx = px1 + px2;
  const double sumPy = py1 + py2;
  const double sumPz = pz1 + pz2;
  const double sumE = e1 + e2;
  const double diffPx = px1 - px2;
  const double diffPy = py1 - py2;
  const double diffPz = pz1 - pz2;
  const double diffE = e1 - e2;

  const double pairPt2 = sumPx * sumPx + sumPy * sumPy;
  const double pairPt = std::sqrt(pairPt2);
  const double pairMt2 = sumE * sumE - sumPz * sumPz;
  const double pairMt = std::sqrt(pairMt2);
  const double pairS = pairMt2 - pairPt2;

  PairQ result;
  result.mass = std::sqrt(std::max(pairS, 0.));

  // invariant definition: Q_inv^2 = -(p1-p2)^2 + (m1^2-m2^2)^2/s
  const double massTerm = mass1 * mass1 - mass2 * mass2;
  const double qInv2 = diffPx * diffPx + diffPy * diffPy + diffPz * diffPz - diffE * diffE + massTerm * massTerm / pairS;
  result.qInv = std::sqrt(std::max(qInv2, 0.));

  // transverse axes; a pair with exactly zero transverse momentum has no out direction, then x is used
  const double outX = pairPt > 0. ? sumPx / pairPt : 1.;
  const double outY = pairPt > 0. ? sumPy / pairPt : 0.;
  const double diffOut = diffPx * outX + diffPy * outY;
  result.qSide = -diffPx * outY + diffPy * outX;

  // boost along the beam to P_z = 0: gamma = E/M_T, gamma*beta = P_z/M_T
  result.qLong = (sumE * diffPz - sumPz * diffE) / pairMt;
  const double diffELcms = (sumE * diffE - sumPz * diffPz) / pairMt;

  // boost along out to the rest frame: gamma = M_T/M, gamma*beta = P_T/M
  result.qOut = (pairMt * diffOut - pairPt * diffELcms) / result.mass;
  return result;
}

} // namespace o2::analysis::qspace

#endif // PWGCF_TWOPARTICLECORRELATIONS_CORE_QSPACEKINEMATICS_H_
