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

/// \file ConversionPointProfile.h
/// \brief radius-profiled conversion point of a photon conversion (alternative to the helix touching point)
/// \author Stefanie Mrozinski, stefanie.mrozinski@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_CONVERSIONPOINTPROFILE_H_
#define PWGEM_PHOTONMESON_UTILS_CONVERSIONPOINTPROFILE_H_

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <DetectorsBase/Propagator.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/TrackParametrization.h>

#include <Math/MatrixFunctions.h>
#include <Math/MatrixRepresentationsStatic.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>

#include <algorithm>
#include <array>
#include <cmath>

namespace o2::pwgem::photonmeson
{

constexpr int NProfPar = 4; // (rphi, z, phi, tgl)
constexpr int NProfR = 75;
constexpr float ProfRMin = 10.f; // trial radii 10 .. 84 cm

/// \brief chi2 of a leg pair vs the trial radius, its minimum and class
struct RadiusProfile {
  std::array<float, NProfR> chi2{};
  std::array<bool, NProfR> ok{};
  int status = 4;
  int iMin = -1;
  float rMin = -1.f;
  float chi2Min = -1.f;
  float halfWidth = -1.f;
  std::array<float, 3> xyz{-999.f, -999.f, -999.f};
};

/// \brief parameters and covariance of one leg at every trial radius
struct LegScan {
  std::array<bool, NProfR> ok{};
  std::array<std::array<float, 4>, NProfR> m{};   // phi of the point, z, phi of the direction, tgl
  std::array<std::array<float, 10>, NProfR> c{};  // lower triangle of the (rphi, z, phi, tgl) covariance
  std::array<std::array<float, 3>, NProfR> pos{}; // global point
};

//_______________________________________________________________________
/// \brief propagate one leg inward from its IU to a radius and keep its parameters there
/// \param trackIU leg at its innermost update
/// \param bz magnetic field
/// \param matCorr material correction of the propagation
/// \param out output: parameters and covariance per radius, ok = false outside the reach of the leg
inline void cpScanLeg(o2::track::TrackParCov const& trackIU, const float bz, const o2::base::Propagator::MatCorrType matCorr, LegScan& out)
{
  auto trk = trackIU;
  const auto iu = trk.getXYZGlo();
  const float rIU = std::hypot(iu.X(), iu.Y());
  for (int k = NProfR - 1; k >= 0; --k) {
    const float r = ProfRMin + k;
    out.ok[k] = false;
    if (r >= rIU - 0.1f) { // o2-linter: disable=magic-number (margin)
      continue;
    }
    float xl = 0.f;
    if (!trk.getXatLabR(r, xl, bz, o2::track::DirInward)) {
      break;
    }
    if (!o2::base::Propagator::Instance()->propagateToX(trk, xl, bz, o2::base::PropagatorImpl<float>::MAX_SIN_PHI, o2::base::PropagatorImpl<float>::MAX_STEP, matCorr)) {
      break;
    }
    const auto g = trk.getXYZGlo();
    const float snp = trk.getSnp();
    const float f = 1.f / std::sqrt(std::max(1.f - snp * snp, 1e-4f)); // o2-linter: disable=magic-number (floor)
    out.m[k] = {std::atan2(g.Y(), g.X()), trk.getZ(), trk.getAlpha() + std::asin(std::clamp(snp, -1.f, 1.f)), trk.getTgl()};
    out.c[k] = {trk.getSigmaY2(), trk.getSigmaZY(), trk.getSigmaZ2(), trk.getSigmaSnpY() * f, trk.getSigmaSnpZ() * f, trk.getSigmaSnp2() * f * f, trk.getSigmaTglY(), trk.getSigmaTglZ(), trk.getSigmaTglSnp() * f, trk.getSigmaTgl2()};
    out.pos[k] = {g.X(), g.Y(), g.Z()};
    out.ok[k] = true;
  }
}

inline void cpProfileFinish(LegScan const& p, LegScan const& n, RadiusProfile& out);

//_______________________________________________________________________
/// \brief chi2 of a common point and a common direction of two legs at every trial radius, then minimum and class
/// \param p scan of the positive leg
/// \param n scan of the negative leg
/// \param zSig2Pos extra z variance of the positive leg (cm^2, legs without ITS)
/// \param zSig2Ele extra z variance of the negative leg (cm^2, legs without ITS)
/// \param out output: profile, minimum, half width and class
inline void cpRadiusProfile(LegScan const& p, LegScan const& n, const float zSig2Pos, const float zSig2Ele, RadiusProfile& out)
{
  out = RadiusProfile{};
  for (int k = 0; k < NProfR; ++k) {
    out.ok[k] = p.ok[k] && n.ok[k];
    if (!out.ok[k]) {
      continue;
    }
    const float r = ProfRMin + k;
    ROOT::Math::SVector<double, NProfPar> d;
    d[0] = r * RecoDecay::constrainAngle(p.m[k][0] - n.m[k][0], -o2::constants::math::PI);
    d[1] = p.m[k][1] - n.m[k][1];
    d[2] = RecoDecay::constrainAngle(p.m[k][2] - n.m[k][2], -o2::constants::math::PI);
    d[3] = p.m[k][3] - n.m[k][3];
    ROOT::Math::SMatrix<double, NProfPar, NProfPar, ROOT::Math::MatRepSym<double, NProfPar>> cov;
    int idx = 0;
    for (int i = 0; i < NProfPar; ++i) {
      for (int j = 0; j <= i; ++j) {
        cov(i, j) = p.c[k][idx] + n.c[k][idx];
        ++idx;
      }
    }
    cov(1, 1) += zSig2Pos + zSig2Ele;
    if (!cov.Invert()) {
      out.ok[k] = false;
      continue;
    }
    out.chi2[k] = static_cast<float>(ROOT::Math::Similarity(d, cov));
    if (out.iMin < 0 || out.chi2[k] < out.chi2Min) {
      out.iMin = k;
      out.chi2Min = out.chi2[k];
    }
  }
  cpProfileFinish(p, n, out);
}
//_______________________________________________________________________
/// \brief minimum, half width of the dchi2 < 1 region, other minima and class of a filled profile
/// \param p scan of the positive leg (point at the minimum)
/// \param n scan of the negative leg (point at the minimum)
/// \param out profile, completed in place
inline void cpProfileFinish(LegScan const& p, LegScan const& n, RadiusProfile& out)
{
  if (out.iMin < 0) {
    return;
  }
  out.rMin = ProfRMin + out.iMin;
  out.xyz = {0.5f * (p.pos[out.iMin][0] + n.pos[out.iMin][0]), 0.5f * (p.pos[out.iMin][1] + n.pos[out.iMin][1]), 0.5f * (p.pos[out.iMin][2] + n.pos[out.iMin][2])};
  // valid range, dchi2 < 1 region around the minimum, other minima
  int kLo = out.iMin;
  int kHi = kLo;
  while (kLo - 1 >= 0 && out.ok[kLo - 1] && out.chi2[kLo - 1] - out.chi2Min < 1.f) {
    --kLo;
  }
  while (kHi + 1 < NProfR && out.ok[kHi + 1] && out.chi2[kHi + 1] - out.chi2Min < 1.f) {
    ++kHi;
  }
  out.halfWidth = 0.5f * (kHi - kLo) + 0.5f;
  const bool edge = (out.iMin == 0 || !out.ok[out.iMin - 1]) || (out.iMin == NProfR - 1 || !out.ok[out.iMin + 1]);
  bool other = false;
  for (int k = 1; k < NProfR - 1 && !other; ++k) {
    if (k >= kLo && k <= kHi) {
      continue;
    }
    if (out.ok[k] && out.ok[k - 1] && out.ok[k + 1] && out.chi2[k] <= out.chi2[k - 1] && out.chi2[k] <= out.chi2[k + 1] && out.chi2[k] - out.chi2Min < 1.f) {
      other = true;
    }
  }
  out.status = edge ? 3 : (other ? 2 : (out.halfWidth > 5.f ? 1 : 0)); // o2-linter: disable=magic-number (broad: half width above 5 cm)
}

} // namespace o2::pwgem::photonmeson

#endif // PWGEM_PHOTONMESON_UTILS_CONVERSIONPOINTPROFILE_H_
