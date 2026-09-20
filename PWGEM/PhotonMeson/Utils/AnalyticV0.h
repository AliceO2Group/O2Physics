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

/// \file AnalyticV0.h
/// \brief Build a photon V0 candidate from two legs without a vertex fitter.
/// \author Stefanie Mrozinski, stefanie.mrozinski@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_ANALYTICV0_H_
#define PWGEM_PHOTONMESON_UTILS_ANALYTICV0_H_

#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h" // CalculateDCAFast, getPropMomentumFromTrackHelix (soa version)

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <MathUtils/Utils.h>
#include <ReconstructionDataFormats/HelixHelper.h>
#include <ReconstructionDataFormats/Track.h>

#include <algorithm>
#include <array>
#include <cmath>

namespace o2::pwgem::photonmeson
{

struct LegParam {
  std::array<float, 3> xyz{}; // conversion point, vertex-relative (cm)
  std::array<float, 3> mom{}; // momentum at the conversion point (GeV/c)
  int charge{0};              // +1 e+, -1 e-
};

inline o2::track::TrackParCov trackFromLeg(LegParam const& leg)
{
  std::array<float, 21> cov{};        // o2-linter: disable=magic-number (lab-frame covariance has 21 elements)
  cov[0] = cov[2] = cov[5] = 0.01f;   // o2-linter: disable=magic-number (placeholder sigma_xyz = 0.1 cm)
  cov[9] = cov[14] = cov[20] = 1e-4f; // o2-linter: disable=magic-number (placeholder sigma_p = 0.01 GeV/c)
  return o2::track::TrackParCov(leg.xyz, leg.mom, cov, leg.charge, true);
}

inline std::array<float, 3> getPropMomentumFromTrackHelix(const float s, const o2::track::TrackParCov& track,
                                                          const o2::track::TrackAuxPar& trHelix, const float bz, float addPhi = 0.f)
{
  const float pt = track.getPt();
  const float phi0 = track.getPhi();                                                    // global azimuth
  const float dphi = -static_cast<float>(track.getSign()) * 0.3f * bz * s / 100.f / pt; // o2-linter: disable=magic-number (0.3 GeV/(T m), s in cm)
  const float dotProd = std::cos(phi0 + addPhi) * pt * trHelix.xC + std::sin(phi0 + addPhi) * pt * trHelix.yC;
  if (dotProd < 0.f) {
    addPhi -= o2::constants::math::PI;
  }
  const float phi = RecoDecay::constrainAngle<float>(phi0 + dphi + addPhi);
  return {std::cos(phi) * pt, std::sin(phi) * pt, track.getTgl() * pt};
}

struct AnalyticV0 {
  bool ok{false};
  std::array<float, 3> vtx{};                   // conversion point (cm), vertex-relative
  std::array<float, 3> mom{};                   // photon momentum at the vertex (GeV/c)
  float pca{999.f};                             // distance of the two legs at the vertex (cm)
  float cospa{-2.f};                            // pointing of mom from the origin (= primary vertex)
  float mee{-1.f};                              // m(e+e-) at the vertex (GeV/c^2)
  float dcaxy{999.f};                           // line through vtx along mom: distance to the origin in xy (cm)
  float dcaz{999.f};                            // line through vtx along mom in z at the xy-closest point (cm)
  float score{999.f};                           // the builder's getScoreV0
  std::array<std::array<float, 3>, 2> legMom{}; // e+ and e- momentum at the vertex (GeV/c)
  void armenteros(float& alpha, float& qt) const
  {
    const float pMag = std::sqrt(mom[0] * mom[0] + mom[1] * mom[1] + mom[2] * mom[2]);
    if (pMag < 1e-9f) { // o2-linter: disable=magic-number (degenerate)
      alpha = 2.f;
      qt = 999.f;
      return;
    }
    const auto& pp = legMom[0];
    const auto& pn = legMom[1];
    const float plP = (pp[0] * mom[0] + pp[1] * mom[1] + pp[2] * mom[2]) / pMag;
    const float plN = (pn[0] * mom[0] + pn[1] * mom[1] + pn[2] * mom[2]) / pMag;
    const float pP2 = pp[0] * pp[0] + pp[1] * pp[1] + pp[2] * pp[2];
    alpha = (plP + plN > 1e-9f) ? (plP - plN) / (plP + plN) : 2.f; // o2-linter: disable=magic-number (degenerate)
    qt = std::sqrt(std::max(0.f, pP2 - plP * plP));
  }
  [[nodiscard]] float rxy() const { return std::hypot(vtx[0], vtx[1]); }
  [[nodiscard]] float pt() const { return std::hypot(mom[0], mom[1]); }
  [[nodiscard]] float eta() const
  {
    const float p = std::sqrt(mom[0] * mom[0] + mom[1] * mom[1] + mom[2] * mom[2]);
    return (p > mom[2]) ? 0.5f * std::log((p + mom[2]) / (p - mom[2])) : 0.f;
  }
  [[nodiscard]] float phi() const { return RecoDecay::constrainAngle<float>(std::atan2(mom[1], mom[0])); }
};

inline float analyticScoreV0(float cospa, float pca, float w)
{
  return w * 60.f * std::acos(cospa) + (1.f - w) * pca / 3.f; // o2-linter: disable=magic-number (getScoreV0 of the builder)
}

inline AnalyticV0 buildAnalyticV0(LegParam const& pos, LegParam const& neg, float bzkG, float scoreWeight)
{
  constexpr float kMe = 0.000510999f; // electron mass, GeV/c^2
  constexpr float kTiny = 1e-6f;
  AnalyticV0 v;
  if (std::fabs(bzkG) < 0.1f) { // o2-linter: disable=magic-number (no field)
    return v;
  }
  const auto tP = trackFromLeg(pos);
  const auto tN = trackFromLeg(neg);
  const o2::track::TrackAuxPar h1(tP, bzkG);
  const o2::track::TrackAuxPar h2(tN, bzkG);

  // 1. xy
  const float dcx = h2.xC - h1.xC, dcy = h2.yC - h1.yC;
  const float d = std::hypot(dcx, dcy);
  if (d < 1e-3f) { // o2-linter: disable=magic-number (concentric circles)
    return v;
  }
  const float ux = dcx / d, uy = dcy / d;
  const float r1 = h1.rC, r2 = h2.rC;

  int nCand = 0;
  float c1x[2], c1y[2], c2x[2], c2y[2];
  float pcaXY = 0.f;
  if (d > r1 + r2) { // external gap: tangent points along the centre line
    c1x[0] = h1.xC + r1 * ux;
    c1y[0] = h1.yC + r1 * uy;
    c2x[0] = h2.xC - r2 * ux;
    c2y[0] = h2.yC - r2 * uy;
    pcaXY = d - r1 - r2;
    nCand = 1;
  } else if (d < std::fabs(r1 - r2)) { // one circle inside the other: closest points on the far side
    const float sgn = (r1 > r2) ? 1.f : -1.f;
    c1x[0] = h1.xC + sgn * r1 * ux;
    c1y[0] = h1.yC + sgn * r1 * uy;
    c2x[0] = h2.xC + sgn * r2 * ux;
    c2y[0] = h2.yC + sgn * r2 * uy;
    pcaXY = std::fabs(r1 - r2) - d;
    nCand = 1;
  } else { // the circles cross: two points with xy distance 0
    const float a = (r1 * r1 - r2 * r2 + d * d) / (2.f * d);
    const float hh = std::sqrt(std::max(0.f, r1 * r1 - a * a));
    const float bx = h1.xC + a * ux, by = h1.yC + a * uy;
    c1x[0] = c2x[0] = bx - hh * uy;
    c1y[0] = c2y[0] = by + hh * ux;
    c1x[1] = c2x[1] = bx + hh * uy;
    c1y[1] = c2y[1] = by - hh * ux;
    pcaXY = 0.f;
    nCand = 2;
  }

  int best = 0;
  float z1 = 0.f, z2 = 0.f, dzBest = 1e9f;
  for (int ic = 0; ic < nCand; ++ic) {
    const o2::math_utils::Point3D<float> xy{0.5f * (c1x[ic] + c2x[ic]), 0.5f * (c1y[ic] + c2y[ic]), 0.f};
    const auto dcaP = CalculateDCAFast(tP, xy, bzkG);
    const auto dcaN = CalculateDCAFast(tN, xy, bzkG);
    if (std::fabs(dcaP[1] - dcaN[1]) < dzBest) {
      dzBest = std::fabs(dcaP[1] - dcaN[1]);
      z1 = dcaP[1];
      z2 = dcaN[1];
      best = ic;
    }
  }
  const float p1x = c1x[best], p1y = c1y[best];
  const float p2x = c2x[best], p2y = c2y[best];
  const float convX = 0.5f * (p1x + p2x), convY = 0.5f * (p1y + p2y);
  const float convZ = 0.5f * (z1 + z2);
  v.vtx = {convX, convY, convZ};
  v.pca = std::hypot(pcaXY, z1 - z2);

  auto arcTo = [](o2::track::TrackAuxPar const& h, float xTrk, float yTrk, float px_, float py_) {
    const float th0 = std::atan2(yTrk - h.yC, xTrk - h.xC);
    const float thv = std::atan2(py_ - h.yC, px_ - h.xC);
    const auto dth = RecoDecay::constrainAngle<float>(thv - th0, -o2::constants::math::PI);
    return std::fabs(h.rC * dth);
  };
  const auto gP = tP.getXYZGlo();
  const auto gN = tN.getXYZGlo();
  const float bzT = bzkG / 10.f; // o2-linter: disable=magic-number (kGauss -> Tesla)
  const auto pP = getPropMomentumFromTrackHelix(arcTo(h1, gP.X(), gP.Y(), p1x, p1y), tP, h1, bzT);
  const auto pN = getPropMomentumFromTrackHelix(arcTo(h2, gN.X(), gN.Y(), p2x, p2y), tN, h2, bzT);
  v.mom = {pP[0] + pN[0], pP[1] + pN[1], pP[2] + pN[2]};
  v.legMom = {pP, pN};

  const float pMag = std::sqrt(v.mom[0] * v.mom[0] + v.mom[1] * v.mom[1] + v.mom[2] * v.mom[2]);
  const float vMag = std::sqrt(convX * convX + convY * convY + convZ * convZ);
  if (pMag < kTiny || vMag < 1e-3f) { // o2-linter: disable=magic-number (degenerate vertex)
    return v;
  }
  v.cospa = std::clamp((convX * v.mom[0] + convY * v.mom[1] + convZ * v.mom[2]) / (pMag * vMag), -1.f, 1.f);
  const float ptSum = std::hypot(v.mom[0], v.mom[1]);
  if (ptSum > kTiny) {
    v.dcaxy = std::fabs(convX * v.mom[1] - convY * v.mom[0]) / ptSum;
    const float t = -(convX * v.mom[0] + convY * v.mom[1]) / (ptSum * ptSum);
    v.dcaz = std::fabs(convZ + t * v.mom[2]);
  }
  const float e0 = std::sqrt(pP[0] * pP[0] + pP[1] * pP[1] + pP[2] * pP[2] + kMe * kMe);
  const float e1 = std::sqrt(pN[0] * pN[0] + pN[1] * pN[1] + pN[2] * pN[2] + kMe * kMe);
  v.mee = std::sqrt(std::max(0.f, (e0 + e1) * (e0 + e1) - pMag * pMag));
  v.score = analyticScoreV0(v.cospa, v.pca, scoreWeight);
  v.ok = true;
  return v;
}

struct PhotonLikeCuts {
  float maxDca{1.5f}, minR{1.f}, maxR{90.f}, minCosPA{0.98f}, maxMee{0.1f}, maxDcaXY{3.f}, maxDcaZ{999.f}; // o2-linter: disable=magic-number (defaults = crosspair_group)
};

inline bool isPhotonLike(AnalyticV0 const& v, PhotonLikeCuts const& c)
{
  return v.ok && v.pca < c.maxDca && v.rxy() > c.minR && v.rxy() < c.maxR &&
         v.cospa > c.minCosPA && v.mee < c.maxMee && v.dcaxy < c.maxDcaXY && v.dcaz < c.maxDcaZ;
}

struct AltPairing {
  bool fitted{false};
  float dca{999.f}, rxy{-1.f}, cospa{-2.f}, mee{-1.f}, dcaz{999.f}, dcaxy{999.f};
};

//_______________________________________________________________________
/// \brief Build the candidate with the DCAFitter
/// \param pos e+ leg
/// \param neg e- leg
/// \param bzT field in Tesla
/// \param scoreWeight w of analyticScoreV0
/// \return the candidate; ok = false when the fitter does not converge
inline AnalyticV0 buildFitterV0(o2::track::TrackParCov tPos, o2::track::TrackParCov tEle, float bzT, float scoreWeight)
{
  constexpr float kMe = 0.000510999f; // electron mass, GeV/c^2
  constexpr float kTiny = 1e-12f;
  static o2::vertexing::DCAFitterN<2> fitter;
  AnalyticV0 v;
  if (std::fabs(bzT) < 1e-9f) { // o2-linter: disable=magic-number (no field)
    return v;
  }
  fitter.setBz(10.f * bzT); // o2-linter: disable=magic-number (Tesla -> kGauss)
  fitter.setPropagateToPCA(true);
  fitter.setMaxR(200.f);            // o2-linter: disable=magic-number (fitter range, cm)
  fitter.setMaxDZIni(8.f);          // o2-linter: disable=magic-number (as the SVertexer photon tune)
  fitter.setMaxDXYIni(8.f);         // o2-linter: disable=magic-number (as the SVertexer photon tune)
  fitter.setMaxChi2(1e9f);          // o2-linter: disable=magic-number (no chi2 cut; the distance is judged afterwards)
  fitter.setMinParamChange(1e-3f);  // o2-linter: disable=magic-number (fitter convergence)
  fitter.setMinRelChi2Change(0.9f); // o2-linter: disable=magic-number (fitter convergence)
  fitter.setUseAbsDCA(true);
  fitter.setCollinear(true);
  if (fitter.process(tPos, tEle) == 0) {
    return v;
  }
  if (!fitter.isPropagateTracksToVertexDone(0) && !fitter.propagateTracksToVertex(0)) {
    return v;
  }
  const auto& vtx = fitter.getPCACandidatePos(0);
  std::array<float, 3> x0{}, x1{}, p0{}, p1{};
  fitter.getTrack(0, 0).getXYZGlo(x0);
  fitter.getTrack(1, 0).getXYZGlo(x1);
  if (!fitter.getTrack(0, 0).getPxPyPzGlo(p0) || !fitter.getTrack(1, 0).getPxPyPzGlo(p1)) {
    return v;
  }
  v.vtx = {static_cast<float>(vtx[0]), static_cast<float>(vtx[1]), static_cast<float>(vtx[2])};
  v.mom = {p0[0] + p1[0], p0[1] + p1[1], p0[2] + p1[2]};
  v.legMom = {p0, p1};
  v.pca = std::sqrt((x0[0] - x1[0]) * (x0[0] - x1[0]) + (x0[1] - x1[1]) * (x0[1] - x1[1]) + (x0[2] - x1[2]) * (x0[2] - x1[2]));
  const float pMag = std::sqrt(v.mom[0] * v.mom[0] + v.mom[1] * v.mom[1] + v.mom[2] * v.mom[2]);
  const float vMag = std::sqrt(v.vtx[0] * v.vtx[0] + v.vtx[1] * v.vtx[1] + v.vtx[2] * v.vtx[2]);
  if (pMag < kTiny || vMag < kTiny) {
    return v;
  }
  v.cospa = std::clamp((v.vtx[0] * v.mom[0] + v.vtx[1] * v.mom[1] + v.vtx[2] * v.mom[2]) / (pMag * vMag), -1.f, 1.f);
  const float ptSum = std::hypot(v.mom[0], v.mom[1]);
  if (ptSum > kTiny) {
    v.dcaxy = std::fabs(v.vtx[0] * v.mom[1] - v.vtx[1] * v.mom[0]) / ptSum;
    const float t = -(v.vtx[0] * v.mom[0] + v.vtx[1] * v.mom[1]) / (ptSum * ptSum);
    v.dcaz = std::fabs(v.vtx[2] + t * v.mom[2]);
  }
  const float e0 = std::sqrt(p0[0] * p0[0] + p0[1] * p0[1] + p0[2] * p0[2] + kMe * kMe);
  const float e1 = std::sqrt(p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2] + kMe * kMe);
  v.mee = std::sqrt(std::max(0.f, (e0 + e1) * (e0 + e1) - pMag * pMag));
  v.score = analyticScoreV0(v.cospa, v.pca, scoreWeight);
  v.ok = true;
  return v;
}

//_______________________________________________________________________
/// \brief buildFitterV0 for two legs given as (point, momentum, charge)
inline AnalyticV0 buildFitterV0(LegParam const& pos, LegParam const& neg, float bzT, float scoreWeight)
{
  constexpr float kFitterMaxR = 200.f;                                        // o2-linter: disable=magic-number (the fitter's own setMaxR)
  const AnalyticV0 seed = buildAnalyticV0(pos, neg, 10.f * bzT, scoreWeight); // o2-linter: disable=magic-number (Tesla -> kGauss)
  if (!seed.ok || seed.rxy() > kFitterMaxR) {
    return AnalyticV0{};
  }
  return buildFitterV0(trackFromLeg(pos), trackFromLeg(neg), bzT, scoreWeight);
}

//_______________________________________________________________________
/// \brief Move a track in the global frame
/// \param t the track; its local frame is rotated by alpha, so a global shift enters x and y rotated
/// \param dx shift in global x (cm)
/// \param dy shift in global y (cm)
/// \param dz shift in global z (cm)
inline void translateTrackParCov(o2::track::TrackParCov& t, float dx, float dy, float dz)
{
  const float ca = std::cos(t.getAlpha()), sa = std::sin(t.getAlpha());
  t.setX(t.getX() + ca * dx + sa * dy);
  t.setY(t.getY() - sa * dx + ca * dy);
  t.setZ(t.getZ() + dz);
}

//_______________________________________________________________________
/// \brief The DCAFitter alternative pairing
/// \param pos e+ leg
/// \param neg e- leg
/// \param bzT field in Tesla
/// \return AltPairing; fitted = false when the fitter does not converge
inline AltPairing evaluateAltPairingFitter(LegParam const& pos, LegParam const& neg, float bzT)
{
  const AnalyticV0 v = buildFitterV0(pos, neg, bzT, 0.5f); // o2-linter: disable=magic-number (score not used here)
  AltPairing r;
  if (!v.ok) {
    return r;
  }
  r.fitted = true;
  r.dca = v.pca;
  r.rxy = v.rxy();
  r.cospa = v.cospa;
  r.mee = v.mee;
  r.dcaxy = v.dcaxy;
  r.dcaz = v.dcaz;
  return r;
}

} // namespace o2::pwgem::photonmeson

#endif // PWGEM_PHOTONMESON_UTILS_ANALYTICV0_H_
