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

/// \commonly used to calculate pair variables
/// \author daiki.sekihata@cern.ch

#ifndef PWGEM_DILEPTON_UTILS_PAIRUTILITIES_H_
#define PWGEM_DILEPTON_UTILS_PAIRUTILITIES_H_

#include <array>
#include <vector>
#include "Math/SMatrix.h"
#include "Common/Core/RecoDecay.h"
#include "ReconstructionDataFormats/TrackFwd.h"

//_______________________________________________________________________
namespace o2::aod::pwgem::dilepton::utils::pairutil
{
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;

template <typename TDCAFitter, typename TCollision, typename TTrack>
bool isSVFound(TDCAFitter fitter, TCollision const& collision, TTrack const& t1, TTrack const& t2, float& pca, float& lxy, float& cosPA)
{
  std::array<float, 3> svpos = {0.}; // secondary vertex position
  std::array<float, 3> pvec0 = {0.};
  std::array<float, 3> pvec1 = {0.};
  auto pTrack = getTrackParCov(t1);
  auto nTrack = getTrackParCov(t2);

  int nCand = fitter.process(pTrack, nTrack);
  if (nCand != 0) {
    fitter.propagateTracksToVertex();
    const auto& vtx = fitter.getPCACandidate();
    for (int i = 0; i < 3; i++) {
      svpos[i] = vtx[i];
    }
    fitter.getTrack(0).getPxPyPzGlo(pvec0); // positive
    fitter.getTrack(1).getPxPyPzGlo(pvec1); // negative
    std::array<float, 3> pvxyz{pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]};
    // float pt = RecoDecay::sqrtSumOfSquares(pvxyz[0], pvxyz[1]);
    // LOGF(info, "pair pT = %f GeV/c at sv", pt);

    pca = std::sqrt(fitter.getChi2AtPCACandidate()); // distance between 2 legs.
    lxy = RecoDecay::sqrtSumOfSquares(svpos[0] - collision.posX(), svpos[1] - collision.posY());
    cosPA = RecoDecay::cpa(std::array{collision.posX(), collision.posY(), collision.posZ()}, svpos, pvxyz);
    return true;
  } else {
    pca = 999.f;
    lxy = 999.f;
    cosPA = 999.f;
    return false;
  }
}

// function call without cosPA
template <typename TDCAFitter, typename TCollision, typename TTrack>
bool isSVFound(TDCAFitter fitter, TCollision const& collision, TTrack const& t1, TTrack const& t2, float& pca, float& lxy)
{
  float cosPA = 999.f;
  return isSVFound(fitter, collision, t1, t2, pca, lxy, cosPA);
}
//_______________________________________________________________________
template <typename TFwdDCAFitter, typename TCollision, typename TTrack>
bool isSVFoundFwd(TFwdDCAFitter fitter, TCollision const& collision, TTrack const& t1, TTrack const& t2, float& pca, float& lxy)
{
  std::array<float, 3> svpos = {0.}; // secondary vertex position
  // std::array<float, 3> pvec0 = {0.};
  // std::array<float, 3> pvec1 = {0.};

  SMatrix5 t1pars(t1.x(), t1.y(), t1.phi(), t1.tgl(), t1.signed1Pt());
  std::vector<double> v1{t1.cXX(), t1.cXY(), t1.cYY(), t1.cPhiX(), t1.cPhiY(),
                         t1.cPhiPhi(), t1.cTglX(), t1.cTglY(), t1.cTglPhi(), t1.cTglTgl(),
                         t1.c1PtX(), t1.c1PtY(), t1.c1PtPhi(), t1.c1PtTgl(), t1.c1Pt21Pt2()};

  SMatrix55 t1covs(v1.begin(), v1.end());
  o2::track::TrackParCovFwd pTrack{t1.z(), t1pars, t1covs, t1.chi2()};

  SMatrix5 t2pars(t2.x(), t2.y(), t2.phi(), t2.tgl(), t2.signed1Pt());
  std::vector<double> v2{t2.cXX(), t2.cXY(), t2.cYY(), t2.cPhiX(), t2.cPhiY(),
                         t2.cPhiPhi(), t2.cTglX(), t2.cTglY(), t2.cTglPhi(), t2.cTglTgl(),
                         t2.c1PtX(), t2.c1PtY(), t2.c1PtPhi(), t2.c1PtTgl(), t2.c1Pt21Pt2()};
  SMatrix55 t2covs(v2.begin(), v2.end());
  o2::track::TrackParCovFwd nTrack{t2.z(), t2pars, t2covs, t2.chi2()};

  // if (abs(t1.cXX()) < 1e-6 || abs(t2.cXX()) < 1e-6) {
  //   pca = 999.f;
  //   lxy = 999.f;
  //   return false;
  // }

  int nCand = fitter.process(pTrack, nTrack);
  if (nCand != 0) {
    fitter.FwdpropagateTracksToVertex();
    const auto& vtx = fitter.getPCACandidate();
    for (int i = 0; i < 3; i++) {
      svpos[i] = vtx[i];
    }
    // fitter.getTrack(0).getPxPyPzGlo(pvec0); // positive
    // fitter.getTrack(1).getPxPyPzGlo(pvec1); // negative
    // std::array<float, 3> pvxyz{pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]};
    // float pt = RecoDecay::sqrtSumOfSquares(pvxyz[0], pvxyz[1]);
    // LOGF(info, "pair pT = %f GeV/c at sv", pt);

    pca = std::sqrt(fitter.getChi2AtPCACandidate()); // distance between 2 legs.
    lxy = RecoDecay::sqrtSumOfSquares(svpos[0] - collision.posX(), svpos[1] - collision.posY());
    return true;
  } else {
    pca = 999.f;
    lxy = 999.f;
    return false;
  }
}
//_______________________________________________________________________
inline float getPhivPair(float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg, int cpos, int cneg, float bz)
{
  // cos(phiv) = w*a /|w||a|
  // with w = u x v
  // and  a = u x z / |u x z|   , unit vector perpendicular to v12 and z-direction (magnetic field)
  // u = v12 / |v12|            , the unit vector of v12
  // v = v1 x v2 / |v1 x v2|    , unit vector perpendicular to v1 and v2

  // momentum of e+ and e- in (ax,ay,az) axis. Note that az=0 by definition.
  // vector product of pep X pem
  // std::array<float, 3> arr_pos{t1.GetPx(), t1.GetPy(), t1.GetPz()};
  // std::array<float, 3> arr_ele{t2.GetPx(), t2.GetPy(), t2.GetPz()};
  std::array<float, 3> arr_pos{pxpos, pypos, pzpos};
  std::array<float, 3> arr_ele{pxneg, pyneg, pzneg};
  std::array<double, 3> pos_x_ele{0, 0, 0};
  // LOGF(info, "Q1 = %d , Q2 = %d", cpos, cneg);

  if (cpos * cneg > 0) { // Like Sign
    if (bz < 0) {
      if (cpos > 0) {
        pos_x_ele = RecoDecay::crossProd(arr_pos, arr_ele);
      } else {
        pos_x_ele = RecoDecay::crossProd(arr_ele, arr_pos);
      }
    } else {
      if (cpos > 0) {
        pos_x_ele = RecoDecay::crossProd(arr_ele, arr_pos);
      } else {
        pos_x_ele = RecoDecay::crossProd(arr_pos, arr_ele);
      }
    }
  } else { // Unlike Sign
    if (bz > 0) {
      if (cpos > 0) {
        pos_x_ele = RecoDecay::crossProd(arr_pos, arr_ele);
      } else {
        pos_x_ele = RecoDecay::crossProd(arr_ele, arr_pos);
      }
    } else {
      if (cpos > 0) {
        pos_x_ele = RecoDecay::crossProd(arr_ele, arr_pos);
      } else {
        pos_x_ele = RecoDecay::crossProd(arr_pos, arr_ele);
      }
    }
  }

  // unit vector of pep X pem
  float vx = pos_x_ele[0] / RecoDecay::sqrtSumOfSquares(pos_x_ele[0], pos_x_ele[1], pos_x_ele[2]);
  float vy = pos_x_ele[1] / RecoDecay::sqrtSumOfSquares(pos_x_ele[0], pos_x_ele[1], pos_x_ele[2]);
  float vz = pos_x_ele[2] / RecoDecay::sqrtSumOfSquares(pos_x_ele[0], pos_x_ele[1], pos_x_ele[2]);

  // unit vector of (pep+pem)
  float ux = (pxpos + pxneg) / RecoDecay::sqrtSumOfSquares(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
  float uy = (pypos + pyneg) / RecoDecay::sqrtSumOfSquares(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);
  float uz = (pzpos + pzneg) / RecoDecay::sqrtSumOfSquares(pxpos + pxneg, pypos + pyneg, pzpos + pzneg);

  float ax = uy / std::sqrt(ux * ux + uy * uy);
  float ay = -ux / std::sqrt(ux * ux + uy * uy);

  // The third axis defined by vector product (ux,uy,uz)X(vx,vy,vz)
  float wx = uy * vz - uz * vy;
  float wy = uz * vx - ux * vz;
  // by construction, (wx,wy,wz) must be a unit vector. Measure angle between (wx,wy,wz) and (ax,ay,0).
  // The angle between them should be small if the pair is conversion. This function then returns values close to pi!
  auto clipToPM1 = [](float x) { return x < -1.f ? -1.f : (x > 1.f ? 1.f : x); };

  // if (!std::isfinite(std::acos(wx * ax + wy * ay))) {
  //   LOGF(info, "pxpos = %f, pypos = %f, pzpos = %f", pxpos, pypos, pzpos);
  //   LOGF(info, "pxneg = %f, pyneg = %f, pzneg = %f", pxneg, pyneg, pzneg);
  //   LOGF(info, "pos_x_ele[0] = %f, pos_x_ele[1] = %f, pos_x_ele[2] = %f", pos_x_ele[0], pos_x_ele[1], pos_x_ele[2]);
  //   LOGF(info, "ux = %f, uy = %f, uz = %f", ux, uy, uz);
  //   LOGF(info, "ax = %f, ay = %f", ax, ay);
  //   LOGF(info, "wx = %f, wy = %f", wx, wy);
  //   LOGF(info, "wx * ax + wy * ay = %f", wx * ax + wy * ay);
  //   LOGF(info, "std::acos(wx * ax + wy * ay) = %f", std::acos(wx * ax + wy * ay));
  //   LOGF(info, "std::acos(clipToPM1(wx * ax + wy * ay)) = %f", std::acos(clipToPM1(wx * ax + wy * ay)));
  // }

  return std::acos(clipToPM1(wx * ax + wy * ay)); // phiv in [0,pi] //cosPhiV = wx * ax + wy * ay;
}
//_______________________________________________________________________
inline float getPsiPair(float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg)
{
  auto clipToPM1 = [](float x) { return x < -1.f ? -1.f : (x > 1.f ? 1.f : x); };
  float ptot2 = RecoDecay::p2(pxpos, pypos, pzpos) * RecoDecay::p2(pxneg, pyneg, pzneg);
  float argcos = RecoDecay::dotProd(std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}) / std::sqrt(ptot2);
  float thetaPos = std::atan2(RecoDecay::sqrtSumOfSquares(pxpos, pypos), pzpos);
  float thetaNeg = std::atan2(RecoDecay::sqrtSumOfSquares(pxneg, pyneg), pzneg);
  float argsin = (thetaNeg - thetaPos) / std::acos(clipToPM1(argcos));
  return std::asin(clipToPM1(argsin));
}
//_______________________________________________________________________
inline float getOpeningAngle(float pxpos, float pypos, float pzpos, float pxneg, float pyneg, float pzneg)
{
  auto clipToPM1 = [](float x) { return x < -1.f ? -1.f : (x > 1.f ? 1.f : x); };
  float ptot2 = RecoDecay::p2(pxpos, pypos, pzpos) * RecoDecay::p2(pxneg, pyneg, pzneg);
  float argcos = RecoDecay::dotProd(std::array{pxpos, pypos, pzpos}, std::array{pxneg, pyneg, pzneg}) / std::sqrt(ptot2);
  return std::acos(clipToPM1(argcos));
}
//_______________________________________________________________________
} // namespace o2::aod::pwgem::dilepton::utils::pairutil
#endif // PWGEM_DILEPTON_UTILS_PAIRUTILITIES_H_
