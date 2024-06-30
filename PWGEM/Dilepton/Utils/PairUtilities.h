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
//_______________________________________________________________________
} // namespace o2::aod::pwgem::dilepton::utils::pairutil
#endif // PWGEM_DILEPTON_UTILS_PAIRUTILITIES_H_
