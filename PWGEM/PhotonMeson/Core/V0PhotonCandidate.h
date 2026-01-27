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

/// \file V0PhotonCandidate.h
/// \brief Getter functions for V0 photon candidate properties.
/// \author Isabel Kantak <isabel.kantak@cern.ch>

#ifndef PWGEM_PHOTONMESON_CORE_V0PHOTONCANDIDATE_H_
#define PWGEM_PHOTONMESON_CORE_V0PHOTONCANDIDATE_H_

#include "PWGEM/Dilepton/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"

#include "Common/Core/RecoDecay.h"

#include <KFParticle.h>

struct V0PhotonCandidate {

 public:
  // Constructor for photonconversionbuilder
  V0PhotonCandidate(const KFParticle& v0, const KFParticle& pos, const KFParticle& ele, const auto& collision, float cospa, float d_bz) : cospa(cospa)
  {
    px = v0.GetPx();
    py = v0.GetPy();
    pz = v0.GetPz();
    pT = RecoDecay::sqrtSumOfSquares(px, py);

    posPx = pos.GetPx();
    posPy = pos.GetPy();
    posPz = pos.GetPz();
    elePx = ele.GetPx();
    elePy = ele.GetPy();
    elePz = ele.GetPz();
    posPT = RecoDecay::sqrtSumOfSquares(posPx, posPy);
    elePT = RecoDecay::sqrtSumOfSquares(elePx, elePy);

    chi2ndf = v0.GetChi2() / v0.GetNDF();
    pca = pos.GetDistanceFromParticle(ele);

    float v0mom = RecoDecay::sqrtSumOfSquares(v0.GetPx(), v0.GetPy(), v0.GetPz());
    float length = RecoDecay::sqrtSumOfSquares(v0.GetX() - collision.posX(), v0.GetY() - collision.posY(), v0.GetZ() - collision.posZ());
    float dcaXV0ToPV = (v0.GetX() - v0.GetPx() * cospa * length / v0mom) - collision.posX();
    float dcaYV0ToPV = (v0.GetY() - v0.GetPy() * cospa * length / v0mom) - collision.posY();
    float tmpSign = (dcaXV0ToPV * dcaYV0ToPV > 0.f) ? +1.f : -1.f;

    dcaXYV0ToPV = RecoDecay::sqrtSumOfSquares(dcaXV0ToPV, dcaYV0ToPV) * tmpSign;
    dcaZV0ToPV = (v0.GetZ() - v0.GetPz() * cospa * length / v0mom) - collision.posZ();

    alpha = v0_alpha(posPx, posPy, posPz, elePx, elePy, elePz);
    qt = v0_qt(posPx, posPy, posPz, elePx, elePy, elePz);
    int posSign = (pos.GetQ() > 0) - (pos.GetQ() < 0);
    int eleSign = (ele.GetQ() > 0) - (ele.GetQ() < 0);
    phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(posPx, posPy, posPz, elePx, elePy, elePz, posSign, eleSign, d_bz);
    psipair = o2::aod::pwgem::dilepton::utils::pairutil::getPsiPair(posPx, posPy, posPz, elePx, elePy, elePz);

    centFT0M = collision.centFT0M();
    centFT0C = collision.centFT0C();
    centFT0A = collision.centFT0A();
  }

  // Constructor for V0PhotonCut
  V0PhotonCandidate(const auto& v0, const auto& pos, const auto& ele, float centFT0A, float centFT0C, float centFT0M, float d_bz) : centFT0A(centFT0A), centFT0C(centFT0C), centFT0M(centFT0M)
  {
    px = v0.px();
    py = v0.py();
    pz = v0.pz();
    pT = v0.pt();

    posPx = pos.px();
    posPy = pos.py();
    posPz = pos.pz();
    elePx = ele.px();
    elePy = ele.py();
    elePz = ele.pz();
    posPT = pos.pt();
    elePT = ele.pt();

    chi2ndf = v0.chiSquareNDF();
    pca = v0.pca();

    dcaXYV0ToPV = v0.dcaXYtopv();
    dcaZV0ToPV = v0.dcaZtopv();

    cospa = v0.cospa();
    alpha = v0.alpha();
    qt = v0.qtarm();

    phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(posPx, posPy, posPz, elePx, elePy, elePz, pos.sign(), ele.sign(), d_bz);
    psipair = o2::aod::pwgem::dilepton::utils::pairutil::getPsiPair(posPx, posPy, posPz, elePx, elePy, elePz);
  }

  // Getter functions
  float getPosPt() const { return posPT; }
  float getElePt() const { return elePT; }
  float getChi2NDF() const { return chi2ndf; }
  float getDcaXYToPV() const { return dcaXYV0ToPV; }
  float getDcaZToPV() const { return dcaZV0ToPV; }
  float getAlpha() const { return alpha; }
  float getQt() const { return qt; }
  float getPhiV() const { return phiv; }
  float getPsiPair() const { return psipair; }
  float getCosPA() const { return cospa; }
  float getPx() const { return px; }
  float getPy() const { return py; }
  float getPz() const { return pz; }
  float getPt() const { return pT; }
  float getPosPx() const { return posPx; }
  float getPosPy() const { return posPy; }
  float getPosPz() const { return posPz; }
  float getElePx() const { return elePx; }
  float getElePy() const { return elePy; }
  float getElePz() const { return elePz; }
  float getCentFT0M() const { return centFT0M; }
  float getCentFT0C() const { return centFT0C; }
  float getCentFT0A() const { return centFT0A; }
  float getPCA() const { return pca; }

 private:
  float px;
  float py;
  float pz;
  float posPx;
  float posPy;
  float posPz;
  float elePx;
  float elePy;
  float elePz;
  float pT;
  float posPT;
  float elePT;
  float dcaXYV0ToPV;
  float dcaZV0ToPV;
  float alpha;
  float qt;
  float phiv;
  float psipair;
  float cospa;
  float chi2ndf;
  float centFT0A;
  float centFT0C;
  float centFT0M;
  float pca;
};

#endif // PWGEM_PHOTONMESON_CORE_V0PHOTONCANDIDATE_H_
