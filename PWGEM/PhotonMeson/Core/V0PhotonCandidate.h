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

#include <cmath>

struct V0PhotonCandidate {

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
  float centrality;
  float pca;

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
    float dca_x_v0_to_pv = (v0.GetX() - v0.GetPx() * cospa * length / v0mom) - collision.posX();
    float dca_y_v0_to_pv = (v0.GetY() - v0.GetPy() * cospa * length / v0mom) - collision.posY();
    float dca_z_v0_to_pv = (v0.GetZ() - v0.GetPz() * cospa * length / v0mom) - collision.posZ();

    float sign_tmp = (dca_x_v0_to_pv * dca_y_v0_to_pv > 0.f) ? +1.f : -1.f;
    dcaXYV0ToPV = RecoDecay::sqrtSumOfSquares(dca_x_v0_to_pv, dca_y_v0_to_pv) * sign_tmp;
    dcaZV0ToPV = dca_z_v0_to_pv;

    alpha = v0_alpha(posPx, posPy, posPz, elePx, elePy, elePz);
    qt = v0_qt(posPx, posPy, posPz, elePx, elePy, elePz);
    int posSign = (pos.GetQ() > 0) - (pos.GetQ() < 0);
    int eleSign = (ele.GetQ() > 0) - (ele.GetQ() < 0);
    phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(posPx, posPy, posPz, elePx, elePy, elePz, posSign, eleSign, d_bz);
    psipair = o2::aod::pwgem::dilepton::utils::pairutil::getPsiPair(posPx, posPy, posPz, elePx, elePy, elePz);

    centrality = collision.centFT0M();
  }

  // Constructor for V0PhotonCut
  V0PhotonCandidate(const auto& v0, const auto& pos, const auto& ele, float cent, float d_bz) : centrality(cent)
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
  float GetPosPt() const { return posPT; }
  float GetElePt() const { return elePT; }
  float GetChi2NDF() const { return chi2ndf; }
  float GetDcaXYToPV() const { return dcaXYV0ToPV; }
  float GetDcaZToPV() const { return dcaZV0ToPV; }
  float GetAlpha() const { return alpha; }
  float GetQt() const { return qt; }
  float GetPhiV() const { return phiv; }
  float GetPsiPair() const { return psipair; }
  float GetCosPA() const { return cospa; }
  float GetPx() const { return px; }
  float GetPy() const { return py; }
  float GetPz() const { return pz; }
  float GetPt() const { return pT; }
  float GetPosPx() const { return posPx; }
  float GetPosPy() const { return posPy; }
  float GetPosPz() const { return posPz; }
  float GetElePx() const { return elePx; }
  float GetElePy() const { return elePy; }
  float GetElePz() const { return elePz; }
  float GetCent() const { return centrality; }
  float GetPCA() const { return pca; }
};

#endif // PWGEM_PHOTONMESON_CORE_V0PHOTONCANDIDATE_H_
