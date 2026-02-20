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

enum CentType : uint8_t {
  CentFT0M = 0,
  CentFT0A = 1,
  CentFT0C = 2
};

struct V0PhotonCandidate {

 public:
  // Empty Constructor
  V0PhotonCandidate() = default;
  // Set method for photonconversionbuilder
  void setPhotonCandidate(const KFParticle& v0, const KFParticle& pos, const KFParticle& ele, const auto& collision, float cospa, float psipair, float phiv, CentType centType)
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

    this->cospa = cospa;
    this->psipair = psipair;
    this->phiv = phiv;
    this->centType = centType;

    switch (centType) {
      case CentType::CentFT0A:
        cent = collision.centFT0A();
        break;
      case CentType::CentFT0C:
        cent = collision.centFT0C();
        break;
      case CentType::CentFT0M:
        cent = collision.centFT0M();
        break;
    }
  }

  // Set-Method for V0PhotonCut
  void setPhoton(const auto& v0, const auto& pos, const auto& ele, float cent, CentType centType)
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
    psipair = 999.f; // default if V0PhotonPhiVPsi table is not included
    phiv = 999.f; // default if V0PhotonPhiVPsi table is not included
    if constexpr( requires{ v0.psipair(); v0.phiv(); } ) {
      psipair = v0.psipair();
      phiv = v0.phiv();
    }
    this->cent = cent;
    this->centType = centType;
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
  float getCent() const { return cent; }
  float getPCA() const { return pca; }
  CentType getCentType() const { return centType; }

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
  float cent;
  float pca;
  CentType centType;
};

#endif // PWGEM_PHOTONMESON_CORE_V0PHOTONCANDIDATE_H_
