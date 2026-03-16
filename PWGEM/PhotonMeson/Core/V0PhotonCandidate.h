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

#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"

#include "Common/Core/RecoDecay.h"

#include <KFParticle.h>

#include <array>
#include <cstdint>

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
  template <class TTrack>
  void setPhotonCandidate(const KFParticle& v0DecayVtx, const KFParticle& v0PV, const TTrack& pos, const KFParticle& posDecayVtx, const TTrack& ele, const KFParticle& eleDecayVtx, const auto& collision, float cospa, float cospaRZ, float cospaXY, float psipair, float phiv, CentType centType, auto posdcaXY, auto eledcaXY, auto posdcaZ, auto eledcaZ)
  {
    conversionPointx = v0DecayVtx.GetX();
    conversionPointy = v0DecayVtx.GetY();
    conversionPointz = v0DecayVtx.GetZ();
    px = v0PV.GetPx();
    py = v0PV.GetPy();
    pz = v0PV.GetPz();
    pT = RecoDecay::sqrtSumOfSquares(px, py);

    posPx = posDecayVtx.GetPx();
    posPy = posDecayVtx.GetPy();
    posPz = posDecayVtx.GetPz();
    elePx = eleDecayVtx.GetPx();
    elePy = eleDecayVtx.GetPy();
    elePz = eleDecayVtx.GetPz();
    posPT = RecoDecay::sqrtSumOfSquares(posPx, posPy);
    elePT = RecoDecay::sqrtSumOfSquares(elePx, elePy);
    posEta = RecoDecay::eta(std::array{posPx, posPy, posPz});
    eleEta = RecoDecay::eta(std::array{elePx, elePy, elePz});

    posTPCNClsShared = pos.tpcNClsShared();
    posTPCNClsFindable = pos.tpcNClsFindable();
    posTPCNClsFindableMinusShared = pos.tpcNClsFindableMinusFound();
    posTPCNClsFindableMinusCrossedRows = pos.tpcNClsFindableMinusCrossedRows();
    posTPCChi2NCl = pos.tpcChi2NCl();
    posTPCSignal = pos.tpcSignal();
    posITSClusterSizes = pos.itsClusterSizes();
    eleTPCNClsShared = ele.tpcNClsShared();
    eleTPCNClsFindable = ele.tpcNClsFindable();
    eleTPCNClsFindableMinusShared = ele.tpcNClsFindableMinusFound();
    eleTPCNClsFindableMinusCrossedRows = ele.tpcNClsFindableMinusCrossedRows();
    eleTPCChi2NCl = ele.tpcChi2NCl();
    eleTPCSignal = ele.tpcSignal();
    eleITSClusterSizes = ele.itsClusterSizes();

    chi2ndf = v0DecayVtx.GetChi2() / v0DecayVtx.GetNDF();
    pca = posDecayVtx.GetDistanceFromParticle(eleDecayVtx);
    eta = RecoDecay::eta(std::array{px, py, pz});
    posEta = RecoDecay::eta(std::array{posPx, posPy, posPz});
    eleEta = RecoDecay::eta(std::array{elePx, elePy, elePz});

    float v0mom = RecoDecay::sqrtSumOfSquares(v0DecayVtx.GetPx(), v0DecayVtx.GetPy(), v0DecayVtx.GetPz());
    float length = RecoDecay::sqrtSumOfSquares(v0DecayVtx.GetX() - collision.posX(), v0DecayVtx.GetY() - collision.posY(), v0DecayVtx.GetZ() - collision.posZ());
    float dcaXV0ToPV = (v0DecayVtx.GetX() - v0DecayVtx.GetPx() * cospa * length / v0mom) - collision.posX();
    float dcaYV0ToPV = (v0DecayVtx.GetY() - v0DecayVtx.GetPy() * cospa * length / v0mom) - collision.posY();
    float tmpSign = (dcaXV0ToPV * dcaYV0ToPV > 0.f) ? +1.f : -1.f;

    dcaXYV0ToPV = RecoDecay::sqrtSumOfSquares(dcaXV0ToPV, dcaYV0ToPV) * tmpSign;
    dcaZV0ToPV = (v0DecayVtx.GetZ() - v0DecayVtx.GetPz() * cospa * length / v0mom) - collision.posZ();

    alpha = v0_alpha(posPx, posPy, posPz, elePx, elePy, elePz);
    qt = v0_qt(posPx, posPy, posPz, elePx, elePy, elePz);

    this->cospa = cospa;
    this->cospaRZ = cospaRZ;
    this->cospaXY = cospaXY;
    this->psipair = psipair;
    this->phiv = phiv;
    this->centType = centType;
    this->posdcaXY = posdcaXY;
    this->eledcaXY = eledcaXY;
    this->posdcaZ = posdcaZ;
    this->eledcaZ = eledcaZ;

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
    conversionPointx = v0.vx();
    conversionPointy = v0.vy();
    conversionPointz = v0.vz();
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
    posEta = pos.eta();
    eleEta = ele.eta();
    posdcaXY = pos.dcaXY();
    posdcaZ = pos.dcaZ();
    eledcaXY = ele.dcaXY();
    eledcaZ = ele.dcaZ();
    posTPCNClsShared = pos.tpcNClsShared();
    posTPCNClsFindable = pos.tpcNClsFindable();
    posTPCNClsFindableMinusShared = pos.tpcNClsFindableMinusFound();
    posTPCNClsFindableMinusCrossedRows = pos.tpcNClsFindableMinusCrossedRows();
    posTPCChi2NCl = pos.tpcChi2NCl();
    posTPCSignal = pos.tpcSignal();
    posITSClusterSizes = pos.itsClusterSizes();
    eleTPCNClsShared = ele.tpcNClsShared();
    eleTPCNClsFindable = ele.tpcNClsFindable();
    eleTPCNClsFindableMinusShared = ele.tpcNClsFindableMinusFound();
    eleTPCNClsFindableMinusCrossedRows = ele.tpcNClsFindableMinusCrossedRows();
    eleTPCChi2NCl = ele.tpcChi2NCl();
    eleTPCSignal = ele.tpcSignal();
    eleITSClusterSizes = ele.itsClusterSizes();

    chi2ndf = v0.chiSquareNDF();
    pca = v0.pca();
    eta = v0.eta();

    dcaXYV0ToPV = v0.dcaXYtopv();
    dcaZV0ToPV = v0.dcaZtopv();

    cospa = v0.cospa();
    cospaRZ = v0.cospaRZ();
    cospaXY = v0.cospaXY();
    alpha = v0.alpha();
    qt = v0.qtarm();
    psipair = 999.f; // default if V0PhotonPhiVPsi table is not included
    phiv = 999.f;    // default if V0PhotonPhiVPsi table is not included
    if constexpr (requires { v0.psipair(); v0.phiv(); }) {
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
  float getCosPARZ() const { return cospaRZ; }
  float getCosPAXY() const { return cospaXY; }
  float getEta() const { return eta; }
  float getPosEta() const { return posEta; }
  float getEleEta() const { return eleEta; }
  float getConversionPointX() const { return conversionPointx; }
  float getConversionPointY() const { return conversionPointy; }
  float getConversionPointZ() const { return conversionPointz; }
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
  float getPosDcaXY() const { return posdcaXY; }
  float getPosDcaZ() const { return posdcaZ; }
  float getEleDcaXY() const { return eledcaXY; }
  float getEleDcaZ() const { return eledcaZ; }
  float getPosTPCNClsShared() const { return posTPCNClsShared; }
  float getPosTPCNClsFindable() const { return posTPCNClsFindable; }
  float getPosTPCNClsFindableMinusShared() const { return posTPCNClsFindableMinusShared; }
  float getPosTPCNClsFindableMinusCrossedRows() const { return posTPCNClsFindableMinusCrossedRows; }
  float getPosTPCChi2NCl() const { return posTPCChi2NCl; }
  float getPosTPCSignal() const { return posTPCSignal; }
  float getPosITSClusterSizes() const { return posITSClusterSizes; }
  float getEleTPCNClsShared() const { return eleTPCNClsShared; }
  float getEleTPCNClsFindable() const { return eleTPCNClsFindable; }
  float getEleTPCNClsFindableMinusShared() const { return eleTPCNClsFindableMinusShared; }
  float getEleTPCNClsFindableMinusCrossedRows() const { return eleTPCNClsFindableMinusCrossedRows; }
  float getEleTPCChi2NCl() const { return eleTPCChi2NCl; }
  float getEleTPCSignal() const { return eleTPCSignal; }
  float getEleITSClusterSizes() const { return eleITSClusterSizes; }
  float getCent() const { return cent; }
  float getPCA() const { return pca; }
  CentType getCentType() const { return centType; }

 private:
  float conversionPointx;
  float conversionPointy;
  float conversionPointz;
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
  float cospaRZ;
  float cospaXY;
  float chi2ndf;
  float cent;
  float pca;
  float eta;
  float posEta;
  float eleEta;
  float posdcaXY;
  float posdcaZ;
  float eledcaXY;
  float eledcaZ;
  float posTPCNClsShared;
  float posTPCNClsFindable;
  float posTPCNClsFindableMinusShared;
  float posTPCNClsFindableMinusCrossedRows;
  float posTPCChi2NCl;
  float posTPCSignal;
  float posITSClusterSizes;
  float eleTPCNClsShared;
  float eleTPCNClsFindable;
  float eleTPCNClsFindableMinusShared;
  float eleTPCNClsFindableMinusCrossedRows;
  float eleTPCChi2NCl;
  float eleTPCSignal;
  float eleITSClusterSizes;
  CentType centType;
};

#endif // PWGEM_PHOTONMESON_CORE_V0PHOTONCANDIDATE_H_
