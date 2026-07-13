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
  void setPhotonCandidate(const KFParticle& v0DecayVtx, const KFParticle& v0PV, const TTrack& pos, const KFParticle& posDecayVtx, const TTrack& ele, const KFParticle& eleDecayVtx, const auto& collision, float cospa_, float cospaRZ_, float cospaXY_, float psipair_, float phiv_, CentType centType_, auto posdcaXY_, auto eledcaXY_, auto posdcaZ_, auto eledcaZ_)
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

    chi2ndf = v0PV.GetChi2() / v0PV.GetNDF();
    pca = posDecayVtx.GetDistanceFromParticle(eleDecayVtx);
    eta = RecoDecay::eta(std::array{px, py, pz});
    posEta = RecoDecay::eta(std::array{posPx, posPy, posPz});
    eleEta = RecoDecay::eta(std::array{elePx, elePy, elePz});

    float v0mom = RecoDecay::sqrtSumOfSquares(v0DecayVtx.GetPx(), v0DecayVtx.GetPy(), v0DecayVtx.GetPz());
    float length = RecoDecay::sqrtSumOfSquares(v0DecayVtx.GetX() - collision.posX(), v0DecayVtx.GetY() - collision.posY(), v0DecayVtx.GetZ() - collision.posZ());
    float dcaXV0ToPV = (v0DecayVtx.GetX() - v0DecayVtx.GetPx() * cospa_ * length / v0mom) - collision.posX();
    float dcaYV0ToPV = (v0DecayVtx.GetY() - v0DecayVtx.GetPy() * cospa_ * length / v0mom) - collision.posY();
    float tmpSign = (dcaXV0ToPV * dcaYV0ToPV > 0.f) ? +1.f : -1.f;

    dcaXYV0ToPV = RecoDecay::sqrtSumOfSquares(dcaXV0ToPV, dcaYV0ToPV) * tmpSign;
    dcaZV0ToPV = (v0DecayVtx.GetZ() - v0DecayVtx.GetPz() * cospa_ * length / v0mom) - collision.posZ();

    alpha = v0_alpha(posPx, posPy, posPz, elePx, elePy, elePz);
    qt = v0_qt(posPx, posPy, posPz, elePx, elePy, elePz);

    cospa = cospa_;
    cospaRZ = cospaRZ_;
    cospaXY = cospaXY_;
    psipair = psipair_;
    phiv = phiv_;
    centType = centType_;
    posdcaXY = posdcaXY_;
    eledcaXY = eledcaXY_;
    posdcaZ = posdcaZ_;
    eledcaZ = eledcaZ_;

    switch (centType_) {
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
  void setPhoton(const auto& v0, const auto& pos, const auto& ele, float cent_, CentType centType_)
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
    cent = cent_;
    centType = centType_;
  }

  // Getter functions
  [[nodiscard]] float getPosPt() const { return posPT; }
  [[nodiscard]] float getElePt() const { return elePT; }
  [[nodiscard]] float getChi2NDF() const { return chi2ndf; }
  [[nodiscard]] float getDcaXYToPV() const { return dcaXYV0ToPV; }
  [[nodiscard]] float getDcaZToPV() const { return dcaZV0ToPV; }
  [[nodiscard]] float getAlpha() const { return alpha; }
  [[nodiscard]] float getQt() const { return qt; }
  [[nodiscard]] float getPhiV() const { return phiv; }
  [[nodiscard]] float getPsiPair() const { return psipair; }
  [[nodiscard]] float getCosPA() const { return cospa; }
  [[nodiscard]] float getCosPARZ() const { return cospaRZ; }
  [[nodiscard]] float getCosPAXY() const { return cospaXY; }
  [[nodiscard]] float getEta() const { return eta; }
  [[nodiscard]] float getPosEta() const { return posEta; }
  [[nodiscard]] float getEleEta() const { return eleEta; }
  [[nodiscard]] float getConversionPointX() const { return conversionPointx; }
  [[nodiscard]] float getConversionPointY() const { return conversionPointy; }
  [[nodiscard]] float getConversionPointZ() const { return conversionPointz; }
  [[nodiscard]] float getPx() const { return px; }
  [[nodiscard]] float getPy() const { return py; }
  [[nodiscard]] float getPz() const { return pz; }
  [[nodiscard]] float getPt() const { return pT; }
  [[nodiscard]] float getPosPx() const { return posPx; }
  [[nodiscard]] float getPosPy() const { return posPy; }
  [[nodiscard]] float getPosPz() const { return posPz; }
  [[nodiscard]] float getElePx() const { return elePx; }
  [[nodiscard]] float getElePy() const { return elePy; }
  [[nodiscard]] float getElePz() const { return elePz; }
  [[nodiscard]] float getPosDcaXY() const { return posdcaXY; }
  [[nodiscard]] float getPosDcaZ() const { return posdcaZ; }
  [[nodiscard]] float getEleDcaXY() const { return eledcaXY; }
  [[nodiscard]] float getEleDcaZ() const { return eledcaZ; }
  [[nodiscard]] float getPosTPCNClsShared() const { return posTPCNClsShared; }
  [[nodiscard]] float getPosTPCNClsFindable() const { return posTPCNClsFindable; }
  [[nodiscard]] float getPosTPCNClsFindableMinusShared() const { return posTPCNClsFindableMinusShared; }
  [[nodiscard]] float getPosTPCNClsFindableMinusCrossedRows() const { return posTPCNClsFindableMinusCrossedRows; }
  [[nodiscard]] float getPosTPCChi2NCl() const { return posTPCChi2NCl; }
  [[nodiscard]] float getPosTPCSignal() const { return posTPCSignal; }
  [[nodiscard]] float getPosITSClusterSizes() const { return posITSClusterSizes; }
  [[nodiscard]] float getEleTPCNClsShared() const { return eleTPCNClsShared; }
  [[nodiscard]] float getEleTPCNClsFindable() const { return eleTPCNClsFindable; }
  [[nodiscard]] float getEleTPCNClsFindableMinusShared() const { return eleTPCNClsFindableMinusShared; }
  [[nodiscard]] float getEleTPCNClsFindableMinusCrossedRows() const { return eleTPCNClsFindableMinusCrossedRows; }
  [[nodiscard]] float getEleTPCChi2NCl() const { return eleTPCChi2NCl; }
  [[nodiscard]] float getEleTPCSignal() const { return eleTPCSignal; }
  [[nodiscard]] float getEleITSClusterSizes() const { return eleITSClusterSizes; }
  [[nodiscard]] float getCent() const { return cent; }
  [[nodiscard]] float getPCA() const { return pca; }
  [[nodiscard]] CentType getCentType() const { return centType; }

 private:
  float conversionPointx{0.f};
  float conversionPointy{0.f};
  float conversionPointz{0.f};
  float px{0.f};
  float py{0.f};
  float pz{0.f};
  float posPx{0.f};
  float posPy{0.f};
  float posPz{0.f};
  float elePx{0.f};
  float elePy{0.f};
  float elePz{0.f};
  float pT{0.f};
  float posPT{0.f};
  float elePT{0.f};
  float dcaXYV0ToPV{0.f};
  float dcaZV0ToPV{0.f};
  float alpha{0.f};
  float qt{0.f};
  float phiv{0.f};
  float psipair{0.f};
  float cospa{0.f};
  float cospaRZ{0.f};
  float cospaXY{0.f};
  float chi2ndf{0.f};
  float cent{0.f};
  float pca{0.f};
  float eta{0.f};
  float posEta{0.f};
  float eleEta{0.f};
  float posdcaXY{0.f};
  float posdcaZ{0.f};
  float eledcaXY{0.f};
  float eledcaZ{0.f};
  float posTPCNClsShared{0.f};
  float posTPCNClsFindable{0.f};
  float posTPCNClsFindableMinusShared{0.f};
  float posTPCNClsFindableMinusCrossedRows{0.f};
  float posTPCChi2NCl{0.f};
  float posTPCSignal{0.f};
  float posITSClusterSizes{0.f};
  float eleTPCNClsShared{0.f};
  float eleTPCNClsFindable{0.f};
  float eleTPCNClsFindableMinusShared{0.f};
  float eleTPCNClsFindableMinusCrossedRows{0.f};
  float eleTPCChi2NCl{0.f};
  float eleTPCSignal{0.f};
  float eleITSClusterSizes{0.f};
  CentType centType{CentFT0C};
};

#endif // PWGEM_PHOTONMESON_CORE_V0PHOTONCANDIDATE_H_
