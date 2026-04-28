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

/// \class to store minimal track info
/// \author daiki.sekihata@cern.ch

#ifndef PWGEM_DILEPTON_UTILS_EMTRACK_H_
#define PWGEM_DILEPTON_UTILS_EMTRACK_H_

#include "Math/Vector4D.h"

namespace o2::aod::pwgem::dilepton::utils
{
class EMTrack
{
 public:
  EMTrack(float pt, float eta, float phi, float mass, int8_t charge = 0, float dcaXY = 0.f, float dcaZ = 0.f, float CYY = 0, float CZY = 0, float CZZ = 0)
  {
    fPt = pt;
    fEta = eta;
    fPhi = phi;
    fMass = mass;
    fCharge = charge;
    fDCAxy = dcaXY;
    fDCAz = dcaZ;
    fCYY = CYY;
    fCZY = CZY;
    fCZZ = CZZ;
  }

  ~EMTrack() {}

  float pt() const { return fPt; }
  float eta() const { return fEta; }
  float phi() const { return fPhi; }
  float mass() const { return fMass; }
  int8_t sign() const { return fCharge; }
  float dcaXY() const { return fDCAxy; }
  float dcaZ() const { return fDCAz; }

  float cYY() const { return fCYY; }
  float cZY() const { return fCZY; }
  float cZZ() const { return fCZZ; }

  float rapidity() const { return std::log((std::sqrt(std::pow(fMass, 2) + std::pow(fPt * std::cosh(fEta), 2)) + fPt * std::sinh(fEta)) / std::sqrt(std::pow(fMass, 2) + std::pow(fPt, 2))); }
  float p() const { return fPt * std::cosh(fEta); }
  float px() const { return fPt * std::cos(fPhi); }
  float py() const { return fPt * std::sin(fPhi); }
  float pz() const { return fPt * std::sinh(fEta); }
  float e() const { return std::hypot(fPt * std::cosh(fEta), fMass); } // e2 = p2 + m2
  float signed1Pt() const { return fCharge * 1.f / fPt; }

 protected:
  float fPt;
  float fEta;
  float fPhi;
  float fMass;
  int8_t fCharge;
  float fDCAxy;
  float fDCAz;
  float fCYY;
  float fCZY;
  float fCZZ;
};

class EMTrackWithCov : public EMTrack
{
 public:
  EMTrackWithCov(float pt, float eta, float phi, float mass, int8_t charge = 0, float dcaXY = 0.f, float dcaZ = 0.f,
                 float CYY = 0.f, float CZY = 0.f, float CZZ = 0.f,
                 float X = 0.f, float Y = 0.f, float Z = 0.f, float Alpha = 0.f, float Snp = 0.f, float Tgl = 0.f,
                 float CSnpY = 0.f, float CSnpZ = 0.f, float CSnpSnp = 0.f,
                 float CTglY = 0.f, float CTglZ = 0.f, float CTglSnp = 0.f, float CTglTgl = 0.f,
                 float C1PtY = 0.f, float C1PtZ = 0.f, float C1PtSnp = 0.f, float C1PtTgl = 0.f, float C1Pt21Pt2 = 0.f) : EMTrack(pt, eta, phi, mass, charge, dcaXY, dcaZ, CYY, CZY, CZZ)
  {
    fX = X;
    fY = Y;
    fZ = Z;
    fAlpha = Alpha;
    fSnp = Snp;
    fTgl = Tgl;
    fCYY = CYY;
    fCZY = CZY;
    fCZZ = CZZ;
    fCSnpY = CSnpY;
    fCSnpZ = CSnpZ;
    fCSnpSnp = CSnpSnp;
    fCTglY = CTglY;
    fCTglZ = CTglZ;
    fCTglSnp = CTglSnp;
    fCTglTgl = CTglTgl;
    fC1PtY = C1PtY;
    fC1PtZ = C1PtZ;
    fC1PtSnp = C1PtSnp;
    fC1PtTgl = C1PtTgl;
    fC1Pt21Pt2 = C1Pt21Pt2;
  }

  ~EMTrackWithCov() {}

  float x() const { return fX; }
  float y() const { return fY; }
  float z() const { return fZ; }
  float alpha() const { return fAlpha; }
  float snp() const { return fSnp; }
  float tgl() const { return fTgl; }

  float cSnpY() const { return fCSnpY; }
  float cSnpZ() const { return fCSnpZ; }
  float cSnpSnp() const { return fCSnpSnp; }
  float cTglY() const { return fCTglY; }
  float cTglZ() const { return fCTglZ; }
  float cTglSnp() const { return fCTglSnp; }
  float cTglTgl() const { return fCTglTgl; }
  float c1PtY() const { return fC1PtY; }
  float c1PtZ() const { return fC1PtZ; }
  float c1PtSnp() const { return fC1PtSnp; }
  float c1PtTgl() const { return fC1PtTgl; }
  float c1Pt21Pt2() const { return fC1Pt21Pt2; }

 protected:
  float fX;
  float fY;
  float fZ;
  float fAlpha;
  float fSnp;
  float fTgl;
  float fCSnpY;
  float fCSnpZ;
  float fCSnpSnp;
  float fCTglY;
  float fCTglZ;
  float fCTglSnp;
  float fCTglTgl;
  float fC1PtY;
  float fC1PtZ;
  float fC1PtSnp;
  float fC1PtTgl;
  float fC1Pt21Pt2;
};

class EMPair : public EMTrack
{
 public:
  EMPair(float pt, float eta, float phi, float mass, int8_t charge = 0) : EMTrack(pt, eta, phi, mass, charge, 0, 0, 0, 0, 0)
  {
    fPairDCA = 999.f;
    fVPos = ROOT::Math::PtEtaPhiMVector(0, 0, 0, 0);
    fVNeg = ROOT::Math::PtEtaPhiMVector(0, 0, 0, 0);
    fVx = 0.f;
    fVy = 0.f;
    fVz = 0.f;
  }

  ~EMPair() {}

  void setPairDCA(float dca) { fPairDCA = dca; }
  float getPairDCA() const { return fPairDCA; }

  void setPositiveLegPtEtaPhiM(float pt, float eta, float phi, float m)
  {
    fVPos.SetPt(pt);
    fVPos.SetEta(eta);
    fVPos.SetPhi(phi);
    fVPos.SetM(m);
  }
  void setNegativeLegPtEtaPhiM(float pt, float eta, float phi, float m)
  {
    fVNeg.SetPt(pt);
    fVNeg.SetEta(eta);
    fVNeg.SetPhi(phi);
    fVNeg.SetM(m);
  }

  void setPositiveLegPxPyPzM(float px, float py, float pz, float m)
  {
    float pt = std::sqrt(px * px + py * py);
    float eta = std::atanh(pz / std::sqrt(std::pow(px, 2) + std::pow(py, 2) + std::pow(pz, 2)));
    float phi = std::atan2(py, px);
    if (phi < 0.f) {
      phi += 2.f * M_PI;
    }

    fVPos.SetPt(pt);
    fVPos.SetEta(eta);
    fVPos.SetPhi(phi);
    fVPos.SetM(m);
  }
  void setNegativeLegPxPyPzM(float px, float py, float pz, float m)
  {
    float pt = std::sqrt(px * px + py * py);
    float eta = std::atanh(pz / std::sqrt(std::pow(px, 2) + std::pow(py, 2) + std::pow(pz, 2)));
    float phi = std::atan2(py, px);
    if (phi < 0.f) {
      phi += 2.f * M_PI;
    }

    fVNeg.SetPt(pt);
    fVNeg.SetEta(eta);
    fVNeg.SetPhi(phi);
    fVNeg.SetM(m);
  }

  ROOT::Math::PtEtaPhiMVector getPositiveLeg() const { return fVPos; }
  ROOT::Math::PtEtaPhiMVector getNegativeLeg() const { return fVNeg; }

  void setConversionPointXYZ(float x, float y, float z)
  {
    fVx = x;
    fVy = y;
    fVz = z;
  }
  float vx() const { return fVx; }
  float vy() const { return fVy; }
  float vz() const { return fVz; }
  float v0radius() const { return std::sqrt(std::pow(fVx, 2) + std::pow(fVy, 2)); }
  float eta_cp() const { return std::atanh(fVz / std::sqrt(std::pow(fVx, 2) + std::pow(fVy, 2) + std::pow(fVz, 2))); }
  float phi_cp() const { return std::atan2(fVy, fVx); }

 protected:
  float fPairDCA;
  ROOT::Math::PtEtaPhiMVector fVPos;
  ROOT::Math::PtEtaPhiMVector fVNeg;

  // only for photon conversion point
  float fVx;
  float fVy;
  float fVz;
};

} // namespace o2::aod::pwgem::dilepton::utils
#endif // PWGEM_DILEPTON_UTILS_EMTRACK_H_
