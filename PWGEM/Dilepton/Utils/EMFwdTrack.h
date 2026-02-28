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

/// \class to store minimal fwdtrack info
/// \author daiki.sekihata@cern.ch

#ifndef PWGEM_DILEPTON_UTILS_EMFWDTRACK_H_
#define PWGEM_DILEPTON_UTILS_EMFWDTRACK_H_

namespace o2::aod::pwgem::dilepton::utils
{
class EMFwdTrack
{
 public:
  EMFwdTrack(float pt, float eta, float phi, float mass, int8_t charge, float dcaX, float dcaY, float cXX, float cXY, float cYY)
  {
    fPt = pt;
    fEta = eta;
    fPhi = phi;
    fMass = mass;
    fCharge = charge;
    fDCAx = dcaX;
    fDCAy = dcaY;
    fCXX = cXX;
    fCXY = cXY;
    fCYY = cYY;
  }

  ~EMFwdTrack() {}

  float pt() const { return fPt; }
  float eta() const { return fEta; }
  float phi() const { return fPhi; }
  float mass() const { return fMass; }
  int8_t sign() const { return fCharge; }
  float fwdDcaX() const { return fDCAx; }
  float fwdDcaY() const { return fDCAy; }
  float fwdDcaXY() const { return std::sqrt(std::pow(fDCAx, 2) + std::pow(fDCAy, 2)); }
  float p() const { return fPt * std::cosh(fEta); }
  float px() const { return fPt * std::cos(fPhi); }
  float py() const { return fPt * std::sin(fPhi); }
  float pz() const { return fPt * std::sinh(fEta); }
  float e() const { return std::hypot(fPt * std::cosh(fEta), fMass); } // e2 = p2 + m2
  float signed1Pt() const { return fCharge * 1.f / fPt; }

  float cXX() const { return fCXX; }
  float cXY() const { return fCXY; }
  float cYY() const { return fCYY; }

 protected:
  float fPt;
  float fEta;
  float fPhi;
  float fMass;
  int8_t fCharge;
  float fDCAx;
  float fDCAy;
  float fCXX;
  float fCXY;
  float fCYY;
};

class EMFwdTrackWithCov : public EMFwdTrack
{
 public:
  EMFwdTrackWithCov(float pt, float eta, float phi, float mass, int8_t charge, float dcaX, float dcaY, float cXX, float cXY, float cYY,
                    float X = 0.f, float Y = 0.f, float Z = 0.f, float tgl = 0.f,
                    float cPhiX = 0.f, float cPhiY = 0.f, float cPhiPhi = 0.f,
                    float cTglX = 0.f, float cTglY = 0.f, float cTglPhi = 0.f, float cTglTgl = 0.f,
                    float c1PtX = 0.f, float c1PtY = 0.f, float c1PtPhi = 0.f, float c1PtTgl = 0.f, float c1Pt21Pt2 = 0.f, float chi2 = 0.f) : EMFwdTrack(pt, eta, phi, mass, charge, dcaX, dcaY, cXX, cXY, cYY)
  {
    fX = X;
    fY = Y;
    fZ = Z;
    fTgl = tgl;
    fCPhiX = cPhiX;
    fCPhiY = cPhiY;
    fCPhiPhi = cPhiPhi;
    fCTglX = cTglX;
    fCTglY = cTglY;
    fCTglPhi = cTglPhi;
    fCTglTgl = cTglTgl;
    fC1PtX = c1PtX;
    fC1PtY = c1PtY;
    fC1PtPhi = c1PtPhi;
    fC1PtTgl = c1PtTgl;
    fC1Pt21Pt2 = c1Pt21Pt2;
    fChi2 = chi2;
  }

  ~EMFwdTrackWithCov() {}

  float x() const { return fX; }
  float y() const { return fY; }
  float z() const { return fZ; }
  float tgl() const { return fTgl; }
  float cPhiX() const { return fCPhiX; }
  float cPhiY() const { return fCPhiY; }
  float cPhiPhi() const { return fCPhiPhi; }
  float cTglX() const { return fCTglX; }
  float cTglY() const { return fCTglY; }
  float cTglPhi() const { return fCTglPhi; }
  float cTglTgl() const { return fCTglTgl; }
  float c1PtX() const { return fC1PtX; }
  float c1PtY() const { return fC1PtY; }
  float c1PtPhi() const { return fC1PtPhi; }
  float c1PtTgl() const { return fC1PtTgl; }
  float c1Pt21Pt2() const { return fC1Pt21Pt2; }
  float chi2() const { return fChi2; }

 protected:
  float fX;
  float fY;
  float fZ;
  float fTgl;
  float fCPhiX;
  float fCPhiY;
  float fCPhiPhi;
  float fCTglX;
  float fCTglY;
  float fCTglPhi;
  float fCTglTgl;
  float fC1PtX;
  float fC1PtY;
  float fC1PtPhi;
  float fC1PtTgl;
  float fC1Pt21Pt2;
  float fChi2; // chi2 not chi2/ndf
};

} // namespace o2::aod::pwgem::dilepton::utils
#endif // PWGEM_DILEPTON_UTILS_EMFWDTRACK_H_
