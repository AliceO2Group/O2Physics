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

#include <vector>

namespace o2::aod::pwgem::dilepton::utils
{
class EMFwdTrack
{
 public:
  EMFwdTrack(int dfId, int globalId, int collisionId, int trackId, float pt, float eta, float phi, float mass, int8_t charge = 0, float dcaX = 0.f, float dcaY = 0.f, std::vector<int> amb_muon_self_ids = {})
  {
    fDFId = dfId;
    fGlobalId = globalId;
    fCollisionId = collisionId;
    fTrackId = trackId;
    fPt = pt;
    fEta = eta;
    fPhi = phi;
    fMass = mass;
    fCharge = charge;
    fDCAx = dcaX;
    fDCAy = dcaY;
    fPairDCAXYinSigmaOTF = 0;

    fAmbMuonSelfIds = amb_muon_self_ids;
    if (fAmbMuonSelfIds.size() > 0) {
      fIsAmbiguous = true;
    } else {
      fIsAmbiguous = false;
    }
  }

  ~EMFwdTrack() {}

  int dfId() const { return fDFId; }
  int globalIndex() const { return fGlobalId; }
  int collisionId() const { return fCollisionId; }
  int fwdtrackId() const { return fTrackId; }
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
  bool has_ambiguousMuons() const { return fIsAmbiguous; }
  std::vector<int> ambiguousMuonsIds() const { return fAmbMuonSelfIds; }
  float signed1Pt() const { return fCharge * 1.f / fPt; }

  float pairDcaXYinSigmaOTF() const { return fPairDCAXYinSigmaOTF; }
  void setPairDcaXYinSigmaOTF(float dca) { fPairDCAXYinSigmaOTF = dca; }

 protected:
  int fDFId;
  int fGlobalId;
  int fCollisionId;
  int fTrackId;
  float fPt;
  float fEta;
  float fPhi;
  float fMass;
  int8_t fCharge;
  float fDCAx;
  float fDCAy;
  float fPairDCAXYinSigmaOTF;
  bool fIsAmbiguous;
  std::vector<int> fAmbMuonSelfIds;
};

class EMFwdTrackWithCov : public EMFwdTrack
{
 public:
  EMFwdTrackWithCov(int dfId, int globalId, int collisionId, int trackId, float pt, float eta, float phi, float mass, int8_t charge = 0, float dcaX = 0.f, float dcaY = 0.f, std::vector<int> amb_muon_self_ids = {},
                    float X = 0.f, float Y = 0.f, float Z = 0.f, float Tgl = 0.f,
                    float CXX = 0.f, float CXY = 0.f, float CYY = 0.f,
                    float CPhiX = 0.f, float CPhiY = 0.f, float CPhiPhi = 0.f,
                    float CTglX = 0.f, float CTglY = 0.f, float CTglPhi = 0.f, float CTglTgl = 0.f,
                    float C1PtX = 0.f, float C1PtY = 0.f, float C1PtPhi = 0.f, float C1PtTgl = 0.f, float C1Pt21Pt2 = 0.f, float chi2 = 0.f) : EMFwdTrack(dfId, globalId, collisionId, trackId, pt, eta, phi, mass, charge, dcaX, dcaY, amb_muon_self_ids)
  {
    fX = X;
    fY = Y;
    fZ = Z;
    fTgl = Tgl;
    fCXX = CXX;
    fCXY = CXY;
    fCYY = CYY;
    fCPhiX = CPhiX;
    fCPhiY = CPhiY;
    fCPhiPhi = CPhiPhi;
    fCTglX = CTglX;
    fCTglY = CTglY;
    fCTglPhi = CTglPhi;
    fCTglTgl = CTglTgl;
    fC1PtX = C1PtX;
    fC1PtY = C1PtY;
    fC1PtPhi = C1PtPhi;
    fC1PtTgl = C1PtTgl;
    fC1Pt21Pt2 = C1Pt21Pt2;
    fChi2 = chi2;
  }

  float x() const { return fX; }
  float y() const { return fY; }
  float z() const { return fZ; }
  float tgl() const { return fTgl; }
  float cXX() const { return fCXX; }
  float cXY() const { return fCXY; }
  float cYY() const { return fCYY; }
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

  void setCXX(float cXX) { fCXX = cXX; }
  void setCXY(float cXY) { fCXY = cXY; }
  void setCYY(float cYY) { fCYY = cYY; }

 protected:
  float fX;
  float fY;
  float fZ;
  float fTgl;
  float fCXX;
  float fCXY;
  float fCYY;
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
  float fChi2;
};

} // namespace o2::aod::pwgem::dilepton::utils
#endif // PWGEM_DILEPTON_UTILS_EMFWDTRACK_H_
