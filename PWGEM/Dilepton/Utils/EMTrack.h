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

#include <vector>

namespace o2::aod::pwgem::dilepton::utils
{
class EMTrack
{
 public:
  EMTrack(int globalId, int collisionId, int trackId, float pt, float eta, float phi, float mass, int8_t charge = 0, float dcaXY = 0.f, float dcaZ = 0.f, std::vector<int> amb_ele_self_ids = {})
  {
    fGlobalId = globalId;
    fCollisionId = collisionId;
    fTrackId = trackId;
    fPt = pt;
    fEta = eta;
    fPhi = phi;
    fMass = mass;
    fCharge = charge;
    fDCAxy = dcaXY;
    fDCAz = dcaZ;
    fPairDCA3DinSigmaOTF = 0;

    fAmbEleSelfIds = amb_ele_self_ids;
    if (fAmbEleSelfIds.size() > 0) {
      fIsAmbiguous = true;
    } else {
      fIsAmbiguous = false;
    }
  }

  ~EMTrack() {}

  int globalIndex() const { return fGlobalId; }
  int collisionId() const { return fCollisionId; }
  int trackId() const { return fTrackId; }
  float pt() const { return fPt; }
  float eta() const { return fEta; }
  float phi() const { return fPhi; }
  float mass() const { return fMass; }
  int8_t sign() const { return fCharge; }
  float dcaXY() const { return fDCAxy; }
  float dcaZ() const { return fDCAz; }
  float px() const { return fPt * std::cos(fPhi); }
  float py() const { return fPt * std::sin(fPhi); }
  float pz() const { return fPt * std::sinh(fEta); }
  bool has_ambiguousElectrons() const { return fIsAmbiguous; }
  std::vector<int> ambiguousElectronsIds() const { return fAmbEleSelfIds; }
  float signed1Pt() const { return fCharge * 1.f / fPt; }

  float pairDca3DinSigmaOTF() const { return fPairDCA3DinSigmaOTF; }
  void setPairDca3DinSigmaOTF(float dca) { fPairDCA3DinSigmaOTF = dca; }

 protected:
  int fGlobalId;
  int fCollisionId;
  int fTrackId;
  float fPt;
  float fEta;
  float fPhi;
  float fMass;
  int8_t fCharge;
  float fDCAxy;
  float fDCAz;
  float fPairDCA3DinSigmaOTF;
  bool fIsAmbiguous;
  std::vector<int> fAmbEleSelfIds;
};

class EMTrackWithCov : public EMTrack
{
 public:
  EMTrackWithCov(int globalId, int collisionId, int trackId, float pt, float eta, float phi, float mass, int8_t charge = 0, float dcaXY = 0.f, float dcaZ = 0.f, std::vector<int> amb_ele_self_ids = {},
                 float X = 0.f, float Y = 0.f, float Z = 0.f, float Alpha = 0.f, float Snp = 0.f, float Tgl = 0.f,
                 float CYY = 0.f, float CZY = 0.f, float CZZ = 0.f,
                 float CSnpY = 0.f, float CSnpZ = 0.f, float CSnpSnp = 0.f,
                 float CTglY = 0.f, float CTglZ = 0.f, float CTglSnp = 0.f, float CTglTgl = 0.f,
                 float C1PtY = 0.f, float C1PtZ = 0.f, float C1PtSnp = 0.f, float C1PtTgl = 0.f, float C1Pt21Pt2 = 0.f) : EMTrack(globalId, collisionId, trackId, pt, eta, phi, mass, charge, dcaXY, dcaZ, amb_ele_self_ids)
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

  float x() const { return fX; }
  float y() const { return fY; }
  float z() const { return fZ; }
  float alpha() const { return fAlpha; }
  float snp() const { return fSnp; }
  float tgl() const { return fTgl; }

  float cYY() const { return fCYY; }
  float cZY() const { return fCZY; }
  float cZZ() const { return fCZZ; }
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

  void setCYY(float cYY) { fCYY = cYY; }
  void setCZY(float cZY) { fCZY = cZY; }
  void setCZZ(float cZZ) { fCZZ = cZZ; }

 protected:
  float fX;
  float fY;
  float fZ;
  float fAlpha;
  float fSnp;
  float fTgl;
  float fCYY;
  float fCZY;
  float fCZZ;
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

} // namespace o2::aod::pwgem::dilepton::utils
#endif // PWGEM_DILEPTON_UTILS_EMTRACK_H_
