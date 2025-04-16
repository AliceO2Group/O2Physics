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
#include "Math/Vector4D.h"

namespace o2::aod::pwgem::dilepton::utils
{
class EMTrack
{
 public:
  EMTrack(int dfId, int globalId, int collisionId, int trackId, float pt, float eta, float phi, float mass, int8_t charge = 0, float dcaXY = 0.f, float dcaZ = 0.f, float cYY = 0.f, float cZY = 0.f, float cZZ = 0.f, std::vector<int> amb_ele_self_ids = {})
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
    fDCAxy = dcaXY;
    fDCAz = dcaZ;
    fCYY = cYY;
    fCZY = cZY;
    fCZZ = cZZ;

    fPairDCA3DinSigmaOTF = 0;

    fAmbEleSelfIds = amb_ele_self_ids;
    if (fAmbEleSelfIds.size() > 0) {
      fIsAmbiguous = true;
    } else {
      fIsAmbiguous = false;
    }
    fVx = 0.f;
    fVy = 0.f;
    fVz = 0.f;
    fVPos = ROOT::Math::PtEtaPhiMVector(0, 0, 0, 0);
    fVNeg = ROOT::Math::PtEtaPhiMVector(0, 0, 0, 0);
    fAmbPosLegSelfIds.clear();
    fAmbNegLegSelfIds.clear();
    fAmbPosLegSelfIds.shrink_to_fit();
    fAmbNegLegSelfIds.shrink_to_fit();
    fGlobalPosId = 0;
    fGlobalNegId = 0;
  }

  ~EMTrack()
  {
    fAmbEleSelfIds.clear();
    fAmbEleSelfIds.shrink_to_fit();
    fAmbPosLegSelfIds.clear();
    fAmbNegLegSelfIds.clear();
    fAmbPosLegSelfIds.shrink_to_fit();
    fAmbNegLegSelfIds.shrink_to_fit();
  }

  int dfId() const { return fDFId; }
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
  float cYY() const { return fCYY; }
  float cZY() const { return fCZY; }
  float cZZ() const { return fCZZ; }
  float p() const { return fPt * std::cosh(fEta); }
  float px() const { return fPt * std::cos(fPhi); }
  float py() const { return fPt * std::sin(fPhi); }
  float pz() const { return fPt * std::sinh(fEta); }
  bool has_ambiguousElectrons() const { return fIsAmbiguous; }
  std::vector<int> ambiguousElectronsIds() const { return fAmbEleSelfIds; }
  float signed1Pt() const { return fCharge * 1.f / fPt; }

  float pairDca3DinSigmaOTF() const { return fPairDCA3DinSigmaOTF; }
  void setPairDca3DinSigmaOTF(float dca) { fPairDCA3DinSigmaOTF = dca; }

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
  float eta_cp() const { return std::atanh(fVz / sqrt(pow(fVx, 2) + pow(fVy, 2) + pow(fVz, 2))); }
  float phi_cp() const { return std::atan2(fVy, fVx); }

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

  ROOT::Math::PtEtaPhiMVector getPositiveLeg() const { return fVPos; }
  ROOT::Math::PtEtaPhiMVector getNegativeLeg() const { return fVNeg; }

  void setGlobalPosId(int id) { fGlobalPosId = id; }
  void setGlobalNegId(int id) { fGlobalNegId = id; }
  int globalIndexPos() const { return fGlobalPosId; }
  int globalIndexNeg() const { return fGlobalNegId; }

  void setAmbPosLegSelfIds(std::vector<int> selfIds) { fAmbPosLegSelfIds = selfIds; }
  void setAmbNegLegSelfIds(std::vector<int> selfIds) { fAmbNegLegSelfIds = selfIds; }
  std::vector<int> ambiguousPosLegIds() const { return fAmbPosLegSelfIds; }
  std::vector<int> ambiguousNegLegIds() const { return fAmbNegLegSelfIds; }

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
  float fDCAxy;
  float fDCAz;
  float fCYY;
  float fCZY;
  float fCZZ;
  float fPairDCA3DinSigmaOTF;
  bool fIsAmbiguous;
  std::vector<int> fAmbEleSelfIds;

  int fGlobalPosId;
  int fGlobalNegId;
  ROOT::Math::PtEtaPhiMVector fVPos;
  ROOT::Math::PtEtaPhiMVector fVNeg;
  std::vector<int> fAmbPosLegSelfIds; // for dileptons
  std::vector<int> fAmbNegLegSelfIds; // for dileptons

  // only for photon conversion point
  float fVx;
  float fVy;
  float fVz;
};

} // namespace o2::aod::pwgem::dilepton::utils
#endif // PWGEM_DILEPTON_UTILS_EMTRACK_H_
