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
  EMFwdTrack(int dfId, int globalId, int collisionId, int trackId, float pt, float eta, float phi, float mass, int8_t charge, float dcaX, float dcaY, std::vector<int> amb_muon_self_ids, float cXX, float cXY, float cYY)
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

    fCXX = cXX;
    fCXY = cXY;
    fCYY = cYY;
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

  float cXXatDCA() const { return fCXX; }
  float cXYatDCA() const { return fCXY; }
  float cYYatDCA() const { return fCYY; }

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
  float fCXX;
  float fCXY;
  float fCYY;
};

} // namespace o2::aod::pwgem::dilepton::utils
#endif // PWGEM_DILEPTON_UTILS_EMFWDTRACK_H_
