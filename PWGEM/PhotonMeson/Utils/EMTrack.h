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

#ifndef PWGEM_PHOTONMESON_UTILS_EMTRACK_H_
#define PWGEM_PHOTONMESON_UTILS_EMTRACK_H_

#include <vector>

class EMTrack
{
 public:
  EMTrack(int collisionId, int trackId, float pt, float eta, float phi, float mass, int8_t charge, float dcaXY, float dcaZ, float cYY, float cZZ, float cZY, std::vector<int> amb_ele_self_ids)
  {
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
    fCZZ = cZZ;
    fCZY = cZY;

    fAmbEleSelfIds = amb_ele_self_ids;
    if (fAmbEleSelfIds.size() > 0) {
      fIsAmbiguous = true;
    } else {
      fIsAmbiguous = false;
    }
  }

  ~EMTrack() {}

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
  float cZZ() const { return fCZZ; }
  float cZY() const { return fCZY; }
  float px() const { return fPt * std::cos(fPhi); }
  float py() const { return fPt * std::sin(fPhi); }
  float pz() const { return fPt * std::sinh(fEta); }
  bool has_ambiguousElectrons() const { return fIsAmbiguous; }
  std::vector<int> ambiguousElectronsIds() const { return fAmbEleSelfIds; }

 private:
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
  float fCZZ;
  float fCZY;
  bool fIsAmbiguous;
  std::vector<int> fAmbEleSelfIds;
};

#endif // PWGEM_PHOTONMESON_UTILS_EMTRACK_H_
