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

/// \class to store minimal photon info
/// \author daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_EMPHOTON_H_
#define PWGEM_PHOTONMESON_UTILS_EMPHOTON_H_

#include <cmath>

#include <math.h>

namespace o2::aod::pwgem::photonmeson::utils
{
class EMPhoton
{
 public:
  EMPhoton(float pt, float eta, float phi, float mass)
  {
    fPt = pt;
    fEta = eta;
    fPhi = phi;
    fMass = mass;
  }

  ~EMPhoton() {}

  float pt() const { return fPt; }
  float eta() const { return fEta; }
  float phi() const { return fPhi; }

  float p() const { return fPt * std::cosh(fEta); }
  float px() const { return fPt * std::cos(fPhi); }
  float py() const { return fPt * std::sin(fPhi); }
  float pz() const { return fPt * std::sinh(fEta); }
  float mass() const { return fMass; }
  float rapidity() const { return std::log((std::sqrt(std::pow(fMass, 2) + std::pow(fPt * std::cosh(fEta), 2)) + fPt * std::sinh(fEta)) / std::sqrt(std::pow(fMass, 2) + std::pow(fPt, 2))); }

 protected:
  float fPt;
  float fEta;
  float fPhi;
  float fMass;
};

} // namespace o2::aod::pwgem::photonmeson::utils
#endif // PWGEM_PHOTONMESON_UTILS_EMPHOTON_H_
