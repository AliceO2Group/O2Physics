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

#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"

#include <cmath>

namespace o2::aod::pwgem::photonmeson::utils
{
class EMPhoton
{
 public:
  EMPhoton(float pt, float eta, float phi, float mass) : fPt(pt), fEta(eta), fPhi(phi), fMass(mass)
  {
  }

  ~EMPhoton() = default;

  [[nodiscard]] float pt() const { return fPt; }
  [[nodiscard]] float eta() const { return fEta; }
  [[nodiscard]] float phi() const { return fPhi; }

  [[nodiscard]] float p() const { return fPt * std::cosh(fEta); }
  [[nodiscard]] float px() const { return fPt * std::cos(fPhi); }
  [[nodiscard]] float py() const { return fPt * std::sin(fPhi); }
  [[nodiscard]] float pz() const { return fPt * std::sinh(fEta); }
  [[nodiscard]] float mass() const { return fMass; }
  [[nodiscard]] float rapidity() const { return std::log((std::sqrt(std::pow(fMass, 2) + std::pow(fPt * std::cosh(fEta), 2)) + fPt * std::sinh(fEta)) / std::sqrt(std::pow(fMass, 2) + std::pow(fPt, 2))); }

  // optional leg track-composition info for pair-level photon-class selections (PCM only)
  void setLegCounts(o2::aod::pwgem::photonmeson::utils::pairutil::V0PhotonLegCounts const& c)
  {
    fLegCounts = c;
    fHasLegCounts = true;
  }
  [[nodiscard]] o2::aod::pwgem::photonmeson::utils::pairutil::V0PhotonLegCounts const& legCounts() const { return fLegCounts; }
  [[nodiscard]] bool hasLegCounts() const { return fHasLegCounts; }

 protected:
  float fPt{0.f};
  float fEta{0.f};
  float fPhi{0.f};
  float fMass{0.f};
  o2::aod::pwgem::photonmeson::utils::pairutil::V0PhotonLegCounts fLegCounts{};
  bool fHasLegCounts{false};
};

} // namespace o2::aod::pwgem::photonmeson::utils
#endif // PWGEM_PHOTONMESON_UTILS_EMPHOTON_H_
