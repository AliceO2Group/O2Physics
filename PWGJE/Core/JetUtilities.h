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

/// \file JetUtilities.h
/// \brief Jet related utilities
///
/// \author Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_CORE_JETUTILITIES_H_
#define PWGJE_CORE_JETUTILITIES_H_

#include "Common/Core/RecoDecay.h"

#include <TKDTree.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <math.h>

namespace jetutilities
{

template <typename T, typename U>
float deltaR(T const& A, U const& B)
{
  float dPhi = RecoDecay::constrainAngle(A.phi() - B.phi(), -M_PI);
  float dEta = A.eta() - B.eta();

  return std::sqrt(dEta * dEta + dPhi * dPhi);
}
// same as deltaR but explicit specification of the eta and phi components
template <typename T, typename U, typename V, typename W>
float deltaR(T const& eta1, U const& phi1, V const& eta2, W const& phi2)
{
  float dPhi = RecoDecay::constrainAngle(phi1 - phi2, -M_PI);
  float dEta = eta1 - eta2;

  return std::sqrt(dEta * dEta + dPhi * dPhi);
}

/// @brief Background estimator using the perpendicular cone method
/// @param inputParticles
/// @param jet
/// @return Rho, RhoM the underlying event density

template <typename T, typename U, typename V>
std::tuple<double, double> estimateRhoPerpCone(const T& inputParticles, const U& jet, V perpConeR)
{

  if (inputParticles.size() == 0) {
    return std::make_tuple(0.0, 0.0);
  }

  double perpPtDensity1 = 0;
  double perpPtDensity2 = 0;
  double perpMdDensity1 = 0;
  double perpMdDensity2 = 0;

  const double jetPhi = RecoDecay::constrainAngle<double, double>(jet.phi(), -M_PI);
  const double jetEta = jet.eta();
  const double radius = static_cast<double>(perpConeR);

  // build 2 perp cones in phi around the leading jet (right and left of the jet)
  double PerpendicularConeAxisPhi1 = RecoDecay::constrainAngle<double, double>(jetPhi + (M_PI / 2.), -M_PI); // This will contrain the angel between -pi & Pi
  double PerpendicularConeAxisPhi2 = RecoDecay::constrainAngle<double, double>(jetPhi - (M_PI / 2.), -M_PI); // This will contrain the angel between -pi & Pi

  for (const auto& particle : inputParticles) {
    // sum the momentum of all paricles that fill the two cones
    const double phi = RecoDecay::constrainAngle<double, double>(particle.phi(), -M_PI);
    double dPhi1 = RecoDecay::constrainAngle<double, double>(phi - PerpendicularConeAxisPhi1, -M_PI); // This will contrain the angel between -pi & Pi
    double dPhi2 = RecoDecay::constrainAngle<double, double>(phi - PerpendicularConeAxisPhi2, -M_PI); // This will contrain the angel between -pi & Pi
    double dEta = jetEta - particle.eta();                                                            // The perp cone eta is the same as the leading jet since the cones are perpendicular only in phi
    if (TMath::Sqrt(dPhi1 * dPhi1 + dEta * dEta) <= static_cast<double>(radius)) {
      perpPtDensity1 += particle.pt();
      perpMdDensity1 += TMath::Sqrt(particle.m() * particle.m() + particle.pt() * particle.pt()) - particle.pt();
    }

    if (TMath::Sqrt(dPhi2 * dPhi2 + dEta * dEta) <= static_cast<double>(radius)) {
      perpPtDensity2 += particle.pt();
      perpMdDensity2 += TMath::Sqrt(particle.m() * particle.m() + particle.pt() * particle.pt()) - particle.pt();
    }
  }

  // Caculate rho as the ratio of average pT of the two cones / the cone area
  double perpPtDensity = (perpPtDensity1 + perpPtDensity2) / (2 * M_PI * static_cast<double>(radius) * static_cast<double>(radius));
  double perpMdDensity = (perpMdDensity1 + perpMdDensity2) / (2 * M_PI * static_cast<double>(radius) * static_cast<double>(radius));

  return std::make_tuple(perpPtDensity, perpMdDensity);
}

}; // namespace jetutilities

#endif // PWGJE_CORE_JETUTILITIES_H_
