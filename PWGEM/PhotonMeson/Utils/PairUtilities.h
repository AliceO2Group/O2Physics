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

/// \commonly used for pair analyses.
/// \author daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_PAIRUTILITIES_H_
#define PWGEM_PHOTONMESON_UTILS_PAIRUTILITIES_H_

#include <TVector2.h>
#include <cmath>

namespace o2::aod::photonpair
{
enum PairType {
  kPCMPCM = 0,
  kPHOSPHOS = 1,
  kEMCEMC = 2,
  kPCMPHOS = 3,
  kPCMEMC = 4,
  kPCMDalitzEE = 5,
  kPCMDalitzMuMu = 6,
  kPHOSEMC = 7,
  kDalitzEEDalitzEE = 8,
  kNpair,
};

template <typename U1, typename U2, typename TG1, typename TG2, typename TCut1, typename TCut2>
bool IsSelectedPair(TG1 const& g1, TG2 const& g2, TCut1 const& cut1, TCut2 const& cut2)
{
  bool is_g1_selected = false;
  bool is_g2_selected = false;
  is_g1_selected = cut1.template IsSelected<U1>(g1);
  is_g2_selected = cut2.template IsSelected<U2>(g2);
  return (is_g1_selected && is_g2_selected);
}

template <typename TV0Leg, typename TCluster>
bool DoesV0LegMatchWithCluster(TV0Leg const& v0leg, TCluster const& cluster, const float max_deta, const float max_dphi, const float max_Ep_width)
{
  float deta = v0leg.eta() - cluster.eta();
  float dphi = TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(v0leg.phi()) - TVector2::Phi_0_2pi(cluster.phi()));
  // float dR = sqrt(deta * deta + dphi * dphi);
  float Ep = cluster.e() / v0leg.p();
  return (pow(deta / max_deta, 2) + pow(dphi / max_dphi, 2) < 1) && (abs(Ep - 1) < max_Ep_width);
}
} // namespace o2::aod::photonpair

#endif // PWGEM_PHOTONMESON_UTILS_PAIRUTILITIES_H_
