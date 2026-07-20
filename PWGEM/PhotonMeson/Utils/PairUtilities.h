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

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>

#include <cmath>

namespace o2::aod::pwgem::photonmeson::utils::pairutil
{
enum class PhotonPrefilterBitDerived : int {
  kPhotonFromPi0gg = 0,  // photon from pi0->gg
  kPhotonFromPi0eeg = 1, // photon from pi0->eeg
};
enum class ElectronPrefilterBitDerived : int {
  kElectronFromPi0eeg = 0, // electron from pi0->eeg
  kElectronFromFakePC = 1, // electron from photon->ee, misidentified photon conversion as virtual photon
};
} // namespace o2::aod::pwgem::photonmeson::utils::pairutil
namespace o2::aod::pwgem::photonmeson::photonpair
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
  kEEEE = 8, // dielectron-dielectron
  kNpair,
};

template <typename U1, typename U2, o2::soa::is_iterator TG1, o2::soa::is_iterator TG2, typename TCut1, typename TCut2>
bool IsSelectedPair(TG1 const& g1, TG2 const& g2, TCut1 const& cut1, TCut2 const& cut2)
{
  bool is_g1_selected = false;
  bool is_g2_selected = false;
  is_g1_selected = cut1.template IsSelected<TG1, U1>(g1);
  is_g2_selected = cut2.template IsSelected<TG2, U2>(g2);
  return (is_g1_selected && is_g2_selected);
}

template <o2::soa::is_iterator TV0Leg, o2::soa::is_iterator TCluster>
bool DoesV0LegMatchWithCluster(TV0Leg const& v0leg, TCluster const& cluster, const float max_deta, const float max_dphi, const float max_Ep_width)
{
  float deta = v0leg.eta() - cluster.eta();
  float dphi = RecoDecay::constrainAngle(RecoDecay::constrainAngle(v0leg.phi()) - RecoDecay::constrainAngle(cluster.phi()), -o2::constants::math::PI);
  // float dR = sqrt(deta * deta + dphi * dphi);
  float Ep = cluster.e() / v0leg.p();
  return (std::pow(deta / max_deta, 2.f) + std::pow(dphi / max_dphi, 2.f) < 1.f) && (std::abs(Ep - 1.f) < max_Ep_width);
}
} // namespace o2::aod::pwgem::photonmeson::photonpair

namespace o2::aod::pwgem::photonmeson::utils::pairutil
{
// ─── photon-class selection by leg track composition ──────────────────────

struct V0PhotonLegCounts {
  int nITSTPC{0};
  int nITSOnly{0};
  int nTPCOnly{0};
  int nTRD{0};
  int nTOF{0};
};

template <o2::soa::is_iterator TLeg>
V0PhotonLegCounts getV0PhotonLegCounts(TLeg const& pos, TLeg const& ele)
{
  V0PhotonLegCounts c;
  auto countLeg = [&c](auto const& l) {
    const bool its = l.hasITS();
    const bool tpc = l.hasTPC();
    if (its && tpc) {
      c.nITSTPC++;
    } else if (its) {
      c.nITSOnly++;
    } else {
      c.nTPCOnly++;
    }
    if (l.hasTRD()) {
      c.nTRD++;
    }
    if (l.hasTOF()) {
      c.nTOF++;
    }
  };
  countLeg(pos);
  countLeg(ele);
  return c;
}

struct V0PhotonClassSelection {
  int minITSTPC{0}, maxITSTPC{2};
  int minITSOnly{0}, maxITSOnly{2};
  int minTPCOnly{0}, maxTPCOnly{2};
  int minTRD{0}, maxTRD{2};
  int minTOF{0}, maxTOF{2};

  [[nodiscard]] bool isSelected(V0PhotonLegCounts const& c) const
  {
    return c.nITSTPC >= minITSTPC && c.nITSTPC <= maxITSTPC &&
           c.nITSOnly >= minITSOnly && c.nITSOnly <= maxITSOnly &&
           c.nTPCOnly >= minTPCOnly && c.nTPCOnly <= maxTPCOnly &&
           c.nTRD >= minTRD && c.nTRD <= maxTRD &&
           c.nTOF >= minTOF && c.nTOF <= maxTOF;
  }
};

inline bool isPairPhotonClassSelected(V0PhotonLegCounts const& c1, V0PhotonLegCounts const& c2,
                                      V0PhotonClassSelection const& selA, V0PhotonClassSelection const& selB)
{
  return (selA.isSelected(c1) && selB.isSelected(c2)) ||
         (selA.isSelected(c2) && selB.isSelected(c1));
}

template <typename TGroup>
V0PhotonClassSelection buildV0PhotonClassSelection(TGroup const& g)
{
  V0PhotonClassSelection s;
  s.minITSTPC = g.cfgMinNLegsITSTPC.value;
  s.maxITSTPC = g.cfgMaxNLegsITSTPC.value;
  s.minITSOnly = g.cfgMinNLegsITSOnly.value;
  s.maxITSOnly = g.cfgMaxNLegsITSOnly.value;
  s.minTPCOnly = g.cfgMinNLegsTPCOnly.value;
  s.maxTPCOnly = g.cfgMaxNLegsTPCOnly.value;
  s.minTRD = g.cfgMinNLegsTRD.value;
  s.maxTRD = g.cfgMaxNLegsTRD.value;
  s.minTOF = g.cfgMinNLegsTOF.value;
  s.maxTOF = g.cfgMaxNLegsTOF.value;
  return s;
}
} // namespace o2::aod::pwgem::photonmeson::utils::pairutil

#endif // PWGEM_PHOTONMESON_UTILS_PAIRUTILITIES_H_
