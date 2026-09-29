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
#include <Framework/Concepts.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h>
#include <Math/Vector4Dfwd.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

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
  kPHOSPHOS,
  kEMCEMC,
  kPCMPHOS,
  kPCMEMC,
  kPCMDalitzEE,
  kPCMDalitzMuMu,
  kPHOSEMC,
  kEEEE, // dielectron-dielectron
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

//_______________________________________________________________________
/// \brief Pair observables of a photon pair: kT, q_inv and the LCMS decomposition
struct PairQ {
  float kt{0.f};
  float qinv{0.f};
  float qabsLcms{0.f};
  float qout{0.f}, qside{0.f}, qlong{0.f};
};

//_______________________________________________________________________
/// \brief Cheap q_inv gate before any boost or sparse fill
/// \param pt1 pT of the first photon
/// \param eta1 eta of the first photon
/// \param phi1 phi of the first photon
/// \param pt2 pT of the second photon
/// \param eta2 eta of the second photon
/// \param phi2 phi of the second photon
/// \param qMax upper limit on q_inv
/// \return true if q_inv <= qMax, using q_inv^2 = 2 pT1 pT2 (cosh(deta) - cos(dphi)), exact for massless particles
[[nodiscard]] inline bool pairBelowQmax(float pt1, float eta1, float phi1,
                                        float pt2, float eta2, float phi2, float qMax)
{
  const float q2 = 2.f * pt1 * pt2 * (std::cosh(eta1 - eta2) - std::cos(phi1 - phi2));
  return q2 <= qMax * qMax;
}

//_______________________________________________________________________
/// \brief Computes kT, q_inv and (q_out, q_side, q_long) in the LCMS
/// \param v1 four-vector of the first photon
/// \param v2 four-vector of the second photon
/// \return PairQ;
inline PairQ computePairQ(ROOT::Math::PtEtaPhiMVector const& v1, ROOT::Math::PtEtaPhiMVector const& v2)
{
  PairQ r;
  const auto k12 = 0.5 * (v1 + v2);
  r.kt = static_cast<float>(k12.Pt());
  r.qinv = std::fabs(static_cast<float>((v1 - v2).M()));
  const double ktMag = k12.Pt();
  const ROOT::Math::XYZVector uvOut = (ktMag > 1e-9) // o2-linter: disable=magic-number (guard against a zero-pT pair)
                                        ? ROOT::Math::XYZVector(k12.Px() / ktMag, k12.Py() / ktMag, 0.)
                                        : ROOT::Math::XYZVector(1., 0., 0.);
  const ROOT::Math::XYZVector uvLong(0., 0., 1.);
  const ROOT::Math::XYZVector uvSide = uvOut.Cross(uvLong);
  const ROOT::Math::PxPyPzEVector v1c(v1), v2c(v2);
  const double betaZ = (v1 + v2).Beta() * std::cos((v1 + v2).Theta());
  const ROOT::Math::Boost bstZ(0., 0., -betaZ);
  const auto q3Lcms = bstZ(v1c - v2c).Vect();
  r.qabsLcms = static_cast<float>(q3Lcms.R());
  r.qout = std::fabs(static_cast<float>(q3Lcms.Dot(uvOut)));
  r.qside = std::fabs(static_cast<float>(q3Lcms.Dot(uvSide)));
  r.qlong = std::fabs(static_cast<float>(q3Lcms.Dot(uvLong)));
  return r;
}

//_______________________________________________________________________
/// \brief computePairQ for anything with pt(), eta() and phi(): V0 photons, MC photons, truth records
/// \param a first photon-like object
/// \param b second photon-like object
/// \return PairQ of the massless pair
template <typename TA, typename TB>
PairQ computePairQ(TA const& a, TB const& b)
{
  return computePairQ(ROOT::Math::PtEtaPhiMVector(a.pt(), a.eta(), a.phi(), 0.f),
                      ROOT::Math::PtEtaPhiMVector(b.pt(), b.eta(), b.phi(), 0.f));
}

//_______________________________________________________________________
/// \brief One leg of a V0 candidat
struct DedupLeg {
  float pt{0.f}, eta{0.f}, phi{0.f};
  float dedx{0.f};       /// TPC dE/dx
  float fracShared{0.f}; /// fraction of shared TPC clusters
  bool hasTpc{false};
  int64_t trackId{-1};
};

//_______________________________________________________________________
/// \brief A selected V0 candidate
struct DedupCand {
  int64_t gi{-1}; /// globalIndex of the V0
  float eta{0.f}, phi{0.f}, pt{0.f}, rConv{0.f}, chi2{1e9f};
  float vx{0.f}, vy{0.f}, vz{0.f};
  int nITSTPC{0};
  std::array<DedupLeg, 2> leg{}; /// 0 = e+, 1 = e-
};

//_______________________________________________________________________
/// \brief Tolerances of the duplicate test
struct DedupConfig {
  float maxLegDEta{0.01f};
  float maxLegDPhi{0.01f};
  float maxLegPtAsym{0.05f};
  float maxLegDeDxAsym{0.05f};
  float maxDVtx3D{999.f};
  bool requireBothLegs{false};
};

//_______________________________________________________________________
/// \brief Distance of two same-charge legs of two candidates
struct LegDistance {
  float dEta{0.f}, dPhi{0.f}, ptAsym{1.f}, dedxAsym{1.f}, fShared{0.f};
  bool dedxValid{false};
};

//_______________________________________________________________________
/// \brief Fills DedupLeg from a V0 leg
/// \param leg V0 leg iterator (aod::V0Legs)
/// \return the leg reduced to the duplicate-test variables
template <o2::soa::is_iterator TLeg>
DedupLeg makeDedupLeg(TLeg const& leg)
{
  DedupLeg l;
  l.pt = static_cast<float>(leg.pt());
  l.eta = static_cast<float>(leg.eta());
  l.phi = static_cast<float>(leg.phi());
  l.hasTpc = leg.hasTPC();
  l.dedx = static_cast<float>(leg.tpcSignal());
  l.fracShared = static_cast<float>(leg.tpcFractionSharedCls());
  l.trackId = leg.trackId();
  return l;
}

//_______________________________________________________________________
/// \brief Fills DedupCand from a V0 photon and its two legs
/// \param v0 V0 photon iterator
/// \param pos positive leg
/// \param ele negative leg
/// \return the candidate reduced to what defines its identity
template <o2::soa::is_iterator TV0, o2::soa::is_iterator TLeg>
DedupCand makeDedupCand(TV0 const& v0, TLeg const& pos, TLeg const& ele)
{
  DedupCand c;
  c.gi = v0.globalIndex();
  c.eta = v0.eta();
  c.phi = v0.phi();
  c.pt = v0.pt();
  c.rConv = std::hypot(v0.vx(), v0.vy());
  c.chi2 = v0.chiSquareNDF();
  c.vx = v0.vx();
  c.vy = v0.vy();
  c.vz = v0.vz();
  c.nITSTPC = getV0PhotonLegCounts(pos, ele).nITSTPC;
  c.leg = {makeDedupLeg(pos), makeDedupLeg(ele)};
  return c;
}

//_______________________________________________________________________
/// \brief Distance of two same-charge legs
/// \param a first candidate
/// \param b second candidate
/// \param il 0 = e+ with e+, 1 = e- with e-
/// \return LegDistance; dedxValid is false when a leg has no TPC
[[nodiscard]] inline LegDistance legDistance(DedupCand const& a, DedupCand const& b, int il)
{
  constexpr float kMinSigma = 1e-9f;
  const auto& la = a.leg[static_cast<size_t>(il)];
  const auto& lb = b.leg[static_cast<size_t>(il)];
  LegDistance d;
  d.dEta = std::fabs(la.eta - lb.eta);
  d.dPhi = std::fabs(RecoDecay::constrainAngle(la.phi - lb.phi, -o2::constants::math::PI));
  const float ptSum = la.pt + lb.pt;
  d.ptAsym = (ptSum > kMinSigma) ? std::fabs(la.pt - lb.pt) / ptSum : 1.f;
  d.dedxValid = la.hasTpc && lb.hasTpc;
  const float dedxSum = la.dedx + lb.dedx;
  d.dedxAsym = (d.dedxValid && dedxSum > kMinSigma) ? std::fabs(la.dedx - lb.dedx) / dedxSum : 1.f;
  d.fShared = std::max(la.fracShared, lb.fracShared);
  return d;
}

//_______________________________________________________________________
/// \brief Test of track identity on reco level: same-charge legs physically the same particle
/// \param a first candidate
/// \param b second candidate
/// \param il 0 = e+, 1 = e-
/// \param cfg tolerances
/// \return true if the two legs agree within all tolerances
[[nodiscard]] inline bool legsIdentical(DedupCand const& a, DedupCand const& b, int il, DedupConfig const& cfg)
{
  if (a.leg[static_cast<size_t>(il)].trackId == b.leg[static_cast<size_t>(il)].trackId) {
    return false;
  }
  const auto d = legDistance(a, b, il);
  if (d.dEta > cfg.maxLegDEta || d.dPhi > cfg.maxLegDPhi) {
    return false;
  }
  if (d.ptAsym > cfg.maxLegPtAsym) {
    return false;
  }
  if (d.dedxValid && d.dedxAsym > cfg.maxLegDeDxAsym) {
    return false;
  }
  return true;
}

//_______________________________________________________________________
/// \brief Test: Duplicate stores
/// \param a first candidate
/// \param b second candidate
/// \param cfg tolerances
/// \return true if at least one same-charge leg pair is identical (both, with cfg.requireBothLegs) and the conversion points are within cfg.maxDVtx3D
[[nodiscard]] inline bool isDuplicatePair(DedupCand const& a, DedupCand const& b, DedupConfig const& cfg)
{
  const float dVtx = std::hypot(a.vx - b.vx, a.vy - b.vy, a.vz - b.vz);
  if (dVtx > cfg.maxDVtx3D) {
    return false;
  }
  const bool posSame = legsIdentical(a, b, 0, cfg);
  const bool negSame = legsIdentical(a, b, 1, cfg);
  return cfg.requireBothLegs ? (posSame && negSame) : (posSame || negSame);
}

//_______________________________________________________________________
/// \brief Quality order among duplicates
/// \param a first candidate
/// \param b second candidate
/// \return true if a should be kept over b
[[nodiscard]] inline bool isBetterCand(DedupCand const& a, DedupCand const& b)
{
  if (a.nITSTPC != b.nITSTPC) {
    return a.nITSTPC > b.nITSTPC;
  }
  return a.chi2 < b.chi2;
}

} // namespace o2::aod::pwgem::photonmeson::utils::pairutil

#endif // PWGEM_PHOTONMESON_UTILS_PAIRUTILITIES_H_
