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
///
/// \file PhotonHbtUtilities.h
/// \brief Helpers of the V0 photon HBT task (photonhbt.cxx): pair kinematics, MC identity of photon pairs, leg geometry.
///        Free functions of photons and pairs; nothing here depends on the task.
/// \author Stefanie Mrozinski, stefanie.mrozinski@cern.ch
///

#ifndef PWGEM_PHOTONMESON_UTILS_PHOTONHBTUTILITIES_H_
#define PWGEM_PHOTONMESON_UTILS_PHOTONHBTUTILITIES_H_

#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Utils/V0DedupUtilities.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <MathUtils/Primitive2D.h>
#include <MathUtils/Utils.h>
#include <ReconstructionDataFormats/Track.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector3D.h> // IWYU pragma: keep
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h> // IWYU pragma: keep
#include <Math/Vector4Dfwd.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <iterator>
#include <tuple>
#include <utility>
#include <vector>

namespace o2::aod::pwgem::photonmeson::hbt
{

inline constexpr float kMinMagnitude = 1e-12f;
inline constexpr float kMinCosine = 1e-12f;
inline constexpr float kMinSigma = 1e-9f;

struct PairQAObservables {
  ROOT::Math::PtEtaPhiMVector v1, v2, k12;
  float x1 = 0.f, y1 = 0.f, z1 = 0.f, x2 = 0.f, y2 = 0.f, z2 = 0.f;
  float r1 = 0.f, r2 = 0.f, dx = 0.f, dy = 0.f, dz = 0.f;
  float deltaR = 0.f, deltaZ = 0.f, deltaR3D = 0.f;
  float opa = 0.f, cosOA = 0.f, drOverCosOA = 0.f;
  float deta = 0.f, dphi = 0.f;
  float kt = 0.f, qinv = 0.f;
  bool valid = true;
};

template <typename TG1, typename TG2>
PairQAObservables buildPairQAObservables(TG1 const& g1, TG2 const& g2)
{
  PairQAObservables o{};
  o.x1 = g1.vx();
  o.y1 = g1.vy();
  o.z1 = g1.vz();
  o.x2 = g2.vx();
  o.y2 = g2.vy();
  o.z2 = g2.vz();
  o.r1 = std::sqrt(o.x1 * o.x1 + o.y1 * o.y1);
  o.r2 = std::sqrt(o.x2 * o.x2 + o.y2 * o.y2);
  o.dx = o.x1 - o.x2;
  o.dy = o.y1 - o.y2;
  o.dz = o.z1 - o.z2;
  o.deltaR = std::fabs(o.r1 - o.r2);
  o.deltaZ = o.dz;
  o.deltaR3D = std::sqrt(o.dx * o.dx + o.dy * o.dy + o.dz * o.dz);
  ROOT::Math::XYZVector cp1(o.x1, o.y1, o.z1), cp2(o.x2, o.y2, o.z2);
  const float mag1 = std::sqrt(cp1.Mag2()), mag2 = std::sqrt(cp2.Mag2());
  if (mag1 < kMinMagnitude || mag2 < kMinMagnitude) {
    o.valid = false;
    return o;
  }
  auto cosPA = static_cast<float>(cp1.Dot(cp2) / (mag1 * mag2));
  cosPA = std::clamp(cosPA, -1.f, 1.f);
  o.opa = std::acos(cosPA);
  o2::math_utils::bringTo02Pi(o.opa);
  if (o.opa > o2::constants::math::PI) {
    o.opa = o2::constants::math::TwoPI - o.opa;
  }
  o.cosOA = std::cos(o.opa / 2.f);
  o.drOverCosOA = (std::fabs(o.cosOA) < kMinCosine) ? 1e12f : (o.deltaR3D / o.cosOA);
  o.v1 = ROOT::Math::PtEtaPhiMVector(g1.pt(), g1.eta(), g1.phi(), 0.f);
  o.v2 = ROOT::Math::PtEtaPhiMVector(g2.pt(), g2.eta(), g2.phi(), 0.f);
  o.k12 = 0.5f * (o.v1 + o.v2);
  o.deta = g1.eta() - g2.eta();
  o.dphi = RecoDecay::constrainAngle(g1.phi() - g2.phi(), -o2::constants::math::PI);
  o.kt = o.k12.Pt();
  o.qinv = std::fabs((o.v1 - o.v2).M());
  return o;
}

inline std::array<float, 3> lcmsQ(ROOT::Math::PtEtaPhiMVector const& v1, ROOT::Math::PtEtaPhiMVector const& v2)
{
  const auto k12 = 0.5 * (v1 + v2);
  const ROOT::Math::XYZVector uvOut(k12.Px() / k12.Pt(), k12.Py() / k12.Pt(), 0), uvLong(0, 0, 1);
  const ROOT::Math::XYZVector uvSide = uvOut.Cross(uvLong);
  const ROOT::Math::PxPyPzEVector v1c(v1), v2c(v2);
  const float betaZ = (v1 + v2).Beta() * std::cos((v1 + v2).Theta());
  const ROOT::Math::Boost bstZ(0, 0, -betaZ);
  const auto q3 = bstZ(v1c - v2c).Vect();
  return {static_cast<float>(std::fabs(q3.Dot(uvOut))), static_cast<float>(std::fabs(q3.Dot(uvSide))), static_cast<float>(std::fabs(q3.Dot(uvLong)))};
}

inline constexpr int kPcGamma = 22, kPcElectron = 11; ///< PDG codes compared with the MC columns of PCMMixLegs

enum PairClass : int {
  kClsSamePhoton = 1,            ///< the same MC photon twice (duplicate)
  kClsSharedLegPhotons = 2,      ///< one MC particle under two track ids, the other two legs both from photons (duplicate)
  kClsSharedLeg = 3,             ///< one MC particle under two track ids, otherwise (duplicate)
  kClsHalfSwap = 4,              ///< the two photons share one true photon mother through one leg each
  kClsBremsFake = 5,             ///< related by MC history (bremsstrahlung, leg daughter of a leg), with a fake photon
  kClsSharedNonPhotonMother = 6, ///< two legs with the same non-photon mother
  kClsBothFake = 7,              ///< both photons fake, otherwise unrelated
  kClsFakeNonPhotonMother = 8,   ///< one fake whose legs share a non-photon mother
  kClsFakeNonElectron = 9,       ///< one fake with a non-electron leg
  kClsFakeTwoPhotons = 10,       ///< one fake built from legs of two different photons
  kClsFakeOther = 11,            ///< one fake, otherwise
  kClsSharedTrack = 12,          ///< the two photons share an AO2D track (duplicate)
  kClsTrueDistinct = 13,         ///< two different true photons
  kClsFullSwap = 14,             ///< the legs of two true photons crossed
  kClsNoLabel = 15,              ///< a leg without MC label
  kClsBrems = 16,                ///< two true photons, one radiated by a leg of the other (bremsstrahlung)
  kClsRelatedSwap = 17,          ///< full swap of a bremsstrahlung pair
  kClsSameMother = 18,           ///< two true photons of the same mother (pi0, eta -> gamma gamma)
  kClsSameMotherSwap = 19,       ///< full swap of two such photons
  kClsSameMotherFake = 20        ///< photons of the same mother, with a fake
};

/// coarse MC truth of a photon pair
enum PairTruth : int {
  kTruthTwoPhotons = 1, ///< two true photons
  kTruthFullSwap = 2,   ///< the legs of two true photons crossed
  kTruthWithFake = 3,   ///< anything with a fake photon
  kTruthNoLabel = 4     ///< a leg without MC mother
};

/// relation of the two photons of a pair through their legs
enum PartnerRelation : int {
  kPartnerOriginals = 0,    ///< both photons are run-1 pairings: the rule does not apply
  kPartnerSharedSource = 1, ///< the photons share a run-1 source photon through their legs
  kPartnerFreeLegs = 2,     ///< both hold a pooled free leg of the same source event, with opposite charge
  kPartnerUnrelated = 3     ///< none of these
};

struct PcLeg {
  int64_t evt = -1;
  int pid = -1, mid = -1, gmid = -1, pdg = 0, mpdg = 0, gpdg = 0, trk = -1, anc = -1, ggid = -1, root = -1, ggpdg = 0, rootPdg = 0;
  float sharedFrac = 0.f;
  float nCls = 0.f;                      // TPC clusters found (PCMMixLegs; analysis-note QA)
  float vx = 0.f, vy = 0.f, vz = 0.f;    // MC production vertex = true conversion point (valid for pid >= 0)
  float px = 0.f, py = 0.f, pz = 0.f;    // reconstructed leg momentum at the fitted vertex, before the mass constraint (V0LegsMix)
  float mpx = 0.f, mpy = 0.f, mpz = 0.f; // MC momentum of the mother (the true photon for a conversion leg)
  bool its = false;                      // leg has ITS hits
  bool free = false;                     // pooled free leg (PCMMixLegs sourceV0Id -2, overlay mixed event)
  int64_t srcPhoton = -1;                // run-1 photon the leg came from (PCMMixLegs sourcePhotonUid; own uid for a free leg)
};

struct PcPhoton {
  float px = 0.f, py = 0.f, pz = 0.f;
  float vx = 0.f, vy = 0.f, vz = 0.f;                     // reconstructed conversion point
  bool original = true;                                   // run-1 pairing (false: crossed pairing built in run 2)
  bool svertexer = true;                                  // run-1 photon was an SVertexer V0 (false: built by the finder)
  float chi2 = 0.f, dcaxy = 0.f, dcaz = 0.f, cospa = 1.f; // KF chi2/ndf, DCA to the PV (cm), cos of the pointing angle to the PV
  float chi2Vtx = -1.f, chi2PV = -1.f, psipair = 0.f;     // KF chi2/ndf of the vertex fit alone, chi2 to the PV, psi_pair
  int nFree = 0;                                          // legs that were pooled free legs (PCMMixLegs sourceV0Id -2, overlay): 0, 1, 2
  bool pass = false;                                      // passes the analysis photon cut (set by the task where it is read)
  int64_t posMix = -1, eleMix = -1;                       // PCMMixLegs rows of the legs (track parameters for the bremsstrahlung pointing tag)
  PcLeg pos, ele;
};

/// The two legs are the same MC particle of the same source event.
inline bool pcSamePart(PcLeg const& x, PcLeg const& y) { return x.evt == y.evt && x.pid >= 0 && x.pid == y.pid; }
/// The two legs have the same MC mother.
inline bool pcSameMother(PcLeg const& x, PcLeg const& y) { return x.evt == y.evt && x.mid >= 0 && x.mid == y.mid; }
/// Leg x is a daughter of leg y.
inline bool pcDaughterOf(PcLeg const& x, PcLeg const& y) { return x.evt == y.evt && x.mid >= 0 && x.mid == y.pid; }
/// The photon of leg x was radiated by leg y (x's grandmother is y).
inline bool pcRadiatedBy(PcLeg const& x, PcLeg const& y) { return x.evt == y.evt && x.gmid >= 0 && x.gmid == y.pid; }
/// The mothers of the two legs are different daughters of the same particle.
inline bool pcSiblings(PcLeg const& x, PcLeg const& y) { return x.evt == y.evt && x.gmid >= 0 && x.gmid == y.gmid && x.mid != y.mid; }
/// Both legs come from the same true photon.
inline bool pcTruePhoton(PcPhoton const& p) { return p.pos.mpdg == kPcGamma && p.ele.mpdg == kPcGamma && pcSameMother(p.pos, p.ele); }

/// Bremsstrahlung: a leg of one photon radiated the other photon (the photon's grandmother is that electron).
inline bool pcRelated(PcPhoton const& a, PcPhoton const& b)
{
  for (const auto* x : {&a.pos, &a.ele}) {
    for (const auto* y : {&b.pos, &b.ele}) {
      if (pcRadiatedBy(*x, *y) || pcRadiatedBy(*y, *x)) {
        return true;
      }
    }
  }
  return false;
}

/// Photons of the same mother (pi0 / eta -> gamma gamma): the legs share the grandmother but not the mother.
inline bool pcSiblingPair(PcPhoton const& a, PcPhoton const& b)
{
  for (const auto* x : {&a.pos, &a.ele}) {
    for (const auto* y : {&b.pos, &b.ele}) {
      if (pcSiblings(*x, *y)) {
        return true;
      }
    }
  }
  return false;
}

/// Truth of a pair (PairTruth): 1 two true photons, 2 the legs of two true photons crossed (full swap),
/// 3 other with a fake, 4 no label.
inline int pcPairTruth(PcPhoton const& a, PcPhoton const& b)
{
  if (a.pos.mid < 0 || a.ele.mid < 0 || b.pos.mid < 0 || b.ele.mid < 0) {
    return kTruthNoLabel;
  }
  auto key = [](PcLeg const& l) { return std::make_pair(l.evt, l.mid); };
  auto isTrue = [](PcLeg const& l) { return l.mid >= 0 && l.mpdg == kPcGamma; };
  const bool allTrue = isTrue(a.pos) && isTrue(a.ele) && isTrue(b.pos) && isTrue(b.ele);
  if (allTrue && key(a.pos) == key(a.ele) && key(b.pos) == key(b.ele) && key(a.pos) != key(b.pos)) {
    return kTruthTwoPhotons;
  }
  if (allTrue && key(a.pos) == key(b.ele) && key(a.ele) == key(b.pos) && key(a.pos) != key(a.ele)) {
    return kTruthFullSwap;
  }
  return kTruthWithFake;
}

/// Class of a pair with a fake (PairClass): 12 shared AO2D track, 1 same MC photon twice, 2/3 shared MC particle under two
/// track ids, 4 shared photon mother (half swap), 5 related by MC history (brems, leg daughter of a leg), 20 photons of the
/// same mother, 6 shared non-photon mother, 7 both fake, 8 fake from one non-photon mother, 9 fake with a non-electron leg,
/// 10 fake from two photons, 11 other. The first test that holds decides.
inline int pcFakeClass(PcPhoton const& a, PcPhoton const& b)
{
  const std::array<const PcLeg*, 2> la{&a.pos, &a.ele}, lb{&b.pos, &b.ele};
  for (const auto* x : la) {
    for (const auto* y : lb) {
      if (x->evt == y->evt && x->trk >= 0 && x->trk == y->trk) {
        return kClsSharedTrack;
      }
    }
  }
  if (pcTruePhoton(a) && pcTruePhoton(b) && pcSameMother(a.pos, b.pos)) {
    return kClsSamePhoton;
  }
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      if (pcSamePart(*la[i], *lb[j])) {
        return (la[1 - i]->mpdg == kPcGamma && lb[1 - j]->mpdg == kPcGamma) ? kClsSharedLegPhotons : kClsSharedLeg;
      }
    }
  }
  for (const auto* x : la) {
    for (const auto* y : lb) {
      if (x->mpdg == kPcGamma && pcSameMother(*x, *y)) {
        return kClsHalfSwap;
      }
    }
  }
  for (const auto* x : la) {
    for (const auto* y : lb) {
      if (pcDaughterOf(*x, *y) || pcDaughterOf(*y, *x)) {
        return kClsBremsFake;
      }
    }
  }
  if (pcRelated(a, b)) {
    return kClsBremsFake;
  }
  if (pcSiblingPair(a, b)) {
    return kClsSameMotherFake;
  }
  for (const auto* x : la) {
    for (const auto* y : lb) {
      if (pcSameMother(*x, *y)) {
        return kClsSharedNonPhotonMother;
      }
    }
  }
  const bool aFake = !pcTruePhoton(a), bFake = !pcTruePhoton(b);
  if (aFake && bFake) {
    return kClsBothFake;
  }
  const auto& f = aFake ? a : b;
  if (pcSameMother(f.pos, f.ele)) {
    return kClsFakeNonPhotonMother;
  }
  if (std::abs(f.pos.pdg) != kPcElectron || std::abs(f.ele.pdg) != kPcElectron) {
    return kClsFakeNonElectron;
  }
  if (f.pos.mpdg == kPcGamma && f.ele.mpdg == kPcGamma) {
    return kClsFakeTwoPhotons;
  }
  return kClsFakeOther;
}

/// MC class of a photon pair (PairClass, 1 to 20): two true photons are split into 13 true distinct, 16 bremsstrahlung and
/// 18 same mother; their full swaps into 14, 17, 19; pairs with a fake are classified by pcFakeClass; 15 no label.
inline int pcPairClass(PcPhoton const& a, PcPhoton const& b)
{
  const int truth = pcPairTruth(a, b);
  if (truth == kTruthWithFake) {
    return pcFakeClass(a, b);
  }
  if (truth == kTruthTwoPhotons || truth == kTruthFullSwap) {
    if (pcRelated(a, b)) {
      return truth == kTruthTwoPhotons ? kClsBrems : kClsRelatedSwap; // bremsstrahlung
    }
    if (pcSiblingPair(a, b)) {
      return truth == kTruthTwoPhotons ? kClsSameMother : kClsSameMotherSwap; // photons of the same mother
    }
    return truth == kTruthTwoPhotons ? kClsTrueDistinct : kClsFullSwap;
  }
  return kClsNoLabel;
}

/// Origin of the two stored photons: 0 both SVertexer originals, 1 one finder original, 2 both finder, 3 a crossed pairing
inline int pcPairOrigin(PcPhoton const& a, PcPhoton const& b)
{
  if (!a.original || !b.original) {
    return 3;
  }
  return (a.svertexer ? 0 : 1) + (b.svertexer ? 0 : 1);
}

/// Partner relation of the two photons (PartnerRelation, 0 to 3); 0 when both are run-1 pairings.
inline int pcPartner(PcPhoton const& a, PcPhoton const& b)
{
  if (a.original && b.original) {
    return kPartnerOriginals;
  }
  const std::array<const PcLeg*, 2> la{&a.pos, &a.ele}, lb{&b.pos, &b.ele};
  for (const auto* x : la) {
    for (const auto* y : lb) {
      if (x->srcPhoton >= 0 && x->srcPhoton == y->srcPhoton) {
        return kPartnerSharedSource;
      }
    }
  }
  if ((a.pos.free && b.ele.free && a.pos.evt == b.ele.evt) || (a.ele.free && b.pos.free && a.ele.evt == b.pos.evt)) {
    return kPartnerFreeLegs;
  }
  return kPartnerUnrelated;
}

inline bool pcPartnerKept(const int partner, const int mode)
{
  return mode <= 0 || partner <= kPartnerSharedSource || (mode >= 2 && partner == kPartnerFreeLegs); // o2-linter: disable=magic-number (partner mode)
}

/// MC relation of a stray leg f of a half swap to the two legs g1, g2 of the shared true photon:
/// 0 conversion e+- of an unrelated photon, 1 other unrelated e+-, 2 same cascade, 3 same mother, 4 same physical primary,
/// 5 not an electron, 6 no MC label.
inline int pcStrayRelation(PcLeg const& f, PcLeg const& g1, PcLeg const& g2)
{
  constexpr int NoLabel = 6, NotElectron = 5, SamePrimary = 4, SameMother = 3, Cascade = 2, OtherElectron = 1, UnrelatedConversion = 0;
  if (f.pid < 0) {
    return NoLabel;
  }
  if (std::abs(f.pdg) != kPcElectron) {
    return NotElectron;
  }
  for (const auto* g : {&g1, &g2}) {
    if (pcSameMother(f, *g) || pcDaughterOf(f, *g) || pcRadiatedBy(f, *g) || pcRadiatedBy(*g, f)) {
      return Cascade;
    }
  }
  if (f.evt == g1.evt && g1.gmid >= 0 && (f.gmid == g1.gmid || f.mid == g1.gmid)) {
    return SameMother;
  }
  if (f.evt == g1.evt && g1.anc >= 0 && f.anc == g1.anc) {
    return SamePrimary;
  }
  return f.mpdg == kPcGamma ? UnrelatedConversion : OtherElectron;
}

/// Largest distance (cm) between the MC production vertices of the given legs, -1 if one has no MC particle.
inline float pcTrueDist(std::initializer_list<const PcLeg*> legs)
{
  float d = 0.f;
  for (auto x = legs.begin(); x != legs.end(); ++x) {
    if ((*x)->pid < 0) {
      return -1.f;
    }
    for (auto y = std::next(x); y != legs.end(); ++y) {
      d = std::max(d, std::hypot((*x)->vx - (*y)->vx, (*x)->vy - (*y)->vy, (*x)->vz - (*y)->vz));
    }
  }
  return d;
}

inline std::array<const PcLeg*, 2> pcTrueMothers(PcPhoton const& a, PcPhoton const& b)
{
  std::vector<const PcLeg*> mothers;
  for (const auto* l : {&a.pos, &a.ele, &b.pos, &b.ele}) {
    if (l->mid < 0 || l->mpdg != kPcGamma) {
      continue;
    }
    bool seen = false;
    for (const auto* m : mothers) {
      seen = seen || pcSameMother(*l, *m);
    }
    if (!seen) {
      mothers.push_back(l);
    }
  }
  if (mothers.size() != 2) { // o2-linter: disable=magic-number (two photons)
    return {nullptr, nullptr};
  }
  return {mothers[0], mothers[1]};
}

/// q_inv (GeV/c) of the two MC photons behind the four legs, -1 unless the four legs come from exactly two photons.
inline float pcTrueQ(PcPhoton const& a, PcPhoton const& b)
{
  std::vector<const PcLeg*> mothers;
  for (const auto* l : {&a.pos, &a.ele, &b.pos, &b.ele}) {
    if (l->mid < 0 || l->mpdg != kPcGamma) {
      continue;
    }
    bool seen = false;
    for (const auto* m : mothers) {
      seen = seen || pcSameMother(*l, *m);
    }
    if (!seen) {
      mothers.push_back(l);
    }
  }
  if (mothers.size() != 2) {
    return -1.f;
  }
  return o2::pwgem::photonmeson::qInv(mothers[0]->mpx, mothers[0]->mpy, mothers[0]->mpz, mothers[1]->mpx, mothers[1]->mpy, mothers[1]->mpz);
}

/// The two pieces of a duplicate: the legs of a and b that are the same AO2D track (sameTrack) or the same MC particle
inline std::tuple<const PcLeg*, const PcLeg*, int> pcSharedPieces(PcPhoton const& a, PcPhoton const& b, const bool sameTrack)
{
  const std::array<const PcLeg*, 2> la{&a.pos, &a.ele}, lb{&b.pos, &b.ele};
  for (int i = 0; i < 2; ++i) {   // o2-linter: disable=magic-number (two legs)
    for (int j = 0; j < 2; ++j) { // o2-linter: disable=magic-number (two legs)
      const bool shares = sameTrack ? (la[i]->evt == lb[j]->evt && la[i]->trk >= 0 && la[i]->trk == lb[j]->trk) : pcSamePart(*la[i], *lb[j]);
      if (shares) {
        return {la[i], lb[j], i == j ? 0 : 1};
      }
    }
  }
  return {nullptr, nullptr, 2};
}

/// The two photons share an AO2D track of the same source event (such pairs are never paired).
inline bool pcSharedTrack(PcPhoton const& a, PcPhoton const& b)
{
  for (const auto* x : {&a.pos, &a.ele}) {
    for (const auto* y : {&b.pos, &b.ele}) {
      if (x->evt == y->evt && x->trk >= 0 && x->trk == y->trk) {
        return true;
      }
    }
  }
  return false;
}

/// Both legs carry an MC label.
inline bool pcHasLabels(PcPhoton const& p) { return p.pos.pid >= 0 && p.ele.pid >= 0; }

/// Signed energy sharing (-1 to 1) of the photon: which charge carries more of the photon momentum.
inline float pcAlpha(PcPhoton const& p)
{
  const float lp = p.pos.px * p.px + p.pos.py * p.py + p.pos.pz * p.pz, le = p.ele.px * p.px + p.ele.py * p.py + p.ele.pz * p.pz;
  return (lp + le) != 0.f ? (lp - le) / (lp + le) : 0.f;
}

/// PcLeg from a V0LegsMix
template <typename TLeg, typename TMixLeg>
PcLeg pcLegFrom(TLeg const& l, TMixLeg const& m)
{
  PcLeg out;
  out.evt = m.sourceEventUid();
  out.pid = m.mcParticleId();
  out.mid = m.mcMotherId();
  out.gmid = m.mcGrandMotherId();
  out.pdg = m.mcPdg();
  out.mpdg = m.mcMotherPdg();
  out.gpdg = m.mcGrandMotherPdg();
  out.trk = m.origTrackId();
  out.anc = m.mcAncestorId();
  out.ggid = m.mcGGMotherId();
  out.root = m.mcRootId();
  out.ggpdg = m.mcGGMotherPdg();
  out.rootPdg = m.mcRootPdg();
  out.sharedFrac = l.tpcFractionSharedCls();
  out.nCls = static_cast<float>(m.tpcNClsFound());
  out.vx = m.mcVx();
  out.vy = m.mcVy();
  out.vz = m.mcVz();
  out.px = l.px();
  out.py = l.py();
  out.pz = l.pz();
  out.mpx = m.mcMotherPx();
  out.mpy = m.mcMotherPy();
  out.mpz = m.mcMotherPz();
  out.its = m.hasITS();
  out.free = m.sourceV0Id() == -2; // o2-linter: disable=magic-number (free-leg stub v0Id of the mixer)
  out.srcPhoton = m.sourcePhotonUid();
  return out;
}

template <typename TPhoton, typename TLegs, typename TMixLegs, typename TKFChi2>
PcPhoton pcPhotonFrom(TPhoton const& photon, TLegs const& legs, TMixLegs const& mixLegs, TKFChi2 const& kfChi2)
{
  PcPhoton p;
  p.px = photon.px();
  p.py = photon.py();
  p.pz = photon.pz();
  p.vx = photon.vx();
  p.vy = photon.vy();
  p.vz = photon.vz();
  p.original = photon.v0Id() >= 0;
  const auto lp = legs.rawIteratorAt(photon.posTrackId());
  const auto ln = legs.rawIteratorAt(photon.negTrackId());
  const auto mp = mixLegs.rawIteratorAt(lp.trackId());
  const auto mn = mixLegs.rawIteratorAt(ln.trackId());
  p.pos = pcLegFrom(lp, mp);
  p.ele = pcLegFrom(ln, mn);
  p.svertexer = p.original && mp.sourceV0Id() >= 0;
  p.nFree = (mp.sourceV0Id() == -2) + (mn.sourceV0Id() == -2); // o2-linter: disable=magic-number (free-leg stub v0Id of the mixer)
  p.chi2 = photon.chiSquareNDF();
  p.dcaxy = photon.dcaXYtopv();
  p.dcaz = photon.dcaZtopv();
  p.cospa = photon.cospa();
  if (photon.globalIndex() < kfChi2.size()) { // row-by-row extension of the photon table
    const auto kc = kfChi2.rawIteratorAt(photon.globalIndex());
    p.chi2Vtx = kc.chi2NDFVtx();
    p.chi2PV = kc.chi2PV();
  }
  p.psipair = photon.psipair();
  p.posMix = lp.trackId();
  p.eleMix = ln.trackId();
  return p;
}

/// MC identity of a stored photon of the classic derived tables (generated-level part)
struct PhotonMCInfo {
  bool hasMC = false;        // both legs have an MC particle
  bool sameMother = false;   // and the same mother
  bool isTruePhoton = false; // and that mother is a photon
  int mcPosId = -1;
  int mcNegId = -1;
  int motherId = -1;
  int motherPdg = 0;
};

/// PhotonMCInfo of photon g from the MC labels of its legs (TLegs carries the MC labels)
template <typename TPhoton, typename TLegs, typename TMCParticles>
PhotonMCInfo buildPhotonMCInfo(TPhoton const& g, TMCParticles const& mcParticles)
{
  PhotonMCInfo info{};
  const auto pos = g.template posTrack_as<TLegs>();
  const auto neg = g.template negTrack_as<TLegs>();
  if (!pos.has_emmcparticle() || !neg.has_emmcparticle()) {
    return info;
  }
  info.hasMC = true;
  info.mcPosId = pos.emmcparticleId();
  info.mcNegId = neg.emmcparticleId();
  const auto mcPos = pos.template emmcparticle_as<TMCParticles>();
  const auto mcNeg = neg.template emmcparticle_as<TMCParticles>();
  if (!mcPos.has_mothers() || !mcNeg.has_mothers()) {
    return info;
  }
  const int mothIdPos = mcPos.mothersIds()[0], mothIdNeg = mcNeg.mothersIds()[0];
  if (mothIdPos != mothIdNeg) {
    return info;
  }
  info.sameMother = true;
  info.motherId = mothIdPos;
  const auto mother = mcParticles.iteratorAt(mothIdPos);
  info.motherPdg = mother.pdgCode();
  info.isTruePhoton = (info.motherPdg == kGamma);
  return info;
}

struct TruthGamma {
  int id = -1, posId = -1, negId = -1;
  float eta = 0.f, phi = 0.f, pt = 0.f;
  float rTrue = -1.f;
  std::array<float, 3> conv{}; // true conversion point (production vertex of the e+)
};

/// radii (m) at which the legs are propagated for the same-sign leg separation
inline constexpr std::array<float, 17> kPhiStarRadiiM = {
  0.85f, 0.95f, 1.05f, 1.15f, 1.25f, 1.35f, 1.45f, 1.55f, 1.65f,
  1.75f, 1.85f, 1.95f, 2.05f, 2.15f, 2.25f, 2.35f, 2.45f};
inline constexpr float kPhiStarInvalid = -999.f;
inline constexpr float kCmPerM = 100.f;

struct LegHelix {
  float phiStar = kPhiStarInvalid; // azimuth of the track at that radius
  float sxy = -1.f;
  int status = 3; // 0 ok, 1 below R_conv, 2 out of reach,3 bad input
};

[[nodiscard]] inline LegHelix legHelixAt(float vxM, float vyM, float phiMom, float pt, int charge, float bzT, float radiusM)
{
  LegHelix out;
  if (pt < kMinSigma || std::fabs(bzT) < kMinSigma) {
    out.status = 3; // o2-linter: disable=magic-number (status code)
    return out;
  }
  const float rStart = std::hypot(vxM, vyM);
  if (rStart < kMinSigma) {
    out.status = 3; // o2-linter: disable=magic-number (status code)
    return out;
  }
  if (radiusM < rStart) {
    out.status = 1;
    return out;
  }
  const float rho = pt / (0.3f * std::fabs(bzT)); // curvature radius, m
  const float sgn = (bzT >= 0.f ? 1.f : -1.f) * static_cast<float>(charge);
  const float cx = vxM + sgn * rho * std::sin(phiMom);
  const float cy = vyM - sgn * rho * std::cos(phiMom);
  const float d = std::hypot(cx, cy);
  if (d < kMinSigma) {
    out.status = 3; // o2-linter: disable=magic-number (status code)
    return out;
  }
  if (radiusM > d + rho || radiusM < std::fabs(d - rho)) {
    out.status = 2; // o2-linter: disable=magic-number (status code)
    return out;
  }
  const float phiC = std::atan2(cy, cx);
  auto clamp1 = [](float x) { return std::max(-1.f, std::min(1.f, x)); };
  auto alphaAt = [&](float r) { return std::acos(clamp1((r * r + d * d - rho * rho) / (2.f * r * d))); };
  auto betaAt = [&](float r) { return std::acos(clamp1((d * d + rho * rho - r * r) / (2.f * d * rho))); };

  const float phiVtx = std::atan2(vyM, vxM);
  const float a0 = alphaAt(rStart);
  const float resPlus = std::fabs(RecoDecay::constrainAngle(phiC + a0 - phiVtx, -o2::constants::math::PI));
  const float resMinus = std::fabs(RecoDecay::constrainAngle(phiC - a0 - phiVtx, -o2::constants::math::PI));
  const float branch = (resPlus < resMinus) ? 1.f : -1.f;

  out.phiStar = RecoDecay::constrainAngle(phiC + branch * alphaAt(radiusM), 0.f);
  out.sxy = rho * std::fabs(betaAt(radiusM) - betaAt(rStart));
  out.status = 0;
  return out;
}

[[nodiscard]] inline float selfTestLegHelix()
{
  float worst = 0.f;
  constexpr float eps = 1.e-5f;

  for (const auto& bz : {0.5f, -0.5f}) {
    for (const auto& q : {+1, -1}) {
      for (const auto& pt : {0.1f, 0.2f, 0.5f}) {
        for (const auto& r : {0.85f, 1.45f, 2.45f}) {
          for (const auto& phi : {0.f, 1.f, -2.f}) {
            const auto h = legHelixAt(eps * std::cos(phi), eps * std::sin(phi), phi, pt, q, bz, r);

            if (h.status != 0) {
              continue;
            }
            const float arg =
              -0.3f * bz * static_cast<float>(q) * r / (2.f * pt);

            if (std::fabs(arg) >= 1.f) {
              continue;
            }
            const float referencePhi =
              RecoDecay::constrainAngle(phi + std::asin(arg), 0.f);

            worst = std::max(
              worst,
              std::fabs(RecoDecay::constrainAngle(
                h.phiStar - referencePhi,
                -o2::constants::math::PI)));
          }
        }
      }
    }
  }

  return worst;
}

/// A photon with its two legs (index 0 = e+, 1 = e-) and, if propagated, phi* and the arc length of each leg at the radii
struct PhotonWithLegs {
  float fPt{0}, fEta{0}, fPhi{0};
  float fVx{0}, fVy{0}, fVz{0};
  float fVtxZ{0};
  std::array<float, 2> fLegPt{}, fLegEta{}, fLegPhi{};
  std::array<std::array<float, kPhiStarRadiiM.size()>, 2> fLegPhiStar{};
  std::array<std::array<float, kPhiStarRadiiM.size()>, 2> fLegSxy{}; // transverse arc length since conversion (m)
  [[nodiscard]] float pt() const { return fPt; }
  [[nodiscard]] float eta() const { return fEta; }
  [[nodiscard]] float phi() const { return fPhi; }
  [[nodiscard]] float vx() const { return fVx; }
  [[nodiscard]] float vy() const { return fVy; }
  [[nodiscard]] float vz() const { return fVz; }
  o2::aod::pwgem::photonmeson::utils::pairutil::V0PhotonLegCounts fLegCounts{};
  [[nodiscard]] o2::aod::pwgem::photonmeson::utils::pairutil::V0PhotonLegCounts const& legCounts() const { return fLegCounts; }
  int fIsTruePhoton{-1};
};

/// PhotonWithLegs of photon g with legs pos, ele in an event with vertex z vtxZ (cm).
template <typename TGamma, typename TLeg>
[[nodiscard]] PhotonWithLegs makePhotonWithLegs(TGamma const& g, TLeg const& pos, TLeg const& ele, float vtxZ, float bzT, bool propagate)
{
  PhotonWithLegs p;
  p.fVtxZ = vtxZ;
  p.fPt = g.pt();
  p.fEta = g.eta();
  p.fPhi = g.phi();
  p.fVx = g.vx();
  p.fVy = g.vy();
  p.fVz = g.vz();
  p.fLegPt = {static_cast<float>(pos.pt()), static_cast<float>(ele.pt())};
  p.fLegEta = {static_cast<float>(pos.eta()), static_cast<float>(ele.eta())};
  p.fLegPhi = {static_cast<float>(pos.phi()), static_cast<float>(ele.phi())};
  p.fLegCounts = o2::aod::pwgem::photonmeson::utils::pairutil::getV0PhotonLegCounts(pos, ele);
  // without propagation no radius is valid: computePairSep finds no common radius (nCommon 0), never a distance of 0
  if (!propagate) {
    for (int il = 0; il < 2; ++il) { // o2-linter: disable=magic-number (two legs)
      p.fLegPhiStar[il].fill(kPhiStarInvalid);
      p.fLegSxy[il].fill(-1.f); // arc length of a radius that was not computed, as in LegHelix
    }
    return p;
  }
  // legHelixAt itself rejects radii below the conversion point
  for (size_t ir = 0; ir < kPhiStarRadiiM.size(); ++ir) {
    for (int il = 0; il < 2; ++il) { // o2-linter: disable=magic-number (two legs)
      const auto h = legHelixAt(p.fVx / kCmPerM, p.fVy / kCmPerM, p.fLegPhi[il], p.fLegPt[il],
                                (il == 0) ? +1 : -1, bzT, kPhiStarRadiiM[ir]);
      p.fLegPhiStar[il][ir] = h.phiStar;
      p.fLegSxy[il][ir] = h.sxy;
    }
  }
  return p;
}

/// Scales of the pair-merging
struct PairSepConfig {
  float sigR = 1.f, sigZ = 1.f, dClose = 1.f;
  std::array<float, 4> closeThr{};
};

/// Same-sign leg separation of a pair
struct PairSep {
  float dMinRPhi{999.f};
  std::array<float, 4> fCloseRPhiThr{};
  int nCommon{0};
  int rphiCharge{-1}; // chosen by d_rphi alone
  std::array<std::array<float, 4>, 2> fCloseRPhiThrPerCharge{};
  float dPairMin{999.f};
  int nCloseScaled{0};
  std::array<int, 2> nPoints{};
};

/// PairSep of two propagated photons (makePhotonWithLegs with propagate).
[[nodiscard]] inline PairSep computePairSep(PhotonWithLegs const& a, PhotonWithLegs const& b, PairSepConfig const& cfg)
{
  PairSep s;

  const float sigR = cfg.sigR;
  const float sigZ = cfg.sigZ;
  const float dClose = cfg.dClose;
  for (int i = 0; i < 2; ++i) { // o2-linter: disable=magic-number (PP and NN)
    int nHere = 0, nSame = 0, nCloseD = 0;
    std::array<int, 4> nCloseThr{};
    float bestR = 999.f, bestD = 999.f;
    for (size_t ir = 0; ir < kPhiStarRadiiM.size(); ++ir) {
      const float p1 = a.fLegPhiStar[i][ir], p2 = b.fLegPhiStar[i][ir];
      if (p1 < -100.f || p2 < -100.f) { // o2-linter: disable=magic-number (invalid phi*)
        continue;
      }
      const float rCm = kPhiStarRadiiM[ir] * kCmPerM;
      const float dphi = RecoDecay::constrainAngle(p1 - p2, -o2::constants::math::PI);
      const float dRPhi = 2.f * rCm * std::fabs(std::sin(0.5f * dphi));
      const float z1 = a.fVz + a.fLegSxy[i][ir] * kCmPerM * std::sinh(a.fLegEta[i]);
      const float z2 = b.fVz + b.fLegSxy[i][ir] * kCmPerM * std::sinh(b.fLegEta[i]);
      const float dZSignedLocal = (z1 - a.fVtxZ) - (z2 - b.fVtxZ);
      const float dZ = std::fabs(dZSignedLocal);
      const bool sameSide = (z1 * z2 >= 0.f);
      const float dPair = std::hypot(dRPhi / sigR, dZ / sigZ);
      ++nHere;
      if (!sameSide) {
        continue;
      }
      ++nSame;
      for (std::size_t t = 0; t < nCloseThr.size(); ++t) {
        if (dRPhi < cfg.closeThr[t]) {
          ++nCloseThr[t];
        }
      }
      if (dPair < dClose) {
        ++nCloseD;
      }
      if (dPair < bestD) {
        bestD = dPair;
      }
      if (dRPhi < bestR) {
        bestR = dRPhi;
      }
    }
    if (nHere == 0) {
      continue;
    }
    s.nPoints[i] = nHere;
    const float invS = (nSame > 0) ? 1.f / static_cast<float>(nSame) : 0.f;
    for (std::size_t t = 0; t < nCloseThr.size(); ++t) {
      s.fCloseRPhiThrPerCharge[i][t] = static_cast<float>(nCloseThr[t]) * invS;
    }
    if (bestD < s.dPairMin) {
      s.dPairMin = bestD;
      s.nCloseScaled = nCloseD;
    }
    if (bestR < s.dMinRPhi) {
      s.dMinRPhi = bestR;
      s.fCloseRPhiThr = s.fCloseRPhiThrPerCharge[i];
      s.rphiCharge = i;
    }
  }
  s.nCommon = (s.rphiCharge >= 0) ? s.nPoints[s.rphiCharge]
                                  : std::max(s.nPoints[0], s.nPoints[1]);
  return s;
}

struct CrossObs {
  std::array<float, 2> mee{999.f, 999.f};
  float meeOverQ{999.f};
  std::array<float, 2> gap{999.f, 999.f};
  std::array<float, 2> dTgl{999.f, 999.f};
  std::array<float, 2> dz{999.f, 999.f};
  bool geomPhotonLike{false};
};

[[nodiscard]] inline CrossObs computeCrossObs(PhotonWithLegs const& a, PhotonWithLegs const& b, float qinv, float bzT, float maxGap, float maxDTgl, float maxDz)
{
  constexpr float kMe = 0.000510999f; // electron mass, GeV/c^2
  CrossObs c;
  auto legVec = [](PhotonWithLegs const& p, int i) {
    return ROOT::Math::PtEtaPhiMVector(p.fLegPt[i], p.fLegEta[i], p.fLegPhi[i], kMe);
  };
  // Combo 0: a's e+ with b's e-.  Combo 1: b's e+ with a's e-.
  for (int combo = 0; combo < 2; ++combo) { // o2-linter: disable=magic-number (the two leg-swap combinations)
    PhotonWithLegs const& pPos = (combo == 0) ? a : b;
    PhotonWithLegs const& pEle = (combo == 0) ? b : a;
    c.mee[combo] = static_cast<float>((legVec(pPos, 0) + legVec(pEle, 1)).M());
  }
  if (qinv > 0.f) {
    c.meeOverQ = std::min(c.mee[0], c.mee[1]) / qinv;
  }
  // geometry: helix circle of leg i (cm) from the conversion point and the leg momentum, as legHelixAt
  auto circle = [&](PhotonWithLegs const& p, int i, float& xc, float& yc, float& r) {
    r = p.fLegPt[i] / (0.3f * std::fabs(bzT)) * kCmPerM; // o2-linter: disable=magic-number (0.3 GeV/c per T m)
    const float sgn = (bzT >= 0.f ? 1.f : -1.f) * (i == 0 ? 1.f : -1.f);
    xc = p.fVx + sgn * r * std::sin(p.fLegPhi[i]);
    yc = p.fVy - sgn * r * std::cos(p.fLegPhi[i]);
  };
  // z of leg i at the point of its circle closest to (xT, yT): signed arc along the motion from the conversion point, local z
  auto zAt = [&](PhotonWithLegs const& p, int i, float xc, float yc, float r, float xT, float yT) {
    const float phiStart = std::atan2(p.fVy - yc, p.fVx - xc);
    const float phiT = std::atan2(yT - yc, xT - xc);
    const bool ccw = std::cos(p.fLegPhi[i]) * (-(p.fVy - yc)) + std::sin(p.fLegPhi[i]) * (p.fVx - xc) > 0.f; // momentum along the counterclockwise tangent
    const float dphi = RecoDecay::constrainAngle(phiT - phiStart, -o2::constants::math::PI);
    return (p.fVz - p.fVtxZ) + r * (ccw ? dphi : -dphi) * std::sinh(p.fLegEta[i]);
  };
  if (std::fabs(bzT) > kMinSigma) {
    for (int combo = 0; combo < 2; ++combo) { // o2-linter: disable=magic-number (the two leg-swap combinations)
      PhotonWithLegs const& pPos = (combo == 0) ? a : b;
      PhotonWithLegs const& pEle = (combo == 0) ? b : a;
      float xP = 0.f, yP = 0.f, rP = 0.f, xE = 0.f, yE = 0.f, rE = 0.f;
      circle(pPos, 0, xP, yP, rP);
      circle(pEle, 1, xE, yE, rE);
      const float dC = std::hypot(xE - xP, yE - yP);
      if (!(dC > kMinSigma) || !(rP + rE > kMinSigma)) {
        continue;
      }
      c.gap[combo] = dC - (rP + rE);
      c.dTgl[combo] = std::fabs(std::sinh(pPos.fLegEta[0]) - std::sinh(pEle.fLegEta[1]));
      const float w = rP / (rP + rE);
      const float xT = xP + w * (xE - xP), yT = yP + w * (yE - yP); // touching point of the circles
      c.dz[combo] = zAt(pPos, 0, xP, yP, rP, xT, yT) - zAt(pEle, 1, xE, yE, rE, xT, yT);
      if (std::fabs(c.gap[combo]) < maxGap && c.dTgl[combo] < maxDTgl && std::fabs(c.dz[combo]) < maxDz) {
        c.geomPhotonLike = true;
      }
    }
  }
  return c;
}

/// Minimum distances (cm) of the same-sign legs of two co-reconstructed photons in the TPC, same-side radii only:
struct PcPairSep {
  float dMin = 999.f, dRPhiMin = 999.f;
};

[[nodiscard]] inline PcPairSep pcPairSep(PcPhoton const& a, PcPhoton const& b, float bzT)
{
  PcPairSep s;
  if (std::fabs(bzT) < kMinSigma) {
    return s;
  }
  const std::array<const PcLeg*, 2> la{&a.pos, &a.ele};
  const std::array<const PcLeg*, 2> lb{&b.pos, &b.ele};
  for (int i = 0; i < 2; ++i) { // o2-linter: disable=magic-number (two same-sign leg pairs)
    const int charge = i == 0 ? +1 : -1;
    const float pt1 = std::hypot(la[i]->px, la[i]->py), pt2 = std::hypot(lb[i]->px, lb[i]->py);
    if (pt1 < kMinSigma || pt2 < kMinSigma) {
      continue;
    }
    const float phi1 = std::atan2(la[i]->py, la[i]->px), phi2 = std::atan2(lb[i]->py, lb[i]->px);
    const float tgl1 = la[i]->pz / pt1, tgl2 = lb[i]->pz / pt2;
    for (size_t ir = 0; ir < kPhiStarRadiiM.size(); ++ir) {
      const auto h1 = legHelixAt(a.vx / kCmPerM, a.vy / kCmPerM, phi1, pt1, charge, bzT, kPhiStarRadiiM[ir]);
      const auto h2 = legHelixAt(b.vx / kCmPerM, b.vy / kCmPerM, phi2, pt2, charge, bzT, kPhiStarRadiiM[ir]);
      if (h1.status != 0 || h2.status != 0) {
        continue;
      }
      const float rCm = kPhiStarRadiiM[ir] * kCmPerM;
      const float dphi = RecoDecay::constrainAngle(h1.phiStar - h2.phiStar, -o2::constants::math::PI);
      const float dRPhi = 2.f * rCm * std::fabs(std::sin(0.5f * dphi));
      const float z1 = a.vz + h1.sxy * kCmPerM * tgl1;
      const float z2 = b.vz + h2.sxy * kCmPerM * tgl2;
      if (z1 * z2 < 0.f) {
        continue;
      }
      s.dRPhiMin = std::min(s.dRPhiMin, dRPhi);
      s.dMin = std::min(s.dMin, std::hypot(dRPhi, std::fabs(z1 - z2)));
    }
  }
  return s;
}

template <typename TMixLegs>
float pcBremsPointingDMin(PcPhoton const& a, PcPhoton const& b, TMixLegs const* mixLegs, float bzT)
{
  constexpr float StepCm = 0.25f, MaxArcCm = 100.f, MarginCm = 2.f;
  float best = 999.f;
  if (mixLegs == nullptr) {
    return best;
  }
  const float bzKG = 10.f * bzT;
  for (int order = 0; order < 2; ++order) {
    PcPhoton const& emitter = order == 0 ? a : b;
    PcPhoton const& brems = order == 0 ? b : a;
    const float pB = std::hypot(brems.px, brems.py, brems.pz);
    if (pB <= 0.f) {
      continue;
    }
    const std::array<float, 3> u{brems.px / pB, brems.py / pB, brems.pz / pB};
    const float rA = std::hypot(emitter.vx, emitter.vy), rB = std::hypot(brems.vx, brems.vy);
    for (const int64_t row : {emitter.posMix, emitter.eleMix}) {
      if (row < 0 || row >= mixLegs->size()) {
        continue;
      }
      const auto leg = mixLegs->rawIteratorAt(row);
      const o2::track::TrackPar par(leg.x(), leg.alpha(), std::array<float, 5>{leg.y(), leg.z(), leg.snp(), leg.tgl(), leg.signed1Pt()});
      o2::math_utils::CircleXYf_t circle;
      float sna = 0.f, csa = 0.f;
      par.getCircleParams(bzKG, circle, sna, csa);
      std::array<float, 3> x0{}, p0{};
      par.getXYZGlo(x0);
      par.getPxPyPzGlo(p0);
      if (circle.rC <= 0.f) {
        continue;
      }
      const float psi0 = std::atan2(x0[1] - circle.yC, x0[0] - circle.xC);
      const float sense = ((x0[0] - circle.xC) * p0[1] - (x0[1] - circle.yC) * p0[0]) > 0.f ? 1.f : -1.f;
      const float tgl = par.getTgl();
      const int nArcSteps = static_cast<int>(std::lround(2.f * MaxArcCm / StepCm));
      for (int iArc = 0; iArc <= nArcSteps; ++iArc) {
        const float sArc = -MaxArcCm + static_cast<float>(iArc) * StepCm;
        const float psi = psi0 + sense * sArc / circle.rC;
        const std::array<float, 3> h{circle.xC + circle.rC * std::cos(psi), circle.yC + circle.rC * std::sin(psi), x0[2] + tgl * sArc};
        const float rH = std::hypot(h[0], h[1]);
        if (rH < rA - MarginCm || rH > rB + MarginCm) {
          continue; // emission between the two conversions
        }
        const std::array<float, 3> d{h[0] - brems.vx, h[1] - brems.vy, h[2] - brems.vz};
        const float t = d[0] * u[0] + d[1] * u[1] + d[2] * u[2];
        if (t > 0.f) {
          continue; // in front of B's conversion point
        }
        best = std::min(best, std::hypot(d[0] - t * u[0], d[1] - t * u[1], d[2] - t * u[2]));
      }
    }
  }
  return best;
}

inline constexpr float kPointTagVeryTight = 0.1f, kPointTagTight = 0.3f, kPointTagLoose = 1.f, kPointTagNone = 100.f;

inline float pcPointingCategory(const float d)
{
  if (d < kPointTagVeryTight) {
    return 0.f;
  }
  if (d < kPointTagTight) {
    return 1.f;
  }
  if (d < kPointTagLoose) {
    return 2.f;
  }
  return d < kPointTagNone ? 3.f : 4.f; // o2-linter: disable=magic-number (category codes)
}

/// Photon-likeness of the two leg-swapped combinations of a pair from the leg track parameters
template <typename TMixLegs>
std::array<float, 2> pcCrossScores(PcPhoton const& a, PcPhoton const& b, TMixLegs const* mixLegs, float bzT, float maxGap, float maxDz, float maxDTgl)
{
  std::array<float, 2> s{999.f, 999.f};
  if (mixLegs == nullptr || std::fabs(bzT) < kMinSigma) {
    return s;
  }
  struct Helix {
    bool ok = false;
    float xC = 0.f, yC = 0.f, rC = 0.f, psi0 = 0.f, sense = 1.f, z0 = 0.f, tgl = 0.f;
  };
  const float bzKG = 10.f * bzT;
  auto helix = [&](const int64_t row) {
    Helix h;
    if (row < 0 || row >= mixLegs->size()) {
      return h;
    }
    const auto leg = mixLegs->rawIteratorAt(row);
    const o2::track::TrackPar par(leg.x(), leg.alpha(), std::array<float, 5>{leg.y(), leg.z(), leg.snp(), leg.tgl(), leg.signed1Pt()});
    o2::math_utils::CircleXYf_t circle;
    float sna = 0.f, csa = 0.f;
    par.getCircleParams(bzKG, circle, sna, csa);
    if (circle.rC <= 0.f) {
      return h;
    }
    std::array<float, 3> x0{}, p0{};
    par.getXYZGlo(x0);
    par.getPxPyPzGlo(p0);
    h.xC = circle.xC;
    h.yC = circle.yC;
    h.rC = circle.rC;
    h.psi0 = std::atan2(x0[1] - circle.yC, x0[0] - circle.xC);
    h.sense = ((x0[0] - circle.xC) * p0[1] - (x0[1] - circle.yC) * p0[0]) > 0.f ? 1.f : -1.f;
    h.z0 = x0[2];
    h.tgl = par.getTgl();
    h.ok = true;
    return h;
  };
  const float gapMax = std::max(maxGap, kMinSigma), dzMax = std::max(maxDz, kMinSigma), dTglMax = std::max(maxDTgl, kMinSigma);
  const std::array<std::pair<int64_t, int64_t>, 2> combos{{{a.posMix, b.eleMix}, {b.posMix, a.eleMix}}};
  for (size_t c = 0; c < combos.size(); ++c) {
    const Helix hP = helix(combos[c].first), hE = helix(combos[c].second);
    const float dC = std::hypot(hE.xC - hP.xC, hE.yC - hP.yC);
    if (!hP.ok || !hE.ok || !(dC > kMinSigma)) {
      continue;
    }
    const float w = hP.rC / (hP.rC + hE.rC);
    const float xT = hP.xC + w * (hE.xC - hP.xC), yT = hP.yC + w * (hE.yC - hP.yC); // touching point of the two circles
    auto zAt = [&](Helix const& h) {
      const float dpsi = RecoDecay::constrainAngle(std::atan2(yT - h.yC, xT - h.xC) - h.psi0, -o2::constants::math::PI);
      return h.z0 + h.tgl * h.sense * dpsi * h.rC;
    };
    s[c] = std::max({std::fabs(dC - (hP.rC + hE.rC)) / gapMax, std::fabs(zAt(hP) - zAt(hE)) / dzMax, std::fabs(hP.tgl - hE.tgl) / dTglMax});
  }
  if (s[1] < s[0]) {
    std::swap(s[0], s[1]);
  }
  return s;
}

} // namespace o2::aod::pwgem::photonmeson::hbt

#endif // PWGEM_PHOTONMESON_UTILS_PHOTONHBTUTILITIES_H_
