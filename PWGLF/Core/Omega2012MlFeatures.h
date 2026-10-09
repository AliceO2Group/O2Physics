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
/// \file Omega2012MlFeatures.h
/// \brief Snapshots, canonical candidates and the two separate feature contracts of the Omega(2012) decay modes
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>
///
/// Mode A (XiK0s): Omega(2012)- -> Xi- K0S. Mode B (Xi1530K): Omega(2012)- -> Xi(1530)0 K- -> Xi- pi+ K-.
/// The two modes have independent, separately versioned contracts (feature_contract_omega2012_v1.json);
/// no feature vector, status or projection is shared between them.

#ifndef PWGLF_CORE_OMEGA2012MLFEATURES_H_
#define PWGLF_CORE_OMEGA2012MLFEATURES_H_

#include "PWGLF/DataModel/LFResonanceTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>

#include <Math/GenVector/Boost.h>
#include <Math/GenVector/VectorUtil.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string_view>
#include <vector>

namespace o2::analysis::omega2012ml
{
inline constexpr std::size_t NItsLayers = 7;
inline constexpr std::size_t NCascadeDaughters = 3;
inline constexpr std::size_t NV0Daughters = 2;
inline constexpr std::size_t NPidSpecies = 3;
inline constexpr std::size_t XiK0sNMasterFeatures = 63;
inline constexpr std::size_t Xi1530KNMasterFeatures = 126;
inline constexpr int8_t SaturatedLow = std::numeric_limits<int8_t>::lowest(); // x10 nSigma code of -inf / missing TOF
inline constexpr int8_t SaturatedHigh = std::numeric_limits<int8_t>::max();   // x10 nSigma code of NaN / +inf
inline constexpr float NSigmaScale = 10.f;
inline constexpr float MicroPidOverflow = 3.5f;

enum class Profile : uint8_t { DetectorV1,
                               RelationalV1,
                               StudyV1 };
enum class BuildStatus : uint8_t {
  Ok,
  InvalidChargePattern,
  ReusedDaughter,
  InvalidMomentum,
  InvalidKinematics,
  InvalidContract
};

/// Charge pattern of a mode-B candidate relative to the Xi charge (metadata, never a feature).
enum ChargePattern : uint8_t {
  kSignalPattern = 0, // pion opposite to the Xi charge, kaon equal to it (Xi- pi+ K- and conjugate)
  kWrongSignPion = 1, // pion with the Xi charge
  kWrongSignKaon = 2, // kaon opposite to the Xi charge
  kWrongSignBoth = 3  // both wrong
};

/// All ResoCascades scalars of one cascade plus derived decay lengths.
/// Daughter arrays: index 0 positive, 1 negative, 2 bachelor; nSigma: [daughter * 3 + species], species pi, ka, pr.
struct CascadeSnapshot {
  int64_t sourceRow = -1; // row of the cascade in the input ResoCascades table
  int8_t sign = 0;
  float px = 0.f, py = 0.f, pz = 0.f;
  float mXi = 0.f, mLambda = 0.f;
  float v0CosPA = 0.f, cascCosPA = 0.f, v0DaughDCA = 0.f, cascDaughDCA = 0.f;
  float dcaPosToPV = 0.f, dcaNegToPV = 0.f, dcaBachToPV = 0.f, dcaV0ToPV = 0.f, dcaXYToPV = 0.f, dcaZToPV = 0.f;
  float v0Radius = 0.f, cascRadius = 0.f;
  float decayVtxX = 0.f, decayVtxY = 0.f, decayVtxZ = 0.f;
  std::array<int8_t, NCascadeDaughters * NPidSpecies> tpcNSigma10{};
  std::array<int8_t, NCascadeDaughters * NPidSpecies> tofNSigma10{};
  std::array<uint8_t, NCascadeDaughters> crossedRows{};
  std::array<int, NCascadeDaughters> daughterIds{};
  float decayLength = 0.f;  // 3D distance from the primary vertex to the decay vertex
  float properLength = 0.f; // decayLength * MassXiMinus / p
};

/// All ResoV0s scalars of one V0 plus derived decay lengths. Daughter arrays: index 0 positive, 1 negative.
struct V0Snapshot {
  int64_t sourceRow = -1; // row of the V0 in the input ResoV0s table
  float px = 0.f, py = 0.f, pz = 0.f;
  float mK0Short = 0.f, mLambda = 0.f, mAntiLambda = 0.f;
  float cosPA = 0.f, daughDCA = 0.f, dcaPosToPV = 0.f, dcaNegToPV = 0.f, dcaV0ToPV = 0.f, radius = 0.f;
  float decayVtxX = 0.f, decayVtxY = 0.f, decayVtxZ = 0.f;
  float alpha = 0.f, qtarm = 0.f;
  std::array<int8_t, NV0Daughters * NPidSpecies> tpcNSigma10{};
  std::array<int8_t, NV0Daughters * NPidSpecies> tofNSigma10{};
  std::array<uint8_t, NV0Daughters> crossedRows{};
  std::array<int, NV0Daughters> daughterIds{};
  float decayLength = 0.f;  // 3D distance from the primary vertex to the decay vertex
  float properLength = 0.f; // decayLength * MassK0Short / p (proper lifetime c*tau)
};

/// ResoMicroTracks_001 row: raw packed bytes and the decoded public accessors.
/// The encoding rules are those of the K1 micro001 contract; the type is separate so the K1 contract stays immutable.
struct TrackSnapshot {
  int64_t sourceTrackId = -1;
  int8_t charge = 0;
  float px = 0.f, py = 0.f, pz = 0.f;
  uint8_t pidNSigmaPiFlag = 0, pidNSigmaKaFlag = 0, pidNSigmaPrFlag = 0;
  uint8_t trackSelectionFlags = 0, trackFlags = 0, tpcNClsCrossedRows = 0, itsClusterMap = 0;
  std::array<float, NPidSpecies> tpcNSigma{}; // pi, ka, pr
  std::array<float, NPidSpecies> tofNSigma{}; // pi, ka, pr
  float dcaXY = std::numeric_limits<float>::quiet_NaN();
  float dcaZ = std::numeric_limits<float>::quiet_NaN();
  bool passedPtDependentDCAxy = false, passedPtDependentDCAz = false;
  bool hasTOF = false, isPVContributor = false;
  std::array<bool, NItsLayers> itsHit{};
};

struct XiK0sCandidateSnapshot {
  CascadeSnapshot xi{};
  V0Snapshot k0s{};
};

struct Xi1530KCandidateSnapshot {
  CascadeSnapshot xi{};
  TrackSnapshot pion{};
  TrackSnapshot kaon{};
  uint8_t chargePattern = kSignalPattern; // metadata
};

template <typename Snapshot>
struct CanonicalizationResult {
  Snapshot candidate{};
  BuildStatus status = BuildStatus::InvalidChargePattern;
  explicit operator bool() const { return status == BuildStatus::Ok; }
};

/// Candidate kinematics from the stored 3-momenta and PDG masses (audit only, never features)
struct KinematicAudit {
  float mass = 0.f;
  float pt = 0.f;
  float eta = 0.f;
  float phi = 0.f;
  float scalarSumPt = 0.f;
};

template <std::size_t N>
struct FeaturePack {
  std::array<float, N> master{};
  BuildStatus status = BuildStatus::InvalidContract;
  KinematicAudit kinematics{};
};
using XiK0sFeaturePack = FeaturePack<XiK0sNMasterFeatures>;
using Xi1530KFeaturePack = FeaturePack<Xi1530KNMasterFeatures>;

namespace detail
{
template <typename Collision, typename Row>
float decayLength(Collision const& collision, Row const& row)
{
  const float dx = row.decayVtxX() - collision.posX();
  const float dy = row.decayVtxY() - collision.posY();
  const float dz = row.decayVtxZ() - collision.posZ();
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

inline float properLength(float length, float px, float py, float pz, double mass)
{
  const float p = std::sqrt(px * px + py * py + pz * pz);
  if (!(p > 0.f)) {
    return std::numeric_limits<float>::quiet_NaN();
  }
  return static_cast<float>(length * mass / p);
}
} // namespace detail

template <typename Collision, typename Cascade>
CascadeSnapshot makeCascadeSnapshot(Collision const& collision, Cascade const& row)
{
  CascadeSnapshot out;
  out.sourceRow = static_cast<int64_t>(row.globalIndex());
  out.sign = static_cast<int8_t>(row.sign());
  out.px = row.px();
  out.py = row.py();
  out.pz = row.pz();
  out.mXi = row.mXi();
  out.mLambda = row.mLambda();
  out.v0CosPA = row.v0CosPA();
  out.cascCosPA = row.cascCosPA();
  out.v0DaughDCA = row.daughDCA();
  out.cascDaughDCA = row.cascDaughDCA();
  out.dcaPosToPV = row.dcapostopv();
  out.dcaNegToPV = row.dcanegtopv();
  out.dcaBachToPV = row.dcabachtopv();
  out.dcaV0ToPV = row.dcav0topv();
  out.dcaXYToPV = row.dcaXYCascToPV();
  out.dcaZToPV = row.dcaZCascToPV();
  out.v0Radius = row.transRadius();
  out.cascRadius = row.cascTransRadius();
  out.decayVtxX = row.decayVtxX();
  out.decayVtxY = row.decayVtxY();
  out.decayVtxZ = row.decayVtxZ();
  out.tpcNSigma10 = {row.daughterTPCNSigmaPosPi10(), row.daughterTPCNSigmaPosKa10(), row.daughterTPCNSigmaPosPr10(),
                     row.daughterTPCNSigmaNegPi10(), row.daughterTPCNSigmaNegKa10(), row.daughterTPCNSigmaNegPr10(),
                     row.daughterTPCNSigmaBachPi10(), row.daughterTPCNSigmaBachKa10(), row.daughterTPCNSigmaBachPr10()};
  out.tofNSigma10 = {row.daughterTOFNSigmaPosPi10(), row.daughterTOFNSigmaPosKa10(), row.daughterTOFNSigmaPosPr10(),
                     row.daughterTOFNSigmaNegPi10(), row.daughterTOFNSigmaNegKa10(), row.daughterTOFNSigmaNegPr10(),
                     row.daughterTOFNSigmaBachPi10(), row.daughterTOFNSigmaBachKa10(), row.daughterTOFNSigmaBachPr10()};
  out.crossedRows = {row.nCrossedRowsPos(), row.nCrossedRowsNeg(), row.nCrossedRowsBach()};
  const auto indices = row.cascadeIndices();
  out.daughterIds = {indices[0], indices[1], indices[2]};
  out.decayLength = detail::decayLength(collision, row);
  out.properLength = detail::properLength(out.decayLength, out.px, out.py, out.pz, o2::constants::physics::MassXiMinus);
  return out;
}

template <typename Collision, typename V0>
V0Snapshot makeV0Snapshot(Collision const& collision, V0 const& row)
{
  V0Snapshot out;
  out.sourceRow = static_cast<int64_t>(row.globalIndex());
  out.px = row.px();
  out.py = row.py();
  out.pz = row.pz();
  out.mK0Short = row.mK0Short();
  out.mLambda = row.mLambda();
  out.mAntiLambda = row.mAntiLambda();
  out.cosPA = row.v0CosPA();
  out.daughDCA = row.daughDCA();
  out.dcaPosToPV = row.dcapostopv();
  out.dcaNegToPV = row.dcanegtopv();
  out.dcaV0ToPV = row.dcav0topv();
  out.radius = row.transRadius();
  out.decayVtxX = row.decayVtxX();
  out.decayVtxY = row.decayVtxY();
  out.decayVtxZ = row.decayVtxZ();
  out.alpha = row.alpha();
  out.qtarm = row.qtarm();
  out.tpcNSigma10 = {row.daughterTPCNSigmaPosPi10(), row.daughterTPCNSigmaPosKa10(), row.daughterTPCNSigmaPosPr10(),
                     row.daughterTPCNSigmaNegPi10(), row.daughterTPCNSigmaNegKa10(), row.daughterTPCNSigmaNegPr10()};
  out.tofNSigma10 = {row.daughterTOFNSigmaPosPi10(), row.daughterTOFNSigmaPosKa10(), row.daughterTOFNSigmaPosPr10(),
                     row.daughterTOFNSigmaNegPi10(), row.daughterTOFNSigmaNegKa10(), row.daughterTOFNSigmaNegPr10()};
  out.crossedRows = {row.nCrossedRowsPos(), row.nCrossedRowsNeg()};
  const auto indices = row.indices();
  out.daughterIds = {indices[0], indices[1]};
  out.decayLength = detail::decayLength(collision, row);
  out.properLength = detail::properLength(out.decayLength, out.px, out.py, out.pz, o2::constants::physics::MassK0Short);
  return out;
}

template <typename MicroTrack>
TrackSnapshot makeTrackSnapshot(MicroTrack const& row)
{
  TrackSnapshot out;
  out.sourceTrackId = static_cast<int64_t>(row.trackId());
  out.charge = static_cast<int8_t>(row.sign());
  out.px = static_cast<float>(row.px());
  out.py = static_cast<float>(row.py());
  out.pz = static_cast<float>(row.pz());
  out.pidNSigmaPiFlag = static_cast<uint8_t>(row.pidNSigmaPiFlag());
  out.pidNSigmaKaFlag = static_cast<uint8_t>(row.pidNSigmaKaFlag());
  out.pidNSigmaPrFlag = static_cast<uint8_t>(row.pidNSigmaPrFlag());
  out.trackSelectionFlags = static_cast<uint8_t>(row.trackSelectionFlags());
  out.trackFlags = static_cast<uint8_t>(row.trackFlags());
  out.tpcNClsCrossedRows = static_cast<uint8_t>(row.tpcNClsCrossedRows());
  out.itsClusterMap = static_cast<uint8_t>(row.itsClusterMap());
  out.tpcNSigma = {static_cast<float>(row.tpcNSigmaPi()), static_cast<float>(row.tpcNSigmaKa()), static_cast<float>(row.tpcNSigmaPr())};
  out.tofNSigma = {static_cast<float>(row.tofNSigmaPi()), static_cast<float>(row.tofNSigmaKa()), static_cast<float>(row.tofNSigmaPr())};
  out.dcaXY = static_cast<float>(row.dcaXY());
  out.dcaZ = static_cast<float>(row.dcaZ());
  out.passedPtDependentDCAxy = row.passedPtDependentDCAxy();
  out.passedPtDependentDCAz = row.passedPtDependentDCAz();
  out.hasTOF = row.hasTOF();
  out.isPVContributor = row.isPVContributor();
  for (std::size_t layer = 0; layer < NItsLayers; ++layer) {
    out.itsHit[layer] = row.hasITSHitInLayer(static_cast<int>(layer));
  }
  return out;
}

inline bool hasFiniteMomentum(float px, float py, float pz)
{
  return std::isfinite(px) && std::isfinite(py) && std::isfinite(pz) && std::hypot(px, py) > 0.f;
}

inline bool sharesDaughter(CascadeSnapshot const& xi, int64_t trackId)
{
  return std::any_of(xi.daughterIds.begin(), xi.daughterIds.end(), [trackId](int id) { return static_cast<int64_t>(id) == trackId; });
}

/// Mode A canonical candidate [Xi, K0s]; the roles are fixed by the object types.
inline CanonicalizationResult<XiK0sCandidateSnapshot> canonicalizeXiK0s(CascadeSnapshot const& xi, V0Snapshot const& k0s)
{
  CanonicalizationResult<XiK0sCandidateSnapshot> result;
  for (const auto& id : k0s.daughterIds) {
    if (sharesDaughter(xi, id)) {
      result.status = BuildStatus::ReusedDaughter;
      return result;
    }
  }
  if (xi.sign != 1 && xi.sign != -1) {
    result.status = BuildStatus::InvalidChargePattern;
    return result;
  }
  if (!hasFiniteMomentum(xi.px, xi.py, xi.pz) || !hasFiniteMomentum(k0s.px, k0s.py, k0s.pz)) {
    result.status = BuildStatus::InvalidMomentum;
    return result;
  }
  result.candidate.xi = xi;
  result.candidate.k0s = k0s;
  result.status = BuildStatus::Ok;
  return result;
}

/// Charge pattern of (Xi, pion, kaon); 0 is the Omega(2012) signal pattern.
inline uint8_t chargePattern(int xiSign, int pionSign, int kaonSign)
{
  uint8_t pattern = kSignalPattern;
  if (pionSign == xiSign) {
    pattern |= kWrongSignPion;
  }
  if (kaonSign != xiSign) {
    pattern |= kWrongSignKaon;
  }
  return pattern;
}

/// Mode B canonical candidate [Xi, pion, kaon]; the roles are fixed by the selected species, the charge pattern is metadata.
inline CanonicalizationResult<Xi1530KCandidateSnapshot> canonicalizeXi1530K(CascadeSnapshot const& xi, TrackSnapshot const& pion, TrackSnapshot const& kaon)
{
  CanonicalizationResult<Xi1530KCandidateSnapshot> result;
  if (pion.sourceTrackId == kaon.sourceTrackId || sharesDaughter(xi, pion.sourceTrackId) || sharesDaughter(xi, kaon.sourceTrackId)) {
    result.status = BuildStatus::ReusedDaughter;
    return result;
  }
  if ((xi.sign != 1 && xi.sign != -1) || (pion.charge != 1 && pion.charge != -1) || (kaon.charge != 1 && kaon.charge != -1)) {
    result.status = BuildStatus::InvalidChargePattern;
    return result;
  }
  if (!hasFiniteMomentum(xi.px, xi.py, xi.pz) || !hasFiniteMomentum(pion.px, pion.py, pion.pz) || !hasFiniteMomentum(kaon.px, kaon.py, kaon.pz)) {
    result.status = BuildStatus::InvalidMomentum;
    return result;
  }
  result.candidate.xi = xi;
  result.candidate.pion = pion;
  result.candidate.kaon = kaon;
  result.candidate.chargePattern = chargePattern(xi.sign, pion.charge, kaon.charge);
  result.status = BuildStatus::Ok;
  return result;
}

namespace detail
{
struct EncodedValue {
  float value;
  float valid;
  float overflow;
};

// int8 x10 nSigma of ResoCascades/ResoV0s: saturated codes carry no usable value
inline EncodedValue encodeNSigma10(int8_t code)
{
  if (code == SaturatedLow || code == SaturatedHigh) {
    return {.value = 0.f, .valid = 0.f, .overflow = 0.f};
  }
  return {.value = static_cast<float>(code) / NSigmaScale, .valid = 1.f, .overflow = 0.f};
}

// micro001 nSigma, as the K1 contract
inline EncodedValue encodePID(float decoded)
{
  if (std::isnan(decoded)) {
    return {.value = 0.f, .valid = 0.f, .overflow = 0.f};
  }
  if (std::isinf(decoded)) {
    return {.value = std::signbit(decoded) ? -MicroPidOverflow : MicroPidOverflow, .valid = 1.f, .overflow = 1.f};
  }
  return {.value = decoded, .valid = 1.f, .overflow = 0.f};
}

// micro001 DCA, as the K1 contract
inline EncodedValue encodeDCA(float decoded)
{
  if (!std::isfinite(decoded)) {
    return {.value = 0.f, .valid = 0.f, .overflow = 0.f};
  }
  const bool overflow = decoded == o2::aod::resomicrodaughter001::DCAEncoding::MaxDCA;
  return {.value = decoded, .valid = 1.f, .overflow = overflow ? 1.f : 0.f};
}

template <std::size_t N>
class Writer
{
 public:
  explicit Writer(std::array<float, N>& output) : out(output) {}
  std::size_t index = 0;
  // A count beyond N is detected by finishPack(); nothing is written past the end.
  void add(float value)
  {
    if (index < N) {
      out[index] = value;
    }
    ++index;
  }
  void add(double value) { add(static_cast<float>(value)); }
  void add(bool value) { add(value ? 1.f : 0.f); }
  void addValid(EncodedValue value)
  {
    add(value.value);
    add(value.valid);
  }
  void addFull(EncodedValue value)
  {
    add(value.value);
    add(value.valid);
    add(value.overflow);
  }

 private:
  std::array<float, N>& out;
};

using LorentzVector = ROOT::Math::PxPyPzMVector;

inline LorentzVector fourMomentum(float px, float py, float pz, double mass)
{
  return {px, py, pz, mass};
}

template <std::size_t N>
void appendCascade(Writer<N>& w, CascadeSnapshot const& xi)
{
  w.add(xi.v0CosPA);
  w.add(xi.cascCosPA);
  w.add(xi.v0DaughDCA);
  w.add(xi.cascDaughDCA);
  w.add(std::abs(xi.dcaPosToPV));
  w.add(std::abs(xi.dcaNegToPV));
  w.add(std::abs(xi.dcaBachToPV));
  w.add(std::abs(xi.dcaV0ToPV));
  w.add(std::abs(xi.dcaXYToPV));
  w.add(std::abs(xi.dcaZToPV));
  w.add(xi.v0Radius);
  w.add(xi.cascRadius);
  w.add(static_cast<double>(xi.mLambda) - o2::constants::physics::MassLambda0);
  w.add(static_cast<double>(xi.mXi) - o2::constants::physics::MassXiMinus);
  w.add(xi.decayLength);
  w.add(xi.properLength);
  // Roles: the (anti)proton is the positive daughter of a Xi-, the negative daughter of a Xi+
  const bool isXiMinus = xi.sign < 0;
  constexpr std::size_t Pi = 0, Pr = 2, Pos = 0, Neg = 1, Bach = 2;
  const std::size_t baryon = isXiMinus ? Pos : Neg;
  const std::size_t meson = isXiMinus ? Neg : Pos;
  w.addValid(encodeNSigma10(xi.tpcNSigma10[baryon * NPidSpecies + Pr]));
  w.addValid(encodeNSigma10(xi.tpcNSigma10[meson * NPidSpecies + Pi]));
  w.addValid(encodeNSigma10(xi.tpcNSigma10[Bach * NPidSpecies + Pi]));
  w.addValid(encodeNSigma10(xi.tofNSigma10[baryon * NPidSpecies + Pr]));
  w.addValid(encodeNSigma10(xi.tofNSigma10[meson * NPidSpecies + Pi]));
  w.addValid(encodeNSigma10(xi.tofNSigma10[Bach * NPidSpecies + Pi]));
  w.add(static_cast<float>(xi.crossedRows[baryon]));
  w.add(static_cast<float>(xi.crossedRows[meson]));
  w.add(static_cast<float>(xi.crossedRows[Bach]));
}

template <std::size_t N>
void appendV0(Writer<N>& w, V0Snapshot const& k0s)
{
  w.add(k0s.cosPA);
  w.add(k0s.daughDCA);
  w.add(std::abs(k0s.dcaPosToPV));
  w.add(std::abs(k0s.dcaNegToPV));
  w.add(std::abs(k0s.dcaV0ToPV));
  w.add(k0s.radius);
  w.add(static_cast<double>(k0s.mK0Short) - o2::constants::physics::MassK0Short);
  w.add(static_cast<double>(k0s.mLambda) - o2::constants::physics::MassLambda0);
  w.add(static_cast<double>(k0s.mAntiLambda) - o2::constants::physics::MassLambda0);
  w.add(k0s.alpha);
  w.add(k0s.qtarm);
  w.add(k0s.decayLength);
  w.add(k0s.properLength);
  constexpr std::size_t Pi = 0, Pos = 0, Neg = 1;
  w.addValid(encodeNSigma10(k0s.tpcNSigma10[Pos * NPidSpecies + Pi]));
  w.addValid(encodeNSigma10(k0s.tpcNSigma10[Neg * NPidSpecies + Pi]));
  w.addValid(encodeNSigma10(k0s.tofNSigma10[Pos * NPidSpecies + Pi]));
  w.addValid(encodeNSigma10(k0s.tofNSigma10[Neg * NPidSpecies + Pi]));
  w.add(static_cast<float>(k0s.crossedRows[Pos]));
  w.add(static_cast<float>(k0s.crossedRows[Neg]));
}

template <std::size_t N>
void appendTrack(Writer<N>& w, TrackSnapshot const& track)
{
  for (const float& decoded : track.tpcNSigma) {
    w.addFull(encodePID(decoded));
  }
  for (const float& decoded : track.tofNSigma) {
    w.addFull(track.hasTOF ? encodePID(decoded) : EncodedValue{.value = 0.f, .valid = 0.f, .overflow = 0.f});
  }
  w.addFull(encodeDCA(track.dcaXY));
  w.addFull(encodeDCA(track.dcaZ));
  w.add(track.passedPtDependentDCAxy);
  w.add(track.passedPtDependentDCAz);
  w.add(track.hasTOF);
  w.add(static_cast<float>(track.tpcNClsCrossedRows));
  for (const bool& hit : track.itsHit) {
    w.add(hit);
  }
  w.add(track.isPVContributor);
}

// delta eta, sin/cos delta phi, delta R, z_pT of two objects
template <std::size_t N>
bool appendPair(Writer<N>& w, LorentzVector const& a, LorentzVector const& b)
{
  const double deta = a.Eta() - b.Eta();
  const double dphi = std::remainder(a.Phi() - b.Phi(), o2::constants::math::TwoPI);
  const double sumPt = a.Pt() + b.Pt();
  if (!std::isfinite(deta) || !std::isfinite(dphi) || !(sumPt > 0.)) {
    return false;
  }
  w.add(deta);
  w.add(std::sin(dphi));
  w.add(std::cos(dphi));
  w.add(std::hypot(deta, dphi));
  w.add(std::min(a.Pt(), b.Pt()) / sumPt);
  return true;
}

// Cosine of the daughter direction in the mother rest frame w.r.t. the mother direction
inline double cosThetaStar(LorentzVector const& daughter, LorentzVector const& mother)
{
  const double p = mother.P();
  if (!(p > 0.) || !(mother.E() > p)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const ROOT::Math::Boost boost{mother.BoostToCM()};
  const auto inRest = boost(daughter);
  const double pStar = inRest.P();
  if (!(pStar > 0.)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return inRest.Vect().Dot(mother.Vect()) / (pStar * p);
}

template <std::size_t N>
bool finishPack(FeaturePack<N>& pack, std::size_t index, LorentzVector const& mother, double scalarSumPt)
{
  const bool allFinite = std::all_of(pack.master.begin(), pack.master.end(), [](float x) { return std::isfinite(x); });
  if (index != N || !allFinite) {
    pack.status = BuildStatus::InvalidContract;
    return false;
  }
  pack.kinematics.mass = static_cast<float>(mother.M());
  pack.kinematics.pt = static_cast<float>(mother.Pt());
  pack.kinematics.eta = static_cast<float>(mother.Eta());
  pack.kinematics.phi = static_cast<float>(mother.Phi());
  pack.kinematics.scalarSumPt = static_cast<float>(scalarSumPt);
  pack.status = BuildStatus::Ok;
  return true;
}

template <std::size_t N>
constexpr bool isStrictlyIncreasingBelow(std::array<std::size_t, N> const& indices, std::size_t bound)
{
  for (std::size_t i = 0; i < N; ++i) {
    if (indices[i] >= bound || (i > 0 && indices[i - 1] >= indices[i])) {
      return false;
    }
  }
  return true;
}
} // namespace detail

// Names, projections and SHA-256 are generated from feature_contract_omega2012_v1.json (one section per mode).
inline constexpr std::array<std::string_view, XiK0sNMasterFeatures> XiK0sMasterFeatureNames{
  "xi.v0_cos_pa",
  "xi.casc_cos_pa",
  "xi.v0_daughter_dca",
  "xi.casc_daughter_dca",
  "xi.abs_dca_pos_to_pv",
  "xi.abs_dca_neg_to_pv",
  "xi.abs_dca_bach_to_pv",
  "xi.abs_dca_v0_to_pv",
  "xi.abs_dca_xy_to_pv",
  "xi.abs_dca_z_to_pv",
  "xi.v0_radius",
  "xi.casc_radius",
  "xi.delta_mass_lambda",
  "xi.delta_mass_xi",
  "xi.decay_length",
  "xi.proper_length",
  "xi.baryon_tpc_nsigma_pr",
  "xi.baryon_tpc_nsigma_pr_valid",
  "xi.meson_tpc_nsigma_pi",
  "xi.meson_tpc_nsigma_pi_valid",
  "xi.bach_tpc_nsigma_pi",
  "xi.bach_tpc_nsigma_pi_valid",
  "xi.baryon_tof_nsigma_pr",
  "xi.baryon_tof_nsigma_pr_valid",
  "xi.meson_tof_nsigma_pi",
  "xi.meson_tof_nsigma_pi_valid",
  "xi.bach_tof_nsigma_pi",
  "xi.bach_tof_nsigma_pi_valid",
  "xi.baryon_crossed_rows",
  "xi.meson_crossed_rows",
  "xi.bach_crossed_rows",
  "k0s.cos_pa",
  "k0s.daughter_dca",
  "k0s.abs_dca_pos_to_pv",
  "k0s.abs_dca_neg_to_pv",
  "k0s.abs_dca_v0_to_pv",
  "k0s.radius",
  "k0s.delta_mass_k0s",
  "k0s.delta_mass_lambda",
  "k0s.delta_mass_antilambda",
  "k0s.arm_alpha",
  "k0s.arm_qt",
  "k0s.decay_length",
  "k0s.proper_length",
  "k0s.pos_tpc_nsigma_pi",
  "k0s.pos_tpc_nsigma_pi_valid",
  "k0s.neg_tpc_nsigma_pi",
  "k0s.neg_tpc_nsigma_pi_valid",
  "k0s.pos_tof_nsigma_pi",
  "k0s.pos_tof_nsigma_pi_valid",
  "k0s.neg_tof_nsigma_pi",
  "k0s.neg_tof_nsigma_pi_valid",
  "k0s.pos_crossed_rows",
  "k0s.neg_crossed_rows",
  "xi.pt_fraction",
  "k0s.pt_fraction",
  "xi__k0s.delta_eta",
  "xi__k0s.sin_delta_phi",
  "xi__k0s.cos_delta_phi",
  "xi__k0s.delta_r",
  "xi__k0s.z_pt",
  "xi__k0s.opening_angle",
  "xi__k0s.cos_theta_star"};
inline constexpr std::array<std::size_t, 54> XiK0sDetectorV1Projection{
  0,
  1,
  2,
  3,
  4,
  5,
  6,
  7,
  8,
  9,
  10,
  11,
  12,
  13,
  14,
  15,
  16,
  17,
  18,
  19,
  20,
  21,
  22,
  23,
  24,
  25,
  26,
  27,
  28,
  29,
  30,
  31,
  32,
  33,
  34,
  35,
  36,
  37,
  38,
  39,
  40,
  41,
  42,
  43,
  44,
  45,
  46,
  47,
  48,
  49,
  50,
  51,
  52,
  53};
inline constexpr std::array<std::size_t, 63> XiK0sRelationalV1Projection{
  0,
  1,
  2,
  3,
  4,
  5,
  6,
  7,
  8,
  9,
  10,
  11,
  12,
  13,
  14,
  15,
  16,
  17,
  18,
  19,
  20,
  21,
  22,
  23,
  24,
  25,
  26,
  27,
  28,
  29,
  30,
  31,
  32,
  33,
  34,
  35,
  36,
  37,
  38,
  39,
  40,
  41,
  42,
  43,
  44,
  45,
  46,
  47,
  48,
  49,
  50,
  51,
  52,
  53,
  54,
  55,
  56,
  57,
  58,
  59,
  60,
  61,
  62};
inline constexpr std::string_view XiK0sFeatureContractSha256 = "b767a5643da563246295f543ec7296c9e21b00fd921fef88acc7b19f3b24942c";

inline constexpr std::array<std::string_view, Xi1530KNMasterFeatures> Xi1530KMasterFeatureNames{
  "xi.v0_cos_pa",
  "xi.casc_cos_pa",
  "xi.v0_daughter_dca",
  "xi.casc_daughter_dca",
  "xi.abs_dca_pos_to_pv",
  "xi.abs_dca_neg_to_pv",
  "xi.abs_dca_bach_to_pv",
  "xi.abs_dca_v0_to_pv",
  "xi.abs_dca_xy_to_pv",
  "xi.abs_dca_z_to_pv",
  "xi.v0_radius",
  "xi.casc_radius",
  "xi.delta_mass_lambda",
  "xi.delta_mass_xi",
  "xi.decay_length",
  "xi.proper_length",
  "xi.baryon_tpc_nsigma_pr",
  "xi.baryon_tpc_nsigma_pr_valid",
  "xi.meson_tpc_nsigma_pi",
  "xi.meson_tpc_nsigma_pi_valid",
  "xi.bach_tpc_nsigma_pi",
  "xi.bach_tpc_nsigma_pi_valid",
  "xi.baryon_tof_nsigma_pr",
  "xi.baryon_tof_nsigma_pr_valid",
  "xi.meson_tof_nsigma_pi",
  "xi.meson_tof_nsigma_pi_valid",
  "xi.bach_tof_nsigma_pi",
  "xi.bach_tof_nsigma_pi_valid",
  "xi.baryon_crossed_rows",
  "xi.meson_crossed_rows",
  "xi.bach_crossed_rows",
  "pion.tpc_nsigma_pi",
  "pion.tpc_nsigma_pi_valid",
  "pion.tpc_nsigma_pi_overflow",
  "pion.tpc_nsigma_ka",
  "pion.tpc_nsigma_ka_valid",
  "pion.tpc_nsigma_ka_overflow",
  "pion.tpc_nsigma_pr",
  "pion.tpc_nsigma_pr_valid",
  "pion.tpc_nsigma_pr_overflow",
  "pion.tof_nsigma_pi",
  "pion.tof_nsigma_pi_valid",
  "pion.tof_nsigma_pi_overflow",
  "pion.tof_nsigma_ka",
  "pion.tof_nsigma_ka_valid",
  "pion.tof_nsigma_ka_overflow",
  "pion.tof_nsigma_pr",
  "pion.tof_nsigma_pr_valid",
  "pion.tof_nsigma_pr_overflow",
  "pion.abs_dca_xy",
  "pion.abs_dca_xy_valid",
  "pion.abs_dca_xy_overflow",
  "pion.abs_dca_z",
  "pion.abs_dca_z_valid",
  "pion.abs_dca_z_overflow",
  "pion.passed_ptdep_dca_xy",
  "pion.passed_ptdep_dca_z",
  "pion.has_tof",
  "pion.tpc_crossed_rows",
  "pion.its_hit_l0",
  "pion.its_hit_l1",
  "pion.its_hit_l2",
  "pion.its_hit_l3",
  "pion.its_hit_l4",
  "pion.its_hit_l5",
  "pion.its_hit_l6",
  "pion.is_pv_contributor",
  "kaon.tpc_nsigma_pi",
  "kaon.tpc_nsigma_pi_valid",
  "kaon.tpc_nsigma_pi_overflow",
  "kaon.tpc_nsigma_ka",
  "kaon.tpc_nsigma_ka_valid",
  "kaon.tpc_nsigma_ka_overflow",
  "kaon.tpc_nsigma_pr",
  "kaon.tpc_nsigma_pr_valid",
  "kaon.tpc_nsigma_pr_overflow",
  "kaon.tof_nsigma_pi",
  "kaon.tof_nsigma_pi_valid",
  "kaon.tof_nsigma_pi_overflow",
  "kaon.tof_nsigma_ka",
  "kaon.tof_nsigma_ka_valid",
  "kaon.tof_nsigma_ka_overflow",
  "kaon.tof_nsigma_pr",
  "kaon.tof_nsigma_pr_valid",
  "kaon.tof_nsigma_pr_overflow",
  "kaon.abs_dca_xy",
  "kaon.abs_dca_xy_valid",
  "kaon.abs_dca_xy_overflow",
  "kaon.abs_dca_z",
  "kaon.abs_dca_z_valid",
  "kaon.abs_dca_z_overflow",
  "kaon.passed_ptdep_dca_xy",
  "kaon.passed_ptdep_dca_z",
  "kaon.has_tof",
  "kaon.tpc_crossed_rows",
  "kaon.its_hit_l0",
  "kaon.its_hit_l1",
  "kaon.its_hit_l2",
  "kaon.its_hit_l3",
  "kaon.its_hit_l4",
  "kaon.its_hit_l5",
  "kaon.its_hit_l6",
  "kaon.is_pv_contributor",
  "xi.pt_fraction",
  "pion.pt_fraction",
  "kaon.pt_fraction",
  "xi__pion.delta_eta",
  "xi__pion.sin_delta_phi",
  "xi__pion.cos_delta_phi",
  "xi__pion.delta_r",
  "xi__pion.z_pt",
  "xi__kaon.delta_eta",
  "xi__kaon.sin_delta_phi",
  "xi__kaon.cos_delta_phi",
  "xi__kaon.delta_r",
  "xi__kaon.z_pt",
  "pion__kaon.delta_eta",
  "pion__kaon.sin_delta_phi",
  "pion__kaon.cos_delta_phi",
  "pion__kaon.delta_r",
  "pion__kaon.z_pt",
  "xi1530__kaon.opening_angle",
  "xi1530__kaon.cos_theta_star",
  "mass_xi_pi",
  "mass_xi_k",
  "mass_pi_k"};
inline constexpr std::array<std::size_t, 103> Xi1530KDetectorV1Projection{
  0,
  1,
  2,
  3,
  4,
  5,
  6,
  7,
  8,
  9,
  10,
  11,
  12,
  13,
  14,
  15,
  16,
  17,
  18,
  19,
  20,
  21,
  22,
  23,
  24,
  25,
  26,
  27,
  28,
  29,
  30,
  31,
  32,
  33,
  34,
  35,
  36,
  37,
  38,
  39,
  40,
  41,
  42,
  43,
  44,
  45,
  46,
  47,
  48,
  49,
  50,
  51,
  52,
  53,
  54,
  55,
  56,
  57,
  58,
  59,
  60,
  61,
  62,
  63,
  64,
  65,
  66,
  67,
  68,
  69,
  70,
  71,
  72,
  73,
  74,
  75,
  76,
  77,
  78,
  79,
  80,
  81,
  82,
  83,
  84,
  85,
  86,
  87,
  88,
  89,
  90,
  91,
  92,
  93,
  94,
  95,
  96,
  97,
  98,
  99,
  100,
  101,
  102};
inline constexpr std::array<std::size_t, 124> Xi1530KRelationalV1Projection{
  0,
  1,
  2,
  3,
  4,
  5,
  6,
  7,
  8,
  9,
  10,
  11,
  12,
  13,
  14,
  15,
  16,
  17,
  18,
  19,
  20,
  21,
  22,
  23,
  24,
  25,
  26,
  27,
  28,
  29,
  30,
  31,
  32,
  33,
  34,
  35,
  36,
  37,
  38,
  39,
  40,
  41,
  42,
  43,
  44,
  45,
  46,
  47,
  48,
  49,
  50,
  51,
  52,
  53,
  54,
  55,
  56,
  57,
  58,
  59,
  60,
  61,
  62,
  63,
  64,
  65,
  66,
  67,
  68,
  69,
  70,
  71,
  72,
  73,
  74,
  75,
  76,
  77,
  78,
  79,
  80,
  81,
  82,
  83,
  84,
  85,
  86,
  87,
  88,
  89,
  90,
  91,
  92,
  93,
  94,
  95,
  96,
  97,
  98,
  99,
  100,
  101,
  102,
  103,
  104,
  105,
  106,
  107,
  108,
  109,
  110,
  111,
  112,
  113,
  114,
  115,
  116,
  117,
  118,
  119,
  120,
  121,
  122,
  123};
inline constexpr std::array<std::size_t, 126> Xi1530KStudyV1Projection{
  0,
  1,
  2,
  3,
  4,
  5,
  6,
  7,
  8,
  9,
  10,
  11,
  12,
  13,
  14,
  15,
  16,
  17,
  18,
  19,
  20,
  21,
  22,
  23,
  24,
  25,
  26,
  27,
  28,
  29,
  30,
  31,
  32,
  33,
  34,
  35,
  36,
  37,
  38,
  39,
  40,
  41,
  42,
  43,
  44,
  45,
  46,
  47,
  48,
  49,
  50,
  51,
  52,
  53,
  54,
  55,
  56,
  57,
  58,
  59,
  60,
  61,
  62,
  63,
  64,
  65,
  66,
  67,
  68,
  69,
  70,
  71,
  72,
  73,
  74,
  75,
  76,
  77,
  78,
  79,
  80,
  81,
  82,
  83,
  84,
  85,
  86,
  87,
  88,
  89,
  90,
  91,
  92,
  93,
  94,
  95,
  96,
  97,
  98,
  99,
  100,
  101,
  102,
  103,
  104,
  105,
  106,
  107,
  108,
  109,
  110,
  111,
  112,
  113,
  114,
  115,
  116,
  117,
  118,
  119,
  120,
  121,
  122,
  123,
  124,
  125};
inline constexpr std::string_view Xi1530KFeatureContractSha256 = "ca244fe167a599f9bc60d3909630edca5c797bb7b84d75d80944e6ad1a44aec5";
static_assert(XiK0sMasterFeatureNames.size() == XiK0sNMasterFeatures, "XiK0s feature name count must match XiK0sNMasterFeatures");
static_assert(Xi1530KMasterFeatureNames.size() == Xi1530KNMasterFeatures, "Xi1530K feature name count must match Xi1530KNMasterFeatures");
static_assert(detail::isStrictlyIncreasingBelow(XiK0sDetectorV1Projection, XiK0sNMasterFeatures), "XiK0s DetectorV1 projection must be strictly increasing");
static_assert(detail::isStrictlyIncreasingBelow(XiK0sRelationalV1Projection, XiK0sNMasterFeatures), "XiK0s RelationalV1 projection must be strictly increasing");
static_assert(detail::isStrictlyIncreasingBelow(Xi1530KDetectorV1Projection, Xi1530KNMasterFeatures), "Xi1530K DetectorV1 projection must be strictly increasing");
static_assert(detail::isStrictlyIncreasingBelow(Xi1530KRelationalV1Projection, Xi1530KNMasterFeatures), "Xi1530K RelationalV1 projection must be strictly increasing");
static_assert(detail::isStrictlyIncreasingBelow(Xi1530KStudyV1Projection, Xi1530KNMasterFeatures), "Xi1530K StudyV1 projection must be strictly increasing");

/// Mode A master features (XiK0s contract)
inline XiK0sFeaturePack buildXiK0sFeatures(XiK0sCandidateSnapshot const& candidate)
{
  using o2::constants::physics::MassK0Short;
  using o2::constants::physics::MassXiMinus;
  XiK0sFeaturePack pack;
  const auto& xi = candidate.xi;
  const auto& k0s = candidate.k0s;
  if (!hasFiniteMomentum(xi.px, xi.py, xi.pz) || !hasFiniteMomentum(k0s.px, k0s.py, k0s.pz)) {
    pack.status = BuildStatus::InvalidMomentum;
    return pack;
  }
  const auto pXi = detail::fourMomentum(xi.px, xi.py, xi.pz, MassXiMinus);
  const auto pK0s = detail::fourMomentum(k0s.px, k0s.py, k0s.pz, MassK0Short);
  const auto mother = pXi + pK0s;
  const double sumPt = pXi.Pt() + pK0s.Pt();
  if (!(sumPt > 0.)) {
    pack.status = BuildStatus::InvalidKinematics;
    return pack;
  }
  detail::Writer<XiK0sNMasterFeatures> w{pack.master};
  detail::appendCascade(w, xi);
  detail::appendV0(w, k0s);
  w.add(pXi.Pt() / sumPt);
  w.add(pK0s.Pt() / sumPt);
  if (!detail::appendPair(w, pXi, pK0s)) {
    pack.status = BuildStatus::InvalidKinematics;
    return pack;
  }
  w.add(ROOT::Math::VectorUtil::Angle(pXi, pK0s));
  w.add(detail::cosThetaStar(pXi, mother));
  detail::finishPack(pack, w.index, mother, sumPt);
  return pack;
}

/// Mode B master features (Xi1530K contract)
inline Xi1530KFeaturePack buildXi1530KFeatures(Xi1530KCandidateSnapshot const& candidate)
{
  using o2::constants::physics::MassKaonCharged;
  using o2::constants::physics::MassPionCharged;
  using o2::constants::physics::MassXiMinus;
  Xi1530KFeaturePack pack;
  const auto& xi = candidate.xi;
  const auto& pion = candidate.pion;
  const auto& kaon = candidate.kaon;
  if (!hasFiniteMomentum(xi.px, xi.py, xi.pz) || !hasFiniteMomentum(pion.px, pion.py, pion.pz) || !hasFiniteMomentum(kaon.px, kaon.py, kaon.pz)) {
    pack.status = BuildStatus::InvalidMomentum;
    return pack;
  }
  if (pion.sourceTrackId == kaon.sourceTrackId || sharesDaughter(xi, pion.sourceTrackId) || sharesDaughter(xi, kaon.sourceTrackId)) {
    pack.status = BuildStatus::ReusedDaughter;
    return pack;
  }
  const auto pXi = detail::fourMomentum(xi.px, xi.py, xi.pz, MassXiMinus);
  const auto pPion = detail::fourMomentum(pion.px, pion.py, pion.pz, MassPionCharged);
  const auto pKaon = detail::fourMomentum(kaon.px, kaon.py, kaon.pz, MassKaonCharged);
  const auto xi1530 = pXi + pPion;
  const auto mother = xi1530 + pKaon;
  const double sumPt = pXi.Pt() + pPion.Pt() + pKaon.Pt();
  if (!(sumPt > 0.)) {
    pack.status = BuildStatus::InvalidKinematics;
    return pack;
  }
  detail::Writer<Xi1530KNMasterFeatures> w{pack.master};
  detail::appendCascade(w, xi);
  detail::appendTrack(w, pion);
  detail::appendTrack(w, kaon);
  w.add(pXi.Pt() / sumPt);
  w.add(pPion.Pt() / sumPt);
  w.add(pKaon.Pt() / sumPt);
  if (!detail::appendPair(w, pXi, pPion) || !detail::appendPair(w, pXi, pKaon) || !detail::appendPair(w, pPion, pKaon)) {
    pack.status = BuildStatus::InvalidKinematics;
    return pack;
  }
  w.add(ROOT::Math::VectorUtil::Angle(xi1530, pKaon));
  w.add(detail::cosThetaStar(xi1530, mother));
  w.add(xi1530.M());
  w.add((pXi + pKaon).M());
  w.add((pPion + pKaon).M());
  detail::finishPack(pack, w.index, mother, sumPt);
  return pack;
}

template <std::size_t N, std::size_t M>
std::vector<float> project(FeaturePack<N> const& pack, std::array<std::size_t, M> const& indices)
{
  std::vector<float> projected;
  projected.reserve(M);
  for (const auto& i : indices) {
    projected.push_back(pack.master[i]);
  }
  return projected;
}

inline std::vector<float> projectFeatures(XiK0sFeaturePack const& pack, Profile profile)
{
  if (pack.status != BuildStatus::Ok) {
    return {};
  }
  switch (profile) {
    case Profile::DetectorV1:
      return project(pack, XiK0sDetectorV1Projection);
    case Profile::RelationalV1:
      return project(pack, XiK0sRelationalV1Projection);
    default:
      return {}; // mode A has no study-only features
  }
}

inline std::vector<float> projectFeatures(Xi1530KFeaturePack const& pack, Profile profile)
{
  if (pack.status != BuildStatus::Ok) {
    return {};
  }
  switch (profile) {
    case Profile::DetectorV1:
      return project(pack, Xi1530KDetectorV1Projection);
    case Profile::RelationalV1:
      return project(pack, Xi1530KRelationalV1Projection);
    case Profile::StudyV1:
      return project(pack, Xi1530KStudyV1Projection);
    default:
      return {};
  }
}
} // namespace o2::analysis::omega2012ml

#endif // PWGLF_CORE_OMEGA2012MLFEATURES_H_
