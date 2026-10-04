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
/// \file K1MlFeatures.h
/// \brief Canonical micro001 K1 candidate and feature contract helpers
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>
///

#ifndef PWGLF_CORE_K1MLFEATURES_H_
#define PWGLF_CORE_K1MLFEATURES_H_

#include "PWGLF/DataModel/LFResonanceTables.h"

#include <CommonConstants/PhysicsConstants.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string_view>
#include <vector>

namespace o2::analysis::k1ml
{
inline constexpr std::size_t NItsLayers = 7;
inline constexpr std::size_t NCandidateTracks = 3;
inline constexpr std::size_t NMasterFeatures = 125;

enum class Role : uint8_t { Kaon,
                            PionSame,
                            PionOpp };
enum class Profile : uint8_t { DetectorV1,
                               RelationalV1,
                               SubstructureV1 };
enum class BuildStatus : uint8_t {
  Ok,
  InvalidChargePattern,
  ReusedTrack,
  InvalidMomentum,
  InvalidKinematics,
  InvalidContract
};

struct TrackSnapshot {
  int64_t sourceTrackId = -1;
  int8_t charge = 0;
  float px = 0.f, py = 0.f, pz = 0.f;
  uint8_t pidNSigmaPiFlag = 0, pidNSigmaKaFlag = 0, pidNSigmaPrFlag = 0;
  uint8_t trackSelectionFlags = 0, trackFlags = 0, tpcNClsCrossedRows = 0, itsClusterMap = 0;
  std::array<float, 3> tpcNSigma{}; // pi, ka, pr
  std::array<float, 3> tofNSigma{}; // pi, ka, pr
  float dcaXY = std::numeric_limits<float>::quiet_NaN();
  float dcaZ = std::numeric_limits<float>::quiet_NaN();
  bool passedPtDependentDCAxy = false, passedPtDependentDCAz = false;
  bool hasTOF = false, isPVContributor = false;
  std::array<bool, NItsLayers> itsHit{};
};

struct CandidateSnapshot {
  std::array<TrackSnapshot, NCandidateTracks> tracks{}; // K, same-charge pion, opposite-charge pion
};

struct CanonicalizationResult {
  CandidateSnapshot candidate{};
  BuildStatus status = BuildStatus::InvalidChargePattern;
  explicit operator bool() const { return status == BuildStatus::Ok; }
};

struct KinematicAudit {
  float massKPiPi = 0.f;
  float massPiPi = 0.f;
  float massKaonPiSame = 0.f;
  float massKaonPiOpp = 0.f;
  float candidatePt = 0.f;
  float candidateEta = 0.f;
  float candidatePhi = 0.f;
  float scalarSumPt = 0.f;
};

struct FeaturePack {
  std::array<float, NMasterFeatures> master{};
  BuildStatus status = BuildStatus::InvalidContract;
  KinematicAudit kinematics{};
};

// Adapt the real ResoMicroTracks_001 public columns and dynamic accessors.
// Keeping the raw bytes alongside decoded accessors makes the encoding auditable.
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

inline bool hasFiniteMomentum(TrackSnapshot const& track)
{
  return std::isfinite(track.px) && std::isfinite(track.py) && std::isfinite(track.pz) &&
         std::hypot(track.px, track.py) > 0.f;
}

inline CanonicalizationResult canonicalizeUS(TrackSnapshot const& kaon,
                                             TrackSnapshot const& pionA,
                                             TrackSnapshot const& pionB)
{
  CanonicalizationResult result;
  if (kaon.sourceTrackId == pionA.sourceTrackId || kaon.sourceTrackId == pionB.sourceTrackId ||
      pionA.sourceTrackId == pionB.sourceTrackId) {
    result.status = BuildStatus::ReusedTrack;
    return result;
  }
  if ((kaon.charge != 1 && kaon.charge != -1) || (pionA.charge != 1 && pionA.charge != -1) ||
      (pionB.charge != 1 && pionB.charge != -1) || pionA.charge == pionB.charge ||
      (pionA.charge != kaon.charge && pionB.charge != kaon.charge)) {
    result.status = BuildStatus::InvalidChargePattern;
    return result;
  }
  if (!hasFiniteMomentum(kaon) || !hasFiniteMomentum(pionA) || !hasFiniteMomentum(pionB)) {
    result.status = BuildStatus::InvalidMomentum;
    return result;
  }
  const auto& same = pionA.charge == kaon.charge ? pionA : pionB;
  const auto& opposite = pionA.charge == kaon.charge ? pionB : pionA;
  result.candidate.tracks = {kaon, same, opposite};
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

inline EncodedValue encodePID(float decoded)
{
  if (std::isnan(decoded)) {
    return {0.f, 0.f, 0.f};
  }
  if (std::isinf(decoded)) {
    return {std::signbit(decoded) ? -3.5f : 3.5f, 1.f, 1.f};
  }
  return {decoded, 1.f, 0.f};
}

inline EncodedValue encodeDCA(float decoded)
{
  if (!std::isfinite(decoded)) {
    return {0.f, 0.f, 0.f};
  }
  const bool overflow = decoded == o2::aod::resomicrodaughter001::DCAEncoding::MaxDCA;
  return {decoded, 1.f, overflow ? 1.f : 0.f};
}

struct Kinematics {
  double pt;
  double eta;
  double phi;
  double energy;
};

inline bool getKinematics(TrackSnapshot const& track, double mass, Kinematics& out)
{
  if (!hasFiniteMomentum(track)) {
    return false;
  }
  const double px = track.px, py = track.py, pz = track.pz;
  out.pt = std::hypot(px, py);
  const double p = std::hypot(out.pt, pz);
  out.eta = std::asinh(pz / out.pt);
  out.phi = std::atan2(py, px);
  out.energy = std::sqrt(p * p + mass * mass);
  return std::isfinite(out.pt) && std::isfinite(out.eta) && std::isfinite(out.phi) && std::isfinite(out.energy);
}

inline float invariantMass(Kinematics const& a, Kinematics const& b)
{
  const auto pxa = a.pt * std::cos(a.phi), pya = a.pt * std::sin(a.phi);
  const auto pxb = b.pt * std::cos(b.phi), pyb = b.pt * std::sin(b.phi);
  const double e = a.energy + b.energy;
  const double px = pxa + pxb, py = pya + pyb;
  // pz is recovered from pt*sinh(eta), matching the stored three-momentum.
  const double pz = a.pt * std::sinh(a.eta) + b.pt * std::sinh(b.eta);
  const double m2 = e * e - px * px - py * py - pz * pz;
  return static_cast<float>(std::sqrt(std::max(0.0, m2)));
}

inline float invariantMass(std::array<Kinematics, NCandidateTracks> const& tracks,
                           std::array<std::size_t, NCandidateTracks> const& indices)
{
  double e = 0., px = 0., py = 0., pz = 0.;
  for (const auto& i : indices) {
    const auto& track = tracks[i];
    e += track.energy;
    px += track.pt * std::cos(track.phi);
    py += track.pt * std::sin(track.phi);
    pz += track.pt * std::sinh(track.eta);
  }
  return static_cast<float>(std::sqrt(std::max(0.0, e * e - px * px - py * py - pz * pz)));
}

inline void append(std::array<float, NMasterFeatures>& out, std::size_t& index, EncodedValue value)
{
  out[index++] = value.value;
  out[index++] = value.valid;
  out[index++] = value.overflow;
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

// Names and projection indices are generated from feature_contract_v1.json.
inline constexpr std::array<std::string_view, NMasterFeatures> MasterFeatureNames{
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
  "kaon.pt_fraction",
  "pion_same.tpc_nsigma_pi",
  "pion_same.tpc_nsigma_pi_valid",
  "pion_same.tpc_nsigma_pi_overflow",
  "pion_same.tpc_nsigma_ka",
  "pion_same.tpc_nsigma_ka_valid",
  "pion_same.tpc_nsigma_ka_overflow",
  "pion_same.tpc_nsigma_pr",
  "pion_same.tpc_nsigma_pr_valid",
  "pion_same.tpc_nsigma_pr_overflow",
  "pion_same.tof_nsigma_pi",
  "pion_same.tof_nsigma_pi_valid",
  "pion_same.tof_nsigma_pi_overflow",
  "pion_same.tof_nsigma_ka",
  "pion_same.tof_nsigma_ka_valid",
  "pion_same.tof_nsigma_ka_overflow",
  "pion_same.tof_nsigma_pr",
  "pion_same.tof_nsigma_pr_valid",
  "pion_same.tof_nsigma_pr_overflow",
  "pion_same.abs_dca_xy",
  "pion_same.abs_dca_xy_valid",
  "pion_same.abs_dca_xy_overflow",
  "pion_same.abs_dca_z",
  "pion_same.abs_dca_z_valid",
  "pion_same.abs_dca_z_overflow",
  "pion_same.passed_ptdep_dca_xy",
  "pion_same.passed_ptdep_dca_z",
  "pion_same.has_tof",
  "pion_same.tpc_crossed_rows",
  "pion_same.its_hit_l0",
  "pion_same.its_hit_l1",
  "pion_same.its_hit_l2",
  "pion_same.its_hit_l3",
  "pion_same.its_hit_l4",
  "pion_same.its_hit_l5",
  "pion_same.its_hit_l6",
  "pion_same.is_pv_contributor",
  "pion_same.pt_fraction",
  "pion_opp.tpc_nsigma_pi",
  "pion_opp.tpc_nsigma_pi_valid",
  "pion_opp.tpc_nsigma_pi_overflow",
  "pion_opp.tpc_nsigma_ka",
  "pion_opp.tpc_nsigma_ka_valid",
  "pion_opp.tpc_nsigma_ka_overflow",
  "pion_opp.tpc_nsigma_pr",
  "pion_opp.tpc_nsigma_pr_valid",
  "pion_opp.tpc_nsigma_pr_overflow",
  "pion_opp.tof_nsigma_pi",
  "pion_opp.tof_nsigma_pi_valid",
  "pion_opp.tof_nsigma_pi_overflow",
  "pion_opp.tof_nsigma_ka",
  "pion_opp.tof_nsigma_ka_valid",
  "pion_opp.tof_nsigma_ka_overflow",
  "pion_opp.tof_nsigma_pr",
  "pion_opp.tof_nsigma_pr_valid",
  "pion_opp.tof_nsigma_pr_overflow",
  "pion_opp.abs_dca_xy",
  "pion_opp.abs_dca_xy_valid",
  "pion_opp.abs_dca_xy_overflow",
  "pion_opp.abs_dca_z",
  "pion_opp.abs_dca_z_valid",
  "pion_opp.abs_dca_z_overflow",
  "pion_opp.passed_ptdep_dca_xy",
  "pion_opp.passed_ptdep_dca_z",
  "pion_opp.has_tof",
  "pion_opp.tpc_crossed_rows",
  "pion_opp.its_hit_l0",
  "pion_opp.its_hit_l1",
  "pion_opp.its_hit_l2",
  "pion_opp.its_hit_l3",
  "pion_opp.its_hit_l4",
  "pion_opp.its_hit_l5",
  "pion_opp.its_hit_l6",
  "pion_opp.is_pv_contributor",
  "pion_opp.pt_fraction",
  "kaon__pion_same.delta_eta",
  "kaon__pion_same.sin_delta_phi",
  "kaon__pion_same.cos_delta_phi",
  "kaon__pion_same.z_pt",
  "kaon__pion_opp.delta_eta",
  "kaon__pion_opp.sin_delta_phi",
  "kaon__pion_opp.cos_delta_phi",
  "kaon__pion_opp.z_pt",
  "pion_same__pion_opp.delta_eta",
  "pion_same__pion_opp.sin_delta_phi",
  "pion_same__pion_opp.cos_delta_phi",
  "pion_same__pion_opp.z_pt",
  "mass_pi_pi",
  "mass_kaon_pion_opp"};
inline constexpr std::array<std::size_t, 108> DetectorV1Projection{
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
  109};
inline constexpr std::array<std::size_t, 123> RelationalV1Projection{
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
  122};
inline constexpr std::array<std::size_t, NMasterFeatures> SubstructureV1Projection{
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
  124};
static_assert(MasterFeatureNames.size() == NMasterFeatures, "master feature name count must match NMasterFeatures");
static_assert(detail::isStrictlyIncreasingBelow(DetectorV1Projection, NMasterFeatures), "DetectorV1 projection indices must be strictly increasing and below NMasterFeatures");
static_assert(detail::isStrictlyIncreasingBelow(RelationalV1Projection, NMasterFeatures), "RelationalV1 projection indices must be strictly increasing and below NMasterFeatures");
static_assert(detail::isStrictlyIncreasingBelow(SubstructureV1Projection, NMasterFeatures), "SubstructureV1 projection indices must be strictly increasing and below NMasterFeatures");
inline constexpr std::string_view FeatureContractSha256 = "39f38ece001581d8ebf57392fad045759a63ba49f57c411b9b566d7c1cc58b8a";

inline FeaturePack buildMasterFeatures(CandidateSnapshot const& candidate)
{
  FeaturePack pack;
  for (auto const& track : candidate.tracks) {
    if ((track.charge != 1 && track.charge != -1) || !hasFiniteMomentum(track)) {
      pack.status = BuildStatus::InvalidMomentum;
      return pack;
    }
  }
  if (candidate.tracks[0].sourceTrackId == candidate.tracks[1].sourceTrackId ||
      candidate.tracks[0].sourceTrackId == candidate.tracks[2].sourceTrackId ||
      candidate.tracks[1].sourceTrackId == candidate.tracks[2].sourceTrackId ||
      candidate.tracks[0].charge != candidate.tracks[1].charge ||
      candidate.tracks[0].charge == candidate.tracks[2].charge) {
    pack.status = BuildStatus::InvalidChargePattern;
    return pack;
  }

  std::array<detail::Kinematics, NCandidateTracks> kin{};
  for (std::size_t i = 0; i < NCandidateTracks; ++i) {
    if (!detail::getKinematics(candidate.tracks[i], i == 0 ? o2::constants::physics::MassKaonCharged : o2::constants::physics::MassPionCharged, kin[i])) {
      pack.status = BuildStatus::InvalidKinematics;
      return pack;
    }
  }
  const double sumPt = kin[0].pt + kin[1].pt + kin[2].pt;
  if (!std::isfinite(sumPt) || sumPt <= 0.) {
    pack.status = BuildStatus::InvalidKinematics;
    return pack;
  }
  const double totalPx = static_cast<double>(candidate.tracks[0].px) + candidate.tracks[1].px + candidate.tracks[2].px;
  const double totalPy = static_cast<double>(candidate.tracks[0].py) + candidate.tracks[1].py + candidate.tracks[2].py;
  const double totalPz = static_cast<double>(candidate.tracks[0].pz) + candidate.tracks[1].pz + candidate.tracks[2].pz;
  const double candidatePt = std::hypot(totalPx, totalPy);
  if (!std::isfinite(candidatePt)) {
    pack.status = BuildStatus::InvalidKinematics;
    return pack;
  }
  pack.kinematics.scalarSumPt = static_cast<float>(sumPt);
  pack.kinematics.candidatePt = static_cast<float>(candidatePt);
  pack.kinematics.candidateEta = candidatePt > 0. ? static_cast<float>(std::asinh(totalPz / candidatePt)) : 0.f;
  pack.kinematics.candidatePhi = static_cast<float>(std::atan2(totalPy, totalPx));
  pack.kinematics.massKPiPi = detail::invariantMass(kin, std::array<std::size_t, NCandidateTracks>{0, 1, 2});
  pack.kinematics.massPiPi = detail::invariantMass(kin[1], kin[2]);
  pack.kinematics.massKaonPiSame = detail::invariantMass(kin[0], kin[1]);
  pack.kinematics.massKaonPiOpp = detail::invariantMass(kin[0], kin[2]);

  std::size_t index = 0;
  for (std::size_t i = 0; i < NCandidateTracks; ++i) {
    const auto& track = candidate.tracks[i];
    for (const float& decoded : track.tpcNSigma) {
      detail::append(pack.master, index, detail::encodePID(decoded));
    }
    for (const float& decoded : track.tofNSigma) {
      detail::append(pack.master, index, track.hasTOF ? detail::encodePID(decoded) : detail::EncodedValue{0.f, 0.f, 0.f});
    }
    detail::append(pack.master, index, detail::encodeDCA(track.dcaXY));
    detail::append(pack.master, index, detail::encodeDCA(track.dcaZ));
    pack.master[index++] = track.passedPtDependentDCAxy ? 1.f : 0.f;
    pack.master[index++] = track.passedPtDependentDCAz ? 1.f : 0.f;
    pack.master[index++] = track.hasTOF ? 1.f : 0.f;
    pack.master[index++] = static_cast<float>(track.tpcNClsCrossedRows);
    for (const bool& hit : track.itsHit) {
      pack.master[index++] = hit ? 1.f : 0.f;
    }
    pack.master[index++] = track.isPVContributor ? 1.f : 0.f;
    pack.master[index++] = static_cast<float>(kin[i].pt / sumPt);
  }

  constexpr std::array<std::array<std::size_t, 2>, 3> PairIndices{{{{0, 1}}, {{0, 2}}, {{1, 2}}}};
  for (auto const& pair : PairIndices) {
    const auto i = pair[0], j = pair[1];
    const double dphi = kin[i].phi - kin[j].phi;
    const double zDenominator = kin[i].pt + kin[j].pt;
    if (!std::isfinite(dphi) || zDenominator <= 0.) {
      pack.status = BuildStatus::InvalidKinematics;
      return pack;
    }
    pack.master[index++] = static_cast<float>(kin[i].eta - kin[j].eta);
    pack.master[index++] = static_cast<float>(std::sin(dphi));
    pack.master[index++] = static_cast<float>(std::cos(dphi));
    pack.master[index++] = static_cast<float>(std::min(kin[i].pt, kin[j].pt) / zDenominator);
  }
  pack.master[index++] = pack.kinematics.massPiPi;
  pack.master[index++] = pack.kinematics.massKaonPiOpp;
  const bool allFinite = std::all_of(pack.master.begin(), pack.master.end(), [](float x) { return std::isfinite(x); });
  if (index != pack.master.size() || !allFinite) {
    pack.status = BuildStatus::InvalidContract;
    return pack;
  }
  pack.status = BuildStatus::Ok;
  return pack;
}

inline std::vector<float> projectFeatures(FeaturePack const& pack, Profile profile)
{
  if (pack.status != BuildStatus::Ok) {
    return {};
  }
  std::vector<float> projected;
  switch (profile) {
    case Profile::DetectorV1:
      projected.reserve(DetectorV1Projection.size());
      for (const auto& i : DetectorV1Projection) {
        projected.push_back(pack.master[i]);
      }
      break;
    case Profile::RelationalV1:
      projected.reserve(RelationalV1Projection.size());
      for (const auto& i : RelationalV1Projection) {
        projected.push_back(pack.master[i]);
      }
      break;
    case Profile::SubstructureV1:
      projected.reserve(SubstructureV1Projection.size());
      for (const auto& i : SubstructureV1Projection) {
        projected.push_back(pack.master[i]);
      }
      break;
    default:
      return {};
  }
  return projected;
}
} // namespace o2::analysis::k1ml

#endif // PWGLF_CORE_K1MLFEATURES_H_
