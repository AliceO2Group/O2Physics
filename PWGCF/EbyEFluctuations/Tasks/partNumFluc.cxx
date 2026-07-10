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

/// \file partNumFluc.cxx
/// \brief Task for particle number fluctuation analysis
/// \author Fan Si <fsi@physi.uni-heidelberg.de>

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/StringHelpers.h>
#include <Framework/runDataProcessing.h>

#include <TFormula.h>
#include <TGraph.h>
#include <TH3.h>
#include <THnBase.h>
#include <TList.h>
#include <TPDGCode.h>
#include <TParticlePDG.h>
#include <TRandom.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <ranges>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#define C_CS(cs) /*NOLINT(cppcoreguidelines-macro-usage)*/                                      \
  []<std::size_t... indices>(std::index_sequence<indices...>) {                                 \
    static_assert(std::is_array_v<std::remove_cvref_t<decltype(cs)>> &&                         \
                  std::same_as<std::remove_extent_t<std::remove_cvref_t<decltype(cs)>>, char>); \
    return ConstStr<(cs)[indices]...>{};                                                        \
  }(std::make_index_sequence<sizeof(cs) - 1>{})
#define C_SV(sv) /*NOLINT(cppcoreguidelines-macro-usage)*/                            \
  []<std::size_t... indices>(std::index_sequence<indices...>) {                       \
    static_assert(std::same_as<std::remove_cvref_t<decltype(sv)>, std::string_view>); \
    return ConstStr<(sv)[indices]...>{};                                              \
  }(std::make_index_sequence<(sv).size()>{})

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
using JoinedCollisions = soa::Join<Collisions, EvSels, PVMults, FT0MultZeqs, CentFT0As, CentFT0Cs, CentFT0Ms>;
using JoinedTracks = soa::Join<Tracks, TracksExtra, TracksDCA, TrackSelection, pidTPCFullPi, pidTPCFullKa, pidTPCFullPr, pidTOFbeta, pidTOFFullPi, pidTOFFullKa, pidTOFFullPr>;
using JoinedCollisionsWithMc = soa::Join<McCollisionLabels, JoinedCollisions>;
using JoinedTracksWithMc = soa::Join<McTrackLabels, JoinedTracks>;
using JoinedMcCollisions = soa::Join<McCollisions, McCollsExtra, MultMCExtras>;

namespace mini_mc_collision
{
DECLARE_SOA_COLUMN(NMcParticlesP, nMcParticlesP, std::uint16_t);
DECLARE_SOA_COLUMN(NMcParticlesM, nMcParticlesM, std::uint16_t);
} // namespace mini_mc_collision

DECLARE_SOA_TABLE(MiniMcCollisions, "AOD", "MINIMCCOLLISION", soa::Index<>, mini_mc_collision::NMcParticlesP, mini_mc_collision::NMcParticlesM);
using MiniMcCollision = MiniMcCollisions::iterator;

namespace mini_collision
{
DECLARE_SOA_COLUMN(Vz, vz, std::int8_t);
DECLARE_SOA_COLUMN(Centrality, centrality, std::uint16_t);
DECLARE_SOA_COLUMN(NTracksP, nTracksP, std::uint16_t);
DECLARE_SOA_COLUMN(NTracksM, nTracksM, std::uint16_t);
} // namespace mini_collision

DECLARE_SOA_TABLE(MiniCollisions, "AOD", "MINICOLLISION", soa::Index<>, mini_collision::Vz, mini_collision::Centrality, mini_collision::NTracksP, mini_collision::NTracksM);
using MiniCollision = MiniCollisions::iterator;

namespace mini_mc_particle
{
DECLARE_SOA_INDEX_COLUMN(MiniMcCollision, miniMcCollision);
DECLARE_SOA_COLUMN(SignedEfficiency, signedEfficiency, std::int16_t);
} // namespace mini_mc_particle

DECLARE_SOA_TABLE(MiniMcParticles, "AOD", "MINIMCPARTICLE", soa::Index<>, mini_mc_particle::MiniMcCollisionId, mini_mc_particle::SignedEfficiency);
using MiniMcParticle = MiniMcParticles::iterator;

namespace mini_track
{
DECLARE_SOA_INDEX_COLUMN(MiniCollision, miniCollision);
DECLARE_SOA_COLUMN(SignedEfficiency, signedEfficiency, std::int16_t);
} // namespace mini_track

DECLARE_SOA_TABLE(MiniTracks, "AOD", "MINITRACK", soa::Index<>, mini_track::MiniCollisionId, mini_track::SignedEfficiency);
using MiniTrack = MiniTracks::iterator;
} // namespace o2::aod

namespace
{
namespace fluctuation_calculator_base
{
inline constexpr std::int8_t MaxOrder{8};
inline constexpr std::int32_t NExponentKeys{MaxOrder * (MaxOrder + 1) / 2};
inline constexpr std::array<std::array<std::int8_t, 2>, NExponentKeys> ExponentKeys{[] {
  std::array<std::array<std::int8_t, 2>, NExponentKeys> result{};
  std::int32_t index{};
  for (std::int32_t const& iExponent : std::views::iota(1, MaxOrder + 1)) {
    for (std::int32_t const& jExponent : std::views::iota(1, iExponent + 1)) {
      result[index++] = {static_cast<std::int8_t>(iExponent), static_cast<std::int8_t>(jExponent)};
    }
  }
  return result;
}()};
inline constexpr std::int32_t NOrderKeys{[] {
  std::array<std::int32_t, MaxOrder + 1> counts{1};
  for (std::array<std::int8_t, 2> const& exponentKey /*o2-linter: disable=const-ref-in-for-loop*/ : ExponentKeys) {
    const std::int32_t weight{exponentKey[0]};
    for (std::int32_t const& sum : std::views::iota(weight, MaxOrder + 1)) {
      counts[sum] += counts[sum - weight];
    }
  }
  std::int32_t result{};
  for (std::int32_t const& count : counts) {
    result += count;
  }
  return result;
}()};
inline constexpr std::array<std::array<std::int8_t, NExponentKeys>, NOrderKeys> OrderKeys{[] {
  std::array<std::array<std::int8_t, NExponentKeys>, NOrderKeys> result{};
  std::array<std::int8_t, NExponentKeys> current{};
  std::int32_t index{};
  const auto fillOrderKeys{[&](const auto& self, const std::int32_t position, const std::int32_t sum, const std::int32_t target) constexpr -> void {
    if (position == NExponentKeys) {
      if (sum == target) {
        result[index++] = current;
      }
      return;
    }

    const std::int32_t weight{ExponentKeys[position][0]};
    for (std::int32_t power{}; sum + power * weight <= target; ++power) {
      current[position] = power;
      self(self, position + 1, sum + power * weight, target);
    }
  }};
  for (std::int32_t const& target : std::views::iota(0, MaxOrder + 1)) {
    fillOrderKeys(fillOrderKeys, 0, 0, target);
  }
  return result;
}()};
} // namespace fluctuation_calculator_base

class FluctuationCalculatorTrack
{
 public:
  FluctuationCalculatorTrack() = default;
  virtual ~FluctuationCalculatorTrack() = default;

  [[nodiscard]] double getProduct(const std::int32_t orderKeyIndex, const double weight = 1.) const
  {
    double product{1.};
    for (std::int32_t const& iExponentKey : std::views::iota(0, static_cast<std::int32_t>(fluctuation_calculator_base::NExponentKeys))) {
      if (fluctuation_calculator_base::OrderKeys[orderKeyIndex][iExponentKey] != 0) {
        product *= std::pow(mQs[iExponentKey], fluctuation_calculator_base::OrderKeys[orderKeyIndex][iExponentKey]);
      }
    }
    return product * weight;
  }
  void init() { mQs.fill(0.); }
  void fill(const double charge, const double efficiency, const double weight = 1.)
  {
    for (std::int32_t const& iExponentKey : std::views::iota(0, static_cast<std::int32_t>(fluctuation_calculator_base::NExponentKeys))) {
      mQs[iExponentKey] += std::pow(charge, fluctuation_calculator_base::ExponentKeys[iExponentKey][0]) / std::pow(efficiency, fluctuation_calculator_base::ExponentKeys[iExponentKey][1]) * weight;
    }
  }

 protected:
  std::array<double, fluctuation_calculator_base::NExponentKeys> mQs{};
};

template <typename... Es>
concept IsValidEnum = ((std::is_enum_v<std::remove_cvref_t<Es>> && requires { std::remove_cvref_t<Es>::N; }) && ...);
template <auto... EValues>
concept IsValid = IsValidEnum<std::remove_cvref_t<decltype(EValues)>...> && ((EValues != std::remove_cvref_t<decltype(EValues)>::N) && ...);

template <typename E>
  requires IsValidEnum<E>
inline constexpr std::int32_t toI(const E e)
{
  return static_cast<std::int32_t>(e);
}

template <typename E>
  requires IsValidEnum<E>
inline constexpr std::int32_t NEs{toI(E::N)};

template <typename E>
  requires IsValidEnum<E>
struct EnumInfo;

enum class NameKind {
  Default = 0,
  Lower,
  Display,
  DisplayLower,
  N
};
template <typename E, NameKind NameKindValue = NameKind::Default>
  requires IsValidEnum<E>
inline constexpr std::string_view getName(const std::int32_t index)
{
  if constexpr (NameKindValue == NameKind::Lower) {
    return EnumInfo<std::remove_cvref_t<E>>::NamesLower.at(index);
  } else if constexpr (NameKindValue == NameKind::Display) {
    return EnumInfo<std::remove_cvref_t<E>>::DisplayNames.at(index);
  } else if constexpr (NameKindValue == NameKind::DisplayLower) {
    return EnumInfo<std::remove_cvref_t<E>>::DisplayNamesLower.at(index);
  } else {
    return EnumInfo<std::remove_cvref_t<E>>::Names.at(index);
  }
}
template <NameKind NameKindValue = NameKind::Default, typename E>
  requires IsValidEnum<E>
inline constexpr std::string_view getName(const E e)
{
  return getName<E, NameKindValue>(toI(e));
}
template <typename E>
  requires IsValidEnum<E>
inline std::vector<std::string> getDisplayNames()
{
  return std::vector<std::string>(EnumInfo<std::remove_cvref_t<E>>::DisplayNames.begin(), EnumInfo<std::remove_cvref_t<E>>::DisplayNames.end());
}

enum class DataMode {
  McMcParticle = 0,
  McTrack,
  RawTrack,
  N
};
enum class CentralityDefinition {
  Ft0a = 0,
  Ft0c,
  Ft0m,
  N
};
enum class DcaKind {
  Mean = 0,
  Sigma,
  N
};
enum class DcaAxis {
  Xy = 0,
  Z,
  N
};
enum class Detector {
  Tpc = 0,
  Tof,
  N
};
enum class PidStrategy {
  Tpc = 0,
  TpcTof,
  N
};
enum class PidStrategyAll {
  Tpc = 0,
  Tof,
  TpcTofSeparated,
  TpcTofCombined,
  N
};
enum class ParticleSpecies {
  Pion = 0,
  Kaon,
  Proton,
  N
};
enum class ParticleSpeciesAll {
  All = 0,
  Pion,
  Kaon,
  Proton,
  N
};
enum class ParticleNumber {
  Charge = 0,
  Kaon,
  Proton,
  N
};
enum class ChargeSpecies {
  Plus = 0,
  Minus,
  N
};
enum class ChargeNumber {
  Plus = 0,
  Minus,
  Total,
  Net,
  N
};

template <>
struct EnumInfo<DcaKind> {
  inline static constexpr std::array<std::string_view, NEs<DcaKind>> Names{"Mean", "Sigma"};
};
template <>
struct EnumInfo<DcaAxis> {
  inline static constexpr std::array<std::string_view, NEs<DcaAxis>> Names{"Xy", "Z"};
  inline static constexpr std::array<std::string_view, NEs<DcaAxis>> DisplayNames{"Xy", "Z"};
};
template <>
struct EnumInfo<Detector> {
  inline static constexpr std::array<std::string_view, NEs<Detector>> Names{"Tpc", "Tof"};
  inline static constexpr std::array<std::string_view, NEs<Detector>> NamesLower{"tpc", "tof"};
  inline static constexpr std::array<std::string_view, NEs<Detector>> DisplayNames{"TPC", "TOF"};
  inline static constexpr std::array<PidStrategyAll, NEs<Detector>> PidStrategyAllValues{PidStrategyAll::Tpc, PidStrategyAll::Tof};
};
template <>
struct EnumInfo<PidStrategy> {
  inline static constexpr std::array<std::string_view, NEs<PidStrategy>> Names{"Tpc", "TpcTof"};
  inline static constexpr std::array<std::string_view, NEs<PidStrategy>> NamesLower{"tpc", "tpcTof"};
  inline static constexpr std::array<std::string_view, NEs<PidStrategy>> DisplayNames{"TPC", "TPC+TOF"};
  inline static constexpr std::array<PidStrategyAll, NEs<PidStrategy>> PidStrategyAllValues{PidStrategyAll::Tpc, PidStrategyAll::TpcTofCombined};
};
template <>
struct EnumInfo<ParticleSpecies> {
  inline static constexpr std::array<std::string_view, NEs<ParticleSpecies>> Names{"Pi", "Ka", "Pr"};
  inline static constexpr std::array<std::string_view, NEs<ParticleSpecies>> DisplayNames{"Pion", "Kaon", "Proton"};
  inline static constexpr std::array<std::string_view, NEs<ParticleSpecies>> DisplayNamesLower{"pion", "kaon", "proton"};
  inline static constexpr std::array<ParticleSpeciesAll, NEs<ParticleSpecies>> ParticleSpeciesAllValues{ParticleSpeciesAll::Pion, ParticleSpeciesAll::Kaon, ParticleSpeciesAll::Proton};
};
template <>
struct EnumInfo<ParticleSpeciesAll> {
  inline static constexpr std::array<std::string_view, NEs<ParticleSpeciesAll>> Names{"", "Pi", "Ka", "Pr"};
  inline static constexpr std::array<std::string_view, NEs<ParticleSpeciesAll>> DisplayNames{"All", "Pion", "Kaon", "Proton"};
  inline static constexpr std::array<std::string_view, NEs<ParticleSpeciesAll>> DisplayNamesLower{"all", "pion", "kaon", "proton"};
  inline static constexpr std::array<ParticleSpecies, NEs<ParticleSpeciesAll>> ParticleSpeciesValues{ParticleSpecies::N, ParticleSpecies::Pion, ParticleSpecies::Kaon, ParticleSpecies::Proton};
  inline static constexpr std::array<std::array<std::int32_t, NEs<ChargeSpecies>>, NEs<ParticleSpeciesAll>> PdgCodes{{{0, 0}, {PDG_t::kPiPlus, PDG_t::kPiMinus}, {PDG_t::kKPlus, PDG_t::kKMinus}, {PDG_t::kProton, PDG_t::kProtonBar}}};
  inline static constexpr std::array<double, NEs<ParticleSpeciesAll>> Masses{0., constants::physics::MassPiPlus, constants::physics::MassKPlus, constants::physics::MassProton};
};
template <>
struct EnumInfo<ParticleNumber> {
  inline static constexpr std::array<std::string_view, NEs<ParticleNumber>> Names{"Ch", "Ka", "Pr"};
  inline static constexpr std::array<std::string_view, NEs<ParticleNumber>> DisplayNames{"Charge", "Kaon", "Proton"};
  inline static constexpr std::array<std::string_view, NEs<ParticleNumber>> DisplayNamesLower{"charge", "kaon", "proton"};
};
template <>
struct EnumInfo<ChargeSpecies> {
  inline static constexpr std::array<std::string_view, NEs<ChargeSpecies>> Names{"P", "M"};
  inline static constexpr std::array<std::string_view, NEs<ChargeSpecies>> NamesLower{"p", "m"};
};
template <>
struct EnumInfo<ChargeNumber> {
  inline static constexpr std::array<std::string_view, NEs<ChargeNumber>> Names{"P", "M", "T", "N"};
};

template <typename To, typename E>
  requires IsValidEnum<To, E>
inline constexpr To getValue(const E e)
{
  if constexpr (std::is_same_v<std::remove_cvref_t<To>, PidStrategyAll>) {
    return EnumInfo<E>::PidStrategyAllValues[toI(e)];
  } else if constexpr (std::is_same_v<std::remove_cvref_t<To>, ParticleSpecies>) {
    return EnumInfo<E>::ParticleSpeciesValues[toI(e)];
  } else if constexpr (std::is_same_v<std::remove_cvref_t<To>, ParticleSpeciesAll>) {
    return EnumInfo<E>::ParticleSpeciesAllValues[toI(e)];
  }
  return To::N;
}
inline constexpr std::int32_t getPdgCode(const ParticleSpeciesAll e, const ChargeSpecies chargeSpecies)
{
  return EnumInfo<ParticleSpeciesAll>::PdgCodes[toI(e)][toI(chargeSpecies)];
}
inline constexpr double getMass(const ParticleSpeciesAll e)
{
  return EnumInfo<ParticleSpeciesAll>::Masses[toI(e)];
}

template <typename T>
  requires std::is_arithmetic_v<T>
inline constexpr std::int32_t nEnabled(const Configurable<LabeledArray<T>>& cfg)
{
  const LabeledArray<T>& la{cfg.value};
  return std::count_if(la[0], la[0] + la.rows() * la.cols(), [](T x) { return x != T{}; });
}
template <typename T>
  requires std::is_arithmetic_v<T>
inline constexpr bool isEnabled(const Configurable<LabeledArray<T>>& cfg)
{
  return nEnabled(cfg) > 0;
}
} // namespace

struct PartNumFluc {
  struct HolderCcdb {
    [[maybe_unused]] inline static constexpr std::int32_t NDimensionsEfficiency{4};

    std::map<std::int32_t, std::pair<std::int32_t, std::int32_t>> runNumbersIndicesGroupIndices;
    std::vector<std::array<std::array<std::array<const TFormula*, NEs<ChargeSpecies>>, NEs<DcaAxis>>, NEs<DcaKind>>> fPtDca;
    std::vector<std::array<std::array<std::array<const TH3*, NEs<ChargeSpecies>>, NEs<ParticleSpecies>>, NEs<Detector>>> hCentralityPtEtaShiftNSigmaPid;
    std::vector<std::array<std::array<std::array<const THnBase*, NEs<ChargeSpecies>>, NEs<ParticleSpecies>>, NEs<PidStrategy>>> hVzCentralityPtEtaEfficiency;

    void clear()
    {
      runNumbersIndicesGroupIndices.clear();
      fPtDca.clear();
      hCentralityPtEtaShiftNSigmaPid.clear();
      hVzCentralityPtEtaEfficiency.clear();
    }
  } holderCcdb{};

  struct HolderMcEvent {
    std::int32_t runNumber{};
    std::int32_t runIndex{};
    std::int32_t runGroupIndex{};
    double vz{};
    std::array<std::array<std::int32_t, NEs<ChargeSpecies>>, NEs<ParticleNumber>> numbers{};
    std::array<std::array<std::int32_t, NEs<ChargeSpecies>>, NEs<ParticleNumber>> numbersEff{};

    void clear() { *this = {}; }
  } holderMcEvent{};

  struct HolderEvent {
    inline static constexpr std::pair<double, double> RangeCentrality{0., 100.};

    std::int32_t runNumber{};
    std::int32_t runIndex{};
    std::int32_t runGroupIndex{};
    double vz{};
    std::array<std::int32_t, NEs<ChargeSpecies>> nGlobalTracks{};
    std::array<std::int32_t, NEs<ChargeSpecies>> nPvContributors{};
    std::array<std::array<std::array<double, NEs<ChargeSpecies>>, NEs<DcaAxis>>, NEs<DcaKind>> dca{};
    std::array<std::int32_t, NEs<ChargeSpecies>> nTofBeta{};
    double centrality{};
    std::int32_t subgroupIndex{};
    std::array<std::array<std::int32_t, NEs<ChargeSpecies>>, NEs<ParticleSpecies>> numbers{};

    void clear() { *this = {}; }
    [[nodiscard]] std::int32_t getNGlobalTracks() const { return std::accumulate(nGlobalTracks.begin(), nGlobalTracks.end(), 0); }
    [[nodiscard]] std::int32_t getNPvContributors() const { return std::accumulate(nPvContributors.begin(), nPvContributors.end(), 0); }
    template <DcaKind DcaKindValue, DcaAxis DcaAxisValue>
      requires IsValid<DcaKindValue, DcaAxisValue>
    [[nodiscard]] double getDca() const
    {
      const std::int32_t sumNGlobalTracks{getNGlobalTracks()};
      if (sumNGlobalTracks == 0) {
        return 0.;
      }

      double sumDca{};
      if constexpr (DcaKindValue == DcaKind::Sigma) {
        const double meanDca{getDca<DcaKind::Mean, DcaAxisValue>()};
        for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
          sumDca += (std::pow(dca[toI(DcaKind::Sigma)][toI(DcaAxisValue)][iChargeSpecies], 2.) + std::pow(dca[toI(DcaKind::Mean)][toI(DcaAxisValue)][iChargeSpecies] - meanDca, 2.)) * nGlobalTracks[iChargeSpecies];
        }
        return std::sqrt(sumDca / sumNGlobalTracks);
      }

      for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
        sumDca += dca[toI(DcaKind::Mean)][toI(DcaAxisValue)][iChargeSpecies] * nGlobalTracks[iChargeSpecies];
      }
      return sumDca / sumNGlobalTracks;
    }
    [[nodiscard]] std::int32_t getNTofBeta() const { return std::accumulate(nTofBeta.begin(), nTofBeta.end(), 0); }
  } holderEvent{};

  struct HolderMcParticle {
    std::int32_t pdgCode{};
    std::int32_t charge{};
    double pt{};
    double eta{};

    void clear() { *this = {}; }
  } holderMcParticle{};

  struct HolderTrack {
    inline static constexpr double TruncationAbsNSigmaPid{999.};
    static constexpr double truncateNSigmaPid(const double value) { return std::abs(value) < TruncationAbsNSigmaPid ? value : -TruncationAbsNSigmaPid; }

    std::array<double, NEs<DcaAxis>> dca{};
    std::int32_t sign{};
    double p{};
    double pt{};
    double eta{};
    double phi{};
    double phiIu{};
    std::array<bool, NEs<Detector>> hasPid{};
    std::array<std::array<double, NEs<ParticleSpecies>>, NEs<PidStrategyAll>> nSigmaPid{[] {
      std::array<std::array<double, NEs<ParticleSpecies>>, NEs<PidStrategyAll>> a{};
      std::array<double, NEs<ParticleSpecies>> b{};
      b.fill(-TruncationAbsNSigmaPid);
      a.fill(b);
      return a;
    }()};

    void clear() { *this = {}; }
  } holderTrack{};

  struct HolderDerivedData {
    template <std::integral T>
    static constexpr T convertFloor(const double value)
    {
      return std::numeric_limits<T>::lowest() <= value && value <= std::numeric_limits<T>::max() ? static_cast<T>(std::floor(value)) : (std::is_signed_v<T> ? std::numeric_limits<T>::lowest() : std::numeric_limits<T>::max());
    }
    template <std::integral T>
    static constexpr T convertRound(const double value)
    {
      return std::numeric_limits<T>::lowest() <= value && value <= std::numeric_limits<T>::max() ? static_cast<T>(std::round(value)) : (std::is_signed_v<T> ? std::numeric_limits<T>::lowest() : std::numeric_limits<T>::max());
    }

    std::array<std::uint16_t, NEs<ChargeSpecies>> nMcParticles{};
    std::array<std::uint16_t, NEs<ChargeSpecies>> nTracks{};
    std::vector<aod::mini_mc_particle::SignedEfficiency::type> signedEfficienciesMcParticle{[] {std::vector<aod::mini_mc_particle::SignedEfficiency::type> v{}; v.reserve(256); return v; }()};
    std::vector<aod::mini_track::SignedEfficiency::type> signedEfficienciesTrack{[] {std::vector<aod::mini_track::SignedEfficiency::type> v{}; v.reserve(256); return v; }()};

    void clear()
    {
      nMcParticles = {};
      nTracks = {};
      signedEfficienciesMcParticle.clear();
      signedEfficienciesTrack.clear();
    }
  } holderDerivedData{};

  std::array<std::array<std::unique_ptr<FluctuationCalculatorTrack>, NEs<ChargeNumber>>, NEs<ParticleNumber>> fluctuationCalculatorTrack{};

  struct : ConfigurableGroup {
    Configurable<std::string> cfgCcdbUrl{"cfgCcdbUrl", "http://ccdb-test.cern.ch:8080", "Url of CCDB"};
    Configurable<std::string> cfgCcdbPath{"cfgCcdbPath", "Users/f/fasi/test", "Path in CCDB"};
    Configurable<std::int64_t> cfgCcdbTimestampLatest{"cfgCcdbTimestampLatest", std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()).time_since_epoch().count(), "Latest timestamp in CCDB"};
  } groupCcdb{};

  struct : ConfigurableGroup {
    Configurable<bool> cfgFlagQaRun{"cfgFlagQaRun", false, "Run QA flag"};
    Configurable<bool> cfgFlagQaEvent{"cfgFlagQaEvent", false, "Event QA flag"};
    Configurable<bool> cfgFlagQaCentrality{"cfgFlagQaCentrality", false, "Centrality QA flag"};
    Configurable<bool> cfgFlagQaTrack{"cfgFlagQaTrack", false, "Track QA flag"};
    Configurable<bool> cfgFlagQaDca{"cfgFlagQaDca", false, "DCA QA flag"};
    Configurable<LabeledArray<std::int32_t>> cfgFlagsQaAcceptance{"cfgFlagsQaAcceptance", {std::array<std::int32_t, NEs<ParticleSpeciesAll>>{0, 0, 0, 0}.data(), NEs<ParticleSpeciesAll>, getDisplayNames<ParticleSpeciesAll>()}, "Acceptance QA flags"};
    Configurable<LabeledArray<std::int32_t>> cfgFlagsQaPhi{"cfgFlagsQaPhi", {std::array<std::int32_t, NEs<ParticleSpeciesAll>>{0, 0, 0, 0}.data(), NEs<ParticleSpeciesAll>, getDisplayNames<ParticleSpeciesAll>()}, "Phi QA flags"};
    Configurable<LabeledArray<std::int32_t>> cfgFlagsQaPid{"cfgFlagsQaPid", {std::array<std::int32_t, NEs<ParticleSpeciesAll>>{0, 0, 0, 0}.data(), NEs<ParticleSpeciesAll>, getDisplayNames<ParticleSpeciesAll>()}, "PID QA flags"};
    Configurable<bool> cfgFlagQaMc{"cfgFlagQaMc", false, "MC QA flag"};
    Configurable<LabeledArray<std::int32_t>> cfgFlagsCalculationYield{"cfgFlagsCalculationYield", {std::array<std::int32_t, NEs<ParticleSpecies>>{0, 0, 0}.data(), NEs<ParticleSpecies>, getDisplayNames<ParticleSpecies>()}, "Yield calculation flags"};
    Configurable<LabeledArray<std::int32_t>> cfgFlagsCalculationPurity{"cfgFlagsCalculationPurity", {std::array<std::int32_t, NEs<ParticleSpecies>>{0, 0, 0}.data(), NEs<ParticleSpecies>, getDisplayNames<ParticleSpecies>()}, "Purity calculation flags"};
    Configurable<LabeledArray<std::int32_t>> cfgFlagsCalculationFractionPrimary{"cfgFlagsCalculationFractionPrimary", {std::array<std::int32_t, NEs<ParticleSpecies>>{0, 0, 0}.data(), NEs<ParticleSpecies>, getDisplayNames<ParticleSpecies>()}, "Primary fraction calculation flags"};
    Configurable<LabeledArray<std::int32_t>> cfgFlagsCalculationFluctuation{"cfgFlagsCalculationFluctuation", {std::array<std::int32_t, NEs<ParticleNumber>>{0, 0, 0}.data(), NEs<ParticleNumber>, getDisplayNames<ParticleNumber>()}, "Fluctuation calculation flags"};
  } groupAnalysis{};

  struct : ConfigurableGroup {
    Configurable<bool> cfgFlagRejectionRunBad{"cfgFlagRejectionRunBad", false, "Bad run rejection flag"};
    Configurable<bool> cfgFlagRejectionRunBadMc{"cfgFlagRejectionRunBadMc", false, "MC bad run rejection flag"};
    Configurable<std::string> cfgLabelFlagsRct{"cfgLabelFlagsRct", "CBT_hadronPID", "RCT flags label"};
    Configurable<LabeledArray<std::int32_t>> cfgFlagsRct{"cfgFlagsRct", {std::array<std::int32_t, 3>{0, 1, 1}.data(), 3, {"ZDC", "Acceptance", "Table"}}, "RCT flags"};
    Configurable<std::uint64_t> cfgBitsSelectionEvent{"cfgBitsSelectionEvent", std::uint64_t{0b10000000001101000000000000000000000000000000000000}, "Event selection bits"};
    Configurable<bool> cfgFlagInelEvent{"cfgFlagInelEvent", true, "Flag of requiring inelastic event"};
    Configurable<bool> cfgFlagInelEventMc{"cfgFlagInelEventMc", false, "Flag of requiring inelastic MC event"};
    Configurable<double> cfgCutMaxAbsVz{"cfgCutMaxAbsVz", 6., "Maximum absolute z-vertex position (cm)"};
    Configurable<bool> cfgFlagCutVzMc{"cfgFlagCutVzMc", false, "Flag of requiring MC z-vertex cut"};
    Configurable<std::int32_t> cfgCutMinDeviationNPvContributors{"cfgCutMinDeviationNPvContributors", -4, "Minimum nPvContributors deviation from nGlobalTracks"};
    Configurable<std::int32_t> cfgIndexDefinitionCentrality{"cfgIndexDefinitionCentrality", 2, "Centrality definition index"};
    ConfigurableAxis cfgAxisCentrality{"cfgAxisCentrality", {18, 0., 90.}, "Centrality axis in fluctuation calculation"};
    Configurable<std::int32_t> cfgNSubgroups{"cfgNSubgroups", 20, "Number of subgroups in fluctuation calculation"};
  } groupEvent{};

  struct : ConfigurableGroup {
    Configurable<bool> cfgFlagPvContributor{"cfgFlagPvContributor", true, "Flag of requiring PV contributor"};
    Configurable<std::int32_t> cfgCutMinItsNCls{"cfgCutMinItsNCls", 5, "Minimum number of clusters ITS"};
    Configurable<double> cfgCutMaxItsChi2NCls{"cfgCutMaxItsChi2NCls", 30., "Maximum chi2 per cluster ITS"};
    Configurable<std::int32_t> cfgCutMinTpcNCls{"cfgCutMinTpcNCls", 55, "Minimum number of clusters TPC"};
    Configurable<double> cfgCutMaxTpcChi2NCls{"cfgCutMaxTpcChi2NCls", 3.5, "Maximum chi2 per cluster TPC"};
    Configurable<double> cfgCutMaxTpcNClsSharedRatio{"cfgCutMaxTpcNClsSharedRatio", 0.25, "Maximum ratio of shared clusters over clusters TPC"};
    Configurable<std::int32_t> cfgCutMinTpcNCrossedRows{"cfgCutMinTpcNCrossedRows", 75, "Minimum number of crossed rows TPC"};
    Configurable<double> cfgCutMinTpcNCrossedRowsRatio{"cfgCutMinTpcNCrossedRowsRatio", 0.8, "Minimum ratio of crossed rows over findable clusters TPC"};
    Configurable<bool> cfgFlagRecalibrationDca{"cfgFlagRecalibrationDca", false, "DCA recalibration flag"};
    Configurable<LabeledArray<double>> cfgCutsMaxAbsNSigmaDca{"cfgCutsMaxAbsNSigmaDca", {std::array<double, NEs<DcaAxis>>{2.5, 2.5}.data(), NEs<DcaAxis>, getDisplayNames<DcaAxis>()}, "Maximum absolute nSigma values of DCA (cm)"};
    Configurable<double> cfgCutMinPt{"cfgCutMinPt", 0.4, "Minimum pT (GeV/c)"};
    Configurable<double> cfgCutMaxPt{"cfgCutMaxPt", 2., "Maximum pT (GeV/c)"};
    Configurable<double> cfgCutMaxAbsEta{"cfgCutMaxAbsEta", 0.8, "Maximum absolute eta"};
    Configurable<LabeledArray<double>> cfgThresholdsPtTofPid{"cfgThresholdsPtTofPid", {std::array<double, NEs<ParticleSpecies>>{0.5, 0.5, 0.8}.data(), NEs<ParticleSpecies>, getDisplayNames<ParticleSpecies>()}, "pT (GeV/c) thresholds for TOF PID"};
    Configurable<LabeledArray<std::int32_t>> cfgFlagsRecalibrationNSigmaPid{"cfgFlagsRecalibrationNSigmaPid", {std::array<std::int32_t, NEs<ParticleSpecies>>{0, 0, 0}.data(), NEs<ParticleSpecies>, getDisplayNames<ParticleSpecies>()}, "nSigma PID recalibration flags"};
    Configurable<bool> cfgFlagRejectionOthers{"cfgFlagRejectionOthers", false, "Other particle species rejection flag"};
    Configurable<LabeledArray<double>> cfgCutsMaxAbsNSigmaPid{"cfgCutsMaxAbsNSigmaPid", {std::array<double, NEs<ParticleSpecies>>{2., 2., 2.}.data(), NEs<ParticleSpecies>, getDisplayNames<ParticleSpecies>()}, "Maximum absolute nSigma values for PID"};
    Configurable<bool> cfgFlagMcParticlePhysicalPrimary{"cfgFlagMcParticlePhysicalPrimary", true, "Flag of requiring physical primary MC particle"};
    Configurable<bool> cfgFlagMcParticleMomentum{"cfgFlagMcParticleMomentum", true, "Flag of using momentum of MC particle"};
  } groupTrack{};

  HistogramRegistry hrCalculationFluctuation{"hrCalculationFluctuation", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry hrCalculationFractionPrimary{"hrCalculationFractionPrimary", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry hrCalculationPurity{"hrCalculationPurity", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry hrCalculationYield{"hrCalculationYield", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry hrQaMc{"hrQaMc", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry hrQaPid{"hrQaPid", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry hrQaPhi{"hrQaPhi", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry hrQaAcceptance{"hrQaAcceptance", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry hrQaDca{"hrQaDca", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry hrQaTrack{"hrQaTrack", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry hrQaCentrality{"hrQaCentrality", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry hrQaEvent{"hrQaEvent", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry hrQaRun{"hrQaRun", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry hrCounter{"hrCounter", {}, OutputObjHandlingPolicy::AnalysisObject};

  aod::rctsel::RCTFlagsChecker rctFlagsChecker;

  Service<framework::O2DatabasePDG> pdg{};
  Service<ccdb::BasicCCDBManager> ccdb{};

  Filter filterCollision{aod::evsel::sel8 == true};
  Filter filterTrack{requireQualityTracksInFilter() && requireTrackCutInFilter(TrackSelectionFlags::kGoldenChi2)};
  Filter filterMcCollision{aod::mccollisionprop::numRecoCollision > 0};

  Preslice<aod::JoinedTracksWithMc> presliceTracksPerCollision{aod::track::collisionId};

  Produces<aod::MiniMcCollisions> miniMcCollision{};
  Produces<aod::MiniCollisions> miniCollision{};
  Produces<aod::MiniMcParticles> miniMcParticle{};
  Produces<aod::MiniTracks> miniTrack{};

  void init(InitContext&)
  {
    gRandom->SetSeed(0);

    if (doprocessRaw.value == doprocessMc.value) {
      LOG(fatal) << "Identical doprocessRaw and doprocessMc!";
    }
    if (doprocessRaw.value) {
      LOG(info) << "Enabling raw data process.";
    } else if (doprocessMc.value) {
      LOG(info) << "Enabling MC data process.";
    }

    rctFlagsChecker.init(groupEvent.cfgLabelFlagsRct.value, static_cast<bool>(groupEvent.cfgFlagsRct.value.get("ZDC")), static_cast<bool>(groupEvent.cfgFlagsRct.value.get("Acceptance")), static_cast<bool>(groupEvent.cfgFlagsRct.value.get("Table")));

    ccdb->setURL(groupCcdb.cfgCcdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(true);
    if (groupCcdb.cfgCcdbTimestampLatest.value >= 0) {
      ccdb->setCreatedNotAfter(groupCcdb.cfgCcdbTimestampLatest.value);
    }
    const TList* const ccdbObject{ccdb->get<TList>(groupCcdb.cfgCcdbPath.value)};
    if (!ccdbObject || ccdbObject->IsA() != TList::Class()) {
      LOG(fatal) << "Invalid CCDB object!";
    }
    holderCcdb.clear();

    std::int32_t nRunsBad{};
    std::int32_t nRunGroups{};
    {
      const TGraph* const gRunNumberGroupIndex{dynamic_cast<const TGraph*>(ccdbObject->FindObject("gRunNumberGroupIndex"))};
      if (!gRunNumberGroupIndex || gRunNumberGroupIndex->IsA() != TGraph::Class()) {
        LOG(fatal) << "Invalid gRunNumberGroupIndex!";
      }
      for (std::int32_t const& iRun : std::views::iota(0, gRunNumberGroupIndex->GetN())) {
        const std::int32_t runGroupIndex{static_cast<std::int32_t>(std::llrint(gRunNumberGroupIndex->GetY()[iRun]))};
        if (runGroupIndex == 0 || (groupEvent.cfgFlagRejectionRunBad.value && runGroupIndex < 0)) {
          ++nRunsBad;
        }
        nRunGroups = std::max(nRunGroups, std::abs(runGroupIndex));
        holderCcdb.runNumbersIndicesGroupIndices[std::llrint(gRunNumberGroupIndex->GetX()[iRun])] = {iRun, runGroupIndex};
      }
    }
    if (groupEvent.cfgFlagRejectionRunBadMc.value) {
      const TGraph* const gRunNumberGroupIndex{dynamic_cast<const TGraph*>(ccdbObject->FindObject("gRunNumberGroupIndex_mc"))};
      if (!gRunNumberGroupIndex || gRunNumberGroupIndex->IsA() != TGraph::Class()) {
        LOG(fatal) << "Invalid gRunNumberGroupIndex_mc!";
      }
      for (std::int32_t const& iRun : std::views::iota(0, gRunNumberGroupIndex->GetN())) {
        if (std::llrint(gRunNumberGroupIndex->GetY()[iRun]) <= 0) {
          auto iter = holderCcdb.runNumbersIndicesGroupIndices.find(std::llrint(gRunNumberGroupIndex->GetX()[iRun]));
          if (iter != holderCcdb.runNumbersIndicesGroupIndices.end() && iter->second.second > 0) {
            iter->second.second = -iter->second.second;
            ++nRunsBad;
          }
        }
      }
    }

    if (holderCcdb.runNumbersIndicesGroupIndices.empty()) {
      LOG(info) << "No run process enabled.";
    } else {
      LOG(info) << "Number of runs: " << holderCcdb.runNumbersIndicesGroupIndices.size();
      if (nRunsBad <= 0) {
        LOG(info) << "No run rejection enabled.";
      } else {
        LOG(info) << "Number of bad runs: " << nRunsBad;
      }
      for (const auto& [runNumber, runIndexGroupIndex] : holderCcdb.runNumbersIndicesGroupIndices) {
        if (runIndexGroupIndex.second == 0 || (groupEvent.cfgFlagRejectionRunBad.value && runIndexGroupIndex.second < 0)) {
          LOG(info) << "Enabling rejecting run: " << runNumber << " (" << runIndexGroupIndex.second << ")";
        } else {
          LOG(info) << "Enabling processing run: " << runNumber << " (" << std::abs(runIndexGroupIndex.second) << ")";
        }
      }
    }

    if (groupEvent.cfgLabelFlagsRct.value.empty()) {
      LOG(info) << "No RCT flags label enabled.";
    } else {
      LOG(info) << "Enabling RCT flags label: " << groupEvent.cfgLabelFlagsRct.value;
    }
    if (static_cast<bool>(groupEvent.cfgFlagsRct.value.get("ZDC"))) {
      LOG(info) << "Enabling RCT flag: ZDC";
    }
    if (static_cast<bool>(groupEvent.cfgFlagsRct.value.get("Acceptance"))) {
      LOG(info) << "Enabling RCT flag: acceptance";
    }
    if (static_cast<bool>(groupEvent.cfgFlagsRct.value.get("Table"))) {
      LOG(info) << "Enabling RCT flag: table";
    }

    if ((groupEvent.cfgBitsSelectionEvent.value & ((std::uint64_t{1} << aod::evsel::EventSelectionFlags::kNsel) - 1)) == 0) {
      LOG(info) << "No event selection bit enabled.";
    } else {
      for (std::int32_t const& iEvSel : std::views::iota(0, aod::evsel::EventSelectionFlags::kNsel)) {
        if (static_cast<bool>((groupEvent.cfgBitsSelectionEvent.value >> iEvSel) & 1)) {
          LOG(info) << "Enabling event selection bit: " << aod::evsel::selectionLabels[iEvSel];
        }
      }
    }

    switch (groupEvent.cfgIndexDefinitionCentrality.value) {
      case toI(CentralityDefinition::Ft0a):
        LOG(info) << "Enabling centrality definition: FT0A";
        break;
      case toI(CentralityDefinition::Ft0c):
        LOG(info) << "Enabling centrality definition: FT0C";
        break;
      default:
        LOG(info) << "Enabling centrality definition: FT0M";
        break;
    }

    const auto readListRunGroup{[&](const std::int32_t runGroupIndex) -> const TList* {
      const char* const name{Form("lRunGroup_%d", runGroupIndex)};
      const TList* const lRunGroup{dynamic_cast<const TList*>(ccdbObject->FindObject(name))};
      if (!lRunGroup) {
        LOG(fatal) << "Invalid " << name << "!";
      }
      return lRunGroup;
    }};

    if (groupTrack.cfgFlagRecalibrationDca.value) {
      LOG(info) << "Enabling DCA recalibration.";

      holderCcdb.fPtDca.resize(nRunGroups);
      for (std::int32_t const& iRunGroup : std::views::iota(0, nRunGroups)) {
        const TList* const lRunGroup{readListRunGroup(iRunGroup + 1)};
        for (std::int32_t const& iDcaKind : std::views::iota(0, NEs<DcaKind>)) {
          for (std::int32_t const& iDcaAxis : std::views::iota(0, NEs<DcaAxis>)) {
            for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
              const char* const name{Form("fPt%sDca%s%s%s_runGroup%d", getName<DcaKind>(iDcaKind).data(), getName<DcaAxis>(iDcaAxis).data(), getName<ChargeSpecies>(iChargeSpecies).data(), doprocessMc.value ? "_mc" : "", iRunGroup + 1)};
              holderCcdb.fPtDca[iRunGroup][iDcaKind][iDcaAxis][iChargeSpecies] = dynamic_cast<const TFormula*>(lRunGroup->FindObject(name));
              if (!holderCcdb.fPtDca[iRunGroup][iDcaKind][iDcaAxis][iChargeSpecies]) {
                LOG(fatal) << "Invalid " << name << "!";
              }
              LOG(info) << "Reading from CCDB: " << name << " \"" << holderCcdb.fPtDca[iRunGroup][iDcaKind][iDcaAxis][iChargeSpecies]->GetExpFormula("clingp") << "\"";
            }
          }
        }
      }
    }

    for (std::int32_t const& iParticleSpecies : std::views::iota(0, NEs<ParticleSpecies>)) {
      if (!static_cast<bool>(groupTrack.cfgFlagsRecalibrationNSigmaPid.value.get(iParticleSpecies))) {
        continue;
      }

      LOG(info) << "Enabling nSigma" << getName<ParticleSpecies>(iParticleSpecies) << " recalibration.";

      holderCcdb.hCentralityPtEtaShiftNSigmaPid.resize(nRunGroups);
      for (std::int32_t const& iRunGroup : std::views::iota(0, nRunGroups)) {
        const TList* const lRunGroup{readListRunGroup(iRunGroup + 1)};
        for (std::int32_t const& iDetector : std::views::iota(0, NEs<Detector>)) {
          for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
            const char* const name{Form("hCentralityPtEtaShift%sNSigma%s%s%s_runGroup%d", getName<Detector>(iDetector).data(), getName<ParticleSpecies>(iParticleSpecies).data(), getName<ChargeSpecies>(iChargeSpecies).data(), doprocessMc.value ? "_mc" : "", iRunGroup + 1)};
            holderCcdb.hCentralityPtEtaShiftNSigmaPid[iRunGroup][iDetector][iParticleSpecies][iChargeSpecies] = dynamic_cast<const TH3*>(lRunGroup->FindObject(name));
            if (!holderCcdb.hCentralityPtEtaShiftNSigmaPid[iRunGroup][iDetector][iParticleSpecies][iChargeSpecies]) {
              LOG(fatal) << "Invalid " << name << "!";
            }
            LOG(info) << "Reading from CCDB: " << name;
          }
        }
      }
    }

    hrCounter.add("hNEvents", ";;No. of Events", {HistType::kTH1D, {{10 + aod::evsel::EventSelectionFlags::kNsel, -0.5, 9.5 + static_cast<double>(aod::evsel::EventSelectionFlags::kNsel), "Selection"}}});
    if (doprocessMc.value) {
      hrCounter.add("hNMcEvents", ";;No. of MC Events", {HistType::kTH1D, {{10, -0.5, 9.5, "Selection"}}});
    }

    if (groupAnalysis.cfgFlagQaRun.value) {
      LOG(info) << "Enabling run QA.";

      const HistogramConfigSpec hcsQaRun(HistType::kTProfile, {{static_cast<std::int32_t>(holderCcdb.runNumbersIndicesGroupIndices.size()), -0.5, holderCcdb.runNumbersIndicesGroupIndices.size() - 0.5, "Run Index"}});

      for (const auto& [name, title, isChargeSpeciesSeparated] : std::to_array<std::tuple<std::string, std::string, bool>>(
             {{"Vx", "#LT#it{V}_{#it{x}}#GT (cm)", false},
              {"Vy", "#LT#it{V}_{#it{y}}#GT (cm)", false},
              {"Vz", "#LT#it{V}_{#it{z}}#GT (cm)", false},
              {"MultiplicityFt0a", "FT0A #LTMultiplicity#GT", false},
              {"MultiplicityFt0c", "FT0C #LTMultiplicity#GT", false},
              {"CentralityFt0a", "FT0A #LTCentrality#GT", false},
              {"CentralityFt0c", "FT0C #LTCentrality#GT", false},
              {"CentralityFt0m", "FT0M #LTCentrality#GT", false},
              {"NGlobalTracks", "#LTnGlobalTracks#GT", true},
              {"NPvContributors", "#LTnPvContributors#GT", true},
              {Form("%sDca%s", getName(DcaKind::Mean).data(), getName(DcaAxis::Xy).data()), "#LT#LTDCA_{#it{xy}}#GT_{event}#GT (cm)", true},
              {Form("%sDca%s", getName(DcaKind::Sigma).data(), getName(DcaAxis::Xy).data()), "#LT#it{#sigma}(DCA_{#it{xy}})_{event}#GT (cm)", true},
              {Form("%sDca%s", getName(DcaKind::Mean).data(), getName(DcaAxis::Z).data()), "#LT#LTDCA_{#it{z}}#GT_{event}#GT (cm)", true},
              {Form("%sDca%s", getName(DcaKind::Sigma).data(), getName(DcaAxis::Z).data()), "#LT#it{#sigma}(DCA_{#it{z}})_{event}#GT (cm)", true},
              {"NTofBeta", "#LTnTofBeta#GT", true},
              {"ItsNCls", "ITS #LTnClusters#GT", true},
              {"ItsChi2NCls", "ITS #LT#it{#chi}^{2}/nClusters#GT", true},
              {"TpcNCls", "TPC #LTnClusters#GT", true},
              {"TpcChi2NCls", "TPC #LT#it{#chi}^{2}/nClusters#GT", true},
              {"TpcNClsSharedRatio", "TPC #LTnSharedClusters/nClusters#GT", true},
              {"TpcNCrossedRows", "TPC #LTnCrossedRows#GT", true},
              {"TpcNCrossedRowsRatio", "TPC #LTnCrossedRows/nFindableClusters#GT", true},
              {"DcaXy", "#LTDCA_{#it{xy}}#GT (cm)", true},
              {"DcaZ", "#LTDCA_{#it{z}}#GT (cm)", true},
              {"Pt", "#LT#it{p}_{T}#GT (GeV/#it{c})", true},
              {"Eta", "#LT#it{#eta}#GT", true},
              {"Phi", "#LT#it{#varphi}#GT", true},
              {"TpcDeDx", "TPC #LTd#it{E}/d#it{x}#GT (a.u.)", true},
              {Form("%sNSigma%s", getName(Detector::Tpc).data(), getName(ParticleSpecies::Pion).data()), Form("%s #LT#it{n}#it{#sigma}_{#pi}#GT", getName<NameKind::Display>(Detector::Tpc).data()), true},
              {Form("%sNSigma%s", getName(Detector::Tpc).data(), getName(ParticleSpecies::Kaon).data()), Form("%s #LT#it{n}#it{#sigma}_{K}#GT", getName<NameKind::Display>(Detector::Tpc).data()), true},
              {Form("%sNSigma%s", getName(Detector::Tpc).data(), getName(ParticleSpecies::Proton).data()), Form("%s #LT#it{n}#it{#sigma}_{p}#GT", getName<NameKind::Display>(Detector::Tpc).data()), true},
              {"TofInverseBeta", "TOF #LT1/#it{#beta}#GT", true},
              {Form("%sNSigma%s", getName(Detector::Tof).data(), getName(ParticleSpecies::Pion).data()), Form("%s #LT#it{n}#it{#sigma}_{#pi}#GT", getName<NameKind::Display>(Detector::Tof).data()), true},
              {Form("%sNSigma%s", getName(Detector::Tof).data(), getName(ParticleSpecies::Kaon).data()), Form("%s #LT#it{n}#it{#sigma}_{K}#GT", getName<NameKind::Display>(Detector::Tof).data()), true},
              {Form("%sNSigma%s", getName(Detector::Tof).data(), getName(ParticleSpecies::Proton).data()), Form("%s #LT#it{n}#it{#sigma}_{p}#GT", getName<NameKind::Display>(Detector::Tof).data()), true}})) {
        if (!isChargeSpeciesSeparated) {
          hrQaRun.add(Form("QaRun/pRunIndex%s", name.data()), Form(";;%s", title.data()), hcsQaRun);
        } else {
          hrQaRun.add(Form("QaRun/pRunIndex%s_%s", name.data(), getName<NameKind::Lower>(ChargeSpecies::Plus).data()), Form(";;%s (#it{q}>0)", title.data()), hcsQaRun);
          hrQaRun.add(Form("QaRun/pRunIndex%s_%s", name.data(), getName<NameKind::Lower>(ChargeSpecies::Minus).data()), Form(";;%s (#it{q}<0)", title.data()), hcsQaRun);
        }
      }
    }

    if (groupAnalysis.cfgFlagQaEvent.value) {
      LOG(info) << "Enabling event QA.";

      const AxisSpec asNTracks(200, -0.5, 199.5);
      const HistogramConfigSpec hcsQaEvent(HistType::kTHnSparseD, {asNTracks, asNTracks});

      hrQaEvent.add("QaEvent/hVxVy", "", {HistType::kTHnSparseD, {{150, -0.15, 0.15, "#it{V}_{#it{x}} (cm)"}, {150, -0.15, 0.15, "#it{V}_{#it{y}} (cm)"}}});
      hrQaEvent.add("QaEvent/hVz", "", {HistType::kTH1D, {{300, -15., 15., "#it{V}_{#it{z}} (cm)"}}});
      hrQaEvent.add("QaEvent/hNPvContributorsNGlobalTracks", ";nPvContributors;nGlobalTracks;", hcsQaEvent);
      hrQaEvent.add("QaEvent/hNGlobalTracksMeanDcaXy", ";nGlobalTracks;", {HistType::kTHnSparseD, {asNTracks, {250, -0.25, 0.25, "#LTDCA_{#it{xy}}#GT_{event} (cm)"}}});
      hrQaEvent.add("QaEvent/hNGlobalTracksMeanDcaXy_nPvContributorsCut", ";nGlobalTracks;", {HistType::kTHnSparseD, {asNTracks, {250, -0.25, 0.25, "#LTDCA_{#it{xy}}#GT_{event} (cm)"}}});
      hrQaEvent.add("QaEvent/hNGlobalTracksMeanDcaZ", ";nGlobalTracks;", {HistType::kTHnSparseD, {asNTracks, {200, -2., 2., "#LTDCA_{#it{z}}#GT_{event} (cm)"}}});
      hrQaEvent.add("QaEvent/hNGlobalTracksMeanDcaZ_nPvContributorsCut", ";nGlobalTracks;", {HistType::kTHnSparseD, {asNTracks, {200, -2., 2., "#LTDCA_{#it{z}}#GT_{event} (cm)"}}});
      hrQaEvent.add("QaEvent/hNTofBetaNGlobalTracks", ";nTofBeta;nGlobalTracks;", hcsQaEvent);
      hrQaEvent.add("QaEvent/hNTofBetaNGlobalTracks_nPvContributorsCut", ";nTofBeta;nGlobalTracks;", hcsQaEvent);
    }

    if (groupAnalysis.cfgFlagQaCentrality.value) {
      LOG(info) << "Enabling centrality QA.";

      hrQaCentrality.add("QaCentrality/hCentralitySelection", "", {HistType::kTHnSparseD, {{100, 0., 100., "Centrality (%)"}, {10 + aod::evsel::EventSelectionFlags::kNsel, -0.5, 9.5 + static_cast<double>(aod::evsel::EventSelectionFlags::kNsel), "Selection"}}});
      hrQaCentrality.add("QaCentrality/hCentralityMultiplicity", "", {HistType::kTHnSparseD, {{100, 0., 100., "Centrality (%)"}, {200, -0.5, 199.5, "Multiplicity"}}});
    }

    if (groupAnalysis.cfgFlagQaTrack.value) {
      LOG(info) << "Enabling track QA.";

      for (const auto& [name, configSpec] : std::to_array<std::pair<std::string_view, HistogramConfigSpec>>(
             {{"ItsNCls", {HistType::kTH1D, {{10, -0.5, 9.5, "ITS nClusters"}}}},
              {"ItsChi2NCls", {HistType::kTH1D, {{80, 0., 40., "ITS #it{#chi}^{2}/nClusters"}}}},
              {"TpcNClsNClsShared", {HistType::kTHnSparseD, {{180, -0.5, 179.5, "TPC nClusters"}, {180, -0.5, 179.5, "TPC nSharedClusters"}}}},
              {"TpcChi2NCls", {HistType::kTH1D, {{100, 0., 5., "TPC #it{#chi}^{2}/nClusters"}}}},
              {"TpcNClsFindableNCrossedRows", {HistType::kTHnSparseD, {{180, -0.5, 179.5, "TPC nFindableClusters"}, {180, -0.5, 179.5, "TPC nCrossedRows"}}}}})) {
        for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
          hrQaTrack.add(Form("QaTrack/h%s_%s", name.data(), getName<ChargeSpecies, NameKind::Lower>(iChargeSpecies).data()), "", configSpec);
        }
      }
    }

    if (groupAnalysis.cfgFlagQaDca.value) {
      LOG(info) << "Enabling DCA QA.";

      const AxisSpec asPt(40, 0., 2., "#it{p}_{T} (GeV/#it{c})");
      const HistogramConfigSpec hcsQaDcaProfile(HistType::kTProfile3D, {asPt, {24, -1.2, 1.2, "#it{#eta}"}, {constants::math::NSectors, 0., constants::math::TwoPI, "#it{#varphi} (rad)"}});

      for (const auto& [name, title, configSpec] : std::to_array<std::tuple<std::string_view, std::string_view, HistogramConfigSpec>>(
             {{"hPtDcaXy", "", {HistType::kTHnSparseD, {asPt, {250, -0.25, 0.25, "DCA_{#it{xy}} (cm)"}}}},
              {"pPtEtaPhiIuDcaXy", ";;#LTDCA_{#it{xy}}#GT (cm)", hcsQaDcaProfile},
              {"hPtDcaZ", "", {HistType::kTHnSparseD, {asPt, {250, -0.5, 0.5, "DCA_{#it{z}} (cm)"}}}},
              {"pPtEtaPhiIuDcaZ", ";;#LTDCA_{#it{z}}#GT (cm)", hcsQaDcaProfile}})) {
        for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
          hrQaDca.add(Form("QaDca/%s_%s", name.data(), getName<ChargeSpecies, NameKind::Lower>(iChargeSpecies).data()), title.data(), configSpec);
        }
      }
    }

    for (std::int32_t const& iParticleSpeciesAll : std::views::iota(0, NEs<ParticleSpeciesAll>)) {
      if (!static_cast<bool>(groupAnalysis.cfgFlagsQaAcceptance.value.get(iParticleSpeciesAll))) {
        continue;
      }

      LOG(info) << "Enabling " << getName<ParticleSpeciesAll, NameKind::DisplayLower>(iParticleSpeciesAll) << " acceptance QA.";

      for (std::int32_t const& iPidStrategy : std::views::iota(0, NEs<PidStrategy>)) {
        for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
          hrQaAcceptance.add(Form("QaAcceptance/h%sPt_%sEdge%s%s", iParticleSpeciesAll == toI(ParticleSpeciesAll::All) ? "Eta" : "Rapidity", getName<PidStrategy, NameKind::Lower>(iPidStrategy).data(), getName<ParticleSpeciesAll>(iParticleSpeciesAll).data(), getName<ChargeSpecies>(iChargeSpecies).data()), "", {HistType::kTHnSparseD, {{300, -1.5, 1.5, iParticleSpeciesAll == toI(ParticleSpeciesAll::All) ? "#it{#eta}" : "#it{y}"}, {250, 0., 2.5, "#it{p}_{T} (GeV/#it{c})"}}});
        }
      }
    }

    for (std::int32_t const& iParticleSpeciesAll : std::views::iota(0, NEs<ParticleSpeciesAll>)) {
      if (!static_cast<bool>(groupAnalysis.cfgFlagsQaPhi.value.get(iParticleSpeciesAll))) {
        continue;
      }

      LOG(info) << "Enabling " << getName<ParticleSpeciesAll, NameKind::DisplayLower>(iParticleSpeciesAll) << " phi QA.";

      const HistogramConfigSpec hcsQaPhi(HistType::kTHnSparseF, {{{0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.}, "Centrality (%)"}, {20, 0., 2., "#it{p}_{T} (GeV/#it{c})"}, {24, -1.2, 1.2, "#it{#eta}"}, {360, 0., constants::math::TwoPI, "#it{#varphi} (rad)"}});

      for (std::int32_t const& iPidStrategy : std::views::iota(0, NEs<PidStrategy>)) {
        for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
          hrQaPhi.add(Form("QaPhi/hCentralityPtEtaPhi_%s%s%s", doprocessMc.value ? Form("mc%s", getName<PidStrategy>(iPidStrategy).data()) : getName<PidStrategy, NameKind::Lower>(iPidStrategy).data(), getName<ParticleSpeciesAll>(iParticleSpeciesAll).data(), getName<ChargeSpecies>(iChargeSpecies).data()), "", hcsQaPhi);
          hrQaPhi.add(Form("QaPhi/hCentralityPtEtaPhiIu_%s%s%s", doprocessMc.value ? Form("mc%s", getName<PidStrategy>(iPidStrategy).data()) : getName<PidStrategy, NameKind::Lower>(iPidStrategy).data(), getName<ParticleSpeciesAll>(iParticleSpeciesAll).data(), getName<ChargeSpecies>(iChargeSpecies).data()), "", hcsQaPhi);
        }
      }
    }

    for (std::int32_t const& iParticleSpeciesAll : std::views::iota(0, NEs<ParticleSpeciesAll>)) {
      if (!static_cast<bool>(groupAnalysis.cfgFlagsQaPid.value.get(iParticleSpeciesAll))) {
        continue;
      }

      LOG(info) << "Enabling " << getName<ParticleSpeciesAll, NameKind::DisplayLower>(iParticleSpeciesAll) << " PID QA.";

      const AxisSpec asCentrality({0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.}, "Centrality (%)");

      if (iParticleSpeciesAll == toI(ParticleSpeciesAll::All)) {
        const AxisSpec asPOverQ(350, -3.5, 3.5, "#it{p}/#it{q} (GeV/#it{c})");
        const AxisSpec asEta(48, -1.2, 1.2, "#it{#eta}");

        hrQaPid.add("QaPid/hCentralityPOverQEtaTpcLnDeDx", "", {HistType::kTHnSparseF, {asCentrality, asPOverQ, asEta, {240, 3., 9., "TPC ln(d#it{E}/d#it{x} (a.u.))"}}});
        hrQaPid.add("QaPid/hCentralityPOverQEtaTofInverseBeta", "", {HistType::kTHnSparseF, {asCentrality, asPOverQ, asEta, {120, 0.5, 3.5, "TOF 1/#it{#beta}"}}});
      } else {
        const HistogramConfigSpec hcsQaPid(HistType::kTHnSparseF, {asCentrality, {40, 0., 2., "#it{p}_{T} (GeV/#it{c})"}, {32, -0.8, 0.8, "#it{#eta}"}, {300, -30., 30.}});

        constexpr std::array<std::string_view, NEs<ParticleSpeciesAll>> ParticleSpeciesAllTitles{"", "#pi", "K", "p"};
        if (doprocessMc.value) {
          for (std::int32_t const& iDetector : std::views::iota(0, NEs<Detector>)) {
            for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
              hrQaPid.add(Form("QaPid/hCentralityPtEta%sNSigma%s_mc%s%s", getName<Detector>(iDetector).data(), getName<ParticleSpeciesAll>(iParticleSpeciesAll).data(), getName<ParticleSpeciesAll>(iParticleSpeciesAll).data(), getName<ChargeSpecies>(iChargeSpecies).data()), Form(";;;;%s #it{n}#it{#sigma}_{%s};", getName<Detector, NameKind::Display>(iDetector).data(), ParticleSpeciesAllTitles[iParticleSpeciesAll].data()), hcsQaPid);
            }
          }
        } else {
          for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
            hrQaPid.add(Form("QaPid/hCentralityPtEta%sNSigma%s_%s", getName(Detector::Tpc).data(), getName<ParticleSpeciesAll>(iParticleSpeciesAll).data(), getName<ChargeSpecies, NameKind::Lower>(iChargeSpecies).data()), Form(";;;;%s #it{n}#it{#sigma}_{%s};", getName<NameKind::Display>(Detector::Tpc).data(), ParticleSpeciesAllTitles[iParticleSpeciesAll].data()), hcsQaPid);
          }
          for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
            hrQaPid.add(Form("QaPid/hCentralityPtEta%sNSigma%s_%s%s%s", getName(Detector::Tpc).data(), getName<ParticleSpeciesAll>(iParticleSpeciesAll).data(), getName<NameKind::Lower>(Detector::Tof).data(), getName<ParticleSpeciesAll>(iParticleSpeciesAll).data(), getName<ChargeSpecies>(iChargeSpecies).data()), Form(";;;;%s #it{n}#it{#sigma}_{%s};", getName<NameKind::Display>(Detector::Tpc).data(), ParticleSpeciesAllTitles[iParticleSpeciesAll].data()), hcsQaPid);
          }
          for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
            hrQaPid.add(Form("QaPid/hCentralityPtEta%sNSigma%s_%s%s%s", getName(Detector::Tof).data(), getName<ParticleSpeciesAll>(iParticleSpeciesAll).data(), getName<NameKind::Lower>(Detector::Tpc).data(), getName<ParticleSpeciesAll>(iParticleSpeciesAll).data(), getName<ChargeSpecies>(iChargeSpecies).data()), Form(";;;;%s #it{n}#it{#sigma}_{%s};", getName<NameKind::Display>(Detector::Tof).data(), ParticleSpeciesAllTitles[iParticleSpeciesAll].data()), hcsQaPid);
          }
          for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
            hrQaPid.add(Form("QaPid/hCentralityPtEta%sNSigma%s_%s", getName(PidStrategy::TpcTof).data(), getName<ParticleSpeciesAll>(iParticleSpeciesAll).data(), getName<ChargeSpecies, NameKind::Lower>(iChargeSpecies).data()), Form(";;;;%s #it{n}#it{#sigma}_{%s};", getName<NameKind::Display>(PidStrategy::TpcTof).data(), ParticleSpeciesAllTitles[iParticleSpeciesAll].data()), hcsQaPid);
          }
        }
      }
    }

    if (doprocessMc.value) {
      if (groupAnalysis.cfgFlagQaMc.value) {
        LOG(info) << "Enabling MC QA.";

        const AxisSpec asCentrality(20, 0., 100., "Centrality (%)");

        hrQaMc.add("QaMc/hCentralityVzMcDeltaVz", "", {HistType::kTHnSparseF, {asCentrality, {static_cast<std::int32_t>(std::llrint(std::ceil(groupEvent.cfgCutMaxAbsVz.value))) * 20, -std::ceil(groupEvent.cfgCutMaxAbsVz.value), std::ceil(groupEvent.cfgCutMaxAbsVz.value), "#it{V}_{#it{z}}^{Gen} (cm)"}, {200, -0.2, 0.2, "#it{V}_{#it{z}}^{Rec}#minus#it{V}_{#it{z}}^{Gen} (cm)"}}});
        hrQaMc.add("QaMc/hCentralityPtMcEtaMcDeltaPt", "", {HistType::kTHnSparseF, {asCentrality, {200, 0., 2., "#it{p}_{T}^{Gen} (GeV/#it{c})"}, {24, -1.2, 1.2, "#it{#eta}_{Gen}"}, {320, -0.8, 0.8, "#it{p}_{T}^{Rec}#minus#it{p}_{T}^{Gen} (GeV/#it{c})"}}});
        hrQaMc.add("QaMc/hCentralityPtMcEtaMcDeltaEta", "", {HistType::kTHnSparseF, {asCentrality, {20, 0., 2., "#it{p}_{T}^{Gen} (GeV/#it{c})"}, {240, -1.2, 1.2, "#it{#eta}_{Gen}"}, {160, -0.4, 0.4, "#it{#eta}_{Rec}#minus#it{#eta}_{Gen}"}}});
      }
    }

    for (std::int32_t const& iParticleSpecies : std::views::iota(0, NEs<ParticleSpecies>)) {
      if (!static_cast<bool>(groupAnalysis.cfgFlagsCalculationYield.value.get(iParticleSpecies))) {
        continue;
      }
      LOG(info) << "Enabling " << getName<ParticleSpecies, NameKind::DisplayLower>(iParticleSpecies) << " yield calculation.";

      const HistogramConfigSpec hcsCalculationYield(HistType::kTHnSparseF, {{static_cast<std::int32_t>(std::llrint(std::ceil(groupEvent.cfgCutMaxAbsVz.value))) * 2, -std::ceil(groupEvent.cfgCutMaxAbsVz.value), std::ceil(groupEvent.cfgCutMaxAbsVz.value), "#it{V}_{#it{z}} (cm)"}, {{0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.}, "Centrality (%)"}, {40, 0., 2., "#it{p}_{T} (GeV/#it{c})"}, {32, -0.8, 0.8, "#it{#eta}"}});

      if (doprocessMc.value) {
        for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
          hrCalculationYield.add(Form("CalculationYield/hVzCentralityPtMcEtaMc_mc%s%s", getName<ParticleSpecies>(iParticleSpecies).data(), getName<ChargeSpecies>(iChargeSpecies).data()), "", hcsCalculationYield);
        }
        for (std::int32_t const& iPidStrategy : std::views::iota(0, NEs<PidStrategy>)) {
          for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
            if (groupTrack.cfgFlagMcParticleMomentum.value) {
              hrCalculationYield.add(Form("CalculationYield/hVzCentralityPtMcEtaMc_mc%s%s%s", getName<PidStrategy>(iPidStrategy).data(), getName<ParticleSpecies>(iParticleSpecies).data(), getName<ChargeSpecies>(iChargeSpecies).data()), "", hcsCalculationYield);
            } else {
              hrCalculationYield.add(Form("CalculationYield/hVzCentralityPtEta_mc%s%s%s", getName<PidStrategy>(iPidStrategy).data(), getName<ParticleSpecies>(iParticleSpecies).data(), getName<ChargeSpecies>(iChargeSpecies).data()), "", hcsCalculationYield);
            }
          }
        }
      } else {
        for (std::int32_t const& iPidStrategy : std::views::iota(0, NEs<PidStrategy>)) {
          for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
            hrCalculationYield.add(Form("CalculationYield/hVzCentralityPtEta_%s%s%s", getName<PidStrategy, NameKind::Lower>(iPidStrategy).data(), getName<ParticleSpecies>(iParticleSpecies).data(), getName<ChargeSpecies>(iChargeSpecies).data()), "", hcsCalculationYield);
          }
        }
      }
    }

    if (doprocessMc.value) {
      for (std::int32_t const& iParticleSpecies : std::views::iota(0, NEs<ParticleSpecies>)) {
        if (!static_cast<bool>(groupAnalysis.cfgFlagsCalculationPurity.value.get(iParticleSpecies))) {
          continue;
        }

        LOG(info) << "Enabling " << getName<ParticleSpecies, NameKind::DisplayLower>(iParticleSpecies) << " purity calculation.";

        const HistogramConfigSpec hcsCalculationPurity(HistType::kTProfile3D, {{{0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.}, "Centrality (%)"}, {20, 0., 2., "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#it{#eta}"}});

        for (std::int32_t const& iPidStrategy : std::views::iota(0, NEs<PidStrategy>)) {
          for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
            hrCalculationPurity.add(Form("CalculationPurity/pCentralityPtEtaPurity%s%s%s", getName<PidStrategy>(iPidStrategy).data(), getName<ParticleSpecies>(iParticleSpecies).data(), getName<ChargeSpecies>(iChargeSpecies).data()), "", hcsCalculationPurity);
          }
        }
      }
    }

    if (doprocessMc.value) {
      for (std::int32_t const& iParticleSpecies : std::views::iota(0, NEs<ParticleSpecies>)) {
        if (!static_cast<bool>(groupAnalysis.cfgFlagsCalculationFractionPrimary.value.get(iParticleSpecies))) {
          continue;
        }

        LOG(info) << "Enabling " << getName<ParticleSpecies, NameKind::DisplayLower>(iParticleSpecies) << " primary fraction calculation.";

        const HistogramConfigSpec hcsCalculationFractionPrimary(HistType::kTProfile3D, {{{0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.}, "Centrality (%)"}, {20, 0., 2., "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#it{#eta}"}});

        for (std::int32_t const& iPidStrategy : std::views::iota(0, NEs<PidStrategy>)) {
          for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
            hrCalculationFractionPrimary.add(Form("CalculationFractionPrimary/pCentralityPtEtaFractionPrimary%s%s%s", getName<PidStrategy>(iPidStrategy).data(), getName<ParticleSpecies>(iParticleSpecies).data(), getName<ChargeSpecies>(iChargeSpecies).data()), "", hcsCalculationFractionPrimary);
          }
        }
      }
    }

    if (nEnabled(groupAnalysis.cfgFlagsCalculationFluctuation) > 1) {
      LOG(fatal) << "Invalid " << groupAnalysis.cfgFlagsCalculationFluctuation.name << "!";
    }
    if (isEnabled(groupAnalysis.cfgFlagsCalculationFluctuation) && groupEvent.cfgNSubgroups.value <= 0) {
      LOG(fatal) << "Invalid " << groupEvent.cfgNSubgroups.name << "!";
    }
    for (std::int32_t const& iParticleNumber : std::views::iota(0, NEs<ParticleNumber>)) {
      if (!static_cast<bool>(groupAnalysis.cfgFlagsCalculationFluctuation.value.get(iParticleNumber))) {
        continue;
      }

      LOG(info) << "Enabling " << getName<ParticleNumber, NameKind::DisplayLower>(iParticleNumber) << " number fluctuation calculation.";

      const AxisSpec asCentrality(groupEvent.cfgAxisCentrality, "Centrality (%)");
      const HistogramConfigSpec hcsDistribution(HistType::kTHnSparseD, {asCentrality, {200, -0.5, 199.5}, {200, -0.5, 199.5}});
      const HistogramConfigSpec hcsFluctuationCalculator(HistType::kTH3D, {asCentrality, {groupEvent.cfgNSubgroups.value, -0.5, groupEvent.cfgNSubgroups.value - 0.5, "Subgroup Index"}, {fluctuation_calculator_base::NOrderKeys, -0.5, fluctuation_calculator_base::NOrderKeys - 0.5, "Order Key Index"}});

      for (std::int32_t const& iChargeNumber : std::views::iota(0, NEs<ChargeNumber>)) {
        fluctuationCalculatorTrack[iParticleNumber][iChargeNumber] = std::make_unique<FluctuationCalculatorTrack>();
      }

      constexpr std::array<std::string_view, NEs<ParticleNumber>> ParticleNumberTitles{"h", "K", "p"};
      constexpr std::array<std::string_view, NEs<ChargeSpecies>> ChargeSpeciesTitles{"+", "#minus"};
      if (doprocessMc.value) {
        hrCalculationFluctuation.add(Form("CalculationFluctuation/hCentralityN%s%sN%s%s_mc", getName<ParticleNumber>(iParticleNumber).data(), getName(ChargeSpecies::Plus).data(), getName<ParticleNumber>(iParticleNumber).data(), getName(ChargeSpecies::Minus).data()), Form(";;#it{N}(%s^{%s});#it{N}(%s^{%s});", ParticleNumberTitles[iParticleNumber].data(), ChargeSpeciesTitles[toI(ChargeSpecies::Plus)].data(), ParticleNumberTitles[iParticleNumber].data(), ChargeSpeciesTitles[toI(ChargeSpecies::Minus)].data()), hcsDistribution);
        hrCalculationFluctuation.add(Form("CalculationFluctuation/hCentralityN%s%sN%s%s_mcEff", getName<ParticleNumber>(iParticleNumber).data(), getName(ChargeSpecies::Plus).data(), getName<ParticleNumber>(iParticleNumber).data(), getName(ChargeSpecies::Minus).data()), Form(";;#it{N}(%s^{%s});#it{N}(%s^{%s});", ParticleNumberTitles[iParticleNumber].data(), ChargeSpeciesTitles[toI(ChargeSpecies::Plus)].data(), ParticleNumberTitles[iParticleNumber].data(), ChargeSpeciesTitles[toI(ChargeSpecies::Minus)].data()), hcsDistribution);
        for (std::int32_t const& iChargeNumber : std::views::iota(0, NEs<ChargeNumber>)) {
          hrCalculationFluctuation.add(Form("CalculationFluctuation/hFluctuationCalculator%s%s_mc", getName<ParticleNumber>(iParticleNumber).data(), getName<ChargeNumber>(iChargeNumber).data()), "", hcsFluctuationCalculator);
        }
      }
      hrCalculationFluctuation.add(Form("CalculationFluctuation/hCentralityN%s%sN%s%s", getName<ParticleNumber>(iParticleNumber).data(), getName(ChargeSpecies::Plus).data(), getName<ParticleNumber>(iParticleNumber).data(), getName(ChargeSpecies::Minus).data()), Form(";;#it{N}(%s^{%s});#it{N}(%s^{%s});", ParticleNumberTitles[iParticleNumber].data(), ChargeSpeciesTitles[toI(ChargeSpecies::Plus)].data(), ParticleNumberTitles[iParticleNumber].data(), ChargeSpeciesTitles[toI(ChargeSpecies::Minus)].data()), hcsDistribution);
      for (std::int32_t const& iChargeNumber : std::views::iota(0, NEs<ChargeNumber>)) {
        hrCalculationFluctuation.add(Form("CalculationFluctuation/hFluctuationCalculator%s%s", getName<ParticleNumber>(iParticleNumber).data(), getName<ChargeNumber>(iChargeNumber).data()), "", hcsFluctuationCalculator);
      }
    }

    for (std::int32_t const& iParticleSpecies : std::views::iota(0, NEs<ParticleSpecies>)) {
      if (!static_cast<bool>(groupAnalysis.cfgFlagsCalculationFluctuation.value.get(toI(ParticleNumber::Charge))) && (iParticleSpecies != toI(ParticleSpecies::Kaon) || !static_cast<bool>(groupAnalysis.cfgFlagsCalculationFluctuation.value.get(toI(ParticleNumber::Kaon)))) && (iParticleSpecies != toI(ParticleSpecies::Proton) || !static_cast<bool>(groupAnalysis.cfgFlagsCalculationFluctuation.value.get(toI(ParticleNumber::Proton))))) {
        continue;
      }

      holderCcdb.hVzCentralityPtEtaEfficiency.resize(nRunGroups);
      for (std::int32_t const& iRunGroup : std::views::iota(0, nRunGroups)) {
        const TList* const lRunGroup{readListRunGroup(iRunGroup + 1)};
        for (std::int32_t const& iPidStrategy : std::views::iota(0, NEs<PidStrategy>)) {
          for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
            const char* const name{Form("hVzCentralityPtEtaEfficiency%s%s%s_runGroup%d", getName<PidStrategy>(iPidStrategy).data(), getName<ParticleSpecies>(iParticleSpecies).data(), getName<ChargeSpecies>(iChargeSpecies).data(), iRunGroup + 1)};
            holderCcdb.hVzCentralityPtEtaEfficiency[iRunGroup][iPidStrategy][iParticleSpecies][iChargeSpecies] = dynamic_cast<const THnBase*>(lRunGroup->FindObject(name));
            if (!holderCcdb.hVzCentralityPtEtaEfficiency[iRunGroup][iPidStrategy][iParticleSpecies][iChargeSpecies] || holderCcdb.hVzCentralityPtEtaEfficiency[iRunGroup][iPidStrategy][iParticleSpecies][iChargeSpecies]->GetNdimensions() != HolderCcdb::NDimensionsEfficiency) {
              LOG(fatal) << "Invalid " << name << "!";
            }
            LOG(info) << "Reading from CCDB: " << name;
          }
        }
      }
    }
  }

  template <PidStrategy PidStrategyValue, ParticleSpecies ParticleSpeciesValue, ChargeSpecies ChargeSpeciesValue>
    requires IsValid<PidStrategyValue, ParticleSpeciesValue, ChargeSpeciesValue>
  double getEfficiency(const bool doUsingMcParticleMomentum)
  {
    const THnBase* const hVzCentralityPtEtaEfficiency{holderCcdb.hVzCentralityPtEtaEfficiency.at(std::abs(holderEvent.runGroupIndex) - 1)[toI(PidStrategyValue)][toI(ParticleSpeciesValue)][toI(ChargeSpeciesValue)]};
    return hVzCentralityPtEtaEfficiency ? hVzCentralityPtEtaEfficiency->GetBinContent(hVzCentralityPtEtaEfficiency->GetBin(std::array<double, HolderCcdb::NDimensionsEfficiency>{holderEvent.vz, holderEvent.centrality, doUsingMcParticleMomentum ? holderMcParticle.pt : holderTrack.pt, doUsingMcParticleMomentum ? holderMcParticle.eta : holderTrack.eta}.data())) : 0.;
  }

  template <Detector DetectorValue, ParticleSpecies ParticleSpeciesValue>
    requires IsValid<ParticleSpeciesValue, DetectorValue>
  double getShiftNSigmaPid()
  {
    if (!groupTrack.cfgFlagsRecalibrationNSigmaPid.value.get(toI(ParticleSpeciesValue))) {
      return 0.;
    }

    static const auto clampInAxis{[](const double value, const TAxis* const axis) {
      const std::int32_t first{std::clamp(axis->GetFirst(), 1, axis->GetNbins())};
      const std::int32_t last{std::clamp(axis->GetLast(), 1, axis->GetNbins())};
      return first == last ? axis->GetBinCenter(first) : std::clamp(value, std::nextafter(axis->GetBinCenter(first), std::numeric_limits<double>::infinity()), std::nextafter(axis->GetBinCenter(last), -std::numeric_limits<double>::infinity()));
    }};

    if (holderTrack.sign == 0) {
      return 0.;
    }
    const TH3* const hCentralityPtEtaShiftNSigmaPid{holderCcdb.hCentralityPtEtaShiftNSigmaPid.at(std::abs(holderEvent.runGroupIndex) - 1)[toI(DetectorValue)][toI(ParticleSpeciesValue)][holderTrack.sign > 0 ? toI(ChargeSpecies::Plus) : toI(ChargeSpecies::Minus)]};
    return hCentralityPtEtaShiftNSigmaPid ? hCentralityPtEtaShiftNSigmaPid->Interpolate(clampInAxis(holderEvent.centrality, hCentralityPtEtaShiftNSigmaPid->GetXaxis()), clampInAxis(holderTrack.pt, hCentralityPtEtaShiftNSigmaPid->GetYaxis()), clampInAxis(holderTrack.eta, hCentralityPtEtaShiftNSigmaPid->GetZaxis())) : 0.;
  }

  template <PidStrategyAll PidStrategyAllValue, ParticleSpeciesAll ParticleSpeciesAllValue>
    requires IsValid<ParticleSpeciesAllValue, PidStrategyAllValue>
  bool isPid(const bool doRejectingOthers)
  {
    if constexpr (ParticleSpeciesAllValue == ParticleSpeciesAll::All) {
      if constexpr (PidStrategyAllValue == PidStrategyAll::Tpc) {
        if (!holderTrack.hasPid[toI(Detector::Tpc)]) {
          return false;
        }
      } else if constexpr (PidStrategyAllValue == PidStrategyAll::Tof) {
        if (!holderTrack.hasPid[toI(Detector::Tof)]) {
          return false;
        }
      } else {
        if (!holderTrack.hasPid[toI(Detector::Tpc)] || !holderTrack.hasPid[toI(Detector::Tof)]) {
          return false;
        }
      }
    } else {
      constexpr std::int32_t ParticleSpeciesIndex{toI(getValue<ParticleSpecies>(ParticleSpeciesAllValue))};
      if constexpr (PidStrategyAllValue == PidStrategyAll::TpcTofSeparated) {
        if (!(std::abs(holderTrack.nSigmaPid[toI(PidStrategyAll::Tpc)][ParticleSpeciesIndex]) < groupTrack.cfgCutsMaxAbsNSigmaPid.value.get(ParticleSpeciesIndex))) {
          return false;
        }
        if (!(std::abs(holderTrack.nSigmaPid[toI(PidStrategyAll::Tof)][ParticleSpeciesIndex]) < groupTrack.cfgCutsMaxAbsNSigmaPid.value.get(ParticleSpeciesIndex))) {
          return false;
        }
        if (doRejectingOthers && !(std::abs(holderTrack.nSigmaPid[toI(PidStrategyAll::Tof)][ParticleSpeciesIndex]) < std::min(std::abs(holderTrack.nSigmaPid[toI(PidStrategyAll::Tof)][(ParticleSpeciesIndex + 1) % NEs<ParticleSpecies>]), std::abs(holderTrack.nSigmaPid[toI(PidStrategyAll::Tof)][(ParticleSpeciesIndex + 2) % NEs<ParticleSpecies>])))) {
          return false;
        }
      } else {
        if (!(std::abs(holderTrack.nSigmaPid[toI(PidStrategyAllValue)][ParticleSpeciesIndex]) < groupTrack.cfgCutsMaxAbsNSigmaPid.value.get(ParticleSpeciesIndex))) {
          return false;
        }
        if (doRejectingOthers && !(std::abs(holderTrack.nSigmaPid[toI(PidStrategyAllValue)][ParticleSpeciesIndex]) < std::min(std::abs(holderTrack.nSigmaPid[toI(PidStrategyAllValue)][(ParticleSpeciesIndex + 1) % NEs<ParticleSpecies>]), std::abs(holderTrack.nSigmaPid[toI(PidStrategyAllValue)][(ParticleSpeciesIndex + 2) % NEs<ParticleSpecies>])))) {
          return false;
        }
      }
    }
    return true;
  }

  template <ParticleSpeciesAll ParticleSpeciesAllValue, ChargeSpecies ChargeSpeciesValue>
    requires IsValid<ParticleSpeciesAllValue, ChargeSpeciesValue>
  bool isPid()
  {
    if constexpr (ParticleSpeciesAllValue == ParticleSpeciesAll::All) {
      if constexpr (ChargeSpeciesValue == ChargeSpecies::Plus) {
        if (holderMcParticle.charge <= 0) {
          return false;
        }
      } else { // ChargeSpeciesValue == ChargeSpecies::Minus
        if (holderMcParticle.charge >= 0) {
          return false;
        }
      }
    } else {
      if (holderMcParticle.pdgCode != getPdgCode(ParticleSpeciesAllValue, ChargeSpeciesValue)) {
        return false;
      }
    }
    return true;
  }

  bool isGoodMomentum(const bool doUsingMcParticleMomentum)
  {
    if (doUsingMcParticleMomentum) {
      if (!(groupTrack.cfgCutMinPt.value < holderMcParticle.pt) || !(holderMcParticle.pt < groupTrack.cfgCutMaxPt.value)) {
        return false;
      }
      if (!(std::abs(holderMcParticle.eta) < groupTrack.cfgCutMaxAbsEta.value)) {
        return false;
      }
    } else {
      if (!(groupTrack.cfgCutMinPt.value < holderTrack.pt) || !(holderTrack.pt < groupTrack.cfgCutMaxPt.value)) {
        return false;
      }
      if (!(std::abs(holderTrack.eta) < groupTrack.cfgCutMaxAbsEta.value)) {
        return false;
      }
    }
    return true;
  }

  bool isGoodDca()
  {
    if (!groupTrack.cfgFlagRecalibrationDca.value) {
      for (std::int32_t const& iDcaAxis : std::views::iota(0, NEs<DcaAxis>)) {
        if (!(std::abs(holderTrack.dca[iDcaAxis]) < groupTrack.cfgCutsMaxAbsNSigmaDca.value.get(iDcaAxis))) {
          return false;
        }
      }
    } else {
      if (holderTrack.sign == 0) {
        return false;
      }
      const std::int32_t chargeSpeciesIndex{holderTrack.sign > 0 ? toI(ChargeSpecies::Plus) : toI(ChargeSpecies::Minus)};
      const std::array<std::array<std::array<const TFormula*, NEs<ChargeSpecies>>, NEs<DcaAxis>>, NEs<DcaKind>>& fPtDcaGroup{holderCcdb.fPtDca.at(std::abs(holderEvent.runGroupIndex) - 1)};
      for (std::int32_t const& iDcaAxis : std::views::iota(0, NEs<DcaAxis>)) {
        if (!fPtDcaGroup[toI(DcaKind::Mean)][iDcaAxis][chargeSpeciesIndex] || !fPtDcaGroup[toI(DcaKind::Sigma)][iDcaAxis][chargeSpeciesIndex] || !(std::abs(holderTrack.dca[iDcaAxis] - fPtDcaGroup[toI(DcaKind::Mean)][iDcaAxis][chargeSpeciesIndex]->Eval(holderTrack.pt)) < groupTrack.cfgCutsMaxAbsNSigmaDca.value.get(iDcaAxis) * fPtDcaGroup[toI(DcaKind::Sigma)][iDcaAxis][chargeSpeciesIndex]->Eval(holderTrack.pt))) {
          return false;
        }
      }
    }
    return true;
  }

  template <typename T>
  bool isGoodTrack(const T& track)
  {
    if (groupTrack.cfgFlagPvContributor.value && !track.isPVContributor()) {
      return false;
    }
    if (!(track.itsNCls() > groupTrack.cfgCutMinItsNCls.value)) {
      return false;
    }
    if (!(track.itsChi2NCl() < groupTrack.cfgCutMaxItsChi2NCls.value)) {
      return false;
    }
    if (!(track.tpcNClsFound() > groupTrack.cfgCutMinTpcNCls.value)) {
      return false;
    }
    if (!(track.tpcChi2NCl() < groupTrack.cfgCutMaxTpcChi2NCls.value)) {
      return false;
    }
    if (!(track.tpcFractionSharedCls() < groupTrack.cfgCutMaxTpcNClsSharedRatio.value)) {
      return false;
    }
    if (!(track.tpcNClsCrossedRows() > groupTrack.cfgCutMinTpcNCrossedRows.value)) {
      return false;
    }
    if (!(track.tpcCrossedRowsOverFindableCls() > groupTrack.cfgCutMinTpcNCrossedRowsRatio.value)) {
      return false;
    }
    return true;
  }

  template <bool IsMc, typename MP>
  bool isGoodMcParticle(const MP& mcParticle)
  {
    if constexpr (IsMc) {
      if (!mcParticle.isPhysicalPrimary()) {
        return false;
      }
    }
    return true;
  }

  template <ChargeSpecies ChargeSpeciesValue, typename T>
    requires IsValid<ChargeSpeciesValue>
  void fillQaRunByTrackByChargeSpecies(const T& track)
  {
    const auto fill{[&](const auto& name, const auto value) {
      hrQaRun.fill(C_CS("QaRun/pRunIndex") + name + C_CS("_") + C_SV(getName<NameKind::Lower>(ChargeSpeciesValue)), holderEvent.runIndex, value);
    }};
    const auto fillNSigmaPidByDetectorParticleSpecies{
      [&]<Detector DetectorValue, ParticleSpecies ParticleSpeciesValue>
        requires IsValid<DetectorValue, getValue<PidStrategyAll>(DetectorValue), ParticleSpeciesValue>
      () {
        const double nSigmaPid{holderTrack.nSigmaPid[toI(getValue<PidStrategyAll>(DetectorValue))][toI(ParticleSpeciesValue)]};
        if (std::abs(nSigmaPid) < HolderTrack::TruncationAbsNSigmaPid) {
          fill(C_SV(getName(DetectorValue)) + C_CS("NSigma") + C_SV(getName(ParticleSpeciesValue)), nSigmaPid);
        }
      }};

    fill(C_CS("ItsNCls"), track.itsNCls());
    fill(C_CS("ItsChi2NCls"), track.itsChi2NCl());
    fill(C_CS("TpcNCls"), track.tpcNClsFound());
    fill(C_CS("TpcChi2NCls"), track.tpcChi2NCl());
    fill(C_CS("TpcNClsSharedRatio"), track.tpcFractionSharedCls());
    fill(C_CS("TpcNCrossedRows"), track.tpcNClsCrossedRows());
    fill(C_CS("TpcNCrossedRowsRatio"), track.tpcCrossedRowsOverFindableCls());
    fill(C_CS("Dca") + C_SV(getName(DcaAxis::Xy)), holderTrack.dca[toI(DcaAxis::Xy)]);
    fill(C_CS("Dca") + C_SV(getName(DcaAxis::Z)), holderTrack.dca[toI(DcaAxis::Z)]);
    fill(C_CS("Pt"), holderTrack.pt);
    fill(C_CS("Eta"), holderTrack.eta);
    fill(C_CS("Phi"), holderTrack.phi);
    if (holderTrack.hasPid[toI(Detector::Tpc)]) {
      fill(C_SV(getName(Detector::Tpc)) + C_CS("DeDx"), track.tpcSignal());
      fillNSigmaPidByDetectorParticleSpecies.template operator()<Detector::Tpc, ParticleSpecies::Pion>();
      fillNSigmaPidByDetectorParticleSpecies.template operator()<Detector::Tpc, ParticleSpecies::Kaon>();
      fillNSigmaPidByDetectorParticleSpecies.template operator()<Detector::Tpc, ParticleSpecies::Proton>();
    }
    if (holderTrack.hasPid[toI(Detector::Tof)]) {
      fill(C_SV(getName(Detector::Tof)) + C_CS("InverseBeta"), 1. / track.beta());
      fillNSigmaPidByDetectorParticleSpecies.template operator()<Detector::Tof, ParticleSpecies::Pion>();
      fillNSigmaPidByDetectorParticleSpecies.template operator()<Detector::Tof, ParticleSpecies::Kaon>();
      fillNSigmaPidByDetectorParticleSpecies.template operator()<Detector::Tof, ParticleSpecies::Proton>();
    }
  }

  template <ChargeSpecies ChargeSpeciesValue>
    requires IsValid<ChargeSpeciesValue>
  void fillQaRunByEventByChargeSpecies()
  {
    const auto fill{[&](const auto& name, const auto value) {
      hrQaRun.fill(C_CS("QaRun/pRunIndex") + name + C_CS("_") + C_SV(getName<NameKind::Lower>(ChargeSpeciesValue)), holderEvent.runIndex, value);
    }};

    fill(C_CS("NGlobalTracks"), holderEvent.nGlobalTracks[toI(ChargeSpeciesValue)]);
    fill(C_CS("NPvContributors"), holderEvent.nPvContributors[toI(ChargeSpeciesValue)]);
    if (holderEvent.nGlobalTracks[toI(ChargeSpeciesValue)] > 0) {
      fill(C_SV(getName(DcaKind::Mean)) + C_CS("Dca") + C_SV(getName(DcaAxis::Xy)), holderEvent.dca[toI(DcaKind::Mean)][toI(DcaAxis::Xy)][toI(ChargeSpeciesValue)]);
      fill(C_SV(getName(DcaKind::Sigma)) + C_CS("Dca") + C_SV(getName(DcaAxis::Xy)), holderEvent.dca[toI(DcaKind::Sigma)][toI(DcaAxis::Xy)][toI(ChargeSpeciesValue)]);
      fill(C_SV(getName(DcaKind::Mean)) + C_CS("Dca") + C_SV(getName(DcaAxis::Z)), holderEvent.dca[toI(DcaKind::Mean)][toI(DcaAxis::Z)][toI(ChargeSpeciesValue)]);
      fill(C_SV(getName(DcaKind::Sigma)) + C_CS("Dca") + C_SV(getName(DcaAxis::Z)), holderEvent.dca[toI(DcaKind::Sigma)][toI(DcaAxis::Z)][toI(ChargeSpeciesValue)]);
    }
    fill(C_CS("NTofBeta"), holderEvent.nTofBeta[toI(ChargeSpeciesValue)]);
  }

  template <ChargeSpecies ChargeSpeciesValue, typename T>
    requires IsValid<ChargeSpeciesValue>
  void fillQaTrackByChargeSpecies(const T& track)
  {
    const auto fill{[&](const auto& name, const auto... positionAndWeight) {
      hrQaTrack.fill(C_CS("QaTrack/h") + name + C_CS("_") + C_SV(getName<NameKind::Lower>(ChargeSpeciesValue)), positionAndWeight...);
    }};

    fill(C_CS("ItsNCls"), track.itsNCls());
    fill(C_CS("ItsChi2NCls"), track.itsChi2NCl());
    fill(C_CS("TpcNClsNClsShared"), track.tpcNClsFound(), track.tpcNClsShared());
    fill(C_CS("TpcChi2NCls"), track.tpcChi2NCl());
    fill(C_CS("TpcNClsFindableNCrossedRows"), track.tpcNClsFindable(), track.tpcNClsCrossedRows());
  }

  template <ChargeSpecies ChargeSpeciesValue>
    requires IsValid<ChargeSpeciesValue>
  void fillQaDcaByChargeSpecies()
  {
    const auto fillByDcaAxis{
      [&]<DcaAxis DcaAxisValue>
        requires IsValid<DcaAxisValue>
      () {
        hrQaDca.fill(C_CS("QaDca/hPtDca") + C_SV(getName(DcaAxisValue)) + C_CS("_") + C_SV(getName<NameKind::Lower>(ChargeSpeciesValue)), holderTrack.pt, holderTrack.dca[toI(DcaAxisValue)]);
        hrQaDca.fill(C_CS("QaDca/pPtEtaPhiIuDca") + C_SV(getName(DcaAxisValue)) + C_CS("_") + C_SV(getName<NameKind::Lower>(ChargeSpeciesValue)), holderTrack.pt, holderTrack.eta, holderTrack.phiIu, holderTrack.dca[toI(DcaAxisValue)]);
      }};

    fillByDcaAxis.template operator()<DcaAxis::Xy>();
    fillByDcaAxis.template operator()<DcaAxis::Z>();
  }

  template <ParticleSpeciesAll ParticleSpeciesAllValue, typename T>
    requires IsValid<ParticleSpeciesAllValue>
  void fillQaAcceptancebyParticleSpeciesAll(const T& track)
  {
    if (!groupAnalysis.cfgFlagsQaAcceptance.value.get(toI(ParticleSpeciesAllValue)) || holderTrack.sign == 0) {
      return;
    }

    const auto fillByChargeSpecies{
      [&]<ChargeSpecies ChargeSpeciesValue>
        requires IsValid<ChargeSpeciesValue>
      (const auto& name, const auto value) {
        const auto fillByPidStrategy{
          [&]<PidStrategy PidStrategyValue>
            requires IsValid<PidStrategyValue>
          () {
            if (isPid<getValue<PidStrategyAll>(PidStrategyValue), ParticleSpeciesAllValue>(false)) {
              hrQaAcceptance.fill(C_CS("QaAcceptance/h") + name + C_CS("Pt_") + C_SV(getName<NameKind::Lower>(PidStrategyValue)) + C_CS("Edge") + C_SV(getName(ParticleSpeciesAllValue)) + C_SV(getName(ChargeSpeciesValue)), value, holderTrack.pt); // NOLINT(clang-analyzer-core.NonNullParamChecker)
            }
          }};

        fillByPidStrategy.template operator()<PidStrategy::Tpc>();
        fillByPidStrategy.template operator()<PidStrategy::TpcTof>();
      }};

    if constexpr (ParticleSpeciesAllValue == ParticleSpeciesAll::All) {
      if (holderTrack.sign > 0) {
        fillByChargeSpecies.template operator()<ChargeSpecies::Plus>(C_CS("Eta"), holderTrack.eta);
      } else {
        fillByChargeSpecies.template operator()<ChargeSpecies::Minus>(C_CS("Eta"), holderTrack.eta);
      }
    } else {
      if (holderTrack.sign > 0) {
        fillByChargeSpecies.template operator()<ChargeSpecies::Plus>(C_CS("Rapidity"), track.rapidity(getMass(ParticleSpeciesAllValue)));
      } else {
        fillByChargeSpecies.template operator()<ChargeSpecies::Minus>(C_CS("Rapidity"), track.rapidity(getMass(ParticleSpeciesAllValue)));
      }
    }
  }

  template <DataMode DataModeValue, ParticleSpeciesAll ParticleSpeciesAllValue>
    requires IsValid<DataModeValue, ParticleSpeciesAllValue> && (DataModeValue != DataMode::McMcParticle)
  void fillQaPhiByParticleSpeciesAll()
  {
    if (!groupAnalysis.cfgFlagsQaPhi.value.get(toI(ParticleSpeciesAllValue))) {
      return;
    }

    if (holderTrack.sign == 0) { // DataModeValue != DataMode::McMcParticle
      return;
    }

    const auto fillByChargeSpecies{
      [&]<ChargeSpecies ChargeSpeciesValue>
        requires IsValid<ChargeSpeciesValue>
      () {
        const auto fillByPidStrategy{
          [&]<PidStrategy PidStrategyValue>
            requires IsValid<PidStrategyValue>
          () {
            if constexpr (DataModeValue == DataMode::McTrack) {
              if (isPid<ParticleSpeciesAllValue, ChargeSpeciesValue>() && isPid<getValue<PidStrategyAll>(PidStrategyValue), ParticleSpeciesAllValue>(false)) {
                hrQaPhi.fill(C_CS("QaPhi/hCentralityPtEtaPhi_mc") + C_SV(getName(PidStrategyValue)) + C_SV(getName(ParticleSpeciesAllValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.centrality, holderTrack.pt, holderTrack.eta, holderTrack.phi);
                hrQaPhi.fill(C_CS("QaPhi/hCentralityPtEtaPhiIu_mc") + C_SV(getName(PidStrategyValue)) + C_SV(getName(ParticleSpeciesAllValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.centrality, holderTrack.pt, holderTrack.eta, holderTrack.phiIu);
              }
            } else { // DataModeValue == DataMode::RawTrack
              if (isPid<getValue<PidStrategyAll>(PidStrategyValue), ParticleSpeciesAllValue>(false)) {
                hrQaPhi.fill(C_CS("QaPhi/hCentralityPtEtaPhi_") + C_SV(getName<NameKind::Lower>(PidStrategyValue)) + C_SV(getName(ParticleSpeciesAllValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.centrality, holderTrack.pt, holderTrack.eta, holderTrack.phi);
                hrQaPhi.fill(C_CS("QaPhi/hCentralityPtEtaPhiIu_") + C_SV(getName<NameKind::Lower>(PidStrategyValue)) + C_SV(getName(ParticleSpeciesAllValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.centrality, holderTrack.pt, holderTrack.eta, holderTrack.phiIu);
              }
            }
          }};

        fillByPidStrategy.template operator()<PidStrategy::Tpc>();
        fillByPidStrategy.template operator()<PidStrategy::TpcTof>();
      }};

    if (holderTrack.sign > 0) {
      fillByChargeSpecies.template operator()<ChargeSpecies::Plus>();
    } else {
      fillByChargeSpecies.template operator()<ChargeSpecies::Minus>();
    }
  }

  template <DataMode DataModeValue, ParticleSpeciesAll ParticleSpeciesAllValue, typename T>
    requires IsValid<DataModeValue, ParticleSpeciesAllValue> && (DataModeValue != DataMode::McMcParticle)
  void fillQaPidByParticleSpeciesAll(const T& track)
  {
    if (!groupAnalysis.cfgFlagsQaPid.value.get(toI(ParticleSpeciesAllValue)) || holderTrack.sign == 0) {
      return;
    }

    if constexpr (ParticleSpeciesAllValue == ParticleSpeciesAll::All) {
      if (isPid<getValue<PidStrategyAll>(PidStrategy::Tpc), ParticleSpeciesAll::All>(false)) {
        hrQaPid.fill(C_CS("QaPid/hCentralityPOverQEtaTpcLnDeDx"), holderEvent.centrality, holderTrack.p / holderTrack.sign, holderTrack.eta, track.tpcSignal());
      }
      if (isPid<getValue<PidStrategyAll>(PidStrategy::TpcTof), ParticleSpeciesAll::All>(false)) {
        hrQaPid.fill(C_CS("QaPid/hCentralityPOverQEtaTofInverseBeta"), holderEvent.centrality, holderTrack.p / holderTrack.sign, holderTrack.eta, 1. / track.beta());
      }
    } else {
      const auto fillByChargeSpecies{
        [&]<ChargeSpecies ChargeSpeciesValue>
          requires IsValid<ChargeSpeciesValue>
        () {
          if constexpr (DataModeValue == DataMode::McTrack) {
            if (isPid<ParticleSpeciesAllValue, ChargeSpeciesValue>()) {
              hrQaPid.fill(C_CS("QaPid/hCentralityPtEta") + C_SV(getName(Detector::Tpc)) + C_CS("NSigma") + C_SV(getName(ParticleSpeciesAllValue)) + C_CS("_mc") + C_SV(getName(ParticleSpeciesAllValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.centrality, holderTrack.pt, holderTrack.eta, holderTrack.nSigmaPid[toI(PidStrategyAll::Tpc)][toI(ParticleSpeciesAllValue)]);
              hrQaPid.fill(C_CS("QaPid/hCentralityPtEta") + C_SV(getName(Detector::Tof)) + C_CS("NSigma") + C_SV(getName(ParticleSpeciesAllValue)) + C_CS("_mc") + C_SV(getName(ParticleSpeciesAllValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.centrality, holderTrack.pt, holderTrack.eta, holderTrack.nSigmaPid[toI(PidStrategyAll::Tof)][toI(ParticleSpeciesAllValue)]);
            }
          } else { // DataModeValue == DataMode::RawTrack
            hrQaPid.fill(C_CS("QaPid/hCentralityPtEta") + C_SV(getName(Detector::Tpc)) + C_CS("NSigma") + C_SV(getName(ParticleSpeciesAllValue)) + C_CS("_") + C_SV(getName<NameKind::Lower>(ChargeSpeciesValue)), holderEvent.centrality, holderTrack.pt, holderTrack.eta, holderTrack.nSigmaPid[toI(PidStrategyAll::Tpc)][toI(ParticleSpeciesAllValue)]);
            if (isPid<PidStrategyAll::Tof, ParticleSpeciesAllValue>(false)) {
              hrQaPid.fill(C_CS("QaPid/hCentralityPtEta") + C_SV(getName(Detector::Tpc)) + C_CS("NSigma") + C_SV(getName(ParticleSpeciesAllValue)) + C_CS("_") + C_SV(getName<NameKind::Lower>(Detector::Tof)) + C_SV(getName(ParticleSpeciesAllValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.centrality, holderTrack.pt, holderTrack.eta, holderTrack.nSigmaPid[toI(PidStrategyAll::Tpc)][toI(ParticleSpeciesAllValue)]);
            }
            if (isPid<PidStrategyAll::Tpc, ParticleSpeciesAllValue>(false)) {
              hrQaPid.fill(C_CS("QaPid/hCentralityPtEta") + C_SV(getName(Detector::Tof)) + C_CS("NSigma") + C_SV(getName(ParticleSpeciesAllValue)) + C_CS("_") + C_SV(getName<NameKind::Lower>(Detector::Tpc)) + C_SV(getName(ParticleSpeciesAllValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.centrality, holderTrack.pt, holderTrack.eta, holderTrack.nSigmaPid[toI(PidStrategyAll::Tof)][toI(ParticleSpeciesAllValue)]);
            }
            hrQaPid.fill(C_CS("QaPid/hCentralityPtEta") + C_SV(getName(PidStrategy::TpcTof)) + C_CS("NSigma") + C_SV(getName(ParticleSpeciesAllValue)) + C_CS("_") + C_SV(getName<NameKind::Lower>(ChargeSpeciesValue)), holderEvent.centrality, holderTrack.pt, holderTrack.eta, holderTrack.nSigmaPid[toI(PidStrategyAll::TpcTofCombined)][toI(ParticleSpeciesAllValue)]);
          }
        }};

      if (holderTrack.sign > 0) {
        fillByChargeSpecies.template operator()<ChargeSpecies::Plus>();
      } else {
        fillByChargeSpecies.template operator()<ChargeSpecies::Minus>();
      }
    }
  }

  template <DataMode DataModeValue, ParticleSpecies ParticleSpeciesValue>
    requires IsValid<DataModeValue, ParticleSpeciesValue>
  void fillCalculationYieldByParticleSpecies()
  {
    if (!groupAnalysis.cfgFlagsCalculationYield.value.get(toI(ParticleSpeciesValue))) {
      return;
    }

    const std::int32_t chargeSign{[&] {
      if constexpr (DataModeValue == DataMode::McMcParticle) {
        return holderMcParticle.charge;
      } else {
        return holderTrack.sign;
      }
    }()};
    if (chargeSign == 0) {
      return;
    }

    const auto fillByChargeSpecies{
      [&]<ChargeSpecies ChargeSpeciesValue>
        requires IsValid<ChargeSpeciesValue>
      () {
        if constexpr (DataModeValue == DataMode::McMcParticle) {
          if (isPid<getValue<ParticleSpeciesAll>(ParticleSpeciesValue), ChargeSpeciesValue>()) {
            hrCalculationYield.fill(C_CS("CalculationYield/hVzCentralityPtMcEtaMc_mc") + C_SV(getName(ParticleSpeciesValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.vz, holderEvent.centrality, holderMcParticle.pt, holderMcParticle.eta);
          }
        } else {
          const auto fillByPidStrategy{
            [&]<PidStrategy PidStrategyValue>
              requires IsValid<PidStrategyValue>
            () {
              if constexpr (DataModeValue == DataMode::McTrack) {
                if (isPid<getValue<ParticleSpeciesAll>(ParticleSpeciesValue), ChargeSpeciesValue>() && isPid<getValue<PidStrategyAll>(PidStrategyValue), getValue<ParticleSpeciesAll>(ParticleSpeciesValue)>(groupTrack.cfgFlagRejectionOthers.value)) {
                  if (groupTrack.cfgFlagMcParticleMomentum.value) {
                    hrCalculationYield.fill(C_CS("CalculationYield/hVzCentralityPtMcEtaMc_mc") + C_SV(getName(PidStrategyValue)) + C_SV(getName(ParticleSpeciesValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.vz, holderEvent.centrality, holderMcParticle.pt, holderMcParticle.eta);
                  } else {
                    hrCalculationYield.fill(C_CS("CalculationYield/hVzCentralityPtEta_mc") + C_SV(getName(PidStrategyValue)) + C_SV(getName(ParticleSpeciesValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.vz, holderEvent.centrality, holderTrack.pt, holderTrack.eta);
                  }
                }
              } else { // DataModeValue == DataMode::RawTrack
                if (isPid<getValue<PidStrategyAll>(PidStrategyValue), getValue<ParticleSpeciesAll>(ParticleSpeciesValue)>(groupTrack.cfgFlagRejectionOthers.value)) {
                  hrCalculationYield.fill(C_CS("CalculationYield/hVzCentralityPtEta_") + C_SV(getName<NameKind::Lower>(PidStrategyValue)) + C_SV(getName(ParticleSpeciesValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.vz, holderEvent.centrality, holderTrack.pt, holderTrack.eta);
                }
              }
            }};

          fillByPidStrategy.template operator()<PidStrategy::Tpc>();
          fillByPidStrategy.template operator()<PidStrategy::TpcTof>();
        }
      }};

    if (chargeSign > 0) {
      fillByChargeSpecies.template operator()<ChargeSpecies::Plus>();
    } else {
      fillByChargeSpecies.template operator()<ChargeSpecies::Minus>();
    }
  }

  template <ParticleSpecies ParticleSpeciesValue>
    requires IsValid<ParticleSpeciesValue>
  void fillCalculationPurityByParticleSpecies()
  {
    if (!groupAnalysis.cfgFlagsCalculationPurity.value.get(toI(ParticleSpeciesValue)) || holderTrack.sign == 0) {
      return;
    }

    const auto fillByChargeSpecies{
      [&]<ChargeSpecies ChargeSpeciesValue>
        requires IsValid<ChargeSpeciesValue>
      () {
        const auto fillByPidStrategy{
          [&]<PidStrategy PidStrategyValue>
            requires IsValid<PidStrategyValue>
          () {
            if (isPid<getValue<PidStrategyAll>(PidStrategyValue), getValue<ParticleSpeciesAll>(ParticleSpeciesValue)>(groupTrack.cfgFlagRejectionOthers.value)) {
              hrCalculationPurity.fill(C_CS("CalculationPurity/pCentralityPtEtaPurity") + C_SV(getName(PidStrategyValue)) + C_SV(getName(ParticleSpeciesValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.centrality, holderTrack.pt, holderTrack.eta, isPid<getValue<ParticleSpeciesAll>(ParticleSpeciesValue), ChargeSpeciesValue>() ? 1. : 0.);
            }
          }};

        fillByPidStrategy.template operator()<PidStrategy::Tpc>();
        fillByPidStrategy.template operator()<PidStrategy::TpcTof>();
      }};

    if (holderTrack.sign > 0) {
      fillByChargeSpecies.template operator()<ChargeSpecies::Plus>();
    } else {
      fillByChargeSpecies.template operator()<ChargeSpecies::Minus>();
    }
  }

  template <ParticleSpecies ParticleSpeciesValue, typename MP>
    requires IsValid<ParticleSpeciesValue>
  void fillCalculationFractionPrimaryByParticleSpecies(const MP& mcParticle)
  {
    if (!groupAnalysis.cfgFlagsCalculationFractionPrimary.value.get(toI(ParticleSpeciesValue)) || holderTrack.sign == 0) {
      return;
    }

    const auto fillByChargeSpecies{
      [&]<ChargeSpecies ChargeSpeciesValue>
        requires IsValid<ChargeSpeciesValue>
      () {
        const auto fillByPidStrategy{
          [&]<PidStrategy PidStrategyValue>
            requires IsValid<PidStrategyValue>
          () {
            if (isPid<getValue<ParticleSpeciesAll>(ParticleSpeciesValue), ChargeSpeciesValue>() && isPid<getValue<PidStrategyAll>(PidStrategyValue), getValue<ParticleSpeciesAll>(ParticleSpeciesValue)>(groupTrack.cfgFlagRejectionOthers.value)) {
              hrCalculationFractionPrimary.fill(C_CS("CalculationFractionPrimary/pCentralityPtEtaFractionPrimary") + C_SV(getName(PidStrategyValue)) + C_SV(getName(ParticleSpeciesValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.centrality, holderTrack.pt, holderTrack.eta, mcParticle.isPhysicalPrimary() ? 1. : 0.);
            }
          }};

        fillByPidStrategy.template operator()<PidStrategy::Tpc>();
        fillByPidStrategy.template operator()<PidStrategy::TpcTof>();
      }};

    if (holderTrack.sign > 0) {
      fillByChargeSpecies.template operator()<ChargeSpecies::Plus>();
    } else {
      fillByChargeSpecies.template operator()<ChargeSpecies::Minus>();
    }
  }

  void initCalculationFluctuation()
  {
    for (std::int32_t const& iParticleNumber : std::views::iota(0, NEs<ParticleNumber>)) {
      if (static_cast<bool>(groupAnalysis.cfgFlagsCalculationFluctuation.value.get(iParticleNumber))) {
        for (std::int32_t const& iChargeNumber : std::views::iota(0, NEs<ChargeNumber>)) {
          fluctuationCalculatorTrack[iParticleNumber][iChargeNumber]->init();
        }
      }
    }
  }

  template <DataMode DataModeValue, ParticleNumber ParticleNumberValue>
    requires IsValid<DataModeValue, ParticleNumberValue>
  void calculateFluctuationByParticleNumber()
  {
    if (!groupAnalysis.cfgFlagsCalculationFluctuation.value.get(toI(ParticleNumberValue))) {
      return;
    }

    const std::int32_t chargeSign{[&] {
      if constexpr (DataModeValue == DataMode::McMcParticle) {
        return holderMcParticle.charge;
      } else {
        return holderTrack.sign;
      }
    }()};
    if (chargeSign == 0) {
      return;
    }

    const bool doUsingMcParticleMomentum{[&] {
      if constexpr (DataModeValue == DataMode::McMcParticle) {
        return true;
      } else if constexpr (DataModeValue == DataMode::McTrack) {
        return groupTrack.cfgFlagMcParticleMomentum.value;
      } else { // DataModeValue == DataMode::RawTrack
        return false;
      }
    }()};
    if (!isGoodMomentum(doUsingMcParticleMomentum)) {
      return;
    }

    if constexpr (DataModeValue == DataMode::McMcParticle) {
      ++holderDerivedData.nMcParticles[toI(chargeSign > 0 ? ChargeSpecies::Plus : ChargeSpecies::Minus)];
    } else {
      ++holderDerivedData.nTracks[toI(chargeSign > 0 ? ChargeSpecies::Plus : ChargeSpecies::Minus)];
    }

    const auto calculateByParticleSpecies{
      [&]<ParticleSpecies ParticleSpeciesValue>
        requires IsValid<ParticleSpeciesValue> && (ParticleNumberValue == ParticleNumber::Charge || (ParticleNumberValue == ParticleNumber::Kaon && ParticleSpeciesValue == ParticleSpecies::Kaon) || (ParticleNumberValue == ParticleNumber::Proton && ParticleSpeciesValue == ParticleSpecies::Proton))
      () {
        const auto calculateByChargeSpecies{
          [&]<ChargeSpecies ChargeSpeciesValue>
            requires IsValid<ChargeSpeciesValue>
          () {
            if constexpr (DataModeValue != DataMode::RawTrack) {
              if (!isPid<getValue<ParticleSpeciesAll>(ParticleSpeciesValue), ChargeSpeciesValue>()) {
                return;
              }
            }

            const bool doUsingTofPid{[&] {
              return [&] {
                if constexpr (DataModeValue == DataMode::McMcParticle) {
                  return holderMcParticle.pt;
                } else if constexpr (DataModeValue == DataMode::McTrack) {
                  return groupTrack.cfgFlagMcParticleMomentum.value ? holderMcParticle.pt : holderTrack.pt;
                } else { // DataModeValue == DataMode::RawTrack
                  return holderTrack.pt;
                }
              }() >= groupTrack.cfgThresholdsPtTofPid.value.get(toI(ParticleSpeciesValue));
            }()};
            if constexpr (DataModeValue != DataMode::McMcParticle) {
              if (!(doUsingTofPid ? isPid<getValue<PidStrategyAll>(PidStrategy::TpcTof), getValue<ParticleSpeciesAll>(ParticleSpeciesValue)>(groupTrack.cfgFlagRejectionOthers.value) : isPid<getValue<PidStrategyAll>(PidStrategy::Tpc), getValue<ParticleSpeciesAll>(ParticleSpeciesValue)>(groupTrack.cfgFlagRejectionOthers.value))) {
                return;
              }
            }

            const double efficiency{doUsingTofPid ? getEfficiency<PidStrategy::TpcTof, ParticleSpeciesValue, ChargeSpeciesValue>(doUsingMcParticleMomentum) : getEfficiency<PidStrategy::Tpc, ParticleSpeciesValue, ChargeSpeciesValue>(doUsingMcParticleMomentum)}; // NOLINT(clang-analyzer-core.NullDereference)
            const auto fill{
              [&] {
                if constexpr (ChargeSpeciesValue == ChargeSpecies::Plus) {
                  fluctuationCalculatorTrack[toI(ParticleNumberValue)][toI(ChargeNumber::Plus)]->fill(1., efficiency);
                  fluctuationCalculatorTrack[toI(ParticleNumberValue)][toI(ChargeNumber::Net)]->fill(1., efficiency);
                } else {
                  fluctuationCalculatorTrack[toI(ParticleNumberValue)][toI(ChargeNumber::Minus)]->fill(1., efficiency);
                  fluctuationCalculatorTrack[toI(ParticleNumberValue)][toI(ChargeNumber::Net)]->fill(-1., efficiency);
                }
                fluctuationCalculatorTrack[toI(ParticleNumberValue)][toI(ChargeNumber::Total)]->fill(1., efficiency);
              }};
            if constexpr (DataModeValue == DataMode::McMcParticle) {
              ++holderMcEvent.numbers[toI(ParticleNumberValue)][toI(ChargeSpeciesValue)];
              if (gRandom->Rndm() < efficiency) {
                ++holderMcEvent.numbersEff[toI(ParticleNumberValue)][toI(ChargeSpeciesValue)];
                fill();
              }
              holderDerivedData.signedEfficienciesMcParticle.push_back(HolderDerivedData::convertRound<aod::mini_mc_particle::SignedEfficiency::type>(std::copysign(std::numeric_limits<aod::mini_mc_particle::SignedEfficiency::type>::max(), chargeSign) * efficiency));
            } else {
              ++holderEvent.numbers[toI(ParticleNumberValue)][toI(ChargeSpeciesValue)];
              fill();
              holderDerivedData.signedEfficienciesTrack.push_back(HolderDerivedData::convertRound<aod::mini_track::SignedEfficiency::type>(std::copysign(std::numeric_limits<aod::mini_track::SignedEfficiency::type>::max(), chargeSign) * efficiency));
            }
          }};

        if (chargeSign > 0) { // NOLINT(clang-analyzer-core.NullDereference)
          calculateByChargeSpecies.template operator()<ChargeSpecies::Plus>();
        } else {
          calculateByChargeSpecies.template operator()<ChargeSpecies::Minus>();
        }
      }};

    if constexpr (ParticleNumberValue == ParticleNumber::Kaon) {
      calculateByParticleSpecies.template operator()<ParticleSpecies::Kaon>();
    } else if constexpr (ParticleNumberValue == ParticleNumber::Proton) {
      calculateByParticleSpecies.template operator()<ParticleSpecies::Proton>();
    } else { // ParticleNumberValue == ParticleNumber::Charge
      calculateByParticleSpecies.template operator()<ParticleSpecies::Pion>();
      calculateByParticleSpecies.template operator()<ParticleSpecies::Kaon>();
      calculateByParticleSpecies.template operator()<ParticleSpecies::Proton>();
    }
  }

  template <DataMode DataModeValue, ParticleNumber ParticleNumberValue>
    requires IsValid<DataModeValue, ParticleNumberValue>
  void fillCalculationFluctuationByParticleNumber()
  {
    if (!groupAnalysis.cfgFlagsCalculationFluctuation.value.get(toI(ParticleNumberValue))) {
      return;
    }

    if constexpr (DataModeValue == DataMode::McMcParticle) {
      hrCalculationFluctuation.fill(C_CS("CalculationFluctuation/hCentralityN") + C_SV(getName(ParticleNumberValue)) + C_SV(getName(ChargeSpecies::Plus)) + C_CS("N") + C_SV(getName(ParticleNumberValue)) + C_SV(getName(ChargeSpecies::Minus)) + C_CS("_mc"), holderEvent.centrality, holderMcEvent.numbers[toI(ParticleNumberValue)][toI(ChargeSpecies::Plus)], holderMcEvent.numbers[toI(ParticleNumberValue)][toI(ChargeSpecies::Minus)]);
      hrCalculationFluctuation.fill(C_CS("CalculationFluctuation/hCentralityN") + C_SV(getName(ParticleNumberValue)) + C_SV(getName(ChargeSpecies::Plus)) + C_CS("N") + C_SV(getName(ParticleNumberValue)) + C_SV(getName(ChargeSpecies::Minus)) + C_CS("_mcEff"), holderEvent.centrality, holderMcEvent.numbersEff[toI(ParticleNumberValue)][toI(ChargeSpecies::Plus)], holderMcEvent.numbersEff[toI(ParticleNumberValue)][toI(ChargeSpecies::Minus)]);
    } else {
      hrCalculationFluctuation.fill(C_CS("CalculationFluctuation/hCentralityN") + C_SV(getName(ParticleNumberValue)) + C_SV(getName(ChargeSpecies::Plus)) + C_CS("N") + C_SV(getName(ParticleNumberValue)) + C_SV(getName(ChargeSpecies::Minus)), holderEvent.centrality, holderEvent.numbers[toI(ParticleNumberValue)][toI(ChargeSpecies::Plus)], holderEvent.numbers[toI(ParticleNumberValue)][toI(ChargeSpecies::Minus)]);
    }

    const auto fillByChargeNumber{
      [&]<ChargeNumber ChargeNumberValue>
        requires IsValid<ChargeNumberValue>
      () {
        for (std::int32_t const& iOrderKey : std::views::iota(0, fluctuation_calculator_base::NOrderKeys)) {
          if constexpr (DataModeValue == DataMode::McMcParticle) {
            hrCalculationFluctuation.fill(C_CS("CalculationFluctuation/hFluctuationCalculator") + C_SV(getName(ParticleNumberValue)) + C_SV(getName(ChargeNumberValue)) + C_CS("_mc"), holderEvent.centrality, holderEvent.subgroupIndex, iOrderKey, fluctuationCalculatorTrack[toI(ParticleNumberValue)][toI(ChargeNumberValue)]->getProduct(iOrderKey));
          } else {
            hrCalculationFluctuation.fill(C_CS("CalculationFluctuation/hFluctuationCalculator") + C_SV(getName(ParticleNumberValue)) + C_SV(getName(ChargeNumberValue)), holderEvent.centrality, holderEvent.subgroupIndex, iOrderKey, fluctuationCalculatorTrack[toI(ParticleNumberValue)][toI(ChargeNumberValue)]->getProduct(iOrderKey));
          }
        }
      }};

    fillByChargeNumber.template operator()<ChargeNumber::Plus>();
    fillByChargeNumber.template operator()<ChargeNumber::Minus>();
    fillByChargeNumber.template operator()<ChargeNumber::Total>();
    fillByChargeNumber.template operator()<ChargeNumber::Net>();
  }

  template <bool DoInitingEvent, typename T, typename TIs>
  bool initTrack(const T& track, const TIs& tracksIu)
  {
    holderTrack.clear();
    holderTrack.dca[toI(DcaAxis::Xy)] = track.dcaXY();
    holderTrack.dca[toI(DcaAxis::Z)] = track.dcaZ();
    holderTrack.sign = track.sign();
    holderTrack.p = track.p();
    holderTrack.pt = track.pt();
    holderTrack.eta = track.eta();
    holderTrack.phi = track.phi();
    {
      const std::int64_t localIndexTrackIu{track.globalIndex() - static_cast<std::int64_t>(tracksIu.offset())};
      if (0 <= localIndexTrackIu && localIndexTrackIu < tracksIu.size()) {
        const auto& trackIu{tracksIu.iteratorAt(localIndexTrackIu)};
        if (track.globalIndex() == trackIu.globalIndex()) {
          holderTrack.phiIu = trackIu.phi();
        } else {
          LOG(warning) << "Mismatched track " << track.globalIndex() << " and trackIu " << trackIu.globalIndex();
        }
      } else {
        LOG(warning) << "Invalid trackIu " << track.globalIndex();
      }
    }
    holderTrack.hasPid[toI(Detector::Tpc)] = (track.hasTPC() && track.tpcSignal() > 0.);
    if (holderTrack.hasPid[toI(Detector::Tpc)]) {
      holderTrack.nSigmaPid[toI(PidStrategyAll::Tpc)] = {HolderTrack::truncateNSigmaPid(track.tpcNSigmaPi() - getShiftNSigmaPid<Detector::Tpc, ParticleSpecies::Pion>()), HolderTrack::truncateNSigmaPid(track.tpcNSigmaKa() - getShiftNSigmaPid<Detector::Tpc, ParticleSpecies::Kaon>()), HolderTrack::truncateNSigmaPid(track.tpcNSigmaPr() - getShiftNSigmaPid<Detector::Tpc, ParticleSpecies::Proton>())};
    }
    holderTrack.hasPid[toI(Detector::Tof)] = (track.hasTOF() && track.beta() > 0.);
    if (holderTrack.hasPid[toI(Detector::Tof)]) {
      holderTrack.nSigmaPid[toI(PidStrategyAll::Tof)] = {HolderTrack::truncateNSigmaPid(track.tofNSigmaPi() - getShiftNSigmaPid<Detector::Tof, ParticleSpecies::Pion>()), HolderTrack::truncateNSigmaPid(track.tofNSigmaKa() - getShiftNSigmaPid<Detector::Tof, ParticleSpecies::Kaon>()), HolderTrack::truncateNSigmaPid(track.tofNSigmaPr() - getShiftNSigmaPid<Detector::Tof, ParticleSpecies::Proton>())};
      if (holderTrack.hasPid[toI(Detector::Tpc)]) {
        for (std::int32_t const& iParticleSpecies : std::views::iota(0, NEs<ParticleSpecies>)) {
          holderTrack.nSigmaPid[toI(PidStrategyAll::TpcTofCombined)][iParticleSpecies] = HolderTrack::truncateNSigmaPid(std::copysign(std::hypot(holderTrack.nSigmaPid[toI(PidStrategyAll::Tpc)][iParticleSpecies], holderTrack.nSigmaPid[toI(PidStrategyAll::Tof)][iParticleSpecies]), holderTrack.nSigmaPid[toI(PidStrategyAll::Tpc)][iParticleSpecies] + holderTrack.nSigmaPid[toI(PidStrategyAll::Tof)][iParticleSpecies]));
        }
      }
    }

    if constexpr (DoInitingEvent) {
      if (track.isPrimaryTrack() && holderTrack.sign != 0) {
        const std::int32_t chargeSpeciesIndex{toI(holderTrack.sign > 0 ? ChargeSpecies::Plus : ChargeSpecies::Minus)};
        ++holderEvent.nGlobalTracks[chargeSpeciesIndex];
        if (track.isPVContributor()) {
          ++holderEvent.nPvContributors[chargeSpeciesIndex];
        }
        for (std::int32_t const& iDcaAxis : std::views::iota(0, NEs<DcaAxis>)) {
          holderEvent.dca[toI(DcaKind::Mean)][iDcaAxis][chargeSpeciesIndex] += holderTrack.dca[iDcaAxis];
          holderEvent.dca[toI(DcaKind::Sigma)][iDcaAxis][chargeSpeciesIndex] += std::pow(holderTrack.dca[iDcaAxis], 2.);
        }
        if (holderTrack.hasPid[toI(Detector::Tof)]) {
          ++holderEvent.nTofBeta[chargeSpeciesIndex];
        }
      }
    }

    if constexpr (DoInitingEvent) {
      if (groupAnalysis.cfgFlagQaRun.value && track.isPrimaryTrack()) {
        if (holderTrack.sign > 0) {
          fillQaRunByTrackByChargeSpecies<ChargeSpecies::Plus>(track);
        } else if (holderTrack.sign < 0) {
          fillQaRunByTrackByChargeSpecies<ChargeSpecies::Minus>(track);
        }
      }
    }

    if constexpr (!DoInitingEvent) {
      if (groupAnalysis.cfgFlagQaTrack.value && track.isPrimaryTrack()) {
        if (holderTrack.sign > 0) {
          fillQaTrackByChargeSpecies<ChargeSpecies::Plus>(track);
        } else if (holderTrack.sign < 0) {
          fillQaTrackByChargeSpecies<ChargeSpecies::Minus>(track);
        }
      }
    }

    if (!isGoodTrack(track)) {
      return false;
    }

    if constexpr (!DoInitingEvent) {
      if (groupAnalysis.cfgFlagQaDca.value) {
        if (holderTrack.sign > 0) {
          fillQaDcaByChargeSpecies<ChargeSpecies::Plus>();
        } else if (holderTrack.sign < 0) {
          fillQaDcaByChargeSpecies<ChargeSpecies::Minus>();
        }
      }
    }

    if (!isGoodDca()) {
      return false;
    }

    if constexpr (!DoInitingEvent) {
      if (isEnabled(groupAnalysis.cfgFlagsQaAcceptance) && (holderTrack.eta * holderEvent.vz > 0. && std::abs(holderEvent.vz) > groupEvent.cfgCutMaxAbsVz.value - 1.)) {
        fillQaAcceptancebyParticleSpeciesAll<ParticleSpeciesAll::All>(track);
        fillQaAcceptancebyParticleSpeciesAll<ParticleSpeciesAll::Pion>(track);
        fillQaAcceptancebyParticleSpeciesAll<ParticleSpeciesAll::Kaon>(track);
        fillQaAcceptancebyParticleSpeciesAll<ParticleSpeciesAll::Proton>(track);
      }
    }

    return true;
  }

  template <bool IsMc, typename MP>
  bool initMcParticle(const MP& mcParticle)
  {
    holderMcParticle.clear();
    holderMcParticle.pdgCode = mcParticle.pdgCode();
    const TParticlePDG* const particlePdg{pdg->GetParticle(mcParticle.pdgCode())};
    if (particlePdg) {
      holderMcParticle.charge = std::llrint(particlePdg->Charge());
    } else {
      switch (std::abs(holderMcParticle.pdgCode) / 100000000) {
        case 10:
          holderMcParticle.charge = holderMcParticle.pdgCode / 10000 % 1000;
          break;
        default:
          break;
      }
    }
    holderMcParticle.pt = mcParticle.pt();
    holderMcParticle.eta = mcParticle.eta();

    return isGoodMcParticle<IsMc>(mcParticle);
  }

  template <typename C, typename Ts, typename TIs>
  bool initEvent(const C& collision, const Ts& tracks, const TIs& tracksIu)
  {
    holderEvent.clear();
    holderEvent.vz = collision.posZ();
    switch (groupEvent.cfgIndexDefinitionCentrality.value) {
      case toI(CentralityDefinition::Ft0a):
        holderEvent.centrality = collision.centFT0A();
        break;
      case toI(CentralityDefinition::Ft0c):
        holderEvent.centrality = collision.centFT0C();
        break;
      default:
        holderEvent.centrality = collision.centFT0M();
        break;
    }

    hrCounter.fill(C_CS("hNEvents"), 0.);
    if (groupAnalysis.cfgFlagQaCentrality.value) {
      hrQaCentrality.fill(C_CS("QaCentrality/hCentralitySelection"), holderEvent.centrality, 0.);
    }

    if (!collision.has_foundBC()) {
      hrCounter.fill(C_CS("hNEvents"), 2.);
      if (groupAnalysis.cfgFlagQaCentrality.value) {
        hrQaCentrality.fill(C_CS("QaCentrality/hCentralitySelection"), holderEvent.centrality, 2.);
      }
      return false;
    }

    const auto& bc{collision.template bc_as<aod::BCsWithTimestamps>()};
    holderEvent.runNumber = bc.runNumber();

    if (!holderCcdb.runNumbersIndicesGroupIndices.contains(holderEvent.runNumber)) {
      hrCounter.fill(C_CS("hNEvents"), 2.);
      if (groupAnalysis.cfgFlagQaCentrality.value) {
        hrQaCentrality.fill(C_CS("QaCentrality/hCentralitySelection"), holderEvent.centrality, 2.);
      }
      return false;
    }

    std::tie(holderEvent.runIndex, holderEvent.runGroupIndex) = holderCcdb.runNumbersIndicesGroupIndices.at(holderEvent.runNumber);

    if (holderEvent.runGroupIndex == 0 || (groupEvent.cfgFlagRejectionRunBad.value && holderEvent.runGroupIndex < 0)) {
      hrCounter.fill(C_CS("hNEvents"), 2.);
      if (groupAnalysis.cfgFlagQaCentrality.value) {
        hrQaCentrality.fill(C_CS("QaCentrality/hCentralitySelection"), holderEvent.centrality, 2.);
      }
      return false;
    }

    if (!rctFlagsChecker.checkTable(collision)) {
      hrCounter.fill(C_CS("hNEvents"), 3.);
      if (groupAnalysis.cfgFlagQaCentrality.value) {
        hrQaCentrality.fill(C_CS("QaCentrality/hCentralitySelection"), holderEvent.centrality, 3.);
      }
      return false;
    }

    for (std::int32_t const& iEvSel : std::views::iota(0, aod::evsel::EventSelectionFlags::kNsel)) {
      if (((groupEvent.cfgBitsSelectionEvent.value >> iEvSel) & 1) && !collision.selection_bit(iEvSel)) {
        hrCounter.fill(C_CS("hNEvents"), 4.);
        hrCounter.fill(C_CS("hNEvents"), 10. + iEvSel);
        if (groupAnalysis.cfgFlagQaCentrality.value) {
          hrQaCentrality.fill(C_CS("QaCentrality/hCentralitySelection"), holderEvent.centrality, 4.);
          hrQaCentrality.fill(C_CS("QaCentrality/hCentralitySelection"), holderEvent.centrality, 10. + iEvSel);
        }
        return false;
      }
    }

    if (groupEvent.cfgFlagInelEvent.value && !collision.isInelGt0()) {
      hrCounter.fill(C_CS("hNEvents"), 5.);
      if (groupAnalysis.cfgFlagQaCentrality.value) {
        hrQaCentrality.fill(C_CS("QaCentrality/hCentralitySelection"), holderEvent.centrality, 5.);
      }
      return false;
    }

    if (groupAnalysis.cfgFlagQaEvent.value) {
      hrQaEvent.fill(C_CS("QaEvent/hVxVy"), collision.posX(), collision.posY());
      hrQaEvent.fill(C_CS("QaEvent/hVz"), holderEvent.vz);
    }

    if (!(std::abs(holderEvent.vz) < groupEvent.cfgCutMaxAbsVz.value)) {
      hrCounter.fill(C_CS("hNEvents"), 6.);
      if (groupAnalysis.cfgFlagQaCentrality.value) {
        hrQaCentrality.fill(C_CS("QaCentrality/hCentralitySelection"), holderEvent.centrality, 6.);
      }
      return false;
    }

    if (groupAnalysis.cfgFlagQaRun.value) {
      hrQaRun.fill(C_CS("QaRun/pRunIndexVx"), holderEvent.runIndex, collision.posX());
      hrQaRun.fill(C_CS("QaRun/pRunIndexVy"), holderEvent.runIndex, collision.posY());
      hrQaRun.fill(C_CS("QaRun/pRunIndexVz"), holderEvent.runIndex, holderEvent.vz);
      hrQaRun.fill(C_CS("QaRun/pRunIndexMultiplicityFt0a"), holderEvent.runIndex, collision.multZeqFT0A());
      hrQaRun.fill(C_CS("QaRun/pRunIndexMultiplicityFt0c"), holderEvent.runIndex, collision.multZeqFT0C());
      if (HolderEvent::RangeCentrality.first <= collision.centFT0A() && collision.centFT0A() <= HolderEvent::RangeCentrality.second) {
        hrQaRun.fill(C_CS("QaRun/pRunIndexCentralityFt0a"), holderEvent.runIndex, collision.centFT0A());
      }
      if (HolderEvent::RangeCentrality.first <= collision.centFT0C() && collision.centFT0C() <= HolderEvent::RangeCentrality.second) {
        hrQaRun.fill(C_CS("QaRun/pRunIndexCentralityFt0c"), holderEvent.runIndex, collision.centFT0C());
      }
      if (HolderEvent::RangeCentrality.first <= collision.centFT0M() && collision.centFT0M() <= HolderEvent::RangeCentrality.second) {
        hrQaRun.fill(C_CS("QaRun/pRunIndexCentralityFt0m"), holderEvent.runIndex, collision.centFT0M());
      }
    }

    for (const auto& track : tracks) {
      if (!track.has_collision() || track.collisionId() != collision.globalIndex()) {
        continue;
      }

      initTrack<true>(track, tracksIu);
    }
    for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
      if (holderEvent.nGlobalTracks[iChargeSpecies] > 0) {
        for (std::int32_t const& iDcaAxis : std::views::iota(0, NEs<DcaAxis>)) {
          holderEvent.dca[toI(DcaKind::Mean)][iDcaAxis][iChargeSpecies] /= holderEvent.nGlobalTracks[iChargeSpecies];
          holderEvent.dca[toI(DcaKind::Sigma)][iDcaAxis][iChargeSpecies] = std::sqrt(std::max(0., holderEvent.dca[toI(DcaKind::Sigma)][iDcaAxis][iChargeSpecies] / holderEvent.nGlobalTracks[iChargeSpecies] - std::pow(holderEvent.dca[toI(DcaKind::Mean)][iDcaAxis][iChargeSpecies], 2.)));
        }
      }
    }

    if (groupAnalysis.cfgFlagQaRun.value) {
      fillQaRunByEventByChargeSpecies<ChargeSpecies::Plus>();
      fillQaRunByEventByChargeSpecies<ChargeSpecies::Minus>();
    }

    if (groupAnalysis.cfgFlagQaEvent.value) {
      hrQaEvent.fill(C_CS("QaEvent/hNPvContributorsNGlobalTracks"), holderEvent.getNPvContributors(), holderEvent.getNGlobalTracks());
      if (holderEvent.getNGlobalTracks() > 0) {
        hrQaEvent.fill(C_CS("QaEvent/hNGlobalTracks") + C_SV(getName(DcaKind::Mean)) + C_CS("Dca") + C_SV(getName(DcaAxis::Xy)), holderEvent.getNGlobalTracks(), holderEvent.getDca<DcaKind::Mean, DcaAxis::Xy>());
        hrQaEvent.fill(C_CS("QaEvent/hNGlobalTracks") + C_SV(getName(DcaKind::Mean)) + C_CS("Dca") + C_SV(getName(DcaAxis::Z)), holderEvent.getNGlobalTracks(), holderEvent.getDca<DcaKind::Mean, DcaAxis::Z>());
      }
      hrQaEvent.fill(C_CS("QaEvent/hNTofBetaNGlobalTracks"), holderEvent.getNTofBeta(), holderEvent.getNGlobalTracks());
    }

    if (!(holderEvent.getNPvContributors() - holderEvent.getNGlobalTracks() > groupEvent.cfgCutMinDeviationNPvContributors.value)) {
      hrCounter.fill(C_CS("hNEvents"), 7.);
      if (groupAnalysis.cfgFlagQaCentrality.value) {
        hrQaCentrality.fill(C_CS("QaCentrality/hCentralitySelection"), holderEvent.centrality, 7.);
      }
      return false;
    }

    hrCounter.fill(C_CS("hNEvents"), 1.);
    if (groupAnalysis.cfgFlagQaCentrality.value) {
      hrQaCentrality.fill(C_CS("QaCentrality/hCentralitySelection"), holderEvent.centrality, 1.);
    }

    if (groupAnalysis.cfgFlagQaEvent.value) {
      if (holderEvent.getNGlobalTracks() > 0) {
        hrQaEvent.fill(C_CS("QaEvent/hNGlobalTracks") + C_SV(getName(DcaKind::Mean)) + C_CS("Dca") + C_SV(getName(DcaAxis::Xy)) + C_CS("_nPvContributorsCut"), holderEvent.getNGlobalTracks(), holderEvent.getDca<DcaKind::Mean, DcaAxis::Xy>());
        hrQaEvent.fill(C_CS("QaEvent/hNGlobalTracks") + C_SV(getName(DcaKind::Mean)) + C_CS("Dca") + C_SV(getName(DcaAxis::Z)) + C_CS("_nPvContributorsCut"), holderEvent.getNGlobalTracks(), holderEvent.getDca<DcaKind::Mean, DcaAxis::Z>());
      }
      hrQaEvent.fill(C_CS("QaEvent/hNTofBetaNGlobalTracks_nPvContributorsCut"), holderEvent.getNTofBeta(), holderEvent.getNGlobalTracks());
    }

    if (groupAnalysis.cfgFlagQaCentrality.value) {
      hrQaCentrality.fill(C_CS("QaCentrality/hCentralityMultiplicity"), holderEvent.centrality, collision.multNTracksPVeta1());
    }

    return true;
  }

  template <typename MC>
  bool initMcEvent(const MC& mcCollision)
  {
    holderMcEvent.clear();
    holderMcEvent.vz = mcCollision.posZ();

    hrCounter.fill(C_CS("hNMcEvents"), 0.);

    if (!mcCollision.has_bc()) {
      hrCounter.fill(C_CS("hNMcEvents"), 2.);
      return false;
    }

    const auto& bc{mcCollision.template bc_as<aod::BCsWithTimestamps>()};
    holderMcEvent.runNumber = bc.runNumber();

    if (!holderCcdb.runNumbersIndicesGroupIndices.contains(holderMcEvent.runNumber)) {
      hrCounter.fill(C_CS("hNMcEvents"), 2.);
      return false;
    }

    std::tie(holderMcEvent.runIndex, holderMcEvent.runGroupIndex) = holderCcdb.runNumbersIndicesGroupIndices.at(holderMcEvent.runNumber);

    if (holderMcEvent.runGroupIndex == 0 || (groupEvent.cfgFlagRejectionRunBad.value && holderMcEvent.runGroupIndex < 0)) {
      hrCounter.fill(C_CS("hNMcEvents"), 2.);
      return false;
    }

    if (groupEvent.cfgFlagInelEventMc.value && !mcCollision.isInelGt0()) {
      hrCounter.fill(C_CS("hNMcEvents"), 3.);
      return false;
    }

    if (groupEvent.cfgFlagCutVzMc.value && !(std::abs(holderMcEvent.vz) < groupEvent.cfgCutMaxAbsVz.value)) {
      hrCounter.fill(C_CS("hNMcEvents"), 4.);
      return false;
    }

    hrCounter.fill(C_CS("hNMcEvents"), 1.);

    return true;
  }

  void processRaw(const soa::Filtered<aod::JoinedCollisions>::iterator& collision, const soa::Filtered<aod::JoinedTracks>& tracks, const aod::TracksIU& tracksIu, const aod::BCsWithTimestamps&)
  {
    if (!initEvent(collision, tracks, tracksIu) || (!groupAnalysis.cfgFlagQaTrack.value && !groupAnalysis.cfgFlagQaDca.value && !isEnabled(groupAnalysis.cfgFlagsQaAcceptance) && !isEnabled(groupAnalysis.cfgFlagsQaPhi) && !isEnabled(groupAnalysis.cfgFlagsQaPid) && !isEnabled(groupAnalysis.cfgFlagsCalculationYield) && !isEnabled(groupAnalysis.cfgFlagsCalculationFluctuation))) {
      return;
    }

    if (isEnabled(groupAnalysis.cfgFlagsCalculationFluctuation)) {
      holderEvent.subgroupIndex = gRandom->Integer(groupEvent.cfgNSubgroups.value);
      initCalculationFluctuation();
      holderDerivedData.clear();
    }

    for (const auto& track : tracks) {
      if (!track.has_collision() || track.collisionId() != collision.globalIndex()) {
        continue;
      }

      if (!initTrack<false>(track, tracksIu)) {
        continue;
      }

      if (isEnabled(groupAnalysis.cfgFlagsQaPhi)) {
        fillQaPhiByParticleSpeciesAll<DataMode::RawTrack, ParticleSpeciesAll::All>();
        fillQaPhiByParticleSpeciesAll<DataMode::RawTrack, ParticleSpeciesAll::Pion>();
        fillQaPhiByParticleSpeciesAll<DataMode::RawTrack, ParticleSpeciesAll::Kaon>();
        fillQaPhiByParticleSpeciesAll<DataMode::RawTrack, ParticleSpeciesAll::Proton>();
      }

      if (isEnabled(groupAnalysis.cfgFlagsQaPid)) {
        fillQaPidByParticleSpeciesAll<DataMode::RawTrack, ParticleSpeciesAll::All>(track);
        fillQaPidByParticleSpeciesAll<DataMode::RawTrack, ParticleSpeciesAll::Pion>(track);
        fillQaPidByParticleSpeciesAll<DataMode::RawTrack, ParticleSpeciesAll::Kaon>(track);
        fillQaPidByParticleSpeciesAll<DataMode::RawTrack, ParticleSpeciesAll::Proton>(track);
      }

      if (isEnabled(groupAnalysis.cfgFlagsCalculationYield)) {
        fillCalculationYieldByParticleSpecies<DataMode::RawTrack, ParticleSpecies::Pion>();
        fillCalculationYieldByParticleSpecies<DataMode::RawTrack, ParticleSpecies::Kaon>();
        fillCalculationYieldByParticleSpecies<DataMode::RawTrack, ParticleSpecies::Proton>();
      }

      if (isEnabled(groupAnalysis.cfgFlagsCalculationFluctuation)) {
        calculateFluctuationByParticleNumber<DataMode::RawTrack, ParticleNumber::Charge>();
        calculateFluctuationByParticleNumber<DataMode::RawTrack, ParticleNumber::Kaon>();
        calculateFluctuationByParticleNumber<DataMode::RawTrack, ParticleNumber::Proton>();
      }
    }

    if (isEnabled(groupAnalysis.cfgFlagsCalculationFluctuation)) {
      fillCalculationFluctuationByParticleNumber<DataMode::RawTrack, ParticleNumber::Charge>();
      fillCalculationFluctuationByParticleNumber<DataMode::RawTrack, ParticleNumber::Kaon>();
      fillCalculationFluctuationByParticleNumber<DataMode::RawTrack, ParticleNumber::Proton>();
      miniCollision(HolderDerivedData::convertFloor<aod::mini_collision::Vz::type>(holderEvent.vz * 10.), HolderDerivedData::convertFloor<aod::mini_collision::Centrality::type>(holderEvent.centrality * 500.), holderDerivedData.nTracks[toI(ChargeSpecies::Plus)], holderDerivedData.nTracks[toI(ChargeSpecies::Minus)]);
      for (auto const& signedEfficiency : holderDerivedData.signedEfficienciesTrack) {
        miniTrack(miniCollision.lastIndex(), signedEfficiency);
      }
    }
  }
  PROCESS_SWITCH(PartNumFluc, processRaw, "Process raw data", true);

  void processMc(const soa::Filtered<aod::JoinedMcCollisions>::iterator& mcCollision, const aod::McParticles& mcParticles, const soa::SmallGroups<aod::JoinedCollisionsWithMc>& collisions, const soa::Filtered<aod::JoinedTracksWithMc>& tracksUngrouped, const aod::TracksIU& tracksIuUngrouped, const aod::BCsWithTimestamps&)
  {
    if (!initMcEvent(mcCollision)) {
      return;
    }

    for (const auto& collision : collisions) {
      if (collision.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }

      const auto& tracks{tracksUngrouped.sliceBy(presliceTracksPerCollision, collision.globalIndex())};
      const auto& tracksIu{tracksIuUngrouped.sliceBy(presliceTracksPerCollision, collision.globalIndex())};

      if (!initEvent(collision, tracks, tracksIu)) {
        continue;
      }

      if (groupAnalysis.cfgFlagQaMc.value) {
        hrQaMc.fill(C_CS("QaMc/hCentralityVzMcDeltaVz"), holderEvent.centrality, holderMcEvent.vz, holderEvent.vz - holderMcEvent.vz);
      }

      if (isEnabled(groupAnalysis.cfgFlagsCalculationYield) || isEnabled(groupAnalysis.cfgFlagsCalculationFluctuation)) {
        if (isEnabled(groupAnalysis.cfgFlagsCalculationFluctuation)) {
          holderEvent.subgroupIndex = gRandom->Integer(groupEvent.cfgNSubgroups.value);
          initCalculationFluctuation();
          holderDerivedData.clear();
        }

        for (const auto& mcParticle : mcParticles) {
          if (!mcParticle.has_mcCollision() || mcParticle.mcCollisionId() != mcCollision.globalIndex()) {
            continue;
          }

          if (!initMcParticle<true>(mcParticle)) {
            continue;
          }

          if (isEnabled(groupAnalysis.cfgFlagsCalculationYield)) {
            fillCalculationYieldByParticleSpecies<DataMode::McMcParticle, ParticleSpecies::Pion>();
            fillCalculationYieldByParticleSpecies<DataMode::McMcParticle, ParticleSpecies::Kaon>();
            fillCalculationYieldByParticleSpecies<DataMode::McMcParticle, ParticleSpecies::Proton>();
          }

          if (isEnabled(groupAnalysis.cfgFlagsCalculationFluctuation)) {
            calculateFluctuationByParticleNumber<DataMode::McMcParticle, ParticleNumber::Charge>();
            calculateFluctuationByParticleNumber<DataMode::McMcParticle, ParticleNumber::Kaon>();
            calculateFluctuationByParticleNumber<DataMode::McMcParticle, ParticleNumber::Proton>();
          }
        }

        if (isEnabled(groupAnalysis.cfgFlagsCalculationFluctuation)) {
          fillCalculationFluctuationByParticleNumber<DataMode::McMcParticle, ParticleNumber::Charge>();
          fillCalculationFluctuationByParticleNumber<DataMode::McMcParticle, ParticleNumber::Kaon>();
          fillCalculationFluctuationByParticleNumber<DataMode::McMcParticle, ParticleNumber::Proton>();
        }
      }

      if (groupAnalysis.cfgFlagQaTrack.value || groupAnalysis.cfgFlagQaDca.value || isEnabled(groupAnalysis.cfgFlagsQaAcceptance) || isEnabled(groupAnalysis.cfgFlagsQaPhi) || isEnabled(groupAnalysis.cfgFlagsQaPid) || isEnabled(groupAnalysis.cfgFlagsCalculationYield) || isEnabled(groupAnalysis.cfgFlagsCalculationPurity) || isEnabled(groupAnalysis.cfgFlagsCalculationFractionPrimary) || isEnabled(groupAnalysis.cfgFlagsCalculationFluctuation)) {
        if (isEnabled(groupAnalysis.cfgFlagsCalculationFluctuation)) {
          initCalculationFluctuation();
        }

        for (const auto& track : tracks) {
          if (!track.has_collision() || track.collisionId() != collision.globalIndex() || !track.has_mcParticle()) {
            continue;
          }

          const auto& mcParticle{track.template mcParticle_as<aod::McParticles>()};
          if (!mcParticle.has_mcCollision() || mcParticle.mcCollisionId() != mcCollision.globalIndex()) {
            continue;
          }

          if (!initTrack<false>(track, tracksIu) || !initMcParticle<false>(mcParticle)) {
            continue;
          }

          if (isEnabled(groupAnalysis.cfgFlagsQaPhi)) {
            fillQaPhiByParticleSpeciesAll<DataMode::McTrack, ParticleSpeciesAll::All>();
            fillQaPhiByParticleSpeciesAll<DataMode::McTrack, ParticleSpeciesAll::Pion>();
            fillQaPhiByParticleSpeciesAll<DataMode::McTrack, ParticleSpeciesAll::Kaon>();
            fillQaPhiByParticleSpeciesAll<DataMode::McTrack, ParticleSpeciesAll::Proton>();
          }

          if (isEnabled(groupAnalysis.cfgFlagsQaPid)) {
            fillQaPidByParticleSpeciesAll<DataMode::McTrack, ParticleSpeciesAll::All>(track);
            fillQaPidByParticleSpeciesAll<DataMode::McTrack, ParticleSpeciesAll::Pion>(track);
            fillQaPidByParticleSpeciesAll<DataMode::McTrack, ParticleSpeciesAll::Kaon>(track);
            fillQaPidByParticleSpeciesAll<DataMode::McTrack, ParticleSpeciesAll::Proton>(track);
          }

          if (groupAnalysis.cfgFlagQaMc.value && (!groupTrack.cfgFlagMcParticlePhysicalPrimary.value || mcParticle.isPhysicalPrimary())) {
            hrQaMc.fill(C_CS("QaMc/hCentralityPtMcEtaMcDeltaPt"), holderEvent.centrality, holderMcParticle.pt, holderMcParticle.eta, holderTrack.pt - holderMcParticle.pt);
            hrQaMc.fill(C_CS("QaMc/hCentralityPtMcEtaMcDeltaEta"), holderEvent.centrality, holderMcParticle.pt, holderMcParticle.eta, holderTrack.eta - holderMcParticle.eta);
          }

          if (isEnabled(groupAnalysis.cfgFlagsCalculationYield) && (!groupTrack.cfgFlagMcParticlePhysicalPrimary.value || mcParticle.isPhysicalPrimary())) {
            fillCalculationYieldByParticleSpecies<DataMode::McTrack, ParticleSpecies::Pion>();
            fillCalculationYieldByParticleSpecies<DataMode::McTrack, ParticleSpecies::Kaon>();
            fillCalculationYieldByParticleSpecies<DataMode::McTrack, ParticleSpecies::Proton>();
          }

          if (isEnabled(groupAnalysis.cfgFlagsCalculationPurity) && (!groupTrack.cfgFlagMcParticlePhysicalPrimary.value || mcParticle.isPhysicalPrimary())) {
            fillCalculationPurityByParticleSpecies<ParticleSpecies::Pion>();
            fillCalculationPurityByParticleSpecies<ParticleSpecies::Kaon>();
            fillCalculationPurityByParticleSpecies<ParticleSpecies::Proton>();
          }

          if (isEnabled(groupAnalysis.cfgFlagsCalculationFractionPrimary)) {
            fillCalculationFractionPrimaryByParticleSpecies<ParticleSpecies::Pion>(mcParticle);
            fillCalculationFractionPrimaryByParticleSpecies<ParticleSpecies::Kaon>(mcParticle);
            fillCalculationFractionPrimaryByParticleSpecies<ParticleSpecies::Proton>(mcParticle);
          }

          if (isEnabled(groupAnalysis.cfgFlagsCalculationFluctuation) && (!groupTrack.cfgFlagMcParticlePhysicalPrimary.value || mcParticle.isPhysicalPrimary())) {
            calculateFluctuationByParticleNumber<DataMode::McTrack, ParticleNumber::Charge>();
            calculateFluctuationByParticleNumber<DataMode::McTrack, ParticleNumber::Kaon>();
            calculateFluctuationByParticleNumber<DataMode::McTrack, ParticleNumber::Proton>();
          }
        }

        if (isEnabled(groupAnalysis.cfgFlagsCalculationFluctuation)) {
          fillCalculationFluctuationByParticleNumber<DataMode::McTrack, ParticleNumber::Charge>();
          fillCalculationFluctuationByParticleNumber<DataMode::McTrack, ParticleNumber::Kaon>();
          fillCalculationFluctuationByParticleNumber<DataMode::McTrack, ParticleNumber::Proton>();
          miniMcCollision(holderDerivedData.nMcParticles[toI(ChargeSpecies::Plus)], holderDerivedData.nMcParticles[toI(ChargeSpecies::Minus)]);
          miniCollision(HolderDerivedData::convertFloor<aod::mini_collision::Vz::type>(holderEvent.vz * 10.), HolderDerivedData::convertFloor<aod::mini_collision::Centrality::type>(holderEvent.centrality * 500.), holderDerivedData.nTracks[toI(ChargeSpecies::Plus)], holderDerivedData.nTracks[toI(ChargeSpecies::Minus)]);
          for (auto const& signedEfficiency : holderDerivedData.signedEfficienciesMcParticle) {
            miniMcParticle(miniMcCollision.lastIndex(), signedEfficiency);
          }
          for (auto const& signedEfficiency : holderDerivedData.signedEfficienciesTrack) {
            miniTrack(miniCollision.lastIndex(), signedEfficiency);
          }
        }
      }

      break;
    }
  }
  PROCESS_SWITCH(PartNumFluc, processMc, "Process MC data", false);
};

WorkflowSpec defineDataProcessing(const ConfigContext& configContext)
{
  return WorkflowSpec{adaptAnalysisTask<PartNumFluc>(configContext)};
}
