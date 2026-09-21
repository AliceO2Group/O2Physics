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

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <ranges>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#define C_CS(cs) /* NOLINT(cppcoreguidelines-macro-usage) */                                                  \
  []<std::size_t... indices>(const std::index_sequence<indices...>) constexpr -> ConstStr<(cs)[indices]...> { \
    static_assert(std::is_array_v<std::remove_cvref_t<decltype(cs)>> &&                                       \
                  std::same_as<std::remove_extent_t<std::remove_cvref_t<decltype(cs)>>, char>);               \
    return ConstStr<(cs)[indices]...>{};                                                                      \
  }(std::make_index_sequence<sizeof(cs) - 1>{})
#define C_SV(sv) /* NOLINT(cppcoreguidelines-macro-usage) */                                                  \
  []<std::size_t... indices>(const std::index_sequence<indices...>) constexpr -> ConstStr<(sv)[indices]...> { \
    static_assert(std::same_as<std::remove_cvref_t<decltype(sv)>, std::string_view>);                         \
    return ConstStr<(sv)[indices]...>{};                                                                      \
  }(std::make_index_sequence<(sv).size()>{})

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
using JoinedMcCollisions = soa::Join<McCollisions, McCollsExtra, MultMCExtras>;
using JoinedCollisions = soa::Join<Collisions, EvSels, PVMults, FT0MultZeqs, CentNTPVs, CentFT0As, CentFT0Cs, CentFT0Ms>;
using JoinedCollisionsWithMc = soa::Join<McCollisionLabels, JoinedCollisions>;
using JoinedTracks = soa::Join<Tracks, TracksExtra, TracksDCA, TrackSelection, pidTPCFullPi, pidTPCFullKa, pidTPCFullPr, pidTOFbeta, pidTOFFullPi, pidTOFFullKa, pidTOFFullPr>;
using JoinedTracksWithMc = soa::Join<McTrackLabels, JoinedTracks>;

namespace mini_collision
{
DECLARE_SOA_COLUMN(Code, code, std::uint16_t);
} // namespace mini_collision

DECLARE_SOA_TABLE(MiniCollisions, "AOD", "MINICOLLISION", soa::Index<>, mini_collision::Code);
using MiniCollision = MiniCollisions::iterator;

namespace mini_mc_particle
{
DECLARE_SOA_INDEX_COLUMN(MiniCollision, miniCollision);
DECLARE_SOA_COLUMN(Code, code, std::uint16_t);
} // namespace mini_mc_particle

DECLARE_SOA_TABLE(MiniMcParticles, "AOD", "MINIMCPARTICLE", soa::Index<>, mini_mc_particle::MiniCollisionId, mini_mc_particle::Code);
using MiniMcParticle = MiniMcParticles::iterator;

namespace mini_track
{
DECLARE_SOA_INDEX_COLUMN(MiniCollision, miniCollision);
DECLARE_SOA_COLUMN(Code, code, std::uint32_t);
} // namespace mini_track

DECLARE_SOA_TABLE(MiniTracks, "AOD", "MINITRACK", soa::Index<>, mini_track::MiniCollisionId, mini_track::Code);
using MiniTrack = MiniTracks::iterator;

namespace tiny_mc_collision
{
DECLARE_SOA_COLUMN(NMcParticlesP, nMcParticlesP, std::uint16_t);
DECLARE_SOA_COLUMN(NMcParticlesM, nMcParticlesM, NMcParticlesP::type);
} // namespace tiny_mc_collision

DECLARE_SOA_TABLE(TinyMcCollisions, "AOD", "TINYMCCOLLISION", soa::Index<>, tiny_mc_collision::NMcParticlesP, tiny_mc_collision::NMcParticlesM);
using TinyMcCollision = TinyMcCollisions::iterator;

namespace tiny_collision
{
DECLARE_SOA_COLUMN(Code, code, mini_collision::Code::type);
DECLARE_SOA_COLUMN(NTracksP, nTracksP, std::uint16_t);
DECLARE_SOA_COLUMN(NTracksM, nTracksM, NTracksP::type);
} // namespace tiny_collision

DECLARE_SOA_TABLE(TinyCollisions, "AOD", "TINYCOLLISION", soa::Index<>, tiny_collision::Code, tiny_collision::NTracksP, tiny_collision::NTracksM);
using TinyCollision = TinyCollisions::iterator;

namespace tiny_mc_particle
{
DECLARE_SOA_INDEX_COLUMN(TinyMcCollision, tinyMcCollision);
DECLARE_SOA_COLUMN(SignedEfficiency, signedEfficiency, std::int16_t);
} // namespace tiny_mc_particle

DECLARE_SOA_TABLE(TinyMcParticles, "AOD", "TINYMCPARTICLE", soa::Index<>, tiny_mc_particle::TinyMcCollisionId, tiny_mc_particle::SignedEfficiency);
using TinyMcParticle = TinyMcParticles::iterator;

namespace tiny_track
{
DECLARE_SOA_INDEX_COLUMN(TinyCollision, tinyCollision);
DECLARE_SOA_COLUMN(SignedEfficiency, signedEfficiency, std::int16_t);
} // namespace tiny_track

DECLARE_SOA_TABLE(TinyTracks, "AOD", "TINYTRACK", soa::Index<>, tiny_track::TinyCollisionId, tiny_track::SignedEfficiency);
using TinyTrack = TinyTracks::iterator;
} // namespace o2::aod

namespace
{
namespace fluctuation_calculator_base
{
inline constexpr std::int8_t MaxOrder{8};
inline constexpr std::int32_t NExponentKeys{MaxOrder * (MaxOrder + 1) / 2};
inline constexpr std::array<std::pair<std::int8_t, std::int8_t>, NExponentKeys> ExponentKeys{[]() consteval noexcept -> std::array<std::pair<std::int8_t, std::int8_t>, NExponentKeys> {
  std::array<std::pair<std::int8_t, std::int8_t>, NExponentKeys> result{};
  std::int32_t index{};
  for (std::int32_t const& iExponent : std::views::iota(1, MaxOrder + 1)) {
    for (std::int32_t const& jExponent : std::views::iota(1, iExponent + 1)) {
      result[index++] = {static_cast<std::int8_t>(iExponent), static_cast<std::int8_t>(jExponent)};
    }
  }
  return result;
}()};
inline constexpr std::int32_t NOrderKeys{[]() consteval noexcept -> std::int32_t {
  std::array<std::int32_t, MaxOrder + 1> counts{1};
  for (std::pair<std::int8_t, std::int8_t> const& exponentKey /* o2-linter: disable=const-ref-in-for-loop */ : ExponentKeys) {
    const std::int32_t weight{exponentKey.first};
    for (std::int32_t const& sum : std::views::iota(weight, MaxOrder + 1)) {
      counts[sum] += counts[sum - weight];
    }
  }
  return std::accumulate(counts.begin(), counts.end(), 0);
}()};
inline constexpr std::array<std::array<std::int8_t, NExponentKeys>, NOrderKeys> OrderKeys{[]() consteval noexcept -> std::array<std::array<std::int8_t, NExponentKeys>, NOrderKeys> {
  constexpr auto FillOrderKeys{[](const auto& self, std::array<std::int8_t, NExponentKeys>& current, std::int32_t& index, std::array<std::array<std::int8_t, NExponentKeys>, NOrderKeys>& result, const std::int32_t position, const std::int32_t sum, const std::int32_t target) consteval noexcept -> void {
    if (position == NExponentKeys) {
      if (sum == target) {
        result[index++] = current;
      }
      return;
    }

    const std::int32_t weight{ExponentKeys[position].first};
    for (std::int32_t const& power : std::views::iota(0, (target - sum) / weight + 1)) {
      current[position] = static_cast<std::int8_t>(power);
      self(self, current, index, result, position + 1, sum + power * weight, target);
    }
  }};

  std::array<std::array<std::int8_t, NExponentKeys>, NOrderKeys> result{};
  std::array<std::int8_t, NExponentKeys> current{};
  std::int32_t index{};
  for (std::int32_t const& target : std::views::iota(0, MaxOrder + 1)) {
    FillOrderKeys(FillOrderKeys, current, index, result, 0, 0, target);
  }
  return result;
}()};
inline constexpr std::array<std::pair<std::int32_t, std::pair<std::int32_t, std::int8_t>>, NOrderKeys> OrderKeyProductLinks{[]() consteval noexcept -> std::array<std::pair<std::int32_t, std::pair<std::int32_t, std::int8_t>>, NOrderKeys> {
  std::array<std::pair<std::int32_t, std::pair<std::int32_t, std::int8_t>>, NOrderKeys> result{};
  for (std::int32_t const& iOrderKey : std::views::iota(1, NOrderKeys)) {
    result[iOrderKey].first = iOrderKey;
    std::array<std::int8_t, NExponentKeys> orderKeyParent{OrderKeys[iOrderKey]};
    for (std::int32_t const& iExponentKey : std::views::iota(0, NExponentKeys) | std::views::reverse) {
      const std::int8_t power{orderKeyParent[iExponentKey]};
      if (power <= 0) {
        continue;
      }

      orderKeyParent[iExponentKey] = {};
      for (std::int32_t const& iOrderKeyParent : std::views::iota(0, iOrderKey)) {
        if (OrderKeys[iOrderKeyParent] == orderKeyParent) {
          result[iOrderKey] = {iOrderKeyParent, {iExponentKey, power}};
          break;
        }
      }
      break;
    }
  }
  return result;
}()};
static_assert([]() consteval noexcept -> bool {
  if (OrderKeys[0] != std::array<std::int8_t, NExponentKeys>{}) {
    return false;
  }
  for (std::int32_t const& iOrderKey : std::views::iota(1, NOrderKeys)) {
    const auto& [iOrderKeyParent, factor]{OrderKeyProductLinks[iOrderKey]};
    const auto& [iExponentKey, power]{factor};
    if (iOrderKeyParent < 0 || iOrderKeyParent >= iOrderKey || iExponentKey < 0 || iExponentKey >= NExponentKeys || power <= 0) {
      return false;
    }
    for (std::int32_t const& iExponentKeyAfter : std::views::iota(iExponentKey + 1, NExponentKeys)) {
      if (OrderKeys[iOrderKey][iExponentKeyAfter] > 0) {
        return false;
      }
    }
    std::array<std::int8_t, NExponentKeys> orderKeyReconstructed{OrderKeys[iOrderKeyParent]};
    if (orderKeyReconstructed[iExponentKey] != 0) {
      return false;
    }
    orderKeyReconstructed[iExponentKey] = power;
    if (orderKeyReconstructed != OrderKeys[iOrderKey]) {
      return false;
    }
  }
  return true;
}());
} // namespace fluctuation_calculator_base

class FluctuationCalculatorTrack
{
 public:
  FluctuationCalculatorTrack() noexcept = default;
  FluctuationCalculatorTrack(const FluctuationCalculatorTrack&) noexcept = default;
  FluctuationCalculatorTrack(FluctuationCalculatorTrack&&) noexcept = default;
  FluctuationCalculatorTrack& operator=(const FluctuationCalculatorTrack&) noexcept = default;
  FluctuationCalculatorTrack& operator=(FluctuationCalculatorTrack&&) noexcept = default;
  virtual ~FluctuationCalculatorTrack() noexcept = default;

  [[nodiscard]] std::array<double, fluctuation_calculator_base::NOrderKeys> getProducts(const double weight = 1.) const noexcept
  {
    std::array<std::array<double, fluctuation_calculator_base::MaxOrder + 1>, fluctuation_calculator_base::NExponentKeys> powersQ{};
    for (std::int32_t const& iExponentKey : std::views::iota(0, fluctuation_calculator_base::NExponentKeys)) {
      powersQ[iExponentKey][1] = mQs[iExponentKey];
      for (std::int32_t const& power : std::views::iota(2, fluctuation_calculator_base::MaxOrder / fluctuation_calculator_base::ExponentKeys[iExponentKey].first + 1)) {
        powersQ[iExponentKey][power] = powersQ[iExponentKey][power - 1] * mQs[iExponentKey];
      }
    }

    std::array<double, fluctuation_calculator_base::NOrderKeys> products{weight};
    for (std::int32_t const& iOrderKey : std::views::iota(1, fluctuation_calculator_base::NOrderKeys)) {
      const auto& [iOrderKeyParent, factor]{fluctuation_calculator_base::OrderKeyProductLinks[iOrderKey]};
      const auto& [iExponentKey, power]{factor};
      products[iOrderKey] = products[iOrderKeyParent] * powersQ[iExponentKey][power];
    }
    return products;
  }
  void addQ(const std::int32_t exponentKeyIndex, const double q) noexcept { mQs[exponentKeyIndex] += q; }
  void clear() noexcept { mQs.fill({}); }
  void fill(const double charge, const double efficiency, const double weight = 1.) noexcept
  {
    const double inverseEfficiency{1. / efficiency};
    std::array<double, fluctuation_calculator_base::MaxOrder + 1> powersCharge{weight};
    std::array<double, fluctuation_calculator_base::MaxOrder + 1> powersInverseEfficiency{1.};
    for (std::int32_t const& exponent : std::views::iota(1, fluctuation_calculator_base::MaxOrder + 1)) {
      powersCharge[exponent] = powersCharge[exponent - 1] * charge;
      powersInverseEfficiency[exponent] = powersInverseEfficiency[exponent - 1] * inverseEfficiency;
    }
    for (std::int32_t const& iExponentKey : std::views::iota(0, fluctuation_calculator_base::NExponentKeys)) {
      const auto& [exponentCharge, exponentEfficiency]{fluctuation_calculator_base::ExponentKeys[iExponentKey]};
      addQ(iExponentKey, powersCharge[exponentCharge] * powersInverseEfficiency[exponentEfficiency]);
    }
  }

 protected:
  std::array<double, fluctuation_calculator_base::NExponentKeys> mQs{};
};

template <typename... Es>
concept IsValidEnum = ((std::is_enum_v<std::remove_cvref_t<Es>> && requires { std::remove_cvref_t<Es>::N; }) && ...);
template <auto... EValues>
concept IsValidEnumValue = IsValidEnum<std::remove_cvref_t<decltype(EValues)>...> && ((EValues != std::remove_cvref_t<decltype(EValues)>::N) && ...);

template <typename E>
  requires IsValidEnum<E>
constexpr std::int32_t toI(const E e) noexcept
{
  return static_cast<std::int32_t>(e);
}
template <typename E>
  requires IsValidEnum<E>
constexpr std::int32_t NEs{toI(E::N)};

enum class NameKind {
  Default = 0,
  Lower,
  Display,
  DisplayLower,
  N
};
enum class DataMode {
  McMcParticle = 0,
  McTrack,
  RawTrack,
  N
};
enum class Selection {
  Loose = 0,
  Default,
  Tight,
  N
};
enum class RangeEdge {
  Min = 0,
  Max,
  N
};
enum class McEventSelection {
  All = 0,
  Good,
  Inel,
  Vz,
  NRecoCollisions,
  N
};
enum class EventSelection {
  All = 0,
  Good,
  Run,
  Rct,
  Inel,
  Bits,
  Centrality,
  Vz,
  Occupancy,
  NPvContributors,
  N
};
enum class CentralityDefinition {
  Ft0a = 0,
  Ft0c,
  Ft0m,
  N
};
enum class DcaMeasure {
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
enum class ChargeSpecies {
  Plus = 0,
  Minus,
  N
};
enum class ParticleNumber {
  Charge = 0,
  Kaon,
  Proton,
  N
};
enum class ChargeNumber {
  Plus = 0,
  Minus,
  Total,
  Net,
  N
};

template <typename E>
  requires IsValidEnum<E>
struct EnumInfo;
template <>
struct EnumInfo<Selection> {
  static constexpr std::array<std::string_view, NEs<Selection>> Names{"Loose", "Default", "Tight"};
};
template <>
struct EnumInfo<RangeEdge> {
  static constexpr std::array<std::string_view, NEs<RangeEdge>> Names{"Min", "Max"};
};
template <>
struct EnumInfo<DcaMeasure> {
  static constexpr std::array<std::string_view, NEs<DcaMeasure>> Names{"Mean", "Sigma"};
};
template <>
struct EnumInfo<DcaAxis> {
  static constexpr std::array<std::string_view, NEs<DcaAxis>> Names{"Xy", "Z"};
};
template <>
struct EnumInfo<Detector> {
  static constexpr std::array<std::string_view, NEs<Detector>> Names{"Tpc", "Tof"};
  static constexpr std::array<std::string_view, NEs<Detector>> NamesLower{"tpc", "tof"};
  static constexpr std::array<std::string_view, NEs<Detector>> DisplayNames{"TPC", "TOF"};
};
template <>
struct EnumInfo<PidStrategy> {
  static constexpr std::array<std::string_view, NEs<PidStrategy>> Names{"Tpc", "TpcTof"};
  static constexpr std::array<std::string_view, NEs<PidStrategy>> NamesLower{"tpc", "tpcTof"};
  static constexpr std::array<std::string_view, NEs<PidStrategy>> DisplayNames{"TPC", "TPC+TOF"};
  static constexpr std::array<PidStrategyAll, NEs<PidStrategy>> PidStrategyAllValues{PidStrategyAll::Tpc, PidStrategyAll::TpcTofCombined};
};
template <>
struct EnumInfo<PidStrategyAll> {
  static constexpr std::array<Detector, NEs<PidStrategyAll>> DetectorValues{Detector::Tpc, Detector::Tof, Detector::N, Detector::N};
};
template <>
struct EnumInfo<ParticleSpecies> {
  static constexpr std::array<std::string_view, NEs<ParticleSpecies>> Names{"Pi", "Ka", "Pr"};
  static constexpr std::array<std::string_view, NEs<ParticleSpecies>> DisplayNames{"Pion", "Kaon", "Proton"};
  static constexpr std::array<std::string_view, NEs<ParticleSpecies>> DisplayNamesLower{"pion", "kaon", "proton"};
  static constexpr std::array<ParticleSpeciesAll, NEs<ParticleSpecies>> ParticleSpeciesAllValues{ParticleSpeciesAll::Pion, ParticleSpeciesAll::Kaon, ParticleSpeciesAll::Proton};
  static constexpr std::array<std::array<std::int32_t, NEs<ChargeSpecies>>, NEs<ParticleSpecies>> PdgCodes{{{PDG_t::kPiPlus, PDG_t::kPiMinus}, {PDG_t::kKPlus, PDG_t::kKMinus}, {PDG_t::kProton, PDG_t::kProtonBar}}};
  static constexpr std::array<double, NEs<ParticleSpecies>> Masses{constants::physics::MassPiPlus, constants::physics::MassKPlus, constants::physics::MassProton};
};
template <>
struct EnumInfo<ParticleSpeciesAll> {
  static constexpr std::array<std::string_view, NEs<ParticleSpeciesAll>> Names{"", "Pi", "Ka", "Pr"};
  static constexpr std::array<std::string_view, NEs<ParticleSpeciesAll>> DisplayNames{"All", "Pion", "Kaon", "Proton"};
  static constexpr std::array<std::string_view, NEs<ParticleSpeciesAll>> DisplayNamesLower{"all", "pion", "kaon", "proton"};
  static constexpr std::array<std::string_view, NEs<ParticleSpeciesAll>> Titles{"", "#pi", "K", "p"};
  static constexpr std::array<ParticleSpecies, NEs<ParticleSpeciesAll>> ParticleSpeciesValues{ParticleSpecies::N, ParticleSpecies::Pion, ParticleSpecies::Kaon, ParticleSpecies::Proton};
};
template <>
struct EnumInfo<ChargeSpecies> {
  static constexpr std::array<std::string_view, NEs<ChargeSpecies>> Names{"P", "M"};
  static constexpr std::array<std::string_view, NEs<ChargeSpecies>> NamesLower{"p", "m"};
  static constexpr std::array<std::string_view, NEs<ChargeSpecies>> Titles{"+", "#minus"};
};
template <>
struct EnumInfo<ParticleNumber> {
  static constexpr std::array<std::string_view, NEs<ParticleNumber>> Names{"Ch", "Ka", "Pr"};
  static constexpr std::array<std::string_view, NEs<ParticleNumber>> DisplayNames{"Charge", "Kaon", "Proton"};
  static constexpr std::array<std::string_view, NEs<ParticleNumber>> DisplayNamesLower{"charge", "kaon", "proton"};
  static constexpr std::array<std::string_view, NEs<ParticleNumber>> Titles{"h", "K", "p"};
  static constexpr std::array<ParticleSpeciesAll, NEs<ParticleNumber>> ParticleSpeciesAllValues{ParticleSpeciesAll::All, ParticleSpeciesAll::Kaon, ParticleSpeciesAll::Proton};
};
template <>
struct EnumInfo<ChargeNumber> {
  static constexpr std::array<std::string_view, NEs<ChargeNumber>> Names{"P", "M", "T", "N"};
};

template <typename E, NameKind NameKindValue = NameKind::Default>
  requires IsValidEnum<E>
constexpr std::string_view getName(const std::int32_t index)
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
constexpr std::string_view getName(const E e)
{
  return getName<E, NameKindValue>(toI(e));
}
template <typename E>
  requires IsValidEnum<E>
std::vector<std::string> getDisplayNames()
{
  if constexpr (requires { EnumInfo<std::remove_cvref_t<E>>::DisplayNames; }) {
    return {EnumInfo<std::remove_cvref_t<E>>::DisplayNames.begin(), EnumInfo<std::remove_cvref_t<E>>::DisplayNames.end()};
  } else {
    return {EnumInfo<std::remove_cvref_t<E>>::Names.begin(), EnumInfo<std::remove_cvref_t<E>>::Names.end()};
  }
}
template <typename E>
  requires IsValidEnum<E>
constexpr std::string_view getTitle(const std::int32_t index)
{
  return EnumInfo<std::remove_cvref_t<E>>::Titles.at(index);
}
template <typename To, typename E>
  requires IsValidEnum<To, E>
constexpr To getValue(const E e) noexcept
{
  if constexpr (std::is_same_v<std::remove_cvref_t<To>, Detector>) {
    return EnumInfo<E>::DetectorValues[toI(e)];
  } else if constexpr (std::is_same_v<std::remove_cvref_t<To>, PidStrategyAll>) {
    return EnumInfo<E>::PidStrategyAllValues[toI(e)];
  } else if constexpr (std::is_same_v<std::remove_cvref_t<To>, ParticleSpecies>) {
    return EnumInfo<E>::ParticleSpeciesValues[toI(e)];
  } else if constexpr (std::is_same_v<std::remove_cvref_t<To>, ParticleSpeciesAll>) {
    return EnumInfo<E>::ParticleSpeciesAllValues[toI(e)];
  } else {
    return To::N;
  }
}
constexpr std::int32_t getPdgCode(const ParticleSpecies particleSpecies, const ChargeSpecies chargeSpecies) noexcept
{
  return EnumInfo<ParticleSpecies>::PdgCodes[toI(particleSpecies)][toI(chargeSpecies)];
}
constexpr double getMass(const ParticleSpecies particleSpecies) noexcept
{
  return EnumInfo<ParticleSpecies>::Masses[toI(particleSpecies)];
}

namespace mini_collision_codec
{
inline constexpr std::int32_t NBinsVz{40};
inline constexpr std::array<double, NEs<RangeEdge>> RangeVz{-10., 10.};
inline constexpr std::int32_t NBinsCentrality{550};
inline constexpr std::array<double, NBinsCentrality + 1> BinEdgesCentrality{[]() consteval noexcept -> std::array<double, NBinsCentrality + 1> {
  constexpr std::array<double, 7> EdgesSelection{0., 0.001, 0.01, 0.1, 1., 10., 100.};
  constexpr std::array<std::int32_t, EdgesSelection.size() - 1> NBinsPerSection{100, 90, 90, 90, 90, 90};
  std::array<double, NBinsCentrality + 1> edges{};
  std::int32_t index{};
  edges[index++] = EdgesSelection.front();
  for (std::int32_t const& iSection : std::views::iota(0, static_cast<std::int32_t>(NBinsPerSection.size()))) {
    const double width{(EdgesSelection[iSection + 1] - EdgesSelection[iSection]) / NBinsPerSection[iSection]};
    for (std::int32_t const& iBin : std::views::iota(1, NBinsPerSection[iSection] + 1)) {
      edges[index++] = EdgesSelection[iSection] + iBin * width;
    }
  }
  return edges;
}()};
inline constexpr std::int32_t NStatesGoodNPvContributors{2};
static_assert(static_cast<std::uint64_t>(NBinsVz) * NBinsCentrality * NStatesGoodNPvContributors <= static_cast<std::uint64_t>(std::numeric_limits<aod::mini_collision::Code::type>::max()) + 1);

std::optional<aod::mini_collision::Code::type> encode(const double vz, const double centrality, const bool isGoodNPvContributors) noexcept
{
  if (!std::isfinite(vz) || !std::isfinite(centrality) || vz < RangeVz[toI(RangeEdge::Min)] || RangeVz[toI(RangeEdge::Max)] <= vz || centrality < BinEdgesCentrality.front() || BinEdgesCentrality.back() <= centrality) {
    return std::nullopt;
  }

  const std::int32_t vzBinIndex{static_cast<std::int32_t>(std::floor((vz - RangeVz[toI(RangeEdge::Min)]) / ((RangeVz[toI(RangeEdge::Max)] - RangeVz[toI(RangeEdge::Min)]) / NBinsVz)))};
  if (vzBinIndex < 0 || NBinsVz <= vzBinIndex) {
    return std::nullopt;
  }

  aod::mini_collision::Code::type code{static_cast<aod::mini_collision::Code::type>(vzBinIndex)};
  aod::mini_collision::Code::type codeSpan{NBinsVz};
  code += (static_cast<std::int32_t>(std::ranges::upper_bound(BinEdgesCentrality, centrality) - BinEdgesCentrality.begin()) - 1) * codeSpan;
  codeSpan *= NBinsCentrality;
  code += static_cast<aod::mini_collision::Code::type>(isGoodNPvContributors) * codeSpan;
  return code;
}
} // namespace mini_collision_codec

namespace mini_track_codec
{
inline constexpr std::int32_t NStatesPvContributor{2};
inline constexpr std::int32_t NBinsPt{18};
inline constexpr std::int32_t NBinsEta{16};
inline constexpr std::array<double, NEs<RangeEdge>> RangePt{0.2, 2.};
inline constexpr std::array<double, NEs<RangeEdge>> RangeEta{-0.8, 0.8};
static_assert(static_cast<std::uint64_t>(NStatesPvContributor) * NEs<Selection> * NEs<Selection> * NEs<Selection> * NEs<Selection> * NEs<Selection> * NEs<Selection> * NEs<Selection> * NEs<Selection> * NBinsPt * NBinsEta * NEs<ChargeSpecies> * NEs<ParticleSpeciesAll> * NEs<Selection> <= static_cast<std::uint64_t>(std::numeric_limits<aod::mini_track::Code::type>::max()) + 1);

struct ConfigSelection {
  std::array<double, NEs<Selection>> minItsNCls{};
  std::array<double, NEs<Selection>> maxItsChi2NCls{};
  std::array<double, NEs<Selection>> minTpcNCls{};
  std::array<double, NEs<Selection>> maxTpcChi2NCls{};
  std::array<double, NEs<Selection>> minTpcNCrossedRows{};
  std::array<double, NEs<Selection>> maxTpcNClsSharedRatio{};
  std::array<std::array<double, NEs<Selection>>, NEs<DcaAxis>> maxAbsNSigmaDca{};
  std::array<std::array<double, NEs<Selection>>, NEs<ParticleSpecies>> maxAbsNSigmaPid{};
};

template <ParticleSpeciesAll ParticleSpeciesAllValue>
  requires IsValidEnumValue<ParticleSpeciesAllValue>
std::optional<aod::mini_track::Code::type> encode(const ConfigSelection& configSelection, const bool isPvContributor, const double itsNCls, const double itsChi2NCls, const double tpcNCls, const double tpcChi2NCls, const double tpcNCrossedRows, const double tpcNClsSharedRatio, const std::array<double, NEs<DcaAxis>>& absNSigmaDca, const double pt, const double eta, const std::int32_t sign, const double absNSigmaPid) noexcept
{
  if (!std::isfinite(pt) || !std::isfinite(eta) || pt < RangePt[toI(RangeEdge::Min)] || RangePt[toI(RangeEdge::Max)] <= pt || eta < RangeEta[toI(RangeEdge::Min)] || RangeEta[toI(RangeEdge::Max)] <= eta) {
    return std::nullopt;
  }

  const std::int32_t ptBinIndex{static_cast<std::int32_t>(std::floor((pt - RangePt[toI(RangeEdge::Min)]) / ((RangePt[toI(RangeEdge::Max)] - RangePt[toI(RangeEdge::Min)]) / NBinsPt)))};
  if (ptBinIndex < 0 || NBinsPt <= ptBinIndex) {
    return std::nullopt;
  }

  const std::int32_t etaBinIndex{static_cast<std::int32_t>(std::floor((eta - RangeEta[toI(RangeEdge::Min)]) / ((RangeEta[toI(RangeEdge::Max)] - RangeEta[toI(RangeEdge::Min)]) / NBinsEta)))};
  if (etaBinIndex < 0 || NBinsEta <= etaBinIndex) {
    return std::nullopt;
  }

  static constexpr auto GetState{
    []<RangeEdge RangeEdgeValue>
      requires IsValidEnumValue<RangeEdgeValue>
    (const double value, const std::array<double, NEs<Selection>>& cuts) noexcept -> std::int32_t {
      for (std::int32_t const& iSelection : std::views::iota(0, NEs<Selection>) | std::views::reverse) {
        if constexpr (RangeEdgeValue == RangeEdge::Min) {
          if (value > cuts[iSelection]) {
            return iSelection;
          }
        } else {
          if (value < cuts[iSelection]) {
            return iSelection;
          }
        }
      }
      return -1;
    }};

  aod::mini_track::Code::type code{static_cast<aod::mini_track::Code::type>(isPvContributor)};
  aod::mini_track::Code::type codeSpan{NStatesPvContributor};
  for (std::int32_t const& state : std::to_array<std::int32_t>(
         // cppcheck-suppress internalAstError
         {GetState.template operator()<RangeEdge::Min>(itsNCls, configSelection.minItsNCls),
          GetState.template operator()<RangeEdge::Max>(itsChi2NCls, configSelection.maxItsChi2NCls),
          GetState.template operator()<RangeEdge::Min>(tpcNCls, configSelection.minTpcNCls),
          GetState.template operator()<RangeEdge::Max>(tpcChi2NCls, configSelection.maxTpcChi2NCls),
          GetState.template operator()<RangeEdge::Min>(tpcNCrossedRows, configSelection.minTpcNCrossedRows),
          GetState.template operator()<RangeEdge::Max>(tpcNClsSharedRatio, configSelection.maxTpcNClsSharedRatio),
          GetState.template operator()<RangeEdge::Max>(absNSigmaDca[toI(DcaAxis::Xy)], configSelection.maxAbsNSigmaDca[toI(DcaAxis::Xy)]),
          GetState.template operator()<RangeEdge::Max>(absNSigmaDca[toI(DcaAxis::Z)], configSelection.maxAbsNSigmaDca[toI(DcaAxis::Z)])})) {
    if (state < 0) {
      return std::nullopt;
    }

    code += state * codeSpan;
    codeSpan *= NEs<Selection>;
  }
  code += ptBinIndex * codeSpan;
  codeSpan *= NBinsPt;
  code += etaBinIndex * codeSpan;
  codeSpan *= NBinsEta;
  code += (sign > 0 ? toI(ChargeSpecies::Plus) : toI(ChargeSpecies::Minus)) * codeSpan;
  codeSpan *= NEs<ChargeSpecies>;
  code += toI(ParticleSpeciesAllValue) * codeSpan;
  codeSpan *= NEs<ParticleSpeciesAll>;
  if constexpr (ParticleSpeciesAllValue != ParticleSpeciesAll::All) {
    // clang-format v20.1.3
    // clang-format off
    const std::int32_t statePid{GetState.template operator()<RangeEdge::Max>(absNSigmaPid, configSelection.maxAbsNSigmaPid[toI(getValue<ParticleSpecies>(ParticleSpeciesAllValue))])};
    // clang-format on
    if (statePid < 0) {
      return std::nullopt;
    }
    code += statePid * codeSpan;
  }
  return code;
}
} // namespace mini_track_codec

namespace mini_mc_particle_codec
{
static_assert(static_cast<std::uint64_t>(mini_track_codec::NBinsPt) * mini_track_codec::NBinsEta * NEs<ChargeSpecies> * NEs<ParticleSpecies> <= static_cast<std::uint64_t>(std::numeric_limits<aod::mini_mc_particle::Code::type>::max()) + 1);

template <ParticleSpecies ParticleSpeciesValue>
  requires IsValidEnumValue<ParticleSpeciesValue>
std::optional<aod::mini_mc_particle::Code::type> encode(const double pt, const double eta, const std::int32_t sign) noexcept
{
  if (!std::isfinite(pt) || !std::isfinite(eta) || pt < mini_track_codec::RangePt[toI(RangeEdge::Min)] || mini_track_codec::RangePt[toI(RangeEdge::Max)] <= pt || eta < mini_track_codec::RangeEta[toI(RangeEdge::Min)] || mini_track_codec::RangeEta[toI(RangeEdge::Max)] <= eta) {
    return std::nullopt;
  }

  const std::int32_t ptBinIndex{static_cast<std::int32_t>(std::floor((pt - mini_track_codec::RangePt[toI(RangeEdge::Min)]) / ((mini_track_codec::RangePt[toI(RangeEdge::Max)] - mini_track_codec::RangePt[toI(RangeEdge::Min)]) / mini_track_codec::NBinsPt)))};
  if (ptBinIndex < 0 || mini_track_codec::NBinsPt <= ptBinIndex) {
    return std::nullopt;
  }

  const std::int32_t etaBinIndex{static_cast<std::int32_t>(std::floor((eta - mini_track_codec::RangeEta[toI(RangeEdge::Min)]) / ((mini_track_codec::RangeEta[toI(RangeEdge::Max)] - mini_track_codec::RangeEta[toI(RangeEdge::Min)]) / mini_track_codec::NBinsEta)))};
  if (etaBinIndex < 0 || mini_track_codec::NBinsEta <= etaBinIndex) {
    return std::nullopt;
  }

  aod::mini_mc_particle::Code::type code{static_cast<aod::mini_mc_particle::Code::type>(ptBinIndex)};
  aod::mini_mc_particle::Code::type codeSpan{mini_track_codec::NBinsPt};
  code += etaBinIndex * codeSpan;
  codeSpan *= mini_track_codec::NBinsEta;
  code += (sign > 0 ? toI(ChargeSpecies::Plus) : toI(ChargeSpecies::Minus)) * codeSpan;
  codeSpan *= NEs<ChargeSpecies>;
  code += toI(ParticleSpeciesValue) * codeSpan;
  return code;
}
} // namespace mini_mc_particle_codec

template <typename T>
  requires std::is_arithmetic_v<T>
constexpr std::int32_t nEnabled(const Configurable<LabeledArray<T>>& cfg)
{
  const LabeledArray<T>& la{cfg.value};
  return std::count_if(la[0], la[0] + la.rows() * la.cols(), [](const T x) constexpr -> bool { return x != T{}; });
}
template <typename T>
  requires std::is_arithmetic_v<T>
constexpr bool isEnabled(const Configurable<LabeledArray<T>>& cfg)
{
  return nEnabled(cfg) > 0;
}

double interpolate(const TH3* const h, const double x, const double y, const double z)
{
  if (!h) {
    return 0.;
  }

  static constexpr auto GetBinIndicesWeights{[](const TAxis* const axis, const double position) -> std::array<std::pair<std::int32_t, double>, NEs<RangeEdge>> {
    if (!axis) {
      return {};
    }

    const std::int32_t nBins{axis->GetNbins()};
    if (nBins == 1 || position <= axis->GetBinCenter(1)) {
      return {{{1, 1.}, {1, 0.}}};
    }
    if (position >= axis->GetBinCenter(nBins)) {
      return {{{nBins, 1.}, {nBins, 0.}}};
    }

    const std::int32_t binIndex{axis->FindFixBin(position)};
    const std::int32_t binIndexLower{position < axis->GetBinCenter(binIndex) ? binIndex - 1 : binIndex};
    const std::int32_t binIndexUpper{binIndexLower + 1};
    const double fraction{(position - axis->GetBinCenter(binIndexLower)) / (axis->GetBinCenter(binIndexUpper) - axis->GetBinCenter(binIndexLower))};

    return {{{binIndexLower, 1. - fraction}, {binIndexUpper, fraction}}};
  }};

  const std::array<std::pair<std::int32_t, double>, NEs<RangeEdge>> binIndicesWeightsX{GetBinIndicesWeights(h->GetXaxis(), x)};
  const std::array<std::pair<std::int32_t, double>, NEs<RangeEdge>> binIndicesWeightsY{GetBinIndicesWeights(h->GetYaxis(), y)};
  const std::array<std::pair<std::int32_t, double>, NEs<RangeEdge>> binIndicesWeightsZ{GetBinIndicesWeights(h->GetZaxis(), z)};

  double result{};
  for (const auto& [iBinX, weightX] : binIndicesWeightsX) {
    for (const auto& [iBinY, weightY] : binIndicesWeightsY) {
      for (const auto& [iBinZ, weightZ] : binIndicesWeightsZ) {
        result += weightX * weightY * weightZ * h->GetBinContent(iBinX, iBinY, iBinZ);
      }
    }
  }
  return result;
}
} // namespace

struct PartNumFluc {
  struct HolderCcdb {
    static constexpr std::int32_t NDimensionsEfficiency{4};

    const TList* lCcdb{};
    std::map<std::int32_t, std::pair<std::int32_t, std::int32_t>> runNumbersIndicesGroupIndices;
    std::int32_t runGroupIndexCurrent{};
    std::array<std::array<std::array<std::pair<const TFormula*, const TH3*>, NEs<ChargeSpecies>>, NEs<DcaAxis>>, NEs<DcaMeasure>> calibrationsPtMeasureDca{};
    std::array<std::array<std::array<const TH3*, NEs<ChargeSpecies>>, NEs<ParticleSpecies>>, NEs<Detector>> hsCentralityPtEtaShiftNSigmaPid{};
    std::array<std::array<std::array<const THnBase*, NEs<ChargeSpecies>>, NEs<ParticleSpecies>>, NEs<PidStrategy>> hsVzCentralityPtEtaEfficiency{};
  } holderCcdb{};

  struct HolderMcEvent {
    double vz{};
    std::array<std::array<std::int32_t, NEs<ChargeSpecies>>, NEs<ParticleNumber>> numbers{};
    std::array<std::array<std::int32_t, NEs<ChargeSpecies>>, NEs<ParticleNumber>> numbersEff{};

    void clear() noexcept { *this = {}; }
  } holderMcEvent{};

  struct HolderEvent {
    static constexpr std::array<double, NEs<RangeEdge>> RangeCentrality{0., 100.};
    static constexpr bool isValidCentrality(const double value) noexcept { return RangeCentrality[toI(RangeEdge::Min)] <= value && value < RangeCentrality[toI(RangeEdge::Max)]; }

    std::int32_t runNumber{};
    std::int32_t runIndex{};
    std::int32_t runGroupIndex{};
    double vz{};
    std::array<std::int32_t, NEs<ChargeSpecies>> nGlobalTracks{};
    std::array<std::int32_t, NEs<ChargeSpecies>> nPvContributors{};
    std::array<std::array<std::array<double, NEs<ChargeSpecies>>, NEs<DcaAxis>>, NEs<DcaMeasure>> measureDca{};
    std::array<std::int32_t, NEs<ChargeSpecies>> nTofBeta{};
    double centralityCalibration{};
    double centrality{};
    std::int32_t subgroupIndex{};
    std::array<std::array<std::int32_t, NEs<ChargeSpecies>>, NEs<ParticleNumber>> numbers{};

    void clear() noexcept { *this = {}; }
    [[nodiscard]] std::int32_t getNGlobalTracks() const noexcept { return std::accumulate(nGlobalTracks.begin(), nGlobalTracks.end(), 0); }
    [[nodiscard]] std::int32_t getNPvContributors() const noexcept { return std::accumulate(nPvContributors.begin(), nPvContributors.end(), 0); }
    template <DcaMeasure DcaMeasureValue, DcaAxis DcaAxisValue>
      requires IsValidEnumValue<DcaMeasureValue, DcaAxisValue>
    [[nodiscard]] double getMeasureDca() const noexcept
    {
      const std::int32_t sumNGlobalTracks{getNGlobalTracks()};
      if (sumNGlobalTracks == 0) {
        return 0.;
      }

      double sumDca{};
      if constexpr (DcaMeasureValue == DcaMeasure::Sigma) {
        const double meanDca{getMeasureDca<DcaMeasure::Mean, DcaAxisValue>()};
        for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
          sumDca += (std::pow(measureDca[toI(DcaMeasure::Sigma)][toI(DcaAxisValue)][iChargeSpecies], 2.) + std::pow(measureDca[toI(DcaMeasure::Mean)][toI(DcaAxisValue)][iChargeSpecies] - meanDca, 2.)) * nGlobalTracks[iChargeSpecies];
        }
        return std::sqrt(sumDca / sumNGlobalTracks);
      } else {
        for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
          sumDca += measureDca[toI(DcaMeasure::Mean)][toI(DcaAxisValue)][iChargeSpecies] * nGlobalTracks[iChargeSpecies];
        }
        return sumDca / sumNGlobalTracks;
      }
    }
    [[nodiscard]] std::int32_t getNTofBeta() const noexcept { return std::accumulate(nTofBeta.begin(), nTofBeta.end(), 0); }
  } holderEvent{};

  struct HolderMcParticle {
    std::int32_t pdgCode{};
    std::int32_t charge{};
    double pt{};
    double eta{};

    void clear() noexcept { *this = {}; }
  } holderMcParticle{};

  struct HolderTrack {
    static constexpr double TruncationAbsNSigmaPid{999.};
    static constexpr double truncateNSigmaPid(const double value, const double shift = {}) noexcept
    {
      const double valueShifted{value - shift};
      return std::abs(value) < TruncationAbsNSigmaPid && std::abs(valueShifted) < TruncationAbsNSigmaPid ? valueShifted : -TruncationAbsNSigmaPid;
    }

    [[nodiscard]] double getNSigmaPidCombined(const std::int32_t particleSpeciesIndex) const noexcept
    {
      return truncateNSigmaPid(std::copysign(std::hypot(nSigmaPid[toI(Detector::Tpc)][particleSpeciesIndex], nSigmaPid[toI(Detector::Tof)][particleSpeciesIndex]), nSigmaPid[toI(Detector::Tpc)][particleSpeciesIndex] + nSigmaPid[toI(Detector::Tof)][particleSpeciesIndex]));
    }

    std::array<double, NEs<DcaAxis>> dca{};
    std::int32_t sign{};
    double pt{};
    double eta{};
    double phi{};
    std::array<bool, NEs<Detector>> hasPid{};
    std::array<std::array<double, NEs<ParticleSpecies>>, NEs<Detector>> nSigmaPid{[]() constexpr -> std::array<std::array<double, NEs<ParticleSpecies>>, NEs<Detector>> {
      std::array<std::array<double, NEs<ParticleSpecies>>, NEs<Detector>> a{};
      std::array<double, NEs<ParticleSpecies>> b{};
      b.fill(-TruncationAbsNSigmaPid);
      a.fill(b);
      return a;
    }()};

    void clear() noexcept { *this = {}; }
  } holderTrack{};

  struct HolderDerivedData {
    template <std::integral T>
    static constexpr T convert(const double value) noexcept
    {
      return std::numeric_limits<T>::lowest() <= value && value <= std::numeric_limits<T>::max() ? static_cast<T>(std::rint(value)) : (std::is_signed_v<T> ? std::numeric_limits<T>::lowest() : std::numeric_limits<T>::max());
    }

    std::array<aod::tiny_mc_collision::NMcParticlesP::type, NEs<ChargeSpecies>> nMcParticles{};
    std::array<aod::tiny_collision::NTracksP::type, NEs<ChargeSpecies>> nTracks{};
    std::vector<aod::tiny_mc_particle::SignedEfficiency::type> signedEfficienciesMcParticle{[]() constexpr -> std::vector<aod::tiny_mc_particle::SignedEfficiency::type> { std::vector<aod::tiny_mc_particle::SignedEfficiency::type> v{}; v.reserve(256); return v; }()};
    std::vector<aod::tiny_track::SignedEfficiency::type> signedEfficienciesTrack{[]() constexpr -> std::vector<aod::tiny_track::SignedEfficiency::type> { std::vector<aod::tiny_track::SignedEfficiency::type> v{}; v.reserve(256); return v; }()};

    void clear() noexcept
    {
      nMcParticles = {};
      nTracks = {};
      signedEfficienciesMcParticle.clear();
      signedEfficienciesTrack.clear();
    }
  } holderDerivedData{};

  struct : ConfigurableGroup {
    std::string prefix{"cgCcdb"};
    Configurable<std::string> cfgUrl{"cfgUrl", "http://ccdb-test.cern.ch:8080", "Url of CCDB"};
    Configurable<std::string> cfgPath{"cfgPath", "Users/f/fasi/test", "Path in CCDB"};
    Configurable<std::int64_t> cfgTimestampLatest{"cfgTimestampLatest", -1, "Latest timestamp in CCDB"};
  } cgCcdb{};

  struct : ConfigurableGroup {
    std::string prefix{"cgAnalysis"};
    Configurable<bool> cfgFlagQaRun{"cfgFlagQaRun", false, "Run QA flag"};
    Configurable<bool> cfgFlagQaEvent{"cfgFlagQaEvent", false, "Event QA flag"};
    Configurable<bool> cfgFlagQaCentrality{"cfgFlagQaCentrality", false, "Centrality QA flag"};
    Configurable<bool> cfgFlagQaTrack{"cfgFlagQaTrack", false, "Track QA flag"};
    Configurable<bool> cfgFlagQaDca{"cfgFlagQaDca", false, "DCA QA flag"};
    Configurable<LabeledArray<std::int32_t>> cfgFlagsQaAcceptance{"cfgFlagsQaAcceptance", {std::array<std::int32_t, NEs<ParticleSpeciesAll>>{}.data(), NEs<ParticleSpeciesAll>, getDisplayNames<ParticleSpeciesAll>()}, "Acceptance QA flags"};
    Configurable<LabeledArray<std::int32_t>> cfgFlagsQaPhi{"cfgFlagsQaPhi", {std::array<std::int32_t, NEs<ParticleSpeciesAll>>{}.data(), NEs<ParticleSpeciesAll>, getDisplayNames<ParticleSpeciesAll>()}, "Phi QA flags"};
    Configurable<LabeledArray<std::int32_t>> cfgFlagsQaPid{"cfgFlagsQaPid", {std::array<std::int32_t, NEs<ParticleSpeciesAll>>{}.data(), NEs<ParticleSpeciesAll>, getDisplayNames<ParticleSpeciesAll>()}, "PID QA flags"};
    Configurable<bool> cfgFlagQaMc{"cfgFlagQaMc", false, "MC QA flag"};
    Configurable<LabeledArray<std::int32_t>> cfgFlagsCalculationYield{"cfgFlagsCalculationYield", {std::array<std::int32_t, NEs<ParticleSpecies>>{}.data(), NEs<ParticleSpecies>, getDisplayNames<ParticleSpecies>()}, "Yield calculation flags"};
    Configurable<LabeledArray<std::int32_t>> cfgFlagsCalculationPurity{"cfgFlagsCalculationPurity", {std::array<std::int32_t, NEs<ParticleSpecies>>{}.data(), NEs<ParticleSpecies>, getDisplayNames<ParticleSpecies>()}, "Purity calculation flags"};
    Configurable<LabeledArray<std::int32_t>> cfgFlagsCalculationFractionPrimary{"cfgFlagsCalculationFractionPrimary", {std::array<std::int32_t, NEs<ParticleSpecies>>{}.data(), NEs<ParticleSpecies>, getDisplayNames<ParticleSpecies>()}, "Primary fraction calculation flags"};
    Configurable<LabeledArray<std::int32_t>> cfgFlagsCalculationFluctuation{"cfgFlagsCalculationFluctuation", {std::array<std::int32_t, NEs<ParticleNumber>>{}.data(), NEs<ParticleNumber>, getDisplayNames<ParticleNumber>()}, "Fluctuation calculation flags"};
    Configurable<LabeledArray<std::int32_t>> cfgFlagsStorageMiniTable{"cfgFlagsStorageMiniTable", {std::array<std::int32_t, NEs<ParticleNumber>>{}.data(), NEs<ParticleNumber>, getDisplayNames<ParticleNumber>()}, "Mini table storage flags"};
  } cgAnalysis{};

  struct : ConfigurableGroup {
    std::string prefix{"cgEvent"};
    Configurable<bool> cfgFlagRejectionRunBad{"cfgFlagRejectionRunBad", false, "Bad run rejection flag"};
    Configurable<bool> cfgFlagRejectionRunBadMc{"cfgFlagRejectionRunBadMc", false, "MC bad run rejection flag"};
    Configurable<std::string> cfgLabelFlagsRct{"cfgLabelFlagsRct", "CBT_hadronPID", "RCT flags label"};
    Configurable<LabeledArray<std::int32_t>> cfgFlagsRct{"cfgFlagsRct", {std::array<std::int32_t, 3>{0, 1, 1}.data(), 3, {"ZDC", "Acceptance", "Table"}}, "RCT flags"};
    Configurable<std::uint64_t> cfgBitsSelection{"cfgBitsSelection", std::uint64_t{0b00000000000001000000000000000000000000000000000000}, "Event selection bits"};
    Configurable<bool> cfgFlagInel{"cfgFlagInel", true, "Flag of requiring INEL > 0"};
    Configurable<bool> cfgFlagInelMc{"cfgFlagInelMc", true, "Flag of requiring MC INEL > 0"};
    Configurable<double> cfgCutMaxAbsVz{"cfgCutMaxAbsVz", 8., "Maximum absolute z-vertex position (cm)"};
    Configurable<double> cfgCutMaxAbsVzMc{"cfgCutMaxAbsVzMc", 8., "Maximum absolute MC z-vertex position (cm)"};
    Configurable<std::int32_t> cfgCutMaxOccupancy{"cfgCutMaxOccupancy", -1, "Maximum occupancy"};
    Configurable<std::int32_t> cfgCutMinDeviationNPvContributors{"cfgCutMinDeviationNPvContributors", -4, "Minimum nPvContributors deviation from nGlobalTracks"};
    Configurable<std::int32_t> cfgIndexDefinitionCentrality{"cfgIndexDefinitionCentrality", 2, "Centrality definition index"};
    Configurable<bool> cfgFlagDefinitionCentralitySameQa{"cfgFlagDefinitionCentralitySameQa", false, "Flag of using the same centrality definition for QA"};
    Configurable<std::int32_t> cfgNMultiplicityBinsAxis{"cfgNMultiplicityBinsAxis", 200, "Number of multiplicity bins for axis"};
    ConfigurableAxis cfgAxisCentralityCalibration{"cfgAxisCentralityCalibration", {VARIABLE_WIDTH, 0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 100.}, "Centrality axis in calibration"};
    ConfigurableAxis cfgAxisCentrality{"cfgAxisCentrality", {20, 0., 100.}, "Centrality axis in fluctuation calculation"};
    Configurable<std::int32_t> cfgNSubgroups{"cfgNSubgroups", 20, "Number of subgroups in fluctuation calculation"};
    Configurable<double> cfgFactorStorageMiniTable{"cfgFactorStorageMiniTable", 1., "Mini table storage inverse probability factor"};
    Configurable<bool> cfgFlagSingleCollisionMc{"cfgFlagSingleCollisionMc", false, "Flag of requiring exactly single collision of MC collision"};
    Configurable<bool> cfgFlagMcCollisionVz{"cfgFlagMcCollisionVz", true, "Flag of using z-vertex position of MC collision"};
  } cgEvent{};

  struct : ConfigurableGroup {
    std::string prefix{"cgTrack"};
    Configurable<bool> cfgFlagPvContributor{"cfgFlagPvContributor", true, "Flag of requiring PV contributor"};
    Configurable<LabeledArray<std::int32_t>> cfgCutsMinItsNCls{"cfgCutsMinItsNCls", {std::array<std::int32_t, NEs<Selection>>{4, 5, 6}.data(), NEs<Selection>, getDisplayNames<Selection>()}, "Minimum numbers of clusters ITS"};
    Configurable<LabeledArray<double>> cfgCutsMaxItsChi2NCls{"cfgCutsMaxItsChi2NCls", {std::array<double, NEs<Selection>>{30., 25., 20.}.data(), NEs<Selection>, getDisplayNames<Selection>()}, "Maximum chi2 per cluster ITS"};
    Configurable<LabeledArray<std::int32_t>> cfgCutsMinTpcNCls{"cfgCutsMinTpcNCls", {std::array<std::int32_t, NEs<Selection>>{45, 55, 65}.data(), NEs<Selection>, getDisplayNames<Selection>()}, "Minimum numbers of clusters TPC"};
    Configurable<double> cfgCutMinTpcChi2NCls{"cfgCutMinTpcChi2NCls", 0., "Minimum chi2 per cluster TPC"};
    Configurable<LabeledArray<double>> cfgCutsMaxTpcChi2NCls{"cfgCutsMaxTpcChi2NCls", {std::array<double, NEs<Selection>>{4., 3.5, 3.}.data(), NEs<Selection>, getDisplayNames<Selection>()}, "Maximum chi2 per cluster TPC"};
    Configurable<LabeledArray<double>> cfgCutsMaxTpcNClsSharedRatio{"cfgCutsMaxTpcNClsSharedRatio", {std::array<double, NEs<Selection>>{0.5, 0.4, 0.3}.data(), NEs<Selection>, getDisplayNames<Selection>()}, "Maximum ratios of shared clusters over clusters TPC"};
    Configurable<LabeledArray<std::int32_t>> cfgCutsMinTpcNCrossedRows{"cfgCutsMinTpcNCrossedRows", {std::array<std::int32_t, NEs<Selection>>{70, 80, 90}.data(), NEs<Selection>, getDisplayNames<Selection>()}, "Minimum numbers of crossed rows TPC"};
    Configurable<double> cfgCutMinTpcNCrossedRowsRatio{"cfgCutMinTpcNCrossedRowsRatio", 0.8, "Minimum ratio of crossed rows over findable clusters TPC"};
    Configurable<bool> cfgFlagRecalibrationDca{"cfgFlagRecalibrationDca", false, "DCA recalibration flag"};
    Configurable<LabeledArray<double>> cfgCutsMaxAbsNSigmaDca{"cfgCutsMaxAbsNSigmaDca", {std::array<double, NEs<DcaAxis> * NEs<Selection>>{3., 2.5, 2., 3., 2.5, 2.}.data(), NEs<DcaAxis>, NEs<Selection>, getDisplayNames<DcaAxis>(), getDisplayNames<Selection>()}, "Maximum absolute nSigma values of DCA (cm)"};
    Configurable<LabeledArray<double>> cfgCutsRangePt{"cfgCutsRangePt", {std::array<double, NEs<ParticleSpecies> * NEs<RangeEdge>>{0.2, 2., 0.3, 2., 0.4, 2.}.data(), NEs<ParticleSpecies>, NEs<RangeEdge>, getDisplayNames<ParticleSpecies>(), getDisplayNames<RangeEdge>()}, "pT ranges (GeV/c)"};
    Configurable<double> cfgCutMaxAbsEta{"cfgCutMaxAbsEta", 0.8, "Maximum absolute eta"};
    Configurable<LabeledArray<double>> cfgThresholdsPtTofPid{"cfgThresholdsPtTofPid", {std::array<double, NEs<ParticleSpecies>>{0.4, 0.4, 0.8}.data(), NEs<ParticleSpecies>, getDisplayNames<ParticleSpecies>()}, "pT (GeV/c) thresholds for TOF PID"};
    Configurable<LabeledArray<std::int32_t>> cfgFlagsRecalibrationNSigmaPid{"cfgFlagsRecalibrationNSigmaPid", {std::array<std::int32_t, NEs<ParticleSpecies>>{}.data(), NEs<ParticleSpecies>, getDisplayNames<ParticleSpecies>()}, "nSigma PID recalibration flags"};
    Configurable<bool> cfgFlagRejectionOthers{"cfgFlagRejectionOthers", false, "Other particle species rejection flag"};
    Configurable<LabeledArray<double>> cfgCutsMaxAbsNSigmaPid{"cfgCutsMaxAbsNSigmaPid", {std::array<double, NEs<ParticleSpecies> * NEs<Selection>>{2.5, 2., 1.5, 2.5, 2., 1.5, 2.5, 2., 1.5}.data(), NEs<ParticleSpecies>, NEs<Selection>, getDisplayNames<ParticleSpecies>(), getDisplayNames<Selection>()}, "Maximum absolute nSigma values for PID"};
    Configurable<bool> cfgFlagMcParticlePhysicalPrimary{"cfgFlagMcParticlePhysicalPrimary", true, "Flag of requiring physical primary MC particle"};
    Configurable<bool> cfgFlagMcParticleMomentum{"cfgFlagMcParticleMomentum", true, "Flag of using momentum of MC particle"};
  } cgTrack{};

  Service<framework::O2DatabasePDG> pdg{};
  Service<ccdb::BasicCCDBManager> ccdb{};

  aod::rctsel::RCTFlagsChecker rctFlagsChecker;

  Filter filterCollision{aod::evsel::sel8 == true};
  Filter filterTrack{requireQualityTracksInFilter() && requireTrackCutInFilter(TrackSelectionFlags::kGoldenChi2)};
  Filter filterMcCollision{aod::mccollisionprop::numRecoCollision > 0};

  Preslice<aod::JoinedTracksWithMc> presliceTracksPerCollision{aod::track::collisionId};
  PresliceUnsorted<aod::JoinedTracksWithMc> presliceTracksPerMcParticle{aod::mctracklabel::mcParticleId};

  Produces<aod::MiniCollisions> miniCollision{};
  Produces<aod::MiniMcParticles> miniMcParticle{};
  Produces<aod::MiniTracks> miniTrack{};
  Produces<aod::TinyMcCollisions> tinyMcCollision{};
  Produces<aod::TinyCollisions> tinyCollision{};
  Produces<aod::TinyMcParticles> tinyMcParticle{};
  Produces<aod::TinyTracks> tinyTrack{};

  HistogramRegistry hrCalculationFluctuation{"CalculationFluctuation", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry hrCalculationFractionPrimary{"CalculationFractionPrimary", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry hrCalculationPurity{"CalculationPurity", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry hrCalculationYield{"CalculationYield", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry hrQaMc{"QaMc", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry hrQaPid{"QaPid", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry hrQaPhi{"QaPhi", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry hrQaAcceptance{"QaAcceptance", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry hrQaDca{"QaDca", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry hrQaTrack{"QaTrack", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry hrQaCentrality{"QaCentrality", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry hrQaEvent{"QaEvent", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry hrQaRun{"QaRun", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry hrCounter{"Counter", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  bool doQaAcceptance{};
  bool doQaPhi{};
  bool doQaPid{};
  bool doCalculationYield{};
  bool doCalculationPurity{};
  bool doCalculationFractionPrimary{};
  bool doCalculationFluctuation{};
  bool doStorageMiniTable{};
  std::array<double, NEs<RangeEdge>> rangePtAll{std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()};
  std::unique_ptr<mini_track_codec::ConfigSelection> configSelection;

  std::array<std::array<std::unique_ptr<FluctuationCalculatorTrack>, NEs<ChargeNumber>>, NEs<ParticleNumber>> fluctuationCalculatorsTrackMcParticle{};
  std::array<std::array<std::unique_ptr<FluctuationCalculatorTrack>, NEs<ChargeNumber>>, NEs<ParticleNumber>> fluctuationCalculatorsTrackTrack{};

  void init(const InitContext&)
  {
    gRandom->SetSeed(0);

    doQaAcceptance = isEnabled(cgAnalysis.cfgFlagsQaAcceptance);
    doQaPhi = isEnabled(cgAnalysis.cfgFlagsQaPhi);
    doQaPid = isEnabled(cgAnalysis.cfgFlagsQaPid);
    doCalculationYield = isEnabled(cgAnalysis.cfgFlagsCalculationYield);
    doCalculationPurity = isEnabled(cgAnalysis.cfgFlagsCalculationPurity);
    doCalculationFractionPrimary = isEnabled(cgAnalysis.cfgFlagsCalculationFractionPrimary);
    doCalculationFluctuation = isEnabled(cgAnalysis.cfgFlagsCalculationFluctuation);
    doStorageMiniTable = isEnabled(cgAnalysis.cfgFlagsStorageMiniTable);

    if (doProcessRaw.value == doProcessMc.value) {
      LOG(fatal) << "Identical values of doProcessRaw and doProcessMc!";
    }
    if (doProcessMc.value) {
      LOG(info) << "Enabling MC data process.";
    } else {
      LOG(info) << "Enabling raw data process.";
    }
    if (nEnabled(cgAnalysis.cfgFlagsCalculationFluctuation) > 1) {
      LOG(fatal) << "Invalid " << cgAnalysis.cfgFlagsCalculationFluctuation.name << "!";
    }
    if (nEnabled(cgAnalysis.cfgFlagsStorageMiniTable) > 1) {
      LOG(fatal) << "Invalid " << cgAnalysis.cfgFlagsStorageMiniTable.name << "!";
    }
    if ((cgAnalysis.cfgFlagQaEvent.value || cgAnalysis.cfgFlagQaCentrality.value || doCalculationFluctuation) && cgEvent.cfgNMultiplicityBinsAxis.value <= 0) {
      LOG(fatal) << "Invalid " << cgEvent.cfgNMultiplicityBinsAxis.name << "!";
    }
    if (doCalculationFluctuation && cgEvent.cfgNSubgroups.value <= 0) {
      LOG(fatal) << "Invalid " << cgEvent.cfgNSubgroups.name << "!";
    }

    if (doStorageMiniTable) {
      configSelection = std::make_unique<mini_track_codec::ConfigSelection>();
      for (std::int32_t const& iSelection : std::views::iota(0, NEs<Selection>)) {
        configSelection->minItsNCls[iSelection] = cgTrack.cfgCutsMinItsNCls.value.get(iSelection);
        configSelection->maxItsChi2NCls[iSelection] = cgTrack.cfgCutsMaxItsChi2NCls.value.get(iSelection);
        configSelection->minTpcNCls[iSelection] = cgTrack.cfgCutsMinTpcNCls.value.get(iSelection);
        configSelection->maxTpcChi2NCls[iSelection] = cgTrack.cfgCutsMaxTpcChi2NCls.value.get(iSelection);
        configSelection->minTpcNCrossedRows[iSelection] = cgTrack.cfgCutsMinTpcNCrossedRows.value.get(iSelection);
        configSelection->maxTpcNClsSharedRatio[iSelection] = cgTrack.cfgCutsMaxTpcNClsSharedRatio.value.get(iSelection);
      }
      for (std::int32_t const& iDcaAxis : std::views::iota(0, NEs<DcaAxis>)) {
        for (std::int32_t const& iSelection : std::views::iota(0, NEs<Selection>)) {
          configSelection->maxAbsNSigmaDca[iDcaAxis][iSelection] = cgTrack.cfgCutsMaxAbsNSigmaDca.value.get(getName<DcaAxis>(iDcaAxis).data(), getName<Selection>(iSelection).data());
        }
      }
      for (std::int32_t const& iParticleSpecies : std::views::iota(0, NEs<ParticleSpecies>)) {
        for (std::int32_t const& iSelection : std::views::iota(0, NEs<Selection>)) {
          configSelection->maxAbsNSigmaPid[iParticleSpecies][iSelection] = cgTrack.cfgCutsMaxAbsNSigmaPid.value.get(getName<ParticleSpecies, NameKind::Display>(iParticleSpecies).data(), getName<Selection>(iSelection).data());
        }
      }

      static constexpr auto ValidateCuts{
        []<RangeEdge RangeEdgeValue>
          requires IsValidEnumValue<RangeEdgeValue>
        (const auto& cuts, const auto& name) constexpr -> void {
          for (std::int32_t const& iSelection : std::views::iota(0, NEs<Selection>)) {
            if (!std::isfinite(cuts[iSelection])) {
              LOG(fatal) << "Non-finite value in " << name << "!";
            }
            if (iSelection == 0) {
              continue;
            }
            if constexpr (RangeEdgeValue == RangeEdge::Min) {
              if (!(cuts[iSelection - 1] < cuts[iSelection])) {
                LOG(fatal) << "Values in " << name << " must satisfy " << getName(Selection::Loose) << " < " << getName(Selection::Default) << " < " << getName(Selection::Tight) << "!";
              }
            } else {
              if (!(cuts[iSelection - 1] > cuts[iSelection])) {
                LOG(fatal) << "Values in " << name << " must satisfy " << getName(Selection::Loose) << " > " << getName(Selection::Default) << " > " << getName(Selection::Tight) << "!";
              }
            }
          }
        }};
      ValidateCuts.template operator()<RangeEdge::Min>(configSelection->minItsNCls, cgTrack.cfgCutsMinItsNCls.name);
      ValidateCuts.template operator()<RangeEdge::Max>(configSelection->maxItsChi2NCls, cgTrack.cfgCutsMaxItsChi2NCls.name);
      ValidateCuts.template operator()<RangeEdge::Min>(configSelection->minTpcNCls, cgTrack.cfgCutsMinTpcNCls.name);
      ValidateCuts.template operator()<RangeEdge::Max>(configSelection->maxTpcChi2NCls, cgTrack.cfgCutsMaxTpcChi2NCls.name);
      ValidateCuts.template operator()<RangeEdge::Min>(configSelection->minTpcNCrossedRows, cgTrack.cfgCutsMinTpcNCrossedRows.name);
      ValidateCuts.template operator()<RangeEdge::Max>(configSelection->maxTpcNClsSharedRatio, cgTrack.cfgCutsMaxTpcNClsSharedRatio.name);
      for (std::int32_t const& iDcaAxis : std::views::iota(0, NEs<DcaAxis>)) {
        ValidateCuts.template operator()<RangeEdge::Max>(configSelection->maxAbsNSigmaDca[iDcaAxis], cgTrack.cfgCutsMaxAbsNSigmaDca.name);
      }
      for (std::int32_t const& iParticleSpecies : std::views::iota(0, NEs<ParticleSpecies>)) {
        ValidateCuts.template operator()<RangeEdge::Max>(configSelection->maxAbsNSigmaPid[iParticleSpecies], cgTrack.cfgCutsMaxAbsNSigmaPid.name);
      }
    }

    if (doStorageMiniTable || static_cast<bool>(cgAnalysis.cfgFlagsCalculationFluctuation.value.get(toI(ParticleNumber::Charge)))) {
      for (std::int32_t const& iParticleSpecies : std::views::iota(0, NEs<ParticleSpecies>)) {
        rangePtAll[toI(RangeEdge::Min)] = std::min(rangePtAll[toI(RangeEdge::Min)], cgTrack.cfgCutsRangePt.value.get(iParticleSpecies, toI(RangeEdge::Min)));
        rangePtAll[toI(RangeEdge::Max)] = std::max(rangePtAll[toI(RangeEdge::Max)], cgTrack.cfgCutsRangePt.value.get(iParticleSpecies, toI(RangeEdge::Max)));
      }
    }

    rctFlagsChecker.init(cgEvent.cfgLabelFlagsRct.value, static_cast<bool>(cgEvent.cfgFlagsRct.value.get("ZDC")), static_cast<bool>(cgEvent.cfgFlagsRct.value.get("Acceptance")), static_cast<bool>(cgEvent.cfgFlagsRct.value.get("Table")));

    ccdb->setURL(cgCcdb.cfgUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(true);
    if (cgCcdb.cfgTimestampLatest.value >= 0) {
      ccdb->setCreatedNotAfter(cgCcdb.cfgTimestampLatest.value);
    }

    readCcdb<true>();
    std::int32_t nRunsBad{};
    for (const auto& [runNumber, runIndexGroupIndex] : holderCcdb.runNumbersIndicesGroupIndices) {
      const std::int32_t runGroupIndex{runIndexGroupIndex.second};
      if (runGroupIndex == 0 || (cgEvent.cfgFlagRejectionRunBad.value && runGroupIndex < 0)) {
        ++nRunsBad;
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
        if (runIndexGroupIndex.second == 0 || (cgEvent.cfgFlagRejectionRunBad.value && runIndexGroupIndex.second < 0)) {
          LOG(info) << "Enabling rejecting run: " << runNumber << " (" << runIndexGroupIndex.second << ")";
        } else {
          LOG(info) << "Enabling processing run: " << runNumber << " (" << std::abs(runIndexGroupIndex.second) << ")";
        }
      }
    }

    if (cgEvent.cfgLabelFlagsRct.value.empty()) {
      LOG(info) << "No RCT flags label enabled.";
    } else {
      LOG(info) << "Enabling RCT flags label: " << cgEvent.cfgLabelFlagsRct.value;
    }
    if (static_cast<bool>(cgEvent.cfgFlagsRct.value.get("ZDC"))) {
      LOG(info) << "Enabling RCT flag: ZDC";
    }
    if (static_cast<bool>(cgEvent.cfgFlagsRct.value.get("Acceptance"))) {
      LOG(info) << "Enabling RCT flag: acceptance";
    }
    if (static_cast<bool>(cgEvent.cfgFlagsRct.value.get("Table"))) {
      LOG(info) << "Enabling RCT flag: table";
    }

    if ((cgEvent.cfgBitsSelection.value & ((std::uint64_t{1} << aod::evsel::EventSelectionFlags::kNsel) - 1)) == 0) {
      LOG(info) << "No event selection bit enabled.";
    } else {
      for (std::int32_t const& iBit : std::views::iota(0, aod::evsel::EventSelectionFlags::kNsel)) {
        if (static_cast<bool>((cgEvent.cfgBitsSelection.value >> iBit) & 1)) {
          LOG(info) << "Enabling event selection bit: " << aod::evsel::selectionLabels[iBit];
        }
      }
    }

    switch (cgEvent.cfgIndexDefinitionCentrality.value) {
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

    hrCounter.add("hNEvents", ";;No. of Events", {HistType::kTH1D, {{NEs<EventSelection> + aod::evsel::EventSelectionFlags::kNsel, -0.5, static_cast<double>(NEs<EventSelection> + aod::evsel::EventSelectionFlags::kNsel) - 0.5, "Selection"}}});
    if (doProcessMc.value) {
      hrCounter.add("hNEventsMc", ";;No. of MC Events", {HistType::kTH1D, {{NEs<McEventSelection>, -0.5, static_cast<double>(NEs<McEventSelection>) - 0.5, "Selection"}}});
    }

    if (cgAnalysis.cfgFlagQaRun.value) {
      LOG(info) << "Enabling run QA.";

      const HistogramConfigSpec hcsQaRun{HistType::kTProfile, {{static_cast<std::int32_t>(holderCcdb.runNumbersIndicesGroupIndices.size()), -0.5, holderCcdb.runNumbersIndicesGroupIndices.size() - 0.5, "Run Index"}}};

      for (const auto& [name, title, isChargeSpeciesSeparated] : std::to_array<std::tuple<std::string, std::string, bool>>(
             {{"Vx", "#LT#it{V}_{#it{x}}#GT (cm)", false},
              {"Vy", "#LT#it{V}_{#it{y}}#GT (cm)", false},
              {"Vz", "#LT#it{V}_{#it{z}}#GT (cm)", false},
              {"MultiplicityFt0a", "FT0A #LTMultiplicity#GT", false},
              {"MultiplicityFt0c", "FT0C #LTMultiplicity#GT", false},
              {"CentralityNtpv", "NTPV #LTCentrality#GT", false},
              {"CentralityFt0a", "FT0A #LTCentrality#GT", false},
              {"CentralityFt0c", "FT0C #LTCentrality#GT", false},
              {"CentralityFt0m", "FT0M #LTCentrality#GT", false},
              {"NGlobalTracks", "#LTnGlobalTracks#GT", true},
              {"NPvContributors", "#LTnPvContributors#GT", true},
              {std::format("{}Dca{}", getName(DcaMeasure::Mean), getName(DcaAxis::Xy)), "#LT#LTDCA_{#it{xy}}#GT_{event}#GT (cm)", true},
              {std::format("{}Dca{}", getName(DcaMeasure::Sigma), getName(DcaAxis::Xy)), "#LT#it{#sigma}(DCA_{#it{xy}})_{event}#GT (cm)", true},
              {std::format("{}Dca{}", getName(DcaMeasure::Mean), getName(DcaAxis::Z)), "#LT#LTDCA_{#it{z}}#GT_{event}#GT (cm)", true},
              {std::format("{}Dca{}", getName(DcaMeasure::Sigma), getName(DcaAxis::Z)), "#LT#it{#sigma}(DCA_{#it{z}})_{event}#GT (cm)", true},
              {"NTofBeta", "#LTnTofBeta#GT", true},
              {"ItsNCls", "ITS #LTnClusters#GT", true},
              {"ItsChi2NCls", "ITS #LT#it{#chi}^{2}/nClusters#GT", true},
              {"TpcNCls", "TPC #LTnClusters#GT", true},
              {"TpcChi2NCls", "TPC #LT#it{#chi}^{2}/nClusters#GT", true},
              {"TpcNClsSharedRatio", "TPC #LTnSharedClusters/nClusters#GT", true},
              {"TpcNCrossedRows", "TPC #LTnCrossedRows#GT", true},
              {"TpcNCrossedRowsRatio", "TPC #LTnCrossedRows/nFindableClusters#GT", true},
              {"Pt", "#LT#it{p}_{T}#GT (GeV/#it{c})", true},
              {"Eta", "#LT#it{#eta}#GT", true},
              {"Phi", "#LT#it{#varphi}#GT", true},
              {"TpcDeDx", "TPC #LTd#it{E}/d#it{x}#GT (a.u.)", true},
              {std::format("{}NSigma{}", getName(Detector::Tpc), getName(ParticleSpecies::Pion)), std::format("{} #LT#it{{n}}#it{{#sigma}}_{{#pi}}#GT", getName<NameKind::Display>(Detector::Tpc)), true},
              {std::format("{}NSigma{}", getName(Detector::Tpc), getName(ParticleSpecies::Kaon)), std::format("{} #LT#it{{n}}#it{{#sigma}}_{{K}}#GT", getName<NameKind::Display>(Detector::Tpc)), true},
              {std::format("{}NSigma{}", getName(Detector::Tpc), getName(ParticleSpecies::Proton)), std::format("{} #LT#it{{n}}#it{{#sigma}}_{{p}}#GT", getName<NameKind::Display>(Detector::Tpc)), true},
              {"TofInverseBeta", "TOF #LT1/#it{#beta}#GT", true},
              {std::format("{}NSigma{}", getName(Detector::Tof), getName(ParticleSpecies::Pion)), std::format("{} #LT#it{{n}}#it{{#sigma}}_{{#pi}}#GT", getName<NameKind::Display>(Detector::Tof)), true},
              {std::format("{}NSigma{}", getName(Detector::Tof), getName(ParticleSpecies::Kaon)), std::format("{} #LT#it{{n}}#it{{#sigma}}_{{K}}#GT", getName<NameKind::Display>(Detector::Tof)), true},
              {std::format("{}NSigma{}", getName(Detector::Tof), getName(ParticleSpecies::Proton)), std::format("{} #LT#it{{n}}#it{{#sigma}}_{{p}}#GT", getName<NameKind::Display>(Detector::Tof)), true}})) {
        if (!isChargeSpeciesSeparated) {
          hrQaRun.add(std::format("pRunIndex{}", name).c_str(), std::format(";;{}", title).c_str(), hcsQaRun);
        } else {
          hrQaRun.add(std::format("pRunIndex{}_{}", name, getName<NameKind::Lower>(ChargeSpecies::Plus)).c_str(), std::format(";;{} (#it{{q}}>0)", title).c_str(), hcsQaRun);
          hrQaRun.add(std::format("pRunIndex{}_{}", name, getName<NameKind::Lower>(ChargeSpecies::Minus)).c_str(), std::format(";;{} (#it{{q}}<0)", title).c_str(), hcsQaRun);
        }
      }
    }

    if (cgAnalysis.cfgFlagQaEvent.value) {
      LOG(info) << "Enabling event QA.";

      const AxisSpec asMultiplicity{cgEvent.cfgNMultiplicityBinsAxis.value, -0.5, cgEvent.cfgNMultiplicityBinsAxis.value - 0.5};
      const HistogramConfigSpec hcsQaEvent{HistType::kTHnSparseD, {asMultiplicity, asMultiplicity}};

      hrQaEvent.add("hVxVy", "", {HistType::kTHnSparseD, {{150, -0.15, 0.15, "#it{V}_{#it{x}} (cm)"}, {150, -0.15, 0.15, "#it{V}_{#it{y}} (cm)"}}});
      hrQaEvent.add("hVz", "", {HistType::kTH1D, {{300, -15., 15., "#it{V}_{#it{z}} (cm)"}}});
      hrQaEvent.add("hNPvContributorsNGlobalTracks", ";nPvContributors;nGlobalTracks;", hcsQaEvent);
      hrQaEvent.add("hNGlobalTracksMeanDcaXy", ";nGlobalTracks;", {HistType::kTHnSparseD, {asMultiplicity, {250, -0.25, 0.25, "#LTDCA_{#it{xy}}#GT_{event} (cm)"}}});
      hrQaEvent.add("hNGlobalTracksMeanDcaXy_nPvContributorsCut", ";nGlobalTracks;", {HistType::kTHnSparseD, {asMultiplicity, {250, -0.25, 0.25, "#LTDCA_{#it{xy}}#GT_{event} (cm)"}}});
      hrQaEvent.add("hNGlobalTracksMeanDcaZ", ";nGlobalTracks;", {HistType::kTHnSparseD, {asMultiplicity, {200, -2., 2., "#LTDCA_{#it{z}}#GT_{event} (cm)"}}});
      hrQaEvent.add("hNGlobalTracksMeanDcaZ_nPvContributorsCut", ";nGlobalTracks;", {HistType::kTHnSparseD, {asMultiplicity, {200, -2., 2., "#LTDCA_{#it{z}}#GT_{event} (cm)"}}});
      hrQaEvent.add("hNTofBetaNGlobalTracks", ";nTofBeta;nGlobalTracks;", hcsQaEvent);
      hrQaEvent.add("hNTofBetaNGlobalTracks_nPvContributorsCut", ";nTofBeta;nGlobalTracks;", hcsQaEvent);
    }

    if (cgAnalysis.cfgFlagQaCentrality.value) {
      LOG(info) << "Enabling centrality QA.";

      hrQaCentrality.add("hCentralitySelection", "", {HistType::kTHnSparseD, {{100, 0., 100., "Centrality (%)"}, {NEs<EventSelection> + aod::evsel::EventSelectionFlags::kNsel, -0.5, static_cast<double>(NEs<EventSelection> + aod::evsel::EventSelectionFlags::kNsel) - 0.5, "Selection"}}});
      hrQaCentrality.add("hCentralityMultiplicity", "", {HistType::kTHnSparseD, {{100, 0., 100., "Centrality (%)"}, {cgEvent.cfgNMultiplicityBinsAxis.value, -0.5, cgEvent.cfgNMultiplicityBinsAxis.value - 0.5, "Multiplicity"}}});
    }

    if (cgAnalysis.cfgFlagQaTrack.value) {
      LOG(info) << "Enabling track QA.";

      for (const auto& [name, hcs] : std::to_array<std::pair<std::string_view, HistogramConfigSpec>>(
             {{"ItsNClsChi2NCls", {HistType::kTHnSparseD, {{10, -0.5, 9.5, "ITS nClusters"}, {80, 0., 40., "ITS #it{#chi}^{2}/nClusters"}}}},
              {"TpcNClsChi2NCls", {HistType::kTHnSparseD, {{180, -0.5, 179.5, "TPC nClusters"}, {100, 0., 5., "TPC #it{#chi}^{2}/nClusters"}}}},
              {"TpcNClsNClsShared", {HistType::kTHnSparseD, {{180, -0.5, 179.5, "TPC nClusters"}, {180, -0.5, 179.5, "TPC nSharedClusters"}}}},
              {"TpcNClsFindableNCrossedRows", {HistType::kTHnSparseD, {{180, -0.5, 179.5, "TPC nFindableClusters"}, {180, -0.5, 179.5, "TPC nCrossedRows"}}}}})) {
        for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
          hrQaTrack.add(std::format("h{}_{}", name, getName<ChargeSpecies, NameKind::Lower>(iChargeSpecies)).c_str(), "", hcs);
        }
      }
    }

    if (cgAnalysis.cfgFlagQaDca.value) {
      LOG(info) << "Enabling DCA QA.";

      const AxisSpec asPt{40, 0., 2., "#it{p}_{T} (GeV/#it{c})"};
      const HistogramConfigSpec hcsQaDcaProfile{HistType::kTProfile3D, {{cgEvent.cfgAxisCentralityCalibration, "Centrality (%)"}, asPt, {24, -1.2, 1.2, "#it{#eta}"}}};

      for (const auto& [name, title, hcs] : std::to_array<std::tuple<std::string_view, std::string_view, HistogramConfigSpec>>(
             {{"hPtDcaXy", "", {HistType::kTHnSparseD, {asPt, {250, -0.25, 0.25, "DCA_{#it{xy}} (cm)"}}}},
              {"pCentralityPtEtaDcaXy", ";;#LTDCA_{#it{xy}}#GT (cm)", hcsQaDcaProfile},
              {"hPtDcaZ", "", {HistType::kTHnSparseD, {asPt, {250, -0.5, 0.5, "DCA_{#it{z}} (cm)"}}}},
              {"pCentralityPtEtaDcaZ", ";;#LTDCA_{#it{z}}#GT (cm)", hcsQaDcaProfile}})) {
        for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
          hrQaDca.add(std::format("{}_{}", name, getName<ChargeSpecies, NameKind::Lower>(iChargeSpecies)).c_str(), title.data(), hcs);
        }
      }
    }

    for (std::int32_t const& iParticleSpeciesAll : std::views::iota(0, NEs<ParticleSpeciesAll>)) {
      if (!static_cast<bool>(cgAnalysis.cfgFlagsQaAcceptance.value.get(iParticleSpeciesAll))) {
        continue;
      }

      LOG(info) << "Enabling " << getName<ParticleSpeciesAll, NameKind::DisplayLower>(iParticleSpeciesAll) << " acceptance QA.";

      for (std::int32_t const& iPidStrategy : std::views::iota(0, NEs<PidStrategy>)) {
        for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
          hrQaAcceptance.add(std::format("h{}Pt_{}Edge{}{}", iParticleSpeciesAll == toI(ParticleSpeciesAll::All) ? "Eta" : "Rapidity", getName<PidStrategy, NameKind::Lower>(iPidStrategy), getName<ParticleSpeciesAll>(iParticleSpeciesAll), getName<ChargeSpecies>(iChargeSpecies)).c_str(), "", {HistType::kTHnSparseD, {{300, -1.5, 1.5, iParticleSpeciesAll == toI(ParticleSpeciesAll::All) ? "#it{#eta}" : "#it{y}"}, {250, 0., 2.5, "#it{p}_{T} (GeV/#it{c})"}}});
        }
      }
    }

    for (std::int32_t const& iParticleSpeciesAll : std::views::iota(0, NEs<ParticleSpeciesAll>)) {
      if (!static_cast<bool>(cgAnalysis.cfgFlagsQaPhi.value.get(iParticleSpeciesAll))) {
        continue;
      }

      LOG(info) << "Enabling " << getName<ParticleSpeciesAll, NameKind::DisplayLower>(iParticleSpeciesAll) << " phi QA.";

      const HistogramConfigSpec hcsQaPhi{HistType::kTHnSparseF, {{cgEvent.cfgAxisCentralityCalibration, "Centrality (%)"}, {20, 0., 2., "#it{p}_{T} (GeV/#it{c})"}, {24, -1.2, 1.2, "#it{#eta}"}, {360, 0., constants::math::TwoPI, "#it{#varphi} (rad)"}}};

      for (std::int32_t const& iPidStrategy : std::views::iota(0, NEs<PidStrategy>)) {
        for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
          hrQaPhi.add(std::format("hCentralityPtEtaPhi_{}{}{}{}", doProcessMc.value ? "mc" : "", doProcessMc.value ? getName<PidStrategy>(iPidStrategy) : getName<PidStrategy, NameKind::Lower>(iPidStrategy), getName<ParticleSpeciesAll>(iParticleSpeciesAll), getName<ChargeSpecies>(iChargeSpecies)).c_str(), "", hcsQaPhi);
        }
      }
    }

    for (std::int32_t const& iParticleSpeciesAll : std::views::iota(0, NEs<ParticleSpeciesAll>)) {
      if (!static_cast<bool>(cgAnalysis.cfgFlagsQaPid.value.get(iParticleSpeciesAll))) {
        continue;
      }

      LOG(info) << "Enabling " << getName<ParticleSpeciesAll, NameKind::DisplayLower>(iParticleSpeciesAll) << " PID QA.";

      const AxisSpec asCentrality{cgEvent.cfgAxisCentralityCalibration, "Centrality (%)"};

      if (iParticleSpeciesAll == toI(ParticleSpeciesAll::All)) {
        const AxisSpec asPOverQ{350, -3.5, 3.5, "#it{p}/#it{q} (GeV/#it{c})"};
        const AxisSpec asEta{48, -1.2, 1.2, "#it{#eta}"};

        hrQaPid.add("hCentralityPOverQEtaTpcLnDeDx", "", {HistType::kTHnSparseF, {asCentrality, asPOverQ, asEta, {240, 3., 9., "TPC ln(d#it{E}/d#it{x} (a.u.))"}}});
        hrQaPid.add("hCentralityPOverQEtaTofInverseBeta", "", {HistType::kTHnSparseF, {asCentrality, asPOverQ, asEta, {120, 0.5, 3.5, "TOF 1/#it{#beta}"}}});
      } else {
        const HistogramConfigSpec hcsQaPid{HistType::kTHnSparseF, {asCentrality, {40, 0., 2., "#it{p}_{T} (GeV/#it{c})"}, {32, -0.8, 0.8, "#it{#eta}"}, {300, -30., 30.}}};

        if (doProcessMc.value) {
          for (std::int32_t const& iDetector : std::views::iota(0, NEs<Detector>)) {
            for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
              hrQaPid.add(std::format("hCentralityPtEta{}NSigma{}_mc{}{}", getName<Detector>(iDetector), getName<ParticleSpeciesAll>(iParticleSpeciesAll), getName<ParticleSpeciesAll>(iParticleSpeciesAll), getName<ChargeSpecies>(iChargeSpecies)).c_str(), std::format(";;;;{} #it{{n}}#it{{#sigma}}_{{{}}};", getName<Detector, NameKind::Display>(iDetector), getTitle<ParticleSpeciesAll>(iParticleSpeciesAll)).c_str(), hcsQaPid);
            }
          }
        } else {
          for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
            hrQaPid.add(std::format("hCentralityPtEta{}NSigma{}_{}", getName(Detector::Tpc), getName<ParticleSpeciesAll>(iParticleSpeciesAll), getName<ChargeSpecies, NameKind::Lower>(iChargeSpecies)).c_str(), std::format(";;;;{} #it{{n}}#it{{#sigma}}_{{{}}};", getName<NameKind::Display>(Detector::Tpc), getTitle<ParticleSpeciesAll>(iParticleSpeciesAll)).c_str(), hcsQaPid);
          }
          for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
            hrQaPid.add(std::format("hCentralityPtEta{}NSigma{}_{}{}{}", getName(Detector::Tpc), getName<ParticleSpeciesAll>(iParticleSpeciesAll), getName<NameKind::Lower>(Detector::Tof), getName<ParticleSpeciesAll>(iParticleSpeciesAll), getName<ChargeSpecies>(iChargeSpecies)).c_str(), std::format(";;;;{} #it{{n}}#it{{#sigma}}_{{{}}};", getName<NameKind::Display>(Detector::Tpc), getTitle<ParticleSpeciesAll>(iParticleSpeciesAll)).c_str(), hcsQaPid);
          }
          for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
            hrQaPid.add(std::format("hCentralityPtEta{}NSigma{}_{}{}{}", getName(Detector::Tof), getName<ParticleSpeciesAll>(iParticleSpeciesAll), getName<NameKind::Lower>(Detector::Tpc), getName<ParticleSpeciesAll>(iParticleSpeciesAll), getName<ChargeSpecies>(iChargeSpecies)).c_str(), std::format(";;;;{} #it{{n}}#it{{#sigma}}_{{{}}};", getName<NameKind::Display>(Detector::Tof), getTitle<ParticleSpeciesAll>(iParticleSpeciesAll)).c_str(), hcsQaPid);
          }
          for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
            hrQaPid.add(std::format("hCentralityPtEta{}NSigma{}_{}", getName(PidStrategy::TpcTof), getName<ParticleSpeciesAll>(iParticleSpeciesAll), getName<ChargeSpecies, NameKind::Lower>(iChargeSpecies)).c_str(), std::format(";;;;{} #it{{n}}#it{{#sigma}}_{{{}}};", getName<NameKind::Display>(PidStrategy::TpcTof), getTitle<ParticleSpeciesAll>(iParticleSpeciesAll)).c_str(), hcsQaPid);
          }
        }
      }
    }

    if (doProcessMc.value) {
      if (cgAnalysis.cfgFlagQaMc.value) {
        LOG(info) << "Enabling MC QA.";

        const double maxAbsVz{std::ceil(cgEvent.cfgFlagMcCollisionVz.value ? cgEvent.cfgCutMaxAbsVzMc.value : cgEvent.cfgCutMaxAbsVz.value)};
        const AxisSpec asCentrality{20, 0., 100., "Centrality (%)"};

        hrQaMc.add("hNCollisionsPerMcCollision", "", {HistType::kTH1D, {{20, -0.5, 19.5, "No. of collisions per MC collision"}}});
        hrQaMc.add("hCentralityNTracksPerMcParticle", "", {HistType::kTHnSparseF, {asCentrality, {20, -0.5, 19.5, "No. of tracks per MC particle"}}});
        hrQaMc.add("hCentralityVzMcDeltaVz", "", {HistType::kTHnSparseF, {asCentrality, {static_cast<std::int32_t>(maxAbsVz) * 20, -maxAbsVz, maxAbsVz, "#it{V}_{#it{z}}^{Gen} (cm)"}, {200, -0.2, 0.2, "#it{V}_{#it{z}}^{Rec}#minus#it{V}_{#it{z}}^{Gen} (cm)"}}});
        hrQaMc.add("hCentralityPtMcEtaMcDeltaPt", "", {HistType::kTHnSparseF, {asCentrality, {200, 0., 2., "#it{p}_{T}^{Gen} (GeV/#it{c})"}, {24, -1.2, 1.2, "#it{#eta}_{Gen}"}, {320, -0.8, 0.8, "#it{p}_{T}^{Rec}#minus#it{p}_{T}^{Gen} (GeV/#it{c})"}}});
        hrQaMc.add("hCentralityPtMcEtaMcDeltaEta", "", {HistType::kTHnSparseF, {asCentrality, {20, 0., 2., "#it{p}_{T}^{Gen} (GeV/#it{c})"}, {240, -1.2, 1.2, "#it{#eta}_{Gen}"}, {160, -0.4, 0.4, "#it{#eta}_{Rec}#minus#it{#eta}_{Gen}"}}});
      }
    }

    for (std::int32_t const& iParticleSpecies : std::views::iota(0, NEs<ParticleSpecies>)) {
      if (!static_cast<bool>(cgAnalysis.cfgFlagsCalculationYield.value.get(iParticleSpecies))) {
        continue;
      }
      LOG(info) << "Enabling " << getName<ParticleSpecies, NameKind::DisplayLower>(iParticleSpecies) << " yield calculation.";

      const double maxAbsVz{std::ceil(doProcessMc.value && cgEvent.cfgFlagMcCollisionVz.value ? cgEvent.cfgCutMaxAbsVzMc.value : cgEvent.cfgCutMaxAbsVz.value)};
      const HistogramConfigSpec hcsCalculationYield{HistType::kTHnSparseF, {{static_cast<std::int32_t>(maxAbsVz) * 2, -maxAbsVz, maxAbsVz, "#it{V}_{#it{z}} (cm)"}, {cgEvent.cfgAxisCentralityCalibration, "Centrality (%)"}, {40, 0., 2., "#it{p}_{T} (GeV/#it{c})"}, {32, -0.8, 0.8, "#it{#eta}"}}};

      if (doProcessMc.value) {
        for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
          hrCalculationYield.add(std::format("hVzCentralityPtMcEtaMc_mc{}{}", getName<ParticleSpecies>(iParticleSpecies), getName<ChargeSpecies>(iChargeSpecies)).c_str(), "", hcsCalculationYield);
        }
        for (std::int32_t const& iPidStrategy : std::views::iota(0, NEs<PidStrategy>)) {
          for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
            if (cgTrack.cfgFlagMcParticleMomentum.value) {
              hrCalculationYield.add(std::format("hVzCentralityPtMcEtaMc_mc{}{}{}", getName<PidStrategy>(iPidStrategy), getName<ParticleSpecies>(iParticleSpecies), getName<ChargeSpecies>(iChargeSpecies)).c_str(), "", hcsCalculationYield);
            } else {
              hrCalculationYield.add(std::format("hVzCentralityPtEta_mc{}{}{}", getName<PidStrategy>(iPidStrategy), getName<ParticleSpecies>(iParticleSpecies), getName<ChargeSpecies>(iChargeSpecies)).c_str(), "", hcsCalculationYield);
            }
          }
        }
      } else {
        for (std::int32_t const& iPidStrategy : std::views::iota(0, NEs<PidStrategy>)) {
          for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
            hrCalculationYield.add(std::format("hVzCentralityPtEta_{}{}{}", getName<PidStrategy, NameKind::Lower>(iPidStrategy), getName<ParticleSpecies>(iParticleSpecies), getName<ChargeSpecies>(iChargeSpecies)).c_str(), "", hcsCalculationYield);
          }
        }
      }
    }

    if (doProcessMc.value) {
      for (std::int32_t const& iParticleSpecies : std::views::iota(0, NEs<ParticleSpecies>)) {
        if (!static_cast<bool>(cgAnalysis.cfgFlagsCalculationPurity.value.get(iParticleSpecies))) {
          continue;
        }

        LOG(info) << "Enabling " << getName<ParticleSpecies, NameKind::DisplayLower>(iParticleSpecies) << " purity calculation.";

        const HistogramConfigSpec hcsCalculationPurity{HistType::kTProfile3D, {{cgEvent.cfgAxisCentralityCalibration, "Centrality (%)"}, {20, 0., 2., "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#it{#eta}"}}};

        for (std::int32_t const& iPidStrategy : std::views::iota(0, NEs<PidStrategy>)) {
          for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
            hrCalculationPurity.add(std::format("pCentralityPtEtaPurity{}{}{}", getName<PidStrategy>(iPidStrategy), getName<ParticleSpecies>(iParticleSpecies), getName<ChargeSpecies>(iChargeSpecies)).c_str(), "", hcsCalculationPurity);
          }
        }
      }
    }

    if (doProcessMc.value) {
      for (std::int32_t const& iParticleSpecies : std::views::iota(0, NEs<ParticleSpecies>)) {
        if (!static_cast<bool>(cgAnalysis.cfgFlagsCalculationFractionPrimary.value.get(iParticleSpecies))) {
          continue;
        }

        LOG(info) << "Enabling " << getName<ParticleSpecies, NameKind::DisplayLower>(iParticleSpecies) << " primary fraction calculation.";

        const HistogramConfigSpec hcsCalculationFractionPrimary{HistType::kTProfile3D, {{cgEvent.cfgAxisCentralityCalibration, "Centrality (%)"}, {20, 0., 2., "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#it{#eta}"}}};

        for (std::int32_t const& iPidStrategy : std::views::iota(0, NEs<PidStrategy>)) {
          for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
            hrCalculationFractionPrimary.add(std::format("pCentralityPtEtaFractionPrimary{}{}{}", getName<PidStrategy>(iPidStrategy), getName<ParticleSpecies>(iParticleSpecies), getName<ChargeSpecies>(iChargeSpecies)).c_str(), "", hcsCalculationFractionPrimary);
          }
        }
      }
    }

    for (std::int32_t const& iParticleNumber : std::views::iota(0, NEs<ParticleNumber>)) {
      if (!static_cast<bool>(cgAnalysis.cfgFlagsCalculationFluctuation.value.get(iParticleNumber))) {
        continue;
      }

      LOG(info) << "Enabling " << getName<ParticleNumber, NameKind::DisplayLower>(iParticleNumber) << " number fluctuation calculation.";

      const AxisSpec asCentrality{cgEvent.cfgAxisCentrality, "Centrality (%)"};
      const AxisSpec asMultiplicity{cgEvent.cfgNMultiplicityBinsAxis.value, -0.5, cgEvent.cfgNMultiplicityBinsAxis.value - 0.5};
      const HistogramConfigSpec hcsDistribution{HistType::kTHnSparseD, {asCentrality, asMultiplicity, asMultiplicity}};
      const HistogramConfigSpec hcsFluctuationCalculator{HistType::kTH3D, {asCentrality, {cgEvent.cfgNSubgroups.value, -0.5, cgEvent.cfgNSubgroups.value - 0.5, "Subgroup Index"}, {fluctuation_calculator_base::NOrderKeys, -0.5, fluctuation_calculator_base::NOrderKeys - 0.5, "Order Key Index"}}};

      for (std::int32_t const& iChargeNumber : std::views::iota(0, NEs<ChargeNumber>)) {
        if (doProcessMc.value) {
          fluctuationCalculatorsTrackMcParticle[iParticleNumber][iChargeNumber] = std::make_unique<FluctuationCalculatorTrack>();
        }
        fluctuationCalculatorsTrackTrack[iParticleNumber][iChargeNumber] = std::make_unique<FluctuationCalculatorTrack>();
      }

      if (doProcessMc.value) {
        hrCalculationFluctuation.add(std::format("hCentralityN{}{}N{}{}_mc", getName<ParticleNumber>(iParticleNumber), getName(ChargeSpecies::Plus), getName<ParticleNumber>(iParticleNumber), getName(ChargeSpecies::Minus)).c_str(), std::format(";;#it{{N}}({}^{{{}}});#it{{N}}({}^{{{}}});", getTitle<ParticleNumber>(iParticleNumber), getTitle<ChargeSpecies>(toI(ChargeSpecies::Plus)), getTitle<ParticleNumber>(iParticleNumber), getTitle<ChargeSpecies>(toI(ChargeSpecies::Minus))).c_str(), hcsDistribution);
        hrCalculationFluctuation.add(std::format("hCentralityN{}{}N{}{}_mcEff", getName<ParticleNumber>(iParticleNumber), getName(ChargeSpecies::Plus), getName<ParticleNumber>(iParticleNumber), getName(ChargeSpecies::Minus)).c_str(), std::format(";;#it{{N}}({}^{{{}}});#it{{N}}({}^{{{}}});", getTitle<ParticleNumber>(iParticleNumber), getTitle<ChargeSpecies>(toI(ChargeSpecies::Plus)), getTitle<ParticleNumber>(iParticleNumber), getTitle<ChargeSpecies>(toI(ChargeSpecies::Minus))).c_str(), hcsDistribution);
        for (std::int32_t const& iChargeNumber : std::views::iota(0, NEs<ChargeNumber>)) {
          hrCalculationFluctuation.add(std::format("hFluctuationCalculator{}{}_mc", getName<ParticleNumber>(iParticleNumber), getName<ChargeNumber>(iChargeNumber)).c_str(), "", hcsFluctuationCalculator);
        }
      }
      hrCalculationFluctuation.add(std::format("hCentralityN{}{}N{}{}", getName<ParticleNumber>(iParticleNumber), getName(ChargeSpecies::Plus), getName<ParticleNumber>(iParticleNumber), getName(ChargeSpecies::Minus)).c_str(), std::format(";;#it{{N}}({}^{{{}}});#it{{N}}({}^{{{}}});", getTitle<ParticleNumber>(iParticleNumber), getTitle<ChargeSpecies>(toI(ChargeSpecies::Plus)), getTitle<ParticleNumber>(iParticleNumber), getTitle<ChargeSpecies>(toI(ChargeSpecies::Minus))).c_str(), hcsDistribution);
      for (std::int32_t const& iChargeNumber : std::views::iota(0, NEs<ChargeNumber>)) {
        hrCalculationFluctuation.add(std::format("hFluctuationCalculator{}{}", getName<ParticleNumber>(iParticleNumber), getName<ChargeNumber>(iChargeNumber)).c_str(), "", hcsFluctuationCalculator);
      }
    }
  }

  template <bool DoInit>
  void readCcdb()
  {
    if constexpr (DoInit) {
      holderCcdb.lCcdb = ccdb->get<TList>(cgCcdb.cfgPath.value);
      if (!holderCcdb.lCcdb || holderCcdb.lCcdb->IsA() != TList::Class()) {
        LOG(fatal) << "Invalid CCDB object!";
      }

      const TGraph* const gRunNumberGroupIndex{dynamic_cast<const TGraph*>(holderCcdb.lCcdb->FindObject("gRunNumberGroupIndex"))};
      if (!gRunNumberGroupIndex || gRunNumberGroupIndex->IsA() != TGraph::Class()) {
        LOG(fatal) << "Invalid gRunNumberGroupIndex!";
      }
      for (std::int32_t const& iRun : std::views::iota(0, gRunNumberGroupIndex->GetN())) {
        holderCcdb.runNumbersIndicesGroupIndices[static_cast<std::int32_t>(std::rint(gRunNumberGroupIndex->GetX()[iRun]))] = {iRun, static_cast<std::int32_t>(std::rint(gRunNumberGroupIndex->GetY()[iRun]))};
      }

      if (cgEvent.cfgFlagRejectionRunBadMc.value) {
        const TGraph* const gRunNumberGroupIndexMc{dynamic_cast<const TGraph*>(holderCcdb.lCcdb->FindObject("gRunNumberGroupIndex_mc"))};
        if (!gRunNumberGroupIndexMc || gRunNumberGroupIndexMc->IsA() != TGraph::Class()) {
          LOG(fatal) << "Invalid gRunNumberGroupIndex_mc!";
        }
        for (std::int32_t const& iRun : std::views::iota(0, gRunNumberGroupIndexMc->GetN())) {
          if (static_cast<std::int32_t>(std::rint(gRunNumberGroupIndexMc->GetY()[iRun])) <= 0) {
            if (const auto iter{holderCcdb.runNumbersIndicesGroupIndices.find(static_cast<std::int32_t>(std::rint(gRunNumberGroupIndexMc->GetX()[iRun])))}; iter != holderCcdb.runNumbersIndicesGroupIndices.end() && iter->second.second > 0) {
              iter->second.second = -iter->second.second;
            }
          }
        }
      }
    } else {
      const std::int32_t runGroupIndex{std::abs(holderEvent.runGroupIndex)};
      if (holderCcdb.runGroupIndexCurrent == runGroupIndex) {
        return;
      }

      holderCcdb.runGroupIndexCurrent = runGroupIndex;
      holderCcdb.calibrationsPtMeasureDca = {};
      holderCcdb.hsCentralityPtEtaShiftNSigmaPid = {};
      holderCcdb.hsVzCentralityPtEtaEfficiency = {};
      if (!cgTrack.cfgFlagRecalibrationDca.value && !isEnabled(cgTrack.cfgFlagsRecalibrationNSigmaPid) && !doCalculationFluctuation) {
        return;
      }

      const std::string nameList{std::format("lRunGroup_{}", runGroupIndex)};
      const TList* const lRunGroup{dynamic_cast<const TList*>(holderCcdb.lCcdb->FindObject(nameList.c_str()))};
      if (!lRunGroup) {
        LOG(fatal) << "Invalid " << nameList << "!";
      }

      if (cgTrack.cfgFlagRecalibrationDca.value) {
        for (std::int32_t const& iDcaMeasure : std::views::iota(0, NEs<DcaMeasure>)) {
          for (std::int32_t const& iDcaAxis : std::views::iota(0, NEs<DcaAxis>)) {
            for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
              std::pair<const TFormula*, const TH3*>& calibration{holderCcdb.calibrationsPtMeasureDca[iDcaMeasure][iDcaAxis][iChargeSpecies]};
              const std::string nameFormula{std::format("fPt{}Dca{}{}{}_runGroup{}", getName<DcaMeasure>(iDcaMeasure), getName<DcaAxis>(iDcaAxis), getName<ChargeSpecies>(iChargeSpecies), doProcessMc.value ? "_mc" : "", runGroupIndex)};
              calibration.first = dynamic_cast<const TFormula*>(lRunGroup->FindObject(nameFormula.c_str()));
              if (!calibration.first || calibration.first->GetNdim() != 1 || calibration.first->GetNpar() <= 0) {
                LOG(fatal) << "Invalid " << nameFormula << "!";
              }
              LOG(info) << "Reading from CCDB: " << nameFormula << " \"" << calibration.first->GetExpFormula() << "\"";
              const std::int32_t nParameters{calibration.first->GetNpar()};

              const std::string nameHistogram{std::format("hCentralityEtaParameterPt{}Dca{}{}{}_runGroup{}", getName<DcaMeasure>(iDcaMeasure), getName<DcaAxis>(iDcaAxis), getName<ChargeSpecies>(iChargeSpecies), doProcessMc.value ? "_mc" : "", runGroupIndex)};
              calibration.second = dynamic_cast<const TH3*>(lRunGroup->FindObject(nameHistogram.c_str()));
              if (calibration.second == nullptr || calibration.second->GetNbinsZ() != nParameters || std::ranges::any_of(std::views::iota(0, nParameters), [zAxis = calibration.second->GetZaxis()](const std::int32_t binIndex) -> bool { return zAxis->GetBinCenter(binIndex + 1) != binIndex; })) {
                LOG(fatal) << "Invalid " << nameHistogram << "!";
              }
              LOG(info) << "Reading from CCDB: " << nameHistogram;
            }
          }
        }
      }

      for (std::int32_t const& iParticleSpecies : std::views::iota(0, NEs<ParticleSpecies>)) {
        if (!static_cast<bool>(cgTrack.cfgFlagsRecalibrationNSigmaPid.value.get(iParticleSpecies))) {
          continue;
        }

        for (std::int32_t const& iDetector : std::views::iota(0, NEs<Detector>)) {
          for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
            const std::string name{std::format("hCentralityPtEtaShift{}NSigma{}{}{}_runGroup{}", getName<Detector>(iDetector), getName<ParticleSpecies>(iParticleSpecies), getName<ChargeSpecies>(iChargeSpecies), doProcessMc.value ? "_mc" : "", runGroupIndex)};
            const TH3*& histogram{holderCcdb.hsCentralityPtEtaShiftNSigmaPid[iDetector][iParticleSpecies][iChargeSpecies]};
            histogram = dynamic_cast<const TH3*>(lRunGroup->FindObject(name.c_str()));
            if (!histogram) {
              LOG(fatal) << "Invalid " << name << "!";
            }
            LOG(info) << "Reading from CCDB: " << name;
          }
        }
      }

      for (std::int32_t const& iParticleSpecies : std::views::iota(0, NEs<ParticleSpecies>)) {
        if (!static_cast<bool>(cgAnalysis.cfgFlagsCalculationFluctuation.value.get(toI(ParticleNumber::Charge))) && (iParticleSpecies != toI(ParticleSpecies::Kaon) || !static_cast<bool>(cgAnalysis.cfgFlagsCalculationFluctuation.value.get(toI(ParticleNumber::Kaon)))) && (iParticleSpecies != toI(ParticleSpecies::Proton) || !static_cast<bool>(cgAnalysis.cfgFlagsCalculationFluctuation.value.get(toI(ParticleNumber::Proton))))) {
          continue;
        }

        for (std::int32_t const& iPidStrategy : std::views::iota(0, NEs<PidStrategy>)) {
          for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
            const std::string name{std::format("hVzCentralityPtEtaEfficiency{}{}{}_runGroup{}", getName<PidStrategy>(iPidStrategy), getName<ParticleSpecies>(iParticleSpecies), getName<ChargeSpecies>(iChargeSpecies), runGroupIndex)};
            const THnBase*& histogram{holderCcdb.hsVzCentralityPtEtaEfficiency[iPidStrategy][iParticleSpecies][iChargeSpecies]};
            histogram = dynamic_cast<const THnBase*>(lRunGroup->FindObject(name.c_str()));
            if (!histogram || histogram->GetNdimensions() != HolderCcdb::NDimensionsEfficiency) {
              LOG(fatal) << "Invalid " << name << "!";
            }
            LOG(info) << "Reading from CCDB: " << name;
          }
        }
      }
    }
  }

  template <PidStrategy PidStrategyValue, ParticleSpecies ParticleSpeciesValue, ChargeSpecies ChargeSpeciesValue>
    requires IsValidEnumValue<PidStrategyValue, ParticleSpeciesValue, ChargeSpeciesValue>
  double getEfficiency(const bool doUseMcParticleMomentum) const
  {
    const THnBase* const hsVzCentralityPtEtaEfficiency{holderCcdb.hsVzCentralityPtEtaEfficiency[toI(PidStrategyValue)][toI(ParticleSpeciesValue)][toI(ChargeSpeciesValue)]};
    return hsVzCentralityPtEtaEfficiency ? hsVzCentralityPtEtaEfficiency->GetBinContent(hsVzCentralityPtEtaEfficiency->GetBin(std::array<double, HolderCcdb::NDimensionsEfficiency>{doProcessMc.value && cgEvent.cfgFlagMcCollisionVz.value ? holderMcEvent.vz : holderEvent.vz, holderEvent.centrality, doUseMcParticleMomentum ? holderMcParticle.pt : holderTrack.pt, doUseMcParticleMomentum ? holderMcParticle.eta : holderTrack.eta}.data())) : 0.;
  }

  template <bool DoRecalibrate, Detector DetectorValue, ParticleSpecies ParticleSpeciesValue>
    requires IsValidEnumValue<ParticleSpeciesValue, DetectorValue>
  double getShiftNSigmaPid() const
  {
    if constexpr (DoRecalibrate) {
      if (cgTrack.cfgFlagsRecalibrationNSigmaPid.value.get(toI(ParticleSpeciesValue))) {
        return interpolate(holderCcdb.hsCentralityPtEtaShiftNSigmaPid[toI(DetectorValue)][toI(ParticleSpeciesValue)][holderTrack.sign > 0 ? toI(ChargeSpecies::Plus) : toI(ChargeSpecies::Minus)], holderEvent.centralityCalibration, holderTrack.pt, holderTrack.eta);
      }
    }
    return 0.;
  }

  double getMeasureDca(const std::pair<const TFormula*, const TH3*>& calibration) const
  {
    static thread_local std::vector<double> parametersPtMeasureDcaScratch{};

    const TFormula* const fPtMeasureDca{calibration.first};
    const TH3* const hCentralityEtaParameterPtMeasureDca{calibration.second};
    const std::int32_t nParameters{fPtMeasureDca->GetNpar()};
    if (std::cmp_less(parametersPtMeasureDcaScratch.size(), nParameters)) {
      parametersPtMeasureDcaScratch.resize(nParameters);
    }
    const std::int32_t centralityBinIndex{std::clamp(hCentralityEtaParameterPtMeasureDca->GetXaxis()->FindFixBin(holderEvent.centralityCalibration), 1, hCentralityEtaParameterPtMeasureDca->GetNbinsX())};
    const std::int32_t etaBinIndex{std::clamp(hCentralityEtaParameterPtMeasureDca->GetYaxis()->FindFixBin(holderTrack.eta), 1, hCentralityEtaParameterPtMeasureDca->GetNbinsY())};
    for (std::int32_t const& iParameter : std::views::iota(0, nParameters)) {
      parametersPtMeasureDcaScratch[iParameter] = hCentralityEtaParameterPtMeasureDca->GetBinContent(centralityBinIndex, etaBinIndex, iParameter + 1);
    }
    return fPtMeasureDca->EvalPar(&holderTrack.pt, parametersPtMeasureDcaScratch.data());
  }

  template <DcaAxis DcaAxisValue>
    requires IsValidEnumValue<DcaAxisValue>
  double getAbsNSigmaDca() const
  {
    if (!cgTrack.cfgFlagRecalibrationDca.value) {
      return std::abs(holderTrack.dca[toI(DcaAxisValue)]);
    }

    const std::int32_t chargeSpeciesIndex{holderTrack.sign > 0 ? toI(ChargeSpecies::Plus) : toI(ChargeSpecies::Minus)};
    const double sigma{getMeasureDca(holderCcdb.calibrationsPtMeasureDca[toI(DcaMeasure::Sigma)][toI(DcaAxisValue)][chargeSpeciesIndex])};
    return sigma > 0. ? std::abs((holderTrack.dca[toI(DcaAxisValue)] - getMeasureDca(holderCcdb.calibrationsPtMeasureDca[toI(DcaMeasure::Mean)][toI(DcaAxisValue)][chargeSpeciesIndex])) / sigma) : std::numeric_limits<double>::infinity();
  }

  template <PidStrategyAll PidStrategyAllValue, ParticleSpeciesAll ParticleSpeciesAllValue, Selection SelectionValue = Selection::Default>
    requires IsValidEnumValue<ParticleSpeciesAllValue, PidStrategyAllValue, SelectionValue>
  bool isPid(const bool doRejectOthers) const
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
        if (!(std::abs(holderTrack.nSigmaPid[toI(Detector::Tpc)][ParticleSpeciesIndex]) < cgTrack.cfgCutsMaxAbsNSigmaPid.value.get(ParticleSpeciesIndex, toI(SelectionValue)))) {
          return false;
        }
        if (!(std::abs(holderTrack.nSigmaPid[toI(Detector::Tof)][ParticleSpeciesIndex]) < cgTrack.cfgCutsMaxAbsNSigmaPid.value.get(ParticleSpeciesIndex, toI(SelectionValue)))) {
          return false;
        }
        if (doRejectOthers && !(std::abs(holderTrack.nSigmaPid[toI(Detector::Tof)][ParticleSpeciesIndex]) < std::min(std::abs(holderTrack.nSigmaPid[toI(Detector::Tof)][(ParticleSpeciesIndex + 1) % NEs<ParticleSpecies>]), std::abs(holderTrack.nSigmaPid[toI(Detector::Tof)][(ParticleSpeciesIndex + 2) % NEs<ParticleSpecies>])))) {
          return false;
        }
      } else if constexpr (PidStrategyAllValue == PidStrategyAll::TpcTofCombined) {
        const double absNSigmaPidCombined{std::abs(holderTrack.getNSigmaPidCombined(ParticleSpeciesIndex))};
        if (!(absNSigmaPidCombined < cgTrack.cfgCutsMaxAbsNSigmaPid.value.get(ParticleSpeciesIndex, toI(SelectionValue)))) {
          return false;
        }
        if (doRejectOthers && !(absNSigmaPidCombined < std::min(std::abs(holderTrack.getNSigmaPidCombined((ParticleSpeciesIndex + 1) % NEs<ParticleSpecies>)), std::abs(holderTrack.getNSigmaPidCombined((ParticleSpeciesIndex + 2) % NEs<ParticleSpecies>))))) {
          return false;
        }
      } else {
        constexpr std::int32_t DetectorIndex{toI(getValue<Detector>(PidStrategyAllValue))};
        if (!(std::abs(holderTrack.nSigmaPid[DetectorIndex][ParticleSpeciesIndex]) < cgTrack.cfgCutsMaxAbsNSigmaPid.value.get(ParticleSpeciesIndex, toI(SelectionValue)))) {
          return false;
        }
        if (doRejectOthers && !(std::abs(holderTrack.nSigmaPid[DetectorIndex][ParticleSpeciesIndex]) < std::min(std::abs(holderTrack.nSigmaPid[DetectorIndex][(ParticleSpeciesIndex + 1) % NEs<ParticleSpecies>]), std::abs(holderTrack.nSigmaPid[DetectorIndex][(ParticleSpeciesIndex + 2) % NEs<ParticleSpecies>])))) {
          return false;
        }
      }
    }
    return true;
  }

  template <ParticleSpeciesAll ParticleSpeciesAllValue, ChargeSpecies ChargeSpeciesValue>
    requires IsValidEnumValue<ParticleSpeciesAllValue, ChargeSpeciesValue>
  bool isPid() const
  {
    if constexpr (ParticleSpeciesAllValue == ParticleSpeciesAll::All) {
      return ChargeSpeciesValue == ChargeSpecies::Plus ? holderMcParticle.charge > 0 : holderMcParticle.charge < 0;
    } else {
      return holderMcParticle.pdgCode == getPdgCode(getValue<ParticleSpecies>(ParticleSpeciesAllValue), ChargeSpeciesValue);
    }
  }

  template <ParticleSpeciesAll ParticleSpeciesAllValue>
    requires IsValidEnumValue<ParticleSpeciesAllValue>
  bool isGoodMomentum(const bool doUseMcParticleMomentum) const
  {
    const double pt{doUseMcParticleMomentum ? holderMcParticle.pt : holderTrack.pt};
    if constexpr (ParticleSpeciesAllValue == ParticleSpeciesAll::All) {
      if (!(rangePtAll[toI(RangeEdge::Min)] < pt) || !(pt < rangePtAll[toI(RangeEdge::Max)])) {
        return false;
      }
    } else {
      const std::int32_t particleSpeciesIndex{toI(getValue<ParticleSpecies>(ParticleSpeciesAllValue))};
      if (!(cgTrack.cfgCutsRangePt.value.get(particleSpeciesIndex, toI(RangeEdge::Min)) < pt) || !(pt < cgTrack.cfgCutsRangePt.value.get(particleSpeciesIndex, toI(RangeEdge::Max)))) {
        return false;
      }
    }
    return std::abs(doUseMcParticleMomentum ? holderMcParticle.eta : holderTrack.eta) < cgTrack.cfgCutMaxAbsEta.value;
  }

  bool isGoodDca() const
  {
    return getAbsNSigmaDca<DcaAxis::Xy>() < cgTrack.cfgCutsMaxAbsNSigmaDca.value.get(toI(DcaAxis::Xy), toI(Selection::Default)) && getAbsNSigmaDca<DcaAxis::Z>() < cgTrack.cfgCutsMaxAbsNSigmaDca.value.get(toI(DcaAxis::Z), toI(Selection::Default));
  }

  template <typename T>
  bool isGoodTrack(const T& track) const
  {
    if (cgTrack.cfgFlagPvContributor.value && !track.isPVContributor()) {
      return false;
    }
    if (!(track.itsNCls() > cgTrack.cfgCutsMinItsNCls.value.get(toI(Selection::Default)))) {
      return false;
    }
    if (!(track.itsChi2NCl() < cgTrack.cfgCutsMaxItsChi2NCls.value.get(toI(Selection::Default)))) {
      return false;
    }
    if (!(track.tpcNClsFound() > cgTrack.cfgCutsMinTpcNCls.value.get(toI(Selection::Default)))) {
      return false;
    }
    if (!(cgTrack.cfgCutMinTpcChi2NCls.value < track.tpcChi2NCl()) || !(track.tpcChi2NCl() < cgTrack.cfgCutsMaxTpcChi2NCls.value.get(toI(Selection::Default)))) {
      return false;
    }
    if (!(track.tpcFractionSharedCls() < cgTrack.cfgCutsMaxTpcNClsSharedRatio.value.get(toI(Selection::Default)))) {
      return false;
    }
    if (!(track.tpcNClsCrossedRows() > cgTrack.cfgCutsMinTpcNCrossedRows.value.get(toI(Selection::Default)))) {
      return false;
    }
    if (!(track.tpcCrossedRowsOverFindableCls() > cgTrack.cfgCutMinTpcNCrossedRowsRatio.value)) {
      return false;
    }
    return true;
  }

  template <ChargeSpecies ChargeSpeciesValue, typename T>
    requires IsValidEnumValue<ChargeSpeciesValue>
  void fillQaRunByTrackByChargeSpecies(const T& track)
  {
    const auto fill{[this](const auto& name, const auto value) -> void {
      hrQaRun.fill(C_CS("pRunIndex") + name + C_CS("_") + C_SV(getName<NameKind::Lower>(ChargeSpeciesValue)), holderEvent.runIndex, value);
    }};
    const auto fillNSigmaPidByDetectorParticleSpecies{
      [this, &fill]<Detector DetectorValue, ParticleSpecies ParticleSpeciesValue>
        requires IsValidEnumValue<DetectorValue, ParticleSpeciesValue>
      () -> void {
        const double nSigmaPid{holderTrack.nSigmaPid[toI(DetectorValue)][toI(ParticleSpeciesValue)]};
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
    requires IsValidEnumValue<ChargeSpeciesValue>
  void fillQaRunByEventByChargeSpecies()
  {
    const auto fill{[this](const auto& name, const auto value) -> void {
      hrQaRun.fill(C_CS("pRunIndex") + name + C_CS("_") + C_SV(getName<NameKind::Lower>(ChargeSpeciesValue)), holderEvent.runIndex, value);
    }};

    fill(C_CS("NGlobalTracks"), holderEvent.nGlobalTracks[toI(ChargeSpeciesValue)]);
    fill(C_CS("NPvContributors"), holderEvent.nPvContributors[toI(ChargeSpeciesValue)]);
    if (holderEvent.nGlobalTracks[toI(ChargeSpeciesValue)] > 0) {
      fill(C_SV(getName(DcaMeasure::Mean)) + C_CS("Dca") + C_SV(getName(DcaAxis::Xy)), holderEvent.measureDca[toI(DcaMeasure::Mean)][toI(DcaAxis::Xy)][toI(ChargeSpeciesValue)]);
      fill(C_SV(getName(DcaMeasure::Sigma)) + C_CS("Dca") + C_SV(getName(DcaAxis::Xy)), holderEvent.measureDca[toI(DcaMeasure::Sigma)][toI(DcaAxis::Xy)][toI(ChargeSpeciesValue)]);
      fill(C_SV(getName(DcaMeasure::Mean)) + C_CS("Dca") + C_SV(getName(DcaAxis::Z)), holderEvent.measureDca[toI(DcaMeasure::Mean)][toI(DcaAxis::Z)][toI(ChargeSpeciesValue)]);
      fill(C_SV(getName(DcaMeasure::Sigma)) + C_CS("Dca") + C_SV(getName(DcaAxis::Z)), holderEvent.measureDca[toI(DcaMeasure::Sigma)][toI(DcaAxis::Z)][toI(ChargeSpeciesValue)]);
    }
    fill(C_CS("NTofBeta"), holderEvent.nTofBeta[toI(ChargeSpeciesValue)]);
  }

  template <ChargeSpecies ChargeSpeciesValue, typename T>
    requires IsValidEnumValue<ChargeSpeciesValue>
  void fillQaTrackByChargeSpecies(const T& track)
  {
    const auto fill{[this](const auto& name, const auto... positionAndWeight) -> void {
      hrQaTrack.fill(C_CS("h") + name + C_CS("_") + C_SV(getName<NameKind::Lower>(ChargeSpeciesValue)), positionAndWeight...);
    }};

    fill(C_CS("ItsNClsChi2NCls"), track.itsNCls(), track.itsChi2NCl());
    fill(C_CS("TpcNClsChi2NCls"), track.tpcNClsFound(), track.tpcChi2NCl());
    fill(C_CS("TpcNClsNClsShared"), track.tpcNClsFound(), track.tpcNClsShared());
    fill(C_CS("TpcNClsFindableNCrossedRows"), track.tpcNClsFindable(), track.tpcNClsCrossedRows());
  }

  template <ChargeSpecies ChargeSpeciesValue>
    requires IsValidEnumValue<ChargeSpeciesValue>
  void fillQaDcaByChargeSpecies()
  {
    const auto fillByDcaAxis{
      [this]<DcaAxis DcaAxisValue>
        requires IsValidEnumValue<DcaAxisValue>
      () -> void {
        hrQaDca.fill(C_CS("hPtDca") + C_SV(getName(DcaAxisValue)) + C_CS("_") + C_SV(getName<NameKind::Lower>(ChargeSpeciesValue)), holderTrack.pt, holderTrack.dca[toI(DcaAxisValue)]);
        hrQaDca.fill(C_CS("pCentralityPtEtaDca") + C_SV(getName(DcaAxisValue)) + C_CS("_") + C_SV(getName<NameKind::Lower>(ChargeSpeciesValue)), holderEvent.centralityCalibration, holderTrack.pt, holderTrack.eta, holderTrack.dca[toI(DcaAxisValue)]);
      }};

    fillByDcaAxis.template operator()<DcaAxis::Xy>();
    fillByDcaAxis.template operator()<DcaAxis::Z>();
  }

  template <ParticleSpeciesAll ParticleSpeciesAllValue, typename T>
    requires IsValidEnumValue<ParticleSpeciesAllValue>
  void fillQaAcceptanceByParticleSpeciesAll(const T& track)
  {
    if (!cgAnalysis.cfgFlagsQaAcceptance.value.get(toI(ParticleSpeciesAllValue))) {
      return;
    }

    const auto fillByChargeSpecies{
      [this]<ChargeSpecies ChargeSpeciesValue>
        requires IsValidEnumValue<ChargeSpeciesValue>
      (const auto& name, const double value) -> void {
        const auto fillByPidStrategy{
          [this, &name, value]<PidStrategy PidStrategyValue>
            requires IsValidEnumValue<PidStrategyValue>
          () -> void {
            if (isPid<getValue<PidStrategyAll>(PidStrategyValue), ParticleSpeciesAllValue>(false)) {
              hrQaAcceptance.fill(C_CS("h") + name + C_CS("Pt_") + C_SV(getName<NameKind::Lower>(PidStrategyValue)) + C_CS("Edge") + C_SV(getName(ParticleSpeciesAllValue)) + C_SV(getName(ChargeSpeciesValue)), value, holderTrack.pt); // NOLINT(clang-analyzer-core.NonNullParamChecker)
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
        fillByChargeSpecies.template operator()<ChargeSpecies::Plus>(C_CS("Rapidity"), track.rapidity(getMass(getValue<ParticleSpecies>(ParticleSpeciesAllValue))));
      } else {
        fillByChargeSpecies.template operator()<ChargeSpecies::Minus>(C_CS("Rapidity"), track.rapidity(getMass(getValue<ParticleSpecies>(ParticleSpeciesAllValue))));
      }
    }
  }

  template <DataMode DataModeValue, ParticleSpeciesAll ParticleSpeciesAllValue>
    requires IsValidEnumValue<DataModeValue, ParticleSpeciesAllValue> && (DataModeValue != DataMode::McMcParticle)
  void fillQaPhiByParticleSpeciesAll()
  {
    if (!cgAnalysis.cfgFlagsQaPhi.value.get(toI(ParticleSpeciesAllValue))) {
      return;
    }

    const auto fillByChargeSpecies{
      [this]<ChargeSpecies ChargeSpeciesValue>
        requires IsValidEnumValue<ChargeSpeciesValue>
      () -> void {
        const auto fillByPidStrategy{
          [this]<PidStrategy PidStrategyValue>
            requires IsValidEnumValue<PidStrategyValue>
          () -> void {
            if constexpr (DataModeValue == DataMode::McTrack) {
              if (isPid<ParticleSpeciesAllValue, ChargeSpeciesValue>() && isPid<getValue<PidStrategyAll>(PidStrategyValue), ParticleSpeciesAllValue>(false)) {
                hrQaPhi.fill(C_CS("hCentralityPtEtaPhi_mc") + C_SV(getName(PidStrategyValue)) + C_SV(getName(ParticleSpeciesAllValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.centralityCalibration, holderTrack.pt, holderTrack.eta, holderTrack.phi);
              }
            } else { // DataModeValue == DataMode::RawTrack
              if (isPid<getValue<PidStrategyAll>(PidStrategyValue), ParticleSpeciesAllValue>(false)) {
                hrQaPhi.fill(C_CS("hCentralityPtEtaPhi_") + C_SV(getName<NameKind::Lower>(PidStrategyValue)) + C_SV(getName(ParticleSpeciesAllValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.centralityCalibration, holderTrack.pt, holderTrack.eta, holderTrack.phi);
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
    requires IsValidEnumValue<DataModeValue, ParticleSpeciesAllValue> && (DataModeValue != DataMode::McMcParticle)
  void fillQaPidByParticleSpeciesAll(const T& track)
  {
    if (!cgAnalysis.cfgFlagsQaPid.value.get(toI(ParticleSpeciesAllValue))) {
      return;
    }

    if constexpr (ParticleSpeciesAllValue == ParticleSpeciesAll::All) {
      if (isPid<getValue<PidStrategyAll>(PidStrategy::Tpc), ParticleSpeciesAll::All>(false)) {
        hrQaPid.fill(C_CS("hCentralityPOverQEtaTpcLnDeDx"), holderEvent.centralityCalibration, track.p() / holderTrack.sign, holderTrack.eta, std::log(track.tpcSignal()));
      }
      if (isPid<getValue<PidStrategyAll>(PidStrategy::TpcTof), ParticleSpeciesAll::All>(false)) {
        hrQaPid.fill(C_CS("hCentralityPOverQEtaTofInverseBeta"), holderEvent.centralityCalibration, track.p() / holderTrack.sign, holderTrack.eta, 1. / track.beta());
      }
    } else {
      constexpr std::int32_t ParticleSpeciesIndex{toI(getValue<ParticleSpecies>(ParticleSpeciesAllValue))};
      const auto fillByChargeSpecies{
        [this]<ChargeSpecies ChargeSpeciesValue>
          requires IsValidEnumValue<ChargeSpeciesValue>
        () -> void {
          if constexpr (DataModeValue == DataMode::McTrack) {
            if (isPid<ParticleSpeciesAllValue, ChargeSpeciesValue>()) {
              hrQaPid.fill(C_CS("hCentralityPtEta") + C_SV(getName(Detector::Tpc)) + C_CS("NSigma") + C_SV(getName(ParticleSpeciesAllValue)) + C_CS("_mc") + C_SV(getName(ParticleSpeciesAllValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.centralityCalibration, holderTrack.pt, holderTrack.eta, holderTrack.nSigmaPid[toI(Detector::Tpc)][ParticleSpeciesIndex]);
              hrQaPid.fill(C_CS("hCentralityPtEta") + C_SV(getName(Detector::Tof)) + C_CS("NSigma") + C_SV(getName(ParticleSpeciesAllValue)) + C_CS("_mc") + C_SV(getName(ParticleSpeciesAllValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.centralityCalibration, holderTrack.pt, holderTrack.eta, holderTrack.nSigmaPid[toI(Detector::Tof)][ParticleSpeciesIndex]);
            }
          } else { // DataModeValue == DataMode::RawTrack
            hrQaPid.fill(C_CS("hCentralityPtEta") + C_SV(getName(Detector::Tpc)) + C_CS("NSigma") + C_SV(getName(ParticleSpeciesAllValue)) + C_CS("_") + C_SV(getName<NameKind::Lower>(ChargeSpeciesValue)), holderEvent.centralityCalibration, holderTrack.pt, holderTrack.eta, holderTrack.nSigmaPid[toI(Detector::Tpc)][ParticleSpeciesIndex]);
            if (isPid<PidStrategyAll::Tof, ParticleSpeciesAllValue>(false)) {
              hrQaPid.fill(C_CS("hCentralityPtEta") + C_SV(getName(Detector::Tpc)) + C_CS("NSigma") + C_SV(getName(ParticleSpeciesAllValue)) + C_CS("_") + C_SV(getName<NameKind::Lower>(Detector::Tof)) + C_SV(getName(ParticleSpeciesAllValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.centralityCalibration, holderTrack.pt, holderTrack.eta, holderTrack.nSigmaPid[toI(Detector::Tpc)][ParticleSpeciesIndex]);
            }
            if (isPid<PidStrategyAll::Tpc, ParticleSpeciesAllValue>(false)) {
              hrQaPid.fill(C_CS("hCentralityPtEta") + C_SV(getName(Detector::Tof)) + C_CS("NSigma") + C_SV(getName(ParticleSpeciesAllValue)) + C_CS("_") + C_SV(getName<NameKind::Lower>(Detector::Tpc)) + C_SV(getName(ParticleSpeciesAllValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.centralityCalibration, holderTrack.pt, holderTrack.eta, holderTrack.nSigmaPid[toI(Detector::Tof)][ParticleSpeciesIndex]);
            }
            hrQaPid.fill(C_CS("hCentralityPtEta") + C_SV(getName(PidStrategy::TpcTof)) + C_CS("NSigma") + C_SV(getName(ParticleSpeciesAllValue)) + C_CS("_") + C_SV(getName<NameKind::Lower>(ChargeSpeciesValue)), holderEvent.centralityCalibration, holderTrack.pt, holderTrack.eta, holderTrack.getNSigmaPidCombined(ParticleSpeciesIndex));
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
    requires IsValidEnumValue<DataModeValue, ParticleSpeciesValue>
  void fillCalculationYieldByParticleSpecies()
  {
    if (!cgAnalysis.cfgFlagsCalculationYield.value.get(toI(ParticleSpeciesValue))) {
      return;
    }

    const std::int32_t chargeSign{[](const std::int32_t chargeMcParticle, const std::int32_t signTrack) constexpr -> std::int32_t {
      if constexpr (DataModeValue == DataMode::McMcParticle) {
        return chargeMcParticle;
      } else {
        return signTrack;
      }
    }(holderMcParticle.charge, holderTrack.sign)};

    const auto fillByChargeSpecies{
      [this]<ChargeSpecies ChargeSpeciesValue>
        requires IsValidEnumValue<ChargeSpeciesValue>
      () -> void {
        if constexpr (DataModeValue == DataMode::McMcParticle) {
          if (isPid<getValue<ParticleSpeciesAll>(ParticleSpeciesValue), ChargeSpeciesValue>()) {
            hrCalculationYield.fill(C_CS("hVzCentralityPtMcEtaMc_mc") + C_SV(getName(ParticleSpeciesValue)) + C_SV(getName(ChargeSpeciesValue)), cgEvent.cfgFlagMcCollisionVz.value ? holderMcEvent.vz : holderEvent.vz, holderEvent.centrality, holderMcParticle.pt, holderMcParticle.eta);
          }
        } else {
          const auto fillByPidStrategy{
            [this]<PidStrategy PidStrategyValue>
              requires IsValidEnumValue<PidStrategyValue>
            () -> void {
              if constexpr (DataModeValue == DataMode::McTrack) {
                if (isPid<getValue<ParticleSpeciesAll>(ParticleSpeciesValue), ChargeSpeciesValue>() && isPid<getValue<PidStrategyAll>(PidStrategyValue), getValue<ParticleSpeciesAll>(ParticleSpeciesValue)>(cgTrack.cfgFlagRejectionOthers.value)) {
                  if (cgTrack.cfgFlagMcParticleMomentum.value) {
                    hrCalculationYield.fill(C_CS("hVzCentralityPtMcEtaMc_mc") + C_SV(getName(PidStrategyValue)) + C_SV(getName(ParticleSpeciesValue)) + C_SV(getName(ChargeSpeciesValue)), cgEvent.cfgFlagMcCollisionVz.value ? holderMcEvent.vz : holderEvent.vz, holderEvent.centrality, holderMcParticle.pt, holderMcParticle.eta);
                  } else {
                    hrCalculationYield.fill(C_CS("hVzCentralityPtEta_mc") + C_SV(getName(PidStrategyValue)) + C_SV(getName(ParticleSpeciesValue)) + C_SV(getName(ChargeSpeciesValue)), cgEvent.cfgFlagMcCollisionVz.value ? holderMcEvent.vz : holderEvent.vz, holderEvent.centrality, holderTrack.pt, holderTrack.eta);
                  }
                }
              } else { // DataModeValue == DataMode::RawTrack
                if (isPid<getValue<PidStrategyAll>(PidStrategyValue), getValue<ParticleSpeciesAll>(ParticleSpeciesValue)>(cgTrack.cfgFlagRejectionOthers.value)) {
                  hrCalculationYield.fill(C_CS("hVzCentralityPtEta_") + C_SV(getName<NameKind::Lower>(PidStrategyValue)) + C_SV(getName(ParticleSpeciesValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.vz, holderEvent.centrality, holderTrack.pt, holderTrack.eta);
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
    requires IsValidEnumValue<ParticleSpeciesValue>
  void fillCalculationPurityByParticleSpecies()
  {
    if (!cgAnalysis.cfgFlagsCalculationPurity.value.get(toI(ParticleSpeciesValue))) {
      return;
    }

    const auto fillByChargeSpecies{
      [this]<ChargeSpecies ChargeSpeciesValue>
        requires IsValidEnumValue<ChargeSpeciesValue>
      () -> void {
        const auto fillByPidStrategy{
          [this]<PidStrategy PidStrategyValue>
            requires IsValidEnumValue<PidStrategyValue>
          () -> void {
            if (isPid<getValue<PidStrategyAll>(PidStrategyValue), getValue<ParticleSpeciesAll>(ParticleSpeciesValue)>(cgTrack.cfgFlagRejectionOthers.value)) {
              hrCalculationPurity.fill(C_CS("pCentralityPtEtaPurity") + C_SV(getName(PidStrategyValue)) + C_SV(getName(ParticleSpeciesValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.centrality, holderTrack.pt, holderTrack.eta, isPid<getValue<ParticleSpeciesAll>(ParticleSpeciesValue), ChargeSpeciesValue>() ? 1. : 0.);
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
    requires IsValidEnumValue<ParticleSpeciesValue>
  void fillCalculationFractionPrimaryByParticleSpecies(const MP& mcParticle)
  {
    if (!cgAnalysis.cfgFlagsCalculationFractionPrimary.value.get(toI(ParticleSpeciesValue))) {
      return;
    }

    const auto fillByChargeSpecies{
      [this, &mcParticle]<ChargeSpecies ChargeSpeciesValue>
        requires IsValidEnumValue<ChargeSpeciesValue>
      () -> void {
        const auto fillByPidStrategy{
          [this, &mcParticle]<PidStrategy PidStrategyValue>
            requires IsValidEnumValue<PidStrategyValue>
          () -> void {
            if (isPid<getValue<ParticleSpeciesAll>(ParticleSpeciesValue), ChargeSpeciesValue>() && isPid<getValue<PidStrategyAll>(PidStrategyValue), getValue<ParticleSpeciesAll>(ParticleSpeciesValue)>(cgTrack.cfgFlagRejectionOthers.value)) {
              hrCalculationFractionPrimary.fill(C_CS("pCentralityPtEtaFractionPrimary") + C_SV(getName(PidStrategyValue)) + C_SV(getName(ParticleSpeciesValue)) + C_SV(getName(ChargeSpeciesValue)), holderEvent.centrality, holderTrack.pt, holderTrack.eta, mcParticle.isPhysicalPrimary() ? 1. : 0.);
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
      if (static_cast<bool>(cgAnalysis.cfgFlagsCalculationFluctuation.value.get(iParticleNumber))) {
        for (std::int32_t const& iChargeNumber : std::views::iota(0, NEs<ChargeNumber>)) {
          if (doProcessMc.value) {
            fluctuationCalculatorsTrackMcParticle[iParticleNumber][iChargeNumber]->clear();
          }
          fluctuationCalculatorsTrackTrack[iParticleNumber][iChargeNumber]->clear();
        }
      }
    }
  }

  template <DataMode DataModeValue, ParticleNumber ParticleNumberValue>
    requires IsValidEnumValue<DataModeValue, ParticleNumberValue>
  void calculateFluctuationByParticleNumber()
  {
    if (!cgAnalysis.cfgFlagsCalculationFluctuation.value.get(toI(ParticleNumberValue))) {
      return;
    }

    const std::int32_t chargeSign{[](const std::int32_t chargeMcParticle, const std::int32_t signTrack) constexpr -> std::int32_t {
      if constexpr (DataModeValue == DataMode::McMcParticle) {
        return chargeMcParticle;
      } else {
        return signTrack;
      }
    }(holderMcParticle.charge, holderTrack.sign)};

    const bool doUseMcParticleMomentum{[](const bool flagMcParticleMomentum) constexpr -> bool {
      if constexpr (DataModeValue == DataMode::McMcParticle) {
        return true;
      } else if constexpr (DataModeValue == DataMode::McTrack) {
        return flagMcParticleMomentum;
      } else {
        return false;
      }
    }(cgTrack.cfgFlagMcParticleMomentum.value)};
    const ChargeSpecies chargeSpecies{chargeSign > 0 ? ChargeSpecies::Plus : ChargeSpecies::Minus};
    if (isGoodMomentum<getValue<ParticleSpeciesAll>(ParticleNumberValue)>(doUseMcParticleMomentum)) {
      if constexpr (DataModeValue == DataMode::McMcParticle) {
        ++holderDerivedData.nMcParticles[toI(chargeSpecies)];
      } else {
        ++holderDerivedData.nTracks[toI(chargeSpecies)];
      }
    }

    const auto calculateByParticleSpecies{
      [this, chargeSign, doUseMcParticleMomentum]<ParticleSpecies ParticleSpeciesValue>
        requires IsValidEnumValue<ParticleSpeciesValue> && (getValue<ParticleSpeciesAll>(ParticleNumberValue) == ParticleSpeciesAll::All || getValue<ParticleSpeciesAll>(ParticleNumberValue) == getValue<ParticleSpeciesAll>(ParticleSpeciesValue))
      () -> bool {
        if (!isGoodMomentum<getValue<ParticleSpeciesAll>(ParticleSpeciesValue)>(doUseMcParticleMomentum) || (DataModeValue != DataMode::RawTrack && (chargeSign > 0 ? !isPid<getValue<ParticleSpeciesAll>(ParticleSpeciesValue), ChargeSpecies::Plus>() : !isPid<getValue<ParticleSpeciesAll>(ParticleSpeciesValue), ChargeSpecies::Minus>()))) {
          return false;
        }

        const bool doUseTofPid{[](const bool doUseMcParticleMomentumValue, const double ptMcParticle, const double ptTrack, const double thresholdPtTofPid) constexpr -> bool {
          if constexpr (DataModeValue == DataMode::McMcParticle) {
            return ptMcParticle >= thresholdPtTofPid;
          } else if constexpr (DataModeValue == DataMode::McTrack) {
            return (doUseMcParticleMomentumValue ? ptMcParticle : ptTrack) >= thresholdPtTofPid;
          } else {
            return ptTrack >= thresholdPtTofPid;
          }
        }(doUseMcParticleMomentum, holderMcParticle.pt, holderTrack.pt, cgTrack.cfgThresholdsPtTofPid.value.get(toI(ParticleSpeciesValue)))};
        if constexpr (DataModeValue != DataMode::McMcParticle) {
          if (!(doUseTofPid ? isPid<getValue<PidStrategyAll>(PidStrategy::TpcTof), getValue<ParticleSpeciesAll>(ParticleSpeciesValue)>(cgTrack.cfgFlagRejectionOthers.value) : isPid<getValue<PidStrategyAll>(PidStrategy::Tpc), getValue<ParticleSpeciesAll>(ParticleSpeciesValue)>(cgTrack.cfgFlagRejectionOthers.value))) {
            return false;
          }
        }

        const auto calculateByChargeSpecies{
          [this, chargeSign, doUseMcParticleMomentum, doUseTofPid]<ChargeSpecies ChargeSpeciesValue>
            requires IsValidEnumValue<ChargeSpeciesValue>
          () -> void {
            const double efficiency{doUseTofPid ? getEfficiency<PidStrategy::TpcTof, ParticleSpeciesValue, ChargeSpeciesValue>(doUseMcParticleMomentum) : getEfficiency<PidStrategy::Tpc, ParticleSpeciesValue, ChargeSpeciesValue>(doUseMcParticleMomentum)};
            const auto fill{
              [this, efficiency]() -> void {
                std::array<std::array<std::unique_ptr<FluctuationCalculatorTrack>, NEs<ChargeNumber>>, NEs<ParticleNumber>>& fluctuationCalculatorsTrack{DataModeValue == DataMode::McMcParticle ? fluctuationCalculatorsTrackMcParticle : fluctuationCalculatorsTrackTrack};
                if constexpr (ChargeSpeciesValue == ChargeSpecies::Plus) {
                  fluctuationCalculatorsTrack[toI(ParticleNumberValue)][toI(ChargeNumber::Plus)]->fill(1., efficiency);
                  fluctuationCalculatorsTrack[toI(ParticleNumberValue)][toI(ChargeNumber::Net)]->fill(1., efficiency);
                } else {
                  fluctuationCalculatorsTrack[toI(ParticleNumberValue)][toI(ChargeNumber::Minus)]->fill(1., efficiency);
                  fluctuationCalculatorsTrack[toI(ParticleNumberValue)][toI(ChargeNumber::Net)]->fill(-1., efficiency);
                }
                fluctuationCalculatorsTrack[toI(ParticleNumberValue)][toI(ChargeNumber::Total)]->fill(1., efficiency);
              }};
            if constexpr (DataModeValue == DataMode::McMcParticle) {
              ++holderMcEvent.numbers[toI(ParticleNumberValue)][toI(ChargeSpeciesValue)];
              if (gRandom->Rndm() < efficiency) {
                ++holderMcEvent.numbersEff[toI(ParticleNumberValue)][toI(ChargeSpeciesValue)];
                fill();
              }
              holderDerivedData.signedEfficienciesMcParticle.push_back(HolderDerivedData::convert<aod::tiny_mc_particle::SignedEfficiency::type>(std::copysign(std::numeric_limits<aod::tiny_mc_particle::SignedEfficiency::type>::max(), chargeSign) * efficiency));
            } else {
              ++holderEvent.numbers[toI(ParticleNumberValue)][toI(ChargeSpeciesValue)];
              fill();
              holderDerivedData.signedEfficienciesTrack.push_back(HolderDerivedData::convert<aod::tiny_track::SignedEfficiency::type>(std::copysign(std::numeric_limits<aod::tiny_track::SignedEfficiency::type>::max(), chargeSign) * efficiency));
            }
          }};

        if (chargeSign > 0) {
          calculateByChargeSpecies.template operator()<ChargeSpecies::Plus>();
        } else {
          calculateByChargeSpecies.template operator()<ChargeSpecies::Minus>();
        }
        return true;
      }};

    if constexpr (getValue<ParticleSpeciesAll>(ParticleNumberValue) == ParticleSpeciesAll::All) {
      if (!calculateByParticleSpecies.template operator()<ParticleSpecies::Pion>() && !calculateByParticleSpecies.template operator()<ParticleSpecies::Kaon>()) {
        calculateByParticleSpecies.template operator()<ParticleSpecies::Proton>();
      }
    } else {
      calculateByParticleSpecies.template operator()<getValue<ParticleSpecies>(getValue<ParticleSpeciesAll>(ParticleNumberValue))>();
    }
  }

  template <ParticleNumber ParticleNumberValue>
    requires IsValidEnumValue<ParticleNumberValue>
  void fillCalculationFluctuationByParticleNumber()
  {
    if (!cgAnalysis.cfgFlagsCalculationFluctuation.value.get(toI(ParticleNumberValue))) {
      return;
    }

    if (doProcessMc.value) {
      hrCalculationFluctuation.fill(C_CS("hCentralityN") + C_SV(getName(ParticleNumberValue)) + C_SV(getName(ChargeSpecies::Plus)) + C_CS("N") + C_SV(getName(ParticleNumberValue)) + C_SV(getName(ChargeSpecies::Minus)) + C_CS("_mc"), holderEvent.centrality, holderMcEvent.numbers[toI(ParticleNumberValue)][toI(ChargeSpecies::Plus)], holderMcEvent.numbers[toI(ParticleNumberValue)][toI(ChargeSpecies::Minus)]);
      hrCalculationFluctuation.fill(C_CS("hCentralityN") + C_SV(getName(ParticleNumberValue)) + C_SV(getName(ChargeSpecies::Plus)) + C_CS("N") + C_SV(getName(ParticleNumberValue)) + C_SV(getName(ChargeSpecies::Minus)) + C_CS("_mcEff"), holderEvent.centrality, holderMcEvent.numbersEff[toI(ParticleNumberValue)][toI(ChargeSpecies::Plus)], holderMcEvent.numbersEff[toI(ParticleNumberValue)][toI(ChargeSpecies::Minus)]);
    }
    hrCalculationFluctuation.fill(C_CS("hCentralityN") + C_SV(getName(ParticleNumberValue)) + C_SV(getName(ChargeSpecies::Plus)) + C_CS("N") + C_SV(getName(ParticleNumberValue)) + C_SV(getName(ChargeSpecies::Minus)), holderEvent.centrality, holderEvent.numbers[toI(ParticleNumberValue)][toI(ChargeSpecies::Plus)], holderEvent.numbers[toI(ParticleNumberValue)][toI(ChargeSpecies::Minus)]);

    const auto fillByChargeNumber{
      [this]<ChargeNumber ChargeNumberValue>
        requires IsValidEnumValue<ChargeNumberValue>
      () -> void {
        if (doProcessMc.value) {
          const std::array<double, fluctuation_calculator_base::NOrderKeys> products{fluctuationCalculatorsTrackMcParticle[toI(ParticleNumberValue)][toI(ChargeNumberValue)]->getProducts()};
          for (std::int32_t const& iOrderKey : std::views::iota(0, fluctuation_calculator_base::NOrderKeys)) {
            hrCalculationFluctuation.fill(C_CS("hFluctuationCalculator") + C_SV(getName(ParticleNumberValue)) + C_SV(getName(ChargeNumberValue)) + C_CS("_mc"), holderEvent.centrality, holderEvent.subgroupIndex, iOrderKey, products[iOrderKey]);
          }
        }

        const std::array<double, fluctuation_calculator_base::NOrderKeys> products{fluctuationCalculatorsTrackTrack[toI(ParticleNumberValue)][toI(ChargeNumberValue)]->getProducts()};
        for (std::int32_t const& iOrderKey : std::views::iota(0, fluctuation_calculator_base::NOrderKeys)) {
          hrCalculationFluctuation.fill(C_CS("hFluctuationCalculator") + C_SV(getName(ParticleNumberValue)) + C_SV(getName(ChargeNumberValue)), holderEvent.centrality, holderEvent.subgroupIndex, iOrderKey, products[iOrderKey]);
        }
      }};

    fillByChargeNumber.template operator()<ChargeNumber::Plus>();
    fillByChargeNumber.template operator()<ChargeNumber::Minus>();
    fillByChargeNumber.template operator()<ChargeNumber::Total>();
    fillByChargeNumber.template operator()<ChargeNumber::Net>();
  }

  template <bool DoRecalibrate, typename T>
  bool setTrack(const T& track)
  {
    if (std::abs(track.sign()) != 1) {
      return false;
    }

    holderTrack.clear();
    holderTrack.dca[toI(DcaAxis::Xy)] = track.dcaXY();
    holderTrack.dca[toI(DcaAxis::Z)] = track.dcaZ();
    holderTrack.sign = track.sign();
    holderTrack.pt = track.pt();
    holderTrack.eta = track.eta();
    holderTrack.phi = track.phi();
    holderTrack.hasPid[toI(Detector::Tpc)] = (track.hasTPC() && track.tpcSignal() > 0.);
    holderTrack.hasPid[toI(Detector::Tof)] = (track.hasTOF() && track.beta() > 0.);
    if (holderTrack.hasPid[toI(Detector::Tpc)]) {
      holderTrack.nSigmaPid[toI(Detector::Tpc)][toI(ParticleSpecies::Pion)] = HolderTrack::truncateNSigmaPid(track.tpcNSigmaPi(), getShiftNSigmaPid<DoRecalibrate, Detector::Tpc, ParticleSpecies::Pion>());
      holderTrack.nSigmaPid[toI(Detector::Tpc)][toI(ParticleSpecies::Kaon)] = HolderTrack::truncateNSigmaPid(track.tpcNSigmaKa(), getShiftNSigmaPid<DoRecalibrate, Detector::Tpc, ParticleSpecies::Kaon>());
      holderTrack.nSigmaPid[toI(Detector::Tpc)][toI(ParticleSpecies::Proton)] = HolderTrack::truncateNSigmaPid(track.tpcNSigmaPr(), getShiftNSigmaPid<DoRecalibrate, Detector::Tpc, ParticleSpecies::Proton>());
    }
    if (holderTrack.hasPid[toI(Detector::Tof)]) {
      holderTrack.nSigmaPid[toI(Detector::Tof)][toI(ParticleSpecies::Pion)] = HolderTrack::truncateNSigmaPid(track.tofNSigmaPi(), getShiftNSigmaPid<DoRecalibrate, Detector::Tof, ParticleSpecies::Pion>());
      holderTrack.nSigmaPid[toI(Detector::Tof)][toI(ParticleSpecies::Kaon)] = HolderTrack::truncateNSigmaPid(track.tofNSigmaKa(), getShiftNSigmaPid<DoRecalibrate, Detector::Tof, ParticleSpecies::Kaon>());
      holderTrack.nSigmaPid[toI(Detector::Tof)][toI(ParticleSpecies::Proton)] = HolderTrack::truncateNSigmaPid(track.tofNSigmaPr(), getShiftNSigmaPid<DoRecalibrate, Detector::Tof, ParticleSpecies::Proton>());
    }

    return true;
  }

  template <typename MP>
  bool initMcParticle(const MP& mcParticle)
  {
    holderMcParticle.clear();
    holderMcParticle.pdgCode = mcParticle.pdgCode();
    const TParticlePDG* const particlePdg{pdg->GetParticle(mcParticle.pdgCode())};
    if (particlePdg) {
      holderMcParticle.charge = static_cast<std::int32_t>(std::rint(particlePdg->Charge()));
    } else {
      switch (std::abs(holderMcParticle.pdgCode) / 100000000) {
        case 10:
          holderMcParticle.charge = holderMcParticle.pdgCode / 10000 % 1000;
          break;
        default:
          break;
      }
    }

    if (holderMcParticle.charge == 0) {
      return false;
    }

    holderMcParticle.pt = mcParticle.pt();
    holderMcParticle.eta = mcParticle.eta();

    return true;
  }

  template <bool DoInitEvent, typename T>
  bool initTrack(const T& track)
  {
    if (!setTrack<!DoInitEvent>(track)) {
      return false;
    }

    if constexpr (DoInitEvent) {
      if (track.isPrimaryTrack()) {
        const std::int32_t chargeSpeciesIndex{toI(holderTrack.sign > 0 ? ChargeSpecies::Plus : ChargeSpecies::Minus)};
        ++holderEvent.nGlobalTracks[chargeSpeciesIndex];
        if (track.isPVContributor()) {
          ++holderEvent.nPvContributors[chargeSpeciesIndex];
        }
        for (std::int32_t const& iDcaAxis : std::views::iota(0, NEs<DcaAxis>)) {
          holderEvent.measureDca[toI(DcaMeasure::Mean)][iDcaAxis][chargeSpeciesIndex] += holderTrack.dca[iDcaAxis];
          holderEvent.measureDca[toI(DcaMeasure::Sigma)][iDcaAxis][chargeSpeciesIndex] += std::pow(holderTrack.dca[iDcaAxis], 2.);
        }
        if (holderTrack.hasPid[toI(Detector::Tof)]) {
          ++holderEvent.nTofBeta[chargeSpeciesIndex];
        }
      }

      if (cgAnalysis.cfgFlagQaRun.value && track.isPrimaryTrack()) {
        if (holderTrack.sign > 0) {
          fillQaRunByTrackByChargeSpecies<ChargeSpecies::Plus>(track);
        } else {
          fillQaRunByTrackByChargeSpecies<ChargeSpecies::Minus>(track);
        }
      }
    } else {
      if (cgAnalysis.cfgFlagQaTrack.value && track.isPrimaryTrack()) {
        if (holderTrack.sign > 0) {
          fillQaTrackByChargeSpecies<ChargeSpecies::Plus>(track);
        } else {
          fillQaTrackByChargeSpecies<ChargeSpecies::Minus>(track);
        }
      }

      if (!isGoodTrack(track)) {
        return false;
      }

      if (cgAnalysis.cfgFlagQaDca.value) {
        if (holderTrack.sign > 0) {
          fillQaDcaByChargeSpecies<ChargeSpecies::Plus>();
        } else {
          fillQaDcaByChargeSpecies<ChargeSpecies::Minus>();
        }
      }

      if (!isGoodDca()) {
        return false;
      }

      if (doQaAcceptance) {
        const double vz{doProcessMc.value && cgEvent.cfgFlagMcCollisionVz.value ? holderMcEvent.vz : holderEvent.vz};
        if (holderTrack.eta * vz > 0. && std::abs(vz) > (doProcessMc.value && cgEvent.cfgFlagMcCollisionVz.value ? cgEvent.cfgCutMaxAbsVzMc.value : cgEvent.cfgCutMaxAbsVz.value) - 1.) {
          fillQaAcceptanceByParticleSpeciesAll<ParticleSpeciesAll::All>(track);
          fillQaAcceptanceByParticleSpeciesAll<ParticleSpeciesAll::Pion>(track);
          fillQaAcceptanceByParticleSpeciesAll<ParticleSpeciesAll::Kaon>(track);
          fillQaAcceptanceByParticleSpeciesAll<ParticleSpeciesAll::Proton>(track);
        }
      }
    }

    return true;
  }

  template <ParticleNumber ParticleNumberValue, typename MPs, typename Ts>
    requires IsValidEnumValue<ParticleNumberValue>
  void fillMiniTableMc(const MPs& mcParticles, const Ts& tracks, const bool isGoodNPvContributors)
  {
    if (!static_cast<bool>(cgAnalysis.cfgFlagsStorageMiniTable.value.get(toI(ParticleNumberValue)))) {
      return;
    }

    const std::optional<aod::mini_collision::Code::type> codeCollision{mini_collision_codec::encode(cgEvent.cfgFlagMcCollisionVz.value ? holderMcEvent.vz : holderEvent.vz, holderEvent.centrality, isGoodNPvContributors)};
    if (!codeCollision.has_value()) {
      return;
    }

    miniCollision(*codeCollision);

    for (const auto& mcParticle : mcParticles) {
      if (!initMcParticle(mcParticle)) {
        continue;
      }

      if (mcParticle.isPhysicalPrimary()) {
        const auto fillByParticleSpecies{
          [this]<ParticleSpecies ParticleSpeciesValue>
            requires IsValidEnumValue<ParticleSpeciesValue> && (getValue<ParticleSpeciesAll>(ParticleNumberValue) == ParticleSpeciesAll::All || getValue<ParticleSpeciesAll>(ParticleNumberValue) == getValue<ParticleSpeciesAll>(ParticleSpeciesValue))
          () -> bool {
            if (!isGoodMomentum<getValue<ParticleSpeciesAll>(ParticleSpeciesValue)>(true) || (holderMcParticle.charge > 0 ? !isPid<getValue<ParticleSpeciesAll>(ParticleSpeciesValue), ChargeSpecies::Plus>() : !isPid<getValue<ParticleSpeciesAll>(ParticleSpeciesValue), ChargeSpecies::Minus>())) {
              return false;
            }

            const std::optional<aod::mini_mc_particle::Code::type> codeMcParticle{mini_mc_particle_codec::encode<ParticleSpeciesValue>(holderMcParticle.pt, holderMcParticle.eta, holderMcParticle.charge)};
            if (codeMcParticle.has_value()) {
              miniMcParticle(miniCollision.lastIndex(), *codeMcParticle);
            }
            return true;
          }};

        if constexpr (ParticleNumberValue == ParticleNumber::Charge) {
          if (!fillByParticleSpecies.template operator()<ParticleSpecies::Pion>() && !fillByParticleSpecies.template operator()<ParticleSpecies::Kaon>()) {
            fillByParticleSpecies.template operator()<ParticleSpecies::Proton>();
          }
        } else {
          fillByParticleSpecies.template operator()<getValue<ParticleSpecies>(getValue<ParticleSpeciesAll>(ParticleNumberValue))>();
        }
      }

      if (!cgTrack.cfgFlagMcParticlePhysicalPrimary.value || mcParticle.isPhysicalPrimary()) {
        const auto& tracksMatched{tracks.sliceBy(presliceTracksPerMcParticle, mcParticle.globalIndex())};
        for (const auto& track : tracksMatched) {
          if (!setTrack<true>(track) || !(cgTrack.cfgCutMinTpcChi2NCls.value < track.tpcChi2NCl()) || !(track.tpcCrossedRowsOverFindableCls() > cgTrack.cfgCutMinTpcNCrossedRowsRatio.value)) {
            continue;
          }

          const auto fillByParticleSpeciesAll{
            [this, &track]<ParticleSpeciesAll ParticleSpeciesAllValue>
              requires IsValidEnumValue<ParticleSpeciesAllValue> && (ParticleSpeciesAllValue != ParticleSpeciesAll::All)
            () -> bool {
              const bool doUseMcParticleMomentum{cgTrack.cfgFlagMcParticleMomentum.value};
              if (!isGoodMomentum<ParticleSpeciesAllValue>(doUseMcParticleMomentum) || (holderTrack.sign > 0 ? !isPid<ParticleSpeciesAllValue, ChargeSpecies::Plus>() : !isPid<ParticleSpeciesAllValue, ChargeSpecies::Minus>())) {
                return false;
              }

              constexpr std::int32_t ParticleSpeciesIndex{toI(getValue<ParticleSpecies>(ParticleSpeciesAllValue))};
              const double pt{doUseMcParticleMomentum ? holderMcParticle.pt : holderTrack.pt};
              const bool doUseTofPid{pt >= cgTrack.cfgThresholdsPtTofPid.value.get(ParticleSpeciesIndex)};
              if (!(doUseTofPid ? isPid<getValue<PidStrategyAll>(PidStrategy::TpcTof), ParticleSpeciesAllValue, Selection::Loose>(cgTrack.cfgFlagRejectionOthers.value) : isPid<getValue<PidStrategyAll>(PidStrategy::Tpc), ParticleSpeciesAllValue, Selection::Loose>(cgTrack.cfgFlagRejectionOthers.value))) {
                return false;
              }

              const std::optional<aod::mini_track::Code::type> codeTrack{mini_track_codec::encode<ParticleSpeciesAllValue>(*configSelection, track.isPVContributor(), track.itsNCls(), track.itsChi2NCl(), track.tpcNClsFound(), track.tpcChi2NCl(), track.tpcNClsCrossedRows(), track.tpcFractionSharedCls(), {getAbsNSigmaDca<DcaAxis::Xy>(), getAbsNSigmaDca<DcaAxis::Z>()}, pt, doUseMcParticleMomentum ? holderMcParticle.eta : holderTrack.eta, holderTrack.sign, doUseTofPid ? std::abs(holderTrack.getNSigmaPidCombined(ParticleSpeciesIndex)) : std::abs(holderTrack.nSigmaPid[toI(Detector::Tpc)][ParticleSpeciesIndex]))};
              if (codeTrack.has_value()) {
                miniTrack(miniCollision.lastIndex(), *codeTrack);
              }
              return true;
            }};

          if constexpr (ParticleNumberValue == ParticleNumber::Charge) {
            if (!fillByParticleSpeciesAll.template operator()<ParticleSpeciesAll::Pion>() && !fillByParticleSpeciesAll.template operator()<ParticleSpeciesAll::Kaon>()) {
              fillByParticleSpeciesAll.template operator()<ParticleSpeciesAll::Proton>();
            }
          } else {
            fillByParticleSpeciesAll.template operator()<getValue<ParticleSpeciesAll>(ParticleNumberValue)>();
          }
        }
      }
    }
  }

  template <ParticleNumber ParticleNumberValue, typename T>
    requires IsValidEnumValue<ParticleNumberValue>
  void fillMiniTableRaw(const T& tracks, const bool isGoodNPvContributors)
  {
    if (!static_cast<bool>(cgAnalysis.cfgFlagsStorageMiniTable.value.get(toI(ParticleNumberValue)))) {
      return;
    }

    const std::optional<aod::mini_collision::Code::type> codeCollision{mini_collision_codec::encode(holderEvent.vz, holderEvent.centrality, isGoodNPvContributors)};
    if (!codeCollision.has_value()) {
      return;
    }

    miniCollision(*codeCollision);

    for (const auto& track : tracks) {
      if (!setTrack<true>(track) || !(cgTrack.cfgCutMinTpcChi2NCls.value < track.tpcChi2NCl()) || !(track.tpcCrossedRowsOverFindableCls() > cgTrack.cfgCutMinTpcNCrossedRowsRatio.value)) {
        continue;
      }

      const auto fill{
        [this, &track]<ParticleSpeciesAll ParticleSpeciesAllValue>
          requires IsValidEnumValue<ParticleSpeciesAllValue>
        (const double absNSigmaPid) -> void {
          const std::optional<aod::mini_track::Code::type> codeTrack{mini_track_codec::encode<ParticleSpeciesAllValue>(*configSelection, track.isPVContributor(), track.itsNCls(), track.itsChi2NCl(), track.tpcNClsFound(), track.tpcChi2NCl(), track.tpcNClsCrossedRows(), track.tpcFractionSharedCls(), {getAbsNSigmaDca<DcaAxis::Xy>(), getAbsNSigmaDca<DcaAxis::Z>()}, holderTrack.pt, holderTrack.eta, holderTrack.sign, absNSigmaPid)};
          if (codeTrack.has_value()) {
            miniTrack(miniCollision.lastIndex(), *codeTrack);
          }
        }};
      const auto fillByParticleSpeciesAll{
        [this, &fill]<ParticleSpeciesAll ParticleSpeciesAllValue>
          requires IsValidEnumValue<ParticleSpeciesAllValue>
        () -> bool {
          if constexpr (ParticleSpeciesAllValue == ParticleSpeciesAll::All) {
            if (!isGoodMomentum<getValue<ParticleSpeciesAll>(ParticleNumberValue)>(false)) {
              return false;
            }

            fill.template operator()<ParticleSpeciesAll::All>(0.);
          } else {
            if (!isGoodMomentum<ParticleSpeciesAllValue>(false)) {
              return false;
            }

            constexpr std::int32_t ParticleSpeciesIndex{toI(getValue<ParticleSpecies>(ParticleSpeciesAllValue))};
            const bool doUseTofPid{holderTrack.pt >= cgTrack.cfgThresholdsPtTofPid.value.get(ParticleSpeciesIndex)};
            if (!(doUseTofPid ? isPid<getValue<PidStrategyAll>(PidStrategy::TpcTof), ParticleSpeciesAllValue, Selection::Loose>(cgTrack.cfgFlagRejectionOthers.value) : isPid<getValue<PidStrategyAll>(PidStrategy::Tpc), ParticleSpeciesAllValue, Selection::Loose>(cgTrack.cfgFlagRejectionOthers.value))) {
              return false;
            }

            fill.template operator()<ParticleSpeciesAllValue>(doUseTofPid ? std::abs(holderTrack.getNSigmaPidCombined(ParticleSpeciesIndex)) : std::abs(holderTrack.nSigmaPid[toI(Detector::Tpc)][ParticleSpeciesIndex]));
          }
          return true;
        }};

      if constexpr (ParticleNumberValue == ParticleNumber::Charge) {
        if (!fillByParticleSpeciesAll.template operator()<ParticleSpeciesAll::Pion>() && !fillByParticleSpeciesAll.template operator()<ParticleSpeciesAll::Kaon>() && !fillByParticleSpeciesAll.template operator()<ParticleSpeciesAll::Proton>()) {
          fillByParticleSpeciesAll.template operator()<ParticleSpeciesAll::All>();
        }
      } else {
        if (!fillByParticleSpeciesAll.template operator()<getValue<ParticleSpeciesAll>(ParticleNumberValue)>()) {
          fillByParticleSpeciesAll.template operator()<ParticleSpeciesAll::All>();
        }
      }
    }
  }

  template <typename MC>
  McEventSelection initMcEvent(const MC& mcCollision)
  {
    holderMcEvent.clear();
    holderMcEvent.vz = mcCollision.posZ();

    const auto fillMcEventSelection{
      [this](const McEventSelection selection) -> McEventSelection {
        hrCounter.fill(C_CS("hNEventsMc"), toI(selection));
        return selection;
      }};

    fillMcEventSelection(McEventSelection::All);

    if (cgEvent.cfgFlagInelMc.value && !mcCollision.isInelGt0()) {
      return fillMcEventSelection(McEventSelection::Inel);
    }

    if (!(std::abs(holderMcEvent.vz) < cgEvent.cfgCutMaxAbsVzMc.value)) {
      return fillMcEventSelection(McEventSelection::Vz);
    }

    if (cgAnalysis.cfgFlagQaMc.value) {
      hrQaMc.fill(C_CS("hNCollisionsPerMcCollision"), mcCollision.numRecoCollision());
    }

    if (cgEvent.cfgFlagSingleCollisionMc.value && mcCollision.numRecoCollision() != 1) {
      return fillMcEventSelection(McEventSelection::NRecoCollisions);
    }

    return fillMcEventSelection(McEventSelection::Good);
  }

  template <typename C, typename Ts>
  EventSelection initEvent(const C& collision, const Ts& tracks)
  {
    holderEvent.clear();
    holderEvent.vz = collision.posZ();
    switch (cgEvent.cfgIndexDefinitionCentrality.value) {
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
    holderEvent.centralityCalibration = (cgEvent.cfgFlagDefinitionCentralitySameQa.value ? holderEvent.centrality : collision.centNTPV());

    const auto fillEventSelection{
      [this](const EventSelection selection, const auto... selectionsAdditional) -> EventSelection
        requires((std::is_arithmetic_v<std::remove_cvref_t<decltype(selectionsAdditional)>> && ...))
      {
        hrCounter.fill(C_CS("hNEvents"), toI(selection));
        (hrCounter.fill(C_CS("hNEvents"), selectionsAdditional), ...);
        if (cgAnalysis.cfgFlagQaCentrality.value) {
          hrQaCentrality.fill(C_CS("hCentralitySelection"), holderEvent.centrality, toI(selection));
          (hrQaCentrality.fill(C_CS("hCentralitySelection"), holderEvent.centrality, selectionsAdditional), ...);
        }
        return selection;
      }};

    fillEventSelection(EventSelection::All);

    if (!collision.has_foundBC()) {
      return fillEventSelection(EventSelection::Run);
    }

    const auto& foundBc{collision.template foundBC_as<aod::BCsWithTimestamps>()};
    holderEvent.runNumber = foundBc.runNumber();

    {
      const auto iter{holderCcdb.runNumbersIndicesGroupIndices.find(holderEvent.runNumber)};
      if (iter == holderCcdb.runNumbersIndicesGroupIndices.end()) {
        return fillEventSelection(EventSelection::Run);
      }
      std::tie(holderEvent.runIndex, holderEvent.runGroupIndex) = iter->second;
    }

    if (holderEvent.runGroupIndex == 0 || (cgEvent.cfgFlagRejectionRunBad.value && holderEvent.runGroupIndex < 0)) {
      return fillEventSelection(EventSelection::Run);
    }

    if (rctFlagsChecker.any() && !rctFlagsChecker.checkTable(collision)) {
      return fillEventSelection(EventSelection::Rct);
    }

    if (cgEvent.cfgFlagInel.value && !collision.isInelGt0()) {
      return fillEventSelection(EventSelection::Inel);
    }

    for (std::int32_t const& iBit : std::views::iota(0, aod::evsel::EventSelectionFlags::kNsel)) {
      if (((cgEvent.cfgBitsSelection.value >> iBit) & 1) && !collision.selection_bit(iBit)) {
        return fillEventSelection(EventSelection::Bits, NEs<EventSelection> + iBit);
      }
    }

    if (!HolderEvent::isValidCentrality(holderEvent.centrality) || !HolderEvent::isValidCentrality(holderEvent.centralityCalibration)) {
      return fillEventSelection(EventSelection::Centrality);
    }

    if (cgAnalysis.cfgFlagQaEvent.value) {
      hrQaEvent.fill(C_CS("hVxVy"), collision.posX(), collision.posY());
      hrQaEvent.fill(C_CS("hVz"), holderEvent.vz);
    }

    if (!(std::abs(holderEvent.vz) < cgEvent.cfgCutMaxAbsVz.value)) {
      return fillEventSelection(EventSelection::Vz);
    }

    if (cgEvent.cfgCutMaxOccupancy.value >= 0) {
      const std::int32_t occupancy{collision.trackOccupancyInTimeRange()};
      if (occupancy < 0 || occupancy >= cgEvent.cfgCutMaxOccupancy.value) {
        return fillEventSelection(EventSelection::Occupancy);
      }
    }

    if (cgAnalysis.cfgFlagQaRun.value) {
      hrQaRun.fill(C_CS("pRunIndexVx"), holderEvent.runIndex, collision.posX());
      hrQaRun.fill(C_CS("pRunIndexVy"), holderEvent.runIndex, collision.posY());
      hrQaRun.fill(C_CS("pRunIndexVz"), holderEvent.runIndex, holderEvent.vz);
      hrQaRun.fill(C_CS("pRunIndexMultiplicityFt0a"), holderEvent.runIndex, collision.multZeqFT0A());
      hrQaRun.fill(C_CS("pRunIndexMultiplicityFt0c"), holderEvent.runIndex, collision.multZeqFT0C());
      if (const double centrality{collision.centNTPV()}; HolderEvent::isValidCentrality(centrality)) {
        hrQaRun.fill(C_CS("pRunIndexCentralityNtpv"), holderEvent.runIndex, centrality);
      }
      if (const double centrality{collision.centFT0A()}; HolderEvent::isValidCentrality(centrality)) {
        hrQaRun.fill(C_CS("pRunIndexCentralityFt0a"), holderEvent.runIndex, centrality);
      }
      if (const double centrality{collision.centFT0C()}; HolderEvent::isValidCentrality(centrality)) {
        hrQaRun.fill(C_CS("pRunIndexCentralityFt0c"), holderEvent.runIndex, centrality);
      }
      if (const double centrality{collision.centFT0M()}; HolderEvent::isValidCentrality(centrality)) {
        hrQaRun.fill(C_CS("pRunIndexCentralityFt0m"), holderEvent.runIndex, centrality);
      }
    }

    for (const auto& track : tracks) {
      initTrack<true>(track);
    }
    for (std::int32_t const& iChargeSpecies : std::views::iota(0, NEs<ChargeSpecies>)) {
      if (holderEvent.nGlobalTracks[iChargeSpecies] > 0) {
        for (std::int32_t const& iDcaAxis : std::views::iota(0, NEs<DcaAxis>)) {
          holderEvent.measureDca[toI(DcaMeasure::Mean)][iDcaAxis][iChargeSpecies] /= holderEvent.nGlobalTracks[iChargeSpecies];
          holderEvent.measureDca[toI(DcaMeasure::Sigma)][iDcaAxis][iChargeSpecies] = std::sqrt(std::max(0., holderEvent.measureDca[toI(DcaMeasure::Sigma)][iDcaAxis][iChargeSpecies] / holderEvent.nGlobalTracks[iChargeSpecies] - std::pow(holderEvent.measureDca[toI(DcaMeasure::Mean)][iDcaAxis][iChargeSpecies], 2.)));
        }
      }
    }

    if (cgAnalysis.cfgFlagQaRun.value) {
      fillQaRunByEventByChargeSpecies<ChargeSpecies::Plus>();
      fillQaRunByEventByChargeSpecies<ChargeSpecies::Minus>();
    }

    if (cgAnalysis.cfgFlagQaEvent.value) {
      hrQaEvent.fill(C_CS("hNPvContributorsNGlobalTracks"), holderEvent.getNPvContributors(), holderEvent.getNGlobalTracks());
      if (holderEvent.getNGlobalTracks() > 0) {
        hrQaEvent.fill(C_CS("hNGlobalTracks") + C_SV(getName(DcaMeasure::Mean)) + C_CS("Dca") + C_SV(getName(DcaAxis::Xy)), holderEvent.getNGlobalTracks(), holderEvent.getMeasureDca<DcaMeasure::Mean, DcaAxis::Xy>());
        hrQaEvent.fill(C_CS("hNGlobalTracks") + C_SV(getName(DcaMeasure::Mean)) + C_CS("Dca") + C_SV(getName(DcaAxis::Z)), holderEvent.getNGlobalTracks(), holderEvent.getMeasureDca<DcaMeasure::Mean, DcaAxis::Z>());
      }
      hrQaEvent.fill(C_CS("hNTofBetaNGlobalTracks"), holderEvent.getNTofBeta(), holderEvent.getNGlobalTracks());
    }

    if (!(holderEvent.getNPvContributors() - holderEvent.getNGlobalTracks() > cgEvent.cfgCutMinDeviationNPvContributors.value)) {
      return fillEventSelection(EventSelection::NPvContributors);
    }

    readCcdb<false>();

    if (cgAnalysis.cfgFlagQaEvent.value) {
      if (holderEvent.getNGlobalTracks() > 0) {
        hrQaEvent.fill(C_CS("hNGlobalTracks") + C_SV(getName(DcaMeasure::Mean)) + C_CS("Dca") + C_SV(getName(DcaAxis::Xy)) + C_CS("_nPvContributorsCut"), holderEvent.getNGlobalTracks(), holderEvent.getMeasureDca<DcaMeasure::Mean, DcaAxis::Xy>());
        hrQaEvent.fill(C_CS("hNGlobalTracks") + C_SV(getName(DcaMeasure::Mean)) + C_CS("Dca") + C_SV(getName(DcaAxis::Z)) + C_CS("_nPvContributorsCut"), holderEvent.getNGlobalTracks(), holderEvent.getMeasureDca<DcaMeasure::Mean, DcaAxis::Z>());
      }
      hrQaEvent.fill(C_CS("hNTofBetaNGlobalTracks_nPvContributorsCut"), holderEvent.getNTofBeta(), holderEvent.getNGlobalTracks());
    }

    if (cgAnalysis.cfgFlagQaCentrality.value) {
      hrQaCentrality.fill(C_CS("hCentralityMultiplicity"), holderEvent.centrality, collision.multNTracksPVeta1());
    }

    return fillEventSelection(EventSelection::Good);
  }

  void processMc(const soa::Filtered<aod::JoinedMcCollisions>::iterator& mcCollision, const aod::McParticles& mcParticles, const soa::SmallGroups<aod::JoinedCollisionsWithMc>& collisions, const soa::Filtered<aod::JoinedTracksWithMc>& tracksUngrouped, const aod::BCsWithTimestamps&)
  {
    if (initMcEvent(mcCollision) != McEventSelection::Good) {
      return;
    }

    for (const auto& collision : collisions) {
      if (collision.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }

      const auto& tracks{tracksUngrouped.sliceBy(presliceTracksPerCollision, collision.globalIndex())};

      const EventSelection eventSelection{initEvent(collision, tracks)};
      if ((eventSelection == EventSelection::Good || toI(eventSelection) >= toI(EventSelection::NPvContributors)) && doStorageMiniTable && gRandom->Rndm() * cgEvent.cfgFactorStorageMiniTable.value < 1.) {
        const bool isGoodNPvContributors{eventSelection != EventSelection::NPvContributors};
        readCcdb<false>();
        fillMiniTableMc<ParticleNumber::Charge>(mcParticles, tracks, isGoodNPvContributors);
        fillMiniTableMc<ParticleNumber::Kaon>(mcParticles, tracks, isGoodNPvContributors);
        fillMiniTableMc<ParticleNumber::Proton>(mcParticles, tracks, isGoodNPvContributors);
      }
      if (eventSelection != EventSelection::Good) {
        continue;
      }

      if (cgAnalysis.cfgFlagQaMc.value) {
        hrQaMc.fill(C_CS("hCentralityVzMcDeltaVz"), holderEvent.centralityCalibration, holderMcEvent.vz, holderEvent.vz - holderMcEvent.vz);
      }

      if (cgAnalysis.cfgFlagQaTrack.value || cgAnalysis.cfgFlagQaDca.value || doQaAcceptance || doQaPhi || doQaPid || cgAnalysis.cfgFlagQaMc.value || doCalculationYield || doCalculationPurity || doCalculationFractionPrimary || doCalculationFluctuation) {
        if (doCalculationFluctuation) {
          holderEvent.subgroupIndex = gRandom->Integer(cgEvent.cfgNSubgroups.value);
          initCalculationFluctuation();
          holderDerivedData.clear();
        }

        for (const auto& mcParticle : mcParticles) {
          if (!initMcParticle(mcParticle)) {
            continue;
          }

          const auto& tracksMatched{tracks.sliceBy(presliceTracksPerMcParticle, mcParticle.globalIndex())};

          if (cgAnalysis.cfgFlagQaMc.value && (!cgTrack.cfgFlagMcParticlePhysicalPrimary.value || mcParticle.isPhysicalPrimary())) {
            hrQaMc.fill(C_CS("hCentralityNTracksPerMcParticle"), holderEvent.centralityCalibration, tracksMatched.size());
          }

          if (mcParticle.isPhysicalPrimary()) {
            if (doCalculationYield) {
              fillCalculationYieldByParticleSpecies<DataMode::McMcParticle, ParticleSpecies::Pion>();
              fillCalculationYieldByParticleSpecies<DataMode::McMcParticle, ParticleSpecies::Kaon>();
              fillCalculationYieldByParticleSpecies<DataMode::McMcParticle, ParticleSpecies::Proton>();
            }

            if (doCalculationFluctuation) {
              calculateFluctuationByParticleNumber<DataMode::McMcParticle, ParticleNumber::Charge>();
              calculateFluctuationByParticleNumber<DataMode::McMcParticle, ParticleNumber::Kaon>();
              calculateFluctuationByParticleNumber<DataMode::McMcParticle, ParticleNumber::Proton>();
            }
          }

          for (const auto& track : tracksMatched) {
            if (!initTrack<false>(track)) {
              continue;
            }

            if (doQaPhi) {
              fillQaPhiByParticleSpeciesAll<DataMode::McTrack, ParticleSpeciesAll::All>();
              fillQaPhiByParticleSpeciesAll<DataMode::McTrack, ParticleSpeciesAll::Pion>();
              fillQaPhiByParticleSpeciesAll<DataMode::McTrack, ParticleSpeciesAll::Kaon>();
              fillQaPhiByParticleSpeciesAll<DataMode::McTrack, ParticleSpeciesAll::Proton>();
            }

            if (doQaPid) {
              fillQaPidByParticleSpeciesAll<DataMode::McTrack, ParticleSpeciesAll::All>(track);
              fillQaPidByParticleSpeciesAll<DataMode::McTrack, ParticleSpeciesAll::Pion>(track);
              fillQaPidByParticleSpeciesAll<DataMode::McTrack, ParticleSpeciesAll::Kaon>(track);
              fillQaPidByParticleSpeciesAll<DataMode::McTrack, ParticleSpeciesAll::Proton>(track);
            }

            if (cgAnalysis.cfgFlagQaMc.value && (!cgTrack.cfgFlagMcParticlePhysicalPrimary.value || mcParticle.isPhysicalPrimary())) {
              hrQaMc.fill(C_CS("hCentralityPtMcEtaMcDeltaPt"), holderEvent.centralityCalibration, holderMcParticle.pt, holderMcParticle.eta, holderTrack.pt - holderMcParticle.pt);
              hrQaMc.fill(C_CS("hCentralityPtMcEtaMcDeltaEta"), holderEvent.centralityCalibration, holderMcParticle.pt, holderMcParticle.eta, holderTrack.eta - holderMcParticle.eta);
            }

            if (doCalculationYield && (!cgTrack.cfgFlagMcParticlePhysicalPrimary.value || mcParticle.isPhysicalPrimary())) {
              fillCalculationYieldByParticleSpecies<DataMode::McTrack, ParticleSpecies::Pion>();
              fillCalculationYieldByParticleSpecies<DataMode::McTrack, ParticleSpecies::Kaon>();
              fillCalculationYieldByParticleSpecies<DataMode::McTrack, ParticleSpecies::Proton>();
            }

            if (doCalculationPurity && (!cgTrack.cfgFlagMcParticlePhysicalPrimary.value || mcParticle.isPhysicalPrimary())) {
              fillCalculationPurityByParticleSpecies<ParticleSpecies::Pion>();
              fillCalculationPurityByParticleSpecies<ParticleSpecies::Kaon>();
              fillCalculationPurityByParticleSpecies<ParticleSpecies::Proton>();
            }

            if (doCalculationFractionPrimary) {
              fillCalculationFractionPrimaryByParticleSpecies<ParticleSpecies::Pion>(mcParticle);
              fillCalculationFractionPrimaryByParticleSpecies<ParticleSpecies::Kaon>(mcParticle);
              fillCalculationFractionPrimaryByParticleSpecies<ParticleSpecies::Proton>(mcParticle);
            }

            if (doCalculationFluctuation && (!cgTrack.cfgFlagMcParticlePhysicalPrimary.value || mcParticle.isPhysicalPrimary())) {
              calculateFluctuationByParticleNumber<DataMode::McTrack, ParticleNumber::Charge>();
              calculateFluctuationByParticleNumber<DataMode::McTrack, ParticleNumber::Kaon>();
              calculateFluctuationByParticleNumber<DataMode::McTrack, ParticleNumber::Proton>();
            }
          }
        }

        if (doCalculationFluctuation) {
          fillCalculationFluctuationByParticleNumber<ParticleNumber::Charge>();
          fillCalculationFluctuationByParticleNumber<ParticleNumber::Kaon>();
          fillCalculationFluctuationByParticleNumber<ParticleNumber::Proton>();
          if (const std::optional<aod::tiny_collision::Code::type> codeCollision{mini_collision_codec::encode(cgEvent.cfgFlagMcCollisionVz.value ? holderMcEvent.vz : holderEvent.vz, holderEvent.centrality, true)}; codeCollision.has_value()) {
            tinyMcCollision(holderDerivedData.nMcParticles[toI(ChargeSpecies::Plus)], holderDerivedData.nMcParticles[toI(ChargeSpecies::Minus)]);
            tinyCollision(*codeCollision, holderDerivedData.nTracks[toI(ChargeSpecies::Plus)], holderDerivedData.nTracks[toI(ChargeSpecies::Minus)]);
            for (aod::tiny_mc_particle::SignedEfficiency::type const& signedEfficiency : holderDerivedData.signedEfficienciesMcParticle) {
              tinyMcParticle(tinyMcCollision.lastIndex(), signedEfficiency);
            }
            for (aod::tiny_track::SignedEfficiency::type const& signedEfficiency : holderDerivedData.signedEfficienciesTrack) {
              tinyTrack(tinyCollision.lastIndex(), signedEfficiency);
            }
          }
        }
      }
    }
  }

  void processRaw(const soa::Filtered<aod::JoinedCollisions>::iterator& collision, const soa::Filtered<aod::JoinedTracks>& tracks, const aod::BCsWithTimestamps&)
  {
    const EventSelection eventSelection{initEvent(collision, tracks)};
    if ((eventSelection == EventSelection::Good || toI(eventSelection) >= toI(EventSelection::NPvContributors)) && doStorageMiniTable && gRandom->Rndm() * cgEvent.cfgFactorStorageMiniTable.value < 1.) {
      const bool isGoodNPvContributors{eventSelection != EventSelection::NPvContributors};
      readCcdb<false>();
      fillMiniTableRaw<ParticleNumber::Charge>(tracks, isGoodNPvContributors);
      fillMiniTableRaw<ParticleNumber::Kaon>(tracks, isGoodNPvContributors);
      fillMiniTableRaw<ParticleNumber::Proton>(tracks, isGoodNPvContributors);
    }
    if (eventSelection != EventSelection::Good) {
      return;
    }

    if (!cgAnalysis.cfgFlagQaTrack.value && !cgAnalysis.cfgFlagQaDca.value && !doQaAcceptance && !doQaPhi && !doQaPid && !doCalculationYield && !doCalculationFluctuation) {
      return;
    }

    if (doCalculationFluctuation) {
      holderEvent.subgroupIndex = gRandom->Integer(cgEvent.cfgNSubgroups.value);
      initCalculationFluctuation();
      holderDerivedData.clear();
    }

    for (const auto& track : tracks) {
      if (!initTrack<false>(track)) {
        continue;
      }

      if (doQaPhi) {
        fillQaPhiByParticleSpeciesAll<DataMode::RawTrack, ParticleSpeciesAll::All>();
        fillQaPhiByParticleSpeciesAll<DataMode::RawTrack, ParticleSpeciesAll::Pion>();
        fillQaPhiByParticleSpeciesAll<DataMode::RawTrack, ParticleSpeciesAll::Kaon>();
        fillQaPhiByParticleSpeciesAll<DataMode::RawTrack, ParticleSpeciesAll::Proton>();
      }

      if (doQaPid) {
        fillQaPidByParticleSpeciesAll<DataMode::RawTrack, ParticleSpeciesAll::All>(track);
        fillQaPidByParticleSpeciesAll<DataMode::RawTrack, ParticleSpeciesAll::Pion>(track);
        fillQaPidByParticleSpeciesAll<DataMode::RawTrack, ParticleSpeciesAll::Kaon>(track);
        fillQaPidByParticleSpeciesAll<DataMode::RawTrack, ParticleSpeciesAll::Proton>(track);
      }

      if (doCalculationYield) {
        fillCalculationYieldByParticleSpecies<DataMode::RawTrack, ParticleSpecies::Pion>();
        fillCalculationYieldByParticleSpecies<DataMode::RawTrack, ParticleSpecies::Kaon>();
        fillCalculationYieldByParticleSpecies<DataMode::RawTrack, ParticleSpecies::Proton>();
      }

      if (doCalculationFluctuation) {
        calculateFluctuationByParticleNumber<DataMode::RawTrack, ParticleNumber::Charge>();
        calculateFluctuationByParticleNumber<DataMode::RawTrack, ParticleNumber::Kaon>();
        calculateFluctuationByParticleNumber<DataMode::RawTrack, ParticleNumber::Proton>();
      }
    }

    if (doCalculationFluctuation) {
      fillCalculationFluctuationByParticleNumber<ParticleNumber::Charge>();
      fillCalculationFluctuationByParticleNumber<ParticleNumber::Kaon>();
      fillCalculationFluctuationByParticleNumber<ParticleNumber::Proton>();
      if (const std::optional<aod::tiny_collision::Code::type> codeCollision{mini_collision_codec::encode(holderEvent.vz, holderEvent.centrality, true)}; codeCollision.has_value()) {
        tinyCollision(*codeCollision, holderDerivedData.nTracks[toI(ChargeSpecies::Plus)], holderDerivedData.nTracks[toI(ChargeSpecies::Minus)]);
        for (aod::tiny_track::SignedEfficiency::type const& signedEfficiency : holderDerivedData.signedEfficienciesTrack) {
          tinyTrack(tinyCollision.lastIndex(), signedEfficiency);
        }
      }
    }
  }

  PROCESS_SWITCH_FULL(PartNumFluc, processMc, ProcessMc, "Flag of processing MC data", false);
  PROCESS_SWITCH_FULL(PartNumFluc, processRaw, ProcessRaw, "Flag of processing raw data", true);
};

WorkflowSpec defineDataProcessing(const ConfigContext& configContext)
{
  return WorkflowSpec{adaptAnalysisTask<PartNumFluc>(configContext)};
}
