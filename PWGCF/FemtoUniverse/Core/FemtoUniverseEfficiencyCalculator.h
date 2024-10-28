#ifndef PWGCF_FEMTOUNIVERSE_CORE_EFFICIENCY_CALCULATOR_H_
#define PWGCF_FEMTOUNIVERSE_CORE_EFFICIENCY_CALCULATOR_H_

#include <concepts>
#include <type_traits>

#include "Framework/Configurable.h"
#include "FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

namespace o2::analysis::femto_universe
{

template <typename T, typename FieldType, auto FieldPtr>
concept HasField = requires(T t) {
  requires std::is_same_v<decltype(FieldPtr), FieldType T::*>;
  { t.*FieldPtr } -> std::convertible_to<FieldType>;
};

template <uint8_t T>
concept IsOneOrTwo = T == 1 || T == 2;

template <typename T>
concept IsConfGroup = std::is_base_of_v<o2::framework::ConfigurableGroup, T>;

using HistGenPart = FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 1>;

template <typename Task>
struct EfficiencyCalculator {
  Task& task;

  std::array<HistGenPart, 2> histGenPart;

  EfficiencyCalculator(Task& task);

  template <uint8_t Index>
    requires IsOneOrTwo<Index>
  auto doMCGen(auto particle) -> void;
};

} // namespace o2::analysis::femto_universe

#endif
