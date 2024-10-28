#ifndef PWGCF_FEMTOUNIVERSE_CORE_EFFICIENCY_CALCULATOR_H_
#define PWGCF_FEMTOUNIVERSE_CORE_EFFICIENCY_CALCULATOR_H_

#include <concepts>
#include <type_traits>

#include "Framework/Configurable.h"
#include "FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

namespace o2::analysis::femto_universe
{

template <uint8_t T>
concept IsOneOrTwo = T == 1 || T == 2;

using MCTruthHist = FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 1>;

template <typename Task>
struct EfficiencyCalculator {
  Task& task;
  std::array<MCTruthHist, 2> mcTruthHist;

  EfficiencyCalculator(Task& task) : task(task)
  {
    // assume the name of registry in specified task
    const auto& registry{task.efficiencyRegistry};

    registry.add(
      "MCTruth/hPt",
      " ;#it{p}_{T} (GeV/c); Entries",
      {HistType::kTH1F, {{100, 0, 4}}} //
    );

    for (uint8_t idx{1}; idx <= 2; idx++) {
      mcTruthHist[idx].init(
        registry,
        task.ConfTempFitVarpTBins,
        task.ConfTempFitVarPDGBins,
        0,
        task.ConfPDGCodePartOne,
        false //
      );

      registry.add(
        HIST("part") + std::to_string(idx),
        "Efficiency origin/generated ; p_{T} (GeV/c); Efficiency",
        {HistType::kTH1F, {{100, 0, 4}}} //
      );
    }
  }

  template <uint8_t Index>
    requires IsOneOrTwo<Index>
  auto doMCGen(auto particle) -> void
  {
    if (!task.ConfNoPDGPartOne && (static_cast<int>(particle.pidcut()) != task.ConfPDGCodePartOne)) {
      return;
    }
    histGenPart[Index].fillQA<false, false>(particle);
  }
};

} // namespace o2::analysis::femto_universe

#endif
