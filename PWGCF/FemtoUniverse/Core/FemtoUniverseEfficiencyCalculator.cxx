#include "FemtoUniverseEfficiencyCalculator.h"

using namespace o2::analysis::femtoUniverse;

template <typename Task>
EfficiencyCalculator<Task>::EfficiencyCalculator(Task& task) : task(task)
{
  // assume the name of registry in specified task
  const auto& registry{task.efficiencyRegistry};

  registry.add("MCTruth/hPt", " ;#it{p}_{T} (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 4}}});

  for (uint8_t idx{1}; idx <= 2; idx++) {
    histGenPart[idx].init(registry, task.ConfTempFitVarpTBins, task.ConfTempFitVarPDGBins, 0, task.ConfPDGCodePartOne, false);
    registry.add(
      HIST("part") + std::to_string(idx),
      "Efficiency origin/generated ; p_{T} (GeV/c); Efficiency",
      {HistType::kTH1F, {{100, 0, 4}}} //
    );
  }
}

/// This function processes the same event and takes care of all the histogramming
template <typename Task>
template <uint8_t Index>
  requires IsOneOrTwo<Index>
auto EfficiencyCalculator<Task>::doMCGen(auto particle) -> void
{
  if (!task.ConfNoPDGPartOne && (static_cast<int>(particle.pidcut()) != task.ConfPDGCodePartOne)) {
    return;
  }
  histGenPart[Index].fillQA<false, false>(particle);
}
