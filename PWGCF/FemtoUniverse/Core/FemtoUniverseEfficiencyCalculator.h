#ifndef PWGCF_FEMTOUNIVERSE_CORE_EFFICIENCY_CALCULATOR_H_
#define PWGCF_FEMTOUNIVERSE_CORE_EFFICIENCY_CALCULATOR_H_

#include <CCDB/BasicCCDBManager.h>
#include <Framework/AnalysisHelpers.h>
#include <iostream>

#include "FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

namespace o2::analysis::femto_universe
{

template <uint8_t T>
concept IsOneOrTwo = T == 1 || T == 2;

template <int N>
  requires IsOneOrTwo<N>
using MCTruthHist = FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, N>;

template <int N>
  requires IsOneOrTwo<N>
using MCRecoHist = FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, N>;

template <typename Task>
struct EfficiencyCalculator {
  Task* task;
  HistogramRegistry& registry;

  int PDG1;
  int PDG2;

  const std::string_view folderPrefix = "Efficiency";

  Service<ccdb::BasicCCDBManager> ccdb;

  MCTruthHist<1> mcTruthHist1;
  MCRecoHist<1>& mcRecoHist1;

  MCTruthHist<2> mcTruthHist2;
  MCRecoHist<2>& mcRecoHist2;

  EfficiencyCalculator(Task* task, int PDG1, int PDG2, MCRecoHist<1>& mcRecoHist1, MCRecoHist<2>& mcRecoHist2)
    : task(task), PDG1(PDG1), PDG2(PDG2), mcRecoHist1(mcRecoHist1), mcRecoHist2(mcRecoHist2), registry(task->efficiencyRegistry)
  {
    mcTruthHist1.init(&registry, task->ConfTempFitVarpTBins, task->ConfTempFitVarPDGBins, 0, PDG1, false);
    mcTruthHist2.init(&registry, task->ConfTempFitVarpTBins, task->ConfTempFitVarPDGBins, 0, PDG2, false);

    registry.add("Efficiency/part1", "Efficiency origin/generated ; p_{T} (GeV/c); Efficiency", HistType::kTH1F, {{100, 0, 4}});
    registry.add("Efficiency/part2", "Efficiency origin/generated ; p_{T} (GeV/c); Efficiency", HistType::kTH1F, {{100, 0, 4}});
  }

  ~EfficiencyCalculator()
  {
    std::cout << "hello world!\n";
  }

  auto calculate() -> void
  {
    std::shared_ptr<TH1> efficiencyHist1{registry.get<TH1>(HIST("Efficiency/part1"))};

    std::shared_ptr<TH1> truth1{registry.get<TH1>(HIST("MCTruthTracks_one/hPt"))};
    std::shared_ptr<TH1> reco1{task->qaRegistry.template get<TH1>(HIST("Tracks_one_MC/hPt"))};

    for (int bin{0}; bin < efficiencyHist1->GetNbinsX(); bin++) {
      auto denom{truth1->GetBinContent(bin)};
      efficiencyHist1->SetBinContent(bin, denom == 0 ? 0 : reco1->GetBinContent(bin) / denom);
    }
  }

  auto save() -> void
  {
  }

  template <uint8_t N>
    requires IsOneOrTwo<N>
  auto fillMCTruth(auto particle) -> void
  {
    if constexpr (N == 1) {
      mcTruthHist1.fillQA<false, false>(particle);
    } else if constexpr (N == 2) {
      mcTruthHist2.fillQA<false, false>(particle);
    }
  }
};

} // namespace o2::analysis::femto_universe

#endif
