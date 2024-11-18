#ifndef PWGCF_FEMTOUNIVERSE_CORE_EFFICIENCY_CALCULATOR_H_
#define PWGCF_FEMTOUNIVERSE_CORE_EFFICIENCY_CALCULATOR_H_

#include <iostream>

#include "CCDB/BasicCCDBManager.h"
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

  const std::string& folderPrefix{"Efficiency"};

  o2::ccdb::CcdbApi ccdbApi;
  o2::ccdb::BasicCCDBManager& ccdb{o2::ccdb::BasicCCDBManager::instance()};

  // const std::string ccdbUrl{"http://alice-ccdb.cern.ch"};
  const std::string& ccdbUrl{"http://ccdb-test.cern.ch:8080"};
  const std::string& ccdbPath{"Users/d/dkarpins/efficiency"};

  MCTruthHist<1> mcTruthHist1;
  MCRecoHist<1>& mcRecoHist1;

  MCTruthHist<2> mcTruthHist2;
  MCRecoHist<2>& mcRecoHist2;

  EfficiencyCalculator(Task* task, int PDG1, int PDG2, MCRecoHist<1>& mcRecoHist1, MCRecoHist<2>& mcRecoHist2)
    : task(task), PDG1(PDG1), PDG2(PDG2), mcRecoHist1(mcRecoHist1), mcRecoHist2(mcRecoHist2), registry(task->efficiencyRegistry)
  {
    mcTruthHist1.init(&registry, task->ConfTempFitVarpTBins, task->ConfTempFitVarPDGBins, false, PDG1, false);
    mcTruthHist2.init(&registry, task->ConfTempFitVarpTBins, task->ConfTempFitVarPDGBins, false, PDG2, false);

    registry.add("Efficiency/part1", "Efficiency origin/generated ; p_{T} (GeV/c); Efficiency", HistType::kTH1F, {{100, 0, 4}});
    registry.add("Efficiency/part2", "Efficiency origin/generated ; p_{T} (GeV/c); Efficiency", HistType::kTH1F, {{100, 0, 4}});

    ccdbApi.init(ccdbUrl);
  }

  // ~EfficiencyCalculator()
  // {
  //   std::cout << "hello world!\n";
  // }

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

  auto getEfficiency() -> TH1*
  {
    auto fullPath{std::format("%s/%s", ccdbUrl, ccdbPath)};
    auto now{std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()};

    auto* efficiency{ccdb.getForTimeStamp<TH1>(fullPath, now)};
    if (!efficiency || efficiency->IsZombie()) {
      LOGF(fatal, "Could not load efficiency protoneff histogram from %s", fullPath);
    }

    // float weight{1.0f};
    // weight = protoneff->GetBinContent(protoneff->FindBin(track.pt(), track.eta())) * phieff->GetBinContent(phieff->FindBin(phicandidate.pt(), phicandidate.eta()));
    return efficiency;
  }

  auto
    save() -> void
  {
    std::cout << "inside save" << std::endl;
    std::shared_ptr<TH1> efficiencyHist1{registry.get<TH1>(HIST("Efficiency/part1"))};
    std::map<string, string> metadata;
    // long now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    // long oneYear = now + (365LL * 24 * 60 * 60 * 1000);
    // metadata["runNumber"] = runNumber;
    ccdbApi.storeAsTFileAny(efficiencyHist1.get(), ccdbPath, metadata);
  }

  template <uint8_t N>
    requires IsOneOrTwo<N>
  auto doMCTruth(int pdg, auto particles) -> void
  {
    for (const auto& particle : particles) {
      if (static_cast<int>(particle.pidcut()) != pdg) {
        continue;
      }

      if constexpr (N == 1) {
        mcTruthHist1.fillQA<false, false>(particle);
      } else if constexpr (N == 2) {
        mcTruthHist2.fillQA<false, false>(particle);
      }
    }
  }
};

} // namespace o2::analysis::femto_universe

#endif
