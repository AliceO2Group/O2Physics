#ifndef PWGCF_FEMTOUNIVERSE_CORE_EFFICIENCY_CALCULATOR_H_
#define PWGCF_FEMTOUNIVERSE_CORE_EFFICIENCY_CALCULATOR_H_

#include "CCDB/BasicCCDBManager.h"
#include "FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

namespace o2::analysis::femto_universe
{

using MCTruthHist = FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 1>;

template <typename Task>
struct EfficiencyCalculator {
  Task* task;
  HistogramRegistry& registry{task->efficiencyRegistry};

  int PDG;

  const std::string& folderPrefix{"Efficiency"};

  o2::ccdb::CcdbApi ccdbApi;
  o2::ccdb::BasicCCDBManager& ccdb{o2::ccdb::BasicCCDBManager::instance()};

  // const std::string ccdbUrl{"http://alice-ccdb.cern.ch"};
  const std::string& ccdbUrl{"http://ccdb-test.cern.ch:8080"};
  const std::string& ccdbPath{"Users/d/dkarpins/efficiency"};

  MCTruthHist mcTruthHist;

  EfficiencyCalculator(Task* task, int PDG) : task(task), PDG(PDG)
  {
    mcTruthHist.init(&registry, task->ConfTempFitVarpTBins, task->ConfTempFitVarPDGBins, false, PDG, false);
    registry.add("Efficiency/part1", "Efficiency origin/generated ; p_{T} (GeV/c); Efficiency", HistType::kTH1F, {{100, 0, 4}});

    ccdbApi.init(ccdbUrl);
  }

  // ~EfficiencyCalculator()
  // {
  //   std::cout << "hello world!\n";
  // }

  auto getEfficiency() -> TH1*
  {
    auto now{std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()};

    auto* efficiency{ccdb.getForTimeStamp<TH1>(ccdbPath, now)};
    if (!efficiency || efficiency->IsZombie()) {
      LOGF(fatal, "Could not load efficiency protoneff histogram from %s", ccdbPath);
    }

    LOG(info) << efficiency->GetTitle();

    // float weight{1.0f};
    // weight = protoneff->GetBinContent(protoneff->FindBin(track.pt(), track.eta())) * phieff->GetBinContent(phieff->FindBin(phicandidate.pt(), phicandidate.eta()));
    return efficiency;
  }

  auto save() -> void
  {
    LOG(info) << "saving...";
    std::shared_ptr<TH1> efficiencyHist{registry.get<TH1>(HIST("Efficiency/part1"))};
    std::map<string, string> metadata;
    // long now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    // long oneYear = now + (365LL * 24 * 60 * 60 * 1000);
    // metadata["runNumber"] = runNumber;
    ccdbApi.storeAsTFileAny(efficiencyHist.get(), ccdbPath, metadata);
  }

  auto doMCTruth(int pdg, auto particles) -> void
  {
    for (const auto& particle : particles) {
      if (static_cast<int>(particle.pidcut()) != pdg) {
        continue;
      }

      mcTruthHist.fillQA<false, false>(particle);
    }
  }
};

} // namespace o2::analysis::femto_universe

#endif
