#ifndef PWGCF_FEMTOUNIVERSE_CORE_EFFICIENCY_CALCULATOR_H_
#define PWGCF_FEMTOUNIVERSE_CORE_EFFICIENCY_CALCULATOR_H_

#include "Framework/EndOfStreamContext.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/InitContext.h"
#include "CCDB/BasicCCDBManager.h"
#include "FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

namespace o2::analysis::femto_universe
{

template <typename Task>
class EfficiencyCalculator
{
 private:
  Task* task{};
  HistogramRegistry* registry{};

  int PDG;
  TH1* output{};

  const std::string folderPrefix{"Efficiency"};

  o2::ccdb::CcdbApi ccdbApi;
  o2::ccdb::BasicCCDBManager& ccdb{o2::ccdb::BasicCCDBManager::instance()};

  // const std::string& ccdbUrl{"http://alice-ccdb.cern.ch"};
  const std::string ccdbUrl{"http://ccdb-test.cern.ch:8080"};
  const std::string ccdbPath{"Users/d/dkarpins/" + folderPrefix};

  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 1> hMCTruth;

 public:
  EfficiencyCalculator(Task* task, HistogramRegistry* registry, int PDG, TH1F* output) : task(task), registry(registry), PDG(PDG), output(output)
  {
    assert(task && "Task cannot be null");

    hMCTruth.init(registry, task->ConfTempFitVarpTBins, task->ConfTempFitVarPDGBins, false, PDG, false);

    ccdbApi.init(ccdbUrl);

    ccdb.setURL(ccdbUrl);
    ccdb.setCaching(true);
    ccdb.setLocalObjectValidityChecking();
  }

  auto getWeight(auto const& particle) const -> float
  {
    auto now{std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()};

    auto* hEff{ccdb.getForTimeStamp<TH1>(ccdbPath, now)};
    if (!hEff || hEff->IsZombie()) {
      LOGF(fatal, "[EFFICIENCY] Could not load histogram from %s", ccdbPath);
      return 1.0f;
    }

    LOG(info) << "[EFFICIENCY] Histogram " << output->GetTitle() << "loaded from " << ccdbPath;

    auto eff = hEff->GetBinContent(hEff->FindBin(particle.pt()));
    return eff == 0 ? 0 : 1 / eff;
  }

  auto save() -> void
  {
    LOG(debug) << "[EFFICIENCY] Attempting save to CCDB...";

    if (!output) {
      LOG(error) << "[EFFICIENCY] No histogram to save";
      return;
    }

    LOG(debug) << "[EFFICIENCY] Found histogram: " << output->GetTitle();

    // long now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    // long oneYear = now + (365LL * 24 * 60 * 60 * 1000);
    // metadata["runNumber"] = runNumber;

    std::map<string, string> metadata;
    if (ccdbApi.storeAsTFileAny(output, ccdbPath, metadata) != 0) {
      LOG(fatal) << "[EFFICIENCY] Histogram save failed";
      return;
    }

    LOG(info) << "[EFFICIENCY] Histogram saved successfully";
  }

  auto doMCTruth(int pdg, auto particles) -> void
  {
    for (const auto& particle : particles) {
      if (static_cast<int>(particle.pidcut()) != pdg) {
        continue;
      }

      hMCTruth.fillQA<false, false>(particle);
    }
  }

  auto saveAfterProcess(InitContext& ic) -> void
  {
    auto& callbacks = ic.services().get<CallbackService>();
    callbacks.set<o2::framework::CallbackService::Id::EndOfStream>([this](EndOfStreamContext&) {
      LOG(info) << "### callback: EndOfStream";
      save();
    });
  }
};

} // namespace o2::analysis::femto_universe

#endif
