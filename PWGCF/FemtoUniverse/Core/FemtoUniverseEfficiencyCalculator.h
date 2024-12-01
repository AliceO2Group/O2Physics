#ifndef PWGCF_FEMTOUNIVERSE_CORE_EFFICIENCY_CALCULATOR_H_
#define PWGCF_FEMTOUNIVERSE_CORE_EFFICIENCY_CALCULATOR_H_

#include "Framework/CallbackService.h"
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
 public:
  EfficiencyCalculator(Task* task, HistogramRegistry* registry, int PDG, TH1F* output) : task(task), registry(registry), PDG(PDG), output(output)
  {
    assert(task && "Task cannot be null");

    hMCTruth.init(registry, task->ConfTempFitVarpTBins, task->ConfTempFitVarPDGBins, false, PDG, false);

    ccdbApi.init(ccdbUrl);

    ccdb.setURL(ccdbUrl);
    ccdb.setCaching(true);
    ccdb.setLocalObjectValidityChecking();

    hEff = loadEfficiencyFromCCDB(1732467686456);
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

  auto saveOnStop(InitContext& ic) -> void
  {
    if (!shouldSaveOnStop) {
      shouldSaveOnStop = true;

      auto& callbacks = ic.services().get<CallbackService>();
      callbacks.set<o2::framework::CallbackService::Id::Stop>([this]() {
        save();
      });
    } else {
      LOG(warn) << "[EFFICIENCY] Save on stop is already set up";
    }
  }

  auto getWeight(auto const& particle) const -> float
  {
    if (auto h{*hEff}) {
      auto eff = h->GetBinContent(h->FindBin(particle.pt()));
      return eff > 0 ? 1.0f / eff : 0.0f;
    }
    return 1.0f;
  }

  auto save() -> void
  {
    LOG(debug) << "[EFFICIENCY] Attempting save to CCDB...";

    if (!output) {
      LOG(error) << "[EFFICIENCY] No histogram to save";
      return;
    }

    LOG(debug) << "[EFFICIENCY] Found histogram: " << output->GetTitle();

    std::map<string, string> metadata{};
    // metadata["runNumber"] = runNumber;

    long now{std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()};
    long oneYear{365LL * 24 * 60 * 60 * 1000};

    if (ccdbApi.storeAsTFileAny(output, ccdbPath, metadata, now, now + oneYear) != 0) {
      LOG(fatal) << "[EFFICIENCY] Histogram save failed";
      return;
    }

    LOG(info) << "[EFFICIENCY] Histogram saved successfully";
  }

 private:
  Task* task{};
  HistogramRegistry* registry{};

  int PDG;
  TH1* output{};
  std::optional<TH1*> hEff{};

  const std::string folderPrefix{"Efficiency"};

  o2::ccdb::CcdbApi ccdbApi;
  o2::ccdb::BasicCCDBManager& ccdb{o2::ccdb::BasicCCDBManager::instance()};

  // const std::string& ccdbUrl{"http://alice-ccdb.cern.ch"};
  const std::string ccdbUrl{"http://ccdb-test.cern.ch:8080"};
  const std::string ccdbPath{"Users/d/dkarpins/" + folderPrefix};

  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 1> hMCTruth;

  bool shouldSaveOnStop{false};

  auto loadEfficiencyFromCCDB(long timestamp) -> std::optional<TH1*>
  {
    auto hEff{ccdb.getForTimeStamp<TH1>(ccdbPath, timestamp)};
    if (!hEff || hEff->IsZombie()) {
      LOGF(fatal, "[EFFICIENCY] Could not load histogram from %s", ccdbPath);
      return {};
    }

    LOG(info) << "[EFFICIENCY] Histogram " << output->GetTitle() << "loaded from " << ccdbPath;
    return hEff;
  }
};

} // namespace o2::analysis::femto_universe

#endif
