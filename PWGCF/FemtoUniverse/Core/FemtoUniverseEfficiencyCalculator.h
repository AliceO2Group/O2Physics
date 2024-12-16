#ifndef PWGCF_FEMTOUNIVERSE_CORE_EFFICIENCY_CALCULATOR_H_
#define PWGCF_FEMTOUNIVERSE_CORE_EFFICIENCY_CALCULATOR_H_

#include "Framework/Configurable.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/CallbackService.h"
#include "Framework/InitContext.h"
#include "CCDB/BasicCCDBManager.h"
#include "FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

namespace o2::analysis::femto_universe::efficiency
{
template <uint8_t T>
concept IsOneOrTwo = T == 1 || T == 2;

struct EfficiencyConfigurableGroup : ConfigurableGroup {
  Configurable<bool> shouldCalculate{"ConfEfficiencyCalculate", false, "Should calculate efficiency"};
  Configurable<bool> shouldApplyCorrections{"ConfEfficiencyApplyCorrections", false, "Should apply corrections from efficiency"};

  Configurable<std::vector<std::string>> ccdbLabels{"ConfCCDBLabels", {"", ""}, "Labels for efficiency objects in CCDB"};

  OutputObj<TH1F> hEfficiency1{"Efficiency part1"};
  OutputObj<TH1F> hEfficiency2{"Efficiency part2"};

  // TODO: move to separate struct?
  Configurable<std::string> ccdbUrl{"ConfCCDBUrl", "http://alice-ccdb.cern.ch", "CCDB URL to be used"};
  Configurable<std::string> ccdbPath{"ConfCCDBPath", "", "CCDB base path to where to upload objects"};
  Configurable<long> ccdbTimestamp{"ConfCCDBTimestamp", -1, "Timestamp from which to query CCDB objects"};

  // TODO: should this be decalred in task directly?
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 1> hMCTruth1;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 2> hMCTruth2;
};

class EfficiencyCalculator
{
 public:
  o2::ccdb::BasicCCDBManager& ccdb{o2::ccdb::BasicCCDBManager::instance()};

  explicit EfficiencyCalculator(EfficiencyConfigurableGroup* config) : config(config)
  {
  }

  auto init() -> void
  {
    ccdbApi.init(config->ccdbUrl);
    ccdb.setURL(config->ccdbUrl);
    ccdb.setLocalObjectValidityChecking();
    ccdb.setFatalWhenNull(false);

    shouldCalculate = config->shouldCalculate;
    shouldApplyCorrections = config->shouldApplyCorrections;

    ccdbFullPath = std::format("{}/{}", config->ccdbPath.value, folderName);

    if (config->shouldCalculate) {
      hOutput = {config->hEfficiency1.object, config->hEfficiency2.object};
    }

    if (config->shouldApplyCorrections) {
      hLoaded = {
        loadEfficiencyFromCCDB<1>(config->ccdbTimestamp),
        loadEfficiencyFromCCDB<2>(config->ccdbTimestamp) //
      };
    }
  }

  template <uint8_t N>
    requires IsOneOrTwo<N>
  auto doMCTruth(auto particles) const -> void
  {
    for (const auto& particle : particles) {
      if constexpr (N == 1) {
        config->hMCTruth1.fillQA<false, false>(particle);
      } else if constexpr (N == 2) {
        config->hMCTruth2.fillQA<false, false>(particle);
      }
    }
  }

  auto uploadOnStop(InitContext& ic) -> EfficiencyCalculator&
  {
    if (!shouldUploadOnStop && config->shouldCalculate) {
      shouldUploadOnStop = true;

      auto& callbacks = ic.services().get<CallbackService>();

      callbacks.set<o2::framework::CallbackService::Id::Stop>([this]() {
        for (auto i = 0UL; i < hOutput.size(); i++) {
          const auto& output = hOutput[i];

          if (isHistogramEmpty(output.get())) {
            LOGF(error, log("Histogram %d is empty - save aborted"), i + 1);
            return;
          }
          LOGF(debug, log("Found histogram %d: %s"), i + 1, output->GetTitle());

          long now = duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
          long oneYear = 365LL * 24 * 60 * 60 * 1000;

          if (ccdbApi.storeAsTFileAny(output.get(), ccdbFullPath, createMetadata(i), now, now + oneYear) == 0) {
            LOGF(info, log("Histogram %d saved successfully"), i + 1);
          } else {
            LOGF(fatal, log("Histogram %d save failed"), i + 1);
          }
        }
      });
    } else {
      LOGF(warn, log("Uploading on stop callback is already set up"));
    }

    return *this;
  }

  template <uint8_t N>
    requires IsOneOrTwo<N>
  auto getWeight(auto const& particle) const -> float
  {
    auto weight = 1.0f;
    auto hEff = hLoaded[N - 1];

    if (shouldApplyCorrections && hEff) {
      auto bin = hEff->FindBin(particle.pt());
      auto eff = hEff->GetBinContent(bin);
      weight /= eff > 0 ? eff : 1.0f;
    }

    return weight;
  }

  template <uint8_t N>
    requires IsOneOrTwo<N>
  auto calculate(std::shared_ptr<TH1> truth, std::shared_ptr<TH1> reco) const
  {
    if (!shouldCalculate) {
      return;
    }

    if (!truth || !reco) {
      LOGF(error, log("MC Truth & MC Reco histograms cannot be null"));
      return;
    }

    auto hEff = hOutput[N - 1];
    if (!hEff) {
      LOGF(error, log("No OutputObj specified for particle %d histogram"), N);
      return;
    }

    for (auto bin = 0; bin < hEff->GetNbinsX(); bin++) {
      auto denom = truth->GetBinContent(bin);
      hEff->SetBinContent(bin, denom == 0 ? 0 : reco->GetBinContent(bin) / denom);
    }
  }

 private:
  static inline auto log(const std::string& msg) -> const std::string
  {
    return std::format("[EFFICIENCY] {}", msg);
  }

  static auto isHistogramEmpty(TH1* hist) -> bool
  {
    if (!hist) {
      return true;
    }

    // check overflow bins as well
    for (auto idx = 0; idx <= hist->GetNbinsX() + 1; idx++) {
      if (hist->GetBinContent(idx) != 0) {
        return false;
      }
    }

    return true;
  }

  auto createMetadata(uint8_t partNo) const -> std::map<std::string, std::string>
  {
    if (config->ccdbLabels->size() != 2) {
      LOGF(fatal, log("CCDB labels configurable should be exactly of size 2"));
    }
    return std::map<std::string, std::string>{
      {"label", config->ccdbLabels.value[partNo]} //
    };
  }

  template <uint8_t N>
    requires IsOneOrTwo<N>
  auto loadEfficiencyFromCCDB(long timestamp) const -> TH1*
  {
    auto hEff = ccdb.getSpecific<TH1>(ccdbFullPath, timestamp, createMetadata(N - 1));
    if (!hEff || hEff->IsZombie()) {
      LOGF(error, log("Could not load histogram from %s"), config->ccdbPath.value);
      return nullptr;
    }

    LOGF(info, log("Histogram \"%s\" loaded from \"%s\""), hEff->GetTitle(), config->ccdbPath.value);
    return hEff;
  }

  EfficiencyConfigurableGroup* config{};

  bool shouldUploadOnStop = false;
  bool shouldCalculate = false;
  bool shouldApplyCorrections = false;

  std::array<std::shared_ptr<TH1>, 2> hOutput{};
  std::array<TH1*, 2> hLoaded{};

  o2::ccdb::CcdbApi ccdbApi{};
  std::string ccdbFullPath{};

  static constexpr std::string folderName{"Efficiency"};
  static constexpr std::array<std::string, 3> histSuffix{"", "_one", "_two"};
};

} // namespace o2::analysis::femto_universe::efficiency

#endif
