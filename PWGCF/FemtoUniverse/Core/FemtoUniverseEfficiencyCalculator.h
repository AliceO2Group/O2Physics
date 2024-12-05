#ifndef PWGCF_FEMTOUNIVERSE_CORE_EFFICIENCY_CALCULATOR_H_
#define PWGCF_FEMTOUNIVERSE_CORE_EFFICIENCY_CALCULATOR_H_

#include "Framework/Configurable.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/CallbackService.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/InitContext.h"
#include "CCDB/BasicCCDBManager.h"
#include "FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

namespace o2::analysis::femto_universe::efficiency
{
#define EFFICIENCY_CONFIGURABLES(name)   \
  struct : EfficiencyConfigurableGroup { \
  } name;

template <uint8_t T>
concept IsOneOrTwo = T == 1 || T == 2;

template <uint8_t N>
  requires IsOneOrTwo<N>
using MCTruthTrackHisto = FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, N>;

struct EfficiencyConfigurableGroup : ConfigurableGroup {
  Configurable<bool> shouldCalculate{"ConfEfficiencyCalculate", false, "Should calculate efficiency"};
  Configurable<bool> shouldApplyCorrections{"ConfEfficiencyApplyCorrections", false, "Should apply corrections from efficinecy"};
  Configurable<std::string> ccdbUrl{"ConfEfficiencyCCDBUrl", "http://alice-ccdb.cern.ch", "CCDB URL to be used"};
  Configurable<std::string> ccdbPath{"ConfEfficiencyCCDBPath", "Users/", "CCDB base path where objects should be put"};

  MCTruthTrackHisto<1> hMCTruth1;
  MCTruthTrackHisto<2> hMCTruth2;

  // TODO: move this to other struct?
  OutputObj<TH1F> hEff1{"Efficiency part1"};
  OutputObj<TH1F> hEff2{"Efficiency part2"};
};

class EfficiencyCalculator
{
 public:
  o2::ccdb::BasicCCDBManager& ccdb{o2::ccdb::BasicCCDBManager::instance()};
  EfficiencyConfigurableGroup config{};

  EfficiencyCalculator(EfficiencyConfigurableGroup& config)
    : shouldCalculate(config.shouldCalculate.value),
      shouldApplyCorrections(config.shouldApplyCorrections.value),
      hOutput({config.hEff1.object, config.hEff2.object}), // TODO: figure out better way
      hMCTruth1(&config.hMCTruth1),
      hMCTruth2(&config.hMCTruth2),
      ccdbUrl(config.ccdbUrl.value),
      ccdbPath(config.ccdbPath.value)
  {
  }

  auto init() -> void
  {
    if (!hOutput[0]) {
      LOG(fatal) << "!!hOutput[0];" << !!hOutput[0];
    }
    if (!hOutput[1]) {
      LOG(fatal) << "!!hOutput[1];" << !!hOutput[1];
    }

    ccdbApi.init(ccdbUrl);
    ccdb.setURL(ccdbUrl);
    ccdb.setLocalObjectValidityChecking();
    ccdb.setFatalWhenNull(false);

    // TODO: upload and load the actual hist for second particle
    hLoaded[0] = loadEfficiencyFromCCDB(1733060890131);
  }

  auto setIsTest(bool value) -> EfficiencyCalculator&
  {
    isTest = value;
    if (value) {
      ccdbUrl = "http://ccdb-test.cern.ch:8080";
    }
    return *this;
  }

  auto withRegistry(HistogramRegistry* reg) -> EfficiencyCalculator&
  {
    if (!reg) {
      LOGF(error, log("Registry value cannot be null"));
    }

    registry = reg;
    return *this;
  }

  auto withCCDBPath(const std::string& path) -> EfficiencyCalculator&
  {
    if (path.empty()) {
      LOGF(fatal, log("CCDB path cannot be empty"));
    }

    // TODO: add more distinction in path
    ccdbPath = std::format("{}/{}", path, folderPrefix);
    return *this;
  }

  template <uint8_t N>
    requires IsOneOrTwo<N>
  auto doMCTruth(auto particles) -> void
  {
    for (const auto& particle : particles) {
      if constexpr (N == 1) {
        hMCTruth1->fillQA<false, false>(particle);
      } else if constexpr (N == 2) {
        hMCTruth2->fillQA<false, false>(particle);
      }
    }
  }

  auto uploadOnStop(InitContext& ic) -> EfficiencyCalculator&
  {
    if (!shouldUploadOnStop) {
      shouldUploadOnStop = true;

      auto& callbacks = ic.services().get<CallbackService>();
      callbacks.set<o2::framework::CallbackService::Id::Stop>([this]() {
        for (uint8_t i{0}; i < hOutput.size(); i++) {
          const auto& output{hOutput[i]};
          if (isHistogramEmpty(output.get())) {
            LOGF(error, log("Histogram %d is empty - save aborted"), i + 1);
            return;
          }
          LOGF(debug, log("Found histogram %d: %s"), i + 1, output->GetTitle());

          std::map<string, string> metadata{};
          // metadata["runNumber"] = runNumber;

          long now{std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()};
          long oneYear{365LL * 24 * 60 * 60 * 1000};

          if (ccdbApi.storeAsTFileAny(output.get(), ccdbPath, metadata, now, now + oneYear) == 0) {
            LOGF(info, log("Histogram saved successfully"));
          } else {
            LOGF(fatal, log("Histogram save failed"));
          }
        }
      });
    } else {
      LOGF(warn, log("Save on stop is already set up"));
    }

    return *this;
  }

  template <uint8_t N>
    requires IsOneOrTwo<N>
  auto getWeight(auto const& particle) -> float
  {
    auto weight{1.0f};
    auto hEff{hLoaded[N - 1]};

    if (shouldApplyCorrections && hEff) {
      auto bin{hEff->FindBin(particle.pt())};
      auto eff{hEff->GetBinContent(bin)};
      weight /= eff > 0 ? eff : 1.0f;
    }
    return weight;
  }

  template <uint8_t N>
    requires IsOneOrTwo<N>
  auto calculate()
  {
    if (!shouldCalculate) {
      return;
    }

    std::shared_ptr<TH1> truth{registry->get<TH1>(HIST("MCTruthTracks_one/hPt"))};
    std::shared_ptr<TH1> reco{registry->get<TH1>(HIST("Tracks_one_MC/hPt"))};
    if (!truth || !reco) {
      LOGF(error, log("MC Truth & MC Reco histograms cannot be null"));
      return;
    }

    auto hEff{hOutput[N - 1]};
    if (!hEff) {
      LOGF(error, log("No OutputObj specified for particle %d histogram"), N);
      return;
    }

    for (int bin{0}; bin < hEff->GetNbinsX(); bin++) {
      auto denom{truth->GetBinContent(bin)};
      hEff->SetBinContent(bin, denom == 0 ? 0 : reco->GetBinContent(bin) / denom);
    }
  }

 private:
  auto log(const std::string& msg) -> const std::string
  {
    return std::format("[EFFICIENCY] {}", msg);
  }

  auto isHistogramEmpty(TH1* hist) -> bool
  {
    if (!hist) {
      return true;
    }

    // check overflow bins as well
    for (int idx{0}; idx <= hist->GetNbinsX() + 1; idx++) {
      if (hist->GetBinContent(idx) != 0) {
        return false;
      }
    }

    return true;
  }

  auto loadEfficiencyFromCCDB(long timestamp) -> TH1*
  {
    // TODO: ccdb.getSpecific<TH1>(ccdbPath, timestamp, metadata);
    auto hEff{ccdb.getForTimeStamp<TH1>(ccdbPath, timestamp)};
    if (!hEff || hEff->IsZombie()) {
      LOGF(fatal, log("Could not load histogram from %s"), ccdbPath);
      return nullptr;
    }

    LOGF(info, log("Histogram \"%s\" loaded from \"%s\""), hEff->GetTitle(), ccdbPath);
    return hEff;
  }

  bool isTest{false};

  bool shouldUploadOnStop{false};
  bool shouldCalculate{false};
  bool shouldApplyCorrections{false};

  HistogramRegistry* registry{};

  std::array<std::shared_ptr<TH1>, 2> hOutput{};
  std::array<TH1*, 2> hLoaded{};

  MCTruthTrackHisto<1>* hMCTruth1;
  MCTruthTrackHisto<2>* hMCTruth2;

  std::string ccdbUrl{};
  std::string ccdbPath{};
  o2::ccdb::CcdbApi ccdbApi;

  const std::string folderPrefix{"Efficiency"};
};

} // namespace o2::analysis::femto_universe::efficiency

#endif
