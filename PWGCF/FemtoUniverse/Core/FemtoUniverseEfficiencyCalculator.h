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
template <uint8_t T>
concept IsOneOrTwo = T == 1 || T == 2;

template <uint8_t N>
  requires IsOneOrTwo<N>
using MCTruthTrackHisto = FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, N>;

struct EfficiencyConfigurableGroup : ConfigurableGroup {
  Configurable<bool> shouldCalculate{"ConfEfficiencyCalculate", false, "Should calculate efficiency"};
  Configurable<bool> shouldApplyCorrections{"ConfEfficiencyApplyCorrections", false, "Should apply corrections from efficinecy"};
  Configurable<std::string> ccdbUrl{"ConfEfficiencyCCDBUrl", "http://alice-ccdb.cern.ch", "CCDB URL to be used"};
  Configurable<std::string> ccdbPath{"ConfEfficiencyCCDBPath", "", "CCDB base path where objects should be put"};

  MCTruthTrackHisto<1> hMCTruth1;
  MCTruthTrackHisto<2> hMCTruth2;

  OutputObj<TH1F> hEff1{"Efficiency part1"};
  OutputObj<TH1F> hEff2{"Efficiency part2"};
};

class EfficiencyCalculator
{
 public:
  o2::ccdb::BasicCCDBManager& ccdb{o2::ccdb::BasicCCDBManager::instance()};
  EfficiencyConfigurableGroup* config{};

  explicit EfficiencyCalculator(EfficiencyConfigurableGroup* config)
    : config(config),
      shouldCalculate(config->shouldCalculate.value),
      shouldApplyCorrections(config->shouldApplyCorrections.value)
  {
  }

  auto init() -> void
  {
    ccdbApi.init(config->ccdbUrl);
    ccdb.setURL(config->ccdbUrl);
    ccdb.setLocalObjectValidityChecking();
    ccdb.setFatalWhenNull(false);

    hOutput = {config->hEff1.object, config->hEff2.object};

    // TODO: upload and load the actual hist for second particle
    // hLoaded = {loadEfficiencyFromCCDB(1733060890131), loadEfficiencyFromCCDB(1733060890131)};
    hLoaded[0] = loadEfficiencyFromCCDB(1733060890131);
  }

  auto withRegistry(HistogramRegistry* reg) -> EfficiencyCalculator&
  {
    if (!reg) {
      LOGF(error, log("Registry value cannot be null"));
    }

    registry = reg;
    return *this;
  }

  template <uint8_t N>
    requires IsOneOrTwo<N>
  auto doMCTruth(auto particles) -> void
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

          if (ccdbApi.storeAsTFileAny(output.get(), config->ccdbPath, metadata, now, now + oneYear) == 0) {
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

  std::string path(const std::string& base, const std::initializer_list<std::string>& segments)
  {
    std::string url = base;
    if (!url.empty() && url.back() == '/') {
      url.pop_back();
    }

    for (const auto& segment : segments) {
      if (!segment.empty()) {
        if (url.back() != '/') {
          url += '/';
        }
        url += segment.front() == '/' ? segment.substr(1) : segment;
      }
    }

    return url;
  }

  auto loadEfficiencyFromCCDB(long timestamp) -> TH1*
  {
    // TODO: ccdb.getSpecific<TH1>(ccdbPath, timestamp, metadata);
    auto fullPath{std::format("{}/{}", config->ccdbPath.value, folderPrefix)};

    auto hEff{ccdb.getForTimeStamp<TH1>(fullPath, timestamp)};
    if (!hEff || hEff->IsZombie()) {
      LOGF(fatal, log("Could not load histogram from %s"), config->ccdbPath.value);
      return nullptr;
    }

    LOGF(info, log("Histogram \"%s\" loaded from \"%s\""), hEff->GetTitle(), config->ccdbPath.value);
    return hEff;
  }

  bool shouldUploadOnStop{false};
  bool shouldCalculate{false};
  bool shouldApplyCorrections{false};

  HistogramRegistry* registry{};

  std::array<std::shared_ptr<TH1>, 2> hOutput{};
  std::array<TH1*, 2> hLoaded{};

  o2::ccdb::CcdbApi ccdbApi;

  const std::string folderPrefix{"Efficiency"};
};

} // namespace o2::analysis::femto_universe::efficiency

#endif
