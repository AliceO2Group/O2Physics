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

template <uint8_t T>
concept IsOneOrTwo = T == 1 || T == 2;

namespace ParticleNo
{
enum ParticleNo {
  ONE,
  TWO
};
}

class EfficiencyCalculator
{
 public:
  o2::ccdb::BasicCCDBManager& ccdb{o2::ccdb::BasicCCDBManager::instance()};

  ConfigurableAxis ConfTempFitVarpTBins{"ConfTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTempFitVarPDGBins{"ConfTempFitVarPDGBins", {6000, -2300, 2300}, "Binning of the PDG code in the pT vs. TempFitVar plot"};

  auto init() -> void
  {
    assert(registry && "Registry has to be set");
    assert(ccdbUrl && "CCDB URL has to be set");

    hMCTruth1.init(registry, ConfTempFitVarpTBins, ConfTempFitVarPDGBins, false, pdgCode[ParticleNo::ONE], false);
    hMCTruth2.init(registry, ConfTempFitVarpTBins, ConfTempFitVarPDGBins, false, pdgCode[ParticleNo::TWO], false);

    ccdbApi.init(ccdbUrl);
    ccdb.setURL(ccdbUrl);
    ccdb.setLocalObjectValidityChecking();
    ccdb.setFatalWhenNull(false);

    hLoaded[ParticleNo::ONE] = loadEfficiencyFromCCDB(1733060890131);
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
    registry = reg;
    return *this;
  }

  template <uint8_t N>
    requires IsOneOrTwo<N>
  auto setParticle(int pdg, TH1F* outputHist = nullptr) -> EfficiencyCalculator&
  {
    pdgCode[N - 1] = pdg;
    hOutput[N - 1] = outputHist;
    return *this;
  }

  auto withCCDBPath(const std::string& path) -> EfficiencyCalculator&
  {
    ccdbPath = path + folderPrefix;
    return *this;
  }

  template <uint8_t N>
    requires IsOneOrTwo<N>
  auto doMCTruth(auto particles) -> void
  {
    for (const auto& particle : particles) {
      if (static_cast<int>(particle.pidcut()) != pdgCode[N - 1]) {
        continue;
      }

      if constexpr (N == 1) {
        hMCTruth1.fillQA<false, false>(particle);
      } else if constexpr (N == 2) {
        hMCTruth2.fillQA<false, false>(particle);
      }
    }
  }

  auto saveOnStop(InitContext& ic) -> void
  {
    if (!shouldSaveOnStop) {
      shouldSaveOnStop = true;

      auto& callbacks = ic.services().get<CallbackService>();
      callbacks.set<o2::framework::CallbackService::Id::Stop>([this]() {
        for (uint8_t i{0}; i < hOutput.size(); i++) {
          const auto& output{hOutput[i]};
          if (isHistogramEmpty(output)) {
            LOGF(error, log("Histogram %s is empty - save aborted"), i + 1);
            return;
          }
          LOGF(debug, log("Found histogram %d: %s"), i + 1, output->GetTitle());

          std::map<string, string> metadata{};
          // metadata["runNumber"] = runNumber;

          long now{std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()};
          long oneYear{365LL * 24 * 60 * 60 * 1000};

          if (ccdbApi.storeAsTFileAny(output, ccdbPath, metadata, now, now + oneYear) == 0) {
            LOGF(info, log("Histogram saved successfully"));
          } else {
            LOGF(fatal, log("Histogram save failed"));
          }
        }
      });
    } else {
      LOGF(warn, log("Save on stop is already set up"));
    }
  }

  template <uint8_t N>
    requires IsOneOrTwo<N>
  auto getWeight(auto const& particle) -> float
  {
    auto weight{1.0f};
    if (auto hEff{hLoaded[N - 1]}) {
      auto bin{hEff->FindBin(particle.pt())};
      auto eff{hEff->GetBinContent(bin)};
      weight /= eff > 0 ? eff : 1.0f;
    } else {
      LOGF(error, log("Called `getWeight` with empty histogram from CCDB"));
    }
    return weight;
  }

  template <uint8_t N>
    requires IsOneOrTwo<N>
  auto calculate()
  {
    std::shared_ptr<TH1> truth{registry->get<TH1>(HIST("MCTruthTracks_one/hPt"))};
    std::shared_ptr<TH1> reco{registry->get<TH1>(HIST("Tracks_one_MC/hPt"))};

    auto hEff{hOutput[N - 1]};
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
    auto hEff{ccdb.getForTimeStamp<TH1>(ccdbPath, timestamp)};
    if (!hEff || hEff->IsZombie()) {
      LOGF(fatal, log("Could not load histogram from %s"), ccdbPath);
      return nullptr;
    }

    LOGF(info, log("Histogram \"%s\" loaded from \"%s\""), hEff->GetTitle(), ccdbPath);
    return hEff;
  }

  bool isTest{false};

  HistogramRegistry* registry{};

  std::array<int, 2> pdgCode;
  std::array<TH1*, 2> hOutput{};
  std::array<TH1*, 2> hLoaded{};

  std::string ccdbUrl{"http://alice-ccdb.cern.ch"};
  std::string ccdbPath{};

  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 1> hMCTruth1;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 2> hMCTruth2;

  const std::string folderPrefix{"Efficiency"};
  o2::ccdb::CcdbApi ccdbApi;

  bool shouldSaveOnStop{false};
};

} // namespace o2::analysis::femto_universe

#endif
