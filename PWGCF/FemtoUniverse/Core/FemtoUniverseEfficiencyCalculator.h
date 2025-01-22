// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoUniverseEfficiencyCalculator.h
/// \brief Abstraction for calculating efficiency and applying corrections with the help of CCDB
/// \author Dawid Karpiński, WUT Warsaw, dawid.karpinski@cern.ch

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEEFFICIENCYCALCULATOR_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEEFFICIENCYCALCULATOR_H_

#include <vector>
#include <memory>
#include <map>
#include <string>

#include "Framework/Configurable.h"
#include "Framework/CallbackService.h"
#include "Framework/InitContext.h"
#include "CCDB/BasicCCDBManager.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "FemtoUniverseParticleHisto.h"

namespace o2::analysis::femto_universe::efficiency
{
template <uint8_t T>
concept isOneOrTwo = T == 1 || T == 2;

struct EfficiencyConfigurableGroup : ConfigurableGroup {
  Configurable<bool> confEfficiencyCalculate{"confEfficiencyCalculate", false, "Should calculate efficiency"};
  Configurable<bool> confEfficiencyApplyCorrections{"confEfficiencyApplyCorrections", false, "Should apply corrections from efficiency"};
  Configurable<std::vector<std::string>> confEfficiencyCCDBLabels{"confEfficiencyCCDBLabels", {}, "Labels for efficiency objects in CCDB"};
  ConfigurableAxis confEfficiencyCCDBTimestamps{"confEfficiencyCCDBTimestamps", {-1, -1}, "Timestamps from which to query CCDB objects (default: -1 for both)"};

  // NOTE: in the future we might move the below configurables to a separate struct, eg. CCDBConfigurableGroup
  Configurable<std::string> confCCDBUrl{"confCCDBUrl", "http://alice-ccdb.cern.ch", "CCDB URL to be used"};
  Configurable<std::string> confCCDBPath{"confCCDBPath", "", "CCDB base path to where to upload objects"};
  Configurable<int64_t> confCCDBLifetime{"confCCDBLifetime", 365LL * 24 * 60 * 60 * 1000, "Lifetime of uploaded objects (default: 1 year)"};
};

class EfficiencyCalculator
{
 public:
  o2::ccdb::BasicCCDBManager& ccdb{o2::ccdb::BasicCCDBManager::instance()};

  explicit EfficiencyCalculator(EfficiencyConfigurableGroup* config) : config(config) // o2-linter: disable=name/function-variable
  {
  }

  auto init() -> void
  {
    ccdbApi.init(config->confCCDBUrl);
    ccdb.setURL(config->confCCDBUrl);
    ccdb.setLocalObjectValidityChecking();
    ccdb.setFatalWhenNull(false);

    int64_t now = duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb.setCreatedNotAfter(now);

    shouldCalculate = config->confEfficiencyCalculate;
    shouldApplyCorrections = config->confEfficiencyApplyCorrections;

    ccdbFullPath = fmt::format("{}/{}", config->confCCDBPath.value, folderName);

    if (config->confEfficiencyApplyCorrections) {
      if (config->confEfficiencyCCDBTimestamps->size() != 2) {
        LOGF(fatal, notify("CCDB timestamps configurable should be exactly of size 2"));
      }
      hLoaded = {
        loadEfficiencyFromCCDB<1>(config->confEfficiencyCCDBTimestamps.value[0]),
        loadEfficiencyFromCCDB<2>(config->confEfficiencyCCDBTimestamps.value[1]) //
      };
    }
  }

  template <uint8_t N>
    requires isOneOrTwo<N>
  auto doMCTruth(FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, N> hMCTruth, auto particles) const -> void
  {
    for (const auto& particle : particles) {
      hMCTruth.template fillQA<false, false>(particle);
    }
  }

  auto uploadOnStop(InitContext& ic) -> void
  {
    if (!shouldCalculate) {
      return;
    }

    if (!shouldUploadOnStop) {
      shouldUploadOnStop = true;

      auto& callbacks = ic.services().get<CallbackService>();
      callbacks.set<o2::framework::CallbackService::Id::EndOfStream>([this](EndOfStreamContext&) {
        for (auto i = 0UL; i < hOutput.size(); i++) {
          const auto& output = hOutput[i];

          if (isHistogramEmpty(output.get())) {
            LOGF(error, notify("Histogram %d is empty - save aborted"), i + 1);
            return;
          }
          LOGF(debug, notify("Found histogram %d: %s"), i + 1, output->GetTitle());

          int64_t now = duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

          if (ccdbApi.storeAsTFileAny(output.get(), ccdbFullPath, createMetadata(i), now, now + config->confCCDBLifetime) == 0) {
            LOGF(info, notify("Histogram %d saved successfully"), i + 1);
          } else {
            LOGF(fatal, notify("Histogram %d save failed"), i + 1);
          }
        }
      });
    } else {
      LOGF(warn, notify("Uploading on stop callback is already set up"));
    }
  }

  template <uint8_t N>
    requires isOneOrTwo<N>
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
    requires isOneOrTwo<N>
  auto calculate(std::shared_ptr<TH1> hEff, std::shared_ptr<TH1> hTruth, std::shared_ptr<TH1> hReco)
  {
    if (!shouldCalculate) {
      return;
    }

    if (!hTruth || !hReco) {
      LOGF(error, notify("MC Truth & MC Reco histograms cannot be null"));
      return;
    }

    if (!hEff) {
      LOGF(error, notify("No target histogram to fill specified for particle %d"), N);
      return;
    }

    for (auto bin = 0; bin < hEff->GetNbinsX(); bin++) {
      auto denom = hTruth->GetBinContent(bin);
      hEff->SetBinContent(bin, denom == 0 ? 0 : hReco->GetBinContent(bin) / denom);
    }

    hOutput[N - 1] = hEff;
  }

 private:
  static inline auto notify(const std::string& msg) -> const std::string
  {
    return fmt::format("[EFFICIENCY] {}", msg);
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
    if (config->confEfficiencyCCDBLabels->size() != 2) {
      LOGF(fatal, notify("CCDB labels configurable should be exactly of size 2"));
    }
    return std::map<std::string, std::string>{
      {"label", config->confEfficiencyCCDBLabels.value[partNo]} //
    };
  }

  template <uint8_t N>
    requires isOneOrTwo<N>
  auto loadEfficiencyFromCCDB(int64_t timestamp) const -> TH1*
  {
    auto hEff = ccdb.getSpecific<TH1>(ccdbFullPath, timestamp, createMetadata(N - 1));
    if (!hEff || hEff->IsZombie()) {
      LOGF(error, notify("Could not load histogram from %s"), config->confCCDBPath.value);
      return nullptr;
    }

    LOGF(info, notify("Histogram \"%s\" loaded from \"%s\""), hEff->GetTitle(), config->confCCDBPath.value);
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

  static constexpr std::string_view folderName{"Efficiency"};
  static constexpr std::array<std::string_view, 3> histSuffix{"", "_one", "_two"};
};

} // namespace o2::analysis::femto_universe::efficiency

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEEFFICIENCYCALCULATOR_H_
