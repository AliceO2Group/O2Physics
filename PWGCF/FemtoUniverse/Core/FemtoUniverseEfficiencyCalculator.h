// Copyright 2025 CERN and copyright holders of ALICE O2.
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
concept isOneOrTwo = T == 1 || T == 2;

struct EfficiencyConfigurableGroup : ConfigurableGroup {
  Configurable<bool> confEfficiencyCalculate{"confEfficiencyCalculate", false, "Should calculate efficiency"};
  Configurable<bool> confEfficiencyApplyCorrections{"confEfficiencyApplyCorrections", false, "Should apply corrections from efficiency"};

  Configurable<std::vector<std::string>> confCCDBLabels{"confCCDBLabels", std::vector<std::string>{"label1", "label2"}, "Labels for efficiency objects in CCDB"};

  OutputObj<TH1F> hEfficiency1{"Efficiency part1"};
  OutputObj<TH1F> hEfficiency2{"Efficiency part2"};

  // TODO: move to separate struct?
  Configurable<std::string> confCCDBUrl{"confCCDBUrl", "http://alice-ccdb.cern.ch", "CCDB URL to be used"};
  Configurable<std::string> confCCDBPath{"confCCDBPath", "", "CCDB base path to where to upload objects"};
  Configurable<long> confCCDBTimestamp{"confCCDBTimestamp", -1, "Timestamp from which to query CCDB objects"};

  // TODO: declare this in task directly?
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 1> hMCTruth1;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 2> hMCTruth2;
};

class EfficiencyCalculator
{
 public:
  o2::ccdb::BasicCCDBManager& ccdb{o2::ccdb::BasicCCDBManager::instance()};

  EfficiencyCalculator(EfficiencyConfigurableGroup* config) : config(config)
  {
  }

  auto init() -> void
  {
    ccdbApi.init(config->confCCDBUrl);
    ccdb.setURL(config->confCCDBUrl);
    ccdb.setLocalObjectValidityChecking();
    ccdb.setFatalWhenNull(false);

    shouldCalculate = config->confEfficiencyCalculate;
    shouldApplyCorrections = config->confEfficiencyApplyCorrections;

    ccdbFullPath = std::format("{}/{}", config->confCCDBPath.value, folderName);

    if (config->confEfficiencyCalculate) {
      hOutput = {config->hEfficiency1.object, config->hEfficiency2.object};
    }

    if (config->confEfficiencyApplyCorrections) {
      hLoaded = {
        loadEfficiencyFromCCDB<1>(config->confCCDBTimestamp),
        loadEfficiencyFromCCDB<2>(config->confCCDBTimestamp) //
      };
    }
  }

  template <uint8_t N>
    requires isOneOrTwo<N>
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

  auto uploadOnStop(InitContext& ic) -> void
  {
    if (!shouldCalculate) {
      return;
    }

    if (!shouldUploadOnStop) {
      shouldUploadOnStop = true;

      auto& callbacks = ic.services().get<CallbackService>();
      callbacks.set<o2::framework::CallbackService::Id::Stop>([this]() {
        for (auto i = 0UL; i < hOutput.size(); i++) {
          const auto& output = hOutput[i];

          if (isHistogramEmpty(output.get())) {
            LOGF(error, notify("Histogram %d is empty - save aborted"), i + 1);
            return;
          }
          LOGF(debug, notify("Found histogram %d: %s"), i + 1, output->GetTitle());

          long now = duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
          long oneYear = 365LL * 24 * 60 * 60 * 1000;

          if (ccdbApi.storeAsTFileAny(output.get(), ccdbFullPath, createMetadata(i), now, now + oneYear) == 0) {
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
  auto calculate(std::shared_ptr<TH1> truth, std::shared_ptr<TH1> reco) const
  {
    if (!shouldCalculate) {
      return;
    }

    if (!truth || !reco) {
      LOGF(error, notify("MC Truth & MC Reco histograms cannot be null"));
      return;
    }

    auto hEff = hOutput[N - 1];
    if (!hEff) {
      LOGF(error, notify("No OutputObj specified for particle %d histogram"), N);
      return;
    }

    for (auto bin = 0; bin < hEff->GetNbinsX(); bin++) {
      auto denom = truth->GetBinContent(bin);
      hEff->SetBinContent(bin, denom == 0 ? 0 : reco->GetBinContent(bin) / denom);
    }
  }

 private:
  static inline auto notify(const std::string& msg) -> const std::string
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
    if (config->confCCDBLabels->size() != 2) {
      LOGF(fatal, notify("CCDB labels configurable should be exactly of size 2"));
    }
    return std::map<std::string, std::string>{
      {"label", config->confCCDBLabels.value[partNo]} //
    };
  }

  template <uint8_t N>
    requires isOneOrTwo<N>
  auto loadEfficiencyFromCCDB(long timestamp) const -> TH1*
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

  static constexpr std::string folderName{"Efficiency"};
  static constexpr std::array<std::string, 3> histSuffix{"", "_one", "_two"};
};

} // namespace o2::analysis::femto_universe::efficiency

#endif
