// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
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
/// \brief Abstraction for applying corrections based on efficiency from CCDB
/// \author Dawid Karpiński, WUT Warsaw, dawid.karpinski@cern.ch

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEEFFICIENCYCALCULATOR_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEEFFICIENCYCALCULATOR_H_

#include "FemtoUniverseParticleHisto.h"

#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/Configurable.h"

#include <TH1.h>

#include <algorithm>
#include <map>
#include <string>
#include <vector>

namespace o2::analysis::femto_universe::efficiency
{
enum ParticleNo : size_t {
  ONE = 1,
  TWO,
};

template <size_t T>
concept IsOneOrTwo = T == ParticleNo::ONE || T == ParticleNo::TWO;

template <typename T>
consteval auto getHistDim() -> int
{
  if (std::is_same_v<T, TH1>)
    return 1;
  else if (std::is_same_v<T, TH2>)
    return 2;
  else if (std::is_same_v<T, TH3>)
    return 3;
  else
    return -1;
}

struct EfficiencyConfigurableGroup : ConfigurableGroup {
  Configurable<bool> confEfficiencyApplyCorrections{"confEfficiencyApplyCorrections", false, "Should apply corrections from efficiency"};
  Configurable<int> confEfficiencyCCDBTrainNumber{"confEfficiencyCCDBTrainNumber", 0, "Train number for which to query CCDB objects (set 0 to ignore)"};
  Configurable<std::vector<std::string>> confEfficiencyCCDBTimestamps{"confEfficiencyCCDBTimestamps", {}, "Timestamps of efficiency histograms in CCDB, to query for specific objects (default: [], set 0 to ignore, useful when running subwagons)"};

  // NOTE: in the future we might move the below configurables to a separate struct, eg. CCDBConfigurableGroup
  Configurable<std::string> confCCDBUrl{"confCCDBUrl", "http://alice-ccdb.cern.ch", "CCDB URL to be used"};
  Configurable<std::string> confCCDBPath{"confCCDBPath", "", "CCDB base path to where to upload objects"};
};

template <typename HistType>
  requires std::is_base_of_v<TH1, HistType>
class EfficiencyCalculator
{
 public:
  explicit EfficiencyCalculator(EfficiencyConfigurableGroup* config) : config(config) // o2-linter: disable=name/function-variable
  {
  }

  auto init() -> void
  {
    ccdb.setURL(config->confCCDBUrl);
    ccdb.setLocalObjectValidityChecking();
    ccdb.setFatalWhenNull(false);

    auto now = duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb.setCreatedNotAfter(now);

    shouldApplyCorrections = config->confEfficiencyApplyCorrections;

    if (config->confEfficiencyApplyCorrections && !config->confEfficiencyCCDBTimestamps.value.empty()) {
      for (auto idx = 0UL; idx < config->confEfficiencyCCDBTimestamps.value.size(); idx++) {
        auto timestamp = 0L;
        try {
          timestamp = std::max(0L, std::stol(config->confEfficiencyCCDBTimestamps.value[idx]));
        } catch (const std::exception&) {
          LOGF(error, notify("Could not parse CCDB timestamp \"%s\""), config->confEfficiencyCCDBTimestamps.value[idx]);
          continue;
        }

        hLoaded[idx] = timestamp > 0 ? loadEfficiencyFromCCDB(timestamp) : nullptr;
      }
    }
  }

  template <typename... BinVars>
    requires(sizeof...(BinVars) == getHistDim<HistType>())
  auto getWeight(ParticleNo partNo, const BinVars&... binVars) const -> float
  {
    auto weight = 1.0f;
    auto hEff = hLoaded[partNo - 1];

    if (shouldApplyCorrections && hEff) {
      auto bin = hEff->FindBin(static_cast<double>(binVars)...);
      auto eff = hEff->GetBinContent(bin);
      weight /= eff > 0 ? eff : 1.0f;
    }

    return weight;
  }

 private:
  static inline auto notify(const std::string& msg) -> const std::string
  {
    return fmt::format("[EFFICIENCY] {}", msg);
  }

  static auto isHistEmpty(HistType* hist) -> bool
  {
    if (!hist) {
      return true;
    }
    for (auto idx = 0; idx <= hist->GetNbinsX() + 1; idx++) {
      if (hist->GetBinContent(idx) > 0) {
        return false;
      }
    }
    return true;
  }

  auto loadEfficiencyFromCCDB(const int64_t timestamp) const -> HistType*
  {
    std::map<std::string, std::string> metadata{};

    if (config->confEfficiencyCCDBTrainNumber > 0) {
      metadata["trainNumber"] = std::to_string(config->confEfficiencyCCDBTrainNumber);
    }

    auto hEff = ccdb.getSpecific<HistType>(config->confCCDBPath, timestamp, metadata);
    if (!hEff || hEff->IsZombie()) {
      LOGF(error, notify("Could not load histogram \"%s/%ld\""), config->confCCDBPath.value, timestamp);
      return nullptr;
    }

    if (isHistEmpty(hEff)) {
      LOGF(warn, notify("Histogram \"%s/%ld\" has been loaded, but it is empty"), config->confCCDBPath.value, timestamp);
    }

    auto clonedEffHist = static_cast<HistType*>(hEff->Clone());
    clonedEffHist->SetDirectory(nullptr);

    LOGF(info, notify("Successfully loaded %ld"), timestamp);
    return clonedEffHist;
  }

  EfficiencyConfigurableGroup* config{};

  bool shouldApplyCorrections = false;

  o2::ccdb::BasicCCDBManager& ccdb{o2::ccdb::BasicCCDBManager::instance()};
  std::array<HistType*, 2> hLoaded{nullptr, nullptr};
};

} // namespace o2::analysis::femto_universe::efficiency

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEEFFICIENCYCALCULATOR_H_
