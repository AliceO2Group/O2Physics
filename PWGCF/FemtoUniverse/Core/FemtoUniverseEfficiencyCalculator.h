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
/// \brief Abstraction for applying corrections based on efficiency from CCDB
/// \author Dawid Karpiński, WUT Warsaw, dawid.karpinski@cern.ch

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEEFFICIENCYCALCULATOR_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEEFFICIENCYCALCULATOR_H_

#include <vector>
#include <map>
#include <string>

#include "Framework/Configurable.h"
#include "CCDB/BasicCCDBManager.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "FemtoUniverseParticleHisto.h"

namespace o2::analysis::femto_universe::efficiency
{
enum ParticleNo : size_t {
  ONE = 1,
  TWO
};

template <size_t T>
concept isOneOrTwo = T == ParticleNo::ONE || T == ParticleNo::TWO;

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
  Configurable<int> confEfficiencyCCDBTrainNumber{"confEfficiencyCCDBTrainNumber", -1, "Train number for which to query CCDB objects (set to -1 to ignore)"};
  Configurable<std::vector<std::string>> confEfficiencyCCDBTimestamps{"confEfficiencyCCDBTimestamps", {}, "Timestamps of efficiency histograms in CCDB, to query for specific objects (default: ['-1', '-1'], gets the latest valid objects for both)"};

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

    int64_t now = duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb.setCreatedNotAfter(now);

    shouldApplyCorrections = config->confEfficiencyApplyCorrections;

    if (config->confEfficiencyApplyCorrections && !config->confEfficiencyCCDBTimestamps.value.empty()) {
      for (const auto& timestamp : config->confEfficiencyCCDBTimestamps.value) {
        hLoaded.push_back(loadEfficiencyFromCCDB(std::stol(timestamp)));
      }

      LOGF(info, notify("Successfully loaded %d efficiency histogram(s)"), hLoaded.size());
    }
  }

  template <typename... BinVars>
    requires(sizeof...(BinVars) == getHistDim<HistType>())
  auto getWeight(ParticleNo partNo, const BinVars&... binVars) const -> float
  {
    auto weight = 1.0f;

    if (partNo - 1 < config->confEfficiencyCCDBTimestamps.value.size()) {
      auto hEff = hLoaded[partNo - 1];

      if (shouldApplyCorrections && hEff) {
        auto bin = hEff->FindBin(binVars...);
        auto eff = hEff->GetBinContent(bin);
        weight /= eff > 0 ? eff : 1.0f;
      }
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

    return hEff;
  }

  EfficiencyConfigurableGroup* config{};

  bool shouldApplyCorrections = false;

  o2::ccdb::BasicCCDBManager& ccdb{o2::ccdb::BasicCCDBManager::instance()};
  std::vector<HistType*> hLoaded{};
};

} // namespace o2::analysis::femto_universe::efficiency

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEEFFICIENCYCALCULATOR_H_
