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

/// \file FemtoUniverseEfficiencyCorrection.h
/// \brief Abstraction for applying efficiency corrections based on weights from CCDB
/// \author Dawid Karpiński, WUT Warsaw, dawid.karpinski@cern.ch

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEEFFICIENCYCORRECTION_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEEFFICIENCYCORRECTION_H_

#include <vector>
#include <string>
#include <algorithm>

#include "Framework/Configurable.h"
#include "CCDB/BasicCCDBManager.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

namespace o2::analysis::femto_universe::efficiency_correction
{
enum ParticleNo : size_t {
  ONE = 1,
  TWO,
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

struct EffCorConfigurableGroup : framework::ConfigurableGroup {
  framework::Configurable<bool> confEffCorApply{"confEffCorApply", false, "[Efficiency Correction] Should apply efficiency corrections"};
  framework::Configurable<std::string> confEffCorCCDBUrl{"confEffCorCCDBUrl", "http://alice-ccdb.cern.ch", "[Efficiency Correction] CCDB URL to use"};
  framework::Configurable<std::string> confEffCorCCDBPath{"confEffCorCCDBPath", "", "[Efficiency Correction] CCDB path to histograms"};
  framework::Configurable<std::vector<std::string>> confEffCorCCDBTimestamps{"confEffCorCCDBTimestamps", {}, "[Efficiency Correction] Timestamps of histograms in CCDB (0 can be used as a placeholder, e.g. when running subwagons)"};
};

template <typename HistType>
  requires std::is_base_of_v<TH1, HistType>
class EfficiencyCorrection
{
 public:
  explicit EfficiencyCorrection(EffCorConfigurableGroup* config) : config(config) // o2-linter: disable=name/function-variable
  {
  }

  auto init() -> void
  {
    ccdb.setURL(config->confEffCorCCDBUrl);
    ccdb.setLocalObjectValidityChecking();
    ccdb.setFatalWhenNull(false);

    auto now = duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb.setCreatedNotAfter(now);

    shouldApplyCorrection = config->confEffCorApply;

    if (shouldApplyCorrection && !config->confEffCorCCDBTimestamps.value.empty()) {
      for (auto idx = 0UL; idx < config->confEffCorCCDBTimestamps.value.size(); idx++) {
        auto timestamp = 0L;
        try {
          timestamp = std::max(0L, std::stol(config->confEffCorCCDBTimestamps.value[idx]));
        } catch (const std::exception&) {
          LOGF(error, notify("Could not parse CCDB timestamp \"%s\""), config->confEffCorCCDBTimestamps.value[idx]);
          continue;
        }

        hLoaded[idx] = timestamp > 0 ? loadHistFromCCDB(timestamp) : nullptr;
      }
    }
  }

  template <typename... BinVars>
    requires(sizeof...(BinVars) == getHistDim<HistType>())
  auto getWeight(ParticleNo partNo, const BinVars&... binVars) const -> float
  {
    auto weight = 1.0f;
    auto hWeights = hLoaded[partNo - 1];

    if (shouldApplyCorrection && hWeights) {
      auto bin = hWeights->FindBin(binVars...);
      weight = hWeights->GetBinContent(bin);
    }

    return weight;
  }

 private:
  static inline auto notify(const std::string& msg) -> const std::string
  {
    return fmt::format("[EFFICIENCY CORRECTION] {}", msg);
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

  auto loadHistFromCCDB(const int64_t timestamp) const -> HistType*
  {
    auto hWeights = ccdb.getForTimeStamp<HistType>(config->confEffCorCCDBPath, timestamp);
    if (!hWeights || hWeights->IsZombie()) {
      LOGF(error, notify("Could not load histogram \"%s/%ld\""), config->confEffCorCCDBPath.value, timestamp);
      return nullptr;
    }

    if (isHistEmpty(hWeights)) {
      LOGF(warn, notify("Histogram \"%s/%ld\" has been loaded, but it is empty"), config->confEffCorCCDBUrl.value, timestamp);
    }

    LOGF(info, notify("Successfully loaded %ld"), timestamp);
    return hWeights;
  }

  EffCorConfigurableGroup* config{};

  bool shouldApplyCorrection = false;

  o2::ccdb::BasicCCDBManager& ccdb{o2::ccdb::BasicCCDBManager::instance()};
  std::array<HistType*, 2> hLoaded{};
};

} // namespace o2::analysis::femto_universe::efficiency_correction

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEEFFICIENCYCORRECTION_H_
