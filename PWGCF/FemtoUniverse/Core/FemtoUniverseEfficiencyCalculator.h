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

struct EfficiencyConfigurableGroup : ConfigurableGroup {
  Configurable<bool> confEfficiencyDoMCTruth{"confEfficiencyDoMCTruth", false, "Should fill MC Truth histogram"};
  Configurable<bool> confEfficiencyApplyCorrections{"confEfficiencyApplyCorrections", false, "Should apply corrections from efficiency"};
  Configurable<std::vector<std::string>> confEfficiencyCCDBLabels{"confEfficiencyCCDBLabels", {}, "Custom labels for efficiency objects in CCDB"};
  Configurable<int> confCCDBTrainNumber{"confCCDBTrainNumber", -1, "Train number for which to query CCDB objects (set to -1 to ignore)"};
  Configurable<std::vector<std::string>> confEfficiencyCCDBTimestamps{"confEfficiencyCCDBTimestamps", {"-1", "-1"}, "Timestamps in CCDB, to query for specific objects (default: -1 for both, the latest valid object)"};

  // NOTE: in the future we might move the below configurables to a separate struct, eg. CCDBConfigurableGroup
  Configurable<std::string> confCCDBUrl{"confCCDBUrl", "http://alice-ccdb.cern.ch", "CCDB URL to be used"};
  Configurable<std::string> confCCDBPath{"confCCDBPath", "", "CCDB base path to where to upload objects"};
};

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

    shouldDoTruth = config->confEfficiencyDoMCTruth;
    shouldApplyCorrections = config->confEfficiencyApplyCorrections;

    if (config->confEfficiencyApplyCorrections) {
      hLoaded = {
        loadEfficiencyFromCCDB(ParticleNo::ONE),
        loadEfficiencyFromCCDB(ParticleNo::TWO), //
      };
    }
  }

  template <size_t N>
    requires isOneOrTwo<N>
  auto doMCTruth(
    FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, N>& hMCTruth,
    const auto& particles) const -> void
  {
    if (shouldDoTruth) {
      for (const auto& particle : particles) {
        hMCTruth.template fillQA<false, false>(particle);
      }
    }
  }

  auto getWeight(const size_t partNo, const auto& particle) const -> float
  {
    auto weight = 1.0f;
    auto hEff = hLoaded[partNo - 1];

    if (shouldApplyCorrections && hEff) {
      auto bin = hEff->FindBin(particle.pt());
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

  auto loadEfficiencyFromCCDB(const size_t partNo) const -> TH1*
  {
    std::map<std::string, std::string> metadata{};

    if (partNo - 1 < config->confEfficiencyCCDBLabels->size()) {
      metadata["label"] = config->confEfficiencyCCDBLabels.value[partNo - 1];
    }
    if (config->confCCDBTrainNumber > 0) {
      metadata["trainNumber"] = std::to_string(config->confCCDBTrainNumber);
    }

    auto timestamp = partNo - 1 < config->confEfficiencyCCDBTimestamps->size()
                       ? std::stoll(config->confEfficiencyCCDBTimestamps.value[partNo - 1])
                       : -1;

    auto hEff = ccdb.getSpecific<TH1>(config->confCCDBPath, timestamp, metadata);
    if (!hEff || hEff->IsZombie()) {
      LOGF(error, notify("Could not load histogram from %s for particle %d"), config->confCCDBPath.value, partNo);
      return nullptr;
    }

    LOGF(info, notify("Histogram for particle %d loaded from \"%s\""), partNo, config->confCCDBPath.value);
    return hEff;
  }

  EfficiencyConfigurableGroup* config{};

  bool shouldDoTruth = false;
  bool shouldApplyCorrections = false;

  o2::ccdb::BasicCCDBManager& ccdb{o2::ccdb::BasicCCDBManager::instance()};
  std::array<TH1*, 2> hLoaded{};
};

} // namespace o2::analysis::femto_universe::efficiency

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEEFFICIENCYCALCULATOR_H_
