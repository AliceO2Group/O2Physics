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

#include <CCDB/BasicCCDBManager.h>
#include <Framework/Configurable.h>
#include <Framework/ConfigurableKinds.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>

#include <Framework/O2DatabasePDGPlugin.h>
#include <TH1.h>
#include <TH2.h>

#include <vector>
#include <string>
#include <algorithm>
#include <ranges>

#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

namespace o2::analysis::femto_universe::efficiency_correction
{
enum ParticleNo : size_t {
  ONE = 1,
  TWO,
};

template <size_t T>
concept IsOneOrTwo = T == ParticleNo::ONE || T == ParticleNo::TWO;

struct EffCorConfigurableGroup : framework::ConfigurableGroup {
  framework::Configurable<bool> confEffCorApply{"confEffCorApply", false, "[Efficiency Correction] Should apply efficiency corrections"};
  framework::Configurable<bool> confEffCorFillHist{"confEffCorFillHist", false, "[Efficiency Correction] Should fill histograms for efficiency corrections"};
  framework::Configurable<std::string> confEffCorCCDBUrl{"confEffCorCCDBUrl", "http://alice-ccdb.cern.ch", "[Efficiency Correction] CCDB URL to use"};
  framework::Configurable<std::string> confEffCorCCDBPath{"confEffCorCCDBPath", "", "[Efficiency Correction] CCDB path to histograms"};
  framework::Configurable<std::vector<std::string>> confEffCorCCDBTimestamps{"confEffCorCCDBTimestamps", {}, "[Efficiency Correction] Timestamps of histograms in CCDB (0 can be used as a placeholder, e.g. when running subwagons)"};
  framework::Configurable<std::string> confEffCorVariables{"confEffCorVariables", "pt", "[Efficiency Correction] Variables for efficiency correction histogram dimensions (available: 'pt'; 'pt,eta'; 'pt,mult'; 'pt,eta,mult')"};
};

class EfficiencyCorrection
{
 public:
  explicit EfficiencyCorrection(EffCorConfigurableGroup* config) : config(config) // o2-linter: disable=name/function-variable
  {
  }

  auto init(framework::HistogramRegistry* registry, std::vector<framework::AxisSpec> axisSpecs) -> void
  {
    shouldFillHistograms = config->confEffCorFillHist;

    histRegistry = registry;
    if (shouldFillHistograms) {
      for (const auto& suffix : histSuffix) {
        auto path = std::format("{}/{}", histDirectory, suffix);
        registry->add((path + "/hMCTruth").c_str(), "MCTruth; #it{p}_{T} (GeV/#it{c}); #it{#eta}; Mult", framework::kTH3F, axisSpecs);
        registry->add((path + "/hPrimary").c_str(), "Primary; #it{p}_{T} (GeV/#it{c}); #it{#eta}; Mult", framework::kTH3F, axisSpecs);
        registry->add((path + "/hSecondary").c_str(), "Secondary; #it{p}_{T} (GeV/#it{c}); #it{#eta}; Mult", framework::kTH3F, axisSpecs);
        registry->add((path + "/hMaterial").c_str(), "Material; #it{p}_{T} (GeV/#it{c}); #it{#eta}; Mult", framework::kTH3F, axisSpecs);
        registry->add((path + "/hFake").c_str(), "Fake; #it{p}_{T} (GeV/#it{c}); #it{#eta}; Mult", framework::kTH3F, axisSpecs);
        registry->add((path + "/hOther").c_str(), "Other; #it{p}_{T} (GeV/#it{c}); #it{#eta}; Mult", framework::kTH3F, axisSpecs);
      }
    }

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

        if (timestamp > 0) {
          switch (getDimensionFromVariables()) {
            case 1:
              hLoaded[idx] = loadHistFromCCDB<TH1>(timestamp);
              break;
            case 2:
              hLoaded[idx] = loadHistFromCCDB<TH2>(timestamp);
              break;
            case 3:
              hLoaded[idx] = loadHistFromCCDB<TH3>(timestamp);
              break;
            default:
              LOGF(fatal, notify("Unknown configuration for efficiency variables"));
              break;
          }
        } else {
          hLoaded[idx] = nullptr;
        }
      }
    }
  }

  template <uint8_t N>
    requires IsOneOrTwo<N>
  void fillTruthHist(auto particle)
  {
    if (!shouldFillHistograms) {
      return;
    }

    histRegistry->fill(HIST(histDirectory) + HIST("/") + HIST(histSuffix[N - 1]) + HIST("/hMCTruth"),
                       particle.pt(),
                       particle.eta(),
                       particle.fdCollision().multV0M());
  }

  template <uint8_t N>
    requires IsOneOrTwo<N>
  void fillRecoHist(auto particle, int particlePDG)
  {
    if (!shouldFillHistograms) {
      return;
    }

    if (!particle.has_fdMCParticle()) {
      return;
    }

    auto mcParticle = particle.fdMCParticle();

    if (mcParticle.pdgMCTruth() == particlePDG) {
      switch (mcParticle.partOriginMCTruth()) {
        case (o2::aod::femtouniverse_mc_particle::kPrimary):
          histRegistry->fill(HIST(histDirectory) + HIST("/") + HIST(histSuffix[N - 1]) + HIST("/hPrimary"),
                             mcParticle.pt(),
                             mcParticle.eta(),
                             particle.fdCollision().multV0M());
          break;

        case (o2::aod::femtouniverse_mc_particle::kDaughter):
        case (o2::aod::femtouniverse_mc_particle::kDaughterLambda):
        case (o2::aod::femtouniverse_mc_particle::kDaughterSigmaplus):
          histRegistry->fill(HIST(histDirectory) + HIST("/") + HIST(histSuffix[N - 1]) + HIST("/hSecondary"),
                             mcParticle.pt(),
                             mcParticle.eta(),
                             particle.fdCollision().multV0M());
          break;

        case (o2::aod::femtouniverse_mc_particle::kMaterial):
          histRegistry->fill(HIST(histDirectory) + HIST("/") + HIST(histSuffix[N - 1]) + HIST("/hMaterial"),
                             mcParticle.pt(),
                             mcParticle.eta(),
                             particle.fdCollision().multV0M());
          break;

        case (o2::aod::femtouniverse_mc_particle::kFake):
          histRegistry->fill(HIST(histDirectory) + HIST("/") + HIST(histSuffix[N - 1]) + HIST("/hFake"),
                             mcParticle.pt(),
                             mcParticle.eta(),
                             particle.fdCollision().multV0M());
          break;

        default:
          histRegistry->fill(HIST(histDirectory) + HIST("/") + HIST(histSuffix[N - 1]) + HIST("/hOther"),
                             mcParticle.pt(),
                             mcParticle.eta(),
                             particle.fdCollision().multV0M());
          break;
      }
    }
  }

  auto getWeight(ParticleNo partNo, auto particle) -> float
  {
    auto weight = 1.0f;
    auto hWeights = hLoaded[partNo - 1];

    if (shouldApplyCorrection && hWeights) {
      auto dim = static_cast<size_t>(hWeights->GetDimension());
      if (dim != getDimensionFromVariables()) {
        LOGF(fatal, notify("Histogram \"%s\" has wrong dimension %d != %d"), config->confEffCorCCDBPath.value, dim, config->confEffCorVariables.value.size());
        return weight;
      }

      auto bin = -1;
      if (config->confEffCorVariables.value == "pt") {
        bin = hLoaded[partNo - 1]->FindBin(particle.pt());
      } else if (config->confEffCorVariables.value == "pt,eta") {
        bin = hLoaded[partNo - 1]->FindBin(particle.pt(), particle.eta());
      } else if (config->confEffCorVariables.value == "pt,mult") {
        bin = hLoaded[partNo - 1]->FindBin(particle.pt(), particle.fdCollision().multV0M());
      } else if (config->confEffCorVariables.value == "pt,eta,mult") {
        bin = hLoaded[partNo - 1]->FindBin(particle.pt(), particle.eta(), particle.fdCollision().multV0M());
      } else {
        LOGF(fatal, notify("Unknown configuration for efficiency variables"));
        return weight;
      }
      weight = hWeights->GetBinContent(bin);
    }

    return weight;
  }

 private:
  static inline auto notify(const std::string& msg) -> const std::string
  {
    return fmt::format("[EFFICIENCY CORRECTION] {}", msg);
  }

  static auto isHistEmpty(TH1* hist) -> bool
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

  template <typename H>
  auto loadHistFromCCDB(const int64_t timestamp) const -> H*
  {
    auto hWeights = ccdb.getForTimeStamp<H>(config->confEffCorCCDBPath, timestamp);
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

  auto getDimensionFromVariables() -> size_t
  {
    auto parts = std::views::split(config->confEffCorVariables.value, ',');
    return std::ranges::distance(parts);
  }

  EffCorConfigurableGroup* config{};

  bool shouldApplyCorrection{false};
  bool shouldFillHistograms{false};

  o2::ccdb::BasicCCDBManager& ccdb{o2::ccdb::BasicCCDBManager::instance()};
  std::array<TH1*, 2> hLoaded{};

  framework::HistogramRegistry* histRegistry{};
  static constexpr std::string_view histDirectory{"EfficiencyCorrection"};
  static constexpr std::string_view histSuffix[2]{"one", "two"};
};

} // namespace o2::analysis::femto_universe::efficiency_correction

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEEFFICIENCYCORRECTION_H_
