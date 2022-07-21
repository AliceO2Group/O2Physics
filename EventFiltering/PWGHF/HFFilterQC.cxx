// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
// O2 includes

/// \file HFFilterQC.cxx
/// \brief task for the quality assurance of the event selection with HFFilter.cxx
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
/// \author Biao Zhang <biao.zhang@cern.ch>, CCNU
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "EventFiltering/filterTables.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::filtering;

namespace
{

enum HfTriggers {
  kHighPt = 0,
  kBeauty,
  kFemto,
  kDoubleCharm,
  kNtriggersHF
};

enum particles {
  kD0 = 0,
  kDplus,
  kDs,
  kLc,
  kXic,
  kNParticles
};

static const std::array<std::string, kNtriggersHF> HfTriggerNames{"HighPt", "Beauty", "Femto", "DoubleCharm"};
static const std::array<std::string, kNParticles> particleNames{"D0", "Dplus", "Ds", "Lc", "Xic"};
static const std::array<int, kNParticles> pdgCodes{421, 411, 431, 4122, 4232};
static const std::tuple pdgDaughters{
  std::array{-321, 211},        // D0
  std::array{-321, 211, 211},   // Dplus
  std::array{321, -321, 211},   // Ds
  std::array{2212, -321, 211},  // Lc
  std::array{2212, -321, 211}}; // Xic

}; // namespace

struct HfFilterQc { // Main struct for HF trigger QC

  std::array<std::shared_ptr<TH2>, kNtriggersHF + 2> hPartPerEvent{};
  std::array<std::shared_ptr<TH2>, kNtriggersHF + 2> hPtDistr{};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(o2::framework::InitContext&)
  {
    // Initialize the histograms
    hPartPerEvent[0] = registry.add<TH2>("hPartPerEventAll", "All events;;number of particles", HistType::kTH2F, {{kNParticles, -0.5, kNParticles - 0.5}, {11, -0.5, 10.5}});
    hPartPerEvent[1] = registry.add<TH2>("hPartPerEventTriggered", "HF triggered events;;number of particles", HistType::kTH2F, {{kNParticles, -0.5, kNParticles - 0.5}, {11, -0.5, 10.5}});
    hPtDistr[0] = registry.add<TH2>("hPtDistrAll", "All events;;#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{kNParticles, -0.5, kNParticles - 0.5}, {50, -0.5, 10.5}});
    hPtDistr[1] = registry.add<TH2>("hPtDistrTriggered", "HF triggered events;;#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{kNParticles, -0.5, kNParticles - 0.5}, {11, -0.5, 10.5}});
    for (auto iTrig = 0; iTrig < kNtriggersHF; ++iTrig) {
      hPartPerEvent[iTrig + 2] = registry.add<TH2>(Form("hPartPerEvent%s", HfTriggerNames[iTrig].data()), Form("%s Filtered events;;number of particles", HfTriggerNames[iTrig].data()), HistType::kTH2F, {{kNParticles, -0.5, kNParticles - 0.5}, {11, -0.5, 10.5}});
      hPtDistr[iTrig + 2] = registry.add<TH2>(Form("hPtDistr%s", HfTriggerNames[iTrig].data()), Form("%s Filtered events;;#it{p}_{T} (GeV/#it{c})", HfTriggerNames[iTrig].data()), HistType::kTH2F, {{kNParticles, -0.5, kNParticles - 0.5}, {11, -0.5, 10.5}});
    }
    for (auto iTrig = 0; iTrig < kNtriggersHF + 1; ++iTrig) {
      for (auto iBin = 0; iBin < kNParticles; ++iBin) {
        hPartPerEvent[iTrig]->GetXaxis()->SetBinLabel(iBin + 1, particleNames[iBin].data());
        hPtDistr[iTrig]->GetXaxis()->SetBinLabel(iBin + 1, particleNames[iBin].data());
      }
    }
  }

  /// Loops over particle species and checks whether the analysed particle is the correct one
  /// \param pdgDau  tuple with PDG daughter codes for the desired decay
  /// \param particlesMC  table with MC particles
  /// \param particle  MC particle
  /// \param nParticles  array with number of particles found for each particle species
  /// \param triggerDecision  array with trigger decision
  template <size_t I = 0, typename... Ts, typename T, typename U, typename A>
  constexpr void checkParticleDecay(std::tuple<Ts...> pdgDau,
                                    const T& particlesMC,
                                    const U& particle,
                                    A& nParticles,
                                    const std::array<bool, kNtriggersHF + 1>& triggerDecision)
  {
    // If we have iterated through all elements
    if constexpr (I == sizeof...(Ts)) {
      // Last case, if nothing is left to
      // iterate, then exit the function
      return;
    } else {
      int8_t sign = 0;
      if (RecoDecay::isMatchedMCGen(particlesMC, particle, pdgCodes[I], std::get<I>(pdgDau), true, &sign, 2)) {
        nParticles[I]++;
        hPtDistr[0]->Fill(static_cast<double>(I), static_cast<double>(particle.pt()));
        for (auto iTrig = 0; iTrig < kNtriggersHF; ++iTrig) {
          if (triggerDecision[iTrig]) {
            hPtDistr[iTrig + 1]->Fill(static_cast<double>(I), static_cast<double>(particle.pt()));
          }
        }
      }

      // Going for next element.
      checkParticleDecay<I + 1>(pdgDau, particlesMC, particle, nParticles, triggerDecision);
    }
  }

  void process(HfFilter const& filterDecision,
               McParticles const& particlesMC)
  {
    bool hasHighPt = filterDecision.hasHfHighPt();
    bool hasBeauty = filterDecision.hasHfBeauty();
    bool hasFemto = filterDecision.hasHfFemto();
    bool hasDoubleCharm = filterDecision.hasHfDoubleCharm();
    bool isTriggered = hasHighPt || hasBeauty || hasFemto || hasDoubleCharm;
    auto triggerDecision = std::array{isTriggered, hasHighPt, hasBeauty, hasFemto, hasDoubleCharm};

    std::array<int, kNParticles> nPart{0};
    // Loop over the MC particles
    for (auto const& particle : particlesMC) {
      // Check if the particle is of interest
      checkParticleDecay(pdgDaughters, particlesMC, particle, nPart, triggerDecision);
    }

    for (auto iPart = 0; iPart < kNParticles; ++iPart) {
      hPartPerEvent[0]->Fill(iPart, nPart[iPart]);
      for (auto iTrig = 0; iTrig < kNtriggersHF; ++iTrig) {
        if (triggerDecision[iTrig]) {
          hPartPerEvent[iTrig + 1]->Fill(iPart, nPart[iPart]);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfFilterQc>(cfg));

  return workflow;
}
