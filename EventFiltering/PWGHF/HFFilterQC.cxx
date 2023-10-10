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

#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h> // needed for HFFilterHelpers, to be fixed

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"

#include "EventFiltering/filterTables.h"
#include "EventFiltering/PWGHF/HFFilterHelpers.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hffilters;
using namespace filtering;

struct HfFilterQc { // Main struct for HF trigger QC

  std::array<std::shared_ptr<TH2>, kNtriggersHF + 2> hPartPerEvent{};
  std::array<std::shared_ptr<TH2>, kNtriggersHF + 2> hPtDistr{};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(InitContext&)
  {
    // Initialize the histograms
    hPartPerEvent[0] = registry.add<TH2>("hPartPerEventAll", "All events;;number of particles", HistType::kTH2F, {{kNCharmParticles, -0.5, kNCharmParticles - 0.5}, {11, -0.5, 10.5}});
    hPartPerEvent[1] = registry.add<TH2>("hPartPerEventTriggered", "HF triggered events;;number of particles", HistType::kTH2F, {{kNCharmParticles, -0.5, kNCharmParticles - 0.5}, {11, -0.5, 10.5}});
    hPtDistr[0] = registry.add<TH2>("hPtDistrAll", "All events;;#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{kNCharmParticles, -0.5, kNCharmParticles - 0.5}, {50, -0.5, 10.5}});
    hPtDistr[1] = registry.add<TH2>("hPtDistrTriggered", "HF triggered events;;#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{kNCharmParticles, -0.5, kNCharmParticles - 0.5}, {11, -0.5, 10.5}});
    for (auto iTrig = 0; iTrig < kNtriggersHF; ++iTrig) {
      hPartPerEvent[iTrig + 2] = registry.add<TH2>(Form("hPartPerEvent%s", hfTriggerNames[iTrig].data()), Form("%s Filtered events;;number of particles", hfTriggerNames[iTrig].data()), HistType::kTH2F, {{kNCharmParticles, -0.5, kNCharmParticles - 0.5}, {11, -0.5, 10.5}});
      hPtDistr[iTrig + 2] = registry.add<TH2>(Form("hPtDistr%s", hfTriggerNames[iTrig].data()), Form("%s Filtered events;;#it{p}_{T} (GeV/#it{c})", hfTriggerNames[iTrig].data()), HistType::kTH2F, {{kNCharmParticles, -0.5, kNCharmParticles - 0.5}, {11, -0.5, 10.5}});
    }
    for (auto iTrig = 0; iTrig < kNtriggersHF + 1; ++iTrig) {
      for (auto iBin = 0; iBin < kNCharmParticles; ++iBin) {
        hPartPerEvent[iTrig]->GetXaxis()->SetBinLabel(iBin + 1, charmParticleNames[iBin].data());
        hPtDistr[iTrig]->GetXaxis()->SetBinLabel(iBin + 1, charmParticleNames[iBin].data());
      }
    }
  }

  /// Loops over particle species and checks whether the analysed particle is the correct one
  /// \param pdgDau  tuple with PDG daughter codes for the desired decay
  /// \param mcParticles  table with MC particles
  /// \param particle  MC particle
  /// \param nParticles  array with number of particles found for each particle species
  /// \param triggerDecision  array with trigger decision
  template <size_t I = 0, typename... Ts, typename T, typename U, typename A>
  constexpr void checkParticleDecay(std::tuple<Ts...> pdgDau,
                                    const T& mcParticles,
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
      if (RecoDecay::isMatchedMCGen(mcParticles, particle, pdgCodesCharm[I], std::get<I>(pdgDau), true, &sign, 2)) {
        nParticles[I]++;
        hPtDistr[0]->Fill(static_cast<double>(I), static_cast<double>(particle.pt()));
        for (auto iTrig = 0; iTrig < kNtriggersHF; ++iTrig) {
          if (triggerDecision[iTrig]) {
            hPtDistr[iTrig + 1]->Fill(static_cast<double>(I), static_cast<double>(particle.pt()));
          }
        }
      }

      // Going for next element.
      checkParticleDecay<I + 1>(pdgDau, mcParticles, particle, nParticles, triggerDecision);
    }
  }

  void process(HfFilter const& filterDecision,
               McParticles const& mcParticles)
  {
    bool hasHighPt2P = filterDecision.hasHfHighPt2P();
    bool hasHighPt3P = filterDecision.hasHfHighPt3P();
    bool hasBeauty3P = filterDecision.hasHfBeauty3P();
    bool hasBeauty4P = filterDecision.hasHfBeauty4P();
    bool hasFemto2P = filterDecision.hasHfFemto2P();
    bool hasFemto3P = filterDecision.hasHfFemto3P();
    bool hasDoubleCharm2P = filterDecision.hasHfDoubleCharm2P();
    bool hasDoubleCharm3P = filterDecision.hasHfDoubleCharm3P();
    bool hasDoubleCharmMix = filterDecision.hasHfDoubleCharmMix();
    bool hasHfV02P = filterDecision.hasHfV0Charm2P();
    bool hasHfV03P = filterDecision.hasHfV0Charm3P();
    bool hasCharmBarToXiBach = filterDecision.hasHfCharmBarToXiBach();
    bool isTriggered = hasHighPt2P || hasHighPt3P || hasBeauty3P || hasBeauty4P || hasFemto2P || hasFemto3P || hasDoubleCharm2P || hasDoubleCharm3P || hasDoubleCharmMix || hasHfV02P || hasHfV03P || hasCharmBarToXiBach;
    auto triggerDecision = std::array{isTriggered, hasHighPt2P, hasHighPt3P, hasBeauty3P, hasBeauty4P, hasFemto2P, hasFemto3P, hasDoubleCharm2P, hasDoubleCharm3P, hasDoubleCharmMix, hasHfV02P, hasHfV03P, hasCharmBarToXiBach};

    std::array<int, kNCharmParticles> nPart{0};
    // Loop over the MC particles
    for (const auto& particle : mcParticles) {
      // Check if the particle is of interest
      checkParticleDecay(pdgCharmDaughters, mcParticles, particle, nPart, triggerDecision);
    }

    for (auto iPart = 0; iPart < kNCharmParticles; ++iPart) {
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
