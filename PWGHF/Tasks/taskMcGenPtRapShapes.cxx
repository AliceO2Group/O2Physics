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

/// \file taskMcGenPtRapShapes.cxx
/// \brief Task to check generated pt and y distributions of charm and beauty hadrons in MC productions

/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include <array>
#include <vector>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/RecoDecay.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

namespace
{
const int nCharmHadrons = 10;
static constexpr std::array<int, nCharmHadrons> pdgCodesCharm = {Pdg::kD0, Pdg::kDPlus, Pdg::kDS, Pdg::kDStar, Pdg::kLambdaCPlus, Pdg::kSigmaC0, Pdg::kSigmaCPlusPlus, Pdg::kXiC0, Pdg::kXiCPlus, Pdg::kOmegaC0};

const int nBeautyHadrons = 4;
static constexpr std::array<int, nBeautyHadrons> pdgCodesBeauty = {Pdg::kB0, Pdg::kBPlus, Pdg::kBS, Pdg::kLambdaB0};
} // namespace

struct HfTaskMcGenPtRapShapes {

  ConfigurableAxis axisPtCharm{"axisPtCharm", {1000, 0.f, 100.f}, "Binning for the pT axis of charm hadrons"};
  ConfigurableAxis axisPtBeauty{"axisPtBeauty", {3000, 0.f, 300.f}, "Binning for the pT axis of beauty hadrons"};
  ConfigurableAxis axisRapCharm{"axisRapCharm", {100, -1.f, 1.f}, "Binning for the y axis of charm hadrons"};
  ConfigurableAxis axisRapBeauty{"axisRapBeauty", {100, -1.f, 1.f}, "Binning for the y axis of beauty hadrons"};

  std::array<std::shared_ptr<TH2>, nCharmHadrons> histRapVsPtCharmPrompt{};
  std::array<std::shared_ptr<TH2>, nCharmHadrons> histRapVsPtCharmNonPrompt{};
  std::array<std::shared_ptr<TH2>, nCharmHadrons> histPtCharmVsPtBeautyNonPrompt{};
  std::array<std::shared_ptr<TH2>, nBeautyHadrons> histRapVsPtBeauty{};

  HistogramRegistry registry{};

  void init(InitContext&)
  {

    for (auto iCharmHad{0}; iCharmHad < nCharmHadrons; ++iCharmHad) {
      histRapVsPtCharmPrompt[iCharmHad] = registry.add<TH2>(Form("CharmHadrons/hRapVsPtPrompt%d", pdgCodesCharm[iCharmHad]), Form("Prompt %d;#it{p}_{T} (GeV/#it{c});#it{y}", pdgCodesCharm[iCharmHad]), {HistType::kTH2F, {axisPtCharm, axisRapCharm}});
      histRapVsPtCharmNonPrompt[iCharmHad] = registry.add<TH2>(Form("CharmHadrons/hRapVsPtNonPrompt%d", pdgCodesCharm[iCharmHad]), Form("Non-prompt %d;#it{p}_{T} (GeV/#it{c});#it{y}", pdgCodesCharm[iCharmHad]), {HistType::kTH2F, {axisPtCharm, axisRapCharm}});
      histPtCharmVsPtBeautyNonPrompt[iCharmHad] = registry.add<TH2>(Form("CharmHadrons/hPtCharmVsPtBeautyNonPrompt%d", pdgCodesCharm[iCharmHad]), Form("Non-prompt %d;#it{p}_{T}(b-had) (GeV/#it{c});#it{p}_{T}(c-had) (GeV/#it{c})", pdgCodesCharm[iCharmHad]), {HistType::kTH2F, {axisPtBeauty, axisPtCharm}});
    }

    for (auto iBeautyHad{0}; iBeautyHad < nBeautyHadrons; ++iBeautyHad) {
      histRapVsPtBeauty[iBeautyHad] = registry.add<TH2>(Form("BeautyHadrons/hRapVsPt%d", pdgCodesBeauty[iBeautyHad]), Form("%d;#it{p}_{T} (GeV/#it{c});#it{y}", pdgCodesBeauty[iBeautyHad]), {HistType::kTH2F, {axisPtBeauty, axisRapBeauty}});
    }
  }

  void process(aod::McParticles const& mcParticles)
  {
    for (auto const& mcParticle : mcParticles) {
      int absPdgCode = std::abs(mcParticle.pdgCode());
      float pt = mcParticle.pt();
      float rap = mcParticle.y();
      auto itCharm = std::find(pdgCodesCharm.begin(), pdgCodesCharm.end(), absPdgCode);
      auto itBeauty = std::find(pdgCodesBeauty.begin(), pdgCodesBeauty.end(), absPdgCode);
      if (itCharm != pdgCodesCharm.end()) {
        auto idxCharm = std::distance(pdgCodesCharm.begin(), itCharm);
        std::vector<int> idxBhadMothers{};
        auto origin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle, false, &idxBhadMothers);
        if (origin == RecoDecay::OriginType::Prompt) {
          histRapVsPtCharmPrompt[idxCharm]->Fill(pt, rap);
        } else if (origin == RecoDecay::OriginType::NonPrompt) {
          histRapVsPtCharmNonPrompt[idxCharm]->Fill(pt, rap);
          if (std::abs(rap) < 0.5) {
            auto bMother = mcParticles.rawIteratorAt(idxBhadMothers[0]);
            histPtCharmVsPtBeautyNonPrompt[idxCharm]->Fill(bMother.pt(), pt);
          }
        }
      } else if (itBeauty != pdgCodesBeauty.end()) {
        auto idxBeauty = std::distance(pdgCodesBeauty.begin(), itBeauty);
        histRapVsPtBeauty[idxBeauty]->Fill(pt, rap);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskMcGenPtRapShapes>(cfgc)};
}
