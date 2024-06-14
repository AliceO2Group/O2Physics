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

/// \file taskMcValidation.cxx
/// \brief Monte Carlo validation task
/// \note generated and reconstructed level validation
///
/// \author Antonio Palasciano <antonio.palasciano@cern.ch>, Università degli Studi di Bari & INFN, Sezione di Bari
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
/// \author Fabio Catalano <fabio.catalano@cern.ch>, CERN

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/StaticFor.h"

#include "Common/DataModel/CollisionAssociationTables.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
enum DecayChannels { DzeroToKPi = 0,
                     DstarToDzeroPi,
                     DplusToPiKPi,
                     DplusToPhiPiToKKPi, // resonant channel with Phi, Dplus -> PhiPi -> KKPi
                     DsToPhiPiToKKPi,    // resonant channel with Phi, Ds -> PhiPi -> KKPi
                     DsToK0starKToKKPi,  // resonant channel with K0*, Ds -> K0*K -> KKPi
                     LcToPKPi,
                     LcToPiK0s,
                     XiCplusToPKPi,
                     XiCplusToXiPiPi,
                     XiCzeroToXiPi,
                     OmegaCToOmegaPi,
                     OmegaCToXiPi,
                     nChannels }; // always keep nChannels at the end

static constexpr int nMesonChannels = 6; // number of meson channels
static constexpr int nOriginTypes = 2;   // number of origin types (prompt, non-prompt)
static constexpr std::array<int, nChannels> PDGArrayParticle = {o2::constants::physics::Pdg::kD0, o2::constants::physics::Pdg::kDStar, o2::constants::physics::Pdg::kDPlus, o2::constants::physics::Pdg::kDPlus,
                                                                o2::constants::physics::Pdg::kDS, o2::constants::physics::Pdg::kDS, o2::constants::physics::Pdg::kLambdaCPlus, o2::constants::physics::Pdg::kLambdaCPlus,
                                                                o2::constants::physics::Pdg::kXiCPlus, o2::constants::physics::Pdg::kXiCPlus, o2::constants::physics::Pdg::kXiC0, o2::constants::physics::Pdg::kOmegaC0, o2::constants::physics::Pdg::kOmegaC0};
static constexpr std::array<unsigned int, nChannels> nDaughters = {2, 3, 3, 3, 3, 3, 3, 3, 3, 5, 4, 4, 4};
static constexpr std::array<int, nChannels> maxDepthForSearch = {1, 2, 2, 2, 2, 2, 2, 3, 2, 4, 3, 3, 3};
// keep coherent indexing with PDGArrayParticle
// FIXME: look for a better solution
static constexpr std::array<std::array<int, 2>, nChannels> arrPDGFinal2Prong = {{{+kPiPlus, -kKPlus}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}};
static constexpr std::array<std::array<int, 3>, nChannels> arrPDGFinal3Prong = {{{}, {+kPiPlus, -kKPlus, +kPiPlus}, {+kPiPlus, -kKPlus, +kPiPlus}, {+kKPlus, -kKPlus, +kPiPlus}, {+kKPlus, -kKPlus, +kPiPlus}, {+kKPlus, -kKPlus, +kPiPlus}, {+kProton, -kKPlus, +kPiPlus}, {+kProton, -kPiPlus, +kPiPlus}, {+kProton, -kKPlus, +kPiPlus}, {}, {}, {}, {}}};
static constexpr std::array<std::array<int, 4>, nChannels> arrPDGFinal4Prong = {{{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {+kPiPlus, -kPiPlus, -kPiPlus, +kProton}, {+kPiPlus, -kKPlus, -kPiPlus, +kProton}, {+kPiPlus, -kPiPlus, -kPiPlus, +kProton}}};
static constexpr std::array<std::array<int, 5>, nChannels> arrPDGFinal5Prong = {{{}, {}, {}, {}, {}, {}, {}, {}, {}, {+kPiPlus, +kPiPlus, -kPiPlus, -kPiPlus, +kProton}, {}, {}, {}}};
static constexpr std::string_view labels[nChannels] = {"D^{0} #rightarrow K#pi", "D*^{+} #rightarrow D^{0}#pi", "D^{+} #rightarrow K#pi#pi", "D^{+} #rightarrow KK#pi", "D_{s}^{+} #rightarrow #Phi#pi #rightarrow KK#pi",
                                                       "D_{s}^{+} #rightarrow #bar{K}^{*0}K #rightarrow KK#pi", "#Lambda_{c}^{+} #rightarrow pK#pi", "#Lambda_{c}^{+} #rightarrow pK^{0}_{s}", "#Xi_{c}^{+} #rightarrow pK#pi",
                                                       "#Xi_{c}^{+} #rightarrow #Xi#pi#pi", "#Xi_{c}^{0} #rightarrow #Xi#pi", "#Omega_{c}^{0} #rightarrow #Omega#pi", "#Omega_{c}^{0} #rightarrow #Xi#pi"};
static constexpr std::string_view particleNames[nChannels] = {"DzeroToKPi", "DstarToDzeroPi", "DplusToPiKPi", "DplusToPhiPiToKKPi", "DsToPhiPiToKKPi", "DsToK0starKToKKPi",
                                                              "LcToPKPi", "LcToPiK0s", "XiCplusToPKPi", "XiCplusToXiPiPi", "XiCzeroToXiPi", "OmegaCToOmegaPi", "OmegaCToXiPi"};
static constexpr std::string_view originNames[nOriginTypes] = {"Prompt", "NonPrompt"};
} // namespace

/// Generated Level Validation
///
/// - Number of HF quarks produced per collision
/// - Number of candidates per collision
/// - Momentum Conservation for these particles
struct HfTaskMcValidationGen {
  Configurable<double> xVertexMin{"xVertexMin", -100., "min. x of generated primary vertex [cm]"};
  Configurable<double> xVertexMax{"xVertexMax", 100., "max. x of generated primary vertex [cm]"};
  Configurable<double> yVertexMin{"yVertexMin", -100., "min. y of generated primary vertex [cm]"};
  Configurable<double> yVertexMax{"yVertexMax", 100., "max. y of generated primary vertex [cm]"};
  Configurable<double> zVertexMin{"zVertexMin", -100., "min. z of generated primary vertex [cm]"};
  Configurable<double> zVertexMax{"zVertexMax", 100., "max. z of generated primary vertex [cm]"};
  Configurable<int> eventGeneratorType{"eventGeneratorType", -1, "If positive, enable event selection using subGeneratorId information. The value indicates which events to keep (0 = MB, 4 = charm triggered, 5 = beauty triggered)"};

  AxisSpec axisNhadrons{10, -0.5, 9.5};
  AxisSpec axisNquarks{20, -0.5, 19.5};
  AxisSpec axisResiduals{100, -0.01, 0.01};
  AxisSpec axisPt{100, 0., 50.};
  AxisSpec axisY{100, -5., 5.};
  AxisSpec axisMesonSpecies{nMesonChannels, -0.5, static_cast<float>(nMesonChannels) - 0.5};
  AxisSpec axisBaryonSpecies{nChannels - nMesonChannels, -0.5, static_cast<float>(nChannels - nMesonChannels) - 0.5};
  AxisSpec axisDecLen{100, 0., 10000.};

  HistogramRegistry registry{
    "registry",
    {{"hMomentumCheck", "Mom. Conservation (1 = true, 0 = false) (#it{#epsilon} = 1 MeV/#it{c}); Mom. Conservation result; entries", {HistType::kTH1F, {{2, -0.5, +1.5}}}},
     {"hPtDiffMotherDaughterGen", "Pt Difference Mother-Daughters; #Delta#it{p}_{T}^{gen} (GeV/#it{c}); entries", {HistType::kTH1F, {axisResiduals}}},
     {"hPxDiffMotherDaughterGen", "Px Difference Mother-Daughters; #Delta#it{p}_{x}^{gen} (GeV/#it{c}); entries", {HistType::kTH1F, {axisResiduals}}},
     {"hPyDiffMotherDaughterGen", "Py Difference Mother-Daughters; #Delta#it{p}_{y}^{gen} (GeV/#it{c}); entries", {HistType::kTH1F, {axisResiduals}}},
     {"hPzDiffMotherDaughterGen", "Pz Difference Mother-Daughters; #Delta#it{p}_{z}^{gen} (GeV/#it{c}); entries", {HistType::kTH1F, {axisResiduals}}},
     {"hPDiffMotherDaughterGen", "P  Difference Mother-Daughters; #Delta#it{p}^{gen} (GeV/#it{c}); entries", {HistType::kTH1F, {axisResiduals}}},
     {"Quarks/hCountC", "Event counter - Number of charm quarks; Events Per Collision; entries", {HistType::kTH1F, {axisNquarks}}},
     {"Quarks/hCountB", "Event counter - Number of beauty quarks; Events Per Collision; entries", {HistType::kTH1F, {axisNquarks}}},
     {"Quarks/hCountCbar", "Event counter - Number of anti-charm quarks; Events Per Collision; entries", {HistType::kTH1F, {axisNquarks}}},
     {"Quarks/hCountBbar", "Event counter - Number of anti-beauty quarks; Events Per Collision; entries", {HistType::kTH1F, {axisNquarks}}},
     {"Quarks/hPtVsYCharmQuark", "Y vs. Pt - charm quarks ; #it{p}_{T}^{gen} (GeV/#it{c}); #it{y}^{gen}", {HistType::kTH2F, {axisPt, axisY}}},
     {"Quarks/hPtVsYBeautyQuark", "Y vs. Pt - beauty quarks ; #it{p}_{T}^{gen} (GeV/#it{c}); #it{y}^{gen}", {HistType::kTH2F, {axisPt, axisY}}},
     {"PromptCharmMesons/hPromptMesonsPtDistr", "Pt distribution vs prompt charm meson in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c})", {HistType::kTH2F, {axisMesonSpecies, axisPt}}},
     {"PromptCharmMesons/hPromptMesonsYDistr", "Y distribution vs prompt charm meson; ; #it{y}^{gen}", {HistType::kTH2F, {axisMesonSpecies, axisY}}},
     {"PromptCharmMesons/hPromptMesonsDecLenDistr", "Decay length distribution vs prompt charm meson; ; decay length (#mum)", {HistType::kTH2F, {axisMesonSpecies, axisDecLen}}},
     {"NonPromptCharmMesons/hNonPromptMesonsPtDistr", "Pt distribution vs non-prompt charm meson in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c})", {HistType::kTH2F, {axisMesonSpecies, axisPt}}},
     {"NonPromptCharmMesons/hNonPromptMesonsYDistr", "Y distribution vs non-prompt charm meson; ; #it{y}^{gen}", {HistType::kTH2F, {axisMesonSpecies, axisY}}},
     {"NonPromptCharmMesons/hNonPromptMesonsDecLenDistr", "Decay length distribution vs non-prompt charm meson; ; decay length (#mum)", {HistType::kTH2F, {axisMesonSpecies, axisDecLen}}},
     {"PromptCharmBaryons/hPromptBaryonsPtDistr", "Pt distribution vs prompt charm baryon in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c})", {HistType::kTH2F, {axisBaryonSpecies, axisPt}}},
     {"PromptCharmBaryons/hPromptBaryonsYDistr", "Y distribution vs prompt charm baryon; ; #it{y}^{gen}", {HistType::kTH2F, {axisBaryonSpecies, axisY}}},
     {"PromptCharmBaryons/hPromptBaryonsDecLenDistr", "Decay length distribution vs prompt charm baryon; ; decay length (#mum)", {HistType::kTH2F, {axisBaryonSpecies, axisDecLen}}},
     {"NonPromptCharmBaryons/hNonPromptBaryonsPtDistr", "Pt distribution vs non-prompt charm baryon in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c})", {HistType::kTH2F, {axisBaryonSpecies, axisPt}}},
     {"NonPromptCharmBaryons/hNonPromptBaryonsYDistr", "Y distribution vs non-prompt charm baryon; ; #it{y}^{gen}", {HistType::kTH2F, {axisBaryonSpecies, axisY}}},
     {"NonPromptCharmBaryons/hNonPromptBaryonsDecLenDistr", "Decay length distribution vs non-prompt charm baryon; ; decay length (#mum)", {HistType::kTH2F, {axisBaryonSpecies, axisDecLen}}}}};

  void init(InitContext&)
  {
    // add per species histograms
    for (size_t iOrigin = 0; iOrigin < nOriginTypes; iOrigin++) {
      for (int iChannel = 0; iChannel < nMesonChannels; iChannel++) {
        registry.add(Form("%sCharmMesons/hCount%s%s", originNames[iOrigin].data(), originNames[iOrigin].data(), particleNames[iChannel].data()),
                     Form("Event counter - %s %s; Events Per Collision; entries", originNames[iOrigin].data(), labels[iChannel].data()), {HistType::kTH1F, {axisNhadrons}});
      }
      for (int iChannel = nMesonChannels; iChannel < nChannels; iChannel++) {
        registry.add(Form("%sCharmBaryons/hCount%s%s", originNames[iOrigin].data(), originNames[iOrigin].data(), particleNames[iChannel].data()),
                     Form("Event counter - %s %s; Events Per Collision; entries", originNames[iOrigin].data(), labels[iChannel].data()), {HistType::kTH1F, {axisNhadrons}});
      }
    }

    for (auto iBin = 1; iBin <= nMesonChannels; ++iBin) {
      registry.get<TH2>(HIST("PromptCharmMesons/hPromptMesonsPtDistr"))->GetXaxis()->SetBinLabel(iBin, labels[iBin - 1].data());
      registry.get<TH2>(HIST("PromptCharmMesons/hPromptMesonsYDistr"))->GetXaxis()->SetBinLabel(iBin, labels[iBin - 1].data());
      registry.get<TH2>(HIST("PromptCharmMesons/hPromptMesonsDecLenDistr"))->GetXaxis()->SetBinLabel(iBin, labels[iBin - 1].data());
      registry.get<TH2>(HIST("NonPromptCharmMesons/hNonPromptMesonsPtDistr"))->GetXaxis()->SetBinLabel(iBin, labels[iBin - 1].data());
      registry.get<TH2>(HIST("NonPromptCharmMesons/hNonPromptMesonsYDistr"))->GetXaxis()->SetBinLabel(iBin, labels[iBin - 1].data());
      registry.get<TH2>(HIST("NonPromptCharmMesons/hNonPromptMesonsDecLenDistr"))->GetXaxis()->SetBinLabel(iBin, labels[iBin - 1].data());
    }
    for (auto iBin = 1; iBin <= nChannels - nMesonChannels; ++iBin) {
      registry.get<TH2>(HIST("PromptCharmBaryons/hPromptBaryonsPtDistr"))->GetXaxis()->SetBinLabel(iBin, labels[iBin + nMesonChannels - 1].data());
      registry.get<TH2>(HIST("PromptCharmBaryons/hPromptBaryonsYDistr"))->GetXaxis()->SetBinLabel(iBin, labels[iBin + nMesonChannels - 1].data());
      registry.get<TH2>(HIST("PromptCharmBaryons/hPromptBaryonsDecLenDistr"))->GetXaxis()->SetBinLabel(iBin, labels[iBin + nMesonChannels - 1].data());
      registry.get<TH2>(HIST("NonPromptCharmBaryons/hNonPromptBaryonsPtDistr"))->GetXaxis()->SetBinLabel(iBin, labels[iBin + nMesonChannels - 1].data());
      registry.get<TH2>(HIST("NonPromptCharmBaryons/hNonPromptBaryonsYDistr"))->GetXaxis()->SetBinLabel(iBin, labels[iBin + nMesonChannels - 1].data());
      registry.get<TH2>(HIST("NonPromptCharmBaryons/hNonPromptBaryonsDecLenDistr"))->GetXaxis()->SetBinLabel(iBin, labels[iBin + nMesonChannels - 1].data());
    }
  }

  /// Primary-vertex selection
  /// \param collision  mcCollision table row
  template <typename Col>
  bool selectVertex(const Col& collision)
  {
    // x position
    if (collision.posX() < xVertexMin || collision.posX() > xVertexMax) {
      return false;
    }
    // y position
    if (collision.posY() < yVertexMin || collision.posY() > yVertexMax) {
      return false;
    }
    // z position
    if (collision.posZ() < zVertexMin || collision.posZ() > zVertexMax) {
      return false;
    }
    return true;
  }

  void process(aod::McCollision const& mcCollision,
               aod::McParticles const& mcParticles)
  {
    if (eventGeneratorType >= 0 && mcCollision.getSubGeneratorId() != eventGeneratorType) {
      return;
    }
    if (!selectVertex(mcCollision)) {
      return;
    }

    int cPerCollision = 0;
    int cBarPerCollision = 0;
    int bPerCollision = 0;
    int bBarPerCollision = 0;
    std::array<int, nChannels> counterPrompt{0}, counterNonPrompt{0};

    for (const auto& particle : mcParticles) {
      if (!particle.has_mothers()) {
        continue;
      }

      int particlePdgCode = particle.pdgCode();
      bool isDiffFromMothers = true;
      for (const auto& mother : particle.mothers_as<aod::McParticles>()) {
        if (particlePdgCode == mother.pdgCode()) {
          isDiffFromMothers = false;
          break;
        }
      }
      if (isDiffFromMothers) {
        switch (particlePdgCode) {
          case kCharm:
            cPerCollision++;
            registry.fill(HIST("Quarks/hPtVsYCharmQuark"), particle.pt(), particle.y());
            break;
          case kCharmBar:
            cBarPerCollision++;
            registry.fill(HIST("Quarks/hPtVsYCharmQuark"), particle.pt(), particle.y());
            break;
          case kBottom:
            bPerCollision++;
            registry.fill(HIST("Quarks/hPtVsYBeautyQuark"), particle.pt(), particle.y());
            break;
          case kBottomBar:
            bBarPerCollision++;
            registry.fill(HIST("Quarks/hPtVsYBeautyQuark"), particle.pt(), particle.y());
            break;
        }
      }

      // Checking the decay of the particles and the momentum conservation
      for (std::size_t iD = 0; iD < PDGArrayParticle.size(); ++iD) {
        if (std::abs(particlePdgCode) != PDGArrayParticle[iD]) {
          continue;
        }

        std::vector<int> listDaughters{};

        // Check that the decay channel is correct and retrieve the daughters
        if (nDaughters[iD] == 2) {
          if (!RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], arrPDGFinal2Prong[iD], true, nullptr, maxDepthForSearch[iD], &listDaughters)) {
            continue;
          }
        }

        if (nDaughters[iD] == 3) {
          if (!RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, PDGArrayParticle[iD], arrPDGFinal3Prong[iD], true, nullptr, maxDepthForSearch[iD], &listDaughters)) {
            continue;
          }
          if (iD == DstarToDzeroPi &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], std::array{+o2::constants::physics::Pdg::kD0, +kPiPlus}, true)) {
            continue;
          }
          if ((iD == DplusToPhiPiToKKPi || iD == DsToPhiPiToKKPi) &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], std::array{+o2::constants::physics::Pdg::kPhi, +kPiPlus}, false) &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, -PDGArrayParticle[iD], std::array{+o2::constants::physics::Pdg::kPhi, -kPiPlus}, false)) {
            continue;
          }
          if (iD == DsToK0starKToKKPi &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], std::array{-313, +kKPlus}, true)) {
            continue;
          }
          if (iD == LcToPiK0s &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], std::array{+kK0Short, +kProton}, false, nullptr, 2) &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, -PDGArrayParticle[iD], std::array{+kK0Short, -kProton}, false, nullptr, 2)) {
            continue;
          }
        }

        if (nDaughters[iD] == 4) {
          if (!RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, PDGArrayParticle[iD], arrPDGFinal4Prong[iD], true, nullptr, maxDepthForSearch[iD], &listDaughters)) {
            continue;
          }
          if ((iD == XiCzeroToXiPi || iD == OmegaCToXiPi) &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], std::array{+kXiMinus, +kPiPlus}, true)) {
            continue;
          }
          if (iD == OmegaCToOmegaPi &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], std::array{+kOmegaMinus, +kPiPlus}, true)) {
            continue;
          }
        }

        if (nDaughters[iD] == 5) {
          if (!RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, PDGArrayParticle[iD], arrPDGFinal5Prong[iD], true, nullptr, maxDepthForSearch[iD], &listDaughters)) {
            continue;
          }
          if (iD == XiCplusToXiPiPi &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], std::array{+kXiMinus, +kPiPlus, +kPiPlus}, true, nullptr, 2)) {
            continue;
          }
        }

        // Check momentum conservation
        double sumPxDau = 0.;
        double sumPyDau = 0.;
        double sumPzDau = 0.;
        bool momentumCheck = true;
        for (std::size_t iDau = 0; iDau < listDaughters.size(); ++iDau) {
          auto daughter = mcParticles.rawIteratorAt(listDaughters.at(iDau) - mcParticles.offset());
          sumPxDau += daughter.px();
          sumPyDau += daughter.py();
          sumPzDau += daughter.pz();
        }
        auto pxDiff = particle.px() - sumPxDau;
        auto pyDiff = particle.py() - sumPyDau;
        auto pzDiff = particle.pz() - sumPzDau;
        if (std::abs(pxDiff) > 0.001 || std::abs(pyDiff) > 0.001 || std::abs(pzDiff) > 0.001) {
          momentumCheck = false;
        }
        double pDiff = RecoDecay::p(pxDiff, pyDiff, pzDiff);
        double ptDiff = RecoDecay::pt(pxDiff, pyDiff);
        registry.fill(HIST("hMomentumCheck"), static_cast<float>(momentumCheck));
        registry.fill(HIST("hPxDiffMotherDaughterGen"), pxDiff);
        registry.fill(HIST("hPyDiffMotherDaughterGen"), pyDiff);
        registry.fill(HIST("hPzDiffMotherDaughterGen"), pzDiff);
        registry.fill(HIST("hPDiffMotherDaughterGen"), pDiff);
        registry.fill(HIST("hPtDiffMotherDaughterGen"), ptDiff);

        int origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle);
        if (origin == RecoDecay::OriginType::Prompt) {
          counterPrompt[iD]++;
        } else if (origin == RecoDecay::OriginType::NonPrompt) {
          counterNonPrompt[iD]++;
        }

        auto daughter0 = particle.daughters_as<aod::McParticles>().begin();
        double vertexDau[3] = {daughter0.vx(), daughter0.vy(), daughter0.vz()};
        double vertexPrimary[3] = {mcCollision.posX(), mcCollision.posY(), mcCollision.posZ()};
        auto decayLength = RecoDecay::distance(vertexPrimary, vertexDau);
        if (iD < nMesonChannels) {
          if (origin == RecoDecay::OriginType::Prompt) {
            if (std::abs(particle.y()) < 0.5) {
              registry.fill(HIST("PromptCharmMesons/hPromptMesonsPtDistr"), iD, particle.pt());
            }
            registry.fill(HIST("PromptCharmMesons/hPromptMesonsYDistr"), iD, particle.y());
            registry.fill(HIST("PromptCharmMesons/hPromptMesonsDecLenDistr"), iD, decayLength * 10000);
          } else if (origin == RecoDecay::OriginType::NonPrompt) {
            if (std::abs(particle.y()) < 0.5) {
              registry.fill(HIST("NonPromptCharmMesons/hNonPromptMesonsPtDistr"), iD, particle.pt());
            }
            registry.fill(HIST("NonPromptCharmMesons/hNonPromptMesonsYDistr"), iD, particle.y());
            registry.fill(HIST("NonPromptCharmMesons/hNonPromptMesonsDecLenDistr"), iD, decayLength * 10000);
          }
        } else {
          if (origin == RecoDecay::OriginType::Prompt) {
            if (std::abs(particle.y()) < 0.5) {
              registry.fill(HIST("PromptCharmBaryons/hPromptBaryonsPtDistr"), iD - nMesonChannels, particle.pt());
            }
            registry.fill(HIST("PromptCharmBaryons/hPromptBaryonsYDistr"), iD - nMesonChannels, particle.y());
            registry.fill(HIST("PromptCharmBaryons/hPromptBaryonsDecLenDistr"), iD - nMesonChannels, decayLength * 10000);
          } else if (origin == RecoDecay::OriginType::NonPrompt) {
            if (std::abs(particle.y()) < 0.5) {
              registry.fill(HIST("NonPromptCharmBaryons/hNonPromptBaryonsPtDistr"), iD - nMesonChannels, particle.pt());
            }
            registry.fill(HIST("NonPromptCharmBaryons/hNonPromptBaryonsYDistr"), iD - nMesonChannels, particle.y());
            registry.fill(HIST("NonPromptCharmBaryons/hNonPromptBaryonsDecLenDistr"), iD - nMesonChannels, decayLength * 10000);
          }
        }
      }
    } // end particles

    registry.fill(HIST("Quarks/hCountC"), cPerCollision);
    registry.fill(HIST("Quarks/hCountB"), bPerCollision);
    registry.fill(HIST("Quarks/hCountCbar"), cBarPerCollision);
    registry.fill(HIST("Quarks/hCountBbar"), bBarPerCollision);
    static_for<0, nMesonChannels - 1>([&](auto i) {
      constexpr int index = i.value;
      registry.fill(HIST("PromptCharmMesons/hCountPrompt") + HIST(particleNames[index]), counterPrompt[index]);
      registry.fill(HIST("NonPromptCharmMesons/hCountNonPrompt") + HIST(particleNames[index]), counterNonPrompt[index]);
    });
    static_for<nMesonChannels, nChannels - 1>([&](auto i) {
      constexpr int index = i.value;
      registry.fill(HIST("PromptCharmBaryons/hCountPrompt") + HIST(particleNames[index]), counterPrompt[index]);
      registry.fill(HIST("NonPromptCharmBaryons/hCountNonPrompt") + HIST(particleNames[index]), counterNonPrompt[index]);
    });
  };
};

/// Reconstruction Level Validation
///
///   - Gen-Rec Level Momentum Difference per component;
///   - Gen-Rec Level Difference for secondary-vertex coordinates and decay length;
struct HfTaskMcValidationRec {
  SliceCache cache;
  Preslice<aod::Tracks> perCol = aod::track::collisionId;

  Configurable<int> eventGeneratorType{"eventGeneratorType", -1, "If positive, enable event selection using subGeneratorId information. The value indicates which events to keep (0 = MB, 4 = charm triggered, 5 = beauty triggered)"};

  std::array<std::shared_ptr<TH1>, nChannels> histDeltaPt, histDeltaPx, histDeltaPy, histDeltaPz, histDeltaSecondaryVertexX, histDeltaSecondaryVertexY, histDeltaSecondaryVertexZ, histDeltaDecayLength;
  std::array<std::array<std::array<std::shared_ptr<TH1>, 5>, 2>, nChannels> histPtDau, histEtaDau, histImpactParameterDau;
  std::array<std::array<std::shared_ptr<TH1>, 2>, nChannels> histPtReco;
  std::array<std::shared_ptr<THnSparse>, 4> histOriginTracks;
  std::shared_ptr<TH2> histAmbiguousTracks, histTracks;
  std::shared_ptr<TH1> histContributors;

  using HfCand2ProngWithMCRec = soa::Join<aod::HfCand2Prong, aod::HfCand2ProngMcRec>;
  using HfCand3ProngWithMCRec = soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec>;
  using CollisionsWithMCLabels = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::HfSelCollision>;
  using TracksWithSel = soa::Join<aod::TracksWMc, aod::TracksExtra, aod::TrackSelection, aod::TrackCompColls>;

  Partition<TracksWithSel> tracksFilteredGlobalTrackWoDCA = requireGlobalTrackWoDCAInFilter();
  Partition<TracksWithSel> tracksInAcc = requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks);

  AxisSpec axisDeltaMom{2000, -1., 1.};
  AxisSpec axisOrigin{4, -1.5, 2.5};
  AxisSpec axisEta{40, -1., 1.};
  AxisSpec axisPt{50, 0., 10.};
  AxisSpec axisPtD{100, 0., 50.};
  AxisSpec axisDeltaVtx{200, -1, 1.};
  AxisSpec axisDecision{2, -0.5, 1.5};
  AxisSpec axisITShits{8, -0.5, 7.5};
  AxisSpec axisMult{200, 0., 200.};
  AxisSpec axisR{100, 0., 0.5};
  AxisSpec axisSmallNum{20, -0.5, 19.5};

  HistogramRegistry registry{
    "registry",
    {{"histNtracks", "Number of global tracks w/o DCA requirement;#it{N}_{tracks};entries", {HistType::kTH1F, {axisMult}}},
     {"histXvtxReco", "Position of reco PV in #it{X};#it{X}^{reco} (cm);entries", {HistType::kTH1F, {axisDeltaVtx}}},
     {"histYvtxReco", "Position of reco PV in #it{Y};#it{Y}^{reco} (cm);entries", {HistType::kTH1F, {axisDeltaVtx}}},
     {"histZvtxReco", "Position of reco PV in #it{Z};#it{Z}^{reco} (cm);entries", {HistType::kTH1F, {{200, -20, 20.}}}},
     {"histDeltaZvtx", "Residual distribution of PV in #it{Z} as a function of number of contributors;number of contributors;#it{Z}^{reco} - #it{Z}^{gen} (cm);entries", {HistType::kTH2F, {{100, -0.5, 99.5}, {1000, -0.5, 0.5}}}},
     {"TrackToCollChecks/histAmbiguousTrackNumCollisions", "Number of collisions associated to an ambiguous track;number of collisions;entries", {HistType::kTH1F, {{30, -0.5, 29.5}}}},
     {"TrackToCollChecks/histAmbiguousTrackZvtxRMS", "RMS of #it{Z}^{reco} of collisions associated to a track;RMS(#it{Z}^{reco}) (cm);entries", {HistType::kTH1F, {{100, 0., 0.5}}}},
     {"TrackToCollChecks/histFracGoodContributors", "Fraction of PV contributors originating from the correct collision;fraction;entries", {HistType::kTH1F, {{101, 0., 1.01}}}},
     {"TrackToCollChecks/histCollisionsSameBC", "Collisions in same BC;number of contributors collision 1;number of contributors collision 2;#it{R}_{xy} collision 1 (cm);#it{R}_{xy} collision 2 (cm);number of contributors from beauty collision 1;number of contributors from beauty collision 2;", {HistType::kTHnSparseF, {axisMult, axisMult, axisR, axisR, axisSmallNum, axisSmallNum}}}}};

  /// RMS calculation
  /// \param vec  vector of values to compute RMS
  template <typename T>
  T computeRMS(std::vector<T>& vec)
  {
    T sum = std::accumulate(vec.begin(), vec.end(), 0.0);
    T mean = sum / vec.size();

    std::vector<T> diff(vec.size());
    std::transform(vec.begin(), vec.end(), diff.begin(), [mean](T x) { return x - mean; });
    T sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    T stdev = std::sqrt(sq_sum / vec.size());

    return stdev;
  }

  /// Fill MC histograms with decay properties at reconstruction level
  /// \param candidate is candidate
  /// \param mother is mother particle
  /// \param whichHad int indicating charm-hadron and decay channel, see enum DecayChannels
  /// \param whichOrigin int indicating origin: prompt or non-prompt
  template <typename T, typename U>
  void fillHisto(const T& candidate, const U& mother, int whichHad, int whichOrigin)
  {
    histDeltaPt[whichHad]->Fill(candidate.pt() - mother.pt());
    histDeltaPx[whichHad]->Fill(candidate.px() - mother.px());
    histDeltaPy[whichHad]->Fill(candidate.py() - mother.py());
    histDeltaPz[whichHad]->Fill(candidate.pz() - mother.pz());
    // Compare Secondary vertex and decay length with MC
    auto daughter0 = mother.template daughters_as<aod::McParticles>().begin();
    double vertexDau[3] = {daughter0.vx(), daughter0.vy(), daughter0.vz()};
    double vertexMoth[3] = {mother.vx(), mother.vy(), mother.vz()};
    auto decayLength = RecoDecay::distance(vertexMoth, vertexDau);

    histDeltaSecondaryVertexX[whichHad]->Fill(candidate.xSecondaryVertex() - vertexDau[0]);
    histDeltaSecondaryVertexY[whichHad]->Fill(candidate.ySecondaryVertex() - vertexDau[1]);
    histDeltaSecondaryVertexZ[whichHad]->Fill(candidate.zSecondaryVertex() - vertexDau[2]);
    histDeltaDecayLength[whichHad]->Fill(candidate.decayLength() - decayLength);
    std::array<double, 3> momDau0 = {candidate.pxProng0(),
                                     candidate.pyProng0(),
                                     candidate.pzProng0()};
    std::array<double, 3> momDau1 = {candidate.pxProng1(),
                                     candidate.pyProng1(),
                                     candidate.pzProng1()};
    histPtReco[whichHad][whichOrigin]->Fill(candidate.pt());
    histPtDau[whichHad][whichOrigin][0]->Fill(RecoDecay::pt(momDau0));
    histEtaDau[whichHad][whichOrigin][0]->Fill(RecoDecay::eta(momDau0));
    histImpactParameterDau[whichHad][whichOrigin][0]->Fill(candidate.impactParameter0());
    histPtDau[whichHad][whichOrigin][1]->Fill(RecoDecay::pt(momDau1));
    histEtaDau[whichHad][whichOrigin][1]->Fill(RecoDecay::eta(momDau1));
    histImpactParameterDau[whichHad][whichOrigin][1]->Fill(candidate.impactParameter1());
  }

  void init(InitContext&)
  {
    histOriginTracks[0] = registry.add<THnSparse>("TrackToCollChecks/histOriginNonAssociatedTracks", ";origin;#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); is PV contributor; has TOF; number of ITS hits", HistType::kTHnSparseF, {axisOrigin, axisPt, axisEta, axisDeltaVtx, axisDecision, axisDecision, axisITShits});           // tracks not associated to any collision
    histOriginTracks[1] = registry.add<THnSparse>("TrackToCollChecks/histOriginAssociatedTracks", ";origin;#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); is PV contributor; has TOF; number of ITS hits", HistType::kTHnSparseF, {axisOrigin, axisPt, axisEta, axisDeltaVtx, axisDecision, axisDecision, axisITShits});              // tracks associasted to a collision
    histOriginTracks[2] = registry.add<THnSparse>("TrackToCollChecks/histOriginGoodAssociatedTracks", ";origin;#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); is PV contributor; has TOF; number of ITS hits", HistType::kTHnSparseF, {axisOrigin, axisPt, axisEta, axisDeltaVtx, axisDecision, axisDecision, axisITShits});          // tracks associated to the correct collision considering only first reco collision (based on the MC collision index)
    histOriginTracks[3] = registry.add<THnSparse>("TrackToCollChecks/histOriginGoodAssociatedTracksAmbiguous", ";origin;#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); is PV contributor; has TOF; number of ITS hits", HistType::kTHnSparseF, {axisOrigin, axisPt, axisEta, axisDeltaVtx, axisDecision, axisDecision, axisITShits}); // tracks associated to the correct collision considering all ambiguous reco collisions (based on the MC collision index)
    for (std::size_t iHist{0}; iHist < histOriginTracks.size(); ++iHist) {
      histOriginTracks[iHist]->GetAxis(0)->SetBinLabel(1, "no MC particle");
      histOriginTracks[iHist]->GetAxis(0)->SetBinLabel(2, "no quark");
      histOriginTracks[iHist]->GetAxis(0)->SetBinLabel(3, "charm");
      histOriginTracks[iHist]->GetAxis(0)->SetBinLabel(4, "beauty");
    }
    histAmbiguousTracks = registry.add<TH2>("TrackToCollChecks/histAmbiguousTracks", "Tracks that are ambiguous vs. origin;#it{p}_{T}^{reco} (GeV/#it{c});entries", HistType::kTH2F, {axisOrigin, axisPt});
    histTracks = registry.add<TH2>("histTracks", "Tracks vs. origin;#it{p}_{T}^{reco} (GeV/#it{c});entries", HistType::kTH2F, {axisOrigin, axisPt});
    histTracks->GetXaxis()->SetBinLabel(1, "no MC particle");
    histTracks->GetXaxis()->SetBinLabel(2, "no quark");
    histTracks->GetXaxis()->SetBinLabel(3, "charm");
    histTracks->GetXaxis()->SetBinLabel(4, "beauty");
    histAmbiguousTracks->GetXaxis()->SetBinLabel(1, "no MC particle");
    histAmbiguousTracks->GetXaxis()->SetBinLabel(2, "no quark");
    histAmbiguousTracks->GetXaxis()->SetBinLabel(3, "charm");
    histAmbiguousTracks->GetXaxis()->SetBinLabel(4, "beauty");
    for (auto iHad = 0; iHad < nChannels; ++iHad) {
      histDeltaPt[iHad] = registry.add<TH1>(Form("%s/histDeltaPt", particleNames[iHad].data()), Form("Pt difference reco - MC %s; #it{p}_{T}^{reco} - #it{p}_{T}^{gen} (GeV/#it{c}); entries", labels[iHad].data()), HistType::kTH1F, {axisDeltaMom});
      histDeltaPx[iHad] = registry.add<TH1>(Form("%s/histDeltaPx", particleNames[iHad].data()), Form("Px difference reco - MC %s; #it{p}_{x}^{reco} - #it{p}_{x}^{gen} (GeV/#it{c}); entries", labels[iHad].data()), HistType::kTH1F, {axisDeltaMom});
      histDeltaPy[iHad] = registry.add<TH1>(Form("%s/histDeltaPy", particleNames[iHad].data()), Form("Py difference reco - MC %s; #it{p}_{y}^{reco} - #it{p}_{y}^{gen} (GeV/#it{c}); entries", labels[iHad].data()), HistType::kTH1F, {axisDeltaMom});
      histDeltaPz[iHad] = registry.add<TH1>(Form("%s/histDeltaPz", particleNames[iHad].data()), Form("Pz difference reco - MC %s; #it{p}_{z}^{reco} - #it{p}_{z}^{gen} (GeV/#it{c}); entries", labels[iHad].data()), HistType::kTH1F, {axisDeltaMom});
      histDeltaSecondaryVertexX[iHad] = registry.add<TH1>(Form("%s/histDeltaSecondaryVertexX", particleNames[iHad].data()), Form("Sec. Vertex difference reco - MC (MC matched) - %s; #Delta x (cm); entries", labels[iHad].data()), HistType::kTH1F, {axisDeltaVtx});
      histDeltaSecondaryVertexY[iHad] = registry.add<TH1>(Form("%s/histDeltaSecondaryVertexY", particleNames[iHad].data()), Form("Sec. Vertex difference reco - MC (MC matched) - %s; #Delta y (cm); entries", labels[iHad].data()), HistType::kTH1F, {axisDeltaVtx});
      histDeltaSecondaryVertexZ[iHad] = registry.add<TH1>(Form("%s/histDeltaSecondaryVertexZ", particleNames[iHad].data()), Form("Sec. Vertex difference reco - MC (MC matched) - %s; #Delta z (cm); entries", labels[iHad].data()), HistType::kTH1F, {axisDeltaVtx});
      histDeltaDecayLength[iHad] = registry.add<TH1>(Form("%s/histDeltaDecayLength", particleNames[iHad].data()), Form("Decay length difference reco - MC (%s); #Delta L (cm); entries", labels[iHad].data()), HistType::kTH1F, {axisDeltaVtx});
      for (auto iOrigin = 0; iOrigin < 2; ++iOrigin) {
        histPtReco[iHad][iOrigin] = registry.add<TH1>(Form("%s/histPtReco%s", particleNames[iHad].data(), originNames[iOrigin].data()), Form("Pt reco %s %s; #it{p}_{T}^{reco} (GeV/#it{c}); entries", originNames[iOrigin].data(), labels[iHad].data()), HistType::kTH1F, {axisPtD});
        for (unsigned int iDau = 0; iDau < nDaughters[iHad]; ++iDau) {
          histPtDau[iHad][iOrigin][iDau] = registry.add<TH1>(Form("%s/histPtDau%d%s", particleNames[iHad].data(), iDau, originNames[iOrigin].data()), Form("Daughter %d Pt reco - %s %s; #it{p}_{T}^{dau, reco} (GeV/#it{c}); entries", iDau, originNames[iOrigin].data(), labels[iHad].data()), HistType::kTH1F, {axisPt});
          histEtaDau[iHad][iOrigin][iDau] = registry.add<TH1>(Form("%s/histEtaDau%d%s", particleNames[iHad].data(), iDau, originNames[iOrigin].data()), Form("Daughter %d Eta reco - %s %s; #it{#eta}^{dau, reco}; entries", iDau, originNames[iOrigin].data(), labels[iHad].data()), HistType::kTH1F, {{100, -1., 1.}});
          histImpactParameterDau[iHad][iOrigin][iDau] = registry.add<TH1>(Form("%s/histImpactParameterDau%d%s", particleNames[iHad].data(), iDau, originNames[iOrigin].data()), Form("Daughter %d DCAxy reco - %s %s; DCAxy (cm); entries", iDau, originNames[iOrigin].data(), labels[iHad].data()), HistType::kTH1F, {axisDeltaVtx});
        }
      }
    }
    histContributors = registry.add<TH1>("TrackToCollChecks/histContributors", "PV contributors from correct/wrong MC collision;;entries", HistType::kTH1F, {axisDecision});
    histContributors->GetXaxis()->SetBinLabel(1, "correct MC collision");
    histContributors->GetXaxis()->SetBinLabel(2, "wrong MC collision");
  }

  void process(CollisionsWithMCLabels::iterator const& collision,
               aod::McCollisions const&)
  {
    // check that collision is selected by hf-track-index-skim-creator-tag-sel-collisions
    if (collision.whyRejectColl() != 0) {
      return;
    }
    if (!collision.has_mcCollision()) {
      return;
    }
    auto mcCollision = collision.mcCollision_as<aod::McCollisions>();
    if (eventGeneratorType >= 0 && mcCollision.getSubGeneratorId() != eventGeneratorType) {
      return;
    }

    registry.fill(HIST("histXvtxReco"), collision.posX());
    registry.fill(HIST("histYvtxReco"), collision.posY());
    registry.fill(HIST("histZvtxReco"), collision.posZ());
    registry.fill(HIST("histDeltaZvtx"), collision.numContrib(), collision.posZ() - mcCollision.posZ());
  } // end process

  void processCollAssoc(CollisionsWithMCLabels const& collisions,
                        TracksWithSel const&,
                        aod::McParticles const& mcParticles,
                        aod::McCollisions const&,
                        aod::BCs const&)
  {
    // loop over collisions
    for (auto collision = collisions.begin(); collision != collisions.end(); ++collision) {
      // check that collision is selected by hf-track-index-skim-creator-tag-sel-collisions
      if (collision.whyRejectColl() != 0) {
        continue;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto mcCollision = collision.mcCollision_as<aod::McCollisions>();
      if (eventGeneratorType >= 0 && mcCollision.getSubGeneratorId() != eventGeneratorType) {
        continue;
      }
      auto tracksGlobalWoDCAColl1 = tracksFilteredGlobalTrackWoDCA->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      registry.fill(HIST("histNtracks"), tracksGlobalWoDCAColl1.size());
      auto tracksColl1 = tracksInAcc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      int nContributors = 0, nGoodContributors = 0;
      for (const auto& track : tracksColl1) {
        if (!track.isPVContributor()) {
          continue;
        }
        if (!track.has_mcParticle()) {
          continue;
        }
        nContributors++;
        auto particle = track.mcParticle();
        if (collision.mcCollisionId() == particle.mcCollisionId()) {
          nGoodContributors++;
        }
      }
      float frac = (nContributors > 0) ? static_cast<float>(nGoodContributors) / nContributors : 1.;
      registry.fill(HIST("TrackToCollChecks/histFracGoodContributors"), frac);
      uint64_t mostProbableBC = collision.bc().globalBC();
      for (auto collision2 = collision + 1; collision2 != collisions.end(); ++collision2) {
        uint64_t mostProbableBC2 = collision2.bc().globalBC();
        if (mostProbableBC2 == mostProbableBC) {
          float radColl1 = std::sqrt(collision.posX() * collision.posX() + collision.posY() * collision.posY());
          float radColl2 = std::sqrt(collision2.posX() * collision2.posX() + collision2.posY() * collision2.posY());
          int nFromBeautyColl1 = 0, nFromBeautyColl2 = 0;
          for (const auto& trackColl1 : tracksColl1) {
            if (trackColl1.has_mcParticle() && trackColl1.isPVContributor()) {
              auto particleColl1 = trackColl1.mcParticle();
              auto origin = RecoDecay::getCharmHadronOrigin(mcParticles, particleColl1, true);
              if (origin == RecoDecay::NonPrompt) {
                nFromBeautyColl1++;
              }
            }
          }
          auto tracksColl2 = tracksInAcc->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);
          for (const auto& trackColl2 : tracksColl2) {
            if (trackColl2.has_mcParticle() && trackColl2.isPVContributor()) {
              auto particleColl2 = trackColl2.mcParticle();
              auto origin = RecoDecay::getCharmHadronOrigin(mcParticles, particleColl2, true);
              if (origin == RecoDecay::NonPrompt) {
                nFromBeautyColl2++;
              }
            }
          }
          registry.fill(HIST("TrackToCollChecks/histCollisionsSameBC"), collision.numContrib(), collision2.numContrib(), radColl1, radColl2, nFromBeautyColl1, nFromBeautyColl2);
          break;
        }
      }
    }

    // loop over tracks
    for (const auto& track : tracksFilteredGlobalTrackWoDCA) {
      // check number of ITS hits
      int nITSlayers = 0;
      uint8_t ITSHitMap = track.itsClusterMap();
      for (int iLayer = 0; iLayer < 7; ++iLayer) {
        if (TESTBIT(ITSHitMap, iLayer)) {
          nITSlayers++;
        }
      }
      uint index = uint(track.collisionId() >= 0);
      if (track.has_mcParticle()) {
        auto particle = track.mcParticle(); // get corresponding MC particle to check origin
        auto mcCollision = particle.mcCollision_as<aod::McCollisions>();
        if (eventGeneratorType >= 0 && mcCollision.getSubGeneratorId() != eventGeneratorType) {
          continue;
        }
        auto origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, true);
        histTracks->Fill(origin, track.pt());
        bool isAmbiguous = (track.compatibleCollIds().size() != 1);
        if (isAmbiguous) {
          registry.fill(HIST("TrackToCollChecks/histAmbiguousTrackNumCollisions"), track.compatibleCollIds().size());
          histAmbiguousTracks->Fill(origin, track.pt());
          std::vector<double> ambCollPosZ{};
          for (const auto& collIdx : track.compatibleCollIds()) {
            auto ambCollision = collisions.rawIteratorAt(collIdx);
            ambCollPosZ.push_back(ambCollision.posZ());
          }
          // here we are only interested to tracks associated to multiple vertices
          if (ambCollPosZ.size() > 0) {
            registry.fill(HIST("TrackToCollChecks/histAmbiguousTrackZvtxRMS"), computeRMS(ambCollPosZ));
          }
        }
        float deltaZ = -999.f;
        if (index) {
          auto collision = track.collision_as<CollisionsWithMCLabels>();
          auto mcCollision = particle.mcCollision_as<aod::McCollisions>();
          deltaZ = collision.posZ() - mcCollision.posZ();
          if (collision.has_mcCollision() && collision.mcCollisionId() == particle.mcCollisionId()) {
            histOriginTracks[index + 1]->Fill(origin, track.pt(), track.eta(), deltaZ, track.isPVContributor(), track.hasTOF(), nITSlayers);
          } else { // if the default associated collision is not the good one, check if the tracks is ambiguous
            if (isAmbiguous) {
              for (const auto& collIdx : track.compatibleCollIds()) {
                auto ambCollision = collisions.rawIteratorAt(collIdx);

                if (ambCollision.has_mcCollision() && ambCollision.mcCollisionId() == particle.mcCollisionId()) {
                  histOriginTracks[index + 2]->Fill(origin, track.pt(), track.eta(), deltaZ, track.isPVContributor(), track.hasTOF(), nITSlayers);
                  break;
                }
              }
            }
          }
          if (track.isPVContributor()) {
            if (collision.has_mcCollision() && collision.mcCollisionId() == particle.mcCollisionId()) {
              histContributors->Fill(0);
            } else {
              histContributors->Fill(1);
            }
          }
        }
        histOriginTracks[index]->Fill(origin, track.pt(), track.eta(), deltaZ, track.isPVContributor(), track.hasTOF(), nITSlayers);
      } else {
        histOriginTracks[index]->Fill(-1.f, track.pt(), track.eta(), -999.f, track.isPVContributor(), track.hasTOF(), nITSlayers);
      }
    }
  }
  PROCESS_SWITCH(HfTaskMcValidationRec, processCollAssoc, "Process collision-association information, requires extra table from TrackToCollisionAssociation task (fillTableOfCollIdsPerTrack=true)", false);

  void processEff(HfCand2ProngWithMCRec const& cand2Prongs,
                  HfCand3ProngWithMCRec const& cand3Prongs,
                  aod::TracksWMc const&,
                  aod::McParticles const& mcParticles,
                  aod::McCollisions const&,
                  CollisionsWithMCLabels const& /*collisions*/)
  {
    // loop over 2-prong candidates
    for (const auto& cand2Prong : cand2Prongs) {

      if (cand2Prong.collision_as<CollisionsWithMCLabels>().has_mcCollision()) {
        auto mcCollision = cand2Prong.collision_as<CollisionsWithMCLabels>().mcCollision_as<aod::McCollisions>();
        if (eventGeneratorType >= 0 && mcCollision.getSubGeneratorId() != eventGeneratorType) {
          continue;
        }
      }

      // determine which kind of candidate it is
      bool isD0Sel = TESTBIT(cand2Prong.hfflag(), o2::aod::hf_cand_2prong::DecayType::D0ToPiK);
      if (!isD0Sel) {
        continue;
      }
      int whichHad = -1;
      if (isD0Sel && TESTBIT(std::abs(cand2Prong.flagMcMatchRec()), hf_cand_2prong::DecayType::D0ToPiK)) {
        whichHad = DzeroToKPi;
      }
      int whichOrigin;
      if (cand2Prong.originMcRec() == RecoDecay::OriginType::Prompt) {
        whichOrigin = 0;
      } else {
        whichOrigin = 1;
      }

      if (whichHad >= 0) {
        int indexParticle = -1;
        if (cand2Prong.prong0_as<aod::TracksWMc>().has_mcParticle()) {
          indexParticle = RecoDecay::getMother(mcParticles, cand2Prong.prong0_as<aod::TracksWMc>().mcParticle(), PDGArrayParticle[whichHad], true);
        }
        if (indexParticle < 0) {
          continue;
        }
        auto mother = mcParticles.rawIteratorAt(indexParticle);
        fillHisto(cand2Prong, mother, whichHad, whichOrigin);
      }
    } // end loop on 2-prong candidates

    // loop over 3-prong candidates
    for (const auto& cand3Prong : cand3Prongs) {

      if (cand3Prong.collision_as<CollisionsWithMCLabels>().has_mcCollision()) {
        auto mcCollision = cand3Prong.collision_as<CollisionsWithMCLabels>().mcCollision_as<aod::McCollisions>();
        if (eventGeneratorType >= 0 && mcCollision.getSubGeneratorId() != eventGeneratorType) {
          continue;
        }
      }

      // determine which kind of candidate it is
      // FIXME: add D* and decays with cascades
      bool isDPlusSel = TESTBIT(cand3Prong.hfflag(), hf_cand_3prong::DecayType::DplusToPiKPi);
      bool isDsSel = TESTBIT(cand3Prong.hfflag(), hf_cand_3prong::DecayType::DsToKKPi);
      bool isLcSel = TESTBIT(cand3Prong.hfflag(), hf_cand_3prong::DecayType::LcToPKPi);
      bool isXicSel = TESTBIT(cand3Prong.hfflag(), hf_cand_3prong::DecayType::XicToPKPi);
      if (!isDPlusSel && !isDsSel && !isLcSel && !isXicSel) {
        continue;
      }
      int whichHad = -1;
      if (isDPlusSel && TESTBIT(std::abs(cand3Prong.flagMcMatchRec()), hf_cand_3prong::DecayType::DplusToPiKPi)) {
        whichHad = DplusToPiKPi;
      } else if (isDsSel && TESTBIT(std::abs(cand3Prong.flagMcMatchRec()), hf_cand_3prong::DecayType::DsToKKPi)) {
        if (cand3Prong.flagMcDecayChanRec() == hf_cand_3prong::DecayChannelDToKKPi::DsToPhiPi) {
          whichHad = DsToPhiPiToKKPi;
        }
        if (cand3Prong.flagMcDecayChanRec() == hf_cand_3prong::DecayChannelDToKKPi::DsToK0starK) {
          whichHad = DsToK0starKToKKPi;
        }
        if (cand3Prong.flagMcDecayChanRec() == hf_cand_3prong::DecayChannelDToKKPi::DplusToPhiPi) {
          whichHad = DplusToPhiPiToKKPi;
        }
      } else if (isLcSel && TESTBIT(std::abs(cand3Prong.flagMcMatchRec()), hf_cand_3prong::DecayType::LcToPKPi)) {
        whichHad = LcToPKPi;
      } else if (isXicSel && TESTBIT(std::abs(cand3Prong.flagMcMatchRec()), hf_cand_3prong::DecayType::XicToPKPi)) {
        whichHad = XiCplusToPKPi;
      }
      int whichOrigin;
      if (cand3Prong.originMcRec() == RecoDecay::OriginType::Prompt) {
        whichOrigin = 0;
      } else {
        whichOrigin = 1;
      }

      if (whichHad >= 0) {
        int indexParticle = -1;
        if (cand3Prong.prong0_as<aod::TracksWMc>().has_mcParticle()) {
          indexParticle = RecoDecay::getMother(mcParticles, cand3Prong.prong0_as<aod::TracksWMc>().mcParticle(), PDGArrayParticle[whichHad], true);
        }
        if (indexParticle < 0) {
          continue;
        }
        auto mother = mcParticles.rawIteratorAt(indexParticle);
        fillHisto(cand3Prong, mother, whichHad, whichOrigin);
        std::array<double, 3> momDau2 = {cand3Prong.pxProng2(),
                                         cand3Prong.pyProng2(),
                                         cand3Prong.pzProng2()};
        histPtDau[whichHad][whichOrigin][2]->Fill(RecoDecay::pt(momDau2));
        histEtaDau[whichHad][whichOrigin][2]->Fill(RecoDecay::eta(momDau2));
        histImpactParameterDau[whichHad][whichOrigin][2]->Fill(cand3Prong.impactParameter2());
      }
    } // end loop on 3-prong candidates
  }
  PROCESS_SWITCH(HfTaskMcValidationRec, processEff, "Compute charm-hadron efficiencies (not all of them are implemented), requires HF candidate creators", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HfTaskMcValidationGen>(cfgc),
    adaptAnalysisTask<HfTaskMcValidationRec>(cfgc)};
  return workflow;
}
