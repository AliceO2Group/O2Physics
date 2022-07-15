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

/// \file HFMCValidation.cxx
/// \brief Monte Carlo validation task
/// \note gen and rec. level validation
///
/// \author Antonio Palasciano <antonio.palasciano@cern.ch>, Università degli Studi di Bari & INFN, Sezione di Bari
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::aod::hf_cand_prong3;

namespace
{
static const int nCharmHadrons = 7;
static const std::array<int, nCharmHadrons> PDGArrayParticle = {pdg::Code::kDPlus, 413, pdg::Code::kD0, 431, pdg::Code::kLambdaCPlus, pdg::Code::kXiCPlus, pdg::Code::kJpsi};
static const std::array<int, nCharmHadrons> nDaughters = {3, 3, 2, 3, 3, 3, 2};
static const std::array<std::string, nCharmHadrons> labels = {"D^{#plus}", "D*^{#plus}", "D^{0}", "D_{s}^{#plus}", "#Lambda_{c}^{#plus} #rightarrow pK^{#minus}#pi^{#plus}", "#Xi_{c}^{#plus} #rightarrow pK^{#minus}#pi^{#plus}", "J/#psi #rightarrow e^{#plus}e^{#minus}"};
static const std::array<std::string, nCharmHadrons> particleNames = {"Dplus", "Dstar", "D0", "Ds", "Lc2pKpi", "Xic2pKpi", "Jpsi2ee"};
static const std::array<std::string, 2> originNames = {"Prompt", "NonPrompt"};
} // namespace

/// Gen Level Validation
///
/// - Number of HF quarks produced per collision
/// - Number of D±      → π± K∓ π±        per collision
///             D*±     → π± K∓ π±,
///             D0(bar) → π± K∓,
///             Ds±     → K± K∓ π±,
///             Λc±     → p(bar) K∓ π±
///             Ξc±     → p(bar) K∓ π±
///             J/psi   → e∓ e±
/// - Momentum Conservation for these particles

struct ValidationGenLevel {
  std::shared_ptr<TH1> hPromptCharmHadronsPtDistr, hPromptCharmHadronsYDistr, hNonPromptCharmHadronsPtDistr, hNonPromptCharmHadronsYDistr, hPromptCharmHadronsDecLenDistr, hNonPromptCharmHadronsDecLenDistr, hQuarkPerEvent;

  HistogramRegistry registry{
    "registry",
    {{"hMomentumCheck", "Mom. Conservation (1 = true, 0 = false) (#it{#epsilon} = 1 MeV/#it{c}); Mom. Conservation result; entries", {HistType::kTH1F, {{2, -0.5, +1.5}}}},
     {"hPtDiffMotherDaughterGen", "Pt Difference Mother-Daughters; #Delta#it{p}_{T}^{gen} (GeV/#it{c}); entries", {HistType::kTH1F, {{100, -0.01, 0.01}}}},
     {"hPxDiffMotherDaughterGen", "Px Difference Mother-Daughters; #Delta#it{p}_{x}^{gen} (GeV/#it{c}); entries", {HistType::kTH1F, {{100, -0.01, 0.01}}}},
     {"hPyDiffMotherDaughterGen", "Py Difference Mother-Daughters; #Delta#it{p}_{y}^{gen} (GeV/#it{c}); entries", {HistType::kTH1F, {{100, -0.01, 0.01}}}},
     {"hPzDiffMotherDaughterGen", "Pz Difference Mother-Daughters; #Delta#it{p}_{z}^{gen} (GeV/#it{c}); entries", {HistType::kTH1F, {{100, -0.01, 0.01}}}},
     {"hPDiffMotherDaughterGen", "P  Difference Mother-Daughters; #Delta#it{p}^{gen} (GeV/#it{c}); entries", {HistType::kTH1F, {{100, -0.01, 0.01}}}},
     {"hCountAverageC", "Event counter - Average Number Charm quark; Events Per Collision; entries", {HistType::kTH1F, {{20, 0., 20.}}}},
     {"hCountAverageB", "Event counter - Average Number Beauty quark; Events Per Collision; entries", {HistType::kTH1F, {{20, 0., 20.}}}},
     {"hCountAverageCbar", "Event counter - Average Number Anti-Charm quark; Events Per Collision; entries", {HistType::kTH1F, {{20, 0., 20.}}}},
     {"hCountAverageBbar", "Event counter - Average Number Anti-Beauty quark; Events Per Collision; entries", {HistType::kTH1F, {{20, 0., 20.}}}},
     {"hCounterPerCollisionPromptDzero", "Event counter - prompt D0; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hCounterPerCollisionPromptDplus", "Event counter - prompt DPlus; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hCounterPerCollisionPromptDs", "Event counter - prompt Ds; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hCounterPerCollisionPromptDstar", "Event counter - prompt Dstar; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hCounterPerCollisionPromptLambdaC", "Event counter - prompt LambdaC; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hCounterPerCollisionPromptXiC", "Event counter - prompt XiC; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hCounterPerCollisionPromptJPsi", "Event counter - prompt JPsi; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hCounterPerCollisionNonPromptDzero", "Event counter - non-prompt D0; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hCounterPerCollisionNonPromptDplus", "Event counter - non-prompt DPlus; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hCounterPerCollisionNonPromptDs", "Event counter - non-prompt Ds; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hCounterPerCollisionNonPromptDstar", "Event counter - non-prompt Dstar; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hCounterPerCollisionNonPromptLambdaC", "Event counter - non-prompt LambdaC; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hCounterPerCollisionNonPromptXiC", "Event counter - non-prompt XiC; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hCounterPerCollisionNonPromptJPsi", "Event counter - non-prompt JPsi; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hPtVsYCharmQuark", "Y vs. Pt - charm quarks ; #it{p}_{T}^{gen} (GeV/#it{c}); #it{y}^{gen}", {HistType::kTH2F, {{100, 0., 50.}, {100, -5., 5.}}}},
     {"hPtVsYBeautyQuark", "Y vs. Pt - beauty quarks ; #it{p}_{T}^{gen} (GeV/#it{c}); #it{y}^{gen}", {HistType::kTH2F, {{100, 0., 50.}, {100, -5., 5.}}}}}};

  void init(o2::framework::InitContext&)
  {
    hPromptCharmHadronsPtDistr = registry.add<TH2>("hPromptCharmHadronsPtDistr", "Pt distribution vs prompt charm hadron in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c})", HistType::kTH2F, {{7, -0.5, 6.5}, {100, 0., 50.}});
    hPromptCharmHadronsYDistr = registry.add<TH2>("hPromptCharmHadronsYDistr", "Y distribution vs prompt charm hadron; ; #it{y}^{gen}", HistType::kTH2F, {{7, -0.5, 6.5}, {100, -5., 5.}});
    hPromptCharmHadronsDecLenDistr = registry.add<TH2>("hPromptCharmHadronsDecLDistr", "Decay length distribution vs prompt charm hadron; ; decay length (#mum)", HistType::kTH2F, {{7, -0.5, 6.5}, {100, 0., 10000.}});
    hNonPromptCharmHadronsPtDistr = registry.add<TH2>("hNonPromptCharmHadronsPtDistr", "Pt distribution vs non-prompt charm hadron in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c})", HistType::kTH2F, {{7, -0.5, 6.5}, {100, 0., 50.}});
    hNonPromptCharmHadronsYDistr = registry.add<TH2>("hNonPromptCharmHadronsYDistr", "Y distribution vs non-prompt charm hadron; ; #it{y}^{gen}", HistType::kTH2F, {{7, -0.5, 6.5}, {100, -5., 5.}});
    hNonPromptCharmHadronsDecLenDistr = registry.add<TH2>("hNonPromptCharmHadronsDecLenDistr", "Decay length distribution vs non-prompt charm hadron; ; decay length (#mum)", HistType::kTH2F, {{7, -0.5, 6.5}, {100, 0., 10000.}});
    for (auto iBin = 1; iBin <= nCharmHadrons; ++iBin) {
      hPromptCharmHadronsPtDistr->GetXaxis()->SetBinLabel(iBin, labels[iBin - 1].data());
      hPromptCharmHadronsYDistr->GetXaxis()->SetBinLabel(iBin, labels[iBin - 1].data());
      hPromptCharmHadronsDecLenDistr->GetXaxis()->SetBinLabel(iBin, labels[iBin - 1].data());
      hNonPromptCharmHadronsPtDistr->GetXaxis()->SetBinLabel(iBin, labels[iBin - 1].data());
      hNonPromptCharmHadronsYDistr->GetXaxis()->SetBinLabel(iBin, labels[iBin - 1].data());
      hNonPromptCharmHadronsDecLenDistr->GetXaxis()->SetBinLabel(iBin, labels[iBin - 1].data());
    }
  }

  void process(aod::McCollision const& mccollision, aod::McParticles const& particlesMC)
  {
    int cPerCollision = 0;
    int cBarPerCollision = 0;
    int bPerCollision = 0;
    int bBarPerCollision = 0;

    // Particles and their decay checked in the second part of the task
    std::array<int, nCharmHadrons> PDGArrayParticle = {pdg::Code::kDPlus, 413, pdg::Code::kD0, 431, pdg::Code::kLambdaCPlus, pdg::Code::kXiCPlus, pdg::Code::kJpsi};
    std::array<std::array<int, 3>, nCharmHadrons> arrPDGFinal = {{{kPiPlus, kPiPlus, -kKPlus}, {kPiPlus, kPiPlus, -kKPlus}, {-kKPlus, kPiPlus, 0}, {kPiPlus, kKPlus, -kKPlus}, {kProton, -kKPlus, kPiPlus}, {kProton, -kKPlus, kPiPlus}, {kElectron, -kElectron, 0}}};
    std::array<int, nCharmHadrons> counterPrompt{0}, counterNonPrompt{0};
    std::vector<int> listDaughters{};

    for (auto& particle : particlesMC) {
      int particlePdgCode = particle.pdgCode();
      if (!particle.has_mothers()) {
        continue;
      }
      auto mother = particle.mothers_as<aod::McParticles>().front();
      if (particlePdgCode != mother.pdgCode()) {
        switch (particlePdgCode) {
          case kCharm:
            cPerCollision++;
            registry.fill(HIST("hPtVsYCharmQuark"), particle.pt(), particle.y());
            break;
          case kCharmBar:
            cBarPerCollision++;
            registry.fill(HIST("hPtVsYCharmQuark"), particle.pt(), particle.y());
            break;
          case kBottom:
            bPerCollision++;
            registry.fill(HIST("hPtVsYBeautyQuark"), particle.pt(), particle.y());
            break;
          case kBottomBar:
            bBarPerCollision++;
            registry.fill(HIST("hPtVsYBeautyQuark"), particle.pt(), particle.y());
            break;
        }
      }

      double sumPxDau = 0.;
      double sumPyDau = 0.;
      double sumPzDau = 0.;
      bool momentumCheck = true;
      listDaughters.clear();

      // Checking the decay of the particles and the momentum conservation
      for (std::size_t iD = 0; iD < PDGArrayParticle.size(); ++iD) {
        int whichHadron = -1;
        if (std::abs(particlePdgCode) == PDGArrayParticle[iD]) {
          whichHadron = iD;
          RecoDecay::getDaughters(particle, &listDaughters, arrPDGFinal[iD], -1);
          std::size_t arrayPDGsize = arrPDGFinal[iD].size() - std::count(arrPDGFinal[iD].begin(), arrPDGFinal[iD].end(), 0);
          int origin = -1;
          if (listDaughters.size() == arrayPDGsize) {
            origin = RecoDecay::getCharmHadronOrigin(particlesMC, particle);
            if (origin == RecoDecay::OriginType::Prompt) {
              counterPrompt[iD]++;
            } else if (origin == RecoDecay::OriginType::NonPrompt) {
              counterNonPrompt[iD]++;
            }
          }
          for (std::size_t iDau = 0; iDau < listDaughters.size(); ++iDau) {
            auto daughter = particlesMC.rawIteratorAt(listDaughters.at(iDau) - particlesMC.offset());
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
          auto daughter0 = particle.daughters_as<aod::McParticles>().begin();
          double vertexDau[3] = {daughter0.vx(), daughter0.vy(), daughter0.vz()};
          double vertexPrimary[3] = {mccollision.posX(), mccollision.posY(), mccollision.posZ()};

          auto decayLength = RecoDecay::distance(vertexPrimary, vertexDau);
          // Filling histograms with per-component momentum conservation
          registry.fill(HIST("hMomentumCheck"), float(momentumCheck));
          registry.fill(HIST("hPxDiffMotherDaughterGen"), pxDiff);
          registry.fill(HIST("hPyDiffMotherDaughterGen"), pyDiff);
          registry.fill(HIST("hPzDiffMotherDaughterGen"), pzDiff);
          registry.fill(HIST("hPDiffMotherDaughterGen"), pDiff);
          registry.fill(HIST("hPtDiffMotherDaughterGen"), ptDiff);
          if (origin == RecoDecay::OriginType::Prompt) {
            if (std::abs(particle.y()) < 0.5) {
              hPromptCharmHadronsPtDistr->Fill(whichHadron, particle.pt());
            }
            hPromptCharmHadronsYDistr->Fill(whichHadron, particle.y());
            hPromptCharmHadronsDecLenDistr->Fill(whichHadron, decayLength * 10000);
          } else if (origin == RecoDecay::OriginType::NonPrompt) {
            if (std::abs(particle.y()) < 0.5) {
              hNonPromptCharmHadronsPtDistr->Fill(whichHadron, particle.pt());
            }
            hNonPromptCharmHadronsYDistr->Fill(whichHadron, particle.y());
            hNonPromptCharmHadronsDecLenDistr->Fill(whichHadron, decayLength * 10000);
          }
        }
      }
    } // end particles
    registry.fill(HIST("hCountAverageC"), cPerCollision);
    registry.fill(HIST("hCountAverageB"), bPerCollision);
    registry.fill(HIST("hCountAverageCbar"), cBarPerCollision);
    registry.fill(HIST("hCountAverageBbar"), bBarPerCollision);
    registry.fill(HIST("hCounterPerCollisionPromptDplus"), counterPrompt[0]);
    registry.fill(HIST("hCounterPerCollisionPromptDstar"), counterPrompt[1]);
    registry.fill(HIST("hCounterPerCollisionPromptDzero"), counterPrompt[2]);
    registry.fill(HIST("hCounterPerCollisionPromptDs"), counterPrompt[3]);
    registry.fill(HIST("hCounterPerCollisionPromptLambdaC"), counterPrompt[4]);
    registry.fill(HIST("hCounterPerCollisionPromptXiC"), counterPrompt[5]);
    registry.fill(HIST("hCounterPerCollisionPromptJPsi"), counterPrompt[6]);
    registry.fill(HIST("hCounterPerCollisionNonPromptDplus"), counterNonPrompt[0]);
    registry.fill(HIST("hCounterPerCollisionNonPromptDstar"), counterNonPrompt[1]);
    registry.fill(HIST("hCounterPerCollisionNonPromptDzero"), counterNonPrompt[2]);
    registry.fill(HIST("hCounterPerCollisionNonPromptDs"), counterNonPrompt[3]);
    registry.fill(HIST("hCounterPerCollisionNonPromptLambdaC"), counterNonPrompt[4]);
    registry.fill(HIST("hCounterPerCollisionNonPromptXiC"), counterNonPrompt[5]);
    registry.fill(HIST("hCounterPerCollisionNonPromptJPsi"), counterNonPrompt[6]);
  }
};

/// Rec Level Validation
///
/// D±      → π± K∓ π±
/// Ds±     → K± K∓ π±,
/// D0(bar) → π± K∓,
/// Λc±     → p(bar) K∓ π±
/// Ξc±     → p(bar) K∓ π±
/// J/psi   → e∓ e±
///   - Gen-Rec Level Momentum Difference per component;
///   - Gen-Rec Level Difference for secondary Vertex coordinates and decay length;
struct ValidationRecLevel {

  static const int nCharmHadrons = 7;
  std::array<std::shared_ptr<TH1>, nCharmHadrons> histDeltaPt, histDeltaPx, histDeltaPy, histDeltaPz, histDeltaSecondaryVertexX, histDeltaSecondaryVertexY, histDeltaSecondaryVertexZ, histDeltaDecayLength;
  std::array<std::array<std::array<std::shared_ptr<TH1>, 3>, 2>, nCharmHadrons> histPtDau, histEtaDau, histImpactParameterDau;
  std::array<std::array<std::shared_ptr<TH1>, 2>, nCharmHadrons> histPtReco;
  std::array<std::shared_ptr<TH1>, 2> histOriginTracks;

  HistogramRegistry registry{"registry", {}};
  void init(o2::framework::InitContext&)
  {
    histOriginTracks[0] = registry.add<TH1>("histOriginNonAssociatedTracks", ";origin;entries", HistType::kTH1F, {{4, -1.5, 2.5}});
    histOriginTracks[1] = registry.add<TH1>("histOriginAssociatedTracks", ";origin;entries", HistType::kTH1F, {{4, -1.5, 2.5}});
    for (std::size_t iHist{0}; iHist < histOriginTracks.size(); ++iHist) {
      histOriginTracks[iHist]->GetXaxis()->SetBinLabel(1, "no MC particle");
      histOriginTracks[iHist]->GetXaxis()->SetBinLabel(2, "no quark");
      histOriginTracks[iHist]->GetXaxis()->SetBinLabel(3, "charm");
      histOriginTracks[iHist]->GetXaxis()->SetBinLabel(4, "beauty");
    }
    for (auto iHad = 0; iHad < nCharmHadrons; ++iHad) {
      histDeltaPt[iHad] = registry.add<TH1>(Form("histDeltaPt%s", particleNames[iHad].data()), Form("Pt difference reco - MC %s; #it{p}_{T}^{reco} - #it{p}_{T}^{gen} (GeV/#it{c}); entries", labels[iHad].data()), HistType::kTH1F, {{2000, -1., 1.}});
      histDeltaPx[iHad] = registry.add<TH1>(Form("histDeltaPx%s", particleNames[iHad].data()), Form("Px difference reco - MC %s; #it{p}_{x}^{reco} - #it{p}_{x}^{gen} (GeV/#it{c}); entries", labels[iHad].data()), HistType::kTH1F, {{2000, -1., 1.}});
      histDeltaPy[iHad] = registry.add<TH1>(Form("histDeltaPy%s", particleNames[iHad].data()), Form("Py difference reco - MC %s; #it{p}_{y}^{reco} - #it{p}_{y}^{gen} (GeV/#it{c}); entries", labels[iHad].data()), HistType::kTH1F, {{2000, -1., 1.}});
      histDeltaPz[iHad] = registry.add<TH1>(Form("histDeltaPz%s", particleNames[iHad].data()), Form("Pz difference reco - MC %s; #it{p}_{z}^{reco} - #it{p}_{z}^{gen} (GeV/#it{c}); entries", labels[iHad].data()), HistType::kTH1F, {{2000, -1., 1.}});
      histDeltaSecondaryVertexX[iHad] = registry.add<TH1>(Form("histDeltaSecondaryVertexX%s", particleNames[iHad].data()), Form("Sec. Vertex difference reco - MC (MC matched) - %s; #Delta x (cm); entries", labels[iHad].data()), HistType::kTH1F, {{200, -1., 1.}});
      histDeltaSecondaryVertexY[iHad] = registry.add<TH1>(Form("histDeltaSecondaryVertexY%s", particleNames[iHad].data()), Form("Sec. Vertex difference reco - MC (MC matched) - %s; #Delta y (cm); entries", labels[iHad].data()), HistType::kTH1F, {{200, -1., 1.}});
      histDeltaSecondaryVertexZ[iHad] = registry.add<TH1>(Form("histDeltaSecondaryVertexZ%s", particleNames[iHad].data()), Form("Sec. Vertex difference reco - MC (MC matched) - %s; #Delta z (cm); entries", labels[iHad].data()), HistType::kTH1F, {{200, -1., 1.}});
      histDeltaDecayLength[iHad] = registry.add<TH1>(Form("histDeltaDecayLength%s", particleNames[iHad].data()), Form("Decay length difference reco - MC (%s); #Delta L (cm); entries", labels[iHad].data()), HistType::kTH1F, {{200, -1., 1.}});
      for (auto iOrigin = 0; iOrigin < 2; ++iOrigin) {
        histPtReco[iHad][iOrigin] = registry.add<TH1>(Form("histPtReco%s%s", originNames[iOrigin].data(), particleNames[iHad].data()), Form("Pt reco %s %s; #it{p}_{T}^{reco} (GeV/#it{c}); entries", originNames[iOrigin].data(), labels[iHad].data()), HistType::kTH1F, {{100, 0., 50.}});
        for (auto iDau = 0; iDau < nDaughters[iHad]; ++iDau) {
          histPtDau[iHad][iOrigin][iDau] = registry.add<TH1>(Form("histPtDau%d%s%s", iDau, originNames[iOrigin].data(), particleNames[iHad].data()), Form("Daughter %d Pt reco - %s %s; #it{p}_{T}^{dau, reco} (GeV/#it{c}); entries", iDau, originNames[iOrigin].data(), labels[iHad].data()), HistType::kTH1F, {{50, 0., 25.}});
          histEtaDau[iHad][iOrigin][iDau] = registry.add<TH1>(Form("histEtaDau%d%s%s", iDau, originNames[iOrigin].data(), particleNames[iHad].data()), Form("Daughter %d Eta reco - %s %s; #it{#eta}^{dau, reco}; entries", iDau, originNames[iOrigin].data(), labels[iHad].data()), HistType::kTH1F, {{100, -1., 1.}});
          histImpactParameterDau[iHad][iOrigin][iDau] = registry.add<TH1>(Form("histImpactParameterDau%d%s%s", iDau, originNames[iOrigin].data(), particleNames[iHad].data()), Form("Daughter %d DCAxy reco - %s %s; DCAxy (cm); entries", iDau, originNames[iOrigin].data(), labels[iHad].data()), HistType::kTH1F, {{200, -1., 1.}});
        }
      }
    }
  }

  using HfCandProng2WithMCRec = soa::Join<aod::HfCandProng2, aod::HfCandProng2MCRec>;
  using HfCandProng3WithMCRec = soa::Join<aod::HfCandProng3, aod::HfCandProng3MCRec>;

  void process(HfCandProng2WithMCRec const& cand2Prongs, HfCandProng3WithMCRec const& cand3Prongs, aod::BigTracksMC const& tracks, aod::McParticles const& particlesMC)
  {
    // loop over tracks
    for (auto& track : tracks) {
      uint index = uint(track.collisionId() >= 0);
      if (track.has_mcParticle()) {
        auto particle = track.mcParticle(); // get corresponding MC particle to check origin
        auto origin = RecoDecay::getCharmHadronOrigin(particlesMC, particle, true);
        histOriginTracks[index]->Fill(origin);
      } else {
        histOriginTracks[index]->Fill(-1.);
      }
    }

    // loop over 2-prong candidates
    for (auto& cand2Prong : cand2Prongs) {

      // determine which kind of candidate it is
      bool isD0Sel = TESTBIT(cand2Prong.hfflag(), o2::aod::hf_cand_prong2::DecayType::D0ToPiK);
      bool isJPsiSel = TESTBIT(cand2Prong.hfflag(), o2::aod::hf_cand_prong2::DecayType::JpsiToEE);
      if (!isD0Sel && !isJPsiSel) {
        continue;
      }
      int whichHad = -1;
      if (isD0Sel && TESTBIT(std::abs(cand2Prong.flagMCMatchRec()), hf_cand_prong2::DecayType::D0ToPiK)) {
        whichHad = 2;
      } else if (isJPsiSel && TESTBIT(std::abs(cand2Prong.flagMCMatchRec()), hf_cand_prong2::DecayType::JpsiToEE)) {
        whichHad = 6;
      }
      int whichOrigin = -1;
      if (cand2Prong.originMCRec() == RecoDecay::OriginType::Prompt) {
        whichOrigin = 0;
      } else {
        whichOrigin = 1;
      }

      if (whichHad >= 0 && whichOrigin >= 0) {
        int indexParticle = 0;
        if (cand2Prong.index0_as<aod::BigTracksMC>().has_mcParticle()) {
          indexParticle = RecoDecay::getMother(particlesMC, cand2Prong.index0_as<aod::BigTracksMC>().mcParticle(), PDGArrayParticle[whichHad], true);
        }
        auto mother = particlesMC.rawIteratorAt(indexParticle);
        histDeltaPt[whichHad]->Fill(cand2Prong.pt() - mother.pt());
        histDeltaPx[whichHad]->Fill(cand2Prong.px() - mother.px());
        histDeltaPy[whichHad]->Fill(cand2Prong.py() - mother.py());
        histDeltaPz[whichHad]->Fill(cand2Prong.pz() - mother.pz());
        // Compare Secondary vertex and decay length with MC
        auto daughter0 = mother.daughters_as<aod::McParticles>().begin();
        double vertexDau[3] = {daughter0.vx(), daughter0.vy(), daughter0.vz()};
        double vertexMoth[3] = {mother.vx(), mother.vy(), mother.vz()};
        auto decayLength = RecoDecay::distance(vertexMoth, vertexDau);

        histDeltaSecondaryVertexX[whichHad]->Fill(cand2Prong.xSecondaryVertex() - vertexDau[0]);
        histDeltaSecondaryVertexY[whichHad]->Fill(cand2Prong.ySecondaryVertex() - vertexDau[1]);
        histDeltaSecondaryVertexZ[whichHad]->Fill(cand2Prong.zSecondaryVertex() - vertexDau[2]);
        histDeltaDecayLength[whichHad]->Fill(cand2Prong.decayLength() - decayLength);
        std::array<double, 3> momDau0 = {cand2Prong.pxProng0(),
                                         cand2Prong.pyProng0(),
                                         cand2Prong.pzProng0()};
        std::array<double, 3> momDau1 = {cand2Prong.pxProng1(),
                                         cand2Prong.pyProng1(),
                                         cand2Prong.pzProng1()};
        histPtReco[whichHad][whichOrigin]->Fill(cand2Prong.pt());
        histPtDau[whichHad][whichOrigin][0]->Fill(RecoDecay::pt(momDau0));
        histEtaDau[whichHad][whichOrigin][0]->Fill(RecoDecay::eta(momDau0));
        histImpactParameterDau[whichHad][whichOrigin][0]->Fill(cand2Prong.impactParameter0());
        histPtDau[whichHad][whichOrigin][1]->Fill(RecoDecay::pt(momDau1));
        histEtaDau[whichHad][whichOrigin][1]->Fill(RecoDecay::eta(momDau1));
        histImpactParameterDau[whichHad][whichOrigin][1]->Fill(cand2Prong.impactParameter1());
      }
    } // end loop on 2-prong candidates

    // loop over 3-prong candidates
    for (auto& cand3Prong : cand3Prongs) {

      // determine which kind of candidate it is
      bool isDPlusSel = TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_prong3::DecayType::DPlusToPiKPi);
      bool isDStarSel = false; // FIXME: add proper check when D* will be added in HF vertexing
      bool isDsSel = TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_prong3::DecayType::DsToPiKK);
      bool isLcSel = TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_prong3::DecayType::LcToPKPi);
      bool isXicSel = TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_prong3::DecayType::XicToPKPi);
      if (!isDPlusSel && !isDStarSel && !isDsSel && !isLcSel && !isXicSel) {
        continue;
      }
      int whichHad = -1;
      if (isDPlusSel && TESTBIT(std::abs(cand3Prong.flagMCMatchRec()), hf_cand_prong3::DecayType::DPlusToPiKPi)) {
        whichHad = 0;
      } else if (isDsSel && TESTBIT(std::abs(cand3Prong.flagMCMatchRec()), hf_cand_prong3::DecayType::DsToPiKK)) {
        whichHad = 3;
      } else if (isLcSel && TESTBIT(std::abs(cand3Prong.flagMCMatchRec()), hf_cand_prong3::DecayType::LcToPKPi)) {
        whichHad = 4;
      } else if (isXicSel && TESTBIT(std::abs(cand3Prong.flagMCMatchRec()), hf_cand_prong3::DecayType::XicToPKPi)) {
        whichHad = 5;
      }
      int whichOrigin = -1;
      if (cand3Prong.originMCRec() == RecoDecay::OriginType::Prompt) {
        whichOrigin = 0;
      } else {
        whichOrigin = 1;
      }

      if (whichHad >= 0) {
        int indexParticle = 0;
        if (cand3Prong.index0_as<aod::BigTracksMC>().has_mcParticle()) {
          indexParticle = RecoDecay::getMother(particlesMC, cand3Prong.index0_as<aod::BigTracksMC>().mcParticle(), PDGArrayParticle[whichHad], true);
        }
        auto mother = particlesMC.rawIteratorAt(indexParticle);
        histDeltaPt[whichHad]->Fill(cand3Prong.pt() - mother.pt());
        histDeltaPx[whichHad]->Fill(cand3Prong.px() - mother.px());
        histDeltaPy[whichHad]->Fill(cand3Prong.py() - mother.py());
        histDeltaPz[whichHad]->Fill(cand3Prong.pz() - mother.pz());
        // Compare Secondary vertex and decay length with MC
        auto daughter0 = mother.daughters_as<aod::McParticles>().begin();
        double vertexDau[3] = {daughter0.vx(), daughter0.vy(), daughter0.vz()};
        double vertexMoth[3] = {mother.vx(), mother.vy(), mother.vz()};
        auto decayLength = RecoDecay::distance(vertexMoth, vertexDau);

        histDeltaSecondaryVertexX[whichHad]->Fill(cand3Prong.xSecondaryVertex() - vertexDau[0]);
        histDeltaSecondaryVertexY[whichHad]->Fill(cand3Prong.ySecondaryVertex() - vertexDau[1]);
        histDeltaSecondaryVertexZ[whichHad]->Fill(cand3Prong.zSecondaryVertex() - vertexDau[2]);
        histDeltaDecayLength[whichHad]->Fill(cand3Prong.decayLength() - decayLength);
        std::array<double, 3> momDau0 = {cand3Prong.pxProng0(),
                                         cand3Prong.pyProng0(),
                                         cand3Prong.pzProng0()};
        std::array<double, 3> momDau1 = {cand3Prong.pxProng1(),
                                         cand3Prong.pyProng1(),
                                         cand3Prong.pzProng1()};
        std::array<double, 3> momDau2 = {cand3Prong.pxProng2(),
                                         cand3Prong.pyProng2(),
                                         cand3Prong.pzProng2()};
        histPtReco[whichHad][whichOrigin]->Fill(cand3Prong.pt());
        histPtDau[whichHad][whichOrigin][0]->Fill(RecoDecay::pt(momDau0));
        histEtaDau[whichHad][whichOrigin][0]->Fill(RecoDecay::eta(momDau0));
        histImpactParameterDau[whichHad][whichOrigin][0]->Fill(cand3Prong.impactParameter0());
        histPtDau[whichHad][whichOrigin][1]->Fill(RecoDecay::pt(momDau1));
        histEtaDau[whichHad][whichOrigin][1]->Fill(RecoDecay::eta(momDau1));
        histImpactParameterDau[whichHad][whichOrigin][1]->Fill(cand3Prong.impactParameter1());
        histPtDau[whichHad][whichOrigin][2]->Fill(RecoDecay::pt(momDau2));
        histEtaDau[whichHad][whichOrigin][2]->Fill(RecoDecay::eta(momDau2));
        histImpactParameterDau[whichHad][whichOrigin][2]->Fill(cand3Prong.impactParameter2());
      }
    } // end loop on 3-prong candidates
  }   // end process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<ValidationGenLevel>(cfgc, TaskName{"hf-mc-validation-gen"}),
    adaptAnalysisTask<ValidationRecLevel>(cfgc, TaskName{"hf-mc-validation-rec"})};
  return workflow;
}
