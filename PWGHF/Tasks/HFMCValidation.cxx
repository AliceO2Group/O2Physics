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
/// \note Gen. and rec. level validation
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

namespace
{
static const int nCharmHad = 7;
static const std::array<int, nCharmHad> PDGArrayParticle = {pdg::Code::kDPlus, 413, pdg::Code::kD0, 431, pdg::Code::kLambdaCPlus, pdg::Code::kXiCPlus, pdg::Code::kJpsi};
static const std::array<std::string, nCharmHad> labels = {"D^{#plus}", "D*^{#plus}", "D^{0}", "D_{s}^{#plus}", "#Lambda_{c}^{#plus} #rightarrow pK^{#minus}#pi^{#plus}", "#Xi_{c}^{#plus} #rightarrow pK^{#minus}#pi^{#plus}", "J/#psi #rightarrow e^{#plus}e^{#minus}"};
static const std::array<std::string, nCharmHad> particleNames = {"Dplus", "Dstar", "D0", "Ds", "Lc2pKpi", "Xic2pKpi", "Jpsi2ee"};
} // namespace

/// Gen Level Validation
///
/// - Number of HF quarks produced per collision
/// - Number of D±      → π± K∓ π±        per collision
///             D*±     → π± K∓ π±,
///             Ds±     → K± K∓ π±,
///             D0(bar) → π± K∓,
///             Λc±     → p(bar) K∓ π±
///             Ξc±     → p(bar) K∓ π±
///             J/psi   → e∓ e±
/// - Momentum Conservation for these particles

struct ValidationGenLevel {
  std::shared_ptr<TH1> hCharmHaronsPtDistr, hCharmHaronsYDistr;

  HistogramRegistry registry{
    "registry",
    {{"hMomentumCheck", "Mom. Conservation (1 = true, 0 = false) (#it{#epsilon} = 1 MeV/#it{c}); Mom. Conservation result; entries", {HistType::kTH1F, {{2, -0.5, +1.5}}}},
     {"hPtDiffMotherDaughterGen", "Pt Difference Mother-Daughters; #Delta#it{p}_{T}^{gen.} (GeV/#it{c}); entries", {HistType::kTH1F, {{100, -0.01, 0.01}}}},
     {"hPxDiffMotherDaughterGen", "Px Difference Mother-Daughters; #Delta#it{p}_{x}^{gen.} (GeV/#it{c}); entries", {HistType::kTH1F, {{100, -0.01, 0.01}}}},
     {"hPyDiffMotherDaughterGen", "Py Difference Mother-Daughters; #Delta#it{p}_{y}^{gen.} (GeV/#it{c}); entries", {HistType::kTH1F, {{100, -0.01, 0.01}}}},
     {"hPzDiffMotherDaughterGen", "Pz Difference Mother-Daughters; #Delta#it{p}_{z}^{gen.} (GeV/#it{c}); entries", {HistType::kTH1F, {{100, -0.01, 0.01}}}},
     {"hPdiffMotherDaughterGen", "P  Difference Mother-Daughters; #Delta#it{p}^{gen.} (GeV/#it{c}); entries", {HistType::kTH1F, {{100, -0.01, 0.01}}}},
     {"hCountAverageC", "Event counter - Average Number Charm quark; Events Per Collision; entries", {HistType::kTH1F, {{20, 0., 20.}}}},
     {"hCountAverageB", "Event counter - Average Number Beauty quark; Events Per Collision; entries", {HistType::kTH1F, {{20, 0., 20.}}}},
     {"hCountAverageCbar", "Event counter - Average Number Anti-Charm quark; Events Per Collision; entries", {HistType::kTH1F, {{20, 0., 20.}}}},
     {"hCountAverageBbar", "Event counter - Average Number Anti-Beauty quark; Events Per Collision; entries", {HistType::kTH1F, {{20, 0., 20.}}}},
     {"hCounterPerCollisionDzero", "Event counter - D0; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hCounterPerCollisionDplus", "Event counter - DPlus; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hCounterPerCollisionDs", "Event counter - Ds; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hCounterPerCollisionDstar", "Event counter - Dstar; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hCounterPerCollisionLambdaC", "Event counter - LambdaC; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hCounterPerCollisionXiC", "Event counter - XiC; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hCounterPerCollisionJPsi", "Event counter - JPsi; Events Per Collision; entries", {HistType::kTH1F, {{10, -0.5, +9.5}}}},
     {"hPtVsYCharmQuark", "Y vs. Pt - charm quarks ; #it{p}_{T}^{gen} (GeV/#it{c}); #it{y}^{gen}", {HistType::kTH2F, {{100, 0., 50.}, {100, -5., 5.}}}},
     {"hPtVsYBeautyQuark", "Y vs. Pt - beauty quarks ; #it{p}_{T}^{gen} (GeV/#it{c}); #it{y}^{gen}", {HistType::kTH2F, {{100, 0., 50.}, {100, -5., 5.}}}}}};

  void init(o2::framework::InitContext&)
  {
    hCharmHaronsPtDistr = registry.add<TH2>("hCharmHaronsPtDistr", "Pt distribution vs charm hadron; ; #it{p}_{T}^{gen} (GeV/#it{c})", HistType::kTH2F, {{7, -0.5, 6.5}, {100, 0., 50.}});
    hCharmHaronsYDistr = registry.add<TH2>("hCharmHaronsYDistr", "Y distribution vs charm hadron; ; #it{y}^{gen}", HistType::kTH2F, {{7, -0.5, 6.5}, {100, 0., 50.}});
    for (auto iBin = 1; iBin <= nCharmHad; ++iBin) {
      hCharmHaronsPtDistr->GetXaxis()->SetBinLabel(iBin, labels[iBin - 1].data());
      hCharmHaronsYDistr->GetXaxis()->SetBinLabel(iBin, labels[iBin - 1].data());
    }
  }

  void process(aod::McCollision const& mccollision, aod::McParticles const& particlesMC)
  {
    int cPerCollision = 0;
    int cBarPerCollision = 0;
    int bPerCollision = 0;
    int bBarPerCollision = 0;

    //Particles and their decay checked in the second part of the task
    std::array<int, nCharmHad> PDGArrayParticle = {pdg::Code::kDPlus, 413, pdg::Code::kD0, 431, pdg::Code::kLambdaCPlus, pdg::Code::kXiCPlus, pdg::Code::kJpsi};
    std::array<std::array<int, 3>, nCharmHad> arrPDGFinal = {{{kPiPlus, kPiPlus, -kKPlus}, {kPiPlus, kPiPlus, -kKPlus}, {-kKPlus, kPiPlus, 0}, {kPiPlus, kKPlus, -kKPlus}, {kProton, -kKPlus, kPiPlus}, {kProton, -kKPlus, kPiPlus}, {kElectron, -kElectron, 0}}};
    std::array<int, nCharmHad> counter{0};
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
      bool momentumCheck = 1;
      listDaughters.clear();

      // Checking the decay of the particles and the momentum conservation
      for (std::size_t iD = 0; iD < PDGArrayParticle.size(); ++iD) {
        int whichHadron = -1;
        if (std::abs(particlePdgCode) == PDGArrayParticle[iD]) {
          whichHadron = iD;
          RecoDecay::getDaughters(particlesMC, particle, &listDaughters, arrPDGFinal[iD], -1);
          std::size_t arrayPDGsize = arrPDGFinal[iD].size() - std::count(arrPDGFinal[iD].begin(), arrPDGFinal[iD].end(), 0);
          if (listDaughters.size() == arrayPDGsize) {
            counter[iD]++;
          }
          for (std::size_t iDau = 0; iDau < listDaughters.size(); ++iDau) {
            auto daughter = particlesMC.rawIteratorAt(listDaughters.at(iDau));
            sumPxDau += daughter.px();
            sumPyDau += daughter.py();
            sumPzDau += daughter.pz();
          }
          auto pxDiff = particle.px() - sumPxDau;
          auto pyDiff = particle.py() - sumPyDau;
          auto pzDiff = particle.pz() - sumPzDau;
          if (std::abs(pxDiff) > 0.001 || std::abs(pyDiff) > 0.001 || std::abs(pzDiff) > 0.001) {
            momentumCheck = 0;
          }
          double pDiff = RecoDecay::P(pxDiff, pyDiff, pzDiff);
          double ptDiff = RecoDecay::Pt(pxDiff, pyDiff);
          //Filling histograms with per-component momentum conservation
          registry.fill(HIST("hMomentumCheck"), momentumCheck);
          registry.fill(HIST("hPxDiffMotherDaughterGen"), pxDiff);
          registry.fill(HIST("hPyDiffMotherDaughterGen"), pyDiff);
          registry.fill(HIST("hPzDiffMotherDaughterGen"), pzDiff);
          registry.fill(HIST("hPdiffMotherDaughterGen"), pDiff);
          registry.fill(HIST("hPtDiffMotherDaughterGen"), ptDiff);
          hCharmHaronsPtDistr->Fill(whichHadron, particle.pt());
          hCharmHaronsYDistr->Fill(whichHadron, particle.y());
        }
      }
    } //end particles
    registry.fill(HIST("hCountAverageC"), cPerCollision);
    registry.fill(HIST("hCountAverageB"), bPerCollision);
    registry.fill(HIST("hCountAverageCbar"), cBarPerCollision);
    registry.fill(HIST("hCountAverageBbar"), bBarPerCollision);
    registry.fill(HIST("hCounterPerCollisionDplus"), counter[0]);
    registry.fill(HIST("hCounterPerCollisionDstar"), counter[1]);
    registry.fill(HIST("hCounterPerCollisionDzero"), counter[2]);
    registry.fill(HIST("hCounterPerCollisionDs"), counter[3]);
    registry.fill(HIST("hCounterPerCollisionLambdaC"), counter[4]);
    registry.fill(HIST("hCounterPerCollisionXiC"), counter[5]);
    registry.fill(HIST("hCounterPerCollisionJPsi"), counter[6]);
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

  static const int nCharmHad = 7;
  std::array<std::shared_ptr<TH1>, nCharmHad> histPt, histPx, histPy, histPz, histSecVx, histSecVy, histSecVz, histDecLen;

  HistogramRegistry registry{"registry", {}};
  void init(o2::framework::InitContext&)
  {
    for (auto iHad = 0; iHad < nCharmHad; ++iHad) {
      histPt[iHad] = registry.add<TH1>(Form("histPt%s", particleNames[iHad].data()), Form("Pt difference reco - MC %s; #it{p}_{T}^{reco} - #it{p}_{T}^{gen.} (GeV/#it{c}); entries", labels[iHad].data()), HistType::kTH1F, {{2000, -1, 1}});
      histPx[iHad] = registry.add<TH1>(Form("histPx%s", particleNames[iHad].data()), Form("Px difference reco - MC %s; #it{p}_{x}^{reco} - #it{p}_{x}^{gen.} (GeV/#it{c}); entries", labels[iHad].data()), HistType::kTH1F, {{2000, -1, 1}});
      histPy[iHad] = registry.add<TH1>(Form("histPy%s", particleNames[iHad].data()), Form("Py difference reco - MC %s; #it{p}_{y}^{reco} - #it{p}_{y}^{gen.} (GeV/#it{c}); entries", labels[iHad].data()), HistType::kTH1F, {{2000, -1, 1}});
      histPz[iHad] = registry.add<TH1>(Form("histPz%s", particleNames[iHad].data()), Form("Pz difference reco - MC %s; #it{p}_{z}^{reco} - #it{p}_{z}^{gen.} (GeV/#it{c}); entries", labels[iHad].data()), HistType::kTH1F, {{2000, -1, 1}});
      histSecVx[iHad] = registry.add<TH1>(Form("histPt%s", particleNames[iHad].data()), Form("Sec. Vertex difference reco - MC (MC matched) - %s; #Delta x (cm); entries", labels[iHad].data()), HistType::kTH1F, {{200, -1, 1}});
      histSecVy[iHad] = registry.add<TH1>(Form("histPx%s", particleNames[iHad].data()), Form("Sec. Vertex difference reco - MC (MC matched) - %s; #Delta y (cm); entries", labels[iHad].data()), HistType::kTH1F, {{200, -1, 1}});
      histSecVz[iHad] = registry.add<TH1>(Form("histPy%s", particleNames[iHad].data()), Form("Sec. Vertex difference reco - MC (MC matched) - %s; #Delta z (cm); entries", labels[iHad].data()), HistType::kTH1F, {{200, -1, 1}});
      histDecLen[iHad] = registry.add<TH1>(Form("histPz%s", particleNames[iHad].data()), Form("Pz difference reco - MC (%s); #Delta L (cm); entries", labels[iHad].data()), HistType::kTH1F, {{200, -1, 1}});
    }
  }

  using HfCandProng2WithMCRec = soa::Join<aod::HfCandProng2, aod::HfCandProng2MCRec>;
  using HfCandProng3WithMCRec = soa::Join<aod::HfCandProng3, aod::HfCandProng3MCRec>;

  void process(HfCandProng2WithMCRec const& cand2Prongs, HfCandProng3WithMCRec const& cand3Prongs, aod::BigTracksMC const& tracks, aod::McParticles const& particlesMC)
  {
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

      if (whichHad >= 0) {
        int indexParticle = 0;
        if (cand2Prong.index0_as<aod::BigTracksMC>().has_mcParticle()) {
          indexParticle = RecoDecay::getMother(particlesMC, cand2Prong.index0_as<aod::BigTracksMC>().mcParticle(), PDGArrayParticle[whichHad], true);
        }
        auto mother = particlesMC.rawIteratorAt(indexParticle);
        histPt[whichHad]->Fill(cand2Prong.pt() - mother.pt());
        histPx[whichHad]->Fill(cand2Prong.px() - mother.px());
        histPy[whichHad]->Fill(cand2Prong.py() - mother.py());
        histPz[whichHad]->Fill(cand2Prong.pz() - mother.pz());
        //Compare Secondary vertex and decay length with MC
        auto daughter0 = mother.daughters_as<aod::McParticles>().begin();
        double vertexDau[3] = {daughter0.vx(), daughter0.vy(), daughter0.vz()};
        double vertexMoth[3] = {mother.vx(), mother.vy(), mother.vz()};
        auto decayLength = RecoDecay::distance(vertexMoth, vertexDau);

        histSecVx[whichHad]->Fill(cand2Prong.xSecondaryVertex() - vertexDau[0]);
        histSecVy[whichHad]->Fill(cand2Prong.ySecondaryVertex() - vertexDau[1]);
        histSecVz[whichHad]->Fill(cand2Prong.zSecondaryVertex() - vertexDau[2]);
        histDecLen[whichHad]->Fill(cand2Prong.decayLength() - decayLength);
      }
    } //end loop on 2-prong candidates

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

      if (whichHad >= 0) {
        int indexParticle = 0;
        if (cand3Prong.index0_as<aod::BigTracksMC>().has_mcParticle()) {
          indexParticle = RecoDecay::getMother(particlesMC, cand3Prong.index0_as<aod::BigTracksMC>().mcParticle(), PDGArrayParticle[whichHad], true);
        }
        auto mother = particlesMC.rawIteratorAt(indexParticle);
        histPt[whichHad]->Fill(cand3Prong.pt() - mother.pt());
        histPx[whichHad]->Fill(cand3Prong.px() - mother.px());
        histPy[whichHad]->Fill(cand3Prong.py() - mother.py());
        histPz[whichHad]->Fill(cand3Prong.pz() - mother.pz());
        //Compare Secondary vertex and decay length with MC
        auto daughter0 = mother.daughters_as<aod::McParticles>().begin();
        double vertexDau[3] = {daughter0.vx(), daughter0.vy(), daughter0.vz()};
        double vertexMoth[3] = {mother.vx(), mother.vy(), mother.vz()};
        auto decayLength = RecoDecay::distance(vertexMoth, vertexDau);

        histSecVx[whichHad]->Fill(cand3Prong.xSecondaryVertex() - vertexDau[0]);
        histSecVy[whichHad]->Fill(cand3Prong.ySecondaryVertex() - vertexDau[1]);
        histSecVz[whichHad]->Fill(cand3Prong.zSecondaryVertex() - vertexDau[2]);
        histDecLen[whichHad]->Fill(cand3Prong.decayLength() - decayLength);
      }
    } //end loop on 3-prong candidates
  }   //end process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<ValidationGenLevel>(cfgc, TaskName{"hf-mc-validation-gen"}),
    adaptAnalysisTask<ValidationRecLevel>(cfgc, TaskName{"hf-mc-validation-rec"})};
  return workflow;
}
