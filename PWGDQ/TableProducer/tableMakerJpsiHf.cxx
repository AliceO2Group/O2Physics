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

/// \file tableMakerJpsiHf.cxx
/// \brief Task for the production of the derived data of J/psi - open HF associate production
/// \author Luca Micheletti <luca.micheletti@to.infn.it>, INFN
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include <string>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::aod::hf_cand_2prong;

// Declarations of various short names
using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyPairCandidatesSelected = soa::Join<aod::Dileptons, aod::DileptonsExtra, aod::DileptonsInfo>;
using MyD0CandidatesSelected = soa::Join<aod::HfCand2Prong, aod::HfSelD0>;
using MyD0CandidatesSelectedWithBdt = soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0>;

constexpr float cutsBdt[1][3] = {{1., 0., 0.}}; // background, prompt, nonprompt
static const std::vector<std::string> labelsBdt = {"Background", "Prompt", "Nonprompt"};
static const std::vector<std::string> labelsEmpty{};

struct tableMakerJpsiHf {
  //
  // This task combines dilepton candidates with a open charm hadron
  //

  // Produce derived tables
  Produces<RedJpDmColls> redCollisions;
  Produces<RedJpDmDmesons> redDmesons;
  Produces<RedJpDmDmesBdts> redDmesBdts;
  Produces<RedJpDmDileptons> redDileptons;

  // TODO: For now this is only used to determine the position in the filter bit map for the hadron cut
  Configurable<std::string> fConfigTrackCuts{"cfgLeptonCuts", "jpsiO2MCdebugCuts2", "Comma separated list of barrel track cuts"};
  // comment: add list of subgroups (must define subgroups under )
  Configurable<std::string> fConfigAddDileptonHadHistogram{"cfgAddDileptonHadHistogram", "", "Comma separated list of histograms"};

  // HF configurables
  // cuts on BDT output scores to be applied only for the histograms
  Configurable<LabeledArray<float>> bdtCutsForHistos{"bdtCutsForHistos", {cutsBdt[0], 1, 3, labelsEmpty, labelsBdt}, "Additional bdt cut values only for histograms"};
  Configurable<double> yCandDmesonMax{"yCandDmesonMax", -1., "max. cand. rapidity"};
  // DQ configurables
  Configurable<double> massDileptonCandMin{"massDileptonCandMin", 1, "minimum dilepton mass"};
  Configurable<double> massDileptonCandMax{"massDileptonCandMax", 5, "maximum dilepton mass"};
  // General configurables
  Configurable<bool> configDebug{"configDebug", true, "If true, fill D0 - J/psi histograms separately"};

  SliceCache cache;
  Partition<MyPairCandidatesSelected> selectedDileptonCandidates = aod::reducedpair::mass > 1.0f && aod::reducedpair::mass < 5.0f && aod::reducedpair::sign == 0;
  Partition<MyD0CandidatesSelected> selectedD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= 1 || aod::hf_sel_candidate_d0::isSelD0bar >= 1;
  Partition<MyD0CandidatesSelectedWithBdt> selectedD0CandidatesWithBdt = aod::hf_sel_candidate_d0::isSelD0 >= 1 || aod::hf_sel_candidate_d0::isSelD0bar >= 1;

  Preslice<MyD0CandidatesSelected> perCollisionDmeson = aod::hf_cand::collisionId;
  Preslice<MyPairCandidatesSelected> perCollisionDilepton = aod::reducedpair::collisionId;

  // Define histograms
  AxisSpec axisPt{100, 0.f, 50.f};
  AxisSpec axisMassDmeson{200, 1.7f, 2.1f}; // TODO: make it dependent on the D-meson species
  AxisSpec axisMassJPsi{300, 2.f, 5.f};
  AxisSpec axisMidY{60, -1.5f, 1.5f};
  AxisSpec axisFwdY{50, -4.5f, -2.0f};
  AxisSpec axisDeltaY{90, 1.f, 5.5f};
  AxisSpec axisPhi{180, 0., 2 * constants::math::PI}; // same for delta phi

  HistogramRegistry registry{"registry",
                             {{"JPsiDmeson/hMassJPsiWithDmeson", ";#it{M}(J/#psi) (GeV/#it{c}^{2});counts", {HistType::kTH1F, {axisMassJPsi}}},
                              {"JPsiDmeson/hPtJPsiWithDmeson", ";#it{p}_{T}(J/#psi) (GeV/#it{c});counts", {HistType::kTH1F, {axisPt}}},
                              {"JPsiDmeson/hRapJPsiWithDmeson", ";#it{y}(J/#psi);counts", {HistType::kTH1F, {axisFwdY}}},
                              {"JPsiDmeson/hPhiJPsiWithDmeson", ";#it{#varphi}(J/#psi);counts", {HistType::kTH1F, {axisPhi}}},
                              {"JPsiDmeson/hMassDmesonWithJPsi", ";#it{M}(D) (GeV/#it{c}^{2});counts", {HistType::kTH1F, {axisMassDmeson}}},
                              {"JPsiDmeson/hPtDmesonWithJPsi", ";#it{p}_{T}(D) (GeV/#it{c});counts", {HistType::kTH1F, {axisPt}}},
                              {"JPsiDmeson/hRapDmesonWithJPsi", ";#it{y}(D);counts", {HistType::kTH1F, {axisMidY}}},
                              {"JPsiDmeson/hPhiDmesonWithJPsi", ";#it{#varphi}(D);counts", {HistType::kTH1F, {axisPhi}}}}};

  void init(o2::framework::InitContext& context)
  {
    if (configDebug) {
      registry.add("JPsi/hMassVsPtJPsi", ";#it{p}_{T}(J/#psi) (GeV/#it{c});#it{M}(J/#psi) (GeV/#it{c}^{2});counts", {HistType::kTH2F, {axisPt, axisMassJPsi}});
      registry.add("JPsi/hRapVsPtJPsi", ";#it{p}_{T}(J/#psi) (GeV/#it{c});#it{y}(J/#psi);counts", {HistType::kTH2F, {axisPt, axisFwdY}});
      registry.add("JPsi/hPhiJPsi", ";#it{#varphi}(J/#psi);counts", {HistType::kTH1F, {axisPhi}});
      registry.add("Dmeson/hMassVsPtDmeson", ";#it{p}_{T}(D) (GeV/#it{c});#it{M}(D) (GeV/#it{c}^{2});counts", {HistType::kTH2F, {axisPt, axisMassDmeson}});
      registry.add("Dmeson/hRapVsPtDmeson", ";#it{p}_{T}(D) (GeV/#it{c});#it{y}(D);counts", {HistType::kTH2F, {axisPt, axisMidY}});
      registry.add("Dmeson/hPhiDmeson", ";#it{#varphi}(D);counts", {HistType::kTH1F, {axisPhi}});
    }
  }

  // Template function to run pair - hadron combinations
  // TODO: generalise to all charm-hadron species
  template <bool withBdt, typename TDqTrack, typename THfTrack>
  void runDileptonDmeson(TDqTrack const& dileptons, THfTrack const& dmesons, MyEvents::iterator const& collision)
  {
    std::vector<float> scores{-1., 2., 2.};
    bool isCollSel{false};
    if (configDebug) {
      for (auto& dmeson : dmesons) {
        if (!TESTBIT(dmeson.hfflag(), DecayType::D0ToPiK)) {
          continue;
        }

        if constexpr (withBdt) {
          scores[0] = dmeson.mlProbD0()[0];
          scores[1] = dmeson.mlProbD0()[1];
          scores[2] = dmeson.mlProbD0()[2];
        }

        auto rapD0 = yD0(dmeson);

        if (yCandDmesonMax >= 0. && std::abs(yD0(dmeson)) > yCandDmesonMax) {
          continue;
        }

        auto ptD0 = dmeson.pt();
        auto phiD0 = dmeson.phi();
        auto massD0 = -1.;
        auto massD0bar = -1.;

        if (scores[0] < bdtCutsForHistos->get(0u, 0u) && scores[1] > bdtCutsForHistos->get(0u, 1u) && scores[2] > bdtCutsForHistos->get(0u, 2u)) {
          if (dmeson.isSelD0() >= 1) {
            massD0 = invMassD0ToPiK(dmeson);
            registry.fill(HIST("Dmeson/hMassVsPtDmeson"), ptD0, massD0);
            registry.fill(HIST("Dmeson/hRapVsPtDmeson"), ptD0, rapD0);
            registry.fill(HIST("Dmeson/hPhiDmeson"), phiD0);
          }

          if (dmeson.isSelD0bar() >= 1) {
            massD0bar = invMassD0ToPiK(dmeson);
            registry.fill(HIST("Dmeson/hMassVsPtDmeson"), ptD0, massD0bar);
            registry.fill(HIST("Dmeson/hRapVsPtDmeson"), ptD0, rapD0);
            registry.fill(HIST("Dmeson/hPhiDmeson"), phiD0);
          }
        }
      }
    }

    // loop over dileptons
    for (auto dilepton : dileptons) {
      auto massJPsi = dilepton.mass();
      auto ptJPsi = dilepton.pt();
      auto rapJPsi = dilepton.rap();
      auto phiJPsi = dilepton.phi() + constants::math::PI; // TODO: check conventions!

      if (massJPsi < massDileptonCandMin || massJPsi > massDileptonCandMax) {
        continue;
      }

      if (configDebug) {
        registry.fill(HIST("JPsi/hMassVsPtJPsi"), ptJPsi, massJPsi);
        registry.fill(HIST("JPsi/hRapVsPtJPsi"), ptJPsi, rapJPsi);
        registry.fill(HIST("JPsi/hPhiJPsi"), phiJPsi);
      }

      // loop over D mesons
      for (auto& dmeson : dmesons) {
        if (!TESTBIT(dmeson.hfflag(), DecayType::D0ToPiK)) {
          continue;
        }

        auto rapD0 = yD0(dmeson);

        if (yCandDmesonMax >= 0. && std::abs(rapD0) > yCandDmesonMax) {
          continue;
        }

        auto ptD0 = dmeson.pt();
        auto phiD0 = dmeson.phi();
        auto massD0 = -1.;
        auto massD0bar = -1.;

        if (dmeson.isSelD0() >= 1 || dmeson.isSelD0bar() >= 1) {
          if (!isCollSel) {
            redCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.numContrib());
            isCollSel = true;
          }
          auto indexRed = redCollisions.lastIndex();
          redDileptons(indexRed, dilepton.px(), dilepton.py(), dilepton.pz(), dilepton.sign(), dilepton.mcDecision(), dilepton.tauz(), dilepton.lz(), dilepton.lxy());
          redDmesons(indexRed, dmeson.px(), dmeson.py(), dmeson.pz(), dmeson.xSecondaryVertex(), dmeson.ySecondaryVertex(), dmeson.zSecondaryVertex(), 0, 0);
          if constexpr (withBdt) {
            auto scores = dmeson.mlProbD0();
            redDmesBdts(scores[0], scores[1], scores[2]);
          }

          if (dmeson.isSelD0() >= 1 && scores[0] < bdtCutsForHistos->get(0u, 0u) && scores[1] > bdtCutsForHistos->get(0u, 1u) && scores[2] > bdtCutsForHistos->get(0u, 2u)) {
            massD0 = invMassD0ToPiK(dmeson);
            registry.fill(HIST("JPsiDmeson/hMassJPsiWithDmeson"), massJPsi);
            registry.fill(HIST("JPsiDmeson/hPtJPsiWithDmeson"), ptJPsi);
            registry.fill(HIST("JPsiDmeson/hRapJPsiWithDmeson"), rapJPsi);
            registry.fill(HIST("JPsiDmeson/hPhiJPsiWithDmeson"), phiJPsi);
            registry.fill(HIST("JPsiDmeson/hMassDmesonWithJPsi"), massD0);
            registry.fill(HIST("JPsiDmeson/hPtDmesonWithJPsi"), ptD0);
            registry.fill(HIST("JPsiDmeson/hRapDmesonWithJPsi"), rapD0);
            registry.fill(HIST("JPsiDmeson/hPhiDmesonWithJPsi"), phiD0);
          }
          if (dmeson.isSelD0bar() >= 1 && scores[0] < bdtCutsForHistos->get(0u, 0u) && scores[1] > bdtCutsForHistos->get(0u, 1u) && scores[2] > bdtCutsForHistos->get(0u, 2u)) {
            massD0bar = invMassD0barToKPi(dmeson);
            registry.fill(HIST("JPsiDmeson/hMassJPsiWithDmeson"), massJPsi);
            registry.fill(HIST("JPsiDmeson/hPtJPsiWithDmeson"), ptJPsi);
            registry.fill(HIST("JPsiDmeson/hRapJPsiWithDmeson"), rapJPsi);
            registry.fill(HIST("JPsiDmeson/hPhiJPsiWithDmeson"), phiJPsi);
            registry.fill(HIST("JPsiDmeson/hMassDmesonWithJPsi"), massD0bar);
            registry.fill(HIST("JPsiDmeson/hPtDmesonWithJPsi"), ptD0);
            registry.fill(HIST("JPsiDmeson/hRapDmesonWithJPsi"), rapD0);
            registry.fill(HIST("JPsiDmeson/hPhiDmesonWithJPsi"), phiD0);
          }
        }
      }
    }
  }

  // process J/psi - D0
  void processJspiD0(MyEvents const& collisions, MyPairCandidatesSelected const& dileptons, MyD0CandidatesSelected const& dmesons)
  {
    for (auto& collision : collisions) {
      auto groupedDmesonCandidates = selectedD0Candidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      auto groupedDileptonCandidates = selectedDileptonCandidates->sliceByCached(aod::reducedpair::collisionId, collision.globalIndex(), cache);
      runDileptonDmeson<false>(groupedDileptonCandidates, groupedDmesonCandidates, collision);
    }
  }

  // process J/psi - D0 adding the BDT output scores to the D0 table
  void processJspiD0WithBdt(MyEvents const& collisions, MyPairCandidatesSelected const& dileptons, MyD0CandidatesSelectedWithBdt const& dmesons)
  {
    for (auto& collision : collisions) {
      auto groupedDmesonCandidates = selectedD0CandidatesWithBdt->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      auto groupedDileptonCandidates = selectedDileptonCandidates->sliceByCached(aod::reducedpair::collisionId, collision.globalIndex(), cache);
      runDileptonDmeson<true>(groupedDileptonCandidates, groupedDmesonCandidates, collision);
    }
  }

  PROCESS_SWITCH(tableMakerJpsiHf, processJspiD0, "Process J/psi - D0", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<tableMakerJpsiHf>(cfgc)};
}
