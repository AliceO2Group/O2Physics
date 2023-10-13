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

/// \file taskJPsiHf.cxx
/// \brief Task for the analysis of J/psi - open HF associate production
/// \author Luca Micheletti <luca.micheletti@to.infn.it>, INFN
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include <string>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

// Declarations of various short names
using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyPairCandidatesSelected = soa::Join<aod::Dileptons, aod::DileptonsExtra, aod::DileptonsInfo>;
using MyD0CandidatesSelected = soa::Join<aod::HfCand2Prong, aod::HfSelD0>;

struct taskJPsiHf {
  //
  // This task combines dilepton candidates with a open charm hadron
  //

  // TODO: For now this is only used to determine the position in the filter bit map for the hadron cut
  Configurable<std::string> fConfigTrackCuts{"cfgLeptonCuts", "jpsiO2MCdebugCuts2", "Comma separated list of barrel track cuts"};
  // comment: add list of subgroups (must define subgroups under )
  Configurable<std::string> fConfigAddDileptonHadHistogram{"cfgAddDileptonHadHistogram", "", "Comma separated list of histograms"};

  // HF configurables
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<int> selectionFlagHf{"selectionFlagHf", 1, "Selection Flag for HF flagged candidates"};
  Configurable<int> selectionTopol{"selectionTopol", 1, "Selection Flag for topologically selected candidates"};
  Configurable<int> selectionCand{"selectionCand", 1, "Selection Flag for conj. topol. selected candidates"};
  Configurable<int> selectionPid{"selectionPid", 1, "Selection Flag for reco PID candidates"};
  // DQ configurables
  Configurable<double> massDileptonCandMin{"massDileptonCandMin", 1, "minimum dilepton mass"};
  Configurable<double> massDileptonCandMax{"massDileptonCandMax", 5, "maximum dilepton mass"};
  // General configurables
  Configurable<bool> configDebug{"configDebug", true, "If true, fill D0 - J/psi histograms separately"};

  SliceCache cache;
  HfHelper hfHelper;

  Partition<MyPairCandidatesSelected> selectedDileptonCandidates = aod::reducedpair::mass > 1.0f && aod::reducedpair::mass < 5.0f && aod::reducedpair::sign == 0;
  Partition<MyD0CandidatesSelected> selectedD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;

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
                             {{"JPsiDmeson/hSparseJPsiDmeson", ";#it{p}_{T}(J/#psi) (GeV/#it{c});#it{p}_{T}(D) (GeV/#it{c};#it{M}(J/#psi) (GeV/#it{c}^{2});#it{M}(D) (GeV/#it{c}^{2});#Delta#it{y};#Delta#varphi", {HistType::kTHnSparseF, {axisPt, axisPt, axisMassJPsi, axisMassDmeson, axisDeltaY, axisPhi}}},
                              {"JPsiDmeson/hMassJPsiWithDmeson", ";#it{M}(J/#psi) (GeV/#it{c}^{2});counts", {HistType::kTH1F, {axisMassJPsi}}},
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
      registry.add("JPsi/hMassJPsi", ";#it{M}(J/#psi) (GeV/#it{c}^{2});counts", {HistType::kTH1F, {axisMassJPsi}});
      registry.add("JPsi/hPtJPsi", ";#it{p}_{T}(J/#psi) (GeV/#it{c});counts", {HistType::kTH1F, {axisPt}});
      registry.add("JPsi/hRapJPsi", ";#it{y}(J/#psi);counts", {HistType::kTH1F, {axisFwdY}});
      registry.add("JPsi/hPhiJPsi", ";#it{#varphi}(J/#psi);counts", {HistType::kTH1F, {axisPhi}});
      registry.add("Dmeson/hMassDmeson", ";#it{M}(D) (GeV/#it{c}^{2});counts", {HistType::kTH1F, {axisMassDmeson}});
      registry.add("Dmeson/hPtDmeson", ";#it{p}_{T}(D) (GeV/#it{c});counts", {HistType::kTH1F, {axisPt}});
      registry.add("Dmeson/hRapDmeson", ";#it{y}(D);counts", {HistType::kTH1F, {axisMidY}});
      registry.add("Dmeson/hPhiDmeson", ";#it{#varphi}(D);counts", {HistType::kTH1F, {axisPhi}});
    }
  }

  // Template function to run pair - hadron combinations
  // TODO: generalise to all charm-hadron species
  template <typename TDqTrack, typename THfTrack>
  void runDileptonDmeson(TDqTrack const& dileptons, THfTrack const& dmesons)
  {
    if (configDebug) {
      for (auto& dmeson : dmesons) {
        if (!TESTBIT(dmeson.hfflag(), hf_cand_2prong::DecayType::D0ToPiK)) {
          continue;
        }

        auto rapD0 = hfHelper.yD0(dmeson);

        if (yCandMax >= 0. && std::abs(hfHelper.yD0(dmeson)) > yCandMax) {
          continue;
        }

        auto ptD0 = dmeson.pt();
        auto phiD0 = dmeson.phi();
        auto massD0 = -1.;
        auto massD0bar = -1.;

        if (dmeson.isSelD0() >= selectionFlagD0) {
          massD0 = hfHelper.invMassD0ToPiK(dmeson);
          registry.fill(HIST("Dmeson/hMassDmeson"), massD0);
          registry.fill(HIST("Dmeson/hPtDmeson"), ptD0);
          registry.fill(HIST("Dmeson/hRapDmeson"), rapD0);
          registry.fill(HIST("Dmeson/hPhiDmeson"), phiD0);
        }

        if (dmeson.isSelD0bar() >= selectionFlagD0bar) {
          massD0bar = hfHelper.invMassD0ToPiK(dmeson);
          registry.fill(HIST("Dmeson/hMassDmeson"), massD0bar);
          registry.fill(HIST("Dmeson/hPtDmeson"), ptD0);
          registry.fill(HIST("Dmeson/hRapDmeson"), rapD0);
          registry.fill(HIST("Dmeson/hPhiDmeson"), phiD0);
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
        registry.fill(HIST("JPsi/hMassJPsi"), massJPsi);
        registry.fill(HIST("JPsi/hPtJPsi"), ptJPsi);
        registry.fill(HIST("JPsi/hRapJPsi"), rapJPsi);
        registry.fill(HIST("JPsi/hPhiJPsi"), phiJPsi);
      }

      // loop over D mesons
      for (auto& dmeson : dmesons) {
        if (!TESTBIT(dmeson.hfflag(), hf_cand_2prong::DecayType::D0ToPiK)) {
          continue;
        }

        auto rapD0 = hfHelper.yD0(dmeson);

        if (yCandMax >= 0. && std::abs(rapD0) > yCandMax) {
          continue;
        }

        auto ptD0 = dmeson.pt();
        auto phiD0 = dmeson.phi();
        auto massD0 = -1.;
        auto massD0bar = -1.;
        auto rapDelta = rapD0 - rapJPsi;
        auto phiDelta = std::abs(phiJPsi - phiD0);

        if (dmeson.isSelD0() >= selectionFlagD0) {
          massD0 = hfHelper.invMassD0ToPiK(dmeson);
          registry.fill(HIST("JPsiDmeson/hSparseJPsiDmeson"), ptJPsi, ptD0, massJPsi, massD0, rapDelta, phiDelta);
          registry.fill(HIST("JPsiDmeson/hMassJPsiWithDmeson"), massJPsi);
          registry.fill(HIST("JPsiDmeson/hPtJPsiWithDmeson"), ptJPsi);
          registry.fill(HIST("JPsiDmeson/hRapJPsiWithDmeson"), rapJPsi);
          registry.fill(HIST("JPsiDmeson/hPhiJPsiWithDmeson"), phiJPsi);
          registry.fill(HIST("JPsiDmeson/hMassDmesonWithJPsi"), massD0);
          registry.fill(HIST("JPsiDmeson/hPtDmesonWithJPsi"), ptD0);
          registry.fill(HIST("JPsiDmeson/hRapDmesonWithJPsi"), rapD0);
          registry.fill(HIST("JPsiDmeson/hPhiDmesonWithJPsi"), phiD0);
        }
        if (dmeson.isSelD0bar() >= selectionFlagD0bar) {
          massD0bar = hfHelper.invMassD0barToKPi(dmeson);
          registry.fill(HIST("JPsiDmeson/hSparseJPsiDmeson"), ptJPsi, ptD0, massJPsi, massD0bar, rapDelta, phiDelta);
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

  // process J/psi - D0
  void processJspiD0(MyEvents const& collisions, MyPairCandidatesSelected const& dileptons, MyD0CandidatesSelected const& dmesons)
  {
    for (auto& collision : collisions) {
      auto groupedDmesonCandidates = selectedD0Candidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      auto groupedDileptonCandidates = selectedDileptonCandidates->sliceByCached(aod::reducedpair::collisionId, collision.globalIndex(), cache);
      runDileptonDmeson(groupedDileptonCandidates, groupedDmesonCandidates);
    }
  }

  PROCESS_SWITCH(taskJPsiHf, processJspiD0, "Process J/psi - D0", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<taskJPsiHf>(cfgc)};
}
