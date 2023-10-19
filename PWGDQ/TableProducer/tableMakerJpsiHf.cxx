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

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/HistogramsLibrary.h"

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

HfHelper hfHelper;

struct tableMakerJpsiHf {
  //
  // This task combines dilepton candidates with a open charm hadron
  //

  // Produce derived tables
  Produces<RedJpDmColls> redCollisions;
  Produces<RedJpDmDmesons> redDmesons;
  Produces<RedJpDmDmesBdts> redDmesBdts;
  Produces<RedJpDmD0Masss> redD0Masses;
  Produces<RedJpDmDileptons> redDileptons;
  Produces<RedJpDmColCounts> redCollCounter;

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

  // Define histograms manager
  float* fValuesDileptonCharmHadron{};
  HistogramManager* fHistMan{};
  OutputObj<THashList> fOutputList{"output"};

  void init(o2::framework::InitContext& context)
  {
    fValuesDileptonCharmHadron = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(true);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
    fHistMan->AddHistClass("Dmeson");
    fHistMan->AddHistClass("JPsi");
    fHistMan->AddHistClass("JPsiDmeson");
    dqhistograms::DefineHistograms(fHistMan, "Dmeson", "dilepton-charmhadron", "dmeson");
    dqhistograms::DefineHistograms(fHistMan, "JPsi", "dilepton-charmhadron", "jpsi");
    dqhistograms::DefineHistograms(fHistMan, "JPsiDmeson", "dilepton-charmhadron", "jpsidmeson");
    VarManager::SetUseVars(fHistMan->GetUsedVars());
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  // Template function to run pair - hadron combinations
  // TODO: generalise to all charm-hadron species
  template <bool withBdt, typename TDqTrack, typename THfTrack>
  void runDileptonDmeson(TDqTrack const& dileptons, THfTrack const& dmesons, MyEvents::iterator const& collision)
  {
    VarManager::ResetValues(0, VarManager::kNVars, fValuesDileptonCharmHadron);

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

        auto rapD0 = hfHelper.yD0(dmeson);

        if (yCandDmesonMax >= 0. && std::abs(rapD0) > yCandDmesonMax) {
          continue;
        }

        if (scores[0] < bdtCutsForHistos->get(0u, 0u) && scores[1] > bdtCutsForHistos->get(0u, 1u) && scores[2] > bdtCutsForHistos->get(0u, 2u)) {
          if (dmeson.isSelD0() >= 1) {
            VarManager::FillSingleDileptonCharmHadron<VarManager::kD0ToPiK>(dmeson, hfHelper, fValuesDileptonCharmHadron);
            fHistMan->FillHistClass("Dmeson", fValuesDileptonCharmHadron);
            VarManager::ResetValues(0, VarManager::kNVars, fValuesDileptonCharmHadron);
          }
          if (dmeson.isSelD0bar() >= 1) {
            VarManager::FillSingleDileptonCharmHadron<VarManager::kD0barToKPi>(dmeson, hfHelper, fValuesDileptonCharmHadron);
            fHistMan->FillHistClass("Dmeson", fValuesDileptonCharmHadron);
            VarManager::ResetValues(0, VarManager::kNVars, fValuesDileptonCharmHadron);
          }
        }
      }
    }

    // loop over dileptons
    for (auto dilepton : dileptons) {
      auto massJPsi = dilepton.mass();

      if (massJPsi < massDileptonCandMin || massJPsi > massDileptonCandMax) {
        continue;
      }

      if (configDebug) {
        VarManager::FillSingleDileptonCharmHadron<VarManager::kJPsiToMuMu>(dilepton, hfHelper, fValuesDileptonCharmHadron);
        fHistMan->FillHistClass("JPsi", fValuesDileptonCharmHadron);
        VarManager::ResetValues(0, VarManager::kNVars, fValuesDileptonCharmHadron);
      }

      // loop over D mesons
      for (auto& dmeson : dmesons) {
        if (!TESTBIT(dmeson.hfflag(), DecayType::D0ToPiK)) {
          continue;
        }

        auto rapD0 = hfHelper.yD0(dmeson);

        if (yCandDmesonMax >= 0. && std::abs(rapD0) > yCandDmesonMax) {
          continue;
        }

        auto massD0 = -1.;
        auto massD0bar = -1.;

        if (dmeson.isSelD0() >= 1 || dmeson.isSelD0bar() >= 1) {
          if (!isCollSel) {
            redCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.numContrib());
            isCollSel = true;
          }
          auto indexRed = redCollisions.lastIndex();
          redDileptons(indexRed, dilepton.px(), dilepton.py(), dilepton.pz(), dilepton.mass(), dilepton.sign(), dilepton.mcDecision(), dilepton.tauz(), dilepton.lz(), dilepton.lxy());
          redDmesons(indexRed, dmeson.px(), dmeson.py(), dmeson.pz(), dmeson.xSecondaryVertex(), dmeson.ySecondaryVertex(), dmeson.zSecondaryVertex(), 0, 0);
          if constexpr (withBdt) {
            auto scores = dmeson.mlProbD0();
            redDmesBdts(scores[0], scores[1], scores[2]);
          }

          if (dmeson.isSelD0() >= 1) {
            massD0 = hfHelper.invMassD0ToPiK(dmeson);
            if (scores[0] < bdtCutsForHistos->get(0u, 0u) && scores[1] > bdtCutsForHistos->get(0u, 1u) && scores[2] > bdtCutsForHistos->get(0u, 2u)) {
              VarManager::FillDileptonCharmHadron<VarManager::kD0ToPiK>(dilepton, dmeson, hfHelper, fValuesDileptonCharmHadron);
              fHistMan->FillHistClass("JPsiDmeson", fValuesDileptonCharmHadron);
              VarManager::ResetValues(0, VarManager::kNVars, fValuesDileptonCharmHadron);
            }
          }
          if (dmeson.isSelD0bar() >= 1) {
            massD0bar = hfHelper.invMassD0barToKPi(dmeson);
            if (scores[0] < bdtCutsForHistos->get(0u, 0u) && scores[1] > bdtCutsForHistos->get(0u, 1u) && scores[2] > bdtCutsForHistos->get(0u, 2u)) {
              VarManager::FillDileptonCharmHadron<VarManager::kD0barToKPi>(dilepton, dmeson, hfHelper, fValuesDileptonCharmHadron);
              fHistMan->FillHistClass("JPsiDmeson", fValuesDileptonCharmHadron);
              VarManager::ResetValues(0, VarManager::kNVars, fValuesDileptonCharmHadron);
            }
          }
          redD0Masses(massD0, massD0bar);
        }
      }
    }
  }

  // process J/psi - D0
  void processJspiD0(MyEvents const& collisions, MyPairCandidatesSelected const& dileptons, MyD0CandidatesSelected const& dmesons)
  {
    redCollCounter(collisions.size());
    for (auto& collision : collisions) {
      auto groupedDmesonCandidates = selectedD0Candidates->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      auto groupedDileptonCandidates = selectedDileptonCandidates->sliceByCached(aod::reducedpair::collisionId, collision.globalIndex(), cache);
      runDileptonDmeson<false>(groupedDileptonCandidates, groupedDmesonCandidates, collision);
    }
  }
  PROCESS_SWITCH(tableMakerJpsiHf, processJspiD0, "Process J/psi - D0", true);

  // process J/psi - D0 adding the BDT output scores to the D0 table
  void processJspiD0WithBdt(MyEvents const& collisions, MyPairCandidatesSelected const& dileptons, MyD0CandidatesSelectedWithBdt const& dmesons)
  {
    redCollCounter(collisions.size());
    for (auto& collision : collisions) {
      auto groupedDmesonCandidates = selectedD0CandidatesWithBdt->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
      auto groupedDileptonCandidates = selectedDileptonCandidates->sliceByCached(aod::reducedpair::collisionId, collision.globalIndex(), cache);
      runDileptonDmeson<true>(groupedDileptonCandidates, groupedDmesonCandidates, collision);
    }
  }
  PROCESS_SWITCH(tableMakerJpsiHf, processJspiD0WithBdt, "Process J/psi - D0", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<tableMakerJpsiHf>(cfgc)};
}
