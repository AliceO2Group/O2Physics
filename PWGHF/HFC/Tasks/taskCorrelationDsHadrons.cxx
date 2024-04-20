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

/// \file taskCorrelationDsHadrons.cxx
/// \brief Ds-Hadrons azimuthal correlations analysis task - data-like, MC-reco and MC-Gen analyses
/// \author Grazia Luparello <grazia.luparello@cern.ch>
/// \author Samuele Cattaruzzi <samuele.cattaruzzi@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"

using namespace o2;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Ds-Hadron correlation pair filling task, from pair tables - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
struct HfTaskCorrelationDsHadrons {
  Configurable<bool> applyEfficiency{"applyEfficiency", false, "Flag for applying efficiency weights"};
  Configurable<int> nTpcCrossedRaws{"nTpcCrossedRaws", 70, "Number of crossed TPC Rows"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 1., "max. DCA_xy of tracks"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1., "max. DCA_z of tracks"};
  Configurable<std::vector<double>> binsPtD{"binsPtD", std::vector<double>{o2::analysis::hf_cuts_ds_to_k_k_pi::vecBinsPt}, "pT bin limits for candidate mass plots and efficiency"};
  Configurable<std::vector<double>> binsPtHadron{"binsPtHadron", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "pT bin limits for assoc particle efficiency"};
  Configurable<std::vector<double>> mlOutputPrompt{"mlScorePrompt", {0.5, 0.5, 0.5, 0.5}, "Machine learning scores for prompt"};
  Configurable<std::vector<double>> mlOutputBkg{"mlScoreBkg", {0.5, 0.5, 0.5, 0.5}, "Machine learning scores for bkg"};
  Configurable<std::vector<double>> binsPtEfficiencyD{"binsPtEfficiencyD", std::vector<double>{o2::analysis::hf_cuts_ds_to_k_k_pi::vecBinsPt}, "pT bin limits for D-meson efficiency"};
  Configurable<std::vector<double>> binsPtEfficiencyHad{"binsPtEfficiencyHad", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "pT bin limits for associated particle efficiency"};
  Configurable<std::vector<double>> efficiencyD{"efficiencyD", {1., 1., 1., 1., 1., 1.}, "efficiency values for Ds meson"};
  Configurable<std::vector<double>> efficiencyHad{"efficiencyHad", {1., 1., 1., 1., 1., 1.}, "efficiency values for associated particles"};
  Configurable<std::vector<double>> signalRegionInner{"signalRegionInner", {1.9440, 1.9440, 1.9440, 1.9440, 1.9440, 1.9440, 1.9440, 1.9440}, "Inner values of signal region vs pT"};
  Configurable<std::vector<double>> signalRegionOuter{"signalRegionOuter", {1.9920, 1.9920, 1.9920, 1.9920, 1.9920, 1.9920, 1.9920, 1.9920}, "Outer values of signal region vs pT"};
  Configurable<std::vector<double>> sidebandLeftInner{"sidebandLeftInner", {1.9360, 1.9360, 1.9360, 1.9360, 1.9360, 1.9360, 1.9360, 1.9360}, "Inner values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandLeftOuter{"sidebandLeftOuter", {1.9040, 1.9040, 1.9040, 1.9040, 1.9040, 1.9040, 1.9040, 1.9040}, "Outer values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightInner{"sidebandRightInner", {2.0000, 2.0000, 2.0000, 2.0000, 2.0000, 2.0000, 2.0000, 2.0000}, "Inner values of right sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightOuter{"sidebandRightOuter", {2.0320, 2.0320, 2.0320, 2.0320, 2.0320, 2.0320, 2.0320, 2.0320}, "Outer values of right sideband vs pT"};
  ConfigurableAxis binsMassD{"binsMassD", {200, 1.7, 2.25}, "inv. mass (K^{#pm}K^{-}#pi^{+}) (GeV/#it{c}^{2})"};
  ConfigurableAxis binsBdtScore{"binsBdtScore", {100, 0., 1.}, "Bdt output scores"};
  ConfigurableAxis binsEta{"binsEta", {100, -2., 2.}, "#it{#eta}"};
  ConfigurableAxis binsPhi{"binsPhi", {64, -PIHalf, 3. * PIHalf}, "#it{#varphi}"};
  ConfigurableAxis binsPoolBin{"binsPoolBin", {9, 0., 9.}, "PoolBin"};

  using DsHadronPairFull = soa::Join<aod::DsHadronPair, aod::DsHadronRecoInfo, aod::DsHadronGenInfo>;
  using DsHadronPairFullWithMl = soa::Join<aod::DsHadronPair, aod::DsHadronRecoInfo, aod::DsHadronGenInfo, aod::DsHadronMlInfo, aod::TrackRecoInfo>;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    AxisSpec axisMassD = {binsMassD, "inv. mass (K^{#pm}K^{-}#pi^{+}) (GeV/#it{c}^{2})"};
    AxisSpec axisDetlaEta = {binsEta, "#it{#eta}^{Hadron}-#it{#eta}^{D}"};
    AxisSpec axisDetlaPhi = {binsPhi, "#it{#varphi}^{Hadron}-#it{#varphi}^{D} (rad)"};
    AxisSpec axisPtD = {(std::vector<double>)binsPtD, "#it{p}_{T}^{D} (GeV/#it{c})"};
    AxisSpec axisPtHadron = {(std::vector<double>)binsPtHadron, "#it{p}_{T}^{Hadron} (GeV/#it{c})"};
    AxisSpec axisPoolBin = {binsPoolBin, "poolBin"};
    AxisSpec axisBdtScore = {binsBdtScore, "Bdt score"};
    AxisSpec axisDsPrompt = {2, -0.5, 1.5, "Prompt Ds"};

    // Histograms for data analysis
    registry.add("hMassDsVsPt", "Ds candidates massVsPt", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
    registry.add("hBdtScorePrompt", "Ds BDT prompt score", {HistType::kTH1F, {axisBdtScore}});
    registry.add("hBdtScoreBkg", "Ds BDT bkg score", {HistType::kTH1F, {axisBdtScore}});
    registry.add("hDeltaEtaPtIntSignalRegion", "Ds-h deltaEta signal region", {HistType::kTH1F, {axisDetlaEta}});
    registry.add("hDeltaPhiPtIntSignalRegion", "Ds-h deltaPhi signal region", {HistType::kTH1F, {axisDetlaPhi}});
    registry.add("hCorrel2DVsPtSignalRegion", "Ds-h correlations signal region", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hDeltaEtaPtIntSidebands", "Ds-h deltaEta sideband region", {HistType::kTH1F, {axisDetlaEta}});
    registry.add("hDeltaPhiPtIntSidebands", "Ds-h deltaPhi sideband region", {HistType::kTH1F, {axisDetlaPhi}});
    registry.add("hCorrel2DVsPtSidebands", "Ds-h correlations sideband region", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    // Histograms for MC Reco analysis
    registry.add("hDeltaEtaPtIntSignalRegionMcRec", "Ds-h deltaEta signal region MC reco", {HistType::kTH1F, {axisDetlaEta}});
    registry.add("hDeltaPhiPtIntSignalRegionMcRec", "Ds-h deltaPhi signal region MC reco", {HistType::kTH1F, {axisDetlaPhi}});
    registry.add("hCorrel2DVsPtSignalRegionMcRec", "Ds-h correlations signal region MC reco", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hCorrel2DVsPtPhysicalPrimaryMcRec", "Ds-h correlations signal region (only true primary particles) MC reco", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}}});
    registry.add("hCorrel2DVsPtSignalRegionMcRecPromptDivision", "Ds-h correlations signal region Prompt-NonPrompt MC reco", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisDsPrompt}, {axisPoolBin}}});
    registry.add("hDeltaEtaPtIntSidebandsMcRec", "Ds-h deltaEta sideband region MC reco", {HistType::kTH1F, {axisDetlaEta}});
    registry.add("hDeltaPhiPtIntSidebandsMcRec", "Ds-h deltaPhi sideband region MC reco", {HistType::kTH1F, {axisDetlaPhi}});
    registry.add("hCorrel2DVsPtSidebandsMcRec", "Ds-h correlations sideband region MC reco", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    // Histograms for MC Gen analysis
    registry.add("hDeltaEtaPtIntMcGen", "Ds-h deltaEta MC Gen", {HistType::kTH1F, {axisDetlaEta}});
    registry.add("hDeltaPhiPtIntMcGen", "Ds-h deltaPhi MC Gen", {HistType::kTH1F, {axisDetlaPhi}});
    registry.add("hCorrel2DVsPtMcGen", "Ds-h correlations MC Gen", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hCorrel2DVsPtMcGenPromptDivision", "Ds-h correlations Prompt-NonPrompt MC Gen", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisDsPrompt}, {axisPoolBin}}});

    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegion"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebands"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMcRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMcRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMcRecPromptDivision"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGen"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGenPromptDivision"))->Sumw2();
  }

  void processData(DsHadronPairFullWithMl const& pairEntries,
                   aod::DsCandRecoInfo const& candidates)
  {
    for (const auto& candidate : candidates) {
      float massD = candidate.mD();
      float ptD = candidate.ptD();
      float bdtScorePrompt = candidate.mlScorePrompt();
      float bdtScoreBkg = candidate.mlScoreBkg();
      int ptBinD = o2::analysis::findBin(binsPtD, ptD);

      if (bdtScorePrompt < mlOutputPrompt->at(ptBinD) || bdtScoreBkg > mlOutputBkg->at(ptBinD)) {
        continue;
      }
      double efficiencyWeightD = 1.;
      if (applyEfficiency) {
        efficiencyWeightD = 1. / efficiencyD->at(o2::analysis::findBin(binsPtEfficiencyD, ptD));
      }
      registry.fill(HIST("hMassDsVsPt"), massD, ptD, efficiencyWeightD);
      registry.fill(HIST("hBdtScorePrompt"), bdtScorePrompt);
      registry.fill(HIST("hBdtScoreBkg"), bdtScoreBkg);
    }

    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float deltaPhi = pairEntry.deltaPhi();
      float deltaEta = pairEntry.deltaEta();
      float ptD = pairEntry.ptD();
      float ptHadron = pairEntry.ptHadron();
      float massD = pairEntry.mD();
      float bdtScorePrompt = pairEntry.mlScorePrompt();
      float bdtScoreBkg = pairEntry.mlScoreBkg();
      float trackDcaXY = pairEntry.trackDcaXY();
      float trackDcaZ = pairEntry.trackDcaZ();
      int trackTpcCrossedRows = pairEntry.trackTPCNClsCrossedRows();
      int poolBin = pairEntry.poolBin();
      int ptBinD = o2::analysis::findBin(binsPtD, ptD);

      if (bdtScorePrompt < mlOutputPrompt->at(ptBinD) || bdtScoreBkg > mlOutputBkg->at(ptBinD)) {
        continue;
      }
      if (trackDcaXY > dcaXYTrackMax || trackDcaZ > dcaZTrackMax || trackTpcCrossedRows < nTpcCrossedRaws) {
        continue;
      }
      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyD->at(o2::analysis::findBin(binsPtEfficiencyD, ptD)) * efficiencyHad->at(o2::analysis::findBin(binsPtEfficiencyHad, ptHadron)));
      }

      // in signal region
      if (massD > signalRegionInner->at(ptBinD) && massD < signalRegionOuter->at(ptBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
      }
      // in sideband region
      if ((massD > sidebandLeftOuter->at(ptBinD) && massD < sidebandLeftInner->at(ptBinD)) ||
          (massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD))) {
        registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi, efficiencyWeight);
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDsHadrons, processData, "Process data", true);

  /// D-Hadron correlation pair filling task, from pair tables - for MC reco-level analysis (candidates matched to true signal only, but also bkg sources are studied)
  void processMcRec(DsHadronPairFullWithMl const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float deltaPhi = pairEntry.deltaPhi();
      float deltaEta = pairEntry.deltaEta();
      float ptD = pairEntry.ptD();
      float ptHadron = pairEntry.ptHadron();
      float massD = pairEntry.mD();
      float bdtScorePrompt = pairEntry.mlScorePrompt();
      float bdtScoreBkg = pairEntry.mlScoreBkg();
      float trackDcaXY = pairEntry.trackDcaXY();
      float trackDcaZ = pairEntry.trackDcaZ();
      int trackTpcCrossedRows = pairEntry.trackTPCNClsCrossedRows();
      int poolBin = pairEntry.poolBin();
      int statusDsPrompt = static_cast<int>(pairEntry.isPrompt());
      int ptBinD = o2::analysis::findBin(binsPtD, ptD);
      bool isPhysicalPrimary = pairEntry.isPhysicalPrimary();

      if (bdtScorePrompt < mlOutputPrompt->at(ptBinD) || bdtScoreBkg > mlOutputBkg->at(ptBinD)) {
        continue;
      }
      if (trackDcaXY > dcaXYTrackMax || trackDcaZ > dcaZTrackMax || trackTpcCrossedRows < nTpcCrossedRaws) {
        continue;
      }
      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyD->at(o2::analysis::findBin(binsPtEfficiencyD, ptD)) * efficiencyHad->at(o2::analysis::findBin(binsPtEfficiencyHad, ptHadron)));
      }
      // in signal region
      if (massD > signalRegionInner->at(ptBinD) && massD < signalRegionOuter->at(ptBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSignalRegionMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionMcRec"), deltaPhi, efficiencyWeight);
        // prompt and non-prompt division
        if (pairEntry.isSignal()) {
          registry.fill(HIST("hCorrel2DVsPtSignalRegionMcRecPromptDivision"), deltaPhi, deltaEta, ptD, ptHadron, statusDsPrompt, poolBin, efficiencyWeight);
          if (isPhysicalPrimary) {
            registry.fill(HIST("hCorrel2DVsPtPhysicalPrimaryMcRec"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
          }
        }
      }
      // in sideband region
      if (((massD > sidebandLeftOuter->at(ptBinD)) && (massD < sidebandLeftInner->at(ptBinD))) ||
          ((massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD)))) {
        registry.fill(HIST("hCorrel2DVsPtSidebandsMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsMcRec"), deltaPhi, efficiencyWeight);
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDsHadrons, processMcRec, "Process MC Reco mode", false);

  /// D-Hadron correlation pair filling task, from pair tables - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGen(DsHadronPairFull const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float deltaPhi = pairEntry.deltaPhi();
      float deltaEta = pairEntry.deltaEta();
      float ptD = pairEntry.ptD();
      float ptHadron = pairEntry.ptHadron();
      int poolBin = pairEntry.poolBin();
      int statusDsPrompt = static_cast<int>(pairEntry.isPrompt());

      registry.fill(HIST("hCorrel2DVsPtMcGen"), deltaPhi, deltaEta, ptD, ptHadron, poolBin);
      registry.fill(HIST("hCorrel2DVsPtMcGenPromptDivision"), deltaPhi, deltaEta, ptD, ptHadron, statusDsPrompt, poolBin);
      registry.fill(HIST("hDeltaEtaPtIntMcGen"), deltaEta);
      registry.fill(HIST("hDeltaPhiPtIntMcGen"), deltaPhi);
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDsHadrons, processMcGen, "Process MC Gen mode", false);

  void processDataME(DsHadronPairFullWithMl const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float deltaPhi = pairEntry.deltaPhi();
      float deltaEta = pairEntry.deltaEta();
      float ptD = pairEntry.ptD();
      float ptHadron = pairEntry.ptHadron();
      float massD = pairEntry.mD();
      float bdtScorePrompt = pairEntry.mlScorePrompt();
      float bdtScoreBkg = pairEntry.mlScoreBkg();
      float trackDcaXY = pairEntry.trackDcaXY();
      float trackDcaZ = pairEntry.trackDcaZ();
      int trackTpcCrossedRows = pairEntry.trackTPCNClsCrossedRows();
      int poolBin = pairEntry.poolBin();
      int ptBinD = o2::analysis::findBin(binsPtD, ptD);

      if (bdtScorePrompt < mlOutputPrompt->at(ptBinD) || bdtScoreBkg > mlOutputBkg->at(ptBinD)) {
        continue;
      }
      if (trackDcaXY > dcaXYTrackMax || trackDcaZ > dcaZTrackMax || trackTpcCrossedRows < nTpcCrossedRaws) {
        continue;
      }
      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyD->at(o2::analysis::findBin(binsPtEfficiencyD, ptD)) * efficiencyHad->at(o2::analysis::findBin(binsPtEfficiencyHad, ptHadron)));
      }

      // in signal region
      if (massD > signalRegionInner->at(ptBinD) && massD < signalRegionOuter->at(ptBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
      }
      // in sideband region
      if ((massD > sidebandLeftOuter->at(ptBinD) && massD < sidebandLeftInner->at(ptBinD)) ||
          (massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD))) {
        registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi, efficiencyWeight);
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDsHadrons, processDataME, "Process data ME", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationDsHadrons>(cfgc)};
}
