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
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
double getDeltaPhi(double phiD, double phiHadron)
{
  return RecoDecay::constrainAngle(phiHadron - phiD, -o2::constants::math::PIHalf);
}

/// Returns phi of candidate/particle evaluated from x and y components of segment connecting primary and secondary vertices
double evaluatePhiByVertex(double xVertex1, double xVertex2, double yVertex1, double yVertex2)
{
  return RecoDecay::phi(xVertex2 - xVertex1, yVertex2 - yVertex1);
}

// string definitions, used for histogram axis labels
const TString stringPtD = "#it{p}_{T}^{D} (GeV/#it{c});";
const TString stringPtHadron = "#it{p}_{T}^{Hadron} (GeV/#it{c});";
const TString stringPoolBin = "poolBin;";
const TString stringDeltaEta = "#it{#eta}^{Hadron}-#it{#eta}^{D};";
const TString stringDeltaPhi = "#it{#varphi}^{Hadron}-#it{#varphi}^{D} (rad);";
const TString stringDsPrompt = "Prompt Ds;";
const TString stringDHadron = "D,Hadron candidates ";
const TString stringSignal = "signal region;";
const TString stringSideband = "sidebands;";
const TString stringMCParticles = "MC gen - D,Hadron particles;";
const TString stringMCReco = "MC reco - D,Hadron candidates ";

// histogram axes definition
AxisSpec axisDetlaEta = {100, -2., 2., ""};
AxisSpec axisDetlaPhi = {64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf, ""};
AxisSpec axisPtD = {10, 0., 10., ""};
AxisSpec axisPtHadron = {11, 0., 11., ""};
AxisSpec axisDsPrompt = {2, -0.5, 1.5, ""};
AxisSpec axisPoolBin = {9, 0., 9., ""};

const int nBinsPtCorrelations = 8;
const double pTBinsCorrelations[nBinsPtCorrelations + 1] = {0., 2., 4., 6., 8., 12., 16., 24., 99.};
auto vecBinsPtCorrelations = std::vector<double>{pTBinsCorrelations, pTBinsCorrelations + nBinsPtCorrelations + 1};
const double signalRegionInnerDefault[nBinsPtCorrelations] = {1.9440, 1.9440, 1.9440, 1.9440, 1.9440, 1.9440, 1.9440, 1.9440};
const double signalRegionOuterDefault[nBinsPtCorrelations] = {1.9920, 1.9920, 1.9920, 1.9920, 1.9920, 1.9920, 1.9920, 1.9920};
const double sidebandLeftOuterDefault[nBinsPtCorrelations] = {1.9040, 1.9040, 1.9040, 1.9040, 1.9040, 1.9040, 1.9040, 1.9040};
const double sidebandLeftInnerDefault[nBinsPtCorrelations] = {1.9360, 1.9360, 1.9360, 1.9360, 1.9360, 1.9360, 1.9360, 1.9360};
const double sidebandRightInnerDefault[nBinsPtCorrelations] = {2.0000, 2.0000, 2.0000, 2.0000, 2.0000, 2.0000, 2.0000, 2.0000};
const double sidebandRightOuterDefault[nBinsPtCorrelations] = {2.0320, 2.0320, 2.0320, 2.0320, 2.0320, 2.0320, 2.0320, 2.0320};
auto vecSignalRegionInner = std::vector<double>{signalRegionInnerDefault, signalRegionInnerDefault + nBinsPtCorrelations};
auto vecSignalRegionOuter = std::vector<double>{signalRegionOuterDefault, signalRegionOuterDefault + nBinsPtCorrelations};
auto vecSidebandLeftInner = std::vector<double>{sidebandLeftInnerDefault, sidebandLeftInnerDefault + nBinsPtCorrelations};
auto vecSidebandLeftOuter = std::vector<double>{sidebandLeftOuterDefault, sidebandLeftOuterDefault + nBinsPtCorrelations};
auto vecSidebandRightInner = std::vector<double>{sidebandRightInnerDefault, sidebandRightInnerDefault + nBinsPtCorrelations};
auto vecSidebandRightOuter = std::vector<double>{sidebandRightOuterDefault, sidebandRightOuterDefault + nBinsPtCorrelations};
const int nBinsPtEfficiency = o2::analysis::hf_cuts_ds_to_k_k_pi::nBinsPt;
const double efficiencyDmesonDefault[nBinsPtEfficiency] = {};
auto vecEfficiencyDmeson = std::vector<double>{efficiencyDmesonDefault, efficiencyDmesonDefault + nBinsPtEfficiency};

/// Ds-Hadron correlation pair filling task, from pair tables - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
struct HfTaskCorrelationDsHadrons {
  Configurable<bool> applyEfficiency{"applyEfficiency", true, "Flag for applying efficiency weights"};
  // pT ranges for correlation plots: the default values are those embedded in hf_cuts_ds_to_k_k_pi (i.e. the mass pT bins), but can be redefined via json files
  Configurable<std::vector<double>> binsPtCorrelations{"binsPtCorrelations", std::vector<double>{vecBinsPtCorrelations}, "pT bin limits for correlation plots"};
  Configurable<std::vector<double>> binsPtEfficiency{"binsPtEfficiency", std::vector<double>{o2::analysis::hf_cuts_ds_to_k_k_pi::vecBinsPt}, "pT bin limits for efficiency"};
  Configurable<std::vector<double>> signalRegionInner{"signalRegionInner", std::vector<double>{vecSignalRegionInner}, "Inner values of signal region vs pT"};
  Configurable<std::vector<double>> signalRegionOuter{"signalRegionOuter", std::vector<double>{vecSignalRegionOuter}, "Outer values of signal region vs pT"};
  Configurable<std::vector<double>> sidebandLeftInner{"sidebandLeftInner", std::vector<double>{vecSidebandLeftInner}, "Inner values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandLeftOuter{"sidebandLeftOuter", std::vector<double>{vecSidebandLeftOuter}, "Outer values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightInner{"sidebandRightInner", std::vector<double>{vecSidebandRightInner}, "Inner values of right sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightOuter{"sidebandRightOuter", std::vector<double>{vecSidebandRightOuter}, "Outer values of right sideband vs pT"};
  Configurable<std::vector<double>> efficiencyD{"efficiencyD", std::vector<double>{vecEfficiencyDmeson}, "Efficiency values for D meson specie under study"};

  using DsHadronPairFull = soa::Join<aod::DsHadronPair, aod::DsHadronRecoInfo, aod::DsHadronGenInfo>;

  HistogramRegistry registry{
    "registry",
    {{"hDeltaEtaPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {axisDetlaEta}}},
     {"hDeltaPhiPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDetlaPhi}}},
     {"hCorrel2DPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDetlaPhi}, {axisDetlaEta}}}},
     {"hCorrel2DVsPtSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + stringPoolBin + "entries", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hDeltaEtaPtIntSidebands", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {axisDetlaEta}}},
     {"hDeltaPhiPtIntSidebands", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDetlaPhi}}},
     {"hCorrel2DPtIntSidebands", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDetlaPhi}, {axisDetlaEta}}}},
     {"hCorrel2DVsPtSidebands", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + stringPoolBin + "entries", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hDeltaEtaPtIntSignalRegionMcRec", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {axisDetlaEta}}},
     {"hDeltaPhiPtIntSignalRegionMcRec", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDetlaPhi}}},
     {"hDeltaEtaPtIntSidebandsMcRec", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {axisDetlaEta}}},
     {"hCorrel2DPtIntSignalRegionMcRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDetlaPhi}, {axisDetlaEta}}}},
     {"hCorrel2DVsPtSignalRegionMcRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + stringPoolBin + "entries", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}},
     {"hCorrel2DVsPtSignalRegionMcRecPromptDivision", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + stringDsPrompt + stringPoolBin + "entries", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisDsPrompt}, {axisPoolBin}}}},
     {"hCorrel2DVsPtSignalMcRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + stringPoolBin + "entries", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}},
     {"hCorrel2DVsPtBkgMcRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + stringPoolBin + "entries", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}},
     {"hDeltaPhiPtIntSidebandsMcRec", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDetlaPhi}}},
     {"hCorrel2DPtIntSidebandsMcRec", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDetlaPhi}, {axisDetlaEta}}}},
     {"hCorrel2DVsPtSidebandsMcRec", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + stringPoolBin + "entries", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}},
     {"hDeltaEtaPtIntMcGen", stringMCParticles + stringDeltaEta + "entries", {HistType::kTH1F, {axisDetlaEta}}},
     {"hDeltaPhiPtIntMcGen", stringMCParticles + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDetlaPhi}}},
     {"hCorrel2DPtIntMcGen", stringMCParticles + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDetlaPhi}, {axisDetlaEta}}}},
     {"hCorrel2DVsPtMcGen", stringMCParticles + stringDeltaPhi + stringDeltaEta + stringPtD + stringPoolBin + "entries", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}}}, // note: axes 3 and 4 (the pT) are updated in the init()
     {"hCorrel2DVsPtMcGenPromptDivision", stringMCParticles + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + stringDsPrompt + stringPoolBin + "entries", {HistType::kTHnSparseD, {{axisDetlaPhi}, {axisDetlaEta}, {axisPtD}, {axisPtHadron}, {axisDsPrompt}, {axisPoolBin}}}}}};

  void init(InitContext&)
  {
    // redefinition of pT axes for THnSparse holding correlation entries
    int nBinsPtAxis = binsPtCorrelations->size() - 1;
    const double* valuesPtAxis = binsPtCorrelations->data();

    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegion"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebands"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegion"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebands"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMcRec"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMcRec"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMcRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMcRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalMcRec"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalMcRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtBkgMcRec"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtBkgMcRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGen"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGen"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGenPromptDivision"))->GetAxis(2)->Set(nBinsPtAxis, valuesPtAxis);
    registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGenPromptDivision"))->Sumw2();
  }

  void processData(DsHadronPairFull const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      int poolBin = pairEntry.poolBin();
      double massD = pairEntry.mD();
      int effBinD = o2::analysis::findBin(binsPtEfficiency, ptD);
      int pTBinD = o2::analysis::findBin(binsPtCorrelations, ptD);
      // reject entries outside pT ranges of interest
      if (pTBinD < 0 || effBinD < 0) {
        continue;
      }
      double efficiencyWeight = 1.;
      double efficiencyHadron = 1.; // Note: To be implemented later on
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyD->at(effBinD) * efficiencyHadron);
      }
      // in signal region
      if (massD > signalRegionInner->at(pTBinD) && massD < signalRegionOuter->at(pTBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegion"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
      }
      // in sideband region
      if ((massD > sidebandLeftOuter->at(pTBinD) && massD < sidebandLeftInner->at(pTBinD)) ||
          (massD > sidebandRightInner->at(pTBinD) && massD < sidebandRightOuter->at(pTBinD))) {
        registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebands"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi, efficiencyWeight);
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDsHadrons, processData, "Process data", true);

  /// D-Hadron correlation pair filling task, from pair tables - for MC reco-level analysis (candidates matched to true signal only, but also bkg sources are studied)
  void processMcRec(DsHadronPairFull const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      int poolBin = pairEntry.poolBin();
      double massD = pairEntry.mD();
      int statusDsPrompt = static_cast<int>(pairEntry.isPrompt());
      int effBinD = o2::analysis::findBin(binsPtEfficiency, ptD);
      int pTBinD = o2::analysis::findBin(binsPtCorrelations, ptD);
      // reject entries outside pT ranges of interest
      if (pTBinD < 0 || effBinD < 0) {
        continue;
      }
      double efficiencyWeight = 1.;
      double efficiencyHadron = 1.; // Note: To be implemented later on
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyD->at(effBinD) * efficiencyHadron);
      }
      // fill correlation plots for signal/bagkground correlations
      if (pairEntry.isSignal()) {
        registry.fill(HIST("hCorrel2DVsPtSignalMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
      } else {
        registry.fill(HIST("hCorrel2DVsPtBkgMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
      }
      // in signal region
      if (massD > signalRegionInner->at(pTBinD) && massD < signalRegionOuter->at(pTBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSignalRegionMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionMcRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionMcRec"), deltaPhi, efficiencyWeight);
        // prompt and non-prompt division
        if (pairEntry.isSignal()) {
          registry.fill(HIST("hCorrel2DVsPtSignalRegionMcRecPromptDivision"), deltaPhi, deltaEta, ptD, ptHadron, statusDsPrompt, poolBin, efficiencyWeight);
        }
      }
      // in sideband region
      if (((massD > sidebandLeftOuter->at(pTBinD)) && (massD < sidebandLeftInner->at(pTBinD))) ||
          ((massD > sidebandRightInner->at(pTBinD) && massD < sidebandRightOuter->at(pTBinD)))) {
        registry.fill(HIST("hCorrel2DVsPtSidebandsMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsMcRec"), deltaPhi, deltaEta, efficiencyWeight);
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
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      int poolBin = pairEntry.poolBin();
      int statusDsPrompt = static_cast<int>(pairEntry.isPrompt());
      // reject entries outside pT ranges of interest
      if (o2::analysis::findBin(binsPtCorrelations, ptD) < 0) {
        continue;
      }

      registry.fill(HIST("hCorrel2DVsPtMcGen"), deltaPhi, deltaEta, ptD, ptHadron, poolBin);
      registry.fill(HIST("hCorrel2DVsPtMcGenPromptDivision"), deltaPhi, deltaEta, ptD, ptHadron, statusDsPrompt, poolBin);
      registry.fill(HIST("hCorrel2DPtIntMcGen"), deltaPhi, deltaEta);
      registry.fill(HIST("hDeltaEtaPtIntMcGen"), deltaEta);
      registry.fill(HIST("hDeltaPhiPtIntMcGen"), deltaPhi);
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDsHadrons, processMcGen, "Process MC Gen mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationDsHadrons>(cfgc)};
}
