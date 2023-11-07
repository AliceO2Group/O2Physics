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

/// \file taskCorrelationLcHadrons.cxx
/// \brief Lc-Hadrons azimuthal correlations analysis task - data-like, Mc-Reco and Mc-Gen analyses
/// \author Marianna Mazzilli <marianna.mazzilli@cern.ch>
/// \author Zhen Zhang <zhenz@cern.ch>

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

///
/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
///
double getDeltaPhi(double phiLc, double phiHadron)
{
  return RecoDecay::constrainAngle(phiHadron - phiLc, -o2::constants::math::PIHalf);
}

///
/// Returns phi of candidate/particle evaluated from x and y components of segment connecting primary and secondary vertices
///
double evaluatePhiByVertex(double xVertex1, double xVertex2, double yVertex1, double yVertex2)
{
  return RecoDecay::phi(xVertex2 - xVertex1, yVertex2 - yVertex1);
}

// string definitions, used for histogram axis labels
const TString stringPtLc = "#it{p}_{T}^{#Lambda_c} (GeV/#it{c});";
const TString stringPtHadron = "#it{p}_{T}^{Hadron} (GeV/#it{c});";
const TString stringPoolBin = "poolBin;";
const TString stringDeltaEta = "#it{#eta}^{Hadron}-#it{#eta}^{#Lambda_c};";
const TString stringDeltaPhi = "#it{#varphi}^{Hadron}-#it{#varphi}^{#Lambda_c} (rad);";
const TString stringLcHadron = "#Lambda_c,Hadron candidates ";
const TString stringSignal = "signal region;";
const TString stringSideband = "sidebands;";
const TString stringMcParticles = "Mc gen - #Lambda_c,Hadron particles;";
const TString stringMcReco = "Mc reco - #Lambda_c,Hadron candidates ";

// histogram axes definition
AxisSpec axisDeltaEta = {100, -2., 2.};
AxisSpec axisDeltaPhi = {64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf};
AxisSpec axisPtLc = {10, 0., 10.};
AxisSpec axisPtHadron = {11, 0., 11.};
AxisSpec axisPoolBin = {9, 0., 9.};

// definition of vectors for standard ptbin and invariant mass configurables
const int nPtBinsCorrelations = 8;
const double pTBinsCorrelations[nPtBinsCorrelations + 1] = {0., 2., 4., 6., 8., 12., 16., 24., 99.};
auto vecBinsPtCorrelations = std::vector<double>{pTBinsCorrelations, pTBinsCorrelations + nPtBinsCorrelations + 1};
const double signalRegionInnerDefault[nPtBinsCorrelations] = {2.269, 2.269, 2.269, 2.269, 2.269, 2.269, 2.269, 2.269};
const double signalRegionOuterDefault[nPtBinsCorrelations] = {2.309, 2.309, 2.309, 2.309, 2.309, 2.309, 2.309, 2.309};
const double sidebandLeftOuterDefault[nPtBinsCorrelations] = {2.209, 2.209, 2.209, 2.209, 2.209, 2.209, 2.209, 2.209};
const double sidebandLeftInnerDefault[nPtBinsCorrelations] = {2.249, 2.249, 2.249, 2.249, 2.249, 2.249, 2.249, 2.249};
const double sidebandRightInnerDefault[nPtBinsCorrelations] = {2.329, 2.329, 2.329, 2.329, 2.329, 2.329, 2.329, 2.329};
const double sidebandRightOuterDefault[nPtBinsCorrelations] = {2.369, 2.369, 2.369, 2.369, 2.369, 2.369, 2.369, 2.369};
auto vecSignalRegionInner = std::vector<double>{signalRegionInnerDefault, signalRegionInnerDefault + nPtBinsCorrelations};
auto vecSignalRegionOuter = std::vector<double>{signalRegionOuterDefault, signalRegionOuterDefault + nPtBinsCorrelations};
auto vecSidebandLeftInner = std::vector<double>{sidebandLeftInnerDefault, sidebandLeftInnerDefault + nPtBinsCorrelations};
auto vecSidebandLeftOuter = std::vector<double>{sidebandLeftOuterDefault, sidebandLeftOuterDefault + nPtBinsCorrelations};
auto vecSidebandRightInner = std::vector<double>{sidebandRightInnerDefault, sidebandRightInnerDefault + nPtBinsCorrelations};
auto vecSidebandRightOuter = std::vector<double>{sidebandRightOuterDefault, sidebandRightOuterDefault + nPtBinsCorrelations};
const int nPtBinsEfficiency = o2::analysis::hf_cuts_lc_to_p_k_pi::nBinsPt;
const double efficiencyLcDefault[nPtBinsEfficiency] = {};
auto vecEfficiencyLc = std::vector<double>{efficiencyLcDefault, efficiencyLcDefault + nPtBinsEfficiency};

/// Lc-Hadron correlation pair filling task, from pair tables - for real data and data-like analysis (i.e. reco-level w/o matching request via Mc truth)
struct HfTaskCorrelationLcHadrons {
  // Pt ranges for correlation plots: the default values are those embedded in hf_cuts_lc_to_p_k_pi (i.e. the mass Pt bins), but can be redefined via json files
  Configurable<int> applyEfficiency{"applyEfficiency", 1, "Flag for applying efficiency weights"};
  Configurable<std::vector<double>> binsPtCorrelations{"binsPtCorrelations", std::vector<double>{vecBinsPtCorrelations}, "Pt bin limits for correlation plots"};
  Configurable<std::vector<double>> binsPtEfficiency{"binsPtEfficiency", std::vector<double>{o2::analysis::hf_cuts_lc_to_p_k_pi::vecBinsPt}, "Pt bin limits for efficiency"};
  // signal and sideband region edges, to be defined via json file (initialised to empty)
  Configurable<std::vector<double>> signalRegionInner{"signalRegionInner", std::vector<double>{vecSignalRegionInner}, "Inner values of signal region vs Pt"};
  Configurable<std::vector<double>> signalRegionOuter{"signalRegionOuter", std::vector<double>{vecSignalRegionOuter}, "Outer values of signal region vs Pt"};
  Configurable<std::vector<double>> sidebandLeftInner{"sidebandLeftInner", std::vector<double>{vecSidebandLeftInner}, "Inner values of left sideband vs Pt"};
  Configurable<std::vector<double>> sidebandLeftOuter{"sidebandLeftOuter", std::vector<double>{vecSidebandLeftOuter}, "Outer values of left sideband vs Pt"};
  Configurable<std::vector<double>> sidebandRightInner{"sidebandRightInner", std::vector<double>{vecSidebandRightInner}, "Inner values of right sideband vs Pt"};
  Configurable<std::vector<double>> sidebandRightOuter{"sidebandRightOuter", std::vector<double>{vecSidebandRightOuter}, "Outer values of right sideband vs Pt"};
  Configurable<std::vector<double>> efficiencyLc{"efficiencyLc", std::vector<double>{vecEfficiencyLc}, "Efficiency values for Lc "};

  using LcHadronPairFull = soa::Join<aod::LcHadronPair, aod::LcHadronRecoInfo>;

  HistogramRegistry registry{
    "registry",
    {
      {"hDeltaEtaPtIntSignalRegion", stringLcHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
      {"hDeltaPhiPtIntSignalRegion", stringLcHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
      {"hCorrel2DPtIntSignalRegion", stringLcHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
      {"hCorrel2DVsPtSignalRegion", stringLcHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + stringPoolBin + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisPoolBin}}}}, // note: axes 3 and 4 (the Pt) are updated in the init()
      {"hDeltaEtaPtIntSidebands", stringLcHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
      {"hDeltaPhiPtIntSidebands", stringLcHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
      {"hCorrel2DPtIntSidebands", stringLcHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
      {"hCorrel2DVsPtSidebands", stringLcHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + stringPoolBin + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisPoolBin}}}}, // note: axes 3 and 4 (the Pt) are updated in the init()
      {"hDeltaEtaPtIntSignalRegionMcRec", stringLcHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
      {"hDeltaPhiPtIntSignalRegionMcRec", stringLcHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
      {"hCorrel2DPtIntSignalRegionMcRec", stringLcHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
      {"hCorrel2DVsPtSignalRegionMcRec", stringLcHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + stringPoolBin + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisPoolBin}}}},
      {"hCorrel2DVsPtSignalMcRec", stringLcHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + stringPoolBin + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisPoolBin}}}},
      {"hCorrel2DVsPtBkgMcRec", stringLcHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + stringPoolBin + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisPoolBin}}}},
      {"hDeltaEtaPtIntSidebandsMcRec", stringLcHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
      {"hDeltaPhiPtIntSidebandsMcRec", stringLcHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
      {"hCorrel2DPtIntSidebandsMcRec", stringLcHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
      {"hCorrel2DVsPtSidebandsMcRec", stringLcHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + stringPoolBin + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisPoolBin}}}},
      {"hDeltaEtaPtIntMcGen", stringMcParticles + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}}},
      {"hDeltaPhiPtIntMcGen", stringMcParticles + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}}},
      {"hCorrel2DPtIntMcGen", stringMcParticles + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}}},
      {"hCorrel2DVsPtMcGen", stringMcParticles + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPoolBin + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisPoolBin}}}}, // note: axes 3 and 4 (the Pt) are updated in the init()
    }};

  void init(InitContext&)
  {
    // redefinition of Pt axes for THnSparse holding correlation entries
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
  }

  void processData(LcHadronPairFull const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptLc = pairEntry.ptLc();
      double ptHadron = pairEntry.ptHadron();
      double massLc = pairEntry.mLc();
      int effBinLc = o2::analysis::findBin(binsPtEfficiency, ptLc);
      int ptBinLc = o2::analysis::findBin(binsPtCorrelations, ptLc);
      int poolBin = pairEntry.poolBin();
      // reject entries outside Pt ranges of interest
      if (ptBinLc < 0 || effBinLc < 0) {
        continue;
      }
      if (ptHadron > 10.0) {
        ptHadron = 10.5;
      }
      double efficiencyWeight = 1.;
      double efficiencyHadron = 1.; // Note: To be implemented later on
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyLc->at(effBinLc) * efficiencyHadron);
      }
      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if (massLc > signalRegionInner->at(ptBinLc) && massLc < signalRegionOuter->at(ptBinLc)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegion"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
      }

      if ((massLc > sidebandLeftOuter->at(ptBinLc) && massLc < sidebandLeftInner->at(ptBinLc)) ||
          (massLc > sidebandRightInner->at(ptBinLc) && massLc < sidebandRightOuter->at(ptBinLc))) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebands"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi, efficiencyWeight);
      }
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationLcHadrons, processData, "Process data", true);

  /// Lc-Hadron correlation pair filling task, from pair tables - for Mc reco-level analysis (candidates matched to true signal only, but also bkg sources are studied)
  void processMcRec(LcHadronPairFull const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptLc = pairEntry.ptLc();
      double ptHadron = pairEntry.ptHadron();
      double massLc = pairEntry.mLc();
      double efficiencyWeight = 1.;
      double efficiencyHadron = 1.;
      int effBinLc = o2::analysis::findBin(binsPtEfficiency, ptLc);
      int ptBinLc = o2::analysis::findBin(binsPtCorrelations, ptLc);
      int poolBin = pairEntry.poolBin();
      if (ptBinLc < 0 || effBinLc < 0) {
        continue;
      }
      if (ptHadron > 10.0) {
        ptHadron = 10.5;
      }
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyLc->at(effBinLc) * efficiencyHadron);
      }
      // fill correlation plots for signal/bagkground correlations
      if (pairEntry.signalStatus()) {
        registry.fill(HIST("hCorrel2DVsPtSignalMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, efficiencyWeight);
      } else {
        registry.fill(HIST("hCorrel2DVsPtBkgMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, efficiencyWeight);
      }
      // reject entries outside Pt ranges of interest

      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if (massLc > signalRegionInner->at(ptBinLc) && massLc < signalRegionOuter->at(ptBinLc)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegionMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionMcRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionMcRec"), deltaPhi, efficiencyWeight);
      }

      if (((massLc > sidebandLeftOuter->at(ptBinLc)) && (massLc < sidebandLeftInner->at(ptBinLc))) ||
          ((massLc > sidebandRightInner->at(ptBinLc) && massLc < sidebandRightOuter->at(ptBinLc)))) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtSidebandsMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsMcRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsMcRec"), deltaPhi, efficiencyWeight);
      }
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationLcHadrons, processMcRec, "Process Mc Reco mode", false);

  /// Lc-Hadron correlation pair filling task, from pair tables - for Mc gen-level analysis (no filter/selection, only true signal)
  void processMcGen(aod::LcHadronPair const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptLc = pairEntry.ptLc();
      double ptHadron = pairEntry.ptHadron();
      int poolBin = pairEntry.poolBin();
      // reject entries outside Pt ranges of interest
      if (o2::analysis::findBin(binsPtCorrelations, ptLc) < 0) {
        continue;
      }
      if (ptHadron > 10.0) {
        ptHadron = 10.5;
      }

      registry.fill(HIST("hCorrel2DVsPtMcGen"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin);
      registry.fill(HIST("hCorrel2DPtIntMcGen"), deltaPhi, deltaEta);
      registry.fill(HIST("hDeltaEtaPtIntMcGen"), deltaEta);
      registry.fill(HIST("hDeltaPhiPtIntMcGen"), deltaPhi);
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationLcHadrons, processMcGen, "Process Mc Gen mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationLcHadrons>(cfgc)};
}
