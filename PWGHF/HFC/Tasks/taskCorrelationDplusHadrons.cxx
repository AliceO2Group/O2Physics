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

/// \file taskCorrelationDplusHadrons.cxx
/// \author Shyam Kumar <shyam.kumar@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"

using namespace o2;
using namespace o2::constants::math;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
using DplusHadronPairFull = soa::Join<aod::DplusHadronPair, aod::DplusHadronRecoInfo>;
} // namespace o2::aod

// string definitions, used for histogram axis labels
const TString stringPtD = "#it{p}_{T}^{D} (GeV/#it{c});";
const TString stringPtHadron = "#it{p}_{T}^{Hadron} (GeV/#it{c});";
const TString stringDeltaEta = "#it{#eta}^{Hadron}-#it{#eta}^{D};";
const TString stringDeltaPhi = "#it{#varphi}^{Hadron}-#it{#varphi}^{D} (rad);";
const TString stringDHadron = "D,Hadron candidates ";
const TString stringSignal = "signal region;";
const TString stringSideband = "sidebands;";
const TString stringMCParticles = "MC gen - D,Hadron particles;";
const TString stringMCReco = "MC reco - D,Hadron candidates ";
const TString stringMCRecoDPrompt = "MC reco, prompt D+;";
const TString stringMCGenDPrompt = "MC gen, prompt D+;";
const TString stringMCRecoDFd = "MC reco, non-prompt D+;";
const TString stringMCGenDFd = "MC gen, non-prompt D+;";

const int npTBinsCorrelations = 8;
const double pTBinsCorrelations[npTBinsCorrelations + 1] = {0., 2., 4., 6., 8., 12., 16., 24., 99.};
auto pTBinsCorrelations_v = std::vector<double>{pTBinsCorrelations, pTBinsCorrelations + npTBinsCorrelations + 1};
const double signalRegionInnerDefault[npTBinsCorrelations] = {1.8490, 1.8490, 1.8490, 1.8490, 1.8490, 1.8490, 1.8490, 1.8490};
const double signalRegionOuterDefault[npTBinsCorrelations] = {1.8890, 1.8890, 1.8890, 1.8890, 1.8890, 1.8890, 1.8890, 1.8890};
const double sidebandLeftOuterDefault[npTBinsCorrelations] = {1.7690, 1.7690, 1.7690, 1.7690, 1.7690, 1.7690, 1.7690, 1.7690};
const double sidebandLeftInnerDefault[npTBinsCorrelations] = {1.8250, 1.8250, 1.8250, 1.8250, 1.8250, 1.8250, 1.8250, 1.8250};
const double sidebandRightInnerDefault[npTBinsCorrelations] = {1.9130, 1.9130, 1.9130, 1.9130, 1.9130, 1.9130, 1.9130, 1.9130};
const double sidebandRightOuterDefault[npTBinsCorrelations] = {1.9690, 1.9690, 1.9690, 1.9690, 1.9690, 1.9690, 1.9690, 1.9690};
auto signalRegionInner_v = std::vector<double>{signalRegionInnerDefault, signalRegionInnerDefault + npTBinsCorrelations};
auto signalRegionOuter_v = std::vector<double>{signalRegionOuterDefault, signalRegionOuterDefault + npTBinsCorrelations};
auto sidebandLeftInner_v = std::vector<double>{sidebandLeftInnerDefault, sidebandLeftInnerDefault + npTBinsCorrelations};
auto sidebandLeftOuter_v = std::vector<double>{sidebandLeftOuterDefault, sidebandLeftOuterDefault + npTBinsCorrelations};
auto sidebandRightInner_v = std::vector<double>{sidebandRightInnerDefault, sidebandRightInnerDefault + npTBinsCorrelations};
auto sidebandRightOuter_v = std::vector<double>{sidebandRightOuterDefault, sidebandRightOuterDefault + npTBinsCorrelations};
const int npTBinsEfficiency = o2::analysis::hf_cuts_dplus_to_pi_k_pi::nBinsPt;
std::vector<double> efficiencyDmeson(npTBinsEfficiency + 1);

/// Dplus-Hadron correlation pair filling task, from pair tables - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
struct HfTaskCorrelationDplusHadrons {
  Configurable<int> applyEfficiency{"applyEfficiency", 1, "Flag for applying efficiency weights"};
  Configurable<int> selectionFlagDplus{"selectionFlagDplus", 7, "Selection Flag for D+"}; // 7 corresponds to topo+PID cuts
  // pT ranges for correlation plots: the default values are those embedded in hf_cuts_dplus_to_pi_k_pi (i.e. the mass pT bins), but can be redefined via json files
  Configurable<std::vector<double>> binsPtCorrelations{"binsPtCorrelations", std::vector<double>{pTBinsCorrelations_v}, "pT bin limits for correlation plots"};
  Configurable<std::vector<double>> binsPtHadron{"binsPtHadron", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "pT bin limits for assoc particle efficiency"};
  Configurable<std::vector<double>> binsPtEfficiencyD{"binsPtEfficiencyD", std::vector<double>{o2::analysis::hf_cuts_dplus_to_pi_k_pi::vecBinsPt}, "pT bin limits for efficiency"};
  Configurable<std::vector<double>> binsPtEfficiencyHad{"binsPtEfficiencyHad", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "pT bin limits for associated particle efficiency"};
  Configurable<std::vector<float>> efficiencyD{"efficiencyD", {1., 1., 1., 1., 1., 1.}, "efficiency values for D+ meson"};
  Configurable<std::vector<float>> efficiencyHad{"efficiencyHad", {1., 1., 1., 1., 1., 1.}, "efficiency values for associated particles"};
  // signal and sideband region edges, to be defined via json file (initialised to empty)
  Configurable<std::vector<double>> signalRegionInner{"signalRegionInner", std::vector<double>{signalRegionInner_v}, "Inner values of signal region vs pT"};
  Configurable<std::vector<double>> signalRegionOuter{"signalRegionOuter", std::vector<double>{signalRegionOuter_v}, "Outer values of signal region vs pT"};
  Configurable<std::vector<double>> sidebandLeftInner{"sidebandLeftInner", std::vector<double>{sidebandLeftInner_v}, "Inner values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandLeftOuter{"sidebandLeftOuter", std::vector<double>{sidebandLeftOuter_v}, "Outer values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightInner{"sidebandRightInner", std::vector<double>{sidebandRightInner_v}, "Inner values of right sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightOuter{"sidebandRightOuter", std::vector<double>{sidebandRightOuter_v}, "Outer values of right sideband vs pT"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. eta of tracks"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};
  Configurable<float> ptCandMax{"ptCandMax", 50., "max. cand pT"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.3, "min. track pT"};
  Configurable<float> ptTrackMax{"ptTrackMax", 50., "max. track pT"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> yCandGenMax{"yCandGenMax", 0.5, "max. gen. cand. rapidity"};
  Configurable<float> ptDaughterMin{"ptDaughterMin", 0.1, "min. daughter pT"};
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable debug histogram"};

  ConfigurableAxis binsEta{"binsEta", {100, -2., 2.}, "#it{#eta}"};
  ConfigurableAxis binsPhi{"binsPhi", {64, -PIHalf, 3. * PIHalf}, "#it{#varphi}"};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {600, 0., 8000.}, "Multiplicity as FT0M signal amplitude"};
  ConfigurableAxis binsPoolBin{"binsPoolBin", {9, 0., 9.}, "PoolBin"};

  HfHelper hfHelper;

  enum CandidateStep { kCandidateStepMcGenAll = 0,
                       kCandidateStepMcGenDplusToPiKPi,
                       kCandidateStepMcCandInAcceptance,
                       kCandidateStepMcDaughtersInAcceptance,
                       kCandidateStepMcReco,
                       kCandidateStepMcRecoInAcceptance,
                       kCandidateNSteps };

  using CandDplusMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfCand3ProngMcRec>>;
  using CandDplusMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

  Filter dplusFilter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << aod::hf_cand_3prong::DecayType::DplusToPiKPi)) != static_cast<uint8_t>(0)) && aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    // Axis definition
    AxisSpec axisPtCorr = {(std::vector<double>)binsPtCorrelations, "#it{p}_{T}^{D} (GeV/#it{c})"};
    AxisSpec axisPtD = {(std::vector<double>)binsPtEfficiencyD, "#it{p}_{T}^{D} (GeV/#it{c})"};
    AxisSpec axisMultFT0M = {binsMultFT0M, "MultiplicityFT0M"};
    AxisSpec axisDeltaEta = {binsEta, "#it{#eta}^{Hadron}-#it{#eta}^{D}"};
    AxisSpec axisDeltaPhi = {binsPhi, "#it{#varphi}^{Hadron}-#it{#varphi}^{D} (rad)"};
    AxisSpec axisPtHadron = {(std::vector<double>)binsPtHadron, "#it{p}_{T}^{Hadron} (GeV/#it{c})"};
    AxisSpec axisPoolBin = {binsPoolBin, "poolBin"};
    AxisSpec axisDplusPrompt = {2, -0.5, 1.5, "Prompt D+"};

    registry.add("hDeltaEtaPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hDeltaEtaPtIntSidebands", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntSidebands", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntSidebands", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtSidebands", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisPoolBin}}});
    // Histograms for MC Reco analysis
    registry.add("hDeltaEtaPtIntSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hDeltaEtaPtIntSidebandsMCRec", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hCorrel2DPtIntSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtSignalRegionMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hCorrel2DVsPtSignalMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hCorrel2DVsPtBkgMCRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hDeltaPhiPtIntSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtSidebandsMCRec", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisPoolBin}}});
    // Histograms for MC Gen analysis
    registry.add("hDeltaEtaPtIntMCGen", stringMCParticles + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntMCGen", stringMCParticles + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntMCGen", stringMCParticles + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtMCGen", stringMCParticles + stringDeltaPhi + stringDeltaEta + stringPtD + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisPoolBin}}});
    // Histograms for efficiencies
    registry.add("Efficiency/hPtCandMcRecPrompt", stringMCRecoDPrompt + stringPtD, {HistType::kTH1F, {axisPtD}});
    registry.add("Efficiency/hPtCandMcGenPrompt", stringMCGenDPrompt + stringPtD, {HistType::kTH1F, {axisPtD}});
    registry.add("Efficiency/hPtCandMcRecNonPrompt", stringMCRecoDFd + stringPtD, {HistType::kTH1F, {axisPtD}});
    registry.add("Efficiency/hPtCandMcGenNonPrompt", stringMCGenDFd + stringPtD, {HistType::kTH1F, {axisPtD}});
    registry.add("Efficiency/hPtCandMcGenDaughterInAcc", stringMCGenDFd + stringPtD, {HistType::kTH1F, {axisPtD}});

    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegion"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebands"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMCRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMCRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalMCRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtBkgMCRec"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtMCGen"))->Sumw2();

    auto hCandidates = registry.add<StepTHn>("hCandidates", "Candidate count at different steps", {HistType::kStepTHnF, {axisPtD, axisMultFT0M, {RecoDecay::OriginType::NonPrompt + 1, +RecoDecay::OriginType::None - 0.5, +RecoDecay::OriginType::NonPrompt + 0.5}}, kCandidateNSteps});
    hCandidates->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hCandidates->GetAxis(1)->SetTitle("multiplicity");
    hCandidates->GetAxis(2)->SetTitle("Charm hadron origin");

    if (activateQA) {
      const int regionLimits = 6;
      std::string labels[regionLimits] = {"SigReg Left", "SigReg Right", "Left SB Low", "Left SB Up", "Right SB Low", "Right SB Up"};
      static const AxisSpec axisSidebandLimits = {regionLimits, 0.5, 6.5, ""};
      auto hSigSidebandLimits = registry.add<TH2>("Inputs/hSigSidebandLimits", "Signal and Sideband Limits;;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {axisSidebandLimits, {(std::vector<double>)binsPtCorrelations, "#it{p}_{T} (GeV/#it{c})"}}});
      for (int iLim = 0; iLim < regionLimits; iLim++) {
        hSigSidebandLimits->GetXaxis()->SetBinLabel(iLim + 1, labels[iLim].data());
      }
      for (int iPtD = 0; iPtD < binsPtCorrelations->size() - 1; iPtD++) {
        hSigSidebandLimits->SetBinContent(1, iPtD + 1, signalRegionInner->at(iPtD));
        hSigSidebandLimits->SetBinContent(2, iPtD + 1, signalRegionOuter->at(iPtD));
        hSigSidebandLimits->SetBinContent(3, iPtD + 1, sidebandLeftOuter->at(iPtD));
        hSigSidebandLimits->SetBinContent(4, iPtD + 1, sidebandLeftInner->at(iPtD));
        hSigSidebandLimits->SetBinContent(5, iPtD + 1, sidebandRightInner->at(iPtD));
        hSigSidebandLimits->SetBinContent(6, iPtD + 1, sidebandRightOuter->at(iPtD));
      }
    }
  }

  void processData(aod::DplusHadronPairFull const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      int poolBin = pairEntry.poolBin();
      double massD = pairEntry.mD();
      int effBinD = o2::analysis::findBin(binsPtEfficiencyD, ptD);
      int pTBinD = o2::analysis::findBin(binsPtCorrelations, ptD);

      // reject entries outside pT ranges of interest
      if (pTBinD < 0 || effBinD < 0) {
        continue;
      }
      if (ptHadron > 10.0) {
        ptHadron = 10.5;
      }
      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyD->at(effBinD) * efficiencyHad->at(o2::analysis::findBin(binsPtEfficiencyHad, ptHadron)));
      }
      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if (massD > signalRegionInner->at(pTBinD) && massD < signalRegionOuter->at(pTBinD)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegion"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
      }

      if ((massD > sidebandLeftOuter->at(pTBinD) && massD < sidebandLeftInner->at(pTBinD)) ||
          (massD > sidebandRightInner->at(pTBinD) && massD < sidebandRightOuter->at(pTBinD))) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebands"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi, efficiencyWeight);
      }
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationDplusHadrons, processData, "Process data", false);

  /// D-Hadron correlation pair filling task, from pair tables - for MC reco-level analysis (candidates matched to true signal only, but also bkg sources are studied)
  void processMcRec(aod::DplusHadronPairFull const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      int poolBin = pairEntry.poolBin();
      double massD = pairEntry.mD();
      int effBinD = o2::analysis::findBin(binsPtEfficiencyD, ptD);
      int pTBinD = o2::analysis::findBin(binsPtCorrelations, ptD);
      if (pTBinD < 0 || effBinD < 0) {
        continue;
      }
      if (ptHadron > 10.0) {
        ptHadron = 10.5;
      }
      double efficiencyWeight = 1.;
      if (applyEfficiency) {
        efficiencyWeight = 1. / (efficiencyD->at(effBinD) * efficiencyHad->at(o2::analysis::findBin(binsPtEfficiencyHad, ptHadron)));
      }
      // fill correlation plots for signal/bagkground correlations
      if (pairEntry.signalStatus()) {
        registry.fill(HIST("hCorrel2DVsPtSignalMCRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
      } else {
        registry.fill(HIST("hCorrel2DVsPtBkgMCRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
      }
      // reject entries outside pT ranges of interest

      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if (massD > signalRegionInner->at(pTBinD) && massD < signalRegionOuter->at(pTBinD)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegionMCRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionMCRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionMCRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionMCRec"), deltaPhi, efficiencyWeight);
      }

      if (((massD > sidebandLeftOuter->at(pTBinD)) && (massD < sidebandLeftInner->at(pTBinD))) ||
          ((massD > sidebandRightInner->at(pTBinD) && massD < sidebandRightOuter->at(pTBinD)))) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtSidebandsMCRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsMCRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsMCRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsMCRec"), deltaPhi, efficiencyWeight);
      }
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationDplusHadrons, processMcRec, "Process MC Reco mode", true);

  /// D-Hadron correlation pair filling task, from pair tables - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGen(aod::DplusHadronPair const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double deltaPhi = pairEntry.deltaPhi();
      double deltaEta = pairEntry.deltaEta();
      double ptD = pairEntry.ptD();
      double ptHadron = pairEntry.ptHadron();
      int poolBin = pairEntry.poolBin();
      // reject entries outside pT ranges of interest
      if (o2::analysis::findBin(binsPtCorrelations, ptD) < 0) {
        continue;
      }
      if (ptHadron > 10.0) {
        ptHadron = 10.5;
      }

      registry.fill(HIST("hCorrel2DVsPtMCGen"), deltaPhi, deltaEta, ptD, ptHadron, poolBin);
      registry.fill(HIST("hCorrel2DPtIntMCGen"), deltaPhi, deltaEta);
      registry.fill(HIST("hDeltaEtaPtIntMCGen"), deltaEta);
      registry.fill(HIST("hDeltaPhiPtIntMCGen"), deltaPhi);
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationDplusHadrons, processMcGen, "Process MC Gen mode", false);

  /// D-Hadron correlation - reconstruction and selection efficiency
  void processMcCandEfficiency(soa::Join<aod::Collisions, aod::FT0Mults> const&,
                               soa::Join<aod::McCollisions, aod::MultsExtraMC> const&,
                               CandDplusMcGen const& mcParticles,
                               CandDplusMcReco const& candidates,
                               aod::TracksWMc const&)
  {
    auto hCandidates = registry.get<StepTHn>(HIST("hCandidates"));

    /// Gen loop
    float multiplicity = -1.;
    for (const auto& mcParticle : mcParticles) {
      // generated candidates
      if (std::abs(mcParticle.pdgCode()) == Pdg::kDPlus) {
        auto mcCollision = mcParticle.template mcCollision_as<soa::Join<aod::McCollisions, aod::MultsExtraMC>>();
        multiplicity = mcCollision.multMCFT0A() + mcCollision.multMCFT0C(); // multFT0M = multFt0A + multFT0C
        hCandidates->Fill(kCandidateStepMcGenAll, mcParticle.pt(), multiplicity, mcParticle.originMcGen());
        if (std::abs(mcParticle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::DplusToPiKPi) {
          hCandidates->Fill(kCandidateStepMcGenDplusToPiKPi, mcParticle.pt(), multiplicity, mcParticle.originMcGen());
          auto yDplus = RecoDecay::y(mcParticle.pVector(), o2::constants::physics::MassDPlus);
          if (std::abs(yDplus) <= yCandGenMax) {
            hCandidates->Fill(kCandidateStepMcCandInAcceptance, mcParticle.pt(), multiplicity, mcParticle.originMcGen());
            if (mcParticle.originMcGen() == RecoDecay::OriginType::Prompt) {
              registry.fill(HIST("Efficiency/hPtCandMcGenPrompt"), mcParticle.pt());
            }
            if (mcParticle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
              registry.fill(HIST("Efficiency/hPtCandMcGenNonPrompt"), mcParticle.pt());
            }
          }
          bool isDaughterInAcceptance = true;
          auto daughters = mcParticle.template daughters_as<CandDplusMcGen>();
          for (const auto& daughter : daughters) {
            if (daughter.pt() < ptDaughterMin || std::abs(daughter.eta()) > etaTrackMax) {
              isDaughterInAcceptance = false;
            }
          }
          if (isDaughterInAcceptance) {
            hCandidates->Fill(kCandidateStepMcDaughtersInAcceptance, mcParticle.pt(), multiplicity, mcParticle.originMcGen());
            registry.fill(HIST("Efficiency/hPtCandMcGenDaughterInAcc"), mcParticle.pt());
          }
        }
      }
    }

    // recontructed candidates loop
    for (const auto& candidate : candidates) {
      if (candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
        continue;
      }
      if (candidate.isSelDplusToPiKPi() < selectionFlagDplus) {
        continue;
      }
      auto collision = candidate.template collision_as<soa::Join<aod::Collisions, aod::FT0Mults>>();
      multiplicity = collision.multFT0M();
      if (std::abs(candidate.flagMcMatchRec()) == 1 << aod::hf_cand_3prong::DecayType::DplusToPiKPi) {
        hCandidates->Fill(kCandidateStepMcReco, candidate.pt(), multiplicity, candidate.originMcRec());
        if (std::abs(hfHelper.yDplus(candidate)) <= yCandMax) {
          hCandidates->Fill(kCandidateStepMcRecoInAcceptance, candidate.pt(), multiplicity, candidate.originMcRec());
          if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
            registry.fill(HIST("Efficiency/hPtCandMcRecPrompt"), candidate.pt());
          }
          if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
            registry.fill(HIST("Efficiency/hPtCandMcRecNonPrompt"), candidate.pt());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDplusHadrons, processMcCandEfficiency, "Process MC for calculating candidate reconstruction efficiency", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationDplusHadrons>(cfgc)};
}
