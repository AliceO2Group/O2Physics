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
/// \brief D+-Hadrons azimuthal correlations analysis task - data-like, MC-reco and MC-Gen analyses
/// \author Shyam Kumar <shyam.kumar@cern.ch>

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/StepTHn.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TString.h>

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <memory> // std::shared_ptr
#include <string>
#include <vector>

using namespace o2;
using namespace o2::constants::math;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

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
const auto ptBinsCorrelationsVec = std::vector<double>{pTBinsCorrelations, pTBinsCorrelations + npTBinsCorrelations + 1};
const double signalRegionInnerDefault[npTBinsCorrelations] = {1.8490, 1.8490, 1.8490, 1.8490, 1.8490, 1.8490, 1.8490, 1.8490};
const double signalRegionOuterDefault[npTBinsCorrelations] = {1.8890, 1.8890, 1.8890, 1.8890, 1.8890, 1.8890, 1.8890, 1.8890};
const double sidebandLeftOuterDefault[npTBinsCorrelations] = {1.7690, 1.7690, 1.7690, 1.7690, 1.7690, 1.7690, 1.7690, 1.7690};
const double sidebandLeftInnerDefault[npTBinsCorrelations] = {1.8250, 1.8250, 1.8250, 1.8250, 1.8250, 1.8250, 1.8250, 1.8250};
const double sidebandRightInnerDefault[npTBinsCorrelations] = {1.9130, 1.9130, 1.9130, 1.9130, 1.9130, 1.9130, 1.9130, 1.9130};
const double sidebandRightOuterDefault[npTBinsCorrelations] = {1.9690, 1.9690, 1.9690, 1.9690, 1.9690, 1.9690, 1.9690, 1.9690};
const auto signalRegionInnerVec = std::vector<double>{signalRegionInnerDefault, signalRegionInnerDefault + npTBinsCorrelations};
const auto signalRegionOuterVec = std::vector<double>{signalRegionOuterDefault, signalRegionOuterDefault + npTBinsCorrelations};
const auto sidebandLeftInnerVec = std::vector<double>{sidebandLeftInnerDefault, sidebandLeftInnerDefault + npTBinsCorrelations};
const auto sidebandLeftOuterVec = std::vector<double>{sidebandLeftOuterDefault, sidebandLeftOuterDefault + npTBinsCorrelations};
const auto sidebandRightInnerVec = std::vector<double>{sidebandRightInnerDefault, sidebandRightInnerDefault + npTBinsCorrelations};
const auto sidebandRightOuterVec = std::vector<double>{sidebandRightOuterDefault, sidebandRightOuterDefault + npTBinsCorrelations};
const int npTBinsEfficiency = o2::analysis::hf_cuts_dplus_to_pi_k_pi::NBinsPt;
const std::vector<double> efficiencyDmeson(npTBinsEfficiency + 1);

/// Dplus-Hadron correlation pair filling task, from pair tables - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
struct HfTaskCorrelationDplusHadrons {
  Configurable<bool> isPromptAnalysis{"isPromptAnalysis", true, "Flag for prompt D+-hadron correlations"};
  Configurable<bool> fillHistoData{"fillHistoData", true, "Flag for filling histograms in data processes"};
  Configurable<bool> fillHistoMcRec{"fillHistoMcRec", true, "Flag for filling histograms in MC Rec processes"};
  Configurable<bool> fillHistoMcGen{"fillHistoMcGen", true, "Flag for filling histograms in MC Gen processes"};
  Configurable<bool> fillHistoMcEff{"fillHistoMcEff", true, "Flag for filling histograms in efficiency processes"};
  Configurable<int> applyEfficiency{"applyEfficiency", 1, "Flag for applying efficiency weights"};
  Configurable<bool> loadAccXEffFromCCDB{"loadAccXEffFromCCDB", false, "Flag for loading efficiency distributions from CCDB"};
  Configurable<int> selectionFlagDplus{"selectionFlagDplus", 7, "Selection Flag for D+"}; // 7 corresponds to topo+PID cuts
  Configurable<bool> selNoSameBunchPileUpColl{"selNoSameBunchPileUpColl", true, "Flag for rejecting the collisions associated with the same bunch crossing"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 1, 2}, "Indexes of ML scores to be stored. Three indexes max."};
  Configurable<std::vector<double>> mlScorePromptOrNonPromptMin{"mlScorePromptOrNonPromptMin", {0.5, 0.5, 0.5, 0.5}, "Minimum Machine learning scores for prompt or Feed-down"};
  Configurable<std::vector<double>> mlScorePromptOrNonPromptMax{"mlScorePromptOrNonPromptMax", {1.0, 1.0, 1.0, 1.0}, "Maximum Machine learning scores for prompt or Feed-down"};
  Configurable<std::vector<double>> mlScoreBkg{"mlScoreBkg", {0.5, 0.5, 0.5, 0.5}, "Machine learning scores for bkg"};
  // pT ranges for correlation plots: the default values are those embedded in hf_cuts_dplus_to_pi_k_pi (i.e. the mass pT bins), but can be redefined via json files
  Configurable<std::vector<double>> binsPtCorrelations{"binsPtCorrelations", std::vector<double>{ptBinsCorrelationsVec}, "pT bin limits for correlation plots"};
  Configurable<std::vector<double>> binsPtHadron{"binsPtHadron", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "pT bin limits for assoc particle efficiency"};
  Configurable<std::vector<double>> binsPtEfficiencyD{"binsPtEfficiencyD", std::vector<double>{o2::analysis::hf_cuts_dplus_to_pi_k_pi::vecBinsPt}, "pT bin limits for efficiency"};
  Configurable<std::vector<double>> binsPtEfficiencyHad{"binsPtEfficiencyHad", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "pT bin limits for associated particle efficiency"};
  Configurable<std::vector<float>> efficiencyD{"efficiencyD", {1., 1., 1., 1., 1., 1.}, "efficiency values for prompt D+ meson"};
  Configurable<std::vector<float>> efficiencyFdD{"efficiencyFdD", {1., 1., 1., 1., 1., 1.}, "efficiency values for beauty feed-down D+ meson"};
  Configurable<std::vector<float>> efficiencyHad{"efficiencyHad", {1., 1., 1., 1., 1., 1.}, "efficiency values for associated particles"};
  // signal and sideband region edges, to be defined via json file (initialised to empty)
  Configurable<std::vector<double>> signalRegionInner{"signalRegionInner", std::vector<double>{signalRegionInnerVec}, "Inner values of signal region vs pT"};
  Configurable<std::vector<double>> signalRegionOuter{"signalRegionOuter", std::vector<double>{signalRegionOuterVec}, "Outer values of signal region vs pT"};
  Configurable<std::vector<double>> sidebandLeftInner{"sidebandLeftInner", std::vector<double>{sidebandLeftInnerVec}, "Inner values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandLeftOuter{"sidebandLeftOuter", std::vector<double>{sidebandLeftOuterVec}, "Outer values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightInner{"sidebandRightInner", std::vector<double>{sidebandRightInnerVec}, "Inner values of right sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightOuter{"sidebandRightOuter", std::vector<double>{sidebandRightOuterVec}, "Outer values of right sideband vs pT"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 1., "max. DCA_xy of tracks"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1., "max. DCA_z of tracks"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. eta of tracks"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};
  Configurable<float> ptCandMax{"ptCandMax", 50., "max. cand pT"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.3, "min. track pT"};
  Configurable<float> ptTrackMax{"ptTrackMax", 50., "max. track pT"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> yCandGenMax{"yCandGenMax", 0.5, "max. gen. cand. rapidity"};
  Configurable<float> ptDaughterMin{"ptDaughterMin", 0.1, "min. daughter pT"};
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable debug histogram"};
  Configurable<int> nTpcCrossedRaws{"nTpcCrossedRaws", 70, "Number of crossed TPC Rows"};
  Configurable<float> cutCollPosZMc{"cutCollPosZMc", 10., "max z-vertex position for collision acceptance"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> associatedEffCcdbPath{"associatedEffCcdbPath", "", "CCDB path for associated efficiency"};
  Configurable<std::string> promptEffCcdbPath{"promptEffCcdbPath", "", "CCDB path for trigger efficiency"};
  Configurable<std::string> fdEffCcdbPath{"fdEffCcdbPath", "", "CCDB path for trigger efficiency"};
  Configurable<int64_t> timestampCcdb{"timestampCcdb", -1, "timestamp of the efficiency files used to query in CCDB"};
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  Service<ccdb::BasicCCDBManager> ccdb{};
  std::shared_ptr<TH1> mEfficiencyPrompt = nullptr;
  std::shared_ptr<TH1> mEfficiencyFD = nullptr;
  std::shared_ptr<TH1> mEfficiencyAssociated = nullptr;
  std::shared_ptr<TH1> effD = nullptr;
  int idxBdtScore = 1; // Index BDTScore 1 for Prompt and 2 for FD Analysis

  enum CandidateStep { kCandidateStepMcGenAll = 0,
                       kCandidateStepMcGenDplusToPiKPi,
                       kCandidateStepMcCandInAcceptance,
                       kCandidateStepMcDaughtersInAcceptance,
                       kCandidateStepMcReco,
                       kCandidateStepMcRecoInAcceptance,
                       kCandidateNSteps };

  using DplusHadronPair = soa::Join<aod::DplusHadronPair, aod::DplusHadronRecoInfo, aod::DplusHadronGenInfo>;
  using DplusHadronPairFullWithMl = soa::Join<aod::DplusHadronPair, aod::DplusHadronRecoInfo, aod::DplusHadronGenInfo, aod::DplusHadronMlInfo, aod::TrkRecInfoDplus>;
  using CandDplusMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi, aod::HfCand3ProngMcRec>>;
  using CandDplusMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;
  using TracksWithMc = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra, o2::aod::McTrackLabels>>; // trackFilter applied

  Filter dplusFilter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << aod::hf_cand_3prong::DecayType::DplusToPiKPi)) != static_cast<uint8_t>(0)) && aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;
  Filter trackFilter = (nabs(aod::track::eta) < etaTrackMax) && (aod::track::pt > ptTrackMin) && (aod::track::pt < ptTrackMax) && (nabs(aod::track::dcaXY) < dcaXYTrackMax) && (nabs(aod::track::dcaZ) < dcaZTrackMax);
  // configurable axis definition
  ConfigurableAxis binsMassD{"binsMassD", {200, 1.7, 2.10}, "inv. mass (#pi^{+}K^{-}#pi^{+}) (GeV/#it{c}^{2})"};
  ConfigurableAxis binsBdtScore{"binsBdtScore", {100, 0., 1.}, "Bdt output scores"};
  ConfigurableAxis binsEta{"binsEta", {100, -2., 2.}, "#it{#eta}"};
  ConfigurableAxis binsPhi{"binsPhi", {64, -PIHalf, 3. * PIHalf}, "#it{#varphi}"};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {600, 0., 8000.}, "Multiplicity as FT0M signal amplitude"};
  ConfigurableAxis binsPoolBin{"binsPoolBin", {9, 0., 9.}, "PoolBin"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    // Axis definition
    AxisSpec axisMassD = {binsMassD, "inv. mass (#pi^{+}K^{-}#pi^{+}) (GeV/#it{c}^{2})"};
    AxisSpec axisPtCorr = {(std::vector<double>)binsPtCorrelations, "#it{p}_{T}^{D} (GeV/#it{c})"};
    AxisSpec axisPtD = {(std::vector<double>)binsPtEfficiencyD, "#it{p}_{T}^{D} (GeV/#it{c})"};
    AxisSpec const axisMultFT0M = {binsMultFT0M, "MultiplicityFT0M"};
    AxisSpec axisDeltaEta = {binsEta, "#it{#eta}^{Hadron}-#it{#eta}^{D}"};
    AxisSpec axisDeltaPhi = {binsPhi, "#it{#varphi}^{Hadron}-#it{#varphi}^{D} (rad)"};
    AxisSpec axisPtHadron = {(std::vector<double>)binsPtHadron, "#it{p}_{T}^{Hadron} (GeV/#it{c})"};
    AxisSpec axisPoolBin = {binsPoolBin, "poolBin"};
    AxisSpec axisDplusPrompt = {2, -0.5, 1.5, "Prompt D+"};
    AxisSpec const axisBdtScore = {binsBdtScore, "Bdt score"};

    // Histograms for data analysis
    registry.add("hBdtScorePrompt", "D+ BDT prompt score", {HistType::kTH1F, {axisBdtScore}});
    registry.add("hBdtScoreBkg", "D+ BDT bkg score", {HistType::kTH1F, {axisBdtScore}});
    registry.add("hMassDplusVsPt", "D+ candidates massVsPt", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
    registry.add("hMassDplusVsPtWoEff", "D+ candidates massVsPt without efficiency", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
    if (fillHistoData) {
      registry.add("hDeltaEtaPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}});
      registry.add("hDeltaPhiPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}});
      registry.add("hCorrel2DPtIntSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
      registry.add("hCorrel2DVsPtSignalRegion", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hDeltaEtaPtIntSidebands", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}});
      registry.add("hDeltaPhiPtIntSidebands", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}});
      registry.add("hCorrel2DPtIntSidebands", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
      registry.add("hCorrel2DVsPtSidebands", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hDeltaEtaPtIntSidebandLeft", stringDHadron + "Left" + stringSideband + stringDeltaEta, {HistType::kTH1F, {axisDeltaEta}});
      registry.add("hDeltaPhiPtIntSidebandLeft", stringDHadron + "Left" + stringSideband + stringDeltaPhi, {HistType::kTH1F, {axisDeltaPhi}});
      registry.add("hDeltaEtaPtIntSidebandRight", stringDHadron + "Right" + stringSideband + stringDeltaEta, {HistType::kTH1F, {axisDeltaEta}});
      registry.add("hDeltaPhiPtIntSidebandRight", stringDHadron + "Right" + stringSideband + stringDeltaPhi, {HistType::kTH1F, {axisDeltaPhi}});
      registry.add("hCorrel2DVsPtSidebandLeft", stringDHadron + "Left" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtSidebandRight", stringDHadron + "Right" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});

      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegion"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebands"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandLeft"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandRight"))->Sumw2();
    }
    // Histograms for MC Reco analysis
    if (fillHistoMcRec) {
      registry.add("hMassPromptDplusVsPt", "D+ prompt candidates mass Vs Pt", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
      registry.add("hMassNonPromptDplusVsPt", "D+ non prompt candidates mass Vs Pt", {HistType::kTH2F, {{axisMassD}, {axisPtD}}});
      registry.add("hDeltaEtaPtIntSignalRegionMcRec", stringDHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}});
      registry.add("hDeltaPhiPtIntSignalRegionMcRec", stringDHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}});
      registry.add("hDeltaEtaPtIntSidebandsMcRec", stringDHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}});
      registry.add("hCorrel2DPtIntSignalRegionMcRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
      registry.add("hCorrel2DVsPtSignalRegionMcRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisDplusPrompt}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtSignalMcRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtBkgMcRec", stringDHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hDeltaPhiPtIntSidebandsMcRec", stringDHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}});
      registry.add("hCorrel2DPtIntSidebandsMcRec", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
      registry.add("hCorrel2DVsPtSidebandsMcRec", stringDHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtPhysicalPrimaryMcRec", stringDHadron + "(only true primary particles)" + stringSignal, {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisDplusPrompt}, {axisPoolBin}}});
      registry.add("hDeltaEtaPtIntSidebandLeftMcRec", stringDHadron + "Left" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTH1F, {axisDeltaEta}});
      registry.add("hDeltaPhiPtIntSidebandLeftMcRec", stringDHadron + "Left" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTH1F, {axisDeltaPhi}});
      registry.add("hCorrel2DVsPtSidebandLeftMcRec", stringDHadron + "Left" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hDeltaEtaPtIntSidebandRightMcRec", stringDHadron + "Right" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTH1F, {axisDeltaEta}});
      registry.add("hDeltaPhiPtIntSidebandRightMcRec", stringDHadron + "Right" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTH1F, {axisDeltaPhi}});
      registry.add("hCorrel2DVsPtSidebandRightMcRec", stringDHadron + "Right" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtD + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtSignalRegionPromptDplusPromptHadronMcRec", stringDHadron + "signal region PromptD - Prompt Track MC reco", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtSignalRegionNonPromptDplusNonPromptHadronMcRec", stringDHadron + " signal region PromptD - NonPrompt Track MC reco", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});

      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMcRec"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMcRec"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalMcRec"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtBkgMcRec"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandLeftMcRec"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandRightMcRec"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMcRec"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtPhysicalPrimaryMcRec"))->Sumw2();
    }
    // Histograms for MC Gen analysis
    if (fillHistoMcGen) {
      registry.add("hDeltaEtaPtIntMcGen", stringMCParticles + stringDeltaEta + "entries", {HistType::kTH1F, {axisDeltaEta}});
      registry.add("hDeltaPhiPtIntMcGen", stringMCParticles + stringDeltaPhi + "entries", {HistType::kTH1F, {axisDeltaPhi}});
      registry.add("hCorrel2DPtIntMcGen", stringMCParticles + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
      registry.add("hCorrel2DVsPtMcGen", stringMCParticles + stringDeltaPhi + stringDeltaEta + stringPtD + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtMcGenPrompt", stringDHadron + " Prompt MC Gen", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtMcGenPromptDPromptHadron", stringDHadron + "prompt D prompt h MC Gen", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtMcGenNonPromptDNonPromptHadron", stringDHadron + " non prompt D non prompt h MC Gen", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtMcGenNonPrompt", stringDHadron + " NonPrompt MC Gen", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});

      registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGen"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGenPrompt"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGenNonPrompt"))->Sumw2();
    }
    // Histograms for efficiencies
    registry.add("Efficiency/hPtCandMcRecPrompt", stringMCRecoDPrompt + stringPtD, {HistType::kTH1F, {axisPtD}});
    registry.add("Efficiency/hPtCandMcGenPrompt", stringMCGenDPrompt + stringPtD, {HistType::kTH1F, {axisPtD}});
    registry.add("Efficiency/hPtCandMcRecNonPrompt", stringMCRecoDFd + stringPtD, {HistType::kTH1F, {axisPtD}});
    registry.add("Efficiency/hPtCandMcGenNonPrompt", stringMCGenDFd + stringPtD, {HistType::kTH1F, {axisPtD}});
    registry.add("Efficiency/hPtCandMcGenDaughterInAcc", stringMCGenDFd + stringPtD, {HistType::kTH1F, {axisPtD}});

    auto hCandidates = registry.add<StepTHn>("hCandidates", "Candidate count at different steps", {HistType::kStepTHnF, {axisPtD, axisMultFT0M, {RecoDecay::OriginType::NonPrompt + 1, +RecoDecay::OriginType::None - 0.5, +RecoDecay::OriginType::NonPrompt + 0.5}}, kCandidateNSteps});
    hCandidates->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hCandidates->GetAxis(1)->SetTitle("multiplicity");
    hCandidates->GetAxis(2)->SetTitle("Charm hadron origin");

    // Loading efficiency histograms from CCDB
    if ((applyEfficiency != 0) && loadAccXEffFromCCDB) {
      ccdb->setURL(ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);

      mEfficiencyPrompt = std::shared_ptr<TH1>(ccdb->getForTimeStamp<TH1F>(promptEffCcdbPath, timestampCcdb));
      if (mEfficiencyPrompt == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", promptEffCcdbPath.value.c_str());
      }
      LOGF(info, "Loaded trigger efficiency (prompt D) histogram from %s", promptEffCcdbPath.value.c_str());

      mEfficiencyFD = std::shared_ptr<TH1>(ccdb->getForTimeStamp<TH1F>(fdEffCcdbPath, timestampCcdb));
      if (mEfficiencyFD == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", fdEffCcdbPath.value.c_str());
      }
      LOGF(info, "Loaded feed-down D meson efficiency histogram from %s", fdEffCcdbPath.value.c_str());

      mEfficiencyAssociated = std::shared_ptr<TH1>(ccdb->getForTimeStamp<TH1F>(associatedEffCcdbPath, timestampCcdb));
      if (mEfficiencyAssociated == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for associated particles from %s", associatedEffCcdbPath.value.c_str());
      }
      LOGF(info, "Loaded associated efficiency histogram from %s", associatedEffCcdbPath.value.c_str());
    }
    auto effD = isPromptAnalysis ? mEfficiencyPrompt : mEfficiencyFD;
    idxBdtScore = isPromptAnalysis ? 1 : 2;

    if (activateQA) {
      const int regionLimits = 6;
      std::string labels[regionLimits] = {"SigReg Left", "SigReg Right", "Left SB Low", "Left SB Up", "Right SB Low", "Right SB Up"};
      static const AxisSpec axisSidebandLimits = {regionLimits, 0.5, 6.5, ""};
      auto hSigSidebandLimits = registry.add<TH2>("Inputs/hSigSidebandLimits", "Signal and Sideband Limits;;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {axisSidebandLimits, {(std::vector<double>)binsPtCorrelations, "#it{p}_{T} (GeV/#it{c})"}}});
      for (int iLim = 0; iLim < regionLimits; iLim++) {
        hSigSidebandLimits->GetXaxis()->SetBinLabel(iLim + 1, labels[iLim].data());
      }
      for (size_t iPtD = 0; iPtD < binsPtCorrelations->size() - 1; iPtD++) {
        hSigSidebandLimits->SetBinContent(1, iPtD + 1, signalRegionInner->at(iPtD));
        hSigSidebandLimits->SetBinContent(2, iPtD + 1, signalRegionOuter->at(iPtD));
        hSigSidebandLimits->SetBinContent(3, iPtD + 1, sidebandLeftOuter->at(iPtD));
        hSigSidebandLimits->SetBinContent(4, iPtD + 1, sidebandLeftInner->at(iPtD));
        hSigSidebandLimits->SetBinContent(5, iPtD + 1, sidebandRightInner->at(iPtD));
        hSigSidebandLimits->SetBinContent(6, iPtD + 1, sidebandRightOuter->at(iPtD));
      }
    }
  }

  void processData(DplusHadronPairFullWithMl const& pairEntries, aod::DplusRecoInfo const& candidates)
  {
    for (const auto& candidate : candidates) {
      float const massD = candidate.mD();
      float const ptD = candidate.ptD();
      float const bdtScorePrompt = candidate.mlScorePrompt();
      float const bdtScoreNonPrompt = candidate.mlScoreNonPrompt();
      float const bdtScoreBkg = candidate.mlScoreBkg();
      int const effBinD = o2::analysis::findBin(binsPtEfficiencyD, ptD);
      float const bdtScorePromptOrNonPrompt = isPromptAnalysis ? bdtScorePrompt : bdtScoreNonPrompt;

      // reject entries outside pT ranges of interest
      if (ptD < binsPtEfficiencyD->front() || ptD > binsPtEfficiencyD->back()) {
        continue;
      }

      if (bdtScorePromptOrNonPrompt < mlScorePromptOrNonPromptMin->at(effBinD) || bdtScorePromptOrNonPrompt > mlScorePromptOrNonPromptMax->at(effBinD) || bdtScoreBkg > mlScoreBkg->at(effBinD)) {
        continue;
      }
      double efficiencyWeightD = 1.;
      if (applyEfficiency != 0) {
        efficiencyWeightD = 1. / efficiencyD->at(o2::analysis::findBin(binsPtEfficiencyD, ptD));
        if (loadAccXEffFromCCDB) {
          efficiencyWeightD = 1. / effD->GetBinContent(effD->FindBin(ptD));
        }
      }
      registry.fill(HIST("hMassDplusVsPt"), massD, ptD, efficiencyWeightD);
      registry.fill(HIST("hMassDplusVsPtWoEff"), massD, ptD);
      registry.fill(HIST("hBdtScorePrompt"), bdtScorePrompt);
      registry.fill(HIST("hBdtScoreBkg"), bdtScoreBkg);
    }

    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float const deltaPhi = pairEntry.deltaPhi();
      float const deltaEta = pairEntry.deltaEta();
      float const ptD = pairEntry.ptD();
      float const ptHadron = pairEntry.ptHadron();
      float const bdtScorePrompt = pairEntry.mlScorePrompt();
      float const bdtScoreNonPrompt = pairEntry.mlScoreNonPrompt();
      float const bdtScoreBkg = pairEntry.mlScoreBkg();
      float const trackDcaXY = pairEntry.trackDcaXY();
      float const trackDcaZ = pairEntry.trackDcaZ();
      int const trackTpcCrossedRows = pairEntry.trackTPCNClsCrossedRows();
      int const poolBin = pairEntry.poolBin();
      double const massD = pairEntry.mD();
      int const effBinD = o2::analysis::findBin(binsPtEfficiencyD, ptD);
      int const pTBinD = o2::analysis::findBin(binsPtCorrelations, ptD);
      float const bdtScorePromptOrNonPrompt = isPromptAnalysis ? bdtScorePrompt : bdtScoreNonPrompt;

      // reject entries outside pT ranges of interest
      if (ptD < binsPtEfficiencyD->front() || ptD > binsPtEfficiencyD->back()) {
        continue;
      }

      if (bdtScorePromptOrNonPrompt < mlScorePromptOrNonPromptMin->at(effBinD) || bdtScorePromptOrNonPrompt > mlScorePromptOrNonPromptMax->at(effBinD) || bdtScoreBkg > mlScoreBkg->at(effBinD)) {
        continue;
      }
      if (trackDcaXY > dcaXYTrackMax || trackDcaZ > dcaZTrackMax || trackTpcCrossedRows < nTpcCrossedRaws) {
        continue;
      }
      double efficiencyWeight = 1.;
      if (applyEfficiency != 0) {
        efficiencyWeight = 1. / (efficiencyD->at(effBinD) * efficiencyHad->at(o2::analysis::findBin(binsPtEfficiencyHad, ptHadron)));
        if (loadAccXEffFromCCDB) {
          efficiencyWeight = 1. / (effD->GetBinContent(effD->FindBin(ptD)) * mEfficiencyAssociated->GetBinContent(mEfficiencyAssociated->FindBin(ptHadron)));
        }
      }
      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if (massD > signalRegionInner->at(pTBinD) && massD < signalRegionOuter->at(pTBinD)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegion"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
      }
      // in sideband left region
      if (massD > sidebandLeftOuter->at(pTBinD) && massD < sidebandLeftInner->at(pTBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSidebandLeft"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandLeft"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandLeft"), deltaPhi, efficiencyWeight);
        registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebands"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi, efficiencyWeight);
      }
      // in sideband right region
      if (massD > sidebandRightInner->at(pTBinD) && massD < sidebandRightOuter->at(pTBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSidebandRight"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandRight"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandRight"), deltaPhi, efficiencyWeight);
        registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebands"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi, efficiencyWeight);
      }
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationDplusHadrons, processData, "Process data", false);

  /// D-Hadron correlation pair filling task, from pair tables - for MC reco-level analysis (candidates matched to true signal only, but also bkg sources are studied)
  void processMcRec(DplusHadronPairFullWithMl const& pairEntries,
                    soa::Join<aod::DplusRecoInfo, aod::DplusGenInfo> const& candidates)
  {
    for (const auto& candidate : candidates) {
      float const massD = candidate.mD();
      float const ptD = candidate.ptD();
      float const bdtScorePrompt = candidate.mlScorePrompt();
      float const bdtScoreNonPrompt = candidate.mlScoreNonPrompt();
      float const bdtScoreBkg = candidate.mlScoreBkg();
      int const effBinD = o2::analysis::findBin(binsPtEfficiencyD, ptD);
      bool const isDplusPrompt = candidate.isPrompt();
      float const bdtScorePromptOrNonPrompt = isPromptAnalysis ? bdtScorePrompt : bdtScoreNonPrompt;

      // reject entries outside pT ranges of interest
      if (ptD < binsPtEfficiencyD->front() || ptD > binsPtEfficiencyD->back()) {
        continue;
      }

      if (bdtScorePromptOrNonPrompt < mlScorePromptOrNonPromptMin->at(effBinD) || bdtScorePromptOrNonPrompt > mlScorePromptOrNonPromptMax->at(effBinD) || bdtScoreBkg > mlScoreBkg->at(effBinD)) {
        continue;
      }
      double efficiencyWeightD = 1.;
      if (applyEfficiency != 0) {
        if (isDplusPrompt) {
          efficiencyWeightD = 1. / efficiencyD->at(effBinD);
          if (loadAccXEffFromCCDB) {
            efficiencyWeightD = 1. / mEfficiencyPrompt->GetBinContent(mEfficiencyPrompt->FindBin(ptD));
          }
          registry.fill(HIST("hMassDplusVsPt"), massD, ptD, efficiencyWeightD);
          registry.fill(HIST("hMassDplusVsPtWoEff"), massD, ptD);
          registry.fill(HIST("hMassPromptDplusVsPt"), massD, ptD, efficiencyWeightD);
          registry.fill(HIST("hBdtScorePrompt"), bdtScorePrompt);
          registry.fill(HIST("hBdtScoreBkg"), bdtScoreBkg);
        } else {
          efficiencyWeightD = 1. / efficiencyFdD->at(effBinD);
          if (loadAccXEffFromCCDB) {
            efficiencyWeightD = 1. / mEfficiencyFD->GetBinContent(mEfficiencyFD->FindBin(ptD));
          }
          registry.fill(HIST("hMassDplusVsPt"), massD, ptD, efficiencyWeightD);
          registry.fill(HIST("hMassDplusVsPtWoEff"), massD, ptD);
          registry.fill(HIST("hMassNonPromptDplusVsPt"), massD, ptD, efficiencyWeightD);
          registry.fill(HIST("hBdtScorePrompt"), bdtScorePrompt);
          registry.fill(HIST("hBdtScoreBkg"), bdtScoreBkg);
        }
      }
    }

    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float const deltaPhi = pairEntry.deltaPhi();
      float const deltaEta = pairEntry.deltaEta();
      float const ptD = pairEntry.ptD();
      float const ptHadron = pairEntry.ptHadron();
      float const massD = pairEntry.mD();
      float const bdtScorePrompt = pairEntry.mlScorePrompt();
      float const bdtScoreNonPrompt = pairEntry.mlScoreNonPrompt();
      float const bdtScoreBkg = pairEntry.mlScoreBkg();
      bool const isPhysicalPrimary = pairEntry.isPhysicalPrimary();
      float const trackDcaXY = pairEntry.trackDcaXY();
      float const trackDcaZ = pairEntry.trackDcaZ();
      int const trackTpcCrossedRows = pairEntry.trackTPCNClsCrossedRows();
      bool const isDplusPrompt = pairEntry.isPrompt();
      int const originHadron = pairEntry.trackOrigin();
      int const poolBin = pairEntry.poolBin();
      int const effBinD = o2::analysis::findBin(binsPtEfficiencyD, ptD);
      int const pTBinD = o2::analysis::findBin(binsPtCorrelations, ptD);
      float const bdtScorePromptOrNonPrompt = isPromptAnalysis ? bdtScorePrompt : bdtScoreNonPrompt;

      // reject entries outside pT ranges of interest
      if (ptD < binsPtEfficiencyD->front() || ptD > binsPtEfficiencyD->back()) {
        continue;
      }

      if (bdtScorePromptOrNonPrompt < mlScorePromptOrNonPromptMin->at(effBinD) || bdtScorePromptOrNonPrompt > mlScorePromptOrNonPromptMax->at(effBinD) || bdtScoreBkg > mlScoreBkg->at(effBinD)) {
        continue;
      }
      if (trackDcaXY > dcaXYTrackMax || trackDcaZ > dcaZTrackMax || trackTpcCrossedRows < nTpcCrossedRaws) {
        continue;
      }
      double efficiencyWeight = 1.;

      if (applyEfficiency != 0) {
        if (isDplusPrompt) {
          efficiencyWeight = 1. / (efficiencyD->at(effBinD) * efficiencyHad->at(o2::analysis::findBin(binsPtEfficiencyHad, ptHadron)));
          if (loadAccXEffFromCCDB) {
            efficiencyWeight = 1. / (mEfficiencyPrompt->GetBinContent(mEfficiencyPrompt->FindBin(ptD)) * mEfficiencyAssociated->GetBinContent(mEfficiencyAssociated->FindBin(ptHadron)));
          }
        } else {
          efficiencyWeight = 1. / (efficiencyFdD->at(effBinD) * efficiencyHad->at(o2::analysis::findBin(binsPtEfficiencyHad, ptHadron)));
          if (loadAccXEffFromCCDB) {
            efficiencyWeight = 1. / (mEfficiencyFD->GetBinContent(mEfficiencyFD->FindBin(ptD)) * mEfficiencyAssociated->GetBinContent(mEfficiencyAssociated->FindBin(ptHadron)));
          }
        }
      }

      // fill correlation plots for signal/bagkground correlations
      if (pairEntry.signalStatus()) {
        registry.fill(HIST("hCorrel2DVsPtSignalMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
      } else {
        registry.fill(HIST("hCorrel2DVsPtBkgMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
      }
      // reject entries outside pT ranges of interest

      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if (massD > signalRegionInner->at(pTBinD) && massD < signalRegionOuter->at(pTBinD)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegionMcRec"), deltaPhi, deltaEta, ptD, ptHadron, static_cast<int>(isDplusPrompt), poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionMcRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionMcRec"), deltaPhi, efficiencyWeight);
        if (isPhysicalPrimary) {
          registry.fill(HIST("hCorrel2DVsPtPhysicalPrimaryMcRec"), deltaPhi, deltaEta, ptD, ptHadron, static_cast<int>(isDplusPrompt), poolBin, efficiencyWeight);
          if (isDplusPrompt && originHadron == RecoDecay::OriginType::Prompt) {
            registry.fill(HIST("hCorrel2DVsPtSignalRegionPromptDplusPromptHadronMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
          } else if (!isDplusPrompt && originHadron == RecoDecay::OriginType::NonPrompt) {
            registry.fill(HIST("hCorrel2DVsPtSignalRegionNonPromptDplusNonPromptHadronMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
          }
        }
      }
      // in sideband left region
      if (massD > sidebandLeftOuter->at(pTBinD) && massD < sidebandLeftInner->at(pTBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSidebandLeftMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandLeftMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandLeftMcRec"), deltaPhi, efficiencyWeight);
        registry.fill(HIST("hCorrel2DVsPtSidebandsMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsMcRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsMcRec"), deltaPhi, efficiencyWeight);
      }
      // in sideband right region
      if (massD > sidebandRightInner->at(pTBinD) && massD < sidebandRightOuter->at(pTBinD)) {
        registry.fill(HIST("hCorrel2DVsPtSidebandRightMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandRightMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandRightMcRec"), deltaPhi, efficiencyWeight);
        registry.fill(HIST("hCorrel2DVsPtSidebandsMcRec"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsMcRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsMcRec"), deltaPhi, efficiencyWeight);
      }
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationDplusHadrons, processMcRec, "Process MC Reco mode", true);

  /// D-Hadron correlation pair filling task, from pair tables - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGen(DplusHadronPair const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float const deltaPhi = pairEntry.deltaPhi();
      float const deltaEta = pairEntry.deltaEta();
      float const ptD = pairEntry.ptD();
      float const ptHadron = pairEntry.ptHadron();
      int const poolBin = pairEntry.poolBin();
      int const originHadron = pairEntry.trackOrigin();
      bool const isDplusPrompt = pairEntry.isPrompt();

      registry.fill(HIST("hCorrel2DVsPtMcGen"), deltaPhi, deltaEta, ptD, ptHadron, poolBin);
      registry.fill(HIST("hDeltaEtaPtIntMcGen"), deltaEta);
      registry.fill(HIST("hDeltaPhiPtIntMcGen"), deltaPhi);
      if (isDplusPrompt) {
        registry.fill(HIST("hCorrel2DVsPtMcGenPrompt"), deltaPhi, deltaEta, ptD, ptHadron, poolBin);
        if (originHadron == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("hCorrel2DVsPtMcGenPromptDPromptHadron"), deltaPhi, deltaEta, ptD, ptHadron, poolBin);
        }
      } else {
        registry.fill(HIST("hCorrel2DVsPtMcGenNonPrompt"), deltaPhi, deltaEta, ptD, ptHadron, poolBin);
        if (originHadron == RecoDecay::OriginType::NonPrompt) {
          registry.fill(HIST("hCorrel2DVsPtMcGenNonPromptDNonPromptHadron"), deltaPhi, deltaEta, ptD, ptHadron, poolBin);
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationDplusHadrons, processMcGen, "Process MC Gen mode", false);

  /// D-Hadron correlation - reconstruction and selection efficiency
  void processMcCandEfficiency(soa::Join<aod::Collisions, aod::FT0Mults, aod::EvSels> const&,
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
        if (std::abs(mcParticle.flagMcMatchGen()) == hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi) {
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
      std::vector<float> outputMl = {-1., -1., -1.};
      if (candidate.isSelDplusToPiKPi() < selectionFlagDplus) {
        continue;
      }
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
        outputMl[iclass] = candidate.mlProbDplusToPiKPi()[classMl->at(iclass)];
      }
      if (outputMl[0] > mlScoreBkg->at(o2::analysis::findBin(binsPtEfficiencyD, candidate.pt())) || outputMl[idxBdtScore] < mlScorePromptOrNonPromptMin->at(o2::analysis::findBin(binsPtEfficiencyD, candidate.pt())) || outputMl[idxBdtScore] > mlScorePromptOrNonPromptMax->at(o2::analysis::findBin(binsPtEfficiencyD, candidate.pt()))) {
        continue;
      }
      auto collision = candidate.template collision_as<soa::Join<aod::Collisions, aod::FT0Mults, aod::EvSels>>();
      if (selNoSameBunchPileUpColl && !(collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))) {
        continue;
      }
      multiplicity = collision.multFT0M();
      if (std::abs(candidate.flagMcMatchRec()) == hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi) {
        hCandidates->Fill(kCandidateStepMcReco, candidate.pt(), multiplicity, candidate.originMcRec());
        if (std::abs(HfHelper::yDplus(candidate)) <= yCandMax) {
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
