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

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/HFC/Utils/utilsCorrelations.h"
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
using namespace o2::analysis::hf_correlations;

// string definitions, used for histogram axis labels
const TString stringPtLc = "#it{p}_{T}^{#Lambda_c} (GeV/#it{c});";
const TString stringPtHadron = "#it{p}_{T}^{Hadron} (GeV/#it{c});";
const TString stringSign = "pairSign;";
const TString stringMass = "M_{pK#pi} (GeV/#it{c^2});";
const TString stringDeltaEta = "#it{#eta}^{Hadron}-#it{#eta}^{#Lambda_c};";
const TString stringDeltaPhi = "#it{#varphi}^{Hadron}-#it{#varphi}^{#Lambda_c} (rad);";
const TString stringLcHadron = "#Lambda_c,Hadron candidates ";
const TString stringSignal = "signal region;";
const TString stringSignMass = "sign and invMass;";
const TString stringSideband = "sidebands;";
const TString stringMcParticles = "MC gen - #Lambda_c,Hadron particles;";
const TString stringMcReco = "MC reco - #Lambda_c,Hadron candidates ";
const TString stringMcRecoLcPrompt = "MC reco, prompt #Lambda_c;";
const TString stringMcGenLcPrompt = "MC gen, prompt #Lambda_c;";
const TString stringMcRecoLcFd = "MC reco, non-prompt #Lambda_c;";
const TString stringMcGenLcFd = "MC gen, non-prompt #Lambda_c;";

// definition of vectors for standard ptbin and invariant mass configurables
const int nPtBinsCorrelations = 8;
const double pTBinsCorrelations[nPtBinsCorrelations + 1] = {0., 2., 4., 6., 8., 12., 16., 24., 99.};
const auto vecBinsPtCorrelations = std::vector<double>{pTBinsCorrelations, pTBinsCorrelations + nPtBinsCorrelations + 1};
const double signalRegionInnerDefault[nPtBinsCorrelations] = {2.269, 2.269, 2.269, 2.269, 2.269, 2.269, 2.269, 2.269};
const double signalRegionOuterDefault[nPtBinsCorrelations] = {2.309, 2.309, 2.309, 2.309, 2.309, 2.309, 2.309, 2.309};
const double sidebandLeftOuterDefault[nPtBinsCorrelations] = {2.209, 2.209, 2.209, 2.209, 2.209, 2.209, 2.209, 2.209};
const double sidebandLeftInnerDefault[nPtBinsCorrelations] = {2.249, 2.249, 2.249, 2.249, 2.249, 2.249, 2.249, 2.249};
const double sidebandRightInnerDefault[nPtBinsCorrelations] = {2.329, 2.329, 2.329, 2.329, 2.329, 2.329, 2.329, 2.329};
const double sidebandRightOuterDefault[nPtBinsCorrelations] = {2.369, 2.369, 2.369, 2.369, 2.369, 2.369, 2.369, 2.369};
const auto vecSignalRegionInner = std::vector<double>{signalRegionInnerDefault, signalRegionInnerDefault + nPtBinsCorrelations};
const auto vecSignalRegionOuter = std::vector<double>{signalRegionOuterDefault, signalRegionOuterDefault + nPtBinsCorrelations};
const auto vecSidebandLeftInner = std::vector<double>{sidebandLeftInnerDefault, sidebandLeftInnerDefault + nPtBinsCorrelations};
const auto vecSidebandLeftOuter = std::vector<double>{sidebandLeftOuterDefault, sidebandLeftOuterDefault + nPtBinsCorrelations};
const auto vecSidebandRightInner = std::vector<double>{sidebandRightInnerDefault, sidebandRightInnerDefault + nPtBinsCorrelations};
const auto vecSidebandRightOuter = std::vector<double>{sidebandRightOuterDefault, sidebandRightOuterDefault + nPtBinsCorrelations};

/// Lc-Hadron correlation pair filling task, from pair tables - for real data and data-like analysis (i.e. reco-level w/o matching request via Mc truth)
struct HfTaskCorrelationLcHadrons {
  Configurable<bool> fillHistoData{"fillHistoData", true, "Flag for filling histograms in data processes"};
  Configurable<bool> fillHistoMcRec{"fillHistoMcRec", true, "Flag for filling histograms in MC Rec processes"};
  Configurable<bool> fillHistoMcGen{"fillHistoMcGen", true, "Flag for filling histograms in MC Gen processes"};
  Configurable<bool> fillHistoMcEff{"fillHistoMcEff", true, "Flag for filling histograms in efficiency processes"};
  Configurable<int> storeMass{"storeMass", 1, "Flag for storing mass information"};
  Configurable<int> applyEfficiency{"applyEfficiency", 1, "Flag for applying efficiency weights"};
  Configurable<bool> loadAccXEffFromCCDB{"loadAccXEffFromCCDB", false, "Flag for loading efficiency distributions from CCDB"};
  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<bool> selNoSameBunchPileUpColl{"selNoSameBunchPileUpColl", true, "Flag for rejecting the collisions associated with the same bunch crossing"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 1, 2}, "Indexes of ML scores to be stored. Three indexes max."};
  Configurable<std::vector<double>> mlOutputPrompt{"mlOutputPrompt", {0.5, 0.5, 0.5, 0.5}, "Machine learning scores for prompt"};
  Configurable<std::vector<double>> mlOutputBkg{"mlOutputBkg", {0.5, 0.5, 0.5, 0.5}, "Machine learning scores for bkg"};
  // Pt ranges for correlation plots: the default values are those embedded in hf_cuts_lc_to_p_k_pi (i.e. the mass Pt bins), but can be redefined via json files
  Configurable<std::vector<double>> binsPtCorrelations{"binsPtCorrelations", std::vector<double>{vecBinsPtCorrelations}, "Pt bin limits for correlation plots"};
  Configurable<std::vector<double>> binsPtHadron{"binsPtHadron", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "Pt bin limits for assoc particle efficiency"};
  Configurable<std::vector<double>> binsPtEfficiencyLc{"binsPtEfficiencyLc", std::vector<double>{o2::analysis::hf_cuts_lc_to_p_k_pi::vecBinsPt}, "Pt bin limits for efficiency"};
  Configurable<std::vector<double>> binsPtEfficiencyHad{"binsPtEfficiencyHad", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "pT bin limits for associated particle efficiency"};
  Configurable<std::vector<float>> efficiencyLc{"efficiencyLc", {1., 1., 1., 1., 1., 1.}, "efficiency values for prompt Lc"};
  Configurable<std::vector<float>> efficiencyFdLc{"efficiencyFdLc", {1., 1., 1., 1., 1., 1.}, "efficiency values for beauty feed-down Lc"};
  Configurable<std::vector<float>> efficiencyHad{"efficiencyHad", {1., 1., 1., 1., 1., 1.}, "efficiency values for associated particles"};
  // Signal and sideband region edges, to be defined via json file (initialised to empty)
  Configurable<std::vector<double>> signalRegionInner{"signalRegionInner", std::vector<double>{vecSignalRegionInner}, "Inner values of signal region vs Pt"};
  Configurable<std::vector<double>> signalRegionOuter{"signalRegionOuter", std::vector<double>{vecSignalRegionOuter}, "Outer values of signal region vs Pt"};
  Configurable<std::vector<double>> sidebandLeftInner{"sidebandLeftInner", std::vector<double>{vecSidebandLeftInner}, "Inner values of left sideband vs Pt"};
  Configurable<std::vector<double>> sidebandLeftOuter{"sidebandLeftOuter", std::vector<double>{vecSidebandLeftOuter}, "Outer values of left sideband vs Pt"};
  Configurable<std::vector<double>> sidebandRightInner{"sidebandRightInner", std::vector<double>{vecSidebandRightInner}, "Inner values of right sideband vs Pt"};
  Configurable<std::vector<double>> sidebandRightOuter{"sidebandRightOuter", std::vector<double>{vecSidebandRightOuter}, "Outer values of right sideband vs Pt"};
  Configurable<bool> isTowardTransverseAway{"isTowardTransverseAway", false, "Divide into three regions: toward, transverse, and away"};
  Configurable<double> leadingParticlePtMin{"leadingParticlePtMin", 0., "Min for leading particle pt"};
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
  Configurable<bool> useCentrality{"useCentrality", false, "Flag for centrality dependent analyses"};
  // sign and invMasss
  Configurable<bool> fillSignAndMass{"fillSignAndMass", false, "flag to select Lc-h corr with Lc invarient mass and sign of pairs"};
  Configurable<bool> calSign{"calSign", false, "flag to calculate sign of pairs"};
  Configurable<bool> fillSign{"fillSign", false, "flag to fill sign of pairs in ThnSparse"};

  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> associatedEffCcdbPath{"associatedEffCcdbPath", "", "CCDB path for associated efficiency"};
  Configurable<std::string> promptEffCcdbPath{"promptEffCcdbPath", "", "CCDB path for trigger efficiency"};
  Configurable<std::string> fdEffCcdbPath{"fdEffCcdbPath", "", "CCDB path for trigger efficiency"};
  Configurable<int64_t> timestampCcdb{"timestampCcdb", -1, "timestamp of the efficiency files used to query in CCDB"};
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  std::shared_ptr<TH1> mEfficiencyPrompt = nullptr;
  std::shared_ptr<TH1> mEfficiencyFD = nullptr;
  std::shared_ptr<TH1> mEfficiencyAssociated = nullptr;

  Service<ccdb::BasicCCDBManager> ccdb;

  enum CandidateStep { kCandidateStepMcGenAll = 0,
                       kCandidateStepMcGenLcToPKPi,
                       kCandidateStepMcCandInAcceptance,
                       kCandidateStepMcDaughtersInAcceptance,
                       kCandidateStepMcReco,
                       kCandidateStepMcRecoInAcceptance,
                       kCandidateNSteps };

  using LcHadronPair = soa::Join<aod::LcHadronPair, aod::LcHadronRecoInfo, aod::LcHadronGenInfo>;
  using LcHadronPairFullWithMl = soa::Join<aod::LcHadronPair, aod::LcHadronRecoInfo, aod::LcHadronGenInfo, aod::LcHadronMlInfo, aod::TrkRecInfoLc>;
  using CandLcMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfMlLcToPKPi, aod::HfCand3ProngMcRec>>;
  using CandLcMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;
  using TracksWithMc = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra, o2::aod::McTrackLabels>>; // trackFilter applied

  Filter lcFilter = ((o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) != static_cast<uint8_t>(0)) && (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLc || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLc);
  Filter trackFilter = (nabs(aod::track::eta) < etaTrackMax) && (aod::track::pt > ptTrackMin) && (aod::track::pt < ptTrackMax) && (nabs(aod::track::dcaXY) < dcaXYTrackMax) && (nabs(aod::track::dcaZ) < dcaZTrackMax);

  // configurable axis definition
  ConfigurableAxis binsMassLc{"binsMassLc", {200, 1.98, 2.58}, "inv. mass (p K #pi) (GeV/#it{c}^{2})"};
  ConfigurableAxis binsBdtScore{"binsBdtScore", {100, 0., 1.}, "Bdt output scores"};
  ConfigurableAxis binsEta{"binsEta", {100, -2., 2.}, "#it{#eta}"};
  ConfigurableAxis binsPhi{"binsPhi", {64, -PIHalf, 3. * PIHalf}, "#it{#varphi}"};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {600, 0., 8000.}, "Multiplicity as FT0M signal amplitude"};
  ConfigurableAxis binsPoolBin{"binsPoolBin", {9, 0., 9.}, "PoolBin"};
  ConfigurableAxis binsCentFt0m{"binsCentFt0m", {100, 0., 100.}, "Centrality percentile (FT0M)"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    // Axis definition
    AxisSpec axisMassLc = {binsMassLc, "inv. mass (p K #pi) (GeV/#it{c}^{2})"};
    AxisSpec axisPtCorr = {(std::vector<double>)binsPtCorrelations, "#it{p}_{T}^{#Lambda_c} (GeV/#it{c})"};
    AxisSpec axisPtLc = {(std::vector<double>)binsPtEfficiencyLc, "#it{p}_{T}^{#Lambda_c} (GeV/#it{c})"};
    AxisSpec const axisMultFT0M = {binsMultFT0M, "MultiplicityFT0M"};
    AxisSpec axisDeltaEta = {binsEta, "#it{#eta}^{Hadron}-#it{#eta}^{#Lambda_c}"};
    AxisSpec axisDeltaPhi = {binsPhi, "#it{#varphi}^{Hadron}-#it{#varphi}^{#Lambda_c} (rad)"};
    AxisSpec axisPtHadron = {(std::vector<double>)binsPtHadron, "#it{p}_{T}^{Hadron} (GeV/#it{c})"};
    AxisSpec axisPoolBin = {binsPoolBin, "poolBin"};
    AxisSpec axisLcPrompt = {2, -0.5, 1.5, "Prompt #Lambda_c"};
    AxisSpec const axisBdtScore = {binsBdtScore, "Bdt score"};
    AxisSpec axisCorrelationState = {2, 0., 2., ""};
    AxisSpec axisSignPair = {4, 1., 5.};
    AxisSpec axisCentFT0M = {binsCentFt0m, "Centrality percentile (FT0M)"};

    // Histograms for data analysis
    registry.add("hBdtScorePrompt", "Lc BDT prompt score", {HistType::kTH1D, {axisBdtScore}});
    registry.add("hBdtScoreBkg", "Lc BDT bkg score", {HistType::kTH1D, {axisBdtScore}});
    registry.add("hMassLcVsPt", "Lc candidates massVsPt", {HistType::kTH2F, {{axisMassLc}, {axisPtLc}}});
    registry.add("hMassLcVsPtWoEff", "Lc candidates massVsPt without efficiency", {HistType::kTH2F, {{axisMassLc}, {axisPtLc}}});
    if (fillHistoData) {
      registry.add("hDeltaEtaPtIntSignalRegion", stringLcHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1D, {axisDeltaEta}});
      registry.add("hDeltaPhiPtIntSignalRegion", stringLcHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1D, {axisDeltaPhi}});
      registry.add("hCorrel2DPtIntSignalRegion", stringLcHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
      registry.add("hDeltaEtaPtIntSidebands", stringLcHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1D, {axisDeltaEta}});
      registry.add("hDeltaPhiPtIntSidebands", stringLcHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1D, {axisDeltaPhi}});
      registry.add("hCorrel2DPtIntSidebands", stringLcHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
      registry.add("hCorrel2DVsPtSidebands", stringLcHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisPoolBin}, {axisCentFT0M}}});
      registry.add("hDeltaEtaPtIntSidebandLeft", stringLcHadron + "Left" + stringSideband + stringDeltaEta, {HistType::kTH1D, {axisDeltaEta}});
      registry.add("hDeltaPhiPtIntSidebandLeft", stringLcHadron + "Left" + stringSideband + stringDeltaPhi, {HistType::kTH1D, {axisDeltaPhi}});
      registry.add("hDeltaEtaPtIntSidebandRight", stringLcHadron + "Right" + stringSideband + stringDeltaEta, {HistType::kTH1D, {axisDeltaEta}});
      registry.add("hDeltaPhiPtIntSidebandRight", stringLcHadron + "Right" + stringSideband + stringDeltaPhi, {HistType::kTH1D, {axisDeltaPhi}});

      if (!fillSign) {
        registry.add("hCorrel2DVsPtSidebandLeft", stringLcHadron + "Left" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisPoolBin}, {axisCentFT0M}}});
        registry.add("hCorrel2DVsPtSidebandRight", stringLcHadron + "Right" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisPoolBin}, {axisCentFT0M}}});
        registry.add("hCorrel2DVsPtSignalRegion", stringLcHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisPoolBin}, {axisCentFT0M}}});
        registry.add("hCorrel2DVsPtGlobalRegion", stringLcHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisPoolBin}, {axisCentFT0M}, {axisMassLc}}});
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandLeft"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandRight"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegion"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtGlobalRegion"))->Sumw2();
      } else {
        registry.add("hCorrel2DVsPtSignSidebandLeft", stringLcHadron + "Left" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + stringSign + "entries", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisSignPair}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtSignSidebandRight", stringLcHadron + "Right" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + stringSign + "entries", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisSignPair}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtSignSignalRegion", stringLcHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + stringSign + "entries", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisSignPair}, {axisPoolBin}}});
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignSidebandLeft"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignSidebandRight"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignSignalRegion"))->Sumw2();
      }
      // Toward Transverse Away
      registry.add("hToward", "Toward invmass; ptLc; correlationState;entries", {HistType::kTH3F, {{axisMassLc}, {axisPtLc}, {axisCorrelationState}}});
      registry.add("hTransverse", "Transverse invmass; ptLc; correlationState;entries", {HistType::kTH3F, {{axisMassLc}, {axisPtLc}, {axisCorrelationState}}});
      registry.add("hAway", "Away invmass; ptLc; correlationState;entries", {HistType::kTH3F, {{axisMassLc}, {axisPtLc}, {axisCorrelationState}}});

      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebands"))->Sumw2();

      if (fillSignAndMass) {
        registry.add("hCorrel2DVsPtSignMass", stringLcHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + stringSignMass + "entries", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisMassLc}, {axisSignPair}, {axisPoolBin}}});
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignMass"))->Sumw2();
      }
    }
    // Histograms for MC Reco analysis
    if (fillHistoMcRec) {
      registry.add("hMassPromptLcVsPt", "Lc prompt candidates mass Vs Pt", {HistType::kTH2F, {{axisMassLc}, {axisPtLc}}});
      registry.add("hMassNonPromptLcVsPt", "Lc non prompt candidates mass Vs Pt", {HistType::kTH2F, {{axisMassLc}, {axisPtLc}}});
      registry.add("hDeltaEtaPtIntSignalRegionMcRec", stringLcHadron + stringSignal + stringDeltaEta + "entries", {HistType::kTH1D, {axisDeltaEta}});
      registry.add("hDeltaPhiPtIntSignalRegionMcRec", stringLcHadron + stringSignal + stringDeltaPhi + "entries", {HistType::kTH1D, {axisDeltaPhi}});
      registry.add("hDeltaEtaPtIntSidebandsMcRec", stringLcHadron + stringSideband + stringDeltaEta + "entries", {HistType::kTH1D, {axisDeltaEta}});
      registry.add("hCorrel2DPtIntSignalRegionMcRec", stringLcHadron + stringSignal + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
      registry.add("hDeltaPhiPtIntSidebandsMcRec", stringLcHadron + stringSideband + stringDeltaPhi + "entries", {HistType::kTH1D, {axisDeltaPhi}});
      registry.add("hCorrel2DPtIntSidebandsMcRec", stringLcHadron + stringSideband + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
      registry.add("hCorrel2DVsPtSidebandsMcRec", stringLcHadron + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisPoolBin}}});
      registry.add("hCorrel2DVsPtPhysicalPrimaryMcRec", stringLcHadron + "(only true primary particles)" + stringSignal, {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisLcPrompt}, {axisPoolBin}}});
      registry.add("hDeltaEtaPtIntSidebandLeftMcRec", stringLcHadron + "Left" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTH1D, {axisDeltaEta}});
      registry.add("hDeltaPhiPtIntSidebandLeftMcRec", stringLcHadron + "Left" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTH1D, {axisDeltaPhi}});
      registry.add("hDeltaEtaPtIntSidebandRightMcRec", stringLcHadron + "Right" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTH1D, {axisDeltaEta}});
      registry.add("hDeltaPhiPtIntSidebandRightMcRec", stringLcHadron + "Right" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTH1D, {axisDeltaPhi}});

      if (!fillSign) {
        registry.add("hCorrel2DVsPtSidebandLeftMcRec", stringLcHadron + "Left" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtSidebandRightMcRec", stringLcHadron + "Right" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtSignalRegionPromptLcPromptHadronMcRec", stringLcHadron + "signal region PromptLc - Prompt Track MC reco", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtSignalRegionNonPromptLcNonPromptHadronMcRec", stringLcHadron + " signal region PromptLc - NonPrompt Track MC reco", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtSignalMcRec", stringLcHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtSignalRegionMcRec", stringLcHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisLcPrompt}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtBkgMcRec", stringLcHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisPoolBin}}});

        registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandLeftMcRec"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandRightMcRec"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalMcRec"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionPromptLcPromptHadronMcRec"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionNonPromptLcNonPromptHadronMcRec"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionMcRec"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtBkgMcRec"))->Sumw2();

      } else {
        registry.add("hCorrel2DVsPtSignSidebandLeftMcRec", stringLcHadron + "Left" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + stringSign + "entries", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisSignPair}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtSignSidebandRightMcRec", stringLcHadron + "Right" + stringSideband + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + stringSign + "entries", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisSignPair}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtSignSignalRegionPromptLcPromptHadronMcRec", stringLcHadron + "signal region PromptLc - Prompt Track MC reco", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisSignPair}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtSignSignalRegionNonPromptLcNonPromptHadronMcRec", stringLcHadron + " signal region PromptLc - NonPrompt Track MC reco", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisSignPair}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtSignSignalMcRec", stringLcHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + stringSign + "entries", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisSignPair}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtSignSignalRegionMcRec", stringLcHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + stringSign + "entries", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisSignPair}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtSignBkgMcRec", stringLcHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + stringSign + "entries", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisSignPair}, {axisPoolBin}}});

        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignSidebandLeftMcRec"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignSidebandRightMcRec"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignSignalMcRec"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignSignalRegionPromptLcPromptHadronMcRec"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignSignalRegionNonPromptLcNonPromptHadronMcRec"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignSignalRegionMcRec"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignBkgMcRec"))->Sumw2();
      }

      if (fillSignAndMass) {
        registry.add("hCorrel2DVsPtSignMassMcRec", stringLcHadron + stringSignal + stringDeltaPhi + stringDeltaEta + stringPtLc + stringPtHadron + stringSignMass + "entries", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisMassLc}, {axisSignPair}, {axisPoolBin}}});
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignMassMcRec"))->Sumw2();
      }

      // Toward Transverse Away for McRec
      registry.add("hTowardRec", "Toward invmass; ptLc; correlationState;entries", {HistType::kTH3F, {{axisMassLc}, {axisPtLc}, {axisCorrelationState}}});
      registry.add("hTransverseRec", "Transverse invmass; ptLc; correlationState;entries", {HistType::kTH3F, {{axisMassLc}, {axisPtLc}, {axisCorrelationState}}});
      registry.add("hAwayRec", "Away invmass; ptLc; correlationState;entries", {HistType::kTH3F, {{axisMassLc}, {axisPtLc}, {axisCorrelationState}}});

      registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsMcRec"))->Sumw2();
      registry.get<THnSparse>(HIST("hCorrel2DVsPtPhysicalPrimaryMcRec"))->Sumw2();
    }
    // Histograms for MC Gen analysis
    if (fillHistoMcGen) {
      registry.add("hDeltaEtaPtIntMcGen", stringMcParticles + stringDeltaEta + "entries", {HistType::kTH1D, {axisDeltaEta}});
      registry.add("hDeltaPhiPtIntMcGen", stringMcParticles + stringDeltaPhi + "entries", {HistType::kTH1D, {axisDeltaPhi}});
      registry.add("hCorrel2DPtIntMcGen", stringMcParticles + stringDeltaPhi + stringDeltaEta + "entries", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});

      if (!fillSign) {
        registry.add("hCorrel2DVsPtMcGen", stringMcParticles + stringDeltaPhi + stringDeltaEta + stringPtLc + "entries", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtMcGenPrompt", stringLcHadron + " Prompt MC Gen", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtMcGenPromptLcPromptHadron", stringLcHadron + "prompt Lc prompt h MC Gen", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtMcGenNonPromptLcNonPromptHadron", stringLcHadron + " non prompt Lc non prompt h MC Gen", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtMcGenNonPrompt", stringLcHadron + " NonPrompt MC Gen", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisPoolBin}}});

        registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGenPrompt"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGenPromptLcPromptHadron"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGenNonPromptLcNonPromptHadron"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGenNonPrompt"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtMcGen"))->Sumw2();
      } else {
        registry.add("hCorrel2DVsPtSignMcGenPrompt", stringLcHadron + " Prompt MC Gen", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisSignPair}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtSignMcGenPromptLcPromptHadron", stringLcHadron + "prompt Lc prompt h MC Gen", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisSignPair}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtSignMcGenNonPromptLcNonPromptHadron", stringLcHadron + " non prompt Lc non prompt h MC Gen", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisSignPair}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtSignMcGenNonPrompt", stringLcHadron + " NonPrompt MC Gen", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtLc}, {axisPtHadron}, {axisSignPair}, {axisPoolBin}}});
        registry.add("hCorrel2DVsPtSignMcGen", stringMcParticles + stringDeltaPhi + stringDeltaEta + stringPtLc + stringSign + "entries", {HistType::kTHnSparseF, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtCorr}, {axisPtHadron}, {axisSignPair}, {axisPoolBin}}});

        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignMcGenPrompt"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignMcGenPromptLcPromptHadron"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignMcGenNonPromptLcNonPromptHadron"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignMcGenNonPrompt"))->Sumw2();
        registry.get<THnSparse>(HIST("hCorrel2DVsPtSignMcGen"))->Sumw2();
      }

      // Toward Transverse Away for McGen
      registry.add("hTowardGen", "Toward invmass; ptLc; correlationState;entries", {HistType::kTH3F, {{axisMassLc}, {axisPtLc}, {axisCorrelationState}}});
      registry.add("hTransverseGen", "Transverse invmass; ptLc; correlationState;entries", {HistType::kTH3F, {{axisMassLc}, {axisPtLc}, {axisCorrelationState}}});
      registry.add("hAwayGen", "Away invmass; ptLc; correlationState;entries", {HistType::kTH3F, {{axisMassLc}, {axisPtLc}, {axisCorrelationState}}});
    }
    // Histograms for efficiencies
    registry.add("Efficiency/hPtCandMcRecPrompt", stringMcRecoLcPrompt + stringPtLc, {HistType::kTH1D, {axisPtLc}});
    registry.add("Efficiency/hPtCandMcGenPrompt", stringMcGenLcPrompt + stringPtLc, {HistType::kTH1D, {axisPtLc}});
    registry.add("Efficiency/hPtCandMcRecNonPrompt", stringMcRecoLcFd + stringPtLc, {HistType::kTH1D, {axisPtLc}});
    registry.add("Efficiency/hPtCandMcGenNonPrompt", stringMcGenLcFd + stringPtLc, {HistType::kTH1D, {axisPtLc}});
    registry.add("Efficiency/hPtCandMcGenDaughterInAcc", stringMcGenLcFd + stringPtLc, {HistType::kTH1D, {axisPtLc}});

    auto hCandidates = registry.add<StepTHn>("hCandidates", "Candidate count at different steps", {HistType::kStepTHnF, {axisPtLc, axisMultFT0M, {RecoDecay::OriginType::NonPrompt + 1, +RecoDecay::OriginType::None - 0.5, +RecoDecay::OriginType::NonPrompt + 0.5}}, kCandidateNSteps});
    hCandidates->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hCandidates->GetAxis(1)->SetTitle("multiplicity");
    hCandidates->GetAxis(2)->SetTitle("Charm hadron origin");

    // Loading efficiency histograms from CCDB
    if ((applyEfficiency != 0) && loadAccXEffFromCCDB) {
      ccdb->setURL(ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);

      mEfficiencyPrompt = std::shared_ptr<TH1>(ccdb->getForTimeStamp<TH1D>(promptEffCcdbPath, timestampCcdb));
      if (mEfficiencyPrompt == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", promptEffCcdbPath.value.c_str());
      }
      LOGF(info, "Loaded trigger efficiency (prompt Lc) histogram from %s", promptEffCcdbPath.value.c_str());

      mEfficiencyFD = std::shared_ptr<TH1>(ccdb->getForTimeStamp<TH1D>(fdEffCcdbPath, timestampCcdb));
      if (mEfficiencyFD == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", fdEffCcdbPath.value.c_str());
      }
      LOGF(info, "Loaded feed-down Lc efficiency histogram from %s", fdEffCcdbPath.value.c_str());

      mEfficiencyAssociated = std::shared_ptr<TH1>(ccdb->getForTimeStamp<TH1D>(associatedEffCcdbPath, timestampCcdb));
      if (mEfficiencyAssociated == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for associated particles from %s", associatedEffCcdbPath.value.c_str());
      }
      LOGF(info, "Loaded associated efficiency histogram from %s", associatedEffCcdbPath.value.c_str());
    }

    if (activateQA) {
      const int regionLimits = 6;
      std::string labels[regionLimits] = {"SigReg Left", "SigReg Right", "Left SB Low", "Left SB Up", "Right SB Low", "Right SB Up"};
      static const AxisSpec axisSidebandLimits = {regionLimits, 0.5, 6.5, ""};
      auto hSigSidebandLimits = registry.add<TH2>("Inputs/hSigSidebandLimits", "Signal and Sideband Limits;;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {axisSidebandLimits, {(std::vector<double>)binsPtCorrelations, "#it{p}_{T} (GeV/#it{c})"}}});
      for (int iLim = 0; iLim < regionLimits; iLim++) {
        hSigSidebandLimits->GetXaxis()->SetBinLabel(iLim + 1, labels[iLim].data());
      }
      for (size_t iPtLc = 0; iPtLc < binsPtCorrelations->size() - 1; iPtLc++) {
        hSigSidebandLimits->SetBinContent(1, iPtLc + 1, signalRegionInner->at(iPtLc));
        hSigSidebandLimits->SetBinContent(2, iPtLc + 1, signalRegionOuter->at(iPtLc));
        hSigSidebandLimits->SetBinContent(3, iPtLc + 1, sidebandLeftOuter->at(iPtLc));
        hSigSidebandLimits->SetBinContent(4, iPtLc + 1, sidebandLeftInner->at(iPtLc));
        hSigSidebandLimits->SetBinContent(5, iPtLc + 1, sidebandRightInner->at(iPtLc));
        hSigSidebandLimits->SetBinContent(6, iPtLc + 1, sidebandRightOuter->at(iPtLc));
      }
    }
  }

  void processData(LcHadronPairFullWithMl const& pairEntries, aod::LcRecoInfo const& candidates)
  {
    for (const auto& candidate : candidates) {
      float const massLc = candidate.mLc();
      float const ptLc = std::abs(candidate.ptLc());
      float const bdtScorePrompt = candidate.mlScorePrompt();
      float const bdtScoreBkg = candidate.mlScoreBkg();
      int const effBinLc = o2::analysis::findBin(binsPtEfficiencyLc, ptLc);

      // reject entries outside Pt ranges of interest
      if (ptLc < binsPtEfficiencyLc->front() || ptLc > binsPtEfficiencyLc->back()) {
        continue;
      }

      if (bdtScorePrompt < mlOutputPrompt->at(effBinLc) || bdtScoreBkg > mlOutputBkg->at(effBinLc)) {
        continue;
      }
      double efficiencyWeightLc = 1.;
      if (applyEfficiency != 0) {
        efficiencyWeightLc = 1. / efficiencyLc->at(o2::analysis::findBin(binsPtEfficiencyLc, ptLc));
        if (loadAccXEffFromCCDB) {
          efficiencyWeightLc = 1. / mEfficiencyPrompt->GetBinContent(mEfficiencyPrompt->FindBin(ptLc));
        }
      }
      registry.fill(HIST("hMassLcVsPt"), massLc, ptLc, efficiencyWeightLc);
      registry.fill(HIST("hMassLcVsPtWoEff"), massLc, ptLc);
      registry.fill(HIST("hBdtScorePrompt"), bdtScorePrompt);
      registry.fill(HIST("hBdtScoreBkg"), bdtScoreBkg);
    }

    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float const deltaPhi = pairEntry.deltaPhi();
      float cent = 0.;
      if (useCentrality) {
        cent = pairEntry.cent();
      }
      float const deltaEta = pairEntry.deltaEta();
      double const ptLc = std::abs(pairEntry.ptLc());
      double const ptHadron = std::abs(pairEntry.ptHadron());
      float const bdtScorePrompt = pairEntry.mlScorePrompt();
      float const bdtScoreBkg = pairEntry.mlScoreBkg();
      float const trackDcaXY = pairEntry.trackDcaXY();
      float const trackDcaZ = pairEntry.trackDcaZ();
      int const trackTpcCrossedRows = pairEntry.trackTPCNClsCrossedRows();
      int const poolBin = pairEntry.poolBin();
      double const massLc = pairEntry.mLc();
      int const effBinLc = o2::analysis::findBin(binsPtEfficiencyLc, ptLc);
      int const ptBinLc = o2::analysis::findBin(binsPtCorrelations, ptLc);
      bool const isAutoCorrelated = pairEntry.isAutoCorrelated();
      int signPair = 0;
      // reject entries outside Pt ranges of interest
      if (ptBinLc < 0 || effBinLc < 0) {
        continue;
      }

      if (bdtScorePrompt < mlOutputPrompt->at(effBinLc) || bdtScoreBkg > mlOutputBkg->at(effBinLc)) {
        continue;
      }
      if (trackDcaXY > dcaXYTrackMax || trackDcaZ > dcaZTrackMax || trackTpcCrossedRows < nTpcCrossedRaws) {
        continue;
      }

      double efficiencyWeight = 1.;
      if (applyEfficiency != 0) {
        efficiencyWeight = 1. / (efficiencyLc->at(effBinLc) * efficiencyHad->at(o2::analysis::findBin(binsPtEfficiencyHad, ptHadron)));
        if (loadAccXEffFromCCDB) {
          efficiencyWeight = 1. / (mEfficiencyPrompt->GetBinContent(mEfficiencyPrompt->FindBin(ptLc)) * mEfficiencyAssociated->GetBinContent(mEfficiencyAssociated->FindBin(ptHadron)));
        }
      }

      // Divide into three regions: toward, transverse, and away
      if (isTowardTransverseAway) {
        if (ptHadron < leadingParticlePtMin) {
          continue;
        }
        Region const region = getRegion(deltaPhi);
        switch (region) {
          case Toward:
            registry.fill(HIST("hToward"), massLc, ptLc, isAutoCorrelated, efficiencyWeight);
            break;
          case Away:
            registry.fill(HIST("hAway"), massLc, ptLc, isAutoCorrelated, efficiencyWeight);
            break;
          case Transverse:
            registry.fill(HIST("hTransverse"), massLc, ptLc, isAutoCorrelated, efficiencyWeight);
            break;
          default:
            break;
        }
      }

      if (calSign) {
        signPair = signCalulation(pairEntry.ptLc(), pairEntry.ptHadron());
      }
      if (fillSignAndMass) {
        registry.fill(HIST("hCorrel2DVsPtSignMass"), deltaPhi, deltaEta, ptLc, ptHadron, massLc, signPair, poolBin, efficiencyWeight);
      }
      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if (storeMass) {
        registry.fill(HIST("hCorrel2DVsPtGlobalRegion"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, cent, massLc, efficiencyWeight);
        continue;
      }
      if (massLc > signalRegionInner->at(ptBinLc) && massLc < signalRegionOuter->at(ptBinLc)) {
        // in signal region
        if (fillSign) {
          registry.fill(HIST("hCorrel2DVsPtSignSignalRegion"), deltaPhi, deltaEta, ptLc, ptHadron, signPair, poolBin, efficiencyWeight);
        } else {
          registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, cent, efficiencyWeight);
        }
        registry.fill(HIST("hCorrel2DPtIntSignalRegion"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
      }
      // in sideband left region
      if (massLc > sidebandLeftOuter->at(ptBinLc) && massLc < sidebandLeftInner->at(ptBinLc)) {
        if (fillSign) {
          registry.fill(HIST("hCorrel2DVsPtSignSidebandLeft"), deltaPhi, deltaEta, ptLc, ptHadron, signPair, poolBin, efficiencyWeight);
        } else {
          registry.fill(HIST("hCorrel2DVsPtSidebandLeft"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, cent, efficiencyWeight);
        }
        registry.fill(HIST("hDeltaEtaPtIntSidebandLeft"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandLeft"), deltaPhi, efficiencyWeight);
        registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, cent, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebands"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi, efficiencyWeight);
      }
      // in sideband right region
      if (massLc > sidebandRightInner->at(ptBinLc) && massLc < sidebandRightOuter->at(ptBinLc)) {
        if (fillSign) {
          registry.fill(HIST("hCorrel2DVsPtSignSidebandRight"), deltaPhi, deltaEta, ptLc, ptHadron, signPair, poolBin, efficiencyWeight);
        } else {
          registry.fill(HIST("hCorrel2DVsPtSidebandRight"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, cent, efficiencyWeight);
        }
        registry.fill(HIST("hDeltaEtaPtIntSidebandRight"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandRight"), deltaPhi, efficiencyWeight);
        registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, cent, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebands"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi, efficiencyWeight);
      }
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationLcHadrons, processData, "Process data", true);

  /// Lc-Hadron correlation pair filling task, from pair tables - for Mc reco-level analysis (candidates matched to true signal only, but also bkg sources are studied)
  void processMcRec(LcHadronPairFullWithMl const& pairEntries,
                    soa::Join<aod::LcRecoInfo, aod::LcGenInfo> const& candidates)
  {
    for (const auto& candidate : candidates) {
      float const massLc = candidate.mLc();
      float const ptLc = std::abs(candidate.ptLc());
      float const bdtScorePrompt = candidate.mlScorePrompt();
      float const bdtScoreBkg = candidate.mlScoreBkg();
      int const effBinLc = o2::analysis::findBin(binsPtEfficiencyLc, ptLc);
      bool const isLcPrompt = candidate.isPrompt();

      // reject entries outside pT ranges of interest
      if (ptLc < binsPtEfficiencyLc->front() || ptLc > binsPtEfficiencyLc->back()) {
        continue;
      }

      if (bdtScorePrompt < mlOutputPrompt->at(effBinLc) || bdtScoreBkg > mlOutputBkg->at(effBinLc)) {
        continue;
      }
      double efficiencyWeightLc = 1.;
      if (applyEfficiency != 0) {
        if (isLcPrompt) {
          efficiencyWeightLc = 1. / efficiencyLc->at(effBinLc);
          if (loadAccXEffFromCCDB) {
            efficiencyWeightLc = 1. / mEfficiencyPrompt->GetBinContent(mEfficiencyPrompt->FindBin(ptLc));
          }
          registry.fill(HIST("hMassLcVsPt"), massLc, ptLc, efficiencyWeightLc);
          registry.fill(HIST("hMassLcVsPtWoEff"), massLc, ptLc);
          registry.fill(HIST("hMassPromptLcVsPt"), massLc, ptLc, efficiencyWeightLc);
          registry.fill(HIST("hBdtScorePrompt"), bdtScorePrompt);
          registry.fill(HIST("hBdtScoreBkg"), bdtScoreBkg);
        } else {
          efficiencyWeightLc = 1. / efficiencyFdLc->at(effBinLc);
          if (loadAccXEffFromCCDB) {
            efficiencyWeightLc = 1. / mEfficiencyFD->GetBinContent(mEfficiencyFD->FindBin(ptLc));
          }
          registry.fill(HIST("hMassLcVsPt"), massLc, ptLc, efficiencyWeightLc);
          registry.fill(HIST("hMassLcVsPtWoEff"), massLc, ptLc);
          registry.fill(HIST("hMassNonPromptLcVsPt"), massLc, ptLc, efficiencyWeightLc);
          registry.fill(HIST("hBdtScorePrompt"), bdtScorePrompt);
          registry.fill(HIST("hBdtScoreBkg"), bdtScoreBkg);
        }
      }
    }

    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float const deltaPhi = pairEntry.deltaPhi();
      float const deltaEta = pairEntry.deltaEta();
      float const ptLc = std::abs(pairEntry.ptLc());
      float const ptHadron = std::abs(pairEntry.ptHadron());
      float const massLc = pairEntry.mLc();
      float const bdtScorePrompt = pairEntry.mlScorePrompt();
      float const bdtScoreBkg = pairEntry.mlScoreBkg();
      bool const isPhysicalPrimary = pairEntry.isPhysicalPrimary();
      float const trackDcaXY = pairEntry.trackDcaXY();
      float const trackDcaZ = pairEntry.trackDcaZ();
      int const trackTpcCrossedRows = pairEntry.trackTPCNClsCrossedRows();
      int const statusLcPrompt = static_cast<int>(pairEntry.isPrompt());
      int const statusPromptHadron = pairEntry.trackOrigin();
      int const poolBin = pairEntry.poolBin();
      int const effBinLc = o2::analysis::findBin(binsPtEfficiencyLc, ptLc);
      int const ptBinLc = o2::analysis::findBin(binsPtCorrelations, ptLc);
      bool const isAutoCorrelated = pairEntry.isAutoCorrelated();
      int signPair = 0;

      // reject entries outside pT ranges of interest
      if (ptLc < binsPtEfficiencyLc->front() || ptLc > binsPtEfficiencyLc->back()) {
        continue;
      }

      if (bdtScorePrompt < mlOutputPrompt->at(effBinLc) || bdtScoreBkg > mlOutputBkg->at(effBinLc)) {
        continue;
      }
      if (trackDcaXY > dcaXYTrackMax || trackDcaZ > dcaZTrackMax || trackTpcCrossedRows < nTpcCrossedRaws) {
        continue;
      }
      double efficiencyWeight = 1.;

      if (applyEfficiency != 0) {
        if (statusLcPrompt != 0) {
          efficiencyWeight = 1. / (efficiencyLc->at(effBinLc) * efficiencyHad->at(o2::analysis::findBin(binsPtEfficiencyHad, ptHadron)));
          if (loadAccXEffFromCCDB) {
            efficiencyWeight = 1. / (mEfficiencyPrompt->GetBinContent(mEfficiencyPrompt->FindBin(ptLc)) * mEfficiencyAssociated->GetBinContent(mEfficiencyAssociated->FindBin(ptHadron)));
          }
        } else {
          efficiencyWeight = 1. / (efficiencyFdLc->at(effBinLc) * efficiencyHad->at(o2::analysis::findBin(binsPtEfficiencyHad, ptHadron)));
          if (loadAccXEffFromCCDB) {
            efficiencyWeight = 1. / (mEfficiencyFD->GetBinContent(mEfficiencyFD->FindBin(ptLc)) * mEfficiencyAssociated->GetBinContent(mEfficiencyAssociated->FindBin(ptHadron)));
          }
        }
      }

      // Divide into three regions: toward, transverse, and away
      if (isTowardTransverseAway) {
        if (ptHadron < leadingParticlePtMin) {
          continue;
        }
        Region const region = getRegion(deltaPhi);
        switch (region) {
          case Toward:
            registry.fill(HIST("hTowardRec"), massLc, ptLc, isAutoCorrelated, efficiencyWeight);
            break;
          case Away:
            registry.fill(HIST("hAwayRec"), massLc, ptLc, isAutoCorrelated, efficiencyWeight);
            break;
          case Transverse:
            registry.fill(HIST("hTransverseRec"), massLc, ptLc, isAutoCorrelated, efficiencyWeight);
            break;
          default:
            break;
        }
      }

      if (calSign) {
        signPair = signCalulation(pairEntry.ptLc(), pairEntry.ptHadron());
      }
      if (fillSignAndMass) {
        registry.fill(HIST("hCorrel2DVsPtSignMassMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, massLc, signPair, poolBin, efficiencyWeight);
      }

      // fill correlation plots for signal/bagkground correlations
      if (pairEntry.signalStatus() != 0) {
        if (fillSign) {
          registry.fill(HIST("hCorrel2DVsPtSignSignalMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, signPair, poolBin, efficiencyWeight);
        } else {
          registry.fill(HIST("hCorrel2DVsPtSignalMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, efficiencyWeight);
        }
      } else {
        if (fillSign) {
          registry.fill(HIST("hCorrel2DVsPtSignBkgMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, signPair, poolBin, efficiencyWeight);
        } else {
          registry.fill(HIST("hCorrel2DVsPtBkgMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, efficiencyWeight);
        }
      }
      // reject entries outside Pt ranges of interest

      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if (massLc > signalRegionInner->at(ptBinLc) && massLc < signalRegionOuter->at(ptBinLc)) {
        // in signal region
        if (fillSign) {
          registry.fill(HIST("hCorrel2DVsPtSignSignalRegionMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, signPair, poolBin, efficiencyWeight);
        } else {
          registry.fill(HIST("hCorrel2DVsPtSignalRegionMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, statusLcPrompt, poolBin, efficiencyWeight);
        }
        registry.fill(HIST("hCorrel2DPtIntSignalRegionMcRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionMcRec"), deltaPhi, efficiencyWeight);
        if (isPhysicalPrimary) {
          registry.fill(HIST("hCorrel2DVsPtPhysicalPrimaryMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, statusLcPrompt, poolBin, efficiencyWeight);
          if (statusLcPrompt == 1 && statusPromptHadron == RecoDecay::OriginType::Prompt) {
            if (fillSign) {
              registry.fill(HIST("hCorrel2DVsPtSignSignalRegionPromptLcPromptHadronMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, signPair, poolBin, efficiencyWeight);
            } else {
              registry.fill(HIST("hCorrel2DVsPtSignalRegionPromptLcPromptHadronMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, efficiencyWeight);
            }
          } else if (statusLcPrompt == 0 && statusPromptHadron == RecoDecay::OriginType::NonPrompt) {
            if (fillSign) {
              registry.fill(HIST("hCorrel2DVsPtSignSignalRegionNonPromptLcNonPromptHadronMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, signPair, poolBin, efficiencyWeight);
            } else {
              registry.fill(HIST("hCorrel2DVsPtSignalRegionNonPromptLcNonPromptHadronMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, efficiencyWeight);
            }
          }
        }
      }
      // in sideband left region
      if (massLc > sidebandLeftOuter->at(ptBinLc) && massLc < sidebandLeftInner->at(ptBinLc)) {
        if (fillSign) {
          registry.fill(HIST("hCorrel2DVsPtSignSidebandLeftMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, signPair, poolBin, efficiencyWeight);
        } else {
          registry.fill(HIST("hCorrel2DVsPtSidebandLeftMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, efficiencyWeight);
        }
        registry.fill(HIST("hDeltaEtaPtIntSidebandLeftMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandLeftMcRec"), deltaPhi, efficiencyWeight);
        registry.fill(HIST("hCorrel2DVsPtSidebandsMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsMcRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsMcRec"), deltaPhi, efficiencyWeight);
      }
      // in sideband right region
      if (massLc > sidebandRightInner->at(ptBinLc) && massLc < sidebandRightOuter->at(ptBinLc)) {
        if (fillSign) {
          registry.fill(HIST("hCorrel2DVsPtSignSidebandRightMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, signPair, poolBin, efficiencyWeight);
        } else {
          registry.fill(HIST("hCorrel2DVsPtSidebandRightMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, efficiencyWeight);
        }
        registry.fill(HIST("hDeltaEtaPtIntSidebandRightMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandRightMcRec"), deltaPhi, efficiencyWeight);
        registry.fill(HIST("hCorrel2DVsPtSidebandsMcRec"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsMcRec"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsMcRec"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsMcRec"), deltaPhi, efficiencyWeight);
      }
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationLcHadrons, processMcRec, "Process Mc Reco mode", false);

  /// Lc-Hadron correlation pair filling task, from pair tables - for Mc gen-level analysis (no filter/selection, only true signal)
  void processMcGen(LcHadronPair const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      float const deltaPhi = pairEntry.deltaPhi();
      float const deltaEta = pairEntry.deltaEta();
      float const ptLc = std::abs(pairEntry.ptLc());
      float const ptHadron = std::abs(pairEntry.ptHadron());
      int const poolBin = pairEntry.poolBin();
      int const statusPromptHadron = pairEntry.trackOrigin();
      bool const isLcPrompt = pairEntry.isPrompt();
      bool const isAutoCorrelated = pairEntry.isAutoCorrelated();
      int signPair = 0;

      if (isTowardTransverseAway) {
        // Divide into three regions: toward, transverse, and away
        if (ptHadron < leadingParticlePtMin) {
          continue;
        }
        Region const region = getRegion(deltaPhi);
        switch (region) {
          case Toward:
            registry.fill(HIST("hTowardRec"), o2::constants::physics::MassLambdaCPlus, ptLc, isAutoCorrelated);
            break;
          case Away:
            registry.fill(HIST("hAwayRec"), o2::constants::physics::MassLambdaCPlus, ptLc, isAutoCorrelated);
            break;
          case Transverse:
            registry.fill(HIST("hTransverseRec"), o2::constants::physics::MassLambdaCPlus, ptLc, isAutoCorrelated);
            break;
          default:
            break;
        }
      }
      if (calSign) {
        signPair = signCalulation(pairEntry.ptLc(), pairEntry.ptHadron());
      }
      if (fillSign) {
        registry.fill(HIST("hCorrel2DVsPtSignMcGen"), deltaPhi, deltaEta, ptLc, ptHadron, signPair, poolBin);
      } else {
        registry.fill(HIST("hCorrel2DVsPtMcGen"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin);
      }
      registry.fill(HIST("hCorrel2DPtIntMcGen"), deltaPhi, deltaEta);
      registry.fill(HIST("hDeltaEtaPtIntMcGen"), deltaEta);
      registry.fill(HIST("hDeltaPhiPtIntMcGen"), deltaPhi);
      if (isLcPrompt) {
        if (fillSign) {
          registry.fill(HIST("hCorrel2DVsPtSignMcGenPrompt"), deltaPhi, deltaEta, ptLc, ptHadron, signPair, poolBin);
        } else {
          registry.fill(HIST("hCorrel2DVsPtMcGenPrompt"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin);
        }
        if (statusPromptHadron == RecoDecay::OriginType::Prompt) {
          if (fillSign) {
            registry.fill(HIST("hCorrel2DVsPtSignMcGenPromptLcPromptHadron"), deltaPhi, deltaEta, ptLc, ptHadron, signPair, poolBin);
          } else {
            registry.fill(HIST("hCorrel2DVsPtMcGenPromptLcPromptHadron"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin);
          }
        }
      } else {
        if (fillSign) {
          registry.fill(HIST("hCorrel2DVsPtSignMcGenNonPrompt"), deltaPhi, deltaEta, ptLc, ptHadron, signPair, poolBin);
        } else {
          registry.fill(HIST("hCorrel2DVsPtMcGenNonPrompt"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin);
        }
        if (statusPromptHadron == RecoDecay::OriginType::NonPrompt) {
          if (fillSign) {
            registry.fill(HIST("hCorrel2DVsPtSignMcGenNonPromptLcNonPromptHadron"), deltaPhi, deltaEta, ptLc, ptHadron, signPair, poolBin);
          } else {
            registry.fill(HIST("hCorrel2DVsPtMcGenNonPromptLcNonPromptHadron"), deltaPhi, deltaEta, ptLc, ptHadron, poolBin);
          }
        }
      }
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationLcHadrons, processMcGen, "Process Mc Gen mode", false);

  /// Lc-Hadron correlation - reconstruction and selection efficiency
  void processMcCandEfficiency(soa::Join<aod::Collisions, aod::FT0Mults, aod::EvSels> const&,
                               soa::Join<aod::McCollisions, aod::MultsExtraMC> const&,
                               CandLcMcGen const& mcParticles,
                               CandLcMcReco const& candidates,
                               aod::TracksWMc const&)
  {
    auto hCandidates = registry.get<StepTHn>(HIST("hCandidates"));

    /// Gen loop
    float multiplicity = -1.;
    for (const auto& mcParticle : mcParticles) {
      // generated candidates
      if (std::abs(mcParticle.pdgCode()) == Pdg::kLambdaCPlus) {
        auto mcCollision = mcParticle.template mcCollision_as<soa::Join<aod::McCollisions, aod::MultsExtraMC>>();
        multiplicity = mcCollision.multMCFT0A() + mcCollision.multMCFT0C(); // multFT0M = multFt0A + multFT0C
        hCandidates->Fill(kCandidateStepMcGenAll, mcParticle.pt(), multiplicity, mcParticle.originMcGen());
        if (std::abs(mcParticle.flagMcMatchGen()) == hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi) {
          hCandidates->Fill(kCandidateStepMcGenLcToPKPi, mcParticle.pt(), multiplicity, mcParticle.originMcGen());
          auto yL = RecoDecay::y(mcParticle.pVector(), o2::constants::physics::MassLambdaCPlus);
          if (std::abs(yL) <= yCandGenMax) {
            hCandidates->Fill(kCandidateStepMcCandInAcceptance, mcParticle.pt(), multiplicity, mcParticle.originMcGen());
            if (mcParticle.originMcGen() == RecoDecay::OriginType::Prompt) {
              registry.fill(HIST("Efficiency/hPtCandMcGenPrompt"), mcParticle.pt());
            }
            if (mcParticle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
              registry.fill(HIST("Efficiency/hPtCandMcGenNonPrompt"), mcParticle.pt());
            }
          }
          bool isDaughterInAcceptance = true;
          auto daughters = mcParticle.template daughters_as<CandLcMcGen>();
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
      if (candidate.isSelLcToPKPi() >= selectionFlagLc) {
        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
          outputMl[iclass] = candidate.mlProbLcToPKPi()[classMl->at(iclass)];
        }
      }
      if (candidate.isSelLcToPiKP() >= selectionFlagLc) {
        for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
          outputMl[iclass] = candidate.mlProbLcToPiKP()[classMl->at(iclass)];
        }
      }
      if (outputMl[0] > mlOutputBkg->at(o2::analysis::findBin(binsPtEfficiencyLc, candidate.pt())) || outputMl[1] < mlOutputPrompt->at(o2::analysis::findBin(binsPtEfficiencyLc, candidate.pt()))) {
        continue;
      }
      auto collision = candidate.template collision_as<soa::Join<aod::Collisions, aod::FT0Mults, aod::EvSels>>();
      if (selNoSameBunchPileUpColl && !(collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))) {
        continue;
      }
      multiplicity = collision.multFT0M();
      if (std::abs(candidate.flagMcMatchRec()) == hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi) {
        hCandidates->Fill(kCandidateStepMcReco, candidate.pt(), multiplicity, candidate.originMcRec());
        if (std::abs(HfHelper::yLc(candidate)) <= yCandMax) {
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
  PROCESS_SWITCH(HfTaskCorrelationLcHadrons, processMcCandEfficiency, "Process MC for calculating candidate reconstruction efficiency", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationLcHadrons>(cfgc)};
}
