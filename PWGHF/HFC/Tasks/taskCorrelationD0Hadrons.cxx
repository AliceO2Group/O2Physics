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

/// \file taskCorrelationD0Hadrons.cxx
/// \brief D0-Hadron correlator task - data-like, MC-reco and MC-kine analyses.
/// \note Extended from taskCorrelationDDbar
///
/// \author Samrangy Sadhu <samrangy.sadhu@cern.ch>, INFN Bari
/// \author Swapnesh Santosh Khade <swapnesh.santosh.khade@cern.ch>, IIT Indore

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/CorrelationTables.h"
#include "PWGHF/HFC/Utils/utilsCorrelations.h"
#include "PWGHF/Utils/utilsAnalysis.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/StepTHn.h>
#include <Framework/runDataProcessing.h>

#include <THnSparse.h>

#include <Rtypes.h>

#include <cstdlib>
#include <vector>

using namespace o2;
using namespace o2::constants::physics;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_correlation_d0_hadron;
using namespace o2::analysis::hf_correlations;

// definition of vectors for standard ptbin and invariant mass configurables
const int nPtBinsCorrelations = 12;
const double pTBinsCorrelations[nPtBinsCorrelations + 1] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 12., 16., 24., 99.};
const auto vecPtBinsCorrelations = std::vector<double>{pTBinsCorrelations, pTBinsCorrelations + nPtBinsCorrelations + 1};
const double signalRegionLeftDefault[nPtBinsCorrelations] = {1.7948, 1.8198, 1.8198, 1.8148, 1.8148, 1.8048, 1.8048, 1.7948, 1.7948, 1.7898, 1.7848, 1.7598};
const double signalRegionRightDefault[nPtBinsCorrelations] = {1.9098, 1.8998, 1.9048, 1.9048, 1.9148, 1.9248, 1.9298, 1.9348, 1.9398, 1.9298, 1.9398, 1.9198};
const double sidebandLeftInnerDefault[nPtBinsCorrelations] = {1.7398, 1.7748, 1.7798, 1.7698, 1.7648, 1.7448, 1.7448, 1.7198, 1.7198, 1.7198, 1.7048, 1.6798};
const double sidebandLeftOuterDefault[nPtBinsCorrelations] = {1.6298, 1.6898, 1.6948, 1.6748, 1.6648, 1.6248, 1.6198, 1.5748, 1.5748, 1.5798, 1.5448, 1.5198};
const double sidebandRightInnerDefault[nPtBinsCorrelations] = {1.9648, 1.9448, 1.9448, 1.9548, 1.9648, 1.9848, 1.9948, 2.0098, 2.0148, 1.9998, 2.0248, 1.9998};
const double sidebandRightOuterDefault[nPtBinsCorrelations] = {2.0748, 2.0248, 2.0298, 2.0448, 2.0648, 2.1048, 2.1148, 2.1548, 2.1648, 2.1398, 2.1848, 2.1598};
const auto vecsignalRegionLeft = std::vector<double>{signalRegionLeftDefault, signalRegionLeftDefault + nPtBinsCorrelations};
const auto vecsignalRegionRight = std::vector<double>{signalRegionRightDefault, signalRegionRightDefault + nPtBinsCorrelations};
const auto vecSidebandLeftInner = std::vector<double>{sidebandLeftInnerDefault, sidebandLeftInnerDefault + nPtBinsCorrelations};
const auto vecSidebandLeftOuter = std::vector<double>{sidebandLeftOuterDefault, sidebandLeftOuterDefault + nPtBinsCorrelations};
const auto vecSidebandRightInner = std::vector<double>{sidebandRightInnerDefault, sidebandRightInnerDefault + nPtBinsCorrelations};
const auto vecSidebandRightOuter = std::vector<double>{sidebandRightOuterDefault, sidebandRightOuterDefault + nPtBinsCorrelations};
const int nPtbinsPtEfficiencyD = o2::analysis::hf_cuts_d0_to_pi_k::NBinsPt;
const double efficiencyDmesonDefault[nPtbinsPtEfficiencyD] = {};
const auto vecEfficiencyDmeson = std::vector<double>{efficiencyDmesonDefault, efficiencyDmesonDefault + nPtbinsPtEfficiencyD};

struct HfTaskCorrelationD0Hadrons {

  // pT ranges: the default values are those embedded in hf_cuts_d0_to_pi_k (i.e. the mass pT bins), but can be redefined via json files
  Configurable<std::vector<double>> binsPtD{"binsPtD", std::vector<double>{o2::analysis::hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits for candidate mass plots and efficiency"};
  Configurable<std::vector<double>> binsPtHadron{"binsPtHadron", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "pT bin limits for assoc particle efficiency"};
  Configurable<std::vector<double>> binsCorrelations{"binsCorrelations", std::vector<double>{vecPtBinsCorrelations}, "pT bin limits for correlation plots"};
  Configurable<std::vector<double>> binsPtEfficiencyD{"binsPtEfficiencyD", std::vector<double>{o2::analysis::hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits for efficiency"};
  Configurable<std::vector<double>> binsPtEfficiencyHad{"binsPtEfficiencyHad", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "pT bin limits for associated particle efficiency"};
  Configurable<std::vector<double>> efficiencyDmeson{"efficiencyDmeson", std::vector<double>{vecEfficiencyDmeson}, "Efficiency values for D meson specie under study"};
  Configurable<std::vector<double>> efficiencyHad{"efficiencyHad", {1., 1., 1., 1., 1., 1.}, "efficiency values for associated particles"};
  Configurable<int> applyEfficiency{"applyEfficiency", 1, "Flag for applying efficiency weights"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 1, 2}, "Indexes of ML scores to be stored. Three indexes max."};
  Configurable<std::vector<double>> mlOutputPromptD0{"mlOutputPromptD0", {0.5, 0.5, 0.5, 0.5}, "Machine learning scores for prompt"};
  Configurable<std::vector<double>> mlOutputBkgD0{"mlOutputBkgD0", {0.5, 0.5, 0.5, 0.5}, "Machine learning scores for bkg"};
  Configurable<std::vector<double>> mlOutputPromptD0bar{"mlOutputPromptD0bar", {0.5, 0.5, 0.5, 0.5}, "Machine learning scores for prompt"};
  Configurable<std::vector<double>> mlOutputBkgD0bar{"mlOutputBkgD0bar", {0.5, 0.5, 0.5, 0.5}, "Machine learning scores for bkg"};
  Configurable<float> ptCandMin{"ptCandMin", 1., "min. cand. pT"};
  Configurable<float> ptCandMax{"ptCandMax", 50., "max. cand pT"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.3, "min. track pT"};
  Configurable<float> ptTrackMax{"ptTrackMax", 99., "max. track pT"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. eta of tracks"};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> yCandGenMax{"yCandGenMax", 0.5, "max. gen. cand. rapidity"};
  Configurable<float> ptDaughterMin{"ptDaughterMin", 0.1, "min. daughter pT"};
  Configurable<int> nTpcCrossedRaws{"nTpcCrossedRaws", 70, "Number of crossed TPC Rows"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 1., "max. DCA_xy of tracks"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1., "max. DCA_z of tracks"};
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0 (bar)"};
  Configurable<bool> selNoSameBunchPileUpColl{"selNoSameBunchPileUpColl", true, "Flag for rejecting the collisions associated with the same bunch crossing"};
  Configurable<std::vector<double>> signalRegionLeft{"signalRegionLeft", std::vector<double>{vecsignalRegionLeft}, "Inner values of signal region vs pT"};
  Configurable<std::vector<double>> signalRegionRight{"signalRegionRight", std::vector<double>{vecsignalRegionRight}, "Outer values of signal region vs pT"};
  Configurable<std::vector<double>> sidebandLeftInner{"sidebandLeftInner", std::vector<double>{vecSidebandLeftInner}, "Inner values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandLeftOuter{"sidebandLeftOuter", std::vector<double>{vecSidebandLeftOuter}, "Outer values of left sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightInner{"sidebandRightInner", std::vector<double>{vecSidebandRightInner}, "Inner values of right sideband vs pT"};
  Configurable<std::vector<double>> sidebandRightOuter{"sidebandRightOuter", std::vector<double>{vecSidebandRightOuter}, "Outer values of right sideband vs pT"};
  Configurable<bool> isTowardTransverseAway{"isTowardTransverseAway", false, "Divide into three regions: toward, transverse, and away"};
  Configurable<double> leadingParticlePtMin{"leadingParticlePtMin", 0., "Min for leading particle pt"};

  HfHelper hfHelper;

  enum CandidateStep { kCandidateStepMcGenAll = 0,
                       kCandidateStepMcGenD0ToPiKPi,
                       kCandidateStepMcCandInAcceptance,
                       kCandidateStepMcDaughtersInAcceptance,
                       kCandidateStepMcReco,
                       kCandidateStepMcRecoInAcceptance,
                       kCandidateNSteps };

  using D0HadronPairFull = soa::Join<aod::D0HadronPair, aod::D0HadronRecoInfo, aod::D0HadronGenInfo>;
  using D0HadronPairFullMl = soa::Join<aod::D0HadronPair, aod::D0HadronRecoInfo, aod::D0HadronGenInfo, aod::D0HadronMlInfo, aod::D0TrackRecoInfo>;
  using CandD0McReco = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0, aod::HfCand2ProngMcRec>>;
  using CandD0McGen = soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>;
  using TracksWithMc = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra, o2::aod::McTrackLabels>>; // trackFilter applied

  Filter d0Filter = (aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0) || (aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0);
  Filter trackFilter = (nabs(aod::track::eta) < etaTrackMax) && (aod::track::pt > ptTrackMin) && (aod::track::pt < ptTrackMax) && (nabs(aod::track::dcaXY) < dcaXYTrackMax) && (nabs(aod::track::dcaZ) < dcaZTrackMax);

  ConfigurableAxis binsMassD{"binsMassD", {200, 1.3848, 2.3848}, "inv. mass (#pi K) (GeV/#it{c}^{2});entries"};
  ConfigurableAxis binsBdtScore{"binsBdtScore", {100, 0., 1.}, "Bdt output scores"};
  ConfigurableAxis binsEta{"binsEta", {100, -2., 2.}, "#it{#eta}"};
  ConfigurableAxis binsPhi{"binsPhi", {64, -PIHalf, 3. * PIHalf}, "#it{#varphi}"};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {600, 0., 8000.}, "Multiplicity as FT0M signal amplitude"};
  ConfigurableAxis binsPosZ{"binsPosZ", {100, -10., 10.}, "primary vertex z coordinate"};
  ConfigurableAxis binsPoolBin{"binsPoolBin", {9, 0., 9.}, "PoolBin"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    AxisSpec axisMassD = {binsMassD, "inv. mass (#pi K) (GeV/#it{c}^{2})"};
    AxisSpec axisDeltaEta = {binsEta, "#it{#eta}^{Hadron}-#it{#eta}^{D}"};
    AxisSpec const axisEta = {binsEta, "#it{#eta}"};
    AxisSpec axisDeltaPhi = {binsPhi, "#it{#varphi}^{Hadron}-#it{#varphi}^{D} (rad)"};
    AxisSpec axisPtD = {(std::vector<double>)binsPtD, "#it{p}_{T}^{D} (GeV/#it{c})"};
    AxisSpec axisPtHadron = {(std::vector<double>)binsPtHadron, "#it{p}_{T}^{Hadron} (GeV/#it{c})"};
    AxisSpec axisPoolBin = {binsPoolBin, "poolBin"};
    AxisSpec const axisBdtScore = {binsBdtScore, "Bdt score"};
    AxisSpec const axisMultFT0M = {binsMultFT0M, "MultiplicityFT0M"};
    AxisSpec const axisPosZ = {binsPosZ, "PosZ"};
    AxisSpec axisD0Prompt = {2, -0.5, 1.5, "Prompt D0"};
    AxisSpec axisCorrelationState = {2, 0., 2., "correlationState"};

    // Histograms for data
    registry.add("hDeltaEtaPtIntSignalRegion", "D0-h deltaEta signal region", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntSignalRegion", "D0-h deltaPhi signal region", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntSignalRegion", "D0-h deltaPhi vs deltaEta signal region", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtSignalRegion", "D0-h correlations signal region", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hDeltaEtaPtIntSignalRegionSoftPi", "D0-h deltaEta signal region soft pi only", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntSignalRegionSoftPi", "D0-h deltaPhi signal region soft pi only", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntSignalRegionSoftPi", "D0-h deltaPhi vs deltaEta signal region soft pi only", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtSignalRegionSoftPi", "D0-h correlations signal region soft pi only", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hDeltaEtaPtIntSidebands", "D0-h deltaEta sidebands", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntSidebands", "D0-h deltaPhi sidebands", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntSidebands", "D0-h deltaPhi vs deltaEta sidebands", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtSidebands", "D0-h correlations sidebands", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hDeltaEtaPtIntSidebandsSoftPi", "D0-h deltaEta sidebands soft pi only", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntSidebandsSoftPi", "D0-h deltaPhi sidebands soft pi only", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntSidebandsSoftPi", "D0-h deltaPhi vs deltaEta sidebands soft pi only", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtSidebandsSoftPi", "D0-h correlations sidebands soft pi only", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegion"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebands"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionSoftPi"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsSoftPi"))->Sumw2();
    // Toward Transverse Away for Data
    registry.add("hToward", "Toward", {HistType::kTH3F, {{axisMassD}, {axisPtD}, {axisCorrelationState}}});
    registry.add("hTransverse", "Transverse", {HistType::kTH3F, {{axisMassD}, {axisPtD}, {axisCorrelationState}}});
    registry.add("hAway", "Away", {HistType::kTH3F, {{axisMassD}, {axisPtD}, {axisCorrelationState}}});

    // Histograms for MC reco
    registry.add("hCorrel2DVsPtRecSig", "D0-h correlations MC reco signal", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.get<THnSparse>(HIST("hCorrel2DVsPtRecSig"))->Sumw2();
    registry.add("hCorrel2DVsPtRecBg", "D0-h correlations MC reco background", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.get<THnSparse>(HIST("hCorrel2DVsPtRecBg"))->Sumw2();
    // MC reco D0, D0bar signal case
    registry.add("hDeltaEtaPtIntSignalRegionRecSig", "D0-h deltaEta MC reco signal region, signal", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntSignalRegionRecSig", "D0-h deltaPhi MC reco signal region, signal", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntSignalRegionRecSig", "D0-h deltaPhi vs deltaEta MC reco signal region, signal", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtSignalRegionRecSig", "D0-h correlations MC reco signal region, signal", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisD0Prompt}, {axisPoolBin}}});
    registry.add("hCorrel2DVsPtPhysicalPrimaryRecSig", "D0-h correlations signal region (only true primary particles) MC reco", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisD0Prompt}, {axisPoolBin}}});
    registry.add("hCorrel2DVsPtSignalRegionPromptD0PromptHadronRecSig", "D0-h correlations signal region Prompt-NonPrompt MC reco", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hCorrel2DVsPtSignalRegionNonPromptD0NonPromptHadronRecSig", "D0-h correlations signal region NonPrompt-NonPrompt MC reco", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hDeltaEtaPtIntSignalRegionSoftPiRecSig", "D0-h deltaEta MC reco signal region, signal", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntSignalRegionSoftPiRecSig", "D0-h deltaPhi MC reco signal region, signal", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntSignalRegionSoftPiRecSig", "D0-h deltaPhi vs deltaEta MC reco signal region, signal", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtSignalRegionSoftPiRecSig", "D0-h correlations MC reco signal region, signal", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hDeltaEtaPtIntSidebandsRecSig", "D0-h deltaEta MC reco sidebands, signal", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntSidebandsRecSig", "D0-h deltaPhi MC reco sidebands, signal", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntSidebandsRecSig", "D0-h deltaPhi vs deltaEta MC reco sidebands, signal", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtSidebandsRecSig", "D0-h correlations MC reco sidebands, signal", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hDeltaEtaPtIntSidebandsSoftPiRecSig", "D0-h deltaEta MC reco sidebands, signal", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntSidebandsSoftPiRecSig", "D0-h deltaPhi MC reco sidebands, signal", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntSidebandsSoftPiRecSig", "D0-h deltaPhi vs deltaEta MC reco sidebands, signal", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtSidebandsSoftPiRecSig", "D0-h correlations MC reco sidebands, signal", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionRecSig"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtPhysicalPrimaryRecSig"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionSoftPiRecSig"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsRecSig"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsSoftPiRecSig"))->Sumw2();
    // MC reco D0, D0bar reflection case
    registry.add("hDeltaEtaPtIntSignalRegionRecRef", "D0-h deltaEta MC reco signal region, refelction", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntSignalRegionRecRef", "D0-h deltaPhi MC reco signal region, refelction", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntSignalRegionRecRef", "D0-h deltaPhi vs deltaEta MC reco signal region, refelction", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtSignalRegionRecRef", "D0-h correlations MC reco signal region, refelction", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hDeltaEtaPtIntSignalRegionSoftPiRecRef", "D0-h deltaEta MC reco signal region, refelction", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntSignalRegionSoftPiRecRef", "D0-h deltaPhi MC reco signal region, refelction", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntSignalRegionSoftPiRecRef", "D0-h deltaPhi vs deltaEta MC reco signal region, refelction", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtSignalRegionSoftPiRecRef", "D0-h correlations MC reco signal region, refelction", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hDeltaEtaPtIntSidebandsRecRef", "D0-h deltaEta MC reco sidebands, refelction", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntSidebandsRecRef", "D0-h deltaPhi MC reco sidebands, refelction", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntSidebandsRecRef", "D0-h deltaPhi vs deltaEta MC reco sidebands, refelction", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtSidebandsRecRef", "D0-h correlations MC reco sidebands, refelction", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hDeltaEtaPtIntSidebandsSoftPiRecRef", "D0-h deltaEta MC reco sidebands, refelction", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntSidebandsSoftPiRecRef", "D0-h deltaPhi MC reco sidebands, refelction", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntSidebandsSoftPiRecRef", "D0-h deltaPhi vs deltaEta MC reco sidebands, refelction", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtSidebandsSoftPiRecRef", "D0-h correlations MC reco sidebands, refelction", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionRecRef"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsRecRef"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionSoftPiRecRef"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsSoftPiRecRef"))->Sumw2();
    // MC reco D0, D0bar background case
    registry.add("hDeltaEtaPtIntSignalRegionRecBg", "D0-h deltaEta MC reco signal region, background", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntSignalRegionRecBg", "D0-h deltaPhi MC reco signal region, background", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntSignalRegionRecBg", "D0-h deltaPhi vs deltaEta MC reco signal region, background", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtSignalRegionRecBg", "D0-h correlations MC reco signal region, background", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hDeltaEtaPtIntSignalRegionSoftPiRecBg", "D0-h deltaEta MC reco signal region, background", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntSignalRegionSoftPiRecBg", "D0-h deltaPhi MC reco signal region, background", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntSignalRegionSoftPiRecBg", "D0-h deltaPhi vs deltaEta MC reco signal region, background", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtSignalRegionSoftPiRecBg", "D0-h correlations MC reco signal region, background", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hDeltaEtaPtIntSidebandsRecBg", "D0-h deltaEta MC reco sidebands, background", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntSidebandsRecBg", "D0-h deltaPhi MC reco sidebands, background", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntSidebandsRecBg", "D0-h deltaPhi vs deltaEta MC reco sidebands, background", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtSidebandsRecBg", "D0-h correlations MC reco sidebands, background", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hDeltaEtaPtIntSidebandsSoftPiRecBg", "D0-h deltaEta MC reco sidebands, background", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntSidebandsSoftPiRecBg", "D0-h deltaPhi MC reco sidebands, background", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntSidebandsSoftPiRecBg", "D0-h deltaPhi vs deltaEta MC reco sidebands, background", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtSidebandsSoftPiRecBg", "D0-h correlations MC reco sidebands, background", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionRecBg"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsRecBg"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSignalRegionSoftPiRecBg"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtSidebandsSoftPiRecBg"))->Sumw2();
    // Toward Transverse Away for McRec
    registry.add("hTowardRec", "Toward", {HistType::kTH3F, {{axisMassD}, {axisPtD}, {axisCorrelationState}}});
    registry.add("hTransverseRec", "Transverse", {HistType::kTH3F, {{axisMassD}, {axisPtD}, {axisCorrelationState}}});
    registry.add("hAwayRec", "Away", {HistType::kTH3F, {{axisMassD}, {axisPtD}, {axisCorrelationState}}});

    // Histograms for MC gen
    registry.add("hDeltaEtaPtIntGen", "D0-h deltaEta MC gen", {HistType::kTH1F, {axisDeltaEta}});
    registry.add("hDeltaPhiPtIntGen", "D0-h deltaPhi MC gen", {HistType::kTH1F, {axisDeltaPhi}});
    registry.add("hCorrel2DPtIntGen", "D0-h deltaPhi vs deltaEta MC gen", {HistType::kTH2F, {{axisDeltaPhi}, {axisDeltaEta}}});
    registry.add("hCorrel2DVsPtGenPrompt", "D0-h correlations MC gen, prompt", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hCorrel2DVsPtGenNonPrompt", "D0-h correlations MC gen, non prompt", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hCorrel2DVsPtGenPromptD0PromptHadron", "D0-h correlations MC gen, prompt D0, non-prompt hadrons", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.add("hCorrel2DVsPtGenNonPromptD0NonPromptHadron", "D0-h correlations MC gen, prompt D0, non-prompt hadrons", {HistType::kTHnSparseD, {{axisDeltaPhi}, {axisDeltaEta}, {axisPtD}, {axisPtHadron}, {axisPoolBin}}});
    registry.get<THnSparse>(HIST("hCorrel2DVsPtGenPrompt"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtGenNonPrompt"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtGenPromptD0PromptHadron"))->Sumw2();
    registry.get<THnSparse>(HIST("hCorrel2DVsPtGenNonPromptD0NonPromptHadron"))->Sumw2();
    // Toward Transverse Away for MC gen
    registry.add("hTowardGen", "Toward", {HistType::kTH3F, {{axisMassD}, {axisPtD}, {axisCorrelationState}}});
    registry.add("hTransverseGen", "Transverse", {HistType::kTH3F, {{axisMassD}, {axisPtD}, {axisCorrelationState}}});
    registry.add("hAwayGen", "Away", {HistType::kTH3F, {{axisMassD}, {axisPtD}, {axisCorrelationState}}});

    // Common histograms
    registry.add("hBdtScorePromptD0", "D0 BDT prompt score", {HistType::kTH1F, {axisBdtScore}});
    registry.add("hBdtScoreBkgD0", "D0 BDT bkg score", {HistType::kTH1F, {axisBdtScore}});
    registry.add("hBdtScorePromptD0bar", "D0bar BDT prompt score", {HistType::kTH1F, {axisBdtScore}});
    registry.add("hBdtScoreBkgD0bar", "D0bar BDT bkg score", {HistType::kTH1F, {axisBdtScore}});

    // Efficiency histograms
    registry.add("hPtCandMcRecPrompt", "D0 prompt candidates pt", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtCandMcRecNonPrompt", "D0 non prompt candidates pt", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtCandMcGenPrompt", "D0,Hadron particles prompt - MC Gen", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtCandMcGenNonPrompt", "D0,Hadron particles non prompt - MC Gen", {HistType::kTH1F, {axisPtD}});
    registry.add("hPtCandMcGenDaughterInAcc", "D0,Hadron particles non prompt - MC Gen", {HistType::kTH1F, {axisPtD}});

    auto hCandidates = registry.add<StepTHn>("hCandidates", "Candidate count at different steps", {HistType::kStepTHnF, {axisPtD, axisMultFT0M, {RecoDecay::OriginType::NonPrompt + 1, +RecoDecay::OriginType::None - 0.5, +RecoDecay::OriginType::NonPrompt + 0.5}}, kCandidateNSteps});
    hCandidates->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hCandidates->GetAxis(1)->SetTitle("multiplicity");
    hCandidates->GetAxis(2)->SetTitle("Charm hadron origin");
  }

  /// D-h correlation pair filling task, from pair tables - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  /// Works on both USL and LS analyses pair tables
  void processData(D0HadronPairFullMl const& pairEntries,
                   aod::D0CandRecoInfo const& candidates)
  {
    for (const auto& candidate : candidates) {
      float const ptD = candidate.ptD();
      float const bdtScorePromptD0 = candidate.mlScorePromptD0();
      float const bdtScoreBkgD0 = candidate.mlScoreBkgD0();
      float const bdtScorePromptD0bar = candidate.mlScorePromptD0bar();
      float const bdtScoreBkgD0bar = candidate.mlScoreBkgD0bar();
      int const effBinD = o2::analysis::findBin(binsPtEfficiencyD, ptD);

      registry.fill(HIST("hBdtScorePromptD0"), bdtScorePromptD0);
      registry.fill(HIST("hBdtScoreBkgD0"), bdtScoreBkgD0);
      registry.fill(HIST("hBdtScorePromptD0bar"), bdtScorePromptD0bar);
      registry.fill(HIST("hBdtScoreBkgD0bar"), bdtScoreBkgD0bar);

      if (bdtScorePromptD0 < mlOutputPromptD0->at(effBinD) || bdtScoreBkgD0 > mlOutputBkgD0->at(effBinD) ||
          bdtScorePromptD0bar < mlOutputPromptD0bar->at(effBinD) || bdtScoreBkgD0bar > mlOutputBkgD0bar->at(effBinD)) {
        continue;
      }
    }

    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double const deltaPhi = pairEntry.deltaPhi();
      double const deltaEta = pairEntry.deltaEta();
      double const ptD = pairEntry.ptD();
      double const ptHadron = pairEntry.ptHadron();
      double const massD = pairEntry.mD();
      double const massDbar = pairEntry.mDbar();
      int const signalStatus = pairEntry.signalStatus();
      int const effBinD = o2::analysis::findBin(binsPtEfficiencyD, ptD);
      int const ptBinD = o2::analysis::findBin(binsCorrelations, ptD);
      int const poolBin = pairEntry.poolBin();
      float const trackDcaXY = pairEntry.trackDcaXY();
      float const trackDcaZ = pairEntry.trackDcaZ();
      int const trackTpcCrossedRows = pairEntry.trackTPCNClsCrossedRows();
      bool const isAutoCorrelated = pairEntry.isAutoCorrelated();
      float const bdtScorePromptD0 = pairEntry.mlScorePromptD0();
      float const bdtScoreBkgD0 = pairEntry.mlScoreBkgD0();
      float const bdtScorePromptD0bar = pairEntry.mlScorePromptD0bar();
      float const bdtScoreBkgD0bar = pairEntry.mlScoreBkgD0bar();
      // reject entries outside pT ranges of interest
      if (ptBinD < 0 || effBinD < 0) {
        continue;
      }
      if (bdtScorePromptD0 < mlOutputPromptD0->at(ptBinD) || bdtScoreBkgD0 > mlOutputBkgD0->at(ptBinD) ||
          bdtScorePromptD0bar < mlOutputPromptD0bar->at(ptBinD) || bdtScoreBkgD0bar > mlOutputBkgD0bar->at(ptBinD)) {
        continue;
      }
      if (trackDcaXY > dcaXYTrackMax || trackDcaZ > dcaZTrackMax || trackTpcCrossedRows < nTpcCrossedRaws) {
        continue;
      }
      double efficiencyWeight = 1.;
      if (applyEfficiency != 0) {
        efficiencyWeight = 1. / (efficiencyDmeson->at(o2::analysis::findBin(binsPtEfficiencyD, ptD)) * efficiencyHad->at(o2::analysis::findBin(binsPtEfficiencyHad, ptHadron)));
      }
      // reject entries outside pT ranges of interest
      if (ptBinD == -1) { // at least one particle outside accepted pT range
        continue;
      }
      //==============================================================================================================
      if (isTowardTransverseAway) {
        // Divide into three regions: toward, transverse, and away
        if (ptHadron < leadingParticlePtMin) {
          continue;
        }
        Region const region = getRegion(deltaPhi);

        switch (region) {
          case Toward:
            registry.fill(HIST("hToward"), massD, ptD, isAutoCorrelated, efficiencyWeight);
            break;
          case Away:
            registry.fill(HIST("hAway"), massD, ptD, isAutoCorrelated, efficiencyWeight);
            break;
          case Transverse:
            registry.fill(HIST("hTransverse"), massD, ptD, isAutoCorrelated, efficiencyWeight);
            break;
          default:
            break;
        }
      }
      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots
      if ((massD > signalRegionLeft->at(ptBinD) && massD < signalRegionRight->at(ptBinD)) && (signalStatus == ParticleTypeData::D0Only)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegion"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
      }

      if ((massD > signalRegionLeft->at(ptBinD) && massD < signalRegionRight->at(ptBinD)) && (signalStatus == ParticleTypeData::D0OnlySoftPi)) {
        // in signal region, fills for soft pion only in ME
        registry.fill(HIST("hCorrel2DVsPtSignalRegionSoftPi"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionSoftPi"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionSoftPi"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionSoftPi"), deltaPhi, efficiencyWeight);
      }

      if ((massDbar > signalRegionLeft->at(ptBinD) && massDbar < signalRegionRight->at(ptBinD)) && (signalStatus == ParticleTypeData::D0barOnly)) {
        // in signal region
        registry.fill(HIST("hCorrel2DVsPtSignalRegion"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegion"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegion"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegion"), deltaPhi, efficiencyWeight);
      }

      if ((massDbar > signalRegionLeft->at(ptBinD) && massDbar < signalRegionRight->at(ptBinD)) && (signalStatus >= ParticleTypeData::D0barOnlySoftPi)) {
        // in signal region, fills for soft pion only in ME
        registry.fill(HIST("hCorrel2DVsPtSignalRegionSoftPi"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionSoftPi"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionSoftPi"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionSoftPi"), deltaPhi, efficiencyWeight);
      }

      if (((massD > sidebandLeftOuter->at(ptBinD) && massD < sidebandLeftInner->at(ptBinD)) ||
           (massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD))) &&
          (signalStatus == ParticleTypeData::D0Only)) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebands"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi, efficiencyWeight);
      }

      if (((massD > sidebandLeftOuter->at(ptBinD) && massD < sidebandLeftInner->at(ptBinD)) ||
           (massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD))) &&
          (signalStatus == ParticleTypeData::D0OnlySoftPi)) {
        // in sideband region, fills for soft pion only in ME
        registry.fill(HIST("hCorrel2DVsPtSidebandsSoftPi"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsSoftPi"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsSoftPi"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsSoftPi"), deltaPhi, efficiencyWeight);
      }

      if (((massDbar > sidebandLeftOuter->at(ptBinD) && massDbar < sidebandLeftInner->at(ptBinD)) ||
           (massDbar > sidebandRightInner->at(ptBinD) && massDbar < sidebandRightOuter->at(ptBinD))) &&
          (signalStatus == ParticleTypeData::D0barOnly)) {
        // in sideband region
        registry.fill(HIST("hCorrel2DVsPtSidebands"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebands"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebands"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebands"), deltaPhi, efficiencyWeight);
      }

      if (((massDbar > sidebandLeftOuter->at(ptBinD) && massDbar < sidebandLeftInner->at(ptBinD)) ||
           (massDbar > sidebandRightInner->at(ptBinD) && massDbar < sidebandRightOuter->at(ptBinD))) &&
          (signalStatus >= ParticleTypeData::D0barOnlySoftPi)) {
        // in sideband region, fills for soft pion only in ME
        registry.fill(HIST("hCorrel2DVsPtSidebandsSoftPi"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsSoftPi"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsSoftPi"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsSoftPi"), deltaPhi, efficiencyWeight);
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationD0Hadrons, processData, "Process data", false);

  void processMcRec(D0HadronPairFullMl const& pairEntries,
                    soa::Join<aod::D0CandRecoInfo, aod::D0CandGenInfo> const& candidates)
  {
    for (const auto& candidate : candidates) {
      float const ptD = candidate.ptD();
      float const bdtScorePromptD0 = candidate.mlScorePromptD0();
      float const bdtScoreBkgD0 = candidate.mlScoreBkgD0();
      float const bdtScorePromptD0bar = candidate.mlScorePromptD0bar();
      float const bdtScoreBkgD0bar = candidate.mlScoreBkgD0bar();
      int const effBinD = o2::analysis::findBin(binsPtEfficiencyD, ptD);

      registry.fill(HIST("hBdtScorePromptD0"), bdtScorePromptD0);
      registry.fill(HIST("hBdtScoreBkgD0"), bdtScoreBkgD0);
      registry.fill(HIST("hBdtScorePromptD0bar"), bdtScorePromptD0bar);
      registry.fill(HIST("hBdtScoreBkgD0bar"), bdtScoreBkgD0bar);

      if (bdtScorePromptD0 < mlOutputPromptD0->at(effBinD) || bdtScoreBkgD0 > mlOutputBkgD0->at(effBinD) ||
          bdtScorePromptD0bar < mlOutputPromptD0bar->at(effBinD) || bdtScoreBkgD0bar > mlOutputBkgD0bar->at(effBinD)) {
        continue;
      }
    }

    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double const deltaPhi = pairEntry.deltaPhi();
      double const deltaEta = pairEntry.deltaEta();
      double const ptD = pairEntry.ptD();
      double const ptHadron = pairEntry.ptHadron();
      double const massD = pairEntry.mD();
      double const massDbar = pairEntry.mDbar();
      int const signalStatus = pairEntry.signalStatus();
      int const ptBinD = o2::analysis::findBin(binsCorrelations, ptD);
      int const poolBin = pairEntry.poolBin();
      float const trackDcaXY = pairEntry.trackDcaXY();
      float const trackDcaZ = pairEntry.trackDcaZ();
      int const trackTpcCrossedRows = pairEntry.trackTPCNClsCrossedRows();
      bool const isAutoCorrelated = pairEntry.isAutoCorrelated();
      float const bdtScorePromptD0 = pairEntry.mlScorePromptD0();
      float const bdtScoreBkgD0 = pairEntry.mlScoreBkgD0();
      float const bdtScorePromptD0bar = pairEntry.mlScorePromptD0bar();
      float const bdtScoreBkgD0bar = pairEntry.mlScoreBkgD0bar();
      bool const isPhysicalPrimary = pairEntry.isPhysicalPrimary();
      bool const isD0Prompt = pairEntry.isPrompt();
      int const statusPromptHadron = pairEntry.trackOrigin();

      if (bdtScorePromptD0 < mlOutputPromptD0->at(ptBinD) || bdtScoreBkgD0 > mlOutputBkgD0->at(ptBinD) ||
          bdtScorePromptD0bar < mlOutputPromptD0bar->at(ptBinD) || bdtScoreBkgD0bar > mlOutputBkgD0bar->at(ptBinD)) {
        continue;
      }
      if (trackDcaXY > dcaXYTrackMax || trackDcaZ > dcaZTrackMax || trackTpcCrossedRows < nTpcCrossedRaws) {
        continue;
      }
      double efficiencyWeight = 1.;
      if (applyEfficiency != 0) {
        efficiencyWeight = 1. / (efficiencyDmeson->at(o2::analysis::findBin(binsPtEfficiencyD, ptD)) * efficiencyHad->at(o2::analysis::findBin(binsPtEfficiencyHad, ptHadron)));
      }
      if (isTowardTransverseAway) {
        // Divide into three regions: toward, transverse, and away
        if (ptHadron < leadingParticlePtMin) {
          continue;
        }
        Region const region = getRegion(deltaPhi);
        switch (region) {
          case Toward:
            registry.fill(HIST("hTowardRec"), massD, ptD, isAutoCorrelated, efficiencyWeight);
            break;
          case Away:
            registry.fill(HIST("hAwayRec"), massD, ptD, isAutoCorrelated, efficiencyWeight);
            break;
          case Transverse:
            registry.fill(HIST("hTransverseRec"), massD, ptD, isAutoCorrelated, efficiencyWeight);
            break;
          default:
            break;
        }
      }
      // fill correlation plots for signal/bagkground correlations
      if (pairEntry.signalStatus() != 0) {
        registry.fill(HIST("hCorrel2DVsPtRecSig"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);

      } else {
        registry.fill(HIST("hCorrel2DVsPtRecBg"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
      }

      // check if correlation entry belongs to signal region, sidebands or is outside both, and fill correlation plots

      // ---------------------- Fill plots for signal case, D0 ->1, D0bar ->8 ---------------------------------------------
      if ((massD > signalRegionLeft->at(ptBinD) && massD < signalRegionRight->at(ptBinD)) && (TESTBIT(signalStatus, ParticleTypeMcRec::D0Sig))) {
        // in signal region, tests bit ParticleTypeMcRec::D0Sig, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSignalRegionRecSig"), deltaPhi, deltaEta, ptD, ptHadron, static_cast<int>(isD0Prompt), poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionRecSig"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionRecSig"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionRecSig"), deltaPhi, efficiencyWeight);
        if (isPhysicalPrimary) {
          registry.fill(HIST("hCorrel2DVsPtPhysicalPrimaryRecSig"), deltaPhi, deltaEta, ptD, ptHadron, static_cast<int>(isD0Prompt), poolBin, efficiencyWeight);
          if (isD0Prompt && statusPromptHadron == RecoDecay::OriginType::Prompt) {
            registry.fill(HIST("hCorrel2DVsPtSignalRegionPromptD0PromptHadronRecSig"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
          } else if (!isD0Prompt && statusPromptHadron == RecoDecay::OriginType::NonPrompt) {
            registry.fill(HIST("hCorrel2DVsPtSignalRegionNonPromptD0NonPromptHadronRecSig"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
          }
        }
      }

      if ((massDbar > signalRegionLeft->at(ptBinD) && massDbar < signalRegionRight->at(ptBinD)) && (TESTBIT(signalStatus, ParticleTypeMcRec::D0barSig))) {
        // in signal region, tests bit ParticleTypeMcRec::D0barSig, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSignalRegionRecSig"), deltaPhi, deltaEta, ptD, ptHadron, static_cast<int>(isD0Prompt), poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionRecSig"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionRecSig"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionRecSig"), deltaPhi, efficiencyWeight);
        if (isPhysicalPrimary) {
          registry.fill(HIST("hCorrel2DVsPtPhysicalPrimaryRecSig"), deltaPhi, deltaEta, ptD, ptHadron, static_cast<int>(isD0Prompt), poolBin, efficiencyWeight);
          if (isD0Prompt && statusPromptHadron == RecoDecay::OriginType::Prompt) {
            registry.fill(HIST("hCorrel2DVsPtSignalRegionPromptD0PromptHadronRecSig"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
          } else if (!isD0Prompt && statusPromptHadron == RecoDecay::OriginType::NonPrompt) {
            registry.fill(HIST("hCorrel2DVsPtSignalRegionNonPromptD0NonPromptHadronRecSig"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
          }
        }
      }

      if ((massDbar > signalRegionLeft->at(ptBinD) && massDbar < signalRegionRight->at(ptBinD)) && (signalStatus >= static_cast<int>(BIT(ParticleTypeMcRec::SoftPi)))) {
        // in signal region, fills for soft pion only for event mixing, tests bit at least ParticleTypeMcRec::SoftPi
        registry.fill(HIST("hCorrel2DVsPtSignalRegionSoftPiRecSig"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionSoftPiRecSig"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionSoftPiRecSig"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionSoftPiRecSig"), deltaPhi, efficiencyWeight);
      }

      if ((((massD > sidebandLeftOuter->at(ptBinD)) && (massD < sidebandLeftInner->at(ptBinD))) ||
           ((massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD)))) &&
          ((TESTBIT(signalStatus, ParticleTypeMcRec::D0Sig)))) {
        // in sideband region, tests bit ParticleTypeMcRec::D0Sig, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSidebandsRecSig"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsRecSig"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsRecSig"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsRecSig"), deltaPhi, efficiencyWeight);
      }

      if ((((massDbar > sidebandLeftOuter->at(ptBinD)) && (massDbar < sidebandLeftInner->at(ptBinD))) ||
           ((massDbar > sidebandRightInner->at(ptBinD) && massDbar < sidebandRightOuter->at(ptBinD)))) &&
          (TESTBIT(signalStatus, ParticleTypeMcRec::D0barSig))) {
        // in sideband region, tests bit ParticleTypeMcRec::D0barSig, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSidebandsRecSig"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsRecSig"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsRecSig"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsRecSig"), deltaPhi, efficiencyWeight);
      }

      if ((((massDbar > sidebandLeftOuter->at(ptBinD)) && (massDbar < sidebandLeftInner->at(ptBinD))) ||
           ((massDbar > sidebandRightInner->at(ptBinD) && massDbar < sidebandRightOuter->at(ptBinD)))) &&
          (signalStatus >= static_cast<int>(BIT(ParticleTypeMcRec::SoftPi)))) {
        // in sideband region, fills for soft pion only for event mixing, tests bit at least ParticleTypeMcRec::SoftPi
        registry.fill(HIST("hCorrel2DVsPtSidebandsSoftPiRecSig"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsSoftPiRecSig"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsSoftPiRecSig"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsSoftPiRecSig"), deltaPhi, efficiencyWeight);
      }

      // ---------------------- Fill plots for reflection case, D0 ->2, D0bar ->16 ---------------------------------------------
      if ((massD > signalRegionLeft->at(ptBinD) && massD < signalRegionRight->at(ptBinD)) && TESTBIT(signalStatus, ParticleTypeMcRec::D0Ref)) {
        // in signal region, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSignalRegionRecRef"), deltaPhi, deltaEta, ptD, ptHadron, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionRecRef"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionRecRef"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionRecRef"), deltaPhi, efficiencyWeight);
      }

      if ((massDbar > signalRegionLeft->at(ptBinD) && massDbar < signalRegionRight->at(ptBinD)) && TESTBIT(signalStatus, ParticleTypeMcRec::D0barRef)) {
        // in signal region, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSignalRegionRecRef"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionRecRef"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionRecRef"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionRecRef"), deltaPhi, efficiencyWeight);
      }

      if ((massDbar > signalRegionLeft->at(ptBinD) && massDbar < signalRegionRight->at(ptBinD)) && (signalStatus >= static_cast<int>(BIT(ParticleTypeMcRec::SoftPi)))) {
        // in signal region, fills for soft pion only for event mixing, tests bit at least ParticleTypeMcRec::SoftPi
        registry.fill(HIST("hCorrel2DVsPtSignalRegionSoftPiRecRef"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionSoftPiRecRef"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionSoftPiRecRef"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionSoftPiRecRef"), deltaPhi, efficiencyWeight);
      }

      if ((((massD > sidebandLeftOuter->at(ptBinD)) && (massD < sidebandLeftInner->at(ptBinD))) ||
           ((massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD)))) &&
          TESTBIT(signalStatus, ParticleTypeMcRec::D0Ref)) {
        // in sideband region, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSidebandsRecRef"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsRecRef"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsRecRef"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsRecRef"), deltaPhi, efficiencyWeight);
      }

      if ((((massDbar > sidebandLeftOuter->at(ptBinD)) && (massDbar < sidebandLeftInner->at(ptBinD))) ||
           ((massDbar > sidebandRightInner->at(ptBinD) && massDbar < sidebandRightOuter->at(ptBinD)))) &&
          TESTBIT(signalStatus, ParticleTypeMcRec::D0barRef)) {
        // in sideband region, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSidebandsRecRef"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsRecRef"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsRecRef"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsRecRef"), deltaPhi, efficiencyWeight);
      }

      if ((((massDbar > sidebandLeftOuter->at(ptBinD)) && (massDbar < sidebandLeftInner->at(ptBinD))) ||
           ((massDbar > sidebandRightInner->at(ptBinD) && massDbar < sidebandRightOuter->at(ptBinD)))) &&
          (signalStatus >= static_cast<int>(BIT(ParticleTypeMcRec::SoftPi)))) {
        // in sideband region, fills for soft pion only for event mixing, tests bit at least ParticleTypeMcRec::SoftPi
        registry.fill(HIST("hCorrel2DVsPtSidebandsSoftPiRecRef"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsSoftPiRecRef"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsSoftPiRecRef"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsSoftPiRecRef"), deltaPhi, efficiencyWeight);
      }

      // ---------------------- Fill plots for background case, D0 ->4, D0bar ->32 ---------------------------------------------
      if ((massD > signalRegionLeft->at(ptBinD) && massD < signalRegionRight->at(ptBinD)) && TESTBIT(signalStatus, ParticleTypeMcRec::D0Bg)) {
        // in signal region, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSignalRegionRecBg"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionRecBg"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionRecBg"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionRecBg"), deltaPhi, efficiencyWeight);
      }

      if ((massDbar > signalRegionLeft->at(ptBinD) && massDbar < signalRegionRight->at(ptBinD)) && TESTBIT(signalStatus, ParticleTypeMcRec::D0barBg)) {
        // in signal region, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSignalRegionRecBg"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionRecBg"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionRecBg"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionRecBg"), deltaPhi, efficiencyWeight);
      }

      if ((massDbar > signalRegionLeft->at(ptBinD) && massDbar < signalRegionRight->at(ptBinD)) && (signalStatus >= static_cast<int>(BIT(ParticleTypeMcRec::SoftPi)))) {
        // in signal region, fills for soft pion only for event mixing, tests bit at least ParticleTypeMcRec::SoftPi
        registry.fill(HIST("hCorrel2DVsPtSignalRegionSoftPiRecBg"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSignalRegionSoftPiRecBg"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSignalRegionSoftPiRecBg"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSignalRegionSoftPiRecBg"), deltaPhi, efficiencyWeight);
      }

      if ((((massD > sidebandLeftOuter->at(ptBinD)) && (massD < sidebandLeftInner->at(ptBinD))) ||
           ((massD > sidebandRightInner->at(ptBinD) && massD < sidebandRightOuter->at(ptBinD)))) &&
          TESTBIT(signalStatus, ParticleTypeMcRec::D0Bg)) {
        // in sideband region, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSidebandsRecBg"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsRecBg"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsRecBg"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsRecBg"), deltaPhi, efficiencyWeight);
      }

      if ((((massDbar > sidebandLeftOuter->at(ptBinD)) && (massDbar < sidebandLeftInner->at(ptBinD))) ||
           ((massDbar > sidebandRightInner->at(ptBinD) && massDbar < sidebandRightOuter->at(ptBinD)))) &&
          TESTBIT(signalStatus, ParticleTypeMcRec::D0barBg)) {
        // in sideband region, SE-> softpi removed, ME-> inclusive
        registry.fill(HIST("hCorrel2DVsPtSidebandsRecBg"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsRecBg"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsRecBg"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsRecBg"), deltaPhi, efficiencyWeight);
      }

      if ((((massDbar > sidebandLeftOuter->at(ptBinD)) && (massDbar < sidebandLeftInner->at(ptBinD))) ||
           ((massDbar > sidebandRightInner->at(ptBinD) && massDbar < sidebandRightOuter->at(ptBinD)))) &&
          (signalStatus >= static_cast<int>(BIT(ParticleTypeMcRec::SoftPi)))) {
        // in sideband region, fills for soft pion only for event mixing, tests bit at least ParticleTypeMcRec::SoftPi
        registry.fill(HIST("hCorrel2DVsPtSidebandsSoftPiRecBg"), deltaPhi, deltaEta, ptD, ptHadron, poolBin, efficiencyWeight);
        registry.fill(HIST("hCorrel2DPtIntSidebandsSoftPiRecBg"), deltaPhi, deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaEtaPtIntSidebandsSoftPiRecBg"), deltaEta, efficiencyWeight);
        registry.fill(HIST("hDeltaPhiPtIntSidebandsSoftPiRecBg"), deltaPhi, efficiencyWeight);
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationD0Hadrons, processMcRec, "Process MC Reco mode", true);

  /// D-Hadron correlation pair filling task, from pair tables - for MC gen-level analysis (no filter/selection, only true signal)
  void processMcGen(D0HadronPairFull const& pairEntries)
  {
    for (const auto& pairEntry : pairEntries) {
      // define variables for widely used quantities
      double const deltaPhi = pairEntry.deltaPhi();
      double const deltaEta = pairEntry.deltaEta();
      double const ptD = pairEntry.ptD();
      double const ptHadron = pairEntry.ptHadron();
      int const poolBin = pairEntry.poolBin();
      double const massD = pairEntry.mD();
      bool const isAutoCorrelated = pairEntry.isAutoCorrelated();
      int const statusPromptHadron = pairEntry.trackOrigin();
      bool const isD0Prompt = pairEntry.isPrompt();

      // reject entries outside pT ranges of interest
      if (o2::analysis::findBin(binsCorrelations, ptD) < 0) {
        continue;
      }
      if (isTowardTransverseAway) {
        // Divide into three regions: toward, transverse, and away
        if (ptHadron < leadingParticlePtMin) {
          continue;
        }
        Region const region = getRegion(deltaPhi);
        switch (region) {
          case Toward:
            registry.fill(HIST("hTowardGen"), massD, ptD, isAutoCorrelated);
            break;
          case Away:
            registry.fill(HIST("hAwayGen"), massD, ptD, isAutoCorrelated);
            break;
          case Transverse:
            registry.fill(HIST("hTransverseGen"), massD, ptD, isAutoCorrelated);
            break;
          default:
            break;
        }
      }
      if (isD0Prompt) {
        registry.fill(HIST("hCorrel2DVsPtGenPrompt"), deltaPhi, deltaEta, ptD, ptHadron, poolBin);
        if (statusPromptHadron == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("hCorrel2DVsPtGenPromptD0PromptHadron"), deltaPhi, deltaEta, ptD, ptHadron, poolBin);
        }
      } else {
        registry.fill(HIST("hCorrel2DVsPtGenNonPrompt"), deltaPhi, deltaEta, ptD, ptHadron, poolBin);
        if (statusPromptHadron == RecoDecay::OriginType::NonPrompt) {
          registry.fill(HIST("hCorrel2DVsPtGenNonPromptD0NonPromptHadron"), deltaPhi, deltaEta, ptD, ptHadron, poolBin);
        }
      }
      registry.fill(HIST("hCorrel2DPtIntGen"), deltaPhi, deltaEta);
      registry.fill(HIST("hDeltaEtaPtIntGen"), deltaEta);
      registry.fill(HIST("hDeltaPhiPtIntGen"), deltaPhi);
    } // end loop
  }
  PROCESS_SWITCH(HfTaskCorrelationD0Hadrons, processMcGen, "Process MC Gen mode", false);

  /// D0 reconstruction and selection efficiency
  void processMcCandEfficiency(soa::Join<aod::Collisions, aod::FT0Mults, aod::EvSels> const&,
                               soa::Join<aod::McCollisions, aod::MultsExtraMC> const&,
                               CandD0McGen const& mcParticles,
                               CandD0McReco const& candidates,
                               aod::TracksWMc const&)
  {
    auto hCandidates = registry.get<StepTHn>(HIST("hCandidates"));

    /// Gen loop
    float multiplicity = -1.;
    for (const auto& mcParticle : mcParticles) {
      // generated candidates
      if (std::abs(mcParticle.pdgCode()) == Pdg::kD0) {
        auto mcCollision = mcParticle.template mcCollision_as<soa::Join<aod::McCollisions, aod::MultsExtraMC>>();
        multiplicity = mcCollision.multMCFT0A() + mcCollision.multMCFT0C(); // multFT0M = multFt0A + multFT0C
        hCandidates->Fill(kCandidateStepMcGenAll, mcParticle.pt(), multiplicity, mcParticle.originMcGen());
        if (std::abs(mcParticle.flagMcMatchGen()) == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
          hCandidates->Fill(kCandidateStepMcGenD0ToPiKPi, mcParticle.pt(), multiplicity, mcParticle.originMcGen());
          auto yD0 = RecoDecay::y(mcParticle.pVector(), o2::constants::physics::MassD0);
          if (std::abs(yD0) <= yCandGenMax) {
            hCandidates->Fill(kCandidateStepMcCandInAcceptance, mcParticle.pt(), multiplicity, mcParticle.originMcGen());
            if (mcParticle.originMcGen() == RecoDecay::OriginType::Prompt) {
              registry.fill(HIST("hPtCandMcGenPrompt"), mcParticle.pt());
            }
            if (mcParticle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
              registry.fill(HIST("hPtCandMcGenNonPrompt"), mcParticle.pt());
            }
          }
          bool isDaughterInAcceptance = true;
          auto daughters = mcParticle.template daughters_as<CandD0McGen>();
          for (const auto& daughter : daughters) {
            if (daughter.pt() < ptDaughterMin || std::abs(daughter.eta()) > etaTrackMax) {
              isDaughterInAcceptance = false;
            }
          }
          if (isDaughterInAcceptance) {
            hCandidates->Fill(kCandidateStepMcDaughtersInAcceptance, mcParticle.pt(), multiplicity, mcParticle.originMcGen());
            registry.fill(HIST("hPtCandMcGenDaughterInAcc"), mcParticle.pt());
          }
        }
      }
    }

    // recontructed candidates loop
    for (const auto& candidate : candidates) {
      if (candidate.pt() < ptCandMin || candidate.pt() > ptCandMax) {
        continue;
      }
      std::vector<float> outputMlD0 = {-1., -1., -1.};
      std::vector<float> outputMlD0bar = {-1., -1., -1.};
      if (candidate.isSelD0() < selectionFlagD0 || candidate.isSelD0bar() < selectionFlagD0) {
        continue;
      }
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
        outputMlD0[iclass] = candidate.mlProbD0()[classMl->at(iclass)];
        outputMlD0bar[iclass] = candidate.mlProbD0bar()[classMl->at(iclass)];
      }
      if (outputMlD0[0] > mlOutputBkgD0->at(o2::analysis::findBin(binsPtEfficiencyD, candidate.pt())) || outputMlD0[1] < mlOutputPromptD0->at(o2::analysis::findBin(binsPtEfficiencyD, candidate.pt()))) {
        continue;
      }
      if (outputMlD0bar[0] > mlOutputBkgD0bar->at(o2::analysis::findBin(binsPtEfficiencyD, candidate.pt())) || outputMlD0bar[1] < mlOutputPromptD0bar->at(o2::analysis::findBin(binsPtEfficiencyD, candidate.pt()))) {
        continue;
      }
      auto collision = candidate.template collision_as<soa::Join<aod::Collisions, aod::FT0Mults, aod::EvSels>>();
      if (selNoSameBunchPileUpColl && !(collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))) {
        continue;
      }
      multiplicity = collision.multFT0M();
      if (std::abs(candidate.flagMcMatchRec()) == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
        hCandidates->Fill(kCandidateStepMcReco, candidate.pt(), multiplicity, candidate.originMcRec());
        if (std::abs(hfHelper.yD0(candidate)) <= yCandMax) {
          hCandidates->Fill(kCandidateStepMcRecoInAcceptance, candidate.pt(), multiplicity, candidate.originMcRec());
          if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
            registry.fill(HIST("hPtCandMcRecPrompt"), candidate.pt());
          }
          if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
            registry.fill(HIST("hPtCandMcRecNonPrompt"), candidate.pt());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationD0Hadrons, processMcCandEfficiency, "Process MC for calculating candidate reconstruction efficiency", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationD0Hadrons>(cfgc)};
}
