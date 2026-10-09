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

/// \file pidLongRange.cxx
/// \brief Identified particles long-range correlations using forward FIT detectors and TPC, with focus on particle
/// \author Thor Jensen (thor.kjaersgaard.jensen@cern.ch), Preet Bhanjan Pati (preet.bhanjan.pati@cern.ch)

#include "PWGCF/TwoParticleCorrelations/Core/DihadronContainer.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <DetectorsCommonDataFormats/AlignParam.h>
#include <FT0Base/Geometry.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/StepTHn.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>

#include <TF1.h>
#include <TFile.h>
#include <TH3.h>
#include <TRandom3.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;
using namespace constants::math;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, (DEFAULT), (HELP)}; // NOLINT(bugprone-macro-parentheses)

template <typename T, typename P>
auto readMatrix(Array2D<T> const& mat, P& array)
{
  for (auto i = 0; i < static_cast<int>(mat.rows); ++i) {
    for (auto j = 0; j < static_cast<int>(mat.cols); ++j) {
      array[i][j] = mat(i, j);
    }
  }

  return;
}

static constexpr std::array<std::array<float, 3>, 20> LongArrayFloat = {{{{1.1, 2.1, 3.1}}, {{1.2, 2.2, 3.2}}, {{1.3, 2.3, 3.3}}, {{-1.1, -2.1, -3.1}}, {{-1.2, -2.2, -3.2}}, {{-1.3, -2.3, -3.3}}, {{1.1, 1.1, 1.1}}, {{1.2, 1.2, 1.2}}, {{1.3, 1.3, 1.3}}, {{-1.1, -1.1, -1.1}}, {{-1.2, -1.2, -1.2}}, {{-1.3, -1.3, -1.3}}, {{1.1, 1.1, 1.1}}, {{1.2, 1.2, 1.2}}, {{1.3, 1.3, 1.3}}, {{-1.1, -1.1, -1.1}}, {{-1.2, -1.2, -1.2}}, {{-1.3, -1.3, -1.3}}, {{1.1, 1.1, 1.1}}, {{1.2, 1.2, 1.2}}}};
static constexpr std::array<std::array<int, 3>, 20> LongArrayInt = {{{{1, 2, 3}}, {{1, 2, 3}}, {{1, 2, 3}}, {{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}, {{1, 1, 1}}, {{1, 1, 1}}, {{1, 1, 1}}, {{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}, {{1, 1, 1}}, {{1, 1, 1}}, {{1, 1, 1}}, {{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}, {{1, 1, 1}}, {{1, 1, 1}}}};

struct PidLongRange {
  o2::aod::ITSResponse itsResponse;
  Service<ccdb::BasicCCDBManager> ccdb{};

  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgQABasic, bool, true, "Enable QA histograms for event and track selection")
    O2_DEFINE_CONFIGURABLE(cfgStrictTrackCounter, bool, false, "Strict track counter for multiplicity correlation cut, counts only tracks that pass all cuts and are used in the correlation")
    O2_DEFINE_CONFIGURABLE(cfgSampleSize, int, 10, "Sample size for mixed event")
    O2_DEFINE_CONFIGURABLE(cfgAnalyzeTPCFT0A, bool, true, "Switch for doing TPC-FT0A correlations")
    O2_DEFINE_CONFIGURABLE(cfgAnalyzeTPCFT0C, bool, true, "Switch for doing TPC-FT0C correlations")
    O2_DEFINE_CONFIGURABLE(cfgMinMixEventNum, int, 5, "Minimum number of events to mix")
    Configurable<std::vector<int>> cfgRunRemoveList{"cfgRunRemoveList", std::vector<int>{-1}, "excluded run numbers"};
    Configurable<std::string> cfgGainEqPath{"cfgGainEqPath", "Analysis/EventPlane/GainEq", "CCDB path for gain equalization constants"};
    Configurable<int> cfgCorrLevel{"cfgCorrLevel", 0, "calibration step: 0 = no corr, 1 = gain corr"};
  } cfgGeneral;

  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
    O2_DEFINE_CONFIGURABLE(cfgEfficiencyNch, std::string, "", "CCDB path to multiplicity dependent efficiency object")
    O2_DEFINE_CONFIGURABLE(cfgUseEventWeights, bool, false, "Use event weights for mixed event")
    O2_DEFINE_CONFIGURABLE(cfgCentralityWeight, std::string, "", "CCDB path to centrality weight object")
    O2_DEFINE_CONFIGURABLE(cfgLocalEfficiency, bool, false, "Use local efficiency object")
    O2_DEFINE_CONFIGURABLE(cfgLocalEfficiencyNch, bool, false, "Use local multiplicity dependent efficiency object");
  } cfgEffWeights;

  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgGlobalPtMin, float, 0.2f, "minimum accepted track pT")
    O2_DEFINE_CONFIGURABLE(cfgGlobalPtMax, float, 10.0f, "maximum accepted track pT")
    O2_DEFINE_CONFIGURABLE(cfgGlobalEta, float, 0.8f, "Eta cut")
    O2_DEFINE_CONFIGURABLE(cfgGlobalChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
    O2_DEFINE_CONFIGURABLE(cfgGlobalDCAz, float, 2.0f, "DCAz range for tracks")
    O2_DEFINE_CONFIGURABLE(cfgCentEstimator, int, 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FT0A")
    Configurable<std::vector<std::string>> cfgTrackCutsDCAxy{"cfgTrackCutsDCAxy", std::vector<std::string>{"(0.0026+0.005/(x^1.01))", "(0.0026+0.005/(x^1.01))"}, "Functional form of pt-dependent DCAxy cut"};
    Configurable<std::vector<std::string>> cfgTrackCutsDCAz{"cfgTrackCutsDCAz", std::vector<std::string>{"", ""}, "Functional form of pt-dependent DCAz cut"};
    Configurable<LabeledArray<float>> cfgTrackCuts{"cfgTrackCuts", {LongArrayFloat.front().data(), 10, 2, {"PtMin", "PtMax", "EtaCut", "Chi2PrTpcCls", "TpcCluster", "TpcCrossedRows", "ItsCluster", "DCAz", "DCAzNsigma", "DCAxyNsigma"}, {"General", "SelNch"}}, "Labeled array for track selections"};
  } cfgTrkSel;

  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgGlobalCutVertex, float, 10.0f, "Accepted z-vertex range")
    O2_DEFINE_CONFIGURABLE(cfgCutOccupancyHigh, int, 2000, "High cut on TPC occupancy")
    O2_DEFINE_CONFIGURABLE(cfgCutOccupancyLow, int, 0, "Low cut on TPC occupancy")
    O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, true, "Use additional event cut on mult correlations")
    O2_DEFINE_CONFIGURABLE(cfgMinMultForCorrelations, int, 0, "minimum multiplicity for correlations")
    O2_DEFINE_CONFIGURABLE(cfgMaxMultForCorrelations, int, 20, "maximum multiplicity for correlations")
    O2_DEFINE_CONFIGURABLE(cfgSelCollByNch, bool, true, "Select collisions by Nch or centrality")
    O2_DEFINE_CONFIGURABLE(cfgMinCentForCorrelations, int, 0, "minimum centrality for correlations")
    O2_DEFINE_CONFIGURABLE(cfgMaxCentForCorrelations, int, 20, "maximum centrality for correlations")
    O2_DEFINE_CONFIGURABLE(cfgEvSelRCTflags, std::string, "", "keep empty to disable, usage: 'CentralBarrelTracking (CBT)', 'CBT_hadronPID' ")
    Configurable<LabeledArray<int>> cfgUseEventCuts{"cfgUseEventCuts", {LongArrayInt.front().data(), 14, 1, {"Filtered Events", "Sel8", "kNoTimeFrameBorder", "kNoITSROFrameBorder", "kNoSameBunchPileup", "kIsGoodZvtxFT0vsPV", "kNoCollInTimeRangeStandard", "kIsGoodITSLayersAll", "kIsGoodITSLayer0123", "kNoCollInRofStandard", "kNoHighMultCollInPrevRof", "Occupancy", "Multcorrelation", "T0AV0ACut"}, {"EvCuts"}}, "Labeled array (int) for various cuts on resonances"};
  } cfgEvSel;

  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgMultCentHighCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 10.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultCentLowCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultT0CCutEnabled, bool, false, "Enable Global multiplicity vs T0C centrality cut")
    Configurable<std::vector<double>> cfgMultT0CCutPars{"cfgMultT0CCutPars", std::vector<double>{143.04, -4.58368, 0.0766055, -0.000727796, 2.86153e-06, 23.3108, -0.36304, 0.00437706, -4.717e-05, 1.98332e-07}, "Global multiplicity vs T0C centrality cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgMultPVT0CCutEnabled, bool, false, "Enable PV multiplicity vs T0C centrality cut")
    Configurable<std::vector<double>> cfgMultPVT0CCutPars{"cfgMultPVT0CCutPars", std::vector<double>{195.357, -6.15194, 0.101313, -0.000955828, 3.74793e-06, 30.0326, -0.43322, 0.00476265, -5.11206e-05, 2.13613e-07}, "PV multiplicity vs T0C centrality cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgMultMultPVHighCutFunction, std::string, "[0]+[1]*x + 5.*([2]+[3]*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultMultPVLowCutFunction, std::string, "[0]+[1]*x - 5.*([2]+[3]*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultGlobalPVCutEnabled, bool, false, "Enable global multiplicity vs PV multiplicity cut")
    Configurable<std::vector<double>> cfgMultGlobalPVCutPars{"cfgMultGlobalPVCutPars", std::vector<double>{-0.140809, 0.734344, 2.77495, 0.0165935}, "PV multiplicity vs T0C centrality cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgMultMultV0AHighCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 4.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultMultV0ALowCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultMultV0ACutEnabled, bool, false, "Enable global multiplicity vs V0A multiplicity cut")
    Configurable<std::vector<double>> cfgMultMultV0ACutPars{"cfgMultMultV0ACutPars", std::vector<double>{534.893, 184.344, 0.423539, -0.00331436, 5.34622e-06, 871.239, 53.3735, -0.203528, 0.000122758, 5.41027e-07}, "Global multiplicity vs V0A multiplicity cut parameter values"};
    std::vector<double> multT0CCutPars;
    std::vector<double> multPVT0CCutPars;
    std::vector<double> multGlobalPVCutPars;
    std::vector<double> multMultV0ACutPars;
    std::unique_ptr<TF1> fMultPVT0CCutLow = nullptr;
    std::unique_ptr<TF1> fMultPVT0CCutHigh = nullptr;
    std::unique_ptr<TF1> fMultT0CCutLow = nullptr;
    std::unique_ptr<TF1> fMultT0CCutHigh = nullptr;
    std::unique_ptr<TF1> fMultGlobalPVCutLow = nullptr;
    std::unique_ptr<TF1> fMultGlobalPVCutHigh = nullptr;
    std::unique_ptr<TF1> fMultMultV0ACutLow = nullptr;
    std::unique_ptr<TF1> fMultMultV0ACutHigh = nullptr;
    std::unique_ptr<TF1> fT0AV0AMean = nullptr;
    std::unique_ptr<TF1> fT0AV0ASigma = nullptr;
    std::unique_ptr<TF1> fPtDepDCAxy = nullptr;
    std::unique_ptr<TF1> fPtDepDCAxyForNch = nullptr;
    std::unique_ptr<TF1> fPtDepDCAz = nullptr;
    std::unique_ptr<TF1> fPtDepDCAzForNch = nullptr;
    O2_DEFINE_CONFIGURABLE(cfgV0AT0Acut, int, 5, "V0AT0A cut")
  } cfgFuncParas;

  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgTofPtCut, float, 0.4f, "Minimum pt to use TOF N-sigma")
    O2_DEFINE_CONFIGURABLE(cfgUseItsPID, bool, false, "Use ITS PID for particle identification")
    O2_DEFINE_CONFIGURABLE(cfgQANsigma, bool, false, "Get QA histograms for selection of pions, kaons, and protons")
    O2_DEFINE_CONFIGURABLE(cfgQAdEdx, bool, false, "Get dEdx histograms for pions, kaons, and protons")
    O2_DEFINE_CONFIGURABLE(cfgPIDUseRejection, bool, true, "True: use exclusion exclusion criteria for PID determination, false: don't use exclusion")
    Configurable<LabeledArray<float>> nSigmas{"nSigmas", {LongArrayFloat.front().data(), 6, 3, {"UpCut_pi", "UpCut_ka", "UpCut_pr", "LowCut_pi", "LowCut_ka", "LowCut_pr"}, {"TPC", "TOF", "ITS"}}, "Labeled array for n-sigma values for TPC, TOF, ITS for pions, kaons, protons (positive and negative)"};
  } cfgPIDConfigs;

  SliceCache cache;

  // Axes for general QA
  ConfigurableAxis axisVertex{"axisVertex", {10, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisMult{"axisMult", {10, 0, 100}, "multiplicity axis for QA histograms"};
  ConfigurableAxis axisCent{"axisCent", {100, 0, 100}, "centrality axis for QA histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {72, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisPtFiner{"axisPtFiner", {98, 0.2, 10.0}, "pt axis for histograms"};
  ConfigurableAxis axisDCAz{"axisDCAz", {200, -2, 2}, "DCA_{z} (cm)"};
  ConfigurableAxis axisDCAxy{"axisDCAxy", {200, -1, 1}, "DCA_{xy} (cm)"};

  // Axes for Correlations
  ConfigurableAxis axisMultForCorr{"axisMultForCorr", {VARIABLE_WIDTH, 0, 5, 10, 20, 50, 100, 150, 200, 300, 500}, "multiplicity axis for correlation histograms"};
  ConfigurableAxis axisCentForCorr{"axisCentForCorr", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, "centrality axis for correlation histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEtaTpcFt0a{"axisDeltaEtaTpcFt0a", {32, -5.8, -2.6}, "delta eta axis, -5.8~-2.6 for TPC-FT0A,"};
  ConfigurableAxis axisDeltaEtaTpcFt0c{"axisDeltaEtaTpcFt0c", {32, 1.2, 4.2}, "delta eta axis, 1.2~4.2 for TPC-FT0C"};
  ConfigurableAxis axisDeltaEtaFt0aFt0c{"axisDeltaEtaFt0aFt0c", {32, -1.5, 3.0}, "delta eta axis"};
  ConfigurableAxis axisVtxMix{"axisVtxMix", {VARIABLE_WIDTH, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, "vertex axis for mixed event histograms"};
  ConfigurableAxis axisMultMix{"axisMultMix", {VARIABLE_WIDTH, 0, 10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260}, "multiplicity / centrality axis for mixed event histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt trigger axis for histograms"};
  ConfigurableAxis axisParticle{"axisParticle", {4, 0, 4}, "particle axis for correlation container"};
  ConfigurableAxis axisTpcNcls{"axisTpcNcls", {500, 0, 500}, "number of TPC clusters used for PID"};

  // Axes for long-range correlations QA
  ConfigurableAxis axisAmplitudeFt0{"axisAmplitudeFt0a", {5000, 0, 1000}, "FT0A amplitude"};
  ConfigurableAxis axisChannelFt0aAxis{"axisChannelFt0aAxis", {96, 0.0, 96.0}, "FT0A channel"};
  ConfigurableAxis axisChannelFt0cAxis{"axisChannelFt0cAxis", {96, 0.0, 96.0}, "FT0C channel"};
  ConfigurableAxis cfgFITamp{"cfgFITamp", {1000, 0, 5000}, "FT0 amplitude"};

  // Axes for PID QA
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {80, -5, 5}, "nsigmaTPC axis"};
  ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {80, -5, 5}, "nsigmaTOF axis"};
  ConfigurableAxis axisNsigmaITS{"axisNsigmaITS", {80, -5, 5}, "nsigmaITS axis"};
  ConfigurableAxis axisTpcSignal{"axisTpcSignal", {250, 0, 250}, "dEdx axis for TPC"};
  ConfigurableAxis axisSigma{"axisSigma", {200, 0, 20}, "sigma axis for TPC"};

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgEvSel.cfgGlobalCutVertex);
  Filter trackFilter = (nabs(aod::track::eta) < cfgTrkSel.cfgGlobalEta) && (aod::track::pt > cfgTrkSel.cfgGlobalPtMin) && (aod::track::pt < cfgTrkSel.cfgGlobalPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t)true)) && (nabs(aod::track::dcaZ) < cfgTrkSel.cfgGlobalDCAz) && (aod::track::tpcChi2NCl < cfgTrkSel.cfgGlobalChi2prTPCcls);

  using FilteredCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSel, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::Mults>>;
  using FilteredTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;

  // FT0 geometry
  o2::ft0::Geometry ft0Det;
  static constexpr uint64_t Ft0IndexA = 96;
  std::vector<o2::detectors::AlignParam>* offsetFT0 = nullptr;
  std::vector<float> cstFT0RelGain;

  // Corrections
  TH3D* mEfficiency = nullptr;
  TH1D* mEfficiencyNch = nullptr;
  TH1D* mCentralityWeight = nullptr;
  bool correctionsLoaded = false;

  // Define the outputs
  OutputObj<DihadronContainer> sameTpcFt0a{"sameEvent_TPC_FT0A"};
  OutputObj<DihadronContainer> mixedTpcFt0a{"mixedEvent_TPC_FT0A"};
  OutputObj<DihadronContainer> sameTpcFt0c{"sameEvent_TPC_FT0C"};
  OutputObj<DihadronContainer> mixedTpcFt0c{"mixedEvent_TPC_FT0C"};
  OutputObj<DihadronContainer> sameFt0aFt0c{"sameEvent_FT0A_FT0C"};
  OutputObj<DihadronContainer> mixedFt0aFt0c{"mixedEvent_FT0A_FT0C"};

  HistogramRegistry histos{"histos"};

  // For Labelled array value containers
  std::array<std::array<float, 2>, 10> trackCuts{};
  std::array<std::array<int, 1>, 14> eventCuts{};
  std::array<std::array<float, 3>, 6> nSigmaVals{};

  // define global variables
  TRandom3 fRandom{0};

  enum EventType {
    SameEvent = 1,
    MixedEvent = 3
  };

  enum FITIndex {
    IndexFT0A = 0,
    IndexFT0C = 1
  };

  enum PIDIndex {
    UseCharged = 0,
    UsePions,
    UseKaons,
    UseProtons,
    UseK0,
    UseLambda,
    UsePhi
  };
  enum PiKpArrayIndex {
    IndexPionUp = 0,
    IndexKaonUp,
    IndexProtonUp,
    IndexPionLow,
    IndexKaonLow,
    IndexProtonLow
  };
  enum DetectorType {
    UseTPC = 0,
    UseTof,
    UseITS
  };
  enum Stage {
    Before = 0,
    After
  };

  enum EventCutTypes {
    FilteredEvents = 0,
    AfterSel8,
    UseNoTimeFrameBorder,
    UseNoITSROFrameBorder,
    UseNoSameBunchPileup,
    UseGoodZvtxFT0vsPV,
    UseNoCollInTimeRangeStandard,
    UseGoodITSLayersAll,
    UseGoodITSLayer0123,
    UseNoCollInRofStandard,
    UseNoHighMultCollInPrevRof,
    UseOccupancy,
    UseMultCorrCut,
    UseT0AV0ACut,
    HaveFT0Cut,
    NumEventCuts
  };

  enum EventCutType {
    EvCut1 = 0,
    NumEvCutTypes = 1
  };

  enum TrackCuts {
    TrkCutPtMin = 0,
    TrkCutPtMax,
    TrkCutEtaCut,
    TrkCutChi2PrTpcCls,
    TrkCutTpcCluster,
    TrkCutTpcCrossedRows,
    TrkCutItsCluster,
    TrkCutDCAz,
    TrkCutDCAzNsigma,
    TrkCutDCAxyNsigma
  };

  enum TrackCutGroup {
    UseGenTrkCuts = 0,
    UseNchSelCuts,
    NumTrackCutTypes
  };

  enum CentEstimators {
    UseCentFT0C = 0,
    UseCentFT0CVariant1,
    UseCentFT0M,
    UseCentFV0A,
    // Count the total number of enum
    NumCentEstimators
  };

  RCTFlagsChecker rctChecker{"CBT"};

  void init(InitContext&)
  {
    // ----------------------------------------------------------------------------------------------------------------------
    // The way the code is initialized, TPC-FT0 correlations are done together and FT0A-FT0C correlations are done separately
    // ----------------------------------------------------------------------------------------------------------------------

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    LOGF(info, "Starting init");

    readMatrix(cfgTrkSel.cfgTrackCuts->getData(), trackCuts);
    readMatrix(cfgEvSel.cfgUseEventCuts->getData(), eventCuts);
    readMatrix(cfgPIDConfigs.nSigmas->getData(), nSigmaVals);

    if (doprocessSameFt0aFt0c || doprocessSameTpcFt0 || doprocessQA) {
      histos.add("hEventCountRct", "Number of Event;; Count", {HistType::kTH1D, {{2, 0, 2}}});
      histos.get<TH1>(HIST("hEventCountRct"))->GetXaxis()->SetBinLabel(1, "rct fail");
      histos.get<TH1>(HIST("hEventCountRct"))->GetXaxis()->SetBinLabel(2, "rct pass");
      histos.add("hEventCount", "Number of Event;; Count", {HistType::kTH1D, {{15, -0.5, 14.5}}});
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(FilteredEvents + 1, "Filtered events");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(AfterSel8 + 1, "After sel8");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseNoTimeFrameBorder + 1, "kNoTimeFrameBorder");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseNoITSROFrameBorder + 1, "kNoITSROFrameBorder");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseNoSameBunchPileup + 1, "kNoSameBunchPileup");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseGoodZvtxFT0vsPV + 1, "kIsGoodZvtxFT0vsPV");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseNoCollInTimeRangeStandard + 1, "kNoCollInTimeRangeStandard");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseGoodITSLayersAll + 1, "kIsGoodITSLayersAll");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseGoodITSLayer0123 + 1, "kIsGoodITSLayer0123");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseNoCollInRofStandard + 1, "kNoCollInRofStandard");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseNoHighMultCollInPrevRof + 1, "kNoHighMultCollInPrevRof");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseOccupancy + 1, "Occupancy Cut");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseMultCorrCut + 1, "MultCorrelation Cut");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseT0AV0ACut + 1, "T0AV0A cut");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(HaveFT0Cut + 1, "Event has FT0 cut");
    }

    if ((doprocessSameFt0aFt0c || doprocessSameTpcFt0 || doprocessQA) && cfgGeneral.cfgQABasic) {
      histos.add("hPassedEventSelection", "Number of Event;; Count", {HistType::kTH1D, {{12, -0.5, 11.5}}});
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(FilteredEvents + 1, "Filtered events");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(AfterSel8 + 1, "After sel8");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(UseNoTimeFrameBorder + 1, "kNoTimeFrameBorder");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(UseNoITSROFrameBorder + 1, "kNoITSROFrameBorder");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(UseNoSameBunchPileup + 1, "kNoSameBunchPileup");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(UseGoodZvtxFT0vsPV + 1, "kIsGoodZvtxFT0vsPV");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(UseNoCollInTimeRangeStandard + 1, "kNoCollInTimeRangeStandard");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(UseGoodITSLayersAll + 1, "kIsGoodITSLayersAll");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(UseGoodITSLayer0123 + 1, "kIsGoodITSLayer0123");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(UseNoCollInRofStandard + 1, "kNoCollInRofStandard");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(UseNoHighMultCollInPrevRof + 1, "kNoHighMultCollInPrevRof");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(UseOccupancy + 1, "Occupancy Cut");
    }

    // Multiplicity correlation cuts
    if (eventCuts[UseMultCorrCut][EvCut1] != 0) {
      cfgFuncParas.multT0CCutPars = cfgFuncParas.cfgMultT0CCutPars;
      cfgFuncParas.multPVT0CCutPars = cfgFuncParas.cfgMultPVT0CCutPars;
      cfgFuncParas.multGlobalPVCutPars = cfgFuncParas.cfgMultGlobalPVCutPars;
      cfgFuncParas.multMultV0ACutPars = cfgFuncParas.cfgMultMultV0ACutPars;

      cfgFuncParas.fMultPVT0CCutLow = std::make_unique<TF1>("fMultPVT0CCutLow", cfgFuncParas.cfgMultCentLowCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultPVT0CCutLow->SetParameters(cfgFuncParas.multPVT0CCutPars.data());
      cfgFuncParas.fMultPVT0CCutHigh = std::make_unique<TF1>("fMultPVT0CCutHigh", cfgFuncParas.cfgMultCentHighCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultPVT0CCutHigh->SetParameters(cfgFuncParas.multPVT0CCutPars.data());

      cfgFuncParas.fMultT0CCutLow = std::make_unique<TF1>("fMultT0CCutLow", cfgFuncParas.cfgMultCentLowCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultT0CCutLow->SetParameters(cfgFuncParas.multT0CCutPars.data());
      cfgFuncParas.fMultT0CCutHigh = std::make_unique<TF1>("fMultT0CCutHigh", cfgFuncParas.cfgMultCentHighCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultT0CCutHigh->SetParameters(cfgFuncParas.multT0CCutPars.data());

      cfgFuncParas.fMultGlobalPVCutLow = std::make_unique<TF1>("fMultGlobalPVCutLow", cfgFuncParas.cfgMultMultPVLowCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultGlobalPVCutLow->SetParameters(cfgFuncParas.multGlobalPVCutPars.data());
      cfgFuncParas.fMultGlobalPVCutHigh = std::make_unique<TF1>("fMultGlobalPVCutHigh", cfgFuncParas.cfgMultMultPVHighCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultGlobalPVCutHigh->SetParameters(cfgFuncParas.multGlobalPVCutPars.data());

      cfgFuncParas.fMultMultV0ACutLow = std::make_unique<TF1>("fMultMultV0ACutLow", cfgFuncParas.cfgMultMultV0ALowCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultMultV0ACutLow->SetParameters(cfgFuncParas.multMultV0ACutPars.data());
      cfgFuncParas.fMultMultV0ACutHigh = std::make_unique<TF1>("fMultMultV0ACutHigh", cfgFuncParas.cfgMultMultV0AHighCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultMultV0ACutHigh->SetParameters(cfgFuncParas.multMultV0ACutPars.data());
    }
    if (eventCuts[UseT0AV0ACut][EvCut1] != 0) {
      cfgFuncParas.fT0AV0AMean = std::make_unique<TF1>("fT0AV0AMean", "[0]+[1]*x", 0, 200000);
      cfgFuncParas.fT0AV0AMean->SetParameters(-1601.0581, 9.417652e-01);
      cfgFuncParas.fT0AV0ASigma = std::make_unique<TF1>("fT0AV0ASigma", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 200000);
      cfgFuncParas.fT0AV0ASigma->SetParameters(463.4144, 6.796509e-02, -9.097136e-07, 7.971088e-12, -2.600581e-17);
    }

    if (!cfgTrkSel.cfgTrackCutsDCAxy.value[UseGenTrkCuts].empty()) {
      cfgFuncParas.fPtDepDCAxy = std::make_unique<TF1>("ptDepDCAxy", cfgTrkSel.cfgTrackCutsDCAxy.value[UseGenTrkCuts].c_str(), 0.001, 1000);
      cfgFuncParas.fPtDepDCAxy->SetParameter(0, trackCuts[TrkCutDCAxyNsigma][UseGenTrkCuts]);
      LOGF(info, "DCAxy pt-dependence function: %0.1f * %s", trackCuts[TrkCutDCAxyNsigma][UseGenTrkCuts], cfgTrkSel.cfgTrackCutsDCAxy.value[UseGenTrkCuts].c_str());
    }
    if (!cfgTrkSel.cfgTrackCutsDCAxy.value[UseNchSelCuts].empty()) {
      cfgFuncParas.fPtDepDCAxyForNch = std::make_unique<TF1>("ptDepDCAxyForNch", cfgTrkSel.cfgTrackCutsDCAxy.value[UseNchSelCuts].c_str(), 0.001, 1000);
      cfgFuncParas.fPtDepDCAxyForNch->SetParameter(0, trackCuts[TrkCutDCAxyNsigma][UseNchSelCuts]);
      LOGF(info, "DCAxy pt-dependence function for Nch: %0.1f * %s", trackCuts[TrkCutDCAxyNsigma][UseNchSelCuts], cfgTrkSel.cfgTrackCutsDCAxy.value[UseNchSelCuts].c_str());
    }

    if (!cfgTrkSel.cfgTrackCutsDCAz.value[UseGenTrkCuts].empty()) {
      cfgFuncParas.fPtDepDCAz = std::make_unique<TF1>("ptDepDCAz", cfgTrkSel.cfgTrackCutsDCAz.value[UseGenTrkCuts].c_str(), 0.001, 1000);
      cfgFuncParas.fPtDepDCAz->SetParameter(0, trackCuts[TrkCutDCAzNsigma][UseGenTrkCuts]);
      LOGF(info, "DCAz pt-dependence function: %0.1f * %s", trackCuts[TrkCutDCAzNsigma][UseGenTrkCuts], cfgTrkSel.cfgTrackCutsDCAz.value[UseGenTrkCuts].c_str());
    }
    if (!cfgTrkSel.cfgTrackCutsDCAz.value[UseNchSelCuts].empty()) {
      cfgFuncParas.fPtDepDCAzForNch = std::make_unique<TF1>("ptDepDCAzForNch", cfgTrkSel.cfgTrackCutsDCAz.value[UseNchSelCuts].c_str(), 0.001, 1000);
      cfgFuncParas.fPtDepDCAzForNch->SetParameter(0, trackCuts[TrkCutDCAzNsigma][UseNchSelCuts]);
      LOGF(info, "DCAz pt-dependence function for Nch: %0.1f * %s", trackCuts[TrkCutDCAzNsigma][UseNchSelCuts], cfgTrkSel.cfgTrackCutsDCAz.value[UseNchSelCuts].c_str());
    }

    const AxisSpec axisT0C{70, 0, 70000, "N_{ch} (T0C)"};
    const AxisSpec axisT0A{200, 0, 200000, "N_{ch} (T0A)"};
    const AxisSpec axisChi2{100, 0., 10.};
    const AxisSpec axisChID = {220, 0, 220};

    auto maxSample = static_cast<double>(cfgGeneral.cfgSampleSize);
    AxisSpec axisSample{cfgGeneral.cfgSampleSize, 0, maxSample, "Sample"};

    // Choose if it is Nch selection or Centrality selection
    AxisSpec axisEventClass = {axisMultForCorr, "N_{ch}"};

    if (!cfgEvSel.cfgSelCollByNch) {
      axisEventClass = {axisCentForCorr, "Centrality"};
    }

    if (doprocessQA) {
      histos.add("h_globalTracks_centT0C_before", "before cut;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {axisCent, axisMult}});
      histos.add("h_PVTracks_centT0C_before", "before cut;Centrality T0C;mulplicity PV tracks", {HistType::kTH2D, {axisCent, axisMult}});
      histos.add("h_globalTracks_PVTracks_before", "before cut;mulplicity PV tracks;mulplicity global tracks", {HistType::kTH2D, {axisMult, axisMult}});
      histos.add("h_globalTracks_multT0A_before", "before cut;mulplicity T0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisMult}});
      histos.add("h_globalTracks_multV0A_before", "before cut;mulplicity V0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisMult}});
      histos.add("h_multV0A_multT0A_before", "before cut;mulplicity T0A;mulplicity V0A", {HistType::kTH2D, {axisT0A, axisT0A}});
      histos.add("h_multT0C_centT0C_before", "before cut;Centrality T0C;mulplicity T0C", {HistType::kTH2D, {axisCent, axisT0C}});

      histos.add("h_globalTracks_centT0C_after", "after cut;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {axisCent, axisMult}});
      histos.add("h_PVTracks_centT0C_after", "after cut;Centrality T0C;mulplicity PV tracks", {HistType::kTH2D, {axisCent, axisMult}});
      histos.add("h_globalTracks_PVTracks_after", "after cut;mulplicity PV tracks;mulplicity global tracks", {HistType::kTH2D, {axisMult, axisMult}});
      histos.add("h_globalTracks_multT0A_after", "after cut;mulplicity T0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisMult}});
      histos.add("h_globalTracks_multV0A_after", "after cut;mulplicity V0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisMult}});
      histos.add("h_multV0A_multT0A_after", "after cut;mulplicity T0A;mulplicity V0A", {HistType::kTH2D, {axisT0A, axisT0A}});
      histos.add("h_multT0C_centT0C_after", "after cut;Centrality T0C;mulplicity T0C", {HistType::kTH2D, {axisCent, axisT0C}});
      histos.add("h_centFT0M_centFT0C", "after cut;Centrality T0C;Centrality T0M", {HistType::kTH2D, {axisCent, axisCent}});
      histos.add("h_centFV0A_centFT0C", "after cut;Centrality T0C;Centrality V0A", {HistType::kTH2D, {axisCent, axisCent}});

      // Add track cuts table
      histos.add("hDCAz_before", "DCAz before cuts; DCAz (cm); Pt", {HistType::kTH2D, {axisDCAz, axisPtFiner}});
      histos.add("hDCAxy_before", "DCAxy before cuts; DCAxy (cm); Pt", {HistType::kTH2D, {axisDCAxy, axisPtFiner}});
      histos.add("hTPCNclsPID_before", "hTPCNclsPID_before", {HistType::kTH1D, {axisTpcNcls}});
      histos.add("hTPCNclsFound_before", "hTPCNclsFound_before", {HistType::kTH1D, {axisTpcNcls}});
      histos.add("hTPCCrossedRows_before", "hTPCCrossedRows_before", {HistType::kTH1D, {axisTpcNcls}});

      histos.add("hDCAz_after", "DCAz after cuts; DCAz (cm); Pt", {HistType::kTH3D, {axisDCAz, axisPtFiner, axisParticle}});
      histos.add("hDCAxy_after", "DCAxy after cuts; DCAxy (cm); Pt", {HistType::kTH3D, {axisDCAxy, axisPtFiner, axisParticle}});
      histos.add("hTPCNclsPID_after", "hTPCNclsPID_after", {HistType::kTH1D, {axisTpcNcls}});
      histos.add("hTPCNclsFound_after", "hTPCNclsFound_after", {HistType::kTH1D, {axisTpcNcls}});
      histos.add("hTPCCrossedRows_after", "hTPCCrossedRows_after", {HistType::kTH1D, {axisTpcNcls}});

      histos.add("hChi2prTPCcls", "#chi^{2}/cluster for the TPC track segment", {HistType::kTH1D, {axisChi2}});
      histos.add("hChi2prITScls", "#chi^{2}/cluster for the ITS track", {HistType::kTH1D, {axisChi2}});
      histos.add("hITSNclsFound", "Number of found ITS clusters", {HistType::kTH1D, {axisTpcNcls}});
    }

    if (doprocessSameFt0aFt0c || doprocessSameTpcFt0 || doprocessQA) {
      histos.add("zVtx", "zVtx", {HistType::kTH1D, {axisVertex}});
      histos.add("Nch", "N_{ch}", {HistType::kTH1D, {axisMult}});
      histos.add("Nch_corrected", "N_{ch} corrected", {HistType::kTH1D, {axisMult}});
      histos.add("Centrality", "Centrality", {HistType::kTH1D, {axisCent}});
      histos.add("CentralityWeighted", "Centrality (weighted)", {HistType::kTH1D, {axisCent}});
      histos.add("hTrackCorrection2d", "Correlation table for number of tracks table; uncorrected track; corrected track", {HistType::kTH2D, {axisMult, axisMult}});
    }

    if (doprocessSameTpcFt0 || doprocessQA) {
      if (cfgGeneral.cfgQABasic) {
        histos.add("Phi", "Phi", {HistType::kTH2D, {axisPhi, axisParticle}});
        histos.add("Eta", "Eta", {HistType::kTH2D, {axisEta, axisParticle}});
        histos.add("EtaCorrected", "EtaCorrected", {HistType::kTH2D, {axisEta, axisParticle}});
        histos.add("pTFiner", "pTFiner", {HistType::kTHnSparseF, {axisPtFiner, axisParticle, axisEventClass}});
        histos.add("pTFinerCorrected", "pTFinerCorrected", {HistType::kTHnSparseF, {axisPtFiner, axisParticle, axisEventClass}});
      }
      // PID nSigma histograms
      if (cfgPIDConfigs.cfgQANsigma) {
        if (!cfgPIDConfigs.cfgUseItsPID) {
          histos.add("TofTpcNsigma_beforeCut", "", {HistType::kTHnSparseD, {{axisNsigmaTPC, axisNsigmaTOF, axisPtTrigger, axisParticle}}});
          histos.add("TofTpcNsigma_afterCut", "", {HistType::kTHnSparseD, {{axisNsigmaTPC, axisNsigmaTOF, axisPtTrigger, axisParticle}}});
        } // TPC-TOF PID QA hists
        if (cfgPIDConfigs.cfgUseItsPID) {
          histos.add("TofItsNsigma_beforeCut", "", {HistType::kTHnSparseD, {{axisNsigmaITS, axisNsigmaTOF, axisPtTrigger, axisParticle}}});
          histos.add("TofItsNsigma_afterCut", "", {HistType::kTHnSparseD, {{axisNsigmaITS, axisNsigmaTOF, axisPtTrigger, axisParticle}}});
        } // ITS-TOF PID QA hists
      } // end of PID QA hists

      // PID dEdx histograms
      if (cfgPIDConfigs.cfgQAdEdx) {
        histos.add("TpcdEdx_ptwise_beforeCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF, axisParticle}}});
        histos.add("ExpTpcdEdx_ptwise_beforeCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF, axisParticle}}});
        histos.add("ExpSigma_ptwise_beforeCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisSigma, axisNsigmaTOF, axisParticle}}});

        histos.add("TpcdEdx_ptwise_afterCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF, axisParticle}}});
        histos.add("ExpTpcdEdx_ptwise_afterCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF, axisParticle}}});
        histos.add("ExpSigma_ptwise_afterCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisSigma, axisNsigmaTOF, axisParticle}}});
      }
    }
    if (doprocessSameTpcFt0) { // QA plots are based on TPC tracks, so they are only included in TPC-FT0A process and not in TPC-FT0C process
      // TPC-FT0A correlation histograms
      if (cfgGeneral.cfgAnalyzeTPCFT0A) {
        histos.add("deltaEta_deltaPhi_same_TPC_FT0A", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0a}}); // check to see the delta eta and delta phi distribution
        histos.add("deltaEta_deltaPhi_mixed_TPC_FT0A", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0a}});
        histos.add("Assoc_amp_same_TPC_FT0A", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0}});
        histos.add("Assoc_amp_mixed_TPC_FT0A", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0}});
        histos.add("Trig_hist_TPC_FT0A", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisEventClass, axisPtTrigger, axisParticle}}});
      }
      // TPC-FT0C correlation histograms
      if (cfgGeneral.cfgAnalyzeTPCFT0C) {
        histos.add("deltaEta_deltaPhi_same_TPC_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0c}}); // check to see the delta eta and delta phi distribution
        histos.add("deltaEta_deltaPhi_mixed_TPC_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0c}});
        histos.add("Assoc_amp_same_TPC_FT0C", "", {HistType::kTH2D, {axisChannelFt0cAxis, axisAmplitudeFt0}});
        histos.add("Assoc_amp_mixed_TPC_FT0C", "", {HistType::kTH2D, {axisChannelFt0cAxis, axisAmplitudeFt0}});
        histos.add("Trig_hist_TPC_FT0C", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisEventClass, axisPtTrigger, axisParticle}}});
      }

      histos.add("FT0Amp", "", {HistType::kTH2F, {axisChID, cfgFITamp}});
      histos.add("FT0AmpCorrect", "", {HistType::kTH2F, {axisChID, cfgFITamp}});
      histos.add("eventcount", "bin", {HistType::kTH1F, {{4, 0, 4, "bin"}}}); // histogram to see how many events are in the same and mixed event
    }

    if (doprocessSameFt0aFt0c) {
      histos.add("Phi_FT0A_FT0C", "Phi_FT0A_FT0C", {HistType::kTH1D, {axisPhi}});
      histos.add("Eta_FT0A_FT0C", "Eta_FT0A_FT0C", {HistType::kTH1D, {axisEta}});
      histos.add("EtaCorrected_FT0A_FT0C", "EtaCorrected_FT0A_FT0C", {HistType::kTH1D, {axisEta}});
      histos.add("pTFiner_FT0A_FT0C", "pTFiner_FT0A_FT0C", {HistType::kTH1D, {axisPtFiner}});
      histos.add("pTFinerCorrected_FT0A_FT0C", "pTFinerCorrected_FT0A_FT0C", {HistType::kTH1D, {axisPtFiner}});
      histos.add("deltaEta_deltaPhi_same_FT0A_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaFt0aFt0c}}); // check to see the delta eta and delta phi distribution
      histos.add("deltaEta_deltaPhi_mixed_FT0A_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaFt0aFt0c}});
      histos.add("Trig_hist_FT0A_FT0C", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisEventClass}}});

      histos.add("FT0Amp", "", {HistType::kTH2F, {axisChID, cfgFITamp}});
      histos.add("FT0AmpCorrect", "", {HistType::kTH2F, {axisChID, cfgFITamp}});
      histos.add("eventcount", "bin", {HistType::kTH1F, {{4, 0, 4, "bin"}}}); // histogram to see how many events are in the same and mixed event
    }

    LOGF(info, "Initializing correlation container");

    // Initialize Nch-related histograms and containers
    std::vector<AxisSpec> corrAxisTpcFt0a = {{axisSample},
                                             {axisVertex, "z-vtx (cm)"},
                                             {axisEventClass},
                                             {axisDeltaPhi, "#Delta#varphi (rad)"},
                                             {axisDeltaEtaTpcFt0a, "#Delta#eta"},
                                             {axisPtTrigger, "p_{T} (GeV/c)"},
                                             {axisParticle, "Particle, 0 = charged, 1 = pion, 2 = kaon, 3 = proton"}};

    std::vector<AxisSpec> corrAxisTpcFt0c = {{axisSample},
                                             {axisVertex, "z-vtx (cm)"},
                                             {axisEventClass},
                                             {axisDeltaPhi, "#Delta#varphi (rad)"},
                                             {axisDeltaEtaTpcFt0c, "#Delta#eta"},
                                             {axisPtTrigger, "p_{T} (GeV/c)"},
                                             {axisParticle, "Particle, 0 = charged, 1 = pion, 2 = kaon, 3 = proton"}};

    std::vector<AxisSpec> corrAxisFt0aFt0c = {{axisSample},
                                              {axisVertex, "z-vtx (cm)"},
                                              {axisEventClass},
                                              {axisDeltaPhi, "#Delta#varphi (rad)"},
                                              {axisDeltaEtaFt0aFt0c, "#Delta#eta"}};

    if (doprocessSameTpcFt0) {
      if (cfgGeneral.cfgAnalyzeTPCFT0A) {
        sameTpcFt0a.setObject(new DihadronContainer("sameEvent_TPC_FT0A", "sameEvent_TPC_FT0A", corrAxisTpcFt0a));
        mixedTpcFt0a.setObject(new DihadronContainer("mixedEvent_TPC_FT0A", "mixedEvent_TPC_FT0A", corrAxisTpcFt0a));
      }
      if (cfgGeneral.cfgAnalyzeTPCFT0C) {
        sameTpcFt0c.setObject(new DihadronContainer("sameEvent_TPC_FT0C", "sameEvent_TPC_FT0C", corrAxisTpcFt0c));
        mixedTpcFt0c.setObject(new DihadronContainer("mixedEvent_TPC_FT0C", "mixedEvent_TPC_FT0C", corrAxisTpcFt0c));
      }
    }

    if (doprocessSameFt0aFt0c) {
      sameFt0aFt0c.setObject(new DihadronContainer("sameEvent_FT0A_FT0C", "sameEvent_FT0A_FT0C", corrAxisFt0aFt0c));
      mixedFt0aFt0c.setObject(new DihadronContainer("mixedEvent_FT0A_FT0C", "mixedEvent_FT0A_FT0C", corrAxisFt0aFt0c));
    }

    if (!cfgEvSel.cfgEvSelRCTflags.value.empty()) {
      rctChecker.init(cfgEvSel.cfgEvSelRCTflags.value); // override initialzation
    }

    LOGF(info, "End of init");
  }

  template <typename TCollision>
  double getCentrality(TCollision const& collision)
  {
    double cent = 0.0;
    switch (cfgTrkSel.cfgCentEstimator) {
      case UseCentFT0C:
        cent = collision.centFT0C();
        break;
      case UseCentFT0CVariant1:
        cent = collision.centFT0CVariant1();
        break;
      case UseCentFT0M:
        cent = collision.centFT0M();
        break;
      case UseCentFV0A:
        cent = collision.centFV0A();
        break;
      default:
        cent = collision.centFT0C();
    }
    return cent;
  }

  template <typename TCollision>
  bool eventRct(TCollision const& collision, const bool fillCounter)
  {
    if (!rctChecker(collision)) {
      if (fillCounter) {
        histos.fill(HIST("hEventCountRct"), 0.5);
      }

      return false;
    }
    if (fillCounter) {
      histos.fill(HIST("hEventCountRct"), 1.5);
    }

    return true;
  }

  template <typename TCollision>
  bool eventSelected(TCollision const& collision, const int mult, const double cent, const bool fillCounter)
  {
    if (fillCounter) {
      histos.fill(HIST("hEventCount"), FilteredEvents);
    }
    if (!collision.sel8()) {
      return false;
    }
    if (fillCounter) {
      histos.fill(HIST("hEventCount"), AfterSel8);
    }

    if (eventCuts[UseNoTimeFrameBorder][EvCut1] && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (fillCounter && eventCuts[UseNoTimeFrameBorder][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseNoTimeFrameBorder);
    }

    if (eventCuts[UseNoITSROFrameBorder][EvCut1] && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (fillCounter && eventCuts[UseNoITSROFrameBorder][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseNoITSROFrameBorder);
    }

    if (eventCuts[UseNoSameBunchPileup][EvCut1] && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return false;
    }
    if (fillCounter && eventCuts[UseNoSameBunchPileup][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseNoSameBunchPileup);
    }

    if (eventCuts[UseGoodZvtxFT0vsPV][EvCut1] && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return false;
    }
    if (fillCounter && eventCuts[UseGoodZvtxFT0vsPV][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseGoodZvtxFT0vsPV);
    }

    if (eventCuts[UseNoCollInTimeRangeStandard][EvCut1] && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      // no collisions in specified time range
      return false;
    }

    if (fillCounter && eventCuts[UseNoCollInTimeRangeStandard][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseNoCollInTimeRangeStandard);
    }

    if (eventCuts[UseGoodITSLayersAll][EvCut1] && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      // from Jan 9 2025 AOT meeting
      // cut time intervals with dead ITS staves
      return false;
    }

    if (fillCounter && eventCuts[UseGoodITSLayersAll][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseGoodITSLayersAll);
    }

    if (eventCuts[UseGoodITSLayer0123][EvCut1] && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
      return false;
    }
    if (fillCounter && eventCuts[UseGoodITSLayer0123][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseGoodITSLayer0123);
    }

    if (eventCuts[UseNoCollInRofStandard][EvCut1] && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      // no other collisions in this Readout Frame with per-collision multiplicity above threshold
      return false;
    }

    if (fillCounter && eventCuts[UseNoCollInRofStandard][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseNoCollInRofStandard);
    }

    if (eventCuts[UseNoHighMultCollInPrevRof][EvCut1] && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      // veto an event if FT0C amplitude in previous ITS ROF is above threshold
      return false;
    }
    if (fillCounter && eventCuts[UseNoHighMultCollInPrevRof][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseNoHighMultCollInPrevRof);
    }

    auto multNTracksPV = collision.multNTracksPV();
    auto occupancy = collision.trackOccupancyInTimeRange();

    if (eventCuts[UseOccupancy][EvCut1] && (occupancy < cfgEvSel.cfgCutOccupancyLow || occupancy > cfgEvSel.cfgCutOccupancyHigh)) {
      return false;
    }
    if (fillCounter && eventCuts[UseOccupancy][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseOccupancy);
    }

    if (eventCuts[UseMultCorrCut][EvCut1]) {
      if (cfgFuncParas.cfgMultPVT0CCutEnabled) {
        if (multNTracksPV < cfgFuncParas.fMultPVT0CCutLow->Eval(cent)) {
          return false;
        }
        if (multNTracksPV > cfgFuncParas.fMultPVT0CCutHigh->Eval(cent)) {
          return false;
        }
      }
      if (cfgFuncParas.cfgMultT0CCutEnabled) {
        if (mult < cfgFuncParas.fMultT0CCutLow->Eval(cent)) {
          return false;
        }
        if (mult > cfgFuncParas.fMultT0CCutHigh->Eval(cent)) {
          return false;
        }
      }
      if (cfgFuncParas.cfgMultGlobalPVCutEnabled) {
        if (mult < cfgFuncParas.fMultGlobalPVCutLow->Eval(multNTracksPV)) {
          return false;
        }
        if (mult > cfgFuncParas.fMultGlobalPVCutHigh->Eval(multNTracksPV)) {
          return false;
        }
      }
      if (cfgFuncParas.cfgMultMultV0ACutEnabled) {
        if (collision.multFV0A() < cfgFuncParas.fMultMultV0ACutLow->Eval(mult)) {
          return false;
        }
        if (collision.multFV0A() > cfgFuncParas.fMultMultV0ACutHigh->Eval(mult)) {
          return false;
        }
      }
    }

    if (fillCounter && eventCuts[UseMultCorrCut][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseMultCorrCut);
    }

    // V0A T0A 5 sigma cut
    if (eventCuts[UseT0AV0ACut][EvCut1] && (std::fabs(collision.multFV0A() - cfgFuncParas.fT0AV0AMean->Eval(collision.multFT0A())) > cfgFuncParas.cfgV0AT0Acut * cfgFuncParas.fT0AV0ASigma->Eval(collision.multFT0A()))) {
      return false;
    }
    if (fillCounter && eventCuts[UseT0AV0ACut][EvCut1]) {
      histos.fill(HIST("hEventCount"), UseT0AV0ACut);
    }

    return true;
  }

  template <typename TCollision>
  void eventSelectedIndividually(TCollision const& collision)
  {
    histos.fill(HIST("hPassedEventSelection"), FilteredEvents);

    if (collision.sel8()) {
      histos.fill(HIST("hPassedEventSelection"), AfterSel8);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      histos.fill(HIST("hPassedEventSelection"), UseNoTimeFrameBorder);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      histos.fill(HIST("hPassedEventSelection"), UseNoITSROFrameBorder);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      histos.fill(HIST("hPassedEventSelection"), UseNoSameBunchPileup);
    }

    if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      histos.fill(HIST("hPassedEventSelection"), UseGoodZvtxFT0vsPV);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      histos.fill(HIST("hPassedEventSelection"), UseNoCollInTimeRangeStandard);
    }

    if (collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      histos.fill(HIST("hPassedEventSelection"), UseGoodITSLayersAll);
    }

    if (collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
      histos.fill(HIST("hPassedEventSelection"), UseGoodITSLayer0123);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      histos.fill(HIST("hPassedEventSelection"), UseNoCollInRofStandard);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      histos.fill(HIST("hPassedEventSelection"), UseNoHighMultCollInPrevRof);
    }

    auto occupancy = collision.trackOccupancyInTimeRange();
    if (eventCuts[UseOccupancy][EvCut1] && (occupancy > cfgEvSel.cfgCutOccupancyLow || occupancy < cfgEvSel.cfgCutOccupancyHigh)) {
      histos.fill(HIST("hPassedEventSelection"), UseOccupancy);
    }
  }

  double getPhiFT0(uint64_t chno, int i)
  {
    // offsetFT0[0]: FT0A, offsetFT0[1]: FT0C
    if (i > 1 || i < 0) {
      LOGF(fatal, "kFIT Index %d out of range", i);
    }

    ft0Det.calculateChannelCenter();
    auto chPos = ft0Det.getChannelCenter(chno);
    return RecoDecay::phi(chPos.X() + (*offsetFT0)[i].getX(), chPos.Y() + (*offsetFT0)[i].getY());
  }

  double getEtaFT0(uint64_t chno, int i)
  {
    // offsetFT0[0]: FT0A, offsetFT0[1]: FT0C
    if (i > 1 || i < 0) {
      LOGF(fatal, "kFIT Index %d out of range", i);
    }
    ft0Det.calculateChannelCenter();
    auto chPos = ft0Det.getChannelCenter(chno);
    auto x = chPos.X() + (*offsetFT0)[i].getX();
    auto y = chPos.Y() + (*offsetFT0)[i].getY();
    auto z = chPos.Z() + (*offsetFT0)[i].getZ();
    if (chno >= Ft0IndexA) {
      z = -z;
    }
    auto r = std::sqrt(x * x + y * y);
    auto theta = std::atan2(r, z);
    return -std::log(std::tan(0.5 * theta));
  }

  void loadAlignParam(uint64_t timestamp)
  {
    offsetFT0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FT0/Calib/Align", timestamp);
    if (offsetFT0 == nullptr) {
      LOGF(fatal, "Could not load FT0/Calib/Align for timestamp %d", timestamp);
    }
  }

  template <typename TTrack>
  bool trackSelected(TTrack const& track)
  {
    if (!cfgTrkSel.cfgTrackCutsDCAxy.value[UseGenTrkCuts].empty() && std::fabs(track.dcaXY()) > cfgFuncParas.fPtDepDCAxy->Eval(track.pt())) {
      return false;
    }

    if (!cfgTrkSel.cfgTrackCutsDCAz.value[UseGenTrkCuts].empty()) {
      if (std::fabs(track.dcaZ()) > cfgFuncParas.fPtDepDCAz->Eval(track.pt())) {
        return false;
      }
    } else {
      if (std::fabs(track.dcaZ()) > trackCuts[TrkCutDCAz][UseGenTrkCuts]) {
        return false;
      }
    }

    return ((track.pt() > trackCuts[TrkCutPtMin][UseGenTrkCuts]) && (track.pt() < trackCuts[TrkCutPtMax][UseGenTrkCuts]) && (std::abs(track.eta()) < trackCuts[TrkCutEtaCut][UseGenTrkCuts]) && (track.tpcChi2NCl() < trackCuts[TrkCutChi2PrTpcCls][UseGenTrkCuts]) && (track.tpcNClsFound() >= trackCuts[TrkCutTpcCluster][UseGenTrkCuts]) && (track.tpcNClsCrossedRows() >= trackCuts[TrkCutTpcCrossedRows][UseGenTrkCuts]) && (track.itsNCls() >= trackCuts[TrkCutItsCluster][UseGenTrkCuts]));
  }

  template <typename TTrack>
  bool trackSelectedForNch(TTrack const& track)
  {
    if (!cfgTrkSel.cfgTrackCutsDCAxy.value[UseNchSelCuts].empty() && (std::fabs(track.dcaXY()) > cfgFuncParas.fPtDepDCAxyForNch->Eval(track.pt()))) {
      return false;
    }
    if (!cfgTrkSel.cfgTrackCutsDCAz.value[UseNchSelCuts].empty()) {
      if (std::fabs(track.dcaZ()) > cfgFuncParas.fPtDepDCAzForNch->Eval(track.pt())) {
        return false;
      }
    } else {
      if (std::fabs(track.dcaZ()) > trackCuts[TrkCutDCAz][UseNchSelCuts]) {
        return false;
      }
    }
    return ((track.pt() > trackCuts[TrkCutPtMin][UseNchSelCuts]) && (track.pt() < trackCuts[TrkCutPtMax][UseNchSelCuts]) && (std::abs(track.eta()) < trackCuts[TrkCutEtaCut][UseNchSelCuts]) && (track.tpcChi2NCl() < trackCuts[TrkCutChi2PrTpcCls][UseNchSelCuts]) && (track.tpcNClsFound() >= trackCuts[TrkCutTpcCluster][UseNchSelCuts]) && (track.tpcNClsCrossedRows() >= trackCuts[TrkCutTpcCrossedRows][UseNchSelCuts]) && (track.itsNCls() >= trackCuts[TrkCutItsCluster][UseNchSelCuts]));
  }

  void loadGain(aod::BCsWithTimestamps::iterator const& bc)
  {
    cstFT0RelGain.clear();
    std::string fullPath;

    auto timestamp = bc.timestamp();
    constexpr int ChannelsFT0 = 208;
    if (cfgGeneral.cfgCorrLevel == 0) {
      for (auto i{0u}; i < ChannelsFT0; i++) {
        cstFT0RelGain.push_back(1.);
      }
    } else {
      fullPath = cfgGeneral.cfgGainEqPath;
      fullPath += "/FT0";
      const auto objft0Gain = ccdb->getForTimeStamp<std::vector<float>>(fullPath, timestamp);
      if (!objft0Gain) {
        for (auto i{0u}; i < ChannelsFT0; i++) {
          cstFT0RelGain.push_back(1.);
        }
      } else {
        cstFT0RelGain = *(objft0Gain);
      }
    }
  }

  template <typename TFT0s>
  void getChannel(TFT0s const& ft0, std::size_t const& iCh, int& id, float& ampl, int fitType, int system)
  {
    if (fitType == IndexFT0C) {
      id = ft0.channelC()[iCh];
      id = id + Ft0IndexA;
      ampl = ft0.amplitudeC()[iCh];
      if (system == SameEvent) {
        histos.fill(HIST("FT0Amp"), id, ampl);
      }
      ampl = ampl / cstFT0RelGain[id];
      if (system == SameEvent) {
        histos.fill(HIST("FT0AmpCorrect"), id, ampl);
      }
    } else if (fitType == IndexFT0A) {
      id = ft0.channelA()[iCh];
      ampl = ft0.amplitudeA()[iCh];
      if (system == SameEvent) {
        histos.fill(HIST("FT0Amp"), id, ampl);
      }
      ampl = ampl / cstFT0RelGain[id];
      if (system == SameEvent) {
        histos.fill(HIST("FT0AmpCorrect"), id, ampl);
      }
    } else {
      LOGF(fatal, "Cor Index %d out of range", fitType);
    }
  }

  template <typename TTrack>
  int getNsigmaPID(const TTrack& track)
  {
    // Computing Nsigma arrays for pion, kaon, and protons
    std::array<float, 3> nSigmaTPC = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::array<float, 3> nSigmaTOF = {track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};
    std::array<float, 3> nSigmaITS = {itsResponse.nSigmaITS<o2::track::PID::Pion>(track), itsResponse.nSigmaITS<o2::track::PID::Kaon>(track), itsResponse.nSigmaITS<o2::track::PID::Proton>(track)};
    int pid = -1; // -1 = not identified, 1 = pion, 2 = kaon, 3 = proton

    std::array<float, 3> nSigmaToUse = cfgPIDConfigs.cfgUseItsPID ? nSigmaITS : nSigmaTPC; // Choose which nSigma to use: TPC or ITS
    int useIndexDetector = UseTPC;
    useIndexDetector = cfgPIDConfigs.cfgUseItsPID ? UseITS : UseTPC; // Choose which nSigma to use: TPC or ITS

    bool isPion = false;
    bool isKaon = false;
    bool isProton = false;
    bool isDetectedPion = false;
    bool isDetectedKaon = false;
    bool isDetectedProton = false;
    bool isTofPion = false;
    bool isTofKaon = false;
    bool isTofProton = false;

    isDetectedPion = nSigmaToUse[IndexPionUp] < nSigmaVals[IndexPionUp][useIndexDetector] && nSigmaToUse[IndexPionUp] > nSigmaVals[IndexPionLow][useIndexDetector];
    isDetectedKaon = nSigmaToUse[IndexKaonUp] < nSigmaVals[IndexKaonUp][useIndexDetector] && nSigmaToUse[IndexKaonUp] > nSigmaVals[IndexKaonLow][useIndexDetector];
    isDetectedProton = nSigmaToUse[IndexProtonUp] < nSigmaVals[IndexProtonUp][useIndexDetector] && nSigmaToUse[IndexProtonUp] > nSigmaVals[IndexProtonLow][useIndexDetector];

    isTofPion = nSigmaTOF[IndexPionUp] < nSigmaVals[IndexPionUp][UseTof] && nSigmaTOF[IndexPionUp] > nSigmaVals[IndexPionLow][UseTof];
    isTofKaon = nSigmaTOF[IndexKaonUp] < nSigmaVals[IndexKaonUp][UseTof] && nSigmaTOF[IndexKaonUp] > nSigmaVals[IndexKaonLow][UseTof];
    isTofProton = nSigmaTOF[IndexProtonUp] < nSigmaVals[IndexProtonUp][UseTof] && nSigmaTOF[IndexProtonUp] > nSigmaVals[IndexProtonLow][UseTof];

    if (track.pt() > cfgPIDConfigs.cfgTofPtCut && !track.hasTOF()) {
      return -1;
    }
    if (track.pt() > cfgPIDConfigs.cfgTofPtCut && track.hasTOF()) {
      isPion = isTofPion && isDetectedPion;
      isKaon = isTofKaon && isDetectedKaon;
      isProton = isTofProton && isDetectedProton;
    } else {
      isPion = isDetectedPion;
      isKaon = isDetectedKaon;
      isProton = isDetectedProton;
    }

    if (cfgPIDConfigs.cfgPIDUseRejection && ((isPion && isKaon) || (isPion && isProton) || (isKaon && isProton))) {
      return -1; // more than one particle satisfy the criteria
    }

    if (isPion) {
      pid = UsePions;
    } else if (isKaon) {
      pid = UseKaons;
    } else if (isProton) {
      pid = UseProtons;
    } else {
      return -1; // no particle satisfies the criteria
    }

    return pid; // -1 = not identified, 1 = pion, 2 = kaon, 3 = proton
  }

  void loadCorrection(uint64_t timestamp)
  {
    if (correctionsLoaded) {
      return;
    }
    if (!cfgEffWeights.cfgEfficiency.value.empty()) {
      if (cfgEffWeights.cfgLocalEfficiency) {
        TFile* fEfficiencyTrigger = nullptr;
        fEfficiencyTrigger = TFile::Open(cfgEffWeights.cfgEfficiency.value.c_str(), "READ");
        mEfficiency = dynamic_cast<TH3D*>(fEfficiencyTrigger->Get("ccdb_object"));

      } else {
        mEfficiency = ccdb->getForTimeStamp<TH3D>(cfgEffWeights.cfgEfficiency, timestamp);
      }
      if (mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEffWeights.cfgEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEffWeights.cfgEfficiency.value.c_str(), static_cast<void*>(mEfficiency));
    }
    if (!cfgEffWeights.cfgEfficiencyNch.value.empty()) {
      if (cfgEffWeights.cfgLocalEfficiencyNch) {
        TFile* fEfficiencyTrigger = nullptr;
        fEfficiencyTrigger = TFile::Open(cfgEffWeights.cfgEfficiencyNch.value.c_str(), "READ");
        mEfficiencyNch = dynamic_cast<TH1D*>(fEfficiencyTrigger->Get("ccdb_object"));

      } else {
        mEfficiencyNch = ccdb->getForTimeStamp<TH1D>(cfgEffWeights.cfgEfficiencyNch, timestamp);
      }
      if (!mEfficiencyNch) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEffWeights.cfgEfficiencyNch.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEffWeights.cfgEfficiencyNch.value.c_str(), static_cast<void*>(mEfficiencyNch));
    }
    if (!cfgEffWeights.cfgCentralityWeight.value.empty()) {
      mCentralityWeight = ccdb->getForTimeStamp<TH1D>(cfgEffWeights.cfgCentralityWeight, timestamp);
      if (mCentralityWeight == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEffWeights.cfgCentralityWeight.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEffWeights.cfgCentralityWeight.value.c_str(), static_cast<void*>(mCentralityWeight));
    }
    correctionsLoaded = true;
  }

  bool getEfficiencyCorrectionNch(float& weightNch, float pt)
  {
    float effNch = 1.;
    if (mEfficiencyNch) {
      int ptBin = mEfficiencyNch->FindBin(pt);
      effNch = mEfficiencyNch->GetBinContent(ptBin);
    } else {
      effNch = 1.0;
    }
    if (effNch == 0) {
      return false;
    }
    weightNch = 1. / effNch;
    return true;
  }

  bool getEfficiencyCorrection(float& weight, float pt, float eta, float vertex)
  {
    float eff = 1.;
    if (mEfficiency) {
      int etaBin = mEfficiency->GetXaxis()->FindBin(eta); // use the eta bin corresponding to eta=0 for the trigger particle efficiency
      int ptBin = mEfficiency->GetYaxis()->FindBin(pt);
      int vertexBin = mEfficiency->GetZaxis()->FindBin(vertex); // use the vertex bin corresponding to z=0 for the trigger particle efficiency
      eff = mEfficiency->GetBinContent(etaBin, ptBin, vertexBin);
    } else {
      eff = 1.0;
    }
    if (eff == 0) {
      return false;
    }
    weight = 1. / eff;
    return true;
  }

  bool getCentralityWeight(float& weightCent, const double centrality)
  {
    float weight = 1.;
    if (mCentralityWeight) {
      weight = mCentralityWeight->GetBinContent(mCentralityWeight->FindBin(centrality));
    } else {
      weight = 1.0;
    }
    if (weight == 0) {
      return false;
    }
    weightCent = weight;
    return true;
  }

  template <typename TTracks>
  void nchCounter(const TTracks& tracks, double& multiplicity, bool fillHistogram) // function to count the number of tracks in the event and fill the histogram
  {
    double nTracksCorrected = 0;
    double nTracksUncorrected = 0;
    float weightNch = 1.0f;
    for (auto const& track : tracks) {
      if (!trackSelectedForNch(track)) {
        continue;
      }

      if (!getEfficiencyCorrectionNch(weightNch, track.pt())) {
        continue;
      }

      nTracksUncorrected += 1.0;
      nTracksCorrected += weightNch;
    }
    if (fillHistogram) {
      histos.fill(HIST("hTrackCorrection2d"), nTracksUncorrected, nTracksCorrected);
    }
    multiplicity = nTracksCorrected;
  }

  template <typename TCollision, typename TTracks>
  void fillYieldFt0aFt0c(const TCollision& collision, const TTracks& tracks) // function to fill the yield and etaphi histograms.
  {

    float weff1 = 1.0;
    float zvtx = collision.posZ();

    for (auto const& track1 : tracks) {

      if (!trackSelected(track1)) {
        continue;
      }

      if (!getEfficiencyCorrection(weff1, track1.pt(), track1.eta(), zvtx)) {
        continue;
      }

      histos.fill(HIST("Phi_FT0A_FT0C"), RecoDecay::constrainAngle(track1.phi(), 0.0));
      histos.fill(HIST("Eta_FT0A_FT0C"), track1.eta());
      histos.fill(HIST("EtaCorrected_FT0A_FT0C"), track1.eta(), weff1);
      histos.fill(HIST("pTFiner_FT0A_FT0C"), track1.pt());
      histos.fill(HIST("pTFinerCorrected_FT0A_FT0C"), track1.pt(), weff1);
    }
  }

  template <typename TTrack>
  void fillTrackQA(TTrack const& track, int pid, double eventClass, Stage stage, float weff) // function to fill the QA after Nsigma selection
  {
    if (cfgGeneral.cfgQABasic) {
      if (stage == Before && pid == UseCharged) {
        histos.fill(HIST("Phi"), RecoDecay::constrainAngle(track.phi(), 0.0), UseCharged);
        histos.fill(HIST("Eta"), track.eta(), UseCharged);
        histos.fill(HIST("EtaCorrected"), track.eta(), UseCharged, weff);
        histos.fill(HIST("pTFiner"), track.pt(), UseCharged, eventClass);
        histos.fill(HIST("pTFinerCorrected"), track.pt(), UseCharged, eventClass, weff);

      } else if (stage == After && pid > UseCharged) {
        histos.fill(HIST("Phi"), RecoDecay::constrainAngle(track.phi(), 0.0), pid);
        histos.fill(HIST("Eta"), track.eta(), pid);
        histos.fill(HIST("EtaCorrected"), track.eta(), pid, weff);
        histos.fill(HIST("pTFiner"), track.pt(), pid, eventClass);
        histos.fill(HIST("pTFinerCorrected"), track.pt(), pid, eventClass, weff);
      }
    }

    if (pid > UseCharged) {
      const bool useITS = cfgPIDConfigs.cfgUseItsPID;
      double tpcNSigma = 0.0;
      double tpcExpSigma = 0.0;
      double tofNSigma = 0.0;
      double itsNSigma = 0.0;
      switch (pid) {
        case UsePions: // For Pions
          tofNSigma = track.tofNSigmaPi();
          if (!useITS) {
            tpcNSigma = track.tpcNSigmaPi();
            tpcExpSigma = track.tpcExpSigmaPi();
          } else {
            itsNSigma = itsResponse.nSigmaITS<o2::track::PID::Pion>(track);
          }
          break;
        case UseKaons: // For Kaons
          tofNSigma = track.tofNSigmaKa();
          if (!useITS) {
            tpcNSigma = track.tpcNSigmaKa();
            tpcExpSigma = track.tpcExpSigmaKa();
          } else {
            itsNSigma = itsResponse.nSigmaITS<o2::track::PID::Kaon>(track);
          }
          break;

        case UseProtons: // For Protons
          tofNSigma = track.tofNSigmaPr();
          if (!useITS) {
            tpcNSigma = track.tpcNSigmaPr();
            tpcExpSigma = track.tpcExpSigmaPr();
          } else {
            itsNSigma = itsResponse.nSigmaITS<o2::track::PID::Proton>(track);
          }
          break;
        default:
          return;
      }

      if (stage == Before) {
        if (!useITS) {
          histos.fill(HIST("TofTpcNsigma_beforeCut"), tpcNSigma, tofNSigma, track.pt(), pid);
        } else {
          histos.fill(HIST("TofItsNsigma_beforeCut"), itsNSigma, tofNSigma, track.pt(), pid);
        }
      } else if (stage == After) {
        if (!useITS) {
          histos.fill(HIST("TofTpcNsigma_afterCut"), tpcNSigma, tofNSigma, track.pt(), pid);
        } else {
          histos.fill(HIST("TofItsNsigma_afterCut"), itsNSigma, tofNSigma, track.pt(), pid);
        }
      }

      if (cfgPIDConfigs.cfgQAdEdx && !useITS) {
        double tpcExpSignal = track.tpcSignal() - (tpcNSigma * tpcExpSigma);

        if (stage == Before) {
          histos.fill(HIST("TpcdEdx_ptwise_beforeCut"), track.pt(), track.tpcSignal(), tofNSigma, pid);
          histos.fill(HIST("ExpTpcdEdx_ptwise_beforeCut"), track.pt(), tpcExpSignal, tofNSigma, pid);
          histos.fill(HIST("ExpSigma_ptwise_beforeCut"), track.pt(), tpcExpSigma, tofNSigma, pid);
        } else if (stage == After) {
          histos.fill(HIST("TpcdEdx_ptwise_afterCut"), track.pt(), track.tpcSignal(), tofNSigma, pid);
          histos.fill(HIST("ExpTpcdEdx_ptwise_afterCut"), track.pt(), tpcExpSignal, tofNSigma, pid);
          histos.fill(HIST("ExpSigma_ptwise_afterCut"), track.pt(), tpcExpSigma, tofNSigma, pid);
        }
      }
    } // if not a charged particle
  }

  template <typename TTracks, typename TFT0s>
  void fillCorrelationsTPCFT0(const TTracks& tracks1, TFT0s const& ft0, float posZ, int system, double eventClass, int corType, float eventWeight) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {
    int fSampleIndex = static_cast<int>(fRandom.Uniform(0.0, cfgGeneral.cfgSampleSize));

    float triggerWeight = 1.0f;
    // loop over all tracks
    for (auto const& track1 : tracks1) {

      if (!trackSelected(track1)) {
        continue;
      }
      if (!getEfficiencyCorrection(triggerWeight, track1.pt(), track1.eta(), posZ)) {
        continue;
      }

      fillTrackQA(track1, UseCharged, eventClass, Before, triggerWeight);
      fillTrackQA(track1, UsePions, eventClass, Before, triggerWeight);
      fillTrackQA(track1, UseKaons, eventClass, Before, triggerWeight);
      fillTrackQA(track1, UseProtons, eventClass, Before, triggerWeight);

      int pidIndex = getNsigmaPID(track1);

      if (pidIndex > UseCharged) {
        fillTrackQA(track1, pidIndex, eventClass, After, triggerWeight);
      }

      if (system == SameEvent) {
        if (corType == IndexFT0C) {
          histos.fill(HIST("Trig_hist_TPC_FT0C"), fSampleIndex, posZ, eventClass, track1.pt(), UseCharged, eventWeight * triggerWeight);
          if (pidIndex > UseCharged) {
            histos.fill(HIST("Trig_hist_TPC_FT0C"), fSampleIndex, posZ, eventClass, track1.pt(), pidIndex, eventWeight * triggerWeight);
          }
        } else if (corType == IndexFT0A) {
          histos.fill(HIST("Trig_hist_TPC_FT0A"), fSampleIndex, posZ, eventClass, track1.pt(), UseCharged, eventWeight * triggerWeight);
          if (pidIndex > UseCharged) {
            histos.fill(HIST("Trig_hist_TPC_FT0A"), fSampleIndex, posZ, eventClass, track1.pt(), pidIndex, eventWeight * triggerWeight);
          }
        }
      }

      std::size_t channelSize = 0;
      if (corType == IndexFT0C) {
        channelSize = ft0.channelC().size();
      } else if (corType == IndexFT0A) {
        channelSize = ft0.channelA().size();
      } else {
        LOGF(fatal, "Cor Index %d out of range", corType);
      }
      for (std::size_t iCh = 0; iCh < channelSize; iCh++) {
        int chanelid = 0;
        float ampl = 0.;
        getChannel(ft0, iCh, chanelid, ampl, corType, system);

        auto phi = getPhiFT0(chanelid, corType);
        auto eta = getEtaFT0(chanelid, corType);

        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - phi, -PIHalf);
        float deltaEta = track1.eta() - eta;
        // fill the right sparse and histograms
        if (system == SameEvent) {
          if (corType == IndexFT0A) {
            if (cfgGeneral.cfgQABasic) {
              histos.fill(HIST("Assoc_amp_same_TPC_FT0A"), chanelid, ampl);
              histos.fill(HIST("deltaEta_deltaPhi_same_TPC_FT0A"), deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
            }
            sameTpcFt0a->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, deltaPhi, deltaEta, track1.pt(), UseCharged, ampl * eventWeight * triggerWeight);
            if (pidIndex > UseCharged) {
              sameTpcFt0a->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, deltaPhi, deltaEta, track1.pt(), pidIndex, ampl * eventWeight * triggerWeight);
            }
          } else if (corType == IndexFT0C) {
            if (cfgGeneral.cfgQABasic) {
              histos.fill(HIST("Assoc_amp_same_TPC_FT0C"), chanelid, ampl);
              histos.fill(HIST("deltaEta_deltaPhi_same_TPC_FT0C"), deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
            }
            sameTpcFt0c->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, deltaPhi, deltaEta, track1.pt(), UseCharged, ampl * eventWeight * triggerWeight);
            if (pidIndex > UseCharged) {
              sameTpcFt0c->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, deltaPhi, deltaEta, track1.pt(), pidIndex, ampl * eventWeight * triggerWeight);
            }
          }
        } else if (system == MixedEvent) {
          if (corType == IndexFT0A) {
            if (cfgGeneral.cfgQABasic) {
              histos.fill(HIST("Assoc_amp_mixed_TPC_FT0A"), chanelid, ampl);
              histos.fill(HIST("deltaEta_deltaPhi_mixed_TPC_FT0A"), deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
            }
            mixedTpcFt0a->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, deltaPhi, deltaEta, track1.pt(), UseCharged, ampl * eventWeight * triggerWeight);
            if (pidIndex > UseCharged) {
              mixedTpcFt0a->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, deltaPhi, deltaEta, track1.pt(), pidIndex, ampl * eventWeight * triggerWeight);
            }
          } else if (corType == IndexFT0C) {
            if (cfgGeneral.cfgQABasic) {
              histos.fill(HIST("Assoc_amp_mixed_TPC_FT0C"), chanelid, ampl);
              histos.fill(HIST("deltaEta_deltaPhi_mixed_TPC_FT0C"), deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
            }
            mixedTpcFt0c->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, deltaPhi, deltaEta, track1.pt(), UseCharged, ampl * eventWeight * triggerWeight);
            if (pidIndex > UseCharged) {
              mixedTpcFt0c->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, deltaPhi, deltaEta, track1.pt(), pidIndex, ampl * eventWeight * triggerWeight);
            }
          }
        }
      }
    }
  }

  template <typename TFT0s>
  void fillCorrelationsFT0AFT0C(TFT0s const& ft0Col1, TFT0s const& ft0Col2, float posZ, int system, double eventClass, float eventWeight) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {
    int fSampleIndex = static_cast<int>(fRandom.Uniform(0.0, cfgGeneral.cfgSampleSize));

    float triggerWeight = 1.0f;
    std::size_t channelASize = ft0Col1.channelA().size();
    std::size_t channelCSize = ft0Col2.channelC().size();
    // loop over all tracks
    for (std::size_t iChA = 0; iChA < channelASize; iChA++) {

      int chanelAid = 0;
      float amplA = 0.;
      getChannel(ft0Col1, iChA, chanelAid, amplA, IndexFT0A, system);
      auto phiA = getPhiFT0(chanelAid, IndexFT0A);
      auto etaA = getEtaFT0(chanelAid, IndexFT0A);

      if (system == SameEvent) {
        histos.fill(HIST("Trig_hist_FT0A_FT0C"), fSampleIndex, posZ, eventClass, eventWeight * amplA);
      }

      for (std::size_t iChC = 0; iChC < channelCSize; iChC++) {
        int chanelCid = 0;
        float amplC = 0.;
        getChannel(ft0Col2, iChC, chanelCid, amplC, IndexFT0C, system);
        auto phiC = getPhiFT0(chanelCid, IndexFT0C);
        auto etaC = getEtaFT0(chanelCid, IndexFT0C);
        float deltaPhi = RecoDecay::constrainAngle(phiA - phiC, -PIHalf);
        float deltaEta = etaA - etaC;

        // fill the right sparse and histograms
        if (system == SameEvent) {
          if (cfgGeneral.cfgQABasic) {
            histos.fill(HIST("deltaEta_deltaPhi_same_FT0A_FT0C"), deltaPhi, deltaEta, amplA * amplC * eventWeight * triggerWeight);
          }
          sameFt0aFt0c->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, deltaPhi, deltaEta, amplA * amplC * eventWeight * triggerWeight);
        } else if (system == MixedEvent) {
          if (cfgGeneral.cfgQABasic) {
            histos.fill(HIST("deltaEta_deltaPhi_mixed_FT0A_FT0C"), deltaPhi, deltaEta, amplA * amplC * eventWeight * triggerWeight);
          }
          mixedFt0aFt0c->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, deltaPhi, deltaEta, amplA * amplC * eventWeight * triggerWeight);
        }
      }
    }
  }

  bool isGoodRun(int runNumber)
  {
    for (const auto& excludedRun : cfgGeneral.cfgRunRemoveList.value) {
      if (runNumber == excludedRun) {
        return false;
      }
    }

    return true;
  }

  void processSameTpcFt0(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
  {
    if (cfgGeneral.cfgQABasic) {
      eventSelectedIndividually(collision);
    }

    if (!eventRct(collision, true)) {
      return;
    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    int currentRunNumber = bc.runNumber();
    if (!cfgGeneral.cfgRunRemoveList.value.empty()) {
      if (!isGoodRun(currentRunNumber)) { // Rejects runs if bad run number
        return;
      }
    }

    auto cent = getCentrality(collision);
    if (cfgEvSel.cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), cent, true)) {
      return;
    }

    if (!collision.has_foundFT0()) {
      return;
    }

    histos.fill(HIST("hEventCount"), HaveFT0Cut);

    loadAlignParam(bc.timestamp());
    loadGain(bc);
    loadCorrection(bc.timestamp());
    float weightCent = 1.0f;

    if (!getCentralityWeight(weightCent, cent)) {
      return; // or continue in a loop
    }

    histos.fill(HIST("Centrality"), cent);
    histos.fill(HIST("CentralityWeighted"), cent, weightCent);

    histos.fill(HIST("zVtx"), collision.posZ());

    histos.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin

    auto multiplicity = static_cast<double>(tracks.size());

    if (cfgGeneral.cfgQABasic) {
      histos.fill(HIST("Nch"), multiplicity);
    }

    if (cfgGeneral.cfgStrictTrackCounter) {
      nchCounter(tracks, multiplicity, false);
    }

    if (cfgGeneral.cfgQABasic) {
      histos.fill(HIST("Nch_corrected"), multiplicity);
    }

    if (cfgEvSel.cfgSelCollByNch && (multiplicity > cfgEvSel.cfgMaxMultForCorrelations || multiplicity < cfgEvSel.cfgMinMultForCorrelations)) {
      return;
    }
    if (!cfgEvSel.cfgSelCollByNch && (cent > cfgEvSel.cfgMaxCentForCorrelations || cent < cfgEvSel.cfgMinCentForCorrelations)) {
      return;
    }

    auto eventClass = (cfgEvSel.cfgSelCollByNch) ? multiplicity : cent;

    const auto& ft0 = collision.foundFT0();
    if (cfgGeneral.cfgAnalyzeTPCFT0A) {
      fillCorrelationsTPCFT0(tracks, ft0, collision.posZ(), SameEvent, eventClass, IndexFT0A, weightCent);
    }
    if (cfgGeneral.cfgAnalyzeTPCFT0C) {
      fillCorrelationsTPCFT0(tracks, ft0, collision.posZ(), SameEvent, eventClass, IndexFT0C, weightCent);
    }
  }
  PROCESS_SWITCH(PidLongRange, processSameTpcFt0, "Process same event for TPC-FT0 correlation", false);

  void processMixedTpcFt0(FilteredCollisions const& collisions, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
  {

    auto getTracksSize = [&tracks, this](FilteredCollisions::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}, true};

    auto tracksTuple = std::make_tuple(tracks, tracks);
    Pair<FilteredCollisions, FilteredTracks, FilteredTracks, MixedBinning> pairs{binningOnVtxAndMult, cfgGeneral.cfgMinMixEventNum, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;

      if (!collision1.sel8() || !collision2.sel8()) { // Doing this at the top of mixed events loop to reduce computation (This cut is also there inside eventSelected)
        continue;
      }

      if (!eventRct(collision1, false) || !eventRct(collision2, false)) {
        continue;
      }

      auto cent1 = getCentrality(collision1);
      auto cent2 = getCentrality(collision2);

      if (cfgEvSel.cfgUseAdditionalEventCut && !eventSelected(collision1, tracks1.size(), cent1, false)) {
        continue;
      }
      if (cfgEvSel.cfgUseAdditionalEventCut && !eventSelected(collision2, tracks2.size(), cent2, false)) {
        continue;
      }

      if (!(collision1.has_foundFT0() && collision2.has_foundFT0())) {
        continue;
      }

      histos.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      int currentRunNumber = bc.runNumber();
      if (!cfgGeneral.cfgRunRemoveList.value.empty()) {
        if (!isGoodRun(currentRunNumber)) { // Rejects runs if bad run number
          continue;
        }
      }

      loadAlignParam(bc.timestamp());
      loadCorrection(bc.timestamp());

      auto multiplicity = static_cast<double>(tracks1.size());

      if (cfgGeneral.cfgStrictTrackCounter) {
        nchCounter(tracks1, multiplicity, false);
      }

      if (cfgEvSel.cfgSelCollByNch && (multiplicity > cfgEvSel.cfgMaxMultForCorrelations || multiplicity < cfgEvSel.cfgMinMultForCorrelations)) {
        continue;
      }
      if (!cfgEvSel.cfgSelCollByNch && (cent1 > cfgEvSel.cfgMaxCentForCorrelations || cent1 < cfgEvSel.cfgMinCentForCorrelations)) {
        continue;
      }
      if (!cfgEvSel.cfgSelCollByNch && (cent2 > cfgEvSel.cfgMaxCentForCorrelations || cent2 < cfgEvSel.cfgMinCentForCorrelations)) {
        continue;
      }

      auto eventClass = (cfgEvSel.cfgSelCollByNch) ? multiplicity : cent1;

      float eventWeight = 1.0f;

      if (cfgEffWeights.cfgUseEventWeights) {
        eventWeight = 1.0f / it.currentWindowNeighbours();
      }
      float weightCent = 1.0f;
      getCentralityWeight(weightCent, cent1);

      const auto& ft0 = collision2.foundFT0();
      if (cfgGeneral.cfgAnalyzeTPCFT0A) {
        fillCorrelationsTPCFT0(tracks1, ft0, collision1.posZ(), MixedEvent, eventClass, IndexFT0A, eventWeight * weightCent);
      }
      if (cfgGeneral.cfgAnalyzeTPCFT0C) {
        fillCorrelationsTPCFT0(tracks1, ft0, collision1.posZ(), MixedEvent, eventClass, IndexFT0C, eventWeight * weightCent);
      }
    }
  }
  PROCESS_SWITCH(PidLongRange, processMixedTpcFt0, "Process mixed events for TPC-FT0A correlation", false);

  void processSameFt0aFt0c(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
  {

    if (cfgGeneral.cfgQABasic) {
      eventSelectedIndividually(collision);
    }

    if (!eventRct(collision, true)) {
      return;
    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int currentRunNumber = bc.runNumber();
    if (!cfgGeneral.cfgRunRemoveList.value.empty()) {
      if (!isGoodRun(currentRunNumber)) { // Rejects runs if bad run number
        return;
      }
    }

    auto cent = getCentrality(collision);
    if (cfgEvSel.cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), cent, true)) {
      return;
    }

    if (!collision.has_foundFT0()) {
      return;
    }
    histos.fill(HIST("hEventCount"), HaveFT0Cut);

    loadAlignParam(bc.timestamp());
    loadGain(bc);
    loadCorrection(bc.timestamp());
    float weightCent = 1.0f;

    getCentralityWeight(weightCent, cent);
    histos.fill(HIST("Centrality"), cent);
    histos.fill(HIST("CentralityWeighted"), cent, weightCent);

    histos.fill(HIST("zVtx"), collision.posZ());
    histos.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin

    fillYieldFt0aFt0c(collision, tracks);

    auto multiplicity = static_cast<double>(tracks.size());

    if (cfgGeneral.cfgQABasic) {
      histos.fill(HIST("Nch"), multiplicity);
    }

    if (cfgGeneral.cfgStrictTrackCounter) {
      nchCounter(tracks, multiplicity, false);
    }

    if (cfgGeneral.cfgQABasic) {
      histos.fill(HIST("Nch_corrected"), multiplicity);
    }

    if (cfgEvSel.cfgSelCollByNch && (multiplicity > cfgEvSel.cfgMaxMultForCorrelations || multiplicity < cfgEvSel.cfgMinMultForCorrelations)) {
      return;
    }
    if (!cfgEvSel.cfgSelCollByNch && (cent > cfgEvSel.cfgMaxCentForCorrelations || cent < cfgEvSel.cfgMinCentForCorrelations)) {
      return;
    }

    auto eventClass = (cfgEvSel.cfgSelCollByNch) ? multiplicity : cent;

    const auto& ft0 = collision.foundFT0();

    fillCorrelationsFT0AFT0C(ft0, ft0, collision.posZ(), SameEvent, eventClass, weightCent);
  }
  PROCESS_SWITCH(PidLongRange, processSameFt0aFt0c, "Process same event for FT0A-FT0C correlation", true);

  void processMixedFt0aFt0c(FilteredCollisions const& collisions, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
  {

    auto getTracksSize = [&tracks, this](FilteredCollisions::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}, true};

    auto tracksTuple = std::make_tuple(tracks, tracks);
    Pair<FilteredCollisions, FilteredTracks, FilteredTracks, MixedBinning> pairs{binningOnVtxAndMult, cfgGeneral.cfgMinMixEventNum, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;

      // should have the same event to TPC-FT0A/C correlations
      if (!collision1.sel8() || !collision2.sel8()) {
        continue;
      }

      if (!eventRct(collision1, false) || !eventRct(collision2, false)) {
        continue;
      }

      auto cent1 = getCentrality(collision1);
      auto cent2 = getCentrality(collision2);

      if (cfgEvSel.cfgUseAdditionalEventCut && !eventSelected(collision1, tracks1.size(), cent1, false)) {
        continue;
      }
      if (cfgEvSel.cfgUseAdditionalEventCut && !eventSelected(collision2, tracks2.size(), cent2, false)) {
        continue;
      }

      if (!(collision1.has_foundFT0() && collision2.has_foundFT0())) {
        continue;
      }

      histos.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      int currentRunNumber = bc.runNumber();
      if (!cfgGeneral.cfgRunRemoveList.value.empty()) {
        if (!isGoodRun(currentRunNumber)) { // Rejects runs if bad run number
          continue;
        }
      }
      loadAlignParam(bc.timestamp());
      loadCorrection(bc.timestamp());

      auto multiplicity = static_cast<double>(tracks1.size());

      if (cfgGeneral.cfgStrictTrackCounter) {
        nchCounter(tracks1, multiplicity, false);
      }

      if (cfgEvSel.cfgSelCollByNch && (multiplicity > cfgEvSel.cfgMaxMultForCorrelations || multiplicity < cfgEvSel.cfgMinMultForCorrelations)) {
        continue;
      }
      if (!cfgEvSel.cfgSelCollByNch && (cent1 > cfgEvSel.cfgMaxCentForCorrelations || cent1 < cfgEvSel.cfgMinCentForCorrelations)) {
        continue;
      }
      if (!cfgEvSel.cfgSelCollByNch && (cent2 > cfgEvSel.cfgMaxCentForCorrelations || cent2 < cfgEvSel.cfgMinCentForCorrelations)) {
        continue;
      }

      auto eventClass = (cfgEvSel.cfgSelCollByNch) ? multiplicity : cent1;

      float eventWeight = 1.0f;

      if (cfgEffWeights.cfgUseEventWeights) {
        eventWeight = 1.0f / it.currentWindowNeighbours();
      }
      float weightCent = 1.0f;
      getCentralityWeight(weightCent, cent1);

      const auto& ft0Col1 = collision1.foundFT0();
      const auto& ft0Col2 = collision2.foundFT0();

      fillCorrelationsFT0AFT0C(ft0Col1, ft0Col2, collision1.posZ(), MixedEvent, eventClass, eventWeight * weightCent);
    }
  }
  PROCESS_SWITCH(PidLongRange, processMixedFt0aFt0c, "Process mixed events for FT0A-FT0C correlation", true);

  void processQA(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
  {
    eventSelectedIndividually(collision);

    if (!eventRct(collision, true)) {
      return;
    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int currentRunNumber = bc.runNumber();
    if (!cfgGeneral.cfgRunRemoveList.value.empty()) {
      if (!isGoodRun(currentRunNumber)) { // Rejects runs if bad run number
        return;
      }
    }

    histos.fill(HIST("h_globalTracks_centT0C_before"), collision.centFT0C(), tracks.size());
    histos.fill(HIST("h_PVTracks_centT0C_before"), collision.centFT0C(), collision.multNTracksPV());
    histos.fill(HIST("h_globalTracks_PVTracks_before"), collision.multNTracksPV(), tracks.size());
    histos.fill(HIST("h_globalTracks_multT0A_before"), collision.multFT0A(), tracks.size());
    histos.fill(HIST("h_globalTracks_multV0A_before"), collision.multFV0A(), tracks.size());
    histos.fill(HIST("h_multV0A_multT0A_before"), collision.multFT0A(), collision.multFV0A());
    histos.fill(HIST("h_multT0C_centT0C_before"), collision.centFT0C(), collision.multFT0C());

    auto cent = getCentrality(collision);
    if (cfgEvSel.cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), cent, true)) {
      return;
    }

    histos.fill(HIST("h_globalTracks_centT0C_after"), collision.centFT0C(), tracks.size());
    histos.fill(HIST("h_PVTracks_centT0C_after"), collision.centFT0C(), collision.multNTracksPV());
    histos.fill(HIST("h_globalTracks_PVTracks_after"), collision.multNTracksPV(), tracks.size());
    histos.fill(HIST("h_globalTracks_multT0A_after"), collision.multFT0A(), tracks.size());
    histos.fill(HIST("h_globalTracks_multV0A_after"), collision.multFV0A(), tracks.size());
    histos.fill(HIST("h_multV0A_multT0A_after"), collision.multFT0A(), collision.multFV0A());
    histos.fill(HIST("h_multT0C_centT0C_after"), collision.centFT0C(), collision.multFT0C());
    histos.fill(HIST("h_centFT0M_centFT0C"), collision.centFT0C(), collision.centFT0M());
    histos.fill(HIST("h_centFV0A_centFT0C"), collision.centFT0C(), collision.centFV0A());

    if (!collision.has_foundFT0()) {
      return;
    }
    histos.fill(HIST("hEventCount"), HaveFT0Cut);

    loadAlignParam(bc.timestamp());
    loadGain(bc);
    loadCorrection(bc.timestamp());
    float weightCent = 1.0f;

    getCentralityWeight(weightCent, cent);
    histos.fill(HIST("Centrality"), cent);
    histos.fill(HIST("CentralityWeighted"), cent, weightCent);

    histos.fill(HIST("zVtx"), collision.posZ());

    auto multiplicity = static_cast<double>(tracks.size());

    histos.fill(HIST("Nch"), multiplicity);

    if (cfgGeneral.cfgStrictTrackCounter) {
      nchCounter(tracks, multiplicity, true);
    }

    histos.fill(HIST("Nch_corrected"), multiplicity);

    auto eventClass = (cfgEvSel.cfgSelCollByNch) ? multiplicity : cent;

    float triggerWeight = 1.0f;
    // loop over all tracks
    for (auto const& track : tracks) {

      histos.fill(HIST("hDCAz_before"), track.dcaZ(), track.pt());
      histos.fill(HIST("hDCAxy_before"), track.dcaXY(), track.pt());
      histos.fill(HIST("hTPCNclsPID_before"), track.tpcNClsPID());
      histos.fill(HIST("hTPCNclsFound_before"), track.tpcNClsFound());
      histos.fill(HIST("hTPCCrossedRows_before"), track.tpcNClsCrossedRows());

      if (!trackSelected(track)) {
        continue;
      }

      histos.fill(HIST("hTPCNclsPID_after"), track.tpcNClsPID());
      histos.fill(HIST("hTPCNclsFound_after"), track.tpcNClsFound());
      histos.fill(HIST("hTPCCrossedRows_after"), track.tpcNClsCrossedRows());
      histos.fill(HIST("hChi2prTPCcls"), track.tpcChi2NCl());
      histos.fill(HIST("hChi2prITScls"), track.itsChi2NCl());
      histos.fill(HIST("hITSNclsFound"), track.itsNCls());

      histos.fill(HIST("hDCAz_after"), track.dcaZ(), track.pt(), UseCharged);
      histos.fill(HIST("hDCAxy_after"), track.dcaXY(), track.pt(), UseCharged);

      fillTrackQA(track, UseCharged, eventClass, Before, triggerWeight);
      fillTrackQA(track, UsePions, eventClass, Before, triggerWeight);
      fillTrackQA(track, UseKaons, eventClass, Before, triggerWeight);
      fillTrackQA(track, UseProtons, eventClass, Before, triggerWeight);

      int pidIndex = getNsigmaPID(track);

      if (pidIndex > UseCharged) {
        fillTrackQA(track, pidIndex, eventClass, After, triggerWeight);
        histos.fill(HIST("hDCAz_after"), track.dcaZ(), track.pt(), pidIndex);
        histos.fill(HIST("hDCAxy_after"), track.dcaXY(), track.pt(), pidIndex);
      }

    } // end of track loop
  }
  PROCESS_SWITCH(PidLongRange, processQA, "Process for storing the QA for events and tracks", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PidLongRange>(cfgc),
  };
}
