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

/// \file corrFit.cxx
/// \brief Ultra long range correlation using forward FIT detectors and TPC, with focus on multiplicity dependence
/// \author Thor Jensen (thor.kjaersgaard.jensen@cern.ch)

#include "PWGCF/Core/CorrelationContainer.h"

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
#include <DataFormatsParameters/GRPMagField.h>
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

struct CorrFit {
  o2::aod::ITSResponse itsResponse;
  Service<ccdb::BasicCCDBManager> ccdb{};

  O2_DEFINE_CONFIGURABLE(cfgUseTransverseMomentum, bool, false, "Use transverse momentum for correlation container")
  O2_DEFINE_CONFIGURABLE(cfgQaCheck, bool, true, "Enable QA histograms for event selection")
  O2_DEFINE_CONFIGURABLE(cfgStrictTrackCounter, bool, false, "Strict track counter for multiplicity correlation cut, counts only tracks that pass all cuts and are used in the correlation")
  O2_DEFINE_CONFIGURABLE(cfgUseRefpTCut, bool, false, "Apply upper pT cut on reference tracks")
  O2_DEFINE_CONFIGURABLE(cfgRefpTMax, float, 3.0f, "maximum pT for reference tracks if cfgUseRefpTCut is true")
  O2_DEFINE_CONFIGURABLE(cfgRefpTMin, float, 0.2f, "minimum pT for reference tracks if cfgUseRefpTCut is true")

  Configurable<std::vector<int>> cfgRunRemoveList{"cfgRunRemoveList", std::vector<int>{-1}, "excluded run numbers"};
  O2_DEFINE_CONFIGURABLE(cfgEvSelRCTflags, std::string, "", "keep empty to disable, usage: 'CentralBarrelTracking (CBT)', 'CBT_hadronPID' ")

  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgGlobalPtMin, float, 0.2f, "minimum accepted track pT")
    O2_DEFINE_CONFIGURABLE(cfgGlobalPtMax, float, 10.0f, "maximum accepted track pT")
    O2_DEFINE_CONFIGURABLE(cfgGlobalEta, float, 0.8f, "Eta cut")
    O2_DEFINE_CONFIGURABLE(cfgGlobalChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
    O2_DEFINE_CONFIGURABLE(cfgGlobalDCAz, float, 2.0f, "DCAz range for tracks")
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
    O2_DEFINE_CONFIGURABLE(cfgCentEstimator, int, 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FT0A")
    Configurable<LabeledArray<int>> cfgUseEventCuts{"cfgUseEventCuts", {LongArrayInt.front().data(), 14, 1, {"Filtered Events", "Sel8", "kNoTimeFrameBorder", "kNoITSROFrameBorder", "kNoSameBunchPileup", "kIsGoodZvtxFT0vsPV", "kNoCollInTimeRangeStandard", "kIsGoodITSLayersAll", "kIsGoodITSLayer0123", "kNoCollInRofStandard", "kNoHighMultCollInPrevRof", "Occupancy", "Multcorrelation", "T0AV0ACut"}, {"EvCuts"}}, "Labeled array (int) for various cuts on resonances"};
  } cfgEvSel;

  O2_DEFINE_CONFIGURABLE(cfgMinMixEventNum, int, 5, "Minimum number of events to mix")
  O2_DEFINE_CONFIGURABLE(cfgMergingCut, float, 0.02, "Merging cut on track merge")
  O2_DEFINE_CONFIGURABLE(cfgApplyTwoTrackEfficiency, bool, true, "Apply two track efficiency for tpc tpc")
  O2_DEFINE_CONFIGURABLE(cfgRadiusLow, float, 0.8, "Low radius for merging cut")
  O2_DEFINE_CONFIGURABLE(cfgRadiusHigh, float, 2.5, "High radius for merging cut")
  O2_DEFINE_CONFIGURABLE(cfgSampleSize, double, 10, "Sample size for mixed event")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgEfficiencyNch, std::string, "", "CCDB path to multiplicity dependent efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgCentralityWeight, std::string, "", "CCDB path to centrality weight object")
  O2_DEFINE_CONFIGURABLE(cfgLocalEfficiency, bool, false, "Use local efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgLocalEfficiencyNch, bool, false, "Use local multiplicity dependent efficiency object");
  O2_DEFINE_CONFIGURABLE(cfgUseEventWeights, bool, false, "Use event weights for mixed event")

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
    O2_DEFINE_CONFIGURABLE(cfgPIDParticle, int, 0, "1 = pion, 2 = kaon, 3 = proton, 4 = kshort, 5 = lambda, 6 = phi, 0 for no PID")
    O2_DEFINE_CONFIGURABLE(cfgGetNsigmaQA, bool, false, "Get QA histograms for selection of pions, kaons, and protons")
    O2_DEFINE_CONFIGURABLE(cfgGetdEdx, bool, false, "Get dEdx histograms for pions, kaons, and protons")
    O2_DEFINE_CONFIGURABLE(cfgPIDUseRejection, bool, true, "True: use exclusion exclusion criteria for PID determination, false: don't use exclusion")
    Configurable<LabeledArray<float>> nSigmas{"nSigmas", {LongArrayFloat.front().data(), 6, 3, {"UpCut_pi", "UpCut_ka", "UpCut_pr", "LowCut_pi", "LowCut_ka", "LowCut_pr"}, {"TPC", "TOF", "ITS"}}, "Labeled array for n-sigma values for TPC, TOF, ITS for pions, kaons, protons (positive and negative)"};
  } cfgPIDConfigs;

  Configurable<float> cfgCutFV0{"cfgCutFV0", 50., "FV0A threshold"};
  Configurable<float> cfgCutFT0A{"cfgCutFT0A", 150., "FT0A threshold"};
  Configurable<float> cfgCutFT0C{"cfgCutFT0C", 50., "FT0C threshold"};
  Configurable<float> cfgCutZDC{"cfgCutZDC", 10., "ZDC threshold"};

  SliceCache cache;

  ConfigurableAxis axisVertex{"axisVertex", {10, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisMult{"axisMult", {10, 0, 100}, "multiplicity axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {72, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisPtFiner{"axisPtFiner", {98, 0.2, 10.0}, "pt axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEtaTpcFt0a{"axisDeltaEtaTpcFt0a", {32, -5.8, -2.6}, "delta eta axis, -5.8~-2.6 for TPC-FT0A,"};
  ConfigurableAxis axisDeltaEtaTpcFt0c{"axisDeltaEtaTpcFt0c", {32, 1.2, 4.2}, "delta eta axis, 1.2~4.2 for TPC-FT0C"};
  ConfigurableAxis axisDeltaEtaFt0aFt0c{"axisDeltaEtaFt0aFt0c", {32, -1.5, 3.0}, "delta eta axis"};
  ConfigurableAxis axisDeltaEtaTpcTpc{"axisDeltaEtaTpcTpc", {32, -1.6, 1.6}, "delta eta axis for TPC-TPC"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt associated axis for histograms"};
  ConfigurableAxis axisVtxMix{"axisVtxMix", {VARIABLE_WIDTH, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, "vertex axis for mixed event histograms"};
  ConfigurableAxis axisMultMix{"axisMultMix", {VARIABLE_WIDTH, 0, 10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260}, "multiplicity / centrality axis for mixed event histograms"};
  ConfigurableAxis axisSample{"axisSample", {cfgSampleSize, 0, cfgSampleSize}, "sample axis for histograms"};
  ConfigurableAxis axisParticle{"axisParticle", {4, 0, 4}, "particle axis for correlation container"};

  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt axis for efficiency histograms"};
  ConfigurableAxis axisAmplitudeFt0a{"axisAmplitudeFt0a", {5000, 0, 1000}, "FT0A amplitude"};
  ConfigurableAxis axisChannelFt0aAxis{"axisChannelFt0aAxis", {96, 0.0, 96.0}, "FT0A channel"};

  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {80, -5, 5}, "nsigmaTPC axis"};
  ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {80, -5, 5}, "nsigmaTOF axis"};
  ConfigurableAxis axisNsigmaITS{"axisNsigmaITS", {80, -5, 5}, "nsigmaITS axis"};
  ConfigurableAxis axisTpcSignal{"axisTpcSignal", {250, 0, 250}, "dEdx axis for TPC"};
  ConfigurableAxis axisSigma{"axisSigma", {200, 0, 20}, "sigma axis for TPC"};

  Configurable<std::string> cfgGainEqPath{"cfgGainEqPath", "Analysis/EventPlane/GainEq", "CCDB path for gain equalization constants"};
  Configurable<int> cfgCorrLevel{"cfgCorrLevel", 0, "calibration step: 0 = no corr, 1 = gain corr"};
  ConfigurableAxis cfgaxisFITamp{"cfgaxisFITamp", {1000, 0, 5000}, ""};
  AxisSpec axisFit{cfgaxisFITamp, "fit amplitude"};
  AxisSpec axisChID = {220, 0, 220};

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
  OutputObj<CorrelationContainer> sameTpcFt0a{"sameEvent_TPC_FT0A"};
  OutputObj<CorrelationContainer> mixedTpcFt0a{"mixedEvent_TPC_FT0A"};
  OutputObj<CorrelationContainer> sameTpcFt0c{"sameEvent_TPC_FT0C"};
  OutputObj<CorrelationContainer> mixedTpcFt0c{"mixedEvent_TPC_FT0C"};
  OutputObj<CorrelationContainer> sameFt0aFt0c{"sameEvent_FT0A_FT0C"};
  OutputObj<CorrelationContainer> mixedFt0aFt0c{"mixedEvent_FT0A_FT0C"};
  OutputObj<CorrelationContainer> sameTPC{"sameEvent_TPC"};
  OutputObj<CorrelationContainer> mixedTPC{"mixedEvent_TPC"};

  HistogramRegistry registry{"registry"};
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

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    LOGF(info, "Starting init");
    readMatrix(cfgTrkSel.cfgTrackCuts->getData(), trackCuts);
    readMatrix(cfgEvSel.cfgUseEventCuts->getData(), eventCuts);
    readMatrix(cfgPIDConfigs.nSigmas->getData(), nSigmaVals);
    if (doprocessSameFt0aFt0c || doprocessSameTpcFt0a || doprocessSameTpcFt0c || doprocessSameTPC) {
      registry.add("hEventCountRct", "Number of Event;; Count", {HistType::kTH1D, {{2, 0, 2}}});
      registry.get<TH1>(HIST("hEventCountRct"))->GetXaxis()->SetBinLabel(1, "rct fail");
      registry.get<TH1>(HIST("hEventCountRct"))->GetXaxis()->SetBinLabel(2, "rct pass");
      registry.add("hEventCount", "Number of Event;; Count", {HistType::kTH1D, {{14, -0.5, 13.5}}});
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(FilteredEvents + 1, "Filtered events");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(AfterSel8 + 1, "After sel8");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseNoTimeFrameBorder + 1, "kNoTimeFrameBorder");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseNoITSROFrameBorder + 1, "kNoITSROFrameBorder");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseNoSameBunchPileup + 1, "kNoSameBunchPileup");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseGoodZvtxFT0vsPV + 1, "kIsGoodZvtxFT0vsPV");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseNoCollInTimeRangeStandard + 1, "kNoCollInTimeRangeStandard");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseGoodITSLayersAll + 1, "kIsGoodITSLayersAll");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseGoodITSLayer0123 + 1, "kIsGoodITSLayer0123");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseNoCollInRofStandard + 1, "kNoCollInRofStandard");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseNoHighMultCollInPrevRof + 1, "kNoHighMultCollInPrevRof");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseOccupancy + 1, "Occupancy Cut");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseMultCorrCut + 1, "MultCorrelation Cut");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(UseT0AV0ACut + 1, "T0AV0A cut");
    }

    if ((doprocessSameFt0aFt0c || doprocessSameTpcFt0a || doprocessSameTpcFt0c || doprocessSameTPC) && cfgQaCheck) {
      registry.add("hPassedEventSelection", "Number of Event;; Count", {HistType::kTH1D, {{12, -0.5, 11.5}}});
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(FilteredEvents + 1, "Filtered events");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(AfterSel8 + 1, "After sel8");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(UseNoTimeFrameBorder + 1, "kNoTimeFrameBorder");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(UseNoITSROFrameBorder + 1, "kNoITSROFrameBorder");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(UseNoSameBunchPileup + 1, "kNoSameBunchPileup");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(UseGoodZvtxFT0vsPV + 1, "kIsGoodZvtxFT0vsPV");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(UseNoCollInTimeRangeStandard + 1, "kNoCollInTimeRangeStandard");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(UseGoodITSLayersAll + 1, "kIsGoodITSLayersAll");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(UseGoodITSLayer0123 + 1, "kIsGoodITSLayer0123");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(UseNoCollInRofStandard + 1, "kNoCollInRofStandard");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(UseNoHighMultCollInPrevRof + 1, "kNoHighMultCollInPrevRof");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(UseOccupancy + 1, "Occupancy Cut");
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

    if (doprocessSameTPC || doprocessSameFt0aFt0c || doprocessSameTpcFt0a || doprocessSameTpcFt0c) {
      registry.add("Phi", "Phi", {HistType::kTH1D, {axisPhi}});
      registry.add("Eta", "Eta", {HistType::kTH1D, {axisEta}});
      registry.add("EtaCorrected", "EtaCorrected", {HistType::kTH1D, {axisEta}});
      registry.add("hParticleCounts", "hParticleCounts", {HistType::kTH1D, {axisParticle}});
      registry.add("hParticleSelected", "hParticleSelected", {HistType::kTH2D, {axisParticle, axisPtTrigger}});
      registry.add("pTFiner", "pTFiner", {HistType::kTH1D, {axisPtFiner}});
      registry.add("pTFinerCorrected", "pTFinerCorrected", {HistType::kTH1D, {axisPtFiner}});
      registry.add("Nch", "N_{ch}", {HistType::kTH1D, {axisMult}});
      registry.add("zVtx", "zVtx", {HistType::kTH1D, {axisVertex}});
      if (doprocessSameFt0aFt0c || doprocessSameTpcFt0a || doprocessSameTpcFt0c) {
        registry.add("FT0Amp", "", {HistType::kTH2F, {axisChID, axisFit}});
        registry.add("FT0AmpCorrect", "", {HistType::kTH2F, {axisChID, axisFit}});
      }
      if (cfgPIDConfigs.cfgGetNsigmaQA) {
        if (!cfgPIDConfigs.cfgUseItsPID) {
          registry.add("TofTpcNsigma_before", "", {HistType::kTHnSparseD, {{axisNsigmaTPC, axisNsigmaTOF, axisPtTrigger}}});
          registry.add("TofTpcNsigma_after", "", {HistType::kTHnSparseD, {{axisNsigmaTPC, axisNsigmaTOF, axisPtTrigger}}});
        } // TPC-TOF PID QA hists
        if (cfgPIDConfigs.cfgUseItsPID) {
          registry.add("TofItsNsigma_before", "", {HistType::kTHnSparseD, {{axisNsigmaITS, axisNsigmaTOF, axisPtTrigger}}});
          registry.add("TofItsNsigma_after", "", {HistType::kTHnSparseD, {{axisNsigmaITS, axisNsigmaTOF, axisPtTrigger}}});
        } // ITS-TOF PID QA hists
      } // end of PID QA hists

      if (cfgPIDConfigs.cfgGetdEdx) {
        registry.add("TpcdEdx_ptwise_beforeCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF}}});
        registry.add("ExpTpcdEdx_ptwise_beforeCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF}}});
        registry.add("ExpSigma_ptwise_beforeCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisSigma, axisNsigmaTOF}}});

        registry.add("TpcdEdx_ptwise_afterCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF}}});
        registry.add("ExpTpcdEdx_ptwise_afterCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF}}});
        registry.add("ExpSigma_ptwise_afterCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisSigma, axisNsigmaTOF}}});
      }

    } // end of single particle distribution histograms

    if (cfgQaCheck) {
      registry.add("Nch_corrected", "N_{ch} corrected", {HistType::kTH1D, {axisMult}});
    }

    if (doprocessSameTpcFt0a) {
      registry.add("deltaEta_deltaPhi_same_TPC_FT0A", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0a}}); // check to see the delta eta and delta phi distribution
      registry.add("deltaEta_deltaPhi_mixed_TPC_FT0A", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0a}});
      registry.add("Assoc_amp_same_TPC_FT0A", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
      registry.add("Assoc_amp_mixed_TPC_FT0A", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
      registry.add("Trig_hist_TPC_FT0A", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger, axisMult}}});
    }
    if (doprocessSameTpcFt0c) {
      registry.add("deltaEta_deltaPhi_same_TPC_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0c}}); // check to see the delta eta and delta phi distribution
      registry.add("deltaEta_deltaPhi_mixed_TPC_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0c}});
      registry.add("Assoc_amp_same_TPC_FT0C", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
      registry.add("Assoc_amp_mixed_TPC_FT0C", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
      registry.add("Trig_hist_TPC_FT0C", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger, axisMult}}});
    }
    if (doprocessSameFt0aFt0c) {
      registry.add("deltaEta_deltaPhi_same_FT0A_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaFt0aFt0c}}); // check to see the delta eta and delta phi distribution
      registry.add("deltaEta_deltaPhi_mixed_FT0A_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaFt0aFt0c}});
      registry.add("Trig_hist_FT0A_FT0C", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger, axisMult}}});
    }
    if (doprocessSameTPC) {
      registry.add("deltaEta_deltaPhi_same_TPC", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcTpc}}); // check to see the delta eta and delta phi distribution
      registry.add("deltaEta_deltaPhi_mixed_TPC", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcTpc}});
      registry.add("Trig_hist_TPC", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger, axisMult}}});
    }

    registry.add("eventcount", "bin", {HistType::kTH1F, {{4, 0, 4, "bin"}}}); // histogram to see how many events are in the same and mixed event

    LOGF(info, "Initializing correlation container");

    // Initialize Nch-related histograms and containers

    std::vector<AxisSpec> corrAxisTpcFt0a = {{axisSample, "Sample"},
                                             {axisVertex, "z-vtx (cm)"},
                                             {axisPtTrigger, "p_{T} (GeV/c)"},
                                             {axisMult, "N_{ch}"},
                                             {axisDeltaPhi, "#Delta#varphi (rad)"},
                                             {axisDeltaEtaTpcFt0a, "#Delta#eta"}};
    std::vector<AxisSpec> effAxis = {
      {axisEtaEfficiency, "#eta"},
      {axisPtEfficiency, "p_{T} (GeV/c)"},
      {axisVertexEfficiency, "z-vtx (cm)"},
    };
    std::vector<AxisSpec> userAxis;

    std::vector<AxisSpec> corrAxisTpcFt0c = {{axisSample, "Sample"},
                                             {axisVertex, "z-vtx (cm)"},
                                             {axisPtTrigger, "p_{T} (GeV/c)"},
                                             {axisMult, "N_{ch}"},
                                             {axisDeltaPhi, "#Delta#varphi (rad)"},
                                             {axisDeltaEtaTpcFt0c, "#Delta#eta"}};

    std::vector<AxisSpec> corrAxisFt0aFt0c = {{axisSample, "Sample"},
                                              {axisVertex, "z-vtx (cm)"},
                                              {axisPtTrigger, "p_{T} (GeV/c)"},
                                              {axisMult, "N_{ch}"},
                                              {axisDeltaPhi, "#Delta#varphi (rad)"},
                                              {axisDeltaEtaFt0aFt0c, "#Delta#eta"}};

    std::vector<AxisSpec> corrAxisTPC = {{axisSample, "Sample"},
                                         {axisVertex, "z-vtx (cm)"},
                                         {axisPtTrigger, "p_{T} (GeV/c)"},
                                         {axisMult, "N_{ch}"},
                                         {axisDeltaPhi, "#Delta#varphi (rad)"},
                                         {axisDeltaEtaTpcTpc, "#Delta#eta"}}; // use the same delta eta axis for TPC-TPC correlation

    if (doprocessSameTpcFt0a) {
      sameTpcFt0a.setObject(new CorrelationContainer("sameEvent_TPC_FT0A", "sameEvent_TPC_FT0A", corrAxisTpcFt0a, effAxis, userAxis));
      mixedTpcFt0a.setObject(new CorrelationContainer("mixedEvent_TPC_FT0A", "mixedEvent_TPC_FT0A", corrAxisTpcFt0a, effAxis, userAxis));
    }
    if (doprocessSameTpcFt0c) {
      sameTpcFt0c.setObject(new CorrelationContainer("sameEvent_TPC_FT0C", "sameEvent_TPC_FT0C", corrAxisTpcFt0c, effAxis, userAxis));
      mixedTpcFt0c.setObject(new CorrelationContainer("mixedEvent_TPC_FT0C", "mixedEvent_TPC_FT0C", corrAxisTpcFt0c, effAxis, userAxis));
    }
    if (doprocessSameFt0aFt0c) {
      sameFt0aFt0c.setObject(new CorrelationContainer("sameEvent_FT0A_FT0C", "sameEvent_FT0A_FT0C", corrAxisFt0aFt0c, effAxis, userAxis));
      mixedFt0aFt0c.setObject(new CorrelationContainer("mixedEvent_FT0A_FT0C", "mixedEvent_FT0A_FT0C", corrAxisFt0aFt0c, effAxis, userAxis));
    }
    if (doprocessSameTPC) {
      sameTPC.setObject(new CorrelationContainer("sameEvent_TPC", "sameEvent_TPC", corrAxisTPC, effAxis, userAxis));
      mixedTPC.setObject(new CorrelationContainer("mixedEvent_TPC", "mixedEvent_TPC", corrAxisTPC, effAxis, userAxis));
    }

    if (!cfgEvSelRCTflags.value.empty()) {
      rctChecker.init(cfgEvSelRCTflags.value); // override initialzation
    }

    LOGF(info, "End of init");
  }

  template <typename TCollision>
  double getCentrality(TCollision const& collision)
  {
    double cent = 0.0;
    switch (cfgEvSel.cfgCentEstimator) {
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
        registry.fill(HIST("hEventCountRct"), 0.5);
      }

      return false;
    }
    if (fillCounter) {
      registry.fill(HIST("hEventCountRct"), 1.5);
    }

    return true;
  }

  template <typename TCollision>
  bool eventSelected(TCollision const& collision, const int mult, const bool fillCounter)
  {
    if (eventCuts[UseNoTimeFrameBorder][EvCut1] && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (fillCounter && eventCuts[UseNoTimeFrameBorder][EvCut1]) {
      registry.fill(HIST("hEventCount"), UseNoTimeFrameBorder);
    }

    if (eventCuts[UseNoITSROFrameBorder][EvCut1] && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (fillCounter && eventCuts[UseNoITSROFrameBorder][EvCut1]) {
      registry.fill(HIST("hEventCount"), UseNoITSROFrameBorder);
    }

    if (eventCuts[UseNoSameBunchPileup][EvCut1] && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return false;
    }
    if (fillCounter && eventCuts[UseNoSameBunchPileup][EvCut1]) {
      registry.fill(HIST("hEventCount"), UseNoSameBunchPileup);
    }

    if (eventCuts[UseGoodZvtxFT0vsPV][EvCut1] && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return false;
    }
    if (fillCounter && eventCuts[UseGoodZvtxFT0vsPV][EvCut1]) {
      registry.fill(HIST("hEventCount"), UseGoodZvtxFT0vsPV);
    }

    if (eventCuts[UseNoCollInTimeRangeStandard][EvCut1] && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      // no collisions in specified time range
      return false;
    }

    if (fillCounter && eventCuts[UseNoCollInTimeRangeStandard][EvCut1]) {
      registry.fill(HIST("hEventCount"), UseNoCollInTimeRangeStandard);
    }

    if (eventCuts[UseGoodITSLayersAll][EvCut1] && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      // from Jan 9 2025 AOT meeting
      // cut time intervals with dead ITS staves
      return false;
    }

    if (fillCounter && eventCuts[UseGoodITSLayersAll][EvCut1]) {
      registry.fill(HIST("hEventCount"), UseGoodITSLayersAll);
    }

    if (eventCuts[UseGoodITSLayer0123][EvCut1] && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
      return false;
    }
    if (fillCounter && eventCuts[UseGoodITSLayer0123][EvCut1]) {
      registry.fill(HIST("hEventCount"), UseGoodITSLayer0123);
    }

    if (eventCuts[UseNoCollInRofStandard][EvCut1] && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      // no other collisions in this Readout Frame with per-collision multiplicity above threshold
      return false;
    }

    if (fillCounter && eventCuts[UseNoCollInRofStandard][EvCut1]) {
      registry.fill(HIST("hEventCount"), UseNoCollInRofStandard);
    }

    if (eventCuts[UseNoHighMultCollInPrevRof][EvCut1] && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      // veto an event if FT0C amplitude in previous ITS ROF is above threshold
      return false;
    }
    if (fillCounter && eventCuts[UseNoHighMultCollInPrevRof][EvCut1]) {
      registry.fill(HIST("hEventCount"), UseNoHighMultCollInPrevRof);
    }

    auto multNTracksPV = collision.multNTracksPV();
    auto occupancy = collision.trackOccupancyInTimeRange();

    if (eventCuts[UseOccupancy][EvCut1] && (occupancy < cfgEvSel.cfgCutOccupancyLow || occupancy > cfgEvSel.cfgCutOccupancyHigh)) {
      return false;
    }
    if (fillCounter && eventCuts[UseOccupancy][EvCut1]) {
      registry.fill(HIST("hEventCount"), UseOccupancy);
    }

    if (eventCuts[UseMultCorrCut][EvCut1]) {
      float cent = getCentrality(collision);
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
      registry.fill(HIST("hEventCount"), UseMultCorrCut);
    }

    // V0A T0A 5 sigma cut
    if (eventCuts[UseT0AV0ACut][EvCut1] && (std::fabs(collision.multFV0A() - cfgFuncParas.fT0AV0AMean->Eval(collision.multFT0A())) > cfgFuncParas.cfgV0AT0Acut * cfgFuncParas.fT0AV0ASigma->Eval(collision.multFT0A()))) {
      return false;
    }
    if (fillCounter && eventCuts[UseT0AV0ACut][EvCut1]) {
      registry.fill(HIST("hEventCount"), UseT0AV0ACut);
    }

    return true;
  }

  template <typename TCollision>
  void eventSelectedIndividually(TCollision const& collision)
  {
    registry.fill(HIST("hPassedEventSelection"), FilteredEvents);

    if (collision.sel8()) {
      registry.fill(HIST("hPassedEventSelection"), AfterSel8);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      registry.fill(HIST("hPassedEventSelection"), UseNoTimeFrameBorder);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      registry.fill(HIST("hPassedEventSelection"), UseNoITSROFrameBorder);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      registry.fill(HIST("hPassedEventSelection"), UseNoSameBunchPileup);
    }

    if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      registry.fill(HIST("hPassedEventSelection"), UseGoodZvtxFT0vsPV);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      registry.fill(HIST("hPassedEventSelection"), UseNoCollInTimeRangeStandard);
    }

    if (collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      registry.fill(HIST("hPassedEventSelection"), UseGoodITSLayersAll);
    }

    if (collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
      registry.fill(HIST("hPassedEventSelection"), UseGoodITSLayer0123);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      registry.fill(HIST("hPassedEventSelection"), UseNoCollInRofStandard);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      registry.fill(HIST("hPassedEventSelection"), UseNoHighMultCollInPrevRof);
    }

    auto occupancy = collision.trackOccupancyInTimeRange();
    if (eventCuts[UseOccupancy][EvCut1] && (occupancy > cfgEvSel.cfgCutOccupancyLow && occupancy < cfgEvSel.cfgCutOccupancyHigh)) {
      registry.fill(HIST("hPassedEventSelection"), UseOccupancy);
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

    return ((std::abs(track.eta()) < trackCuts[TrkCutEtaCut][UseGenTrkCuts]) && (track.tpcChi2NCl() < trackCuts[TrkCutChi2PrTpcCls][UseGenTrkCuts]) && (track.tpcNClsFound() >= trackCuts[TrkCutTpcCluster][UseGenTrkCuts]) && (track.tpcNClsCrossedRows() >= trackCuts[TrkCutTpcCrossedRows][UseGenTrkCuts]) && (track.itsNCls() >= trackCuts[TrkCutItsCluster][UseGenTrkCuts]));
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

  template <typename TTrack, typename TTrackAssoc>
  float getDPhiStar(TTrack const& track1, TTrackAssoc const& track2, float radius, int magField)
  {
    float charge1 = track1.sign();
    float charge2 = track2.sign();

    float phi1 = track1.phi();
    float phi2 = track2.phi();

    float pt1 = track1.pt();
    float pt2 = track2.pt();

    int fbSign = (magField > 0) ? 1 : -1;

    float dPhiStar = phi1 - phi2 - charge1 * fbSign * std::asin(0.075 * radius / pt1) + charge2 * fbSign * std::asin(0.075 * radius / pt2);

    if (dPhiStar > constants::math::PI) {
      dPhiStar = constants::math::TwoPI - dPhiStar;
    }
    if (dPhiStar < -constants::math::PI) {
      dPhiStar = -constants::math::TwoPI - dPhiStar;
    }

    return dPhiStar;
  }

  void loadGain(aod::BCsWithTimestamps::iterator const& bc)
  {
    cstFT0RelGain.clear();
    std::string fullPath;

    auto timestamp = bc.timestamp();
    constexpr int ChannelsFT0 = 208;
    if (cfgCorrLevel == 0) {
      for (auto i{0u}; i < ChannelsFT0; i++) {
        cstFT0RelGain.push_back(1.);
      }
    } else {
      fullPath = cfgGainEqPath;
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
        registry.fill(HIST("FT0Amp"), id, ampl);
      }
      ampl = ampl / cstFT0RelGain[id];
      if (system == SameEvent) {
        registry.fill(HIST("FT0AmpCorrect"), id, ampl);
      }
    } else if (fitType == IndexFT0A) {
      id = ft0.channelA()[iCh];
      ampl = ft0.amplitudeA()[iCh];
      if (system == SameEvent) {
        registry.fill(HIST("FT0Amp"), id, ampl);
      }
      ampl = ampl / cstFT0RelGain[id];
      if (system == SameEvent) {
        registry.fill(HIST("FT0AmpCorrect"), id, ampl);
      }
    } else {
      LOGF(fatal, "Cor Index %d out of range", fitType);
    }
  }

  template <typename TTrack>
  int getNsigmaPID(const TTrack& track, bool fillYields)
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

    if (fillYields) {
      registry.fill(HIST("hParticleCounts"), UseCharged);
    }

    if (isPion) {
      pid = UsePions;
      if (fillYields) {
        registry.fill(HIST("hParticleCounts"), UsePions);
      }
    } else if (isKaon) {
      pid = UseKaons;
      if (fillYields) {
        registry.fill(HIST("hParticleCounts"), UseKaons);
      }
    } else if (isProton) {
      pid = UseProtons;
      if (fillYields) {
        registry.fill(HIST("hParticleCounts"), UseProtons);
      }
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
    if (!cfgEfficiency.value.empty()) {
      if (cfgLocalEfficiency) {
        TFile* fEfficiencyTrigger = TFile::Open(cfgEfficiency.value.c_str(), "READ");
        mEfficiency = dynamic_cast<TH3D*>(fEfficiencyTrigger->Get("ccdb_object"));

      } else {
        mEfficiency = ccdb->getForTimeStamp<TH3D>(cfgEfficiency, timestamp);
      }
      if (mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEfficiency.value.c_str(), static_cast<void*>(mEfficiency));
    }
    if (!cfgEfficiencyNch.value.empty()) {
      if (cfgLocalEfficiencyNch) {
        TFile* fEfficiencyTrigger = TFile::Open(cfgEfficiencyNch.value.c_str(), "READ");
        mEfficiencyNch = dynamic_cast<TH1D*>(fEfficiencyTrigger->Get("ccdb_object"));

      } else {
        mEfficiencyNch = ccdb->getForTimeStamp<TH1D>(cfgEfficiencyNch, timestamp);
      }
      if (!mEfficiencyNch) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiencyNch.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEfficiencyNch.value.c_str(), static_cast<void*>(mEfficiencyNch));
    }
    if (!cfgCentralityWeight.value.empty()) {
      mCentralityWeight = ccdb->getForTimeStamp<TH1D>(cfgCentralityWeight, timestamp);
      if (mCentralityWeight == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgCentralityWeight.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgCentralityWeight.value.c_str(), static_cast<void*>(mCentralityWeight));
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

  int getMagneticField(uint64_t timestamp)
  {
    // Get the magnetic field
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("/GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  template <typename TTracks>
  void trackCounter(const TTracks& tracks, double& multiplicity) // function to count the number of tracks in the event and fill the histogram
  {
    double nTracksCorrected = 0;
    float weightNch = 1.0f;
    for (auto const& track : tracks) {
      if (!trackSelectedForNch(track)) {
        continue;
      }

      if (!getEfficiencyCorrectionNch(weightNch, track.pt())) {
        continue;
      }

      nTracksCorrected += weightNch;
    }
    multiplicity = nTracksCorrected;
  }

  template <typename TCollision, typename TTracks>
  void fillYield(const TCollision& collision, const TTracks& tracks) // function to fill the yield and etaphi histograms.
  {

    float weff1 = 1.0;
    float zvtx = collision.posZ();

    for (auto const& track1 : tracks) {

      if (!trackSelected(track1)) {
        continue;
      }
      if (cfgPIDConfigs.cfgPIDParticle && getNsigmaPID(track1, true) != cfgPIDConfigs.cfgPIDParticle) {
        continue; // if PID is selected, check if the track has the right PID
      }

      if (!getEfficiencyCorrection(weff1, track1.eta(), track1.pt(), zvtx)) {
        continue;
      }

      registry.fill(HIST("Phi"), RecoDecay::constrainAngle(track1.phi(), 0.0));
      registry.fill(HIST("Eta"), track1.eta());
      registry.fill(HIST("EtaCorrected"), track1.eta(), weff1);
      registry.fill(HIST("pTFiner"), track1.pt());
      registry.fill(HIST("pTFinerCorrected"), track1.pt(), weff1);
      registry.fill(HIST("hParticleSelected"), cfgPIDConfigs.cfgPIDParticle.value, track1.pt());
    }
  }

  template <typename TTrack>
  void fillNsigmaBeforeCut(const TTrack& track1, int pid) // function to fill the QA before Nsigma selection
  {
    switch (pid) {
      case UsePions: // For Pions
        if (!cfgPIDConfigs.cfgUseItsPID) {
          if (cfgPIDConfigs.cfgGetNsigmaQA) {
            registry.fill(HIST("TofTpcNsigma_before"), track1.tpcNSigmaPi(), track1.tofNSigmaPi(), track1.pt());
          }
          if (cfgPIDConfigs.cfgGetdEdx) {
            double tpcExpSignalPi = track1.tpcSignal() - (track1.tpcNSigmaPi() * track1.tpcExpSigmaPi());

            registry.fill(HIST("TpcdEdx_ptwise_beforeCut"), track1.pt(), track1.tpcSignal(), track1.tofNSigmaPi());
            registry.fill(HIST("ExpTpcdEdx_ptwise_beforeCut"), track1.pt(), tpcExpSignalPi, track1.tofNSigmaPi());
            registry.fill(HIST("ExpSigma_ptwise_beforeCut"), track1.pt(), track1.tpcExpSigmaPi(), track1.tofNSigmaPi());
          }
        }
        if (cfgPIDConfigs.cfgGetNsigmaQA && cfgPIDConfigs.cfgUseItsPID) {
          registry.fill(HIST("TofItsNsigma_before"), itsResponse.nSigmaITS<o2::track::PID::Pion>(track1), track1.tofNSigmaPi(), track1.pt());
        }
        break;
      case UseKaons: // For Kaons
        if (!cfgPIDConfigs.cfgUseItsPID) {
          if (cfgPIDConfigs.cfgGetNsigmaQA) {
            registry.fill(HIST("TofTpcNsigma_before"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
          }
          if (cfgPIDConfigs.cfgGetdEdx) {
            double tpcExpSignalKa = track1.tpcSignal() - (track1.tpcNSigmaKa() * track1.tpcExpSigmaKa());

            registry.fill(HIST("TpcdEdx_ptwise_beforeCut"), track1.pt(), track1.tpcSignal(), track1.tofNSigmaKa());
            registry.fill(HIST("ExpTpcdEdx_ptwise_beforeCut"), track1.pt(), tpcExpSignalKa, track1.tofNSigmaKa());
            registry.fill(HIST("ExpSigma_ptwise_beforeCut"), track1.pt(), track1.tpcExpSigmaKa(), track1.tofNSigmaKa());
          }
        }
        if (cfgPIDConfigs.cfgGetNsigmaQA && cfgPIDConfigs.cfgUseItsPID) {
          registry.fill(HIST("TofItsNsigma_before"), itsResponse.nSigmaITS<o2::track::PID::Kaon>(track1), track1.tofNSigmaKa(), track1.pt());
        }
        break;

      case UseProtons: // For Protons
        if (!cfgPIDConfigs.cfgUseItsPID) {
          if (cfgPIDConfigs.cfgGetNsigmaQA) {
            registry.fill(HIST("TofTpcNsigma_before"), track1.tpcNSigmaPr(), track1.tofNSigmaPr(), track1.pt());
          }
          if (cfgPIDConfigs.cfgGetdEdx) {
            double tpcExpSignalPr = track1.tpcSignal() - (track1.tpcNSigmaPr() * track1.tpcExpSigmaPr());

            registry.fill(HIST("TpcdEdx_ptwise_beforeCut"), track1.pt(), track1.tpcSignal(), track1.tofNSigmaPr());
            registry.fill(HIST("ExpTpcdEdx_ptwise_beforeCut"), track1.pt(), tpcExpSignalPr, track1.tofNSigmaPr());
            registry.fill(HIST("ExpSigma_ptwise_beforeCut"), track1.pt(), track1.tpcExpSigmaPr(), track1.tofNSigmaPr());
          }
        }
        if (cfgPIDConfigs.cfgGetNsigmaQA && cfgPIDConfigs.cfgUseItsPID) {
          registry.fill(HIST("TofItsNsigma_before"), itsResponse.nSigmaITS<o2::track::PID::Proton>(track1), track1.tofNSigmaPr(), track1.pt());
        }
        break;
      default:
        LOGF(info, "PID %d not recognized, skipping nSigma QA filling", pid);
        break;
    } // end of switch
  } // end of fillNsigmaBeforeCut

  template <typename TTrack>
  void fillNsigmaAfterCut(const TTrack& track1, int pid) // function to fill the QA after Nsigma selection
  {
    switch (pid) {
      case UsePions: // For Pions
        if (!cfgPIDConfigs.cfgUseItsPID) {
          if (cfgPIDConfigs.cfgGetdEdx) {
            double tpcExpSignalPi = track1.tpcSignal() - (track1.tpcNSigmaPi() * track1.tpcExpSigmaPi());

            registry.fill(HIST("TpcdEdx_ptwise_afterCut"), track1.pt(), track1.tpcSignal(), track1.tofNSigmaPi());
            registry.fill(HIST("ExpTpcdEdx_ptwise_afterCut"), track1.pt(), tpcExpSignalPi, track1.tofNSigmaPi());
            registry.fill(HIST("ExpSigma_ptwise_afterCut"), track1.pt(), track1.tpcExpSigmaPi(), track1.tofNSigmaPi());
          }
          if (cfgPIDConfigs.cfgGetNsigmaQA) {
            registry.fill(HIST("TofTpcNsigma_after"), track1.tpcNSigmaPi(), track1.tofNSigmaPi(), track1.pt());
          }
        }
        if (cfgPIDConfigs.cfgUseItsPID) {
          registry.fill(HIST("TofItsNsigma_after"), itsResponse.nSigmaITS<o2::track::PID::Pion>(track1), track1.tofNSigmaPi(), track1.pt());
        }
        break;
      case UseKaons: // For Kaons
        if (!cfgPIDConfigs.cfgUseItsPID) {
          if (cfgPIDConfigs.cfgGetdEdx) {
            double tpcExpSignalKa = track1.tpcSignal() - (track1.tpcNSigmaKa() * track1.tpcExpSigmaKa());

            registry.fill(HIST("TpcdEdx_ptwise_afterCut"), track1.pt(), track1.tpcSignal(), track1.tofNSigmaKa());
            registry.fill(HIST("ExpTpcdEdx_ptwise_afterCut"), track1.pt(), tpcExpSignalKa, track1.tofNSigmaKa());
            registry.fill(HIST("ExpSigma_ptwise_afterCut"), track1.pt(), track1.tpcExpSigmaKa(), track1.tofNSigmaKa());
          }
          if (cfgPIDConfigs.cfgGetNsigmaQA) {
            registry.fill(HIST("TofTpcNsigma_after"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
          }
        }
        if (cfgPIDConfigs.cfgUseItsPID) {
          registry.fill(HIST("TofItsNsigma_after"), itsResponse.nSigmaITS<o2::track::PID::Kaon>(track1), track1.tofNSigmaKa(), track1.pt());
        }
        break;
      case UseProtons: // For Protons
        if (!cfgPIDConfigs.cfgUseItsPID) {
          if (cfgPIDConfigs.cfgGetdEdx) {
            double tpcExpSignalPr = track1.tpcSignal() - (track1.tpcNSigmaPr() * track1.tpcExpSigmaPr());

            registry.fill(HIST("TpcdEdx_ptwise_afterCut"), track1.pt(), track1.tpcSignal(), track1.tofNSigmaPr());
            registry.fill(HIST("ExpTpcdEdx_ptwise_afterCut"), track1.pt(), tpcExpSignalPr, track1.tofNSigmaPr());
            registry.fill(HIST("ExpSigma_ptwise_afterCut"), track1.pt(), track1.tpcExpSigmaPr(), track1.tofNSigmaPr());
          }
          if (cfgPIDConfigs.cfgGetNsigmaQA) {
            registry.fill(HIST("TofTpcNsigma_after"), track1.tpcNSigmaPr(), track1.tofNSigmaPr(), track1.pt());
          }
        }
        if (cfgPIDConfigs.cfgUseItsPID && cfgPIDConfigs.cfgGetNsigmaQA) {
          registry.fill(HIST("TofItsNsigma_after"), itsResponse.nSigmaITS<o2::track::PID::Proton>(track1), track1.tofNSigmaPr(), track1.pt());
        }
        break;
      default:
        LOGF(info, "PID %d not recognized, skipping nSigma QA filling, breaking", pid);
        break;
    } // end of switch
  }

  template <CorrelationContainer::CFStep step, typename TTracks, typename TFT0s>
  void fillCorrelationsTPCFT0(const TTracks& tracks1, TFT0s const& ft0, float posZ, int system, double multiplicity, int corType, float eventWeight) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {
    auto fSampleIndex = static_cast<int>(fRandom.Uniform(0.0, cfgSampleSize));

    float triggerWeight = 1.0f;
    // loop over all tracks
    for (auto const& track1 : tracks1) {

      if (!trackSelected(track1)) {
        continue;
      }
      if (cfgPIDConfigs.cfgPIDParticle > UseCharged) {
        fillNsigmaBeforeCut(track1, cfgPIDConfigs.cfgPIDParticle);
      }

      if (cfgPIDConfigs.cfgPIDParticle && getNsigmaPID(track1, false) != cfgPIDConfigs.cfgPIDParticle) {
        continue; // if PID is selected, check if the track has the right PID
      }

      if (cfgPIDConfigs.cfgPIDParticle > UseCharged) {
        fillNsigmaAfterCut(track1, cfgPIDConfigs.cfgPIDParticle);
      }

      if (!getEfficiencyCorrection(triggerWeight, track1.pt(), track1.eta(), posZ)) {
        continue;
      }

      if (system == SameEvent) {
        if (corType == IndexFT0C) {
          registry.fill(HIST("Trig_hist_TPC_FT0C"), fSampleIndex, posZ, track1.pt(), multiplicity, eventWeight * triggerWeight);
        } else if (corType == IndexFT0A) {
          registry.fill(HIST("Trig_hist_TPC_FT0A"), fSampleIndex, posZ, track1.pt(), multiplicity, eventWeight * triggerWeight);
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
            if (cfgQaCheck) {
              registry.fill(HIST("Assoc_amp_same_TPC_FT0A"), chanelid, ampl);
              registry.fill(HIST("deltaEta_deltaPhi_same_TPC_FT0A"), deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
            }
            sameTpcFt0a->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), multiplicity, deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
          } else if (corType == IndexFT0C) {
            if (cfgQaCheck) {
              registry.fill(HIST("Assoc_amp_same_TPC_FT0C"), chanelid, ampl);
              registry.fill(HIST("deltaEta_deltaPhi_same_TPC_FT0C"), deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
            }
            sameTpcFt0c->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), multiplicity, deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
          }
        } else if (system == MixedEvent) {
          if (corType == IndexFT0A) {
            if (cfgQaCheck) {
              registry.fill(HIST("Assoc_amp_mixed_TPC_FT0A"), chanelid, ampl);
              registry.fill(HIST("deltaEta_deltaPhi_mixed_TPC_FT0A"), deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
            }
            mixedTpcFt0a->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), multiplicity, deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
          } else if (corType == IndexFT0C) {
            if (cfgQaCheck) {
              registry.fill(HIST("Assoc_amp_mixed_TPC_FT0C"), chanelid, ampl);
              registry.fill(HIST("deltaEta_deltaPhi_mixed_TPC_FT0C"), deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
            }
            mixedTpcFt0c->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), multiplicity, deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
          }
        }
      }
    }
  }

  template <CorrelationContainer::CFStep step, typename TFT0s>
  void fillCorrelationsFT0AFT0C(TFT0s const& ft0Col1, TFT0s const& ft0Col2, float posZ, int system, double multiplicity, float eventWeight) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {
    auto fSampleIndex = static_cast<int>(fRandom.Uniform(0.0, cfgSampleSize));

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
        registry.fill(HIST("Trig_hist_FT0A_FT0C"), fSampleIndex, posZ, 0.5, multiplicity, eventWeight * amplA);
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
          if (cfgQaCheck) {
            registry.fill(HIST("deltaEta_deltaPhi_same_FT0A_FT0C"), deltaPhi, deltaEta, amplA * amplC * eventWeight * triggerWeight);
          }
          sameFt0aFt0c->getPairHist()->Fill(step, fSampleIndex, posZ, 0.5, multiplicity, deltaPhi, deltaEta, amplA * amplC * eventWeight * triggerWeight);
        } else if (system == MixedEvent) {
          if (cfgQaCheck) {
            registry.fill(HIST("deltaEta_deltaPhi_mixed_FT0A_FT0C"), deltaPhi, deltaEta, amplA * amplC * eventWeight * triggerWeight);
          }
          mixedFt0aFt0c->getPairHist()->Fill(step, fSampleIndex, posZ, 0.5, multiplicity, deltaPhi, deltaEta, amplA * amplC * eventWeight * triggerWeight);
        }
      }
    }
  }

  template <CorrelationContainer::CFStep step, typename TTracks, typename TTracksAssoc>
  void fillCorrelations(const TTracks& tracks1, const TTracksAssoc& tracks2, float posZ, int system, double multiplicity, int magneticField, float eventWeight) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {

    auto fSampleIndex = static_cast<int>(fRandom.Uniform(0.0, cfgSampleSize));

    float triggerWeight = 1.0f;

    float associateWeight = 1.0f;

    // loop over all tracks
    for (auto const& track1 : tracks1) {

      if (!trackSelected(track1)) {
        continue;
      }

      // Fill Nsigma QA
      if (cfgPIDConfigs.cfgPIDParticle > UseCharged) {
        fillNsigmaBeforeCut(track1, cfgPIDConfigs.cfgPIDParticle);
      }

      if (cfgPIDConfigs.cfgPIDParticle && getNsigmaPID(track1, false) != cfgPIDConfigs.cfgPIDParticle) {
        continue; // if PID is selected, check if the track has the right PID
      }

      if (cfgPIDConfigs.cfgPIDParticle > UseCharged) {
        fillNsigmaAfterCut(track1, cfgPIDConfigs.cfgPIDParticle);
      }

      if (!getEfficiencyCorrectionNch(triggerWeight, track1.pt())) {
        continue;
      }

      if (system == SameEvent) {
        registry.fill(HIST("Trig_hist_TPC"), fSampleIndex, posZ, track1.pt(), multiplicity, eventWeight * triggerWeight);
      }

      for (auto const& track2 : tracks2) {

        if (!trackSelected(track2)) {
          continue;
        }

        if (!getEfficiencyCorrectionNch(associateWeight, track2.pt())) {
          continue;
        }

        if (cfgUseRefpTCut) {
          if (track2.pt() > cfgRefpTMax || track2.pt() < cfgRefpTMin) {
            continue;
          }
        }
        if (track1.pt() <= track2.pt()) {
          continue; // skip if the trigger pt is less than the associate pt
        }

        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -PIHalf);
        float deltaEta = track1.eta() - track2.eta();

        if (cfgApplyTwoTrackEfficiency && std::abs(deltaEta) < cfgMergingCut) {

          double dPhiStarHigh = getDPhiStar(track1, track2, cfgRadiusHigh, magneticField);
          double dPhiStarLow = getDPhiStar(track1, track2, cfgRadiusLow, magneticField);

          const double kLimit = 3.0 * cfgMergingCut;

          bool bIsBelow = false;

          constexpr double loopIncrement = 0.01;
          const int nLoops = static_cast<int>((cfgRadiusHigh - cfgRadiusLow) / loopIncrement);

          if (std::abs(dPhiStarLow) < kLimit || std::abs(dPhiStarHigh) < kLimit || dPhiStarLow * dPhiStarHigh < 0) {
            for (int i = 0; i < nLoops; i++) {
              double rad = cfgRadiusLow + i * loopIncrement;
              double dPhiStar = getDPhiStar(track1, track2, rad, magneticField);
              if (std::abs(dPhiStar) < kLimit) {
                bIsBelow = true;
                break;
              }
            }
            if (bIsBelow) {
              continue;
            }
          }
        }

        // fill the right sparse and histograms
        if (system == SameEvent) {

          sameTPC->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), multiplicity, deltaPhi, deltaEta, triggerWeight * associateWeight);
          if (cfgQaCheck) {
            registry.fill(HIST("deltaEta_deltaPhi_same_TPC"), deltaPhi, deltaEta, triggerWeight * associateWeight);
          }
        } else if (system == MixedEvent) {

          mixedTPC->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), multiplicity, deltaPhi, deltaEta, triggerWeight * associateWeight);

          if (cfgQaCheck) {
            registry.fill(HIST("deltaEta_deltaPhi_mixed_TPC"), deltaPhi, deltaEta, triggerWeight * associateWeight);
          }
        }
      }
    }
  }

  bool isGoodRun(int runNumber)
  {
    for (const auto& excludedRun : cfgRunRemoveList.value) {
      if (runNumber == excludedRun) {
        return false;
      }
    }

    return true;
  }

  void processSameTpcFt0a(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
  {
    if (cfgQaCheck) {
      eventSelectedIndividually(collision);
    }

    if (!collision.sel8()) {
      return;
    }

    if (!eventRct(collision, true)) {
      return;
    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    int currentRunNumber = bc.runNumber();
    if (!cfgRunRemoveList.value.empty()) {
      if (!isGoodRun(currentRunNumber)) { // Rejects runs if bad run number
        return;
      }
    }

    if (cfgEvSel.cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), true)) {
      return;
    }

    if (!collision.has_foundFT0()) {
      return;
    }

    loadAlignParam(bc.timestamp());
    loadGain(bc);
    loadCorrection(bc.timestamp());
    float eventWeight = 1.0f;

    registry.fill(HIST("zVtx"), collision.posZ());

    registry.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin

    fillYield(collision, tracks);

    auto multiplicity = static_cast<double>(tracks.size());

    if (cfgQaCheck) {
      registry.fill(HIST("Nch"), multiplicity);
    }

    if (cfgStrictTrackCounter) {
      trackCounter(tracks, multiplicity);
    }

    if (cfgQaCheck) {
      registry.fill(HIST("Nch_corrected"), multiplicity);
    }

    if (multiplicity > cfgEvSel.cfgMaxMultForCorrelations || multiplicity < cfgEvSel.cfgMinMultForCorrelations) {
      return;
    }

    const auto& ft0 = collision.foundFT0();
    fillCorrelationsTPCFT0<CorrelationContainer::kCFStepReconstructed>(tracks, ft0, collision.posZ(), SameEvent, multiplicity, IndexFT0A, eventWeight);
  }
  PROCESS_SWITCH(CorrFit, processSameTpcFt0a, "Process same event for TPC-FT0 correlation", false);

  void processMixedTpcFt0a(FilteredCollisions const& collisions, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
  {

    auto getTracksSize = [&tracks, this](FilteredCollisions::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}};

    auto tracksTuple = std::make_tuple(tracks, tracks);
    Pair<FilteredCollisions, FilteredTracks, FilteredTracks, MixedBinning> pairs{binningOnVtxAndMult, cfgMinMixEventNum, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;

      if (!collision1.sel8() || !collision2.sel8()) {
        continue;
      }

      if (!eventRct(collision1, false) || !eventRct(collision2, false)) {
        continue;
      }

      if (cfgEvSel.cfgUseAdditionalEventCut && !eventSelected(collision1, tracks1.size(), false)) {
        continue;
      }
      if (cfgEvSel.cfgUseAdditionalEventCut && !eventSelected(collision2, tracks2.size(), false)) {
        continue;
      }

      if (!(collision1.has_foundFT0() && collision2.has_foundFT0())) {
        continue;
      }

      registry.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      int currentRunNumber = bc.runNumber();
      if (!cfgRunRemoveList.value.empty()) {
        if (!isGoodRun(currentRunNumber)) { // Rejects runs if bad run number
          continue;
        }
      }

      loadAlignParam(bc.timestamp());
      loadCorrection(bc.timestamp());
      float eventWeight = 1.0f;

      auto multiplicity = static_cast<double>(tracks1.size());

      if (cfgStrictTrackCounter) {
        trackCounter(tracks1, multiplicity);
      }

      if (multiplicity > cfgEvSel.cfgMaxMultForCorrelations || multiplicity < cfgEvSel.cfgMinMultForCorrelations) {
        continue;
      }

      const auto& ft0 = collision2.foundFT0();
      fillCorrelationsTPCFT0<CorrelationContainer::kCFStepReconstructed>(tracks1, ft0, collision1.posZ(), MixedEvent, multiplicity, IndexFT0A, eventWeight);
    }
  }
  PROCESS_SWITCH(CorrFit, processMixedTpcFt0a, "Process mixed events for TPC-FT0A correlation", false);

  void processSameTpcFt0c(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
  {

    if (cfgQaCheck) {
      eventSelectedIndividually(collision);
    }

    if (!collision.sel8()) {
      return;
    }

    if (!eventRct(collision, true)) {
      return;
    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int currentRunNumber = bc.runNumber();
    if (!cfgRunRemoveList.value.empty()) {
      if (!isGoodRun(currentRunNumber)) { // Rejects runs if bad run number
        return;
      }
    }

    if (cfgEvSel.cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), true)) {
      return;
    }
    if (!collision.has_foundFT0()) {
      return;
    }
    loadAlignParam(bc.timestamp());
    loadGain(bc);
    loadCorrection(bc.timestamp());

    registry.fill(HIST("zVtx"), collision.posZ());

    registry.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin

    const auto& ft0 = collision.foundFT0();

    auto multiplicity = static_cast<double>(tracks.size());

    if (cfgQaCheck) {
      registry.fill(HIST("Nch"), multiplicity);
    }

    if (cfgStrictTrackCounter) {
      trackCounter(tracks, multiplicity);
    }

    if (cfgQaCheck) {
      registry.fill(HIST("Nch_corrected"), multiplicity);
    }

    if (multiplicity > cfgEvSel.cfgMaxMultForCorrelations || multiplicity < cfgEvSel.cfgMinMultForCorrelations) {
      return;
    }

    fillCorrelationsTPCFT0<CorrelationContainer::kCFStepReconstructed>(tracks, ft0, collision.posZ(), SameEvent, multiplicity, IndexFT0C, 1.0f);
  }
  PROCESS_SWITCH(CorrFit, processSameTpcFt0c, "Process same event for TPC-FT0C correlation", false);

  void processMixedTpcFt0c(FilteredCollisions const& collisions, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
  {

    auto getTracksSize = [&tracks, this](FilteredCollisions::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}};

    auto tracksTuple = std::make_tuple(tracks, tracks);
    Pair<FilteredCollisions, FilteredTracks, FilteredTracks, MixedBinning> pairs{binningOnVtxAndMult, cfgMinMixEventNum, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;
      if (!collision1.sel8() || !collision2.sel8()) {
        continue;
      }

      if (!eventRct(collision1, false) || !eventRct(collision2, false)) {
        continue;
      }

      if (cfgEvSel.cfgUseAdditionalEventCut && !eventSelected(collision1, tracks1.size(), false)) {
        continue;
      }
      if (cfgEvSel.cfgUseAdditionalEventCut && !eventSelected(collision2, tracks2.size(), false)) {
        continue;
      }

      if (!(collision1.has_foundFT0() && collision2.has_foundFT0())) {
        continue;
      }

      registry.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin
      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      int currentRunNumber = bc.runNumber();
      if (!cfgRunRemoveList.value.empty()) {
        if (!isGoodRun(currentRunNumber)) { // Rejects runs if bad run number
          continue;
        }
      }
      loadAlignParam(bc.timestamp());
      loadCorrection(bc.timestamp());
      float eventWeight = 1.0f;

      const auto& ft0 = collision2.foundFT0();
      auto multiplicity = static_cast<double>(tracks1.size());

      if (cfgStrictTrackCounter) {
        trackCounter(tracks1, multiplicity);
      }

      if (multiplicity > cfgEvSel.cfgMaxMultForCorrelations || multiplicity < cfgEvSel.cfgMinMultForCorrelations) {
        continue;
      }

      fillCorrelationsTPCFT0<CorrelationContainer::kCFStepReconstructed>(tracks1, ft0, collision1.posZ(), MixedEvent, multiplicity, IndexFT0C, eventWeight);
    }
  }
  PROCESS_SWITCH(CorrFit, processMixedTpcFt0c, "Process mixed events for TPC-FT0C correlation", false);

  void processSameFt0aFt0c(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
  {

    if (cfgQaCheck) {
      eventSelectedIndividually(collision);
    }

    if (!collision.sel8()) {
      return;
    }

    if (!eventRct(collision, true)) {
      return;
    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int currentRunNumber = bc.runNumber();
    if (!cfgRunRemoveList.value.empty()) {
      if (!isGoodRun(currentRunNumber)) { // Rejects runs if bad run number
        return;
      }
    }

    if (cfgEvSel.cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), true)) {
      return;
    }

    if (!collision.has_foundFT0()) {
      return;
    }

    loadAlignParam(bc.timestamp());
    loadGain(bc);
    loadCorrection(bc.timestamp());
    float eventWeight = 1.0f;

    const auto& ft0 = collision.foundFT0();

    fillYield(collision, tracks);

    registry.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin

    auto multiplicity = static_cast<double>(tracks.size());

    if (cfgQaCheck) {
      registry.fill(HIST("Nch"), multiplicity);
    }

    if (cfgStrictTrackCounter) {
      trackCounter(tracks, multiplicity);
    }

    if (cfgQaCheck) {
      registry.fill(HIST("Nch_corrected"), multiplicity);
    }

    if (multiplicity > cfgEvSel.cfgMaxMultForCorrelations || multiplicity < cfgEvSel.cfgMinMultForCorrelations) {
      return;
    }

    fillCorrelationsFT0AFT0C<CorrelationContainer::kCFStepReconstructed>(ft0, ft0, collision.posZ(), SameEvent, multiplicity, eventWeight);
  }
  PROCESS_SWITCH(CorrFit, processSameFt0aFt0c, "Process same event for FT0A-FT0C correlation", true);

  void processMixedFt0aFt0c(FilteredCollisions const& collisions, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
  {

    auto getTracksSize = [&tracks, this](FilteredCollisions::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}};

    auto tracksTuple = std::make_tuple(tracks, tracks);
    Pair<FilteredCollisions, FilteredTracks, FilteredTracks, MixedBinning> pairs{binningOnVtxAndMult, cfgMinMixEventNum, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;

      // should have the same event to TPC-FT0A/C correlations
      if (!collision1.sel8() || !collision2.sel8()) {
        continue;
      }

      if (!eventRct(collision1, false) || !eventRct(collision2, false)) {
        continue;
      }

      if (cfgEvSel.cfgUseAdditionalEventCut && !eventSelected(collision1, tracks1.size(), false)) {
        continue;
      }
      if (cfgEvSel.cfgUseAdditionalEventCut && !eventSelected(collision2, tracks2.size(), false)) {
        continue;
      }

      if (!(collision1.has_foundFT0() && collision2.has_foundFT0())) {
        continue;
      }

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      int currentRunNumber = bc.runNumber();
      if (!cfgRunRemoveList.value.empty()) {
        if (!isGoodRun(currentRunNumber)) { // Rejects runs if bad run number
          continue;
        }
      }
      loadAlignParam(bc.timestamp());
      loadCorrection(bc.timestamp());
      float eventWeight = 1.0f;

      const auto& ft0Col1 = collision1.foundFT0();
      const auto& ft0Col2 = collision2.foundFT0();

      auto multiplicity = static_cast<double>(tracks1.size());

      if (cfgStrictTrackCounter) {
        trackCounter(tracks1, multiplicity);
      }

      if (multiplicity > cfgEvSel.cfgMaxMultForCorrelations || multiplicity < cfgEvSel.cfgMinMultForCorrelations) {
        continue;
      }

      registry.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin

      fillCorrelationsFT0AFT0C<CorrelationContainer::kCFStepReconstructed>(ft0Col1, ft0Col2, collision1.posZ(), MixedEvent, multiplicity, eventWeight);
    }
  }
  PROCESS_SWITCH(CorrFit, processMixedFt0aFt0c, "Process mixed events for FT0A-FT0C correlation", true);

  void processSameTPC(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::BCsWithTimestamps const&)
  {

    if (cfgQaCheck) {
      eventSelectedIndividually(collision);
    }
    if (!eventRct(collision, true)) {
      return;
    }

    if (!collision.sel8()) {
      return;
    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int currentRunNumber = bc.runNumber();
    if (!cfgRunRemoveList.value.empty()) {
      if (!isGoodRun(currentRunNumber)) { // Rejects runs if bad run number
        return;
      }
    }

    if (cfgEvSel.cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), true)) {
      return;
    }

    registry.fill(HIST("zVtx"), collision.posZ());

    registry.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin
    loadCorrection(bc.timestamp());

    fillYield(collision, tracks);

    auto multiplicity = static_cast<double>(tracks.size());

    if (cfgQaCheck) {
      registry.fill(HIST("Nch"), multiplicity);
    }

    if (cfgStrictTrackCounter) {
      trackCounter(tracks, multiplicity);
    }

    if (cfgQaCheck) {
      registry.fill(HIST("Nch_corrected"), multiplicity);
    }

    if (multiplicity > cfgEvSel.cfgMaxMultForCorrelations || multiplicity < cfgEvSel.cfgMinMultForCorrelations) {
      return;
    }

    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(tracks, tracks, collision.posZ(), SameEvent, multiplicity, getMagneticField(bc.timestamp()), 1.0f);
  }
  PROCESS_SWITCH(CorrFit, processSameTPC, "Process same event for TPC-TPC correlation", false);

  void processMixedTPC(FilteredCollisions const& collisions, FilteredTracks const& tracks, aod::BCsWithTimestamps const&)
  {

    auto getTracksSize = [&tracks, this](FilteredCollisions::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}};

    auto tracksTuple = std::make_tuple(tracks, tracks);
    Pair<FilteredCollisions, FilteredTracks, FilteredTracks, MixedBinning> pairs{binningOnVtxAndMult, cfgMinMixEventNum, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;

      if (!eventRct(collision1, false) || !eventRct(collision2, false)) {
        continue;
      }

      if (!collision1.sel8() || !collision2.sel8()) {
        continue;
      }

      if (cfgEvSel.cfgUseAdditionalEventCut && !eventSelected(collision1, tracks1.size(), false)) {
        continue;
      }
      if (cfgEvSel.cfgUseAdditionalEventCut && !eventSelected(collision2, tracks2.size(), false)) {
        continue;
      }

      registry.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      int currentRunNumber = bc.runNumber();
      if (!cfgRunRemoveList.value.empty()) {
        if (!isGoodRun(currentRunNumber)) { // Rejects runs if bad run number
          continue;
        }
      }

      loadCorrection(bc.timestamp());

      auto multiplicity = static_cast<double>(tracks1.size());

      if (cfgStrictTrackCounter) {
        trackCounter(tracks1, multiplicity);
      }

      if (multiplicity > cfgEvSel.cfgMaxMultForCorrelations || multiplicity < cfgEvSel.cfgMinMultForCorrelations) {
        continue;
      }

      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(tracks1, tracks2, collision1.posZ(), MixedEvent, multiplicity, getMagneticField(collision1.bc_as<aod::BCsWithTimestamps>().timestamp()), 1.0f);
    }
  }
  PROCESS_SWITCH(CorrFit, processMixedTPC, "Process mixed events for TPC-TPC correlation", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CorrFit>(cfgc),
  };
}
