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

  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, true, "Use additional event cut on mult correlations")
  O2_DEFINE_CONFIGURABLE(cfgZVtxCut, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgQABasic, bool, true, "Enable QA histograms for event selection")
  O2_DEFINE_CONFIGURABLE(cfgStrictTrackCounter, bool, false, "Strict track counter for multiplicity correlation cut, counts only tracks that pass all cuts and are used in the correlation")
  O2_DEFINE_CONFIGURABLE(cfgRefpTt, bool, false, "Apply upper pT cut on reference tracks")
  O2_DEFINE_CONFIGURABLE(cfgRefpTMax, float, 3.0f, "maximum pT for reference tracks if cfgRefpTt is true")
  O2_DEFINE_CONFIGURABLE(cfgMinMultForCorrelations, int, 0, "minimum multiplicity for correlations")
  O2_DEFINE_CONFIGURABLE(cfgMaxMultForCorrelations, int, 20, "maximum multiplicity for correlations")
  O2_DEFINE_CONFIGURABLE(cfgSelCollByNch, bool, true, "Select collisions by Nch or centrality")
  O2_DEFINE_CONFIGURABLE(cfgMinCentForCorrelations, int, 0, "minimum centrality for correlations")
  O2_DEFINE_CONFIGURABLE(cfgMaxCentForCorrelations, int, 20, "maximum centrality for correlations")
  O2_DEFINE_CONFIGURABLE(cfgRefMultiplicity, bool, false, "Use multiplicity of reference tracks for multiplicity correlation cut instead of Nch")
  Configurable<std::vector<int>> cfgRunRemoveList{"cfgRunRemoveList", std::vector<int>{-1}, "excluded run numbers"};
  O2_DEFINE_CONFIGURABLE(cfgEvSelRCTflags, std::string, "", "keep empty to disable, usage: 'CentralBarrelTracking (CBT)', 'CBT_hadronPID' ")
  O2_DEFINE_CONFIGURABLE(cfgAnalyzeTPCFT0A, bool, true, "Switch for doing TPC-FT0A correlations")
  O2_DEFINE_CONFIGURABLE(cfgAnalyzeTPCFT0C, bool, true, "Switch for doing TPC-FT0C correlations")
  O2_DEFINE_CONFIGURABLE(cfgUseEventWeights, bool, false, "Use event weights for mixed event")

  struct : ConfigurableGroup{
             O2_DEFINE_CONFIGURABLE(cfgPtCutMin, float, 0.2f, "minimum accepted track pT")
               O2_DEFINE_CONFIGURABLE(cfgPtCutMax, float, 10.0f, "maximum accepted track pT")
                 O2_DEFINE_CONFIGURABLE(cfgEtaCut, float, 0.8f, "Eta cut")
                   O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5f, "max chi2 per TPC clusters")
                     O2_DEFINE_CONFIGURABLE(cfgCutTPCclu, float, 50.0f, "minimum TPC clusters")
                       O2_DEFINE_CONFIGURABLE(cfgCutTPCCrossedRows, float, 70.0f, "minimum TPC crossed rows")
                         O2_DEFINE_CONFIGURABLE(cfgCutITSclu, float, 5.0f, "minimum ITS clusters")
                           O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "max DCA to vertex z")} cfgTrackCuts;

  struct : ConfigurableGroup{
             O2_DEFINE_CONFIGURABLE(cfgCutOccupancyHigh, int, 2000, "High cut on TPC occupancy")
               O2_DEFINE_CONFIGURABLE(cfgCutOccupancyLow, int, 0, "Low cut on TPC occupancy")} cfgEventSelection;

  Configurable<LabeledArray<int>> cfgUseEventCuts{"cfgUseEventCuts", {LongArrayInt.front().data(), 14, 1, {"Filtered Events", "Sel8", "kNoTimeFrameBorder", "kNoITSROFrameBorder", "kNoSameBunchPileup", "kIsGoodZvtxFT0vsPV", "kNoCollInTimeRangeStandard", "kIsGoodITSLayersAll", "kIsGoodITSLayer0123", "kNoCollInRofStandard", "kNoHighMultCollInPrevRof", "Occupancy", "Multcorrelation", "T0AV0ACut"}, {"EvCuts"}}, "Labeled array (int) for various cuts on resonances"};

  O2_DEFINE_CONFIGURABLE(cfgMinMixEventNum, int, 5, "Minimum number of events to mix")

  O2_DEFINE_CONFIGURABLE(cfgSampleSize, double, 10, "Sample size for mixed event")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgEfficiencyNch, std::string, "", "CCDB path to multiplicity dependent efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgCentralityWeight, std::string, "", "CCDB path to centrality weight object")
  O2_DEFINE_CONFIGURABLE(cfgLocalEfficiency, bool, false, "Use local efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgLocalEfficiencyNch, bool, false, "Use local multiplicity dependent efficiency object");
  O2_DEFINE_CONFIGURABLE(cfgCentEstimator, int, 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FT0A")

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
    TF1* fMultPVT0CCutLow = nullptr;
    TF1* fMultPVT0CCutHigh = nullptr;
    TF1* fMultT0CCutLow = nullptr;
    TF1* fMultT0CCutHigh = nullptr;
    TF1* fMultGlobalPVCutLow = nullptr;
    TF1* fMultGlobalPVCutHigh = nullptr;
    TF1* fMultMultV0ACutLow = nullptr;
    TF1* fMultMultV0ACutHigh = nullptr;
    TF1* fT0AV0AMean = nullptr;
    TF1* fT0AV0ASigma = nullptr;
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

  Configurable<float> cfgCutFV0{"cfgCutFV0", 50., "FV0A threshold"};
  Configurable<float> cfgCutFT0A{"cfgCutFT0A", 150., "FT0A threshold"};
  Configurable<float> cfgCutFT0C{"cfgCutFT0C", 50., "FT0C threshold"};
  Configurable<float> cfgCutZDC{"cfgCutZDC", 10., "ZDC threshold"};

  SliceCache cache;

  // Axes for general QA
  ConfigurableAxis axisVertex{"axisVertex", {10, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisMult{"axisMult", {10, 0, 100}, "multiplicity axis for QA histograms"};
  ConfigurableAxis axisCent{"axisCent", {100, 0, 100}, "centrality axis for QA histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {72, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisPtFiner{"axisPtFiner", {98, 0.2, 10.0}, "pt axis for histograms"};

  // Axes for Correlations
  ConfigurableAxis axisMultForCorr{"axisMultForCorr", {VARIABLE_WIDTH, 0, 5, 10, 20, 50, 100, 150, 200, 300, 500}, "multiplicity axis for correlation histograms"};
  ConfigurableAxis axisCentForCorr{"axisCentForCorr", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, "centrality axis for correlation histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEtaTpcFt0a{"axisDeltaEtaTpcFt0a", {32, -5.8, -2.6}, "delta eta axis, -5.8~-2.6 for TPC-FT0A,"};
  ConfigurableAxis axisDeltaEtaTpcFt0c{"axisDeltaEtaTpcFt0c", {32, 1.2, 4.2}, "delta eta axis, 1.2~4.2 for TPC-FT0C"};
  ConfigurableAxis axisDeltaEtaFt0aFt0c{"axisDeltaEtaFt0aFt0c", {32, -1.5, 3.0}, "delta eta axis"};
  ConfigurableAxis axisVtxMix{"axisVtxMix", {VARIABLE_WIDTH, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, "vertex axis for mixed event histograms"};
  ConfigurableAxis axisMultMix{"axisMultMix", {VARIABLE_WIDTH, 0, 10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260}, "multiplicity / centrality axis for mixed event histograms"};
  ConfigurableAxis axisSample{"axisSample", {cfgSampleSize, 0, cfgSampleSize}, "sample axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt trigger axis for histograms"};
  ConfigurableAxis axisParticle{"axisParticle", {4, 0, 4}, "particle axis for correlation container"};
  ConfigurableAxis axisTpcNcls{"axisTpcNcls", {500, 0, 500}, "number of TPC clusters used for PID"};

  // Axes for long-range correlations QA
  ConfigurableAxis axisAmplitudeFt0{"axisAmplitudeFt0a", {5000, 0, 1000}, "FT0A amplitude"};
  ConfigurableAxis axisChannelFt0aAxis{"axisChannelFt0aAxis", {96, 0.0, 96.0}, "FT0A channel"};
  ConfigurableAxis axisChannelFt0cAxis{"axisChannelFt0cAxis", {96, 0.0, 96.0}, "FT0C channel"};
  Configurable<std::string> cfgGainEqPath{"cfgGainEqPath", "Analysis/EventPlane/GainEq", "CCDB path for gain equalization constants"};

  Configurable<int> cfgCorrLevel{"cfgCorrLevel", 0, "calibration step: 0 = no corr, 1 = gain corr"};
  ConfigurableAxis cfgaxisFITamp{"cfgaxisFITamp", {1000, 0, 5000}, ""};
  AxisSpec axisFit{cfgaxisFITamp, "fit amplitude"};
  AxisSpec axisChID = {220, 0, 220};

  // Axes for PID QA
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {80, -5, 5}, "nsigmaTPC axis"};
  ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {80, -5, 5}, "nsigmaTOF axis"};
  ConfigurableAxis axisNsigmaITS{"axisNsigmaITS", {80, -5, 5}, "nsigmaITS axis"};
  ConfigurableAxis axisTpcSignal{"axisTpcSignal", {250, 0, 250}, "dEdx axis for TPC"};
  ConfigurableAxis axisSigma{"axisSigma", {200, 0, 20}, "sigma axis for TPC"};

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgZVtxCut);
  Filter trackFilter = (nabs(aod::track::eta) < cfgTrackCuts.cfgEtaCut) && (cfgTrackCuts.cfgPtCutMin < aod::track::pt) && (cfgTrackCuts.cfgPtCutMax > aod::track::pt) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t)true));

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
  std::array<std::array<int, 1>, 14> eventCuts{};
  std::array<std::array<float, 3>, 6> nSigmaVals{};

  // define global variables
  TRandom3* gRandom = new TRandom3();

  enum EventType {
    SameEvent = 1,
    MixedEvent = 3
  };

  enum FITIndex {
    kFT0A = 0,
    kFT0C = 1
  };

  enum PIDIndex {
    kCharged = 0,
    kPions,
    kKaons,
    kProtons,
    kK0,
    kLambda,
    kPhi
  };
  enum PiKpArrayIndex {
    iPionUp = 0,
    iKaonUp,
    iProtonUp,
    iPionLow,
    iKaonLow,
    iProtonLow
  };
  enum DetectorType {
    kTPC = 0,
    kTOF,
    kITS
  };
  enum Stage {
    Before = 0,
    After
  };

  enum EventCutTypes {
    kFilteredEvents = 0,
    kAfterSel8,
    kUseNoTimeFrameBorder,
    kUseNoITSROFrameBorder,
    kUseNoSameBunchPileup,
    kUseGoodZvtxFT0vsPV,
    kUseNoCollInTimeRangeStandard,
    kUseGoodITSLayersAll,
    kUseGoodITSLayer0123,
    kUseNoCollInRofStandard,
    kUseNoHighMultCollInPrevRof,
    kUseOccupancy,
    kUseMultCorrCut,
    kUseT0AV0ACut,
    kHaveFT0Cut,
    kNEventCuts
  };

  enum EventCutType {
    kEvCut1 = 0,
    kNEvCutTypes = 1
  };

  enum CentEstimators {
    kCentFT0C = 0,
    kCentFT0CVariant1,
    kCentFT0M,
    kCentFV0A,
    // Count the total number of enum
    kCount_CentEstimators
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

    readMatrix(cfgUseEventCuts->getData(), eventCuts);
    readMatrix(cfgPIDConfigs.nSigmas->getData(), nSigmaVals);

    if (doprocessSameFt0aFt0c || doprocessSameTpcFt0) {
      histos.add("hEventCountRct", "Number of Event;; Count", {HistType::kTH1D, {{2, 0, 2}}});
      histos.get<TH1>(HIST("hEventCountRct"))->GetXaxis()->SetBinLabel(1, "rct fail");
      histos.get<TH1>(HIST("hEventCountRct"))->GetXaxis()->SetBinLabel(2, "rct pass");
      histos.add("hEventCount", "Number of Event;; Count", {HistType::kTH1D, {{15, -0.5, 14.5}}});
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kFilteredEvents + 1, "Filtered events");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kAfterSel8 + 1, "After sel8");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoTimeFrameBorder + 1, "kNoTimeFrameBorder");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoITSROFrameBorder + 1, "kNoITSROFrameBorder");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoSameBunchPileup + 1, "kNoSameBunchPileup");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseGoodZvtxFT0vsPV + 1, "kIsGoodZvtxFT0vsPV");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoCollInTimeRangeStandard + 1, "kNoCollInTimeRangeStandard");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseGoodITSLayersAll + 1, "kIsGoodITSLayersAll");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseGoodITSLayer0123 + 1, "kIsGoodITSLayer0123");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoCollInRofStandard + 1, "kNoCollInRofStandard");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoHighMultCollInPrevRof + 1, "kNoHighMultCollInPrevRof");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseOccupancy + 1, "Occupancy Cut");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseMultCorrCut + 1, "MultCorrelation Cut");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseT0AV0ACut + 1, "T0AV0A cut");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kHaveFT0Cut + 1, "Event has FT0 cut");
    }

    if ((doprocessSameFt0aFt0c || doprocessSameTpcFt0) && cfgQABasic) {
      histos.add("hPassedEventSelection", "Number of Event;; Count", {HistType::kTH1D, {{12, -0.5, 11.5}}});
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kFilteredEvents + 1, "Filtered events");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kAfterSel8 + 1, "After sel8");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kUseNoTimeFrameBorder + 1, "kNoTimeFrameBorder");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kUseNoITSROFrameBorder + 1, "kNoITSROFrameBorder");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kUseNoSameBunchPileup + 1, "kNoSameBunchPileup");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kUseGoodZvtxFT0vsPV + 1, "kIsGoodZvtxFT0vsPV");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kUseNoCollInTimeRangeStandard + 1, "kNoCollInTimeRangeStandard");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kUseGoodITSLayersAll + 1, "kIsGoodITSLayersAll");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kUseGoodITSLayer0123 + 1, "kIsGoodITSLayer0123");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kUseNoCollInRofStandard + 1, "kNoCollInRofStandard");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kUseNoHighMultCollInPrevRof + 1, "kNoHighMultCollInPrevRof");
      histos.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kUseOccupancy + 1, "Occupancy Cut");
    }

    // Multiplicity correlation cuts
    if (eventCuts[kUseMultCorrCut][kEvCut1] != 0) {
      cfgFuncParas.multT0CCutPars = cfgFuncParas.cfgMultT0CCutPars;
      cfgFuncParas.multPVT0CCutPars = cfgFuncParas.cfgMultPVT0CCutPars;
      cfgFuncParas.multGlobalPVCutPars = cfgFuncParas.cfgMultGlobalPVCutPars;
      cfgFuncParas.multMultV0ACutPars = cfgFuncParas.cfgMultMultV0ACutPars;
      cfgFuncParas.fMultPVT0CCutLow = new TF1("fMultPVT0CCutLow", cfgFuncParas.cfgMultCentLowCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultPVT0CCutLow->SetParameters(cfgFuncParas.multPVT0CCutPars.data());
      cfgFuncParas.fMultPVT0CCutHigh = new TF1("fMultPVT0CCutHigh", cfgFuncParas.cfgMultCentHighCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultPVT0CCutHigh->SetParameters(cfgFuncParas.multPVT0CCutPars.data());

      cfgFuncParas.fMultT0CCutLow = new TF1("fMultT0CCutLow", cfgFuncParas.cfgMultCentLowCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultT0CCutLow->SetParameters(cfgFuncParas.multT0CCutPars.data());
      cfgFuncParas.fMultT0CCutHigh = new TF1("fMultT0CCutHigh", cfgFuncParas.cfgMultCentHighCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultT0CCutHigh->SetParameters(cfgFuncParas.multT0CCutPars.data());

      cfgFuncParas.fMultGlobalPVCutLow = new TF1("fMultGlobalPVCutLow", cfgFuncParas.cfgMultMultPVLowCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultGlobalPVCutLow->SetParameters(cfgFuncParas.multGlobalPVCutPars.data());
      cfgFuncParas.fMultGlobalPVCutHigh = new TF1("fMultGlobalPVCutHigh", cfgFuncParas.cfgMultMultPVHighCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultGlobalPVCutHigh->SetParameters(cfgFuncParas.multGlobalPVCutPars.data());

      cfgFuncParas.fMultMultV0ACutLow = new TF1("fMultMultV0ACutLow", cfgFuncParas.cfgMultMultV0ALowCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultMultV0ACutLow->SetParameters(cfgFuncParas.multMultV0ACutPars.data());
      cfgFuncParas.fMultMultV0ACutHigh = new TF1("fMultMultV0ACutHigh", cfgFuncParas.cfgMultMultV0AHighCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultMultV0ACutHigh->SetParameters(cfgFuncParas.multMultV0ACutPars.data());
    }
    if (eventCuts[kUseT0AV0ACut][kEvCut1] != 0) {
      cfgFuncParas.fT0AV0AMean = new TF1("fT0AV0AMean", "[0]+[1]*x", 0, 200000);
      cfgFuncParas.fT0AV0AMean->SetParameters(-1601.0581, 9.417652e-01);
      cfgFuncParas.fT0AV0ASigma = new TF1("fT0AV0ASigma", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 200000);
      cfgFuncParas.fT0AV0ASigma->SetParameters(463.4144, 6.796509e-02, -9.097136e-07, 7.971088e-12, -2.600581e-17);
    }

    // Choose if it is Nch selection or Centrality selection
    AxisSpec axisEventClass = {axisMultForCorr, "N_{ch}"};

    if (!cfgSelCollByNch) {
      axisEventClass = {axisCentForCorr, "Centrality"};
    }

    if (doprocessSameFt0aFt0c || doprocessSameTpcFt0) {
      histos.add("zVtx", "zVtx", {HistType::kTH1D, {axisVertex}});
      histos.add("Nch", "N_{ch}", {HistType::kTH1D, {axisMult}});
      histos.add("Nch_corrected", "N_{ch} corrected", {HistType::kTH1D, {axisMult}});
      histos.add("Centrality", "Centrality", {HistType::kTH1D, {axisCent}});
      histos.add("CentralityWeighted", "Centrality (weighted)", {HistType::kTH1D, {axisCent}});

      histos.add("FT0Amp", "", {HistType::kTH2F, {axisChID, axisFit}});
      histos.add("FT0AmpCorrect", "", {HistType::kTH2F, {axisChID, axisFit}});
    }

    if (doprocessSameTpcFt0) { // QA plots are based on TPC tracks, so they are only included in TPC-FT0A process and not in TPC-FT0C process
      histos.add("Phi", "Phi", {HistType::kTH2D, {axisPhi, axisParticle}});
      histos.add("Eta", "Eta", {HistType::kTH2D, {axisEta, axisParticle}});
      histos.add("EtaCorrected", "EtaCorrected", {HistType::kTH2D, {axisEta, axisParticle}});
      histos.add("pTFiner", "pTFiner", {HistType::kTH2D, {axisPtFiner, axisParticle}});
      histos.add("pTFinerCorrected", "pTFinerCorrected", {HistType::kTH2D, {axisPtFiner, axisParticle}});

      histos.add("hTPCNclsPID", "hTPCNclsPID", {HistType::kTH1D, {axisTpcNcls}});
      histos.add("hTPCNclsFound", "hTPCNclsFound", {HistType::kTH1D, {axisTpcNcls}});
      histos.add("hTPCCrossedRows", "hTPCCrossedRows", {HistType::kTH1D, {axisTpcNcls}});

      // PID nSigma histograms
      if (cfgPIDConfigs.cfgQANsigma) {
        if (!cfgPIDConfigs.cfgUseItsPID) {
          histos.add("TofTpcNsigma_before", "", {HistType::kTHnSparseD, {{axisNsigmaTPC, axisNsigmaTOF, axisPtTrigger}}});
          histos.add("TofTpcNsigma_after", "", {HistType::kTHnSparseD, {{axisNsigmaTPC, axisNsigmaTOF, axisPtTrigger}}});
        } // TPC-TOF PID QA hists
        if (cfgPIDConfigs.cfgUseItsPID) {
          histos.add("TofItsNsigma_before", "", {HistType::kTHnSparseD, {{axisNsigmaITS, axisNsigmaTOF, axisPtTrigger}}});
          histos.add("TofItsNsigma_after", "", {HistType::kTHnSparseD, {{axisNsigmaITS, axisNsigmaTOF, axisPtTrigger}}});
        } // ITS-TOF PID QA hists
      } // end of PID QA hists

      // PID dEdx histograms
      if (cfgPIDConfigs.cfgQAdEdx) {
        histos.add("TpcdEdx_ptwise_beforeCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF}}});
        histos.add("ExpTpcdEdx_ptwise_beforeCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF}}});
        histos.add("ExpSigma_ptwise_beforeCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisSigma, axisNsigmaTOF}}});

        histos.add("TpcdEdx_ptwise_afterCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF}}});
        histos.add("ExpTpcdEdx_ptwise_afterCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF}}});
        histos.add("ExpSigma_ptwise_afterCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisSigma, axisNsigmaTOF}}});
      }

      // TPC-FT0A correlation histograms
      if (cfgAnalyzeTPCFT0A) {
        histos.add("deltaEta_deltaPhi_same_TPC_FT0A", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0a}}); // check to see the delta eta and delta phi distribution
        histos.add("deltaEta_deltaPhi_mixed_TPC_FT0A", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0a}});
        histos.add("Assoc_amp_same_TPC_FT0A", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0}});
        histos.add("Assoc_amp_mixed_TPC_FT0A", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0}});
        histos.add("Trig_hist_TPC_FT0A", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisEventClass, axisPtTrigger, axisParticle}}});
      }
      // TPC-FT0C correlation histograms
      if (cfgAnalyzeTPCFT0C) {
        histos.add("deltaEta_deltaPhi_same_TPC_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0c}}); // check to see the delta eta and delta phi distribution
        histos.add("deltaEta_deltaPhi_mixed_TPC_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0c}});
        histos.add("Assoc_amp_same_TPC_FT0C", "", {HistType::kTH2D, {axisChannelFt0cAxis, axisAmplitudeFt0}});
        histos.add("Assoc_amp_mixed_TPC_FT0C", "", {HistType::kTH2D, {axisChannelFt0cAxis, axisAmplitudeFt0}});
        histos.add("Trig_hist_TPC_FT0C", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisEventClass, axisPtTrigger, axisParticle}}});
      }
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
    }

    histos.add("eventcount", "bin", {HistType::kTH1F, {{4, 0, 4, "bin"}}}); // histogram to see how many events are in the same and mixed event

    LOGF(info, "Initializing correlation container");

    // Initialize Nch-related histograms and containers

    std::vector<AxisSpec> corrAxisTpcFt0a = {{axisSample, "Sample"},
                                             {axisVertex, "z-vtx (cm)"},
                                             {axisEventClass},
                                             {axisDeltaPhi, "#Delta#varphi (rad)"},
                                             {axisDeltaEtaTpcFt0a, "#Delta#eta"},
                                             {axisPtTrigger, "p_{T} (GeV/c)"},
                                             {axisParticle, "Particle, 0 = charged, 1 = pion, 2 = kaon, 3 = proton"}};

    std::vector<AxisSpec> corrAxisTpcFt0c = {{axisSample, "Sample"},
                                             {axisVertex, "z-vtx (cm)"},
                                             {axisEventClass},
                                             {axisDeltaPhi, "#Delta#varphi (rad)"},
                                             {axisDeltaEtaTpcFt0c, "#Delta#eta"},
                                             {axisPtTrigger, "p_{T} (GeV/c)"},
                                             {axisParticle, "Particle, 0 = charged, 1 = pion, 2 = kaon, 3 = proton"}};

    std::vector<AxisSpec> corrAxisFt0aFt0c = {{axisSample, "Sample"},
                                              {axisVertex, "z-vtx (cm)"},
                                              {axisEventClass},
                                              {axisDeltaPhi, "#Delta#varphi (rad)"},
                                              {axisDeltaEtaFt0aFt0c, "#Delta#eta"}};

    if (doprocessSameTpcFt0) {
      if (cfgAnalyzeTPCFT0A) {
        sameTpcFt0a.setObject(new DihadronContainer("sameEvent_TPC_FT0A", "sameEvent_TPC_FT0A", corrAxisTpcFt0a));
        mixedTpcFt0a.setObject(new DihadronContainer("mixedEvent_TPC_FT0A", "mixedEvent_TPC_FT0A", corrAxisTpcFt0a));
      }
      if (cfgAnalyzeTPCFT0C) {
        sameTpcFt0c.setObject(new DihadronContainer("sameEvent_TPC_FT0C", "sameEvent_TPC_FT0C", corrAxisTpcFt0c));
        mixedTpcFt0c.setObject(new DihadronContainer("mixedEvent_TPC_FT0C", "mixedEvent_TPC_FT0C", corrAxisTpcFt0c));
      }
    }

    if (doprocessSameFt0aFt0c) {
      sameFt0aFt0c.setObject(new DihadronContainer("sameEvent_FT0A_FT0C", "sameEvent_FT0A_FT0C", corrAxisFt0aFt0c));
      mixedFt0aFt0c.setObject(new DihadronContainer("mixedEvent_FT0A_FT0C", "mixedEvent_FT0A_FT0C", corrAxisFt0aFt0c));
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
    switch (cfgCentEstimator) {
      case kCentFT0C:
        cent = collision.centFT0C();
        break;
      case kCentFT0CVariant1:
        cent = collision.centFT0CVariant1();
        break;
      case kCentFT0M:
        cent = collision.centFT0M();
        break;
      case kCentFV0A:
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
      histos.fill(HIST("hEventCount"), kFilteredEvents);
    }
    if (!collision.sel8()) {
      return false;
    }
    if (fillCounter) {
      histos.fill(HIST("hEventCount"), kAfterSel8);
    }

    if (eventCuts[kUseNoTimeFrameBorder][kEvCut1] && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (fillCounter && eventCuts[kUseNoTimeFrameBorder][kEvCut1]) {
      histos.fill(HIST("hEventCount"), kUseNoTimeFrameBorder);
    }

    if (eventCuts[kUseNoITSROFrameBorder][kEvCut1] && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (fillCounter && eventCuts[kUseNoITSROFrameBorder][kEvCut1]) {
      histos.fill(HIST("hEventCount"), kUseNoITSROFrameBorder);
    }

    if (eventCuts[kUseNoSameBunchPileup][kEvCut1] && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return false;
    }
    if (fillCounter && eventCuts[kUseNoSameBunchPileup][kEvCut1]) {
      histos.fill(HIST("hEventCount"), kUseNoSameBunchPileup);
    }

    if (eventCuts[kUseGoodZvtxFT0vsPV][kEvCut1] && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return false;
    }
    if (fillCounter && eventCuts[kUseGoodZvtxFT0vsPV][kEvCut1]) {
      histos.fill(HIST("hEventCount"), kUseGoodZvtxFT0vsPV);
    }

    if (eventCuts[kUseNoCollInTimeRangeStandard][kEvCut1] && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      // no collisions in specified time range
      return false;
    }

    if (fillCounter && eventCuts[kUseNoCollInTimeRangeStandard][kEvCut1]) {
      histos.fill(HIST("hEventCount"), kUseNoCollInTimeRangeStandard);
    }

    if (eventCuts[kUseGoodITSLayersAll][kEvCut1] && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      // from Jan 9 2025 AOT meeting
      // cut time intervals with dead ITS staves
      return false;
    }

    if (fillCounter && eventCuts[kUseGoodITSLayersAll][kEvCut1]) {
      histos.fill(HIST("hEventCount"), kUseGoodITSLayersAll);
    }

    if (eventCuts[kUseGoodITSLayer0123][kEvCut1] && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
      return false;
    }
    if (fillCounter && eventCuts[kUseGoodITSLayer0123][kEvCut1]) {
      histos.fill(HIST("hEventCount"), kUseGoodITSLayer0123);
    }

    if (eventCuts[kUseNoCollInRofStandard][kEvCut1] && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      // no other collisions in this Readout Frame with per-collision multiplicity above threshold
      return false;
    }

    if (fillCounter && eventCuts[kUseNoCollInRofStandard][kEvCut1]) {
      histos.fill(HIST("hEventCount"), kUseNoCollInRofStandard);
    }

    if (eventCuts[kUseNoHighMultCollInPrevRof][kEvCut1] && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      // veto an event if FT0C amplitude in previous ITS ROF is above threshold
      return false;
    }
    if (fillCounter && eventCuts[kUseNoHighMultCollInPrevRof][kEvCut1]) {
      histos.fill(HIST("hEventCount"), kUseNoHighMultCollInPrevRof);
    }

    auto multNTracksPV = collision.multNTracksPV();
    auto occupancy = collision.trackOccupancyInTimeRange();

    if (eventCuts[kUseOccupancy][kEvCut1] && (occupancy < cfgEventSelection.cfgCutOccupancyLow || occupancy > cfgEventSelection.cfgCutOccupancyHigh)) {
      return false;
    }
    if (fillCounter && eventCuts[kUseOccupancy][kEvCut1]) {
      histos.fill(HIST("hEventCount"), kUseOccupancy);
    }

    if (eventCuts[kUseMultCorrCut][kEvCut1]) {
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

    if (fillCounter && eventCuts[kUseMultCorrCut][kEvCut1]) {
      histos.fill(HIST("hEventCount"), kUseMultCorrCut);
    }

    // V0A T0A 5 sigma cut
    if (eventCuts[kUseT0AV0ACut][kEvCut1] && (std::fabs(collision.multFV0A() - cfgFuncParas.fT0AV0AMean->Eval(collision.multFT0A())) > cfgFuncParas.cfgV0AT0Acut * cfgFuncParas.fT0AV0ASigma->Eval(collision.multFT0A()))) {
      return false;
    }
    if (fillCounter && eventCuts[kUseT0AV0ACut][kEvCut1]) {
      histos.fill(HIST("hEventCount"), kUseT0AV0ACut);
    }

    return true;
  }

  template <typename TCollision>
  void eventSelectedIndividually(TCollision const& collision)
  {
    histos.fill(HIST("hPassedEventSelection"), kFilteredEvents);

    if (collision.sel8()) {
      histos.fill(HIST("hPassedEventSelection"), kAfterSel8);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      histos.fill(HIST("hPassedEventSelection"), kUseNoTimeFrameBorder);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      histos.fill(HIST("hPassedEventSelection"), kUseNoITSROFrameBorder);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      histos.fill(HIST("hPassedEventSelection"), kUseNoSameBunchPileup);
    }

    if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      histos.fill(HIST("hPassedEventSelection"), kUseGoodZvtxFT0vsPV);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      histos.fill(HIST("hPassedEventSelection"), kUseNoCollInTimeRangeStandard);
    }

    if (collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      histos.fill(HIST("hPassedEventSelection"), kUseGoodITSLayersAll);
    }

    if (collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
      histos.fill(HIST("hPassedEventSelection"), kUseGoodITSLayer0123);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      histos.fill(HIST("hPassedEventSelection"), kUseNoCollInRofStandard);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      histos.fill(HIST("hPassedEventSelection"), kUseNoHighMultCollInPrevRof);
    }

    auto occupancy = collision.trackOccupancyInTimeRange();
    if (eventCuts[kUseOccupancy][kEvCut1] && (occupancy < cfgEventSelection.cfgCutOccupancyLow || occupancy > cfgEventSelection.cfgCutOccupancyHigh)) {
      histos.fill(HIST("hPassedEventSelection"), kUseOccupancy);
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
    return ((track.tpcNClsFound() >= cfgTrackCuts.cfgCutTPCclu) && (track.tpcNClsCrossedRows() >= cfgTrackCuts.cfgCutTPCCrossedRows) && (track.itsNCls() >= cfgTrackCuts.cfgCutITSclu) && (track.tpcChi2NCl() < cfgTrackCuts.cfgCutChi2prTPCcls) && (track.dcaZ() < cfgTrackCuts.cfgCutDCAz));
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
    if (fitType == kFT0C) {
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
    } else if (fitType == kFT0A) {
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
    int kIndexDetector = cfgPIDConfigs.cfgUseItsPID ? kITS : kTPC;                         // Choose which nSigma to use: TPC or ITS

    bool isPion = false;
    bool isKaon = false;
    bool isProton = false;
    bool isDetectedPion = nSigmaToUse[iPionUp] < nSigmaVals[iPionUp][kIndexDetector] && nSigmaToUse[iPionUp] > nSigmaVals[iPionLow][kIndexDetector];
    bool isDetectedKaon = nSigmaToUse[iKaonUp] < nSigmaVals[iKaonUp][kIndexDetector] && nSigmaToUse[iKaonUp] > nSigmaVals[iKaonLow][kIndexDetector];
    bool isDetectedProton = nSigmaToUse[iProtonUp] < nSigmaVals[iProtonUp][kIndexDetector] && nSigmaToUse[iProtonUp] > nSigmaVals[iProtonLow][kIndexDetector];

    bool isTofPion = nSigmaTOF[iPionUp] < nSigmaVals[iPionUp][kTOF] && nSigmaTOF[iPionUp] > nSigmaVals[iPionLow][kTOF];
    bool isTofKaon = nSigmaTOF[iKaonUp] < nSigmaVals[iKaonUp][kTOF] && nSigmaTOF[iKaonUp] > nSigmaVals[iKaonLow][kTOF];
    bool isTofProton = nSigmaTOF[iProtonUp] < nSigmaVals[iProtonUp][kTOF] && nSigmaTOF[iProtonUp] > nSigmaVals[iProtonLow][kTOF];

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
      pid = kPions;
    } else if (isKaon) {
      pid = kKaons;
    } else if (isProton) {
      pid = kProtons;
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
  void trackCounter(const TTracks& tracks, double& multiplicity) // function to count the number of tracks in the event and fill the histogram
  {
    double nTracksCorrected = 0;
    float weightNch = 1.0f;
    for (auto const& track : tracks) {

      if (cfgRefMultiplicity) {
        if (track.pt() > cfgRefpTMax) {
          continue;
        }
      }

      if (!getEfficiencyCorrectionNch(weightNch, track.pt())) {
        continue;
      }

      nTracksCorrected += weightNch;
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
  void fillTrackQA(TTrack const& track, Int_t pid, Stage stage, float weff) // function to fill the QA after Nsigma selection
  {
    if (cfgQABasic) {
      if (stage == Before && pid == kCharged) {
        histos.fill(HIST("Phi"), RecoDecay::constrainAngle(track.phi(), 0.0), kCharged);
        histos.fill(HIST("Eta"), track.eta(), kCharged);
        histos.fill(HIST("EtaCorrected"), track.eta(), kCharged, weff);
        histos.fill(HIST("pTFiner"), track.pt(), kCharged);
        histos.fill(HIST("pTFinerCorrected"), track.pt(), kCharged, weff);
        histos.fill(HIST("hTPCNclsPID"), track.tpcNClsPID());
        histos.fill(HIST("hTPCNclsFound"), track.tpcNClsFound());
        histos.fill(HIST("hTPCCrossedRows"), track.tpcNClsCrossedRows());
      } else if (stage == After && pid > kCharged) {
        histos.fill(HIST("Phi"), RecoDecay::constrainAngle(track.phi(), 0.0), pid);
        histos.fill(HIST("Eta"), track.eta(), pid);
        histos.fill(HIST("EtaCorrected"), track.eta(), pid, weff);
        histos.fill(HIST("pTFiner"), track.pt(), pid);
        histos.fill(HIST("pTFinerCorrected"), track.pt(), pid, weff);
      }
    }

    if (pid > kCharged) {
      const bool useITS = cfgPIDConfigs.cfgUseItsPID;
      double tpcNSigma, tpcExpSigma, tofNSigma, itsNSigma;
      switch (pid) {
        case kPions: // For Pions
          tofNSigma = track.tofNSigmaPi();
          if (!useITS) {
            tpcNSigma = track.tpcNSigmaPi();
            tpcExpSigma = track.tpcExpSigmaPi();
          } else {
            itsNSigma = itsResponse.nSigmaITS<o2::track::PID::Pion>(track);
          }
          break;
        case kKaons: // For Kaons
          tofNSigma = track.tofNSigmaKa();
          if (!useITS) {
            tpcNSigma = track.tpcNSigmaKa();
            tpcExpSigma = track.tpcExpSigmaKa();
          } else {
            itsNSigma = itsResponse.nSigmaITS<o2::track::PID::Kaon>(track);
          }
          break;

        case kProtons: // For Protons
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

      if (cfgPIDConfigs.cfgQANsigma) {
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
    int fSampleIndex = static_cast<int>(gRandom->Uniform(0.0, cfgSampleSize));

    float triggerWeight = 1.0f;
    // loop over all tracks
    for (auto const& track1 : tracks1) {

      if (!trackSelected(track1)) {
        continue;
      }

      fillTrackQA(track1, kCharged, Before, triggerWeight);
      fillTrackQA(track1, kPions, Before, triggerWeight);
      fillTrackQA(track1, kKaons, Before, triggerWeight);
      fillTrackQA(track1, kProtons, Before, triggerWeight);

      int pidIndex = getNsigmaPID(track1);

      if (pidIndex > kCharged) {
        fillTrackQA(track1, pidIndex, After, triggerWeight);
      }

      if (!getEfficiencyCorrection(triggerWeight, track1.pt(), track1.eta(), posZ)) {
        continue;
      }

      if (system == SameEvent) {
        if (corType == kFT0C) {
          histos.fill(HIST("Trig_hist_TPC_FT0C"), fSampleIndex, posZ, eventClass, track1.pt(), kCharged, eventWeight * triggerWeight);
          if (pidIndex > kCharged) {
            histos.fill(HIST("Trig_hist_TPC_FT0C"), fSampleIndex, posZ, eventClass, track1.pt(), pidIndex, eventWeight * triggerWeight);
          }
        } else if (corType == kFT0A) {
          histos.fill(HIST("Trig_hist_TPC_FT0A"), fSampleIndex, posZ, eventClass, track1.pt(), kCharged, eventWeight * triggerWeight);
          if (pidIndex > kCharged) {
            histos.fill(HIST("Trig_hist_TPC_FT0A"), fSampleIndex, posZ, eventClass, track1.pt(), pidIndex, eventWeight * triggerWeight);
          }
        }
      }

      std::size_t channelSize = 0;
      if (corType == kFT0C) {
        channelSize = ft0.channelC().size();
      } else if (corType == kFT0A) {
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
          if (corType == kFT0A) {
            if (cfgQABasic) {
              histos.fill(HIST("Assoc_amp_same_TPC_FT0A"), chanelid, ampl);
              histos.fill(HIST("deltaEta_deltaPhi_same_TPC_FT0A"), deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
            }
            sameTpcFt0a->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, deltaPhi, deltaEta, track1.pt(), kCharged, ampl * eventWeight * triggerWeight);
            if (pidIndex > kCharged) {
              sameTpcFt0a->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, deltaPhi, deltaEta, track1.pt(), pidIndex, ampl * eventWeight * triggerWeight);
            }
          } else if (corType == kFT0C) {
            if (cfgQABasic) {
              histos.fill(HIST("Assoc_amp_same_TPC_FT0C"), chanelid, ampl);
              histos.fill(HIST("deltaEta_deltaPhi_same_TPC_FT0C"), deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
            }
            sameTpcFt0c->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, deltaPhi, deltaEta, track1.pt(), kCharged, ampl * eventWeight * triggerWeight);
            if (pidIndex > kCharged) {
              sameTpcFt0c->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, deltaPhi, deltaEta, track1.pt(), pidIndex, ampl * eventWeight * triggerWeight);
            }
          }
        } else if (system == MixedEvent) {
          if (corType == kFT0A) {
            if (cfgQABasic) {
              histos.fill(HIST("Assoc_amp_mixed_TPC_FT0A"), chanelid, ampl);
              histos.fill(HIST("deltaEta_deltaPhi_mixed_TPC_FT0A"), deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
            }
            mixedTpcFt0a->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, deltaPhi, deltaEta, track1.pt(), kCharged, ampl * eventWeight * triggerWeight);
            if (pidIndex > kCharged) {
              mixedTpcFt0a->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, deltaPhi, deltaEta, track1.pt(), pidIndex, ampl * eventWeight * triggerWeight);
            }
          } else if (corType == kFT0C) {
            if (cfgQABasic) {
              histos.fill(HIST("Assoc_amp_mixed_TPC_FT0C"), chanelid, ampl);
              histos.fill(HIST("deltaEta_deltaPhi_mixed_TPC_FT0C"), deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
            }
            mixedTpcFt0c->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, deltaPhi, deltaEta, track1.pt(), kCharged, ampl * eventWeight * triggerWeight);
            if (pidIndex > kCharged) {
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
    int fSampleIndex = static_cast<int>(gRandom->Uniform(0.0, cfgSampleSize));

    float triggerWeight = 1.0f;
    std::size_t channelASize = ft0Col1.channelA().size();
    std::size_t channelCSize = ft0Col2.channelC().size();
    // loop over all tracks
    for (std::size_t iChA = 0; iChA < channelASize; iChA++) {

      int chanelAid = 0;
      float amplA = 0.;
      getChannel(ft0Col1, iChA, chanelAid, amplA, kFT0A, system);
      auto phiA = getPhiFT0(chanelAid, kFT0A);
      auto etaA = getEtaFT0(chanelAid, kFT0A);

      if (system == SameEvent) {
        histos.fill(HIST("Trig_hist_FT0A_FT0C"), fSampleIndex, posZ, eventClass, eventWeight * amplA);
      }

      for (std::size_t iChC = 0; iChC < channelCSize; iChC++) {
        int chanelCid = 0;
        float amplC = 0.;
        getChannel(ft0Col2, iChC, chanelCid, amplC, kFT0C, system);
        auto phiC = getPhiFT0(chanelCid, kFT0C);
        auto etaC = getEtaFT0(chanelCid, kFT0C);
        float deltaPhi = RecoDecay::constrainAngle(phiA - phiC, -PIHalf);
        float deltaEta = etaA - etaC;

        // fill the right sparse and histograms
        if (system == SameEvent) {
          if (cfgQABasic) {
            histos.fill(HIST("deltaEta_deltaPhi_same_FT0A_FT0C"), deltaPhi, deltaEta, amplA * amplC * eventWeight * triggerWeight);
          }
          sameFt0aFt0c->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, deltaPhi, deltaEta, amplA * amplC * eventWeight * triggerWeight);
        } else if (system == MixedEvent) {
          if (cfgQABasic) {
            histos.fill(HIST("deltaEta_deltaPhi_mixed_FT0A_FT0C"), deltaPhi, deltaEta, amplA * amplC * eventWeight * triggerWeight);
          }
          mixedFt0aFt0c->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, deltaPhi, deltaEta, amplA * amplC * eventWeight * triggerWeight);
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

  void processSameTpcFt0(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
  {
    if (cfgQABasic) {
      eventSelectedIndividually(collision);
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

    auto cent = getCentrality(collision);
    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), cent, true)) {
      return;
    }

    if (!collision.has_foundFT0()) {
      return;
    }

    histos.fill(HIST("hEventCount"), kHaveFT0Cut);

    loadAlignParam(bc.timestamp());
    loadGain(bc);
    loadCorrection(bc.timestamp());
    float weightCent = 1.0f;

    getCentralityWeight(weightCent, cent);
    histos.fill(HIST("Centrality"), cent);
    histos.fill(HIST("CentralityWeighted"), cent, weightCent);

    histos.fill(HIST("zVtx"), collision.posZ());

    histos.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin

    auto multiplicity = static_cast<double>(tracks.size());

    if (cfgQABasic) {
      histos.fill(HIST("Nch"), multiplicity);
    }

    if (cfgStrictTrackCounter) {
      trackCounter(tracks, multiplicity);
    }

    if (cfgQABasic) {
      histos.fill(HIST("Nch_corrected"), multiplicity);
    }

    if (cfgSelCollByNch && (multiplicity > cfgMaxMultForCorrelations || multiplicity < cfgMinMultForCorrelations)) {
      return;
    }
    if (!cfgSelCollByNch && (cent > cfgMaxCentForCorrelations || cent < cfgMinCentForCorrelations)) {
      return;
    }

    auto eventClass = (cfgSelCollByNch) ? multiplicity : cent;

    const auto& ft0 = collision.foundFT0();
    if (cfgAnalyzeTPCFT0A) {
      fillCorrelationsTPCFT0(tracks, ft0, collision.posZ(), SameEvent, eventClass, kFT0A, weightCent);
    }
    if (cfgAnalyzeTPCFT0C) {
      fillCorrelationsTPCFT0(tracks, ft0, collision.posZ(), SameEvent, eventClass, kFT0C, weightCent);
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
    Pair<FilteredCollisions, FilteredTracks, FilteredTracks, MixedBinning> pairs{binningOnVtxAndMult, cfgMinMixEventNum, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
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

      if (cfgUseAdditionalEventCut && !eventSelected(collision1, tracks1.size(), cent1, false)) {
        continue;
      }
      if (cfgUseAdditionalEventCut && !eventSelected(collision2, tracks2.size(), cent2, false)) {
        continue;
      }

      if (!(collision1.has_foundFT0() && collision2.has_foundFT0())) {
        continue;
      }

      histos.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      int currentRunNumber = bc.runNumber();
      if (!cfgRunRemoveList.value.empty()) {
        if (!isGoodRun(currentRunNumber)) { // Rejects runs if bad run number
          continue;
        }
      }

      loadAlignParam(bc.timestamp());
      loadCorrection(bc.timestamp());

      auto multiplicity = static_cast<double>(tracks1.size());

      if (cfgStrictTrackCounter) {
        trackCounter(tracks1, multiplicity);
      }

      if (cfgSelCollByNch && (multiplicity > cfgMaxMultForCorrelations || multiplicity < cfgMinMultForCorrelations)) {
        continue;
      }
      if (!cfgSelCollByNch && (cent1 > cfgMaxCentForCorrelations || cent1 < cfgMinCentForCorrelations)) {
        continue;
      }

      auto eventClass = (cfgSelCollByNch) ? multiplicity : cent1;

      float eventWeight = 1.0f;

      if (cfgUseEventWeights) {
        eventWeight = 1.0f / it.currentWindowNeighbours();
      }
      float weightCent = 1.0f;
      getCentralityWeight(weightCent, cent1);

      const auto& ft0 = collision2.foundFT0();
      if (cfgAnalyzeTPCFT0A) {
        fillCorrelationsTPCFT0(tracks1, ft0, collision1.posZ(), MixedEvent, eventClass, kFT0A, eventWeight * weightCent);
      }
      if (cfgAnalyzeTPCFT0C) {
        fillCorrelationsTPCFT0(tracks1, ft0, collision1.posZ(), MixedEvent, eventClass, kFT0C, eventWeight * weightCent);
      }
    }
  }
  PROCESS_SWITCH(PidLongRange, processMixedTpcFt0, "Process mixed events for TPC-FT0A correlation", false);

  void processSameFt0aFt0c(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
  {

    if (cfgQABasic) {
      eventSelectedIndividually(collision);
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

    auto cent = getCentrality(collision);
    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), cent, true)) {
      return;
    }

    if (!collision.has_foundFT0()) {
      return;
    }
    histos.fill(HIST("hEventCount"), kHaveFT0Cut);

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

    if (cfgQABasic) {
      histos.fill(HIST("Nch"), multiplicity);
    }

    if (cfgStrictTrackCounter) {
      trackCounter(tracks, multiplicity);
    }

    if (cfgQABasic) {
      histos.fill(HIST("Nch_corrected"), multiplicity);
    }

    if (cfgSelCollByNch && (multiplicity > cfgMaxMultForCorrelations || multiplicity < cfgMinMultForCorrelations)) {
      return;
    }
    if (!cfgSelCollByNch && (cent > cfgMaxCentForCorrelations || cent < cfgMinCentForCorrelations)) {
      return;
    }

    auto eventClass = (cfgSelCollByNch) ? multiplicity : cent;

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

      auto cent1 = getCentrality(collision1);
      auto cent2 = getCentrality(collision2);

      if (cfgUseAdditionalEventCut && !eventSelected(collision1, tracks1.size(), cent1, false)) {
        continue;
      }
      if (cfgUseAdditionalEventCut && !eventSelected(collision2, tracks2.size(), cent2, false)) {
        continue;
      }

      if (!(collision1.has_foundFT0() && collision2.has_foundFT0())) {
        continue;
      }

      histos.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      int currentRunNumber = bc.runNumber();
      if (!cfgRunRemoveList.value.empty()) {
        if (!isGoodRun(currentRunNumber)) { // Rejects runs if bad run number
          continue;
        }
      }
      loadAlignParam(bc.timestamp());
      loadCorrection(bc.timestamp());

      auto multiplicity = static_cast<double>(tracks1.size());

      if (cfgStrictTrackCounter) {
        trackCounter(tracks1, multiplicity);
      }

      if (cfgSelCollByNch && (multiplicity > cfgMaxMultForCorrelations || multiplicity < cfgMinMultForCorrelations)) {
        continue;
      }
      if (!cfgSelCollByNch && (cent1 > cfgMaxCentForCorrelations || cent1 < cfgMinCentForCorrelations)) {
        continue;
      }

      auto eventClass = (cfgSelCollByNch) ? multiplicity : cent1;

      float eventWeight = 1.0f;

      if (cfgUseEventWeights) {
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
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PidLongRange>(cfgc),
  };
}
