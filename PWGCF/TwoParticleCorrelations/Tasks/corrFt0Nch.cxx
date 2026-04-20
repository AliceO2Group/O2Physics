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

/// \file corrFt0Nch.cxx
/// \brief Ultra long range correlation using forward FIT detectors and TPC, with foxus on multiplicity dependence
/// \author Thor Jensen (thor.kjaersgaard.jensen@cern.ch)

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGMM/Mult/DataModel/bestCollisionTable.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsCommonDataFormats/AlignParam.h>
#include <FT0Base/Geometry.h>
#include <FV0Base/Geometry.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/StepTHn.h>
#include <Framework/runDataProcessing.h>

#include <TF1.h>
#include <TFile.h>
#include <TH3.h>
#include <TRandom3.h>

#include <algorithm>
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
using namespace constants::math;

// define the filtered collisions and tracks
#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct CorrFt0Nch {

  Service<ccdb::BasicCCDBManager> ccdb;

  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, false, "Use additional event cut on mult correlations")
  O2_DEFINE_CONFIGURABLE(cfgZVtxCut, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgUseTransverseMomentum, bool, false, "Use transverse momentum for correlation container")
  struct : ConfigurableGroup{
             O2_DEFINE_CONFIGURABLE(cfgPtCutMin, float, 0.2f, "minimum accepted track pT")
               O2_DEFINE_CONFIGURABLE(cfgPtCutMax, float, 10.0f, "maximum accepted track pT")
                 O2_DEFINE_CONFIGURABLE(cfgEtaCut, float, 0.8f, "Eta cut")
                   O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5f, "max chi2 per TPC clusters")
                     O2_DEFINE_CONFIGURABLE(cfgCutTPCclu, float, 50.0f, "minimum TPC clusters")
                       O2_DEFINE_CONFIGURABLE(cfgCutTPCCrossedRows, float, 70.0f, "minimum TPC crossed rows")
                         O2_DEFINE_CONFIGURABLE(cfgCutITSclu, float, 5.0f, "minimum ITS clusters")
                           O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "max DCA to vertex z")

           } cfgTrackCuts;

  struct : ConfigurableGroup{
             O2_DEFINE_CONFIGURABLE(cfgMinMult, int, 0, "Minimum multiplicity for collision")
               O2_DEFINE_CONFIGURABLE(cfgMaxMult, int, 10, "Maximum multiplicity for collision")
                 O2_DEFINE_CONFIGURABLE(cfgEvSelkNoSameBunchPileup, bool, false, "rejects collisions which are associated with the same found-by-T0 bunch crossing")
                   O2_DEFINE_CONFIGURABLE(cfgEvSelkNoITSROFrameBorder, bool, false, "reject events at ITS ROF border")
                     O2_DEFINE_CONFIGURABLE(cfgEvSelkNoTimeFrameBorder, bool, false, "reject events at TF border")
                       O2_DEFINE_CONFIGURABLE(cfgEvSelkIsGoodZvtxFT0vsPV, bool, false, "removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference, use this cut at low multiplicities with caution")
                         O2_DEFINE_CONFIGURABLE(cfgEvSelkNoCollInTimeRangeStandard, bool, false, "no collisions in specified time range")
                           O2_DEFINE_CONFIGURABLE(cfgEvSelkIsGoodITSLayer0123, bool, true, "cut time intervals with dead ITS layers 0,1,2,3")
                             O2_DEFINE_CONFIGURABLE(cfgEvSelkIsGoodITSLayersAll, bool, true, "cut time intervals with dead ITS staves")
                               O2_DEFINE_CONFIGURABLE(cfgEvSelkNoCollInRofStandard, bool, false, "no other collisions in this Readout Frame with per-collision multiplicity above threshold")
                                 O2_DEFINE_CONFIGURABLE(cfgEvSelkNoHighMultCollInPrevRof, bool, false, "veto an event if FT0C amplitude in previous ITS ROF is above threshold")
                                   O2_DEFINE_CONFIGURABLE(cfgEvSelMultCorrelation, bool, true, "Multiplicity correlation cut")
                                     O2_DEFINE_CONFIGURABLE(cfgEvSelV0AT0ACut, bool, true, "V0A T0A 5 sigma cut")
                                       O2_DEFINE_CONFIGURABLE(cfgEvSelOccupancy, bool, true, "Occupancy cut")
                                         O2_DEFINE_CONFIGURABLE(cfgCutOccupancyHigh, int, 2000, "High cut on TPC occupancy")
                                           O2_DEFINE_CONFIGURABLE(cfgCutOccupancyLow, int, 0, "Low cut on TPC occupancy")

           } cfgEventSelection;

  struct : ConfigurableGroup{
             O2_DEFINE_CONFIGURABLE(cfgSystematicsVariation, bool, false, "Enable systematics variation for track cuts")
               O2_DEFINE_CONFIGURABLE(cfgSystematicsCutChi2prTPCcls, float, 3.0f, "max chi2 per TPC clusters for systematics variation")
                 O2_DEFINE_CONFIGURABLE(cfgSystematicsCutTPCclu, float, 40.0f, "minimum TPC clusters for systematics variation")
                   O2_DEFINE_CONFIGURABLE(cfgSystematicsCutTPCCrossedRows, float, 60.0f, "minimum TPC crossed rows for systematics variation")
                     O2_DEFINE_CONFIGURABLE(cfgSystematicsCutITSclu, float, 4.0f, "minimum ITS clusters for systematics variation")
                       O2_DEFINE_CONFIGURABLE(cfgSystematicsCutDCAz, float, 1.5f, "max DCA to vertex z for systematics variation")} cfgSystematics;

  O2_DEFINE_CONFIGURABLE(cfgMinMixEventNum, int, 5, "Minimum number of events to mix")
  O2_DEFINE_CONFIGURABLE(cfgMergingCut, float, 0.02, "Merging cut on track merge")
  O2_DEFINE_CONFIGURABLE(cfgApplyTwoTrackEfficiency, bool, true, "Apply two track efficiency for tpc tpc")
  O2_DEFINE_CONFIGURABLE(cfgRadiusLow, float, 0.8, "Low radius for merging cut")
  O2_DEFINE_CONFIGURABLE(cfgRadiusHigh, float, 2.5, "High radius for merging cut")
  O2_DEFINE_CONFIGURABLE(cfgSampleSize, double, 10, "Sample size for mixed event")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgCentralityWeight, std::string, "", "CCDB path to centrality weight object")
  O2_DEFINE_CONFIGURABLE(cfgLocalEfficiency, bool, false, "Use local efficiency object")
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
  } cfgFuncParas;

  Configurable<float> cfgCutFV0{"cfgCutFV0", 50., "FV0A threshold"};
  Configurable<float> cfgCutFT0A{"cfgCutFT0A", 150., "FT0A threshold"};
  Configurable<float> cfgCutFT0C{"cfgCutFT0C", 50., "FT0C threshold"};
  Configurable<float> cfgCutZDC{"cfgCutZDC", 10., "ZDC threshold"};

  SliceCache cache;

  ConfigurableAxis axisVertex{"axisVertex", {10, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisMult{"axisMult", {10, 0, 100}, "multiplicity axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {72, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEtaTpcFt0a{"axisDeltaEtaTpcFt0a", {32, -5.8, -2.6}, "delta eta axis, -5.8~-2.6 for TPC-FT0A,"};
  ConfigurableAxis axisDeltaEtaTpcFt0c{"axisDeltaEtaTpcFt0c", {32, 1.2, 4.2}, "delta eta axis, 1.2~4.2 for TPC-FT0C"};
  ConfigurableAxis axisDeltaEtaFt0aFt0c{"axisDeltaEtaFt0aFt0c", {32, -1.5, 3.0}, "delta eta axis"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt associated axis for histograms"};
  ConfigurableAxis axisVtxMix{"axisVtxMix", {VARIABLE_WIDTH, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, "vertex axis for mixed event histograms"};
  ConfigurableAxis axisMultMix{"axisMultMix", {VARIABLE_WIDTH, 0, 10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260}, "multiplicity / centrality axis for mixed event histograms"};
  ConfigurableAxis axisSample{"axisSample", {cfgSampleSize, 0, cfgSampleSize}, "sample axis for histograms"};

  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt axis for efficiency histograms"};
  ConfigurableAxis axisAmplitudeFt0a{"axisAmplitudeFt0a", {5000, 0, 1000}, "FT0A amplitude"};
  ConfigurableAxis axisChannelFt0aAxis{"axisChannelFt0aAxis", {96, 0.0, 96.0}, "FT0A channel"};

  Configurable<std::string> cfgGainEqPath{"cfgGainEqPath", "Analysis/EventPlane/GainEq", "CCDB path for gain equalization constants"};
  Configurable<int> cfgCorrLevel{"cfgCorrLevel", 0, "calibration step: 0 = no corr, 1 = gain corr"};
  ConfigurableAxis cfgaxisFITamp{"cfgaxisFITamp", {1000, 0, 5000}, ""};
  AxisSpec axisFit{cfgaxisFITamp, "fit amplitude"};
  AxisSpec axisChID = {220, 0, 220};

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgZVtxCut);
  Filter trackFilter = (nabs(aod::track::eta) < cfgTrackCuts.cfgEtaCut) && (cfgTrackCuts.cfgPtCutMin < aod::track::pt) && (cfgTrackCuts.cfgPtCutMax > aod::track::pt) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t)true)) && (aod::track::tpcChi2NCl < cfgTrackCuts.cfgCutChi2prTPCcls) && (aod::track::dcaZ < cfgTrackCuts.cfgCutDCAz);

  using FilteredCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSel, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::Mults>>;
  using FilteredTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA>>;

  // FT0 geometry
  o2::ft0::Geometry ft0Det;
  static constexpr uint64_t Ft0IndexA = 96;
  std::vector<o2::detectors::AlignParam>* offsetFT0;
  std::vector<float> cstFT0RelGain{};

  // Corrections
  TH3D* mEfficiency = nullptr;
  TH1D* mCentralityWeight = nullptr;
  bool correctionsLoaded = false;

  // Define the outputs
  OutputObj<CorrelationContainer> sameTpcFt0a{"sameEvent_TPC_FT0A"};
  OutputObj<CorrelationContainer> mixedTpcFt0a{"mixedEvent_TPC_FT0A"};
  OutputObj<CorrelationContainer> sameTpcFt0c{"sameEvent_TPC_FT0C"};
  OutputObj<CorrelationContainer> mixedTpcFt0c{"mixedEvent_TPC_FT0C"};
  OutputObj<CorrelationContainer> sameFt0aFt0c{"sameEvent_FT0A_FT0C"};
  OutputObj<CorrelationContainer> mixedFt0aFt0c{"mixedEvent_FA_FT0C"};
  OutputObj<CorrelationContainer> sameTPC{"sameEvent_TPC"};
  OutputObj<CorrelationContainer> mixedTPC{"mixedEvent_TPC"};

  HistogramRegistry registry{"registry"};

  // define global variables
  TRandom3* gRandom = new TRandom3();

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
    kUseVertexITSTPC,
    kUseTVXinTRD,
    kNEventCuts
  };

  enum EventType {
    SameEvent = 1,
    MixedEvent = 3
  };

  enum FITIndex {
    kFT0A = 0,
    kFT0C = 1
  };

  enum DetectorChannels {
    kFT0AInnerRingMin = 0,
    kFT0AInnerRingMax = 31,
    kFT0AOuterRingMin = 32,
    kFT0AOuterRingMax = 95,
    kFT0CInnerRingMin = 96,
    kFT0CInnerRingMax = 143,
    kFT0COuterRingMin = 144,
    kFT0COuterRingMax = 207
  };

  std::array<std::array<int, 1>, 16> eventCuts;

  void init(InitContext&)
  {

    const AxisSpec axisPhi{72, 0.0, constants::math::TwoPI, "#varphi"};
    const AxisSpec axisEta{40, -1., 1., "#eta"};
    const AxisSpec axisEtaFull{90, -4., 5., "#eta"};

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    LOGF(info, "Starting init");

    if (doprocessSameFt0aFt0c || doprocessSameTpcFt0a || doprocessSameTpcFt0c || doprocessSameTPC) {
      registry.add("hEventCountSpecific", "Number of Event;; Count", {HistType::kTH1D, {{13, 0, 13}}});
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(1, "after sel8");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(2, "kNoSameBunchPileup");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(3, "kNoITSROFrameBorder");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(4, "kNoTimeFrameBorder");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(5, "kIsGoodZvtxFT0vsPV");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(6, "kNoCollInTimeRangeStandard");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(7, "kIsGoodITSLayer0123");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(8, "kIsGoodITSLayersAll");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(9, "kNoCollInRofStandard");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(10, "kNoHighMultCollInPrevRof");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(11, "occupancy");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(12, "MultCorrelation");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(13, "cfgEvSelV0AT0ACut");
    }

    if (doprocessSameTPC || doprocessSameFt0aFt0c || doprocessSameTpcFt0a || doprocessSameTpcFt0c) {
      registry.add("Phi", "Phi", {HistType::kTH1D, {axisPhi}});
      registry.add("Eta", "Eta", {HistType::kTH1D, {axisEta}});
      registry.add("EtaCorrected", "EtaCorrected", {HistType::kTH1D, {axisEta}});
      registry.add("pT", "pT", {HistType::kTH1D, {axisPtTrigger}});
      registry.add("pTCorrected", "pTCorrected", {HistType::kTH1D, {axisPtTrigger}});
      registry.add("Nch", "N_{ch}", {HistType::kTH1D, {axisMult}});
      registry.add("zVtx", "zVtx", {HistType::kTH1D, {axisVertex}});
      if (doprocessSameFt0aFt0c || doprocessSameTpcFt0a || doprocessSameTpcFt0c) {
        registry.add("FT0Amp", "", {HistType::kTH2F, {axisChID, axisFit}});
        registry.add("FT0AmpCorrect", "", {HistType::kTH2F, {axisChID, axisFit}});
      }
    }

    if (doprocessSameTpcFt0a) {
      registry.add("deltaEta_deltaPhi_same_TPC_FT0A", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0a}}); // check to see the delta eta and delta phi distribution
      registry.add("deltaEta_deltaPhi_mixed_TPC_FT0A", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0a}});
      registry.add("Assoc_amp_same_TPC_FT0A", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
      registry.add("Assoc_amp_mixed_TPC_FT0A", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
      registry.add("Trig_hist_TPC_FT0A", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
    }
    if (doprocessSameTpcFt0c) {
      registry.add("deltaEta_deltaPhi_same_TPC_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0c}}); // check to see the delta eta and delta phi distribution
      registry.add("deltaEta_deltaPhi_mixed_TPC_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0c}});
      registry.add("Assoc_amp_same_TPC_FT0C", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
      registry.add("Assoc_amp_mixed_TPC_FT0C", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
      registry.add("Trig_hist_TPC_FT0C", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
    }
    if (doprocessSameFt0aFt0c) {
      registry.add("deltaEta_deltaPhi_same_FT0A_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaFt0aFt0c}}); // check to see the delta eta and delta phi distribution
      registry.add("deltaEta_deltaPhi_mixed_FT0A_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaFt0aFt0c}});
      registry.add("Trig_hist_FT0A_FT0C", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
    }
    if (doprocessSameTPC) {
      registry.add("deltaEta_deltaPhi_same_TPC", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0a}}); // check to see the delta eta and delta phi distribution
      registry.add("deltaEta_deltaPhi_mixed_TPC", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0a}});
      registry.add("Trig_hist_TPC", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
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
                                         {axisDeltaEtaTpcFt0a, "#Delta#eta"}}; // use the same delta eta axis for TPC-TPC correlation

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

    LOGF(info, "End of init");
  }

  template <typename TCollision>
  bool eventSelected(TCollision collision, const int multTrk, const bool fillCounter)
  {
    registry.fill(HIST("hEventCountSpecific"), 0.5);
    if (cfgEventSelection.cfgEvSelkNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return 0;
    }
    if (fillCounter && cfgEventSelection.cfgEvSelkNoSameBunchPileup)
      registry.fill(HIST("hEventCountSpecific"), 1.5);
    if (cfgEventSelection.cfgEvSelkNoITSROFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return 0;
    }
    if (fillCounter && cfgEventSelection.cfgEvSelkNoITSROFrameBorder)
      registry.fill(HIST("hEventCountSpecific"), 2.5);
    if (cfgEventSelection.cfgEvSelkNoTimeFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return 0;
    }
    if (fillCounter && cfgEventSelection.cfgEvSelkNoTimeFrameBorder)
      registry.fill(HIST("hEventCountSpecific"), 3.5);
    if (cfgEventSelection.cfgEvSelkIsGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return 0;
    }
    if (fillCounter && cfgEventSelection.cfgEvSelkIsGoodZvtxFT0vsPV)
      registry.fill(HIST("hEventCountSpecific"), 4.5);
    if (cfgEventSelection.cfgEvSelkNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      // no collisions in specified time range
      return 0;
    }
    if (fillCounter && cfgEventSelection.cfgEvSelkNoCollInTimeRangeStandard)
      registry.fill(HIST("hEventCountSpecific"), 5.5);

    if (cfgEventSelection.cfgEvSelkIsGoodITSLayer0123 && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
      // from Jan 9 2025 AOT meeting
      // cut time intervals with dead ITS staves
      return 0;
    }
    if (fillCounter && cfgEventSelection.cfgEvSelkIsGoodITSLayer0123)
      registry.fill(HIST("hEventCountSpecific"), 6.5);

    if (cfgEventSelection.cfgEvSelkIsGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      // from Jan 9 2025 AOT meeting
      // cut time intervals with dead ITS staves
      return 0;
    }

    if (fillCounter && cfgEventSelection.cfgEvSelkIsGoodITSLayersAll)
      registry.fill(HIST("hEventCountSpecific"), 7.5);

    if (cfgEventSelection.cfgEvSelkNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      // no other collisions in this Readout Frame with per-collision multiplicity above threshold
      return 0;
    }
    if (fillCounter && cfgEventSelection.cfgEvSelkNoCollInRofStandard)
      registry.fill(HIST("hEventCountSpecific"), 8.5);
    if (cfgEventSelection.cfgEvSelkNoHighMultCollInPrevRof && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      // veto an event if FT0C amplitude in previous ITS ROF is above threshold
      return 0;
    }
    if (fillCounter && cfgEventSelection.cfgEvSelkNoHighMultCollInPrevRof)
      registry.fill(HIST("hEventCountSpecific"), 9.5);
    auto occupancy = collision.trackOccupancyInTimeRange();
    if (cfgEventSelection.cfgEvSelOccupancy && (occupancy < cfgEventSelection.cfgCutOccupancyLow || occupancy > cfgEventSelection.cfgCutOccupancyHigh))
      return 0;
    if (fillCounter && cfgEventSelection.cfgEvSelOccupancy)
      registry.fill(HIST("hEventCountSpecific"), 10.5);

    auto multNTracksPV = collision.multNTracksPV();

    if (cfgFuncParas.cfgMultGlobalPVCutEnabled) {
      if (multTrk < cfgFuncParas.fMultGlobalPVCutLow->Eval(multNTracksPV))
        return 0;
      if (multTrk > cfgFuncParas.fMultGlobalPVCutHigh->Eval(multNTracksPV))
        return 0;
    }
    if (cfgFuncParas.cfgMultMultV0ACutEnabled) {
      if (collision.multFV0A() < cfgFuncParas.fMultMultV0ACutLow->Eval(multTrk))
        return 0;
      if (collision.multFV0A() > cfgFuncParas.fMultMultV0ACutHigh->Eval(multTrk))
        return 0;
    }

    if (fillCounter && cfgEventSelection.cfgEvSelMultCorrelation)
      registry.fill(HIST("hEventCountSpecific"), 11.5);

    // V0A T0A 5 sigma cut
    float sigma = 5.0;
    if (cfgEventSelection.cfgEvSelV0AT0ACut && (std::fabs(collision.multFV0A() - cfgFuncParas.fT0AV0AMean->Eval(collision.multFT0A())) > sigma * cfgFuncParas.fT0AV0ASigma->Eval(collision.multFT0A())))
      return 0;
    if (fillCounter && cfgEventSelection.cfgEvSelV0AT0ACut)
      registry.fill(HIST("hEventCountSpecific"), 12.5);

    return 1;
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
  bool trackSelected(TTrack track)
  {
    return ((track.tpcNClsFound() >= cfgTrackCuts.cfgCutTPCclu) && (track.tpcNClsCrossedRows() >= cfgTrackCuts.cfgCutTPCCrossedRows) && (track.itsNCls() >= cfgTrackCuts.cfgCutITSclu));
  }

  template <typename TTrack>
  bool trackSelectedSystematics(TTrack track)
  {
    return ((track.tpcNClsFound() >= cfgSystematics.cfgSystematicsCutTPCclu) && (track.tpcNClsCrossedRows() >= cfgSystematics.cfgSystematicsCutTPCCrossedRows) && (track.itsNCls() >= cfgSystematics.cfgSystematicsCutITSclu));
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

    if (dPhiStar > constants::math::PI)
      dPhiStar = constants::math::TwoPI - dPhiStar;
    if (dPhiStar < -constants::math::PI)
      dPhiStar = -constants::math::TwoPI - dPhiStar;

    return dPhiStar;
  }

  void loadGain(aod::BCsWithTimestamps::iterator const& bc)
  {
    cstFT0RelGain.clear();
    cstFT0RelGain = {};
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
      if (system == SameEvent)
        registry.fill(HIST("FT0Amp"), id, ampl);
      ampl = ampl / cstFT0RelGain[id];
      if (system == SameEvent) {
        registry.fill(HIST("FT0AmpCorrect"), id, ampl);
      }
    } else if (fitType == kFT0A) {
      id = ft0.channelA()[iCh];
      ampl = ft0.amplitudeA()[iCh];
      if (system == SameEvent)
        registry.fill(HIST("FT0Amp"), id, ampl);
      ampl = ampl / cstFT0RelGain[id];
      if (system == SameEvent) {
        registry.fill(HIST("FT0AmpCorrect"), id, ampl);
      }
    } else {
      LOGF(fatal, "Cor Index %d out of range", fitType);
    }
  }

  void loadCorrection(uint64_t timestamp)
  {
    if (correctionsLoaded) {
      return;
    }
    if (cfgEfficiency.value.empty() == false) {
      if (cfgLocalEfficiency > 0) {
        TFile* fEfficiencyTrigger = TFile::Open(cfgEfficiency.value.c_str(), "READ");
        mEfficiency = reinterpret_cast<TH3D*>(fEfficiencyTrigger->Get("ccdb_object"));
      } else {
        mEfficiency = ccdb->getForTimeStamp<TH3D>(cfgEfficiency, timestamp);
      }
      if (mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEfficiency.value.c_str(), (void*)mEfficiency);
    }
    if (cfgCentralityWeight.value.empty() == false) {
      mCentralityWeight = ccdb->getForTimeStamp<TH1D>(cfgCentralityWeight, timestamp);
      if (mCentralityWeight == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgCentralityWeight.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgCentralityWeight.value.c_str(), (void*)mCentralityWeight);
    }
    correctionsLoaded = true;
  }

  bool getEfficiencyCorrection(float& weight_nue, float eta, float pt, float posZ)
  {
    float eff = 1.;
    if (mEfficiency) {
      int etaBin = mEfficiency->GetXaxis()->FindBin(eta);
      int ptBin = mEfficiency->GetYaxis()->FindBin(pt);
      int zBin = mEfficiency->GetZaxis()->FindBin(posZ);
      eff = mEfficiency->GetBinContent(etaBin, ptBin, zBin);
    } else {
      eff = 1.0;
    }
    if (eff == 0)
      return false;
    weight_nue = 1. / eff;
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
  void trackCounter(TTracks tracks, int& multiplicity) // function to count the number of tracks in the event and fill the histogram
  {
    int mult = 0;
    for (auto const& track : tracks) {
      if (!trackSelected(track))
        continue;
      mult++;
    }
    multiplicity = mult;
  }

  template <CorrelationContainer::CFStep step, typename TTracks, typename TFT0s>
  void fillCorrelationsTPCFT0(TTracks tracks1, TFT0s const& ft0, float posZ, int system, int corType, float multiplicity, float eventWeight) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {

    if (system == SameEvent) {
      registry.fill(HIST("Nch"), multiplicity);
    }

    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);

    float triggerWeight = 1.0f;
    // loop over all tracks
    for (auto const& track1 : tracks1) {

      if (!trackSelected(track1))
        continue;
      if (!getEfficiencyCorrection(triggerWeight, track1.eta(), track1.pt(), posZ))
        continue;

      if (system == SameEvent) {
        if (corType == kFT0C) {
          registry.fill(HIST("Trig_hist_TPC_FT0C"), fSampleIndex, posZ, track1.pt(), eventWeight * triggerWeight);
        } else if (corType == kFT0A) {
          registry.fill(HIST("Trig_hist_TPC_FT0A"), fSampleIndex, posZ, track1.pt(), eventWeight * triggerWeight);
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
            registry.fill(HIST("Assoc_amp_same_TPC_FT0A"), chanelid, ampl);
            registry.fill(HIST("deltaEta_deltaPhi_same_TPC_FT0A"), deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
            sameTpcFt0a->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), multiplicity, deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
          } else if (corType == kFT0C) {
            registry.fill(HIST("Assoc_amp_same_TPC_FT0C"), chanelid, ampl);
            registry.fill(HIST("deltaEta_deltaPhi_same_TPC_FT0C"), deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
            sameTpcFt0c->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), multiplicity, deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
          }
        } else if (system == MixedEvent) {
          if (corType == kFT0A) {
            registry.fill(HIST("Assoc_amp_mixed_TPC_FT0A"), chanelid, ampl);
            registry.fill(HIST("deltaEta_deltaPhi_mixed_TPC_FT0A"), deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
            mixedTpcFt0a->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), multiplicity, deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
          } else if (corType == kFT0C) {
            registry.fill(HIST("Assoc_amp_mixed_TPC_FT0C"), chanelid, ampl);
            registry.fill(HIST("deltaEta_deltaPhi_mixed_TPC_FT0C"), deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
            mixedTpcFt0c->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), multiplicity, deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
          }
        }
      }
    }
  }

  template <typename TCollision, typename TTracks>
  void fillYield(TCollision collision, TTracks tracks) // function to fill the yield and etaphi histograms.
  {

    float weff1 = 1.0;
    float zvtx = collision.posZ();

    for (auto const& track1 : tracks) {

      if (!trackSelected(track1)) {
        continue;
      }
      if (!getEfficiencyCorrection(weff1, track1.eta(), track1.pt(), zvtx)) {
        continue;
      }

      registry.fill(HIST("Phi"), RecoDecay::constrainAngle(track1.phi(), 0.0));
      registry.fill(HIST("Eta"), track1.eta());
      registry.fill(HIST("EtaCorrected"), track1.eta(), weff1);
      registry.fill(HIST("pT"), track1.pt());
      registry.fill(HIST("pTCorrected"), track1.pt(), weff1);
    }
  }

  template <CorrelationContainer::CFStep step, typename TFT0s>
  void fillCorrelationsFT0AFT0C(TFT0s const& ft0Col1, TFT0s const& ft0Col2, float posZ, int system, float multiplicity, float eventWeight) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {
    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);

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
        registry.fill(HIST("Trig_hist_FT0A_FT0C"), fSampleIndex, posZ, 0.5, eventWeight * amplA);
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
          registry.fill(HIST("deltaEta_deltaPhi_same_FT0A_FT0C"), deltaPhi, deltaEta, amplA * amplC * eventWeight * triggerWeight);
          sameFt0aFt0c->getPairHist()->Fill(step, fSampleIndex, posZ, 0.5, multiplicity, deltaPhi, deltaEta, amplA * amplC * eventWeight * triggerWeight);
        } else if (system == MixedEvent) {
          registry.fill(HIST("deltaEta_deltaPhi_mixed_FT0A_FT0C"), deltaPhi, deltaEta, amplA * amplC * eventWeight * triggerWeight);
          mixedFt0aFt0c->getPairHist()->Fill(step, fSampleIndex, posZ, 0.5, multiplicity, deltaPhi, deltaEta, amplA * amplC * eventWeight * triggerWeight);
        }
      }
    }
  }

  template <CorrelationContainer::CFStep step, typename TTracks, typename TTracksAssoc>
  void fillCorrelations(TTracks tracks1, TTracksAssoc tracks2, float posZ, float multiplicity, int system, int magneticField) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {

    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);

    float triggerWeight = 1.0f;

    // loop over all tracks
    for (auto const& track1 : tracks1) {

      if (!trackSelected(track1))
        continue;

      if (!getEfficiencyCorrection(triggerWeight, track1.eta(), track1.pt(), posZ))
        continue;

      if (system == SameEvent) {
        registry.fill(HIST("Trig_hist"), fSampleIndex, posZ, track1.pt(), triggerWeight);
      }

      for (auto const& track2 : tracks2) {

        if (!trackSelected(track2))
          continue;

        if (track1.pt() <= track2.pt())
          continue; // skip if the trigger pt is less than the associate pt

        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -PIHalf);
        float deltaEta = track1.eta() - track2.eta();

        if (cfgApplyTwoTrackEfficiency && std::abs(deltaEta) < cfgMergingCut) {

          double dPhiStarHigh = getDPhiStar(track1, track2, cfgRadiusHigh, magneticField);
          double dPhiStarLow = getDPhiStar(track1, track2, cfgRadiusLow, magneticField);

          const double kLimit = 3.0 * cfgMergingCut;

          bool bIsBelow = false;

          if (std::abs(dPhiStarLow) < kLimit || std::abs(dPhiStarHigh) < kLimit || dPhiStarLow * dPhiStarHigh < 0) {
            for (double rad(cfgRadiusLow); rad < cfgRadiusHigh; rad += 0.01) {
              double dPhiStar = getDPhiStar(track1, track2, rad, magneticField);
              if (std::abs(dPhiStar) < kLimit) {
                bIsBelow = true;
                break;
              }
            }
            if (bIsBelow)
              continue;
          }
        }

        // fill the right sparse and histograms
        if (system == SameEvent) {

          sameTPC->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), multiplicity, deltaPhi, deltaEta);
          registry.fill(HIST("deltaEta_deltaPhi_same_TPC"), deltaPhi, deltaEta);

        } else if (system == MixedEvent) {

          mixedTPC->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), multiplicity, deltaPhi, deltaEta);
          registry.fill(HIST("deltaEta_deltaPhi_mixed_TPC"), deltaPhi, deltaEta);
        }
      }
    }
  }

  void processSameTpcFt0a(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
  {
    if (!collision.sel8())
      return;

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), true))
      return;

    if (!collision.has_foundFT0())
      return;

    loadAlignParam(bc.timestamp());
    loadGain(bc);
    loadCorrection(bc.timestamp());
    float eventWeight = 1.0f;

    registry.fill(HIST("zVtx"), collision.posZ());

    registry.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin

    fillYield(collision, tracks);

    const auto& ft0 = collision.foundFT0();
    fillCorrelationsTPCFT0<CorrelationContainer::kCFStepReconstructed>(tracks, ft0, collision.posZ(), SameEvent, kFT0A, tracks.size(), eventWeight);
  }
  PROCESS_SWITCH(CorrFt0Nch, processSameTpcFt0a, "Process same event for TPC-FT0 correlation", false);

  void processMixedTpcFt0a(FilteredCollisions const& collisions, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
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
      if (!collision1.sel8() || !collision2.sel8())
        continue;

      if (cfgUseAdditionalEventCut && !eventSelected(collision1, tracks1.size(), false))
        continue;
      if (cfgUseAdditionalEventCut && !eventSelected(collision2, tracks2.size(), false))
        continue;

      if (!(collision1.has_foundFT0() && collision2.has_foundFT0()))
        continue;

      registry.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin
      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      loadAlignParam(bc.timestamp());
      loadCorrection(bc.timestamp());
      float eventWeight = 1.0f;

      const auto& ft0 = collision2.foundFT0();
      fillCorrelationsTPCFT0<CorrelationContainer::kCFStepReconstructed>(tracks1, ft0, collision1.posZ(), MixedEvent, kFT0A, tracks1.size(), eventWeight);
    }
  }
  PROCESS_SWITCH(CorrFt0Nch, processMixedTpcFt0a, "Process mixed events for TPC-FT0A correlation", false);

  void processSameTpcFt0c(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
  {
    if (!collision.sel8())
      return;
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), true))
      return;
    if (!collision.has_foundFT0())
      return;
    loadAlignParam(bc.timestamp());
    loadGain(bc);
    loadCorrection(bc.timestamp());

    registry.fill(HIST("zVtx"), collision.posZ());

    registry.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin

    const auto& ft0 = collision.foundFT0();

    int multiplicity = 0;
    trackCounter(tracks, multiplicity);

    fillCorrelationsTPCFT0<CorrelationContainer::kCFStepReconstructed>(tracks, ft0, collision.posZ(), SameEvent, kFT0C, multiplicity, 1.0f);
  }
  PROCESS_SWITCH(CorrFt0Nch, processSameTpcFt0c, "Process same event for TPC-FT0C correlation", true);

  void processMixedTpcFt0c(FilteredCollisions const& collisions, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
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
      if (!collision1.sel8() || !collision2.sel8())
        continue;

      if (cfgUseAdditionalEventCut && !eventSelected(collision1, tracks1.size(), false))
        continue;

      if (cfgUseAdditionalEventCut && !eventSelected(collision2, tracks2.size(), false))
        continue;

      if (!(collision1.has_foundFT0() && collision2.has_foundFT0()))
        continue;

      registry.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin
      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      loadAlignParam(bc.timestamp());
      loadCorrection(bc.timestamp());
      float eventWeight = 1.0f;

      const auto& ft0 = collision2.foundFT0();

      int multiplicity = 0;
      trackCounter(tracks1, multiplicity);

      fillCorrelationsTPCFT0<CorrelationContainer::kCFStepReconstructed>(tracks1, ft0, collision1.posZ(), MixedEvent, kFT0C, multiplicity, eventWeight);
    }
  }
  PROCESS_SWITCH(CorrFt0Nch, processMixedTpcFt0c, "Process mixed events for TPC-FT0C correlation", true);

  void processSameFt0aFt0c(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
  {
    if (!collision.sel8())
      return;
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), true))
      return;
    if (!collision.has_foundFT0())
      return;
    loadAlignParam(bc.timestamp());
    loadGain(bc);
    loadCorrection(bc.timestamp());
    float eventWeight = 1.0f;

    // should have the same event to TPC-FT0A/C correlations

    const auto& ft0 = collision.foundFT0();
    int multiplicity = 0;
    trackCounter(tracks, multiplicity);

    fillYield(collision, tracks);

    fillCorrelationsFT0AFT0C<CorrelationContainer::kCFStepReconstructed>(ft0, ft0, collision.posZ(), SameEvent, multiplicity, eventWeight);
  }
  PROCESS_SWITCH(CorrFt0Nch, processSameFt0aFt0c, "Process same event for FT0A-FT0C correlation", false);

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
      if (!collision1.sel8() || !collision2.sel8())
        continue;

      if (cfgUseAdditionalEventCut && !eventSelected(collision1, tracks1.size(), false))
        continue;
      if (cfgUseAdditionalEventCut && !eventSelected(collision2, tracks2.size(), false))
        continue;

      if (!(collision1.has_foundFT0() && collision2.has_foundFT0()))
        continue;

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      loadAlignParam(bc.timestamp());
      loadCorrection(bc.timestamp());
      float eventWeight = 1.0f;

      const auto& ft0Col1 = collision1.foundFT0();
      const auto& ft0Col2 = collision2.foundFT0();

      int multiplicity = 0;

      trackCounter(tracks, multiplicity);

      fillCorrelationsFT0AFT0C<CorrelationContainer::kCFStepReconstructed>(ft0Col1, ft0Col2, collision1.posZ(), MixedEvent, multiplicity, eventWeight);
    }
  }
  PROCESS_SWITCH(CorrFt0Nch, processMixedFt0aFt0c, "Process mixed events for FT0A-FT0C correlation", false);

  void processSameTPC(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::BCsWithTimestamps const&)
  {
    if (!collision.sel8())
      return;

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), true))
      return;

    registry.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin

    int multiplicity = 0;
    trackCounter(tracks, multiplicity);

    fillYield(collision, tracks);

    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(tracks, tracks, collision.posZ(), multiplicity, SameEvent, getMagneticField(bc.timestamp()));
  }
  PROCESS_SWITCH(CorrFt0Nch, processSameTPC, "Process same event for TPC-TPC correlation", false);

  void processMixedTPC(FilteredCollisions const& collisions, FilteredTracks const& tracks, aod::BCsWithTimestamps const&)
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
      if (!collision1.sel8() || !collision2.sel8())
        continue;

      if (cfgUseAdditionalEventCut && !eventSelected(collision1, tracks1.size(), false))
        continue;
      if (cfgUseAdditionalEventCut && !eventSelected(collision2, tracks2.size(), false))
        continue;

      registry.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin

      int multiplicity = 0;
      trackCounter(tracks1, multiplicity);

      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(tracks1, tracks2, collision1.posZ(), multiplicity, MixedEvent, getMagneticField(collision1.bc_as<aod::BCsWithTimestamps>().timestamp()));
    }
  }
  PROCESS_SWITCH(CorrFt0Nch, processMixedTPC, "Process mixed events for TPC-TPC correlation", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CorrFt0Nch>(cfgc),
  };
}
