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

/// \file corrReso.cxx
/// \brief Ultra long range correlation using forward FIT detectors and TPC, with focus on resonances
/// \author Thor Jensen (thor.kjaersgaard.jensen@cern.ch)

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

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
#include <CommonConstants/PhysicsConstants.h>
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

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

static constexpr float LongArrayFloat[3][20] = {{1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}, {2.1, 2.2, 2.3, -2.1, -2.2, -2.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}, {3.1, 3.2, 3.3, -3.1, -3.2, -3.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}};
static constexpr int LongArrayInt[3][20] = {{1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1}, {2, 2, 2, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1}, {3, 3, 3, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1}};

struct CorrReso {
  o2::aod::ITSResponse itsResponse;
  Service<ccdb::BasicCCDBManager> ccdb;

  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, true, "Use additional event cut on mult correlations")
  O2_DEFINE_CONFIGURABLE(cfgZVtxCut, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgUseTransverseMomentum, bool, false, "Use transverse momentum for correlation container")
  O2_DEFINE_CONFIGURABLE(cfgQaCheck, bool, true, "Enable QA histograms for event selection")
  O2_DEFINE_CONFIGURABLE(cfgStrictTrackCounter, bool, false, "Strict track counter for multiplicity correlation cut, counts only tracks that pass all cuts and are used in the correlation")
  O2_DEFINE_CONFIGURABLE(cfgRefpTt, bool, false, "Apply upper pT cut on reference tracks")
  O2_DEFINE_CONFIGURABLE(cfgRefpTMax, float, 3.0f, "maximum pT for reference tracks if cfgRefpTt is true")
  O2_DEFINE_CONFIGURABLE(cfgMinMultForCorrelations, int, 0, "minimum multiplicity for correlations")
  O2_DEFINE_CONFIGURABLE(cfgMaxMultForCorrelations, int, 20, "maximum multiplicity for correlations")
  O2_DEFINE_CONFIGURABLE(cfgRefMultiplicity, bool, false, "Use multiplicity of reference tracks for multiplicity correlation cut instead of Nch")
  Configurable<std::vector<int>> cfgRunRemoveList{"cfgRunRemoveList", std::vector<int>{-1}, "excluded run numbers"};

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

  Configurable<LabeledArray<int>> cfgUseEventCuts{"cfgUseEventCuts", {LongArrayInt[0], 14, 1, {"Filtered Events", "Sel8", "kNoTimeFrameBorder", "kNoITSROFrameBorder", "kNoSameBunchPileup", "kIsGoodZvtxFT0vsPV", "kNoCollInTimeRangeStandard", "kIsGoodITSLayersAll", "kIsGoodITSLayer0123", "kNoCollInRofStandard", "kNoHighMultCollInPrevRof", "Occupancy", "Multcorrelation", "T0AV0ACut"}, {"EvCuts"}}, "Labeled array (int) for various cuts on resonances"};

  O2_DEFINE_CONFIGURABLE(cfgMinMixEventNum, int, 5, "Minimum number of events to mix")
  O2_DEFINE_CONFIGURABLE(cfgSampleSize, double, 10, "Sample size for mixed event")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgEfficiencyNch, std::string, "", "CCDB path to multiplicity dependent efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgCentralityWeight, std::string, "", "CCDB path to centrality weight object")
  O2_DEFINE_CONFIGURABLE(cfgLocalEfficiency, bool, false, "Use local efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgLocalEfficiencyNch, bool, false, "Use local multiplicity dependent efficiency object");
  O2_DEFINE_CONFIGURABLE(cfgUseEventWeights, bool, false, "Use event weights for mixed event")
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
    O2_DEFINE_CONFIGURABLE(cfgUseOnlyTPC, bool, true, "Use only TPC PID for daughter selection")
    O2_DEFINE_CONFIGURABLE(cfgUseAntiLambda, bool, true, "Use AntiLambda candidates for analysis")
    O2_DEFINE_CONFIGURABLE(cfgPIDUseRejection, bool, true, "True: use exclusion exclusion criteria for PID determination, false: don't use exclusion")
    O2_DEFINE_CONFIGURABLE(cfgTpcCut, float, 3.0f, "TPC N-sigma cut for pions, kaons, protons")
    O2_DEFINE_CONFIGURABLE(cfgPIDParticle, int, 0, "4 = kshort, 5 = lambda, 6 = phi, 0 for no PID")
    O2_DEFINE_CONFIGURABLE(cfgUseItsPID, bool, true, "Use ITS PID for particle identification")
    O2_DEFINE_CONFIGURABLE(cfgTofPtCut, float, 0.4f, "Minimum pt to use TOF N-sigma")
    Configurable<LabeledArray<float>> nSigmas{"nSigmas", {LongArrayFloat[0], 6, 3, {"UpCut_pi", "UpCut_ka", "UpCut_pr", "LowCut_pi", "LowCut_ka", "LowCut_pr"}, {"TPC", "TOF", "ITS"}}, "Labeled array for n-sigma values for TPC, TOF, ITS for pions, kaons, protons (positive and negative)"};
    Configurable<LabeledArray<float>> cfgResoCuts{"cfgResoCuts", {LongArrayFloat[0], 12, 3, {"cos_PAs", "massMin", "massMax", "PosTrackPt", "NegTrackPt", "DCAPosToPVMin", "DCANegToPVMin", "Lifetime", "RadiusMin", "RadiusMax", "Rapidity", "ArmPodMinVal"}, {"K0", "Lambda", "Phi"}}, "Labeled array (float) for various cuts on resonances"};
    Configurable<LabeledArray<int>> cfgResoSwitches{"cfgResoSwitches", {LongArrayInt[0], 6, 3, {"UseCosPA", "NMassBins", "DCABetDaug", "UseProperLifetime", "UseV0Radius", "UseArmPodCut"}, {"K0", "Lambda", "Phi"}}, "Labeled array (int) for various cuts on resonances"};
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
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt trigger axis for histograms"};
  ConfigurableAxis axisVtxMix{"axisVtxMix", {VARIABLE_WIDTH, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, "vertex axis for mixed event histograms"};
  ConfigurableAxis axisMultMix{"axisMultMix", {VARIABLE_WIDTH, 0, 10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260}, "multiplicity / centrality axis for mixed event histograms"};
  ConfigurableAxis axisSample{"axisSample", {cfgSampleSize, 0, cfgSampleSize}, "sample axis for histograms"};
  ConfigurableAxis axisNch{"axisNch", {VARIABLE_WIDTH, 0, 10, 50, 70, 100}, "multiplicity axis for correlation container"};

  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {80, -5, 5}, "nsigmaTPC axis"};
  ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {80, -5, 5}, "nsigmaTOF axis"};

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
  Filter trackFilter = (nabs(aod::track::eta) < cfgTrackCuts.cfgEtaCut) && (cfgTrackCuts.cfgPtCutMin < aod::track::pt) && (cfgTrackCuts.cfgPtCutMax > aod::track::pt) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t)true));

  using FilteredCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSel, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::Mults>>;
  using FilteredTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;
  using V0TrackCandidate = aod::V0Datas;

  // FT0 geometry
  o2::ft0::Geometry ft0Det;
  static constexpr uint64_t Ft0IndexA = 96;
  std::vector<o2::detectors::AlignParam>* offsetFT0;
  std::vector<float> cstFT0RelGain{};

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

  HistogramRegistry registry{"registry"};

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
  enum ResoArrayIndex {
    iK0 = 0,
    iLambda = 1,
    iPhi = 2,
    NResoParticles = 3
  };
  enum ResoParticleCuts {
    kCosPA = 0,
    kMassMin,
    kMassMax,
    kPosTrackPt,
    kNegTrackPt,
    kDCAPosToPVMin,
    kDCANegToPVMin,
    kLifeTime,
    kRadiusMin,
    kRadiusMax,
    kRapidity,
    kArmPodMinVal,
    kNParticleCuts
  };
  enum ResoParticleSwitches {
    kUseCosPA = 0,
    kMassBins,
    kDCABetDaug,
    kUseProperLifetime,
    kUseV0Radius,
    kUseArmPodCut,
    kNParticleSwitches
  };
  enum DetectorType {
    kTPC = 0,
    kTOF,
    kITS
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
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    LOGF(info, "Starting init");

    if (cfgPIDConfigs.cfgPIDParticle < kK0) {
      LOGF(fatal, "The input number does not correspond to any particle");
      return;
    }

    // Creating mass axis depending on particle - 4 = kshort, 5 = lambda, 6 = phi
    AxisSpec axisInvMass = {10, 0, 1, "mass"};
    if (cfgPIDConfigs.cfgPIDParticle == kK0)
      axisInvMass = {cfgPIDConfigs.cfgResoSwitches->getData()[kMassBins][iK0], cfgPIDConfigs.cfgResoCuts->getData()[kMassMin][iK0], cfgPIDConfigs.cfgResoCuts->getData()[kMassMax][iK0], "M_{#pi^{+}#pi^{-}} (GeV/c^{2})"};
    if (cfgPIDConfigs.cfgPIDParticle == kLambda)
      axisInvMass = {cfgPIDConfigs.cfgResoSwitches->getData()[kMassBins][iLambda], cfgPIDConfigs.cfgResoCuts->getData()[kMassMin][iLambda], cfgPIDConfigs.cfgResoCuts->getData()[kMassMax][iLambda], "M_{p#pi^{-}} (GeV/c^{2})"};
    if (cfgPIDConfigs.cfgPIDParticle == kPhi)
      axisInvMass = {cfgPIDConfigs.cfgResoSwitches->getData()[kMassBins][iPhi], cfgPIDConfigs.cfgResoCuts->getData()[kMassMin][iPhi], cfgPIDConfigs.cfgResoCuts->getData()[kMassMax][iPhi], "M_{K^{+}K^{-}} (GeV/c^{2})"};

    if (doprocessSameTpcFt0a || doprocessSameTpcFt0c) {
      if (cfgPIDConfigs.cfgPIDParticle == kK0) { // For K0
        registry.add("PiPlusTPC_K0", "", {HistType::kTH2D, {{axisPtFiner, axisNsigmaTPC}}});
        registry.add("PiMinusTPC_K0", "", {HistType::kTH2D, {{axisPtFiner, axisNsigmaTPC}}});
        registry.add("PiPlusTOF_K0", "", {HistType::kTH2D, {{axisPtFiner, axisNsigmaTOF}}});
        registry.add("PiMinusTOF_K0", "", {HistType::kTH2D, {{axisPtFiner, axisNsigmaTOF}}});
        registry.add("hK0Phi", "", {HistType::kTH1D, {axisPhi}});
        registry.add("hK0Eta", "", {HistType::kTH1D, {axisEta}});

        registry.add("hK0Count", "Number of K0;; Count", {HistType::kTH1D, {{11, 0, 11}}});
        registry.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(1, "K0 candidates");
        registry.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(2, "Daughter pt");
        registry.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(3, "Mass cut");
        registry.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(4, "Rapidity cut");
        registry.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(5, "DCA to PV");
        registry.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(6, "DCA between daughters");
        registry.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(7, "V0radius");
        registry.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(8, "CosPA");
        registry.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(9, "Proper lifetime");
        registry.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(10, "ArmenterosPod");
        registry.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(11, "Daughter track selection");
      }
      if (cfgPIDConfigs.cfgPIDParticle == kLambda) { // For Lambda
        registry.add("PrPlusTPC_La", "", {HistType::kTH2D, {{axisPtFiner, axisNsigmaTPC}}});
        registry.add("PiMinusTPC_La", "", {HistType::kTH2D, {{axisPtFiner, axisNsigmaTPC}}});
        registry.add("PrPlusTOF_La", "", {HistType::kTH2D, {{axisPtFiner, axisNsigmaTOF}}});
        registry.add("PiMinusTOF_La", "", {HistType::kTH2D, {{axisPtFiner, axisNsigmaTOF}}});

        if (cfgPIDConfigs.cfgUseAntiLambda) {
          registry.add("PrMinusTPC_Al", "", {HistType::kTH2D, {{axisPtFiner, axisNsigmaTPC}}});
          registry.add("PiPlusTPC_Al", "", {HistType::kTH2D, {{axisPtFiner, axisNsigmaTPC}}});
          registry.add("PrMinusTOF_Al", "", {HistType::kTH2D, {{axisPtFiner, axisNsigmaTOF}}});
          registry.add("PiPlusTOF_Al", "", {HistType::kTH2D, {{axisPtFiner, axisNsigmaTOF}}});
        }

        registry.add("hLambdaPhi", "", {HistType::kTH1D, {axisPhi}});
        registry.add("hLambdaEta", "", {HistType::kTH1D, {axisEta}});

        registry.add("hLambdaCount", "Number of Lambda;; Count", {HistType::kTH1D, {{10, 0, 10}}});
        registry.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(1, "Lambda candidates");
        registry.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(2, "Daughter pt");
        registry.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(3, "Mass cut");
        registry.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(4, "Rapidity cut");
        registry.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(5, "DCA to PV");
        registry.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(6, "DCA between daughters");
        registry.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(7, "V0radius");
        registry.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(8, "CosPA");
        registry.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(9, "Proper lifetime");
        registry.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(10, "Daughter track selection");
      }
    }
    if (doprocessSameFt0aFt0c || doprocessSameTpcFt0a || doprocessSameTpcFt0c) {
      registry.add("hEventCountRct", "Number of Event;; Count", {HistType::kTH1D, {{2, 0, 2}}});
      registry.get<TH1>(HIST("hEventCountRct"))->GetXaxis()->SetBinLabel(1, "rct fail");
      registry.get<TH1>(HIST("hEventCountRct"))->GetXaxis()->SetBinLabel(2, "rct pass");
      registry.add("hEventCount", "Number of Event;; Count", {HistType::kTH1D, {{14, -0.5, 13.5}}});
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kFilteredEvents + 1, "Filtered events");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kAfterSel8 + 1, "After sel8");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoTimeFrameBorder + 1, "kNoTimeFrameBorder");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoITSROFrameBorder + 1, "kNoITSROFrameBorder");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoSameBunchPileup + 1, "kNoSameBunchPileup");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseGoodZvtxFT0vsPV + 1, "kIsGoodZvtxFT0vsPV");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoCollInTimeRangeStandard + 1, "kNoCollInTimeRangeStandard");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseGoodITSLayersAll + 1, "kIsGoodITSLayersAll");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseGoodITSLayer0123 + 1, "kIsGoodITSLayer0123");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoCollInRofStandard + 1, "kNoCollInRofStandard");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoHighMultCollInPrevRof + 1, "kNoHighMultCollInPrevRof");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseOccupancy + 1, "Occupancy Cut");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseMultCorrCut + 1, "MultCorrelation Cut");
      registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseT0AV0ACut + 1, "T0AV0A cut");
    }

    if ((doprocessSameFt0aFt0c || doprocessSameTpcFt0a || doprocessSameTpcFt0c) && cfgQaCheck) {
      registry.add("hPassedEventSelection", "Number of Event;; Count", {HistType::kTH1D, {{12, -0.5, 11.5}}});
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kFilteredEvents + 1, "Filtered events");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kAfterSel8 + 1, "After sel8");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kUseNoTimeFrameBorder + 1, "kNoTimeFrameBorder");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kUseNoITSROFrameBorder + 1, "kNoITSROFrameBorder");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kUseNoSameBunchPileup + 1, "kNoSameBunchPileup");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kUseGoodZvtxFT0vsPV + 1, "kIsGoodZvtxFT0vsPV");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kUseNoCollInTimeRangeStandard + 1, "kNoCollInTimeRangeStandard");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kUseGoodITSLayersAll + 1, "kIsGoodITSLayersAll");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kUseGoodITSLayer0123 + 1, "kIsGoodITSLayer0123");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kUseNoCollInRofStandard + 1, "kNoCollInRofStandard");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kUseNoHighMultCollInPrevRof + 1, "kNoHighMultCollInPrevRof");
      registry.get<TH1>(HIST("hPassedEventSelection"))->GetXaxis()->SetBinLabel(kUseOccupancy + 1, "Occupancy Cut");
    }

    // Multiplicity correlation cuts
    if (cfgUseEventCuts->getData()[kUseMultCorrCut][kEvCut1]) {
      cfgFuncParas.multT0CCutPars = cfgFuncParas.cfgMultT0CCutPars;
      cfgFuncParas.multPVT0CCutPars = cfgFuncParas.cfgMultPVT0CCutPars;
      cfgFuncParas.multGlobalPVCutPars = cfgFuncParas.cfgMultGlobalPVCutPars;
      cfgFuncParas.multMultV0ACutPars = cfgFuncParas.cfgMultMultV0ACutPars;
      cfgFuncParas.fMultPVT0CCutLow = new TF1("fMultPVT0CCutLow", cfgFuncParas.cfgMultCentLowCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultPVT0CCutLow->SetParameters(&(cfgFuncParas.multPVT0CCutPars[0]));
      cfgFuncParas.fMultPVT0CCutHigh = new TF1("fMultPVT0CCutHigh", cfgFuncParas.cfgMultCentHighCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultPVT0CCutHigh->SetParameters(&(cfgFuncParas.multPVT0CCutPars[0]));

      cfgFuncParas.fMultT0CCutLow = new TF1("fMultT0CCutLow", cfgFuncParas.cfgMultCentLowCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultT0CCutLow->SetParameters(&(cfgFuncParas.multT0CCutPars[0]));
      cfgFuncParas.fMultT0CCutHigh = new TF1("fMultT0CCutHigh", cfgFuncParas.cfgMultCentHighCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultT0CCutHigh->SetParameters(&(cfgFuncParas.multT0CCutPars[0]));

      cfgFuncParas.fMultGlobalPVCutLow = new TF1("fMultGlobalPVCutLow", cfgFuncParas.cfgMultMultPVLowCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultGlobalPVCutLow->SetParameters(&(cfgFuncParas.multGlobalPVCutPars[0]));
      cfgFuncParas.fMultGlobalPVCutHigh = new TF1("fMultGlobalPVCutHigh", cfgFuncParas.cfgMultMultPVHighCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultGlobalPVCutHigh->SetParameters(&(cfgFuncParas.multGlobalPVCutPars[0]));

      cfgFuncParas.fMultMultV0ACutLow = new TF1("fMultMultV0ACutLow", cfgFuncParas.cfgMultMultV0ALowCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultMultV0ACutLow->SetParameters(&(cfgFuncParas.multMultV0ACutPars[0]));
      cfgFuncParas.fMultMultV0ACutHigh = new TF1("fMultMultV0ACutHigh", cfgFuncParas.cfgMultMultV0AHighCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultMultV0ACutHigh->SetParameters(&(cfgFuncParas.multMultV0ACutPars[0]));
    }
    if (cfgUseEventCuts->getData()[kUseT0AV0ACut][kEvCut1]) {
      cfgFuncParas.fT0AV0AMean = new TF1("fT0AV0AMean", "[0]+[1]*x", 0, 200000);
      cfgFuncParas.fT0AV0AMean->SetParameters(-1601.0581, 9.417652e-01);
      cfgFuncParas.fT0AV0ASigma = new TF1("fT0AV0ASigma", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 200000);
      cfgFuncParas.fT0AV0ASigma->SetParameters(463.4144, 6.796509e-02, -9.097136e-07, 7.971088e-12, -2.600581e-17);
    }

    if (doprocessSameTpcFt0a || doprocessSameTpcFt0c || doprocessSameFt0aFt0c) {
      registry.add("Phi", "Phi", {HistType::kTH1D, {axisPhi}});
      registry.add("Eta", "Eta", {HistType::kTH1D, {axisEta}});
      registry.add("EtaCorrected", "EtaCorrected", {HistType::kTH1D, {axisEta}});
      registry.add("pTFiner", "pTFiner", {HistType::kTH1D, {axisPtFiner}});
      registry.add("pTFinerCorrected", "pTFinerCorrected", {HistType::kTH1D, {axisPtFiner}});
      registry.add("Nch", "N_{ch}", {HistType::kTH1D, {axisMult}});
      registry.add("zVtx", "zVtx", {HistType::kTH1D, {axisVertex}});

      registry.add("FT0Amp", "", {HistType::kTH2F, {axisChID, axisFit}});
      registry.add("FT0AmpCorrect", "", {HistType::kTH2F, {axisChID, axisFit}});
    } // end of single particle distribution histograms

    if (cfgQaCheck) {
      registry.add("Nch_corrected", "N_{ch} corrected", {HistType::kTH1D, {axisMult}});
    }

    if (doprocessSameTpcFt0a) {
      registry.add("deltaEta_deltaPhi_same_TPC_FT0A", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0a}}); // check to see the delta eta and delta phi distribution
      registry.add("deltaEta_deltaPhi_mixed_TPC_FT0A", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0a}});
      registry.add("Assoc_amp_same_TPC_FT0A", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
      registry.add("Assoc_amp_mixed_TPC_FT0A", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
      registry.add("Trig_hist_TPC_FT0A", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger, axisInvMass}}});
    }
    if (doprocessSameTpcFt0c) {
      registry.add("deltaEta_deltaPhi_same_TPC_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0c}}); // check to see the delta eta and delta phi distribution
      registry.add("deltaEta_deltaPhi_mixed_TPC_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0c}});
      registry.add("Assoc_amp_same_TPC_FT0C", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
      registry.add("Assoc_amp_mixed_TPC_FT0C", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
      registry.add("Trig_hist_TPC_FT0C", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger, axisInvMass}}});
    }
    if (doprocessSameFt0aFt0c) {
      registry.add("deltaEta_deltaPhi_same_FT0A_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaFt0aFt0c}}); // check to see the delta eta and delta phi distribution
      registry.add("deltaEta_deltaPhi_mixed_FT0A_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaFt0aFt0c}});
      registry.add("Trig_hist_FT0A_FT0C", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
    }

    registry.add("eventcount", "bin", {HistType::kTH1F, {{4, 0, 4, "bin"}}}); // histogram to see how many events are in the same and mixed event

    LOGF(info, "Initializing correlation container");

    // Initialize Nch-related histograms and containers

    std::vector<AxisSpec> corrAxisTpcFt0a = {{axisSample, "Sample"},
                                             {axisVertex, "z-vtx (cm)"},
                                             {axisPtTrigger, "p_{T} (GeV/c)"},
                                             {axisInvMass},
                                             {axisDeltaPhi, "#Delta#varphi (rad)"},
                                             {axisDeltaEtaTpcFt0a, "#Delta#eta"}};

    std::vector<AxisSpec> corrAxisTpcFt0c = {{axisSample, "Sample"},
                                             {axisVertex, "z-vtx (cm)"},
                                             {axisPtTrigger, "p_{T} (GeV/c)"},
                                             {axisInvMass},
                                             {axisDeltaPhi, "#Delta#varphi (rad)"},
                                             {axisDeltaEtaTpcFt0c, "#Delta#eta"}};

    std::vector<AxisSpec> corrAxisFt0aFt0c = {{axisSample, "Sample"},
                                              {axisVertex, "z-vtx (cm)"},
                                              {axisPtTrigger, "p_{T} (GeV/c)"},
                                              {axisNch, "N_{ch}"},
                                              {axisDeltaPhi, "#Delta#varphi (rad)"},
                                              {axisDeltaEtaFt0aFt0c, "#Delta#eta"}};

    std::vector<AxisSpec> effAxis = {
      {axisEtaEfficiency, "#eta"},
      {axisPtEfficiency, "p_{T} (GeV/c)"},
      {axisVertexEfficiency, "z-vtx (cm)"},
    };
    std::vector<AxisSpec> userAxis;

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

    LOGF(info, "End of init");
  }

  template <typename TCollision>
  float getCentrality(TCollision const& collision)
  {
    float cent;
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
      if (fillCounter)
        registry.fill(HIST("hEventCountRct"), 0.5);

      return 0;
    }
    if (fillCounter)
      registry.fill(HIST("hEventCountRct"), 1.5);

    return 1;
  }

  template <typename TCollision>
  bool eventSelected(TCollision const& collision, const int mult, const bool fillCounter)
  {
    if (cfgUseEventCuts->getData()[kUseNoTimeFrameBorder][kEvCut1] && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return 0;
    }
    if (fillCounter && cfgUseEventCuts->getData()[kUseNoTimeFrameBorder][kEvCut1])
      registry.fill(HIST("hEventCount"), kUseNoTimeFrameBorder);

    if (cfgUseEventCuts->getData()[kUseNoITSROFrameBorder][kEvCut1] && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return 0;
    }
    if (fillCounter && cfgUseEventCuts->getData()[kUseNoITSROFrameBorder][kEvCut1])
      registry.fill(HIST("hEventCount"), kUseNoITSROFrameBorder);

    if (cfgUseEventCuts->getData()[kUseNoSameBunchPileup][kEvCut1] && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return 0;
    }
    if (fillCounter && cfgUseEventCuts->getData()[kUseNoSameBunchPileup][kEvCut1])
      registry.fill(HIST("hEventCount"), kUseNoSameBunchPileup);

    if (cfgUseEventCuts->getData()[kUseGoodZvtxFT0vsPV][kEvCut1] && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return 0;
    }
    if (fillCounter && cfgUseEventCuts->getData()[kUseGoodZvtxFT0vsPV][kEvCut1])
      registry.fill(HIST("hEventCount"), kUseGoodZvtxFT0vsPV);

    if (cfgUseEventCuts->getData()[kUseNoCollInTimeRangeStandard][kEvCut1] && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      // no collisions in specified time range
      return 0;
    }

    if (fillCounter && cfgUseEventCuts->getData()[kUseNoCollInTimeRangeStandard][kEvCut1])
      registry.fill(HIST("hEventCount"), kUseNoCollInTimeRangeStandard);

    if (cfgUseEventCuts->getData()[kUseGoodITSLayersAll][kEvCut1] && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      // from Jan 9 2025 AOT meeting
      // cut time intervals with dead ITS staves
      return 0;
    }

    if (fillCounter && cfgUseEventCuts->getData()[kUseGoodITSLayersAll][kEvCut1])
      registry.fill(HIST("hEventCount"), kUseGoodITSLayersAll);

    if (cfgUseEventCuts->getData()[kUseGoodITSLayer0123][kEvCut1] && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
      return 0;
    }
    if (fillCounter && cfgUseEventCuts->getData()[kUseGoodITSLayer0123][kEvCut1])
      registry.fill(HIST("hEventCount"), kUseGoodITSLayer0123);

    if (cfgUseEventCuts->getData()[kUseNoCollInRofStandard][kEvCut1] && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      // no other collisions in this Readout Frame with per-collision multiplicity above threshold
      return 0;
    }

    if (fillCounter && cfgUseEventCuts->getData()[kUseNoCollInRofStandard][kEvCut1])
      registry.fill(HIST("hEventCount"), kUseNoCollInRofStandard);

    if (cfgUseEventCuts->getData()[kUseNoHighMultCollInPrevRof][kEvCut1] && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      // veto an event if FT0C amplitude in previous ITS ROF is above threshold
      return 0;
    }
    if (fillCounter && cfgUseEventCuts->getData()[kUseNoHighMultCollInPrevRof][kEvCut1])
      registry.fill(HIST("hEventCount"), kUseNoHighMultCollInPrevRof);

    auto multNTracksPV = collision.multNTracksPV();
    auto occupancy = collision.trackOccupancyInTimeRange();

    if (cfgUseEventCuts->getData()[kUseOccupancy][kEvCut1] && (occupancy < cfgEventSelection.cfgCutOccupancyLow || occupancy > cfgEventSelection.cfgCutOccupancyHigh)) {
      return 0;
    }
    if (fillCounter && cfgUseEventCuts->getData()[kUseOccupancy][kEvCut1])
      registry.fill(HIST("hEventCount"), kUseOccupancy);

    if (cfgUseEventCuts->getData()[kUseMultCorrCut][kEvCut1]) {
      float cent = getCentrality(collision);
      if (cfgFuncParas.cfgMultPVT0CCutEnabled) {
        if (multNTracksPV < cfgFuncParas.fMultPVT0CCutLow->Eval(cent))
          return 0;
        if (multNTracksPV > cfgFuncParas.fMultPVT0CCutHigh->Eval(cent))
          return 0;
      }
      if (cfgFuncParas.cfgMultT0CCutEnabled) {
        if (mult < cfgFuncParas.fMultT0CCutLow->Eval(cent))
          return 0;
        if (mult > cfgFuncParas.fMultT0CCutHigh->Eval(cent))
          return 0;
      }
      if (cfgFuncParas.cfgMultGlobalPVCutEnabled) {
        if (mult < cfgFuncParas.fMultGlobalPVCutLow->Eval(multNTracksPV))
          return 0;
        if (mult > cfgFuncParas.fMultGlobalPVCutHigh->Eval(multNTracksPV))
          return 0;
      }
      if (cfgFuncParas.cfgMultMultV0ACutEnabled) {
        if (collision.multFV0A() < cfgFuncParas.fMultMultV0ACutLow->Eval(mult))
          return 0;
        if (collision.multFV0A() > cfgFuncParas.fMultMultV0ACutHigh->Eval(mult))
          return 0;
      }
    }

    if (fillCounter && cfgUseEventCuts->getData()[kUseMultCorrCut][kEvCut1])
      registry.fill(HIST("hEventCount"), kUseMultCorrCut);

    // V0A T0A 5 sigma cut
    if (cfgUseEventCuts->getData()[kUseT0AV0ACut][kEvCut1] && (std::fabs(collision.multFV0A() - cfgFuncParas.fT0AV0AMean->Eval(collision.multFT0A())) > cfgFuncParas.cfgV0AT0Acut * cfgFuncParas.fT0AV0ASigma->Eval(collision.multFT0A())))
      return 0;
    if (fillCounter && cfgUseEventCuts->getData()[kUseT0AV0ACut][kEvCut1])
      registry.fill(HIST("hEventCount"), kUseT0AV0ACut);

    return 1;
  }

  template <typename TCollision>
  void eventSelectedIndividually(TCollision const& collision)
  {
    registry.fill(HIST("hPassedEventSelection"), kFilteredEvents);

    if (collision.sel8()) {
      registry.fill(HIST("hPassedEventSelection"), kAfterSel8);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      registry.fill(HIST("hPassedEventSelection"), kUseNoTimeFrameBorder);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      registry.fill(HIST("hPassedEventSelection"), kUseNoITSROFrameBorder);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      registry.fill(HIST("hPassedEventSelection"), kUseNoSameBunchPileup);
    }

    if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      registry.fill(HIST("hPassedEventSelection"), kUseGoodZvtxFT0vsPV);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      registry.fill(HIST("hPassedEventSelection"), kUseNoCollInTimeRangeStandard);
    }

    if (collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      registry.fill(HIST("hPassedEventSelection"), kUseGoodITSLayersAll);
    }

    if (collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
      registry.fill(HIST("hPassedEventSelection"), kUseGoodITSLayer0123);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      registry.fill(HIST("hPassedEventSelection"), kUseNoCollInRofStandard);
    }

    if (collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      registry.fill(HIST("hPassedEventSelection"), kUseNoHighMultCollInPrevRof);
    }

    auto occupancy = collision.trackOccupancyInTimeRange();
    if (cfgUseEventCuts->getData()[kUseOccupancy][kEvCut1] && (occupancy < cfgEventSelection.cfgCutOccupancyLow || occupancy > cfgEventSelection.cfgCutOccupancyHigh)) {
      registry.fill(HIST("hPassedEventSelection"), kUseOccupancy);
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

  template <typename TTrack>
  int getNsigmaPID(TTrack track)
  {
    // Computing Nsigma arrays for pion, kaon, and protons
    std::array<float, 3> nSigmaTPC = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::array<float, 3> nSigmaTOF = {track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};
    std::array<float, 3> nSigmaITS = {itsResponse.nSigmaITS<o2::track::PID::Pion>(track), itsResponse.nSigmaITS<o2::track::PID::Kaon>(track), itsResponse.nSigmaITS<o2::track::PID::Proton>(track)};
    int pid = -1; // -1 = not identified, 1 = pion, 2 = kaon, 3 = proton

    std::array<float, 3> nSigmaToUse = cfgPIDConfigs.cfgUseItsPID ? nSigmaITS : nSigmaTPC; // Choose which nSigma to use: TPC or ITS
    int kIndexDetector = cfgPIDConfigs.cfgUseItsPID ? kITS : kTPC;                         // Choose which nSigma to use: TPC or ITS

    bool isPion, isKaon, isProton;
    bool isDetectedPion = nSigmaToUse[iPionUp] < cfgPIDConfigs.nSigmas->getData()[iPionUp][kIndexDetector] && nSigmaToUse[iPionUp] > cfgPIDConfigs.nSigmas->getData()[iPionLow][kIndexDetector];
    bool isDetectedKaon = nSigmaToUse[iKaonUp] < cfgPIDConfigs.nSigmas->getData()[iKaonUp][kIndexDetector] && nSigmaToUse[iKaonUp] > cfgPIDConfigs.nSigmas->getData()[iKaonLow][kIndexDetector];
    bool isDetectedProton = nSigmaToUse[iProtonUp] < cfgPIDConfigs.nSigmas->getData()[iProtonUp][kIndexDetector] && nSigmaToUse[iProtonUp] > cfgPIDConfigs.nSigmas->getData()[iProtonLow][kIndexDetector];

    bool isTofPion = nSigmaTOF[iPionUp] < cfgPIDConfigs.nSigmas->getData()[iPionUp][kTOF] && nSigmaTOF[iPionUp] > cfgPIDConfigs.nSigmas->getData()[iPionLow][kTOF];
    bool isTofKaon = nSigmaTOF[iKaonUp] < cfgPIDConfigs.nSigmas->getData()[iKaonUp][kTOF] && nSigmaTOF[iKaonUp] > cfgPIDConfigs.nSigmas->getData()[iKaonLow][kTOF];
    bool isTofProton = nSigmaTOF[iProtonUp] < cfgPIDConfigs.nSigmas->getData()[iProtonUp][kTOF] && nSigmaTOF[iProtonUp] > cfgPIDConfigs.nSigmas->getData()[iProtonLow][kTOF];

    if (track.pt() > cfgPIDConfigs.cfgTofPtCut && !track.hasTOF()) {
      return -1;
    } else if (track.pt() > cfgPIDConfigs.cfgTofPtCut && track.hasTOF()) {
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

  template <typename TTrack>
  bool selectionV0Daughter(TTrack const& track, int pid)
  {
    if (!(track.itsNCls() > cfgTrackCuts.cfgCutITSclu))
      return false;
    if (!track.hasTPC())
      return false;
    if (!(track.tpcNClsFound() >= cfgTrackCuts.cfgCutTPCclu))
      return false;
    if (!(track.tpcNClsCrossedRows() >= cfgTrackCuts.cfgCutTPCCrossedRows))
      return false;
    if (!(track.dcaZ() < cfgTrackCuts.cfgCutDCAz))
      return false;

    if (cfgPIDConfigs.cfgUseOnlyTPC) {
      if (pid == kPions && std::abs(track.tpcNSigmaPi()) > cfgPIDConfigs.cfgTpcCut)
        return false;
      if (pid == kKaons && std::abs(track.tpcNSigmaKa()) > cfgPIDConfigs.cfgTpcCut)
        return false;
      if (pid == kProtons && std::abs(track.tpcNSigmaPr()) > cfgPIDConfigs.cfgTpcCut)
        return false;
    } else {
      int partIndex = getNsigmaPID(track);
      int pidIndex = partIndex; // 1 = pion, 2 = kaon, 3 = proton
      if (pidIndex != pid)
        return false;
    }

    return true;
  }

  void loadCorrection(uint64_t timestamp)
  {
    if (correctionsLoaded) {
      return;
    }
    if (cfgEfficiency.value.empty() == false) {
      if (cfgLocalEfficiency) {
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
    if (cfgEfficiencyNch.value.empty() == false) {
      if (cfgLocalEfficiencyNch) {
        TFile* fEfficiencyTrigger = TFile::Open(cfgEfficiencyNch.value.c_str(), "READ");
        mEfficiencyNch = reinterpret_cast<TH1D*>(fEfficiencyTrigger->Get("ccdb_object"));

      } else {
        mEfficiencyNch = ccdb->getForTimeStamp<TH1D>(cfgEfficiencyNch, timestamp);
      }
      if (!mEfficiencyNch) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiencyNch.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEfficiencyNch.value.c_str(), (void*)mEfficiencyNch);
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

  bool getEfficiencyCorrectionNch(float& weightNch, float pt)
  {
    float effNch = 1.;
    if (mEfficiencyNch) {

      int ptBin = mEfficiencyNch->FindBin(pt);
      effNch = mEfficiencyNch->GetBinContent(ptBin);

    } else {
      effNch = 1.0;
    }
    if (effNch == 0)
      return false;
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
    if (eff == 0)
      return false;
    weight = 1. / eff;
    return true;
  }

  template <typename TTracks>
  void trackCounter(TTracks tracks, double& multiplicity) // function to count the number of tracks in the event and fill the histogram
  {
    double nTracksCorrected = 0;
    float weightNch = 1.0f;
    for (auto const& track : tracks) {

      if (cfgRefMultiplicity) {
        if (track.pt() > cfgRefpTMax)
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
      registry.fill(HIST("pTFiner"), track1.pt());
      registry.fill(HIST("pTFinerCorrected"), track1.pt(), weff1);
    }
  }

  template <typename V0>
  bool isSelectedK0(V0 const& candidate, float posZ, float posY, float posX)
  {
    double mk0 = candidate.mK0Short();

    // separate the positive and negative V0 daughters
    auto postrack = candidate.template posTrack_as<FilteredTracks>();
    auto negtrack = candidate.template negTrack_as<FilteredTracks>();

    registry.fill(HIST("hK0Count"), 0.5);
    if (postrack.pt() < cfgPIDConfigs.cfgResoCuts->getData()[kPosTrackPt][iK0] || negtrack.pt() < cfgPIDConfigs.cfgResoCuts->getData()[kNegTrackPt][iK0])
      return false;
    registry.fill(HIST("hK0Count"), 1.5);
    if (mk0 < cfgPIDConfigs.cfgResoCuts->getData()[kMassMin][iK0] && mk0 > cfgPIDConfigs.cfgResoCuts->getData()[kMassMax][iK0])
      return false;
    registry.fill(HIST("hK0Count"), 2.5);
    // Rapidity correction
    if (candidate.yK0Short() > cfgPIDConfigs.cfgResoCuts->getData()[kRapidity][iK0])
      return false;
    registry.fill(HIST("hK0Count"), 3.5);
    // DCA cuts for K0short
    if (std::abs(candidate.dcapostopv()) < cfgPIDConfigs.cfgResoCuts->getData()[kDCAPosToPVMin][iK0] || std::abs(candidate.dcanegtopv()) < cfgPIDConfigs.cfgResoCuts->getData()[kDCANegToPVMin][iK0])
      return false;
    registry.fill(HIST("hK0Count"), 4.5);
    if (std::abs(candidate.dcaV0daughters()) > cfgPIDConfigs.cfgResoSwitches->getData()[kDCABetDaug][iK0])
      return false;
    registry.fill(HIST("hK0Count"), 5.5);
    // v0 radius cuts
    if (cfgPIDConfigs.cfgResoSwitches->getData()[kUseV0Radius][iK0] && (candidate.v0radius() < cfgPIDConfigs.cfgResoCuts->getData()[kRadiusMin][iK0] || candidate.v0radius() > cfgPIDConfigs.cfgResoCuts->getData()[kRadiusMax][iK0]))
      return false;
    registry.fill(HIST("hK0Count"), 6.5);
    // cosine pointing angle cuts
    if (candidate.v0cosPA() < cfgPIDConfigs.cfgResoCuts->getData()[kCosPA][iK0])
      return false;
    registry.fill(HIST("hK0Count"), 7.5);
    // Proper lifetime
    float cTauK0 = candidate.distovertotmom(posX, posY, posZ) * massK0Short;
    if (cfgPIDConfigs.cfgResoSwitches->getData()[kUseProperLifetime][iK0] && std::abs(cTauK0) > cfgPIDConfigs.cfgResoCuts->getData()[kLifeTime][iK0])
      return false;
    registry.fill(HIST("hK0Count"), 8.5);
    // ArmenterosPodolanskiCut
    if (cfgPIDConfigs.cfgResoSwitches->getData()[kUseArmPodCut][iK0] && (candidate.qtarm() / std::abs(candidate.alpha())) < cfgPIDConfigs.cfgResoCuts->getData()[kArmPodMinVal][iK0])
      return false;
    registry.fill(HIST("hK0Count"), 9.5);
    // Selection on V0 daughters
    if (!selectionV0Daughter(postrack, kPions) || !selectionV0Daughter(negtrack, kPions))
      return false;
    registry.fill(HIST("hK0Count"), 10.5);

    registry.fill(HIST("hK0Phi"), candidate.phi());
    registry.fill(HIST("hK0Eta"), candidate.eta());
    registry.fill(HIST("PiPlusTPC_K0"), postrack.pt(), postrack.tpcNSigmaPi());
    registry.fill(HIST("PiPlusTOF_K0"), postrack.pt(), postrack.tofNSigmaPi());
    registry.fill(HIST("PiMinusTPC_K0"), negtrack.pt(), negtrack.tpcNSigmaPi());
    registry.fill(HIST("PiMinusTOF_K0"), negtrack.pt(), negtrack.tofNSigmaPi());

    return true;
  }

  template <typename V0>
  bool isSelectedLambda(V0 const& candidate, float posZ, float posY, float posX)
  {
    bool isL = false;  // Is lambda candidate
    bool isAL = false; // Is anti-lambda candidate

    double mlambda = candidate.mLambda();
    double mantilambda = candidate.mAntiLambda();

    // separate the positive and negative V0 daughters
    auto postrack = candidate.template posTrack_as<FilteredTracks>();
    auto negtrack = candidate.template negTrack_as<FilteredTracks>();

    registry.fill(HIST("hLambdaCount"), 0.5);
    if (postrack.pt() < cfgPIDConfigs.cfgResoCuts->getData()[kPosTrackPt][iLambda] || negtrack.pt() < cfgPIDConfigs.cfgResoCuts->getData()[kNegTrackPt][iLambda])
      return false;

    registry.fill(HIST("hLambdaCount"), 1.5);
    if (mlambda > cfgPIDConfigs.cfgResoCuts->getData()[kMassMin][iLambda] && mlambda < cfgPIDConfigs.cfgResoCuts->getData()[kMassMax][iLambda])
      isL = true;
    if (mantilambda > cfgPIDConfigs.cfgResoCuts->getData()[kMassMin][iLambda] && mantilambda < cfgPIDConfigs.cfgResoCuts->getData()[kMassMax][iLambda])
      isAL = true;

    if (!isL && !isAL) {
      return false;
    }
    registry.fill(HIST("hLambdaCount"), 2.5);

    // Rapidity correction
    if (candidate.yLambda() > cfgPIDConfigs.cfgResoCuts->getData()[kRapidity][iLambda])
      return false;
    registry.fill(HIST("hLambdaCount"), 3.5);
    // DCA cuts for lambda and antilambda
    if (isL) {
      if (std::abs(candidate.dcapostopv()) < cfgPIDConfigs.cfgResoCuts->getData()[kDCAPosToPVMin][iLambda] || std::abs(candidate.dcanegtopv()) < cfgPIDConfigs.cfgResoCuts->getData()[kDCANegToPVMin][iLambda])
        return false;
    }
    if (isAL) {
      if (std::abs(candidate.dcapostopv()) < cfgPIDConfigs.cfgResoCuts->getData()[kDCANegToPVMin][iLambda] || std::abs(candidate.dcanegtopv()) < cfgPIDConfigs.cfgResoCuts->getData()[kDCAPosToPVMin][iLambda])
        return false;
    }
    registry.fill(HIST("hLambdaCount"), 4.5);
    if (std::abs(candidate.dcaV0daughters()) > cfgPIDConfigs.cfgResoSwitches->getData()[kDCABetDaug][iLambda])
      return false;
    registry.fill(HIST("hLambdaCount"), 5.5);
    // v0 radius cuts
    if (cfgPIDConfigs.cfgResoSwitches->getData()[kUseV0Radius][iLambda] && (candidate.v0radius() < cfgPIDConfigs.cfgResoCuts->getData()[kRadiusMin][iLambda] || candidate.v0radius() > cfgPIDConfigs.cfgResoCuts->getData()[kRadiusMax][iLambda]))
      return false;
    registry.fill(HIST("hLambdaCount"), 6.5);
    // cosine pointing angle cuts
    if (candidate.v0cosPA() < cfgPIDConfigs.cfgResoCuts->getData()[kCosPA][iLambda])
      return false;
    registry.fill(HIST("hLambdaCount"), 7.5);
    // Proper lifetime
    float cTauLambda = candidate.distovertotmom(posX, posY, posZ) * massLambda;
    if (cfgPIDConfigs.cfgResoSwitches->getData()[kUseProperLifetime][iLambda] && cTauLambda > cfgPIDConfigs.cfgResoCuts->getData()[kLifeTime][iLambda])
      return false;
    registry.fill(HIST("hLambdaCount"), 8.5);
    if (isL) {
      if (!selectionV0Daughter(postrack, kProtons) || !selectionV0Daughter(negtrack, kPions))
        return false;
    }
    if (isAL) {
      if (!selectionV0Daughter(postrack, kPions) || !selectionV0Daughter(negtrack, kProtons))
        return false;
    }
    registry.fill(HIST("hLambdaCount"), 9.5);

    if (!cfgPIDConfigs.cfgUseAntiLambda && isAL) { // Reject the track if it is antilambda
      return false;
    }

    registry.fill(HIST("hLambdaPhi"), candidate.phi());
    registry.fill(HIST("hLambdaEta"), candidate.eta());
    if (isL) {
      registry.fill(HIST("PrPlusTPC_La"), postrack.pt(), postrack.tpcNSigmaPr());
      registry.fill(HIST("PrPlusTOF_La"), postrack.pt(), postrack.tofNSigmaPr());
      registry.fill(HIST("PiMinusTPC_La"), negtrack.pt(), negtrack.tpcNSigmaPi());
      registry.fill(HIST("PiMinusTOF_La"), negtrack.pt(), negtrack.tofNSigmaPi());
    }
    if (cfgPIDConfigs.cfgUseAntiLambda && isAL) {
      registry.fill(HIST("PrMinusTPC_Al"), negtrack.pt(), negtrack.tpcNSigmaPr());
      registry.fill(HIST("PrMinusTOF_Al"), negtrack.pt(), negtrack.tofNSigmaPr());
      registry.fill(HIST("PiPlusTPC_Al"), postrack.pt(), postrack.tpcNSigmaPi());
      registry.fill(HIST("PiPlusTOF_Al"), postrack.pt(), postrack.tofNSigmaPi());
    }

    return true;
  }

  template <CorrelationContainer::CFStep step, typename TV0Tracks, typename TFT0s>
  void fillCorrelationsTPCFT0(TV0Tracks tracks1, TFT0s const& ft0, float posZ, float posY, float posX, int system, int corType, float eventWeight) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {
    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);

    float triggerWeight = 1.0f;
    // loop over all tracks
    for (auto const& track1 : tracks1) {
      double resoMass = -1;

      // 4 = kshort, 5 = lambda, 6 = phi
      if (cfgPIDConfigs.cfgPIDParticle == kK0) {
        if (!isSelectedK0(track1, posZ, posY, posX))
          continue; // Reject if called for K0 but V0 is not K0

        resoMass = track1.mK0Short();
      }

      if (cfgPIDConfigs.cfgPIDParticle == kLambda) {
        if (!isSelectedLambda(track1, posZ, posY, posX))
          continue; // Reject if called for Lambda but V0 is not lambda

        resoMass = track1.mLambda();
      }

      if (system == SameEvent) {
        if (corType == kFT0C) {
          registry.fill(HIST("Trig_hist_TPC_FT0C"), fSampleIndex, posZ, track1.pt(), resoMass, eventWeight * triggerWeight);
        } else if (corType == kFT0A) {
          registry.fill(HIST("Trig_hist_TPC_FT0A"), fSampleIndex, posZ, track1.pt(), resoMass, eventWeight * triggerWeight);
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
            if (cfgQaCheck) {
              registry.fill(HIST("Assoc_amp_same_TPC_FT0A"), chanelid, ampl);
              registry.fill(HIST("deltaEta_deltaPhi_same_TPC_FT0A"), deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
            }
            sameTpcFt0a->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), resoMass, deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
          } else if (corType == kFT0C) {
            if (cfgQaCheck) {
              registry.fill(HIST("Assoc_amp_same_TPC_FT0C"), chanelid, ampl);
              registry.fill(HIST("deltaEta_deltaPhi_same_TPC_FT0C"), deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
            }
            sameTpcFt0c->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), resoMass, deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
          }
        } else if (system == MixedEvent) {
          if (corType == kFT0A) {
            if (cfgQaCheck) {
              registry.fill(HIST("Assoc_amp_mixed_TPC_FT0A"), chanelid, ampl);
              registry.fill(HIST("deltaEta_deltaPhi_mixed_TPC_FT0A"), deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
            }
            mixedTpcFt0a->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), resoMass, deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
          } else if (corType == kFT0C) {
            if (cfgQaCheck) {
              registry.fill(HIST("Assoc_amp_mixed_TPC_FT0C"), chanelid, ampl);
              registry.fill(HIST("deltaEta_deltaPhi_mixed_TPC_FT0C"), deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
            }
            mixedTpcFt0c->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), resoMass, deltaPhi, deltaEta, ampl * eventWeight * triggerWeight);
          }
        }
      }
    }
  }

  template <CorrelationContainer::CFStep step, typename TFT0s>
  void fillCorrelationsFT0AFT0C(TFT0s const& ft0Col1, TFT0s const& ft0Col2, float posZ, int system, int multiplicity, float eventWeight) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
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

  bool isGoodRun(int runNumber)
  {
    for (const auto& ExcludedRun : cfgRunRemoveList.value) {
      if (runNumber == ExcludedRun) {
        return false;
      }
    }

    return true;
  }

  double massKaPlus = o2::constants::physics::MassKPlus;    // same as MassKMinus
  double massLambda = o2::constants::physics::MassLambda;   // same as MassLambda0 and MassLambda0Bar
  double massK0Short = o2::constants::physics::MassK0Short; // same as o2::constants::physics::MassK0 and o2::constants::physics::MassK0Bar

  void processSameTpcFt0a(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&, aod::V0Datas const& V0s)
  {
    registry.fill(HIST("hEventCount"), kFilteredEvents);

    if (cfgQaCheck) {
      eventSelectedIndividually(collision);
    }

    if (!collision.sel8())
      return;

    registry.fill(HIST("hEventCount"), kAfterSel8);

    if (!eventRct(collision, true))
      return;

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    int currentRunNumber = bc.runNumber();
    if (!cfgRunRemoveList.value.empty()) {
      if (!isGoodRun(currentRunNumber)) // Rejects runs if bad run number
        return;
    }

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

    double multiplicity = tracks.size();

    if (cfgQaCheck)
      registry.fill(HIST("Nch"), multiplicity);

    if (cfgStrictTrackCounter) {
      trackCounter(tracks, multiplicity);
    }

    if (cfgQaCheck) {
      registry.fill(HIST("Nch_corrected"), multiplicity);
    }

    if (multiplicity > cfgMaxMultForCorrelations || multiplicity < cfgMinMultForCorrelations) {
      return;
    }

    const auto& ft0 = collision.foundFT0();
    fillCorrelationsTPCFT0<CorrelationContainer::kCFStepReconstructed>(V0s, ft0, collision.posZ(), collision.posY(), collision.posX(), SameEvent, kFT0A, eventWeight);
  }
  PROCESS_SWITCH(CorrReso, processSameTpcFt0a, "Process same event for TPC-FT0 correlation", false);

  void processMixedTpcFt0a(FilteredCollisions const& collisions, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&, aod::V0Datas const& V0s)
  {

    auto getTracksSize = [&tracks, this](FilteredCollisions::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}, true};

    auto tracksTuple = std::make_tuple(V0s, tracks);
    Pair<FilteredCollisions, aod::V0Datas, FilteredTracks, MixedBinning> pairs{binningOnVtxAndMult, cfgMinMixEventNum, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, v0s1, collision2, tracks2] = *it;

      if (!collision1.sel8() || !collision2.sel8())
        continue;

      if (!eventRct(collision1, false) || !eventRct(collision2, false))
        continue;

      auto tracks1 = tracks.sliceByCached(o2::aod::track::collisionId, collision1.globalIndex(), this->cache);

      if (cfgUseAdditionalEventCut && !eventSelected(collision1, tracks1.size(), false))
        continue;
      if (cfgUseAdditionalEventCut && !eventSelected(collision2, tracks2.size(), false))
        continue;

      if (!(collision1.has_foundFT0() && collision2.has_foundFT0()))
        continue;

      registry.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      int currentRunNumber = bc.runNumber();
      if (!cfgRunRemoveList.value.empty()) {
        if (!isGoodRun(currentRunNumber)) // Rejects runs if bad run number
          return;
      }

      loadAlignParam(bc.timestamp());
      loadCorrection(bc.timestamp());
      float eventWeight = 1.0f;

      double multiplicity = tracks1.size();

      if (cfgStrictTrackCounter) {
        trackCounter(tracks1, multiplicity);
      }

      if (multiplicity > cfgMaxMultForCorrelations || multiplicity < cfgMinMultForCorrelations) {
        return;
      }

      const auto& ft0 = collision2.foundFT0();
      fillCorrelationsTPCFT0<CorrelationContainer::kCFStepReconstructed>(v0s1, ft0, collision1.posZ(), collision1.posY(), collision1.posX(), MixedEvent, kFT0A, eventWeight);
    }
  }
  PROCESS_SWITCH(CorrReso, processMixedTpcFt0a, "Process mixed events for TPC-FT0A correlation", false);

  void processSameTpcFt0c(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&, aod::V0Datas const& V0s)
  {
    registry.fill(HIST("hEventCount"), kFilteredEvents);

    if (cfgQaCheck) {
      eventSelectedIndividually(collision);
    }

    if (!collision.sel8())
      return;

    registry.fill(HIST("hEventCount"), kAfterSel8);

    if (!eventRct(collision, true))
      return;

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int currentRunNumber = bc.runNumber();
    if (!cfgRunRemoveList.value.empty()) {
      if (!isGoodRun(currentRunNumber)) // Rejects runs if bad run number
        return;
    }

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

    double multiplicity = tracks.size();

    if (cfgQaCheck)
      registry.fill(HIST("Nch"), multiplicity);

    if (cfgStrictTrackCounter) {
      trackCounter(tracks, multiplicity);
    }

    if (cfgQaCheck) {
      registry.fill(HIST("Nch_corrected"), multiplicity);
    }

    if (multiplicity > cfgMaxMultForCorrelations || multiplicity < cfgMinMultForCorrelations) {
      return;
    }

    fillCorrelationsTPCFT0<CorrelationContainer::kCFStepReconstructed>(V0s, ft0, collision.posZ(), collision.posY(), collision.posX(), SameEvent, kFT0C, 1.0f);
  }
  PROCESS_SWITCH(CorrReso, processSameTpcFt0c, "Process same event for TPC-FT0C correlation", false);

  void processMixedTpcFt0c(FilteredCollisions const& collisions, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&, aod::V0Datas const& V0s)
  {

    auto getTracksSize = [&tracks, this](FilteredCollisions::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}, true};

    auto tracksTuple = std::make_tuple(V0s, tracks);
    Pair<FilteredCollisions, aod::V0Datas, FilteredTracks, MixedBinning> pairs{binningOnVtxAndMult, cfgMinMixEventNum, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, v0s1, collision2, tracks2] = *it;
      if (!collision1.sel8() || !collision2.sel8())
        continue;

      if (!eventRct(collision1, false) || !eventRct(collision2, false))
        continue;

      auto tracks1 = tracks.sliceByCached(o2::aod::track::collisionId, collision1.globalIndex(), this->cache);

      if (cfgUseAdditionalEventCut && !eventSelected(collision1, tracks1.size(), false))
        continue;

      if (cfgUseAdditionalEventCut && !eventSelected(collision2, tracks2.size(), false))
        continue;

      if (!(collision1.has_foundFT0() && collision2.has_foundFT0()))
        continue;

      registry.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin
      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      int currentRunNumber = bc.runNumber();
      if (!cfgRunRemoveList.value.empty()) {
        if (!isGoodRun(currentRunNumber)) // Rejects runs if bad run number
          return;
      }
      loadAlignParam(bc.timestamp());
      loadCorrection(bc.timestamp());
      float eventWeight = 1.0f;

      const auto& ft0 = collision2.foundFT0();
      double multiplicity = tracks1.size();

      if (cfgStrictTrackCounter) {
        trackCounter(tracks, multiplicity);
      }

      if (multiplicity > cfgMaxMultForCorrelations || multiplicity < cfgMinMultForCorrelations) {
        return;
      }

      fillCorrelationsTPCFT0<CorrelationContainer::kCFStepReconstructed>(v0s1, ft0, collision1.posZ(), collision1.posY(), collision1.posX(), MixedEvent, kFT0C, eventWeight);
    }
  }
  PROCESS_SWITCH(CorrReso, processMixedTpcFt0c, "Process mixed events for TPC-FT0C correlation", false);

  void processSameFt0aFt0c(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::FT0s const&, aod::BCsWithTimestamps const&)
  {
    registry.fill(HIST("hEventCount"), kFilteredEvents);

    if (cfgQaCheck) {
      eventSelectedIndividually(collision);
    }

    if (!collision.sel8())
      return;

    registry.fill(HIST("hEventCount"), kAfterSel8);

    if (!eventRct(collision, true))
      return;

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int currentRunNumber = bc.runNumber();
    if (!cfgRunRemoveList.value.empty()) {
      if (!isGoodRun(currentRunNumber)) // Rejects runs if bad run number
        return;
    }

    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), true))
      return;

    if (!collision.has_foundFT0())
      return;

    loadAlignParam(bc.timestamp());
    loadGain(bc);
    loadCorrection(bc.timestamp());
    float eventWeight = 1.0f;

    const auto& ft0 = collision.foundFT0();

    fillYield(collision, tracks);

    registry.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin

    double multiplicity = tracks.size();

    if (cfgQaCheck)
      registry.fill(HIST("Nch"), multiplicity);

    if (cfgStrictTrackCounter) {
      trackCounter(tracks, multiplicity);
    }

    if (cfgQaCheck) {
      registry.fill(HIST("Nch_corrected"), multiplicity);
    }

    if (multiplicity > cfgMaxMultForCorrelations || multiplicity < cfgMinMultForCorrelations) {
      return;
    }

    fillCorrelationsFT0AFT0C<CorrelationContainer::kCFStepReconstructed>(ft0, ft0, collision.posZ(), SameEvent, multiplicity, eventWeight);
  }
  PROCESS_SWITCH(CorrReso, processSameFt0aFt0c, "Process same event for FT0A-FT0C correlation", true);

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

      if (!eventRct(collision1, false) || !eventRct(collision2, false))
        continue;

      if (cfgUseAdditionalEventCut && !eventSelected(collision1, tracks1.size(), false))
        continue;
      if (cfgUseAdditionalEventCut && !eventSelected(collision2, tracks2.size(), false))
        continue;

      if (!(collision1.has_foundFT0() && collision2.has_foundFT0()))
        continue;

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      int currentRunNumber = bc.runNumber();
      if (!cfgRunRemoveList.value.empty()) {
        if (!isGoodRun(currentRunNumber)) // Rejects runs if bad run number
          return;
      }
      loadAlignParam(bc.timestamp());
      loadCorrection(bc.timestamp());
      float eventWeight = 1.0f;

      const auto& ft0Col1 = collision1.foundFT0();
      const auto& ft0Col2 = collision2.foundFT0();

      double multiplicity = tracks1.size();

      if (cfgStrictTrackCounter) {
        trackCounter(tracks1, multiplicity);
      }

      if (multiplicity > cfgMaxMultForCorrelations || multiplicity < cfgMinMultForCorrelations) {
        return;
      }

      registry.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin

      fillCorrelationsFT0AFT0C<CorrelationContainer::kCFStepReconstructed>(ft0Col1, ft0Col2, collision1.posZ(), MixedEvent, multiplicity, eventWeight);
    }
  }
  PROCESS_SWITCH(CorrReso, processMixedFt0aFt0c, "Process mixed events for FT0A-FT0C correlation", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CorrReso>(cfgc),
  };
}
