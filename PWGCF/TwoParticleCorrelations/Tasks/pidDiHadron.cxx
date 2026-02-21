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

/// \file pidDiHadron.cxx
/// \brief di-hadron correlation of PID for O-O, Pb-Pb collisions
/// \author Preet Bhanjan Pati (preet.bhanjan.pati@cern.ch), Zhiyong Lu (zhiyong.lu@cern.ch)
/// \since  July/29/2025

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"
#include <CCDB/BasicCCDBManager.h>

#include "TF1.h"
#include "TRandom3.h"
#include <TPDGCode.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// define the filtered collisions and tracks
#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

template <typename T, typename P>
auto readMatrix(Array2D<T> const& mat, P& array, int rowStart, int rowEnd, int colStart, int colEnd)
{
  for (auto i = rowStart; i < rowEnd; ++i) {
    for (auto j = colStart; j < colEnd; ++j) {
      array[i][j] = mat(i, j);
    }
  }

  return;
}

static constexpr float LongArrayFloat[3][20] = {{1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}, {2.1, 2.2, 2.3, -2.1, -2.2, -2.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}, {3.1, 3.2, 3.3, -3.1, -3.2, -3.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}};
static constexpr int LongArrayInt[3][20] = {{1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1}, {2, 2, 2, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1}, {3, 3, 3, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1}};

struct PidDiHadron {
  o2::aod::ITSResponse itsResponse;
  Service<ccdb::BasicCCDBManager> ccdb;

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "minimum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "maximum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta cut")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5f, "max chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgTpcCluster, float, 50.0f, "minimum TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgTpcCrossRows, float, 70.0f, "minimum TPC crossed rows")
  O2_DEFINE_CONFIGURABLE(cfgTofPtCut, float, 0.5f, "Minimum pt to use TOF N-sigma")
  O2_DEFINE_CONFIGURABLE(cfgTpcCut, float, 3.0f, "TPC N-sigma cut for pions, kaons, protons")
  O2_DEFINE_CONFIGURABLE(cfgITScluster, float, 5.0f, "minimum ITS clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyMax, int, 2000, "Minimum occupancy cut")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyMin, int, 0, "Maximum occupancy cut")
  O2_DEFINE_CONFIGURABLE(cfgFakeKaonCut, float, 0.1f, "Maximum difference in measured momentum and TPC inner ring momentum of particle")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "max DCA to vertex z")
  O2_DEFINE_CONFIGURABLE(cfgCutMerging, float, 0.0, "Merging cut on track merge")
  O2_DEFINE_CONFIGURABLE(cfgSelCollByNch, bool, true, "Select collisions by Nch or centrality")
  O2_DEFINE_CONFIGURABLE(cfgCutMultMin, int, 0, "Minimum multiplicity for collision")
  O2_DEFINE_CONFIGURABLE(cfgCutMultMax, int, 10, "Maximum multiplicity for collision")
  O2_DEFINE_CONFIGURABLE(cfgCutCentMin, float, 60.0f, "Minimum centrality for collision")
  O2_DEFINE_CONFIGURABLE(cfgCutCentMax, float, 80.0f, "Maximum centrality for collision")
  O2_DEFINE_CONFIGURABLE(cfgMixEventNumMin, int, 5, "Minimum number of events to mix")
  O2_DEFINE_CONFIGURABLE(cfgRadiusLow, float, 0.8, "Low radius for merging cut")
  O2_DEFINE_CONFIGURABLE(cfgRadiusHigh, float, 2.5, "High radius for merging cut")
  O2_DEFINE_CONFIGURABLE(cfgSampleSize, double, 10, "Sample size for mixed event")
  O2_DEFINE_CONFIGURABLE(cfgCentEstimator, int, 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FT0A")
  O2_DEFINE_CONFIGURABLE(cfgCentTableUnavailable, bool, false, "if a dataset does not provide centrality information")
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, false, "Use additional event cut on mult correlations")
  O2_DEFINE_CONFIGURABLE(cfgV0AT0Acut, int, 5, "V0AT0A cut")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgCentralityWeight, std::string, "", "CCDB path to centrality weight object")
  O2_DEFINE_CONFIGURABLE(cfgLocalEfficiency, bool, false, "Use local efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgVerbosity, bool, false, "Verbose output")
  O2_DEFINE_CONFIGURABLE(cfgUseEventWeights, bool, false, "Use event weights for mixed event")
  O2_DEFINE_CONFIGURABLE(cfgUsePtOrder, bool, false, "enable trigger pT < associated pT cut")
  O2_DEFINE_CONFIGURABLE(cfgUsePtOrderInMixEvent, bool, false, "enable trigger pT < associated pT cut in mixed event")
  O2_DEFINE_CONFIGURABLE(cfgUseItsPID, bool, true, "Use ITS PID for particle identification")
  O2_DEFINE_CONFIGURABLE(cfgUseOnlyTPC, bool, true, "Use only TPC PID for daughter selection")
  O2_DEFINE_CONFIGURABLE(cfgPIDParticle, int, 0, "1 = pion, 2 = kaon, 3 = proton, 4 = kshort, 5 = lambda, 6 = phi, 0 for no PID")
  O2_DEFINE_CONFIGURABLE(cfgGetNsigmaQA, bool, true, "Get QA histograms for selection of pions, kaons, and protons")
  O2_DEFINE_CONFIGURABLE(cfgGetdEdx, bool, true, "Get QA histograms for selection of pions, kaons, and protons")
  O2_DEFINE_CONFIGURABLE(cfgUseAntiLambda, bool, true, "Use AntiLambda candidates for analysis")
  O2_DEFINE_CONFIGURABLE(cfgPIDUseRejection, bool, true, "Turn off and on the exclusion criteria for PID determination")

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
    O2_DEFINE_CONFIGURABLE(cfgDCAxyNSigma, float, 7, "Cut on number of sigma deviations from expected DCA in the transverse direction");
    O2_DEFINE_CONFIGURABLE(cfgDCAxy, std::string, "(0.0026+0.005/(x^1.01))", "Functional form of pt-dependent DCAxy cut");
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
    TF1* fPtDepDCAxy = nullptr;
  } cfgFuncParas;

  SliceCache cache;

  ConfigurableAxis axisVertex{"axisVertex", {10, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260}, "multiplicity axis for histograms"};
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, "centrality axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {48, -2.4, 2.4}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt associated axis for histograms"};
  ConfigurableAxis axisVtxMix{"axisVtxMix", {VARIABLE_WIDTH, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, "vertex axis for mixed event histograms"};
  ConfigurableAxis axisMultMix{"axisMultMix", {VARIABLE_WIDTH, 0, 10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260}, "multiplicity / centrality axis for mixed event histograms"};
  ConfigurableAxis axisSample{"axisSample", {cfgSampleSize, 0, cfgSampleSize}, "sample axis for histograms"};
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {80, -5, 5}, "nsigmaTPC axis"};
  ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {80, -5, 5}, "nsigmaTOF axis"};
  ConfigurableAxis axisNsigmaITS{"axisNsigmaITS", {80, -5, 5}, "nsigmaITS axis"};
  ConfigurableAxis axisTpcSignal{"axisTpcSignal", {250, 0, 250}, "dEdx axis for TPC"};
  ConfigurableAxis axisSigma{"axisSigma", {200, 0, 20}, "sigma axis for TPC"};

  Configurable<LabeledArray<int>> cfgUseEventCuts{"cfgUseEventCuts", {LongArrayInt[0], 15, 1, {"Filtered Events", "Sel8", "kNoTimeFrameBorder", "kNoITSROFrameBorder", "kNoSameBunchPileup", "kIsGoodZvtxFT0vsPV", "kNoCollInTimeRangeStandard", "kIsGoodITSLayersAll", "kNoCollInRofStandard", "kNoHighMultCollInPrevRof", "Occupancy", "Multcorrelation", "T0AV0ACut", "kIsVertexITSTPC", "kTVXinTRD"}, {"EvCuts"}}, "Labeled array (int) for various cuts on resonances"};
  Configurable<LabeledArray<float>> nSigmas{"nSigmas", {LongArrayFloat[0], 6, 3, {"UpCut_pi", "UpCut_ka", "UpCut_pr", "LowCut_pi", "LowCut_ka", "LowCut_pr"}, {"TPC", "TOF", "ITS"}}, "Labeled array for n-sigma values for TPC, TOF, ITS for pions, kaons, protons (positive and negative)"};
  Configurable<LabeledArray<float>> resonanceCuts{"resonanceCuts", {LongArrayFloat[0], 12, 3, {"cos_PAs", "massMin", "massMax", "PosTrackPt", "NegTrackPt", "DCAPosToPVMin", "DCANegToPVMin", "Lifetime", "RadiusMin", "RadiusMax", "Rapidity", "ArmPodMinVal"}, {"K0", "Lambda", "Phi"}}, "Labeled array (float) for various cuts on resonances"};
  Configurable<LabeledArray<int>> resonanceSwitches{"resonanceSwitches", {LongArrayInt[0], 6, 3, {"UseCosPA", "NMassBins", "DCABetDaug", "UseProperLifetime", "UseV0Radius", "UseArmPodCut"}, {"K0", "Lambda", "Phi"}}, "Labeled array (int) for various cuts on resonances"};

  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {1, 0, 1}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {1, 0, 1}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {1, 0, 1}, "pt axis for efficiency histograms"};

  // make the filters and cuts.
  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgCutVertex);
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtPOIMin) && (aod::track::pt < cfgCutPtPOIMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t)true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using FilteredCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSel, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::Mults>>;
  using FilteredTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;
  using FilteredMcTracks = soa::Filtered<soa::Join<aod::Tracks, aod::McTrackLabels, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;
  using V0TrackCandidate = aod::V0Datas;

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  Preslice<aod::V0Datas> perCollisionV0 = aod::v0::collisionId;

  // Corrections
  TH3D* mEfficiency = nullptr;
  TH1D* mCentralityWeight = nullptr;
  bool correctionsLoaded = false;

  // Define the outputs
  OutputObj<CorrelationContainer> same{"sameEvent"};
  OutputObj<CorrelationContainer> mixed{"mixedEvent"};
  OutputObj<CorrelationContainer> sameReso{"sameEventReso"};
  OutputObj<CorrelationContainer> mixedReso{"mixedEventReso"};
  HistogramRegistry histos{"histos"};

  // define global variables
  TRandom3* gRandom = new TRandom3();
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
  enum EventCutTypes {
    kFilteredEvents = 0,
    kAfterSel8,
    kUseNoTimeFrameBorder,
    kUseNoITSROFrameBorder,
    kUseNoSameBunchPileup,
    kUseGoodZvtxFT0vsPV,
    kUseNoCollInTimeRangeStandard,
    kUseGoodITSLayersAll,
    kUseNoCollInRofStandard,
    kUseNoHighMultCollInPrevRof,
    kUseOccupancy,
    kUseMultCorrCut,
    kUseT0AV0ACut,
    kUseVertexITSTPC,
    kUseTVXinTRD,
    kNEventCuts
  };
  enum CentEstimators {
    kCentFT0C = 0,
    kCentFT0CVariant1,
    kCentFT0M,
    kCentFV0A,
    // Count the total number of enum
    kCount_CentEstimators
  };
  enum EventType {
    SameEvent = 1,
    MixedEvent = 3
  };
  enum DetectorType {
    kTPC = 0,
    kTOF,
    kITS
  };
  enum EventCutType {
    kEvCut1 = 0,
    kNEvCutTypes = 1
  };

  std::array<std::array<int, 1>, 15> eventCuts;
  std::array<std::array<float, 3>, 12> resoCutVals;
  std::array<std::array<int, 3>, 7> resoSwitchVals;
  std::array<float, 6> tofNsigmaCut;
  std::array<float, 6> itsNsigmaCut;
  std::array<float, 6> tpcNsigmaCut;

  // persistent caches
  std::vector<float> efficiencyAssociatedCache;

  void init(InitContext&)
  {
    readMatrix(resonanceCuts->getData(), resoCutVals, kCosPA, kNParticleCuts, iK0, NResoParticles);
    readMatrix(resonanceSwitches->getData(), resoSwitchVals, kUseCosPA, kNParticleSwitches, iK0, NResoParticles);
    readMatrix(cfgUseEventCuts->getData(), eventCuts, kFilteredEvents, kNEventCuts, kEvCut1, kNEvCutTypes);

    tpcNsigmaCut[iPionUp] = nSigmas->getData()[iPionUp][kTPC];
    tpcNsigmaCut[iKaonUp] = nSigmas->getData()[iKaonUp][kTPC];
    tpcNsigmaCut[iProtonUp] = nSigmas->getData()[iProtonUp][kTPC];
    tpcNsigmaCut[iPionLow] = nSigmas->getData()[iPionLow][kTPC];
    tpcNsigmaCut[iKaonLow] = nSigmas->getData()[iKaonLow][kTPC];
    tpcNsigmaCut[iProtonLow] = nSigmas->getData()[iProtonLow][kTPC];

    tofNsigmaCut[iPionUp] = nSigmas->getData()[iPionUp][kTOF];
    tofNsigmaCut[iKaonUp] = nSigmas->getData()[iKaonUp][kTOF];
    tofNsigmaCut[iProtonUp] = nSigmas->getData()[iProtonUp][kTOF];
    tofNsigmaCut[iPionLow] = nSigmas->getData()[iPionLow][kTOF];
    tofNsigmaCut[iKaonLow] = nSigmas->getData()[iKaonLow][kTOF];
    tofNsigmaCut[iProtonLow] = nSigmas->getData()[iProtonLow][kTOF];

    itsNsigmaCut[iPionUp] = nSigmas->getData()[iPionUp][kITS];
    itsNsigmaCut[iKaonUp] = nSigmas->getData()[iKaonUp][kITS];
    itsNsigmaCut[iProtonUp] = nSigmas->getData()[iProtonUp][kITS];
    itsNsigmaCut[iPionLow] = nSigmas->getData()[iPionLow][kITS];
    itsNsigmaCut[iKaonLow] = nSigmas->getData()[iKaonLow][kITS];
    itsNsigmaCut[iProtonLow] = nSigmas->getData()[iProtonLow][kITS];

    AxisSpec axisK0Mass = {resoSwitchVals[kMassBins][iK0], resoCutVals[kMassMin][iK0], resoCutVals[kMassMax][iK0]};
    AxisSpec axisLambdaMass = {resoSwitchVals[kMassBins][iLambda], resoCutVals[kMassMin][iLambda], resoCutVals[kMassMax][iLambda]};
    AxisSpec axisPhiMass = {resoSwitchVals[kMassBins][iPhi], resoCutVals[kMassMin][iPhi], resoCutVals[kMassMax][iPhi]};

    if (cfgCentTableUnavailable && !cfgSelCollByNch) {
      LOGF(fatal, "Centrality table is unavailable, cannot select collisions by centrality");
    }
    const AxisSpec axisPhi{72, 0.0, constants::math::TwoPI, "#varphi"};
    const AxisSpec axisEta{40, -1., 1., "#eta"};

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    LOGF(info, "Starting init");

    // Creating mass axis depending on particle - 4 = kshort, 5 = lambda, 6 = phi
    AxisSpec massAxisReso = {10, 0, 1, "mass"};
    if (cfgPIDParticle == kK0)
      massAxisReso = {resoSwitchVals[kMassBins][iK0], resoCutVals[kMassMin][iK0], resoCutVals[kMassMax][iK0], "M_{#pi^{+}#pi^{-}} (GeV/c^{2})"};
    if (cfgPIDParticle == kLambda)
      massAxisReso = {resoSwitchVals[kMassBins][iLambda], resoCutVals[kMassMin][iLambda], resoCutVals[kMassMax][iLambda], "M_{p#pi^{-}} (GeV/c^{2})"};
    if (cfgPIDParticle == kPhi)
      massAxisReso = {resoSwitchVals[kMassBins][iPhi], resoCutVals[kMassMin][iPhi], resoCutVals[kMassMax][iPhi], "M_{K^{+}K^{-}} (GeV/c^{2})"};

    // Event Counter
    if ((doprocessSame || doprocessSameReso || doprocessMC) && cfgUseAdditionalEventCut) {
      histos.add("hEventCount", "Number of Events;; Count", {HistType::kTH1D, {{15, -0.5, 14.5}}});
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kFilteredEvents + 1, "Filtered event");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kAfterSel8 + 1, "After sel8");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoTimeFrameBorder + 1, "kNoTimeFrameBorder");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoITSROFrameBorder + 1, "kNoITSROFrameBorder");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoSameBunchPileup + 1, "kNoSameBunchPileup");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseGoodZvtxFT0vsPV + 1, "kIsGoodZvtxFT0vsPV");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoCollInTimeRangeStandard + 1, "kNoCollInTimeRangeStandard");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseGoodITSLayersAll + 1, "kIsGoodITSLayersAll");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoCollInRofStandard + 1, "kNoCollInRofStandard");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoHighMultCollInPrevRof + 1, "kNoHighMultCollInPrevRof");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseOccupancy + 1, "Occupancy Cut");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseMultCorrCut + 1, "MultCorrelation Cut");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseT0AV0ACut + 1, "T0AV0A cut");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseVertexITSTPC + 1, "kIsVertexITSTPC");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseTVXinTRD + 1, "kTVXinTRD");
    }

    if (cfgPIDParticle == kK0) { // For K0
      histos.add("PiPlusTPC_K0", "", {HistType::kTH2D, {{axisPt, axisNsigmaTPC}}});
      histos.add("PiMinusTPC_K0", "", {HistType::kTH2D, {{axisPt, axisNsigmaTPC}}});
      histos.add("PiPlusTOF_K0", "", {HistType::kTH2D, {{axisPt, axisNsigmaTOF}}});
      histos.add("PiMinusTOF_K0", "", {HistType::kTH2D, {{axisPt, axisNsigmaTOF}}});
      histos.add("hK0Phi", "", {HistType::kTH1D, {axisPhi}});
      histos.add("hK0Eta", "", {HistType::kTH1D, {axisEta}});

      histos.add("hK0Count", "Number of K0;; Count", {HistType::kTH1D, {{11, 0, 11}}});
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(1, "K0 candidates");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(2, "Daughter pt");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(3, "Mass cut");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(4, "Rapidity cut");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(5, "DCA to PV");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(6, "DCA between daughters");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(7, "V0radius");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(8, "CosPA");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(9, "Proper lifetime");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(10, "ArmenterosPod");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(11, "Daughter track selection");
    }
    if (cfgPIDParticle == kLambda) { // For Lambda
      histos.add("PrPlusTPC_L", "", {HistType::kTH2D, {{axisPt, axisNsigmaTPC}}});
      histos.add("PiMinusTPC_L", "", {HistType::kTH2D, {{axisPt, axisNsigmaTPC}}});
      histos.add("PrPlusTOF_L", "", {HistType::kTH2D, {{axisPt, axisNsigmaTOF}}});
      histos.add("PiMinusTOF_L", "", {HistType::kTH2D, {{axisPt, axisNsigmaTOF}}});
      histos.add("hLambdaPhi", "", {HistType::kTH1D, {axisPhi}});
      histos.add("hLambdaEta", "", {HistType::kTH1D, {axisEta}});

      histos.add("hLambdaCount", "Number of Lambda;; Count", {HistType::kTH1D, {{10, 0, 10}}});
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(1, "Lambda candidates");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(2, "Daughter pt");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(3, "Mass cut");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(4, "Rapidity cut");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(5, "DCA to PV");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(6, "DCA between daughters");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(7, "V0radius");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(8, "CosPA");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(9, "Proper lifetime");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(10, "Daughter track selection");
    }
    if (cfgPIDParticle == kPhi) { // For Phi
      histos.add("KaPlusTPC", "", {HistType::kTH2D, {{axisPt, axisNsigmaTPC}}});
      histos.add("KaMinusTPC", "", {HistType::kTH2D, {{axisPt, axisNsigmaTPC}}});
      histos.add("KaPlusTOF", "", {HistType::kTH2D, {{axisPt, axisNsigmaTOF}}});
      histos.add("KaMinusTOF", "", {HistType::kTH2D, {{axisPt, axisNsigmaTOF}}});
      histos.add("hPhiPhi", "", {HistType::kTH1D, {axisPhi}});
      histos.add("hPhiEta", "", {HistType::kTH1D, {axisEta}});
      histos.add("hPhiMass_sparse", "", {HistType::kTHnSparseD, {{axisPhiMass, axisPt, axisMultiplicity}}});

      histos.add("hPhiCount", "Number of Phi;; Count", {HistType::kTH1D, {{5, 0, 5}}});
      histos.get<TH1>(HIST("hPhiCount"))->GetXaxis()->SetBinLabel(1, "Phi candidates");
      histos.get<TH1>(HIST("hPhiCount"))->GetXaxis()->SetBinLabel(2, "Daughter track selection");
      histos.get<TH1>(HIST("hPhiCount"))->GetXaxis()->SetBinLabel(3, "Fake Kaon");
      histos.get<TH1>(HIST("hPhiCount"))->GetXaxis()->SetBinLabel(4, "CosPA");
      histos.get<TH1>(HIST("hPhiCount"))->GetXaxis()->SetBinLabel(5, "Rapidity cut");
    }

    // Multiplicity correlation cuts
    if (eventCuts[kUseMultCorrCut][kEvCut1]) {
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
    if (eventCuts[kUseT0AV0ACut][kEvCut1]) {
      cfgFuncParas.fT0AV0AMean = new TF1("fT0AV0AMean", "[0]+[1]*x", 0, 200000);
      cfgFuncParas.fT0AV0AMean->SetParameters(-1601.0581, 9.417652e-01);
      cfgFuncParas.fT0AV0ASigma = new TF1("fT0AV0ASigma", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 200000);
      cfgFuncParas.fT0AV0ASigma->SetParameters(463.4144, 6.796509e-02, -9.097136e-07, 7.971088e-12, -2.600581e-17);
    }

    cfgFuncParas.fPtDepDCAxy = new TF1("ptDepDCAxy", Form("[0]*%s", cfgFuncParas.cfgDCAxy->c_str()), 0.001, 100);
    cfgFuncParas.fPtDepDCAxy->SetParameter(0, cfgFuncParas.cfgDCAxyNSigma);
    LOGF(info, "DCAxy pt-dependence function: %s", Form("[0]*%s", cfgFuncParas.cfgDCAxy->c_str()));

    std::string hCentTitle = "Centrality distribution, Estimator " + std::to_string(cfgCentEstimator);
    // Make histograms to check the distributions after cuts
    if (doprocessSame || doprocessSameReso) {
      histos.add("deltaEta_deltaPhi_same", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}}); // check to see the delta eta and delta phi distribution
      histos.add("Nch", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}});
      histos.add("Nch_used", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}}); // histogram to see how many events are in the same and mixed event
      histos.add("Centrality", hCentTitle.c_str(), {HistType::kTH1D, {axisCentrality}});
      histos.add("CentralityWeighted", hCentTitle.c_str(), {HistType::kTH1D, {{100, 0, 100}}});
      histos.add("Centrality_used", hCentTitle.c_str(), {HistType::kTH1D, {axisCentrality}}); // histogram to see how many events are in the same and mixed event
      histos.add("zVtx", "zVtx", {HistType::kTH1D, {axisVertex}});
      histos.add("zVtx_used", "zVtx_used", {HistType::kTH1D, {axisVertex}});

      if (cfgPIDParticle == kCharged || cfgPIDParticle == kPions || cfgPIDParticle == kKaons || cfgPIDParticle == kProtons) {
        histos.add("Phi", "Phi", {HistType::kTH1D, {axisPhi}});
        histos.add("Eta", "Eta", {HistType::kTH1D, {axisEta}});
        histos.add("EtaCorrected", "EtaCorrected", {HistType::kTH1D, {axisEta}});
        histos.add("pT", "pT", {HistType::kTH1D, {axisPtTrigger}});
        histos.add("pTCorrected", "pTCorrected", {HistType::kTH1D, {axisPtTrigger}});
        histos.add("Trig_hist", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});

        if (cfgGetNsigmaQA) {
          if (!cfgUseItsPID) {
            histos.add("TofTpcNsigma_before", "", {HistType::kTHnSparseD, {{axisNsigmaTPC, axisNsigmaTOF, axisPt}}});
            histos.add("TofTpcNsigma_after", "", {HistType::kTHnSparseD, {{axisNsigmaTPC, axisNsigmaTOF, axisPt}}});

            if (cfgGetdEdx) {
              histos.add("TpcdEdx_ptwise", "", {HistType::kTHnSparseD, {{axisPt, axisTpcSignal, axisNsigmaTOF}}});
              histos.add("ExpTpcdEdx_ptwise", "", {HistType::kTHnSparseD, {{axisPt, axisTpcSignal, axisNsigmaTOF}}});
              histos.add("ExpSigma_ptwise", "", {HistType::kTHnSparseD, {{axisPt, axisSigma, axisNsigmaTOF}}});
            }
          } // TPC-TOF PID QA hists
          if (cfgUseItsPID) {
            histos.add("TofItsNsigma_before", "", {HistType::kTHnSparseD, {{axisNsigmaITS, axisNsigmaTOF, axisPt}}});
            histos.add("TofItsNsigma_after", "", {HistType::kTHnSparseD, {{axisNsigmaITS, axisNsigmaTOF, axisPt}}});
          } // ITS-TOF PID QA hists
        } // end of PID QA hists
      }

      if (cfgPIDParticle == kK0 || cfgPIDParticle == kLambda || cfgPIDParticle == kPhi) {
        histos.add("Trig_histReso", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger, massAxisReso}}});
      }
    }
    if (doprocessMixed || doprocessMixedReso) {
      histos.add("deltaEta_deltaPhi_mixed", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}});
    }
    if (doprocessMC) {
      histos.add("hNsigmaPionTruePositives", "hNsigmaPionTruePositives", {HistType::kTH1D, {axisPt}});     // Fraction of particles that are pions and selected as pions
      histos.add("hNsigmaKaonTruePositives", "hNsigmaKaonTruePositives", {HistType::kTH1D, {axisPt}});     // Fraction of particles that are kaons and selected as kaons
      histos.add("hNsigmaProtonTruePositives", "hNsigmaProtonTruePositives", {HistType::kTH1D, {axisPt}}); // Fraction of particles that are protons and selected as protons

      histos.add("hNsigmaPionSelected", "hNsigmaPionSelected", {HistType::kTH1D, {axisPt}});
      histos.add("hNsigmaKaonSelected", "hNsigmaKaonSelected", {HistType::kTH1D, {axisPt}});
      histos.add("hNsigmaProtonSelected", "hNsigmaProtonSelected", {HistType::kTH1D, {axisPt}});

      histos.add("hNsigmaPionTrue", "hNsigmaPionTrue", {HistType::kTH1D, {axisPt}});     // All true pions from MC
      histos.add("hNsigmaKaonTrue", "hNsigmaKaonTrue", {HistType::kTH1D, {axisPt}});     // All true kaons from MC
      histos.add("hNsigmaProtonTrue", "hNsigmaProtonTrue", {HistType::kTH1D, {axisPt}}); // All true protons from MC

      histos.add("TpcdEdx_ptwise_pi", "", {HistType::kTHnSparseD, {{axisPt, axisTpcSignal, axisNsigmaTOF}}});
      histos.add("TpcdEdx_ptwise_ka", "", {HistType::kTHnSparseD, {{axisPt, axisTpcSignal, axisNsigmaTOF}}});
      histos.add("TpcdEdx_ptwise_pr", "", {HistType::kTHnSparseD, {{axisPt, axisTpcSignal, axisNsigmaTOF}}});

      histos.add("ExpTpcdEdx_ptwise_pi", "", {HistType::kTHnSparseD, {{axisPt, axisTpcSignal, axisNsigmaTOF}}});
      histos.add("ExpTpcdEdx_ptwise_ka", "", {HistType::kTHnSparseD, {{axisPt, axisTpcSignal, axisNsigmaTOF}}});
      histos.add("ExpTpcdEdx_ptwise_pr", "", {HistType::kTHnSparseD, {{axisPt, axisTpcSignal, axisNsigmaTOF}}});

    }

    histos.add("eventcount", "bin", {HistType::kTH1F, {{4, 0, 4, "bin"}}}); // histogram to see how many events are in the same and mixed event

    LOGF(info, "Initializing correlation container");
    std::vector<AxisSpec> corrAxis = {{axisSample, "Sample"},
                                      {axisVertex, "z-vtx (cm)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisDeltaEta, "#Delta#eta"}};

    std::vector<AxisSpec> corrAxisReso = {{axisSample, "Sample"},
                                          {axisVertex, "z-vtx (cm)"},
                                          {axisPtTrigger, "p_{T} (GeV/c)"},
                                          {massAxisReso},
                                          {axisDeltaPhi, "#Delta#varphi (rad)"},
                                          {axisDeltaEta, "#Delta#eta"}};
    std::vector<AxisSpec> effAxis = {
      {axisEtaEfficiency, "#eta"},
      {axisPtEfficiency, "p_{T} (GeV/c)"},
      {axisVertexEfficiency, "z-vtx (cm)"},
    };
    std::vector<AxisSpec> userAxis;

    if (cfgPIDParticle == kCharged || cfgPIDParticle == kPions || cfgPIDParticle == kKaons || cfgPIDParticle == kProtons) {
      same.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, userAxis));
      mixed.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, userAxis));
    }

    if (cfgPIDParticle == kK0 || cfgPIDParticle == kLambda || cfgPIDParticle == kPhi) {
      sameReso.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxisReso, effAxis, userAxis));
      mixedReso.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxisReso, effAxis, userAxis));
    }
    LOGF(info, "End of init");
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

  template <typename TTrack>
  bool trackSelected(TTrack track)
  {
    if (cfgFuncParas.cfgDCAxyNSigma && (std::fabs(track.dcaXY()) > cfgFuncParas.fPtDepDCAxy->Eval(track.pt())))
      return false;
    return ((track.tpcNClsFound() >= cfgTpcCluster) && (track.tpcNClsCrossedRows() >= cfgTpcCrossRows) && (track.itsNCls() >= cfgITScluster));
  }

  template <typename TTrack>
  bool selectionV0Daughter(TTrack const& track, int pid)
  {
    if (!(track.itsNCls() > cfgITScluster))
      return 0;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < cfgTpcCluster)
      return false;
    if (!(track.tpcNClsCrossedRows() > cfgTpcCrossRows))
      return 0;

    if (cfgUseOnlyTPC) {
      if (pid == kPions && std::abs(track.tpcNSigmaPi()) > cfgTpcCut)
        return false;
      if (pid == kKaons && std::abs(track.tpcNSigmaKa()) > cfgTpcCut)
        return false;
      if (pid == kProtons && std::abs(track.tpcNSigmaPr()) > cfgTpcCut)
        return false;
    } else {
      int partIndex = getNsigmaPID(track);
      int pidIndex = partIndex; // 1 = pion, 2 = kaon, 3 = proton
      if (pidIndex != pid)
        return false;
    }

    return true;
  }

  template <typename TTrack>
  int getNsigmaPID(TTrack track)
  {
    // Computing Nsigma arrays for pion, kaon, and protons
    std::array<float, 3> nSigmaTPC = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::array<float, 3> nSigmaTOF = {track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};
    std::array<float, 3> nSigmaITS = {itsResponse.nSigmaITS<o2::track::PID::Pion>(track), itsResponse.nSigmaITS<o2::track::PID::Kaon>(track), itsResponse.nSigmaITS<o2::track::PID::Proton>(track)};
    int pid = -1; // -1 = not identified, 1 = pion, 2 = kaon, 3 = proton

    std::array<float, 3> nSigmaToUse = cfgUseItsPID ? nSigmaITS : nSigmaTPC;             // Choose which nSigma to use: TPC or ITS
    std::array<float, 6> detectorNsigmaCut = cfgUseItsPID ? itsNsigmaCut : tpcNsigmaCut; // Choose which nSigma to use: TPC or ITS

    bool isPion, isKaon, isProton;
    bool isDetectedPion = nSigmaToUse[iPionUp] < detectorNsigmaCut[iPionUp] && nSigmaToUse[iPionUp] > detectorNsigmaCut[iPionLow];
    bool isDetectedKaon = nSigmaToUse[iKaonUp] < detectorNsigmaCut[iKaonUp] && nSigmaToUse[iKaonUp] > detectorNsigmaCut[iKaonLow];
    bool isDetectedProton = nSigmaToUse[iProtonUp] < detectorNsigmaCut[iProtonUp] && nSigmaToUse[iProtonUp] > detectorNsigmaCut[iProtonLow];

    bool isTofPion = nSigmaTOF[iPionUp] < tofNsigmaCut[iPionUp] && nSigmaTOF[iPionUp] > tofNsigmaCut[iPionLow];
    bool isTofKaon = nSigmaTOF[iKaonUp] < tofNsigmaCut[iKaonUp] && nSigmaTOF[iKaonUp] > tofNsigmaCut[iKaonLow];
    bool isTofProton = nSigmaTOF[iProtonUp] < tofNsigmaCut[iProtonUp] && nSigmaTOF[iProtonUp] > tofNsigmaCut[iProtonLow];

    if (track.pt() > cfgTofPtCut && !track.hasTOF()) {
      return -1;
    } else if (track.pt() > cfgTofPtCut && track.hasTOF()) {
      isPion = isTofPion && isDetectedPion;
      isKaon = isTofKaon && isDetectedKaon;
      isProton = isTofProton && isDetectedProton;
    } else {
      isPion = isDetectedPion;
      isKaon = isDetectedKaon;
      isProton = isDetectedProton;
    }

    if (cfgPIDUseRejection && ((isPion && isKaon) || (isPion && isProton) || (isKaon && isProton))) {
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

  bool getCentralityWeight(float& weightCent, const float centrality)
  {
    float weight = 1.;
    if (mCentralityWeight)
      weight = mCentralityWeight->GetBinContent(mCentralityWeight->FindBin(centrality));
    else
      weight = 1.0;
    if (weight == 0)
      return false;
    weightCent = weight;
    return true;
  }

  // fill multiple histograms
  template <typename TCollision, typename TTracks>
  void fillYield(TCollision collision, TTracks tracks) // function to fill the yield and etaphi histograms.
  {
    float weff1 = 1;
    float vtxz = collision.posZ();
    for (auto const& track1 : tracks) {
      if (!trackSelected(track1))
        continue;
      if (cfgPIDParticle && getNsigmaPID(track1) != cfgPIDParticle)
        continue; // if PID is selected, check if the track has the right PID
      if (!getEfficiencyCorrection(weff1, track1.eta(), track1.pt(), vtxz))
        continue;
      histos.fill(HIST("Phi"), RecoDecay::constrainAngle(track1.phi(), 0.0));
      histos.fill(HIST("Eta"), track1.eta());
      histos.fill(HIST("EtaCorrected"), track1.eta(), weff1);
      histos.fill(HIST("pT"), track1.pt());
      histos.fill(HIST("pTCorrected"), track1.pt(), weff1);
    }
  }

  template <typename TTrack>
  void fillNsigmaAfterCut(TTrack track1, Int_t pid) // function to fill the QA after Nsigma selection
  {
    switch (pid) {
      case 1: // For Pions
        if (!cfgUseItsPID)
          histos.fill(HIST("TofTpcNsigma_after"), track1.tpcNSigmaPi(), track1.tofNSigmaPi(), track1.pt());
        if (cfgUseItsPID)
          histos.fill(HIST("TofItsNsigma_after"), itsResponse.nSigmaITS<o2::track::PID::Pion>(track1), track1.tofNSigmaPi(), track1.pt());
        break;
      case 2: // For Kaons
        if (!cfgUseItsPID)
          histos.fill(HIST("TofTpcNsigma_after"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
        if (cfgUseItsPID)
          histos.fill(HIST("TofItsNsigma_after"), itsResponse.nSigmaITS<o2::track::PID::Kaon>(track1), track1.tofNSigmaKa(), track1.pt());
        break;
      case 3: // For Protons
        if (!cfgUseItsPID)
          histos.fill(HIST("TofTpcNsigma_after"), track1.tpcNSigmaPr(), track1.tofNSigmaPr(), track1.pt());
        if (cfgUseItsPID)
          histos.fill(HIST("TofItsNsigma_after"), itsResponse.nSigmaITS<o2::track::PID::Proton>(track1), track1.tofNSigmaPr(), track1.pt());
        break;
    } // end of switch
  }

  float getDPhiStar(float charge1, float charge2, float phi1, float phi2, float pt1, float pt2, float radius, int magField)
  {
    int fbSign = (magField > 0) ? 1 : -1;

    float dPhiStar = phi1 - phi2 - charge1 * fbSign * std::asin(0.075 * radius / pt1) + charge2 * fbSign * std::asin(0.075 * radius / pt2);

    if (dPhiStar > constants::math::PI)
      dPhiStar = constants::math::TwoPI - dPhiStar;
    if (dPhiStar < -constants::math::PI)
      dPhiStar = -constants::math::TwoPI - dPhiStar;

    return dPhiStar;
  }

  template <CorrelationContainer::CFStep step, typename TTracks, typename TTracksAssoc>
  void fillCorrelations(TTracks tracks1, TTracksAssoc tracks2, float posZ, int system, int magneticField, float cent, float eventWeight) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {
    // Cache efficiency for particles (too many FindBin lookups)
    if (mEfficiency) {
      efficiencyAssociatedCache.clear();
      efficiencyAssociatedCache.reserve(tracks2.size());
      for (const auto& track2 : tracks2) {
        float weff = 1.;
        getEfficiencyCorrection(weff, track2.eta(), track2.pt(), posZ);
        efficiencyAssociatedCache.push_back(weff);
      }
    }

    if (system == SameEvent) {
      if (!cfgCentTableUnavailable)
        histos.fill(HIST("Centrality_used"), cent);
      histos.fill(HIST("Nch_used"), tracks1.size());
    }

    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);

    float triggerWeight = 1.0f;
    float associatedWeight = 1.0f;
    // loop over all tracks
    for (auto const& track1 : tracks1) {

      if (!trackSelected(track1))
        continue;

      // Fill Nsigma QA
      if (cfgGetNsigmaQA && !cfgUseItsPID) {
        if (cfgPIDParticle == kPions) {
          histos.fill(HIST("TofTpcNsigma_before"), track1.tpcNSigmaPi(), track1.tofNSigmaPi(), track1.pt());
          if (cfgGetdEdx) {
            double tpcExpSignalPi = track1.tpcSignal() - (track1.tpcNSigmaPi() * track1.tpcExpSigmaPi());

            histos.fill(HIST("TpcdEdx_ptwise"), track1.pt(), track1.tpcSignal(), track1.tofNSigmaPi());
            histos.fill(HIST("ExpTpcdEdx_ptwise"), track1.pt(), tpcExpSignalPi, track1.tofNSigmaPi());
            histos.fill(HIST("ExpSigma_ptwise"), track1.pt(), track1.tpcExpSigmaPi(), track1.tofNSigmaPi());
          }
        }
        if (cfgPIDParticle == kKaons) {
          histos.fill(HIST("TofTpcNsigma_before"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
          if (cfgGetdEdx) {
            double tpcExpSignalKa = track1.tpcSignal() - (track1.tpcNSigmaKa() * track1.tpcExpSigmaKa());

            histos.fill(HIST("TpcdEdx_ptwise"), track1.pt(), track1.tpcSignal(), track1.tofNSigmaKa());
            histos.fill(HIST("ExpTpcdEdx_ptwise"), track1.pt(), tpcExpSignalKa, track1.tofNSigmaKa());
            histos.fill(HIST("ExpSigma_ptwise"), track1.pt(), track1.tpcExpSigmaKa(), track1.tofNSigmaKa());
          }
        }
        if (cfgPIDParticle == kProtons) {
          histos.fill(HIST("TofTpcNsigma_before"), track1.tpcNSigmaPr(), track1.tofNSigmaPr(), track1.pt());
          if (cfgGetdEdx) {
            double tpcExpSignalPr = track1.tpcSignal() - (track1.tpcNSigmaPr() * track1.tpcExpSigmaPr());

            histos.fill(HIST("TpcdEdx_ptwise"), track1.pt(), track1.tpcSignal(), track1.tofNSigmaPr());
            histos.fill(HIST("ExpTpcdEdx_ptwise"), track1.pt(), tpcExpSignalPr, track1.tofNSigmaPr());
            histos.fill(HIST("ExpSigma_ptwise"), track1.pt(), track1.tpcExpSigmaPr(), track1.tofNSigmaPr());
          }
        }
        
      }
      if (cfgGetNsigmaQA && cfgUseItsPID) {
        if (cfgPIDParticle == kPions)
          histos.fill(HIST("TofItsNsigma_before"), itsResponse.nSigmaITS<o2::track::PID::Pion>(track1), track1.tofNSigmaPi(), track1.pt());
        if (cfgPIDParticle == kKaons)
          histos.fill(HIST("TofItsNsigma_before"), itsResponse.nSigmaITS<o2::track::PID::Kaon>(track1), track1.tofNSigmaKa(), track1.pt());
        if (cfgPIDParticle == kProtons)
          histos.fill(HIST("TofItsNsigma_before"), itsResponse.nSigmaITS<o2::track::PID::Proton>(track1), track1.tofNSigmaPr(), track1.pt());
      }

      if (cfgPIDParticle && getNsigmaPID(track1) != cfgPIDParticle)
        continue; // if PID is selected, check if the track has the right PID

      if (cfgGetNsigmaQA)
        fillNsigmaAfterCut(track1, cfgPIDParticle);

      if (!getEfficiencyCorrection(triggerWeight, track1.eta(), track1.pt(), posZ))
        continue;
      if (system == SameEvent) {
        histos.fill(HIST("Trig_hist"), fSampleIndex, posZ, track1.pt(), eventWeight * triggerWeight);
      }

      for (auto const& track2 : tracks2) {

        if (!trackSelected(track2))
          continue;
        if (mEfficiency) {
          associatedWeight = efficiencyAssociatedCache[track2.filteredIndex()];
        }

        if (!cfgUsePtOrder && track1.globalIndex() == track2.globalIndex())
          continue; // For pt-differential correlations, skip if the trigger and associate are the same track
        if (cfgUsePtOrder && system == SameEvent && track1.pt() <= track2.pt())
          continue; // Without pt-differential correlations, skip if the trigger pt is less than the associate pt
        if (cfgUsePtOrder && system == MixedEvent && cfgUsePtOrderInMixEvent && track1.pt() <= track2.pt())
          continue; // For pt-differential correlations in mixed events, skip if the trigger pt is less than the associate pt

        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -PIHalf);
        float deltaEta = track1.eta() - track2.eta();

        if (std::abs(deltaEta) < cfgCutMerging) {

          double dPhiStarHigh = getDPhiStar(track1.sign(), track2.sign(), track1.phi(), track2.phi(), track1.pt(), track2.pt(), cfgRadiusHigh, magneticField);
          double dPhiStarLow = getDPhiStar(track1.sign(), track2.sign(), track1.phi(), track2.phi(), track1.pt(), track2.pt(), cfgRadiusLow, magneticField);

          const double kLimit = 3.0 * cfgCutMerging;

          bool bIsBelow = false;

          if (std::abs(dPhiStarLow) < kLimit || std::abs(dPhiStarHigh) < kLimit || dPhiStarLow * dPhiStarHigh < 0) {
            for (double rad(cfgRadiusLow); rad < cfgRadiusHigh; rad += 0.01) {
              double dPhiStar = getDPhiStar(track1.sign(), track2.sign(), track1.phi(), track2.phi(), track1.pt(), track2.pt(), rad, magneticField);
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

          same->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track2.pt(), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          histos.fill(HIST("deltaEta_deltaPhi_same"), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
        } else if (system == MixedEvent) {

          mixed->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track2.pt(), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          histos.fill(HIST("deltaEta_deltaPhi_mixed"), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
        }
      }
    }
  }

  template <CorrelationContainer::CFStep step, typename TV0Tracks, typename TTracksAssoc>
  void fillCorrelationsReso(TV0Tracks tracks1, TTracksAssoc tracks2, float posZ, float posY, float posX, int system, int magneticField, float cent, float eventWeight) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {
    // Cache efficiency for particles (too many FindBin lookups)
    if (mEfficiency) {
      efficiencyAssociatedCache.clear();
      efficiencyAssociatedCache.reserve(tracks2.size());
      for (const auto& track2 : tracks2) {
        float weff = 1.;
        getEfficiencyCorrection(weff, track2.eta(), track2.pt(), posZ);
        efficiencyAssociatedCache.push_back(weff);
      }
    }

    if (system == SameEvent) {
      if (!cfgCentTableUnavailable)
        histos.fill(HIST("Centrality_used"), cent);
      histos.fill(HIST("Nch_used"), tracks2.size()); // Taking Nch from tracks2 since tracks1 are V0s
    }

    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);

    float triggerWeight = 1.0f;
    float associatedWeight = 1.0f;
    // loop over all tracks
    for (auto const& track1 : tracks1) {

      double resoMass = -1;

      // 4 = kshort, 5 = lambda, 6 = phi
      if (cfgPIDParticle == kK0) {
        if (!isSelectedK0(track1, posZ, posY, posX))
          continue; // Reject if called for K0 but V0 is not K0

        resoMass = track1.mK0Short();
      }

      if (cfgPIDParticle == kLambda) {
        if (!isSelectedLambda(track1, posZ, posY, posX))
          continue; // Reject if called for Lambda but V0 is not lambda

        resoMass = track1.mLambda();
      }

      if (system == SameEvent) {
        histos.fill(HIST("Trig_histReso"), fSampleIndex, posZ, track1.pt(), resoMass, eventWeight * triggerWeight);
      }

      for (auto const& track2 : tracks2) {

        if (!trackSelected(track2))
          continue;
        if (track2.pt() < cfgCutPtMin || track2.pt() > cfgCutPtMax) // Select associated particles in the pt range 0.2 - 3.0 GeV/c
          continue;
        if (mEfficiency) {
          associatedWeight = efficiencyAssociatedCache[track2.filteredIndex()];
        }

        if (cfgUsePtOrder && system == SameEvent && track1.pt() <= track2.pt())
          continue; // Without pt-differential correlations, skip if the trigger pt is less than the associate pt
        if (cfgUsePtOrder && system == MixedEvent && cfgUsePtOrderInMixEvent && track1.pt() <= track2.pt())
          continue; // For pt-differential correlations in mixed events, skip if the trigger pt is less than the associate pt

        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -PIHalf);
        float deltaEta = track1.eta() - track2.eta();

        if (std::abs(deltaEta) < cfgCutMerging) {

          double dPhiStarHigh = getDPhiStar(0, track2.sign(), track1.phi(), track2.phi(), track1.pt(), track2.pt(), cfgRadiusHigh, magneticField);
          double dPhiStarLow = getDPhiStar(0, track2.sign(), track1.phi(), track2.phi(), track1.pt(), track2.pt(), cfgRadiusLow, magneticField);

          const double kLimit = 3.0 * cfgCutMerging;

          bool bIsBelow = false;

          if (std::abs(dPhiStarLow) < kLimit || std::abs(dPhiStarHigh) < kLimit || dPhiStarLow * dPhiStarHigh < 0) {
            for (double rad(cfgRadiusLow); rad < cfgRadiusHigh; rad += 0.01) {
              double dPhiStar = getDPhiStar(0, track2.sign(), track1.phi(), track2.phi(), track1.pt(), track2.pt(), rad, magneticField);
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
          sameReso->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), resoMass, deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          histos.fill(HIST("deltaEta_deltaPhi_same"), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
        } else if (system == MixedEvent) {
          mixedReso->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), resoMass, deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          histos.fill(HIST("deltaEta_deltaPhi_mixed"), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
        }
      }
    }
  }

  template <typename TCollision>
  bool selectionEvent(TCollision collision, const int mult, const float cent, const bool fillCounter)
  {
    if (fillCounter)
      histos.fill(HIST("hEventCount"), kFilteredEvents);
    if (!collision.sel8()) {
      return 0;
    }
    if (fillCounter)
      histos.fill(HIST("hEventCount"), kAfterSel8);

    if (eventCuts[kUseNoTimeFrameBorder][kEvCut1] && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return 0;
    }
    if (fillCounter && eventCuts[kUseNoTimeFrameBorder][kEvCut1])
      histos.fill(HIST("hEventCount"), kUseNoTimeFrameBorder);

    if (eventCuts[kUseNoITSROFrameBorder][kEvCut1] && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return 0;
    }
    if (fillCounter && eventCuts[kUseNoITSROFrameBorder][kEvCut1])
      histos.fill(HIST("hEventCount"), kUseNoITSROFrameBorder);

    if (eventCuts[kUseNoSameBunchPileup][kEvCut1] && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return 0;
    }
    if (fillCounter && eventCuts[kUseNoSameBunchPileup][kEvCut1])
      histos.fill(HIST("hEventCount"), kUseNoSameBunchPileup);

    if (eventCuts[kUseGoodZvtxFT0vsPV][kEvCut1] && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return 0;
    }
    if (fillCounter && eventCuts[kUseGoodZvtxFT0vsPV][kEvCut1])
      histos.fill(HIST("hEventCount"), kUseGoodZvtxFT0vsPV);

    if (eventCuts[kUseNoCollInTimeRangeStandard][kEvCut1] && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      // no collisions in specified time range
      return 0;
    }
    if (fillCounter && eventCuts[kUseNoCollInTimeRangeStandard][kEvCut1])
      histos.fill(HIST("hEventCount"), kUseNoCollInTimeRangeStandard);

    if (eventCuts[kUseGoodITSLayersAll][kEvCut1] && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      // from Jan 9 2025 AOT meeting
      // cut time intervals with dead ITS staves
      return 0;
    }
    if (fillCounter && eventCuts[kUseGoodITSLayersAll][kEvCut1])
      histos.fill(HIST("hEventCount"), kUseGoodITSLayersAll);

    if (eventCuts[kUseNoCollInRofStandard][kEvCut1] && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      // no other collisions in this Readout Frame with per-collision multiplicity above threshold
      return 0;
    }
    if (fillCounter && eventCuts[kUseNoCollInRofStandard][kEvCut1])
      histos.fill(HIST("hEventCount"), kUseNoCollInRofStandard);

    if (eventCuts[kUseNoHighMultCollInPrevRof][kEvCut1] && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      // veto an event if FT0C amplitude in previous ITS ROF is above threshold
      return 0;
    }
    if (fillCounter && eventCuts[kUseNoHighMultCollInPrevRof][kEvCut1])
      histos.fill(HIST("hEventCount"), kUseNoHighMultCollInPrevRof);

    auto multNTracksPV = collision.multNTracksPV();
    auto occupancy = collision.trackOccupancyInTimeRange();

    if (eventCuts[kUseOccupancy][kEvCut1] && (occupancy < cfgCutOccupancyMin || occupancy > cfgCutOccupancyMax)) {
      return 0;
    }
    if (fillCounter && eventCuts[kUseOccupancy][kEvCut1])
      histos.fill(HIST("hEventCount"), kUseOccupancy);

    if (eventCuts[kUseMultCorrCut][kEvCut1]) {
      if (cfgFuncParas.cfgMultPVT0CCutEnabled && !cfgCentTableUnavailable) {
        if (multNTracksPV < cfgFuncParas.fMultPVT0CCutLow->Eval(cent))
          return 0;
        if (multNTracksPV > cfgFuncParas.fMultPVT0CCutHigh->Eval(cent))
          return 0;
      }
      if (cfgFuncParas.cfgMultT0CCutEnabled && !cfgCentTableUnavailable) {
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

    if (fillCounter && eventCuts[kUseMultCorrCut][kEvCut1])
      histos.fill(HIST("hEventCount"), kUseMultCorrCut);

    // V0A T0A 5 sigma cut
    if (eventCuts[kUseT0AV0ACut][kEvCut1] && (std::fabs(collision.multFV0A() - cfgFuncParas.fT0AV0AMean->Eval(collision.multFT0A())) > cfgV0AT0Acut * cfgFuncParas.fT0AV0ASigma->Eval(collision.multFT0A())))
      return 0;
    if (fillCounter && eventCuts[kUseT0AV0ACut][kEvCut1])
      histos.fill(HIST("hEventCount"), kUseT0AV0ACut);

    if (eventCuts[kUseVertexITSTPC][kEvCut1] && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
      return 0;
    if (fillCounter && eventCuts[kUseVertexITSTPC][kEvCut1])
      histos.fill(HIST("hEventCount"), kUseVertexITSTPC);

    if (eventCuts[kUseTVXinTRD][kEvCut1] && collision.alias_bit(kTVXinTRD)) {
      return 0;
    }
    if (fillCounter && eventCuts[kUseTVXinTRD][kEvCut1])
      histos.fill(HIST("hEventCount"), kUseTVXinTRD);

    return 1;
  }

  void processSame(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float cent = -1.;
    float weightCent = 1.0f;
    if (!cfgCentTableUnavailable)
      cent = getCentrality(collision);
    if (cfgUseAdditionalEventCut && !selectionEvent(collision, tracks.size(), cent, true))
      return;

    loadCorrection(bc.timestamp());
    if (!cfgCentTableUnavailable) {
      getCentralityWeight(weightCent, cent);
      histos.fill(HIST("Centrality"), cent);
      histos.fill(HIST("CentralityWeighted"), cent, weightCent);
    }
    histos.fill(HIST("Nch"), tracks.size());
    histos.fill(HIST("zVtx"), collision.posZ());

    if (cfgSelCollByNch && (tracks.size() < cfgCutMultMin || tracks.size() >= cfgCutMultMax)) {
      return;
    }
    if (!cfgSelCollByNch && !cfgCentTableUnavailable && (cent < cfgCutCentMin || cent >= cfgCutCentMax)) {
      return;
    }
    histos.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin
    fillYield(collision, tracks);

    same->fillEvent(tracks.size(), CorrelationContainer::kCFStepReconstructed);
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(tracks, tracks, collision.posZ(), SameEvent, getMagneticField(bc.timestamp()), cent, weightCent);
  }
  PROCESS_SWITCH(PidDiHadron, processSame, "Process same event", true);

  // the process for filling the mixed events
  void processMixed(FilteredCollisions const& collisions, FilteredTracks const& tracks, aod::BCsWithTimestamps const&)
  {

    auto getTracksSize = [&tracks, this](FilteredCollisions::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}, true};

    auto tracksTuple = std::make_tuple(tracks, tracks);
    Pair<FilteredCollisions, FilteredTracks, FilteredTracks, MixedBinning> pairs{binningOnVtxAndMult, cfgMixEventNumMin, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;
      if (!collision1.sel8() || !collision2.sel8())
        continue;

      if (cfgSelCollByNch && (tracks1.size() < cfgCutMultMin || tracks1.size() >= cfgCutMultMax))
        continue;

      if (cfgSelCollByNch && (tracks2.size() < cfgCutMultMin || tracks2.size() >= cfgCutMultMax))
        continue;

      float cent1 = -1;
      float cent2 = -1;
      if (!cfgCentTableUnavailable) {
        cent1 = getCentrality(collision1);
        cent2 = getCentrality(collision2);
      }
      if (cfgUseAdditionalEventCut && !selectionEvent(collision1, tracks1.size(), cent1, false))
        continue;
      if (cfgUseAdditionalEventCut && !selectionEvent(collision2, tracks2.size(), cent2, false))
        continue;

      if (!cfgSelCollByNch && !cfgCentTableUnavailable && (cent1 < cfgCutCentMin || cent1 >= cfgCutCentMax))
        continue;

      if (!cfgSelCollByNch && !cfgCentTableUnavailable && (cent2 < cfgCutCentMin || cent2 >= cfgCutCentMax))
        continue;

      histos.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin
      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      loadCorrection(bc.timestamp());
      float eventWeight = 1.0f;
      if (cfgUseEventWeights) {
        eventWeight = 1.0f / it.currentWindowNeighbours();
      }
      float weightCent = 1.0f;
      if (!cfgCentTableUnavailable)
        getCentralityWeight(weightCent, cent1);

      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(tracks1, tracks2, collision1.posZ(), MixedEvent, getMagneticField(bc.timestamp()), cent1, eventWeight * weightCent);
    }
  }
  PROCESS_SWITCH(PidDiHadron, processMixed, "Process mixed events", true);

  template <typename V0>
  bool isSelectedK0(V0 const& candidate, float posZ, float posY, float posX)
  {
    double mk0 = candidate.mK0Short();

    // separate the positive and negative V0 daughters
    auto postrack = candidate.template posTrack_as<FilteredTracks>();
    auto negtrack = candidate.template negTrack_as<FilteredTracks>();

    histos.fill(HIST("hK0Count"), 0.5);
    if (postrack.pt() < resoCutVals[kPosTrackPt][iK0] || negtrack.pt() < resoCutVals[kNegTrackPt][iK0])
      return false;
    histos.fill(HIST("hK0Count"), 1.5);
    if (mk0 < resoCutVals[kMassMin][iK0] && mk0 > resoCutVals[kMassMax][iK0])
      return false;
    histos.fill(HIST("hK0Count"), 2.5);
    // Rapidity correction
    if (candidate.yK0Short() > resoCutVals[kRapidity][iK0])
      return false;
    histos.fill(HIST("hK0Count"), 3.5);
    // DCA cuts for K0short
    if (std::abs(candidate.dcapostopv()) < resoCutVals[kDCAPosToPVMin][iK0] || std::abs(candidate.dcanegtopv()) < resoCutVals[kDCANegToPVMin][iK0])
      return false;
    histos.fill(HIST("hK0Count"), 4.5);
    if (std::abs(candidate.dcaV0daughters()) > resoSwitchVals[kDCABetDaug][iK0])
      return false;
    histos.fill(HIST("hK0Count"), 5.5);
    // v0 radius cuts
    if (resoSwitchVals[kUseV0Radius][iK0] && (candidate.v0radius() < resoCutVals[kRadiusMin][iK0] || candidate.v0radius() > resoCutVals[kRadiusMax][iK0]))
      return false;
    histos.fill(HIST("hK0Count"), 6.5);
    // cosine pointing angle cuts
    if (candidate.v0cosPA() < resoCutVals[kCosPA][iK0])
      return false;
    histos.fill(HIST("hK0Count"), 7.5);
    // Proper lifetime
    float cTauK0 = candidate.distovertotmom(posX, posY, posZ) * massK0Short;
    if (resoSwitchVals[kUseProperLifetime][iK0] && std::abs(cTauK0) > resoCutVals[kLifeTime][iK0])
      return false;
    histos.fill(HIST("hK0Count"), 8.5);
    // ArmenterosPodolanskiCut
    if (resoSwitchVals[kUseArmPodCut][iK0] && (candidate.qtarm() / std::abs(candidate.alpha())) < resoCutVals[kArmPodMinVal][iK0])
      return false;
    histos.fill(HIST("hK0Count"), 9.5);
    // Selection on V0 daughters
    if (!selectionV0Daughter(postrack, kPions) || !selectionV0Daughter(negtrack, kPions))
      return false;
    histos.fill(HIST("hK0Count"), 10.5);

    histos.fill(HIST("hK0Phi"), candidate.phi());
    histos.fill(HIST("hK0Eta"), candidate.eta());
    histos.fill(HIST("PiPlusTPC_K0"), postrack.pt(), postrack.tpcNSigmaPi());
    histos.fill(HIST("PiPlusTOF_K0"), postrack.pt(), postrack.tofNSigmaPi());
    histos.fill(HIST("PiMinusTPC_K0"), negtrack.pt(), negtrack.tpcNSigmaPi());
    histos.fill(HIST("PiMinusTOF_K0"), negtrack.pt(), negtrack.tofNSigmaPi());

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

    histos.fill(HIST("hLambdaCount"), 0.5);
    if (postrack.pt() < resoCutVals[kPosTrackPt][iLambda] || negtrack.pt() < resoCutVals[kNegTrackPt][iLambda])
      return false;

    histos.fill(HIST("hLambdaCount"), 1.5);
    if (mlambda > resoCutVals[kMassMin][iLambda] && mlambda < resoCutVals[kMassMax][iLambda])
      isL = true;
    if (mantilambda > resoCutVals[kMassMin][iLambda] && mantilambda < resoCutVals[kMassMax][iLambda])
      isAL = true;

    if (!isL && !isAL) {
      return false;
    }
    histos.fill(HIST("hLambdaCount"), 2.5);

    // Rapidity correction
    if (candidate.yLambda() > resoCutVals[kRapidity][iLambda])
      return false;
    histos.fill(HIST("hLambdaCount"), 3.5);
    // DCA cuts for lambda and antilambda
    if (isL) {
      if (std::abs(candidate.dcapostopv()) < resoCutVals[kDCAPosToPVMin][iLambda] || std::abs(candidate.dcanegtopv()) < resoCutVals[kDCANegToPVMin][iLambda])
        return false;
    }
    if (isAL) {
      if (std::abs(candidate.dcapostopv()) < resoCutVals[kDCANegToPVMin][iLambda] || std::abs(candidate.dcanegtopv()) < resoCutVals[kDCAPosToPVMin][iLambda])
        return false;
    }
    histos.fill(HIST("hLambdaCount"), 4.5);
    if (std::abs(candidate.dcaV0daughters()) > resoSwitchVals[kDCABetDaug][iLambda])
      return false;
    histos.fill(HIST("hLambdaCount"), 5.5);
    // v0 radius cuts
    if (resoSwitchVals[kUseV0Radius][iLambda] && (candidate.v0radius() < resoCutVals[kRadiusMin][iLambda] || candidate.v0radius() > resoCutVals[kRadiusMax][iLambda]))
      return false;
    histos.fill(HIST("hLambdaCount"), 6.5);
    // cosine pointing angle cuts
    if (candidate.v0cosPA() < resoCutVals[kCosPA][iLambda])
      return false;
    histos.fill(HIST("hLambdaCount"), 7.5);
    // Proper lifetime
    float cTauLambda = candidate.distovertotmom(posX, posY, posZ) * massLambda;
    if (resoSwitchVals[kUseProperLifetime][iLambda] && cTauLambda > resoCutVals[kLifeTime][iLambda])
      return false;
    histos.fill(HIST("hLambdaCount"), 8.5);
    if (isL) {
      if (!selectionV0Daughter(postrack, kProtons) || !selectionV0Daughter(negtrack, kPions))
        return false;
    }
    if (isAL) {
      if (!selectionV0Daughter(postrack, kPions) || !selectionV0Daughter(negtrack, kProtons))
        return false;
    }
    histos.fill(HIST("hLambdaCount"), 9.5);

    if (!cfgUseAntiLambda && isAL) { // Reject the track if it is antilambda
      return false;
    }

    histos.fill(HIST("hLambdaPhi"), candidate.phi());
    histos.fill(HIST("hLambdaEta"), candidate.eta());
    histos.fill(HIST("PrPlusTPC_L"), postrack.pt(), postrack.tpcNSigmaPr());
    histos.fill(HIST("PrPlusTOF_L"), postrack.pt(), postrack.tofNSigmaPr());
    histos.fill(HIST("PiMinusTPC_L"), negtrack.pt(), negtrack.tpcNSigmaPi());
    histos.fill(HIST("PiMinusTOF_L"), negtrack.pt(), negtrack.tofNSigmaPi());

    return true;
  }

  double massKaPlus = o2::constants::physics::MassKPlus;    // same as MassKMinus
  double massLambda = o2::constants::physics::MassLambda;   // same as MassLambda0 and MassLambda0Bar
  double massK0Short = o2::constants::physics::MassK0Short; // same as o2::constants::physics::MassK0 and o2::constants::physics::MassK0Bar

  void processSameReso(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::BCsWithTimestamps const&, aod::V0Datas const& V0s)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float cent = -1.;
    float weightCent = 1.0f;
    if (!cfgCentTableUnavailable)
      cent = getCentrality(collision);
    if (cfgUseAdditionalEventCut && !selectionEvent(collision, tracks.size(), cent, true))
      return;

    loadCorrection(bc.timestamp());
    if (!cfgCentTableUnavailable) {
      getCentralityWeight(weightCent, cent);
      histos.fill(HIST("Centrality"), cent);
      histos.fill(HIST("CentralityWeighted"), cent, weightCent);
    }
    histos.fill(HIST("Nch"), tracks.size());
    histos.fill(HIST("zVtx"), collision.posZ());

    if (cfgSelCollByNch && (tracks.size() < cfgCutMultMin || tracks.size() >= cfgCutMultMax)) {
      return;
    }
    if (!cfgSelCollByNch && !cfgCentTableUnavailable && (cent < cfgCutCentMin || cent >= cfgCutCentMax)) {
      return;
    }

    histos.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin

    sameReso->fillEvent(tracks.size(), CorrelationContainer::kCFStepReconstructed);
    fillCorrelationsReso<CorrelationContainer::kCFStepReconstructed>(V0s, tracks, collision.posZ(), collision.posY(), collision.posX(), SameEvent, getMagneticField(bc.timestamp()), cent, weightCent);
  }
  PROCESS_SWITCH(PidDiHadron, processSameReso, "Process same event for resonances", true);

  void processMixedReso(FilteredCollisions const& collisions, FilteredTracks const& tracks, aod::BCsWithTimestamps const&, aod::V0Datas const& V0s)
  {
    auto getTracksSize = [&tracks, this](FilteredCollisions::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}, true};

    auto tracksTuple = std::make_tuple(V0s, tracks);
    Pair<FilteredCollisions, aod::V0Datas, FilteredTracks, MixedBinning> pairs{binningOnVtxAndMult, cfgMixEventNumMin, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, v0s1, collision2, tracks2] = *it;
      if (!collision1.sel8() || !collision2.sel8())
        continue;

      auto tracks1 = tracks.sliceByCached(o2::aod::track::collisionId, collision1.globalIndex(), this->cache);
      if (cfgSelCollByNch && (tracks1.size() < cfgCutMultMin || tracks1.size() >= cfgCutMultMax))
        continue;

      if (cfgSelCollByNch && (tracks2.size() < cfgCutMultMin || tracks2.size() >= cfgCutMultMax))
        continue;

      float cent1 = -1;
      float cent2 = -1;
      if (!cfgCentTableUnavailable) {
        cent1 = getCentrality(collision1);
        cent2 = getCentrality(collision2);
      }
      if (cfgUseAdditionalEventCut && !selectionEvent(collision1, tracks1.size(), cent1, false))
        continue;
      if (cfgUseAdditionalEventCut && !selectionEvent(collision2, tracks2.size(), cent2, false))
        continue;

      if (!cfgSelCollByNch && !cfgCentTableUnavailable && (cent1 < cfgCutCentMin || cent1 >= cfgCutCentMax))
        continue;

      if (!cfgSelCollByNch && !cfgCentTableUnavailable && (cent2 < cfgCutCentMin || cent2 >= cfgCutCentMax))
        continue;

      histos.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin
      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      loadCorrection(bc.timestamp());
      float eventWeight = 1.0f;
      if (cfgUseEventWeights) {
        eventWeight = 1.0f / it.currentWindowNeighbours();
      }
      float weightCent = 1.0f;
      if (!cfgCentTableUnavailable)
        getCentralityWeight(weightCent, cent1);

      fillCorrelationsReso<CorrelationContainer::kCFStepReconstructed>(v0s1, tracks2, collision1.posZ(), collision1.posY(), collision1.posX(), MixedEvent, getMagneticField(bc.timestamp()), cent1, eventWeight * weightCent);
    }
  }
  PROCESS_SWITCH(PidDiHadron, processMixedReso, "Process mixed events", true);

  void processMC(FilteredCollisions::iterator const& collision, FilteredMcTracks const& tracksmc, aod::BCsWithTimestamps const&, aod::McParticles const&)
  {
    float cent = -1.;
    if (!cfgCentTableUnavailable)
      cent = getCentrality(collision);
    if (cfgUseAdditionalEventCut && !selectionEvent(collision, tracksmc.size(), cent, true))
      return;

    if (cfgSelCollByNch && (tracksmc.size() < cfgCutMultMin || tracksmc.size() >= cfgCutMultMax)) {
      return;
    }
    if (!cfgSelCollByNch && !cfgCentTableUnavailable && (cent < cfgCutCentMin || cent >= cfgCutCentMax)) {
      return;
    }

    // loop over all tracks
    for (auto const& track : tracksmc) {

      if (!trackSelected(track))
        continue;

      int pidIndex = getNsigmaPID(track);

      // Fill Counts for selection through Nsigma cuts
      if (pidIndex == kPions)
        histos.fill(HIST("hNsigmaPionSelected"), track.pt());
      if (pidIndex == kKaons)
        histos.fill(HIST("hNsigmaKaonSelected"), track.pt());
      if (pidIndex == kProtons)
        histos.fill(HIST("hNsigmaProtonSelected"), track.pt());

      // Check the PDG code for the particles (MC truth) and match with analysed Nsigma PID
      if (std::abs(track.mcParticle().pdgCode()) == PDG_t::kPiPlus) {
        histos.fill(HIST("hNsigmaPionTrue"), track.pt());
        histos.fill(HIST("TpcdEdx_ptwise_pi"), track.pt(), track.tpcSignal(), track.tofNSigmaPi());

        double tpcExpSignalPi = track.tpcSignal() - (track.tpcNSigmaPi() * track.tpcExpSigmaPi());
        histos.fill(HIST("ExpTpcdEdx_ptwise_pi"), track.pt(), tpcExpSignalPi, track.tofNSigmaPi());

        if (pidIndex == kPions) {
          histos.fill(HIST("hNsigmaPionTruePositives"), track.pt());
        }
      } // Pion condition

      if (std::abs(track.mcParticle().pdgCode()) == PDG_t::kKPlus) {
        histos.fill(HIST("hNsigmaKaonTrue"), track.pt());
        histos.fill(HIST("TpcdEdx_ptwise_ka"), track.pt(), track.tpcSignal(), track.tofNSigmaKa());

        double tpcExpSignalKa = track.tpcSignal() - (track.tpcNSigmaKa() * track.tpcExpSigmaKa());
        histos.fill(HIST("ExpTpcdEdx_ptwise_ka"), track.pt(), tpcExpSignalKa, track.tofNSigmaKa());

        if (pidIndex == kKaons) {
          histos.fill(HIST("hNsigmaKaonTruePositives"), track.pt());
        }
      } // Kaon condition

      if (std::abs(track.mcParticle().pdgCode()) == PDG_t::kProton) {
        histos.fill(HIST("hNsigmaProtonTrue"), track.pt());
        histos.fill(HIST("TpcdEdx_ptwise_pr"), track.pt(), track.tpcSignal(), track.tofNSigmaPr());

        double tpcExpSignalPr = track.tpcSignal() - (track.tpcNSigmaPr() * track.tpcExpSigmaPr());
        histos.fill(HIST("ExpTpcdEdx_ptwise_pr"), track.pt(), tpcExpSignalPr, track.tofNSigmaPr());

        if (pidIndex == kProtons) {
          histos.fill(HIST("hNsigmaProtonTruePositives"), track.pt());
        }
      } // Proton condition

    } // end of tracks MC loop
  } // end of process MC
  PROCESS_SWITCH(PidDiHadron, processMC, "Process Monte Carlo", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PidDiHadron>(cfgc),
  };
}
