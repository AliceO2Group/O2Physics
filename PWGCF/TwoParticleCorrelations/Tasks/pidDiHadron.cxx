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

#include "PWGCF/TwoParticleCorrelations/Core/DihadronContainer.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
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
#include <DataFormatsParameters/GRPMagField.h>
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
#include <TPDGCode.h>
#include <TRandom3.h>
#include <TString.h>

#include <RtypesCore.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;

// define the filtered collisions and tracks
#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

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

struct PidDiHadron {
  o2::aod::ITSResponse itsResponse;
  Service<ccdb::BasicCCDBManager> ccdb;

  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgQABasic, bool, true, "Enable QA histograms for event and track selection")
    O2_DEFINE_CONFIGURABLE(cfgStrictTrackCounter, bool, false, "Strict track counter for multiplicity correlation cut, counts only tracks that pass all cuts and are used in the correlation")
    O2_DEFINE_CONFIGURABLE(cfgSampleSize, int, 10, "Sample size for mixed event")
    O2_DEFINE_CONFIGURABLE(cfgAsscPtMin, float, 0.2f, "Minimal pT for ref tracks")
    O2_DEFINE_CONFIGURABLE(cfgAsscPtMax, float, 3.0f, "Maximal pT for ref tracks")
    O2_DEFINE_CONFIGURABLE(cfgCutMerging, float, 0.0, "Merging cut on track merge")
    O2_DEFINE_CONFIGURABLE(cfgUseCutMerging, bool, false, "use two track cut");
    O2_DEFINE_CONFIGURABLE(cfgMixEventNumMin, int, 5, "Minimum number of events to mix")
    O2_DEFINE_CONFIGURABLE(cfgRadiusLow, float, 0.8, "Low radius for merging cut")
    O2_DEFINE_CONFIGURABLE(cfgRadiusHigh, float, 2.5, "High radius for merging cut")
    O2_DEFINE_CONFIGURABLE(cfgUsePtOrder, bool, false, "enable trigger pT < associated pT cut")
    O2_DEFINE_CONFIGURABLE(cfgUsePtOrderInMixEvent, bool, false, "enable trigger pT < associated pT cut in mixed event")
  } cfgGeneral;

  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
    O2_DEFINE_CONFIGURABLE(cfgEfficiencyNch, std::string, "", "CCDB path to multiplicity dependent efficiency object")
    O2_DEFINE_CONFIGURABLE(cfgUseEventWeights, bool, false, "Use event weights for mixed event")
    O2_DEFINE_CONFIGURABLE(cfgCentralityWeight, std::string, "", "CCDB path to centrality weight object")
    O2_DEFINE_CONFIGURABLE(cfgLocalEfficiency, bool, false, "Use local efficiency object")
    O2_DEFINE_CONFIGURABLE(cfgLocalEfficiencyNch, bool, false, "Use local multiplicity dependent efficiency object");
    O2_DEFINE_CONFIGURABLE(cfgUseTrigWeight, bool, false, "enable trigger efficiency correction")
  } cfgEffWeights;

  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgGlobalCutVertex, float, 10.0f, "Accepted z-vertex range")
    O2_DEFINE_CONFIGURABLE(cfgCutOccupancyMax, int, 2000, "Maximum occupancy cut")
    O2_DEFINE_CONFIGURABLE(cfgCutOccupancyMin, int, 0, "Minimum occupancy cut")
    O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, true, "Use additional event cut on mult correlations")
    O2_DEFINE_CONFIGURABLE(cfgCutMultMin, int, 0, "Minimum multiplicity for collision")
    O2_DEFINE_CONFIGURABLE(cfgCutMultMax, int, 10, "Maximum multiplicity for collision")
    O2_DEFINE_CONFIGURABLE(cfgSelCollByNch, bool, true, "Select collisions by Nch or centrality")
    O2_DEFINE_CONFIGURABLE(cfgCentEstimator, int, 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FT0A")
    O2_DEFINE_CONFIGURABLE(cfgCentTableUnavailable, bool, false, "if a dataset does not provide centrality information")
    O2_DEFINE_CONFIGURABLE(cfgCutCentMin, float, 60.0f, "Minimum centrality for collision")
    O2_DEFINE_CONFIGURABLE(cfgCutCentMax, float, 80.0f, "Maximum centrality for collision")
    O2_DEFINE_CONFIGURABLE(cfgUseCol2ForMixedEvSel, bool, true, "Use Mult/Cent cut on collision 2 during mixed event");
    Configurable<LabeledArray<int>> cfgUseEventCuts{"cfgUseEventCuts", {LongArrayInt.front().data(), 14, 1, {"Filtered Events", "Sel8", "kNoTimeFrameBorder", "kNoITSROFrameBorder", "kNoSameBunchPileup", "kIsGoodZvtxFT0vsPV", "kNoCollInTimeRangeStandard", "kIsGoodITSLayersAll", "kIsGoodITSLayer0123", "kNoCollInRofStandard", "kNoHighMultCollInPrevRof", "Occupancy", "Multcorrelation", "T0AV0ACut"}, {"EvCuts"}}, "Labeled array (int) for various cuts on resonances"};
  } cfgEvSel;

  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgGlobalPtMin, float, 0.2f, "minimum accepted track pT")
    O2_DEFINE_CONFIGURABLE(cfgGlobalPtMax, float, 10.0f, "maximum accepted track pT")
    O2_DEFINE_CONFIGURABLE(cfgGlobalEta, float, 0.8f, "Eta cut")
    O2_DEFINE_CONFIGURABLE(cfgGlobalChi2prTPCcls, float, 2.5f, "Maximum chi2 per TPC clusters")
    O2_DEFINE_CONFIGURABLE(cfgGlobalDCAz, float, 2.0f, "Maximum DCA to vertex z")
    Configurable<std::vector<std::string>> cfgTrackCutsDCAxy{"cfgTrackCutsDCAxy", std::vector<std::string>{"[0]*(0.0026+0.005/(x^1.01))", "[0]*(0.0026+0.005/(x^1.01))"}, "Functional form of pt-dependent DCAxy cut"};
    Configurable<std::vector<std::string>> cfgTrackCutsDCAz{"cfgTrackCutsDCAz", std::vector<std::string>{"", ""}, "Functional form of pt-dependent DCAz cut"};
    Configurable<LabeledArray<float>> cfgTrackCuts{"cfgTrackCuts", {LongArrayFloat.front().data(), 10, 2, {"PtMin", "PtMax", "EtaCut", "Chi2PrTpcCls", "TpcCluster", "TpcCrossedRows", "ItsCluster", "DCAz", "DCAzNsigma", "DCAxyNsigma"}, {"General", "SelNch"}}, "Labeled array for track selections"};
  } cfgTrkSel;

  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgQAdEdx, bool, false, "Get dEdx histograms for pions, kaons, and protons")
    O2_DEFINE_CONFIGURABLE(cfgQANsigma, bool, false, "Get QA histograms for selection of pions, kaons, and protons")
    O2_DEFINE_CONFIGURABLE(cfgPIDTofPtCut, float, 0.4f, "Minimum pt to use TOF N-sigma")
    O2_DEFINE_CONFIGURABLE(cfgPIDUseIts, bool, false, "Use ITS PID for particle identification")
    O2_DEFINE_CONFIGURABLE(cfgPIDUseRejection, bool, true, "True: use exclusion exclusion criteria for PID determination, false: don't use exclusion")
    O2_DEFINE_CONFIGURABLE(cfgPIDOnlyTpcCut, float, 3.0f, "TPC N-sigma cut for pions, kaons, protons")
    O2_DEFINE_CONFIGURABLE(cfgPIDUseOnlyTpc, bool, true, "Use only TPC PID for daughter selection")
    O2_DEFINE_CONFIGURABLE(cfgPIDParticle, int, 0, "1 = pion, 2 = kaon, 3 = proton, 4 = kshort, 5 = lambda, 0 for no PID")
    O2_DEFINE_CONFIGURABLE(cfgPIDUseReso, bool, false, "Switch to turn on resonance analysis")
    Configurable<LabeledArray<float>> cfgPIDnSigmas{"cfgPIDnSigmas", {LongArrayFloat.front().data(), 6, 3, {"UpCut_pi", "UpCut_ka", "UpCut_pr", "LowCut_pi", "LowCut_ka", "LowCut_pr"}, {"TPC", "TOF", "ITS"}}, "Labeled array for n-sigma values for TPC, TOF, ITS for pions, kaons, protons (positive and negative)"};
    Configurable<LabeledArray<float>> cfgPIDResoCuts{"cfgPIDResoCuts", {LongArrayFloat.front().data(), 15, 2, {"cos_PAs", "massMin", "massMax", "PosTrackPt", "NegTrackPt", "DCAPosToPVMin", "DCANegToPVMin", "Lifetime", "RadiusMin", "RadiusMax", "RapidityMin", "RapidityMax", "ArmPodMinVal", "DCABetDaug", "CompMassRej"}, {"K0", "Lambda"}}, "Labeled array (float) for various cuts on resonances"};
    Configurable<LabeledArray<int>> cfgPIDResoSwitches{"cfgPIDResoSwitches", {LongArrayInt.front().data(), 7, 2, {"UseCosPA", "NMassBins", "UseDCABetDaug", "UseProperLifetime", "UseV0Radius", "UseArmPodCut", "UseCompMassRej"}, {"K0", "Lambda"}}, "Labeled array (int) for various cuts on resonances"};
  } cfgPIDConfigs;

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

  SliceCache cache;

  // Axes for general QA
  ConfigurableAxis axisVertex{"axisVertex", {10, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {100, 0, 1000}, "multiplicity axis for histograms"};
  ConfigurableAxis axisCentrality{"axisCentrality", {100, 0, 100}, "centrality axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {98, 0.2, 10.0}, "pt axis for histograms"};

  // Axes for Correlations
  ConfigurableAxis axisMultForCorr{"axisMultForCorr", {VARIABLE_WIDTH, 0, 5, 10, 20, 50, 100, 150, 200, 300, 500}, "multiplicity axis for correlation histograms"};
  ConfigurableAxis axisCentForCorr{"axisCentForCorr", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, "centrality axis for correlation histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {48, -2.4, 2.4}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt trigger axis for histograms"};
  ConfigurableAxis axisParticle{"axisParticle", {4, 0, 4}, "particle axis for correlation container"};
  ConfigurableAxis axisVtxMix{"axisVtxMix", {VARIABLE_WIDTH, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, "vertex axis for mixed event histograms"};
  ConfigurableAxis axisMultMix{"axisMultMix", {VARIABLE_WIDTH, 0, 10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260}, "multiplicity / centrality axis for mixed event histograms"};

  // Axes for PID QA
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {80, -5, 5}, "nsigmaTPC axis"};
  ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {80, -5, 5}, "nsigmaTOF axis"};
  ConfigurableAxis axisNsigmaITS{"axisNsigmaITS", {80, -5, 5}, "nsigmaITS axis"};
  ConfigurableAxis axisTpcSignal{"axisTpcSignal", {250, 0, 250}, "dEdx axis for TPC"};
  ConfigurableAxis axisSigma{"axisSigma", {200, 0, 20}, "sigma axis for TPC"};

  // make the filters and cuts.
  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgEvSel.cfgGlobalCutVertex);
  Filter trackFilter = (nabs(aod::track::eta) < cfgTrkSel.cfgGlobalEta) && (aod::track::pt > cfgTrkSel.cfgGlobalPtMin) && (aod::track::pt < cfgTrkSel.cfgGlobalPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t)true)) && (aod::track::tpcChi2NCl < cfgTrkSel.cfgGlobalChi2prTPCcls) && (nabs(aod::track::dcaZ) < cfgTrkSel.cfgGlobalDCAz);

  using FilteredCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSel, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::Mults>>;
  using FilteredTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;
  using FilteredMcTracks = soa::Filtered<soa::Join<aod::Tracks, aod::McTrackLabels, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;
  using V0TrackCandidate = aod::V0Datas;

  Partition<FilteredTracks> asscTracks = (aod::track::pt > cfgGeneral.cfgAsscPtMin) && (aod::track::pt < cfgGeneral.cfgAsscPtMax);

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  Preslice<aod::V0Datas> perCollisionV0 = aod::v0::collisionId;

  // Corrections
  TH3D* mEfficiency = nullptr;
  TH1D* mEfficiencyNch = nullptr;
  TH1D* mCentralityWeight = nullptr;
  bool correctionsLoaded = false;

  // Define the outputs
  OutputObj<DihadronContainer> same{"sameEvent"};
  OutputObj<DihadronContainer> mixed{"mixedEvent"};
  OutputObj<DihadronContainer> sameReso{"sameEventReso"};
  OutputObj<DihadronContainer> mixedReso{"mixedEventReso"};
  HistogramRegistry histos{"histos"};

  // define global variables
  TRandom3 fRandom{0};
  enum PIDIndex {
    UseCharged = 0,
    UsePions,
    UseKaons,
    UseProtons,
    UseK0,
    UseLambda
  };
  enum PiKpArrayIndex {
    IndexPionUp = 0,
    IndexKaonUp,
    IndexProtonUp,
    IndexPionLow,
    IndexKaonLow,
    IndexProtonLow
  };
  enum ResoArrayIndex {
    IndexK0 = 0,
    IndexLambda,
    NumResoParticles
  };
  enum ResoParticleCuts {
    CutCosPA = 0,
    CutMassMin,
    CutMassMax,
    CutPosTrackPt,
    CutNegTrackPt,
    CutDCAPosToPVMin,
    CutDCANegToPVMin,
    CutLifeTime,
    CutRadiusMin,
    CutRadiusMax,
    CutRapidityMin,
    CutRapidityMax,
    CutArmPodMinVal,
    CutDCABetDaug,
    CutCompMassRej,
    NumParticleCuts
  };
  enum ResoParticleSwitches {
    UseCosPA = 0,
    UseMassBins,
    UseDCABetDaug,
    UseProperLifetime,
    UseV0Radius,
    UseArmPodCut,
    UseCompMassRej,
    NumParticleSwitches
  };

  enum ResoSelectionFlags {
    FillCandidate = 0,
    FillDaughterPt,
    FillMassCut,
    FillRapidityCut,
    FillDCAtoPV,
    FillDCAxDaughters,
    FillV0Radius,
    FillCosPA,
    FillProperLifetime,
    FillArmPodCut,
    FillCompetingMass,
    FillDaughterTrackSelected,
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
    NumEventCuts
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
    NumCentEstimators
  };
  enum EventType {
    SameEvent = 1,
    MixedEvent = 3
  };
  enum DetectorType {
    UseTPC = 0,
    UseTOF,
    UseITS
  };
  enum EventCutType {
    EvCut1 = 0,
    NumEvCutTypes = 1
  };
  enum Stage {
    Before = 0,
    After
  };

  std::array<std::array<float, 2>, 10> trackCuts{};
  std::array<std::array<int, 1>, 14> eventCuts{};
  std::array<std::array<float, 3>, 6> nSigmaVals{};
  std::array<std::array<float, 2>, 15> resoCutVals{};
  std::array<std::array<int, 2>, 7> resoSwitchVals{};
  // persistent caches
  std::vector<float> efficiencyAssociatedCache;
  std::vector<int> useAssociatedCache;

  void init(InitContext&)
  {
    readMatrix(cfgTrkSel.cfgTrackCuts->getData(), trackCuts);
    readMatrix(cfgEvSel.cfgUseEventCuts->getData(), eventCuts);
    readMatrix(cfgPIDConfigs.cfgPIDnSigmas->getData(), nSigmaVals);
    readMatrix(cfgPIDConfigs.cfgPIDResoCuts->getData(), resoCutVals);
    readMatrix(cfgPIDConfigs.cfgPIDResoSwitches->getData(), resoSwitchVals);

    if (cfgEvSel.cfgCentTableUnavailable && !cfgEvSel.cfgSelCollByNch) {
      LOGF(fatal, "Centrality table is unavailable, cannot select collisions by centrality");
    }
    const AxisSpec axisPhi{72, 0.0, constants::math::TwoPI, "#varphi"};
    const AxisSpec axisEta{40, -1., 1., "#eta"};

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    LOGF(info, "Starting init");

    // Creating mass axis depending on particle - 4 = kshort, 5 = lambda
    AxisSpec axisInvMass = {10, 0, 1, "mass"};
    if (cfgPIDConfigs.cfgPIDParticle == UseK0) {
      const auto bins = resoSwitchVals[UseMassBins][IndexK0];
      const auto min = resoCutVals[CutMassMin][IndexK0];
      const auto max = resoCutVals[CutMassMax][IndexK0];

      axisInvMass = {bins, min, max, "M_{#pi^{+}#pi^{-}} (GeV/c^{2})"};
    }
    if (cfgPIDConfigs.cfgPIDParticle == UseLambda) {
      const auto bins = resoSwitchVals[UseMassBins][IndexLambda];
      const auto min = resoCutVals[CutMassMin][IndexLambda];
      const auto max = resoCutVals[CutMassMax][IndexLambda];

      axisInvMass = {bins, min, max, "M_{p#pi} (GeV/c^{2})"};
    }

    // Set the sample axis
    auto maxSample = static_cast<double>(cfgGeneral.cfgSampleSize);
    AxisSpec axisSample{cfgGeneral.cfgSampleSize, 0, maxSample, "Sample"};

    // Choose if it is Nch selection or Centrality selection
    AxisSpec axisEventClass = {axisMultForCorr, "N_{ch}"};

    if (!cfgEvSel.cfgSelCollByNch) {
      axisEventClass = {axisCentForCorr, "Centrality"};
    }

    // Event Counter
    if ((doprocessSame || doprocessSameReso || doprocessMC) && cfgEvSel.cfgUseAdditionalEventCut) {
      histos.add("hEventCount", "Number of Events;; Count", {HistType::kTH1D, {{15, -0.5, 14.5}}});
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
    }

    if (cfgPIDConfigs.cfgPIDUseReso && cfgPIDConfigs.cfgPIDParticle == UseK0) { // For K0
      histos.add("PiPlusTPC_K0", "", {HistType::kTH2D, {{axisPtTrigger, axisNsigmaTPC}}});
      histos.add("PiMinusTPC_K0", "", {HistType::kTH2D, {{axisPtTrigger, axisNsigmaTPC}}});
      histos.add("PiPlusTOF_K0", "", {HistType::kTH2D, {{axisPtTrigger, axisNsigmaTOF}}});
      histos.add("PiMinusTOF_K0", "", {HistType::kTH2D, {{axisPtTrigger, axisNsigmaTOF}}});
      histos.add("hK0Phi", "", {HistType::kTH1D, {axisPhi}});
      histos.add("hK0Eta", "", {HistType::kTH1D, {axisEta}});

      histos.add("hK0Count", "Number of K0;; Count", {HistType::kTH1D, {{12, -0.5, 11.5}}});
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(FillCandidate + 1, "K0 candidates");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(FillDaughterPt + 1, "Daughter pt");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(FillMassCut + 1, "Mass cut");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(FillRapidityCut + 1, "Rapidity cut");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(FillDCAtoPV + 1, "DCA to PV");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(FillDCAxDaughters + 1, "DCA between daughters");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(FillV0Radius + 1, "V0radius");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(FillCosPA + 1, "CosPA");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(FillProperLifetime + 1, "Proper lifetime");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(FillArmPodCut + 1, "ArmenterosPod");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(FillCompetingMass + 1, "CompetingMassRejection");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(FillDaughterTrackSelected + 1, "Daughter track selection");
    }
    if (cfgPIDConfigs.cfgPIDUseReso && cfgPIDConfigs.cfgPIDParticle == UseLambda) { // For Lambda
      histos.add("PrPlusTPC_L", "", {HistType::kTH2D, {{axisPtTrigger, axisNsigmaTPC}}});
      histos.add("PiMinusTPC_L", "", {HistType::kTH2D, {{axisPtTrigger, axisNsigmaTPC}}});
      histos.add("PrPlusTOF_L", "", {HistType::kTH2D, {{axisPtTrigger, axisNsigmaTOF}}});
      histos.add("PiMinusTOF_L", "", {HistType::kTH2D, {{axisPtTrigger, axisNsigmaTOF}}});
      histos.add("hLambdaPhi", "", {HistType::kTH1D, {axisPhi}});
      histos.add("hLambdaEta", "", {HistType::kTH1D, {axisEta}});

      histos.add("hLambdaCount", "Number of Lambda;; Count", {HistType::kTH1D, {{12, -0.5, 11.5}}});
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(FillCandidate + 1, "Lambda candidates");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(FillDaughterPt + 1, "Daughter pt");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(FillMassCut + 1, "Mass cut");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(FillRapidityCut + 1, "Rapidity cut");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(FillDCAtoPV + 1, "DCA to PV");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(FillDCAxDaughters + 1, "DCA between daughters");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(FillV0Radius + 1, "V0radius");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(FillCosPA + 1, "CosPA");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(FillProperLifetime + 1, "Proper lifetime");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(FillArmPodCut + 1, "ArmenterosPod");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(FillCompetingMass + 1, "CompetingMassRejection");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(FillDaughterTrackSelected + 1, "Daughter track selection");
    }

    // Multiplicity correlation cuts
    if (eventCuts[UseMultCorrCut][EvCut1]) {
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

    std::string hCentTitle = "Centrality distribution, Estimator " + std::to_string(cfgEvSel.cfgCentEstimator);

    // Make histograms to check the distributions after cuts
    if (doprocessSame || doprocessSameReso) {
      histos.add("Nch", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}});
      histos.add("Nch_corrected", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}});
      histos.add("Centrality", hCentTitle.c_str(), {HistType::kTH1D, {axisCentrality}});
      histos.add("CentralityWeighted", hCentTitle.c_str(), {HistType::kTH1D, {{100, 0, 100}}});
      histos.add("zVtx", "zVtx", {HistType::kTH1D, {axisVertex}});
      histos.add("deltaEta_deltaPhi_same", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}}); // check to see the delta eta and delta phi distribution
    }

    if (doprocessSame) {
      if (cfgGeneral.cfgQABasic) {
        histos.add("Phi", "Phi", {HistType::kTH2D, {axisPhi, axisParticle}});
        histos.add("Eta", "Eta", {HistType::kTH2D, {axisEta, axisParticle}});
        histos.add("EtaCorrected", "EtaCorrected", {HistType::kTH2D, {axisEta, axisParticle}});
        histos.add("pTFiner", "pTFiner", {HistType::kTHnSparseF, {axisPt, axisParticle, axisEventClass}});
        histos.add("pTFinerCorrected", "pTFinerCorrected", {HistType::kTHnSparseF, {axisPt, axisParticle, axisEventClass}});
      }

      if (cfgPIDConfigs.cfgQANsigma) {
        if (!cfgPIDConfigs.cfgPIDUseIts) { // TPC-TOF PID QA hists
          histos.add("TofTpcNsigma_beforeCut", "", {HistType::kTHnSparseD, {{axisNsigmaTPC, axisNsigmaTOF, axisPtTrigger, axisParticle}}});
          histos.add("TofTpcNsigma_afterCut", "", {HistType::kTHnSparseD, {{axisNsigmaTPC, axisNsigmaTOF, axisPtTrigger, axisParticle}}});
        }
        if (cfgPIDConfigs.cfgPIDUseIts) { // ITS-TOF PID QA hists
          histos.add("TofItsNsigma_beforeCut", "", {HistType::kTHnSparseD, {{axisNsigmaITS, axisNsigmaTOF, axisPtTrigger, axisParticle}}});
          histos.add("TofItsNsigma_afterCut", "", {HistType::kTHnSparseD, {{axisNsigmaITS, axisNsigmaTOF, axisPtTrigger, axisParticle}}});
        }
      } // end of PID nsigma QA hists
      if (cfgPIDConfigs.cfgQAdEdx) {
        histos.add("TpcdEdx_ptwise_beforeCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF, axisParticle}}});
        histos.add("ExpTpcdEdx_ptwise_beforeCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF, axisParticle}}});
        histos.add("ExpSigma_ptwise_beforeCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisSigma, axisNsigmaTOF, axisParticle}}});

        histos.add("TpcdEdx_ptwise_afterCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF, axisParticle}}});
        histos.add("ExpTpcdEdx_ptwise_afterCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF, axisParticle}}});
        histos.add("ExpSigma_ptwise_afterCut", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisSigma, axisNsigmaTOF, axisParticle}}});
      } // end of PID dEdx QA hists

      // Trigger histogram
      histos.add("Trig_hist", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisEventClass, axisPtTrigger, axisParticle}}});
    }

    if (doprocessSameReso) {
      if (cfgPIDConfigs.cfgPIDParticle == UseK0 || cfgPIDConfigs.cfgPIDParticle == UseLambda) {
        histos.add("Trig_histReso", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisEventClass, axisPtTrigger, axisInvMass}}});
      }
    }

    if (doprocessMixed || doprocessMixedReso) {
      histos.add("deltaEta_deltaPhi_mixed", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}});
    }
    if (doprocessMC) {
      histos.add("hNsigmaPionTruePositives", "hNsigmaPionTruePositives", {HistType::kTH1D, {axisPtTrigger}});     // Fraction of particles that are pions and selected as pions
      histos.add("hNsigmaKaonTruePositives", "hNsigmaKaonTruePositives", {HistType::kTH1D, {axisPtTrigger}});     // Fraction of particles that are kaons and selected as kaons
      histos.add("hNsigmaProtonTruePositives", "hNsigmaProtonTruePositives", {HistType::kTH1D, {axisPtTrigger}}); // Fraction of particles that are protons and selected as protons

      histos.add("hNsigmaPionSelected", "hNsigmaPionSelected", {HistType::kTH1D, {axisPtTrigger}});
      histos.add("hNsigmaKaonSelected", "hNsigmaKaonSelected", {HistType::kTH1D, {axisPtTrigger}});
      histos.add("hNsigmaProtonSelected", "hNsigmaProtonSelected", {HistType::kTH1D, {axisPtTrigger}});

      histos.add("hNsigmaPionTrue", "hNsigmaPionTrue", {HistType::kTH1D, {axisPtTrigger}});     // All true pions from MC
      histos.add("hNsigmaKaonTrue", "hNsigmaKaonTrue", {HistType::kTH1D, {axisPtTrigger}});     // All true kaons from MC
      histos.add("hNsigmaProtonTrue", "hNsigmaProtonTrue", {HistType::kTH1D, {axisPtTrigger}}); // All true protons from MC

      histos.add("TpcdEdx_ptwise_pi", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF}}});
      histos.add("TpcdEdx_ptwise_ka", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF}}});
      histos.add("TpcdEdx_ptwise_pr", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF}}});

      histos.add("ExpTpcdEdx_ptwise_pi", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF}}});
      histos.add("ExpTpcdEdx_ptwise_ka", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF}}});
      histos.add("ExpTpcdEdx_ptwise_pr", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF}}});

      histos.add("TpcdEdx_ptwise_pi_physprim", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF}}});
      histos.add("TpcdEdx_ptwise_ka_physprim", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF}}});
      histos.add("TpcdEdx_ptwise_pr_physprim", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF}}});

      histos.add("TpcdEdx_ptwise_all_physprim", "", {HistType::kTHnSparseD, {{axisPtTrigger, axisTpcSignal, axisNsigmaTOF}}});
    }

    histos.add("eventcount", "bin", {HistType::kTH1F, {{4, 0, 4, "bin"}}}); // histogram to see how many events are in the same and mixed event

    LOGF(info, "Initializing correlation container");
    std::vector<AxisSpec> corrAxis = {{axisSample},
                                      {axisVertex, "z-vtx (cm)"},
                                      {axisEventClass},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisParticle, "Particle, 0 = charged, 1 = pion, 2 = kaon, 3 = proton"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisDeltaEta, "#Delta#eta"}};

    std::vector<AxisSpec> corrAxisReso = {{axisSample},
                                          {axisVertex, "z-vtx (cm)"},
                                          {axisEventClass},
                                          {axisPtTrigger, "p_{T} (GeV/c)"},
                                          {axisInvMass},
                                          {axisDeltaPhi, "#Delta#varphi (rad)"},
                                          {axisDeltaEta, "#Delta#eta"}};

    if (!cfgPIDConfigs.cfgPIDUseReso) {
      same.setObject(new DihadronContainer("sameEvent", "sameEvent", corrAxis));
      mixed.setObject(new DihadronContainer("mixedEvent", "mixedEvent", corrAxis));
    }

    if (cfgPIDConfigs.cfgPIDUseReso && (cfgPIDConfigs.cfgPIDParticle == UseK0 || cfgPIDConfigs.cfgPIDParticle == UseLambda)) {
      sameReso.setObject(new DihadronContainer("sameEvent", "sameEvent", corrAxisReso));
      mixedReso.setObject(new DihadronContainer("mixedEvent", "mixedEvent", corrAxisReso));
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

  template <typename TTrack>
  bool selectionV0Daughter(TTrack const& track, int pid)
  {
    if (!(track.itsNCls() >= trackCuts[TrkCutItsCluster][UseGenTrkCuts])) {
      return false;
    }
    if (!track.hasTPC()) {
      return false;
    }
    if (!(track.tpcNClsFound() >= trackCuts[TrkCutTpcCluster][UseGenTrkCuts])) {
      return false;
    }
    if (!(track.tpcNClsCrossedRows() >= trackCuts[TrkCutTpcCrossedRows][UseGenTrkCuts])) {
      return false;
    }

    if (cfgPIDConfigs.cfgPIDUseOnlyTpc) {
      if (pid == UsePions && std::abs(track.tpcNSigmaPi()) > cfgPIDConfigs.cfgPIDOnlyTpcCut) {
        return false;
      }
      if (pid == UseKaons && std::abs(track.tpcNSigmaKa()) > cfgPIDConfigs.cfgPIDOnlyTpcCut) {
        return false;
      }
      if (pid == UseProtons && std::abs(track.tpcNSigmaPr()) > cfgPIDConfigs.cfgPIDOnlyTpcCut) {
        return false;
      }
    } else {
      int partIndex = getNsigmaPID(track);
      int pidIndex = partIndex; // 1 = pion, 2 = kaon, 3 = proton
      if (pidIndex != pid)
        return false;
    }

    return true;
  }

  template <typename TTrack>
  int getNsigmaPID(const TTrack& track)
  {
    // Computing Nsigma arrays for pion, kaon, and protons
    std::array<float, 3> nSigmaTPC = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::array<float, 3> nSigmaTOF = {track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};
    std::array<float, 3> nSigmaITS = {itsResponse.nSigmaITS<o2::track::PID::Pion>(track), itsResponse.nSigmaITS<o2::track::PID::Kaon>(track), itsResponse.nSigmaITS<o2::track::PID::Proton>(track)};
    int pid = -1; // -1 = not identified, 1 = pion, 2 = kaon, 3 = proton

    std::array<float, 3> nSigmaToUse = cfgPIDConfigs.cfgPIDUseIts ? nSigmaITS : nSigmaTPC; // Choose which nSigma to use: TPC or ITS
    int useIndexDetector = UseTPC;
    useIndexDetector = cfgPIDConfigs.cfgPIDUseIts ? UseITS : UseTPC; // Choose which nSigma to use: TPC or ITS

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

    isTofPion = nSigmaTOF[IndexPionUp] < nSigmaVals[IndexPionUp][UseTOF] && nSigmaTOF[IndexPionUp] > nSigmaVals[IndexPionLow][UseTOF];
    isTofKaon = nSigmaTOF[IndexKaonUp] < nSigmaVals[IndexKaonUp][UseTOF] && nSigmaTOF[IndexKaonUp] > nSigmaVals[IndexKaonLow][UseTOF];
    isTofProton = nSigmaTOF[IndexProtonUp] < nSigmaVals[IndexProtonUp][UseTOF] && nSigmaTOF[IndexProtonUp] > nSigmaVals[IndexProtonLow][UseTOF];

    if (track.pt() > cfgPIDConfigs.cfgPIDTofPtCut && !track.hasTOF()) {
      return -1;
    } else if (track.pt() > cfgPIDConfigs.cfgPIDTofPtCut && track.hasTOF()) {
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
    if (eff == 0) {
      return false;
    }
    weight_nue = 1. / eff;
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
  void nchCounter(const TTracks& tracks, double& multiplicity) // function to count the number of tracks in the event and fill the histogram
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

    if (pid > UseCharged && (cfgPIDConfigs.cfgQANsigma || cfgPIDConfigs.cfgQAdEdx)) {
      const bool useITS = cfgPIDConfigs.cfgPIDUseIts;
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

      if (cfgPIDConfigs.cfgQAdEdx) {
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

  template <typename TTracks, typename TTracksAssoc>
  void fillCorrelations(const TTracks& triggers, const TTracksAssoc& assocs, float posZ, int system, int magneticField, double eventClass, float eventWeight) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {
    const int kNullIndex = -44;
    efficiencyAssociatedCache.clear();
    efficiencyAssociatedCache.reserve(assocs.size());
    useAssociatedCache.clear();
    useAssociatedCache.reserve(assocs.size());

    for (const auto& assc : assocs) { // loading efficiencies for associated particles and simultaneously applying track selection
      if (trackSelected(assc)) {
        useAssociatedCache.push_back(1);
        if (mEfficiency) {
          float weff = 1.;
          getEfficiencyCorrection(weff, assc.eta(), assc.pt(), posZ);
          efficiencyAssociatedCache.push_back(weff);
        }
      } else {
        useAssociatedCache.push_back(0);
        if (mEfficiency) {
          efficiencyAssociatedCache.push_back(kNullIndex);
        }
      }
    } // end of efficieny requirement

    int fSampleIndex = static_cast<int>(fRandom.Uniform(0.0, cfgGeneral.cfgSampleSize));

    float triggerWeight = 1.0f;
    float associatedWeight = 1.0f;
    // loop over all tracks
    for (auto const& trig : triggers) {

      if (!trackSelected(trig)) {
        continue;
      }

      if (cfgEffWeights.cfgUseTrigWeight && !getEfficiencyCorrection(triggerWeight, trig.eta(), trig.pt(), posZ)) {
        continue;
      }

      if (system == SameEvent) {
        fillTrackQA(trig, UseCharged, eventClass, Before, triggerWeight);
        fillTrackQA(trig, UsePions, eventClass, Before, triggerWeight);
        fillTrackQA(trig, UseKaons, eventClass, Before, triggerWeight);
        fillTrackQA(trig, UseProtons, eventClass, Before, triggerWeight);
      }

      int pidIndex = getNsigmaPID(trig);

      if (system == SameEvent) {
        histos.fill(HIST("Trig_hist"), fSampleIndex, posZ, eventClass, trig.pt(), UseCharged, eventWeight * triggerWeight);
        if (pidIndex > UseCharged) {
          fillTrackQA(trig, pidIndex, eventClass, After, triggerWeight);
          histos.fill(HIST("Trig_hist"), fSampleIndex, posZ, eventClass, trig.pt(), pidIndex, eventWeight * triggerWeight);
        }
      }

      for (auto const& assc : assocs) {

        if (useAssociatedCache[assc.filteredIndex()] == 0) { // rejecting associated tracks
          continue;
        }

        if (mEfficiency && efficiencyAssociatedCache[assc.filteredIndex()] != kNullIndex) { // loading associated efficiency
          associatedWeight = efficiencyAssociatedCache[assc.filteredIndex()];
        }

        if (!cfgGeneral.cfgUsePtOrder && trig.globalIndex() == assc.globalIndex()) {
          continue; // For pt-differential correlations, skip if the trigger and associate are the same track
        }
        if (cfgGeneral.cfgUsePtOrder && system == SameEvent && trig.pt() <= assc.pt()) {
          continue; // Without pt-differential correlations, skip if the trigger pt is less than the associate pt
        }
        if (cfgGeneral.cfgUsePtOrder && system == MixedEvent && cfgGeneral.cfgUsePtOrderInMixEvent && trig.pt() <= assc.pt()) {
          continue; // For pt-differential correlations in mixed events, skip if the trigger pt is less than the associate pt
        }

        float deltaPhi = RecoDecay::constrainAngle(trig.phi() - assc.phi(), -PIHalf);
        float deltaEta = trig.eta() - assc.eta();

        if (std::abs(deltaEta) < cfgGeneral.cfgCutMerging && cfgGeneral.cfgUseCutMerging) {

          double dPhiStarHigh = getDPhiStar(trig.sign(), assc.sign(), trig.phi(), assc.phi(), trig.pt(), assc.pt(), cfgGeneral.cfgRadiusHigh, magneticField);
          double dPhiStarLow = getDPhiStar(trig.sign(), assc.sign(), trig.phi(), assc.phi(), trig.pt(), assc.pt(), cfgGeneral.cfgRadiusLow, magneticField);

          const double kLimit = 3.0 * cfgGeneral.cfgCutMerging;

          bool bIsBelow = false;

          if (std::abs(dPhiStarLow) < kLimit || std::abs(dPhiStarHigh) < kLimit || dPhiStarLow * dPhiStarHigh < 0) {
            for (double rad(cfgGeneral.cfgRadiusLow); rad < cfgGeneral.cfgRadiusHigh; rad += 0.01) {
              double dPhiStar = getDPhiStar(trig.sign(), assc.sign(), trig.phi(), assc.phi(), trig.pt(), assc.pt(), rad, magneticField);
              if (std::abs(dPhiStar) < kLimit) {
                bIsBelow = true;
                break;
              }
            }
            if (bIsBelow) {
              continue;
            }
          }
        } // End of merging cut condition

        // fill the right sparse and histograms
        if (system == SameEvent) {
          same->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, trig.pt(), UseCharged, deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          if (pidIndex > UseCharged) {
            same->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, trig.pt(), pidIndex, deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          }
          histos.fill(HIST("deltaEta_deltaPhi_same"), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
        } else if (system == MixedEvent) {
          mixed->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, trig.pt(), UseCharged, deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          if (pidIndex > UseCharged) {
            mixed->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, trig.pt(), pidIndex, deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          }
          histos.fill(HIST("deltaEta_deltaPhi_mixed"), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
        }
      }
    }
  }

  template <typename TV0Tracks, typename TTracksAssoc>
  void fillCorrelationsReso(const TV0Tracks& triggers, const TTracksAssoc& assocs, float posZ, float posY, float posX, int system, int magneticField, double eventClass, float eventWeight) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {
    const int kNullIndex = -44;
    efficiencyAssociatedCache.clear();
    efficiencyAssociatedCache.reserve(assocs.size());
    useAssociatedCache.clear();
    useAssociatedCache.reserve(assocs.size());

    for (const auto& assc : assocs) { // loading efficiencies for associated particles and simultaneously applying track selection
      if (trackSelected(assc)) {
        useAssociatedCache.push_back(1);
        if (mEfficiency) {
          float weff = 1.;
          getEfficiencyCorrection(weff, assc.eta(), assc.pt(), posZ);
          efficiencyAssociatedCache.push_back(weff);
        }
      } else {
        useAssociatedCache.push_back(0);
        if (mEfficiency) {
          efficiencyAssociatedCache.push_back(kNullIndex);
        }
      }
    } // end of efficieny requirement

    int fSampleIndex = static_cast<int>(fRandom.Uniform(0.0, cfgGeneral.cfgSampleSize));

    float triggerWeight = 1.0f;
    float associatedWeight = 1.0f;
    // loop over all tracks
    for (auto const& trig : triggers) {

      double resoMass = -1;

      // 4 = kshort, 5 = lambda
      if (cfgPIDConfigs.cfgPIDParticle == UseK0) {
        if (!isSelectedK0(trig, posZ, posY, posX, system)) {
          continue; // Reject if called for K0 but V0 is not K0
        }

        resoMass = trig.mK0Short();
      }

      if (cfgPIDConfigs.cfgPIDParticle == UseLambda) {
        if (!isSelectedLambda(trig, posZ, posY, posX, system)) {
          continue; // Reject if called for Lambda but V0 is not lambda
        }
        resoMass = trig.mLambda();
      }

      if (system == SameEvent) {
        histos.fill(HIST("Trig_histReso"), fSampleIndex, posZ, eventClass, trig.pt(), resoMass, eventWeight * triggerWeight);
      }

      for (auto const& assc : assocs) {

        if (useAssociatedCache[assc.filteredIndex()] == 0) {
          continue; // rejecting associated tracks
        }

        if (mEfficiency && efficiencyAssociatedCache[assc.filteredIndex()] != kNullIndex) {
          associatedWeight = efficiencyAssociatedCache[assc.filteredIndex()]; // loading associated efficiency
        }

        if (cfgGeneral.cfgUsePtOrder && system == SameEvent && trig.pt() <= assc.pt()) {
          continue; // Without pt-differential correlations, skip if the trigger pt is less than the associate pt
        }
        if (cfgGeneral.cfgUsePtOrder && system == MixedEvent && cfgGeneral.cfgUsePtOrderInMixEvent && trig.pt() <= assc.pt()) {
          continue; // For pt-differential correlations in mixed events, skip if the trigger pt is less than the associate pt
        }

        float deltaPhi = RecoDecay::constrainAngle(trig.phi() - assc.phi(), -PIHalf);
        float deltaEta = trig.eta() - assc.eta();

        if (std::abs(deltaEta) < cfgGeneral.cfgCutMerging && cfgGeneral.cfgUseCutMerging) {

          double dPhiStarHigh = getDPhiStar(0, assc.sign(), trig.phi(), assc.phi(), trig.pt(), assc.pt(), cfgGeneral.cfgRadiusHigh, magneticField);
          double dPhiStarLow = getDPhiStar(0, assc.sign(), trig.phi(), assc.phi(), trig.pt(), assc.pt(), cfgGeneral.cfgRadiusLow, magneticField);

          const double kLimit = 3.0 * cfgGeneral.cfgCutMerging;

          bool bIsBelow = false;

          if (std::abs(dPhiStarLow) < kLimit || std::abs(dPhiStarHigh) < kLimit || dPhiStarLow * dPhiStarHigh < 0) {
            for (double rad(cfgGeneral.cfgRadiusLow); rad < cfgGeneral.cfgRadiusHigh; rad += 0.01) {
              double dPhiStar = getDPhiStar(0, assc.sign(), trig.phi(), assc.phi(), trig.pt(), assc.pt(), rad, magneticField);
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
          sameReso->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, trig.pt(), resoMass, deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          histos.fill(HIST("deltaEta_deltaPhi_same"), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
        } else if (system == MixedEvent) {
          mixedReso->getCorrHist()->Fill(0, fSampleIndex, posZ, eventClass, trig.pt(), resoMass, deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          histos.fill(HIST("deltaEta_deltaPhi_mixed"), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
        }
      }
    }
  }

  template <typename TCollision>
  bool selectionEvent(TCollision const& collision, const int mult, const double cent, const bool fillCounter)
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

    if (eventCuts[UseOccupancy][EvCut1] && (occupancy < cfgEvSel.cfgCutOccupancyMin || occupancy > cfgEvSel.cfgCutOccupancyMax)) {
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

  void processSame(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float weightCent = 1.0f;
    float cent = cfgEvSel.cfgCentTableUnavailable ? -1.f : getCentrality(collision);

    if (cfgEvSel.cfgUseAdditionalEventCut && !selectionEvent(collision, tracks.size(), cent, true)) {
      return;
    }

    loadCorrection(bc.timestamp());
    if (!cfgEvSel.cfgCentTableUnavailable) {
      getCentralityWeight(weightCent, cent);
      histos.fill(HIST("Centrality"), cent);
      histos.fill(HIST("CentralityWeighted"), cent, weightCent);
    }
    histos.fill(HIST("Nch"), tracks.size());
    histos.fill(HIST("zVtx"), collision.posZ());

    auto multiplicity = static_cast<double>(tracks.size());

    if (cfgGeneral.cfgStrictTrackCounter) {
      nchCounter(tracks, multiplicity);
    }
    histos.fill(HIST("Nch_corrected"), multiplicity);

    if (cfgEvSel.cfgSelCollByNch && (multiplicity < cfgEvSel.cfgCutMultMin || multiplicity >= cfgEvSel.cfgCutMultMax)) {
      return;
    }
    if (!cfgEvSel.cfgSelCollByNch && !cfgEvSel.cfgCentTableUnavailable && (cent < cfgEvSel.cfgCutCentMin || cent >= cfgEvSel.cfgCutCentMax)) {
      return;
    }
    histos.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin

    auto eventClass = (cfgEvSel.cfgSelCollByNch) ? multiplicity : cent;

    auto assocs = asscTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    fillCorrelations(tracks, assocs, collision.posZ(), SameEvent, getMagneticField(bc.timestamp()), eventClass, weightCent);
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

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}};

    auto tracksTuple = std::make_tuple(tracks, tracks);
    Pair<FilteredCollisions, FilteredTracks, FilteredTracks, MixedBinning> pairs{binningOnVtxAndMult, cfgGeneral.cfgMixEventNumMin, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;
      if (!collision1.sel8() || !collision2.sel8()) {
        continue;
      }

      float cent1 = cfgEvSel.cfgCentTableUnavailable ? -1.f : getCentrality(collision1);
      float cent2 = cfgEvSel.cfgCentTableUnavailable ? -1.f : getCentrality(collision2);

      if (cfgEvSel.cfgUseAdditionalEventCut && !selectionEvent(collision1, tracks1.size(), cent1, false)) {
        continue;
      }
      if (cfgEvSel.cfgUseAdditionalEventCut && !selectionEvent(collision2, tracks2.size(), cent2, false)) {
        continue;
      }

      if (!cfgEvSel.cfgSelCollByNch && !cfgEvSel.cfgCentTableUnavailable && (cent1 < cfgEvSel.cfgCutCentMin || cent1 >= cfgEvSel.cfgCutCentMax)) {
        continue;
      }

      if (!cfgEvSel.cfgSelCollByNch && !cfgEvSel.cfgCentTableUnavailable && (cent2 < cfgEvSel.cfgCutCentMin || cent2 >= cfgEvSel.cfgCutCentMax)) {
        continue;
      }

      auto multiplicity1 = static_cast<double>(tracks1.size());
      auto multiplicity2 = static_cast<double>(tracks2.size());

      if (cfgGeneral.cfgStrictTrackCounter) {
        nchCounter(tracks1, multiplicity1);
        if (cfgEvSel.cfgUseCol2ForMixedEvSel) {
          nchCounter(tracks2, multiplicity2);
        }
      }

      if (cfgEvSel.cfgSelCollByNch && (multiplicity1 < cfgEvSel.cfgCutMultMin || multiplicity1 >= cfgEvSel.cfgCutMultMax)) {
        continue;
      }

      if (cfgEvSel.cfgSelCollByNch && (multiplicity2 < cfgEvSel.cfgCutMultMin || multiplicity2 >= cfgEvSel.cfgCutMultMax)) {
        continue;
      }

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      loadCorrection(bc.timestamp());
      float eventWeight = 1.0f;
      if (cfgEffWeights.cfgUseEventWeights) {
        eventWeight = 1.0f / it.currentWindowNeighbours();
      }

      histos.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin
      auto eventClass = (cfgEvSel.cfgSelCollByNch) ? multiplicity1 : cent1;
      float weightCent = 1.0f;
      if (!cfgEvSel.cfgCentTableUnavailable) {
        getCentralityWeight(weightCent, cent1);
      }

      auto assocs = asscTracks->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);
      fillCorrelations(tracks1, assocs, collision1.posZ(), MixedEvent, getMagneticField(bc.timestamp()), eventClass, eventWeight * weightCent);
    }
  }
  PROCESS_SWITCH(PidDiHadron, processMixed, "Process mixed events", true);

  template <typename V0>
  bool isSelectedK0(V0 const& candidate, float posZ, float posY, float posX, int system)
  {
    double mk0 = candidate.mK0Short();

    // separate the positive and negative V0 daughters
    auto postrack = candidate.template posTrack_as<FilteredTracks>();
    auto negtrack = candidate.template negTrack_as<FilteredTracks>();

    if (system == SameEvent) {
      histos.fill(HIST("hK0Count"), FillCandidate);
    }
    if (postrack.pt() < resoCutVals[CutPosTrackPt][IndexK0] || negtrack.pt() < resoCutVals[CutNegTrackPt][IndexK0]) {
      return false;
    }
    if (system == SameEvent) {
      histos.fill(HIST("hK0Count"), FillDaughterPt);
    }
    if (mk0 < resoCutVals[CutMassMin][IndexK0] || mk0 > resoCutVals[CutMassMax][IndexK0]) {
      return false;
    }
    if (system == SameEvent) {
      histos.fill(HIST("hK0Count"), FillMassCut);
    }
    // Rapidity correction
    if (candidate.yK0Short() > resoCutVals[CutRapidityMax][IndexK0] || candidate.yK0Short() < resoCutVals[CutRapidityMin][IndexK0]) {
      return false;
    }
    if (system == SameEvent) {
      histos.fill(HIST("hK0Count"), FillRapidityCut);
    }
    // DCA cuts for K0short
    if (std::abs(candidate.dcapostopv()) < resoCutVals[CutDCAPosToPVMin][IndexK0] || std::abs(candidate.dcanegtopv()) < resoCutVals[CutDCANegToPVMin][IndexK0]) {
      return false;
    }
    if (system == SameEvent) {
      histos.fill(HIST("hK0Count"), FillDCAtoPV);
    }
    if (resoSwitchVals[UseDCABetDaug][IndexK0] && std::abs(candidate.dcaV0daughters()) > resoCutVals[CutDCABetDaug][IndexK0]) {
      return false;
    }
    if (system == SameEvent) {
      histos.fill(HIST("hK0Count"), FillDCAxDaughters);
    }
    // v0 radius cuts
    if (resoSwitchVals[UseV0Radius][IndexK0] && (candidate.v0radius() < resoCutVals[CutRadiusMin][IndexK0] || candidate.v0radius() > resoCutVals[CutRadiusMax][IndexK0])) {
      return false;
    }
    if (system == SameEvent) {
      histos.fill(HIST("hK0Count"), FillV0Radius);
    }
    // cosine pointing angle cuts
    if (resoSwitchVals[UseCosPA][IndexK0] && candidate.v0cosPA() < resoCutVals[CutCosPA][IndexK0]) {
      return false;
    }
    if (system == SameEvent) {
      histos.fill(HIST("hK0Count"), FillCosPA);
    }
    // Proper lifetime
    float cTauK0 = candidate.distovertotmom(posX, posY, posZ) * massK0Short;
    if (resoSwitchVals[UseProperLifetime][IndexK0] && cTauK0 > resoCutVals[CutLifeTime][IndexK0]) {
      return false;
    }
    if (system == SameEvent) {
      histos.fill(HIST("hK0Count"), FillProperLifetime);
    }
    // ArmenterosPodolanskiCut
    if (resoSwitchVals[UseArmPodCut][IndexK0] && (candidate.qtarm() / std::abs(candidate.alpha())) < resoCutVals[CutArmPodMinVal][IndexK0]) {
      return false;
    }
    if (system == SameEvent) {
      histos.fill(HIST("hK0Count"), FillArmPodCut);
    }

    if (resoSwitchVals[UseCompMassRej][IndexK0]) {
      if (std::abs(candidate.mLambda() - o2::constants::physics::MassLambda) < resoCutVals[CutCompMassRej][IndexK0]) {
        return false;
      }
      if (std::abs(candidate.mAntiLambda() - o2::constants::physics::MassLambda) < resoCutVals[CutCompMassRej][IndexK0]) {
        return false;
      }
    }
    if (system == SameEvent) {
      histos.fill(HIST("hK0Count"), FillCompetingMass);
    }

    // Selection on V0 daughters
    if (!selectionV0Daughter(postrack, UsePions) || !selectionV0Daughter(negtrack, UsePions)) {
      return false;
    }
    if (system == SameEvent) {
      histos.fill(HIST("hK0Count"), FillDaughterTrackSelected);
      histos.fill(HIST("hK0Phi"), candidate.phi());
      histos.fill(HIST("hK0Rap"), candidate.yK0Short());
      histos.fill(HIST("hK0Eta"), candidate.eta());
      histos.fill(HIST("PiPlusTPC_K0"), postrack.pt(), postrack.tpcNSigmaPi());
      histos.fill(HIST("PiPlusTOF_K0"), postrack.pt(), postrack.tofNSigmaPi());
      histos.fill(HIST("PiMinusTPC_K0"), negtrack.pt(), negtrack.tpcNSigmaPi());
      histos.fill(HIST("PiMinusTOF_K0"), negtrack.pt(), negtrack.tofNSigmaPi());
    }

    return true;
  }

  template <typename V0>
  bool isSelectedLambda(V0 const& candidate, float posZ, float posY, float posX, int system)
  {
    bool isL = false;  // Is lambda candidate
    bool isAL = false; // Is anti-lambda candidate

    double mlambda = candidate.mLambda();
    double mantilambda = candidate.mAntiLambda();

    // separate the positive and negative V0 daughters
    auto postrack = candidate.template posTrack_as<FilteredTracks>();
    auto negtrack = candidate.template negTrack_as<FilteredTracks>();

    if (system == SameEvent) {
      histos.fill(HIST("hLambdaCount"), FillCandidate);
    }
    if (postrack.pt() < resoCutVals[CutPosTrackPt][IndexLambda] || negtrack.pt() < resoCutVals[CutNegTrackPt][IndexLambda]) {
      return false;
    }
    if (system == SameEvent) {
      histos.fill(HIST("hLambdaCount"), FillDaughterPt);
    }
    if (mlambda > resoCutVals[CutMassMin][IndexLambda] && mlambda < resoCutVals[CutMassMax][IndexLambda]) {
      isL = true;
    }
    if (mantilambda > resoCutVals[CutMassMin][IndexLambda] && mantilambda < resoCutVals[CutMassMax][IndexLambda]) {
      isAL = true;
    }

    if (!isL && !isAL) {
      return false;
    }
    if (system == SameEvent) {
      histos.fill(HIST("hLambdaCount"), FillMassCut);
    }

    // Rapidity correction
    if (candidate.yLambda() > resoCutVals[CutRapidityMax][IndexLambda] || candidate.yLambda() < resoCutVals[CutRapidityMin][IndexLambda]) {
      return false;
    }
    if (system == SameEvent) {
      histos.fill(HIST("hLambdaCount"), FillRapidityCut);
    }

    // DCA cuts for lambda and antilambda
    if (isL && (std::abs(candidate.dcapostopv()) < resoCutVals[CutDCAPosToPVMin][IndexLambda] || std::abs(candidate.dcanegtopv()) < resoCutVals[CutDCANegToPVMin][IndexLambda])) {
      isL = false;
    }
    if (isAL && (std::abs(candidate.dcapostopv()) < resoCutVals[CutDCANegToPVMin][IndexLambda] || std::abs(candidate.dcanegtopv()) < resoCutVals[CutDCAPosToPVMin][IndexLambda])) {
      isAL = false;
    }
    if (!isL && !isAL) {
      return false;
    }
    if (system == SameEvent) {
      histos.fill(HIST("hLambdaCount"), FillDCAtoPV);
    }
    if (resoSwitchVals[UseDCABetDaug][IndexLambda] && std::abs(candidate.dcaV0daughters()) > resoCutVals[CutDCABetDaug][IndexLambda]) {
      return false;
    }
    if (system == SameEvent) {
      histos.fill(HIST("hLambdaCount"), FillDCAxDaughters);
    }
    // v0 radius cuts
    if (resoSwitchVals[UseV0Radius][IndexLambda] && (candidate.v0radius() < resoCutVals[CutRadiusMin][IndexLambda] || candidate.v0radius() > resoCutVals[CutRadiusMax][IndexLambda])) {
      return false;
    }
    if (system == SameEvent) {
      histos.fill(HIST("hLambdaCount"), FillV0Radius);
    }
    // cosine pointing angle cuts
    if (resoSwitchVals[UseCosPA][IndexLambda] && candidate.v0cosPA() < resoCutVals[CutCosPA][IndexLambda]) {
      return false;
    }
    if (system == SameEvent) {
      histos.fill(HIST("hLambdaCount"), FillCosPA);
    }
    // Proper lifetime
    float cTauLambda = candidate.distovertotmom(posX, posY, posZ) * massLambda;
    if (resoSwitchVals[UseProperLifetime][IndexLambda] && cTauLambda > resoCutVals[CutLifeTime][IndexLambda]) {
      return false;
    }
    if (system == SameEvent) {
      histos.fill(HIST("hLambdaCount"), FillProperLifetime);
    }

    // ArmenterosPodolanskiCut
    if (resoSwitchVals[UseArmPodCut][IndexLambda] && (candidate.qtarm() / std::abs(candidate.alpha())) < resoCutVals[CutArmPodMinVal][IndexLambda]) {
      return false;
    }
    if (system == SameEvent) {
      histos.fill(HIST("hLambdaCount"), FillArmPodCut);
    }

    if (resoSwitchVals[UseCompMassRej][IndexLambda]) {
      if (std::abs(candidate.mK0Short() - o2::constants::physics::MassK0Short) < resoCutVals[CutCompMassRej][IndexLambda]) {
        return false;
      }
    }
    if (system == SameEvent) {
      histos.fill(HIST("hLambdaCount"), FillCompetingMass);
    }

    if (isL) {
      if (!selectionV0Daughter(postrack, UseProtons) || !selectionV0Daughter(negtrack, UsePions)) {
        isL = false;
      }
    }
    if (isAL) {
      if (!selectionV0Daughter(postrack, UsePions) || !selectionV0Daughter(negtrack, UseProtons)) {
        isAL = false;
      }
    }
    if (!isL && !isAL) {
      return false;
    }
    if (system == SameEvent) {
      histos.fill(HIST("hLambdaCount"), FillDaughterTrackSelected);
      histos.fill(HIST("hLambdaPhi"), candidate.phi());
      histos.fill(HIST("hLambdaEta"), candidate.eta());
      histos.fill(HIST("hLambdaRap"), candidate.yLambda());
    }
    if (isL && system == SameEvent) {
      histos.fill(HIST("PrPlusTPC_La"), postrack.pt(), postrack.tpcNSigmaPr());
      histos.fill(HIST("PrPlusTOF_La"), postrack.pt(), postrack.tofNSigmaPr());
      histos.fill(HIST("PiMinusTPC_La"), negtrack.pt(), negtrack.tpcNSigmaPi());
      histos.fill(HIST("PiMinusTOF_La"), negtrack.pt(), negtrack.tofNSigmaPi());
    }
    if (isAL && system == SameEvent) {
      histos.fill(HIST("PrMinusTPC_Al"), negtrack.pt(), negtrack.tpcNSigmaPr());
      histos.fill(HIST("PrMinusTOF_Al"), negtrack.pt(), negtrack.tofNSigmaPr());
      histos.fill(HIST("PiPlusTPC_Al"), postrack.pt(), postrack.tpcNSigmaPi());
      histos.fill(HIST("PiPlusTOF_Al"), postrack.pt(), postrack.tofNSigmaPi());
    }

    return true;
  }

  double massKaPlus = o2::constants::physics::MassKPlus;    // same as MassKMinus
  double massLambda = o2::constants::physics::MassLambda;   // same as MassLambda0 and MassLambda0Bar
  double massK0Short = o2::constants::physics::MassK0Short; // same as o2::constants::physics::MassK0 and o2::constants::physics::MassK0Bar

  void processSameReso(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::BCsWithTimestamps const&, aod::V0Datas const& V0s)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float weightCent = 1.0f;
    float cent = cfgEvSel.cfgCentTableUnavailable ? -1.f : getCentrality(collision);

    if (cfgEvSel.cfgUseAdditionalEventCut && !selectionEvent(collision, tracks.size(), cent, true)) {
      return;
    }

    loadCorrection(bc.timestamp());
    if (!cfgEvSel.cfgCentTableUnavailable) {
      getCentralityWeight(weightCent, cent);
      histos.fill(HIST("Centrality"), cent);
      histos.fill(HIST("CentralityWeighted"), cent, weightCent);
    }
    histos.fill(HIST("Nch"), tracks.size());
    histos.fill(HIST("zVtx"), collision.posZ());

    auto multiplicity = static_cast<double>(tracks.size());

    if (cfgGeneral.cfgStrictTrackCounter) {
      nchCounter(tracks, multiplicity);
    }

    histos.fill(HIST("Nch_corrected"), multiplicity);

    if (cfgEvSel.cfgSelCollByNch && (multiplicity < cfgEvSel.cfgCutMultMin || multiplicity >= cfgEvSel.cfgCutMultMax)) {
      return;
    }
    if (!cfgEvSel.cfgSelCollByNch && !cfgEvSel.cfgCentTableUnavailable && (cent < cfgEvSel.cfgCutCentMin || cent >= cfgEvSel.cfgCutCentMax)) {
      return;
    }

    histos.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin

    auto eventClass = (cfgEvSel.cfgSelCollByNch) ? multiplicity : cent;

    auto assocs = asscTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    fillCorrelationsReso(V0s, assocs, collision.posZ(), collision.posY(), collision.posX(), SameEvent, getMagneticField(bc.timestamp()), eventClass, weightCent);
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

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}};

    auto tracksTuple = std::make_tuple(V0s, tracks);
    Pair<FilteredCollisions, aod::V0Datas, FilteredTracks, MixedBinning> pairs{binningOnVtxAndMult, cfgGeneral.cfgMixEventNumMin, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, v0s1, collision2, tracks2] = *it;
      if (!collision1.sel8() || !collision2.sel8()) {
        continue;
      }

      auto tracks1 = tracks.sliceByCached(o2::aod::track::collisionId, collision1.globalIndex(), this->cache);

      float cent2 = cfgEvSel.cfgCentTableUnavailable ? -1.f : getCentrality(collision2);
      float cent1 = cfgEvSel.cfgCentTableUnavailable ? -1.f : getCentrality(collision1);

      if (cfgEvSel.cfgUseAdditionalEventCut && !selectionEvent(collision1, tracks1.size(), cent1, false)) {
        continue;
      }
      if (cfgEvSel.cfgUseAdditionalEventCut && !selectionEvent(collision2, tracks2.size(), cent2, false)) {
        continue;
      }

      if (!cfgEvSel.cfgSelCollByNch && !cfgEvSel.cfgCentTableUnavailable && (cent1 < cfgEvSel.cfgCutCentMin || cent1 >= cfgEvSel.cfgCutCentMax)) {
        continue;
      }

      if (!cfgEvSel.cfgSelCollByNch && !cfgEvSel.cfgCentTableUnavailable && (cent2 < cfgEvSel.cfgCutCentMin || cent2 >= cfgEvSel.cfgCutCentMax)) {
        continue;
      }

      auto multiplicity1 = static_cast<double>(tracks1.size());
      auto multiplicity2 = static_cast<double>(tracks2.size());

      if (cfgGeneral.cfgStrictTrackCounter) {
        nchCounter(tracks1, multiplicity1);

        if (cfgEvSel.cfgUseCol2ForMixedEvSel) {
          nchCounter(tracks2, multiplicity2);
        }
      }

      if (cfgEvSel.cfgSelCollByNch && (multiplicity1 < cfgEvSel.cfgCutMultMin || multiplicity1 >= cfgEvSel.cfgCutMultMax)) {
        continue;
      }

      if (cfgEvSel.cfgSelCollByNch && (multiplicity2 < cfgEvSel.cfgCutMultMin || multiplicity2 >= cfgEvSel.cfgCutMultMax)) {
        continue;
      }

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      loadCorrection(bc.timestamp());
      float eventWeight = 1.0f;
      if (cfgEffWeights.cfgUseEventWeights) {
        eventWeight = 1.0f / it.currentWindowNeighbours();
      }

      histos.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin
      auto eventClass = (cfgEvSel.cfgSelCollByNch) ? multiplicity1 : cent1;
      float weightCent = 1.0f;
      if (!cfgEvSel.cfgCentTableUnavailable) {
        getCentralityWeight(weightCent, cent1);
      }

      auto assocs = asscTracks->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);

      fillCorrelationsReso(v0s1, assocs, collision1.posZ(), collision1.posY(), collision1.posX(), MixedEvent, getMagneticField(bc.timestamp()), eventClass, eventWeight * weightCent);
    }
  }
  PROCESS_SWITCH(PidDiHadron, processMixedReso, "Process mixed events", true);

  void processMC(FilteredCollisions::iterator const& collision, FilteredMcTracks const& tracksmc, aod::BCsWithTimestamps const&, aod::McParticles const&)
  {
    float cent = -1.;
    if (!cfgEvSel.cfgCentTableUnavailable) {
      cent = getCentrality(collision);
    }
    if (cfgEvSel.cfgUseAdditionalEventCut && !selectionEvent(collision, tracksmc.size(), cent, true)) {
      return;
    }

    if (cfgEvSel.cfgSelCollByNch && (tracksmc.size() < cfgEvSel.cfgCutMultMin || tracksmc.size() >= cfgEvSel.cfgCutMultMax)) {
      return;
    }
    if (!cfgEvSel.cfgSelCollByNch && !cfgEvSel.cfgCentTableUnavailable && (cent < cfgEvSel.cfgCutCentMin || cent >= cfgEvSel.cfgCutCentMax)) {
      return;
    }

    // loop over all tracks
    for (auto const& track : tracksmc) {

      if (!trackSelected(track))
        continue;

      int pidIndex = getNsigmaPID(track);

      // Fill Counts for selection through Nsigma cuts
      if (pidIndex == UsePions) {
        histos.fill(HIST("hNsigmaPionSelected"), track.pt());
      }
      if (pidIndex == UseKaons) {
        histos.fill(HIST("hNsigmaKaonSelected"), track.pt());
      }
      if (pidIndex == UseProtons) {
        histos.fill(HIST("hNsigmaProtonSelected"), track.pt());
      }

      if (track.mcParticle().isPhysicalPrimary()) {
        if (cfgPIDConfigs.cfgPIDParticle == UsePions) {
          histos.fill(HIST("TpcdEdx_ptwise_all_physprim"), track.pt(), track.tpcSignal(), track.tofNSigmaPi());
        }
        if (cfgPIDConfigs.cfgPIDParticle == UseKaons) {
          histos.fill(HIST("TpcdEdx_ptwise_all_physprim"), track.pt(), track.tpcSignal(), track.tofNSigmaKa());
        }
        if (cfgPIDConfigs.cfgPIDParticle == UseProtons) {
          histos.fill(HIST("TpcdEdx_ptwise_all_physprim"), track.pt(), track.tpcSignal(), track.tofNSigmaPr());
        }
      }

      // Check the PDG code for the particles (MC truth) and match with analysed Nsigma PID
      if (std::abs(track.mcParticle().pdgCode()) == PDG_t::kPiPlus) {
        histos.fill(HIST("hNsigmaPionTrue"), track.pt());
        histos.fill(HIST("TpcdEdx_ptwise_pi"), track.pt(), track.tpcSignal(), track.tofNSigmaPi());
        if (track.mcParticle().isPhysicalPrimary()) {
          histos.fill(HIST("TpcdEdx_ptwise_pi_physprim"), track.pt(), track.tpcSignal(), track.tofNSigmaPi());
        }

        double tpcExpSignalPi = track.tpcSignal() - (track.tpcNSigmaPi() * track.tpcExpSigmaPi());
        histos.fill(HIST("ExpTpcdEdx_ptwise_pi"), track.pt(), tpcExpSignalPi, track.tofNSigmaPi());

        if (pidIndex == UsePions) {
          histos.fill(HIST("hNsigmaPionTruePositives"), track.pt());
        }
      } // Pion condition

      if (std::abs(track.mcParticle().pdgCode()) == PDG_t::kKPlus) {
        histos.fill(HIST("hNsigmaKaonTrue"), track.pt());
        histos.fill(HIST("TpcdEdx_ptwise_ka"), track.pt(), track.tpcSignal(), track.tofNSigmaKa());
        if (track.mcParticle().isPhysicalPrimary()) {
          histos.fill(HIST("TpcdEdx_ptwise_ka_physprim"), track.pt(), track.tpcSignal(), track.tofNSigmaKa());
        }

        double tpcExpSignalKa = track.tpcSignal() - (track.tpcNSigmaKa() * track.tpcExpSigmaKa());
        histos.fill(HIST("ExpTpcdEdx_ptwise_ka"), track.pt(), tpcExpSignalKa, track.tofNSigmaKa());

        if (pidIndex == UseKaons) {
          histos.fill(HIST("hNsigmaKaonTruePositives"), track.pt());
        }
      } // Kaon condition

      if (std::abs(track.mcParticle().pdgCode()) == PDG_t::kProton) {
        histos.fill(HIST("hNsigmaProtonTrue"), track.pt());
        histos.fill(HIST("TpcdEdx_ptwise_pr"), track.pt(), track.tpcSignal(), track.tofNSigmaPr());
        if (track.mcParticle().isPhysicalPrimary()) {
          histos.fill(HIST("TpcdEdx_ptwise_pr_physprim"), track.pt(), track.tpcSignal(), track.tofNSigmaPr());
        }

        double tpcExpSignalPr = track.tpcSignal() - (track.tpcNSigmaPr() * track.tpcExpSigmaPr());
        histos.fill(HIST("ExpTpcdEdx_ptwise_pr"), track.pt(), tpcExpSignalPr, track.tofNSigmaPr());

        if (pidIndex == UseProtons) {
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
