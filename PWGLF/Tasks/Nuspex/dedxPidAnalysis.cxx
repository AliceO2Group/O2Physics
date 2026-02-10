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
///
/// \author Paola Vargas Torres  (paola.vargas.torres@cern.ch)
/// \since January 8, 2025
/// \file dedxPidAnalysis.cxx
/// \brief  Analysis to do PID

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "TF1.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::physics;

using PIDTracks = soa::Join<
  aod::Tracks, aod::TracksExtra, aod::TrackSelectionExtension, aod::TracksDCA, aod::TrackSelection,
  aod::pidTOFFullPi, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFbeta, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCEl>;

using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, o2::aod::BarrelMults>;
using BCsRun3 = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

static constexpr int NCentHists{10};
std::array<std::shared_ptr<TH3>, NCentHists> hDedxVsMomentumVsCentPos{};
std::array<std::shared_ptr<TH3>, NCentHists> hDedxVsMomentumVsCentNeg{};

struct DedxPidAnalysis {

  // dE/dx for all charged particles
  HistogramRegistry registryDeDx{
    "registryDeDx",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  // Constant values
  static constexpr int EtaIntervals = 8;
  static constexpr int ParticlesType = 4;
  static constexpr int CentralityClasses = 10;
  float tpcCut = 0.6;
  float pionMin = 0.35;
  float pionMax = 0.45;
  float elTofCut = 0.1;
  float pionTofCut = 1.0;
  float pTcut = 2.0;

  bool fillHist = false;

  enum V0SelectionMode {
    V0TPC = 1,
    V0TOF = 2,
    V0TPCTOF = 3

  };

  enum NINELSelectionMode {
    NoSelINEL = 1,
    SelINELgt0 = 2,
    SelINELgt1 = 3

  };

  enum MomentumMode {
    TpcInnerParam = 1,
    TotalMomentum = 2
  };

  // Event cut labels
  enum EvCutLabel {
    AllEv = 1,
    SelEigth,
    ZVtxCut,
    NoSameBunchPileup,
    GoodZvtxFT0vsPV,
    INELgt
  };

  // Track primary label
  enum TrkPriCutLabel {
    AllPri = 1,
    SelectionPrim,
    PhiVarCutPri,
    NClTPCPIDCutPri,
    NClTPCFoundCutPri
  };

  // Track secondary lebel
  enum TrkSecCutLabel {
    AllSec = 1,
    V0Type,
    V0CosPA,
    V0DecayRadius,
    V0Daughters,
    TPCRefit,
    PhiVarCutSec,
    NClTPCFoundCutSec,
    NClTPCPIDCutSec,
    AllK0s,
    V0RapidityK0s,
    V0ProperLifetimeK0s,
    MassCutK0s,
    AllLambda,
    V0RapidityLambda,
    V0ProperLifetimeLambda,
    MassCutLambda,
    AllAntiLambda,
    V0RapidityAntiLambda,
    V0ProperLifetimeAntiLambda,
    MassCutAntiLambda,
    AllGamma,
    V0RapidityGamma,
    MassCutGamma
  };
  // Configurable Parameters
  // Tracks cuts
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 70.0f,
                                      "min number of found TPC clusters"};
  Configurable<float> minTPCnClsPID{"minTPCnClsPID", 70.0f,
                                    "min number of PID TPC clusters"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.0f, "min number of found TPC crossed rows"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f,
                                 "max chi2 per cluster TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f,
                                 "max chi2 per cluster ITS"};
  Configurable<float> maxZDistanceToIP{"maxZDistanceToIP", 10.0f,
                                       "max z distance to IP"};
  Configurable<float> etaMin{"etaMin", -0.8f, "etaMin"};
  Configurable<float> etaMax{"etaMax", +0.8f, "etaMax"};
  Configurable<float> minNCrossedRowsOverFindableClustersTPC{"minNCrossedRowsOverFindableClustersTPC", 0.8f, "Additional cut on the minimum value of the ratio between crossed rows and findable clusters in the TPC"};
  Configurable<float> maxDCAz{"maxDCAz", 0.1f, "maxDCAz"};
  // v0 cuts
  Configurable<float> v0cospaMin{"v0cospaMin", 0.998f, "Minimum V0 CosPA"};
  Configurable<float> minimumV0Radius{"minimumV0Radius", 0.5f,
                                      "Minimum V0 Radius"};
  Configurable<float> maximumV0Radius{"maximumV0Radius", 100.0f,
                                      "Maximum V0 Radius"};
  Configurable<float> dcaV0DaughtersMax{"dcaV0DaughtersMax", 0.5f,
                                        "Maximum DCA Daughters"};
  Configurable<float> v0rapidityCut{"v0rapidityCut", 0.5f, "V0 rapidity cut"};
  Configurable<float> v0ProperLifetimeCutK0s{"v0ProperLifetimeCutK0s", 20.f, "V0 proper lifetime cut for K0s"};
  Configurable<float> v0ProperLifetimeCutLambda{"v0ProperLifetimeCutLambda", 30.f, "V0 proper lifetime cut for Lambda"};
  Configurable<float> nsigmaTOFmax{"nsigmaTOFmax", 3.0f, "Maximum nsigma TOF"};
  Configurable<float> invMassCutK0s{"invMassCutK0s", 0.2f, "invariant Mass Cut for K0s"};
  Configurable<float> invMassCutLambda{"invMassCutLambda", 0.1f, "invariant Mass Cut for Lambda"};
  Configurable<float> invMassCutGamma{"invMassCutGamma", 0.1f, "invariant Mass Cut for Gamma"};
  Configurable<bool> calibrationMode{"calibrationMode", false, "calibration mode"};
  Configurable<bool> phiVarCut{"phiVarCut", true, "phi var cut"};
  Configurable<bool> nClTPCFoundCut{"nClTPCFoundCut", false, "number of found clusters in TPC cut"};
  Configurable<bool> nClTPCPIDCut{"nClTPCPIDCut", true, "number of PID clusters in TPC cut"};
  Configurable<bool> nGoodZvtx{"nGoodZvtx", true, "Rejects events with no vertex match between FT0 and PV"}; //
  Configurable<bool> nPileUp{"nPileUp", true, "Rejects events with pileup in the same bunch crossing"};
  Configurable<int> nINELSelectionMode{"nINELSelectionMode", 2, "INEL event selection: 1 no sel, 2 INEL>0, 3 INEL>1"};
  Configurable<int> v0SelectionMode{"v0SelectionMode", 3, "V0 Selection base on TPC: 1, TOF:2 ,Both:3"};
  Configurable<int> momentumMode{"momentumMode", 2, "1: TPC inner param, 2: Total momentum p"};
  Configurable<uint8_t> v0TypeSelection{"v0TypeSelection", 1, "select on a certain V0 type (leave negative if no selection desired)"};
  Configurable<double> lowParam1{"lowParam1", 0.119297, "First parameter for low phi cut"};
  Configurable<double> lowParam2{"lowParam2", 0.000379693, "Second parameter for low phi cut"};
  Configurable<double> highParam1{"highParam1", 0.16685, "First parameter for high phi cut"};
  Configurable<double> highParam2{"highParam2", 0.00981942, "Second parameter for high phi cut"};
  // Histograms names
  static constexpr std::string_view DedxvsMomentumPos[ParticlesType] = {"dEdx_vs_Momentum_all_Pos", "dEdx_vs_Momentum_Pi_v0_Pos", "dEdx_vs_Momentum_Pr_v0_Pos", "dEdx_vs_Momentum_El_v0_Pos"};
  static constexpr std::string_view DedxvsMomentumNeg[ParticlesType] = {"dEdx_vs_Momentum_all_Neg", "dEdx_vs_Momentum_Pi_v0_Neg", "dEdx_vs_Momentum_Pr_v0_Neg", "dEdx_vs_Momentum_El_v0_Neg"};
  static constexpr std::string_view DedxvsMomentumvsCentPos[CentralityClasses] = {"dEdx_vs_Momentum_Cent0_1_Pos", "dEdx_vs_Momentum_Cent1_5_Pos", "dEdx_vs_Momentum_Cent5_10_Pos", "dEdx_vs_Momentum_Cent10_15_Pos", "dEdx_vs_Momentum_Cent15_20_Pos", "dEdx_vs_Momentum_Cent20_30_Pos", "dEdx_vs_Momentum_Cent30_40_Pos", "dEdx_vs_Momentum_Cent40_50_Pos", "dEdx_vs_Momentum_Cent50_70_Pos", "dEdx_vs_Momentum_Cent70_100_Pos"};
  static constexpr std::string_view DedxvsMomentumvsCentNeg[CentralityClasses] = {"dEdx_vs_Momentum_Cent0_1_Neg", "dEdx_vs_Momentum_Cent1_5_Neg", "dEdx_vs_Momentum_Cent5_10_Neg", "dEdx_vs_Momentum_Cent10_15_Neg", "dEdx_vs_Momentum_Cent15_20_Neg", "dEdx_vs_Momentum_Cent20_30_Neg", "dEdx_vs_Momentum_Cent30_40_Neg", "dEdx_vs_Momentum_Cent40_50_Neg", "dEdx_vs_Momentum_Cent50_70_Neg", "dEdx_vs_Momentum_Cent70_100_Neg"};
  // Ncl TPC
  static constexpr std::string_view NclTPCDedxMomentumNegBefore[EtaIntervals] = {"Ncl_FoundTPC_vs_dEdx_vs_Momentum_Neg_1_Before", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Neg_2_Before", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Neg_3_Before", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Neg_4_Before", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Neg_5_Before", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Neg_6_Before", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Neg_7_Before", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Neg_8_Before"};
  static constexpr std::string_view NclTPCDedxMomentumPosBefore[EtaIntervals] = {"Ncl_FoundTPC_vs_dEdx_vs_Momentum_Pos_1_Before", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Pos_2_Before", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Pos_3_Before", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Pos_4_Before", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Pos_5_Before", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Pos_6_Before", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Pos_7_Before", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Pos_8_Before"};
  static constexpr std::string_view NclTPCDedxMomentumNegAfter[EtaIntervals] = {"Ncl_TFoundPC_vs_dEdx_vs_Momentum_Neg_1_After", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Neg_2_After", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Neg_3_After", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Neg_4_After", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Neg_5_After", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Neg_6_After", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Neg_7_After", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Neg_8_After"};
  static constexpr std::string_view NclTPCDedxMomentumPosAfter[EtaIntervals] = {"Ncl_FoundTPC_vs_dEdx_vs_Momentum_Pos_1_After", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Pos_2_After", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Pos_3_After", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Pos_4_After", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Pos_5_After", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Pos_6_After", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Pos_7_After", "Ncl_FoundTPC_vs_dEdx_vs_Momentum_Pos_8_After"};
  // Ncl PID TPC
  static constexpr std::string_view NclPIDTPCDedxMomentumNegBefore[EtaIntervals] = {"Ncl_PIDTPC_vs_dEdx_vs_Momentum_Neg_1_Before", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Neg_2_Before", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Neg_3_Before", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Neg_4_Before", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Neg_5_Before", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Neg_6_Before", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Neg_7_Before", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Neg_8_Before"};
  static constexpr std::string_view NclPIDTPCDedxMomentumPosBefore[EtaIntervals] = {"Ncl_PIDTPC_vs_dEdx_vs_Momentum_Pos_1_Before", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Pos_2_Before", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Pos_3_Before", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Pos_4_Before", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Pos_5_Before", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Pos_6_Before", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Pos_7_Before", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Pos_8_Before"};
  static constexpr std::string_view NclPIDTPCDedxMomentumNegAfter[EtaIntervals] = {"Ncl_PIDTPC_vs_dEdx_vs_Momentum_Neg_1_After", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Neg_2_After", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Neg_3_After", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Neg_4_After", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Neg_5_After", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Neg_6_After", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Neg_7_After", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Neg_8_After"};
  static constexpr std::string_view NclPIDTPCDedxMomentumPosAfter[EtaIntervals] = {"Ncl_PIDTPC_vs_dEdx_vs_Momentum_Pos_1_After", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Pos_2_After", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Pos_3_After", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Pos_4_After", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Pos_5_After", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Pos_6_After", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Pos_7_After", "Ncl_PIDTPC_vs_dEdx_vs_Momentum_Pos_8_After"};
  static constexpr double EtaCut[EtaIntervals + 1] = {-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8};
  static constexpr double CentClasses[CentralityClasses + 1] = {0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};
  Configurable<std::vector<float>> calibrationFactorNeg{"calibrationFactorNeg", {50.4011, 50.4764, 50.186, 49.2955, 48.8222, 49.4273, 49.9292, 50.0556}, "negative calibration factors"};
  Configurable<std::vector<float>> calibrationFactorPos{"calibrationFactorPos", {50.5157, 50.6359, 50.3198, 49.3345, 48.9197, 49.4931, 50.0188, 50.1406}, "positive calibration factors"};
  ConfigurableAxis binP{"binP", {VARIABLE_WIDTH, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0}, ""};

  // phi cut fits
  TF1* fphiCutHigh = nullptr;
  TF1* fphiCutLow = nullptr;
  Service<ccdb::BasicCCDBManager> ccdb;

  TrackSelection myTrackSelection()
  {
    TrackSelection selectedTracks;
    selectedTracks.SetPtRange(0.1f, 1e10f);
    selectedTracks.SetEtaRange(etaMin, etaMax);
    selectedTracks.SetRequireITSRefit(true);
    selectedTracks.SetRequireTPCRefit(true);
    selectedTracks.SetMinNCrossedRowsTPC(minNCrossedRowsTPC);
    selectedTracks.SetMinNCrossedRowsOverFindableClustersTPC(minNCrossedRowsOverFindableClustersTPC);
    selectedTracks.SetMaxChi2PerClusterTPC(maxChi2TPC);
    selectedTracks.SetRequireHitsInITSLayers(1, {0, 1, 2});
    selectedTracks.SetMaxChi2PerClusterITS(maxChi2ITS);
    selectedTracks.SetMaxDcaXYPtDep([](float pt) { return 0.0105f + 0.0350f / std::pow(pt, 1.1f); });
    selectedTracks.SetMaxDcaZ(maxDCAz);
    selectedTracks.SetRequireGoldenChi2(true);

    return selectedTracks;
  }

  TrackSelection mySelectionPrim;

  void init(InitContext const&)
  {
    const char* label = "No INEL selection";
    if (nPileUp) {
      LOGF(info, "Applying NoSameBunchPileup cut");
    } else {
      LOGF(info, "NoSameBunchPileup cut disabled");
    }
    if (nGoodZvtx) {
      LOGF(info, "Applying GoodZvtxFT0vsPV cut");
    } else {
      LOGF(info, "GoodZvtxFT0vsPV cut disabled");
    }
    if (nINELSelectionMode == NoSelINEL) {
      LOGF(info, "INEL cut disabled");
    } else if (nINELSelectionMode == SelINELgt0) {
      LOGF(info, "Applying INEL > 0 cut");
      label = "INEL > 0";
    } else if (nINELSelectionMode == SelINELgt1) {
      LOGF(info, "Applying INEL > 1 cut");
      label = "INEL > 1";
    }
    if (v0SelectionMode == V0TPC) {
      LOGF(info, "V0 seleccion using TPC only");
    } else if (v0SelectionMode == V0TOF) {
      LOGF(info, "V0 seleccion using TOF only");
    } else if (v0SelectionMode == V0TPCTOF) {
      LOGF(info, "V0 seleccion using TOF + TPC");
    }
    if (calibrationMode) {
      LOGF(info, "Calibration mode activated");
    } else {
      LOGF(info, "Running data applying calibration values");
    }

    if (phiVarCut) {
      LOGF(info, "Applying phi prime cut");
    } else {
      LOGF(info, "Phi prime cut disabled");
    }

    if (nClTPCFoundCut) {
      LOGF(info, "Applying TPC found clusters cut");
    } else {
      LOGF(info, "TPC found clusters cut disabled");
    }

    if (nClTPCPIDCut) {
      LOGF(info, "Applying TPC clusters for PID cut");
    } else {
      LOGF(info, "TPC clusters for PID cut disabled");
    }
    if (momentumMode == TpcInnerParam) {
      LOGF(info, "Using TPC inner parameter for momentum calculation");
    } else {
      LOGF(info, "Using total momentum (p) for calculation");
    }
    LOGF(info, "Centrality clases between %.1f - %.1f", CentClasses[0], CentClasses[10]);

    AxisSpec dedxAxis{100, 0.0, 100.0, "dE/dx (a. u.)"};
    AxisSpec ptAxis = {binP, "pT (GeV/c)"};
    AxisSpec etaAxis{8, -0.8, 0.8, "#eta"};
    AxisSpec pAxis = {binP, "#it{p}/Z (GeV/c)"};
    AxisSpec pAxisTrack = {binP, "#it{p} (GeV/c)"};
    fphiCutLow = new TF1("StandardPhiCutLow",
                         Form("%f/x/x+pi/18.0-%f", lowParam1.value, lowParam2.value),
                         0, 50);
    fphiCutHigh = new TF1("StandardPhiCutHigh",
                          Form("%f/x+pi/18.0+%f", highParam1.value, highParam2.value),
                          0, 50);
    LOGF(info, "=== Phi Cut Parameters ===");
    LOGF(info, "Low cut: %.6f/x² + pi/18 - %.6f", lowParam1.value, lowParam2.value);
    LOGF(info, "High cut: %.6f/x + pi/18 + %.6f", highParam1.value, highParam2.value);

    if (calibrationMode) {
      // MIP for pions
      registryDeDx.add(
        "hdEdx_vs_eta_Neg_Pi", "dE/dx", HistType::kTH2F,
        {{etaAxis}, {dedxAxis}});
      registryDeDx.add(
        "hdEdx_vs_eta_Pos_Pi", "dE/dx", HistType::kTH2F,
        {{etaAxis}, {dedxAxis}});
      // MIP for electrons
      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Neg_El", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});
      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Pos_El", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});
      // Pions from TOF
      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Neg_TOF", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});
      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Pos_TOF", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});

    } else {
      // MIP for pions
      registryDeDx.add(
        "hdEdx_vs_eta_Neg_calibrated_Pi", "dE/dx", HistType::kTH2F,
        {{etaAxis}, {dedxAxis}});

      registryDeDx.add(
        "hdEdx_vs_eta_Pos_calibrated_Pi", "dE/dx", HistType::kTH2F,
        {{etaAxis}, {dedxAxis}});

      // MIP for electrons
      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Neg_calibrated_El", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});

      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Pos_calibrated_El", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});

      // Pions from TOF
      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Neg_calibrated_TOF", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});

      registryDeDx.add(
        "hdEdx_vs_eta_vs_p_Pos_calibrated_TOF", "dE/dx", HistType::kTH3F,
        {{etaAxis}, {dedxAxis}, {pAxis}});

      // pt vs p
      registryDeDx.add(
        "heta_vs_pt_vs_p_all_Neg", "eta_vs_pT_vs_p", HistType::kTH3F,
        {{etaAxis}, {ptAxis}, {pAxisTrack}});
      registryDeDx.add(
        "heta_vs_pt_vs_p_all_Pos", "eta_vs_pT_vs_p", HistType::kTH3F,
        {{etaAxis}, {ptAxis}, {pAxisTrack}});

      // De/Dx for ch and v0 particles
      for (int i = 0; i < ParticlesType; ++i) {
        registryDeDx.add(DedxvsMomentumPos[i].data(), "dE/dx", HistType::kTH3F,
                         {{pAxisTrack}, {dedxAxis}, {etaAxis}});
        registryDeDx.add(DedxvsMomentumNeg[i].data(), "dE/dx", HistType::kTH3F,
                         {{pAxisTrack}, {dedxAxis}, {etaAxis}});
      }

      for (int i = 0; i < CentralityClasses; ++i) {
        hDedxVsMomentumVsCentPos[i] = registryDeDx.add<TH3>(DedxvsMomentumvsCentPos[i].data(), "dE/dx", HistType::kTH3F, {{pAxisTrack}, {dedxAxis}, {etaAxis}});
        hDedxVsMomentumVsCentNeg[i] = registryDeDx.add<TH3>(DedxvsMomentumvsCentNeg[i].data(), "dE/dx", HistType::kTH3F, {{pAxisTrack}, {dedxAxis}, {etaAxis}});
      }
    }

    registryDeDx.add(
      "hdEdx_vs_phi", "dE/dx", HistType::kTH2F,
      {{100, 0.0, 6.4, "#phi"}, {dedxAxis}});
    // phi cut
    if (phiVarCut) {
      // pt for found
      registryDeDx.add(
        "hpt_vs_phi_NclFound_TPC_After", "phi cut vs pt", HistType::kTH3F,
        {{ptAxis}, {100, 0.0, 0.4, "#varphi^{'}"}, {100, 0, 160, "N_{cl, found}"}});

      registryDeDx.add(
        "hpt_vs_phi_NclFound_TPC_Before", "phi cut vs pt", HistType::kTH3F,
        {{ptAxis}, {100, 0.0, 0.4, "#varphi^{'}"}, {100, 0, 160, "N_{cl, found}"}});
      // p
      registryDeDx.add(
        "hp_vs_phi_NclFound_TPC_After", "phi cut vs p", HistType::kTH3F,
        {{pAxis}, {100, 0.0, 0.4, "#varphi^{'}"}, {100, 0, 160, "N_{cl, found}"}});

      registryDeDx.add(
        "hp_vs_phi_NclFound_TPC_Before", "phi cut vs p", HistType::kTH3F,
        {{pAxis}, {100, 0.0, 0.4, "#varphi^{'}"}, {100, 0, 160, "N_{cl, found}"}});

      // pt for PID
      registryDeDx.add(
        "hpt_vs_phi_NclPID_TPC_After", "phi cut vs pt", HistType::kTH3F,
        {{ptAxis}, {100, 0.0, 0.4, "#varphi^{'}"}, {100, 0, 160, "N_{cl, PID}"}});

      registryDeDx.add(
        "hpt_vs_phi_NclPID_TPC_Before", "phi cut vs pt", HistType::kTH3F,
        {{ptAxis}, {100, 0.0, 0.4, "#varphi^{'}"}, {100, 0, 160, "N_{cl, PID}"}});
      // p
      registryDeDx.add(
        "hp_vs_phi_NclPID_TPC_After", "phi cut vs p", HistType::kTH3F,
        {{pAxis}, {100, 0.0, 0.4, "#varphi^{'}"}, {100, 0, 160, "N_{cl, PID}"}});

      registryDeDx.add(
        "hp_vs_phi_NclPID_TPC_Before", "phi cut vs p", HistType::kTH3F,
        {{pAxis}, {100, 0.0, 0.4, "#varphi^{'}"}, {100, 0, 160, "N_{cl, PID}"}});
    }
    // Ncl vs de/dx TPC
    if (nClTPCFoundCut) {
      for (int i = 0; i < EtaIntervals; ++i) {
        registryDeDx.add(NclTPCDedxMomentumPosBefore[i].data(), "Ncl found TPC vs dE/dx vs Momentum Positive before", HistType::kTH3F,
                         {{100, 0, 160, "N_{cl, found}^{TPC}"}, {dedxAxis}, {pAxis}});
        registryDeDx.add(NclTPCDedxMomentumNegBefore[i].data(), "Ncl found TPC vs dE/dx vs Momentum Negative before", HistType::kTH3F,
                         {{100, 0, 160, "N_{cl, found}^{TPC}"}, {dedxAxis}, {pAxis}});

        registryDeDx.add(NclTPCDedxMomentumPosAfter[i].data(), "Ncl found TPC vs dE/dx vs Momentum Positive after", HistType::kTH3F,
                         {{100, 0, 160, "N_{cl, found}^{TPC}"}, {dedxAxis}, {pAxis}});
        registryDeDx.add(NclTPCDedxMomentumNegAfter[i].data(), "Ncl found TPC vs dE/dx vs Momentum Negative after", HistType::kTH3F,
                         {{100, 0, 160, "N_{cl, found}^{TPC}"}, {dedxAxis}, {pAxis}});
      }
    }

    // Ncl vs de/dx ITS
    if (nClTPCPIDCut) {
      for (int i = 0; i < EtaIntervals; ++i) {
        registryDeDx.add(NclPIDTPCDedxMomentumPosBefore[i].data(), "Ncl PID TPC vs dE/dx vs Momentum Positive before", HistType::kTH3F,
                         {{100, 0, 160, "N_{cl, PID}^{TPC}"}, {dedxAxis}, {pAxis}});
        registryDeDx.add(NclPIDTPCDedxMomentumNegBefore[i].data(), "Ncl PID TPC vs dE/dx vs Momentum Negative before", HistType::kTH3F,
                         {{100, 0, 160, "N_{cl, PID}^{TPC}"}, {dedxAxis}, {pAxis}});

        registryDeDx.add(NclPIDTPCDedxMomentumPosAfter[i].data(), "Ncl PID TPC vs dE/dx vs Momentum Positive after", HistType::kTH3F,
                         {{100, 0, 160, "N_{cl, PID}^{TPC}"}, {dedxAxis}, {pAxis}});
        registryDeDx.add(NclPIDTPCDedxMomentumNegAfter[i].data(), "Ncl PID TPC vs dE/dx vs Momentum Negative after", HistType::kTH3F,
                         {{100, 0, 160, "N_{cl, PID}^{TPC}"}, {dedxAxis}, {pAxis}});
      }
    }
    // eta
    registryDeDx.add(
      "heta_vs_NclFound_TPC_Before_Primary", "eta and N_{cl}", HistType::kTH2F,
      {{100, -0.8, 0.8, "#eta"}, {100, 0, 160, "N_{cl, found}"}});

    registryDeDx.add(
      "heta_vs_NclFound_TPC_After_Primary", "eta and N_{cl} ", HistType::kTH2F,
      {{100, -0.8, 0.8, "#eta"}, {100, 0, 160, "N_{cl, found}"}});

    registryDeDx.add(
      "heta_vs_NclPID_TPC_Before_Primary", "eta and N_{cl, PID}", HistType::kTH2F,
      {{100, -0.8, 0.8, "#eta"}, {100, 0, 160, "N_{cl, PID}"}});

    registryDeDx.add(
      "heta_vs_NclPID_TPC_After_Primary", "eta and N_{cl, PID} ", HistType::kTH2F,
      {{100, -0.8, 0.8, "#eta"}, {100, 0, 160, "N_{cl, PID}"}});
    // momentum for primaries
    registryDeDx.add(
      "hp_vs_NclPID_TPC_Before_Primary", "p and N_{cl, PID}", HistType::kTH2F,
      {{pAxisTrack}, {100, 0, 160, "N_{cl, PID}"}});

    registryDeDx.add(
      "hp_vs_NclPID_TPC_After_Primary", "p and N_{cl, PID} ", HistType::kTH2F,
      {{pAxisTrack}, {100, 0, 160, "N_{cl, PID}"}});

    // eta for secondaries
    registryDeDx.add(
      "heta_vs_NclPID_TPC_Before_PionsK0s", "eta and N_{cl, PID}", HistType::kTH2F,
      {{100, -0.8, 0.8, "#eta"}, {100, 0, 160, "N_{cl, PID}"}});

    registryDeDx.add(
      "heta_vs_NclPID_TPC_After_PionsK0s", "eta and N_{cl, PID} ", HistType::kTH2F,
      {{100, -0.8, 0.8, "#eta"}, {100, 0, 160, "N_{cl, PID}"}});

    registryDeDx.add(
      "heta_vs_NclPID_TPC_Before_PionsLambda", "eta and N_{cl, PID}", HistType::kTH2F,
      {{100, -0.8, 0.8, "#eta"}, {100, 0, 160, "N_{cl, PID}"}});

    registryDeDx.add(
      "heta_vs_NclPID_TPC_After_PionsLambda", "eta and N_{cl, PID} ", HistType::kTH2F,
      {{100, -0.8, 0.8, "#eta"}, {100, 0, 160, "N_{cl, PID}"}});

    registryDeDx.add(
      "heta_vs_NclPID_TPC_Before_ProtonsLambda", "eta and N_{cl, PID}", HistType::kTH2F,
      {{100, -0.8, 0.8, "#eta"}, {100, 0, 160, "N_{cl, PID}"}});

    registryDeDx.add(
      "heta_vs_NclPID_TPC_After_ProtonsLambda", "eta and N_{cl, PID} ", HistType::kTH2F,
      {{100, -0.8, 0.8, "#eta"}, {100, 0, 160, "N_{cl, PID}"}});

    registryDeDx.add(
      "heta_vs_NclPID_TPC_Before_ElectronsGamma", "eta and N_{cl, PID}", HistType::kTH2F,
      {{100, -0.8, 0.8, "#eta"}, {100, 0, 160, "N_{cl, PID}"}});

    registryDeDx.add(
      "heta_vs_NclPID_TPC_After_ElectronsGamma", "eta and N_{cl, PID} ", HistType::kTH2F,
      {{100, -0.8, 0.8, "#eta"}, {100, 0, 160, "N_{cl, PID}"}});

    // momentum for secondaries
    registryDeDx.add(
      "hp_vs_NclPID_TPC_Before_PionsK0s", "p and N_{cl, PID}", HistType::kTH2F,
      {{pAxisTrack}, {100, 0, 160, "N_{cl, PID}"}});

    registryDeDx.add(
      "hp_vs_NclPID_TPC_After_PionsK0s", "p and N_{cl, PID} ", HistType::kTH2F,
      {{pAxisTrack}, {100, 0, 160, "N_{cl, PID}"}});

    registryDeDx.add(
      "hp_vs_NclPID_TPC_Before_PionsLambda", "p and N_{cl, PID}", HistType::kTH2F,
      {{pAxisTrack}, {100, 0, 160, "N_{cl, PID}"}});

    registryDeDx.add(
      "hp_vs_NclPID_TPC_After_PionsLambda", "p and N_{cl, PID} ", HistType::kTH2F,
      {{pAxisTrack}, {100, 0, 160, "N_{cl, PID}"}});

    registryDeDx.add(
      "hp_vs_NclPID_TPC_Before_ProtonsLambda", "p and N_{cl, PID}", HistType::kTH2F,
      {{pAxisTrack}, {100, 0, 160, "N_{cl, PID}"}});

    registryDeDx.add(
      "hp_vs_NclPID_TPC_After_ProtonsLambda", "p and N_{cl, PID} ", HistType::kTH2F,
      {{pAxisTrack}, {100, 0, 160, "N_{cl, PID}"}});

    registryDeDx.add(
      "hp_vs_NclPID_TPC_Before_ElectronsGamma", "p and N_{cl, PID}", HistType::kTH2F,
      {{pAxisTrack}, {100, 0, 160, "N_{cl, PID}"}});

    registryDeDx.add(
      "hp_vs_NclPID_TPC_After_ElectronsGamma", "p and N_{cl, PID} ", HistType::kTH2F,
      {{pAxisTrack}, {100, 0, 160, "N_{cl, PID}"}});

    // beta plot
    registryDeDx.add(
      "hbeta_vs_p_Neg", "beta", HistType::kTH2F,
      {{pAxis}, {100, 0.0, 1.1, "#beta"}});

    registryDeDx.add(
      "hbeta_vs_p_Pos", "beta", HistType::kTH2F,
      {{pAxis}, {100, 0.0, 1.1, "#beta"}});

    // Event Counter
    registryDeDx.add("histRecVtxZData", "collision z position", HistType::kTH1F, {{100, -20.0, +20.0, "z_{vtx} (cm)"}});

    // Event Counter by centrality
    registryDeDx.add("histCentrality", "collision centrality", HistType::kTH1F, {{100, 0.0, 100, "cent"}});

    // Event Counter
    registryDeDx.add("evsel", "events selected", HistType::kTH1F, {{6, 0.5, 6.5, ""}});
    auto hstat = registryDeDx.get<TH1>(HIST("evsel"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(AllEv, "AllEv");
    x->SetBinLabel(SelEigth, "SelEigth");
    x->SetBinLabel(ZVtxCut, "ZVtxCut");
    x->SetBinLabel(NoSameBunchPileup, "NoSameBunchPileup");
    x->SetBinLabel(GoodZvtxFT0vsPV, "GoodZvtxFT0vsPV");
    x->SetBinLabel(INELgt, label);

    // Track Prim Counter
    registryDeDx.add("trackselAll", "track selected all particles", HistType::kTH1F, {{5, 0.5, 5.5, ""}});
    auto htrackAll = registryDeDx.get<TH1>(HIST("trackselAll"));
    auto* xAll = htrackAll->GetXaxis();
    xAll->SetBinLabel(AllPri, "AllPri");
    xAll->SetBinLabel(SelectionPrim, "SelectionPrim");
    xAll->SetBinLabel(PhiVarCutPri, "PhiVarCutPri");
    xAll->SetBinLabel(NClTPCPIDCutPri, "NClTPCPIDCutPri");
    xAll->SetBinLabel(NClTPCFoundCutPri, "NClTPCFoundCutPri");

    registryDeDx.add("trackselSec", "track selected sec particles", HistType::kTH1F, {{24, 0.5, 24.5, ""}});
    auto htrackSec = registryDeDx.get<TH1>(HIST("trackselSec"));
    auto* xSec = htrackSec->GetXaxis();
    xSec->SetBinLabel(AllSec, "AllSec");
    xSec->SetBinLabel(V0Type, "V0Type");
    xSec->SetBinLabel(V0CosPA, "V0CosPA");
    xSec->SetBinLabel(V0DecayRadius, "V0DecayRadius");
    xSec->SetBinLabel(V0Daughters, "V0Daughters");
    xSec->SetBinLabel(TPCRefit, "TPCRefit");
    xSec->SetBinLabel(PhiVarCutSec, "PhiVarCutSec");
    xSec->SetBinLabel(NClTPCFoundCutSec, "NClTPCFoundCutSec");
    xSec->SetBinLabel(NClTPCPIDCutSec, "NClTPCPIDCutSec");
    xSec->SetBinLabel(AllK0s, "AllK0s");
    xSec->SetBinLabel(V0RapidityK0s, "V0RapidityK0s");
    xSec->SetBinLabel(V0ProperLifetimeK0s, "V0ProperLifetimeK0s");
    xSec->SetBinLabel(MassCutK0s, "MassCutK0s");
    xSec->SetBinLabel(AllLambda, "AllLambda");
    xSec->SetBinLabel(V0RapidityLambda, "V0RapidityLambda");
    xSec->SetBinLabel(V0ProperLifetimeLambda, "V0ProperLifetimeLambda");
    xSec->SetBinLabel(MassCutLambda, "MassCutLambda");
    xSec->SetBinLabel(AllAntiLambda, "AllAntiLambda");
    xSec->SetBinLabel(V0RapidityAntiLambda, "V0RapidityAntiLambda");
    xSec->SetBinLabel(V0ProperLifetimeAntiLambda, "V0ProperLifetimeAntiLambda");
    xSec->SetBinLabel(MassCutAntiLambda, "MassCutAntiLambda");
    xSec->SetBinLabel(AllGamma, "AllGamma");
    xSec->SetBinLabel(V0RapidityGamma, "V0RapidityGamma");
    xSec->SetBinLabel(MassCutGamma, "MassCutGamma");

    mySelectionPrim = myTrackSelection();
  }

  // Single-Track Selection
  template <typename T1, typename C>
  bool passedSingleTrackSelection(const T1& track, const C& /*collision*/)
  {
    // Single-Track Selections
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if (track.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (track.eta() < etaMin || track.eta() > etaMax)
      return false;

    return true;
  }

  // Momentum
  template <typename T1>
  float getMomentum(const T1& track)
  {
    return (momentumMode == TpcInnerParam) ? track.tpcInnerParam() : track.p();
  };

  // General V0 Selections
  template <typename T1, typename C>
  bool passedV0Selection(const T1& v0, const C& /*collision*/)
  {
    if (v0.v0cosPA() < v0cospaMin)
      return false;
    registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::V0CosPA);

    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;
    registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::V0DecayRadius);

    if (v0.dcaV0daughters() > dcaV0DaughtersMax)
      return false;
    registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::V0Daughters);

    return true;
  }

  // K0s Selections
  template <typename T1, typename T2, typename C>
  bool passedK0Selection(const T1& v0, const T2& ntrack, const T2& ptrack,
                         const C& collision)
  {
    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack, collision))
      return false;
    if (!passedSingleTrackSelection(ntrack, collision))
      return false;
    double sigmap = 0.0;
    double sigman = 0.0;

    if (v0SelectionMode == V0TPC) {
      sigmap = ptrack.tpcNSigmaPi();
      sigman = ntrack.tpcNSigmaPi();
    } else if (v0SelectionMode == V0TOF) {
      sigmap = ptrack.tofNSigmaPi();
      sigman = ntrack.tofNSigmaPi();
    } else if (v0SelectionMode == V0TPCTOF) {
      sigmap = std::hypot(ptrack.tpcNSigmaPi(), ptrack.tofNSigmaPi());
      sigman = std::hypot(ntrack.tpcNSigmaPi(), ntrack.tofNSigmaPi());
    }

    if (ptrack.tpcInnerParam() > tpcCut) {
      if (!ptrack.hasTOF())
        return false;
      if (std::abs(sigmap) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > tpcCut) {
      if (!ntrack.hasTOF())
        return false;
      if (std::abs(sigman) > nsigmaTOFmax)
        return false;
    }
    if (fillHist)
      registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::AllK0s);

    if (std::abs(v0.yK0Short()) > v0rapidityCut)
      return false;
    if (fillHist)
      registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::V0RapidityK0s);

    float properLifetime = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * MassK0Short;

    if (properLifetime > v0ProperLifetimeCutK0s)
      return false;
    if (fillHist)
      registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::V0ProperLifetimeK0s);

    if (std::abs(v0.mK0Short() - MassK0Short) > invMassCutK0s)
      return false;
    if (fillHist)
      registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::MassCutK0s);

    return true;
  }

  // Lambda Selections
  template <typename T1, typename T2, typename C>
  bool passedLambdaSelection(const T1& v0, const T2& ntrack, const T2& ptrack,
                             const C& collision)
  {
    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack, collision))
      return false;
    if (!passedSingleTrackSelection(ntrack, collision))
      return false;

    double sigmap = 0.0;
    double sigman = 0.0;

    if (v0SelectionMode == V0TPC) {
      sigmap = ptrack.tpcNSigmaPr();
      sigman = ntrack.tpcNSigmaPi();
    } else if (v0SelectionMode == V0TOF) {
      sigmap = ptrack.tofNSigmaPr();
      sigman = ntrack.tofNSigmaPi();
    } else if (v0SelectionMode == V0TPCTOF) {
      sigmap = std::hypot(ptrack.tpcNSigmaPr(), ptrack.tofNSigmaPr());
      sigman = std::hypot(ntrack.tpcNSigmaPi(), ntrack.tofNSigmaPi());
    }

    if (ptrack.tpcInnerParam() > tpcCut) {
      if (!ptrack.hasTOF())
        return false;
      if (std::abs(sigmap) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > tpcCut) {
      if (!ntrack.hasTOF())
        return false;
      if (std::abs(sigman) > nsigmaTOFmax)
        return false;
    }
    if (fillHist)
      registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::AllLambda);

    if (std::abs(v0.yLambda()) > v0rapidityCut)
      return false;
    if (fillHist)
      registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::V0RapidityLambda);

    float properLifetime = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * MassLambda;

    if (properLifetime > v0ProperLifetimeCutLambda)
      return false;

    if (fillHist)
      registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::V0ProperLifetimeLambda);

    if (std::abs(v0.mLambda() - MassLambda) > invMassCutLambda) {
      return false;
    }
    if (fillHist)
      registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::MassCutLambda);

    return true;
  }

  // AntiLambda Selections
  template <typename T1, typename T2, typename C>
  bool passedAntiLambdaSelection(const T1& v0, const T2& ntrack,
                                 const T2& ptrack, const C& collision)
  {

    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack, collision))
      return false;
    if (!passedSingleTrackSelection(ntrack, collision))
      return false;
    double sigmap = 0.0;
    double sigman = 0.0;

    if (v0SelectionMode == V0TPC) {
      sigmap = ptrack.tpcNSigmaPi();
      sigman = ntrack.tpcNSigmaPr();
    } else if (v0SelectionMode == V0TOF) {
      sigmap = ptrack.tofNSigmaPi();
      sigman = ntrack.tofNSigmaPr();
    } else if (v0SelectionMode == V0TPCTOF) {
      sigmap = std::hypot(ptrack.tpcNSigmaPi(), ptrack.tofNSigmaPi());
      sigman = std::hypot(ntrack.tpcNSigmaPr(), ntrack.tofNSigmaPr());
    }
    if (ptrack.tpcInnerParam() > tpcCut) {
      if (!ptrack.hasTOF())
        return false;
      if (std::abs(sigmap) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > tpcCut) {
      if (!ntrack.hasTOF())
        return false;
      if (std::abs(sigman) > nsigmaTOFmax)
        return false;
    }
    if (fillHist)
      registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::AllAntiLambda);

    if (std::abs(v0.yLambda()) > v0rapidityCut)
      return false;

    if (fillHist)
      registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::V0RapidityAntiLambda);

    float properLifetime = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * MassLambda;

    if (properLifetime > v0ProperLifetimeCutLambda)
      return false;

    if (fillHist)
      registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::V0ProperLifetimeAntiLambda);

    if (std::abs(v0.mAntiLambda() - MassLambda) > invMassCutLambda)
      return false;

    if (fillHist)
      registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::MassCutAntiLambda);

    return true;
  }

  // Gamma Selections
  template <typename T1, typename T2, typename C>
  bool passedGammaSelection(const T1& v0, const T2& ntrack, const T2& ptrack,
                            const C& collision)
  {
    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack, collision))
      return false;
    if (!passedSingleTrackSelection(ntrack, collision))
      return false;

    double sigmap = 0.0;
    double sigman = 0.0;

    if (v0SelectionMode == V0TPC) {
      sigmap = ptrack.tpcNSigmaEl();
      sigman = ntrack.tpcNSigmaEl();
    } else if (v0SelectionMode == V0TOF) {
      sigmap = ptrack.tofNSigmaEl();
      sigman = ntrack.tofNSigmaEl();
    } else if (v0SelectionMode == V0TPCTOF) {
      sigmap = std::hypot(ptrack.tpcNSigmaEl(), ptrack.tofNSigmaEl());
      sigman = std::hypot(ntrack.tpcNSigmaEl(), ntrack.tofNSigmaEl());
    }

    if (ptrack.tpcInnerParam() > tpcCut) {
      if (!ptrack.hasTOF())
        return false;
      if (std::abs(sigmap) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > tpcCut) {
      if (!ntrack.hasTOF())
        return false;
      if (std::abs(sigman) > nsigmaTOFmax)
        return false;
    }
    const float gammaMass = 2 * MassElectron; // GeV/c^2

    if (fillHist)
      registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::AllGamma);

    const float yGamma = RecoDecay::y(std::array{v0.px(), v0.py(), v0.pz()}, MassGamma);

    if (std::abs(yGamma) > v0rapidityCut)
      return false;

    if (fillHist)
      registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::V0RapidityGamma);

    if (std::abs(v0.mGamma() - gammaMass) > invMassCutGamma)
      return false;

    if (fillHist)
      registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::MassCutGamma);

    return true;
  }
  // Magnetic field
  int getMagneticField(uint64_t timestamp)
  {
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }
  // Phi cut
  template <typename T>
  bool passedPhiCutPri(const T& trk, float magField, const TF1& fphiCutLow, const TF1& fphiCutHigh)
  {
    float p = trk.p();
    float pt = trk.pt();
    float phi = trk.phi();
    int charge = trk.sign();
    auto nTPCCl = trk.tpcNClsFound();
    auto nTPCPIDCl = trk.tpcNClsPID();

    if (pt < pTcut)
      return true;

    if (magField < 0) // for negatve polarity field
      phi = o2::constants::math::TwoPI - phi;
    if (charge < 0) // for negatve charge
      phi = o2::constants::math::TwoPI - phi;

    // to center gap in the middle
    phi += o2::constants::math::PI / 18.0f;
    phi = std::fmod(phi, o2::constants::math::PI / 9.0f);

    registryDeDx.fill(HIST("hpt_vs_phi_NclFound_TPC_Before"), pt, phi, nTPCCl);
    registryDeDx.fill(HIST("hp_vs_phi_NclFound_TPC_Before"), p, phi, nTPCCl);
    registryDeDx.fill(HIST("hpt_vs_phi_NclPID_TPC_Before"), pt, phi, nTPCPIDCl);
    registryDeDx.fill(HIST("hp_vs_phi_NclPID_TPC_Before"), p, phi, nTPCPIDCl);

    // cut phi
    if (phi < fphiCutHigh.Eval(pt) && phi > fphiCutLow.Eval(pt))
      return false; // reject track

    registryDeDx.fill(HIST("hpt_vs_phi_NclFound_TPC_After"), pt, phi, nTPCCl);
    registryDeDx.fill(HIST("hp_vs_phi_NclFound_TPC_After"), p, phi, nTPCCl);
    registryDeDx.fill(HIST("hpt_vs_phi_NclPID_TPC_After"), pt, phi, nTPCPIDCl);
    registryDeDx.fill(HIST("hp_vs_phi_NclPID_TPC_After"), p, phi, nTPCPIDCl);

    return true;
  }

  // NclCutTPC
  template <typename T>
  bool passedNClTPCFoundCutPri(const T& trk)
  {
    float eta = trk.eta();
    float sigP = trk.sign() * getMomentum(trk);
    auto nTPCCl = trk.tpcNClsFound();

    if (eta >= EtaCut[0] && eta < EtaCut[1]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclTPCDedxMomentumNegBefore[0]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclTPCDedxMomentumPosBefore[0]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[1] && eta < EtaCut[2]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclTPCDedxMomentumNegBefore[1]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclTPCDedxMomentumPosBefore[1]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[2] && eta < EtaCut[3]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclTPCDedxMomentumNegBefore[2]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclTPCDedxMomentumPosBefore[2]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[3] && eta < EtaCut[4]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclTPCDedxMomentumNegBefore[3]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclTPCDedxMomentumPosBefore[3]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[4] && eta < EtaCut[5]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclTPCDedxMomentumNegBefore[4]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclTPCDedxMomentumPosBefore[4]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[5] && eta < EtaCut[6]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclTPCDedxMomentumNegBefore[5]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclTPCDedxMomentumPosBefore[5]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[6] && eta < EtaCut[7]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclTPCDedxMomentumNegBefore[6]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclTPCDedxMomentumPosBefore[6]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[7] && eta < EtaCut[8]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclTPCDedxMomentumNegBefore[7]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclTPCDedxMomentumPosBefore[7]), nTPCCl, trk.tpcSignal(), sigP);
      }
    }

    if (nTPCCl < minTPCnClsFound)
      return false;

    if (eta >= EtaCut[0] && eta < EtaCut[1]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclTPCDedxMomentumNegAfter[0]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclTPCDedxMomentumPosAfter[0]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[1] && eta < EtaCut[2]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclTPCDedxMomentumNegAfter[1]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclTPCDedxMomentumPosAfter[1]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[2] && eta < EtaCut[3]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclTPCDedxMomentumNegAfter[2]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclTPCDedxMomentumPosAfter[2]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[3] && eta < EtaCut[4]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclTPCDedxMomentumNegAfter[3]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclTPCDedxMomentumPosAfter[3]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[4] && eta < EtaCut[5]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclTPCDedxMomentumNegAfter[4]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclTPCDedxMomentumPosAfter[4]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[5] && eta < EtaCut[6]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclTPCDedxMomentumNegAfter[5]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclTPCDedxMomentumPosAfter[5]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[6] && eta < EtaCut[7]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclTPCDedxMomentumNegAfter[6]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclTPCDedxMomentumPosAfter[6]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[7] && eta < EtaCut[8]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclTPCDedxMomentumNegAfter[7]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclTPCDedxMomentumPosAfter[7]), nTPCCl, trk.tpcSignal(), sigP);
      }
    }

    return true;
  }

  // NclPIDCutTPC
  template <typename T>
  bool passedNClTPCPIDCutPri(const T& trk)
  {
    float eta = trk.eta();
    float sigP = trk.sign() * getMomentum(trk);
    auto nTPCCl = trk.tpcNClsPID();

    if (eta >= EtaCut[0] && eta < EtaCut[1]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumNegBefore[0]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumPosBefore[0]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[1] && eta < EtaCut[2]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumNegBefore[1]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumPosBefore[1]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[2] && eta < EtaCut[3]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumNegBefore[2]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumPosBefore[2]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[3] && eta < EtaCut[4]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumNegBefore[3]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumPosBefore[3]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[4] && eta < EtaCut[5]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumNegBefore[4]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumPosBefore[4]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[5] && eta < EtaCut[6]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumNegBefore[5]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumPosBefore[5]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[6] && eta < EtaCut[7]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumNegBefore[6]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumPosBefore[6]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[7] && eta < EtaCut[8]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumNegBefore[7]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumPosBefore[7]), nTPCCl, trk.tpcSignal(), sigP);
      }
    }

    if (nTPCCl < minTPCnClsPID)
      return false;

    if (eta >= EtaCut[0] && eta < EtaCut[1]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumNegAfter[0]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumPosAfter[0]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[1] && eta < EtaCut[2]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumNegAfter[1]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumPosAfter[1]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[2] && eta < EtaCut[3]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumNegAfter[2]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumPosAfter[2]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[3] && eta < EtaCut[4]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumNegAfter[3]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumPosAfter[3]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[4] && eta < EtaCut[5]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumNegAfter[4]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumPosAfter[4]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[5] && eta < EtaCut[6]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumNegAfter[5]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumPosAfter[5]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[6] && eta < EtaCut[7]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumNegAfter[6]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumPosAfter[6]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta >= EtaCut[7] && eta < EtaCut[8]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumNegAfter[7]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(NclPIDTPCDedxMomentumPosAfter[7]), nTPCCl, trk.tpcSignal(), sigP);
      }
    }

    return true;
  }

  // Phi cut Secondaries
  template <typename T>
  bool passedPhiCutSecondaries(const T& trk, float magField, const TF1& fphiCutLow, const TF1& fphiCutHigh)
  {
    float pt = trk.pt();
    float phi = trk.phi();
    int charge = trk.sign();

    if (pt < pTcut)
      return true;

    if (magField < 0) // for negatve polarity field
      phi = o2::constants::math::TwoPI - phi;
    if (charge < 0) // for negatve charge
      phi = o2::constants::math::TwoPI - phi;

    // to center gap in the middle
    phi += o2::constants::math::PI / 18.0f;
    phi = std::fmod(phi, o2::constants::math::PI / 9.0f);

    // cut phi
    if (phi < fphiCutHigh.Eval(pt) && phi > fphiCutLow.Eval(pt))
      return false; // reject track

    return true;
  }

  // NclCutTPC
  template <typename T>
  bool passedNClTPCFoundCutSecondaries(const T& trk)
  {
    auto nTPCCl = trk.tpcNClsFound();

    if (nTPCCl < minTPCnClsFound)
      return false;

    return true;
  }

  // NclCutPIDTPC secondary
  template <typename T>
  bool passedNClTPCPIDCutSecondaries(const T& trk)
  {
    auto nTPCCl = trk.tpcNClsPID();

    if (nTPCCl < minTPCnClsPID)
      return false;

    return true;
  }

  // Process Data
  void process(SelectedCollisions::iterator const& collision, BCsRun3 const& /**/,
               aod::V0Datas const& fullV0s, PIDTracks const& tracks)
  {
    registryDeDx.fill(HIST("evsel"), EvCutLabel::AllEv);
    // Event Selection
    if (!collision.sel8())
      return;

    registryDeDx.fill(HIST("evsel"), EvCutLabel::SelEigth);

    if (std::abs(collision.posZ()) > maxZDistanceToIP)
      return;

    registryDeDx.fill(HIST("evsel"), EvCutLabel::ZVtxCut);

    if (nPileUp) {
      if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
        return;
      registryDeDx.fill(HIST("evsel"), EvCutLabel::NoSameBunchPileup);
    }

    if (nGoodZvtx) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
        return;
      registryDeDx.fill(HIST("evsel"), EvCutLabel::GoodZvtxFT0vsPV);
    }

    if (nINELSelectionMode == NoSelINEL) {
    } else if (nINELSelectionMode == SelINELgt0) {
      if (!collision.isInelGt0())
        return;
      registryDeDx.fill(HIST("evsel"), EvCutLabel::INELgt);
    } else if (nINELSelectionMode == SelINELgt1) {
      if (!collision.isInelGt1())
        return;
      registryDeDx.fill(HIST("evsel"), EvCutLabel::INELgt);
    }

    // Event Counter
    registryDeDx.fill(HIST("histRecVtxZData"), collision.posZ());

    // For magnetic field
    const auto& foundBC = collision.foundBC_as<BCsRun3>();
    const uint64_t timeStamp{foundBC.timestamp()};
    const int magField{getMagneticField(timeStamp)};

    float centrality = collision.centFT0C();
    if (centrality < CentClasses[0] || centrality > CentClasses[10])
      return;

    // Event Counter by cent
    registryDeDx.fill(HIST("histCentrality"), centrality);

    for (const auto& trk : tracks) {
      registryDeDx.fill(HIST("trackselAll"), TrkPriCutLabel::AllPri);

      // track Selection
      if (!mySelectionPrim.IsSelected(trk))
        continue;

      registryDeDx.fill(HIST("trackselAll"), TrkPriCutLabel::SelectionPrim);

      // phi and Ncl cut
      if (phiVarCut) {
        if (!passedPhiCutPri(trk, magField, *fphiCutLow, *fphiCutHigh))
          continue;
        registryDeDx.fill(HIST("trackselAll"), TrkPriCutLabel::PhiVarCutPri);
      }

      // NCl cut PID TPC
      registryDeDx.fill(HIST("heta_vs_NclPID_TPC_Before_Primary"), trk.eta(), trk.tpcNClsPID());
      if (nClTPCPIDCut) {
        if (!passedNClTPCPIDCutPri(trk))
          continue;
        registryDeDx.fill(HIST("trackselAll"), TrkPriCutLabel::NClTPCPIDCutPri);
      }

      registryDeDx.fill(HIST("heta_vs_NclPID_TPC_After_Primary"), trk.eta(), trk.tpcNClsPID());

      // NCl cut TPC
      registryDeDx.fill(HIST("heta_vs_NclFound_TPC_Before_Primary"), trk.eta(), trk.tpcNClsFound());
      if (nClTPCFoundCut) {
        if (!passedNClTPCFoundCutPri(trk))
          continue;
        registryDeDx.fill(HIST("trackselAll"), TrkPriCutLabel::NClTPCFoundCutPri);
      }

      registryDeDx.fill(HIST("heta_vs_NclFound_TPC_After_Primary"), trk.eta(), trk.tpcNClsFound());

      float signedP = trk.sign() * getMomentum(trk);

      // MIP calibration for pions
      if (getMomentum(trk) >= pionMin && getMomentum(trk) <= pionMax) {
        if (calibrationMode) {
          if (signedP < 0) {
            registryDeDx.fill(HIST("hdEdx_vs_eta_Neg_Pi"), trk.eta(), trk.tpcSignal());
          } else {
            registryDeDx.fill(HIST("hdEdx_vs_eta_Pos_Pi"), trk.eta(), trk.tpcSignal());
          }

        } else {
          for (int i = 0; i < EtaIntervals; ++i) {
            if (trk.eta() >= EtaCut[i] && trk.eta() < EtaCut[i + 1]) {
              if (signedP < 0) {
                registryDeDx.fill(HIST("hdEdx_vs_eta_Neg_calibrated_Pi"), trk.eta(), trk.tpcSignal() * 50 / calibrationFactorNeg->at(i));
              } else {
                registryDeDx.fill(HIST("hdEdx_vs_eta_Pos_calibrated_Pi"), trk.eta(), trk.tpcSignal() * 50 / calibrationFactorPos->at(i));
              }
            }
          }
        }
      }
      // Beta from TOF
      if (signedP < 0) {
        registryDeDx.fill(HIST("hbeta_vs_p_Neg"), std::abs(signedP), trk.beta());
      } else {
        registryDeDx.fill(HIST("hbeta_vs_p_Pos"), signedP, trk.beta());
      }
      // Electrons from TOF
      if (std::abs(trk.beta() - 1) < elTofCut) { // beta cut
        if (calibrationMode) {
          if (signedP < 0) {
            registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Neg_El"), trk.eta(), trk.tpcSignal(), std::abs(signedP));
          } else {
            registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Pos_El"), trk.eta(), trk.tpcSignal(), signedP);
          }
        } else {
          for (int i = 0; i < EtaIntervals; ++i) {
            if (trk.eta() >= EtaCut[i] && trk.eta() < EtaCut[i + 1]) {
              if (signedP < 0) {
                registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Neg_calibrated_El"), trk.eta(), trk.tpcSignal() * 50 / calibrationFactorNeg->at(i), std::abs(signedP));
              } else {
                registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Pos_calibrated_El"), trk.eta(), trk.tpcSignal() * 50 / calibrationFactorPos->at(i), signedP);
              }
            }
          }
        }
      }
      // pions from TOF
      if (trk.beta() > pionTofCut && trk.beta() < pionTofCut + 0.05) { // beta cut
        if (calibrationMode) {
          if (signedP < 0) {
            registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Neg_TOF"), trk.eta(), trk.tpcSignal(), std::abs(signedP));
          } else {
            registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Pos_TOF"), trk.eta(), trk.tpcSignal(), signedP);
          }
        } else {
          for (int i = 0; i < EtaIntervals; ++i) {
            if (trk.eta() >= EtaCut[i] && trk.eta() < EtaCut[i + 1]) {
              if (signedP < 0) {
                registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Neg_calibrated_TOF"), trk.eta(), trk.tpcSignal() * 50 / calibrationFactorNeg->at(i), std::abs(signedP));
              } else {
                registryDeDx.fill(HIST("hdEdx_vs_eta_vs_p_Pos_calibrated_TOF"), trk.eta(), trk.tpcSignal() * 50 / calibrationFactorPos->at(i), signedP);
              }
            }
          }
        }
      }

      registryDeDx.fill(HIST("hdEdx_vs_phi"), trk.phi(), trk.tpcSignal());

      if (!calibrationMode) {
        int centIndex = -1;
        for (int j = 0; j < CentralityClasses; ++j) {
          if (centrality >= CentClasses[j] && centrality < CentClasses[j + 1]) {
            centIndex = j;
            break;
          }
        }
        if (centIndex == -1)
          continue;

        for (int i = 0; i < EtaIntervals; ++i) {
          if (trk.eta() >= EtaCut[i] && trk.eta() < EtaCut[i + 1]) {
            if (signedP > 0) {
              registryDeDx.fill(HIST(DedxvsMomentumPos[0]), signedP, trk.tpcSignal() * 50 / calibrationFactorPos->at(i), trk.eta());
              registryDeDx.fill(HIST("heta_vs_pt_vs_p_all_Pos"), trk.eta(), trk.pt(), trk.p());
              hDedxVsMomentumVsCentPos[centIndex]->Fill(signedP, trk.tpcSignal() * 50 / calibrationFactorPos->at(i), trk.eta());
            } else {
              registryDeDx.fill(HIST(DedxvsMomentumNeg[0]), std::abs(signedP), trk.tpcSignal() * 50 / calibrationFactorNeg->at(i), trk.eta());
              registryDeDx.fill(HIST("heta_vs_pt_vs_p_all_Neg"), trk.eta(), trk.pt(), trk.p());
              hDedxVsMomentumVsCentNeg[centIndex]->Fill(std::abs(signedP), trk.tpcSignal() * 50 / calibrationFactorNeg->at(i), trk.eta());
            }
          }
        }
      }
    }

    // Loop over Reconstructed V0s
    if (!calibrationMode) {
      for (const auto& v0 : fullV0s) {

        // Standard V0 Selections
        registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::AllSec);

        // Select V0 type
        if (v0.v0Type() != v0TypeSelection)
          continue;

        registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::V0Type);

        if (!passedV0Selection(v0, collision)) {
          continue;
        }

        // Positive and Negative Tracks
        const auto& posTrack = v0.posTrack_as<PIDTracks>();
        const auto& negTrack = v0.negTrack_as<PIDTracks>();

        if (!posTrack.passedTPCRefit())
          continue;
        if (!negTrack.passedTPCRefit())
          continue;

        registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::TPCRefit);
        // phi and Ncl cut
        if (phiVarCut) {
          if (!passedPhiCutSecondaries(posTrack, magField, *fphiCutLow, *fphiCutHigh))
            continue;

          if (!passedPhiCutSecondaries(negTrack, magField, *fphiCutLow, *fphiCutHigh))
            continue;

          registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel ::PhiVarCutSec);
        }

        fillHist = false;
        // K0s Selection
        if (passedK0Selection(v0, negTrack, posTrack, collision)) {
          registryDeDx.fill(HIST("heta_vs_NclPID_TPC_Before_PionsK0s"), posTrack.eta(), posTrack.tpcNClsPID());
          registryDeDx.fill(HIST("heta_vs_NclPID_TPC_Before_PionsK0s"), negTrack.eta(), negTrack.tpcNClsPID());
          registryDeDx.fill(HIST("hp_vs_NclPID_TPC_Before_PionsK0s"), posTrack.p(), posTrack.tpcNClsPID());
          registryDeDx.fill(HIST("hp_vs_NclPID_TPC_Before_PionsK0s"), negTrack.p(), negTrack.tpcNClsPID());
        }
        // Lambda Selection
        if (passedLambdaSelection(v0, negTrack, posTrack, collision)) {
          registryDeDx.fill(HIST("heta_vs_NclPID_TPC_Before_ProtonsLambda"), posTrack.eta(), posTrack.tpcNClsPID());
          registryDeDx.fill(HIST("heta_vs_NclPID_TPC_Before_PionsLambda"), negTrack.eta(), negTrack.tpcNClsPID());
          registryDeDx.fill(HIST("hp_vs_NclPID_TPC_Before_ProtonsLambda"), posTrack.p(), posTrack.tpcNClsPID());
          registryDeDx.fill(HIST("hp_vs_NclPID_TPC_Before_PionsLambda"), negTrack.p(), negTrack.tpcNClsPID());
        }
        // AntiLambda Selection
        if (passedAntiLambdaSelection(v0, negTrack, posTrack, collision)) {
          registryDeDx.fill(HIST("heta_vs_NclPID_TPC_Before_PionsLambda"), posTrack.eta(), posTrack.tpcNClsPID());
          registryDeDx.fill(HIST("heta_vs_NclPID_TPC_Before_ProtonsLambda"), negTrack.eta(), negTrack.tpcNClsPID());
          registryDeDx.fill(HIST("hp_vs_NclPID_TPC_Before_PionsLambda"), posTrack.p(), posTrack.tpcNClsPID());
          registryDeDx.fill(HIST("hp_vs_NclPID_TPC_Before_ProtonsLambda"), negTrack.p(), negTrack.tpcNClsPID());
        }
        // Gamma Selection
        if (passedGammaSelection(v0, negTrack, posTrack, collision)) {
          registryDeDx.fill(HIST("heta_vs_NclPID_TPC_Before_ElectronsGamma"), posTrack.eta(), posTrack.tpcNClsPID());
          registryDeDx.fill(HIST("heta_vs_NclPID_TPC_Before_ElectronsGamma"), negTrack.eta(), negTrack.tpcNClsPID());
          registryDeDx.fill(HIST("hp_vs_NclPID_TPC_Before_ElectronsGamma"), posTrack.p(), posTrack.tpcNClsPID());
          registryDeDx.fill(HIST("hp_vs_NclPID_TPC_Before_ElectronsGamma"), negTrack.p(), negTrack.tpcNClsPID());
        }

        if (nClTPCFoundCut) {
          if (!passedNClTPCFoundCutSecondaries(posTrack))
            continue;

          if (!passedNClTPCFoundCutSecondaries(negTrack))
            continue;

          registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel ::NClTPCFoundCutSec);
        }

        if (nClTPCPIDCut) {
          if (!passedNClTPCPIDCutSecondaries(posTrack))
            continue;

          if (!passedNClTPCPIDCutSecondaries(negTrack))
            continue;

          registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel ::NClTPCPIDCutSec);
        }

        float signedPpos = posTrack.sign() * getMomentum(posTrack);
        float signedPneg = negTrack.sign() * getMomentum(negTrack);

        fillHist = true;

        // K0s Selection
        if (passedK0Selection(v0, negTrack, posTrack, collision)) {
          registryDeDx.fill(HIST("heta_vs_NclPID_TPC_After_PionsK0s"), posTrack.eta(), posTrack.tpcNClsPID());
          registryDeDx.fill(HIST("heta_vs_NclPID_TPC_After_PionsK0s"), negTrack.eta(), negTrack.tpcNClsPID());
          registryDeDx.fill(HIST("hp_vs_NclPID_TPC_After_PionsK0s"), posTrack.p(), posTrack.tpcNClsPID());
          registryDeDx.fill(HIST("hp_vs_NclPID_TPC_After_PionsK0s"), negTrack.p(), negTrack.tpcNClsPID());
          for (int i = 0; i < EtaIntervals; ++i) {
            if (negTrack.eta() > EtaCut[i] && negTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(DedxvsMomentumNeg[1]), std::abs(signedPneg), negTrack.tpcSignal() * 50 / calibrationFactorNeg->at(i), negTrack.eta());
            }
            if (posTrack.eta() > EtaCut[i] && posTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(DedxvsMomentumPos[1]), signedPpos, posTrack.tpcSignal() * 50 / calibrationFactorPos->at(i), posTrack.eta());
            }
          }
        }

        // Lambda Selection
        if (passedLambdaSelection(v0, negTrack, posTrack, collision)) {
          registryDeDx.fill(HIST("heta_vs_NclPID_TPC_After_ProtonsLambda"), posTrack.eta(), posTrack.tpcNClsPID());
          registryDeDx.fill(HIST("heta_vs_NclPID_TPC_After_PionsLambda"), negTrack.eta(), negTrack.tpcNClsPID());
          registryDeDx.fill(HIST("hp_vs_NclPID_TPC_After_ProtonsLambda"), posTrack.p(), posTrack.tpcNClsPID());
          registryDeDx.fill(HIST("hp_vs_NclPID_TPC_After_PionsLambda"), negTrack.p(), negTrack.tpcNClsPID());
          for (int i = 0; i < EtaIntervals; ++i) {
            if (negTrack.eta() > EtaCut[i] && negTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(DedxvsMomentumNeg[1]), std::abs(signedPneg), negTrack.tpcSignal() * 50 / calibrationFactorNeg->at(i), negTrack.eta());
            }
            if (posTrack.eta() > EtaCut[i] && posTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(DedxvsMomentumPos[2]), signedPpos, posTrack.tpcSignal() * 50 / calibrationFactorPos->at(i), posTrack.eta());
            }
          }
        }

        // AntiLambda Selection
        if (passedAntiLambdaSelection(v0, negTrack, posTrack, collision)) {
          registryDeDx.fill(HIST("heta_vs_NclPID_TPC_After_PionsLambda"), posTrack.eta(), posTrack.tpcNClsPID());
          registryDeDx.fill(HIST("heta_vs_NclPID_TPC_After_ProtonsLambda"), negTrack.eta(), negTrack.tpcNClsPID());
          registryDeDx.fill(HIST("hp_vs_NclPID_TPC_After_PionsLambda"), posTrack.p(), posTrack.tpcNClsPID());
          registryDeDx.fill(HIST("hp_vs_NclPID_TPC_After_ProtonsLambda"), negTrack.p(), negTrack.tpcNClsPID());
          for (int i = 0; i < EtaIntervals; ++i) {
            if (negTrack.eta() > EtaCut[i] && negTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(DedxvsMomentumNeg[2]), std::abs(signedPneg), negTrack.tpcSignal() * 50 / calibrationFactorNeg->at(i), negTrack.eta());
            }
            if (posTrack.eta() > EtaCut[i] && posTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(DedxvsMomentumPos[1]), signedPpos, posTrack.tpcSignal() * 50 / calibrationFactorPos->at(i), posTrack.eta());
            }
          }
        }

        // Gamma Selection
        if (passedGammaSelection(v0, negTrack, posTrack, collision)) {
          registryDeDx.fill(HIST("heta_vs_NclPID_TPC_After_ElectronsGamma"), posTrack.eta(), posTrack.tpcNClsPID());
          registryDeDx.fill(HIST("heta_vs_NclPID_TPC_After_ElectronsGamma"), negTrack.eta(), negTrack.tpcNClsPID());
          registryDeDx.fill(HIST("hp_vs_NclPID_TPC_After_ElectronsGamma"), posTrack.p(), posTrack.tpcNClsPID());
          registryDeDx.fill(HIST("hp_vs_NclPID_TPC_After_ElectronsGamma"), negTrack.p(), negTrack.tpcNClsPID());
          for (int i = 0; i < EtaIntervals; ++i) {
            if (negTrack.eta() > EtaCut[i] && negTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(DedxvsMomentumNeg[3]), std::abs(signedPneg), negTrack.tpcSignal() * 50 / calibrationFactorNeg->at(i), negTrack.eta());
            }
            if (posTrack.eta() > EtaCut[i] && posTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(DedxvsMomentumPos[3]), signedPpos, posTrack.tpcSignal() * 50 / calibrationFactorPos->at(i), posTrack.eta());
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<DedxPidAnalysis>(cfgc)};
}
