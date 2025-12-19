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
/// \file dedxAnalysis.cxx
/// \brief  Analysis to do PID

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Logger.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/V0.h"

#include "TF1.h"

using namespace o2;
using namespace o2::framework;
using namespace constants::physics;

using PIDTracks = soa::Join<
  aod::Tracks, aod::TracksExtra, aod::TrackSelectionExtension, aod::TracksDCA, aod::TrackSelection,
  aod::pidTOFFullPi, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFbeta>;

using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using BCsRun3 = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

struct DedxAnalysis {

  // dE/dx for all charged particles
  HistogramRegistry registryDeDx{
    "registryDeDx",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  // Constant values
  static constexpr int kEtaIntervals = 8;
  static constexpr int kParticlesType = 4;
  float tpcCut = 0.6;
  float pionMin = 0.35;
  float pionMax = 0.45;
  float elTofCut = 0.1;
  float pionTofCut = 1.0;
  float invMassCut = 0.01;
  float invMassCutGamma = 0.0015;
  float pTcut = 2.0;

  // Event cut labels
  enum EvCutLabel {
    AllEv = 1,
    SelEigth,
    ZVtxCut
  };

  // Track primary label
  enum TrkPriCutLabel {
    AllPri = 1,
    SelectionPrim,
    PhiVarCutPri,
    NTPCClCutPri,
    NITSClCutPri
  };

  // Track secondary lebel
  enum TrkSecCutLabel {
    AllSec = 1,
    V0CosPA,
    V0DecayRadius,
    V0Daughters,
    TPCRefit,
    PhiVarCutSec,
    NTPCClCutSec,
    NITSClCutSec,
    V0RapidityK0s,
    V0ProperLifetimeK0s,
    V0RapidityLambda,
    V0ProperLifetimeLambda,
    V0RapidityAntiLambda,
    V0ProperLifetimeAntiLambda
  };
  // Configurable Parameters
  // Tracks cuts
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 70.0f,
                                      "min number of found TPC clusters"};
  Configurable<float> minITSnCls{"minITSnCls", 70.0f,
                                 "min number of ITS clusters"};
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
  Configurable<float> minMassK0s{"minMassK0s", 0.4f, "Minimum Mass K0s"};
  Configurable<float> maxMassK0s{"maxMassK0s", 0.6f, "Maximum Mass K0s"};
  Configurable<float> minMassLambda{"minMassLambda", 1.1f,
                                    "Minimum Mass Lambda"};
  Configurable<float> maxMassLambda{"maxMassLambda", 1.2f,
                                    "Maximum Mass Lambda"};
  Configurable<float> minMassGamma{"minMassGamma", 0.000922f,
                                   "Minimum Mass Gamma"};
  Configurable<float> maxMassGamma{"maxMassGamma", 0.002022f,
                                   "Maximum Mass Gamma"};
  Configurable<bool> calibrationMode{"calibrationMode", false, "calibration mode"};
  Configurable<bool> phiVarCut{"phiVarCut", true, "phi var cut"};
  Configurable<bool> nTPCClCut{"nTPCClCut", true, "number of clusters in TPC cut"};
  Configurable<bool> nITSClCut{"nITSClCut", true, "number of clusters in ITS cut"};
  // Histograms names
  static constexpr std::string_view kDedxvsMomentumPos[kParticlesType] = {"dEdx_vs_Momentum_all_Pos", "dEdx_vs_Momentum_Pi_v0_Pos", "dEdx_vs_Momentum_Pr_v0_Pos", "dEdx_vs_Momentum_El_v0_Pos"};
  static constexpr std::string_view kDedxvsMomentumNeg[kParticlesType] = {"dEdx_vs_Momentum_all_Neg", "dEdx_vs_Momentum_Pi_v0_Neg", "dEdx_vs_Momentum_Pr_v0_Neg", "dEdx_vs_Momentum_El_v0_Neg"};
  // Ncl TPC
  static constexpr std::string_view kNclTPCDedxMomentumNegBefore[kEtaIntervals] = {"Ncl_TPC_vs_dEdx_vs_Momentum_Neg_1_Before", "Ncl_TPC_vs_dEdx_vs_Momentum_Neg_2_Before", "Ncl_TPC_vs_dEdx_vs_Momentum_Neg_3_Before", "Ncl_TPC_vs_dEdx_vs_Momentum_Neg_4_Before", "Ncl_TPC_vs_dEdx_vs_Momentum_Neg_5_Before", "Ncl_TPC_vs_dEdx_vs_Momentum_Neg_6_Before", "Ncl_TPC_vs_dEdx_vs_Momentum_Neg_7_Before", "Ncl_TPC_vs_dEdx_vs_Momentum_Neg_8_Before"};
  static constexpr std::string_view kNclTPCDedxMomentumPosBefore[kEtaIntervals] = {"Ncl_TPC_vs_dEdx_vs_Momentum_Pos_1_Before", "Ncl_TPC_vs_dEdx_vs_Momentum_Pos_2_Before", "Ncl_TPC_vs_dEdx_vs_Momentum_Pos_3_Before", "Ncl_TPC_vs_dEdx_vs_Momentum_Pos_4_Before", "Ncl_TPC_vs_dEdx_vs_Momentum_Pos_5_Before", "Ncl_TPC_vs_dEdx_vs_Momentum_Pos_6_Before", "Ncl_TPC_vs_dEdx_vs_Momentum_Pos_7_Before", "Ncl_TPC_vs_dEdx_vs_Momentum_Pos_8_Before"};
  static constexpr std::string_view kNclTPCDedxMomentumNegAfter[kEtaIntervals] = {"Ncl_TPC_vs_dEdx_vs_Momentum_Neg_1_After", "Ncl_TPC_vs_dEdx_vs_Momentum_Neg_2_After", "Ncl_TPC_vs_dEdx_vs_Momentum_Neg_3_After", "Ncl_TPC_vs_dEdx_vs_Momentum_Neg_4_After", "Ncl_TPC_vs_dEdx_vs_Momentum_Neg_5_After", "Ncl_TPC_vs_dEdx_vs_Momentum_Neg_6_After", "Ncl_TPC_vs_dEdx_vs_Momentum_Neg_7_After", "Ncl_TPC_vs_dEdx_vs_Momentum_Neg_8_After"};
  static constexpr std::string_view kNclTPCDedxMomentumPosAfter[kEtaIntervals] = {"Ncl_TPC_vs_dEdx_vs_Momentum_Pos_1_After", "Ncl_TPC_vs_dEdx_vs_Momentum_Pos_2_After", "Ncl_TPC_vs_dEdx_vs_Momentum_Pos_3_After", "Ncl_TPC_vs_dEdx_vs_Momentum_Pos_4_After", "Ncl_TPC_vs_dEdx_vs_Momentum_Pos_5_After", "Ncl_TPC_vs_dEdx_vs_Momentum_Pos_6_After", "Ncl_TPC_vs_dEdx_vs_Momentum_Pos_7_After", "Ncl_TPC_vs_dEdx_vs_Momentum_Pos_8_After"};
  // Ncl TPC
  static constexpr std::string_view kNclITSDedxMomentumNegBefore[kEtaIntervals] = {"Ncl_ITS_vs_dEdx_vs_Momentum_Neg_1_Before", "Ncl_ITS_vs_dEdx_vs_Momentum_Neg_2_Before", "Ncl_ITS_vs_dEdx_vs_Momentum_Neg_3_Before", "Ncl_ITS_vs_dEdx_vs_Momentum_Neg_4_Before", "Ncl_ITS_vs_dEdx_vs_Momentum_Neg_5_Before", "Ncl_ITS_vs_dEdx_vs_Momentum_Neg_6_Before", "Ncl_ITS_vs_dEdx_vs_Momentum_Neg_7_Before", "Ncl_ITS_vs_dEdx_vs_Momentum_Neg_8_Before"};
  static constexpr std::string_view kNclITSDedxMomentumPosBefore[kEtaIntervals] = {"Ncl_ITS_vs_dEdx_vs_Momentum_Pos_1_Before", "Ncl_ITS_vs_dEdx_vs_Momentum_Pos_2_Before", "Ncl_ITS_vs_dEdx_vs_Momentum_Pos_3_Before", "Ncl_ITS_vs_dEdx_vs_Momentum_Pos_4_Before", "Ncl_ITS_vs_dEdx_vs_Momentum_Pos_5_Before", "Ncl_ITS_vs_dEdx_vs_Momentum_Pos_6_Before", "Ncl_ITS_vs_dEdx_vs_Momentum_Pos_7_Before", "Ncl_ITS_vs_dEdx_vs_Momentum_Pos_8_Before"};
  static constexpr std::string_view kNclITSDedxMomentumNegAfter[kEtaIntervals] = {"Ncl_ITS_vs_dEdx_vs_Momentum_Neg_1_After", "Ncl_ITS_vs_dEdx_vs_Momentum_Neg_2_After", "Ncl_ITS_vs_dEdx_vs_Momentum_Neg_3_After", "Ncl_ITS_vs_dEdx_vs_Momentum_Neg_4_After", "Ncl_ITS_vs_dEdx_vs_Momentum_Neg_5_After", "Ncl_ITS_vs_dEdx_vs_Momentum_Neg_6_After", "Ncl_ITS_vs_dEdx_vs_Momentum_Neg_7_After", "Ncl_ITS_vs_dEdx_vs_Momentum_Neg_8_After"};
  static constexpr std::string_view kNclITSDedxMomentumPosAfter[kEtaIntervals] = {"Ncl_ITS_vs_dEdx_vs_Momentum_Pos_1_After", "Ncl_ITS_vs_dEdx_vs_Momentum_Pos_2_After", "Ncl_ITS_vs_dEdx_vs_Momentum_Pos_3_After", "Ncl_ITS_vs_dEdx_vs_Momentum_Pos_4_After", "Ncl_ITS_vs_dEdx_vs_Momentum_Pos_5_After", "Ncl_ITS_vs_dEdx_vs_Momentum_Pos_6_After", "Ncl_ITS_vs_dEdx_vs_Momentum_Pos_7_After", "Ncl_ITS_vs_dEdx_vs_Momentum_Pos_8_After"};
  static constexpr double EtaCut[kEtaIntervals + 1] = {-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8};
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
    AxisSpec dedxAxis{100, 0.0, 100.0, "dE/dx (a. u.)"};
    AxisSpec ptAxis = {binP, "pT (GeV/c)"};
    AxisSpec etaAxis{8, -0.8, 0.8, "#eta"};
    AxisSpec pAxis = {binP, "#it{p}/Z (GeV/c)"};
    fphiCutLow = new TF1("StandardPhiCutLow", "0.119297/x/x+pi/18.0-0.000379693", 0, 50);
    fphiCutHigh = new TF1("StandardPhiCutHigh", "0.16685/x+pi/18.0+0.00981942", 0, 50);
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
        "heta_vs_p_vs_pt_all_Neg", "eta_vs_p_vs_pT", HistType::kTH3F,
        {{etaAxis}, {ptAxis}, {pAxis}});
      registryDeDx.add(
        "heta_vs_p_vs_pt_all_Pos", "eta_vs_p_vs_pT", HistType::kTH3F,
        {{etaAxis}, {ptAxis}, {pAxis}});

      // De/Dx for ch and v0 particles
      for (int i = 0; i < kParticlesType; ++i) {
        registryDeDx.add(kDedxvsMomentumPos[i].data(), "dE/dx", HistType::kTH3F,
                         {{pAxis}, {dedxAxis}, {etaAxis}});
        registryDeDx.add(kDedxvsMomentumNeg[i].data(), "dE/dx", HistType::kTH3F,
                         {{pAxis}, {dedxAxis}, {etaAxis}});
      }
    }

    registryDeDx.add(
      "hdEdx_vs_phi", "dE/dx", HistType::kTH2F,
      {{100, 0.0, 6.4, "#phi"}, {dedxAxis}});
    // phi cut
    if (phiVarCut) {
      registryDeDx.add(
        "hpt_vs_phi_Ncl_TPC_After", "phi cut", HistType::kTH3F,
        {{ptAxis}, {100, 0.0, 0.4, "#varphi^{'}"}, {100, 0, 160, "N_{cl}"}});

      registryDeDx.add(
        "hpt_vs_phi_Ncl_TPC_Before", "phi cut", HistType::kTH3F,
        {{ptAxis}, {100, 0.0, 0.4, "#varphi^{'}"}, {100, 0, 160, "N_{cl}"}});

      // Ncl vs de/dx TPC

      for (int i = 0; i < kEtaIntervals; ++i) {
        registryDeDx.add(kNclTPCDedxMomentumPosBefore[i].data(), "Ncl TPC vs dE/dx vs Momentum Positive before", HistType::kTH3F,
                         {{100, 0, 160, "N_{cl}^{TPC}"}, {dedxAxis}, {pAxis}});
        registryDeDx.add(kNclTPCDedxMomentumNegBefore[i].data(), "Ncl TPC vs dE/dx vs Momentum Negative before", HistType::kTH3F,
                         {{100, 0, 160, "N_{cl}^{TPC}"}, {dedxAxis}, {pAxis}});
        if (nTPCClCut) {
          registryDeDx.add(kNclTPCDedxMomentumPosAfter[i].data(), "Ncl TPC vs dE/dx vs Momentum Positive after", HistType::kTH3F,
                           {{100, 0, 160, "N_{cl}^{TPC}"}, {dedxAxis}, {pAxis}});
          registryDeDx.add(kNclTPCDedxMomentumNegAfter[i].data(), "Ncl TPC vs dE/dx vs Momentum Negative after", HistType::kTH3F,
                           {{100, 0, 160, "N_{cl}^{TPC}"}, {dedxAxis}, {pAxis}});
        }
      }
    }

    // Ncl vs de/dx ITS
    if (nITSClCut) {
      for (int i = 0; i < kEtaIntervals; ++i) {
        registryDeDx.add(kNclITSDedxMomentumPosBefore[i].data(), "Ncl ITS vs dE/dx vs Momentum Positive before", HistType::kTH3F,
                         {{100, 0, 160, "N_{cl}^{ITS}"}, {dedxAxis}, {pAxis}});
        registryDeDx.add(kNclITSDedxMomentumNegBefore[i].data(), "Ncl ITS vs dE/dx vs Momentum Negative before", HistType::kTH3F,
                         {{100, 0, 160, "N_{cl}^{ITS}"}, {dedxAxis}, {pAxis}});

        registryDeDx.add(kNclITSDedxMomentumPosAfter[i].data(), "Ncl ITS vs dE/dx vs Momentum Positive after", HistType::kTH3F,
                         {{100, 0, 160, "N_{cl}^{ITS}"}, {dedxAxis}, {pAxis}});
        registryDeDx.add(kNclITSDedxMomentumNegAfter[i].data(), "Ncl ITS vs dE/dx vs Momentum Negative after", HistType::kTH3F,
                         {{100, 0, 160, "N_{cl}^{ITS}"}, {dedxAxis}, {pAxis}});
      }
    }

    // beta plot
    registryDeDx.add(
      "hbeta_vs_p_Neg", "beta", HistType::kTH2F,
      {{pAxis}, {100, 0.0, 1.1, "#beta"}});

    registryDeDx.add(
      "hbeta_vs_p_Pos", "beta", HistType::kTH2F,
      {{pAxis}, {100, 0.0, 1.1, "#beta"}});

    // Event Counter
    registryDeDx.add("histRecVtxZData", "collision z position", HistType::kTH1F, {{100, -20.0, +20.0, "z_{vtx} (cm)"}});

    // Event Counter
    registryDeDx.add("evsel", "events selected", HistType::kTH1F, {{3, 0.5, 3.5, ""}});
    auto hstat = registryDeDx.get<TH1>(HIST("evsel"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "AllEv");
    x->SetBinLabel(2, "SelEigth");
    x->SetBinLabel(3, "ZVtxCut");

    // Track Counter
    registryDeDx.add("trackselAll", "track selected all particles", HistType::kTH1F, {{5, 0.5, 5.5, ""}});
    auto htrackAll = registryDeDx.get<TH1>(HIST("trackselAll"));
    auto* xAll = htrackAll->GetXaxis();
    xAll->SetBinLabel(1, "AllPri");
    xAll->SetBinLabel(2, "SelectionPrim");
    xAll->SetBinLabel(3, "PhiVarCutPri");
    xAll->SetBinLabel(4, "NTPCClCutPri");
    xAll->SetBinLabel(5, "NITSClCutPri");

    registryDeDx.add("trackselSec", "track selected sec particles", HistType::kTH1F, {{13, 0.5, 13.5, ""}});
    auto htrackSec = registryDeDx.get<TH1>(HIST("trackselSec"));
    auto* xSec = htrackSec->GetXaxis();
    xSec->SetBinLabel(1, "AllSec");
    xSec->SetBinLabel(2, "V0CosPA");
    xSec->SetBinLabel(3, "V0DecayRadius");
    xSec->SetBinLabel(4, "V0Daughters");
    xSec->SetBinLabel(5, "TPCRefit");
    xSec->SetBinLabel(6, "PhiVarCutSec");
    xSec->SetBinLabel(7, "NTPCClCutSec");
    xSec->SetBinLabel(8, "NITSClCutSec");
    xSec->SetBinLabel(9, "V0RapidityK0s");
    xSec->SetBinLabel(10, "V0ProperLifetimeK0s");
    xSec->SetBinLabel(11, "V0RapidityLambda");
    xSec->SetBinLabel(12, "V0ProperLifetimeLambda");
    xSec->SetBinLabel(13, "V0RapidityAntiLambda");
    xSec->SetBinLabel(14, "V0ProperLifetimeAntiLambda");

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

    if (ptrack.tpcInnerParam() > tpcCut) {
      if (!ptrack.hasTOF())
        return false;
      if (std::abs(ptrack.tofNSigmaPi()) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > tpcCut) {
      if (!ntrack.hasTOF())
        return false;
      if (std::abs(ntrack.tofNSigmaPi()) > nsigmaTOFmax)
        return false;
    }

    // Invariant-Mass Selection
    if (v0.mK0Short() < minMassK0s || v0.mK0Short() > maxMassK0s)
      return false;

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

    if (ptrack.tpcInnerParam() > tpcCut) {
      if (!ptrack.hasTOF())
        return false;
      if (std::abs(ptrack.tofNSigmaPr()) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > tpcCut) {
      if (!ntrack.hasTOF())
        return false;
      if (std::abs(ntrack.tofNSigmaPi()) > nsigmaTOFmax)
        return false;
    }

    // Invariant-Mass Selection
    if (v0.mLambda() < minMassLambda || v0.mLambda() > maxMassLambda)
      return false;

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

    if (ptrack.tpcInnerParam() > tpcCut) {
      if (!ptrack.hasTOF())
        return false;
      if (std::abs(ptrack.tofNSigmaPi()) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > tpcCut) {
      if (!ntrack.hasTOF())
        return false;
      if (std::abs(ntrack.tofNSigmaPr()) > nsigmaTOFmax)
        return false;
    }

    // Invariant-Mass Selection
    if (v0.mAntiLambda() < minMassLambda || v0.mAntiLambda() > maxMassLambda)
      return false;

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

    if (ptrack.tpcInnerParam() > tpcCut) {
      if (!ptrack.hasTOF())
        return false;
      if (std::abs(ptrack.tofNSigmaEl()) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > tpcCut) {
      if (!ntrack.hasTOF())
        return false;
      if (std::abs(ntrack.tofNSigmaEl()) > nsigmaTOFmax)
        return false;
    }

    // Invariant-Mass Selection
    if (v0.mGamma() < minMassGamma || v0.mGamma() > maxMassGamma)
      return false;

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
  bool passedPhiVarCut(const T& trk, float magField, const TF1& fphiCutLow, const TF1& fphiCutHigh)
  {
    float pt = trk.pt();
    float phi = trk.phi();
    int charge = trk.sign();
    float eta = trk.eta();
    float sigP = trk.sign() * trk.tpcInnerParam();
    auto nTPCCl = trk.tpcNClsFound();

    if (pt < pTcut)
      return true;

    if (magField < 0) // for negatve polarity field
      phi = o2::constants::math::TwoPI - phi;
    if (charge < 0) // for negatve charge
      phi = o2::constants::math::TwoPI - phi;

    // to center gap in the middle
    phi += o2::constants::math::PI / 18.0f;
    phi = std::fmod(phi, o2::constants::math::PI / 9.0f);

    registryDeDx.fill(HIST("hpt_vs_phi_Ncl_TPC_Before"), pt, phi, nTPCCl);

    // cut phi
    if (phi < fphiCutHigh.Eval(pt) && phi > fphiCutLow.Eval(pt))
      return false; // reject track

    if (eta > EtaCut[0] && eta < EtaCut[1]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclTPCDedxMomentumNegBefore[0]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclTPCDedxMomentumPosBefore[0]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[1] && eta < EtaCut[2]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclTPCDedxMomentumNegBefore[1]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclTPCDedxMomentumPosBefore[1]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[2] && eta < EtaCut[3]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclTPCDedxMomentumNegBefore[2]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclTPCDedxMomentumPosBefore[2]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[3] && eta < EtaCut[4]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclTPCDedxMomentumNegBefore[3]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclTPCDedxMomentumPosBefore[3]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[4] && eta < EtaCut[5]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclTPCDedxMomentumNegBefore[4]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclTPCDedxMomentumPosBefore[4]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[5] && eta < EtaCut[6]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclTPCDedxMomentumNegBefore[5]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclTPCDedxMomentumPosBefore[5]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[6] && eta < EtaCut[7]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclTPCDedxMomentumNegBefore[6]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclTPCDedxMomentumPosBefore[6]), nTPCCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[7] && eta < EtaCut[8]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclTPCDedxMomentumNegBefore[7]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclTPCDedxMomentumPosBefore[7]), nTPCCl, trk.tpcSignal(), sigP);
      }
    }
    registryDeDx.fill(HIST("trackselAll"), TrkPriCutLabel::PhiVarCutPri);

    if (nTPCClCut) {
      // cut Ncl
      if (nTPCCl < minTPCnClsFound)
        return false;

      registryDeDx.fill(HIST("hpt_vs_phi_Ncl_TPC_After"), pt, phi, nTPCCl);

      if (eta > EtaCut[0] && eta < EtaCut[1]) {
        if (sigP < 0) {
          registryDeDx.fill(HIST(kNclTPCDedxMomentumNegAfter[0]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
        } else {
          registryDeDx.fill(HIST(kNclTPCDedxMomentumPosAfter[0]), nTPCCl, trk.tpcSignal(), sigP);
        }
      } else if (eta > EtaCut[1] && eta < EtaCut[2]) {
        if (sigP < 0) {
          registryDeDx.fill(HIST(kNclTPCDedxMomentumNegAfter[1]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
        } else {
          registryDeDx.fill(HIST(kNclTPCDedxMomentumPosAfter[1]), nTPCCl, trk.tpcSignal(), sigP);
        }
      } else if (eta > EtaCut[2] && eta < EtaCut[3]) {
        if (sigP < 0) {
          registryDeDx.fill(HIST(kNclTPCDedxMomentumNegAfter[2]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
        } else {
          registryDeDx.fill(HIST(kNclTPCDedxMomentumPosAfter[2]), nTPCCl, trk.tpcSignal(), sigP);
        }
      } else if (eta > EtaCut[3] && eta < EtaCut[4]) {
        if (sigP < 0) {
          registryDeDx.fill(HIST(kNclTPCDedxMomentumNegAfter[3]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
        } else {
          registryDeDx.fill(HIST(kNclTPCDedxMomentumPosAfter[3]), nTPCCl, trk.tpcSignal(), sigP);
        }
      } else if (eta > EtaCut[4] && eta < EtaCut[5]) {
        if (sigP < 0) {
          registryDeDx.fill(HIST(kNclTPCDedxMomentumNegAfter[4]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
        } else {
          registryDeDx.fill(HIST(kNclTPCDedxMomentumPosAfter[4]), nTPCCl, trk.tpcSignal(), sigP);
        }
      } else if (eta > EtaCut[5] && eta < EtaCut[6]) {
        if (sigP < 0) {
          registryDeDx.fill(HIST(kNclTPCDedxMomentumNegAfter[5]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
        } else {
          registryDeDx.fill(HIST(kNclTPCDedxMomentumPosAfter[5]), nTPCCl, trk.tpcSignal(), sigP);
        }
      } else if (eta > EtaCut[6] && eta < EtaCut[7]) {
        if (sigP < 0) {
          registryDeDx.fill(HIST(kNclTPCDedxMomentumNegAfter[6]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
        } else {
          registryDeDx.fill(HIST(kNclTPCDedxMomentumPosAfter[6]), nTPCCl, trk.tpcSignal(), sigP);
        }
      } else if (eta > EtaCut[7] && eta < EtaCut[8]) {
        if (sigP < 0) {
          registryDeDx.fill(HIST(kNclTPCDedxMomentumNegAfter[7]), nTPCCl, trk.tpcSignal(), std::abs(sigP));
        } else {
          registryDeDx.fill(HIST(kNclTPCDedxMomentumPosAfter[7]), nTPCCl, trk.tpcSignal(), sigP);
        }
      }
    }
    registryDeDx.fill(HIST("trackselAll"), TrkPriCutLabel::NTPCClCutPri);
    return true;
  }

  // NclCutITS
  template <typename T>
  bool passedNITSClCut(const T& trk)
  {
    float eta = trk.eta();
    float sigP = trk.sign() * trk.tpcInnerParam();
    auto nITSCl = trk.itsNCls();

    if (eta > EtaCut[0] && eta < EtaCut[1]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclITSDedxMomentumNegBefore[0]), nITSCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclITSDedxMomentumPosBefore[0]), nITSCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[1] && eta < EtaCut[2]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclITSDedxMomentumNegBefore[1]), nITSCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclITSDedxMomentumPosBefore[1]), nITSCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[2] && eta < EtaCut[3]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclITSDedxMomentumNegBefore[2]), nITSCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclITSDedxMomentumPosBefore[2]), nITSCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[3] && eta < EtaCut[4]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclITSDedxMomentumNegBefore[3]), nITSCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclITSDedxMomentumPosBefore[3]), nITSCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[4] && eta < EtaCut[5]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclITSDedxMomentumNegBefore[4]), nITSCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclITSDedxMomentumPosBefore[4]), nITSCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[5] && eta < EtaCut[6]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclITSDedxMomentumNegBefore[5]), nITSCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclITSDedxMomentumPosBefore[5]), nITSCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[6] && eta < EtaCut[7]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclITSDedxMomentumNegBefore[6]), nITSCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclITSDedxMomentumPosBefore[6]), nITSCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[7] && eta < EtaCut[8]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclITSDedxMomentumNegBefore[7]), nITSCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclITSDedxMomentumPosBefore[7]), nITSCl, trk.tpcSignal(), sigP);
      }
    }

    if (nITSCl < minITSnCls)
      return false;

    if (eta > EtaCut[0] && eta < EtaCut[1]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclITSDedxMomentumNegAfter[0]), nITSCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclITSDedxMomentumPosAfter[0]), nITSCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[1] && eta < EtaCut[2]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclITSDedxMomentumNegAfter[1]), nITSCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclITSDedxMomentumPosAfter[1]), nITSCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[2] && eta < EtaCut[3]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclITSDedxMomentumNegAfter[2]), nITSCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclITSDedxMomentumPosAfter[2]), nITSCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[3] && eta < EtaCut[4]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclITSDedxMomentumNegAfter[3]), nITSCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclITSDedxMomentumPosAfter[3]), nITSCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[4] && eta < EtaCut[5]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclITSDedxMomentumNegAfter[4]), nITSCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclITSDedxMomentumPosAfter[4]), nITSCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[5] && eta < EtaCut[6]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclITSDedxMomentumNegAfter[5]), nITSCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclITSDedxMomentumPosAfter[5]), nITSCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[6] && eta < EtaCut[7]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclITSDedxMomentumNegAfter[6]), nITSCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclITSDedxMomentumPosAfter[6]), nITSCl, trk.tpcSignal(), sigP);
      }
    } else if (eta > EtaCut[7] && eta < EtaCut[8]) {
      if (sigP < 0) {
        registryDeDx.fill(HIST(kNclITSDedxMomentumNegAfter[7]), nITSCl, trk.tpcSignal(), std::abs(sigP));
      } else {
        registryDeDx.fill(HIST(kNclITSDedxMomentumPosAfter[7]), nITSCl, trk.tpcSignal(), sigP);
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
  bool passedNTPCClCutSecondaries(const T& trk)
  {
    auto nTPCCl = trk.tpcNClsFound();

    if (nTPCCl < minTPCnClsFound)
      return false;

    return true;
  }

  // NclCutITS primary
  template <typename T>
  bool passedNITSClCutSecondaries(const T& trk)
  {
    auto nITSCl = trk.itsNCls();

    if (nITSCl < minITSnCls)
      return false;

    return true;
  }

  // Process Data
  void process(SelectedCollisions::iterator const& collision,
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

    // Event Counter
    registryDeDx.fill(HIST("histRecVtxZData"), collision.posZ());

    // For magnetic field
    const auto& foundBC = collision.foundBC_as<BCsRun3>();
    const uint64_t timeStamp{foundBC.timestamp()};
    const int magField{getMagneticField(timeStamp)};

    for (const auto& trk : tracks) {
      registryDeDx.fill(HIST("trackselAll"), TrkPriCutLabel::AllPri);
      // track Selection
      if (!mySelectionPrim.IsSelected(trk))
        continue;

      registryDeDx.fill(HIST("trackselAll"), TrkPriCutLabel::SelectionPrim);

      // phi and Ncl cut
      if (phiVarCut) {
        if (!passedPhiVarCut(trk, magField, *fphiCutLow, *fphiCutHigh))
          continue;
      }

      // NCl cut ITS
      if (nITSClCut) {
        if (!passedNITSClCut(trk))
          continue;
      }
      registryDeDx.fill(HIST("trackselAll"), TrkPriCutLabel::NITSClCutPri);

      float signedP = trk.sign() * trk.tpcInnerParam();

      // MIP calibration for pions
      if (trk.tpcInnerParam() >= pionMin && trk.tpcInnerParam() <= pionMax) {
        if (calibrationMode) {
          if (signedP < 0) {
            registryDeDx.fill(HIST("hdEdx_vs_eta_Neg_Pi"), trk.eta(), trk.tpcSignal());
          } else {
            registryDeDx.fill(HIST("hdEdx_vs_eta_Pos_Pi"), trk.eta(), trk.tpcSignal());
          }

        } else {
          for (int i = 0; i < kEtaIntervals; ++i) {
            if (trk.eta() > EtaCut[i] && trk.eta() < EtaCut[i + 1]) {
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
          for (int i = 0; i < kEtaIntervals; ++i) {
            if (trk.eta() > EtaCut[i] && trk.eta() < EtaCut[i + 1]) {
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
          for (int i = 0; i < kEtaIntervals; ++i) {
            if (trk.eta() > EtaCut[i] && trk.eta() < EtaCut[i + 1]) {
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
        for (int i = 0; i < kEtaIntervals; ++i) {
          if (trk.eta() > EtaCut[i] && trk.eta() < EtaCut[i + 1]) {
            if (signedP > 0) {
              registryDeDx.fill(HIST(kDedxvsMomentumPos[0]), signedP, trk.tpcSignal() * 50 / calibrationFactorPos->at(i), trk.eta());
              registryDeDx.fill(HIST("heta_vs_p_vs_pt_all_Pos"), trk.eta(), trk.pt(), signedP);
            } else {
              registryDeDx.fill(HIST(kDedxvsMomentumNeg[0]), std::abs(signedP), trk.tpcSignal() * 50 / calibrationFactorNeg->at(i), trk.eta());
              registryDeDx.fill(HIST("heta_vs_p_vs_pt_all_Neg"), trk.eta(), trk.pt(), std::abs(signedP));
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
        }
        registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel ::PhiVarCutSec);

        if (nTPCClCut) {
          if (!passedNTPCClCutSecondaries(posTrack))
            continue;

          if (!passedNTPCClCutSecondaries(negTrack))
            continue;
        }

        registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel ::NTPCClCutSec);

        if (nITSClCut) {
          if (!passedNITSClCutSecondaries(posTrack))
            continue;

          if (!passedNITSClCutSecondaries(negTrack))
            continue;
        }

        registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel ::NITSClCutSec);

        float signedPpos = posTrack.sign() * posTrack.tpcInnerParam();
        float signedPneg = negTrack.sign() * negTrack.tpcInnerParam();

        float pxPos = posTrack.px();
        float pyPos = posTrack.py();
        float pzPos = posTrack.pz();

        float pxNeg = negTrack.px();
        float pyNeg = negTrack.py();
        float pzNeg = negTrack.pz();

        const float gammaMass = 2 * MassElectron; // GeV/c^2

        // K0s Selection
        if (passedK0Selection(v0, negTrack, posTrack, collision)) {

          if (std::abs(v0.rapidity(MassK0Short)) > v0rapidityCut)
            continue;

          registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::V0RapidityK0s);
          float properLifetime = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * MassK0Short;

          if (properLifetime > v0ProperLifetimeCutK0s)
            continue;

          registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::V0ProperLifetimeK0s);

          float ePosPi = posTrack.energy(MassPionCharged);
          float eNegPi = negTrack.energy(MassPionCharged);

          float invMass = std::sqrt((eNegPi + ePosPi) * (eNegPi + ePosPi) - ((pxNeg + pxPos) * (pxNeg + pxPos) + (pyNeg + pyPos) * (pyNeg + pyPos) + (pzNeg + pzPos) * (pzNeg + pzPos)));

          if (std::abs(invMass - MassK0Short) > invMassCut) {
            continue;
          }

          for (int i = 0; i < kEtaIntervals; ++i) {
            if (negTrack.eta() > EtaCut[i] && negTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumNeg[1]), std::abs(signedPneg), negTrack.tpcSignal() * 50 / calibrationFactorNeg->at(i), negTrack.eta());
            }
            if (posTrack.eta() > EtaCut[i] && posTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumPos[1]), signedPpos, posTrack.tpcSignal() * 50 / calibrationFactorPos->at(i), posTrack.eta());
            }
          }
        }

        // Lambda Selection
        if (passedLambdaSelection(v0, negTrack, posTrack, collision)) {

          if (std::abs(v0.rapidity(MassLambda)) > v0rapidityCut)
            continue;

          registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::V0RapidityLambda);
          float properLifetime = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * MassLambda;

          if (properLifetime > v0ProperLifetimeCutLambda)
            continue;

          registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::V0ProperLifetimeLambda);

          float ePosPr = posTrack.energy(MassProton);
          float eNegPi = negTrack.energy(MassPionCharged);

          float invMass = std::sqrt((eNegPi + ePosPr) * (eNegPi + ePosPr) - ((pxNeg + pxPos) * (pxNeg + pxPos) + (pyNeg + pyPos) * (pyNeg + pyPos) + (pzNeg + pzPos) * (pzNeg + pzPos)));

          if (std::abs(invMass - MassLambda) > invMassCut) {
            continue;
          }

          for (int i = 0; i < kEtaIntervals; ++i) {
            if (negTrack.eta() > EtaCut[i] && negTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumNeg[1]), std::abs(signedPneg), negTrack.tpcSignal() * 50 / calibrationFactorNeg->at(i), negTrack.eta());
            }
            if (posTrack.eta() > EtaCut[i] && posTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumPos[2]), signedPpos, posTrack.tpcSignal() * 50 / calibrationFactorPos->at(i), posTrack.eta());
            }
          }
        }

        // AntiLambda Selection
        if (passedAntiLambdaSelection(v0, negTrack, posTrack, collision)) {

          if (std::abs(v0.rapidity(MassLambda)) > v0rapidityCut)
            continue;

          registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::V0RapidityAntiLambda);
          float properLifetime = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * MassLambda;

          if (properLifetime > v0ProperLifetimeCutLambda)
            continue;

          registryDeDx.fill(HIST("trackselSec"), TrkSecCutLabel::V0ProperLifetimeAntiLambda);

          float ePosPi = posTrack.energy(MassPionCharged);
          float eNegPr = negTrack.energy(MassProton);

          float invMass = std::sqrt((eNegPr + ePosPi) * (eNegPr + ePosPi) - ((pxNeg + pxPos) * (pxNeg + pxPos) + (pyNeg + pyPos) * (pyNeg + pyPos) + (pzNeg + pzPos) * (pzNeg + pzPos)));

          if (std::abs(invMass - MassLambda) > invMassCut) {
            continue;
          }

          for (int i = 0; i < kEtaIntervals; ++i) {
            if (negTrack.eta() > EtaCut[i] && negTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumNeg[2]), std::abs(signedPneg), negTrack.tpcSignal() * 50 / calibrationFactorNeg->at(i), negTrack.eta());
            }
            if (posTrack.eta() > EtaCut[i] && posTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumPos[1]), signedPpos, posTrack.tpcSignal() * 50 / calibrationFactorPos->at(i), posTrack.eta());
            }
          }
        }

        // Gamma Selection
        if (passedGammaSelection(v0, negTrack, posTrack, collision)) {

          float ePosEl = posTrack.energy(MassElectron);
          float eNegEl = negTrack.energy(MassElectron);

          float invMass = std::sqrt((eNegEl + ePosEl) * (eNegEl + ePosEl) - ((pxNeg + pxPos) * (pxNeg + pxPos) + (pyNeg + pyPos) * (pyNeg + pyPos) + (pzNeg + pzPos) * (pzNeg + pzPos)));

          if (std::abs(invMass - gammaMass) > invMassCutGamma) {
            continue;
          }

          for (int i = 0; i < kEtaIntervals; ++i) {
            if (negTrack.eta() > EtaCut[i] && negTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumNeg[3]), std::abs(signedPneg), negTrack.tpcSignal() * 50 / calibrationFactorNeg->at(i), negTrack.eta());
            }
            if (posTrack.eta() > EtaCut[i] && posTrack.eta() < EtaCut[i + 1]) {
              registryDeDx.fill(HIST(kDedxvsMomentumPos[3]), signedPpos, posTrack.tpcSignal() * 50 / calibrationFactorPos->at(i), posTrack.eta());
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<DedxAnalysis>(cfgc)};
}
