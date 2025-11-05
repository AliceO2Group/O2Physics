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

/// \file   flowTask.cxx
/// \author Zhiyong Lu (zhiyong.lu@cern.ch)
/// \since  Dec/10/2023
/// \brief  jira: PWGCF-254, task to measure flow observables with cumulant method

#include "FlowContainer.h"
#include "FlowPtContainer.h"
#include "GFW.h"
#include "GFWConfig.h"
#include "GFWCumulant.h"
#include "GFWPowerArray.h"
#include "GFWWeights.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPMagField.h>

#include "TList.h"
#include <TF1.h>
#include <TObjArray.h>
#include <TProfile.h>
#include <TRandom3.h>

#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowTask {

  // Basic event&track selections
  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCentEstimator, int, 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FT0A")
  O2_DEFINE_CONFIGURABLE(cfgCentFT0CMin, float, 0.0f, "Minimum centrality (FT0C) to cut events in filter")
  O2_DEFINE_CONFIGURABLE(cfgCentFT0CMax, float, 100.0f, "Maximum centrality (FT0C) to cut events in filter")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtRefMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtRefMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for all tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 10.0f, "Maximal pT for all tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgEtaPtPt, float, 0.4, "eta cut for pt-pt correlations");
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5f, "max chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCclu, float, 50.0f, "minimum TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCCrossedRows, float, 70.0f, "minimum TPC crossed rows")
  O2_DEFINE_CONFIGURABLE(cfgCutITSclu, float, 5.0f, "minimum ITS clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "max DCA to vertex z")
  // Additional events selection flags
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, false, "Use additional event cut on mult correlations")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoSameBunchPileup, bool, false, "rejects collisions which are associated with the same found-by-T0 bunch crossing")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoITSROFrameBorder, bool, false, "reject events at ITS ROF border")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoTimeFrameBorder, bool, false, "reject events at TF border")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkIsGoodZvtxFT0vsPV, bool, false, "removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference, use this cut at low multiplicities with caution")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoCollInTimeRangeStandard, bool, false, "no collisions in specified time range")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkIsGoodITSLayersAll, bool, true, "cut time intervals with dead ITS staves")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoCollInRofStandard, bool, false, "no other collisions in this Readout Frame with per-collision multiplicity above threshold")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoHighMultCollInPrevRof, bool, false, "veto an event if FT0C amplitude in previous ITS ROF is above threshold")
  O2_DEFINE_CONFIGURABLE(cfgEvSelMultCorrelation, bool, true, "Multiplicity correlation cut")
  O2_DEFINE_CONFIGURABLE(cfgEvSelV0AT0ACut, bool, true, "V0A T0A 5 sigma cut")
  O2_DEFINE_CONFIGURABLE(cfgGetInteractionRate, bool, false, "Get interaction rate from CCDB")
  O2_DEFINE_CONFIGURABLE(cfgUseInteractionRateCut, bool, false, "Use events with low interaction rate")
  O2_DEFINE_CONFIGURABLE(cfgCutMaxIR, float, 50.0f, "maximum interaction rate (kHz)")
  O2_DEFINE_CONFIGURABLE(cfgCutMinIR, float, 0.0f, "minimum interaction rate (kHz)")
  O2_DEFINE_CONFIGURABLE(cfgEvSelOccupancy, bool, true, "Occupancy cut")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyHigh, int, 500, "High cut on TPC occupancy")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyLow, int, 0, "Low cut on TPC occupancy")
  // User configuration for ananlysis
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Use Nch for flow observables")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 30, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeights, bool, false, "Fill and output NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeightsRefPt, bool, false, "NUA weights are filled in ref pt bins")
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeightsRunbyRun, bool, false, "NUA weights are filled run-by-run")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgUseSmallMemory, bool, false, "Use small memory mode")
  O2_DEFINE_CONFIGURABLE(cfgTrackDensityCorrUse, bool, false, "Use track density efficiency correction")
  O2_DEFINE_CONFIGURABLE(cfgUseCentralMoments, bool, true, "Use central moments in vn-pt calculations")
  O2_DEFINE_CONFIGURABLE(cfgUsePtRef, bool, true, "Use refernce pt range for pt container (if you are checking the spectra, you need to extent it)")
  O2_DEFINE_CONFIGURABLE(cfgMpar, int, 4, "Highest order of pt-pt correlations")
  Configurable<std::vector<std::string>> cfgUserDefineGFWCorr{"cfgUserDefineGFWCorr", std::vector<std::string>{"refN02 {2} refP02 {-2}", "refN12 {2} refP12 {-2}"}, "User defined GFW CorrelatorConfig"};
  Configurable<std::vector<std::string>> cfgUserDefineGFWName{"cfgUserDefineGFWName", std::vector<std::string>{"Ch02Gap22", "Ch12Gap22"}, "User defined GFW Name"};
  Configurable<GFWCorrConfigs> cfgUserPtVnCorrConfig{"cfgUserPtVnCorrConfig", {{"refP {2} refN {-2}", "refP {3} refN {-3}"}, {"ChGap22", "ChGap32"}, {0, 0}, {3, 3}}, "Configurations for vn-pt correlations"};
  Configurable<std::vector<int>> cfgRunRemoveList{"cfgRunRemoveList", std::vector<int>{-1}, "excluded run numbers"};
  struct : ConfigurableGroup {
    Configurable<std::vector<double>> cfgTrackDensityP0{"cfgTrackDensityP0", std::vector<double>{0.7217476707, 0.7384792571, 0.7542625668, 0.7640680200, 0.7701951667, 0.7755299053, 0.7805901710, 0.7849446786, 0.7957356586, 0.8113039262, 0.8211968966, 0.8280558878, 0.8329342135}, "parameter 0 for track density efficiency correction"};
    Configurable<std::vector<double>> cfgTrackDensityP1{"cfgTrackDensityP1", std::vector<double>{-2.169488e-05, -2.191913e-05, -2.295484e-05, -2.556538e-05, -2.754463e-05, -2.816832e-05, -2.846502e-05, -2.843857e-05, -2.705974e-05, -2.477018e-05, -2.321730e-05, -2.203315e-05, -2.109474e-05}, "parameter 1 for track density efficiency correction"};
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
    // for TPC sector boundary
    O2_DEFINE_CONFIGURABLE(cfgShowTPCsectorOverlap, bool, true, "Draw TPC sector overlap")
    O2_DEFINE_CONFIGURABLE(cfgRejectionTPCsectorOverlap, bool, false, "rejection for TPC sector overlap")
    O2_DEFINE_CONFIGURABLE(cfgMagnetField, std::string, "GLO/Config/GRPMagField", "CCDB path to Magnet field object")
    ConfigurableAxis axisPhiMod{"axisPhiMod", {100, 0, constants::math::PI / 9}, "fmod(#varphi,#pi/9)"};
    O2_DEFINE_CONFIGURABLE(cfgTPCPhiCutLowCutFunction, std::string, "0.1/x-0.005", "Function for TPC mod phi-pt cut");
    O2_DEFINE_CONFIGURABLE(cfgTPCPhiCutHighCutFunction, std::string, "0.1/x+0.01", "Function for TPC mod phi-pt cut");
    O2_DEFINE_CONFIGURABLE(cfgTPCPhiCutPtMin, float, 2.0f, "start point of phi-pt cut")
    TF1* fPhiCutLow = nullptr;
    TF1* fPhiCutHigh = nullptr;
    // for deltaPt/<pT> vs c22 vs centrality
    O2_DEFINE_CONFIGURABLE(cfgDptDisEnable, bool, false, "Produce deltaPt/meanPt vs c22 vs centrality")
    O2_DEFINE_CONFIGURABLE(cfgDptDisCorrConfigIndex, int, 7, "Index in CorrConfig, this decide which correlation is filled")
    TH1D* hEvAvgMeanPt = nullptr;
    O2_DEFINE_CONFIGURABLE(cfgDptDishEvAvgMeanPt, std::string, "", "CCDB path to hMeanPt object")
    ConfigurableAxis cfgDptDisAxisCorr{"cfgDptDisAxisCorr", {50, 0., 0.005}, "pt axis for histograms"};
    ConfigurableAxis cfgDptDisAxisNormal{"cfgDptDisAxisNormal", {200, -1., 1.}, "normalized axis"};
  } cfgFuncParas;

  ConfigurableAxis axisPtHist{"axisPtHist", {100, 0., 10.}, "pt axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10}, "pt axis for histograms"};
  ConfigurableAxis axisIndependent{"axisIndependent", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "X axis for histograms"};
  ConfigurableAxis axisNch{"axisNch", {4000, 0, 4000}, "N_{ch}"};
  ConfigurableAxis axisDCAz{"axisDCAz", {200, -2, 2}, "DCA_{z} (cm)"};
  ConfigurableAxis axisDCAxy{"axisDCAxy", {200, -1, 1}, "DCA_{xy} (cm)"};

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgCutVertex) && (aod::cent::centFT0C > cfgCentFT0CMin) && (aod::cent::centFT0C < cfgCentFT0CMax);
  Filter trackFilter = ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  using FilteredCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::Mults>>;
  using FilteredTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA>>;
  // Filter for MCcollisions
  Filter mccollisionFilter = nabs(aod::mccollision::posZ) < cfgCutVertex;
  using FilteredMcCollisions = soa::Filtered<aod::McCollisions>;
  // Filter for MCParticle
  Filter particleFilter = (nabs(aod::mcparticle::eta) < cfgCutEta) && (aod::mcparticle::pt > cfgCutPtMin) && (aod::mcparticle::pt < cfgCutPtMax);
  using FilteredMcParticles = soa::Filtered<aod::McParticles>;

  using FilteredSmallGroupMcCollisions = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSel, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::Mults>>;

  // Corrections
  TH1D* mEfficiency = nullptr;
  GFWWeights* mAcceptance = nullptr;
  bool correctionsLoaded = false;

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  // Define output
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  OutputObj<FlowContainer> fFCgen{FlowContainer("FlowContainer_gen")};
  OutputObj<FlowPtContainer> fFCpt{FlowPtContainer("FlowPtContainer")};
  OutputObj<FlowPtContainer> fFCptgen{FlowPtContainer("FlowPtContainer_gen")};
  OutputObj<GFWWeights> fWeights{GFWWeights("weights")};
  HistogramRegistry registry{"registry"};

  // define global variables
  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  GFWCorrConfigs gfwConfigs;
  std::vector<GFW::CorrConfig> corrconfigsPtVn;
  TAxis* fPtAxis;
  TRandom3* fRndm = new TRandom3(0);
  int lastRunNumber = -1;
  std::vector<int> runNumbers;
  std::map<int, std::shared_ptr<TH3>> th3sPerRun; // map of TH3 histograms for all runs
  enum CentEstimators {
    kCentFT0C = 0,
    kCentFT0CVariant1,
    kCentFT0M,
    kCentFV0A,
    // Count the total number of enum
    kCount_CentEstimators
  };
  enum DataType {
    kReco,
    kGen
  };
  int mRunNumber{-1};
  uint64_t mSOR{0};
  double mMinSeconds{-1.};
  std::unordered_map<int, TH2*> gHadronicRate;
  ctpRateFetcher mRateFetcher;
  TH2* gCurrentHadronicRate;

  // phi-EP correction
  std::vector<TF1*> funcEff;
  TH1D* hFindPtBin;
  TF1* funcV2;
  TF1* funcV3;
  TF1* funcV4;

  void init(InitContext const&)
  {
    const AxisSpec axisVertex{40, -20, 20, "Vtxz (cm)"};
    const AxisSpec axisPhi{60, 0.0, constants::math::TwoPI, "#varphi"};
    const AxisSpec axisEta{40, -1., 1., "#eta"};
    const AxisSpec axisCentForQA{100, 0, 100, "centrality (%)"};
    const AxisSpec axisT0C{70, 0, 70000, "N_{ch} (T0C)"};
    const AxisSpec axisT0A{200, 0, 200000, "N_{ch} (T0A)"};

    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    // Add some output objects to the histogram registry
    // Event QA
    registry.add("hEventCount", "Number of Event;; Count", {HistType::kTH1D, {{5, 0, 5}}});
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(1, "Filtered event");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(2, "after sel8");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(3, "after supicious Runs removal");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(4, "after additional event cut");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(5, "after correction loads");
    registry.add("hEventCountSpecific", "Number of Event;; Count", {HistType::kTH1D, {{12, 0, 12}}});
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(1, "after sel8");
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(2, "kNoSameBunchPileup");
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(3, "kNoITSROFrameBorder");
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(4, "kNoTimeFrameBorder");
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(5, "kIsGoodZvtxFT0vsPV");
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(6, "kNoCollInTimeRangeStandard");
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(7, "kIsGoodITSLayersAll");
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(8, "kNoCollInRofStandard");
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(9, "kNoHighMultCollInPrevRof");
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(10, "occupancy");
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(11, "MultCorrelation");
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(12, "cfgEvSelV0AT0ACut");
    registry.add("hVtxZ", "Vexter Z distribution", {HistType::kTH1D, {axisVertex}});
    registry.add("hMult", "Multiplicity distribution", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    std::string hCentTitle = "Centrality distribution, Estimator " + std::to_string(cfgCentEstimator);
    registry.add("hCent", hCentTitle.c_str(), {HistType::kTH1D, {{100, 0, 100}}});
    if (doprocessMCGen) {
      registry.add("MCGen/MChVtxZ", "Vexter Z distribution", {HistType::kTH1D, {axisVertex}});
      registry.add("MCGen/MChMult", "Multiplicity distribution", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
      registry.add("MCGen/MChCent", hCentTitle.c_str(), {HistType::kTH1D, {{100, 0, 100}}});
    }
    if (!cfgUseSmallMemory) {
      registry.add("BeforeSel8_globalTracks_centT0C", "before sel8;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {axisCentForQA, axisNch}});
      registry.add("BeforeCut_globalTracks_centT0C", "before cut;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {axisCentForQA, axisNch}});
      registry.add("BeforeCut_PVTracks_centT0C", "before cut;Centrality T0C;mulplicity PV tracks", {HistType::kTH2D, {axisCentForQA, axisNch}});
      registry.add("BeforeCut_globalTracks_PVTracks", "before cut;mulplicity PV tracks;mulplicity global tracks", {HistType::kTH2D, {axisNch, axisNch}});
      registry.add("BeforeCut_globalTracks_multT0A", "before cut;mulplicity T0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisNch}});
      registry.add("BeforeCut_globalTracks_multV0A", "before cut;mulplicity V0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisNch}});
      registry.add("BeforeCut_multV0A_multT0A", "before cut;mulplicity T0A;mulplicity V0A", {HistType::kTH2D, {axisT0A, axisT0A}});
      registry.add("BeforeCut_multT0C_centT0C", "before cut;Centrality T0C;mulplicity T0C", {HistType::kTH2D, {axisCentForQA, axisT0C}});
      registry.add("globalTracks_centT0C", "after cut;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {axisCentForQA, axisNch}});
      registry.add("PVTracks_centT0C", "after cut;Centrality T0C;mulplicity PV tracks", {HistType::kTH2D, {axisCentForQA, axisNch}});
      registry.add("globalTracks_PVTracks", "after cut;mulplicity PV tracks;mulplicity global tracks", {HistType::kTH2D, {axisNch, axisNch}});
      registry.add("globalTracks_multT0A", "after cut;mulplicity T0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisNch}});
      registry.add("globalTracks_multV0A", "after cut;mulplicity V0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisNch}});
      registry.add("multV0A_multT0A", "after cut;mulplicity T0A;mulplicity V0A", {HistType::kTH2D, {axisT0A, axisT0A}});
      registry.add("multT0C_centT0C", "after cut;Centrality T0C;mulplicity T0C", {HistType::kTH2D, {axisCentForQA, axisT0C}});
      registry.add("centFT0CVar_centFT0C", "after cut;Centrality T0C;Centrality T0C Var", {HistType::kTH2D, {axisCentForQA, axisCentForQA}});
      registry.add("centFT0M_centFT0C", "after cut;Centrality T0C;Centrality T0M", {HistType::kTH2D, {axisCentForQA, axisCentForQA}});
      registry.add("centFV0A_centFT0C", "after cut;Centrality T0C;Centrality V0A", {HistType::kTH2D, {axisCentForQA, axisCentForQA}});
    }
    // Track QA
    registry.add("hPhi", "#phi distribution", {HistType::kTH1D, {axisPhi}});
    registry.add("hPhiWeighted", "corrected #phi distribution", {HistType::kTH1D, {axisPhi}});
    registry.add("hEta", "#eta distribution", {HistType::kTH1D, {axisEta}});
    registry.add("hPt", "p_{T} distribution before cut", {HistType::kTH1D, {axisPtHist}});
    registry.add("hPtRef", "p_{T} distribution after cut", {HistType::kTH1D, {axisPtHist}});
    registry.add("pt_phi_bef", "before cut;p_{T};#phi_{modn}", {HistType::kTH2D, {axisPt, cfgFuncParas.axisPhiMod}});
    registry.add("pt_phi_aft", "after cut;p_{T};#phi_{modn}", {HistType::kTH2D, {axisPt, cfgFuncParas.axisPhiMod}});
    registry.add("hChi2prTPCcls", "#chi^{2}/cluster for the TPC track segment", {HistType::kTH1D, {{100, 0., 5.}}});
    registry.add("hChi2prITScls", "#chi^{2}/cluster for the ITS track", {HistType::kTH1D, {{100, 0., 50.}}});
    registry.add("hnTPCClu", "Number of found TPC clusters", {HistType::kTH1D, {{100, 40, 180}}});
    registry.add("hnITSClu", "Number of found ITS clusters", {HistType::kTH1D, {{100, 0, 20}}});
    registry.add("hnTPCCrossedRow", "Number of crossed TPC Rows", {HistType::kTH1D, {{100, 40, 180}}});
    registry.add("hDCAz", "DCAz after cuts; DCAz (cm); Pt", {HistType::kTH2D, {{200, -0.5, 0.5}, {200, 0, 5}}});
    registry.add("hDCAxy", "DCAxy after cuts; DCAxy (cm); Pt", {HistType::kTH2D, {{200, -0.5, 0.5}, {200, 0, 5}}});
    registry.add("hTrackCorrection2d", "Correlation table for number of tracks table; uncorrected track; corrected track", {HistType::kTH2D, {axisNch, axisNch}});
    registry.add("hMeanPt", "", {HistType::kTProfile, {axisIndependent}});
    registry.add("hMeanPtWithinGap08", "", {HistType::kTProfile, {axisIndependent}});
    registry.add("c22_gap08_Weff", "", {HistType::kTProfile, {axisIndependent}});
    registry.add("c22_gap08_trackMeanPt", "", {HistType::kTProfile, {axisIndependent}});
    registry.add("PtVariance_partA_WithinGap08", "", {HistType::kTProfile, {axisIndependent}});
    registry.add("PtVariance_partB_WithinGap08", "", {HistType::kTProfile, {axisIndependent}});
    if (cfgFuncParas.cfgDptDisEnable) {
      registry.add("hNormDeltaPt_Corr_X", " #delta p_{T}/[p_{T}]; Corr; X", {HistType::kTH3D, {cfgFuncParas.cfgDptDisAxisNormal, cfgFuncParas.cfgDptDisAxisCorr, axisIndependent}});
    }
    if (doprocessMCGen) {
      registry.add("MCGen/MChPhi", "#phi distribution", {HistType::kTH1D, {axisPhi}});
      registry.add("MCGen/MChEta", "#eta distribution", {HistType::kTH1D, {axisEta}});
      registry.add("MCGen/MChPtRef", "p_{T} distribution after cut", {HistType::kTH1D, {axisPtHist}});
    }

    o2::framework::AxisSpec axis = axisPt;
    int nPtBins = axis.binEdges.size() - 1;
    double* ptBins = &(axis.binEdges)[0];
    fPtAxis = new TAxis(nPtBins, ptBins);

    if (cfgOutputNUAWeights) {
      fWeights->setPtBins(nPtBins, ptBins);
      fWeights->init(true, false);
    }

    // add in FlowContainer to Get boostrap sample automatically
    TObjArray* oba = new TObjArray();
    oba->Add(new TNamed("ChGap22", "ChGap22"));
    oba->Add(new TNamed("ChFull22", "ChFull22"));
    oba->Add(new TNamed("ChFull32", "ChFull32"));
    oba->Add(new TNamed("ChFull42", "ChFull42"));
    oba->Add(new TNamed("ChFull24", "ChFull24"));
    oba->Add(new TNamed("ChFull26", "ChFull26"));
    for (auto i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("ChFull22_pt_%i", i + 1), "ChFull22_pTDiff"));
    for (auto i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("ChFull24_pt_%i", i + 1), "ChFull24_pTDiff"));
    for (auto i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("ChFull26_pt_%i", i + 1), "ChFull26_pTDiff"));
    oba->Add(new TNamed("Ch04Gap22", "Ch04Gap22"));
    oba->Add(new TNamed("Ch06Gap22", "Ch06Gap22"));
    oba->Add(new TNamed("Ch08Gap22", "Ch08Gap22"));
    oba->Add(new TNamed("Ch10Gap22", "Ch10Gap22"));
    for (auto i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("Ch10Gap22_pt_%i", i + 1), "Ch10Gap22_pTDiff"));
    oba->Add(new TNamed("Ch12Gap22", "Ch12Gap22"));
    oba->Add(new TNamed("Ch04Gap32", "Ch04Gap32"));
    oba->Add(new TNamed("Ch06Gap32", "Ch06Gap32"));
    oba->Add(new TNamed("Ch08Gap32", "Ch08Gap32"));
    oba->Add(new TNamed("Ch10Gap32", "Ch10Gap32"));
    for (auto i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("Ch10Gap32_pt_%i", i + 1), "Ch10Gap32_pTDiff"));
    oba->Add(new TNamed("Ch12Gap32", "Ch12Gap32"));
    oba->Add(new TNamed("Ch04Gap42", "Ch04Gap42"));
    oba->Add(new TNamed("Ch06Gap42", "Ch06Gap42"));
    oba->Add(new TNamed("Ch08Gap42", "Ch08Gap42"));
    oba->Add(new TNamed("Ch10Gap42", "Ch10Gap42"));
    for (auto i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("Ch10Gap42_pt_%i", i + 1), "Ch10Gap42_pTDiff"));
    oba->Add(new TNamed("Ch12Gap42", "Ch12Gap42"));
    oba->Add(new TNamed("ChFull422", "ChFull422"));
    oba->Add(new TNamed("Ch04GapA422", "Ch04GapA422"));
    oba->Add(new TNamed("Ch04GapB422", "Ch04GapB422"));
    oba->Add(new TNamed("Ch10GapA422", "Ch10GapA422"));
    oba->Add(new TNamed("Ch10GapB422", "Ch10GapB422"));
    oba->Add(new TNamed("ChFull3232", "ChFull3232"));
    oba->Add(new TNamed("ChFull4242", "ChFull4242"));
    oba->Add(new TNamed("Ch04Gap3232", "Ch04Gap3232"));
    oba->Add(new TNamed("Ch04Gap4242", "Ch04Gap4242"));
    oba->Add(new TNamed("Ch04Gap24", "Ch04Gap24"));
    oba->Add(new TNamed("Ch10Gap3232", "Ch10Gap3232"));
    oba->Add(new TNamed("Ch10Gap4242", "Ch10Gap4242"));
    oba->Add(new TNamed("Ch10Gap24", "Ch10Gap24"));
    for (auto i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("Ch10Gap24_pt_%i", i + 1), "Ch10Gap24_pTDiff"));
    std::vector<std::string> userDefineGFWCorr = cfgUserDefineGFWCorr;
    std::vector<std::string> userDefineGFWName = cfgUserDefineGFWName;
    if (!userDefineGFWCorr.empty() && !userDefineGFWName.empty()) {
      for (uint i = 0; i < userDefineGFWName.size(); i++) {
        oba->Add(new TNamed(userDefineGFWName.at(i).c_str(), userDefineGFWName.at(i).c_str()));
      }
    }
    fFC->SetName("FlowContainer");
    fFC->SetXAxis(fPtAxis);
    fFC->Initialize(oba, axisIndependent, cfgNbootstrap);
    if (doprocessMCGen) {
      fFCgen->SetName("FlowContainer_gen");
      fFCgen->SetXAxis(fPtAxis);
      fFCgen->Initialize(oba, axisIndependent, cfgNbootstrap);
    }
    delete oba;

    // eta region
    fGFW->AddRegion("full", -0.8, 0.8, 1, 1);
    fGFW->AddRegion("refN00", -0.8, 0., 1, 1);   // gap0 negative region
    fGFW->AddRegion("refP00", 0., 0.8, 1, 1);    // gap0 positve region
    fGFW->AddRegion("refN02", -0.8, -0.1, 1, 1); // gap2 negative region
    fGFW->AddRegion("refP02", 0.1, 0.8, 1, 1);   // gap2 positve region
    fGFW->AddRegion("refN04", -0.8, -0.2, 1, 1); // gap4 negative region
    fGFW->AddRegion("refP04", 0.2, 0.8, 1, 1);   // gap4 positve region
    fGFW->AddRegion("refN06", -0.8, -0.3, 1, 1); // gap6 negative region
    fGFW->AddRegion("refP06", 0.3, 0.8, 1, 1);   // gap6 positve region
    fGFW->AddRegion("refN08", -0.8, -0.4, 1, 1);
    fGFW->AddRegion("refP08", 0.4, 0.8, 1, 1);
    fGFW->AddRegion("refN10", -0.8, -0.5, 1, 1);
    fGFW->AddRegion("refP10", 0.5, 0.8, 1, 1);
    fGFW->AddRegion("refN12", -0.8, -0.6, 1, 1);
    fGFW->AddRegion("refP12", 0.6, 0.8, 1, 1);
    fGFW->AddRegion("refN14", -0.8, -0.7, 1, 1);
    fGFW->AddRegion("refP14", 0.7, 0.8, 1, 1);
    fGFW->AddRegion("refN", -0.8, -0.4, 1, 1);
    fGFW->AddRegion("refP", 0.4, 0.8, 1, 1);
    fGFW->AddRegion("refM", -0.4, 0.4, 1, 1);
    fGFW->AddRegion("poiN", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 2);
    fGFW->AddRegion("poiN10", -0.8, -0.5, 1 + fPtAxis->GetNbins(), 2);
    fGFW->AddRegion("poifull", -0.8, 0.8, 1 + fPtAxis->GetNbins(), 2);
    fGFW->AddRegion("olN", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 4);
    fGFW->AddRegion("olN10", -0.8, -0.5, 1 + fPtAxis->GetNbins(), 4);
    fGFW->AddRegion("olfull", -0.8, 0.8, 1 + fPtAxis->GetNbins(), 4);

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {3 -3}", "ChFull32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {4 -4}", "ChFull42", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 -2 -2}", "ChFull24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 2 -2 -2 -2}", "ChFull26", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN04 {2} refP04 {-2}", "Ch04Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN06 {2} refP06 {-2}", "Ch06Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "Ch08Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {2} refP10 {-2}", "Ch10Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN12 {2} refP12 {-2}", "Ch12Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN04 {3} refP04 {-3}", "Ch04Gap32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN06 {3} refP06 {-3}", "Ch06Gap32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {3} refP08 {-3}", "Ch08Gap32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {3} refP10 {-3}", "Ch10Gap32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN12 {3} refP12 {-3}", "Ch12Gap32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN04 {4} refP04 {-4}", "Ch04Gap42", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN06 {4} refP06 {-4}", "Ch06Gap42", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {4} refP08 {-4}", "Ch08Gap42", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {4} refP10 {-4}", "Ch10Gap42", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN12 {4} refP12 {-4}", "Ch12Gap42", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN {2} refP {-2}", "ChGap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poifull full | olfull {2 -2}", "ChFull22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poifull full | olfull {2 2 -2 -2}", "ChFull24", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poifull full | olfull {2 2 2 -2 -2 -2}", "ChFull26", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiN10 refN10 | olN10 {2} refP10 {-2}", "Ch10Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiN10 refN10 | olN10 {3} refP10 {-3}", "Ch10Gap32", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiN10 refN10 | olN10 {4} refP10 {-4}", "Ch10Gap42", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {4 -2 -2}", "ChFull422", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN04 {-2 -2} refP04 {4}", "Ch04GapA422", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN04 {4} refP04 {-2 -2}", "Ch04GapB422", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {-2 -2} refP10 {4}", "Ch10GapA422", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {4} refP10 {-2 -2}", "Ch10GapB422", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {3 2 -3 -2}", "ChFull3232", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {4 2 -4 -2}", "ChFull4242", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN04 {3 2} refP04 {-3 -2}", "Ch04Gap3232", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN04 {4 2} refP04 {-4 -2}", "Ch04Gap4242", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN04 {2 2} refP04 {-2 -2}", "Ch04Gap24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {3 2} refP10 {-3 -2}", "Ch10Gap3232", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {4 2} refP10 {-4 -2}", "Ch10Gap4242", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {2 2} refP10 {-2 -2}", "Ch10Gap24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiN10 refN10 | olN10 {2 2} refP10 {-2 -2}", "Ch10Gap24", kTRUE));
    if (!userDefineGFWCorr.empty() && !userDefineGFWName.empty()) {
      LOGF(info, "User adding GFW CorrelatorConfig:");
      // attentaion: here we follow the index of cfgUserDefineGFWCorr
      for (uint i = 0; i < userDefineGFWCorr.size(); i++) {
        if (i >= userDefineGFWName.size()) {
          LOGF(fatal, "The names you provided are more than configurations. userDefineGFWName.size(): %d > userDefineGFWCorr.size(): %d", userDefineGFWName.size(), userDefineGFWCorr.size());
          break;
        }
        LOGF(info, "%d: %s %s", i, userDefineGFWCorr.at(i).c_str(), userDefineGFWName.at(i).c_str());
        corrconfigs.push_back(fGFW->GetCorrelatorConfig(userDefineGFWCorr.at(i).c_str(), userDefineGFWName.at(i).c_str(), kFALSE));
      }
    }

    gfwConfigs.SetCorrs(cfgUserPtVnCorrConfig->GetCorrs());
    gfwConfigs.SetHeads(cfgUserPtVnCorrConfig->GetHeads());
    gfwConfigs.SetpTDifs(cfgUserPtVnCorrConfig->GetpTDifs());
    // Mask 1: vn-[pT], 3: vn-[pT^2], 7: vn-[pT^3], 15: vn-[pT^4]
    gfwConfigs.SetpTCorrMasks(cfgUserPtVnCorrConfig->GetpTCorrMasks());
    gfwConfigs.Print();
    fFCpt->setUseCentralMoments(cfgUseCentralMoments);
    fFCpt->setUseGapMethod(true);
    fFCpt->initialise(axisIndependent, cfgMpar, gfwConfigs, cfgNbootstrap);
    for (auto i = 0; i < gfwConfigs.GetSize(); ++i) {
      corrconfigsPtVn.push_back(fGFW->GetCorrelatorConfig(gfwConfigs.GetCorrs()[i], gfwConfigs.GetHeads()[i], gfwConfigs.GetpTDifs()[i]));
    }
    fGFW->CreateRegions();

    if (cfgEvSelMultCorrelation) {
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

      cfgFuncParas.fT0AV0AMean = new TF1("fT0AV0AMean", "[0]+[1]*x", 0, 200000);
      cfgFuncParas.fT0AV0AMean->SetParameters(-1601.0581, 9.417652e-01);
      cfgFuncParas.fT0AV0ASigma = new TF1("fT0AV0ASigma", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 200000);
      cfgFuncParas.fT0AV0ASigma->SetParameters(463.4144, 6.796509e-02, -9.097136e-07, 7.971088e-12, -2.600581e-17);
    }

    if (cfgFuncParas.cfgShowTPCsectorOverlap) {
      cfgFuncParas.fPhiCutLow = new TF1("fPhiCutLow", cfgFuncParas.cfgTPCPhiCutLowCutFunction->c_str(), 0, 100);
      cfgFuncParas.fPhiCutHigh = new TF1("fPhiCutHigh", cfgFuncParas.cfgTPCPhiCutHighCutFunction->c_str(), 0, 100);
    }

    if (cfgTrackDensityCorrUse) {
      std::vector<double> pTEffBins = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0};
      hFindPtBin = new TH1D("hFindPtBin", "hFindPtBin", pTEffBins.size() - 1, &pTEffBins[0]);
      funcEff.resize(pTEffBins.size() - 1);
      // LHC24g3 Eff
      std::vector<double> f1p0 = cfgFuncParas.cfgTrackDensityP0;
      std::vector<double> f1p1 = cfgFuncParas.cfgTrackDensityP1;
      for (uint ifunc = 0; ifunc < pTEffBins.size() - 1; ifunc++) {
        funcEff[ifunc] = new TF1(Form("funcEff%i", ifunc), "[0]+[1]*x", 0, 3000);
        funcEff[ifunc]->SetParameters(f1p0[ifunc], f1p1[ifunc]);
      }
      funcV2 = new TF1("funcV2", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      funcV2->SetParameters(0.0186111, 0.00351907, -4.38264e-05, 1.35383e-07, -3.96266e-10);
      funcV3 = new TF1("funcV3", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      funcV3->SetParameters(0.0174056, 0.000703329, -1.45044e-05, 1.91991e-07, -1.62137e-09);
      funcV4 = new TF1("funcV4", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      funcV4->SetParameters(0.008845, 0.000259668, -3.24435e-06, 4.54837e-08, -6.01825e-10);
    }
  }

  void createOutputObjectsForRun(int runNumber)
  {
    const AxisSpec axisPhi{60, 0.0, constants::math::TwoPI, "#varphi"};
    std::shared_ptr<TH3> histPhiEtaVtxz = registry.add<TH3>(Form("%d/hPhiEtaVtxz", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, {40, -10, 10}}});
    th3sPerRun.insert(std::make_pair(runNumber, histPhiEtaVtxz));
  }

  template <char... chars>
  void fillProfile(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& cent)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (std::fabs(val) < 1)
        registry.fill(tarName, cent, val, dnx);
      return;
    }
    return;
  }

  template <char... chars, char... chars2>
  void fillpTvnProfile(const GFW::CorrConfig& corrconf, const double& sum_pt, const double& WeffEvent, const ConstStr<chars...>& vnWeff, const ConstStr<chars2...>& vnpT, const double& cent)
  {
    double meanPt = sum_pt / WeffEvent;
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (std::fabs(val) < 1) {
        registry.fill(vnWeff, cent, val, dnx * WeffEvent);
        registry.fill(vnpT, cent, val * meanPt, dnx * WeffEvent);
      }
      return;
    }
    return;
  }

  template <char... chars>
  void fillDeltaPtvsCorr(const GFW::CorrConfig& corrconf, const double& sum_pt, const double& WeffEvent, const ConstStr<chars...>& histo, const double& cent)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (std::fabs(val) < 1) {
        double meanPt = sum_pt / WeffEvent;
        double deltaPt = meanPt - cfgFuncParas.hEvAvgMeanPt->GetBinContent(cfgFuncParas.hEvAvgMeanPt->FindBin(cent));
        registry.fill(histo, deltaPt / meanPt, val, cent, dnx * WeffEvent);
      }
      return;
    }
    return;
  }

  template <DataType dt>
  void fillFC(const GFW::CorrConfig& corrconf, const double& cent, const double& rndm)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (!corrconf.pTDif) {
      if (dnx == 0)
        return;
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (std::fabs(val) < 1) {
        (dt == kGen) ? fFCgen->FillProfile(corrconf.Head.c_str(), cent, val, dnx, rndm) : fFC->FillProfile(corrconf.Head.c_str(), cent, val, dnx, rndm);
      }
      return;
    }
    for (auto i = 1; i <= fPtAxis->GetNbins(); i++) {
      dnx = fGFW->Calculate(corrconf, i - 1, kTRUE).real();
      if (dnx == 0)
        continue;
      val = fGFW->Calculate(corrconf, i - 1, kFALSE).real() / dnx;
      if (std::fabs(val) < 1) {
        (dt == kGen) ? fFCgen->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i), cent, val, dnx, rndm) : fFC->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i), cent, val, dnx, rndm);
      }
    }
    return;
  }

  template <DataType dt, typename TTrack>
  inline void fillPtSums(TTrack track, float weff)
  {
    if (std::abs(track.eta()) < cfgEtaPtPt) {
      (dt == kGen) ? fFCptgen->fill(1., track.pt()) : fFCpt->fill(weff, track.pt());
    }
  }

  template <DataType dt>
  void fillPtContainers(const float& centmult, const double& rndm)
  {
    (dt == kGen) ? fFCptgen->calculateCorrelations() : fFCpt->calculateCorrelations();
    (dt == kGen) ? fFCptgen->fillPtProfiles(centmult, rndm) : fFCpt->fillPtProfiles(centmult, rndm);
    (dt == kGen) ? fFCptgen->fillCMProfiles(centmult, rndm) : fFCpt->fillCMProfiles(centmult, rndm);
    for (uint l_ind = 0; l_ind < corrconfigsPtVn.size(); ++l_ind) {
      if (!corrconfigsPtVn.at(l_ind).pTDif) {
        auto dnx = fGFW->Calculate(corrconfigsPtVn.at(l_ind), 0, kTRUE).real();
        if (dnx == 0)
          continue;
        auto val = fGFW->Calculate(corrconfigsPtVn.at(l_ind), 0, kFALSE).real() / dnx;
        if (std::abs(val) < 1) {
          (dt == kGen) ? fFCptgen->fillVnPtProfiles(l_ind, centmult, val, dnx, rndm, gfwConfigs.GetpTCorrMasks()[l_ind]) : fFCpt->fillVnPtProfiles(l_ind, centmult, val, dnx, rndm, gfwConfigs.GetpTCorrMasks()[l_ind]);
        }
        continue;
      }
    }
    return;
  }

  void loadCorrections(uint64_t timestamp)
  {
    if (correctionsLoaded)
      return;
    if (cfgAcceptance.value.empty() == false) {
      mAcceptance = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance, timestamp);
      if (mAcceptance)
        LOGF(info, "Loaded acceptance weights from %s (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance);
      else
        LOGF(warning, "Could not load acceptance weights from %s (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance);
    }
    if (cfgEfficiency.value.empty() == false) {
      mEfficiency = ccdb->getForTimeStamp<TH1D>(cfgEfficiency, timestamp);
      if (mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEfficiency.value.c_str(), (void*)mEfficiency);
    }
    if (cfgFuncParas.cfgDptDisEnable && cfgFuncParas.cfgDptDishEvAvgMeanPt.value.empty() == false) {
      cfgFuncParas.hEvAvgMeanPt = ccdb->getForTimeStamp<TH1D>(cfgFuncParas.cfgDptDishEvAvgMeanPt, timestamp);
      if (cfgFuncParas.hEvAvgMeanPt == nullptr) {
        LOGF(fatal, "Could not load mean pT histogram from %s", cfgFuncParas.cfgDptDishEvAvgMeanPt.value.c_str());
      }
      LOGF(info, "Loaded mean pT histogram from %s (%p)", cfgFuncParas.cfgDptDishEvAvgMeanPt.value.c_str(), (void*)cfgFuncParas.hEvAvgMeanPt);
    }
    correctionsLoaded = true;
  }

  bool setCurrentParticleWeights(float& weight_nue, float& weight_nua, float phi, float eta, float pt, float vtxz)
  {
    float eff = 1.;
    if (mEfficiency)
      eff = mEfficiency->GetBinContent(mEfficiency->FindBin(pt));
    else
      eff = 1.0;
    if (eff == 0)
      return false;
    weight_nue = 1. / eff;

    if (mAcceptance)
      weight_nua = mAcceptance->getNUA(phi, eta, vtxz);
    else
      weight_nua = 1;
    return true;
  }

  template <typename TCollision>
  bool eventSelected(TCollision collision, const int multTrk, const float centrality)
  {
    registry.fill(HIST("hEventCountSpecific"), 0.5);
    if (cfgEvSelkNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return 0;
    }
    if (cfgEvSelkNoSameBunchPileup)
      registry.fill(HIST("hEventCountSpecific"), 1.5);
    if (cfgEvSelkNoITSROFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return 0;
    }
    if (cfgEvSelkNoITSROFrameBorder)
      registry.fill(HIST("hEventCountSpecific"), 2.5);
    if (cfgEvSelkNoTimeFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return 0;
    }
    if (cfgEvSelkNoTimeFrameBorder)
      registry.fill(HIST("hEventCountSpecific"), 3.5);
    if (cfgEvSelkIsGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return 0;
    }
    if (cfgEvSelkIsGoodZvtxFT0vsPV)
      registry.fill(HIST("hEventCountSpecific"), 4.5);
    if (cfgEvSelkNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      // no collisions in specified time range
      return 0;
    }
    if (cfgEvSelkNoCollInTimeRangeStandard)
      registry.fill(HIST("hEventCountSpecific"), 5.5);
    if (cfgEvSelkIsGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      // from Jan 9 2025 AOT meeting
      // cut time intervals with dead ITS staves
      return 0;
    }
    if (cfgEvSelkIsGoodITSLayersAll)
      registry.fill(HIST("hEventCountSpecific"), 6.5);
    if (cfgEvSelkNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      // no other collisions in this Readout Frame with per-collision multiplicity above threshold
      return 0;
    }
    if (cfgEvSelkNoCollInRofStandard)
      registry.fill(HIST("hEventCountSpecific"), 7.5);
    if (cfgEvSelkNoHighMultCollInPrevRof && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      // veto an event if FT0C amplitude in previous ITS ROF is above threshold
      return 0;
    }
    if (cfgEvSelkNoHighMultCollInPrevRof)
      registry.fill(HIST("hEventCountSpecific"), 8.5);
    auto multNTracksPV = collision.multNTracksPV();
    auto occupancy = collision.trackOccupancyInTimeRange();
    if (cfgEvSelOccupancy && (occupancy < cfgCutOccupancyLow || occupancy > cfgCutOccupancyHigh))
      return 0;
    if (cfgEvSelOccupancy)
      registry.fill(HIST("hEventCountSpecific"), 9.5);

    if (cfgEvSelMultCorrelation) {
      if (cfgFuncParas.cfgMultPVT0CCutEnabled) {
        if (multNTracksPV < cfgFuncParas.fMultPVT0CCutLow->Eval(centrality))
          return 0;
        if (multNTracksPV > cfgFuncParas.fMultPVT0CCutHigh->Eval(centrality))
          return 0;
      }
      if (cfgFuncParas.cfgMultT0CCutEnabled) {
        if (multTrk < cfgFuncParas.fMultT0CCutLow->Eval(centrality))
          return 0;
        if (multTrk > cfgFuncParas.fMultT0CCutHigh->Eval(centrality))
          return 0;
      }
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
    }
    if (cfgEvSelMultCorrelation)
      registry.fill(HIST("hEventCountSpecific"), 10.5);

    // V0A T0A 5 sigma cut
    float sigma = 5.0;
    if (cfgEvSelV0AT0ACut && (std::fabs(collision.multFV0A() - cfgFuncParas.fT0AV0AMean->Eval(collision.multFT0A())) > sigma * cfgFuncParas.fT0AV0ASigma->Eval(collision.multFT0A())))
      return 0;
    if (cfgEvSelV0AT0ACut)
      registry.fill(HIST("hEventCountSpecific"), 11.5);

    return 1;
  }

  int getMagneticField(uint64_t timestamp)
  {
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(cfgFuncParas.cfgMagnetField, timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found in %s for timestamp %llu", cfgFuncParas.cfgMagnetField.value.c_str(), timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP from %s for timestamp %llu with magnetic field of %d kG", cfgFuncParas.cfgMagnetField.value.c_str(), timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  template <typename TTrack>
  bool trackSelected(TTrack track)
  {
    return ((track.tpcNClsFound() >= cfgCutTPCclu) && (track.tpcNClsCrossedRows() >= cfgCutTPCCrossedRows) && (track.itsNCls() >= cfgCutITSclu));
  }

  template <typename TTrack>
  bool rejectionTPCoverlap(TTrack track, const int field)
  {
    double phimodn = track.phi();
    if (field < 0) // for negative polarity field
      phimodn = o2::constants::math::TwoPI - phimodn;
    if (track.sign() < 0) // for negative charge
      phimodn = o2::constants::math::TwoPI - phimodn;
    if (phimodn < 0)
      LOGF(warning, "phi < 0: %g", phimodn);

    float middle = o2::constants::math::TwoPI / 18.0;
    phimodn += middle; // to center gap in the middle
    phimodn = fmod(phimodn, o2::constants::math::TwoPI / 9.0);
    registry.fill(HIST("pt_phi_bef"), track.pt(), phimodn);
    if (cfgFuncParas.cfgRejectionTPCsectorOverlap) {
      if (track.pt() >= cfgFuncParas.cfgTPCPhiCutPtMin && phimodn < cfgFuncParas.fPhiCutHigh->Eval(track.pt()) && phimodn > cfgFuncParas.fPhiCutLow->Eval(track.pt()))
        return false; // reject track
    }
    registry.fill(HIST("pt_phi_aft"), track.pt(), phimodn);
    return true;
  }

  void initHadronicRate(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    if (gHadronicRate.find(mRunNumber) == gHadronicRate.end()) {
      auto runDuration = ccdb->getRunDuration(mRunNumber);
      mSOR = runDuration.first;
      mMinSeconds = std::floor(mSOR * 1.e-3);                /// round tsSOR to the highest integer lower than tsSOR
      double maxSec = std::ceil(runDuration.second * 1.e-3); /// round tsEOR to the lowest integer higher than tsEOR
      const AxisSpec axisSeconds{static_cast<int>((maxSec - mMinSeconds) / 20.f), 0, maxSec - mMinSeconds, "Seconds since SOR"};
      gHadronicRate[mRunNumber] = registry.add<TH2>(Form("HadronicRate/%i", mRunNumber), ";Time since SOR (s);Hadronic rate (kHz)", kTH2D, {axisSeconds, {510, 0., 51.}}).get();
    }
    gCurrentHadronicRate = gHadronicRate[mRunNumber];
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

  void processData(FilteredCollisions::iterator const& collision, aod::BCsWithTimestamps const&, FilteredTracks const& tracks)
  {
    registry.fill(HIST("hEventCount"), 0.5);
    if (!cfgUseSmallMemory && tracks.size() >= 1) {
      registry.fill(HIST("BeforeSel8_globalTracks_centT0C"), collision.centFT0C(), tracks.size());
    }
    if (!collision.sel8())
      return;
    if (tracks.size() < 1)
      return;
    registry.fill(HIST("hEventCount"), 1.5);
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int currentRunNumber = bc.runNumber();
    for (const auto& ExcludedRun : cfgRunRemoveList.value) {
      if (currentRunNumber == ExcludedRun) {
        return;
      }
    }
    if (cfgOutputNUAWeightsRunbyRun && currentRunNumber != lastRunNumber) {
      lastRunNumber = currentRunNumber;
      if (std::find(runNumbers.begin(), runNumbers.end(), currentRunNumber) == runNumbers.end()) {
        // if run number is not in the preconfigured list, create new output histograms for this run
        createOutputObjectsForRun(currentRunNumber);
        runNumbers.push_back(currentRunNumber);
      }

      if (th3sPerRun.find(currentRunNumber) == th3sPerRun.end()) {
        LOGF(fatal, "RunNumber %d not found in th3sPerRun", currentRunNumber);
        return;
      }
    }

    registry.fill(HIST("hEventCount"), 2.5);
    if (!cfgUseSmallMemory) {
      registry.fill(HIST("BeforeCut_globalTracks_centT0C"), collision.centFT0C(), tracks.size());
      registry.fill(HIST("BeforeCut_PVTracks_centT0C"), collision.centFT0C(), collision.multNTracksPV());
      registry.fill(HIST("BeforeCut_globalTracks_PVTracks"), collision.multNTracksPV(), tracks.size());
      registry.fill(HIST("BeforeCut_globalTracks_multT0A"), collision.multFT0A(), tracks.size());
      registry.fill(HIST("BeforeCut_globalTracks_multV0A"), collision.multFV0A(), tracks.size());
      registry.fill(HIST("BeforeCut_multV0A_multT0A"), collision.multFT0A(), collision.multFV0A());
      registry.fill(HIST("BeforeCut_multT0C_centT0C"), collision.centFT0C(), collision.multFT0C());
    }
    float cent = getCentrality(collision);

    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), cent))
      return;
    registry.fill(HIST("hEventCount"), 3.5);
    float lRandom = fRndm->Rndm();
    float vtxz = collision.posZ();
    registry.fill(HIST("hVtxZ"), vtxz);
    registry.fill(HIST("hMult"), tracks.size());
    registry.fill(HIST("hCent"), cent);
    fGFW->Clear();
    fFCpt->clearVector();
    if (cfgGetInteractionRate) {
      initHadronicRate(bc);
      double hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), mRunNumber, "ZNC hadronic") * 1.e-3; //
      double seconds = bc.timestamp() * 1.e-3 - mMinSeconds;
      if (cfgUseInteractionRateCut && (hadronicRate < cfgCutMinIR || hadronicRate > cfgCutMaxIR)) // cut on hadronic rate
        return;
      gCurrentHadronicRate->Fill(seconds, hadronicRate);
    }
    loadCorrections(bc.timestamp());
    registry.fill(HIST("hEventCount"), 4.5);

    // fill event QA
    if (!cfgUseSmallMemory) {
      registry.fill(HIST("globalTracks_centT0C"), collision.centFT0C(), tracks.size());
      registry.fill(HIST("PVTracks_centT0C"), collision.centFT0C(), collision.multNTracksPV());
      registry.fill(HIST("globalTracks_PVTracks"), collision.multNTracksPV(), tracks.size());
      registry.fill(HIST("globalTracks_multT0A"), collision.multFT0A(), tracks.size());
      registry.fill(HIST("globalTracks_multV0A"), collision.multFV0A(), tracks.size());
      registry.fill(HIST("multV0A_multT0A"), collision.multFT0A(), collision.multFV0A());
      registry.fill(HIST("multT0C_centT0C"), collision.centFT0C(), collision.multFT0C());
      registry.fill(HIST("centFT0CVar_centFT0C"), collision.centFT0C(), collision.centFT0CVariant1());
      registry.fill(HIST("centFT0M_centFT0C"), collision.centFT0C(), collision.centFT0M());
      registry.fill(HIST("centFV0A_centFT0C"), collision.centFT0C(), collision.centFV0A());
    }

    // track weights
    float weff = 1, wacc = 1;
    double weffEvent = 0;
    double ptSum = 0., ptSum_Gap08 = 0.;
    double weffEventWithinGap08 = 0., weffEventSquareWithinGap08 = 0.;
    double sumPtsquareWsquareWithinGap08 = 0., sumPtWsquareWithinGap08 = 0.;
    double nTracksCorrected = 0;
    int magnetfield = 0;
    float independent = cent;
    if (cfgUseNch)
      independent = static_cast<float>(tracks.size());
    if (cfgFuncParas.cfgShowTPCsectorOverlap) {
      // magnet field dependence cut
      magnetfield = getMagneticField(bc.timestamp());
    }

    double psi2Est = 0, psi3Est = 0, psi4Est = 0;
    float wEPeff = 1;
    double v2 = 0, v3 = 0, v4 = 0;
    // be cautious, this only works for Pb-Pb
    // esimate the Event plane and vn for this event
    if (cfgTrackDensityCorrUse) {
      double q2x = 0, q2y = 0;
      double q3x = 0, q3y = 0;
      double q4x = 0, q4y = 0;
      for (const auto& track : tracks) {
        bool withinPtRef = (cfgCutPtRefMin < track.pt()) && (track.pt() < cfgCutPtRefMax); // within RF pT rang
        if (withinPtRef) {
          q2x += std::cos(2 * track.phi());
          q2y += std::sin(2 * track.phi());
          q3x += std::cos(3 * track.phi());
          q3y += std::sin(3 * track.phi());
          q4x += std::cos(4 * track.phi());
          q4y += std::sin(4 * track.phi());
        }
      }
      psi2Est = std::atan2(q2y, q2x) / 2.;
      psi3Est = std::atan2(q3y, q3x) / 3.;
      psi4Est = std::atan2(q4y, q4x) / 4.;
      v2 = funcV2->Eval(cent);
      v3 = funcV3->Eval(cent);
      v4 = funcV4->Eval(cent);
    }

    for (const auto& track : tracks) {
      if (!trackSelected(track))
        continue;
      if (cfgFuncParas.cfgShowTPCsectorOverlap && !rejectionTPCoverlap(track, magnetfield))
        continue;
      bool withinPtPOI = (cfgCutPtPOIMin < track.pt()) && (track.pt() < cfgCutPtPOIMax); // within POI pT range
      bool withinPtRef = (cfgCutPtRefMin < track.pt()) && (track.pt() < cfgCutPtRefMax); // within RF pT range
      bool withinEtaGap08 = (std::abs(track.eta()) < cfgEtaPtPt);
      if (cfgOutputNUAWeights) {
        if (cfgOutputNUAWeightsRefPt) {
          if (withinPtRef) {
            fWeights->fill(track.phi(), track.eta(), vtxz, track.pt(), cent, 0);
            if (cfgOutputNUAWeightsRunbyRun)
              th3sPerRun[currentRunNumber]->Fill(track.phi(), track.eta(), collision.posZ());
          }
        } else {
          fWeights->fill(track.phi(), track.eta(), vtxz, track.pt(), cent, 0);
          if (cfgOutputNUAWeightsRunbyRun)
            th3sPerRun[currentRunNumber]->Fill(track.phi(), track.eta(), collision.posZ());
        }
      }
      if (!setCurrentParticleWeights(weff, wacc, track.phi(), track.eta(), track.pt(), vtxz))
        continue;
      if (cfgTrackDensityCorrUse && withinPtRef) {
        double fphi = v2 * std::cos(2 * (track.phi() - psi2Est)) + v3 * std::cos(3 * (track.phi() - psi3Est)) + v4 * std::cos(4 * (track.phi() - psi4Est));
        fphi = (1 + 2 * fphi);
        int pTBinForEff = hFindPtBin->FindBin(track.pt());
        if (pTBinForEff >= 1 && pTBinForEff <= hFindPtBin->GetNbinsX()) {
          wEPeff = funcEff[pTBinForEff - 1]->Eval(fphi * tracks.size());
          if (wEPeff > 0.) {
            wEPeff = 1. / wEPeff;
            weff *= wEPeff;
          }
        }
      }
      registry.fill(HIST("hPt"), track.pt());
      if (withinPtRef) {
        registry.fill(HIST("hPhi"), track.phi());
        registry.fill(HIST("hPhiWeighted"), track.phi(), wacc);
        registry.fill(HIST("hEta"), track.eta());
        registry.fill(HIST("hPtRef"), track.pt());
        registry.fill(HIST("hChi2prTPCcls"), track.tpcChi2NCl());
        registry.fill(HIST("hChi2prITScls"), track.itsChi2NCl());
        registry.fill(HIST("hnTPCClu"), track.tpcNClsFound());
        registry.fill(HIST("hnITSClu"), track.itsNCls());
        registry.fill(HIST("hnTPCCrossedRow"), track.tpcNClsCrossedRows());
        registry.fill(HIST("hDCAz"), track.dcaZ(), track.pt());
        registry.fill(HIST("hDCAxy"), track.dcaXY(), track.pt());
        weffEvent += weff;
        ptSum += weff * track.pt();
        nTracksCorrected += weff;
        if (withinEtaGap08) {
          ptSum_Gap08 += weff * track.pt();
          sumPtWsquareWithinGap08 += weff * weff * track.pt();
          sumPtsquareWsquareWithinGap08 += weff * weff * track.pt() * track.pt();
          weffEventWithinGap08 += weff;
          weffEventSquareWithinGap08 += weff * weff;
        }
      }
      if (withinPtRef)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), wacc * weff, 1);
      if (withinPtPOI)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), wacc * weff, 2);
      if (withinPtPOI && withinPtRef)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), wacc * weff, 4);
      if (cfgUsePtRef && withinPtRef)
        fillPtSums<kReco>(track, weff);
      if (!cfgUsePtRef && withinPtPOI)
        fillPtSums<kReco>(track, weff);
    }
    registry.fill(HIST("hTrackCorrection2d"), tracks.size(), nTracksCorrected);

    if (cfgFuncParas.cfgDptDisEnable) {
      fillDeltaPtvsCorr(corrconfigs.at(cfgFuncParas.cfgDptDisCorrConfigIndex), ptSum, weffEvent, HIST("hNormDeltaPt_Corr_X"), independent);
    }

    double weffEventDiffWithGap08 = weffEventWithinGap08 * weffEventWithinGap08 - weffEventSquareWithinGap08;
    // MeanPt
    if (weffEvent) {
      registry.fill(HIST("hMeanPt"), independent, ptSum / weffEvent, weffEvent);
    }
    if (weffEventWithinGap08)
      registry.fill(HIST("hMeanPtWithinGap08"), independent, ptSum_Gap08 / weffEventWithinGap08, weffEventWithinGap08);
    // c22_gap8 * pt_withGap8
    if (weffEventWithinGap08)
      fillpTvnProfile(corrconfigs.at(7), ptSum_Gap08, weffEventWithinGap08, HIST("c22_gap08_Weff"), HIST("c22_gap08_trackMeanPt"), independent);
    // PtVariance
    if (weffEventDiffWithGap08) {
      registry.fill(HIST("PtVariance_partA_WithinGap08"), independent,
                    (ptSum_Gap08 * ptSum_Gap08 - sumPtsquareWsquareWithinGap08) / weffEventDiffWithGap08,
                    weffEventDiffWithGap08);
      registry.fill(HIST("PtVariance_partB_WithinGap08"), independent,
                    (weffEventWithinGap08 * ptSum_Gap08 - sumPtWsquareWithinGap08) / weffEventDiffWithGap08,
                    weffEventDiffWithGap08);
    }

    // Filling Flow Container
    for (uint l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
      fillFC<kReco>(corrconfigs.at(l_ind), independent, lRandom);
    }
    // Filling pt Container
    fillPtContainers<kReco>(independent, lRandom);
  }
  PROCESS_SWITCH(FlowTask, processData, "Process analysis for non-derived data", true);

  void processMCGen(FilteredMcCollisions::iterator const& mcCollision, FilteredSmallGroupMcCollisions const& collisions, FilteredMcParticles const& mcParticles)
  {
    if (collisions.size() != 1)
      return;

    float cent = -1.;
    for (const auto& collision : collisions) {
      cent = getCentrality(collision);
    }

    float lRandom = fRndm->Rndm();
    float vtxz = mcCollision.posZ();
    registry.fill(HIST("MCGen/MChVtxZ"), vtxz);
    registry.fill(HIST("MCGen/MChMult"), mcParticles.size());
    registry.fill(HIST("MCGen/MChCent"), cent);
    float independent = cent;
    if (cfgUseNch)
      independent = static_cast<float>(mcParticles.size());

    fGFW->Clear();
    fFCptgen->clearVector();

    for (const auto& mcParticle : mcParticles) {
      if (!mcParticle.isPhysicalPrimary())
        continue;
      bool withinPtPOI = (cfgCutPtPOIMin < mcParticle.pt()) && (mcParticle.pt() < cfgCutPtPOIMax); // within POI pT range
      bool withinPtRef = (cfgCutPtRefMin < mcParticle.pt()) && (mcParticle.pt() < cfgCutPtRefMax); // within RF pT range

      if (withinPtRef) {
        registry.fill(HIST("MCGen/MChPhi"), mcParticle.phi());
        registry.fill(HIST("MCGen/MChEta"), mcParticle.eta());
        registry.fill(HIST("MCGen/MChPtRef"), mcParticle.pt());
      }
      if (withinPtRef)
        fGFW->Fill(mcParticle.eta(), fPtAxis->FindBin(mcParticle.pt()) - 1, mcParticle.phi(), 1., 1);
      if (withinPtPOI)
        fGFW->Fill(mcParticle.eta(), fPtAxis->FindBin(mcParticle.pt()) - 1, mcParticle.phi(), 1., 2);
      if (withinPtPOI && withinPtRef)
        fGFW->Fill(mcParticle.eta(), fPtAxis->FindBin(mcParticle.pt()) - 1, mcParticle.phi(), 1., 4);
      if (cfgUsePtRef && withinPtRef)
        fillPtSums<kGen>(mcParticle, 1.);
      if (!cfgUsePtRef && withinPtPOI)
        fillPtSums<kGen>(mcParticle, 1.);
    }

    // Filling Flow Container
    for (uint l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
      fillFC<kGen>(corrconfigs.at(l_ind), independent, lRandom);
    }
    // Filling pt Container
    fillPtContainers<kGen>(independent, lRandom);
  }
  PROCESS_SWITCH(FlowTask, processMCGen, "Process analysis for MC generated events", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowTask>(cfgc)};
}
