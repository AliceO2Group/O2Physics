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

/// \file   flowSmallSystems.cxx
/// \author Vojtech Machacek (vojtech.machacek@cern.ch)
/// \since  Oct/05/2026
/// \brief  Task to measure flow observables with cumulants method in small collision systems.

#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWConfig.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TNamed.h>
#include <TRandom3.h>
#include <TString.h>

#include <sys/types.h>

#include <RtypesCore.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::genericframework;
using namespace o2::aod::rctsel;
#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP}; // NOLINT(bugprone-macro-parentheses)
static constexpr std::array<std::array<float, 2>, 7> LongArrayFloat = {{{{2.5, 2.5}}, {{50.0, 50.0}}, {{70.0, 70.0}}, {{5.0, 5.0}}, {{7.0, 7.0}}, {{7.0, 7.0}}, {{0., 0.}}}};
static constexpr std::array<std::array<int, 2>, 9> LongArrayInt = {{{{1, 2}}, {{1, 2}}, {{1, 2}}, {{1, 2}}, {{1, 2}}, {{1, 2}}, {{1, 2}}, {{1, 2}}, {{1, 2}}}};
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
struct FlowSmallSystems {
  struct : ConfigurableGroup{
             O2_DEFINE_CONFIGURABLE(cfgCentEstimator, int, 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FV0A; 4:NTPV; 5:NGlobal; 6:MFT")
               O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
                 O2_DEFINE_CONFIGURABLE(cfgEvSelRCTflags, std::string, "CBT", "keep empty to disable, usage: 'CBT', 'CBT_hadronPID' ")
                   O2_DEFINE_CONFIGURABLE(cfgEvSelkIsGoodITSLayersFlag, int, 0, "0: Disabled; 1: kIsGoodITSLayersAll; 2: kIsGoodITSLayer0123; 3: kIsGoodITSLayer3")
                     O2_DEFINE_CONFIGURABLE(cfgEvSelUsage, bool, true, "Use event selection cuts: (0/1)")
                       O2_DEFINE_CONFIGURABLE(cfgEvSel8, bool, true, "rejects collisions which are not selected by sel8")
                         O2_DEFINE_CONFIGURABLE(cfgEvSelkNoSameBunchPileup, bool, true, "rejects collisions which are associated with the same found-by-T0 bunch crossing")
                           O2_DEFINE_CONFIGURABLE(cfgEvSelkNoITSROFrameBorder, bool, true, "reject events at ITS ROF border")
                             O2_DEFINE_CONFIGURABLE(cfgEvSelkNoTimeFrameBorder, bool, true, "reject events at TF border")
                               O2_DEFINE_CONFIGURABLE(cfgEvSelkIsGoodZvtxFT0vsPV, bool, true, "removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference, use this cut at low multiplicities with caution")
                                 O2_DEFINE_CONFIGURABLE(cfgEvSelkNoCollInTimeRangeStandard, bool, true, "no collisions in specified time range")
                                   O2_DEFINE_CONFIGURABLE(cfgEvSelkNoCollInRofStandard, bool, true, "no other collisions in this Readout Frame with per-collision multiplicity above threshold")
                                     O2_DEFINE_CONFIGURABLE(cfgEvSelkNoHighMultCollInPrevRof, bool, true, "veto an event if FT0C amplitude in previous ITS ROF is above threshold")} cfgEventSelection;
  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgEvSelMultCorrelation, bool, true, "Multiplicity correlation cut")
    O2_DEFINE_CONFIGURABLE(cfgGlobalMultCentCutEnabled, bool, true, "Enable global multiplicity vs centrality cut")
    O2_DEFINE_CONFIGURABLE(cfgGlobalMultCentLowCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for global multiplicity vs centrality lower cut");
    O2_DEFINE_CONFIGURABLE(cfgGlobalMultCentHighCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 10.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for global multiplicity vs centrality upper cut");
    Configurable<std::vector<double>> cfgGlobalMultCentCutPars{"cfgGlobalMultCentCutPars", std::vector<double>{177.409, -5.06256, 0.0695392, -0.000539089, 0.00000177447, 27.3677, -0.427405, 0.00452115, -0.0000382143, 1.13518e-7}, "Global multiplicity vs centrality cut values"};
    O2_DEFINE_CONFIGURABLE(cfgGlobalMultPVMultCutEnabled, bool, true, "Enable global multiplicity vs PV multiplicity cut");
    O2_DEFINE_CONFIGURABLE(cfgGlobalMultPVMultLowCutFunction, std::string, "[0]+[1]*x - 5.*([2]+[3]*x)", "Functional for global multiplicity vs PV multiplicity correlation lower cut");
    O2_DEFINE_CONFIGURABLE(cfgGlobalMultPVMultHighCutFunction, std::string, "[0]+[1]*x + 5.*([2]+[3]*x)", "Functional for global multiplicity vs PV multiplicity correlation upper cut");
    Configurable<std::vector<double>> cfgGlobalMultPVMultCutPars{"cfgGlobalMultPVMultCutPars", std::vector<double>{-0.0893236, 0.721511, 2.05725, 0.0244317}, "global multiplicity vs PV multiplicity cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgGlobalMultV0AMultCutEnabled, bool, true, "Enable global multiplicity vs V0A multiplicity cut");
    O2_DEFINE_CONFIGURABLE(cfgGlobalMultV0AMultLowCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for global multiplicity vs V0A multiplicity correlation lower cut");
    O2_DEFINE_CONFIGURABLE(cfgGlobalMultV0AMultHighCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 4.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for global multiplicity vs V0A multiplicity correlation upper cut");
    Configurable<std::vector<double>> cfgGlobalMultV0AMultCutPars{"cfgGlobalMultV0AMultCutPars", std::vector<double>{485.385, 191.561, 0.296174, -0.00282594, 0.00000497112, 370.073, 102.096, -1.00597, 0.00405269, -0.00000520521}, "Global multiplicity vs V0A multiplicity cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgPVMultCentCutEnabled, bool, true, "Enable primary vertex multiplicity vs centrality cut");
    O2_DEFINE_CONFIGURABLE(cfgPVMultCentLowCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for primary vertex multiplicity vs centrality upper cut");
    O2_DEFINE_CONFIGURABLE(cfgPVMultCentHighCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 10.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for primary vertex multiplicity vs centrality lower cut");
    Configurable<std::vector<double>> cfgPVMultCentCutPars{"cfgPVMultCentCutPars", std::vector<double>{245.882, -6.90255, 0.0930817, -0.000710842, 0.00000231095, 35.9647, -0.530632, 0.00513954, -0.0000411738, 1.10106e-7}, "Primary vertex multiplicity vs centrality cut parameter values"};
  } cfgMultiplicityCorrelationCut;
  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgUseGFW, bool, true, "Enable analysis using GFW");
    O2_DEFINE_CONFIGURABLE(cfgUseMultWeights, bool, true, "Enable multiplicity weights in flow calculation");
    O2_DEFINE_CONFIGURABLE(cfgBootstrapSamples, int, 10, "Number of bootstrap subsamples");
    O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, 0, "Use Nch for cumulant analysis instead of centrality");
    O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
    O2_DEFINE_CONFIGURABLE(cfgEfficiencyNch, std::string, "", "CCDB path to efficiency object, only for Nch correction")
    O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")
    O2_DEFINE_CONFIGURABLE(cfgCutPtRefMax, float, 3.0f, "Maximal pT for poi tracks")
    O2_DEFINE_CONFIGURABLE(cfgCutPtRefMin, float, 0.2f, "Minimal pT for poi tracks")
    O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 100.0f, "Maximal pT for poi tracks")
    O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
    Configurable<GFWRegions> cfgRegions{"cfgRegions", {{"refN", "refP", "refFull"}, {-0.8, 0.4, -0.8}, {-0.4, 0.8, 0.8}, {0, 0, 0}, {1, 1, 1}}, "Extra configurations for GFW regions (Name, Eta max, Eta min, pT-diff, BitMask)"};
    Configurable<LabeledArray<int>> cfgDefineGFWRegions{"cfgDefineGFWRegions", {LongArrayInt.front().data(), 9, 1, {"Standard", "DeltaEta>0", "DeltaEta>0.2", "DeltaEta>0.4", "DeltaEta>0.6", "DeltaEta>0.8", "DeltaEta>1.0", "DeltaEta>1.2", "DeltaEta>1.4"}, {"Enabled"}}, "Switch (0/1) for various GFW region initialzation"};
    Configurable<GFWCorrConfigs> cfgCorrConfig{"cfgCorrConfig", {{"refP {2} refN {-2}", "refP {3} refN {-3}", "refP {4} refN {-4}", "refFull {2 -2}", "refFull {2 2 -2 -2}"}, {"ChGap22", "ChGap32", "ChGap42", "ChFull22", "ChFull24"}, {0, 0, 0, 0, 0}, {15, 1, 1, 0, 0}}, "Configurations for each correlation to calculate"};
  } cfgGFWSettings;
  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgUseHighMultTrigger, bool, false, "Enable HM trigger");
    O2_DEFINE_CONFIGURABLE(cfgHighMultiplicityTriggerMenu, int, 0, "0: HighTrackMult; 1: HighMultFv0; 2: HighFt0Mult; 3: HighFt0cFv0Mult; 4: HighFt0Flat; 5: HighFt0cFv0Flat");
  } cfgTriggeredData;
  O2_DEFINE_CONFIGURABLE(cfgNUAWeightsRunByRun, bool, true, "Fill NUA weights run by run")
  O2_DEFINE_CONFIGURABLE(cfgGetNumberOfCorr, bool, false, "Get histograms with number of correlations")
  Configurable<std::vector<int>> cfgFirstFillRunsList{"cfgFirstFillRunsList", std::vector<int>{}, "First run of the fill for time-dependent analysis, keep empty to disable"};
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 10.0f, "Maximal pT for all tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for all tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxy, std::string, "(0.0026+0.005/(x^1.01))", "function for pT dependent DCAxy cut, keep empty for const. cut based on DCAxy (Nsigma)")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, std::string, "(0.0026+0.005/(x^1.01))", "function for pT dependent DCAz cut, keep empty for const. cut based on DCAz (Nsigma)")
  Configurable<std::vector<int>> cfgRemovedRunList{"cfgRemovedRunList", std::vector<int>{-1}, "excluded run numbers"};
  Configurable<LabeledArray<float>> cfgTrackCuts{"cfgTrackCuts", {LongArrayFloat.front().data(), 7, 2, {"chi2 per TPCcls", "TPC cluster", "TPC crossed rows", "ITS cluster", "DCAxy (Nsigma)", "DCAz (Nsigma)", "DCAz (cm) override"}, {"Nch", "Observable"}}, "separate Nch and observable track selections"};
  ConfigurableAxis binningNch{"binningNch", {VARIABLE_WIDTH, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 460, 470, 480, 490, 500}, "binning for N_{ch} histograms"};
  ConfigurableAxis binningTime{"binningTime", {12, 0, 12}, "binning for time-dependent histograms QA"};
  ConfigurableAxis binningPt{"binningPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10}, "binning for pT histograms"};
  HistogramRegistry registryQA{"registryQA"};
  HistogramRegistry registryCorr{"registryCorr"};
  HistogramRegistry registryTrigger{"registryTrigger"};
  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgEventSelection.cfgCutVertex);
  Filter collisionFilterMC = (nabs(aod::mccollision::posZ) < cfgEventSelection.cfgCutVertex);
  Filter trackFilter = ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == static_cast<uint8_t>(true))) && (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax);
  Service<ccdb::BasicCCDBManager> ccdb{};
  Service<o2::framework::O2DatabasePDG> pdgDB{};
  RCTFlagsChecker rctChecker{"CBT"};
  enum TrackSelectionGroup {
    TrkSelNch = 0,
    TrkSelObs,
    TrkReco,
    TrkCounter
  };
  enum TrackCut {
    Chi2prTPCcls = 0, // max chi2 per TPC clusters
    TPCclu,           // minimum TPC found clusters
    TPCCrossedRows,   // minimum TPC crossed rows
    ITSclu,           // minimum ITS found clusters
    DCAxyNsigma,      // cut on number of sigma deviations from expected DCA in the transverse direction, nsigma=7 is the same with global track; set cfgCutDCAxy to "1" to change the cut to a fixed one
    DCAzNsigma,       // cut on number of sigma deviations from expected DCA in the transverse direction, nsigma=7 is the same with global track; set cfgCutDCAxy to "1" to change the cut to a fixed one
    DCAz,             // constant DCAz cut, set to 0 to override with pT dependent one
    Count_TrackCut
  };
  enum CentEstimators {
    CentFT0C = 0,
    CentFT0CVariant1,
    CentFT0M,
    CentFV0A,
    CentNTPV,
    CentNGlobal,
    CentMFT,
    // Count the total number of enum
    Count_CentEstimators
  };
  enum GoodITSLayersFlag {
    Disabled = 0,
    ITSLayersAll,
    ITSLayer0123,
    ITSLayer3,
    Count_ITSLayersFlag
  };
  enum EventSelectionCuts {
    FilteredEvent = 1,
    Sel8,
    RCTflags,
    IsGoodITSLayers,
    NoSameBunchPileup,
    NoITSROFrameBorder,
    NoTimeFrameBorder,
    IsGoodZvtxFT0vsPV,
    NoCollInTimeRangeStandard,
    NoCollInRofStandard,
    NoHighMultCollInPrevRof,
    MultCorrelation,
    BadRunsRemoved
  };
  enum EtaGapsSettings {
    NoGap = -1,
    Gap00,
    Gap02,
    Gap04,
    Gap06,
    Gap08,
    Gap10,
    Gap12,
    Gap14,
    EtaGapCounter
  };
  enum QAFillTime {
    BeforeSel,
    AfterSel
  };
  enum MinParticlesInCorrelation {
    Min2PC = 0,
    Min4PC = 1,
    Min6PC = 2,
    Min8PC = 3,

  };
  enum HighMultTriggerMenu {
    HighTrackMult = 0,
    HighMultFv0,
    HighFt0Mult,
    HighFt0cFv0Mult,
    HighFt0Flat,
    HighFt0cFv0Flat,
    TriggerMenuCounter
  };
  struct EventCut {
    bool enabled;
    int histBin;
    int flag; // just store the enum
  };
  struct TrackRegion {
    float forward;
    float backward;
    float gap;
  };
  struct Weights {
    TH1D* hEfficiencyNch = nullptr;
    TH1D* hEfficiency = nullptr;
    GFWWeights* hAcceptance = nullptr;
  } weight;
  static constexpr std::array<std::string_view, 2> FillTimeName = {"beforeEventSelQA/", "afterEventSelQA/"};
  static constexpr std::array<std::string_view, 3> FillTrackName = {"cutNch/", "cutObs/", "Reco/"};
  std::array<std::string, TriggerMenuCounter> triggerName = {"fHighTrackMult", "HighMultFv0", "fHighFt0Mult", "fHighFt0cFv0Mult", "fHighFt0Flat", "fHighFt0cFv0Flat"};

  std::array<EventCut, 8> eventSelectionCuts{};
  std::array<std::array<float, 2>, 7> trackSel{};
  std::array<std::array<int, 1>, 9> gfwRegions{};

  TF1* fGlobalMultPVMultLowCut = nullptr;
  TF1* fGlobalMultPVMultHighCut = nullptr;
  TF1* fGlobalMultV0AMultLowCut = nullptr;
  TF1* fGlobalMultV0AMultHighCut = nullptr;
  TF1* fGlobalMultCentLowCut = nullptr;
  TF1* fGlobalMultCentHighCut = nullptr;
  TF1* fPVMultCentLowCut = nullptr;
  TF1* fPVMultCentHighCut = nullptr;

  std::array<TF1*, 2> fDCAzCut{};
  std::array<TF1*, 2> fDCAxyCut{};
  std::vector<int> runNumbers;
  std::map<int, std::shared_ptr<TH3>> th3sPerRun;
  std::map<int, std::shared_ptr<TH2>> th2sPerRun;

  int firstRunOfFill = -1;

  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrConfigs;
  GFWRegions cacheRegions;
  GFWCorrConfigs cacheCorrConfig;
  TAxis* axisPtGFW = nullptr;
  OutputObj<FlowContainer> fFC{"FlowContainer"};

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  TH1D* hEff = nullptr;
  TH1D* mEfficiencyForNch = nullptr;
  GFWWeights* mAcceptance = nullptr;
  TRandom3* fRndm = new TRandom3(0);

  int getCorrectITSLayersFlag(const int& flag)
  {
    switch (flag) {
      case ITSLayersAll:
        return o2::aod::evsel::kIsGoodITSLayersAll;
      case ITSLayer0123:
        return o2::aod::evsel::kIsGoodITSLayer0123;
      case ITSLayer3:
        return o2::aod::evsel::kIsGoodITSLayer3;
      default:
        return 0; // Disabled
    }
  }
  inline void setGFWRegionPredefined(int flagRegion)
  {
    switch (flagRegion) {
      case NoGap:
        fGFW->AddRegion("full", -0.8, 0.8, 1, 1);
        LOGF(info, "cfgRegion added: full, Eta Min: -0.8, Eta Max: 0.8, pT diff: 0, bitmask: 1");
        break;
      case Gap00:
        fGFW->AddRegion("refN00", -0.8, 0., 1, 1);
        fGFW->AddRegion("refP00", 0., 0.8, 1, 1);
        LOGF(info, "cfgRegion added: refN00, Eta Min: -0.8, Eta Max: 0., pT diff: 0, bitmask: 1");
        LOGF(info, "cfgRegion added: refP00, Eta Min: 0., Eta Max: 0.8, pT diff: 0, bitmask: 1");
        break;
      case Gap02:
        fGFW->AddRegion("refN02", -0.8, -0.1, 1, 1);
        fGFW->AddRegion("refP02", 0.1, 0.8, 1, 1);
        LOGF(info, "cfgRegion added: refN02, Eta Min: -0.8, Eta Max: -0.1, pT diff: 0, bitmask: 1");
        LOGF(info, "cfgRegion added: refP02, Eta Min: 0.1, Eta Max: 0.8, pT diff: 0, bitmask: 1");
        break;
      case Gap04:
        fGFW->AddRegion("refN04", -0.8, -0.2, 1, 1);
        fGFW->AddRegion("refP04", 0.2, 0.8, 1, 1);
        LOGF(info, "cfgRegion added: refN04, Eta Min: -0.8, Eta Max: -0.2, pT diff: 0, bitmask: 1");
        LOGF(info, "cfgRegion added: refP04, Eta Min: 0.2, Eta Max: 0.8, pT diff: 0, bitmask: 1");
        break;
      case Gap06:
        fGFW->AddRegion("refN06", -0.8, -0.3, 1, 1);
        fGFW->AddRegion("refP06", 0.3, 0.8, 1, 1);
        LOGF(info, "cfgRegion added: refN06, Eta Min: -0.8, Eta Max: -0.3, pT diff: 0, bitmask: 1");
        LOGF(info, "cfgRegion added: refP06, Eta Min: 0.3, Eta Max: 0.8, pT diff: 0, bitmask: 1");
        break;
      case Gap08:
        fGFW->AddRegion("refN08", -0.8, -0.4, 1, 1);
        fGFW->AddRegion("refP08", 0.4, 0.8, 1, 1);
        LOGF(info, "cfgRegion added: refN08, Eta Min: -0.8, Eta Max: -0.4, pT diff: 0, bitmask: 1");
        LOGF(info, "cfgRegion added: refP08, Eta Min: 0.4, Eta Max: 0.8, pT diff: 0, bitmask: 1");
        break;
      case Gap10:
        fGFW->AddRegion("refN10", -0.8, -0.5, 1, 1);
        fGFW->AddRegion("refP10", 0.5, 0.8, 1, 1);
        LOGF(info, "cfgRegion added: refN10, Eta Min: -0.8, Eta Max: -0.5, pT diff: 0, bitmask: 1");
        LOGF(info, "cfgRegion added: refP10, Eta Min: 0.5, Eta Max: 0.8, pT diff: 0, bitmask: 1");
        break;
      case Gap12:
        fGFW->AddRegion("refN12", -0.8, -0.6, 1, 1);
        fGFW->AddRegion("refP12", 0.6, 0.8, 1, 1);
        LOGF(info, "cfgRegion added: refN12, Eta Min: -0.8, Eta Max: -0.6, pT diff: 0, bitmask: 1");
        LOGF(info, "cfgRegion added: refP12, Eta Min: 0.6, Eta Max: 0.8, pT diff: 0, bitmask: 1");
        break;
      case Gap14:
        fGFW->AddRegion("refN14", -0.8, -0.7, 1, 1);
        fGFW->AddRegion("refP14", 0.7, 0.8, 1, 1);
        LOGF(info, "cfgRegion added: refN14, Eta Min: -0.8, Eta Max: -0.7, pT diff: 0, bitmask: 1");
        LOGF(info, "cfgRegion added: refP14, Eta Min: 0.7, Eta Max: 0.8, pT diff: 0, bitmask: 1");
        break;
      default:
        LOGF(warning, Form("Used flag: %i is not defined! No region has been set.", flagRegion));
    }
  }
  void init(InitContext const&)
  {
    double binningNchMax = binningNch.value.back();
    std::vector<double> binningCent = std::vector<double>(101);
    iota(binningCent.begin(), binningCent.end(), 0);
    AxisSpec axisVertex = {40, -20, 20, "Vtxz (cm)"};
    AxisSpec axisPhi = {60, 0.0, constants::math::TwoPI, "#varphi"};
    AxisSpec axisEta = {40, -1., 1., "#eta"};
    AxisSpec axisCent = {binningCent, "Centrality"};
    AxisSpec axisT0C = {1000, 0, 50000, "N_{ch} (T0C)"};
    AxisSpec axisT0A = {1800, 0, 180000, "N_{ch} (T0A)"};
    AxisSpec axisV0A = {1800, 0, 180000, "N_{ch} (V0A)"};
    AxisSpec axisNch = {static_cast<int>(binningNchMax), 0., binningNchMax, "N_{ch}"};
    AxisSpec axisPt = {binningPt, "p_{T} (GeV/c)"};
    AxisSpec axisGFW = {(cfgGFWSettings.cfgUseNch) ? binningNch : binningCent};
    int nPtBins = axisPt.binEdges.size() - 1;
    axisPtGFW = new TAxis(nPtBins, binningPt.value.data());
    if (doprocessData) {
      registryQA.add("cutNch/hTrackPhi", "Track #varphi distribution;#varphi;Count", {HistType::kTH1D, {axisPhi}});
      registryQA.add("cutNch/hTrackEta", "Track #eta distribution;#eta;Count", {HistType::kTH1D, {axisEta}});
      registryQA.add("cutNch/hPt", "p_{T} distribution for reconstructed particles in a given pT range;p_{T};Counts", {HistType::kTH1D, {axisPt}});
      registryQA.add("cutNch/hChi2prTPCcls", "#chi^{2}/cluster for the TPC track segment;#chi^{2}/cluster;Counts", {HistType::kTH1D, {{100, 0., 5.}}});
      registryQA.add("cutNch/hChi2prITScls", "#chi^{2}/cluster for the ITS track;#chi^{2}/cluster;Counts", {HistType::kTH1D, {{100, 0., 50.}}});
      registryQA.add("cutNch/hnTPCClu", "Number of found TPC clusters", {HistType::kTH1D, {{100, 40, 180}}});
      registryQA.add("cutNch/hnITSClu", "Number of found ITS clusters", {HistType::kTH1D, {{100, 0, 20}}});
      registryQA.add("cutNch/hnTPCCrossedRow", "Number of crossed TPC Rows", {HistType::kTH1D, {{100, 40, 180}}});
      registryQA.add("cutNch/hDCAz", "DCAz after cuts; DCAz (cm); p_{T}", {HistType::kTH2D, {{200, -1., 1.}, {200, 0, 5}}});
      registryQA.add("cutNch/hDCAxy", "DCAxy after cuts; DCAxy (cm); p_{T}", {HistType::kTH2D, {{200, -1., 1.}, {200, 0, 5}}});
      registryQA.add("cutNch/hPhiEtaVtxZ", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, {32, -0.8, 0.8}, {40, -10, 10}}});
      registryQA.addClone("cutNch/", "cutObs/");
      registryQA.add("hVtxZ", "Vertex Z distribution;Vtxz (cm);Count", {HistType::kTH1D, {axisVertex}});
      registryQA.add("hCentrality", "Centrality distribution;Centrality;Count", {HistType::kTH1D, {axisCent}});
      registryQA.add("hNch", "N_{ch} distribution;N_{ch}^{raw};Count", {HistType::kTH1D, {axisNch}});
      if (!cfgGFWSettings.cfgEfficiencyNch.value.empty()) {
        registryQA.add("hNchCorrected", "N_{ch} distribution;N_{ch}^{corrected};Count", {HistType::kTH1D, {axisNch}});
      }
    } else if (doprocessMC) {
      registryQA.add("Reco/hTrackPhi", "Track #varphi distribution;#varphi;Count", {HistType::kTH1D, {axisPhi}});
      registryQA.add("Reco/hTrackEta", "Track #eta distribution;#eta;Count", {HistType::kTH1D, {axisEta}});
      registryQA.add("Reco/hPt", "p_{T} distribution for reconstructed particles in a given pT range;p_{T};Counts", {HistType::kTH1D, {axisPt}});
      registryQA.add("Reco/hPtNch", "p_{T} distribution correlated with Nch distribution for reconstructed particles in a given pT range;p_{T};N_{ch};Counts", {HistType::kTH2D, {axisPt, axisNch}});
      registryQA.add("Reco/hPtCent", "p_{T} distribution correlated with Nch distribution for reconstructed particles in a given pT range;p_{T};Centrality (%);Counts", {HistType::kTH2D, {axisPt, axisCent}});
      registryQA.add("Reco/hChi2prTPCcls", "#chi^{2}/cluster for the TPC track segment;#chi^{2}/cluster;Counts", {HistType::kTH1D, {{100, 0., 5.}}});
      registryQA.add("Reco/hChi2prITScls", "#chi^{2}/cluster for the ITS track;#chi^{2}/cluster;Counts", {HistType::kTH1D, {{100, 0., 50.}}});
      registryQA.add("Reco/hnTPCClu", "Number of found TPC clusters", {HistType::kTH1D, {{100, 40, 180}}});
      registryQA.add("Reco/hnITSClu", "Number of found ITS clusters", {HistType::kTH1D, {{100, 0, 20}}});
      registryQA.add("Reco/hnTPCCrossedRow", "Number of crossed TPC Rows", {HistType::kTH1D, {{100, 40, 180}}});
      registryQA.add("Reco/hDCAz", "DCAz after cuts; DCAz (cm); p_{T}", {HistType::kTH2D, {{200, -1., 1.}, {200, 0, 5}}});
      registryQA.add("Reco/hDCAxy", "DCAxy after cuts; DCAxy (cm); p_{T}", {HistType::kTH2D, {{200, -1., 1.}, {200, 0, 5}}});
      registryQA.add("Reco/hPhiEtaVtxZ", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, {32, -0.8, 0.8}, {40, -10, 10}}});
      registryQA.add("Gen/hTrackPhi", "Track #varphi distribution;#varphi;Count", {HistType::kTH1D, {axisPhi}});
      registryQA.add("Gen/hTrackEta", "Track #eta distribution;#eta;Count", {HistType::kTH1D, {axisEta}});
      registryQA.add("Gen/hPt", "p_{T} distribution for generated particles in a given pT range;p_{T};Counts", {HistType::kTH1D, {axisPt}});
      registryQA.add("Gen/hPtNch", "p_{T} distribution correlated with Nch distribution for generated particles in a given pT range;p_{T};N_{ch};Counts", {HistType::kTH2D, {axisPt, axisNch}});
      registryQA.add("Gen/hPtCent", "p_{T} distribution correlated with Nch distribution for reconstructed particles in a given pT range;p_{T};Centrality (%);Counts", {HistType::kTH2D, {axisPt, axisCent}});
      registryQA.add("hVtxZ", "Vertex Z distribution;Vtxz (cm);Count", {HistType::kTH1D, {axisVertex}});
      registryQA.add("hCentrality", "Centrality distribution;Centrality;Count", {HistType::kTH1D, {axisCent}});
      registryQA.add("hNch", "N_{ch} distribution;N_{ch}^{raw};Count", {HistType::kTH1D, {axisNch}});
      registryQA.add("hGenVsReco", "Correlation between N_{ch} Gen and N_{ch} Reco; N_{ch} Gen; N_{ch} Reco", {HistType::kTH2D, {axisNch, axisNch}});
      registryQA.add("hEventTrueVsGen", "Number of generated and reconstructed events used", {HistType::kTH1D, {{2, 0., 2.}}});
      registryQA.get<TH1>(HIST("hEventTrueVsGen"))->GetXaxis()->SetBinLabel(1, "Gen Events");
      registryQA.get<TH1>(HIST("hEventTrueVsGen"))->GetXaxis()->SetBinLabel(2, "Reco Events");
    }
    registryQA.add("hEventCount", "Number of Event;; Count", {HistType::kTH1D, {{13, 0.5, 13.5}}});
    registryQA.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(FilteredEvent, (cfgTriggeredData.cfgUseHighMultTrigger) ? "Filtered event + HM Trigger" : "Filtered event");
    registryQA.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(Sel8, "after sel8");
    registryQA.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(RCTflags, "RCTflags");
    std::array<std::string, 4> itsLayersFlag = {"ITSLayers-Disabled", "ITSLayersAll", "ITSLayers0123", "ITSLayer3"};
    registryQA.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(IsGoodITSLayers, Form("kIsGood%s", itsLayersFlag[cfgEventSelection.cfgEvSelkIsGoodITSLayersFlag].c_str()));
    registryQA.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(NoSameBunchPileup, "kNoSameBunchPileup");
    registryQA.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(NoITSROFrameBorder, "kNoITSROFrameBorder");
    registryQA.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(NoTimeFrameBorder, "kNoTimeFrameBorder");
    registryQA.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(IsGoodZvtxFT0vsPV, "kIsGoodZvtxFT0vsPV");
    registryQA.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(NoCollInTimeRangeStandard, "kNoCollInTimeRangeStandard");
    registryQA.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(NoCollInRofStandard, "kNoCollInRofStandard");
    registryQA.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(NoHighMultCollInPrevRof, "kNoHighMultCollInPrevRof");
    registryQA.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(MultCorrelation, "MultCorrelation");
    registryQA.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(BadRunsRemoved, "BadRunsRemoved");
    registryQA.add("beforeEventSelQA/globalTracks_centT0C", "before event selection;FT0C centrality (%);mulplicity global tracks", {HistType::kTH2D, {axisCent, axisNch}});
    registryQA.add("beforeEventSelQA/PVTracks_centT0C", "before event selection;FT0C centrality (%);mulplicity PV tracks", {HistType::kTH2D, {axisCent, axisNch}});
    registryQA.add("beforeEventSelQA/globalTracks_PVTracks", "before event selection;mulplicity PV tracks;mulplicity global tracks", {HistType::kTH2D, {axisNch, axisNch}});
    registryQA.add("beforeEventSelQA/globalTracks_multT0A", "before event selection;mulplicity T0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisNch}});
    registryQA.add("beforeEventSelQA/globalTracks_multV0A", "before event selection;mulplicity V0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisNch}});
    registryQA.add("beforeEventSelQA/multV0A_multT0A", "before event selection;mulplicity T0A;mulplicity V0A", {HistType::kTH2D, {axisT0A, axisT0A}});
    registryQA.add("beforeEventSelQA/multT0C_centT0C", "before event selection;FT0C centrality (%);mulplicity T0C", {HistType::kTH2D, {axisCent, axisT0C}});
    if (!cfgGFWSettings.cfgAcceptance.value.empty()) {
      registryQA.add("cutObs/hPhiCorrected", ";#varphi;Count", {HistType::kTH1D, {axisPhi}});
    }
    if (!cfgGFWSettings.cfgEfficiencyNch.value.empty()) {
      registryQA.add("hMultCorrected_MultUncorrected", "Multiplicity correlation matrix; N_ch^corected; N_ch^uncorected", {HistType::kTH2D, {axisNch, axisNch}});
    }
    if (cfgEventSelection.cfgEvSelUsage) {
      registryQA.addClone("beforeEventSelQA/", "afterEventSelQA/");
    }
    eventSelectionCuts = {
      EventCut({.enabled = (cfgEventSelection.cfgEvSelkIsGoodITSLayersFlag != Disabled), .histBin = IsGoodITSLayers, .flag = getCorrectITSLayersFlag(cfgEventSelection.cfgEvSelkIsGoodITSLayersFlag)}),
      EventCut({.enabled = cfgEventSelection.cfgEvSelkNoSameBunchPileup, .histBin = NoSameBunchPileup, .flag = o2::aod::evsel::kNoSameBunchPileup}),
      EventCut({.enabled = cfgEventSelection.cfgEvSelkNoITSROFrameBorder, .histBin = NoITSROFrameBorder, .flag = o2::aod::evsel::kNoITSROFrameBorder}),
      EventCut({.enabled = cfgEventSelection.cfgEvSelkNoTimeFrameBorder, .histBin = NoTimeFrameBorder, .flag = o2::aod::evsel::kNoTimeFrameBorder}),
      EventCut({.enabled = cfgEventSelection.cfgEvSelkIsGoodZvtxFT0vsPV, .histBin = IsGoodZvtxFT0vsPV, .flag = o2::aod::evsel::kIsGoodZvtxFT0vsPV}),
      EventCut({.enabled = cfgEventSelection.cfgEvSelkNoCollInTimeRangeStandard, .histBin = NoCollInTimeRangeStandard, .flag = o2::aod::evsel::kNoCollInTimeRangeStandard}),
      EventCut({.enabled = cfgEventSelection.cfgEvSelkNoCollInRofStandard, .histBin = NoCollInRofStandard, .flag = o2::aod::evsel::kNoCollInRofStandard}),
      EventCut({.enabled = cfgEventSelection.cfgEvSelkNoHighMultCollInPrevRof, .histBin = NoHighMultCollInPrevRof, .flag = o2::aod::evsel::kNoHighMultCollInPrevRof})};
    if (cfgEventSelection.cfgEvSelUsage) {
      LOGF(info, "Event selection cuts are enabled");
      LOGF(info, "Sel8 is %senabled", (cfgEventSelection.cfgEvSel8) ? "" : "not ");
      if (!cfgEventSelection.cfgEvSelRCTflags.value.empty()) {
        rctChecker.init(cfgEventSelection.cfgEvSelRCTflags.value);
        LOGF(info, "RCT is %senabled", (!cfgEventSelection.cfgEvSelRCTflags.value.empty()) ? "" : "not ");
        LOGF(info, "RCT selection flag: %s", cfgEventSelection.cfgEvSelRCTflags.value.c_str());
      }
      for (const auto& cut : eventSelectionCuts) {
        LOGF(info, "Flag %d is %senabled", cut.histBin, (cut.enabled) ? "" : "not ");
      }
      if (cfgMultiplicityCorrelationCut.cfgEvSelMultCorrelation) {
        LOGF(info, "Multiplicity correlation cuts are enabled");
        if (cfgMultiplicityCorrelationCut.cfgGlobalMultPVMultCutEnabled) {
          LOGF(info, "Global multiplicity vs PV multiplicity cut is enabled: %f", axisNch.binEdges.back());
          fGlobalMultPVMultHighCut = new TF1("fGlobalMultPVMultHighCut", cfgMultiplicityCorrelationCut.cfgGlobalMultPVMultHighCutFunction->c_str(), 0, axisNch.binEdges.back());
          fGlobalMultPVMultHighCut->SetParameters(cfgMultiplicityCorrelationCut.cfgGlobalMultPVMultCutPars.value.data());
          fGlobalMultPVMultLowCut = new TF1("fGlobalMultPVMultLowCut", cfgMultiplicityCorrelationCut.cfgGlobalMultPVMultLowCutFunction->c_str(), 0, axisNch.binEdges.back());
          fGlobalMultPVMultLowCut->SetParameters(cfgMultiplicityCorrelationCut.cfgGlobalMultPVMultCutPars.value.data());
          LOGF(info, "Global multiplicity vs PV multiplicity Upper Cut: %s", cfgMultiplicityCorrelationCut.cfgGlobalMultPVMultHighCutFunction->c_str());
          LOGF(info, "Global multiplicity vs PV multiplicity Lower Cut: %s", cfgMultiplicityCorrelationCut.cfgGlobalMultPVMultLowCutFunction->c_str());
          LOGF(info, "Global multiplicity vs PV multiplicity cut parameters:");
          for (std::size_t i = 0; i < cfgMultiplicityCorrelationCut.cfgGlobalMultPVMultCutPars.value.size(); ++i) {
            LOGF(info, "fGlobalMultPVMultLowCut par %zu = %g", i, fGlobalMultPVMultLowCut->GetParameter(i));
            LOGF(info, "fGlobalMultPVMultHighCut par %zu = %g", i, fGlobalMultPVMultHighCut->GetParameter(i));
          }
        }
        if (cfgMultiplicityCorrelationCut.cfgGlobalMultV0AMultCutEnabled) {
          LOGF(info, Form("Global multiplicity vs V0A multiplicity cut is enabled: %f ", axisV0A.binEdges.back()));

          fGlobalMultV0AMultHighCut = new TF1("fGlobalMultV0AMultHighCut", cfgMultiplicityCorrelationCut.cfgGlobalMultV0AMultHighCutFunction->c_str(), 0, axisV0A.binEdges.back());
          fGlobalMultV0AMultHighCut->SetParameters(cfgMultiplicityCorrelationCut.cfgGlobalMultV0AMultCutPars.value.data());
          fGlobalMultV0AMultLowCut = new TF1("fGlobalMultV0AMultLowCut", cfgMultiplicityCorrelationCut.cfgGlobalMultV0AMultLowCutFunction->c_str(), 0, axisV0A.binEdges.back());
          fGlobalMultV0AMultLowCut->SetParameters(cfgMultiplicityCorrelationCut.cfgGlobalMultV0AMultCutPars.value.data());
          LOGF(info, "Global multiplicity vs V0A multiplicity Upper Cut: %s", cfgMultiplicityCorrelationCut.cfgGlobalMultV0AMultHighCutFunction->c_str());
          LOGF(info, "Global multiplicity vs V0A multiplicity Lower Cut: %s", cfgMultiplicityCorrelationCut.cfgGlobalMultV0AMultLowCutFunction->c_str());
          LOGF(info, "Global multiplicity vs V0A multiplicity cut parameters:");
          for (std::size_t i = 0; i < cfgMultiplicityCorrelationCut.cfgGlobalMultV0AMultCutPars.value.size(); ++i) {
            LOGF(info, "fGlobalMultV0AMultLowCut par %zu = %g", i, fGlobalMultV0AMultLowCut->GetParameter(i));
            LOGF(info, "fGlobalMultV0AMultHighCut par %zu = %g", i, fGlobalMultV0AMultHighCut->GetParameter(i));
          }
        }
        if (cfgMultiplicityCorrelationCut.cfgGlobalMultCentCutEnabled) {
          LOGF(info, "Global multiplicity vs centrality cut is enabled");
          fGlobalMultCentHighCut = new TF1("fGlobalMultCentHighCut", cfgMultiplicityCorrelationCut.cfgGlobalMultCentHighCutFunction->c_str(), 0, 100);
          fGlobalMultCentHighCut->SetParameters(cfgMultiplicityCorrelationCut.cfgGlobalMultCentCutPars.value.data());
          fGlobalMultCentLowCut = new TF1("fGlobalMultCentLowCut", cfgMultiplicityCorrelationCut.cfgGlobalMultCentLowCutFunction->c_str(), 0, 100);
          fGlobalMultCentLowCut->SetParameters(cfgMultiplicityCorrelationCut.cfgGlobalMultCentCutPars.value.data());
          LOGF(info, "Global multiplicity vs centrality Upper Cut: %s", cfgMultiplicityCorrelationCut.cfgGlobalMultCentHighCutFunction->c_str());
          LOGF(info, "Global multiplicity vs centrality Lower Cut: %s", cfgMultiplicityCorrelationCut.cfgGlobalMultCentLowCutFunction->c_str());
          LOGF(info, "Global multiplicity vs centrality cut parameters:");
          for (std::size_t i = 0; i < cfgMultiplicityCorrelationCut.cfgGlobalMultCentCutPars.value.size(); ++i) {
            LOGF(info, "fGlobalMultCentLowCut par %zu = %g", i, fGlobalMultCentLowCut->GetParameter(i));
            LOGF(info, "fGlobalMultCentHighCut par %zu = %g", i, fGlobalMultCentHighCut->GetParameter(i));
          }
        }
        if (cfgMultiplicityCorrelationCut.cfgPVMultCentCutEnabled) {
          LOGF(info, "PV multiplicity vs centrality cut is enabled");
          fPVMultCentHighCut = new TF1("fPVMultCentHighCut", cfgMultiplicityCorrelationCut.cfgPVMultCentHighCutFunction->c_str(), 0, 100);
          fPVMultCentHighCut->SetParameters(cfgMultiplicityCorrelationCut.cfgPVMultCentCutPars.value.data());
          fPVMultCentLowCut = new TF1("fPVMultCentLowCut", cfgMultiplicityCorrelationCut.cfgPVMultCentLowCutFunction->c_str(), 0, 100);
          fPVMultCentLowCut->SetParameters(cfgMultiplicityCorrelationCut.cfgPVMultCentCutPars.value.data());
          LOGF(info, "PV multiplicity vs centrality Upper Cut: %s", cfgMultiplicityCorrelationCut.cfgPVMultCentHighCutFunction->c_str());
          LOGF(info, "PV multiplicity vs centrality Lower Cut: %s", cfgMultiplicityCorrelationCut.cfgPVMultCentLowCutFunction->c_str());
          LOGF(info, "PV multiplicity vs centrality cut parameters:");
          for (std::size_t i = 0; i < cfgMultiplicityCorrelationCut.cfgPVMultCentCutPars.value.size(); ++i) {
            LOGF(info, "fPVMultCentLowCut par %zu = %g", i, fPVMultCentLowCut->GetParameter(i));
            LOGF(info, "fPVMultCentHighCut par %zu = %g", i, fPVMultCentHighCut->GetParameter(i));
          }
        }
      }
    } else {
      LOGF(warning, "Event selection cuts are disabled");
    }
    readMatrix(cfgTrackCuts->getData(), trackSel);
    if (trackSel[DCAz][TrkSelNch] == 0) {
      fDCAzCut[TrkSelNch] = new TF1("pT-DCAzNch", Form("[0]*%s", cfgCutDCAz->c_str()), 0.001, 100);
      fDCAzCut[TrkSelNch]->SetParameter(0, trackSel[DCAzNsigma][TrkSelNch]);
      LOGF(info, Form("pT-dependent DCAz cut has been selcted: [0]*%s", cfgCutDCAz->c_str()));
      LOGF(info, Form("pT-dependent DCAz cut has been selcted: %f*%s", fDCAzCut[TrkSelNch]->GetParameter(0), cfgCutDCAz->c_str()));
    }
    if (trackSel[DCAz][TrkSelObs] == 0) {
      fDCAzCut[TrkSelObs] = new TF1("pT-DCAzObs", Form("[0]*%s", cfgCutDCAz->c_str()), 0.001, 100);
      fDCAzCut[TrkSelObs]->SetParameter(0, trackSel[DCAzNsigma][TrkSelObs]);
    }
    fDCAxyCut[TrkSelNch] = new TF1("pT-DCAxyNch", Form("[0]*%s", cfgCutDCAxy->c_str()), 0.001, 100);
    fDCAxyCut[TrkSelNch]->SetParameter(0, trackSel[DCAxyNsigma][TrkSelNch]);
    fDCAxyCut[TrkSelObs] = new TF1("pT-DCAxyObs", Form("[0]*%s", cfgCutDCAxy->c_str()), 0.001, 100);
    fDCAxyCut[TrkSelObs]->SetParameter(0, trackSel[DCAxyNsigma][TrkSelObs]);
    LOGF(info, Form("pT-dependent DCAxy cut for Nch : %f*%s", fDCAxyCut[TrkSelObs]->GetParameter(0), cfgCutDCAxy->c_str()));
    std::array<std::string, 8> etaGapsFlag = {"2-sub |#Delta#eta|>0", "2-sub |#Delta#eta|>0.2", "2-sub |#Delta#eta|>0.4", "2-sub |#Delta#eta|>0.6", "2-sub |#Delta#eta|>0.8", "2-sub |#Delta#eta|>1.0", "2-sub |#Delta#eta|>1.2", "2-sub |#Delta#eta|>1.4"};
    if (cfgGetNumberOfCorr) {
      registryCorr.add("hTwoParticleCorrelations", "2 particle correlations; Nch;", {HistType::kTH2F, {axisNch, {8, -0.5, 7.5}}});
      registryCorr.add("hFourParticleCorrelations", "4 particle correlations; Nch;", {HistType::kTH2F, {axisNch, {8, -0.5, 7.5}}});
      registryCorr.add("hSixParticleCorrelations", "6 particle correlations; Nch;", {HistType::kTH2F, {axisNch, {8, -0.5, 7.5}}});
      registryCorr.add("hEightParticleCorrelations", "8 particle correlations; Nch;", {HistType::kTH2F, {axisNch, {8, -0.5, 7.5}}});
      for (uint iSet = 0; iSet < static_cast<uint>(EtaGapCounter); iSet++) {
        registryCorr.get<TH2>(HIST("hTwoParticleCorrelations"))->GetYaxis()->SetBinLabel(iSet + 1, etaGapsFlag[iSet].c_str());
        registryCorr.get<TH2>(HIST("hFourParticleCorrelations"))->GetYaxis()->SetBinLabel(iSet + 1, etaGapsFlag[iSet].c_str());
        registryCorr.get<TH2>(HIST("hSixParticleCorrelations"))->GetYaxis()->SetBinLabel(iSet + 1, etaGapsFlag[iSet].c_str());
        registryCorr.get<TH2>(HIST("hEightParticleCorrelations"))->GetYaxis()->SetBinLabel(iSet + 1, etaGapsFlag[iSet].c_str());
      }
    }
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
    if (!cfgFirstFillRunsList.value.empty()) {
      registryQA.add("multT0A_multT0C", "FT0A-C multiplicity correlation;N_{ch}(FT0A);N_{ch}(FT0C);t[h]", {HistType::kTH3D, {axisT0A, axisT0C, binningTime}});
      firstRunOfFill = *std::min_element(cfgFirstFillRunsList.value.begin(), cfgFirstFillRunsList.value.end());
    }
    std::array<std::string, 9> predefRegions = {"Standard", "Gap00", "Gap02", "Gap04", "Gap06", "Gap08", "Gap10", "Gap12", "Gap14"};
    auto* oba = new TObjArray();
    if (cfgGFWSettings.cfgUseGFW) {
      readMatrix(cfgGFWSettings.cfgDefineGFWRegions->getData(), gfwRegions); // read matrix to get initial settings for GFW regions
      LOGF(info, "GFW regions Defined:");
      for (int iRegion = -1; iRegion < EtaGapCounter; iRegion++) {
        if (gfwRegions[iRegion + 1][0] != 0) {
          setGFWRegionPredefined(iRegion);
        }
      }
      if (cfgGFWSettings.cfgRegions->GetSize() == -1) { // define regions for correlation
        LOGF(warning, "cfgRegions: All columns used for GFW regions are not the same size!");
      } else if (cfgGFWSettings.cfgRegions->GetSize() == 0) {
        LOGF(info, "cfgRegions: No extra GFW regions to be added");
      } else {
        cacheRegions.SetNames(cfgGFWSettings.cfgRegions->GetNames());
        cacheRegions.SetEtaMin(cfgGFWSettings.cfgRegions->GetEtaMin());
        cacheRegions.SetEtaMax(cfgGFWSettings.cfgRegions->GetEtaMax());
        cacheRegions.SetpTDifs(cfgGFWSettings.cfgRegions->GetpTDifs());
        cacheRegions.SetBitmasks(cfgGFWSettings.cfgRegions->GetBitmasks());
        for (int i = 0; i < cacheRegions.GetSize(); i++) {
          fGFW->AddRegion(cacheRegions.GetNames().at(i), cacheRegions.GetEtaMin().at(i), cacheRegions.GetEtaMax().at(i), (cacheRegions.GetpTDifs().at(i) != 0) ? 1 + axisPtGFW->GetNbins() : 1, cacheRegions.GetBitmasks().at(i));
          LOGF(info, "cfgRegion added: %s, Eta Min: %f, Eta Max: %f, pT diff: %i, bitmask: %i", cacheRegions.GetNames().at(i), cacheRegions.GetEtaMin().at(i), cacheRegions.GetEtaMax().at(i), cacheRegions.GetpTDifs().at(i), cacheRegions.GetBitmasks().at(i));
        }
      }
      if (cfgGFWSettings.cfgCorrConfig->GetSize() == -1) { // define N particle correlators
        LOGF(fatal, "cfgCorrConfig: All columns used for GFW correlations are not the same size!");
      } else if (cfgGFWSettings.cfgCorrConfig->GetSize() == 0) {
        LOGF(info, "cfgCorrConfig: No GFW regions have been added");
      } else {
        cacheCorrConfig.SetCorrs(cfgGFWSettings.cfgCorrConfig->GetCorrs());
        cacheCorrConfig.SetHeads(cfgGFWSettings.cfgCorrConfig->GetHeads());
        cacheCorrConfig.SetpTDifs(cfgGFWSettings.cfgCorrConfig->GetpTDifs());
        cacheCorrConfig.SetpTCorrMasks(cfgGFWSettings.cfgCorrConfig->GetpTCorrMasks());
        std::string suffix = "_ptDiff";
        for (int i = 0; i < cacheCorrConfig.GetSize(); i++) {
          corrConfigs.push_back(fGFW->GetCorrelatorConfig(cacheCorrConfig.GetCorrs().at(i), cacheCorrConfig.GetHeads().at(i), cacheCorrConfig.GetpTDifs().at(i) != 0));
          if (cacheCorrConfig.GetpTDifs().at(i) != 0) {
            for (auto iBin = 1; iBin <= axisPtGFW->GetNbins(); iBin++) {
              std::string index = Form("_pt_%i", iBin);
              oba->Add(new TNamed(cacheCorrConfig.GetHeads().at(i) + index, cacheCorrConfig.GetHeads().at(i) + suffix));
              std::string sName = cacheCorrConfig.GetHeads().at(i) + index;
              LOGF(info, Form("Correlator to be calculated: %s", sName.c_str()));
            }
          } else {
            oba->Add(new TNamed(cacheCorrConfig.GetHeads().at(i), cacheCorrConfig.GetHeads().at(i)));
            LOGF(info, Form("Correlator to be calculated: %s", cacheCorrConfig.GetHeads().at(i).c_str()));
          }
        }
      }
      fGFW->CreateRegions();
      fFC.setObject(new FlowContainer("FlowContainer"));
      fFC->SetXAxis(axisPtGFW);
      fFC->Initialize(oba, axisGFW, cfgGFWSettings.cfgBootstrapSamples);
    }
    delete oba;
    if (cfgTriggeredData.cfgUseHighMultTrigger) {
      zorroSummary.setObject(zorro.getZorroSummary());
    }
  }
  template <typename TCollision>
  float getCentrality(const TCollision& collision)
  {
    switch (cfgEventSelection.cfgCentEstimator) {
      case CentFT0C:
        return collision.centFT0C();
      case CentFT0CVariant1:
        return collision.centFT0CVariant1();
      case CentFT0M:
        return collision.centFT0M();
      case CentFV0A:
        return collision.centFV0A();
      case CentNTPV:
        return collision.centNTPV();
      case CentNGlobal:
        return collision.centNGlobal();
      case CentMFT:
        return collision.centMFT();
      default:
        return collision.centFT0C();
    }
  }
  template <typename TCollision>
  bool selectMultiplicityCorrelation(const TCollision& collision, const int& multTrk, const float& centrality)
  {
    auto multNTracksPV = collision.multNTracksPV();
    if (cfgMultiplicityCorrelationCut.cfgGlobalMultPVMultCutEnabled) {
      if (multTrk > fGlobalMultPVMultHighCut->Eval(multNTracksPV) || multTrk < fGlobalMultPVMultLowCut->Eval(multNTracksPV)) {
        return false;
      }
    }
    if (cfgMultiplicityCorrelationCut.cfgGlobalMultV0AMultCutEnabled) {
      if (collision.multFV0A() > fGlobalMultV0AMultHighCut->Eval(multTrk) || collision.multFV0A() < fGlobalMultV0AMultLowCut->Eval(multTrk)) {
        return false;
      }
    }
    if (cfgMultiplicityCorrelationCut.cfgGlobalMultCentCutEnabled) {
      if (multTrk > fGlobalMultCentHighCut->Eval(centrality) || multTrk < fGlobalMultCentLowCut->Eval(centrality)) {
        return false;
      }
    }
    if (cfgMultiplicityCorrelationCut.cfgPVMultCentCutEnabled) {
      if (multNTracksPV > fPVMultCentHighCut->Eval(centrality) || multNTracksPV < fPVMultCentLowCut->Eval(centrality)) {
        return false;
      }
    }
    return true;
  }
  template <typename TCollision>
  bool eventSelected(const TCollision& collision, const int& multTrk, const float& centrality, const int& runNumber)
  {
    if (!cfgEventSelection.cfgEvSelUsage) {
      return true;
    }
    if (cfgEventSelection.cfgEvSel8) {
      if (!collision.sel8()) {
        return false;
      }
      registryQA.fill(HIST("hEventCount"), Sel8);
    }
    if (!cfgEventSelection.cfgEvSelRCTflags.value.empty()) {
      if (!rctChecker(collision)) {
        return false;
      }
      registryQA.fill(HIST("hEventCount"), RCTflags);
    }
    for (const auto& cut : eventSelectionCuts) {
      if (!cut.enabled) {
        continue;
      }
      if (!collision.selection_bit(cut.flag)) {
        return false;
      }
      registryQA.fill(HIST("hEventCount"), cut.histBin);
    }
    if (cfgMultiplicityCorrelationCut.cfgEvSelMultCorrelation) {
      if (!selectMultiplicityCorrelation(collision, multTrk, centrality)) {
        return false;
      }
      registryQA.fill(HIST("hEventCount"), MultCorrelation);
    }
    if (std::find(cfgRemovedRunList.value.begin(), cfgRemovedRunList.value.end(), runNumber) != cfgRemovedRunList.value.end()) {
      return false;
    }
    registryQA.fill(HIST("hEventCount"), BadRunsRemoved);
    return true;
  }
  template <TrackSelectionGroup sg, typename TTracks>
  bool trackSelected(const TTracks& track)
  {
    if (trackSel[DCAz][sg] == 0 && std::fabs(track.dcaZ()) > fDCAzCut[sg]->Eval(track.pt())) {
      return false;
    }
    if (trackSel[DCAz][sg] != 0 && std::fabs(track.dcaZ()) > trackSel[DCAz][sg]) {
      return false;
    }
    return ((track.tpcChi2NCl() < trackSel[Chi2prTPCcls][sg]) && (track.tpcNClsFound() >= trackSel[TPCclu][sg]) && (track.itsNCls() >= trackSel[ITSclu][sg]) && (track.tpcNClsCrossedRows() >= trackSel[TPCCrossedRows][sg]) && (std::fabs(track.dcaXY()) < fDCAxyCut[sg]->Eval(track.pt())));
  }
  template <QAFillTime ft, typename TCollision, typename TTracks>
  inline void fillEventQA(const TCollision& collision, const TTracks& tracks)
  {
    registryQA.fill(HIST(FillTimeName[ft]) + HIST("globalTracks_centT0C"), collision.centFT0C(), tracks.size());
    registryQA.fill(HIST(FillTimeName[ft]) + HIST("PVTracks_centT0C"), collision.centFT0C(), collision.multNTracksPV());
    registryQA.fill(HIST(FillTimeName[ft]) + HIST("globalTracks_PVTracks"), collision.multNTracksPV(), tracks.size());
    registryQA.fill(HIST(FillTimeName[ft]) + HIST("globalTracks_multT0A"), collision.multFT0A(), tracks.size());
    registryQA.fill(HIST(FillTimeName[ft]) + HIST("globalTracks_multV0A"), collision.multFV0A(), tracks.size());
    registryQA.fill(HIST(FillTimeName[ft]) + HIST("multV0A_multT0A"), collision.multFT0A(), collision.multFV0A());
    registryQA.fill(HIST(FillTimeName[ft]) + HIST("multT0C_centT0C"), collision.centFT0C(), collision.multFT0C());
  }
  template <QAFillTime ft, typename TCollision>
  inline void fillMCEventQA(const TCollision& collision, const int& ntrk)
  {
    registryQA.fill(HIST(FillTimeName[ft]) + HIST("globalTracks_centT0C"), collision.centFT0C(), ntrk);
    registryQA.fill(HIST(FillTimeName[ft]) + HIST("PVTracks_centT0C"), collision.centFT0C(), collision.multNTracksPV());
    registryQA.fill(HIST(FillTimeName[ft]) + HIST("globalTracks_PVTracks"), collision.multNTracksPV(), ntrk);
    registryQA.fill(HIST(FillTimeName[ft]) + HIST("globalTracks_multT0A"), collision.multFT0A(), ntrk);
    registryQA.fill(HIST(FillTimeName[ft]) + HIST("globalTracks_multV0A"), collision.multFV0A(), ntrk);
    registryQA.fill(HIST(FillTimeName[ft]) + HIST("multV0A_multT0A"), collision.multFT0A(), collision.multFV0A());
    registryQA.fill(HIST(FillTimeName[ft]) + HIST("multT0C_centT0C"), collision.centFT0C(), collision.multFT0C());
  }
  template <TrackSelectionGroup sg, typename TTracks>
  inline void fillTrackQA(const TTracks& track, const float& vtxz)
  {
    registryQA.fill(HIST(FillTrackName[sg]) + HIST("hTrackPhi"), track.phi());
    registryQA.fill(HIST(FillTrackName[sg]) + HIST("hTrackEta"), track.eta());
    registryQA.fill(HIST(FillTrackName[sg]) + HIST("hPhiEtaVtxZ"), track.phi(), track.eta(), vtxz);
    registryQA.fill(HIST(FillTrackName[sg]) + HIST("hPt"), track.pt());
    registryQA.fill(HIST(FillTrackName[sg]) + HIST("hChi2prTPCcls"), track.tpcChi2NCl());
    registryQA.fill(HIST(FillTrackName[sg]) + HIST("hChi2prITScls"), track.itsChi2NCl());
    registryQA.fill(HIST(FillTrackName[sg]) + HIST("hnTPCClu"), track.tpcNClsFound());
    registryQA.fill(HIST(FillTrackName[sg]) + HIST("hnITSClu"), track.itsNCls());
    registryQA.fill(HIST(FillTrackName[sg]) + HIST("hnTPCCrossedRow"), track.tpcNClsCrossedRows());
    registryQA.fill(HIST(FillTrackName[sg]) + HIST("hDCAz"), track.dcaZ(), track.pt());
    registryQA.fill(HIST(FillTrackName[sg]) + HIST("hDCAxy"), track.dcaXY(), track.pt());
  }
  template <typename TTrack, typename T>
  inline void fillCorrProfile(const TTrack& track, std::array<T, EtaGapCounter>& trackRegion)
  {
    float eta = track.eta();
    for (uint iRegion = 0; iRegion < static_cast<int>(EtaGapCounter); iRegion++) {
      if (std::fabs(eta) > trackRegion[iRegion].gap) {
        if (eta > 0) {
          trackRegion[iRegion].forward++;
        } else {
          trackRegion[iRegion].backward++;
        }
      }
    }
  }
  void loadCorrection(const int& runNumber)
  {
    if (!cfgGFWSettings.cfgAcceptance.value.empty()) {
      weight.hAcceptance = ccdb->getForRun<GFWWeights>(cfgGFWSettings.cfgAcceptance, runNumber);
      if (weight.hAcceptance == nullptr) {
        LOGF(fatal, "Could not load acceptance weight for run %i from %s", runNumber, cfgGFWSettings.cfgAcceptance.value.c_str());
      } else {
        LOGF(info, "Loaded acceptance weight for run %i from %s", runNumber, cfgGFWSettings.cfgAcceptance.value.c_str());
      }
    }
    if (!cfgGFWSettings.cfgEfficiency.value.empty()) {
      weight.hEfficiency = ccdb->getForRun<TH1D>(cfgGFWSettings.cfgEfficiency, runNumber);
      if (weight.hEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency weight for run %i from %s", runNumber, cfgGFWSettings.cfgEfficiency.value.c_str());
      } else {
        LOGF(info, "Loaded efficiency weight for run %i from %s", runNumber, cfgGFWSettings.cfgEfficiency.value.c_str());
      }
    }
    if (!cfgGFWSettings.cfgEfficiencyNch.value.empty() && cfgGFWSettings.cfgUseNch) {
      weight.hEfficiencyNch = ccdb->getForRun<TH1D>(cfgGFWSettings.cfgEfficiencyNch, runNumber);
      if (weight.hEfficiencyNch == nullptr) {
        LOGF(fatal, "Could not load Nch-efficiency weight for run %i from %s - No correction on x-axis is applied", runNumber, cfgGFWSettings.cfgEfficiencyNch.value.c_str());
      } else {
        LOGF(info, "Loaded Nch-efficiency weight for run %i from %s", runNumber, cfgGFWSettings.cfgEfficiencyNch.value.c_str());
      }
    }
  }
  inline void fillFlowContainer(const GFW::CorrConfig& nCorrConfig, const double& xAxis, const int& iBin, const double& rndm, const std::string& corrName)
  {
    auto dnx = fGFW->Calculate(nCorrConfig, iBin, kTRUE).real();
    if (dnx == 0) {
      return;
    }
    auto val = fGFW->Calculate(nCorrConfig, iBin, kFALSE).real() / dnx;
    if (std::abs(val) < 1) {
      fFC->FillProfile(corrName.c_str(), xAxis, val, cfgGFWSettings.cfgUseMultWeights ? dnx : 1.0, rndm);
    }
  }
  template <typename T>
  inline void fillCorelationCounter(const int& nchRaw, const std::array<T, EtaGapCounter>& trackCorr)
  {
    float binContent = 0;
    float nCorrelations = 0;
    for (uint iRegion = 0; iRegion < static_cast<int>(EtaGapCounter); iRegion++) {
      if (trackCorr[iRegion].forward > static_cast<float>(Min2PC) && trackCorr[iRegion].backward > static_cast<float>(Min2PC)) {
        binContent = registryCorr.get<TH2>(HIST("hTwoParticleCorrelations"))->GetBinContent(nchRaw, iRegion + 1);
        nCorrelations = trackCorr[iRegion].forward * trackCorr[iRegion].backward;
        registryCorr.get<TH2>(HIST("hTwoParticleCorrelations"))->SetBinContent(nchRaw, iRegion + 1, binContent + nCorrelations);
      }
      if (trackCorr[iRegion].forward > static_cast<float>(Min4PC) && trackCorr[iRegion].backward > static_cast<float>(Min4PC)) {
        binContent = registryCorr.get<TH2>(HIST("hFourParticleCorrelations"))->GetBinContent(nchRaw, iRegion + 1);
        nCorrelations = trackCorr[iRegion].forward * (trackCorr[iRegion].forward - 1) * trackCorr[iRegion].backward * (trackCorr[iRegion].backward - 1) / 4.;
        registryCorr.get<TH2>(HIST("hFourParticleCorrelations"))->SetBinContent(nchRaw, iRegion + 1, binContent + nCorrelations);
      }
      if (trackCorr[iRegion].forward > static_cast<float>(Min6PC) && trackCorr[iRegion].backward > static_cast<float>(Min6PC)) {
        binContent = registryCorr.get<TH2>(HIST("hSixParticleCorrelations"))->GetBinContent(nchRaw, iRegion + 1);
        nCorrelations = trackCorr[iRegion].forward * (trackCorr[iRegion].forward - 1) * (trackCorr[iRegion].forward - 2) * trackCorr[iRegion].backward * (trackCorr[iRegion].backward - 1) * (trackCorr[iRegion].backward - 2) / 36.;
        registryCorr.get<TH2>(HIST("hSixParticleCorrelations"))->SetBinContent(nchRaw, iRegion + 1, binContent + nCorrelations);
      }
      if (trackCorr[iRegion].forward > static_cast<float>(Min8PC) && trackCorr[iRegion].backward > static_cast<float>(Min8PC)) {
        binContent = registryCorr.get<TH2>(HIST("hEightParticleCorrelations"))->GetBinContent(nchRaw, iRegion + 1);
        nCorrelations = trackCorr[iRegion].forward * (trackCorr[iRegion].forward - 1) * (trackCorr[iRegion].forward - 2) * (trackCorr[iRegion].forward - 3) * trackCorr[iRegion].backward * (trackCorr[iRegion].backward - 1) * (trackCorr[iRegion].backward - 2) * (trackCorr[iRegion].backward - 3) / 576;
        registryCorr.get<TH2>(HIST("hEightParticleCorrelations"))->SetBinContent(nchRaw, iRegion + 1, binContent + nCorrelations);
      }
    }
  }
  template <typename T>
  void fillOutputObjects(const float& cent, const int& nchRaw, const float& nchCorrected, const double& rndm, const std::array<T, EtaGapCounter>& trackCorr)
  {
    if (cfgGetNumberOfCorr) {
      fillCorelationCounter(nchRaw, trackCorr);
    }
    if (cfgGFWSettings.cfgUseGFW) {
      double xAxis = (cfgGFWSettings.cfgUseNch) ? nchCorrected : cent;
      for (const auto& corr : corrConfigs) {
        if (!corr.pTDif) {
          fillFlowContainer(corr, xAxis, 0, rndm, corr.Head);
        } else {
          std::string nameFC;
          for (int iBin = 1; iBin <= axisPtGFW->GetNbins(); iBin++) {
            nameFC = corr.Head + "_pt_" + iBin;
            fillFlowContainer(corr, xAxis, iBin, rndm, nameFC);
          }
        }
      }
    }
  }
  double getTimeSinceStartOfFill(const uint64_t& timestamp, const int& firstRun)
  {
    auto runDuration = ccdb->getRunDuration(firstRun);
    uint64_t tsSOF = runDuration.first;
    uint64_t diff = timestamp - tsSOF;
    return static_cast<double>(diff) / 3600000.0;
  }
  template <typename TTrack>
  inline void fillGFW(const TTrack& track, const float& vtxZ)
  {
    bool isRef = (track.pt() < cfgGFWSettings.cfgCutPtRefMax && track.pt() > cfgGFWSettings.cfgCutPtRefMin);
    bool isPOI = (track.pt() < cfgGFWSettings.cfgCutPtPOIMax && track.pt() > cfgGFWSettings.cfgCutPtPOIMin);

    if (!isRef && !isPOI) {
      return;
    }
    double wNUE = (weight.hEfficiency != nullptr) ? 1. / weight.hEfficiency->GetBinContent(weight.hEfficiency->FindBin(track.pt())) : 1.0;
    double wNUA = (weight.hAcceptance != nullptr) ? weight.hAcceptance->getNUA(track.phi(), track.eta(), vtxZ) : 1.0;
    if (weight.hAcceptance != nullptr) {
      registryQA.fill(HIST("cutObs/hPhiCorrected"), track.phi(), wNUA);
    }
    if (isRef) {
      fGFW->Fill(track.eta(), axisPtGFW->FindBin(track.pt()) - 1, track.phi(), wNUA * wNUE, 1);
    }
    if (isPOI) {
      fGFW->Fill(track.eta(), axisPtGFW->FindBin(track.pt()) - 1, track.phi(), wNUA * wNUE, 2);
    }
    if (isRef && isPOI) {
      fGFW->Fill(track.eta(), axisPtGFW->FindBin(track.pt()) - 1, track.phi(), wNUA * wNUE, 4);
    }
  }
  using FilteredTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA>>;
  using FilteredCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs, aod::CentNGlobals, aod::CentMFTs, aod::Mults>>;
  void processData(FilteredCollisions::iterator const& collision, aod::BCsWithTimestamps const&, FilteredTracks const& tracks)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int run = bc.runNumber();
    float nch = 0;
    int nchCorrCounter = 0;
    float nchCorrected = 0;
    double lRandom = fRndm->Rndm();
    std::array<TrackRegion, 8> trackCorr{TrackRegion{.forward = 0, .backward = 0, .gap = 0.}, TrackRegion{.forward = 0, .backward = 0, .gap = 0.1}, TrackRegion{.forward = 0, .backward = 0, .gap = 0.2}, TrackRegion{.forward = 0, .backward = 0, .gap = 0.3}, TrackRegion{.forward = 0, .backward = 0, .gap = 0.4}, TrackRegion{.forward = 0, .backward = 0, .gap = 0.5}, TrackRegion{.forward = 0, .backward = 0, .gap = 0.6}, TrackRegion{.forward = 0, .backward = 0, .gap = 0.7}};
    if (cfgTriggeredData.cfgUseHighMultTrigger) {
      zorro.initCCDB(ccdb.service, run, bc.timestamp(), triggerName[cfgTriggeredData.cfgHighMultiplicityTriggerMenu]);
      zorro.populateHistRegistry(registryTrigger, run);
      if (!zorro.isSelected(bc.globalBC())) {
        return;
      }
    }
    registryQA.fill(HIST("hEventCount"), FilteredEvent);
    if (tracks.size() < 1) {
      return;
    }
    fillEventQA<BeforeSel>(collision, tracks);
    const float cent = getCentrality(collision);
    const float vtxz = collision.posZ();
    if (!eventSelected(collision, tracks.size(), cent, run)) {
      return;
    }
    if (cfgEventSelection.cfgEvSelUsage) {
      fillEventQA<AfterSel>(collision, tracks);
    }
    if (!cfgFirstFillRunsList.value.empty()) {
      if (std::find(cfgFirstFillRunsList.value.begin(), cfgFirstFillRunsList.value.end(), run) != cfgFirstFillRunsList.value.end() && firstRunOfFill != run) {
        firstRunOfFill = run;
      }
      registryQA.fill(HIST("multT0A_multT0C"), collision.multFT0A(), collision.multFT0C(), getTimeSinceStartOfFill(bc.timestamp(), firstRunOfFill));
    }
    registryQA.fill(HIST("hVtxZ"), vtxz);
    registryQA.fill(HIST("hCentrality"), cent);
    if (std::find(runNumbers.begin(), runNumbers.end(), run) == runNumbers.end()) {
      if (cfgNUAWeightsRunByRun) {
        std::shared_ptr<TH3> histNUA = registryQA.add<TH3>(Form("NUA/hPhiEtaVtxZ_%i", run), "", {HistType::kTH3D, {{60, 0.0, constants::math::TwoPI}, {32, -0.8, 0.8}, {40, -10, 10}}});
        th3sPerRun.insert(std::make_pair(run, histNUA));
      }
      if (cfgGFWSettings.cfgUseGFW) {
        loadCorrection(run);
      }
      runNumbers.push_back(run);
    }
    if (cfgGFWSettings.cfgUseGFW) {
      fGFW->Clear();
    }
    for (const auto& track : tracks) {
      if (trackSelected<TrkSelNch>(track)) {
        fillTrackQA<TrkSelNch>(track, vtxz);
        if (weight.hEfficiencyNch != nullptr) {
          nchCorrected += 1. / weight.hEfficiency->GetBinContent(weight.hEfficiencyNch->FindBin(track.pt()));
        } else {
          nchCorrected++;
        }
        nch++;
      }
      if (!trackSelected<TrkSelObs>(track)) {
        continue;
      }
      fillTrackQA<TrkSelObs>(track, vtxz);
      if (cfgGetNumberOfCorr) {
        fillCorrProfile(track, trackCorr);
        nchCorrCounter++;
      }
      if (cfgNUAWeightsRunByRun) {
        th3sPerRun[run]->Fill(track.phi(), track.eta(), vtxz);
      }
      if (cfgGFWSettings.cfgUseGFW) {
        fillGFW(track, vtxz);
      }
    }
    if (!cfgGFWSettings.cfgEfficiencyNch.value.empty() && weight.hEfficiencyNch != nullptr) {
      registryQA.fill(HIST("hMultCorrected_MultUncorrected"), nch, nchCorrected);
      registryQA.fill(HIST("hNchCorrected"), nchCorrected);
    }
    registryQA.fill(HIST("hNch"), nch);
    fillOutputObjects(cent, nchCorrCounter, nchCorrected, lRandom, trackCorr);
  }
  PROCESS_SWITCH(FlowSmallSystems, processData, "Process analysis for non-derived data", false);
  using FilteredTracksMC = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels>>;
  using CollisionsReco = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs, aod::CentNGlobals, aod::CentMFTs, aod::Mults>>;
  using FilteredCollisionsMC = soa::Filtered<aod::McCollisions>;
  void processMC(aod::McCollisions::iterator const&, CollisionsReco const& collisions, FilteredTracksMC const& tracks, aod::McParticles const& particlesMC, aod::BCsWithTimestamps const&)
  {
    int run = 0;
    int nchReco = 0;
    int nchGen = 0;
    float cent = 0;
    float vtxz = 0;
    if (collisions.size() != 1) {
      return;
    }
    for (const auto& particleMC : particlesMC) {
      if (!particleMC.isPhysicalPrimary()) {
        continue;
      }
      auto pdgParticle = pdgDB->GetParticle(particleMC.pdgCode());
      if (pdgParticle == 0 || pdgParticle->Charge() == 0) {
        continue;
      }
      if (std::fabs(particleMC.eta()) > cfgCutEta || particleMC.pt() < cfgCutPtMin || particleMC.pt() > cfgCutPtMax) {
        continue;
      }
      nchGen++;
    }
    if (nchGen == 0) {
      return;
    }
    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      run = bc.runNumber();
      cent = getCentrality(collision);
      vtxz = collision.posZ();
      registryQA.fill(HIST("hEventCount"), FilteredEvent);
      if (tracks.size() < 1) {
        LOGF(info, "Size of the event is: %i", tracks.size());
        return;
      }
      for (const auto& track : tracks) {
        if (track.collisionId() != collision.globalIndex()) {
          continue;
        }
        if (!trackSelected<TrkSelObs>(track)) {
          continue;
        }
        if (track.mcParticleId() < 0 || !track.has_mcParticle()) {
          continue;
        }
        auto mcParticle = track.mcParticle();
        auto pdgParticle = pdgDB->GetParticle(mcParticle.pdgCode());
        if (pdgParticle == 0) {
          continue;
        }
        if (pdgParticle->Charge() == 0.) {
          continue;
        }
        if (!mcParticle.isPhysicalPrimary()) {
          continue;
        }
        if (std::fabs(mcParticle.eta()) > cfgCutEta || mcParticle.pt() < cfgCutPtMin || mcParticle.pt() > cfgCutPtMax) {
          continue;
        }
        fillTrackQA<TrkReco>(track, vtxz);
        nchReco++;
      }
      if (nchReco == 0) {
        return;
      }
      fillMCEventQA<BeforeSel>(collision, nchReco);
      if (cfgEventSelection.cfgEvSelUsage) {
        if (!eventSelected(collision, nchReco, cent, run)) {
          return;
        }
        fillMCEventQA<AfterSel>(collision, nchReco);
      }
      registryQA.fill(HIST("hEventTrueVsGen"), 0.5);
      registryQA.fill(HIST("hVtxZ"), vtxz);
      registryQA.fill(HIST("hCentrality"), cent);
      for (const auto& track : tracks) {
        if (track.collisionId() != collision.globalIndex()) {
          continue;
        }
        if (!trackSelected<TrkSelObs>(track)) {
          continue;
        }
        if (track.mcParticleId() < 0 || !track.has_mcParticle()) {
          continue;
        }
        auto mcParticle = track.mcParticle();
        auto pdgParticle = pdgDB->GetParticle(mcParticle.pdgCode());
        if (pdgParticle == 0) {
          continue;
        }
        if (pdgParticle->Charge() == 0.) {
          continue;
        }
        if (!mcParticle.isPhysicalPrimary()) {
          continue;
        }
        if (std::fabs(mcParticle.eta()) > cfgCutEta || mcParticle.pt() < cfgCutPtMin || mcParticle.pt() > cfgCutPtMax) {
          continue;
        }
        registryQA.fill(HIST("Reco/hPtCent"), track.pt(), cent);
        registryQA.fill(HIST("Reco/hPtNch"), track.pt(), nchReco);
      }
      registryQA.fill(HIST("hNch"), nchReco);
    }
    for (const auto& particleMC : particlesMC) {
      if (!particleMC.isPhysicalPrimary()) {
        continue;
      }
      auto pdgParticle = pdgDB->GetParticle(particleMC.pdgCode());
      if (!pdgParticle || pdgParticle->Charge() == 0) {
        continue;
      }
      if (std::fabs(particleMC.eta()) > cfgCutEta || particleMC.pt() < cfgCutPtMin || particleMC.pt() > cfgCutPtMax) {
        continue;
      }
      registryQA.fill(HIST("Gen/hTrackPhi"), particleMC.phi());
      registryQA.fill(HIST("Gen/hTrackEta"), particleMC.eta());
      registryQA.fill(HIST("Gen/hPt"), particleMC.pt());
      registryQA.fill(HIST("Gen/hPtNch"), particleMC.pt(), nchReco);
      registryQA.fill(HIST("Gen/hPtCent"), particleMC.pt(), cent);
    }
    registryQA.fill(HIST("hGenVsReco"), nchGen, nchReco);
    registryQA.fill(HIST("hEventTrueVsGen"), 1.5);
    if (std::find(runNumbers.begin(), runNumbers.end(), run) == runNumbers.end()) {
      std::shared_ptr<TH2> hGenVsReco = registryQA.add<TH2>(Form("GenVsReco/%i", run), ";N_{ch} Gen; N_{ch} Reco", {HistType::kTH2D, {{static_cast<int>(binningNch.value.back()), 0., binningNch.value.back()}, {static_cast<int>(binningNch.value.back()), 0., binningNch.value.back()}}});
      th2sPerRun.insert(std::make_pair(run, hGenVsReco));
      runNumbers.push_back(run);
    }
    th2sPerRun[run]->Fill(nchGen, nchReco);
  }
  PROCESS_SWITCH(FlowSmallSystems, processMC, "Process analysis for MC gen and reco levels", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowSmallSystems>(cfgc)};
}
