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

/// \file   flowTransmutation.cxx
/// \author Vojtech Machacek (vojtech.machacek@cern.ch)
/// \since  Oct/05/2026
/// \brief  Task to measure flow observables with cumulants method in transmuted events.

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
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
#include <Framework/runDataProcessing.h>

#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TString.h>

#include <sys/types.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;
#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP}; // NOLINT(bugprone-macro-parentheses)
static constexpr std::array<std::array<float, 2>, 7> LongArrayFloat = {{{{2.5, 2.5}}, {{50.0, 50.0}}, {{70.0, 70.0}}, {{5.0, 5.0}}, {{7.0, 7.0}}, {{7.0, 7.0}}, {{0., 0.}}}};
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
struct FlowTransmutation {
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
    Configurable<std::vector<double>> cfgPVMultCentCutCutPars{"cfgPVMultCentCutPars", std::vector<double>{245.882, -6.90255, 0.0930817, -0.000710842, 0.00000231095, 35.9647, -0.530632, 0.00513954, -0.0000411738, 1.10106e-7}, "Primary vertex multiplicity vs centrality cut parameter values"};
  } cfgMultiplicityCorrelationCut;
  O2_DEFINE_CONFIGURABLE(cfgNUAWeightsRunByRun, bool, true, "Fill NUA weights run by run")
  O2_DEFINE_CONFIGURABLE(cfgGetNumberOfCorr, bool, false, "Get histograms with number of correlations")
  Configurable<std::vector<int>> cfgFirstFillRunsList{"cfgFirstFillRunsList", std::vector<int>{}, "First run of the fill for time-dependent analysis, keep empty to disable"};
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 10.0f, "Maximal pT for all tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for all tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 3.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxy, std::string, "(0.0026+0.005/(x^1.01))", "function for pT dependent DCAxy cut, keep empty for const. cut based on DCAxy (Nsigma)")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, std::string, "(0.0026+0.005/(x^1.01))", "function for pT dependent DCAz cut, keep empty for const. cut based on DCAz (Nsigma)")
  Configurable<std::vector<int>> cfgRemovedRunList{"cfgRemovedRunList", std::vector<int>{-1}, "excluded run numbers"};
  Configurable<LabeledArray<float>> cfgTrackCuts{"cfgTrackCuts", {LongArrayFloat.front().data(), 7, 2, {"chi2 per TPCcls", "TPC cluster", "TPC crossed rows", "ITS cluster", "DCAxy (Nsigma)", "DCAz (Nsigma)", "DCAz (cm) override"}, {"Nch", "Observable"}}, "separate Nch and observable track selections"};
  ConfigurableAxis binningNch{"binningNch", {500, 0, 500}, "binning for N_{ch} histograms"};
  ConfigurableAxis binningTime{"binningTime", {12, 0, 12}, "binning for time-dependent histograms transmutation QA"};
  ConfigurableAxis binningPt{"binningPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10}, "binning for pT histograms"};
  HistogramRegistry registryQA{"registryQA"};
  HistogramRegistry registryCorr{"registryCorr"};
  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgEventSelection.cfgCutVertex);
  Filter trackFilter = ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t)true)) && (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax);
  Service<ccdb::BasicCCDBManager> ccdb{};
  RCTFlagsChecker rctChecker{"CBT"};
  enum TrackSelectionGroup {
    TrkSelNch = 0,
    TrkSelObs,
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
  static constexpr std::array<std::string_view, 2> FillTimeName = {"beforeEventSelQA/", "afterEventSelQA/"};
  static constexpr std::array<std::string_view, 2> FillTrackName = {"cutNch/", "cutObs/"};
  std::array<EventCut, 8> eventSelectionCuts{};
  std::array<std::array<float, 2>, 7> trackSel{};
  TF1* fGlobalMultPVMultLowCut = nullptr;
  TF1* fGlobalMultPVMultHighCut = nullptr;
  TF1* fGlobalMultV0AMultLowCut = nullptr;
  TF1* fGlobalMultV0AMultHighCut = nullptr;
  TF1* fGlobalMultCentLowCut = nullptr;
  TF1* fGlobalMultCentHighCut = nullptr;
  TF1* fPVMultCentLowCut = nullptr;
  TF1* fPVMultCentHighCut = nullptr;
  TH1D* hEff = nullptr;
  std::array<TF1*, 2> fDCAzCut{};
  std::array<TF1*, 2> fDCAxyCut{};
  std::vector<int> runNumbers;
  std::map<int, std::shared_ptr<TH3>> th3sPerRun;
  int firstRunOfFill = -1;
  int getCorrectITSLayersFlag(int flag)
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
  void init(InitContext const&)
  {
    AxisSpec axisVertex{40, -20, 20, "Vtxz (cm)"};
    AxisSpec axisPhi{60, 0.0, constants::math::TwoPI, "#varphi"};
    AxisSpec axisEta{40, -1., 1., "#eta"};
    AxisSpec axisCent{100, 0., 100., "Centrality"};
    AxisSpec axisT0C = {1000, 0, 50000, "N_{ch} (T0C)"};
    AxisSpec axisT0A = {1800, 0, 180000, "N_{ch} (T0A)"};
    AxisSpec axisV0A = {1800, 0, 180000, "N_{ch} (V0A)"};
    AxisSpec axisNch = {binningNch, "N_{ch}"};
    AxisSpec axisPt = {binningPt, "p_{T} (GeV/c)"};
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
    registryQA.add("hEventCount", "Number of Event;; Count", {HistType::kTH1D, {{13, 0.5, 13.5}}});
    registryQA.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(FilteredEvent, "Filtered event");
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
          fPVMultCentHighCut->SetParameters(cfgMultiplicityCorrelationCut.cfgPVMultCentCutCutPars.value.data());
          fPVMultCentLowCut = new TF1("fPVMultCentLowCut", cfgMultiplicityCorrelationCut.cfgPVMultCentLowCutFunction->c_str(), 0, 100);
          fPVMultCentLowCut->SetParameters(cfgMultiplicityCorrelationCut.cfgPVMultCentCutCutPars.value.data());
          LOGF(info, "PV multiplicity vs centrality Upper Cut: %s", cfgMultiplicityCorrelationCut.cfgPVMultCentHighCutFunction->c_str());
          LOGF(info, "PV multiplicity vs centrality Lower Cut: %s", cfgMultiplicityCorrelationCut.cfgPVMultCentLowCutFunction->c_str());
          LOGF(info, "PV multiplicity vs centrality cut parameters:");
          for (std::size_t i = 0; i < cfgMultiplicityCorrelationCut.cfgPVMultCentCutCutPars.value.size(); ++i) {
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
  bool selectMultiplicityCorrelation(const TCollision& collision, const int multTrk, const float& centrality)
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
  bool eventSelected(const TCollision& collision, const int multTrk, const float centrality, const int runNumber)
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
    for (auto& cut : eventSelectionCuts) {
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
  double getTimeSinceStartOfFill(uint64_t timestamp, int firstRun)
  {
    auto runDuration = ccdb->getRunDuration(firstRun);
    uint64_t tsSOF = runDuration.first;
    uint64_t diff = timestamp - tsSOF;
    return static_cast<double>(diff) / 3600000.0;
  }
  using FilteredTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA>>;
  using FilteredCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs, aod::CentNGlobals, aod::CentMFTs, aod::Mults>>;
  void process(FilteredCollisions::iterator const& collision, aod::BCsWithTimestamps const&, FilteredTracks const& tracks)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int run = bc.runNumber();
    int nch = 0;
    std::array<TrackRegion, 8> trackCorr{TrackRegion{.forward = 0, .backward = 0, .gap = 0.}, TrackRegion{.forward = 0, .backward = 0, .gap = 0.1}, TrackRegion{.forward = 0, .backward = 0, .gap = 0.2}, TrackRegion{.forward = 0, .backward = 0, .gap = 0.3}, TrackRegion{.forward = 0, .backward = 0, .gap = 0.4}, TrackRegion{.forward = 0, .backward = 0, .gap = 0.5}, TrackRegion{.forward = 0, .backward = 0, .gap = 0.6}, TrackRegion{.forward = 0, .backward = 0, .gap = 0.7}};
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
      runNumbers.push_back(run);
    }
    for (const auto& track : tracks) {

      if (trackSelected<TrkSelNch>(track)) {
        fillTrackQA<TrkSelNch>(track, vtxz);
      }
      if (!trackSelected<TrkSelObs>(track)) {
        continue;
      }
      fillTrackQA<TrkSelObs>(track, vtxz);
      if (cfgGetNumberOfCorr) {
        fillCorrProfile(track, trackCorr);

        nch++;
      }
      if (cfgNUAWeightsRunByRun) {
        th3sPerRun[run]->Fill(track.phi(), track.eta(), vtxz);
      }
    }
    float binContent = 0;
    float nCorrelations = 0;
    if (cfgGetNumberOfCorr) {
      for (uint iRegion = 0; iRegion < static_cast<int>(EtaGapCounter); iRegion++) {
        if (trackCorr[iRegion].forward > 0 && trackCorr[iRegion].backward > 0) {
          binContent = registryCorr.get<TH2>(HIST("hTwoParticleCorrelations"))->GetBinContent(nch, iRegion + 1);
          nCorrelations = trackCorr[iRegion].forward * trackCorr[iRegion].backward;
          registryCorr.get<TH2>(HIST("hTwoParticleCorrelations"))->SetBinContent(nch, iRegion + 1, binContent + nCorrelations);
        }
        if (trackCorr[iRegion].forward > 1 && trackCorr[iRegion].backward > 1) {
          binContent = registryCorr.get<TH2>(HIST("hFourParticleCorrelations"))->GetBinContent(nch, iRegion + 1);
          nCorrelations = trackCorr[iRegion].forward * (trackCorr[iRegion].forward - 1) * trackCorr[iRegion].backward * (trackCorr[iRegion].backward - 1) / 4.;
          registryCorr.get<TH2>(HIST("hFourParticleCorrelations"))->SetBinContent(nch, iRegion + 1, binContent + nCorrelations);
        }
        if (trackCorr[iRegion].forward > 2 && trackCorr[iRegion].backward > 2) {
          binContent = registryCorr.get<TH2>(HIST("hSixParticleCorrelations"))->GetBinContent(nch, iRegion + 1);
          nCorrelations = trackCorr[iRegion].forward * (trackCorr[iRegion].forward - 1) * (trackCorr[iRegion].forward - 2) * trackCorr[iRegion].backward * (trackCorr[iRegion].backward - 1) * (trackCorr[iRegion].backward - 2) / 36.;
          registryCorr.get<TH2>(HIST("hSixParticleCorrelations"))->SetBinContent(nch, iRegion + 1, binContent + nCorrelations);
        }
        if (trackCorr[iRegion].forward > 3 && trackCorr[iRegion].backward > 3) {
          binContent = registryCorr.get<TH2>(HIST("hEightParticleCorrelations"))->GetBinContent(nch, iRegion + 1);
          nCorrelations = trackCorr[iRegion].forward * (trackCorr[iRegion].forward - 1) * (trackCorr[iRegion].forward - 2) * (trackCorr[iRegion].forward - 3) * trackCorr[iRegion].backward * (trackCorr[iRegion].backward - 1) * (trackCorr[iRegion].backward - 2) * (trackCorr[iRegion].backward - 3) / 576;
          registryCorr.get<TH2>(HIST("hEightParticleCorrelations"))->SetBinContent(nch, iRegion + 1, binContent + nCorrelations);
        }
      }
    }
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowTransmutation>(cfgc)};
}
