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

#include <CCDB/BasicCCDBManager.h>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <string>
#include <memory>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/CCDB/ctpRateFetcher.h"

#include "GFWPowerArray.h"
#include "GFW.h"
#include "GFWCumulant.h"
#include "GFWWeights.h"
#include "FlowContainer.h"
#include "TList.h"
#include <TProfile.h>
#include <TRandom3.h>
#include <TF1.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowTask {

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
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5f, "max chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCclu, float, 70.0f, "minimum TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutITSclu, float, 5.0f, "minimum ITS clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "max DCA to vertex z")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxyppPass3Enabled, bool, false, "switch of ppPass3 DCAxy pt dependent cut")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAzPtDepEnabled, bool, false, "switch of DCAz pt dependent cut")
  O2_DEFINE_CONFIGURABLE(cfgTrkSelSwitch, bool, false, "switch for self-defined track selection")
  O2_DEFINE_CONFIGURABLE(cfgTrkSelRun3ITSMatch, bool, false, "GlobalTrackRun3ITSMatching::Run3ITSall7Layers selection")
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, false, "Use additional event cut on mult correlations")
  O2_DEFINE_CONFIGURABLE(cfgUseTentativeEventCounter, bool, false, "After sel8(), count events regardless of real event selection")
  O2_DEFINE_CONFIGURABLE(cfgEvSelkNoSameBunchPileup, bool, false, "rejects collisions which are associated with the same found-by-T0 bunch crossing")
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
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Use Nch for flow observables")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 30, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeights, bool, false, "Fill and output NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeightsRefPt, bool, false, "NUA weights are filled in ref pt bins")
  O2_DEFINE_CONFIGURABLE(cfgOutputGroupNUAWeights, bool, false, "Fill and output group NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptanceGroup, std::string, "", "CCDB path to group acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptanceGroupUse, bool, false, "Apply group acceptance, this option overrides cfgAcceptance")
  O2_DEFINE_CONFIGURABLE(cfgEvSelOccupancy, bool, true, "Occupancy cut")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyHigh, int, 500, "High cut on TPC occupancy")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyLow, int, 0, "Low cut on TPC occupancy")
  O2_DEFINE_CONFIGURABLE(cfgUseSmallMemory, bool, false, "Use small memory mode")
  Configurable<std::vector<std::string>> cfgUserDefineGFWCorr{"cfgUserDefineGFWCorr", std::vector<std::string>{"refN02 {2} refP02 {-2}", "refN12 {2} refP12 {-2}"}, "User defined GFW CorrelatorConfig"};
  Configurable<std::vector<std::string>> cfgUserDefineGFWName{"cfgUserDefineGFWName", std::vector<std::string>{"Ch02Gap22", "Ch12Gap22"}, "User defined GFW Name"};
  Configurable<std::vector<int>> cfgRunRemoveList{"cfgRunRemoveList", std::vector<int>{-1}, "excluded run numbers"};
  Configurable<std::vector<int>> cfgGroupSplitRunNumber{"cfgGroupSplitRunNumber", std::vector<int>{544510, 544653}, "runnumbers for group splitting (suppose run numbers are increasing monotonically) "};

  ConfigurableAxis axisVertex{"axisVertex", {40, -20, 20}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPtHist{"axisPtHist", {100, 0., 10.}, "pt axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10}, "pt axis for histograms"};
  ConfigurableAxis axisIndependent{"axisIndependent", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "X axis for histograms"};
  ConfigurableAxis axisCentForQA{"axisCentForQA", {100, 0, 100}, "centrality for QA"};
  ConfigurableAxis axisNch{"axisNch", {4000, 0, 4000}, "N_{ch}"};
  ConfigurableAxis axisT0C{"axisT0C", {70, 0, 70000}, "N_{ch} (T0C)"};
  ConfigurableAxis axisT0A{"axisT0A", {200, 0, 200000}, "N_{ch} (T0A)"};
  ConfigurableAxis axisNchPV{"axisNchPV", {4000, 0, 4000}, "N_{ch} (PV)"};
  ConfigurableAxis axisDCAz{"axisDCAz", {200, -2, 2}, "DCA_{z} (cm)"};
  ConfigurableAxis axisDCAxy{"axisDCAxy", {200, -1, 1}, "DCA_{xy} (cm)"};

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgCutVertex) && (aod::cent::centFT0C > cfgCentFT0CMin) && (aod::cent::centFT0C < cfgCentFT0CMax);
  Filter trackFilter = ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  // Corrections
  TH1D* mEfficiency = nullptr;
  GFWWeights* mAcceptance = nullptr;
  TList* mGroupAcceptanceList = nullptr;
  bool correctionsLoaded = false;

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  // Define output
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  OutputObj<GFWWeights> fWeights{GFWWeights("weights")};
  HistogramRegistry registry{"registry"};
  OutputObj<TList> fGroupNUAList{"GroupNUAList", OutputObjHandlingPolicy::AnalysisObject, OutputObjSourceType::OutputObjSource};

  // define global variables
  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  TAxis* fPtAxis;
  TRandom3* fRndm = new TRandom3(0);
  enum CentEstimators {
    kCentFT0C = 0,
    kCentFT0CVariant1,
    kCentFT0M,
    kCentFV0A,
    // Count the total number of enum
    kCount_CentEstimators
  };
  int mRunNumber{-1};
  uint64_t mSOR{0};
  double mMinSeconds{-1.};
  std::unordered_map<int, TH2*> gHadronicRate;
  ctpRateFetcher mRateFetcher;
  TH2* gCurrentHadronicRate;
  std::vector<std::shared_ptr<GFWWeights>> groupNUAWeightPtr;

  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::Mults>>;
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA>>;

  // Track selection
  TrackSelection myTrackSel;
  TF1* fPhiCutLow = nullptr;
  TF1* fPhiCutHigh = nullptr;
  // Additional Event selection cuts - Copy from flowGenericFramework.cxx
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;
  TF1* fT0AV0AMean = nullptr;
  TF1* fT0AV0ASigma = nullptr;

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);

    // Add some output objects to the histogram registry
    // Event QA
    registry.add("hEventCount", "Number of Event;; Count", {HistType::kTH1D, {{5, 0, 5}}});
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(1, "Filtered event");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(2, "after sel8");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(3, "after supicious Runs removal");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(4, "after additional event cut");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(5, "after correction loads");
    registry.add("hEventCountSpecific", "Number of Event;; Count", {HistType::kTH1D, {{10, 0, 10}}});
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(1, "after sel8");
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(2, "kNoSameBunchPileup");
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(3, "kIsGoodZvtxFT0vsPV");
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(4, "kNoCollInTimeRangeStandard");
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(5, "kIsGoodITSLayersAll");
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(6, "kNoCollInRofStandard");
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(7, "kNoHighMultCollInPrevRof");
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(8, "occupancy");
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(9, "MultCorrelation");
    registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(10, "cfgEvSelV0AT0ACut");
    if (cfgUseTentativeEventCounter) {
      registry.add("hEventCountTentative", "Number of Event;; Count", {HistType::kTH1D, {{10, 0, 10}}});
      registry.get<TH1>(HIST("hEventCountTentative"))->GetXaxis()->SetBinLabel(1, "after sel8");
      registry.get<TH1>(HIST("hEventCountTentative"))->GetXaxis()->SetBinLabel(2, "kNoSameBunchPileup");
      registry.get<TH1>(HIST("hEventCountTentative"))->GetXaxis()->SetBinLabel(3, "kIsGoodZvtxFT0vsPV");
      registry.get<TH1>(HIST("hEventCountTentative"))->GetXaxis()->SetBinLabel(4, "kNoCollInTimeRangeStandard");
      registry.get<TH1>(HIST("hEventCountTentative"))->GetXaxis()->SetBinLabel(5, "kIsGoodITSLayersAll");
      registry.get<TH1>(HIST("hEventCountTentative"))->GetXaxis()->SetBinLabel(6, "kNoCollInRofStandard");
      registry.get<TH1>(HIST("hEventCountTentative"))->GetXaxis()->SetBinLabel(7, "kNoHighMultCollInPrevRof");
      registry.get<TH1>(HIST("hEventCountTentative"))->GetXaxis()->SetBinLabel(8, "occupancy");
      registry.get<TH1>(HIST("hEventCountTentative"))->GetXaxis()->SetBinLabel(9, "MultCorrelation");
      registry.get<TH1>(HIST("hEventCountTentative"))->GetXaxis()->SetBinLabel(10, "cfgEvSelV0AT0ACut");
    }
    registry.add("hVtxZ", "Vexter Z distribution", {HistType::kTH1D, {axisVertex}});
    registry.add("hMult", "Multiplicity distribution", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    std::string hCentTitle = "Centrality distribution, Estimator " + std::to_string(cfgCentEstimator);
    registry.add("hCent", hCentTitle.c_str(), {HistType::kTH1D, {{90, 0, 90}}});
    if (!cfgUseSmallMemory) {
      registry.add("BeforeSel8_globalTracks_centT0C", "before sel8;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {axisCentForQA, axisNch}});
      registry.add("BeforeCut_globalTracks_centT0C", "before cut;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {axisCentForQA, axisNch}});
      registry.add("BeforeCut_PVTracks_centT0C", "before cut;Centrality T0C;mulplicity PV tracks", {HistType::kTH2D, {axisCentForQA, axisNchPV}});
      registry.add("BeforeCut_globalTracks_PVTracks", "before cut;mulplicity PV tracks;mulplicity global tracks", {HistType::kTH2D, {axisNchPV, axisNch}});
      registry.add("BeforeCut_globalTracks_multT0A", "before cut;mulplicity T0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisNch}});
      registry.add("BeforeCut_globalTracks_multV0A", "before cut;mulplicity V0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisNch}});
      registry.add("BeforeCut_multV0A_multT0A", "before cut;mulplicity T0A;mulplicity V0A", {HistType::kTH2D, {axisT0A, axisT0A}});
      registry.add("BeforeCut_multT0C_centT0C", "before cut;Centrality T0C;mulplicity T0C", {HistType::kTH2D, {axisCentForQA, axisT0C}});
      registry.add("globalTracks_centT0C", "after cut;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {axisCentForQA, axisNch}});
      registry.add("PVTracks_centT0C", "after cut;Centrality T0C;mulplicity PV tracks", {HistType::kTH2D, {axisCentForQA, axisNchPV}});
      registry.add("globalTracks_PVTracks", "after cut;mulplicity PV tracks;mulplicity global tracks", {HistType::kTH2D, {axisNchPV, axisNch}});
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
    registry.add("hChi2prTPCcls", "#chi^{2}/cluster for the TPC track segment", {HistType::kTH1D, {{100, 0., 5.}}});
    registry.add("hChi2prITScls", "#chi^{2}/cluster for the ITS track", {HistType::kTH1D, {{100, 0., 50.}}});
    registry.add("hnTPCClu", "Number of found TPC clusters", {HistType::kTH1D, {{100, 40, 180}}});
    registry.add("hnITSClu", "Number of found ITS clusters", {HistType::kTH1D, {{100, 0, 20}}});
    registry.add("hnTPCCrossedRow", "Number of crossed TPC Rows", {HistType::kTH1D, {{100, 40, 180}}});
    registry.add("hDCAz", "DCAz after cuts; DCAz (cm); Pt", {HistType::kTH2D, {{200, -0.5, 0.5}, {200, 0, 5}}});
    registry.add("hDCAxy", "DCAxy after cuts; DCAxy (cm); Pt", {HistType::kTH2D, {{200, -0.5, 0.5}, {200, 0, 5}}});
    registry.add("hTrackCorrection2d", "Correlation table for number of tracks table; uncorrected track; corrected track", {HistType::kTH2D, {axisNch, axisNch}});

    o2::framework::AxisSpec axis = axisPt;
    int nPtBins = axis.binEdges.size() - 1;
    double* ptBins = &(axis.binEdges)[0];
    fPtAxis = new TAxis(nPtBins, ptBins);

    if (cfgOutputNUAWeights) {
      fWeights->setPtBins(nPtBins, ptBins);
      fWeights->init(true, false);
    }

    TList* groupNUAWeightlist = new TList();
    groupNUAWeightlist->SetOwner(true);
    fGroupNUAList.setObject(groupNUAWeightlist);

    if (cfgOutputGroupNUAWeights) {
      groupNUAWeightPtr.resize(cfgGroupSplitRunNumber.value.size() + 1);
      for (uint i = 0; i < cfgGroupSplitRunNumber.value.size() + 1; i++) {
        GFWWeights* groupweight = nullptr;
        if (i < cfgGroupSplitRunNumber.value.size())
          groupweight = new GFWWeights(Form("groupweight_%d", cfgGroupSplitRunNumber.value[i]));
        else
          groupweight = new GFWWeights(Form("groupweight_last"));

        groupweight->setPtBins(nPtBins, ptBins);
        groupweight->init(true, false);
        groupNUAWeightlist->Add(groupweight);
        std::shared_ptr<GFWWeights> sharePtrGroupWeight(groupweight);
        groupNUAWeightPtr[i] = sharePtrGroupWeight;
      }
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
    fGFW->CreateRegions();

    if (cfgUseAdditionalEventCut) {
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(3257.29, -121.848, 1.98492, -0.0172128, 6.47528e-05, 154.756, -1.86072, -0.0274713, 0.000633499, -3.37757e-06);
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(3257.29, -121.848, 1.98492, -0.0172128, 6.47528e-05, 154.756, -1.86072, -0.0274713, 0.000633499, -3.37757e-06);

      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutLow->SetParameters(1654.46, -47.2379, 0.449833, -0.0014125, 150.773, -3.67334, 0.0530503, -0.000614061, 3.15956e-06);
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutHigh->SetParameters(1654.46, -47.2379, 0.449833, -0.0014125, 150.773, -3.67334, 0.0530503, -0.000614061, 3.15956e-06);

      fT0AV0AMean = new TF1("fT0AV0AMean", "[0]+[1]*x", 0, 200000);
      fT0AV0AMean->SetParameters(-1601.0581, 9.417652e-01);
      fT0AV0ASigma = new TF1("fT0AV0ASigma", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 200000);
      fT0AV0ASigma->SetParameters(463.4144, 6.796509e-02, -9.097136e-07, 7.971088e-12, -2.600581e-17);
    }

    if (cfgTrkSelRun3ITSMatch) {
      myTrackSel = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSall7Layers, TrackSelection::GlobalTrackRun3DCAxyCut::Default);
    } else {
      myTrackSel = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, TrackSelection::GlobalTrackRun3DCAxyCut::Default);
    }
    myTrackSel.SetMinNClustersTPC(cfgCutTPCclu);
    myTrackSel.SetMinNClustersITS(cfgCutITSclu);
    if (cfgCutDCAxyppPass3Enabled)
      myTrackSel.SetMaxDcaXYPtDep([](float pt) { return 0.004f + 0.013f / pt; }); // Tuned on the LHC22f anchored MC LHC23d1d on primary pions. 7 Sigmas of the resolution
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

  void fillFC(const GFW::CorrConfig& corrconf, const double& cent, const double& rndm)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (std::fabs(val) < 1)
        fFC->FillProfile(corrconf.Head.c_str(), cent, val, dnx, rndm);
      return;
    }
    for (auto i = 1; i <= fPtAxis->GetNbins(); i++) {
      dnx = fGFW->Calculate(corrconf, i - 1, kTRUE).real();
      if (dnx == 0)
        continue;
      val = fGFW->Calculate(corrconf, i - 1, kFALSE).real() / dnx;
      if (std::fabs(val) < 1)
        fFC->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i), cent, val, dnx, rndm);
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
    if (cfgAcceptanceGroup.value.empty() == false) {
      mGroupAcceptanceList = ccdb->getForTimeStamp<TList>(cfgAcceptance, timestamp);
      if (mGroupAcceptanceList == nullptr) {
        LOGF(fatal, "Could not load grouped acceptance weights from %s", cfgAcceptanceGroup.value.c_str());
      }
      LOGF(info, "Loaded grouped acceptance weights from %s (%p)", cfgAcceptanceGroup.value.c_str(), (void*)mGroupAcceptanceList);
    }
    if (cfgEfficiency.value.empty() == false) {
      mEfficiency = ccdb->getForTimeStamp<TH1D>(cfgEfficiency, timestamp);
      if (mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEfficiency.value.c_str(), (void*)mEfficiency);
    }
    correctionsLoaded = true;
  }

  bool setCurrentParticleWeights(float& weight_nue, float& weight_nua, float phi, float eta, float pt, float vtxz, int groupNUAIndex = 0)
  {
    float eff = 1.;
    if (mEfficiency)
      eff = mEfficiency->GetBinContent(mEfficiency->FindBin(pt));
    else
      eff = 1.0;
    if (eff == 0)
      return false;
    weight_nue = 1. / eff;
    if (cfgAcceptanceGroupUse) {
      if (mGroupAcceptanceList && mGroupAcceptanceList->At(groupNUAIndex)) {
        weight_nua = reinterpret_cast<GFWWeights*>(mGroupAcceptanceList->At(groupNUAIndex))->getNUA(phi, eta, vtxz);
      } else {
        weight_nua = 1;
      }
      return true;
    }
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
    if (cfgEvSelkIsGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return 0;
    }
    if (cfgEvSelkIsGoodZvtxFT0vsPV)
      registry.fill(HIST("hEventCountSpecific"), 2.5);
    if (cfgEvSelkNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      // no collisions in specified time range
      return 0;
    }
    if (cfgEvSelkNoCollInTimeRangeStandard)
      registry.fill(HIST("hEventCountSpecific"), 3.5);
    if (cfgEvSelkIsGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      // from Jan 9 2025 AOT meeting
      // cut time intervals with dead ITS staves
      return 0;
    }
    if (cfgEvSelkIsGoodITSLayersAll)
      registry.fill(HIST("hEventCountSpecific"), 4.5);
    if (cfgEvSelkNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      // no other collisions in this Readout Frame with per-collision multiplicity above threshold
      return 0;
    }
    if (cfgEvSelkNoCollInRofStandard)
      registry.fill(HIST("hEventCountSpecific"), 5.5);
    if (cfgEvSelkNoHighMultCollInPrevRof && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      // veto an event if FT0C amplitude in previous ITS ROF is above threshold
      return 0;
    }
    if (cfgEvSelkNoHighMultCollInPrevRof)
      registry.fill(HIST("hEventCountSpecific"), 6.5);
    auto multNTracksPV = collision.multNTracksPV();
    auto occupancy = collision.trackOccupancyInTimeRange();
    if (cfgEvSelOccupancy && (occupancy < cfgCutOccupancyLow || occupancy > cfgCutOccupancyHigh))
      return 0;
    if (cfgEvSelOccupancy)
      registry.fill(HIST("hEventCountSpecific"), 7.5);

    if (cfgEvSelMultCorrelation) {
      if (multNTracksPV < fMultPVCutLow->Eval(centrality))
        return 0;
      if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
        return 0;
      if (multTrk < fMultCutLow->Eval(centrality))
        return 0;
      if (multTrk > fMultCutHigh->Eval(centrality))
        return 0;
    }
    if (cfgEvSelMultCorrelation)
      registry.fill(HIST("hEventCountSpecific"), 8.5);

    // V0A T0A 5 sigma cut
    if (cfgEvSelV0AT0ACut && (std::fabs(collision.multFV0A() - fT0AV0AMean->Eval(collision.multFT0A())) > 5 * fT0AV0ASigma->Eval(collision.multFT0A())))
      return 0;
    if (cfgEvSelV0AT0ACut)
      registry.fill(HIST("hEventCountSpecific"), 9.5);

    return 1;
  }

  template <typename TCollision>
  void eventCounterQA(TCollision collision, const int multTrk, const float centrality)
  {
    registry.fill(HIST("hEventCountTentative"), 0.5);
    // Regradless of the event selection, fill the event counter histograms
    if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      registry.fill(HIST("hEventCountTentative"), 1.5);
    if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      registry.fill(HIST("hEventCountTentative"), 2.5);
    if (collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))
      registry.fill(HIST("hEventCountTentative"), 3.5);
    if (collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll))
      registry.fill(HIST("hEventCountTentative"), 4.5);
    if (collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard))
      registry.fill(HIST("hEventCountTentative"), 5.5);
    if (collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof))
      registry.fill(HIST("hEventCountTentative"), 6.5);
    auto multNTracksPV = collision.multNTracksPV();
    auto occupancy = collision.trackOccupancyInTimeRange();
    if (!(occupancy < cfgCutOccupancyLow || occupancy > cfgCutOccupancyHigh))
      registry.fill(HIST("hEventCountTentative"), 7.5);
    if (!((multNTracksPV < fMultPVCutLow->Eval(centrality)) || (multNTracksPV > fMultPVCutHigh->Eval(centrality)) || (multTrk < fMultCutLow->Eval(centrality)) || (multTrk > fMultCutHigh->Eval(centrality))))
      registry.fill(HIST("hEventCountTentative"), 8.5);
    if (!(std::fabs(collision.multFV0A() - fT0AV0AMean->Eval(collision.multFT0A())) > 5 * fT0AV0ASigma->Eval(collision.multFT0A())))
      registry.fill(HIST("hEventCountTentative"), 9.5);
  }

  template <typename TTrack>
  bool trackSelected(TTrack track)
  {
    if (cfgCutDCAzPtDepEnabled && (std::fabs(track.dcaZ()) > (0.004f + 0.013f / track.pt())))
      return false;

    if (cfgTrkSelSwitch) {
      return myTrackSel.IsSelected(track);
    } else {
      return ((track.tpcNClsFound() >= cfgCutTPCclu) && (track.itsNCls() >= cfgCutITSclu));
    }
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

  void process(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracks const& tracks)
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
    if (cfgUseTentativeEventCounter)
      eventCounterQA(collision, tracks.size(), cent);
    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), cent))
      return;
    registry.fill(HIST("hEventCount"), 3.5);
    float lRandom = fRndm->Rndm();
    float vtxz = collision.posZ();
    registry.fill(HIST("hVtxZ"), vtxz);
    registry.fill(HIST("hMult"), tracks.size());
    registry.fill(HIST("hCent"), cent);
    fGFW->Clear();
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
    double nTracksCorrected = 0;
    float independent = cent;
    if (cfgUseNch)
      independent = static_cast<float>(tracks.size());

    int groupNUAIndex = 0;
    if (cfgOutputGroupNUAWeights || cfgAcceptanceGroupUse) {
      for (uint i = 0; i < cfgGroupSplitRunNumber.value.size(); i++) {
        if (currentRunNumber < cfgGroupSplitRunNumber.value.at(i)) {
          break;
        } else {
          groupNUAIndex++;
        }
      }
    }

    for (const auto& track : tracks) {
      if (!trackSelected(track))
        continue;
      bool withinPtPOI = (cfgCutPtPOIMin < track.pt()) && (track.pt() < cfgCutPtPOIMax); // within POI pT range
      bool withinPtRef = (cfgCutPtRefMin < track.pt()) && (track.pt() < cfgCutPtRefMax); // within RF pT range
      if (cfgOutputNUAWeights) {
        if (cfgOutputNUAWeightsRefPt) {
          if (withinPtRef)
            fWeights->fill(track.phi(), track.eta(), vtxz, track.pt(), cent, 0);
        } else {
          fWeights->fill(track.phi(), track.eta(), vtxz, track.pt(), cent, 0);
        }
      }
      if (cfgOutputGroupNUAWeights) {
        if (cfgOutputNUAWeightsRefPt) {
          if (withinPtRef) {
            groupNUAWeightPtr[groupNUAIndex]->fill(track.phi(), track.eta(), vtxz, track.pt(), cent, 0);
          }
        } else {
          groupNUAWeightPtr[groupNUAIndex]->fill(track.phi(), track.eta(), vtxz, track.pt(), cent, 0);
        }
      }
      if (!setCurrentParticleWeights(weff, wacc, track.phi(), track.eta(), track.pt(), vtxz, groupNUAIndex))
        continue;
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
        nTracksCorrected += weff;
      }
      if (withinPtRef)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), wacc * weff, 1);
      if (withinPtPOI)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), wacc * weff, 2);
      if (withinPtPOI && withinPtRef)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), wacc * weff, 4);
    }
    registry.fill(HIST("hTrackCorrection2d"), tracks.size(), nTracksCorrected);

    // Filling Flow Container
    for (uint l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
      fillFC(corrconfigs.at(l_ind), independent, lRandom);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowTask>(cfgc)};
}
