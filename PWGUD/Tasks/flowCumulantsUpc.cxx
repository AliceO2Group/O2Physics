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

/// \file   flowCumulantsUpc.cxx
/// \author Yongxi Du (yongxi.du@cern.ch), Mingrui Zhao (mingrui.zhao@mail.labz0.org, mingrui.zhao@cern.ch), Zhiyong Lu (zhiyong.lu@cern.ch)
/// \since  Mar/2025
/// \brief  jira: PWGUD-45, task to measure flow observables with cumulant method

#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/StringHelpers.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PxPyPzE4D.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TString.h>

#include <sys/types.h>

#include <RtypesCore.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowCumulantsUpc {

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgIfVertex, bool, false, "choose vertex or not")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.1f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtRefMin, float, 0.1f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtRefMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "max Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCCrossedRows, float, 70.0f, "minimum number of crossed TPC Rows")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCclu, float, 50.0f, "minimum number of found TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutITSclu, float, 5.0f, "minimum number of ITS clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCChi2NCl, int, 4, "max chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 30, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeights, bool, false, "Fill and output NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeightsRefPt, bool, false, "NUA weights are filled in ref pt bins")
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeightsRunbyRun, bool, false, "NUA weights are filled run-by-run")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgEvSelOccupancy, bool, false, "Occupancy cut")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyHigh, int, 1000, "High cut on TPC occupancy")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyLow, int, 0, "Low cut on TPC occupancy")
  O2_DEFINE_CONFIGURABLE(cfgGapSideA, bool, true, "only pass gapside A, otherwise pass A+C")
  O2_DEFINE_CONFIGURABLE(cfgGapSideC, bool, false, "only pass gapside C, otherwise pass A+C")
  O2_DEFINE_CONFIGURABLE(cfgGlobalTrack, bool, true, "require TPC+ITS track")
  O2_DEFINE_CONFIGURABLE(cfgUseNchCorrected, bool, true, "use corrected Nch for X axis")
  O2_DEFINE_CONFIGURABLE(cfgDcaxy, bool, true, "choose dcaxy")
  O2_DEFINE_CONFIGURABLE(cfgDcaz, bool, false, "choose dcaz")
  O2_DEFINE_CONFIGURABLE(cfgDcazCut, float, 10.0, "dcaz cut")
  O2_DEFINE_CONFIGURABLE(cfgConsistentEventFlag, int, 0, "Flag to select consistent events - 0: off, 1: v2{2} gap calculable, 2: v2{4} full calculable, 4: v2{4} gap calculable, 8: v2{4} 3sub calculable")

  Configurable<std::vector<std::string>> cfgUserDefineGFWCorr{"cfgUserDefineGFWCorr", std::vector<std::string>{"refN02 {2} refP02 {-2}", "refN12 {2} refP12 {-2}"}, "User defined GFW CorrelatorConfig"};
  Configurable<std::vector<std::string>> cfgUserDefineGFWName{"cfgUserDefineGFWName", std::vector<std::string>{"Ch02Gap22", "Ch12Gap22"}, "User defined GFW Name"};
  Configurable<std::vector<float>> cfgConsistentEventVector{"cfgConsistentEventVector", std::vector<float>{-0.8, -0.5, -0.4, 0.4, 0.5, 0.8}, "eta regions: left(min,max), mid(min,max), right(min,max)"};
  struct AcceptedTracks {
    int nNeg;
    int nMid;
    int nPos;
    int nFull;
  };

  ConfigurableAxis axisPtHist{"axisPtHist", {100, 0., 10.}, "pt axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10}, "pt axis for histograms"};
  ConfigurableAxis axisIndependent{"axisIndependent", {VARIABLE_WIDTH, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60}, "X axis for histograms"};
  ConfigurableAxis axisNch{"axisNch", {300, 0, 300}, "N_{ch}"};
  ConfigurableAxis axisDCAz{"axisDCAz", {200, -2, 2}, "DCA_{z} (cm)"};
  ConfigurableAxis axisDCAxy{"axisDCAxy", {200, -1, 1}, "DCA_{xy} (cm)"};

  // Added UPC Cuts
  SGSelector sgSelector;
  Configurable<float> cfgCutFV0{"cfgCutFV0", 50., "FV0A threshold"};
  Configurable<float> cfgCutFT0A{"cfgCutFT0A", 150., "FT0A threshold"};
  Configurable<float> cfgCutFT0C{"cfgCutFT0C", 50., "FT0C threshold"};
  Configurable<float> cfgCutZDC{"cfgCutZDC", 10., "ZDC threshold"};

  // Corrections
  TH1D* mEfficiency = nullptr;
  GFWWeights* mAcceptance = nullptr;
  bool correctionsLoaded = false;

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  // Define output
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  OutputObj<GFWWeights> fWeights{GFWWeights("weights")};
  HistogramRegistry registry{"registry"};
  OutputObj<FlowContainer> fFCMc{FlowContainer("FlowContainerMC")};
  OutputObj<GFWWeights> fWeightsMc{GFWWeights("weightsMC")};

  // define global variables
  GFW* fGFW = new GFW();
  GFW* fGFWMC = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  std::vector<GFW::CorrConfig> corrconfigsmc;
  TAxis* fPtAxis;
  TRandom3* fRndm = new TRandom3(0);
  TRandom3* fRndmMc = new TRandom3(0);
  int lastRunNumber = -1;
  std::vector<int> runNumbers;
  std::map<int, std::shared_ptr<TH3>> th3sPerRun; // map of TH3 histograms for all runs

  using UdTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using UdTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced, aod::UDCollisionSelExtras>;

  void init(InitContext const&)
  {
    const AxisSpec axisVertex{40, -20, 20, "Vtxz (cm)"};
    const AxisSpec axisPhi{60, 0.0, constants::math::TwoPI, "#varphi"};
    const AxisSpec axisEta{40, -1., 1., "#eta"};

    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    // Add some output objects to the histogram registry
    // Event QA
    registry.add("hEventCount", "Number of Event;; Count", {HistType::kTH1D, {{6, 0, 6}}});
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(1, "Filtered event");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(2, "after gapside selection");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(3, "after vertex selection");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(4, "after occupancy");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(5, "after loadcorrection");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(6, "after consistency check");

    registry.add("hTrackCount", "Number of tracks;; Count", {HistType::kTH1D, {{9, 0, 9}}});
    registry.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(1, "default UD tracks");
    registry.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(2, "after PVContributor");
    registry.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(3, "after ITS+TPC");
    registry.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(4, "after TPC CrossRow");
    registry.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(5, "after TPC Clu");
    registry.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(6, "after TPC chi2");
    registry.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(7, "after ITS Clu");
    registry.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(8, "after DCAz");
    registry.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(9, "after DCAxy");
    registry.add("hVtxZ", "Vexter Z distribution", {HistType::kTH1D, {axisVertex}});
    registry.add("hMultWoSel", "Multiplicity distribution", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    registry.add("hMult", "Multiplicity distribution", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    // Reco Track QA
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
    // Mc track QA
    registry.add("hPhiMC", "#phi distribution", {HistType::kTH1D, {axisPhi}});
    registry.add("hPhiWeightedMC", "corrected #phi distribution", {HistType::kTH1D, {axisPhi}});
    registry.add("hEtaMC", "#eta distribution", {HistType::kTH1D, {axisEta}});
    registry.add("hPtMC", "p_{T} distribution before cut", {HistType::kTH1D, {axisPtHist}});
    registry.add("hPtRefMC", "p_{T} distribution after cut", {HistType::kTH1D, {axisPtHist}});
    registry.add("hTrackCorrection2dMC", "Correlation table for number of tracks table; uncorrected track; corrected track", {HistType::kTH2D, {axisNch, axisNch}});

    // // MC event QA histograms
    if (doprocessSim) {
      registry.add("eventCounterMC", "Number of MC Events;; Count", {HistType::kTH1D, {{5, 0, 5}}});
      registry.add("hVtxZMC", "Vexter Z distribution (MC)", {HistType::kTH1D, {axisVertex}});
      registry.add("hMultMC", "Multiplicity distribution (MC)", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
      registry.add("numberOfTracksMC", "Number of MC tracks;; Count", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    }

    o2::framework::AxisSpec axis = axisPt;
    int nPtBins = axis.binEdges.size() - 1;
    double* ptBins = &(axis.binEdges)[0];
    fPtAxis = new TAxis(nPtBins, ptBins);

    if (cfgOutputNUAWeights) {
      fWeights->setPtBins(nPtBins, ptBins);
      fWeights->init(true, false);
      fWeightsMc->setPtBins(nPtBins, ptBins);
      fWeightsMc->init(true, false);
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
    fFCMc->SetName("FlowContainerMC");
    fFCMc->SetXAxis(fPtAxis);
    fFCMc->Initialize(oba, axisIndependent, cfgNbootstrap);
    delete oba;

    // eta region
    float maxEta = std::fabs(cfgCutEta);
    float minEta = -1. * maxEta;
    fGFW->AddRegion("full", minEta, maxEta, 1, 1);
    fGFW->AddRegion("refN00", minEta, 0., 1, 1);   // gap0 negative region
    fGFW->AddRegion("refP00", 0., maxEta, 1, 1);   // gap0 positve region
    fGFW->AddRegion("refN02", minEta, -0.1, 1, 1); // gap2 negative region
    fGFW->AddRegion("refP02", 0.1, maxEta, 1, 1);  // gap2 positve region
    fGFW->AddRegion("refN04", minEta, -0.2, 1, 1); // gap4 negative region
    fGFW->AddRegion("refP04", 0.2, maxEta, 1, 1);  // gap4 positve region
    fGFW->AddRegion("refN06", minEta, -0.3, 1, 1); // gap6 negative region
    fGFW->AddRegion("refP06", 0.3, maxEta, 1, 1);  // gap6 positve region
    fGFW->AddRegion("refN08", minEta, -0.4, 1, 1);
    fGFW->AddRegion("refP08", 0.4, maxEta, 1, 1);
    fGFW->AddRegion("refN10", minEta, -0.5, 1, 1);
    fGFW->AddRegion("refP10", 0.5, maxEta, 1, 1);
    fGFW->AddRegion("refN12", minEta, -0.6, 1, 1);
    fGFW->AddRegion("refP12", 0.6, maxEta, 1, 1);
    fGFW->AddRegion("refN14", minEta, -0.7, 1, 1);
    fGFW->AddRegion("refP14", 0.7, maxEta, 1, 1);
    fGFW->AddRegion("refN", minEta, -0.4, 1, 1);
    fGFW->AddRegion("refP", 0.4, maxEta, 1, 1);
    fGFW->AddRegion("refM", -0.4, 0.4, 1, 1);
    fGFW->AddRegion("poiN", minEta, -0.4, 1 + fPtAxis->GetNbins(), 2);
    fGFW->AddRegion("poiN10", minEta, -0.5, 1 + fPtAxis->GetNbins(), 2);
    fGFW->AddRegion("poifull", minEta, maxEta, 1 + fPtAxis->GetNbins(), 2);
    fGFW->AddRegion("olN", minEta, -0.4, 1 + fPtAxis->GetNbins(), 4);
    fGFW->AddRegion("olN10", minEta, -0.5, 1 + fPtAxis->GetNbins(), 4);
    fGFW->AddRegion("olfull", minEta, maxEta, 1 + fPtAxis->GetNbins(), 4);

    // eta region for MC, can be different from data to study the effect of acceptance
    fGFWMC->AddRegion("full", minEta, maxEta, 1, 1);
    fGFWMC->AddRegion("refN00", minEta, 0., 1, 1);   // gap0 negative region
    fGFWMC->AddRegion("refP00", 0., maxEta, 1, 1);   // gap0 positve region
    fGFWMC->AddRegion("refN02", minEta, -0.1, 1, 1); // gap2 negative region
    fGFWMC->AddRegion("refP02", 0.1, maxEta, 1, 1);  // gap2 positve region
    fGFWMC->AddRegion("refN04", minEta, -0.2, 1, 1); // gap4 negative region
    fGFWMC->AddRegion("refP04", 0.2, maxEta, 1, 1);  // gap4 positve region
    fGFWMC->AddRegion("refN06", minEta, -0.3, 1, 1); // gap6 negative region
    fGFWMC->AddRegion("refP06", 0.3, maxEta, 1, 1);  // gap6 positve region
    fGFWMC->AddRegion("refN08", minEta, -0.4, 1, 1);
    fGFWMC->AddRegion("refP08", 0.4, maxEta, 1, 1);
    fGFWMC->AddRegion("refN10", minEta, -0.5, 1, 1);
    fGFWMC->AddRegion("refP10", 0.5, maxEta, 1, 1);
    fGFWMC->AddRegion("refN12", minEta, -0.6, 1, 1);
    fGFWMC->AddRegion("refP12", 0.6, maxEta, 1, 1);
    fGFWMC->AddRegion("refN14", minEta, -0.7, 1, 1);
    fGFWMC->AddRegion("refP14", 0.7, maxEta, 1, 1);
    fGFWMC->AddRegion("refN", minEta, -0.4, 1, 1);
    fGFWMC->AddRegion("refP", 0.4, maxEta, 1, 1);
    fGFWMC->AddRegion("refM", -0.4, 0.4, 1, 1);
    fGFWMC->AddRegion("poiN", minEta, -0.4, 1 + fPtAxis->GetNbins(), 2);
    fGFWMC->AddRegion("poiN10", minEta, -0.5, 1 + fPtAxis->GetNbins(), 2);
    fGFWMC->AddRegion("poifull", minEta, maxEta, 1 + fPtAxis->GetNbins(), 2);
    fGFWMC->AddRegion("olN", minEta, -0.4, 1 + fPtAxis->GetNbins(), 4);
    fGFWMC->AddRegion("olN10", minEta, -0.5, 1 + fPtAxis->GetNbins(), 4);
    fGFWMC->AddRegion("olfull", minEta, maxEta, 1 + fPtAxis->GetNbins(), 4);

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
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("full {3 -3}", "ChFull32", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("full {4 -4}", "ChFull42", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("full {2 2 -2 -2}", "ChFull24", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("full {2 2 2 -2 -2 -2}", "ChFull26", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN04 {2} refP04 {-2}", "Ch04Gap22", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN06 {2} refP06 {-2}", "Ch06Gap22", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "Ch08Gap22", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN10 {2} refP10 {-2}", "Ch10Gap22", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN12 {2} refP12 {-2}", "Ch12Gap22", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN04 {3} refP04 {-3}", "Ch04Gap32", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN06 {3} refP06 {-3}", "Ch06Gap32", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN08 {3} refP08 {-3}", "Ch08Gap32", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN10 {3} refP10 {-3}", "Ch10Gap32", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN12 {3} refP12 {-3}", "Ch12Gap32", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN04 {4} refP04 {-4}", "Ch04Gap42", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN06 {4} refP06 {-4}", "Ch06Gap42", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN08 {4} refP08 {-4}", "Ch08Gap42", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN10 {4} refP10 {-4}", "Ch10Gap42", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN12 {4} refP12 {-4}", "Ch12Gap42", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN {2} refP {-2}", "ChGap22", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("poifull full | olfull {2 -2}", "ChFull22", kTRUE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("poifull full | olfull {2 2 -2 -2}", "ChFull24", kTRUE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("poifull full | olfull {2 2 2 -2 -2 -2}", "ChFull26", kTRUE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("poiN10 refN10 | olN10 {2} refP10 {-2}", "Ch10Gap22", kTRUE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("poiN10 refN10 | olN10 {3} refP10 {-3}", "Ch10Gap32", kTRUE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("poiN10 refN10 | olN10 {4} refP10 {-4}", "Ch10Gap42", kTRUE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("full {4 -2 -2}", "ChFull422", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN04 {-2 -2} refP04 {4}", "Ch04GapA422", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN04 {4} refP04 {-2 -2}", "Ch04GapB422", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN10 {-2 -2} refP10 {4}", "Ch10GapA422", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN10 {4} refP10 {-2 -2}", "Ch10GapB422", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("full {3 2 -3 -2}", "ChFull3232", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("full {4 2 -4 -2}", "ChFull4242", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN04 {3 2} refP04 {-3 -2}", "Ch04Gap3232", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN04 {4 2} refP04 {-4 -2}", "Ch04Gap4242", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN04 {2 2} refP04 {-2 -2}", "Ch04Gap24", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN10 {3 2} refP10 {-3 -2}", "Ch10Gap3232", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN10 {4 2} refP10 {-4 -2}", "Ch10Gap4242", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("refN10 {2 2} refP10 {-2 -2}", "Ch10Gap24", kFALSE));
    corrconfigsmc.push_back(fGFWMC->GetCorrelatorConfig("poiN10 refN10 | olN10 {2 2} refP10 {-2 -2}", "Ch10Gap24", kTRUE));
    LOGF(info, "Embedded GFW CorrelatorConfig:");
    for (auto icorr = 0; icorr < static_cast<int>(corrconfigs.size()); icorr++) {
      LOGF(info, "corrconfigs.at(%d): %s", icorr, corrconfigs.at(icorr).Head.c_str());
    }
    if (!userDefineGFWCorr.empty() && !userDefineGFWName.empty()) {
      LOGF(info, "User adding GFW CorrelatorConfig:");
      // attentaion: here we follow the index of cfgUserDefineGFWCorr
      for (uint i = 0; i < userDefineGFWCorr.size(); i++) {
        if (i >= userDefineGFWName.size()) {
          LOGF(fatal, "The names you provided are more than configurations. userDefineGFWName.size(): %d > userDefineGFWCorr.size(): %d", userDefineGFWName.size(), userDefineGFWCorr.size());
          break;
        }
        if (userDefineGFWCorr.at(i).find("poi") != std::string::npos) {
          corrconfigs.push_back(fGFW->GetCorrelatorConfig(userDefineGFWCorr.at(i).c_str(), userDefineGFWName.at(i).c_str(), kTRUE));
          LOGF(info, "corrconfigs.at(%d): enable pt-Diff for %s %s", corrconfigs.size() - 1, userDefineGFWCorr.at(i).c_str(), userDefineGFWName.at(i).c_str());
        } else {
          corrconfigs.push_back(fGFW->GetCorrelatorConfig(userDefineGFWCorr.at(i).c_str(), userDefineGFWName.at(i).c_str(), kFALSE));
          LOGF(info, "corrconfigs.at(%d): %s %s", corrconfigs.size() - 1, userDefineGFWCorr.at(i).c_str(), userDefineGFWName.at(i).c_str());
        }
      }
    }
    fGFW->CreateRegions();
  }

  void createOutputObjectsForRun(int runNumber)
  {
    const AxisSpec axisPhi{60, 0.0, constants::math::TwoPI, "#varphi"};
    std::shared_ptr<TH3> histPhiEtaVtxz = registry.add<TH3>(Form("%d/hPhiEtaVtxz", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, {40, -10, 10}}});
    th3sPerRun.insert(std::make_pair(runNumber, histPhiEtaVtxz));
  }

  void fillFC(const GFW::CorrConfig& corrconf, const double& cent, const double& rndm)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (!corrconf.pTDif) {
      if (dnx == 0) {
        return;
      }
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (std::fabs(val) < 1) {
        fFC->FillProfile(corrconf.Head.c_str(), cent, val, dnx, rndm);
      }
      return;
    }
    for (auto i = 1; i <= fPtAxis->GetNbins(); i++) {
      dnx = fGFW->Calculate(corrconf, i - 1, kTRUE).real();
      if (dnx == 0) {
        continue;
      }
      val = fGFW->Calculate(corrconf, i - 1, kFALSE).real() / dnx;
      if (std::fabs(val) < 1) {
        fFC->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i), cent, val, dnx, rndm);
      }
    }
    return;
  }

  void fillFCMC(const GFW::CorrConfig& corrconf, const double& cent, const double& rndm)
  {
    double dnx, val;
    dnx = fGFWMC->Calculate(corrconf, 0, kTRUE).real();
    if (!corrconf.pTDif) {
      if (dnx == 0) {
        return;
      }
      val = fGFWMC->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (std::fabs(val) < 1) {
        fFCMc->FillProfile(corrconf.Head.c_str(), cent, val, dnx, rndm);
      }
      return;
    }
    for (auto i = 1; i <= fPtAxis->GetNbins(); i++) {
      dnx = fGFWMC->Calculate(corrconf, i - 1, kTRUE).real();
      if (dnx == 0) {
        continue;
      }
      val = fGFWMC->Calculate(corrconf, i - 1, kFALSE).real() / dnx;
      if (std::fabs(val) < 1) {
        fFCMc->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i), cent, val, dnx, rndm);
      }
    }
    return;
  }

  void loadCorrections(uint64_t timestamp)
  {
    if (correctionsLoaded) {
      return;
    }
    if (cfgAcceptance.value.empty() == false) {
      mAcceptance = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance, timestamp);
      if (mAcceptance) {
        LOGF(info, "Loaded acceptance weights from %s (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance);
      } else {
        LOGF(warning, "Could not load acceptance weights from %s (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance);
      }
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

  template <typename TTrack>
  bool trackSelected(TTrack track)
  {
    registry.fill(HIST("hTrackCount"), 0.5);
    // UPC selection
    if (!track.isPVContributor()) {
      return false;
    }
    registry.fill(HIST("hTrackCount"), 1.5);
    if (cfgGlobalTrack && !(track.hasITS() && track.hasTPC())) {
      return false;
    }
    registry.fill(HIST("hTrackCount"), 2.5);
    if (track.tpcNClsCrossedRows() < cfgCutTPCCrossedRows) {
      return false;
    }
    registry.fill(HIST("hTrackCount"), 3.5);
    auto tpcClu = track.tpcNClsFindable() - track.tpcNClsFindableMinusFound();
    if (tpcClu < cfgCutTPCclu) {
      return false;
    }
    registry.fill(HIST("hTrackCount"), 4.5);
    if (track.tpcChi2NCl() >= cfgCutTPCChi2NCl) {
      return false;
    }
    registry.fill(HIST("hTrackCount"), 5.5);
    if (track.itsNCls() < cfgCutITSclu) {
      return false;
    }
    registry.fill(HIST("hTrackCount"), 6.5);
    if (cfgDcaz && !(std::fabs(track.dcaZ()) < cfgDcazCut)) {
      return false;
    }
    registry.fill(HIST("hTrackCount"), 7.5);
    double dcaLimit = 0.0105 + 0.035 / std::pow(track.pt(), 1.1);
    if (cfgDcaxy && !(std::fabs(track.dcaXY()) < dcaLimit)) {
      return false;
    }
    registry.fill(HIST("hTrackCount"), 8.5);

    return true;
  }

  void processData(UDCollisionsFull::iterator const& collision, UdTracksFull const& tracks)
  {
    registry.fill(HIST("hEventCount"), 0.5);
    int gapSide = collision.gapSide();
    if (gapSide == 0) {
      if (!cfgGapSideA) {
        return;
      }
    }
    if (gapSide == 1) {
      if (!cfgGapSideC) {
        return;
      }
    }
    if (gapSide != 0 && gapSide != 1) {
      return;
    }
    int trueGapSide = sgSelector.trueGap(collision, cfgCutFV0, cfgCutFT0A, cfgCutFT0C, cfgCutZDC);
    gapSide = trueGapSide;
    if (gapSide == 0) {
      if (!cfgGapSideA) {
        return;
      }
    }
    if (gapSide == 1) {
      if (!cfgGapSideC) {
        return;
      }
    }
    if (gapSide != 0 && gapSide != 1) {
      return;
    }
    registry.fill(HIST("hEventCount"), 1.5);
    float lRandom = fRndm->Rndm();
    float vtxz = collision.posZ();
    registry.fill(HIST("hVtxZ"), vtxz);
    registry.fill(HIST("hMultWoSel"), tracks.size());
    fGFW->Clear();
    if (cfgIfVertex && std::abs(vtxz) > cfgCutVertex) {
      return;
    }
    registry.fill(HIST("hEventCount"), 2.5);
    int occupancy = collision.occupancyInTime();
    if (cfgEvSelOccupancy && (occupancy < cfgCutOccupancyLow || occupancy > cfgCutOccupancyHigh)) {
      return;
    }
    registry.fill(HIST("hEventCount"), 3.5);

    auto currentRunNumber = collision.runNumber();
    auto runDuration = ccdb->getRunDuration(currentRunNumber);
    loadCorrections(runDuration.first);
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
    registry.fill(HIST("hEventCount"), 4.5);

    // // track weights
    float weff = 1, wacc = 1;
    double nTracksRaw = 0.;
    double nTracksCorrected = 0.;
    AcceptedTracks acceptedTracks{0, 0, 0, 0};
    std::vector<float> consistentEventVector = cfgConsistentEventVector;

    for (const auto& track : tracks) {
      if (!trackSelected(track)) {
        continue;
      }
      auto momentum = std::array<double, 3>{track.px(), track.py(), track.pz()};
      double pt = RecoDecay::pt(momentum);
      double phi = RecoDecay::phi(momentum);
      double eta = RecoDecay::eta(momentum);
      bool withinPtPOI = (cfgCutPtPOIMin < pt) && (pt < cfgCutPtPOIMax); // within POI pT range
      bool withinPtRef = (cfgCutPtRefMin < pt) && (pt < cfgCutPtRefMax); // within RF pT range
      if (cfgOutputNUAWeights) {
        if (cfgOutputNUAWeightsRefPt) {
          if (withinPtRef) {
            fWeights->fill(phi, eta, vtxz, pt, 0., 0);
            if (cfgOutputNUAWeightsRunbyRun)
              th3sPerRun[currentRunNumber]->Fill(phi, eta, vtxz);
          }
        } else {
          fWeights->fill(phi, eta, vtxz, pt, 0., 0);
          if (cfgOutputNUAWeightsRunbyRun)
            th3sPerRun[currentRunNumber]->Fill(phi, eta, vtxz);
        }
      }
      if (!setCurrentParticleWeights(weff, wacc, phi, eta, pt, vtxz)) {
        continue;
      }
      registry.fill(HIST("hPt"), track.pt());

      if (cfgConsistentEventFlag && consistentEventVector.size() == 6) { // o2-linter: disable=magic-number (size match)
        acceptedTracks.nFull += 1;
        if (eta > consistentEventVector[0] && eta < consistentEventVector[1])
          acceptedTracks.nNeg += 1;
        if (eta > consistentEventVector[2] && eta < consistentEventVector[3])
          acceptedTracks.nMid += 1;
        if (eta > consistentEventVector[4] && eta < consistentEventVector[5])
          acceptedTracks.nPos += 1;
      }
      if (withinPtRef) {
        registry.fill(HIST("hChi2prTPCcls"), track.tpcChi2NCl());
        registry.fill(HIST("hChi2prITScls"), track.itsChi2NCl());
        auto tpcClu = track.tpcNClsFindable() - track.tpcNClsFindableMinusFound();
        registry.fill(HIST("hnTPCClu"), tpcClu);
        registry.fill(HIST("hnITSClu"), track.itsNCls());
        registry.fill(HIST("hnTPCCrossedRow"), track.tpcNClsCrossedRows());
        registry.fill(HIST("hPhi"), phi);
        registry.fill(HIST("hPhiWeighted"), phi, wacc);
        registry.fill(HIST("hEta"), eta);
        registry.fill(HIST("hPtRef"), pt);
        registry.fill(HIST("hDCAz"), track.dcaZ(), track.pt());
        registry.fill(HIST("hDCAxy"), track.dcaXY(), track.pt());
        nTracksRaw += 1.;
        nTracksCorrected += weff;
      }
      if (withinPtRef) {
        fGFW->Fill(eta, fPtAxis->FindBin(pt) - 1, phi, wacc * weff, 1);
      }
      if (withinPtPOI) {
        fGFW->Fill(eta, fPtAxis->FindBin(pt) - 1, phi, wacc * weff, 2);
      }
      if (withinPtPOI && withinPtRef) {
        fGFW->Fill(eta, fPtAxis->FindBin(pt) - 1, phi, wacc * weff, 4);
      }
    }
    registry.fill(HIST("hTrackCorrection2d"), nTracksRaw, nTracksCorrected);
    registry.fill(HIST("hMult"), nTracksRaw);
    float independent = nTracksRaw;
    if (cfgUseNchCorrected) {
      independent = nTracksCorrected;
    }

    if (cfgConsistentEventFlag) {
      if (cfgConsistentEventFlag & 1) {
        if (!acceptedTracks.nPos || !acceptedTracks.nNeg)
          return;
      } else if (cfgConsistentEventFlag & 2) {
        if (acceptedTracks.nFull < 4) // o2-linter: disable=magic-number (at least four tracks in full acceptance)
          return;
      } else if (cfgConsistentEventFlag & 4) {
        if (acceptedTracks.nPos < 2 || acceptedTracks.nNeg < 2) // o2-linter: disable=magic-number (at least two tracks in each subevent)
          return;
      }
      if (cfgConsistentEventFlag & 8) {
        if (acceptedTracks.nPos < 2 || acceptedTracks.nMid < 2 || acceptedTracks.nNeg < 2) // o2-linter: disable=magic-number (at least two tracks in all three subevents)
          return;
      }
    }
    registry.fill(HIST("hEventCount"), 5.5);

    // Filling Flow Container
    for (uint l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
      fillFC(corrconfigs.at(l_ind), independent, lRandom);
    }
  }
  PROCESS_SWITCH(FlowCumulantsUpc, processData, "processData", true);

  //-----------------------------------------------------------------------------------------------------------------------
  void processSim(aod::UDMcCollision const& mcCollision, aod::UDMcParticles const& mcParticles)
  {
    registry.fill(HIST("eventCounterMC"), 0.5);

    registry.fill(HIST("hEventCount"), 1.5);
    float vtxz = mcCollision.posZ();
    registry.fill(HIST("hVtxZMC"), vtxz);
    registry.fill(HIST("hMultMC"), mcParticles.size());

    auto massPion = o2::constants::physics::MassPionCharged;
    registry.fill(HIST("numberOfTracksMC"), mcParticles.size());
    // LOGF(info, "New event! mcParticles.size() = %d", mcParticles.size());

    float lRandomMc = fRndmMc->Rndm();
    fGFWMC->Clear();

    // // track weights
    float weff = 1, wacc = 1;
    double nTracksCorrected = 0;
    float independent = static_cast<float>(mcParticles.size());

    for (const auto& mcParticle : mcParticles) {
      if (!mcParticle.isPhysicalPrimary())
        continue;
      std::array<double, 3> momentum = {mcParticle.px(), mcParticle.py(), mcParticle.pz()};
      double energy = std::sqrt(momentum[0] * momentum[0] + momentum[1] * momentum[1] + momentum[2] * momentum[2] + massPion * massPion);
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> protoMC(momentum[0], momentum[1], momentum[2], energy);
      double kEtaCut = 0.8;
      double kPtCut = 0.1;
      if (!(std::fabs(protoMC.Eta()) < kEtaCut && protoMC.Pt() > kPtCut)) {
        continue;
      }
      // auto momentum = std::array<double, 3>{mcParticle.px(), mcParticle.py(), mcParticle.pz()};
      double pt = RecoDecay::pt(momentum);
      double phi = RecoDecay::phi(momentum);
      double eta = RecoDecay::eta(momentum);
      bool withinPtPOI = (cfgCutPtPOIMin < pt) && (pt < cfgCutPtPOIMax); // within POI pT range
      bool withinPtRef = (cfgCutPtRefMin < pt) && (pt < cfgCutPtRefMax); // within RF pT range
      if (cfgOutputNUAWeights) {
        if (cfgOutputNUAWeightsRefPt) {
          if (withinPtRef) {
            fWeightsMc->fill(phi, eta, vtxz, pt, 0., 0);
          }
        } else {
          fWeightsMc->fill(phi, eta, vtxz, pt, 0., 0);
        }
      }
      if (!setCurrentParticleWeights(weff, wacc, phi, eta, pt, vtxz)) {
        continue;
      }
      registry.fill(HIST("hPtMC"), pt);
      if (withinPtRef) {
        registry.fill(HIST("hPhiMC"), phi);
        registry.fill(HIST("hPhiWeightedMC"), phi, wacc);
        registry.fill(HIST("hEtaMC"), eta);
        registry.fill(HIST("hPtRefMC"), pt);
        nTracksCorrected += weff;
      }
      if (withinPtRef) {
        fGFWMC->Fill(eta, fPtAxis->FindBin(pt) - 1, phi, wacc * weff, 1);
      }
      if (withinPtPOI) {
        fGFWMC->Fill(eta, fPtAxis->FindBin(pt) - 1, phi, wacc * weff, 2);
      }
      if (withinPtPOI && withinPtRef) {
        fGFWMC->Fill(eta, fPtAxis->FindBin(pt) - 1, phi, wacc * weff, 4);
      }
    }
    registry.fill(HIST("hTrackCorrection2dMC"), mcParticles.size(), nTracksCorrected);

    // Filling Flow Container
    for (uint l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
      fillFCMC(corrconfigs.at(l_ind), independent, lRandomMc);
    }
  }
  PROCESS_SWITCH(FlowCumulantsUpc, processSim, "processSim", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowCumulantsUpc>(cfgc),
  };
}
