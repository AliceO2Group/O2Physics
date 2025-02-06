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
//

/// \file   flowRunbyRun.cxx
/// \author Zhiyong Lu (zhiyong.lu@cern.ch)
/// \since  Oct/30/2024
/// \brief  jira: PWGCF-254, produce Run-by-Run QA plots and flow analysis for Run3

#include <CCDB/BasicCCDBManager.h>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <map>
#include <string>
#include <memory>
#include <utility>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/CCDB/ctpRateFetcher.h"

#include "GFWPowerArray.h"
#include "GFW.h"
#include "GFWCumulant.h"
#include "GFWWeights.h"
#include "GFWWeightsList.h"
#include "FlowContainer.h"
#include "TList.h"
#include <TProfile.h>
#include <TRandom3.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowRunbyRun {

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtRefMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtRefMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for all tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 10.0f, "Maximal pT for all tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "max DCA to vertex z")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Use Nch for flow observables")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 30, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeights, bool, false, "NUA weights are filled in ref pt bins")
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeightsRefPt, bool, false, "NUA weights are filled in ref pt bins")
  O2_DEFINE_CONFIGURABLE(cfgDynamicRunNumber, bool, false, "Add runNumber during runtime")
  O2_DEFINE_CONFIGURABLE(cfgGetInteractionRate, bool, false, "Get interaction rate from CCDB")
  O2_DEFINE_CONFIGURABLE(cfgUseInteractionRateCut, bool, false, "Use events with low interaction rate")
  O2_DEFINE_CONFIGURABLE(cfgCutMaxIR, float, 50.0f, "maximum interaction rate (kHz)")
  O2_DEFINE_CONFIGURABLE(cfgCutMinIR, float, 0.0f, "minimum interaction rate (kHz)")
  Configurable<std::vector<int>> cfgRunNumbers{"cfgRunNumbers", std::vector<int>{544095, 544098, 544116, 544121, 544122, 544123, 544124}, "Preconfigured run numbers"};
  Configurable<std::vector<std::string>> cfgUserDefineGFWCorr{"cfgUserDefineGFWCorr", std::vector<std::string>{"refN10 {2} refP10 {-2}"}, "User defined GFW CorrelatorConfig"};
  Configurable<std::vector<std::string>> cfgUserDefineGFWName{"cfgUserDefineGFWName", std::vector<std::string>{"Ch10Gap22"}, "User defined GFW Name"};

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10}, "pt axis for histograms"};
  ConfigurableAxis axisIndependent{"axisIndependent", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "X axis for histograms"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  // Define output
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  OutputObj<GFWWeightsList> fGFWWeightsList{GFWWeightsList("GFWWeightsList")};
  HistogramRegistry registry{"registry"};

  // define global variables
  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  std::vector<GFW::CorrConfig> corrconfigsFC;
  TAxis* fPtAxis;
  TRandom3* fRndm = new TRandom3(0);
  int lastRunNumer = -1;
  std::vector<int> runNumbers;                                        // vector of run numbers
  std::map<int, std::vector<std::shared_ptr<TH1>>> th1sList;          // map of histograms for all runs
  std::map<int, std::vector<std::shared_ptr<TProfile>>> profilesList; // map of profiles for all runs
  enum OutputTH1Names {
    // here are TProfiles for vn-pt correlations that are not implemented in GFW
    hPhi = 0,
    hEta,
    hVtxZ,
    hMult,
    hCent,
    kCount_TH1Names
  };
  enum OutputTProfileNames {
    c22 = 0,
    c22_gap10,
    kCount_TProfileNames
  };
  int mRunNumber{-1};
  uint64_t mSOR{0};
  double mMinSeconds{-1.};
  std::unordered_map<int, TH2*> gHadronicRate;
  ctpRateFetcher mRateFetcher;
  TH2* gCurrentHadronicRate;

  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults>>;
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA>>;

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);

    fGFWWeightsList->init("weightList");

    // Add output histograms to the registry
    runNumbers = cfgRunNumbers;
    for (const auto& runNumber : runNumbers) {
      createOutputObjectsForRun(runNumber);
    }

    o2::framework::AxisSpec axis = axisPt;
    int nPtBins = axis.binEdges.size() - 1;
    double* ptBins = &(axis.binEdges)[0];
    fPtAxis = new TAxis(nPtBins, ptBins);

    // Create FlowContainer
    TObjArray* oba = new TObjArray();
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

    fGFW->AddRegion("full", -0.8, 0.8, 1, 1);
    fGFW->AddRegion("refN10", -0.8, -0.5, 1, 1);
    fGFW->AddRegion("refP10", 0.5, 0.8, 1, 1);
    corrconfigs.resize(kCount_TProfileNames);
    corrconfigs[c22] = fGFW->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE);
    corrconfigs[c22_gap10] = fGFW->GetCorrelatorConfig("refN10 {2} refP10 {-2}", "Ch10Gap22", kFALSE);
    if (!userDefineGFWCorr.empty() && !userDefineGFWName.empty()) {
      LOGF(info, "User adding GFW CorrelatorConfig:");
      // attentaion: here we follow the index of cfgUserDefineGFWCorr
      for (uint i = 0; i < userDefineGFWCorr.size(); i++) {
        if (i >= userDefineGFWName.size()) {
          LOGF(fatal, "The names you provided are more than configurations. userDefineGFWName.size(): %d > userDefineGFWCorr.size(): %d", userDefineGFWName.size(), userDefineGFWCorr.size());
          break;
        }
        LOGF(info, "%d: %s %s", i, userDefineGFWCorr.at(i).c_str(), userDefineGFWName.at(i).c_str());
        corrconfigsFC.push_back(fGFW->GetCorrelatorConfig(userDefineGFWCorr.at(i).c_str(), userDefineGFWName.at(i).c_str(), kFALSE));
      }
    }
    fGFW->CreateRegions();
  }

  template <char... chars>
  void fillProfile(const GFW::CorrConfig& corrconf, std::shared_ptr<TProfile> profile, const double& cent)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (std::fabs(val) < 1)
        profile->Fill(cent, val, dnx);
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

  void createOutputObjectsForRun(int runNumber)
  {
    std::vector<std::shared_ptr<TH1>> histos(kCount_TH1Names);
    histos[hPhi] = registry.add<TH1>(Form("%d/hPhi", runNumber), "", {HistType::kTH1D, {axisPhi}});
    histos[hEta] = registry.add<TH1>(Form("%d/hEta", runNumber), "", {HistType::kTH1D, {axisEta}});
    histos[hVtxZ] = registry.add<TH1>(Form("%d/hVtxZ", runNumber), "", {HistType::kTH1D, {axisVertex}});
    histos[hMult] = registry.add<TH1>(Form("%d/hMult", runNumber), "", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    histos[hCent] = registry.add<TH1>(Form("%d/hCent", runNumber), "", {HistType::kTH1D, {{90, 0, 90}}});
    th1sList.insert(std::make_pair(runNumber, histos));

    std::vector<std::shared_ptr<TProfile>> profiles(kCount_TProfileNames);
    profiles[c22] = registry.add<TProfile>(Form("%d/c22", runNumber), "", {HistType::kTProfile, {axisIndependent}});
    profiles[c22_gap10] = registry.add<TProfile>(Form("%d/c22_gap10", runNumber), "", {HistType::kTProfile, {axisIndependent}});
    profilesList.insert(std::make_pair(runNumber, profiles));

    if (cfgOutputNUAWeights) {
      // weightsList
      o2::framework::AxisSpec axis = axisPt;
      int nPtBins = axis.binEdges.size() - 1;
      double* ptBins = &(axis.binEdges)[0];
      fGFWWeightsList->addGFWWeightsByRun(runNumber, nPtBins, ptBins, true, false);
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
    if (!collision.sel8())
      return;
    if (tracks.size() < 1)
      return;
    // detect run number
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (cfgGetInteractionRate) {
      initHadronicRate(bc);
      double hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), mRunNumber, "ZNC hadronic") * 1.e-3; //
      double seconds = bc.timestamp() * 1.e-3 - mMinSeconds;
      if (cfgUseInteractionRateCut && (hadronicRate < cfgCutMinIR || hadronicRate > cfgCutMaxIR)) // cut on hadronic rate
        return;
      gCurrentHadronicRate->Fill(seconds, hadronicRate);
    }
    int runNumber = bc.runNumber();
    float lRandom = fRndm->Rndm();
    if (runNumber != lastRunNumer) {
      lastRunNumer = runNumber;
      if (cfgDynamicRunNumber && std::find(runNumbers.begin(), runNumbers.end(), runNumber) == runNumbers.end()) {
        // if run number is not in the preconfigured list, create new output histograms for this run
        createOutputObjectsForRun(runNumber);
        runNumbers.push_back(runNumber);
      }

      if (th1sList.find(runNumber) == th1sList.end()) {
        LOGF(fatal, "RunNumber %d not found in th1sList", runNumber);
        return;
      }
    }

    th1sList[runNumber][hVtxZ]->Fill(collision.posZ());
    th1sList[runNumber][hMult]->Fill(tracks.size());
    th1sList[runNumber][hCent]->Fill(collision.centFT0C());

    fGFW->Clear();
    const auto cent = collision.centFT0C();
    float weff = 1, wacc = 1;
    for (const auto& track : tracks) {
      th1sList[runNumber][hPhi]->Fill(track.phi());
      th1sList[runNumber][hEta]->Fill(track.eta());
      // bool WithinPtPOI = (cfgCutPtPOIMin < track.pt()) && (track.pt() < cfgCutPtPOIMax); // within POI pT range
      bool withinPtRef = (cfgCutPtRefMin < track.pt()) && (track.pt() < cfgCutPtRefMax); // within RF pT range
      if (withinPtRef) {
        fGFW->Fill(track.eta(), 1, track.phi(), wacc * weff, 1);
      }
      if (cfgOutputNUAWeights) {
        if (cfgOutputNUAWeightsRefPt) {
          if (withinPtRef) {
            GFWWeights* weight = fGFWWeightsList->getGFWWeightsByRun(runNumber);
            if (!weight) {
              LOGF(fatal, "Could not find the weight for run %d", runNumber);
              return;
            }
            weight->fill(track.phi(), track.eta(), collision.posZ(), track.pt(), cent, 0);
          }
        } else {
          GFWWeights* weight = fGFWWeightsList->getGFWWeightsByRun(runNumber);
          if (!weight) {
            LOGF(fatal, "Could not find the weight for run %d", runNumber);
            return;
          }
          weight->fill(track.phi(), track.eta(), collision.posZ(), track.pt(), cent, 0);
        }
      }
    }

    // Filling TProfile
    for (uint i = 0; i < kCount_TProfileNames; ++i) {
      fillProfile(corrconfigs[i], profilesList[runNumber][i], cent);
    }
    // Filling Flow Container
    for (uint l_ind = 0; l_ind < corrconfigsFC.size(); l_ind++) {
      fillFC(corrconfigsFC.at(l_ind), cent, lRandom);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowRunbyRun>(cfgc)};
}
