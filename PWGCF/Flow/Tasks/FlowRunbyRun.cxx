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
// code author: Zhiyong Lu (zhiyong.lu@cern.ch)
// jira: PWGCF-254
// Produce Run-by-Run QA plots and flow analysis for Run3

#include <CCDB/BasicCCDBManager.h>
#include <cmath>
#include <vector>
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

#include "GFWPowerArray.h"
#include "GFW.h"
#include "GFWCumulant.h"
#include "GFWWeights.h"
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
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  // Define output
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  HistogramRegistry registry{"registry"};

  // define global variables
  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  std::vector<GFW::CorrConfig> corrconfigsFC;
  TAxis* fPtAxis;
  TRandom3* fRndm = new TRandom3(0);
  std::vector<int> RunNumbers;                                        // vector of run numbers
  std::map<int, std::vector<std::shared_ptr<TH1>>> TH1sList;          // map of histograms for all runs
  std::map<int, std::vector<std::shared_ptr<TProfile>>> ProfilesList; // map of profiles for all runs
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

  using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults>>;
  using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA>>;

  void init(InitContext const&)
  {
    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(nolaterthan.value);

    // Add output histograms to the registry
    RunNumbers = cfgRunNumbers;
    for (auto& runNumber : RunNumbers) {
      CreateOutputObjectsForRun(runNumber);
    }

    o2::framework::AxisSpec axis = axisPt;
    int nPtBins = axis.binEdges.size() - 1;
    double* PtBins = &(axis.binEdges)[0];
    fPtAxis = new TAxis(nPtBins, PtBins);

    // Create FlowContainer
    TObjArray* oba = new TObjArray();
    std::vector<std::string> UserDefineGFWCorr = cfgUserDefineGFWCorr;
    std::vector<std::string> UserDefineGFWName = cfgUserDefineGFWName;
    if (!UserDefineGFWCorr.empty() && !UserDefineGFWName.empty()) {
      for (uint i = 0; i < UserDefineGFWName.size(); i++) {
        oba->Add(new TNamed(UserDefineGFWName.at(i).c_str(), UserDefineGFWName.at(i).c_str()));
      }
    }
    fFC->SetName("FlowContainer");
    fFC->SetXAxis(fPtAxis);
    fFC->Initialize(oba, axisIndependent, 1);
    delete oba;

    fGFW->AddRegion("full", -0.8, 0.8, 1, 1);
    fGFW->AddRegion("refN10", -0.8, -0.5, 1, 1);
    fGFW->AddRegion("refP10", 0.5, 0.8, 1, 1);
    corrconfigs.resize(kCount_TProfileNames);
    corrconfigs[c22] = fGFW->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE);
    corrconfigs[c22_gap10] = fGFW->GetCorrelatorConfig("refN10 {2} refP10 {-2}", "Ch10Gap22", kFALSE);
    if (!UserDefineGFWCorr.empty() && !UserDefineGFWName.empty()) {
      LOGF(info, "User adding GFW CorrelatorConfig:");
      // attentaion: here we follow the index of cfgUserDefineGFWCorr
      for (uint i = 0; i < UserDefineGFWCorr.size(); i++) {
        if (i >= UserDefineGFWName.size()) {
          LOGF(fatal, "The names you provided are more than configurations. UserDefineGFWName.size(): %d > UserDefineGFWCorr.size(): %d", UserDefineGFWName.size(), UserDefineGFWCorr.size());
          break;
        }
        LOGF(info, "%d: %s %s", i, UserDefineGFWCorr.at(i).c_str(), UserDefineGFWName.at(i).c_str());
        corrconfigsFC.push_back(fGFW->GetCorrelatorConfig(UserDefineGFWCorr.at(i).c_str(), UserDefineGFWName.at(i).c_str(), kFALSE));
      }
    }
    fGFW->CreateRegions();
  }

  template <char... chars>
  void FillProfile(const GFW::CorrConfig& corrconf, std::shared_ptr<TProfile> profile, const double& cent)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1)
        profile->Fill(cent, val, dnx);
      return;
    }
    return;
  }

  void FillFC(const GFW::CorrConfig& corrconf, const double& cent, const double& rndm)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1)
        fFC->FillProfile(corrconf.Head.c_str(), cent, val, dnx, rndm);
      return;
    }
    for (Int_t i = 1; i <= fPtAxis->GetNbins(); i++) {
      dnx = fGFW->Calculate(corrconf, i - 1, kTRUE).real();
      if (dnx == 0)
        continue;
      val = fGFW->Calculate(corrconf, i - 1, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1)
        fFC->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i), cent, val, dnx, rndm);
    }
    return;
  }

  void CreateOutputObjectsForRun(int runNumber)
  {
    std::vector<std::shared_ptr<TH1>> histos(kCount_TH1Names);
    histos[hPhi] = registry.add<TH1>(Form("%d/hPhi", runNumber), "", {HistType::kTH1D, {axisPhi}});
    histos[hEta] = registry.add<TH1>(Form("%d/hEta", runNumber), "", {HistType::kTH1D, {axisEta}});
    histos[hVtxZ] = registry.add<TH1>(Form("%d/hVtxZ", runNumber), "", {HistType::kTH1D, {axisVertex}});
    histos[hMult] = registry.add<TH1>(Form("%d/hMult", runNumber), "", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    histos[hCent] = registry.add<TH1>(Form("%d/hCent", runNumber), "", {HistType::kTH1D, {{90, 0, 90}}});
    TH1sList.insert(std::make_pair(runNumber, histos));

    std::vector<std::shared_ptr<TProfile>> profiles(kCount_TProfileNames);
    profiles[c22] = registry.add<TProfile>(Form("%d/c22", runNumber), "", {HistType::kTProfile, {axisIndependent}});
    profiles[c22_gap10] = registry.add<TProfile>(Form("%d/c22_gap10", runNumber), "", {HistType::kTProfile, {axisIndependent}});
    ProfilesList.insert(std::make_pair(runNumber, profiles));
  }

  void process(aodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, aodTracks const& tracks)
  {
    if (!collision.sel8())
      return;
    if (tracks.size() < 1)
      return;
    // detect run number
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int runNumber = bc.runNumber();
    float l_Random = fRndm->Rndm();
    if (std::find(RunNumbers.begin(), RunNumbers.end(), runNumber) == RunNumbers.end()) {
      // if run number is not in the preconfigured list, create new output histograms for this run
      CreateOutputObjectsForRun(runNumber);
      RunNumbers.push_back(runNumber);
    }

    if (TH1sList.find(runNumber) == TH1sList.end()) {
      LOGF(fatal, "RunNumber %d not found in TH1sList", runNumber);
      return;
    }

    TH1sList[runNumber][hVtxZ]->Fill(collision.posZ());
    TH1sList[runNumber][hMult]->Fill(tracks.size());
    TH1sList[runNumber][hCent]->Fill(collision.centFT0C());

    fGFW->Clear();
    const auto cent = collision.centFT0C();
    float weff = 1, wacc = 1;
    for (auto& track : tracks) {
      TH1sList[runNumber][hPhi]->Fill(track.phi());
      TH1sList[runNumber][hEta]->Fill(track.eta());
      bool WithinPtRef = (cfgCutPtRefMin < track.pt()) && (track.pt() < cfgCutPtRefMax); // within RF pT range
      if (WithinPtRef) {
        fGFW->Fill(track.eta(), 1, track.phi(), wacc * weff, 1);
      }
    }

    // Filling TProfile
    for (uint i = 0; i < kCount_TProfileNames; ++i) {
      FillProfile(corrconfigs[i], ProfilesList[runNumber][i], cent);
    }
    // Filling Flow Container
    for (uint l_ind = 0; l_ind < corrconfigsFC.size(); l_ind++) {
      FillFC(corrconfigsFC.at(l_ind), cent, l_Random);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowRunbyRun>(cfgc)};
}
