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

#include <CCDB/BasicCCDBManager.h>
#include <cmath>
#include <vector>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"

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

struct FlowPbPbTask {

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Use Nch for flow observables")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeights, bool, false, "Fill and output NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPtHist{"axisPtHist", {100, 0., 10.}, "pt axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10}, "pt axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "centrality axis for histograms"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);

  // Corrections
  TH1D* mEfficiency = nullptr;
  GFWWeights* mAcceptance = nullptr;
  bool correctionsLoaded = false;

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  // Define output
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  OutputObj<GFWWeights> fWeights{GFWWeights("weights")};
  HistogramRegistry registry{"registry"};

  // define global variables
  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  TAxis* fPtAxis;
  TRandom3* fRndm = new TRandom3(0);
  std::vector<std::vector<std::shared_ptr<TProfile>>> BootstrapArray;
  enum ExtraProfile {
    // here are TProfiles for vn-pt correlations that are not implemented in GFW
    kMeanPt_InGap08 = 0,
    kC22_Gap08_Weff,
    kC22_Gap08_MeanPt,
    kPtVarParA_InGap08,
    kPtVarParB_InGap08,
    // Count the total number of enum
    kCount_ExtraProfile
  };

  using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>;
  using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra>>;

  void init(InitContext const&)
  {
    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(nolaterthan.value);

    // Add some output objects to the histogram registry
    registry.add("hPhi", "", {HistType::kTH1D, {axisPhi}});
    registry.add("hEta", "", {HistType::kTH1D, {axisEta}});
    registry.add("hPt", "", {HistType::kTH1D, {axisPtHist}});
    registry.add("hVtxZ", "", {HistType::kTH1D, {axisVertex}});
    registry.add("hMult", "", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    registry.add("hCent", "", {HistType::kTH1D, {{90, 0, 90}}});
    registry.add("hMeanPt", "", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("hMeanPtWithinGap08", "", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c22_gap08_Weff", "", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c22_gap08_trackMeanPt", "", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("PtVariance_partA_WithinGap08", "", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("PtVariance_partB_WithinGap08", "", {HistType::kTProfile, {axisMultiplicity}});

    // initial array
    BootstrapArray.resize(cfgNbootstrap);
    for (int i = 0; i < cfgNbootstrap; i++) {
      BootstrapArray[i].resize(kCount_ExtraProfile);
    }
    for (int i = 0; i < cfgNbootstrap; i++) {
      BootstrapArray[i][kMeanPt_InGap08] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("BootstrapContainer_%d/hMeanPtWithinGap08", i), "", {HistType::kTProfile, {axisMultiplicity}}));
      BootstrapArray[i][kC22_Gap08_Weff] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("BootstrapContainer_%d/c22_gap08_Weff", i), "", {HistType::kTProfile, {axisMultiplicity}}));
      BootstrapArray[i][kC22_Gap08_MeanPt] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("BootstrapContainer_%d/c22_gap08_trackMeanPt", i), "", {HistType::kTProfile, {axisMultiplicity}}));
      BootstrapArray[i][kPtVarParA_InGap08] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("BootstrapContainer_%d/PtVariance_partA_WithinGap08", i), "", {HistType::kTProfile, {axisMultiplicity}}));
      BootstrapArray[i][kPtVarParB_InGap08] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("BootstrapContainer_%d/PtVariance_partB_WithinGap08", i), "", {HistType::kTProfile, {axisMultiplicity}}));
    }

    o2::framework::AxisSpec axis = axisPt;
    int nPtBins = axis.binEdges.size() - 1;
    double* PtBins = &(axis.binEdges)[0];
    fPtAxis = new TAxis(nPtBins, PtBins);

    if (cfgOutputNUAWeights) {
      fWeights->SetPtBins(nPtBins, PtBins);
      fWeights->Init(true, false);
    }

    // add in FlowContainer to Get boostrap sample automatically
    TObjArray* oba = new TObjArray();
    oba->Add(new TNamed("ChGap22", "ChGap22"));
    for (Int_t i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("ChGap22_pt_%i", i + 1), "ChGap22_pTDiff"));
    oba->Add(new TNamed("ChFull22", "ChFull22"));
    oba->Add(new TNamed("ChFull32", "ChFull32"));
    oba->Add(new TNamed("ChFull42", "ChFull42"));
    oba->Add(new TNamed("ChFull24", "ChFull24"));
    oba->Add(new TNamed("ChFull26", "ChFull26"));
    oba->Add(new TNamed("Ch04Gap22", "Ch04Gap22"));
    oba->Add(new TNamed("Ch06Gap22", "Ch06Gap22"));
    oba->Add(new TNamed("Ch08Gap22", "Ch08Gap22"));
    oba->Add(new TNamed("Ch10Gap22", "Ch10Gap22"));
    oba->Add(new TNamed("Ch04Gap32", "Ch04Gap32"));
    oba->Add(new TNamed("Ch06Gap32", "Ch06Gap32"));
    oba->Add(new TNamed("Ch08Gap32", "Ch08Gap32"));
    oba->Add(new TNamed("Ch10Gap32", "Ch10Gap32"));
    oba->Add(new TNamed("Ch04Gap42", "Ch04Gap42"));
    oba->Add(new TNamed("Ch06Gap42", "Ch06Gap42"));
    oba->Add(new TNamed("Ch08Gap42", "Ch08Gap42"));
    oba->Add(new TNamed("Ch10Gap42", "Ch10Gap42"));
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
    fFC->SetName("FlowContainer");
    fFC->SetXAxis(fPtAxis);
    fFC->Initialize(oba, axisMultiplicity, cfgNbootstrap);
    delete oba;

    // eta region
    fGFW->AddRegion("full", -0.8, 0.8, 1, 1);
    fGFW->AddRegion("refN04", -0.8, -0.2, 1, 1); // gap4 negative region
    fGFW->AddRegion("refP04", 0.2, 0.8, 1, 1);   // gap4 positve region
    fGFW->AddRegion("refN06", -0.8, -0.3, 1, 1); // gap6 negative region
    fGFW->AddRegion("refP06", 0.3, 0.8, 1, 1);   // gap6 positve region
    fGFW->AddRegion("refN08", -0.8, -0.4, 1, 1);
    fGFW->AddRegion("refP08", 0.4, 0.8, 1, 1);
    fGFW->AddRegion("refN10", -0.8, -0.5, 1, 1);
    fGFW->AddRegion("refP10", 0.5, 0.8, 1, 1);
    fGFW->AddRegion("refP", 0.4, 0.8, 1, 1);
    fGFW->AddRegion("refN", -0.8, -0.4, 1, 1);
    fGFW->AddRegion("poiN", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 2);
    fGFW->AddRegion("olN", -0.8, -0.4, 1, 4);

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {3 -3}", "ChFull32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {4 -4}", "ChFull42", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 -2 -2}", "ChFull24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 2 -2 -2 -2}", "ChFull26", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN04 {2} refP04 {-2}", "Ch04Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN06 {2} refP06 {-2}", "Ch06Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "Ch08Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {2} refP10 {-2}", "Ch10Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN04 {3} refP04 {-3}", "Ch04Gap32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN06 {3} refP06 {-3}", "Ch06Gap32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {3} refP08 {-3}", "Ch08Gap32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {3} refP10 {-3}", "Ch10Gap32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN04 {4} refP04 {-4}", "Ch04Gap42", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN06 {4} refP06 {-4}", "Ch06Gap42", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {4} refP08 {-4}", "Ch08Gap42", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {4} refP10 {-4}", "Ch10Gap42", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN {2} refP {-2}", "ChGap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiN refN | olN {2} refP {-2}", "ChGap22", kTRUE));
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
    fGFW->CreateRegions();
  }

  template <char... chars>
  void FillProfile(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& cent)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1)
        registry.fill(tarName, cent, val, dnx);
      return;
    }
    return;
  }

  template <char... chars, char... chars2>
  void FillpTvnProfile(const GFW::CorrConfig& corrconf, const double& sum_pt, const double& WeffEvent, const ConstStr<chars...>& vnWeff, const ConstStr<chars2...>& vnpT, const double& cent)
  {
    double meanPt = sum_pt / WeffEvent;
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1) {
        registry.fill(vnWeff, cent, val, dnx * WeffEvent);
        registry.fill(vnpT, cent, val * meanPt, dnx * WeffEvent);
      }
      return;
    }
    return;
  }

  void FillpTvnProfile(const GFW::CorrConfig& corrconf, const double& sum_pt, const double& WeffEvent, std::shared_ptr<TProfile> vnWeff, std::shared_ptr<TProfile> vnpT, const double& cent)
  {
    double meanPt = sum_pt / WeffEvent;
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1) {
        vnWeff->Fill(cent, val, dnx * WeffEvent);
        vnpT->Fill(cent, val * meanPt, dnx * WeffEvent);
      }
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
    correctionsLoaded = true;
  }

  bool setCurrentParticleWeights(float& weight_nue, float& weight_nua, const float& phi, const float& eta, const float& pt, const float& vtxz)
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
      weight_nua = mAcceptance->GetNUA(phi, eta, vtxz);
    else
      weight_nua = 1;
    return true;
  }

  void process(aodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, aodTracks const& tracks)
  {
    if (!collision.sel8())
      return;
    int Ntot = tracks.size();
    if (Ntot < 1)
      return;
    float l_Random = fRndm->Rndm();
    float vtxz = collision.posZ();
    registry.fill(HIST("hVtxZ"), vtxz);
    registry.fill(HIST("hMult"), Ntot);
    registry.fill(HIST("hCent"), collision.centFT0C());
    fGFW->Clear();
    const auto cent = collision.centFT0C();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    loadCorrections(bc.timestamp());

    float weff = 1, wacc = 1;
    double weffEvent = 0, waccEvent = 0;
    int TrackNum = 0;
    double ptSum = 0., ptSum_Gap08 = 0.;
    double weffEvent_WithinGap08 = 0., weffEventSquare_WithinGap08 = 0.;
    double sum_ptSquare_wSquare_WithinGap08 = 0., sum_pt_wSquare_WithinGap08 = 0.;

    for (auto& track : tracks) {
      double pt = track.pt();
      double eta = track.eta();
      double phi = track.phi();
      if (cfgOutputNUAWeights)
        fWeights->Fill(phi, eta, vtxz, pt, cent, 0);
      if (!setCurrentParticleWeights(weff, wacc, phi, eta, pt, vtxz))
        continue;
      bool WithinPtPOI = (cfgCutPtPOIMin < pt) && (pt < cfgCutPtPOIMax); // within POI pT range
      bool WithinPtRef = (cfgCutPtMin < pt) && (pt < cfgCutPtMax);       // within RF pT range
      bool WithinEtaGap08 = (eta >= -0.4) && (eta <= 0.4);
      if (WithinPtRef) {
        registry.fill(HIST("hPhi"), phi);
        registry.fill(HIST("hEta"), track.eta());
        registry.fill(HIST("hPt"), pt);
        weffEvent += weff;
        waccEvent += wacc;
        ptSum += weff * pt;
        TrackNum++;
        if (WithinEtaGap08) {
          ptSum_Gap08 += weff * pt;
          sum_pt_wSquare_WithinGap08 += weff * weff * pt;
          sum_ptSquare_wSquare_WithinGap08 += weff * weff * pt * pt;
          weffEvent_WithinGap08 += weff;
          weffEventSquare_WithinGap08 += weff * weff;
        }
      }
      if (WithinPtRef)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, phi, wacc * weff, 1);
      if (WithinPtPOI)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, phi, wacc * weff, 2);
      if (WithinPtPOI && WithinPtRef)
        fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, phi, wacc * weff, 4);
    }

    double WeffEvent_diff_WithGap08 = weffEvent_WithinGap08 * weffEvent_WithinGap08 - weffEventSquare_WithinGap08;
    // Filling TProfile
    // MeanPt
    registry.fill(HIST("hMeanPt"), cent, ptSum / weffEvent, weffEvent);
    if (weffEvent_WithinGap08 > 1e-6)
      registry.fill(HIST("hMeanPtWithinGap08"), cent, ptSum_Gap08 / weffEvent_WithinGap08, weffEvent_WithinGap08);
    // v22-Pt
    // c22_gap8 * pt_withGap8
    if (weffEvent_WithinGap08 > 1e-6)
      FillpTvnProfile(corrconfigs.at(7), ptSum_Gap08, weffEvent_WithinGap08, HIST("c22_gap08_Weff"), HIST("c22_gap08_trackMeanPt"), cent);
    // PtVariance
    if (WeffEvent_diff_WithGap08 > 1e-6) {
      registry.fill(HIST("PtVariance_partA_WithinGap08"), cent,
                    (ptSum_Gap08 * ptSum_Gap08 - sum_ptSquare_wSquare_WithinGap08) / WeffEvent_diff_WithGap08,
                    WeffEvent_diff_WithGap08);
      registry.fill(HIST("PtVariance_partB_WithinGap08"), cent,
                    (weffEvent_WithinGap08 * ptSum_Gap08 - sum_pt_wSquare_WithinGap08) / WeffEvent_diff_WithGap08,
                    WeffEvent_diff_WithGap08);
    }

    // Filling Bootstrap Samples
    int SampleIndex = static_cast<int>(cfgNbootstrap * l_Random);
    if (weffEvent_WithinGap08 > 1e-6)
      BootstrapArray[SampleIndex][kMeanPt_InGap08]->Fill(cent, ptSum_Gap08 / weffEvent_WithinGap08, weffEvent_WithinGap08);
    if (weffEvent_WithinGap08 > 1e-6)
      FillpTvnProfile(corrconfigs.at(7), ptSum_Gap08, weffEvent_WithinGap08, BootstrapArray[SampleIndex][kC22_Gap08_Weff], BootstrapArray[SampleIndex][kC22_Gap08_MeanPt], cent);
    if (WeffEvent_diff_WithGap08 > 1e-6) {
      BootstrapArray[SampleIndex][kPtVarParA_InGap08]->Fill(cent,
                                                            (ptSum_Gap08 * ptSum_Gap08 - sum_ptSquare_wSquare_WithinGap08) / WeffEvent_diff_WithGap08,
                                                            WeffEvent_diff_WithGap08);
      BootstrapArray[SampleIndex][kPtVarParB_InGap08]->Fill(cent,
                                                            (weffEvent_WithinGap08 * ptSum_Gap08 - sum_pt_wSquare_WithinGap08) / WeffEvent_diff_WithGap08,
                                                            WeffEvent_diff_WithGap08);
    }

    // Filling Flow Container
    for (uint l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
      FillFC(corrconfigs.at(l_ind), cent, l_Random);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowPbPbTask>(cfgc)};
}
