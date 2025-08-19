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

/// \file   flowPtEfficiency.cxx
/// \author Mingrui Zhao (mingrui.zhao@cern.ch), Zhiyong Lu (zhiyong.lu@cern.ch), Tao Jiang (tao.jiang@cern.ch)
/// \since  Jun/08/2023
/// \brief  a task to calculate the pt efficiency

#include "FlowContainer.h"
#include "GFW.h"
#include "GFWCumulant.h"
#include "GFWPowerArray.h"
#include "GFWWeights.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include <CCDB/BasicCCDBManager.h>

#include <TF1.h>
#include <TPDGCode.h>
#include <TProfile.h>
#include <TRandom3.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowPtEfficiency {

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 1000.0f, "Maximal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgkIsTrackGlobal, bool, false, "GlobalTrack requirement for tracks")
  O2_DEFINE_CONFIGURABLE(cfgTrkSelRun3ITSMatch, bool, false, "GlobalTrackRun3ITSMatching::Run3ITSall7Layers selection")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5f, "max chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCclu, float, 70.0f, "minimum TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutITSclu, float, 5.0f, "minimum ITS clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCcrossedrows, float, 70.0f, "minimum TPC crossed rows")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxy, float, 0.2f, "DCAxy cut for tracks")
  O2_DEFINE_CONFIGURABLE(cfgDCAxyNSigma, float, 7, "Cut on number of sigma deviations from expected DCA in the transverse direction");
  O2_DEFINE_CONFIGURABLE(cfgDCAxyFunction, std::string, "(0.0015+0.005/(x^1.1))", "Functional form of pt-dependent DCAxy cut");
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "DCAz cut for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxyppPass3Enabled, bool, false, "switch of ppPass3 DCAxy pt dependent cut")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAzPtDepEnabled, bool, false, "switch of DCAz pt dependent cut")
  O2_DEFINE_CONFIGURABLE(cfgEnableITSCuts, bool, true, "switch of enabling ITS based track selection cuts")
  O2_DEFINE_CONFIGURABLE(cfgSelRunNumberEnabled, bool, false, "switch of run number selection")
  O2_DEFINE_CONFIGURABLE(cfgFlowEnabled, bool, false, "switch of calculating flow")
  O2_DEFINE_CONFIGURABLE(cfgFlowNbootstrap, int, 30, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgFlowCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgFlowCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgFlowCutPtRefMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgFlowCutPtRefMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCentVsIPTruth, std::string, "", "CCDB path to centrality vs IP truth")
  O2_DEFINE_CONFIGURABLE(cfgCentVsIPReco, std::string, "", "CCDB path to centrality vs IP reco")
  O2_DEFINE_CONFIGURABLE(cfgFlowAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgFlowEfficiency, std::string, "", "CCDB path to efficiency object")
  Configurable<std::vector<int>> cfgRunNumberList{"cfgRunNumberList", std::vector<int>{-1}, "runnumber list in consideration for analysis"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10}, "pt axis for histograms"};
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "X axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {100, 0.0f, constants::math::TwoPI}, ""};
  ConfigurableAxis axisB{"axisB", {100, 0.0f, 20.0f}, "b (fm)"};
  ConfigurableAxis axisNch{"axisNch", {6000, 0, 6000}, "N_{ch}"};

  // Filter the tracks
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);
  using MyTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::McTrackLabels>>;

  // Filter for collisions
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  using MyCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>>;

  // Filter for MCParticle
  Filter particleFilter = (nabs(aod::mcparticle::eta) < cfgCutEta) && (aod::mcparticle::pt > cfgCutPtMin) && (aod::mcparticle::pt < cfgCutPtMax);
  using MyMcParticles = soa::Filtered<aod::McParticles>;

  // Filter for MCcollisions
  Filter mccollisionFilter = nabs(aod::mccollision::posZ) < cfgCutVertex;
  using MyMcCollisions = soa::Filtered<aod::McCollisions>;

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  // Additional filters for tracks
  TrackSelection myTrackSel;

  // Cent vs IP
  TH1D* mCentVsIPTruth = nullptr;
  bool centVsIPTruthLoaded = false;
  TH1D* mCentVsIPReco = nullptr;
  bool centVsIPRecoLoaded = false;

  // corrections
  TH1D* mEfficiency = nullptr;
  GFWWeights* mAcceptance = nullptr;
  bool correctionsLoaded = false;

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  // Define the output
  HistogramRegistry registry{"registry"};
  OutputObj<FlowContainer> fFCTrue{FlowContainer("FlowContainerTrue")};
  OutputObj<FlowContainer> fFCReco{FlowContainer("FlowContainerReco")};
  OutputObj<GFWWeights> fWeights{GFWWeights("weights")};
  GFW* fGFWTrue = new GFW();
  GFW* fGFWReco = new GFW();
  TAxis* fPtAxis;
  std::vector<GFW::CorrConfig> corrconfigsTruth;
  std::vector<GFW::CorrConfig> corrconfigsReco;
  TRandom3* fRndm = new TRandom3(0);
  TF1* fPtDepDCAxy = nullptr;

  bool isStable(int pdg)
  {
    if (std::abs(pdg) == PDG_t::kPiPlus)
      return true;
    if (std::abs(pdg) == PDG_t::kKPlus)
      return true;
    if (std::abs(pdg) == PDG_t::kProton)
      return true;
    if (std::abs(pdg) == PDG_t::kElectron)
      return true;
    if (std::abs(pdg) == PDG_t::kMuonMinus)
      return true;
    return false;
  }

  void init(InitContext const&)
  {
    const AxisSpec axisVertex{20, -10, 10, "Vtxz (cm)"};
    const AxisSpec axisEta{20, -1., 1., "#eta"};
    const AxisSpec axisCounter{1, 0, +1, ""};
    // create histograms
    registry.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    registry.add("hPtMCRec", "Monte Carlo Reco", {HistType::kTH1D, {axisPt}});
    registry.add("hPtNchMCRec", "Reco production; pT (GeV/c); multiplicity", {HistType::kTH2D, {axisPt, axisNch}});
    registry.add("hBVsPtVsPhiRec", "hBVsPtVsPhiRec", HistType::kTH3D, {axisB, axisPhi, axisPt});
    registry.add("hEtaPtVzRec", "hEtaPtVz Reconstructed", HistType::kTH3D, {axisEta, axisPt, axisVertex});

    registry.add("mcEventCounter", "Monte Carlo Truth EventCounter", kTH1F, {axisCounter});
    registry.add("hPtMCGen", "Monte Carlo Truth", {HistType::kTH1D, {axisPt}});
    registry.add("hPtNchMCGen", "Truth production; pT (GeV/c); multiplicity", {HistType::kTH2D, {axisPt, axisNch}});
    registry.add("numberOfRecoCollisions", "numberOfRecoCollisions", kTH1F, {{10, -0.5f, 9.5f}});
    registry.add("hBVsPtVsPhiTrue", "hBVsPtVsPhiTrue", HistType::kTH3D, {axisB, axisPhi, axisPt});
    registry.add("hEtaPtVzTrue", "hEtaPtVz True", HistType::kTH3D, {axisEta, axisPt, axisVertex});

    if (cfgFlowEnabled) {
      registry.add("hImpactParameterReco", "hImpactParameterReco", {HistType::kTH1D, {axisB}});
      registry.add("hImpactParameterTruth", "hImpactParameterTruth", {HistType::kTH1D, {axisB}});
      registry.add("hPhi", "#phi distribution", {HistType::kTH1D, {axisPhi}});
      registry.add("hPhiMCTruth", "#phi distribution", {HistType::kTH1D, {axisPhi}});
      registry.add("hPhiWeighted", "corrected #phi distribution", {HistType::kTH1D, {axisPhi}});

      o2::framework::AxisSpec axis = axisPt;
      int nPtBins = axis.binEdges.size() - 1;
      double* ptBins = &(axis.binEdges)[0];
      fPtAxis = new TAxis(nPtBins, ptBins);

      fWeights->setPtBins(nPtBins, ptBins);
      fWeights->init(true, false);

      TObjArray* oba = new TObjArray();
      oba->Add(new TNamed("ChFull22", "ChFull22"));
      for (auto i = 0; i < fPtAxis->GetNbins(); i++)
        oba->Add(new TNamed(Form("ChFull22_pt_%i", i + 1), "ChFull22_pTDiff"));
      oba->Add(new TNamed("Ch10Gap22", "Ch10Gap22"));
      for (auto i = 0; i < fPtAxis->GetNbins(); i++)
        oba->Add(new TNamed(Form("Ch10Gap22_pt_%i", i + 1), "Ch10Gap22_pTDiff"));
      fFCTrue->SetName("FlowContainerTrue");
      fFCTrue->SetXAxis(fPtAxis);
      fFCTrue->Initialize(oba, axisCentrality, cfgFlowNbootstrap);
      fFCReco->SetName("FlowContainerReco");
      fFCReco->SetXAxis(fPtAxis);
      fFCReco->Initialize(oba, axisCentrality, cfgFlowNbootstrap);
      delete oba;

      fGFWTrue->AddRegion("full", -0.8, 0.8, 1, 1);
      fGFWTrue->AddRegion("refN10", -0.8, -0.5, 1, 1);
      fGFWTrue->AddRegion("refP10", 0.5, 0.8, 1, 1);
      fGFWTrue->AddRegion("poiN10", -0.8, -0.5, 1 + fPtAxis->GetNbins(), 2);
      fGFWTrue->AddRegion("poifull", -0.8, 0.8, 1 + fPtAxis->GetNbins(), 2);
      fGFWTrue->AddRegion("olN10", -0.8, -0.5, 1 + fPtAxis->GetNbins(), 4);
      fGFWTrue->AddRegion("olfull", -0.8, 0.8, 1 + fPtAxis->GetNbins(), 4);
      corrconfigsTruth.push_back(fGFWTrue->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE));
      corrconfigsTruth.push_back(fGFWTrue->GetCorrelatorConfig("poifull full | olfull {2 -2}", "ChFull22", kTRUE));
      corrconfigsTruth.push_back(fGFWTrue->GetCorrelatorConfig("refN10 {2} refP10 {-2}", "Ch10Gap22", kFALSE));
      corrconfigsTruth.push_back(fGFWTrue->GetCorrelatorConfig("poiN10 refN10 | olN10 {2} refP10 {-2}", "Ch10Gap22", kTRUE));
      fGFWTrue->CreateRegions();

      fGFWReco->AddRegion("full", -0.8, 0.8, 1, 1);
      fGFWReco->AddRegion("refN10", -0.8, -0.5, 1, 1);
      fGFWReco->AddRegion("refP10", 0.5, 0.8, 1, 1);
      fGFWReco->AddRegion("poiN10", -0.8, -0.5, 1 + fPtAxis->GetNbins(), 2);
      fGFWReco->AddRegion("poifull", -0.8, 0.8, 1 + fPtAxis->GetNbins(), 2);
      fGFWReco->AddRegion("olN10", -0.8, -0.5, 1 + fPtAxis->GetNbins(), 4);
      fGFWReco->AddRegion("olfull", -0.8, 0.8, 1 + fPtAxis->GetNbins(), 4);
      corrconfigsReco.push_back(fGFWReco->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE));
      corrconfigsReco.push_back(fGFWReco->GetCorrelatorConfig("poifull full | olfull {2 -2}", "ChFull22", kTRUE));
      corrconfigsReco.push_back(fGFWReco->GetCorrelatorConfig("refN10 {2} refP10 {-2}", "Ch10Gap22", kFALSE));
      corrconfigsReco.push_back(fGFWReco->GetCorrelatorConfig("poiN10 refN10 | olN10 {2} refP10 {-2}", "Ch10Gap22", kTRUE));
      fGFWReco->CreateRegions();
    }

    if (cfgEnableITSCuts) {
      if (cfgTrkSelRun3ITSMatch) {
        myTrackSel = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSall7Layers, TrackSelection::GlobalTrackRun3DCAxyCut::Default);
      } else {
        myTrackSel = getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, TrackSelection::GlobalTrackRun3DCAxyCut::Default);
      }
    }
    if (cfgCutDCAxyppPass3Enabled) {
      myTrackSel.SetMaxDcaXYPtDep([](float pt) { return 0.004f + 0.013f / pt; });
    } else {
      if (cfgCutDCAxy != 0.0) {
        myTrackSel.SetMaxDcaXY(cfgCutDCAxy);
      } else {
        fPtDepDCAxy = new TF1("ptDepDCAxy", Form("[0]*%s", cfgDCAxyFunction->c_str()), 0.001, 100);
        fPtDepDCAxy->SetParameter(0, cfgDCAxyNSigma);
        LOGF(info, "DCAxy pt-dependence function: %s", Form("[0]*%s", cfgDCAxyFunction->c_str()));
        myTrackSel.SetMaxDcaXYPtDep([fPtDepDCAxy = this->fPtDepDCAxy](float pt) { return fPtDepDCAxy->Eval(pt); });
      }
    }
    myTrackSel.SetMinNClustersTPC(cfgCutTPCclu);
    myTrackSel.SetMinNCrossedRowsTPC(cfgCutTPCcrossedrows);
    if (cfgEnableITSCuts)
      myTrackSel.SetMinNClustersITS(cfgCutITSclu);
    if (!cfgCutDCAzPtDepEnabled)
      myTrackSel.SetMaxDcaZ(cfgCutDCAz);
  }

  template <char... chars>
  void fillProfile(GFW* fGFW, const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& cent)
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

  void fillFC(GFW* fGFW, bool isMCTruth, const GFW::CorrConfig& corrconf, const double& cent, const double& rndm)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (std::fabs(val) < 1) {
        if (isMCTruth)
          fFCTrue->FillProfile(corrconf.Head.c_str(), cent, val, dnx, rndm);
        else
          fFCReco->FillProfile(corrconf.Head.c_str(), cent, val, dnx, rndm);
      }
      return;
    }
    for (auto i = 1; i <= fPtAxis->GetNbins(); i++) {
      dnx = fGFW->Calculate(corrconf, i - 1, kTRUE).real();
      if (dnx == 0)
        continue;
      val = fGFW->Calculate(corrconf, i - 1, kFALSE).real() / dnx;
      if (std::fabs(val) < 1) {
        if (isMCTruth)
          fFCTrue->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i), cent, val, dnx, rndm);
        else
          fFCReco->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i), cent, val, dnx, rndm);
      }
    }
    return;
  }

  void loadCentVsIPTruth(uint64_t timestamp)
  {
    if (centVsIPTruthLoaded)
      return;
    if (cfgCentVsIPTruth.value.empty() == false) {
      mCentVsIPTruth = ccdb->getForTimeStamp<TH1D>(cfgCentVsIPTruth, timestamp);
      if (mCentVsIPTruth)
        LOGF(info, "Loaded CentVsIPTruth weights from %s (%p)", cfgCentVsIPTruth.value.c_str(), (void*)mCentVsIPTruth);
      else
        LOGF(fatal, "Failed to load CentVsIPTruth weights from %s", cfgCentVsIPTruth.value.c_str());

      centVsIPTruthLoaded = true;
    } else {
      LOGF(fatal, "when calculate flow, Cent Vs IP distribution must be provided");
    }
  }

  void loadCentVsIPReco(uint64_t timestamp)
  {
    if (centVsIPRecoLoaded)
      return;
    if (cfgCentVsIPReco.value.empty() == false) {
      mCentVsIPReco = ccdb->getForTimeStamp<TH1D>(cfgCentVsIPReco, timestamp);
      if (mCentVsIPReco)
        LOGF(info, "Loaded CentVsIPReco weights from %s (%p)", cfgCentVsIPReco.value.c_str(), (void*)mCentVsIPReco);
      else
        LOGF(fatal, "Failed to load CentVsIPReco weights from %s", cfgCentVsIPReco.value.c_str());

      centVsIPRecoLoaded = true;
    } else {
      LOGF(fatal, "when calculate flow, Cent Vs IP distribution must be provided");
    }
  }

  void loadCorrections(uint64_t timestamp)
  {
    if (correctionsLoaded)
      return;
    if (cfgFlowAcceptance.value.empty() == false) {
      mAcceptance = ccdb->getForTimeStamp<GFWWeights>(cfgFlowAcceptance, timestamp);
      if (mAcceptance)
        LOGF(info, "Loaded acceptance weights from %s (%p)", cfgFlowAcceptance.value.c_str(), (void*)mAcceptance);
      else
        LOGF(warning, "Could not load acceptance weights from %s (%p)", cfgFlowAcceptance.value.c_str(), (void*)mAcceptance);
    }
    if (cfgFlowEfficiency.value.empty() == false) {
      mEfficiency = ccdb->getForTimeStamp<TH1D>(cfgFlowEfficiency, timestamp);
      if (mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgFlowEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgFlowEfficiency.value.c_str(), (void*)mEfficiency);
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
    if (cfgkIsTrackGlobal && !track.isGlobalTrack()) {
      return false;
    }
    if (cfgCutDCAzPtDepEnabled && (track.dcaZ() > (0.004f + 0.013f / track.pt()))) {
      return false;
    }
    return myTrackSel.IsSelected(track);
  }

  void processReco(MyCollisions::iterator const& collision, aod::BCsWithTimestamps const&, MyTracks const& tracks, aod::McParticles const&, aod::McCollisions const&)
  {
    registry.fill(HIST("eventCounter"), 0.5);
    if (!collision.sel8())
      return;
    if (tracks.size() < 1)
      return;
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int runNumber = bc.runNumber();
    if (cfgSelRunNumberEnabled) {
      if (!std::count(cfgRunNumberList.value.begin(), cfgRunNumberList.value.end(), runNumber))
        return;
    }
    float imp = 0;
    bool impFetched = false;
    float evPhi = 0;
    float centrality = 0.;
    float lRandom = fRndm->Rndm();
    float vtxz = collision.posZ();
    float wacc = 1.0f;
    float weff = 1.0f;
    if (cfgFlowEnabled) {
      loadCentVsIPReco(bc.timestamp());
      loadCorrections(bc.timestamp());

      fGFWReco->Clear();
    }
    for (const auto& track : tracks) {
      if (!trackSelected(track))
        continue;
      if (track.has_mcParticle()) {
        auto mcParticle = track.mcParticle();
        if (cfgFlowEnabled && !impFetched) {
          auto mcCollision = mcParticle.mcCollision();
          imp = mcCollision.impactParameter();
          registry.fill(HIST("hImpactParameterReco"), imp);
          centrality = mCentVsIPReco->GetBinContent(mCentVsIPReco->GetXaxis()->FindBin(imp));
          evPhi = RecoDecay::constrainAngle(mcCollision.eventPlaneAngle());
          impFetched = true;
        }
        if (isStable(mcParticle.pdgCode())) {
          registry.fill(HIST("hPtMCRec"), track.pt());
          registry.fill(HIST("hPtNchMCRec"), track.pt(), tracks.size());
          registry.fill(HIST("hEtaPtVzRec"), track.eta(), track.pt(), vtxz);

          if (cfgFlowEnabled) {
            float deltaPhi = RecoDecay::constrainAngle(track.phi() - evPhi);
            registry.fill(HIST("hBVsPtVsPhiRec"), imp, deltaPhi, track.pt());
            bool withinPtPOI = (cfgFlowCutPtPOIMin < track.pt()) && (track.pt() < cfgFlowCutPtPOIMax); // within POI pT range
            bool withinPtRef = (cfgFlowCutPtRefMin < track.pt()) && (track.pt() < cfgFlowCutPtRefMax); // within RF pT range
            if (withinPtRef)
              fWeights->fill(track.phi(), track.eta(), vtxz, track.pt(), centrality, 0);
            if (!setCurrentParticleWeights(weff, wacc, track.phi(), track.eta(), track.pt(), vtxz))
              continue;
            if (withinPtRef) {
              registry.fill(HIST("hPhi"), track.phi());
              registry.fill(HIST("hPhiWeighted"), track.phi(), wacc);
            }
            if (withinPtRef)
              fGFWReco->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), wacc * weff, 1);
            if (withinPtPOI)
              fGFWReco->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), wacc * weff, 2);
            if (withinPtPOI && withinPtRef)
              fGFWReco->Fill(track.eta(), fPtAxis->FindBin(track.pt()) - 1, track.phi(), wacc * weff, 4);
          }
        }
      }
    }
    if (cfgFlowEnabled) {
      // Filling Flow Container
      for (uint l_ind = 0; l_ind < corrconfigsReco.size(); l_ind++) {
        fillFC(fGFWReco, false, corrconfigsReco.at(l_ind), centrality, lRandom);
      }
    }
  }
  PROCESS_SWITCH(FlowPtEfficiency, processReco, "process reconstructed information", true);

  void processSim(MyMcCollisions::iterator const& mcCollision, aod::BCsWithTimestamps const&, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions, MyMcParticles const& mcParticles, MyTracks const& tracks)
  {
    if (cfgSelRunNumberEnabled) {
      auto bc = mcCollision.bc_as<aod::BCsWithTimestamps>();
      int runNumber = bc.runNumber();
      if (!std::count(cfgRunNumberList.value.begin(), cfgRunNumberList.value.end(), runNumber))
        return;
    }

    float imp = mcCollision.impactParameter();
    float evPhi = RecoDecay::constrainAngle(mcCollision.eventPlaneAngle());
    float centrality = 0.;
    if (cfgFlowEnabled) {
      registry.fill(HIST("hImpactParameterTruth"), imp);
      auto bc = mcCollision.bc_as<aod::BCsWithTimestamps>();
      loadCentVsIPTruth(bc.timestamp());
      centrality = mCentVsIPTruth->GetBinContent(mCentVsIPTruth->GetXaxis()->FindBin(imp));

      fGFWTrue->Clear();
    }
    float lRandom = fRndm->Rndm();
    float wacc = 1.0f;
    float weff = 1.0f;
    float vtxz = mcCollision.posZ();

    if (collisions.size() > -1) {
      registry.fill(HIST("mcEventCounter"), 0.5);

      registry.fill(HIST("numberOfRecoCollisions"), collisions.size()); // number of times coll was reco-ed

      std::vector<int> numberOfTracks;
      for (auto const& collision : collisions) {
        auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());
        numberOfTracks.emplace_back(groupedTracks.size());
      }

      for (const auto& mcParticle : mcParticles) {
        if (mcParticle.isPhysicalPrimary() && isStable(mcParticle.pdgCode())) {
          registry.fill(HIST("hPtMCGen"), mcParticle.pt());
          if (collisions.size() > 0) {
            registry.fill(HIST("hPtNchMCGen"), mcParticle.pt(), numberOfTracks[0]);
          }
          registry.fill(HIST("hEtaPtVzTrue"), mcParticle.eta(), mcParticle.pt(), vtxz);

          if (cfgFlowEnabled) {
            float deltaPhi = RecoDecay::constrainAngle(mcParticle.phi() - evPhi);
            registry.fill(HIST("hBVsPtVsPhiTrue"), imp, deltaPhi, mcParticle.pt());
            bool withinPtPOI = (cfgFlowCutPtPOIMin < mcParticle.pt()) && (mcParticle.pt() < cfgFlowCutPtPOIMax); // within POI pT range
            bool withinPtRef = (cfgFlowCutPtRefMin < mcParticle.pt()) && (mcParticle.pt() < cfgFlowCutPtRefMax); // within RF pT range
            if (withinPtRef) {
              registry.fill(HIST("hPhiMCTruth"), mcParticle.phi());
            }
            if (withinPtRef)
              fGFWTrue->Fill(mcParticle.eta(), fPtAxis->FindBin(mcParticle.pt()) - 1, mcParticle.phi(), wacc * weff, 1);
            if (withinPtPOI)
              fGFWTrue->Fill(mcParticle.eta(), fPtAxis->FindBin(mcParticle.pt()) - 1, mcParticle.phi(), wacc * weff, 2);
            if (withinPtPOI && withinPtRef)
              fGFWTrue->Fill(mcParticle.eta(), fPtAxis->FindBin(mcParticle.pt()) - 1, mcParticle.phi(), wacc * weff, 4);
          }
        }
      }
      if (cfgFlowEnabled) {
        // Filling Flow Container
        for (uint l_ind = 0; l_ind < corrconfigsTruth.size(); l_ind++) {
          fillFC(fGFWTrue, true, corrconfigsTruth.at(l_ind), centrality, lRandom);
        }
      }
    }
  }
  PROCESS_SWITCH(FlowPtEfficiency, processSim, "process pure simulation information", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowPtEfficiency>(cfgc)};
}
