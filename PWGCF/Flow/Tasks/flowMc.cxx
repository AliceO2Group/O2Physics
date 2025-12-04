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

/// \file   flowMc.cxx
/// \author Zhiyong Lu (zhiyong.lu@cern.ch)
/// \since  Feb/5/2025
/// \brief  QC of synthetic flow exercise

#include <CCDB/BasicCCDBManager.h>
#include <vector>
#include <string>
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGMM/Mult/DataModel/Index.h" // for Particles2Tracks table
#include "GFWPowerArray.h"
#include "GFW.h"
#include "GFWCumulant.h"
#include "GFWWeights.h"
#include "FlowContainer.h"
#include <TProfile.h>
#include <TRandom3.h>
#include <TPDGCode.h>
#include <TF1.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowMc {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> minB{"minB", 0.0f, "min impact parameter"};
  Configurable<float> maxB{"maxB", 20.0f, "max impact parameter"};
  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 1000.0f, "Maximal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeights, bool, false, "Fill and output NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgCutPtRefMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtRefMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCclu, float, 50.0f, "minimum TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutITSclu, float, 5.0f, "minimum ITS clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5f, "max chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutTPCcrossedrows, float, 70.0f, "minimum TPC crossed rows")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxy, float, 0.2f, "DCAxy cut for tracks")
  O2_DEFINE_CONFIGURABLE(cfgDCAxyNSigma, float, 7, "Cut on number of sigma deviations from expected DCA in the transverse direction");
  O2_DEFINE_CONFIGURABLE(cfgDCAxyFunction, std::string, "(0.0015+0.005/(x^1.1))", "Functional form of pt-dependent DCAxy cut");
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "DCAz cut for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAzPtDepEnabled, bool, false, "switch of DCAz pt dependent cut")
  O2_DEFINE_CONFIGURABLE(cfgEnableITSCuts, bool, true, "switch of enabling ITS based track selection cuts")
  O2_DEFINE_CONFIGURABLE(cfgTrkSelRun3ITSMatch, bool, false, "GlobalTrackRun3ITSMatching::Run3ITSall7Layers selection")
  O2_DEFINE_CONFIGURABLE(cfgFlowAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgFlowEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgCentVsIPTruth, std::string, "", "CCDB path to centrality vs IP truth")
  O2_DEFINE_CONFIGURABLE(cfgIsGlobalTrack, bool, false, "Use global tracks instead of hasTPC&&hasITS")
  O2_DEFINE_CONFIGURABLE(cfgK0Lambda0Enabled, bool, false, "Add K0 and Lambda0")
  O2_DEFINE_CONFIGURABLE(cfgFlowCumulantEnabled, bool, false, "switch of calculating flow")
  O2_DEFINE_CONFIGURABLE(cfgFlowCumulantNbootstrap, int, 30, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgTrackDensityCorrUse, bool, false, "Use track density efficiency correction")
  O2_DEFINE_CONFIGURABLE(cfgTrackDensityCorrSlopeFactor, float, 1.0f, "A factor to scale the track density efficiency slope")
  O2_DEFINE_CONFIGURABLE(cfgRecoEvRejectMC, bool, false, "reject both MC and Reco events when reco do not pass")
  O2_DEFINE_CONFIGURABLE(cfgRecoEvSel8, bool, false, "require sel8 for reconstruction events")
  O2_DEFINE_CONFIGURABLE(cfgRecoEvkIsGoodITSLayersAll, bool, false, "require kIsGoodITSLayersAll for reconstruction events")
  O2_DEFINE_CONFIGURABLE(cfgRecoEvkNoSameBunchPileup, bool, false, "require kNoSameBunchPileup for reconstruction events")
  O2_DEFINE_CONFIGURABLE(cfgRecoEvSelkIsGoodZvtxFT0vsPV, bool, false, "removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference, use this cut at low multiplicities with caution")
  O2_DEFINE_CONFIGURABLE(cfgRecoEvSelkNoITSROFrameBorder, bool, false, "reject events at ITS ROF border")
  O2_DEFINE_CONFIGURABLE(cfgRecoEvSelkNoTimeFrameBorder, bool, false, "reject events at TF border")
  O2_DEFINE_CONFIGURABLE(cfgRecoEvSelkNoCollInTimeRangeStandard, bool, false, "no collisions in specified time range")
  O2_DEFINE_CONFIGURABLE(cfgRecoEvSelkNoCollInRofStandard, bool, false, "no other collisions in this Readout Frame with per-collision multiplicity above threshold")
  O2_DEFINE_CONFIGURABLE(cfgRecoEvSelkNoHighMultCollInPrevRof, bool, false, "veto an event if FT0C amplitude in previous ITS ROF is above threshold")

  Configurable<std::vector<double>> cfgTrackDensityP0{"cfgTrackDensityP0", std::vector<double>{0.6003720411, 0.6152630970, 0.6288860646, 0.6360694031, 0.6409494798, 0.6450540203, 0.6482117301, 0.6512592056, 0.6640008690, 0.6862631416, 0.7005738691, 0.7106567432, 0.7170728333}, "parameter 0 for track density efficiency correction"};
  Configurable<std::vector<double>> cfgTrackDensityP1{"cfgTrackDensityP1", std::vector<double>{-1.007592e-05, -8.932635e-06, -9.114538e-06, -1.054818e-05, -1.220212e-05, -1.312304e-05, -1.376433e-05, -1.412813e-05, -1.289562e-05, -1.050065e-05, -8.635725e-06, -7.380821e-06, -6.201250e-06}, "parameter 1 for track density efficiency correction"};
  float maxEta = 0.8;

  ConfigurableAxis axisB{"axisB", {100, 0.0f, 20.0f}, ""};
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "X axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {100, 0.0f, constants::math::TwoPI}, ""};
  ConfigurableAxis axisNch{"axisNch", {300, 0.0f, 3000.0f}, "Nch in |eta|<0.8"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f}, "pt axis"};

  // Filter for MCcollisions
  Filter mccollisionFilter = nabs(aod::mccollision::posZ) < cfgCutVertex;
  using FilteredMcCollisions = soa::Filtered<aod::McCollisions>;
  // Filter for MCParticle
  Filter particleFilter = (nabs(aod::mcparticle::eta) < cfgCutEta) && (aod::mcparticle::pt > cfgCutPtMin) && (aod::mcparticle::pt < cfgCutPtMax);
  using FilteredMcParticles = soa::Filtered<soa::Join<aod::McParticles, aod::ParticlesToTracks>>;
  // Filter for reco tracks
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);
  using FilteredTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::McTrackLabels>>;

  // using FilteredTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TrackSelection>;

  // Additional filters for tracks
  TrackSelection myTrackSel;
  TF1* fPtDepDCAxy = nullptr;

  // Cent vs IP
  TH1D* mCentVsIPTruth = nullptr;
  bool centVsIPTruthLoaded = false;

  // Corrections
  TH1D* mEfficiency = nullptr;
  GFWWeights* mAcceptance = nullptr;
  bool correctionsLoaded = false;

  std::vector<TF1*> funcEff;
  TH1D* hFindPtBin;
  TF1* funcV2;
  TF1* funcV3;
  TF1* funcV4;

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  OutputObj<GFWWeights> fWeights{GFWWeights("weights")};
  OutputObj<FlowContainer> fFCTrue{FlowContainer("FlowContainerTrue")};
  OutputObj<FlowContainer> fFCReco{FlowContainer("FlowContainerReco")};
  GFW* fGFWTrue = new GFW();
  GFW* fGFWReco = new GFW();
  TAxis* fPtAxis;
  std::vector<GFW::CorrConfig> corrconfigsTruth;
  std::vector<GFW::CorrConfig> corrconfigsReco;
  TRandom3* fRndm = new TRandom3(0);
  double epsilon = 1e-6;

  void init(InitContext&)
  {
    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    const AxisSpec axisVertex{20, -10, 10, "Vtxz (cm)"};
    const AxisSpec axisEta{20, -1., 1., "#eta"};
    const AxisSpec axisCounter{1, 0, +1, ""};
    // QA histograms
    histos.add<TH1>("mcEventCounter", "Monte Carlo Truth EventCounter", HistType::kTH1F, {axisCounter});
    histos.add<TH1>("numberOfRecoCollisions", "numberOfRecoCollisions", HistType::kTH1F, {{10, -0.5f, 9.5f}});
    histos.add<TH1>("RecoEventCounter", "Reconstruction EventCounter", HistType::kTH1F, {axisCounter});
    histos.add<TH1>("hnTPCClu", "Number of found TPC clusters", HistType::kTH1D, {{100, 40, 180}});
    histos.add<TH1>("hnITSClu", "Number of found ITS clusters", HistType::kTH1D, {{100, 0, 20}});
    // pT histograms
    histos.add<TH1>("hImpactParameter", "hImpactParameter", HistType::kTH1D, {axisB});
    histos.add<TH2>("hNchVsImpactParameter", "hNchVsImpactParameter", HistType::kTH2D, {axisB, axisNch});
    histos.add<TH1>("hEventPlaneAngle", "hEventPlaneAngle", HistType::kTH1D, {axisPhi});
    histos.add<TH2>("hPtVsPhiGenerated", "hPtVsPhiGenerated", HistType::kTH2D, {axisPhi, axisPt});
    histos.add<TH2>("hPtVsPhiGlobal", "hPtVsPhiGlobal", HistType::kTH2D, {axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiGenerated", "hBVsPtVsPhiGenerated", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiGlobal", "hBVsPtVsPhiGlobal", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiAny", "hBVsPtVsPhiAny", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiTPCTrack", "hBVsPtVsPhiTPCTrack", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiITSTrack", "hBVsPtVsPhiITSTrack", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiITSABTrack", "hBVsPtVsPhiITSABTrack", HistType::kTH3D, {axisB, axisPhi, axisPt});

    histos.add<TH3>("hBVsPtVsPhiGeneratedK0Short", "hBVsPtVsPhiGeneratedK0Short", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiGlobalK0Short", "hBVsPtVsPhiGlobalK0Short", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiGeneratedLambda", "hBVsPtVsPhiGeneratedLambda", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiGlobalLambda", "hBVsPtVsPhiGlobalLambda", HistType::kTH3D, {axisB, axisPhi, axisPt});

    histos.add<TH3>("hBVsPtVsPhiGeneratedXi", "hBVsPtVsPhiGeneratedXi", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiGlobalXi", "hBVsPtVsPhiGlobalXi", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiGeneratedOmega", "hBVsPtVsPhiGeneratedOmega", HistType::kTH3D, {axisB, axisPhi, axisPt});
    histos.add<TH3>("hBVsPtVsPhiGlobalOmega", "hBVsPtVsPhiGlobalOmega", HistType::kTH3D, {axisB, axisPhi, axisPt});

    histos.add<TH1>("hPhi", "#phi distribution", HistType::kTH1D, {axisPhi});
    histos.add<TH1>("hPhiWeighted", "corrected #phi distribution", HistType::kTH1D, {axisPhi});
    histos.add<TH2>("hEPVsPhiMC", "hEPVsPhiMC;Event Plane Angle; #varphi", HistType::kTH2D, {axisPhi, axisPhi});
    histos.add<TH2>("hEPVsPhi", "hEPVsPhi;Event Plane Angle; #varphi", HistType::kTH2D, {axisPhi, axisPhi});
    histos.add<TH2>("hPtNchGenerated", "Reco production; pT (GeV/c); multiplicity", HistType::kTH2D, {axisPt, axisNch});
    histos.add<TH2>("hPtNchGlobal", "Global production; pT (GeV/c); multiplicity", HistType::kTH2D, {axisPt, axisNch});
    histos.add<TH2>("hPtNchGeneratedPion", "Reco production; pT (GeV/c); multiplicity", HistType::kTH2D, {axisPt, axisNch});
    histos.add<TH2>("hPtNchGlobalPion", "Global production; pT (GeV/c); multiplicity", HistType::kTH2D, {axisPt, axisNch});
    histos.add<TH2>("hPtNchGeneratedKaon", "Reco production; pT (GeV/c); multiplicity", HistType::kTH2D, {axisPt, axisNch});
    histos.add<TH2>("hPtNchGlobalKaon", "Global production; pT (GeV/c); multiplicity", HistType::kTH2D, {axisPt, axisNch});
    histos.add<TH2>("hPtNchGeneratedProton", "Reco production; pT (GeV/c); multiplicity", HistType::kTH2D, {axisPt, axisNch});
    histos.add<TH2>("hPtNchGlobalProton", "Global production; pT (GeV/c); multiplicity", HistType::kTH2D, {axisPt, axisNch});
    histos.add<TH2>("hPtNchGeneratedK0", "Reco production; pT (GeV/c); multiplicity", HistType::kTH2D, {axisPt, axisNch});
    histos.add<TH2>("hPtNchGlobalK0", "Global production; pT (GeV/c); multiplicity", HistType::kTH2D, {axisPt, axisNch});
    histos.add<TH2>("hPtNchGeneratedLambda", "Reco production; pT (GeV/c); multiplicity", HistType::kTH2D, {axisPt, axisNch});
    histos.add<TH2>("hPtNchGlobalLambda", "Global production; pT (GeV/c); multiplicity", HistType::kTH2D, {axisPt, axisNch});
    histos.add<TH1>("hPtMCGen", "Monte Carlo Truth; pT (GeV/c);", {HistType::kTH1D, {axisPt}});
    histos.add<TH3>("hEtaPtVtxzMCGen", "Monte Carlo Truth; #eta; p_{T} (GeV/c); V_{z} (cm);", {HistType::kTH3D, {axisEta, axisPt, axisVertex}});
    histos.add<TH1>("hPtMCGlobal", "Monte Carlo Global; pT (GeV/c);", {HistType::kTH1D, {axisPt}});
    histos.add<TH3>("hEtaPtVtxzMCGlobal", "Monte Carlo Global; #eta; p_{T} (GeV/c); V_{z} (cm);", {HistType::kTH3D, {axisEta, axisPt, axisVertex}});
    histos.add<TH1>("hPhiWeightedTrDen", "corrected #phi distribution, considering track density", {HistType::kTH1D, {axisPhi}});

    o2::framework::AxisSpec axis = axisPt;
    int nPtBins = axis.binEdges.size() - 1;
    double* ptBins = &(axis.binEdges)[0];
    fPtAxis = new TAxis(nPtBins, ptBins);

    if (cfgOutputNUAWeights) {
      fWeights->setPtBins(nPtBins, ptBins);
      fWeights->init(true, false);
    }

    if (cfgFlowCumulantEnabled) {
      TObjArray* oba = new TObjArray();
      oba->Add(new TNamed("ChFull22", "ChFull22"));
      for (auto i = 0; i < fPtAxis->GetNbins(); i++)
        oba->Add(new TNamed(Form("ChFull22_pt_%i", i + 1), "ChFull22_pTDiff"));
      oba->Add(new TNamed("Ch10Gap22", "Ch10Gap22"));
      for (auto i = 0; i < fPtAxis->GetNbins(); i++)
        oba->Add(new TNamed(Form("Ch10Gap22_pt_%i", i + 1), "Ch10Gap22_pTDiff"));
      fFCTrue->SetName("FlowContainerTrue");
      fFCTrue->SetXAxis(fPtAxis);
      fFCTrue->Initialize(oba, axisCentrality, cfgFlowCumulantNbootstrap);
      fFCReco->SetName("FlowContainerReco");
      fFCReco->SetXAxis(fPtAxis);
      fFCReco->Initialize(oba, axisCentrality, cfgFlowCumulantNbootstrap);
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
    if (cfgCutDCAxy != 0.0) {
      myTrackSel.SetMaxDcaXY(cfgCutDCAxy);
    } else {
      fPtDepDCAxy = new TF1("ptDepDCAxy", Form("[0]*%s", cfgDCAxyFunction->c_str()), 0.001, 100);
      fPtDepDCAxy->SetParameter(0, cfgDCAxyNSigma);
      LOGF(info, "DCAxy pt-dependence function: %s", Form("[0]*%s", cfgDCAxyFunction->c_str()));
      myTrackSel.SetMaxDcaXYPtDep([fPtDepDCAxy = this->fPtDepDCAxy](float pt) { return fPtDepDCAxy->Eval(pt); });
    }
    myTrackSel.SetMinNClustersTPC(cfgCutTPCclu);
    myTrackSel.SetMinNCrossedRowsTPC(cfgCutTPCcrossedrows);
    if (cfgEnableITSCuts)
      myTrackSel.SetMinNClustersITS(cfgCutITSclu);
    if (!cfgCutDCAzPtDepEnabled)
      myTrackSel.SetMaxDcaZ(cfgCutDCAz);

    if (cfgTrackDensityCorrUse) {
      std::vector<double> pTEffBins = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0};
      hFindPtBin = new TH1D("hFindPtBin", "hFindPtBin", pTEffBins.size() - 1, &pTEffBins[0]);
      funcEff.resize(pTEffBins.size() - 1);
      // LHC24g3 Eff
      std::vector<double> f1p0 = cfgTrackDensityP0;
      std::vector<double> f1p1 = cfgTrackDensityP1;
      for (uint ifunc = 0; ifunc < pTEffBins.size() - 1; ifunc++) {
        funcEff[ifunc] = new TF1(Form("funcEff%i", ifunc), "[0]+[1]*x", 0, 3000);
        funcEff[ifunc]->SetParameters(f1p0[ifunc], f1p1[ifunc] * cfgTrackDensityCorrSlopeFactor);
      }
      funcV2 = new TF1("funcV2", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      funcV2->SetParameters(0.0186111, 0.00351907, -4.38264e-05, 1.35383e-07, -3.96266e-10);
      funcV3 = new TF1("funcV3", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      funcV3->SetParameters(0.0174056, 0.000703329, -1.45044e-05, 1.91991e-07, -1.62137e-09);
      funcV4 = new TF1("funcV4", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      funcV4->SetParameters(0.008845, 0.000259668, -3.24435e-06, 4.54837e-08, -6.01825e-10);
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

  void fillFC(GFW* fGFW, bool isMCTruth, const GFW::CorrConfig& corrconf, const double& cent, const double& rndm)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (!corrconf.pTDif) {
      if (dnx == 0)
        return;
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

  template <typename TCollision>
  bool eventSelected(TCollision collision)
  {
    if (std::fabs(collision.posZ()) > cfgCutVertex) {
      return 0;
    }
    if (cfgRecoEvSel8 && !collision.sel8()) {
      return 0;
    }
    if (cfgRecoEvkNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return 0;
    }
    if (cfgRecoEvkIsGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      // from Jan 9 2025 AOT meeting
      // cut time intervals with dead ITS staves
      return 0;
    }
    if (cfgRecoEvSelkIsGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return 0;
    }
    if (cfgRecoEvSelkNoITSROFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return 0;
    }
    if (cfgRecoEvSelkNoTimeFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return 0;
    }
    if (cfgRecoEvSelkNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      // no collisions in specified time range
      return 0;
    }
    if (cfgRecoEvSelkNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      // no other collisions in this Readout Frame with per-collision multiplicity above threshold
      return 0;
    }
    if (cfgRecoEvSelkNoHighMultCollInPrevRof && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      // veto an event if FT0C amplitude in previous ITS ROF is above threshold
      return 0;
    }
    return 1;
  }

  template <typename TTrack>
  bool trackSelected(TTrack track)
  {
    if (cfgCutDCAzPtDepEnabled && (track.dcaZ() > (0.004f + 0.013f / track.pt()))) {
      return false;
    }
    return myTrackSel.IsSelected(track);
  }

  void process(FilteredMcCollisions::iterator const& mcCollision, aod::BCsWithTimestamps const&, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels>> const& collisions, FilteredMcParticles const& mcParticles, FilteredTracks const&)
  {

    float imp = mcCollision.impactParameter();
    float evPhi = mcCollision.eventPlaneAngle();
    float vtxz = mcCollision.posZ();
    evPhi = RecoDecay::constrainAngle(evPhi);

    int64_t nCh = 0;
    int64_t nChGlobal = 0;
    float centrality = 0;
    float lRandom = fRndm->Rndm();
    float weff = 1.;
    float wacc = 1.;
    auto bc = mcCollision.bc_as<aod::BCsWithTimestamps>();
    loadCorrections(bc.timestamp());

    if (collisions.size() > -1) {
      histos.fill(HIST("mcEventCounter"), 0.5);
      histos.fill(HIST("numberOfRecoCollisions"), collisions.size()); // number of times coll was reco-ed
      if (cfgRecoEvRejectMC) {
        if (collisions.size() != 1) { // only pass those have one reconstruction event
          return;
        }
        for (auto const& collision : collisions) {
          if (!eventSelected(collision))
            return;
        }
      }
      histos.fill(HIST("RecoEventCounter"), 0.5);
    }

    if (imp > minB && imp < maxB) {
      // event within range
      histos.fill(HIST("hImpactParameter"), imp);
      histos.fill(HIST("hEventPlaneAngle"), evPhi);
      if (cfgFlowCumulantEnabled) {
        loadCentVsIPTruth(bc.timestamp());
        centrality = mCentVsIPTruth->GetBinContent(mCentVsIPTruth->GetXaxis()->FindBin(imp));
        fGFWTrue->Clear();
        fGFWReco->Clear();
      }

      double psi2Est = 0, psi3Est = 0, psi4Est = 0;
      float wEPeff = 1;
      double v2 = 0, v3 = 0, v4 = 0;
      double q2x = 0, q2y = 0;
      double q3x = 0, q3y = 0;
      double q4x = 0, q4y = 0;
      for (auto const& mcParticle : mcParticles) {
        int pdgCode = std::abs(mcParticle.pdgCode());
        if (pdgCode != PDG_t::kElectron && pdgCode != PDG_t::kMuonMinus && pdgCode != PDG_t::kPiPlus && pdgCode != kKPlus && pdgCode != PDG_t::kProton)
          continue;
        if (!mcParticle.isPhysicalPrimary())
          continue;
        if (std::fabs(mcParticle.eta()) > maxEta) // main acceptance
          continue;
        if (mcParticle.has_tracks()) {
          auto const& tracks = mcParticle.tracks_as<FilteredTracks>();
          for (auto const& track : tracks) {
            if (!trackSelected(track)) {
              continue;
            }
            if (cfgIsGlobalTrack && track.isGlobalTrack()) {
              nChGlobal++;
            }
            if (!cfgIsGlobalTrack && track.hasTPC() && track.hasITS()) {
              nChGlobal++;
            }
            if (cfgTrackDensityCorrUse && cfgFlowCumulantEnabled) {
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
          }
        }
      }
      if (cfgTrackDensityCorrUse && cfgFlowCumulantEnabled) {
        psi2Est = std::atan2(q2y, q2x) / 2.;
        psi3Est = std::atan2(q3y, q3x) / 3.;
        psi4Est = std::atan2(q4y, q4x) / 4.;
        v2 = funcV2->Eval(centrality);
        v3 = funcV3->Eval(centrality);
        v4 = funcV4->Eval(centrality);
      }

      for (auto const& mcParticle : mcParticles) {
        // focus on bulk: e, mu, pi, k, p
        int pdgCode = std::abs(mcParticle.pdgCode());
        bool extraPDGType = true;
        if (cfgK0Lambda0Enabled) {
          extraPDGType = (pdgCode != PDG_t::kK0Short && pdgCode != PDG_t::kLambda0);
        }
        if (extraPDGType && pdgCode != PDG_t::kElectron && pdgCode != PDG_t::kMuonMinus && pdgCode != PDG_t::kPiPlus && pdgCode != kKPlus && pdgCode != PDG_t::kProton)
          continue;

        if (!mcParticle.isPhysicalPrimary())
          continue;
        if (std::fabs(mcParticle.eta()) > maxEta) // main acceptance
          continue;

        float deltaPhi = mcParticle.phi() - mcCollision.eventPlaneAngle();
        deltaPhi = RecoDecay::constrainAngle(deltaPhi);
        histos.fill(HIST("hPtVsPhiGenerated"), deltaPhi, mcParticle.pt());
        histos.fill(HIST("hBVsPtVsPhiGenerated"), imp, deltaPhi, mcParticle.pt());
        histos.fill(HIST("hPtNchGenerated"), mcParticle.pt(), nChGlobal);
        histos.fill(HIST("hPtMCGen"), mcParticle.pt());
        histos.fill(HIST("hEtaPtVtxzMCGen"), mcParticle.eta(), mcParticle.pt(), vtxz);
        if (pdgCode == PDG_t::kPiPlus)
          histos.fill(HIST("hPtNchGeneratedPion"), mcParticle.pt(), nChGlobal);
        if (pdgCode == PDG_t::kKPlus)
          histos.fill(HIST("hPtNchGeneratedKaon"), mcParticle.pt(), nChGlobal);
        if (pdgCode == PDG_t::kProton)
          histos.fill(HIST("hPtNchGeneratedProton"), mcParticle.pt(), nChGlobal);
        if (pdgCode == PDG_t::kK0Short)
          histos.fill(HIST("hPtNchGeneratedK0"), mcParticle.pt(), nChGlobal);
        if (pdgCode == PDG_t::kLambda0)
          histos.fill(HIST("hPtNchGeneratedLambda"), mcParticle.pt(), nChGlobal);

        nCh++;

        bool validGlobal = false;
        bool validTrack = false;
        bool validTPCTrack = false;
        bool validITSTrack = false;
        bool validITSABTrack = false;
        if (mcParticle.has_tracks()) {
          auto const& tracks = mcParticle.tracks_as<FilteredTracks>();
          for (auto const& track : tracks) {
            if (!trackSelected(track)) {
              continue;
            }
            histos.fill(HIST("hnTPCClu"), track.tpcNClsFound());
            histos.fill(HIST("hnITSClu"), track.itsNCls());
            if (cfgIsGlobalTrack && track.isGlobalTrack()) {
              validGlobal = true;
            }
            if (!cfgIsGlobalTrack && track.hasTPC() && track.hasITS()) {
              validGlobal = true;
            }
            if (track.hasTPC() || track.hasITS()) {
              validTrack = true;
            }
            if (track.hasTPC()) {
              validTPCTrack = true;
            }
            if (track.hasITS() && track.itsChi2NCl() > -1. * epsilon) {
              validITSTrack = true;
            }
            if (track.hasITS() && track.itsChi2NCl() < -1. * epsilon) {
              validITSABTrack = true;
            }
          }
        }

        bool withinPtRef = (cfgCutPtRefMin < mcParticle.pt()) && (mcParticle.pt() < cfgCutPtRefMax); // within RF pT range
        bool withinPtPOI = (cfgCutPtPOIMin < mcParticle.pt()) && (mcParticle.pt() < cfgCutPtPOIMax); // within POI pT range
        if (cfgOutputNUAWeights && withinPtRef)
          fWeights->fill(mcParticle.phi(), mcParticle.eta(), vtxz, mcParticle.pt(), 0, 0);
        if (!setCurrentParticleWeights(weff, wacc, mcParticle.phi(), mcParticle.eta(), mcParticle.pt(), vtxz))
          continue;
        if (cfgTrackDensityCorrUse && cfgFlowCumulantEnabled && withinPtRef) {
          double fphi = v2 * std::cos(2 * (mcParticle.phi() - psi2Est)) + v3 * std::cos(3 * (mcParticle.phi() - psi3Est)) + v4 * std::cos(4 * (mcParticle.phi() - psi4Est));
          fphi = (1 + 2 * fphi);
          int pTBinForEff = hFindPtBin->FindBin(mcParticle.pt());
          if (pTBinForEff >= 1 && pTBinForEff <= hFindPtBin->GetNbinsX()) {
            wEPeff = funcEff[pTBinForEff - 1]->Eval(fphi * nChGlobal);
            if (wEPeff > 0.) {
              wEPeff = 1. / wEPeff;
              weff *= wEPeff;
              histos.fill(HIST("hPhiWeightedTrDen"), mcParticle.phi(), wacc * wEPeff);
            }
          }
        }

        if (cfgFlowCumulantEnabled) {
          if (withinPtRef)
            fGFWTrue->Fill(mcParticle.eta(), fPtAxis->FindBin(mcParticle.pt()) - 1, mcParticle.phi(), wacc * weff, 1);
          if (withinPtPOI)
            fGFWTrue->Fill(mcParticle.eta(), fPtAxis->FindBin(mcParticle.pt()) - 1, mcParticle.phi(), wacc * weff, 2);
          if (withinPtPOI && withinPtRef)
            fGFWTrue->Fill(mcParticle.eta(), fPtAxis->FindBin(mcParticle.pt()) - 1, mcParticle.phi(), wacc * weff, 4);

          if (validGlobal) {
            if (withinPtRef)
              fGFWReco->Fill(mcParticle.eta(), fPtAxis->FindBin(mcParticle.pt()) - 1, mcParticle.phi(), wacc * weff, 1);
            if (withinPtPOI)
              fGFWReco->Fill(mcParticle.eta(), fPtAxis->FindBin(mcParticle.pt()) - 1, mcParticle.phi(), wacc * weff, 2);
            if (withinPtPOI && withinPtRef)
              fGFWReco->Fill(mcParticle.eta(), fPtAxis->FindBin(mcParticle.pt()) - 1, mcParticle.phi(), wacc * weff, 4);
          }
        }

        if (withinPtRef) {
          histos.fill(HIST("hEPVsPhiMC"), evPhi, mcParticle.phi());
        }

        if (validGlobal && withinPtRef) {
          histos.fill(HIST("hPhi"), mcParticle.phi());
          histos.fill(HIST("hPhiWeighted"), mcParticle.phi(), wacc);
          histos.fill(HIST("hEPVsPhi"), evPhi, mcParticle.phi());
        }

        // if valid global, fill
        if (validGlobal) {
          histos.fill(HIST("hPtVsPhiGlobal"), deltaPhi, mcParticle.pt(), wacc * weff);
          histos.fill(HIST("hBVsPtVsPhiGlobal"), imp, deltaPhi, mcParticle.pt(), wacc * weff);
          histos.fill(HIST("hPtNchGlobal"), mcParticle.pt(), nChGlobal);
          histos.fill(HIST("hPtMCGlobal"), mcParticle.pt());
          histos.fill(HIST("hEtaPtVtxzMCGlobal"), mcParticle.eta(), mcParticle.pt(), vtxz);
          if (pdgCode == PDG_t::kPiPlus)
            histos.fill(HIST("hPtNchGlobalPion"), mcParticle.pt(), nChGlobal);
          if (pdgCode == PDG_t::kKPlus)
            histos.fill(HIST("hPtNchGlobalKaon"), mcParticle.pt(), nChGlobal);
          if (pdgCode == PDG_t::kProton)
            histos.fill(HIST("hPtNchGlobalProton"), mcParticle.pt(), nChGlobal);
          if (pdgCode == PDG_t::kK0Short)
            histos.fill(HIST("hPtNchGlobalK0"), mcParticle.pt(), nChGlobal);
          if (pdgCode == PDG_t::kLambda0)
            histos.fill(HIST("hPtNchGlobalLambda"), mcParticle.pt(), nChGlobal);
        }
        // if any track present, fill
        if (validTrack)
          histos.fill(HIST("hBVsPtVsPhiAny"), imp, deltaPhi, mcParticle.pt(), wacc * weff);
        if (validTPCTrack)
          histos.fill(HIST("hBVsPtVsPhiTPCTrack"), imp, deltaPhi, mcParticle.pt(), wacc * weff);
        if (validITSTrack)
          histos.fill(HIST("hBVsPtVsPhiITSTrack"), imp, deltaPhi, mcParticle.pt(), wacc * weff);
        if (validITSABTrack)
          histos.fill(HIST("hBVsPtVsPhiITSABTrack"), imp, deltaPhi, mcParticle.pt(), wacc * weff);
      }

      if (cfgFlowCumulantEnabled) {
        for (uint l_ind = 0; l_ind < corrconfigsTruth.size(); l_ind++) {
          fillFC(fGFWTrue, true, corrconfigsTruth.at(l_ind), centrality, lRandom);
        }
        for (uint l_ind = 0; l_ind < corrconfigsReco.size(); l_ind++) {
          fillFC(fGFWReco, false, corrconfigsReco.at(l_ind), centrality, lRandom);
        }
      }
    }
    histos.fill(HIST("hNchVsImpactParameter"), imp, nCh);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowMc>(cfgc)};
}
