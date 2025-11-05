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
/// \file eventMeanPtId.cxx
/// \brief Analysis task to study Mean pT Fluctuations using two particle correlator using Cumulant Method
/// \author Sweta Singh (sweta.singh@cern.ch)

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"

#include "TF1.h"
#include <TPDGCode.h>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

double massPi = o2::constants::physics::MassPionCharged;
double massKa = o2::constants::physics::MassKaonCharged;
double massPr = o2::constants::physics::MassProton;

using namespace o2::constants::physics;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;
using o2::constants::physics::Pdg;

struct EventMeanPtId {

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};
  // Configurables
  Configurable<float> cVtxZcut{"cVtxZcut", 10.0f, "Vertex Z"};
  Configurable<float> cEtacut{"cEtacut", 0.8f, "Eta cut"};
  Configurable<float> cPtmincut{"cPtmincut", 0.15f, "Pt min cut"};
  Configurable<float> cPtmaxcut{"cPtmaxcut", 2.0f, "Pt max cut"};
  Configurable<float> cPtmincut1{"cPtmincut1", 0.15f, " Pt min cut"};
  Configurable<float> cPtmaxcut1{"cPtmaxcut1", 2.0f, " Pt max cut"};
  Configurable<float> cDcaXYcut{"cDcaXYcut", 0.3f, "DCA XY cut"};
  Configurable<float> cDcaZcut{"cDcaZcut", 2.0f, "DCA Z cut"};
  Configurable<float> cCentmincut{"cCentmincut", 0.0, "Min cent cut"};
  Configurable<float> cCentmaxcut{"cCentmaxcut", 90.0, "Max cent cut"};
  Configurable<int> csyTPCcrosscut{"csyTPCcrosscut", 70, "TPC crossrows cut"};
  Configurable<int> csysItsChiCut{"csysItsChiCut", 36, "ITS chi2 cluster cut"};
  Configurable<int> csysTpcChiCut{"csysTpcChiCut", 4, "TPC chi2 cluster cut"};
  Configurable<int> csysnITSClustersCut{"csysnITSClustersCut", 5, "Number of ITS clusters cut"};
  Configurable<int> csystpcNClsCut{"csystpcNClsCut", 80, "No. of TPC clusters cut"};
  Configurable<double> threshold{"threshold", 1e-6, "Delta eta bin count"};
  Configurable<float> ptMax{"ptMax", 2.0, "maximum pT"};
  Configurable<float> ptMin{"ptMin", 0.15, "minimum pT"};
  Configurable<std::vector<double>> ptBins{"ptBins", {0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00}, "p_{T} bins"};
  // Event selections
  Configurable<bool> cSel8Trig{"cSel8Trig", true, "Sel8 (T0A + T0C) Selection Run3"};
  Configurable<bool> cpileupTFBorder{"cpileupTFBorder", true, "Timeframe Border Selection"};
  Configurable<bool> cpileupNoItsROBorder{"cpileupNoItsROBorder", true, "No ITSRO Border Cut"};
  Configurable<bool> cpileupItsTpcVtx{"cpileupItsTpcVtx", true, "ITS+TPC Vertex Selection"};
  Configurable<bool> cpileupSameBunch{"cpileupSameBunch", true, "Pileup rejection"};
  Configurable<bool> cpileupZVtxTimeDiff{"cpileupZVtxTimeDiff", true, "z-vtx time diff selection"};
  Configurable<bool> cIsGoodITSLayers{"cIsGoodITSLayers", true, "Good ITS Layers All"};
  Configurable<bool> cpileupItslayerall{"cpileupItslayerall", true, "dead staves of ITS removed"};
  Configurable<bool> cpileupvtxtofmatched{"cpileupvtxtofmatched", true, "TOF vertex matched"};
  Configurable<bool> citsNCluster{"citsNCluster", false, "Enable Number of ITS clusters"};
  Configurable<bool> ctpcNClusterFound{"ctpcNClusterFound", false, "Enable Number of TPC clusters"};
  Configurable<bool> cPVContributor{"cPVContributor", false, "Enable Primary Vertex Contributor"};
  Configurable<bool> csyDCAxy{"csyDCAxy", true, "DCAxy cut"};
  Configurable<bool> csyDCAz{"csyDCAz", true, "DCAz cut"};
  Configurable<bool> csyTPCcr{"csyTPCcr", true, "tpc crossed rows"};
  Configurable<bool> csyITSchi{"csyITSchi", true, "ITS chi2"};
  Configurable<bool> csyTPCchi{"csyTPCchi", true, "TPC chi2"};
  Configurable<bool> ccentFT0C{"ccentFT0C", true, "Use FT0C centraity"};
  Configurable<bool> pidSwitch{"pidSwitch", false, "pid calculations"};
  Configurable<bool> pidSwitchHistoFill{"pidSwitchHistoFill", false, "pid histogram filling"};
  Configurable<bool> effSwitch{"effSwitch", false, "efficiency calculations"};
  Configurable<bool> effSwitchHistoFill{"effSwitchHistoFill", false, "efficiency histogram filling"};
  // PID selection configurables
  Configurable<float> cpidPionPmincut{"cpidPionPmincut", 0.15, "pion min cut of pion"};
  Configurable<float> cpidKaonPmincut{"cpidKaonPmincut", 0.15, "kaon min cut of kaon"};
  Configurable<float> cpidProtonPmincut{"cpidProtonPmincut", 0.15, "proton min cut of proton"};
  Configurable<float> cpidPionPmaxcut{"cpidPionPmaxcut", 2.0, "pion min cut of pion"};
  Configurable<float> cpidKaonPmaxcut{"cpidKaonPmaxcut", 2.0, "kaon min cut of kaon"};
  Configurable<float> cpidProtonPmaxcut{"cpidProtonPmaxcut", 2.0, "proton min cut of proton"};
  Configurable<float> cpidPionPthcut{"cpidPionPthcut", 0.65, "pion threshold cut of pion"};
  Configurable<float> cpidKaonPthcut{"cpidKaonPthcut", 0.65, "kaon threshold cut of kaon"};
  Configurable<float> cpidProtonPthcut{"cpidProtonPthcut", 1.0, "proton threshold cut of proton"};
  Configurable<float> cNSigCut2{"cNSigCut2", 2.0, "nSigma cut (2)"};
  Configurable<float> cNSigCut3{"cNSigCut3", 3.0, "nSigma cut (3)"};
  Configurable<float> cElMinCut{"cElMinCut", -3.0, "electron min cut"};
  Configurable<float> cElMaxCut{"cElMaxCut", 5.0, "electron max cut"};
  Configurable<float> cTwoPtlCut2{"cTwoPtlCut2", 2.0, "n2ptl cut"};
  Configurable<float> cRapidityCut05{"cRapidityCut05", 0.5, "rapidity cut"};
  Configurable<int> nchBins{"nchBins", 4000, "Number of bins for nch axis"};
  Configurable<float> nchMin{"nchMin", 0.0, "Minimum value for nch axis"};
  Configurable<float> nchMax{"nchMax", 4000.0, "Maximum value for nch axis"};
  Configurable<float> cSigmaLowHighcut{"cSigmaLowHighcut", 3.0f, "lower and upper sigma cut"};
  O2_DEFINE_CONFIGURABLE(cfgEvSelMultCorrelation, bool, true, "Multiplicity correlation cut")
  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgMultPVFT0CCutEnabled, bool, true, "Enable PV multiplicity vs FT0C centrality cut")
    O2_DEFINE_CONFIGURABLE(cfgMultGlobalFT0CCutEnabled, bool, true, "Enable globalTracks vs FT0C centrality cut")
    O2_DEFINE_CONFIGURABLE(cfgMultGlobalPVCutEnabled, bool, true, "Enable globalTracks vs PV multiplicity cut")
    Configurable<std::vector<double>> cfgMultPVFT0CCutPars{"cfgMultPVFT0CCutPars",
                                                           std::vector<double>{3303.11, -121.316, 1.90207, -0.0152644, 5.10121e-05, 190.633, -4.32972, 0.0340001, -5.83261e-05, -3.19566e-07},
                                                           "PV multiplicity vs T0C centrality cut parameter values"};
    Configurable<std::vector<double>> cfgMultGlobalFT0CCutPars{"cfgMultGlobalFT0CCutPars",
                                                               std::vector<double>{1893.97, -61.3423, 0.790664, -0.00507208, 1.41683e-05, 167.997, -5.29125, 0.0840145, -0.000748102, 2.75743e-06},
                                                               "globalTracks vs FT0C cut parameter values"};
    Configurable<std::vector<double>> cfgMultGlobalPVCutPars{"cfgMultGlobalPVCutPars",
                                                             std::vector<double>{65.0322, 0.557725, -0.772828, 0.059224, -1.96379e-05, 4.46295e-09},
                                                             "globalTracks vs PV cut parameter values"};
    std::vector<double> multPVFT0CCutPars;
    std::vector<double> multGlobalFT0CPars;
    std::vector<double> multGlobalPVCutPars;
    TF1* fMultPVFT0CCutLow = nullptr;
    TF1* fMultPVFT0CCutHigh = nullptr;
    TF1* fMultGlobalFT0CCutLow = nullptr;
    TF1* fMultGlobalFT0CCutHigh = nullptr;
    TF1* fMultGlobalPVCutLow = nullptr;
    TF1* fMultGlobalPVCutHigh = nullptr;
  } cfgFunCoeff;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;

  Filter collisionFilter = nabs(aod::collision::posZ) <= cVtxZcut;
  Filter trackFilter = (nabs(aod::track::eta) < cEtacut) && (aod::track::pt > ptMin) && (aod::track::pt < ptMax) && (requireGlobalTrackInFilter());

  using MyCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFV0As>>;
  using MyCollision = MyCollisions::iterator;
  using MyTracks = soa::Filtered<soa::Join<aod::FullTracks,
                                           aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                           aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
                                           aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::StoredTracks,
                                           aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFbeta, aod::TOFSignal, aod::TracksExtra, aod::TracksIU, aod::TracksDCA, aod::TrackSelection>>;
  using MyTrack = MyTracks::iterator;
  using MyMCRecoCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::McCollisionLabels, aod::CentFT0Ms, aod::CentFT0As, aod::CentFV0As>>;
  using MyMCRecoCollision = MyMCRecoCollisions::iterator;
  using MyMCRecoTracks = soa::Filtered<soa::Join<aod::FullTracks,
                                                 aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                                 aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
                                                 aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::StoredTracks,
                                                 aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFbeta, aod::TOFSignal, aod::TracksExtra, aod::TracksIU, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>>;
  using MyMCRecoTrack = MyMCRecoTracks::iterator;
  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Cs, aod::CentFT0Ms, aod::Mults>;

  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> cfgUrlCCDB{"cfgUrlCCDB", "http://alice-ccdb.cern.ch", "url of ccdb"};
  Configurable<std::string> cfgPathCCDB{"cfgPathCCDB", "Users/s/swsingh/My/Object/GlobalRun3DCAcuts", "Path for ccdb-object"};
  Configurable<bool> cfgLoadEff{"cfgLoadEff", true, "Load efficiency"};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  TH2D* ptHistogramAllchargeRec = nullptr;

  void init(o2::framework::InitContext&)
  {
    if (cfgLoadEff) {
      // Set CCDB url
      ccdb->setURL(cfgUrlCCDB.value);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      // ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);
      // LOGF(info, "Getting object %s", ccdbPath.value.data());
      TList* lst = ccdb->getForTimeStamp<TList>(cfgPathCCDB.value, -1);
      ptHistogramAllchargeRec = reinterpret_cast<TH2D*>(lst->FindObject("hPtEta_rec"));
    }
    std::vector<double> ptBinning = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0};
    AxisSpec vtxZAxis = {100, -20.0, 20.0, "Z (cm)"};
    AxisSpec dcaAxis = {1002, -5.01, 5.01, "DCA_{xy} (cm)"};
    AxisSpec dcazAxis = {1002, -5.01, 5.01, "DCA_{z} (cm)"};
    AxisSpec ptAxis = {600, 0.0, 6.0, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pAxis = {400, 0.0, 4.0, "#it{p} (GeV/#it{c})"};
    AxisSpec betaAxis = {200, 0.0, 2.0, "TOF_{#beta} (GeV/#it{c})"};
    AxisSpec dEdxAxis = {2000, 0.0, 200.0, "dE/dx (GeV/#it{c})"};
    AxisSpec etaAxis = {300, -1.5, 1.5, "#eta"}; // 300, -1.5, 1.5
    AxisSpec nSigmaTPCAxis = {170, -8.5, 8.5, "n#sigma_{TPC}^{proton}"};
    AxisSpec nSigmaTPCAxispid = {170, -8.5, 8.5, "n#sigma_{TPC}"};
    AxisSpec nSigmaTOFAxispid = {170, -8.5, 8.5, "n#sigma_{TOF}"};
    AxisSpec centAxis = {100, 0., 100., "centrality"};
    AxisSpec subAxis = {30, 0., 30., "sample"};
    AxisSpec tnchAxis = {40, 0., 4000., "nch"};
    AxisSpec nchAxis = {nchBins, nchMin, nchMax, "nch"};
    AxisSpec varAxis1 = {400, 0., 4., "var1"};
    AxisSpec varAxis2 = {400, 0., 4., "var2"};
    AxisSpec tpcchi2Axis = {700, 0., 7., "tpc Chi2"};
    AxisSpec itschi2Axis = {400, 0., 40., "its Chi2"};
    AxisSpec crossedRowTpcAxis = {1600, 0., 160., "TPC Crossed rows"};
    AxisSpec counter = {10, 0., 10., "events"};
    // QA Plots
    histos.add("hEventcounter", "event counts", kTH1D, {counter});
    auto h = histos.add<TH1>("tracksel", "tracksel", HistType::kTH1D, {{15, 0.5, 15.5}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "Global track passed");
    h->GetXaxis()->SetBinLabel(3, "DCAxy passed");
    h->GetXaxis()->SetBinLabel(4, "DCAz passed");
    h->GetXaxis()->SetBinLabel(5, "Eta-cut passed");
    h->GetXaxis()->SetBinLabel(6, "pT-cut passed");
    h->GetXaxis()->SetBinLabel(7, "TPC crossed rows passed");
    h->GetXaxis()->SetBinLabel(8, "TPC Chai2cluster passed");
    h->GetXaxis()->SetBinLabel(9, "ITS Chai2cluster passed");
    h->GetXaxis()->SetBinLabel(10, "No. of ITS cluster 5 passed");
    h->GetXaxis()->SetBinLabel(11, "No. of TPC cluster 80 passed");
    auto hRec = histos.add<TH1>("trackSelRec", "trackSelRec", HistType::kTH1D, {{10, 0.5, 10.5}});
    hRec->GetXaxis()->SetBinLabel(1, "has_mcCollision() read");
    hRec->GetXaxis()->SetBinLabel(2, "Vertex Z > 10cm passed");
    hRec->GetXaxis()->SetBinLabel(3, "sel 8 passed");
    hRec->GetXaxis()->SetBinLabel(4, "kNoSameBunchPileup passed");
    hRec->GetXaxis()->SetBinLabel(5, "kNoITSROFrameBorder passed");
    hRec->GetXaxis()->SetBinLabel(6, "klsGoodZvtxFT0vsPV passed");
    hRec->GetXaxis()->SetBinLabel(7, "klsVertexITSTPC passed");
    histos.add("Data/hZvtx_before_sel", "hZvtx_before_sel", kTH1D, {vtxZAxis});
    histos.add("Data/hZvtx_after_sel8", "hZvtx_after_sel8", kTH1D, {vtxZAxis});
    histos.add("Data/hP", "hP", kTH1D, {pAxis});
    histos.add("Data/hEta", ";hEta", kTH1D, {etaAxis});
    histos.add("Data/hPt", ";#it{p}_{T} (GeV/#it{c})", kTH1D, {ptAxis});
    histos.add("Data/hPtvar", ";#it{p}_{T} (GeV/#it{c})", kTH1D, {ptAxis});
    histos.add("Data/hDCAxy", "hDCAxy", kTH1D, {dcaAxis});
    histos.add("Data/hDCAz", "hDCAz", kTH1D, {dcazAxis});
    histos.add("Data/hPtDCAxy", "hPtDCAxy", kTH2D, {ptAxis, dcaAxis});
    histos.add("Data/hPtDCAz", "hPtDCAz", kTH2D, {ptAxis, dcazAxis});
    histos.add("Data/hVar1", "hVar1", kTH2D, {subAxis, centAxis});
    histos.add("Data/hVar2", "hVar2", kTH2D, {subAxis, centAxis});
    histos.add("Data/hVar2meanpt", "hVar2meanpt", kTH2D, {centAxis, varAxis2});
    histos.add("Data/hVarc", "hVarc", kTH2D, {subAxis, centAxis});
    histos.add("Data/hnchAll", ";hnchAll", kTH1D, {nchAxis});
    histos.add("Data/hnchAll_bf_cut", ";hnchAll_bf_cut", kTH1D, {nchAxis});
    histos.add("Data/hnch", ";hnch", kTH1D, {nchAxis});
    histos.add("Data/hnchTrue", ";hnchTrue", kTH1D, {nchAxis});
    histos.add("Data/hnchTrue_pt", ";hnchTrue_pt", kTH1D, {nchAxis});
    histos.add("Data/hVar1x", "hVar1x", kTH2D, {subAxis, nchAxis});
    histos.add("Data/hVar2x", "hVar2x", kTH2D, {subAxis, nchAxis});
    histos.add("Data/hVarx", "hVarx", kTH2D, {subAxis, nchAxis});
    histos.add("Data/hdiffVar1x", "hdiffVar1x", kTH2D, {subAxis, nchAxis});
    histos.add("Data/hdiffVar2x", "hdiffVar2x", kTH2D, {subAxis, nchAxis});
    histos.add("Data/hdiffVarx", "hdiffVarx", kTH2D, {subAxis, nchAxis});
    histos.add("Data/hVar2meanptx", "hVar2meanptx", kTH2D, {nchAxis, varAxis2});
    histos.add("Data/hCentrality", "hCentrality", kTH1D, {centAxis});
    histos.add("Data/hPEta", "hPEta", kTH2D, {pAxis, etaAxis});
    histos.add("Data/hPtEta", "hPtEta", kTH2D, {ptAxis, etaAxis});
    histos.add("Data/hTPCchi2perCluster_before", "TPC #Chi^{2}/Cluster", kTH1D, {tpcchi2Axis});
    histos.add("Data/hITSchi2perCluster_before", "ITS #Chi^{2}/Cluster", kTH1D, {itschi2Axis});
    histos.add("Data/hTPCCrossedrows_before", "Crossed TPC rows", kTH1D, {crossedRowTpcAxis});
    histos.add("Data/hTPCchi2perCluster_after", "TPC #Chi^{2}/Cluster", kTH1D, {tpcchi2Axis});
    histos.add("Data/hITSchi2perCluster_after", "ITS #Chi^{2}/Cluster", kTH1D, {itschi2Axis});
    histos.add("Data/hTPCCrossedrows_after", "Crossed TPC rows", kTH1D, {crossedRowTpcAxis});
    histos.add("Data/hcent_nacc", "hcent_nacc", kTH2D, {centAxis, nchAxis});
    histos.add("Data/hcentFT0A_nacc", "hcentFT0A_nacc", kTH2D, {centAxis, nchAxis});
    histos.add("Data/hcentFT0M_nacc", "hcentFT0M_nacc", kTH2D, {centAxis, nchAxis});
    histos.add("Data/hcentFV0A_nacc", "hcentFV0A_nacc", kTH2D, {centAxis, nchAxis});
    histos.add("Data/hNchPV_NchGlobal_before", "hNchPV_NchGlobal_before", kTH2D, {nchAxis, nchAxis});
    histos.add("Data/hcentFT0C_GlobalNch_before", "hcentFT0C_GlobalNch_before", kTH2D, {centAxis, nchAxis});
    histos.add("Data/hcentFT0C_NchPV_before", "hcentFT0C_NchPV_before", kTH2D, {centAxis, nchAxis});
    histos.add("Data/hNchPV_NchGlobal_after", "hNchPV_NchGlobal_after", kTH2D, {nchAxis, nchAxis});
    histos.add("Data/hcentFT0C_GlobalNch_after", "hcentFT0C_GlobalNch_after", kTH2D, {centAxis, nchAxis});
    histos.add("Data/hcentFT0C_NchPV_after", "hcentFT0C_NchPV_after", kTH2D, {centAxis, nchAxis});
    histos.add("ptHistogramAllchargeRec", "ptHistogramAllchargeRec", kTH1D, {ptAxis});
    histos.add("hEta_rec", "", kTH1F, {etaAxis});
    histos.add("hPt_rec", "", kTH1F, {ptAxis});
    histos.add("hPtEta_rec", "hPtEta_rec", kTH2D, {ptAxis, etaAxis});
    histos.add("hNch_vs_Nch", "hNch_vs_Nch", kTH2D, {subAxis, nchAxis});
    histos.add("hterm1", "hterm1", kTProfile, {tnchAxis});
    histos.add("hterm2", "hterm2", kTProfile, {tnchAxis});
    histos.add("hCentrality_rec_before", "hCentrality_rec_before", kTH1D, {centAxis});
    histos.add("hEta1", ";hEta1", kTH1D, {etaAxis});
    if (effSwitchHistoFill) {
      histos.add("hEffVar1x_data", "hEffVar1x_data", kTH2D, {subAxis, nchAxis});
      histos.add("hEffVar2x_data", "hEffVar2x_data", kTH2D, {subAxis, nchAxis});
      histos.add("hEffVarx_data", "hEffVarx_data", kTH2D, {subAxis, nchAxis});
      histos.add("hEffVar2Meanptx_data", "hEffVar2Meanptx_data", kTH2D, {nchAxis, varAxis2});
      histos.add("hEffVar1x_Naccorr_data", "hEffVar1x_Naccorr_data", kTH2D, {subAxis, nchAxis});
      histos.add("hEffVar2x_Naccorr_data", "hEffVar2x_Naccorr_data", kTH2D, {subAxis, nchAxis});
      histos.add("hEffVarx_Naccorr_data", "hEffVarx_Naccorr_data", kTH2D, {subAxis, nchAxis});
      histos.add("hEffVar1x_Naccorr_xaxis_data", "hEffVar1x_Naccorr_xaxis_data", kTH2D, {subAxis, nchAxis});
      histos.add("hEffVar2x_Naccorr_xaxis_data", "hEffVar2x_Naccorr_xaxis_data", kTH2D, {subAxis, nchAxis});
      histos.add("hEta_rec_corr", "", kTH1F, {etaAxis});
      histos.add("hPt_rec_corr", "", kTH1F, {ptAxis});
      histos.add("hcent_nacc_corr", "hcent_nacc_corr", kTH2D, {centAxis, nchAxis});
      histos.add("hNch_vs_corr", "hNch_vs_corr", kTH2D, {subAxis, nchAxis});
    }
    if (pidSwitchHistoFill) {
      histos.add("Data/NSigamaTPCpion", "NSigamaTPCpion", kTH2D, {ptAxis, nSigmaTPCAxispid});
      histos.add("Data/NSigamaTPCkaon", "NSigamaTPCkaon", kTH2D, {ptAxis, nSigmaTPCAxispid});
      histos.add("Data/NSigamaTPCproton", "NSigamaTPCproton", kTH2D, {ptAxis, nSigmaTPCAxispid});
      histos.add("Data/NSigamaTOFpion", "NSigamaTOFpion", kTH2D, {ptAxis, nSigmaTOFAxispid});
      histos.add("Data/NSigamaTOFkaon", "NSigamaTOFkaon", kTH2D, {ptAxis, nSigmaTOFAxispid});
      histos.add("Data/NSigamaTOFproton", "NSigamaTOFproton", kTH2D, {ptAxis, nSigmaTOFAxispid});
      histos.add("Data/NSigamaTPCTOFpion", "NSigamaTPCTOFpion", kTH2D, {nSigmaTPCAxispid, nSigmaTOFAxispid});
      histos.add("Data/NSigamaTPCTOFkaon", "NSigamaTPCTOFkaon", kTH2D, {nSigmaTPCAxispid, nSigmaTOFAxispid});
      histos.add("Data/NSigamaTPCTOFproton", "NSigamaTPCTOFproton", kTH2D, {nSigmaTPCAxispid, nSigmaTOFAxispid});
      histos.add("Data/hPtPion", ";#it{p}_{T} (GeV/#it{c})", kTH1D, {ptAxis});
      histos.add("Data/hPtKaon", ";#it{p}_{T} (GeV/#it{c})", kTH1D, {ptAxis});
      histos.add("Data/hPtProton", ";#it{p}_{T} (GeV/#it{c})", kTH1D, {ptAxis});
      histos.add("Data/hEtaPion", ";hEta", kTH1D, {etaAxis});
      histos.add("Data/hEtaKaon", ";hEta", kTH1D, {etaAxis});
      histos.add("Data/hEtaProton", ";hEta", kTH1D, {etaAxis});
      histos.add("Data/hyPion", ";hyPion", kTH1D, {etaAxis});
      histos.add("Data/hyKaon", ";hyKaon", kTH1D, {etaAxis});
      histos.add("Data/hyProton", ";hyProton", kTH1D, {etaAxis});
      histos.add("Data/hVar1pix", "hVar1pix", kTH2D, {subAxis, nchAxis});
      histos.add("Data/hVar2pix", "hVar2pix", kTH2D, {subAxis, nchAxis});
      histos.add("Data/hVarpix", "hVarpix", kTH2D, {subAxis, nchAxis});
      histos.add("Data/hVar2meanptpix", "hVar2meanptpix", kTH2D, {nchAxis, varAxis2});
      histos.add("Data/hVar1kx", "hVar1kx", kTH2D, {subAxis, nchAxis});
      histos.add("Data/hVar2kx", "hVar2kx", kTH2D, {subAxis, nchAxis});
      histos.add("Data/hVarkx", "hVarkx", kTH2D, {subAxis, nchAxis});
      histos.add("Data/hVar2meanptkx", "hVar2meanptkx", kTH2D, {nchAxis, varAxis2});
      histos.add("Data/hVar1px", "hVar1px", kTH2D, {subAxis, nchAxis});
      histos.add("Data/hVar2px", "hVar2px", kTH2D, {subAxis, nchAxis});
      histos.add("Data/hVarpx", "hVarpx", kTH2D, {subAxis, nchAxis});
      histos.add("Data/hVar2meanptpx", "hVar2meanptpx", kTH2D, {nchAxis, varAxis2});
      histos.add("Data/hPtyPion", "hPtyPion", kTH2D, {ptAxis, etaAxis});
      histos.add("Data/hPtyKaon", "hPtyKaon", kTH2D, {ptAxis, etaAxis});
      histos.add("Data/hPtyProton", "hPtyProton", kTH2D, {ptAxis, etaAxis});
      histos.add("Data/hTOFbeta", "hTOFbeta", kTH2D, {pAxis, betaAxis});
      histos.add("Data/hdEdx", "hdEdx", kTH2D, {pAxis, dEdxAxis});
      histos.add("Data/hTOFbeta_afterselection", "hTOFbeta_afterselection", kTH2D, {pAxis, betaAxis});
      histos.add("Data/hdEdx_afterselection", "hdEdx_afterselection", kTH2D, {pAxis, dEdxAxis});
      histos.add("Data/hTOFbeta_afterselection1", "hTOFbeta_afterselection1", kTH2D, {pAxis, betaAxis});
      histos.add("Data/hdEdx_afterselection1", "hdEdx_afterselection1", kTH2D, {pAxis, dEdxAxis});
      histos.add("Data/hVar1pi", "hVar1pi", kTH2D, {subAxis, centAxis});
      histos.add("Data/hVar2pi", "hVar2pi", kTH2D, {subAxis, centAxis});
      histos.add("Data/hVar2meanptpi", "hVar2meanptpi", kTH2D, {centAxis, varAxis2});
      histos.add("Data/hVar1k", "hVar1k", kTH2D, {subAxis, centAxis});
      histos.add("Data/hVar2k", "hVar2k", kTH2D, {subAxis, centAxis});
      histos.add("Data/hVar2meanptk", "hVar2meanptk", kTH2D, {centAxis, varAxis2});
      histos.add("Data/hVar1p", "hVar1p", kTH2D, {subAxis, centAxis});
      histos.add("Data/hVar2p", "hVar2p", kTH2D, {subAxis, centAxis});
      histos.add("Data/hVarp", "hVarp", kTH2D, {subAxis, centAxis});
      histos.add("Data/hVar2meanptp", "hVar2meanptp", kTH2D, {centAxis, varAxis2});
      //===============reco level pid histograms==============================//
      histos.add("NSigamaTPCpion_rec", "NSigamaTPCpion_rec", kTH2D, {pAxis, nSigmaTPCAxispid});
      histos.add("NSigamaTPCkaon_rec", "NSigamaTPCkaon_rec", kTH2D, {pAxis, nSigmaTPCAxispid});
      histos.add("NSigamaTPCproton_rec", "NSigamaTPCproton_rec", kTH2D, {pAxis, nSigmaTPCAxispid});
      histos.add("NSigamaTOFpion_rec", "NSigamaTOFpion_rec", kTH2D, {pAxis, nSigmaTOFAxispid});
      histos.add("NSigamaTOFkaon_rec", "NSigamaTOFkaon_rec", kTH2D, {pAxis, nSigmaTOFAxispid});
      histos.add("NSigamaTOFproton_rec", "NSigamaTOFproton_rec", kTH2D, {pAxis, nSigmaTOFAxispid});
      histos.add("NSigamaTPCTOFpion_rec", "NSigamaTPCTOFpion_rec", kTH2D, {nSigmaTPCAxispid, nSigmaTOFAxispid});
      histos.add("NSigamaTPCTOFkaon_rec", "NSigamaTPCTOFkaon_rec", kTH2D, {nSigmaTPCAxispid, nSigmaTOFAxispid});
      histos.add("NSigamaTPCTOFproton_rec", "NSigamaTPCTOFproton_rec", kTH2D, {nSigmaTPCAxispid, nSigmaTOFAxispid});
      histos.add("hPtyPion_rec", "hPtyPion_rec", kTH2D, {ptAxis, etaAxis});
      histos.add("hPtyKaon_rec", "hPtyKaon_rec", kTH2D, {ptAxis, etaAxis});
      histos.add("hPtyProton_rec", "hPtyProton_rec", kTH2D, {ptAxis, etaAxis});
      histos.add("hPyPion_rec", "hPyPion_rec", kTH2D, {pAxis, etaAxis});
      histos.add("hPyKaon_rec", "hPyKaon_rec", kTH2D, {pAxis, etaAxis});
      histos.add("hPyProton_rec", "hPyProton_rec", kTH2D, {pAxis, etaAxis});
      histos.add("hTOFbeta_afterselection_rec_beforepidcut", "hTOFbeta_afterselection_rec_beforepidcut", kTH2D, {pAxis, betaAxis});
      histos.add("hdEdx_afterselection_rec_beforepidcut", "hdEdx_afterselection_rec_beforepidcut", kTH2D, {pAxis, dEdxAxis});
      histos.add("ptHistogramPionrec", "ptHistogramPionrec", kTH1D, {ptAxis});
      histos.add("ptHistogramKaonrec", "ptHistogramKaonrec", kTH1D, {ptAxis});
      histos.add("ptHistogramProtonrec", "ptHistogramProtonrec", kTH1D, {ptAxis});
      histos.add("ptHistogramPionrec_purity", "ptHistogramPionrec_purity", kTH1D, {ptAxis});
      histos.add("ptHistogramKaonrec_purity", "ptHistogramKaonrec_purity", kTH1D, {ptAxis});
      histos.add("ptHistogramProtonrec_purity", "ptHistogramProtonrec_purity", kTH1D, {ptAxis});
      histos.add("ptHistogramPionrec_pdg", "ptHistogramPionrec_pdg", kTH1D, {ptAxis});
      histos.add("ptHistogramKaonrec_pdg", "ptHistogramKaonrec_pdg", kTH1D, {ptAxis});
      histos.add("ptHistogramProtonrec_pdg", "ptHistogramProtonrec_pdg", kTH1D, {ptAxis});
      histos.add("hPtEta_pi_rec", "hPtEta_pi_rec", kTH2D, {ptAxis, etaAxis});
      histos.add("hPtEta_ka_rec", "hPtEta_ka_rec", kTH2D, {ptAxis, etaAxis});
      histos.add("hPtEta_pr_rec", "hPtEta_pr_rec", kTH2D, {ptAxis, etaAxis});
      //===============generated level pid histograms==============================//
      histos.add("hVar1pix_gen", "hVar1pix_gen", kTH2D, {subAxis, nchAxis});
      histos.add("hVar2pix_gen", "hVar2pix_gen", kTH2D, {subAxis, nchAxis});
      histos.add("hVarpix_gen", "hVarpix_gen", kTH2D, {subAxis, nchAxis});
      histos.add("hVar2meanptpix_gen", "hVar2meanptpix_gen", kTH2D, {nchAxis, varAxis2});
      histos.add("hVar1kx_gen", "hVar1kx_gen", kTH2D, {subAxis, nchAxis});
      histos.add("hVar2kx_gen", "hVar2kx_gen", kTH2D, {subAxis, nchAxis});
      histos.add("hVarkx_gen", "hVarkx_gen", kTH2D, {subAxis, nchAxis});
      histos.add("hVar2meanptkx_gen", "hVar2meanptkx_gen", kTH2D, {nchAxis, varAxis2});
      histos.add("hVar1px_gen", "hVar1px_gen", kTH2D, {subAxis, nchAxis});
      histos.add("hVar2px_gen", "hVar2px_gen", kTH2D, {subAxis, nchAxis});
      histos.add("hVarpx_gen", "hVarpx_gen", kTH2D, {subAxis, nchAxis});
      histos.add("hVar2meanptpx_gen", "hVar2meanptpx_gen", kTH2D, {nchAxis, varAxis2});
      histos.add("hPty_pi_gen", "hPty_pi_gen", kTH2D, {ptAxis, etaAxis});
      histos.add("hPty_ka_gen", "hPty_ka_gen", kTH2D, {ptAxis, etaAxis});
      histos.add("hPty_pr_gen", "hPty_pr_gen", kTH2D, {ptAxis, etaAxis});
      histos.add("hPtEta_pi_gen", "hPtEta_pi_gen", kTH2D, {ptAxis, etaAxis});
      histos.add("hPtEta_ka_gen", "hPtEta_ka_gen", kTH2D, {ptAxis, etaAxis});
      histos.add("hPtEta_pr_gen", "hPtEta_pr_gen", kTH2D, {ptAxis, etaAxis});
      histos.add("ptHistogramPion", "ptHistogramPion", kTH1D, {ptAxis});
      histos.add("ptHistogramKaon", "ptHistogramKaon", kTH1D, {ptAxis});
      histos.add("ptHistogramProton", "ptHistogramProton", kTH1D, {ptAxis});
      histos.add("hnch_pi", ";hnch_pi", kTH1D, {nchAxis});
      histos.add("hnch_ka", ";hnch_ka", kTH1D, {nchAxis});
      histos.add("hnch_pr", ";hnch_pr", kTH1D, {nchAxis});
    }
    //===============generated level allharons histograms==============================//
    histos.add("ptHistogram_allcharge_gen", "ptHistogram_allcharge_gen", kTH1D, {ptAxis});
    histos.add("hnch_gen_all", ";hnch_gen_all", kTH1D, {nchAxis});
    histos.add("hnch_gen_after_etacut", ";hnch_gen_after_etacut", kTH1D, {nchAxis});
    histos.add("hnch_afterPhysPrimary", ";hnch_afterPhysPrimary", kTH1D, {nchAxis});
    histos.add("hnch_gen", ";hnch_gen", kTH1D, {nchAxis});
    histos.add("hPtvar_gen", ";#it{p}_{T} (GeV/#it{c})", kTH1D, {ptAxis});
    histos.add("hVar1x_gen", "hVar1x_gen", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2x_gen", "hVar2x_gen", kTH2D, {subAxis, nchAxis});
    histos.add("hVarx_gen", "hVarx_gen", kTH2D, {subAxis, nchAxis});
    histos.add("hdiffVar1x_gen", "hdiffVar1x_gen", kTH2D, {subAxis, nchAxis});
    histos.add("hdiffVar2x_gen", "hdiffVar2x_gen", kTH2D, {subAxis, nchAxis});
    histos.add("hdiffVarx_gen", "hdiffVarx_gen", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2meanptx_gen", "hVar2meanptx_gen", kTH2D, {nchAxis, varAxis2});
    histos.add("hcent_nacc_gen", "hcent_nacc_gen", kTH2D, {centAxis, nchAxis});
    histos.add("hVtxZ_before_gen", "", kTH1F, {vtxZAxis});
    histos.add("hVtxZ_after_gensim", "", kTH1F, {vtxZAxis});
    histos.add("hEta_gen", "", kTH1F, {etaAxis});
    histos.add("hPt_gen", "", kTH1F, {ptAxis});
    histos.add("hPtEta_gen", "hPtEta_gen", kTH2D, {ptAxis, etaAxis});
    histos.add("hVar1_gen", "hVar1_gen", kTH2D, {subAxis, centAxis});
    histos.add("hVar2_gen", "hVar2_gen", kTH2D, {subAxis, centAxis});
    histos.add("hVarc_gen", "hVarc_gen", kTH2D, {subAxis, centAxis});
    histos.add("hterm1_gen", "hterm1_gen", kTProfile, {tnchAxis});
    histos.add("hterm2_gen", "hterm2_gen", kTProfile, {tnchAxis});
    cfgFunCoeff.multPVFT0CCutPars = cfgFunCoeff.cfgMultPVFT0CCutPars;
    cfgFunCoeff.multGlobalFT0CPars = cfgFunCoeff.cfgMultGlobalFT0CCutPars;
    cfgFunCoeff.multGlobalPVCutPars = cfgFunCoeff.cfgMultGlobalPVCutPars;
    cfgFunCoeff.fMultPVFT0CCutLow =
      new TF1("fMultPVFT0CCutLow",
              "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.0*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)",
              0, 100);
    cfgFunCoeff.fMultPVFT0CCutLow->SetParameters(&(cfgFunCoeff.multPVFT0CCutPars[0]));
    cfgFunCoeff.fMultPVFT0CCutHigh =
      new TF1("fMultPVFT0CCutHigh",
              "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.0*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)",
              0, 100);
    cfgFunCoeff.fMultPVFT0CCutHigh->SetParameters(&(cfgFunCoeff.multPVFT0CCutPars[0]));
    cfgFunCoeff.fMultGlobalFT0CCutLow =
      new TF1("fMultGlobalFT0CCutLow",
              "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.0*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)",
              0, 100);
    cfgFunCoeff.fMultGlobalFT0CCutLow->SetParameters(&(cfgFunCoeff.multGlobalFT0CPars[0]));
    cfgFunCoeff.fMultGlobalFT0CCutHigh =
      new TF1("fMultGlobalFT0CCutHigh",
              "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.0*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)",
              0, 100);
    cfgFunCoeff.fMultGlobalFT0CCutHigh->SetParameters(&(cfgFunCoeff.multGlobalFT0CPars[0]));
    cfgFunCoeff.fMultGlobalPVCutLow =
      new TF1("fMultGlobalPVCutLow",
              "[0] + [1]*x - 3.0*([2] + [3]*x + [4]*x*x + [5]*x*x*x)",
              0, 100);
    cfgFunCoeff.fMultGlobalPVCutLow->SetParameters(&(cfgFunCoeff.multGlobalPVCutPars[0]));
    cfgFunCoeff.fMultGlobalPVCutHigh =
      new TF1("fMultGlobalPVCutHigh",
              "[0] + [1]*x + 3.0*([2] + [3]*x + [4]*x*x + [5]*x*x*x)",
              0, 100);
    cfgFunCoeff.fMultGlobalPVCutHigh->SetParameters(&(cfgFunCoeff.multGlobalPVCutPars[0]));
    LOG(info) << "Printing Stored Registry Information";
    histos.print();
  }

  bool eventSelected(const float& globalNch, const float& pvTrack, const float& centrality)
  {
    if (cfgFunCoeff.cfgMultPVFT0CCutEnabled) {
      if (pvTrack < cfgFunCoeff.fMultPVFT0CCutLow->Eval(centrality))
        return false;
      if (pvTrack > cfgFunCoeff.fMultPVFT0CCutHigh->Eval(centrality))
        return false;
    }
    if (cfgFunCoeff.cfgMultGlobalFT0CCutEnabled) {
      if (globalNch < cfgFunCoeff.fMultGlobalFT0CCutLow->Eval(centrality))
        return false;
      if (globalNch > cfgFunCoeff.fMultGlobalFT0CCutHigh->Eval(centrality))
        return false;
    }
    if (cfgFunCoeff.cfgMultGlobalPVCutEnabled) {
      if (globalNch < cfgFunCoeff.fMultGlobalPVCutLow->Eval(pvTrack))
        return false;
      if (globalNch > cfgFunCoeff.fMultGlobalPVCutHigh->Eval(pvTrack))
        return false;
    }
    return true;
  }

  template <typename C>
  bool selCollision(C const& coll, float& cent)
  {
    if (std::abs(coll.posZ()) >= cVtxZcut) {
      return false;
    } // Reject the collisions with large vertex-z
    histos.fill(HIST("hEventcounter"), 2.);

    if (ccentFT0C) {
      cent = coll.centFT0C(); // centrality from FT0C
    } else {
      cent = coll.centFT0M(); // centrality from FT0M
    }

    if (cSel8Trig && !coll.sel8()) {
      return false;
    } // require min bias trigger
    histos.fill(HIST("hEventcounter"), 3.);

    if (cpileupTFBorder && !coll.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (cpileupNoItsROBorder && !coll.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    histos.fill(HIST("trackSelRec"), 4);

    if (cpileupSameBunch && !coll.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    histos.fill(HIST("trackSelRec"), 5);

    if (cpileupZVtxTimeDiff && !coll.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    histos.fill(HIST("trackSelRec"), 6);

    if (cpileupItsTpcVtx && !coll.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    histos.fill(HIST("trackSelRec"), 7);

    // if (cpileupItslayerall && !coll.selection_bit(aod::evsel::kIsGoodITSLayersAll))         {return false;}
    histos.fill(HIST("trackSelRec"), 8);

    if (cpileupvtxtofmatched && !coll.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    histos.fill(HIST("trackSelRec"), 9);

    return true; // if all checks pass, accept the collision
  }

  template <typename T>
  bool selTrack(T const& track)
  {
    if (!track.isGlobalTrack()) {
      return false;
    } // accept only global tracks
    histos.fill(HIST("tracksel"), 2);

    if (csyDCAxy && std::fabs(track.dcaXY()) > cDcaXYcut) {
      return false;
    }
    histos.fill(HIST("tracksel"), 3);

    if (csyDCAz && std::fabs(track.dcaZ()) > cDcaZcut) {
      return false;
    }
    histos.fill(HIST("tracksel"), 4);

    if (std::fabs(track.eta()) >= cEtacut) {
      return false;
    }
    histos.fill(HIST("tracksel"), 5);

    if (csyTPCcr && track.tpcNClsCrossedRows() < csyTPCcrosscut) {
      return false;
    }
    histos.fill(HIST("tracksel"), 6);

    if (csyITSchi && track.itsChi2NCl() >= csysItsChiCut) {
      return false;
    }
    histos.fill(HIST("tracksel"), 7);

    if (csyTPCchi && track.tpcChi2NCl() >= csysTpcChiCut) {
      return false;
    }
    histos.fill(HIST("tracksel"), 8);

    if (track.sign() == 0) {
      return false;
    }

    if (cPVContributor) {
      if (!(track.isPVContributor())) {
        return false;
      }
      histos.fill(HIST("tracksel"), 9);
    }

    if (citsNCluster) {
      if (track.itsNCls() < csysnITSClustersCut) {
        return false;
      }
      histos.fill(HIST("tracksel"), 10);
    }

    if (ctpcNClusterFound) {
      if (track.tpcNClsFound() < csystpcNClsCut) {
        return false;
      }
      histos.fill(HIST("tracksel"), 11);
    }

    return true; // if all checks pass, accept the collision
  }

  template <typename T>
  bool rejEl(T const& track)
  {
    if (track.tpcNSigmaEl() > cElMinCut && track.tpcNSigmaEl() < cElMaxCut && std::fabs(track.tpcNSigmaPi()) > cNSigCut3 && std::fabs(track.tpcNSigmaKa()) > cNSigCut3 && std::fabs(track.tpcNSigmaPr()) > cNSigCut3) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool selProton(T const& track)
  {
    //! if pt < threshold (For tracks without TOF information)
    if (track.p() > cpidProtonPmincut && track.p() <= cpidProtonPthcut) {
      if (track.hasTPC() && std::fabs(track.tpcNSigmaPr()) < cNSigCut2 && std::fabs(track.tpcNSigmaPi()) > cNSigCut2 && std::fabs(track.tpcNSigmaKa()) > cNSigCut2) {
        return true;
      }
    }

    //! if pt < threshold (For tracks with TOF information)
    if (track.p() > cpidProtonPmincut && track.p() <= cpidProtonPthcut) {
      if (track.hasTOF() && std::fabs(track.tpcNSigmaPr()) < cNSigCut2 && std::fabs(track.tofNSigmaPr()) < cNSigCut2 && std::fabs(track.tpcNSigmaPi()) > cNSigCut2 && std::fabs(track.tpcNSigmaKa()) > cNSigCut2) {
        return true;
      }
    }

    //! if pt > threshold (For tracks with TOF information)
    if (track.p() > cpidProtonPthcut && track.p() <= cpidProtonPmaxcut) {
      if (track.hasTPC() && track.hasTOF() && std::fabs(track.tpcNSigmaPr()) < cNSigCut2 && std::fabs(track.tofNSigmaPr()) < cNSigCut2 && std::hypot(track.tofNSigmaPi(), track.tpcNSigmaPi()) > cNSigCut2 && std::hypot(track.tofNSigmaKa(), track.tpcNSigmaKa()) > cNSigCut2) {
        return true;
      }
    }

    return false;
  }

  template <typename T>
  bool selKaon(T const& track)
  {
    //! if pt < threshold (For tracks without TOF information)
    if (track.p() > cpidKaonPmincut && track.p() <= cpidKaonPthcut) {
      if (track.hasTPC() && std::fabs(track.tpcNSigmaKa()) < cNSigCut2 && std::fabs(track.tpcNSigmaPi()) > cNSigCut2 && std::fabs(track.tpcNSigmaPr()) > cNSigCut2) {
        return true;
      }
    }

    //! if pt < threshold (For tracks with TOF information)
    if (track.p() > cpidKaonPmincut && track.p() <= cpidKaonPthcut) {
      if (track.hasTOF() && std::fabs(track.tpcNSigmaKa()) < cNSigCut2 && std::fabs(track.tofNSigmaKa()) < cNSigCut2 && std::fabs(track.tpcNSigmaPi()) > cNSigCut2 && std::fabs(track.tpcNSigmaPr()) > cNSigCut2) {
        return true;
      }
    }

    //! if pt > threshold (For tracks with TOF information)
    if (track.p() > cpidKaonPthcut && track.p() <= cpidKaonPmaxcut) {
      if (track.hasTPC() && track.hasTOF() && std::fabs(track.tpcNSigmaKa()) < cNSigCut2 && std::fabs(track.tofNSigmaKa()) < cNSigCut2 && std::hypot(track.tofNSigmaPi(), track.tpcNSigmaPi()) > cNSigCut2 && std::hypot(track.tofNSigmaPr(), track.tpcNSigmaPr()) > cNSigCut2) {
        return true;
      }
    }

    return false;
  }

  template <typename T>
  bool selPion(T const& track)
  {
    //! if pt < threshold (For tracks without TOF information)
    if (track.p() > cpidPionPmincut && track.p() <= cpidPionPthcut) {
      if (track.hasTPC() && std::fabs(track.tpcNSigmaPi()) < cNSigCut2 && std::fabs(track.tpcNSigmaKa()) > cNSigCut2 && std::fabs(track.tpcNSigmaPr()) > cNSigCut2) {
        return true;
      }
    }

    //! if pt < threshold (For tracks with TOF information)
    if (track.p() > cpidPionPmincut && track.p() <= cpidPionPthcut) {
      if (track.hasTOF() && std::fabs(track.tpcNSigmaPi()) < cNSigCut2 && std::fabs(track.tofNSigmaPi()) < cNSigCut2 && std::fabs(track.tpcNSigmaKa()) > cNSigCut2 && std::fabs(track.tpcNSigmaPr()) > cNSigCut2) {
        return true;
      }
    }

    //! if pt > threshold (For tracks with TOF information)
    if (track.p() > cpidPionPthcut && track.p() <= cpidPionPmaxcut) {
      if (track.hasTPC() && track.hasTOF() && std::fabs(track.tpcNSigmaPi()) < cNSigCut2 && std::fabs(track.tofNSigmaPi()) < cNSigCut2 && std::hypot(track.tofNSigmaKa(), track.tpcNSigmaKa()) > cNSigCut2 && std::hypot(track.tofNSigmaPr(), track.tpcNSigmaPr()) > cNSigCut2) {
        return true;
      }
    }

    return false;
  }

  double getEfficiency(double pt, double eta, TH2D* ptHistogramAllchargeRec)
  {
    int xbin = ptHistogramAllchargeRec->GetXaxis()->FindBin(pt);
    int ybin = ptHistogramAllchargeRec->GetYaxis()->FindBin(eta);

    if (xbin < 1 || xbin > ptHistogramAllchargeRec->GetNbinsX() || ybin < 1 || ybin > ptHistogramAllchargeRec->GetNbinsY()) {
      LOGF(warn, "pt or eta out of histograms bounds : %f, eta = %f", pt, eta);
      return 1e-6;
    }
    double eff = ptHistogramAllchargeRec->GetBinContent(xbin, ybin);
    return (eff > 0) ? eff : 1e-6; // Avoid division by zero
  }

  //++++++++++++++++++++++++++++++++++++DATA CALCULATION +++++++++++++++++++++++++++++++++++++++++++++++++++++//
  void processData(MyCollision const& coll, MyTracks const& inputTracks)

  {
    float cent = -1;
    histos.fill(HIST("hEventcounter"), 1.);
    histos.fill(HIST("Data/hZvtx_before_sel"), coll.posZ());

    if (!selCollision(coll, cent))
      return;
    {
      histos.fill(HIST("Data/hZvtx_after_sel8"), coll.posZ());
    }

    histos.fill(HIST("Data/hCentrality"), cent);

    float globalNch = inputTracks.size();
    float pvTrack = coll.multNTracksPV();

    histos.fill(HIST("Data/hNchPV_NchGlobal_before"), pvTrack, globalNch);
    histos.fill(HIST("Data/hcentFT0C_GlobalNch_before"), coll.centFT0C(), globalNch);
    histos.fill(HIST("Data/hcentFT0C_NchPV_before"), coll.centFT0C(), pvTrack);

    if (cfgEvSelMultCorrelation && !eventSelected(globalNch, pvTrack, cent)) {
      return;
    }

    histos.fill(HIST("Data/hNchPV_NchGlobal_after"), pvTrack, globalNch);
    histos.fill(HIST("Data/hcentFT0C_GlobalNch_after"), coll.centFT0C(), globalNch);
    histos.fill(HIST("Data/hcentFT0C_NchPV_after"), coll.centFT0C(), pvTrack);

    double nchAll = 0., nchAllBfCut = 0., nchEta = 0., nchPt = 0., nch = 0., nchPi = 0., nchKa = 0., nchPr = 0.;
    double q1 = 0., q2 = 0., var1 = 0., var2 = 0.;
    double sumPtWeight = 0., sumWeight = 0., sumPtPtWeight = 0., var1Eff = 0., var2Eff = 0.;
    double q1Pi = 0., q2Pi = 0., var1Pi = 0., var2Pi = 0.;
    double q1Ka = 0., q2Ka = 0., var1Ka = 0., var2Ka = 0.;
    double q1Pr = 0., q2Pr = 0., var1Pr = 0., var2Pr = 0.;

    int sample = histos.get<TH1>(HIST("Data/hZvtx_after_sel8"))->GetEntries();
    sample = sample % 30;

    for (const auto& track : inputTracks) {
      nchAllBfCut += 1.;
      histos.fill(HIST("Data/hnchAll_bf_cut"), nchAllBfCut);

      histos.fill(HIST("tracksel"), 1);
      histos.fill(HIST("Data/hTPCchi2perCluster_before"), track.tpcChi2NCl());
      histos.fill(HIST("Data/hITSchi2perCluster_before"), track.itsChi2NCl());
      histos.fill(HIST("Data/hTPCCrossedrows_before"), track.tpcNClsCrossedRows());

      if (std::fabs(track.eta()) <= cEtacut) {
        nchEta++;
        histos.fill(HIST("Data/hnchTrue"), nchEta);
      }
      if (track.pt() >= cPtmincut && track.pt() <= cPtmaxcut) {
        nchPt += 1.;
        histos.fill(HIST("Data/hnchTrue_pt"), nchPt);
      }
      if (!selTrack(track))
        continue;

      if (track.pt() >= cPtmincut1 && track.pt() <= cPtmaxcut1) {
        nch += 1.;
        histos.fill(HIST("Data/hnch"), nch);
        histos.fill(HIST("Data/hPtvar"), track.pt());
      }

      if (track.pt() < cPtmincut || track.pt() > cPtmaxcut)
        continue;

      nchAll += 1.;
      q1 += track.pt();
      q2 += (track.pt() * track.pt());

      histos.fill(HIST("Data/hnchAll"), nchAll);
      histos.fill(HIST("Data/hPt"), track.pt());
      histos.fill(HIST("Data/hEta"), track.eta());
      histos.fill(HIST("Data/hDCAxy"), track.dcaXY());
      histos.fill(HIST("Data/hDCAz"), track.dcaZ());
      histos.fill(HIST("Data/hTPCCrossedrows_after"), track.tpcNClsCrossedRows());
      histos.fill(HIST("Data/hTPCchi2perCluster_after"), track.tpcChi2NCl());
      histos.fill(HIST("Data/hITSchi2perCluster_after"), track.itsChi2NCl());
      histos.fill(HIST("Data/hP"), track.p());
      histos.fill(HIST("Data/hPtDCAxy"), track.pt(), track.dcaXY());
      histos.fill(HIST("Data/hPtDCAz"), track.pt(), track.dcaZ());
      histos.fill(HIST("Data/hPtEta"), track.pt(), track.eta());
      histos.fill(HIST("Data/hPEta"), track.p(), track.eta());

      if (effSwitch) {
        double eff = getEfficiency(track.pt(), track.eta(), ptHistogramAllchargeRec);
        if (eff < threshold)
          continue;
        double weight = 1. / eff;
        sumPtWeight += track.pt() / eff;
        sumPtPtWeight += (track.pt() * track.pt()) / (eff * eff);
        sumWeight += weight;
      }

      if (pidSwitch) {
        // only TPC tracks: Pion, Kaon, Proton
        if (track.hasTPC() && std::abs(track.tpcNSigmaPi()) < cNSigCut3)
          histos.fill(HIST("Data/NSigamaTPCpion"), track.pt(), track.tpcNSigmaPi());
        if (track.hasTPC() && std::abs(track.tpcNSigmaKa()) < cNSigCut3)
          histos.fill(HIST("Data/NSigamaTPCkaon"), track.pt(), track.tpcNSigmaKa());
        if (track.hasTPC() && std::abs(track.tpcNSigmaPr()) < cNSigCut3)
          histos.fill(HIST("Data/NSigamaTPCproton"), track.pt(), track.tpcNSigmaPr());

        // only TOF tracks: Pion, Kaon, Proton
        if (track.hasTOF() && std::abs(track.tofNSigmaPi()) < cNSigCut3)
          histos.fill(HIST("Data/NSigamaTOFpion"), track.pt(), track.tofNSigmaPi());
        if (track.hasTOF() && std::abs(track.tofNSigmaKa()) < cNSigCut3)
          histos.fill(HIST("Data/NSigamaTOFkaon"), track.pt(), track.tofNSigmaKa());
        if (track.hasTOF() && std::abs(track.tofNSigmaPr()) < cNSigCut3)
          histos.fill(HIST("Data/NSigamaTOFproton"), track.pt(), track.tofNSigmaPr());

        if (track.hasTPC())
          histos.fill(HIST("Data/hdEdx"), track.p(), track.tpcSignal());
        if (track.hasTOF())
          histos.fill(HIST("Data/hTOFbeta"), track.p(), track.beta());

        //===================================pion===========================================================
        // only TPC+TOF tracks: Pion, Kaon, Proton
        if ((track.hasTPC() && std::abs(track.tpcNSigmaPi()) < cNSigCut3) && (track.hasTOF() && std::abs(track.tofNSigmaPi()) < cNSigCut3)) {
          histos.fill(HIST("Data/NSigamaTPCTOFpion"), track.tpcNSigmaPi(), track.tofNSigmaPi());
          histos.fill(HIST("Data/hdEdx_afterselection"), track.p(), track.tpcSignal());
          histos.fill(HIST("Data/hTOFbeta_afterselection"), track.p(), track.beta());
        }
        if (selPion(track)) {
          histos.fill(HIST("Data/hPtPion"), track.pt());
          histos.fill(HIST("Data/hEtaPion"), track.eta());
          histos.fill(HIST("Data/hyPion"), track.rapidity(massPi));
          histos.fill(HIST("Data/hPtyPion"), track.pt(), track.rapidity(massPi));
          nchPi += 1.;
          q1Pi += track.pt();
          q2Pi += (track.pt() * track.pt());
          if (track.beta() > 1)
            continue;
          histos.fill(HIST("Data/hdEdx_afterselection1"), track.p(), track.tpcSignal());
          histos.fill(HIST("Data/hTOFbeta_afterselection1"), track.p(), track.beta());
        }

        //===========================kaon===============================================================
        if ((track.hasTPC() && std::abs(track.tpcNSigmaKa()) < cNSigCut3) && (track.hasTOF() && std::abs(track.tofNSigmaKa()) < cNSigCut3)) {
          histos.fill(HIST("Data/NSigamaTPCTOFkaon"), track.tpcNSigmaKa(), track.tofNSigmaKa());
          histos.fill(HIST("Data/hdEdx_afterselection"), track.p(), track.tpcSignal());
          histos.fill(HIST("Data/hTOFbeta_afterselection"), track.p(), track.beta());
        }
        if (selKaon(track)) {
          histos.fill(HIST("Data/hPtKaon"), track.pt());
          histos.fill(HIST("Data/hEtaKaon"), track.eta());
          histos.fill(HIST("Data/hyKaon"), track.rapidity(massKa));
          histos.fill(HIST("Data/hPtyKaon"), track.pt(), track.rapidity(massKa));
          nchKa += 1.;
          q1Ka += track.pt();
          q2Ka += (track.pt() * track.pt());
          if (track.beta() > 1)
            continue;
          histos.fill(HIST("Data/hdEdx_afterselection1"), track.p(), track.tpcSignal());
          histos.fill(HIST("Data/hTOFbeta_afterselection1"), track.p(), track.beta());
        }

        //============================proton===========================================================
        if ((track.hasTPC() && std::abs(track.tpcNSigmaPr()) < cNSigCut3) && (track.hasTOF() && std::abs(track.tofNSigmaPr()) < cNSigCut3)) {
          histos.fill(HIST("Data/NSigamaTPCTOFproton"), track.tpcNSigmaPr(), track.tofNSigmaPr());
          histos.fill(HIST("Data/hdEdx_afterselection"), track.p(), track.tpcSignal());
          histos.fill(HIST("Data/hTOFbeta_afterselection"), track.p(), track.beta());
        }
        if (selProton(track)) {
          histos.fill(HIST("Data/hPtProton"), track.pt());
          histos.fill(HIST("Data/hEtaProton"), track.eta());
          histos.fill(HIST("Data/hyProton"), track.rapidity(massPr));
          histos.fill(HIST("Data/hPtyProton"), track.pt(), track.rapidity(massPr));
          nchPr += 1.;
          q1Pr += track.pt();
          q2Pr += (track.pt() * track.pt());
          if (track.beta() > 1)
            continue;
          histos.fill(HIST("Data/hdEdx_afterselection1"), track.p(), track.tpcSignal());
          histos.fill(HIST("Data/hTOFbeta_afterselection1"), track.p(), track.beta());
        }
      }

    } // Track loop ends!
    histos.fill(HIST("Data/hcentFV0A_nacc"), coll.multFV0A(), nchAll);
    histos.fill(HIST("Data/hcentFT0A_nacc"), coll.multFT0A(), nchAll);
    histos.fill(HIST("Data/hcentFT0M_nacc"), coll.centFT0M(), nchAll);
    histos.fill(HIST("Data/hcent_nacc"), cent, nchAll);

    if (nchAll < cTwoPtlCut2)
      return;
    var1 = (q1 * q1 - q2) / (nchAll * (nchAll - 1));
    var2 = (q1 / nchAll);

    //---------------------- pions ----------------------------------------
    if (nchPi >= cTwoPtlCut2) {
      var1Pi = (q1Pi * q1Pi - q2Pi) / (nchPi * (nchPi - 1));
      var2Pi = (q1Pi / nchPi);
    }
    //----------------------- kaons ---------------------------------------
    if (nchKa >= cTwoPtlCut2) {
      var1Ka = (q1Ka * q1Ka - q2Ka) / (nchKa * (nchKa - 1));
      var2Ka = (q1Ka / nchKa);
    }
    //---------------------------- protons ----------------------------------
    if (nchPr >= cTwoPtlCut2) {
      var1Pr = (q1Pr * q1Pr - q2Pr) / (nchPr * (nchPr - 1));
      var2Pr = (q1Pr / nchPr);
    }

    //------------------ all charges-------------------------------------
    histos.fill(HIST("Data/hVar1"), sample, cent, var1);
    histos.fill(HIST("Data/hVar2"), sample, cent, var2);
    histos.fill(HIST("Data/hVarc"), sample, cent);
    histos.fill(HIST("Data/hVar2meanpt"), cent, var2);

    //-----------------------nch-------------------------------------
    histos.fill(HIST("Data/hVar1x"), sample, nchAll, var1);
    histos.fill(HIST("Data/hVar2x"), sample, nchAll, var2);
    histos.fill(HIST("Data/hVarx"), sample, nchAll);
    histos.fill(HIST("Data/hVar2meanptx"), nchAll, var2);
    histos.fill(HIST("Data/hdiffVar1x"), sample, nch, var1);
    histos.fill(HIST("Data/hdiffVar2x"), sample, nch, var2);
    histos.fill(HIST("Data/hdiffVarx"), sample, nch);

    if (effSwitchHistoFill) {
      //------------------ Efficiency corrected histograms ---------------
      var1Eff = (sumPtWeight * sumPtWeight - sumPtPtWeight) / (sumWeight * (sumWeight - 1));
      var2Eff = (sumPtWeight / sumWeight);

      histos.fill(HIST("hEffVar1x_data"), sample, nchAll, var1Eff);
      histos.fill(HIST("hEffVar2x_data"), sample, nchAll, var2Eff);
      histos.fill(HIST("hEffVarx_data"), sample, nchAll);
      histos.fill(HIST("hEffVar2Meanptx_data"), nchAll, var2Eff);
      histos.fill(HIST("hEffVar1x_Naccorr_data"), sample, sumWeight, var1Eff);
      histos.fill(HIST("hEffVar2x_Naccorr_data"), sample, sumWeight, var2Eff);
      histos.fill(HIST("hEffVarx_Naccorr_data"), sample, sumWeight);
      histos.fill(HIST("hEffVar1x_Naccorr_xaxis_data"), sample, sumWeight, var1);
      histos.fill(HIST("hEffVar2x_Naccorr_xaxis_data"), sample, sumWeight, var2);
    }

    if (pidSwitchHistoFill) {
      histos.fill(HIST("Data/hVar1pi"), sample, cent, var1Pi);
      histos.fill(HIST("Data/hVar2pi"), sample, cent, var2Pi);
      histos.fill(HIST("Data/hVar2meanptpi"), cent, var2Pi);
      histos.fill(HIST("Data/hVar1k"), sample, cent, var1Ka);
      histos.fill(HIST("Data/hVar2k"), sample, cent, var2Ka);
      histos.fill(HIST("Data/hVar2meanptk"), cent, var2Ka);
      histos.fill(HIST("Data/hVar1p"), sample, cent, var1Pr);
      histos.fill(HIST("Data/hVar2p"), sample, cent, var2Pr);
      histos.fill(HIST("Data/hVar2meanptp"), cent, var2Pr);
      histos.fill(HIST("Data/hVar1pix"), sample, nchAll, var1Pi);
      histos.fill(HIST("Data/hVar2pix"), sample, nchAll, var2Pi);
      histos.fill(HIST("Data/hVarpix"), sample, nchPi);
      histos.fill(HIST("Data/hVar2meanptpix"), nchAll, var2Pi);
      histos.fill(HIST("Data/hVar1kx"), sample, nchAll, var1Ka);
      histos.fill(HIST("Data/hVar2kx"), sample, nchAll, var2Ka);
      histos.fill(HIST("Data/hVarkx"), sample, nchKa);
      histos.fill(HIST("Data/hVar2meanptkx"), nchAll, var2Ka);
      histos.fill(HIST("Data/hVar1px"), sample, nchAll, var1Pr);
      histos.fill(HIST("Data/hVar2px"), sample, nchAll, var2Pr);
      histos.fill(HIST("Data/hVarpx"), sample, nchPr);
      histos.fill(HIST("Data/hVar2meanptpx"), nchAll, var2Pr);
    }

  } // event loop ends!

  PROCESS_SWITCH(EventMeanPtId, processData, "process real data information", true);

  //++++++++++++++++++++++++++++++++++++MC Reconstructed +++++++++++++++++++++++++++++++++++++++++++++++++++++//
  SliceCache cache;
  Preslice<aod::McParticles> mcTrack = o2::aod::mcparticle::mcCollisionId;
  void processMcReco(MyMCRecoCollision const& coll, MyMCRecoTracks const& inputTracks, aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    float cent = -1;
    (void)mcCollisions;
    if (!coll.has_mcCollision()) {
      return;
    }

    histos.fill(HIST("Data/hZvtx_before_sel"), coll.posZ());
    histos.fill(HIST("hVtxZ_before_gen"), coll.mcCollision().posZ());
    histos.fill(HIST("hCentrality_rec_before"), cent);

    if (!selCollision(coll, cent))
      return;

    histos.fill(HIST("Data/hZvtx_after_sel8"), coll.posZ());
    histos.fill(HIST("Data/hCentrality"), cent);

    float globalNch = inputTracks.size();
    float pvTrack = coll.multNTracksPV();

    histos.fill(HIST("Data/hNchPV_NchGlobal_before"), pvTrack, globalNch);
    histos.fill(HIST("Data/hcentFT0C_GlobalNch_before"), coll.centFT0C(), globalNch);
    histos.fill(HIST("Data/hcentFT0C_NchPV_before"), coll.centFT0C(), pvTrack);

    if (cfgEvSelMultCorrelation && !eventSelected(globalNch, pvTrack, cent)) {
      return;
    }

    histos.fill(HIST("Data/hNchPV_NchGlobal_after"), pvTrack, globalNch);
    histos.fill(HIST("Data/hcentFT0C_GlobalNch_after"), coll.centFT0C(), globalNch);
    histos.fill(HIST("Data/hcentFT0C_NchPV_after"), coll.centFT0C(), pvTrack);

    double nch = 0., nchPi = 0., nchKa = 0., nchPr = 0., nchAll = 0., nchAllBfCut = 0., nchEta = 0., nchPt = 0.;
    double q1 = 0., q2 = 0.;
    double q1Pi = 0., q2Pi = 0., q1Ka = 0., q2Ka = 0., q1Pr = 0., q2Pr = 0.;
    double var1 = 0., var2 = 0.;
    double var1Pi = 0., var2Pi = 0., var1Ka = 0., var2Ka = 0., var1Pr = 0., var2Pr = 0.;
    double sumPtWeight = 0., sumWeight = 0., sumPtPtWeight = 0., var1Eff = 0., var2Eff = 0.;

    int sample = histos.get<TH1>(HIST("Data/hZvtx_after_sel8"))->GetEntries();
    sample = sample % 30;

    for (const auto& track : inputTracks) {
      nchAllBfCut += 1.;
      histos.fill(HIST("Data/hnchAll_bf_cut"), nchAllBfCut);
      histos.fill(HIST("Data/hTPCchi2perCluster_before"), track.tpcChi2NCl());
      histos.fill(HIST("Data/hITSchi2perCluster_before"), track.itsChi2NCl());
      histos.fill(HIST("Data/hTPCCrossedrows_before"), track.tpcNClsCrossedRows());

      if (std::fabs(track.eta()) <= cEtacut) {
        nchEta++;
        histos.fill(HIST("Data/hnchTrue"), nchEta);
      }
      if (track.pt() >= cPtmincut && track.pt() <= cPtmaxcut) {
        nchPt += 1.;
        histos.fill(HIST("Data/hnchTrue_pt"), nchPt);
      }

      if (!selTrack(track))
        continue;

      if (track.pt() >= cPtmincut1 && track.pt() <= cPtmaxcut1) {
        nch += 1.;
        histos.fill(HIST("Data/hnch"), nch);
        histos.fill(HIST("Data/hPtvar"), track.pt());
      }
      if (track.pt() < cPtmincut || track.pt() > cPtmaxcut)
        continue;

      // if (std::fabs(track.y()) > 0.5) continue;
      histos.fill(HIST("hPt_rec"), track.pt());
      histos.fill(HIST("hEta_rec"), track.eta());

      if (effSwitch) {
        double eff = getEfficiency(track.pt(), track.eta(), ptHistogramAllchargeRec);
        if (eff < threshold)
          continue;
        double weight = 1.0 / eff;
        sumPtWeight += track.pt() * weight;
        sumPtPtWeight += (track.pt() * track.pt() * weight * weight);
        sumWeight += weight;

        histos.fill(HIST("hPt_rec_corr"), track.pt(), weight);
        histos.fill(HIST("hEta_rec_corr"), track.eta(), weight);
      }

      auto mcParticle = track.mcParticle();
      nchAll += 1.;
      q1 += track.pt();
      q2 += (track.pt() * track.pt());

      histos.fill(HIST("Data/hnchAll"), nchAll);
      histos.fill(HIST("ptHistogramAllchargeRec"), track.pt());
      histos.fill(HIST("Data/hDCAxy"), track.dcaXY());
      histos.fill(HIST("Data/hDCAz"), track.dcaZ());
      histos.fill(HIST("Data/hTPCCrossedrows_after"), track.tpcNClsCrossedRows());
      histos.fill(HIST("Data/hTPCchi2perCluster_after"), track.tpcChi2NCl());
      histos.fill(HIST("Data/hITSchi2perCluster_after"), track.itsChi2NCl());
      histos.fill(HIST("Data/hP"), track.p());
      histos.fill(HIST("Data/hPt"), track.pt());
      histos.fill(HIST("Data/hEta"), track.eta());
      histos.fill(HIST("Data/hPtDCAxy"), track.pt(), track.dcaXY());
      histos.fill(HIST("Data/hPtDCAz"), track.pt(), track.dcaZ());
      histos.fill(HIST("Data/hPtEta"), track.pt(), track.eta());
      histos.fill(HIST("Data/hPEta"), track.p(), track.eta());
      histos.fill(HIST("hPtEta_rec"), track.pt(), track.eta());

      if (pidSwitch) {
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus)
          histos.fill(HIST("ptHistogramPionrec_pdg"), track.pt());
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus)
          histos.fill(HIST("ptHistogramKaonrec_pdg"), track.pt());
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kProton)
          histos.fill(HIST("ptHistogramProtonrec_pdg"), track.pt());

        // only TPC tracks: Pion, Kaon, Proton
        if (track.hasTPC() && std::abs(track.tpcNSigmaPi()) < cNSigCut3)
          histos.fill(HIST("Data/NSigamaTPCpion"), track.pt(), track.tpcNSigmaPi());
        if (track.hasTPC() && std::abs(track.tpcNSigmaKa()) < cNSigCut3)
          histos.fill(HIST("Data/NSigamaTPCkaon"), track.pt(), track.tpcNSigmaKa());
        if (track.hasTPC() && std::abs(track.tpcNSigmaPr()) < cNSigCut3)
          histos.fill(HIST("Data/NSigamaTPCproton"), track.pt(), track.tpcNSigmaPr());

        // only TOF tracks: Pion, Kaon, Proton
        if (track.hasTOF() && std::abs(track.tofNSigmaPi()) < cNSigCut3)
          histos.fill(HIST("Data/NSigamaTOFpion"), track.pt(), track.tofNSigmaPi());
        if (track.hasTOF() && std::abs(track.tofNSigmaKa()) < cNSigCut3)
          histos.fill(HIST("Data/NSigamaTOFkaon"), track.pt(), track.tofNSigmaKa());
        if (track.hasTOF() && std::abs(track.tofNSigmaPr()) < cNSigCut3)
          histos.fill(HIST("Data/NSigamaTOFproton"), track.pt(), track.tofNSigmaPr());

        if (track.hasTPC())
          histos.fill(HIST("Data/hdEdx"), track.p(), track.tpcSignal());
        if (track.hasTOF())
          histos.fill(HIST("Data/hTOFbeta"), track.p(), track.beta());
        if (track.hasTPC())
          histos.fill(HIST("hdEdx_afterselection_rec_beforepidcut"), track.p(), track.tpcSignal());
        if (track.hasTOF())
          histos.fill(HIST("hTOFbeta_afterselection_rec_beforepidcut"), track.p(), track.beta());

        //===================================pion==============================================================
        if ((track.hasTPC() && std::abs(track.tpcNSigmaPi()) < cNSigCut3) && (track.hasTOF() && std::abs(track.tofNSigmaPi()) < cNSigCut3)) {
          histos.fill(HIST("Data/NSigamaTPCTOFpion"), track.tpcNSigmaPi(), track.tofNSigmaPi());

          histos.fill(HIST("Data/hdEdx_afterselection"), track.p(), track.tpcSignal());
          histos.fill(HIST("Data/hTOFbeta_afterselection"), track.p(), track.beta());
        }

        if (selPion(track)) {
          if (std::fabs(track.y()) > cRapidityCut05)
            continue;
          if (track.beta() > 1)
            continue;
          histos.fill(HIST("ptHistogramPionrec"), track.pt());
          histos.fill(HIST("hPtEta_pi_rec"), track.pt(), track.eta());
          histos.fill(HIST("Data/hPtPion"), track.pt());
          histos.fill(HIST("Data/hEtaPion"), track.eta());
          histos.fill(HIST("Data/hyPion"), track.rapidity(massPi));
          histos.fill(HIST("Data/hPtyPion"), track.pt(), track.rapidity(massPi));
          histos.fill(HIST("NSigamaTPCpion_rec"), track.p(), track.tpcNSigmaPi());
          histos.fill(HIST("NSigamaTOFpion_rec"), track.p(), track.tofNSigmaPi());
          histos.fill(HIST("NSigamaTPCTOFpion_rec"), track.tpcNSigmaPi(), track.tofNSigmaPi());
          histos.fill(HIST("Data/hdEdx_afterselection1"), track.p(), track.tpcSignal());
          histos.fill(HIST("Data/hTOFbeta_afterselection1"), track.p(), track.beta());
          if (std::abs(track.mcParticle().pdgCode()) == PDG_t::kPiPlus) {
            histos.fill(HIST("ptHistogramPionrec_purity"), track.pt());
          }
          nchPi += 1.;
          q1Pi += track.pt();
          q2Pi += (track.pt() * track.pt());

          histos.fill(HIST("hPyPion_rec"), track.p(), track.rapidity(massPi));
          histos.fill(HIST("hPtyPion_rec"), track.pt(), track.rapidity(massPi));
        }

        //===========================kaon===============================================================

        if ((track.hasTPC() && std::abs(track.tpcNSigmaKa()) < cNSigCut3) && (track.hasTOF() && std::abs(track.tofNSigmaKa()) < cNSigCut3)) {
          histos.fill(HIST("Data/NSigamaTPCTOFkaon"), track.tpcNSigmaKa(), track.tofNSigmaKa());
          histos.fill(HIST("Data/hdEdx_afterselection"), track.p(), track.tpcSignal());
          histos.fill(HIST("Data/hTOFbeta_afterselection"), track.p(), track.beta());
        }

        if (selKaon(track)) {
          if (std::fabs(track.y()) > cRapidityCut05)
            continue;
          if (track.beta() > 1)
            continue;
          histos.fill(HIST("ptHistogramKaonrec"), track.pt());
          histos.fill(HIST("hPtEta_ka_rec"), track.pt(), track.eta());
          histos.fill(HIST("Data/hPtKaon"), track.pt());
          histos.fill(HIST("Data/hEtaKaon"), track.eta());
          histos.fill(HIST("Data/hyKaon"), track.rapidity(massKa));
          histos.fill(HIST("Data/hPtyKaon"), track.pt(), track.rapidity(massKa));
          histos.fill(HIST("NSigamaTPCkaon_rec"), track.p(), track.tpcNSigmaKa());
          histos.fill(HIST("NSigamaTOFkaon_rec"), track.p(), track.tofNSigmaKa());
          histos.fill(HIST("NSigamaTPCTOFkaon_rec"), track.tpcNSigmaKa(), track.tofNSigmaKa());
          histos.fill(HIST("Data/hdEdx_afterselection1"), track.p(), track.tpcSignal());
          histos.fill(HIST("Data/hTOFbeta_afterselection1"), track.p(), track.beta());
          if (std::abs(track.mcParticle().pdgCode()) == PDG_t::kKPlus) {
            histos.fill(HIST("ptHistogramKaonrec_purity"), track.pt());
          }
          nchKa += 1.;
          q1Ka += track.pt();
          q2Ka += (track.pt() * track.pt());

          histos.fill(HIST("hPyKaon_rec"), track.p(), track.rapidity(massKa));
          histos.fill(HIST("hPtyKaon_rec"), track.pt(), track.rapidity(massKa));
        }

        //============================proton===========================================================

        if ((track.hasTPC() && std::abs(track.tpcNSigmaPr()) < cNSigCut3) && (track.hasTOF() && std::abs(track.tofNSigmaPr()) < cNSigCut3)) {
          histos.fill(HIST("Data/NSigamaTPCTOFproton"), track.tpcNSigmaPr(), track.tofNSigmaPr());
          histos.fill(HIST("Data/hdEdx_afterselection"), track.p(), track.tpcSignal());
          histos.fill(HIST("Data/hTOFbeta_afterselection"), track.p(), track.beta());
        }

        if (selProton(track)) {
          if (std::fabs(track.y()) > cRapidityCut05)
            continue;
          if (track.beta() > 1)
            continue;
          histos.fill(HIST("ptHistogramProtonrec"), track.pt());
          histos.fill(HIST("hPtEta_pr_rec"), track.pt(), track.eta());
          histos.fill(HIST("Data/hPtProton"), track.pt());
          histos.fill(HIST("Data/hEtaProton"), track.eta());
          histos.fill(HIST("Data/hyProton"), track.rapidity(massPr));
          histos.fill(HIST("Data/hPtyProton"), track.pt(), track.rapidity(massPr));
          histos.fill(HIST("NSigamaTPCproton_rec"), track.p(), track.tpcNSigmaPr());
          histos.fill(HIST("NSigamaTOFproton_rec"), track.p(), track.tofNSigmaPr());
          histos.fill(HIST("NSigamaTPCTOFproton_rec"), track.tpcNSigmaPr(), track.tofNSigmaPr());
          histos.fill(HIST("Data/hdEdx_afterselection1"), track.p(), track.tpcSignal());
          histos.fill(HIST("Data/hTOFbeta_afterselection1"), track.p(), track.beta());
          if (std::abs(track.mcParticle().pdgCode()) == PDG_t::kProton) {
            histos.fill(HIST("ptHistogramProtonrec_purity"), track.pt());
          }
          nchPr += 1.;
          q1Pr += track.pt();
          q2Pr += (track.pt() * track.pt());

          histos.fill(HIST("hPyProton_rec"), track.p(), track.rapidity(massPr));
          histos.fill(HIST("hPtyProton_rec"), track.pt(), track.rapidity(massPr));
        }
      }

    } // loop over tracks
    histos.fill(HIST("Data/hcent_nacc"), cent, nchAll);
    histos.fill(HIST("hNch_vs_Nch"), sample, nchAll, nchAll);

    if (nchAll < cTwoPtlCut2)
      return;
    var1 = (q1 * q1 - q2) / (nchAll * (nchAll - 1));
    var2 = (q1 / nchAll);

    histos.fill(HIST("hterm1"), nchAll, var1);
    histos.fill(HIST("hterm2"), nchAll, var2);

    histos.fill(HIST("Data/hVar1"), sample, cent, var1);
    histos.fill(HIST("Data/hVar2"), sample, cent, var2);
    histos.fill(HIST("Data/hVarc"), sample, cent);
    histos.fill(HIST("Data/hVar2meanpt"), cent, var2);

    //---------------------- pions ----------------------------------------
    if (nchPi >= cTwoPtlCut2) {
      var1Pi = (q1Pi * q1Pi - q2Pi) / (nchPi * (nchPi - 1));
      var2Pi = (q1Pi / nchPi);
    }
    //----------------------- kaons ---------------------------------------
    if (nchKa >= cTwoPtlCut2) {
      var1Ka = (q1Ka * q1Ka - q2Ka) / (nchKa * (nchKa - 1));
      var2Ka = (q1Ka / nchKa);
    }
    //---------------------------- protons ----------------------------------
    if (nchPr >= cTwoPtlCut2) {
      var1Pr = (q1Pr * q1Pr - q2Pr) / (nchPr * (nchPr - 1));
      var2Pr = (q1Pr / nchPr);
    }

    //-----------------------nch-------------------------------------
    histos.fill(HIST("Data/hVar1x"), sample, nchAll, var1);
    histos.fill(HIST("Data/hVar2x"), sample, nchAll, var2);
    histos.fill(HIST("Data/hVarx"), sample, nchAll);
    histos.fill(HIST("Data/hdiffVar1x"), sample, nch, var1);
    histos.fill(HIST("Data/hdiffVar2x"), sample, nch, var2);
    histos.fill(HIST("Data/hdiffVarx"), sample, nch);
    histos.fill(HIST("Data/hVar2meanptx"), nchAll, var2);

    if (pidSwitchHistoFill) {
      histos.fill(HIST("Data/hVar1pi"), sample, cent, var1Pi);
      histos.fill(HIST("Data/hVar2pi"), sample, cent, var2Pi);
      histos.fill(HIST("Data/hVar2meanptpi"), cent, var2Pi);
      histos.fill(HIST("Data/hVar1k"), sample, cent, var1Ka);
      histos.fill(HIST("Data/hVar2k"), sample, cent, var2Ka);
      histos.fill(HIST("Data/hVar2meanptk"), cent, var2Ka);
      histos.fill(HIST("Data/hVar1p"), sample, cent, var1Pr);
      histos.fill(HIST("Data/hVar2p"), sample, cent, var2Pr);
      histos.fill(HIST("Data/hVar2meanptp"), cent, var2Pr);
      histos.fill(HIST("Data/hVar1pix"), sample, nchAll, var1Pi);
      histos.fill(HIST("Data/hVar2pix"), sample, nchAll, var2Pi);
      histos.fill(HIST("Data/hVarpix"), sample, nchPi);
      histos.fill(HIST("Data/hVar2meanptpix"), nchAll, var2Pi);
      histos.fill(HIST("Data/hVar1kx"), sample, nchAll, var1Ka);
      histos.fill(HIST("Data/hVar2kx"), sample, nchAll, var2Ka);
      histos.fill(HIST("Data/hVarkx"), sample, nchKa);
      histos.fill(HIST("Data/hVar2meanptkx"), nchAll, var2Ka);
      histos.fill(HIST("Data/hVar1px"), sample, nchAll, var1Pr);
      histos.fill(HIST("Data/hVar2px"), sample, nchAll, var2Pr);
      histos.fill(HIST("Data/hVarpx"), sample, nchPr);
      histos.fill(HIST("Data/hVar2meanptpx"), nchAll, var2Pr);
    }

    if (effSwitchHistoFill) {
      var1Eff = (sumPtWeight * sumPtWeight - sumPtPtWeight) / (sumWeight * (sumWeight - 1));
      var2Eff = (sumPtWeight / sumWeight);

      histos.fill(HIST("hEffVar1x_data"), sample, nchAll, var1Eff);
      histos.fill(HIST("hEffVar2x_data"), sample, nchAll, var2Eff);
      histos.fill(HIST("hEffVarx_data"), sample, nchAll);
      histos.fill(HIST("hEffVar2Meanptx_data"), nchAll, var2Eff);
      histos.fill(HIST("hEffVar1x_Naccorr_data"), sample, sumWeight, var1Eff);
      histos.fill(HIST("hEffVar2x_Naccorr_data"), sample, sumWeight, var2Eff);
      histos.fill(HIST("hEffVarx_Naccorr_data"), sample, sumWeight);
      histos.fill(HIST("hEffVar1x_Naccorr_xaxis_data"), sample, sumWeight, var1);
      histos.fill(HIST("hEffVar2x_Naccorr_xaxis_data"), sample, sumWeight, var2);
      histos.fill(HIST("hcent_nacc_corr"), cent, sumWeight);
      histos.fill(HIST("hNch_vs_corr"), sample, nchAll, sumWeight);
    }
    //================= generated level==============================

    const auto& mccolgen = coll.mcCollision_as<aod::McCollisions>();
    if (std::abs(mccolgen.posZ()) > cVtxZcut) {
      return;
    }
    const auto& mcpartgen = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mccolgen.globalIndex(), cache);
    histos.fill(HIST("hVtxZ_after_gensim"), mccolgen.posZ());

    double nchGenAll = 0., nchGenTrue = 0., nch1 = 0., nchgen = 0.;
    double nchPiGen = 0., nchKaGen = 0., nchPrGen = 0.;
    double q1AllGen = 0, q2AllGen = 0.;
    double q1PiGen = 0, q2PiGen = 0, q1KaGen = 0, q2KaGen = 0, q1PrGen = 0, q2PrGen = 0;
    double var1AllGen = 0, var2AllGen = 0.;
    double var1PiGen = 0, var2PiGen = 0, var1KaGen = 0, var2KaGen = 0, var1PrGen = 0, var2PrGen = 0;

    int sampleGen = histos.get<TH1>(HIST("hVtxZ_after_gensim"))->GetEntries();
    sampleGen = sampleGen % 30;

    for (const auto& mcpart : mcpartgen) {
      // auto  pdgcode = std::abs(mcpart.pdgCode());
      if (!mcpart.isPhysicalPrimary()) {
        continue;
      }
      nch1++;
      histos.fill(HIST("hnch_afterPhysPrimary"), nch1);

      int pid = mcpart.pdgCode();
      auto sign = 0;
      auto* pd = pdg->GetParticle(pid);
      if (pd != nullptr) {
        sign = pd->Charge() / 3.;
      }
      if (sign == 0) {
        continue;
      }
      //    histos.fill(HIST("gen_hSign"), sign);
      if (std::fabs(mcpart.eta()) > cEtacut)
        continue;
      nchGenTrue++;
      histos.fill(HIST("hnch_gen_after_etacut"), nchGenTrue);

      if (mcpart.pt() >= cPtmincut1 && mcpart.pt() <= cPtmaxcut1) {
        nchgen += 1.;
        histos.fill(HIST("hnch_gen"), nchgen);
        histos.fill(HIST("hPtvar_gen"), mcpart.pt());
      }

      if ((mcpart.pt() < cPtmincut) || (mcpart.pt() > cPtmaxcut))
        continue;
      histos.fill(HIST("hPt_gen"), mcpart.pt());
      histos.fill(HIST("hEta_gen"), mcpart.eta());
      histos.fill(HIST("ptHistogram_allcharge_gen"), mcpart.pt());
      nchGenAll += 1.;
      q1AllGen += mcpart.pt();
      q2AllGen += (mcpart.pt() * mcpart.pt());
      histos.fill(HIST("hnch_gen_all"), nchGenAll);
      histos.fill(HIST("hPtEta_gen"), mcpart.pt(), mcpart.eta());

      if (pidSwitch) {
        if (std::fabs(mcpart.y()) < cRapidityCut05) {

          if (mcpart.pdgCode() == PDG_t::kPiPlus || mcpart.pdgCode() == PDG_t::kPiMinus) {
            histos.fill(HIST("ptHistogramPion"), mcpart.pt());
            histos.fill(HIST("hPtEta_pi_gen"), mcpart.pt(), mcpart.eta());
            histos.fill(HIST("hPty_pi_gen"), mcpart.pt(), mcpart.y());
            nchPiGen += 1.;
            q1PiGen += mcpart.pt();
            q2PiGen += (mcpart.pt() * mcpart.pt());
            histos.fill(HIST("hnch_pi"), nchPiGen);
          }

          if (mcpart.pdgCode() == PDG_t::kKPlus || mcpart.pdgCode() == PDG_t::kKMinus) {
            histos.fill(HIST("ptHistogramKaon"), mcpart.pt());
            histos.fill(HIST("hPtEta_ka_gen"), mcpart.pt(), mcpart.eta());
            histos.fill(HIST("hPty_ka_gen"), mcpart.pt(), mcpart.y());
            nchKaGen += 1.;
            q1KaGen += mcpart.pt();
            q2KaGen += (mcpart.pt() * mcpart.pt());
            histos.fill(HIST("hnch_ka"), nchKaGen);
          }

          if (mcpart.pdgCode() == PDG_t::kProton || mcpart.pdgCode() == PDG_t::kProtonBar) {
            histos.fill(HIST("ptHistogramProton"), mcpart.pt());
            histos.fill(HIST("hPtEta_pr_gen"), mcpart.pt(), mcpart.eta());
            histos.fill(HIST("hPty_pr_gen"), mcpart.pt(), mcpart.y());
            nchPrGen += 1.;
            q1PrGen += mcpart.pt();
            q2PrGen += (mcpart.pt() * mcpart.pt());
            histos.fill(HIST("hnch_pr"), nchPrGen);
          }

        } //|y| < 0.5 cut ends!
      } // pid flag
    } // track loop ends!
    histos.fill(HIST("hcent_nacc_gen"), cent, nchGenAll);

    if (nchGenAll < cTwoPtlCut2)
      return;
    var1AllGen = (q1AllGen * q1AllGen - q2AllGen) / (nchGenAll * (nchGenAll - 1));
    var2AllGen = (q1AllGen / nchGenAll);

    histos.fill(HIST("hVar1_gen"), sampleGen, cent, var1AllGen);
    histos.fill(HIST("hVar2_gen"), sampleGen, cent, var2AllGen);
    histos.fill(HIST("hVarc_gen"), sampleGen, cent);

    histos.fill(HIST("hterm1_gen"), nchGenAll, var1AllGen);
    histos.fill(HIST("hterm2_gen"), nchGenAll, var2AllGen);
    //-----------------------nch-------------------------------------
    histos.fill(HIST("hVar1x_gen"), sampleGen, nchGenAll, var1AllGen);
    histos.fill(HIST("hVar2x_gen"), sampleGen, nchGenAll, var2AllGen);
    histos.fill(HIST("hVarx_gen"), sampleGen, nchGenAll);
    histos.fill(HIST("hdiffVar1x_gen"), sampleGen, nchgen, var1AllGen);
    histos.fill(HIST("hdiffVar2x_gen"), sampleGen, nchgen, var2AllGen);
    histos.fill(HIST("hdiffVarx_gen"), sampleGen, nchgen);
    histos.fill(HIST("hVar2meanptx_gen"), nchGenAll, var2AllGen);

    if (pidSwitchHistoFill) {
      //--------------------------Pions-------------------------------------------
      if (nchPiGen >= cTwoPtlCut2) {
        var1PiGen = (q1PiGen * q1PiGen - q2PiGen) / (nchPiGen * (nchPiGen - 1));
        var2PiGen = (q1PiGen / nchPiGen);
      }
      //----------------------- kaons ---------------------------------------
      if (nchKaGen >= cTwoPtlCut2) {
        var1KaGen = (q1KaGen * q1KaGen - q2KaGen) / (nchKaGen * (nchKaGen - 1));
        var2KaGen = (q1KaGen / nchKaGen);
      }
      //---------------------------- protons ----------------------------------
      if (nchPrGen >= cTwoPtlCut2) {
        var1PrGen = (q1PrGen * q1PrGen - q2PrGen) / (nchPrGen * (nchPrGen - 1));
        var2PrGen = (q1PrGen / nchPrGen);
      }
      histos.fill(HIST("hVar1pix_gen"), sampleGen, nchGenAll, var1PiGen);
      histos.fill(HIST("hVar2pix_gen"), sampleGen, nchGenAll, var2PiGen);
      histos.fill(HIST("hVarpix_gen"), sampleGen, nchPiGen);
      histos.fill(HIST("hVar2meanptpix_gen"), nchGenAll, var2PiGen);
      histos.fill(HIST("hVar1kx_gen"), sampleGen, nchGenAll, var1KaGen);
      histos.fill(HIST("hVar2kx_gen"), sampleGen, nchGenAll, var2KaGen);
      histos.fill(HIST("hVarkx_gen"), sampleGen, nchKaGen);
      histos.fill(HIST("hVar2meanptkx_gen"), nchGenAll, var2KaGen);
      histos.fill(HIST("hVar1px_gen"), sampleGen, nchGenAll, var1PrGen);
      histos.fill(HIST("hVar2px_gen"), sampleGen, nchGenAll, var2PrGen);
      histos.fill(HIST("hVarpx_gen"), sampleGen, nchPrGen);
      histos.fill(HIST("hVar2meanptpx_gen"), nchGenAll, var2PrGen);
    }

  } // void process
  PROCESS_SWITCH(EventMeanPtId, processMcReco, "Process reconstructed", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<EventMeanPtId>(cfgc)};
  return workflow;
}
