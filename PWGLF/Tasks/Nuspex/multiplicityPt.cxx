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

/// \author Dushmanta Sahu (dushmanta.sahu@cern.ch)
/// \file multiplicityPt.cxx
/// \brief Analysis to do PID with MC - Full correction factors for pions, kaons, protons

#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/DataModel/spectraTOF.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>

#include <TF1.h>
#include <TMCProcess.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels,
                          aod::Run3MatchedToBCSparse>;

// Reconstructed collisions joined with MC label + centrality
using ColEvSelsMC = soa::Join<aod::Collisions, aod::EvSels,
                              aod::McCollisionLabels,
                              aod::CentFT0Ms,
                              aod::TPCMults, aod::PVMults>;

// Data collision table (for processData)
using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels,
                                     aod::CentFT0Ms,
                                     aod::TPCMults, aod::PVMults>;

// Track tables - with TPC PID only
using TrackTableData = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                 aod::TrackSelection,
                                 aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;

using TracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                           aod::TrackSelection, aod::McTrackLabels,
                           aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;

static constexpr int NCentHists{10};
static constexpr int NPartHists{5};
std::array<std::shared_ptr<TH3>, NCentHists> hDedxVsMomentumVsCentPos{};
std::array<std::shared_ptr<TH3>, NCentHists> hDedxVsMomentumVsCentNeg{};
std::array<std::shared_ptr<TH3>, NCentHists + 1> hDedxVspTMomentumVsCent{};
std::array<std::shared_ptr<TH2>, NCentHists + 1> hMomentumVsEtaPos{};
std::array<std::shared_ptr<TH2>, NCentHists + 1> hMomentumVsEtaNeg{};
std::array<std::shared_ptr<TH2>, NCentHists + 1> hpTVsEtaPos{};
std::array<std::shared_ptr<TH2>, NCentHists + 1> hpTVsEtaNeg{};
// Total counts
std::array<std::shared_ptr<TH2>, NCentHists + 1> hTotalMomPosCent{};
std::array<std::shared_ptr<TH2>, NCentHists + 1> hTotalMomNegCent{};
std::array<std::shared_ptr<TH2>, NCentHists + 1> hTotalPtPosCent{};
std::array<std::shared_ptr<TH2>, NCentHists + 1> hTotalPtNegCent{};
// Counts for particles
std::array<std::array<std::shared_ptr<TH2>, NCentHists + 1>, NPartHists> hFracMomPosCent{};
std::array<std::array<std::shared_ptr<TH2>, NCentHists + 1>, NPartHists> hFracMomNegCent{};
std::array<std::array<std::shared_ptr<TH2>, NCentHists + 1>, NPartHists> hFracPtPosCent{};
std::array<std::array<std::shared_ptr<TH2>, NCentHists + 1>, NPartHists> hFracPtNegCent{};

struct MultiplicityPt {

  // ── Services ──────────────────────────────────────────────
  Service<o2::framework::O2DatabasePDG> pdg;
  Service<ccdb::BasicCCDBManager> ccdb;

  // ── Constant values ──────────────────────────────────
  static constexpr int CentBinMax = 100;
  static constexpr int NEvLabel = 15;
  static constexpr int ResponseMatrixTypes = 7;

  // ── Configurables: event ──────────────────────────────────
  Configurable<bool> isRun3{"isRun3", true, "is Run3 dataset"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<int> cfgINELCut{"cfgINELCut", 0, "INEL event selection: 0 no sel, 1 INEL>0, 2 INEL>1"};
  Configurable<bool> askForCustomTVX{"askForCustomTVX", false, "Ask for custom TVX rather than sel8"};
  Configurable<bool> removeITSROFrameBorder{"removeITSROFrameBorder", false, "Remove ITS Read-Out Frame border"};
  Configurable<bool> removeNoSameBunchPileup{"removeNoSameBunchPileup", false, "Remove no same bunch pileup"};
  Configurable<bool> requireIsGoodZvtxFT0vsPV{"requireIsGoodZvtxFT0vsPV", false, "Require good Z vertex FT0 vs PV"};
  Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "Require vertex ITSTPC"};
  Configurable<bool> removeNoTimeFrameBorder{"removeNoTimeFrameBorder", false, "Remove no time frame border"};

  // Gen-level event selection
  Configurable<bool> selTVXMC{"selTVXMC", true, "Require TVX-equivalent at gen level"};
  Configurable<bool> isZvtxPosSelMC{"isZvtxPosSelMC", true, "Require |Zvtx|<cut at gen level"};

  // ── Configurables: track ──────────────────────────────────
  Configurable<float> cfgCutEtaMax{"cfgCutEtaMax", 0.8f, "Max eta range for tracks"};
  Configurable<float> cfgCutEtaMin{"cfgCutEtaMin", -0.8f, "Min eta range for tracks"};
  Configurable<float> cfgCutY{"cfgCutY", 0.5f, "Y range for tracks"};
  Configurable<float> cfgCutNsigma{"cfgCutNsigma", 3.0f, "nsigma cut range for tracks"};
  Configurable<int> lastRequiredTrdCluster{"lastRequiredTrdCluster", -1, "Last cluster to require in TRD"};
  Configurable<bool> requireTrdOnly{"requireTrdOnly", false, "Require only tracks from TRD"};
  Configurable<bool> requireNoTrd{"requireNoTrd", false, "Require tracks without TRD"};
  Configurable<int> multiplicityEstimator{"multiplicityEstimator", 6,
                                          "Multiplicity estimator: 0=NoMult, 1=MultFV0M, 2=MultFT0M, 3=MultFDDM, 4=MultTracklets, 5=MultTPC, 6=MultNTracksPV, 7=MultNTracksPVeta1, 8=CentFT0C, 9=CentFT0M, 10=CentFV0A"};

  // Track cuts
  Configurable<bool> enableDCAHistograms{"enableDCAHistograms", false, "Enable DCA histograms"};
  Configurable<bool> enablePIDHistograms{"enablePIDHistograms", true, "Enable PID histograms"};
  Configurable<bool> useCustomTrackCuts{"useCustomTrackCuts", true, "Flag to use custom track cuts"};
  Configurable<int> itsPattern{"itsPattern", 0, "0 = Run3ITSibAny, 1 = Run3ITSallAny, 2 = Run3ITSall7Layers, 3 = Run3ITSibTwo"};
  Configurable<bool> requireITS{"requireITS", true, "Additional cut on the ITS requirement"};
  Configurable<bool> requireTPC{"requireTPC", true, "Additional cut on the TPC requirement"};
  Configurable<bool> requireGoldenChi2{"requireGoldenChi2", true, "Additional cut on the GoldenChi2"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.f, "Additional cut on the minimum number of crossed rows in the TPC"};
  Configurable<float> minNCrossedRowsOverFindableClustersTPC{"minNCrossedRowsOverFindableClustersTPC", 0.8f, "Additional cut on the minimum value of the ratio between crossed rows and findable clusters in the TPC"};
  Configurable<float> maxChi2PerClusterTPC{"maxChi2PerClusterTPC", 4.f, "Additional cut on the maximum value of the chi2 per cluster in the TPC"};
  Configurable<float> minChi2PerClusterTPC{"minChi2PerClusterTPC", 0.5f, "Additional cut on the minimum value of the chi2 per cluster in the TPC"};
  Configurable<float> maxChi2PerClusterITS{"maxChi2PerClusterITS", 36.f, "Additional cut on the maximum value of the chi2 per cluster in the ITS"};
  Configurable<float> nSigmaDCAxy{"nSigmaDCAxy", 1.f, "Additional cut on the maximum value of the DCA xy (multiplicative factor)"};
  Configurable<float> dcaXYp0{"dcaXYp0", 0.0105f, "DCAxy formula: p0 + p1/pt^p2"};
  Configurable<float> dcaXYp1{"dcaXYp1", 0.0350f, "DCAxy p1 parameter"};
  Configurable<float> dcaXYp2{"dcaXYp2", 1.1f, "DCAxy p2 parameter"};
  Configurable<float> nSigmaDCAz{"nSigmaDCAz", 1.f, "Additional cut on the maximum value of the DCA z (multiplicative factor)"};
  Configurable<float> dcaZp0{"dcaZp0", 0.0105f, "DCAz formula: p0 + p1/pt^p2"};
  Configurable<float> dcaZp1{"dcaZp1", 0.0350f, "DCAz p1 parameter"};
  Configurable<float> dcaZp2{"dcaZp2", 1.1f, "DCAz p2 parameter"};
  // Configurable<float> maxDcaZ{"maxDcaZ", 2.0f, "Additional cut on the maximum value of the DCA z"};
  Configurable<float> minTPCNClsFound{"minTPCNClsFound", 70.0f, "min number of found TPC clusters"};
  Configurable<float> minTPCNClsPID{"minTPCNClsPID", 130.0f, "min number of PID TPC clusters"};
  Configurable<bool> nClTPCFoundCut{"nClTPCFoundCut", false, "Apply TPC found clusters cut"};
  Configurable<bool> nClTPCPIDCut{"nClTPCPIDCut", true, "Apply TPC clusters for PID cut"};
  Configurable<int> minITSnClusters{"minITSnClusters", 5, "minimum number of found ITS clusters"};
  Configurable<float> cfgTrkLowPtCut{"cfgTrkLowPtCut", 0.15f, "Minimum constituent pT"};

  // PID selection
  Configurable<float> cfgCutNsigmaPi{"cfgCutNsigmaPi", 3.0f, "nsigma cut for pions"};
  Configurable<float> cfgCutNsigmaKa{"cfgCutNsigmaKa", 2.5f, "nsigma cut for kaons"};
  Configurable<float> cfgCutNsigmaPr{"cfgCutNsigmaPr", 2.5f, "nsigma cut for protons"};

  // Phi cut parameters
  Configurable<bool> applyPhiCut{"applyPhiCut", true, "Apply phi sector cut"};
  Configurable<float> pTthresholdPhiCut{"pTthresholdPhiCut", 2.0f, "pT threshold above which to apply phi cut"};
  Configurable<double> phiCutLowParam1{"phiCutLowParam1", 0.119297, "First parameter for low phi cut"};
  Configurable<double> phiCutLowParam2{"phiCutLowParam2", 0.000379693, "Second parameter for low phi cut"};
  Configurable<double> phiCutHighParam1{"phiCutHighParam1", 0.16685, "First parameter for high phi cut"};
  Configurable<double> phiCutHighParam2{"phiCutHighParam2", 0.00981942, "Second parameter for high phi cut"};

  // ── Nch window for multiplicity axis ──────────────────────
  Configurable<float> tpcNchAcceptance{"tpcNchAcceptance", 0.8f,
                                       "|eta| window for counting gen. Nch (multiplicity axis)"};

  // ── Histogram binning ────────────────────────────────────
  Configurable<int> nBinsNch{"nBinsNch", 200, "Bins on gen-Nch axis"};
  Configurable<float> maxNch{"maxNch", 200.0f, "Max gen Nch on histogram axis"};
  Configurable<int> nBinsNPV{"nBinsNPV", 600, "N bins ITS tracks"};
  Configurable<float> minNpv{"minNpv", 0, "Min NPV"};
  Configurable<float> maxNpv{"maxNpv", 600, "Max NPV"};
  ConfigurableAxis ptBinning{"ptBinning", {VARIABLE_WIDTH, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0}, "pT bin limits"};
  ConfigurableAxis pFineBins{"pFineBins", {1995, 0.1, 40}, "Binning for momentum"};
  ConfigurableAxis dedxBins{"dedxBins", {100, 0, 100}, "Binning for dedx"};
  std::vector<double> centBinningStd = {0., 1., 5., 10., 15., 20., 30., 40., 50., 70., 100.};
  ConfigurableAxis dcaBins{"dcaBins", {200, -0.5, 0.5}, "Binning for DCA plots"};
  // ── Custom track-selection object ────────────────────────
  TrackSelection customTrackCuts;

  // ── TF1 pointers for phi cuts ────────────────────────────
  TF1* fphiCutLow = nullptr;
  TF1* fphiCutHigh = nullptr;

  // ── Histogram registry ───────────────────────────────────
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryFrac{"registryFrac", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // ── Preslice: group tracks by reco collision ──────────────
  Preslice<TracksMC> perCollision = aod::track::collisionId;

  // FT0 acceptance for TVX-equivalent gen-level selection
  static constexpr float MinFT0A = 3.5f;
  static constexpr float MaxFT0A = 4.9f;
  static constexpr float MinFT0C = -3.3f;
  static constexpr float MaxFT0C = -2.1f;

  static constexpr float MinCharge = 3.0f;
  static constexpr int CentralityClasses = 10;

  static constexpr int ParticleTypes = 4;
  // Response Matrix histogram names
  static constexpr std::string_view EtavspvspTPosPart[ResponseMatrixTypes] = {"heta_vs_pt_vs_p_all_Pos", "heta_vs_pt_vs_p_all_Pos_Pri", "heta_vs_pt_vs_p_all_Pos_Pri_MC", "heta_vs_pt_vs_p_all_Pos_Pri_MC_Part", "heta_vs_pt_vs_p_Pi_Pos", "heta_vs_pt_vs_p_K_Pos", "heta_vs_pt_vs_p_Pr_Pos"};
  static constexpr std::string_view EtavspvspTNegPart[ResponseMatrixTypes] = {"heta_vs_pt_vs_p_all_Neg", "heta_vs_pt_vs_p_all_Neg_Pri", "heta_vs_pt_vs_p_all_Neg_Pri_MC", "heta_vs_pt_vs_p_all_Neg_Pri_MC_Part", "heta_vs_pt_vs_p_Pi_Neg", "heta_vs_pt_vs_p_K_Neg", "heta_vs_pt_vs_p_Pr_Neg"};

  // Event counter bins
  enum EvCutLabel {
    kAllGen = 1,
    kTVXequiv,
    kVtxZ,
    kINELgt0,
    kRecoColl,
    kRecoSelected
  };

  // Particle species enum
  enum ParticleSpecies : int {
    kPion = 0,
    kKaon = 1,
    kProton = 2,
    kAllCharged = 3,
    kNSpecies = 4
  };

  enum INELCutSelection : int {
    INEL = 0,
    INELgt0 = 1,
    INELgt1 = 2
  };

  void init(InitContext const&)
  {
    // Setup custom track cuts
    if (useCustomTrackCuts.value) {
      customTrackCuts = getGlobalTrackSelectionRun3ITSMatch(itsPattern.value);
      customTrackCuts.SetRequireITSRefit(requireITS.value);
      customTrackCuts.SetRequireTPCRefit(requireTPC.value);
      customTrackCuts.SetMinNClustersITS(minITSnClusters.value);
      customTrackCuts.SetRequireGoldenChi2(requireGoldenChi2.value);
      customTrackCuts.SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC.value);
      customTrackCuts.SetMaxChi2PerClusterITS(maxChi2PerClusterITS.value);
      customTrackCuts.SetMinNCrossedRowsTPC(minNCrossedRowsTPC.value);
      customTrackCuts.SetMinNCrossedRowsOverFindableClustersTPC(minNCrossedRowsOverFindableClustersTPC.value);
      customTrackCuts.SetMaxDcaXYPtDep([](float /*pt*/) { return 10000.f; });
      // customTrackCuts.SetMaxDcaZ(maxDcaZ.value);
    }

    // Initialize phi cut functions if enabled
    if (applyPhiCut.value) {
      fphiCutLow = new TF1("StandardPhiCutLow",
                           Form("%f/x/x+pi/18.0-%f",
                                phiCutLowParam1.value, phiCutLowParam2.value),
                           0, 50);
      fphiCutHigh = new TF1("StandardPhiCutHigh",
                            Form("%f/x+pi/18.0+%f",
                                 phiCutHighParam1.value, phiCutHighParam2.value),
                            0, 50);
      LOGF(info, "Phi cut ENABLED for pT > %.1f GeV/c", pTthresholdPhiCut.value);
    }

    // Define axes
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec dedxAxis = {dedxBins, "dE/dx (a. u.)"};
    AxisSpec etaAxis{8, -0.8, 0.8, "#eta"};
    AxisSpec pAxis = {ptBinning, "#it{p} (GeV/#it{c})"};
    AxisSpec pFineAxis{pFineBins, "#it{p} (GeV/c)"};
    AxisSpec pTFineAxis{pFineBins, "#it{p}_{T} (GeV/c)"};

    const AxisSpec centAxis{centBinningStd, "FT0M Centrality (%)"};

    // Fine centrality binning (100 bins)
    std::vector<double> centBinningFine;
    for (int i = 0; i <= CentBinMax; i++) {
      centBinningFine.push_back(static_cast<double>(i));
    }
    const AxisSpec centFineAxis{centBinningFine, "FT0M Centrality (%)"};

    const AxisSpec nchAxis{nBinsNch.value, 0.f, maxNch.value,
                           Form("Gen. N_{ch} (|#eta|<%.1f)", tpcNchAcceptance.value)};

    const AxisSpec npvAxis{nBinsNPV.value, minNpv.value, maxNpv.value, "NPV"};
    const AxisSpec dcaXYAxis{105, -1.05f, 1.05f, "DCA_{xy} (cm)"};
    const AxisSpec zvtxAxis{60, -30.0, 30.0, "Vtx_{z} (cm)"};
    const AxisSpec nclAxis{161, -0.5, 160.5, "N_{cl} TPC"};
    AxisSpec dcaAxis{dcaBins, ""};
    // ========================================================================
    // EVENT COUNTER AND BASIC HISTOGRAMS
    // ========================================================================
    registry.add("EventCounter", ";;Events", kTH1F, {{8, 0.5, 8.5}});
    {
      auto h = registry.get<TH1>(HIST("EventCounter"));
      h->GetXaxis()->SetBinLabel(kAllGen, "All gen.");
      h->GetXaxis()->SetBinLabel(kTVXequiv, "TVX-equiv.");
      h->GetXaxis()->SetBinLabel(kVtxZ, "|Zvtx|<cut");
      h->GetXaxis()->SetBinLabel(kINELgt0, "INEL>0");
      h->GetXaxis()->SetBinLabel(kRecoColl, ">=1 reco coll.");
      h->GetXaxis()->SetBinLabel(kRecoSelected, ">=1 reco+sel.");
    }

    registry.add("NumberOfRecoCollisions", "Reco collisions per gen. collision;N_{reco};Entries",
                 kTH1F, {{10, -0.5, 9.5}});
    registry.add("zPosMC", "Gen. Vtx_{z} (evts w/ >=1 reco+sel. coll.);Vtx_{z} (cm);Entries",
                 kTH1F, {zvtxAxis});
    registry.add("zPosReco", "Reco. Vtx_{z} (selected);Vtx_{z} (cm);Entries",
                 kTH1F, {zvtxAxis});
    registry.add("T0Ccent", "FT0M centrality (selected);Centrality (%);Entries",
                 kTH1F, {centAxis});
    registry.add("T0CcentVsFoundFT0", "Found(1.5) NOT Found(0.5);;Status;", kTH2F, {{centAxis, {2, 0, 2}}});
    registry.add("T0CcentVsFoundFT0AndTVX", "Found(1.5) NOT Found(0.5);;Status;", kTH2F, {{centAxis, {2, 0, 2}}});
    registry.add("hvtxZ", "Vertex Z (data);Vertex Z (cm);Events", kTH1F, {{40, -20.0, 20.0}});
    registry.add("hvtxZmc", "MC vertex Z;Vertex Z (cm);Events", kTH1F, {{40, -20.0, 20.0}});

    // ========================================================================
    // CENTRALITY DIAGNOSTIC HISTOGRAMS
    // ========================================================================
    registry.add("Centrality/hCentRaw", "Raw FT0M Centrality (no cuts);Centrality (%);Counts",
                 kTH1D, {centFineAxis});
    registry.add("Centrality/hCentAfterVtx", "Centrality after vertex cut;Centrality (%);Counts",
                 kTH1D, {centFineAxis});
    registry.add("Centrality/hCentAfterINEL", "Centrality after INEL cut;Centrality (%);Counts",
                 kTH1D, {centFineAxis});
    registry.add("Centrality/hCentAfterAll", "Centrality after all cuts;Centrality (%);Counts",
                 kTH1D, {centFineAxis});
    registry.add("Centrality/hCentVsMult", "Centrality vs Generated Multiplicity;Centrality (%);N_{ch}^{gen}",
                 kTH2D, {centFineAxis, nchAxis});
    registry.add("Centrality/hMultVsCent", "Generated Multiplicity vs Centrality;N_{ch}^{gen};Centrality (%)",
                 kTH2D, {nchAxis, centFineAxis});
    registry.add("Centrality/hCentVsVz", "Centrality vs Vertex Z;Centrality (%);V_{z} (cm)",
                 kTH2D, {centFineAxis, {40, -20, 20}});
    registry.add("Centrality/hRecoMultVsCent", "Reconstructed Track Multiplicity vs Centrality;Centrality (%);N_{tracks}^{reco}",
                 kTH2D, {centFineAxis, {100, 0, 100}});
    registry.add("Centrality/hVertexResVsCent", "Vertex Resolution vs Centrality;Centrality (%);V_{z} resolution (cm)",
                 kTH2D, {centFineAxis, {100, -1, 1}});
    registry.add("NchVsNPV", ";Nch; NPV;", kTH2D, {{npvAxis, nchAxis}});
    registry.add("ExcludedEvtVsNch", ";Nch;Entries;", kTH1F, {nchAxis});
    registry.add("ExcludedEvtVsNPV", ";NPV;Entries;", kTH1F, {npvAxis});

    // ========================================================================
    // INEL CLASS HISTOGRAMS
    // ========================================================================
    registry.add("INEL/hINELClass", "INEL Class for MC Collisions;INEL Class;Counts",
                 kTH1D, {{3, 0.5, 3.5}});
    auto hINEL = registry.get<TH1>(HIST("INEL/hINELClass"));
    hINEL->GetXaxis()->SetBinLabel(1, "INEL0");
    hINEL->GetXaxis()->SetBinLabel(2, "INEL>0");
    hINEL->GetXaxis()->SetBinLabel(3, "INEL>1");
    registry.add("INEL/hINELVsCent", "INEL Class vs Centrality;Centrality (%);INEL Class",
                 kTH2D, {centFineAxis, {3, 0.5, 3.5}});

    // ========================================================================
    // CUT FLOW HISTOGRAMS
    // ========================================================================
    registry.add("CutFlow/hCutStats", "Cut Statistics;Cut Stage;Counts",
                 kTH1D, {{6, 0.5, 6.5}});
    auto hCut = registry.get<TH1>(HIST("CutFlow/hCutStats"));
    hCut->GetXaxis()->SetBinLabel(1, "All reco events");
    hCut->GetXaxis()->SetBinLabel(2, "Has MC match");
    hCut->GetXaxis()->SetBinLabel(3, "Has centrality");
    hCut->GetXaxis()->SetBinLabel(4, "Pass vertex");
    hCut->GetXaxis()->SetBinLabel(5, "Pass INEL");
    hCut->GetXaxis()->SetBinLabel(6, "Selected");
    registry.add("CutFlow/hCentPerCut", "Centrality Distribution at Each Cut;Cut Stage;Centrality (%)",
                 kTH2D, {{6, 0.5, 6.5}, centFineAxis});
    registry.add("evsel", "Event selection", HistType::kTH1D, {{20, 0.5, 20.5}});
    auto hEvSel = registry.get<TH1>(HIST("evsel"));
    for (int i = 1; i <= NEvLabel; i++) {
      hEvSel->GetXaxis()->SetBinLabel(i, Form("Step %d", i));
    }

    // ========================================================================
    // MC COLLISION HISTOGRAMS
    // ========================================================================
    registry.add("MC/GenRecoCollisions", "Generated and Reconstructed MC Collisions",
                 kTH1D, {{10, 0.5, 10.5}});
    auto hColl = registry.get<TH1>(HIST("MC/GenRecoCollisions"));
    hColl->GetXaxis()->SetBinLabel(1, "Collisions generated");
    hColl->GetXaxis()->SetBinLabel(2, "Collisions reconstructed");
    hColl->GetXaxis()->SetBinLabel(3, "INEL>0");
    hColl->GetXaxis()->SetBinLabel(4, "INEL>1");
    registry.add("NchMCcentVsTVX", ";Passed(1.5) NOT Passed(0.5);", kTH2F, {{nchAxis, {2, 0, 2}}});
    registry.add("Centrality_AllRecoEvt", "Generated Events Irrespective of number of reconstructions;Centrality;Entries", kTH1F, {centAxis});
    registry.add("Centrality_WRecoEvt", "Generated Events With at least One Rec. Collision;Centrality;Entries", kTH1F, {centAxis});
    registry.add("Centrality_WRecoEvtWSelCri", "Generated Events With at least One Rec. Collision + Sel. criteria;Centrality;Entries", kTH1F, {centAxis});

    // ========================================================================
    // 2D and 3D CORRELATION HISTOGRAMS
    // ========================================================================
    registry.add("NchMCVsCent", "Gen. N_{ch} vs FT0M centrality (reco+sel. evts);FT0M Centrality (%);Gen. N_{ch}",
                 kTH2F, {{centAxis, nchAxis}});
    registry.add("NchMCVsCentVsPt", "pT vs Gen. N_{ch} vs FT0M centrality;FT0M Centrality (%);Gen. N_{ch};#it{p}_{T}",
                 kTH3F, {{centAxis, nchAxis, ptAxis}});
    registry.add("hEta", "Track eta;#eta;Counts", kTH1D, {{20, -0.8, 0.8}});
    registry.add("hPhi", "Track phi;#varphi (rad);Counts", kTH1D, {{64, 0, o2::constants::math::TwoPI}});
    registry.add("EtaVsPhi", ";#eta;#varphi;", kTH2F, {{50, -1.0, 1.0}, {100, 0, o2::constants::math::TwoPI}});

    // ========================================================================
    // EVENT LOSS HISTOGRAMS
    // ========================================================================
    registry.add("NchMC_AllGen", "EVENT LOSS denom.;Gen. N_{ch};Entries", kTH1F, {nchAxis});
    registry.add("NchMC_WithRecoEvt", "EVENT LOSS numer.;Gen. N_{ch};Entries", kTH1F, {nchAxis});
    registry.add("MC/EventLoss/NchGenerated", "Generated charged multiplicity;N_{ch}^{gen};Counts", kTH1D, {nchAxis});
    registry.add("MC/EventLoss/NchGenerated_PhysicsSelected", "Generated charged multiplicity (physics selected);N_{ch}^{gen};Counts", kTH1D, {nchAxis});
    registry.add("MC/EventLoss/NchGenerated_Reconstructed", "Generated charged multiplicity (reconstructed);N_{ch}^{gen};Counts", kTH1D, {nchAxis});
    registry.add("MC/EventLoss/GenMultVsCent", "Generated charged particles vs FT0M centrality;FT0M Centrality (%);N_{ch}^{gen}", kTH2D, {centAxis, nchAxis});
    registry.add("MC/EventLoss/GenMultVsCent_Selected", "Generated vs FT0M centrality (selected);FT0M Centrality (%);N_{ch}^{gen}", kTH2D, {centAxis, nchAxis});
    registry.add("MC/EventLoss/GenMultVsCent_Rejected", "Generated vs FT0M centrality (rejected);FT0M Centrality (%);N_{ch}^{gen}", kTH2D, {centAxis, nchAxis});
    registry.add("hEventLossBreakdown", "Event loss breakdown", kTH1D, {{4, 0.5, 4.5}});
    auto hLoss = registry.get<TH1>(HIST("hEventLossBreakdown"));
    hLoss->GetXaxis()->SetBinLabel(1, "Physics selected");
    hLoss->GetXaxis()->SetBinLabel(2, "Reconstructed");
    hLoss->GetXaxis()->SetBinLabel(3, "Selected");
    hLoss->GetXaxis()->SetBinLabel(4, "Final efficiency");

    // ========================================================================
    // SIGNAL LOSS HISTOGRAMS
    // ========================================================================
    const std::vector<std::string> speciesNames = {"Pi", "Ka", "Pr", "All"};
    for (int i = 0; i < ParticleTypes; ++i) {
      const std::string& name = speciesNames[i];
      registry.add(Form("Pt%sVsNchMC_AllGen", name.c_str()),
                   Form("SIGNAL LOSS denom. (%s): all gen. evts.;#it{p}_{T};Gen. N_{ch}", name.c_str()),
                   kTH2F, {{ptAxis, nchAxis}});
      registry.add(Form("Pt%sVsNchMC_WithRecoEvt", name.c_str()),
                   Form("SIGNAL LOSS numer. (%s): gen. evts. w/ >=1 reco+sel.;#it{p}_{T};Gen. N_{ch}", name.c_str()),
                   kTH2F, {{ptAxis, nchAxis}});
    }

    // ========================================================================
    // TRACKING EFFICIENCY HISTOGRAMS
    // ========================================================================
    for (int i = 0; i < ParticleTypes; ++i) {
      const std::string& name = speciesNames[i];
      registry.add(Form("Pt%sVsCentMC_WithRecoEvt", name.c_str()),
                   Form("EFF denom. (%s): gen pT in reco+sel. evts.;#it{p}_{T}^{gen};FT0M Cent. (%s)", name.c_str(), ")"),
                   kTH2F, {{ptAxis, centAxis}});
      registry.add(Form("Pt%sVsCent_WithRecoEvt", name.c_str()),
                   Form("EFF numer. (%s): reco primaries (gen pT);#it{p}_{T}^{gen};FT0M Cent. (%s)", name.c_str(), ")"),
                   kTH2F, {{ptAxis, centAxis}});
      registry.add(Form("PtGen%sVsNchMC_WithRecoEvt", name.c_str()),
                   Form("EFF denom. (%s) vs gen Nch;#it{p}_{T}^{gen};Gen. N_{ch}", name.c_str()),
                   kTH2F, {{ptAxis, nchAxis}});
      registry.add(Form("PtGen%sVsNchMC_RecoTrk", name.c_str()),
                   Form("EFF numer. (%s) vs gen Nch;#it{p}_{T}^{gen};Gen. N_{ch}", name.c_str()),
                   kTH2F, {{ptAxis, nchAxis}});
    }

    // ========================================================================
    // MC CLOSURE HISTOGRAMS
    // ========================================================================
    registry.add("MCclosure_PtMCPiVsNchMC", "MC closure Pi (gen);#it{p}_{T}^{gen};Gen. N_{ch}", kTH2F, {{ptAxis, nchAxis}});
    registry.add("MCclosure_PtMCKaVsNchMC", "MC closure Ka (gen);#it{p}_{T}^{gen};Gen. N_{ch}", kTH2F, {{ptAxis, nchAxis}});
    registry.add("MCclosure_PtMCPrVsNchMC", "MC closure Pr (gen);#it{p}_{T}^{gen};Gen. N_{ch}", kTH2F, {{ptAxis, nchAxis}});
    registry.add("MCclosure_PtPiVsNchMC", "MC closure Pi (reco);#it{p}_{T}^{reco};Gen. N_{ch}", kTH2F, {{ptAxis, nchAxis}});
    registry.add("MCclosure_PtKaVsNchMC", "MC closure Ka (reco);#it{p}_{T}^{reco};Gen. N_{ch}", kTH2F, {{ptAxis, nchAxis}});
    registry.add("MCclosure_PtPrVsNchMC", "MC closure Pr (reco);#it{p}_{T}^{reco};Gen. N_{ch}", kTH2F, {{ptAxis, nchAxis}});
    registry.add("MC/GenPtVsNch", "Generated pT vs Multiplicity;#it{p}_{T};N_{ch}^{gen}", kTH2D, {ptAxis, nchAxis});
    registry.add("MC/GenPtVsNch_PhysicsSelected", "Generated pT vs Multiplicity (physics selected);#it{p}_{T};N_{ch}^{gen}", kTH2D, {ptAxis, nchAxis});

    // ========================================================================
    // DCA HISTOGRAMS FOR PRIMARY FRACTION
    // ========================================================================
    registry.add("dcaVsPtPi", "PRIMARY FRAC. (Pi) - primaries;#it{p}_{T};DCA_{xy};FT0M Cent.",
                 kTH3F, {{ptAxis, dcaXYAxis, centAxis}});
    registry.add("dcaVsPtPr", "PRIMARY FRAC. (Pr) - primaries;#it{p}_{T};DCA_{xy};FT0M Cent.",
                 kTH3F, {{ptAxis, dcaXYAxis, centAxis}});
    registry.add("dcaVsPtPiDec", "PRIMARY FRAC. (Pi) - sec. from decays;#it{p}_{T};DCA_{xy};FT0M Cent.",
                 kTH3F, {{ptAxis, dcaXYAxis, centAxis}});
    registry.add("dcaVsPtPrDec", "PRIMARY FRAC. (Pr) - sec. from decays;#it{p}_{T};DCA_{xy};FT0M Cent.",
                 kTH3F, {{ptAxis, dcaXYAxis, centAxis}});
    registry.add("dcaVsPtPiMat", "PRIMARY FRAC. (Pi) - sec. from material;#it{p}_{T};DCA_{xy};FT0M Cent.",
                 kTH3F, {{ptAxis, dcaXYAxis, centAxis}});
    registry.add("dcaVsPtPrMat", "PRIMARY FRAC. (Pr) - sec. from material;#it{p}_{T};DCA_{xy};FT0M Cent.",
                 kTH3F, {{ptAxis, dcaXYAxis, centAxis}});

    // ========================================================================
    // MEASURED SPECTRA
    // ========================================================================
    for (int i = 0; i < ParticleTypes; ++i) {
      const std::string& name = speciesNames[i];
      registry.add(Form("Pt%sMeasuredVsCent", name.c_str()),
                   Form("Measured %s (PID);#it{p}_{T};FT0M Centrality (%s)", name.c_str(), ")"),
                   kTH2F, {{ptAxis, centAxis}});
      registry.add(Form("Pt%sMeasuredVsNch", name.c_str()),
                   Form("Measured %s (PID) vs gen Nch;#it{p}_{T};Gen. N_{ch}", name.c_str()),
                   kTH2F, {{ptAxis, nchAxis}});
    }

    // ========================================================================
    // TPC CLUSTER HISTOGRAMS
    // ========================================================================
    registry.add("hNclFoundTPC", "Number of TPC found clusters;N_{cl, found};Counts", kTH1D, {nclAxis});
    registry.add("hNclPIDTPC", "Number of TPC PID clusters;N_{cl, PID};Counts", kTH1D, {nclAxis});
    registry.add("hNclFoundTPCvsPt", "TPC found clusters vs pT;#it{p}_{T};N_{cl,found}", kTH2D, {ptAxis, nclAxis});
    registry.add("hNclPIDTPCvsPt", "TPC PID clusters vs pT;#it{p}_{T};N_{cl,PID}", kTH2D, {ptAxis, nclAxis});

    // ========================================================================
    // INCLUSIVE HISTOGRAMS
    // ========================================================================
    registry.add("Inclusive/hPtPrimGenAll", "All generated primaries (no cuts);#it{p}_{T};Counts", kTH1D, {ptAxis});
    registry.add("Inclusive/hPtPrimBadVertex", "Generated primaries (bad vertex);#it{p}_{T};Counts", kTH1D, {ptAxis});
    registry.add("Inclusive/hPtPrimGen", "Generated primaries (after physics selection);#it{p}_{T};Counts", kTH1D, {ptAxis});
    registry.add("Inclusive/hPtPrimRecoEv", "Generated primaries (reco events);#it{p}_{T};Counts", kTH1D, {ptAxis});
    registry.add("Inclusive/hPtPrimGoodEv", "Generated primaries (good events);#it{p}_{T};Counts", kTH1D, {ptAxis});
    registry.add("Inclusive/hPtNumEff", "Tracking efficiency numerator;#it{p}_{T};Counts", kTH1D, {ptAxis});
    registry.add("Inclusive/hPtDenEff", "Tracking efficiency denominator;#it{p}_{T};Counts", kTH1D, {ptAxis});
    registry.add("Inclusive/hPtAllReco", "All reconstructed tracks;#it{p}_{T};Counts", kTH1D, {ptAxis});
    registry.add("Inclusive/hPtPrimReco", "Reconstructed primaries;#it{p}_{T};Counts", kTH1D, {ptAxis});
    registry.add("Inclusive/hPtSecReco", "Reconstructed secondaries;#it{p}_{T};Counts", kTH1D, {ptAxis});
    registry.add("Inclusive/hPtMeasured", "All measured tracks;#it{p}_{T};Counts", kTH1D, {ptAxis});
    registry.add("Inclusive/hPtMeasuredVsCent", "All measured tracks (PID) vs centrality;#it{p}_{T};FT0M Centrality (%s)", kTH2D, {ptAxis, centAxis});
    registry.add("Inclusive/hPtMeasuredVsMult", "All measured tracks vs mult;#it{p}_{T};Mult Class (%)", kTH2D, {ptAxis, nchAxis});

    // Inclusive vs Multiplicity
    registry.add("Inclusive/hPtPrimGenAllVsMult", "All generated primaries vs mult;#it{p}_{T};Mult Class (%)", kTH2D, {ptAxis, nchAxis});
    registry.add("Inclusive/hPtPrimBadVertexVsMult", "Generated primaries (bad vertex) vs mult;#it{p}_{T};Mult Class (%)", kTH2D, {ptAxis, nchAxis});
    registry.add("Inclusive/hPtPrimGenVsMult", "Generated primaries (after phys sel) vs mult;#it{p}_{T};Mult Class (%)", kTH2D, {ptAxis, nchAxis});
    registry.add("Inclusive/hPtPrimRecoEvVsMult", "Generated primaries (reco events) vs mult;#it{p}_{T};Mult Class (%)", kTH2D, {ptAxis, nchAxis});
    registry.add("Inclusive/hPtPrimGoodEvVsMult", "Generated primaries (good events) vs mult;#it{p}_{T};Mult Class (%)", kTH2D, {ptAxis, nchAxis});
    registry.add("Inclusive/hPtNumEffVsMult", "Tracking efficiency numerator vs mult;#it{p}_{T};Mult Class (%)", kTH2D, {ptAxis, nchAxis});
    registry.add("Inclusive/hPtDenEffVsMult", "Tracking efficiency denominator vs mult;#it{p}_{T};Mult Class (%)", kTH2D, {ptAxis, nchAxis});
    registry.add("Inclusive/hPtAllRecoVsMult", "All reconstructed tracks vs mult;#it{p}_{T};Mult Class (%)", kTH2D, {ptAxis, nchAxis});
    registry.add("Inclusive/hPtPrimRecoVsMult", "Reconstructed primaries vs mult;#it{p}_{T};Mult Class (%)", kTH2D, {ptAxis, nchAxis});
    registry.add("Inclusive/hPtSecRecoVsMult", "Reconstructed secondaries vs mult;#it{p}_{T};Mult Class (%)", kTH2D, {ptAxis, nchAxis});

    // ========================================================================
    // PER-SPECIES INCLUSIVE HISTOGRAMS
    // ========================================================================
    const std::array<std::string, 3> particleNames = {"Pion", "Kaon", "Proton"};
    const std::array<std::string, 3> particleSymbols = {"#pi^{#pm}", "K^{#pm}", "p+#bar{p}"};

    for (int iSpecies = 0; iSpecies < ParticleTypes - 1; ++iSpecies) {
      const auto& name = particleNames[iSpecies];
      const auto& symbol = particleSymbols[iSpecies];

      registry.add(Form("%s/hPtPrimGenAll", name.c_str()),
                   Form("All generated %s (no cuts);#it{p}_{T};Counts", symbol.c_str()), kTH1D, {ptAxis});
      registry.add(Form("%s/hPtPrimBadVertex", name.c_str()),
                   Form("Generated %s (bad vertex);#it{p}_{T};Counts", symbol.c_str()), kTH1D, {ptAxis});
      registry.add(Form("%s/hPtPrimGen", name.c_str()),
                   Form("Generated %s (after physics selection);#it{p}_{T};Counts", symbol.c_str()), kTH1D, {ptAxis});
      registry.add(Form("%s/hPtPrimRecoEv", name.c_str()),
                   Form("Generated %s (reco events);#it{p}_{T};Counts", symbol.c_str()), kTH1D, {ptAxis});
      registry.add(Form("%s/hPtPrimGoodEv", name.c_str()),
                   Form("Generated %s (good events);#it{p}_{T};Counts", symbol.c_str()), kTH1D, {ptAxis});
      registry.add(Form("%s/hPtNumEff", name.c_str()),
                   Form("%s tracking efficiency numerator;#it{p}_{T};Counts", symbol.c_str()), kTH1D, {ptAxis});
      registry.add(Form("%s/hPtDenEff", name.c_str()),
                   Form("%s tracking efficiency denominator;#it{p}_{T};Counts", symbol.c_str()), kTH1D, {ptAxis});
      registry.add(Form("%s/hPtAllReco", name.c_str()),
                   Form("All reconstructed %s;#it{p}_{T};Counts", symbol.c_str()), kTH1D, {ptAxis});
      registry.add(Form("%s/hPtPrimReco", name.c_str()),
                   Form("Reconstructed primary %s;#it{p}_{T};Counts", symbol.c_str()), kTH1D, {ptAxis});
      registry.add(Form("%s/hPtSecReco", name.c_str()),
                   Form("Reconstructed secondary %s;#it{p}_{T};Counts", symbol.c_str()), kTH1D, {ptAxis});
      registry.add(Form("%s/hPtMeasured", name.c_str()),
                   Form("Measured %s;#it{p}_{T};Counts", symbol.c_str()), kTH1D, {ptAxis});

      // Per-species vs multiplicity
      registry.add(Form("%s/hPtPrimGenAllVsMult", name.c_str()),
                   Form("All generated %s vs mult;#it{p}_{T};Mult Class (%s)", symbol.c_str(), ")"), kTH2D, {ptAxis, nchAxis});
      registry.add(Form("%s/hPtPrimBadVertexVsMult", name.c_str()),
                   Form("Generated %s (bad vertex) vs mult;#it{p}_{T};Mult Class (%s)", symbol.c_str(), ")"), kTH2D, {ptAxis, nchAxis});
      registry.add(Form("%s/hPtPrimGenVsMult", name.c_str()),
                   Form("Generated %s (after phys sel) vs mult;#it{p}_{T};Mult Class (%s)", symbol.c_str(), ")"), kTH2D, {ptAxis, nchAxis});
      registry.add(Form("%s/hPtPrimRecoEvVsMult", name.c_str()),
                   Form("Generated %s (reco events) vs mult;#it{p}_{T};Mult Class (%s)", symbol.c_str(), ")"), kTH2D, {ptAxis, nchAxis});
      registry.add(Form("%s/hPtPrimGoodEvVsMult", name.c_str()),
                   Form("Generated %s (good events) vs mult;#it{p}_{T};Mult Class (%s)", symbol.c_str(), ")"), kTH2D, {ptAxis, nchAxis});
      registry.add(Form("%s/hPtNumEffVsMult", name.c_str()),
                   Form("%s tracking eff numerator vs mult;#it{p}_{T};Mult Class (%s)", symbol.c_str(), ")"), kTH2D, {ptAxis, nchAxis});
      registry.add(Form("%s/hPtDenEffVsMult", name.c_str()),
                   Form("%s tracking eff denominator vs mult;#it{p}_{T};Mult Class (%s)", symbol.c_str(), ")"), kTH2D, {ptAxis, nchAxis});
      registry.add(Form("%s/hPtAllRecoVsMult", name.c_str()),
                   Form("All reconstructed %s vs mult;#it{p}_{T};Mult Class (%s)", symbol.c_str(), ")"), kTH2D, {ptAxis, nchAxis});
      registry.add(Form("%s/hPtPrimRecoVsMult", name.c_str()),
                   Form("Reconstructed primary %s vs mult;#it{p}_{T};Mult Class (%s)", symbol.c_str(), ")"), kTH2D, {ptAxis, nchAxis});
      registry.add(Form("%s/hPtSecRecoVsMult", name.c_str()),
                   Form("Reconstructed secondary %s vs mult;#it{p}_{T};Mult Class (%s)", symbol.c_str(), ")"), kTH2D, {ptAxis, nchAxis});
      registry.add(Form("%s/hPtMeasuredVsMult", name.c_str()),
                   Form("Measured %s vs mult;#it{p}_{T};Mult Class (%s)", symbol.c_str(), ")"), kTH2D, {ptAxis, nchAxis});
    }

    // ========================================================================
    // PID HISTOGRAMS
    // ========================================================================
    registry.add("PtResolution", "pT resolution;#it{p}_{T}^{gen};(p_{T}^{reco}-p_{T}^{gen})/p_{T}^{gen}",
                 kTH2F, {{ptAxis, {100, -1.0, 1.0}}});

    if (enablePIDHistograms) {
      registry.add("Pion/hNsigmaTPC", "Pion TPC n#sigma;#it{p}_{T};n#sigma_{TPC}",
                   kTH2D, {{ptAxis, {200, -10, 10}}});
      registry.add("Kaon/hNsigmaTPC", "Kaon TPC n#sigma;#it{p}_{T};n#sigma_{TPC}",
                   kTH2D, {{ptAxis, {200, -10, 10}}});
      registry.add("Proton/hNsigmaTPC", "Proton TPC n#sigma;#it{p}_{T};n#sigma_{TPC}",
                   kTH2D, {{ptAxis, {200, -10, 10}}});
    }

    // ========================================================================
    // PHI CUT MONITORING
    // ========================================================================
    if (applyPhiCut.value) {
      registry.add("PhiCut/hPtVsPhiPrimeBefore", "pT vs #phi' before cut;p_{T};#phi'",
                   kTH2F, {{100, 0, 10}, {100, 0, 0.4}});
      registry.add("PhiCut/hPtVsPhiPrimeAfter", "pT vs #phi' after cut;p_{T};#phi'",
                   kTH2F, {{100, 0, 10}, {100, 0, 0.4}});
    }
    // ========================================================================
    // DCA CUT MONITORING
    // ========================================================================

    registry.add("hDCAxyVsPt_before", "DCAxy vs pT before cut;#it{p}_{T} (GeV/c);DCA_{xy} (cm)",
                 HistType::kTH2F, {{ptAxis}, {dcaAxis}});
    registry.add("hDCAzVsPt_before", "DCAz vs pT before cut;#it{p}_{T} (GeV/c);DCA_{z} (cm)",
                 HistType::kTH2F, {{ptAxis}, {dcaAxis}});
    registry.add("hDCAxyVsPt_after", "DCAxy vs pT after cut;#it{p}_{T} (GeV/c);DCA_{xy} (cm)",
                 HistType::kTH2F, {{ptAxis}, {dcaAxis}});
    registry.add("hDCAzVsPt_after", "DCAz vs pT after cut;#it{p}_{T} (GeV/c);DCA_{z} (cm)",
                 HistType::kTH2F, {{ptAxis}, {dcaAxis}});

    // ========================================================================
    // CALIBRATION HISTOGRAMS
    // ========================================================================
    registry.add("Calibration/hRawMultiplicity", "Raw multiplicity distribution;N_{ch};Events",
                 kTH1D, {{150, 0, 150}});

    // ========================================================================
    // DEDX VS MOMENTUM HISTOGRAMS
    // ========================================================================
    const std::array<std::string, NCentHists + 1> centNames = {
      "Cent0_1", "Cent1_5", "Cent5_10", "Cent10_15", "Cent15_20",
      "Cent20_30", "Cent30_40", "Cent40_50", "Cent50_70", "Cent70_100", "MB"};
    const std::array<std::string, ParticleTypes> v0Names = {
      "all", "Pi_v0", "Pr_v0", "El_v0"};
    for (int i = 0; i < ParticleTypes; ++i) {
      const auto& part = v0Names[i];
      registry.add(Form("DedxVsMomentum/dEdx_vs_Momentum_%s_Pos", part.c_str()),
                   "dE/dx vs Momentum Positive", kTH3F, {{pAxis}, {dedxAxis}, {etaAxis}});
      registry.add(Form("DedxVsMomentum/dEdx_vs_Momentum_%s_Neg", part.c_str()),
                   "dE/dx vs Momentum Negative", kTH3F, {{pAxis}, {dedxAxis}, {etaAxis}});
    }
    for (int i = 0; i < CentralityClasses; ++i) {
      const auto& cent = centNames[i];
      hDedxVsMomentumVsCentPos[i] = registry.add<TH3>(Form("DedxVsMomentum/dEdx_vs_Momentum_%s_Pos", cent.c_str()), "dE/dx vs Momentum Positive", HistType::kTH3F, {{ptAxis}, {dedxAxis}, {etaAxis}});
      hDedxVsMomentumVsCentNeg[i] = registry.add<TH3>(Form("DedxVsMomentum/dEdx_vs_Momentum_%s_Neg", cent.c_str()), "dE/dx vs Momentum Negative", HistType::kTH3F, {{ptAxis}, {dedxAxis}, {etaAxis}});
    }
    for (int i = 0; i < CentralityClasses + 1; ++i) {
      const auto& cent = centNames[i];
      hDedxVspTMomentumVsCent[i] = registry.add<TH3>(Form("DedxVsMomentum/dEdx_vs_pT_%s", cent.c_str()), "dE/dx vs pT", HistType::kTH3F, {{ptAxis}, {dedxAxis}, {etaAxis}});
    }
    // ========================================================================
    // RESPONSE MATRIX HISTOGRAMS
    // ========================================================================
    for (int i = 0; i < ResponseMatrixTypes; ++i) {
      registry.add(("ResponseMatrix/" + std::string(EtavspvspTPosPart[i])).c_str(),
                   "eta vs pT vs p Positive", HistType::kTH3F,
                   {{etaAxis}, {ptAxis}, {pAxis}});
      registry.add(("ResponseMatrix/" + std::string(EtavspvspTNegPart[i])).c_str(),
                   "eta vs pT vs p Negative", HistType::kTH3F,
                   {{etaAxis}, {ptAxis}, {pAxis}});
    }
    // ========================================================================
    // FINNER BINNING HISTOGRAMS
    // ========================================================================
    for (int i = 0; i < CentralityClasses + 1; ++i) {
      const auto& cent = centNames[i];
      hMomentumVsEtaPos[i] = registry.add<TH2>(Form("Binning/p_vs_eta_%s_Pos", cent.c_str()), "p vs eta", HistType::kTH2F, {{etaAxis}, {pFineAxis}});
      hMomentumVsEtaNeg[i] = registry.add<TH2>(Form("Binning/p_vs_eta_%s_Neg", cent.c_str()), "p vs eta", HistType::kTH2F, {{etaAxis}, {pFineAxis}});
      hpTVsEtaPos[i] = registry.add<TH2>(Form("Binning/pT_vs_eta_%s_Pos", cent.c_str()), "pT vs eta", HistType::kTH2F, {{etaAxis}, {pTFineAxis}});
      hpTVsEtaNeg[i] = registry.add<TH2>(Form("Binning/pT_vs_eta_%s_Neg", cent.c_str()), "pT vs eta", HistType::kTH2F, {{etaAxis}, {pTFineAxis}});
    }

    // ========================================================================
    // PARTICLE FRACTIONS HISTOGRAMS
    // ========================================================================
    const std::array<std::string, NPartHists> partName = {"Pion", "Kaon", "Proton", "Electron", "Muon"};
    for (int ic = 0; ic < NCentHists + 1; ++ic) {
      const auto& cent = centNames[ic];
      hTotalMomPosCent[ic] = registryFrac.add<TH2>(
        Form("ParticleFractions/hTotalCountsVsMomentumPos_%s", cent.c_str()),
        "Total counts vs momentum", HistType::kTH2D, {{etaAxis}, {pAxis}});
      hTotalMomNegCent[ic] = registryFrac.add<TH2>(
        Form("ParticleFractions/hTotalCountsVsMomentumNeg_%s", cent.c_str()),
        "Total counts vs momentum", HistType::kTH2D, {{etaAxis}, {pAxis}});
      hTotalPtPosCent[ic] = registryFrac.add<TH2>(
        Form("ParticleFractions/hTotalCountsVsPtPos_%s", cent.c_str()),
        "Total counts vs pT", HistType::kTH2D, {{etaAxis}, {ptAxis}});
      hTotalPtNegCent[ic] = registryFrac.add<TH2>(
        Form("ParticleFractions/hTotalCountsVsPtNeg_%s", cent.c_str()),
        "Total counts vs pT", HistType::kTH2D, {{etaAxis}, {ptAxis}});

      for (int ip = 0; ip < NPartHists; ++ip) {
        const auto& part = partName[ip];
        hFracMomPosCent[ip][ic] = registryFrac.add<TH2>(
          Form("ParticleFractions/hFractionVsMomentum_%s_Pos_%s", part.c_str(), cent.c_str()),
          "Fraction vs momentum", HistType::kTH2D, {{etaAxis}, {pAxis}});
        hFracMomNegCent[ip][ic] = registryFrac.add<TH2>(
          Form("ParticleFractions/hFractionVsMomentum_%s_Neg_%s", part.c_str(), cent.c_str()),
          "Fraction vs momentum", HistType::kTH2D, {{etaAxis}, {pAxis}});
        hFracPtPosCent[ip][ic] = registryFrac.add<TH2>(
          Form("ParticleFractions/hFractionVsPt_%s_Pos_%s", part.c_str(), cent.c_str()),
          "Fraction vs pT", HistType::kTH2D, {{etaAxis}, {ptAxis}});
        hFracPtNegCent[ip][ic] = registryFrac.add<TH2>(
          Form("ParticleFractions/hFractionVsPt_%s_Neg_%s", part.c_str(), cent.c_str()),
          "Fraction vs pT", HistType::kTH2D, {{etaAxis}, {ptAxis}});
      }
    }

    LOG(info) << "=== MultiplicityPt initialized with ALL histograms (including dE/dx) ===";
    LOG(info) << "tpcNchAcceptance = " << tpcNchAcceptance.value;
    LOG(info) << "cfgINELCut       = " << cfgINELCut.value;
    LOG(info) << "selTVXMC         = " << selTVXMC.value;
    LOG(info) << "applyPhiCut      = " << applyPhiCut.value;
    // LOG(info) << "maxDcaZ          = " << maxDcaZ.value;
  }

  // Get magnetic field from CCDB
  int getMagneticField(uint64_t timestamp)
  {
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
    }
    return grpo->getNominalL3Field();
  }

  // Transform phi for phi cut
  float getTransformedPhi(const float phi, const int charge, const float magField) const
  {
    float transformedPhi = phi;
    if (magField < 0) {
      transformedPhi = o2::constants::math::TwoPI - transformedPhi;
    }
    if (charge < 0) {
      transformedPhi = o2::constants::math::TwoPI - transformedPhi;
    }
    transformedPhi += o2::constants::math::PI / 18.0f;
    transformedPhi = std::fmod(transformedPhi, o2::constants::math::PI / 9.0f);
    return transformedPhi;
  }

  // Check phi cut
  template <typename T>
  bool passedPhiCut(const T& track, float magField) const
  {
    if (!applyPhiCut.value)
      return true;
    if (track.pt() < pTthresholdPhiCut.value)
      return true;

    float phiPrime = getTransformedPhi(track.phi(), track.sign(), magField);

    if (phiPrime < fphiCutHigh->Eval(track.pt()) && phiPrime > fphiCutLow->Eval(track.pt())) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool passesTrackSelectionNoDCA(const T& track) const
  {
    if (track.eta() < cfgCutEtaMin.value || track.eta() > cfgCutEtaMax.value)
      return false;
    if (track.pt() < cfgTrkLowPtCut.value)
      return false;
    if (track.tpcChi2NCl() < minChi2PerClusterTPC.value || track.tpcChi2NCl() > maxChi2PerClusterTPC.value)
      return false;

    if (useCustomTrackCuts.value) {
      for (int i = 0; i < static_cast<int>(TrackSelection::TrackCuts::kNCuts); i++) {
        if (i == static_cast<int>(TrackSelection::TrackCuts::kDCAxy))
          continue;
        if (!customTrackCuts.IsSelected(track, static_cast<TrackSelection::TrackCuts>(i)))
          return false;
      }
    } else {
      if (!track.isGlobalTrackWoDCA())
        return false;
    }

    if (nClTPCFoundCut.value && track.tpcNClsFound() < minTPCNClsFound.value)
      return false;
    if (nClTPCPIDCut.value && track.tpcNClsPID() < minTPCNClsPID.value)
      return false;

    return true;
  }

  // DCA xy cut
  template <typename T1>
  bool passesDCAxyCut(const T1& track) const
  {
    const float maxDcaXY = nSigmaDCAxy.value * (dcaXYp0.value + dcaXYp1.value / std::pow(track.pt(), dcaXYp2.value)) / 3.0;
    return std::abs(track.dcaXY()) < maxDcaXY;
  }
  // DCA z cut
  template <typename T1>
  bool passesDCAzCut(const T1& track) const
  {
    const float maxiDcaZ = nSigmaDCAz.value * (dcaZp0.value + dcaZp1.value / std::pow(track.pt(), dcaZp2.value)) / 3.0;
    return std::abs(track.dcaZ()) < maxiDcaZ;
  }

  // Full track selection
  template <typename T>
  bool passesTrackSelection(const T& track) const
  {
    return passesTrackSelectionNoDCA(track) && passesDCAxyCut(track) && passesDCAzCut(track);
  }

  template <typename C>
  bool isEventSelected(const C& col) const
  {
    if (askForCustomTVX.value) {
      if (!col.selection_bit(aod::evsel::kIsTriggerTVX))
        return false;
    } else {
      if (!col.sel8())
        return false;
    }
    if (removeITSROFrameBorder.value && !col.selection_bit(aod::evsel::kNoITSROFrameBorder))
      return false;
    if (removeNoSameBunchPileup.value && !col.selection_bit(aod::evsel::kNoSameBunchPileup))
      return false;
    if (requireIsGoodZvtxFT0vsPV.value && !col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
      return false;
    if (requireIsVertexITSTPC.value && !col.selection_bit(aod::evsel::kIsVertexITSTPC))
      return false;
    if (removeNoTimeFrameBorder.value && !col.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return false;
    if (std::abs(col.posZ()) > cfgCutVertex.value)
      return false;
    return true;
  }

  template <typename C>
  bool isEventSelectedMC(const C& col) const
  {
    if (std::abs(col.posZ()) > cfgCutVertex.value)
      return false;
    return true;
  }

  template <typename T>
  int bestPIDHypothesis(const T& track) const
  {
    const float nsPi = std::abs(track.tpcNSigmaPi());
    const float nsKa = std::abs(track.tpcNSigmaKa());
    const float nsPr = std::abs(track.tpcNSigmaPr());
    float best = 999.f;
    int id = -1;
    if (nsPi < cfgCutNsigmaPi.value && nsPi < best) {
      best = nsPi;
      id = kPion;
    }
    if (nsKa < cfgCutNsigmaKa.value && nsKa < best) {
      best = nsKa;
      id = kKaon;
    }
    if (nsPr < cfgCutNsigmaPr.value && nsPr < best) {
      best = nsPr;
      id = kProton;
    }
    return id;
  }

  int countGeneratedChargedPrimaries(const aod::McParticles& particles, float etaMax) const
  {
    int count = 0;
    for (const auto& particle : particles) {
      auto* pdgPart = pdg->GetParticle(particle.pdgCode());
      if (!pdgPart || std::abs(pdgPart->Charge()) < MinCharge)
        continue;
      if (!particle.isPhysicalPrimary())
        continue;
      if (std::abs(particle.eta()) > etaMax)
        continue;
      if (particle.pt() < cfgTrkLowPtCut.value)
        continue;
      count++;
    }
    return count;
  }

  void processSim(aod::McCollisions::iterator const& mcCollision,
                  soa::SmallGroups<ColEvSelsMC> const& collisions,
                  aod::McParticles const& mcParticles,
                  TracksMC const& tracksMC,
                  BCsRun3 const& /*bcs*/)
  {
    registry.fill(HIST("EventCounter"), kAllGen);

    int nChFT0A = 0, nChFT0C = 0;
    int nChINEL = 0;
    int nChMCEta = 0;
    std::vector<float> particlePtBySpecies[4]; // Pi, Ka, Pr, All
    std::vector<float> particlePtAll;

    // Store particle pT for MC closure
    std::vector<float> mcPiPt, mcKaPt, mcPrPt;

    for (const auto& particle : mcParticles) {
      auto* pdgPart = pdg->GetParticle(particle.pdgCode());
      if (!pdgPart || std::abs(pdgPart->Charge()) < MinCharge)
        continue;
      if (!particle.isPhysicalPrimary())
        continue;

      const float eta = particle.eta();
      const float pt = particle.pt();

      if (eta > MinFT0A && eta < MaxFT0A)
        nChFT0A++;
      if (eta > MinFT0C && eta < MaxFT0C)
        nChFT0C++;
      if (std::abs(eta) < 1.0f)
        nChINEL++;

      if (std::abs(eta) < tpcNchAcceptance.value) {
        nChMCEta++;

        const int absPDG = std::abs(particle.pdgCode());
        if (absPDG == PDG_t::kPiPlus) {
          particlePtBySpecies[kPion].push_back(pt);
          mcPiPt.push_back(pt);
        } else if (absPDG == PDG_t::kKPlus) {
          particlePtBySpecies[kKaon].push_back(pt);
          mcKaPt.push_back(pt);
        } else if (absPDG == PDG_t::kProton) {
          particlePtBySpecies[kProton].push_back(pt);
          mcPrPt.push_back(pt);
        }
        particlePtAll.push_back(pt);
      }
    }

    // Fill NchMCcentVsTVX before TVX selection
    registry.fill(HIST("NchMCcentVsTVX"), nChMCEta, 0.5);

    if (selTVXMC.value && !(nChFT0A > 0 && nChFT0C > 0))
      return;
    registry.fill(HIST("NchMCcentVsTVX"), nChMCEta, 1.5);
    registry.fill(HIST("EventCounter"), kTVXequiv);

    if (isZvtxPosSelMC.value && std::abs(mcCollision.posZ()) > cfgCutVertex.value)
      return;
    registry.fill(HIST("EventCounter"), kVtxZ);

    if (cfgINELCut.value == 1 && nChINEL == 0)
      return;
    if (cfgINELCut.value == INELgt1 && nChINEL < INELgt1)
      return;
    registry.fill(HIST("EventCounter"), kINELgt0);

    const float nchF = static_cast<float>(nChMCEta);

    // Fill event loss denominator and MC closure
    registry.fill(HIST("NchMC_AllGen"), nchF);
    registry.fill(HIST("MC/EventLoss/NchGenerated"), nchF);

    for (const float& pt : mcPiPt) {
      registry.fill(HIST("MCclosure_PtMCPiVsNchMC"), pt, nchF);
      registry.fill(HIST("MC/GenPtVsNch"), pt, nchF);
    }
    for (const float& pt : mcKaPt) {
      registry.fill(HIST("MCclosure_PtMCKaVsNchMC"), pt, nchF);
      registry.fill(HIST("MC/GenPtVsNch"), pt, nchF);
    }
    for (const float& pt : mcPrPt) {
      registry.fill(HIST("MCclosure_PtMCPrVsNchMC"), pt, nchF);
      registry.fill(HIST("MC/GenPtVsNch"), pt, nchF);
    }

    // Fill physics-selected histograms (after vertex and INEL cuts)
    registry.fill(HIST("MC/EventLoss/NchGenerated_PhysicsSelected"), nchF);
    for (const float& pt : mcPiPt) {
      registry.fill(HIST("MC/GenPtVsNch_PhysicsSelected"), pt, nchF);
    }
    for (const float& pt : mcKaPt) {
      registry.fill(HIST("MC/GenPtVsNch_PhysicsSelected"), pt, nchF);
    }
    for (const float& pt : mcPrPt) {
      registry.fill(HIST("MC/GenPtVsNch_PhysicsSelected"), pt, nchF);
    }

    // Fill signal loss denominators
    for (const float& pt : particlePtBySpecies[kPion]) {
      registry.fill(HIST("PtPiVsNchMC_AllGen"), pt, nchF);
    }
    for (const float& pt : particlePtBySpecies[kKaon]) {
      registry.fill(HIST("PtKaVsNchMC_AllGen"), pt, nchF);
    }
    for (const float& pt : particlePtBySpecies[kProton]) {
      registry.fill(HIST("PtPrVsNchMC_AllGen"), pt, nchF);
    }
    for (const float& pt : particlePtAll) {
      registry.fill(HIST("PtAllVsNchMC_AllGen"), pt, nchF);
    }

    // Fill inclusive histograms - all generated
    for (const float& pt : particlePtAll) {
      registry.fill(HIST("Inclusive/hPtPrimGenAll"), pt);
      registry.fill(HIST("Inclusive/hPtPrimGenAllVsMult"), pt, nchF);
    }
    for (const float& pt : mcPiPt) {
      registry.fill(HIST("Pion/hPtPrimGenAll"), pt);
      registry.fill(HIST("Pion/hPtPrimGenAllVsMult"), pt, nchF);
    }
    for (const float& pt : mcKaPt) {
      registry.fill(HIST("Kaon/hPtPrimGenAll"), pt);
      registry.fill(HIST("Kaon/hPtPrimGenAllVsMult"), pt, nchF);
    }
    for (const float& pt : mcPrPt) {
      registry.fill(HIST("Proton/hPtPrimGenAll"), pt);
      registry.fill(HIST("Proton/hPtPrimGenAllVsMult"), pt, nchF);
    }

    for (const float& pt : particlePtAll) {
      registry.fill(HIST("Inclusive/hPtPrimGen"), pt);
      registry.fill(HIST("Inclusive/hPtPrimGenVsMult"), pt, nchF);
    }

    const int nRecColls = collisions.size();
    registry.fill(HIST("NumberOfRecoCollisions"), nRecColls);

    if (nRecColls == 0)
      return;
    registry.fill(HIST("EventCounter"), kRecoColl);

    int biggestNContribs = -1;
    int bestCollisionIndex = -1;
    for (const auto& col : collisions) {
      if (col.numContrib() > biggestNContribs) {
        biggestNContribs = col.numContrib();
        bestCollisionIndex = col.globalIndex();
      }
    }

    for (const auto& collision : collisions) {
      if (collision.globalIndex() != bestCollisionIndex)
        continue;
      if (!isEventSelectedMC(collision))
        continue;

      registry.fill(HIST("EventCounter"), kRecoSelected);

      const float centrality = collision.centFT0M();

      float magField = 0.f;
      if (applyPhiCut.value) {
        const auto& bc = collision.foundBC_as<BCsRun3>();
        magField = static_cast<float>(getMagneticField(bc.timestamp()));
      }

      // Fill centrality and correlation histograms
      registry.fill(HIST("Centrality/hCentRaw"), centrality);
      registry.fill(HIST("NchMCVsCent"), centrality, nchF);
      registry.fill(HIST("Centrality/hCentVsMult"), centrality, nchF);
      registry.fill(HIST("Centrality/hMultVsCent"), nchF, centrality);
      registry.fill(HIST("Centrality/hCentVsVz"), centrality, collision.posZ());
      registry.fill(HIST("Centrality_WRecoEvt"), centrality);
      registry.fill(HIST("Centrality_WRecoEvtWSelCri"), centrality);

      for (const float& pt : particlePtAll) {
        registry.fill(HIST("NchMCVsCentVsPt"), centrality, nchF, pt);
      }

      registry.fill(HIST("NchMC_WithRecoEvt"), nchF);
      registry.fill(HIST("MC/EventLoss/NchGenerated_Reconstructed"), nchF);
      registry.fill(HIST("MC/EventLoss/GenMultVsCent"), centrality, nchF);
      registry.fill(HIST("zPosMC"), mcCollision.posZ());
      registry.fill(HIST("zPosReco"), collision.posZ());
      registry.fill(HIST("T0Ccent"), centrality);

      if (collision.has_foundFT0()) {
        registry.fill(HIST("T0CcentVsFoundFT0"), centrality, 1.5);
      } else {
        registry.fill(HIST("T0CcentVsFoundFT0"), centrality, 0.5);
      }
      if (collision.has_foundFT0() && collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
        registry.fill(HIST("T0CcentVsFoundFT0AndTVX"), centrality, 1.5);
      } else {
        registry.fill(HIST("T0CcentVsFoundFT0AndTVX"), centrality, 0.5);
      }

      // Fill inclusive histograms for reco events
      for (const float& pt : particlePtAll) {
        registry.fill(HIST("Inclusive/hPtPrimRecoEv"), pt);
        registry.fill(HIST("Inclusive/hPtPrimRecoEvVsMult"), pt, nchF);
        registry.fill(HIST("Inclusive/hPtPrimGoodEv"), pt);
        registry.fill(HIST("Inclusive/hPtPrimGoodEvVsMult"), pt, nchF);
      }

      // Fill signal loss numerators and efficiency denominators
      for (const float& pt : particlePtBySpecies[kPion]) {
        registry.fill(HIST("PtPiVsNchMC_WithRecoEvt"), pt, nchF);
        registry.fill(HIST("PtPiVsCentMC_WithRecoEvt"), pt, centrality);
        registry.fill(HIST("PtGenPiVsNchMC_WithRecoEvt"), pt, nchF);
      }
      for (const float& pt : particlePtBySpecies[kKaon]) {
        registry.fill(HIST("PtKaVsNchMC_WithRecoEvt"), pt, nchF);
        registry.fill(HIST("PtKaVsCentMC_WithRecoEvt"), pt, centrality);
        registry.fill(HIST("PtGenKaVsNchMC_WithRecoEvt"), pt, nchF);
      }
      for (const float& pt : particlePtBySpecies[kProton]) {
        registry.fill(HIST("PtPrVsNchMC_WithRecoEvt"), pt, nchF);
        registry.fill(HIST("PtPrVsCentMC_WithRecoEvt"), pt, centrality);
        registry.fill(HIST("PtGenPrVsNchMC_WithRecoEvt"), pt, nchF);
      }
      for (const float& pt : particlePtAll) {
        registry.fill(HIST("PtAllVsNchMC_WithRecoEvt"), pt, nchF);
        registry.fill(HIST("PtAllVsCentMC_WithRecoEvt"), pt, centrality);
        registry.fill(HIST("PtGenAllVsNchMC_WithRecoEvt"), pt, nchF);
      }

      // Fill efficiency denominator histograms
      for (const float& pt : mcPiPt) {
        registry.fill(HIST("Pion/hPtDenEff"), pt);
        registry.fill(HIST("Pion/hPtDenEffVsMult"), pt, nchF);
      }
      for (const float& pt : mcKaPt) {
        registry.fill(HIST("Kaon/hPtDenEff"), pt);
        registry.fill(HIST("Kaon/hPtDenEffVsMult"), pt, nchF);
      }
      for (const float& pt : mcPrPt) {
        registry.fill(HIST("Proton/hPtDenEff"), pt);
        registry.fill(HIST("Proton/hPtDenEffVsMult"), pt, nchF);
      }
      for (const float& pt : particlePtAll) {
        registry.fill(HIST("Inclusive/hPtDenEff"), pt);
        registry.fill(HIST("Inclusive/hPtDenEffVsMult"), pt, nchF);
      }

      const auto& groupedTracks = tracksMC.sliceBy(perCollision, collision.globalIndex());

      for (const auto& track : groupedTracks) {
        if (track.eta() < cfgCutEtaMin.value || track.eta() > cfgCutEtaMax.value)
          continue;
        if (track.pt() < cfgTrkLowPtCut.value)
          continue;
        if (!track.has_mcParticle())
          continue;

        // Before DCA cuts
        registry.fill(HIST("hDCAxyVsPt_before"), track.pt(), track.dcaXY());
        registry.fill(HIST("hDCAzVsPt_before"), track.pt(), track.dcaZ());

        if (applyPhiCut.value && track.pt() >= pTthresholdPhiCut.value) {
          float phiPrime = getTransformedPhi(track.phi(), track.sign(), magField);
          registry.fill(HIST("PhiCut/hPtVsPhiPrimeBefore"), track.pt(), phiPrime);
        }

        const auto& particle = track.mcParticle();
        auto* pdgPart = pdg->GetParticle(particle.pdgCode());
        if (!pdgPart || std::abs(pdgPart->Charge()) < MinCharge)
          continue;

        const bool isPrimary = particle.isPhysicalPrimary();
        const bool isDecay = (!isPrimary) && (particle.getProcess() == TMCProcess::kPDecay);
        const bool isMaterial = (!isPrimary) && (!isDecay);
        (void)isMaterial;
        const int absPDG = std::abs(particle.pdgCode());
        const bool isPi = (absPDG == PDG_t::kPiPlus);
        const bool isKa = (absPDG == PDG_t::kKPlus);
        const bool isPr = (absPDG == PDG_t::kProton);

        float momentum = track.p();
        float tpcSignal = track.tpcSignal();
        float eta = track.eta();
        int charge = track.sign();

        int centIndex = -1;
        for (int j = 0; j < CentralityClasses; ++j) {
          if (centrality >= centBinningStd[j] && centrality < centBinningStd[j + 1]) {
            centIndex = j;
            break;
          }
        }
        if (centIndex == -1)
          continue;

        registry.fill(HIST("hEta"), track.eta());
        registry.fill(HIST("hPhi"), track.phi());
        registry.fill(HIST("EtaVsPhi"), track.eta(), track.phi());

        // Fill TPC cluster histograms
        registry.fill(HIST("hNclFoundTPC"), track.tpcNClsFound());
        registry.fill(HIST("hNclPIDTPC"), track.tpcNClsPID());
        registry.fill(HIST("hNclFoundTPCvsPt"), track.pt(), track.tpcNClsFound());
        registry.fill(HIST("hNclPIDTPCvsPt"), track.pt(), track.tpcNClsPID());

        // Fill inclusive reconstructed histograms
        registry.fill(HIST("Inclusive/hPtAllReco"), track.pt());
        registry.fill(HIST("Inclusive/hPtAllRecoVsMult"), track.pt(), nchF);

        if (isPi) {
          registry.fill(HIST("Pion/hPtAllReco"), track.pt());
          registry.fill(HIST("Pion/hPtAllRecoVsMult"), track.pt(), nchF);
        } else if (isKa) {
          registry.fill(HIST("Kaon/hPtAllReco"), track.pt());
          registry.fill(HIST("Kaon/hPtAllRecoVsMult"), track.pt(), nchF);
        } else if (isPr) {
          registry.fill(HIST("Proton/hPtAllReco"), track.pt());
          registry.fill(HIST("Proton/hPtAllRecoVsMult"), track.pt(), nchF);
        }

        if (passesTrackSelectionNoDCA(track)) {
          if (isPrimary) {
            if (isPi)
              registry.fill(HIST("dcaVsPtPi"), track.pt(), track.dcaXY(), centrality);
            if (isPr)
              registry.fill(HIST("dcaVsPtPr"), track.pt(), track.dcaXY(), centrality);
          } else if (isDecay) {
            if (isPi)
              registry.fill(HIST("dcaVsPtPiDec"), track.pt(), track.dcaXY(), centrality);
            if (isPr)
              registry.fill(HIST("dcaVsPtPrDec"), track.pt(), track.dcaXY(), centrality);
          } else {
            if (isPi)
              registry.fill(HIST("dcaVsPtPiMat"), track.pt(), track.dcaXY(), centrality);
            if (isPr)
              registry.fill(HIST("dcaVsPtPrMat"), track.pt(), track.dcaXY(), centrality);
          }
        }

        if (!passesTrackSelection(track))
          continue;

        // After Trk cuts
        registry.fill(HIST("hDCAxyVsPt_after"), track.pt(), track.dcaXY());
        registry.fill(HIST("hDCAzVsPt_after"), track.pt(), track.dcaZ());

        if (applyPhiCut.value && !passedPhiCut(track, magField))
          continue;

        if (applyPhiCut.value && track.pt() >= pTthresholdPhiCut.value) {
          float phiPrime = getTransformedPhi(track.phi(), track.sign(), magField);
          registry.fill(HIST("PhiCut/hPtVsPhiPrimeAfter"), track.pt(), phiPrime);
        }

        // Fill efficiency numerators (reconstructed primaries)
        if (isPrimary) {
          registry.fill(HIST("PtResolution"), particle.pt(),
                        (track.pt() - particle.pt()) / particle.pt());

          registry.fill(HIST("Inclusive/hPtNumEff"), particle.pt());
          registry.fill(HIST("Inclusive/hPtNumEffVsMult"), particle.pt(), nchF);
          registry.fill(HIST("Inclusive/hPtPrimReco"), track.pt());
          registry.fill(HIST("Inclusive/hPtPrimRecoVsMult"), track.pt(), nchF);

          if (isPi) {
            registry.fill(HIST("Pion/hPtNumEff"), particle.pt());
            registry.fill(HIST("Pion/hPtNumEffVsMult"), particle.pt(), nchF);
            registry.fill(HIST("Pion/hPtPrimReco"), track.pt());
            registry.fill(HIST("Pion/hPtPrimRecoVsMult"), track.pt(), nchF);
          } else if (isKa) {
            registry.fill(HIST("Kaon/hPtNumEff"), particle.pt());
            registry.fill(HIST("Kaon/hPtNumEffVsMult"), particle.pt(), nchF);
            registry.fill(HIST("Kaon/hPtPrimReco"), track.pt());
            registry.fill(HIST("Kaon/hPtPrimRecoVsMult"), track.pt(), nchF);
          } else if (isPr) {
            registry.fill(HIST("Proton/hPtNumEff"), particle.pt());
            registry.fill(HIST("Proton/hPtNumEffVsMult"), particle.pt(), nchF);
            registry.fill(HIST("Proton/hPtPrimReco"), track.pt());
            registry.fill(HIST("Proton/hPtPrimRecoVsMult"), track.pt(), nchF);
          }

          // Fill efficiency numerator histograms
          if (isPi) {
            registry.fill(HIST("PtPiVsCent_WithRecoEvt"), track.pt(), centrality);
            registry.fill(HIST("PtGenPiVsNchMC_RecoTrk"), particle.pt(), nchF);
            registry.fill(HIST("MCclosure_PtPiVsNchMC"), track.pt(), nchF);
          } else if (isKa) {
            registry.fill(HIST("PtKaVsCent_WithRecoEvt"), track.pt(), centrality);
            registry.fill(HIST("PtGenKaVsNchMC_RecoTrk"), particle.pt(), nchF);
            registry.fill(HIST("MCclosure_PtKaVsNchMC"), track.pt(), nchF);
          } else if (isPr) {
            registry.fill(HIST("PtPrVsCent_WithRecoEvt"), track.pt(), centrality);
            registry.fill(HIST("PtGenPrVsNchMC_RecoTrk"), particle.pt(), nchF);
            registry.fill(HIST("MCclosure_PtPrVsNchMC"), track.pt(), nchF);
          }
          registry.fill(HIST("PtAllVsCent_WithRecoEvt"), track.pt(), centrality);
          registry.fill(HIST("PtGenAllVsNchMC_RecoTrk"), particle.pt(), nchF);
        } else {
          // Fill secondary particle histograms
          registry.fill(HIST("Inclusive/hPtSecReco"), track.pt());
          registry.fill(HIST("Inclusive/hPtSecRecoVsMult"), track.pt(), nchF);
          if (isPi) {
            registry.fill(HIST("Pion/hPtSecReco"), track.pt());
            registry.fill(HIST("Pion/hPtSecRecoVsMult"), track.pt(), nchF);
          } else if (isKa) {
            registry.fill(HIST("Kaon/hPtSecReco"), track.pt());
            registry.fill(HIST("Kaon/hPtSecRecoVsMult"), track.pt(), nchF);
          } else if (isPr) {
            registry.fill(HIST("Proton/hPtSecReco"), track.pt());
            registry.fill(HIST("Proton/hPtSecRecoVsMult"), track.pt(), nchF);
          }
        }

        // ====================================================================
        // DEDX VS MOMENTUM HISTOGRAMS FILLING - ALL TRACKS
        // ====================================================================
        hDedxVspTMomentumVsCent[10]->Fill(track.pt(), tpcSignal, eta);
        if (charge > 0) {
          registry.fill(HIST("DedxVsMomentum/dEdx_vs_Momentum_all_Pos"), momentum, tpcSignal, eta);
          hDedxVsMomentumVsCentPos[centIndex]->Fill(momentum, tpcSignal, eta);
          hDedxVspTMomentumVsCent[centIndex]->Fill(track.pt(), tpcSignal, eta);
          hMomentumVsEtaPos[centIndex]->Fill(eta, momentum);
          hMomentumVsEtaPos[10]->Fill(eta, momentum);
          hpTVsEtaPos[centIndex]->Fill(eta, track.pt());
          hpTVsEtaPos[10]->Fill(eta, track.pt());
          registry.fill(HIST("ResponseMatrix/heta_vs_pt_vs_p_all_Pos"), eta, track.pt(), momentum);
        } else {
          registry.fill(HIST("DedxVsMomentum/dEdx_vs_Momentum_all_Neg"), momentum, tpcSignal, eta);
          hDedxVsMomentumVsCentNeg[centIndex]->Fill(momentum, tpcSignal, eta);
          hDedxVspTMomentumVsCent[centIndex]->Fill(track.pt(), tpcSignal, eta);
          hMomentumVsEtaNeg[centIndex]->Fill(eta, momentum);
          hMomentumVsEtaNeg[10]->Fill(eta, momentum);
          hpTVsEtaNeg[centIndex]->Fill(eta, track.pt());
          hpTVsEtaNeg[10]->Fill(eta, track.pt());
          registry.fill(HIST("ResponseMatrix/heta_vs_pt_vs_p_all_Neg"), eta, track.pt(), momentum);
        }

        if (isPrimary) {
          if (charge > 0) {
            registry.fill(HIST("ResponseMatrix/heta_vs_pt_vs_p_all_Pos_Pri"), eta, track.pt(), momentum);
          } else {
            registry.fill(HIST("ResponseMatrix/heta_vs_pt_vs_p_all_Neg_Pri"), eta, track.pt(), momentum);
          }
        }

        // ====================================================================
        // DEDX VS MOMENTUM HISTOGRAMS FILLING - PARTICLE SPECIFIC
        // ====================================================================
        if (track.has_mcParticle() && isPrimary) {
          int pdgCode = std::abs(particle.pdgCode());

          if (charge > 0) {
            registry.fill(HIST("ResponseMatrix/heta_vs_pt_vs_p_all_Pos_Pri_MC"), eta, track.pt(), momentum);
          } else {
            registry.fill(HIST("ResponseMatrix/heta_vs_pt_vs_p_all_Neg_Pri_MC"), eta, track.pt(), momentum);
          }

          if (pdgCode == PDG_t::kPiPlus || pdgCode == PDG_t::kKPlus || pdgCode == PDG_t::kProton ||
              pdgCode == PDG_t::kElectron || pdgCode == PDG_t::kMuonPlus) {
            if (charge > 0) {
              hTotalMomPosCent[centIndex]->Fill(eta, momentum);
              hTotalMomPosCent[10]->Fill(eta, momentum);
              hTotalPtPosCent[centIndex]->Fill(eta, track.pt());
              hTotalPtPosCent[10]->Fill(eta, track.pt());
              registry.fill(HIST("ResponseMatrix/heta_vs_pt_vs_p_all_Pos_Pri_MC_Part"), eta, track.pt(), momentum);
            } else {
              hTotalMomNegCent[centIndex]->Fill(eta, momentum);
              hTotalMomNegCent[10]->Fill(eta, momentum);
              hTotalPtNegCent[centIndex]->Fill(eta, track.pt());
              hTotalPtNegCent[10]->Fill(eta, track.pt());
              registry.fill(HIST("ResponseMatrix/heta_vs_pt_vs_p_all_Neg_Pri_MC_Part"), eta, track.pt(), momentum);
            }
          }

          if (pdgCode == PDG_t::kPiPlus) {
            if (charge > 0) {
              hFracMomPosCent[0][centIndex]->Fill(eta, momentum);
              hFracMomPosCent[0][10]->Fill(eta, momentum);
              hFracPtPosCent[0][centIndex]->Fill(eta, track.pt());
              hFracPtPosCent[0][10]->Fill(eta, track.pt());
              registry.fill(HIST("DedxVsMomentum/dEdx_vs_Momentum_Pi_v0_Pos"), momentum, tpcSignal, eta);
              registry.fill(HIST("ResponseMatrix/heta_vs_pt_vs_p_Pi_Pos"), eta, track.pt(), momentum);
            } else {
              hFracMomNegCent[0][centIndex]->Fill(eta, momentum);
              hFracMomNegCent[0][10]->Fill(eta, momentum);
              hFracPtNegCent[0][centIndex]->Fill(eta, track.pt());
              hFracPtNegCent[0][10]->Fill(eta, track.pt());
              registry.fill(HIST("DedxVsMomentum/dEdx_vs_Momentum_Pi_v0_Neg"), momentum, tpcSignal, eta);
              registry.fill(HIST("ResponseMatrix/heta_vs_pt_vs_p_Pi_Neg"), eta, track.pt(), momentum);
            }
          } else if (pdgCode == PDG_t::kKPlus) {
            if (charge > 0) {
              hFracMomPosCent[1][centIndex]->Fill(eta, momentum);
              hFracMomPosCent[1][10]->Fill(eta, momentum);
              hFracPtPosCent[1][centIndex]->Fill(eta, track.pt());
              hFracPtPosCent[1][10]->Fill(eta, track.pt());
              registry.fill(HIST("ResponseMatrix/heta_vs_pt_vs_p_K_Pos"), eta, track.pt(), momentum);
            } else {
              hFracMomNegCent[1][centIndex]->Fill(eta, momentum);
              hFracMomNegCent[1][10]->Fill(eta, momentum);
              hFracPtNegCent[1][centIndex]->Fill(eta, track.pt());
              hFracPtNegCent[1][10]->Fill(eta, track.pt());
              registry.fill(HIST("ResponseMatrix/heta_vs_pt_vs_p_K_Neg"), eta, track.pt(), momentum);
            }
          } else if (pdgCode == PDG_t::kProton) {
            if (charge > 0) {
              hFracMomPosCent[2][centIndex]->Fill(eta, momentum);
              hFracMomPosCent[2][10]->Fill(eta, momentum);
              hFracPtPosCent[2][centIndex]->Fill(eta, track.pt());
              hFracPtPosCent[2][10]->Fill(eta, track.pt());
              registry.fill(HIST("DedxVsMomentum/dEdx_vs_Momentum_Pr_v0_Pos"), momentum, tpcSignal, eta);
              registry.fill(HIST("ResponseMatrix/heta_vs_pt_vs_p_Pr_Pos"), eta, track.pt(), momentum);
            } else {
              hFracMomNegCent[2][centIndex]->Fill(eta, momentum);
              hFracMomNegCent[2][10]->Fill(eta, momentum);
              hFracPtNegCent[2][centIndex]->Fill(eta, track.pt());
              hFracPtNegCent[2][10]->Fill(eta, track.pt());
              registry.fill(HIST("DedxVsMomentum/dEdx_vs_Momentum_Pr_v0_Neg"), momentum, tpcSignal, eta);
              registry.fill(HIST("ResponseMatrix/heta_vs_pt_vs_p_Pr_Neg"), eta, track.pt(), momentum);
            }
          } else if (pdgCode == PDG_t::kElectron) {
            if (charge > 0) {
              hFracMomPosCent[3][centIndex]->Fill(eta, momentum);
              hFracMomPosCent[3][10]->Fill(eta, momentum);
              hFracPtPosCent[3][centIndex]->Fill(eta, track.pt());
              hFracPtPosCent[3][10]->Fill(eta, track.pt());
              registry.fill(HIST("DedxVsMomentum/dEdx_vs_Momentum_El_v0_Pos"), momentum, tpcSignal, eta);
            } else {
              hFracMomNegCent[3][centIndex]->Fill(eta, momentum);
              hFracMomNegCent[3][10]->Fill(eta, momentum);
              hFracPtNegCent[3][centIndex]->Fill(eta, track.pt());
              hFracPtNegCent[3][10]->Fill(eta, track.pt());
              registry.fill(HIST("DedxVsMomentum/dEdx_vs_Momentum_El_v0_Neg"), momentum, tpcSignal, eta);
            }
          } else if (pdgCode == PDG_t::kMuonPlus) {
            if (charge > 0) {
              hFracMomPosCent[4][centIndex]->Fill(eta, momentum);
              hFracMomPosCent[4][10]->Fill(eta, momentum);
              hFracPtPosCent[4][centIndex]->Fill(eta, track.pt());
              hFracPtPosCent[4][10]->Fill(eta, track.pt());
            } else {
              hFracMomNegCent[4][centIndex]->Fill(eta, momentum);
              hFracMomNegCent[4][10]->Fill(eta, momentum);
              hFracPtNegCent[4][centIndex]->Fill(eta, track.pt());
              hFracPtNegCent[4][10]->Fill(eta, track.pt());
            }
          }
        }

        // PID and measured spectra
        const int species = bestPIDHypothesis(track);

        registry.fill(HIST("Inclusive/hPtMeasured"), track.pt());
        registry.fill(HIST("Inclusive/hPtMeasuredVsCent"), track.pt(), centrality);
        registry.fill(HIST("Inclusive/hPtMeasuredVsMult"), track.pt(), nchF);

        if (species == kPion) {
          registry.fill(HIST("PtPiMeasuredVsCent"), track.pt(), centrality);
          registry.fill(HIST("PtPiMeasuredVsNch"), track.pt(), nchF);
          registry.fill(HIST("Pion/hPtMeasured"), track.pt());
          registry.fill(HIST("Pion/hPtMeasuredVsMult"), track.pt(), nchF);
          if (enablePIDHistograms) {
            registry.fill(HIST("Pion/hNsigmaTPC"), track.pt(), track.tpcNSigmaPi());
          }
        } else if (species == kKaon) {
          registry.fill(HIST("PtKaMeasuredVsCent"), track.pt(), centrality);
          registry.fill(HIST("PtKaMeasuredVsNch"), track.pt(), nchF);
          registry.fill(HIST("Kaon/hPtMeasured"), track.pt());
          registry.fill(HIST("Kaon/hPtMeasuredVsMult"), track.pt(), nchF);
          if (enablePIDHistograms) {
            registry.fill(HIST("Kaon/hNsigmaTPC"), track.pt(), track.tpcNSigmaKa());
          }
        } else if (species == kProton) {
          registry.fill(HIST("PtPrMeasuredVsCent"), track.pt(), centrality);
          registry.fill(HIST("PtPrMeasuredVsNch"), track.pt(), nchF);
          registry.fill(HIST("Proton/hPtMeasured"), track.pt());
          registry.fill(HIST("Proton/hPtMeasuredVsMult"), track.pt(), nchF);
          if (enablePIDHistograms) {
            registry.fill(HIST("Proton/hNsigmaTPC"), track.pt(), track.tpcNSigmaPr());
          }
        }

        registry.fill(HIST("PtAllMeasuredVsCent"), track.pt(), centrality);
        registry.fill(HIST("PtAllMeasuredVsNch"), track.pt(), nchF);
      }
      break;
    }
  }
  PROCESS_SWITCH(MultiplicityPt, processSim, "Process MC simulation", true);

  void processData(CollisionTableData::iterator const& collision,
                   TrackTableData const& tracks,
                   BCsRun3 const& bcs)
  {
    (void)bcs;

    if (!isEventSelected(collision))
      return;

    const float centrality = collision.centFT0M();

    float magField = 0;
    if (applyPhiCut.value) {
      const auto& bc = collision.bc_as<BCsRun3>();
      magField = getMagneticField(bc.timestamp());
    }

    registry.fill(HIST("hvtxZ"), collision.posZ());

    for (const auto& track : tracks) {
      if (track.eta() < cfgCutEtaMin.value || track.eta() > cfgCutEtaMax.value)
        continue;
      if (track.pt() < cfgTrkLowPtCut.value)
        continue;

      if (applyPhiCut.value && track.pt() >= pTthresholdPhiCut.value) {
        float phiPrime = getTransformedPhi(track.phi(), track.sign(), magField);
        registry.fill(HIST("PhiCut/hPtVsPhiPrimeBefore"), track.pt(), phiPrime);
      }

      if (!passesTrackSelection(track))
        continue;

      if (applyPhiCut.value && !passedPhiCut(track, magField))
        continue;

      if (applyPhiCut.value && track.pt() >= pTthresholdPhiCut.value) {
        float phiPrime = getTransformedPhi(track.phi(), track.sign(), magField);
        registry.fill(HIST("PhiCut/hPtVsPhiPrimeAfter"), track.pt(), phiPrime);
      }

      registry.fill(HIST("hEta"), track.eta());
      registry.fill(HIST("hPhi"), track.phi());
      registry.fill(HIST("Inclusive/hPtMeasured"), track.pt());
      registry.fill(HIST("Inclusive/hPtMeasuredVsCent"), track.pt(), centrality);

      const int species = bestPIDHypothesis(track);
      if (species == kPion) {
        registry.fill(HIST("PtPiMeasuredVsCent"), track.pt(), centrality);
        registry.fill(HIST("Pion/hPtMeasured"), track.pt());
        if (enablePIDHistograms) {
          registry.fill(HIST("Pion/hNsigmaTPC"), track.pt(), track.tpcNSigmaPi());
        }
      } else if (species == kKaon) {
        registry.fill(HIST("PtKaMeasuredVsCent"), track.pt(), centrality);
        registry.fill(HIST("Kaon/hPtMeasured"), track.pt());
        if (enablePIDHistograms) {
          registry.fill(HIST("Kaon/hNsigmaTPC"), track.pt(), track.tpcNSigmaKa());
        }
      } else if (species == kProton) {
        registry.fill(HIST("PtPrMeasuredVsCent"), track.pt(), centrality);
        registry.fill(HIST("Proton/hPtMeasured"), track.pt());
        if (enablePIDHistograms) {
          registry.fill(HIST("Proton/hNsigmaTPC"), track.pt(), track.tpcNSigmaPr());
        }
      }

      registry.fill(HIST("PtAllMeasuredVsCent"), track.pt(), centrality);
    }
  }
  PROCESS_SWITCH(MultiplicityPt, processData, "Process data", false);
};

// ============================================================
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiplicityPt>(cfgc)};
}
