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

/// \file   flowSP.cxx
/// \author Noor Koster
/// \since  01/12/2024
/// \brief  task to evaluate flow with respect to spectator plane.

#include "GFWWeights.h"

#include "PWGCF/DataModel/SPTableZDC.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"

#include "TF1.h"
#include "TPDGCode.h"

#include <algorithm>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
// using namespace o2::analysis;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowSP {
  // QA Plots
  O2_DEFINE_CONFIGURABLE(cfgFillEventQA, bool, false, "Fill histograms for event QA");
  O2_DEFINE_CONFIGURABLE(cfgFillTrackQA, bool, false, "Fill histograms for track QA");
  O2_DEFINE_CONFIGURABLE(cfgFillPIDQA, bool, false, "Fill histograms for PID QA");
  O2_DEFINE_CONFIGURABLE(cfgFillEventPlaneQA, bool, false, "Fill histograms for Event Plane QA");
  O2_DEFINE_CONFIGURABLE(cfgFillQABefore, bool, false, "Fill QA histograms before cuts, only for processData");
  // Flags to make and fill histograms
  O2_DEFINE_CONFIGURABLE(cfgFillGeneralV1Histos, bool, true, "Fill histograms for vn analysis");
  O2_DEFINE_CONFIGURABLE(cfgFillMixedHarmonics, bool, true, "Flag to make and fill histos for mixed harmonics");
  O2_DEFINE_CONFIGURABLE(cfgFillEventPlane, bool, false, "Flag to make and fill histos with Event Plane");
  O2_DEFINE_CONFIGURABLE(cfgFillXandYterms, bool, false, "Flag to make and fill histos for with separate x and y terms for SPM");
  O2_DEFINE_CONFIGURABLE(cfgFillChargeDependence, bool, true, "Flag to make and fill histos for charge dependent flow");
  O2_DEFINE_CONFIGURABLE(cfgFillPID, bool, false, "Flag to make and fill histos for PID flow");
  // Centrality Estimators -> standard is FT0C
  O2_DEFINE_CONFIGURABLE(cfgCentFT0Cvariant1, bool, false, "Set centrality estimator to cfgCentFT0Cvariant1");
  O2_DEFINE_CONFIGURABLE(cfgCentFT0M, bool, false, "Set centrality estimator to cfgCentFT0M");
  O2_DEFINE_CONFIGURABLE(cfgCentFV0A, bool, false, "Set centrality estimator to cfgCentFV0A");
  O2_DEFINE_CONFIGURABLE(cfgCentNGlobal, bool, false, "Set centrality estimator to cfgCentNGlobal");
  // Standard selections
  O2_DEFINE_CONFIGURABLE(cfgTrackSelsDCAxy, float, 0.2, "Cut on DCA in the transverse direction (cm)");
  O2_DEFINE_CONFIGURABLE(cfgTrackSelsDCAz, float, 2, "Cut on DCA in the longitudinal direction (cm)");
  O2_DEFINE_CONFIGURABLE(cfgTrackSelsNcls, float, 70, "Cut on number of TPC clusters found");
  O2_DEFINE_CONFIGURABLE(cfgTrackSelsFshcls, float, 0.4, "Cut on fraction of shared TPC clusters found");
  O2_DEFINE_CONFIGURABLE(cfgTrackSelsPtmin, float, 0.2, "minimum pt (GeV/c)");
  O2_DEFINE_CONFIGURABLE(cfgTrackSelsPtmax, float, 10, "maximum pt (GeV/c)");
  O2_DEFINE_CONFIGURABLE(cfgTrackSelsEta, float, 0.8, "eta cut");
  O2_DEFINE_CONFIGURABLE(cfgEvSelsVtxZ, float, 10, "vertex cut (cm)");
  O2_DEFINE_CONFIGURABLE(cfgMagField, float, 99999, "Configurable magnetic field;default CCDB will be queried");
  O2_DEFINE_CONFIGURABLE(cfgCentMin, float, 0, "Minimum cenrality for selected events");
  O2_DEFINE_CONFIGURABLE(cfgCentMax, float, 90, "Maximum cenrality for selected events");
  // NUA and NUE weights
  O2_DEFINE_CONFIGURABLE(cfgFillWeights, bool, true, "Fill NUA weights");
  O2_DEFINE_CONFIGURABLE(cfgFillWeightsPOS, bool, true, "Fill NUA weights only for positive charges");
  O2_DEFINE_CONFIGURABLE(cfgFillWeightsNEG, bool, true, "Fill NUA weights only for negative charges");
  O2_DEFINE_CONFIGURABLE(cfguseNUA1D, bool, false, "Use 1D NUA weights (only phi)");
  O2_DEFINE_CONFIGURABLE(cfguseNUA2D, bool, true, "Use 2D NUA weights (phi and eta)");
  // Additional track Selections
  O2_DEFINE_CONFIGURABLE(cfgTrackSelsUseAdditionalTrackCut, bool, false, "Bool to enable Additional Track Cut");
  O2_DEFINE_CONFIGURABLE(cfgTrackSelsDoDCApt, bool, false, "Apply Pt dependent DCAz cut");
  O2_DEFINE_CONFIGURABLE(cfgTrackSelsDCApt1, float, 0.1, "DcaZ < a * b / pt^1.1 -> this sets a");
  O2_DEFINE_CONFIGURABLE(cfgTrackSelsDCApt2, float, 0.035, "DcaZ < a * b / pt^1.1 -> this sets b");
  O2_DEFINE_CONFIGURABLE(cfgTrackSelsPIDNsigma, float, 2.0, "nSigma cut for PID");
  O2_DEFINE_CONFIGURABLE(cfgTrackSelDoTrackQAvsCent, bool, true, "Do track selection QA plots as function of centrality");
  // Additional event selections
  O2_DEFINE_CONFIGURABLE(cfgEvSelsUseAdditionalEventCut, bool, true, "Bool to enable Additional Event Cut");
  O2_DEFINE_CONFIGURABLE(cfgEvSelsMaxOccupancy, int, 10000, "Maximum occupancy of selected events");
  O2_DEFINE_CONFIGURABLE(cfgEvSelsNoSameBunchPileupCut, bool, true, "kNoSameBunchPileupCut");
  O2_DEFINE_CONFIGURABLE(cfgEvSelsIsGoodZvtxFT0vsPV, bool, true, "kIsGoodZvtxFT0vsPV");
  O2_DEFINE_CONFIGURABLE(cfgEvSelsNoCollInTimeRangeStandard, bool, true, "kNoCollInTimeRangeStandard");
  O2_DEFINE_CONFIGURABLE(cfgEvSelsDoOccupancySel, bool, true, "Bool for event selection on detector occupancy");
  O2_DEFINE_CONFIGURABLE(cfgEvSelsTVXinTRD, bool, false, "Use kTVXinTRD (reject TRD triggered events)");
  O2_DEFINE_CONFIGURABLE(cfgEvSelsIsVertexITSTPC, bool, true, "Selects collisions with at least one ITS-TPC track");
  O2_DEFINE_CONFIGURABLE(cfgEvSelsIsGoodITSLayersAll, bool, true, "Cut time intervals with dead ITS staves");
  // harmonics for v coefficients
  O2_DEFINE_CONFIGURABLE(cfgHarm, int, 1, "Flow harmonic n for ux and uy: (Cos(n*phi), Sin(n*phi))");
  O2_DEFINE_CONFIGURABLE(cfgHarmMixed1, int, 2, "Flow harmonic n for ux and uy in mixed harmonics (MH): (Cos(n*phi), Sin(n*phi))");
  O2_DEFINE_CONFIGURABLE(cfgHarmMixed2, int, 3, "Flow harmonic n for ux and uy in mixed harmonics (MH): (Cos(n*phi), Sin(n*phi))");
  // settings for CCDB data
  O2_DEFINE_CONFIGURABLE(cfgCCDBdir_QQ, std::string, "Users/c/ckoster/ZDC/LHC23_PbPb_pass5/meanQQ/Default", "ccdb dir for average QQ values in 1% centrality bins");
  O2_DEFINE_CONFIGURABLE(cfgCCDBdir_SP, std::string, "", "ccdb dir for average event plane resolution in 1% centrality bins");
  O2_DEFINE_CONFIGURABLE(cfgCCDB_NUA, std::string, "Users/c/ckoster/flowSP/LHC23_PbPb_pass5/Default", "ccdb dir for NUA corrections");
  O2_DEFINE_CONFIGURABLE(cfgCCDB_NUE, std::string, "Users/c/ckoster/flowSP/LHC23_PbPb_pass5/NUE/Default", "ccdb dir for NUE corrections");
  O2_DEFINE_CONFIGURABLE(cfgCCDBdir_centrality, std::string, "", "ccdb dir for Centrality corrections");
  // Confogirable axis
  ConfigurableAxis axisCentrality{"axisCentrality", {20, 0, 100}, "Centrality bins for vn "};
  ConfigurableAxis axisNch = {"axisNch", {400, 0, 4000}, "Global N_{ch}"};
  ConfigurableAxis axisMultpv = {"axisMultpv", {400, 0, 4000}, "N_{ch} (PV)"};
  // Configurables containing vector
  Configurable<std::vector<double>> cfgEvSelsMultPv{"cfgEvSelsMultPv", std::vector<double>{2228.05, -75.5988, 0.976695, -0.00585275, 1.40738e-05, 3795.65, -136.988, 2.12393, -0.017028, 5.78679e-05}, "Multiplicity cuts (PV) first 5 parameters cutLOW last 5 cutHIGH (Default is +-2sigma pass5) "};
  Configurable<std::vector<double>> cfgEvSelsMult{"cfgEvSelsMult", std::vector<double>{1308.86, -41.9314, 0.488423, -0.00248178, 4.71554e-06, 2973.55, -103.092, 1.47673, -0.0106685, 3.29348e-05}, "Multiplicity cuts (Global) first 5 parameters cutLOW last 5 cutHIGH (Default is +-2sigma pass5) "};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgEvSelsVtxZ;
  Filter trackFilter = nabs(aod::track::eta) < cfgTrackSelsEta && aod::track::pt > cfgTrackSelsPtmin&& aod::track::pt < cfgTrackSelsPtmax && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && nabs(aod::track::dcaXY) < cfgTrackSelsDCAxy&& nabs(aod::track::dcaZ) < cfgTrackSelsDCAz;
  Filter trackFilterMC = nabs(aod::mcparticle::eta) < cfgTrackSelsEta && aod::mcparticle::pt > cfgTrackSelsPtmin&& aod::mcparticle::pt < cfgTrackSelsPtmax;
  using GeneralCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::CentNGlobals>;
  using UnfilteredTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;

  using UsedTracks = soa::Filtered<UnfilteredTracks>;
  using ZDCCollisions = soa::Filtered<soa::Join<GeneralCollisions, aod::SPTableZDC>>;

  // For MC Reco and Gen
  using CCs = soa::Filtered<soa::Join<GeneralCollisions, aod::McCollisionLabels>>; // without SPTableZDC
  using CC = CCs::iterator;
  using TCs = soa::Join<UnfilteredTracks, aod::McTrackLabels>;
  using FilteredTCs = soa::Filtered<TCs>;
  using TC = TCs::iterator;
  using MCs = soa::Filtered<aod::McParticles>;

  Preslice<aod::McParticles> partPerMcCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<CCs> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<TCs> trackPerMcParticle = aod::mctracklabel::mcParticleId;
  Preslice<TCs> trackPerCollision = aod::track::collisionId;

  //  Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;

  // struct to hold the correction histos/
  struct Config {
    std::vector<TH1D*> mEfficiency = {};
    std::vector<GFWWeights*> mAcceptance = {};
    std::vector<TH3D*> mAcceptance2D = {};
    bool correctionsLoaded = false;
    int lastRunNumber = 0;

    TProfile* hcorrQQ = nullptr;
    TProfile* hcorrQQx = nullptr;
    TProfile* hcorrQQy = nullptr;
    TProfile* hEvPlaneRes = nullptr;
    TH1D* hCentrality = nullptr;

    bool clQQ = false;
    bool clEvPlaneRes = false;
    bool clCentrality = false;

  } cfg;

  OutputObj<GFWWeights> fWeights{GFWWeights("weights")};
  OutputObj<GFWWeights> fWeightsPOS{GFWWeights("weights_positive")};
  OutputObj<GFWWeights> fWeightsNEG{GFWWeights("weights_negative")};

  HistogramRegistry registry{"registry"};
  HistogramRegistry histos{"QAhistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  // Event selection cuts
  TF1* fPhiCutLow = nullptr;
  TF1* fPhiCutHigh = nullptr;
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  enum SelectionCriteria {
    evSel_FilteredEvent,
    evSel_sel8,
    evSel_occupancy,
    evSel_kTVXinTRD,
    evSel_kNoSameBunchPileup,
    evSel_kIsGoodZvtxFT0vsPV,
    evSel_kNoCollInTimeRangeStandard,
    evSel_kIsVertexITSTPC,
    evSel_MultCuts,
    evSel_kIsGoodITSLayersAll,
    evSel_isSelectedZDC,
    evSel_CentCuts,
    nEventSelections
  };

  enum TrackSelections {
    trackSel_ZeroCharge,
    trackSel_Eta,
    trackSel_Pt,
    trackSel_DCAxy,
    trackSel_DCAz,
    trackSel_GlobalTracks,
    trackSel_NCls,
    trackSel_FshCls,
    trackSel_TPCBoundary,
    trackSel_ParticleWeights,
    nTrackSelections
  };

  enum ChargeType {
    kInclusive,
    kPositive,
    kNegative
  };

  enum FillType {
    kBefore,
    kAfter
  };

  enum ModeType {
    kGen,
    kReco
  };

  enum ParticleType {
    kUnidentified,
    kPions,
    kKaons,
    kProtons
  };

  static constexpr std::string_view Charge[] = {"incl/", "pos/", "neg/"};
  static constexpr std::string_view Species[] = {"", "pion/", "kaon/", "proton/"};
  static constexpr std::string_view Time[] = {"before/", "after/"};

  void init(InitContext const&)
  {

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    AxisSpec axisDCAz = {100, -.5, .5, "DCA_{z} (cm)"};
    AxisSpec axisDCAxy = {100, -.5, .5, "DCA_{xy} (cm)"};
    AxisSpec axisPhiMod = {100, 0, constants::math::PI / 9, "fmod(#varphi,#pi/9)"};
    AxisSpec axisPhi = {60, 0, constants::math::TwoPI, "#varphi"};
    AxisSpec axisEta = {64, -1.6, 1.6, "#eta"};
    AxisSpec axisEtaVn = {8, -.8, .8, "#eta"};
    AxisSpec axisVx = {40, -0.01, 0.01, "v_{x}"};
    AxisSpec axisVy = {40, -0.01, 0.01, "v_{y}"};
    AxisSpec axisVz = {40, -10, 10, "v_{z}"};
    AxisSpec axisCent = {90, 0, 90, "Centrality(%)"};
    AxisSpec axisPhiPlane = {100, -constants::math::PI, constants::math::PI, "#Psi"};
    AxisSpec axisT0c = {70, 0, 100000, "N_{ch} (T0C)"};
    AxisSpec axisT0a = {70, 0, 200000, "N_{ch} (T0A)"};
    AxisSpec axisV0a = {70, 0, 200000, "N_{ch} (V0A)"};
    AxisSpec axisShCl = {40, 0, 1, "Fraction shared cl. TPC"};
    AxisSpec axisCl = {80, 0, 160, "Number of cl. TPC"};
    AxisSpec axisNsigma = {100, -10, 10, "Nsigma for TPC and TOF"};
    AxisSpec axisdEdx = {300, 0, 300, "dEdx for PID"};
    AxisSpec axisBeta = {150, 0, 1.5, "Beta for PID"};

    std::vector<double> ptbinning = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10};
    AxisSpec axisPt = {ptbinning, "#it{p}_{T} GeV/#it{c}"};

    int ptbins = ptbinning.size() - 1;

    histos.add("hEventCount", "Number of Event; Cut; #Events Passed Cut", {HistType::kTH1D, {{nEventSelections, 0, nEventSelections}}});
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_FilteredEvent + 1, "Filtered event");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_sel8 + 1, "Sel8");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_occupancy + 1, "kOccupancy");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kTVXinTRD + 1, "kTVXinTRD");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kNoSameBunchPileup + 1, "kNoSameBunchPileup");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kIsGoodZvtxFT0vsPV + 1, "kIsGoodZvtxFT0vsPV");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kNoCollInTimeRangeStandard + 1, "kNoCollInTimeRangeStandard");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kIsVertexITSTPC + 1, "kIsVertexITSTPC");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_MultCuts + 1, "Multiplicity cuts");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_kIsGoodITSLayersAll + 1, "kkIsGoodITSLayersAll");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_isSelectedZDC + 1, "isSelected");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(evSel_CentCuts + 1, "Cenrality range");

    histos.add("hTrackCount", "Number of Tracks; Cut; #Tracks Passed Cut", {HistType::kTH1D, {{nTrackSelections, 0, nTrackSelections}}});
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(trackSel_Eta + 1, "Eta");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(trackSel_Pt + 1, "Pt");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(trackSel_DCAxy + 1, "DCAxy");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(trackSel_DCAz + 1, "DCAz");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(trackSel_GlobalTracks + 1, "GlobalTracks");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(trackSel_NCls + 1, "nClusters TPC");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(trackSel_FshCls + 1, "Frac. sh. Cls TPC");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(trackSel_TPCBoundary + 1, "TPC Boundary");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(trackSel_ZeroCharge + 1, "Only charged");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(trackSel_ParticleWeights + 1, "Apply weights");

    if (cfgFillWeights) {
      registry.add<TH3>("weights2D/hPhi_Eta_vz", "", kTH3D, {axisPhi, axisEta, axisVz});
      registry.add<TH3>("weights2D/hPhi_Eta_vz_positive", "", kTH3D, {axisPhi, axisEta, axisVz});
      registry.add<TH3>("weights2D/hPhi_Eta_vz_negative", "", kTH3D, {axisPhi, axisEta, axisVz});

      // define output objects
      fWeights->setPtBins(ptbins, &ptbinning[0]);
      fWeights->init(true, false);

      fWeightsPOS->setPtBins(ptbins, &ptbinning[0]);
      fWeightsPOS->init(true, false);

      fWeightsNEG->setPtBins(ptbins, &ptbinning[0]);
      fWeightsNEG->init(true, false);
    }

    if (cfgFillEventQA) {
      histos.add("QA/after/hCentFT0C", " ; Cent FT0C (%); ", {HistType::kTH1D, {axisCent}});
      histos.add("QA/after/hCentFT0M", "; Cent FT0M (%); ", {HistType::kTH1D, {axisCent}});
      histos.add("QA/after/hCentFV0A", "; Cent FV0A (%); ", {HistType::kTH1D, {axisCent}});
      histos.add("QA/after/hCentNGlobal", "; Cent NGlobal (%); ", {HistType::kTH1D, {axisCent}});
      histos.add("QA/after/globalTracks_centT0C", "", {HistType::kTH2D, {axisCent, axisNch}});
      histos.add("QA/after/PVTracks_centT0C", "", {HistType::kTH2D, {axisCent, axisMultpv}});
      histos.add("QA/after/globalTracks_PVTracks", "", {HistType::kTH2D, {axisMultpv, axisNch}});
      histos.add("QA/after/globalTracks_multT0A", "", {HistType::kTH2D, {axisT0a, axisNch}});
      histos.add("QA/after/globalTracks_multV0A", "", {HistType::kTH2D, {axisV0a, axisNch}});
      histos.add("QA/after/multV0A_multT0A", "", {HistType::kTH2D, {axisT0a, axisV0a}});
      histos.add("QA/after/multT0C_centT0C", "", {HistType::kTH2D, {axisCent, axisT0c}});
      histos.add("QA/after/CentFT0C_vs_CentFT0Cvariant1", " ; Cent FT0C (%); Cent FT0Cvariant1 (%) ", {HistType::kTH2D, {axisCent, axisCent}});
      histos.add("QA/after/CentFT0C_vs_CentFT0M", " ; Cent FT0C (%); Cent FT0M (%) ", {HistType::kTH2D, {axisCent, axisCent}});
      histos.add("QA/after/CentFT0C_vs_CentFV0A", " ; Cent FT0C (%); Cent FV0A (%) ", {HistType::kTH2D, {axisCent, axisCent}});
      histos.add("QA/after/CentFT0C_vs_CentNGlobal", " ; Cent FT0C (%); Cent NGlobal (%) ", {HistType::kTH2D, {axisCent, axisCent}});

      if (cfgFillEventPlaneQA && doprocessData) {
        histos.add("QA/after/PsiA_vs_Cent", "", {HistType::kTH2D, {axisPhiPlane, axisCent}});
        histos.add("QA/after/PsiC_vs_Cent", "", {HistType::kTH2D, {axisPhiPlane, axisCent}});
        histos.add("QA/after/PsiFull_vs_Cent", "", {HistType::kTH2D, {axisPhiPlane, axisCent}});
        histos.add("QA/after/PsiA_vs_Vx", "", {HistType::kTH2D, {axisPhiPlane, axisVx}});
        histos.add("QA/after/PsiC_vs_Vx", "", {HistType::kTH2D, {axisPhiPlane, axisVx}});
        histos.add("QA/after/PsiFull_vs_Vx", "", {HistType::kTH2D, {axisPhiPlane, axisVx}});
        histos.add("QA/after/PsiA_vs_Vy", "", {HistType::kTH2D, {axisPhiPlane, axisVy}});
        histos.add("QA/after/PsiC_vs_Vy", "", {HistType::kTH2D, {axisPhiPlane, axisVy}});
        histos.add("QA/after/PsiFull_vs_Vy", "", {HistType::kTH2D, {axisPhiPlane, axisVy}});
        histos.add("QA/after/PsiA_vs_Vz", "", {HistType::kTH2D, {axisPhiPlane, axisVz}});
        histos.add("QA/after/PsiC_vs_Vz", "", {HistType::kTH2D, {axisPhiPlane, axisVz}});
        histos.add("QA/after/PsiFull_vs_Vz", "", {HistType::kTH2D, {axisPhiPlane, axisVz}});
      }

      if (cfgFillQABefore) {
        histos.addClone("QA/after/", "QA/before/");
      }
    }

    if (doprocessData || doprocessMCReco) {

      if (cfgFillTrackQA) {
        histos.add("incl/QA/after/pt_phi", "", {HistType::kTH2D, {axisPt, axisPhiMod}});
        histos.add<TH3>("incl/QA/after/hPhi_Eta_vz", "", kTH3D, {axisPhi, axisEta, axisVz});
        histos.add<TH3>("incl/QA/after/hPhi_Eta_vz_corrected", "", kTH3D, {axisPhi, axisEta, axisVz});
        histos.add<TH2>("incl/QA/after/hDCAxy_pt", "", kTH2D, {axisPt, axisDCAxy});
        histos.add<TH2>("incl/QA/after/hDCAz_pt", "", kTH2D, {axisPt, axisDCAz});
        histos.add("incl/QA/after/hSharedClusters_pt", "", {HistType::kTH2D, {axisPt, axisShCl}});
        histos.add("incl/QA/after/hCrossedRows_pt", "", {HistType::kTH2D, {axisPt, axisCl}});
        histos.add("incl/QA/after/hCrossedRows_vs_SharedClusters", "", {HistType::kTH2D, {axisCl, axisShCl}});

        if (cfgTrackSelDoTrackQAvsCent) {
          histos.add<TH3>("incl/QA/after/hPt_Eta", "", kTH3D, {axisPt, axisEta, axisCent});
          histos.add<TH3>("incl/QA/after/hPt_Eta_uncorrected", "", kTH3D, {axisPt, axisEta, axisCent});
          histos.add<TH3>("incl/QA/after/hPhi_Eta", "", kTH3D, {axisPhi, axisEta, axisCent});
          histos.add<TH3>("incl/QA/after/hPhi_Eta_uncorrected", "", kTH3D, {axisPhi, axisEta, axisCent});
        } else {
          histos.add<TH3>("incl/QA/after/hPhi_Eta_Pt", "", kTH3D, {axisPhi, axisEta, axisPt});
          histos.add<TH3>("incl/QA/after/hPhi_Eta_Pt_corrected", "", kTH3D, {axisPhi, axisEta, axisPt});
        }

        if (cfgFillQABefore)
          histos.addClone("incl/QA/after/", "incl/QA/before/");
      }

      if (cfgFillPIDQA) {
        histos.add<TH2>("hPIDcounts", "", kTH2D, {{{4, 0, 4}, axisPt}});
        histos.get<TH2>(HIST("hPIDcounts"))->GetXaxis()->SetBinLabel(1, "UFO");
        histos.get<TH2>(HIST("hPIDcounts"))->GetXaxis()->SetBinLabel(2, "Pion");
        histos.get<TH2>(HIST("hPIDcounts"))->GetXaxis()->SetBinLabel(3, "Kaon");
        histos.get<TH2>(HIST("hPIDcounts"))->GetXaxis()->SetBinLabel(4, "Proton");

        histos.add("incl/QA/after/hdEdxTPC_pt", "", {HistType::kTH2D, {axisPt, axisdEdx}});
        histos.add("incl/QA/after/hBetaTOF_pt", "", {HistType::kTH2D, {axisPt, axisBeta}});
        histos.add("incl/QA/before/hdEdxTPC_pt", "", {HistType::kTH2D, {axisPt, axisdEdx}});
        histos.add("incl/QA/before/hBetaTOF_pt", "", {HistType::kTH2D, {axisPt, axisBeta}});

        histos.add("incl/pion/QA/after/hNsigmaTPC_pt", "", {HistType::kTH2D, {axisPt, axisNsigma}});
        histos.add("incl/pion/QA/after/hNsigmaTOF_pt", "", {HistType::kTH2D, {axisPt, axisNsigma}});

        histos.add<TH3>("incl/pion/QA/after/hPhi_Eta_vz", "", kTH3D, {axisPhi, axisEta, axisVz});
        histos.add<TH3>("incl/pion/QA/after/hPhi_Eta_vz_corrected", "", kTH3D, {axisPhi, axisEta, axisVz});
        histos.add<TH3>("incl/pion/QA/after/hPhi_Eta_Pt", "", kTH3D, {axisPhi, axisEta, axisPt});
        histos.add<TH3>("incl/pion/QA/after/hPhi_Eta_Pt_corrected", "", kTH3D, {axisPhi, axisEta, axisPt});
        histos.add<TH2>("incl/pion/QA/after/hDCAxy_pt", "", kTH2D, {axisPt, axisDCAxy});
        histos.add<TH2>("incl/pion/QA/after/hDCAz_pt", "", kTH2D, {axisPt, axisDCAz});
        histos.add("incl/pion/QA/after/hSharedClusters_pt", "", {HistType::kTH2D, {axisPt, axisShCl}});
        histos.add("incl/pion/QA/after/hCrossedRows_pt", "", {HistType::kTH2D, {axisPt, axisCl}});
        histos.add("incl/pion/QA/after/hCrossedRows_vs_SharedClusters", "", {HistType::kTH2D, {axisCl, axisShCl}});

        if (cfgFillQABefore) {
          histos.addClone("incl/pion/QA/after/", "incl/pion/QA/before/");
        }

        histos.addClone("incl/pion/", "incl/kaon/");
        histos.addClone("incl/pion/", "incl/proton/");
      }

      if (doprocessMCReco) {
        registry.add("trackMCReco/after/hIsPhysicalPrimary", "", {HistType::kTH2D, {{2, 0, 2}, axisCentrality}});
        registry.add("trackMCReco/hTrackSize_unFiltered", "", {HistType::kTH2D, {{100, 0, 200000}, axisCentrality}});
        registry.add("trackMCReco/hTrackSize_Filtered", "", {HistType::kTH2D, {{100, 0, 20000}, axisCentrality}});
        registry.get<TH2>(HIST("trackMCReco/after/hIsPhysicalPrimary"))->GetXaxis()->SetBinLabel(1, "Secondary");
        registry.get<TH2>(HIST("trackMCReco/after/hIsPhysicalPrimary"))->GetXaxis()->SetBinLabel(2, "Primary");
        registry.add("trackMCReco/after/incl/hPt_hadron", "", {HistType::kTH3D, {axisPt, axisEta, axisCentrality}});
        registry.add("trackMCReco/after/incl/hPt_proton", "", {HistType::kTH3D, {axisPt, axisEta, axisCentrality}});
        registry.add("trackMCReco/after/incl/hPt_pion", "", {HistType::kTH3D, {axisPt, axisEta, axisCentrality}});
        registry.add("trackMCReco/after/incl/hPt_kaon", "", {HistType::kTH3D, {axisPt, axisEta, axisCentrality}});
        // Clone into particles and before/after
        registry.addClone("trackMCReco/after/incl/", "trackMCReco/after/pos/");
        registry.addClone("trackMCReco/after/incl/", "trackMCReco/after/neg/");
        registry.addClone("trackMCReco/after/", "trackMCReco/before/");
      }
      if (doprocessData) {
        registry.add<TProfile>("QQCorrelations/qAqCX", "", kTProfile, {axisCent});
        registry.add<TProfile>("QQCorrelations/qAqCY", "", kTProfile, {axisCent});
        registry.add<TProfile>("QQCorrelations/qAqCXY", "", kTProfile, {axisCent});
        registry.add<TProfile>("QQCorrelations/qAXqCY", "", kTProfile, {axisCent});
        registry.add<TProfile>("QQCorrelations/qAYqCX", "", kTProfile, {axisCent});
        registry.add<TProfile>("QQCorrelations/qAXYqCXY", "", kTProfile, {axisCent});

        if (cfgFillGeneralV1Histos) {
          // track properties per centrality and per eta, pt bin
          registry.add<TProfile2D>("incl/vnC_eta", "", kTProfile2D, {axisEtaVn, axisCentrality});
          registry.add<TProfile2D>("incl/vnA_eta", "", kTProfile2D, {axisEtaVn, axisCentrality});
          registry.add<TProfile2D>("incl/vnC_pt", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile2D>("incl/vnA_pt", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile2D>("incl/vnC_pt_odd", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile2D>("incl/vnA_pt_odd", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile>("incl/vnC_cent_minEta", "", kTProfile, {axisCent});
          registry.add<TProfile>("incl/vnA_cent_minEta", "", kTProfile, {axisCent});
          registry.add<TProfile>("incl/vnC_cent_plusEta", "", kTProfile, {axisCent});
          registry.add<TProfile>("incl/vnA_cent_plusEta", "", kTProfile, {axisCent});
        }
        if (cfgFillPID) {
          registry.add<TProfile2D>("incl/pion/vnC_eta", "", kTProfile2D, {axisEtaVn, axisCentrality});
          registry.add<TProfile2D>("incl/pion/vnA_eta", "", kTProfile2D, {axisEtaVn, axisCentrality});
          registry.add<TProfile2D>("incl/pion/vnC_pt", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile2D>("incl/pion/vnA_pt", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile2D>("incl/pion/vnC_pt_odd", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile2D>("incl/pion/vnA_pt_odd", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile>("incl/pion/vnC_cent_minEta", "", kTProfile, {axisCent});
          registry.add<TProfile>("incl/pion/vnA_cent_minEta", "", kTProfile, {axisCent});
          registry.add<TProfile>("incl/pion/vnC_cent_plusEta", "", kTProfile, {axisCent});
          registry.add<TProfile>("incl/pion/vnA_cent_plusEta", "", kTProfile, {axisCent});
        }
        if (cfgFillXandYterms) {
          registry.add<TProfile2D>("incl/vnAx_eta", "", kTProfile2D, {axisEtaVn, axisCentrality});
          registry.add<TProfile2D>("incl/vnAy_eta", "", kTProfile2D, {axisEtaVn, axisCentrality});
          registry.add<TProfile2D>("incl/vnCx_eta", "", kTProfile2D, {axisEtaVn, axisCentrality});
          registry.add<TProfile2D>("incl/vnCy_eta", "", kTProfile2D, {axisEtaVn, axisCentrality});
          registry.add<TProfile2D>("incl/vnAx_pt", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile2D>("incl/vnAy_pt", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile2D>("incl/vnCx_pt", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile2D>("incl/vnCy_pt", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile2D>("incl/vnCx_pt_odd", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile2D>("incl/vnAx_pt_odd", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile2D>("incl/vnCy_pt_odd", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile2D>("incl/vnAy_pt_odd", "", kTProfile2D, {axisPt, axisCentrality});
          if (cfgFillPID) {
            registry.add<TProfile2D>("incl/pion/vnAx_eta", "", kTProfile2D, {axisEtaVn, axisCentrality});
            registry.add<TProfile2D>("incl/pion/vnAy_eta", "", kTProfile2D, {axisEtaVn, axisCentrality});
            registry.add<TProfile2D>("incl/pion/vnCx_eta", "", kTProfile2D, {axisEtaVn, axisCentrality});
            registry.add<TProfile2D>("incl/pion/vnCy_eta", "", kTProfile2D, {axisEtaVn, axisCentrality});
            registry.add<TProfile2D>("incl/pion/vnAx_pt", "", kTProfile2D, {axisPt, axisCentrality});
            registry.add<TProfile2D>("incl/pion/vnAy_pt", "", kTProfile2D, {axisPt, axisCentrality});
            registry.add<TProfile2D>("incl/pion/vnCx_pt", "", kTProfile2D, {axisPt, axisCentrality});
            registry.add<TProfile2D>("incl/pion/vnCy_pt", "", kTProfile2D, {axisPt, axisCentrality});
            registry.add<TProfile2D>("incl/pion/vnCx_pt_odd", "", kTProfile2D, {axisPt, axisCentrality});
            registry.add<TProfile2D>("incl/pion/vnAx_pt_odd", "", kTProfile2D, {axisPt, axisCentrality});
            registry.add<TProfile2D>("incl/pion/vnCy_pt_odd", "", kTProfile2D, {axisPt, axisCentrality});
            registry.add<TProfile2D>("incl/pion/vnAy_pt_odd", "", kTProfile2D, {axisPt, axisCentrality});
          }
        }
        if (cfgFillMixedHarmonics) {
          registry.add<TProfile2D>("incl/v2/vnAxCxUx_pt_MH", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile2D>("incl/v2/vnAxCyUx_pt_MH", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile2D>("incl/v2/vnAxCyUy_pt_MH", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile2D>("incl/v2/vnAyCxUy_pt_MH", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile>("incl/v2/vnAxCxUx_cent_MH", "", kTProfile, {axisCent});
          registry.add<TProfile>("incl/v2/vnAxCyUx_cent_MH", "", kTProfile, {axisCent});
          registry.add<TProfile>("incl/v2/vnAxCyUy_cent_MH", "", kTProfile, {axisCent});
          registry.add<TProfile>("incl/v2/vnAyCxUy_cent_MH", "", kTProfile, {axisCent});
          registry.add<TProfile2D>("incl/v2/vnAxCxUx_eta_MH", "", kTProfile2D, {axisEtaVn, axisCentrality});
          registry.add<TProfile2D>("incl/v2/vnAxCyUx_eta_MH", "", kTProfile2D, {axisEtaVn, axisCentrality});
          registry.add<TProfile2D>("incl/v2/vnAxCyUy_eta_MH", "", kTProfile2D, {axisEtaVn, axisCentrality});
          registry.add<TProfile2D>("incl/v2/vnAyCxUy_eta_MH", "", kTProfile2D, {axisEtaVn, axisCentrality});
          registry.addClone("incl/v2/", "incl/v3/");
          if (cfgFillPID) {
            registry.add<TProfile2D>("incl/pion/v2/vnAxCxUx_pt_MH", "", kTProfile2D, {axisPt, axisCentrality});
            registry.add<TProfile2D>("incl/pion/v2/vnAxCyUx_pt_MH", "", kTProfile2D, {axisPt, axisCentrality});
            registry.add<TProfile2D>("incl/pion/v2/vnAxCyUy_pt_MH", "", kTProfile2D, {axisPt, axisCentrality});
            registry.add<TProfile2D>("incl/pion/v2/vnAyCxUy_pt_MH", "", kTProfile2D, {axisPt, axisCentrality});
            registry.add<TProfile>("incl/pion/v2/vnAxCxUx_cent_MH", "", kTProfile, {axisCent});
            registry.add<TProfile>("incl/pion/v2/vnAxCyUx_cent_MH", "", kTProfile, {axisCent});
            registry.add<TProfile>("incl/pion/v2/vnAxCyUy_cent_MH", "", kTProfile, {axisCent});
            registry.add<TProfile>("incl/pion/v2/vnAyCxUy_cent_MH", "", kTProfile, {axisCent});
            registry.add<TProfile2D>("incl/pion/v2/vnAxCxUx_eta_MH", "", kTProfile2D, {axisEtaVn, axisCentrality});
            registry.add<TProfile2D>("incl/pion/v2/vnAxCyUx_eta_MH", "", kTProfile2D, {axisEtaVn, axisCentrality});
            registry.add<TProfile2D>("incl/pion/v2/vnAxCyUy_eta_MH", "", kTProfile2D, {axisEtaVn, axisCentrality});
            registry.add<TProfile2D>("incl/pion/v2/vnAyCxUy_eta_MH", "", kTProfile2D, {axisEtaVn, axisCentrality});
            registry.addClone("incl/pion/v2/", "incl/pion/v3/");
          }
        }
        if (cfgFillEventPlane) {
          registry.add<TProfile>("incl/vnA_cent_EP", "", kTProfile, {axisCent});
          registry.add<TProfile>("incl/vnC_cent_EP", "", kTProfile, {axisCent});
          registry.add<TProfile>("incl/vnFull_cent_EP", "", kTProfile, {axisCent});
          registry.add<TProfile2D>("incl/vnA_pt_EP", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile2D>("incl/vnC_pt_EP", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile2D>("incl/vnFull_pt_EP", "", kTProfile2D, {axisPt, axisCentrality});
          registry.add<TProfile2D>("incl/vnA_eta_EP", "", kTProfile2D, {axisEtaVn, axisCentrality});
          registry.add<TProfile2D>("incl/vnC_eta_EP", "", kTProfile2D, {axisEtaVn, axisCentrality});
          registry.add<TProfile2D>("incl/vnFull_eta_EP", "", kTProfile2D, {axisEtaVn, axisCentrality});
        }
        if (cfgFillEventPlaneQA) {
          histos.add<TH1>("QA/hSPplaneA", "hSPplaneA", kTH1D, {axisPhiPlane});
          histos.add<TH1>("QA/hSPplaneC", "hSPplaneC", kTH1D, {axisPhiPlane});
          histos.add<TH1>("QA/hSPplaneFull", "hSPplaneFull", kTH1D, {axisPhiPlane});
          histos.add<TProfile>("QA/hCosPhiACosPhiC", "hCosPhiACosPhiC; Centrality(%); #LT Cos(#Psi^{A})Cos(#Psi^{C})#GT", kTProfile, {axisCent});
          histos.add<TProfile>("QA/hSinPhiASinPhiC", "hSinPhiASinPhiC; Centrality(%); #LT Sin(#Psi^{A})Sin(#Psi^{C})#GT", kTProfile, {axisCent});
          histos.add<TProfile>("QA/hSinPhiACosPhiC", "hSinPhiACosPhiC; Centrality(%); #LT Sin(#Psi^{A})Cos(#Psi^{C})#GT", kTProfile, {axisCent});
          histos.add<TProfile>("QA/hCosPhiASinsPhiC", "hCosPhiASinsPhiC; Centrality(%); #LT Cos(#Psi^{A})Sin(#Psi^{C})#GT", kTProfile, {axisCent});
          histos.add<TProfile>("QA/hFullEvPlaneRes", "hFullEvPlaneRes; Centrality(%); -#LT Cos(#Psi^{A} - #Psi^{C})#GT ", kTProfile, {axisCent});
        }
        if (cfgFillEventQA) {
          histos.add("QA/hCentFull", " ; Centrality (%); ", {HistType::kTH1D, {axisCent}});
        }
      } // end of doprocessData
      if (cfgFillChargeDependence || cfgFillPID) {
        registry.addClone("incl/pion/", "incl/proton/");
        registry.addClone("incl/pion/", "incl/kaon/");
        registry.addClone("incl/", "pos/");
        registry.addClone("incl/", "neg/");
      }
      if (cfgFillPIDQA) {
        histos.addClone("incl/", "pos/");
        histos.addClone("incl/", "neg/");
      }

    } else if (doprocessMCGen) {
      registry.add("trackMCGen/nCollReconstructedPerMcCollision", "", {HistType::kTH1D, {{10, -5, 5}}});
      registry.add("trackMCGen/after/incl/hPt_hadron", "", {HistType::kTH3D, {axisPt, axisEta, axisCentrality}});
      registry.add("trackMCGen/after/incl/hPt_proton", "", {HistType::kTH3D, {axisPt, axisEta, axisCentrality}});
      registry.add("trackMCGen/after/incl/hPt_pion", "", {HistType::kTH3D, {axisPt, axisEta, axisCentrality}});
      registry.add("trackMCGen/after/incl/hPt_kaon", "", {HistType::kTH3D, {axisPt, axisEta, axisCentrality}});
      registry.add("trackMCGen/after/incl/phi_eta_vtxZ_gen", "", {HistType::kTH3D, {axisPhi, axisEta, axisVz}});
      registry.addClone("trackMCGen/after/incl/", "trackMCGen/after/pos/");
      registry.addClone("trackMCGen/after/incl/", "trackMCGen/after/neg/");
      registry.addClone("trackMCGen/after/", "trackMCGen/before/");
    }

    if (cfgEvSelsUseAdditionalEventCut) {
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);

      std::vector<double> paramsMultPVCut = cfgEvSelsMultPv;
      std::vector<double> paramsMultCut = cfgEvSelsMult;

      // number of parameters required in cfgEvSelsMultPv and cfgEvSelsMult.  (5 Low + 5 High)
      uint64_t nParams = 10;

      if (paramsMultPVCut.size() < nParams) {
        LOGF(fatal, "cfgEvSelsMultPv not set properly.. size = %d (should be 10) --> Check your config files!", paramsMultPVCut.size());
      } else if (paramsMultCut.size() < nParams) {
        LOGF(fatal, "cfgEvSelsMult not set properly.. size = %d (should be 10) --> Check your config files!", paramsMultCut.size());
      } else {
        fMultPVCutLow->SetParameters(paramsMultPVCut[0], paramsMultPVCut[1], paramsMultPVCut[2], paramsMultPVCut[3], paramsMultPVCut[4]);
        fMultPVCutHigh->SetParameters(paramsMultPVCut[5], paramsMultPVCut[6], paramsMultPVCut[7], paramsMultPVCut[8], paramsMultPVCut[9]);
        fMultCutLow->SetParameters(paramsMultCut[0], paramsMultCut[1], paramsMultCut[2], paramsMultCut[3], paramsMultCut[4]);
        fMultCutHigh->SetParameters(paramsMultCut[5], paramsMultCut[6], paramsMultCut[7], paramsMultCut[8], paramsMultCut[9]);
      }
    }

    if (cfgTrackSelsUseAdditionalTrackCut) {
      fPhiCutLow = new TF1("fPhiCutLow", "0.06/x+pi/18.0-0.06", 0, 100);
      fPhiCutHigh = new TF1("fPhiCutHigh", "0.1/x+pi/18.0+0.06", 0, 100);
    }
  } // end of init

  float getNUA2D(TH3D* hNUA, float eta, float phi, float vtxz)
  {
    int xind = hNUA->GetXaxis()->FindBin(phi);
    int etaind = hNUA->GetYaxis()->FindBin(eta);
    int vzind = hNUA->GetZaxis()->FindBin(vtxz);
    float weight = hNUA->GetBinContent(xind, etaind, vzind);
    if (weight != 0)
      return 1. / weight;
    return 1;
  }

  template <typename TrackObject>
  int getTrackPID(TrackObject track)
  {

    float usedNSigmaPi = -1;
    float usedNSigmaKa = -1;
    float usedNSigmaPr = -1;

    if (track.hasTOF() && track.hasTPC()) {
      usedNSigmaPi = std::hypot(track.tofNSigmaPi(), track.tpcNSigmaPi());
      usedNSigmaKa = std::hypot(track.tofNSigmaKa(), track.tpcNSigmaKa());
      usedNSigmaPr = std::hypot(track.tofNSigmaPr(), track.tpcNSigmaPr());
    } else if (track.hasTOF()) {
      usedNSigmaPi = track.tofNSigmaPi();
      usedNSigmaKa = track.tofNSigmaKa();
      usedNSigmaPr = track.tofNSigmaPr();
    } else if (track.hasTPC()) {
      usedNSigmaPi = track.tpcNSigmaPi();
      usedNSigmaKa = track.tpcNSigmaKa();
      usedNSigmaPr = track.tpcNSigmaPr();
    } else {
      return kUnidentified; // No PID information available
    }

    std::unordered_map<float, int> usedNSigma = {{usedNSigmaPi, kPions}, {usedNSigmaKa, kKaons}, {usedNSigmaPr, kProtons}};

    int nIdentified = 0;
    int valPID = 0;

    for (const auto& nsigma : usedNSigma) {
      if (std::abs(nsigma.first) < cfgTrackSelsPIDNsigma) {
        valPID = nsigma.second;
        nIdentified++;
      }
    }

    if (nIdentified == 0) {
      return kUnidentified; // No PID match found
    } else if (nIdentified == 1) {
      return valPID;
    } else {
      return kUnidentified; // Multiple PID matches found
    }

    return -1;
  }

  int getMagneticField(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  std::pair<float, uint16_t> getCrossingAngleCCDB(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    auto grpo = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", timestamp);
    if (grpo == nullptr) {
      LOGF(fatal, "GRP object for Crossing Angle not found for timestamp %llu", timestamp);
      return {0, 0};
    }
    float crossingAngle = grpo->getCrossingAngle();
    uint16_t crossingAngleTime = grpo->getCrossingAngleTime();
    return {crossingAngle, crossingAngleTime};
  }

  // From Generic Framework
  void loadCorrections(uint64_t timestamp)
  {
    // corrections saved on CCDB as TList {incl, pos, neg} of GFWWeights (acc) TH1D (eff) objects!
    if (cfg.correctionsLoaded)
      return;

    int nWeights = 3;

    if (cfguseNUA1D) {
      if (cfgCCDB_NUA.value.empty() == false) {
        TList* listCorrections = ccdb->getForTimeStamp<TList>(cfgCCDB_NUA, timestamp);
        cfg.mAcceptance.push_back(reinterpret_cast<GFWWeights*>(listCorrections->FindObject("weights")));
        cfg.mAcceptance.push_back(reinterpret_cast<GFWWeights*>(listCorrections->FindObject("weights_positive")));
        cfg.mAcceptance.push_back(reinterpret_cast<GFWWeights*>(listCorrections->FindObject("weights_negative")));
        int sizeAcc = cfg.mAcceptance.size();
        if (sizeAcc < nWeights)
          LOGF(fatal, "Could not load acceptance weights from %s", cfgCCDB_NUA.value.c_str());
        else
          LOGF(info, "Loaded acceptance weights from %s", cfgCCDB_NUA.value.c_str());
      } else {
        LOGF(info, "cfgCCDB_NUA empty! No corrections loaded");
      }
    } else if (cfguseNUA2D) {
      if (cfgCCDB_NUA.value.empty() == false) {
        TH3D* hNUA2D = ccdb->getForTimeStamp<TH3D>(cfgCCDB_NUA, timestamp);
        if (!hNUA2D) {
          LOGF(fatal, "Could not load acceptance weights from %s", cfgCCDB_NUA.value.c_str());
        } else {
          LOGF(info, "Loaded acceptance weights from %s", cfgCCDB_NUA.value.c_str());
          cfg.mAcceptance2D.push_back(hNUA2D);
        }
      } else {
        LOGF(info, "cfgCCDB_NUA empty! No corrections loaded");
      }
    }
    // Get Efficiency correction
    if (cfgCCDB_NUE.value.empty() == false) {
      TList* listCorrections = ccdb->getForTimeStamp<TList>(cfgCCDB_NUE, timestamp);
      cfg.mEfficiency.push_back(reinterpret_cast<TH1D*>(listCorrections->FindObject("Efficiency")));
      cfg.mEfficiency.push_back(reinterpret_cast<TH1D*>(listCorrections->FindObject("Efficiency_pos")));
      cfg.mEfficiency.push_back(reinterpret_cast<TH1D*>(listCorrections->FindObject("Efficiency_neg")));
      int sizeEff = cfg.mEfficiency.size();
      if (sizeEff < nWeights)
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgCCDB_NUE.value.c_str());
      else
        LOGF(info, "Loaded efficiency histogram from %s", cfgCCDB_NUE.value.c_str());
    } else {
      LOGF(info, "cfgCCDB_NUE empty! No corrections loaded");
    }
    cfg.correctionsLoaded = true;
  }

  // From Generic Framework
  bool setCurrentParticleWeights(int pID, float& weight_nue, float& weight_nua, const float& phi, const float& eta, const float& pt, const float& vtxz)
  {
    float eff = 1.;
    int sizeEff = cfg.mEfficiency.size();
    if (sizeEff > pID)
      eff = cfg.mEfficiency[pID]->GetBinContent(cfg.mEfficiency[pID]->FindBin(pt));
    else
      eff = 1.0;
    if (eff == 0)
      return false;
    weight_nue = 1. / eff;

    if (cfguseNUA1D) {
      int sizeAcc = cfg.mAcceptance.size();
      if (sizeAcc > pID) {
        weight_nua = cfg.mAcceptance[pID]->getNUA(phi, eta, vtxz);
      } else {
        weight_nua = 1;
      }
    } else if (cfguseNUA2D) {
      if (cfg.mAcceptance2D.size() > 0) {
        weight_nua = getNUA2D(cfg.mAcceptance2D[0], eta, phi, vtxz);
      } else {
        weight_nua = 1;
      }
    }
    return true;
  }

  template <typename TCollision>
  bool eventSelected(TCollision collision, const int& multTrk)
  {
    if (!collision.sel8())
      return 0;
    histos.fill(HIST("hEventCount"), evSel_sel8);

    // Occupancy
    if (cfgEvSelsDoOccupancySel) {
      auto occupancy = collision.trackOccupancyInTimeRange();
      if (occupancy > cfgEvSelsMaxOccupancy) {
        return 0;
      }
      histos.fill(HIST("hEventCount"), evSel_occupancy);
    }

    if (cfgEvSelsTVXinTRD) {
      if (collision.alias_bit(kTVXinTRD)) {
        // TRD triggered
        // "CMTVX-B-NOPF-TRD,minbias_TVX"
        return 0;
      }
      histos.fill(HIST("hEventCount"), evSel_kTVXinTRD);
    }

    if (cfgEvSelsNoSameBunchPileupCut) {
      if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        // rejects collisions which are associated with the same "found-by-T0" bunch crossing
        // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
        return 0;
      }
      histos.fill(HIST("hEventCount"), evSel_kNoSameBunchPileup);
    }
    if (cfgEvSelsIsGoodZvtxFT0vsPV) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
        // use this cut at low multiplicities with caution
        return 0;
      }
      histos.fill(HIST("hEventCount"), evSel_kIsGoodZvtxFT0vsPV);
    }
    if (cfgEvSelsNoCollInTimeRangeStandard) {
      if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        //  Rejection of the collisions which have other events nearby
        return 0;
      }
      histos.fill(HIST("hEventCount"), evSel_kNoCollInTimeRangeStandard);
    }

    if (cfgEvSelsIsVertexITSTPC) {
      if (!collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        // selects collisions with at least one ITS-TPC track, and thus rejects vertices built from ITS-only tracks
        return 0;
      }
      histos.fill(HIST("hEventCount"), evSel_kIsVertexITSTPC);
    }

    if (cfgEvSelsUseAdditionalEventCut) {
      float vtxz = -999;
      if (collision.numContrib() > 1) {
        vtxz = collision.posZ();
        float zRes = std::sqrt(collision.covZZ());
        float minzRes = 0.25;
        int maxNumContrib = 20;
        if (zRes > minzRes && collision.numContrib() < maxNumContrib)
          vtxz = -999;
      }

      auto multNTracksPV = collision.multNTracksPV();

      if (vtxz > cfgEvSelsVtxZ || vtxz < -cfgEvSelsVtxZ)
        return 0;
      if (multNTracksPV < fMultPVCutLow->Eval(collision.centFT0C()))
        return 0;
      if (multNTracksPV > fMultPVCutHigh->Eval(collision.centFT0C()))
        return 0;
      if (multTrk < fMultCutLow->Eval(collision.centFT0C()))
        return 0;
      if (multTrk > fMultCutHigh->Eval(collision.centFT0C()))
        return 0;

      histos.fill(HIST("hEventCount"), evSel_MultCuts);
    }

    if (cfgEvSelsIsGoodITSLayersAll) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
        // New event selection bits to cut time intervals with dead ITS staves
        // https://indico.cern.ch/event/1493023/ (09-01-2025)
        return 0;
      }
      histos.fill(HIST("hEventCount"), evSel_kIsGoodITSLayersAll);
    }

    return 1;
  }

  template <typename TrackObject>
  bool trackSelected(TrackObject track, const int& field)
  {
    if (std::fabs(track.eta()) > cfgTrackSelsEta)
      return false;
    histos.fill(HIST("hTrackCount"), trackSel_Eta);

    if (track.pt() < cfgTrackSelsPtmin || track.pt() > cfgTrackSelsPtmax)
      return false;

    histos.fill(HIST("hTrackCount"), trackSel_Pt);

    if (track.dcaXY() > cfgTrackSelsDCAxy)
      return false;

    histos.fill(HIST("hTrackCount"), trackSel_DCAxy);

    if (track.dcaZ() > cfgTrackSelsDCAz)
      return false;

    if (cfgTrackSelsDoDCApt && std::fabs(track.dcaZ()) > (cfgTrackSelsDCApt1 * cfgTrackSelsDCApt2) / (std::pow(track.pt(), 1.1)))
      return false;

    histos.fill(HIST("hTrackCount"), trackSel_DCAz);

    if (track.tpcNClsFound() < cfgTrackSelsNcls)
      return false;
    histos.fill(HIST("hTrackCount"), trackSel_NCls);

    if (track.tpcFractionSharedCls() > cfgTrackSelsFshcls)
      return false;
    histos.fill(HIST("hTrackCount"), trackSel_FshCls);

    double phimodn = track.phi();
    if (field < 0) // for negative polarity field
      phimodn = o2::constants::math::TwoPI - phimodn;
    if (track.sign() < 0) // for negative charge
      phimodn = o2::constants::math::TwoPI - phimodn;
    if (phimodn < 0)
      LOGF(warning, "phi < 0: %g", phimodn);

    phimodn += o2::constants::math::PI / 18.0; // to center gap in the middle
    phimodn = fmod(phimodn, o2::constants::math::PI / 9.0);
    if (cfgFillTrackQA && cfgFillQABefore)
      histos.fill(HIST("incl/QA/before/pt_phi"), track.pt(), phimodn);

    if (cfgTrackSelsUseAdditionalTrackCut) {
      if (phimodn < fPhiCutHigh->Eval(track.pt()) && phimodn > fPhiCutLow->Eval(track.pt()))
        return false; // reject track
    }
    if (cfgFillTrackQA)
      histos.fill(HIST("incl/QA/after/pt_phi"), track.pt(), phimodn);
    histos.fill(HIST("hTrackCount"), trackSel_TPCBoundary);
    return true;
  }

  template <FillType ft, typename CollisionObject, typename TracksObject>
  inline void fillEventQA(CollisionObject collision, TracksObject tracks, double centWeight = 1.0)
  {
    if (!cfgFillEventQA)
      return;

    static constexpr std::string_view Time[] = {"before", "after"};

    histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hCentFT0C"), collision.centFT0C(), centWeight);
    histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hCentNGlobal"), collision.centNGlobal(), centWeight);
    histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hCentFT0M"), collision.centFT0M(), centWeight);
    histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/hCentFV0A"), collision.centFV0A(), centWeight);
    histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/globalTracks_centT0C"), collision.centFT0C(), tracks.size(), centWeight);
    histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/PVTracks_centT0C"), collision.centFT0C(), collision.multNTracksPV(), centWeight);
    histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/globalTracks_PVTracks"), collision.multNTracksPV(), tracks.size(), centWeight);
    histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/globalTracks_multT0A"), collision.multFT0A(), tracks.size(), centWeight);
    histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/globalTracks_multV0A"), collision.multFV0A(), tracks.size(), centWeight);
    histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/multV0A_multT0A"), collision.multFT0A(), collision.multFV0A(), centWeight);
    histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/multT0C_centT0C"), collision.centFT0C(), collision.multFT0C(), centWeight);
    histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/CentFT0C_vs_CentFT0Cvariant1"), collision.centFT0C(), collision.centFT0CVariant1(), centWeight);
    histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/CentFT0C_vs_CentFT0M"), collision.centFT0C(), collision.centFT0M(), centWeight);
    histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/CentFT0C_vs_CentFV0A"), collision.centFT0C(), collision.centFV0A(), centWeight);
    histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/CentFT0C_vs_CentNGlobal"), collision.centFT0C(), collision.centNGlobal(), centWeight);

    if (cfgFillEventPlaneQA) {
      if constexpr (o2::framework::has_type_v<aod::sptablezdc::Vx, typename CollisionObject::all_columns>) {
        double psiA = 1.0 * std::atan2(collision.qyA(), collision.qxA());
        double psiC = 1.0 * std::atan2(collision.qyC(), collision.qxC());
        double psiFull = 1.0 * std::atan2(collision.qyA() + collision.qyC(), collision.qxA() + collision.qxC());

        histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/PsiA_vs_Cent"), psiA, collision.centFT0C(), centWeight);
        histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/PsiC_vs_Cent"), psiC, collision.centFT0C(), centWeight);
        histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/PsiFull_vs_Cent"), psiFull, collision.centFT0C(), centWeight);
        histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/PsiA_vs_Vx"), psiA, collision.vx(), centWeight);
        histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/PsiC_vs_Vx"), psiC, collision.vx(), centWeight);
        histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/PsiFull_vs_Vx"), psiFull, collision.vx(), centWeight);
        histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/PsiA_vs_Vy"), psiA, collision.vy(), centWeight);
        histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/PsiC_vs_Vy"), psiC, collision.vy(), centWeight);
        histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/PsiFull_vs_Vy"), psiFull, collision.vy(), centWeight);
        histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/PsiA_vs_Vz"), psiA, collision.posZ(), centWeight);
        histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/PsiC_vs_Vz"), psiC, collision.posZ(), centWeight);
        histos.fill(HIST("QA/") + HIST(Time[ft]) + HIST("/PsiFull_vs_Vz"), psiFull, collision.posZ(), centWeight);
      }
    }
    return;
  }

  template <ChargeType ct, ParticleType pt, typename TrackObject>
  inline void fillHistograms(TrackObject track, float wacc, float weff, double centWeight, double ux, double uy, double uxMH, double uyMH, double uxMH2, double uyMH2, double qxA, double qyA, double qxC, double qyC, double corrQQx, double corrQQy, double corrQQ, double vnA, double vnC, double vnFull, double centrality)
  {

    double weight = wacc * weff * centWeight;

    if (cfgFillGeneralV1Histos) {
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnA_eta"), track.eta(), centrality, (uy * qyA + ux * qxA) / std::sqrt(std::fabs(corrQQ)), weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnC_eta"), track.eta(), centrality, (uy * qyC + ux * qxC) / std::sqrt(std::fabs(corrQQ)), weight);

      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnA_pt"), track.pt(), centrality, (uy * qyA + ux * qxA) / std::sqrt(std::fabs(corrQQ)), weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnC_pt"), track.pt(), centrality, (uy * qyC + ux * qxC) / std::sqrt(std::fabs(corrQQ)), weight);
    }

    if (cfgFillMixedHarmonics) {
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v2/vnAxCxUx_eta_MH"), track.eta(), centrality, (uxMH * qxA * qxC) / corrQQx, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v2/vnAxCyUx_eta_MH"), track.eta(), centrality, (uxMH * qyA * qyC) / corrQQy, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v2/vnAxCyUy_eta_MH"), track.eta(), centrality, (uyMH * qxA * qyC) / corrQQx, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v2/vnAyCxUy_eta_MH"), track.eta(), centrality, (uyMH * qyA * qxC) / corrQQy, weight);

      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v2/vnAxCxUx_pt_MH"), track.pt(), centrality, (uxMH * qxA * qxC) / corrQQx, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v2/vnAxCyUx_pt_MH"), track.pt(), centrality, (uxMH * qyA * qyC) / corrQQy, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v2/vnAxCyUy_pt_MH"), track.pt(), centrality, (uyMH * qxA * qyC) / corrQQx, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v2/vnAyCxUy_pt_MH"), track.pt(), centrality, (uyMH * qyA * qxC) / corrQQy, weight);

      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v2/vnAxCxUx_cent_MH"), centrality, (uxMH * qxA * qxC) / corrQQx, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v2/vnAxCyUx_cent_MH"), centrality, (uxMH * qyA * qyC) / corrQQy, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v2/vnAxCyUy_cent_MH"), centrality, (uyMH * qxA * qyC) / corrQQx, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v2/vnAyCxUy_cent_MH"), centrality, (uyMH * qyA * qxC) / corrQQy, weight);

      // -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v3/vnAxCxUx_eta_MH"), track.eta(), centrality, (uxMH2 * qxA * qxC) / corrQQx, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v3/vnAxCyUx_eta_MH"), track.eta(), centrality, (uxMH2 * qyA * qyC) / corrQQy, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v3/vnAxCyUy_eta_MH"), track.eta(), centrality, (uyMH2 * qxA * qyC) / corrQQx, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v3/vnAyCxUy_eta_MH"), track.eta(), centrality, (uyMH2 * qyA * qxC) / corrQQy, weight);

      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v3/vnAxCxUx_pt_MH"), track.pt(), centrality, (uxMH2 * qxA * qxC) / corrQQx, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v3/vnAxCyUx_pt_MH"), track.pt(), centrality, (uxMH2 * qyA * qyC) / corrQQy, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v3/vnAxCyUy_pt_MH"), track.pt(), centrality, (uyMH2 * qxA * qyC) / corrQQx, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v3/vnAyCxUy_pt_MH"), track.pt(), centrality, (uyMH2 * qyA * qxC) / corrQQy, weight);

      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v3/vnAxCxUx_cent_MH"), centrality, (uxMH2 * qxA * qxC) / corrQQx, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v3/vnAxCyUx_cent_MH"), centrality, (uxMH2 * qyA * qyC) / corrQQy, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v3/vnAxCyUy_cent_MH"), centrality, (uyMH2 * qxA * qyC) / corrQQx, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("v3/vnAyCxUy_cent_MH"), centrality, (uyMH2 * qyA * qxC) / corrQQy, weight);
    }

    if (cfgFillXandYterms) {
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnAx_eta"), track.eta(), centrality, (ux * qxA) / std::sqrt(std::fabs(corrQQx)), weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnAy_eta"), track.eta(), centrality, (uy * qyA) / std::sqrt(std::fabs(corrQQy)), weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnCx_eta"), track.eta(), centrality, (ux * qxC) / std::sqrt(std::fabs(corrQQx)), weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnCy_eta"), track.eta(), centrality, (uy * qyC) / std::sqrt(std::fabs(corrQQy)), weight);

      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnAx_pt"), track.pt(), centrality, (ux * qxA) / std::sqrt(std::fabs(corrQQx)), weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnAy_pt"), track.pt(), centrality, (uy * qyA) / std::sqrt(std::fabs(corrQQy)), weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnCx_pt"), track.pt(), centrality, (ux * qxC) / std::sqrt(std::fabs(corrQQx)), weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnCy_pt"), track.pt(), centrality, (uy * qyC) / std::sqrt(std::fabs(corrQQy)), weight);
    }

    if (cfgFillEventPlane && pt == 0) { // only fill for inclusive!
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnA_eta_EP"), track.eta(), centrality, vnA, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnC_eta_EP"), track.eta(), centrality, vnC, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnFull_eta_EP"), track.eta(), centrality, vnFull, weight);

      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnA_pt_EP"), track.pt(), centrality, vnA, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnC_pt_EP"), track.pt(), centrality, vnC, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnFull_pt_EP"), track.pt(), centrality, vnFull, weight);

      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnA_cent_EP"), centrality, vnA, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnC_cent_EP"), centrality, vnC, weight);
      registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnFull_cent_EP"), centrality, vnFull, weight);
    }

    // For integrated v1 take only tracks from eta>0.
    // Following https://arxiv.org/pdf/1306.4145
    if (cfgFillGeneralV1Histos) {
      if (track.eta() < 0 && cfgHarm == 1) {
        registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnA_cent_minEta"), centrality, -1.0 * (uy * qyA + ux * qxA) / std::sqrt(std::fabs(corrQQ)), weight);
        registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnC_cent_minEta"), centrality, -1.0 * (uy * qyC + ux * qxC) / std::sqrt(std::fabs(corrQQ)), weight);

        registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnA_pt_odd"), track.pt(), centrality, -1.0 * (uy * qyA + ux * qxA) / std::sqrt(std::fabs(corrQQ)), weight);
        registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnC_pt_odd"), track.pt(), centrality, -1.0 * (uy * qyC + ux * qxC) / std::sqrt(std::fabs(corrQQ)), weight);

        if (cfgFillXandYterms) {
          registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnAx_pt_odd"), track.pt(), centrality, -1.0 * (ux * qxA) / std::sqrt(std::fabs(corrQQx)), weight);
          registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnAy_pt_odd"), track.pt(), centrality, -1.0 * (uy * qyA) / std::sqrt(std::fabs(corrQQy)), weight);
          registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnCx_pt_odd"), track.pt(), centrality, -1.0 * (ux * qxC) / std::sqrt(std::fabs(corrQQx)), weight);
          registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnCy_pt_odd"), track.pt(), centrality, -1.0 * (uy * qyC) / std::sqrt(std::fabs(corrQQy)), weight);
        }
      } else {
        registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnA_cent_plusEta"), centrality, (uy * qyA + ux * qxA) / std::sqrt(std::fabs(corrQQ)), weight);
        registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnC_cent_plusEta"), centrality, (uy * qyC + ux * qxC) / std::sqrt(std::fabs(corrQQ)), weight);

        registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnA_pt_odd"), track.pt(), centrality, (uy * qyA + ux * qxA) / std::sqrt(std::fabs(corrQQ)), weight);
        registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnC_pt_odd"), track.pt(), centrality, (uy * qyC + ux * qxC) / std::sqrt(std::fabs(corrQQ)), weight);

        if (cfgFillXandYterms) {
          registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnAx_pt_odd"), track.pt(), centrality, (ux * qxA) / std::sqrt(std::fabs(corrQQx)), weight);
          registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnAy_pt_odd"), track.pt(), centrality, (uy * qyA) / std::sqrt(std::fabs(corrQQy)), weight);
          registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnCx_pt_odd"), track.pt(), centrality, (ux * qxC) / std::sqrt(std::fabs(corrQQx)), weight);
          registry.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("vnCy_pt_odd"), track.pt(), centrality, (uy * qyC) / std::sqrt(std::fabs(corrQQy)), weight);
        }
      }
    }
  }

  template <FillType ft, ChargeType ct, ParticleType pt, typename TrackObject>
  inline void fillTrackQA(TrackObject track, double vz, float wacc = 1, float weff = 1)
  {
    if (!cfgFillTrackQA)
      return;

    static constexpr std::string_view Time[] = {"before/", "after/"};
    // NOTE: species[kUnidentified] = "" (when no PID)
    histos.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("QA/") + HIST(Time[ft]) + HIST("hPhi_Eta_Pt"), track.phi(), track.eta(), track.pt());
    histos.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("QA/") + HIST(Time[ft]) + HIST("hPhi_Eta_Pt_corrected"), track.phi(), track.eta(), track.pt(), wacc * weff);

    histos.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("QA/") + HIST(Time[ft]) + HIST("hPhi_Eta_vz"), track.phi(), track.eta(), vz);
    histos.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("QA/") + HIST(Time[ft]) + HIST("hPhi_Eta_vz_corrected"), track.phi(), track.eta(), vz, wacc * weff);
    histos.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("QA/") + HIST(Time[ft]) + HIST("hDCAxy_pt"), track.pt(), track.dcaXY(), wacc * weff);
    histos.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("QA/") + HIST(Time[ft]) + HIST("hDCAz_pt"), track.pt(), track.dcaZ(), wacc * weff);
    histos.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("QA/") + HIST(Time[ft]) + HIST("hSharedClusters_pt"), track.pt(), track.tpcFractionSharedCls(), wacc * weff);
    histos.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("QA/") + HIST(Time[ft]) + HIST("hCrossedRows_pt"), track.pt(), track.tpcNClsFound(), wacc * weff);
    histos.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("QA/") + HIST(Time[ft]) + HIST("hCrossedRows_vs_SharedClusters"), track.tpcNClsFound(), track.tpcFractionSharedCls(), wacc * weff);
  }

  template <FillType ft, ChargeType ct, ParticleType pt, typename TrackObject>
  inline void fillTrackQA(TrackObject track, double vz, double centrality, float wacc = 1, float weff = 1)
  {
    if (!cfgFillTrackQA)
      return;

    static constexpr std::string_view Time[] = {"before/", "after/"};
    // NOTE: species[kUnidentified] = "" (when no PID)
    histos.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("QA/") + HIST(Time[ft]) + HIST("hPt_Eta"), track.pt(), track.eta(), centrality, wacc * weff);
    histos.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("QA/") + HIST(Time[ft]) + HIST("hPt_Eta_uncorrected"), track.pt(), track.eta(), centrality);
    histos.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("QA/") + HIST(Time[ft]) + HIST("hPhi_Eta"), track.phi(), track.eta(), centrality, wacc * weff);
    histos.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("QA/") + HIST(Time[ft]) + HIST("hPhi_Eta_uncorrected"), track.phi(), track.eta(), centrality);

    histos.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("QA/") + HIST(Time[ft]) + HIST("hPhi_Eta_vz"), track.phi(), track.eta(), vz);
    histos.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("QA/") + HIST(Time[ft]) + HIST("hPhi_Eta_vz_corrected"), track.phi(), track.eta(), vz, wacc);
    histos.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("QA/") + HIST(Time[ft]) + HIST("hDCAxy_pt"), track.pt(), track.dcaXY(), wacc * weff);
    histos.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("QA/") + HIST(Time[ft]) + HIST("hDCAz_pt"), track.pt(), track.dcaZ(), wacc * weff);
    histos.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("QA/") + HIST(Time[ft]) + HIST("hSharedClusters_pt"), track.pt(), track.tpcFractionSharedCls(), wacc * weff);
    histos.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("QA/") + HIST(Time[ft]) + HIST("hCrossedRows_pt"), track.pt(), track.tpcNClsFound(), wacc * weff);
    histos.fill(HIST(Charge[ct]) + HIST(Species[pt]) + HIST("QA/") + HIST(Time[ft]) + HIST("hCrossedRows_vs_SharedClusters"), track.tpcNClsFound(), track.tpcFractionSharedCls(), wacc * weff);
  }

  template <FillType ft, ChargeType ct, typename TrackObject>
  inline void fillPIDQA(TrackObject track)
  {
    if (!cfgFillPIDQA || !cfgFillTrackQA)
      return;

    histos.fill(HIST(Charge[ct]) + HIST("pion/") + HIST("QA/") + HIST(Time[ft]) + HIST("hNsigmaTOF_pt"), track.pt(), track.tofNSigmaPi());
    histos.fill(HIST(Charge[ct]) + HIST("pion/") + HIST("QA/") + HIST(Time[ft]) + HIST("hNsigmaTPC_pt"), track.pt(), track.tpcNSigmaPi());
    histos.fill(HIST(Charge[ct]) + HIST("kaon/") + HIST("QA/") + HIST(Time[ft]) + HIST("hNsigmaTOF_pt"), track.pt(), track.tofNSigmaKa());
    histos.fill(HIST(Charge[ct]) + HIST("kaon/") + HIST("QA/") + HIST(Time[ft]) + HIST("hNsigmaTPC_pt"), track.pt(), track.tpcNSigmaKa());
    histos.fill(HIST(Charge[ct]) + HIST("proton/") + HIST("QA/") + HIST(Time[ft]) + HIST("hNsigmaTOF_pt"), track.pt(), track.tofNSigmaPr());
    histos.fill(HIST(Charge[ct]) + HIST("proton/") + HIST("QA/") + HIST(Time[ft]) + HIST("hNsigmaTPC_pt"), track.pt(), track.tpcNSigmaPr());

    histos.fill(HIST(Charge[ct]) + HIST("proton/") + HIST("QA/") + HIST(Time[ft]) + HIST("hPhi_Eta_Pt"), track.phi(), track.eta(), track.pt());

    histos.fill(HIST(Charge[ct]) + HIST("QA/") + HIST(Time[ft]) + HIST("hdEdxTPC_pt"), track.pt(), track.tpcSignal());
    histos.fill(HIST(Charge[ct]) + HIST("QA/") + HIST(Time[ft]) + HIST("hBetaTOF_pt"), track.pt(), track.beta());
  }

  template <FillType ft, ModeType md, typename TrackObject>
  inline void fillMCPtHistos(TrackObject track, int pdgCode, float centrality)
  {
    static constexpr std::string_view Time[] = {"before/", "after/"};
    static constexpr std::string_view Mode[] = {"Gen/", "Reco/"};

    registry.fill(HIST("trackMC") + HIST(Mode[md]) + HIST(Time[ft]) + HIST("incl/hPt_hadron"), track.pt(), track.eta(), centrality);
    if (pdgCode > 0) {
      registry.fill(HIST("trackMC") + HIST(Mode[md]) + HIST(Time[ft]) + HIST("pos/hPt_hadron"), track.pt(), track.eta(), centrality);
    } else {
      registry.fill(HIST("trackMC") + HIST(Mode[md]) + HIST(Time[ft]) + HIST("neg/hPt_hadron"), track.pt(), track.eta(), centrality);
    }

    if (pdgCode == kPiPlus || pdgCode == kPiMinus) {
      registry.fill(HIST("trackMC") + HIST(Mode[md]) + HIST(Time[ft]) + HIST("incl/hPt_pion"), track.pt(), track.eta(), centrality);
      if (pdgCode == kPiPlus) {
        registry.fill(HIST("trackMC") + HIST(Mode[md]) + HIST(Time[ft]) + HIST("pos/hPt_pion"), track.pt(), track.eta(), centrality);
      } else {
        registry.fill(HIST("trackMC") + HIST(Mode[md]) + HIST(Time[ft]) + HIST("neg/hPt_pion"), track.pt(), track.eta(), centrality);
      }
    } else if (pdgCode == kKPlus || pdgCode == kKMinus) {
      registry.fill(HIST("trackMC") + HIST(Mode[md]) + HIST(Time[ft]) + HIST("incl/hPt_kaon"), track.pt(), track.eta(), centrality);
      if (pdgCode == kKPlus) {
        registry.fill(HIST("trackMC") + HIST(Mode[md]) + HIST(Time[ft]) + HIST("pos/hPt_kaon"), track.pt(), track.eta(), centrality);
      } else {
        registry.fill(HIST("trackMC") + HIST(Mode[md]) + HIST(Time[ft]) + HIST("neg/hPt_kaon"), track.pt(), track.eta(), centrality);
      }
    } else if (pdgCode == kProton || pdgCode == kProtonBar) {
      registry.fill(HIST("trackMC") + HIST(Mode[md]) + HIST(Time[ft]) + HIST("incl/hPt_proton"), track.pt(), track.eta(), centrality);
      if (pdgCode == kProton) {
        registry.fill(HIST("trackMC") + HIST(Mode[md]) + HIST(Time[ft]) + HIST("pos/hPt_proton"), track.pt(), track.eta(), centrality);
      } else {
        registry.fill(HIST("trackMC") + HIST(Mode[md]) + HIST(Time[ft]) + HIST("neg/hPt_proton"), track.pt(), track.eta(), centrality);
      }
    }
  }

  template <FillType ft, typename McParticleObject>
  inline void fillPrimaryHistos(McParticleObject mcparticle, float centrality)
  {
    static constexpr std::string_view Time[] = {"/before", "/after"};

    if (!mcparticle.isPhysicalPrimary()) {
      registry.fill(HIST("trackMCReco") + HIST(Time[ft]) + HIST("/hIsPhysicalPrimary"), 0, centrality);
    } else {
      registry.fill(HIST("trackMCReco") + HIST(Time[ft]) + HIST("/hIsPhysicalPrimary"), 1, centrality);
    }
  }

  template <FillType ft, ParticleType ct, typename TrackObject>
  void fillAllQA(TrackObject track, double vtxz, double centrality, bool pos, float wacc = 1, float weff = 1, float waccP = 1, float weffP = 1, float waccN = 1, float weffN = 1)
  {
    if (!cfgTrackSelDoTrackQAvsCent) {
      fillTrackQA<ft, kInclusive, ct>(track, vtxz, wacc, weff);
    } else {
      fillTrackQA<ft, kInclusive, ct>(track, vtxz, centrality, wacc, weff);
    }
    fillPIDQA<ft, kInclusive>(track);
    if (pos) {
      if (!cfgTrackSelDoTrackQAvsCent) {
        fillTrackQA<ft, kPositive, ct>(track, vtxz, waccP, weffP);
      } else {
        fillTrackQA<ft, kPositive, ct>(track, vtxz, centrality, waccP, weffP);
      }
      fillPIDQA<ft, kPositive>(track);
    } else {
      if (!cfgTrackSelDoTrackQAvsCent) {
        fillTrackQA<ft, kNegative, ct>(track, vtxz, waccN, weffN);
      } else {
        fillTrackQA<ft, kNegative, ct>(track, vtxz, centrality, waccN, weffN);
      }
      fillPIDQA<ft, kNegative>(track);
    }
  }

  void processData(ZDCCollisions::iterator const& collision, aod::BCsWithTimestamps const&, UsedTracks const& tracks)
  {
    histos.fill(HIST("hEventCount"), evSel_FilteredEvent);
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int standardMagField = 99999;
    auto field = (cfgMagField == standardMagField) ? getMagneticField(bc.timestamp()) : cfgMagField;

    if (bc.runNumber() != cfg.lastRunNumber) {
      cfg.correctionsLoaded = false;
      cfg.clCentrality = false;
      cfg.lastRunNumber = bc.runNumber();
      cfg.mAcceptance.clear();
      LOGF(info, "Size of mAcceptance: %i (should be 0)", (int)cfg.mAcceptance.size());
    }

    if (cfgFillQABefore)
      fillEventQA<kBefore>(collision, tracks);

    loadCorrections(bc.timestamp());

    float centrality = collision.centFT0C();

    if (cfgCentFT0Cvariant1)
      centrality = collision.centFT0CVariant1();
    if (cfgCentFT0M)
      centrality = collision.centFT0M();
    if (cfgCentFV0A)
      centrality = collision.centFV0A();
    if (cfgCentNGlobal)
      centrality = collision.centNGlobal();

    if (!eventSelected(collision, tracks.size()))
      return;

    if (collision.isSelected()) {

      histos.fill(HIST("hEventCount"), evSel_isSelectedZDC);

      double qxA = collision.qxA();
      double qyA = collision.qyA();
      double qxC = collision.qxC();
      double qyC = collision.qyC();

      double vtxz = collision.posZ();
      double psiA = 1.0 * std::atan2(qyA, qxA);
      double psiC = 1.0 * std::atan2(qyC, qxC);

      // https://twiki.cern.ch/twiki/pub/ALICE/DirectedFlowAnalysisNote/vn_ZDC_ALICE_INT_NOTE_version02.pdf
      double psiFull = 1.0 * std::atan2(qyA + qyC, qxA + qxC);

      // always fill these histograms!
      registry.fill(HIST("QQCorrelations/qAqCXY"), centrality, qxA * qxC + qyA * qyC);
      registry.fill(HIST("QQCorrelations/qAXqCY"), centrality, qxA * qyC);
      registry.fill(HIST("QQCorrelations/qAYqCX"), centrality, qyA * qxC);
      registry.fill(HIST("QQCorrelations/qAXYqCXY"), centrality, qyA * qxC + qxA * qyC);
      registry.fill(HIST("QQCorrelations/qAqCX"), centrality, qxA * qxC);
      registry.fill(HIST("QQCorrelations/qAqCY"), centrality, qyA * qyC);

      if (cfgFillEventQA) {
        histos.fill(HIST("QA/hCentFull"), centrality, 1);
      }
      if (cfgFillEventPlaneQA) {
        histos.fill(HIST("QA/hSPplaneA"), psiA, 1);
        histos.fill(HIST("QA/hSPplaneC"), psiC, 1);
        histos.fill(HIST("QA/hSPplaneFull"), psiFull, 1);
        histos.fill(HIST("QA/hCosPhiACosPhiC"), centrality, std::cos(psiA) * std::cos(psiC));
        histos.fill(HIST("QA/hSinPhiASinPhiC"), centrality, std::sin(psiA) * std::sin(psiC));
        histos.fill(HIST("QA/hSinPhiACosPhiC"), centrality, std::sin(psiA) * std::cos(psiC));
        histos.fill(HIST("QA/hCosPhiASinsPhiC"), centrality, std::cos(psiA) * std::sin(psiC));
        histos.fill(HIST("QA/hFullEvPlaneRes"), centrality, -1 * std::cos(psiA - psiC));
      }

      if (centrality > cfgCentMax || centrality < cfgCentMin)
        return;

      histos.fill(HIST("hEventCount"), evSel_CentCuts);

      // Load correlations and SP resolution needed for Scalar Product and event plane methods.
      // Only load once!
      // If not loaded set to 1
      double corrQQ = 1., corrQQx = 1., corrQQy = 1.;
      if (cfgCCDBdir_QQ.value.empty() == false) {
        if (!cfg.clQQ) {
          TList* hcorrList = ccdb->getForTimeStamp<TList>(cfgCCDBdir_QQ.value, bc.timestamp());
          cfg.hcorrQQ = reinterpret_cast<TProfile*>(hcorrList->FindObject("qAqCXY"));
          cfg.hcorrQQx = reinterpret_cast<TProfile*>(hcorrList->FindObject("qAqCX"));
          cfg.hcorrQQy = reinterpret_cast<TProfile*>(hcorrList->FindObject("qAqCY"));
          cfg.clQQ = true;
        }
        corrQQ = cfg.hcorrQQ->GetBinContent(cfg.hcorrQQ->FindBin(centrality));
        corrQQx = cfg.hcorrQQx->GetBinContent(cfg.hcorrQQx->FindBin(centrality));
        corrQQy = cfg.hcorrQQy->GetBinContent(cfg.hcorrQQy->FindBin(centrality));
      }

      double evPlaneRes = 1.;
      if (cfgCCDBdir_SP.value.empty() == false) {
        if (!cfg.clEvPlaneRes) {
          cfg.hEvPlaneRes = ccdb->getForTimeStamp<TProfile>(cfgCCDBdir_SP.value, bc.timestamp());
          cfg.clEvPlaneRes = true;
        }
        evPlaneRes = cfg.hEvPlaneRes->GetBinContent(cfg.hEvPlaneRes->FindBin(centrality));
        if (evPlaneRes < 0)
          LOGF(fatal, "<Cos(PsiA-PsiC)> > 0 for centrality %.2f! Cannot determine resolution.. Change centrality ranges!!!", centrality);
        evPlaneRes = std::sqrt(evPlaneRes);
      }

      double centWeight = 1.0;
      if (cfgCCDBdir_centrality.value.empty() == false) {
        if (!cfg.clCentrality) {
          cfg.hCentrality = ccdb->getForTimeStamp<TH1D>(cfgCCDBdir_centrality.value, bc.timestamp());
          cfg.clCentrality = true;
        }
        centWeight = cfg.hCentrality->GetBinContent(cfg.hCentrality->FindBin(centrality));
        if (centWeight < 0)
          LOGF(fatal, "Centrality weight cannot be negative.. abort for (%.2f)", centrality);
      }

      fillEventQA<kAfter>(collision, tracks, centWeight);

      for (const auto& track : tracks) {

        int trackPID = (cfgFillPID || cfgFillPIDQA) ? getTrackPID(track) : kUnidentified;

        if (cfgFillPIDQA)
          histos.fill(HIST("hPIDcounts"), trackPID, track.pt());

        float weff = 1., wacc = 1.;
        float weffP = 1., waccP = 1.;
        float weffN = 1., waccN = 1.;

        if (track.sign() == 0.0)
          continue;

        histos.fill(HIST("hTrackCount"), trackSel_ZeroCharge);
        bool pos = (track.sign() > 0) ? true : false;

        if (cfgFillQABefore) {
          switch (trackPID) {
            case kUnidentified:
              fillAllQA<kBefore, kUnidentified>(track, vtxz, centrality, pos);
              break;
            case kPions:
              fillAllQA<kBefore, kPions>(track, vtxz, centrality, pos);
              break;
            case kKaons:
              fillAllQA<kBefore, kKaons>(track, vtxz, centrality, pos);
              break;
            case kProtons:
              fillAllQA<kBefore, kProtons>(track, vtxz, centrality, pos);
              break;
          }
        }

        if (!trackSelected(track, field))
          continue;

        // constrain angle to 0 -> [0,0+2pi]
        auto phi = RecoDecay::constrainAngle(track.phi(), 0);

        // Fill NUA weights (last 0 is for Data see GFWWeights class (not a weight))
        // ToDo: Add pi, ka, proton here!
        if (cfgFillWeights) {
          fWeights->fill(phi, track.eta(), vtxz, track.pt(), centrality, 0);
          registry.fill(HIST("weights2D/hPhi_Eta_vz"), phi, track.eta(), vtxz, 1);
        }
        if (cfgFillWeightsPOS && pos) {
          fWeightsPOS->fill(phi, track.eta(), vtxz, track.pt(), centrality, 0);
          registry.fill(HIST("weights2D/hPhi_Eta_vz_positive"), phi, track.eta(), vtxz, 1);
        }
        if (cfgFillWeightsNEG && !pos) {
          fWeightsNEG->fill(phi, track.eta(), vtxz, track.pt(), centrality, 0);
          registry.fill(HIST("weights2D/hPhi_Eta_vz_negative"), phi, track.eta(), vtxz, 1);
        }

        // Set weff and wacc for inclusive, negative and positive hadrons
        if (!setCurrentParticleWeights(kInclusive, weff, wacc, phi, track.eta(), track.pt(), vtxz))
          continue;
        if (pos && !setCurrentParticleWeights(kPositive, weffP, waccP, phi, track.eta(), track.pt(), vtxz))
          continue;
        if (!pos && !setCurrentParticleWeights(kNegative, weffN, waccN, phi, track.eta(), track.pt(), vtxz))
          continue;

        histos.fill(HIST("hTrackCount"), trackSel_ParticleWeights);

        switch (trackPID) {
          case kUnidentified:
            fillAllQA<kAfter, kUnidentified>(track, vtxz, centrality, pos, wacc, weff, waccP, weffP, waccN, weffN);
            break;
          case kPions:
            fillAllQA<kAfter, kPions>(track, vtxz, centrality, pos, wacc, weff, waccP, weffP, waccN, weffN);
            break;
          case kKaons:
            fillAllQA<kAfter, kKaons>(track, vtxz, centrality, pos, wacc, weff, waccP, weffP, waccN, weffN);
            break;
          case kProtons:
            fillAllQA<kAfter, kProtons>(track, vtxz, centrality, pos, wacc, weff, waccP, weffP, waccN, weffN);
            break;
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        auto ux = std::cos(cfgHarm * phi);
        auto uy = std::sin(cfgHarm * phi);

        auto uxMH = std::cos(cfgHarmMixed1 * phi);
        auto uyMH = std::sin(cfgHarmMixed1 * phi);

        auto uxMH2 = std::cos(cfgHarmMixed2 * phi);
        auto uyMH2 = std::sin(cfgHarmMixed2 * phi);
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        double vnA = std::cos(cfgHarm * (phi - psiA)) / evPlaneRes;
        double vnC = std::cos(cfgHarm * (phi - psiC)) / evPlaneRes;
        double vnFull = std::cos(cfgHarm * (phi - psiFull)) / evPlaneRes;

        fillHistograms<kInclusive, kUnidentified>(track, wacc, weff, centWeight, ux, uy, uxMH, uyMH, uxMH2, uyMH2, qxA, qyA, qxC, qyC, corrQQx, corrQQy, corrQQ, vnA, vnC, vnFull, centrality);

        // NOTE: these are not the right wacc and weff!! Need to do for each particle species! This is just to test for now.
        if (cfgFillPID) {
          if (trackPID == kPions) {
            fillHistograms<kInclusive, kPions>(track, wacc, weff, centWeight, ux, uy, uxMH, uyMH, uxMH2, uyMH2, qxA, qyA, qxC, qyC, corrQQx, corrQQy, corrQQ, vnA, vnC, vnFull, centrality);
            if (pos)
              fillHistograms<kPositive, kPions>(track, waccP, weffP, centWeight, ux, uy, uxMH, uyMH, uxMH2, uyMH2, qxA, qyA, qxC, qyC, corrQQx, corrQQy, corrQQ, vnA, vnC, vnFull, centrality);
            else
              fillHistograms<kNegative, kPions>(track, waccN, weffN, centWeight, ux, uy, uxMH, uyMH, uxMH2, uyMH2, qxA, qyA, qxC, qyC, corrQQx, corrQQy, corrQQ, vnA, vnC, vnFull, centrality);
          } else if (trackPID == kKaons) {
            fillHistograms<kInclusive, kKaons>(track, wacc, weff, centWeight, ux, uy, uxMH, uyMH, uxMH2, uyMH2, qxA, qyA, qxC, qyC, corrQQx, corrQQy, corrQQ, vnA, vnC, vnFull, centrality);
            if (pos)
              fillHistograms<kPositive, kKaons>(track, waccP, weffP, centWeight, ux, uy, uxMH, uyMH, uxMH2, uyMH2, qxA, qyA, qxC, qyC, corrQQx, corrQQy, corrQQ, vnA, vnC, vnFull, centrality);
            else
              fillHistograms<kNegative, kKaons>(track, waccN, weffN, centWeight, ux, uy, uxMH, uyMH, uxMH2, uyMH2, qxA, qyA, qxC, qyC, corrQQx, corrQQy, corrQQ, vnA, vnC, vnFull, centrality);
          } else if (trackPID == kProtons) {
            fillHistograms<kInclusive, kProtons>(track, wacc, weff, centWeight, ux, uy, uxMH, uyMH, uxMH2, uyMH2, qxA, qyA, qxC, qyC, corrQQx, corrQQy, corrQQ, vnA, vnC, vnFull, centrality);
            if (pos)
              fillHistograms<kPositive, kProtons>(track, waccP, weffP, centWeight, ux, uy, uxMH, uyMH, uxMH2, uyMH2, qxA, qyA, qxC, qyC, corrQQx, corrQQy, corrQQ, vnA, vnC, vnFull, centrality);
            else
              fillHistograms<kNegative, kProtons>(track, waccN, weffN, centWeight, ux, uy, uxMH, uyMH, uxMH2, uyMH2, qxA, qyA, qxC, qyC, corrQQx, corrQQy, corrQQ, vnA, vnC, vnFull, centrality);
          }
        }

        if (cfgFillChargeDependence) {
          if (pos) {
            fillHistograms<kPositive, kUnidentified>(track, waccP, weffP, centWeight, ux, uy, uxMH, uyMH, uxMH2, uyMH2, qxA, qyA, qxC, qyC, corrQQx, corrQQy, corrQQ, vnA, vnC, vnFull, centrality);
          } else {
            fillHistograms<kNegative, kUnidentified>(track, waccN, weffN, centWeight, ux, uy, uxMH, uyMH, uxMH2, uyMH2, qxA, qyA, qxC, qyC, corrQQx, corrQQy, corrQQ, vnA, vnC, vnFull, centrality);
          }
        }
      } // end of track loop
    } // end of collision isSelected loop
  }

  PROCESS_SWITCH(FlowSP, processData, "Process analysis for non-derived data", true);

  void processMCReco(CC const& collision, aod::BCsWithTimestamps const&, TCs const& tracks, FilteredTCs const& filteredTracks, aod::McParticles const&)
  {
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    int standardMagField = 99999;
    auto field = (cfgMagField == standardMagField) ? getMagneticField(bc.timestamp()) : cfgMagField;

    double vtxz = collision.posZ();
    float centrality = collision.centFT0C();
    if (cfgCentFT0Cvariant1)
      centrality = collision.centFT0CVariant1();
    if (cfgCentFT0M)
      centrality = collision.centFT0M();
    if (cfgCentFV0A)
      centrality = collision.centFV0A();
    if (cfgCentNGlobal)
      centrality = collision.centNGlobal();

    if (cfgFillQABefore)
      fillEventQA<kBefore>(collision, filteredTracks);

    if (!eventSelected(collision, filteredTracks.size()))
      return;

    if (centrality > cfgCentMax || centrality < cfgCentMin)
      return;

    histos.fill(HIST("hEventCount"), evSel_CentCuts);

    if (!collision.has_mcCollision()) {
      LOGF(info, "No mccollision found for this collision");
      return;
    }

    fillEventQA<kAfter>(collision, filteredTracks);

    registry.fill(HIST("trackMCReco/hTrackSize_unFiltered"), tracks.size(), centrality);
    registry.fill(HIST("trackMCReco/hTrackSize_Filtered"), filteredTracks.size(), centrality);

    for (const auto& track : filteredTracks) {
      auto mcParticle = track.mcParticle();
      if (track.sign() == 0.0)
        continue;
      histos.fill(HIST("hTrackCount"), trackSel_ZeroCharge);

      fillMCPtHistos<kBefore, kReco>(track, mcParticle.pdgCode(), centrality);
      bool pos = (track.sign() > 0) ? true : false;

      fillPrimaryHistos<kBefore>(mcParticle, centrality);

      // This neglects PID (for now) later use getPID like in data.
      if (cfgFillQABefore) {
        switch (std::abs(mcParticle.pdgCode())) {
          case kPions:
            fillAllQA<kBefore, kPions>(track, vtxz, centrality, pos);
            break;
          case kKaons:
            fillAllQA<kBefore, kKaons>(track, vtxz, centrality, pos);
            break;
          case kProtons:
            fillAllQA<kBefore, kProtons>(track, vtxz, centrality, pos);
            break;
          default:
            fillAllQA<kBefore, kUnidentified>(track, vtxz, centrality, pos);
        }
      }

      if (!trackSelected(track, field))
        continue;

      fillMCPtHistos<kAfter, kReco>(track, mcParticle.pdgCode(), centrality);

      fillTrackQA<kAfter, kInclusive, kUnidentified>(track, vtxz);
      fillPrimaryHistos<kAfter>(mcParticle, centrality);

    } // end of track loop
  }
  PROCESS_SWITCH(FlowSP, processMCReco, "Process analysis for MC reconstructed events", false);

  void processMCGen(aod::McCollisions const& mcCollisions, CCs const& collisions, TCs const& tracks, FilteredTCs const& filteredTracks, MCs const& McParts)
  {

    for (const auto& mcCollision : mcCollisions) {
      float centrality = -1;
      bool colSelected = true;

      // get McParticles which belong to mccollision
      auto partSlice = McParts.sliceBy(partPerMcCollision, mcCollision.globalIndex());

      // get reconstructed collision which belongs to mccollision
      auto colSlice = collisions.sliceBy(colPerMcCollision, mcCollision.globalIndex());
      registry.fill(HIST("trackMCGen/nCollReconstructedPerMcCollision"), colSlice.size());
      if (colSlice.size() != 1) { // check if MC collision is only reconstructed once! (https://indico.cern.ch/event/1425820/contributions/6170879/attachments/2947721/5180548/DDChinellato-O2AT4-HandsOn-03a.pdf)
        continue;
      }

      for (const auto& col : colSlice) {
        // get tracks that belong to reconstructed collision
        auto trackSlice = tracks.sliceBy(trackPerCollision, col.globalIndex());

        auto filteredTrackSlice = filteredTracks.sliceBy(trackPerCollision, col.globalIndex());

        centrality = col.centFT0C();
        if (cfgCentFT0Cvariant1)
          centrality = col.centFT0CVariant1();
        if (cfgCentFT0M)
          centrality = col.centFT0M();
        if (cfgCentFV0A)
          centrality = col.centFV0A();
        if (cfgCentNGlobal)
          centrality = col.centNGlobal();

        if (cfgFillQABefore)
          fillEventQA<kBefore>(col, filteredTrackSlice);

        if (trackSlice.size() < 1) {
          colSelected = false;
          continue;
        }
        if (!eventSelected(col, filteredTrackSlice.size())) {
          colSelected = false;
          continue;
        }

        if (centrality > cfgCentMax || centrality < cfgCentMin) {
          colSelected = false;
          continue;
        }
        histos.fill(HIST("hEventCount"), evSel_CentCuts);

        fillEventQA<kAfter>(col, filteredTrackSlice);

      } // leave reconstructed collision loop

      if (!colSelected)
        continue;

      float vtxz = mcCollision.posZ();

      for (const auto& particle : partSlice) {
        if (!particle.isPhysicalPrimary())
          continue;

        int charge = 0;

        auto pdgCode = particle.pdgCode();
        auto pdgInfo = pdg->GetParticle(pdgCode);
        if (pdgInfo != nullptr) {
          charge = pdgInfo->Charge();
        }

        if (std::fabs(charge) < 1)
          continue;

        bool pos = (charge > 0) ? true : false;

        fillMCPtHistos<kBefore, kGen>(particle, pdgCode, centrality);

        registry.fill(HIST("trackMCGen/before/incl/phi_eta_vtxZ_gen"), particle.phi(), particle.eta(), vtxz);

        if (pos) {
          registry.fill(HIST("trackMCGen/before/pos/phi_eta_vtxZ_gen"), particle.phi(), particle.eta(), vtxz);
        } else {
          registry.fill(HIST("trackMCGen/before/neg/phi_eta_vtxZ_gen"), particle.phi(), particle.eta(), vtxz);
        }

        if (particle.eta() < -cfgTrackSelsEta || particle.eta() > cfgTrackSelsEta || particle.pt() < cfgTrackSelsPtmin || particle.pt() > cfgTrackSelsPtmax)
          continue;

        fillMCPtHistos<kAfter, kGen>(particle, pdgCode, centrality);

        registry.fill(HIST("trackMCGen/after/incl/phi_eta_vtxZ_gen"), particle.phi(), particle.eta(), vtxz);

        if (pos) {
          registry.fill(HIST("trackMCGen/after/pos/phi_eta_vtxZ_gen"), particle.phi(), particle.eta(), vtxz);
        } else {
          registry.fill(HIST("trackMCGen/after/neg/phi_eta_vtxZ_gen"), particle.phi(), particle.eta(), vtxz);
        }
      }
    }
  }
  PROCESS_SWITCH(FlowSP, processMCGen, "Process analysis for MC generated events", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowSP>(cfgc)};
}
