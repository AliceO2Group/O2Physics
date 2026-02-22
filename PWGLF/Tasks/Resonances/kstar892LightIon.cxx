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

/// \file kstar892LightIon.cxx
/// \brief Code for K*0(892) resonance without resonance initializer in Light Ion collisions
/// \author Subhadeep Mandal <subhadeep.mandal@cern.ch>
/// \since 22/11/2025

#include "PWGLF/Utils/inelGt.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TRandom3.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>

#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using std::array;
using namespace o2::aod::rctsel;

struct Kstar892LightIon {
  SliceCache cache;

  // Histograms are defined with HistogramRegistry
  HistogramRegistry hEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hInvMass{"InvMass", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hMC{"MC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hPID{"PID", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hOthers{"Others", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  struct : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", false, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCut;
  RCTFlagsChecker rctChecker;

  struct : ConfigurableGroup {
    // Configurables for event selections
    Configurable<float> cfgVrtxZCut{"cfgVrtxZCut", 10.0f, "Accepted z-vertex range (cm)"};
    Configurable<bool> isApplysel8{"isApplysel8", true, "Apply sel8 event selection"};
    Configurable<bool> isApplyINELgt0{"isApplyINELgt0", true, "INEL>0 selection"};
    Configurable<bool> isTriggerTVX{"isTriggerTVX", true, "TriggerTVX"};
    Configurable<bool> isGoodZvtxFT0vsPV{"isGoodZvtxFT0vsPV", true, "IsGoodZvtxFT0vsPV"};
    Configurable<bool> isApplyOccCut{"isApplyOccCut", false, "Apply occupancy cut"};
    Configurable<float> cfgOccCut{"cfgOccCut", 1000., "Occupancy cut"};
    Configurable<bool> isNoSameBunchPileup{"isNoSameBunchPileup", true, "kNoSameBunchPileup"};
    Configurable<bool> isGoodITSLayersAll{"isGoodITSLayersAll", false, "Require all ITS layers to be good"};
    Configurable<bool> isNoTimeFrameBorder{"isNoTimeFrameBorder", true, "kNoTimeFrameBorder"};
    Configurable<bool> isNoITSROFrameBorder{"isNoITSROFrameBorder", true, "kNoITSROFrameBorder"};
    Configurable<bool> isApplyDeepAngle{"isApplyDeepAngle", false, "Deep Angle cut"};
    Configurable<bool> isNoCollInTimeRangeStandard{"isNoCollInTimeRangeStandard", false, "No collision in time range standard"};
    Configurable<bool> isVertexITSTPC{"isVertexITSTPC", false, "Vertex ITS TPC"};
    Configurable<bool> isVertexTOFMatched{"isVertexTOFMatched", false, "Vertex TOF Matched"};

    // check
    Configurable<bool> isApplyMCGenInelgt0{"isApplyMCGenInelgt0", true, "Apply INEL>0 cut in MC Gen Collisions"};
    Configurable<bool> isApplyMCGenTVX{"isApplyMCGenTVX", true, "Apply TVX cut in MC Gen Collisions"};
    Configurable<bool> isApplyMCGenVz{"isApplyMCGenVz", true, "Apply Vz cut in MC Gen Collisions"};

    // Configurables for track selections
    Configurable<bool> isPVContributor{"isPVContributor", true, "PV contributor track selection"}; // PV Contriuibutor
    Configurable<bool> isPrimaryTrack{"isPrimaryTrack", true, "Primary track selection"};          // kGoldenChi2 | kDCAxy | kDCAz
    Configurable<bool> isGlobalTracks{"isGlobalTracks", true, "isGlobalTracks"};

    Configurable<float> cfgCutPT{"cfgCutPT", 0.1f, "PT cut on daughter track"};
    Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta cut on daughter track"};
    Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 0.1f, "DCAxy range for tracks"};
    Configurable<float> cfgCutDCAz{"cfgCutDCAz", 0.1f, "DCAz range for tracks"};
    Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 15, "Number of mixed events per event"};
    Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
    Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
    Configurable<float> cfgRCRFC{"cfgRCRFC", 0.8f, "Crossed Rows to Findable Clusters"};
    Configurable<float> cfgITSChi2NCl{"cfgITSChi2NCl", 36.0, "ITS Chi2/NCl"};
    Configurable<float> cfgTPCChi2NClMax{"cfgTPCChi2NClMax", 4.0, "TPC Chi2/NCl"};
    Configurable<float> cfgTPCChi2NClMin{"cfgTPCChi2NClMin", 0.0, "TPC Chi2/NCl"};
    Configurable<bool> isUseITSTPCRefit{"isUseITSTPCRefit", false, "Require ITS Refit"};
    Configurable<bool> isApplyPtDepDCAxyCut{"isApplyPtDepDCAxyCut", false, "Apply pT dependent DCAxy cut"};
    Configurable<bool> isGoldenChi2{"isGoldenChi2", false, "Apply golden chi2 cut"};
    Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
    // Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", false, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
    Configurable<float> cfgBetaCutTOF{"cfgBetaCutTOF", 0.0, "cut TOF beta"};

    Configurable<bool> isAvoidsplitrackMC{"isAvoidsplitrackMC", true, "avoid split track in MC"};

    // cuts on mother
    // Configurable<bool> isApplyCutsOnMother{"isApplyCutsOnMother", false, "Enable additional cuts on Kstar mother"};
    // Configurable<float> cMaxPtMotherCut{"cMaxPtMotherCut", 15.0, "Maximum pt of mother cut"};
    // Configurable<float> cMaxMinvMotherCut{"cMaxMinvMotherCut", 1.5, "Maximum mass of mother cut"};
    Configurable<float> motherRapidityMax{"motherRapidityMax", 0.5, "Maximum rapidity of mother"};
    Configurable<float> motherRapidityMin{"motherRapidityMin", -0.5, "Minimum rapidity of mother"};

    // PID selections
    Configurable<bool> onlyTOF{"onlyTOF", false, "only TOF tracks"};
    Configurable<bool> onlyTOFHIT{"onlyTOFHIT", false, "accept only TOF hit tracks at high pt"};
    Configurable<bool> onlyTPC{"onlyTPC", false, "only TPC tracks"};
    Configurable<bool> isApplypTdepPID{"isApplypTdepPID", false, "Apply pT dependent PID"};
    Configurable<bool> isApplypTdepPIDwTOF{"isApplypTdepPIDwTOF", false, "Apply pT dependent PID with compulsory TOF condition in a pT range"};
    Configurable<bool> isApplyMID{"isApplyMID", false, "Apply particle MID"};
    Configurable<bool> isApplypTdepMID{"isApplypTdepMID", false, "Apply pT dependent particle MID"};

    Configurable<float> nsigmaCutTPCPi{"nsigmaCutTPCPi", 3.0, "TPC Nsigma cut for pions"};
    Configurable<float> nsigmaCutTPCKa{"nsigmaCutTPCKa", 3.0, "TPC Nsigma cut for kaons"};
    Configurable<float> nsigmaCutTOFPi{"nsigmaCutTOFPi", 3.0, "TOF Nsigma cut for pions"};
    Configurable<float> nsigmaCutTOFKa{"nsigmaCutTOFKa", 3.0, "TOF Nsigma cut for kaons"};
    Configurable<float> nsigmaCutCombinedKa{"nsigmaCutCombinedKa", 3.0, "Combined Nsigma cut for kaon"};
    Configurable<float> nsigmaCutCombinedPi{"nsigmaCutCombinedPi", 3.0, "Combined Nsigma cut for pion"};

    Configurable<float> nsigmaCutCombinedMID{"nsigmaCutCombinedMID", 3.0, "Combined Nsigma cut for pion in MID"};
    Configurable<float> nsigmaCutTPCMID{"nsigmaCutTPCMID", 1.0, "MID Nsigma cut for pion in TPC"};

    // Fixed variables
    float lowPtCutPID = 0.5;

    Configurable<bool> selHasFT0{"selHasFT0", true, "Has FT0?"};
    Configurable<bool> isZvtxPosSelMC{"isZvtxPosSelMC", true, "Zvtx position selection for MC events?"};
    Configurable<bool> selTVXMC{"selTVXMC", true, "apply TVX selection in MC?"};
    Configurable<bool> selINELgt0{"selINELgt0", true, "Select INEL > 0?"};
  } selectionConfig;

  Configurable<bool> calcLikeSign{"calcLikeSign", true, "Calculate Like Sign"};
  Configurable<bool> calcRotational{"calcRotational", true, "Calculate Rotational"};
  Configurable<int> cRotations{"cRotations", 3, "Number of random rotations in the rotational background"};
  Configurable<int> rotationalCut{"rotationalCut", 10, "Cut value (Rotation angle pi - pi/cut and pi + pi/cut)"};

  // Confugrable for QA histograms
  Configurable<bool> cQAplots{"cQAplots", true, "cQAplots"};
  Configurable<bool> cQAevents{"cQAevents", true, "centrality dist, DCAxy, DCAz"};

  Configurable<int> selectCentEstimator{"selectCentEstimator", 0, "Select centrality estimator: 0 - FT0M, 1 - FT0A, 2 - FT0C, 3 - FV0A"};

  Configurable<int> reflectionType{"reflectionType", 0, "Reflection: 0=Rho, 1=Omega, 2=Phi, 3=Kstar (for processRecReflection)"};

  Configurable<float> nchAcceptance{"nchAcceptance", 0.5, "Eta window to measure Nch MC for Nch vs Cent distribution"};

  Configurable<int> nBinsNch{"nBinsNch", 400, "N bins Nch (|eta|<0.8)"};
  Configurable<float> minNch{"minNch", 0, "Min Nch (|eta|<0.8)"};
  Configurable<float> maxNch{"maxNch", 400, "Max Nch (|eta|<0.8)"};

  // Configurable for histograms
  ConfigurableAxis binsCentPlot{"binsCentPlot", {110, 0.0, 110}, "Centrality axis"};
  ConfigurableAxis axisdEdx{"axisdEdx", {1, 0.0f, 200.0f}, "dE/dx (a.u.)"};
  ConfigurableAxis axisPtfordEbydx{"axisPtfordEbydx", {1, 0, 20}, "pT (GeV/c)"};
  ConfigurableAxis invMassKstarAxis{"invMassKstarAxis", {300, 0.7f, 1.3f}, "Kstar invariant mass axis"};
  ConfigurableAxis ptAxisKstar{"ptAxisKstar", {200, 0.0f, 20.0f}, "Kstar pT axis"};
  ConfigurableAxis binsImpactPar{"binsImpactPar", {100, 0, 25}, "Binning of the impact parameter axis"};
  ConfigurableAxis axisNch{"axisNch", {100, 0.0f, 100.0f}, "Number of charged particles in |y| < 0.5"};

  enum CentEstimator {
    kFT0M,
    kFT0A,
    kFT0C,
    kFV0A,
    kFV0C,
    kFV0M,
    kNEstimators // useful if you want to iterate or size things
  };

  enum PIDParticle {
    kPion,
    kKaon,
    kProton
  };

  enum PartReflection {
    kRho,
    kOmega,
    kPhi,
    kKstar
  };

  int noOfDaughters = 2;

  double pionPIDpTlow = 1, pionPIDpThigh = 2.5, kaonPIDpTlow = 0.7, kaonPIDpThigh = 2.5;

  TRandom* rn = new TRandom();

  void init(InitContext const&)
  {
    rctChecker.init(rctCut.cfgEvtRCTFlagCheckerLabel, rctCut.cfgEvtRCTFlagCheckerZDCCheck, rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);
    // Axes
    AxisSpec vertexZAxis = {60, -15., 15., "V_{Z} (cm) for plots"};
    AxisSpec ptAxis = {ptAxisKstar, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invmassAxis = {invMassKstarAxis, "Invariant mass (GeV/#it{c}^{2})"};
    AxisSpec centralityAxis = {binsCentPlot, "Centrality Axis"};
    AxisSpec impactParAxis = {binsImpactPar, "Impact Parameter (cm)"};

    // Histograms
    // Event selection
    hEventSelection.add("hVertexZ", "hVertexZ", kTH1F, {vertexZAxis});
    hEventSelection.add("hCentrality", "Centrality percentile", kTH1F, {{110, 0, 110}});
    hEventSelection.add("hOccupancy", "Occupancy distribution", kTH1F, {{1000, 0, 15000}});

    hEventSelection.add("hEventCut", "No. of event after cuts", kTH1D, {{20, 0, 20}});
    auto check = [](bool enabled) { return enabled ? "" : " #otimes"; }; // check if a cut is enabled and put #otimes beside that label if not enabled
    std::vector<std::string> eveCutLabels = {
      "All Events",
      Form("|Vz| < %.1f", selectionConfig.cfgVrtxZCut.value),
      std::string("sel8") + check(selectionConfig.isApplysel8.value),
      std::string("kNoTimeFrameBorder") + check(selectionConfig.isNoTimeFrameBorder.value),
      std::string("kNoITSROFrameBorder") + check(selectionConfig.isNoITSROFrameBorder.value),
      std::string("kIsTriggerTVX") + check(selectionConfig.isTriggerTVX.value),
      std::string("kNoSameBunchPileup") + check(selectionConfig.isNoSameBunchPileup.value),
      std::string("kIsGoodITSLayersAll") + check(selectionConfig.isGoodITSLayersAll.value),
      std::string("kNoCollInTimeRangeStandard") + check(selectionConfig.isNoCollInTimeRangeStandard.value),
      Form("Occupancy < %.0f%s", selectionConfig.cfgOccCut.value, check(selectionConfig.isApplyOccCut.value)),
      std::string("rctChecker") + check(rctCut.requireRCTFlagChecker.value),
      std::string("kIsGoodZvtxFT0vsPV") + check(selectionConfig.isGoodZvtxFT0vsPV.value),
      std::string("isVertexITSTPC") + check(selectionConfig.isVertexITSTPC.value),
      std::string("isVertexTOFMatched") + check(selectionConfig.isVertexTOFMatched.value),
      std::string("INEL > 0") + check(selectionConfig.isApplyINELgt0.value)};
    // assign labels
    for (size_t i = 0; i < eveCutLabels.size(); ++i) {
      hEventSelection.get<TH1>(HIST("hEventCut"))->GetXaxis()->SetBinLabel(i + 1, eveCutLabels[i].c_str());
    }

    // for primary tracksbinsCentPlot
    if (cQAplots) {
      hOthers.add("dE_by_dx_TPC", "dE/dx signal in the TPC as a function of pT", kTH2F, {axisPtfordEbydx, axisdEdx});
      hOthers.add("hEta_after", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
      hOthers.add("hCRFC_after", "CRFC after distribution", kTH1F, {{100, 0.0f, 10.0f}});
      hOthers.add("hCRFC_before", "CRFC before distribution", kTH1F, {{100, 0.0f, 10.0f}});

      hOthers.add("hKstar_rap_pt", "Pair rapidity distribution; y; p_{T}; Counts", kTH2F, {{400, -2.0f, 2.0f}, ptAxis});
      hOthers.add("hKstar_eta_pt", "Pair eta distribution; #eta; p_{T}; Counts", kTH2F, {{400, -2.0f, 2.0f}, ptAxis});

      hOthers.add("hDcaxyPi", "Dcaxy distribution of selected Pions", kTH1F, {{200, -1.0f, 1.0f}});
      hOthers.add("hDcaxyKa", "Dcaxy distribution of selected Kaons", kTH1F, {{200, -1.0f, 1.0f}});
      hOthers.add("hDcazPi", "Dcaz distribution of selected Pions", kTH1F, {{200, -1.0f, 1.0f}});
      hOthers.add("hDcazKa", "Dcaz distribution of selected Kaons", kTH1F, {{200, -1.0f, 1.0f}});

      hOthers.add("hDcaxy_cent_pt", "Dcaxy distribution before PID", kTH3F, {{200, -1.0f, 1.0f}, centralityAxis, ptAxis});
      hOthers.add("hDcaz_cent_pt", "Dcaz distribution before PID", kTH3F, {{200, -1.0f, 1.0f}, centralityAxis, ptAxis});

      hPID.add("Before/hNsigma_TPC_TOF_Ka_pt", "N #sigma Kaon TPC TOF before", kTH3F, {{50, -5.0f, 5.0f}, {50, -5.0f, 5.0f}, ptAxis});
      hPID.add("Before/hNsigma_TPC_TOF_Pi_pt", "N #sigma Pion TPC TOF before", kTH3F, {{50, -5.0f, 5.0f}, {50, -5.0f, 5.0f}, ptAxis});

      hPID.add("Before/hTPCnsigKa_Neg_mult_pt", "TPC nsigma of K^{-} before PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("Before/hTPCnsigPi_Neg_mult_pt", "TPC nsigma of #pi^{-} before PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("Before/hTOFnsigKa_Neg_mult_pt", "TOF nsigma of K^{-} before PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("Before/hTOFnsigPi_Neg_mult_pt", "TOF nsigma of #pi^{-} before PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});

      hPID.add("Before/hTPCnsigKa_Pos_mult_pt", "TPC nsigma of K^{+} before PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("Before/hTPCnsigPi_Pos_mult_pt", "TPC nsigma of #pi^{+} before PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("Before/hTOFnsigKa_Pos_mult_pt", "TOF nsigma of K^{+} before PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("Before/hTOFnsigPi_Pos_mult_pt", "TOF nsigma of #pi^{+} before PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});

      hPID.add("After/hNsigma_TPC_TOF_Ka_pt", "N #sigma Kaon TPC TOF after", kTH3F, {{50, -5.0f, 5.0f}, {50, -5.0f, 5.0f}, ptAxis});
      hPID.add("After/hNsigma_TPC_TOF_Pi_pt", "N #sigma Pion TPC TOF after", kTH3F, {{50, -5.0f, 5.0f}, {50, -5.0f, 5.0f}, ptAxis});

      hPID.add("After/hTPCnsigKa_Neg_mult_pt", "TPC nsigma of K^{-} after PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("After/hTPCnsigPi_Neg_mult_pt", "TPC nsigma of #pi^{-} after PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("After/hTOFnsigKa_Neg_mult_pt", "TOF nsigma of K^{-} after PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("After/hTOFnsigPi_Neg_mult_pt", "TOF nsigma of #pi^{-} after PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});

      hPID.add("After/hTPCnsigKa_Pos_mult_pt", "TPC nsigma of K^{+} after PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("After/hTPCnsigPi_Pos_mult_pt", "TPC nsigma of #pi^{+} after PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("After/hTOFnsigKa_Pos_mult_pt", "TOF nsigma of K^{+} after PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("After/hTOFnsigPi_Pos_mult_pt", "TOF nsigma of #pi^{+} after PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
    }

    // KStar histograms
    hInvMass.add("h3KstarInvMassUnlikeSign", "kstar Unlike Sign", kTH3F, {centralityAxis, ptAxis, invmassAxis});
    hInvMass.add("h3KstarInvMassMixed", "kstar Mixed", kTH3F, {centralityAxis, ptAxis, invmassAxis});
    if (calcLikeSign) {
      hInvMass.add("h3KstarInvMasslikeSignPP", "kstar like Sign", kTH3F, {centralityAxis, ptAxis, invmassAxis});
      hInvMass.add("h3KstarInvMasslikeSignMM", "kstar like Sign", kTH3F, {centralityAxis, ptAxis, invmassAxis});
    }
    if (calcRotational)
      hInvMass.add("h3KstarInvMassRotated", "kstar rotated", kTH3F, {centralityAxis, ptAxis, invmassAxis});

    // MC histograms
    if (doprocessGen) {
      hMC.add("Gen/hGenNo", "MC Event statistics", kTH1F, {{10, 0.0f, 10.0f}});
      hMC.add("Gen/hk892GenpT", "pT distribution of True MC K(892)0", kTH2F, {ptAxis, centralityAxis});
      hMC.add("Gen/hk892GenpT2", "pT distribution of True MC K(892)0", kTH2F, {ptAxis, centralityAxis});
      hMC.add("Gen/h1genmass", "Invariant mass of generated kstar meson", kTH1F, {invmassAxis});
      hMC.add("Gen/h1GenCent", "centrality generated", kTH1F, {centralityAxis});
      hMC.add("Gen/hAllGenCollisions", "All generated events", kTH1F, {centralityAxis});
      hMC.add("Gen/hAllGenCollisions1Rec", "All gen events with at least one rec event", kTH1F, {centralityAxis});
      hMC.add("Gen/hAllKstarGenCollisisons", "All generated Kstar in events with rapidity in 0.5", kTH2F, {ptAxis, centralityAxis});
      hMC.add("Gen/hAllKstarGenCollisisons1Rec", "All generated Kstar in events with at least one rec event in rapidity in 0.5", kTH2F, {ptAxis, centralityAxis});
    }

    if (doprocessRec) {
      hMC.add("Rec/hAllRecCollisions", "All reconstructed events", kTH1F, {centralityAxis});
      hMC.add("Rec/h1KstarRecMass", "Invariant mass of kstar meson", kTH1F, {invmassAxis});
      hMC.add("Rec/h2KstarRecpt1", "pT of kstar meson", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("Rec/h2KstarRecpt2", "pT of kstar meson", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("Rec/h1RecCent", "centrality reconstructed", kTH1F, {centralityAxis});
      hMC.add("Rec/h1KSRecsplit", "KS meson Rec split", kTH1F, {{100, 0.0f, 10.0f}});
    }

    // Signal Loss & Event Loss
    if (doprocessEvtLossSigLossMC) {
      hMC.add("ImpactCorr/hImpactParameterGen", "Impact parameter of generated MC events", kTH1F, {impactParAxis});
      hMC.add("ImpactCorr/hImpactParameterRec", "Impact parameter of selected MC events", kTH1F, {impactParAxis});
      hMC.add("ImpactCorr/hImpactParvsCentrRec", "Impact parameter of selected MC events vs centrality", kTH2F, {{centralityAxis}, impactParAxis});
      hMC.add("ImpactCorr/hKstarGenBeforeEvtSel", "K*0 before event selections", kTH2F, {ptAxis, impactParAxis});
      hMC.add("ImpactCorr/hKstarGenAfterEvtSel", "K*0 after event selections", kTH2F, {ptAxis, impactParAxis});
    }

    if (doprocessLossMCMultiplicity) {
      hMC.add("LossMult/hMultMC", "Charged Paticle multiplicity in generated MC before event selection", kTH1F, {axisNch});
      hMC.add("LossMult/hCentVsMultMC", "Centrality vs Charged Particle Multiplicity", kTH2F, {centralityAxis, axisNch});
      hMC.add("LossMult/hCentVsMultMC_EvtSel", "Centrality vs Charged Particle Multiplicity after event selection", kTH2F, {centralityAxis, axisNch});
      hMC.add("LossMult/hGenEvt_vs_multMC", "Charged Paticle multiplicity in generated MC after event selection", kTH1F, {axisNch});
      hMC.add("LossMult/hGenEvtRecoEvt_vs_multMC", "Charged Paticle multiplicity in generated MC before event selection with reconstruction", kTH1F, {axisNch});
      hMC.add("LossMult/hGenKstar_vs_pt_vs_multMC", "pT vs Charged particle multiplicity", kTH2F, {ptAxis, axisNch});
      hMC.add("LossMult/hGenKstarRecoEvt_vs_pt_vs_multMC", "pT vs Charged particle multiplicity with reconstruction", kTH2F, {ptAxis, axisNch});
    }

    if (doprocessAllLossMC) {
      hMC.add("AllLoss/hImpactParameterGen", "Impact parameter of generated MC events", kTH1F, {impactParAxis});
      hMC.add("AllLoss/hMultEta05Gen", "Charged Paticle multiplicity in generated MC before event selection in #eta 0.5", kTH1F, {axisNch});
      hMC.add("AllLoss/hMultEta08Gen", "Charged Paticle multiplicity in generated MC before event selection in #eta 0.8", kTH1F, {axisNch});
      hMC.add("AllLoss/hImpactParvsCentr", "Impact parameter of all MC events vs centrality", kTH2F, {{centralityAxis}, impactParAxis});
      hMC.add("AllLoss/hMultEta05vsCentr", "Centrality vs Charged Particle Multiplicity in #eta 0.5 for all MC", kTH2F, {centralityAxis, axisNch});
      hMC.add("AllLoss/hMultEta08vsCentr", "Centrality vs Charged Particle Multiplicity in #eta 0.8 for all MC", kTH2F, {centralityAxis, axisNch});
      hMC.add("AllLoss/hImpactParameterRec", "Impact parameter of reconstructed MC events", kTH1F, {impactParAxis});
      hMC.add("AllLoss/hImpactParvsCentrRec", "Impact parameter of selected MC events vs centrality", kTH2F, {{centralityAxis}, impactParAxis});
      hMC.add("AllLoss/hMultEta05Rec", "Charged Paticle multiplicity in reconstructed MC before event selection in #eta 0.5", kTH1F, {axisNch});
      hMC.add("AllLoss/hMultEta05vsCentrRec", "Centrality vs Charged Particle Multiplicity in #eta 0.5 after event selection", kTH2F, {centralityAxis, axisNch});
      hMC.add("AllLoss/hMultEta08Rec", "Charged Paticle multiplicity in reconstructed MC before event selection in #eta 0.8", kTH1F, {axisNch});
      hMC.add("AllLoss/hMultEta08vsCentrRec", "Centrality vs Charged Particle Multiplicity in #eta 0.8 after event selection", kTH2F, {centralityAxis, axisNch});
      hMC.add("AllLoss/hKstarpTGenVsImpactParBeforeEvtSel", "K*0 before event selections", kTH2F, {ptAxis, impactParAxis});
      hMC.add("AllLoss/hKstarpTGenVsMultEta05BeforeEvtSel", "pT vs Charged particle multiplicity in #eta 0.5 before event selection", kTH2F, {ptAxis, axisNch});
      hMC.add("AllLoss/hKstarpTGenVsMultEta08BeforeEvtSel", "pT vs Charged particle multiplicity in #eta 0.8 before event selection", kTH2F, {ptAxis, axisNch});
      hMC.add("AllLoss/hKstarpTGenVsImpactParAfterEvtSel", "K*0 after event selections", kTH2F, {ptAxis, impactParAxis});
      hMC.add("AllLoss/hKstarpTGenVsMultEta05AfterEvtSel", "pT vs Charged particle multiplicity in #eta 0.5 after event selection", kTH2F, {ptAxis, axisNch});
      hMC.add("AllLoss/hKstarpTGenVsMultEta08AfterEvtSel", "pT vs Charged particle multiplicity in #eta 0.8 after event selection", kTH2F, {ptAxis, axisNch});
    }

    if (doprocessRecMisID) {
      hMC.add("RecMisID/hMassMisID", "Reconstruction misidentification", kTH3F, {ptAxis, centralityAxis, invmassAxis});
    }

    if (doprocessRecReflection) {
      hMC.add("Reflections/hReflection", "Refelction template of Rho", kTH3F, {ptAxis, centralityAxis, invmassAxis});
    }

    if (doprocessMCCheck) {
      hMC.add("MCCheck/CentVsFoundFT0", "Found(=1.5) NOT Found(=0.5);;Status;", kTH2F, {{{centralityAxis}, {2, 0, 2}}});
      hMC.add("MCCheck/zPosMC", "Generated Events With at least One Rec. Collision + Sel. criteria;;Entries;", kTH1F, {vertexZAxis});
      hMC.add("MCCheck/zPos", "With Event Selection;;Entries;", kTH1F, {vertexZAxis});
      hMC.add("MCCheck/Cent", ";;Entries", kTH1F, {centralityAxis});
      hMC.add("MCCheck/NchVsCent", "Measured Nch v.s. Centrality (At least Once Rec. Coll. + Sel. criteria);;Nch", kTH2F, {{centralityAxis, {nBinsNch, minNch, maxNch}}});

      // MC events passing the TVX requirement
      hMC.add("MCCheck/NchMCcentVsTVX", ";Passed(=1.5) NOT Passed(=0.5);", kTH2F, {{{nBinsNch, minNch, maxNch}, {2, 0, 2}}});

      hMC.add("MCCheck/NumberOfRecoCollisions", "Number of times Gen. Coll.are reconstructed;N;Entries", kTH1F, {{10, -0.5, 9.5}});

      // Needed for the Gen. Nch to Centrality conversion
      hMC.add("MCCheck/NchMCVsCent", "Generated Nch v.s. Centrality (At least Once Rec. Coll. + Sel. criteria);;Gen. Nch MC (|#eta|<0.8)", kTH2F, {{centralityAxis, {nBinsNch, minNch, maxNch}}});

      // Needed to measure Event Loss
      hMC.add("MCCheck/NchMC_WithRecoEvt", "Generated Nch of Evts With at least one Rec. Coll. + Sel. criteria;Gen. Nch MC (|#eta|<0.8);Entries", kTH1F, {{nBinsNch, minNch, maxNch}});
      hMC.add("MCCheck/NchMC_AllGen", "Generated Nch of All Gen. Evts.;Gen. Nch;Entries", kTH1F, {{nBinsNch, minNch, maxNch}});

      // Needed to measure Event Splitting
      hMC.add("MCCheck/Centrality_WRecoEvt", "Generated Events With at least One Rec. Collision And NO Sel. criteria;;Entries", kTH1F, {centralityAxis});
      hMC.add("MCCheck/Centrality_WRecoEvtWSelCri", "Generated Events With at least One Rec. Collision + Sel. criteria;;Entries", kTH1F, {centralityAxis});
      hMC.add("MCCheck/Centrality_AllRecoEvt", "Generated Events Irrespective of the number of times it was reconstructed + Evt. Selections;;Entries", kTH1F, {centralityAxis});

      hMC.add("MCCheck/PtKstarVsCentMC_WithRecoEvt", "Generated Events With at least One Rec. Collision;;;", kTH2F, {ptAxis, centralityAxis});

      // Needed to calculate the numerator of the Signal Loss correction
      hMC.add("MCCheck/PtKstarVsNchMC_WithRecoEvt", "Generated Events With at least One Rec. Collision;;Gen. Nch (|#eta|<0.8);", kTH2F, {{ptAxis, {nBinsNch, minNch, maxNch}}});

      // Needed to calculate the denominator of the Signal Loss correction
      hMC.add("MCCheck/PtKstarVsNchMC_AllGen", "All Generated Events;;Gen. Nch (|#eta|<0.8);", kTH2F, {{ptAxis, {nBinsNch, minNch, maxNch}}});
    }
  }

  double massPi = o2::constants::physics::MassPiPlus;
  double massKa = o2::constants::physics::MassKPlus;

  template <typename Coll>
  bool selectionEvent(const Coll& collision, bool fillHist = false) // default to false
  {
    if (fillHist)
      hEventSelection.fill(HIST("hEventCut"), 0);

    if (std::abs(collision.posZ()) > selectionConfig.cfgVrtxZCut)
      return false;
    if (fillHist)
      hEventSelection.fill(HIST("hEventCut"), 1);

    if (selectionConfig.isApplysel8 && !collision.sel8())
      return false;
    if (fillHist)
      hEventSelection.fill(HIST("hEventCut"), 2);

    if (selectionConfig.isNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return false;
    if (fillHist)
      hEventSelection.fill(HIST("hEventCut"), 3);

    if (selectionConfig.isNoITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))
      return false;
    if (fillHist)
      hEventSelection.fill(HIST("hEventCut"), 4);

    if (selectionConfig.isTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX))
      return false;
    if (fillHist)
      hEventSelection.fill(HIST("hEventCut"), 5);

    if (selectionConfig.isNoSameBunchPileup && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)))
      return false;
    if (fillHist)
      hEventSelection.fill(HIST("hEventCut"), 6);

    if (selectionConfig.isGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll))
      return false;
    if (fillHist)
      hEventSelection.fill(HIST("hEventCut"), 7);

    if (selectionConfig.isNoCollInTimeRangeStandard && (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)))
      return false;
    if (fillHist)
      hEventSelection.fill(HIST("hEventCut"), 8);

    if (selectionConfig.isApplyOccCut && (std::abs(collision.trackOccupancyInTimeRange()) > selectionConfig.cfgOccCut))
      return false;
    if (fillHist)
      hEventSelection.fill(HIST("hEventCut"), 9);

    if (rctCut.requireRCTFlagChecker && !rctChecker(collision))
      return false;
    if (fillHist)
      hEventSelection.fill(HIST("hEventCut"), 10);

    if (selectionConfig.isGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
      return false;
    if (fillHist)
      hEventSelection.fill(HIST("hEventCut"), 11);

    if (selectionConfig.isVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (fillHist)
      hEventSelection.fill(HIST("hEventCut"), 12);

    if (selectionConfig.isVertexTOFMatched && !collision.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if (fillHist)
      hEventSelection.fill(HIST("hEventCut"), 13);

    if (selectionConfig.isApplyINELgt0 && !collision.isInelGt0()) {
      return false;
    }
    if (fillHist)
      hEventSelection.fill(HIST("hEventCut"), 14);

    return true;
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (selectionConfig.isGlobalTracks) {
      if (!candidate.isGlobalTrackWoDCA())
        return false;
      if (std::abs(candidate.pt()) < selectionConfig.cfgCutPT)
        return false;

      if (std::abs(candidate.eta()) > selectionConfig.cfgCutEta)
        return false;
      if (!selectionConfig.isApplyPtDepDCAxyCut) {
        if (std::abs(candidate.dcaXY()) > selectionConfig.cfgCutDCAxy)
          return false;
      } else {
        if (std::abs(candidate.dcaXY()) > (0.0105 + 0.035 / std::pow(candidate.pt(), 1.1)))
          return false;
      }
      if (selectionConfig.isGoldenChi2 && !candidate.passedGoldenChi2())
        return false;
      if (std::abs(candidate.dcaZ()) > selectionConfig.cfgCutDCAz)
        return false;
      if (candidate.tpcCrossedRowsOverFindableCls() < selectionConfig.cfgRCRFC)
        return false;
      if (candidate.itsNCls() < selectionConfig.cfgITScluster)
        return false;
      if (candidate.tpcNClsFound() < selectionConfig.cfgTPCcluster)
        return false;
      if (candidate.itsChi2NCl() >= selectionConfig.cfgITSChi2NCl)
        return false;
      if (candidate.tpcChi2NCl() >= selectionConfig.cfgTPCChi2NClMax || candidate.tpcChi2NCl() < selectionConfig.cfgTPCChi2NClMin)
        return false;
      if (selectionConfig.isPVContributor && !candidate.isPVContributor())
        return false;
      if (selectionConfig.isUseITSTPCRefit && (!(o2::aod::track::ITSrefit) || !(o2::aod::track::TPCrefit)))
        return false;
    } else if (!selectionConfig.isGlobalTracks) {
      if (std::abs(candidate.pt()) < selectionConfig.cfgCutPT)
        return false;
      // if (std::abs(candidate.eta()) > selectionConfig.cfgCutEta || std::abs(candidate.eta()) < selectionConfig.cfgCutEtaMin)
      if (std::abs(candidate.eta()) > selectionConfig.cfgCutEta)
        return false;
      // if (std::abs(candidate.dcaXY()) > selectionConfig.cfgCutDCAxy || std::abs(candidate.dcaXY()) < selectionConfig.cfgCutDCAxyMin)
      if (std::abs(candidate.dcaXY()) > selectionConfig.cfgCutDCAxy)
        return false;
      if (std::abs(candidate.dcaZ()) > selectionConfig.cfgCutDCAz)
        return false;
      if (candidate.tpcCrossedRowsOverFindableCls() < selectionConfig.cfgRCRFC)
        return false;
      if (candidate.itsNCls() < selectionConfig.cfgITScluster)
        return false;
      if (candidate.tpcNClsFound() < selectionConfig.cfgTPCcluster)
        return false;
      if (candidate.itsChi2NCl() >= selectionConfig.cfgITSChi2NCl)
        return false;
      if (candidate.tpcChi2NCl() >= selectionConfig.cfgTPCChi2NClMax || candidate.tpcChi2NCl() < selectionConfig.cfgTPCChi2NClMin)
        return false;
      if (selectionConfig.isPVContributor && !candidate.isPVContributor())
        return false;
      if (selectionConfig.isPrimaryTrack && !candidate.isPrimaryTrack())
        return false;
    }

    return true;
  }

  // deep angle cut on pair to remove photon conversion
  template <typename T1, typename T2>
  bool selectionPair(const T1& candidate1, const T2& candidate2)
  {
    double pt1, pt2, pz1, pz2, p1, p2, angle;
    pt1 = candidate1.pt();
    pt2 = candidate2.pt();
    pz1 = candidate1.pz();
    pz2 = candidate2.pz();
    p1 = candidate1.p();
    p2 = candidate2.p();
    angle = std::acos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    if (selectionConfig.isApplyDeepAngle && angle < selectionConfig.cfgDeepAngle) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool selectionPID(const T& candidate, int PID)
  {
    if (PID == PIDParticle::kPion) {
      if (selectionConfig.onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < selectionConfig.nsigmaCutTOFPi && candidate.beta() > selectionConfig.cfgBetaCutTOF) {
          return true;
        }
      } else if (selectionConfig.onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < selectionConfig.nsigmaCutTOFPi && candidate.beta() > selectionConfig.cfgBetaCutTOF) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < selectionConfig.nsigmaCutTPCPi) {
          return true;
        }
      } else if (selectionConfig.onlyTPC) {
        if (std::abs(candidate.tpcNSigmaPi()) < selectionConfig.nsigmaCutTPCPi) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < (selectionConfig.nsigmaCutCombinedPi * selectionConfig.nsigmaCutCombinedPi) && candidate.beta() > selectionConfig.cfgBetaCutTOF) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < selectionConfig.nsigmaCutTPCPi) {
          return true;
        }
      }
    } else if (PID == PIDParticle::kKaon) {
      if (selectionConfig.onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < selectionConfig.nsigmaCutTOFKa && candidate.beta() > selectionConfig.cfgBetaCutTOF) {
          return true;
        }
      } else if (selectionConfig.onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < selectionConfig.nsigmaCutTOFKa && candidate.beta() > selectionConfig.cfgBetaCutTOF) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < selectionConfig.nsigmaCutTPCKa) {
          return true;
        }
      } else if (selectionConfig.onlyTPC) {
        if (std::abs(candidate.tpcNSigmaKa()) < selectionConfig.nsigmaCutTPCKa) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (selectionConfig.nsigmaCutCombinedKa * selectionConfig.nsigmaCutCombinedKa) && candidate.beta() > selectionConfig.cfgBetaCutTOF) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < selectionConfig.nsigmaCutTPCKa) {
          return true;
        }
      }
    }
    return false;
  }

  template <typename T>
  bool selectionPIDpTdep(const T& candidate, int PID)
  {
    if (PID == PIDParticle::kPion) {
      if (candidate.pt() < selectionConfig.lowPtCutPID && std::abs(candidate.tpcNSigmaPi()) < selectionConfig.nsigmaCutTPCPi) {
        return true;
      }
      if (candidate.pt() >= selectionConfig.lowPtCutPID && std::abs(candidate.tpcNSigmaPi()) < selectionConfig.nsigmaCutTPCPi && candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < selectionConfig.nsigmaCutTOFPi) {
        return true;
      }
      if (candidate.pt() >= selectionConfig.lowPtCutPID && std::abs(candidate.tpcNSigmaPi()) < selectionConfig.nsigmaCutTPCPi && !candidate.hasTOF()) {
        return true;
      }
    } else if (PID == PIDParticle::kKaon) {
      if (candidate.pt() < selectionConfig.lowPtCutPID && std::abs(candidate.tpcNSigmaKa()) < selectionConfig.nsigmaCutTPCKa) {
        return true;
      }
      if (candidate.pt() >= selectionConfig.lowPtCutPID && std::abs(candidate.tpcNSigmaKa()) < selectionConfig.nsigmaCutTPCKa && candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < selectionConfig.nsigmaCutTOFKa) {
        return true;
      }
      if (candidate.pt() >= selectionConfig.lowPtCutPID && std::abs(candidate.tpcNSigmaKa()) < selectionConfig.nsigmaCutTPCKa && !candidate.hasTOF()) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selectionPIDpTdepTOF(const T& candidate, int PID)
  {
    if (PID == PIDParticle::kPion) {
      if (candidate.pt() < pionPIDpTlow || candidate.pt() > pionPIDpThigh) {
        if (candidate.pt() < selectionConfig.lowPtCutPID && std::abs(candidate.tpcNSigmaPi()) < selectionConfig.nsigmaCutTPCPi) {
          return true;
        }
        if (candidate.pt() >= selectionConfig.lowPtCutPID && std::abs(candidate.tpcNSigmaPi()) < selectionConfig.nsigmaCutTPCPi && candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < selectionConfig.nsigmaCutTOFPi) {
          return true;
        }
        if (candidate.pt() >= selectionConfig.lowPtCutPID && std::abs(candidate.tpcNSigmaPi()) < selectionConfig.nsigmaCutTPCPi && !candidate.hasTOF()) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < (selectionConfig.nsigmaCutCombinedPi * selectionConfig.nsigmaCutCombinedPi)) {
          return true;
        }
      }
    } else if (PID == PIDParticle::kKaon) {
      if (candidate.pt() < kaonPIDpTlow || candidate.pt() > kaonPIDpThigh) {
        if (candidate.pt() < selectionConfig.lowPtCutPID && std::abs(candidate.tpcNSigmaKa()) < selectionConfig.nsigmaCutTPCKa) {
          return true;
        }
        if (candidate.pt() >= selectionConfig.lowPtCutPID && std::abs(candidate.tpcNSigmaKa()) < selectionConfig.nsigmaCutTPCKa && candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < selectionConfig.nsigmaCutTOFKa) {
          return true;
        }
        if (candidate.pt() >= selectionConfig.lowPtCutPID && std::abs(candidate.tpcNSigmaKa()) < selectionConfig.nsigmaCutTPCKa && !candidate.hasTOF()) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (selectionConfig.nsigmaCutCombinedKa * selectionConfig.nsigmaCutCombinedKa)) {
          return true;
        }
      }
    }
    return false;
  }

  template <typename T>
  bool selectionMID(const T& candidate, int PID)
  {
    if (PID == PIDParticle::kPion) {
      if (selectionConfig.onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < selectionConfig.nsigmaCutTPCMID) {
          return true;
        }
      } else if (selectionConfig.onlyTPC) {
        if (std::abs(candidate.tpcNSigmaPi()) < selectionConfig.nsigmaCutTPCMID) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < (selectionConfig.nsigmaCutCombinedMID * selectionConfig.nsigmaCutCombinedMID)) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < selectionConfig.nsigmaCutTPCMID) {
          return true;
        }
      }
    } else if (PID == PIDParticle::kKaon) {
      if (selectionConfig.onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < selectionConfig.nsigmaCutTPCMID) {
          return true;
        }
      } else if (selectionConfig.onlyTPC) {
        if (std::abs(candidate.tpcNSigmaKa()) < selectionConfig.nsigmaCutTPCMID) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (selectionConfig.nsigmaCutCombinedMID * selectionConfig.nsigmaCutCombinedMID)) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < selectionConfig.nsigmaCutTPCMID) {
          return true;
        }
      }
    }
    return false;
  }

  template <typename T>
  bool selectionMIDpTdep(const T& candidate, int PID)
  {
    if (PID == PIDParticle::kPion) {
      if (candidate.pt() >= pionPIDpTlow && candidate.pt() < pionPIDpThigh && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < selectionConfig.nsigmaCutTPCMID) {
        return true;
      }
    } else if (PID == PIDParticle::kKaon) {
      if (candidate.pt() >= kaonPIDpTlow && candidate.pt() < kaonPIDpThigh && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < selectionConfig.nsigmaCutTPCMID) {
        return true;
      }
    }
    return false;
  }

  //*********Varibles declaration***************
  float centrality{-1.0}, theta2;
  ROOT::Math::PxPyPzMVector daughter1, daughter2, daughterRot, mother, motherRot;
  bool isMix = false;

  template <typename T1, typename T2>
  void fillInvMass(const T1& daughter1, const T1& daughter2, const T1& mother, float centrality, bool isMix, const T2& track1, const T2& track2)
  {
    if (track1.sign() * track2.sign() < 0) {
      if (!isMix) {
        if (mother.Rapidity() > selectionConfig.motherRapidityMin && mother.Rapidity() < selectionConfig.motherRapidityMax) {
          hInvMass.fill(HIST("h3KstarInvMassUnlikeSign"), centrality, mother.Pt(), mother.M());
        }
        for (int i = 0; i < cRotations; i++) {
          theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / rotationalCut, o2::constants::math::PI + o2::constants::math::PI / rotationalCut);

          daughterRot = ROOT::Math::PxPyPzMVector(daughter1.Px() * std::cos(theta2) - daughter1.Py() * std::sin(theta2), daughter1.Px() * std::sin(theta2) + daughter1.Py() * std::cos(theta2), daughter1.Pz(), daughter1.M());
          motherRot = daughterRot + daughter2;

          if (calcRotational && (motherRot.Rapidity() > selectionConfig.motherRapidityMin && motherRot.Rapidity() < selectionConfig.motherRapidityMax))
            hInvMass.fill(HIST("h3KstarInvMassRotated"), centrality, motherRot.Pt(), motherRot.M());
        }
      } else if (isMix && (mother.Rapidity() > selectionConfig.motherRapidityMin && mother.Rapidity() < selectionConfig.motherRapidityMax)) {
        hInvMass.fill(HIST("h3KstarInvMassMixed"), centrality, mother.Pt(), mother.M());
      }
    } else {
      if (!isMix) {
        if (calcLikeSign && (mother.Rapidity() > selectionConfig.motherRapidityMin && mother.Rapidity() < selectionConfig.motherRapidityMax)) {
          if (track1.sign() > 0 && track2.sign() > 0) {
            hInvMass.fill(HIST("h3KstarInvMasslikeSignPP"), centrality, mother.Pt(), mother.M());
          } else if (track1.sign() < 0 && track2.sign() < 0) {
            hInvMass.fill(HIST("h3KstarInvMasslikeSignMM"), centrality, mother.Pt(), mother.M());
          }
        }
      }
    }
  }

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection
  // requirements
  // Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < selectionConfig.cfgVrtxZCut);
  Filter acceptanceFilter = (nabs(aod::track::eta) < selectionConfig.cfgCutEta && nabs(aod::track::pt) > selectionConfig.cfgCutPT);
  Filter fDCAcutFilter = (nabs(aod::track::dcaXY) < selectionConfig.cfgCutDCAxy) && (nabs(aod::track::dcaZ) < selectionConfig.cfgCutDCAz);

  // using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::CentFV0As, aod::PVMults>>; // aod::CentNGlobals, aod::CentNTPVs, aod::CentMFTs
  //  using EventCandidatesMC = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::FT0Mults, aod::PVMults, aod::CentFV0As>>;
  // using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTOFbeta, aod::TrackSelectionExtension>>;
  // using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::McTrackLabels, aod::pidTOFbeta, aod::TrackSelectionExtension>>;
  using EventCandidatesMix = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::CentFV0As>>; // aod::CentNGlobals, aod::CentNTPVs, aod::CentMFTs

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::CentFV0As, aod::PVMults>; // aod::CentNGlobals, aod::CentNTPVs, aod::CentMFTs
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTOFbeta, aod::TrackSelectionExtension>;

  using EventMCGenerated = soa::Join<aod::McCollisions, aod::MultMCExtras>; // aod::CentNGlobals, aod::CentNTPVs, aod::CentMFTs
  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::MultZeqs, aod::FT0Mults, aod::PVMults, aod::CentFV0As>;
  using TrackCandidatesMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::McTrackLabels, aod::pidTOFbeta, aod::TrackSelectionExtension>;

  void processSE(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {
    int occupancy = collision.trackOccupancyInTimeRange();
    hEventSelection.fill(HIST("hOccupancy"), occupancy);

    if (!selectionEvent(collision, true)) { // fill data event cut histogram
      return;
    }

    centrality = -1;

    if (selectCentEstimator == kFT0M) {
      centrality = collision.centFT0M();
    } else if (selectCentEstimator == kFT0A) {
      centrality = collision.centFT0A();
    } else if (selectCentEstimator == kFT0C) {
      centrality = collision.centFT0C();
    } else if (selectCentEstimator == kFV0A) {
      centrality = collision.centFV0A();
    } else {
      centrality = collision.centFT0M(); // default
    }

    /* else if (selectCentEstimator == 4) {
      centrality = collision.centMFT();
    } */
    /* else if (selectCentEstimator == 5) {
      centrality = collision.centNGlobal();
    } */
    /* else if (selectCentEstimator == 6) {
      centrality = collision.centNTPV();
    } */

    // Fill the event counter
    if (cQAevents) {
      hEventSelection.fill(HIST("hVertexZ"), collision.posZ());
      hEventSelection.fill(HIST("hCentrality"), centrality);
    }

    for (const auto& [track1, track2] : combinations(CombinationsFullIndexPolicy(tracks, tracks))) {
      if (!selectionTrack(track1) || !selectionTrack(track2)) {
        continue;
      }

      if (track1.globalIndex() == track2.globalIndex())
        continue;

      if (!selectionPair(track1, track2)) {
        continue;
      }

      if (cQAplots) {
        hOthers.fill(HIST("hCRFC_before"), track1.tpcCrossedRowsOverFindableCls());
        hOthers.fill(HIST("dE_by_dx_TPC"), track1.p(), track1.tpcSignal());

        if (track1.sign() < 0) {
          hPID.fill(HIST("Before/hTPCnsigKa_Neg_mult_pt"), track1.tpcNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTPCnsigPi_Neg_mult_pt"), track1.tpcNSigmaPi(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTOFnsigKa_Neg_mult_pt"), track1.tofNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTOFnsigPi_Neg_mult_pt"), track1.tofNSigmaPi(), centrality, track1.pt());
        } else if (track1.sign() > 0) {
          hPID.fill(HIST("Before/hTPCnsigKa_Pos_mult_pt"), track1.tpcNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTPCnsigPi_Pos_mult_pt"), track1.tpcNSigmaPi(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTOFnsigKa_Pos_mult_pt"), track1.tofNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTOFnsigPi_Pos_mult_pt"), track1.tofNSigmaPi(), centrality, track1.pt());
        }

        hPID.fill(HIST("Before/hNsigma_TPC_TOF_Ka_pt"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
        hPID.fill(HIST("Before/hNsigma_TPC_TOF_Pi_pt"), track1.tpcNSigmaPi(), track1.tofNSigmaPi(), track1.pt());
      }

      if (cQAevents) {
        hOthers.fill(HIST("hDcaxy_cent_pt"), track1.dcaXY(), centrality, track1.pt());
        hOthers.fill(HIST("hDcaz_cent_pt"), track1.dcaZ(), centrality, track1.pt());
      }

      // since we are using combinations full index policy, so repeated pairs are allowed, so we can check one with Kaon and other with pion
      if ((!selectionConfig.isApplypTdepPID && !selectionConfig.isApplypTdepPIDwTOF) && !selectionPID(track1, 1)) // Track 1 is checked with Kaon
        continue;
      if ((!selectionConfig.isApplypTdepPID && !selectionConfig.isApplypTdepPIDwTOF) && !selectionPID(track2, 0)) // Track 2 is checked with Pion
        continue;

      if (selectionConfig.isApplypTdepPID && !selectionPIDpTdep(track1, 1)) // Track 1 is checked with Kaon
        continue;
      if (selectionConfig.isApplypTdepPID && !selectionPIDpTdep(track2, 0)) // Track 2 is checked with Pion
        continue;

      if (selectionConfig.isApplypTdepPIDwTOF && !selectionPIDpTdepTOF(track1, 1)) // Track 1 is checked with Kaon
        continue;
      if (selectionConfig.isApplypTdepPIDwTOF && !selectionPIDpTdepTOF(track2, 0)) // Track 2 is checked with Pion
        continue;

      if (selectionConfig.isApplyMID && (selectionMID(track1, 0) || selectionMID(track2, 1)))
        continue;

      if (selectionConfig.isApplypTdepMID && (selectionMIDpTdep(track1, 0) || selectionMIDpTdep(track2, 1)))
        continue;

      if (cQAplots) {
        hOthers.fill(HIST("hEta_after"), track1.eta());
        hOthers.fill(HIST("hCRFC_after"), track1.tpcCrossedRowsOverFindableCls());
        hOthers.fill(HIST("hDcaxyPi"), track2.dcaXY());
        hOthers.fill(HIST("hDcaxyKa"), track1.dcaXY());
        hOthers.fill(HIST("hDcazPi"), track2.dcaZ());
        hOthers.fill(HIST("hDcazKa"), track1.dcaZ());

        if (track1.sign() < 0) {
          hPID.fill(HIST("After/hTPCnsigKa_Neg_mult_pt"), track1.tpcNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("After/hTOFnsigKa_Neg_mult_pt"), track1.tofNSigmaKa(), centrality, track1.pt());
        } else if (track1.sign() > 0) {
          hPID.fill(HIST("After/hTPCnsigKa_Pos_mult_pt"), track1.tpcNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("After/hTOFnsigKa_Pos_mult_pt"), track1.tofNSigmaKa(), centrality, track1.pt());
        }

        if (track2.sign() < 0) {
          hPID.fill(HIST("After/hTPCnsigPi_Neg_mult_pt"), track2.tpcNSigmaPi(), centrality, track2.pt());
          hPID.fill(HIST("After/hTOFnsigPi_Neg_mult_pt"), track2.tofNSigmaPi(), centrality, track2.pt());
        } else if (track2.sign() > 0) {
          hPID.fill(HIST("After/hTPCnsigPi_Pos_mult_pt"), track2.tpcNSigmaPi(), centrality, track2.pt());
          hPID.fill(HIST("After/hTOFnsigPi_Pos_mult_pt"), track2.tofNSigmaPi(), centrality, track2.pt());
        }

        hPID.fill(HIST("After/hNsigma_TPC_TOF_Ka_pt"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
        hPID.fill(HIST("After/hNsigma_TPC_TOF_Pi_pt"), track2.tpcNSigmaPi(), track2.tofNSigmaPi(), track2.pt());
      }

      daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
      daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
      mother = daughter1 + daughter2; // Kstar meson

      /* if (selectionConfig.isApplyCutsOnMother) {
        if (mother.Pt() >= selectionConfig.cMaxPtMotherCut) // excluding candidates in overflow
          continue;
        if (mother.M() >= selectionConfig.cMaxMinvMotherCut) // excluding candidates in overflow
          continue;
      } */

      hOthers.fill(HIST("hKstar_rap_pt"), mother.Rapidity(), mother.Pt());
      hOthers.fill(HIST("hKstar_eta_pt"), mother.Eta(), mother.Pt());

      isMix = false;
      fillInvMass(daughter1, daughter2, mother, centrality, isMix, track1, track2);
    }
  }
  PROCESS_SWITCH(Kstar892LightIon, processSE, "Process Same event", true);

  void processSEMC(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&, aod::McCollisions const& /*mcCollisions*/)
  {
    // auto oldindex = -999;
    if (!collision.has_mcCollision()) {
      return;
    }
    int occupancy = collision.trackOccupancyInTimeRange();
    hEventSelection.fill(HIST("hOccupancy"), occupancy);

    if (!selectionEvent(collision, false)) { // don't fill event cut histogram
      return;
    }

    centrality = -1;

    if (selectCentEstimator == kFT0M) {
      centrality = collision.centFT0M();
    } else if (selectCentEstimator == kFT0A) {
      centrality = collision.centFT0A();
    } else if (selectCentEstimator == kFT0C) {
      centrality = collision.centFT0C();
    } else if (selectCentEstimator == kFV0A) {
      centrality = collision.centFV0A();
    } else {
      centrality = collision.centFT0M(); // default
    }

    // Fill the event counter
    if (cQAevents) {
      hEventSelection.fill(HIST("hVertexZ"), collision.posZ());
      hEventSelection.fill(HIST("hCentrality"), centrality);
    }

    for (const auto& [track1, track2] : combinations(CombinationsFullIndexPolicy(tracks, tracks))) {
      if (!selectionTrack(track1) || !selectionTrack(track2))
        continue;

      const auto mctrack1 = track1.mcParticle();
      const auto mctrack2 = track2.mcParticle();

      if (!track1.has_mcParticle() || !track2.has_mcParticle())
        continue; // skip if no MC particle associated

      if (!mctrack1.isPhysicalPrimary() || !mctrack2.isPhysicalPrimary())
        continue;

      if (track1.globalIndex() == track2.globalIndex())
        continue;

      if (cQAplots) {
        hOthers.fill(HIST("hCRFC_before"), track1.tpcCrossedRowsOverFindableCls());
        hOthers.fill(HIST("dE_by_dx_TPC"), track1.p(), track1.tpcSignal());

        if (track1.sign() < 0) {
          hPID.fill(HIST("Before/hTPCnsigKa_Neg_mult_pt"), track1.tpcNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTPCnsigPi_Neg_mult_pt"), track1.tpcNSigmaPi(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTOFnsigKa_Neg_mult_pt"), track1.tofNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTOFnsigPi_Neg_mult_pt"), track1.tofNSigmaPi(), centrality, track1.pt());
        } else {
          hPID.fill(HIST("Before/hTPCnsigKa_Pos_mult_pt"), track1.tpcNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTPCnsigPi_Pos_mult_pt"), track1.tpcNSigmaPi(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTOFnsigKa_Pos_mult_pt"), track1.tofNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTOFnsigPi_Pos_mult_pt"), track1.tofNSigmaPi(), centrality, track1.pt());
        }

        hPID.fill(HIST("Before/hNsigma_TPC_TOF_Ka_pt"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
        hPID.fill(HIST("Before/hNsigma_TPC_TOF_Pi_pt"), track1.tpcNSigmaPi(), track1.tofNSigmaPi(), track1.pt());
      }

      if (cQAevents) {
        hOthers.fill(HIST("hDcaxy_cent_pt"), track1.dcaXY(), centrality, track1.pt());
        hOthers.fill(HIST("hDcaz_cent_pt"), track1.dcaZ(), centrality, track1.pt());
      }

      // since we are using combinations full index policy, so repeated pairs are allowed, so we can check one with Kaon and other with pion
      if ((!selectionConfig.isApplypTdepPID && !selectionConfig.isApplypTdepPIDwTOF) && (!selectionPID(track1, 1) || !selectionPID(track2, 0))) // Track 1 is checked with Kaon, track 2 is checked with Pion
        continue;

      if (selectionConfig.isApplypTdepPID && (!selectionPIDpTdep(track1, 1) || !selectionPIDpTdep(track2, 0))) // Track 1 is checked with Kaon, track 2 is checked with Pion
        continue;

      if (selectionConfig.isApplypTdepPIDwTOF && (!selectionPIDpTdepTOF(track1, 1) || !selectionPIDpTdepTOF(track2, 0))) // Track 1 is checked with Kaon, track 2 is checked with Pion
        continue;

      if (selectionConfig.isApplyMID && (selectionMID(track1, 0) || selectionMID(track2, 1)))
        continue;

      if (selectionConfig.isApplypTdepMID && (selectionMIDpTdep(track1, 0) || selectionMIDpTdep(track2, 1)))
        continue;

      if (cQAplots) {
        hOthers.fill(HIST("hEta_after"), track1.eta());
        hOthers.fill(HIST("hCRFC_after"), track1.tpcCrossedRowsOverFindableCls());
        hOthers.fill(HIST("hDcaxyPi"), track2.dcaXY());
        hOthers.fill(HIST("hDcaxyKa"), track1.dcaXY());
        hOthers.fill(HIST("hDcazPi"), track2.dcaZ());
        hOthers.fill(HIST("hDcazKa"), track1.dcaZ());

        if (track1.sign() < 0) {
          hPID.fill(HIST("After/hTPCnsigKa_Neg_mult_pt"), track1.tpcNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("After/hTOFnsigKa_Neg_mult_pt"), track1.tofNSigmaKa(), centrality, track1.pt());
        } else if (track1.sign() > 0) {
          hPID.fill(HIST("After/hTPCnsigKa_Pos_mult_pt"), track1.tpcNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("After/hTOFnsigKa_Pos_mult_pt"), track1.tofNSigmaKa(), centrality, track1.pt());
        }

        if (track2.sign() < 0) {
          hPID.fill(HIST("After/hTPCnsigPi_Neg_mult_pt"), track2.tpcNSigmaPi(), centrality, track2.pt());
          hPID.fill(HIST("After/hTOFnsigPi_Neg_mult_pt"), track2.tofNSigmaPi(), centrality, track2.pt());
        } else if (track2.sign() > 0) {
          hPID.fill(HIST("After/hTPCnsigPi_Pos_mult_pt"), track2.tpcNSigmaPi(), centrality, track2.pt());
          hPID.fill(HIST("After/hTOFnsigPi_Pos_mult_pt"), track2.tofNSigmaPi(), centrality, track2.pt());
        }

        hPID.fill(HIST("After/hNsigma_TPC_TOF_Ka_pt"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
        hPID.fill(HIST("After/hNsigma_TPC_TOF_Pi_pt"), track2.tpcNSigmaPi(), track2.tofNSigmaPi(), track2.pt());
      }
      daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
      daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
      mother = daughter1 + daughter2; // Kstar meson

      /* if (selectionConfig.isApplyCutsOnMother) {
        if (mother.Pt() >= selectionConfig.cMaxPtMotherCut) // excluding candidates in overflow
          continue;
        if (mother.M() >= selectionConfig.cMaxMinvMotherCut) // excluding candidates in overflow
          continue;
      } */

      hOthers.fill(HIST("hKstar_rap_pt"), mother.Rapidity(), mother.Pt());
      hOthers.fill(HIST("hKstar_eta_pt"), mother.Eta(), mother.Pt());

      isMix = false;
      fillInvMass(daughter1, daughter2, mother, centrality, isMix, track1, track2);
    }
  }
  PROCESS_SWITCH(Kstar892LightIon, processSEMC, "Process same event in MC", false);

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for ME mixing"};
  // ConfigurableAxis axisCentralityClass{"axisCentralityClass", {10, 0, 100}, "centrality percentile for ME mixing"};
  ConfigurableAxis axisCentrality{"axisCentrality", {2000, 0, 10000}, "TPC centrality axis for ME mixing"};

  // using BinningTypeTPCcentrality = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  using BinningTypeFT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypeFT0A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0A>;
  using BinningTypeFT0C = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  using BinningTypeFV0A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFV0A>;

  BinningTypeFT0M binningOnFT0M{{axisVertex, axisCentrality}, true};
  BinningTypeFT0A binningOnFT0A{{axisVertex, axisCentrality}, true};
  BinningTypeFT0C binningOnFT0C{{axisVertex, axisCentrality}, true};
  BinningTypeFV0A binningOnFV0A{{axisVertex, axisCentrality}, true};

  SameKindPair<EventCandidates, TrackCandidates, BinningTypeFT0M> pair1{binningOnFT0M, selectionConfig.cfgNoMixedEvents, -1, &cache};
  SameKindPair<EventCandidates, TrackCandidates, BinningTypeFT0A> pair2{binningOnFT0A, selectionConfig.cfgNoMixedEvents, -1, &cache};
  SameKindPair<EventCandidates, TrackCandidates, BinningTypeFT0C> pair3{binningOnFT0C, selectionConfig.cfgNoMixedEvents, -1, &cache};
  SameKindPair<EventCandidates, TrackCandidates, BinningTypeFV0A> pair4{binningOnFV0A, selectionConfig.cfgNoMixedEvents, -1, &cache};

  void processME(EventCandidatesMix const&, TrackCandidates const&)
  {
    // Map estimator to pair and centrality accessor
    auto runMixing = [&](auto& pair, auto centralityGetter) {
      for (const auto& [c1, tracks1, c2, tracks2] : pair) {

        if (!selectionEvent(c1, false) || !selectionEvent(c2, false)) { // don't fill event cut histogram
          continue;
        }

        centrality = centralityGetter(c1);

        for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
          if (!selectionTrack(t1) || !selectionTrack(t2))
            continue;

          if (!selectionPair(t1, t2)) {
            continue;
          }

          if ((!selectionConfig.isApplypTdepPID && !selectionConfig.isApplypTdepPIDwTOF) && (!selectionPID(t1, 1) || !selectionPID(t2, 0))) // Track 1 is checked with Kaon, track 2 is checked with Pion
            continue;

          if (selectionConfig.isApplypTdepPID && (!selectionPIDpTdep(t1, 1) || !selectionPIDpTdep(t2, 0))) // Track 1 is checked with Kaon, track 2 is checked with Pion
            continue;

          if (selectionConfig.isApplypTdepPIDwTOF && (!selectionPIDpTdepTOF(t1, 1) || !selectionPIDpTdepTOF(t2, 0))) // Track 1 is checked with Kaon, track 2 is checked with Pion
            continue;

          if (selectionConfig.isApplyMID && (selectionMID(t1, 0) || selectionMID(t2, 1)))
            continue;

          if (selectionConfig.isApplypTdepMID && (selectionMIDpTdep(t1, 0) || selectionMIDpTdep(t2, 1)))
            continue;

          daughter1 = ROOT::Math::PxPyPzMVector(t1.px(), t1.py(), t1.pz(), massKa);
          daughter2 = ROOT::Math::PxPyPzMVector(t2.px(), t2.py(), t2.pz(), massPi);
          mother = daughter1 + daughter2;

          isMix = true;

          fillInvMass(daughter1, daughter2, mother, centrality, isMix, t1, t2);
        }
      }
    };

    // Call mixing based on selected estimator
    if (selectCentEstimator == kFT0M) {
      runMixing(pair1, [](const auto& c) { return c.centFT0M(); });
    } else if (selectCentEstimator == kFT0A) {
      runMixing(pair2, [](const auto& c) { return c.centFT0A(); });
    } else if (selectCentEstimator == kFT0C) {
      runMixing(pair3, [](const auto& c) { return c.centFT0C(); });
    } else if (selectCentEstimator == kFV0A) {
      runMixing(pair4, [](const auto& c) { return c.centFV0A(); });
    }
  }
  PROCESS_SWITCH(Kstar892LightIon, processME, "Process Mixed event", true);

  using BinningTypeMCFT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypeMCFT0A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0A>;
  using BinningTypeMCFT0C = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  using BinningTypeMCFV0A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFV0A>;

  BinningTypeMCFT0M binningOnMCFT0M{{axisVertex, axisCentrality}, true};
  BinningTypeMCFT0A binningOnMCFT0A{{axisVertex, axisCentrality}, true};
  BinningTypeMCFT0C binningOnMCFT0C{{axisVertex, axisCentrality}, true};
  BinningTypeMCFV0A binningOnMCFV0A{{axisVertex, axisCentrality}, true};

  SameKindPair<EventCandidatesMC, TrackCandidatesMC, BinningTypeMCFT0M> pairmc1{binningOnMCFT0M, selectionConfig.cfgNoMixedEvents, -1, &cache};
  SameKindPair<EventCandidatesMC, TrackCandidatesMC, BinningTypeMCFT0A> pairmc2{binningOnMCFT0A, selectionConfig.cfgNoMixedEvents, -1, &cache};
  SameKindPair<EventCandidatesMC, TrackCandidatesMC, BinningTypeMCFT0C> pairmc3{binningOnMCFT0C, selectionConfig.cfgNoMixedEvents, -1, &cache};
  SameKindPair<EventCandidatesMC, TrackCandidatesMC, BinningTypeMCFV0A> pairmc4{binningOnMCFV0A, selectionConfig.cfgNoMixedEvents, -1, &cache};

  void processMEMC(EventCandidatesMC const&, TrackCandidatesMC const&, aod::McParticles const&, aod::McCollisions const&)
  {
    auto runMixing = [&](auto& pair, auto centralityGetter) {
      for (const auto& [c1, tracks1, c2, tracks2] : pair) {

        if (!selectionEvent(c1, false) || !selectionEvent(c2, false)) { // don't fill event cut histogram
          continue;
        }

        if (!c1.has_mcCollision() || !c2.has_mcCollision()) {
          continue; // skip if no MC collision associated
        }

        centrality = centralityGetter(c1);

        for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
          if (!selectionTrack(t1) || !selectionTrack(t2))
            continue;

          if (!selectionPair(t1, t2)) {
            continue;
          }

          if ((!selectionConfig.isApplypTdepPID && !selectionConfig.isApplypTdepPIDwTOF) && (!selectionPID(t1, 1) || !selectionPID(t2, 0))) // Track 1 is checked with Kaon, track 2 is checked with Pion
            continue;

          if (selectionConfig.isApplypTdepPID && (!selectionPIDpTdep(t1, 1) || !selectionPIDpTdep(t2, 0))) // Track 1 is checked with Kaon, track 2 is checked with Pion
            continue;

          if (selectionConfig.isApplypTdepPIDwTOF && (!selectionPIDpTdepTOF(t1, 1) || !selectionPIDpTdepTOF(t2, 0))) // Track 1 is checked with Kaon, track 2 is checked with Pion
            continue;

          if (selectionConfig.isApplyMID && (selectionMID(t1, 0) || selectionMID(t2, 1)))
            continue;

          if (selectionConfig.isApplypTdepMID && (selectionMIDpTdep(t1, 0) || selectionMIDpTdep(t2, 1)))
            continue;

          if (!t1.has_mcParticle() || !t2.has_mcParticle()) {
            continue; // skip if no MC particle associated
          }

          const auto mctrack1 = t1.mcParticle();
          const auto mctrack2 = t2.mcParticle();

          if (!mctrack1.isPhysicalPrimary() || !mctrack2.isPhysicalPrimary()) {
            continue;
          }

          daughter1 = ROOT::Math::PxPyPzMVector(t1.px(), t1.py(), t1.pz(), massKa);
          daughter2 = ROOT::Math::PxPyPzMVector(t2.px(), t2.py(), t2.pz(), massPi);
          mother = daughter1 + daughter2;

          isMix = true;

          fillInvMass(daughter1, daughter2, mother, centrality, isMix, t1, t2);
        }
      }
    };
    // Call mixing based on selected estimator
    if (selectCentEstimator == kFT0M) {
      runMixing(pairmc1, [](const auto& c) { return c.centFT0M(); });
    } else if (selectCentEstimator == kFT0A) {
      runMixing(pairmc2, [](const auto& c) { return c.centFT0A(); });
    } else if (selectCentEstimator == kFT0C) {
      runMixing(pairmc3, [](const auto& c) { return c.centFT0C(); });
    } else if (selectCentEstimator == kFV0A) {
      runMixing(pairmc4, [](const auto& c) { return c.centFV0A(); });
    }
  }
  PROCESS_SWITCH(Kstar892LightIon, processMEMC, "Process mixed-event in MC", false);

  void processGen(EventMCGenerated::iterator const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<EventCandidatesMC>& collisions)
  {

    if (selectionConfig.isApplyMCGenInelgt0 && !mcCollision.isInelGt0()) {
      return;
    }

    if (selectionConfig.isApplyMCGenTVX && !(mcCollision.multMCFT0C() > 0 && mcCollision.multMCFT0A() > 0)) {
      return;
    }

    if (selectionConfig.isApplyMCGenVz && std::abs(mcCollision.posZ()) >= selectionConfig.cfgVrtxZCut) {
      return;
    }

    std::vector<int64_t> selectedEvents(collisions.size());
    int nevts = 0;
    centrality = -1.0;

    for (const auto& collision : collisions) {

      hMC.fill(HIST("Gen/hGenNo"), 0.5);

      if (!selectionEvent(collision, false)) { // don't fill event cut histogram
        continue;
      }
      centrality = collision.centFT0M();

      if (selectCentEstimator == kFT0M) {
        centrality = collision.centFT0M();
      } else if (selectCentEstimator == kFT0A) {
        centrality = collision.centFT0A();
      } else if (selectCentEstimator == kFT0C) {
        centrality = collision.centFT0C();
      } else if (selectCentEstimator == kFV0A) {
        centrality = collision.centFV0A();
      } else {
        centrality = collision.centFT0M(); // default
      }
      hMC.fill(HIST("Gen/h1GenCent"), centrality);

      int occupancy = collision.trackOccupancyInTimeRange();
      hEventSelection.fill(HIST("hOccupancy"), occupancy);

      selectedEvents[nevts++] = collision.mcCollision_as<EventMCGenerated>().globalIndex();
    }
    selectedEvents.resize(nevts);

    for (const auto& mcParticle : mcParticles) {
      if ((mcParticle.y() > selectionConfig.motherRapidityMin && mcParticle.y() < selectionConfig.motherRapidityMax) && std::abs(mcParticle.pdgCode()) == o2::constants::physics::kK0Star892) {
        hMC.fill(HIST("Gen/hAllKstarGenCollisisons"), mcParticle.pt(), centrality);
      }
    }

    const auto evtReconstructedAndSelected = std::find(selectedEvents.begin(), selectedEvents.end(), mcCollision.globalIndex()) != selectedEvents.end();
    hMC.fill(HIST("Gen/hAllGenCollisions"), centrality);
    if (!evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      return;
    }
    hMC.fill(HIST("Gen/hGenNo"), 1.5);

    hMC.fill(HIST("Gen/hAllGenCollisions1Rec"), centrality);

    for (const auto& mcParticle : mcParticles) {

      if (mcParticle.y() < selectionConfig.motherRapidityMin || mcParticle.y() > selectionConfig.motherRapidityMax) {
        continue;
      }

      /* if (selectionConfig.isApplyCutsOnMother) {
        if (mcParticle.pt() >= selectionConfig.cMaxPtMotherCut) // excluding candidates in overflow
          continue;
        if ((std::sqrt(mcParticle.e() * mcParticle.e() - mcParticle.p() * mcParticle.p())) >= selectionConfig.cMaxMinvMotherCut) // excluding candidates in overflow
          continue;
      } */

      if (std::abs(mcParticle.pdgCode()) != o2::constants::physics::kK0Star892) {
        continue;
      }
      hMC.fill(HIST("Gen/hAllKstarGenCollisisons1Rec"), mcParticle.pt(), centrality);

      auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
      if (kDaughters.size() != noOfDaughters) {
        continue;
      }

      bool hasPos = false;
      bool hasNeg = false;

      auto passkaon = false;
      auto passpion = false;

      for (const auto& kCurrentDaughter : kDaughters) {
        if (!kCurrentDaughter.isPhysicalPrimary()) {
          continue;
        }

        int pdgDau = kCurrentDaughter.pdgCode();
        int charge = (pdgDau > 0) - (pdgDau < 0);

        if (charge > 0)
          hasPos = true;
        if (charge < 0)
          hasNeg = true;

        if (std::abs(pdgDau) == PDG_t::kKPlus) {
          passkaon = true;
        } else if (std::abs(pdgDau) == PDG_t::kPiPlus) {
          passpion = true;
        }
      }

      if ((passkaon && passpion) && (hasPos && hasNeg)) {
        hMC.fill(HIST("Gen/hk892GenpT"), mcParticle.pt(), centrality);
      }
    }
  }
  PROCESS_SWITCH(Kstar892LightIon, processGen, "Process Generated", false);

  void processRec(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&, EventMCGenerated const&)
  {
    if (!collision.has_mcCollision()) {
      return;
    }

    centrality = collision.centFT0M();

    if (selectCentEstimator == kFT0M) {
      centrality = collision.centFT0M();
    } else if (selectCentEstimator == kFT0A) {
      centrality = collision.centFT0A();
    } else if (selectCentEstimator == kFT0C) {
      centrality = collision.centFT0C();
    } else if (selectCentEstimator == kFV0A) {
      centrality = collision.centFV0A();
    } else {
      centrality = collision.centFT0M(); // default
    }

    hMC.fill(HIST("Rec/hAllRecCollisions"), centrality);

    if (!selectionEvent(collision, true)) { // fill MC event cut histogram
      return;
    }

    hMC.fill(HIST("Rec/h1RecCent"), centrality);

    if (cQAevents) {
      hEventSelection.fill(HIST("hVertexZ"), collision.posZ());
    }

    auto oldindex = -999;

    for (const auto& [track1, track2] : combinations(CombinationsFullIndexPolicy(tracks, tracks))) {
      if (!selectionTrack(track1) || !selectionTrack(track2)) {
        continue;
      }

      if (!track1.has_mcParticle() || !track2.has_mcParticle()) {
        continue;
      }

      if (track1.index() <= track2.index())
        continue;

      if (!selectionPair(track1, track2)) {
        continue;
      }

      if (cQAevents) {
        hOthers.fill(HIST("hDcaxy_cent_pt"), track1.dcaXY(), centrality, track1.pt());
        hOthers.fill(HIST("hDcaz_cent_pt"), track1.dcaZ(), centrality, track1.pt());
      }

      if (track1.sign() * track2.sign() >= 0) {
        continue;
      }

      const auto mctrack1 = track1.mcParticle();
      const auto mctrack2 = track2.mcParticle();
      if (!mctrack1.isPhysicalPrimary() || !mctrack2.isPhysicalPrimary()) {
        continue;
      }

      int track1PDG = std::abs(mctrack1.pdgCode());
      int track2PDG = std::abs(mctrack2.pdgCode());

      if (cQAplots) {
        if (mctrack2.pdgCode() == PDG_t::kPiPlus) { // pion
          hPID.fill(HIST("Before/hTPCnsigPi_Pos_mult_pt"), track2.tpcNSigmaPi(), centrality, track2.pt());
          hPID.fill(HIST("Before/hTOFnsigPi_Pos_mult_pt"), track2.tofNSigmaPi(), centrality, track2.pt());
        }
        if (mctrack2.pdgCode() == PDG_t::kKPlus) { // kaon
          hPID.fill(HIST("Before/hTPCnsigKa_Pos_mult_pt"), track2.tpcNSigmaKa(), centrality, track2.pt());
          hPID.fill(HIST("Before/hTOFnsigKa_Pos_mult_pt"), track2.tofNSigmaKa(), centrality, track2.pt());
        }
        if (mctrack2.pdgCode() == PDG_t::kPiMinus) { // negative track pion
          hPID.fill(HIST("Before/hTPCnsigPi_Neg_mult_pt"), track2.tpcNSigmaPi(), centrality, track2.pt());
          hPID.fill(HIST("Before/hTOFnsigPi_Neg_mult_pt"), track2.tofNSigmaPi(), centrality, track2.pt());
        }
        if (mctrack2.pdgCode() == PDG_t::kKMinus) { // negative track kaon
          hPID.fill(HIST("Before/hTPCnsigKa_Neg_mult_pt"), track2.tpcNSigmaKa(), centrality, track2.pt());
          hPID.fill(HIST("Before/hTOFnsigKa_Neg_mult_pt"), track2.tofNSigmaKa(), centrality, track2.pt());
        }
        if (std::abs(mctrack1.pdgCode()) == PDG_t::kKPlus && std::abs(mctrack2.pdgCode()) == PDG_t::kPiPlus) {
          hPID.fill(HIST("Before/hNsigma_TPC_TOF_Ka_pt"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
          hPID.fill(HIST("Before/hNsigma_TPC_TOF_Pi_pt"), track2.tpcNSigmaPi(), track2.tofNSigmaPi(), track2.pt());
        }
      }

      if (!(track1PDG == PDG_t::kPiPlus && track2PDG == PDG_t::kKPlus) && !(track1PDG == PDG_t::kKPlus && track2PDG == PDG_t::kPiPlus)) {
        continue;
      }

      for (const auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
        for (const auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
          if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
            continue;
          }

          if (mothertrack1.globalIndex() != mothertrack2.globalIndex()) {
            continue;
          }

          if (!mothertrack1.producedByGenerator()) {
            continue;
          }

          if (mothertrack1.y() < selectionConfig.motherRapidityMin || mothertrack1.y() > selectionConfig.motherRapidityMax) {
            continue;
          }

          if (std::abs(mothertrack1.pdgCode()) != o2::constants::physics::kK0Star892) {
            continue;
          }

          if (track1PDG == PDG_t::kPiPlus) {
            if ((!selectionConfig.isApplypTdepPID && !selectionConfig.isApplypTdepPIDwTOF) && !(selectionPID(track1, 0) && selectionPID(track2, 1))) { // pion and kaon
              continue;
            } else if (selectionConfig.isApplypTdepPID && !(selectionPIDpTdep(track1, 0) && selectionPIDpTdep(track2, 1))) { // pion and kaon
              continue;
            } else if (selectionConfig.isApplypTdepPIDwTOF && !(selectionPIDpTdepTOF(track1, 0) && selectionPIDpTdepTOF(track2, 1))) {
              continue;
            }

            if (selectionConfig.isApplyMID && (selectionMID(track1, 1) || selectionMID(track2, 0)))
              continue;

            if (selectionConfig.isApplypTdepMID && (selectionMIDpTdep(track1, 1) || selectionMIDpTdep(track2, 0)))
              continue;

            if (cQAplots) {
              if (track1.sign() < 0 && track2.sign() > 0) {
                hPID.fill(HIST("After/hTPCnsigPi_Neg_mult_pt"), track1.tpcNSigmaPi(), centrality, track1.pt());
                hPID.fill(HIST("After/hTOFnsigPi_Neg_mult_pt"), track1.tofNSigmaPi(), centrality, track1.pt());
                hPID.fill(HIST("After/hTPCnsigKa_Pos_mult_pt"), track2.tpcNSigmaKa(), centrality, track2.pt());
                hPID.fill(HIST("After/hTOFnsigKa_Pos_mult_pt"), track2.tofNSigmaKa(), centrality, track2.pt());
              } else {
                hPID.fill(HIST("After/hTPCnsigPi_Pos_mult_pt"), track1.tpcNSigmaPi(), centrality, track1.pt());
                hPID.fill(HIST("After/hTOFnsigPi_Pos_mult_pt"), track1.tofNSigmaPi(), centrality, track1.pt());
                hPID.fill(HIST("After/hTPCnsigKa_Neg_mult_pt"), track2.tpcNSigmaKa(), centrality, track2.pt());
                hPID.fill(HIST("After/hTOFnsigKa_Neg_mult_pt"), track2.tofNSigmaKa(), centrality, track2.pt());
              }
            }

          } else if (track1PDG == PDG_t::kKPlus) {
            if ((!selectionConfig.isApplypTdepPID && !selectionConfig.isApplypTdepPIDwTOF) && !(selectionPID(track1, 1) && selectionPID(track2, 0))) { // kaon and pion
              continue;
            } else if (selectionConfig.isApplypTdepPID && !(selectionPIDpTdep(track1, 1) && selectionPIDpTdep(track2, 0))) { // kaon and pion
              continue;
            } else if (selectionConfig.isApplypTdepPIDwTOF && !(selectionPIDpTdepTOF(track1, 1) && selectionPIDpTdepTOF(track2, 0))) {
              continue;
            }

            if (selectionConfig.isApplyMID && (selectionMID(track1, 0) || selectionMID(track2, 1)))
              continue;

            if (selectionConfig.isApplypTdepMID && (selectionMIDpTdep(track1, 0) || selectionMIDpTdep(track2, 1)))
              continue;

            if (cQAplots) {
              if (track1.sign() < 0 && track2.sign() > 0) {
                hPID.fill(HIST("After/hTPCnsigKa_Neg_mult_pt"), track1.tpcNSigmaKa(), centrality, track1.pt());
                hPID.fill(HIST("After/hTOFnsigKa_Neg_mult_pt"), track1.tofNSigmaKa(), centrality, track1.pt());
                hPID.fill(HIST("After/hTPCnsigPi_Pos_mult_pt"), track2.tpcNSigmaPi(), centrality, track2.pt());
                hPID.fill(HIST("After/hTOFnsigPi_Pos_mult_pt"), track2.tofNSigmaPi(), centrality, track2.pt());
              } else {
                hPID.fill(HIST("After/hTPCnsigKa_Pos_mult_pt"), track1.tpcNSigmaKa(), centrality, track1.pt());
                hPID.fill(HIST("After/hTOFnsigKa_Pos_mult_pt"), track1.tofNSigmaKa(), centrality, track1.pt());
                hPID.fill(HIST("After/hTPCnsigPi_Neg_mult_pt"), track2.tpcNSigmaPi(), centrality, track2.pt());
                hPID.fill(HIST("After/hTOFnsigPi_Neg_mult_pt"), track2.tofNSigmaPi(), centrality, track2.pt());
              }
            }
          }

          /* if (selectionConfig.isApplyCutsOnMother) {
            if (mothertrack1.pt() >= selectionConfig.cMaxPtMotherCut) // excluding candidates in overflow
              continue;
            if ((std::sqrt(mothertrack1.e() * mothertrack1.e() - mothertrack1.p() * mothertrack1.p())) >= selectionConfig.cMaxMinvMotherCut) // excluding candidates in overflow
              continue;
          } */

          if (selectionConfig.isAvoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
            hMC.fill(HIST("Rec/h1KSRecsplit"), mothertrack1.pt());
            continue;
          }

          oldindex = mothertrack1.globalIndex();

          if (track1PDG == PDG_t::kPiPlus) {
            daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massPi);
            daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
          } else if (track1PDG == PDG_t::kKPlus) {
            daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
            daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
          }

          mother = daughter1 + daughter2; // Kstar meson

          hMC.fill(HIST("Rec/h2KstarRecpt2"), mothertrack1.pt(), centrality, std::sqrt(mothertrack1.e() * mothertrack1.e() - mothertrack1.p() * mothertrack1.p()));

          if (mother.Rapidity() < selectionConfig.motherRapidityMin || mother.Rapidity() > selectionConfig.motherRapidityMax) {
            continue;
          }

          hMC.fill(HIST("Rec/h1KstarRecMass"), mother.M());
          hMC.fill(HIST("Rec/h2KstarRecpt1"), mother.Pt(), centrality, mother.M());
        }
      }
    }
  }
  PROCESS_SWITCH(Kstar892LightIon, processRec, "Process Reconstructed", false);

  void processEvtLossSigLossMC(EventMCGenerated::iterator const& mcCollision, const soa::SmallGroups<EventCandidatesMC>& recCollisions, aod::McParticles const& mcParticles)
  {
    if (selectionConfig.isApplyMCGenInelgt0 && !mcCollision.isInelGt0()) {
      return;
    }

    if (selectionConfig.isApplyMCGenTVX && !(mcCollision.multMCFT0C() > 0 && mcCollision.multMCFT0A() > 0)) {
      return;
    }

    if (selectionConfig.isApplyMCGenVz && std::abs(mcCollision.posZ()) >= selectionConfig.cfgVrtxZCut) {
      return;
    }

    auto impactPar = mcCollision.impactParameter();
    hMC.fill(HIST("ImpactCorr/hImpactParameterGen"), impactPar);

    bool isSelectedEvent = false;
    centrality = -1.f;

    for (const auto& RecCollision : recCollisions) {
      if (!RecCollision.has_mcCollision())
        continue;
      if (!selectionEvent(RecCollision, false)) // don't fill event cut histogram
        continue;

      if (selectCentEstimator == kFT0M) {
        centrality = RecCollision.centFT0M();
      } else if (selectCentEstimator == kFT0A) {
        centrality = RecCollision.centFT0A();
      } else if (selectCentEstimator == kFT0C) {
        centrality = RecCollision.centFT0C();
      } else if (selectCentEstimator == kFV0A) {
        centrality = RecCollision.centFV0A();
      } else {
        centrality = RecCollision.centFT0M(); // default
      }

      isSelectedEvent = true;
    }

    if (isSelectedEvent) {
      hMC.fill(HIST("ImpactCorr/hImpactParameterRec"), impactPar);
      hMC.fill(HIST("ImpactCorr/hImpactParvsCentrRec"), centrality, impactPar);
    }

    // Generated MC
    for (const auto& mcPart : mcParticles) {
      if ((mcPart.y() < selectionConfig.motherRapidityMin || mcPart.y() > selectionConfig.motherRapidityMax) || std::abs(mcPart.pdgCode()) != o2::constants::physics::kK0Star892)
        continue;

      // signal loss estimation
      hMC.fill(HIST("ImpactCorr/hKstarGenBeforeEvtSel"), mcPart.pt(), impactPar);
      if (isSelectedEvent) {
        // signal loss estimation
        hMC.fill(HIST("ImpactCorr/hKstarGenAfterEvtSel"), mcPart.pt(), impactPar);
      }
    } // end loop on gen particles
  }
  PROCESS_SWITCH(Kstar892LightIon, processEvtLossSigLossMC, "Process Signal Loss, Event Loss using impact parameter", false);

  void processLossMCMultiplicity(EventMCGenerated::iterator const& mcCollision, aod::McParticles const& mcParticles, soa::SmallGroups<EventCandidatesMC> const& recCollisions)
  {

    if (selectionConfig.isApplyMCGenInelgt0 && !mcCollision.isInelGt0()) {
      return;
    }

    if (selectionConfig.isApplyMCGenTVX && !(mcCollision.multMCFT0C() > 0 && mcCollision.multMCFT0A() > 0)) {
      return;
    }

    if (selectionConfig.isApplyMCGenVz && std::abs(mcCollision.posZ()) >= selectionConfig.cfgVrtxZCut) {
      return;
    }

    const int multMC = mcCollision.multMCNParticlesEta05();
    hMC.fill(HIST("LossMult/hMultMC"), multMC);

    bool isSelectedEvent = false;
    centrality = -1.f;

    for (auto const& collision : recCollisions) {

      if (!selectionEvent(collision, false))
        continue;

      if (selectCentEstimator == kFT0M) {
        centrality = collision.centFT0M();
      } else if (selectCentEstimator == kFT0A) {
        centrality = collision.centFT0A();
      } else if (selectCentEstimator == kFT0C) {
        centrality = collision.centFT0C();
      } else if (selectCentEstimator == kFV0A) {
        centrality = collision.centFV0A();
      } else {
        centrality = collision.centFT0M(); // default
      }

      isSelectedEvent = true;
    }

    hMC.fill(HIST("LossMult/hCentVsMultMC"), centrality, multMC);

    // Event loss histograms
    hMC.fill(HIST("LossMult/hGenEvt_vs_multMC"), multMC);

    if (isSelectedEvent) {
      hMC.fill(HIST("LossMult/hCentVsMultMC_EvtSel"), centrality, multMC);
      hMC.fill(HIST("LossMult/hGenEvtRecoEvt_vs_multMC"), multMC);
    }

    // Signal loss histograms
    for (auto const& mcPart : mcParticles) {

      if ((mcPart.y() < selectionConfig.motherRapidityMin || mcPart.y() > selectionConfig.motherRapidityMax) || std::abs(mcPart.pdgCode()) != o2::constants::physics::kK0Star892)
        continue;

      const float pt = mcPart.pt();

      hMC.fill(HIST("LossMult/hGenKstar_vs_pt_vs_multMC"), pt, multMC);

      if (isSelectedEvent) {
        hMC.fill(HIST("LossMult/hGenKstarRecoEvt_vs_pt_vs_multMC"), pt, multMC);
      }
    }
  }
  PROCESS_SWITCH(Kstar892LightIon, processLossMCMultiplicity, "Signal + Event loss (using MC multiplicity)", false);

  void processAllLossMC(EventMCGenerated::iterator const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<EventCandidatesMC>& recCollisions)
  {
    if (selectionConfig.isApplyMCGenInelgt0 && !mcCollision.isInelGt0()) {
      return;
    }

    if (selectionConfig.isApplyMCGenTVX && !(mcCollision.multMCFT0C() > 0 && mcCollision.multMCFT0A() > 0)) {
      return;
    }

    if (selectionConfig.isApplyMCGenVz && std::abs(mcCollision.posZ()) >= selectionConfig.cfgVrtxZCut) {
      return;
    }

    // Event loss estimation
    auto impactPar = mcCollision.impactParameter();
    auto mult05 = mcCollision.multMCNParticlesEta05();
    auto mult08 = mcCollision.multMCNParticlesEta08();
    hMC.fill(HIST("AllLoss/hImpactParameterGen"), impactPar);
    hMC.fill(HIST("AllLoss/hMultEta05Gen"), mult05);
    hMC.fill(HIST("AllLoss/hMultEta08Gen"), mult08);
    bool isSelectedEvent = false;
    float centrality = -1.f;

    for (auto const& collision : recCollisions) {

      if (!selectionEvent(collision, false))
        continue;

      if (selectCentEstimator == kFT0M) {
        centrality = collision.centFT0M();
      } else if (selectCentEstimator == kFT0A) {
        centrality = collision.centFT0A();
      } else if (selectCentEstimator == kFT0C) {
        centrality = collision.centFT0C();
      } else if (selectCentEstimator == kFV0A) {
        centrality = collision.centFV0A();
      } else {
        centrality = collision.centFT0M(); // default
      }

      isSelectedEvent = true;
    }

    hMC.fill(HIST("AllLoss/hImpactParvsCentr"), centrality, impactPar);
    hMC.fill(HIST("AllLoss/hMultEta05vsCentr"), centrality, mult05);
    hMC.fill(HIST("AllLoss/hMultEta08vsCentr"), centrality, mult08);

    if (isSelectedEvent) {
      hMC.fill(HIST("AllLoss/hImpactParameterRec"), impactPar);
      hMC.fill(HIST("AllLoss/hImpactParvsCentrRec"), centrality, impactPar);
      hMC.fill(HIST("AllLoss/hMultEta05Rec"), mult05);
      hMC.fill(HIST("AllLoss/hMultEta05vsCentrRec"), centrality, mult05);
      hMC.fill(HIST("AllLoss/hMultEta08Rec"), mult08);
      hMC.fill(HIST("AllLoss/hMultEta08vsCentrRec"), centrality, mult08);
    }

    // Generated MC
    for (const auto& mcPart : mcParticles) {

      if ((mcPart.y() < selectionConfig.motherRapidityMin || mcPart.y() > selectionConfig.motherRapidityMax) || std::abs(mcPart.pdgCode()) != o2::constants::physics::kK0Star892)
        continue;

      const float pt = mcPart.pt();

      // signal loss estimation
      hMC.fill(HIST("AllLoss/hKstarpTGenVsImpactParBeforeEvtSel"), pt, impactPar);
      hMC.fill(HIST("AllLoss/hKstarpTGenVsMultEta05BeforeEvtSel"), pt, mult05);
      hMC.fill(HIST("AllLoss/hKstarpTGenVsMultEta08BeforeEvtSel"), pt, mult08);
      if (isSelectedEvent) {
        // signal loss estimation
        hMC.fill(HIST("AllLoss/hKstarpTGenVsImpactParAfterEvtSel"), pt, impactPar);
        hMC.fill(HIST("AllLoss/hKstarpTGenVsMultEta05AfterEvtSel"), pt, mult05);
        hMC.fill(HIST("AllLoss/hKstarpTGenVsMultEta08AfterEvtSel"), pt, mult08);
      }
    } // end loop on gen particles
  }
  PROCESS_SWITCH(Kstar892LightIon, processAllLossMC, "Process All Signal Loss, Event Loss", false);

  void processRecMisID(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&, EventMCGenerated const&)
  {
    if (!collision.has_mcCollision()) {
      return;
    }

    if (selectCentEstimator == kFT0M) {
      centrality = collision.centFT0M();
    } else if (selectCentEstimator == kFT0A) {
      centrality = collision.centFT0A();
    } else if (selectCentEstimator == kFT0C) {
      centrality = collision.centFT0C();
    } else if (selectCentEstimator == kFV0A) {
      centrality = collision.centFV0A();
    } else {
      centrality = collision.centFT0M(); // default
    }

    if (!selectionEvent(collision, false)) {
      return;
    }

    for (const auto& [track1, track2] : combinations(CombinationsFullIndexPolicy(tracks, tracks))) {

      if (!selectionTrack(track1) || !selectionTrack(track2))
        continue;

      if (track1.index() >= track2.index())
        continue;

      if (track1.sign() * track2.sign() >= 0)
        continue;

      if (!track1.has_mcParticle() || !track2.has_mcParticle())
        continue;

      const auto mc1 = track1.mcParticle();
      const auto mc2 = track2.mcParticle();

      if (!mc1.isPhysicalPrimary() || !mc2.isPhysicalPrimary())
        continue;

      int pdg1 = std::abs(mc1.pdgCode());
      int pdg2 = std::abs(mc2.pdgCode());

      bool ok1 = (pdg1 == PDG_t::kPiPlus || pdg1 == PDG_t::kKPlus);
      bool ok2 = (pdg2 == PDG_t::kPiPlus || pdg2 == PDG_t::kKPlus);
      if (!ok1 || !ok2)
        continue;

      // pi-pi misidentification
      if (pdg1 == PDG_t::kPiPlus && pdg2 == PDG_t::kPiPlus) {
        ROOT::Math::PxPyPzMVector p1Fake(track1.px(), track1.py(), track1.pz(), massKa);
        ROOT::Math::PxPyPzMVector p2True(track2.px(), track2.py(), track2.pz(), massPi);

        auto misIDMother = p1Fake + p2True;
        if (misIDMother.Rapidity() > selectionConfig.motherRapidityMin && misIDMother.Rapidity() < selectionConfig.motherRapidityMax) {
          hMC.fill(HIST("RecMisID/hMassMisID"), misIDMother.Pt(), centrality, misIDMother.M());
        }
      }

      // KK misidentification
      if (pdg1 == PDG_t::kKPlus && pdg2 == PDG_t::kKPlus) {
        ROOT::Math::PxPyPzMVector p1Fake(track1.px(), track1.py(), track1.pz(), massPi);
        ROOT::Math::PxPyPzMVector p2True(track2.px(), track2.py(), track2.pz(), massKa);

        auto misIDMother = p1Fake + p2True;
        if (misIDMother.Rapidity() > selectionConfig.motherRapidityMin && misIDMother.Rapidity() < selectionConfig.motherRapidityMax) {
          hMC.fill(HIST("RecMisID/hMassMisID"), misIDMother.Pt(), centrality, misIDMother.M());
        }
      }
    }
  }
  PROCESS_SWITCH(Kstar892LightIon, processRecMisID, "Process Reconstructed MisID Background", false);

  void processRecReflection(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&, EventMCGenerated const&)
  {
    if (!collision.has_mcCollision())
      return;

    if (!selectionEvent(collision, false))
      return;

    if (selectCentEstimator == kFT0M) {
      centrality = collision.centFT0M();
    } else if (selectCentEstimator == kFT0A) {
      centrality = collision.centFT0A();
    } else if (selectCentEstimator == kFT0C) {
      centrality = collision.centFT0C();
    } else if (selectCentEstimator == kFV0A) {
      centrality = collision.centFV0A();
    } else {
      centrality = collision.centFT0M(); // default
    }

    for (const auto& [track1, track2] : combinations(CombinationsFullIndexPolicy(tracks, tracks))) {

      if (!selectionTrack(track1) || !selectionTrack(track2))
        continue;

      if (track1.index() >= track2.index())
        continue;

      if (track1.sign() * track2.sign() >= 0)
        continue;

      if (!track1.has_mcParticle() || !track2.has_mcParticle())
        continue;

      const auto mc1 = track1.mcParticle();
      const auto mc2 = track2.mcParticle();

      if (!mc1.isPhysicalPrimary() || !mc2.isPhysicalPrimary())
        continue;

      bool sameMother = false;
      int motherPDG = 0;

      for (const auto& m1 : mc1.mothers_as<aod::McParticles>()) {
        for (const auto& m2 : mc2.mothers_as<aod::McParticles>()) {
          if (m1.globalIndex() == m2.globalIndex()) {
            motherPDG = std::abs(m1.pdgCode());
            sameMother = true;
            break;
          }
        }
        if (sameMother)
          break;
      }

      if (!sameMother)
        continue;

      if (reflectionType == kRho) { // Rho0 (770) -> pi pi -> K pi
        if (motherPDG != PDG_t::kRho770_0)
          continue;

        if (std::abs(mc1.pdgCode()) != PDG_t::kPiPlus ||
            std::abs(mc2.pdgCode()) != PDG_t::kPiPlus)
          continue;

        // ---- permutation 1: track1 -> K
        ROOT::Math::PxPyPzMVector p1K(track1.px(), track1.py(), track1.pz(), massKa);
        ROOT::Math::PxPyPzMVector p2Pi(track2.px(), track2.py(), track2.pz(), massPi);

        auto fake1 = p1K + p2Pi;

        if (fake1.Rapidity() > selectionConfig.motherRapidityMin && fake1.Rapidity() < selectionConfig.motherRapidityMax)
          hMC.fill(HIST("Reflections/hReflection"), fake1.Pt(), centrality, fake1.M());

        // ---- permutation 2: track2 -> K
        ROOT::Math::PxPyPzMVector p1Pi(track1.px(), track1.py(), track1.pz(), massPi);
        ROOT::Math::PxPyPzMVector p2K(track2.px(), track2.py(), track2.pz(), massKa);

        auto fake2 = p1Pi + p2K;

        if (fake2.Rapidity() > selectionConfig.motherRapidityMin && fake2.Rapidity() < selectionConfig.motherRapidityMax)
          hMC.fill(HIST("Reflections/hReflection"), fake2.Pt(), centrality, fake2.M());

      } else if (reflectionType == kOmega) { // Omega (782) -> pi pi (pi0) -> K pi
        if (motherPDG != o2::constants::physics::kOmega)
          continue;

        if (std::abs(mc1.pdgCode()) != PDG_t::kPiPlus ||
            std::abs(mc2.pdgCode()) != PDG_t::kPiPlus)
          continue;

        // same two permutations as rho
        ROOT::Math::PxPyPzMVector p1K(track1.px(), track1.py(), track1.pz(), massKa);
        ROOT::Math::PxPyPzMVector p2Pi(track2.px(), track2.py(), track2.pz(), massPi);

        auto fake1 = p1K + p2Pi;

        if (fake1.Rapidity() > selectionConfig.motherRapidityMin && fake1.Rapidity() < selectionConfig.motherRapidityMax)
          hMC.fill(HIST("Reflections/hReflection"), fake1.Pt(), centrality, fake1.M());

        ROOT::Math::PxPyPzMVector p1Pi(track1.px(), track1.py(), track1.pz(), massPi);
        ROOT::Math::PxPyPzMVector p2K(track2.px(), track2.py(), track2.pz(), massKa);

        auto fake2 = p1Pi + p2K;

        if (fake2.Rapidity() > selectionConfig.motherRapidityMin && fake2.Rapidity() < selectionConfig.motherRapidityMax)
          hMC.fill(HIST("Reflections/hReflection"), fake2.Pt(), centrality, fake2.M());

      } else if (reflectionType == kPhi) { // Phi (1020) -> K K -> K pi
        if (motherPDG != o2::constants::physics::kPhi)
          continue;

        if (std::abs(mc1.pdgCode()) != PDG_t::kKPlus || std::abs(mc2.pdgCode()) != PDG_t::kKPlus)
          continue;

        // ---- permutation 1: track1 -> π
        ROOT::Math::PxPyPzMVector p1Pi(track1.px(), track1.py(), track1.pz(), massPi);
        ROOT::Math::PxPyPzMVector p2K(track2.px(), track2.py(), track2.pz(), massKa);

        auto fake1 = p1Pi + p2K;

        if (fake1.Rapidity() > selectionConfig.motherRapidityMin && fake1.Rapidity() < selectionConfig.motherRapidityMax)
          hMC.fill(HIST("Reflections/hReflection"), fake1.Pt(), centrality, fake1.M());

        // ---- permutation 2: track2 -> π
        ROOT::Math::PxPyPzMVector p1K(track1.px(), track1.py(), track1.pz(), massKa);
        ROOT::Math::PxPyPzMVector p2Pi(track2.px(), track2.py(), track2.pz(), massPi);

        auto fake2 = p1K + p2Pi;

        if (fake2.Rapidity() > selectionConfig.motherRapidityMin && fake2.Rapidity() < selectionConfig.motherRapidityMax)
          hMC.fill(HIST("Reflections/hReflection"), fake2.Pt(), centrality, fake2.M());
      } else if (reflectionType == kKstar) { //  K*0 (892) Self-Reflection (swap)

        if (motherPDG != o2::constants::physics::kK0Star892)
          continue;

        if (!((std::abs(mc1.pdgCode()) == PDG_t::kPiPlus && std::abs(mc2.pdgCode()) == PDG_t::kKPlus) || (std::abs(mc1.pdgCode()) == PDG_t::kKPlus && std::abs(mc2.pdgCode()) == PDG_t::kPiPlus)))
          continue;

        ROOT::Math::PxPyPzMVector p1Swap(track1.px(), track1.py(), track1.pz(), std::abs(mc1.pdgCode()) == PDG_t::kKPlus ? massPi : massKa);

        ROOT::Math::PxPyPzMVector p2Swap(track2.px(), track2.py(), track2.pz(), std::abs(mc2.pdgCode()) == PDG_t::kKPlus ? massPi : massKa);

        auto fake = p1Swap + p2Swap;

        if (fake.Rapidity() > selectionConfig.motherRapidityMin && fake.Rapidity() < selectionConfig.motherRapidityMax)
          hMC.fill(HIST("Reflections/hReflection"), fake.Pt(), centrality, fake.M());
      }
    }
  }
  PROCESS_SWITCH(Kstar892LightIon, processRecReflection, "Process particle reflection", false);

  Service<o2::framework::O2DatabasePDG> pdg;

  void processMCCheck(aod::McCollisions::iterator const& mccollision, soa::SmallGroups<EventCandidatesMC> const& collisions, aod::McParticles const& mcParticles, TrackCandidatesMC const&)
  {

    //---------------------------
    // Only INEL > 0 generated collisions
    // By counting number of primary charged particles in |eta| < 1
    //---------------------------
    int nChMC{0};
    int nChMCEta08{0};
    int nChFT0A{0};
    int nChFT0C{0};
    static constexpr float MinCharge{3.f};
    static constexpr float MinFT0A{3.5f};
    static constexpr float MaxFT0A{4.9f};
    static constexpr float MinFT0C{-3.3f};
    static constexpr float MaxFT0C{-2.1f};
    static constexpr float One{1.0f};
    static constexpr int ZeroInt{0};

    for (const auto& particle : mcParticles) {

      auto charge{0.};
      // Get the MC particle
      const auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (pdgParticle != nullptr) {
        charge = pdgParticle->Charge();
      } else {
        continue;
      }

      // Is it a charged particle?
      if (std::abs(charge) < MinCharge)
        continue;

      // Is it a primary particle?
      if (!particle.isPhysicalPrimary())
        continue;

      const float eta{particle.eta()};

      // TVX requirement
      if (eta > MinFT0A && eta < MaxFT0A) {
        nChFT0A++;
      }

      if (eta > MinFT0C && eta < MaxFT0C) {
        nChFT0C++;
      }

      if (std::abs(eta) < nchAcceptance) {
        nChMCEta08++;
      }

      // INEL > 0
      if (std::abs(eta) > One)
        continue;

      nChMC++;
    }

    //---------------------------
    // Only events with at least one charged particle in the FT0A and FT0C acceptances
    //---------------------------
    if (selectionConfig.selTVXMC) {
      if (!(nChFT0A > ZeroInt && nChFT0C > ZeroInt)) {
        hMC.fill(HIST("MCCheck/NchMCcentVsTVX"), nChMC, 0.5);
        return;
      }
      hMC.fill(HIST("MCCheck/NchMCcentVsTVX"), nChMC, 1.5);
    }

    //---------------------------
    // Only MC events with |Vtx Z| < 10 cm
    //---------------------------
    if (selectionConfig.isZvtxPosSelMC && (std::fabs(mccollision.posZ()) > selectionConfig.cfgVrtxZCut)) {
      return;
    }

    //---------------------------
    // Only INEL > 0 generated events
    //---------------------------
    if (selectionConfig.selINELgt0) {
      if (!(nChMC > ZeroInt)) {
        return;
      }
    }

    const auto& nRecColls{collisions.size()};
    hMC.fill(HIST("MCCheck/NumberOfRecoCollisions"), nRecColls);

    //---------------------------
    // Only Generated evets with at least one reconstrued collision
    //---------------------------
    if (nRecColls > ZeroInt) {

      // Finds the collisions with the largest number of contributors
      // in case nRecColls is larger than One
      int biggestNContribs{-1};
      int bestCollisionIndex{-1};
      centrality = -1.f;
      for (const auto& collision : collisions) {

        if (selectCentEstimator == kFT0M) {
          centrality = collision.centFT0M();
        } else if (selectCentEstimator == kFT0A) {
          centrality = collision.centFT0A();
        } else if (selectCentEstimator == kFT0C) {
          centrality = collision.centFT0C();
        } else if (selectCentEstimator == kFV0A) {
          centrality = collision.centFV0A();
        } else {
          centrality = collision.centFT0M(); // default
        }

        if (selectionConfig.selHasFT0 && !collision.has_foundFT0()) {
          continue;
        }

        if (biggestNContribs < collision.numContrib()) {
          biggestNContribs = collision.numContrib();
          bestCollisionIndex = collision.globalIndex();
        }

        // Needed to calculate denominator of the Event Splitting correction
        if (selectionEvent(collision, false)) {
          hMC.fill(HIST("MCCheck/Centrality_AllRecoEvt"), centrality);
        }
      }

      //---------------------------
      // Loop over the reconstructed collisions
      // Only that one with the largest number of contributors is considered
      //---------------------------
      centrality = -1.f;
      for (const auto& collision : collisions) {

        if (selectCentEstimator == kFT0M) {
          centrality = collision.centFT0M();
        } else if (selectCentEstimator == kFT0A) {
          centrality = collision.centFT0A();
        } else if (selectCentEstimator == kFT0C) {
          centrality = collision.centFT0C();
        } else if (selectCentEstimator == kFV0A) {
          centrality = collision.centFV0A();
        } else {
          centrality = collision.centFT0M(); // default
        }

        //---------------------------
        // Reject collisions if has_foundFT0() returns false
        //---------------------------
        if (selectionConfig.selHasFT0 && !collision.has_foundFT0()) {
          hMC.fill(HIST("MCCheck/CentVsFoundFT0"), centrality, 0.5);
          continue;
        }
        hMC.fill(HIST("MCCheck/CentVsFoundFT0"), centrality, 1.5);

        //---------------------------
        // Pick the collisions with the largest number of contributors
        //---------------------------
        if (bestCollisionIndex != collision.globalIndex()) {
          continue;
        }

        //---------------------------
        // Needed to construct the correlation between MC Nch v.s. centrality
        //---------------------------

        hMC.fill(HIST("MCCheck/Centrality_WRecoEvt"), centrality);
        hMC.fill(HIST("MCCheck/zPosMC"), mccollision.posZ());

        //---------------------------
        // Event selection
        // for reconstructed collisions
        //---------------------------
        if (!selectionEvent(collision, false)) {
          continue;
        }

        hMC.fill(HIST("MCCheck/Centrality_WRecoEvtWSelCri"), centrality);
        hMC.fill(HIST("MCCheck/NchMCVsCent"), centrality, nChMCEta08);
        hMC.fill(HIST("MCCheck/NchMC_WithRecoEvt"), nChMCEta08); // Numerator of event loss correction
        hMC.fill(HIST("MCCheck/zPos"), collision.posZ());
        hMC.fill(HIST("MCCheck/Cent"), centrality);

        //---------------------------
        // All Generated events with at least one associated reconstructed collision
        // The Generated events are not subjected to any selection criteria
        // However, the associated reconstructed collisions pass the selection criteria
        // This histograms are used for the denominator of the tracking efficiency
        //---------------------------
        for (const auto& mcPart : mcParticles) {
          if ((mcPart.y() < selectionConfig.motherRapidityMin || mcPart.y() > selectionConfig.motherRapidityMax) || std::abs(mcPart.pdgCode()) != o2::constants::physics::kK0Star892)
            continue;

          hMC.fill(HIST("MCCheck/PtKstarVsCentMC_WithRecoEvt"), mcPart.pt(), centrality);
          hMC.fill(HIST("MCCheck/PtKstarVsNchMC_WithRecoEvt"), mcPart.pt(), nChMCEta08); // Numerator of signal loss
        } // Loop over generated particles per generated collision
        // hMC.fill(HIST("MCCheck/NchVsCent"), centrality, nCh);
      } // Loop over Reco. Collisions: Only the collisions with the largest number of contributors
    } // If condition: Only simulated evets with at least one reconstrued collision

    //---------------------------
    // All Generated events irrespective of whether there is an associated reconstructed collision
    // Consequently, the centrality being a reconstructed quantity, might not always be available
    // Therefore it is expressed as a function of the generated pT and the generated Nch in ∣eta∣ < 0.8
    // This is used for the denominator of the signal loss correction
    //---------------------------
    for (const auto& mcPart : mcParticles) {
      if ((mcPart.y() < selectionConfig.motherRapidityMin || mcPart.y() > selectionConfig.motherRapidityMax) || std::abs(mcPart.pdgCode()) != o2::constants::physics::kK0Star892)
        continue;

      hMC.fill(HIST("MCCheck/PtKstarVsNchMC_AllGen"), mcPart.pt(), nChMCEta08);
    } // Loop over Generated Particles

    //---------------------------
    //  This is used for the denominator of the event loss correction
    //---------------------------
    hMC.fill(HIST("MCCheck/NchMC_AllGen"), nChMCEta08);
  }
  PROCESS_SWITCH(Kstar892LightIon, processMCCheck, "Cross-check MC analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Kstar892LightIon>(cfgc)};
}
