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

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <THn.h>
#include <TPDGCode.h>
#include <TRandom.h>
#include <TString.h>

#include <algorithm>
// #include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
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

    Configurable<bool> isApplyhasFT0{"isApplyhasFT0", false, "Apply has_foundFT0 event selection"};

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
    Configurable<float> cfgITSChi2NCl{"cfgITSChi2NCl", 36.0, "ITS Chi2/NCl"};
    Configurable<float> cfgTPCChi2NClMax{"cfgTPCChi2NClMax", 4.0, "TPC Chi2/NCl"};
    Configurable<float> cfgTPCChi2NClMin{"cfgTPCChi2NClMin", 0.0, "TPC Chi2/NCl"};
    // Configurable<bool> isUseITSTPCRefit{"isUseITSTPCRefit", false, "Require ITS Refit"};
    Configurable<bool> isApplyPtDepDCAxyCut{"isApplyPtDepDCAxyCut", false, "Apply pT dependent DCAxy cut"};
    Configurable<bool> isGoldenChi2{"isGoldenChi2", false, "Apply golden chi2 cut"};
    Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
    // Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", false, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
    Configurable<float> cfgTOFBetaCut{"cfgTOFBetaCut", 0.0, "cut TOF beta"};

    Configurable<bool> isAvoidsplitrackMC{"isAvoidsplitrackMC", true, "avoid split track in MC"};

    // cuts on mother
    // Configurable<bool> isApplyCutsOnMother{"isApplyCutsOnMother", false, "Enable additional cuts on Kstar mother"};
    // Configurable<float> cMaxPtMotherCut{"cMaxPtMotherCut", 15.0, "Maximum pt of mother cut"};
    // Configurable<float> cMaxMinvMotherCut{"cMaxMinvMotherCut", 1.5, "Maximum mass of mother cut"};
    Configurable<float> motherRapidityMax{"motherRapidityMax", 0.5, "Maximum rapidity of mother"};
    Configurable<float> motherRapidityMin{"motherRapidityMin", -0.5, "Minimum rapidity of mother"};

    // PID selections
    Configurable<int> pidStrategy{"pidStrategy", 0, "0=Standard, 1=pTDependent, 2=pTDependentTOF, 3=ThreePtDependent, 4=PIDCompare"};
    Configurable<int> pidMode{"pidMode", 0, "0=Combined,1=TPC,2=TOF,3=TOFHIT,4=TOFVeto"};
    Configurable<int> misIdStrategy{"misIdStrategy", -1, "-1=Disabled, 0=Standard, 1=PtDependent, 2=PtDependentCompare"};

    Configurable<float> nsigmaCutTPCPi{"nsigmaCutTPCPi", 3.0, "TPC Nsigma cut for pions"};
    Configurable<float> nsigmaCutTPCKa{"nsigmaCutTPCKa", 3.0, "TPC Nsigma cut for kaons"};
    Configurable<float> nsigmaCutTOFPi{"nsigmaCutTOFPi", 3.0, "TOF Nsigma cut for pions"};
    Configurable<float> nsigmaCutTOFKa{"nsigmaCutTOFKa", 3.0, "TOF Nsigma cut for kaons"};
    Configurable<float> nsigmaCutCombinedKa{"nsigmaCutCombinedKa", 3.0, "Combined Nsigma cut for kaon"};
    Configurable<float> nsigmaCutCombinedPi{"nsigmaCutCombinedPi", 3.0, "Combined Nsigma cut for pion"};

    Configurable<float> nsigmaCutTpcMisId{"nsigmaCutTpcMisId", 1.0, "MID Nsigma cut for pion and kaon in TPC"};

    Configurable<bool> isApplyFakeTrack{"isApplyFakeTrack", false, "Fake track selection"};
    Configurable<float> cfgFakeTrackCutKa{"cfgFakeTrackCutKa", 0.3, "Cut based on momentum difference in global and TPC tracks for kaons"};
    Configurable<float> cfgFakeTrackCutPi{"cfgFakeTrackCutPi", 0.3, "Cut based on momentum difference in global and TPC tracks for pions"};

    Configurable<float> pionMisIdPtLow{"pionMisIdPtLow", 1.0, "Low pT cut for pion in MID"};
    Configurable<float> pionMisIdPtHigh{"pionMisIdPtHigh", 2.5, "High pT cut for pion in MID"};
    Configurable<float> kaonMisIdPtLow{"kaonMisIdPtLow", 0.7, "Low pT cut for kaon in MID"};
    Configurable<float> kaonMisIdPtHigh{"kaonMisIdPtHigh", 2.5, "High pT cut for kaon in MID"};

    Configurable<float> lowPtCutPid{"lowPtCutPid", 0.5, "Low pT cut for PID"};
    Configurable<float> highPtCutPid{"highPtCutPid", 6.0, "High pT cut for PID"};

    Configurable<float> lowPtPid{"lowPtPid", 2.5, "Low pT cut for PID"};
    Configurable<float> midPtPid{"midPtPid", 1.5, "Mid pT cut for PID"};
    Configurable<float> highPtPid{"highPtPid", 2.5, "High pT cut for PID"};

    Configurable<bool> selHasFT0MC{"selHasFT0MC", true, "Has FT0?"};
    Configurable<bool> isZvtxPosSelMC{"isZvtxPosSelMC", true, "Zvtx position selection for MC events?"};
    Configurable<bool> selTVXMC{"selTVXMC", true, "apply TVX selection in MC?"};
    Configurable<bool> selINELgt0MC{"selINELgt0MC", true, "Select INEL > 0?"};
  } selectionConfig;

  Configurable<bool> calcLikeSign{"calcLikeSign", true, "Calculate Like Sign"};
  Configurable<bool> calcRotational{"calcRotational", true, "Calculate Rotational"};
  Configurable<int> cRotations{"cRotations", 3, "Number of random rotations in the rotational background"};
  Configurable<int> rotationalCut{"rotationalCut", 10, "Cut value (Rotation angle pi - pi/cut and pi + pi/cut)"};

  // Confugrable for QA histograms
  Configurable<bool> cQAplots{"cQAplots", true, "cQAplots"};
  Configurable<bool> additionalKin{"additionalKin", false, "Additional kinematics histograms for processMisIdKinematics"};
  Configurable<bool> cQAevents{"cQAevents", true, "centrality dist, Vz, Event Cut hists"};

  Configurable<int> selectCentEstimator{"selectCentEstimator", 0, "Select centrality estimator: 0 - FT0M, 1 - FT0A, 2 - FT0C, 3 - FV0A"};

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

  enum class PIDParticle {
    kPion,
    kKaon
  };

  enum class PIDStrategy {
    Standard,
    PtDependent,
    PtDependentTOF,
    ThreePtDependent,
    PIDCompare
  };

  enum class PIDMode {
    Combined,
    TPC,
    TOF,
    TOFHIT,
    TOFVeto
  };

  enum class MIDStrategy {
    Disabled = -1,
    Standard = 0,
    PtDependent = 1,
    PtDependentCompare = 2
  };

  int noOfDaughters = 2;
  int initialValue = -9999999;
  double massPi = o2::constants::physics::MassPiPlus;
  double massKa = o2::constants::physics::MassKPlus;

  double pionPidPtLow = 1.0, pionPidPtHigh = 2.5, kaonPidPtLow = 0.7, kaonPidPtHigh = 2.5;

  int kKstar14300 = 10311, kK11400Plus = 20323, kKStar1680Plus = 30323, kK21770Plus = 10325, kK21820Plus = 20325, kKstar14100 = 100313, kKstar16800 = 30313, kK218200 = 20315, kK217700 = 10315, kKstar214300 = 315, kKstar1430Plus = 10321;

  // ThreePtDependent PID cuts
  // Regions:
  // 0 : pT < lowPtCutPid
  // 1 : lowPtCutPid <= pT < highPtCutPID
  // 2 : pT >= highPtCutPID
  // static constexpr std::array<float, 3> TPCPiCuts3Pt{3.0f, 1.5f, 2.0f};
  // static constexpr std::array<float, 3> TOFPiCuts3Pt{3.0f, 2.0f, 3.0f};

  // static constexpr std::array<float, 3> TPCKaCuts3Pt{2.0f, 1.5f, 2.0f};
  // static constexpr std::array<float, 3> TOFKaCuts3Pt{3.0f, 2.0f, 3.0f};

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
      std::string("INEL > 0") + check(selectionConfig.isApplyINELgt0.value),
      std::string("hasFT0") + check(selectionConfig.isApplyhasFT0.value)};

    // assign labels
    for (size_t i = 0; i < eveCutLabels.size(); ++i) {
      hEventSelection.get<TH1>(HIST("hEventCut"))->GetXaxis()->SetBinLabel(i + 1, eveCutLabels[i].c_str());
    }

    // for primary tracksbinsCentPlot
    if (cQAplots) {
      hOthers.add("dE_by_dx_TPC", "dE/dx signal in the TPC as a function of pT", kTH2F, {axisPtfordEbydx, axisdEdx});
      hOthers.add("hEta_after", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
      hOthers.add("hTOFBetaKa", "Beta distribution from TOF", kTH1F, {{50, 0.0f, 1.0f}});
      hOthers.add("hTOFBetaPi", "Beta distribution from TOF", kTH1F, {{50, 0.0f, 1.0f}});
      hOthers.add("hTPCpDiffKa", "Momentum difference between global and TPC tracks for Kaons", kTH2F, {{100, -1.0f, 1.0f}, ptAxis});
      hOthers.add("hTPCpDiffPi", "Momentum difference between global and TPC tracks for Kaons", kTH2F, {{100, -1.0f, 1.0f}, ptAxis});

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

      hPID.add("Before/hTPCnsigKa_Neg_mult_p", "TPC nsigma of K^{-} before PID with p and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("Before/hTPCnsigPi_Neg_mult_p", "TPC nsigma of #pi^{-} before PID with p and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("Before/hTOFnsigKa_Neg_mult_p", "TOF nsigma of K^{-} before PID with p and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("Before/hTOFnsigPi_Neg_mult_p", "TOF nsigma of #pi^{-} before PID with p and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});

      hPID.add("Before/hTPCnsigKa_Pos_mult_pt", "TPC nsigma of K^{+} before PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("Before/hTPCnsigPi_Pos_mult_pt", "TPC nsigma of #pi^{+} before PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("Before/hTOFnsigKa_Pos_mult_pt", "TOF nsigma of K^{+} before PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("Before/hTOFnsigPi_Pos_mult_pt", "TOF nsigma of #pi^{+} before PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});

      hPID.add("Before/hTPCnsigKa_Pos_mult_p", "TPC nsigma of K^{-} before PID with p and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("Before/hTPCnsigPi_Pos_mult_p", "TPC nsigma of #pi^{-} before PID with p and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("Before/hTOFnsigKa_Pos_mult_p", "TOF nsigma of K^{-} before PID with p and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("Before/hTOFnsigPi_Pos_mult_p", "TOF nsigma of #pi^{-} before PID with p and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});

      hPID.add("After/hNsigma_TPC_TOF_Ka_pt", "N #sigma Kaon TPC TOF after", kTH3F, {{50, -5.0f, 5.0f}, {50, -5.0f, 5.0f}, ptAxis});
      hPID.add("After/hNsigma_TPC_TOF_Pi_pt", "N #sigma Pion TPC TOF after", kTH3F, {{50, -5.0f, 5.0f}, {50, -5.0f, 5.0f}, ptAxis});

      hPID.add("After/hTPCnsigKa_Neg_mult_pt", "TPC nsigma of K^{-} after PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("After/hTPCnsigPi_Neg_mult_pt", "TPC nsigma of #pi^{-} after PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("After/hTOFnsigKa_Neg_mult_pt", "TOF nsigma of K^{-} after PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("After/hTOFnsigPi_Neg_mult_pt", "TOF nsigma of #pi^{-} after PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});

      hPID.add("After/hTPCnsigKa_Neg_mult_p", "TPC nsigma of K^{-} After PID with p and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("After/hTPCnsigPi_Neg_mult_p", "TPC nsigma of #pi^{-} After PID with p and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("After/hTOFnsigKa_Neg_mult_p", "TOF nsigma of K^{-} After PID with p and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("After/hTOFnsigPi_Neg_mult_p", "TOF nsigma of #pi^{-} After PID with p and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});

      hPID.add("After/hTPCnsigKa_Pos_mult_pt", "TPC nsigma of K^{+} after PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("After/hTPCnsigPi_Pos_mult_pt", "TPC nsigma of #pi^{+} after PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("After/hTOFnsigKa_Pos_mult_pt", "TOF nsigma of K^{+} after PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("After/hTOFnsigPi_Pos_mult_pt", "TOF nsigma of #pi^{+} after PID with pt and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});

      hPID.add("After/hTPCnsigKa_Pos_mult_p", "TPC nsigma of K^{-} After PID with p and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("After/hTPCnsigPi_Pos_mult_p", "TPC nsigma of #pi^{-} After PID with p and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("After/hTOFnsigKa_Pos_mult_p", "TOF nsigma of K^{-} After PID with p and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
      hPID.add("After/hTOFnsigPi_Pos_mult_p", "TOF nsigma of #pi^{-} After PID with p and centrality", kTH3F, {{100, -10.0f, 10.0f}, centralityAxis, ptAxis});
    }

    // KStar histograms
    hInvMass.add("h3KstarInvMassUnlikeSign", "kstar Unlike Sign", kTH3F, {centralityAxis, ptAxis, invmassAxis});
    hInvMass.add("h3KstarInvMassMixed", "kstar Mixed", kTH3F, {centralityAxis, ptAxis, invmassAxis});
    if (calcLikeSign) {
      hInvMass.add("h3KstarInvMasslikeSignPP", "kstar like Sign", kTH3F, {centralityAxis, ptAxis, invmassAxis});
      hInvMass.add("h3KstarInvMasslikeSignMM", "kstar like Sign", kTH3F, {centralityAxis, ptAxis, invmassAxis});
    }
    if (calcRotational) {
      hInvMass.add("h3KstarInvMassRotated", "kstar rotated", kTH3F, {centralityAxis, ptAxis, invmassAxis});
    }

    // MC histograms
    if (doprocessGen) {
      hMC.add("Gen/hGenNo", "MC Event statistics", kTH1F, {{10, 0.0f, 10.0f}});
      hMC.add("Gen/h3KstarPtCentMassDirect", "pT distribution of generated K*0 with centrality and mass (direct from generator)", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("Gen/h3KstarPtCentMassRec", "pT distribution of generated K*0 with centrality and mass (reconstructed from daughter tracks)", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("Gen/h1GenCent", "centrality generated", kTH1F, {centralityAxis});
      hMC.add("Gen/hAllGenCollisions", "All generated events", kTH1F, {centralityAxis});
      hMC.add("Gen/hAllGenCollisions1Rec", "All gen events with at least one rec event", kTH1F, {centralityAxis});
      hMC.add("Gen/hAllKstarGenCollisisons", "All generated Kstar in events with rapidity in 0.5", kTH2F, {ptAxis, centralityAxis});
      hMC.add("Gen/hAllKstarGenCollisisons1Rec", "All generated Kstar in events with at least one rec event in rapidity in 0.5", kTH2F, {ptAxis, centralityAxis});
    }

    if (doprocessRec) {
      hMC.add("Rec/hAllRecCollisions", "All reconstructed events", kTH1F, {centralityAxis});
      hMC.add("Rec/h3KstarPtCentMassrec", "pT of reconstructed K*0 with centrality and mass (reconstructed from rec daughter tracks)", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("Rec/h3KstarPtCentMassgen", "pT of reconstructed K*0 with centrality and mass (reconstructed from gen daughter tracks)", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("Rec/h1RecCent", "centrality reconstructed", kTH1F, {centralityAxis});
      hMC.add("Rec/h1KSRecsplit", "K*0 Rec split", kTH1F, {{100, 0.0f, 10.0f}});

      hMC.add("Rec/hMassShift", "#Delta M = m_{rec} - m_{gen}; #it{p}_{T}_{gen}; #it{p}_{T}_{rec}; #Delta M", kTH3F, {ptAxis, ptAxis, {2000, -0.1, 0.1}});
    }

    // Signal Loss & Event Loss
    if (doprocessEvtLossSigLossMC) {
      hMC.add("ImpactCorr/hImpactParameterGen", "Impact parameter of generated MC events", kTH1F, {impactParAxis});
      hMC.add("ImpactCorr/hImpactParameterRec", "Impact parameter of selected MC events", kTH1F, {impactParAxis});
      hMC.add("ImpactCorr/hImpactParvsCentrRec", "Impact parameter of selected MC events vs centrality", kTH2F, {centralityAxis, impactParAxis});
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
      hMC.add("RecMisID/hMassMisIDPiPi", "Reconstruction misidentification", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("RecMisID/hMassMisIDKK", "Reconstruction misidentification", kTH3F, {ptAxis, centralityAxis, invmassAxis});
    }

    if (doprocessRecReflection) {
      hMC.add("Reflections/hRhoToKpi", "Refelction template of Rho", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("Reflections/hOmegaToKpi", "Refelction template of Omega", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("Reflections/hPhiToKpi", "Refelction template of Phi", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("Reflections/hKstarSelf", "Refelction template of Kstar", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("Reflections/hEtaToKpi", "Refelction template of Eta", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("Reflections/hEtaPrimeToKpi", "Refelction template of Eta'", kTH3F, {ptAxis, centralityAxis, invmassAxis});

      hMC.add("FeedDown/hK1_1270", "Distribution for K1(1270)", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("FeedDown/hK1_1400", "Distribution for K1(1270)", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("FeedDown/hKstar1410", "Distribution for K*(1410)", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("FeedDown/hKstar0_1430_0", "Distribution for K*0(1430)", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("FeedDown/hKstar0_1430_ch", "Distribution for K*+-(1430)", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("FeedDown/hKstar2_1430", "Distribution for K*2(1430)", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("FeedDown/hK2_1770", "Distribution for K2(1770)", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("FeedDown/hK2_1820", "Distribution for K2(1820)", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("FeedDown/hKstar1680", "Distribution for K2(1680)", kTH3F, {ptAxis, centralityAxis, invmassAxis});
    }

    if (doprocessRecCorrelatedBackground) {
      hMC.add("CorrelatedBG/hKstar1410", "Wrong pair distribution for K*(1410)", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("CorrelatedBG/hKstar0_1430", "Wrong pair distribution for K*(1410)", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("CorrelatedBG/hK1_1270", "Wrong pair distribution for K*(1410)", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("CorrelatedBG/hK1_1400", "Wrong pair distribution for K*(1410)", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("CorrelatedBG/hKstar1680", "Wrong pair distribution for K*(1410)", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("CorrelatedBG/hK2_1770", "Wrong pair distribution for K*(1410)", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("CorrelatedBG/hK2_1820", "Wrong pair distribution for K*(1410)", kTH3F, {ptAxis, centralityAxis, invmassAxis});
    }
    if (doprocessTemplateMC) {
      hMC.add("Template/hSignal", "True K*0 signal", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("Template/hKstarReflection", "K*0 signal due to mis-identification", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("Template/hSameMotherOther", "Kpi pair from same mother other than K*0", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("Template/hDifferentMother", "Kpi pair from different mothers", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("Template/hPrimaryResonance", "Kpi pair one primary another from decay", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("Template/hPrimaryPrimary", "Primary Kpi pair", kTH3F, {ptAxis, centralityAxis, invmassAxis});
    }

    if (doprocessRecKinematics) {
      hMC.add("Kinematics/h1RecCent", "centrality reconstructed", kTH1F, {centralityAxis});

      hMC.add("Kinematics/hPDGMotherKstar", "Mother PDG of Ka and Pi from K*0", kTH2F, {{100000, -50000.0, 50000.0}, invmassAxis});
      hMC.add("Kinematics/hOpenAngleKstar", "Opening angle vs M(K#pi);Opening angle (rad);M(K#pi) (GeV/c^{2})", kTH2F, {{180, 0., o2::constants::math::PI}, invmassAxis});
      hMC.add("Kinematics/hDeltaPhiKstar", "#Delta#phi vs M(K#pi);#Delta#phi (rad);M(K#pi) (GeV/c^{2})", kTH2F, {{180, -o2::constants::math::PI, o2::constants::math::PI}, invmassAxis});
      hMC.add("Kinematics/hDeltaEtaKstar", "#Delta#eta vs M(K#pi);#Delta#eta;M(K#pi) (GeV/c^{2})", kTH2F, {{200, -2.0, 2.0}, invmassAxis});
      hMC.add("Kinematics/hDeltaRKstar", "#DeltaR vs M(K#pi);#DeltaR;M(K#pi) (GeV/c^{2})", kTH2F, {{200, 0.0, 5.0}, invmassAxis});
      hMC.add("Kinematics/hKaonPtKstar", "Kaon p_{T} vs M(K#pi);p_{T}^{K} (GeV/c);M(K#pi) (GeV/c^{2})", kTH2F, {ptAxis, invmassAxis});
      hMC.add("Kinematics/hPionPtKstar", "Pion p_{T} vs M(K#pi);p_{T}^{#pi} (GeV/c);M(K#pi) (GeV/c^{2})", kTH2F, {ptAxis, invmassAxis});

      hMC.add("Kinematics/hPtCentMassKstar", "p_{T} vs Centrality vs M(K#pi);p_{T} (GeV/c);Centrality (%);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, centralityAxis, invmassAxis});

      hMC.add("Kinematics/hOpenAngleOther", "Opening angle vs M(K#pi) (non-K*);Opening angle (rad);M(K#pi) (GeV/c^{2})", kTH2F, {{180, 0., o2::constants::math::PI}, invmassAxis});
      hMC.add("Kinematics/hDeltaPhiOther", "#Delta#phi vs M(K#pi) (non-K*);#Delta#phi (rad);M(K#pi) (GeV/c^{2})", kTH2F, {{180, -o2::constants::math::PI, o2::constants::math::PI}, invmassAxis});
      hMC.add("Kinematics/hDeltaEtaOther", "#Delta#eta vs M(K#pi) (non-K*);#Delta#eta;M(K#pi) (GeV/c^{2})", kTH2F, {{200, -2.0, 2.0}, invmassAxis});
      hMC.add("Kinematics/hDeltaROther", "#DeltaR vs M(K#pi) (non-K*);#DeltaR;M(K#pi) (GeV/c^{2})", kTH2F, {{200, 0.0, 5.0}, invmassAxis});
      hMC.add("Kinematics/hKaonPtOther", "Kaon p_{T} vs M(K#pi) (non-K*);p_{T}^{K} (GeV/c);M(K#pi) (GeV/c^{2})", kTH2F, {ptAxis, invmassAxis});
      hMC.add("Kinematics/hPionPtOther", "Pion p_{T} vs M(K#pi) (non-K*);p_{T}^{#pi} (GeV/c);M(K#pi) (GeV/c^{2})", kTH2F, {ptAxis, invmassAxis});
      hMC.add("Kinematics/hPtCentMassOther", "p_{T} vs Centrality vs M(K#pi) (non-K*);p_{T} (GeV/c);Centrality (%);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, centralityAxis, invmassAxis});
      hMC.add("Kinematics/hOtherMotherPDG", "Common mother PDG of non-K* pairs;PDG Code;M(K#pi) (GeV/c^{2})", kTH2F, {{100000, -50000.0, 50000.0}, invmassAxis});
    }

    if (doprocessMisIdKinematics) {
      hMC.add("KinematicsMisId/hRecCent", "centrality reconstructed", kTH1F, {centralityAxis});

      hMC.add("KinematicsMisId/hTrueKstar_PtDeltaRMass", "#DeltaR vs M(K#pi);#DeltaR;M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {200, 0.0, 5.0}, invmassAxis});
      hMC.add("KinematicsMisId/hTrueKstar_PtCentMass", "p_{T} vs Centrality vs M(K#pi);p_{T} (GeV/c);Centrality (%);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, centralityAxis, invmassAxis});

      hMC.add("KinematicsMisId/hOmega_PtDeltaRMass", "#DeltaR vs M(K#pi);#DeltaR;M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {200, 0.0, 5.0}, invmassAxis});
      hMC.add("KinematicsMisId/hOmega_PtCentMass", "p_{T} vs Centrality vs M(K#pi);p_{T} (GeV/c);Centrality (%);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, centralityAxis, invmassAxis});

      hMC.add("KinematicsMisId/hRho_PtDeltaRMass", "#DeltaR vs M(K#pi);#DeltaR;M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {200, 0.0, 5.0}, invmassAxis});
      hMC.add("KinematicsMisId/hRho_PtCentMass", "p_{T} vs Centrality vs M(K#pi);p_{T} (GeV/c);Centrality (%);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, centralityAxis, invmassAxis});

      hMC.add("KinematicsMisId/hEta_PtDeltaRMass", "#DeltaR vs M(K#pi);#DeltaR;M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {200, 0.0, 5.0}, invmassAxis});
      hMC.add("KinematicsMisId/hEta_PtCentMass", "p_{T} vs Centrality vs M(K#pi);p_{T} (GeV/c);Centrality (%);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, centralityAxis, invmassAxis});

      hMC.add("KinematicsMisId/hEtaP_PtDeltaRMass", "#DeltaR vs M(K#pi);#DeltaR;M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {200, 0.0, 5.0}, invmassAxis});
      hMC.add("KinematicsMisId/hEtaP_PtCentMass", "p_{T} vs Centrality vs M(K#pi);p_{T} (GeV/c);Centrality (%);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, centralityAxis, invmassAxis});

      hMC.add("KinematicsMisId/hPhi_PtDeltaRMass", "#DeltaR vs M(K#pi);#DeltaR;M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {200, 0.0, 5.0}, invmassAxis});
      hMC.add("KinematicsMisId/hPhi_PtCentMass", "p_{T} vs Centrality vs M(K#pi);p_{T} (GeV/c);Centrality (%);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, centralityAxis, invmassAxis});

      hMC.add("KinematicsMisId/hKstarMisId_PtDeltaRMass", "#DeltaR vs M(K#pi);#DeltaR;M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {200, 0.0, 5.0}, invmassAxis});
      hMC.add("KinematicsMisId/hKstarMisId_PtCentMass", "p_{T} vs Centrality vs M(K#pi);p_{T} (GeV/c);Centrality (%);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, centralityAxis, invmassAxis});

      hMC.add("KinematicsMisId/hOtherMisId_PtDeltaRMass", "#DeltaR vs M(K#pi);#DeltaR;M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {200, 0.0, 5.0}, invmassAxis});
      hMC.add("KinematicsMisId/hOtherMisId_PtCentMass", "p_{T} vs Centrality vs M(K#pi);p_{T} (GeV/c);Centrality (%);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, centralityAxis, invmassAxis});

      if (additionalKin) {
        hMC.add("KinematicsMisId/hTrueKstar_PtOpenAngleMass", "Opening angle vs M(K#pi);Opening angle (rad);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {180, 0., o2::constants::math::PI}, invmassAxis});
        hMC.add("KinematicsMisId/hTrueKstar_PtDeltaPhiMass", "#Delta#phi vs M(K#pi);#Delta#phi (rad);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {180, -o2::constants::math::PI, o2::constants::math::PI}, invmassAxis});
        hMC.add("KinematicsMisId/hTrueKstar_PtDeltaEtaMass", "#Delta#eta vs M(K#pi);#Delta#eta;M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {200, -2.0, 2.0}, invmassAxis});
        hMC.add("KinematicsMisId/hTrueKstar_PtKaonPtMass", "Kaon p_{T} vs M(K#pi);p_{T}^{K} (GeV/c);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, ptAxis, invmassAxis});
        hMC.add("KinematicsMisId/hTrueKstar_PtPionPtMass", "Pion p_{T} vs M(K#pi);p_{T}^{#pi} (GeV/c);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, ptAxis, invmassAxis});

        hMC.add("KinematicsMisId/hOmega_PtOpenAngleMass", "Opening angle vs M(K#pi);Opening angle (rad);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {180, 0., o2::constants::math::PI}, invmassAxis});
        hMC.add("KinematicsMisId/hOmega_PtDeltaPhiMass", "#Delta#phi vs M(K#pi);#Delta#phi (rad);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {180, -o2::constants::math::PI, o2::constants::math::PI}, invmassAxis});
        hMC.add("KinematicsMisId/hOmega_PtDeltaEtaMass", "#Delta#eta vs M(K#pi);#Delta#eta;M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {200, -2.0, 2.0}, invmassAxis});
        hMC.add("KinematicsMisId/hOmega_PtKaonPtMass", "Kaon p_{T} vs M(K#pi);p_{T}^{K} (GeV/c);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, ptAxis, invmassAxis});
        hMC.add("KinematicsMisId/hOmega_PtPionPtMass", "Pion p_{T} vs M(K#pi);p_{T}^{#pi} (GeV/c);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, ptAxis, invmassAxis});

        hMC.add("KinematicsMisId/hRho_PtOpenAngleMass", "Opening angle vs M(K#pi);Opening angle (rad);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {180, 0., o2::constants::math::PI}, invmassAxis});
        hMC.add("KinematicsMisId/hRho_PtDeltaPhiMass", "#Delta#phi vs M(K#pi);#Delta#phi (rad);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {180, -o2::constants::math::PI, o2::constants::math::PI}, invmassAxis});
        hMC.add("KinematicsMisId/hRho_PtDeltaEtaMass", "#Delta#eta vs M(K#pi);#Delta#eta;M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {200, -2.0, 2.0}, invmassAxis});
        hMC.add("KinematicsMisId/hRho_PtKaonPtMass", "Kaon p_{T} vs M(K#pi);p_{T}^{K} (GeV/c);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, ptAxis, invmassAxis});
        hMC.add("KinematicsMisId/hRho_PtPionPtMass", "Pion p_{T} vs M(K#pi);p_{T}^{#pi} (GeV/c);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, ptAxis, invmassAxis});

        hMC.add("KinematicsMisId/hEta_PtOpenAngleMass", "Opening angle vs M(K#pi);Opening angle (rad);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {180, 0., o2::constants::math::PI}, invmassAxis});
        hMC.add("KinematicsMisId/hEta_PtDeltaPhiMass", "#Delta#phi vs M(K#pi);#Delta#phi (rad);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {180, -o2::constants::math::PI, o2::constants::math::PI}, invmassAxis});
        hMC.add("KinematicsMisId/hEta_PtDeltaEtaMass", "#Delta#eta vs M(K#pi);#Delta#eta;M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {200, -2.0, 2.0}, invmassAxis});
        hMC.add("KinematicsMisId/hEta_PtKaonPtMass", "Kaon p_{T} vs M(K#pi);p_{T}^{K} (GeV/c);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, ptAxis, invmassAxis});
        hMC.add("KinematicsMisId/hEta_PtPionPtMass", "Pion p_{T} vs M(K#pi);p_{T}^{#pi} (GeV/c);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, ptAxis, invmassAxis});

        hMC.add("KinematicsMisId/hEtaP_PtOpenAngleMass", "Opening angle vs M(K#pi);Opening angle (rad);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {180, 0., o2::constants::math::PI}, invmassAxis});
        hMC.add("KinematicsMisId/hEtaP_PtDeltaPhiMass", "#Delta#phi vs M(K#pi);#Delta#phi (rad);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {180, -o2::constants::math::PI, o2::constants::math::PI}, invmassAxis});
        hMC.add("KinematicsMisId/hEtaP_PtDeltaEtaMass", "#Delta#eta vs M(K#pi);#Delta#eta;M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {200, -2.0, 2.0}, invmassAxis});
        hMC.add("KinematicsMisId/hEtaP_PtKaonPtMass", "Kaon p_{T} vs M(K#pi);p_{T}^{K} (GeV/c);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, ptAxis, invmassAxis});
        hMC.add("KinematicsMisId/hEtaP_PtPionPtMass", "Pion p_{T} vs M(K#pi);p_{T}^{#pi} (GeV/c);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, ptAxis, invmassAxis});

        hMC.add("KinematicsMisId/hPhi_PtOpenAngleMass", "Opening angle vs M(K#pi);Opening angle (rad);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {180, 0., o2::constants::math::PI}, invmassAxis});
        hMC.add("KinematicsMisId/hPhi_PtDeltaPhiMass", "#Delta#phi vs M(K#pi);#Delta#phi (rad);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {180, -o2::constants::math::PI, o2::constants::math::PI}, invmassAxis});
        hMC.add("KinematicsMisId/hPhi_PtDeltaEtaMass", "#Delta#eta vs M(K#pi);#Delta#eta;M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {200, -2.0, 2.0}, invmassAxis});
        hMC.add("KinematicsMisId/hPhi_PtKaonPtMass", "Kaon p_{T} vs M(K#pi);p_{T}^{K} (GeV/c);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, ptAxis, invmassAxis});
        hMC.add("KinematicsMisId/hPhi_PtPionPtMass", "Pion p_{T} vs M(K#pi);p_{T}^{#pi} (GeV/c);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, ptAxis, invmassAxis});

        hMC.add("KinematicsMisId/hKstarMisId_PtOpenAngleMass", "Opening angle vs M(K#pi);Opening angle (rad);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {180, 0., o2::constants::math::PI}, invmassAxis});
        hMC.add("KinematicsMisId/hKstarMisId_PtDeltaPhiMass", "#Delta#phi vs M(K#pi);#Delta#phi (rad);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {180, -o2::constants::math::PI, o2::constants::math::PI}, invmassAxis});
        hMC.add("KinematicsMisId/hKstarMisId_PtDeltaEtaMass", "#Delta#eta vs M(K#pi);#Delta#eta;M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {200, -2.0, 2.0}, invmassAxis});
        hMC.add("KinematicsMisId/hKstarMisId_PtKaonPtMass", "Kaon p_{T} vs M(K#pi);p_{T}^{K} (GeV/c);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, ptAxis, invmassAxis});
        hMC.add("KinematicsMisId/hKstarMisId_PtPionPtMass", "Pion p_{T} vs M(K#pi);p_{T}^{#pi} (GeV/c);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, ptAxis, invmassAxis});

        hMC.add("KinematicsMisId/hOtherMisId_PtOpenAngleMass", "Opening angle vs M(K#pi);Opening angle (rad);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {180, 0., o2::constants::math::PI}, invmassAxis});
        hMC.add("KinematicsMisId/hOtherMisId_PtDeltaPhiMass", "#Delta#phi vs M(K#pi);#Delta#phi (rad);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {180, -o2::constants::math::PI, o2::constants::math::PI}, invmassAxis});
        hMC.add("KinematicsMisId/hOtherMisId_PtDeltaEtaMass", "#Delta#eta vs M(K#pi);#Delta#eta;M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, {200, -2.0, 2.0}, invmassAxis});
        hMC.add("KinematicsMisId/hOtherMisId_PtKaonPtMass", "Kaon p_{T} vs M(K#pi);p_{T}^{K} (GeV/c);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, ptAxis, invmassAxis});
        hMC.add("KinematicsMisId/hOtherMisId_PtPionPtMass", "Pion p_{T} vs M(K#pi);p_{T}^{#pi} (GeV/c);M(K#pi) (GeV/c^{2})", kTH3F, {ptAxis, ptAxis, invmassAxis});
      }
    }
  }

  template <typename Coll>
  bool selectionEvent(const Coll& collision, bool fillHist = false) // default to false
  {
    if (fillHist) {
      hEventSelection.fill(HIST("hEventCut"), 0);
    }

    if (std::abs(collision.posZ()) > selectionConfig.cfgVrtxZCut) {
      return false;
    }
    if (fillHist) {
      hEventSelection.fill(HIST("hEventCut"), 1);
    }

    if (selectionConfig.isApplysel8 && !collision.sel8()) {
      return false;
    }
    if (fillHist) {
      hEventSelection.fill(HIST("hEventCut"), 2);
    }

    if (selectionConfig.isNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (fillHist) {
      hEventSelection.fill(HIST("hEventCut"), 3);
    }

    if (selectionConfig.isNoITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (fillHist) {
      hEventSelection.fill(HIST("hEventCut"), 4);
    }

    if (selectionConfig.isTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (fillHist) {
      hEventSelection.fill(HIST("hEventCut"), 5);
    }

    if (selectionConfig.isNoSameBunchPileup && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup))) {
      return false;
    }
    if (fillHist) {
      hEventSelection.fill(HIST("hEventCut"), 6);
    }

    if (selectionConfig.isGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }
    if (fillHist) {
      hEventSelection.fill(HIST("hEventCut"), 7);
    }

    if (selectionConfig.isNoCollInTimeRangeStandard && (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))) {
      return false;
    }
    if (fillHist) {
      hEventSelection.fill(HIST("hEventCut"), 8);
    }

    if (selectionConfig.isApplyOccCut && (std::abs(collision.trackOccupancyInTimeRange()) > selectionConfig.cfgOccCut)) {
      return false;
    }
    if (fillHist) {
      hEventSelection.fill(HIST("hEventCut"), 9);
    }

    if (rctCut.requireRCTFlagChecker && !rctChecker(collision)) {
      return false;
    }
    if (fillHist) {
      hEventSelection.fill(HIST("hEventCut"), 10);
    }

    if (selectionConfig.isGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (fillHist) {
      hEventSelection.fill(HIST("hEventCut"), 11);
    }

    if (selectionConfig.isVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (fillHist) {
      hEventSelection.fill(HIST("hEventCut"), 12);
    }

    if (selectionConfig.isVertexTOFMatched && !collision.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if (fillHist) {
      hEventSelection.fill(HIST("hEventCut"), 13);
    }

    if (selectionConfig.isApplyINELgt0 && !collision.isInelGt0()) {
      return false;
    }
    if (fillHist) {
      hEventSelection.fill(HIST("hEventCut"), 14);
    }

    if (selectionConfig.isApplyhasFT0 && !collision.has_foundFT0()) {
      return false;
    }
    if (fillHist) {
      hEventSelection.fill(HIST("hEventCut"), 15);
    }

    return true;
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (selectionConfig.isGlobalTracks) {
      if (!candidate.isGlobalTrackWoDCA()) {
        return false;
      }
      if (std::abs(candidate.pt()) < selectionConfig.cfgCutPT) {
        return false;
      }

      if (std::abs(candidate.eta()) > selectionConfig.cfgCutEta) {
        return false;
      }
      if (!selectionConfig.isApplyPtDepDCAxyCut) {
        if (std::abs(candidate.dcaXY()) > selectionConfig.cfgCutDCAxy) {
          return false;
        }
      } else {
        if (std::abs(candidate.dcaXY()) > (0.0105 + 0.035 / std::pow(candidate.pt(), 1.1))) {
          return false;
        }
      }
      if (selectionConfig.isGoldenChi2 && !candidate.passedGoldenChi2()) {
        return false;
      }
      if (std::abs(candidate.dcaZ()) > selectionConfig.cfgCutDCAz) {
        return false;
      }
      if (candidate.itsNCls() < selectionConfig.cfgITScluster) {
        return false;
      }
      if (candidate.tpcNClsFound() < selectionConfig.cfgTPCcluster) {
        return false;
      }
      if (candidate.itsChi2NCl() >= selectionConfig.cfgITSChi2NCl) {
        return false;
      }
      if (candidate.tpcChi2NCl() >= selectionConfig.cfgTPCChi2NClMax || candidate.tpcChi2NCl() < selectionConfig.cfgTPCChi2NClMin) {
        return false;
      }
      if (selectionConfig.isPVContributor && !candidate.isPVContributor()) {
        return false;
      }
      // if (selectionConfig.isUseITSTPCRefit && (!(o2::aod::track::ITSrefit) || !(o2::aod::track::TPCrefit)))
      //   return false;
    } else {
      if (std::abs(candidate.pt()) < selectionConfig.cfgCutPT) {
        return false;
      }
      // if (std::abs(candidate.eta()) > selectionConfig.cfgCutEta || std::abs(candidate.eta()) < selectionConfig.cfgCutEtaMin)
      if (std::abs(candidate.eta()) > selectionConfig.cfgCutEta) {
        return false;
      }
      // if (std::abs(candidate.dcaXY()) > selectionConfig.cfgCutDCAxy || std::abs(candidate.dcaXY()) < selectionConfig.cfgCutDCAxyMin)
      if (std::abs(candidate.dcaXY()) > selectionConfig.cfgCutDCAxy) {
        return false;
      }
      if (std::abs(candidate.dcaZ()) > selectionConfig.cfgCutDCAz) {
        return false;
      }
      if (candidate.itsNCls() < selectionConfig.cfgITScluster) {
        return false;
      }
      if (candidate.tpcNClsFound() < selectionConfig.cfgTPCcluster) {
        return false;
      }
      if (candidate.itsChi2NCl() >= selectionConfig.cfgITSChi2NCl) {
        return false;
      }
      if (candidate.tpcChi2NCl() >= selectionConfig.cfgTPCChi2NClMax || candidate.tpcChi2NCl() < selectionConfig.cfgTPCChi2NClMin) {
        return false;
      }
      if (selectionConfig.isPVContributor && !candidate.isPVContributor()) {
        return false;
      }
      if (selectionConfig.isPrimaryTrack && !candidate.isPrimaryTrack()) {
        return false;
      }
    }

    return true;
  }

  // deep angle cut on pair to remove photon conversion
  template <typename T1, typename T2>
  bool selectionPair(const T1& candidate1, const T2& candidate2)
  {
    if (!selectionConfig.isApplyDeepAngle) {
      return true;
    }

    const double cosAngle = std::clamp(static_cast<double>((candidate1.pt() * candidate2.pt() + candidate1.pz() * candidate2.pz())) / (candidate1.p() * candidate2.p()), -1.0, 1.0);
    const double cosDeepAngle = std::cos(selectionConfig.cfgDeepAngle);
    return cosAngle <= cosDeepAngle;
  }

  template <typename T>
  bool isFakeTrack(const T& track, PIDParticle PID)
  {
    const auto pglobal = track.p();
    const auto ptpc = track.tpcInnerParam();
    if (PID == PIDParticle::kPion) {
      if (std::abs(pglobal - ptpc) > selectionConfig.cfgFakeTrackCutPi) {
        return true;
      }
    } else if (PID == PIDParticle::kKaon) {
      if (std::abs(pglobal - ptpc) > selectionConfig.cfgFakeTrackCutKa) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  float tpcSigma(const T& c, PIDParticle pid)
  {
    return pid == PIDParticle::kPion ? c.tpcNSigmaPi() : c.tpcNSigmaKa();
  }

  template <typename T>
  float tofSigma(const T& c, PIDParticle pid)
  {
    return pid == PIDParticle::kPion ? c.tofNSigmaPi() : c.tofNSigmaKa();
  }

  template <typename T>
  float combinedNSigma2(const T& c, PIDParticle pid)
  {
    const float tpc = tpcSigma(c, pid);
    const float tof = tofSigma(c, pid);

    return tpc * tpc + tof * tof;
  }

  float tpcCut(PIDParticle pid)
  {
    return pid == PIDParticle::kPion ? selectionConfig.nsigmaCutTPCPi : selectionConfig.nsigmaCutTPCKa;
  }

  float tofCut(PIDParticle pid)
  {
    return pid == PIDParticle::kPion ? selectionConfig.nsigmaCutTOFPi : selectionConfig.nsigmaCutTOFKa;
  }

  float combinedCut(PIDParticle pid)
  {
    return pid == PIDParticle::kPion ? selectionConfig.nsigmaCutCombinedPi : selectionConfig.nsigmaCutCombinedKa;
  }

  template <typename T>
  bool passTPC(const T& c, PIDParticle pid, float cut)
  {
    return std::abs(tpcSigma(c, pid)) < cut;
  }

  template <typename T>
  bool passTOF(const T& c, PIDParticle pid, float cut)
  {
    return c.hasTOF() && c.beta() > selectionConfig.cfgTOFBetaCut && std::abs(tofSigma(c, pid)) < cut;
  }

  template <typename T>
  bool passCombined(const T& c, PIDParticle pid, float cut)
  {
    if (!c.hasTOF() || c.beta() <= selectionConfig.cfgTOFBetaCut) {
      return false;
    }

    return combinedNSigma2(c, pid) < cut * cut;
  }

  template <typename T>
  bool selectionPID(const T& candidate, PIDParticle pid)
  {
    const auto strategy = static_cast<PIDStrategy>(selectionConfig.pidStrategy.value);
    const auto mode = static_cast<PIDMode>(selectionConfig.pidMode.value);

    switch (strategy) {
      case PIDStrategy::Standard: {
        switch (mode) {
          case PIDMode::TPC: // Apply only TPC cut
            return passTPC(candidate, pid, tpcCut(pid));

          case PIDMode::TOF: // Apply only TOF cut
            return passTOF(candidate, pid, tofCut(pid));

          case PIDMode::TOFVeto: // Require TOF hit; apply both TPC and TOF cuts separately. Tracks without TOF are rejected.
            return passTPC(candidate, pid, tpcCut(pid)) && passTOF(candidate, pid, tofCut(pid));

          case PIDMode::TOFHIT: // Apply only TOF cut that has tof hits, apply TPC cut for the rest
            if (candidate.hasTOF()) {
              return passTOF(candidate, pid, tofCut(pid));
            }

            return passTPC(candidate, pid, tpcCut(pid));

          case PIDMode::Combined: // Apply combined cut if TOF is available, apply TPC cut for the rest
            if (candidate.hasTOF()) {
              return passCombined(candidate, pid, combinedCut(pid));
            }

            return passTPC(candidate, pid, tpcCut(pid));

          default:
            LOGF(fatal, "Unknown PID mode!");
            return false;
        }
      }

      case PIDStrategy::PtDependent: // Apply TPC and TOF cut above lowPtCut, apply TPC cut for the rest
      {
        if (candidate.hasTOF() && candidate.pt() >= selectionConfig.lowPtCutPid) {
          return passTPC(candidate, pid, tpcCut(pid)) && passTOF(candidate, pid, tofCut(pid));
        }

        return passTPC(candidate, pid, tpcCut(pid));
      }

      case PIDStrategy::PtDependentTOF: // Apply combined cut (requires TOF) in specific pT window, apply TPC+TOF cut for tracks with TOF hit above lowPtCut, apply TPC cut for the rest.
      {
        const bool inPtWindow = (pid == PIDParticle::kPion) ? (candidate.pt() >= pionPidPtLow && candidate.pt() <= pionPidPtHigh) : (candidate.pt() >= kaonPidPtLow && candidate.pt() <= kaonPidPtHigh);

        if (inPtWindow) {
          return passCombined(candidate, pid, combinedCut(pid));
        }

        if (candidate.hasTOF() && candidate.pt() >= selectionConfig.lowPtCutPid) {
          return passTPC(candidate, pid, tpcCut(pid)) && passTOF(candidate, pid, tofCut(pid));
        }

        return passTPC(candidate, pid, tpcCut(pid));
      }

        /* case PIDStrategy::ThreePtDependent: // Apply pT-dependent TPC and TOF cuts using three pT regions
        {
          const int region = (candidate.pt() < selectionConfig.lowPtCutPid) ? 0 : (candidate.pt() < selectionConfig.highPtCutPid) ? 1
                                                                                                                                  : 2;

          const float regionTPCCut = (pid == PIDParticle::kPion) ? TPCPiCuts3Pt[region] : TPCKaCuts3Pt[region];

          const float regionTOFCut = (pid == PIDParticle::kPion) ? TOFPiCuts3Pt[region] : TOFKaCuts3Pt[region];

          if (candidate.hasTOF()) {
            return (std::abs(tpcSigma(candidate, pid)) < regionTPCCut) && (std::abs(tofSigma(candidate, pid)) < regionTOFCut) && (candidate.beta() > selectionConfig.cfgTOFBetaCut);
          }

          return (std::abs(tpcSigma(candidate, pid)) < regionTPCCut);
        } */

      case PIDStrategy::ThreePtDependent: // Apply region-dependent PID cuts using one PID cut per pT region. The middle pT region supports two TOF modes.
      {
        const float pidCut = (candidate.pt() < selectionConfig.lowPtCutPid) ? selectionConfig.lowPtPid : (candidate.pt() < selectionConfig.highPtCutPid) ? selectionConfig.midPtPid
                                                                                                                                                         : selectionConfig.highPtPid;

        // Low-pT region
        if (candidate.pt() < selectionConfig.lowPtCutPid) {

          if (!candidate.hasTOF()) {
            return passTPC(candidate, pid, pidCut);
          }

          return passCombined(candidate, pid, pidCut);
        }

        // Mid-pT region
        if (candidate.pt() < selectionConfig.highPtCutPid) {
          switch (mode) {
            case PIDMode::TOFVeto: // Mid-pT region (TOF mandatory with combined cut)
              if (!candidate.hasTOF() || candidate.beta() <= selectionConfig.cfgTOFBetaCut) {
                return false;
              }
              return passCombined(candidate, pid, pidCut);

            case PIDMode::Combined: // Mid-pT region (If TOF available use combined cut, otherwise use TPC cut)
              if (!candidate.hasTOF()) {
                return passTPC(candidate, pid, pidCut);
              }

              return passCombined(candidate, pid, pidCut);

            default:
              LOGF(fatal, "ThreePtDependent supports only Combined and TOFVeto modes.");
              return false;
          }
        }

        // High-pT region
        if (!candidate.hasTOF()) {
          return passTPC(candidate, pid, pidCut);
        }

        return passCombined(candidate, pid, pidCut);
      }

      case PIDStrategy::PIDCompare: // Apply independent TPC/TOF cuts below lowPtCutPid, within lowPtCutPid and highPtCutPid, require TOF and accept the track only if its combined nsigma is smaller than that of the alternate PID hypothesis and above highPtCutPid if TOF is available, accept if the track passes combined cut and if TOF is not available, check if it passes TPC cut
      {
        if (candidate.pt() < selectionConfig.lowPtCutPid) {
          if (candidate.hasTOF()) {
            return passTPC(candidate, pid, tpcCut(pid)) && passTOF(candidate, pid, tofCut(pid));
          }
          return passTPC(candidate, pid, tpcCut(pid));
        }

        if (candidate.pt() < selectionConfig.highPtCutPid) {
          if (!candidate.hasTOF() || candidate.beta() <= selectionConfig.cfgTOFBetaCut) {
            return false;
          }

          const float sigmaComb2 = combinedNSigma2(candidate, pid);

          const float cut = combinedCut(pid);
          if (sigmaComb2 >= cut * cut) {
            return false;
          }

          const PIDParticle otherPID = (pid == PIDParticle::kPion) ? PIDParticle::kKaon : PIDParticle::kPion;

          return sigmaComb2 < combinedNSigma2(candidate, otherPID);
        }

        if (candidate.hasTOF()) {
          return passCombined(candidate, pid, combinedCut(pid));
        }

        return passTPC(candidate, pid, tpcCut(pid));
      }

      default:
        LOGF(fatal, "Unknown PID strategy!");
    }
    return false;
  }

  template <typename T>
  bool isMisidentified(const T& candidate, PIDParticle pid)
  {
    const auto strategy = static_cast<MIDStrategy>(selectionConfig.misIdStrategy.value);
    const PIDParticle misPID = (pid == PIDParticle::kPion) ? PIDParticle::kKaon : PIDParticle::kPion;

    switch (strategy) {
      case MIDStrategy::Disabled:
        return false;

      case MIDStrategy::Standard: // Reject if track also satisfies TPC cut for the opposite hypothesis
        return std::abs(tpcSigma(candidate, misPID)) < selectionConfig.nsigmaCutTpcMisId;

      case MIDStrategy::PtDependent: // Reject if track is in MID pT window, has no TOF, and passes opposite TPC cut
      {
        const float lowPt = (pid == PIDParticle::kPion) ? selectionConfig.pionMisIdPtLow : selectionConfig.kaonMisIdPtLow;
        const float highPt = (pid == PIDParticle::kPion) ? selectionConfig.pionMisIdPtHigh : selectionConfig.kaonMisIdPtHigh;

        return candidate.pt() >= lowPt && candidate.pt() < highPt && !candidate.hasTOF() && std::abs(tpcSigma(candidate, misPID)) < selectionConfig.nsigmaCutTpcMisId;
      }

      case MIDStrategy::PtDependentCompare: // Reject if track is in MID pT window, has no TOF, and looks more like the opposite hypothesis than the signal hypothesis
      {
        const float lowPt = (pid == PIDParticle::kPion) ? selectionConfig.pionMisIdPtLow : selectionConfig.kaonMisIdPtLow;
        const float highPt = (pid == PIDParticle::kPion) ? selectionConfig.pionMisIdPtHigh : selectionConfig.kaonMisIdPtHigh;

        return candidate.pt() >= lowPt && candidate.pt() < highPt && !candidate.hasTOF() && std::abs(tpcSigma(candidate, misPID)) < std::abs(tpcSigma(candidate, pid));
      }

      default:
        LOGF(fatal, "Unknown MID strategy!");
    }
    return false;
  }

  template <typename MCParticle>
  std::vector<std::pair<int, int>> getAncestors(const MCParticle& particle)
  {
    std::vector<std::pair<int, int>> ancestors;

    std::function<void(const MCParticle&)> collect = [&](const auto& p) {
      for (const auto& mother : p.template mothers_as<aod::McParticles>()) {

        ancestors.emplace_back(mother.globalIndex(), std::abs(mother.pdgCode()));
        collect(mother);
      }
    };

    collect(particle);

    return ancestors;
  }

  //*********Varibles declaration***************
  ROOT::Math::PxPyPzMVector daughter1, daughter2, daughterRot, genDaughter1, genDaughter2, mother, motherRot, genMother;

  float centrality = -1.0f, theta2 = 0.0f;
  double genMass = 0.0, recMass = 0.0, recPt = 0.0, genPt = 0.0;

  bool isMix = false;

  template <typename T1, typename T2>
  void fillInvMass(const T1& argDaughter1, const T1& argDaughter2, const T1& argMother, float argCentrality, bool argIsMix, const T2& track1, const T2& track2)
  {
    if (track1.sign() * track2.sign() < 0) {
      if (!argIsMix) {
        if (argMother.Rapidity() > selectionConfig.motherRapidityMin && argMother.Rapidity() < selectionConfig.motherRapidityMax) {
          hInvMass.fill(HIST("h3KstarInvMassUnlikeSign"), argCentrality, argMother.Pt(), argMother.M());
        }
        for (int i = 0; i < cRotations; i++) {
          theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / rotationalCut, o2::constants::math::PI + o2::constants::math::PI / rotationalCut);

          daughterRot = ROOT::Math::PxPyPzMVector(argDaughter1.Px() * std::cos(theta2) - argDaughter1.Py() * std::sin(theta2), argDaughter1.Px() * std::sin(theta2) + argDaughter1.Py() * std::cos(theta2), argDaughter1.Pz(), argDaughter1.M());
          motherRot = daughterRot + argDaughter2;

          if (calcRotational && (motherRot.Rapidity() > selectionConfig.motherRapidityMin && motherRot.Rapidity() < selectionConfig.motherRapidityMax)) {
            hInvMass.fill(HIST("h3KstarInvMassRotated"), argCentrality, motherRot.Pt(), motherRot.M());
          }
        }
      } else if (argMother.Rapidity() > selectionConfig.motherRapidityMin && argMother.Rapidity() < selectionConfig.motherRapidityMax) {
        hInvMass.fill(HIST("h3KstarInvMassMixed"), argCentrality, argMother.Pt(), argMother.M());
      }
    } else {
      if (!argIsMix) {
        if (calcLikeSign && (argMother.Rapidity() > selectionConfig.motherRapidityMin && argMother.Rapidity() < selectionConfig.motherRapidityMax)) {
          if (track1.sign() > 0 && track2.sign() > 0) {
            hInvMass.fill(HIST("h3KstarInvMasslikeSignPP"), argCentrality, argMother.Pt(), argMother.M());
          } else if (track1.sign() < 0 && track2.sign() < 0) {
            hInvMass.fill(HIST("h3KstarInvMasslikeSignMM"), argCentrality, argMother.Pt(), argMother.M());
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

    if (selectCentEstimator == kFT0A) {
      centrality = collision.centFT0A();
    } else if (selectCentEstimator == kFT0C) {
      centrality = collision.centFT0C();
    } else if (selectCentEstimator == kFV0A) {
      centrality = collision.centFV0A();
    } else {
      centrality = collision.centFT0M(); // default includes kFT0M
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

      if (track1.globalIndex() == track2.globalIndex()) {
        continue;
      }

      if (!selectionPair(track1, track2)) {
        continue;
      }

      if (cQAplots) {
        hOthers.fill(HIST("dE_by_dx_TPC"), track1.p(), track1.tpcSignal());

        if (track1.sign() < 0) {
          hPID.fill(HIST("Before/hTPCnsigKa_Neg_mult_pt"), track1.tpcNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTPCnsigPi_Neg_mult_pt"), track1.tpcNSigmaPi(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTOFnsigKa_Neg_mult_pt"), track1.tofNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTOFnsigPi_Neg_mult_pt"), track1.tofNSigmaPi(), centrality, track1.pt());

          hPID.fill(HIST("Before/hTPCnsigKa_Neg_mult_p"), track1.tpcNSigmaKa(), centrality, track1.p());
          hPID.fill(HIST("Before/hTPCnsigPi_Neg_mult_p"), track1.tpcNSigmaPi(), centrality, track1.p());
          hPID.fill(HIST("Before/hTOFnsigKa_Neg_mult_p"), track1.tofNSigmaKa(), centrality, track1.p());
          hPID.fill(HIST("Before/hTOFnsigPi_Neg_mult_p"), track1.tofNSigmaPi(), centrality, track1.p());
        } else if (track1.sign() > 0) {
          hPID.fill(HIST("Before/hTPCnsigKa_Pos_mult_pt"), track1.tpcNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTPCnsigPi_Pos_mult_pt"), track1.tpcNSigmaPi(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTOFnsigKa_Pos_mult_pt"), track1.tofNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTOFnsigPi_Pos_mult_pt"), track1.tofNSigmaPi(), centrality, track1.pt());

          hPID.fill(HIST("Before/hTPCnsigKa_Pos_mult_p"), track1.tpcNSigmaKa(), centrality, track1.p());
          hPID.fill(HIST("Before/hTPCnsigPi_Pos_mult_p"), track1.tpcNSigmaPi(), centrality, track1.p());
          hPID.fill(HIST("Before/hTOFnsigKa_Pos_mult_p"), track1.tofNSigmaKa(), centrality, track1.p());
          hPID.fill(HIST("Before/hTOFnsigPi_Pos_mult_p"), track1.tofNSigmaPi(), centrality, track1.p());
        }

        hPID.fill(HIST("Before/hNsigma_TPC_TOF_Ka_pt"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
        hPID.fill(HIST("Before/hNsigma_TPC_TOF_Pi_pt"), track1.tpcNSigmaPi(), track1.tofNSigmaPi(), track1.pt());

        hOthers.fill(HIST("hDcaxy_cent_pt"), track1.dcaXY(), centrality, track1.pt());
        hOthers.fill(HIST("hDcaz_cent_pt"), track1.dcaZ(), centrality, track1.pt());
      }

      // since we are using combinations full index policy, so repeated pairs are allowed, so we can check one with Kaon and other with pion
      if (!selectionPID(track1, PIDParticle::kKaon)) { // Track 1 is checked with Kaon
        continue;
      }

      if (!selectionPID(track2, PIDParticle::kPion)) { // Track 2 is checked with Pion
        continue;
      }

      if (isMisidentified(track1, PIDParticle::kKaon)) {
        continue;
      }
      if (isMisidentified(track2, PIDParticle::kPion)) {
        continue;
      }

      if (selectionConfig.isApplyFakeTrack && (isFakeTrack(track1, PIDParticle::kKaon) || isFakeTrack(track2, PIDParticle::kPion))) {
        continue;
      }

      if (cQAplots) {
        hOthers.fill(HIST("hEta_after"), track1.eta());
        hOthers.fill(HIST("hDcaxyPi"), track2.dcaXY());
        hOthers.fill(HIST("hDcaxyKa"), track1.dcaXY());
        hOthers.fill(HIST("hDcazPi"), track2.dcaZ());
        hOthers.fill(HIST("hDcazKa"), track1.dcaZ());

        if (track1.sign() < 0) {
          hPID.fill(HIST("After/hTPCnsigKa_Neg_mult_pt"), track1.tpcNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("After/hTOFnsigKa_Neg_mult_pt"), track1.tofNSigmaKa(), centrality, track1.pt());

          hPID.fill(HIST("After/hTPCnsigKa_Neg_mult_p"), track1.tpcNSigmaKa(), centrality, track1.p());
          hPID.fill(HIST("After/hTOFnsigKa_Neg_mult_p"), track1.tofNSigmaKa(), centrality, track1.p());
        } else if (track1.sign() > 0) {
          hPID.fill(HIST("After/hTPCnsigKa_Pos_mult_pt"), track1.tpcNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("After/hTOFnsigKa_Pos_mult_pt"), track1.tofNSigmaKa(), centrality, track1.pt());

          hPID.fill(HIST("After/hTPCnsigKa_Pos_mult_p"), track1.tpcNSigmaKa(), centrality, track1.p());
          hPID.fill(HIST("After/hTOFnsigKa_Pos_mult_p"), track1.tofNSigmaKa(), centrality, track1.p());
        }

        if (track2.sign() < 0) {
          hPID.fill(HIST("After/hTPCnsigPi_Neg_mult_pt"), track2.tpcNSigmaPi(), centrality, track2.pt());
          hPID.fill(HIST("After/hTOFnsigPi_Neg_mult_pt"), track2.tofNSigmaPi(), centrality, track2.pt());

          hPID.fill(HIST("After/hTPCnsigPi_Neg_mult_p"), track2.tpcNSigmaPi(), centrality, track2.p());
          hPID.fill(HIST("After/hTOFnsigPi_Neg_mult_p"), track2.tofNSigmaPi(), centrality, track2.p());
        } else if (track2.sign() > 0) {
          hPID.fill(HIST("After/hTPCnsigPi_Pos_mult_pt"), track2.tpcNSigmaPi(), centrality, track2.pt());
          hPID.fill(HIST("After/hTOFnsigPi_Pos_mult_pt"), track2.tofNSigmaPi(), centrality, track2.pt());

          hPID.fill(HIST("After/hTPCnsigPi_Pos_mult_p"), track2.tpcNSigmaPi(), centrality, track2.p());
          hPID.fill(HIST("After/hTOFnsigPi_Pos_mult_p"), track2.tofNSigmaPi(), centrality, track2.p());
        }

        hPID.fill(HIST("After/hNsigma_TPC_TOF_Ka_pt"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
        hPID.fill(HIST("After/hNsigma_TPC_TOF_Pi_pt"), track2.tpcNSigmaPi(), track2.tofNSigmaPi(), track2.pt());

        if (track1.hasTOF()) {
          hOthers.fill(HIST("hTOFBetaKa"), track1.beta());
        }
        if (track2.hasTOF()) {
          hOthers.fill(HIST("hTOFBetaPi"), track2.beta());
        }
        hOthers.fill(HIST("hTPCpDiffKa"), track1.p() - track1.tpcInnerParam(), track1.pt());
        hOthers.fill(HIST("hTPCpDiffPi"), track2.p() - track2.tpcInnerParam(), track2.pt());
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
      if (cQAplots) {
        hOthers.fill(HIST("hKstar_rap_pt"), mother.Rapidity(), mother.Pt());
        hOthers.fill(HIST("hKstar_eta_pt"), mother.Eta(), mother.Pt());
      }

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

    if (selectCentEstimator == kFT0A) {
      centrality = collision.centFT0A();
    } else if (selectCentEstimator == kFT0C) {
      centrality = collision.centFT0C();
    } else if (selectCentEstimator == kFV0A) {
      centrality = collision.centFV0A();
    } else {
      centrality = collision.centFT0M(); // default includes kFT0M
    }

    // Fill the event counter
    if (cQAevents) {
      hEventSelection.fill(HIST("hVertexZ"), collision.posZ());
      hEventSelection.fill(HIST("hCentrality"), centrality);
    }

    for (const auto& [track1, track2] : combinations(CombinationsFullIndexPolicy(tracks, tracks))) {
      if (!selectionTrack(track1) || !selectionTrack(track2)) {
        continue;
      }

      const auto mctrack1 = track1.mcParticle();
      const auto mctrack2 = track2.mcParticle();

      if (!track1.has_mcParticle() || !track2.has_mcParticle()) {
        continue; // skip if no MC particle associated
      }

      if (!mctrack1.isPhysicalPrimary() || !mctrack2.isPhysicalPrimary()) {
        continue;
      }

      if (track1.globalIndex() == track2.globalIndex()) {
        continue;
      }

      if (cQAplots) {
        hOthers.fill(HIST("dE_by_dx_TPC"), track1.p(), track1.tpcSignal());

        if (track1.sign() < 0) {
          hPID.fill(HIST("Before/hTPCnsigKa_Neg_mult_pt"), track1.tpcNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTPCnsigPi_Neg_mult_pt"), track1.tpcNSigmaPi(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTOFnsigKa_Neg_mult_pt"), track1.tofNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTOFnsigPi_Neg_mult_pt"), track1.tofNSigmaPi(), centrality, track1.pt());

          hPID.fill(HIST("Before/hTPCnsigKa_Neg_mult_p"), track1.tpcNSigmaKa(), centrality, track1.p());
          hPID.fill(HIST("Before/hTPCnsigPi_Neg_mult_p"), track1.tpcNSigmaPi(), centrality, track1.p());
          hPID.fill(HIST("Before/hTOFnsigKa_Neg_mult_p"), track1.tofNSigmaKa(), centrality, track1.p());
          hPID.fill(HIST("Before/hTOFnsigPi_Neg_mult_p"), track1.tofNSigmaPi(), centrality, track1.p());
        } else {
          hPID.fill(HIST("Before/hTPCnsigKa_Pos_mult_pt"), track1.tpcNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTPCnsigPi_Pos_mult_pt"), track1.tpcNSigmaPi(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTOFnsigKa_Pos_mult_pt"), track1.tofNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("Before/hTOFnsigPi_Pos_mult_pt"), track1.tofNSigmaPi(), centrality, track1.pt());

          hPID.fill(HIST("Before/hTPCnsigKa_Pos_mult_p"), track1.tpcNSigmaKa(), centrality, track1.p());
          hPID.fill(HIST("Before/hTPCnsigPi_Pos_mult_p"), track1.tpcNSigmaPi(), centrality, track1.p());
          hPID.fill(HIST("Before/hTOFnsigKa_Pos_mult_p"), track1.tofNSigmaKa(), centrality, track1.p());
          hPID.fill(HIST("Before/hTOFnsigPi_Pos_mult_p"), track1.tofNSigmaPi(), centrality, track1.p());
        }

        hPID.fill(HIST("Before/hNsigma_TPC_TOF_Ka_pt"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
        hPID.fill(HIST("Before/hNsigma_TPC_TOF_Pi_pt"), track1.tpcNSigmaPi(), track1.tofNSigmaPi(), track1.pt());

        hOthers.fill(HIST("hDcaxy_cent_pt"), track1.dcaXY(), centrality, track1.pt());
        hOthers.fill(HIST("hDcaz_cent_pt"), track1.dcaZ(), centrality, track1.pt());
      }

      // since we are using combinations full index policy, so repeated pairs are allowed, so we can check one with Kaon and other with pion
      if (!selectionPID(track1, PIDParticle::kKaon)) { // Track 1 is checked with Kaon
        continue;
      }

      if (!selectionPID(track2, PIDParticle::kPion)) { // Track 2 is checked with Pion
        continue;
      }

      if (isMisidentified(track1, PIDParticle::kKaon)) {
        continue;
      }
      if (isMisidentified(track2, PIDParticle::kPion)) {
        continue;
      }

      if (selectionConfig.isApplyFakeTrack && (isFakeTrack(track1, PIDParticle::kKaon) || isFakeTrack(track2, PIDParticle::kPion))) {
        continue;
      }

      if (cQAplots) {
        hOthers.fill(HIST("hEta_after"), track1.eta());
        hOthers.fill(HIST("hDcaxyPi"), track2.dcaXY());
        hOthers.fill(HIST("hDcaxyKa"), track1.dcaXY());
        hOthers.fill(HIST("hDcazPi"), track2.dcaZ());
        hOthers.fill(HIST("hDcazKa"), track1.dcaZ());

        if (track1.sign() < 0) {
          hPID.fill(HIST("After/hTPCnsigKa_Neg_mult_pt"), track1.tpcNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("After/hTOFnsigKa_Neg_mult_pt"), track1.tofNSigmaKa(), centrality, track1.pt());

          hPID.fill(HIST("After/hTPCnsigKa_Neg_mult_p"), track1.tpcNSigmaKa(), centrality, track1.p());
          hPID.fill(HIST("After/hTOFnsigKa_Neg_mult_p"), track1.tofNSigmaKa(), centrality, track1.p());
        } else if (track1.sign() > 0) {
          hPID.fill(HIST("After/hTPCnsigKa_Pos_mult_pt"), track1.tpcNSigmaKa(), centrality, track1.pt());
          hPID.fill(HIST("After/hTOFnsigKa_Pos_mult_pt"), track1.tofNSigmaKa(), centrality, track1.pt());

          hPID.fill(HIST("After/hTPCnsigKa_Pos_mult_p"), track1.tpcNSigmaKa(), centrality, track1.p());
          hPID.fill(HIST("After/hTOFnsigKa_Pos_mult_p"), track1.tofNSigmaKa(), centrality, track1.p());
        }

        if (track2.sign() < 0) {
          hPID.fill(HIST("After/hTPCnsigPi_Neg_mult_pt"), track2.tpcNSigmaPi(), centrality, track2.pt());
          hPID.fill(HIST("After/hTOFnsigPi_Neg_mult_pt"), track2.tofNSigmaPi(), centrality, track2.pt());

          hPID.fill(HIST("After/hTPCnsigPi_Neg_mult_p"), track2.tpcNSigmaPi(), centrality, track2.p());
          hPID.fill(HIST("After/hTOFnsigPi_Neg_mult_p"), track2.tofNSigmaPi(), centrality, track2.p());
        } else if (track2.sign() > 0) {
          hPID.fill(HIST("After/hTPCnsigPi_Pos_mult_pt"), track2.tpcNSigmaPi(), centrality, track2.pt());
          hPID.fill(HIST("After/hTOFnsigPi_Pos_mult_pt"), track2.tofNSigmaPi(), centrality, track2.pt());

          hPID.fill(HIST("After/hTPCnsigPi_Pos_mult_p"), track2.tpcNSigmaPi(), centrality, track2.p());
          hPID.fill(HIST("After/hTOFnsigPi_Pos_mult_p"), track2.tofNSigmaPi(), centrality, track2.p());
        }

        hPID.fill(HIST("After/hNsigma_TPC_TOF_Ka_pt"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
        hPID.fill(HIST("After/hNsigma_TPC_TOF_Pi_pt"), track2.tpcNSigmaPi(), track2.tofNSigmaPi(), track2.pt());

        if (track1.hasTOF()) {
          hOthers.fill(HIST("hTOFBetaKa"), track1.beta());
        }
        if (track2.hasTOF()) {
          hOthers.fill(HIST("hTOFBetaPi"), track2.beta());
        }
        hOthers.fill(HIST("hTPCpDiffKa"), track1.p() - track1.tpcInnerParam(), track1.pt());
        hOthers.fill(HIST("hTPCpDiffPi"), track2.p() - track2.tpcInnerParam(), track2.pt());
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
      if (cQAplots) {
        hOthers.fill(HIST("hKstar_rap_pt"), mother.Rapidity(), mother.Pt());
        hOthers.fill(HIST("hKstar_eta_pt"), mother.Eta(), mother.Pt());
      }
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
    auto runMixing = [&](const auto& pair, auto centralityGetter) {
      for (const auto& [c1, tracks1, c2, tracks2] : pair) {

        if (!selectionEvent(c1, false) || !selectionEvent(c2, false)) { // don't fill event cut histogram
          continue;
        }

        centrality = centralityGetter(c1);

        for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
          if (!selectionTrack(t1) || !selectionTrack(t2)) {
            continue;
          }

          if (!selectionPair(t1, t2)) {
            continue;
          }

          // t1 checked as kaon and t2 checked as pion
          if (!selectionPID(t1, PIDParticle::kKaon) || !selectionPID(t2, PIDParticle::kPion)) {
            continue;
          }

          if (isMisidentified(t1, PIDParticle::kKaon) || isMisidentified(t2, PIDParticle::kPion)) {
            continue;
          }

          if (selectionConfig.isApplyFakeTrack && (isFakeTrack(t1, PIDParticle::kKaon) || isFakeTrack(t2, PIDParticle::kPion))) {
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
    auto runMixing = [&](const auto& pair, auto centralityGetter) {
      for (const auto& [c1, tracks1, c2, tracks2] : pair) {

        if (!selectionEvent(c1, false) || !selectionEvent(c2, false)) { // don't fill event cut histogram
          continue;
        }

        if (!c1.has_mcCollision() || !c2.has_mcCollision()) {
          continue; // skip if no MC collision associated
        }

        centrality = centralityGetter(c1);

        for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
          if (!selectionTrack(t1) || !selectionTrack(t2)) {
            continue;
          }

          if (!selectionPair(t1, t2)) {
            continue;
          }

          // t1 checked as kaon and t2 checked as pion
          if (!selectionPID(t1, PIDParticle::kKaon) || !selectionPID(t2, PIDParticle::kPion)) {
            continue;
          }

          if (isMisidentified(t1, PIDParticle::kKaon) || isMisidentified(t2, PIDParticle::kPion)) {
            continue;
          }

          if (selectionConfig.isApplyFakeTrack && (isFakeTrack(t1, PIDParticle::kKaon) || isFakeTrack(t2, PIDParticle::kPion))) {
            continue;
          }

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

    if (selectionConfig.isApplyMCGenTVX && (mcCollision.multMCFT0C() <= 0 || mcCollision.multMCFT0A() <= 0)) {
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

      if (selectCentEstimator == kFT0A) {
        centrality = collision.centFT0A();
      } else if (selectCentEstimator == kFT0C) {
        centrality = collision.centFT0C();
      } else if (selectCentEstimator == kFV0A) {
        centrality = collision.centFV0A();
      } else {
        centrality = collision.centFT0M(); // default includes kFT0M
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
        int charge = static_cast<int>(pdgDau > 0) - static_cast<int>(pdgDau < 0);

        if (charge > 0) {
          hasPos = true;
        }
        if (charge < 0) {
          hasNeg = true;
        }

        if (std::abs(pdgDau) == PDG_t::kKPlus) {
          passkaon = true;
          daughter1 = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
        } else if (std::abs(pdgDau) == PDG_t::kPiPlus) {
          daughter2 = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massPi);
          passpion = true;
        }
      }

      if ((passkaon && passpion) && (hasPos && hasNeg)) {
        mother = daughter1 + daughter2;
        hMC.fill(HIST("Gen/h3KstarPtCentMassDirect"), mcParticle.pt(), centrality, std::sqrt(mcParticle.e() * mcParticle.e() - mcParticle.p() * mcParticle.p()));
        hMC.fill(HIST("Gen/h3KstarPtCentMassRec"), mother.Pt(), centrality, mother.M());
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

    if (selectCentEstimator == kFT0A) {
      centrality = collision.centFT0A();
    } else if (selectCentEstimator == kFT0C) {
      centrality = collision.centFT0C();
    } else if (selectCentEstimator == kFV0A) {
      centrality = collision.centFV0A();
    } else {
      centrality = collision.centFT0M(); // default includes kFT0M
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

      if (track1.index() <= track2.index()) {
        continue;
      }

      if (!selectionPair(track1, track2)) {
        continue;
      }

      if (cQAplots) {
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

          hPID.fill(HIST("Before/hTPCnsigPi_Pos_mult_p"), track2.tpcNSigmaPi(), centrality, track2.p());
          hPID.fill(HIST("Before/hTOFnsigPi_Pos_mult_p"), track2.tofNSigmaPi(), centrality, track2.p());
        }
        if (mctrack2.pdgCode() == PDG_t::kKPlus) { // kaon
          hPID.fill(HIST("Before/hTPCnsigKa_Pos_mult_pt"), track2.tpcNSigmaKa(), centrality, track2.pt());
          hPID.fill(HIST("Before/hTOFnsigKa_Pos_mult_pt"), track2.tofNSigmaKa(), centrality, track2.pt());

          hPID.fill(HIST("Before/hTPCnsigKa_Pos_mult_p"), track2.tpcNSigmaKa(), centrality, track2.p());
          hPID.fill(HIST("Before/hTOFnsigKa_Pos_mult_p"), track2.tofNSigmaKa(), centrality, track2.p());
        }
        if (mctrack2.pdgCode() == PDG_t::kPiMinus) { // negative track pion
          hPID.fill(HIST("Before/hTPCnsigPi_Neg_mult_pt"), track2.tpcNSigmaPi(), centrality, track2.pt());
          hPID.fill(HIST("Before/hTOFnsigPi_Neg_mult_pt"), track2.tofNSigmaPi(), centrality, track2.pt());

          hPID.fill(HIST("Before/hTPCnsigPi_Neg_mult_p"), track2.tpcNSigmaPi(), centrality, track2.p());
          hPID.fill(HIST("Before/hTOFnsigPi_Neg_mult_p"), track2.tofNSigmaPi(), centrality, track2.p());
        }
        if (mctrack2.pdgCode() == PDG_t::kKMinus) { // negative track kaon
          hPID.fill(HIST("Before/hTPCnsigKa_Neg_mult_pt"), track2.tpcNSigmaKa(), centrality, track2.pt());
          hPID.fill(HIST("Before/hTOFnsigKa_Neg_mult_pt"), track2.tofNSigmaKa(), centrality, track2.pt());

          hPID.fill(HIST("Before/hTPCnsigKa_Neg_mult_p"), track2.tpcNSigmaKa(), centrality, track2.p());
          hPID.fill(HIST("Before/hTOFnsigKa_Neg_mult_p"), track2.tofNSigmaKa(), centrality, track2.p());
        }
        if (std::abs(mctrack1.pdgCode()) == PDG_t::kKPlus && std::abs(mctrack2.pdgCode()) == PDG_t::kPiPlus) {
          hPID.fill(HIST("Before/hNsigma_TPC_TOF_Ka_pt"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
          hPID.fill(HIST("Before/hNsigma_TPC_TOF_Pi_pt"), track2.tpcNSigmaPi(), track2.tofNSigmaPi(), track2.pt());
        }
      }

      if ((track1PDG != PDG_t::kPiPlus || track2PDG != PDG_t::kKPlus) && (track1PDG != PDG_t::kKPlus || track2PDG != PDG_t::kPiPlus)) {
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
            if (!selectionPID(track1, PIDParticle::kPion) || !selectionPID(track2, PIDParticle::kKaon)) { // Treat track1 as the pion candidate and track2 as the kaon candidate
              continue;
            }

            if (isMisidentified(track1, PIDParticle::kPion) || isMisidentified(track2, PIDParticle::kKaon)) {
              continue;
            }

            if (selectionConfig.isApplyFakeTrack && (isFakeTrack(track1, PIDParticle::kPion) || isFakeTrack(track2, PIDParticle::kKaon))) {
              continue;
            }

            if (cQAplots) {
              if (track1.sign() < 0 && track2.sign() > 0) {
                hPID.fill(HIST("After/hTPCnsigPi_Neg_mult_pt"), track1.tpcNSigmaPi(), centrality, track1.pt());
                hPID.fill(HIST("After/hTOFnsigPi_Neg_mult_pt"), track1.tofNSigmaPi(), centrality, track1.pt());
                hPID.fill(HIST("After/hTPCnsigKa_Pos_mult_pt"), track2.tpcNSigmaKa(), centrality, track2.pt());
                hPID.fill(HIST("After/hTOFnsigKa_Pos_mult_pt"), track2.tofNSigmaKa(), centrality, track2.pt());

                hPID.fill(HIST("After/hTPCnsigPi_Neg_mult_p"), track1.tpcNSigmaPi(), centrality, track1.p());
                hPID.fill(HIST("After/hTOFnsigPi_Neg_mult_p"), track1.tofNSigmaPi(), centrality, track1.p());
                hPID.fill(HIST("After/hTPCnsigKa_Pos_mult_p"), track2.tpcNSigmaKa(), centrality, track2.p());
                hPID.fill(HIST("After/hTOFnsigKa_Pos_mult_p"), track2.tofNSigmaKa(), centrality, track2.p());
              } else {
                hPID.fill(HIST("After/hTPCnsigPi_Pos_mult_pt"), track1.tpcNSigmaPi(), centrality, track1.pt());
                hPID.fill(HIST("After/hTOFnsigPi_Pos_mult_pt"), track1.tofNSigmaPi(), centrality, track1.pt());
                hPID.fill(HIST("After/hTPCnsigKa_Neg_mult_pt"), track2.tpcNSigmaKa(), centrality, track2.pt());
                hPID.fill(HIST("After/hTOFnsigKa_Neg_mult_pt"), track2.tofNSigmaKa(), centrality, track2.pt());

                hPID.fill(HIST("After/hTPCnsigPi_Pos_mult_p"), track1.tpcNSigmaPi(), centrality, track1.p());
                hPID.fill(HIST("After/hTOFnsigPi_Pos_mult_p"), track1.tofNSigmaPi(), centrality, track1.p());
                hPID.fill(HIST("After/hTPCnsigKa_Neg_mult_p"), track2.tpcNSigmaKa(), centrality, track2.p());
                hPID.fill(HIST("After/hTOFnsigKa_Neg_mult_p"), track2.tofNSigmaKa(), centrality, track2.p());
              }

              if (track1.hasTOF()) {
                hOthers.fill(HIST("hTOFBetaPi"), track1.beta());
              }
              if (track2.hasTOF()) {
                hOthers.fill(HIST("hTOFBetaKa"), track2.beta());
              }
              hOthers.fill(HIST("hTPCpDiffPi"), track1.p() - track1.tpcInnerParam(), track1.pt());
              hOthers.fill(HIST("hTPCpDiffKa"), track2.p() - track2.tpcInnerParam(), track2.pt());
            }

          } else if (track1PDG == PDG_t::kKPlus) {
            if (!selectionPID(track1, PIDParticle::kKaon) || !selectionPID(track2, PIDParticle::kPion)) { // Treat track1 as the kaon candidate and track2 as the pion candidate
              continue;
            }

            if (isMisidentified(track1, PIDParticle::kKaon) || isMisidentified(track2, PIDParticle::kPion)) {
              continue;
            }

            if (selectionConfig.isApplyFakeTrack && (isFakeTrack(track1, PIDParticle::kKaon) || isFakeTrack(track2, PIDParticle::kPion))) {
              continue;
            }

            if (cQAplots) {
              if (track1.sign() < 0 && track2.sign() > 0) {
                hPID.fill(HIST("After/hTPCnsigKa_Neg_mult_pt"), track1.tpcNSigmaKa(), centrality, track1.pt());
                hPID.fill(HIST("After/hTOFnsigKa_Neg_mult_pt"), track1.tofNSigmaKa(), centrality, track1.pt());
                hPID.fill(HIST("After/hTPCnsigPi_Pos_mult_pt"), track2.tpcNSigmaPi(), centrality, track2.pt());
                hPID.fill(HIST("After/hTOFnsigPi_Pos_mult_pt"), track2.tofNSigmaPi(), centrality, track2.pt());

                hPID.fill(HIST("After/hTPCnsigKa_Neg_mult_p"), track1.tpcNSigmaKa(), centrality, track1.p());
                hPID.fill(HIST("After/hTOFnsigKa_Neg_mult_p"), track1.tofNSigmaKa(), centrality, track1.p());
                hPID.fill(HIST("After/hTPCnsigPi_Pos_mult_p"), track2.tpcNSigmaPi(), centrality, track2.p());
                hPID.fill(HIST("After/hTOFnsigPi_Pos_mult_p"), track2.tofNSigmaPi(), centrality, track2.p());
              } else {
                hPID.fill(HIST("After/hTPCnsigKa_Pos_mult_pt"), track1.tpcNSigmaKa(), centrality, track1.pt());
                hPID.fill(HIST("After/hTOFnsigKa_Pos_mult_pt"), track1.tofNSigmaKa(), centrality, track1.pt());
                hPID.fill(HIST("After/hTPCnsigPi_Neg_mult_pt"), track2.tpcNSigmaPi(), centrality, track2.pt());
                hPID.fill(HIST("After/hTOFnsigPi_Neg_mult_pt"), track2.tofNSigmaPi(), centrality, track2.pt());

                hPID.fill(HIST("After/hTPCnsigKa_Pos_mult_p"), track1.tpcNSigmaKa(), centrality, track1.p());
                hPID.fill(HIST("After/hTOFnsigKa_Pos_mult_p"), track1.tofNSigmaKa(), centrality, track1.p());
                hPID.fill(HIST("After/hTPCnsigPi_Neg_mult_p"), track2.tpcNSigmaPi(), centrality, track2.p());
                hPID.fill(HIST("After/hTOFnsigPi_Neg_mult_p"), track2.tofNSigmaPi(), centrality, track2.p());
              }

              if (track1.hasTOF()) {
                hOthers.fill(HIST("hTOFBetaKa"), track1.beta());
              }
              if (track2.hasTOF()) {
                hOthers.fill(HIST("hTOFBetaPi"), track2.beta());
              }
              hOthers.fill(HIST("hTPCpDiffKa"), track1.p() - track1.tpcInnerParam(), track1.pt());
              hOthers.fill(HIST("hTPCpDiffPi"), track2.p() - track2.tpcInnerParam(), track2.pt());
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

            genDaughter1 = ROOT::Math::PxPyPzMVector(mctrack1.px(), mctrack1.py(), mctrack1.pz(), massPi);
            genDaughter2 = ROOT::Math::PxPyPzMVector(mctrack2.px(), mctrack2.py(), mctrack2.pz(), massKa);

          } else if (track1PDG == PDG_t::kKPlus) {
            daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
            daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);

            genDaughter1 = ROOT::Math::PxPyPzMVector(mctrack1.px(), mctrack1.py(), mctrack1.pz(), massKa);
            genDaughter2 = ROOT::Math::PxPyPzMVector(mctrack2.px(), mctrack2.py(), mctrack2.pz(), massPi);
          }

          mother = daughter1 + daughter2;          // Rec Kstar meson
          genMother = genDaughter1 + genDaughter2; // Gen Kstar from MC daughters

          if (mother.Rapidity() < selectionConfig.motherRapidityMin || mother.Rapidity() > selectionConfig.motherRapidityMax) {
            continue;
          }

          recMass = mother.M();
          recPt = mother.Pt();

          genMass = genMother.M();
          genPt = mothertrack1.pt();

          hMC.fill(HIST("Rec/hMassShift"), genPt, recPt, recMass - genMass);
          hMC.fill(HIST("Rec/h3KstarPtCentMassgen"), genPt, centrality, genMass);
          hMC.fill(HIST("Rec/h3KstarPtCentMassrec"), recPt, centrality, recMass);
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

    if (selectionConfig.isApplyMCGenTVX && (mcCollision.multMCFT0C() <= 0 || mcCollision.multMCFT0A() <= 0)) {
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
      if (!RecCollision.has_mcCollision()) {
        continue;
      }
      if (!selectionEvent(RecCollision, false)) { // don't fill event cut histogram
        continue;
      }

      if (selectCentEstimator == kFT0A) {
        centrality = RecCollision.centFT0A();
      } else if (selectCentEstimator == kFT0C) {
        centrality = RecCollision.centFT0C();
      } else if (selectCentEstimator == kFV0A) {
        centrality = RecCollision.centFV0A();
      } else {
        centrality = RecCollision.centFT0M(); // default includes kFT0M
      }

      isSelectedEvent = true;
    }

    if (isSelectedEvent) {
      hMC.fill(HIST("ImpactCorr/hImpactParameterRec"), impactPar);
      hMC.fill(HIST("ImpactCorr/hImpactParvsCentrRec"), centrality, impactPar);
    }

    // Generated MC
    for (const auto& mcPart : mcParticles) {
      if ((mcPart.y() < selectionConfig.motherRapidityMin || mcPart.y() > selectionConfig.motherRapidityMax) || std::abs(mcPart.pdgCode()) != o2::constants::physics::kK0Star892) {
        continue;
      }

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

    if (selectionConfig.isApplyMCGenTVX && (mcCollision.multMCFT0C() <= 0 || mcCollision.multMCFT0A() <= 0)) {
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

      if (!selectionEvent(collision, false)) {
        continue;
      }

      if (selectCentEstimator == kFT0A) {
        centrality = collision.centFT0A();
      } else if (selectCentEstimator == kFT0C) {
        centrality = collision.centFT0C();
      } else if (selectCentEstimator == kFV0A) {
        centrality = collision.centFV0A();
      } else {
        centrality = collision.centFT0M(); // default includes kFT0M
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

      if ((mcPart.y() < selectionConfig.motherRapidityMin || mcPart.y() > selectionConfig.motherRapidityMax) || std::abs(mcPart.pdgCode()) != o2::constants::physics::kK0Star892) {
        continue;
      }

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

    if (selectionConfig.isApplyMCGenTVX && (mcCollision.multMCFT0C() <= 0 || mcCollision.multMCFT0A() <= 0)) {
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
    centrality = -1.f;

    for (auto const& collision : recCollisions) {

      if (!selectionEvent(collision, false)) {
        continue;
      }

      if (selectCentEstimator == kFT0A) {
        centrality = collision.centFT0A();
      } else if (selectCentEstimator == kFT0C) {
        centrality = collision.centFT0C();
      } else if (selectCentEstimator == kFV0A) {
        centrality = collision.centFV0A();
      } else {
        centrality = collision.centFT0M(); // default includes kFT0M
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

      if ((mcPart.y() < selectionConfig.motherRapidityMin || mcPart.y() > selectionConfig.motherRapidityMax) || std::abs(mcPart.pdgCode()) != o2::constants::physics::kK0Star892) {
        continue;
      }

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

    if (selectCentEstimator == kFT0A) {
      centrality = collision.centFT0A();
    } else if (selectCentEstimator == kFT0C) {
      centrality = collision.centFT0C();
    } else if (selectCentEstimator == kFV0A) {
      centrality = collision.centFV0A();
    } else {
      centrality = collision.centFT0M(); // default includes kFT0M
    }

    if (!selectionEvent(collision, false)) {
      return;
    }

    for (const auto& [track1, track2] : combinations(CombinationsFullIndexPolicy(tracks, tracks))) {

      if (!selectionTrack(track1) || !selectionTrack(track2)) {
        continue;
      }

      if (track1.index() >= track2.index()) {
        continue;
      }

      if (track1.sign() * track2.sign() >= 0) {
        continue;
      }

      if (!track1.has_mcParticle() || !track2.has_mcParticle()) {
        continue;
      }

      const auto mc1 = track1.mcParticle();
      const auto mc2 = track2.mcParticle();

      if (!mc1.isPhysicalPrimary() || !mc2.isPhysicalPrimary()) {
        continue;
      }

      int pdg1 = std::abs(mc1.pdgCode());
      int pdg2 = std::abs(mc2.pdgCode());

      bool ok1 = (pdg1 == PDG_t::kPiPlus || pdg1 == PDG_t::kKPlus);
      bool ok2 = (pdg2 == PDG_t::kPiPlus || pdg2 == PDG_t::kKPlus);
      if (!ok1 || !ok2) {
        continue;
      }

      // pi-pi misidentification
      if (pdg1 == PDG_t::kPiPlus && pdg2 == PDG_t::kPiPlus) {
        ROOT::Math::PxPyPzMVector p1Fake(track1.px(), track1.py(), track1.pz(), massKa);
        ROOT::Math::PxPyPzMVector p2True(track2.px(), track2.py(), track2.pz(), massPi);

        auto misIDMother = p1Fake + p2True;
        if (misIDMother.Rapidity() > selectionConfig.motherRapidityMin && misIDMother.Rapidity() < selectionConfig.motherRapidityMax) {
          hMC.fill(HIST("RecMisID/hMassMisIDPiPi"), misIDMother.Pt(), centrality, misIDMother.M());
        }

        ROOT::Math::PxPyPzMVector p1True(track1.px(), track1.py(), track1.pz(), massPi);
        ROOT::Math::PxPyPzMVector p2Fake(track2.px(), track2.py(), track2.pz(), massKa);

        misIDMother = p1True + p2Fake;
        if (misIDMother.Rapidity() > selectionConfig.motherRapidityMin && misIDMother.Rapidity() < selectionConfig.motherRapidityMax) {
          hMC.fill(HIST("RecMisID/hMassMisIDPiPi"), misIDMother.Pt(), centrality, misIDMother.M());
        }
      }

      // KK misidentification
      if (pdg1 == PDG_t::kKPlus && pdg2 == PDG_t::kKPlus) {
        ROOT::Math::PxPyPzMVector p1Fake(track1.px(), track1.py(), track1.pz(), massPi);
        ROOT::Math::PxPyPzMVector p2True(track2.px(), track2.py(), track2.pz(), massKa);

        auto misIDMother = p1Fake + p2True;
        if (misIDMother.Rapidity() > selectionConfig.motherRapidityMin && misIDMother.Rapidity() < selectionConfig.motherRapidityMax) {
          hMC.fill(HIST("RecMisID/hMassMisIDKK"), misIDMother.Pt(), centrality, misIDMother.M());
        }

        ROOT::Math::PxPyPzMVector p1True(track1.px(), track1.py(), track1.pz(), massKa);
        ROOT::Math::PxPyPzMVector p2Fake(track2.px(), track2.py(), track2.pz(), massPi);

        misIDMother = p1True + p2Fake;
        if (misIDMother.Rapidity() > selectionConfig.motherRapidityMin && misIDMother.Rapidity() < selectionConfig.motherRapidityMax) {
          hMC.fill(HIST("RecMisID/hMassMisIDKK"), misIDMother.Pt(), centrality, misIDMother.M());
        }
      }
    }
  }
  PROCESS_SWITCH(Kstar892LightIon, processRecMisID, "Process Reconstructed MisID Background", false);

  void processRecReflection(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&, EventMCGenerated const&)
  {
    if (!collision.has_mcCollision()) {
      return;
    }

    if (!selectionEvent(collision, false)) {
      return;
    }

    centrality = -1.f;

    if (selectCentEstimator == kFT0A) {
      centrality = collision.centFT0A();
    } else if (selectCentEstimator == kFT0C) {
      centrality = collision.centFT0C();
    } else if (selectCentEstimator == kFV0A) {
      centrality = collision.centFV0A();
    } else {
      centrality = collision.centFT0M(); // default includes kFT0M
    }

    for (const auto& [track1, track2] : combinations(CombinationsFullIndexPolicy(tracks, tracks))) {

      if (!selectionTrack(track1) || !selectionTrack(track2)) {
        continue;
      }

      if (track1.index() >= track2.index()) {
        continue;
      }

      if (track1.sign() * track2.sign() >= 0) {
        continue;
      }

      if (!track1.has_mcParticle() || !track2.has_mcParticle()) {
        continue;
      }

      const auto mc1 = track1.mcParticle();
      const auto mc2 = track2.mcParticle();

      if (!mc1.isPhysicalPrimary() || !mc2.isPhysicalPrimary()) {
        continue;
      }

      bool sameMother = false;
      int motherPDG = 0;

      for (const auto& m1 : mc1.mothers_as<aod::McParticles>()) {
        for (const auto& m2 : mc2.mothers_as<aod::McParticles>()) {
          if (m1.globalIndex() == m2.globalIndex()) {
            sameMother = true;
            motherPDG = std::abs(m1.pdgCode());
            break;
          }
        }
        if (sameMother) {
          break;
        }
      }

      if (!sameMother) {
        continue;
      }

      int pdg1 = std::abs(mc1.pdgCode());
      int pdg2 = std::abs(mc2.pdgCode());

      // =====================================================
      // Rho0 (770) -> pi pi -> K pi
      // =====================================================
      if (motherPDG == PDG_t::kRho770_0 && pdg1 == PDG_t::kPiPlus && pdg2 == PDG_t::kPiPlus) {

        // track 1 -> K
        ROOT::Math::PxPyPzMVector p1K(track1.px(), track1.py(), track1.pz(), massKa);
        ROOT::Math::PxPyPzMVector p2Pi(track2.px(), track2.py(), track2.pz(), massPi);
        auto fake1 = p1K + p2Pi;

        if (fake1.Rapidity() > selectionConfig.motherRapidityMin && fake1.Rapidity() < selectionConfig.motherRapidityMax) {
          hMC.fill(HIST("Reflections/hRhoToKpi"), fake1.Pt(), centrality, fake1.M());
        }

        // track 2 -> K
        ROOT::Math::PxPyPzMVector p1Pi(track1.px(), track1.py(), track1.pz(), massPi);
        ROOT::Math::PxPyPzMVector p2K(track2.px(), track2.py(), track2.pz(), massKa);
        auto fake2 = p1Pi + p2K;

        if (fake2.Rapidity() > selectionConfig.motherRapidityMin && fake2.Rapidity() < selectionConfig.motherRapidityMax) {
          hMC.fill(HIST("Reflections/hRhoToKpi"), fake2.Pt(), centrality, fake2.M());
        }
      }

      // =====================================================
      // Omega (782) -> pi pi(pi0) -> K pi
      // =====================================================
      if (motherPDG == o2::constants::physics::kOmega && pdg1 == PDG_t::kPiPlus && pdg2 == PDG_t::kPiPlus) {

        // track 1 -> K
        ROOT::Math::PxPyPzMVector p1K(track1.px(), track1.py(), track1.pz(), massKa);
        ROOT::Math::PxPyPzMVector p2Pi(track2.px(), track2.py(), track2.pz(), massPi);
        auto fake1 = p1K + p2Pi;

        if (fake1.Rapidity() > selectionConfig.motherRapidityMin && fake1.Rapidity() < selectionConfig.motherRapidityMax) {
          hMC.fill(HIST("Reflections/hOmegaToKpi"), fake1.Pt(), centrality, fake1.M());
        }

        // track 2 -> K
        ROOT::Math::PxPyPzMVector p1Pi(track1.px(), track1.py(), track1.pz(), massPi);
        ROOT::Math::PxPyPzMVector p2K(track2.px(), track2.py(), track2.pz(), massKa);
        auto fake2 = p1Pi + p2K;

        if (fake2.Rapidity() > selectionConfig.motherRapidityMin && fake2.Rapidity() < selectionConfig.motherRapidityMax) {
          hMC.fill(HIST("Reflections/hOmegaToKpi"), fake2.Pt(), centrality, fake2.M());
        }
      }

      // =====================================================
      //  Eta -> pi pi(pi0) -> K pi
      // =====================================================
      if (motherPDG == o2::constants::physics::kEta && pdg1 == PDG_t::kPiPlus && pdg2 == PDG_t::kPiPlus) {

        // track 1 -> K
        ROOT::Math::PxPyPzMVector p1K(track1.px(), track1.py(), track1.pz(), massKa);
        ROOT::Math::PxPyPzMVector p2Pi(track2.px(), track2.py(), track2.pz(), massPi);
        auto fake1 = p1K + p2Pi;

        if (fake1.Rapidity() > selectionConfig.motherRapidityMin && fake1.Rapidity() < selectionConfig.motherRapidityMax) {
          hMC.fill(HIST("Reflections/hEtaToKpi"), fake1.Pt(), centrality, fake1.M());
        }

        // track 2 -> K
        ROOT::Math::PxPyPzMVector p1Pi(track1.px(), track1.py(), track1.pz(), massPi);
        ROOT::Math::PxPyPzMVector p2K(track2.px(), track2.py(), track2.pz(), massKa);
        auto fake2 = p1Pi + p2K;

        if (fake2.Rapidity() > selectionConfig.motherRapidityMin && fake2.Rapidity() < selectionConfig.motherRapidityMax) {
          hMC.fill(HIST("Reflections/hEtaToKpi"), fake2.Pt(), centrality, fake2.M());
        }
      }

      // =====================================================
      //  Eta' (958) -> pi pi(eta) -> K pi
      // =====================================================
      if (motherPDG == o2::constants::physics::kEtaPrime && pdg1 == PDG_t::kPiPlus && pdg2 == PDG_t::kPiPlus) {

        // track 1 -> K
        ROOT::Math::PxPyPzMVector p1K(track1.px(), track1.py(), track1.pz(), massKa);
        ROOT::Math::PxPyPzMVector p2Pi(track2.px(), track2.py(), track2.pz(), massPi);
        auto fake1 = p1K + p2Pi;

        if (fake1.Rapidity() > selectionConfig.motherRapidityMin && fake1.Rapidity() < selectionConfig.motherRapidityMax) {
          hMC.fill(HIST("Reflections/hEtaPrimeToKpi"), fake1.Pt(), centrality, fake1.M());
        }

        // track 2 -> K
        ROOT::Math::PxPyPzMVector p1Pi(track1.px(), track1.py(), track1.pz(), massPi);
        ROOT::Math::PxPyPzMVector p2K(track2.px(), track2.py(), track2.pz(), massKa);
        auto fake2 = p1Pi + p2K;

        if (fake2.Rapidity() > selectionConfig.motherRapidityMin && fake2.Rapidity() < selectionConfig.motherRapidityMax) {
          hMC.fill(HIST("Reflections/hEtaPrimeToKpi"), fake2.Pt(), centrality, fake2.M());
        }
      }

      // =====================================================
      // Phi (1020) -> KK -> K pi
      // =====================================================
      if (motherPDG == o2::constants::physics::kPhi && pdg1 == PDG_t::kKPlus && pdg2 == PDG_t::kKPlus) {

        // track 1 -> pi
        ROOT::Math::PxPyPzMVector p1Pi(track1.px(), track1.py(), track1.pz(), massPi);
        ROOT::Math::PxPyPzMVector p2K(track2.px(), track2.py(), track2.pz(), massKa);
        auto fake1 = p1Pi + p2K;

        if (fake1.Rapidity() > selectionConfig.motherRapidityMin && fake1.Rapidity() < selectionConfig.motherRapidityMax) {
          hMC.fill(HIST("Reflections/hPhiToKpi"), fake1.Pt(), centrality, fake1.M());
        }

        // track 2 -> pi
        ROOT::Math::PxPyPzMVector p1K(track1.px(), track1.py(), track1.pz(), massKa);
        ROOT::Math::PxPyPzMVector p2Pi(track2.px(), track2.py(), track2.pz(), massPi);
        auto fake2 = p1K + p2Pi;

        if (fake2.Rapidity() > selectionConfig.motherRapidityMin && fake2.Rapidity() < selectionConfig.motherRapidityMax) {
          hMC.fill(HIST("Reflections/hPhiToKpi"), fake2.Pt(), centrality, fake2.M());
        }
      }

      // =====================================================
      // K*0 Self-reflection
      // =====================================================
      if (motherPDG == o2::constants::physics::kK0Star892 && ((pdg1 == PDG_t::kPiPlus && pdg2 == PDG_t::kKPlus) || (pdg1 == PDG_t::kKPlus && pdg2 == PDG_t::kPiPlus))) {

        ROOT::Math::PxPyPzMVector p1Swap(track1.px(), track1.py(), track1.pz(), pdg1 == PDG_t::kKPlus ? massPi : massKa);

        ROOT::Math::PxPyPzMVector p2Swap(track2.px(), track2.py(), track2.pz(), pdg2 == PDG_t::kKPlus ? massPi : massKa);

        auto fake = p1Swap + p2Swap;

        if (fake.Rapidity() > selectionConfig.motherRapidityMin && fake.Rapidity() < selectionConfig.motherRapidityMax) {
          hMC.fill(HIST("Reflections/hKstarSelf"), fake.Pt(), centrality, fake.M());
        }
      }

      // =========================================
      // Higher resonance feed-down
      // =========================================

      if ((pdg1 == PDG_t::kKPlus && pdg2 == PDG_t::kPiPlus) || (pdg1 == PDG_t::kPiPlus && pdg2 == PDG_t::kKPlus)) {

        ROOT::Math::PxPyPzMVector p1(track1.px(), track1.py(), track1.pz(), pdg1 == PDG_t::kKPlus ? massKa : massPi);

        ROOT::Math::PxPyPzMVector p2(track2.px(), track2.py(), track2.pz(), pdg2 == PDG_t::kKPlus ? massKa : massPi);

        auto pair = p1 + p2;

        if (pair.Rapidity() < selectionConfig.motherRapidityMin || pair.Rapidity() > selectionConfig.motherRapidityMax) {
          continue;
        }

        // ==========================================
        // K1(1270)
        // ==========================================
        if (motherPDG == o2::constants::physics::kK1_1270Plus) {
          hMC.fill(HIST("FeedDown/hK1_1270"), pair.Pt(), centrality, pair.M());
        }

        // ==========================================
        // K1(1400)
        // ==========================================
        if (motherPDG == kK11400Plus) {
          hMC.fill(HIST("FeedDown/hK1_1400"), pair.Pt(), centrality, pair.M());
        }

        // ==========================================
        // K*(1410)
        // ==========================================
        if (motherPDG == kKstar14100) {
          hMC.fill(HIST("FeedDown/hKstar1410"), pair.Pt(), centrality, pair.M());
        }

        // ==========================================
        // K*0(1430)
        // ==========================================
        if (motherPDG == kKstar14300) {
          hMC.fill(HIST("FeedDown/hKstar0_1430_0"), pair.Pt(), centrality, pair.M());
        }

        // ==========================================
        // K*±(1430)
        // ==========================================
        if (motherPDG == kKstar1430Plus) {
          hMC.fill(HIST("FeedDown/hKstar0_1430_ch"), pair.Pt(), centrality, pair.M());
        }

        // ==========================================
        // K*2(1430)
        // ==========================================
        if (motherPDG == kKstar214300) {
          hMC.fill(HIST("FeedDown/hKstar2_1430"), pair.Pt(), centrality, pair.M());
        }

        // ==========================================
        // K2(1770)
        // ==========================================
        if (motherPDG == kK217700) {
          hMC.fill(HIST("FeedDown/hK2_1770"), pair.Pt(), centrality, pair.M());
        }

        // ==========================================
        // K2(1820)
        // ==========================================
        if (motherPDG == kK218200) {
          hMC.fill(HIST("FeedDown/hK2_1820"), pair.Pt(), centrality, pair.M());
        }

        // ==========================================
        // K*(1680)
        // ==========================================
        if (motherPDG == kKstar16800) {
          hMC.fill(HIST("FeedDown/hKstar1680"), pair.Pt(), centrality, pair.M());
        }
      }
    }
  }
  PROCESS_SWITCH(Kstar892LightIon, processRecReflection, "Process reconstructed reflections", false);

  void processRecCorrelatedBackground(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&, EventMCGenerated const&)
  {
    if (!collision.has_mcCollision()) {
      return;
    }

    if (!selectionEvent(collision, false)) {
      return;
    }

    centrality = -1.f;

    if (selectCentEstimator == kFT0A) {
      centrality = collision.centFT0A();
    } else if (selectCentEstimator == kFT0C) {
      centrality = collision.centFT0C();
    } else if (selectCentEstimator == kFV0A) {
      centrality = collision.centFV0A();
    } else {
      centrality = collision.centFT0M(); // default includes kFT0M
    }

    for (const auto& [track1, track2] : combinations(CombinationsFullIndexPolicy(tracks, tracks))) {

      if (!selectionTrack(track1)) {
        continue;
      }

      if (!selectionTrack(track2)) {
        continue;
      }

      if (track1.index() >= track2.index()) {
        continue;
      }

      if (track1.sign() * track2.sign() >= 0) {
        continue;
      }

      if (!track1.has_mcParticle()) {
        continue;
      }

      if (!track2.has_mcParticle()) {
        continue;
      }

      auto mc1 = track1.mcParticle();
      auto mc2 = track2.mcParticle();

      // if (!mc1.isPhysicalPrimary() || !mc2.isPhysicalPrimary())
      // continue;

      int pdg1 = std::abs(mc1.pdgCode());
      int pdg2 = std::abs(mc2.pdgCode());

      if ((pdg1 != PDG_t::kKPlus || pdg2 != PDG_t::kPiPlus) && (pdg1 != PDG_t::kPiPlus || pdg2 != PDG_t::kKPlus)) {
        continue;
      }

      ROOT::Math::PxPyPzMVector p1(track1.px(), track1.py(), track1.pz(), pdg1 == PDG_t::kKPlus ? massKa : massPi);
      ROOT::Math::PxPyPzMVector p2(track2.px(), track2.py(), track2.pz(), pdg2 == PDG_t::kKPlus ? massKa : massPi);
      auto pair = p1 + p2;

      if (pair.Rapidity() < selectionConfig.motherRapidityMin || pair.Rapidity() > selectionConfig.motherRapidityMax) {
        continue;
      }

      bool trueKstar = false;

      for (const auto& m1 : mc1.mothers_as<aod::McParticles>()) {
        for (const auto& m2 : mc2.mothers_as<aod::McParticles>()) {

          if (m1.globalIndex() == m2.globalIndex() && std::abs(m1.pdgCode()) == o2::constants::physics::kK0Star892) {
            trueKstar = true;
            break;
          }
        }
        if (trueKstar) {
          break;
        }
      }
      if (trueKstar) {
        continue;
      }

      auto ancestors1 = getAncestors(mc1);
      auto ancestors2 = getAncestors(mc2);

      bool fromK11270 = false;
      bool fromK11400 = false;
      bool fromKstar1680 = false;
      bool fromK21770 = false;
      bool fromK21820 = false;
      bool fromKstar1410 = false;
      bool fromKstar01430 = false;

      for (auto const& a1 : ancestors1) {

        for (auto const& a2 : ancestors2) {

          if (a1.first != a2.first) {
            continue;
          }

          int ancestorPdg = a1.second;

          if (ancestorPdg == o2::constants::physics::kK1_1270Plus) {
            fromK11270 = true;
          }

          if (ancestorPdg == kK11400Plus) {
            fromK11400 = true;
          }

          if (ancestorPdg == kKStar1680Plus) {
            fromKstar1680 = true;
          }

          if (ancestorPdg == kK21770Plus) {
            fromK21770 = true;
          }

          if (ancestorPdg == kK21820Plus) {
            fromK21820 = true;
          }

          if (ancestorPdg == kKstar14100) {
            fromKstar1410 = true;
          }

          if (ancestorPdg == kKstar14300) {
            fromKstar01430 = true;
          }
        }
      }

      if (fromKstar1410) {
        hMC.fill(HIST("CorrelatedBG/hKstar1410"), pair.Pt(), centrality, pair.M());
      }

      if (fromKstar01430) {
        hMC.fill(HIST("CorrelatedBG/hKstar0_1430"), pair.Pt(), centrality, pair.M());
      }

      if (fromK11270) {
        hMC.fill(HIST("CorrelatedBG/hK1_1270"), pair.Pt(), centrality, pair.M());
      }

      if (fromK11400) {
        hMC.fill(HIST("CorrelatedBG/hK1_1400"), pair.Pt(), centrality, pair.M());
      }

      if (fromKstar1680) {
        hMC.fill(HIST("CorrelatedBG/hKstar1680"), pair.Pt(), centrality, pair.M());
      }

      if (fromK21770) {
        hMC.fill(HIST("CorrelatedBG/hK2_1770"), pair.Pt(), centrality, pair.M());
      }

      if (fromK21820) {
        hMC.fill(HIST("CorrelatedBG/hK2_1820"), pair.Pt(), centrality, pair.M());
      }
    }
  }

  PROCESS_SWITCH(Kstar892LightIon, processRecCorrelatedBackground, "Process correlated background", false);

  void processTemplateMC(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&, EventMCGenerated const&)
  {
    if (!collision.has_mcCollision()) {
      return;
    }

    if (!selectionEvent(collision, false)) {
      return;
    }

    centrality = -1.f;
    if (selectCentEstimator == kFT0A) {
      centrality = collision.centFT0A();
    } else if (selectCentEstimator == kFT0C) {
      centrality = collision.centFT0C();
    } else if (selectCentEstimator == kFV0A) {
      centrality = collision.centFV0A();
    } else {
      centrality = collision.centFT0M(); // default includes kFT0M
    }

    auto classifyOrigin = [](const aod::McParticle& part, int64_t& motherIdx) -> int {
      motherIdx = -1;
      if (!part.has_mothers()) {
        return 0;
      }
      for (const auto& mom : part.mothers_as<aod::McParticles>()) {
        if (!mom.producedByGenerator()) {
          continue;
        }
        motherIdx = mom.globalIndex();
        if (std::abs(mom.pdgCode()) == o2::constants::physics::kK0Star892) {
          return 1;
        }
        return 2;
      }
      return 0;
    };

    for (const auto& [track1, track2] : combinations(CombinationsFullIndexPolicy(tracks, tracks))) {

      if (!selectionTrack(track1) || !selectionTrack(track2)) {
        continue;
      }

      if (track1.globalIndex() == track2.globalIndex()) {
        continue;
      }

      if (!selectionPair(track1, track2)) {
        continue;
      }

      if (!selectionPID(track1, PIDParticle::kKaon)) {
        continue;
      }
      if (!selectionPID(track2, PIDParticle::kPion)) {
        continue;
      }

      if (isMisidentified(track1, PIDParticle::kKaon)) {
        continue;
      }
      if (isMisidentified(track2, PIDParticle::kPion)) {
        continue;
      }

      if (selectionConfig.isApplyFakeTrack && (isFakeTrack(track1, PIDParticle::kKaon) || isFakeTrack(track2, PIDParticle::kPion))) {
        continue;
      }

      if (!track1.has_mcParticle() || !track2.has_mcParticle()) {
        continue;
      }

      daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
      daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
      mother = daughter1 + daughter2;

      if (mother.Rapidity() < selectionConfig.motherRapidityMin || mother.Rapidity() > selectionConfig.motherRapidityMax) {
        continue;
      }

      recMass = mother.M();
      recPt = mother.Pt();

      const auto mcKa = track1.mcParticle();
      const auto mcPi = track2.mcParticle();

      if (selectionConfig.isPrimaryTrack && (!mcKa.isPhysicalPrimary() || !mcPi.isPhysicalPrimary())) {
        continue;
      }

      int64_t kaMotherIdx = -1;
      int64_t piMotherIdx = -1;
      const int kaOrigin = classifyOrigin(mcKa, kaMotherIdx);
      const int piOrigin = classifyOrigin(mcPi, piMotherIdx);

      const bool kaIsKstar = (kaOrigin == 1);
      const bool piIsKstar = (piOrigin == 1);
      const bool kaHasParent = (kaOrigin != 0);
      const bool piHasParent = (piOrigin != 0);
      const bool sameMother = (kaMotherIdx >= 0) && (piMotherIdx >= 0) && (kaMotherIdx == piMotherIdx);

      if (kaIsKstar && piIsKstar && sameMother) {
        const bool correctAssignment = (std::abs(mcKa.pdgCode()) == PDG_t::kKPlus) && (std::abs(mcPi.pdgCode()) == PDG_t::kPiPlus);
        if (correctAssignment) {
          hMC.fill(HIST("Template/hSignal"), recPt, centrality, recMass);
        } else {
          hMC.fill(HIST("Template/hKstarReflection"), recPt, centrality, recMass);
        }

      } else if (sameMother && kaHasParent && piHasParent) {
        hMC.fill(HIST("Template/hSameMotherOther"), recPt, centrality, recMass);

      } else if (!sameMother && kaHasParent && piHasParent) {
        hMC.fill(HIST("Template/hDifferentMother"), recPt, centrality, recMass);

      } else if (kaHasParent != piHasParent) {
        hMC.fill(HIST("Template/hPrimaryResonance"), recPt, centrality, recMass);

      } else {
        hMC.fill(HIST("Template/hPrimaryPrimary"), recPt, centrality, recMass);
      }
    }
  }
  PROCESS_SWITCH(Kstar892LightIon, processTemplateMC, "Process MC Template Background", false);

  void processRecKinematics(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&, EventMCGenerated const&)
  {
    if (!collision.has_mcCollision()) {
      return;
    }

    centrality = collision.centFT0M();

    if (selectCentEstimator == kFT0A) {
      centrality = collision.centFT0A();
    } else if (selectCentEstimator == kFT0C) {
      centrality = collision.centFT0C();
    } else if (selectCentEstimator == kFV0A) {
      centrality = collision.centFV0A();
    } else {
      centrality = collision.centFT0M(); // default includes kFT0M
    }

    if (!selectionEvent(collision, true)) { // fill MC event cut histogram
      return;
    }

    hMC.fill(HIST("Kinematics/h1RecCent"), centrality);

    // auto oldindex = -999;

    for (const auto& [track1, track2] : combinations(CombinationsFullIndexPolicy(tracks, tracks))) {
      if (!selectionTrack(track1) || !selectionTrack(track2)) {
        continue;
      }

      if (!track1.has_mcParticle() || !track2.has_mcParticle()) {
        continue;
      }

      if (track1.index() <= track2.index()) {
        continue;
      }

      if (!selectionPair(track1, track2)) {
        continue;
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

      if ((track1PDG != PDG_t::kPiPlus || track2PDG != PDG_t::kKPlus) && (track1PDG != PDG_t::kKPlus || track2PDG != PDG_t::kPiPlus)) {
        continue;
      }

      ROOT::Math::PxPyPzMVector kVec;
      ROOT::Math::PxPyPzMVector piVec;

      double phiK = 0.0, phiPi = 0.0;
      double etaK = 0.0, etaPi = 0.0;
      double ptK = 0.0, ptPi = 0.0;

      if (track1PDG == PDG_t::kKPlus) {

        if (!selectionPID(track1, PIDParticle::kKaon) || !selectionPID(track2, PIDParticle::kPion)) { // Treat track1 as the kaon candidate and track2 as the pion candidate
          continue;
        }

        if (isMisidentified(track1, PIDParticle::kKaon) || isMisidentified(track2, PIDParticle::kPion)) {
          continue;
        }

        if (selectionConfig.isApplyFakeTrack && (isFakeTrack(track1, PIDParticle::kKaon) || isFakeTrack(track2, PIDParticle::kPion))) {
          continue;
        }

        kVec = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        piVec = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);

        phiK = track1.phi();
        etaK = track1.eta();
        ptK = track1.pt();

        phiPi = track2.phi();
        etaPi = track2.eta();
        ptPi = track2.pt();

      } else {

        if (!selectionPID(track1, PIDParticle::kPion) || !selectionPID(track2, PIDParticle::kKaon)) { // Treat track1 as the pion candidate and track2 as the kaon candidate
          continue;
        }

        if (isMisidentified(track1, PIDParticle::kPion) || isMisidentified(track2, PIDParticle::kKaon)) {
          continue;
        }

        if (selectionConfig.isApplyFakeTrack && (isFakeTrack(track1, PIDParticle::kPion) || isFakeTrack(track2, PIDParticle::kKaon))) {
          continue;
        }

        kVec = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        piVec = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massPi);

        phiK = track2.phi();
        etaK = track2.eta();
        ptK = track2.pt();

        phiPi = track1.phi();
        etaPi = track1.eta();
        ptPi = track1.pt();
      }

      auto pair = kVec + piVec;

      double mass = pair.M();

      double cosAngle = (kVec.Px() * piVec.Px() + kVec.Py() * piVec.Py() + kVec.Pz() * piVec.Pz()) / (kVec.P() * piVec.P());

      cosAngle = std::clamp(cosAngle, -1.0, 1.0);

      double openAngle = std::acos(cosAngle);

      double dPhi = RecoDecay::constrainAngle(phiK - phiPi, -o2::constants::math::PI);

      double dEta = etaK - etaPi;

      double dR = std::sqrt(dPhi * dPhi + dEta * dEta);

      bool isTrueKstar = false;
      bool hasCommonMother = false;

      int commonMotherPDG = initialValue;

      for (const auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
        for (const auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {

          if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
            continue;
          }

          if (mothertrack1.globalIndex() != mothertrack2.globalIndex()) {
            continue;
          }

          hasCommonMother = true;
          commonMotherPDG = mothertrack1.pdgCode();

          if (std::abs(commonMotherPDG) == o2::constants::physics::kK0Star892) {
            isTrueKstar = true;
          }
        }
      }

      if (pair.Rapidity() < selectionConfig.motherRapidityMin || pair.Rapidity() > selectionConfig.motherRapidityMax) {
        continue;
      }

      if (isTrueKstar) {

        hMC.fill(HIST("Kinematics/hPDGMotherKstar"), commonMotherPDG, mass);

        hMC.fill(HIST("Kinematics/hOpenAngleKstar"), openAngle, mass);

        hMC.fill(HIST("Kinematics/hDeltaPhiKstar"), dPhi, mass);

        hMC.fill(HIST("Kinematics/hDeltaEtaKstar"), dEta, mass);

        hMC.fill(HIST("Kinematics/hDeltaRKstar"), dR, mass);

        hMC.fill(HIST("Kinematics/hKaonPtKstar"), ptK, mass);

        hMC.fill(HIST("Kinematics/hPionPtKstar"), ptPi, mass);

        hMC.fill(HIST("Kinematics/hPtCentMassKstar"), pair.Pt(), centrality, mass);

      } else {

        hMC.fill(HIST("Kinematics/hOpenAngleOther"), openAngle, mass);

        hMC.fill(HIST("Kinematics/hDeltaPhiOther"), dPhi, mass);

        hMC.fill(HIST("Kinematics/hDeltaEtaOther"), dEta, mass);

        hMC.fill(HIST("Kinematics/hDeltaROther"), dR, mass);

        hMC.fill(HIST("Kinematics/hKaonPtOther"), ptK, mass);

        hMC.fill(HIST("Kinematics/hPionPtOther"), ptPi, mass);

        hMC.fill(HIST("Kinematics/hPtCentMassOther"), pair.Pt(), centrality, mass);

        if (hasCommonMother) {
          hMC.fill(HIST("Kinematics/hOtherMotherPDG"), commonMotherPDG, mass);
        }
      }
    }
  }
  PROCESS_SWITCH(Kstar892LightIon, processRecKinematics, "Process Reconstructed Kinematics", false);

  void processMisIdKinematics(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&, EventMCGenerated const&)
  {
    if (!collision.has_mcCollision()) {
      return;
    }

    if (!selectionEvent(collision, true)) {
      return;
    }

    centrality = -1;
    if (selectCentEstimator == kFT0A) {
      centrality = collision.centFT0A();
    } else if (selectCentEstimator == kFT0C) {
      centrality = collision.centFT0C();
    } else if (selectCentEstimator == kFV0A) {
      centrality = collision.centFV0A();
    } else {
      centrality = collision.centFT0M();
    }

    hMC.fill(HIST("KinematicsMisId/hRecCent"), centrality);

    for (const auto& [track1, track2] : combinations(CombinationsFullIndexPolicy(tracks, tracks))) {

      if (!selectionTrack(track1) || !selectionTrack(track2)) {
        continue;
      }
      if (track1.globalIndex() == track2.globalIndex()) {
        continue;
      }
      if (!selectionPair(track1, track2)) {
        continue;
      }

      if (!selectionPID(track1, PIDParticle::kKaon)) {
        continue;
      }
      if (!selectionPID(track2, PIDParticle::kPion)) {
        continue;
      }

      if (isMisidentified(track1, PIDParticle::kKaon)) {
        continue;
      }
      if (isMisidentified(track2, PIDParticle::kPion)) {
        continue;
      }

      if (selectionConfig.isApplyFakeTrack && (isFakeTrack(track1, PIDParticle::kKaon) || isFakeTrack(track2, PIDParticle::kPion))) {
        continue;
      }

      daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
      daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
      mother = daughter1 + daughter2;

      double ptK = daughter1.Pt();
      double ptPi = daughter2.Pt();
      double ptPair = mother.Pt();
      double mass = mother.M();

      double cosAngle = (daughter1.Px() * daughter2.Px() + daughter1.Py() * daughter2.Py() + daughter1.Pz() * daughter2.Pz()) / (daughter1.P() * daughter2.P());
      cosAngle = std::clamp(cosAngle, -1.0, 1.0);
      double openAngle = std::acos(cosAngle);

      double dPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -o2::constants::math::PI);
      double dEta = track1.eta() - track2.eta();
      double dR = std::sqrt(dPhi * dPhi + dEta * dEta);

      const auto mctrack1 = track1.mcParticle();
      const auto mctrack2 = track2.mcParticle();

      if (!mctrack1.isPhysicalPrimary() || !mctrack2.isPhysicalPrimary()) {
        continue;
      }

      int truePdg1 = std::abs(mctrack1.pdgCode());
      int truePdg2 = std::abs(mctrack2.pdgCode());

      bool isFakeKaon = (truePdg1 != PDG_t::kKPlus);
      bool isFakePion = (truePdg2 != PDG_t::kPiPlus);

      bool hasCommonMother = false;
      int commonMotherPDG = initialValue;

      for (const auto& mother1 : mctrack1.mothers_as<aod::McParticles>()) {
        for (const auto& mother2 : mctrack2.mothers_as<aod::McParticles>()) {
          if (mother1.pdgCode() == mother2.pdgCode() && mother1.globalIndex() == mother2.globalIndex()) {
            hasCommonMother = true;
            commonMotherPDG = std::abs(mother1.pdgCode());
            break;
          }
        }
        if (hasCommonMother) {
          break;
        }
      }

      if (!isFakeKaon && !isFakePion && hasCommonMother && commonMotherPDG == o2::constants::physics::kK0Star892) { // True Signal (Correct identification and from a K*0)
        if (additionalKin) {
          hMC.fill(HIST("KinematicsMisId/hTrueKstar_PtOpenAngleMass"), ptPair, openAngle, mass);
          hMC.fill(HIST("KinematicsMisId/hTrueKstar_PtDeltaPhiMass"), ptPair, dPhi, mass);
          hMC.fill(HIST("KinematicsMisId/hTrueKstar_PtDeltaEtaMass"), ptPair, dEta, mass);
          hMC.fill(HIST("KinematicsMisId/hTrueKstar_PtKaonPtMass"), ptPair, ptK, mass);
          hMC.fill(HIST("KinematicsMisId/hTrueKstar_PtPionPtMass"), ptPair, ptPi, mass);
        }
        hMC.fill(HIST("KinematicsMisId/hTrueKstar_PtDeltaRMass"), ptPair, dR, mass);
        hMC.fill(HIST("KinematicsMisId/hTrueKstar_PtCentMass"), ptPair, centrality, mass);
      } else if (isFakeKaon || isFakePion) { // Correct identification but from a different mother
        if (hasCommonMother) {
          if (commonMotherPDG == o2::constants::physics::kOmega) { // Omega reflection
            if (additionalKin) {
              hMC.fill(HIST("KinematicsMisId/hOmega_PtOpenAngleMass"), ptPair, openAngle, mass);
              hMC.fill(HIST("KinematicsMisId/hOmega_PtDeltaPhiMass"), ptPair, dPhi, mass);
              hMC.fill(HIST("KinematicsMisId/hOmega_PtDeltaEtaMass"), ptPair, dEta, mass);
              hMC.fill(HIST("KinematicsMisId/hOmega_PtKaonPtMass"), ptPair, ptK, mass);
              hMC.fill(HIST("KinematicsMisId/hOmega_PtPionPtMass"), ptPair, ptPi, mass);
            }
            hMC.fill(HIST("KinematicsMisId/hOmega_PtDeltaRMass"), ptPair, dR, mass);
            hMC.fill(HIST("KinematicsMisId/hOmega_PtCentMass"), ptPair, centrality, mass);
          } else if (commonMotherPDG == PDG_t::kRho770_0) { // Rho reflection
            if (additionalKin) {
              hMC.fill(HIST("KinematicsMisId/hRho_PtOpenAngleMass"), ptPair, openAngle, mass);
              hMC.fill(HIST("KinematicsMisId/hRho_PtDeltaPhiMass"), ptPair, dPhi, mass);
              hMC.fill(HIST("KinematicsMisId/hRho_PtDeltaEtaMass"), ptPair, dEta, mass);
              hMC.fill(HIST("KinematicsMisId/hRho_PtKaonPtMass"), ptPair, ptK, mass);
              hMC.fill(HIST("KinematicsMisId/hRho_PtPionPtMass"), ptPair, ptPi, mass);
            }
            hMC.fill(HIST("KinematicsMisId/hRho_PtDeltaRMass"), ptPair, dR, mass);
            hMC.fill(HIST("KinematicsMisId/hRho_PtCentMass"), ptPair, centrality, mass);
          } else if (commonMotherPDG == o2::constants::physics::kEta) { // Eta reflection
            if (additionalKin) {
              hMC.fill(HIST("KinematicsMisId/hEta_PtOpenAngleMass"), ptPair, openAngle, mass);
              hMC.fill(HIST("KinematicsMisId/hEta_PtDeltaPhiMass"), ptPair, dPhi, mass);
              hMC.fill(HIST("KinematicsMisId/hEta_PtDeltaEtaMass"), ptPair, dEta, mass);
              hMC.fill(HIST("KinematicsMisId/hEta_PtKaonPtMass"), ptPair, ptK, mass);
              hMC.fill(HIST("KinematicsMisId/hEta_PtPionPtMass"), ptPair, ptPi, mass);
            }
            hMC.fill(HIST("KinematicsMisId/hEta_PtDeltaRMass"), ptPair, dR, mass);
            hMC.fill(HIST("KinematicsMisId/hEta_PtCentMass"), ptPair, centrality, mass);
          } else if (commonMotherPDG == o2::constants::physics::kEtaPrime) { // Eta' reflection
            if (additionalKin) {
              hMC.fill(HIST("KinematicsMisId/hEtaP_PtOpenAngleMass"), ptPair, openAngle, mass);
              hMC.fill(HIST("KinematicsMisId/hEtaP_PtDeltaPhiMass"), ptPair, dPhi, mass);
              hMC.fill(HIST("KinematicsMisId/hEtaP_PtDeltaEtaMass"), ptPair, dEta, mass);
              hMC.fill(HIST("KinematicsMisId/hEtaP_PtKaonPtMass"), ptPair, ptK, mass);
              hMC.fill(HIST("KinematicsMisId/hEtaP_PtPionPtMass"), ptPair, ptPi, mass);
            }
            hMC.fill(HIST("KinematicsMisId/hEtaP_PtDeltaRMass"), ptPair, dR, mass);
            hMC.fill(HIST("KinematicsMisId/hEtaP_PtCentMass"), ptPair, centrality, mass);
          } else if (commonMotherPDG == o2::constants::physics::kPhi) { // Phi reflection
            if (additionalKin) {
              hMC.fill(HIST("KinematicsMisId/hPhi_PtOpenAngleMass"), ptPair, openAngle, mass);
              hMC.fill(HIST("KinematicsMisId/hPhi_PtDeltaPhiMass"), ptPair, dPhi, mass);
              hMC.fill(HIST("KinematicsMisId/hPhi_PtDeltaEtaMass"), ptPair, dEta, mass);
              hMC.fill(HIST("KinematicsMisId/hPhi_PtKaonPtMass"), ptPair, ptK, mass);
              hMC.fill(HIST("KinematicsMisId/hPhi_PtPionPtMass"), ptPair, ptPi, mass);
            }
            hMC.fill(HIST("KinematicsMisId/hPhi_PtDeltaRMass"), ptPair, dR, mass);
            hMC.fill(HIST("KinematicsMisId/hPhi_PtCentMass"), ptPair, centrality, mass);
          } else if (commonMotherPDG == o2::constants::physics::kK0Star892) { // K*0 self reflection
            if (additionalKin) {
              hMC.fill(HIST("KinematicsMisId/hKstarMisId_PtOpenAngleMass"), ptPair, openAngle, mass);
              hMC.fill(HIST("KinematicsMisId/hKstarMisId_PtDeltaPhiMass"), ptPair, dPhi, mass);
              hMC.fill(HIST("KinematicsMisId/hKstarMisId_PtDeltaEtaMass"), ptPair, dEta, mass);
              hMC.fill(HIST("KinematicsMisId/hKstarMisId_PtKaonPtMass"), ptPair, ptK, mass);
              hMC.fill(HIST("KinematicsMisId/hKstarMisId_PtPionPtMass"), ptPair, ptPi, mass);
            }
            hMC.fill(HIST("KinematicsMisId/hKstarMisId_PtDeltaRMass"), ptPair, dR, mass);
            hMC.fill(HIST("KinematicsMisId/hKstarMisId_PtCentMass"), ptPair, centrality, mass);
          }
        } else { // Correct indentification but no common mother
          if (additionalKin) {
            hMC.fill(HIST("KinematicsMisId/hOtherMisId_PtOpenAngleMass"), ptPair, openAngle, mass);
            hMC.fill(HIST("KinematicsMisId/hOtherMisId_PtDeltaPhiMass"), ptPair, dPhi, mass);
            hMC.fill(HIST("KinematicsMisId/hOtherMisId_PtDeltaEtaMass"), ptPair, dEta, mass);
            hMC.fill(HIST("KinematicsMisId/hOtherMisId_PtKaonPtMass"), ptPair, ptK, mass);
            hMC.fill(HIST("KinematicsMisId/hOtherMisId_PtPionPtMass"), ptPair, ptPi, mass);
          }
          hMC.fill(HIST("KinematicsMisId/hOtherMisId_PtDeltaRMass"), ptPair, dR, mass);
          hMC.fill(HIST("KinematicsMisId/hOtherMisId_PtCentMass"), ptPair, centrality, mass);
        }
      }
    }
  }
  PROCESS_SWITCH(Kstar892LightIon, processMisIdKinematics, "Process Signal vs Misid Background Kinematics", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Kstar892LightIon>(cfgc)};
}
