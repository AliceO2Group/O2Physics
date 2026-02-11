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

/// \file kstarqa.cxx
/// \brief Code for Kstar resonance without resonance initializer
/// \author prottay das, sawan
/// \since 13/03/2024

// #include <TDatabasePDG.h>
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
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
#include <TVector2.h>

#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using std::array;
using namespace o2::aod::rctsel;

struct Kstarqa {

  SliceCache cache;

  struct : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCut;
  RCTFlagsChecker rctChecker;

  struct : ConfigurableGroup {
    // Configurables for event selections
    Configurable<bool> isINELgt0{"isINELgt0", true, "INEL>0 selection"};
    Configurable<bool> isSel8{"isSel8", false, "Event selection sel8"};
    Configurable<bool> isTriggerTVX{"isTriggerTVX", false, "TriggerTVX"};
    Configurable<bool> isGoodZvtxFT0vsPV{"isGoodZvtxFT0vsPV", false, "IsGoodZvtxFT0vsPV"};
    Configurable<bool> isApplyOccCut{"isApplyOccCut", true, "Apply occupancy cut"};
    Configurable<bool> isNoSameBunchPileup{"isNoSameBunchPileup", true, "kNoSameBunchPileup"};
    Configurable<bool> isGoodITSLayersAll{"isGoodITSLayersAll", true, "Require all ITS layers to be good"};
    Configurable<bool> isNoTimeFrameBorder{"isNoTimeFrameBorder", true, "kNoTimeFrameBorder"};
    Configurable<bool> isNoITSROFrameBorder{"isNoITSROFrameBorder", true, "kNoITSROFrameBorder"};
    Configurable<bool> isapplypTdepPID{"isapplypTdepPID", false, "Apply pT dependent PID"};
    Configurable<bool> isapplypTdepPIDTOF{"isapplypTdepPIDTOF", false, "Apply pT dependent PID for TOF"};
    Configurable<bool> isApplyParticleMID{"isApplyParticleMID", false, "Apply particle misidentification"};
    Configurable<bool> isApplyParticleMIDPtDep{"isApplyParticleMIDPtDep", false, "Apply pT dependent MID selection"};

    Configurable<bool> checkVzEvSigLoss{"checkVzEvSigLoss", false, "Check Vz event signal loss"};
    Configurable<bool> isApplyDeepAngle{"isApplyDeepAngle", false, "Deep Angle cut"};
    Configurable<bool> isApplyMCchecksClosure{"isApplyMCchecksClosure", true, "Apply MC checks for closure test"};
    Configurable<float> deltaRCut{"deltaRCut", 0.0f, "Apply deltaR cut between two daughters"};

    // Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
    Configurable<float> configOccCut{"configOccCut", 1000., "Occupancy cut"};

    // Configurables for track selections
    Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"}; // PV Contriuibutor
    Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", false, "Primary track selection"};          // kGoldenChi2 | kDCAxy | kDCAz
    Configurable<bool> isGlobalTracks{"isGlobalTracks", true, "isGlobalTracks"};

    Configurable<int> rotationalCut{"rotationalCut", 10, "Cut value (Rotation angle pi - pi/cut and pi + pi/cut)"};
    // Configurable<float> cfgCutPT{"cfgCutPT", 0.2f, "PT cut on daughter track"};
    // Configurable<float> cfgCutEtaMax{"cfgCutEtaMax", 0.8f, "Eta cut on daughter track"};
    // Configurable<float> cfgCutEtaMin{"cfgCutEtaMin", 0.0f, "Eta cut on daughter track"};
    // Configurable<float> cfgCutDCAxyMin{"cfgCutDCAxyMin", 0.0f, "DCAxy range for tracks"};
    Configurable<float> cfgCutDCAxyMax{"cfgCutDCAxyMax", 2.0f, "DCAxy range for tracks"};
    Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
    Configurable<float> ctrackRapidity{"ctrackRapidity", 0.3f, "Cut on track rapidity"};
    Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};
    Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
    Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
    Configurable<float> cfgRCRFC{"cfgRCRFC", 0.8f, "Crossed Rows to Findable Clusters"};
    Configurable<float> cfgITSChi2NCl{"cfgITSChi2NCl", 36.0, "ITS Chi2/NCl"};
    Configurable<float> cfgTPCChi2NClMax{"cfgTPCChi2NClMax", 4.0, "TPC Chi2/NCl"};
    Configurable<float> cfgTPCChi2NClMin{"cfgTPCChi2NClMin", 0.0, "TPC Chi2/NCl"};
    Configurable<bool> hasITS{"hasITS", true, "Required ITS"};
    Configurable<bool> isITSTPCRefit{"isITSTPCRefit", false, "Require ITS Refit"};
    Configurable<bool> isVertexITSTPC{"isVertexITSTPC", false, "Vertex ITS TPC"};
    Configurable<bool> isVertexTOFMatched{"isVertexTOFMatched", false, "Vertex TOF Matched"};
    Configurable<bool> isNoCollInTimeRangeStandard{"isNoCollInTimeRangeStandard", false, "No collision in time range standard"};
    Configurable<bool> isApplyPtDepDCAxyCut{"isApplyPtDepDCAxyCut", false, "Apply pT dependent DCAxy cut"};
    // Configurable<bool> isGoldenChi2{"isGoldenChi2", false, "Apply golden chi2 cut"};
    // Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};

    // cuts on mother
    // Configurable<bool> isApplyCutsOnMother{"isApplyCutsOnMother", false, "Enable additional cuts on Kstar mother"};
    // Configurable<float> cMaxPtMotherCut{"cMaxPtMotherCut", 15.0, "Maximum pt of mother cut"};
    // Configurable<float> cMaxMinvMotherCut{"cMaxMinvMotherCut", 1.5, "Maximum mass of mother cut"};
    // Configurable<float> rapidityMotherData{"rapidityMotherData", 0.5, "Maximum rapidity of mother"};
    Configurable<bool> isPDGCheckMC{"isPDGCheckMC", true, "Check PDG code in MC (false for MC closure test)"};

    // PID selections
    Configurable<float> nsigmaCutTPCPi{"nsigmaCutTPCPi", 3.0, "TPC Nsigma cut for pions"};
    Configurable<float> nsigmaCutTPCKa{"nsigmaCutTPCKa", 3.0, "TPC Nsigma cut for kaons"};
    Configurable<float> nsigmaCutTOFPi{"nsigmaCutTOFPi", 3.0, "TOF Nsigma cut for pions"};
    Configurable<float> nsigmaCutTOFKa{"nsigmaCutTOFKa", 3.0, "TOF Nsigma cut for kaons"};
    Configurable<float> nsigmaCutCombinedKa{"nsigmaCutCombinedKa", 3.0, "Combined Nsigma cut for kaon"};
    Configurable<float> nsigmaCutCombinedPi{"nsigmaCutCombinedPi", 3.0, "Combined Nsigma cut for pion"};
    // Configurable<float> nsigmaCutCombinedMIDKa{"nsigmaCutCombinedMIDKa", 3.0, "Combined Nsigma cut for kaon in MID"};
    Configurable<float> nsigmaCutCombinedMID{"nsigmaCutCombinedMID", 3.0, "Combined Nsigma cut for pion in MID"};
    Configurable<float> nsigmaCutTPCMID{"nsigmaCutTPCMID", 1.0, "MID Nsigma cut for pion in TPC"};
    // Configurable<float> nsigmaCutTPCMIDKa{"nsigmaCutTPCMIDKa", 1.0, "MID Nsigma cut for kaon in TPC"};
    // Configurable<float> nsigmaCutTOFMIDPi{"nsigmaCutTOFMIDPi", 1.0, "MID Nsigma cut for pion in TOF"};
    // Configurable<float> nsigmaCutTOFMIDKa{"nsigmaCutTOFMIDKa", 1.0, "MID Nsigma cut for kaon in TOF"};
    // Configurable<float> nsigmaCutTPCMIDPr{"nsigmaCutTPCMIDPr", 3.0, "TPC Nsigma cut for protons (for MID)"};
    // Configurable<float> nsigmaCutTOFMIDPr{"nsigmaCutTOFMIDPr", 3.0, "TOF Nsigma cut for protons (for MID)"};

    // Other fixed variables
    float lowPtCutPID = 0.5;
    int noOfDaughters = 2;
    // float rapidityMotherData = 0.5;
    float cutzvertex = 10.0f;
    float cfgCutEtaMax = 0.8f;
    float cfgCutPT = 0.2f;
    float cfgDeepAngle = 0.04;
    float rapidityMotherData = 0.5;
    float pidCutTemp = 2.0f;

  } selectionConfig;

  enum MultEstimator {
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

  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hInvMass{"hInvMass", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hPID{"hPID", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry hOthers{"hOthers", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Confugrable for QA histograms
  Configurable<bool> calcLikeSign{"calcLikeSign", true, "Calculate Like Sign"};
  Configurable<bool> calcRotational{"calcRotational", false, "Calculate Rotational"};
  Configurable<bool> cQAplots{"cQAplots", true, "cQAplots"};
  Configurable<bool> cPIDAdditionalChecks{"cPIDAdditionalChecks", true, "cPIDAdditionalChecks"};
  Configurable<bool> cQAevents{"cQAevents", true, "Multiplicity dist, DCAxy, DCAz"};
  Configurable<bool> onlyTOF{"onlyTOF", false, "only TOF tracks"};
  // Configurable<bool> onlyTOFHIT{"onlyTOFHIT", false, "accept only TOF hit tracks at high pt"};
  Configurable<bool> onlyTPC{"onlyTPC", true, "only TPC tracks"};
  Configurable<int> cRotations{"cRotations", 3, "Number of random rotations in the rotational background"};
  Configurable<int> cSelectMultEstimator{"cSelectMultEstimator", 0, "Select multiplicity estimator: 0 - FT0M, 1 - FT0A, 2 - FT0C, 3 - FV0A"};
  Configurable<bool> applyRecMotherRapidity{"applyRecMotherRapidity", true, "Apply rapidity cut on reconstructed mother track"};

  // Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", false, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<float> cBetaCutTOF{"cBetaCutTOF", 0.0, "cut TOF beta"};
  // Configurable<bool> cFakeTrack{"cFakeTrack", true, "Fake track selection"};
  // Configurable<float> cFakeTrackCutKa{"cFakeTrackCutKa", 0.5, "Cut based on momentum difference in global and TPC tracks for Kaons"};
  // Configurable<float> cFakeTrackCutPi{"cFakeTrackCutPi", 0.5, "Cut based on momentum difference in global and TPC tracks for Pions"};

  // Configurable for histograms
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", true, "avoid split track in MC"};
  // Configurable<bool> cAllGenCollisions{"cAllGenCollisions", false, "To fill all generated collisions for the signal loss calculations"};
  ConfigurableAxis binsMultPlot{"binsMultPlot", {110, 0.0, 110}, "THnSpare multiplicity axis"};
  // ConfigurableAxis axisdEdx{"axisdEdx", {1, 0.0f, 200.0f}, "dE/dx (a.u.)"};
  // ConfigurableAxis axisPtfordEbydx{"axisPtfordEbydx", {1, 0, 20}, "pT (GeV/c)"};
  // ConfigurableAxis axisMultdist{"axisMultdist", {1, 0, 70000}, "Multiplicity distribution"};
  ConfigurableAxis configThnAxisPOL{"configThnAxisPOL", {20, -1.0, 1.0}, "Costheta axis"};
  ConfigurableAxis invMassKstarAxis{"invMassKstarAxis", {300, 0.7f, 1.3f}, "Kstar invariant mass axis"};
  ConfigurableAxis ptAxisKstar{"ptAxisKstar", {200, 0.0f, 20.0f}, "Kstar pT axis"};
  ConfigurableAxis binsImpactPar{"binsImpactPar", {100, 0, 25}, "Binning of the impact parameter axis"};
  ConfigurableAxis axisNch{"axisNch", {100, 0.0f, 100.0f}, "Number of charged particles in |y| < 0.5"};

  // Event plane configurables
  Configurable<bool> boostDaugter1{"boostDaugter1", false, "Boost daughter Kaon in the COM frame"};
  Configurable<bool> boostDaugter2{"boostDaugter2", true, "Boost daughter Pion in the COM frame"};
  Configurable<bool> activateTHnSparseCosThStarHelicity{"activateTHnSparseCosThStarHelicity", true, "Activate the THnSparse with cosThStar w.r.t. helicity axis"};
  Configurable<bool> activateTHnSparseCosThStarProduction{"activateTHnSparseCosThStarProduction", false, "Activate the THnSparse with cosThStar w.r.t. production axis"};
  Configurable<bool> activateTHnSparseCosThStarBeam{"activateTHnSparseCosThStarBeam", false, "Activate the THnSparse with cosThStar w.r.t. beam axis (Gottified jackson frame)"};
  Configurable<bool> activateTHnSparseCosThStarRandom{"activateTHnSparseCosThStarRandom", false, "Activate the THnSparse with cosThStar w.r.t. random axis"};

  TRandom* rn = new TRandom();

  void init(InitContext const&)
  {
    rctChecker.init(rctCut.cfgEvtRCTFlagCheckerLabel, rctCut.cfgEvtRCTFlagCheckerZDCCheck, rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);
    // Axes
    AxisSpec vertexZAxis = {60, -15., 15., "vrtx_{Z} [cm] for plots"};
    AxisSpec ptAxis = {ptAxisKstar, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invmassAxis = {invMassKstarAxis, "Invariant mass (GeV/#it{c}^{2})"};
    AxisSpec thnAxisPOL{configThnAxisPOL, "cos(#theta)"};
    AxisSpec multiplicityAxis = {binsMultPlot, "Multiplicity Axis"};
    AxisSpec impactParAxis = {binsImpactPar, "Impact Parameter (cm)"};

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    rEventSelection.add("hMultiplicity", "Multiplicity percentile", kTH1F, {{110, 0, 110}});

    rEventSelection.add("hEventCut", "No. of event after cuts", kTH1D, {{20, 0, 20}});
    std::shared_ptr<TH1> hCutFlow = rEventSelection.get<TH1>(HIST("hEventCut"));

    auto check = [](bool enabled) { return enabled ? "" : " #otimes"; }; // check if a cut is enabled and put #otimes beside that label if not enabled

    std::vector<std::string> eveCutLabels = {
      "All Events",
      Form("|Vz| < %.1f", selectionConfig.cutzvertex),
      "sel8",
      std::string("kNoTimeFrameBorder") + check(selectionConfig.isNoTimeFrameBorder.value),
      std::string("kNoITSROFrameBorder") + check(selectionConfig.isNoITSROFrameBorder.value),
      std::string("kNoSameBunchPileup") + check(selectionConfig.isNoSameBunchPileup.value),
      std::string("kIsGoodITSLayersAll") + check(selectionConfig.isGoodITSLayersAll.value),
      std::string("kNoCollInTimeRangeStandard") + check(selectionConfig.isNoCollInTimeRangeStandard.value),
      Form("Occupancy < %.0f%s", selectionConfig.configOccCut.value, check(selectionConfig.isApplyOccCut.value)),
      std::string("rctChecker") + check(rctCut.requireRCTFlagChecker.value),
      std::string("kIsTriggerTVX") + check(selectionConfig.isTriggerTVX.value),
      std::string("kIsGoodZvtxFT0vsPV") + check(selectionConfig.isGoodZvtxFT0vsPV.value),
      std::string("IsINELgt0") + check(selectionConfig.isINELgt0.value),
      std::string("isVertexITSTPC") + check(selectionConfig.isVertexITSTPC.value),
      std::string("isVertexTOFMatched") + check(selectionConfig.isVertexTOFMatched.value)};
    // assign labels
    for (size_t i = 0; i < eveCutLabels.size(); ++i) {
      hCutFlow->GetXaxis()->SetBinLabel(i + 1, eveCutLabels[i].c_str());
    }

    // for primary tracksbinsMultPlot
    if (cQAplots) {
      // hOthers.add("dE_by_dx_TPC", "dE/dx signal in the TPC as a function of pT", kTH2F, {axisPtfordEbydx, axisdEdx});
      hOthers.add("hphi", "Phi distribution", kTH1F, {{65, 0, 6.5}});
      hOthers.add("hEta_after", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
      hOthers.add("hCRFC_after", "CRFC after distribution", kTH1F, {{100, 0.0f, 10.0f}});
      hOthers.add("hCRFC_before", "CRFC before distribution", kTH1F, {{100, 0.0f, 10.0f}});

      hOthers.add("hKstar_rap_pt", "Pair rapidity distribution; y; p_{T}; Counts", kTH2F, {{400, -2.0f, 2.0f}, ptAxis});
      hOthers.add("hKstar_eta_pt", "Pair eta distribution; #eta; p_{T}; Counts", kTH2F, {{400, -2.0f, 2.0f}, ptAxis});

      hPID.add("Before/hNsigma_TPC_TOF_Ka_before", "N #sigma Kaon TOF before", kTH3F, {{50, -5.0f, 5.0f}, {50, -5.0f, 5.0f}, ptAxis});
      hPID.add("Before/hNsigma_TPC_TOF_Pi_before", "N #sigma Pion TOF before", kTH3F, {{50, -5.0f, 5.0f}, {50, -5.0f, 5.0f}, ptAxis});

      hPID.add("Before/hTPCnsigKa_mult_pt", "TPC nsigma of Kaon brfore PID with pt and multiplicity", kTH3F, {{100, -10.0f, 10.0f}, multiplicityAxis, ptAxis});
      hPID.add("Before/hTPCnsigPi_mult_pt", "TPC nsigma of Pion brfore PID with pt and multiplicity", kTH3F, {{100, -10.0f, 10.0f}, multiplicityAxis, ptAxis});
      hPID.add("Before/hTOFnsigKa_mult_pt", "TOF nsigma of Kaon brfore PID with pt and multiplicity", kTH3F, {{100, -10.0f, 10.0f}, multiplicityAxis, ptAxis});
      hPID.add("Before/hTOFnsigPi_mult_pt", "TOF nsigma of Pion brfore PID with pt and multiplicity", kTH3F, {{100, -10.0f, 10.0f}, multiplicityAxis, ptAxis});

      hPID.add("After/hNsigma_TPC_TOF_Ka_after", "N #sigma Kaon TOF after", kTH3F, {{50, -5.0f, 5.0f}, {50, -5.0f, 5.0f}, ptAxis});
      hPID.add("After/hNsigma_TPC_TOF_Pi_after", "N #sigma Pion TOF after", kTH3F, {{50, -5.0f, 5.0f}, {50, -5.0f, 5.0f}, ptAxis});

      hPID.add("After/hDcaxyPi", "Dcaxy distribution of selected Pions", kTH1F, {{200, -1.0f, 1.0f}});
      hPID.add("After/hDcaxyKa", "Dcaxy distribution of selected Kaons", kTH1F, {{200, -1.0f, 1.0f}});
      hPID.add("After/hDcazPi", "Dcaz distribution of selected Pions", kTH1F, {{200, -1.0f, 1.0f}});
      hPID.add("After/hDcazKa", "Dcaz distribution of selected Kaons", kTH1F, {{200, -1.0f, 1.0f}});

      hPID.add("After/hTPCnsigKa_mult_pt", "TPC nsigma of Kaon after PID with pt and multiplicity", kTH3F, {{100, -10.0f, 10.0f}, multiplicityAxis, ptAxis});
      hPID.add("After/hTPCnsigPi_mult_pt", "TPC nsigma of Pion after PID with pt and multiplicity", kTH3F, {{100, -10.0f, 10.0f}, multiplicityAxis, ptAxis});
      hPID.add("After/hTOFnsigKa_mult_pt", "TOF nsigma of Kaon after PID with pt and multiplicity", kTH3F, {{100, -10.0f, 10.0f}, multiplicityAxis, ptAxis});
      hPID.add("After/hTOFnsigPi_mult_pt", "TOF nsigma of Pion after PID with pt and multiplicity", kTH3F, {{100, -10.0f, 10.0f}, multiplicityAxis, ptAxis});
    }

    if (cPIDAdditionalChecks) {
      hPID.add("AdditionalChecks/hTPCnsigKa_pt", "TPC nsigma of Kaon with pt", kTH2F, {{100, -10.0f, 10.0f}, ptAxis});
      hPID.add("AdditionalChecks/hTPCnsigPi_pt", "TPC nsigma of Pion with pt", kTH2F, {{100, -10.0f, 10.0f}, ptAxis});

      hPID.add("AdditionalChecks/hTOFnsigKa_ptPrCut", "TOF nsigma of Kaon with pt", kTH2F, {{100, -10.0f, 10.0f}, ptAxis});
      hPID.add("AdditionalChecks/hTOFnsigPi_ptPrCut", "TOF nsigma of Pion with pt", kTH2F, {{100, -10.0f, 10.0f}, ptAxis});
    }

    // KStar histograms
    hInvMass.add("h3KstarInvMassUnlikeSign", "kstar Unlike Sign", kTHnSparseF, {multiplicityAxis, ptAxis, invmassAxis, thnAxisPOL});
    hInvMass.add("h3KstarInvMassMixed", "kstar Mixed", kTHnSparseF, {multiplicityAxis, ptAxis, invmassAxis, thnAxisPOL});
    if (calcLikeSign) {
      hInvMass.add("h3KstarInvMasslikeSignPP", "kstar like Sign", kTHnSparseF, {multiplicityAxis, ptAxis, invmassAxis, thnAxisPOL});
      hInvMass.add("h3KstarInvMasslikeSignMM", "kstar like Sign", kTHnSparseF, {multiplicityAxis, ptAxis, invmassAxis, thnAxisPOL});
    }
    if (calcRotational)
      hInvMass.add("h3KstarInvMassRotated", "kstar rotated", kTHnSparseF, {multiplicityAxis, ptAxis, invmassAxis, thnAxisPOL});

    // MC histograms
    hInvMass.add("hk892GenpT", "pT distribution of True MC K(892)0", kTHnSparseF, {ptAxis, multiplicityAxis});
    hInvMass.add("hk892GenpT2", "pT distribution of True MC K(892)0", kTHnSparseF, {ptAxis, multiplicityAxis});
    hInvMass.add("hk892GenpTCalib1", "pT distribution of True MC K(892)0", kTHnSparseF, {ptAxis, multiplicityAxis});
    hInvMass.add("hk892GenpTCalib2", "pT distribution of True MC K(892)0", kTHnSparseF, {ptAxis, multiplicityAxis});
    hInvMass.add("h1KstarRecMass", "Invariant mass of kstar meson", kTH1F, {invmassAxis});
    hInvMass.add("h2KstarRecpt1", "pT of kstar meson", kTHnSparseF, {ptAxis, multiplicityAxis, invmassAxis});
    hInvMass.add("h2KstarRecpt2", "pT of kstar meson", kTHnSparseF, {ptAxis, multiplicityAxis, invmassAxis});
    hInvMass.add("h2KstarRecptCalib1", "pT of kstar meson", kTHnSparseF, {ptAxis, multiplicityAxis, invmassAxis});
    hInvMass.add("h2KstarRecptCalib2", "pT of kstar meson", kTHnSparseF, {ptAxis, multiplicityAxis, invmassAxis});
    hInvMass.add("h1genmass", "Invariant mass of generated kstar meson", kTH1F, {invmassAxis});
    hInvMass.add("h1GenMult", "Multiplicity generated", kTH1F, {multiplicityAxis});
    hInvMass.add("h1GenMult2", "Multiplicity generated (direct)", kTH1F, {multiplicityAxis});
    hInvMass.add("h1RecMult", "Multiplicity reconstructed", kTH1F, {multiplicityAxis});
    hInvMass.add("h1RecMult2", "Multiplicity reconstructed", kTH1F, {multiplicityAxis});
    hInvMass.add("h1KSRecsplit", "KS meson Rec split", kTH1F, {{100, 0.0f, 10.0f}});
    // hInvMass.add("hAllGenCollisionsImpact", "All generated collisions vs impact parameter", kTH1F, {multiplicityAxis});
    hInvMass.add("hAllGenCollisions", "All generated events", kTH1F, {multiplicityAxis});
    hInvMass.add("hAllGenCollisions1Rec", "All gen events with at least one rec event", kTH1F, {multiplicityAxis});
    hInvMass.add("hAllKstarGenCollisisons", "All generated Kstar in events with rapidity in 0.5", kTH2F, {{multiplicityAxis}, {ptAxis}});
    hInvMass.add("hAllKstarGenCollisisons1Rec", "All generated Kstar in events with at least one rec event in rapidity in 0.5", kTH2F, {{multiplicityAxis}, {ptAxis}});
    hInvMass.add("hAllRecCollisions", "All reconstructed events", kTH1F, {multiplicityAxis});
    hInvMass.add("hAllRecCollisionsCalib", "All reconstructed events", kTH1F, {multiplicityAxis});
    hInvMass.add("sigEvLossFromGenRec/MultiplicityGen", "Multiplicity in generated MC", kTH1F, {multiplicityAxis});
    hInvMass.add("sigEvLossFromGenRec/MultiplicityRec", "Multiplicity in generated MC with at least 1 reconstruction", kTH1F, {multiplicityAxis});
    hInvMass.add("sigEvLossFromGenRec/hSignalLossDenominator", "Kstar generated before event selection", kTH2F, {{ptAxis}, {multiplicityAxis}});
    hInvMass.add("sigEvLossFromGenRec/hSignalLossNumerator", "Kstar generated after event selection", kTH2F, {{ptAxis}, {multiplicityAxis}});

    if (doprocessEvtLossSigLossMC || doprocessEvtLossSigLossMCPhi) {
      hInvMass.add("MCcorrections/hSignalLossDenominator", "Kstar generated before event selection", kTH2F, {{ptAxis}, {multiplicityAxis}});
      hInvMass.add("MCcorrections/hSignalLossNumerator", "Kstar generated after event selection", kTH2F, {{ptAxis}, {multiplicityAxis}});
      hInvMass.add("MCcorrections/hImpactParameterRec", "Impact parameter in reconstructed MC", kTH1F, {impactParAxis});
      hInvMass.add("MCcorrections/MultiplicityRec", "Multiplicity in generated MC with at least 1 reconstruction", kTH1F, {multiplicityAxis});
      hInvMass.add("MCcorrections/hImpactParameterGen", "Impact parameter in generated MC", kTH1F, {impactParAxis});
      hInvMass.add("MCcorrections/MultiplicityGen", "Multiplicity in generated MC", kTH1F, {multiplicityAxis});
      hInvMass.add("MCcorrections/hImpactParametervsMultiplicity", "Impact parameter vs multiplicity in reconstructed MC", kTH2F, {{impactParAxis}, {multiplicityAxis}});
    }

    // Signal Loss & Event Loss in Light Ion Collisions
    if (doprocessEvtLossSigLossLightIonMC) {
      hInvMass.add("MCcorrections/hImpactParameterGen", "Impact parameter of generated MC events", kTH1F, {impactParAxis});
      hInvMass.add("MCcorrections/hImpactParameterRec", "Impact parameter of selected MC events", kTH1F, {impactParAxis});
      hInvMass.add("MCcorrections/hImpactParvsCentrRec", "Impact parameter of selected MC events vs centrality", kTH2F, {{multiplicityAxis}, impactParAxis});
      hInvMass.add("MCcorrections/hKstarGenBeforeEvtSel", "K*0 before event selections", kTH2F, {ptAxis, impactParAxis});
      hInvMass.add("MCcorrections/hKstarGenAfterEvtSel", "K*0 after event selections", kTH2F, {ptAxis, impactParAxis});
    }

    if (doprocessEvSigLossFactors || doprocessEvSigLossFactorsPhi) {
      hInvMass.add("CorrFactors/hMultiplicityVsMultMC", "Event centrality vs MC centrality", kTH2F, {{101, 0.0f, 101.0f}, axisNch});
      hInvMass.add("CorrFactors/hEventCentrality", "Event centrality", kTH1F, {{101, 0, 101}});
      hInvMass.add("CorrFactors/hNrecInGen", "Number of collisions in MC", kTH1F, {{4, -0.5, 3.5}});
      hInvMass.add("CorrFactors/hGenEvents", "Generated events", kTH2F, {{axisNch}, {4, 0, 4}});
      auto hGenEvents = hInvMass.get<TH2>(HIST("CorrFactors/hGenEvents"));
      hGenEvents->GetYaxis()->SetBinLabel(1, "All generated events");
      hGenEvents->GetYaxis()->SetBinLabel(2, "Generated events with Mc collision V_{z} cut");
      hGenEvents->GetYaxis()->SetBinLabel(3, "Generated events with Mc INEL>0");
      hGenEvents->GetYaxis()->SetBinLabel(4, "Generated events with at least one reconstructed event");

      hInvMass.add("CorrFactors/h2dGenKstar", "Centrality vs p_{T}", kTH2D, {{101, 0.0f, 101.0f}, ptAxis});
      hInvMass.add("CorrFactors/h3dGenKstarVsMultMCVsMultiplicity", "MC centrality vs centrality vs p_{T}", kTH3D, {axisNch, {101, 0.0f, 101.0f}, ptAxis});
      hInvMass.add("CorrFactors/hSignalLoss1", "Kstar generated before event selection", kTH2F, {{ptAxis}, {axisNch}});
      hInvMass.add("CorrFactors/hSignalLoss2", "Kstar generated after event selection", kTH2F, {{ptAxis}, {axisNch}});
      hInvMass.add("CorrFactors/MultiplicityRec", "Multiplicity in generated MC with at least 1 reconstruction", kTH1F, {axisNch});
      hInvMass.add("CorrFactors/MultiplicityGen", "Multiplicity in generated MC", kTH1F, {axisNch});
    }

    rEventSelection.add("tracksCheckData", "No. of events in the data", kTH1I, {{10, 0, 10}});
    rEventSelection.add("eventsCheckGen", "No. of events in the generated MC", kTH1I, {{10, 0, 10}});
    rEventSelection.add("recMCparticles", "No. of events in the reconstructed MC", kTH1I, {{20, 0, 20}});
    rEventSelection.add("hOccupancy", "Occupancy distribution", kTH1F, {{1000, 0, 15000}});

    std::shared_ptr<TH1> hrecLabel = rEventSelection.get<TH1>(HIST("recMCparticles"));
    hrecLabel->GetXaxis()->SetBinLabel(1, "All tracks");
    hrecLabel->GetXaxis()->SetBinLabel(2, "has_MC");
    hrecLabel->GetXaxis()->SetBinLabel(3, "Track selection");
    hrecLabel->GetXaxis()->SetBinLabel(4, "StrictlyUpperIndex");
    hrecLabel->GetXaxis()->SetBinLabel(5, "Unlike Sign");
    hrecLabel->GetXaxis()->SetBinLabel(6, "Physical Primary");
    hrecLabel->GetXaxis()->SetBinLabel(7, "TrackPDG Check1");
    hrecLabel->GetXaxis()->SetBinLabel(8, "TrackPDG Check2");
    hrecLabel->GetXaxis()->SetBinLabel(9, "Global Index");
    hrecLabel->GetXaxis()->SetBinLabel(10, "Generator");
    hrecLabel->GetXaxis()->SetBinLabel(11, "Mother y");
    hrecLabel->GetXaxis()->SetBinLabel(12, "Mother PDG");
    hrecLabel->GetXaxis()->SetBinLabel(13, "Track PID");
    hrecLabel->GetXaxis()->SetBinLabel(14, "Track MID");
    hrecLabel->GetXaxis()->SetBinLabel(15, "Track y");
    hrecLabel->GetXaxis()->SetBinLabel(16, "Split tracks");
    hrecLabel->GetXaxis()->SetBinLabel(17, "DeepAngle Cut");

    std::shared_ptr<TH1> hDataTracks = rEventSelection.get<TH1>(HIST("tracksCheckData"));
    hDataTracks->GetXaxis()->SetBinLabel(1, "All tracks");
    hDataTracks->GetXaxis()->SetBinLabel(2, "Track selection");
    hDataTracks->GetXaxis()->SetBinLabel(3, "PID Cut");
    hDataTracks->GetXaxis()->SetBinLabel(4, "Remove Fake Tracks");
    hDataTracks->GetXaxis()->SetBinLabel(5, "Rapidity Cut");
    hDataTracks->GetXaxis()->SetBinLabel(6, "MID");
    hDataTracks->GetXaxis()->SetBinLabel(7, "DeepAngle Cut");

    std::shared_ptr<TH1> hGenTracks = rEventSelection.get<TH1>(HIST("eventsCheckGen"));
    hGenTracks->GetXaxis()->SetBinLabel(1, "All events");
    hGenTracks->GetXaxis()->SetBinLabel(2, "INELgt0+vtz");
    hGenTracks->GetXaxis()->SetBinLabel(3, "INELgt0");
    hGenTracks->GetXaxis()->SetBinLabel(4, "All Collisions");
    hGenTracks->GetXaxis()->SetBinLabel(5, "Event Reconstructed");

    // Multplicity distribution
    if (cQAevents) {
      // rEventSelection.add("multdist_FT0M", "FT0M Multiplicity distribution", kTH1F, {axisMultdist});
      // hInvMass.add("multdist_FT0A", "FT0A Multiplicity distribution", kTH1F, {axisMultdist});
      // hInvMass.add("multdist_FT0C", "FT0C Multiplicity distribution", kTH1F, {axisMultdist});
      // hInvMass.add("hNcontributor", "Number of primary vertex contributor", kTH1F, {{2000, 0.0f, 10000.0f}});
      rEventSelection.add("hDcaxy_cent_pt", "Dcaxy distribution", kTH3F, {{200, -1.0f, 1.0f}, multiplicityAxis, ptAxis});
      rEventSelection.add("hDcaz_cent_pt", "Dcaz distribution", kTH3F, {{200, -1.0f, 1.0f}, multiplicityAxis, ptAxis});
    }
  }

  // double massPi = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass(); // FIXME: Get from the common header
  double massPi = o2::constants::physics::MassPiPlus;
  double massKa = o2::constants::physics::MassKPlus;

  template <typename Coll>
  bool selectionEvent(const Coll& collision, bool fillHist = false) // default to false
  {
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 0);

    if (std::abs(collision.posZ()) > selectionConfig.cutzvertex)
      return false;
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 1);

    if (selectionConfig.isSel8 && !collision.sel8())
      return false;
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 2);

    if (selectionConfig.isNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return false;
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 3);

    if (selectionConfig.isNoITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))
      return false;
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 4);

    if (selectionConfig.isNoSameBunchPileup && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)))
      return false;
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 5);

    if (selectionConfig.isGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll))
      return false;
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 6);

    if (selectionConfig.isNoCollInTimeRangeStandard && (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)))
      return false;
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 7);

    if (selectionConfig.isApplyOccCut && (std::abs(collision.trackOccupancyInTimeRange()) > selectionConfig.configOccCut))
      return false;
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 8);

    if (rctCut.requireRCTFlagChecker && !rctChecker(collision))
      return false;
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 9);

    if (selectionConfig.isTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX))
      return false;
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 10);

    if (selectionConfig.isGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
      return false;
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 11);

    if (selectionConfig.isINELgt0 && !collision.isInelGt0()) {
      return false;
    }
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 12);

    if (selectionConfig.isVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 13);

    if (selectionConfig.isVertexTOFMatched && !collision.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if (fillHist)
      rEventSelection.fill(HIST("hEventCut"), 14);

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
      // if (std::abs(candidate.eta()) > selectionConfig.cfgCutEtaMax || std::abs(candidate.eta()) < selectionConfig.cfgCutEtaMin)
      if (std::abs(candidate.eta()) > selectionConfig.cfgCutEtaMax)
        return false;
      if (!selectionConfig.isApplyPtDepDCAxyCut) {
        // if (std::abs(candidate.dcaXY()) > selectionConfig.cfgCutDCAxyMax || std::abs(candidate.dcaXY()) < selectionConfig.cfgCutDCAxyMin)
        if (std::abs(candidate.dcaXY()) > selectionConfig.cfgCutDCAxyMax)
          return false;
      } else {
        if (std::abs(candidate.dcaXY()) > (0.0105 + 0.035 / std::pow(candidate.pt(), 1.1)))
          return false;
      }
      // if (selectionConfig.isGoldenChi2 && !candidate.passedGoldenChi2())
      //   return false;
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
      if (selectionConfig.cfgPVContributor && !candidate.isPVContributor())
        return false;
      if (selectionConfig.hasITS && !candidate.hasITS())
        return false;
      if (selectionConfig.isITSTPCRefit && (!(o2::aod::track::ITSrefit) || !(o2::aod::track::TPCrefit)))
        return false;
    } else if (!selectionConfig.isGlobalTracks) {
      if (std::abs(candidate.pt()) < selectionConfig.cfgCutPT)
        return false;
      // if (std::abs(candidate.eta()) > selectionConfig.cfgCutEtaMax || std::abs(candidate.eta()) < selectionConfig.cfgCutEtaMin)
      if (std::abs(candidate.eta()) > selectionConfig.cfgCutEtaMax)
        return false;
      // if (std::abs(candidate.dcaXY()) > selectionConfig.cfgCutDCAxyMax || std::abs(candidate.dcaXY()) < selectionConfig.cfgCutDCAxyMin)
      if (std::abs(candidate.dcaXY()) > selectionConfig.cfgCutDCAxyMax)
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
      if (selectionConfig.cfgPVContributor && !candidate.isPVContributor())
        return false;
      if (selectionConfig.cfgPrimaryTrack && !candidate.isPrimaryTrack())
        return false;
    }

    return true;
  }

  // template <typename T>
  // bool isFakeTrack(const T& track, int PID)
  // {
  //   const auto pglobal = track.p();
  //   const auto ptpc = track.tpcInnerParam();
  //   if (PID == 0 && std::abs(pglobal - ptpc) > cFakeTrackCutPi) {
  //     return true;
  //   }
  //   if (PID == 1 && std::abs(pglobal - ptpc) > cFakeTrackCutKa) {
  //     return true;
  //   }
  //   return false;
  // }

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
    double deltaRvalue = std::sqrt(TVector2::Phi_mpi_pi(candidate1.phi() - candidate2.phi()) * TVector2::Phi_mpi_pi(candidate1.phi() - candidate2.phi()) + (candidate1.eta() - candidate2.eta()) * (candidate1.eta() - candidate2.eta()));
    if (deltaRvalue < selectionConfig.deltaRCut) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool selectionPID(const T& candidate, int PID)
  {
    if (PID == PIDParticle::kPion) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < selectionConfig.nsigmaCutTOFPi && candidate.beta() > cBetaCutTOF) {
          return true;
        }
      } else if (onlyTPC) {
        if (std::abs(candidate.tpcNSigmaPi()) < selectionConfig.nsigmaCutTPCPi) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < (selectionConfig.nsigmaCutCombinedPi * selectionConfig.nsigmaCutCombinedPi) && candidate.beta() > cBetaCutTOF) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < selectionConfig.nsigmaCutTPCPi) {
          return true;
        }
      }
    } else if (PID == PIDParticle::kKaon) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < selectionConfig.nsigmaCutTOFKa && candidate.beta() > cBetaCutTOF) {
          return true;
        }
      } else if (onlyTPC) {
        if (std::abs(candidate.tpcNSigmaKa()) < selectionConfig.nsigmaCutTPCKa) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (selectionConfig.nsigmaCutCombinedKa * selectionConfig.nsigmaCutCombinedKa) && candidate.beta() > cBetaCutTOF) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < selectionConfig.nsigmaCutTPCKa) {
          return true;
        }
      }
    }
    /* else if (PID == PIDParticle::kProton) { // for proton
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPr()) < selectionConfig.nsigmaCutTOFPr && candidate.beta() > cBetaCutTOF) {
          return true;
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPr()) < selectionConfig.nsigmaCutTOFPr && candidate.beta() > cBetaCutTOF) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < selectionConfig.nsigmaCutTPCPr) {
          return true;
        }
      } else if (onlyTPC) {
        if (std::abs(candidate.tpcNSigmaPr()) < selectionConfig.nsigmaCutTPCPr) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaPr() * candidate.tofNSigmaPr() + candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr()) < (selectionConfig.nsigmaCutTOFPr * selectionConfig.nsigmaCutTOFPr) && candidate.beta() > cBetaCutTOF) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < selectionConfig.nsigmaCutTPCPr) {
          return true;
        }
      }
    } */
    return false;
  }

  template <typename T>
  bool selectionMID(const T& candidate, int PID)
  {
    if (PID == PIDParticle::kPion) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < selectionConfig.nsigmaCutTPCMID) {
          return true;
        }
      } else if (onlyTPC) {
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
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < selectionConfig.nsigmaCutTPCMID) {
          return true;
        }
      } else if (onlyTPC) {
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
    // else if (PID == PIDParticle::kProton) { // for proton
    //   if (onlyTOF) {
    //     if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPr()) < selectionConfig.nsigmaCutTOFPr && candidate.beta() > cBetaCutTOF) {
    //       return true;
    //     }
    //   } else if (onlyTOFHIT) {
    //     if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPr()) < selectionConfig.nsigmaCutTOFPr && candidate.beta() > cBetaCutTOF) {
    //       return true;
    //     }
    //     if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < selectionConfig.nsigmaCutTPCPr) {
    //       return true;
    //     }
    //   } else if (onlyTPC) {
    //     if (std::abs(candidate.tpcNSigmaPr()) < selectionConfig.nsigmaCutTPCPr) {
    //       return true;
    //     }
    //   } else {
    //     if (candidate.hasTOF() && (candidate.tofNSigmaPr() * candidate.tofNSigmaPr() + candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr()) < (selectionConfig.nsigmaCutTOFPr * selectionConfig.nsigmaCutTOFPr) && candidate.beta() > cBetaCutTOF) {
    //       return true;
    //     }
    //     if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < selectionConfig.nsigmaCutTPCPr) {
    //       return true;
    //     }
    //   }
    // }
    return false;
  }

  template <typename T>
  bool selectionMIDPtDep(const T& candidate, int PID)
  {
    if (PID == PIDParticle::kPion) {
      if (candidate.pt() >= 1.0 && candidate.pt() < 2.5 && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < selectionConfig.nsigmaCutTPCMID) {
        return true;
      }
    } else if (PID == PIDParticle::kKaon) {
      if (candidate.pt() >= 0.7 && candidate.pt() < 2.5 && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < selectionConfig.nsigmaCutTPCMID) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selectionPIDPtDep(const T& candidate, int PID)
  {
    if (PID == PIDParticle::kPion) {
      if (candidate.pt() < selectionConfig.lowPtCutPID && std::abs(candidate.tpcNSigmaPi()) < selectionConfig.nsigmaCutTPCPi) {
        return true;
      }
      if (candidate.pt() >= selectionConfig.lowPtCutPID && !candidate.hasTOF()) {
        if (candidate.pt() < 1.0 && std::abs(candidate.tpcNSigmaPi()) < selectionConfig.nsigmaCutTPCPi) {
          return true;
        }
        if (candidate.pt() >= 1.0 && candidate.pt() < 2.0 && candidate.tpcNSigmaPi() > -1.5 && candidate.tpcNSigmaPi() < selectionConfig.nsigmaCutTPCPi) {
          return true;
        }
        if (candidate.pt() >= 2.0 && candidate.pt() < 5.0 && candidate.tpcNSigmaPi() > -1.0 && candidate.tpcNSigmaPi() < selectionConfig.nsigmaCutTPCPi) {
          return true;
        }
        if (candidate.pt() >= 5.0 && candidate.tpcNSigmaPi() > -0.5 && candidate.tpcNSigmaPi() < selectionConfig.nsigmaCutTPCPi) {
          return true;
        }
      }
      if (candidate.pt() >= selectionConfig.lowPtCutPID && std::abs(candidate.tpcNSigmaPi()) < selectionConfig.nsigmaCutTPCPi && candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < (selectionConfig.nsigmaCutCombinedPi * selectionConfig.nsigmaCutCombinedPi)) {
        return true;
      }
    } else if (PID == PIDParticle::kKaon) {
      if (candidate.pt() < selectionConfig.lowPtCutPID && candidate.tpcNSigmaKa() > -2.0 && candidate.tpcNSigmaKa() < selectionConfig.nsigmaCutTPCKa) {
        return true;
      }
      if (candidate.pt() >= selectionConfig.lowPtCutPID && !candidate.hasTOF()) {
        if (candidate.pt() >= 0.5 && candidate.pt() < 0.6 && candidate.tpcNSigmaKa() > -1.5 && candidate.tpcNSigmaKa() < selectionConfig.nsigmaCutTPCKa) {
          return true;
        }
        if (candidate.pt() >= 0.6 && candidate.pt() < 0.7 && candidate.tpcNSigmaKa() > -1.0 && candidate.tpcNSigmaKa() < selectionConfig.nsigmaCutTPCKa) {
          return true;
        }
        if (candidate.pt() >= 0.7 && candidate.pt() < 0.8 && candidate.tpcNSigmaKa() > -0.5 && candidate.tpcNSigmaKa() < selectionConfig.nsigmaCutTPCKa) {
          return true;
        }
        if (candidate.pt() >= 0.8 && candidate.pt() < 1.0 && candidate.tpcNSigmaKa() > 0.0 && candidate.tpcNSigmaKa() < selectionConfig.nsigmaCutTPCKa) {
          return true;
        }
        if (candidate.pt() >= 1.0 && candidate.pt() < 2.0 && candidate.tpcNSigmaKa() > -selectionConfig.nsigmaCutTPCKa && candidate.tpcNSigmaKa() < 0.0) {
          return true;
        }
        if (candidate.pt() >= 2.0 && candidate.pt() < 4.0 && candidate.tpcNSigmaKa() > -selectionConfig.nsigmaCutTPCKa && candidate.tpcNSigmaKa() < 1.0) {
          return true;
        }
        if (candidate.pt() >= 4.0 && candidate.tpcNSigmaKa() > -selectionConfig.nsigmaCutTPCKa && candidate.tpcNSigmaKa() < 1.5) {
          return true;
        }
      }
      if (candidate.pt() >= selectionConfig.lowPtCutPID && candidate.hasTOF()) {
        if (candidate.pt() < 2.0 && std::abs(candidate.tpcNSigmaKa()) < selectionConfig.nsigmaCutTPCKa && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (selectionConfig.nsigmaCutCombinedKa * selectionConfig.nsigmaCutCombinedKa)) {
          return true;
        }
        if (candidate.pt() >= 2.0 && std::abs(candidate.tpcNSigmaKa()) < selectionConfig.nsigmaCutTPCKa && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (1.0 * 1.0)) {
          return true;
        }
      }
    }
    return false;
  }

  template <typename T>
  bool selectionPIDptDepTOF(const T& candidate, int PID)
  {
    if (PID == PIDParticle::kPion) {
      if (candidate.pt() < 1 || candidate.pt() > 2.5) {
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
      if (candidate.pt() < 0.7 || candidate.pt() > 2.5) {
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

  std::array<float, 3> pvec0;
  std::array<float, 3> pvec1;

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection
  // requirements
  // Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < selectionConfig.cutzvertex);

  // Filter acceptanceFilter = (nabs(aod::track::eta) < selectionConfig.cfgCutEtaMax && nabs(aod::track::pt) > selectionConfig.cfgCutPT) && (nabs(aod::track::eta) > selectionConfig.cfgCutEtaMin);
  Filter acceptanceFilter = (nabs(aod::track::eta) < selectionConfig.cfgCutEtaMax && nabs(aod::track::pt) > selectionConfig.cfgCutPT);
  // Filter fDCAcutFilter = (nabs(aod::track::dcaXY) < selectionConfig.cfgCutDCAxyMax) && (nabs(aod::track::dcaZ) < selectionConfig.cfgCutDCAz) && (nabs(aod::track::dcaXY) > selectionConfig.cfgCutDCAxyMin);
  Filter fDCAcutFilter = (nabs(aod::track::dcaXY) < selectionConfig.cfgCutDCAxyMax) && (nabs(aod::track::dcaZ) < selectionConfig.cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::CentFV0As, aod::PVMults>>; // aod::CentNGlobals, aod::CentNTPVs, aod::CentMFTs
  using EventCandidatesMix = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::CentFV0As>>;            // aod::CentNGlobals, aod::CentNTPVs, aod::CentMFTs
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTOFbeta, aod::TrackSelectionExtension>>;
  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::MultZeqs, aod::FT0Mults, aod::PVMults, aod::CentFV0As>;
  // using EventCandidatesMC = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::FT0Mults, aod::PVMults, aod::CentFV0As>>;

  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::McTrackLabels, aod::pidTOFbeta, aod::TrackSelectionExtension>>;
  // using EventMCGenerated = soa::Join<aod::McCollisions, aod::MultsExtraMC>; // aod::CentNGlobals, aod::CentNTPVs, aod::CentMFTs
  using EventMCGenerated = soa::Join<aod::McCollisions, aod::McCentFT0Ms, aod::MultsExtraMC>;

  using McCollisionMults = soa::Join<aod::McCollisions, aod::MultMCExtras>;
  using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;

  //*********Varibles declaration***************
  float multiplicity{-1.0}, theta2;
  ROOT::Math::PxPyPzMVector daughter1, daughter2, daughterRot, mother, motherRot, daughterSelected, fourVecDauCM, daughterRotCM;
  ROOT::Math::XYZVector randomVec, beamVec, normalVec;
  bool isMix = false;

  template <typename T1, typename T2>
  void fillInvMass(const T1& daughter1, const T1& daughter2, const T1& mother, float multiplicity, bool isMix, const T2& track1, const T2& track2)
  {
    daughterSelected = (boostDaugter1) ? daughter1 : daughter2; // polarization calculations
    ROOT::Math::Boost boost{mother.BoostToCM()};                // boost mother to center of mass frame
    fourVecDauCM = boost(daughterSelected);                     // boost the frame of daughter same as mother

    // if (std::abs(mother.Rapidity()) < selectionConfig.rapidityMotherData) {
    if (activateTHnSparseCosThStarHelicity) {
      auto cosThetaStarHelicity = mother.Vect().Dot(fourVecDauCM.Vect()) / (std::sqrt(fourVecDauCM.Vect().Mag2()) * std::sqrt(mother.Vect().Mag2()));

      if (track1.sign() * track2.sign() < 0) {
        if (!isMix) {
          if (std::abs(mother.Rapidity()) < selectionConfig.rapidityMotherData) {
            hInvMass.fill(HIST("h3KstarInvMassUnlikeSign"), multiplicity, mother.Pt(), mother.M(), cosThetaStarHelicity);
          }

          for (int i = 0; i < cRotations; i++) {
            theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / selectionConfig.rotationalCut, o2::constants::math::PI + o2::constants::math::PI / selectionConfig.rotationalCut);

            daughterRot = ROOT::Math::PxPyPzMVector(daughter1.Px() * std::cos(theta2) - daughter1.Py() * std::sin(theta2), daughter1.Px() * std::sin(theta2) + daughter1.Py() * std::cos(theta2), daughter1.Pz(), daughter1.M());

            motherRot = daughterRot + daughter2;

            ROOT::Math::Boost boost2{motherRot.BoostToCM()};
            daughterRotCM = boost2(daughterRot);

            auto cosThetaStarHelicityRot = motherRot.Vect().Dot(daughterRotCM.Vect()) / (std::sqrt(daughterRotCM.Vect().Mag2()) * std::sqrt(motherRot.Vect().Mag2()));

            if (calcRotational && std::abs(motherRot.Rapidity()) < selectionConfig.rapidityMotherData)
              hInvMass.fill(HIST("h3KstarInvMassRotated"), multiplicity, motherRot.Pt(), motherRot.M(), cosThetaStarHelicityRot);
          }
        } else if (isMix && std::abs(mother.Rapidity()) < selectionConfig.rapidityMotherData) {
          hInvMass.fill(HIST("h3KstarInvMassMixed"), multiplicity, mother.Pt(), mother.M(), cosThetaStarHelicity);
        }
      } else {
        if (!isMix) {
          if (calcLikeSign && std::abs(mother.Rapidity()) < selectionConfig.rapidityMotherData) {
            if (track1.sign() > 0 && track2.sign() > 0) {
              hInvMass.fill(HIST("h3KstarInvMasslikeSignPP"), multiplicity, mother.Pt(), mother.M(), cosThetaStarHelicity);
            } else if (track1.sign() < 0 && track2.sign() < 0) {
              hInvMass.fill(HIST("h3KstarInvMasslikeSignMM"), multiplicity, mother.Pt(), mother.M(), cosThetaStarHelicity);
            }
          }
        }
      }

    } else if (activateTHnSparseCosThStarProduction) {
      normalVec = ROOT::Math::XYZVector(mother.Py(), -mother.Px(), 0.f);
      auto cosThetaStarProduction = normalVec.Dot(fourVecDauCM.Vect()) / (std::sqrt(fourVecDauCM.Vect().Mag2()) * std::sqrt(normalVec.Mag2()));

      if (track1.sign() * track2.sign() < 0) {
        if (!isMix) {
          if (std::abs(mother.Rapidity()) < selectionConfig.rapidityMotherData) {
            hInvMass.fill(HIST("h3KstarInvMassUnlikeSign"), multiplicity, mother.Pt(), mother.M(), cosThetaStarProduction);
          }
          for (int i = 0; i < cRotations; i++) {
            theta2 = rn->Uniform(0, o2::constants::math::PI);
            daughterRot = ROOT::Math::PxPyPzMVector(daughter1.Px() * std::cos(theta2) - daughter1.Py() * std::sin(theta2), daughter1.Px() * std::sin(theta2) + daughter1.Py() * std::cos(theta2), daughter1.Pz(), daughter1.M());

            motherRot = daughterRot + daughter2;
            if (calcRotational && std::abs(motherRot.Rapidity()) < selectionConfig.rapidityMotherData)
              hInvMass.fill(HIST("h3KstarInvMassRotated"), multiplicity, motherRot.Pt(), motherRot.M(), cosThetaStarProduction);
          }
        } else if (isMix && std::abs(mother.Rapidity()) < selectionConfig.rapidityMotherData) {
          hInvMass.fill(HIST("h3KstarInvMassMixed"), multiplicity, mother.Pt(), mother.M(), cosThetaStarProduction);
        }
      } else {
        if (!isMix) {
          if (calcLikeSign && std::abs(mother.Rapidity()) < selectionConfig.rapidityMotherData) {
            if (track1.sign() > 0 && track2.sign() > 0) {
              hInvMass.fill(HIST("h3KstarInvMasslikeSignPP"), multiplicity, mother.Pt(), mother.M(), cosThetaStarProduction);
            } else if (track1.sign() < 0 && track2.sign() < 0) {
              hInvMass.fill(HIST("h3KstarInvMasslikeSignMM"), multiplicity, mother.Pt(), mother.M(), cosThetaStarProduction);
            }
          }
        }
      }
    } else if (activateTHnSparseCosThStarBeam) {
      beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
      auto cosThetaStarBeam = beamVec.Dot(fourVecDauCM.Vect()) / std::sqrt(fourVecDauCM.Vect().Mag2());

      if (track1.sign() * track2.sign() < 0) {
        if (!isMix) {
          if (std::abs(mother.Rapidity()) < selectionConfig.rapidityMotherData) {
            hInvMass.fill(HIST("h3KstarInvMassUnlikeSign"), multiplicity, mother.Pt(), mother.M(), cosThetaStarBeam);
          }
          for (int i = 0; i < cRotations; i++) {
            theta2 = rn->Uniform(0, o2::constants::math::PI);
            daughterRot = ROOT::Math::PxPyPzMVector(daughter1.Px() * std::cos(theta2) - daughter1.Py() * std::sin(theta2), daughter1.Px() * std::sin(theta2) + daughter1.Py() * std::cos(theta2), daughter1.Pz(), daughter1.M());

            motherRot = daughterRot + daughter2;
            if (calcRotational && std::abs(motherRot.Rapidity()) < selectionConfig.rapidityMotherData)
              hInvMass.fill(HIST("h3KstarInvMassRotated"), multiplicity, motherRot.Pt(), motherRot.M(), cosThetaStarBeam);
          }
        } else if (isMix && std::abs(mother.Rapidity()) < selectionConfig.rapidityMotherData) {
          hInvMass.fill(HIST("h3KstarInvMassMixed"), multiplicity, mother.Pt(), mother.M(), cosThetaStarBeam);
        }
      } else {
        if (calcLikeSign && std::abs(mother.Rapidity()) < selectionConfig.rapidityMotherData) {
          if (track1.sign() > 0 && track2.sign() > 0) {
            hInvMass.fill(HIST("h3KstarInvMasslikeSignPP"), multiplicity, mother.Pt(), mother.M(), cosThetaStarBeam);
          } else if (track1.sign() < 0 && track2.sign() < 0) {
            hInvMass.fill(HIST("h3KstarInvMasslikeSignMM"), multiplicity, mother.Pt(), mother.M(), cosThetaStarBeam);
          }
        }
      }
    } else if (activateTHnSparseCosThStarRandom) {
      auto phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
      auto thetaRandom = gRandom->Uniform(0.f, constants::math::PI);

      randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
      auto cosThetaStarRandom = randomVec.Dot(fourVecDauCM.Vect()) / std::sqrt(fourVecDauCM.Vect().Mag2());

      if (track1.sign() * track2.sign() < 0) {
        if (!isMix) {
          if (std::abs(mother.Rapidity()) < selectionConfig.rapidityMotherData) {
            hInvMass.fill(HIST("h3KstarInvMassUnlikeSign"), multiplicity, mother.Pt(), mother.M(), cosThetaStarRandom);
          }
          for (int i = 0; i < cRotations; i++) {
            theta2 = rn->Uniform(0, o2::constants::math::PI);
            daughterRot = ROOT::Math::PxPyPzMVector(daughter1.Px() * std::cos(theta2) - daughter1.Py() * std::sin(theta2), daughter1.Px() * std::sin(theta2) + daughter1.Py() * std::cos(theta2), daughter1.Pz(), daughter1.M());

            motherRot = daughterRot + daughter2;
            if (calcRotational && std::abs(motherRot.Rapidity()) < selectionConfig.rapidityMotherData)
              hInvMass.fill(HIST("h3KstarInvMassRotated"), multiplicity, motherRot.Pt(), motherRot.M(), cosThetaStarRandom);
          }
        } else if (isMix && std::abs(mother.Rapidity()) < selectionConfig.rapidityMotherData) {
          hInvMass.fill(HIST("h3KstarInvMassMixed"), multiplicity, mother.Pt(), mother.M(), cosThetaStarRandom);
        }
      } else {
        if (!isMix) {
          if (calcLikeSign && std::abs(mother.Rapidity()) < selectionConfig.rapidityMotherData) {
            if (track1.sign() > 0 && track2.sign() > 0) {
              hInvMass.fill(HIST("h3KstarInvMasslikeSignPP"), multiplicity, mother.Pt(), mother.M(), cosThetaStarRandom);
            } else if (track1.sign() < 0 && track2.sign() < 0) {
              hInvMass.fill(HIST("h3KstarInvMasslikeSignMM"), multiplicity, mother.Pt(), mother.M(), cosThetaStarRandom);
            }
          }
        }
      }
    }
    // }
  }

  // int counter = 0;
  void processSE(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)

  {
    int occupancy = collision.trackOccupancyInTimeRange();
    rEventSelection.fill(HIST("hOccupancy"), occupancy);

    if (!selectionEvent(collision, true)) { // fill data event cut histogram
      return;
    }

    multiplicity = -1;

    if (cSelectMultEstimator == kFT0M) {
      multiplicity = collision.centFT0M();
    } else if (cSelectMultEstimator == kFT0A) {
      multiplicity = collision.centFT0A();
    } else if (cSelectMultEstimator == kFT0C) {
      multiplicity = collision.centFT0C();
    } else if (cSelectMultEstimator == kFV0A) {
      multiplicity = collision.centFV0A();
    } else {
      multiplicity = collision.centFT0M(); // default
    }

    // Fill the event counter
    if (cQAevents) {
      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
      rEventSelection.fill(HIST("hMultiplicity"), multiplicity);
    }

    for (const auto& [track1, track2] : combinations(CombinationsFullIndexPolicy(tracks, tracks))) {
      rEventSelection.fill(HIST("tracksCheckData"), 0.5);
      if (!selectionTrack(track1)) {
        continue;
      }
      if (!selectionTrack(track2)) {
        continue;
      }
      rEventSelection.fill(HIST("tracksCheckData"), 1.5);

      if (track1.globalIndex() == track2.globalIndex())
        continue;

      if (cQAplots) {
        hPID.fill(HIST("Before/hNsigma_TPC_TOF_Ka_before"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
        hPID.fill(HIST("Before/hNsigma_TPC_TOF_Pi_before"), track1.tpcNSigmaPi(), track1.tofNSigmaPi(), track1.pt());

        hPID.fill(HIST("Before/hTPCnsigKa_mult_pt"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
        hPID.fill(HIST("Before/hTPCnsigPi_mult_pt"), track1.tpcNSigmaPi(), multiplicity, track1.pt());
        hPID.fill(HIST("Before/hTOFnsigKa_mult_pt"), track1.tofNSigmaKa(), multiplicity, track1.pt());
        hPID.fill(HIST("Before/hTOFnsigPi_mult_pt"), track1.tofNSigmaPi(), multiplicity, track1.pt());

        hOthers.fill(HIST("hCRFC_before"), track1.tpcCrossedRowsOverFindableCls());
        // hOthers.fill(HIST("dE_by_dx_TPC"), track1.p(), track1.tpcSignal());
        hOthers.fill(HIST("hphi"), track1.phi());
      }

      if (cPIDAdditionalChecks) {
        if (std::abs(track1.tpcNSigmaPi()) > selectionConfig.pidCutTemp)
          hPID.fill(HIST("AdditionalChecks/hTPCnsigKa_pt"), track1.tpcNSigmaKa(), track1.pt());
        if (std::abs(track2.tpcNSigmaKa()) > selectionConfig.pidCutTemp)
          hPID.fill(HIST("AdditionalChecks/hTPCnsigPi_pt"), track1.tpcNSigmaPi(), track1.pt());
        if (std::abs(track1.tpcNSigmaPr()) > selectionConfig.pidCutTemp) {
          hPID.fill(HIST("AdditionalChecks/hTOFnsigKa_ptPrCut"), track1.tpcNSigmaKa(), track1.pt());
          hPID.fill(HIST("AdditionalChecks/hTOFnsigPi_ptPrCut"), track1.tpcNSigmaPi(), track1.pt());
        }
      }

      if (cQAevents) {
        rEventSelection.fill(HIST("hDcaxy_cent_pt"), track1.dcaXY(), multiplicity, track1.pt());
        rEventSelection.fill(HIST("hDcaz_cent_pt"), track1.dcaZ(), multiplicity, track1.pt());
      }

      // since we are using combinations full index policy, so repeated pairs are allowed, so we can check one with Kaon and other with pion
      if (!selectionConfig.isapplypTdepPID && !selectionConfig.isapplypTdepPIDTOF && !selectionPID(track1, 1)) // Track 1 is checked with Kaon
        continue;
      if (!selectionConfig.isapplypTdepPID && !selectionConfig.isapplypTdepPIDTOF && !selectionPID(track2, 0)) // Track 2 is checked with Pion
        continue;

      if (selectionConfig.isapplypTdepPID && !selectionConfig.isapplypTdepPIDTOF && !selectionPIDPtDep(track1, 1)) // Track 1 is checked with Kaon
        continue;
      if (selectionConfig.isapplypTdepPID && !selectionConfig.isapplypTdepPIDTOF && !selectionPIDPtDep(track2, 0)) // Track 2 is checked with Pion
        continue;

      if (!selectionConfig.isapplypTdepPID && selectionConfig.isapplypTdepPIDTOF && !selectionPIDptDepTOF(track1, 1)) // Track 1 is checked with Kaon
        continue;
      if (!selectionConfig.isapplypTdepPID && selectionConfig.isapplypTdepPIDTOF && !selectionPIDptDepTOF(track2, 0)) // Track 2 is checked with Pion
        continue;

      rEventSelection.fill(HIST("tracksCheckData"), 2.5);

      // if (cFakeTrack && isFakeTrack(track1, 1)) // Kaon
      //   continue;
      // if (cFakeTrack && isFakeTrack(track2, 0)) // Pion
      //   continue;
      rEventSelection.fill(HIST("tracksCheckData"), 3.5);

      if (std::abs(track1.rapidity(o2::track::PID::getMass(o2::track::PID::Kaon))) > selectionConfig.ctrackRapidity)
        continue;

      if (std::abs(track2.rapidity(o2::track::PID::getMass(o2::track::PID::Pion))) > selectionConfig.ctrackRapidity)
        continue;
      rEventSelection.fill(HIST("tracksCheckData"), 4.5);

      if (selectionConfig.isApplyParticleMID) {
        if (selectionMID(track1, 0)) // Kaon misidentified as pion
          continue;
        // if (selectionMID(track1, 2)) // Kaon misidentified as proton
        //   continue;
        if (selectionMID(track2, 1)) // Pion misidentified as kaon
          continue;
        // if (selectionMID(track2, 2)) // Pion misidentified as proton
        //   continue;
      }

      if (selectionConfig.isApplyParticleMIDPtDep) {
        if (selectionMIDPtDep(track1, 0)) // Kaon misidentified as pion
          continue;
        if (selectionMIDPtDep(track2, 1)) // Pion misidentified as kaon
          continue;
      }

      rEventSelection.fill(HIST("tracksCheckData"), 5.5);

      if (cQAplots) {
        hPID.fill(HIST("After/hDcaxyPi"), track2.dcaXY());
        hPID.fill(HIST("After/hDcaxyKa"), track1.dcaXY());
        hPID.fill(HIST("After/hDcazPi"), track2.dcaZ());
        hPID.fill(HIST("After/hDcazKa"), track1.dcaZ());

        hPID.fill(HIST("After/hTPCnsigKa_mult_pt"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
        hPID.fill(HIST("After/hTPCnsigPi_mult_pt"), track2.tpcNSigmaPi(), multiplicity, track2.pt());
        hPID.fill(HIST("After/hTOFnsigKa_mult_pt"), track1.tofNSigmaKa(), multiplicity, track1.pt());
        hPID.fill(HIST("After/hTOFnsigPi_mult_pt"), track2.tofNSigmaPi(), multiplicity, track2.pt());
        hOthers.fill(HIST("hEta_after"), track1.eta());
        hOthers.fill(HIST("hCRFC_after"), track1.tpcCrossedRowsOverFindableCls());
        hPID.fill(HIST("After/hNsigma_TPC_TOF_Ka_after"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
        hPID.fill(HIST("After/hNsigma_TPC_TOF_Pi_after"), track2.tpcNSigmaPi(), track2.tofNSigmaPi(), track2.pt());
      }

      if (!selectionPair(track1, track2)) {
        continue;
      }
      rEventSelection.fill(HIST("tracksCheckData"), 6.5);

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
      fillInvMass(daughter1, daughter2, mother, multiplicity, isMix, track1, track2);
    }
  }

  PROCESS_SWITCH(Kstarqa, processSE, "Process Same event", true);

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for ME mixing"};
  // ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {10, 0, 100}, "multiplicity percentile for ME mixing"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {2000, 0, 10000}, "TPC multiplicity axis for ME mixing"};

  // using BinningTypeTPCMultiplicity = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  using BinningTypeFT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypeFT0A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0A>;
  using BinningTypeFT0C = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  using BinningTypeFV0A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFV0A>;

  using BinningTypeMCFT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypeMCFT0A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0A>;
  using BinningTypeMCFT0C = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  using BinningTypeMCFV0A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFV0A>;

  BinningTypeFT0M binningOnFT0M{{axisVertex, axisMultiplicity}, true};
  BinningTypeFT0A binningOnFT0A{{axisVertex, axisMultiplicity}, true};
  BinningTypeFT0C binningOnFT0C{{axisVertex, axisMultiplicity}, true};
  BinningTypeFV0A binningOnFV0A{{axisVertex, axisMultiplicity}, true};

  BinningTypeMCFT0M binningOnMCFT0M{{axisVertex, axisMultiplicity}, true};
  BinningTypeMCFT0A binningOnMCFT0A{{axisVertex, axisMultiplicity}, true};
  BinningTypeMCFT0C binningOnMCFT0C{{axisVertex, axisMultiplicity}, true};
  BinningTypeMCFV0A binningOnMCFV0A{{axisVertex, axisMultiplicity}, true};

  SameKindPair<EventCandidates, TrackCandidates, BinningTypeFT0M> pair1{binningOnFT0M, selectionConfig.cfgNoMixedEvents, -1, &cache};
  SameKindPair<EventCandidates, TrackCandidates, BinningTypeFT0A> pair2{binningOnFT0A, selectionConfig.cfgNoMixedEvents, -1, &cache};
  SameKindPair<EventCandidates, TrackCandidates, BinningTypeFT0C> pair3{binningOnFT0C, selectionConfig.cfgNoMixedEvents, -1, &cache};
  SameKindPair<EventCandidates, TrackCandidates, BinningTypeFV0A> pair4{binningOnFV0A, selectionConfig.cfgNoMixedEvents, -1, &cache};

  SameKindPair<EventCandidatesMC, TrackCandidatesMC, BinningTypeMCFT0M> pairmc1{binningOnMCFT0M, selectionConfig.cfgNoMixedEvents, -1, &cache};
  SameKindPair<EventCandidatesMC, TrackCandidatesMC, BinningTypeMCFT0A> pairmc2{binningOnMCFT0A, selectionConfig.cfgNoMixedEvents, -1, &cache};
  SameKindPair<EventCandidatesMC, TrackCandidatesMC, BinningTypeMCFT0C> pairmc3{binningOnMCFT0C, selectionConfig.cfgNoMixedEvents, -1, &cache};
  SameKindPair<EventCandidatesMC, TrackCandidatesMC, BinningTypeMCFV0A> pairmc4{binningOnMCFV0A, selectionConfig.cfgNoMixedEvents, -1, &cache};

  void processME(EventCandidatesMix const&, TrackCandidates const&)
  {
    // Map estimator to pair and multiplicity accessor
    auto runMixing = [&](auto& pair, auto multiplicityGetter) {
      for (const auto& [c1, tracks1, c2, tracks2] : pair) {
        // if (!c1.sel8() || !c2.sel8())
        //   continue;

        if (!selectionEvent(c1, false) || !selectionEvent(c2, false)) { // don't fill event cut histogram
          continue;
        }

        multiplicity = multiplicityGetter(c1);

        for (const auto& [t1, t2] : o2::soa::combinations(
               o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
          if (!selectionTrack(t1) || !selectionTrack(t2))
            continue;
          // if (!selectionPID(t1, 1) || !selectionPID(t2, 0))
          //   continue;
          if (!selectionConfig.isapplypTdepPID && !selectionConfig.isapplypTdepPIDTOF && !selectionPID(t1, 1)) // Track 1 is checked with Kaon
            continue;
          if (!selectionConfig.isapplypTdepPID && !selectionConfig.isapplypTdepPIDTOF && !selectionPID(t2, 0)) // Track 2 is checked with Pion
            continue;

          if (selectionConfig.isapplypTdepPID && !selectionConfig.isapplypTdepPIDTOF && !selectionPIDPtDep(t1, 1)) // Track 1 is checked with Kaon
            continue;
          if (selectionConfig.isapplypTdepPID && !selectionConfig.isapplypTdepPIDTOF && !selectionPIDPtDep(t2, 0)) // Track 2 is checked with Pion
            continue;

          if (!selectionConfig.isapplypTdepPID && selectionConfig.isapplypTdepPIDTOF && !selectionPIDptDepTOF(t1, 1)) // Track 1 is checked with Kaon
            continue;
          if (!selectionConfig.isapplypTdepPID && selectionConfig.isapplypTdepPIDTOF && !selectionPIDptDepTOF(t2, 0)) // Track 2 is checked with Pion
            continue;

          if (std::abs(t1.rapidity(o2::track::PID::getMass(o2::track::PID::Kaon))) > selectionConfig.ctrackRapidity)
            continue;

          if (std::abs(t2.rapidity(o2::track::PID::getMass(o2::track::PID::Pion))) > selectionConfig.ctrackRapidity)
            continue;

          // if (cFakeTrack && isFakeTrack(t1, 1)) // Kaon
          //   continue;
          // if (cFakeTrack && isFakeTrack(t2, 0)) // Pion
          //   continue;

          if (selectionConfig.isApplyParticleMID) {
            if (selectionMID(t1, 0)) // Kaon misidentified as pion
              continue;
            // if (selectionMID(t1, 2)) // Kaon misidentified as proton
            //   continue;
            if (selectionMID(t2, 1)) // Pion misidentified as kaon
              continue;
            // if (selectionMID(t2, 2)) // Pion misidentified as proton
            //   continue;
          }

          if (selectionConfig.isApplyParticleMIDPtDep) {
            if (selectionMIDPtDep(t1, 0)) // Kaon misidentified as pion
              continue;
            if (selectionMIDPtDep(t2, 1)) // Pion misidentified as kaon
              continue;
          }

          if (!selectionPair(t1, t2)) {
            continue;
          }

          daughter1 = ROOT::Math::PxPyPzMVector(t1.px(), t1.py(), t1.pz(), massKa);
          daughter2 = ROOT::Math::PxPyPzMVector(t2.px(), t2.py(), t2.pz(), massPi);
          mother = daughter1 + daughter2;

          isMix = true;

          if (std::abs(mother.Rapidity()) < selectionConfig.rapidityMotherData) {
            fillInvMass(daughter1, daughter2, mother, multiplicity, isMix, t1, t2);
          }
        }
      }
    };

    // Call mixing based on selected estimator
    if (cSelectMultEstimator == kFT0M) {
      runMixing(pair1, [](const auto& c) { return c.centFT0M(); });
    } else if (cSelectMultEstimator == kFT0A) {
      runMixing(pair2, [](const auto& c) { return c.centFT0A(); });
    } else if (cSelectMultEstimator == kFT0C) {
      runMixing(pair3, [](const auto& c) { return c.centFT0C(); });
    } else if (cSelectMultEstimator == kFV0A) {
      runMixing(pair4, [](const auto& c) { return c.centFV0A(); });
    }
  }
  PROCESS_SWITCH(Kstarqa, processME, "Process Mixed event", true);

  void processMEMC(EventCandidatesMC const&, TrackCandidatesMC const&, aod::McParticles const&, aod::McCollisions const&)
  {
    auto runMixing = [&](auto& pair, auto multiplicityGetter) {
      for (const auto& [c1, tracks1, c2, tracks2] : pair) {

        if (!selectionEvent(c1, false) || !selectionEvent(c2, false)) { // don't fill event cut histogram
          continue;
        }

        if (selectionConfig.isApplyMCchecksClosure && (!c1.has_mcCollision() || !c2.has_mcCollision())) {
          continue; // skip if no MC collision associated
        }

        multiplicity = multiplicityGetter(c1);
        // multiplicity = c1.centFT0M(); // default, can be changed later

        for (const auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
          if (!selectionTrack(t1) || !selectionTrack(t2))
            continue;
          // if (!selectionPID(t1, 1) || !selectionPID(t2, 0))
          //   continue;
          if (!selectionConfig.isapplypTdepPID && !selectionPID(t1, 1)) // Track 1 is checked with Kaon
            continue;
          if (!selectionConfig.isapplypTdepPID && !selectionPID(t2, 0)) // Track 2 is checked with Pion
            continue;

          if (selectionConfig.isapplypTdepPID && !selectionPIDPtDep(t1, 1)) // Track 1 is checked with Kaon
            continue;
          if (selectionConfig.isapplypTdepPID && !selectionPIDPtDep(t2, 0)) // Track 2 is checked with Pion
            continue;

          if (selectionConfig.isApplyParticleMID) {
            if (selectionMID(t1, 0)) // Kaon misidentified as pion
              continue;
            // if (selectionMID(t1, 2)) // Kaon misidentified as proton
            //   continue;
            if (selectionMID(t2, 1)) // Pion misidentified as kaon
              continue;
            // if (selectionMID(t2, 2)) // Pion misidentified as proton
            //   continue;
          }

          if (selectionConfig.isApplyParticleMIDPtDep) {
            if (selectionMIDPtDep(t1, 0)) // Kaon misidentified as pion
              continue;
            if (selectionMIDPtDep(t2, 1)) // Pion misidentified as kaon
              continue;
          }

          if (selectionConfig.isApplyMCchecksClosure && (!t1.has_mcParticle() || !t2.has_mcParticle())) {
            continue; // skip if no MC particle associated
          }

          if (std::abs(t1.rapidity(o2::track::PID::getMass(o2::track::PID::Kaon))) > selectionConfig.ctrackRapidity)
            continue;

          if (std::abs(t2.rapidity(o2::track::PID::getMass(o2::track::PID::Pion))) > selectionConfig.ctrackRapidity)
            continue;

          const auto mctrack1 = t1.mcParticle();
          const auto mctrack2 = t2.mcParticle();

          if (selectionConfig.isApplyMCchecksClosure && !mctrack1.isPhysicalPrimary()) {
            continue;
          }

          if (selectionConfig.isApplyMCchecksClosure && !mctrack2.isPhysicalPrimary()) {
            continue;
          }

          daughter1 = ROOT::Math::PxPyPzMVector(t1.px(), t1.py(), t1.pz(), massKa);
          daughter2 = ROOT::Math::PxPyPzMVector(t2.px(), t2.py(), t2.pz(), massPi);
          mother = daughter1 + daughter2;

          isMix = true;

          if (std::abs(mother.Rapidity()) < selectionConfig.rapidityMotherData) {
            fillInvMass(daughter1, daughter2, mother, multiplicity, isMix, t1, t2);
          }
        }
      }
    };
    // Call mixing based on selected estimator
    if (cSelectMultEstimator == kFT0M) {
      runMixing(pairmc1, [](const auto& c) { return c.centFT0M(); });
    } else if (cSelectMultEstimator == kFT0A) {
      runMixing(pairmc2, [](const auto& c) { return c.centFT0A(); });
    } else if (cSelectMultEstimator == kFT0C) {
      runMixing(pairmc3, [](const auto& c) { return c.centFT0C(); });
    } else if (cSelectMultEstimator == kFV0A) {
      runMixing(pairmc4, [](const auto& c) { return c.centFV0A(); });
    }
  }
  PROCESS_SWITCH(Kstarqa, processMEMC, "Process mixed-event in MC", false);

  void processSEMC(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&, aod::McCollisions const& /*mcCollisions*/)
  {
    auto oldindex = -999;
    if (!collision.has_mcCollision()) {
      return;
    }
    int occupancy = collision.trackOccupancyInTimeRange();
    rEventSelection.fill(HIST("hOccupancy"), occupancy);

    if (!selectionEvent(collision, false)) { // don't fill event cut histogram
      return;
    }

    multiplicity = -1;

    if (cSelectMultEstimator == kFT0M) {
      multiplicity = collision.centFT0M();
    } else if (cSelectMultEstimator == kFT0A) {
      multiplicity = collision.centFT0A();
    } else if (cSelectMultEstimator == kFT0C) {
      multiplicity = collision.centFT0C();
    } else if (cSelectMultEstimator == kFV0A) {
      multiplicity = collision.centFV0A();
    } else {
      multiplicity = collision.centFT0M(); // default
    }

    // Fill the event counter
    if (cQAevents) {
      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
      rEventSelection.fill(HIST("hMultiplicity"), multiplicity);
    }

    for (const auto& [track1, track2] : combinations(CombinationsFullIndexPolicy(tracks, tracks))) {
      rEventSelection.fill(HIST("tracksCheckData"), 0.5);
      if (!selectionTrack(track1)) {
        continue;
      }
      if (!selectionTrack(track2)) {
        continue;
      }

      const auto mctrack1 = track1.mcParticle();
      const auto mctrack2 = track2.mcParticle();

      if (selectionConfig.isApplyMCchecksClosure && (!track1.has_mcParticle() || !track2.has_mcParticle())) {
        continue; // skip if no MC particle associated
      }

      if (selectionConfig.isApplyMCchecksClosure && !mctrack1.isPhysicalPrimary()) {
        continue;
      }

      if (selectionConfig.isApplyMCchecksClosure && !mctrack2.isPhysicalPrimary()) {
        continue;
      }
      rEventSelection.fill(HIST("tracksCheckData"), 1.5);

      if (cQAplots) {
        hPID.fill(HIST("Before/hNsigma_TPC_TOF_Ka_before"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
        hPID.fill(HIST("Before/hNsigma_TPC_TOF_Pi_before"), track1.tpcNSigmaPi(), track1.tofNSigmaPi(), track1.pt());

        hPID.fill(HIST("Before/hTPCnsigKa_mult_pt"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
        hPID.fill(HIST("Before/hTPCnsigPi_mult_pt"), track1.tpcNSigmaPi(), multiplicity, track1.pt());
        hPID.fill(HIST("Before/hTOFnsigKa_mult_pt"), track1.tofNSigmaKa(), multiplicity, track1.pt());
        hPID.fill(HIST("Before/hTOFnsigPi_mult_pt"), track1.tofNSigmaPi(), multiplicity, track1.pt());

        hOthers.fill(HIST("hCRFC_before"), track1.tpcCrossedRowsOverFindableCls());
        hOthers.fill(HIST("hphi"), track1.phi());
      }

      if (cQAevents) {
        rEventSelection.fill(HIST("hDcaxy_cent_pt"), track1.dcaXY(), multiplicity, track1.pt());
        rEventSelection.fill(HIST("hDcaz_cent_pt"), track1.dcaZ(), multiplicity, track1.pt());
      }

      // since we are using combinations full index policy, so repeated pairs are allowed, so we can check one with Kaon and other with pion
      if (!selectionConfig.isapplypTdepPID && !selectionPID(track1, 1)) // Track 1 is checked with Kaon
        continue;
      if (!selectionConfig.isapplypTdepPID && !selectionPID(track2, 0)) // Track 2 is checked with Pion
        continue;
      rEventSelection.fill(HIST("tracksCheckData"), 2.5);

      if (selectionConfig.isapplypTdepPID && !selectionPIDPtDep(track1, 1)) // Track 1 is checked with Kaon
        continue;
      if (selectionConfig.isapplypTdepPID && !selectionPIDPtDep(track2, 0)) // Track 2 is checked with Pion
        continue;
      rEventSelection.fill(HIST("tracksCheckData"), 3.5);

      if (selectionConfig.isApplyParticleMID) {
        if (selectionMID(track1, 0)) // Kaon misidentified as pion
          continue;
        // if (selectionMID(track1, 2)) // Kaon misidentified as proton
        //   continue;
        if (selectionMID(track2, 1)) // Pion misidentified as kaon
          continue;
        // if (selectionMID(track2, 2)) // Pion misidentified as proton
        //   continue;
      }

      if (selectionConfig.isApplyParticleMIDPtDep) {
        if (selectionMIDPtDep(track1, 0)) // Kaon misidentified as pion
          continue;
        if (selectionMIDPtDep(track2, 1)) // Pion misidentified as kaon
          continue;
      }

      rEventSelection.fill(HIST("tracksCheckData"), 4.5);
      if (std::abs(track1.rapidity(o2::track::PID::getMass(o2::track::PID::Kaon))) > selectionConfig.ctrackRapidity)
        continue;

      if (std::abs(track2.rapidity(o2::track::PID::getMass(o2::track::PID::Pion))) > selectionConfig.ctrackRapidity)
        continue;

      // if (cFakeTrack && isFakeTrack(track1, 1)) // Kaon
      //   continue;
      // if (cFakeTrack && isFakeTrack(track2, 0)) // Pion
      //   continue;
      rEventSelection.fill(HIST("tracksCheckData"), 5.5);

      if (track1.globalIndex() == track2.globalIndex())
        continue;

      rEventSelection.fill(HIST("tracksCheckData"), 6.5);

      if (cQAplots) {
        hPID.fill(HIST("After/hDcaxyPi"), track2.dcaXY());
        hPID.fill(HIST("After/hDcaxyKa"), track1.dcaXY());
        hPID.fill(HIST("After/hDcazPi"), track2.dcaZ());
        hPID.fill(HIST("After/hDcazKa"), track1.dcaZ());

        hPID.fill(HIST("After/hTPCnsigKa_mult_pt"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
        hPID.fill(HIST("After/hTPCnsigPi_mult_pt"), track2.tpcNSigmaPi(), multiplicity, track2.pt());
        hPID.fill(HIST("After/hTOFnsigKa_mult_pt"), track1.tofNSigmaKa(), multiplicity, track1.pt());
        hPID.fill(HIST("After/hTOFnsigPi_mult_pt"), track2.tofNSigmaPi(), multiplicity, track2.pt());
        hOthers.fill(HIST("hEta_after"), track1.eta());
        hOthers.fill(HIST("hCRFC_after"), track1.tpcCrossedRowsOverFindableCls());
        // hPID.fill(HIST("After/hNsigmaKaonTPC_after"), track1.pt(), track1.tpcNSigmaKa());
        // hPID.fill(HIST("After/hNsigmaKaonTOF_after"), track1.pt(), track1.tofNSigmaKa());
        // hPID.fill(HIST("After/hNsigmaPionTPC_after"), track2.pt(), track2.tpcNSigmaPi());
        // hPID.fill(HIST("After/hNsigmaPionTOF_after"), track2.pt(), track2.tofNSigmaPi());
        hPID.fill(HIST("After/hNsigma_TPC_TOF_Ka_after"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
        hPID.fill(HIST("After/hNsigma_TPC_TOF_Pi_after"), track2.tpcNSigmaPi(), track2.tofNSigmaPi(), track2.pt());
      }

      if (selectionConfig.isApplyMCchecksClosure) {
        for (const auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
          for (const auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {

            if (mothertrack1.globalIndex() != mothertrack2.globalIndex()) {
              continue;
            }

            if (selectionConfig.isApplyMCchecksClosure && !mothertrack1.producedByGenerator()) {
              continue;
            }

            if (selectionConfig.isApplyMCchecksClosure && avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 11.5);

            oldindex = mothertrack1.globalIndex();

            daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
            daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
            mother = daughter1 + daughter2; // Kstar meson

            hOthers.fill(HIST("hKstar_rap_pt"), mother.Rapidity(), mother.Pt());
            hOthers.fill(HIST("hKstar_eta_pt"), mother.Eta(), mother.Pt());

            isMix = false;
            fillInvMass(daughter1, daughter2, mother, multiplicity, isMix, track1, track2);
          }
        }
      } else {
        daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
        mother = daughter1 + daughter2; // Kstar meson

        hOthers.fill(HIST("hKstar_rap_pt"), mother.Rapidity(), mother.Pt());
        hOthers.fill(HIST("hKstar_eta_pt"), mother.Eta(), mother.Pt());

        isMix = false;
        fillInvMass(daughter1, daughter2, mother, multiplicity, isMix, track1, track2);
      }
    }
  }
  PROCESS_SWITCH(Kstarqa, processSEMC, "Process same event in MC", false);

  Service<o2::framework::O2DatabasePDG> pdgDB;

  void processGen(EventMCGenerated::iterator const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<EventCandidatesMC>& collisions)
  // void processGen(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<EventCandidatesMC>& collisions)
  {
    rEventSelection.fill(HIST("eventsCheckGen"), 0.5);

    int nChInel = 0;
    for (const auto& mcParticle : mcParticles) {
      auto pdgcode = std::abs(mcParticle.pdgCode());
      if (mcParticle.isPhysicalPrimary() && (pdgcode == PDG_t::kPiPlus || pdgcode == PDG_t::kKPlus || pdgcode == PDG_t::kProton || pdgcode == std::abs(PDG_t::kElectron) || pdgcode == std::abs(PDG_t::kMuonMinus))) {
        if (std::abs(mcParticle.eta()) < 1.0) {
          nChInel = nChInel + 1;
        }
      }
    }
    if (nChInel > 0 && std::abs(mcCollision.posZ()) < selectionConfig.cutzvertex)
      rEventSelection.fill(HIST("eventsCheckGen"), 1.5);

    std::vector<int64_t> selectedEvents(collisions.size());
    int nevts = 0;
    multiplicity = -1.0;
    // float impactParameter = mcCollision.impactParameter();

    bool isINELgt0true = false;

    if (pwglf::isINELgtNmc(mcParticles, 0, pdgDB)) {
      isINELgt0true = true;
    }
    if (selectionConfig.isINELgt0 && !isINELgt0true) {
      return;
    }

    // if (selectionConfig.isINELgt0 && !mcCollision.isInelGt0()) {
    //   return;
    // }
    rEventSelection.fill(HIST("eventsCheckGen"), 2.5);

    for (const auto& collision : collisions) {
      rEventSelection.fill(HIST("eventsCheckGen"), 3.5);

      if (!selectionEvent(collision, false)) { // don't fill event cut histogram
        continue;
      }
      multiplicity = collision.centFT0M();

      if (cSelectMultEstimator == kFT0M) {
        multiplicity = collision.centFT0M();
      } else if (cSelectMultEstimator == kFT0A) {
        multiplicity = collision.centFT0A();
      } else if (cSelectMultEstimator == kFT0C) {
        multiplicity = collision.centFT0C();
      } else if (cSelectMultEstimator == kFV0A) {
        multiplicity = collision.centFV0A();
      } else {
        multiplicity = collision.centFT0M(); // default
      }
      hInvMass.fill(HIST("h1GenMult"), multiplicity);

      int occupancy = collision.trackOccupancyInTimeRange();
      rEventSelection.fill(HIST("hOccupancy"), occupancy);

      selectedEvents[nevts++] = collision.mcCollision_as<EventMCGenerated>().globalIndex();
    }
    selectedEvents.resize(nevts);

    for (const auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.y()) < selectionConfig.rapidityMotherData && std::abs(mcParticle.pdgCode()) == o2::constants::physics::kK0Star892) {
        hInvMass.fill(HIST("hAllKstarGenCollisisons"), multiplicity, mcParticle.pt());
      }
    }

    if (std::abs(mcCollision.posZ()) > selectionConfig.cutzvertex) {
      return;
    }

    const auto evtReconstructedAndSelected = std::find(selectedEvents.begin(), selectedEvents.end(), mcCollision.globalIndex()) != selectedEvents.end();
    hInvMass.fill(HIST("hAllGenCollisions"), multiplicity);
    // if (!evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
    //   return;
    // }
    double genMultiplicity = mcCollision.centFT0M();
    hInvMass.fill(HIST("sigEvLossFromGenRec/MultiplicityGen"), genMultiplicity);
    if (evtReconstructedAndSelected) {
      hInvMass.fill(HIST("sigEvLossFromGenRec/MultiplicityRec"), genMultiplicity);
      hInvMass.fill(HIST("h1GenMult2"), genMultiplicity);
      hInvMass.fill(HIST("hAllGenCollisions1Rec"), multiplicity);
      rEventSelection.fill(HIST("eventsCheckGen"), 4.5);
    }

    for (const auto& mcParticle : mcParticles) {

      if (std::abs(mcParticle.y()) >= selectionConfig.rapidityMotherData) {
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
      if (evtReconstructedAndSelected)
        hInvMass.fill(HIST("hAllKstarGenCollisisons1Rec"), multiplicity, mcParticle.pt());

      auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
      if (kDaughters.size() != selectionConfig.noOfDaughters) {
        continue;
      }

      auto passkaon = false;
      auto passpion = false;
      for (const auto& kCurrentDaughter : kDaughters) {
        if (!kCurrentDaughter.isPhysicalPrimary()) {
          continue;
        }

        if (std::abs(kCurrentDaughter.pdgCode()) == PDG_t::kKPlus) {
          passkaon = true;
          daughter1 = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);

        } else if (std::abs(kCurrentDaughter.pdgCode()) == PDG_t::kPiPlus) {
          passpion = true;
          daughter2 = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massPi);
        }
      }
      if (passkaon && passpion) {
        mother = daughter1 + daughter2; // Kstar meson
        hInvMass.fill(HIST("sigEvLossFromGenRec/hSignalLossDenominator"), mother.Pt(), genMultiplicity);
        if (evtReconstructedAndSelected) {
          hInvMass.fill(HIST("hk892GenpT"), mcParticle.pt(), multiplicity);
          hInvMass.fill(HIST("hk892GenpT2"), mother.Pt(), multiplicity);
          hInvMass.fill(HIST("hk892GenpTCalib1"), mcParticle.pt(), genMultiplicity);
          hInvMass.fill(HIST("hk892GenpTCalib2"), mother.Pt(), genMultiplicity);
          hInvMass.fill(HIST("h1genmass"), mother.M());
          hInvMass.fill(HIST("sigEvLossFromGenRec/hSignalLossNumerator"), mother.Pt(), genMultiplicity);
        }
      }
    }
  }
  PROCESS_SWITCH(Kstarqa, processGen, "Process Generated", false);

  void processEvtLossSigLossMC(EventMCGenerated::iterator const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<EventCandidatesMC>& recCollisions)
  // void processEvtLossSigLossMC(aod::McCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<EventCandidatesMC>& recCollisions)
  {
    // auto multiplicityRec = -1;
    bool isSelectedEvent = false;
    auto multiplicity1 = -999.;

    for (const auto& RecCollision : recCollisions) {
      if (!RecCollision.has_mcCollision())
        continue;
      if (!selectionEvent(RecCollision, false)) // don't fill event cut histogram
        continue;

      // const auto& mcCollisionRec = RecCollision.mcCollision_as<EventMCGenerated>();
      // multiplicityRec = mcCollisionRec.centFT0M();

      if (cSelectMultEstimator == kFT0M) {
        multiplicity1 = RecCollision.centFT0M();
      } else if (cSelectMultEstimator == kFT0A) {
        multiplicity1 = RecCollision.centFT0A();
      } else if (cSelectMultEstimator == kFT0C) {
        multiplicity1 = RecCollision.centFT0C();
      } else if (cSelectMultEstimator == kFV0A) {
        multiplicity1 = RecCollision.centFV0A();
      } else {
        multiplicity1 = RecCollision.centFT0M(); // default
      }

      isSelectedEvent = true;
    }

    bool isINELgt0true = false;

    if (pwglf::isINELgtNmc(mcParticles, 0, pdgDB)) {
      isINELgt0true = true;
    }
    if (selectionConfig.isINELgt0 && !isINELgt0true) {
      return;
    }

    if (selectionConfig.checkVzEvSigLoss && (std::abs(mcCollision.posZ()) > selectionConfig.cutzvertex)) {
      return;
    }

    auto impactPar = mcCollision.impactParameter();
    auto multiplicityGen = -1;
    multiplicityGen = mcCollision.centFT0M();
    hInvMass.fill(HIST("MCcorrections/hImpactParameterGen"), impactPar);
    hInvMass.fill(HIST("MCcorrections/MultiplicityGen"), multiplicityGen);

    // Event loss
    if (isSelectedEvent) {
      hInvMass.fill(HIST("MCcorrections/hImpactParameterRec"), impactPar);
      hInvMass.fill(HIST("MCcorrections/MultiplicityRec"), multiplicityGen);
      hInvMass.fill(HIST("MCcorrections/hImpactParametervsMultiplicity"), impactPar, multiplicity1);
    }

    // Generated MC
    for (const auto& mcPart : mcParticles) {
      if (std::abs(mcPart.y()) >= selectionConfig.rapidityMotherData || std::abs(mcPart.pdgCode()) != o2::constants::physics::kK0Star892)
        continue;

      hInvMass.fill(HIST("MCcorrections/hSignalLossDenominator"), mcPart.pt(), multiplicityGen);
      if (isSelectedEvent) {
        hInvMass.fill(HIST("MCcorrections/hSignalLossNumerator"), mcPart.pt(), multiplicityGen);
      }
    } // end loop on gen particles
  }
  PROCESS_SWITCH(Kstarqa, processEvtLossSigLossMC, "Process Signal Loss, Event Loss", false);

  void processEvSigLossFactors(McCollisionMults::iterator const& mcCollision, soa::SmallGroups<EventCandidatesMC> const& collisions, LabeledTracks const&, aod::McParticles const& mcParticles)
  {
    auto multiplicityNch = -1;
    multiplicityNch = mcCollision.multMCNParticlesEta05();
    hInvMass.fill(HIST("CorrFactors/hGenEvents"), multiplicityNch, 0.5);

    if (selectionConfig.checkVzEvSigLoss && std::abs(mcCollision.posZ()) > selectionConfig.cutzvertex)
      return;

    hInvMass.fill(HIST("CorrFactors/hGenEvents"), multiplicityNch, 1.5);

    if (selectionConfig.isINELgt0 && !mcCollision.isInelGt0()) {
      return;
    }
    hInvMass.fill(HIST("CorrFactors/hGenEvents"), multiplicityNch, 2.5);

    float multiplicity = -1.0;
    bool isSelectedEvent = false;

    for (auto const& collision : collisions) {
      if (!collision.has_mcCollision())
        continue;
      if (!selectionEvent(collision, false)) // don't fill event cut histogram
        continue;

      if (cSelectMultEstimator == kFT0M) {
        multiplicity = collision.centFT0M();
      } else if (cSelectMultEstimator == kFT0A) {
        multiplicity = collision.centFT0A();
      } else if (cSelectMultEstimator == kFT0C) {
        multiplicity = collision.centFT0C();
      } else if (cSelectMultEstimator == kFV0A) {
        multiplicity = collision.centFV0A();
      } else {
        multiplicity = collision.centFT0M(); // default
      }
      isSelectedEvent = true;
    }

    // auto multiplicityGen = -1;
    // multiplicityGen = mcCollision.centFT0M();

    hInvMass.fill(HIST("CorrFactors/hMultiplicityVsMultMC"), multiplicity, multiplicityNch);
    hInvMass.fill(HIST("CorrFactors/hNrecInGen"), collisions.size());
    hInvMass.fill(HIST("CorrFactors/MultiplicityGen"), multiplicityNch);
    if (isSelectedEvent) {
      hInvMass.fill(HIST("CorrFactors/MultiplicityRec"), multiplicityNch);
    }

    for (const auto& mcParticle : mcParticles) {

      if (std::abs(mcParticle.y()) >= selectionConfig.rapidityMotherData)
        continue;

      if (std::abs(mcParticle.pdgCode()) == o2::constants::physics::kK0Star892) {

        auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
        if (kDaughters.size() != selectionConfig.noOfDaughters) {
          continue;
        }

        bool hasPos = false;
        bool hasNeg = false;

        auto passkaon = false;
        auto passpion = false;
        for (const auto& kCurrentDaughter : kDaughters) {
          // if (!kCurrentDaughter.isPhysicalPrimary())
          //   continue;

          int pdgDau = kCurrentDaughter.pdgCode();
          int sign = (pdgDau > 0) - (pdgDau < 0);

          if (sign > 0)
            hasPos = true;
          if (sign < 0)
            hasNeg = true;

          if (std::abs(kCurrentDaughter.pdgCode()) == PDG_t::kKPlus) {
            passkaon = true;
            daughter1 = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);

          } else if (std::abs(kCurrentDaughter.pdgCode()) == PDG_t::kPiPlus) {
            passpion = true;
            daughter2 = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massPi);
          }
        }

        if ((passkaon && passpion) && (hasPos && hasNeg)) {
          mother = daughter1 + daughter2; // Kstar meson

          hInvMass.fill(HIST("CorrFactors/h2dGenKstar"), multiplicity, mother.Pt());
          hInvMass.fill(HIST("CorrFactors/h3dGenKstarVsMultMCVsMultiplicity"), multiplicityNch, multiplicity, mother.Pt());
          hInvMass.fill(HIST("CorrFactors/hSignalLoss1"), mother.pt(), multiplicityNch);
          if (isSelectedEvent) {
            hInvMass.fill(HIST("CorrFactors/hSignalLoss2"), mother.pt(), multiplicityNch);
          }
        }
      }
    }

    if (collisions.size() == 0)
      return;

    hInvMass.fill(HIST("CorrFactors/hGenEvents"), multiplicityNch, 3.5);
  }
  PROCESS_SWITCH(Kstarqa, processEvSigLossFactors, "Process Event and Signal loss", false);

  void processRec(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&, EventMCGenerated const&)
  {

    if (!collision.has_mcCollision()) {
      return;
    }

    double multiplicityRec = -1.0;
    // multiplicityRec = collision.mcCollision_as<EventMCGenerated>().centFT0M();
    const auto& mcCollisionRec = collision.mcCollision_as<EventMCGenerated>();
    multiplicityRec = mcCollisionRec.centFT0M();

    if (selectionConfig.isINELgt0 && !collision.isInelGt0()) {
      return;
    }
    // multiplicity = collision.centFT0M();

    if (cSelectMultEstimator == kFT0M) {
      multiplicity = collision.centFT0M();
    } else if (cSelectMultEstimator == kFT0A) {
      multiplicity = collision.centFT0A();
    } else if (cSelectMultEstimator == kFT0C) {
      multiplicity = collision.centFT0C();
    } else if (cSelectMultEstimator == kFV0A) {
      multiplicity = collision.centFV0A();
    } else {
      multiplicity = collision.centFT0M(); // default
    }

    hInvMass.fill(HIST("hAllRecCollisions"), multiplicity);
    hInvMass.fill(HIST("hAllRecCollisionsCalib"), multiplicityRec);

    if (!selectionEvent(collision, true)) { // fill MC event cut histogram
      return;
    }

    hInvMass.fill(HIST("h1RecMult"), multiplicity);
    hInvMass.fill(HIST("h1RecMult2"), multiplicityRec);

    if (cQAevents) {
      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
    }

    auto oldindex = -999;
    for (const auto& track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }

      if (!track1.has_mcParticle()) {
        continue;
      }

      if (cQAevents) {
        rEventSelection.fill(HIST("hDcaxy_cent_pt"), track1.dcaXY(), multiplicity, track1.pt());
        rEventSelection.fill(HIST("hDcaz_cent_pt"), track1.dcaZ(), multiplicity, track1.pt());
      }

      auto track1ID = track1.index();
      for (const auto& track2 : tracks) {
        rEventSelection.fill(HIST("recMCparticles"), 0.5);
        if (!track2.has_mcParticle()) {
          continue;
        }
        rEventSelection.fill(HIST("recMCparticles"), 1.5);

        if (!selectionTrack(track2)) {
          continue;
        }
        rEventSelection.fill(HIST("recMCparticles"), 2.5);

        auto track2ID = track2.index();
        if (track2ID <= track1ID) {
          continue;
        }
        rEventSelection.fill(HIST("recMCparticles"), 3.5);

        if (track1.sign() * track2.sign() >= 0) {
          continue;
        }
        rEventSelection.fill(HIST("recMCparticles"), 4.5);

        const auto mctrack1 = track1.mcParticle();
        const auto mctrack2 = track2.mcParticle();
        int track1PDG = std::abs(mctrack1.pdgCode());
        int track2PDG = std::abs(mctrack2.pdgCode());

        if (cQAplots && (std::abs(mctrack1.pdgCode()) == PDG_t::kKPlus && std::abs(mctrack2.pdgCode()) == PDG_t::kPiPlus)) {
          hPID.fill(HIST("Before/hNsigma_TPC_TOF_Ka_before"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
          hPID.fill(HIST("Before/hNsigma_TPC_TOF_Pi_before"), track1.tpcNSigmaPi(), track1.tofNSigmaPi(), track1.pt());
          hPID.fill(HIST("Before/hTPCnsigKa_mult_pt"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
          hPID.fill(HIST("Before/hTPCnsigPi_mult_pt"), track1.tpcNSigmaPi(), multiplicity, track1.pt());
          hPID.fill(HIST("Before/hTOFnsigKa_mult_pt"), track1.tofNSigmaKa(), multiplicity, track1.pt());
          hPID.fill(HIST("Before/hTOFnsigPi_mult_pt"), track1.tofNSigmaPi(), multiplicity, track1.pt());
        }

        if (!mctrack1.isPhysicalPrimary()) {
          continue;
        }

        if (!mctrack2.isPhysicalPrimary()) {
          continue;
        }
        rEventSelection.fill(HIST("recMCparticles"), 5.5);
        rEventSelection.fill(HIST("recMCparticles"), 6.5);

        if (selectionConfig.isPDGCheckMC && !(track1PDG == PDG_t::kPiPlus && track2PDG == PDG_t::kKPlus) && !(track1PDG == PDG_t::kKPlus && track2PDG == PDG_t::kPiPlus)) {
          continue;
        }
        rEventSelection.fill(HIST("recMCparticles"), 7.5);

        for (const auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
          for (const auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
            if (selectionConfig.isPDGCheckMC && (mothertrack1.pdgCode() != mothertrack2.pdgCode())) {
              continue;
            }

            if (mothertrack1.globalIndex() != mothertrack2.globalIndex()) {
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 8.5);

            if (!mothertrack1.producedByGenerator()) {
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 9.5);

            if (std::abs(mothertrack1.y()) >= selectionConfig.rapidityMotherData) {
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 10.5);

            if (selectionConfig.isPDGCheckMC && (std::abs(mothertrack1.pdgCode()) != o2::constants::physics::kK0Star892)) {
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 11.5);

            if (selectionConfig.isPDGCheckMC && (track1PDG == PDG_t::kPiPlus)) {
              if (!selectionConfig.isapplypTdepPID && !selectionConfig.isapplypTdepPIDTOF && !(selectionPID(track1, 0) && selectionPID(track2, 1))) { // pion and kaon
                continue;
              } else if (selectionConfig.isapplypTdepPID && !selectionConfig.isapplypTdepPIDTOF && !(selectionPIDPtDep(track1, 0) && selectionPIDPtDep(track2, 1))) { // pion and kaon
                continue;
              } else if (!selectionConfig.isapplypTdepPID && selectionConfig.isapplypTdepPIDTOF && !(selectionPIDptDepTOF(track1, 0) && selectionPIDptDepTOF(track2, 1))) { // pion and kaon
                continue;
              }

              rEventSelection.fill(HIST("recMCparticles"), 12.5);
              if (selectionConfig.isApplyParticleMID) {
                if (selectionMID(track2, 0)) // Kaon misidentified as pion
                  continue;
                // if (selectionMID(track2, 2)) // Kaon misidentified as proton
                //   continue;
                if (selectionMID(track1, 1)) // Pion misidentified as kaon
                  continue;
                // if (selectionMID(track1, 2)) // Pion misidentified as proton
                //   continue;
              }

              if (selectionConfig.isApplyParticleMIDPtDep) {
                if (selectionMIDPtDep(track1, 1)) // Kaon misidentified as pion
                  continue;
                if (selectionMIDPtDep(track2, 0)) // Pion misidentified as kaon
                  continue;
              }
              rEventSelection.fill(HIST("recMCparticles"), 13.5);

              if (std::abs(track1.rapidity(o2::track::PID::getMass(o2::track::PID::Pion))) > selectionConfig.ctrackRapidity)
                continue;

              if (std::abs(track2.rapidity(o2::track::PID::getMass(o2::track::PID::Kaon))) > selectionConfig.ctrackRapidity)
                continue;
              rEventSelection.fill(HIST("recMCparticles"), 14.5);

              if (cQAplots) {
                hPID.fill(HIST("After/hTPCnsigPi_mult_pt"), track1.tpcNSigmaPi(), multiplicity, track1.pt());
                hPID.fill(HIST("After/hTPCnsigKa_mult_pt"), track2.tpcNSigmaKa(), multiplicity, track2.pt());
                hPID.fill(HIST("After/hTOFnsigPi_mult_pt"), track1.tofNSigmaPi(), multiplicity, track1.pt());
                hPID.fill(HIST("After/hTOFnsigKa_mult_pt"), track2.tofNSigmaKa(), multiplicity, track2.pt());
              }

            } else if (selectionConfig.isPDGCheckMC && (track1PDG == PDG_t::kKPlus)) {
              if (!selectionConfig.isapplypTdepPID && !selectionConfig.isapplypTdepPIDTOF && !(selectionPID(track1, 1) && selectionPID(track2, 0))) { // pion and kaon
                continue;
              } else if (selectionConfig.isapplypTdepPID && !selectionConfig.isapplypTdepPIDTOF && !(selectionPIDPtDep(track1, 1) && selectionPIDPtDep(track2, 0))) { // pion and kaon
                continue;
              } else if (!selectionConfig.isapplypTdepPID && selectionConfig.isapplypTdepPIDTOF && !(selectionPIDptDepTOF(track1, 1) && selectionPIDptDepTOF(track2, 0))) { // pion and kaon
                continue;
              }
              rEventSelection.fill(HIST("recMCparticles"), 12.5);

              if (selectionConfig.isApplyParticleMID) {
                if (selectionMID(track2, 1)) // Pion misidentified as kaon
                  continue;
                // if (selectionMID(track2, 2)) // Kaon misidentified as proton
                //   continue;
                if (selectionMID(track1, 0)) // Kaon misidentified as pion
                  continue;
                // if (selectionMID(track1, 2)) // Pion misidentified as proton
                //   continue;
              }

              if (selectionConfig.isApplyParticleMIDPtDep) {
                if (selectionMIDPtDep(track1, 0)) // Kaon misidentified as pion
                  continue;
                if (selectionMIDPtDep(track2, 1)) // Pion misidentified as kaon
                  continue;
              }

              if (std::abs(track1.rapidity(o2::track::PID::getMass(o2::track::PID::Kaon))) > selectionConfig.ctrackRapidity)
                continue;

              if (std::abs(track2.rapidity(o2::track::PID::getMass(o2::track::PID::Pion))) > selectionConfig.ctrackRapidity)
                continue;
              rEventSelection.fill(HIST("recMCparticles"), 14.5);

              if (cQAplots) {
                hPID.fill(HIST("After/hTPCnsigKa_mult_pt"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
                hPID.fill(HIST("After/hTPCnsigPi_mult_pt"), track2.tpcNSigmaPi(), multiplicity, track2.pt());
                hPID.fill(HIST("After/hTOFnsigKa_mult_pt"), track1.tofNSigmaKa(), multiplicity, track1.pt());
                hPID.fill(HIST("After/hTOFnsigPi_mult_pt"), track2.tofNSigmaPi(), multiplicity, track2.pt());
              }
            }

            /* if (selectionConfig.isApplyCutsOnMother) {
              if (mothertrack1.pt() >= selectionConfig.cMaxPtMotherCut) // excluding candidates in overflow
                continue;
              if ((std::sqrt(mothertrack1.e() * mothertrack1.e() - mothertrack1.p() * mothertrack1.p())) >= selectionConfig.cMaxMinvMotherCut) // excluding candidates in overflow
                continue;
            } */

            if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
              hInvMass.fill(HIST("h1KSRecsplit"), mothertrack1.pt());
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 15.5);

            if (!selectionPair(track1, track2)) {
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 16.5);

            oldindex = mothertrack1.globalIndex();
            if (track1PDG == PDG_t::kPiPlus) {
              daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massPi);
              daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
            } else if (track1PDG == PDG_t::kKPlus) {
              daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
              daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
            }
            mother = daughter1 + daughter2; // Kstar meson

            hInvMass.fill(HIST("h2KstarRecpt2"), mothertrack1.pt(), multiplicity, std::sqrt(mothertrack1.e() * mothertrack1.e() - mothertrack1.p() * mothertrack1.p()));
            hInvMass.fill(HIST("h2KstarRecptCalib2"), mothertrack1.pt(), multiplicityRec, std::sqrt(mothertrack1.e() * mothertrack1.e() - mothertrack1.p() * mothertrack1.p()));

            if (applyRecMotherRapidity && mother.Rapidity() >= selectionConfig.rapidityMotherData) {
              continue;
            }

            hInvMass.fill(HIST("h1KstarRecMass"), mother.M());
            hInvMass.fill(HIST("h2KstarRecpt1"), mother.Pt(), multiplicity, mother.M());
            hInvMass.fill(HIST("h2KstarRecptCalib1"), mother.Pt(), multiplicityRec, mother.M());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(Kstarqa, processRec, "Process Reconstructed", false);

  void processRec2(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&, EventMCGenerated const&)
  {

    if (!collision.has_mcCollision()) {
      return;
    }

    double multiplicityRec = -1.0;
    const auto& mcCollisionRec = collision.mcCollision_as<EventMCGenerated>();
    multiplicityRec = mcCollisionRec.centFT0M();

    if (selectionConfig.isINELgt0 && !collision.isInelGt0()) {
      return;
    }
    // multiplicity = collision.centFT0M();

    if (cSelectMultEstimator == kFT0M) {
      multiplicity = collision.centFT0M();
    } else if (cSelectMultEstimator == kFT0A) {
      multiplicity = collision.centFT0A();
    } else if (cSelectMultEstimator == kFT0C) {
      multiplicity = collision.centFT0C();
    } else if (cSelectMultEstimator == kFV0A) {
      multiplicity = collision.centFV0A();
    } else {
      multiplicity = collision.centFT0M(); // default
    }

    hInvMass.fill(HIST("hAllRecCollisions"), multiplicity);
    hInvMass.fill(HIST("hAllRecCollisionsCalib"), multiplicityRec);

    if (!selectionEvent(collision, true)) { // don't fill event cut histogram
      return;
    }

    hInvMass.fill(HIST("h1RecMult"), multiplicity);
    hInvMass.fill(HIST("h1RecMult2"), multiplicityRec);

    if (cQAevents) {
      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
    }

    auto oldindex = -999;
    for (const auto& track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }

      if (!track1.has_mcParticle()) {
        continue;
      }

      if (cQAevents) {
        rEventSelection.fill(HIST("hDcaxy_cent_pt"), track1.dcaXY(), multiplicity, track1.pt());
        rEventSelection.fill(HIST("hDcaz_cent_pt"), track1.dcaZ(), multiplicity, track1.pt());
      }

      auto track1ID = track1.index();
      for (const auto& track2 : tracks) {
        rEventSelection.fill(HIST("recMCparticles"), 0.5);
        if (!track2.has_mcParticle()) {
          continue;
        }
        rEventSelection.fill(HIST("recMCparticles"), 1.5);

        if (!selectionTrack(track2)) {
          continue;
        }
        rEventSelection.fill(HIST("recMCparticles"), 2.5);

        auto track2ID = track2.index();
        if (track2ID == track1ID) {
          continue;
        }
        rEventSelection.fill(HIST("recMCparticles"), 3.5);

        if (track1.sign() * track2.sign() >= 0) {
          continue;
        }
        rEventSelection.fill(HIST("recMCparticles"), 4.5);

        const auto mctrack1 = track1.mcParticle();
        const auto mctrack2 = track2.mcParticle();
        int track1PDG = std::abs(mctrack1.pdgCode());
        int track2PDG = std::abs(mctrack2.pdgCode());

        if (cQAplots && track1PDG == PDG_t::kKPlus && track2PDG == PDG_t::kPiPlus) {
          hPID.fill(HIST("Before/hNsigma_TPC_TOF_Ka_before"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
          hPID.fill(HIST("Before/hNsigma_TPC_TOF_Pi_before"), track2.tpcNSigmaPi(), track2.tofNSigmaPi(), track2.pt());
          hPID.fill(HIST("Before/hTPCnsigKa_mult_pt"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
          hPID.fill(HIST("Before/hTPCnsigPi_mult_pt"), track2.tpcNSigmaPi(), multiplicity, track2.pt());
          hPID.fill(HIST("Before/hTOFnsigKa_mult_pt"), track1.tofNSigmaKa(), multiplicity, track1.pt());
          hPID.fill(HIST("Before/hTOFnsigPi_mult_pt"), track2.tofNSigmaPi(), multiplicity, track2.pt());
        }

        if (!mctrack1.isPhysicalPrimary()) {
          continue;
        }

        if (!mctrack2.isPhysicalPrimary()) {
          continue;
        }
        rEventSelection.fill(HIST("recMCparticles"), 5.5);

        // if (!(track1PDG == PDG_t::kKPlus && track2PDG == PDG_t::kPiPlus)) {
        //   continue;
        // }
        if (selectionConfig.isPDGCheckMC && (track1PDG != PDG_t::kKPlus)) {
          continue;
        }
        if (selectionConfig.isPDGCheckMC && (track2PDG != PDG_t::kPiPlus)) {
          continue;
        }
        rEventSelection.fill(HIST("recMCparticles"), 6.5);

        for (const auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
          for (const auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
            if (selectionConfig.isPDGCheckMC && (mothertrack1.pdgCode() != mothertrack2.pdgCode())) {
              continue;
            }

            if (mothertrack1.globalIndex() != mothertrack2.globalIndex()) {
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 7.5);

            if (!mothertrack1.producedByGenerator()) {
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 8.5);

            if (std::abs(mothertrack1.y()) >= selectionConfig.rapidityMotherData) {
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 9.5);

            if (selectionConfig.isPDGCheckMC && (std::abs(mothertrack1.pdgCode()) != o2::constants::physics::kK0Star892)) {
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 10.5);

            if (!selectionConfig.isapplypTdepPID && !selectionConfig.isapplypTdepPIDTOF && !(selectionPID(track1, 0) && selectionPID(track2, 1))) { // pion and kaon
              continue;
            } else if (selectionConfig.isapplypTdepPID && !selectionConfig.isapplypTdepPIDTOF && !(selectionPIDPtDep(track1, 0) && selectionPIDPtDep(track2, 1))) { // pion and kaon
              continue;
            } else if (!selectionConfig.isapplypTdepPID && selectionConfig.isapplypTdepPIDTOF && !(selectionPIDptDepTOF(track1, 0) && selectionPIDptDepTOF(track2, 1))) { // pion and kaon
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 11.5);

            if (selectionConfig.isApplyParticleMID) {
              if (selectionMID(track1, 1)) // Pion misidentified as kaon
                continue;
              // if (selectionMID(track1, 2)) // Kaon misidentified as proton
              //   continue;
              if (selectionMID(track2, 0)) // Kaon misidentified as pion
                continue;
              // if (selectionMID(track2, 2)) // Pion misidentified as proton
              //   continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 12.5);

            if (selectionConfig.isApplyParticleMIDPtDep) {
              if (selectionMIDPtDep(track1, 1)) // Pion misidentified as kaon
                continue;
              if (selectionMIDPtDep(track2, 0)) // Kaon misidentified as pion
                continue;
            }

            if (std::abs(track1.rapidity(o2::track::PID::getMass(o2::track::PID::Kaon))) > selectionConfig.ctrackRapidity)
              continue;

            if (std::abs(track2.rapidity(o2::track::PID::getMass(o2::track::PID::Pion))) > selectionConfig.ctrackRapidity)
              continue;

            rEventSelection.fill(HIST("recMCparticles"), 13.5);

            /* if (selectionConfig.isApplyCutsOnMother) {
              if (mothertrack1.pt() >= selectionConfig.cMaxPtMotherCut) // excluding candidates in overflow
                continue;
              if ((std::sqrt(mothertrack1.e() * mothertrack1.e() - mothertrack1.p() * mothertrack1.p())) >= selectionConfig.cMaxMinvMotherCut) // excluding candidates in overflow
                continue;
            } */

            if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
              hInvMass.fill(HIST("h1KSRecsplit"), mothertrack1.pt());
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 14.5);

            if (!selectionPair(track1, track2)) {
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 15.5);

            oldindex = mothertrack1.globalIndex();
            if (track1.sign() * track2.sign() < 0) {
              daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
              daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
              mother = daughter1 + daughter2; // Kstar meson

              hInvMass.fill(HIST("h2KstarRecpt2"), mothertrack1.pt(), multiplicity, std::sqrt(mothertrack1.e() * mothertrack1.e() - mothertrack1.p() * mothertrack1.p()));
              hInvMass.fill(HIST("h2KstarRecptCalib2"), mothertrack1.pt(), multiplicityRec, std::sqrt(mothertrack1.e() * mothertrack1.e() - mothertrack1.p() * mothertrack1.p()));

              if (applyRecMotherRapidity && mother.Rapidity() >= selectionConfig.rapidityMotherData) {
                continue;
              }

              hInvMass.fill(HIST("h1KstarRecMass"), mother.M());
              hInvMass.fill(HIST("h2KstarRecpt1"), mother.Pt(), multiplicity, mother.M());
              hInvMass.fill(HIST("h2KstarRecptCalib1"), mother.Pt(), multiplicityRec, mother.M());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(Kstarqa, processRec2, "Process Reconstructed 2", false);

  void processSEPhi(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)

  {

    int occupancy = collision.trackOccupancyInTimeRange();
    rEventSelection.fill(HIST("hOccupancy"), occupancy);

    if (!selectionEvent(collision, true)) { // fill data event cut histogram
      return;
    }

    multiplicity = -1;

    if (cSelectMultEstimator == kFT0M) {
      multiplicity = collision.centFT0M();
    } else if (cSelectMultEstimator == kFT0A) {
      multiplicity = collision.centFT0A();
    } else if (cSelectMultEstimator == kFT0C) {
      multiplicity = collision.centFT0C();
    } else if (cSelectMultEstimator == kFV0A) {
      multiplicity = collision.centFV0A();
    } else {
      multiplicity = collision.centFT0M(); // default
    }

    // Fill the event counter
    if (cQAevents) {
      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
      rEventSelection.fill(HIST("hMultiplicity"), multiplicity);
      // rEventSelection.fill(HIST("multdist_FT0M"), collision.multFT0M());
    }

    for (const auto& [track1, track2] : combinations(CombinationsStrictlyUpperIndexPolicy(tracks, tracks))) {
      rEventSelection.fill(HIST("tracksCheckData"), 0.5);
      if (!selectionTrack(track1)) {
        continue;
      }
      if (!selectionTrack(track2)) {
        continue;
      }
      rEventSelection.fill(HIST("tracksCheckData"), 1.5);

      if (track1.globalIndex() == track2.globalIndex())
        continue;

      if (cQAplots) {
        // hPID.fill(HIST("Before/hNsigmaTPC_Ka_before"), track1.pt(), track1.tpcNSigmaKa());
        // hPID.fill(HIST("Before/hNsigmaTOF_Ka_before"), track1.pt(), track1.tofNSigmaKa());
        // hPID.fill(HIST("Before/hNsigmaTPC_Pi_before"), track2.pt(), track2.tpcNSigmaPi());
        // hPID.fill(HIST("Before/hNsigmaTOF_Pi_before"), track2.pt(), track2.tofNSigmaPi());
        hPID.fill(HIST("Before/hNsigma_TPC_TOF_Ka_before"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
        hPID.fill(HIST("Before/hNsigma_TPC_TOF_Pi_before"), track2.tpcNSigmaPi(), track2.tofNSigmaPi(), track2.pt());

        hPID.fill(HIST("Before/hTPCnsigKa_mult_pt"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
        // hPID.fill(HIST("Before/hTPCnsigPi_mult_pt"), track2.tpcNSigmaPi(), multiplicity, track2.pt());
        hPID.fill(HIST("Before/hTOFnsigKa_mult_pt"), track1.tofNSigmaKa(), multiplicity, track1.pt());
        hPID.fill(HIST("Before/hTOFnsigPi_mult_pt"), track2.tofNSigmaPi(), multiplicity, track2.pt());

        hOthers.fill(HIST("hCRFC_before"), track1.tpcCrossedRowsOverFindableCls());
        // hOthers.fill(HIST("dE_by_dx_TPC"), track1.p(), track1.tpcSignal());
        hOthers.fill(HIST("hphi"), track1.phi());

        // if (track1.sign() < 0) {
        //   hPID.fill(HIST("Before/h1PID_TPC_neg_kaon"), track1.tpcNSigmaKa());
        //   // hPID.fill(HIST("Before/h1PID_TPC_neg_pion"), track2.tpcNSigmaPi());
        //   hPID.fill(HIST("Before/h1PID_TOF_neg_kaon"), track1.tofNSigmaKa());
        //   // hPID.fill(HIST("Before/h1PID_TOF_neg_pion"), track2.tofNSigmaPi());
        // } else {
        //   hPID.fill(HIST("Before/h1PID_TPC_pos_kaon"), track1.tpcNSigmaKa());
        //   // hPID.fill(HIST("Before/h1PID_TPC_pos_pion"), track2.tpcNSigmaPi());
        //   hPID.fill(HIST("Before/h1PID_TOF_pos_kaon"), track1.tofNSigmaKa());
        //   // hPID.fill(HIST("Before/h1PID_TOF_pos_pion"), track2.tofNSigmaPi());
        // }
      }

      if (cQAevents) {
        rEventSelection.fill(HIST("hDcaxy_cent_pt"), track1.dcaXY(), multiplicity, track1.pt());
        rEventSelection.fill(HIST("hDcaz_cent_pt"), track1.dcaZ(), multiplicity, track1.pt());
      }

      // since we are using combinations full index policy, so repeated pairs are allowed, so we can check one with Kaon and other with kaon
      if (!selectionConfig.isapplypTdepPID && !selectionPID(track1, 1)) // Track 1 is checked with Kaon
        continue;
      if (!selectionConfig.isapplypTdepPID && !selectionPID(track2, 1)) // Track 2 is checked with kaon
        continue;

      if (selectionConfig.isapplypTdepPID && !selectionPIDPtDep(track1, 1)) // Track 1 is checked with Kaon
        continue;
      if (selectionConfig.isapplypTdepPID && !selectionPIDPtDep(track2, 1)) // Track 2 is checked with kaon
        continue;

      rEventSelection.fill(HIST("tracksCheckData"), 2.5);

      // if (cFakeTrack && isFakeTrack(track1, 1)) // Kaon
      //   continue;
      // if (cFakeTrack && isFakeTrack(track2, 1)) // Pion
      //   continue;
      rEventSelection.fill(HIST("tracksCheckData"), 3.5);

      if (std::abs(track1.rapidity(o2::track::PID::getMass(o2::track::PID::Kaon))) > selectionConfig.ctrackRapidity)
        continue;

      if (std::abs(track2.rapidity(o2::track::PID::getMass(o2::track::PID::Kaon))) > selectionConfig.ctrackRapidity)
        continue;
      rEventSelection.fill(HIST("tracksCheckData"), 4.5);

      /* if (selectionConfig.isApplyParticleMID) {
        if (selectionMID(track1, 0)) // Kaon misidentified as pion
          continue;
        if (selectionMID(track1, 2)) // Kaon misidentified as proton
          continue;
        if (selectionMID(track2, 0)) // Pion misidentified as pion
          continue;
        if (selectionMID(track2, 2)) // Pion misidentified as proton
          continue;
      } */

      rEventSelection.fill(HIST("tracksCheckData"), 5.5);

      if (cQAplots) {
        // hPID.fill(HIST("After/hDcaxyPi"), track2.dcaXY());
        hPID.fill(HIST("After/hDcaxyKa"), track1.dcaXY());
        // hPID.fill(HIST("After/hDcazPi"), track2.dcaZ());
        hPID.fill(HIST("After/hDcazKa"), track1.dcaZ());

        hPID.fill(HIST("After/hTPCnsigKa_mult_pt"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
        // hPID.fill(HIST("After/hTPCnsigPi_mult_pt"), track2.tpcNSigmaPi(), multiplicity, track2.pt());
        hPID.fill(HIST("After/hTOFnsigKa_mult_pt"), track1.tofNSigmaKa(), multiplicity, track1.pt());
        hPID.fill(HIST("After/hTOFnsigPi_mult_pt"), track2.tofNSigmaPi(), multiplicity, track2.pt());
        hOthers.fill(HIST("hEta_after"), track1.eta());
        hOthers.fill(HIST("hCRFC_after"), track1.tpcCrossedRowsOverFindableCls());
        // hPID.fill(HIST("After/hNsigmaKaonTPC_after"), track1.pt(), track1.tpcNSigmaKa());
        // hPID.fill(HIST("After/hNsigmaKaonTOF_after"), track1.pt(), track1.tofNSigmaKa());
        // hPID.fill(HIST("After/hNsigmaPionTPC_after"), track2.pt(), track2.tpcNSigmaPi());
        // hPID.fill(HIST("After/hNsigmaPionTOF_after"), track2.pt(), track2.tofNSigmaPi());
        hPID.fill(HIST("After/hNsigma_TPC_TOF_Ka_after"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
        hPID.fill(HIST("After/hNsigma_TPC_TOF_Pi_after"), track2.tpcNSigmaPi(), track2.tofNSigmaPi(), track2.pt());
      }

      if (!selectionPair(track1, track2)) {
        continue;
      }
      rEventSelection.fill(HIST("tracksCheckData"), 6.5);

      daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
      daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
      mother = daughter1 + daughter2; // Phi meson

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
      fillInvMass(daughter1, daughter2, mother, multiplicity, isMix, track1, track2);
    }
  }

  PROCESS_SWITCH(Kstarqa, processSEPhi, "Process Same event for Phi", false);

  void processMEPhi(EventCandidatesMix const&, TrackCandidates const&)
  {
    // Map estimator to pair and multiplicity accessor
    auto runMixing = [&](auto& pair, auto multiplicityGetter) {
      for (const auto& [c1, tracks1, c2, tracks2] : pair) {
        // if (!c1.sel8() || !c2.sel8())
        //   continue;

        if (!selectionEvent(c1, false) || !selectionEvent(c2, false)) { // don't fill event cut histogram
          continue;
        }

        multiplicity = multiplicityGetter(c1);

        for (const auto& [t1, t2] : o2::soa::combinations(
               o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
          if (!selectionTrack(t1) || !selectionTrack(t2))
            continue;
          if (!selectionPID(t1, 1) || !selectionPID(t2, 1))
            continue;

          if (std::abs(t1.rapidity(o2::track::PID::getMass(o2::track::PID::Kaon))) > selectionConfig.ctrackRapidity)
            continue;

          if (std::abs(t2.rapidity(o2::track::PID::getMass(o2::track::PID::Kaon))) > selectionConfig.ctrackRapidity)
            continue;

          // if (cFakeTrack && isFakeTrack(t1, 1)) // Kaon
          //   continue;
          // if (cFakeTrack && isFakeTrack(t2, 1)) // Kaon
          //   continue;

          /* if (selectionConfig.isApplyParticleMID) {
            if (selectionMID(t1, 0)) // Kaon misidentified as pion
              continue;
            if (selectionMID(t1, 2)) // Kaon misidentified as proton
              continue;
            if (selectionMID(t2, 0)) // Pion misidentified as pion
              continue;
            if (selectionMID(t2, 2)) // Pion misidentified as proton
              continue;
          } */

          if (!selectionPair(t1, t2)) {
            continue;
          }

          daughter1 = ROOT::Math::PxPyPzMVector(t1.px(), t1.py(), t1.pz(), massKa);
          daughter2 = ROOT::Math::PxPyPzMVector(t2.px(), t2.py(), t2.pz(), massKa);
          mother = daughter1 + daughter2;

          isMix = true;

          if (std::abs(mother.Rapidity()) < selectionConfig.rapidityMotherData) {
            fillInvMass(daughter1, daughter2, mother, multiplicity, isMix, t1, t2);
          }
        }
      }
    };

    // Call mixing based on selected estimator
    if (cSelectMultEstimator == kFT0M) {
      runMixing(pair1, [](const auto& c) { return c.centFT0M(); });
    } else if (cSelectMultEstimator == kFT0A) {
      runMixing(pair2, [](const auto& c) { return c.centFT0A(); });
    } else if (cSelectMultEstimator == kFT0C) {
      runMixing(pair3, [](const auto& c) { return c.centFT0C(); });
    } else if (cSelectMultEstimator == kFV0A) {
      runMixing(pair4, [](const auto& c) { return c.centFV0A(); });
    }
  }
  PROCESS_SWITCH(Kstarqa, processMEPhi, "Process Mixed event for Phi", false);

  void processGenPhi(EventMCGenerated::iterator const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<EventCandidatesMC>& collisions)
  // void processGenPhi(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<EventCandidatesMC>& collisions)
  {
    rEventSelection.fill(HIST("eventsCheckGen"), 0.5);

    int nChInel = 0;
    for (const auto& mcParticle : mcParticles) {
      auto pdgcode = std::abs(mcParticle.pdgCode());
      if (mcParticle.isPhysicalPrimary() && (pdgcode == PDG_t::kPiPlus || pdgcode == PDG_t::kKPlus || pdgcode == PDG_t::kProton || pdgcode == std::abs(PDG_t::kElectron) || pdgcode == std::abs(PDG_t::kMuonMinus))) {
        if (std::abs(mcParticle.eta()) < 1.0) {
          nChInel = nChInel + 1;
        }
      }
    }
    if (nChInel > 0 && std::abs(mcCollision.posZ()) < selectionConfig.cutzvertex)
      rEventSelection.fill(HIST("eventsCheckGen"), 1.5);

    std::vector<int64_t> selectedEvents(collisions.size());
    int nevts = 0;
    multiplicity = -1.0;
    // float impactParameter = mcCollision.impactParameter();

    bool isINELgt0true = false;

    if (pwglf::isINELgtNmc(mcParticles, 0, pdgDB)) {
      isINELgt0true = true;
    }
    if (selectionConfig.isINELgt0 && !isINELgt0true) {
      return;
    }

    // if (selectionConfig.isINELgt0 && !mcCollision.isInelGt0()) {
    //   return;
    // }
    rEventSelection.fill(HIST("eventsCheckGen"), 2.5);

    for (const auto& collision : collisions) {
      rEventSelection.fill(HIST("eventsCheckGen"), 3.5);

      if (!selectionEvent(collision, false)) { // don't fill event cut histogram
        continue;
      }
      multiplicity = collision.centFT0M();

      if (cSelectMultEstimator == kFT0M) {
        multiplicity = collision.centFT0M();
      } else if (cSelectMultEstimator == kFT0A) {
        multiplicity = collision.centFT0A();
      } else if (cSelectMultEstimator == kFT0C) {
        multiplicity = collision.centFT0C();
      } else if (cSelectMultEstimator == kFV0A) {
        multiplicity = collision.centFV0A();
      } else {
        multiplicity = collision.centFT0M(); // default
      }
      hInvMass.fill(HIST("h1GenMult"), multiplicity);

      int occupancy = collision.trackOccupancyInTimeRange();
      rEventSelection.fill(HIST("hOccupancy"), occupancy);

      selectedEvents[nevts++] = collision.mcCollision_as<EventMCGenerated>().globalIndex();
    }
    selectedEvents.resize(nevts);

    for (const auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.y()) < selectionConfig.rapidityMotherData && std::abs(mcParticle.pdgCode()) == o2::constants::physics::kPhi) {
        hInvMass.fill(HIST("hAllKstarGenCollisisons"), multiplicity, mcParticle.pt());
      }
    }

    const auto evtReconstructedAndSelected = std::find(selectedEvents.begin(), selectedEvents.end(), mcCollision.globalIndex()) != selectedEvents.end();
    hInvMass.fill(HIST("hAllGenCollisions"), multiplicity);
    if (!evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      return;
    }
    double genMultiplicity = mcCollision.centFT0M();
    hInvMass.fill(HIST("h1GenMult2"), genMultiplicity);
    hInvMass.fill(HIST("hAllGenCollisions1Rec"), multiplicity);
    rEventSelection.fill(HIST("eventsCheckGen"), 4.5);

    for (const auto& mcParticle : mcParticles) {

      if (std::abs(mcParticle.y()) >= selectionConfig.rapidityMotherData) {
        continue;
      }

      /* if (selectionConfig.isApplyCutsOnMother) {
        if (mcParticle.pt() >= selectionConfig.cMaxPtMotherCut) // excluding candidates in overflow
          continue;
        if ((std::sqrt(mcParticle.e() * mcParticle.e() - mcParticle.p() * mcParticle.p())) >= selectionConfig.cMaxMinvMotherCut) // excluding candidates in overflow
          continue;
      } */

      if (std::abs(mcParticle.pdgCode()) != o2::constants::physics::kPhi) {
        continue;
      }
      hInvMass.fill(HIST("hAllKstarGenCollisisons1Rec"), multiplicity, mcParticle.pt());

      auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
      if (kDaughters.size() != selectionConfig.noOfDaughters) {
        continue;
      }

      auto passPosKaon = false;
      auto passNegKaon = false;
      for (const auto& kCurrentDaughter : kDaughters) {
        if (!kCurrentDaughter.isPhysicalPrimary()) {
          continue;
        }

        if (kCurrentDaughter.pdgCode() == PDG_t::kKPlus) {
          passPosKaon = true;
          daughter1 = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);

        } else if (kCurrentDaughter.pdgCode() == PDG_t::kKMinus) {
          passNegKaon = true;
          daughter2 = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
        }
      }
      if (passPosKaon && passNegKaon) {
        mother = daughter1 + daughter2; // Phi meson
        hInvMass.fill(HIST("hk892GenpT"), mcParticle.pt(), multiplicity);
        hInvMass.fill(HIST("hk892GenpT2"), mother.Pt(), multiplicity);
        hInvMass.fill(HIST("hk892GenpTCalib1"), mcParticle.pt(), genMultiplicity);
        hInvMass.fill(HIST("hk892GenpTCalib2"), mother.Pt(), genMultiplicity);
        hInvMass.fill(HIST("h1genmass"), mother.M());
      }
    }
  }
  PROCESS_SWITCH(Kstarqa, processGenPhi, "Process Generated for Phi meson", false);

  void processEvSigLossFactorsPhi(McCollisionMults::iterator const& mcCollision, soa::SmallGroups<EventCandidatesMC> const& collisions, LabeledTracks const&, aod::McParticles const& mcParticles)
  {
    auto multiplicityNch = -1;
    multiplicityNch = mcCollision.multMCNParticlesEta05();
    hInvMass.fill(HIST("CorrFactors/hGenEvents"), multiplicityNch, 0.5);

    if (selectionConfig.checkVzEvSigLoss && std::abs(mcCollision.posZ()) > selectionConfig.cutzvertex)
      return;

    hInvMass.fill(HIST("CorrFactors/hGenEvents"), multiplicityNch, 1.5);

    if (selectionConfig.isINELgt0 && !mcCollision.isInelGt0()) {
      return;
    }
    hInvMass.fill(HIST("CorrFactors/hGenEvents"), multiplicityNch, 2.5);

    float multiplicity = -1.0;
    bool isSelectedEvent = false;

    for (auto const& collision : collisions) {
      if (!collision.has_mcCollision())
        continue;
      if (!selectionEvent(collision, false)) // don't fill event cut histogram
        continue;

      if (cSelectMultEstimator == kFT0M) {
        multiplicity = collision.centFT0M();
      } else if (cSelectMultEstimator == kFT0A) {
        multiplicity = collision.centFT0A();
      } else if (cSelectMultEstimator == kFT0C) {
        multiplicity = collision.centFT0C();
      } else if (cSelectMultEstimator == kFV0A) {
        multiplicity = collision.centFV0A();
      } else {
        multiplicity = collision.centFT0M(); // default
      }
      isSelectedEvent = true;
    }

    // auto multiplicityGen = -1;
    // multiplicityGen = mcCollision.centFT0M();

    hInvMass.fill(HIST("CorrFactors/hMultiplicityVsMultMC"), multiplicity, multiplicityNch);
    hInvMass.fill(HIST("CorrFactors/hNrecInGen"), collisions.size());
    hInvMass.fill(HIST("CorrFactors/MultiplicityGen"), multiplicityNch);
    if (isSelectedEvent) {
      hInvMass.fill(HIST("CorrFactors/MultiplicityRec"), multiplicityNch);
    }

    for (const auto& mcParticle : mcParticles) {

      if (std::abs(mcParticle.y()) >= selectionConfig.rapidityMotherData)
        continue;

      if (std::abs(mcParticle.pdgCode()) == o2::constants::physics::kPhi) {

        auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
        if (kDaughters.size() != selectionConfig.noOfDaughters) {
          continue;
        }

        bool hasPos = false;
        bool hasNeg = false;

        auto passkaon = false;
        auto passpion = false;
        for (const auto& kCurrentDaughter : kDaughters) {
          // if (!kCurrentDaughter.isPhysicalPrimary())
          //   continue;

          int pdgDau = kCurrentDaughter.pdgCode();
          int sign = (pdgDau > 0) - (pdgDau < 0);

          if (sign > 0)
            hasPos = true;
          if (sign < 0)
            hasNeg = true;

          if (kCurrentDaughter.pdgCode() == PDG_t::kKPlus) {
            passkaon = true;
            daughter1 = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);

          } else if (kCurrentDaughter.pdgCode() == PDG_t::kKMinus) {
            passpion = true;
            daughter2 = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
          }
        }

        if ((passkaon && passpion) && (hasPos && hasNeg)) {
          mother = daughter1 + daughter2; // Kstar meson

          hInvMass.fill(HIST("CorrFactors/h2dGenKstar"), multiplicity, mother.Pt());
          hInvMass.fill(HIST("CorrFactors/h3dGenKstarVsMultMCVsMultiplicity"), multiplicityNch, multiplicity, mother.Pt());
          hInvMass.fill(HIST("CorrFactors/hSignalLoss1"), mother.pt(), multiplicityNch);
          if (isSelectedEvent) {
            hInvMass.fill(HIST("CorrFactors/hSignalLoss2"), mother.pt(), multiplicityNch);
          }
        }
      }
    }

    if (collisions.size() == 0)
      return;

    hInvMass.fill(HIST("CorrFactors/hGenEvents"), multiplicityNch, 3.5);
  }
  PROCESS_SWITCH(Kstarqa, processEvSigLossFactorsPhi, "Process Event and Signal loss", false);

  void processRecPhi(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&, EventMCGenerated const&)
  {

    if (!collision.has_mcCollision()) {
      return;
    }

    double multiplicityRec = -1.0;
    // multiplicityRec = collision.mcCollision_as<EventMCGenerated>().centFT0M();
    const auto& mcCollisionRec = collision.mcCollision_as<EventMCGenerated>();
    multiplicityRec = mcCollisionRec.centFT0M();

    if (selectionConfig.isINELgt0 && !collision.isInelGt0()) {
      return;
    }
    // multiplicity = collision.centFT0M();

    if (cSelectMultEstimator == kFT0M) {
      multiplicity = collision.centFT0M();
    } else if (cSelectMultEstimator == kFT0A) {
      multiplicity = collision.centFT0A();
    } else if (cSelectMultEstimator == kFT0C) {
      multiplicity = collision.centFT0C();
    } else if (cSelectMultEstimator == kFV0A) {
      multiplicity = collision.centFV0A();
    } else {
      multiplicity = collision.centFT0M(); // default
    }

    hInvMass.fill(HIST("hAllRecCollisions"), multiplicity);
    hInvMass.fill(HIST("hAllRecCollisionsCalib"), multiplicityRec);

    if (!selectionEvent(collision, true)) { // fill MC event cut histogram
      return;
    }

    hInvMass.fill(HIST("h1RecMult"), multiplicity);
    hInvMass.fill(HIST("h1RecMult2"), multiplicityRec);

    if (cQAevents) {
      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
    }

    auto oldindex = -999;
    for (const auto& track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }

      if (!track1.has_mcParticle()) {
        continue;
      }

      if (cQAevents) {
        rEventSelection.fill(HIST("hDcaxy_cent_pt"), track1.dcaXY(), multiplicity, track1.pt());
        rEventSelection.fill(HIST("hDcaz_cent_pt"), track1.dcaZ(), multiplicity, track1.pt());
      }

      auto track1ID = track1.index();
      for (const auto& track2 : tracks) {
        rEventSelection.fill(HIST("recMCparticles"), 0.5);
        if (!track2.has_mcParticle()) {
          continue;
        }
        rEventSelection.fill(HIST("recMCparticles"), 1.5);

        if (!selectionTrack(track2)) {
          continue;
        }
        rEventSelection.fill(HIST("recMCparticles"), 2.5);

        auto track2ID = track2.index();
        if (track2ID <= track1ID) {
          continue;
        }
        rEventSelection.fill(HIST("recMCparticles"), 3.5);

        if (track1.sign() * track2.sign() >= 0) {
          continue;
        }
        rEventSelection.fill(HIST("recMCparticles"), 4.5);

        const auto mctrack1 = track1.mcParticle();
        const auto mctrack2 = track2.mcParticle();
        int track1PDG = std::abs(mctrack1.pdgCode());
        int track2PDG = std::abs(mctrack2.pdgCode());
        if (cQAplots) {
          hPID.fill(HIST("Before/hTPCnsigKa_mult_pt"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
          // hPID.fill(HIST("Before/hTPCnsigPi_mult_pt"), track2.tpcNSigmaPi(), multiplicity, track2.pt());
          hPID.fill(HIST("Before/hTOFnsigKa_mult_pt"), track1.tofNSigmaKa(), multiplicity, track1.pt());
          hPID.fill(HIST("Before/hTOFnsigPi_mult_pt"), track2.tofNSigmaPi(), multiplicity, track2.pt());
        }
        if (cQAplots && (track1PDG == PDG_t::kKPlus && track2PDG == PDG_t::kKMinus)) {
          hPID.fill(HIST("Before/hNsigma_TPC_TOF_Ka_before"), track1.tpcNSigmaKa(), track1.tofNSigmaKa(), track1.pt());
        }

        if (!mctrack1.isPhysicalPrimary()) {
          continue;
        }

        if (!mctrack2.isPhysicalPrimary()) {
          continue;
        }
        rEventSelection.fill(HIST("recMCparticles"), 5.5);

        if (selectionConfig.isPDGCheckMC && ((track1PDG != PDG_t::kKPlus) || (track2PDG != PDG_t::kKPlus))) {
          continue;
        }
        rEventSelection.fill(HIST("recMCparticles"), 6.5);
        rEventSelection.fill(HIST("recMCparticles"), 7.5);

        for (const auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
          for (const auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
            if (selectionConfig.isPDGCheckMC && (mothertrack1.pdgCode() != mothertrack2.pdgCode())) {
              continue;
            }

            if (mothertrack1.globalIndex() != mothertrack2.globalIndex()) {
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 8.5);

            if (!mothertrack1.producedByGenerator()) {
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 9.5);

            if (std::abs(mothertrack1.y()) >= selectionConfig.rapidityMotherData) {
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 10.5);

            if (selectionConfig.isPDGCheckMC && (std::abs(mothertrack1.pdgCode()) != o2::constants::physics::kPhi)) {
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 11.5);

            if (!selectionConfig.isapplypTdepPID && !(selectionPID(track1, 1) && selectionPID(track2, 1))) { // kaon and kaon
              continue;
            } else if (selectionConfig.isapplypTdepPID && !(selectionPIDPtDep(track1, 1) && selectionPIDPtDep(track2, 1))) { // kaon and kaon
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 12.5);

            /* if (selectionConfig.isApplyParticleMID) {
              if (selectionMID(track1, 0)) // Kaon misidentified as pion
                continue;
              if (selectionMID(track1, 2)) // Kaon misidentified as proton
                continue;
              if (selectionMID(track2, 0)) // Kaon misidentified as pion
                continue;
              if (selectionMID(track2, 2)) // Kaon misidentified as proton
                continue;
            } */
            rEventSelection.fill(HIST("recMCparticles"), 13.5);

            if (std::abs(track1.rapidity(o2::track::PID::getMass(o2::track::PID::Kaon))) > selectionConfig.ctrackRapidity)
              continue;

            if (std::abs(track2.rapidity(o2::track::PID::getMass(o2::track::PID::Kaon))) > selectionConfig.ctrackRapidity)
              continue;
            rEventSelection.fill(HIST("recMCparticles"), 14.5);
            if (cQAplots) {
              hPID.fill(HIST("After/hTPCnsigKa_mult_pt"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
              // hPID.fill(HIST("After/hTPCnsigPi_mult_pt"), track2.tpcNSigmaPi(), multiplicity, track2.pt());
              hPID.fill(HIST("After/hTOFnsigKa_mult_pt"), track1.tofNSigmaKa(), multiplicity, track1.pt());
              // hPID.fill(HIST("After/hTOFnsigPi_mult_pt"), track2.tofNSigmaPi(), multiplicity, track2.pt());
            }

            /* if (selectionConfig.isApplyCutsOnMother) {
              if (mothertrack1.pt() >= selectionConfig.cMaxPtMotherCut) // excluding candidates in overflow
                continue;
              if ((std::sqrt(mothertrack1.e() * mothertrack1.e() - mothertrack1.p() * mothertrack1.p())) >= selectionConfig.cMaxMinvMotherCut) // excluding candidates in overflow
                continue;
            } */

            if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
              hInvMass.fill(HIST("h1KSRecsplit"), mothertrack1.pt());
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 15.5);

            if (!selectionPair(track1, track2)) {
              continue;
            }
            rEventSelection.fill(HIST("recMCparticles"), 16.5);

            oldindex = mothertrack1.globalIndex();
            if (track1.sign() * track2.sign() < 0) {
              daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
              daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
              mother = daughter1 + daughter2; // Phi meson

              hInvMass.fill(HIST("h2KstarRecpt2"), mothertrack1.pt(), multiplicity, std::sqrt(mothertrack1.e() * mothertrack1.e() - mothertrack1.p() * mothertrack1.p()));
              hInvMass.fill(HIST("h2KstarRecptCalib2"), mothertrack1.pt(), multiplicityRec, std::sqrt(mothertrack1.e() * mothertrack1.e() - mothertrack1.p() * mothertrack1.p()));

              if (applyRecMotherRapidity && mother.Rapidity() >= selectionConfig.rapidityMotherData) {
                continue;
              }

              hInvMass.fill(HIST("h1KstarRecMass"), mother.M());
              hInvMass.fill(HIST("h2KstarRecpt1"), mother.Pt(), multiplicity, mother.M());
              hInvMass.fill(HIST("h2KstarRecptCalib1"), mother.Pt(), multiplicityRec, mother.M());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(Kstarqa, processRecPhi, "Process Reconstructed", false);

  void processEvtLossSigLossLightIonMC(aod::McCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<EventCandidatesMC>& recCollisions)
  {

    // Event loss estimation
    auto impactPar = mcCollision.impactParameter();
    hInvMass.fill(HIST("MCcorrections/hImpactParameterGen"), impactPar);

    bool isSelectedEvent = false;
    auto centrality = -999.;
    for (const auto& RecCollision : recCollisions) {
      if (!RecCollision.has_mcCollision())
        continue;
      if (!selectionEvent(RecCollision, false)) // don't fill event cut histogram
        continue;

      if (cSelectMultEstimator == kFT0M) {
        centrality = RecCollision.centFT0M();
      } else if (cSelectMultEstimator == kFT0A) {
        centrality = RecCollision.centFT0A();
      } else if (cSelectMultEstimator == kFT0C) {
        centrality = RecCollision.centFT0C();
      } else if (cSelectMultEstimator == kFV0A) {
        centrality = RecCollision.centFV0A();
      } else {
        centrality = RecCollision.centFT0M(); // default
      }

      isSelectedEvent = true;
    }

    if (isSelectedEvent) {
      hInvMass.fill(HIST("MCcorrections/hImpactParameterRec"), impactPar);
      hInvMass.fill(HIST("MCcorrections/hImpactParvsCentrRec"), centrality, impactPar);
    }

    // Generated MC
    for (const auto& mcPart : mcParticles) {

      if (std::abs(mcPart.y()) >= selectionConfig.rapidityMotherData || std::abs(mcPart.pdgCode()) != o2::constants::physics::kK0Star892)
        continue;

      // signal loss estimation
      hInvMass.fill(HIST("MCcorrections/hKstarGenBeforeEvtSel"), mcPart.pt(), impactPar);
      if (isSelectedEvent) {
        // signal loss estimation
        hInvMass.fill(HIST("MCcorrections/hKstarGenAfterEvtSel"), mcPart.pt(), impactPar);
      }
    } // end loop on gen particles
  }
  PROCESS_SWITCH(Kstarqa, processEvtLossSigLossLightIonMC, "Process Signal Loss, Event Loss for Kstar in Light Ion", false);

  void processEvtLossSigLossMCPhi(EventMCGenerated::iterator const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<EventCandidatesMC>& recCollisions)
  // void processEvtLossSigLossMCPhi(aod::McCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<EventCandidatesMC>& recCollisions)
  {
    // if (selectionConfig.isINELgt0 && !mcCollision.isInelGt0()) {
    //   return;
    // }

    bool isINELgt0true = false;

    if (pwglf::isINELgtNmc(mcParticles, 0, pdgDB)) {
      isINELgt0true = true;
    }
    if (selectionConfig.isINELgt0 && !isINELgt0true) {
      return;
    }

    if (selectionConfig.checkVzEvSigLoss && (std::abs(mcCollision.posZ()) > selectionConfig.cutzvertex)) {
      return;
    }

    auto impactPar = mcCollision.impactParameter();
    // auto multiplicityRec = -1;
    auto multiplicityGen = -1;
    multiplicityGen = mcCollision.centFT0M();
    hInvMass.fill(HIST("MCcorrections/hImpactParameterGen"), impactPar);
    hInvMass.fill(HIST("MCcorrections/MultiplicityGen"), multiplicityGen);

    bool isSelectedEvent = false;
    auto multiplicity1 = -999.;
    for (const auto& RecCollision : recCollisions) {
      if (!RecCollision.has_mcCollision())
        continue;
      if (!selectionEvent(RecCollision, false)) // don't fill event cut histogram
        continue;
      // multiplicity1 = RecCollision.centFT0M();
      // const auto& mcCollisionRec = RecCollision.mcCollision_as<EventMCGenerated>();
      // multiplicityRec = mcCollisionRec.centFT0M();

      if (cSelectMultEstimator == kFT0M) {
        multiplicity1 = RecCollision.centFT0M();
      } else if (cSelectMultEstimator == kFT0A) {
        multiplicity1 = RecCollision.centFT0A();
      } else if (cSelectMultEstimator == kFT0C) {
        multiplicity1 = RecCollision.centFT0C();
      } else if (cSelectMultEstimator == kFV0A) {
        multiplicity1 = RecCollision.centFV0A();
      } else {
        multiplicity1 = RecCollision.centFT0M(); // default
      }
      isSelectedEvent = true;
    }

    // Event loss
    if (isSelectedEvent) {
      hInvMass.fill(HIST("MCcorrections/hImpactParameterRec"), impactPar);
      hInvMass.fill(HIST("MCcorrections/MultiplicityRec"), multiplicityGen);
      hInvMass.fill(HIST("MCcorrections/hImpactParametervsMultiplicity"), impactPar, multiplicity1);
    }

    // Generated MC
    for (const auto& mcPart : mcParticles) {
      if (std::abs(mcPart.y()) >= selectionConfig.rapidityMotherData || std::abs(mcPart.pdgCode()) != o2::constants::physics::kPhi)
        continue;

      // signal loss estimation
      hInvMass.fill(HIST("MCcorrections/hSignalLossDenominator"), mcPart.pt(), multiplicityGen);
      if (isSelectedEvent) {
        hInvMass.fill(HIST("MCcorrections/hSignalLossNumerator"), mcPart.pt(), multiplicityGen);
      }
    } // end loop on gen particles
  }
  PROCESS_SWITCH(Kstarqa, processEvtLossSigLossMCPhi, "Process Signal Loss, Event Loss for Phi", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Kstarqa>(cfgc)};
}
