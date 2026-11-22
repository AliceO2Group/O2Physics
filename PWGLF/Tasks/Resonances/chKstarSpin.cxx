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
///
/// \file chKstarSpin.cxx
/// \brief Charged K*(892) spin alignment on AO2D
/// \author Prottay Das <prottay.das@cern.ch>, Subhadeep Mandal <subhadeep.mandal@cern.ch>
/// \since 22/11/2026

#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector3D.h> // IWYU pragma: keep
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h> // IWYU pragma: keep
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <THn.h>
#include <TPDGCode.h>
#include <TRandom3.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <cinttypes>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <iterator>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;

struct chKstarSpin {

  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  // ------------------------------------------------------------------
  // Event selection
  // ------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "eventSelection";
    Configurable<float> cfgCutVertexZ{"cfgCutVertexZ", 10.0f, "Accepted |z-vertex| range (cm)"};
    Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 0.0f, "Accepted minimum FT0C centrality"};
    Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 80.0f, "Accepted maximum FT0C centrality"};
    Configurable<int> cfgOccCut{"cfgOccCut", 2000, "Maximum track occupancy in time range"};
    Configurable<bool> useNoSameBunchPileup{"useNoSameBunchPileup", true, "Apply kNoSameBunchPileup"};
    Configurable<bool> useGoodITSLayersAll{"useGoodITSLayersAll", true, "Apply kIsGoodITSLayersAll"};
    Configurable<bool> useNoCollInTimeRangeStandard{"useNoCollInTimeRangeStandard", false, "Apply kNoCollInTimeRangeStandard"};
    Configurable<bool> useOccCut{"useOccCut", false, "Apply Occupancy Cut"};
    Configurable<bool> useIsGoodZvtxFT0vsPV{"useIsGoodZvtxFT0vsPV", true, "Apply kIsGoodZvtxFT0vsPV"};
    Configurable<bool> useIsVertexITSTPC{"useIsVertexITSTPC", false, "Apply kIsVertexITSTPC"};
    Configurable<bool> useIsVertexTOFmatched{"useIsVertexTOFmatched", false, "Apply kIsVertexTOFmatched"};
    Configurable<bool> useInelGt0{"useInelGt0", false, "Apply isInelGt0"};
    Configurable<bool> useHasFT0{"useHasFT0", false, "Apply has_foundFT0"};
    Configurable<bool> useHasBC{"useHasBC", false, "Apply has_foundBC"};
  } eventCuts;

  // ------------------------------------------------------------------
  // RCT event-quality selection
  // ------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "rctSelection";
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCut;

  // ------------------------------------------------------------------
  // Bachelor-track selection
  // ------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "trackSelection";
    Configurable<float> cfgCutPT{"cfgCutPT", 0.2f, "PT cut on daughter track"};
    Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta cut on daughter track"};
    Configurable<int> cfgITScluster{"cfgITScluster", 3, "Minimum number of ITS clusters"};
    Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Minimum number of TPC crossed rows"};
    Configurable<bool> useGlobalTrack{"useGlobalTrack", true, "Require isGlobalTrack()"};
  } trackCuts;

  // ------------------------------------------------------------------
  // Bachelor PID
  // ------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "bachelorPID";
    Configurable<float> cfgPIDSwitchPt{"cfgPIDSwitchPt", 0.5f, "pT at which TOF is required/used for bachelor PID (GeV/c)"};
    Configurable<bool> cfgRequireTOFAboveSwitch{"cfgRequireTOFAboveSwitch", true, "Require TOF above cfgPIDSwitchPt"};
    Configurable<bool> cfgUseExclusivePionPID{"cfgUseExclusivePionPID", true, "Reject bachelor tracks also compatible with K or p"};
    Configurable<float> nsigmaCutTPCPi{"nsigmaCutTPCPi", 3.0f, "TPC n-sigma cut for pion"};
    Configurable<float> nsigmaCutTOFPicomb{"nsigmaCutTOFPicomb", 3.0f, "Combined TPC+TOF n-sigma cut for pion"};
    Configurable<float> nsigmaCutTPCKa{"nsigmaCutTPCKa", 3.0f, "TPC n-sigma cut for kaon (exclusive PID only)"};
    Configurable<float> nsigmaCutTOFKacomb{"nsigmaCutTOFKacomb", 3.0f, "Combined TPC+TOF n-sigma cut for kaon (exclusive PID only)"};
    Configurable<float> nsigmaCutTPCPr{"nsigmaCutTPCPr", 3.0f, "TPC n-sigma cut for proton (exclusive PID only)"};
    Configurable<float> nsigmaCutTOFPrcomb{"nsigmaCutTOFPrcomb", 3.0f, "Combined TPC+TOF n-sigma cut for proton (exclusive PID only)"};
  } pidCuts;

  // ------------------------------------------------------------------
  // K0S selection (signal window)
  // ------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "k0sSelection";
    Configurable<float> confV0PtMin{"confV0PtMin", 0.0f, "Minimum K0S pT (GeV/c)"};
    Configurable<float> confV0PtMax{"confV0PtMax", 1000.0f, "Maximum K0S pT (GeV/c)"};
    Configurable<float> confV0Rap{"confV0Rap", 0.8f, "Maximum |y| of K0S"};
    Configurable<float> confV0DCADaughMax{"confV0DCADaughMax", 1.0f, "Maximum DCA between V0 daughters (cm)"};
    Configurable<double> confV0CPAMin{"confV0CPAMin", 0.9998, "Minimum K0S cosine pointing angle"};
    Configurable<float> confV0TranRadV0Min{"confV0TranRadV0Min", 1.5f, "Minimum K0S transverse radius (cm)"};
    Configurable<float> confV0TranRadV0Max{"confV0TranRadV0Max", 100.0f, "Maximum K0S transverse radius (cm)"};
    Configurable<double> cMaxV0DCA{"cMaxV0DCA", 1.2, "Maximum K0S DCA to PV (cm)"};
    Configurable<float> cMinV0DCAPi{"cMinV0DCAPi", 0.05f, "Minimum V0-daughter DCA to PV (cm)"};
    Configurable<float> cMaxV0LifeTime{"cMaxV0LifeTime", 50.0f, "Maximum K0S proper decay length c*tau (cm)"};
    Configurable<float> qtArmenterosMin{"qtArmenterosMin", 0.2f, "Minimum qT/|alpha| for K0S"};
    Configurable<float> confDaughEta{"confDaughEta", 0.8f, "Maximum |eta| of V0 daughters"};
    Configurable<float> cfgDaughPiPt{"cfgDaughPiPt", 0.2f, "Minimum pT of V0 daughters (GeV/c)"};
    Configurable<float> confDaughTPCnclsMin{"confDaughTPCnclsMin", 50.0f, "Minimum V0-daughter TPC clusters"};
    Configurable<float> confDaughTPCncrwsMin{"confDaughTPCncrwsMin", 70.0f, "Minimum V0-daughter TPC crossed rows"};
    Configurable<float> confDaughPIDCuts{"confDaughPIDCuts", 3.0f, "TPC pion PID cut for V0 daughters"};
    Configurable<float> cfgK0sMassMin{"cfgK0sMassMin", 0.45f, "Minimum K0S signal invariant mass (GeV/c^2)"};
    Configurable<float> cfgK0sMassMax{"cfgK0sMassMax", 0.55f, "Maximum K0S signal invariant mass (GeV/c^2)"};
  } k0sCuts;

  // ------------------------------------------------------------------
  // K* selection
  // ------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "kstarSelection";
    Configurable<float> cfgKstarRapidityMax{"cfgKstarRapidityMax", 0.5f, "Maximum |y| of K*(892)^+-"};
  } kstarCuts;

  // ------------------------------------------------------------------
  // Spin-alignment coordinate systems & rotational background
  // (Adapted from kstarpbpb.cxx)
  // ------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "spinAlignment";
    Configurable<int> cfgSAFrame{"cfgSAFrame", 0, "SA quantization axis: 0 = event plane (FT0C), 1 = production plane, 2 = random event plane"};
    Configurable<bool> usepolar{"usepolar", true, "true = fill cos(theta*), false = fill cos(2*(phi* - Psi))"};
    Configurable<int> cfgRndSeed{"cfgRndSeed", 0, "Seed for random event plane (0 = unique seed per job)"};
    Configurable<float> cfgEPNormalHarmonic{"cfgEPNormalHarmonic", 1.0f, "EP normal = (sin(n*Psi), -cos(n*Psi), 0); n = 1 is normal to Psi2 plane"};
    Configurable<bool> fillRotation{"fillRotation", true, "Fill rotational combinatorial background in SE"};
    Configurable<bool> cfgRotatePion{"cfgRotatePion", false, "false = rotate K0S (default), true = rotate bachelor pion"};
    Configurable<float> confMinRot{"confMinRot", 5.0f * o2::constants::math::PI / 6.0f, "Minimum rotation angle (rad)"};
    Configurable<float> confMaxRot{"confMaxRot", 7.0f * o2::constants::math::PI / 6.0f, "Maximum rotation angle (rad)"};
    Configurable<int> nBkgRotations{"nBkgRotations", 9, "Number of rotated copies per candidate"};
    Configurable<bool> additionalQAplots{"additionalQAplots", true, "Fill additional EP/subevent resolution and track QA"};
  } saConfig;

  // ------------------------------------------------------------------
  // MC efficiency & EP-frame settings (processMC)
  // ------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "mcEffEP";
    Configurable<bool> pdgcheck{"pdgcheck", true, "Check mother PDG == +-323 (K*+-)"};
    Configurable<bool> genacceptancecut{"genacceptancecut", true, "Apply pT and eta cuts on generated daughters"};
    Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", true, "Avoid split reconstructed candidates from same MC mother"};
    Configurable<int> psiSource{"psiSource", 1, "Psi in MC efficiency: 0 = fixed Psi = 0, 1 = Psi2 from gen charged primaries in [etaMin, etaMax], 2 = mcCollision.eventPlaneAngle()"};
    Configurable<float> etaMin{"etaMin", -3.3f, "eta min of generated-track Q2 vector (FT0C acceptance)"};
    Configurable<float> etaMax{"etaMax", -2.1f, "eta max of generated-track Q2 vector (FT0C acceptance)"};
    Configurable<bool> fillPhiVMinusPsi{"fillPhiVMinusPsi", false, "Fill gen/rec Acc x Eff maps vs pT, cos(theta*), and phi_V - Psi"};
    ConfigurableAxis axisPhiVMinusPsi{"axisPhiVMinusPsi", {6, 0.0, o2::constants::math::PI}, "#varphi_{V} - #Psi"};
  } mcEffEP;

  // ------------------------------------------------------------------
  // Event loss & signal loss in MC (processEvtLossSigLossMC)
  // ------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "evtSigLoss";
    Configurable<bool> cutVzGen{"cutVzGen", true, "Cut |vz| < cfgCutVertexZ on generated collision"};
    Configurable<bool> isApplyInelgt0{"isApplyInelgt0", false, "Require INEL > 0 for generated collision"};
    Configurable<bool> isApplyTVX{"isApplyTVX", false, "Require TVX (FT0A and FT0C MC mult > 0) for generated collision"};
    Configurable<bool> requireSingleReco{"requireSingleReco", true, "Require exactly one selected reco collision per MC collision"};
    ConfigurableAxis axisImpactPar{"axisImpactPar", {200, 0.0, 20.0}, "Impact parameter (fm)"};
    ConfigurableAxis axisMultMC{"axisMultMC", {500, 0.0, 2000.0}, "Generated charged particles in |eta| < 0.5"};
  } evtSigLoss;

  // ------------------------------------------------------------------
  // MC closure test settings (processMCClosure*, processMCClosure*RecoEP)
  // ------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "closure";
    Configurable<bool> doClosureTest{"doClosureTest", false, "Fill inline closure histograms in SE/NormalME/MC"};
    Configurable<bool> injectPolarization{"injectPolarization", false, "Reweight true K*+- with W(cos theta*) for rho00 = rho00Inj"};
    Configurable<float> rho00Inj{"rho00Inj", 1.0f / 3.0f, "Injected rho00 in [0,1] (1/3 = unpolarized)"};
    Configurable<bool> useTrueEP{"useTrueEP", true, "Use MC true event plane angle in processMCClosure (false: Psi = 0)"};
    Configurable<bool> splitSample{"splitSample", true, "Even MC events build Acc x Eff, odd MC events are test sample"};
    Configurable<bool> fillRotBkg{"fillRotBkg", false, "Fill rotated-pair background in MC closure"};
    Configurable<bool> genAcceptanceCut{"genAcceptanceCut", false, "Apply daughter pT/eta cut on generated K*+- in closure"};
    Configurable<bool> requirePrimaryTracks{"requirePrimaryTracks", false, "Require bachelor track to be MC physical primary in closure"};
    Configurable<bool> meRequireSingleReco{"meRequireSingleReco", true, "Require both ME collisions to come from single-reco MC collisions"};
    ConfigurableAxis axisCosThetaStar{"axisCosThetaStar", {10, -1.0, 1.0}, "cos(#vartheta*)"};
  } closure;

  // ------------------------------------------------------------------
  // Standard O2 dataframe event mixing
  // ------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "normalMixing";
    Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of standard mixed events per reference event"};
    Configurable<bool> cfgRequireSameRun{"cfgRequireSameRun", true, "Require both mixed events to belong to the same run"};
    ConfigurableAxis cfgMixVtxBins{"cfgMixVtxBins", {VARIABLE_WIDTH, -10.0f, -8.0f, -6.0f, -4.0f, -2.0f, 0.0f, 2.0f, 4.0f, 6.0f, 8.0f, 10.0f}, "Standard-ME z-vertex bins (cm)"};
    ConfigurableAxis cfgMixCentBins{"cfgMixCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f}, "Standard-ME FT0C centrality bins (%)"};
    ConfigurableAxis cfgMixPsiBins{"cfgMixPsiBins", {VARIABLE_WIDTH, -1.5707963f, -1.3089969f, -1.0471976f, -0.7853982f, -0.5235988f, -0.2617994f, 0.0f, 0.2617994f, 0.5235988f, 0.7853982f, 1.0471976f, 1.3089969f, 1.5707963f}, "Standard-ME FT0C Psi_2 bins (rad)"};
  } normalMixing;

  struct : ConfigurableGroup {
    std::string prefix = "normalMixingQAAxes";
    ConfigurableAxis cfgDeltaVz{"cfgDeltaVz", {VARIABLE_WIDTH, -4.0f, -3.0f, -2.0f, -1.5f, -1.0f, -0.5f, 0.0f, 0.5f, 1.0f, 1.5f, 2.0f, 3.0f, 4.0f}, "standard-ME Delta z_vtx"};
    ConfigurableAxis cfgDeltaCent{"cfgDeltaCent", {VARIABLE_WIDTH, -20.0f, -15.0f, -10.0f, -7.5f, -5.0f, -2.5f, 0.0f, 2.5f, 5.0f, 7.5f, 10.0f, 15.0f, 20.0f}, "standard-ME Delta centrality"};
    ConfigurableAxis cfgDeltaPsi{"cfgDeltaPsi", {VARIABLE_WIDTH, -1.5707963f, -1.3089969f, -1.0471976f, -0.7853982f, -0.5235988f, -0.2617994f, 0.0f, 0.2617994f, 0.5235988f, 0.7853982f, 1.0471976f, 1.3089969f, 1.5707963f}, "standard-ME wrapped Delta Psi_2"};
  } normalMixingQAAxes;

  // ------------------------------------------------------------------
  // Data-driven daughter-replacement mixing
  // ------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "dataDrivenMixing";
    Configurable<bool> cfgDoK0sReplacement{"cfgDoK0sReplacement", true, "Enable K0S(B)+pi(A) replacement branch"};
    Configurable<bool> cfgDoPionReplacement{"cfgDoPionReplacement", true, "Enable K0S(A)+pi(B) replacement branch"};
    Configurable<int> cfgMaxDistinctEventsPerLeg{"cfgMaxDistinctEventsPerLeg", 20, "Maximum distinct replacement events used per SE target and replacement leg"};
    Configurable<int> cfgMaxMatchesPerLeg{"cfgMaxMatchesPerLeg", 20, "Maximum accepted replacements per SE target and replacement leg"};
    Configurable<int> cfgMaxMatchesPerEventPerLeg{"cfgMaxMatchesPerEventPerLeg", 1, "Maximum accepted replacements from one event for one target/leg"};
    Configurable<int> cfgMaxK0sCandidatesScannedPerCell{"cfgMaxK0sCandidatesScannedPerCell", 512, "Maximum newest K0S pool candidates inspected per kinematic cell; <=0 scans full cell"};
    Configurable<int> cfgMaxPionCandidatesScannedPerCell{"cfgMaxPionCandidatesScannedPerCell", 128, "Maximum newest pion pool candidates inspected per kinematic cell; <=0 scans full cell"};
    Configurable<int> cfgK0sPoolDepthPerBin{"cfgK0sPoolDepthPerBin", 512, "Maximum stored K0S candidates per DD pool cell"};
    Configurable<int> cfgPionPoolDepthPerBin{"cfgPionPoolDepthPerBin", 512, "Maximum stored bachelor-pion candidates per DD pool cell"};
    Configurable<int> cfgMaxTotalK0sPoolCandidates{"cfgMaxTotalK0sPoolCandidates", 1000000, "Global K0S rolling-pool cap; <=0 disables global cap"};
    Configurable<int> cfgMaxTotalPionPoolCandidates{"cfgMaxTotalPionPoolCandidates", 1000000, "Global pion rolling-pool cap; <=0 disables global cap"};
    Configurable<int> cfgMaxK0sPerEventPerBin{"cfgMaxK0sPerEventPerBin", 8, "Maximum K0S candidates from one event inserted into one DD pool cell"};
    Configurable<int> cfgMaxPionsPerEventPerBin{"cfgMaxPionsPerEventPerBin", 12, "Maximum bachelor pions from one event inserted into one DD pool cell"};
    Configurable<int> cfgPoolMaxDataFrameAge{"cfgPoolMaxDataFrameAge", 120, "Maximum age of rolling-pool entries in processed dataframes; <=0 disables age pruning"};
    Configurable<int> cfgPoolPruneEveryNDataFrames{"cfgPoolPruneEveryNDataFrames", 5, "Prune stale DD pool cells every N processed dataframes"};
    Configurable<bool> cfgLogPool{"cfgLogPool", false, "Log rolling DD-pool size once per dataframe"};
    Configurable<uint64_t> cfgMixSeed{"cfgMixSeed", 0xdecafbadULL, "Deterministic seed for per-event/per-target downsampling"};
    Configurable<bool> cfgFillMatchedSESparses{"cfgFillMatchedSESparses", true, "Fill branch-specific SE target sparses for exact matched-population closure"};
    Configurable<float> cfgCandidateMassMin{"cfgCandidateMassMin", 0.60f, "Minimum mixed K0S-pi mass kept (GeV/c^2)"};
    Configurable<float> cfgCandidateMassMax{"cfgCandidateMassMax", 1.00f, "Maximum mixed K0S-pi mass kept (GeV/c^2)"};
    Configurable<float> cfgK0sMaxDeltaPt{"cfgK0sMaxDeltaPt", 0.10f, "Maximum |pT(K0S_B)-pT(K0S_A)| (GeV/c)"};
    Configurable<float> cfgK0sMaxDeltaEta{"cfgK0sMaxDeltaEta", 0.05f, "Maximum |eta(K0S_B)-eta(K0S_A)|"};
    Configurable<float> cfgK0sMaxDeltaPhi{"cfgK0sMaxDeltaPhi", 0.05f, "Maximum wrapped |phi(K0S_B)-phi(K0S_A)| (rad)"};
    Configurable<float> cfgPionMaxDeltaPt{"cfgPionMaxDeltaPt", 0.10f, "Maximum |pT(pi_B)-pT(pi_A)| (GeV/c)"};
    Configurable<float> cfgPionMaxDeltaEta{"cfgPionMaxDeltaEta", 0.05f, "Maximum |eta(pi_B)-eta(pi_A)|"};
    Configurable<float> cfgPionMaxDeltaPhi{"cfgPionMaxDeltaPhi", 0.05f, "Maximum wrapped |phi(pi_B)-phi(pi_A)| (rad)"};
    Configurable<int> cfgNeighborK0sPt{"cfgNeighborK0sPt", 1, "Neighbor K0S-pT lookup bins searched on each side"};
    Configurable<int> cfgNeighborK0sEta{"cfgNeighborK0sEta", 1, "Neighbor K0S-eta lookup bins searched on each side"};
    Configurable<int> cfgNeighborK0sPhi{"cfgNeighborK0sPhi", 1, "Neighbor K0S-phi lookup bins searched on each side (periodic)"};
    Configurable<int> cfgNeighborPionPt{"cfgNeighborPionPt", 1, "Neighbor pion-pT lookup bins searched on each side"};
    Configurable<int> cfgNeighborPionEta{"cfgNeighborPionEta", 1, "Neighbor pion-eta lookup bins searched on each side"};
    Configurable<int> cfgNeighborPionPhi{"cfgNeighborPionPhi", 1, "Neighbor pion-phi lookup bins searched on each side (periodic)"};

    ConfigurableAxis cfgEventVtxBins{"cfgEventVtxBins", {VARIABLE_WIDTH, -10.0f, -8.0f, -6.0f, -4.0f, -2.0f, 0.0f, 2.0f, 4.0f, 6.0f, 8.0f, 10.0f}, "DD event-pool z-vertex bins (cm)"};
    ConfigurableAxis cfgEventCentBins{"cfgEventCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f}, "DD event-pool FT0C centrality bins (%)"};
    ConfigurableAxis cfgK0sPtBins{"cfgK0sPtBins", {VARIABLE_WIDTH, 0.0f, 0.25f, 0.5f, 0.75f, 1.0f, 1.25f, 1.5f, 1.75f, 2.0f, 2.25f, 2.5f, 3.0f, 3.5f, 4.0f, 5.0f, 6.0f, 8.0f, 10.0f}, "K0S pT lookup bins (GeV/c)"};
    ConfigurableAxis cfgK0sEtaBins{"cfgK0sEtaBins", {VARIABLE_WIDTH, -0.8f, -0.7f, -0.6f, -0.5f, -0.4f, -0.3f, -0.2f, -0.1f, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f}, "K0S eta lookup bins"};
    ConfigurableAxis cfgK0sPhiBins{"cfgK0sPhiBins", {VARIABLE_WIDTH, 0.0f, 0.2617994f, 0.5235988f, 0.7853982f, 1.047198f, 1.308997f, 1.570796f, 1.832596f, 2.094395f, 2.356194f, 2.617994f, 2.879793f, 3.141593f, 3.403392f, 3.665191f, 3.926991f, 4.18879f, 4.45059f, 4.712389f, 4.974188f, 5.235988f, 5.497787f, 5.759587f, 6.021386f, 6.283185f}, "K0S phi lookup bins (rad)"};
    ConfigurableAxis cfgPionPtBins{"cfgPionPtBins", {VARIABLE_WIDTH, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.75f, 1.0f, 1.25f, 1.5f, 1.75f, 2.0f, 2.5f, 3.0f, 4.0f, 5.0f, 6.0f, 8.0f, 10.0f}, "Bachelor-pion pT lookup bins (GeV/c)"};
    ConfigurableAxis cfgPionEtaBins{"cfgPionEtaBins", {VARIABLE_WIDTH, -0.8f, -0.7f, -0.6f, -0.5f, -0.4f, -0.3f, -0.2f, -0.1f, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f}, "Bachelor-pion eta lookup bins"};
    ConfigurableAxis cfgPionPhiBins{"cfgPionPhiBins", {VARIABLE_WIDTH, 0.0f, 0.2617994f, 0.5235988f, 0.7853982f, 1.047198f, 1.308997f, 1.570796f, 1.832596f, 2.094395f, 2.356194f, 2.617994f, 2.879793f, 3.141593f, 3.403392f, 3.665191f, 3.926991f, 4.18879f, 4.45059f, 4.712389f, 4.974188f, 5.235988f, 5.497787f, 5.759587f, 6.021386f, 6.283185f}, "Bachelor-pion phi lookup bins (rad)"};
  } dataDrivenMixing;

  struct : ConfigurableGroup {
    std::string prefix = "dataDrivenMixingQAAxes";
    ConfigurableAxis cfgDeltaPt{"cfgDeltaPt", {VARIABLE_WIDTH, -0.20f, -0.15f, -0.10f, -0.075f, -0.05f, -0.025f, 0.0f, 0.025f, 0.05f, 0.075f, 0.10f, 0.15f, 0.20f}, "mixed-target Delta pT"};
    ConfigurableAxis cfgDeltaEta{"cfgDeltaEta", {VARIABLE_WIDTH, -0.10f, -0.075f, -0.05f, -0.04f, -0.03f, -0.02f, -0.01f, 0.0f, 0.01f, 0.02f, 0.03f, 0.04f, 0.05f, 0.075f, 0.10f}, "mixed-target Delta eta"};
    ConfigurableAxis cfgDeltaPhi{"cfgDeltaPhi", {VARIABLE_WIDTH, -0.20f, -0.15f, -0.10f, -0.075f, -0.05f, -0.025f, 0.0f, 0.025f, 0.05f, 0.075f, 0.10f, 0.15f, 0.20f}, "wrapped mixed-target Delta phi"};
    ConfigurableAxis cfgDeltaVz{"cfgDeltaVz", {VARIABLE_WIDTH, -4.0f, -3.0f, -2.0f, -1.5f, -1.0f, -0.5f, 0.0f, 0.5f, 1.0f, 1.5f, 2.0f, 3.0f, 4.0f}, "replacement-target Delta z_vtx"};
    ConfigurableAxis cfgDeltaCent{"cfgDeltaCent", {VARIABLE_WIDTH, -20.0f, -15.0f, -10.0f, -7.5f, -5.0f, -2.5f, 0.0f, 2.5f, 5.0f, 7.5f, 10.0f, 15.0f, 20.0f}, "replacement-target Delta centrality"};
    ConfigurableAxis cfgDeltaPsi{"cfgDeltaPsi", {VARIABLE_WIDTH, -1.5707963f, -1.3089969f, -1.0471976f, -0.7853982f, -0.5235988f, -0.2617994f, 0.0f, 0.2617994f, 0.5235988f, 0.7853982f, 1.0471976f, 1.3089969f, 1.5707963f}, "replacement-target wrapped Delta Psi2"};
    ConfigurableAxis cfgNMatches{"cfgNMatches", {VARIABLE_WIDTH, -0.5f, 0.5f, 1.5f, 2.5f, 3.5f, 4.5f, 5.5f, 6.5f, 7.5f, 8.5f, 9.5f, 10.5f, 11.5f, 12.5f, 13.5f, 14.5f, 15.5f, 16.5f, 17.5f, 18.5f, 19.5f, 20.5f}, "accepted replacements per target"};
    ConfigurableAxis cfgPoolCandidates{"cfgPoolCandidates", {VARIABLE_WIDTH, 0.0f, 1000.0f, 5000.0f, 10000.0f, 25000.0f, 50000.0f, 100000.0f, 200000.0f, 500000.0f, 750000.0f, 1000000.0f, 1250000.0f}, "Candidates currently stored in a rolling DD pool"};
    ConfigurableAxis cfgMatchEffPt{"cfgMatchEffPt", {VARIABLE_WIDTH, 0.0f, 0.5f, 1.0f, 1.5f, 2.0f, 2.5f, 3.0f, 4.0f, 5.0f, 6.0f, 8.0f, 10.0f}, "SE K0S-pion pT for DD matching-efficiency QA"};
    ConfigurableAxis cfgMatchEffCent{"cfgMatchEffCent", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f}, "centrality for DD matching-efficiency QA"};
    ConfigurableAxis cfgMatchEffMass{"cfgMatchEffMass", {VARIABLE_WIDTH, 0.60f, 0.65f, 0.70f, 0.75f, 0.80f, 0.82f, 0.84f, 0.86f, 0.88f, 0.90f, 0.92f, 0.94f, 0.96f, 0.98f, 1.00f}, "SE K0S-pion mass for DD matching-efficiency QA"};
  } dataDrivenMixingQAAxes;

  // Main analysis sparse axes
  struct : ConfigurableGroup {
    std::string prefix = "sparseAxes";
    ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {VARIABLE_WIDTH, 0.6f, 0.8f, 0.9f, 1.0f}, "#it{M}_{K^{0}_{S}#pi} (GeV/#it{c}^{2})"};
    ConfigurableAxis configThnAxisPt{"configThnAxisPt", {VARIABLE_WIDTH, 0.0f, 10.0f}, "#it{p}_{T}^{K^{0}_{S}#pi} (GeV/#it{c})"};
    ConfigurableAxis configThnAxisSA{"configThnAxisSA", {VARIABLE_WIDTH, -1.0f, 0.0f, 1.0f}, "cos#it{#theta}^{*}"};
    ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {VARIABLE_WIDTH, 0.0f, 80.0f}, "Centrality (%)"};
    ConfigurableAxis configThnAxisCharge{"configThnAxisCharge", {VARIABLE_WIDTH, -1.5f, 0.0f, 1.5f}, "charge"};
  } sparseAxes;

  // Event QA axes
  struct : ConfigurableGroup {
    std::string prefix = "eventQAAxes";
    ConfigurableAxis cfgAxisEvtSel{"cfgAxisEvtSel", {VARIABLE_WIDTH, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f}, "Event-selection step"};
    ConfigurableAxis cfgAxisCentQA{"cfgAxisCentQA", {80, 0.0f, 80.0f}, "FT0C centrality (%)"};
    ConfigurableAxis cfgAxisVzQA{"cfgAxisVzQA", {40, -10.0f, 10.0f}, "z_{vtx} (cm)"};
    ConfigurableAxis cfgAxisPsiEventQA{"cfgAxisPsiEventQA", {72, -1.6f, 1.6f}, "#Psi_{2}^{FT0C} (rad)"};
  } eventQAAxes;

  // Candidate QA axes
  struct : ConfigurableGroup {
    std::string prefix = "candidateQAAxes";
    ConfigurableAxis cfgAxisCosThetaStarQA{"cfgAxisCosThetaStarQA", {40, -1.0f, 1.0f}, "cos#theta^{*}"};
  } candidateQAAxes;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  RCTFlagsChecker rctChecker;
  TRandom3 rndGen;

  enum SAFrame : int { kEventPlane = 0,
                       kProductionPlane = 1,
                       kRandomEventPlane = 2 };

  enum class PIDSpecies : uint8_t { Pion,
                                    Kaon,
                                    Proton };

  struct StoredK0sCandidate {
    uint64_t eventKey = 0;
    uint64_t dataFrameId = 0;
    int runNumber = -1;
    float psi2 = 0.0f;
    float posZ = 0.0f;
    float centrality = 0.0f;
    float pt = 0.0f;
    float eta = 0.0f;
    float phi = 0.0f;
    float px = 0.0f;
    float py = 0.0f;
    float pz = 0.0f;
  };

  struct StoredPionCandidate {
    uint64_t eventKey = 0;
    uint64_t dataFrameId = 0;
    int runNumber = -1;
    float psi2 = 0.0f;
    float posZ = 0.0f;
    float centrality = 0.0f;
    int8_t charge = 0;
    int64_t trackId = -1;
    float pt = 0.0f;
    float eta = 0.0f;
    float phi = 0.0f;
    float px = 0.0f;
    float py = 0.0f;
    float pz = 0.0f;
  };

  std::unordered_map<uint64_t, std::deque<StoredK0sCandidate>> ddRollingK0sPool;
  std::unordered_map<uint64_t, std::deque<StoredPionCandidate>> ddRollingPionPool;
  uint64_t ddDataFrameCounter = 0;
  size_t ddK0sPoolCandidateCount = 0;
  size_t ddPionPoolCandidateCount = 0;

  // ------------------------------------------------------------------
  // Data and MC Table Type Definitions, Filters, and Preslices
  // ------------------------------------------------------------------
  Filter collisionFilter = nabs(aod::collision::posZ) < eventCuts.cfgCutVertexZ;
  Filter centralityFilter = (aod::cent::centFT0C >= eventCuts.cfgCutCentralityMin) && (aod::cent::centFT0C < eventCuts.cfgCutCentralityMax);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::EPCalibrationTables, aod::FT0Mults, aod::Mults>>;
  using FullTrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr>;
  using ResoV0s = aod::V0Datas;

  // MC tables
  using CollisionMCTrueTable = aod::McCollisions;
  using McCollisionMults = soa::Join<aod::McCollisions, aod::MultMCExtras>;
  using TrackMCTrueTable = aod::McParticles;
  using CollisionMCRecTableCentFT0C = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs, aod::EvSels, aod::FT0Mults, aod::Mults>>;
  using CollisionMCRecTableClosure = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs, aod::EvSels, aod::EPCalibrationTables, aod::FT0Mults, aod::Mults>>;
  using CollisionMCRecMixTable = soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::CentFT0Cs, aod::EvSels, aod::FT0Mults, aod::Mults>>;
  using CollisionMCRecMixTableRecoEP = soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::CentFT0Cs, aod::EvSels, aod::EPCalibrationTables, aod::FT0Mults, aod::Mults>>;
  using TrackMCRecTable = soa::Join<FullTrackCandidates, aod::McTrackLabels>;

  Preslice<FullTrackCandidates> tracksPerCollisionDD = aod::track::collisionId;
  Preslice<TrackMCRecTable> mcTracksPerCollision = aod::track::collisionId;
  Preslice<ResoV0s> v0sPerCollisionDD = aod::v0data::collisionId;

  SliceCache normalMixingCache;
  using NormalMixingBinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, aod::epcalibrationtable::PsiFT0C>;
  NormalMixingBinningType normalMixingBinning{{normalMixing.cfgMixVtxBins, normalMixing.cfgMixCentBins, normalMixing.cfgMixPsiBins}, true};

  Pair<EventCandidates, FullTrackCandidates, ResoV0s, NormalMixingBinningType> normalMixedEventPairs{normalMixingBinning, normalMixing.cfgNoMixedEvents, -1, &normalMixingCache};

  using BinningTypeMCClosure = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  BinningTypeMCClosure closureMixingBinning{{normalMixing.cfgMixVtxBins, normalMixing.cfgMixCentBins}, true};
  Pair<CollisionMCRecMixTable, TrackMCRecTable, ResoV0s, BinningTypeMCClosure> mcClosureMixedEventPairs{closureMixingBinning, normalMixing.cfgNoMixedEvents, -1, &normalMixingCache};

  using BinningTypeMCClosureRecoEP = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, aod::epcalibrationtable::PsiFT0C>;
  BinningTypeMCClosureRecoEP closureMixingBinningRecoEP{{normalMixing.cfgMixVtxBins, normalMixing.cfgMixCentBins, normalMixing.cfgMixPsiBins}, true};
  Pair<CollisionMCRecMixTableRecoEP, TrackMCRecTable, ResoV0s, BinningTypeMCClosureRecoEP> mcClosureRecoEPMixedEventPairs{closureMixingBinningRecoEP, normalMixing.cfgNoMixedEvents, -1, &normalMixingCache};

  // ------------------------------------------------------------------
  // Utility & Coordinate-System Helper Functions
  // ------------------------------------------------------------------
  static float wrapPhiDD(float phi)
  {
    while (phi < 0.0f) {
      phi += o2::constants::math::TwoPI;
    }
    while (phi >= o2::constants::math::TwoPI) {
      phi -= o2::constants::math::TwoPI;
    }
    return phi;
  }

  static double getPhiInRange(double phi)
  {
    double result = phi;
    while (result < 0.0) {
      result += o2::constants::math::PI;
    }
    while (result >= o2::constants::math::PI) {
      result -= o2::constants::math::PI;
    }
    return result;
  }

  static float wrappedDeltaPhiDD(float phi1, float phi2)
  {
    float d = phi1 - phi2;
    while (d > o2::constants::math::PI) {
      d -= o2::constants::math::TwoPI;
    }
    while (d <= -o2::constants::math::PI) {
      d += o2::constants::math::TwoPI;
    }
    return d;
  }

  static float wrappedDeltaPsi2DD(float psi1, float psi2)
  {
    float d = psi1 - psi2;
    while (d > 0.5f * o2::constants::math::PI) {
      d -= o2::constants::math::PI;
    }
    while (d <= -0.5f * o2::constants::math::PI) {
      d += o2::constants::math::PI;
    }
    return d;
  }

  // Event-level reference angle for Spin Alignment: FT0C Psi2 or Random EP angle
  float getSAEventAngle(float psiEP)
  {
    if (saConfig.cfgSAFrame.value == kRandomEventPlane) {
      return static_cast<float>(rndGen.Uniform(-o2::constants::math::PIHalf, o2::constants::math::PIHalf));
    }
    return psiEP;
  }

  // In-plane angle of the quantization plane (used when usepolar == false)
  double getSAPlaneAngle(const ROOT::Math::PxPyPzMVector& mother, double psiSA) const
  {
    if (saConfig.cfgSAFrame.value == kProductionPlane) {
      return mother.Phi();
    }
    return psiSA;
  }

  // Quantization axis (normal to Event Plane, Production Plane, or Random Event Plane)
  ROOT::Math::XYZVector getSAAxis(const ROOT::Math::PxPyPzMVector& mother, double psiSA) const
  {
    if (saConfig.cfgSAFrame.value == kProductionPlane) {
      return ROOT::Math::XYZVector(0.0, 0.0, 1.0).Cross(mother.Vect()).Unit(); // z x p_K*
    }
    const auto n = static_cast<double>(saConfig.cfgEPNormalHarmonic.value);
    return {std::sin(n * psiSA), -std::cos(n * psiSA), 0.0};
  }

  // Computes cos(theta*) (when usepolar=true) or cos(2*(phi* - Psi_plane)) (when usepolar=false)
  // for the K0S daughter in the K*+- rest frame
  bool computeSA(const ROOT::Math::PxPyPzMVector& k0s, const ROOT::Math::PxPyPzMVector& kstar, float psiSA, float& saOut, float* rawCosThetaStarOut = nullptr) const
  {
    ROOT::Math::Boost boost{kstar.BoostToCM()};
    const ROOT::Math::PxPyPzMVector k0sCM = boost(k0s);
    const ROOT::Math::XYZVector k0sCMVec = k0sCM.Vect();
    const ROOT::Math::XYZVector axis = getSAAxis(kstar, psiSA);

    const double denom = std::sqrt(k0sCMVec.Mag2()) * std::sqrt(axis.Mag2());
    if (denom <= 0.0 || !std::isfinite(denom)) {
      return false;
    }

    auto cosThetaStar = static_cast<float>(axis.Dot(k0sCMVec) / denom);
    if (!std::isfinite(cosThetaStar)) {
      return false;
    }
    cosThetaStar = std::clamp(cosThetaStar, -1.0f, 1.0f);
    if (rawCosThetaStarOut != nullptr) {
      *rawCosThetaStarOut = cosThetaStar;
    }

    if (saConfig.usepolar.value) {
      saOut = cosThetaStar;
    } else {
      const double phiStarMinusPsi = getPhiInRange(k0sCM.Phi() - getSAPlaneAngle(kstar, psiSA));
      saOut = static_cast<float>(std::cos(2.0 * phiStarMinusPsi));
    }
    return std::isfinite(saOut);
  }

  // Forward-rapidity generated Q2 event-plane angle (psiSource == 1)
  template <typename TParticles>
  float getGenPsi2(TParticles const& particles) const
  {
    double qx = 0.0, qy = 0.0;
    for (const auto& particle : particles) {
      if (!particle.isPhysicalPrimary() || particle.eta() < mcEffEP.etaMin.value || particle.eta() > mcEffEP.etaMax.value) {
        continue;
      }
      const int absPdg = std::abs(particle.pdgCode());
      if (absPdg != PDG_t::kPiPlus && absPdg != PDG_t::kKPlus && absPdg != PDG_t::kProton && absPdg != PDG_t::kElectron && absPdg != PDG_t::kMuonMinus) {
        continue;
      }
      qx += std::cos(2.0 * particle.phi());
      qy += std::sin(2.0 * particle.phi());
    }
    return static_cast<float>(0.5 * std::atan2(qy, qx));
  }

  template <typename TMcCollision, typename TParticles>
  float getMCEffPsi(TMcCollision const& mcCollision, TParticles const& particles) const
  {
    if (mcEffEP.psiSource.value == 1) {
      return getGenPsi2(particles);
    }
    if (mcEffEP.psiSource.value == 2) {
      return mcCollision.eventPlaneAngle();
    }
    return 0.0f;
  }

  double getInjectionWeight(double cosThetaStar) const
  {
    if (!closure.injectPolarization.value) {
      return 1.0;
    }
    const double rho = closure.rho00Inj.value;
    return 1.5 * ((1.0 - rho) + (3.0 * rho - 1.0) * cosThetaStar * cosThetaStar);
  }

  float getClosurePsi(float trueEventPlaneAngle) { return getSAEventAngle(closure.useTrueEP.value ? trueEventPlaneAngle : 0.0f); }

  static constexpr int64_t SplitSampleModulus = 2;
  bool isClosureEffSample(int64_t mcCollisionIndex) const { return !closure.splitSample.value || (mcCollisionIndex % SplitSampleModulus == 0); }
  bool isClosureTestSample(int64_t mcCollisionIndex) const { return !closure.splitSample.value || (mcCollisionIndex % SplitSampleModulus == 1); }

  std::vector<int> countRecoPerMcCollision(aod::McCollisionLabels const& mcLabels, int64_t nMcCollisions) const
  {
    std::vector<int> nRecPerMc(nMcCollisions, 0);
    for (const auto& label : mcLabels) {
      if (label.mcCollisionId() >= 0 && label.mcCollisionId() < nMcCollisions) {
        nRecPerMc[label.mcCollisionId()]++;
      }
    }
    return nRecPerMc;
  }

  // ------------------------------------------------------------------
  // MC Truth Matching Helpers for K*(892)+- -> K0S + pi+-
  // Handles both direct (323 -> 310 + 211) and 2-step (323 -> 311 (-> 310) + 211)
  // decay trees in generators.
  // ------------------------------------------------------------------
  template <typename V0T, typename TrackT>
  bool matchRecoChKstar(const V0T& v0, const TrackT& bachelor, int64_t& outMotherIndex, ROOT::Math::PxPyPzMVector& outK0sGen, ROOT::Math::PxPyPzMVector& outPionGen) const
  {
    auto posTrack = v0.template posTrack_as<TrackMCRecTable>();
    auto negTrack = v0.template negTrack_as<TrackMCRecTable>();
    if (!posTrack.has_mcParticle() || !negTrack.has_mcParticle() || !bachelor.has_mcParticle()) {
      return false;
    }

    const auto mcPos = posTrack.mcParticle();
    const auto mcNeg = negTrack.mcParticle();
    const auto mcBach = bachelor.mcParticle();

    if (!mcBach.isPhysicalPrimary() || std::abs(mcBach.pdgCode()) != PDG_t::kPiPlus) {
      return false;
    }
    if (std::abs(mcPos.pdgCode()) != PDG_t::kPiPlus || std::abs(mcNeg.pdgCode()) != PDG_t::kPiPlus) {
      return false;
    }

    // Find common K0S (310) mother of V0 daughters
    for (const auto& momPos : mcPos.template mothers_as<aod::McParticles>()) {
      if (momPos.pdgCode() != PDG_t::kK0Short) {
        continue;
      }
      bool sharedK0s = false;
      for (const auto& momNeg : mcNeg.template mothers_as<aod::McParticles>()) {
        if (momNeg == momPos) {
          sharedK0s = true;
          break;
        }
      }
      if (!sharedK0s) {
        continue;
      }

      // Check mothers of K0S (either K*+- directly, or K0 (311) -> K*+-)
      for (const auto& k0sMom : momPos.template mothers_as<aod::McParticles>()) {
        if (std::abs(k0sMom.pdgCode()) == o2::constants::physics::Pdg::kKPlusStar892) {
          for (const auto& bachMom : mcBach.template mothers_as<aod::McParticles>()) {
            if (bachMom == k0sMom && std::abs(k0sMom.y()) < kstarCuts.cfgKstarRapidityMax.value) {
              outMotherIndex = k0sMom.globalIndex();
              outK0sGen = ROOT::Math::PxPyPzMVector(momPos.px(), momPos.py(), momPos.pz(), o2::constants::physics::MassK0Short);
              outPionGen = ROOT::Math::PxPyPzMVector(mcBach.px(), mcBach.py(), mcBach.pz(), o2::constants::physics::MassPionCharged);
              return true;
            }
          }
        } else if (std::abs(k0sMom.pdgCode()) == PDG_t::kK0) {
          for (const auto& k0GrandMom : k0sMom.template mothers_as<aod::McParticles>()) {
            if (std::abs(k0GrandMom.pdgCode()) != o2::constants::physics::Pdg::kKPlusStar892) {
              continue;
            }
            for (const auto& bachMom : mcBach.template mothers_as<aod::McParticles>()) {
              if (bachMom == k0GrandMom && std::abs(k0GrandMom.y()) < kstarCuts.cfgKstarRapidityMax.value) {
                outMotherIndex = k0GrandMom.globalIndex();
                outK0sGen = ROOT::Math::PxPyPzMVector(momPos.px(), momPos.py(), momPos.pz(), o2::constants::physics::MassK0Short);
                outPionGen = ROOT::Math::PxPyPzMVector(mcBach.px(), mcBach.py(), mcBach.pz(), o2::constants::physics::MassPionCharged);
                return true;
              }
            }
          }
        }
      }
    }
    return false;
  }

  // Extract generated K0S and bachelor pi+- daughters of a true K*+- mother
  template <typename McParticleT>
  bool extractGenChKstarDaughters(const McParticleT& mcMother, bool applyAccCut, ROOT::Math::PxPyPzMVector& outK0sGen, ROOT::Math::PxPyPzMVector& outPionGen, float& outCharge) const
  {
    if (mcEffEP.pdgcheck.value && std::abs(mcMother.pdgCode()) != o2::constants::physics::Pdg::kKPlusStar892) {
      return false;
    }
    if (std::abs(mcMother.y()) >= kstarCuts.cfgKstarRapidityMax.value) {
      return false;
    }

    const auto daughters = mcMother.template daughters_as<aod::McParticles>();
    if (daughters.size() != 2) {
      return false;
    }

    bool foundK0s = false;
    bool foundPion = false;
    for (const auto& dau : daughters) {
      const int absPdg = std::abs(dau.pdgCode());
      if (absPdg == PDG_t::kPiPlus) {
        if (!dau.isPhysicalPrimary()) {
          continue;
        }
        if (applyAccCut && (dau.pt() <= trackCuts.cfgCutPT.value || std::abs(dau.eta()) >= trackCuts.cfgCutEta.value)) {
          continue;
        }
        outPionGen = ROOT::Math::PxPyPzMVector(dau.px(), dau.py(), dau.pz(), o2::constants::physics::MassPionCharged);
        outCharge = (dau.pdgCode() > 0) ? 1.0f : -1.0f;
        foundPion = true;
      } else if (dau.pdgCode() == PDG_t::kK0Short) {
        if (applyAccCut && (dau.pt() <= k0sCuts.confV0PtMin.value || dau.pt() >= k0sCuts.confV0PtMax.value || std::abs(dau.y()) >= k0sCuts.confV0Rap.value)) {
          continue;
        }
        outK0sGen = ROOT::Math::PxPyPzMVector(dau.px(), dau.py(), dau.pz(), o2::constants::physics::MassK0Short);
        foundK0s = true;
      } else if (absPdg == PDG_t::kK0) {
        for (const auto& k0Dau : dau.template daughters_as<aod::McParticles>()) {
          if (k0Dau.pdgCode() == PDG_t::kK0Short) {
            if (applyAccCut && (k0Dau.pt() <= k0sCuts.confV0PtMin.value || k0Dau.pt() >= k0sCuts.confV0PtMax.value || std::abs(k0Dau.y()) >= k0sCuts.confV0Rap.value)) {
              continue;
            }
            outK0sGen = ROOT::Math::PxPyPzMVector(k0Dau.px(), k0Dau.py(), k0Dau.pz(), o2::constants::physics::MassK0Short);
            foundK0s = true;
            break;
          }
        }
      }
    }
    return foundK0s && foundPion;
  }

  // ------------------------------------------------------------------
  // Data-driven pool management helpers
  // ------------------------------------------------------------------
  static int findVariableBinDD(float x, const std::vector<double>& edges)
  {
    if (edges.size() < 2 || x < edges.front() || x >= edges.back()) {
      return -1;
    }
    const auto it = std::upper_bound(edges.begin(), edges.end(), static_cast<double>(x));
    return static_cast<int>(std::distance(edges.begin(), it)) - 1;
  }

  static int wrapBinDD(int b, int n)
  {
    if (n <= 0) {
      return -1;
    }
    b %= n;
    if (b < 0) {
      b += n;
    }
    return b;
  }

  static uint64_t splitmix64DD(uint64_t x)
  {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
  }

  static uint64_t linearKeyDD(int z, int c, int p, int e, int ph, int nC, int nP, int nE, int nPh)
  {
    auto k = static_cast<uint64_t>(z);
    k = k * static_cast<uint64_t>(nC) + static_cast<uint64_t>(c);
    k = k * static_cast<uint64_t>(nP) + static_cast<uint64_t>(p);
    k = k * static_cast<uint64_t>(nE) + static_cast<uint64_t>(e);
    k = k * static_cast<uint64_t>(nPh) + static_cast<uint64_t>(ph);
    return k;
  }

  static uint64_t linearKeyDDCharge(int z, int c, int p, int e, int ph, int q, int nC, int nP, int nE, int nPh) { return linearKeyDD(z, c, p, e, ph, nC, nP, nE, nPh) * 2ULL + static_cast<uint64_t>(q); }

  void pruneDDRollingPool(uint64_t currentDataFrame)
  {
    const int maxAge = dataDrivenMixing.cfgPoolMaxDataFrameAge.value;
    if (maxAge <= 0) {
      return;
    }

    for (auto it = ddRollingK0sPool.begin(); it != ddRollingK0sPool.end();) {
      auto& cell = it->second;
      while (!cell.empty() && currentDataFrame > cell.front().dataFrameId && currentDataFrame - cell.front().dataFrameId > static_cast<uint64_t>(maxAge)) {
        cell.pop_front();
        --ddK0sPoolCandidateCount;
      }
      if (cell.empty()) {
        it = ddRollingK0sPool.erase(it);
      } else {
        ++it;
      }
    }

    for (auto it = ddRollingPionPool.begin(); it != ddRollingPionPool.end();) {
      auto& cell = it->second;
      while (!cell.empty() && currentDataFrame > cell.front().dataFrameId && currentDataFrame - cell.front().dataFrameId > static_cast<uint64_t>(maxAge)) {
        cell.pop_front();
        --ddPionPoolCandidateCount;
      }
      if (cell.empty()) {
        it = ddRollingPionPool.erase(it);
      } else {
        ++it;
      }
    }
  }

  void insertDDK0sPoolCandidate(uint64_t key, const StoredK0sCandidate& candidate)
  {
    auto& cell = ddRollingK0sPool[key];
    cell.push_back(candidate);
    ++ddK0sPoolCandidateCount;
    const int maxDepth = dataDrivenMixing.cfgK0sPoolDepthPerBin.value;
    if (maxDepth > 0) {
      while (static_cast<int>(cell.size()) > maxDepth) {
        cell.pop_front();
        --ddK0sPoolCandidateCount;
      }
    }
  }

  void insertDDPionPoolCandidate(uint64_t key, const StoredPionCandidate& candidate)
  {
    auto& cell = ddRollingPionPool[key];
    cell.push_back(candidate);
    ++ddPionPoolCandidateCount;
    const int maxDepth = dataDrivenMixing.cfgPionPoolDepthPerBin.value;
    if (maxDepth > 0) {
      while (static_cast<int>(cell.size()) > maxDepth) {
        cell.pop_front();
        --ddPionPoolCandidateCount;
      }
    }
  }

  template <typename PoolT>
  void trimDDPoolToGlobalCap(PoolT& pool, size_t& totalCount, int configuredCap)
  {
    if (configuredCap <= 0 || totalCount <= static_cast<size_t>(configuredCap)) {
      return;
    }
    const auto cap = static_cast<size_t>(configuredCap);

    while (totalCount > cap) {
      uint64_t oldest = std::numeric_limits<uint64_t>::max();
      for (const auto& [key, cell] : pool) {
        (void)key;
        if (!cell.empty()) {
          oldest = std::min(oldest, cell.front().dataFrameId);
        }
      }
      if (oldest == std::numeric_limits<uint64_t>::max()) {
        break;
      }

      bool removed = false;
      for (auto it = pool.begin(); it != pool.end();) {
        auto& cell = it->second;
        while (!cell.empty() && cell.front().dataFrameId == oldest && totalCount > cap) {
          cell.pop_front();
          --totalCount;
          removed = true;
        }
        if (cell.empty()) {
          it = pool.erase(it);
        } else {
          ++it;
        }
        if (totalCount <= cap) {
          break;
        }
      }
      if (!removed) {
        break;
      }
    }
  }

  void enforceDDGlobalPoolCaps()
  {
    trimDDPoolToGlobalCap(ddRollingK0sPool, ddK0sPoolCandidateCount, dataDrivenMixing.cfgMaxTotalK0sPoolCandidates.value);
    trimDDPoolToGlobalCap(ddRollingPionPool, ddPionPoolCandidateCount, dataDrivenMixing.cfgMaxTotalPionPoolCandidates.value);
  }

  // ------------------------------------------------------------------
  // Selection Helpers (Event, Bachelor Pion, V0 K0S)
  // ------------------------------------------------------------------
  template <typename T>
  bool selectionTrack(const T& candidate) const
  {
    if (candidate.pt() <= trackCuts.cfgCutPT.value || std::abs(candidate.eta()) >= trackCuts.cfgCutEta.value) {
      return false;
    }
    if (!candidate.isPVContributor()) {
      return false;
    }
    if (candidate.itsNCls() <= trackCuts.cfgITScluster.value || candidate.tpcNClsCrossedRows() <= trackCuts.cfgTPCcluster.value) {
      return false;
    }
    if (trackCuts.useGlobalTrack.value && !candidate.isGlobalTrack()) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool selectionPID(const T& candidate, PIDSpecies species) const
  {
    float nSigmaTPC = 999.0f;
    float nSigmaTOF = 999.0f;
    float cutTPC = 0.0f;
    float cutCombined = 0.0f;

    switch (species) {
      case PIDSpecies::Pion:
        nSigmaTPC = candidate.tpcNSigmaPi();
        cutTPC = pidCuts.nsigmaCutTPCPi.value;
        cutCombined = pidCuts.nsigmaCutTOFPicomb.value;
        break;
      case PIDSpecies::Kaon:
        nSigmaTPC = candidate.tpcNSigmaKa();
        cutTPC = pidCuts.nsigmaCutTPCKa.value;
        cutCombined = pidCuts.nsigmaCutTOFKacomb.value;
        break;
      case PIDSpecies::Proton:
        nSigmaTPC = candidate.tpcNSigmaPr();
        cutTPC = pidCuts.nsigmaCutTPCPr.value;
        cutCombined = pidCuts.nsigmaCutTOFPrcomb.value;
        break;
    }

    if (!std::isfinite(nSigmaTPC) || std::abs(nSigmaTPC) >= cutTPC) {
      return false;
    }

    if (candidate.pt() < pidCuts.cfgPIDSwitchPt.value) {
      return true;
    }

    if (!candidate.hasTOF()) {
      return !pidCuts.cfgRequireTOFAboveSwitch.value;
    }

    switch (species) {
      case PIDSpecies::Pion:
        nSigmaTOF = candidate.tofNSigmaPi();
        break;
      case PIDSpecies::Kaon:
        nSigmaTOF = candidate.tofNSigmaKa();
        break;
      case PIDSpecies::Proton:
        nSigmaTOF = candidate.tofNSigmaPr();
        break;
    }

    if (!std::isfinite(nSigmaTOF)) {
      return false;
    }

    return std::hypot(nSigmaTPC, nSigmaTOF) < cutCombined;
  }

  template <typename T>
  bool selectionBachelorPion(const T& candidate) const
  {
    if (!selectionTrack(candidate) || !selectionPID(candidate, PIDSpecies::Pion)) {
      return false;
    }
    if (pidCuts.cfgUseExclusivePionPID.value && (selectionPID(candidate, PIDSpecies::Kaon) || selectionPID(candidate, PIDSpecies::Proton))) {
      return false;
    }
    return candidate.sign() != 0;
  }

  template <typename Collision, typename V0>
  bool selectionV0(const Collision& collision, const V0& candidate) const
  {
    const float pT = candidate.pt();
    const float tranRad = candidate.v0radius();
    const float dcaDaughV0 = std::abs(candidate.dcaV0daughters());
    const float cpaV0 = candidate.v0cosPA();
    const float ctauKShort = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;

    if (std::abs(candidate.dcav0topv()) >= k0sCuts.cMaxV0DCA.value) {
      return false;
    }
    if (pT <= k0sCuts.confV0PtMin.value || pT >= k0sCuts.confV0PtMax.value) {
      return false;
    }
    if (dcaDaughV0 >= k0sCuts.confV0DCADaughMax.value) {
      return false;
    }
    if (cpaV0 <= k0sCuts.confV0CPAMin.value) {
      return false;
    }
    if (tranRad <= k0sCuts.confV0TranRadV0Min.value || tranRad >= k0sCuts.confV0TranRadV0Max.value) {
      return false;
    }
    if (std::abs(ctauKShort) >= k0sCuts.cMaxV0LifeTime.value) {
      return false;
    }
    if (candidate.qtarm() <= k0sCuts.qtArmenterosMin.value * std::abs(candidate.alpha())) {
      return false;
    }
    if (std::abs(candidate.yK0Short()) >= k0sCuts.confV0Rap.value) {
      return false;
    }
    return true;
  }

  template <typename TrackTableT = FullTrackCandidates, typename V0>
  bool isSelectedV0Daughter(const V0& candidate) const
  {
    auto posTrack = candidate.template posTrack_as<TrackTableT>();
    auto negTrack = candidate.template negTrack_as<TrackTableT>();

    if (posTrack.sign() <= 0 || negTrack.sign() >= 0) {
      return false;
    }
    if (posTrack.tpcNClsCrossedRows() <= k0sCuts.confDaughTPCncrwsMin.value || negTrack.tpcNClsCrossedRows() <= k0sCuts.confDaughTPCncrwsMin.value) {
      return false;
    }
    if (posTrack.tpcNClsFound() <= k0sCuts.confDaughTPCnclsMin.value || negTrack.tpcNClsFound() <= k0sCuts.confDaughTPCnclsMin.value) {
      return false;
    }
    if (std::abs(posTrack.tpcNSigmaPi()) >= k0sCuts.confDaughPIDCuts.value || std::abs(negTrack.tpcNSigmaPi()) >= k0sCuts.confDaughPIDCuts.value) {
      return false;
    }
    if (candidate.positivept() <= k0sCuts.cfgDaughPiPt.value || candidate.negativept() <= k0sCuts.cfgDaughPiPt.value) {
      return false;
    }
    if (std::abs(candidate.positiveeta()) >= k0sCuts.confDaughEta.value || std::abs(candidate.negativeeta()) >= k0sCuts.confDaughEta.value) {
      return false;
    }
    if (std::abs(candidate.dcapostopv()) <= k0sCuts.cMinV0DCAPi.value || std::abs(candidate.dcanegtopv()) <= k0sCuts.cMinV0DCAPi.value) {
      return false;
    }
    return true;
  }

  template <typename TrackTableT = FullTrackCandidates, typename Collision, typename V0>
  bool selectionK0s(const Collision& collision, const V0& v0) const
  {
    if (!isSelectedV0Daughter<TrackTableT>(v0) || !selectionV0(collision, v0)) {
      return false;
    }
    return v0.mK0Short() >= k0sCuts.cfgK0sMassMin.value && v0.mK0Short() <= k0sCuts.cfgK0sMassMax.value;
  }

  template <typename Collision>
  bool passEventQuality(const Collision& collision, bool fillHist = false, bool checkRCT = true)
  {
    if (fillHist) {
      histos.fill(HIST("EventQA/hEventCut"), 0);
    }

    if (std::abs(collision.posZ()) > eventCuts.cfgCutVertexZ.value) {
      return false;
    }
    if (fillHist) {
      histos.fill(HIST("EventQA/hEventCut"), 1);
    }

    if (!collision.sel8()) {
      return false;
    }
    if (fillHist) {
      histos.fill(HIST("EventQA/hEventCut"), 2);
    }

    if (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (fillHist) {
      histos.fill(HIST("EventQA/hEventCut"), 3);
    }

    if (!collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (fillHist) {
      histos.fill(HIST("EventQA/hEventCut"), 4);
    }

    if (!collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (fillHist) {
      histos.fill(HIST("EventQA/hEventCut"), 5);
    }

    if (eventCuts.useNoSameBunchPileup.value && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (fillHist) {
      histos.fill(HIST("EventQA/hEventCut"), 6);
    }

    if (eventCuts.useGoodITSLayersAll.value && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }
    if (fillHist) {
      histos.fill(HIST("EventQA/hEventCut"), 7);
    }

    if (eventCuts.useNoCollInTimeRangeStandard.value && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (fillHist) {
      histos.fill(HIST("EventQA/hEventCut"), 8);
    }

    if (eventCuts.useOccCut.value && (std::abs(collision.trackOccupancyInTimeRange()) >= eventCuts.cfgOccCut.value)) {
      return false;
    }
    if (fillHist) {
      histos.fill(HIST("EventQA/hEventCut"), 9);
    }

    if (checkRCT && rctCut.requireRCTFlagChecker.value && !rctChecker(collision)) {
      return false;
    }
    if (fillHist) {
      histos.fill(HIST("EventQA/hEventCut"), 10);
    }

    if (eventCuts.useIsGoodZvtxFT0vsPV.value && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (fillHist) {
      histos.fill(HIST("EventQA/hEventCut"), 11);
    }

    if (eventCuts.useIsVertexITSTPC.value && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (fillHist) {
      histos.fill(HIST("EventQA/hEventCut"), 12);
    }

    if (eventCuts.useIsVertexTOFmatched.value && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if (fillHist) {
      histos.fill(HIST("EventQA/hEventCut"), 13);
    }

    if (eventCuts.useInelGt0.value && !collision.isInelGt0()) {
      return false;
    }
    if (fillHist) {
      histos.fill(HIST("EventQA/hEventCut"), 14);
    }

    if (eventCuts.useHasFT0.value && !collision.has_foundFT0()) {
      return false;
    }
    if (fillHist) {
      histos.fill(HIST("EventQA/hEventCut"), 15);
    }

    if (eventCuts.useHasBC.value && !collision.has_foundBC()) {
      return false;
    }
    if (fillHist) {
      histos.fill(HIST("EventQA/hEventCut"), 16);
    }

    return true;
  }

  template <typename Collision>
  bool passEventQualityMC(const Collision& collision)
  {
    if (collision.centFT0C() < eventCuts.cfgCutCentralityMin.value || collision.centFT0C() >= eventCuts.cfgCutCentralityMax.value) {
      return false;
    }
    return passEventQuality(collision, false, false);
  }

  // ------------------------------------------------------------------
  // init: Histogram and THnSparse Definitions
  // ------------------------------------------------------------------
  void init(o2::framework::InitContext&)
  {
    rctChecker.init(rctCut.cfgEvtRCTFlagCheckerLabel, rctCut.cfgEvtRCTFlagCheckerZDCCheck, rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    if (saConfig.cfgSAFrame.value < kEventPlane || saConfig.cfgSAFrame.value > kRandomEventPlane) {
      LOGF(fatal, "cfgSAFrame = %d not supported (0 = event plane, 1 = production plane, 2 = random event plane)", saConfig.cfgSAFrame.value);
    }
    if (saConfig.cfgSAFrame.value != kProductionPlane && saConfig.cfgEPNormalHarmonic.value != 1.0f) {
      LOGF(warning, "cfgEPNormalHarmonic = %.1f: the SA axis is not the normal to the Psi2 plane, so R2 resolution correction does not close", saConfig.cfgEPNormalHarmonic.value);
    }
    rndGen.SetSeed(saConfig.cfgRndSeed.value);

    const AxisSpec phiAxis = {500, -6.28, 6.28, "phi"};
    const AxisSpec resAxis = {6000, -30, 30, "Res"};
    const AxisSpec centAxis = {8, 0, 80, "FT0C (%)"};

    // Event QA
    histos.add("EventQA/hEvtSelInfo", "Event selection;selection step;events", kTH1F, {eventQAAxes.cfgAxisEvtSel});
    auto hEvtSel = histos.get<TH1>(HIST("EventQA/hEvtSelInfo"));
    hEvtSel->GetXaxis()->SetBinLabel(1, "z/cent prefilter");
    hEvtSel->GetXaxis()->SetBinLabel(2, "event quality");
    hEvtSel->GetXaxis()->SetBinLabel(3, "valid EP");
    hEvtSel->GetXaxis()->SetBinLabel(4, "accepted");

    histos.add("EventQA/hEventCut", "No. of event after cuts", kTH1D, {{20, 0, 20}});
    auto check = [](bool enabled) { return enabled ? "" : " #otimes"; };
    std::vector<std::string> eveCutLabels = {
      "All Events",
      Form("|Vz| < %.1f", eventCuts.cfgCutVertexZ.value),
      std::string("sel8"), std::string("kNoTimeFrameBorder"),
      std::string("kNoITSROFrameBorder"), std::string("kIsTriggerTVX"),
      std::string("kNoSameBunchPileup") + check(eventCuts.useNoSameBunchPileup.value),
      std::string("kIsGoodITSLayersAll") + check(eventCuts.useGoodITSLayersAll.value),
      std::string("kNoCollInTimeRangeStandard") + check(eventCuts.useNoCollInTimeRangeStandard.value),
      Form("Occupancy < %d%s", eventCuts.cfgOccCut.value, check(eventCuts.useOccCut.value)),
      std::string("rctChecker") + check(rctCut.requireRCTFlagChecker.value),
      std::string("kIsGoodZvtxFT0vsPV") + check(eventCuts.useIsGoodZvtxFT0vsPV.value),
      std::string("isVertexITSTPC") + check(eventCuts.useIsVertexITSTPC.value),
      std::string("isVertexTOFMatched") + check(eventCuts.useIsVertexTOFmatched.value),
      std::string("INEL > 0") + check(eventCuts.useInelGt0.value),
      std::string("hasFT0") + check(eventCuts.useHasFT0.value),
      std::string("hasBC") + check(eventCuts.useHasBC.value)};

    for (size_t i = 0; i < eveCutLabels.size(); ++i) {
      histos.get<TH1>(HIST("EventQA/hEventCut"))->GetXaxis()->SetBinLabel(i + 1, eveCutLabels[i].c_str());
    }

    histos.add("EventQA/hCent", "Accepted events;FT0C centrality (%);events", kTH1F, {eventQAAxes.cfgAxisCentQA});
    histos.add("EventQA/hVz", "Accepted events;z_{vtx} (cm);events", kTH1F, {eventQAAxes.cfgAxisVzQA});
    histos.add("EventQA/hPsiFT0C", "FT0C second-order event plane;#Psi_{2}^{FT0C};events", kTH1F, {eventQAAxes.cfgAxisPsiEventQA});
    histos.add("EventQA/pResFT0CFT0A", "FT0C-FT0A correlation;FT0C centrality (%);<#LTcos[2(#Psi_{2}^{C}-#Psi_{2}^{A})]#GT>", kTProfile, {eventQAAxes.cfgAxisCentQA});
    histos.add("EventQA/pResFT0CTPC", "FT0C-TPC correlation;FT0C centrality (%);<#LTcos[2(#Psi_{2}^{C}-#Psi_{2}^{TPC})]#GT>", kTProfile, {eventQAAxes.cfgAxisCentQA});
    histos.add("EventQA/pResFT0ATPC", "FT0A-TPC correlation;FT0C centrality (%);<#LTcos[2(#Psi_{2}^{A}-#Psi_{2}^{TPC})]#GT>", kTProfile, {eventQAAxes.cfgAxisCentQA});

    if (saConfig.additionalQAplots.value) {
      histos.add("EventQA/ResFT0CTPC", "ResFT0CTPC", kTH2F, {centAxis, resAxis});
      histos.add("EventQA/ResFT0CFT0A", "ResFT0CFT0A", kTH2F, {centAxis, resAxis});
      histos.add("EventQA/ResFT0ATPC", "ResFT0ATPC", kTH2F, {centAxis, resAxis});
      histos.add("EventQA/ResFT0CTPCSP", "ResFT0CTPCSP", kTH2F, {centAxis, resAxis});
      histos.add("EventQA/ResFT0CFT0ASP", "ResFT0CFT0ASP", kTH2F, {centAxis, resAxis});
      histos.add("EventQA/ResFT0ATPCSP", "ResFT0ATPCSP", kTH2F, {centAxis, resAxis});
      histos.add("CandidateQA/hK0sMass", "Selected V0 mass;M_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2});candidates", kTH1F, {{200, 0.43f, 0.57f}});
    }
    if (saConfig.cfgSAFrame.value == kRandomEventPlane) {
      histos.add("EventQA/hPsiRandom", "Random event plane angle;FT0C centrality (%);#Psi_{rnd}", kTH2F, {centAxis, phiAxis});
    }

    // Same-Event & Same-Event Rotational Background
    histos.add("hSparsesame", "Same-event K^{0}_{S}#pi^{#pm};M_{K^{0}_{S}#pi};p_{T};centrality;SA;charge", kTHnSparseF, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, sparseAxes.configThnAxisSA, sparseAxes.configThnAxisCharge});

    if (saConfig.fillRotation.value) {
      histos.add("hRotation", "Rotation angle;#Delta#varphi_{rot} (rad);entries", kTH1F, {{360, 0.0, o2::constants::math::TwoPI}});
      histos.add("hSparseRot", "Rotational background K^{0}_{S}#pi^{#pm};M_{K^{0}_{S}#pi};p_{T};centrality;SA;charge", kTHnSparseF, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, sparseAxes.configThnAxisSA, sparseAxes.configThnAxisCharge});
    }

    // Standard Mixed-Event Background
    histos.add("hSparseMENormal", "Standard mixed-event K^{0}_{S}#pi^{#pm};M_{K^{0}_{S}#pi};p_{T};centrality;SA;charge", kTHnSparseF, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, sparseAxes.configThnAxisSA, sparseAxes.configThnAxisCharge});

    histos.add("NormalMixQA/hDeltaVz", "Standard ME;z_{vtx}^{K^{0}_{S} event}-z_{vtx}^{#pi event} (cm);event pairs", kTH1F, {normalMixingQAAxes.cfgDeltaVz});
    histos.add("NormalMixQA/hDeltaCent", "Standard ME;cent_{K^{0}_{S} event}-cent_{#pi event} (%);event pairs", kTH1F, {normalMixingQAAxes.cfgDeltaCent});
    histos.add("NormalMixQA/hDeltaPsi2", "Standard ME;wrapped (#Psi_{2}^{K^{0}_{S} event}-#Psi_{2}^{#pi event});event pairs", kTH1F, {normalMixingQAAxes.cfgDeltaPsi});
    histos.add("NormalMixQA/hPairsVsCent", "Accepted standard mixed-event pairs;FT0C centrality of reference event (%);event pairs", kTH1F, {eventQAAxes.cfgAxisCentQA});

    // Data-Driven Daughter-Replacement Mixed-Event Background
    histos.add("hSparseMEK0sReplaced", "DD ME, K0S replaced;M_{K^{0}_{S}(B)#pi(A)};p_{T,mix};centrality_{A};SA_{mix};charge", kTHnSparseF, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, sparseAxes.configThnAxisSA, sparseAxes.configThnAxisCharge});

    histos.add("hSparseMEPionReplaced", "DD ME, pion replaced;M_{K^{0}_{S}(A)#pi(B)};p_{T,mix};centrality_{A};SA_{mix};charge", kTHnSparseF, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, sparseAxes.configThnAxisSA, sparseAxes.configThnAxisCharge});

    if (dataDrivenMixing.cfgFillMatchedSESparses.value) {
      histos.add("hSparseSEMatchedK0s", "SE targets with K0S replacement;M_{SE};p_{T,SE};centrality;SA_{SE};charge", kTHnSparseF, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, sparseAxes.configThnAxisSA, sparseAxes.configThnAxisCharge});
      histos.add("hSparseSEMatchedPion", "SE targets with pion replacement;M_{SE};p_{T,SE};centrality;SA_{SE};charge", kTHnSparseF, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, sparseAxes.configThnAxisSA, sparseAxes.configThnAxisCharge});
    }

    histos.add("DataDrivenQA/hNMatchesK0s", "K0S-replacement matches per SE target;N_{mix};targets", kTH1F, {dataDrivenMixingQAAxes.cfgNMatches});
    histos.add("DataDrivenQA/hNMatchesPion", "Pion-replacement matches per SE target;N_{mix};targets", kTH1F, {dataDrivenMixingQAAxes.cfgNMatches});
    histos.add("DataDrivenQA/hDeltaPtK0s", "K0S replacement;p_{T}^{B}-p_{T}^{A};weighted replacements", kTH1F, {dataDrivenMixingQAAxes.cfgDeltaPt});
    histos.add("DataDrivenQA/hDeltaEtaK0s", "K0S replacement;#eta^{B}-#eta^{A};weighted replacements", kTH1F, {dataDrivenMixingQAAxes.cfgDeltaEta});
    histos.add("DataDrivenQA/hDeltaPhiK0s", "K0S replacement;wrapped (#varphi^{B}-#varphi^{A});weighted replacements", kTH1F, {dataDrivenMixingQAAxes.cfgDeltaPhi});
    histos.add("DataDrivenQA/hDeltaPtPion", "Pion replacement;p_{T}^{B}-p_{T}^{A};weighted replacements", kTH1F, {dataDrivenMixingQAAxes.cfgDeltaPt});
    histos.add("DataDrivenQA/hDeltaEtaPion", "Pion replacement;#eta^{B}-#eta^{A};weighted replacements", kTH1F, {dataDrivenMixingQAAxes.cfgDeltaEta});
    histos.add("DataDrivenQA/hDeltaPhiPion", "Pion replacement;wrapped (#varphi^{B}-#varphi^{A});weighted replacements", kTH1F, {dataDrivenMixingQAAxes.cfgDeltaPhi});
    histos.add("DataDrivenQA/hDeltaVz", "DD event matching;z_{vtx}^{rep}-z_{vtx}^{target} (cm);weighted replacements", kTH1F, {dataDrivenMixingQAAxes.cfgDeltaVz});
    histos.add("DataDrivenQA/hDeltaCent", "DD event matching;cent_{rep}-cent_{target} (%);weighted replacements", kTH1F, {dataDrivenMixingQAAxes.cfgDeltaCent});
    histos.add("DataDrivenQA/hDeltaPsi2", "DD event-plane decorrelation;wrapped (#Psi_{2}^{rep}-#Psi_{2}^{target});weighted replacements", kTH1F, {dataDrivenMixingQAAxes.cfgDeltaPsi});
    histos.add("DataDrivenQA/pK0sReplacementFoundVsPt", "K0S replacement efficiency;p_{T,SE};fraction with >=1 replacement", kTProfile, {dataDrivenMixingQAAxes.cfgMatchEffPt});
    histos.add("DataDrivenQA/pK0sReplacementFoundVsCosTheta", "K0S replacement efficiency;cos#theta^{*}_{SE};fraction with >=1 replacement", kTProfile, {candidateQAAxes.cfgAxisCosThetaStarQA});
    histos.add("DataDrivenQA/pK0sReplacementFoundVsCent", "K0S replacement efficiency;centrality (%);fraction with >=1 replacement", kTProfile, {dataDrivenMixingQAAxes.cfgMatchEffCent});
    histos.add("DataDrivenQA/pK0sReplacementFoundVsMass", "K0S replacement efficiency;M_{SE} (GeV/#it{c}^{2});fraction with >=1 replacement", kTProfile, {dataDrivenMixingQAAxes.cfgMatchEffMass});
    histos.add("DataDrivenQA/pPionReplacementFoundVsPt", "Pion replacement efficiency;p_{T,SE};fraction with >=1 replacement", kTProfile, {dataDrivenMixingQAAxes.cfgMatchEffPt});
    histos.add("DataDrivenQA/pPionReplacementFoundVsCosTheta", "Pion replacement efficiency;cos#theta^{*}_{SE};fraction with >=1 replacement", kTProfile, {candidateQAAxes.cfgAxisCosThetaStarQA});
    histos.add("DataDrivenQA/pPionReplacementFoundVsCent", "Pion replacement efficiency;centrality (%);fraction with >=1 replacement", kTProfile, {dataDrivenMixingQAAxes.cfgMatchEffCent});
    histos.add("DataDrivenQA/pPionReplacementFoundVsMass", "Pion replacement efficiency;M_{SE} (GeV/#it{c}^{2});fraction with >=1 replacement", kTProfile, {dataDrivenMixingQAAxes.cfgMatchEffMass});
    histos.add("DataDrivenQA/pAnyReplacementFoundVsCent", "Any DD replacement found;centrality (%);fraction with >=1 replacement", kTProfile, {dataDrivenMixingQAAxes.cfgMatchEffCent});
    histos.add("DataDrivenQA/pAnyReplacementFoundVsMass", "Any DD replacement found;M_{SE} (GeV/#it{c}^{2});fraction with >=1 replacement", kTProfile, {dataDrivenMixingQAAxes.cfgMatchEffMass});
    histos.add("DataDrivenQA/hTargetK0sKinematics", "Successful target K0S;p_{T};#eta;#varphi", HistType::kTH3F, {dataDrivenMixing.cfgK0sPtBins, dataDrivenMixing.cfgK0sEtaBins, dataDrivenMixing.cfgK0sPhiBins});
    histos.add("DataDrivenQA/hReplacementK0sKinematics", "Replacement K0S;p_{T};#eta;#varphi", HistType::kTH3F, {dataDrivenMixing.cfgK0sPtBins, dataDrivenMixing.cfgK0sEtaBins, dataDrivenMixing.cfgK0sPhiBins});
    histos.add("DataDrivenQA/hTargetPionKinematics", "Successful target pion;p_{T};#eta;#varphi", HistType::kTH3F, {dataDrivenMixing.cfgPionPtBins, dataDrivenMixing.cfgPionEtaBins, dataDrivenMixing.cfgPionPhiBins});
    histos.add("DataDrivenQA/hReplacementPionKinematics", "Replacement pion;p_{T};#eta;#varphi", HistType::kTH3F, {dataDrivenMixing.cfgPionPtBins, dataDrivenMixing.cfgPionEtaBins, dataDrivenMixing.cfgPionPhiBins});
    histos.add("DataDrivenQA/hK0sPoolCandidates", "Rolling DD pool;stored K^{0}_{S} candidates;dataframes", kTH1F, {dataDrivenMixingQAAxes.cfgPoolCandidates});
    histos.add("DataDrivenQA/hPionPoolCandidates", "Rolling DD pool;stored bachelor-pion candidates;dataframes", kTH1F, {dataDrivenMixingQAAxes.cfgPoolCandidates});

    // MC Efficiency Histograms (processMC)
    if (doprocessMC) {
      histos.add("MC/hMC", "MC Event statistics", kTH1F, {{10, 0.0f, 10.0f}});
      histos.add("MC/CentPercentileMCRecHist", "MC Centrality;FT0C centrality (%);events", kTH1F, {{100, 0.0f, 100.0f}});
      histos.add("MC/hMCEffPsi", "Psi of the MC efficiency;centrality (%);#Psi", kTH2F, {sparseAxes.configThnAxisCentrality, {36, -o2::constants::math::PI, o2::constants::math::PI}});
      histos.add("MC/hSparseChKstarMCGenSA", "Generated K*(892)^{#pm};M_{K^{0}_{S}#pi};p_{T};centrality;SA;charge", kTHnSparseD, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, sparseAxes.configThnAxisSA, sparseAxes.configThnAxisCharge});
      histos.add("MC/hSparseChKstarMCGenCosThetaStar_effy", "Generated K*(892)^{#pm} cos#theta^{*};M_{K^{0}_{S}#pi};p_{T};centrality;cos#theta^{*};charge", kTHnSparseD, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, sparseAxes.configThnAxisSA, sparseAxes.configThnAxisCharge});
      histos.add("MC/hSparseChKstarMCRecSA", "Reconstructed true K*(892)^{#pm};M_{K^{0}_{S}#pi};p_{T};centrality;SA;charge", kTHnSparseD, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, sparseAxes.configThnAxisSA, sparseAxes.configThnAxisCharge});
      histos.add("MC/hSparseChKstarMCRecCosThetaStar_effy", "Reconstructed true K*(892)^{#pm} cos#theta^{*};M_{K^{0}_{S}#pi};p_{T};centrality;cos#theta^{*};charge", kTHnSparseD, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, sparseAxes.configThnAxisSA, sparseAxes.configThnAxisCharge});
      if (mcEffEP.fillPhiVMinusPsi.value) {
        const AxisSpec phiVMinusPsiAxis{mcEffEP.axisPhiVMinusPsi, "#varphi_{V} - #Psi"};
        histos.add("MC/hSparseChKstarMCGenCosThetaStarPhiVPsi", "Generated K*(892)^{#pm}: p_{T}, cos#theta^{*}, #varphi_{V} - #Psi, centrality", kTHnSparseD, {sparseAxes.configThnAxisPt, sparseAxes.configThnAxisSA, phiVMinusPsiAxis, sparseAxes.configThnAxisCentrality});
        histos.add("MC/hSparseChKstarMCRecCosThetaStarPhiVPsi", "Reconstructed K*(892)^{#pm}: p_{T}, cos#theta^{*}, #varphi_{V} - #Psi, centrality", kTHnSparseD, {sparseAxes.configThnAxisPt, sparseAxes.configThnAxisSA, phiVMinusPsiAxis, sparseAxes.configThnAxisCentrality});
      }
    }

    // Inline closure test histograms
    if (closure.doClosureTest.value) {
      histos.add("closure/ChKstarSameEventCosThetaStar", "ChK* closure: same-event;M;p_{T};centrality;cos#theta^{*};charge", kTHnSparseF, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, sparseAxes.configThnAxisSA, sparseAxes.configThnAxisCharge});
      histos.add("closure/ChKstarMixedEventCosThetaStar", "ChK* closure: mixed-event;M;p_{T};centrality;cos#theta^{*};charge", kTHnSparseF, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, sparseAxes.configThnAxisSA, sparseAxes.configThnAxisCharge});
      histos.add("closure/ChKstarMCRecTrueCosThetaStar", "ChK* closure: true matched reco;M;p_{T};centrality;cos#theta^{*};charge", kTHnSparseF, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, sparseAxes.configThnAxisSA, sparseAxes.configThnAxisCharge});
      histos.add("closure/ChKstarMCGenCosThetaStar", "ChK* closure: generated;M;p_{T};centrality;cos#theta^{*};charge", kTHnSparseF, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, sparseAxes.configThnAxisSA, sparseAxes.configThnAxisCharge});
    }

    // Event Loss and Signal Loss Histograms (processEvtLossSigLossMC)
    if (doprocessEvtLossSigLossMC) {
      const AxisSpec impactParAxis{evtSigLoss.axisImpactPar, "Impact parameter (fm)"};
      const AxisSpec multMCAxis{evtSigLoss.axisMultMC, "N_{ch}^{gen} (|#eta| < 0.5)"};
      const AxisSpec centLossAxis{100, 0.0, 100.0, "Centrality FT0C (%)"};
      const AxisSpec ptKstarAxis{sparseAxes.configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};

      histos.add("evtSigLoss/MCEventHist", "MC event statistics", kTH1F, {{4, 0.5, 4.5}});
      auto hstat = histos.get<TH1>(HIST("evtSigLoss/MCEventHist"));
      hstat->GetXaxis()->SetBinLabel(1, "All MC events");
      hstat->GetXaxis()->SetBinLabel(2, "MC events with selected reco event");
      hstat->GetXaxis()->SetBinLabel(3, "MC events with no reco event");
      hstat->GetXaxis()->SetBinLabel(4, "MC events with >1 reco event");

      histos.add("evtSigLoss/hImpactParameterGen", "Impact parameter, all generated events", kTH1F, {impactParAxis});
      histos.add("evtSigLoss/hImpactParameterRec", "Impact parameter, generated events with selected reco event", kTH1F, {impactParAxis});
      histos.add("evtSigLoss/hImpactParameterGenNoReco", "Impact parameter, generated events with no reco event", kTH1F, {impactParAxis});
      histos.add("evtSigLoss/hImpactParVsCentRec", "Impact parameter vs centrality of selected reco event", kTH2F, {centLossAxis, impactParAxis});
      histos.add("evtSigLoss/hMultEta05Gen", "N_{ch}^{gen}, all generated events", kTH1F, {multMCAxis});
      histos.add("evtSigLoss/hMultEta05Rec", "N_{ch}^{gen}, generated events with selected reco event", kTH1F, {multMCAxis});
      histos.add("evtSigLoss/hMultEta05GenNoReco", "N_{ch}^{gen}, generated events with no reco event", kTH1F, {multMCAxis});
      histos.add("evtSigLoss/hMultEta05VsCentRec", "N_{ch}^{gen} vs centrality of selected reco event", kTH2F, {centLossAxis, multMCAxis});

      histos.add("evtSigLoss/hChKstarGenBeforeEvtSel", "Generated K*+-, all events", kTH2F, {ptKstarAxis, impactParAxis});
      histos.add("evtSigLoss/hChKstarGenAfterEvtSel", "Generated K*+-, events with selected reco event", kTH2F, {ptKstarAxis, impactParAxis});
      histos.add("evtSigLoss/hChKstarGenVsMultBeforeEvtSel", "Generated K*+-, all events", kTH2F, {ptKstarAxis, multMCAxis});
      histos.add("evtSigLoss/hChKstarGenVsMultAfterEvtSel", "Generated K*+-, events with selected reco event", kTH2F, {ptKstarAxis, multMCAxis});
      histos.add("evtSigLoss/hChKstarGenAfterEvtSelVsCent", "Generated K*+-, events with selected reco event", kTH2F, {ptKstarAxis, centLossAxis});
    }

    // Dedicated MC Closure (True EP)
    if (doprocessMCClosure || doprocessMCClosureME) {
      if (closure.rho00Inj.value < 0.0f || closure.rho00Inj.value > 1.0f) {
        LOGF(fatal, "closure.rho00Inj = %f is outside [0,1]", closure.rho00Inj.value);
      }
      const AxisSpec cosAxisClosure{closure.axisCosThetaStar, "cos(#vartheta*)"};
      const AxisSpec cosFineAxis{100, -1.0, 1.0, "cos(#vartheta*)"};
      if (doprocessMCClosure) {
        histos.add("closureChKstar/hMC", "MC event statistics", kTH1F, {{5, 0.0f, 5.0f}});
        histos.add("closureChKstar/hSparseGenEff", "Generated K*+-, efficiency sample", kTHnSparseD, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, cosAxisClosure, sparseAxes.configThnAxisCharge}, true);
        histos.add("closureChKstar/hSparseRecTrueEff", "Reconstructed true K*+-, efficiency sample", kTHnSparseD, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, cosAxisClosure, sparseAxes.configThnAxisCharge}, true);
        histos.add("closureChKstar/hSparseGen", "Generated K*+- (input), test sample", kTHnSparseD, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, cosAxisClosure, sparseAxes.configThnAxisCharge}, true);
        histos.add("closureChKstar/hSparseRecTrue", "Reconstructed true K*+-, test sample", kTHnSparseD, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, cosAxisClosure, sparseAxes.configThnAxisCharge}, true);
        histos.add("closureChKstar/hSparseSE", "Same-event K0S+pi+- pairs, test sample", kTHnSparseD, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, cosAxisClosure, sparseAxes.configThnAxisCharge}, true);
        if (closure.fillRotBkg.value) {
          histos.add("closureChKstar/hSparseRot", "Rotated K0S+pi+- pairs, test sample", kTHnSparseD, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, cosAxisClosure, sparseAxes.configThnAxisCharge}, true);
        }
        histos.add("closureChKstar/hCosThetaStarGenVsRec", "cos(#vartheta*) gen vs rec, true K*+-", kTH2F, {cosFineAxis, cosFineAxis});
      }
      if (doprocessMCClosureME) {
        histos.add("closureChKstar/hSparseME", "Mixed-event K0S+pi+- pairs, test sample", kTHnSparseD, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, cosAxisClosure, sparseAxes.configThnAxisCharge}, true);
      }
    }

    // Dedicated MC Closure (Reconstructed FT0C EP)
    if (doprocessMCClosureRecoEP || doprocessMCClosureMERecoEP) {
      if (closure.rho00Inj.value < 0.0f || closure.rho00Inj.value > 1.0f) {
        LOGF(fatal, "closure.rho00Inj = %f is outside [0,1]", closure.rho00Inj.value);
      }
      const AxisSpec cosAxisClosure{closure.axisCosThetaStar, "cos(#vartheta*)"};
      const AxisSpec cosFineAxis{100, -1.0, 1.0, "cos(#vartheta*)"};
      if (doprocessMCClosureRecoEP) {
        histos.add("closureChKstarRecoEP/hMC", "MC event statistics", kTH1F, {{5, 0.0f, 5.0f}});
        histos.add("closureChKstarRecoEP/hSparseGenEff", "Generated K*+-, efficiency sample", kTHnSparseD, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, cosAxisClosure, sparseAxes.configThnAxisCharge}, true);
        histos.add("closureChKstarRecoEP/hSparseRecTrueEff", "Reconstructed true K*+-, efficiency sample", kTHnSparseD, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, cosAxisClosure, sparseAxes.configThnAxisCharge}, true);
        histos.add("closureChKstarRecoEP/hSparseGen", "Generated K*+- (input), test sample", kTHnSparseD, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, cosAxisClosure, sparseAxes.configThnAxisCharge}, true);
        histos.add("closureChKstarRecoEP/hSparseRecTrue", "Reconstructed true K*+-, test sample", kTHnSparseD, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, cosAxisClosure, sparseAxes.configThnAxisCharge}, true);
        histos.add("closureChKstarRecoEP/hSparseSE", "Same-event K0S+pi+- pairs, test sample", kTHnSparseD, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, cosAxisClosure, sparseAxes.configThnAxisCharge}, true);
        if (closure.fillRotBkg.value) {
          histos.add("closureChKstarRecoEP/hSparseRot", "Rotated K0S+pi+- pairs, test sample", kTHnSparseD, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, cosAxisClosure, sparseAxes.configThnAxisCharge}, true);
        }
        histos.add("closureChKstarRecoEP/hCosThetaStarGenVsRec", "cos(#vartheta*) gen vs rec, true K*+-", kTH2F, {cosFineAxis, cosFineAxis});
        histos.add("closureChKstarRecoEP/ResFT0CTPCSP", "ResFT0CTPCSP", kTH2F, {centAxis, resAxis});
        histos.add("closureChKstarRecoEP/ResFT0CFT0ASP", "ResFT0CFT0ASP", kTH2F, {centAxis, resAxis});
        histos.add("closureChKstarRecoEP/ResFT0ATPCSP", "ResFT0ATPCSP", kTH2F, {centAxis, resAxis});
        histos.add("closureChKstarRecoEP/ResFT0CTPC", "ResFT0CTPC", kTH2F, {centAxis, resAxis});
        histos.add("closureChKstarRecoEP/ResFT0CFT0A", "ResFT0CFT0A", kTH2F, {centAxis, resAxis});
        histos.add("closureChKstarRecoEP/ResFT0ATPC", "ResFT0ATPC", kTH2F, {centAxis, resAxis});
        histos.add("closureChKstarRecoEP/ResTrueFT0C", "ResTrueFT0C", kTH2F, {centAxis, resAxis});
      }
      if (doprocessMCClosureMERecoEP) {
        histos.add("closureChKstarRecoEP/hSparseME", "Mixed-event K0S+pi+- pairs, test sample", kTHnSparseD, {sparseAxes.configThnAxisInvMass, sparseAxes.configThnAxisPt, sparseAxes.configThnAxisCentrality, cosAxisClosure, sparseAxes.configThnAxisCharge}, true);
      }
    }
  }
  // ------------------------------------------------------------------
  // Same-event process: signal + rotational background
  // Supports Event Plane, Production Plane, and Random Event Plane
  // ------------------------------------------------------------------
  void processSE(EventCandidates::iterator const& collision, FullTrackCandidates const& tracks, ResoV0s const& v0s, aod::BCs const&)
  {
    const float centrality = collision.centFT0C();
    const float vz = collision.posZ();

    histos.fill(HIST("EventQA/hEvtSelInfo"), 0.5);

    if (!passEventQuality(collision, true)) {
      return;
    }

    histos.fill(HIST("EventQA/hEvtSelInfo"), 1.5);

    const float psiFT0C = collision.psiFT0C();
    if (!collision.triggereventep() || !std::isfinite(psiFT0C)) {
      return;
    }

    histos.fill(HIST("EventQA/hEvtSelInfo"), 2.5);
    histos.fill(HIST("EventQA/hCent"), centrality);
    histos.fill(HIST("EventQA/hVz"), vz);
    histos.fill(HIST("EventQA/hPsiFT0C"), psiFT0C);

    const float psiFT0A = collision.psiFT0A();
    const float psiTPC = collision.psiTPC();
    const float qFT0C = collision.qFT0C();
    const float qFT0A = collision.qFT0A();
    const float qTPC = collision.qTPC();

    if (qFT0A > 0.0f) {
      histos.fill(HIST("EventQA/pResFT0CFT0A"), centrality, std::cos(2.0f * (psiFT0C - psiFT0A)));
    }
    if (qTPC > 0.0f) {
      histos.fill(HIST("EventQA/pResFT0CTPC"), centrality, std::cos(2.0f * (psiFT0C - psiTPC)));
    }
    if (qFT0A > 0.0f && qTPC > 0.0f) {
      histos.fill(HIST("EventQA/pResFT0ATPC"), centrality, std::cos(2.0f * (psiFT0A - psiTPC)));
    }

    if (saConfig.additionalQAplots.value) {
      histos.fill(HIST("EventQA/ResFT0CTPC"), centrality, std::cos(2.0f * (psiFT0C - psiTPC)));
      histos.fill(HIST("EventQA/ResFT0CFT0A"), centrality, std::cos(2.0f * (psiFT0C - psiFT0A)));
      histos.fill(HIST("EventQA/ResFT0ATPC"), centrality, std::cos(2.0f * (psiTPC - psiFT0A)));
      histos.fill(HIST("EventQA/ResFT0CTPCSP"), centrality, qFT0C * qTPC * std::cos(2.0f * (psiFT0C - psiTPC)));
      histos.fill(HIST("EventQA/ResFT0CFT0ASP"), centrality, qFT0C * qFT0A * std::cos(2.0f * (psiFT0C - psiFT0A)));
      histos.fill(HIST("EventQA/ResFT0ATPCSP"), centrality, qTPC * qFT0A * std::cos(2.0f * (psiTPC - psiFT0A)));
    }

    const float psiSA = getSAEventAngle(psiFT0C);
    if (saConfig.cfgSAFrame.value == kRandomEventPlane) {
      histos.fill(HIST("EventQA/hPsiRandom"), centrality, psiSA);
    }

    histos.fill(HIST("EventQA/hEvtSelInfo"), 3.5);

    // Outer loop over V0s so topological selections are evaluated once per K0S
    for (const auto& v0 : v0s) {
      if (!selectionK0s<FullTrackCandidates>(collision, v0)) {
        continue;
      }

      if (saConfig.additionalQAplots.value) {
        histos.fill(HIST("CandidateQA/hK0sMass"), v0.mK0Short());
      }

      auto posTrack = v0.template posTrack_as<FullTrackCandidates>();
      auto negTrack = v0.template negTrack_as<FullTrackCandidates>();
      const auto posID = posTrack.globalIndex();
      const auto negID = negTrack.globalIndex();

      const ROOT::Math::PxPyPzMVector k0s(v0.px(), v0.py(), v0.pz(), o2::constants::physics::MassK0Short);

      for (const auto& bachelor : tracks) {
        if (!selectionBachelorPion(bachelor)) {
          continue;
        }

        const auto bachelorID = bachelor.globalIndex();
        if (bachelorID == posID || bachelorID == negID) {
          continue; // prevent self-pairing with a K0S daughter track
        }

        const ROOT::Math::PxPyPzMVector pionBach(bachelor.px(), bachelor.py(), bachelor.pz(), o2::constants::physics::MassPionCharged);
        const auto kstar = k0s + pionBach;

        const auto rapidity = static_cast<float>(kstar.Rapidity());
        if (!std::isfinite(rapidity) || std::abs(rapidity) >= kstarCuts.cfgKstarRapidityMax.value) {
          continue;
        }

        float saVal = 0.0f;
        float cosThetaStar = 0.0f;
        if (!computeSA(k0s, kstar, psiSA, saVal, &cosThetaStar)) {
          continue;
        }

        const auto mass = static_cast<float>(kstar.M());
        const auto pt = static_cast<float>(kstar.Pt());
        const auto charge = static_cast<float>(bachelor.sign());

        // Same-event signal candidate
        histos.fill(HIST("hSparsesame"), mass, pt, centrality, saVal, charge);

        if (closure.doClosureTest.value) {
          histos.fill(HIST("closure/ChKstarSameEventCosThetaStar"), mass, pt, centrality, cosThetaStar, charge);
        }

        // Same-event rotational combinatorial background
        if (saConfig.fillRotation.value && saConfig.nBkgRotations.value > 0) {
          const int nRot = saConfig.nBkgRotations.value;
          const float angleStart = saConfig.confMinRot.value;
          const float angleEnd = saConfig.confMaxRot.value;
          const float angleStep = (nRot > 1) ? ((angleEnd - angleStart) / static_cast<float>(nRot - 1)) : 0.0f;

          for (int irot = 0; irot < nRot; ++irot) {
            const float rotAngle = angleStart + static_cast<float>(irot) * angleStep;
            histos.fill(HIST("hRotation"), rotAngle);

            const float cosA = std::cos(rotAngle);
            const float sinA = std::sin(rotAngle);

            ROOT::Math::PxPyPzMVector k0sRot = k0s;
            ROOT::Math::PxPyPzMVector pionRot = pionBach;
            if (saConfig.cfgRotatePion.value) {
              const float rotPx = bachelor.px() * cosA - bachelor.py() * sinA;
              const float rotPy = bachelor.px() * sinA + bachelor.py() * cosA;
              pionRot = ROOT::Math::PxPyPzMVector(rotPx, rotPy, bachelor.pz(), o2::constants::physics::MassPionCharged);
            } else {
              const float rotPx = v0.px() * cosA - v0.py() * sinA;
              const float rotPy = v0.px() * sinA + v0.py() * cosA;
              k0sRot = ROOT::Math::PxPyPzMVector(rotPx, rotPy, v0.pz(), o2::constants::physics::MassK0Short);
            }

            const auto kstarRot = k0sRot + pionRot;
            const auto rapRot = static_cast<float>(kstarRot.Rapidity());
            if (!std::isfinite(rapRot) || std::abs(rapRot) >= kstarCuts.cfgKstarRapidityMax.value) {
              continue;
            }

            float saRot = 0.0f;
            if (!computeSA(k0sRot, kstarRot, psiSA, saRot)) {
              continue;
            }

            histos.fill(HIST("hSparseRot"), static_cast<float>(kstarRot.M()), static_cast<float>(kstarRot.Pt()), centrality, saRot, charge);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(chKstarSpin, processSE, "Process Same event", true);

  // ------------------------------------------------------------------
  // Conventional O2 dataframe mixed-event process
  // ------------------------------------------------------------------
  void processNormalME(EventCandidates const& /*collisions*/, FullTrackCandidates const& /*tracks*/, ResoV0s const& /*v0s*/, aod::BCsWithTimestamps const&)
  {
    for (const auto& [collision1, tracks1, collision2, v0s2] : normalMixedEventPairs) {
      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }

      if (!passEventQuality(collision1, false) || !passEventQuality(collision2, false)) {
        continue;
      }

      const float psi1 = collision1.psiFT0C();
      const float psi2 = collision2.psiFT0C();
      if (!collision1.triggereventep() || !collision2.triggereventep() || !std::isfinite(psi1) || !std::isfinite(psi2)) {
        continue;
      }

      if (normalMixing.cfgRequireSameRun.value) {
        const auto bc1 = collision1.template bc_as<aod::BCsWithTimestamps>();
        const auto bc2 = collision2.template bc_as<aod::BCsWithTimestamps>();
        if (bc1.runNumber() != bc2.runNumber()) {
          continue;
        }
      }

      const float centrality = collision1.centFT0C();
      const float psiSA1 = getSAEventAngle(psi1);

      histos.fill(HIST("NormalMixQA/hDeltaVz"), collision2.posZ() - collision1.posZ());
      histos.fill(HIST("NormalMixQA/hDeltaCent"), collision2.centFT0C() - collision1.centFT0C());
      histos.fill(HIST("NormalMixQA/hDeltaPsi2"), wrappedDeltaPsi2DD(psi2, psi1));
      histos.fill(HIST("NormalMixQA/hPairsVsCent"), centrality);

      for (const auto& v0 : v0s2) {
        if (!selectionK0s<FullTrackCandidates>(collision2, v0)) {
          continue;
        }

        const ROOT::Math::PxPyPzMVector k0s(v0.px(), v0.py(), v0.pz(), o2::constants::physics::MassK0Short);

        for (const auto& bachelor : tracks1) {
          if (!selectionBachelorPion(bachelor)) {
            continue;
          }

          const ROOT::Math::PxPyPzMVector pionBach(bachelor.px(), bachelor.py(), bachelor.pz(), o2::constants::physics::MassPionCharged);
          const auto kstar = k0s + pionBach;

          const auto rapidity = static_cast<float>(kstar.Rapidity());
          if (!std::isfinite(rapidity) || std::abs(rapidity) >= kstarCuts.cfgKstarRapidityMax.value) {
            continue;
          }

          float saVal = 0.0f;
          float cosThetaStar = 0.0f;
          if (!computeSA(k0s, kstar, psiSA1, saVal, &cosThetaStar)) {
            continue;
          }

          const auto mass = static_cast<float>(kstar.M());
          const auto pt = static_cast<float>(kstar.Pt());
          const auto charge = static_cast<float>(bachelor.sign());

          histos.fill(HIST("hSparseMENormal"), mass, pt, centrality, saVal, charge);

          if (closure.doClosureTest.value) {
            histos.fill(HIST("closure/ChKstarMixedEventCosThetaStar"), mass, pt, centrality, cosThetaStar, charge);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(chKstarSpin, processNormalME, "Process standard O2 dataframe mixed events", false);

  // ------------------------------------------------------------------
  // Data-driven daughter-replacement mixed-event process
  // ------------------------------------------------------------------
  void processDataDrivenME(EventCandidates const& collisions, FullTrackCandidates const& tracks, ResoV0s const& v0s, aod::BCsWithTimestamps const&)
  {
    const uint64_t currentDataFrame = ++ddDataFrameCounter;
    const int pruneEvery = std::max(1, dataDrivenMixing.cfgPoolPruneEveryNDataFrames.value);
    if ((currentDataFrame % static_cast<uint64_t>(pruneEvery)) == 0) {
      pruneDDRollingPool(currentDataFrame);
    }

    const auto& zEdges = dataDrivenMixing.cfgEventVtxBins.value;
    const auto& centEdges = dataDrivenMixing.cfgEventCentBins.value;
    const auto& k0PtEdges = dataDrivenMixing.cfgK0sPtBins.value;
    const auto& k0EtaEdges = dataDrivenMixing.cfgK0sEtaBins.value;
    const auto& k0PhiEdges = dataDrivenMixing.cfgK0sPhiBins.value;
    const auto& piPtEdges = dataDrivenMixing.cfgPionPtBins.value;
    const auto& piEtaEdges = dataDrivenMixing.cfgPionEtaBins.value;
    const auto& piPhiEdges = dataDrivenMixing.cfgPionPhiBins.value;

    const int nC = static_cast<int>(centEdges.size()) - 1;
    const int nK0Pt = static_cast<int>(k0PtEdges.size()) - 1;
    const int nK0Eta = static_cast<int>(k0EtaEdges.size()) - 1;
    const int nK0Phi = static_cast<int>(k0PhiEdges.size()) - 1;
    const int nPiPt = static_cast<int>(piPtEdges.size()) - 1;
    const int nPiEta = static_cast<int>(piEtaEdges.size()) - 1;
    const int nPiPhi = static_cast<int>(piPhiEdges.size()) - 1;

    if (nC <= 0 || nK0Pt <= 0 || nK0Eta <= 0 || nK0Phi <= 0 || nPiPt <= 0 || nPiEta <= 0 || nPiPhi <= 0) {
      return;
    }

    struct CurrentK0s {
      int64_t posId = -1;
      int64_t negId = -1;
      float px = 0.0f;
      float py = 0.0f;
      float pz = 0.0f;
      float pt = 0.0f;
      float eta = 0.0f;
      float phi = 0.0f;
    };

    struct CurrentPion {
      int64_t trackId = -1;
      int8_t charge = 0;
      float px = 0.0f;
      float py = 0.0f;
      float pz = 0.0f;
      float pt = 0.0f;
      float eta = 0.0f;
      float phi = 0.0f;
    };

    struct CurrentEvent {
      uint64_t eventKey = 0;
      int runNumber = -1;
      float psi2 = 0.0f;
      float psiSA = 0.0f;
      float posZ = 0.0f;
      float centrality = 0.0f;
      int zBin = -1;
      int centBin = -1;
      std::vector<CurrentK0s> k0s;
      std::vector<CurrentPion> pions;
    };

    struct MixedMatch {
      uint64_t eventKey = 0;
      int runNumber = -1;
      float replacementPsi2 = 0.0f;
      float replacementVz = 0.0f;
      float replacementCent = 0.0f;
      float mass = 0.0f;
      float pt = 0.0f;
      float saVal = 0.0f;
      float dPt = 0.0f;
      float dEta = 0.0f;
      float dPhi = 0.0f;
      float repPt = 0.0f;
      float repEta = 0.0f;
      float repPhi = 0.0f;
    };

    auto limitEventContribution = [&](auto& byCell, int maxPerEvent, uint64_t eventKey, uint64_t salt) {
      if (maxPerEvent <= 0) {
        return;
      }
      for (auto& [key, vec] : byCell) {
        if (static_cast<int>(vec.size()) <= maxPerEvent) {
          continue;
        }
        uint64_t state = splitmix64DD(dataDrivenMixing.cfgMixSeed.value ^ eventKey ^ key ^ salt);
        for (size_t i = vec.size(); i > 1; --i) {
          state = splitmix64DD(state);
          const auto j = static_cast<size_t>(state % i);
          std::swap(vec[i - 1], vec[j]);
        }
        vec.resize(static_cast<size_t>(maxPerEvent));
      }
    };

    auto downsampleMatches = [&](std::vector<MixedMatch> matches, uint64_t seed) {
      const int maxMatches = dataDrivenMixing.cfgMaxMatchesPerLeg.value;
      if (maxMatches <= 0 || static_cast<int>(matches.size()) <= maxMatches) {
        return matches;
      }
      uint64_t state = splitmix64DD(seed);
      for (size_t i = matches.size(); i > 1; --i) {
        state = splitmix64DD(state);
        const auto j = static_cast<size_t>(state % i);
        std::swap(matches[i - 1], matches[j]);
      }
      matches.resize(static_cast<size_t>(maxMatches));
      return matches;
    };

    auto keepDistinctEvents = [&](std::vector<MixedMatch>& matches, uint64_t seed) {
      const int maxEvents = dataDrivenMixing.cfgMaxDistinctEventsPerLeg.value;
      const int maxPerEvent = dataDrivenMixing.cfgMaxMatchesPerEventPerLeg.value;
      if (matches.empty()) {
        return;
      }

      uint64_t state = splitmix64DD(seed);
      for (size_t i = matches.size(); i > 1; --i) {
        state = splitmix64DD(state);
        const auto j = static_cast<size_t>(state % i);
        std::swap(matches[i - 1], matches[j]);
      }

      if (maxEvents <= 0 && maxPerEvent <= 0) {
        return;
      }
      std::unordered_map<uint64_t, int> counts;
      std::vector<MixedMatch> filtered;
      filtered.reserve(matches.size());
      for (const auto& m : matches) {
        auto& count = counts[m.eventKey];
        if (maxPerEvent > 0 && count >= maxPerEvent) {
          continue;
        }
        if (maxEvents > 0 && count == 0 && static_cast<int>(counts.size()) >= maxEvents) {
          continue;
        }
        ++count;
        filtered.push_back(m);
      }
      matches.swap(filtered);
    };

    auto newestCellBeginOffset = [](size_t size, int maxScan) -> size_t {
      if (maxScan <= 0 || size <= static_cast<size_t>(maxScan)) {
        return 0;
      }
      return size - static_cast<size_t>(maxScan);
    };

    std::vector<CurrentEvent> events;
    events.reserve(collisions.size());

    // PASS 1: populate rolling pools from current dataframe
    for (const auto& collision : collisions) {
      if (!passEventQuality(collision, false)) {
        continue;
      }
      const float psi2 = collision.psiFT0C();
      if (!collision.triggereventep() || !std::isfinite(psi2)) {
        continue;
      }

      const auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      const int runNumber = bc.runNumber();
      const float centrality = collision.centFT0C();
      const float posZ = collision.posZ();
      const int zBin = findVariableBinDD(posZ, zEdges);
      const int centBin = findVariableBinDD(centrality, centEdges);
      if (zBin < 0 || centBin < 0) {
        continue;
      }

      const uint64_t eventKey = (currentDataFrame << 32) | (static_cast<uint64_t>(collision.globalIndex()) & 0xffffffffULL);

      CurrentEvent evt;
      evt.eventKey = eventKey;
      evt.runNumber = runNumber;
      evt.psi2 = psi2;
      evt.psiSA = getSAEventAngle(psi2);
      evt.posZ = posZ;
      evt.centrality = centrality;
      evt.zBin = zBin;
      evt.centBin = centBin;

      const auto tracksThis = tracks.sliceBy(tracksPerCollisionDD, collision.globalIndex());
      const auto v0sThis = v0s.sliceBy(v0sPerCollisionDD, collision.globalIndex());

      std::unordered_map<uint64_t, std::vector<StoredK0sCandidate>> k0ByCell;
      std::unordered_map<uint64_t, std::vector<StoredPionCandidate>> pionByCell;

      evt.k0s.reserve(v0sThis.size());
      for (const auto& v0 : v0sThis) {
        if (!selectionK0s<FullTrackCandidates>(collision, v0)) {
          continue;
        }
        const ROOT::Math::PxPyPzMVector k0s(v0.px(), v0.py(), v0.pz(), o2::constants::physics::MassK0Short);
        const auto pt = static_cast<float>(k0s.Pt());
        const auto eta = static_cast<float>(k0s.Eta());
        const float phi = wrapPhiDD(static_cast<float>(k0s.Phi()));
        if (!std::isfinite(pt) || !std::isfinite(eta) || !std::isfinite(phi)) {
          continue;
        }

        auto posTrack = v0.template posTrack_as<FullTrackCandidates>();
        auto negTrack = v0.template negTrack_as<FullTrackCandidates>();
        evt.k0s.push_back(CurrentK0s{posTrack.globalIndex(), negTrack.globalIndex(), v0.px(), v0.py(), v0.pz(), pt, eta, phi});

        if (dataDrivenMixing.cfgDoK0sReplacement.value) {
          const int ip = findVariableBinDD(pt, k0PtEdges);
          const int ie = findVariableBinDD(eta, k0EtaEdges);
          const int iph = findVariableBinDD(phi, k0PhiEdges);
          if (ip >= 0 && ie >= 0 && iph >= 0) {
            const uint64_t key = linearKeyDD(zBin, centBin, ip, ie, iph, nC, nK0Pt, nK0Eta, nK0Phi);
            k0ByCell[key].push_back(StoredK0sCandidate{eventKey, currentDataFrame, runNumber, psi2, posZ, centrality, pt, eta, phi, v0.px(), v0.py(), v0.pz()});
          }
        }
      }

      evt.pions.reserve(tracksThis.size());
      for (const auto& track : tracksThis) {
        if (!selectionBachelorPion(track)) {
          continue;
        }
        const auto charge = static_cast<int8_t>(track.sign());
        if (charge == 0) {
          continue;
        }
        const float pt = track.pt();
        const float eta = track.eta();
        const float phi = wrapPhiDD(track.phi());
        if (!std::isfinite(pt) || !std::isfinite(eta) || !std::isfinite(phi)) {
          continue;
        }

        evt.pions.push_back(CurrentPion{track.globalIndex(), charge, track.px(), track.py(), track.pz(), pt, eta, phi});

        if (dataDrivenMixing.cfgDoPionReplacement.value) {
          const int ip = findVariableBinDD(pt, piPtEdges);
          const int ie = findVariableBinDD(eta, piEtaEdges);
          const int iph = findVariableBinDD(phi, piPhiEdges);
          if (ip >= 0 && ie >= 0 && iph >= 0) {
            const int qBin = charge > 0 ? 1 : 0;
            const uint64_t key = linearKeyDDCharge(zBin, centBin, ip, ie, iph, qBin, nC, nPiPt, nPiEta, nPiPhi);
            pionByCell[key].push_back(StoredPionCandidate{eventKey, currentDataFrame, runNumber, psi2, posZ, centrality, charge, track.globalIndex(), pt, eta, phi, track.px(), track.py(), track.pz()});
          }
        }
      }

      limitEventContribution(k0ByCell, dataDrivenMixing.cfgMaxK0sPerEventPerBin.value, eventKey, 0x11111111ULL);
      limitEventContribution(pionByCell, dataDrivenMixing.cfgMaxPionsPerEventPerBin.value, eventKey, 0x22222222ULL);

      for (const auto& [key, vec] : k0ByCell) {
        for (const auto& candidate : vec) {
          insertDDK0sPoolCandidate(key, candidate);
        }
      }
      for (const auto& [key, vec] : pionByCell) {
        for (const auto& candidate : vec) {
          insertDDPionPoolCandidate(key, candidate);
        }
      }

      events.push_back(std::move(evt));
    }

    enforceDDGlobalPoolCaps();

    // PASS 2: construct SE targets and perform daughter replacements
    uint64_t targetOrdinal = 0;
    for (const auto& evt : events) {
      for (const auto& targetK0 : evt.k0s) {
        const ROOT::Math::PxPyPzMVector targetK0s(targetK0.px, targetK0.py, targetK0.pz, o2::constants::physics::MassK0Short);

        const int targetK0PtBin = findVariableBinDD(targetK0.pt, k0PtEdges);
        const int targetK0EtaBin = findVariableBinDD(targetK0.eta, k0EtaEdges);
        const int targetK0PhiBin = findVariableBinDD(targetK0.phi, k0PhiEdges);

        for (const auto& targetPi : evt.pions) {
          if (targetPi.trackId == targetK0.posId || targetPi.trackId == targetK0.negId) {
            continue;
          }

          const ROOT::Math::PxPyPzMVector targetPion(targetPi.px, targetPi.py, targetPi.pz, o2::constants::physics::MassPionCharged);
          const auto targetMother = targetK0s + targetPion;
          const auto targetMass = static_cast<float>(targetMother.M());
          const auto targetPt = static_cast<float>(targetMother.Pt());
          const auto targetY = static_cast<float>(targetMother.Rapidity());
          if (!std::isfinite(targetMass) || !std::isfinite(targetPt) || !std::isfinite(targetY)) {
            continue;
          }
          if (std::abs(targetY) >= kstarCuts.cfgKstarRapidityMax.value || targetMass < dataDrivenMixing.cfgCandidateMassMin.value || targetMass >= dataDrivenMixing.cfgCandidateMassMax.value) {
            continue;
          }

          float targetSA = 0.0f;
          float targetCosThetaStar = 0.0f;
          if (!computeSA(targetK0s, targetMother, evt.psiSA, targetSA, &targetCosThetaStar)) {
            continue;
          }
          const auto charge = static_cast<float>(targetPi.charge);

          // Leg 1: replace K0S
          std::vector<MixedMatch> k0Matches;
          if (dataDrivenMixing.cfgDoK0sReplacement.value && targetK0PtBin >= 0 && targetK0EtaBin >= 0 && targetK0PhiBin >= 0) {
            const int dPtBins = std::max(0, dataDrivenMixing.cfgNeighborK0sPt.value);
            const int dEtaBins = std::max(0, dataDrivenMixing.cfgNeighborK0sEta.value);
            const int dPhiBins = std::max(0, dataDrivenMixing.cfgNeighborK0sPhi.value);

            for (int dp = -dPtBins; dp <= dPtBins; ++dp) {
              const int ip = targetK0PtBin + dp;
              if (ip < 0 || ip >= nK0Pt) {
                continue;
              }
              for (int de = -dEtaBins; de <= dEtaBins; ++de) {
                const int ie = targetK0EtaBin + de;
                if (ie < 0 || ie >= nK0Eta) {
                  continue;
                }
                for (int dph = -dPhiBins; dph <= dPhiBins; ++dph) {
                  const int iph = wrapBinDD(targetK0PhiBin + dph, nK0Phi);
                  const uint64_t key = linearKeyDD(evt.zBin, evt.centBin, ip, ie, iph, nC, nK0Pt, nK0Eta, nK0Phi);
                  const auto poolIt = ddRollingK0sPool.find(key);
                  if (poolIt == ddRollingK0sPool.end()) {
                    continue;
                  }
                  const auto& cell = poolIt->second;
                  const size_t begin = newestCellBeginOffset(cell.size(), dataDrivenMixing.cfgMaxK0sCandidatesScannedPerCell.value);
                  for (size_t idx = begin; idx < cell.size(); ++idx) {
                    const auto& rep = cell[idx];
                    if (rep.eventKey == evt.eventKey || rep.runNumber != evt.runNumber) {
                      continue;
                    }
                    const float dpt = rep.pt - targetK0.pt;
                    const float deta = rep.eta - targetK0.eta;
                    const float dphi = wrappedDeltaPhiDD(rep.phi, targetK0.phi);
                    if (std::abs(dpt) > dataDrivenMixing.cfgK0sMaxDeltaPt.value || std::abs(deta) > dataDrivenMixing.cfgK0sMaxDeltaEta.value || std::abs(dphi) > dataDrivenMixing.cfgK0sMaxDeltaPhi.value) {
                      continue;
                    }

                    const ROOT::Math::PxPyPzMVector repK0s(rep.px, rep.py, rep.pz, o2::constants::physics::MassK0Short);
                    const auto mixedMother = repK0s + targetPion;
                    const auto mass = static_cast<float>(mixedMother.M());
                    const auto pt = static_cast<float>(mixedMother.Pt());
                    const auto y = static_cast<float>(mixedMother.Rapidity());
                    if (!std::isfinite(mass) || !std::isfinite(pt) || !std::isfinite(y)) {
                      continue;
                    }
                    if (std::abs(y) >= kstarCuts.cfgKstarRapidityMax.value || mass < dataDrivenMixing.cfgCandidateMassMin.value || mass >= dataDrivenMixing.cfgCandidateMassMax.value) {
                      continue;
                    }

                    float saMix = 0.0f;
                    if (!computeSA(repK0s, mixedMother, evt.psiSA, saMix)) {
                      continue;
                    }
                    k0Matches.push_back(MixedMatch{rep.eventKey, rep.runNumber, rep.psi2, rep.posZ, rep.centrality, mass, pt, saMix, dpt, deta, dphi, rep.pt, rep.eta, rep.phi});
                  }
                }
              }
            }
          }

          const uint64_t k0Seed = dataDrivenMixing.cfgMixSeed.value ^ evt.eventKey ^ splitmix64DD(++targetOrdinal) ^ 0x4b3053ULL;
          keepDistinctEvents(k0Matches, k0Seed);
          k0Matches = downsampleMatches(std::move(k0Matches), k0Seed ^ 0x1111ULL);

          // Leg 2: replace bachelor pion
          std::vector<MixedMatch> pionMatches;
          const int targetPiPtBin = findVariableBinDD(targetPi.pt, piPtEdges);
          const int targetPiEtaBin = findVariableBinDD(targetPi.eta, piEtaEdges);
          const int targetPiPhiBin = findVariableBinDD(targetPi.phi, piPhiEdges);
          if (dataDrivenMixing.cfgDoPionReplacement.value && targetPiPtBin >= 0 && targetPiEtaBin >= 0 && targetPiPhiBin >= 0) {
            const int dPtBins = std::max(0, dataDrivenMixing.cfgNeighborPionPt.value);
            const int dEtaBins = std::max(0, dataDrivenMixing.cfgNeighborPionEta.value);
            const int dPhiBins = std::max(0, dataDrivenMixing.cfgNeighborPionPhi.value);
            const int qBin = targetPi.charge > 0 ? 1 : 0;

            for (int dp = -dPtBins; dp <= dPtBins; ++dp) {
              const int ip = targetPiPtBin + dp;
              if (ip < 0 || ip >= nPiPt) {
                continue;
              }
              for (int de = -dEtaBins; de <= dEtaBins; ++de) {
                const int ie = targetPiEtaBin + de;
                if (ie < 0 || ie >= nPiEta) {
                  continue;
                }
                for (int dph = -dPhiBins; dph <= dPhiBins; ++dph) {
                  const int iph = wrapBinDD(targetPiPhiBin + dph, nPiPhi);
                  const uint64_t key = linearKeyDDCharge(evt.zBin, evt.centBin, ip, ie, iph, qBin, nC, nPiPt, nPiEta, nPiPhi);
                  const auto poolIt = ddRollingPionPool.find(key);
                  if (poolIt == ddRollingPionPool.end()) {
                    continue;
                  }
                  const auto& cell = poolIt->second;
                  const size_t begin = newestCellBeginOffset(cell.size(), dataDrivenMixing.cfgMaxPionCandidatesScannedPerCell.value);
                  for (size_t idx = begin; idx < cell.size(); ++idx) {
                    const auto& rep = cell[idx];
                    if (rep.eventKey == evt.eventKey || rep.runNumber != evt.runNumber || rep.charge != targetPi.charge) {
                      continue;
                    }
                    const float dpt = rep.pt - targetPi.pt;
                    const float deta = rep.eta - targetPi.eta;
                    const float dphi = wrappedDeltaPhiDD(rep.phi, targetPi.phi);
                    if (std::abs(dpt) > dataDrivenMixing.cfgPionMaxDeltaPt.value || std::abs(deta) > dataDrivenMixing.cfgPionMaxDeltaEta.value || std::abs(dphi) > dataDrivenMixing.cfgPionMaxDeltaPhi.value) {
                      continue;
                    }

                    const ROOT::Math::PxPyPzMVector repPion(rep.px, rep.py, rep.pz, o2::constants::physics::MassPionCharged);
                    const auto mixedMother = targetK0s + repPion;
                    const auto mass = static_cast<float>(mixedMother.M());
                    const auto pt = static_cast<float>(mixedMother.Pt());
                    const auto y = static_cast<float>(mixedMother.Rapidity());
                    if (!std::isfinite(mass) || !std::isfinite(pt) || !std::isfinite(y)) {
                      continue;
                    }
                    if (std::abs(y) >= kstarCuts.cfgKstarRapidityMax.value || mass < dataDrivenMixing.cfgCandidateMassMin.value || mass >= dataDrivenMixing.cfgCandidateMassMax.value) {
                      continue;
                    }

                    float saMix = 0.0f;
                    if (!computeSA(targetK0s, mixedMother, evt.psiSA, saMix)) {
                      continue;
                    }
                    pionMatches.push_back(MixedMatch{rep.eventKey, rep.runNumber, rep.psi2, rep.posZ, rep.centrality, mass, pt, saMix, dpt, deta, dphi, rep.pt, rep.eta, rep.phi});
                  }
                }
              }
            }
          }

          const uint64_t pionSeed = dataDrivenMixing.cfgMixSeed.value ^ evt.eventKey ^ splitmix64DD(++targetOrdinal) ^ 0x50494f4eULL;
          keepDistinctEvents(pionMatches, pionSeed);
          pionMatches = downsampleMatches(std::move(pionMatches), pionSeed ^ 0x2222ULL);

          histos.fill(HIST("DataDrivenQA/hNMatchesK0s"), static_cast<float>(k0Matches.size()));
          histos.fill(HIST("DataDrivenQA/hNMatchesPion"), static_cast<float>(pionMatches.size()));

          const bool hasK0 = !k0Matches.empty();
          const bool hasPi = !pionMatches.empty();
          const bool hasAny = hasK0 || hasPi;
          histos.fill(HIST("DataDrivenQA/pK0sReplacementFoundVsPt"), targetPt, hasK0 ? 1.0f : 0.0f);
          histos.fill(HIST("DataDrivenQA/pK0sReplacementFoundVsCosTheta"), targetCosThetaStar, hasK0 ? 1.0f : 0.0f);
          histos.fill(HIST("DataDrivenQA/pK0sReplacementFoundVsCent"), evt.centrality, hasK0 ? 1.0f : 0.0f);
          histos.fill(HIST("DataDrivenQA/pK0sReplacementFoundVsMass"), targetMass, hasK0 ? 1.0f : 0.0f);
          histos.fill(HIST("DataDrivenQA/pPionReplacementFoundVsPt"), targetPt, hasPi ? 1.0f : 0.0f);
          histos.fill(HIST("DataDrivenQA/pPionReplacementFoundVsCosTheta"), targetCosThetaStar, hasPi ? 1.0f : 0.0f);
          histos.fill(HIST("DataDrivenQA/pPionReplacementFoundVsCent"), evt.centrality, hasPi ? 1.0f : 0.0f);
          histos.fill(HIST("DataDrivenQA/pPionReplacementFoundVsMass"), targetMass, hasPi ? 1.0f : 0.0f);
          histos.fill(HIST("DataDrivenQA/pAnyReplacementFoundVsCent"), evt.centrality, hasAny ? 1.0f : 0.0f);
          histos.fill(HIST("DataDrivenQA/pAnyReplacementFoundVsMass"), targetMass, hasAny ? 1.0f : 0.0f);

          if (!hasAny) {
            continue;
          }

          if (dataDrivenMixing.cfgFillMatchedSESparses.value) {
            if (hasK0) {
              histos.fill(HIST("hSparseSEMatchedK0s"), targetMass, targetPt, evt.centrality, targetSA, charge);
            }
            if (hasPi) {
              histos.fill(HIST("hSparseSEMatchedPion"), targetMass, targetPt, evt.centrality, targetSA, charge);
            }
          }

          const float k0LegNorm = hasK0 ? (hasPi ? 0.5f : 1.0f) : 0.0f;
          const float piLegNorm = hasPi ? (hasK0 ? 0.5f : 1.0f) : 0.0f;

          if (hasK0) {
            histos.fill(HIST("DataDrivenQA/hTargetK0sKinematics"), targetK0.pt, targetK0.eta, targetK0.phi, k0LegNorm);
            const float w = k0LegNorm / static_cast<float>(k0Matches.size());
            for (const auto& m : k0Matches) {
              histos.fill(HIST("hSparseMEK0sReplaced"), m.mass, m.pt, evt.centrality, m.saVal, charge, w);
              histos.fill(HIST("DataDrivenQA/hDeltaPtK0s"), m.dPt, w);
              histos.fill(HIST("DataDrivenQA/hDeltaEtaK0s"), m.dEta, w);
              histos.fill(HIST("DataDrivenQA/hDeltaPhiK0s"), m.dPhi, w);
              histos.fill(HIST("DataDrivenQA/hDeltaVz"), m.replacementVz - evt.posZ, w);
              histos.fill(HIST("DataDrivenQA/hDeltaCent"), m.replacementCent - evt.centrality, w);
              histos.fill(HIST("DataDrivenQA/hDeltaPsi2"), wrappedDeltaPsi2DD(m.replacementPsi2, evt.psi2), w);
              histos.fill(HIST("DataDrivenQA/hReplacementK0sKinematics"), m.repPt, m.repEta, m.repPhi, w);
            }
          }

          if (hasPi) {
            histos.fill(HIST("DataDrivenQA/hTargetPionKinematics"), targetPi.pt, targetPi.eta, targetPi.phi, piLegNorm);
            const float w = piLegNorm / static_cast<float>(pionMatches.size());
            for (const auto& m : pionMatches) {
              histos.fill(HIST("hSparseMEPionReplaced"), m.mass, m.pt, evt.centrality, m.saVal, charge, w);
              histos.fill(HIST("DataDrivenQA/hDeltaPtPion"), m.dPt, w);
              histos.fill(HIST("DataDrivenQA/hDeltaEtaPion"), m.dEta, w);
              histos.fill(HIST("DataDrivenQA/hDeltaPhiPion"), m.dPhi, w);
              histos.fill(HIST("DataDrivenQA/hDeltaVz"), m.replacementVz - evt.posZ, w);
              histos.fill(HIST("DataDrivenQA/hDeltaCent"), m.replacementCent - evt.centrality, w);
              histos.fill(HIST("DataDrivenQA/hDeltaPsi2"), wrappedDeltaPsi2DD(m.replacementPsi2, evt.psi2), w);
              histos.fill(HIST("DataDrivenQA/hReplacementPionKinematics"), m.repPt, m.repEta, m.repPhi, w);
            }
          }
        }
      }
    }

    histos.fill(HIST("DataDrivenQA/hK0sPoolCandidates"), static_cast<float>(ddK0sPoolCandidateCount));
    histos.fill(HIST("DataDrivenQA/hPionPoolCandidates"), static_cast<float>(ddPionPoolCandidateCount));

    if (dataDrivenMixing.cfgLogPool.value) {
      LOGF(info, "DD daughter pools: dataframe=%" PRIu64 " K0S(cells=%zu,candidates=%zu) pion(cells=%zu,candidates=%zu)", currentDataFrame, ddRollingK0sPool.size(), ddK0sPoolCandidateCount, ddRollingPionPool.size(), ddPionPoolCandidateCount);
    }
  }
  PROCESS_SWITCH(chKstarSpin, processDataDrivenME, "Process DD daughter-replacement event mixing", true);

  // ------------------------------------------------------------------
  // MC efficiency process (processMC):
  // Reconstructed and generated K*(892)+- -> K0S + pi+- in the chosen
  // coordinate frame (Event Plane, Production Plane, or Random EP)
  // ------------------------------------------------------------------
  void processMC(CollisionMCTrueTable::iterator const& trueCollision, CollisionMCRecTableCentFT0C const& recCollisions, TrackMCTrueTable const& genParticles, TrackMCRecTable const& recTracks, ResoV0s const& v0s)
  {
    histos.fill(HIST("MC/hMC"), 0.5);
    if (recCollisions.size() == 0) {
      histos.fill(HIST("MC/hMC"), 1.5);
      return;
    }
    if (recCollisions.size() > 1) {
      histos.fill(HIST("MC/hMC"), 2.5);
      return;
    }

    const float psiEP = getMCEffPsi(trueCollision, genParticles);

    for (const auto& recCollision : recCollisions) {
      histos.fill(HIST("MC/hMC"), 3.5);
      if (!recCollision.sel8()) {
        histos.fill(HIST("MC/hMC"), 4.5);
        continue;
      }
      if (!passEventQualityMC(recCollision)) {
        histos.fill(HIST("MC/hMC"), 5.5);
        continue;
      }
      histos.fill(HIST("MC/hMC"), 6.5);

      const float centrality = recCollision.centFT0C();
      histos.fill(HIST("MC/CentPercentileMCRecHist"), centrality);

      const float psiSA = getSAEventAngle(psiEP);
      histos.fill(HIST("MC/hMCEffPsi"), centrality, psiSA);

      const auto recTracksThis = recTracks.sliceBy(mcTracksPerCollision, recCollision.globalIndex());
      const auto v0sThis = v0s.sliceBy(v0sPerCollisionDD, recCollision.globalIndex());

      std::vector<int64_t> seenMothers;

      // Reconstructed true K*(892)+- loop
      for (const auto& v0 : v0sThis) {
        if (!selectionK0s<TrackMCRecTable>(recCollision, v0)) {
          continue;
        }

        auto posTrack = v0.template posTrack_as<TrackMCRecTable>();
        auto negTrack = v0.template negTrack_as<TrackMCRecTable>();
        const auto posID = posTrack.globalIndex();
        const auto negID = negTrack.globalIndex();

        const ROOT::Math::PxPyPzMVector k0sRec(v0.px(), v0.py(), v0.pz(), o2::constants::physics::MassK0Short);

        for (const auto& bachelor : recTracksThis) {
          if (!selectionBachelorPion(bachelor)) {
            continue;
          }
          const auto bachelorID = bachelor.globalIndex();
          if (bachelorID == posID || bachelorID == negID) {
            continue;
          }

          int64_t motherIndex = -1;
          ROOT::Math::PxPyPzMVector k0sGen, pionGen;
          if (!matchRecoChKstar(v0, bachelor, motherIndex, k0sGen, pionGen)) {
            continue;
          }

          if (mcEffEP.avoidsplitrackMC.value) {
            if (std::find(seenMothers.begin(), seenMothers.end(), motherIndex) != seenMothers.end()) {
              continue;
            }
            seenMothers.push_back(motherIndex);
          }

          const ROOT::Math::PxPyPzMVector pionRec(bachelor.px(), bachelor.py(), bachelor.pz(), o2::constants::physics::MassPionCharged);
          const auto kstarRec = k0sRec + pionRec;

          const auto rapRec = static_cast<float>(kstarRec.Rapidity());
          if (!std::isfinite(rapRec) || std::abs(rapRec) >= kstarCuts.cfgKstarRapidityMax.value) {
            continue;
          }

          float saRec = 0.0f;
          float cosThetaStarRec = 0.0f;
          if (!computeSA(k0sRec, kstarRec, psiSA, saRec, &cosThetaStarRec)) {
            continue;
          }

          const auto massRec = static_cast<float>(kstarRec.M());
          const auto ptRec = static_cast<float>(kstarRec.Pt());
          const auto chargeRec = static_cast<float>(bachelor.sign());

          histos.fill(HIST("MC/hSparseChKstarMCRecSA"), massRec, ptRec, centrality, saRec, chargeRec);
          histos.fill(HIST("MC/hSparseChKstarMCRecCosThetaStar_effy"), massRec, ptRec, centrality, cosThetaStarRec, chargeRec);

          if (mcEffEP.fillPhiVMinusPsi.value) {
            histos.fill(HIST("MC/hSparseChKstarMCRecCosThetaStarPhiVPsi"), ptRec, cosThetaStarRec, getPhiInRange(kstarRec.Phi() - psiSA), centrality);
          }
          if (closure.doClosureTest.value) {
            histos.fill(HIST("closure/ChKstarMCRecTrueCosThetaStar"), massRec, ptRec, centrality, cosThetaStarRec, chargeRec);
          }
        }
      }

      // Generated K*(892)+- loop
      for (const auto& mcParticle : genParticles) {
        ROOT::Math::PxPyPzMVector k0sGen, pionGen;
        float chargeGen = 0.0f;
        if (!extractGenChKstarDaughters(mcParticle, mcEffEP.genacceptancecut.value, k0sGen, pionGen, chargeGen)) {
          continue;
        }

        const auto kstarGen = k0sGen + pionGen;
        const auto rapGen = static_cast<float>(kstarGen.Rapidity());
        if (!std::isfinite(rapGen) || std::abs(rapGen) >= kstarCuts.cfgKstarRapidityMax.value) {
          continue;
        }

        float saGen = 0.0f;
        float cosThetaStarGen = 0.0f;
        if (!computeSA(k0sGen, kstarGen, psiSA, saGen, &cosThetaStarGen)) {
          continue;
        }

        const auto massGen = static_cast<float>(kstarGen.M());
        const auto ptGen = static_cast<float>(kstarGen.Pt());

        histos.fill(HIST("MC/hSparseChKstarMCGenSA"), massGen, ptGen, centrality, saGen, chargeGen);
        histos.fill(HIST("MC/hSparseChKstarMCGenCosThetaStar_effy"), massGen, ptGen, centrality, cosThetaStarGen, chargeGen);

        if (mcEffEP.fillPhiVMinusPsi.value) {
          histos.fill(HIST("MC/hSparseChKstarMCGenCosThetaStarPhiVPsi"), ptGen, cosThetaStarGen, getPhiInRange(kstarGen.Phi() - psiSA), centrality);
        }
        if (closure.doClosureTest.value) {
          histos.fill(HIST("closure/ChKstarMCGenCosThetaStar"), massGen, ptGen, centrality, cosThetaStarGen, chargeGen);
        }
      }
    }
  }
  PROCESS_SWITCH(chKstarSpin, processMC, "Process MC efficiency for charged K*(892)", false);

  // ------------------------------------------------------------------
  // Event loss and signal loss in MC for K*(892)+-
  // ------------------------------------------------------------------
  void processEvtLossSigLossMC(McCollisionMults::iterator const& mcCollision, CollisionMCRecTableCentFT0C const& recCollisions, TrackMCTrueTable const& genParticles)
  {
    if (evtSigLoss.cutVzGen.value && std::abs(mcCollision.posZ()) > eventCuts.cfgCutVertexZ.value) {
      return;
    }
    if (evtSigLoss.isApplyInelgt0.value && !mcCollision.isInelGt0()) {
      return;
    }
    if (evtSigLoss.isApplyTVX.value && (mcCollision.multMCFT0C() <= 0 || mcCollision.multMCFT0A() <= 0)) {
      return;
    }

    const float impactPar = mcCollision.impactParameter();
    const float multMC = mcCollision.multMCNParticlesEta05();

    histos.fill(HIST("evtSigLoss/MCEventHist"), 1);
    histos.fill(HIST("evtSigLoss/hImpactParameterGen"), impactPar);
    histos.fill(HIST("evtSigLoss/hMultEta05Gen"), multMC);

    if (recCollisions.size() == 0) {
      histos.fill(HIST("evtSigLoss/MCEventHist"), 3);
      histos.fill(HIST("evtSigLoss/hImpactParameterGenNoReco"), impactPar);
      histos.fill(HIST("evtSigLoss/hMultEta05GenNoReco"), multMC);
    }
    if (recCollisions.size() > 1) {
      histos.fill(HIST("evtSigLoss/MCEventHist"), 4);
    }

    bool selected = false;
    float centrality = -999.0f;
    if (!evtSigLoss.requireSingleReco.value || recCollisions.size() == 1) {
      for (const auto& recCollision : recCollisions) {
        if (!passEventQualityMC(recCollision)) {
          continue;
        }
        selected = true;
        centrality = recCollision.centFT0C();
      }
    }

    if (selected) {
      histos.fill(HIST("evtSigLoss/MCEventHist"), 2);
      histos.fill(HIST("evtSigLoss/hImpactParameterRec"), impactPar);
      histos.fill(HIST("evtSigLoss/hMultEta05Rec"), multMC);
      histos.fill(HIST("evtSigLoss/hImpactParVsCentRec"), centrality, impactPar);
      histos.fill(HIST("evtSigLoss/hMultEta05VsCentRec"), centrality, multMC);
    }

    for (const auto& mcParticle : genParticles) {
      ROOT::Math::PxPyPzMVector k0sGen, pionGen;
      float chargeGen = 0.0f;
      if (!extractGenChKstarDaughters(mcParticle, false, k0sGen, pionGen, chargeGen)) {
        continue;
      }

      const double pt = (k0sGen + pionGen).Pt();
      histos.fill(HIST("evtSigLoss/hChKstarGenBeforeEvtSel"), pt, impactPar);
      histos.fill(HIST("evtSigLoss/hChKstarGenVsMultBeforeEvtSel"), pt, multMC);
      if (selected) {
        histos.fill(HIST("evtSigLoss/hChKstarGenAfterEvtSel"), pt, impactPar);
        histos.fill(HIST("evtSigLoss/hChKstarGenVsMultAfterEvtSel"), pt, multMC);
        histos.fill(HIST("evtSigLoss/hChKstarGenAfterEvtSelVsCent"), pt, centrality);
      }
    }
  }
  PROCESS_SWITCH(chKstarSpin, processEvtLossSigLossMC, "Process event loss and signal loss for charged K*(892)", false);

  // ------------------------------------------------------------------
  // Dedicated MC Closure Test (True EP): SE, RecTrue, Gen, and Rotated
  // ------------------------------------------------------------------
  void processMCClosure(CollisionMCTrueTable::iterator const& trueCollision, CollisionMCRecTableCentFT0C const& recCollisions, TrackMCTrueTable const& genParticles, TrackMCRecTable const& recTracks, ResoV0s const& v0s)
  {
    histos.fill(HIST("closureChKstar/hMC"), 0.5);
    if (recCollisions.size() == 0) {
      histos.fill(HIST("closureChKstar/hMC"), 1.5);
      return;
    }
    if (recCollisions.size() > 1) {
      histos.fill(HIST("closureChKstar/hMC"), 2.5);
      return;
    }

    const bool isEffSample = isClosureEffSample(trueCollision.globalIndex());
    const bool isTestSample = isClosureTestSample(trueCollision.globalIndex());
    const float psi = getClosurePsi(trueCollision.eventPlaneAngle());

    for (const auto& recCollision : recCollisions) {
      if (!passEventQualityMC(recCollision)) {
        continue;
      }
      if (isEffSample) {
        histos.fill(HIST("closureChKstar/hMC"), 3.5);
      }
      if (isTestSample) {
        histos.fill(HIST("closureChKstar/hMC"), 4.5);
      }

      const float centrality = recCollision.centFT0C();
      const auto recTracksThis = recTracks.sliceBy(mcTracksPerCollision, recCollision.globalIndex());
      const auto v0sThis = v0s.sliceBy(v0sPerCollisionDD, recCollision.globalIndex());

      std::vector<int64_t> seenMothersEff;
      std::vector<int64_t> seenMothersTest;

      // Reconstructed K0S + bachelor pi+- pairs
      for (const auto& v0 : v0sThis) {
        if (!selectionK0s<TrackMCRecTable>(recCollision, v0)) {
          continue;
        }

        auto posTrack = v0.template posTrack_as<TrackMCRecTable>();
        auto negTrack = v0.template negTrack_as<TrackMCRecTable>();
        const auto posID = posTrack.globalIndex();
        const auto negID = negTrack.globalIndex();

        const ROOT::Math::PxPyPzMVector k0s(v0.px(), v0.py(), v0.pz(), o2::constants::physics::MassK0Short);

        for (const auto& bachelor : recTracksThis) {
          if (!selectionBachelorPion(bachelor)) {
            continue;
          }
          if (closure.requirePrimaryTracks.value && (!bachelor.has_mcParticle() || !bachelor.mcParticle().isPhysicalPrimary())) {
            continue;
          }

          const auto bachelorID = bachelor.globalIndex();
          if (bachelorID == posID || bachelorID == negID) {
            continue;
          }

          const ROOT::Math::PxPyPzMVector pion(bachelor.px(), bachelor.py(), bachelor.pz(), o2::constants::physics::MassPionCharged);
          const auto kstar = k0s + pion;
          const auto rap = static_cast<float>(kstar.Rapidity());
          if (!std::isfinite(rap) || std::abs(rap) >= kstarCuts.cfgKstarRapidityMax.value) {
            continue;
          }

          float saRec = 0.0f;
          float cosRec = 0.0f;
          if (!computeSA(k0s, kstar, psi, saRec, &cosRec)) {
            continue;
          }

          const auto mass = static_cast<float>(kstar.M());
          const auto pt = static_cast<float>(kstar.Pt());
          const auto charge = static_cast<float>(bachelor.sign());

          double weight = 1.0;
          int64_t motherIndex = -1;
          ROOT::Math::PxPyPzMVector k0sGen, pionGen;
          if (matchRecoChKstar(v0, bachelor, motherIndex, k0sGen, pionGen)) {
            const auto kstarGen = k0sGen + pionGen;
            float saGen = 0.0f;
            float cosGen = 0.0f;
            if (computeSA(k0sGen, kstarGen, psi, saGen, &cosGen)) {
              weight = getInjectionWeight(cosGen);

              if (isEffSample) {
                bool isDupEff = mcEffEP.avoidsplitrackMC.value && (std::find(seenMothersEff.begin(), seenMothersEff.end(), motherIndex) != seenMothersEff.end());
                if (!isDupEff) {
                  if (mcEffEP.avoidsplitrackMC.value) {
                    seenMothersEff.push_back(motherIndex);
                  }
                  histos.fill(HIST("closureChKstar/hSparseRecTrueEff"), mass, pt, centrality, cosRec, charge, 1.0);
                }
              }

              if (isTestSample) {
                bool isDupTest = mcEffEP.avoidsplitrackMC.value && (std::find(seenMothersTest.begin(), seenMothersTest.end(), motherIndex) != seenMothersTest.end());
                if (!isDupTest) {
                  if (mcEffEP.avoidsplitrackMC.value) {
                    seenMothersTest.push_back(motherIndex);
                  }
                  histos.fill(HIST("closureChKstar/hSparseRecTrue"), mass, pt, centrality, cosRec, charge, weight);
                  histos.fill(HIST("closureChKstar/hCosThetaStarGenVsRec"), cosGen, cosRec);
                }
              }
            }
          }

          if (!isTestSample) {
            continue;
          }

          histos.fill(HIST("closureChKstar/hSparseSE"), mass, pt, centrality, cosRec, charge, weight);

          // Optional rotational background in closure test sample
          if (closure.fillRotBkg.value && saConfig.nBkgRotations.value > 1) {
            const int nRot = saConfig.nBkgRotations.value;
            const float angleStart = saConfig.confMinRot.value;
            const float angleStep = (saConfig.confMaxRot.value - angleStart) / static_cast<float>(nRot - 1);
            for (int irot = 0; irot < nRot; ++irot) {
              const float rotAngle = angleStart + static_cast<float>(irot) * angleStep;
              const float cosA = std::cos(rotAngle);
              const float sinA = std::sin(rotAngle);

              ROOT::Math::PxPyPzMVector k0sRot = k0s;
              ROOT::Math::PxPyPzMVector pionRot = pion;
              if (saConfig.cfgRotatePion.value) {
                pionRot = ROOT::Math::PxPyPzMVector(bachelor.px() * cosA - bachelor.py() * sinA, bachelor.px() * sinA + bachelor.py() * cosA, bachelor.pz(), o2::constants::physics::MassPionCharged);
              } else {
                k0sRot = ROOT::Math::PxPyPzMVector(v0.px() * cosA - v0.py() * sinA, v0.px() * sinA + v0.py() * cosA, v0.pz(), o2::constants::physics::MassK0Short);
              }

              const auto kstarRot = k0sRot + pionRot;
              const auto rapRot = static_cast<float>(kstarRot.Rapidity());
              if (!std::isfinite(rapRot) || std::abs(rapRot) >= kstarCuts.cfgKstarRapidityMax.value) {
                continue;
              }

              float saRot = 0.0f;
              float cosRot = 0.0f;
              if (!computeSA(k0sRot, kstarRot, psi, saRot, &cosRot)) {
                continue;
              }
              histos.fill(HIST("closureChKstar/hSparseRot"), static_cast<float>(kstarRot.M()), static_cast<float>(kstarRot.Pt()), centrality, cosRot, charge);
            }
          }
        }
      }

      // Generated K*(892)+- loop for closure denominator (eff sample) and input truth (test sample)
      for (const auto& mcParticle : genParticles) {
        ROOT::Math::PxPyPzMVector k0sGen, pionGen;
        float chargeGen = 0.0f;
        if (!extractGenChKstarDaughters(mcParticle, closure.genAcceptanceCut.value, k0sGen, pionGen, chargeGen)) {
          continue;
        }

        const auto kstarGen = k0sGen + pionGen;
        const auto rapGen = static_cast<float>(kstarGen.Rapidity());
        if (!std::isfinite(rapGen) || std::abs(rapGen) >= kstarCuts.cfgKstarRapidityMax.value) {
          continue;
        }

        float saGen = 0.0f;
        float cosGen = 0.0f;
        if (!computeSA(k0sGen, kstarGen, psi, saGen, &cosGen)) {
          continue;
        }

        const auto massGen = static_cast<float>(kstarGen.M());
        const auto ptGen = static_cast<float>(kstarGen.Pt());

        if (isEffSample) {
          histos.fill(HIST("closureChKstar/hSparseGenEff"), massGen, ptGen, centrality, cosGen, chargeGen, 1.0);
        }
        if (isTestSample) {
          const double weight = getInjectionWeight(cosGen);
          histos.fill(HIST("closureChKstar/hSparseGen"), massGen, ptGen, centrality, cosGen, chargeGen, weight);
        }
      }
    }
  }
  PROCESS_SWITCH(chKstarSpin, processMCClosure, "Process MC closure of cos(theta*) for charged K*: SE, true rec, gen", false);

  // ------------------------------------------------------------------
  // Dedicated MC Closure Test (True EP): Mixed Event
  // ------------------------------------------------------------------
  void processMCClosureME(CollisionMCRecMixTable const& /*collisions*/, TrackMCRecTable const& /*tracks*/, ResoV0s const& /*v0s*/, aod::McCollisions const& mcCollisions, aod::McCollisionLabels const& mcLabels, aod::McParticles const&)
  {
    const auto nRecPerMc = countRecoPerMcCollision(mcLabels, mcCollisions.size());

    for (const auto& [collision1, tracks1, collision2, v0s2] : mcClosureMixedEventPairs) {
      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }
      if (!passEventQualityMC(collision1) || !passEventQualityMC(collision2)) {
        continue;
      }
      if (!collision1.has_mcCollision() || !collision2.has_mcCollision()) {
        continue;
      }
      if (!isClosureTestSample(collision1.mcCollisionId()) || !isClosureTestSample(collision2.mcCollisionId())) {
        continue;
      }
      if (closure.meRequireSingleReco.value && (nRecPerMc[collision1.mcCollisionId()] != 1 || nRecPerMc[collision2.mcCollisionId()] != 1)) {
        continue;
      }

      const float centrality = collision1.centFT0C();
      const float psi = getClosurePsi(collision1.mcCollision().eventPlaneAngle());

      for (const auto& v0 : v0s2) {
        if (!selectionK0s<TrackMCRecTable>(collision2, v0)) {
          continue;
        }
        const ROOT::Math::PxPyPzMVector k0s(v0.px(), v0.py(), v0.pz(), o2::constants::physics::MassK0Short);

        for (const auto& bachelor : tracks1) {
          if (!selectionBachelorPion(bachelor)) {
            continue;
          }
          if (closure.requirePrimaryTracks.value && (!bachelor.has_mcParticle() || !bachelor.mcParticle().isPhysicalPrimary())) {
            continue;
          }

          const ROOT::Math::PxPyPzMVector pion(bachelor.px(), bachelor.py(), bachelor.pz(), o2::constants::physics::MassPionCharged);
          const auto kstar = k0s + pion;
          const auto rap = static_cast<float>(kstar.Rapidity());
          if (!std::isfinite(rap) || std::abs(rap) >= kstarCuts.cfgKstarRapidityMax.value) {
            continue;
          }

          float saRec = 0.0f;
          float cosRec = 0.0f;
          if (!computeSA(k0s, kstar, psi, saRec, &cosRec)) {
            continue;
          }

          histos.fill(HIST("closureChKstar/hSparseME"), static_cast<float>(kstar.M()), static_cast<float>(kstar.Pt()), centrality, cosRec, static_cast<float>(bachelor.sign()));
        }
      }
    }
  }
  PROCESS_SWITCH(chKstarSpin, processMCClosureME, "Process MC closure of cos(theta*) for charged K*: mixed event", false);

  // ------------------------------------------------------------------
  // Dedicated MC Closure Test (Reconstructed FT0C EP): SE, RecTrue, Gen
  // ------------------------------------------------------------------
  void processMCClosureRecoEP(CollisionMCTrueTable::iterator const& trueCollision, CollisionMCRecTableClosure const& recCollisions, TrackMCTrueTable const& genParticles, TrackMCRecTable const& recTracks, ResoV0s const& v0s)
  {
    histos.fill(HIST("closureChKstarRecoEP/hMC"), 0.5);
    if (recCollisions.size() == 0) {
      histos.fill(HIST("closureChKstarRecoEP/hMC"), 1.5);
      return;
    }
    if (recCollisions.size() > 1) {
      histos.fill(HIST("closureChKstarRecoEP/hMC"), 2.5);
      return;
    }

    const bool isEffSample = isClosureEffSample(trueCollision.globalIndex());
    const bool isTestSample = isClosureTestSample(trueCollision.globalIndex());
    const float psiTrue = getSAEventAngle(trueCollision.eventPlaneAngle());

    for (const auto& recCollision : recCollisions) {
      if (!passEventQualityMC(recCollision)) {
        continue;
      }
      const float recoPsiFT0C = recCollision.psiFT0C();
      if (!recCollision.triggereventep() || !std::isfinite(recoPsiFT0C)) {
        continue;
      }

      if (isEffSample) {
        histos.fill(HIST("closureChKstarRecoEP/hMC"), 3.5);
      }
      if (isTestSample) {
        histos.fill(HIST("closureChKstarRecoEP/hMC"), 4.5);
      }

      const float centrality = recCollision.centFT0C();
      // const float psi = getSAEventAngle(recoPsiFT0C);
      const float psi = (saConfig.cfgSAFrame.value == kRandomEventPlane) ? psiTrue : getSAEventAngle(recoPsiFT0C);

      {
        const double psiTPC = recCollision.psiTPC();
        const double psiFT0A = recCollision.psiFT0A();
        const double qFT0C = recCollision.qFT0C();
        const double qTPC = recCollision.qTPC();
        const double qFT0A = recCollision.qFT0A();
        histos.fill(HIST("closureChKstarRecoEP/ResFT0CTPCSP"), centrality, qFT0C * qTPC * std::cos(2.0 * (recoPsiFT0C - psiTPC)));
        histos.fill(HIST("closureChKstarRecoEP/ResFT0CFT0ASP"), centrality, qFT0C * qFT0A * std::cos(2.0 * (recoPsiFT0C - psiFT0A)));
        histos.fill(HIST("closureChKstarRecoEP/ResFT0ATPCSP"), centrality, qTPC * qFT0A * std::cos(2.0 * (psiTPC - psiFT0A)));
        histos.fill(HIST("closureChKstarRecoEP/ResFT0CTPC"), centrality, std::cos(2.0 * (recoPsiFT0C - psiTPC)));
        histos.fill(HIST("closureChKstarRecoEP/ResFT0CFT0A"), centrality, std::cos(2.0 * (recoPsiFT0C - psiFT0A)));
        histos.fill(HIST("closureChKstarRecoEP/ResFT0ATPC"), centrality, std::cos(2.0 * (psiTPC - psiFT0A)));
        histos.fill(HIST("closureChKstarRecoEP/ResTrueFT0C"), centrality, std::cos(2.0 * (recoPsiFT0C - trueCollision.eventPlaneAngle())));
      }

      const auto recTracksThis = recTracks.sliceBy(mcTracksPerCollision, recCollision.globalIndex());
      const auto v0sThis = v0s.sliceBy(v0sPerCollisionDD, recCollision.globalIndex());

      std::vector<int64_t> seenMothersEff;
      std::vector<int64_t> seenMothersTest;

      // Reconstructed pairs
      for (const auto& v0 : v0sThis) {
        if (!selectionK0s<TrackMCRecTable>(recCollision, v0)) {
          continue;
        }

        auto posTrack = v0.template posTrack_as<TrackMCRecTable>();
        auto negTrack = v0.template negTrack_as<TrackMCRecTable>();
        const auto posID = posTrack.globalIndex();
        const auto negID = negTrack.globalIndex();

        const ROOT::Math::PxPyPzMVector k0s(v0.px(), v0.py(), v0.pz(), o2::constants::physics::MassK0Short);

        for (const auto& bachelor : recTracksThis) {
          if (!selectionBachelorPion(bachelor)) {
            continue;
          }
          if (closure.requirePrimaryTracks.value && (!bachelor.has_mcParticle() || !bachelor.mcParticle().isPhysicalPrimary())) {
            continue;
          }

          const auto bachelorID = bachelor.globalIndex();
          if (bachelorID == posID || bachelorID == negID) {
            continue;
          }

          const ROOT::Math::PxPyPzMVector pion(bachelor.px(), bachelor.py(), bachelor.pz(), o2::constants::physics::MassPionCharged);
          const auto kstar = k0s + pion;
          const auto rap = static_cast<float>(kstar.Rapidity());
          if (!std::isfinite(rap) || std::abs(rap) >= kstarCuts.cfgKstarRapidityMax.value) {
            continue;
          }

          float saRec = 0.0f;
          float cosRec = 0.0f;
          if (!computeSA(k0s, kstar, psi, saRec, &cosRec)) {
            continue;
          }

          const auto mass = static_cast<float>(kstar.M());
          const auto pt = static_cast<float>(kstar.Pt());
          const auto charge = static_cast<float>(bachelor.sign());

          double weight = 1.0;
          int64_t motherIndex = -1;
          ROOT::Math::PxPyPzMVector k0sGen, pionGen;
          if (matchRecoChKstar(v0, bachelor, motherIndex, k0sGen, pionGen)) {
            const auto kstarGen = k0sGen + pionGen;
            float saGen = 0.0f;
            float cosGen = 0.0f;
            if (computeSA(k0sGen, kstarGen, psiTrue, saGen, &cosGen)) {
              weight = getInjectionWeight(cosGen);

              if (isEffSample) {
                bool isDupEff = mcEffEP.avoidsplitrackMC.value && (std::find(seenMothersEff.begin(), seenMothersEff.end(), motherIndex) != seenMothersEff.end());
                if (!isDupEff) {
                  if (mcEffEP.avoidsplitrackMC.value) {
                    seenMothersEff.push_back(motherIndex);
                  }
                  histos.fill(HIST("closureChKstarRecoEP/hSparseRecTrueEff"), mass, pt, centrality, cosRec, charge, 1.0);
                }
              }

              if (isTestSample) {
                bool isDupTest = mcEffEP.avoidsplitrackMC.value && (std::find(seenMothersTest.begin(), seenMothersTest.end(), motherIndex) != seenMothersTest.end());
                if (!isDupTest) {
                  if (mcEffEP.avoidsplitrackMC.value) {
                    seenMothersTest.push_back(motherIndex);
                  }
                  histos.fill(HIST("closureChKstarRecoEP/hSparseRecTrue"), mass, pt, centrality, cosRec, charge, weight);
                  histos.fill(HIST("closureChKstarRecoEP/hCosThetaStarGenVsRec"), cosGen, cosRec);
                }
              }
            }
          }

          if (!isTestSample) {
            continue;
          }

          histos.fill(HIST("closureChKstarRecoEP/hSparseSE"), mass, pt, centrality, cosRec, charge, weight);

          // Optional rotational background in RecoEP closure test sample
          if (closure.fillRotBkg.value && saConfig.nBkgRotations.value > 1) {
            const int nRot = saConfig.nBkgRotations.value;
            const float angleStart = saConfig.confMinRot.value;
            const float angleStep = (saConfig.confMaxRot.value - angleStart) / static_cast<float>(nRot - 1);
            for (int irot = 0; irot < nRot; ++irot) {
              const float rotAngle = angleStart + static_cast<float>(irot) * angleStep;
              const float cosA = std::cos(rotAngle);
              const float sinA = std::sin(rotAngle);

              ROOT::Math::PxPyPzMVector k0sRot = k0s;
              ROOT::Math::PxPyPzMVector pionRot = pion;
              if (saConfig.cfgRotatePion.value) {
                pionRot = ROOT::Math::PxPyPzMVector(bachelor.px() * cosA - bachelor.py() * sinA, bachelor.px() * sinA + bachelor.py() * cosA, bachelor.pz(), o2::constants::physics::MassPionCharged);
              } else {
                k0sRot = ROOT::Math::PxPyPzMVector(v0.px() * cosA - v0.py() * sinA, v0.px() * sinA + v0.py() * cosA, v0.pz(), o2::constants::physics::MassK0Short);
              }

              const auto kstarRot = k0sRot + pionRot;
              const auto rapRot = static_cast<float>(kstarRot.Rapidity());
              if (!std::isfinite(rapRot) || std::abs(rapRot) >= kstarCuts.cfgKstarRapidityMax.value) {
                continue;
              }

              float saRot = 0.0f;
              float cosRot = 0.0f;
              if (!computeSA(k0sRot, kstarRot, psi, saRot, &cosRot)) {
                continue;
              }
              histos.fill(HIST("closureChKstarRecoEP/hSparseRot"), static_cast<float>(kstarRot.M()), static_cast<float>(kstarRot.Pt()), centrality, cosRot, charge);
            }
          }
        }
      }

      // Generated K*(892)+- loop (evaluated with respect to psiTrue)
      for (const auto& mcParticle : genParticles) {
        ROOT::Math::PxPyPzMVector k0sGen, pionGen;
        float chargeGen = 0.0f;
        if (!extractGenChKstarDaughters(mcParticle, closure.genAcceptanceCut.value, k0sGen, pionGen, chargeGen)) {
          continue;
        }

        const auto kstarGen = k0sGen + pionGen;
        const auto rapGen = static_cast<float>(kstarGen.Rapidity());
        if (!std::isfinite(rapGen) || std::abs(rapGen) >= kstarCuts.cfgKstarRapidityMax.value) {
          continue;
        }

        float saGen = 0.0f;
        float cosGen = 0.0f;
        if (!computeSA(k0sGen, kstarGen, psiTrue, saGen, &cosGen)) {
          continue;
        }

        const auto massGen = static_cast<float>(kstarGen.M());
        const auto ptGen = static_cast<float>(kstarGen.Pt());

        if (isEffSample) {
          histos.fill(HIST("closureChKstarRecoEP/hSparseGenEff"), massGen, ptGen, centrality, cosGen, chargeGen, 1.0);
        }
        if (isTestSample) {
          const double weight = getInjectionWeight(cosGen);
          histos.fill(HIST("closureChKstarRecoEP/hSparseGen"), massGen, ptGen, centrality, cosGen, chargeGen, weight);
        }
      }
    }
  }
  PROCESS_SWITCH(chKstarSpin, processMCClosureRecoEP, "Process MC closure of cos(theta*) for charged K*, reconstructed Psi_FT0C: SE, true rec, gen", false);

  // ------------------------------------------------------------------
  // Dedicated MC Closure Test (Reconstructed FT0C EP): Mixed Event
  // ------------------------------------------------------------------
  void processMCClosureMERecoEP(CollisionMCRecMixTableRecoEP const& /*collisions*/, TrackMCRecTable const& /*tracks*/, ResoV0s const& /*v0s*/, aod::McCollisions const& mcCollisions, aod::McCollisionLabels const& mcLabels, aod::McParticles const&)
  {
    const auto nRecPerMc = countRecoPerMcCollision(mcLabels, mcCollisions.size());

    for (const auto& [collision1, tracks1, collision2, v0s2] : mcClosureRecoEPMixedEventPairs) {
      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }
      if (!passEventQualityMC(collision1) || !passEventQualityMC(collision2)) {
        continue;
      }
      if (!collision1.triggereventep() || !collision2.triggereventep() || !std::isfinite(collision1.psiFT0C()) || !std::isfinite(collision2.psiFT0C())) {
        continue;
      }
      if (!collision1.has_mcCollision() || !collision2.has_mcCollision()) {
        continue;
      }
      if (!isClosureTestSample(collision1.mcCollisionId()) || !isClosureTestSample(collision2.mcCollisionId())) {
        continue;
      }
      if (closure.meRequireSingleReco.value && (nRecPerMc[collision1.mcCollisionId()] != 1 || nRecPerMc[collision2.mcCollisionId()] != 1)) {
        continue;
      }

      const float centrality = collision1.centFT0C();
      const float psi = getSAEventAngle(collision1.psiFT0C());

      for (const auto& v0 : v0s2) {
        if (!selectionK0s<TrackMCRecTable>(collision2, v0)) {
          continue;
        }
        const ROOT::Math::PxPyPzMVector k0s(v0.px(), v0.py(), v0.pz(), o2::constants::physics::MassK0Short);

        for (const auto& bachelor : tracks1) {
          if (!selectionBachelorPion(bachelor)) {
            continue;
          }
          if (closure.requirePrimaryTracks.value && (!bachelor.has_mcParticle() || !bachelor.mcParticle().isPhysicalPrimary())) {
            continue;
          }

          const ROOT::Math::PxPyPzMVector pion(bachelor.px(), bachelor.py(), bachelor.pz(), o2::constants::physics::MassPionCharged);
          const auto kstar = k0s + pion;
          const auto rap = static_cast<float>(kstar.Rapidity());
          if (!std::isfinite(rap) || std::abs(rap) >= kstarCuts.cfgKstarRapidityMax.value) {
            continue;
          }

          float saRec = 0.0f;
          float cosRec = 0.0f;
          if (!computeSA(k0s, kstar, psi, saRec, &cosRec)) {
            continue;
          }

          histos.fill(HIST("closureChKstarRecoEP/hSparseME"), static_cast<float>(kstar.M()), static_cast<float>(kstar.Pt()), centrality, cosRec, static_cast<float>(bachelor.sign()));
        }
      }
    }
  }
  PROCESS_SWITCH(chKstarSpin, processMCClosureMERecoEP, "Process MC closure of cos(theta*) for charged K*, reconstructed Psi_FT0C: mixed event", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<chKstarSpin>(cfgc)}; }
